
  SUBROUTINE solve_tracer_upwind_sed(ctf, ctfold)
  !VF
  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM
  !a USE g_PARSUP
  IMPLICIT NONE
  real(kind=WP), intent(inout) :: ctf(nsigma-1, nod2D)
  real(kind=WP), intent(inout) :: ctfold(nsigma-1, nod2D)
  integer      :: nod1, nod2, n, nz, ed, el1, el2
  integer      :: i
  real(kind=WP) :: c1, c2, c1_c, c2_c
  real(kind=WP) :: cvert(nsigma),  db, dg
  real(kind=WP) :: dmean, un1, s, c, aux
  real(kind=WP) :: Tx, Ty, Tmean, Kh
  real(kind=WP) :: Vel_nor(nod2D), rho
  real(kind=WP) :: scalearea_inv, area_inv
  real(kind=WP) :: c_aux, aux_c
  real(kind=WP) :: ttxa(nsigma-1,elem2D), ttya(nsigma-1,elem2D)
  real(kind=WP) :: ctrhs(nsigma-1,nod2D),trhs_c(nsigma-1,nod2D)

  ! =================
  ! AB interpolation
  ! =================
call vel_gradients
call h_viscosity
DO n=1, nod2D
   ctfold(:,n) = -(0.5_WP+epsilon)*ctfold(:,n)+(1.5_WP+epsilon)*ctf(:,n)  - c_aver
   ctf(:,n)    = ctf(:,n) - c_aver
END DO
     ! ctfold now contains AB-interpolated tracer! (n+1/2)


   ctrhs  = 0.0_WP
   trhs_c = 0.0_WP
   ttxa = 0.0_WP
   ttya = 0.0_WP

  call tracer_gradient_elements(ctf, ttxa, ttya)
  !call fill_up_dn_grad(ttxa,ttya)



  ! =================
  ! Horizontal advection and diffusion and
  !  Average diffusive flux
  ! =================

  DO ed=1, edge2D_in
     nod1 = edge_nodes(1,ed)
     nod2 = edge_nodes(2,ed)
     el1  = edge_tri(1,ed)
     el2  = edge_tri(2,ed)

     ! ============
     ! Average diffusive flux
     ! ============


     DO nz=1, nsigma-1
       Kh = 0.5_WP* (Visc(nz,nod1)+Visc(nz,nod2))
       Tx=0.5_WP*Kh*(ttxa(nz,el1)+ttxa(nz,el2))
       Ty=0.5_WP*Kh*(ttya(nz,el1)+ttya(nz,el2))

       ! ============
       ! MUSCL type reconstruction
       ! ============
        if (V_n_filt(nz,el1)*edge_cross_dxdy(1,ed) > U_n_filt(nz,el1)*edge_cross_dxdy(2,ed)) then
           Tmean = ctfold(nz, nod2)
        else
           Tmean = ctfold(nz, nod1)
        end if
        c1= (Je(nz,el1)*V_n_filt(nz,el1)*Tmean - Ty)*edge_cross_dxdy(1,ed)  &
          - (Je(nz,el1)*U_n_filt(nz,el1)*Tmean - Tx)*edge_cross_dxdy(2,ed)

     ! Second segment

           ! ============
           ! Linear upwind reconstruction
           ! ============
        if (V_n_filt(nz,el2)*edge_cross_dxdy(3,ed) < U_n_filt(nz,el2)*edge_cross_dxdy(4,ed) ) then
           Tmean= ctfold(nz, nod2)
        else
           Tmean= ctfold(nz, nod1)
        end if
        c2 = - (Je(nz,el2)*V_n_filt(nz,el2)*Tmean -Ty)*edge_cross_dxdy(3,ed) &
             + (Je(nz,el2)*U_n_filt(nz,el2)*Tmean -Tx)*edge_cross_dxdy(4,ed)

        c_aux = c1 + c2

        ! ============
        ! Average diffusive flux
        ! ============
        ! First segment
        c1_c = V_n_filt(nz,el1)*edge_cross_dxdy(1,ed)  &
             - U_n_filt(nz,el1)*edge_cross_dxdy(2,ed)

        ! Second segment
        c2_c = - V_n_filt(nz,el2)*edge_cross_dxdy(3,ed) &
               + U_n_filt(nz,el2)*edge_cross_dxdy(4,ed)

        aux_c = Je(nz,el1)*c1_c + Je(nz,el2)*c2_c


        trhs_c(nz,nod1) = trhs_c(nz,nod1) + aux_c
        ctrhs(nz,nod1)  = ctrhs(nz,nod1) + c_aux

        trhs_c(nz,nod2) = trhs_c(nz,nod2) - aux_c
        ctrhs(nz,nod2)  = ctrhs(nz,nod2) - c_aux

     END DO
  END DO
!write(*,*) 'trhs_1', maxval(trhs_c),minval(trhs_c), maxval(ctrhs),minval(ctrhs)
  DO ed = edge2D_in+1, edge2D
     nod2= edge_nodes(2,ed)
     nod1= edge_nodes(1,ed)
     el1    = edge_tri(1,ed)

     ! ============
     ! First segment only
     ! ============
     DO nz=1, nsigma-1
     Kh = 0.5_WP* (Visc(nz,nod1)+Visc(nz,nod2))
        ! ============
        ! Average diffusive flux
        ! ============
        Tx = Kh*ttxa(nz,el1)
        Ty = Kh*ttya(nz,el1)

        ! ============
        ! MUSCL type reconstruction
        ! ============
        if(V_n_filt(nz,el1)*edge_cross_dxdy(1,ed) > U_n_filt(nz,el1)*edge_cross_dxdy(2,ed)) then
           Tmean = ctfold(nz, nod2)
        else
           Tmean = ctfold(nz, nod1)
        end if
        c_aux = (Je(nz,el1)*V_n_filt(nz,el1)*Tmean - Ty)*edge_cross_dxdy(1,ed)  &
                         - (Je(nz,el1)*U_n_filt(nz,el1)*Tmean - Tx)*edge_cross_dxdy(2,ed)

        ! ============
        ! Average diffusive flux
        ! ============
        aux_c = Je(nz,el1)*(V_n_filt(nz,el1)*edge_cross_dxdy(1,ed) &
                                       - U_n_filt(nz,el1)*edge_cross_dxdy(2,ed))

     trhs_c(nz,nod1) = trhs_c(nz,nod1) + aux_c
     trhs_c(nz,nod2) = trhs_c(nz,nod2) - aux_c
     ctrhs(nz,nod1) = ctrhs(nz,nod1) + c_aux
     ctrhs(nz,nod2) = ctrhs(nz,nod2) - c_aux

     END DO
  END DO

  ! ===================
  ! Vertical advection
  ! ===================

  cvert(1:nsigma)   = 0.0_WP

  DO n=1, nod2D
     ! ===========
     ! Fluxes in the column
     ! ===========
     ! Surface forcing



     cvert(1)= -Wvel(1,n)*ctfold(1,n)
     cvert(nsigma)=0.0_WP

     ! Extract first special case nz = 2
     if((Wvel(2,n)) > 0.0_WP) then
        dg =  (zbar(2,n) - Z(1,n))*(Z(2,n) - zbar(2,n)) /((Z(3,n) - Z(1,n))*(Z(2,n) - Z(3,n)))
        db = -(zbar(2,n) - Z(1,n))*(Z(3,n) - zbar(2,n)) /((Z(2,n) - Z(1,n))*(Z(2,n) - Z(3,n)))

        cvert(2) = - (ctfold(1,n)*(1.0_WP-dg-db) +ctfold(2,n)*db +ctfold(3,n)*dg)*(Wvel(2,n))

     else
        cvert(2) = - ctfold(1,n) *(Wvel(2,n))

     end if


     DO nz=3, nsigma-2

        ! ============
        ! QUICK upwind (3rd order)
        ! ============
        if((Wvel(nz,n))>= 0.0_WP) then
           !a = zbar(nz,n) - Z(nz-1,n)
           !b = Z(nz,n) - zbar(nz,n)
           !c = Z(nz+1,n) - zbar(nz,n)
           dg =  (zbar(nz,n) - Z(nz-1,n))*(Z(nz,n) - zbar(nz,n)) &
                /((Z(nz+1,n) - Z(nz-1,n))*(Z(nz,n) - Z(nz+1,n)))    !
           db = -(zbar(nz,n) - Z(nz-1,n))*(Z(nz+1,n) - zbar(nz,n)) &
                /((Z(nz,n)   - Z(nz-1,n))*(Z(nz,n) - Z(nz+1,n)))    !

           cvert(nz) = -(ctfold(nz-1,n)*(1.0_WP-dg-db) +ctfold(nz,n)*db +ctfold(nz+1,n)*dg)*(Wvel(nz,n))

        else
           !a = Z(nz,n) - zbar(nz,n)
           !b = zbar(nz,n) - Z(nz-1,n)
           !c = zbar(nz,n) - Z(nz-2,n)
           dg =  (Z(nz,n) - zbar(nz,n))*(zbar(nz,n) - Z(nz-1,n)) &
                /((Z(nz,n) - Z(nz-2,n) )*(Z(nz-2,n) - Z(nz-1,n)))    !  a*b/((c+a)*(b-c))
           db = -(Z(nz,n) - zbar(nz,n))*(zbar(nz,n) - Z(nz-2,n)) &
                /((Z(nz,n) - Z(nz-1,n) )*(Z(nz-2,n) - Z(nz-1,n)))    ! -a*c/((b+a)*(b-c))

           cvert(nz) = -(ctfold(nz,n)*(1.0_WP-dg-db) +ctfold(nz-1,n)*db +ctfold(nz-2,n)*dg)*(Wvel(nz,n))
        end if

        ! cvert(nz)= -Tmean*Wvel(nz,n)! + &
        !	      Kvv*(ctf(nz-1,n)-ctf(nz,n))/(Z(nz-1,n)-Z(nz,n)))
     END DO

     ! And the other special case: nz = nsigma-1
     if((Wvel(nsigma-1,n)) >= 0.0_WP) then
        cvert(nsigma-1) = - ctfold(nsigma-1,n) * (Wvel(nsigma-1,n))
     else
        dg =  (Z(nsigma-1,n) - zbar(nsigma-1,n))*(zbar(nsigma-1,n) - Z(nsigma-2,n)) &
            /((Z(nsigma-1,n) - Z(nsigma-3,n) )  *(Z(   nsigma-3,n) - Z(nsigma-2,n)))
        db = -(Z(nsigma-1,n) - zbar(nsigma-1,n))*(zbar(nsigma-1,n) - Z(nsigma-3,n)) &
            /((Z(nsigma-1,n) - Z(nsigma-2,n) )  *(Z(   nsigma-3,n) - Z(nsigma-2,n)))

        cvert(nsigma-1) = -(ctfold(nsigma-1,n)*(1.0_WP-dg-db) +ctfold(nsigma-2,n)*db &
                          + ctfold(nsigma-3,n)*dg ) * (Wvel(nsigma-1,n))
     end if

     area_inv = 1._WP / area(n)
     DO nz=1,nsigma-1
        ctrhs(nz,n) = ( (ctrhs(nz,n) - ctfold(nz,n)*trhs_c(nz,n))*area_inv   &
                     + ctfold(nz,n)*(Wvel(nz,n) - Wvel(nz+1,n)) + (cvert(nz)-cvert(nz+1))) &
                   / Jc(nz,n)
     END DO
  END DO


  !=================
  ! Relax to climatology if needed
  ! =================
  if (cl_relax_sed .or. sed_boun_flux) then

  do n=1,nod2D
     Vel_nor(n) = 0.0_WP
  enddo

  DO ed = 1+edge2D_in, edge2D
     nod1 = edge_nodes(1,ed)
     nod2 = edge_nodes(2,ed)
     el1  = edge_tri(1,ed)

     dmean = max(Dmin, 0.5_WP*(eta_n(nod1)+depth(nod1) + eta_n(nod2)+depth(nod2)))

     un1 = (V_filt_2D(el1)*edge_cross_dxdy(1,ed)- U_filt_2D(el1)*edge_cross_dxdy(2,ed))*dmean

     Vel_nor(nod1) = Vel_nor(nod1) + un1/area(nod1)
     Vel_nor(nod2) = Vel_nor(nod2) - un1/area(nod2)
  END DO

  if (cl_relax_sed) then

     relax2clim_ac = clim_relax/(24.0_WP*10.0_WP)

     DO n=1, nod2D
        if (index_nod2D(n) == 2) then
           if (Vel_nor(n) > 0.0_WP ) then
             DO nz=1,nsigma-1
                 ctrhs(nz,n) = ctrhs(nz,n) +relax2clim_ac*(Cclim(nz,n) -(ctf(nz,n)+c_aver))
             END DO
          else
              DO nz=1,nsigma-1
                 ctrhs(nz,n) = clim_relax*(Cclim(nz,n) -(ctf(nz,n)+c_aver))
              END DO
           endif
        endif
     END DO

  endif

  if (sed_boun_flux) then
     DO n=1, nobn
        if (Vel_nor(n) > 0.0_WP) then
           DO nz=1,nsigma-1
              ctrhs(nz,n) = ctrhs(nz,n) +relax2sed*(cr_distr2(nz,n)-(ctf(nz,n)+c_aver))
           END DO
        else
           DO nz=1,nsigma-1
             ctrhs(nz,n) = sed_relax*(cr_distr2(nz,n)-(ctf(nz,n)+c_aver))
           END DO
        endif
    END DO
  endif
  endif
  ! =================
  ! Update ctfold (to be used on the next time level)
  ! and compute new ctf
  ! =================

  DO n=1,nod2D
     DO nz=1,nsigma-1
        ctfold(nz,n) = ctf(nz,n) + c_aver
        ctf(nz,n)    = ctf(nz,n) + ctrhs(nz,n)*dt + c_aver
     END DO
  END DO

end subroutine solve_tracer_upwind_sed


subroutine tracer_impl_vert_diff_sed
USE o_MESH
USE o_PARAM
USE o_ARRAYS
!USE g_PARSUP
IMPLICIT NONE

real(kind=WP)              ::  a(nsigma), b(nsigma), c(nsigma), cr(nsigma)
real(kind=WP)              ::  sr(nsigma)
real(kind=WP)              ::  cp(nsigma), tp(nsigma), sp(nsigma), aux,C_d_node(nod2D)
real(kind=WP)              ::  t1, s1, dt_inv, m, sed_source(nsigma-1),frw
integer                    ::  nz, n, i

   a = 0.0_WP
   b = 0.0_WP
   c = 0.0_WP
   cr = 0.0_WP
   sr = 0.0_WP
   cp  = 0.0_WP
   tp = 0.0_WP
   sp = 0.0_WP

if (BFC==1)  then
C_d_node=C_d
else
call compute_el2nodes_2D(C_d_el,C_d_node)
endif

!$OMP PARALLEL PRIVATE(n,nz,a,b,c,cr,sr,cp,tp,sp,t1,s1,dt_inv,m,sed_source,frw)
!$OMP DO
do n=1,nod2D
   if (mask_wd_node(n) /= 0.0_WP) then
   CF(:,n) = CF(:,n) - c_aver
   endif
enddo
!$OMP END DO

dt_inv=1.0_WP/dt
sed_source=0.0_WP
!$OMP DO
DO n=1,nod2D

  if (mask_wd_node(n) /= 0.0_WP) then

   ! Regular part of coefficients:
   DO nz=2, nsigma-2
      a(nz) =-dt*Av_node(nz,n)/((Z(nz-1,n)-Z(nz,n))*(zbar(nz,n)-zbar(nz+1,n)))
      c(nz) =-dt*Av_node(nz,n)/((Z(nz,n)-Z(nz+1,n))*(zbar(nz,n)-zbar(nz+1,n)))
      b(nz) = -a(nz)-c(nz)+1.0_WP   !dt_inv
   END DO
   ! The first row
   c(1) = -dt*Av_node(2,n)/((Z(1,n)-Z(2,n))*(zbar(1,n)-zbar(2,n)))
   a(1) = 0.0_WP
   b(1) = -c(1)+1.0_WP !dt_inv
   ! The last row
   a(nsigma-1) = -dt*Av_node(nsigma-1,n)/((Z(nsigma-2,n)-Z(nsigma-1,n))*(zbar(nsigma-1,n)-zbar(nsigma,n)))
   b(nsigma-1) = -a(nsigma-1)+1.0_WP !dt_inv
   c(nsigma-1) = 0.0_WP
   ! ===========================================
     call calc_sed_source(n,C_d_node(n),sed_source)
      cr(1) = c(1)*(-CF(2,n)+CF(1,n))+dt*sed_source(1)

      DO nz=2,nsigma-2
         cr(nz) = -a(nz)*CF(nz-1,n)-c(nz)*CF(nz+1,n)+(a(nz)+c(nz))*CF(nz,n)+dt*sed_source(nz)
      END DO

      cr(nsigma-1) = a(nz)*(CF(nsigma-1,n)-CF(nsigma-2,n))+dt*sed_source(nsigma-1)


   !     nz=nsigma-1
   !     cr(nz)=-a(nz)*CF(nz-1,n)+a(nz)*CF(nz,n)
   !     sr(nz)=-a(nz)*SF(nz-1,n)+a(nz)*SF(nz,n)

   !  The first row contains also surface forcing
   !
   !   cr(1)=c(1)*(CF(1,n)- CF(2,n)) + dt*(-1.e-6*heat_flux(n)/(4.2_WP)+surf_relax_T*(Tsurf(n)-CF(1,n)))
   !   sr(1)=c(1)*(SF(1,n)-SF(2,n)) ! + &
   !      dt*(SF(1,n)*water_flux(n) +surf_relax_S*(Ssurf(n)-SF(1,n)))

   ! =============================================
   ! The sweep algorithm
   ! initialize c-prime and s,t-prime
   cp(1) = c(1)/b(1)
   tp(1) = cr(1)/b(1)
   ! solve for vectors c-prime and t, s-prime
   do nz = 2,nsigma-1
      m = 1._WP/(b(nz)-cp(nz-1)*a(nz))
      cp(nz) = c(nz)*m
      tp(nz) = (cr(nz)-tp(nz-1)*a(nz))*m
   enddo
   ! initialize x
   cr(nsigma-1) = tp(nsigma-1)
   ! solve for x from the vectors c-prime and d-prime
   do nz = nsigma-2, 1, -1
      cr(nz) = tp(nz)-cp(nz)*cr(nz+1)
   end do

   DO nz=1,nsigma-1
      CF(nz,n) = CF(nz,n)+cr(nz) + c_aver
   END DO
     endif  ! end for node mask
END DO   !!! cycle over nodes
!$OMP END DO
!$OMP END PARALLEL
            !write(*,*) '!!!!!!! ', maxval(SF(:,1:myDim_nod2D)), minval(SF(:,1:myDim_nod2D))


    DO n=1, nod2D
        DO nz=1,nsigma-1
         !  if ((abs( CF(nz,n))< 1.d-12).and.(CF(nz,n)<0.0_WP)) CF(nz,n) = 0.0_WP
           if (CF(nz,n)< 1.e-12) CF(nz,n)=0.0_WP
        enddo
     END DO
!write(*,*) '2sed', sum(CF)


 END subroutine tracer_impl_vert_diff_sed


 SUBROUTINE calc_sed_source(n, C_d_local, sed_source)

    !
    USE o_PARAM
    USE o_MESH
    USE o_ARRAYS

    IMPLICIT NONE
    INTEGER, INTENT(in):: n
    REAL(kind=WP),INTENT(in):: C_d_local
    REAL(kind=WP),INTENT(out):: sed_source(nsigma-1)
    INTEGER  :: j,nz
    REAL(kind=WP):: D, E, density0_inv, tb, u_taub2_aux, dst, b,s, teta,f, tcr,c, rho,co, co_1


rho=rho_c(1,n)+density_0

co=CF(1,n)+c_aver

w_s(1,n)=((1.0_WP-co/plop_s)*(1.0_WP-co/(0.6_WP*plop_s))**1.5_WP)&
*g*(plop_s-rho)*d_m**2/(18.0_WP*rho*snu_kin)

sed_source(1)=-w_s(1,n)*co/Jc(1,n)

DO  nz=2,nsigma-2

  co_1=co

  co=CF(nz,n)+c_aver

  rho=rho_c(nz,n)+density_0

  w_s(nz,n)=((1.0_WP-co/plop_s)*(1.0_WP-co/(0.6_WP*plop_s))**1.5_WP)&
  *g*(plop_s-rho)*d_m**2/(18.0_WP*rho*snu_kin)

  sed_source(nz)=-w_s(nz,n)*co/Jc(nz,n) + w_s(nz-1,n)*co_1/Jc(nz-1,n)

enddo

  co_1=co

  co=CF(nsigma-1,n)+c_aver

  rho=rho_c(nsigma-1,n)+density_0

  w_s(nsigma-1,n)=((1.0_WP-co/plop_s)*(1.0_WP-co/(0.6_WP*plop_s))**1.5_WP)&
  *g*(plop_s-rho)*d_m**2/(18.0_WP*rho*snu_kin)
!-------------------------------------------------------------------------------------
!   compute particle parameter dst
!   Leo and van Rijn "Sediment Transport, Part II:
!   Suspended load Transport". p 1613.
!

   b = 1.0_WP/3.0_WP
   s = plop_s/rho - 1.0_WP
   dst = d_m*(g*s/snu_kin**2)**b


! Criteria is from Miller, Cave and Komar, 1977, Threshold of sediment motion under unidirectional current + Shields criterium (dst>=4)

   teta=0.3_WP/(1+1.2_WP*dst)+0.055_WP*(1-exp(-0.022_WP*dst))




    tcr=teta*d_m*g*(plop_s-rho)

      u_taub2_aux =  C_d_local*(Unode(nsigma-1,n)**2  + Vnode(nsigma-1,n)**2)
      tb=u_taub2_aux*rho
      f=1._WP ! volumetric fraction of the sediment class j, right know we have only one class
       E=Er0*f*(1-e_p)*max((tb/tcr-1),0.0_WP)
     ! Erosion rate  from Vera
     ! Ev= (plop_s-density_0)*(1-e_p)*u_taub-(-Wvel(nsigma-1,n)

     ! Simultaneous erosion and deposition are possible
     D=w_s(nsigma-1,n)*co
     Er_Dep(n)=E-D
     call bottom_sed_disc_calc(n, dst, teta, dsqrt(u_taub2_aux))
     sed_source(nsigma-1)= (E-D)/Jc(nsigma-1,n)+w_s(nsigma-2,n)*co_1/Jc(nsigma-2,n)


END SUBROUTINE calc_sed_source

!=============================================
subroutine bottom_sed_disc_calc(n,dst, teta, u_taub)
  use o_MESH
  use o_ARRAYS
  use o_PARAM
  implicit none

  real(kind=WP) :: s,b, u_st_cr, a,uv, Tsp
  integer, INTENT(in):: n
  real(kind=WP), INTENT(in):: dst, teta, u_taub
! The calculation of the bed load transport according to Van Rijn L.C.
! Sediment transport. Pt.1 Bed load transport, 1984
          s = plop_s/density_0 - 1.0_WP
	  b = d_m*g*s
          u_st_cr = sqrt(b*teta)
          a = u_taub**2/(u_st_cr)**2 - 1.0_WP
	  Tsp = max(0.0_WP,a)

          a = 0.053_WP*dsqrt(g*s)*d_m**1.5_WP/dst**0.3_WP
          uv = dsqrt(Unode(nsigma-1,n)**2 + Vnode(nsigma-1,n)**2)

! compute bed-load transport qb

          if (uv > 0.0_WP) then
          qb(1,n) = a*Tsp**2.1_WP*Unode(nsigma-1,n)/uv
          qb(2,n) = a*Tsp**2.1_WP*Vnode(nsigma-1,n)/uv
          else
          qb(1,n) = 0.0_WP
	  qb(2,n) = 0.0_WP
         endif

end subroutine bottom_sed_disc_calc


subroutine bottom_evolution2

  use o_MESH
  use o_ARRAYS
  use o_PARAM
  implicit none

  real(kind=WP) :: c1, c2, a
  integer       :: ed, el(2), enodes(2),elnodes(4), n, elem
  real(kind=WP), allocatable :: cHrhs(:)

  allocate(cHrhs(nod2D))

  cHrhs = 0.0_WP
 ! ==============
 ! internal edges
 ! ==============

    do elem=1,elem2D
         elnodes=elem2D_nodes(:,elem)
         qbu(elem)=sum(w_cv(1:4,elem)*qb(1,elnodes))
         qbv(elem)=sum(w_cv(1:4,elem)*qb(2,elnodes))
    enddo

 DO ed=1,edge2D_in
    enodes=edge_nodes(:,ed)
    el=edge_tri(:,ed)
    c1=-qbv(el(1))*edge_cross_dxdy(1,ed)+qbu(el(1))*edge_cross_dxdy(2,ed)
    c2=qbv(el(2))*edge_cross_dxdy(3,ed)-qbu(el(2))*edge_cross_dxdy(4,ed)
    c1 = c1 + c2
    cHrhs(enodes(1)) = cHrhs(enodes(1)) + c1
    cHrhs(enodes(2)) = cHrhs(enodes(2)) - c1
 END DO
 ! =============
 ! boundary edges
 ! only the left element (1) is available
 ! =============
  DO ed=1+edge2D_in, edge2D
    enodes=edge_nodes(:,ed)
    el=edge_tri(:,ed)
    c1=-qbv(el(1))*edge_cross_dxdy(1,ed)+qbu(el(1))*edge_cross_dxdy(2,ed)
    cHrhs(enodes(1)) = cHrhs(enodes(1)) + c1
    cHrhs(enodes(2)) = cHrhs(enodes(2)) - c1
 END DO

   a = dt_2D/(1.0_WP - e_p)
   do n=1,nod2D
   h_var2(n)=h_var2(n) + ac(n)*a*(cHrhs(n)/area(n) - Er_Dep(n)/plop_s)
   if (h_var2(n)>15) then
    write(*,*) 'h_crap:' , cHrhs(n), Er_Dep(n),h_var2(n),n, qb(:,n),Unode(nsigma-1,n),Vnode(nsigma-1,n)
   stop
   endif
  ENDDO

  deallocate(cHrhs)
  end subroutine bottom_evolution2
