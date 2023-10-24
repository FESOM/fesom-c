 SUBROUTINE compute_ssh_rhs_elem

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM
  USE i_PARAM
  USE i_ARRAYS
  use g_parsup
  use g_comm_auto

  IMPLICIT NONE

  ! -div(Hu) is computed in cycle over edges
  !
  ! The depth is estimated at elements
  ! c1 is -u_n*L_left*d1, c2  -u_n*L_right*d1

  integer        :: ed, el(2), elem, j, n, q, i
  real(kind=WP)  :: c1(myDim_edge2D+eDim_edge2D)
  real(Kind=WP)  :: aux_elem(myDim_elem2D+eDim_elem2D+eXDim_elem2D)
  !SH Test
  real(Kind=WP)  :: ttldep, aux_depAB(4,myDim_elem2D)
  integer :: edglim

  real(Kind=WP)  ::  wdx, wdy, wndm, rtmp, A_el

  do i=1,myDim_edge2D+eDim_edge2D
     C1(i)=0.0_WP
  end do

  if (use_ice) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(elem,wdx,wdy,wndm,rtmp,A_el)
!$OMP DO
    do elem=1,myDim_elem2D
        A_el =   w_cv(1,elem) *a_ice(elem2D_nodes(1,elem)) + w_cv(2,elem) *a_ice(elem2D_nodes(2,elem)) &
            + w_cv(3,elem) *a_ice(elem2D_nodes(3,elem)) + w_cv(4,elem) *a_ice(elem2D_nodes(4,elem))
        if (A_el > 0.0) then
            wdx = U_n_ice(1,elem) - U_n(1,elem) ! icedrift - ocean current ( x direction)
            wdy = U_n_ice(2,elem) - V_n(1,elem) ! icedrift - ocean current ( y direction)
            wndm = SQRT( wdx * wdx + wdy * wdy )
            rtmp = wndm * Cd_oce_ice   ! from fesom2, Cd=5.5*10-3 , https://doi.org/10.1175/JPO-D-18-0185.1
            taux(elem) = rtmp * wdx  ! new taux (stress along x)
            tauy(elem) = rtmp * wdy  ! new tauy (stress along y)
        endif
    enddo
#ifdef USE_MPI
  call exchange_elem(taux)
  call exchange_elem(tauy)
#endif
!$OMP END DO
!$OMP END PARALLEL
  endif

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n,ed,el,elem)
!$OMP DO
  DO n=1,myDim_nod2D+eDim_nod2D
     ssh_rhs(n)=0.0_WP
  ENDDO
!$OMP END DO NOWAIT

  ! ==============
  ! fill mean depth
  ! ==============

  !SH Test intermediate step
  DO elem=1,myDim_elem2D
     do i=1,4
        aux_depAB(i,elem)=Dmin
        ttldep=etaAB(elem2D_nodes(i,elem)) + depth(elem2D_nodes(i,elem))
        if (ttldep>Dmin) aux_depAB(i,elem)=ttldep
     end do
  END DO

!$OMP DO
  DO elem=1,myDim_elem2D
     !NR    elnodes=elem2D_nodes(:,elem)
     !NR    fD=max(Dmin,depth(elnodes) + etaAB(elnodes))
     !NR    aux_elem(q)=sum(w_cv(1:4,elem)*fD)
!!$     aux_elem(elem) = w_cv(1,elem) *max(Dmin, etaAB(elem2D_nodes(1,elem)) + depth(elem2D_nodes(1,elem))) &
!!$          + w_cv(2,elem) *max(Dmin, etaAB(elem2D_nodes(2,elem)) + depth(elem2D_nodes(2,elem))) &
!!$          + w_cv(3,elem) *max(Dmin, etaAB(elem2D_nodes(3,elem)) + depth(elem2D_nodes(3,elem))) &
!!$          + w_cv(4,elem) *max(Dmin, etaAB(elem2D_nodes(4,elem)) + depth(elem2D_nodes(4,elem)))

     aux_elem(elem) = w_cv(1,elem)*aux_depAB(1,elem) +  w_cv(2,elem)*aux_depAB(2,elem) + &
                      w_cv(3,elem)*aux_depAB(3,elem) +  w_cv(4,elem)*aux_depAB(4,elem)
  END DO
!$OMP END DO

#ifdef USE_MPI
  call exchange_elem(aux_elem)
#endif

!SH print *,'AUX_ELEM:',mype,minval(aux_elem),maxval(aux_elem)


  ! ==============
  ! internal edges
  ! ==============

#ifdef USE_MPI
  edglim=myDim_edge2D+eDim_edge2D
#else
  edglim=edge2D_in
#endif

!$OMP DO
  DO ed=1,edglim

#ifdef USE_MPI
     if (myList_edge2D(ed)>edge2D_in) cycle
#endif

     el     = edge_tri(:,ed)

#ifdef DEBUG
        if (edge_tri(1,ed)<=0 .or. edge_tri(2,ed)<=0) print *,'WARNING: edge_tri Problem!',mype,ed
#endif

     c1(ed) = ( UAB(2,el(1))*edge_cross_dxdy(1,ed) - UAB(1,el(1))*edge_cross_dxdy(2,ed)) *aux_elem(el(1))  &
            + (-UAB(2,el(2))*edge_cross_dxdy(3,ed) + UAB(1,el(2))*edge_cross_dxdy(4,ed)) *aux_elem(el(2))

  END DO

!$OMP END DO NOWAIT

  ! =============
  ! boundary edges
  ! only the left element (1) is available
  ! =============

#ifdef USE_MPI

  DO ed=1,edglim
     if (myList_edge2D(ed)>edge2D_in) then

#ifdef DEBUG
        if (edge_tri(1,ed)<=0 .or. edge_tri(2,ed)>0) print *,'WARNING: edge_tri Problem!',mype,ed
#endif

        elem=edge_tri(1,ed)
        c1(ed) = (  UAB(2,elem)*edge_cross_dxdy(1,ed) &
                  - UAB(1,elem)*edge_cross_dxdy(2,ed) ) * aux_elem(edge_tri(1,ed))

 !if (abs(c1(ed))>1.e-20) print *,'HUAA', myList_edge2D(ed),c1(ed)
     end if
  END DO

#else

!$OMP DO
  DO ed=1+edge2D_in, edge2D
     c1(ed) = ( UAB(2,edge_tri(1,ed))*edge_cross_dxdy(1,ed) &
          -UAB(1,edge_tri(1,ed))*edge_cross_dxdy(2,ed) ) * aux_elem(edge_tri(1,ed))
!if (abs(c1(ed))>1.e-20) print *,'HUAA', ed,c1(ed)
  END DO
!$OMP END DO

#endif


!!SH print *,'EDGES TEST',mype,minval(c1),maxval(c1)

!$OMP END PARALLEL

!NROMP  To be parallelized:
  DO ed=1,myDim_edge2D+eDim_edge2D
     ssh_rhs(edge_nodes(1,ed)) = ssh_rhs(edge_nodes(1,ed))+c1(ed)
     ssh_rhs(edge_nodes(2,ed)) = ssh_rhs(edge_nodes(2,ed))-c1(ed)
  END DO

!!$#ifdef USE_MPI
!!$  call exchange_nod(ssh_rhs)
!!$#endif

  !increase SSH level due to rivers (based on river runoff)
  ! new rivers IK
  if (key_rivers) then
     !runoff_rivers m3/s (interpolated in time runoff, see rivers_do, fv_rivers.f90)
     !my_Nnods_rivers number of nodes with rivers
     !my_nod_rivers array of nodes indexes
     do n = 1, my_Nnods_rivers
        ssh_rhs(my_nod_rivers(n)) = ssh_rhs(my_nod_rivers(n))+runoff_rivers(n)
     end do
  endif

  !    open boundary nodes
  ! newly done SH, IK
  !NROMP   Build index array for open boundary nodes, to simplify and parallelize this loop.
  !VF, TF_presence....., ssh_rhs
  if (TF_presence .and. mynobn>0)  then

     ssh_rhs(my_in_obn)=0.0_WP

!!$OMP PARALLEL PRIVATE(n,q)
!!$OMP DO
     do i=1,mynobn
        n=my_in_obn(i)
        !SH q=i
        q=my_in_obn_idx(i)
!!$OMP REDUCTION(+:ssh_rhs(n))
        do j=1,12
           ssh_rhs(n) = ssh_rhs(n)+amp_factor(j)*ampt(q,j)&
              *cos((time_2D-time_jd0*86400.0_WP)*2.0_WP*pi/(Harmonics_tf_period(j)*3600._WP) - (fazt(q,j)+phy_factor(j)))
        end do
        if (.not. restart) then
            if (time_jd-time_jd0<10.0_WP) then
            ssh_rhs(n)=ssh_rhs(n)*(time_jd-time_jd0)/10.0_WP
            end if
        endif
     enddo
!!$OMP END DO
!!$OMP END PARALLEL
!if (mype==2) then
!   write(*,*) ssh_rhs(my_in_obn(43)),ssh_rhs(my_in_obn(44)),ssh_rhs(my_in_obn(44))-ssh_rhs(my_in_obn(43))
!endif
  endif

!SH Check if the exchange is necessary!
!!$#ifdef USE_MPI
!!$  call exchange_nod(ssh_rhs)
!!$#endif

!!SHDB print *,'compute_ssh_rhs_elem: ssh_rhs:',minval(ssh_rhs(:)),maxval(ssh_rhs(:))


END subroutine compute_ssh_rhs_elem
!===========================================================================
SUBROUTINE update_2D_vel(step)

  ! It is intended for AB3

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  use g_parsup
  use g_comm_auto

  IMPLICIT NONE

  integer, intent(in)   :: step
  integer               :: el, elem, elnodes(4), n, rie
  real(kind=WP)    :: dmean, dmean2, ff, friction
  real(kind=WP)    :: InBr           ! inverse barometer
  real(kind=WP)    :: fD(4), fDD(4), eta(myDim_nod2D+eDim_nod2D)
  real(kind=WP), parameter ::  de=0.614_WP
  real(kind=WP), parameter ::  ga=0.088_WP
  real(kind=WP), parameter ::  ep=0.013_WP
  real(kind=WP), parameter ::  fe=1.0_WP-de-ga-ep
  real(kind=WP)    :: a_Cd, density0_inv, C_d2, C_d3, sqfv


  call vel_gradients_2D(UAB)   !AA ! We need them only for biharmonic viscosity    !CHECKSH
                               !NROMP and for something else? It may not be skipped,
                               !NROMP even if filtering viscosity is used.
!!$print *,'vel_grad1',step,minval(vel_grad(1,:)),maxval(vel_grad(1,:))
!!$print *,'vel_grad2',step,minval(vel_grad(2,:)),maxval(vel_grad(2,:))
!!$print *,'vel_grad3',step,minval(vel_grad(3,:)),maxval(vel_grad(3,:))
!!$print *,'vel_grad4',step,minval(vel_grad(4,:)),maxval(vel_grad(4,:))
  C_d_el=C_d_2d_e
  !call bottom_friction_C_d_el(step) !SH not parallel only case BFC==1 works  !CHECKSH
   C_d_el = C_d
  ! ====================
  ! Sea level contribution   -D\nabla\eta
  ! and  the Coriolis force (elemental part)
  ! ====================
  ! Below are the interpolation coefficients for the elevation
  ! as recommended in Shchepetkin and McWilliams (2005).
  ! CFL stabily limit is essentially improved with them.
  ! Interpolate AM4

  density0_inv = 1._WP/(density_0)

  !if (T_potential) call potential           ! AA geopotential now computed in:   oce_timestep


!$OMP PARALLEL PRIVATE(n,density0_inv,el,dmean)
  if (key_atm) then
!$OMP DO
     DO n=1,myDim_nod2D+eDim_nod2D
        eta_n_2(n) = de*ssh_rhs(n) + fe*eta_n(n) + ga*eta_n_1(n) + ep*eta_n_2(n)
        etaAB(n)   = eta_n_2(n)
        eta(n)     = g*(a_tp(n,1)*eta_n_2(n) - emp(n)*density0_inv - a_tp(n,2)*ssh_gp(n)) &
             + (P_ref - mslp(n))*density0_inv ! Use AM4 interpolation
     END DO
!$OMP END DO
  else ! key_atm=False
     !$OMP DO

     DO n=1,myDim_nod2D+eDim_nod2D
        eta_n_2(n) = de*ssh_rhs(n) + fe*eta_n(n) + ga*eta_n_1(n) + ep*eta_n_2(n)
        etaAB(n)   = eta_n_2(n)
        !      eta(n)     = g*(a_tp(n,1)*eta_n_2(n) - emp(n)*density0_inv - a_tp(n,2)*ssh_gp(n))
        eta(n)     = g*(a_tp(n,1)*eta_n_2(n) - a_tp(n,2)*ssh_gp(n))

     END DO
!$OMP END DO
  endif

  !++++++++++++++++++++++++++++++++++
  ! compute mask for wetting/drying
  !++++++++++++++++++++++++++++++++++

if (WET_DRY_ON) call wad_mask

!$OMP DO
  DO el=1,myDim_elem2D
     dmean =  max(Dmin, w_cv(1,el)*(etaAB(elem2D_nodes(1,el)) + depth(elem2D_nodes(1,el))) &
                     +  w_cv(2,el)*(etaAB(elem2D_nodes(2,el)) + depth(elem2D_nodes(2,el))) &
                     +  w_cv(3,el)*(etaAB(elem2D_nodes(3,el)) + depth(elem2D_nodes(3,el))) &
                     +  w_cv(4,el)*(etaAB(elem2D_nodes(4,el)) + depth(elem2D_nodes(4,el)))  )

     U_rhs_2D(1,el) = elem_area(el) * (-dmean*( gradient_sca(1,el)*eta(elem2D_nodes(1,el)) &
          + gradient_sca(2,el)*eta(elem2D_nodes(2,el)) &
          + gradient_sca(3,el)*eta(elem2D_nodes(3,el)) &
          + gradient_sca(4,el)*eta(elem2D_nodes(4,el)) ) &
          + UAB(2,el)*dmean*coriolis(el)                                                 &
          + taux(el)*density0_inv  )

     U_rhs_2D(2,el) = elem_area(el) * (-dmean*( gradient_sca(5,el)*eta(elem2D_nodes(1,el)) &
          + gradient_sca(6,el)*eta(elem2D_nodes(2,el)) &
          + gradient_sca(7,el)*eta(elem2D_nodes(3,el)) &
          + gradient_sca(8,el)*eta(elem2D_nodes(4,el)) ) &
          - UAB(1,el)*dmean*coriolis(el)                                                 &
          + tauy(el)*density0_inv  )

  END DO
!$OMP END DO
!$OMP END PARALLEL

#ifdef USE_MPI
  call exchange_elem(U_rhs_2D) !SH necessary??
#endif


!!$  print *,'U_rhs_2D : ',mype,minval(U_rhs_2D(1,:)),maxval(U_rhs_2D(1,:))
!!$  print *,'V_rhs_2D : ',mype,minval(U_rhs_2D(2,:)),maxval(U_rhs_2D(2,:))

  ! ======================
  ! Momentum advection   -\int div(uDu)dS=-\sum uD(un)l
  ! and viscosity: (assembly over edges)
  ! ======================
!SH SKIPPED FOR NOW  if (mom_adv_2D == 1) call momentum_adv_scalar_2D
  if (mom_adv_2D == 2) call momentum_adv_upwind_2D

  If (filt_2D) call viscosity_filt_2D
  If (filt_bi_2D) call viscosity_filt2x_2D
  IF (bih_2D) call biharmonic_2D

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(elem,ff,dmean, friction,a_Cd,elnodes)
  !AA connection with 3D part of the model. If we compute only barotropic task 2D or 3D this is not needed.

  if (type_task>2) then
!$OMP DO SCHEDULE(STATIC)
     do elem=1,myDim_elem2D + eDim_elem2D+eXDim_elem2D
        U_rhs_2D(1,elem) = U_rhs_2D(1,elem) + U_rhs_2D_3D(1,elem)
        U_rhs_2D(2,elem) = U_rhs_2D(2,elem) + U_rhs_2D_3D(2,elem)
     enddo
!$OMP END DO NOWAIT
  endif

!!$#ifdef USE_MPI
!!$  call exchange_elem(U_rhs_2D) !SH necessary??
!!$#endif

  ! ==============
  ! Final update and bottom drag
  ! The contribution from the bottom drag is taken implicitly.
  ! ==============
  !VF, task, friction
!$OMP DO SCHEDULE(STATIC)

  DO elem=1,myDim_elem2D
     elnodes = elem2D_nodes(:,elem)

     dmean = max(Dmin,sum(w_cv(1:4,elem)*(eta_n(elnodes) + depth(elnodes))))
     a_Cd = C_d_el(elem)*((1.0_WP - sum(w_cv(1:4,elem)*ac(elnodes)))*Cd_bz + 1.0_WP)
     if (type_task == 1) then

        friction = a_Cd*sqrt(U_n_2D(1,elem)**2+U_n_2D(2,elem)**2)


        if (use_wef) then
            ff = 1._WP/(max(Dmin,sum(w_cv(1:4,elem)*(ssh_rhs(elnodes) + depth(elnodes)))) + friction*dt_2D + &
                        dmean*sum(w_cv(1:4,elem)*wef(elnodes))*dt_2D)
        else
            ff = 1._WP/(max(Dmin,sum(w_cv(1:4,elem)*(ssh_rhs(elnodes) + depth(elnodes)))) + friction*dt_2D)
        end if

        U_rhs_2D(1,elem) = mask_wd(elem)*((dmean*U_n_2D(1,elem) &
             + dt_2d*Bar_pr_2D(1,elem) +&
             dt_2D*U_rhs_2D(1,elem)/elem_area(elem))*ff + &
             UV2_rhs(1,elem))

        U_rhs_2D(2,elem) = mask_wd(elem)*((dmean*U_n_2D(2,elem) &
             + dt_2d*Bar_pr_2D(2,elem) +&
             dt_2D*U_rhs_2D(2,elem)/elem_area(elem))*ff + &
             UV2_rhs(2,elem))


     else

        friction = a_Cd*sqrt(U_n(nsigma-1,elem)**2+V_n(nsigma-1,elem)**2)

        ff = 1._WP/(max(Dmin,sum(w_cv(1:4,elem)*(ssh_rhs(elnodes) + depth(elnodes)))) + friction*dt_2D)


        U_rhs_2D(1,elem) = mask_wd(elem)*((dmean*U_n_2D(1,elem) - &
             dt_2D*friction*(U_n(nsigma-1,elem) - U_n_2D(1,elem)) + dt_2d*Bar_pr_2D(1,elem) +&
             dt_2D*U_rhs_2D(1,elem)/elem_area(elem))*ff + &
             UV2_rhs(1,elem))
        U_rhs_2D(2,elem) = mask_wd(elem)*((dmean*U_n_2D(2,elem) - &
             dt_2D*friction*(V_n(nsigma-1,elem) - U_n_2D(2,elem)) + dt_2d*Bar_pr_2D(2,elem) +&
             dt_2D*U_rhs_2D(2,elem)/elem_area(elem))*ff + &
             UV2_rhs(2,elem))

     endif

     ! Now it contains the updated velocity

  END DO

!if (mype==0) print *,'HUGH',step,minval(V_n(nsigma-1,:)),maxval(V_n(nsigma-1,:))
!if (mype==0) print *,'HUGH',step,V_n(nsigma-1,4767),V_n(nsigma-1,4434),V_n(nsigma-1,11595),V_n(nsigma-1,21928)

!$OMP END DO
!$OMP END PARALLEL

#ifdef USE_MPI
  call exchange_elem(U_rhs_2D)
#endif

!!$  print *,'Bar_pr1: ',mype,minval(Bar_pr_2D(1,:)),maxval(Bar_pr_2D(1,:))
!!$  print *,'Bar_pr2: ',mype,minval(Bar_pr_2D(2,:)),maxval(Bar_pr_2D(2,:))
!!$
!!$  print *,'U_n 2: ',mype,minval(U_n(1,:)),maxval(U_n(1,:))
!!$  print *,'V_n 2: ',mype,minval(U_n(2,:)),maxval(U_n(2,:))

!SHDB  print *,'U_rhs_2D 2: ',mype,minval(U_rhs_2D(1,:)),maxval(U_rhs_2D(1,:))
!SHDB  print *,'V_rhs_2D 2: ',mype,minval(U_rhs_2D(2,:)),maxval(U_rhs_2D(2,:))


END SUBROUTINE update_2D_vel
!==========================================================================
subroutine compute_vortex_2D
  use o_MESH
  use o_ARRAYS
  use o_PARAM

  implicit none
  integer               :: ed, el(2), enodes(2) , n
  real(kind=WP)     :: deltaX1, deltaX2, deltaY1, deltaY2, c1

 vorticity_2D=0.0_WP

do ed=1,edge2D
   enodes=edge_nodes(:,ed)
   el=edge_tri(:,ed)
   deltaX1=edge_cross_dxdy(1,ed)
   deltaY1=edge_cross_dxdy(2,ed)
   c1=deltaX1*U_n_2D(1,el(1))+deltaY1*U_n_2D(2,el(1))
   if(el(2)>0) then
    deltaX2=edge_cross_dxdy(3,ed)
    deltaY2=edge_cross_dxdy(4,ed)
    c1=c1-deltaX2*U_n_2D(1,el(2))-deltaY2*U_n_2D(2,el(2))
   endif
   vorticity_2D(enodes(1))=vorticity_2D(enodes(1))+c1
   vorticity_2D(enodes(2))=vorticity_2D(enodes(2))-c1
end do
do n=1,nod2D
 vorticity_2D(n)=vorticity_2D(n)/area(n)
enddo

end subroutine compute_vortex_2D
!===============================================================================================
subroutine energy(eout)

  use o_MESH
  use o_PARAM
  use o_ARRAYS

  use g_parsup
  use g_comm_auto

  implicit none

  integer          :: n, el, elnodes(4)
  real(kind=WP)  :: eout, eout_sum, fD(4), my_area, ttl_area
  real(kind=WP), allocatable :: rbuffer(:), ekin(:)
  integer :: ierror

  allocate(rbuffer(elem2D), ekin(myDim_elem2D))

  !IK print *,'ENERGY: ',minval(eta_n),maxval(eta_n)

  eout=0.0_WP
  my_area=0.0_WP

  ! potential energy
  do n=1,myDim_nod2D
!     if (abs(eta_n(n))<1.e-20) eta_n(n)=0.0_WP  !SH Workaround
     eout=eout + 0.5_WP*g*area(n)*eta_n(n)**2
     my_area=my_area+area(n)
  end do

!IK print *,'ENERGY_U: ',minval(U_n_2D(1,:)),maxval(U_n_2D(1,:))
!IK print *,'ENERGY_V: ',minval(U_n_2D(2,:)),maxval(U_n_2D(2,:))

  ! kinetic energy
  do el=1,myDim_elem2D
     elnodes=elem2D_nodes(:,el)
     !if (elnodes(1)<=myDim_nod2D .and. elnodes(2)<=myDim_nod2D .and. &
     !    elnodes(3)<=myDim_nod2D .and. elnodes(4)<=myDim_nod2D ) then
!     if (U_n_2D(1,el)<1.e-20)  U_n_2D(1,el)=0.0_WP !SH Workaround
!     if (U_n_2D(2,el)<1.0e-20) U_n_2D(2,el)=0.0_WP !SH Workaround
        fD=max(Dmin,depth(elnodes) + eta_n(elnodes))
        ekin(el)=0.5_WP*sum(w_cv(1:4,el)*fD)*elem_area(el)*&
             (U_n_2D(1,el)**2+U_n_2D(2,el)**2)

        !eout=eout+0.5_WP*sum(w_cv(1:4,el)*fD)*elem_area(el)*&
        !          (U_n_2D(1,el)**2+U_n_2D(2,el)**2)
        !my_area=my_area+elem_area(el)
     !end if
  end do

#ifdef USE_MPI
  call gather_elem(ekin,rbuffer)

  call MPI_AllREDUCE(eout, eout_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
          MPI_COMM_FESOM_C, MPIerr)

  if (mype==0) then
  do el=1,elem2D
     eout_sum=eout_sum+rbuffer(el)
  end do
  endif
  call MPI_BCast(eout_sum, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)

  call MPI_AllREDUCE(my_area, ttl_area, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
          MPI_COMM_FESOM_C, MPIerr)
#else
  eout_sum=eout+SUM(ekin)
  ttl_area=my_area
#endif

!IK print *,'mype, ttl_area= ',mype, my_area, ttl_area

!VF, per unit area to whole area

   eout=eout_sum/ttl_area    ! Energy for area

deallocate(rbuffer, ekin)

end subroutine energy

!=============================================================================================
subroutine C_d_el_calc_2d(elem)
  use o_MESH
  use o_PARAM
  use o_ARRAYS
!VF: please, look on Bi and Toorman, 2015 "Mixed sediment transport modelling in Scheldt estuary
!with a physics-based botttom friction law"
  implicit none

  integer,INTENT(IN)            :: elem
  integer                       :: elnodes(4),n
  real(kind=WP)  :: dmean, vel,sqfv, con_el, h_p, A_p, phi, aux, fA, C_ref, h_ref,snu_aux

elnodes = elem2D_nodes(:,elem)
dmean=max(Dmin,sum(w_cv(1:4,elem)*(eta_n(elnodes) + depth(elnodes))))
vel=dsqrt(U_n_2D(1,elem)**2+U_n_2D(2,elem)**2)

If (comp_sediment) then
con_el=sum(w_cv(1:4,elem)*con(elnodes))*plop_s*dmean
else
con_el=0.0_WP
endif
phi=con_el/plop_s !concentr./plop_s
 !version 1: Suspension viscosity calculation, C_ref, h_ref can be calibrated!!!
 !C_ref=0.222_WP  !g/l (kg/m**3)
 !h_ref=0.12_WP   !m
 !snu_aux=snu_kin*(1+con_el*dmean/(C_ref*h_ref))
 !version 2: approach by Fei and Xiangjun, 1982
snu_aux=snu_kin*(1.0_WP-1.35_WP*phi)**(-2.5_WP)
!Non-dimensionalised water depth
 h_p=dsqrt(3*vel*dmean/snu_aux)

!Turbulence dampig factor,A_+(A_p) is an empirical value
A_p=17.0_WP
fA=(1-exp(-h_p/A_p))**2

!Suspension friction, empirical value
beta=0.045_WP

! Volumetric suspended particle concentration

aux=(z0b_min+beta*phi*dmean)*1.5_WP*vel/dmean
! Calculation of the squared friction velocity
sqfv=fA*(cka*vel/(log(dmean/z0b_min)-1.0_WP+z0b_min/dmean))**2 &
+ (dsqrt(aux**2 + 3.0_WP*snu_aux*vel/dmean)+aux)**2

 C_d_el(elem)=min(sqfv/vel**2,0.05)

end subroutine C_d_el_calc_2d

!=============================================================================================
subroutine C_d_el_calc_3d(elem)
  use o_MESH
  use o_PARAM
  use o_ARRAYS
!VF: please, look on Bi and Toorman, 2015 "Mixed sediment transport modelling in Scheldt estuary
!with a physics-based botttom friction law"
  implicit none

  integer,INTENT(IN)            :: elem
  integer                       :: elnodes(4),j
  real(kind=WP)  :: sqfv, vel, con_el, phi, snu_aux,nu, z0, u_taub
  real(kind=WP)  :: v_sqrt, z0b_gotm, rr



! elnodes = elem2D_nodes(:,elem)

 !vel=dsqrt(U_n(nsigma-1,elem)**2+V_n(nsigma-1,elem)**2)

!If (comp_sediment) then
!con_el=sum(w_cv(1:4,elem)*CF(nsigma-1,elnodes))
!else
!con_el=0.0_WP
!endif
!nu=sum(w_cv(1:4,elem)*Z(nsigma-1,elnodes))/dmean
!nu=1-0.5_WP*(sigma(nsigma-1)-sigma(nsigma))
!phi=con_el/plop_s
!snu_aux=snu_kin*(1.0_WP-1.35_WP*phi)**(-2.5_WP)
!sqfv=vel*(Av(nsigma-1,elem)+snu_aux)/((1.0_WP-nu)*0.5_WP*(sigma(nsigma-1)-sigma(nsigma))*dmean)
!C_d_el(elem)= sqfv/(vel**2)
!endif

 nu= max(Je(nsigma-1,elem),dmin)*0.5_WP

 if (comp_sediment) then
    za=2.85_WP*d_m/30.0_WP
 else
    za=0.0_WP
 endif
   !  Iterate bottom roughness length 10 times, see Burchard, 2001; Kagan, 1995

! vel = U_n(nsigma-1,elem)**2+V_n(nsigma-1,elem)**2

! if (vel > 0.0_WP) then
!    v_sqrt = sqrt(vel)*cka
!    z0b_gotm=z0b_min+za
!    u_taub = v_sqrt / log(1._WP + nu/z0b_gotm)

!     DO j=1,2
!     z0b_gotm = 0.1_WP*min(snul/u_taub,1.0_WP)+z0b_min+za
!     u_taub = v_sqrt / log(1._WP + nu/(z0b_gotm))
!     END DO
!     z0b_gotm_el(elem)=z0b_gotm!0.1*snu_kin/max(snu_kin,u_taub)+0.03_WP*(z0b_min+za)
!     C_d_el(elem) = min(u_taub**2/vel,0.05_WP)
!  else
     z0b_gotm_el(elem)=z0b_min+za
     C_d_el(elem) = (cka/(log(1._WP + nu/(z0b_min+za))))**2
! endif

end subroutine C_d_el_calc_3d

subroutine bottom_friction_C_d_el(step)
  use o_MESH
  use o_PARAM
  use o_ARRAYS
!VF: calculation of C_d_el if it is needed (C_d is not fixed constant)
  implicit none
  integer, intent(in)           :: step
  integer                       :: elem, elnodes(4)
  real(kind=WP)                 :: dmean, Cd_crit_depth, Cd_crit_depth_coef


  if (BFC==1) then
     if (use_Cd_2d) then
        C_d_el=C_d_2d_e
     else
        C_d_el=C_d
     end if


  elseif (BFC==3) then
    !IK : increase Cd in the shallow area
     C_d_el = C_d
     Cd_crit_depth      = 3.0  ! critical depth, below this depth Cd will be increased
     Cd_crit_depth_coef = 10.0  !  multiplication coefficient for Cd at depth = 0.0
!$OMP PARALLEL DO private(elem,elnodes,dmean)
     do elem = 1, elem2D
        elnodes = elem2D_nodes(:,elem)
        dmean = sum(w_cv(1:4,elem)*(eta_n(elnodes) + depth(elnodes)))
        if ( dmean <= Cd_crit_depth ) then
!            C_d_el(elem) = C_d*Cd_crit_depth_coef
            C_d_el(elem) = C_d*( dmean * (1-Cd_crit_depth_coef)/Cd_crit_depth+Cd_crit_depth_coef )
        end if
     enddo
!$OMP END PARALLEL DO

  elseif (BFC==2) then
     if (type_task==1) then
!$OMP PARALLEL DO private(elem,elnodes,dmean)
        do elem=1,elem2D
           if ((U_n_2D(1,elem)**2+U_n_2D(2,elem)**2)>0.0_WP) then
              call C_d_el_calc_2d(elem)
           else
              elnodes = elem2D_nodes(:,elem)
              dmean = max(Dmin,sum(w_cv(1:4,elem)*(eta_n(elnodes) + depth(elnodes))))
              C_d_el(elem) = (cka/log(dmean/z0b_min))**2
           endif
        enddo
!$OMP END PARALLEL DO
     elseif (step==1) then !NR As far as I can tell, the result does not change
                           !NR during the barotropic time stepping (result is based on
                           !NR 3D-arrays U_n, V_n)  -> Calculate only once in the first step
!$OMP PARALLEL DO private(elem)
        do elem=1,elem2D
           call C_d_el_calc_3d(elem)
        enddo
!$OMP END PARALLEL DO
     endif
  endif

end subroutine bottom_friction_C_d_el
