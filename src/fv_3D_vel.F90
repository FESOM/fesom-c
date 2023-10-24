SUBROUTINE compute_vel_rhs

  !a modification by Alexey Androsov
  !a (for sigma coordinates)
  !a 09.10.14

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  USE g_PARSUP
  use g_comm_auto

  IMPLICIT NONE

  integer          :: el, elnodes(4), nz, rie , elem
  real(kind=WP)    :: eta(4), ff, gg, mm , pre(4)
  real(kind=WP)    :: Fx, Fy, area_inv
  real(kind=WP)    :: Atm_mean, InBr   ! inverce barometer
  !logical, save    :: lfirst=.true.
  real(kind=WP)    :: selnodes, density0_inv
  real(kind=WP)    :: emp_mean         ! evaporation - precipitation on element

  density0_inv = 1._WP/density_0

  ! ====================
  ! Sea level and pressure contribution   -\nabla(\eta +hpressure/rho_0)
  ! and the Coriolis force + metric terms
  ! ====================
  !write(*,*) 'Urhs1', maxval(U_rhs)
  !VF, change PRIVATE(el,pre,ff...
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(el,pre,ff,Fx,Fy,InBr,Atm_mean)
!$OMP DO

  DO el=1,myDim_elem2D

     pre = -g*(a_tp(elem2D_nodes(:,el),1)*eta_n(elem2D_nodes(:,el)) - &
          a_tp(elem2D_nodes(:,el),2)*ssh_gp(elem2D_nodes(:,el)))

     if (key_atm) then
        Atm_mean = w_cv(1,el)*mslp(elem2D_nodes(1,el)) &
             +  w_cv(2,el)*mslp(elem2D_nodes(2,el)) &
             +  w_cv(3,el)*mslp(elem2D_nodes(3,el)) &
             +  w_cv(4,el)*mslp(elem2D_nodes(4,el))

        ! for emp used only density_0 , but not a density on surface : rho_c(1,n)
        emp_mean = w_cv(1,el)*emp(elem2D_nodes(1,el)) &
             +  w_cv(2,el)*emp(elem2D_nodes(2,el)) &
             +  w_cv(3,el)*emp(elem2D_nodes(3,el)) &
             +  w_cv(4,el)*emp(elem2D_nodes(4,el))

        pre = pre - (P_ref - Atm_mean)*density0_inv + g*emp_mean*density0_inv
     endif


     Fx = sum(gradient_sca(1:4,el)*pre)*elem_area(el)
     Fy = sum(gradient_sca(5:8,el)*pre)*elem_area(el)

     ! =================
     ! Take care of the AB part
     ! =================
     !NR U_rhs(1:nsigma-1,el) = -(0.5_WP+epsilon)*U_rhsAB(1:nsigma-1,el)
     !NR V_rhs(1:nsigma-1,el) = -(0.5_WP+epsilon)*V_rhsAB(1:nsigma-1,el)

     U_rhs(1:nsigma-1,el) =  -(0.5_WP+epsilon)*U_rhsAB(1:nsigma-1,el) &
          + Fx*Je(1:nsigma-1,el)

     V_rhs(1:nsigma-1,el) =  -(0.5_WP+epsilon)*V_rhsAB(1:nsigma-1,el) &
          + Fy*Je(1:nsigma-1,el)

     U_rhsAB(1:nsigma-1,el) = &
          V_n(1:nsigma-1,el)*coriolis(el)*elem_area(el)*Je(1:nsigma-1,el)
     ! +mm*U_n(nz,el)*V_n(nz,el)
     V_rhsAB(1:nsigma-1,el) =&
          -U_n(1:nsigma-1,el)*coriolis(el)*elem_area(el)*Je(1:nsigma-1,el)
     ! -mm*U_n(nz,el)*U_n(nz,el)

  END DO
!$OMP END DO NOWAIT


#ifdef USE_MPI
  call exchange_elem(U_rhs)
  call exchange_elem(V_rhs)
  call exchange_elem(U_rhsAB)
  call exchange_elem(V_rhsAB)
#endif

  ! ====================
  ! Compute velocity gradients
  ! (to be used in viscosity operator
  ! and in flux estimates)
  ! ====================

  call vel_gradients

!$OMP END PARALLEL

  ! ====================
  ! Horizontal advection
  ! ====================

  if(mom_adv_3D == 1) call momentum_adv_upwind
!SH SKIPPED FOR NOW  If(mom_adv_3D == 2) call momentum_adv_p1
!SH SKIPPED FOR NOW  if(mom_adv_3D == 3) call momentum_adv_scalar_3D

  !======
  ! Horizontal viscosity part
  !======

  if (filt_3D) call viscosity_filt_3D
  if (filt_bi_3D) call viscosity_filt2x_3D
  if (bih_3D) call biharmonic_viscosity

  !===============================
  ! Vertical advection
  !===============================
  if (vert_adv == 1) call momentum_vert_adv_upwind
!SH SKIPPED FOR NOW    if (ex_vert_visc) call momentum_vert_expl_visc

  ! =======================
  ! Update the rhs
  ! =======================
  ff=(1.5_WP+epsilon)
  if(lfirst) then
     ff=1.0_WP
     lfirst=.false.
  end if

!$OMP PARALLEL  DEFAULT(SHARED) PRIVATE(el, nz,area_inv)
!$OMP DO

  DO el=1,mydim_elem2D
     area_inv = 1._WP/elem_area(el)
     DO nz=1,nsigma-1
        U_rhs(nz,el) = mask_wd(el)*(U_n(nz,el)  + Bar_pru_3D(nz,el)*dt +&
             (dt*(U_rhs(nz,el)+U_rhsAB(nz,el)*ff)*area_inv +&
             UV_rhs(1,nz,el))/Je(nz,el))
        V_rhs(nz,el) = mask_wd(el)*(V_n(nz,el) + Bar_prv_3D(nz,el)*dt +&
             (dt*(V_rhs(nz,el)+V_rhsAB(nz,el)*ff)*area_inv +&
             UV_rhs(2,nz,el))/Je(nz,el))
     END DO
  END DO

#ifdef USE_MPI
  call exchange_elem(U_rhs)
  call exchange_elem(V_rhs)
#endif

#ifdef USE_MPI
!     call MPI_AllREDUCE(minval(Je(1,:)),mnv, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
!          MPI_COMM_FESOM_C, MPIerr)

!     call MPI_AllREDUCE(maxval(Je(1,:)),mxv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
!          MPI_COMM_FESOM_C, MPIerr)

!IK     if (mype==0) print *,'MIN_VAL rhsAB',mnv
!IK     if (mype==0) print *,'MAX_VAL rhsAB',mxv
#else
!IK     print *,'MIN_VAL rhsAB',minval(Je(1,:))
!IK     print *,'MAX_VAL rhsAB',maxval(Je(1,:))
#endif


#ifdef USE_MPI
!     call MPI_AllREDUCE(minval(U_rhs),mnv, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
!          MPI_COMM_FESOM_C, MPIerr)

!     call MPI_AllREDUCE(maxval(U_rhs),mxv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
!          MPI_COMM_FESOM_C, MPIerr)

!IK     if (mype==0) print *,'MIN_VAL U_rhs',mnv
!IK     if (mype==0) print *,'MAX_VAL U_rhs',mxv
#else
!IK     print *,'MIN_VAL U_rhs',minval(U_rhs)
!IK     print *,'MAX_VAL U_rhs',maxval(U_rhs)
#endif


!$OMP END DO
!$OMP END PARALLEL
!!! comment from here to 0>
  if (im_vert_visc) call momentum_vert_impl_visc
  if (im_vert_visc_fic) then
    call momentum_vert_impl_visc_tmp

  endif
!!! up to here <0
!!! new version of the MOST
!  call vert_impl_MOST

  ! =======================
  ! U_rhs contains all contributions to velocity from old time steps
  ! =======================
  ! call CPU_time(t6)

END SUBROUTINE compute_vel_rhs

!=======================================================================

SUBROUTINE update_3D_vel

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  USE g_PARSUP
  use g_comm_auto

  IMPLICIT NONE

  integer   :: el,  nz, ed
  real(kind=WP) :: a(nsigma-1),  dmean_inv
  real(kind=WP) :: U_2D_p, V_2D_p

!$OMP DO
  do el=1,myDim_elem2D

     a(1:nsigma-1) = w_cv(1,el)*Jc_old(1:nsigma-1,elem2D_nodes(1,el))  &
          + w_cv(2,el)*Jc_old(1:nsigma-1,elem2D_nodes(2,el))  &
          + w_cv(3,el)*Jc_old(1:nsigma-1,elem2D_nodes(3,el))  &
          + w_cv(4,el)*Jc_old(1:nsigma-1,elem2D_nodes(4,el))

     dmean_inv = 1._WP/max(Dmin,  &
          sum(w_cv(1:4,el)*(depth(elem2D_nodes(:,el)) + eta_n(elem2D_nodes(:,el)))))

     U_2D_p = sum(U_rhs(1:nsigma-1,el) *a(1:nsigma-1))*dmean_inv
     V_2D_p = sum(V_rhs(1:nsigma-1,el) *a(1:nsigma-1))*dmean_inv

     U_filt_2D(el) = U_filt_2D(el)*dmean_inv
     V_filt_2D(el) = V_filt_2D(el)*dmean_inv


     a(1:nsigma-1) = a(1:nsigma-1) /Je(1:nsigma-1,el)

     ! velocity for transport equations

     U_n_filt(1:nsigma-1,el)=  a(1:nsigma-1)*U_rhs(1:nsigma-1,el)+ U_filt_2D(el) - U_2D_p
     V_n_filt(1:nsigma-1,el)=  a(1:nsigma-1)*V_rhs(1:nsigma-1,el)+ V_filt_2D(el) - V_2D_p

     U_n(1:nsigma-1,el)=  a(1:nsigma-1)*U_rhs(1:nsigma-1,el) + U_n_2D(1,el) - U_2D_p
     V_n(1:nsigma-1,el)=  a(1:nsigma-1)*V_rhs(1:nsigma-1,el) + U_n_2D(2,el) - V_2D_p

  END DO
!$OMP END DO

  !SH NECESSARY?
#ifdef USE_MPI
  call exchange_elem(U_n_filt)
  call exchange_elem(V_n_filt)
  call exchange_elem(U_n)
  call exchange_elem(V_n)
#endif

end subroutine update_3D_vel
!=================================================================================
SUBROUTINE compute_puls_vel

  !a Alexey Androsov
  !a 15.10.14

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  USE g_PARSUP
  use g_comm_auto

  IMPLICIT NONE

  integer          :: el, elnodes(4)
  real(kind=WP)    :: dmean, dmean_new, a


!SHDB print *,'PULS ssh :',mype,minval(ssh_rhs),maxval(ssh_rhs)
!SHDB print *,'PULS dep :',mype,minval(depth),maxval(depth)
!SHDB print *,'PULS dmn :',mype,minval(dmean_n),maxval(dmean_n)
!SHDB print *,'PULS Dmin:',mype,Dmin

!$OMP DO
  DO el=1,myDim_elem2D

     !NR dmean = max(Dmin,sum(w_cv(1:4,el)*(eta_n(elem2D_nodes(:,el))  +depth(elem2D_nodes(:,el)))))
     !NR was calucated just before in timestep_AB_2D as dmean_n

     !NR dmean_new = max(Dmin,sum(w_cv(1:4,el)*(ssh_rhs(elem2D_nodes(:,el))+depth(elem2D_nodes(:,el)))))
     !NR a = dmean/dmean_new

     a = dmean_n(el) /  max(Dmin,sum(w_cv(1:4,el)*(ssh_rhs(elem2D_nodes(:,el))+depth(elem2D_nodes(:,el)))))
     !a = dmean_n(el) /  max(Dmin,ssh_rhs(elem2D_nodes(:,el)))) !! TEST!!CHECKTHIS

     U_puls(1:nsigma-1,el)= U_n(1:nsigma-1,el) - a*U_n_2D(1,el)
     V_puls(1:nsigma-1,el)= V_n(1:nsigma-1,el) - a*U_n_2D(2,el)

  END DO
!$OMP END DO


#ifdef USE_MPI
  call exchange_elem(U_puls)
  call exchange_elem(V_puls)
#endif

!SHDB print *,'U_puls :',mype,minval(U_puls),maxval(U_puls)
!SHDB print *,'V_puls :',mype,minval(V_puls),maxval(V_puls)

END SUBROUTINE compute_puls_vel
!=========================================================================
SUBROUTINE vert_vel_sigma

  !a Alexey Androsov
  !a 20.10.14

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  use g_parsup
  use g_comm_auto

  IMPLICIT NONE

  integer       :: el(2),el1, n, nz, ed
  real(kind=WP) :: c1(nsigma-1), area_inv,rho_inv, rho
  !$ integer    :: me, nthreads,  chunk
  integer       ::  me_nod1, me_nod2

  integer :: edglim

  me_nod1 = 1
  me_nod2 = myDim_nod2D+eDim_nod2D

#ifdef USE_MPI
  edglim=myDim_edge2D+eDim_edge2D
#else
  edglim=edge2D_in
#endif

  !NR For OpenMP, we compute a domain decomposition for Wvel:
  !$  nthreads = omp_get_num_threads()
  !$  me = omp_get_thread_num()
  !$  chunk = nod2D / nthreads + 1
  !$  me_nod1 = me*chunk +1
  !$  me_nod2 = min((me+1)*chunk, nod2D)

  DO n=me_nod1,me_nod2
     Wvel(1:nsigma-1,n) = 0.0_WP
  ENDDO

  ! ==============
  ! internal edges
  ! ==============

  DO ed=1,edglim

#ifdef USE_MPI
     if (mylist_edge2D(ed)>edge2D_in) cycle
#endif

     !$ if ( (edge_nodes(1,ed) >= me_nod1 .and. edge_nodes(1,ed)<= me_nod2 ) &
     !$  .or.(edge_nodes(2,ed) >= me_nod1 .and. edge_nodes(2,ed)<= me_nod2 ) ) then

     el=edge_tri(:,ed)

     do nz=1,nsigma-1

        c1(nz) = Je(nz,el(1)) * (                                          &
             (V_n_filt(nz,el(1)) - V_filt_2D(el(1)))*edge_cross_dxdy(1,ed) &
             -(U_n_filt(nz,el(1)) - U_filt_2D(el(1)))*edge_cross_dxdy(2,ed))&
             + Je(nz,el(2)) * (                                       &
             -(V_n_filt(nz,el(2)) - V_filt_2D(el(2)))*edge_cross_dxdy(3,ed) &
             +(U_n_filt(nz,el(2)) - U_filt_2D(el(2)))*edge_cross_dxdy(4,ed))
     enddo

     !$ if (edge_nodes(1,ed) >= me_nod1 .and. edge_nodes(1,ed)<= me_nod2 ) then
     Wvel(1:nsigma-1,edge_nodes(1,ed)) = Wvel(1:nsigma-1,edge_nodes(1,ed)) + c1(1:nsigma-1)
     !$ end if
     !$ if (edge_nodes(2,ed) >= me_nod1 .and. edge_nodes(2,ed)<= me_nod2 )  then
     Wvel(1:nsigma-1,edge_nodes(2,ed)) = Wvel(1:nsigma-1,edge_nodes(2,ed)) - c1(1:nsigma-1)
     !$ end if
     !$ endif
  END DO

  ! =============
  ! boundary edges
  ! only the left element (1) is available
  ! =============


#ifdef USE_MPI

  call exchange_nod(Wvel)

  DO ed=1,edglim
     if (mylist_edge2D(ed)>edge2D_in) then
        !$ if ( (edge_nodes(1,ed) >= me_nod1 .and. edge_nodes(1,ed)<= me_nod2 ) &
        !$  .or.(edge_nodes(2,ed) >= me_nod1 .and. edge_nodes(2,ed)<= me_nod2 ) ) then

        el1=edge_tri(1,ed)

        do nz=1,nsigma-1
           c1(nz) = Je(nz,el1)*                           &
                ((V_n_filt(nz,el1) - V_filt_2D(el1))*edge_cross_dxdy(1,ed) &
                -(U_n_filt(nz,el1) - U_filt_2D(el1))*edge_cross_dxdy(2,ed))
        enddo
        !$ if (edge_nodes(1,ed) >= me_nod1 .and. edge_nodes(1,ed)<= me_nod2 ) then
        Wvel(1:nsigma-1,edge_nodes(1,ed)) = Wvel(1:nsigma-1,edge_nodes(1,ed)) + c1(1:nsigma-1)
        !$ end if
        !$ if (edge_nodes(2,ed) >= me_nod1 .and. edge_nodes(2,ed)<= me_nod2 )  then
        Wvel(1:nsigma-1,edge_nodes(2,ed)) = Wvel(1:nsigma-1,edge_nodes(2,ed)) - c1(1:nsigma-1)
        !$ end if
        !$ endif
     end if
  END DO

  call exchange_nod(Wvel)

#else
  DO ed=1+edge2D_in, edge2D
     !$ if ( (edge_nodes(1,ed) >= me_nod1 .and. edge_nodes(1,ed)<= me_nod2 ) &
     !$  .or.(edge_nodes(2,ed) >= me_nod1 .and. edge_nodes(2,ed)<= me_nod2 ) ) then

     el1=edge_tri(1,ed)

     do nz=1,nsigma-1
        c1(nz) = Je(nz,el1)*                           &
             ((V_n_filt(nz,el1) - V_filt_2D(el1))*edge_cross_dxdy(1,ed) &
             -(U_n_filt(nz,el1) - U_filt_2D(el1))*edge_cross_dxdy(2,ed))
     enddo
     !$ if (edge_nodes(1,ed) >= me_nod1 .and. edge_nodes(1,ed)<= me_nod2 ) then
     Wvel(1:nsigma-1,edge_nodes(1,ed)) = Wvel(1:nsigma-1,edge_nodes(1,ed)) + c1(1:nsigma-1)
     !$ end if
     !$ if (edge_nodes(2,ed) >= me_nod1 .and. edge_nodes(2,ed)<= me_nod2 )  then
     Wvel(1:nsigma-1,edge_nodes(2,ed)) = Wvel(1:nsigma-1,edge_nodes(2,ed)) - c1(1:nsigma-1)
     !$ end if
     !$ endif
  END DO
#endif

  ! ===================
  ! Sum up to get W
  ! ===================

  Do n=me_nod1,me_nod2
     area_inv = 1._WP/area(n)
     DO nz=nsigma-1,1,-1
        Wvel(nz,n) = Wvel(nz,n) + Wvel(nz+1,n)
     END DO
     Wvel(1:nsigma-1,n) = Wvel(1:nsigma-1,n) * area_inv
  END DO

#ifdef USE_MPI
  call exchange_nod(Wvel)
#endif

 end subroutine vert_vel_sigma
!==========================================================================
SUBROUTINE vert_vel_cart

USE o_MESH
USE o_ARRAYS
USE o_PARAM
IMPLICIT NONE
!
! -div(Hu) is computed in cycle over edges
!
! The depth is estimated at elements
! c1 is -u_n*L_left*d1, c2  -u_n*L_right*d1

integer                      :: ed, el(2),el1,enodes(2), elnodes(4), elem, n, nz
real(kind=WP)           :: c1, c2, fD(4), delta, a, b, wad_nodes
real(Kind=WP), allocatable   :: auh_elem(:,:), w_rhs(:,:)

allocate(auh_elem(nsigma-1,elem2D), w_rhs(nsigma-1,nod2D))
w_rhs=0.0_WP
W_n=0.0_WP
 ! ==============
 ! fill mean depth
 ! ==============

  DO elem=1,elem2D
    elnodes=elem2D_nodes(:,elem)
    do nz=1,nsigma-1
      delta = 0.5_WP*(sigma(nz) + sigma(nz+1))
     fD=(depth(elnodes) + eta_n(elnodes))*delta - depth(elnodes)
     auh_elem(nz,elem)=sum(w_cv(1:4,elem)*fD)*Je(nz,elem)
    enddo
 END DO
 ! ==============
 ! internal edges
 ! ==============
 DO ed=1,edge2D_in
    enodes=edge_nodes(:,ed)
    el=edge_tri(:,ed)
     do nz=1,nsigma-1
    c1=V_n(nz,el(1))*edge_cross_dxdy(1,ed)-U_n(nz,el(1))*edge_cross_dxdy(2,ed)
    c2=-V_n(nz,el(2))*edge_cross_dxdy(3,ed)+U_n(nz,el(2))*edge_cross_dxdy(4,ed)
    c1=c1*auh_elem(nz,el(1)) + c2*auh_elem(nz,el(2))
    w_rhs(nz,enodes(1))=w_rhs(nz,enodes(1)) + c1
    w_rhs(nz,enodes(2))=w_rhs(nz,enodes(2)) - c1
     enddo
 END DO
 ! =============
 ! boundary edges
 ! only the left element (1) is available
 ! =============
  DO ed=1+edge2D_in, edge2D
    enodes=edge_nodes(:,ed)
    el1=edge_tri(1,ed)
    do nz=1,nsigma-1
    c1=V_n(nz,el1)*edge_cross_dxdy(1,ed)-U_n(nz,el1)*edge_cross_dxdy(2,ed)
    c1=c1*auh_elem(nz,el1)
    w_rhs(nz,enodes(1))=w_rhs(nz,enodes(1)) + c1
    w_rhs(nz,enodes(2))=w_rhs(nz,enodes(2)) - c1
    enddo
 END DO
!++++++++++++++++++++++++++++++++++
! cartesian vertical velocity
!++++++++++++++++++++++++++++++++++
  do n=1,nod2D
   a = depth(n) + eta_n(n)
   b = depth(n) + eta_p(n)
   if (a <= Dmin) wad_nodes = 0.0_WP
   if (a > Dmin)   wad_nodes = 1.0_WP
   do nz=1,nsigma-1
    delta = 0.5_WP*(sigma(nz) + sigma(nz+1))
    c1 = a*delta - depth(n)
    c2 = b*delta - depth(n)
   if (index_nod2D(n) /= 2) W_n(nz,n)=wad_nodes*((c1*Jc(nz,n)-c2*Jc_old(nz,n))/dt +&
                                                       w_rhs(nz,n)/area(n) + &
                                                       (Wvel(nz,n) - Wvel(nz+1,n))*c1)/Jc(nz,n)
   enddo
  ENDDO
deallocate(auh_elem, w_rhs)

 end subroutine vert_vel_cart
!==========================================================================
subroutine compute_vortex_3D
  use o_MESH
  use o_ARRAYS
  use o_PARAM

  implicit none
  integer               :: ed, el(2), enodes(2) , n, nz
  real(kind=WP)     :: deltaX1, deltaX2, deltaY1, deltaY2, c1

 vorticity_3D=0.0_WP
do ed=1,edge2D
   enodes=edge_nodes(:,ed)
   el=edge_tri(:,ed)
   deltaX1=edge_cross_dxdy(1,ed)
   deltaY1=edge_cross_dxdy(2,ed)
   if(el(2)>0) then
    deltaX2=edge_cross_dxdy(3,ed)
    deltaY2=edge_cross_dxdy(4,ed)
   endif
    do nz=1,nsigma-1
   c1=deltaX1*U_n(nz,el(1))+deltaY1*V_n(nz,el(1))
   if(el(2)>0) then
    c1=c1-deltaX2*U_n(nz,el(2))-deltaY2*V_n(nz,el(2))
   endif
   vorticity_3D(nz,enodes(1))=vorticity_3D(nz,enodes(1))+c1
   vorticity_3D(nz,enodes(2))=vorticity_3D(nz,enodes(2))-c1
   end do
end do
do n=1,nod2D
  do nz=1,nsigma-1
   vorticity_3D(nz,n)=vorticity_3D(nz,n)/area(n)
  enddo
enddo
end subroutine compute_vortex_3D
!=======================================================================
