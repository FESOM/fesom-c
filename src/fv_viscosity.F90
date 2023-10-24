! ===========================================================================
SUBROUTINE h_viscosity
  !
  ! Coefficient of horizontal viscosity is a combination of
  ! Smagorinsky contribution (with Smag2)
  ! and background (with A_hor)
  ! Background is scaled as sqrt(A/A_0), others are scaled
  ! in natural way by construction

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  use g_parsup
  use g_comm_auto

  IMPLICIT NONE

  real(kind=WP)  :: Smag2, A_hor_scale, d1, d2, acc
  integer             :: elem, nz

  Smag2=C_Smag*C_Smag

  ! Fill in viscosity:
  DO  elem=1,myDim_elem2D
     acc = sum(w_cv(1:4,elem)*ac(elem2D_nodes(:,elem)))
     A_hor_scale = A_hor*sqrt(elem_area(elem)/scale_area)
     do nz=1,nsigma-1
        d1=vel_grad_ux(nz,elem)-vel_grad_vy(nz,elem)
        d2=vel_grad_uy(nz,elem)+vel_grad_vx(nz,elem)
        Visc(nz,elem)=acc*( elem_area(elem)* &
	              sqrt(Smag2*(d1**2+d2**2))+ A_hor_scale)
     end do
  END DO

#ifdef USE_MPI
  call exchange_elem(Visc)
#endif

END subroutine h_viscosity
! ===========================================================================
SUBROUTINE h_viscosity_2D
  !
  ! Coefficient of horizontal viscosity is a combination of
  ! Smagorinsky contribution (with Smag2)
  ! and background (with A_hor)
  ! Background is scaled as sqrt(A/A_0), others are scaled
  ! in natural way by construction

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  USE g_PARSUP
  use g_comm_auto

  IMPLICIT NONE

  real(kind=WP)  :: d1, d2, acc, A_vel
  integer             :: elem

  ! Fill in viscosity:
  DO  elem=1, myDim_elem2D
     acc = sum(w_cv(1:4,elem)*ac(elem2D_nodes(:,elem)))
     A_vel = A_hor !2D*elem_area(elem)/(2.*dt_2D)  !75.*ddd
     d1=vel_grad(1,elem)-vel_grad(4,elem)
     d2=vel_grad(2,elem)+vel_grad(3,elem)
     Visc2D(elem)=acc*(A_vel+ elem_area(elem)*C_Smag*sqrt(d1**2+d2**2))
  END DO
  ! Here multiplication with depth is just for simplicity,
  ! in reality the contribution from viscosity to the rhs
  ! has to be multiplied.

#ifdef USE_MPI
  call exchange_elem(Visc2D)
#endif


END subroutine h_viscosity_2D

!===========================================================================================
!Biharmonic viscosity routines.
! Contains: biharmonic_viscosity (only biharmonic operator)
!           biPlusharmonic_viscosity (harmonic+biharmonic)
!           vel_lapl_gradiens ==vel_gradients, but simplified (no account for free slip).
!
! True viscosity operator is implemented in laplacian_viscosity (div sigma_ij,
! with sigma_ij=2*A_hor(epsilon_ij-delta_ij*div u/2). If we apply it twice, the resultant
! biharmonic operator does not kill noise if eddies are strong. We use operator
! Lu=nabla nabla u, corrected as suggested by Blazek, 2001, (the correction concerns the
! way nabla u is computed on cell-centered discretization) and apply it twice:
! bh_visc=-L Abh Lu.
!
! If Laplacian and biharmonic are used simultaneously, we assemble
! lapl_visc=nabla A_hor nabla u and bh_visc in the same subroutine to
! increase numerical efficiency. Note, however, that operators introduced in that way are
! not zero under solid body rotations, while the operator in laplacian_viscosity is.
!
! Note also that we ignore metric differentiation terms because our viscosity
! is just a parameterization. There is no problem with including them.
!
! sergey.danilov@awi.de  2012
! With Laplacian regularization
SUBROUTINE biharmonic_viscosity

  !a changes have been made by Alexey Androsov
  !a (for sigma coordinates, and for tri-quad elements)
  !a 09.10.14

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  use g_parsup
  use g_comm_auto

  IMPLICIT NONE

  real(kind=WP)          :: xe, ye, u1, v1, Ah, acc
  real(kind=WP)          :: tx, ty, tt, g1, g2, tl
  real(kind=WP), allocatable  :: U_c(:,:), V_c(:,:)
  integer               :: ed, el(2), nz, elem

  integer :: edglim
  integer            :: elem_size!,node_size

  allocate(U_c(nsigma-1,myDim_elem2D+eDim_elem2D+eXDim_elem2D), &
           V_c(nsigma-1,myDim_elem2D+eDim_elem2D+eXDim_elem2D))

  !Acall h_viscosity

  U_c=0.0_WP
  V_c=0.0_WP

#ifdef USE_MPI
  edglim=myDim_edge2D+eDim_edge2D
#else
  edglim=edge2D_in
#endif

  DO ed=1, edglim

#ifdef USE_MPI
     if (myList_edge2D(ed)>edge2D_in) cycle
#endif

     el=edge_tri(:,ed)
     xe=edge_dxdy(1,ed)*r_earth*0.5_WP*(elem_cos(el(1))+elem_cos(el(2)))
     ye=edge_dxdy(2,ed)*r_earth

     tx=-edge_cross_dxdy(1,ed)+edge_cross_dxdy(3,ed)
     ty=-edge_cross_dxdy(2,ed)+edge_cross_dxdy(4,ed)
     tl=sqrt(tx**2+ty**2)
     tx=tx/tl
     ty=ty/tl

     DO nz=1, nsigma-1
        g1=0.5_WP*(vel_grad_ux(nz,el(1))+vel_grad_ux(nz,el(2)))
        g2=0.5_WP*(vel_grad_uy(nz,el(1))+vel_grad_uy(nz,el(2)))
        tt=g1*tx+g2*ty-(U_n(nz,el(2))-U_n(nz,el(1)))/tl
        g1=g1-tx*tt
        g2=g2-ty*tt
        u1=(g1*ye-g2*xe)
        g1=0.5_WP*(vel_grad_vx(nz,el(1))+vel_grad_vx(nz,el(2)))
        g2=0.5_WP*(vel_grad_vy(nz,el(1))+vel_grad_vy(nz,el(2)))
        tt=g1*tx+g2*ty-(V_n(nz,el(2))-V_n(nz,el(1)))/tl
        g1=g1-tx*tt
        g2=g2-ty*tt
        v1=(g1*ye-g2*xe)

        U_c(nz,el(1))=U_c(nz,el(1))+u1*Jd(nz,ed)
        V_c(nz,el(1))=V_c(nz,el(1))+v1*Jd(nz,ed)
        U_c(nz,el(2))=U_c(nz,el(2))-u1*Jd(nz,ed)
        V_c(nz,el(2))=V_c(nz,el(2))-v1*Jd(nz,ed)

     END DO
  END DO

  ! ============
  ! Contribution from boundary edges
  ! ============

#ifdef USE_MPI

  call exchange_elem(U_c) !SH check this
  call exchange_elem(V_c) !SH check this

  DO ed=1,edglim
     if (myList_edge2D(ed)>edge2D_in) then
        el=edge_tri(:,ed)
        xe=edge_dxdy(1,ed)*r_earth*elem_cos(el(1))
        ye=edge_dxdy(2,ed)*r_earth

        DO nz=1,nsigma-1
           u1=(vel_grad_ux(nz,el(1)))*ye-(vel_grad_uy(nz,el(1)))*xe
           v1=(vel_grad_vx(nz,el(1)))*ye-(vel_grad_vy(nz,el(1)))*xe
           U_c(nz,el(1))=U_c(nz,el(1))+u1*Jd(nz,ed)
           V_c(nz,el(1))=V_c(nz,el(1))+v1*Jd(nz,ed)
        END DO
     end if
  END DO

  call exchange_elem(U_c)
  call exchange_elem(V_c)

#else

  DO ed=1+edge2d_in,edge2D
     el=edge_tri(:,ed)
     xe=edge_dxdy(1,ed)*r_earth*elem_cos(el(1))
     ye=edge_dxdy(2,ed)*r_earth

     DO nz=1,nsigma-1
        u1=(vel_grad_ux(nz,el(1)))*ye-(vel_grad_uy(nz,el(1)))*xe
        v1=(vel_grad_vx(nz,el(1)))*ye-(vel_grad_vy(nz,el(1)))*xe
        U_c(nz,el(1))=U_c(nz,el(1))+u1*Jd(nz,ed)
        V_c(nz,el(1))=V_c(nz,el(1))+v1*Jd(nz,ed)
     END DO
  END DO

#endif

  ! =============
  ! Multiply with biharmonic viscosity
  ! =============

  DO elem=1, myDim_elem2D
     ! we have to divide on elem_area
     ! - takes into account the fact that
     ! biharmonic goes with different sign.
     !     acc = sum(w_cv(1:4,elem)*ac(elem2D_nodes(:,elem)))
     Ah=-Abh0*sqrt(elem_area(elem)/scale_area)/scale_area
     U_c(:,elem)=U_c(:,elem)*(Ah) !-Visc(nz,elem))
     V_c(:,elem)=V_c(:,elem)*(Ah) !-Visc(nz,elem))
  END DO

#ifdef USE_MPI
  call exchange_elem(U_c)
  call exchange_elem(V_c)
#endif

  call vel_lapl_gradients(U_c, V_c)

  ! =============
  ! Apply Laplace operator once more
  ! =============

  DO ed=1, edglim
#ifdef USE_MPI
     if (myList_edge2D(ed)>edge2D_in) cycle
#endif

    el=edge_tri(:,ed)
   xe=edge_dxdy(1,ed)*r_earth*0.5_WP*(elem_cos(el(1))+elem_cos(el(2)))
   ye=edge_dxdy(2,ed)*r_earth

   tx=-edge_cross_dxdy(1,ed)+edge_cross_dxdy(3,ed)
   ty=-edge_cross_dxdy(2,ed)+edge_cross_dxdy(4,ed)
   tl=sqrt(tx**2+ty**2)
   tx=tx/tl
   ty=ty/tl
   acc=sum(ac(edge_nodes(:,ed)))/2.0_WP


   DO nz=1, nsigma-1
   g1=0.5_WP*(vel_grad_ux(nz,el(1))+vel_grad_ux(nz,el(2)))
   g2=0.5_WP*(vel_grad_uy(nz,el(1))+vel_grad_uy(nz,el(2)))
   tt=g1*tx+g2*ty-(U_c(nz,el(2))-U_c(nz,el(1)))/tl
   g1=g1-tx*tt
   g2=g2-ty*tt
   u1=(g1*ye-g2*xe)
   g1=0.5_WP*(vel_grad_vx(nz,el(1))+vel_grad_vx(nz,el(2)))
   g2=0.5_WP*(vel_grad_vy(nz,el(1))+vel_grad_vy(nz,el(2)))
   tt=g1*tx+g2*ty-(V_c(nz,el(2))-V_c(nz,el(1)))/tl
   g1=g1-tx*tt
   g2=g2-ty*tt
   v1=(g1*ye-g2*xe)

   U_rhs(nz,el(1))=U_rhs(nz,el(1))+u1*Jd(nz,ed)*acc
   V_rhs(nz,el(1))=V_rhs(nz,el(1))+v1*Jd(nz,ed)*acc
   U_rhs(nz,el(2))=U_rhs(nz,el(2))-u1*Jd(nz,ed)*acc
   V_rhs(nz,el(2))=V_rhs(nz,el(2))-v1*Jd(nz,ed)*acc

   END DO
END DO

! ============
! Contribution from boundary edges
! ============
#ifdef USE_MPI

  call exchange_elem(U_rhs) !SH check this
  call exchange_elem(V_rhs) !SH check this

  DO ed=1,edglim
        if (myList_edge2D(ed)>edge2D_in) then
            el=edge_tri(:,ed)
            xe=edge_dxdy(1,ed)*r_earth*elem_cos(el(1))
            ye=edge_dxdy(2,ed)*r_earth
            acc = sum(ac(edge_nodes(:,ed)))/2.0_WP


   DO nz=1, nsigma-1
                u1=(vel_grad_ux(nz,el(1)))*ye-(vel_grad_uy(nz,el(1)))*xe
                v1=(vel_grad_vx(nz,el(1)))*ye-(vel_grad_vy(nz,el(1)))*xe
                U_rhs(nz,el(1))=U_rhs(nz,el(1))+u1*Jd(nz,ed)*acc
                V_rhs(nz,el(1))=V_rhs(nz,el(1))+v1*Jd(nz,ed)*acc
   END DO

        end if
  END DO

  call exchange_elem(U_rhs)
  call exchange_elem(V_rhs)

#else
 DO ed=1+edge2D_in, edge2D
   el=edge_tri(:,ed)
   xe=edge_dxdy(1,ed)*r_earth*elem_cos(el(1))
   ye=edge_dxdy(2,ed)*r_earth
   acc = sum(ac(edge_nodes(:,ed)))/2.0_WP


   DO nz=1, nsigma-1
   u1=(vel_grad_ux(nz,el(1)))*ye-(vel_grad_uy(nz,el(1)))*xe
   v1=(vel_grad_vx(nz,el(1)))*ye-(vel_grad_vy(nz,el(1)))*xe
   U_rhs(nz,el(1))=U_rhs(nz,el(1))+u1*Jd(nz,ed)*acc
   V_rhs(nz,el(1))=V_rhs(nz,el(1))+v1*Jd(nz,ed)*acc
   END DO
END DO

#endif
deallocate(V_c, U_c)
END subroutine biharmonic_viscosity
! =========================================================================
! compute biharmonic diffusion for 2D velocity
!
SUBROUTINE biharmonic_2D
USE o_MESH
USE o_ARRAYS
USE o_PARAM

  use g_parsup
  use g_comm_auto

IMPLICIT NONE

real(kind=WP)  :: xe, ye, u1, v1, xed
real(kind=WP)  :: Ah, tl, tt, tx, ty, g1, g2, acc !, fD(4)
real(kind=WP), allocatable :: Ul(:,:)
integer              :: ed, el(2),elem, elnodes(4)

  integer :: edglim

 ! This routine applies twice viscosity2

  allocate(Ul(2,myDim_elem2D+eDim_elem2D+eXDim_elem2D))
 Ul=0.0_WP
 ! =================
 ! Compute Laplacian of velocity
 ! =================

#ifdef USE_MPI
  edglim=myDim_edge2D+eDim_edge2D
#else
  edglim=edge2D_in
#endif

  DO  ed=1,edglim

#ifdef USE_MPI
     if (myList_edge2D(ed)>edge2D_in) cycle
#endif

   el=edge_tri(:,ed)
   xe=edge_dxdy(1,ed)*r_earth*0.5_WP*(elem_cos(el(1))+elem_cos(el(2)))
   ye=edge_dxdy(2,ed)*r_earth

   tx=-edge_cross_dxdy(1,ed)+edge_cross_dxdy(3,ed)
   ty=-edge_cross_dxdy(2,ed)+edge_cross_dxdy(4,ed)
   tl=sqrt(tx**2+ty**2)
   tx=tx/tl
   ty=ty/tl   ! this is the vector of normal to the edge

   g1=0.5_WP*(vel_grad(1,el(1))+vel_grad(1,el(2)))
   g2=0.5_WP*(vel_grad(2,el(1))+vel_grad(2,el(2)))
   tt=g1*tx+g2*ty-(UAB(1,el(2))-UAB(1,el(1)))/tl
   g1=g1-tx*tt
   g2=g2-ty*tt
   u1=(g1*ye-g2*xe)
   g1=0.5_WP*(vel_grad(3,el(1))+vel_grad(3,el(2)))
   g2=0.5_WP*(vel_grad(4,el(1))+vel_grad(4,el(2)))
   tt=g1*tx+g2*ty-(UAB(2,el(2))-UAB(2,el(1)))/tl
   g1=g1-tx*tt
   g2=g2-ty*tt
   v1=(g1*ye-g2*xe)

   Ul(1,el(1))=Ul(1,el(1))+u1
   Ul(2,el(1))=Ul(2,el(1))+v1
   Ul(1,el(2))=Ul(1,el(2))-u1
   Ul(2,el(2))=Ul(2,el(2))-v1
 END DO


#ifdef USE_MPI

  call exchange_elem(Ul) !SH check this

  DO ed=1,edglim
     if (myList_edge2D(ed)>edge2D_in) then

        el=edge_tri(:,ed)
        xe=edge_dxdy(1,ed)*r_earth
        xed=xe*elem_cos(el(1))
        ye=edge_dxdy(2,ed)*r_earth
        u1=(vel_grad(1,el(1)))*ye-(vel_grad(2,el(1)))*xed
        v1=(vel_grad(3,el(1)))*ye-(vel_grad(4,el(1)))*xed

        if (free_slip) then
           !   remove tangent component
           !   Projection; division because xe, ye are not normalized
           u1=ye*(u1*ye-v1*xe)/(xe*xe+ye*ye)
           v1=-xe*(u1*ye-v1*xe)/(xe*xe+ye*ye)
        end if

        Ul(1,el(1))=Ul(1,el(1))+u1
        Ul(2,el(1))=Ul(2,el(1))+v1

     end if
  END DO

  call exchange_elem(Ul)

#else

 DO ed=1+edge2D_in,edge2D
   el=edge_tri(:,ed)
   xe=edge_dxdy(1,ed)*r_earth
   xed=xe*elem_cos(el(1))
   ye=edge_dxdy(2,ed)*r_earth
   u1=(vel_grad(1,el(1)))*ye-(vel_grad(2,el(1)))*xed
   v1=(vel_grad(3,el(1)))*ye-(vel_grad(4,el(1)))*xed

   if (free_slip) then
 !   remove tangent component
 !   Projection; division because xe, ye are not normalized
   u1=ye*(u1*ye-v1*xe)/(xe*xe+ye*ye)
   v1=-xe*(u1*ye-v1*xe)/(xe*xe+ye*ye)
   end if
   Ul(1,el(1))=Ul(1,el(1))+u1
   Ul(2,el(1))=Ul(2,el(1))+v1
 END DO

#endif

  !call h_viscosity
 ! =============
 ! Multiply with biharmonic viscosity
 ! =============
 call h_viscosity_2D

  DO elem=1,myDim_elem2D
     elnodes=elem2D_nodes(:,elem)
     acc = sum(w_cv(1:4,elem)*ac(elnodes))
     Ah=-Abh0*sqrt(elem_area(elem)/scale_area)/scale_area
 !a    Ah=Ah*max(1.0_WP, sum(depth(elnodes))/4.0_WP)
    !  fD=depth(elnodes)
    ! Ah=Ah*max(1.0_8, sum(w_cv(1:4,elem)*fD))
     Ul(:,elem)=Ul(:,elem)*(Ah)*acc !-Visc2D(elem))
 END DO

#ifdef USE_MPI
  call exchange_elem(Ul)
#endif

  call vel_gradients_2D(Ul)
  ! ===============
 ! Apply Laplacian second time
 ! ===============

  DO  ed=1,edglim

#ifdef USE_MPI
     if (myList_edge2D(ed)>edge2D_in) cycle
#endif

   el=edge_tri(:,ed)

   xe=edge_dxdy(1,ed)*r_earth*0.5_WP*(elem_cos(el(1))+elem_cos(el(2)))
   ye=edge_dxdy(2,ed)*r_earth
   tx=-edge_cross_dxdy(1,ed)+edge_cross_dxdy(3,ed)
   ty=-edge_cross_dxdy(2,ed)+edge_cross_dxdy(4,ed)
   tl=sqrt(tx**2+ty**2)
   tx=tx/tl
   ty=ty/tl    ! this is the vector of normal to the edge

   g1=0.5_WP*(vel_grad(1,el(1))+vel_grad(1,el(2)))
   g2=0.5_WP*(vel_grad(2,el(1))+vel_grad(2,el(2)))
   tt=g1*tx+g2*ty-(Ul(1,el(2))-Ul(1,el(1)))/tl
   g1=g1-tx*tt
   g2=g2-ty*tt
   u1=(g1*ye-g2*xe)
   g1=0.5_WP*(vel_grad(3,el(1))+vel_grad(3,el(2)))
   g2=0.5_WP*(vel_grad(4,el(1))+vel_grad(4,el(2)))
   tt=g1*tx+g2*ty-(Ul(2,el(2))-Ul(2,el(1)))/tl
   g1=g1-tx*tt
   g2=g2-ty*tt
   v1=(g1*ye-g2*xe)

   U_rhs_2D(1,el(1))=U_rhs_2D(1,el(1))+u1
   U_rhs_2D(2,el(1))=U_rhs_2D(2,el(1))+v1
   U_rhs_2D(1,el(2))=U_rhs_2D(1,el(2))-u1
   U_rhs_2D(2,el(2))=U_rhs_2D(2,el(2))-v1
 END DO


#ifdef USE_MPI

  call exchange_elem(U_rhs_2D) !SH check this

  DO ed=1,edglim
     if (myList_edge2D(ed)>edge2D_in) then

        el=edge_tri(:,ed)
        xe=edge_dxdy(1,ed)*r_earth
        xed=xe*elem_cos(el(1))
        ye=edge_dxdy(2,ed)*r_earth
        u1=(vel_grad(1,el(1)))*ye-(vel_grad(2,el(1)))*xed
        v1=(vel_grad(3,el(1)))*ye-(vel_grad(4,el(1)))*xed
        if (free_slip) then
           ! remove tangent component
           ! Projection; division because xe, ye are not normalized
           u1=ye*(u1*ye-v1*xe)/(xe*xe+ye*ye)
           v1=-xe*(u1*ye-v1*xe)/(xe*xe+ye*ye)
        end if

        U_rhs_2D(1,el(1))=U_rhs_2D(1,el(1))+u1
        U_rhs_2D(2,el(1))=U_rhs_2D(2,el(1))+v1

     end if
  END DO

  call exchange_elem(U_rhs_2D)

#else

 DO ed=1+edge2D_in,edge2D
   el=edge_tri(:,ed)
   xe=edge_dxdy(1,ed)*r_earth
   xed=xe*elem_cos(el(1))
   ye=edge_dxdy(2,ed)*r_earth
   u1=(vel_grad(1,el(1)))*ye-(vel_grad(2,el(1)))*xed
   v1=(vel_grad(3,el(1)))*ye-(vel_grad(4,el(1)))*xed
   if (free_slip) then
   ! remove tangent component
   ! Projection; division because xe, ye are not normalized
   u1=ye*(u1*ye-v1*xe)/(xe*xe+ye*ye)
   v1=-xe*(u1*ye-v1*xe)/(xe*xe+ye*ye)
   end if

   U_rhs_2D(1,el(1))=U_rhs_2D(1,el(1))+u1
   U_rhs_2D(2,el(1))=U_rhs_2D(2,el(1))+v1
 END DO

#endif

 deallocate(Ul)
end subroutine biharmonic_2D
! ===================================================================
SUBROUTINE viscosity_filt2x_3D
USE o_MESH
USE o_ARRAYS
USE o_PARAM

  use g_parsup
  use g_comm_auto

IMPLICIT NONE

integer         :: el, eln(4), nz
real(kind=WP)   :: tau_inv, K_bi
  real(kind=WP)   :: U_c(nsigma-1,myDim_elem2D+eDim_elem2D+eXDim_elem2D)
  real(kind=WP)   :: V_c(nsigma-1,myDim_elem2D+eDim_elem2D+eXDim_elem2D)

!tau_c = 0.4_WP  ! 0.4_WP
tau_inv=dt*tau_c/3600.0/24     ! SET IT experimentally

 ! Filter is applied twice. It should be approximately
 ! equivalent to biharmonic operator with the coefficient
 ! (tau_c/day)a^3/9. Scaling inside is found to help
 ! with smoothness in places of mesh transition. *(it makes a^3 from a^4)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(el,eln,nz,K_bi)
!$OMP DO
  DO el=1,myDim_elem2D
   K_bi = tau_inv*sqrt(scale_area/elem_area(el))

   eln(1:4) = elem_neighbors(1:4,el)

   if ( eln(1) > 0) then
      DO nz=1,nsigma-1
         U_c(nz,el) =           U_n(nz,el)*Je(nz,el) - U_n(nz,eln(1))*Je(nz,eln(1))
         V_c(nz,el) =           V_n(nz,el)*Je(nz,el) - V_n(nz,eln(1))*Je(nz,eln(1))
      END DO
   else
      U_c(1:nsigma-1,el) = 0._WP
      V_c(1:nsigma-1,el) = 0._WP
   endif
   if ( eln(2) > 0) then
      DO nz=1,nsigma-1
         U_c(nz,el) = U_c(nz,el) + U_n(nz,el)*Je(nz,el) - U_n(nz,eln(2))*Je(nz,eln(2))
         V_c(nz,el) = V_c(nz,el) + V_n(nz,el)*Je(nz,el) - V_n(nz,eln(2))*Je(nz,eln(2))
      END DO
   endif
   if ( eln(3) > 0) then
      DO nz=1,nsigma-1
         U_c(nz,el) = U_c(nz,el) + U_n(nz,el)*Je(nz,el) - U_n(nz,eln(3))*Je(nz,eln(3))
         V_c(nz,el) = V_c(nz,el) + V_n(nz,el)*Je(nz,el) - V_n(nz,eln(3))*Je(nz,eln(3))
      END DO
   endif
   if ( eln(4) > 0  .and. eln(4)/=eln(1)) then
      DO nz=1,nsigma-1
         U_c(nz,el) = U_c(nz,el) + U_n(nz,el)*Je(nz,el) - U_n(nz,eln(4))*Je(nz,eln(4))
         V_c(nz,el) = V_c(nz,el) + V_n(nz,el)*Je(nz,el) - V_n(nz,eln(4))*Je(nz,eln(4))
      END DO
   endif

   U_c(1:nsigma-1,el) = U_c(1:nsigma-1,el)*K_bi
   V_c(1:nsigma-1,el) = V_c(1:nsigma-1,el)*K_bi
ENDDO
!$OMP END DO

#ifdef USE_MPI
  call exchange_elem(U_c)
  call exchange_elem(V_c)
#endif


!$OMP DO
  DO el=1,myDim_elem2D
   eln(1:4) = elem_neighbors(1:4,el)

   if ( eln(1) > 0) then
      UV_rhs(1,1:nsigma-1,el) =               - U_c(1:nsigma-1,el) + U_c(1:nsigma-1,eln(1))
      UV_rhs(2,1:nsigma-1,el) =               - V_c(1:nsigma-1,el) + V_c(1:nsigma-1,eln(1))
   else
      UV_rhs(1:2,1:nsigma-1,el) = 0._WP
   endif
   if ( eln(2) > 0) then
      UV_rhs(1,1:nsigma-1,el) = UV_rhs(1,1:nsigma-1,el) - U_c(1:nsigma-1,el) + U_c(1:nsigma-1,eln(2))
      UV_rhs(2,1:nsigma-1,el) = UV_rhs(2,1:nsigma-1,el) - V_c(1:nsigma-1,el) + V_c(1:nsigma-1,eln(2))
   endif
   if ( eln(3) > 0) then
      UV_rhs(1,1:nsigma-1,el) = UV_rhs(1,1:nsigma-1,el) - U_c(1:nsigma-1,el) + U_c(1:nsigma-1,eln(3))
      UV_rhs(2,1:nsigma-1,el) = UV_rhs(2,1:nsigma-1,el) - V_c(1:nsigma-1,el) + V_c(1:nsigma-1,eln(3))
   endif
   if ( eln(4) > 0  .and. eln(4)/=eln(1)) then
      UV_rhs(1,1:nsigma-1,el) = UV_rhs(1,1:nsigma-1,el) - U_c(1:nsigma-1,el) + U_c(1:nsigma-1,eln(4))
      UV_rhs(2,1:nsigma-1,el) = UV_rhs(2,1:nsigma-1,el) - V_c(1:nsigma-1,el) + V_c(1:nsigma-1,eln(4))
   endif

ENDDO
!$OMP END DO
!$OMP END PARALLEL

#ifdef USE_MPI
  call exchange_elem(UV_rhs) !SH Size of the array correct order?
#endif

end subroutine viscosity_filt2x_3D
!=======================================================================
SUBROUTINE viscosity_filt2x_2D
USE o_MESH
USE o_ARRAYS
USE o_PARAM

  use g_parsup
  use g_comm_auto

IMPLICIT NONE

real(kind=WP) :: tau_inv, K_bi, u1, v1
  real(kind=WP)   :: U_c(myDim_elem2D+eDim_elem2D+eXDim_elem2D)
  real(kind=WP)   :: V_c(myDim_elem2D+eDim_elem2D+eXDim_elem2D)


integer       :: els(2), ed, el

 ! Filter is applied twice. It should be approximately
 ! equivalent to biharmonic operator with the coefficient
 ! (tau_c/day)a^3/9. Scaling inside is found to help
 ! with smoothness in places of mesh transition. *(it makes a^3 from a^4)

! tau_c = 0.4_WP ! 0.4_WP
 tau_inv=dt_2D*tau_c/3600.0/24.0     ! SET IT experimentally

!NROMP  A more compact data structure would be nice with
!NROMP  direct access to the "non-virtual" neighbors.
!NROMP  However, here is not the crucial part for computation time.

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(el,K_bi)
!$OMP DO
  DO el=1,myDim_elem2D

   U_c(el) = 0._WP
   V_c(el) = 0._WP
   K_bi = tau_inv*sqrt(scale_area/elem_area(el))

   if ( elem_neighbors(1,el) > 0) then
      U_c(el) =           U_n_2D(1,el) - U_n_2D(1,elem_neighbors(1,el))
      V_c(el) =           U_n_2D(2,el) - U_n_2D(2,elem_neighbors(1,el))
   endif
   if ( elem_neighbors(2,el) > 0) then
      U_c(el) = U_c(el) + U_n_2D(1,el) - U_n_2D(1,elem_neighbors(2,el))
      V_c(el) = V_c(el) + U_n_2D(2,el) - U_n_2D(2,elem_neighbors(2,el))
   endif
   if ( elem_neighbors(3,el) > 0) then
      U_c(el) = U_c(el) + U_n_2D(1,el) - U_n_2D(1,elem_neighbors(3,el))
      V_c(el) = V_c(el) + U_n_2D(2,el) - U_n_2D(2,elem_neighbors(3,el))
   endif
   if ( elem_neighbors(4,el) > 0  .and. elem_neighbors(4,el)/=elem_neighbors(1,el)) then
      U_c(el) = U_c(el) + U_n_2D(1,el) - U_n_2D(1,elem_neighbors(4,el))
      V_c(el) = V_c(el) + U_n_2D(2,el) - U_n_2D(2,elem_neighbors(4,el))
   endif

   U_c(el) = U_c(el)*K_bi
   V_c(el) = V_c(el)*K_bi
ENDDO
!$OMP END DO

#ifdef USE_MPI
  call exchange_elem(U_c)
  call exchange_elem(V_c)
#endif

!$OMP DO
  DO el=1,myDim_elem2D
   UV2_rhs(1,el) = 0._WP
   UV2_rhs(2,el) = 0._WP

   if ( elem_neighbors(1,el) > 0) then
      UV2_rhs(1,el) =               - U_c(el) + U_c(elem_neighbors(1,el))
      UV2_rhs(2,el) =               - V_c(el) + V_c(elem_neighbors(1,el))
   endif
   if ( elem_neighbors(2,el) > 0) then
      UV2_rhs(1,el) = UV2_rhs(1,el) - U_c(el) + U_c(elem_neighbors(2,el))
      UV2_rhs(2,el) = UV2_rhs(2,el) - V_c(el) + V_c(elem_neighbors(2,el))
   endif
   if ( elem_neighbors(3,el) > 0) then
      UV2_rhs(1,el) = UV2_rhs(1,el) - U_c(el) + U_c(elem_neighbors(3,el))
      UV2_rhs(2,el) = UV2_rhs(2,el) - V_c(el) + V_c(elem_neighbors(3,el))
   endif
   if ( elem_neighbors(4,el) > 0  .and. elem_neighbors(4,el)/=elem_neighbors(1,el)) then
      UV2_rhs(1,el) = UV2_rhs(1,el) - U_c(el) + U_c(elem_neighbors(4,el))
      UV2_rhs(2,el) = UV2_rhs(2,el) - V_c(el) + V_c(elem_neighbors(4,el))
   endif

ENDDO
!$OMP END DO
!$OMP END PARALLEL

#ifdef USE_MPI
  call exchange_elem(UV2_rhs) !SH Size of the array correct order?
#endif


end subroutine viscosity_filt2x_2D
!=======================================================================
SUBROUTINE viscosity_filt_2D

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  use g_parsup
  use g_comm_auto

  IMPLICIT NONE

  real(kind=WP)  :: u1(2), tau_inv
  integer       :: elem, nelem, j, q, elnodes(4), jend

  ! subroutine applies filter to velocity, which on uniform
  ! equilateral mesh will be equivalent to the Laplacian operator
  ! It adds to the rhs(0) tau_filt*(u1*h1+u2*h2+u3*h3-3*u0*h0)/3
  ! The contribution from boundary edges is neglected, it can
  ! be done in function of boundary conditions
!now from fv_run
  tau_inv=2.0_WP/3600.0_WP*2.0_WP     ! it is 0.5/day now SET IT experimentally
  tau_inv= tau_inv_filt
  !if tau_inv=3 Visc/h^2, the result should be roughly equivalent
  !to ordinary viscosity. tau_inv=3*Visc/2/elem_area, and multiplication
  !with the element area furthercan be removed! TO DO if needed

  Do elem=1,myDim_elem2D
     elnodes=elem2D_nodes(:,elem)
     u1=0.0_WP
     jend = 4
     if(elnodes(1)==elnodes(4)) jend = 3
     Do j=1,jend
        nelem=elem_neighbors(j,elem)
        if(nelem>0) u1=u1+UAB(:, nelem) -UAB(:,elem)
     end do
     elnodes=elem2D_nodes(:,elem)
     U_rhs_2D(:,elem)=U_rhs_2D(:,elem)+tau_inv*u1*elem_area(elem) &
          *sum(w_cv(1:4,elem)*depth(elnodes))
  end do

#ifdef USE_MPI
  call exchange_elem(U_rhs_2D)
#endif

end subroutine viscosity_filt_2D
! ===================================================================
SUBROUTINE viscosity_filt_3D

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  use g_parsup
  use g_comm_auto

  IMPLICIT NONE

  real(kind=WP)  :: u1, v1, tau_inv
  integer          :: elem, nelem, j, nz, jend, elnodes(4)

  ! subroutine applies filter to velocity, which on uniform
  ! equilateral mesh will be equivalent to the Laplacian operator
  ! It adds to the rhs(0) tau_filt*(u1*h1+u2*h2+u3*h3-3*u0*h0)/3
  ! The contribution from boundary edges is neglected, it can
  ! be done in function of boundary conditions
!now from fv_run
  tau_inv=2.0_WP/3600.0_WP*2.0_WP   ! it is 0.5/day now SET IT experimentally
  tau_inv=tau_inv_filt
  !if tau_inv=3 Visc/h^2, the result should be roughly equivalent
  !to ordinary viscosity. tau_inv=3*Visc/2/elem_area, and multiplication
  !with the element area furthercan be removed! TO DO if needed

!$OMP PARALLEL DO PRIVATE(elem,elnodes,jend,nz,nelem,u1,v1)
  Do elem=1,myDim_elem2D
     elnodes=elem2D_nodes(:,elem)
     jend = 4
     if(elnodes(1)==elnodes(4)) jend = 3

     do nz=1,nsigma-1
        u1 = 0.0_WP
        v1 = 0.0_WP
        Do j=1,jend
           nelem=elem_neighbors(j,elem)
           if(nelem>0) then
              u1 = u1+U_n(nz, nelem) -U_n(nz,elem)
              v1 = v1+V_n(nz, nelem) -V_n(nz,elem)
           endif
        end do
        U_rhs(nz,elem) = U_rhs(nz,elem)+tau_inv*u1*elem_area(elem)*Je(nz,elem)
        V_rhs(nz,elem) = V_rhs(nz,elem)+tau_inv*v1*elem_area(elem)*Je(nz,elem)
     end do
  end do
!$OMP END PARALLEL DO

#ifdef USE_MPI
  call exchange_elem(U_rhs)
  call exchange_elem(V_rhs)
#endif

end subroutine viscosity_filt_3D
! ===================================================================
SUBROUTINE viscosity_filt_3D_to_2D

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  use g_parsup
  use g_comm_auto

  IMPLICIT NONE

  real(kind=WP)  ::  u1, v1, tau_inv
  real(kind=WP), allocatable  ::  U_3Dfiltr_puls(:,:), V_3Dfiltr_puls(:,:)
  integer          :: elem, nelem, j, elnodes(4), nz, jend, elemsize

  allocate(U_3Dfiltr_puls(nsigma-1,myDim_elem2D), V_3Dfiltr_puls(nsigma-1,myDim_elem2D) )

  U_3Dfiltr_puls = 0.0_WP
  V_3Dfiltr_puls = 0.0_WP

  ! subroutine applies filter to velocity, which on uniform
  ! equilateral mesh will be equivalent to the Laplacian operator
  ! It adds to the rhs(0) tau_filt*(u1*h1+u2*h2+u3*h3-3*u0*h0)/3
  ! The contribution from boundary edges is neglected, it can
  ! be done in function of boundary conditions
!now from fv_run
  tau_inv=2.0_WP/3600.0_WP*2.0_WP      ! it is 1/day now SET IT experimentally
  tau_inv=tau_inv_filt
  !if tau_inv=3 Visc/h^2, the result should be roughly equivalent
  !to ordinary viscosity. tau_inv=3*Visc/2/elem_area, and multiplication
  !with the element area furthercan be removed! TO DO if needed
  Do elem=1,myDim_elem2D
     elnodes=elem2D_nodes(:,elem)
     jend = 4
     if(elnodes(1)==elnodes(4)) jend = 3
     do nz=1,nsigma-1
        u1 = 0.0_WP
        v1 = 0.0_WP
        Do j=1,jend
           nelem=elem_neighbors(j,elem)
           if (nelem>0) then
              u1=u1+U_puls(nz, nelem)*Je(nz,nelem) -U_puls(nz,elem)*Je(nz,elem)
              v1=v1+V_puls(nz, nelem)*Je(nz,nelem) -V_puls(nz,elem)*Je(nz,elem)
           endif
        end do
        U_3Dfiltr_puls(nz,elem)=tau_inv*u1*elem_area(elem)
        V_3Dfiltr_puls(nz,elem)=tau_inv*v1*elem_area(elem)
     end do
  end do

! No exchange needed!
!#ifdef USE_MPI
!  call exchange_elem(U_3Dfiltr_puls)
!  call exchange_elem(V_3Dfiltr_puls)
!#endif

  DO elem=1,myDim_elem2D
     DO nz=1,nsigma-1
        U_rhs_2D_3D(1,elem) = U_rhs_2D_3D(1,elem) + U_3Dfiltr_puls(nz,elem)*Je(nz,elem)
        U_rhs_2D_3D(2,elem) = U_rhs_2D_3D(2,elem) + V_3Dfiltr_puls(nz,elem)*Je(nz,elem)
     END DO
  END DO

  deallocate(U_3Dfiltr_puls, V_3Dfiltr_puls)

#ifdef USE_MPI
  call exchange_elem(U_rhs_2D_3D)
#endif

end subroutine viscosity_filt_3D_to_2D

! ===================================================================
!Biharmonic viscosity routines.
! Contains: biharmonic_viscosity (only biharmonic operator)
!           biPlusharmonic_viscosity (harmonic+biharmonic)
!           vel_lapl_gradiens ==vel_gradients, but simplified (no account for free slip).
!
! True viscosity operator is implemented in laplacian_viscosity (div sigma_ij,
! with sigma_ij=2*A_hor(epsilon_ij-delta_ij*div u/2). If we apply it twice, the resultant
! biharmonic operator does not kill noise if eddies are strong. We use operator
! Lu=nabla nabla u, corrected as suggested by Blazek, 2001, (the correction concerns the
! way nabla u is computed on cell-centered discretization) and apply it twice:
! bh_visc=-L Abh Lu.
!
! If Laplacian and biharmonic are used simultaneously, we assemble
! lapl_visc=nabla A_hor nabla u and bh_visc in the same subroutine to
! increase numerical efficiency. Note, however, that operators introduced in that way are
! not zero under solid body rotations, while the operator in laplacian_viscosity is.
!
! Note also that we ignore metric differentiation terms because our viscosity
! is just a parameterization. There is no problem with including them.
!
! sergey.danilov@awi.de  2012
! With Laplacian regularization
SUBROUTINE biharmonic_viscosity_3D_to_2D
!
!a changes have been made by Alexey Androsov
!a (for sigma coordinates, and for tri-quad elements)
!a 09.10.14
!
USE o_MESH
USE o_ARRAYS
USE o_PARAM

  use g_parsup
  use g_comm_auto

IMPLICIT NONE
 real(kind=WP)          :: xe, ye, u1, v1, Ah, acc !, delta_sigma
 real(kind=WP)          :: tx, ty, tl, tt, g1, g2
 real(kind=WP), allocatable  :: U_c(:,:), V_c(:,:)
 real(kind=WP), allocatable  :: U_3Dbih_puls(:,:), V_3Dbih_puls(:,:)
 integer               :: ed, el(2), nz, elem

  integer :: edglim

  allocate(U_c(nsigma-1,myDim_elem2D+eDim_elem2D+eXDim_elem2D), &
           V_c(nsigma-1, myDim_elem2D+eDim_elem2D+eXDim_elem2D))
  allocate(U_3Dbih_puls(nsigma-1, myDim_elem2D+eDim_elem2D+eXDim_elem2D), &
           V_3Dbih_puls(nsigma-1, myDim_elem2D+eDim_elem2D+eXDim_elem2D))

U_3Dbih_puls = 0.0_WP
V_3Dbih_puls = 0.0_WP

 call vel_gradients_puls

   U_c=0.0_WP
   V_c=0.0_WP

#ifdef USE_MPI
  edglim=myDim_edge2D+eDim_edge2D
#else
  edglim=edge2D_in
#endif

  DO ed=1, edglim

#ifdef USE_MPI
     if (myList_edge2D(ed)>edge2D_in) cycle
#endif

   el=edge_tri(:,ed)
   xe=edge_dxdy(1,ed)*r_earth*0.5_WP*(elem_cos(el(1))+elem_cos(el(2)))
   ye=edge_dxdy(2,ed)*r_earth

   tx=-edge_cross_dxdy(1,ed)+edge_cross_dxdy(3,ed)
   ty=-edge_cross_dxdy(2,ed)+edge_cross_dxdy(4,ed)
   tl=sqrt(tx**2+ty**2)
   tx=tx/tl
   ty=ty/tl

   DO nz=1, nsigma-1
   g1=0.5_WP*(vel_grad_ux(nz,el(1))+vel_grad_ux(nz,el(2)))
   g2=0.5_WP*(vel_grad_uy(nz,el(1))+vel_grad_uy(nz,el(2)))
   tt=g1*tx+g2*ty-(U_puls(nz,el(2))-U_puls(nz,el(1)))/tl
   g1=g1-tx*tt
   g2=g2-ty*tt
   u1=(g1*ye-g2*xe)
   g1=0.5_WP*(vel_grad_vx(nz,el(1))+vel_grad_vx(nz,el(2)))
   g2=0.5_WP*(vel_grad_vy(nz,el(1))+vel_grad_vy(nz,el(2)))
   tt=g1*tx+g2*ty-(V_puls(nz,el(2))-V_puls(nz,el(1)))/tl
   g1=g1-tx*tt
   g2=g2-ty*tt
   v1=(g1*ye-g2*xe)

   U_c(nz,el(1))=U_c(nz,el(1))+u1*Jd(nz,ed)
   V_c(nz,el(1))=V_c(nz,el(1))+v1*Jd(nz,ed)
   U_c(nz,el(2))=U_c(nz,el(2))-u1*Jd(nz,ed)
   V_c(nz,el(2))=V_c(nz,el(2))-v1*Jd(nz,ed)

   END DO
END DO

 ! ============
 ! Contribution from boundary edges
 ! ============


#ifdef USE_MPI

  call exchange_elem(U_c) !SH check this
  call exchange_elem(V_c) !SH check this

  DO ed=1,edglim
     if (myList_edge2D(ed)>edge2D_in) then

        el=edge_tri(:,ed)
        xe=edge_dxdy(1,ed)*r_earth*elem_cos(el(1))
        ye=edge_dxdy(2,ed)*r_earth

        DO nz=1,nsigma-1
           u1=(vel_grad_ux(nz,el(1)))*ye-(vel_grad_uy(nz,el(1)))*xe
           v1=(vel_grad_vx(nz,el(1)))*ye-(vel_grad_vy(nz,el(1)))*xe
           U_c(nz,el(1))=U_c(nz,el(1))+u1*Jd(nz,ed)
           V_c(nz,el(1))=V_c(nz,el(1))+v1*Jd(nz,ed)
        END DO

     end if
  END DO

  call exchange_elem(U_c)
  call exchange_elem(V_c)

#else

 DO ed=1+edge2d_in,edge2D
   el=edge_tri(:,ed)
   xe=edge_dxdy(1,ed)*r_earth*elem_cos(el(1))
   ye=edge_dxdy(2,ed)*r_earth

   DO nz=1,nsigma-1
   u1=(vel_grad_ux(nz,el(1)))*ye-(vel_grad_uy(nz,el(1)))*xe
   v1=(vel_grad_vx(nz,el(1)))*ye-(vel_grad_vy(nz,el(1)))*xe
   U_c(nz,el(1))=U_c(nz,el(1))+u1*Jd(nz,ed)
   V_c(nz,el(1))=V_c(nz,el(1))+v1*Jd(nz,ed)
   END DO
END DO

#endif

 ! =============
 ! Multiply with biharmonic viscosity
 ! =============

  DO elem=1, myDim_elem2D
 ! we have to divide on elem_area
 ! - takes into account the fact that
 ! biharmonic goes with different sign.
         acc = sum(w_cv(1:4,elem)*ac(elem2D_nodes(:,elem)))
   Ah=-acc*Abh0*sqrt(elem_area(elem)/scale_area)/scale_area
   DO nz=1, nsigma-1
     U_c(nz,elem)=U_c(nz,elem)*(Ah-Visc(nz,elem))
     V_c(nz,elem)=V_c(nz,elem)*(Ah-Visc(nz,elem))
   END DO
 END DO

#ifdef USE_MPI
  call exchange_elem(U_c)
  call exchange_elem(V_c)
#endif

 call vel_lapl_gradients(U_c, V_c)


 ! =============
 ! Apply Laplace operator once more
 ! =============

  DO ed=1, edglim

#ifdef USE_MPI
     if (myList_edge2D(ed)>edge2D_in) cycle
#endif

   el=edge_tri(:,ed)

   xe=edge_dxdy(1,ed)*r_earth*0.5_WP*(elem_cos(el(1))+elem_cos(el(2)))
   ye=edge_dxdy(2,ed)*r_earth

   tx=-edge_cross_dxdy(1,ed)+edge_cross_dxdy(3,ed)
   ty=-edge_cross_dxdy(2,ed)+edge_cross_dxdy(4,ed)
   tl=sqrt(tx**2+ty**2)
   tx=tx/tl
   ty=ty/tl

   DO nz=1, nsigma-1
   g1=0.5_WP*(vel_grad_ux(nz,el(1))+vel_grad_ux(nz,el(2)))
   g2=0.5_WP*(vel_grad_uy(nz,el(1))+vel_grad_uy(nz,el(2)))
   tt=g1*tx+g2*ty-(U_c(nz,el(2))-U_c(nz,el(1)))/tl
   g1=g1-tx*tt
   g2=g2-ty*tt
   u1=(g1*ye-g2*xe)
   g1=0.5_WP*(vel_grad_vx(nz,el(1))+vel_grad_vx(nz,el(2)))
   g2=0.5_WP*(vel_grad_vy(nz,el(1))+vel_grad_vy(nz,el(2)))
   tt=g1*tx+g2*ty-(V_c(nz,el(2))-V_c(nz,el(1)))/tl
   g1=g1-tx*tt
   g2=g2-ty*tt
   v1=(g1*ye-g2*xe)

   U_3Dbih_puls(nz,el(1))=U_3Dbih_puls(nz,el(1))+u1 !*Jd(nz,ed)
   V_3Dbih_puls(nz,el(1))=V_3Dbih_puls(nz,el(1))+v1 !*Jd(nz,ed)
   U_3Dbih_puls(nz,el(2))=U_3Dbih_puls(nz,el(2))-u1 !*Jd(nz,ed)
   V_3Dbih_puls(nz,el(2))=V_3Dbih_puls(nz,el(2))-v1 !*Jd(nz,ed)

   END DO
END DO

! ============
! Contribution from boundary edges
! ============

#ifdef USE_MPI

  call exchange_elem(U_3Dbih_puls) !SH Check this
  call exchange_elem(V_3Dbih_puls) !SH Check this

  DO ed=1,edglim
     if (myList_edge2D(ed)>edge2D_in) then

        el=edge_tri(:,ed)
        xe=edge_dxdy(1,ed)*r_earth*elem_cos(el(1))
        ye=edge_dxdy(2,ed)*r_earth

        DO nz=1, nsigma-1
           u1=(vel_grad_ux(nz,el(1)))*ye-(vel_grad_uy(nz,el(1)))*xe
           v1=(vel_grad_vx(nz,el(1)))*ye-(vel_grad_vy(nz,el(1)))*xe
           U_3Dbih_puls(nz,el(1))=U_3Dbih_puls(nz,el(1))+u1 !*Jd(nz,ed)
           V_3Dbih_puls(nz,el(1))=V_3Dbih_puls(nz,el(1))+v1 !*Jd(nz,ed)
        END DO

     end if

  END DO

  call exchange_elem(U_3Dbih_puls)
  call exchange_elem(V_3Dbih_puls)

#else

 DO ed=1+edge2D_in, edge2D

   el=edge_tri(:,ed)
   xe=edge_dxdy(1,ed)*r_earth*elem_cos(el(1))
   ye=edge_dxdy(2,ed)*r_earth

   DO nz=1, nsigma-1
   u1=(vel_grad_ux(nz,el(1)))*ye-(vel_grad_uy(nz,el(1)))*xe
   v1=(vel_grad_vx(nz,el(1)))*ye-(vel_grad_vy(nz,el(1)))*xe
   U_3Dbih_puls(nz,el(1))=U_3Dbih_puls(nz,el(1))+u1 !*Jd(nz,ed)
   V_3Dbih_puls(nz,el(1))=V_3Dbih_puls(nz,el(1))+v1 !*Jd(nz,ed)
   END DO
END DO

#endif

  DO elem=1,myDim_elem2D
   DO nz=1,nsigma-1
 !  delta_sigma = sigma(nz) - sigma(nz+1)
   U_rhs_2D_3D(1,elem) = U_rhs_2D_3D(1,elem) + U_3Dbih_puls(nz,elem)*Je(nz,elem) !delta_sigma
   U_rhs_2D_3D(2,elem) = U_rhs_2D_3D(2,elem) + V_3Dbih_puls(nz,elem)*Je(nz,elem) !delta_sigma
    END DO
  END DO

#ifdef USE_MPI
  call exchange_elem(U_rhs_2D_3D)
#endif

deallocate(U_3Dbih_puls, V_3Dbih_puls)
deallocate(V_c, U_c)
END subroutine biharmonic_viscosity_3D_to_2D
! ===================================================================================
SUBROUTINE momentum_vert_expl_visc
!
! Vertical viscosity, explicit
! 24.11.15
!
USE o_PARAM
USE o_MESH
USE o_ARRAYS
IMPLICIT NONE

integer               :: elem, elnodes(4), nz
real(kind=WP)    :: friction, a, b, zinv, a_Cd, acc
real(kind=WP)    ::  uvert(2,nsigma)

 DO elem=1,elem2D
    elnodes=elem2D_nodes(:,elem)
    acc = sum(w_cv(1:4,elem)*ac(elnodes))
    !VF diff. bottom friction options are used
    a_Cd = C_d_el(elem)*((1.0_WP - acc)*Cd_bz + 1.0_WP)

   DO nz=2, nsigma-1
      a = sum(w_cv(1:4,elem)*Z(nz,elnodes))
      b = sum(w_cv(1:4,elem)*Z(nz-1,elnodes))
      zinv=Av(nz,elem)*elem_area(elem)/(a-b)
      uvert(1,nz)= (U_n(nz-1,elem)-U_n(nz,elem))*zinv
      uvert(2,nz)= (V_n(nz-1,elem)-V_n(nz,elem))*zinv
   END DO
   ! Wind stress and bottom drag
     uvert(1,1)= taux(elem)*elem_area(elem)/density_0
     uvert(2,1)= tauy(elem)*elem_area(elem)/density_0
     friction=a_Cd*sqrt(U_n(nsigma-1,elem)**2+ &
                          V_n(nsigma-1,elem)**2)*elem_area(elem)

     uvert(1,nsigma)=friction*U_n(nsigma-1,elem)
     uvert(2,nsigma)=friction*V_n(nsigma-1,elem)
     ! + sign here because it is subtracted!
   DO nz=1,nsigma-1
      U_rhs(nz,elem)=U_rhs(nz,elem)+(uvert(1,nz)-uvert(1,nz+1))
      V_rhs(nz,elem)=V_rhs(nz,elem)+(uvert(2,nz)-uvert(2,nz+1))
   END DO
 END DO

END SUBROUTINE momentum_vert_expl_visc
!===========================================================================

subroutine momentum_vert_impl_visc

  USE o_MESH
  USE o_PARAM
  USE o_ARRAYS

  use g_parsup
  use g_comm_auto

  IMPLICIT NONE

  real(kind=WP)              ::  a(nsigma-1), b(nsigma-1), c(nsigma-1)
  real(kind=WP)              ::  ur(nsigma-1), vr(nsigma-1)
  real(kind=WP)              ::  cp(nsigma-1), up(nsigma-1), vp(nsigma-1)
  real(kind=WP)              ::  Z_c, Z_u, Z_b
  real(kind=WP)              ::  zinv, m, friction, acc, a_Cd
  integer                    ::  elnodes(4)
  integer                    ::  nz, elem

!$OMP DO
  DO elem=1,myDim_elem2D

     elnodes=elem2D_nodes(:,elem)
     acc = sum(w_cv(1:4,elem)*ac(elnodes))
     !VF diff. bottom friction options are used
     a_Cd = C_d_el(elem)*((1.0_WP - acc)*Cd_bz + 1.0_WP)

     ! ==========================
     ! Operator
     ! ==========================
     ! Regular part of coefficients:
     DO nz=2, nsigma-2
        zinv = dt/Je(nz,elem)
        Z_c = sum(w_cv(1:4,elem)*Z(nz,elnodes))
        Z_u = sum(w_cv(1:4,elem)*Z(nz-1,elnodes))
        Z_b = sum(w_cv(1:4,elem)*Z(nz+1,elnodes))
        a(nz)=-Av(nz,elem)/(Z_c-Z_u)*zinv
        c(nz)=-Av(nz+1,elem)/(Z_b-Z_c)*zinv
        b(nz)=-a(nz)-c(nz)+1.0_WP
     END DO
     ! The last row
     zinv=dt/Je(nsigma-1,elem)
     Z_c = sum(w_cv(1:4,elem)*Z(nsigma-1,elnodes))
     Z_u = sum(w_cv(1:4,elem)*Z(nsigma-2,elnodes))
     a(nsigma-1)=-Av(nsigma-1,elem)/(Z_c-Z_u)*zinv
     b(nsigma-1)=-a(nsigma-1)+1.0_WP
     c(nsigma-1)=0.0_WP
     ! The first row
     zinv=dt/Je(1,elem)
     Z_c = sum(w_cv(1:4,elem)*Z(1,elnodes))
     Z_b = sum(w_cv(1:4,elem)*Z(2,elnodes))
     c(1)=-Av(2,elem)/(Z_b-Z_c)*zinv
     a(1)=0.0_WP
     b(1)=-c(1)+1.0_WP
     ! ===========================
     ! The rhs:
     ! ===========================
     ur(1:nsigma-1)=U_rhs(1:nsigma-1,elem)   ! no dt here, it is already in RHSs
     vr(1:nsigma-1)=V_rhs(1:nsigma-1,elem)
     ! The first row contains surface forcing
     ur(1)= ur(1) +zinv*taux(elem)/density_0
     vr(1)= vr(1) +zinv*tauy(elem)/density_0
     ! The last row contains bottom friction
     zinv=dt/Je(nsigma-1,elem)
     friction=-a_Cd*sqrt(U_n(nsigma-1,elem)**2+ &
          V_n(nsigma-1,elem)**2)
     ur(nsigma-1)=ur(nsigma-1) +zinv*friction*U_n(nsigma-1,elem)
     vr(nsigma-1)=vr(nsigma-1) +zinv*friction*V_n(nsigma-1,elem)
     ! ===========================
     ! The sweep algorithm
     ! ===========================
     ! initialize c-prime and s,t-prime
     cp(1) = c(1)/b(1)
     up(1) = ur(1)/b(1)
     vp(1) = vr(1)/b(1)
     ! solve for vectors c-prime and t, s-prime
     do nz = 2,nsigma-1
        m = b(nz)-cp(nz-1)*a(nz)
        cp(nz) = c(nz)/m
        up(nz) = (ur(nz)-up(nz-1)*a(nz))/m
        vp(nz) = (vr(nz)-vp(nz-1)*a(nz))/m
     enddo
     ! initialize x
     ur(nsigma-1) = up(nsigma-1)
     vr(nsigma-1) = vp(nsigma-1)
     ! solve for x from the vectors c-prime and d-prime
     do nz = nsigma-2, 1, -1
        ur(nz) = up(nz)-cp(nz)*ur(nz+1)
        vr(nz) = vp(nz)-cp(nz)*vr(nz+1)
     end do
     ! ===========================
     ! RHS update
     ! ===========================
     DO nz=1,nsigma-1
        U_rhs(nz,elem)=ur(nz)*mask_wd(elem)
        V_rhs(nz,elem)=vr(nz)*mask_wd(elem)
     END DO
  END DO   !!! cycle over elements
  !$OMP END DO

#ifdef USE_MPI
  call exchange_elem(U_rhs)
  call exchange_elem(V_rhs)
#endif

END subroutine momentum_vert_impl_visc

!=======================================================================
subroutine momentum_vert_impl_visc_tmp
    USE o_MESH
    USE o_ARRAYS
    USE o_PARAM
    USE i_ARRAYS
    USE i_PARAM


    use g_parsup
    use g_comm_auto
IMPLICIT NONE
    integer   :: elem, nz, elnodes(4)
    real(kind=WP) :: Fx, Fy, delta_1, delta_2, delta, A_el
    real(kind=WP) :: zt, a, b, c, d, de, friction, rw, acc, a_Cd
    real(kind=WP), allocatable  ::  xx1(:), yy1_u(:), yy1_v(:)
    real(kind=WP), allocatable  ::  U_tmp(:), V_tmp(:)

    allocate(xx1(nsigma),yy1_u(nsigma), yy1_v(nsigma))
    allocate(U_tmp(nsigma), V_tmp(nsigma))

    U_tmp=0.0_WP
    V_tmp=0.0_WP
    xx1=0.0_WP
    yy1_u=0.0_WP
    yy1_v=0.0_WP

    do elem=1,myDim_elem2D
        if (mask_wd(elem) /= 0.0_WP) then
            elnodes=elem2D_nodes(:,elem)
            xx1(1) = 1.0_WP
            de = Je(1,elem)/(Av(1,elem) + Av(2,elem))
            rw = de/density_0
            if (use_ice) then
                A_el =   w_cv(1,elem) *a_ice(elem2D_nodes(1,elem)) + w_cv(2,elem) *a_ice(elem2D_nodes(2,elem)) &
                        + w_cv(3,elem) *a_ice(elem2D_nodes(3,elem)) + w_cv(4,elem) *a_ice(elem2D_nodes(4,elem))

                yy1_u(1) =  rw*tau_ice(2,elem)*(U_n_ice(1,elem) - U_n(1,elem))* A_el + taux(elem)*rw * (1.0_WP - A_el)
                yy1_v(1) =  rw*tau_ice(2,elem)*(U_n_ice(2,elem) - V_n(1,elem))* A_el + taux(elem)*rw * (1.0_WP - A_el)
            else
                yy1_u(1) =  taux(elem)*rw
                yy1_v(1) =  tauy(elem)*rw
            endif

            do nz=1,nsigma-1
                Fx = U_rhs(nz,elem)
                Fy = V_rhs(nz,elem)
                if (nz == 1) then
                    delta_1 = 0.5_WP*Je(nz,elem)
                else
                    delta_1 = 0.5_WP*(Je(nz,elem) + Je(nz-1,elem))
                endif
                if (nz == nsigma-1) then
                    delta_2= 0.5_WP*Je(nz,elem)
                else
                    delta_2 = 0.5_WP*(Je(nz,elem) + Je(nz+1,elem))
                endif
                zt = dt/Je(nz,elem)
                a = zt*Av(nz,elem)/delta_1
                c = zt*Av(nz+1,elem)/delta_2
                b = a + c + 1.0_WP
                d = b - a*xx1(nz)
                xx1(nz+1) = c/d
                yy1_u(nz+1) = (a*yy1_u(nz) + Fx)/d
                yy1_v(nz+1) = (a*yy1_v(nz) + Fy)/d
            enddo
    !+++++++++++++++++++++++++++++++
    ! bottom boundary conditions
    !+++++++++++++++++++++++++++++++

            acc = sum(w_cv(1:4,elem)*ac(elnodes))
            a_Cd = C_d_el(elem)*((1.0_WP - acc)*Cd_bz + 1.0_WP)

            de = Je(nsigma-1,elem)/(Av(nsigma,elem) + Av(nsigma-1,elem))
            friction=de*a_Cd*sqrt(U_n(nsigma-1,elem)**2+V_n(nsigma-1,elem)**2)

            if (fic_point == 0) then
                U_tmp(nsigma) = 0.0_WP
                V_tmp(nsigma) = 0.0_WP
            elseif (fic_point == 1) then
                a = friction + 1.0_WP - xx1(nsigma)
                U_tmp(nsigma) =yy1_u(nsigma)/a
                V_tmp(nsigma) =yy1_v(nsigma)/a
            elseif (fic_point == 2) then
                a = 1.0_WP/(1.0_WP - friction) - xx1(nsigma)
                U_tmp(nsigma) =yy1_u(nsigma)/a
                V_tmp(nsigma) =yy1_v(nsigma)/a
            endif

            do nz=nsigma-1,1,-1
                U_tmp(nz) = U_tmp(nz+1)*xx1(nz+1) + yy1_u(nz+1)
                V_tmp(nz) = V_tmp(nz+1)*xx1(nz+1) + yy1_v(nz+1)
            enddo

            do nz=1,nsigma-1
                U_rhs(nz,elem) = U_tmp(nz)*mask_wd(elem)
                V_rhs(nz,elem) = V_tmp(nz)*mask_wd(elem)
            enddo

        else
            do nz=1,nsigma-1
                U_rhs(nz,elem) = 0.0_WP
                V_rhs(nz,elem) = 0.0_WP
            enddo
        endif

    enddo
#ifdef USE_MPI
    call exchange_elem(U_rhs)
    call exchange_elem(V_rhs)
#endif
    deallocate(xx1, yy1_u, yy1_v, U_tmp, V_tmp)

end subroutine momentum_vert_impl_visc_tmp
!
!=======================================================================
!aa subroutine momentum_vert_impl_visc_tmp
!aa USE o_MESH
!aa USE o_ARRAYS
!aa USE o_PARAM

!aa IMPLICIT NONE
!aa integer   :: elem, nz, elnodes(4)
!aa real(kind=WP) :: Fx, Fy, delta_1, delta_2, delta
!aa real(kind=WP) :: zt, a, b, c, d, de, friction, rw, acc, a_Cd
!aa real(kind=WP) :: xx1(nsigma),yy1_u(nsigma), yy1_v(nsigma)
!aa real(kind=WP) ::  U_tmp(nsigma), V_tmp(nsigma)

!aa U_tmp=0.0_WP
!aa V_tmp=0.0_WP
!aa xx1=0.0_WP
!aa yy1_u=0.0_WP
!aa yy1_v=0.0_WP

!aa!$OMP DO
!aa  DO elem=1,elem2D

!aa    elnodes=elem2D_nodes(:,elem)

!aa    xx1(1) = 1.0_WP
!aa    rw = Je(1,elem)/((Av(1,elem) + Av(2,elem))*density_0)
!aa    yy1_u(1) = taux(elem)*rw
!aa    yy1_v(1) = tauy(elem)*rw

    ! nz = 1

!aa    a = dt*Av(1,elem)/(0.5_WP*Je(1,elem)*Je(1,elem))
!aa    c = dt*Av(2,elem)/(0.5_WP*(Je(1,elem) + Je(2,elem))*Je(1,elem))
!aa    b = a + c + 1.0_WP
!aa    d = 1._WP/(b - a)
!aa    xx1(2) = c*d
!aa    yy1_u(2) = (a*yy1_u(1) + U_rhs(1,elem))*d
!aa    yy1_v(2) = (a*yy1_v(1) + V_rhs(1,elem))*d

!aa    DO nz=2,nsigma-2
!aa      a = dt*Av(nz  ,elem)/(0.5_WP*(Je(nz,elem) + Je(nz-1,elem))*Je(nz,elem))
!aa      c = dt*Av(nz+1,elem)/(0.5_WP*(Je(nz,elem) + Je(nz+1,elem))*Je(nz,elem))
!aa      b = a + c + 1.0_WP
!aa      d = 1._WP/(b - a*xx1(nz))
!aa      xx1(nz+1) = c*d
!aa      yy1_u(nz+1) = (a*yy1_u(nz) + U_rhs(nz,elem))*d
!aa      yy1_v(nz+1) = (a*yy1_v(nz) + V_rhs(nz,elem))*d

!aa    enddo

    ! nz = nsigma-1
!aa    a = dt*Av(nsigma-1,elem)/(0.5_WP*(Je(nsigma-1,elem) + Je(nsigma-2,elem))*Je(nsigma-1,elem))
!aa    c = dt*Av(nsigma,  elem)/(0.5_WP*Je(nsigma-1,elem)*Je(nsigma-1,elem))
!aa    b = a + c + 1.0_WP
!aa    d = 1._WP/(b - a*xx1(nsigma-1))
!aa    xx1(nsigma)   = c*d
!aa    yy1_u(nsigma) = (a*yy1_u(nsigma-1) + U_rhs(nsigma-1,elem))*d
!aa    yy1_v(nsigma) = (a*yy1_v(nsigma-1) + V_rhs(nsigma-1,elem))*d

!+++++++++++++++++++++++++++++++
! bottom boundary conditions
!+++++++++++++++++++++++++++++++
!aa    acc = sum(w_cv(1:4,elem)*ac(elnodes))
!aa    a_Cd = C_d*((1.0_WP - acc)*Cd_bz + 1.0_WP)
!aa    de = Je(nsigma-1,elem)/(Av(nsigma,elem) + Av(nsigma-1,elem))
!aa    friction=de*a_Cd*sqrt(U_n(nsigma-1,elem)**2+V_n(nsigma-1,elem)**2)

!aa    if (fic_point == 0) then
!aa       U_tmp(nsigma) = 0.0_WP
!aa       V_tmp(nsigma) = 0.0_WP
!aa    elseif (fic_point == 1) then
!aa       a = 1._WP/(friction + 1.0_WP - xx1(nsigma))
!aa       U_tmp(nsigma) =yy1_u(nsigma)*a
!aa       V_tmp(nsigma) =yy1_v(nsigma)*a
!aa    elseif (fic_point == 2) then
!       a = 1._WP/(1.0_WP/(1.0_WP - friction) - xx1(nsigma))
! should be the same as: (only one division)
!aa        a = (1.0_WP - friction)/(1.0_WP - (1.0_WP - friction)*xx1(nsigma))
!aa       U_tmp(nsigma) =yy1_u(nsigma)*a
!aa       V_tmp(nsigma) =yy1_v(nsigma)*a
!aa    endif

!aa    do nz=nsigma-1,1,-1
!aa     U_tmp(nz) = U_tmp(nz+1)*xx1(nz+1) + yy1_u(nz+1)
!aa     V_tmp(nz) = V_tmp(nz+1)*xx1(nz+1) + yy1_v(nz+1)
!aa    enddo

!aa    do nz=1,nsigma-1
!aa     U_rhs(nz,elem) = U_tmp(nz)*mask_wd(elem)
!aa     V_rhs(nz,elem) = V_tmp(nz)*mask_wd(elem)
!aa    enddo

!aa  enddo
!aa !$OMP END DO

!aa end subroutine momentum_vert_impl_visc_tmp
!=======================================================================



!
!=======================================================================
! alexey.androsov@awi.de
! 06.11.18 (original date), introduced to MPI version at 04.11.20 (IK))
! correction bottom BC (MOST)

subroutine vert_impl_MOST
    use o_MESH
    use o_ARRAYS
    use o_PARAM

    use g_parsup
    use g_comm_auto
IMPLICIT NONE
    integer   :: elem, nz, elnodes(4)
    real(kind=WP) :: Sx, Sy, lambda_2, lambda_3, am = 5.0_WP, FX, Fy
    real(kind=WP) :: Z_1, Z_2, Z_3, zx_MOST, zy_MOST, u_st, v_st, bu_st, bv_st
    real(kind=WP) :: delta_1, delta_2, delta, Lu_st, Lv_st, dmean
    real(kind=WP) :: zt, a, b, c, d, de, friction, rw, acc, a_Cd, zx_MOST_max, Zy_MOST_max
    real(kind=WP) :: zx_MOST_min, zy_MOST_min, z0b_new_x, z0b_new_y
    real(kind=WP), allocatable  ::  xx1(:), yy1_u(:), yy1_v(:)
    real(kind=WP), allocatable  ::  U_tmp(:), V_tmp(:)


    allocate(xx1(nsigma),yy1_u(nsigma), yy1_v(nsigma))
    allocate(U_tmp(nsigma), V_tmp(nsigma))

    U_tmp=0.0_WP
    V_tmp=0.0_WP
    xx1=0.0_WP
    yy1_u=0.0_WP
    yy1_v=0.0_WP

    zx_MOST_max = -999999.0_WP
    zx_MOST_min = 9999999.0_WP
    zy_MOST_max = -999999.0_WP
    zy_MOST_min = 9999999.0_WP


    do elem=1,myDim_elem2D
        if (mask_wd(elem) /= 0.0_WP) then
            elnodes=elem2D_nodes(:,elem)
            acc = sum(w_cv(1:4,elem)*ac(elnodes))
            if (acc .gt. 0.9_WP) then
                dmean = max(Dmin,sum(w_cv(1:4,elem)*(depth(elem2D_nodes(:,elem)) + eta_n(elem2D_nodes(:,elem)))))
                a = U_rhs(nsigma-2,elem) - U_rhs(nsigma-1,elem)
                if (a == 0.0_WP) then
                    Sx = 0.0_WP
                else
                    Sx = (U_rhs(nsigma-3,elem) - U_rhs(nsigma-1,elem))/(U_rhs(nsigma-2,elem) - U_rhs(nsigma-1,elem))
                endif
                b = V_rhs(nsigma-2,elem) - V_rhs(nsigma-1,elem)
                if (b == 0.0_WP) then
                    Sy = 0.0_WP
                else
                    Sy = (V_rhs(nsigma-3,elem) - V_rhs(nsigma-1,elem))/(V_rhs(nsigma-2,elem) - V_rhs(nsigma-1,elem))
                endif

                Z_1 = w_cv(1,elem)*Z(nsigma-1,elem2D_nodes(1,elem)) &
                    + w_cv(2,elem)*Z(nsigma-1,elem2D_nodes(2,elem)) &
                    + w_cv(3,elem)*Z(nsigma-1,elem2D_nodes(3,elem)) &
                    + w_cv(4,elem)*Z(nsigma-1,elem2D_nodes(4,elem))
                Z_1 = dmean - Z_1
                Z_2 = w_cv(1,elem)*Z(nsigma-2,elem2D_nodes(1,elem)) &
                    + w_cv(2,elem)*Z(nsigma-2,elem2D_nodes(2,elem)) &
                    + w_cv(3,elem)*Z(nsigma-2,elem2D_nodes(3,elem)) &
                    + w_cv(4,elem)*Z(nsigma-2,elem2D_nodes(4,elem))
                Z_2 = dmean - Z_2
                Z_3 = w_cv(1,elem)*Z(nsigma-3,elem2D_nodes(1,elem)) &
                    + w_cv(2,elem)*Z(nsigma-3,elem2D_nodes(2,elem)) &
                    + w_cv(3,elem)*Z(nsigma-3,elem2D_nodes(3,elem)) &
                    + w_cv(4,elem)*Z(nsigma-3,elem2D_nodes(4,elem))
                Z_3 = dmean - Z_3

!                ZZ_1(elem) = Z_1
!                Vrhs(elem) = V_rhs(nsigma-1,elem)

                lambda_2 = Z_2/Z_1
                lambda_3 = Z_3/Z_1

                zx_MOST = (DLOG(lambda_3) - Sx*DLOG(lambda_2))/(am*(Sx*lambda_2 - Sx - lambda_3 + 1.0_WP))
                zy_MOST = (DLOG(lambda_3) - Sy*DLOG(lambda_2))/(am*(Sy*lambda_2 - Sy - lambda_3 + 1.0_WP))
                !   zx_MOST = abs(zx_MOST)
                !   zy_MOST = abs(zy_MOST)

                if (zx_MOST >= 0.0_WP) zx_MOST = -0.01_WP
                if (zy_MOST >= 0.0_WP) zy_MOST = -0.01_WP

                !zyMOST(elem) = zy_MOST

                 !  if (zx_MOST >= zx_MOST_max) zx_MOST_max = zx_MOST
                 !  if (zy_MOST >= zy_MOST_max) zy_MOST_max = zy_MOST
                 !  if (zx_MOST <= zx_MOST_min) zx_MOST_min = zx_MOST
                 !  if (zy_MOST <= zy_MOST_min) zy_MOST_min = zy_MOST

                u_st = cka*(U_rhs(nsigma-3,elem) - U_rhs(nsigma-1,elem))/(DLOG(lambda_3) + am*zx_MOST*(lambda_3 - 1.0_WP))
                v_st = cka*(V_rhs(nsigma-3,elem) - V_rhs(nsigma-1,elem))/(DLOG(lambda_3) + am*zy_MOST*(lambda_3 - 1.0_WP))


!               if (u_st == 0.0_WP) u_st = cka*0.01_WP/(DLOG(lambda_3) + am*zx_MOST*(lambda_3 - 1.0_WP))

!               if (v_st == 0.0_WP) v_st = cka*0.01_WP/(DLOG(lambda_3) + am*zx_MOST*(lambda_3 - 1.0_WP))



!               if (u_st >= 0.2_WP)u_st = 0.2_WP
!               if (v_st >= 0.2_WP)v_st = 0.2_WP
!               if (u_st < -0.2_WP)u_st = -0.2_WP
!               if (v_st < -0.2_WP)v_st = -0.2_WP

!               if (u_st > 0.0_WP .and. u_st < 0.0001_WP) u_st = 0.0001_WP
!               if (v_st > 0.0_WP .and. v_st < 0.0001_WP) v_st = 0.0001_WP
!               if (u_st < -0.0_WP .and. u_st > -0.0001_WP) u_st = -0.0001_WP
!               if (v_st < 0.0_WP .and. v_st > -0.0001_WP) v_st = -0.0001_WP


                !u_star(elem) = u_st
                !v_star(elem) = v_st


                 !  if (u_st >= zx_MOST_max) zx_MOST_max = u_st
                 !  if (v_st >= zy_MOST_max) zy_MOST_max = v_st
                 !  if (u_st <= zx_MOST_min) zx_MOST_min = u_st
                 !  if (v_st <= zy_MOST_min) zy_MOST_min = v_st

                 !   c = Sx*am*(lambda_2 - lambda_3)
                 !    if (c == 0.0_WP) then
                 !  Lu_st = z0b_min
                 !    else
                 !  bu_st = (DLOG(lambda_3) - Sx*DLOG(lambda_2))/c
                 !  bu_st = (bu_st*u_st**2)/(Z_1*cka)
                 !  Lu_st = u_st**2/(cka*bu_st)
                 ! Lu_st = abs(Sx*am*(lambda_2 - lambda_3)/(DLOG(lambda_3) - Sx*DLOG(lambda_2)))
                 !    endif
                 !  d = Sy*am*(lambda_2 - lambda_3)
                 !    if (d == 0.0_WP) then
                !   Lv_st = z0b_min
                 !    else
                 !  bv_st = (DLOG(lambda_3) - Sy*DLOG(lambda_2))/d
                 !  bv_st = (bv_st*v_st**2)/(Z_1*cka)
                 !  Lv_st = v_st**2/(cka*bv_st)
                 !    endif
                 ! Lv_st = abs(Sy*am*(lambda_2 - lambda_3)/(DLOG(lambda_3) - Sy*DLOG(lambda_2)))

                c = U_rhs(nsigma-1,elem)*cka/u_st - am*zx_MOST
                d = V_rhs(nsigma-1,elem)*cka/v_st - am*zy_MOST
                z0b_new_x = Z_1*exp(-c)
                z0b_new_y = Z_1*exp(-d)

            !    z0b_new_x = 0.01_WP
            !    z0b_new_y = 0.01_WP

                if (z0b_new_x >= Z_1)z0b_new_x = Z_1
                if (z0b_new_y >= Z_1)z0b_new_y = Z_1
                if (z0b_new_x < 0.001_WP)z0b_new_x = 0.001_WP
                if (z0b_new_y < 0.001_WP)z0b_new_y = 0.001_WP

                !bottomfric_x(elem) = z0b_new_x
                !bottomfric_y(elem) = z0b_new_y


                if (z0b_new_x >= zx_MOST_max) zx_MOST_max = z0b_new_x
                if (z0b_new_y >= zy_MOST_max) zy_MOST_max = z0b_new_y
                if (z0b_new_x <= zx_MOST_min) zx_MOST_min = z0b_new_x
                if (z0b_new_y <= zy_MOST_min) zy_MOST_min = z0b_new_y


  ! if (Lu_st <= z0b_min) Lu_st = z0b_min
  ! if (Lv_st <= z0b_min) Lv_st = z0b_min

                U_tmp(nsigma-1) = DLOG(Z_1/z0b_new_x) + am*zx_MOST
                U_tmp(nsigma-1) = U_tmp(nsigma-1)*u_st/cka
                V_tmp(nsigma-1) = DLOG(Z_1/z0b_new_y) + am*zy_MOST
                V_tmp(nsigma-1) = V_tmp(nsigma-1)*v_st/cka

                xx1(2) = 0.0_WP
                yy1_u(2) =  U_rhs(1,elem)
                yy1_v(2) =  V_rhs(1,elem)
                do nz=2,nsigma-2
                    Fx = U_rhs(nz,elem)
                    Fy = V_rhs(nz,elem)

                  !    if (nz == 1) then
                  !    delta_1 = 0.5_WP*Je(nz,elem)
                  !    else
                    delta_1 = 0.5_WP*(Je(nz,elem) + Je(nz-1,elem))
                   !   endif
                    delta_2 = 0.5_WP*(Je(nz,elem) + Je(nz+1,elem))
                    zt = dt/Je(nz,elem)
                    a = 0.5_WP*zt*Av(nz,elem)/delta_1
                    c = 0.5_WP*zt*Av(nz+1,elem)/delta_2
                    b = a + c + 1.0_WP
                    d = b - a*xx1(nz)
                    xx1(nz+1) = c/d
                    yy1_u(nz+1) = (a*yy1_u(nz) + Fx)/d
                    yy1_v(nz+1) = (a*yy1_v(nz) + Fy)/d

                enddo
!+++++++++++++++++++++++++++++++
! bottom boundary conditions
!+++++++++++++++++++++++++++++++

                do nz=nsigma-2,1,-1
                 U_tmp(nz) = U_tmp(nz+1)*xx1(nz+1) + yy1_u(nz+1)
                 V_tmp(nz) = V_tmp(nz+1)*xx1(nz+1) + yy1_v(nz+1)
                enddo

                do nz=1,nsigma-1
                 U_rhs(nz,elem) = U_tmp(nz)*mask_wd(elem)
                 V_rhs(nz,elem) = V_tmp(nz)*mask_wd(elem)
                enddo

            else

                do nz=1,nsigma-1
                 U_rhs(nz,elem) = 0.0_WP
                 V_rhs(nz,elem) = 0.0_WP
                enddo

            endif

        endif

    enddo

 ! write(99,*) zx_MOST_max, zy_MOST_max, zx_MOST_min, zy_MOST_min
#ifdef USE_MPI
    call exchange_elem(U_rhs)
    call exchange_elem(V_rhs)
#endif
    deallocate(xx1, yy1_u, yy1_v, U_tmp, V_tmp)

end subroutine vert_impl_MOST
