!=====================================================================================
SUBROUTINE vel_gradients

  ! Compute derivatives of 3D velocity by least square interpolation.
  ! The interpolation coefficients are already saved
  ! For the no-slip case, it is assumed that velocity at 
  ! the boundary edge == 0. For the free-slip case, there are only 2
  ! neighbours

  !a modification by Alexey Androsov
  !a (for sigma coordinates, and for tri-quad elements)
  !a 09.10.14

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  USE g_PARSUP
  use g_comm_auto

  IMPLICIT NONE

  real(kind=WP)    :: u, v, un
  real(kind=WP)    :: zc(2)
  integer         :: el, el2, j, nz,  jend

!$OMP DO 

  DO el=1, myDim_elem2D

     jend = 4
     if(elem2D_nodes(1,el)==elem2D_nodes(4,el)) jend = 3

      
     vel_grad_ux(1:nsigma-1,el) = 0._WP
     vel_grad_uy(1:nsigma-1,el) = 0._WP
     vel_grad_vx(1:nsigma-1,el) = 0._WP
     vel_grad_vy(1:nsigma-1,el) = 0._WP

     DO j=1,jend
        el2=elem_neighbors(j,el)
      
        if (el2>0) then
           ! ======================
           ! Filling in velocity gradients:
           ! ======================
           DO nz=1, nsigma-1
              vel_grad_ux(nz,el) = vel_grad_ux(nz,el) +gradient_vec(j,  el)*(U_n(nz,el2)-U_n(nz,el))
              vel_grad_uy(nz,el) = vel_grad_uy(nz,el) +gradient_vec(j+4,el)*(U_n(nz,el2)-U_n(nz,el))
              vel_grad_vx(nz,el) = vel_grad_vx(nz,el) +gradient_vec(j,  el)*(V_n(nz,el2)-V_n(nz,el))
              vel_grad_vy(nz,el) = vel_grad_vy(nz,el) +gradient_vec(j+4,el)*(V_n(nz,el2)-V_n(nz,el))
           END DO

        else

           ! ===============
           ! Boundary element
           ! ===============
           !    (Here we do not have place for virtual velocities
           !     in the velocity array so we use auxiliary array)
           ! ======================
           if(free_slip) then
              zc=edge_dxdy(:,elem_edges(j,el))
              zc(1)=zc(1)*elem_cos(el)
              DO nz=1,nsigma-1
                 un=-2*(U_n(nz,el)*zc(2)-V_n(nz,el)*zc(1))/(zc(1)*zc(1)+zc(2)*zc(2))
                 vel_grad_ux(nz,el) = vel_grad_ux(nz,el) + gradient_vec(j,  el)*un*zc(2)
                 vel_grad_uy(nz,el) = vel_grad_uy(nz,el) + gradient_vec(j+4,el)*un*zc(2)
                 vel_grad_vx(nz,el) = vel_grad_vx(nz,el) - gradient_vec(j,  el)*un*zc(1)
                 vel_grad_vy(nz,el) = vel_grad_vy(nz,el) - gradient_vec(j+4,el)*un*zc(1)
              END DO
           else     ! noslip
              DO nz=1,nsigma-1
                 vel_grad_ux(nz,el) = vel_grad_ux(nz,el) -2._WP*gradient_vec(j,  el)*U_n(nz,el)
                 vel_grad_uy(nz,el) = vel_grad_uy(nz,el) -2._WP*gradient_vec(j+4,el)*U_n(nz,el)
                 vel_grad_vx(nz,el) = vel_grad_vx(nz,el) -2._WP*gradient_vec(j,  el)*V_n(nz,el)
                 vel_grad_vy(nz,el) = vel_grad_vy(nz,el) -2._WP*gradient_vec(j+4,el)*V_n(nz,el)
              END DO
           end if

        end if
     end do   ! cycle over neighbor elements
  END DO
!$OMP END DO

#ifdef USE_MPI
  call exchange_elem(vel_grad_ux)
  call exchange_elem(vel_grad_uy)
  call exchange_elem(vel_grad_vx)
  call exchange_elem(vel_grad_vy)
#endif

END SUBROUTINE vel_gradients

!===================================================================================

SUBROUTINE vel_gradients_2D(Uin)

  ! Compute derivatives of velocity by least square interpolation.
  ! The interpolation coefficients are already saved
  ! For the no-slip case, it is assumed that velocity at 
  ! the boundary edge == 0 (opposite on ghost element). For the free-slip case, 
  ! only normal component is reflected on the ghost element

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  use g_parsup
  use g_comm_auto

  IMPLICIT NONE

  real(kind=WP)             :: u, v
  real(kind=WP)             :: zc(2), un
  real(kind=WP), intent(in) :: Uin(2,myDim_elem2D+eDim_elem2D+eXDim_elem2D)
  integer                   :: el, el2, j,  jend, elnodes(4)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(u,v,zc,un,el,el2,j,jend,elnodes)

!!$print *,'grad_vec1',minval(gradient_vec(1,:)),maxval(gradient_vec(1,:))
  DO el=1, myDim_elem2D

     jend = 4
     if(elem2D_nodes(1,el) == elem2D_nodes(4,el) ) jend = 3

     vel_grad(:,el) = 0.0_WP

     DO j=1,jend
        el2=elem_neighbors(j,el)
      
        if (el2 > 0) then
           ! ======================
           ! Filling in velocity gradients:
           ! ======================
           vel_grad(1,el) = vel_grad(1,el) + gradient_vec(j,  el)*(Uin(1,el2)-Uin(1,el))
           vel_grad(2,el) = vel_grad(2,el) + gradient_vec(j+4,el)*(Uin(1,el2)-Uin(1,el))
           vel_grad(3,el) = vel_grad(3,el) + gradient_vec(j,  el)*(Uin(2,el2)-Uin(2,el))
           vel_grad(4,el) = vel_grad(4,el) + gradient_vec(j+4,el)*(Uin(2,el2)-Uin(2,el))
        else
           ! ===============
           ! Boundary element
           ! ===============

           if(free_slip) then


              zc(1) = edge_dxdy(1,elem_edges(j,el))*elem_cos(el)
              zc(2) = edge_dxdy(2,elem_edges(j,el))
              un    =-2.0_WP*( Uin(1,el)*zc(2) -Uin(2,el)*zc(1))/(zc(1)*zc(1) + zc(2)*zc(2))

              vel_grad(1,el) = vel_grad(1,el) + gradient_vec(j,  el)*un*zc(2)
              vel_grad(2,el) = vel_grad(2,el) + gradient_vec(j+4,el)*un*zc(2)
              vel_grad(3,el) = vel_grad(3,el) - gradient_vec(j,  el)*un*zc(1)
              vel_grad(4,el) = vel_grad(4,el) - gradient_vec(j+4,el)*un*zc(1)
             
           else     ! noslip

              vel_grad(1,el) = vel_grad(1,el) -2._WP*gradient_vec(j,  el)*Uin(1,el)
              vel_grad(2,el) = vel_grad(2,el) -2._WP*gradient_vec(j+4,el)*Uin(1,el)
              vel_grad(3,el) = vel_grad(3,el) -2._WP*gradient_vec(j,  el)*Uin(2,el) 
              vel_grad(4,el) = vel_grad(4,el) -2._WP*gradient_vec(j+4,el)*Uin(2,el) 

           end if
	  
        end if
     end do   ! cycle over neighbor elements
    
  END DO
!$OMP END PARALLEL DO

#ifdef USE_MPI
  call exchange_elem(vel_grad)
#endif

END SUBROUTINE vel_gradients_2D

!===========================================================================

SUBROUTINE tracer_gradient_elements (ttf, ttx, tty)

  ! ttf - tracer field
  ! ttx, tty  elemental gradient*Jacobian of tracer 

  USE o_PARAM
  USE o_MESH
  USE o_ARRAYS

  USE g_PARSUP
  use g_comm_auto

  IMPLICIT NONE

  real(kind=WP), intent(in)  :: ttf(nsigma-1,myDim_nod2D+eDim_nod2D)
  real(kind=WP), intent(out) :: ttx(nsigma-1,myDim_elem2D)
  real(kind=WP), intent(out) :: tty(nsigma-1,myDim_elem2D)
  integer           :: el
  integer           :: n, nz

!$OMP DO
  DO el=1,myDim_elem2D
    DO nz=1, nsigma-1
      ttx(nz,el) = sum(gradient_sca(1:4,el)*ttf(nz,elem2D_nodes(1:4,el))) !*Je(nz,elem)
      tty(nz,el) = sum(gradient_sca(5:8,el)*ttf(nz,elem2D_nodes(1:4,el))) !*Je(nz,elem)
    END DO
  END DO
!$OMP END DO

#ifdef USE_MPI
  call exchange_elem(ttx)
  call exchange_elem(tty)
#endif

END SUBROUTINE tracer_gradient_elements

!=======================================================================
SUBROUTINE tracer_gradient_nodes (ttx, tty, ttxnodes, ttynodes)

  ! ttx, tty  elemental gradient of tracer 
  ! ttxnodes, ttynodes nodal gradients obtained by averaging

  USE o_PARAM
  USE o_MESH
  USE o_ARRAYS

  USE g_PARSUP

  IMPLICIT NONE

  real(kind=WP), intent(in) :: ttx(nsigma-1,myDim_elem2D)
  real(kind=WP), intent(in)   :: tty(nsigma-1,myDim_elem2D)
  real(kind=WP), intent(out)   :: ttxnodes(nsigma-1,myDim_nod2D+eDim_nod2D)
  real(kind=WP), intent(out)   :: ttynodes(nsigma-1,myDim_nod2D+eDim_nod2D)
  integer         :: n, nz, el, k, n_elem
  real(kind=WP)   :: tvol, tx(nsigma-1), ty(nsigma-1)

!$OMP DO
  DO n=1, myDim_nod2D+eDim_nod2D

     n_elem = nod_in_elem2D_num(n)
     tvol = 1._WP/sum(elem_area(nod_in_elem2D(1:n_elem,n)))

     tx(1:nsigma-1) = 0.
     ty(1:nsigma-1) = 0.

     DO k = 1, n_elem
        el = nod_in_elem2D(k,n)
        DO nz=1, nsigma-1
           tx(nz) = tx(nz)+ttx(nz,el)*elem_area(el)
	   ty(nz) = ty(nz)+tty(nz,el)*elem_area(el)
        END DO
     END DO

     ttxnodes(1:nsigma-1,n) = tx(1:nsigma-1)*tvol
     ttynodes(1:nsigma-1,n) = ty(1:nsigma-1)*tvol

   END DO 
!$OMP END DO 

END SUBROUTINE tracer_gradient_nodes

!==============================================================================
SUBROUTINE fill_up_dn_grad (ttx, tty)
!
!a Alexey Androsov
!a 20.10.14
!a
! ttx, tty  elemental gradient of tracer 
USE o_PARAM
USE o_MESH
USE o_ARRAYS
!a USE g_PARSUP
IMPLICIT NONE

real(kind=WP)   :: ttx(nsigma-1,elem2D) 
real(kind=WP)   :: tty(nsigma-1,elem2D) 
integer        :: n, nz, elem, k, k1, edge, ednodes(2)
real(kind=WP)   :: tvol, tx, ty
  
 
  DO edge=1,edge2D
     ednodes=edge_nodes(:,edge)
    if((edge_up_dn_tri(1,edge).ne.0).and.(edge_up_dn_tri(2,edge).ne.0)) then
    
     DO nz=1, nsigma-1
        edge_up_dn_grad(1:2,nz,edge)=ttx(nz,edge_up_dn_tri(:,edge))
	edge_up_dn_grad(3:4,nz,edge)=tty(nz,edge_up_dn_tri(:,edge))
     END DO
  
     else
    
      ! Only linear reconstruction part
     DO nz=1,nsigma-1
        tvol=0.0_WP
	tx=0.0_WP
	ty=0.0_WP
	DO k=1, nod_in_elem2D_num(ednodes(1))
           elem=nod_in_elem2D(k,ednodes(1))
           tvol=tvol+elem_area(elem)
           tx=tx+ttx(nz,elem)*elem_area(elem)
	   ty=ty+tty(nz,elem)*elem_area(elem)
        END DO
	edge_up_dn_grad(1,nz,edge)=tx/tvol
	edge_up_dn_grad(3,nz,edge)=ty/tvol
     END DO
     DO nz=1,nsigma-1
        tvol=0.0_WP
	tx=0.0_WP
	ty=0.0_WP
	DO k=1, nod_in_elem2D_num(ednodes(2))
           elem=nod_in_elem2D(k,ednodes(2))
           tvol=tvol+elem_area(elem)
           tx=tx+ttx(nz,elem)*elem_area(elem)
	   ty=ty+tty(nz,elem)*elem_area(elem)
        END DO
	edge_up_dn_grad(2,nz,edge)=tx/tvol
	edge_up_dn_grad(4,nz,edge)=ty/tvol
     END DO
    
    end if  
   END DO 
   
END SUBROUTINE fill_up_dn_grad

!===========================================================================

SUBROUTINE vel_lapl_gradients(U_c,V_c)
  ! Similar to vel_gradients, the difference is  
  ! the argument and the implementation of free slip option
  
  !a changes have been made by Alexey Androsov
  !a (for sigma coordinates, and for tri-quad elements)
  !a 09.10.14

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  USE g_PARSUP
  use g_comm_auto

  IMPLICIT NONE

  real(kind=WP)    :: u, v
  real(kind=WP)    :: grad_aux(4, nsigma)
  real(kind=WP)    :: U_c(nsigma-1,myDim_elem2D+eDim_elem2D+eXDim_elem2D)
  real(kind=WP)    :: V_c(nsigma-1,myDim_elem2D+eDim_elem2D+eXDim_elem2D)
  integer               :: elem, el, j, nz, elnodes(4), jend
  
  DO elem=1,myDim_elem2D
     elnodes=elem2D_nodes(:,elem)
     jend = 4
     if(elnodes(1)==elnodes(4)) jend = 3
     grad_aux=0.0_WP
     DO j=1,jend
        el=elem_neighbors(j,elem)
      
        if (el>0) then
           ! ======================
           ! Filling in velocity gradients:
           ! ======================
           DO nz=1, nsigma-1
              u=U_c(nz,el)-U_c(nz,elem)
              v=V_c(nz,el)-V_c(nz,elem)
              grad_aux(1,nz)=grad_aux(1,nz)+gradient_vec(j,elem)*u
              grad_aux(2,nz)=grad_aux(2,nz)+gradient_vec(j+4,elem)*u
              grad_aux(3,nz)=grad_aux(3,nz)+gradient_vec(j,elem)*v
              grad_aux(4,nz)=grad_aux(4,nz)+gradient_vec(j+4,elem)*v
           END DO
        else
           ! ===============
           ! Boundary element
           ! ===============
           !    (Here we do not have place for virtual velocities
           !     in the velocity array so we use auxiliary array)
           ! ======================
           ! Filling in velocity gradients:
           ! ======================
           DO nz=1, nsigma-1
              u=-2.0_WP*U_c(nz,elem)
              v=-2.0_WP*V_c(nz,elem)
              grad_aux(1,nz)=grad_aux(1,nz)+gradient_vec(j,elem)*u
              grad_aux(2,nz)=grad_aux(2,nz)+gradient_vec(j+4,elem)*u
              grad_aux(3,nz)=grad_aux(3,nz)+gradient_vec(j,elem)*v
              grad_aux(4,nz)=grad_aux(4,nz)+gradient_vec(j+4,elem)*v
           END DO
        end if
     end do   ! cycle over neighbor elements
     vel_grad_ux(1:nsigma-1,elem)=grad_aux(1,1:nsigma-1)
     vel_grad_uy(1:nsigma-1,elem)=grad_aux(2,1:nsigma-1)
     vel_grad_vx(1:nsigma-1,elem)=grad_aux(3,1:nsigma-1)
     vel_grad_vy(1:nsigma-1,elem)=grad_aux(4,1:nsigma-1)	 
  END DO

#ifdef USE_MPI
  call exchange_elem(vel_grad_ux)
  call exchange_elem(vel_grad_uy)
  call exchange_elem(vel_grad_vx)
  call exchange_elem(vel_grad_vy)
#endif

END SUBROUTINE vel_lapl_gradients
!=================================================================================
SUBROUTINE vel_gradients_puls
! Compute derivatives of 3D velocity by least square interpolation.
! The interpolation coefficients are already saved
! For the no-slip case, it is assumed that velocity at 
! the boundary edge == 0. For the free-slip case, there are only 2
! neighbours
!
!a modification by Alexey Androsov
!a (for sigma coordinates, and for tri-quad elements)
!a 09.10.14
!
USE o_MESH
USE o_ARRAYS
USE o_PARAM

  use g_parsup
  use g_comm_auto

IMPLICIT NONE
real(kind=WP)    :: u, v, un
real(kind=WP)    :: zc(2), grad_aux(4,nsigma-1)
real(kind=WP)    :: u_aux(nsigma), v_aux(nsigma)
integer         :: elem, el, j, nz, elnodes(4), jend

  DO elem=1, myDim_elem2D
   elnodes=elem2D_nodes(:,elem)
    jend = 4
    if(elnodes(1)==elnodes(4)) jend = 3
      grad_aux=0.0_WP
      DO j=1,jend
	  el=elem_neighbors(j,elem)
      
	  if (el>0) then
		! ======================
		! Filling in velocity gradients:
		! ======================
	     DO nz=1, nsigma-1
		u=U_puls(nz,el)-U_puls(nz,elem)
                v=V_puls(nz,el)-V_puls(nz,elem)
         	grad_aux(1,nz)=grad_aux(1,nz)+gradient_vec(j,elem)*u
		grad_aux(2,nz)=grad_aux(2,nz)+gradient_vec(j+4,elem)*u
		grad_aux(3,nz)=grad_aux(3,nz)+gradient_vec(j,elem)*v
		grad_aux(4,nz)=grad_aux(4,nz)+gradient_vec(j+4,elem)*v
	     END DO

	  else
	  ! ===============
	  ! Boundary element
	  ! ===============
	  !    (Here we do not have place for virtual velocities
	  !     in the velocity array so we use auxiliary array)
	  ! ======================
	     if(free_slip) then
	       zc=edge_dxdy(:,elem_edges(j,elem))
	       zc(1)=zc(1)*elem_cos(elem)
	       DO nz=1,nsigma-1
		  un=-2*(U_puls(nz,elem)*zc(2)-V_puls(nz,elem)*zc(1))/sum(zc*zc)
                  u_aux(nz)=U_puls(nz,elem)+un*zc(2)
                  v_aux(nz)=V_puls(nz,elem)-un*zc(1)
	       END DO
             else     ! noslip
	       DO nz=1,nsigma-1
		  u_aux(nz)=-U_puls(nz,elem)
		  v_aux(nz)=-V_puls(nz,elem)     
	       END DO
	     end if 
		 ! ======================
		 ! Filling in velocity gradients:
		 ! ======================
	     DO nz=1, nsigma-1
		u=u_aux(nz)-U_puls(nz,elem)
                v=v_aux(nz)-V_puls(nz,elem)
         	grad_aux(1,nz)=grad_aux(1,nz)+gradient_vec(j,elem)*u
		grad_aux(2,nz)=grad_aux(2,nz)+gradient_vec(j+4,elem)*u
		grad_aux(3,nz)=grad_aux(3,nz)+gradient_vec(j,elem)*v
		grad_aux(4,nz)=grad_aux(4,nz)+gradient_vec(j+4,elem)*v
	     END DO
	   end if
	  end do   ! cycle over neighbor elements
         vel_grad_ux(1:nsigma-1,elem)=grad_aux(1,1:nsigma-1)
	 vel_grad_uy(1:nsigma-1,elem)=grad_aux(2,1:nsigma-1)
	 vel_grad_vx(1:nsigma-1,elem)=grad_aux(3,1:nsigma-1)
	 vel_grad_vy(1:nsigma-1,elem)=grad_aux(4,1:nsigma-1)
  END DO	

#ifdef USE_MPI
  call exchange_elem(vel_grad_ux)
  call exchange_elem(vel_grad_uy)
  call exchange_elem(vel_grad_vx)
  call exchange_elem(vel_grad_vy)
#endif

END SUBROUTINE vel_gradients_puls
!===================================================================================
