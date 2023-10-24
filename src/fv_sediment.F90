 Subroutine sediment
!++++++++++++++++++++++++++
! sediment transport model: Ver. 01
! 05.04.2016
! Androsov Alexey
! alexey.androsov@awi.de
!++++++++++++++++++++++++++
  use o_MESH
  use o_ARRAYS
  use o_PARAM
  implicit none
  
  integer             :: elnodes(4), elem, n, k

  real(kind=WP) :: Zw, a, b, c, uv, a_min, c_e, s
  real(kind=WP) :: Q_cr, u_st_cr, u_st, Chez, dmean
  real(kind=WP) :: tvol, tx, ty, hama_c
  real(kind=WP) :: Tsp_el, D_st_el, sk,w_s1
  real(Kind=WP), allocatable   :: U2D_node(:), V2D_node(:)
  real(Kind=WP), allocatable   :: D_st(:), Tsp(:)
  real(kind=WP), Allocatable   :: qsu(:), qsv(:)
  
  real(kind=WP), allocatable   :: test1(:), test2(:), test3(:), test4(:), test5(:)

  Allocate(U2D_node(nod2D),V2D_node(nod2D))
  Allocate(D_st(nod2D),Tsp(nod2D))
  allocate(qsu(nod2D),qsv(nod2D))
  allocate(test1(nod2D),test2(nod2D),test3(nod2D),test4(nod2D),test5(nod2D))
  
!         Kh_c - horizontal diffusion for concentration [1/s]
!         w_s1    --->  settling velocity of a single particle in tranquil water [m/s]
!         D_sed --->  deposition of non-cohesive sediment (k1*l)
!         C_a    --->  near-bed volumetric sediment concentration
!         f_DW ---> Darcy-Weisbach friction factor (const in this simulation = 0.1)
!         u_st   ---> friction velocity
!         u_sur  --->  surface velocity
!         Ren     ---> particle Reynolds number
!         Q_c   ----> critical Shields parameter
!         Q_p   ----> Shields parameter
!
!-------------------------------------------------------------------------------------
!   compute particle parameter D_st
!   Leo and van Rijn "Sediment Transport, Part II:
!   Suspended load Transport". p 1613.
!

    b = 1.0_WP/3.0_WP
    do n=1,nod2D
         s = plop_s/density_0 - 1.0_WP
         D_st(n) = d_m*(g*s/snu_kin**2)**b
    enddo
!    write(*,*) 'particle parameter max, min: D_st = ' , maxval(D_st),minval(D_st)
!-------------------------------------------------------------------------------------
!   compute 2D velocity in the node points
!
    DO n=1, nod2D 
       tvol=0.0_WP
       tx=0.0_WP
       ty=0.0_WP
       DO k=1, nod_in_elem2D_num(n)
           elem=nod_in_elem2D(k,n)
           tvol=tvol+elem_area(elem)
           tx=tx+U_n_2D(1,elem)*elem_area(elem)
	   ty=ty+U_n_2D(2,elem)*elem_area(elem)
        END DO
	U2D_node(n)=tx/tvol
	V2D_node(n)=ty/tvol
     END DO

    do n=1,nod2D

!******************************************************************
!   compute critical velocity (Shields criterium) u_st_cr=f(D_st,s,d_m)
!******************************************************************
         if (D_st(n) <= 4.0_WP) Q_cr = 0.24_WP/D_st(n)
         If (D_st(n) > 4.0_WP .and. D_st(n) <= 10.0_WP) Q_cr = 0.14_WP/D_st(n)**0.64_WP
         If (D_st(n) > 10.0_WP .and. D_st(n) <= 20.0_WP) Q_cr = 0.04_WP/D_st(n)**0.1_WP
         If (D_st(n) > 20.0_WP .and. D_st(n) <= 150.0_WP) Q_cr = 0.013_WP*D_st(n)**0.29_WP
         If (D_st(n) > 150.0_WP) Q_cr = 0.055_WP
!write(*,*), 'Q_cr', Q_cr
          s = plop_s/density_0 - 1.0_WP
	  b = d_m*g*s
          u_st_cr = sqrt(b*Q_cr)
	  test3(n) = u_st_cr

!************************************
!   compute dynamical velocity  (u_st)
!************************************
          a = sqrt(U2D_node(n)**2 + V2D_node(n)**2)
         !VF z0r is replaced by z0b_min+za, where
         !za is roughness caused by suspended sediment
         !z0b_min is bottom roughness height, minimum value

          sk = (z0b_min+za)*30.0_WP
          dmean=eta_n(n) + depth(n)
          dmean=max(Dmin,dmean)
	  c=max(50.0,12.0_WP*dmean/(sk))
          Chez = 18.0_WP*dlog10(c)
     !     write(*,*), 'Chez', Chez
          u_st = sqrt(g)*a/(Chez)
	  test2(n) = u_st
	  
!************************************************************************
! average settling velocity of the material available for transport (w_s1)
! Simpson and Castelltort 2006, Coupled model of surface water flow, 
! sediment transport and morphological evolution. Computers & Geosciences,
! 32, pp. 1600-1614.
! ************************************************************************
          c = 13.95_WP*snu_kin/d_m
          w_s1 = sqrt(c*c + 1.09_WP*b ) - c
	!  write(*,*), 'w_s1', w_s1
!w_s1=0.001
!************************
! compute transport stage
!************************
          a = u_st**2/(u_st_cr)**2 - 1.0_WP
	  Tsp(n) = max(0.0_WP,a)
        !  write(*,*), 'Tsp', Tsp
!******************************
! compute weight parameter Zw
!******************************
          Zw = w_s1/(cka*(u_st + 2.0_WP*w_s1))
	  Zw = min(0.9999_WP,Zw)

          hama_v(n) = 0.980_WP - 0.198_WP*Zw + 0.032_WP*Zw**2
          hama_c = 0.434_WP + 5.975_WP*Zw + 2.888_WP*Zw**2
	  test1(n) = hama_c
	  
!*********************************
! compute weighted sediment transport qs
!*********************************
          dmean=max(5.0_WP,eta_n(n) + depth(n))
          b = hama_v(n)*dmean*con(n)
          qsu(n) = b*U2D_node(n)
          qsv(n) = b*V2D_node(n)
	  
!*********************************
! akkumulation sediment particles
!*********************************

          a_min=0.01_WP*dmean*10.0_WP  ! for test *10
          c_e = 0.015_WP*d_m*Tsp(n)**1.5_WP/(D_st(n)**0.3_WP*a_min*dmean)
   !       c_e = 0.015_WP*d_m*Tsp(n)**1.5_WP
          test4(n) = c_e
          E_sed(n) = hama_c*w_s1*(c_e - con(n))
    !      if (E_sed(n) > 0.0_WP) E_sed(n)=min(0.04_WP,E_sed(n))
    !      if (E_sed(n) < 0.0_WP) E_sed(n)=max(-0.04_WP,E_sed(n))

   enddo
 !    write(*,*) 'akkumulation of sediment part max, min: E_sed = ' , maxval(E_sed),minval(E_sed)
 !    write(*,*) 'sediment transport max: qsu,qsv = ' , maxval(qsu),minval(qsv)
 !    write(*,*) 'hama =                                         ', maxval(hama_v),maxval(test1)
 !   write(*,*) 'c_e, max, min = ', maxval(test4),minval(test4)
 !   write(*,*) 'transport stage, max, min =', maxval(Tsp),minval(Tsp)
 !    write(*,*) 'dynamical velocity, max, min = ', maxval(test2),minval(test2)
 !    write(*,*) 'critical dynamical velocity, max, min = ', maxval(test3),minval(test3)
 !    write(*,*) '.............................'
!*********************************
! compute bed-load transport qb
!*********************************
    do elem=1,elem2D
         elnodes=elem2D_nodes(:,elem)
          Tsp_el=sum(w_cv(1:4,elem)*Tsp(elnodes))
          D_st_el=sum(w_cv(1:4,elem)*D_st(elnodes))
          a = 0.053_WP*dsqrt(g*(plop_s/density_0 - 1.0_WP))*d_m**1.5_WP/D_st_el**0.3_WP
          uv = dsqrt(U_n_2D(1,elem)**2 + U_n_2D(2,elem)**2)
	  if (uv < 1.e-12) then
	   qbu(elem) = 0.0_WP
	   qbv(elem) = 0.0_WP
         else
          qbu(elem) = a*Tsp_el**2.1_WP*U_n_2D(1,elem)/uv
          qbv(elem) = a*Tsp_el**2.1_WP*U_n_2D(2,elem)/uv
         endif
     enddo

!*********************************
! compute equation for concentration
!*********************************
       call concentration_sed
 !  write(*,*) 'concentration of sedimentation max, min: con = ' , maxval(con),minval(con)
   
        call bottom_evolution
 !  write(*,*) 'bottom variation max, min: h_var = ' , maxval(h_var),minval(h_var)
   
     deallocate(U2D_node,V2D_node)
     deallocate(D_st,Tsp)
     deallocate(qsu,qsv)
     deallocate(test1,test2,test3,test4,test5)

END subroutine sediment
!============================================================================
SUBROUTINE concentration_sed
USE o_MESH
USE o_ARRAYS
USE o_PARAM
IMPLICIT NONE

 integer      :: el(2), enodes(2), n, edge, ed, nodes(2)
 
 real(kind=WP) :: c1, c2, deltaX1, deltaY1, deltaX2, deltaY2
 real(kind=WP) :: Cmean, fD, fDold, x1, y1, dmean, un1 !, relax2clim_ac
 
 real(kind=WP), allocatable :: cHrhs(:), Vel_nor(:)

 allocate(cHrhs(nod2D))
 allocate(Vel_nor(nod2D))

! =========================
! momentum equation for
! concentration in conservative form
! =========================
  DO n=1, nod2D
     fDold=max(Dmin,depth(n) + eta_n_1(n))
     con(n) = con(n)*fDold
  END DO
! =================
! Clean the rhs
! =================          
 cHrhs=0.0_WP
! =================
! Horizontal advection
! =================
  DO edge=1, edge2D 
   enodes=edge_nodes(:,edge)   
   el=edge_tri(:,edge)
   c1 = 0.0_WP
   c2 = 0.0_WP
   deltaX1=edge_cross_dxdy(1,edge)
   deltaY1=edge_cross_dxdy(2,edge)
!========== ============
! First segment
! Linear upwind reconstruction
! ======================
   if(U_n_2D(2,el(1))*deltaX1- U_n_2D(1,el(1))*deltaY1>0.0_WP) then   
      Cmean=con(enodes(2))*hama_v(enodes(2))
   else
      Cmean=con(enodes(1))*hama_v(enodes(1))
   endif
   c1=U_n_2D(2,el(1))*Cmean*deltaX1 - U_n_2D(1,el(1))*Cmean*deltaY1
   cHrhs(enodes(1)) = cHrhs(enodes(1)) + c1
   cHrhs(enodes(2)) = cHrhs(enodes(2)) - c1

! ======================
! Second segment
! Linear upwind reconstruction
!========== ============
   if(el(2)>0) then
   deltaX2=edge_cross_dxdy(3,edge)
   deltaY2=edge_cross_dxdy(4,edge)
   if(U_n_2D(2,el(2))*deltaX2- U_n_2D(1,el(2))*deltaY2<0.0_WP) then   
      Cmean= con(enodes(2))*hama_v(enodes(2))
   else
      Cmean= con(enodes(1))*hama_v(enodes(1))
   endif
   c2=-U_n_2D(2,el(2))*Cmean*deltaX2 + U_n_2D(1,el(2))*Cmean*deltaY2
   cHrhs(enodes(1))=cHrhs(enodes(1)) + c2
   cHrhs(enodes(2))=cHrhs(enodes(2)) - c2

  end if
   
 END DO
  DO n=1,nod2D
     cHrhs(n) = cHrhs(n)/area(n) !+ E_sed(n)
  END DO

!=================
! boundary conditions
! =================
Vel_nor = 0.0_WP

DO ed=1+edge2D_in, edge2D 
   nodes=edge_nodes(:,ed)   
   el=edge_tri(:,ed)
   x1=edge_cross_dxdy(1,ed)
   y1=edge_cross_dxdy(2,ed)
   dmean=sum(eta_n(nodes)+depth(nodes))/2.0_WP
   dmean=max(Dmin,dmean)
   un1=(U_n_2D(2,el(1))*x1- U_n_2D(1,el(1))*y1)*dmean
 
   Vel_nor(nodes(1))=Vel_nor(nodes(1))+un1/area(nodes(1))
   Vel_nor(nodes(2))=Vel_nor(nodes(2))-un1/area(nodes(2))
 END DO

  !   relax2clim_ac = relax2clim/24.0_WP/100.0_WP
  DO n=1, nod2D
     if (index_nod2D(n) == 2) then
       fD=max(Dmin,depth(n) + eta_n_1(n))
       if (Vel_nor(n) > 0.0_WP) then
        cHrhs(n)=cHrhs(n) + relax2clim_ac*(con_bc(n)*fD - con(n))/24.0_WP
       else
       cHrhs(n)=clim_relax*(con_bc(n)*fD - con(n))
       endif
      endif
  END DO
  
! =================
! Update ttfold (to be used on the next time level)
! and compute new ttf
! =================          

  DO n=1,nod2D
     fD=max(Dmin,depth(n) + eta_n(n))
     con(n) = ac(n)*(con(n) + (cHrhs(n)+E_sed(n))*dt_2D)/fD 
  END DO
  
deallocate(cHrhs, Vel_nor)
end subroutine concentration_sed
!=============================================
   subroutine bottom_evolution

  use o_MESH
  use o_ARRAYS
  use o_PARAM
  implicit none

  real(kind=WP) :: c1, c2, a
  integer             :: ed, el(2), enodes(2), n

  real(kind=WP), allocatable :: cHrhs(:)

  allocate(cHrhs(nod2D))
 
  cHrhs = 0.0_WP
 ! ==============
 ! internal edges
 ! ==============
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
    h_var(n)=h_var(n) + ac(n)*a*(cHrhs(n)/area(n) + E_sed(n))
  ENDDO

  deallocate(cHrhs)
  end subroutine bottom_evolution
