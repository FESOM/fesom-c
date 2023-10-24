subroutine oce_mixing_pp
    !  Compute Richardson number dependent Av and Kv following
    !  Pacanowski and Philander, 1981
    !  Av = Avmax * factor**2 + Av0,
    !  Kv = Kvmax * factor**3 + Kv0,
    !  factor=1/(1+10Ri),  Ri=N**2/(dU/dz)**2 is the Richardson number
    !                     N is the buoyancy frequency
    !                     dU/dz is the vertical velocity shear
    !  Avmax, Kvmax are tunable

    ! Output: Av(2:nlevels(elem)-1,:)  == vert. visc. coeff. at zbar(2:...)
    !         Kv(2:nlevels_nod2D(node),:) == vert. diff. coeff. at zbar(2:...)
!
! modification by Alexey Androsov
! (for sigma coordinates)
! 08.10.14
!
USE o_PARAM
USE o_MESH
USE o_ARRAYS
!a USE g_PARSUP
IMPLICIT NONE

real(kind=WP)                :: dz_inv, bv, shear, a, rho_up, rho_dn, t, s
integer                           :: node, nz, elem, elnodes(4)

  DO node=1, nod2D
     DO nz=2,nsigma-1
        dz_inv=1.0_WP/(Z(nz,node) - Z(nz-1,node))
        t=TF(nz-1, node)
        s=SF(nz-1, node)
        call densityJM(t, s, -zbar(nz,node), rho_up)
        t=TF(nz, node)
        s=SF(nz, node)
        call densityJM(t, s, -zbar(nz,node), rho_dn)


        ! both densities are referenced to the plane
        ! where flux is computed (zbar(nz)
        bv  = -g*dz_inv*(rho_up-rho_dn)/density_0

        if (bv<=0.0_WP) then
           a=1.0_WP
        else
        shear=(Unode(nz,node)-Unode(nz-1,node))**2 +&
	          (Vnode(nz,node)-Vnode(nz-1,node))**2
        shear=shear*dz_inv*dz_inv
	a=shear/(shear+5.0_WP*bv+1.0d-14)  ! To avoid NaNs at start
        end if
        Kv(nz,node)=a !/10.0_WP
     END DO
   END DO

   DO elem=1,elem2D
     elnodes=elem2D_nodes(:,elem)
     DO nz=2,nsigma-1
  	Av(nz,elem)=mix_coeff_PP*sum(w_cv(1:4,elem)*Kv(nz,elnodes)**2)+A_ver
     END DO
   END DO

   DO node=1,nod2D
     DO nz=2,nsigma-1
        Kv(nz,node)=mix_coeff_PP*Kv(nz,node)**3+K_ver
     END DO
   END DO
end subroutine oce_mixing_pp

!========================================================================================
subroutine d3_end

    ! Output: Av(2:nlevels(elem)-1,:)  == vert. visc. coeff. at zbar(2:...)
    !         Kv(2:nlevels_nod2D(node),:) == vert. diff. coeff. at zbar(2:...)
!
! Alexey Androsov
! 20.01.15
!
USE o_PARAM
USE o_MESH
USE o_ARRAYS
IMPLICIT NONE

real(kind=WP)                :: rho_up, rho_dn, t, s, a, b, c
real(kind=WP)                :: dmean, dmean2, d1, d2, d, delta_b, B1cka2
real(kind=WP)                :: ap1, ap2, c1, zt
real(kind=WP)                :: snu_p, snu_m, fr, am, cm, bm
real(kind=WP)                :: t_inv, bbb_max, bott_R, M_pr, abc
real(kind=WP)                :: aa, dss, Zh, Zd, Z0, scale_t, u_din_sur
real(kind=WP), allocatable  ::  xx1(:), yy1(:), btr(:), btr1(:), Ri(:), b1(:), b2(:), c2(:)

integer                           :: node, nz, elem, elnodes(4), it

 allocate(xx1(nsigma),yy1(nsigma),btr(nsigma), btr1(nsigma),Ri(nsigma), b1(nsigma), b2(nsigma), c2(nsigma))
 xx1 = 0.0_WP
 yy1 = 0.0_WP
 btr  = 0.0_WP
 btr1 = 0.0_WP
 Ri = 0.0_WP

 t_inv = 1.0_WP/dt
 B1cka2 = 16.6_WP*cka**2

do node=1,nod2D

 dmean = max(Dmin,depth(node) + eta_n(node))
 dmean2 = dmean**2

 do nz=1,nsigma
   btr(nz) = bt(nz,node)
 enddo

 do nz=2,nsigma-1
    d1 = sigma(nz-1) - sigma(nz)
    d2 = sigma(nz)   - sigma(nz+1)
    d = 0.5_WP*(d1 + d2)
!+++++++++++++++++++++++++++++++++++++++
! compute vertical shear production
!+++++++++++++++++++++++++++++++++++++++
    ap1 = Unode(nz-1,node) - Unode(nz,node)
    ap2 = Vnode(nz-1,node) - Vnode(nz,node)
    b1(nz)  = snu(nz,node)*(ap1**2 + ap2**2)/d**2
    b1(nz) = b1(nz)/dmean2

!+++++++++++++++++++++++++++++++++++++++
! compute baroclinic production
!+++++++++++++++++++++++++++++++++++++++
    t=TF(nz-1, node)
    s=SF(nz-1, node)
    call densityJM(t, s, -zbar(nz,node), rho_up)
    t=TF(nz, node)
    s=SF(nz, node)
    call densityJM(t, s, -zbar(nz,node), rho_dn)
    ap1 = (rho_up - rho_dn)/d
      if (ap1 > 0.0_WP) then
       b2(nz) = 0.0_WP
      else
       b2(nz) = -0.5_WP*g*Kv(nz,node)*ap1/(density_0*dmean)
      endif
!+++++++++++++++++++++++++++++++++++++++
! compute diffusion
!+++++++++++++++++++++++++++++++++++++++
    snu_p = 0.5_WP*(snu(nz-1,node) + snu(nz,node))
    snu_m = 0.5_WP*(snu(nz,node) + snu(nz+1,node))
    c2(nz) = 0.5_WP*0.73_WP*( snu_p*(bt(nz-1,node) - bt(nz,node))/d1 - &
				       snu_m*(bt(nz,node) - bt(nz+1,node))/d2 )/d
    c2(nz) = c2(nz)/dmean2

 enddo
!++++++++++++++++++++++++++++++++++++++++
!  surface boundary conditions (???)
!++++++++++++++++++++++++++++++++++++++++
    xx1(2) = 1.0_WP
    yy1(2) = 0.0_WP
!++++++++++++++++++++++++++++++++++++++++
!  wind boundary conditions
!++++++++++++++++++++++++++++++++++++++++
!!   for dynamical velocity u_*=0.05*W
!!   W ---> wind speed
!!****************************************
   d1 = 0.5_WP*(sigma(nsigma-1) - sigma(nsigma))
   u_din_sur = 0.05_WP*(windx(node)**2 + windy(node)**2)
   snu_m = 0.5_WP*(snu(nsigma,node) + snu(nsigma-1,node))
   yy1(2) = (0.4d-3*u_din_sur**3*dmean*d1)/(snu_m*0.76_WP)

   do it=1,200

  do nz=2,nsigma-1
    d1 = sigma(nz-1) - sigma(nz)
    d2 = sigma(nz)   - sigma(nz+1)
    d = 0.5_WP*(d1 + d2)

!+++++++++++++++++++++++++++++++++++++++
! compute dissipation
!+++++++++++++++++++++++++++++++++++++++
    c1 = 0.046_WP*btr(nz)**2/snu(nz,node)

!+++++++++++++++++++++++++++++++++++++++
! right hand for budget turbulent eq.
!+++++++++++++++++++++++++++++++++++++++
    fr = - b1(nz) - b2(nz) - c1 - t_inv*bt(nz,node) - c2(nz)

!+++++++++++++++++++++++++++++++++++++++
! 3points Thoma's schema coeff.
!+++++++++++++++++++++++++++++++++++++++
    snu_p = 0.5_WP*(snu(nz-1,node) + snu(nz,node))
    snu_m = 0.5_WP*(snu(nz,node) + snu(nz+1,node))

    am = 0.5_WP*(snu_p*0.73_WP)/(dmean2*d1*d)
    cm = 0.5_WP*(snu_m*0.73_WP)/(dmean2*d2*d)
    bm = am + cm + (2.0_WP*0.046_WP*btr(nz))/snu(nz,node) + t_inv
    a = bm - am*xx1(nz)
    xx1(nz+1) = cm/a
    yy1(nz+1) = (am*yy1(nz) - fr)/a

  enddo
!++++++++++++++++++++++++++++++++++++++++
!  bottom boundary conditions
!++++++++++++++++++++++++++++++++++++++++
   delta_b = Jc(nsigma-1,node) + z0b_min
   bott_R = B1cka2/(dlog(delta_b/z0b_min))**2
   btr1(nsigma) =  bott_R*(Unode(nsigma-1,node)**2 + Vnode(nsigma-1,node)**2)

   do nz = nsigma,2,-1
    btr1(nz-1) = btr1(nz)*xx1(nz) + yy1(nz)
   enddo
!++++++++++++++++++++++++++++++++++++++++
!  convergence of iteration
!  if max diff < 10^(-6)  OK
!++++++++++++++++++++++++++++++++++++++++

      bbb_max = 0.0_WP

	do nz=1,nsigma
	  a = btr(nz) - btr1(nz)
	  b = max(btr1(nz),btr(nz)) + 1.0d-8
	  c = dabs(a/b)
	 if (c > bbb_max) bbb_max = c
	enddo
	if (bbb_max < 1.0d-6) go to 5
!++++++++++++++++++++++++++++++++++++++++
!  next iteration
!++++++++++++++++++++++++++++++++++++++++
       do nz=1,nsigma
        btr(nz) = btr1(nz)
       enddo

       enddo ! end for iteration

5     continue

       do nz=1,nsigma
	bt(nz,node) = max(btr(nz),1.0d-8)

!++++++++++++++++++++++++++++++++++++++++
! scale of turbulence
!++++++++++++++++++++++++++++++++++++++++
       a = max(depth(node), Dmin)
       zt = sigma(nz)*a
       Zh = a + zt + z0b_min
       Zd = abs(eta_n(node) - zt + z0s_min)
       Z0 = 1.0_WP  - beta_scale*zh*zd/dmean**2

!	 if (Ri(nz) > 0.0_WP) then
!	  abc = 1.0_WP - 5.0_WP*Ri(nz)
!	 else
!	  abc = 1.0_WP
!	  endif
         scale_t = cka*Zh*Zd*Z0/dmean
         snu(nz,node) = 1.5_8*scale_t*sqrt(bt(nz,node)) + snul
	 if (snu(nz,node) >= 1.0_WP) snu(nz,node) = 1.0_WP
	 Kv(nz,node) = Pr_num*(snu(nz,node))
        end do

       enddo ! end for nodes

     DO elem=1,elem2D
     elnodes=elem2D_nodes(:,elem)
     DO nz=1,nsigma
  	Av(nz,elem)=sum(w_cv(1:4,elem)*snu(nz,elnodes))
     END DO
   END DO

 deallocate(xx1, yy1, btr, btr1, Ri, b1, b2, c2)

end subroutine d3_end

subroutine GOTM
!==============================================================================|
!  A ROUTINE TO CALL THE GOTM TURBULENCE LIBRARY                               |
!==============================================================================|
!--Variables from our 3D Model------------------------------
USE o_PARAM, only: nsigma, g, z0b_min, z0s_min, snul, WP, charnock, dt, density_0, za, Dmin,comp_sediment
USE o_MESH, only: nod2D, elem2D, depth, sigma, w_cv, elem2D_nodes, C_d_el
USE o_ARRAYS, only: taux_node, tauy_node,taux, tauy, Unode, Vnode, Kv, Av, Av_node, teps, tke, L, eta_n, &
TF, SF, CF, zbar, Jc , z0b_gotm_el

!--Variables from GOTM Modules-------------------------------

USE turbulence, only: do_turbulence
USE turbulence, only: tke1d => tke, eps1d => eps, l1d => L, num, nuh !Shear and buoyancy productions can be included  P1d => P, B1d => B, if it is necessary


!--Local Variables (Temporary)-------------------------------
IMPLICIT NONE
integer                   :: node, elem, nz, J, elnodes(4)
real(kind=WP)             :: rr, ztemp, d1, d2, d, t, s, c,rho_up, rho_dn,  ap1, ap2, txr, tyr
real(kind=WP)             :: den_0(1:nsigma, 1:nod2D)
DOUBLE PRECISION          :: NNsq(1:nsigma),SSsq(1:nsigma)



!--Variables for Interfacing with GOTM-----------------------

   DOUBLE PRECISION          :: u_taus,u_taub,z0s_gotm,z0b_gotm, z0b_gotm_nodes(1:nod2D),dmean
   DOUBLE PRECISION          :: h(0:nsigma-1), C_d_node(1:nod2D)
   DOUBLE PRECISION          :: NN1d(0:nsigma-1),SS1d(0:nsigma-1)

!--Local parameters---------------------------------------------------------------------------------------------

   REAL(kind=WP), PARAMETER       :: KAPPA  = .4          !!Von Karman's Constant
   REAL(kind=WP), PARAMETER       :: CHARNOK_VAL = 1400.  !!Charnok Constant
   INTEGER, PARAMETER             :: MaxItz0b = 10        !!Number of iterations of local bottom roughness length



!  Initialize bottom and surface stress to zero
   u_taub = 0.0_WP
   u_taus = 0.0_WP
   NNsq = 0.0_WP
   SSsq = 0.0_WP

 !  Define den_0, which equals to density_0, but has general dim = [nsigma,nod2D]

  den_0= density_0
  ! there is no need to recalculate taux_node,... it is done in fv_sbc
!   call compute_el2nodes_2D_2fields(taux,tauy,taux_node,tauy_node)
   call compute_el2nodes_2D(C_d_el,C_d_node)
   !VF z0b_gotm_el is calculated in the average dynamics routine now, here only interpolation on nodes takes place.
   call compute_el2nodes_2D(z0b_gotm_el, z0b_gotm_nodes)

  DO node=1,nod2D  ! Main loop over all elements

 !+++++++++++++++++++++++++++++++++++++++
 !Set up Depth, [m]
 !+++++++++++++++++++++++++++++++++++++++
 dmean = max(depth(node) + eta_n(node),Dmin)

  DO nz=2,nsigma-1
   ! d1 = sigma(nz-1) - sigma(nz)
   ! d2 = sigma(nz) - sigma(nz+1)
    d = 0.5_WP*(Jc(nz-1,node) + Jc(nz,node))

 !+++++++++++++++++++++++++++++++++++++++
 !  Calculate Buoyancy Frequency Squared (NNsq)
 !+++++++++++++++++++++++++++++++++++++++

    if (comp_sediment) then
    t=TF(nz-1, node)
    s=SF(nz-1, node)
    c=CF(nz-1, node)
    call densityJM(t, s, -zbar(nz,node), rho_up,c)
    t=TF(nz, node)
    s=SF(nz, node)
    c=CF(nz, node)
    call densityJM(t, s, -zbar(nz,node), rho_dn,c)
    else
    t=TF(nz-1, node)
    s=SF(nz-1, node)
    call densityJM(t, s, -zbar(nz,node), rho_up)
    t=TF(nz, node)
    s=SF(nz, node)
    call densityJM(t, s, -zbar(nz,node), rho_dn)
    endif


    NNsq(nz) = -g*(rho_up - rho_dn)/(d*density_0)

 !+++++++++++++++++++++++++++++++++++++++
 !  Calculate Shear Frequency Squared (SSsq),
 !  This discretisation of vertical shear squared (together with NNsq discretisation) should be checked
 !  in order to guarantee the conservation of kinetic energy when transformed
 !  from mean kinetic energy to turbulent kinetic energy.
 !+++++++++++++++++++++++++++++++++++++++
    ap1 = Unode(nz-1,node) - Unode(nz,node)
    ap2 = Vnode(nz-1,node) - Vnode(nz,node)
	SSsq (nz)= (ap1**2 + ap2**2)/(d)**2
     END DO
 !Set BC's as GOTM Does
   NNsq(1)  = NNsq(2)
   NNsq(nsigma) = NNsq(nsigma-1)
   SSsq(1)  = SSsq(2)
   SSsq(nsigma) = SSsq(nsigma-1)

     !+++++++++++++++++++++++++++++++++++++++
	!Surface and bottom friction Velocities [m/s], GOTM formulation for the bottom friction
	!+++++++++++++++++++++++++++++++++++++++


   !+++++++++++++++++++++++++++++++++++++++
   ! Set Surface Roughness Height Using Charnok Formula[m], if charnock = true

  ! taux and tauy are the surface shear-stresses


   u_taus=sqrt(sqrt(taux_node(node)**2+tauy_node(node)**2)/density_0)


   if (charnock) then
    z0s_gotm = charnok_val*u_taus**2/g

    if(z0s_gotm < z0s_min) z0s_gotm = z0s_min

   else
      z0s_gotm=z0s_min
   end if

   !+++++++++++++++++++++++++++++++++++++++
   !  Iterate bottom roughness length MaxItz0b times

  rr=dsqrt(C_d_node(node))

  u_taub = rr*dsqrt( Unode(nsigma-1,node)**2  + Vnode(nsigma-1,node)**2)
  z0b_gotm=z0b_gotm_nodes(node)

      !Set up Layer Thicknesses [m]

      h(1:nsigma-1) = Jc(nsigma-1:1:-1,node)

      !Set Up 1-D Arrays
      SS1d(0:nsigma-1)  = SSsq(nsigma:1:-1)    !Shear Frequency Squared       [1/s^2]
      NN1d(0:nsigma-1)  = NNsq(nsigma:1:-1)    !Buoyancy Frequency Squared    [1/s^2]

      num(0:nsigma-1)   = Av_node(nsigma:1:-1,node)   !Vertical Kinematic Viscosity  [m^2/s]
      nuh(0:nsigma-1)   = Kv(nsigma:1:-1,node)     !Vertical Kinematic Viscosity  [m^2/s]    [1/s^2]
      tke1d(0:nsigma-1) = tke(nsigma:1:-1,node)  !Turbulent Kinetic Energy      [m^2/s^2]
      eps1d(0:nsigma-1) = teps(nsigma:1:-1,node) !Turbulence Dissipation Rate   [m^2/s^3]
      l1d(0:nsigma-1)   = L(nsigma:1:-1,node)


      !Update Turbulence Model

       call  do_turbulence(nsigma-1,dt,sum(Jc(:,node)),u_taus,u_taub,z0s_gotm,z0b_gotm,h, NN1d,SS1d)

      !Update 3D Fields of tke,eps,Av(KM),Kv(KH) and tlen, you can also include shear and buoyancy productions P,B, if it is necessary

       tke(1:nsigma, node) = tke1d(nsigma-1:0:-1)   ! turbulent kinetic energy
       teps(1:nsigma, node) = eps1d(nsigma-1:0:-1)  ! dissipation rate
	   if (teps(1,node)/= teps(1,node)) then
	   teps(1,node)=0.0
	   endif
       Av_node(1:nsigma, node) = min(num(nsigma-1:0:-1),2.0_WP)+ snul !eddy viscosity
       Kv(1:nsigma, node) =  min(nuh(nsigma-1:0:-1),2.0_WP) + snul    !eddy diffusivity
       L(1:nsigma, node) = l1d(nsigma-1:0:-1)         !integral length scale, calculated based on eps

   END DO


   do elem=1,elem2D
    elnodes=elem2D_nodes(:,elem)
     do nz=1,nsigma
      Av(nz,elem) = sum(w_cv(1:4,elem)*Av_node(nz,elnodes))
     enddo
   enddo

   RETURN
  END subroutine GOTM
