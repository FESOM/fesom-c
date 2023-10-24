!!=========================================================================================
!++++++++++++++++++++++++++++++++++++++++++
!  compute density in the vertical column
!  19.09.2014
!  Androsov Alexey
!++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE pressure

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  use g_parsup
  use g_comm_auto

  IMPLICIT NONE

  INTEGER        :: nz, n
  REAL(kind=WP)  :: g_density0_inv, hpress_clim(nsigma,myDim_nod2D+eDim_nod2D), rho_clim
  REAL(kind=WP)  :: mx_rho, mn_rho

  interface

     SUBROUTINE densityJM(t, s, pz, rho_out,con_ss)

       USE o_MESH
       USE o_ARRAYS
       USE o_PARAM

       IMPLICIT NONE

       real(kind=WP), intent(IN)       :: t,s,pz
       real(kind=WP), intent(OUT) :: rho_out
       real(kind=WP)              :: rhopot, bulk, ss3
       real(kind=WP), intent(IN), optional   :: con_ss

     END SUBROUTINE densityJM

  end interface

  g_density0_inv = g/density_0
  !VF, add if based on mixed_BP_grad_OB

  !$OMP PARALLEL PRIVATE(n,nz,rho_clim)
  if (mixed_BP_grad_OB) then

     if (comp_sediment) then
!$OMP DO
        DO n=1,myDim_nod2D+eDim_nod2D

           hpressure(1,n)   = 0.0_WP
           hpress_clim(1,n) = 0.0_WP

           ! VF, compute density in the vertical column and add it to rho_c - current rho, needed in fv_tracer
           ! VF, add effect of the susp. sediments
           DO nz=1,nsigma-1
              call densityJM(TF(nz,n),   SF(nz,n),    -Z(nz,n), rho_c(nz,n),CF(nz,n))
              call densityJM(Tclim(nz,n),Sclim(nz,n), -Z(nz,n), rho_clim, Cclim(nz,n))
              hpressure(nz+1,n)   = hpressure(nz,n)   + rho_c(nz,n)*Jc(nz,n)*g_density0_inv
              hpress_clim(nz+1,n) = hpress_clim(nz,n) + rho_clim   *Jc(nz,n)*g_density0_inv
           END DO
        END DO
!$OMP END DO
     else
!$OMP DO
        DO n=1,myDim_nod2D

           hpressure(1,n)   = 0.0_WP
           hpress_clim(1,n) = 0.0_WP

           ! VF, compute density in the vertical column and add it to rho_c - current rho, needed in fv_tracer
           ! VF, add effect of the susp. sediments
           DO nz=1,nsigma-1
              call densityJM(TF(nz,n),   SF(nz,n),    -Z(nz,n), rho_c(nz,n))
              call densityJM(Tclim(nz,n),Sclim(nz,n), -Z(nz,n), rho_clim)
              hpressure(nz+1,n)   = hpressure(nz,n)   + rho_c(nz,n)*Jc(nz,n)*g_density0_inv
              hpress_clim(nz+1,n) = hpress_clim(nz,n) + rho_clim   *Jc(nz,n)*g_density0_inv
           END DO
        END DO
!$OMP END DO
     endif

#ifdef USE_MPI
     call exchange_nod(hpressure)
     call exchange_nod(hpress_clim)
#endif

     call baroclinic_pressure_gradient_3D_clim

  else
!$OMP DO
     DO n=1,myDim_nod2D+eDim_nod2D

        hpressure(1, n)=0.0_WP

        ! Compute density in the vertical column
        DO nz=1,nsigma-1

           if (comp_sediment) then
              call densityJM(TF(nz,n), SF(nz,n), -Z(nz,n), rho_c(nz,n), CF(nz,n))
           else
              call densityJM(TF(nz,n), SF(nz,n), -Z(nz,n), rho_c(nz,n))
           endif

           hpressure(nz+1,n) = hpressure(nz,n) + rho_c(nz,n)*Jc(nz,n)*g_density0_inv
        END DO
     END DO
!$OMP END DO

#ifdef USE_MPI
     call exchange_nod(hpressure)
     call exchange_nod(hpress_clim)
#endif


     call baroclinic_pressure_gradient_3D

  endif


#ifdef USE_MPI
     call MPI_AllREDUCE(minval(rho_c),mn_rho, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
          MPI_COMM_FESOM_C, MPIerr)

     call MPI_AllREDUCE(maxval(rho_c),mx_rho, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
          MPI_COMM_FESOM_C, MPIerr)

     !if (mype==0) print *,'MIN_RHO',mn_rho
     !if (mype==0) print *,'MAX_RHO',mx_rho
#else
!     print *,'MIN_RHO',minval(rho_c)
!     print *,'MAX_RHO',maxval(rho_c)
#endif

!$OMP END PARALLEL

contains

!===========================================================================
!++++++++++++++++++++++++++++++++++++++++++
!  compute 3D baroclinic pressure gradient
!  02.12.2015
!  Androsov Alexey
!++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE baroclinic_pressure_gradient_3D

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  use g_parsup
  use g_comm_auto

  IMPLICIT NONE

  INTEGER         :: elem, nz, elnodes(4)
  REAL(KIND=WP)   :: pre(4)
  real(kind=WP)   :: Ax, Ay, Bx, By
  real(kind=WP)   :: qd_up, qd_low
  real(kind=WP)   :: acc, acc_to_Bp_bz, accJe_inv

!$OMP DO
  DO elem=1,myDim_elem2D

     elnodes=elem2D_nodes(:,elem)
     Ax = sum(gradient_sca(1:4,elem)*(depth(elnodes)+eta_n(elnodes)))
     Ay = sum(gradient_sca(5:8,elem)*(depth(elnodes)+eta_n(elnodes)))
     Bx = sum(gradient_sca(1:4,elem)*(depth(elnodes)))
     By = sum(gradient_sca(5:8,elem)*(depth(elnodes)))

     acc = sum(w_cv(1:4,elem)*ac(elnodes))
     acc_to_Bp_bz = acc**Bp_bz*mask_bpr(elem)

     qd_up  = 0._WP

     DO nz=1,nsigma-1
        accJe_inv = -acc_to_Bp_bz/Je(nz,elem)

        qd_up = sum(w_cv(1:4,elem)*hpressure(nz,elnodes))
        qd_low = sum(w_cv(1:4,elem)*hpressure(nz+1,elnodes))

        pre = 0.5_WP*(hpressure(nz,elnodes) + hpressure(nz+1,elnodes))*Jc(nz,elnodes)

        Bar_pru_3D(nz,elem) = accJe_inv*(sum(gradient_sca(1:4,elem)*pre) &
             + qd_low*(sigma(nz+1)*Ax - Bx) - qd_up*(sigma(nz)*Ax - Bx) )

        Bar_prv_3D(nz,elem) = accJe_inv*(sum(gradient_sca(5:8,elem)*pre) &
             + qd_low*(sigma(nz+1)*Ay - By) - qd_up*(sigma(nz)*Ay - By) )

        !qd_up = qd_low ! SH this version produced differences in MPI version
                        ! Replaced it by direct calculation

     END DO
     !Bar_pru_3D=Bar_pru_3D*0.0_WP
     !Bar_prv_3D=Bar_prv_3D*0.0_WP

     if (time_jd-time_jd0<30.0_WP) then
        Bar_pru_3D=Bar_pru_3D*(time_jd-time_jd0)/30.0_WP
        Bar_prv_3D=Bar_prv_3D*(time_jd-time_jd0)/30.0_WP
     end if

     ! baroclinic_pressure_gradient_2D (inlined from subroutine)
     Bar_pr_2D(1,elem) = sum(Bar_pru_3D(1:nsigma-1,elem)*Je(1:nsigma-1,elem))
     Bar_pr_2D(2,elem) = sum(Bar_prv_3D(1:nsigma-1,elem)*Je(1:nsigma-1,elem))
  enddo
!$OMP END DO

  !NR inlined into loop above (Cache!!!)
  ! call baroclinic_pressure_gradient_2D

#ifdef USE_MPI
  call exchange_elem(Bar_pr_2D)
#endif

!SHDB print *,'PRTEST hprs: ',mype,minval(hpressure(nsigma-1,:)),maxval(hpressure(nsigma-1,:))

!SHDB print *,'PRTEST Bar_pr1: ',mype,minval(Bar_pr_2D(1,:)),maxval(Bar_pr_2D(1,:))
!SHDB print *,'PRTEST Bar_pr2: ',mype,minval(Bar_pr_2D(2,:)),maxval(Bar_pr_2D(2,:))

END SUBROUTINE baroclinic_pressure_gradient_3D


!==================================================================

SUBROUTINE baroclinic_pressure_gradient_3D_clim

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  use g_parsup
  use g_comm_auto

  IMPLICIT NONE

  INTEGER         :: elem, nz, elnodes(4)
  REAL(KIND=WP)   :: pre(4), pre_clim(4)
  real(kind=WP)   :: Ax, Ay, Bx, By
  real(kind=WP)   :: qd_up, qd_low, qd_up_clim, qd_low_clim
  real(kind=WP)   :: acc, acc_to_Bp_bz, accJe_inv

!$OMP DO
  DO elem=1,myDim_elem2D

     elnodes=elem2D_nodes(:,elem)
     Ax = sum(gradient_sca(1:4,elem)*(depth(elnodes)+eta_n(elnodes)))
     Ay = sum(gradient_sca(5:8,elem)*(depth(elnodes)+eta_n(elnodes)))
     Bx = sum(gradient_sca(1:4,elem)*(depth(elnodes)))
     By = sum(gradient_sca(5:8,elem)*(depth(elnodes)))

     acc = sum(w_cv(1:4,elem)*ac(elnodes))
     acc_to_Bp_bz = acc**Bp_bz*mask_bpr(elem)

     qd_up       = 0._WP
     qd_up_clim  = 0._WP

     DO nz=1,nsigma-1

        accJe_inv = -acc_to_Bp_bz/Je(nz,elem)

        qd_low      = sum(w_cv(1:4,elem)*hpressure(nz+1,elnodes))
        qd_low_clim = sum(w_cv(1:4,elem)*hpress_clim(nz+1,elnodes))

        pre      = 0.5_WP*(hpressure(nz,elnodes)   + hpressure(nz+1,elnodes))  *Jc(nz,elnodes)
        pre_clim = 0.5_WP*(hpress_clim(nz,elnodes) + hpress_clim(nz+1,elnodes))*Jc(nz,elnodes)

        Bar_pru_3D(nz,elem)      = accJe_inv*(sum(gradient_sca(1:4,elem)*pre) &
             + qd_low*(sigma(nz+1)*Ax - Bx)     - qd_up*(sigma(nz)*Ax - Bx) )

        Bar_prv_3D(nz,elem)      = accJe_inv*(sum(gradient_sca(5:8,elem)*pre) &
             + qd_low*(sigma(nz+1)*Ay - By)     - qd_up*(sigma(nz)*Ay - By) )

        Bar_pru_3D_clim(nz,elem) = accJe_inv*(sum(gradient_sca(1:4,elem)*pre_clim) &
             + qd_low_clim*(sigma(nz+1)*Ax - Bx) - qd_up_clim*(sigma(nz)*Ax - Bx) )

        Bar_prv_3D_clim(nz,elem) = accJe_inv*(sum(gradient_sca(5:8,elem)*pre_clim) &
             + qd_low_clim*(sigma(nz+1)*Ay - By) - qd_up_clim*(sigma(nz)*Ay - By) )

        qd_up      = qd_low
        qd_up_clim = qd_low_clim
     END DO

     DO nz=1,nsigma-1
        Bar_pru_3D(nz,elem) = (1.0_WP - acc)*Bar_pru_3D_clim(nz,elem) &
             + acc *Bar_pru_3D(     nz,elem)
        Bar_prv_3D(nz,elem) = (1.0_WP - acc)*Bar_prv_3D_clim(nz,elem) &
             + acc *Bar_prv_3D(     nz,elem)
     enddo

     ! baroclinic_pressure_gradient_2D (inlined from subroutine)
     Bar_pr_2D(1,elem) = sum(Bar_pru_3D(1:nsigma-1,elem)*Je(1:nsigma-1,elem))
     Bar_pr_2D(2,elem) = sum(Bar_prv_3D(1:nsigma-1,elem)*Je(1:nsigma-1,elem))
  enddo
!$OMP END DO

#ifdef USE_MPI
  call exchange_elem(Bar_pr_2D)
#endif

  !NR inlined into loop above (Cache!!!)
  ! call baroclinic_pressure_gradient_2D

END SUBROUTINE baroclinic_pressure_gradient_3D_clim
!===========================================================================
!++++++++++++++++++++++++++++++++++++++++++
!  compute average baroclinic pressure gradient
!  for 2D average on vertical velocity
!  02.12.2015
!  Androsov Alexey
!++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE baroclinic_pressure_gradient_2D

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  use g_parsup

  IMPLICIT NONE

  INTEGER   :: elem

!$OMP DO
  DO elem=1,myDim_elem2D+eDim_elem2D+eXDim_elem2D
     Bar_pr_2D(1,elem) = sum(Bar_pru_3D(1:nsigma-1,elem)*Je(1:nsigma-1,elem))
     Bar_pr_2D(2,elem) = sum(Bar_prv_3D(1:nsigma-1,elem)*Je(1:nsigma-1,elem))
  END DO
!$OMP END DO

END SUBROUTINE baroclinic_pressure_gradient_2D
!===========================================================================


end SUBROUTINE pressure

!===========================================================================

SUBROUTINE densityJM(t, s, pz, rho_out,con_ss)

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  IMPLICIT NONE

  !
  ! - calculates in-situ density as a function of potential temperature
  !   (relative to the surface)
  !   using the Jackett and McDougall equation of state
  !   (Copyright (c) 1992, CSIRO, Australia)
  ! - has been derived from the SPEM subroutine rhocal
  !
  ! Ralph Timmermann, August 2005
  !---------------------------------------------------------------------------


  real(kind=WP), intent(IN)       :: t,s,pz
  real(kind=WP), intent(OUT) :: rho_out
  real(kind=WP)              :: rhopot, bulk, ss3
  real(kind=WP), intent(IN), optional   :: con_ss

  !compute secant bulk modulus

  ! if (s.lt.0.) then
  !    write (*,*)'s<0 happens!',t,s
  !     stop
  !  endif
  !   rho_out=-0.25_WP*(t-10.0_WP)
  !   return

  ss3 = sqrt(s*s*s)

  bulk = 19092.56_WP + t*(209.8925_WP 			&
       - t*(3.041638_WP - t*(-1.852732e-3			&
       - t*(1.361629e-5))))				&
       + s*(104.4077_WP - t*(6.500517_WP			&
       -  t*(.1553190_WP - t*(-2.326469e-4))))		&
       + ss3*(-5.587545_WP				&
       + t*(0.7390729_WP - t*(1.909078e-2)))		&
       - pz *(4.721788e-1 + t*(1.028859e-2		&
       + t*(-2.512549e-4 - t*(5.939910e-7))))		&
       - pz*s*(-1.571896e-2				&
       - t*(2.598241e-4 + t*(-7.267926e-6)))		&
       - pz*ss3					&
       *2.042967e-3 + pz*pz*(1.045941e-5		&
       - t*(5.782165e-10 - t*(1.296821e-7)))		&
       + pz*pz*s					&
       *(-2.595994e-7					&
       + t*(-1.248266e-9 + t*(-3.508914e-9)))

  rhopot = 999.842594_WP				&
       + t*( 6.793952e-2		        &
       + t*(-9.095290e-3			&
       + t*( 1.001685e-4			&
       + t*(-1.120083e-6			&
       + t* 6.536332e-9))))			&
       + s*( 0.824493_WP				&
       + t *(-4.08990e-3			&
       + t *( 7.64380e-5			&
       + t *(-8.24670e-7			&
       + t * 5.38750e-9))))			&
       + ss3*(-5.72466e-3		&
       + t*( 1.02270e-4			&
       + t*(-1.65460e-6)))			&
       + 4.8314e-4*s*s

  rho_out = rhopot*bulk / (bulk + 0.1_WP*pz) - density_0
  if (comp_sediment) then
     rho_out=rho_out+con_ss*(plop_s-(rho_out+density_0))/plop_s
  endif

end subroutine densityJM
!===================================================================================
function theta(s,t,p,pr)
  ! Compute local potential temperature at pr
  ! using bryden 1973 polynomial for adiabatic lapse rate
  ! and runge-kutta 4-th order integration algorithm.
  ! ref: bryden,h.,1973,deep-sea res.,20,401-408
  ! fofonoff,n.,1977,deep-sea res.,24,489-491
  ! units:
  !       pressure        p        decibars
  !       temperature     t        deg celsius (ipts-68)
  !       salinity        s        (ipss-78)
  !       reference prs   pr       decibars
  !       potential tmp.  theta    deg celsius
  ! checkvalue: theta= 36.89073 c,s=40 (ipss-78),t=40 deg c,
  ! p=10000 decibars,pr=0 decibars
  !
  USE o_PARAM

  implicit none
  real(kind=WP)	            :: theta, s, t, p, pr
  real(kind=WP) 		    :: h, xk, q
  real(kind=WP), external	    :: atg

  h = pr - p
  xk = h*atg(s,t,p)
  t = t + 0.5_WP*xk
  q = xk
  p = p + 0.5_WP*h
  xk = h*atg(s,t,p)
  t = t + 0.29289322_WP*(xk-q)
  q = 0.58578644_WP*xk + 0.121320344_WP*q
  xk = h*atg(s,t,p)
  t = t + 1.707106781_WP*(xk-q)
  q = 3.414213562_WP*xk - 4.121320344_WP*q
  p = p + 0.5_WP*h
  xk = h*atg(s,t,p)
  theta = t + (xk-2.0_WP*q)/6.0_WP
  return
end function theta
!=============================================================
function atg(s,t,p)
  ! adiabatic temperature gradient deg c per decibar
  ! ref: bryden,h.,1973,deep-sea res.,20,401-408
  ! units:
  !       pressure        p        decibars
  !       temperature     t        deg celsius (ipts-68)
  !       salinity        s        (ipss-78)
  !       adiabatic       atg      deg. c/decibar
  ! checkvalue: atg=3.255976e-4 c/dbar for s=40 (ipss-78),
  ! t=40 deg c,p0=10000 decibars
   USE o_PARAM

  implicit none
  real(kind=WP)      ::  atg, s, t, p, ds

  ds = s - 35.0_WP
  atg = (((-2.1687e-16*t+1.8676e-14)*t-4.6206e-13)*p   &
       +((2.7759e-12*t-1.1351e-10)*ds+((-5.4481e-14*t        &
       +8.733e-12)*t-6.7795e-10)*t+1.8741e-8))*p             &
       +(-4.2393e-8*t+1.8932e-6)*ds                          &
       +((6.6228e-10*t-6.836e-8)*t+8.5258e-6)*t+3.5803e-5

  return
end function atg
!=============================================================
