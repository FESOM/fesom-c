module d3_module

    implicit none

contains

    subroutine d3_end

        use o_PARAM
        use o_MESH
        use o_ARRAYS
        use g_PARSUP
    
        implicit none

        real(kind=WP)                :: rho_up, rho_dn, t, s, a, b, c
        real(kind=WP)                :: dmean, dmean2, d1, d2, d, delta_b, B1cka2
        real(kind=WP)                :: ap1, ap2, c1, zt
        real(kind=WP)                :: snu_p, snu_m, fr, am, cm, bm
        real(kind=WP)                :: t_inv, bbb_max, bott_R, M_pr, abc
        real(kind=WP)                :: aa, dss, Zh, Zd, Z0, scale_t, u_din_sur
        real(kind=WP), allocatable  ::  xx1(:), yy1(:), btr(:), btr1(:), b1(:), b2(:), c2(:)
        real(kind=WP)                :: deltt
      
        integer                           :: node, nz, elem, elnodes(4), it
      
        allocate(xx1(nsigma),yy1(nsigma),btr(nsigma), btr1(nsigma), b1(nsigma), b2(nsigma), c2(nsigma))
        xx1 = 0.0_WP
        yy1 = 0.0_WP
        btr  = 0.0_WP
        btr1 = 0.0_WP
        Ri = 0.0_WP
      
        t_inv = 1.0_WP/dt
        B1cka2 = 16.6_WP*cka**2
      
        do node=1,myDim_nod2D+eDim_nod2D
      
           !-----------------------------
           ! !!   Changing Pr_num taking surface and bottom T
           ! deltt = (TF(nsigma-1, node) - TF(1, node) -0.5_WP)*5.0_WP
           ! Pr_num = 0.001_WP + 0.0035_WP * (1.0_WP + tanh( deltt ))
           !----------------------------
      
           dmean = max(Dmin,depth(node) + eta_n(node))
           dmean2 = dmean**2
      
           do nz=1,nsigma
              btr(nz) = bt(nz,node)
           end do
      
           do nz=2,nsigma-1
              d1 = sigma(nz-1) - sigma(nz)
              d2 = sigma(nz)    - sigma(nz+1)
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
              b2(nz) = -0.5_WP*g*Kv(nz,node)*ap1/(density_0*dmean)
      
              !+++++++++++++++++++++++++++++++++++++++
              ! compute diffusion
              !+++++++++++++++++++++++++++++++++++++++
      
              snu_p = 0.5_WP*(snu(nz-1,node) + snu(nz,node))
              snu_m = 0.5_WP*(snu(nz,node) + snu(nz+1,node))
              c2(nz) = 0.5_WP*0.73_WP*( snu_p*(bt(nz-1,node) - bt(nz,node))/d1 - &
                                        snu_m*(bt(nz,node) - bt(nz+1,node))/d2 )/d
              c2(nz) = c2(nz)/dmean2
      
      
              if (b1(nz) > 0.0_WP) then
                  !if snu = 0 then b1 = 0
                  Ri(nz,node) = -g*ap1/(density_0*dmean)/b1(nz)*snu(nz,node)
              else
                  Ri(nz,node) = 10000.0_WP
              end if
      
      
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
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !u_din_sur = 0.05_WP*100_WP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
                 tke_dissip(nz,node) = c1
      
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
      
              if (bbb_max < 1.0d-6) exit  !go to 5
      
              !++++++++++++++++++++++++++++++++++++++++
              !  next iteration
              !++++++++++++++++++++++++++++++++++++++++
      
              do nz=1,nsigma
                 btr(nz) = btr1(nz)
              enddo
      
           enddo ! end for iteration
      
      !!!5    continue
      
           do nz=1,nsigma
      
              bt(nz,node) = max(btr(nz),1.0d-8)
      
              !++++++++++++++++++++++++++++++++++++++++
              ! scale of turbulence
              !++++++++++++++++++++++++++++++++++++++++
              !a = max(depth(node), Dmin)
              zt = sigma(nz)*dmean
              Zh = zt + z0b_min
              Zd = abs(dmean - zt) + z0s_min
              Z0 = 1.0_WP  - beta_scale*zh*zd/dmean**2
      
              !if (Ri(nz) > 0.0_WP) then
              !  abc = 1.0_WP - 5.0_WP*Ri(nz)
              !else
              !  abc = 1.0_WP
          !endif
              scale_t = cka*Zh*Zd*Z0/dmean
              snu(nz,node) = scale_t*sqrt(bt(nz,node)) + snul
              if (snu(nz,node) >= 0.2_WP) snu(nz,node) = 0.2_WP
              Kv(nz,node) = Pr_num*(snu(nz,node))
           end do
      
        enddo ! end for nodes

#ifdef USE_MPI
        call exchange_nod(Kv)
        call exchange_nod(snu)
        call exchange_nod(Ri)
        call exchange_nod(bt)
      
        call exchange_nod(tke_dissip)
#endif
      
        DO elem=1,myDim_elem2D
           elnodes=elem2D_nodes(:,elem)
           DO nz=1,nsigma
            Av(nz,elem)=sum(w_cv(1:4,elem)*snu(nz,elnodes))
           END DO
        END DO
      
#ifdef USE_MPI
        call exchange_elem(Av)
#endif
      
        deallocate(xx1, yy1, btr, btr1, b1, b2, c2)
      
#ifdef USE_MPI
      !!$     call MPI_AllREDUCE(minval(w_cv(1,:)),mnv, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
      !!$          MPI_COMM_FESOM_C, MPIerr)
      !!$
      !!$     call MPI_AllREDUCE(maxval(w_cv(1,:)),mxv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
      !!$          MPI_COMM_FESOM_C, MPIerr)
      !!$
      !!$     if (mype==0) print *,'MIN_VAL w_cv',mnv
      !!$     if (mype==0) print *,'MAX_VAL w_cv',mxv
      !!$
      !!$     call MPI_AllREDUCE(minval(Av(1,:)),mnv, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
      !!$          MPI_COMM_FESOM_C, MPIerr)
      !!$
      !!$     call MPI_AllREDUCE(maxval(Av(1,:)),mxv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
      !!$          MPI_COMM_FESOM_C, MPIerr)
      !!$
      !!$     if (mype==0) print *,'MIN_VAL Av',mnv
      !!$     if (mype==0) print *,'MAX_VAL Av',mxv
#else
      !     print *,'MIN_VAL w_cv',minval(w_cv(1,:))
      !     print *,'MAX_VAL w_cv',maxval(w_cv(1,:))
      !     print *,'MIN_VAL Av',minval(Av(1,:))
      !     print *,'MAX_VAL Av',maxval(Av(1,:))
#endif
      


    end subroutine d3_end


    SUBROUTINE densityJM(t, s, pz, rho_out,con_ss)

        use o_MESH
        use o_ARRAYS
        use o_PARAM
      
        implicit none 
      
        !
        ! - calculates in-situ density as a function of potential temperature
        !   (relative to the surface)
        !   using the Jackett and McDougall equation of state
        !   (Copyright (c) 1992, CSIRO, Australia)
        ! - has been derived from the SPEM subroutine rhocal
        !
        ! Ralph Timmermann, August 2005
        !---------------------------------------------------------------------------
      
      
        real(kind=WP), intent(IN)  :: t,s,pz
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
end module d3_module