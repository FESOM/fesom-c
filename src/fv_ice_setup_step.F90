
subroutine ice_setup()

    use i_param
    use i_arrays

    use i_therm_parms
    !use i_dyn_parms

    use o_mesh
    use o_arrays
    use o_param

    use g_parsup
    use o_UTILIT

    implicit none

    integer   :: err_alloc
    character(len=256)   :: nmlfile
    real(kind=WP)  :: delta_t
    integer            :: node_size, elem_size   ! number of nodes and elements on current PE (all one)

    node_size=myDim_nod2D+eDim_nod2D
    elem_size=myDim_elem2D+eDim_elem2D+eXDim_elem2D

    nmlfile ='namelist.ice'    ! name of ice namelist file
    open(20,file=nmlfile)
    read(20,NML=ice_dyn)
    read(20,NML=ice_therm)
    close(20)

    ! initializing sea ice model
    allocate(U_rhs_ice(2,elem_size), UiceAB(2,elem_size))
    !allocate(U_rhs_ice2(2,elem_size))
    allocate(U_rhs_ice_tmp(elem_size))    
    allocate(U_n_ice(2,elem_size), U_n_1ice(2,elem_size), U_n_2ice(2,elem_size))
    allocate(tau_ice(2,elem_size))
    allocate(grad_ssh_x(elem_size), grad_ssh_y(elem_size))
    allocate(Cdu_lf(elem_size), Cdv_lf(elem_size))


    allocate(delta_stress_e(elem_size),delta_stress_n(node_size),grad_delta_stress_e(elem_size))
    grad_delta_stress_e=0.0_WP
    delta_stress_e=0.0_WP
    delta_stress_n=0.0_WP
    U_rhs_ice=0.0_WP
    U_rhs_ice_tmp=0.0_WP
    UiceAB=0.0_WP
    U_n_ice=0.0_WP
    U_n_1ice=0.0_WP
    U_n_2ice=0.0_WP
    tau_ice=0.0_WP
    grad_ssh_x=0.0_WP
    grad_ssh_y=0.0_WP
    Cdu_lf = 0.0_WP
    Cdv_lf = 0.0_WP

    allocate(ice_pressure(elem_size))
    ice_pressure = 0.0_WP

    

    ! ice reology VPE
    allocate(sigma11(elem_size), sigma12(elem_size), sigma22(elem_size))
    sigma11=0.0_WP
    sigma12=0.0_WP
    sigma22=0.0_WP
    allocate(ice_delta(elem_size), ice_div(elem_size), ice_share(elem_size))
    ice_delta = 0.0_WP
    ice_div = 0.0_WP
    ice_share = 0.0_WP


    
!    allocate(Vel_nor(node_size))
!    allocate(tmp1(3,elem_size),tmp2(3,elem_size))
!    tmp1 = 0.0_WP
!    tmp2 = 0.0_WP
!    Vel_nor = 0.0_WP

    allocate(tmparr(3,node_size),tmparr_e(3,elem_size))
    tmparr = 0.0_WP
    tmparr_e = 0.0_WP

    allocate(cHrhsi(3,node_size),a_ice_RK(node_size), m_ice_RK(node_size), m_snow_RK(node_size))
    a_ice_RK = 0.0_WP
    m_ice_RK = 0.0_WP
    m_snow_RK = 0.0_WP
    cHrhsi = 0.0_WP
    U_rhs_ice = 0.0_WP
    UiceAB = 0.0_WP
    U_n_ice = 0.0_WP
    U_n_1ice = 0.0_WP
    U_n_2ice = 0.0_WP
    tau_ice = 0.0_WP
    allocate( a_ice(node_size), ice_temp(node_size), &
              m_ice(node_size), m_snow(node_size), &
              t_skin(node_size), &
                &      STAT=err_alloc )
!              u_ice(nod2D), v_ice(nod2D), &
    if( err_alloc /= 0 )   STOP 'ice_setup: failed to allocate arrays'
    a_ice = 0.0_WP
    ice_temp = 0.0_WP
    m_ice = 0.0_WP
    m_snow = 0.0_WP
 !   u_ice = 0.0_WP
 !   v_ice = 0.0_WP
    t_skin = 0.0_WP
    Ch_atm_oce_arr = Ch_atm_ice
    Ce_atm_oce_arr = Ce_atm_ice

! ice dynamics as function of eta in river (here is time as eta proxi)
    icedyn_tmask = 1.0_WP
    if (mask_icedyn) then
        icedyn_nmask=48
        allocate(icedyn_time(icedyn_nmask),icedyn_vmask(icedyn_nmask),STAT=err_alloc )
        if( err_alloc /= 0 )   STOP 'ice_setup: failed to allocate arrays'
!        icedyn_time=(/2451684.5, 2451690.5, 2451819.5, 2451849.5, 2452047.5, 2452057.5, &
!            2452198.5, 2452211.5, 2452418.5, 2452435.5, 2452567.5, 2452583.5, &
!            2452782.5, 2452795.5, 2452937.5, 2452950.5, 2453151.5, 2453166.5, &
!            2453295.5, 2453321.5, 2453513.5, 2453519.5, 2453666.5, 2453684.5, &
!            2453880.5, 2453887.5, 2454026.5, 2454040.5, 2454241.5, 2454261.5, &
!            2454398.5, 2454414.5, 2454625.5, 2454634.5, 2454771.5, 2454780.5, &
!            2454985.5, 2454997.5, 2455124.5, 2455145.5, 2455331.5, 2455345.5, &
!            2455484.5, 2455510.5, 2455699.5, 2455710.5, 2455854.5, 2455871.5/)
        icedyn_time=(/2451674.5, 2451690.5, 2451819.5, 2451859.5, 2452037.5, 2452057.5, &
            2452198.5, 2452221.5, 2452408.5, 2452435.5, 2452567.5, 2452573.5, &
            2452772.5, 2452795.5, 2452937.5, 2452960.5, 2453141.5, 2453166.5, &
            2453295.5, 2453331.5, 2453503.5, 2453519.5, 2453666.5, 2453694.5, &
            2453870.5, 2453887.5, 2454026.5, 2454050.5, 2454231.5, 2454261.5, &
            2454398.5, 2454424.5, 2454615.5, 2454634.5, 2454771.5, 2454790.5, &
            2454975.5, 2454997.5, 2455124.5, 2455155.5, 2455321.5, 2455345.5, &
            2455484.5, 2455520.5, 2455689.5, 2455710.5, 2455854.5, 2455881.5/)
        icedyn_vmask=(/0, 1, 1, 0, 0, 1, &
                       1, 0, 0, 1, 1, 0, &
                       0, 1, 1, 0, 0, 1, &
                       1, 0, 0, 1, 1, 0, &
                       0, 1, 1, 0, 0, 1, &
                       1, 0, 0, 1, 1, 0, &
                       0, 1, 1, 0, 0, 1, &
                       1, 0, 0, 1, 1, 0/)

        call binarysearch(icedyn_nmask,icedyn_time,time_jd,icedyn_t_indx)

        if ( icedyn_t_indx < icedyn_nmask ) then
         icedyn_t_indx_p1 = icedyn_t_indx + 1
         delta_t = icedyn_time(icedyn_t_indx_p1) - icedyn_time(icedyn_t_indx)
        else ! NO extrapolation to future
         icedyn_t_indx_p1 = icedyn_t_indx
         delta_t = 1.0_wp
         if (mype==0) then
            write(*,*) "     WARNING:  time lager than last time step in icedyn mask,"
            write(*,*) "        last time step will be used as a constant mask"
         endif
        end if
        if ( time_jd < icedyn_time(1) ) then  ! NO extrapolation back in time
         icedyn_t_indx_p1 = icedyn_t_indx
         delta_t = 1.0_wp
        end if
        ! linear time interpolation ,
        ! calculate new coefficients for interpolations
        ! delta_t was calculated before, in case of extrapolation delta_t=1
        icedyn_coefa = (icedyn_vmask(icedyn_t_indx_p1) - icedyn_vmask(icedyn_t_indx)) / delta_t
        icedyn_coefb = icedyn_vmask(icedyn_t_indx) - icedyn_coefa * icedyn_time(icedyn_t_indx)
    end if


end subroutine

subroutine ice_timestep()
    use i_arrays
    use i_param
    use o_mesh
    use o_arrays
    use o_param

    use g_parsup
    use g_comm_auto
    use o_UTILIT
    implicit none

    real(kind=WP) :: x1, max_m_ice, mx_tmp1,mx_tmp2,mx_tmp3
    integer       :: n, nRK4, nevp
    real(kind=WP)  :: delta_t

    call thermodynamics()

    if (mask_icedyn) then  ! if apply mask ice dynamic as function of eta/time
        if ( time_jd > icedyn_time(icedyn_t_indx_p1) ) then !if time to get new interpolation coefficients
            call binarysearch(icedyn_nmask,icedyn_time,time_jd,icedyn_t_indx)
            if ( icedyn_t_indx < icedyn_nmask ) then
                icedyn_t_indx_p1 = icedyn_t_indx + 1
                delta_t = icedyn_time(icedyn_t_indx_p1) - icedyn_time(icedyn_t_indx)
            else ! NO extrapolation to future
                icedyn_t_indx_p1 = icedyn_t_indx
                delta_t = 1.0_wp
            end if
            if ( time_jd < icedyn_time(1) ) then  ! NO extrapolation back in time
                icedyn_t_indx_p1 = icedyn_t_indx
                delta_t = 1.0_wp
            end if
            icedyn_coefa = (icedyn_vmask(icedyn_t_indx_p1) - icedyn_vmask(icedyn_t_indx)) / delta_t
            icedyn_coefb = icedyn_vmask(icedyn_t_indx) - icedyn_coefa * icedyn_time(icedyn_t_indx)
        end if
        icedyn_tmask = time_jd * icedyn_coefa + icedyn_coefb
     endif

#ifdef USE_MPI
    call MPI_AllREDUCE(maxval(m_ice),max_m_ice, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
             MPI_COMM_FESOM_C, MPIerr)
#else
    max_m_ice=maxval(m_ice)
#endif
    if (max_m_ice>0.05) then
         call EVPdynamics
    endif

end subroutine

!======================= Thermodynamic ========================

!
! For every surface node, this subroutine extracts the information
! needed for computation of thermodydnamics, calls the relevant
! subroutine, and returns the information to the vectors of prognostic
! variables.
!
! Originally programmed by N. Yakovlev/S. Danilov.
!
! Adjusted for upgraded model physics (NCEP forcing data; parameterization
! of ice-ocean heat flux considering friction velocity) by Ralph Timmermann.
!
! Adjusted for general forcing data and NlFs, 13.01.2009

! adjusted for FESOM-C from FESOM1.4
! by Ivan Kuznetsov, 15.12.2021
! adjusted for MPI , Ivan Kuznetsov 23.06.22
!
!----------------------------------------------------------------------------
subroutine thermodynamics()
  !
  ! For every surface node, this subroutine extracts the information
  ! needed for computation of thermodynamics, calls the relevant
  ! subroutine, and returns the information to the vectors of prognostic
  ! variables.
  !------------------------------------------------------------------------

    use i_arrays
    use i_param
    use i_therm_parms
    use o_param
    use o_mesh
    use o_arrays
#ifdef USE_MPI
    use fv_sbcmpi, only: atmdata, key_snow, nm_prec_coef     ! MPI version of surface boundary cond.
#else
    use fv_sbc, only: atmdata, key_snow, nm_prec_coef
#endif
    use g_parsup
    use g_comm_auto

    implicit none

    integer  :: n, elem, k
    real(kind=WP)  :: rsss ! Sea surface Salinity
    real(kind=WP)  :: ustar, uice, vice, tot_area
    real(kind=WP)  :: h, hsn, A, fsh, flo, Ta, qa, rain, snow, runo, ug, T_oc, S_oc
    real(kind=WP)  :: ch, ch_i, ce, ce_i, h_ml, fw, ehf, t, h00, evap, rsf
    real(kind=WP)  :: ithdgr, ithdgrsn, iflice, hflatow, hfsenow, hflwrdout,o2ihf,Sice0

!    real(kind=WP), allocatable  :: ustar_aux(:)

    ! u_ice and v_ice are at nodes
    ! u_w, v_w are at nodes (interpolated from elements)
    ! u_wind and v_wind are always at nodes
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n,h,hsn,A,fsh,flo,Ta,qa,rain,snow,runo,&
!$OMP& ug,tot_area,uice,vice,k,elem,ustar,T_oc,S_oc,rsss,t,ch,ce,ch_i,ce_i,h_ml,fw,&
!$OMP& ehf,h00,evap,rhowat,Sice0,rsf)

    do n=1, myDim_nod2D !myDim_nod2d+eDim_nod2D
        ! if there is a cavity no sea ice thermodynamics is apllied
        ! fesom_c (no cavity)
        ! if(ulevels_nod2d(i)>1) cycle fesom2
        ! if(cavity_flag_nod2d(i)==1) cycle fesom14

        h       = m_ice(n)     ! ice thickness [m]
        hsn     = m_snow(n)    ! snow thickness
        A       = a_ice(n)     ! ice concentation [%]
        fsh     = atmdata(4,n) ! atmdata from fv_sbc , i_qsr   = 4 ! index of solar heat         [W/m2]
        flo     = atmdata(5,n) ! atmdata from fv_sbc , i_qlw   = 5 ! index of Long wave          [W/m2]
        Ta      = atmdata(6,n) -273.15_WP! atmdata from fv_sbc , i_tair  = 6 ! index of 2m air temperature [degK] -273.15 [convert to C]
        qa      = atmdata(3,n) ! atmdata from fv_sbc , i_humi  = 3 ! index of specific humidity  [kg/kg]
!if (n==1) then
!    write(*,*) "ice: ",h,hsn,A,fsh,flo,Ta
!end if
        if (.not. key_snow) then
            if (Ta>=0.0_WP) then
                rain=atmdata(7,n)/nm_prec_coef/1000.0_WP !i_prec  = 7 ! index of total precipitation (rain+snow) [Kg/m^2/s] -> [m/s]
                snow=0.0_WP
            else
                rain=0.0_WP
                snow=atmdata(7,n)/nm_prec_coef/1000.0_WP !i_prec  = 7 ! index of total precipitation (rain+snow) [Kg/m^2/s] -> [m/s]
            endif
        else
            rain=atmdata(7,n)/nm_prec_coef/1000.0_WP !i_prec  = 7 ! index of total precipitation (rain+snow) [Kg/m^2/s] -> [m/s]
            STOP "snow is not implemented in sbc"
            snow=0.0_WP ! do not forget to /1000.0_WP
        end if
        runo    = 0.0_WP !runoff_rivers(n)
        ug      = sqrt(windx(n)**2+windy(n)**2)
        !convert uice on elements to nods
        tot_area=0.0_WP
        uice=0.0_WP
        vice=0.0_WP
        do k=1, nod_in_elem2D_num(n)
            elem=nod_in_elem2D(k,n)
            tot_area=tot_area+elem_area(elem)
            uice=uice+U_n_ice(1,elem)*elem_area(elem)
            vice=vice+U_n_ice(2,elem)*elem_area(elem)
        enddo
        uice=uice/tot_area
        vice=vice/tot_area
        ustar   = sqrt(Cd_oce_ice)*sqrt((uice-Unode(1,n))**2+(vice-Vnode(1,n))**2)
        T_oc    = TF(1,n)
        S_oc    = SF(1,n)
        !if(ref_sss_local) rsss = S_oc
        rsss = S_oc
        t       = t_skin(n)
        ch	     = Ch_atm_oce_arr(n)
        ce	     = Ce_atm_oce_arr(n)
        ch_i    = Ch_atm_ice
        ce_i    = Ce_atm_ice
        h_ml    = Jc(1,n)! 1.0_WP
        fw      = 0.0_WP
        ehf     = 0.0_WP
        h00=ice_h0
        evap = evaporation(n)

        rhowat = 1025.0_WP
        !IK we can use it for fresh body call densityJM(T_oc, S_oc, 0.0_WP, rhowat)
        !IK we can use it for fresh body   rhowat = rhowat + density_0 ! Volumetr. heat cap. of water [J/m**3/K](cc = rhowat*cp_water)
        cc = rhowat * 4190.0_WP
        Sice0 = Sice
        !IK Sice0 = S_oc / 10._WP

        call therm_ice(h,hsn,A,fsh,flo,Ta,qa,rain,snow,runo,rsss, &
                   &   ug,ustar,T_oc,S_oc,h_ml,t,ch,ce,ch_i,ce_i,fw,ehf,evap, &
                   &   rsf, ithdgr, ithdgrsn, iflice, hflatow, hfsenow, hflwrdout, h00, o2ihf,n,rhowat,Sice0)
!if (n==1) then
!    write(*,*) mype,"ice out: ",h,hsn,A,fsh,flo,Ta,t
!end if
        ! change fluxes due to ice presence
        qns(n) = qns(n)*(1.0_WP-A)+ehf! -o2ihf
        qsr(n) = qsr(n)*(1.0_WP-A)
!        emp(n) = (-1.0_WP)*fw*1000.0_WP
!if (n==1) write(*,*) qns(n),qsr(n),ehf
        m_ice(n)         = h
        m_snow(n)        = hsn
        a_ice(n)         = A
        t_skin(n)        = t
        fresh_water_flux(n) = fw
        net_heat_flux(n) = ehf
        evaporation(n)   = evap    !negative up
!        #ifdef use_fullfreesurf
        real_salt_flux(n)= 0.0_WP !rsf
!        #endif
!
!        thdgr(i)         = ithdgr ! growth rate ice
!        thdgrsn(i)       = ithdgrsn  ! growth rate snow
!        flice(i)         = iflice
!        olat_heat(i)     = hflatow
!        osen_heat(i)     = hfsenow
!        olwout(i)        = hflwrdout
    end do
!$OMP END PARALLEL DO
#ifdef USE_MPI
    call exchange_nod(qns)
    call exchange_nod(qsr)
    call exchange_nod(emp)
    call exchange_nod(m_ice)
    call exchange_nod(m_snow)
    call exchange_nod(a_ice)
    call exchange_nod(fresh_water_flux)
    call exchange_nod(net_heat_flux)
!    call exchange_nod(evaporation)
    call exchange_nod(real_salt_flux)
#endif
!    deallocate(ustar_aux)

end subroutine

subroutine therm_ice(h,hsn,A,fsh,flo,Ta,qa,rain,snow,runo,rsss, &
                    & ug,ustar,T_oc,S_oc,H_ML,t,ch,ce,ch_i,ce_i,fw,ehf,evap, &
                    & rsf, dhgrowth, dhsngrowth, iflice, hflatow, hfsenow, hflwrdout, h00, o2ihf, node,rhowat0,Sice0)
  ! Ice Thermodynamic growth model
  !
  ! Input parameters:
  !------------------
  ! h - ice mass [m]
  ! hsn - snow mass [m]
  ! A - ice compactness
  ! fsh - shortwave radiation
  ! flo - longwave radiation
  ! Ta - air temperature
  ! qa - specific humidity
  ! rain - precipitation rain
  ! snow - precipitation snow
  ! runo - runoff
  ! ug - wind speed
  ! ustar - friction velocity
  ! T_oc, S_oc - ocean temperature and salinity beneath the ice (mixed layer)
  ! H_ML - mixed layer depth - should be specified.
  ! t - temperature of snow/ice top surface
  ! dt - time step [s]
  ! ch - transfer coefficient for sensible heat (for open ocean)
  ! ce - transfer coefficient for evaporation   (for open ocean)
  ! ch_i - transfer coefficient for sensible heat (for ice)
  ! ce_i - transfer coefficient for evaporation   (for ice)

  ! Output parameters:
  !-------------------
  ! h - ice mass
  ! hsn - snow mass
  ! A - ice compactness
  ! t - temperature of snow/ice top surface
  ! fw - freshwater flux due to ice melting [m water/dt]
  ! ehf - net heat flux at the ocean surface [W/m2]        !RTnew
  use i_therm_parms

  implicit none

  integer k, node
  real(kind=WP) :: h,hsn,A,fsh,flo,Ta,qa,rain,snow,runo,rsss,rhowat0,Sice0
  real(kind=WP) :: ug,ustar,T_oc,S_oc,H_ML,t,ch,ce,ch_i,ce_i,fw,ehf
  real(kind=WP) :: dhgrowth,dhsngrowth,ahf,prec,subli,subli_i,rsf
  real(kind=WP) :: rhow,show,rhice,shice,sh,thick,thact
  real(kind=WP) :: rh,rA,qhst,sn,hsntmp,o2ihf,evap
  real(kind=WP) :: iflice,hflatow,hfsenow,hflwrdout,h00
  real(kind=WP), external  :: TFrez  ! Sea water freeze temperature.

  ! Store ice thickness at start of growth routine
  dhgrowth=h

  ! determine h(i,j)/a(i,j) = actual ice thickness.
  ! if snow layer is present, add hsn weighted with quotient
  ! of conductivities of ice and snow, according to 0-layer approach
  ! of Semtner (1976).
  ! thickness at the ice covered part
  thick=hsn*(con/consn)/max(A,Armin)    ! Effective snow thickness
  thick=thick+h/max(A,Armin)            ! Effective total snow-ice thickness

!#if defined (__oasis) || defined(__uncplecham6)
!   rhow =-fsh/cl
!   rhice=-flo/cl
!   subli=ta
!#else
  ! Growth rate for ice in open ocean
  rhow=0.0_WP
  evap=0.0_WP
  call obudget(qa,fsh,flo,T_oc,ug,ta,ch,ce,rhow,evap,hflatow,hfsenow,hflwrdout,rhowat0)
  hflatow=hflatow*(1.0_WP-A)
  hfsenow=hfsenow*(1.0_WP-A)
  hflwrdout=hflwrdout*(1.0_WP-A)

  ! add heat loss at open ocean due to melting snow fall
  !rhow=rhow+snow*1000.0/rhoice !qiang
  ! dt and (1-A) will be multiplied afterwards

  ! growth rate of ice in ice covered part
  ! following Hibler 1984
  ! assuming ice thickness has an euqal, 7-level distribution from zero to two times h
  rhice=0.0_WP
  subli=0.0_WP
  if (thick.gt.hmin) then
     do k=1,iclasses
        thact = (2.0_WP*k-1.0_WP)*thick/float(iclasses)  	! Thicknesses of actual ice class
        call budget(thact,hsn,t,ta,qa,fsh,flo,ug,S_oc,ch_i,ce_i,shice,subli_i,rhowat0)
        !Thick ice K-class growth rate
        rhice=rhice+shice/float(iclasses)      	! Add to average heat flux
        subli=subli+subli_i/float(iclasses)
     end do
  end if
!#endif
  ! Convert growth rates [m ice/sec] into growth per time step DT.
  rhow=rhow*dt
  rhice=rhice*dt

  ! Multiply ice growth of open water and ice
  ! with the corresponding areal fractions of grid cell
  show =rhow*(1.0-A)
  shice=rhice*A
  sh   =show+shice

  ! Store atmospheric heat flux, average over grid cell [W/m**2]
  ahf=-cl*sh/dt

  ! precipitation (into the ocean)
  prec=rain+runo+snow*(1.0-A)  	        ! m water/s

  ! snow fall above ice
  hsn=hsn+snow*dt*A*1000.0/rhosno	! Add snow fall to temporary snow thickness    !!!
  dhsngrowth=hsn   		        ! Store snow thickness after snow fall

  ! evap
  evap=evap*(1.0-A)    			! m water/s
  subli=subli*A

  ! If there is atmospheric melting, first melt any snow that is present.
  ! Atmospheric heat flux available for melting
  ! sh = MINUS atm. heat flux / specific latent heat of sea ice
  ! Note: (sh<0) for melting, (sh>0) for freezing
  hsntmp= -min(sh,0.0)*rhoice/rhosno

  ! hsntmp is the decrease in snow thickness due to atmospheric melting
  ! [m/DT]. Do not melt more snow than available
  hsntmp=min(hsntmp,hsn)
  hsn=hsn-hsntmp  ! Update snow thickness after atmospheric snow melt

  ! Negative atmospheric heat flux left after melting of snow
  ! Note: (sh<0) and (hsntmp>0) for melting conditions
  ! hsntmp=0 for non-snow-melting conditions
  rh=sh+hsntmp*rhosno/rhoice
  h=max(h,0.)

  ! Compute heat content qhst of mixed layer - sea ice system
  !
  ! Total heat content is the sum of
  !	h	ice thickness after calculation of dynamic effects
  !	rh	change in ice thickness due to atmospheric forcing
  ! and heat available in mixed layer, with
  !	T_oc	temperature of ocean surface layer
  !	Tfrez	freezing point of sea water
  !	H_ML	thickness of uppermost layer
  !
  !RT:
  ! There are three possibilities to do this.
  ! 1.: Assume an instantaneous adjustment of mixed layer heat content.
  !     Any heat available is then instantaneously used to melt ice.
  !     (so-called ice-bath approach)
  !     This is what used to be used in the Lemke sea ice-mixed layer model.
  ! rh=rh-(T_oc-TFrez(S_oc))*H_ML*cc/cl
  ! qhst=h+rh
  !
  ! 2.: Parameterize the ocean-to-ice heat flux (o2ihf)
  !     as a function of temperature difference. For a first step
  !     we can assume a constant exchange coefficient gamma_t:
  ! o2ihf= (T_oc-TFrez(S_oc))*gamma_t*cc*A     &
  !        +(T_oc-Tfrez(S_oc))*H_ML/dt*cc*(1.0-A) ! [W/m2]
  ! rh=rh-o2ihf*dt/cl
  ! qhst=h+rh		                      	! [m]
  !
  ! 3.  Parameterize the ocean-to-ice heat flux (o2ihf)
  !     as a function of temperature difference and the
  !     friction velocity:
  o2ihf= (T_oc-TFrez(S_oc))*0.006*ustar*cc*A  &
       +(T_oc-Tfrez(S_oc))*H_ML/dt*cc*(1.0-A)  	! [W/m2]
  !o2ihf=min(0.0_WP,o2ihf)
  !if (node==1) then
    !write(*,*) T_oc,Tfrez(S_oc)
    !write(*,*) o2ihf, rh, o2ihf*dt/cl, (T_oc-TFrez(S_oc))*0.006*ustar*cc*A,(T_oc-Tfrez(S_oc))*H_ML/dt*cc*(1.0-A)
  !end if
  rh=rh-o2ihf*dt/cl
  qhst=h+rh		              		! [m]

  ! Melt snow if there is any ML heat content left (qhst<0).
  ! This may be the case if advection moves ice (with snow) to regions
  ! with a warm mixed layer.
  sn=hsn+min(qhst,0.)*rhoice/rhosno

  ! New temporary snow thickness must not be negative:
  sn=max(sn,0.)

  ! Update snow and ice depth
  hsn=sn
  h=max(qhst,0.)
  if (h.lt.1E-6) h=0.        ! Avoid very small ice thicknesses

  ! heat and fresh water fluxes
  dhgrowth=h-dhgrowth        ! Change in ice thickness due to thermodynamic effects
  dhsngrowth=hsn-dhsngrowth  ! Change in snow thickness due to thermodynamic melting
  ! (without snow fall). This is a negative value (MINUS snow melt)

  dhgrowth=dhgrowth/dt       ! Conversion: 'per time step' -> 'per second'
  dhsngrowth=dhsngrowth/dt   ! Conversion: 'per time step' -> 'per second'

  ! (radiation+turbulent) + freezing(-melting) sea-ice&snow
  ehf = ahf + cl*(dhgrowth+(rhosno/rhoice)*dhsngrowth)
  ! (prec+runoff)+evap - freezing(+melting) ice&snow

!#ifdef use_fullfreesurf
  fw= prec+evap - dhgrowth*rhoice/rhowat0 - dhsngrowth*rhosno/rhowat0
  rsf= -dhgrowth*rhoice/rhowat0*Sice0
!  rsf = 0.0_WP
!#else
!  fw= prec+evap - dhgrowth*rhoice/rhowat*(rsss-Sice)/rsss - dhsngrowth*rhosno/rhowat
!#endif

  ! Changes in compactnesses (equation 16 of Hibler 1979)
  rh=-min(h,-rh)   ! Make sure we do not try to melt more ice than is available
  rA= rhow - o2ihf*dt/cl !Qiang: it was -(T_oc-TFrez(S_oc))*H_ML*cc/cl, changed in June 2010
  !rA= rhow - (T_oc-TFrez(S_oc))*H_ML*cc/cl*(1.0-A)
  A=A + 0.5*min(rh,0.)*A/max(h,hmin) + max(rA,0.)*(1.-A)/h00
  !meaning:       melting                  freezing

  A=min(A,h*1.e6)     ! A -> 0 for h -> 0
  A=min(max(A,0.),1.) ! A >= 0, A <= 1

  ! Flooding (snow to ice conversion)
!  iflice=h
!  call flooding(h,hsn,rhowat0)
!  iflice=(h-iflice)/dt

  ! to maintain salt conservation for the current model version
  !(a way to avoid producing net salt from snow-type-ice)

!#ifdef use_fullfreesurf
!  rsf=rsf-iflice*rhoice/rhowat0*Sice0
  rsf = 0.0_WP
!#else
!  fw=fw+iflice*rhoice/rhowat*Sice/rsss
!#endif

  ! add sublimation to evap
  evap=evap+subli  !negative up, m/s


  return
end subroutine therm_ice
!
!=====================================================================================
!
subroutine budget (hice,hsn,t,ta,qa,fsh,flo,ug,S_oc,ch_i,ce_i,fh,subli,rhowat0)
  ! Thick ice growth rate [m ice/sec]
  !
  ! INPUT:
  ! hice - actual ice thickness [m]
  ! hsn - snow thickness, used for albedo parameterization [m]
  ! t - temperature of snow/ice surface [C]
  ! ta - air temperature [C]
  ! qa - specific humidity [Kg/Kg]
  ! fsh - shortwave radiation [W/m**2]
  ! flo - longwave radiation  [W/m**2]
  ! ug - wind speed [m/sec]
  ! S_oc - ocean salinity for the temperature of the ice base calculation [ppt]
  ! ch_i - transfer coefficient for sensible heat (for ice)
  ! ce_i - transfer coefficient for evaporation   (for ice)
  !
  ! OUTPUT: fh - growth rate
  !
  ! qiang: The formular for saturated humidity was modified according to Large/Yeager2004
  ! to allow direct comparison with the CORE results (Griffies et al. 2009). The new
  ! formular does not require sea level pressure.
  ! A similar change was also made for the obudget routine.
  ! It was found through experiments that the results are quite similar to that from the
  ! original code, and the simulated ice volume is only slightly larger after modification.

  use i_therm_parms
  implicit none

  integer iter, imax      ! Number of iterations
  real*8  hice,hsn,t,ta,qa,fsh,flo,ug,S_oc,ch_i,ce_i,fh
  real*8  hfsen,hfrad,hflat,hftot,subli,rhowat0
  real*8  alb             ! Albedo of sea ice
  real*8  q1, q2	  ! coefficevapients for saturated specific humidity
  real*8  A1,A2,A3,B,C, d1, d2, d3
  real*8, external :: TFrez

  data q1 /11637800.0/, q2 /-5897.8/
  data imax /5/

  ! set albedo
  ! ice and snow, freezing and melting conditions are distinguished.
  if (t<0.0) then	        ! freezing condition
     if (hsn.gt.0.0) then	!   snow cover present
        alb=albsn
     else              		!   no snow cover
        alb=albi
     endif
  else			        ! melting condition
     if (hsn.gt.0.0) then	!   snow cover present
        alb=albsnm
     else			!   no snow cover
        alb=albim
     endif
  endif

  d1=rhoair*cpair*Ch_i
  d2=rhoair*Ce_i
  d3=d2*clhi

  ! total incoming atmospheric heat flux
  A1=(1.0-alb)*fsh + flo + d1*ug*ta + d3*ug*qa   ! in LY2004 emiss is multiplied wiht flo

  ! NEWTON-RHAPSON TO GET TEMPERATURE AT THE TOP OF THE ICE LAYER

  do iter=1,imax

     B=q1/rhoair*exp(q2/(t+tmelt))		! (saturated) specific humidity over ice
     A2=-d1*ug*t-d3*ug*B &
          -emiss_ice*boltzmann*((t+tmelt)**4)	! sensible and latent heat,and outward radiation
     A3=-d3*ug*B*q2/((t+tmelt)**2)		! gradient coefficient for the latent heat part
     C=con/hice                     		! gradient coefficient for downward heat conductivity
     A3=A3+C+d1*ug & 			! gradient coefficient for sensible heat and radiation
          +4.0*emiss_ice*boltzmann*((t+tmelt)**3)
     C=C*(TFrez(S_oc)-t)       		! downward conductivity term

     t=t+(A1+A2+C)/A3 		        ! NEW ICE TEMPERATURE AS THE SUM OF ALL COMPONENTS
  end do

  t=min(0.0,t)

  ! heat fluxes [W/m**2]:

  hfrad= (1.0-alb)*fsh &	        ! absorbed short wave radiation
       +flo &           	        ! long wave radiation coming in  ! in LY2004 emiss is multiplied
       -emiss_ice*boltzmann*((t+tmelt)**4) 	! long wave radiation going out

  hfsen=d1*ug*(ta-t) 			! sensible heat
  subli=d2*ug*(qa-B) 			! sublimation
  hflat=clhi*subli                     	! latent heat

  hftot=hfrad+hfsen+hflat               ! total heat

  fh= -hftot/cl                         ! growth rate [m ice/sec]
  !                                      	+: ML gains energy, ice melts
  !                                      	-: ML loses energy, ice grows
  subli=subli/rhowat0                    ! negative upward

  return
end subroutine budget
!
!======================================================================================
!
subroutine obudget (qa,fsh,flo,t,ug,ta,ch,ce,fh,evap,hflatow,hfsenow,hflwrdout,rhowat0)
  ! Ice growth rate for open ocean [m ice/sec]
  !
  ! INPUT:
  ! t - temperature of open water [C]
  ! fsh - shortwave radiation
  ! flo - longwave radiation
  ! ta - air temperature [C]
  ! qa  - specific humidity
  ! ug - wind speed [m/sec]
  ! ch - transfer coefficient for sensible heat
  ! ce - transfer coefficient for evaporation
  !
  ! OUTPUT: fh - growth rate
  !         evap - evaporation

  use i_therm_parms
  implicit none

  real*8 qa,t,Ta,fsh,flo,ug,ch,ce,fh,evap,rhowat0
  real*8 hfsenow,hfradow,hflatow,hftotow,hflwrdout,b
  real*8 q1, q2 		! coefficients for saturated specific humidity
  real*8 c1, c4, c5
  logical :: standard_saturation_shum_formula=.false.

  data c1, c4, c5 /3.8e-3, 17.27, 237.3/
  data q1 /640380./, q2 /-5107.4/

  ! (saturated) surface specific humidity
  if(standard_saturation_shum_formula) then
     b=c1*exp(c4*t/(t+c5))                      ! a standard one
  else
     b=0.98*q1/rhoair*exp(q2/(t+tmelt))         ! LY2004 NCAR version
  end if

  ! heat fluxes [W/m**2]:

  hfradow= (1.0-albw)*fsh &	            	! absorbed short wave radiation
       +flo             	                ! long wave radiation coming in !put emiss/check
  hflwrdout=-emiss_wat*boltzmann*((t+tmelt)**4) 	! long wave radiation going out !in LY2004 emiss=1
  hfradow=hfradow+hflwrdout

  hfsenow=rhoair*cpair*ch*ug*(ta-t)             ! sensible heat
  evap=rhoair*ce*ug*(qa-b)			! evaporation kg/m2/s
  hflatow=clhw*evap                       	! latent heat W/m2

  hftotow=hfradow+hfsenow+hflatow             	! total heat W/m2

  fh= -hftotow/cl                             	! growth rate [m ice/sec]
  !                                           	+: ML gains energy, ice melts
  !                                           	-: ML loses energy, ice grows
  evap=evap/rhowat0	 			! evaporation rate [m water/s],negative up !!!

  return
end subroutine obudget
!
!======================================================================================
!
subroutine flooding (h,hsn,rhowat0)
  use i_therm_parms

  real*8 h,hsn,hdraft,hflood,rhowat0

  hdraft=(rhosno*hsn+h*rhoice)/rhowat0 ! Archimedes: displaced water
  hflood=hdraft-min(hdraft,h)         ! Increase in mean ice thickness due to flooding
  h=h+hflood                          ! Add converted snow to ice volume
  hsn=hsn-hflood*rhoice/rhosno        ! Subtract snow from snow layer

  !RT   This is what all AWI sea ice models do, but
  !RT   I wonder whether it really is correct for the heat budget.
  !RT   I suggest we initially keep it to allow for a comparison with BRIOS results
  !RT   and rethink it at a later stage.

  return
end subroutine flooding
!
!======================================================================================
!
function TFrez(S)
  use o_PARAM
  ! Nonlinear correlation for the water freezing temperature.
  ! Millero (1978) - UNESCO. Reference - See A. Gill, 1982.
  implicit none
  real(kind=WP) S, TFrez

  TFrez= -0.0575*S+1.7105e-3 *sqrt(S**3)-2.155e-4 *S*S

end function TFrez

!
!======================================================================================
subroutine stress_tensor
! EVP rheology. The routine computes stress tensor components based on ice
! velocity field. They are stored as elemental arrays (sigma11, sigma22 and
! sigma12).

    use i_arrays
    use i_param
    use i_therm_parms
    use o_mesh
    use o_arrays
    use o_param

    use g_parsup
    use g_comm_auto

    implicit none

    real(kind=WP)   :: eps11, eps12, eps22, eta, xi, ice_strength, delta, aa
    real(kind=WP)   :: asum, msum, vale, dx(3), dy(3), eps1, eps2, pressure
    real(kind=WP)   :: det1, det2, r1, r2, r3, si1, si2, dte
    real(kind=WP)   :: zeta, delta_inv, Tevp_inv

!
    real(kind=WP)   :: mx_sigma11
!


    integer         :: elem, elnodes(4), nRK4, i

    vale = 1.0_WP/(ellipse**2)

!
! time scale of transition from elastic behavior to the VP reology
!

    det2 = 1.0_WP/(1.0_WP + alpha_evp)
    det1 = alpha_evp*det2

    call vel_gradients_2D(UiceAB)
    do elem=1,myDim_elem2D
        elnodes = elem2D_nodes(:,elem)
!+++++++++++++++++++++++++++++++++++++++++++
!  ice strength (Hunke and Dukowicz c*h*p*)
!+++++++++++++++++++++++++++++++++++++++++++
        msum = w_cv(1,elem)*m_ice(elnodes(1)) + w_cv(2,elem)*m_ice(elnodes(2)) &
             + w_cv(3,elem)*m_ice(elnodes(3)) + w_cv(4,elem)*m_ice(elnodes(4))
        asum = w_cv(1,elem) *a_ice(elnodes(1)) + w_cv(2,elem) *a_ice(elnodes(2)) &
             + w_cv(3,elem) *a_ice(elnodes(3)) + w_cv(4,elem) *a_ice(elnodes(4))

!        asum = w_cv(1,elem) *exp(-c_pressure*(1.0_WP - a_ice(elnodes(1)))) &
!             + w_cv(2,elem) *exp(-c_pressure*(1.0_WP - a_ice(elnodes(2)))) &
!             + w_cv(3,elem) *exp(-c_pressure*(1.0_WP - a_ice(elnodes(3)))) &
!             + w_cv(4,elem) *exp(-c_pressure*(1.0_WP - a_ice(elnodes(4))))


        ice_strength = Pstar*(msum)*asum!exp(-c_pressure*(1.0_WP - asum))

!+++++++++++++++++++++++++++++++++++++++++++
! deformation rate tensor on element elem
!+++++++++++++++++++++++++++++++++++++++++++

        eps11 = vel_grad(1,elem) - UiceAB(2,elem)*metric_factor(elem)
        eps22 = vel_grad(4,elem)
        eps12 = 0.5_WP*(vel_grad(2,elem) + vel_grad(3,elem) + UiceAB(1,elem)*metric_factor(elem))

!+++++++++++++++++++++++++++++++++++++++++++
! switch to eps1,eps2
!+++++++++++++++++++++++++++++++++++++++++++
       eps1 = eps11 + eps22
       eps2 = eps11 - eps22

!+++++++++++++++++++++++++++++++++++++++++++
! moduli
!+++++++++++++++++++++++++++++++++++++++++++
       delta  =eps1**2 + vale*(eps2**2 + 4.0_WP*eps12**2)
       delta = sqrt(delta)
 !aa150623      delta_stress_e(elem) = delta
       pressure = ice_strength/(delta + delta_min)
       !dyagnostic
       ice_pressure(elem) = pressure

       r1 = pressure*(eps1 - delta) 
       r2 = pressure*eps2*vale
       r3 = pressure*eps12*vale
       si1 = sigma11(elem) + sigma22(elem)    
       si2 = sigma11(elem) - sigma22(elem)

       si1 = det1*si1 + det2*r1
       si2 = det1*si2 + det2*r2
       sigma12(elem) = det1*sigma12(elem) + det2*r3     
       sigma11(elem) = 0.5_WP*(si1 + si2)                 
       sigma22(elem) = 0.5_WP*(si1 - si2)  

       !dyagnostic
       ice_delta(elem) = delta
       ice_div(elem) = (eps11+eps12)*0.5_WP
       ice_share(elem) = sqrt(((eps11-eps22)*0.5_WP)**2+eps12**2)

    end do

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Equations solved in terms of si1, si2, eps1, eps2 are (43)-(45) of 
! Boullion et al Ocean Modelling 2013, but in an implicit mode:
! si1_{p+1}=det1*si1_p+det2*r1, where det1=alpha/(1+alpha) and det2=1/(1+alpha),
! and similarly for si2 and sigma12
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef USE_MPI
    call exchange_elem(delta_stress_e)
    call exchange_elem(sigma12)
    call exchange_elem(sigma11)
    call exchange_elem(sigma22)

    call exchange_elem(ice_delta)
    call exchange_elem(ice_div)
    call exchange_elem(ice_share)  ! ice share sqrt(((eps11-eps22)*0.5)^2+eps12^2) (for output of stress )        

    call exchange_elem(ice_pressure)


#endif

!aa150623   call compute_el2nodes_2D(delta_stress_e,delta_stress_n)
!aa150623    do elem=1,myDim_elem2D
!aa150623        elnodes=elem2D_nodes(:,elem)
!aa150623        grad_delta_stress_e(elem) =     sum(gradient_sca(1:4,elem)*delta_stress_e(elnodes(1)) &
!aa150623                              + gradient_sca(2,elem)*delta_stress_e(elnodes(2)) &
!aa150623                              + gradient_sca(3,elem)*delta_stress_e(elnodes(3)) &
!aa150623                              + gradient_sca(4,elem)*delta_stress_e(elnodes(4)) )

!aa150623    end do


!aa150623 #ifdef USE_MPI
!aa150623    call exchange_elem(grad_delta_stress_e)
!aa150623 #endif    


end subroutine stress_tensor

subroutine stress2rhs
    ! EVP implementation:
    ! Computes the divergence of stress tensor and puts the result into the tmp rhs vectors 
        use i_arrays
        use i_param
        use i_therm_parms
        use o_mesh
        use o_arrays
        use o_param
    
        use g_parsup
        use g_comm_auto
    IMPLICIT NONE
    
        INTEGER      :: elem, ed, elnodes(4), el(2), edglim, nodes(2)
        REAL(kind=WP) :: mass, uc, vc, xe, ye, a, b, acc, c, inv_mass

#ifdef USE_MPI
        edglim=myDim_edge2D+eDim_edge2D
#else
        edglim=edge2D
#endif
        do ed=1, edglim
    
            nodes=edge_nodes(:,ed)
            acc = sum(ac(nodes))/2.0_WP
            el = edge_tri(:,ed)
            xe = edge_dxdy(1,ed)
            ye = edge_dxdy(2,ed)         ! xe, ye are in radians!
    
            if (el(2) > 0) then
            !if (myList_edge2D(ed)<=edge2D_in) then    
                uc = 0.5_WP*r_earth*(ye*(sigma11(el(1)) + sigma11(el(2))) - &
                            xe*(elem_cos(el(1))*sigma12(el(1)) + &
                            elem_cos(el(2))*sigma12(el(2))))
                vc = 0.5_WP*r_earth*(ye*(sigma12(el(1)) + sigma12(el(2))) - &
                            xe*(elem_cos(el(1))*sigma22(el(1)) + &
                            elem_cos(el(2))*sigma22(el(2))))
                U_rhs_ice(1,el(1)) = U_rhs_ice(1,el(1)) + uc !*acc
                U_rhs_ice(1,el(2)) = U_rhs_ice(1,el(2)) - uc !*acc
                U_rhs_ice(2,el(1)) = U_rhs_ice(2,el(1)) + vc !*acc
                U_rhs_ice(2,el(2)) = U_rhs_ice(2,el(2)) - vc !*acc
            else
                uc = r_earth*(ye*sigma11(el(1)) - xe*elem_cos(el(1))*sigma12(el(1)))
                vc = r_earth*(ye*sigma12(el(1)) - xe*elem_cos(el(1))*sigma22(el(1)))
                U_rhs_ice(1,el(1)) = U_rhs_ice(1,el(1)) + uc !*acc
                U_rhs_ice(2,el(1)) = U_rhs_ice(2,el(1)) + vc !     *acc
            end if
    
        enddo
#ifdef USE_MPI
        call exchange_elem(U_rhs_ice) 
#endif

    do elem=1, myDim_elem2D
        elnodes=elem2D_nodes(:,elem)
        a = max( hmin, w_cv(1,elem) *m_ice(elnodes(1)) + w_cv(2,elem) *m_ice(elnodes(2)) &
                     + w_cv(3,elem) *m_ice(elnodes(3)) + w_cv(4,elem) *m_ice(elnodes(4)) )
        b = w_cv(1,elem) *m_snow(elnodes(1)) + w_cv(2,elem) *m_snow(elnodes(2)) &
                     + w_cv(3,elem) *m_snow(elnodes(3)) + w_cv(4,elem) *m_snow(elnodes(4))
        c = w_cv(1,elem) *a_ice(elnodes(1)) + w_cv(2,elem) *a_ice(elnodes(2)) &
                     + w_cv(3,elem) *a_ice(elnodes(3)) + w_cv(4,elem) *a_ice(elnodes(4))

        inv_mass = 1.0_WP/max(0.1_WP,rhoice*a + rhosno*b )    

        if (mass > 0.1_WP) then
            U_rhs_ice(1,elem) = U_rhs_ice(1,elem)*inv_mass/elem_area(elem) + grad_ssh_x(elem)
            U_rhs_ice(2,elem) = U_rhs_ice(2,elem)*inv_mass/elem_area(elem) + grad_ssh_y(elem)
        else
            U_rhs_ice(1,elem) = 0.0_WP
            U_rhs_ice(2,elem) = 0.0_WP
        end if
    enddo
#ifdef USE_MPI
    call exchange_elem(U_rhs_ice)
#endif


    end subroutine stress2rhs
    

!==============================
! SSH gradient for ice dynamic
! 16.05.2023 AA, IK
!==============================
SUBROUTINE gradssh_2_RHice
    use i_PARAM
    use i_ARRAYS
    use i_therm_parms

    use o_MESH
    use o_ARRAYS
    use o_PARAM

    use g_parsup
    use g_comm_auto

    implicit none

    integer               :: elem, elnodes(4)

!$OMP DO
    do elem=1,myDim_elem2D

        elnodes=elem2D_nodes(:,elem)
        grad_ssh_x(elem) =     g*(  gradient_sca(1,elem)*eta_n(elnodes(1)) &
                                  + gradient_sca(2,elem)*eta_n(elnodes(2)) &
                                  + gradient_sca(3,elem)*eta_n(elnodes(3)) &
                                  + gradient_sca(4,elem)*eta_n(elnodes(4)) )

        grad_ssh_y(elem) =     g*(  gradient_sca(5,elem)*eta_n(elnodes(1)) &
                                  + gradient_sca(6,elem)*eta_n(elnodes(2)) &
                                  + gradient_sca(7,elem)*eta_n(elnodes(3)) &
                                  + gradient_sca(8,elem)*eta_n(elnodes(4)) )

    end do
!$OMP END DO

#ifdef USE_MPI
    call exchange_elem(grad_ssh_x)
    call exchange_elem(grad_ssh_y)
#endif

END SUBROUTINE gradssh_2_RHice

!--------------------------------------------------------------------
!====================================================================
! Landfast ice. Friction. No-slip boundary conditions.
! 06.07.23
! see: Lemieux et al., 2015. A basal stress parametrization
!                            for modeling landfast ice.
! ===================================================================

subroutine landfast_friction_NS

    use i_PARAM
    use i_ARRAYS
    use i_therm_parms

    use o_MESH
    use o_ARRAYS
    use o_PARAM

    use g_parsup
    use g_comm_auto

    implicit none

  real(kind=WP)  :: dmean, a, b, c, hl, rhowat_lf, hd, hrb, hrt, hh, hc, k2_lf, UiceM, modUVice_inv
  integer        :: elem, elnodes(4)

        rhowat_lf = 1025.0_WP

    do elem=1,myDim_elem2D

        elnodes=elem2D_nodes(:,elem)

        b = max( hmin, w_cv(1,elem) *m_ice(elnodes(1)) + w_cv(2,elem) *m_ice(elnodes(2)) &
                     + w_cv(3,elem) *m_ice(elnodes(3)) + w_cv(4,elem) *m_ice(elnodes(4)) )
        c = w_cv(1,elem) *a_ice(elnodes(1)) + w_cv(2,elem) *a_ice(elnodes(2)) &
                     + w_cv(3,elem) *a_ice(elnodes(3)) + w_cv(4,elem) *a_ice(elnodes(4))


      if (c > 0.01_WP) then
     dmean =  max( Dmin, w_cv(1,elem)*(eta_n(elnodes(1)) + depth(elnodes(1))) &
                     +   w_cv(2,elem)*(eta_n(elnodes(2)) + depth(elnodes(2))) &
                     +   w_cv(3,elem)*(eta_n(elnodes(3)) + depth(elnodes(3))) &
                     +   w_cv(4,elem)*(eta_n(elnodes(4)) + depth(elnodes(4)))  )


!   a = 2.0_WP - 2.0_WP*beta_lf + beta_lf*q_lf * beta_lf*q_lf*hama_lf
 !  hl = 2.0_WP*(b/c)/a   ! Most ice thickness
 !  hd = rhoice*hl/rhowat_lf         ! rhowat_lf - density of the surface water can be VARIABLE! hd - part of ice thickness in the water
 !  hrb = (hama_lf - 1.0_WP)*rhoice*hl/rhowat_lf   ! bottom ridge of ice
!   hrt = (hama_lf - 1.0_WP)*(rhowat_lf-rhoice)*hl/rhowat_lf   ! top ridge of ice
!   hh = c*hl*(1.0_WP - beta_lf + beta_lf*q_lf) + 0.5*c*beta_lf*q_lf*(hrt + hrb) ! ice in a grid cell
!   if (hh > dmean) hh = dmean   ! limit for the ice grid cell thickness
!    hc = dmean*c*rhowat_lf*a/(2.0_WP*rhoice*hama_lf)
!   hc = hrb + hd  ! critical thickness of the ice
    hc = dmean*c/20.0_WP
 !  if (hc > dmean) hc = dmean   ! limit for the ice thickness
   k2_lf = 15.0_WP !mu_lf*(2.0_WP*beta_lf*q_lf*hama_lf*rhoice*g)/a
   UiceM = sqrt(U_n_ice(1,elem)**2 + U_n_ice(2,elem)**2 )


       if (b > hc .and. UiceM > 0.01) then
      
     !   UiceM = sqrt(U_n_ice(1,elem)**2 + U_n_ice(2,elem)**2 + UiceM0**2)
    
        modUVice_inv = 1.0_WP/UiceM
        if (modUVice_inv>10.0) then
                modUVice_inv = 10.0
        endif        
    !    Cdu_lf(elem) = modUVice_inv*k2_lf*U_n_ice(1,elem)*(b - hc)*exp(-c_pressure*(1.0_WP - c))
    !    Cdv_lf(elem) = modUVice_inv*k2_lf*U_n_ice(2,elem)*(b - hc)*exp(-c_pressure*(1.0_WP - c))
        Cdu_lf(elem) = 0.1_WP*modUVice_inv*k2_lf*LOG(exp(10._WP*(b - hc)+1.0_WP))*exp(-c_pressure*(1.0_WP - c))
        Cdv_lf(elem) = 0.1_WP*modUVice_inv*k2_lf*LOG(exp(10._WP*(b - hc)+1.0_WP))*exp(-c_pressure*(1.0_WP - c))
        Cdu_lf(elem) = Cdu_lf(elem)/rhoice
        Cdv_lf(elem) = Cdv_lf(elem)/rhoice
       else

        Cdu_lf(elem) = 0.0_WP
        Cdv_lf(elem) = 0.0_WP

       endif
   
      else

        Cdu_lf(elem) = 0.0_WP
        Cdv_lf(elem) = 0.0_WP

      endif
  enddo
end subroutine landfast_friction_NS

!--------------------------------------------------------------------
!===================================================================
!
subroutine EVPdynamics
  ! assemble rhs and solve for ice velocity
  ! New implementation based on Bouillion et al. Ocean Modelling 2013
  !---------------------------------------------------------

    use i_PARAM
    use i_ARRAYS
    use i_therm_parms

    use o_MESH
    use o_ARRAYS
    use o_PARAM

    use g_parsup
    use g_comm_auto

    implicit none

  integer         :: shortstep, i, j, elem, elnodes(4), n, step_count!, ice_elem, sum_ice_elem
  real(kind=WP)    :: rdt, drag1, drag2, dragw, detu, detv, fc, a, b, c, inv_mass, inv_mass_a, acc, x1
  real(kind=WP)    :: thickness, inv_thickness, umod, umodw, rhsu, rhsv, wx_el, wy_el, ice_elem, sum_ice_elem

!aa67  FOR TEST
  integer ::    elem_size
  real(kind=WP)    :: aa_u, aa_v, max_aauv, dif_aauv, max_uvice, dif_sum, dif_sum_reol,dif_maxr
    elem_size=myDim_elem2D+eDim_elem2D+eXDim_elem2D
!aa67  FOR TEST

  
    do elem=1,myDim_elem2D
  UiceAB(1,elem) = U_n_ice(1,elem)    ! Initialize solver variables
  UiceAB(2,elem) = U_n_ice(2,elem)
    enddo

!   call ice_strength_mass_p1p1_m   ! Compute arrays that are constant during subcycles
   call gradssh_2_RHice
   call landfast_friction_NS             

     U_rhs_ice_tmp = -999999.0_WP
     dif_sum_reol = 0.0_WP
  step_count = 0
  do shortstep=1, evp_rheol_steps            ! Sub-cycling
     step_count = shortstep
     call stress_tensor
     U_rhs_ice = 0.0_WP
     call stress2rhs
      dif_sum = 0.0_WP

    ice_elem = 0.0_WP
    do elem=1,myDim_elem2D

        elnodes=elem2D_nodes(:,elem)

!aa67  FOR TEST
        aa_u = UiceAB(1,elem)
        aa_v = UiceAB(2,elem)
!aa67    

        a = max( hmin, w_cv(1,elem) *m_ice(elnodes(1)) + w_cv(2,elem) *m_ice(elnodes(2)) &
                     + w_cv(3,elem) *m_ice(elnodes(3)) + w_cv(4,elem) *m_ice(elnodes(4)) )
        b = w_cv(1,elem) *m_snow(elnodes(1)) + w_cv(2,elem) *m_snow(elnodes(2)) &
                     + w_cv(3,elem) *m_snow(elnodes(3)) + w_cv(4,elem) *m_snow(elnodes(4))
        c = w_cv(1,elem) *a_ice(elnodes(1)) + w_cv(2,elem) *a_ice(elnodes(2)) &
                     + w_cv(3,elem) *a_ice(elnodes(3)) + w_cv(4,elem) *a_ice(elnodes(4))

        acc = w_cv(1,elem) *ac(elnodes(1)) + w_cv(2,elem) *ac(elnodes(2)) &
                     + w_cv(3,elem) *ac(elnodes(3)) + w_cv(4,elem) *ac(elnodes(4))

        inv_mass = 1.0_WP/max(0.1_WP,rhoice*a + rhosno*b )
        inv_mass_a = inv_mass*c

        umod = sqrt((U_n(1,elem) - UiceAB(1,elem))**2 + (V_n(1,elem) - UiceAB(2,elem))**2)
        drag1 = dt*(Cd_oce_ice + Cdu_lf(elem))*umod*rhowat*inv_mass_a
        drag2 = dt*(Cd_oce_ice + Cdv_lf(elem))*umod*rhowat*inv_mass_a

        wx_el =   w_cv(1,elem) *windx(elnodes(1)) + w_cv(2,elem) *windx(elnodes(2)) &
                + w_cv(3,elem) *windx(elnodes(3)) + w_cv(4,elem) *windx(elnodes(4))

        wy_el =   w_cv(1,elem) *windy(elnodes(1)) + w_cv(2,elem) *windy(elnodes(2)) &
                + w_cv(3,elem) *windy(elnodes(3)) + w_cv(4,elem) *windy(elnodes(4))

   !     umodw = sqrt( (wx_el - UiceAB(1,elem))**2 + (wy_el - UiceAB(2,elem))**2 )
        umodw = sqrt( wx_el**2 + wy_el**2 )

        dragw = dt*Cda*umodw*rhoair*inv_mass_a

!-----------------
! rhs for water stress, air stress, and rhs_u/v (internal stress + ssh)
!-----------------
         rhsu = U_n_ice(1,elem) + drag1*U_n(1,elem) + dragw*wx_el + dt*U_rhs_ice(1,elem)
         rhsv = U_n_ice(2,elem) + drag2*V_n(1,elem) + dragw*wy_el + dt*U_rhs_ice(2,elem)

         rhsu = beta_evp*UiceAB(1,elem) + rhsu  
	 rhsv = beta_evp*UiceAB(2,elem) + rhsv  
!-----------------
! solve (Coriolis and water stress are treated implicitly)
!-----------------

         fc = dt*coriolis(elem)
         detu = (1.0_WP + beta_evp + drag1)**2 + fc**2
         detv = (1.0_WP + beta_evp + drag2)**2 + fc**2

         detu = mask_wd(elem)/detu
         detv = mask_wd(elem)/detv

         UiceAB(1,elem) = detu*((1.0_WP + beta_evp + drag1 )*rhsu + fc*rhsv)
         UiceAB(2,elem) = detv*((1.0_WP + beta_evp + drag2 )*rhsv - fc*rhsu)


         if (c <= 0.01_WP) then
            UiceAB(1,elem) = 0.0_WP
            UiceAB(2,elem) = 0.0_WP
         else   
            ice_elem = ice_elem+1.0_WP
            dif_aauv = sqrt((UiceAB(1,elem) - aa_u)**2 + (UiceAB(2,elem) - aa_v)**2)
            dif_sum = dif_sum + dif_aauv
            if (dif_aauv >= U_rhs_ice_tmp(elem)) U_rhs_ice_tmp(elem) = dif_aauv
         end if

!aa67
!         if (shortstep .eq. evp_rheol_steps) then
!            dif_aauv = sqrt((UiceAB(1,elem) - aa_u)**2 + (UiceAB(2,elem) - aa_v)**2)
!            dif_sum = dif_sum + dif_aauv
!            if (dif_aauv >= U_rhs_ice_tmp(elem)) U_rhs_ice_tmp(elem) = dif_aauv
!         endif
        
!aa67 

      enddo

#ifdef USE_MPI
     call exchange_elem(UiceAB)
!!!!!!! test
   call MPI_AllREDUCE(dif_sum,dif_maxr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                      MPI_COMM_FESOM_C, MPIerr)
   call MPI_AllREDUCE(ice_elem,sum_ice_elem, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                      MPI_COMM_FESOM_C, MPIerr)
#endif
      !if (dif_maxr/float(elem2D) <= 6.0e-9_WP) exit
      if (dif_maxr/sum_ice_elem <= 6.0e-6_WP) then
              if (step_count>24)  exit
      endif

    end do

     call exchange_elem(U_rhs_ice_tmp)


#ifdef USE_MPI
     call MPI_AllREDUCE(maxval(U_rhs_ice_tmp),max_uvice, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
          MPI_COMM_FESOM_C, MPIerr)
!     call MPI_AllREDUCE(dif_sum,dif_maxr, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
!          MPI_COMM_FESOM_C, MPIerr)
#endif
   !  if (mype==0) write(85,*)  max_uvice, dif_maxr/float(evp_rheol_steps),dif_sum = dif_sum/float(elem_size)
     if (mype==0) write(85,*)  max_uvice, dif_maxr/sum_ice_elem, step_count, sum_ice_elem, time_jd-time_jd0

  
!    do elem=1, myDim_elem2D
!        elnodes=elem2D_nodes(:,elem)
!       c = w_cv(1,elem) *a_ice(elnodes(1)) + w_cv(2,elem) *a_ice(elnodes(2)) &
!                     + w_cv(3,elem) *a_ice(elnodes(3)) + w_cv(4,elem) *a_ice(elnodes(4))

!    if (c <= 0.01_WP) then
!    UiceAB(1,elem) = 0.0_WP
!    UiceAB(2,elem) = 0.0_WP
!    end if
!    end do 


    
    do elem=1,myDim_elem2D
  U_n_ice(1,elem) = UiceAB(1,elem)    ! Initialize solver variables
  U_n_ice(2,elem) = UiceAB(2,elem)
    enddo

     call exchange_elem(U_n_ice)   !!!!!!!!!!!!!!!!!!!!!



!aaaaaaaaaaaaaaaaaaaa
 
!     call compute_AIS_rhs

!     do n=1,myDim_nod2D
     !Calculate AIS for step_n+1
!                    x1 = dt*ac(n)/area(n)
!      m_ice(n) = m_ice(n) + x1*ssh_rhs(n)
!  enddo

!#ifdef USE_MPI
!  call exchange_nod(m_ice)
!#endif
!          do n=1,myDim_nod2D
!                    ! cut off
!                    if (a_ice(n) > 1.0_WP) a_ice(n) = 1.0_WP
!                    if (a_ice(n) < 0.1e-8_WP) a_ice(n) = 0.0_WP
!                    if ( m_ice(n) < 0.1e-8_WP) m_ice(n) = 0.0_WP
!                    if ( m_ice(n) > 2.0_WP) m_ice(n) = 2.0_WP
!                enddo

!#ifdef USE_MPI
!                call exchange_nod(a_ice)
!                call exchange_nod(m_ice)
!                call exchange_nod(m_snow)
!#endif

!aaaaaaaaaaaaaaaaa

                 call ice_con2

                do n=1,myDim_nod2D+eDim_nod2D
                    x1 = dt! /evp_rheol_steps
                    a_ice(n) = a_ice(n) + x1*cHrhsi(1,n)
                    m_ice(n) = m_ice(n) + x1*cHrhsi(2,n)
                    m_snow(n) = m_snow(n) + x1*cHrhsi(3,n)
                end do

                do n=1,myDim_nod2D+eDim_nod2D
                    ! cut off
                    if (a_ice(n) > 1.0_WP) then
                        a_ice(n) = 1.0_WP
                    endif
                    if (a_ice(n) < 0.1e-8_WP) then
                        a_ice(n) = 0.0_WP
                    endif
                    if ( m_ice(n) < 0.1e-8_WP) then
                        m_ice(n) = 0.0_WP
                    endif
                    if ( m_ice(n) > 2.0_WP) then
                        m_ice(n) = 2.0_WP
                    endif
                end do
#ifdef USE_MPI
                call exchange_nod(a_ice)
                call exchange_nod(m_ice)
                call exchange_nod(m_snow)
#endif


end subroutine EVPdynamics
! =============================================================================

!--------------------------------------------------------------------
!
! Modified EVP following Bouillion et al. 2013
! 
! Arrays for pressure and ice mass are introduced.
! Compute auxuliary arrays that are not sub-cycled
! ===================================================================

!subroutine ice_strength_mass_p1p1_m
!  ! Compute auxuliary arrays that are not sub-cycled

!    use i_PARAM
!    use i_ARRAYS
!    use i_therm_parms

!    use o_MESH
!    use o_ARRAYS
!    use o_PARAM

!    use g_parsup
!    use g_comm_auto

!    implicit none

!  integer        :: i, n, m, row

!  do i=1,myDim_nod2D+eDim_nod2D
!      n=myList_nod2D(i)
!      inv_mass(n)=1.0_8/max(0.1,rhoice*m_ice(n)+rhosno*m_snow(n))
!      inv_mass_a(n)=inv_mass(n)*a_ice(n)
!  enddo
  ! Attention: Here we in essence require that there is no dynamics
  ! if sea ice is thinner than 1 mm 
!end subroutine ice_strength_mass_p1p1_m


SUBROUTINE ice_con2

    USE i_ARRAYS
    USE i_PARAM
    USE i_therm_parms
    USE o_MESH
    USE o_ARRAYS
    USE o_PARAM

    use g_parsup
    use g_comm_auto
#ifdef USE_MPI
        USE fv_sbcmpi      ! MPI version of surface boundary cond.
#else
        USE fv_sbc
#endif

IMPLICIT NONE

    integer      :: el(2), enodes(2), n, edge, ed, nodes(2), elem, elnodes(4), j, jend, k, nelem, nfilt, i
    integer      :: edglim, node_size, elem_size

    real(kind=WP) :: deltaX1, deltaY1, deltaX2, deltaY2
    real(kind=WP) :: fD, fDold, x1, y1, dmean, un1
    real(kind=WP) :: a1, a2, b1, b2, c1, c2, u1, u2, Aice, Mice, Msnow, tau_inv, tvol, tx, ty, tz

    real(kind=WP) :: mx_tmp1,mx_tmp2,mx_tmp3

    real(kind=WP) :: tmp_n(myDim_nod2D+eDim_nod2D) !used for MPI node exchange

    node_size=myDim_nod2D+eDim_nod2D
    elem_size=myDim_elem2D+eDim_elem2D+eXDim_elem2D

    ! =========================
    ! equation for
    ! ice concentration (a_ice), ice thickness (m_ice)
    ! and snow (m_snow)
    ! 27.12.2021 AA, IK
    ! =========================

    ! ===============
    ! Clean the rhs
    ! ===============
    cHrhsi = 0.0_WP

    ! ===========================================
    ! Horizontal advection:
    ! Upwind with linear velocity reconstruction
    ! ===========================================
!$OMP DO
#ifdef USE_MPI
    edglim=myDim_edge2D+eDim_edge2D
#else
    edglim=edge2D_in
#endif

! first loop over the all edges inside domain
    do edge=1, edglim

#ifdef USE_MPI
     ! jump out if boundary edge 
     if (mylist_edge2D(edge)>edge2D_in) cycle
#endif

        enodes = edge_nodes(:,edge)
        el = edge_tri(:,edge)

        deltaX1 = edge_cross_dxdy(1,edge)
        deltaY1 = edge_cross_dxdy(2,edge)
        !===============================
        ! First segment: linear upwind reconstruction
        ! ===============================
        if (UiceAB(2,el(1))*deltaX1 - UiceAB(1,el(1))*deltaY1 > 0.0_WP) then
            Aice = a_ice(enodes(2))
            Mice = m_ice(enodes(2))
            Msnow = m_snow(enodes(2))
        else
            Aice = a_ice(enodes(1))
            Mice = m_ice(enodes(1))
            Msnow = m_snow(enodes(1))
        endif
        u1 = UiceAB(2,el(1))*deltaX1
        u2 = UiceAB(1,el(1))*deltaY1
        a1 = u1*Aice - u2*Aice
        b1 = u1*Mice - u2*Mice
        c1 = u1*Msnow - u2*Msnow

        cHrhsi(1,enodes(1)) = cHrhsi(1,enodes(1)) + a1
        cHrhsi(1,enodes(2)) = cHrhsi(1,enodes(2))  - a1
        cHrhsi(2,enodes(1)) = cHrhsi(2,enodes(1)) + b1
        cHrhsi(2,enodes(2)) = cHrhsi(2,enodes(2))  - b1
        cHrhsi(3,enodes(1)) = cHrhsi(3,enodes(1)) + c1
        cHrhsi(3,enodes(2)) = cHrhsi(3,enodes(2))  - c1

        ! ================================
        ! Second segment: linear upwind reconstruction
        !========== ======================
!        if (el(2) > 0) then  !we do not need this check, we are inside the domain

            deltaX2 = edge_cross_dxdy(3,edge)
            deltaY2 = edge_cross_dxdy(4,edge)
            if(UiceAB(2,el(2))*deltaX2 - UiceAB(1,el(2))*deltaY2 < 0.0_WP) then
                Aice = a_ice(enodes(2))
                Mice = m_ice(enodes(2))
                Msnow = m_snow(enodes(2))
            else
                Aice = a_ice(enodes(1))
                Mice = m_ice(enodes(1))
                Msnow = m_snow(enodes(1))
            endif

            u1 = UiceAB(2,el(2))*deltaX2
            u2 = UiceAB(1,el(2))*deltaY2
            a2 = -u1*Aice + u2*Aice
            b2 = -u1*Mice + u2*Mice
            c2 = -u1*Msnow + u2*Msnow

            cHrhsi(1,enodes(1)) = cHrhsi(1,enodes(1)) + a2
            cHrhsi(1,enodes(2))=cHrhsi(1,enodes(2)) - a2
            cHrhsi(2,enodes(1)) = cHrhsi(2,enodes(1)) + b2
            cHrhsi(2,enodes(2))=cHrhsi(2,enodes(2)) - b2
            cHrhsi(3,enodes(1)) = cHrhsi(3,enodes(1)) + c2
            cHrhsi(3,enodes(2))=cHrhsi(3,enodes(2)) - c2

!        end if

    end do

#ifdef USE_MPI
    DO ed=1,edglim
       if (myList_edge2D(ed)>edge2D_in) then !loop over boundary edges only
            enodes = edge_nodes(:,edge)
            el = edge_tri(:,edge)

            deltaX1 = edge_cross_dxdy(1,edge)
            deltaY1 = edge_cross_dxdy(2,edge)
            !===============================
            ! First segment: linear upwind reconstruction
            ! ===============================
            if (UiceAB(2,el(1))*deltaX1 - UiceAB(1,el(1))*deltaY1 > 0.0_WP) then
                Aice = a_ice(enodes(2))
                Mice = m_ice(enodes(2))
                Msnow = m_snow(enodes(2))
            else
                Aice = a_ice(enodes(1))
                Mice = m_ice(enodes(1))
                Msnow = m_snow(enodes(1))
            endif
            u1 = UiceAB(2,el(1))*deltaX1
            u2 = UiceAB(1,el(1))*deltaY1
            a1 = u1*Aice - u2*Aice
            b1 = u1*Mice - u2*Mice
            c1 = u1*Msnow - u2*Msnow

            cHrhsi(1,enodes(1)) = cHrhsi(1,enodes(1)) + a1
            cHrhsi(1,enodes(2)) = cHrhsi(1,enodes(2))  - a1
            cHrhsi(2,enodes(1)) = cHrhsi(2,enodes(1)) + b1
            cHrhsi(2,enodes(2)) = cHrhsi(2,enodes(2))  - b1
            cHrhsi(3,enodes(1)) = cHrhsi(3,enodes(1)) + c1
            cHrhsi(3,enodes(2)) = cHrhsi(3,enodes(2))  - c1
       end if
    END DO
  
#else  
! for serial version loop over boundary edges    
  !$OMP DO
    DO ed=1+edge2D_in, edge2D
        enodes = edge_nodes(:,edge)
        el = edge_tri(:,edge)

        deltaX1 = edge_cross_dxdy(1,edge)
        deltaY1 = edge_cross_dxdy(2,edge)
        !===============================
        ! First segment: linear upwind reconstruction
        ! ===============================
        if (UiceAB(2,el(1))*deltaX1 - UiceAB(1,el(1))*deltaY1 > 0.0_WP) then
            Aice = a_ice(enodes(2))
            Mice = m_ice(enodes(2))
            Msnow = m_snow(enodes(2))
        else
            Aice = a_ice(enodes(1))
            Mice = m_ice(enodes(1))
            Msnow = m_snow(enodes(1))
        endif
        u1 = UiceAB(2,el(1))*deltaX1
        u2 = UiceAB(1,el(1))*deltaY1
        a1 = u1*Aice - u2*Aice
        b1 = u1*Mice - u2*Mice
        c1 = u1*Msnow - u2*Msnow

        cHrhsi(1,enodes(1)) = cHrhsi(1,enodes(1)) + a1
        cHrhsi(1,enodes(2)) = cHrhsi(1,enodes(2))  - a1
        cHrhsi(2,enodes(1)) = cHrhsi(2,enodes(1)) + b1
        cHrhsi(2,enodes(2)) = cHrhsi(2,enodes(2))  - b1
        cHrhsi(3,enodes(1)) = cHrhsi(3,enodes(1)) + c1
        cHrhsi(3,enodes(2)) = cHrhsi(3,enodes(2))  - c1        
    END DO
  !$OMP END DO
#endif
    
!$OMP END DO
!$OMP DO
    do n=1,myDim_nod2D
        cHrhsi(1,n) = ac(n)*cHrhsi(1,n)/area(n)
        cHrhsi(2,n) = ac(n)*cHrhsi(2,n)/area(n)
        cHrhsi(3,n) = ac(n)*cHrhsi(3,n)/area(n)
    end do
!$OMP END DO
#ifdef USE_MPI
    do i=1,3
        tmp_n = cHrhsi(i,:)
        call exchange_nod(tmp_n)
        cHrhsi(i,:) = tmp_n
    enddo
#endif
end subroutine ice_con2

