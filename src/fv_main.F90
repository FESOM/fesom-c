! Shallow water code based cell-vertex finite volume discretization
! AB3-AM4 time stepping as suggested in Shchepetkin and McWilliams
! (2005)
! Serial version
! July 2012
! sergey.danilov@awi.de

PROGRAM MAIN

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM
  USE o_UTILIT
  USE i_ARRAYS
  USE i_PARAM
!  USE turbulence, only: init_turbulence    no gotm
!  USE mtridiagonal, only: init_tridiagonal no gotm
#ifdef USE_MPI
  USE pnetcdf ! lib. PnetCDF (read/write netcdf in MPI)
  USE fv_ncoutputmpi ! MPI version of netCDF output
  USE fv_sbcmpi      ! MPI version of surface boundary cond.
#else
  USE fv_ncoutput
  USE fv_sbc
#endif

  USE fv_obc
  USE fv_rivers

  USE g_parsup
  use g_comm_auto

  IMPLICIT NONE

  integer      :: i, k, n, period_m2, out_fft, istep, ist,elnodes(4), nz, rk, rk2, turn_on_riv, n_dt2, flag_riv,vsp(11), nod
  integer      :: fid_ssh_fft=111, fid_vel2D_fft=100, fid_vel3Du_fft=3, fid_vel3Dv_fft=4, Tfile=5, Sfile=6, tracer_sch, elem, nn, sk
  real(kind=WP) :: eout, x, y, tt, ss, pp, pr
  real(kind=WP), external :: theta
  double precision :: L_min
  integer      :: ierr ! return values, error flags

  real(kind=WP) :: mx_t, mn_t, mx_s, mn_s, mx_eta, mn_eta
  real(kind=WP) :: mx_mice,mn_mice,mx_aice,mn_aice,mx_uice,mn_uice,mx_vice,mn_vice
  real(kind=WP) :: mx_bp2, mn_bp2, mx_bpu3, mx_bpv3, mx_un2, mn_un2, mx_un, mx_vn, mx_w, mn_w
  real(kind=WP) :: mx_cd, mn_cd, mx_av, mx_kv

   real(kind=WP) :: mx_cdulf, mx_cdvlf

  real(kind=WP) :: start_time, end_time, t0, t_end, t1,t2,t3,t4
  integer :: c1,c1rate,c1max,c2,c2rate,c2max

  integer :: IREP_NC, ierror, fID, ind

  logical :: enable_output_main_switch
  real(kind=WP)   :: t_0, t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8



  real(kind=WP) ::   tau_inv_filt_0 ! coefficient for reducing clima_relax after initialization (several days))
  logical :: SHnc !netcdf output by SH

  character(len = 256) :: fname,tmpstr ! file name

  SHnc = .FALSE.
  enable_output_main_switch = .TRUE.

#ifdef USE_MPI
  call MPI_INIT(ierr)
  call par_init
#else
  print *,'Run serial version'
  mype=0
  npes=1
#endif


!  print *,'Init phase: mype,npes = ',mype,npes


  !++++++++++++++++++++++++++
  ! Open some output files
  !++++++++++++++++++++++++++

  if (mype==0) then
    open(51,file='max_min_ssh.dat')
    !open(52,file='max_min_vel.dat')
    open(53,file='energy.dat')
    open(54,file='cpu_time.dat')
    open(55,file='MPI_Wtime.dat')
  end if


  !++++++++++++++++++++++++++
  ! Set model parameters
  !++++++++++++++++++++++++++

  call READ_DATA_RUN

  IREP_NC=IREP

  !++++++++++++++++++++++++++
  ! Load the mesh
  !++++++++++++++++++++++++++

#ifdef USE_MPI
  ! taken from subroutine mesh_setup
  call read_mesh_par
  call set_par_support
#else

  call read_mesh_ser

  ! These values are set in read_mesh_ser
  ! myDim_elem2D = elem2D
  ! eDim_elem2D  = 0
  ! eXDim_elem2D = 0
  ! myDim_nod2D  = nod2D
  ! seDim_nod2D   = 0
#endif

  !++++++++++++++++++++++++++
  ! Assemble mesh arrays
  !++++++++++++++++++++++++++

  call test_elem

#ifdef USE_MPI
  call load_edges

#else
  call find_edges
  myDim_edge2D = edge2D
  eDim_edge2D  = 0
#endif

!  if (.not.cartesian) then
!    print *,'PE, Longitude range: ',mype,minval(coord_nod2D(1,:))/rad,maxval(coord_nod2D(1,:))/rad
!    print *,'PE, Latitude range : ',mype,minval(coord_nod2D(2,:))/rad,maxval(coord_nod2D(2,:))/rad
!    print *,'PE, Depth range    : ',mype,minval(depth(:)),maxval(depth(:))
!  end if


  call find_elem_neighbors


  call mesh_arrays1
  call mesh_arrays2


!!$#ifdef USE_MPI
!!$  call par_ex
!!$  stop 'MPI Version: Leaving for now...'
!!$#else
!!$  stop 'Serial Version: Leaving for now...'
!!$#endif

  call array_setup

!IK   call find_up_downwind_triangles
!IK ! edge_up_dn_tri is not calculated in fv_mesh_array, due to the error in the coord_elem(:,:,:) , zeors (0) in coordinates of elements beggire then my_elem2d

  if (type_task>1) call jacobian

  !=====================
  ! Initialize fields
  !=====================

  call initial_state

  if (type_task>2) then

#ifdef USE_MPI
     call MPI_AllREDUCE(maxval(TF),mx_t, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
          MPI_COMM_FESOM_C, MPIerr)
     call MPI_AllREDUCE(minval(TF),mn_t, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
          MPI_COMM_FESOM_C, MPIerr)
     call MPI_AllREDUCE(maxval(SF),mx_s, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
          MPI_COMM_FESOM_C, MPIerr)
     call MPI_AllREDUCE(minval(SF),mn_s, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
          MPI_COMM_FESOM_C, MPIerr)

     T_aver=0.5_WP*(mn_t+mx_t)
     S_aver=0.5_WP*(mn_s+mx_s)
     T_maxini = mx_t
#else
     T_aver=0.5_WP*(minval(TF)+maxval(TF))
     S_aver=0.5_WP*(minval(SF)+maxval(SF))
     T_maxini = maxval(TF)
#endif

     if (mype==0) write(*,*) 'aver TS= ', T_aver, S_aver

  end if

  !=====================
  ! Initialize turbulence related stuff
  !=====================

!IK print *,mype,type_task,ver_mix

  if ((type_task>1).and.(ver_mix == 2)) then
!SH SKIPPED FOR NOW, NOT YET USED     call init_turbulence(namlst,'gotmturb.nml',nsigma-1)
!SH SKIPPED FOR NOW, NOT YET USED     call init_tridiagonal(nsigma-1)
     if (mype==0) print *,'init_turbulence completed!'
     STOP 'GOTM turned off'
  end if


  !=====================
  ! Time stepping
  !=====================


!  write(*,*) ini_time

  !SH open preliminary netcdf output file
  !SH TO BE REPLACED
  allocate(coord_nod2D_glob(2,nod2D), index_nod2D_glob(nod2D))
  allocate(elem2D_nodes_glob(4,elem2D))
  if (SHnc) then
    allocate(eta_n_2_glob(nod2D))
    allocate(TF_glob(nsigma-1,nod2D))
  endif
  !Generate ac file for sponge layer
  ALLOCATE(ac(myDim_nod2D+eDim_nod2D))
  ac(:)=1.0_WP

#ifdef USE_MPI

  call gather_nod(coord_nod2D, coord_nod2D_glob)
  call gather_nod(index_nod2D, index_nod2D_glob)


  if (nobn>0) then
    allocate(X1obn(nobn),X2obn(nobn),in_obn(nobn))

    if ((SL_obn).and.(.not.(SL_obn_user_def))) then

        if (mype==0) call ac_create_Xobn

        !call MPI_BARRIER(MPI_COMM_FESOM_C,ierror)
        call MPI_BCast(X1obn,nobn,MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
        call MPI_BCast(X2obn,nobn,MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
        call MPI_BCast(in_obn,nobn,MPI_INTEGER, 0, MPI_COMM_FESOM_C, ierror)
        call ac_create(ac)
        call exchange_nod(ac)
        if (mype==0) print *,'AC values (open boundaries) are set.'
    end if

    if (mynobn>0) then
     ind=0
     do n=1,myDim_nod2D
        if (index_nod2D(n)==2) then
           ind=ind+1
           do k=1,nobn
              if (myList_nod2D(n)==in_obn(k)) then
                 my_in_obn_idx(ind)=k
                 exit
              end if
           end do
           if (my_in_obn_idx(ind)==0) print *,'WARNING: open bnd index discrepancy'
        endif

     end do
    end if
  endif
!if (mynobn>0) then
!        write(tmpstr,*) mype
!        write(fname,*) 'ampt_',trim(ADJUSTL(trim(tmpstr))),'.out'
!       open(50,file=fname)
!       write(50,*) "mype=",mype,"mynobn=",mynobn
!       write(50,*) "my_in_obn="
!       write(50,*) my_in_obn
!       write(50,*) "my_in_obn_idx="
!       write(50,*) my_in_obn_idx
!       write(50,*) "myList_nod2D(my_in_obn_idx)"
!       do n=1,mynobn
!       write(50,*) myList_nod2D(my_in_obn(n))
!       end do
!endif


#else

  coord_nod2D_glob=coord_nod2D
  elem2D_nodes_glob=elem2D_nodes
  index_nod2D_glob=index_nod2D

  !aaa  if (SL_obn_user_def) then
  !aaa  open (23,file=trim(meshpath)//trim(TITLE)//'_ac.out', status='old')
  !aaa  read(23,*) ac
  !aaa  CLOSE(23)
  !aaa  endif

  if (nobn>0) then
    allocate(X1obn(nobn),X2obn(nobn))
    if ((SL_obn).and.(.not.(SL_obn_user_def))) then
        call ac_create_Xobn
        call ac_create(ac)
        if (mype==0) print *,'AC values (open boundaries) are set.'
    end if
  endif
open(22,file='ac_values_ser')
do n=1,nod2D
   write (22,*) n,ac(n)
end do
close(22)

#endif
  if (SHnc) then
    if (mype==0) call nc_init
  endif

  deallocate(coord_nod2D_glob,elem2D_nodes_glob)

  !++++++++++++++++++++++++++++++++++
  ! compute mask for wetting/drying
  !++++++++++++++++++++++++++++++++++

  if (WET_DRY_ON) then
!$OMP PARALLEL
     call wad_mask
!$OMP END PARALLEL
  endif

  time_jd = time_jd0


  !=================================
  ! inizialization of ice
  !=================================
  if (use_ice) call ice_setup

  !Read restart (after allocation in initialization of ice ...)
  ! sbc,obc,rivers should be after restart
  if(restart) call read_restart_separate
  !=================================
  ! inizialization of ocean forcing
  !=================================
#ifdef USE_MPI
  if (key_atm) call sbcmpi_ini
#else
  if (key_atm) call sbc_ini
#endif
  !=================================
  ! initialization of open boundary, fv_obc.F90, work for serial and MPI
  !=================================
  if (key_obc)    call obc_ini
  !=================================
  ! initialization of rivers, fv_rivers.F90, work for serial and MPI
  !=================================
  if (key_rivers)    call rivers_ini
  !=================================
  ! Initialization of NC output
  !=================================
#ifdef USE_MPI
  if (key_nc_output) call ncoutputmpi_ini
#else
  if (key_nc_output) call output_ini
#endif

  if (type_task>1) call jacobian

  if (comp_sediment) then
     h_var_old = h_var
     h_var_old2 = h_var2
  end if

  tau_inv_filt_0 = tau_inv_filt !initial state
!===============================
!
!        MAIN LOOP
!
!===============================

  n_dt2=0
  sk=1;


!#ifdef USE_MPI
!  t0=MPI_Wtime()
!#else
!  call cpu_time(start_time)
!  call system_clock(c1,c1rate,c1max)
!  print *,'CLOCK:',c1,c1rate,c1max
!#endif

  do n_dt = 1, nsteps

!if (mype==0) then
!print *,'============'
!print *,'STEP:',n_dt
!endif

!     if (mype==0) write(*,*) "time step: ", n_dt
     n_dt2=n_dt2+1

     !$ if (iverbosity >= 2) t1=omp_get_wtime()

     time_jd = time_jd+dt/86400.0_WP
     time = time_jd*86400.0_WP

     ! set date and seconds for new time step
     !     run_date  = time_jd
     !     run_sec   = mod(nint(time),86400)

     ! surface boundary conditions
     !$ if (iverbosity >= 3) t2=omp_get_wtime()
!#ifdef USE_MPI
! t_0=MPI_Wtime()
!#endif
     if (key_rivers) call rivers_do
!#ifdef USE_MPI
! t_1=MPI_Wtime()
!#endif
#ifdef USE_MPI
     if (key_atm) call sbcmpi_do
#else
     if (key_atm) call sbc_do
#endif
!#ifdef USE_MPI
! t_2=MPI_Wtime()
!#endif
     if (use_ice) call ice_timestep
!#ifdef USE_MPI
! t_3=MPI_Wtime()
!#endif
     !$ if (iverbosity >= 3) t3=omp_get_wtime()

     ! call change in Tclim and Sclim before time step, works for MPI and serial
     if (key_obc) call obc_do
!#ifdef USE_MPI
! t_4=MPI_Wtime()
!#endif
!#ifdef USE_MPI
!  t1=MPI_Wtime()
!#endif

     call oce_timestep

!#ifdef USE_MPI
!  t2=MPI_Wtime()
!#endif
!#ifdef USE_MPI
! t_5=MPI_Wtime()
!#endif
!#ifdef USE_MPI
!  #if (mype==0) print *,'ZWISCHENZEIT: oce_timestep: ',t2-t1
!  #if (mype==0) write(54,*) t2-t1

!#endif



     !*****************************************
     ! part of sediment model (every tidal cycle change the bottom due to sediment)
     ! AA
!     If (comp_sediment) then
!
!        if ( ((n/period_m2)*period_m2) == n) then
!           do nn=1,myDim_nod2D
!              depth(nn) = depth(nn) + (h_var(nn) - h_var_old(nn))
!              h_var_old(nn) = h_var(nn)
!              !VF: alternative calculation
!              h_var_old2(nn) = h_var2(nn)
!           enddo
!        endif
!
!     endif
     !AA
     !*************************************

     !make NetCDF output
#ifdef USE_MPI
     if (key_nc_output) call ncoutputmpi_do
#else
     if (key_nc_output) call output_do
#endif
!#ifdef USE_MPI
! t_6=MPI_Wtime()
!#endif
     !$ if (iverbosity >= 3) t4=omp_get_wtime()

!#ifdef USE_MPI
!
!     call MPI_REDUCE(maxval(eta_n),mx_eta, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
!          0, MPI_COMM_FESOM_C, MPIerr)
!     call MPI_REDUCE(minval(eta_n),mn_eta, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
!          0, MPI_COMM_FESOM_C, MPIerr)
!
!     if (mype==0) write(51,'(3e13.5)') (time- time_jd0*86400.0_WP)/3600.0_WP,mx_eta,mn_eta
!!!     if (mype==0) write(*,'(3e13.5)') (time- time_jd0*86400.0_WP)/3600.0_WP,mx_eta,mn_eta
!
!#else
!
!     write(51,'(3e13.5)') (time- time_jd0*86400.0_WP)/3600.0_WP,maxval(eta_n),minval(eta_n)
!!     write(*,'(3e13.5)') (time- time_jd0*86400.0_WP)/3600.0_WP,maxval(eta_n),minval(eta_n)
!
!#endif

     !==================================================
     !VF, update sediments concentration at the ob, if necessary
     !==================================================

!SH skipping for now     if (sed_boun_flux)    call update_info_sed(sk)

     !=======================
     ! ouput only for control
     !=======================


     if ( enable_output_main_switch .AND. mod(n_dt,IREP)==0 ) then

#ifdef USE_MPI
        call MPI_AllREDUCE(maxval(eta_n),mx_eta, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
             MPI_COMM_FESOM_C, MPIerr)
        call MPI_AllREDUCE(minval(eta_n),mn_eta, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
             MPI_COMM_FESOM_C, MPIerr)
        if (type_task>1) then
           call MPI_AllREDUCE(maxval(U_n),mx_un, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                MPI_COMM_FESOM_C, MPIerr)
           call MPI_AllREDUCE(maxval(V_n),mx_vn, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                MPI_COMM_FESOM_C, MPIerr)
        end if
        call MPI_AllREDUCE(maxval(u_n_2D),mx_un2, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
             MPI_COMM_FESOM_C, MPIerr)
        call MPI_AllREDUCE(minval(u_n_2D),mn_un2, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
             MPI_COMM_FESOM_C, MPIerr)
        if (type_task>1) then
           call MPI_AllREDUCE(maxval(Wvel(:,1)),mx_w, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                MPI_COMM_FESOM_C, MPIerr)
           call MPI_AllREDUCE(minval(Wvel(:,1)),mn_w, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
                MPI_COMM_FESOM_C, MPIerr)
        end if
        call MPI_AllREDUCE(maxval(Bar_pr_2D),mx_bp2, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
             MPI_COMM_FESOM_C, MPIerr)
        call MPI_AllREDUCE(minval(Bar_pr_2D),mn_bp2, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
             MPI_COMM_FESOM_C, MPIerr)
        if (type_task>1) then
           call MPI_AllREDUCE(maxval(Bar_pru_3D),mx_bpu3, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                MPI_COMM_FESOM_C, MPIerr)
           call MPI_AllREDUCE(maxval(Bar_prv_3D),mx_bpv3, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                MPI_COMM_FESOM_C, MPIerr)
        end if
        call MPI_AllREDUCE(maxval(C_d_el),mx_cd, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
             MPI_COMM_FESOM_C, MPIerr)
        call MPI_AllREDUCE(minval(C_d_el),mn_cd, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
             MPI_COMM_FESOM_C, MPIerr)
        if (type_task>1) then
           call MPI_AllREDUCE(maxval(Av),mx_av, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                MPI_COMM_FESOM_C, MPIerr)
           call MPI_AllREDUCE(maxval(Kv),mx_kv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                MPI_COMM_FESOM_C, MPIerr)
        end if
        if (allocated(TF)) then
           call MPI_AllREDUCE(maxval(TF),mx_t, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                MPI_COMM_FESOM_C, MPIerr)
           call MPI_AllREDUCE(minval(TF),mn_t, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
                MPI_COMM_FESOM_C, MPIerr)
           call MPI_AllREDUCE(maxval(SF),mx_s, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                MPI_COMM_FESOM_C, MPIerr)
           call MPI_AllREDUCE(minval(SF),mn_s, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
                MPI_COMM_FESOM_C, MPIerr)
        end if
        if (use_ice) then
           call MPI_AllREDUCE(maxval(m_ice),mx_mice, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                MPI_COMM_FESOM_C, MPIerr)
           call MPI_AllREDUCE(minval(m_ice),mn_mice, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
                MPI_COMM_FESOM_C, MPIerr)
           call MPI_AllREDUCE(maxval(a_ice),mx_aice, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                MPI_COMM_FESOM_C, MPIerr)
           call MPI_AllREDUCE(minval(a_ice),mn_aice, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
                MPI_COMM_FESOM_C, MPIerr)
           call MPI_AllREDUCE(minval(U_n_ice(1,:)),mx_uice, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                MPI_COMM_FESOM_C, MPIerr)
           call MPI_AllREDUCE(minval(U_n_ice(1,:)),mn_uice, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
                MPI_COMM_FESOM_C, MPIerr)
           call MPI_AllREDUCE(minval(U_n_ice(2,:)),mx_vice, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                MPI_COMM_FESOM_C, MPIerr)
           call MPI_AllREDUCE(minval(U_n_ice(2,:)),mn_vice, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
                MPI_COMM_FESOM_C, MPIerr)
          call MPI_AllREDUCE(maxval(Cdu_lf),mx_cdulf, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                MPI_COMM_FESOM_C, MPIerr)
          call MPI_AllREDUCE(maxval(Cdv_lf),mx_cdvlf, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                MPI_COMM_FESOM_C, MPIerr)


         endif
#else
        mx_mice = maxval(m_ice)
        mn_mice = minval(m_ice)
        mx_aice = maxval(a_ice)
        mn_aice = minval(a_ice)
        mx_uice = maxval(U_n_ice(1,:))
        mn_uice = minval(U_n_ice(1,:))
        mx_vice = maxval(U_n_ice(2,:))
        mn_vice = minval(U_n_ice(2,:))
#endif

        if (SHnc) then
        if ( mod(n_dt,IREP_NC)==0 ) then
           print *,'OUTPUT',mype
#ifdef USE_MPI
           call gather_nod(eta_n_2, eta_n_2_glob)
           if (allocated(TF)) call gather_nod(TF, TF_glob)
#else
           eta_n_2_glob=eta_n_2
           if (allocated(TF)) TF_glob=TF
#endif
           if (mype==0) call nc_out(n_dt)
        end if
        endif

        call energy(eout)

        if (mype==0) then
           print *,'Energy ENERGY ENERGY : ',eout
           write(*,*) "icedyn_tmask = ",icedyn_tmask
           write(53,'(4e16.7)') time/86400.0_WP-time_jd0,eout
           if (iverbosity >= 1) then
              write(*,*) 'time= ',time/86400.0_WP-time_jd0,'time_all= ',time/86400.0_WP-time_jd0
              write(*,*)  ' energy= ',eout
              write(*,*) 'Mass Conserv. imbalance: max, index_max, mean', ib, im_index, ib_mean


              if (type_task>1) then
#ifdef USE_MPI
                 write(*,*) 'max_min_ssh:' , mx_eta, mn_eta
                 write(*,*) 'max_min_vel 2D:', mx_un2, mn_un2

                 write(*,*) 'max_vel 3D:', mx_un, mx_vn
                 write(*,*) 'vert_vel= ', mx_w, mn_w
                 write(*,*) 'maxmin_2Dpr=',mx_bp2, mn_bp2
                 write(*,*) 'max_3Dpr=',mx_bpu3, mx_bpv3
                 write(*,*) 'C_d_el=',mx_cd, mn_cd
                 write(*,*) 'max_vert_vis_dif= ', mx_av, mx_kv
                 if (type_task>2) then
                    write(*,*) 'T_maxmin= ', mx_t, mn_t
                    write(*,*) 'S_maxmin= ', mx_s, mn_s
                    write(*,*) 'dt', dt
                 endif
#else
                 write(*,*) 'max_min_ssh:' , maxval(eta_n), minval(eta_n)
                 write(*,*) 'max_min_vel 2D:', maxval(U_n_2D),minval(U_n_2D)

                 if (type_task>2) then
                    write(*,*) 'T_maxmin= ', maxval(TF),minval(TF)
                    write(*,*) 'S_maxmin= ', maxval(SF),minval(SF)
                    write(*,*) 'dt', dt
                 endif
#endif
                 if (use_ice) then
                    write(*,*) 'mice_maxmin= ', mx_mice, mn_mice
                    write(*,*) 'aice_maxmin= ', mx_aice, mn_aice
                    write(*,*) 'uice_maxmin= ', mx_uice, mn_uice
                    write(*,*) 'vice_maxmin= ', mx_vice, mn_vice
                    write(*,*) 'CD_LANDFAST_ICE=', mx_cdulf, mx_cdvlf
                 end if
              endif
              write(*,*) 'dt_2D', dt_2D
           endif
        endif
     end if
!#ifdef USE_MPI
! t_7=MPI_Wtime()
!#endif
     !aaa output for LE only!!!!
!SH skipped for now     if ( mod(n_dt, IREC)==0 ) then
!SH skipped for now        if (type_task>1) call cross_sec_LE
!SH skipped for now        write(*,*) 'output cross section LE experiment', time/60.0_WP
!SH skipped for now     endif
     !aaa
     !if ( mod(n_dt, IREC)==0 ) then
     !if (type_task>1) then
     !  call compute_vortex_2D
     !  call compute_vortex_3D
     !  call output_vert_vel
     ! vsp =(/81, 182, 283,384, 485, 586,687,788,889,990,1091/)

     !  If (comp_sediment) write(fid_ssh_fft,'(e14.6)') CF(:,vsp)
     !   write(fid_vel2D_fft,'(e14.6)') Av_node(:,vsp), Kv(:,vsp), tke(:,vsp), teps(:,vsp)

     ! AA
     ! output sediment
     !         open(3,file=cross_file1)
     !         do n=1,nod2D
     !            x = coord_nod2d(1,n)*r_earth
     !           y = coord_nod2d(2,n)*r_earth
     !          write(3,'(2f12.4,2e16.8)') x,y,con(n),h_var(n)
     !	enddo
     !AA
     !          endif

     !        call cross_sec
     !        write(*,*) 'output cross section TF'
     !  endif
     if ( mod(n_dt, IRESTART)==0 ) call write_restart_separate(n_dt)
!#ifdef USE_MPI
! t_8=MPI_Wtime()
!#endif
!#ifdef USE_MPI
!    t_8=t_8-t_7
!    t_7=t_7-t_6
!    t_6=t_6-t_5
!    t_5=t_5-t_4
!    t_4=t_4-t_3
!    t_3=t_3-t_2
!    t_2=t_2-t_1
!    t_1=t_1-t_0
!    call MPI_AllREDUCE(t_1,t1, 1, MPI_DOUBLE_PRECISION, MPI_MAX,MPI_COMM_FESOM_C, MPIerr)
!    t_1=t1
!    call MPI_AllREDUCE(t_2,t1, 1, MPI_DOUBLE_PRECISION, MPI_MAX,MPI_COMM_FESOM_C, MPIerr)
!    t_2=t1
!    call MPI_AllREDUCE(t_3,t1, 1, MPI_DOUBLE_PRECISION, MPI_MAX,MPI_COMM_FESOM_C, MPIerr)
!    t_3=t1
!    call MPI_AllREDUCE(t_4,t1, 1, MPI_DOUBLE_PRECISION, MPI_MAX,MPI_COMM_FESOM_C, MPIerr)
!    t_4=t1
!    call MPI_AllREDUCE(t_5,t1, 1, MPI_DOUBLE_PRECISION, MPI_MAX,MPI_COMM_FESOM_C, MPIerr)
!    t_5=t1
!    call MPI_AllREDUCE(t_6,t1, 1, MPI_DOUBLE_PRECISION, MPI_MAX,MPI_COMM_FESOM_C, MPIerr)
!    t_6=t1
!    call MPI_AllREDUCE(t_7,t1, 1, MPI_DOUBLE_PRECISION, MPI_MAX,MPI_COMM_FESOM_C, MPIerr)
!    t_7=t1
!    call MPI_AllREDUCE(t_8,t1, 1, MPI_DOUBLE_PRECISION, MPI_MAX,MPI_COMM_FESOM_C, MPIerr)
!    t_8=t1
!    if (mype==0) then
!        write(55,'(8e16.7)') t_1,t_2,t_3,t_4,t_5,t_6,t_7,t_8
!    endif
!#endif

     !$ if (iverbosity >= 2) then
     !$    t5=omp_get_wtime()
     !$    if (iverbosity >= 3) then
     !$       write(*,'("-- sbc_do took           ",f10.4," s --")') t3-t2
     !$       write(*,'("-- oce_timestep took     ",f10.4," s --")') t4-t3
     !$       write(*,'("-- output+diagnosis took ",f10.4," s --")') t5-t4
     !$    endif
     !$       write(*, '(" =========== STEP took",f10.4," s =========== " )')  t5-t1
     !$       write(*,*)  " "
     !$ end if

  end do


!===============================
!
!    END OF MAIN LOOP
!
!===============================


#ifdef USE_MPI
  t_end=MPI_Wtime()
  if (mype==0) print *,'Time needed for main loop (parallel) p0: ',t_end-t0
  if (mype==1) print *,'Time needed for main loop (parallel) p1: ',t_end-t0
#else
  call cpu_time(end_time)
  print *,'Time needed for main loop (serial):',end_time-start_time
  call system_clock(c2,c2rate,c2max)
  print *,'CLOCK:',c2,c2rate,c2max
  print *,'CLOCKDIFF:',real(c2-c1,4)/real(c2rate,4)
#endif

  close(fid_ssh_fft)
  close(fid_vel2D_fft)

  if (mype==0) close(51)
  !close(52)
  if (mype==0) close(53,status='KEEP')  ! close file ---> energy.dat
  if (mype==0) close(54,status='KEEP')  ! close file ---> cpu

#ifdef USE_MPI
  call par_ex
  stop 'MPI Version: Leaving for now...'
#else
  stop 'Serial Version: Leaving for now...'
#endif


END PROGRAM MAIN



!===========================================================================
! Read mesh in serial version

subroutine read_mesh_ser

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM
  USE o_UTILIT

  use g_parsup

  IMPLICIT NONE

  INTEGER               :: nq, nt, n_quad, n_tri
  INTEGER               :: n1,n2,n3, n4, nod(4)
  INTEGER               :: n,ind
  REAL(kind=WP)         :: x1, x2
  INTEGER, allocatable  :: elem_data(:)
  INTEGER               :: i_error

  ! Requires new mesh format: elem2D lists 4 nodes in each row.
  ! For triangles the forth node coincides with the first one or just list the numbers of 3 nodes

  open (20,file=trim(meshpath)//'nod2d.out', status='old')
  open (21,file=trim(meshpath)//'elem2d.out', status='old')

  READ(20,*) nod2D
  ALLOCATE(coord_nod2D(2,nod2D),index_nod2D(nod2D))
  nobn=0

  if (cartesian) then
     do n=1,nod2D
        read(20,*) nq, x1, x2, ind
        index_nod2D(nq) = ind
        coord_nod2D(1,nq)=x1/r_earth
        coord_nod2D(2,nq)=x2/r_earth
        if (ind==2) then
           nobn=nobn+1
        endif
     end do
  else
     do n=1,nod2D
        read(20,*) nq, x1, x2, ind
        index_nod2D(nq) = ind
        coord_nod2D(1,nq)=x1*rad
        coord_nod2D(2,nq)=x2*rad
        if (ind==2) then
           nobn=nobn+1
        endif
     end do
  endif

  close(20)

  mynobn=nobn ! Defined for compatibility with MPI version

  allocate(in_obn(nobn), my_in_obn(nobn), my_in_obn_idx(nobn))
  ind=0
  do n=1,nod2D
     if (index_nod2D(n)==2) then
        ind=ind+1
        in_obn(ind)=n
        my_in_obn_idx(ind)=ind
     endif
  end do

  my_in_obn=in_obn ! Defined for compatibility with MPI version

  read(21,*)  elem2D
  ALLOCATE(elem2D_nodes(4,elem2D))
  ALLOCATE(elem_data(4*elem2D))
  elem_data(:)=-1

  ! meshes with quads have 4 columns, but TsunAWI grids may be
  ! purely triangular, with 3 columns each. Test, how many
  ! columns there are!

  read(21,*,iostat=i_error) elem_data(1:4*elem2D)
  write(*,*) i_error
  if (i_error == 0) then
     ! There is a fourth column => quad or mixed mesh
     n_quad = 0
     n_tri = 0
     do n=1,elem2D
        nod(1:4) = elem_data((n-1)*4+1:n*4)
        ! It might be important to have the elements sorted:
        ! triangles first.
        if (nod(1) == nod(4)) then   ! triangle
           n_tri = n_tri+1
           elem2D_nodes(1:4,n_tri) = nod(1:4)
        endif
     enddo
     do n=1,elem2D
        nod(1:4) = elem_data((n-1)*4+1:n*4)
        if (nod(1) /= nod(4)) then   ! quad
           n_quad = n_quad+1
           elem2D_nodes(1:4,n_tri + n_quad) = nod(1:4)
        endif
     enddo

  else
     ! No third column => triangles only
     n_quad=0
     n_tri=elem2D
     do n=1,elem2D
        elem2D_nodes(1:3,n) = elem_data((n-1)*3+1:n*3)
        elem2D_nodes(4,n) = elem_data((n-1)*3+1)
     enddo
  end if

  close(21)

  deallocate(elem_data)

  elem2D_tri = n_tri

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! ATTENTION: The list of elements may now differ from that in elem2d.out.
  !  In order to work with the same lists put triangles first.
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  myDim_elem2D = elem2D
  eDim_elem2D  = 0
  eXDim_elem2D = 0
  myDim_nod2D  = nod2D
  eDim_nod2D   = 0


  allocate(depth(nod2D))

  if (len(trim(title))==2 .AND. title=='LE') then
     depth=20.0_WP
  else
     open (22,file=trim(meshpath)//'depth.out', status='old')
     read(22,*) depth
     CLOSE(22)
  end if

!SHTEST TOPOGRAPHY
!do n=1,nod2D
!  if (depth(n)<10.0) depth(n)=10.0
!end do

  !SH check this!! (what if part of domain is dry)
  !do n=1,nod2D
  !   if(depth(n) < Dcr ) depth(n) = Dcr  ! depth is positive
  !enddo

  if (type_task>1) then
     call SET_SIGMA
  endif

  write(*,*) 'Mesh is read     ', 'nod2D=', nod2D,' elem2D=', elem2D
  write(*,*) 'Mesh includes    ', elem2D_tri, 'triangles'
  write(*,*) 'Mesh includes    ', elem2D-elem2D_tri, 'quads'
  write(*,*) 'Amount of open boundary nodes ', nobn
  write(*,*) 'The list of elements starts from triangles !!!'

  mynobn=nobn

END SUBROUTINE read_mesh_ser

SUBROUTINE array_setup

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  USE i_PARAM
  USE i_ARRAYS

  use g_parsup

  IMPLICIT NONE

  integer            :: sbc_alloc
  integer            :: node_size, elem_size

  node_size=myDim_nod2D+eDim_nod2D
  elem_size=myDim_elem2D+eDim_elem2D+eXDim_elem2D

  if (use_wef)  then
    allocate(wef(node_size), STAT=sbc_alloc )
    if( sbc_alloc /= 0 )   STOP 'array_setup: failed to allocate arrays'
    wef=0.0_WP
    if (type_task>2) then
        write(*,*) "!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!"
        write(*,*) "!!!! type_task > 2 and use_wef = True"
    end if
  end if

  if (key_atm .or. use_ice) then
    allocate(Ch_atm_oce_arr(node_size),Ce_atm_oce_arr(node_size),&
            & Cd_atm_oce_arr(node_size),Cd_atm_ice_arr(node_size), net_heat_flux(node_size),&
            & evaporation(node_size), fresh_water_flux(node_size), real_salt_flux(node_size), &
            & STAT=sbc_alloc )
    if( sbc_alloc /= 0 )   STOP 'array_setup: failed to allocate arrays'
    Ch_atm_oce_arr = 0.0_WP
    Ce_atm_oce_arr = 0.0_WP
    Cd_atm_oce_arr = 0.0_WP
    Cd_atm_ice_arr = 0.0_WP
    evaporation = 0.0_WP
    fresh_water_flux = 0.0_WP
    real_salt_flux = 0.0_WP
    net_heat_flux = 0.0_WP
    allocate(mask_ice_adv_n(node_size),mask_ice_adv_e(elem_size), STAT=sbc_alloc )
    if( sbc_alloc /= 0 )   STOP 'array_setup: failed to allocate arrays mask_ice_adv_n'
    mask_ice_adv_n = 1.0_WP
    mask_ice_adv_e = 1.0_WP
    allocate(mask_ice_adv_big_n(node_size),mask_ice_adv_big_e(elem_size), STAT=sbc_alloc )
    if( sbc_alloc /= 0 )   STOP 'array_setup: failed to allocate arrays mask_ice_adv_big_n'
    mask_ice_adv_big_n = 1.0_WP
    mask_ice_adv_big_e = 1.0_WP

  endif
  allocate(C_d_2d_n(node_size), C_d_2d_e(elem_size), STAT=sbc_alloc )
  if( sbc_alloc /= 0 )   STOP 'array_setup: failed to allocate arrays mask_ice_adv_big_n'
  C_d_2d_n = 0.003_WP
  C_d_2d_e = 0.003_WP

  allocate(mslp(node_size), STAT=sbc_alloc )
  mslp = P_ref
  if( sbc_alloc /= 0 )   STOP 'array_setup: failed to allocate arrays'

  allocate(U_n_2D(2, elem_size), U_n_1(2, elem_size), U_n_2(2, elem_size))
  !NR  Allocate(U_n_2D_old(2,elem_size))  !NR not used
  U_n_2D=0.0_WP
  U_n_1=0.0_WP
  U_n_2=0.0_WP
  !NR  U_n_2D_old=0.0_WP

  allocate(dmean_n(elem_size))
  allocate(eta_p(node_size))
  eta_p=0.0_WP

  allocate(UAB(2, elem_size), U_rhs_2D(2, elem_size))
  UAB=0.0_WP
  U_rhs_2D=0.0_WP
  allocate(U_rhs_2D_3D(2,elem_size))
  U_rhs_2D_3D=0.0_WP
  allocate( UV2_rhs(2,elem_size))
  UV2_rhs = 0.0_WP

  allocate(eta_n(node_size), eta_n_1(node_size), eta_n_2(node_size))
  eta_n=0.0_WP
  eta_n_1=0.0_WP
  eta_n_2=0.0_WP

  allocate(etaAB(node_size),ssh_rhs(node_size), ssh_gp(node_size))
  allocate(vel_grad(4,elem_size))
  U_rhs_2D=0.0_WP
  ssh_rhs=0.0_WP
  vel_grad=0.0_WP
  ssh_gp=0.0_WP

  allocate(taux(elem_size), tauy(elem_size))
  allocate(taux_node(node_size), tauy_node(node_size))
  taux=0.0_WP
  tauy=0.0_WP
  taux_node=0.0_WP
  tauy_node=0.0_WP

  allocate(windx(node_size), windy(node_size), wind(node_size))
  windx=0.0_WP
  windy=0.0_WP
  allocate(qns(node_size), emp(node_size), qsr(node_size))
  qns = 0.0_WP
  emp = 0.0_WP
  qsr = 0.0_WP

  allocate(relax_coef(node_size))
  relax_coef=0.0_WP

  ! add 2D part
  allocate(Bar_pr_2D(2,elem_size),hpre_2D(node_size))
  Bar_pr_2D=0.0_WP
  hpre_2D=0.0_WP

  ! 3D part

  if (type_task>1) then
     allocate(TF(nsigma-1,node_size), SF(nsigma-1,node_size))
     allocate(T_old(nsigma-1,node_size), S_old(nsigma-1,node_size))
     allocate(Tclim(nsigma-1,node_size), Sclim(nsigma-1,node_size))
     TF(:,:)=0.0_WP; SF(:,:)=0.0_WP
  endif
  if (comp_sediment) allocate(CF(nsigma-1,node_size),c_old(nsigma-1,node_size),w_s(nsigma-1,node_size),Cclim(nsigma-1,node_size))
  allocate(rho_c(nsigma-1,node_size))
  rho_c=0.0_WP

  Allocate(C_d_el(elem_size))

  if (type_task>1) then
     allocate(z0b_gotm_el(elem_size))
     z0b_gotm_el=z0b_min
     allocate(U_n(nsigma-1,elem_size), V_n(nsigma-1,elem_size))
     U_n=0.0_WP
     V_n=0.0_WP
     allocate(hpressure(nsigma,node_size))
     hpressure=0.0_WP
     allocate(Bar_pru_3D(nsigma-1,elem_size), Bar_prv_3D(nsigma-1,elem_size))
     Bar_pru_3D=0.0_WP
     Bar_prv_3D=0.0_WP
     allocate(Bar_pru_3D_clim(nsigma-1,elem_size), Bar_prv_3D_clim(nsigma-1,elem_size))
     Bar_pru_3D_clim=0.0_WP
     Bar_prv_3D_clim=0.0_WP
     allocate(snu(nsigma,node_size), bt(nsigma,node_size),Ri(nsigma,node_size),tke_dissip(nsigma,node_size))
     Ri = 0.0_WP
     snu = snul
     bt = 1.0d-8
     tke_dissip = 0.0_WP
     allocate(UV_rhs(2,nsigma-1,elem_size))
     UV_rhs = 0.0_WP

     Allocate(Jc(nsigma-1,node_size),Jc_old(nsigma-1,node_size))
     Allocate(Je(nsigma-1,elem_size),Jd(nsigma-1,myDim_edge2D+eDim_edge2D))
     Jc=0.0_WP
     Jc_old=0.0_WP
     Je=0.0_WP
     Jd=0.0_WP
     Allocate(Unode(nsigma-1,node_size), Vnode(nsigma-1,node_size))
     Unode=0.0_WP
     Vnode=0.0_WP
     allocate(Kv(nsigma,node_size), Av_node(nsigma,node_size),Av(nsigma,elem_size))
     allocate(L(nsigma,node_size), teps(nsigma,node_size), tepsb(nsigma,node_size),tke(nsigma,node_size))
     tke=1.0e-8
     teps=1.0e-12
     Kv=snul
     Av=snul
     Av_node=snul
     allocate(U_rhs(nsigma-1,elem_size), V_rhs(nsigma-1,elem_size))
     Allocate(U_rhsAB(nsigma-1,elem_size), V_rhsAB(nsigma-1,elem_size))
     allocate(vel_grad_ux(nsigma-1,elem_size), vel_grad_uy(nsigma-1,elem_size))
     allocate(vel_grad_vx(nsigma-1,elem_size), vel_grad_vy(nsigma-1,elem_size))
     vel_grad_ux=0.0_WP
     vel_grad_uy=0.0_WP
     vel_grad_vx=0.0_WP
     vel_grad_vy=0.0_WP
     allocate(Visc(nsigma-1,elem_size))
     Visc=0.0_WP

     Allocate(U_puls(nsigma-1,elem_size), V_puls(nsigma-1,elem_size))
     Allocate(Unode_p(nsigma-1,node_size), Vnode_p(nsigma-1,node_size))
     Allocate(U_n_filt(nsigma-1,elem_size), V_n_filt(nsigma-1,elem_size))

     Allocate(W_n(nsigma,node_size))
     Unode_p=0.0_WP
     Vnode_p=0.0_WP
     U_puls=0.0_WP
     V_puls=0.0_WP
     U_n_filt=0.0_WP
     V_n_filt=0.0_WP
     U_rhsAB=0.0_WP
     V_rhsAB=0.0_WP
     U_rhs=0.0_WP
     V_rhs=0.0_WP
     W_n=0.0_WP

     Allocate(vorticity_2D(node_size), Vorticity_3D(nsigma-1,node_size))

     allocate(Z(nsigma-1,node_size),zbar(nsigma,node_size))

     Allocate(Wvel(nsigma,node_size))
     Wvel=0.0_WP

     allocate(edge_up_dn_grad(4,nsigma-1,myDim_edge2D+eDim_edge2D))

     edge_up_dn_grad=0.0_WP

  endif

  allocate(Visc2D(elem_size))
  Visc2D=0.0_WP
  Allocate(edge_up_dn_tri(2,myDim_edge2D+eDim_edge2D))
  edge_up_dn_tri=0
  allocate(relax2clim(node_size))
  relax2clim = 0.0_WP          ! relaxation to climatology
  Allocate(U_filt_2D(elem_size), V_filt_2D(elem_size))
  U_filt_2D=0.0_WP
  V_filt_2D=0.0_WP
  allocate(mask_wd(elem_size))
  mask_wd = 1.0_WP          ! wet elements
  allocate(mask_wd_node(node_size))
  mask_wd_node = 1.0_WP          ! wet elements
  Allocate(mask_ad(node_size))
  mask_ad=1.0_WP            ! mask for order advection scheme for tracer
  ! 0.0  ---> upwind
  ! 1.0 ---> 2 order (MIURA....)
  ! 0.0 ---> 1.0 reduce accuracy near critical depth
  ! H_advtr_crit  ---> if full depth > H_adv_crit ---> high order scheme for tracer
  ! H_advtr_min ---> if full depth < H_advtr_min ---> UPWIND
  ALLOCATE(mask_bpr(elem_size))
  mask_bpr=1.0_WP

  !++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! SEDIMENT
  !++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! remove comment,
!  if (comp_sediment) then
  allocate(qbu(elem_size), qbv(elem_size))
  qbu = 0.0_WP
  qbv = 0.0_WP
  Allocate(hama_v(node_size), E_sed(node_size), h_var(node_size), h_var_old(node_size))
  Allocate(Er_Dep(node_size), h_var2(node_size), h_var_old2(node_size))
  Allocate(qb(2,node_size))
  qb=0.0_WP
  hama_v = 0.0_WP
  E_sed = 0.0_WP
  Er_Dep = 0.0_WP
  h_var = 1.0_WP
  h_var_old = 1.0_WP
  h_var2 = 1.0_WP
  h_var_old2 = 1.0_WP
  allocate(con(node_size), con_bc(node_size))
  con = 0.0_WP
  con_bc = 0.0_WP
!  endif
  if (mype==0) write(*,*) 'Arrays are set up'

END SUBROUTINE array_setup

!====================================================================================

subroutine jacobian

  !++++++++++++++++++++++++++++++++++
  ! Compute vertical Jacobian
  ! 19.09.2014
  ! Androsov Alexey
  !++++++++++++++++++++++++++++++++++

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  USE g_PARSUP
  use g_comm_auto

  IMPLICIT NONE

  INTEGER        :: n, nz, el, ed
  real(kind=WP)  :: a

!$OMP DO
  do n=1,myDim_nod2D+eDim_nod2D

     Jc_old(1:nsigma-1,n) = Jc(1:nsigma-1,n)

     a = max(depth(n) + eta_n(n), Dmin)

     zbar(1:nsigma,n) = (1.0_WP - sigma(1:nsigma))*a

     Jc(1:nsigma-1,n) = a*( sigma(1:nsigma-1) - sigma(2:nsigma))

     Z(1:nsigma-1,n)  = 0.5_WP*(zbar(1:nsigma-1,n) + zbar(2:nsigma,n))

  end do
!$OMP END DO

!$OMP DO
  do el=1,myDim_elem2D

     Je(1:nsigma-1,el) = w_cv(1,el)*Jc(1:nsigma-1,elem2D_nodes(1,el)) &
                       + w_cv(2,el)*Jc(1:nsigma-1,elem2D_nodes(2,el)) &
                       + w_cv(3,el)*Jc(1:nsigma-1,elem2D_nodes(3,el)) &
                       + w_cv(4,el)*Jc(1:nsigma-1,elem2D_nodes(4,el))

  end do
!$OMP END DO NOWAIT

#ifdef USE_MPI
  call exchange_elem(Je)
#endif

!$OMP DO
  DO ed=1, myDim_edge2D+eDim_edge2D
     Jd(1:nsigma-1,ed) = 0.5_WP*(Jc(1:nsigma-1,edge_nodes(1,ed)) + Jc(1:nsigma-1,edge_nodes(2,ed)))
  enddo
!$OMP END DO

!IK  if (mype==0) print *,'SUBROUTINE jacobian COMPLETED'

End subroutine jacobian
!===========================================================================
subroutine wad_mask

  ! IK: Calculate masks (1/0) for wetting and drying on elemens and on nodes.

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  use g_parsup
  use g_comm_auto

  IMPLICIT NONE

  integer         :: el,nod

!$OMP DO
  do el=1,myDim_elem2D

     !NR  .false. = dry element, .true. = wet element

     mask_wd(el) = merge(1.0_WP,0.0_WP,( min(depth(elem2D_nodes(1,el)),depth(elem2D_nodes(2,el)),   &
                                             depth(elem2D_nodes(3,el)),depth(elem2D_nodes(4,el)))   &
                                       + max(etaAB(elem2D_nodes(1,el)),etaAB(elem2D_nodes(2,el)),   &
                                             etaAB(elem2D_nodes(3,el)),etaAB(elem2D_nodes(4,el))) > Dmin ) )

  enddo
!$OMP END DO NOWAIT

#ifdef USE_MPI
  call exchange_elem(mask_wd)
#endif

!$OMP DO
  do nod=1,mydim_nod2D+eDim_nod2D
     mask_wd_node(nod) = merge(1.0_WP,0.0_WP,( depth(nod) + etaAB(nod) > Dmin))
  enddo
!$OMP END DO

end subroutine wad_mask
!===========================================================================
subroutine adv_tracer_mask

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  use g_parsup

  IMPLICIT NONE

  real(kind=WP)   :: dmean,x,y
  integer              :: n

  mask_ad = 1.0_WP

  do n=1,myDim_nod2D+eDim_nod2D
     dmean = max(Dmin,(depth(n) + eta_n(n)))
     if (dmean > H_advtr_crit) then
        mask_ad(n) = 1.0_WP  ! high order advection for tracer
     elseif (dmean < H_advtr_min) then
        mask_ad(n) = 0.0_WP  ! first order advection for tracer
     else
        mask_ad(n) = (dmean - H_advtr_min)/(H_advtr_crit - H_advtr_min)
     endif
  enddo

end subroutine adv_tracer_mask

!===========================================================================

subroutine bpr_mask

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  use g_parsup

  IMPLICIT NONE

  integer         :: elem

  mask_bpr = 1.0_WP

  do elem=1,myDim_elem2D !SH +eDim_elem2D+eXDim_elem2D

     dmean_n(elem) = max(Dmin, w_cv(1,elem)*(eta_n(elem2D_nodes(1,elem)) + depth(elem2D_nodes(1,elem))) &
          +  w_cv(2,elem)*(eta_n(elem2D_nodes(2,elem)) + depth(elem2D_nodes(2,elem))) &
          +  w_cv(3,elem)*(eta_n(elem2D_nodes(3,elem)) + depth(elem2D_nodes(3,elem))) &
          +  w_cv(4,elem)*(eta_n(elem2D_nodes(4,elem)) + depth(elem2D_nodes(4,elem))))

     if (dmean_n(elem) > H_bpr_crit) then
        mask_bpr(elem) = 1.0_WP  ! compute Baroclinic pressure
     elseif (dmean_n(elem) < H_bpr_min) then
        mask_bpr(elem) = 0.0_WP  ! without Baroclinic pressure
     else
        mask_bpr(elem) = (dmean_n(elem) - H_bpr_min)/(H_bpr_crit - H_bpr_min)
     endif
  enddo

end subroutine bpr_mask
!===========================================================================

SUBROUTINE timestep_AB_2D(step)

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  use g_parsup
  use g_comm_auto

  IMPLICIT NONE

  real(kind=WP)   :: dmean, fD(4), fDD(4), ibv, rho_inv
  integer         :: step, elem, i, elnodes(4),k,ed
  integer         :: n, el
  real(kind=WP)   :: t_0, t_1, t_2, t_3, t_4, t_5, t_6
  ! No flux form, only U\nablau

  ! =============
  !  AB3 interpolate
  ! =============

!aa13.02.20 print *,'step, U_n_2D :',mype,step,minval(U_n_2D(2,:)),maxval(U_n_2D(2,:))

#ifdef USE_MPI
 t_0=MPI_Wtime()
#endif


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(n,el)
!$OMP DO
  do n=1,myDim_nod2D+eDim_nod2D
     etaAB(n)  = (1.5_WP+beta)*eta_n(n) - (0.5_WP+2.0_WP*beta)*eta_n_1(n)  + beta*eta_n_2(n)
  enddo
!$OMP END DO NOWAIT
!$OMP DO
  do el=1,myDim_elem2D+eDim_elem2D+eXDim_elem2D
     UAB(1,el) = (1.5_WP+beta)*U_n_2D(1,el) - (0.5_WP+2.0_WP*beta)*U_n_1(1,el) + beta*U_n_2(1,el)
     UAB(2,el) = (1.5_WP+beta)*U_n_2D(2,el) - (0.5_WP+2.0_WP*beta)*U_n_1(2,el) + beta*U_n_2(2,el)
  enddo
!$OMP END DO
!$OMP END PARALLEL

!#ifdef USE_MPI
! t_1=MPI_Wtime()
!#endif

  !==================
  !  compute ssh_rhs
  !==================

  call compute_ssh_rhs_elem

!#ifdef USE_MPI
! t_2=MPI_Wtime()
!#endif


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
!$OMP DO

  do i=1,myDim_nod2D
     !Calculate eta_n+1
     if (index_nod2D(i) < 2) ssh_rhs(i) = eta_n(i) + dt_2D*ssh_rhs(i)/area(i)
  enddo

#ifdef USE_MPI
  call exchange_nod(ssh_rhs)
#endif
!if (mype==12) then
!   write(*,*) ssh_rhs(my_in_obn(1)),ssh_rhs(my_in_obn(10)),ssh_rhs(my_in_obn(1))-ssh_rhs(my_in_obn(10))
!endif


!$OMP END DO
!$OMP END PARALLEL

  !##########################################################
  !##

  IF ( (step==1) .and. (type_task>1) ) THEN
     !================================================
     ! right hand from 3D advection and diffusion to 2D
     !================================================
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(elem,i)
!$OMP DO
     DO elem=1,myDim_elem2D+eDim_elem2D+eXDim_elem2D
        U_rhs_2D_3D(1,elem) =0.0_WP
        U_rhs_2D_3D(2,elem) =0.0_WP
     ENDDO
!$OMP END DO NOWAIT
     !++++++++++++++++++++++++++++++++++
     ! compute depth and ssh
     ! on predictor step
     !++++++++++++++++++++++++++++++++++
!$OMP DO
     DO elem=1,myDim_elem2D

        dmean_n(elem) = max(Dmin, w_cv(1,elem)*(eta_n(elem2D_nodes(1,elem)) + depth(elem2D_nodes(1,elem))) &
             +  w_cv(2,elem)*(eta_n(elem2D_nodes(2,elem)) + depth(elem2D_nodes(2,elem))) &
             +  w_cv(3,elem)*(eta_n(elem2D_nodes(3,elem)) + depth(elem2D_nodes(3,elem))) &
             +  w_cv(4,elem)*(eta_n(elem2D_nodes(4,elem)) + depth(elem2D_nodes(4,elem))))
     END DO
!$OMP END DO NOWAIT

#ifdef USE_MPI
     call exchange_elem(dmean_n)
#endif

!aa13.02.20 print *,'ZWEITERTEST',mype,minval(dmean_n),maxval(dmean_n)
!aa13.02.20 print *,'ERSTERaaa',step,minval(ssh_rhs),maxval(ssh_rhs)


!$OMP DO
     DO i=1,myDim_nod2D+eDim_nod2D
        eta_p(i) = eta_n(i)
     END DO
!$OMP END DO
!aa13.02.20 print *,'ERSTERbbbb',mype,minval(ssh_rhs),maxval(ssh_rhs)

     call compute_puls_vel

     call compute_el2nodes_3D(U_puls,V_puls,Unode_p,Vnode_p)

!$OMP END PARALLEL

     call momentum_adv_P1_3D_to_2D

     if (filt_3D)  call viscosity_filt_3D_to_2D

     if (bih_3D)  call biharmonic_viscosity_3D_to_2D

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(elem)
     DO elem=1,myDim_elem2D
        !NR     elnodes=elem2D_nodes(:,elem)
        !NR     fD=depth(elnodes) + eta_n(elnodes)
        !NR     dmean = max(Dmin,sum(w_cv(1:4,elem)*fD))
        !NR dmean ist already calculated as dmean_n(elem), see loop above

        U_rhs_2D_3D(1,elem) = U_rhs_2D_3D(1,elem)/dmean_n(elem)
        U_rhs_2D_3D(2,elem) = U_rhs_2D_3D(2,elem)/dmean_n(elem)

     END DO
!$OMP END PARALLEL DO

#ifdef USE_MPI
     call exchange_elem(U_rhs_2D_3D)
#endif

  END IF ! part for step==1 only

!#ifdef USE_MPI
! t_3=MPI_Wtime()
!#endif

  !##
  !##########################################################


!aa13.02.20 print *,'HUHU U',step,minval(U_rhs_2D(1,:)),maxval(U_rhs_2D(1,:))
!aa13.02.20 print *,'HUHU V',step,minval(U_rhs_2D(2,:)),maxval(U_rhs_2D(2,:))


  call update_2D_vel(step)
!aa13.02.20 print *,'HUHU2',step,minval(U_rhs_2D(2,:)),maxval(U_rhs_2D(2,:))

!#ifdef USE_MPI
! t_4=MPI_Wtime()
!#endif

  ! =============
  ! Update elevation and velocity
  ! =============

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,n,elem,dmean)
!$OMP DO
  do n=1,myDim_nod2D
     eta_n_2(n) = eta_n_1(n)
     eta_n_1(n) = eta_n(n)
     eta_n(n)   = ssh_rhs(n)
  end do
  !$OMP END DO NOWAIT

#ifdef USE_MPI
  call exchange_nod(eta_n_1)
  call exchange_nod(eta_n_2)
  call exchange_nod(eta_n)
#endif

!aa13.02.20 print *,'ETA_N:',step,minval(eta_n),maxval(eta_n)
!aa13.02.20 print *,'WCV:',step,minval(w_cv),maxval(w_cv)
!aa13.02.20 print *,'U_FILT:',step,minval(U_filt_2D),maxval(U_filt_2D)
!aa13.02.20 print *,'U_RHS:',step,minval(U_rhs_2D(1,:)),maxval(U_rhs_2D(1,:))
!aa13.02.20 print *,'V_RHS:',step,minval(U_rhs_2D(2,:)),maxval(U_rhs_2D(2,:))



!$OMP DO
  do elem=1,myDim_elem2D
     U_n_2(1,elem)      = U_n_1(1,elem)
     U_n_2(2,elem)      = U_n_1(2,elem)
     U_n_1(1,elem)      = U_n_2D(1,elem)
     U_n_1(2,elem)      = U_n_2D(2,elem)
     U_n_2D(1,elem)     = U_rhs_2D(1,elem)
     U_n_2D(2,elem)     = U_rhs_2D(2,elem)
!NR  U_n_2D_old(1,elem) = U_n_2D(1,elem)  !NR never used
!NR  U_n_2D_old(2,elem) = U_n_2D(2,elem)  !NR never used
  enddo
!$OMP END DO NOWAIT

#ifdef USE_MPI
  call exchange_elem(U_n_1)
  call exchange_elem(U_n_2)
  call exchange_elem(U_n_2D)
#endif

!#ifdef USE_MPI
! t_5=MPI_Wtime()
!#endif
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  preparation for compute Filtering 2D velocity
!
!$OMP DO
  do elem=1,myDim_elem2D

    dmean = max(Dmin, w_cv(1,elem)*(eta_n(elem2D_nodes(1,elem)) + depth(elem2D_nodes(1,elem))) &
                   +  w_cv(2,elem)*(eta_n(elem2D_nodes(2,elem)) + depth(elem2D_nodes(2,elem))) &
                   +  w_cv(3,elem)*(eta_n(elem2D_nodes(3,elem)) + depth(elem2D_nodes(3,elem))) &
                   +  w_cv(4,elem)*(eta_n(elem2D_nodes(4,elem)) + depth(elem2D_nodes(4,elem))))

!NR U_rhs_2D allows the loop above to "NOWAIT"
!NR    U_filt_2D(elem) = U_filt_2D(elem) + U_n_2D(1,elem)*dmean
!NR    V_filt_2D(elem) = V_filt_2D(elem) + U_n_2D(2,elem)*dmean
    U_filt_2D(elem) = U_filt_2D(elem) + U_rhs_2D(1,elem)*dmean
    V_filt_2D(elem) = V_filt_2D(elem) + U_rhs_2D(2,elem)*dmean
  enddo
!$OMP END DO
!$OMP END PARALLEL

#ifdef USE_MPI
  call exchange_elem(U_filt_2D)
  call exchange_elem(V_filt_2D)
#endif

!#ifdef USE_MPI
! t_6=MPI_Wtime()
!#endif

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! VF calling of sediment procedure is shifted here,
! so do not spoil parallelization
! AA compute sediment model
!SH SKIPPING FOR NOW    If (comp_sediment) call sediment
! AA

#ifdef USE_MPI
 if(mype==0) then
  !write(91,*) t_1-t_0,t_2-t_1,t_3-t_2,t_4-t_3,t_5-t_4,t_6-t_5
  !write(91,*) t_4-t_3,t_2-t_1,t_3-t_0
 endif
#endif


END SUBROUTINE timestep_AB_2D
!===========================================================================
SUBROUTINE oce_barotropic_timestep

  !a modification by Alexey Androsov
  !a (for sigma coordinates)
  !a 14.10.14

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  USE g_PARSUP
  use g_comm_auto

  IMPLICIT NONE

  integer       :: n, el
  real(kind=WP) :: Mt_inv

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(el)
  do el=1,myDim_elem2D+eDim_elem2D+eXDim_elem2D
     U_filt_2D(el) = 0.0_WP
     V_filt_2D(el) = 0.0_WP
  enddo
!$OMP END PARALLEL DO

!if (mype==0) print *,'Mt =',Mt
  do n=1,Mt
     time_2D=time + dt_2D*n
     call timestep_AB_2D(n)
!aa13.02.20 write (*,*) 'After AB: ',n,minval(U_filt_2D),maxval(U_filt_2D)
  enddo

  ! Filterin velocity

  Mt_inv = 1._WP/real(Mt,WP)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(el)
  do el=1,myDim_elem2D
     U_filt_2D(el) = U_filt_2D(el) * Mt_inv
     V_filt_2D(el) = V_filt_2D(el) * Mt_inv
  enddo
!$OMP END PARALLEL DO

!SH Necessary?
#ifdef USE_MPI
  call exchange_elem(U_filt_2D)
  call exchange_elem(V_filt_2D)
#endif

#ifdef DEBUG

#ifdef USE_MPI
!     call MPI_AllREDUCE(minval(U_filt_2D),mnv, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
!          MPI_COMM_FESOM_C, MPIerr)

!     call MPI_AllREDUCE(maxval(U_filt_2D),mxv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
!          MPI_COMM_FESOM_C, MPIerr)

!aa13.02.20     if (mype==0) print *,'MIN_VAL Ufilt',mnv
!aa13.02.20     if (mype==0) print *,'MAX_VAL Ufilt',mxv
#else
!     print *,'MIN_VAL Ufilt',minval(U_filt_2D)
!     print *,'MAX_VAL Ufilt',maxval(U_filt_2D)
#endif

#endif

end subroutine oce_barotropic_timestep

!===================================================================================

SUBROUTINE solve_tracers

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  use g_parsup

  IMPLICIT NONE

  real(kind=WP) :: T_min, T_max, S_min, S_max
  real(kind=WP) :: T_min_loc, T_max_loc, S_min_loc, S_max_loc
  integer       :: n, nz, node
  real(kind=WP) :: t_0, t_1, t_2, t_3, t_4
  real(kind=WP) :: river_coeff, river_delta

  !a if(tracer_adv==1) then
  !a call solve_tracer_miura(TF, 't')
  !a call solve_tracer_miura(SF, 's')
  !a end if

  !a if(tracer_adv==2) then
  !a call solve_tracer_second_order_rec(TF, 't')
  !a call solve_tracer_second_order_rec(SF, 's')
  !a end if

#ifdef USE_MPI
 t_0=MPI_Wtime()
#endif

  T_min_loc = TF(1,1)
  T_max_loc = TF(1,1)
  S_min_loc = SF(1,1)
  S_max_loc = SF(1,1)


!$OMP PARALLEL REDUCTION(min:T_min, S_min) &
!$OMP&         REDUCTION(max:T_max, S_max)
!$OMP DO
  do n=1,myDim_nod2D
     T_min_loc = min(T_min_loc,minval(TF(:,n)))
     T_max_loc = max(T_max_loc,maxval(TF(:,n)))
  enddo
!$OMP END DO NOWAIT
!$OMP DO

  do n=1,myDim_nod2D
     S_min_loc = min(S_min_loc,minval(SF(:,n)))
     S_max_loc = max(S_max_loc,maxval(SF(:,n)))
  enddo
!$OMP END DO
!$OMP END PARALLEL

#ifdef USE_MPI
  call MPI_AllREDUCE(T_max_loc,T_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
       MPI_COMM_FESOM_C, MPIerr)
  call MPI_AllREDUCE(T_min_loc,T_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
       MPI_COMM_FESOM_C, MPIerr)
  call MPI_AllREDUCE(S_max_loc,S_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
       MPI_COMM_FESOM_C, MPIerr)
  call MPI_AllREDUCE(S_min_loc,S_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
       MPI_COMM_FESOM_C, MPIerr)
#else
  T_max=T_max_loc
  T_min=T_min_loc
  S_max=S_max_loc
  S_min=S_min_loc
#endif

  T_aver=0.5_WP*(T_min + T_max)
  S_aver=0.5_WP*(S_min + S_max)

#ifdef USE_MPI
 t_1=MPI_Wtime()
#endif

  call adv_tracer_mask

#ifdef USE_MPI
 t_2=MPI_Wtime()
#endif

  !write(*,*) 'update T_aver, S_aver'
  if(tracer_adv==2) then
!SH SKIPPED FOR NOW     call solve_tracer_muscl(TF, T_old, 't')
!SH SKIPPED FOR NOW     call solve_tracer_muscl(SF, S_old, 's')
  elseif(tracer_adv==1) then
     call solve_tracer_upwind(TF, T_old, SF, S_old)

  elseif (tracer_adv==4) then

     call solve_tracer_miura(TF, T_old, 't')
!SH SKIPPED FOR NOW     call solve_tracer_miura(SF, S_old, 's')
  end if

#ifdef USE_MPI
 t_3=MPI_Wtime()
#endif
  ! add salinity (remove) due to rivers
  ! runoff_rivers m3/s (interpolated in time runoff, see rivers_do, fv_rivers.f90)
  ! my_Nnods_rivers number of nodes with rivers
  ! my_nod_rivers array of nodes indexes
  !    runoff * delta Sigma for current layer * timestep / area / cell height
  if (key_rivers) then
     ! to optimize: precalculate river_coeff and runoff_rivers(i)*(sigma(nz) - sigma(nz+1)))
     do n=1, my_Nnods_rivers
        do nz=1,nsigma-1
           !0.2 - salinity in the river
           node = my_nod_rivers(n)
           river_coeff = runoff_rivers(n)*(sigma(nz) - sigma(nz+1)) &
                         * dt / (area(node)*Jc(nz,node))
           river_delta = (0.2 - S_old(nz,node))
           SF(nz,node) = SF(nz,node) + river_delta * river_coeff
        enddo
     enddo
  endif

  call tracer_impl_vert_diff

#ifdef USE_MPI
 t_4=MPI_Wtime()
#endif

#ifdef USE_MPI
 if(mype==0) then
  write(92,*) t_1-t_0,t_2-t_1,t_3-t_2,t_4-t_3
 endif
#endif

end subroutine solve_tracers
!==========================================================================
!

!===================================================================================
SUBROUTINE solve_tracers_sed
USE o_MESH
USE o_ARRAYS
USE o_PARAM
IMPLICIT NONE
real(kind=WP) :: c_min, c_max
integer       :: n

c_min = CF(1,1)
c_max = CF(1,1)

!$OMP PARALLEL REDUCTION(min:c_min) &
!$OMP&         REDUCTION(max:c_max)
!$OMP DO
do n=1,nod2D
   c_min = min(c_min,minval(CF(:,n)))
   c_max = max(c_max,maxval(CF(:,n)))
enddo
!$OMP END DO
!$OMP END PARALLEL
c_aver=0.5_WP*(c_min + c_max)


call adv_tracer_mask
! IK adv_tracer_mask does not work after step back (co ver 429) with old/new fv_tracer
STOP "ERROR: check here,fv_main,solve_tracers_sed,call adv_tracer_mask"

 call solve_tracer_upwind_sed(CF, c_old)

 call tracer_impl_vert_diff_sed
 call bottom_evolution2

end subroutine solve_tracers_sed

!==========================================================================

SUBROUTINE oce_timestep

  !a modification by Alexey Androsov
  !a (for sigma coordinates)
  !a 13.10.14

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  USE g_PARSUP

  IMPLICIT NONE

  real(kind=WP)      :: t_0,t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, t_9, t_10

  integer                :: tracer_sch, n, nz


!aa13.02.20  print *,mype,'oce_timestep started!'


  if (T_potential) call potential           ! AA compute tidal geopotential here

  select case (type_task)

  case(1)

!!$     !$ if (iverbosity >= 3) t1=omp_get_wtime()
     call oce_barotropic_timestep

!!$     !$ if (iverbosity >= 3) then
!!$     !$    t2=omp_get_wtime()
!!$     !$    write(*,'("oce_barotropic_timestep took ",f10.4,"s")') t2-t1
!!$     !$ endif

  case(2)

     !$ if (iverbosity >= 3) t1=omp_get_wtime()
!$OMP PARALLEL
!SHDB print *,'ARGH ',minval(U_n),maxval(U_n)
!SHDB print *,'ARGH ',minval(V_n),maxval(V_n)

     call compute_el2nodes_3D(U_n,V_n,Unode,Vnode)
!$OMP END PARALLEL
     !+++++++++++++++++++++++++
     ! vertical mixing scheme
     !+++++++++++++++++++++++++
     !$ if (iverbosity >= 3) t2=omp_get_wtime()

!SH SKIPPED     if (ver_mix == 1) call oce_mixing_PP
!SH SKIPPED     if (ver_mix == 2) call GOTM
     if (ver_mix == 3) then
        if (len(trim(title))==2 .AND. title=='LE') then
           call d3_end_LE
        else
           call d3_end
        end if
     end if

     !$ if (iverbosity >= 3) t3=omp_get_wtime()

     call compute_vel_rhs

     !$ if (iverbosity >= 3) t4=omp_get_wtime()
     call oce_barotropic_timestep
     !$ if (iverbosity >= 3) t5=omp_get_wtime()

!$OMP PARALLEL
     call jacobian

     !$ if (iverbosity >= 3) t6=omp_get_wtime()
     call update_3D_vel

     !$ if (iverbosity >= 3) t7=omp_get_wtime()
     call vert_vel_sigma

!$OMP END PARALLEL

     !$ if (iverbosity >= 3) then
     !$   t8=omp_get_wtime()
     !$    write(*,'("compute_vel_nodes took       ",f10.4," s")') t2-t1
     !$    write(*,'("mixing took                  ",f10.4," s")') t3-t2
     !$    write(*,'("compute_vel_rhs took         ",f10.4," s")') t4-t3
     !$    write(*,'("oce_barotropic_timestep took ",f10.4," s")') t5-t4
     !$    write(*,'("jacobian took                ",f10.4," s")') t6-t5
     !$    write(*,'("update_3D_vel took           ",f10.4," s")') t7-t6
     !$    write(*,'("vert_vel_sigma took          ",f10.4," s")') t8-t7
     !$ endif

     if (comp_sediment) then
        !VF If comp_sediment and no other tracers, the density still changes due to varying concentration
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n, nz)
        DO n=1,nod2D
           ! VF, compute density in the vertical column and add it to rho_c - current rho, needed in fv_tracer
           ! VF, add effect of the susp. sediments
           DO nz=1,nsigma-1
!              call densityJM(T_const, S_const, -Z(nz,n), rho_c(nz,n),CF(nz,n))
           END DO
        END DO
        !$OMP END PARALLEL DO

        call solve_tracers_sed

     endif

  case(3)

!#ifdef USE_MPI
! t_0 = MPI_Wtime()
!#endif
     if (Mask_Bar_pr) call bpr_mask

     !$ if (iverbosity >= 3) t1=omp_get_wtime()

!#ifdef USE_MPI
! t_1 = MPI_Wtime()
!#endif

     call pressure
     ! if ((time/3600.0_WP/24.0_WP)>4) call pressure

     !$ if (iverbosity >= 3) t2=omp_get_wtime()

!#ifdef USE_MPI
! t_2 = MPI_Wtime()
!#endif

!$OMP PARALLEL
     call compute_el2nodes_3D(U_n,V_n,Unode,Vnode)

!#ifdef USE_MPI
! t_3 = MPI_Wtime()
!#endif

!$OMP END PARALLEL

     !$ if (iverbosity >= 3) t3=omp_get_wtime()

     !+++++++++++++++++++++++++
     ! vertical mixing scheme
     !+++++++++++++++++++++++++
!SH SKIPPED FOR NOW     if (ver_mix == 1) call oce_mixing_PP
!SH SKIPPED FOR NOW     if (ver_mix == 2) call GOTM
     if (ver_mix == 3) then
        if (len(trim(title))==2 .AND. title=='LE') then
           call d3_end_LE
        else
           call d3_end
        end if
     end if

!#ifdef USE_MPI
! t_4 = MPI_Wtime()
!#endif

     !$ if (iverbosity >= 3) t4=omp_get_wtime()
     call compute_vel_rhs

!#ifdef USE_MPI
! t_5 = MPI_Wtime()
!#endif

     !$ if (iverbosity >= 3) t5=omp_get_wtime()
     call oce_barotropic_timestep

!#ifdef USE_MPI
! t_6 = MPI_Wtime()
!#endif

     !$ if (iverbosity >= 3) t6=omp_get_wtime()

!$OMP PARALLEL
     call jacobian

!#ifdef USE_MPI
! t_7 = MPI_Wtime()
!#endif

     !$ if (iverbosity >= 3) t7=omp_get_wtime()
     call update_3D_vel

!#ifdef USE_MPI
! t_8 = MPI_Wtime()
!#endif

     !$ if (iverbosity >= 3) t8=omp_get_wtime()
     call vert_vel_sigma
!#ifdef USE_MPI
! t_9 = MPI_Wtime()
!#endif


!$OMP END PARALLEL


     !$ if (iverbosity >= 3) t9=omp_get_wtime()

     !   print *,'================================ CALL TRACERS ===================='
     call solve_tracers
!#ifdef USE_MPI
! t_10 = MPI_Wtime()
!#endif


!SH SKIPPED FOR NOW     if (comp_sediment) call solve_tracers_sed

     !$ if (iverbosity >= 3) then
     !$   t10=omp_get_wtime()
     !$    write(*,'("pressure took                ",f10.4," s")') t2-t1
     !$    write(*,'("compute_vel_nodes took       ",f10.4," s")') t3-t2
     !$    write(*,'("mixing took                  ",f10.4," s")') t4-t3
     !$    write(*,'("compute_vel_rhs took         ",f10.4," s")') t5-t4
     !$    write(*,'("oce_barotropic_timestep took ",f10.4," s")') t6-t5
     !$    write(*,'("jacobian took                ",f10.4," s")') t7-t6
     !$    write(*,'("update_3D_vel took           ",f10.4," s")') t8-t7
     !$    write(*,'("vert_vel_sigma took          ",f10.4," s")') t9-t8
     !$    write(*,'("solve_tracers took           ",f10.4," s")') t10-t9
     !$ endif


  end select

!#ifdef USE_MPI
! if(mype==0) then
!  !write(93,*) t_1-t_0,t_2-t_1,t_3-t_2,t_4-t_3,t_5-t_4,t_6-t_5,t_7-t_6,t_8-t_7
! endif
!#endif

END SUBROUTINE oce_timestep
!==========================================================================
