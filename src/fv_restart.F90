subroutine write_restart(n)
    use o_MESH
    use o_PARAM
    use o_ARRAYS
    use g_parsup
    use g_comm_auto

    implicit none

    integer  :: n
    character(len = 256) :: fname,tmp ! file name
    ! variables of some preliminary fields used for restart
    real(kind=WP), allocatable, dimension(:) :: array1d_glob
    real(kind=WP), allocatable, dimension(:,:) :: array2d_glob

#ifdef USE_MPI
    if (mype==0) then
    !write(*,*) 'Write restart', n
      write(tmp,*) n
      write(fname,*) 'restart_',trim(ADJUSTL(trim(tmp))),'.out'
      fname=ADJUSTL(trim(fname))
      write(*,*) 'Write restart file: ', fname
      open(25,file=fname, form='unformatted')
    endif
    !
    ! barotropic part
    !
!array1d_glob
    if (mype==0) then
        allocate(array1d_glob(nod2D))
    endif
    call gather_nod(eta_n, array1d_glob)
    if (mype==0) then
       write(25) array1d_glob
    endif
    call gather_nod(eta_n_1, array1d_glob)
    if (mype==0) then
       write(25) array1d_glob
    endif
    call gather_nod(eta_n_2, array1d_glob)
    if (mype==0) then
       write(25) array1d_glob
    endif
    if (mype==0) then
        deallocate(array1d_glob)
    endif
!array2d_glob
    if (mype==0) then
        allocate(array2d_glob(2,elem2D))
    endif
    call gather_elem(U_n_2D, array2d_glob)
    if (mype==0) then
       write(25) array2d_glob
    endif
    call gather_elem(U_n_1, array2d_glob)
    if (mype==0) then
       write(25) array2d_glob
    endif
    call gather_elem(U_n_2, array2d_glob)
    if (mype==0) then
       write(25) array2d_glob
    endif
    if (mype==0) then
        deallocate(array2d_glob)
    endif
    !
    ! 3D velocity
    !
    if (mype==0) then
        allocate(array2d_glob(nsigma-1,elem2D))
    endif
    call gather_elem(U_n, array2d_glob)
    if (mype==0) then
       write(25) array2d_glob
    endif
    call gather_elem(V_n, array2d_glob)
    if (mype==0) then
       write(25) array2d_glob
    endif
    call gather_elem(U_rhsAB, array2d_glob)
    if (mype==0) then
       write(25) array2d_glob
    endif
    call gather_elem(V_rhsAB, array2d_glob)
    if (mype==0) then
       write(25) array2d_glob
    endif
    if (mype==0) then
        deallocate(array2d_glob)
        allocate(array2d_glob(nsigma,nod2D))
    endif
    call gather_nod(Wvel, array2d_glob)
    if (mype==0) then
       write(25) array2d_glob
       deallocate(array2d_glob)
    endif

    !
    ! tracer
    if (mype==0) then
       allocate(array2d_glob(nsigma-1,nod2D))
    endif
    if (allocated(TF)) then
        call gather_nod(TF, array2d_glob)
        if (mype==0) then
            write(25) array2d_glob
        endif
    end if
    if (allocated(SF)) then
        call gather_nod(SF, array2d_glob)
        if (mype==0) then
            write(25) array2d_glob
        endif
    end if
    if (allocated(TF)) then
        call gather_nod(T_old, array2d_glob)
        if (mype==0) then
            write(25) array2d_glob
        endif
    end if
    if (allocated(SF)) then
        call gather_nod(S_old, array2d_glob)
        if (mype==0) then
            write(25) array2d_glob
        endif
    end if
    if (mype==0) then
        deallocate(array2d_glob)
    endif
    !
    ! vertical mixing
    !
    if (ver_mix == 2) then
        if (mype==0) then
            allocate(array2d_glob(nsigma,nod2D))
        endif
        call gather_nod(tke, array2d_glob)
        if (mype==0) then
            write(25) array2d_glob
        endif
        call gather_nod(Av_node, array2d_glob)
        if (mype==0) then
            write(25) array2d_glob
        endif
        call gather_nod(teps, array2d_glob)
        if (mype==0) then
            write(25) array2d_glob
        endif
        call gather_nod(Kv, array2d_glob)
        if (mype==0) then
            write(25) array2d_glob
        endif
        if (mype==0) then
            deallocate(array2d_glob)
        endif
    endif

    if (ver_mix == 3) then
        if (mype==0) then
            allocate(array2d_glob(nsigma,nod2D))
        endif
        call gather_nod(bt, array2d_glob)
        if (mype==0) then
            write(25) array2d_glob
        endif
        call gather_nod(snu, array2d_glob)
        if (mype==0) then
            write(25) array2d_glob
        endif
        call gather_nod(KV, array2d_glob)
        if (mype==0) then
            write(25) array2d_glob
        endif
        if (mype==0) then
            deallocate(array2d_glob)
        endif
    endif
    !
    ! time interval (+1)
    !
    if (mype==0) then
       write(25) n
       write(25) time_jd, time_jd0, dt, dt_2D
    !
    ! T_counter (for tidal output)
    !
       write(25) T_counter
    endif
    !
    ! tracer Clim
    if (mype==0) then
       allocate(array2d_glob(nsigma-1,nod2D))
    endif
    if (allocated(TF)) then
        call gather_nod(Tclim, array2d_glob)
        if (mype==0) then
            write(25) array2d_glob
        endif
    end if
    if (allocated(SF)) then
        call gather_nod(Sclim, array2d_glob)
        if (mype==0) then
            write(25) array2d_glob
        endif
    end if
    if (mype==0) then
        deallocate(array2d_glob)
    endif

    if (mype==0) then
       close(25)
    endif

#else
    write(*,*) 'Write restart', n
    !
    ! Write restart
    write(tmp,*) n
    write(fname,*) 'restart_',trim(ADJUSTL(trim(tmp))),'.out'
    fname=ADJUSTL(trim(fname))
    write(*,*) 'Write restart file: ', fname
    open(25,file=fname, form='unformatted')
    !
    ! barotropic part
    !
     write(25) eta_n, eta_n_1, eta_n_2, U_n_2D, U_n_1, U_n_2
    !
    ! 3D velocity
    !
     write(25) U_n, V_n, U_rhsAB, V_rhsAB, Wvel
    !
    ! tracer
    !
     write(25) TF, SF, T_old, S_old
    !
    ! vertical mixing
    !
     if (ver_mix == 2) write(25) tke, Av_node, teps, Kv
     if (ver_mix == 3) write(25) bt, snu, KV
    !
    ! time interval (+1)
    !
     write(25) n
     write(25) time_jd, time_jd0, dt, dt_2D
    !
    ! T_counter (for tidal output)
    !
     write(25) T_counter
     write(25) Tclim, Sclim
     close(25)
#endif
end subroutine write_restart

!=======================================   READ ==========================
subroutine read_restart
    use o_MESH
    use o_PARAM
    use o_ARRAYS
    use g_parsup
    use g_comm_auto

implicit none
    real(kind=WP) :: tmp_r
    real(kind=WP), allocatable, dimension(:) :: eta_n_g, eta_n_1_g, eta_n_2_g
    real(kind=WP), allocatable, dimension(:,:) :: U_n_2D_g, U_n_1_g, U_n_2_g
    real(kind=WP), allocatable, dimension(:,:) :: U_n_g, V_n_g, U_rhsAB_g, V_rhsAB_g
    real(kind=WP), allocatable, dimension(:,:) :: Wvel_g
    real(kind=WP), allocatable, dimension(:,:) :: TF_g, SF_g, T_old_g, S_old_g
    real(kind=WP), allocatable, dimension(:,:) :: tke_g, Av_node_g, teps_g, Kv_g
    real(kind=WP), allocatable, dimension(:,:) :: bt_g, snu_g

#ifdef USE_MPI
    if (mype==0) write(*,*) 'Reading restart'

  ! Read restart
    if (mype==0) then
        allocate(eta_n_g(nod2D),eta_n_1_g(nod2D),eta_n_2_g(nod2D))
        allocate(U_n_2D_g(2,elem2D),U_n_1_g(2,elem2D),U_n_2_g(2,elem2D))

        open(25,file='restart.out', form='unformatted')
        read(25) eta_n_g, eta_n_1_g, eta_n_2_g, U_n_2D_g, U_n_1_g, U_n_2_g
    endif

    if (allocated(eta_n)) call broadcast_nod(eta_n, eta_n_g)
    if (allocated(eta_n_1)) call broadcast_nod(eta_n_1, eta_n_1_g)
    if (allocated(eta_n_2)) call broadcast_nod(eta_n_2, eta_n_2_g)
    if (allocated(U_n_2D)) call broadcast_elem(U_n_2D, U_n_2D_g)
    if (allocated(U_n_1)) call broadcast_elem(U_n_1, U_n_1_g)
    if (allocated(U_n_2)) call broadcast_elem(U_n_2, U_n_2_g)
    if (mype==0) then
        deallocate(eta_n_g,eta_n_1_g,eta_n_2_g)
        deallocate(U_n_2D_g,U_n_1_g,U_n_2_g)
    endif

    if (mype==0) then
        allocate(U_n_g(nsigma-1,elem2D), V_n_g(nsigma-1,elem2D), U_rhsAB_g(nsigma-1,elem2D), V_rhsAB_g(nsigma-1,elem2D))
        allocate(Wvel_g(nsigma,nod2D))
        read(25) U_n_g, V_n_g, U_rhsAB_g, V_rhsAB_g, Wvel_g
    endif

    if (allocated(Wvel)) call broadcast_nod(Wvel, Wvel_g)
    if (allocated(U_n)) call broadcast_elem(U_n, U_n_g)
    if (allocated(V_n)) call broadcast_elem(V_n, V_n_g)
    if (allocated(U_rhsAB)) call broadcast_elem(U_rhsAB, U_rhsAB_g)
    if (allocated(V_rhsAB)) call broadcast_elem(V_rhsAB, V_rhsAB_g)
    if (mype==0) then
        deallocate(U_n_g, V_n_g, U_rhsAB_g, V_rhsAB_g)
        deallocate(Wvel_g)
    endif



    if (mype==0) then
       allocate(TF_g(nsigma-1,nod2D),SF_g(nsigma-1,nod2D),T_old_g(nsigma-1,nod2D),S_old_g(nsigma-1,nod2D))
       read(25) TF_g, SF_g, T_old_g, S_old_g
    endif

    if (allocated(TF)) call broadcast_nod(TF, TF_g)
    if (allocated(SF)) call broadcast_nod(SF, SF_g)
    if (allocated(T_old)) call broadcast_nod(T_old, T_old_g)
    if (allocated(S_old)) call broadcast_nod(S_old, S_old_g)
    if (mype==0) then
        deallocate(TF_g,SF_g,T_old_g,S_old_g)
    endif

    if (ver_mix == 2) then
       if (mype==0) then
          allocate(tke_g(nsigma,nod2D), Av_node_g(nsigma,elem2D), teps_g(nsigma,nod2D), Kv_g(nsigma,nod2D))
          read(25) tke_g, Av_node_g, teps_g, Kv_g
       endif
       if (allocated(tke)) call broadcast_nod(tke, tke_g)
       if (allocated(teps)) call broadcast_nod(teps, teps_g)
       if (allocated(Kv)) call broadcast_nod(Kv, Kv_g)
       if (allocated(Av_node)) call broadcast_elem(Av_node, Av_node_g)
       if (mype==0) then
          deallocate(tke_g, Av_node_g, teps_g, Kv_g)
       endif
    endif
    if (ver_mix == 3) then
       if (mype==0) then
          allocate(bt_g(nsigma,nod2D),snu_g(nsigma,nod2D),KV_g(nsigma,nod2D))
          read(25) bt_g, snu_g, KV_g
       endif
       if (allocated(bt)) call broadcast_nod(bt, bt_g)
       if (allocated(snu)) call broadcast_nod(snu, snu_g)
       if (allocated(KV)) call broadcast_nod(KV, KV_g)
       if (mype==0) then
           deallocate(bt_g,snu_g,KV_g)
       endif
    endif
!
! time interval (+1)
!
    if (mype==0) then
       read(25) ini_time
	   !VF The nsteps changes (decreasing) if you use hot_start
       nsteps=nsteps-ini_time
       ini_time=ini_time+1
       read(25) time_jd, tmp_r, dt_restart, dt_2D_restart !time_jd0 -> tmp_r
       read(25) T_counter
       close(25)

    endif

    call MPI_BCast(time_jd, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, MPIerr)

!    time_jd0=time_jd
    if (mype==0) write(*,*) 'ini_time_jd', time_jd0
    lfirst=.false.
#else

    write(*,*) 'Reading restart'
    !
    ! Read restart
    open(25,file='restart.out', form='unformatted')

    read(25) eta_n, eta_n_1, eta_n_2, U_n_2D, U_n_1, U_n_2
    read(25) U_n, V_n, U_rhsAB, V_rhsAB, Wvel
    read(25) TF, SF, T_old, S_old
    if (ver_mix == 2) read(25) tke, Av_node, teps, Kv
    if (ver_mix == 3) read(25) bt, snu, Kv
    read(25) ini_time
    !VF The nsteps changes (decreasing) if you use hot_start
    nsteps=nsteps-ini_time
    ini_time=ini_time+1
    read(25) time_jd, tmp_r, dt_restart, dt_2D_restart !time_jd0 -> tmp_r
    read(25) T_counter
    close(25)

    write(*,*) 'T_counter', T_counter
    write(*,*) 'ini_time', ini_time
!    time_jd0=time_jd
    write(*,*) 'ini_time_jd', time_jd0
    lfirst=.false.
#endif
end subroutine read_restart


subroutine write_restart_separate(n)
    use o_MESH
    use o_PARAM
    use o_ARRAYS
    use g_parsup
    use g_comm_auto
    use i_ARRAYS
    use i_PARAM

    implicit none

    integer  :: n
    character(len = 256) :: fname,tmp,tmp2 ! file name
    character(len=10)  :: dirname
    real(kind=WP) :: tmp_n(myDim_nod2D+eDim_nod2D)
#ifdef USE_MPI
    if (allocated(eta_n)) call exchange_nod(eta_n)
    if (allocated(eta_n_1)) call exchange_nod(eta_n_1)
    if (allocated(eta_n_2)) call exchange_nod(eta_n_2)
    if (allocated(U_n_2D)) call exchange_elem(U_n_2D)
    if (allocated(U_n_1)) call exchange_elem(U_n_1)
    if (allocated(U_n_2)) call exchange_elem(U_n_2)
    if (allocated(U_n)) call exchange_elem(U_n)
    if (allocated(V_n)) call exchange_elem(V_n)
    if (allocated(U_rhsAB)) call exchange_elem(U_rhsAB)
    if (allocated(Wvel)) call exchange_nod(Wvel)
    if (allocated(TF)) call exchange_nod(TF)
    if (allocated(SF)) call exchange_nod(SF)
    if (allocated(T_old)) call exchange_nod(T_old)
    if (allocated(S_old)) call exchange_nod(S_old)
    if (ver_mix == 2) then
       if (allocated(tke)) call exchange_nod(tke)
       if (allocated(Av_node)) call exchange_nod(Av_node)
       if (allocated(teps)) call exchange_nod(teps)
       if (allocated(Kv)) call exchange_nod(Kv)
    endif
    if (ver_mix == 3) then
       if (allocated(bt)) call exchange_nod(bt)
       if (allocated(snu)) call exchange_nod(snu)
       if (allocated(Kv)) call exchange_nod(Kv)
    endif
    if (allocated(Tclim)) call exchange_nod(Tclim)
    if (allocated(Sclim)) call exchange_nod(Sclim)

    if (use_ice) then
        call exchange_elem(U_n_ice)
        call exchange_elem(U_n_1ice)
        call exchange_elem(U_n_2ice)
        call exchange_elem(UiceAB)
        call exchange_nod(a_ice)
        call exchange_nod(m_ice)
        call exchange_nod(m_snow)
        call exchange_nod(ice_temp)
        call exchange_nod(t_skin)
        !call exchange_nod(cHrhsi)
!        do i=1,3
!            tmp_n = cHrhsi(i,:)
!            call exchange_nod(tmp_n)
!            cHrhsi(i,:) = tmp_n
!        enddo
    endif
#endif

    ! Write restart
    if (mype==0) write(*,*) ' Write Restart at step: ', n
    write(tmp,*) mype
    write (dirname, '(I10.10)')  n
!    if (mype==0) call execute_command_line('mkdir -p restarts/' // adjustl(trim( dirname ) ) )
    write(fname,*) 'restarts/',trim(ADJUSTL(trim(dirname))),'_restart_',trim(ADJUSTL(trim(tmp))),'.out'
    fname=ADJUSTL(trim(fname))
    !write(*,*) 'Write restart file: ', fname
    open(25,file=fname, form='unformatted')
    !
    ! time interval (+1)
    !
    write(25) n
    write(25) time_jd, time_jd0, dt, dt_2D
    !
    ! barotropic part
    !
    write(25) eta_n, eta_n_1, eta_n_2, U_n_2D, U_n_1, U_n_2

    if (type_task>1) then
    !
    ! 3D velocity
    !
        write(25) U_n, V_n, U_rhsAB, V_rhsAB, Wvel
    !
    ! vertical mixing
    !
        if (ver_mix == 2) write(25) tke, Av_node, teps, Kv
        if (ver_mix == 3) write(25) bt, snu, KV

    endif
    !
    ! tracer
    !
    if (type_task>2) then
        write(25) TF, SF, T_old, S_old
        write(25) Tclim, Sclim
    endif

    if (use_ice) then
        write(25) U_n_ice, U_n_1ice,U_n_2ice,UiceAB
        write(25) a_ice, m_ice, m_snow, ice_temp, t_skin
        write(25) cHrhsi
    endif
    !
    ! T_counter (for tidal output)
    !
    ! write(25) T_counter

    close(25)

end subroutine write_restart_separate

subroutine read_restart_separate
    use o_MESH
    use o_PARAM
    use o_ARRAYS
    use g_parsup
    use g_comm_auto
    use i_ARRAYS
    use i_PARAM

implicit none
    real(kind=WP) :: tmp_r
    character(len = 256) :: fname,tmp,tmp2 ! file name
    character(len=10)  :: dirname
    integer :: i
    real(kind=WP) :: tmp_n(myDim_nod2D+eDim_nod2D) !temporary array for exchange between cpu

    if (mype==0) write(*,*) 'Reading restart'

    write(tmp,*) mype
    write(fname,*) 'restarts/restart_',trim(ADJUSTL(trim(tmp))),'.out'
    fname=ADJUSTL(trim(fname))

    ! Read restart
    open(25,file=fname, form='unformatted')

    read(25) ini_time  ! old n_dt in  restart , ridiculous
    nsteps=nsteps-ini_time
    ini_time=ini_time+1
    read(25) time_jd, tmp_r, dt_restart, dt_2D_restart !time_jd0 -> tmp_r

    read(25) eta_n, eta_n_1, eta_n_2, U_n_2D, U_n_1, U_n_2
    call exchange_nod(eta_n)
    call exchange_nod(eta_n_1)
    call exchange_nod(eta_n_2)
    call exchange_elem(U_n_2D)
    call exchange_elem(U_n_1)
    call exchange_elem(U_n_2)

    if (type_task>1) then
        read(25) U_n, V_n, U_rhsAB, V_rhsAB, Wvel
        call exchange_elem(U_n)
        call exchange_elem(V_n)
        call exchange_elem(U_rhsAB)
        call exchange_nod(Wvel)
        if (ver_mix == 2) then
           read(25) tke, Av_node, teps, Kv
           call exchange_nod(tke)
           call exchange_nod(Av_node)
           call exchange_nod(teps)
           call exchange_nod(Kv)
        endif
        if (ver_mix == 3) then
           read(25) bt, snu, Kv
           call exchange_nod(bt)
           call exchange_nod(snu)
           call exchange_nod(Kv)
        endif
    endif

    if (type_task>2) then
        read(25) TF, SF, T_old, S_old
        call exchange_nod(TF)
        call exchange_nod(SF)
        call exchange_nod(T_old)
        call exchange_nod(S_old)
        read(25) Tclim, Sclim
        call exchange_nod(Tclim)
        call exchange_nod(Sclim)
    endif

    if (use_ice) then
        read(25) U_n_ice, U_n_1ice,U_n_2ice,UiceAB
        read(25) a_ice, m_ice, m_snow!, ice_temp, t_skin
!        read(25) cHrhsi
!        call exchange_elem(U_n_ice)
!        call exchange_elem(U_n_1ice)
!        call exchange_elem(U_n_2ice)
!        call exchange_elem(UiceAB)
!        call exchange_nod(a_ice)
!        call exchange_nod(m_ice)
!        call exchange_nod(m_snow)
!        call exchange_nod(ice_temp)
!        call exchange_nod(t_skin)
!        call exchange_nod(cHrhsi)
!        do i=1,3
!            tmp_n = cHrhsi(i,:)
!            call exchange_nod(tmp_n)
!            cHrhsi(i,:) = tmp_n
!        enddo
    endif
    !read(25) T_counter
    close(25)

!    write(*,*) 'T_counter', T_counter
!    write(*,*) 'ini_time', ini_time
!    time_jd0=time_jd
    if (mype==0) write(*,*) 'Restart reading finished, time_jd= ', time_jd, ',time_jd0=',time_jd0
    lfirst=.false.
end subroutine read_restart_separate
