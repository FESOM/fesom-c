subroutine initial_state

  use o_MESH
  use o_ARRAYS
  use o_PARAM
  use fv_ic

  use i_ARRAYS
  use i_PARAM

  use g_parsup
  use g_comm_auto

  implicit none
!  character(len = 256) :: fname,tmpstr ! file name

  integer       :: i,el, elnodes(4), n, m, j, jjj, nl, nz, c_riv, eledges(4), flag, elem
  real(kind=WP) :: x1, x2, y1, y2, amp, a, den_0
  real(kind=WP) :: dx, dy, period, hu, x, y
  real(kind=WP) :: pr, tt, ss, pp, d
  real(kind=WP) :: de, ga, ep, fe
  real(kind=WP), external :: theta
  real(kind=WP), allocatable, dimension(:,:)  ::   aux
  real(kind=WP), allocatable, dimension(:)  ::   aux_buff

  real, allocatable, dimension(:,:)  ::  S_t, T_t
  real(kind=WP) :: dmean

  integer :: node_size
  real(kind=WP) :: mx_eta, mn_eta, mx_dep, mn_dep
  real(kind=WP) :: mx_t, mn_t, mx_s, mn_s
  integer :: ierror

  integer :: ttldim

  real(kind=WP) :: all_elem_area(elem2D),tmp

  real(kind=WP), allocatable, dimension(:) :: Rj_buff,aj_buff, bj_buff
  real(kind=WP)                            :: Rval, aval, bval
  integer, allocatable, dimension(:) :: mapping
  integer                            :: iofs
  logical ::  there
  real(kind=WP), allocatable, dimension(:)  :: wef_buff
  real(kind=WP), allocatable, dimension(:)  :: mask_n_buff
  real(kind=WP), allocatable, dimension(:)  :: mask_n_ice_buff
  real(kind=WP) :: bcast_buff(12)
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
  node_size=myDim_nod2D+eDim_nod2D

  !====================================
  ! Calculating of reference density
  !====================================

#ifdef USE_MPI

  allocate(Rj_buff(nod2D), aj_buff(nod2D), bj_buff(nod2D))
  allocate(mapping(nod2D))
  Rj_buff(:)=0.0_WP; aj_buff(:)=0.0_WP; bj_buff(:)=0.0_WP
  mapping(:)=0

  do n=1, myDim_nod2D+eDim_nod2D
     iofs=myList_nod2D(n)
     mapping(iofs)=n
  end do
  all_elem_area = 0.0_WP
  call gather_elem(elem_area,all_elem_area)
  scale_area=sum(all_elem_area)/elem2D
  call MPI_BCast(scale_area, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
#else
  scale_area=sum(elem_area)/elem2D
#endif

  if (mype==0) write(*,*) 'scale_area = ',scale_area

  !====================================
  ! Calculating coeff. to determine proper eta (see fv_average_dynamic)
  ! in case of absence or presence tidal potential.
  ! Attention, in case of presence tidal potential, a_tp can be also calculating
  ! for each time step to take into account the influence of internal waves presence
  ! to tidal potential energy
  !=====================================

  allocate(a_tp(node_size,2))
  if (T_potential) then
     a_tp(:,1)=0.918_WP
     a_tp(:,2)=0.69_WP
  else
     a_tp(:,1)=1.0_WP
     a_tp(:,2)=0.0_WP
  endif
  do i=1,node_size
     if (index_nod2D(i) == 2) then
        a_tp(i,1)=1.0_WP
        a_tp(i,2)=0.0_WP
     endif
  enddo

  ! =========
  ! initialize with zeros
  ! =========

  eta_n=0.0_WP
  U_n_2D=0.0_WP

  eta_n_1=0.0_WP
  eta_n_2=0.0_WP
  U_n_1=0.0_WP
  U_n_2=0.0_WP

!SHTEST
! Constant wind
!!$print *,'CONSTANT WIND FIELD'
!!$taux(:) = -0.05_WP
!!$tauy(:) = 0.0_WP

  do n=1,node_size
     d = depth(n) + eta_n(n)
     !print *, depth(n), eta_n(n), d
     if (d < 0.0_WP) then
        eta_n(n) = -depth(n)
        eta_n_1(n) = -depth(n)
        eta_n_2(n) = -depth(n)
     endif
     !  print*, eta_n(n)
  enddo

  !+++++++++++++++++++++++++++++
  ! for compute mask WaD
  !     we need a etaAB
  !+++++++++++++++++++++++++++++

  de=0.614_WP
  ga=0.088_WP
  ep=0.013_WP
  fe=1.0_WP-de-ga-ep
  ! Interpolate AM4
  ssh_rhs = eta_n
  etaAB=de*ssh_rhs+fe*eta_n+ga*eta_n_1+ep*eta_n_2


#ifdef USE_MPI
  call MPI_REDUCE(maxval(eta_n), mx_eta, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
       0, MPI_COMM_FESOM_C, MPIerr)
  call MPI_REDUCE(minval(eta_n), mn_eta, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
       0, MPI_COMM_FESOM_C, MPIerr)
  call MPI_REDUCE(maxval(depth), mx_dep, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
       0, MPI_COMM_FESOM_C, MPIerr)
  call MPI_REDUCE(minval(depth), mn_dep, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
       0, MPI_COMM_FESOM_C, MPIerr)

  if (mype==0) then
     write(*,*) 'MAX_MIN_SSH_initial= ', mx_eta, mn_eta
     write(*,*) 'MAX_MIN_depth= ', mx_dep, mn_dep
  end if

#else

  write(*,*) 'MAX_MIN_SSH_initial= ', maxval(eta_n),minval(eta_n)
  write(*,*) 'MAX_MIN_depth= ', maxval(depth),minval(depth)

#endif

  ! AA initialzation for sediment
  do n=1,node_size
     dmean=max(Dmin,eta_n(n) + depth(n))
     con(n) = 0.0_WP/dmean
     con_bc(n) = con(n)
  enddo
  ! AA

  !VF initialization of sediment module and taking into account
  !presence of sediments into the density calculation

  if (comp_sediment) call comp_sediment_ini

  density_0=1000.0_WP

  !VF, calculate density using info about T_const, S_const and Concentration of sus. sediments
  tmp=0.0_WP
  if (comp_sediment) then
!SH skipped for now max/min CF must be properly computed!!
!     call densityJM(T_const, S_const, tmp, den_0,(maxval(CF)+minval(CF))*0.5_WP)
  else
     call densityJM(T_const, S_const, tmp, den_0)
  endif

  if (mype==0) write(*,*) 'den_0 = ', den_0
  density_0=density_0+den_0

  if (mype==0) write(*,*) 'density_0', density_0


  if (TF_presence) then
     if (mype==0) then
        open(50,file=trim(meshpath)//'m2.out', status='old')
        read(50,*) n
     end if

#ifdef USE_MPI
     call MPI_BCast(n, 1, MPI_INTEGER, 0, MPI_COMM_FESOM_C, ierror)
#endif
!write(*,*) "N nodes at OB: ",n

     allocate(ampt(n,12),fazt(n,12))
     ampt=0.0_WP
     fazt=0.0_WP
     allocate(amp_factor(12),phy_factor(12))
     amp_factor = 1.0_WP
     phy_factor = 0.0_WP


     if (mype==0) then
        allocate(aux(n,n_th*2))
        allocate(aux_buff(n_th*2))
        do i=1,n
           read(50,*) jjj,aux_buff
           aux(i,:)=aux_buff
           do j=2,n_th*2,2
              if (aux(i,j) < 0.0_WP) aux(i,j) = 360.0_WP + aux(i,j)
           end do
        enddo
        close(50)
        j=1
        do i=1,12
           if (a_th(i)==1) then
              ampt(:,i)=aux(:,j)
              fazt(:,i)=aux(:,j+1)*rad
              j=j+2
           endif
        end do
        write (*,*) 'max_amplitude', maxval(ampt)
        deallocate(aux)
        deallocate(aux_buff)

        INQUIRE( FILE='ampl_factor.out', EXIST=there )
        if (there) then
	    ! read nodal/satellite correction factors
            open(50,file='ampl_factor.out', status='old')
            !(/'M2 ', 'S2 ', 'N2 ', 'K2 ','K1 ', 'O1 ', 'P1 ', 'Q1 ', 'Mf ', 'Mm ', 'Ssa', 'M4 '/)
            do i = 1, 12
                read(50,*) amp_factor(i),phy_factor(i)
            enddo
            phy_factor = phy_factor*rad  !convert to radians ?
            close(50)
            write(*,*) " amplification factors from ampl_factor.out: "
            do i = 1, 12
                write(*,*) amp_factor(i),phy_factor(i)
            enddo
        else
            write(*,*) " !!!!!!!!!!!!!!!   WARNING  !!!!!!!!!!!!!!!!!!!!"
            write(*,*) " !!!!!  no file ampl_factor.out found      !!!!!"
            write(*,*) " !!!!!  but TF_presence = True, tides at OB!!!!!"
            write(*,*) " !!  All amplification factors set to :!!!!!!!!!"
            write(*,*) " !!!!!  1.0 for amplitudes, and 0,0 for phases!!"
            write(*,*) " !!!!!!!!!!!!!!!   WARNING  !!!!!!!!!!!!!!!!!!!!"
            amp_factor = 1.0_WP
            phy_factor = 0.0_WP
        end if
     end if

#ifdef USE_MPI
     do i=1,n
        bcast_buff=ampt(i,1:12)
        call MPI_BCast(bcast_buff, 12, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
        ampt(i,1:12)=bcast_buff
        bcast_buff=fazt(i,1:12)
        call MPI_BCast(bcast_buff, 12, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
        fazt(i,1:12)=bcast_buff
     end do
     call MPI_BCast(amp_factor, 12, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
     call MPI_BCast(phy_factor, 12, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)

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
!       write(50,*) myList_nod2D(my_in_obn_idx)
!       end do
!
!endif

#endif

  endif ! TF_presence

  if (use_wef) then
     allocate(wef_buff(nod2D))

     if (mype==0) then
        write(*,*) "Read wef"
        open(50,file='wef.out', status='old')
        do n=1,nod2D
            read(50,*) wef_buff(n)
        enddo
        close(50)
     endif
#ifdef USE_MPI
     call MPI_BCast(wef_buff, nod2D, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
     do n=1,nod2D
        if (mapping(n)>0) then
            wef(mapping(n))=wef_buff(n)
        end if
     end do
#else
     wef=wef_buff
#endif
  end if

! allocate buffer for masks read
  allocate(mask_n_ice_buff(nod2D))

!  Read Cd 2D  ======================================
  if (use_Cd_2d) then
     if (mype==0) then
        write(*,*) "Read Cd 2D"
        open(50,file='cd_2d.out', status='old')
        do n=1,nod2D
            read(50,*) mask_n_ice_buff(n)
        enddo
        close(50)
     endif
#ifdef USE_MPI
     call MPI_BCast(mask_n_ice_buff, nod2D, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
     do n=1,nod2D
        if (mapping(n)>0) then
            C_d_2d_n(mapping(n))=mask_n_ice_buff(n)
        end if
     end do
#else
     C_d_2d_n=mask_n_ice_buff
#endif
     if (mype==0) then
        write(*,*) "Cd max min: ",maxval(C_d_2d_n),minval(C_d_2d_n)
     end if
     do elem=1,myDim_elem2D
        elnodes = elem2D_nodes(:,elem)
        C_d_2d_e(elem) = sum(w_cv(1:4,elem)*C_d_2d_n(elnodes))
     enddo
    if (mype==0) then
        write(*,*) "Read Cd 2D: DONE"
    end if
  end if

!   Read Ice masks     ==========================================
!  if (use_ice) then
  if (mask_icedyn) then
     if (mype==0) then
        write(*,*) "Read mask for ice dynamics"
        open(50,file='mask_ice_adv.out', status='old')
        do n=1,nod2D
            read(50,*) mask_n_ice_buff(n)
        enddo
        close(50)
     endif
#ifdef USE_MPI
     call MPI_BCast(mask_n_ice_buff, nod2D, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
     do n=1,nod2D
        if (mapping(n)>0) then
            mask_ice_adv_n(mapping(n))=mask_n_ice_buff(n)
        end if
     end do
#else
     mask_ice_adv_n=mask_n_ice_buff
#endif
     do elem=1,myDim_elem2D
        elnodes = elem2D_nodes(:,elem)
        mask_ice_adv_e(elem) = sum(w_cv(1:4,elem)*mask_ice_adv_n(elnodes))
     enddo
    if (mype==0) then
        write(*,*) "Read mask for ice dynamics: DONE"
    end if
     if (mype==0) then
        write(*,*) "Read 2 mask for ice dynamics"
        open(50,file='mask_ice_adv_big.out', status='old')
        do n=1,nod2D
            read(50,*) mask_n_ice_buff(n)
        enddo
        close(50)
     endif
#ifdef USE_MPI
     call MPI_BCast(mask_n_ice_buff, nod2D, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
     do n=1,nod2D
        if (mapping(n)>0) then
            mask_ice_adv_big_n(mapping(n))=mask_n_ice_buff(n)
        end if
     end do
#else
     mask_ice_adv_big_n=mask_n_ice_buff
#endif
     do elem=1,myDim_elem2D
        elnodes = elem2D_nodes(:,elem)
        mask_ice_adv_big_e(elem) = sum(w_cv(1:4,elem)*mask_ice_adv_big_n(elnodes))
     enddo
    if (mype==0) then
        write(*,*) "Read 2 mask for ice dynamics: DONE"
    end if
  end if

  !
  ! Irradiance in the upper ocean

  !
  if (type_task>2) then

!SH OBACHT INITIALIZE JER
! NOT YET MPI



     if (.not.(jer_const)) then

        allocate(aj(node_size),bj(node_size),Rj(node_size))

#ifdef USE_MPI

        if (mype==0) then

           open(50,file=trim(meshpath)//'jer.out', status='old')

           do j=1,nod2D
              read(50,*) Rj_buff(j),aj_buff(j),bj_buff(j)
           enddo

        end if

        call MPI_BCast(Rj_buff(:), nod2D, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
        call MPI_BCast(aj_buff(:), nod2D, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
        call MPI_BCast(bj_buff(:), nod2D, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)

        do j=1,nod2D

           Rval=Rj_buff(n)
           aval=aj_buff(n)
           bval=bj_buff(n)

           if (mapping(j)>0) then
              Rj(mapping(j))=Rval
              aj(mapping(j))=aval
              bj(mapping(j))=bval
           end if
        end do

        if (mype==0) close(50)

#else
        open(50,file=trim(meshpath)//'jer.out', status='old')

        do j=1,nod2D
           read(50,*) Rj(j),aj(j),bj(j)
        enddo

        close(50)

#endif

     endif

     if (key_ic) then
        call ic_do  ! from fv_ic.f90
     else
        if (TEMP_ON) then
!not MPI
!           open(50,file=trim(meshpath)//trim(TITLE)//'_temp.out', status='old')
!           do n=1,nod2D
!            read(50,*) (TF(nz,n),nz=1,nsigma-1)
!           enddo
!           close(50)
                 !  else
           TF = T_const
           if (len(trim(title))==2 .AND. title=='LE') then
              do n=1,node_size
                 a = coord_nod2D(1,n)
                 if (a < 0.0_WP) then
                    TF(:,n) = 5.0_WP/0.28_WP
                 else
                    TF(:,n) = 10.0_WP/0.28_WP
                 endif
              enddo
           endif
        endif
        if (SALINITY_ON) then
!not MPI
!           open(50,file=trim(meshpath)//trim(TITLE)//'_sal.out', status='old')
!           do n=1,nod2D
!            read(50,*) (SF(nz,n),nz=1,nsigma-1)
!           enddo
!           close(50)
!        else
            SF = S_const
        endif
     endif
!IK     if (TEMP_ON) then
!IK        !aa experiment lock exchange ONLY
!IK        do n=1,node_size
!IK           a = coord_nod2D(1,n)
!IK           do nl=1,nsigma-1
!IK              if (a < 0.0_WP) then
!IK                 TF(nl,n) = 5.0_WP/0.28_WP
!IK              else
!IK                 !  if (a > 1.E-4*rad) then
!IK                 TF(nl,n) = 10.0_WP/0.28_WP
!IK                 !  else
!IK                 !                        TF(nl,n) = 0.5_WP*(5.0_WP + 10.0_WP)/0.28_WP
!IK              endif
!IK           enddo
!IK        enddo
!IK     endif
!IK     !aa END for lock exchange experiment

     !
     ! Climatology
     !
     !   convert in situ temperature into potential temperature
     !
     pr=0.0_WP

!     do n=1,node_size
!        do nz=1,nsigma-1
!           tt = TF(nz,n)
!           ss = SF(nz,n)
!           pp = Z(nz,n)
!           TF(nz,n)=theta(ss,tt,pp,pr)
!        enddo
!     enddo
     T_old= TF
     S_old= SF
     Tclim = TF
     Sclim = SF

#ifdef USE_MPI

     call MPI_REDUCE(maxval(TF), mx_t, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
          0, MPI_COMM_FESOM_C, MPIerr)
     call MPI_REDUCE(minval(TF), mn_t, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
          0, MPI_COMM_FESOM_C, MPIerr)
     call MPI_REDUCE(maxval(SF), mx_s, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
          0, MPI_COMM_FESOM_C, MPIerr)
     call MPI_REDUCE(minval(SF), mn_s, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
          0, MPI_COMM_FESOM_C, MPIerr)

     if (mype==0) then
        write(*,*)'T_max, T_min', mx_t, mn_t
        write(*,*)'S_max, S_min', mx_s, mn_s
     end if

#else

     write(*,*)'T_max, T_min', maxval(TF), minval(TF)
     write(*,*)'S_max, S_min', maxval(SF), minval(SF)

#endif
  endif
  !VF: Initialisation of C_d depends on user choice
  if (BFC==2) then
     DO elem=1,myDim_elem2D
        elnodes=elem2D_nodes(:,elem)
        dmean = max(Dmin,sum(w_cv(1:4,elem)*(eta_n(elnodes) + depth(elnodes))))
        C_d_el(elem)=1.0_WP/(log(dmean/(z0b_min+za))/cka)**2
     enddo
#ifdef USE_MPI
     call exchange_elem(C_d_el)
#endif
  else
     C_d_el=C_d
  endif

  if (mype==0) write(*,*)'finish initial_state'

  !a call solve_tracer_second_order_rec(SF, 's')
  !a end if

end subroutine initial_state

subroutine initial_state_wind

  use o_MESH
  use o_ARRAYS
  use o_PARAM
  !
  implicit none

  integer             :: el, q, elnodes(4), windi, mass
  real(kind=WP) :: y1, y, dst



  windi=0
  mass=0

  ! =========
  ! initialize with zeros
  ! =========
  eta_n=0.0_WP
  U_n_2D=0.0_WP
  if(Tstepping==1) then
  eta_n_1=0.0_WP
  eta_n_1=0.0_WP
  U_n_1=0.0_WP
  U_n_2=0.0_WP
  end if


  IF (windi==1) THEN
  ! Stress divided by density_0

  ! initialize wind
  y=minval(coord_nod2D(2,:))
  y1=maxval(coord_nod2D(2,:))
  DO el=1,elem2D
  elnodes=elem2D_nodes(:,el)
  q=4
  if(elnodes(1)==elnodes(4)) q=3
     taux(el)=0.0001_WP*cos(pi*(sum(coord_nod2d(2,elnodes(1:q)))/dble(q)-y)/(y1-y))
     tauy(el)=0.0_WP
  END DO
    ! Southern ocean forcing
  !DO el=1, elem2D
  !   elnodes=elem2D_nodes(:, el)
  !   q=4
  ! if(elnodes(1)==elnodes(4)) q=3
  !   y=sum(coord_nod2D(2,elnodes(1:q)))/dble(q)
  !   if(y<-32.*rad) then
  !  dst=1.0
  !   if(y<-64.*rad) dst=(y/rad+85)/21.0
  !   taux(el)=0.2*dst*dst*sin(pi*abs(y+32.0*rad)/(32.0*rad))
  !   end if
  !  tauy(el)=0.
  !END DO
  END IF
  if(mass==1) THEN
  ! provide a pattern with relaxation to prescribed elevation
  DO el=1,nod2D
  y=coord_nod2D(2,el)
  y1=coord_nod2D(1,el)
  dst=((y-50.0_WP*rad)**2+y1**2)/rad**2/4.0_WP   ! radius of 2 degrees
  if(dst<1) then
  relax_coef(el)=cos(pi*sqrt(dst)/2.0_WP)/3600.0_WP/24.0_WP/3.0_WP
  end if
  ! Set signal frequency
  !O_ext=2*pi/24.0_WP/3600.0_WP/365.0_WP/2.0_WP

  ENDDO
  END IF

end subroutine initial_state_wind

!===========================================================================

subroutine comp_sediment_ini

  use o_MESH
  use o_ARRAYS
  use o_PARAM

  use g_parsup

  implicit none
  integer       :: n,m,j,i
  real(kind=WP) :: b,s,dst, teta,w_s1,w_s2,w_s3,c

  integer  :: ierror

  CF=0.0_WP
  c_old=CF
  Cclim = CF

  if (sed_boun_flux) then
     if (mype==0) then

        open(50,file=trim(meshpath)//trim(TITLE)//'_sed_ob.out', status='old')
        read(50,*) n
        read(50,*) m

        close(50)
     end if

#ifdef USE_MPI
     call MPI_BCast(n, 1, MPI_INTEGER, 0, MPI_COMM_FESOM_C, ierror)
     call MPI_BCast(m, 1, MPI_INTEGER, 0, MPI_COMM_FESOM_C, ierror)
#endif

     allocate(cr_distr(nsigma-1,nobn,m),cr_distr2(nsigma-1,nobn),sed_vt2(m))

     if (mype==0) then
        open(50,file=trim(meshpath)//trim(TITLE)//'_sed_ob.out', status='old')
        read(50,*) n
        read(50,*) m

        do j=1,m
           read(50,*) sed_vt2(j)
           do i=1,n
              read(50,*)cr_distr(1:nsigma-1,i,j)
           enddo
        enddo
        close(50)

        cr_distr2(:,:)=cr_distr(:,:,1)
     end if

#ifdef USE_MPI
     call MPI_BCast(sed_vt2(:), m, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
     do i=1,nsigma-1
        call MPI_BCast(cr_distr2(i,1:nobn), nobn, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
     end do
     do i=1,nsigma-1
        do j=1,nobn
           call MPI_BCast(cr_distr(i,j,1:m), nobn, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
        end do
     end do
#endif

  end if


  !************************************************************************
  ! average settling velocity of the material available for transport (w_s)
  ! Simpson and Castelltort 2006, Coupled model of surface water flow,
  ! sediment transport and morphological evolution. Computers & Geosciences,
  ! 32, pp. 1600-1614.
  ! ************************************************************************
  b = 1.0_WP/3.0_WP
  s = plop_s/density_0 - 1.0_WP
  dst = d_m*(g*s/snu_kin**2)**b


  ! Criteria is from Miller, Cave and Komar, 1977, Threshold of sediment motion under unidirectional current + Shields criterium (dst>=4)

  if (dst<4) then
     teta=0.115_WP*(dst)**(-0.5_WP)
  else
     if (dst<=10) then
        teta=0.14_WP*(dst)**(-0.64_WP)
     else
        if (dst<=20) then
           teta=0.04_WP*(dst)**(-0.1_WP)
        else
           if (dst<=150) then
              teta=0.013_WP*(dst)**(0.29_WP)
           else
              teta= 0.055_WP
           endif
        endif
     endif
  endif

  c = 13.95_WP*snu_kin/d_m
  w_s1 = sqrt(c*c + 1.09_WP*d_m*g*s) - c
  w_s2 =snu_kin/d_m*((10.36_WP**2+1.049*teta**3)**(0.5)-1-.36_WP)
  w_s3=g*(plop_s-density_0)*d_m**2/(18.0_WP*snu_kin*density_0)

  if (mype==0) write(*,*) 'w_s1, w_s2, w_s3,...', w_s1,w_s2,w_s3

end subroutine comp_sediment_ini

