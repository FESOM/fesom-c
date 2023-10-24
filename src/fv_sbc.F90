MODULE fv_sbc
   !!===========================================================================
   !! Ocean forcing:
   !!===========================================================================
   !! History: 0.1 ! 07/2015 I. Kuznetsov
   !!
   !! WARNING: in getcoeffld;nc_readTimeGrid
   !!   module will flip data in infiles for lat from -90 to 90 (NCEP-DOE Reanalysis 2 standart)
   !!   if you are going to use other coordinates in input files, please rewrite getcoeffld and nc_readTimeGrid functions.
   !!! WARNING : for now reading only files were time counts in hours from 1800 year; nc_readTimeGrid ; done FOR NCEP data
   !!! WARNING : move time forward for dt/2 , last is a last+dt from last -1 ; nc_readTimeGrid ; done FOR NCEP data

   !! Description:
   !!   read and interpolate atmpospheric forcing on model grid,
   !!     or use constants from namelist each time step
   !!
   !!   first initialization before first time step
   !!     model will read namelist, made some checks and prepare first interpolation coeficients
   !!   during model time steping, each time step sb_do is calling
   !!     during run model check if there is a time to read new data and construct new interpolation coefients
   !!     and interpolate data for new time step
   !!
   !!   taux and tuay defined and allocated outside of this module, but will be changed in this module
   !!   qns - Downward Non Solar heat flux over the ocean defined here and changed here, to use it outside:
   !!   USE fv_sbc
   !!   emp - evaporation minus precipitation defined here and changed here, to use it outside:
   !!   USE fv_sbc
   !!
   !! NetCDF:
   !!   we assume that all NetCDF files have identical grid and time variable
   !!   nm_sbc=2  nm_sbc_ftype=2 nm_tauwind=2
   !!
   !! public:
   !!   sbc_ini  -- inizialization atmpospheric forcing
   !!   sbc_do   -- provide a sbc (surface boundary conditions) each time step
   !!
   USE o_ARRAYS
   USE o_MESH
   USE o_PARAM


   IMPLICIT NONE

   include 'netcdf.inc'

   public  sbc_ini  ! routine called before 1st time step (open files, read namelist,...)
   public  sbc_do   ! routine called each time step to provide a sbc fileds (wind,...)
   public  sbc_end  ! routine called after last time step
   public  julday   ! get julian day from date

   private



   ! namelists
   integer, save  :: nm_sbc_unit     = 101       ! unit to open namelist file
  !============== namelistatmdata variables ================
   integer, save  :: nm_sbc       = 1        ! data  1= constant, 2=from file
!   integer, save  :: nm_sbc_ftype = 1        ! input file type if nm_sbc=2 : 1 = ASCII, 2 = netcdf   Deprecated
   integer, save  :: nm_tauwind   = 2        ! 1 = wind stress, 2 = wind 10 m  .
   character(len=256), save   :: nm_xwind_file = 'xwind.dat' ! name of file with winds/stress, if netcdf file then provide only name from "nameyyyy.nc" yyyy.nc will be added by model
   character(len=256), save   :: nm_ywind_file = 'ywind.dat' ! name of file with winds/stress, if netcdf file then provide only name from "nameyyyy.nc" yyyy.nc will be added by model
   character(len=256), save   :: nm_humi_file  = 'humidity.dat' ! name of file with humidity,  if netcdf file then provide only name from "nameyyyy.nc" yyyy.nc will be added by model
   character(len=256), save   :: nm_qsr_file   = 'qsr.dat'   ! name of file with solar heat,   if netcdf file then provide only name from "nameyyyy.nc" yyyy.nc will be added by model
   character(len=256), save   :: nm_qlw_file   = 'qlw.dat'   ! name of file with Long wave,    if netcdf file then provide only name from "nameyyyy.nc" yyyy.nc will be added by model
   character(len=256), save   :: nm_tair_file  = 'tair.dat'  ! name of file with 2m air temperature, if netcdf file then provide only name from "nameyyyy.nc" yyyy.nc will be added by model
   character(len=256), save   :: nm_prec_file  = 'prec.dat'  ! name of file with total precipitation, if netcdf file then provide only name from "nameyyyy.nc" yyyy.nc will be added by model
   character(len=256), save   :: nm_mslp_file  = 'mslp.dat'  ! name of file with mean sea level pressure, if netcdf file then provide only name from "nameyyyy.nc" yyyy.nc will be added by model
   character(len=256), save   :: nm_cloud_file  = 'cloud.dat'  ! name of file with clouds, if netcdf file then provide only name from "nameyyyy.nc" yyyy.nc will be added by model

   character(len=34), save   :: nm_xwind_var = 'uwnd' ! name of variable in file with wind
   character(len=34), save   :: nm_ywind_var = 'vwnd' ! name of variable in file with wind
   character(len=34), save   :: nm_humi_var  = 'shum' ! name of variable in file with humidity
   character(len=34), save   :: nm_qsr_var   = 'dswrf'! name of variable in file with solar heat
   character(len=34), save   :: nm_qlw_var   = 'dlwrf'! name of variable in file with Long wave
   character(len=34), save   :: nm_tair_var  = 'air'  ! name of variable in file with 2m air temperature
   character(len=34), save   :: nm_prec_var  = 'prate'! name of variable in file with total precipitation
   character(len=34), save   :: nm_mslp_var  = 'mslp' ! name of variable in file with mean sea level pressure
   character(len=34), save   :: nm_cloud_var = 'cloud'! name of variable in file with clouds

   real(wp), save :: nm_xwind0 = 0._wp    ! constant 10m. wind value in i-direction/stress, if nm_sbc=1
   real(wp), save :: nm_ywind0 = 0._wp    ! constant 10m. wind value in j-direction, if nm_sbc=1
   real(wp), save :: nm_humi0  = 0._wp    ! constant humidity
   real(wp), save :: nm_qsr0   = 0._wp    ! constant solar heat
   real(wp), save :: nm_qlw0   = 0._wp    ! constant Long wave
   real(wp), save :: nm_tair   = 0._wp    ! constant 2m air temperature
   real(wp), save :: nm_prec   = 0._wp    ! constant total precipitation
   real(wp), save :: nm_mslp   = 0._wp    ! constant mean sea level pressure
   real(wp), save :: nm_cloud  = 0._wp    ! constant clouds,

   real(wp),public, save :: depth_swr = 5._wp    ! depth of swr penetration
   real(wp), save,public :: nm_prec_coef = 1._wp ! precipitation will be devide by this constant (3600 for CoastDat, 1 for NCEP) (3600 - total precipitation per hour)
   integer , save :: nm_net_flux  = 0     ! constant for downward longwave heat over the ocean: 0 - downward, 1 - Net downward

!   integer, save  :: nm_calc_flux = 0 !  =0 use precalculated I0,... (based in NEMO subroutins)
   ! ========== netCDF time param
   integer, save :: nm_nc_iyear = 1948    ! initial year of time axis in netCDF (1948 like CoastDat,1800 NCEP)
   integer, save :: nm_nc_imm = 1         ! initial month of time axis in netCDF
   integer, save :: nm_nc_idd = 1         ! initial day of time axis in netCDF
   real, save :: nm_nc_secstep = 86400.0_WP ! time units coef (86400 CoastDat, 24 NCEP)
   integer,save            :: warn       ! warning switch node/element coordinate out of forcing bounds

   ! ========== interpolation coeficients
   integer,  allocatable, save, dimension(:)     :: bilin_indx_i ! indexs i for interpolation
   integer,  allocatable, save, dimension(:)     :: bilin_indx_j ! indexs j for interpolation

   real(wp), allocatable, save, dimension(:,:)   :: coef_b ! time inerp coef. b (x=a*t+b)
   real(wp), allocatable, save, dimension(:,:)   :: coef_a ! time inerp coef. a (x=a*t+b)

   real(wp), allocatable, save, dimension(:)   :: datawx_f ! wind X data from first time slice
   real(wp), allocatable, save, dimension(:)   :: datawx_s ! wind X data from second time slice
   real(wp), allocatable, save, dimension(:)   :: datawy_f ! wind Y data from first time slice
   real(wp), allocatable, save, dimension(:)   :: datawy_s ! wind Y data from second time slice

   logical, save :: one_field = .false.! only one field used for forcing

   real(wp), allocatable, save, dimension(:,:), public   :: atmdata ! atmosperic data for current time step


!============== NETCDF ==========================================

!   character(len=256) :: tmp_str

   integer, parameter :: i_totfl = 8 ! total number of fluxes
   integer, parameter :: i_xwind = 1 ! index of 10m wind velocity (x-component) [m/s]
   integer, parameter :: i_ywind = 2 ! index of 10m wind velocity (y-component) [m/s]
   integer, parameter :: i_humi  = 3 ! index of specific humidity               [kg/kg]
   integer, parameter :: i_qsr   = 4 ! index of solar heat                      [W/m2]
   integer, parameter :: i_qlw   = 5 ! index of Long wave                       [W/m2]
   integer, parameter :: i_tair  = 6 ! index of 2m air temperature              [degK]
   integer, parameter :: i_prec  = 7 ! index of total precipitation (rain+snow) [Kg/m^2/s]
   integer, parameter :: i_mslp  = 8 ! index of mean sea level pressure         [Pascals]
   integer, parameter :: i_cloud = 9 ! index of clouds         [0-1]

   logical, save, public :: key_snow = .false.! True if snow is in forcing (used in ice model)

   type, public ::   flfi_type    !flux file information
      character(len = 256) :: file_name ! file name
      character(len = 34)  :: var_name  ! variable name in the NetCDF file
   end type flfi_type

  type(flfi_type),save, dimension(i_totfl) :: sbc_flfi  !array for information about flux files

  ! arrays of time, lon and lat in INfiles
   real(wp), allocatable, save, dimension(:)  :: nc_lon
   real(wp), allocatable, save, dimension(:)  :: nc_lat
   real(wp), allocatable, save, dimension(:)  :: nc_time
  ! lenght of arrays in INfiles
   integer,save              :: nc_Nlon
   integer,save              :: nc_Nlat
   integer,save              :: nc_Ntime
   ! time index for NC time array
   integer,save              :: t_indx    ! now time index in nc_time array
   integer,save              :: t_indx_p1 ! now time index +1 in nc_time array

  ! flip latitude from infiles (for example  NCEP-DOE Reanalysis 2 standart)
  integer, save              :: flip_lat ! 1 if we need to flip
!============== NETCDF ==========================================


CONTAINS
   SUBROUTINE nc_readTimeGrid(flf)
   ! Read time array and grid from nc file
      IMPLICIT NONE

      type(flfi_type),intent(in) :: flf

      integer              :: iost !I/O status
      integer              :: ncid      ! netcdf file id
      integer              :: i
      ! ID dimensions and variables:
      integer              :: id_lon
      integer              :: id_lat
      integer              :: id_time
      integer              :: id_lond
      integer              :: id_latd
      integer              :: id_timed
!      integer              :: nf_dims(4) ! dimensions (temporal)
      integer              :: nf_start(4)
      integer              :: nf_edges(4)
      integer              :: zero_year,yyyy,mm,dd
      character(len = 256) :: att_string ! attribute
      integer              :: sbc_alloc                   !: allocation status



      !open file
      iost = nf_open(flf%file_name,NF_NOWRITE,ncid)
      call check_nferr(iost,flf%file_name)

      ! get dimensions
      iost = nf_inq_dimid(ncid, "latitude", id_latd)
      call check_nferr(iost,flf%file_name)
      iost = nf_inq_dimid(ncid, "longitude", id_lond)
      call check_nferr(iost,flf%file_name)
      iost = nf_inq_dimid(ncid, "time", id_timed)
      call check_nferr(iost,flf%file_name)

      ! get variable id
      iost = nf_inq_varid(ncid, "longitude", id_lon)
      call check_nferr(iost,flf%file_name)
      iost = nf_inq_varid(ncid, "latitude", id_lat)
      call check_nferr(iost,flf%file_name)
      iost = nf_inq_varid(ncid, "time", id_time)
      call check_nferr(iost,flf%file_name)
      !  get dimensions size
      iost = nf_inq_dimlen(ncid, id_latd, nc_Nlat)
      call check_nferr(iost,flf%file_name)
      iost = nf_inq_dimlen(ncid, id_lond, nc_Nlon)
      call check_nferr(iost,flf%file_name)
      iost = nf_inq_dimlen(ncid, id_timed, nc_Ntime)
      call check_nferr(iost,flf%file_name)

      ALLOCATE( nc_lon(nc_Nlon), nc_lat(nc_Nlat), nc_time(nc_Ntime),&
                &      STAT=sbc_alloc )
      if( sbc_alloc /= 0 )   STOP 'read_sbc: failed to allocate arrays'

   !read variables from file
   ! coordinates
      nf_start(1)=1
      nf_edges(1)=nc_Nlat
      iost = nf_get_vara_double(ncid, id_lat, nf_start, nf_edges, nc_lat)
      call check_nferr(iost,flf%file_name)
      nf_start(1)=1
      nf_edges(1)=nc_Nlon
      iost = nf_get_vara_double(ncid, id_lon, nf_start, nf_edges, nc_lon)
      call check_nferr(iost,flf%file_name)
      nf_start(1)=1
      nf_edges(1)=nc_Ntime
      iost = nf_get_vara_double(ncid, id_time, nf_start, nf_edges, nc_time)
      call check_nferr(iost,flf%file_name)
      iost = nf_close(ncid)
      call check_nferr(iost,flf%file_name)

      ! convert time to days
      nc_time = nc_time / nm_nc_secstep + julday(nm_nc_iyear,nm_nc_imm,nm_nc_idd) - 0.5_WP
      !flip lat and data in case of lat from -90 to 90
      flip_lat = 0
      if ( nc_Nlat > 1 ) then
         if ( nc_lat(1) > nc_lat(nc_Nlat) ) then
            flip_lat = 1
            nc_lat=nc_lat(nc_Nlat:1:-1)
            write(*,*) "fv_sbc: nc_readTimeGrid: FLIP lat and data while lat from -90 to 90"
         endif
      endif



   END SUBROUTINE nc_readTimeGrid

   SUBROUTINE nc_sbc_ini_fillnames(yyear)
      character(len=4),intent(in)   :: yyear
      character(len=27) :: filenameend

      !! ** Purpose : Fill names of sbc_flfi array (file names and variable names)

      write(filenameend ,*) yyear,'.nc'
      filenameend = trim(filenameend)
      filenameend='       '
      filenameend = trim(filenameend)
      !prepare proper nc file (add year and .nc to the end of the file name from namelist
      write(sbc_flfi(i_xwind)%file_name,*) trim(nm_xwind_file),filenameend
      write(sbc_flfi(i_ywind)%file_name,*) trim(nm_ywind_file),filenameend
      write(sbc_flfi(i_humi)%file_name, *) trim(nm_humi_file),filenameend
      write(sbc_flfi(i_qsr)%file_name, *) trim(nm_qsr_file),filenameend
      write(sbc_flfi(i_qlw)%file_name, *) trim(nm_qlw_file),filenameend
      write(sbc_flfi(i_tair)%file_name, *) trim(nm_tair_file),filenameend
      write(sbc_flfi(i_prec)%file_name, *) trim(nm_prec_file),filenameend
      write(sbc_flfi(i_mslp)%file_name, *) trim(nm_mslp_file),filenameend

!      if (nm_calc_flux==1) then
!         write(sbc_flfi(i_cloud)%file_name, *) trim(nm_cloud_file),filenameend
!      end if

      sbc_flfi(i_xwind)%file_name=ADJUSTL(trim(sbc_flfi(i_xwind)%file_name))
      sbc_flfi(i_ywind)%file_name=ADJUSTL(trim(sbc_flfi(i_ywind)%file_name))
      sbc_flfi(i_humi)%file_name=ADJUSTL(trim(sbc_flfi(i_humi)%file_name))
      sbc_flfi(i_qsr)%file_name=ADJUSTL(trim(sbc_flfi(i_qsr)%file_name))
      sbc_flfi(i_qlw)%file_name=ADJUSTL(trim(sbc_flfi(i_qlw)%file_name))
      sbc_flfi(i_tair)%file_name=ADJUSTL(trim(sbc_flfi(i_tair)%file_name))
      sbc_flfi(i_prec)%file_name=ADJUSTL(trim(sbc_flfi(i_prec)%file_name))
      sbc_flfi(i_mslp)%file_name=ADJUSTL(trim(sbc_flfi(i_mslp)%file_name))
!      if (nm_calc_flux==1) then
!         sbc_flfi(i_cloud)%file_name=ADJUSTL(trim(sbc_flfi(i_cloud)%file_name))
!      end if

      sbc_flfi(i_xwind)%var_name=ADJUSTL(trim(nm_xwind_var))
      sbc_flfi(i_ywind)%var_name=ADJUSTL(trim(nm_ywind_var))
      sbc_flfi(i_humi)%var_name=ADJUSTL(trim(nm_humi_var))
      sbc_flfi(i_qsr)%var_name=ADJUSTL(trim(nm_qsr_var))
      sbc_flfi(i_qlw)%var_name=ADJUSTL(trim(nm_qlw_var))
      sbc_flfi(i_tair)%var_name=ADJUSTL(trim(nm_tair_var))
      sbc_flfi(i_prec)%var_name=ADJUSTL(trim(nm_prec_var))
      sbc_flfi(i_mslp)%var_name=ADJUSTL(trim(nm_mslp_var))
!      if (nm_calc_flux==1) then
!         sbc_flfi(i_cloud)%var_name=ADJUSTL(trim(nm_cloud_var))
!      end if

   END SUBROUTINE nc_sbc_ini_fillnames

   SUBROUTINE nc_sbc_ini(idate,isec)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE nc_sbc_ini ***
      !!
      !! ** Purpose : initialization of ocean forcing from NETCDF file
      !!----------------------------------------------------------------------

      IMPLICIT NONE

      integer,intent(in) :: idate ! initialization date
      integer,intent(in) :: isec ! initialization seconds

      character(len=4)   :: yyear
      integer            :: yyyy,mm,dd

      integer            :: i
      integer            :: sbc_alloc

      integer            :: elnodes(4) !4 nodes from one element
      integer            :: numnodes   ! nu,ber of nodes in elem (3 for triangle, 4 for ... )
      real(wp)           :: x, y       ! coordinates of elements

      ! get ini year; Fill names of sbc_flfi
      i=int(time_jd0+0.5_WP)
      call calendar_date(i,yyyy,dd,mm)
      write(yyear,"(I4)") yyyy
      call nc_sbc_ini_fillnames(yyear)

      ! we assume that all NetCDF files have identical grid and time variable
      call nc_readTimeGrid(sbc_flfi(i_xwind))
      ! prepare nearest coordinates in INfile , save to bilin_indx_i/j
      do i = 1, nod2D
         x  = coord_nod2d(1,i)/rad
         y  = coord_nod2d(2,i)/rad

         ! find nearest
         if ( x < nc_lon(nc_Nlon) .and. x >= nc_lon(1) ) then
            call binarysearch(nc_Nlon, nc_lon, x, bilin_indx_i(i))
         else ! NO extrapolation in space
            if ( x < nc_lon(1) ) then
               bilin_indx_i(i)=-1
            else
               bilin_indx_i(i)=0
            end if
         end if
         if ( y < nc_lat(nc_Nlat) .and. y >= nc_lat(1) ) then
            call binarysearch(nc_Nlat, nc_lat, y, bilin_indx_j(i))
         else ! NO extrapolation in space
            if ( y < nc_lat(1) ) then
               bilin_indx_j(i)=-1
            else
               bilin_indx_j(i)=0
            end if
         end if
         if (bilin_indx_i(i) < 1 .or. bilin_indx_j(i) < 1) then
            WRITE(*,*) '     WARNING:  node/element coordinate out of forcing bounds,'
            WRITE(*,*) '        nearest value will be used as a constant field'
         end if
      end do
      ! get first coefficients for time interpolation on model grid for all data
      call getcoeffld(time_jd)
      ! interpolate in time

      call data_timeinterp(time_jd)
      call apply_atm_fluxes
   END SUBROUTINE nc_sbc_ini

   SUBROUTINE apply_atm_fluxes
      !!----------------------------------------------------------------------
      !! ** Purpose : Change model variables according to atm fluxes
      !! source of original code: NEMO 3.1.1 + NCAR
      !!----------------------------------------------------------------------
      IMPLICIT NONE

      integer             :: i
      integer             :: num, n
      real(wp)            :: rtmp    ! temporal real
      real(wp)            :: wndm    ! delta of wind module and ocean curent module
      real(wp)            :: wdx,wdy ! delta of wind x/y and ocean curent x/y
      real(wp)            :: q_sat   ! sea surface specific humidity         [kg/kg]
      real(wp), parameter :: rhoa = 1.3 !1.22_WP ! air density
      real(wp), parameter :: cpa  = 1005.0_WP         ! specific heat of air
      real(wp), parameter :: Lv   =    2.5e6_WP       ! latent heat of vaporization
      real(wp), parameter :: Stef =    5.67e-8_WP     ! Stefan Boltzmann constant
      real(wp), parameter :: albo =    0.1!0.066_WP       ! ocean albedo assumed to be contant
      real(wp)            :: zst     ! surface temperature in Kelvin
      real(wp)            :: emiss_wat = 0.97

      real(wp)           ::  &
         Cd,       &     ! transfer coefficient for momentum         (tau)
         Ch,       &     ! transfer coefficient for sensible heat (Q_sens)
         Ce,       &     ! transfert coefficient for evaporation   (Q_lat)
         t_zu,     &     ! air temp. shifted at zu                     [K]
         q_zu            ! spec. hum.  shifted at zu               [kg/kg]

      real(wp)           :: zevap, zqsb, zqla, zqlw

      do i = 1,nod2D
         windx(i) = atmdata(i_xwind,i)
         windy(i) = atmdata(i_ywind,i)
          !IF you going to use winds on nodes (nod2d) U_n and V_n should be changed to Unode and Vnode
         if (type_task>1) then
            wdx = atmdata(i_xwind,i) - Unode(1,i) ! wind from data - ocean current ( x direction)
            wdy = atmdata(i_ywind,i) - Vnode(1,i) ! wind from data - ocean current ( y direction)
         else
            ! type_task = 1 ;convert 2D velocity on element to nodes (wdx), than modify wdx by
            n = i
            num = nod_in_elem2D_num(n)
            wdx = sum(U_n_2D(1,nod_in_elem2D(1:num,n))*elem_area(nod_in_elem2D(1:num,n))) &
                                                  /sum(elem_area(nod_in_elem2D(1:num,n)))
            wdy = sum(U_n_2D(2,nod_in_elem2D(1:num,n))*elem_area(nod_in_elem2D(1:num,n))) &
                                                  /sum(elem_area(nod_in_elem2D(1:num,n)))
            wdx = atmdata(i_xwind,i) - wdx
            wdy = atmdata(i_ywind,i) - wdy
         endif

         wndm = SQRT( wdx * wdx + wdy * wdy )

         if (type_task>2) then
            zst = TF(1,i)+273.15_WP
         else
            zst = T_const+273.15_WP
         endif

         q_sat = 0.98_WP * 640380._WP / rhoa * EXP( -5107.4_WP / zst )

         call core_coeff_2z(2.0_wp, 10.0_wp, zst, atmdata(i_tair,i), &
                           q_sat, atmdata(i_humi,i), wndm, Cd, Ch, Ce, t_zu, q_zu)

         Ch_atm_oce_arr(i)=Ch
         Ce_atm_oce_arr(i)=Ce
         Cd_atm_oce_arr(i)=Cd

         rtmp = rhoa * wndm * Cd

         taux_node(i) = rtmp * wdx  ! new taux (stress along x)
         tauy_node(i) = rtmp * wdy  ! new tauy (stress along y)

         zevap = MAX( 0.e0_WP, rhoa    *Ce*( q_sat - q_zu ) * wndm )   ! Evaporation
         zqsb  = rhoa*cpa*Ch*( zst - t_zu ) * wndm     ! Sensible Heat

         if ( nm_net_flux == 0 ) then
            ! NCEP case
!            zqlw  = (  atmdata(i_qlw,i) - Stef * zst*zst*zst*zst  )    ! Long  Wave
            zqlw  = (  atmdata(i_qlw,i) - emiss_wat*Stef * zst*zst*zst*zst  )    ! Long  Wave
         else
            ! CoastDat case with Net flux
            zqlw  = atmdata(i_qlw,i)
         endif
         zqla  = Lv * zevap   ! Latent Heat
!!!================================ data for use in ocean part ================================
         ! Downward Non Solar heat flux over the ocean
         qns(i) = zqlw - zqsb - zqla
    !!!!WARNING  !!!!!!!!!!!!!!!!!!
    ! when we reach temperature less 0 we start to have a problems, why ?
    ! need to check fluxes in qns , could be due to dis-balance of zqlw
         if (.not. use_ice) then
             if (type_task>2) then
                 if (TF(1,i)<-2.0_WP) then
                    if (qns(i) < 0.0_WP) then
                        qns(i) = 0.0_WP
                    end if
                 end if
             end if
         endif
    !!!WARNING   !!!!!!!!!!!!!!!!!!
         ! zqlw   = downward longwave heat over the ocean
         ! - zqsb = downward sensible heat over the ocean
         ! - zqla = downward latent   heat over the ocean
         ! qns    = downward non solar heat over the ocean
         emp(i) = zevap - atmdata(i_prec,i)/nm_prec_coef
         qsr(i) = (1.0_wp - albo ) * atmdata(i_qsr,i)
         mslp(i)= atmdata(i_mslp,i)
!!!================================ data for use in ocean part ================================
      end do
      ! convert tau to elements
      do i = 1, elem2D
         taux(i) = sum(w_cv(1:4,i)*taux_node(elem2D_nodes(:,i)))
         tauy(i) = sum(w_cv(1:4,i)*tauy_node(elem2D_nodes(:,i)))
      end do

   END SUBROUTINE apply_atm_fluxes


   SUBROUTINE core_coeff_2z(zt, zu, sst, T_zt, q_sat, q_zt, dU, Cd, Ch, Ce, T_zu, q_zu)
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  core_coeff_2z  ***
      !!
      !! ** Purpose :   Computes turbulent transfert coefficients of surface
      !!                fluxes according to Large & Yeager (2004).
      !!
      !! ** Method  :   I N E R T I A L   D I S S I P A T I O N   M E T H O D
      !!      Momentum, Latent and sensible heat exchange coefficients
      !!      Caution: this procedure should only be used in cases when air
      !!      temperature (T_air) and air specific humidity (q_air) are at 2m
      !!      whereas wind (dU) is at 10m.
      !!
      !! References :   Large & Yeager, 2004 : ???
      !! code was adopted from NEMO 3.3.1
      !!----------------------------------------------------------------------
      IMPLICIT NONE

      real(wp)            :: dU10        ! dU                             [m/s]
      real(wp)            :: dT          ! air/sea temperature difference   [K]
      real(wp)            :: dq          ! air/sea humidity difference      [K]
      real(wp)            :: Cd_n10      ! 10m neutral drag coefficient
      real(wp)            :: Ce_n10      ! 10m neutral latent coefficient
      real(wp)            :: Ch_n10      ! 10m neutral sensible coefficient
      real(wp)            :: sqrt_Cd_n10 ! root square of Cd_n10
      real(wp)            :: sqrt_Cd     ! root square of Cd
      real(wp)            :: T_vpot      ! virtual potential temperature    [K]
      real(wp)            :: T_star      ! turbulent scale of tem. fluct.
      real(wp)            :: q_star      ! turbulent humidity of temp. fluct.
      real(wp)            :: U_star      ! turb. scale of velocity fluct.
      real(wp)            :: L           ! Monin-Obukov length              [m]
      real(wp)            :: zeta_u      ! stability parameter at height zu
      real(wp)            :: zeta_t      ! stability parameter at height zt
      real(wp)            :: U_n10       ! neutral wind velocity at 10m     [m]
      real(wp)            :: xlogt , xct , zpsi_hu , zpsi_ht , zpsi_m
      real(wp)            :: stab        ! 1st guess stability test integer
      !!
      real(wp), intent(in)   :: &
         zt,      &     ! height for T_zt and q_zt                   [m]
         zu             ! height for dU                              [m]
      real(wp), intent(in)   ::  &
         sst,      &     ! sea surface temperature              [Kelvin]
         T_zt,     &     ! potential air temperature            [Kelvin]
         q_sat,    &     ! sea surface specific humidity         [kg/kg]
         q_zt,     &     ! specific air humidity                 [kg/kg]
         dU              ! relative wind module |U(zu)-U(0)|       [m/s]
      real(wp), intent(out)  ::  &
         Cd,       &     ! transfer coefficient for momentum         (tau)
         Ch,       &     ! transfer coefficient for sensible heat (Q_sens)
         Ce,       &     ! transfert coefficient for evaporation   (Q_lat)
         T_zu,     &     ! air temp. shifted at zu                     [K]
         q_zu            ! spec. hum.  shifted at zu               [kg/kg]

      integer :: j_itt
      integer,  parameter :: nb_itt = 5   ! number of itterations
      real(wp), parameter ::                        &
         grav   = 9.8_WP,      &  ! gravity
         kappa  = 0.4_WP          ! von Karman's constant
      !!----------------------------------------------------------------------
      !!  * Start

      !! Initial air/sea differences
      dU10 = max(0.5_wp, dU)      !  we don't want to fall under 0.5 m/s
      dT = T_zt - sst
      dq = q_zt - q_sat

      !! Neutral Drag Coefficient :
      stab = 0.5_WP + sign(0.5_wp,dT)                 ! stab = 1  if dT > 0  -> STABLE
      Cd_n10  = 1E-3_WP*( 2.7_WP/dU10 + 0.142_WP + dU10/13.09_WP )
      sqrt_Cd_n10 = sqrt(Cd_n10)
      Ce_n10  = 1E-3_WP*( 34.6_WP * sqrt_Cd_n10 )
      Ch_n10  = 1E-3_WP*sqrt_Cd_n10*(18._WP*stab + 32.7_WP*(1._WP - stab))

      !! Initializing transf. coeff. with their first guess neutral equivalents :
      Cd = Cd_n10 ;  Ce = Ce_n10 ;  Ch = Ch_n10 ;  sqrt_Cd = sqrt(Cd)

      !! Initializing z_u values with z_t values :
      T_zu = T_zt ;  q_zu = q_zt

      !!  * Now starting iteration loop
      do j_itt=1, nb_itt
         dT = T_zu - sst ;  dq = q_zu - q_sat ! Updating air/sea differences
         T_vpot = T_zu*(1._WP + 0.608_WP*q_zu)    ! Updating virtual potential temperature at zu
         U_star = sqrt_Cd*dU10                ! Updating turbulent scales :   (L & Y eq. (7))
         T_star  = Ch/sqrt_Cd*dT              !
         q_star  = Ce/sqrt_Cd*dq              !
         !!
         L = (U_star*U_star) &                ! Estimate the Monin-Obukov length at height zu
              & / (kappa*grav/T_vpot*(T_star*(1.+0.608_WP*q_zu) + 0.608_WP*T_zu*q_star))
         !! Stability parameters :
         zeta_u  = zu/L  ;  zeta_u = sign( min(abs(zeta_u),10.0_WP), zeta_u )
         zeta_t  = zt/L  ;  zeta_t = sign( min(abs(zeta_t),10.0_WP), zeta_t )
         zpsi_hu = psi_h(zeta_u)
         zpsi_ht = psi_h(zeta_t)
         zpsi_m  = psi_m(zeta_u)
         !!
         !! Shifting the wind speed to 10m and neutral stability : (L & Y eq.(9a))
!        U_n10 = dU10/(1. + sqrt_Cd_n10/kappa*(log(zu/10.) - psi_m(zeta_u)))
         !   In very rare low-wind conditions, the old way of estimating the
         !   neutral wind speed at 10m leads to a negative value that causes the code
         !   to crash. To prevent this a threshold of 0.25m/s is now imposed.
         U_n10 = max(0.25_WP , dU10/(1._WP + sqrt_Cd_n10/kappa*(log(zu/10._WP) - zpsi_m)))
         !!
         !! Shifting temperature and humidity at zu :          (L & Y eq. (9b-9c))
!        T_zu = T_zt - T_star/kappa*(log(zt/zu) + psi_h(zeta_u) - psi_h(zeta_t))
         T_zu = T_zt - T_star/kappa*(log(zt/zu) + zpsi_hu - zpsi_ht)
!        q_zu = q_zt - q_star/kappa*(log(zt/zu) + psi_h(zeta_u) - psi_h(zeta_t))
         q_zu = q_zt - q_star/kappa*(log(zt/zu) + zpsi_hu - zpsi_ht)
         !!
         !! q_zu cannot have a negative value : forcing 0
         stab = 0.5_WP + sign(0.5_wp,q_zu) ;  q_zu = stab*q_zu
         !!
         !! Updating the neutral 10m transfer coefficients :
         Cd_n10  = 1E-3_WP * (2.7_WP/U_n10 + 0.142_WP + U_n10/13.09_WP)    ! L & Y eq. (6a)
         sqrt_Cd_n10 = sqrt(Cd_n10)
         Ce_n10  = 1E-3_WP * (34.6_WP * sqrt_Cd_n10)                 ! L & Y eq. (6b)
         stab    = 0.5_WP + sign(0.5_wp,zeta_u)
         Ch_n10  = 1E-3_WP*sqrt_Cd_n10*(18._WP*stab + 32.7_WP*(1.0_WP-stab)) ! L & Y eq. (6c-6d)
         !!
         !!
         !! Shifting the neutral 10m transfer coefficients to (zu,zeta_u) :
!        xct = 1. + sqrt_Cd_n10/kappa*(log(zu/10.) - psi_m(zeta_u))
         xct = 1._WP + sqrt_Cd_n10/kappa*(log(zu/10._WP) - zpsi_m)
         Cd = Cd_n10/(xct*xct) ; sqrt_Cd = sqrt(Cd)
         !!
!        xlogt = log(zu/10.) - psi_h(zeta_u)
         xlogt = log(zu/10._WP) - zpsi_hu
         !!
         xct = 1._WP + Ch_n10*xlogt/kappa/sqrt_Cd_n10
         Ch  = Ch_n10*sqrt_Cd/sqrt_Cd_n10/xct
         !!
         xct = 1._WP + Ce_n10*xlogt/kappa/sqrt_Cd_n10
         Ce  = Ce_n10*sqrt_Cd/sqrt_Cd_n10/xct
         !!
         !!
      end do
      !!

      !
   END SUBROUTINE core_coeff_2z

   FUNCTION psi_h( zta )
      !! Psis, L & Y eq. (8c), (8d), (8e)
      !-------------------------------------------------------------------------------
      real(wp)             :: X2
      real(wp)             :: X
      real(wp)             :: stabit
      !
      real(wp), intent(in) ::   zta
      real(wp)             ::   psi_h
      !-------------------------------------------------------------------------------

      X2 = sqrt(abs(1._WP - 16._WP*zta))  ;  X2 = max(X2 , 1._WP) ;  X  = sqrt(X2)
      stabit    = 0.5_WP + sign(0.5_wp,zta)
      psi_h = -5._WP*zta*stabit  &                                       ! Stable
         &    + (1._WP - stabit)*(2._WP*log( (1._WP + X2)/2._WP ))                 ! Unstable
   END FUNCTION psi_h

   FUNCTION psi_m( zta )
   !! Psis, L & Y eq. (8c), (8d), (8e)
   !-------------------------------------------------------------------------------
      real(wp)             :: X2
      real(wp)             :: X
      real(wp)             :: stabit
      !!
      real(wp), intent(in) ::   zta

      real(wp), parameter  :: pi = 3.141592653589793_wp
      real(wp)             :: psi_m
      !-------------------------------------------------------------------------------

      X2 = sqrt(abs(1._WP - 16._WP*zta))  ;  X2 = max(X2 , 1.0_WP) ;  X  = sqrt(X2)
      stabit    = 0.5_WP + sign(0.5_wp,zta)
      psi_m = -5._WP*zta*stabit  &                                                          ! Stable
         &    + (1._WP - stabit)*(2.0_WP*log((1._WP + X)/2_WP) + log((1._WP + X2)/2.0_WP) - 2.0_WP*atan(X) + pi/2.0_WP)  ! Unstable

      !
   END FUNCTION psi_m

   SUBROUTINE getcoeffld(rdate)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE getcoeffld ***
      !!
      !! ** Purpose : read fields from files, interpolate on model mesh and prepare interpolation coefficients
      !! ** Method  :
      !! ** Action  :
      !!----------------------------------------------------------------------
      IMPLICIT NONE
!      integer,intent(in)   :: idate ! initialization date
!      integer,intent(in)   :: isec ! initialization seconds
      real(wp),intent(in)   :: rdate ! initialization date

      integer              :: iost !I/O status
      integer              :: ncid      ! netcdf file id
      ! ID dimensions and variables:
      integer              :: id_data
      integer              :: nf_start(4)
      integer              :: nf_edges(4)
!      integer              :: zero_year,yyyy,mm,dd
!      character(len = 256) :: att_string ! attribute
      integer              :: fld_idx, i,j,ii, ip1, jp1, extrp
      integer              :: sbc_alloc, itot

      real(wp)             :: denom, x1, x2, y1, y2, x, y
      real(wp)             :: now_date

      real(wp), allocatable, dimension(:,:)  :: sbcdata1,sbcdata2
      real(wp)             :: data1,data2
      real(wp)             :: delta_t   ! time(t_indx) - time(t_indx+1)

      integer              :: elnodes(4) !4 nodes from one element
      integer              :: numnodes   ! nu,ber of nodes in elem (3 for triangle, 4 for ... )

      ALLOCATE( sbcdata1(nc_Nlon,nc_Nlat), sbcdata2(nc_Nlon,nc_Nlat),&
                &      STAT=sbc_alloc )
!                data1(elem2D),data2(elem2D), &
      if( sbc_alloc /= 0 )   STOP 'getcoeffld: failed to allocate arrays'

      ! find time index in files
      now_date = rdate
      call binarysearch(nc_Ntime,nc_time,now_date,t_indx)
      if ( now_date < nc_time(1) ) then  ! NO extrapolation back in time
         t_indx = 1
         t_indx_p1 = t_indx
         delta_t = 1.0_wp
      else
      if ( t_indx < nc_Ntime ) then
         t_indx_p1 = t_indx + 1
         delta_t = nc_time(t_indx_p1) - nc_time(t_indx)
      else ! NO extrapolation to future
         t_indx_p1 = t_indx
         delta_t = 1.0_wp
         WRITE(*,*) '     WARNING:  time lager than last time step in forcing file,'
         WRITE(*,*) '        last time step will be used as a constant field'
         one_field = .true.
      end if
      end if
      itot = i_totfl
!            if (nm_calc_flux==0) then
      itot = 8
!            end if

      do fld_idx = 1, itot
         !open file sbc_flfi
         iost = nf_open(sbc_flfi(fld_idx)%file_name,NF_NOWRITE,ncid)
         call check_nferr(iost,sbc_flfi(fld_idx)%file_name)
         ! get variable id
         iost = nf_inq_varid(ncid, sbc_flfi(fld_idx)%var_name, id_data)
         call check_nferr(iost,sbc_flfi(fld_idx)%file_name)
         !read data from file
         nf_start(1)=1
         nf_edges(1)=nc_Nlon
         nf_start(2)=1
         nf_edges(2)=nc_Nlat
         nf_start(3)=t_indx
         nf_edges(3)=1
         iost = nf_get_vara_double(ncid, id_data, nf_start, nf_edges, sbcdata1)
         call check_nferr(iost,sbc_flfi(fld_idx)%file_name)
         ! read next time step in file (check for +1 done before)
         nf_start(3)=t_indx_p1
         nf_edges(3)=1
         iost = nf_get_vara_double(ncid, id_data, nf_start, nf_edges, sbcdata2)
         call check_nferr(iost,sbc_flfi(fld_idx)%file_name)
         !flip data in case of lat from -90 to 90
!!!! WARNING
         if ( flip_lat == 1 ) then
             sbcdata1=sbcdata1(:,nc_Nlat:1:-1)
             sbcdata2=sbcdata2(:,nc_Nlat:1:-1)
         end if


         ! bilinear space interpolation, and time interpolation ,
         ! data is assumed to be sampled on a regular grid
!!$OMP PARALLEL
!!$OMP DO
         do ii = 1, nod2D
            i = bilin_indx_i(ii)
            j = bilin_indx_j(ii)
            ip1 = i + 1
            jp1 = j + 1

!!WARNING !! get coordinates of element (mean of nodes in elem)
!            elnodes = elem2D_nodes(:,ii) !! 4 nodes in element
!            numnodes = 4
!            if(elnodes(1)==elnodes(4)) numnodes = 3  !! set to 3 if we have triangle
!            y = sum(coord_nod2d(2,elnodes(1:numnodes)))/dble(numnodes) /rad
!            x = sum(coord_nod2d(1,elnodes(1:numnodes)))/dble(numnodes) /rad
!  use these lines in case if we use sbc on nodes
            x  = coord_nod2d(1,ii)/rad
            y  = coord_nod2d(2,ii)/rad
            extrp = 0
            if ( i == 0 ) then
               i   = nc_Nlon
               ip1 = i
               extrp = extrp + 1
            end if
            if ( i == -1 ) then
               i   = 1
               ip1 = i
               extrp = extrp + 1
            end if
            if ( j == 0 ) then
               j   = nc_Nlat
               jp1 = j
               extrp = extrp + 2
            end if
            if ( j == -1 ) then
               j   = 1
               jp1 = j
               extrp = extrp + 2
            end if

            x1 = nc_lon(i)
            x2 = nc_lon(ip1)
            y1 = nc_lat(j)
            y2 = nc_lat(jp1)

            if ( extrp == 0 ) then
            ! if point inside forcing domain
               denom = (x2 - x1)*(y2 - y1)
               data1 = ( sbcdata1(i,j)   * (x2-x)*(y2-y)   + sbcdata1(ip1,j)    * (x-x1)*(y2-y) + &
                     sbcdata1(i,jp1) * (x2-x)*(y-y1)   + sbcdata1(ip1, jp1) * (x-x1)*(y-y1)     ) / denom
               data2 = ( sbcdata2(i,j)   * (x2-x)*(y2-y)   + sbcdata2(ip1,j)    * (x-x1)*(y2-y) + &
                     sbcdata2(i,jp1) * (x2-x)*(y-y1)   + sbcdata2(ip1, jp1) * (x-x1)*(y-y1)     ) / denom
            else if ( extrp == 1 ) then !  "extrapolation" in x direction
               denom = (y2 - y1)
               data1 = ( sbcdata1(i,j)   * (y2-y)   + sbcdata1(ip1, jp1) * (y-y1) ) / denom
               data2 = ( sbcdata2(i,j)   * (y2-y)   + sbcdata2(ip1, jp1) * (y-y1) ) / denom
            else if ( extrp == 2 ) then !  "extrapolation" in y direction
               denom = (x2 - x1)
               data1 = ( sbcdata1(i,j)   * (x2-x)   + sbcdata1(ip1, jp1) * (x-x1) ) / denom
               data2 = ( sbcdata2(i,j)   * (x2-x)   + sbcdata2(ip1, jp1) * (x-x1) ) / denom
            else if ( extrp == 3 ) then !  "extrapolation" in x and y direction
               data1 = sbcdata1(i,j)
               data2 = sbcdata2(i,j)
            end if
            ! calculate new coefficients for interpolations
            coef_a(fld_idx, ii) = ( data2 - data1 ) / delta_t !( nc_time(t_indx+1) - nc_time(t_indx) )
            coef_b(fld_idx, ii) = data1 - coef_a(fld_idx, ii) * nc_time(t_indx)

         end do !ii
!!$OMP END DO
!!$OMP END PARALLEL
         iost = nf_close(ncid)
         call check_nferr(iost,sbc_flfi(fld_idx)%file_name)

      end do !fld_idx

      DEALLOCATE( sbcdata1, sbcdata2 )


   END SUBROUTINE getcoeffld

   SUBROUTINE data_timeinterp(rdate)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE data_timeinterp ***
      !!
      !! ** Purpose : interpolation of fields(interpolated on model grid) from IN files in time
      !! ** Method  :
      !! ** Action  :
      !!----------------------------------------------------------------------
      IMPLICIT NONE
!      integer,intent(in) :: ndate ! date
!      integer,intent(in) :: nsec  ! seconds
      real(wp),intent(in)    :: rdate  ! seconds

     ! assign data from interpolation to taux and tauy
      integer            :: fld_idx, i,j,ii
      real(wp)           :: now_date

!      now_date = ndate+nsec/86400.0_wp
      now_date = rdate
!!$OMP PARALLEL
!!$OMP DO
      do i = 1, nod2D
         do fld_idx = 1, i_totfl

            atmdata(fld_idx,i) = now_date * coef_a(fld_idx,i) + coef_b(fld_idx,i)

         end do !fld_idx
      end do !elem2D
!!$OMP END DO
!!$OMP END PARALLEL
   END SUBROUTINE data_timeinterp

   SUBROUTINE sbc_ini
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_ini ***
      !!
      !! ** Purpose : inizialization of ocean forcing
      !! ** Method  :
      !! ** Action  :
      !!----------------------------------------------------------------------
      IMPLICIT NONE

      integer            :: idate ! initialization date
      integer            :: isec  ! initialization seconds
      integer            :: iost  ! I/O status
      integer            :: sbc_alloc                   !: allocation status

      real(wp)           :: tx, ty



      namelist/nam_sbc/ nm_nc_secstep, nm_sbc, nm_tauwind, nm_xwind0, nm_ywind0,  &
                        nm_xwind_file, nm_ywind_file, nm_humi_file, nm_qsr_file, &
                        nm_qlw_file, nm_tair_file, nm_prec_file, nm_humi0, nm_qsr0, &
                        nm_qlw0, nm_tair, nm_prec,nm_xwind_var,nm_ywind_var,nm_humi_var, &
                        nm_qsr_var, nm_qlw_var, nm_tair_var, nm_prec_var, &
                        nm_mslp_var, nm_mslp_file, &
                        nm_cloud_var, nm_cloud_file, &
                        depth_swr, nm_prec_coef, nm_net_flux,  &
                        nm_nc_iyear, nm_nc_imm, nm_nc_idd !nm_calc_flux,nm_sbc_ftype,

      write(*,*) "Start: Ocean forcing inizialization."

      idate = time_jd0
      isec  = 86400.0_WP *(time_jd0 - idate)
! file for debuging
!      open(unit=unit_deb,file='obc_deb.dat')

      ! OPEN and read namelist for SBC
      open( unit=nm_sbc_unit, file='namelist_bc.nml', form='formatted', access='sequential', status='old', iostat=iost )
      if( iost == 0 ) then
         WRITE(*,*) '     file   : ', 'namelist_bc.nml',' open ok'
      else
         WRITE(*,*) 'ERROR: --> bad opening file   : ', 'namelist_bc.nml',' ; iostat=',iost
         STOP 'ERROR: --> sbc_ini'
      endif
      READ( nm_sbc_unit, nml=nam_sbc, iostat=iost )
      close( nm_sbc_unit )
      if( nm_sbc == -1 ) then
         ! IF module not in use
         return
      endif
      write(*,*) "Start: Ocean forcing inizialization."
      write(*,*) "   Surface boundary conditions parameters:"
      write(*,*) "      nm_sbc        = ", nm_sbc,"   ! 1= constant, 2=from file, -1=module not in use"
      write(*,*) "      nm_tauwind    = ", nm_tauwind ," ! if 1 = wind stress, 2 = wind 10 m  .    "
    !      write(*,*) "      nm_sbc_ftype  = ", nm_sbc_ftype ," ! input file type: 1 = ASCII, 2 = netcdf   "
      if (nm_sbc==1) then
         write(*,*) "      nm_xwind0     = ", nm_xwind0 ," ! constant 10m. wind value in i-direction/stress, if nm_sbc=1   "
         write(*,*) "      nm_ywind0     = ", nm_ywind0 ," ! constant 10m. wind value in j-direction, if nm_sbc=1  "
         write(*,*) "      nm_humi0      = ", nm_humi0 ," ! constant humidity "
         write(*,*) "      nm_qsr0       = ", nm_qsr0 ," ! constant solar heat  "
         write(*,*) "      nm_qlw0       = ", nm_qlw0 ," ! constant Long wave "
         write(*,*) "      nm_tair       = ", nm_tair," ! constant 2m air temperature "
         write(*,*) "      nm_prec       = ", nm_prec ," ! constant total precipitation "
         write(*,*) " WARNING:: only wind constants in use, rest set to 0 !!!    "
      else
         write(*,*) "      nm_xwind_file = ", trim(nm_xwind_file) ," ! name of file with winds, if nm_sbc=2 "
         write(*,*) "      nm_ywind_file = ", trim(nm_ywind_file) ," ! name of file with winds, if nm_sbc=2 "
         write(*,*) "      nm_humi_file  = ", trim(nm_humi_file) ," ! name of file with humidity "
         write(*,*) "      nm_qsr_file   = ", trim(nm_qsr_file) ," ! name of file with solar heat "
         write(*,*) "      nm_qlw_file   = ", trim(nm_qlw_file) ," ! name of file with Long wave "
         write(*,*) "      nm_tair_file  = ", trim(nm_tair_file) ," ! name of file with 2m air temperature "
         write(*,*) "      nm_prec_file  = ", trim(nm_prec_file) ," ! name of file with total precipitation "
         write(*,*) "      nm_mslp_file  = ", trim(nm_mslp_file)," !air_pressure_at_sea_level "
         write(*,*) "      nm_cloud_file  = ", trim(nm_cloud_file)," !clouds "
         write(*,*) "      nm_xwind_var  = ", trim(nm_xwind_var) ," ! name of variable in file with wind "
         write(*,*) "      nm_ywind_var  = ", trim(nm_ywind_var) ," ! name of variable in file with wind "
         write(*,*) "      nm_humi_var   = ", trim(nm_humi_var) ," ! name of variable in file with humidity  "
         write(*,*) "      nm_qsr_var    = ", trim(nm_qsr_var) ," ! name of variable in file with solar heat "
         write(*,*) "      nm_qlw_var    = ", trim(nm_qlw_var) ," ! name of variable in file with Long wave "
         write(*,*) "      nm_tair_var   = ", trim(nm_tair_var) ," ! name of variable in file with 2m air temperature "
         write(*,*) "      nm_prec_var   = ", trim(nm_prec_var) ," ! name of variable in file with total precipitation  "
         write(*,*) "      nm_mslp_var   = ", trim(nm_mslp_var) ," ! name of variable in file with air_pressure_at_sea_level "
         write(*,*) "      nm_cloud_var   = ", trim(nm_cloud_var) ," ! name of variable in file with clouds "
      endif
      write(*,*) "      depth_swr     = ", depth_swr ," !  depth of swr penetration"
      write(*,*) "      nm_net_flux   = ", nm_net_flux," !  key for downward longwave heat over the ocean: 0 - , 1 - Net"
          write(*,*) "      nm_prec_coef  = ", nm_prec_coef,&
                                                " !  precipitation will be devide by this constant(3600-CoastDat,1-NCEP)"
      write(*,*) "      nm_nc_secstep = ", nm_nc_secstep ,&
              " !  time units coef (86400 CoastDat, 24 NCEP),netcdf time to [day] by Tnetcdf/nm_nc_secstep"
      write(*,*) "      nm_nc_iyear   = ", nm_nc_iyear ," ! initial year of time axis in netCDF (1948 like CoastDat,1800 NCEP)"
      write(*,*) "      nm_nc_imm     = ", nm_nc_imm ," ! initial month of time axis in netCDF "
      write(*,*) "      nm_nc_idd     = ", nm_nc_idd ," ! initial day of time axis in netCDF "
      ALLOCATE( coef_a(i_totfl,nod2D), coef_b(i_totfl,nod2D), &
              & atmdata(i_totfl,nod2D), &
                   &      STAT=sbc_alloc )
      if( sbc_alloc /= 0 )   STOP 'sbc_ini: failed to allocate arrays'

      ALLOCATE( bilin_indx_i(nod2D),bilin_indx_j(nod2D), &
                   &      STAT=sbc_alloc )
!            & qns(nod2D), emp(nod2D), qsr(nod2D),&!qsr3d(nsigma-1,nod2D),  &
!              & qns_2(nod2D), emp_2(nod2D), qsr_2(nod2D),  &
!              & taux_node_2(nod2D), tauy_node_2(nod2D),  &

      if( sbc_alloc /= 0 )   STOP 'sbc_ini: failed to allocate arrays'

! IF constant wind/stress
      if( nm_sbc == 1 ) then
         STOP "ERROR: sbc_ini: nm_sbc == 1, constant wind/stress is under developing "

         if( nm_tauwind == 1 ) then
            write(*,*) "WARNING:: nm_sbc == 1 and  nm_tauwind == 1, "
            write(*,*) "           to calculate wind at 10m will be used simple Large and Pond (1981) _revers_ parametrization."
            call stress2wind_scalar(nm_xwind0, nm_ywind0, tx, ty)
            taux = nm_xwind0
            tauy = nm_ywind0
            taux_node = nm_xwind0
            tauy_node = nm_ywind0

         else
            write(*,*) "WARNING:: nm_sbc == 1 and  nm_tauwind == 2, "
            write(*,*) "           to calculate wind stress will be used simple Large and Pond (1981) parametrization."
!            STOP "ERROR: sbc_ini"
            call wind2stress_scalar(nm_xwind0, nm_ywind0, tx, ty)
            taux = tx
            tauy = ty
            taux_node = tx
            tauy_node = ty
            windx = nm_xwind0
            windy = nm_ywind0
         endif
      endif

! ========================== NC PART start===============================================================
      if( nm_sbc == 2 ) then
         call nc_sbc_ini(idate,isec)
      endif


      write(*,*) "DONE:  Ocean forcing inizialization."
   END SUBROUTINE sbc_ini

   SUBROUTINE sbc_do
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_do ***
      !!
      !! ** Purpose : provide at each time-step: wind stress, ...
      !! ** Method  :
      !! ** Action  :
      !!----------------------------------------------------------------------
      IMPLICIT NONE

      integer            :: rdate ! date
      integer            :: rsec  ! seconds

      if( nm_sbc == -1 ) then
         ! IF module not in use
         return
      endif

      rdate = time_jd
      rsec  = 86400.0_WP *(time_jd - rdate)


         if( .not. one_field ) then
            ! IF more field available
            if( time_jd > nc_time(t_indx_p1) ) then
               ! get new coefficients for time interpolation on model grid for all data
               call getcoeffld(time_jd)
            endif
         endif
         ! interpolate in time
         call data_timeinterp(time_jd)
         ! change model ocean parameters
         call apply_atm_fluxes

   END SUBROUTINE sbc_do

   SUBROUTINE wind2stress_scalar(U10x,U10y,tx,ty)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE  wind2stress ***
      !!
      !! ** Purpose : convert wind 10m to wind stress (return scalar)
      !! ** Method  :
      !! ** Action  :
      !!----------------------------------------------------------------------
   IMPLICIT NONE
      real(wp), intent(in)  :: U10x ! wind X direction
      real(wp), intent(in)  :: U10y ! wind Y direction
      real(wp), intent(out) :: tx   ! tau X direction
      real(wp), intent(out) :: ty   ! tau X direction

      real(wp)             :: U10

!     Large and Pond (1981), J. Phys. Oceanog., 11, 324-336.
!     Tau = Cd * ro_air * U10^2
!     U10    = wind speed at 10 m above the sea surface
!     ro_air =  1.22 kg m-3
!     Cd     = dimensionless drag coefficient :
!         1.2 x 10^-3 for 4 < U10 < 11 m s-1
!         10^-3 (0.49 + 0.065 U10)     for     11 < U10 < 25 m s-1
      U10 = sqrt( U10x*U10x + U10y*U10y )
      if ( U10 <= 11 ) then
         tx = 1.2d-3 * 1.22 * U10 * U10x
         ty = 1.2d-3 * 1.22 * U10 * U10y
      else
         tx = 1.d-3 * (0.49 + 0.065 * U10) * 1.22 * U10 * U10x
         ty = 1.d-3 * (0.49 + 0.065 * U10) * 1.22 * U10 * U10y
      endif

   END SUBROUTINE wind2stress_scalar

   SUBROUTINE stress2wind_scalar(tx, ty, U10x, U10y)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE  wind2stress ***
      !!
      !! ** Purpose : convert wind stress to wind 10m (return scalar)
      !! ** Method  :  it is very bad way to do it, sorry
      !! ** Action  :
      !!----------------------------------------------------------------------
   IMPLICIT NONE
      real(wp), intent(out)  :: U10x ! wind X direction
      real(wp), intent(out)  :: U10y ! wind Y direction
      real(wp), intent(in) :: tx   ! tau X direction
      real(wp), intent(in) :: ty   ! tau X direction

!     Tau = Cd * ro_air * U10^2
!     U10    = wind speed at 10 m above the sea surface
!     ro_air =  1.22 kg m-3
!     Cd     = dimensionless drag coefficient : 1.2 x 10^-3
!

      U10y = sqrt(ty / (1.22 * 1.2d-3 * sqrt(tx*tx/(ty*ty)+1) ))
      U10x = sqrt(tx / (1.22 * 1.2d-3 * sqrt(ty*ty/(tx*tx)+1) ))


   END SUBROUTINE stress2wind_scalar

   FUNCTION julday(yyyy,mm,dd)

   IMPLICIT NONE
      integer, INTENT(IN) :: mm, dd, yyyy
      integer             :: julday
! In this routine julday returns the Julian Day Number that begins at noon of the calendar
!    date specified by month mm , day dd , and year yyyy , all integer variables. Positive year
!    signifies A.D.; negative, B.C. Remember that the year after 1 B.C. was 1 A.D. (from Num. Rec.)
      integer, PARAMETER  :: IGREG=15+31*(10+12*1582)
! Gregorian Calendar adopted Oct. 15, 1582.
      integer             :: ja,jm,jy

      jy = yyyy
      if (jy == 0) STOP 'julday: there is no year zero'
      if (jy < 0) jy=jy+1
      if (mm > 2) then
         jm=mm+1
      else
         jy=jy-1
         jm=mm+13
      endif
      julday=int(365.25_wp*jy)+int(30.6001_wp*jm)+dd+1720995
!Test whether to change to Gregorian Calendar.
      if (dd+31.0_WP*(mm+12.0_WP*yyyy) >= IGREG) then
         ja=int(0.01_WP*jy)
         julday=julday+2.0_WP-ja+int(0.25_wp*ja)
      end if
   END FUNCTION julday


   SUBROUTINE calendar_date(julian,yyyy,mm,dd)

!  Converts a Julian day to a calendar date (year, month and day). Numerical Recipes
   IMPLICIT NONE
!
      integer,intent(in) :: julian
      integer            :: yyyy,mm,dd

      integer, parameter :: IGREG=2299161
      integer            :: ja,jb,jc,jd,je
      real(wp)           :: x
!
!-----------------------------------------------------------------------
      if (julian >= IGREG ) then
         x = ((julian-1867216)-0.25)/36524.25
         ja = julian+1+int(x)-int(0.25*x)
      else
         ja = julian
      end if

      jb = ja+1524
      jc = int(6680 + ((jb-2439870)-122.1)/365.25)
      jd = int(365*jc+(0.25*jc))
      je = int((jb-jd)/30.6001)

      dd = jb-jd-int(30.6001*je)
      mm = je-1
      if (mm > 12) mm = mm-12
      yyyy = jc - 4715
      if (mm > 2) yyyy = yyyy-1
      if (yyyy <= 0) yyyy = yyyy-1

      return
   END SUBROUTINE calendar_date

   SUBROUTINE sbc_end

      IMPLICIT NONE

      if( nm_sbc == 2 ) then
         DEALLOCATE( nc_lon, nc_lat, nc_time)
      endif

      DEALLOCATE( coef_a, coef_b, atmdata, &
                  &  bilin_indx_i, bilin_indx_j)!,  &
!                  &  qns, emp, qsr)




   END SUBROUTINE sbc_end

   SUBROUTINE check_nferr(iost,fname)
   IMPLICIT NONE
      character(len=256), intent(in) :: fname
      integer, intent(in) :: iost

      if (iost .ne. NF_NOERR) then
         write(*,*) 'ERROR: I/O status= "',trim(nf_strerror(iost)),'";',iost,' file= ',fname
         STOP 'ERROR: stop'
      endif
   END SUBROUTINE

   SUBROUTINE binarysearch(length, array, value, ind)!, delta)
      ! Given an array and a value, returns the index of the element that
      ! is closest to, but less than, the given value.
      ! Uses a binary search algorithm.
      ! "delta" is the tolerance used to determine if two values are equal
      ! if ( abs(x1 - x2) <= delta) then
      !    assume x1 = x2
      ! endif
      !org. source from: https://github.com/cfinch/Shocksolution_Examples/blob/master/FORTRAN/BilinearInterpolation/interpolation.f90

      IMPLICIT NONE
      integer,  intent(in) :: length
      real(wp), dimension(length), intent(in) :: array
      real(wp), intent(in) :: value
   !   real, intent(in), optional :: delta

   !   integer :: binarysearch
      integer, intent(out) :: ind

      integer :: left, middle, right
      real(wp):: d

   !   if (present(delta) .eqv. .true.) then
   !      d = delta
   !   else
      d = 1e-9
   !   endif
      left = 1
      right = length
      do
         if (left > right) then
            exit
         endif
         middle = nint((left+right) / 2.0_WP)
         if ( abs(array(middle) - value) <= d) then
            ind = middle
            return
         else if (array(middle) > value) then
            right = middle - 1
         else
            left = middle + 1
         end if
      end do
      ind = right

   END SUBROUTINE binarysearch

END MODULE fv_sbc
