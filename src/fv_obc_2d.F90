MODULE fv_obc_2d
   !!===========================================================================
   !! Open boundary conditions for 2D variables (ssh):
   !!===========================================================================
   !! History: 0.1 ! 03/2016 I. Kuznetsov
   !!
   !!
   !! public:
   !!   obc_2d_ini  -- inizialization of open boundary
   !!   obc_2d_do   -- provide a obc (open boundary conditions) each time step
   !!
   USE o_ARRAYS
   USE o_MESH
   USE o_PARAM


   IMPLICIT NONE

   include 'netcdf.inc'

   public  obc_2d_ini  ! routine called before 1st time step (open files, read namelist,...)
   public  obc_2d_do   ! routine called each time step to provide a sbc fileds (wind,...)
   public  obc_2d_end  ! routine called after last time step

   private


   ! namelists
   integer, save  :: nm_obc_unit     = 201       ! unit to open namelist file
  !============== namelistatmdata variables ================
   character(len=256), save   :: nm_sshobc_file = 'sshobc' ! name of file with ssh at OB

   character(len=34), save   :: nm_sshobc_var = 'ssh'  ! name of variable in file with ssh

   real(wp), save :: nm_obc_tau_relax = 86400.0 ! relaxation coeficient for OB [sec]
   integer, save  :: nm_obc_2d   ! switch of module if -1

   ! ========== netCDF time param
   integer, save :: nm_nc_iyear = 1948    ! initial year of time axis in netCDF (1948 like CoastDat,1800 NCEP)
   integer, save :: nm_nc_imm = 1         ! initial month of time axis in netCDF
   integer, save :: nm_nc_idd = 1         ! initial day of time axis in netCDF
   real(wp), save :: nm_nc_secstep = 86400.0 ! time units coef (86400 CoastDat, 24 NCEP)

   ! ========== interpolation coeficients

   real(wp), allocatable, save, dimension(:,:)   :: coef_b ! time inerp coef. b (x=a*t+b) (fld_indx,node)
   real(wp), allocatable, save, dimension(:,:)   :: coef_a ! time inerp coef. a (x=a*t+b) (fld_indx,node)

   logical, save :: one_field = .false.! only one time step used for openboundary

   real(wp), allocatable, save, dimension(:,:)   :: obcdata ! OB data for current time step
   integer,  allocatable, save, dimension(:)     :: index_nod2D_ob ! indexes of nods at OB

   real(wp), save        :: intdata ! interpolated OB data for current time step and node (call apply_atm_fluxes)


!============== NETCDF ==========================================
!structure similar to sbc (so fluxes here is concentrations)

   integer, parameter :: i_totc = 1 ! total number of 2D (temperature,...) values (2=T,S)
   integer, parameter :: i_ssh  = 1 ! index of sea level h. [m]


   type, public ::   flfi_type    !flux file informations
      character(len = 256) :: file_name ! file name
      character(len = 34)  :: var_name  ! variable name in the NetCDF file
   end type flfi_type

  type(flfi_type),save, dimension(i_totc) :: obc_flfi  !array for information about flux files

  ! arrays of time, lon and lat in INfiles
   real(wp), allocatable, save, dimension(:)  :: nc_lon  ! coordinates of nodes in infile
   real(wp), allocatable, save, dimension(:)  :: nc_lat  ! coordinates of nodes in infile
   real(wp), allocatable, save, dimension(:)  :: nc_node ! index of nodes in infile corespond to node in mesh
   real(wp), allocatable, save, dimension(:)  :: nc_time ! time array
!!   real(wp), allocatable, save, dimension(:)  :: nc_depth! depth (z) in netcdf file [m] , assume we have z-coordinate system
  ! lenght of arrays in INfiles
   integer,save              :: nc_Nnode
!!   integer,save              :: nc_Ndepth
   integer,save              :: nc_Ntime
   ! time index for NC time array
   integer,save              :: t_indx    ! now time index in nc_time array
   integer,save              :: t_indx_p1 ! now time index +1 in nc_time array



!============== NETCDF ==========================================
CONTAINS
   SUBROUTINE nc_readTimeGrid(flf)
   ! Read time array and grid from nc file
      IMPLICIT NONE

      type(flfi_type),intent(in) :: flf

      integer              :: iost !I/O status
      integer              :: ncid      ! netcdf file id
      integer              :: i, warn
      ! ID dimensions and variables:
!      integer              :: id_depth
      integer              :: id_node
      integer              :: id_time
      integer              :: id_lon
      integer              :: id_lat
      integer              :: id_noded
      integer              :: id_timed
!      integer              :: id_depthd
!      integer              :: nf_dims(4) ! dimensions (temporal)
      integer              :: nf_start(4)
      integer              :: nf_edges(4)
      integer              :: zero_year,yyyy,mm,dd
      character(len = 256) :: att_string ! attribute
      integer              :: obc_alloc                   !: allocation status



      !open file
      iost = nf_open(flf%file_name,NF_NOWRITE,ncid)
      call check_nferr(iost,flf%file_name)

      ! get dimensions
      iost = nf_inq_dimid(ncid, "node", id_noded)
      call check_nferr(iost,flf%file_name)
!      iost = nf_inq_dimid(ncid, "depth", id_depthd)
!      call check_nferr(iost,flf%file_name)
      iost = nf_inq_dimid(ncid, "time", id_timed)
      call check_nferr(iost,flf%file_name)

      ! get variable id
      iost = nf_inq_varid(ncid, "node", id_node)
      call check_nferr(iost,flf%file_name)
!      iost = nf_inq_varid(ncid, "depth", id_depth)
!      call check_nferr(iost,flf%file_name)
      iost = nf_inq_varid(ncid, "time", id_time)
      call check_nferr(iost,flf%file_name)
      iost = nf_inq_varid(ncid, "lon", id_lon)
      call check_nferr(iost,flf%file_name)
      iost = nf_inq_varid(ncid, "lat", id_lat)
      call check_nferr(iost,flf%file_name)
      !  get dimensions size
!       iost = nf_inq_dimlen(ncid, id_depthd, nc_Ndepth)
!       call check_nferr(iost,flf%file_name)
      iost = nf_inq_dimlen(ncid, id_noded, i)
      call check_nferr(iost,flf%file_name)
      if( i /= nc_Nnode ) STOP 'nc_readTimeGrid: size of OB nodes in file .neq. to mesh'
      iost = nf_inq_dimlen(ncid, id_timed, nc_Ntime)
      call check_nferr(iost,flf%file_name)

      ALLOCATE( nc_lon(nc_Nnode), nc_lat(nc_Nnode), &
                nc_time(nc_Ntime), nc_node(nc_Nnode), &
                &      STAT=obc_alloc )
      if( obc_alloc /= 0 )   STOP 'nc_readTimeGrid: failed to allocate arrays'
 !(fld_indx,node,depth)
      ALLOCATE( coef_a(i_totc,nc_Nnode), coef_b(i_totc,nc_Nnode), &
              & obcdata(i_totc,nc_Nnode), &
                   &      STAT=obc_alloc )
      if( obc_alloc /= 0 )   STOP 'nc_readTimeGrid: failed to allocate arrays'
      ALLOCATE( eta_n_obc(nc_Nnode),STAT=obc_alloc )
      if( obc_alloc /= 0 )   STOP 'nc_readTimeGrid: failed to allocate arrays'

   !read variables from file
   !time
      nf_start(1)=1
      nf_edges(1)=nc_Ntime
      iost = nf_get_vara_double(ncid, id_time, nf_start, nf_edges, nc_time)
      call check_nferr(iost,flf%file_name)
      ! convert time to days
      !!NCEP     nc_time = nc_time / 24.0 + julday(1800,1,1)
      nc_time = nc_time / nm_nc_secstep + julday(nm_nc_iyear,nm_nc_imm,nm_nc_idd)
!!! WARNING : move time forward for dt/2 , last is a last+dt from last -1
!      if (nc_Ntime > 1) then
!         do i = 1, nc_Ntime-1
!            nc_time(i) = (nc_time(i+1) + nc_time(i))/2.0
!         end do

!         nc_time(nc_Ntime) = nc_time(nc_Ntime) + (nc_time(nc_Ntime) - nc_time(nc_Ntime-1))/2.0
!      end if
   ! coordinates
      nf_start(1)=1
      nf_edges(1)=nc_Nnode
      iost = nf_get_vara_double(ncid, id_lat, nf_start, nf_edges, nc_lat)
      call check_nferr(iost,flf%file_name)
      iost = nf_get_vara_double(ncid, id_lon, nf_start, nf_edges, nc_lon)
      call check_nferr(iost,flf%file_name)
      iost = nf_get_vara_double(ncid, id_node, nf_start, nf_edges,nc_node)
      call check_nferr(iost,flf%file_name)

      warn = 0
      do i = 1, nc_Nnode
         if ( nc_node(i) /= index_nod2D_ob(i) ) then
            if ( warn == 0 ) then
               write(*,*) 'WARNING:: OB nodes in Netcdf file are not the same as in mesh. node:',index_nod2D_ob(i)
               warn = 1
            endif
         endif
      enddo
      do i = 1, nc_Nnode
         if ( in_obn(i) /= index_nod2D_ob(i) ) then
               write(*,*) 'ERROR:: OB nodes in fv_obc_2d are not the same as in fv_main; node:',index_nod2D_ob(i)
               STOP
         endif
      enddo
!         x  = coord_nod2d(1,i)/rad
!         y  = coord_nod2d(2,i)/rad

      iost = nf_close(ncid)
      call check_nferr(iost,flf%file_name)

   END SUBROUTINE nc_readTimeGrid

   SUBROUTINE nc_obc_ini_fillnames(yyear)
      character(len=4),intent(in)   :: yyear

      !! ** Purpose : Fill names of obc_flfi array (file names and variable names)

      !prepare proper nc file (add year and .nc to the end of the file name from namelist
      write(obc_flfi(i_ssh)%file_name, *) trim(nm_sshobc_file),yyear,'.nc'

      obc_flfi(i_ssh)%file_name=ADJUSTL(trim(obc_flfi(i_ssh)%file_name))

      obc_flfi(i_ssh)%var_name=ADJUSTL(trim(nm_sshobc_var))
   END SUBROUTINE nc_obc_ini_fillnames

   SUBROUTINE nc_obc_ini(rdate)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE nc_obc_ini ***
      !!
      !! ** Purpose : inizialization of OB from NETCDF file
      !!----------------------------------------------------------------------

      IMPLICIT NONE
      integer  :: idate ! initialization date
!      integer,intent(in) :: isec ! initialization seconds
      real(wp),intent(in)   :: rdate ! initialization date

      character(len=4)   :: yyear
      integer            :: yyyy,mm,dd

      integer            :: i
      integer            :: obc_alloc

      integer            :: elnodes(4) !4 nodes from one element





      idate = rdate
      ! get ini year; Fill names of obc_flfi
      call calendar_date(idate,yyyy,dd,mm)
      write(yyear,"(I4)") yyyy
      call nc_obc_ini_fillnames(yyear)

      ! we assume that all NetCDF files have identical grid and time variable
      ! read axis, sizes, allocate memory
      call nc_readTimeGrid(obc_flfi(i_ssh))



      ! get first coeficients for time inerpolation on model grid for all datas
      call getcoeffld(rdate)

      ! interpolate in time

      call data_timeinterp(rdate)

      call apply_atm_fluxes


   END SUBROUTINE nc_obc_ini

   SUBROUTINE apply_atm_fluxes
      !!----------------------------------------------------------------------
      !! ** Purpose : Change model variables according to OB conditions
      !  REMOVED: Relaxation with the time scale nm_obc_tau_relax towards a precribed OB profile.

      !!----------------------------------------------------------------------
      IMPLICIT NONE

!      integer :: i, j !
!      integer               :: obc_alloc
!      integer               :: d_indx, d_indx_p1  ! index of nearest (but less) depth in OB data
!      real(wp)              :: data1, data2, cf_a, cf_b, delta_d


      eta_n_obc = obcdata(i_ssh,:)


   END SUBROUTINE apply_atm_fluxes


   SUBROUTINE getcoeffld(rdate)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE getcoeffld ***
      !!
      !! ** Purpose : read fields from files, prepare interpolation coeffients (in time)
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
      integer              :: obc_alloc

      real(wp)             :: denom, x1, x2, y1, y2, x, y
      real(wp)             :: now_date

      real(wp), allocatable, dimension(:)  :: sbcdata1,sbcdata2
      real(wp)             :: data1,data2
      real(wp)             :: delta_t   ! time(t_indx) - time(t_indx+1)

      integer              :: elnodes(4) !4 nodes from one element
      integer              :: numnodes   ! nu,ber of nodes in elem (3 for triangle, 4 for ... )

      ALLOCATE( sbcdata1(nc_Nnode), sbcdata2(nc_Nnode),&
                &      STAT=obc_alloc )
!                data1(elem2D),data2(elem2D), &
      if( obc_alloc /= 0 )   STOP 'getcoeffld: failed to allocate arrays'
      ! find time index in files
      now_date = rdate
      call binarysearch(nc_Ntime,nc_time,now_date,t_indx)
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

      if ( now_date < nc_time(1) ) then  ! NO extrapolation back in time
         t_indx = 1
         t_indx_p1 = t_indx
         delta_t = 1.0_wp
      end if
      do fld_idx = 1, i_totc
         !open file obc_flfi
         iost = nf_open(obc_flfi(fld_idx)%file_name,NF_NOWRITE,ncid)
         call check_nferr(iost,obc_flfi(fld_idx)%file_name)
         ! get variable id
         iost = nf_inq_varid(ncid, obc_flfi(fld_idx)%var_name, id_data)
         call check_nferr(iost,obc_flfi(fld_idx)%file_name)
         !read data from file
         nf_start(1)=1
         nf_edges(1)=nc_Nnode
         nf_start(2)=t_indx
         nf_edges(2)=1
         iost = nf_get_vara_double(ncid, id_data, nf_start, nf_edges, sbcdata1)
         call check_nferr(iost,obc_flfi(fld_idx)%file_name)
         ! read next time step in file (check for +1 done before)
         nf_start(2)=t_indx_p1
         nf_edges(2)=1
         iost = nf_get_vara_double(ncid, id_data, nf_start, nf_edges, sbcdata2)
         call check_nferr(iost,obc_flfi(fld_idx)%file_name)

         ! bilinear time interpolation ,
         ! data is assumed to be sampled on a regular grid
         ! calculate new coefficients for interpolations
         ! delta_t was calculated before, in case of extrapolation delta_t=1
         coef_a(fld_idx, :) = ( sbcdata2 - sbcdata1 ) / delta_t !( nc_time(t_indx+1) - nc_time(t_indx) )
         coef_b(fld_idx, :) = sbcdata1 - coef_a(fld_idx, :) * nc_time(t_indx)

         iost = nf_close(ncid)
         call check_nferr(iost,obc_flfi(fld_idx)%file_name)

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
      real(wp),intent(in) :: rdate ! date
!      integer,intent(in) :: nsec  ! seconds
     ! assign data from interpolation to taux and tauy
      integer            :: fld_idx, i,j,ii
      real(wp)           :: now_date

      now_date = rdate

      obcdata = now_date * coef_a + coef_b

   END SUBROUTINE data_timeinterp

   SUBROUTINE obc_2d_ini
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE obc_2d_ini ***
      !!
      !! ** Purpose : inizialization of ocean OB 2D
      !! ** Method  :
      !! ** Action  :
      !!----------------------------------------------------------------------
      IMPLICIT NONE

      integer            :: idate ! initialization date
      integer            :: isec  ! initialization seconds
      integer            :: iost  ! I/O status
      integer            :: obc_alloc                   !: allocation status
      integer            :: i,itmp     ! temporary integer



      namelist/nam_obc_2d/ nm_sshobc_file, nm_obc_2d, &
                        nm_sshobc_var, &
                        nm_obc_tau_relax,&
                        nm_nc_iyear, nm_nc_imm, nm_nc_idd, nm_nc_secstep

      write(*,*) "Start: 2D Open Bounday inizialization."

      idate = time_jd0
      isec  = 86400.0 *(time_jd0 - idate)

      ! OPEN and read namelist for OBC
      open( unit=nm_obc_unit, file='namelist_bc.nml', form='formatted', access='sequential', status='old', iostat=iost )
      if( iost == 0 ) then
         WRITE(*,*) '     file   : ', 'namelist_bc.nml',' open ok'
      else
         WRITE(*,*) 'ERROR: --> bad opening file   : ', 'namelist_bc.nml',' ; iostat=',iost
         STOP 'ERROR: --> obc_2d_ini'
      endif
      READ( nm_obc_unit, nml=nam_obc_2d, iostat=iost )
      close( unit=nm_obc_unit)

      if( nm_obc_2d == -1 ) then
         write(*,*) 'obc_2d_ini: nm_obc_2d = -1; no 2D open boundary'
         return
      endif

      !get total number of nodes at open boundary
      nc_Nnode = 0
      do i = 1, nod2D
         if (index_nod2D(i) == 2) nc_Nnode = nc_Nnode + 1
      enddo

      if( nc_Nnode == 0 ) then
         write(*,*) 'obc_2d_ini: NO nodes at open boundary. switch nm_obc_2d to -1 '
         nm_obc_2d = -1
         return
      endif

      ALLOCATE( index_nod2D_ob(nc_Nnode),&
                &      STAT=obc_alloc )
      if( obc_alloc /= 0 )   STOP 'obc_ini: failed to allocate arrays'
      ! write indexs of OB nodes to index_nod2D_ob
      itmp = 0
      do i = 1, nod2D
         if (index_nod2D(i) == 2) then
            itmp = itmp + 1
            index_nod2D_ob(itmp) = i
         endif
      enddo
!    allocation done after reading nc file

      call nc_obc_ini(time_jd0)



      write(*,*) "DONE:  Open Boundary inizialization."
   END SUBROUTINE obc_2d_ini





   SUBROUTINE obc_2d_do(rdate)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE obc_do ***
      !!
      !! ** Purpose : interpolate OB (open boundary) from IN file on model Z and apply
      !! ** Method  :
      !! ** Action  :
      !!----------------------------------------------------------------------
      IMPLICIT NONE

      real(wp),intent(in)           :: rdate ! date

      if( nm_obc_2d == -1 ) then
         ! do not apply open boundary
         return
      endif

      if( nc_Nnode == 0 ) then
         ! NO nodes at open boundary
         return
      endif


!      rdate = time_jd

      if( .not. one_field ) then
         ! IF NOT one time step in IN files
         !assuming that OBC were changed during previews time
         if( rdate > nc_time(t_indx_p1) ) then
            !if rdate later than second date in forcing
            ! get new coeficients for time inerpolation on model grid for all datas
            call getcoeffld(rdate)
         endif

         ! interpolate in time
         call data_timeinterp(rdate)

      endif


      ! change model ocean parameteres
      call apply_atm_fluxes




   END SUBROUTINE obc_2d_do




   SUBROUTINE err_call(iost,fname)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE  err_call ***
      !!
      !! ** Purpose : call Error
      !! ** Method  :
      !! ** Action  :
      !!----------------------------------------------------------------------
   IMPLICIT NONE
      integer, intent(in)            :: iost
      character(len=256), intent(in) :: fname
      write(*,*) 'ERROR: I/O status=',iost,' file= ',fname
      STOP 'ERROR:  stop'


   END SUBROUTINE err_call

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
      if (dd+31*(mm+12*yyyy) >= IGREG) then
         ja=int(0.01*jy)
         julday=julday+2-ja+int(0.25_wp*ja)
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


   SUBROUTINE obc_2d_end

      IMPLICIT NONE


      DEALLOCATE( nc_lon, nc_lat, &
                  nc_time, nc_node, &
                  coef_a, coef_b, &
                  obcdata, index_nod2D_ob,&
                  eta_n_obc)



   END SUBROUTINE obc_2d_end

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
         middle = nint((left+right) / 2.0)
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

END MODULE fv_obc_2d
