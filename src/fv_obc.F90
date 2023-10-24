MODULE fv_obc
   !!===========================================================================
   !! Open boundary conditions:
   !!===========================================================================
   !! History: 0.1 ! 02/2016 I. Kuznetsov
   !!          1.01! 08/2020 I. Kuznetsov
   !!                              MPI version, based on netcdf libs
   !!
   !! public:
   !!   obc_ini  -- inizialization of open boundary
   !!   obc_do   -- provide a obc (open boundary conditions) each time step
   !!
   USE o_ARRAYS
   USE o_MESH
   USE o_PARAM

   USE g_parsup

   IMPLICIT NONE

   include 'netcdf.inc'

   public  obc_ini  ! routine called before 1st time step (open files, read namelist,...)
   public  obc_do   ! routine called each time step to provide a obc fileds
   public  obc_end  ! routine called after last time step


   private


   ! namelists
   integer, save  :: nm_obc_unit     = 201       ! unit to open namelist file
  !============== namelistatmdata variables ================
   character(len=256), save   :: nm_tsobc_file = 'tsobc.' ! name of file with temperature and salinity at OB (open boundary)

   character(len=34), save   :: nm_tobc_var   = 'temp' ! name of variable in file with temperature
   character(len=34), save   :: nm_sobc_var   = 'salt' ! name of variable in file with salinity

   integer, save  :: nm_obc   ! switch of module if -1
   ! ========== netCDF time param
   integer, save :: nm_nc_iyear = 1948    ! initial year of time axis in netCDF (1948 like CoastDat,1800 NCEP)
   integer, save :: nm_nc_imm = 1         ! initial month of time axis in netCDF
   integer, save :: nm_nc_idd = 1         ! initial day of time axis in netCDF
   real,    save :: nm_nc_secstep = 86400.0 ! time units coef (86400 CoastDat, 24 NCEP)

   ! ========== interpolation coeficients

   real(wp), allocatable, save, dimension(:,:,:)   :: coef_b ! time inerp coef. b (x=a*t+b) (fld_indx,node,depth)
   real(wp), allocatable, save, dimension(:,:,:)   :: coef_a ! time inerp coef. a (x=a*t+b) (fld_indx,node,depth)

   logical, save :: one_field = .false.! only one time step used for openboundary

   real(wp), allocatable, save, dimension(:,:,:)   :: obcdata ! OB data for current time step
   integer,  allocatable, save, dimension(:)       :: index_nod2D_ob ! indexes of nods at OB

   real(wp), allocatable, save, dimension(:)       :: intdata ! interpolated OB data for current time step and node (call apply_obc)


!============== NETCDF ==========================================
!structure similar to sbc (so fluxes here is concentrations)

   integer, parameter :: i_totc = 2 ! total number of 3D (temperature,...) values (2=T,S)
   integer, parameter :: i_temp = 1 ! index of temperature [grad C]
   integer, parameter :: i_salt = 2 ! index of salinity []


   type, public ::   flfi_type    !flux file informations
      character(len = 256) :: file_name ! file name
      character(len = 34)  :: var_name  ! variable name in the NetCDF file
   end type flfi_type

  type(flfi_type),save, dimension(i_totc) :: obc_flfi  !array for information about flux files

  ! arrays of time, lon and lat in INfiles
   real(wp), allocatable, save, dimension(:)  :: nc_node ! index of nodes in infile corespond to node in mesh
   real(wp), allocatable, save, dimension(:)  :: nc_time ! time array
   real(wp), allocatable, save, dimension(:)  :: nc_depth! depth (z) in netcdf file [m] , assume we have z-coordinate system
  ! lenght of arrays in INfiles
   integer,save              :: nc_Ndepth
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
      integer              :: id_depth
      integer              :: id_node
      integer              :: id_time
      integer              :: id_lon
      integer              :: id_lat
      integer              :: id_noded
      integer              :: id_timed
      integer              :: id_depthd
!      integer              :: nf_dims(4) ! dimensions (temporal)
      integer              :: nf_start(4)
      integer              :: nf_edges(4)
      integer              :: zero_year,yyyy,mm,dd
      character(len = 256) :: att_string ! attribute
      integer              :: obc_alloc                   !: allocation status


      if (mype==0) then

          !open file
          iost = nf_open(flf%file_name,NF_NOWRITE,ncid)
          call check_nferr(iost,flf%file_name)

          ! get dimensions
          iost = nf_inq_dimid(ncid, "node", id_noded)
          call check_nferr(iost,flf%file_name)
          iost = nf_inq_dimid(ncid, "depth", id_depthd)
          call check_nferr(iost,flf%file_name)
          iost = nf_inq_dimid(ncid, "time", id_timed)
          call check_nferr(iost,flf%file_name)

          ! get variable id
          iost = nf_inq_varid(ncid, "node", id_node)
          call check_nferr(iost,flf%file_name)
          iost = nf_inq_varid(ncid, "depth", id_depth)
          call check_nferr(iost,flf%file_name)
          iost = nf_inq_varid(ncid, "time", id_time)
          call check_nferr(iost,flf%file_name)
          iost = nf_inq_varid(ncid, "lon", id_lon)
          call check_nferr(iost,flf%file_name)
          iost = nf_inq_varid(ncid, "lat", id_lat)
          call check_nferr(iost,flf%file_name)
          !  get dimensions size
          iost = nf_inq_dimlen(ncid, id_depthd, nc_Ndepth)
          call check_nferr(iost,flf%file_name)
          iost = nf_inq_dimlen(ncid, id_noded, i)
          call check_nferr(iost,flf%file_name)
          if( i /= nobn ) STOP 'nc_readTimeGrid: size of OB nodes in file .neq. to mesh'
          iost = nf_inq_dimlen(ncid, id_timed, nc_Ntime)
          call check_nferr(iost,flf%file_name)
      endif
#ifdef USE_MPI
      call MPI_BCast(nc_Ntime, 1, MPI_INTEGER, 0, MPI_COMM_FESOM_C, MPIerr)
      call MPI_BCast(nc_Ndepth, 1, MPI_INTEGER, 0, MPI_COMM_FESOM_C, MPIerr)
#endif

      ALLOCATE( nc_time(nc_Ntime), nc_node(nobn), &
                nc_depth(nc_Ndepth), &
                &      STAT=obc_alloc )
      if( obc_alloc /= 0 )   STOP 'nc_readTimeGrid: failed to allocate arrays'
 !(fld_indx,node,depth)
      ALLOCATE( coef_a(i_totc,nobn,nc_Ndepth), coef_b(i_totc,nobn,nc_Ndepth), &
              & obcdata(i_totc,nobn,nc_Ndepth), &
                   &      STAT=obc_alloc )
      if( obc_alloc /= 0 )   STOP 'nc_readTimeGrid: failed to allocate arrays'
      ALLOCATE( intdata(nsigma-1),&
                &      STAT=obc_alloc )
      if( obc_alloc /= 0 )   STOP 'nc_readTimeGrid: failed to allocate arrays'

      if (mype==0) then
       !read variables from file
       !time
          nf_start(1)=1
          nf_edges(1)=nc_Ntime
          iost = nf_get_vara_double(ncid, id_time, nf_start, nf_edges, nc_time)
          call check_nferr(iost,flf%file_name)
          ! convert time to days
          !!NCEP     nc_time = nc_time / 24.0 + julday(1800,1,1)
    !      write(*,*) 'nc_time(1)=',nc_time(1)
          nc_time = nc_time / nm_nc_secstep + julday(nm_nc_iyear,nm_nc_imm,nm_nc_idd) - 0.5
    !   write(*,*) 'nc_time(1)=',nc_time(1),nm_nc_secstep,nm_nc_iyear,nm_nc_imm,nm_nc_idd,julday(nm_nc_iyear,nm_nc_imm,nm_nc_idd)
       ! coordinates
          nf_start(1)=1
          nf_edges(1)=nobn
          iost = nf_get_vara_double(ncid, id_node, nf_start, nf_edges,nc_node)
          call check_nferr(iost,flf%file_name)
       ! depth
          nf_start(1)=1
          nf_edges(1)=nc_Ndepth
          iost = nf_get_vara_double(ncid, id_depth, nf_start, nf_edges,nc_depth)
          call check_nferr(iost,flf%file_name)

          iost = nf_close(ncid)
          call check_nferr(iost,flf%file_name)
          warn = 0
write(*,*) "nc_node:"          ,nc_node
write(*,*) "in_obn:"          ,in_obn

          do i = 1, nobn
             if ( nc_node(i) /= in_obn(i) ) then
                if ( warn == 0 ) then
                   write(*,*) 'WARNING:: OB nodes in Netcdf file are not the same as in mesh. node:',nc_node(i)
                   warn = 1
                endif
             endif
          enddo
      endif! (mype==0)
#ifdef USE_MPI
      call MPI_BCast(nc_time, nc_Ntime, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, MPIerr)
      call MPI_BCast(nc_depth, nc_Ndepth, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, MPIerr)
#endif
      DEALLOCATE( nc_node)


   END SUBROUTINE nc_readTimeGrid

   SUBROUTINE nc_obc_ini_fillnames(yyear)
      character(len=4),intent(in)   :: yyear

      !! ** Purpose : Fill names of obc_flfi array (file names and variable names)

      !prepare proper nc file (add year and .nc to the end of the file name from namelist
      write(obc_flfi(i_temp)%file_name,*) trim(nm_tsobc_file)!,yyear,'.nc'
      write(obc_flfi(i_salt)%file_name,*) trim(nm_tsobc_file)!,yyear,'.nc'

      obc_flfi(i_temp)%file_name=ADJUSTL(trim(obc_flfi(i_temp)%file_name))
      obc_flfi(i_salt)%file_name=ADJUSTL(trim(obc_flfi(i_salt)%file_name))

      obc_flfi(i_temp)%var_name=ADJUSTL(trim(nm_tobc_var))
      obc_flfi(i_salt)%var_name=ADJUSTL(trim(nm_sobc_var))
   END SUBROUTINE nc_obc_ini_fillnames

   SUBROUTINE nc_obc_ini(rdate)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE nc_obc_ini ***
      !!
      !! ** Purpose : inizialization of OB from NETCDF file
      !!----------------------------------------------------------------------

      IMPLICIT NONE
      integer :: idate ! initialization date
      real(wp),intent(in) :: rdate ! initialization seconds

      character(len=4)   :: yyear
      integer            :: yyyy,mm,dd

      integer            :: i
      integer            :: obc_alloc

      integer            :: elnodes(4) !4 nodes from one element






      ! get ini year; Fill names of obc_flfi
      idate = rdate
      call calendar_date(idate,yyyy,dd,mm)
      write(yyear,"(I4)") yyyy
      call nc_obc_ini_fillnames(yyear)

      ! we assume that all NetCDF files have identical grid and time variable
      ! read axis, sizes, allocate memory
      call nc_readTimeGrid(obc_flfi(i_temp))


!write(*,*) 'idate,isec',idate,isec
      ! get first coeficients for time inerpolation on model grid for all datas
      call getcoeffld(rdate)

      ! interpolate in time

      call data_timeinterp(rdate)

!      call apply_obc
!     no OB apply before first time step


   END SUBROUTINE nc_obc_ini

   SUBROUTINE apply_obc
      !!----------------------------------------------------------------------
      !! ** Purpose : Change model variables according to OB conditions
      !!----------------------------------------------------------------------
      IMPLICIT NONE

      integer :: i, j !
      integer               :: obc_alloc
      integer               :: d_indx, d_indx_p1  ! index of nearest (but less) depth in OB data
      real(wp)              :: data1, data2, cf_a, cf_b, delta_d



      do i = 1, mynobn ! go over all OB nodes saved in index_nod2D_ob
         ! temperature
         ! interpolate vertical profile from OB obcdata
         !   on curent Z (depth) of node (new values for intdata)
         call profile_interp(i_temp,i)
         ! We change Tclim ans Sclim, so adding to T and S will be done in fv_tracer by cl_relax key
         Tclim(:,my_in_obn(i)) = intdata(:)
         ! salinity
         call profile_interp(i_salt,i)
         Sclim(:,my_in_obn(i)) = intdata(:)

      enddo


   END SUBROUTINE apply_obc


   SUBROUTINE profile_interp(fld_idx,nod_idx)
      !!----------------------------------------------------------------------
      !! ** Purpose : profile interpolation  interpolate data from obcdata on Z and write to intdata
      !! output: intdata
      !!----------------------------------------------------------------------
      IMPLICIT NONE
      integer, intent(in)   :: fld_idx ! index of field to interpolate
      integer, intent(in)   :: nod_idx ! index of nod to interpolate
      integer :: i
      integer               :: d_indx, d_indx_p1  ! index of nearest (but less) depth in OB data
      real(wp)              :: data1, data2, cf_a, cf_b, delta_d


      do i= 1, nsigma-1
         call binarysearch(nc_Ndepth,nc_depth,Z(i,my_in_obn(nod_idx)),d_indx)
         if ( d_indx < nc_Ndepth ) then
            d_indx_p1 = d_indx+1
            delta_d = nc_depth(d_indx+1)-nc_depth(d_indx)
         else
            ! extrapolation
            d_indx_p1 = d_indx
            delta_d = 1.0_wp
         endif
         ! values from OB data for nearest depth
         data1 = obcdata(fld_idx,nod_idx,d_indx)
         data2 = obcdata(fld_idx,nod_idx,d_indx_p1)
         ! line a*z+b coefficients calculation
         cf_a  = (data2 - data1)/ delta_d
         ! value of interpolated OB data on Z from model
         cf_b  = data1 - cf_a * nc_depth(d_indx)
         intdata(i) = cf_a * Z(i,my_in_obn(nod_idx)) + cf_b


      enddo

   END SUBROUTINE profile_interp

   SUBROUTINE getcoeffld(rdate)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE getcoeffld ***
      !!
      !! ** Purpose : read fields from files, prepare interpolation coeffients (in time)
      !! ** Method  :
      !! ** Action  :
      !!----------------------------------------------------------------------
      IMPLICIT NONE
      real(wp),intent(in)   :: rdate ! initialization date
!      integer,intent(in)   :: isec ! initialization seconds

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

      real(wp), allocatable, dimension(:,:)  :: sbcdata1,sbcdata2
      real(wp)             :: data1,data2
      real(wp)             :: delta_t   ! time(t_indx) - time(t_indx+1)

      integer              :: elnodes(4) !4 nodes from one element
      integer              :: numnodes   ! nu,ber of nodes in elem (3 for triangle, 4 for ... )

      ALLOCATE( sbcdata1(nobn,nc_Ndepth), sbcdata2(nobn,nc_Ndepth),&
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
         if (mype==0) then
            write(*,*) "     WARNING:  time lager than last time step in forcing file,"
            write(*,*) "        last time step will be used as a constant field"
            write(*,*) "now_date=",now_date,"t_indx=",t_indx,"nc_Ntime=",nc_Ntime,"nc_time(t_indx)=",nc_time(t_indx)
         endif
         one_field = .true.
      end if

      if ( now_date < nc_time(1) ) then  ! NO extrapolation back in time
         t_indx_p1 = t_indx
         delta_t = 1.0_wp
      end if

      do fld_idx = 1, i_totc
         if (mype==0) then
             !open file obc_flfi
             iost = nf_open(obc_flfi(fld_idx)%file_name,NF_NOWRITE,ncid)
             call check_nferr(iost,obc_flfi(fld_idx)%file_name)
             ! get variable id
             iost = nf_inq_varid(ncid, obc_flfi(fld_idx)%var_name, id_data)
             call check_nferr(iost,obc_flfi(fld_idx)%file_name)
             !read data from file
             nf_start(1)=1
             nf_edges(1)=nobn
             nf_start(2)=1
             nf_edges(2)=nc_Ndepth
             nf_start(3)=t_indx
             nf_edges(3)=1
             iost = nf_get_vara_double(ncid, id_data, nf_start, nf_edges, sbcdata1)
             call check_nferr(iost,obc_flfi(fld_idx)%file_name)
             ! read next time step in file (check for +1 done before)
             nf_start(3)=t_indx_p1
             nf_edges(3)=1
             iost = nf_get_vara_double(ncid, id_data, nf_start, nf_edges, sbcdata2)
             call check_nferr(iost,obc_flfi(fld_idx)%file_name)

             iost = nf_close(ncid)
             call check_nferr(iost,obc_flfi(fld_idx)%file_name)
         endif

#ifdef USE_MPI
         call MPI_BCast(sbcdata1(1:nobn,1:nc_Ndepth), nobn*nc_Ndepth, &
                          MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, MPIerr)
         call MPI_BCast(sbcdata2(1:nobn,1:nc_Ndepth), nobn*nc_Ndepth, &
                          MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, MPIerr)

#endif
         ! bilinear time interpolation ,
         ! data is assumed to be sampled on a regular grid
         ! calculate new coeficients for interpolations
         ! delta_t was calculated before, in case of extrapolation delta_t=1
         coef_a(fld_idx, :,:) = ( sbcdata2 - sbcdata1 ) / delta_t !( nc_time(t_indx+1) - nc_time(t_indx) )
         coef_b(fld_idx, :,:) = sbcdata1 - coef_a(fld_idx, :,:) * nc_time(t_indx)


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

      now_date = rdate !+nsec/86400.0_wp

      obcdata = now_date * coef_a + coef_b

   END SUBROUTINE data_timeinterp

   SUBROUTINE obc_ini
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE obc_ini ***
      !!
      !! ** Purpose : inizialization of ocean OB
      !! ** Method  :
      !! ** Action  :
      !!----------------------------------------------------------------------
      IMPLICIT NONE

!      integer            :: rdate ! initialization date
!      integer            :: isec  ! initialization seconds
      integer            :: iost  ! I/O status
      integer            :: obc_alloc                   !: allocation status
      integer            :: i,itmp     ! temporary integer



      namelist/nam_obc/ nm_tsobc_file, nm_obc, &
                        nm_tobc_var, nm_sobc_var, &
                        nm_nc_iyear, nm_nc_imm, nm_nc_idd, nm_nc_secstep

      if (mype==0) write(*,*) "= mype= 0 ================================================"
      if (mype==0) write(*,*) "========= Start: Open Boundary initialization. ==========="
      if (mype==0) write(*,*) "=========================================================="

      if( nobn == 0 ) then
         write(*,*) 'obc_ini: NO nodes at open boundary. switch nm_obc to -1 '
         nm_obc = -1
         return
      endif
!      rdate = time_jd
!      isec  = 86400.0 *(time_jd - idate)

      ! OPEN and read namelist for OBC
      open( unit=nm_obc_unit, file='namelist_bc.nml', form='formatted', access='sequential', status='old', iostat=iost )
      if( iost == 0 ) then
         if (mype==0) WRITE(*,*) '     file   : ', 'namelist_bc.nml',' open ok'
      else
         WRITE(*,*) 'ERROR: --> bad opening file   : ', 'namelist_bc.nml',' ; iostat=',iost
         STOP 'ERROR: --> obc_ini'
      endif
      READ( nm_obc_unit, nml=nam_obc, iostat=iost )
      if( iost == 0 ) then
      else
         WRITE(*,*) 'ERROR: --> bad reading namelist nam_obc',' ; iostat=',iost
         STOP 'ERROR: --> obc_ini'
      endif
      close( unit=nm_obc_unit )

      if (mype==0) then
        write(*,*) "obc_ini: nm_obc_unit:"
        write(*,*) "                    : nm_obc = ",nm_obc
        write(*,*) "                    : nm_tsobc_file = ",nm_tsobc_file
        write(*,*) "                    : nm_tobc_var = ",nm_tobc_var
        write(*,*) "                    : nm_sobc_var = ",nm_sobc_var
        write(*,*) "                    : nm_nc_iyear = ",nm_nc_iyear
        write(*,*) "                    : nm_nc_imm = ",nm_nc_imm
        write(*,*) "                    : nm_nc_idd = ",nm_nc_idd
        write(*,*) "                    : nm_nc_secstep = ",nm_nc_secstep
      endif

      if( nm_obc == -1 ) then
         if (mype==0) write(*,*) "obc_ini: nm_obc = -1; no open boundary"
         return
      endif

      !get total number of nodes at open boundary
!      nc_Nnode = 0
!      do i = 1, nod2D
!         if (index_nod2D(i) == 2) nc_Nnode = nc_Nnode + 1
!      enddo

!      indexes are found in main



!      ALLOCATE( index_nod2D_ob(nc_Nnode),&
!                &      STAT=obc_alloc )
!      if( obc_alloc /= 0 )   STOP 'obc_ini: failed to allocate arrays'
      ! write indexs of OB nodes to index_nod2D_ob
!      itmp = 0
!      do i = 1, nod2D
!         if (index_nod2D(i) == 2) then
!            itmp = itmp + 1
!            index_nod2D_ob(itmp) = i
!         endif
!      enddo
!    allocation done after reading nc file

      call nc_obc_ini(time_jd)



      if (mype==0) write(*,*) "DONE:  Open Boundary initialization."
   END SUBROUTINE obc_ini





   SUBROUTINE obc_do
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE obc_do ***
      !!
      !! ** Purpose : interpolate OB (open boundary) from IN file on model Z and apply
      !! ** Method  :
      !! ** Action  :
      !!----------------------------------------------------------------------
      IMPLICIT NONE

      real(wp)            :: rdate ! date
!      integer            :: rsec  ! seconds

      if( nm_obc == -1 ) then
         ! do not apply open boundary
         return
      endif

      rdate = time_jd
!      rsec  = 86400.0 *(time_jd - rdate)
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


      ! change model ocean parameters (Tclim and Sclim)
      call apply_obc




   END SUBROUTINE obc_do




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


   SUBROUTINE obc_end

      IMPLICIT NONE


      DEALLOCATE( nc_time, &
                  nc_depth, coef_a, coef_b, &
                  obcdata, intdata, index_nod2D_ob)



   END SUBROUTINE obc_end

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

END MODULE fv_obc
