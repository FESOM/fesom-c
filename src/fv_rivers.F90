MODULE fv_rivers
   !!===========================================================================
   !! Rivers:
   !!===========================================================================
   !! History: v1.0 , 10/2021, Ivan Kuznetsov (AWI),   MPI version
   !!
   !! public:
   !!   rivers_ini  -- initialization of rivers
   !!   rivers_do   -- change of rivers runoff at each time step
   !!
   USE o_ARRAYS
   USE o_MESH
   USE o_PARAM

   USE g_parsup

   IMPLICIT NONE

   !include 'netcdf.inc' do not use netcdf in this version

   public  rivers_ini  ! routine called before 1st time step (open files, read namelist,...)
   public  rivers_do   ! routine called each time step to provide a obc fileds
   public  rivers_end  ! routine called after last time step


   private

   ! namelists
   integer, save  :: rivers_unit     = 209       ! unit to open data file
  !============== namelistatmdata variables ================
   character(len=256), save   :: rivers_file = 'rivers.dat' ! name of file with temperature and salinity at OB (open boundary)

   ! ========== interpolation coefficients

   real(wp), allocatable, save, dimension(:)   :: coef_b ! time inerp coef. b (x=a*t+b) (node)
   real(wp), allocatable, save, dimension(:)   :: coef_a ! time inerp coef. a (x=a*t+b) (node)

   logical, save :: one_field = .false.! only one time step used for rivers

   real(wp), allocatable, save, dimension(:,:) :: my_data_rivers ! rivers runoff data for all time step (local PE)
   real(wp), allocatable, save, dimension(:)   :: time_rivers ! rivers times in INfile for all time step
! now in fv_var.f90   integer,  allocatable, save, dimension(:)   :: my_nod_rivers ! indexes of rivers nods (local PE (my))

  ! length of arrays
! now in fv_var.f90   integer,save              :: my_Nnods_rivers !(local)
   integer,save              :: Ntime_rivers
   ! time index for INfiles time array
   integer,save              :: t_indx    ! now time index in INfile time array
   integer,save              :: t_indx_p1 ! now time index +1 in INfile time array



!============== NETCDF ==========================================
CONTAINS

SUBROUTINE rivers_ini
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE obc_ini ***
      !!
      !! ** Purpose : initialization of ocean rivers
      !! ** Method  :
      !!       1. make mapping for MPI
      !!       2. read data from ascii file by pe=0 and distribute to others PE
      !!           ascii file as following:
      !!              Ntime_rivers - number of records (time) (can be 1 for constant runoff)
      !!              Nnods_rivers - number of nodes with runoff
      !!              time_rivers  - array of time step in julian days
      !!              nodes_rivers   - array of nodes indexes (integer)
      !!              data_rivers(:,:) - array with runoff, record by record
      !!                                 line with runoff for 1 time record, .... ,... Ntime_rivers time record
      !!       3. on each PE choose only there on data (nodes, runoff)
      !! ** Action  :
      !!----------------------------------------------------------------------
      IMPLICIT NONE

      integer            :: iost  ! I/O status
      integer            :: err_alloc                   !: allocation status
      integer            :: n,i,itmp,j    ! temporary integer

      integer,  allocatable, dimension(:)     :: nodes_rivers ! indexes of rivers nods   (global)
      real(wp), allocatable, dimension(:,:)   :: data_rivers ! runoff (time,nods)   (global)
      real(wp), allocatable, dimension(:)   :: data_rivers_tmp ! runoff (time,nods)    (global)
      integer                                 :: Nnods_rivers !(global)

      integer, allocatable, dimension(:)    :: mapping


      if (mype==0) write(*,*) "= mype= 0 ================================================"
      if (mype==0) write(*,*) "============= Start: Rivers initialization. =============="
      if (mype==0) write(*,*) "=========================================================="

      rivers_file = trim(meshpath)//'rivers.dat'
#ifdef USE_MPI
      allocate(mapping(nod2D),STAT=err_alloc )
      if( err_alloc /= 0 )   STOP 'failed to allocate arrays'
      mapping(:)=0
      do n=1, myDim_nod2D+eDim_nod2D
        mapping(myList_nod2D(n))=n
      end do
#endif
      !read rivers data
      if (mype==0) then
          open( unit=rivers_unit, FILE=rivers_file, status='old', iostat=iost )
          if( iost == 0 ) then
             write(*,*) '     file   : ', trim(rivers_file),' open ok'
          else
             write(*,*) 'ERROR: --> bad opening file   : ', trim(rivers_file),' ; iostat=',iost
             STOP 'ERROR: --> sbc_ini'
          endif
          !read rivers Ntime and Nnods
          read( unit=rivers_unit, FMT=*,iostat=iost) Ntime_rivers
          if( iost /= 0 ) STOP 'Error while reading rivers'
          read( unit=rivers_unit, FMT=*,iostat=iost) Nnods_rivers
          if( iost /= 0 ) STOP 'Error while reading rivers'
      endif

#ifdef USE_MPI
      call MPI_BCast(Ntime_rivers, 1, MPI_INTEGER, 0, MPI_COMM_FESOM_C, MPIerr)
      call MPI_BCast(Nnods_rivers, 1, MPI_INTEGER, 0, MPI_COMM_FESOM_C, MPIerr)
#endif
      !allocate local array for all rivers data
      allocate( time_rivers(Ntime_rivers), data_rivers(Ntime_rivers,Nnods_rivers), nodes_rivers(Nnods_rivers),&
                & data_rivers_tmp(Nnods_rivers), STAT=err_alloc )
      if( err_alloc /= 0 )   STOP 'Failed to allocate arrays'
      time_rivers = 0
      data_rivers = 0
      nodes_rivers = 0
      if (mype==0) then
          !read rivers time, nods, runoff
          read( unit=rivers_unit, FMT=*,iostat=iost) time_rivers
          if( iost /= 0 ) STOP 'Error while reading rivers'
          read( unit=rivers_unit, FMT=*,iostat=iost) nodes_rivers
          if( iost /= 0 ) STOP 'Error while reading rivers'
          do i = 1, Ntime_rivers
            read( unit=rivers_unit, FMT=*,iostat=iost) data_rivers_tmp
            data_rivers(i,:)=data_rivers_tmp
            if( iost /= 0 ) STOP 'Error while reading rivers'
          end do
          close( rivers_unit )
      endif

      ! send data to all PEs
#ifdef USE_MPI
      call MPI_BCast(nodes_rivers, Nnods_rivers, MPI_INTEGER, 0, MPI_COMM_FESOM_C, MPIerr)
      call MPI_BCast(time_rivers, Ntime_rivers, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, MPIerr)
      do i = 1, Ntime_rivers
         data_rivers_tmp = data_rivers(i,:)
         call MPI_BCast(data_rivers_tmp, Nnods_rivers, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, MPIerr)
         data_rivers(i,:) = data_rivers_tmp
         !call MPI_BCast(data_rivers(i,:), Nnods_rivers, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, MPIerr)
      enddo
      ! find number of river nodes at current PE
      my_Nnods_rivers = 0
      do i = 1, Nnods_rivers
        itmp = nodes_rivers(i)
        if (mapping(itmp) > 0) then
           my_Nnods_rivers = my_Nnods_rivers +1
        end if
      end do

      if (my_Nnods_rivers>0) then
          ! fill my_nod_rivers and my_data_rivers at current PE by glob data
          allocate( my_data_rivers(Ntime_rivers,my_Nnods_rivers), my_nod_rivers(my_Nnods_rivers),&
                    coef_a(my_Nnods_rivers), coef_b(my_Nnods_rivers), runoff_rivers(my_Nnods_rivers), &
                        &      STAT=err_alloc )
          if( err_alloc /= 0 )   STOP 'Failed to allocate arrays'
          my_nod_rivers(:) = 0
          my_data_rivers(:,:) = 0.0_WP
          coef_a(:) = 0.0_WP
          coef_b(:) = 0.0_WP
          runoff_rivers(:) = 0.0_WP
          j=0
          do i = 1, Nnods_rivers
            itmp = nodes_rivers(i)
            if (mapping(itmp)>0) then
               j = j +1
               my_nod_rivers(j) = mapping(itmp)
               my_data_rivers(:,j) = data_rivers(:,i)
            end if
          end do
      endif
#else
      ! in case serial, copy globe to my_ arrays
      my_Nnods_rivers = Nnods_rivers
      allocate( my_data_rivers(Ntime_rivers,my_Nnods_rivers), my_nod_rivers(my_Nnods_rivers),&
                coef_a(my_Nnods_rivers), coef_b(my_Nnods_rivers), runoff_rivers(my_Nnods_rivers),&
                &      STAT=err_alloc )
      my_nod_rivers = nodes_rivers
      my_data_rivers = data_rivers
      coef_a(:) = 0.0_WP
      coef_b(:) = 0.0_WP
      runoff_rivers(:) = 0.0_WP
#endif

      deallocate(mapping,data_rivers,nodes_rivers,&
                       STAT=err_alloc )
      if( err_alloc /= 0 )   STOP 'Failed to deallocate arrays'

      if (my_Nnods_rivers>0) then
         ! get first coefficients (coef_a,coef_b) for time interpolation
         call getcoeffld(time_jd)
         ! interpolate in time
         runoff_rivers = time_jd * coef_a + coef_b
         runoff_rivers = 0.0_WP !at first time step
      endif

      !write(*,*) "rivers: ",mype,my_Nnods_rivers,Ntime_rivers,my_nod_rivers

      if (mype==0) write(*,*) "DONE:  Rivers initialization. ============================"
   END SUBROUTINE rivers_ini

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

      integer              :: i
      real(wp)             :: denom, x1, x2, y1, y2, x, y
      real(wp)             :: now_date
      real(wp)             :: delta_t   ! time(t_indx) - time(t_indx+1)

      integer              :: elnodes(4) !4 nodes from one element
      integer              :: numnodes   ! number of nodes in elem (3 for triangle, 4 for ... )

      ! find time index of now time
      now_date = rdate
      call binarysearch(Ntime_rivers,time_rivers,now_date,t_indx)
      if ( t_indx < Ntime_rivers ) then
         t_indx_p1 = t_indx + 1
         delta_t = time_rivers(t_indx_p1) - time_rivers(t_indx)
      else ! NO extrapolation to future
         t_indx_p1 = t_indx
         delta_t = 1.0_wp
         if (mype==0) then
            write(*,*) "     WARNING:  time lager than last time step in rivers file,"
            write(*,*) "        last time step will be used as a constant runoff for rivers"
         endif
         one_field = .true.
      end if

      if ( now_date < time_rivers(1) ) then  ! NO extrapolation back in time
         t_indx_p1 = t_indx
         delta_t = 1.0_wp
      end if

      ! linear time interpolation ,
      ! calculate new coefficients for interpolations
      ! delta_t was calculated before, in case of extrapolation delta_t=1
      ! coef_a,coef_b,
      coef_a(:) = (my_data_rivers(t_indx_p1,:) - my_data_rivers(t_indx,:)) / delta_t
      coef_b(:) = my_data_rivers(t_indx,:) - coef_a(:) * time_rivers(t_indx)


   END SUBROUTINE getcoeffld


   SUBROUTINE rivers_do
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE rivers_do ***
      !!
      !! ** Purpose : interpolate  rivers runoff from data in memory on local PE nodes and time
      !! ** Method  :
      !! ** Action  :
      !!----------------------------------------------------------------------
      use o_PARAM
      IMPLICIT NONE

      real(wp)            :: rdate ! date
!      integer            :: rsec  ! seconds

      if (my_Nnods_rivers==0) then
        return
      end if

      rdate = time_jd

      if (.not. one_field ) then
         ! IF NOT one  time record in IN files

         if ( rdate > time_rivers(t_indx_p1) ) then
            !if rdate later than second date in forcing
            ! get new coeficients for time inerpolation on model grid for all datas
            call getcoeffld(rdate)
         endif

         ! interpolate in time
         runoff_rivers = rdate * coef_a + coef_b
         if (time_jd-time_jd0<10.0_WP) then
           runoff_rivers=runoff_rivers*(time_jd-time_jd0)/10.0_WP
         end if
      endif

   END SUBROUTINE rivers_do

   SUBROUTINE rivers_end

      IMPLICIT NONE

      !deallocate()

   END SUBROUTINE rivers_end

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

END MODULE fv_rivers
