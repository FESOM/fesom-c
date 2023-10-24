MODULE fv_ncoutputmpi
#ifdef USE_MPI


   !!===========================================================================
   !! FVcoastal netcdf OUTPUT module:
   !!===========================================================================
   !!
   !! History: 0.1 ! 07/2015 I. Kuznetsov
   !!           basic configuration to write snapshot to nc file
   !!          0.9 ! 12/2015 I. Kuznetsov
   !!           new output structures, multiple files, averaging and snapshots, automatic control
   !!          0.95 ! 05/2016 I. Kuznetsov
   !!           improved output structures, now every file has possibility to save on specific nodes (stations, transaction, zoom field) (but cannot write yet zooms),
   !!           filelist(1)%ftype and filelist(1)%pos?_?? positions added to outputfilelist
   !!          0.99 ! 12/2018 I. Kuznetsov
   !!           bug fix: w_cv, x and y on elements
   !!           add:  additional info to mesh: nod_in_elem2d, nod_in_elem2d_num, depth_elem...
   !!           mesh writes correct on both Cartesian and non Cartesian coordinates
   !!           looking for nearest stations in both cases.
   !!           checked averaging function
   !!          2.01 ! 08/2020 I.K. MPI version.
   !!            mpi version with PnetCDF lib
   !!            https://parallel-netcdf.github.io
   !!            Jianwei Li, Wei-keng Liao, Alok Choudhary, Robert Ross, Rajeev Thakur, William Gropp, Rob Latham, Andrew Siegel, Brad Gallagher, and Michael Zingale.
   !!            Parallel netCDF: A Scientific High-Performance I/O Interface. In the Proceedings of ACM/IEEE conference on Supercomputing, pp. 39, November, 2003.
   !!
   !! Description:
   !!     * restart done every N steps, if run was from restart, then counter will start again.
   !!       to avoid problems with combining files use Restart N / output step  = int
   !!    first initialization of netcdf file before first time step should be done.
   !!    at each time step subroutine output_do should be called
   !!    Warning: in case of variable dt (time step) averaging will be incorrect
   !!     to control output file 'outputfilelist.dat' should be modified:
   !!     each file could be or snapshots or averaged values for fixed frequency in time steps.
   !!     for example: write to nc file "output100avg.nc" mean values of one variable(ilelist(1)%numvars): eta_n
   !!                  averaged for 100 time steps (filelist(1)%freq),
   !!                  in the list of variabes (filelist(1)%vars) should be indicated wich variabes to save:
   !!           !total number of files:
   !!           1
   !!           !filelist(1)%name=
   !!           output100avg.nc
   !!           !filelist(1)%snapshot=
   !!           0
   !!           !filelist(1)%ftype=
   !!           0
   !!           !filelist(1)%pos?_?? positions posx_ll, posy_ll, posx_ur, posy_ur
   !!           4.0, 53.0, 4.0, 53.0
   !!           !filelist(1)%freq=
   !!           100
   !!           !filelist(1)%specfreq=
   !!           0
   !!           !filelist(1)%numvars=
   !!           1
   !!           !filelist(1)%vars(1)=
   !!           1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
   !!
   !!
   !!  list of possible variabes is in  SUBROUTINE define_varlist
   !!  to add more variables subroutines define_varlist and one of
   !!  getvar0d/getvar1d/getvar2d (depens on number of variable dimensions) should be modified
   !!  add description and parameters of variable in define_varlist
   !!  and extand select case in getvar?d .
   !!  example how to extend output list with Jc(nsigma-1,nod2D) :
   !!    1. extend subroutines define_varlist
   !!          at the end of variabes list add:
   !!          nvar=nvar+1     ! increase number of possible variables
   !!          varlist(nvar)%name =           'temperature' ! name how it will be in netCDF file
   !!          varlist(nvar)%modelname =      'TF' ! name of variable in model
   !!          varlist(nvar)%standard_name =  'temperature'! full variable name in the NetCDF file
   !!          varlist(nvar)%units =          'grad C'! variable units in the NetCDF file
   !!          varlist(nvar)%onnode =         1 ! 1 for temperature while it is on node , if on nod =1 , if on elements =0, other =-1
   !!          varlist(nvar)%timedim =        1 ! 1 while it resolved in time, dimension 0 (not in time), 1 (in time)
   !!          varlist(nvar)%vertdim =        2 ! 2 while it is allocate(TF(nsigma-1,nod2D)) vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
   !!          varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
   !!          varlist(nvar)%varid =          nvar !8 ! ID of variable
   !!    2. TF is 3D variable, however x-y plane is only one dimension so for us it is 2d variable
   !!       extend getvar2d subroutine, add one more case were nvar should be from define_varlist (arlist(nvar)%varid = nvar):
   !!         case(8)
   !!            varbas2d(i)%snapshot =  TF
   !!    3. now it is possible to use 8 in !filelist(1)%vars(1)= from outputfilelist.dat

   !! public:
   !!   output_ini  -- inizialization of outfiles
   !!   output_do   -- call each time step

   USE o_MESH
   USE o_ARRAYS
   USE o_PARAM
!   USE fv_sbc
   USE g_parsup
   USE pnetcdf
   USE i_ARRAYS
   USE i_PARAM
   IMPLICIT NONE


   public  ncoutputmpi_ini    ! routine called before 1st time step (open files, read namelist,...)
   public  ncoutputmpi_do     ! routine called each time step

   private
      integer      :: iost !I/O status
      integer      :: ncid ! netcdf file ID for temporal
      integer,save :: ncid_m ! netcdf file ID for output (Main (snapshots))
      integer,save :: nc_rec ! record in NC file (=0 at inizialization)

      ! temporal arrays
      integer      :: nf_dims(4) ! dimensions (temporal)
      integer(kind=MPI_OFFSET_KIND)      :: nf_start(4)
      integer(kind=MPI_OFFSET_KIND)      :: nf_edges(4)
      type  ::   var_vect0    ! allocatable vector for data , used by var_type1d and var_type2d to save avg values
          real(kind=WP)                                :: v ! name
          real(kind=WP)                                :: sn ! name
      end type var_vect0
      type  ::   var_vect1    ! allocatable vector for data , used by var_type1d and var_type2d to save avg values
          real(kind=WP),allocatable,dimension(:)       :: v ! name
          real(kind=WP),allocatable,dimension(:)       :: sn ! name
          integer                                      :: sz
      end type var_vect1
      type  ::   var_vect2    ! allocatable vector for data , used by var_type1d and var_type2d to save avg values
          real(kind=WP),allocatable,dimension(:,:)     :: v ! name
          real(kind=WP),allocatable,dimension(:,:)     :: sn ! name
          integer                                      :: sz(2)
      end type var_vect2
      type  ::   var_type0d    ! type for temporal values for output
         real(kind=WP)                                :: snapshot ! last value
         type(var_vect0),allocatable,dimension(:)     :: avg ! name
         integer                          :: varid  ! ID of variable
         integer                          :: nfid    ! number of files used avg  ; avg(fid)
      end type var_type0d
      type  ::   var_type1d    ! type for temporal values for output (for example surface, eta_n)
         real(kind=WP),allocatable,dimension(:)       :: snapshot ! name
         type(var_vect1),allocatable,dimension(:)     :: avg ! name
         integer                          :: varid  ! ID of variable
         integer                          :: nfid    ! number of files used avg  ; avg(:,fid)
      end type var_type1d
      type  ::   var_type2d    ! type for temporal values for output (for example 3d temperature, TF)
         real(kind=WP),allocatable,dimension(:,:)     :: snapshot ! name
         type(var_vect2),allocatable,dimension(:)   :: avg ! name
         integer                          :: varid  ! ID of variable
         integer                          :: nfid    ! number of files used avg  ; avg(:,:,fid)
      end type var_type2d

      type  ::   varlist_type    ! information about variables for output
         character(len = 100) :: name ! name
         character(len = 100) :: modelname ! name of variable in model
         character(len = 256) :: standard_name  ! full variable name in the NetCDF file (attribute)
         character(len = 256) :: units  ! variable units in the NetCDF file
         integer              :: onnode ! if on nod =1 , if on elements =0, other =-1
         integer              :: timedim! time dimension 0 (not in time), 1 (in time)
         integer              :: vertdim! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
         integer              :: varid  ! ID of variable
         integer              :: type_task ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
         integer              :: ndim ! number of dimensions
         integer              :: varsize(2)  ! size of dimensions (1-if time then time, else nod/elem, 2 - if time then nod/elem, 3 - real size (or
      end type varlist_type

      type  ::   dimid_type    ! list of dimensions ID in netCDF files
         integer              :: node, nelem,maxnod, time, sigma, sigmam1, maxnodelem
      end type dimid_type
      type  ::   varid_type    ! list of variables ID in netCDF files
         integer  :: nv, mesh, lon, lat, lon_elem, lat_elem, time, &
                                 sigma, sigmam1, elem_area, w_cv, area, &
                                 nod_in_elem2D, nod_in_elem2D_num, depth, depth_elem,&
                                 node, elem
      end type varid_type
      integer,parameter      ::   maxnumbervars = 30 ! Maximum number of variables in one file
      integer,parameter      ::   maxnumberallvars = 100 ! Maximum number of declared variables


      type  ::   filelist_type    ! type for descripbing of output files
         character(len = 256) :: name     ! file name
         integer              :: snapshot ! if no averaging/snapshots =1 ,  averaging ocer freq(sec) =0
         integer              :: freq     ! frequency of output in seconds or steps (check freq units) , for now only steps
!         integer              :: frequnits! units of frequency: seconds =0, model steps =1
         integer              :: specfreq ! special frequency: read freq =0, daily=1, monthly=2, seasons=3, yearly=4
         integer              :: fileid   ! file ID
         integer              :: numvars  ! number of variables for output in this file
         integer,dimension(maxnumbervars):: vars     ! list of variables id from varlist (varid) which should be put in this file
         integer,dimension(maxnumbervars):: varsid   ! ID of each variable in netCDF file
         integer,dimension(maxnumbervars):: fid      ! index in var?d array
         integer,dimension(maxnumbervars):: avgid    ! index in var?d%avg array
         type(dimid_type)     :: dimid    ! list of dimensions ID
         type(varid_type)     :: varid    ! list of variables ID after dimension
         integer              :: nc_rec   ! record in NC file (=0 at inizialization)
         integer              :: ftype    ! type of file: 0 - ggeneral (full size outpu), 1 - 2d zoom to certan region, 2 - transaction between two points, 3 - station
         real(kind=WP)        :: posx_ll  ! position lower left X (longitude or x in cartesian) of the file if it is station/transection/field
         real(kind=WP)        :: posy_ll  ! position lower left Y (longitude or x in cartesian) of the file if it is station/transection/field
         real(kind=WP)        :: posx_ur  ! position upper right X (longitude or x in cartesian) of the file if it is station/transection/field
         real(kind=WP)        :: posy_ur  ! position upper right Y (longitude or x in cartesian) of the file if it is station/transection/field
         integer,allocatable,dimension(:) :: nods ! list of nodes indexes in global distribution  (one for station,...)
         integer,allocatable,dimension(:) :: elems ! list of nodes indexes in global distribution (one for station,...)
         integer,allocatable,dimension(:) :: mynods ! list of nodes indexes in local PE distribution  (one for station,...)
         integer,allocatable,dimension(:) :: myelems ! list of nodes indexes in local PE distribution (one for station,...)
         integer(kind=MPI_OFFSET_KIND) :: nnods ! number of nodes of my proc to write (eq 1 for station,...)
         integer(kind=MPI_OFFSET_KIND) :: nelems ! number of elements of my proc to write (eq 1 for station,...)
         integer(kind=MPI_OFFSET_KIND) :: nfilenods ! total number of nodes in file (one for station,...)
         integer(kind=MPI_OFFSET_KIND) :: nfileelems ! total number of elements in file (one for station,...)
         integer :: mynodstart ! starting position where in file to write nodes (=1 for station)
         integer :: myelemstart ! starting position where in file to write nodes (=1 for station)
         integer :: mystation ! 1 if station at my PE, else 0 (file will be used only by one PE)

       end type filelist_type

       type(filelist_type), allocatable, dimension(:), save   :: filelist
       type(varlist_type),  dimension(maxnumberallvars), save  :: varlist
       type(var_type0d),  allocatable, dimension(:), save  :: varbas0d
       type(var_type1d),  allocatable, dimension(:), save  :: varbas1d
       type(var_type2d),  allocatable, dimension(:), save  :: varbas2d
       integer, save                             :: nvar0d, nvar1d, nvar2d ! number of variables in varbas?d
       integer, save                             :: numfiles ! number of files for output

CONTAINS
   SUBROUTINE define_varlist
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE define_varlist ***
      !!
      !! ** Purpose : define variables in varlist
      !! ** Method  :
      !! ** Action  :
      !!----------------------------------------------------------------------
      IMPLICIT NONE
      integer              :: nvar  ! temporal variable ID
      integer              :: i, ndim
      integer              :: varsize(2)  ! size of dimensions
      character(len = 100) :: tmp,vname ! name
      character(len = 3)   :: vnamef ! name


      nvar=0
      nvar=nvar+1
      varlist(nvar)%name =           'eta' ! name
      varlist(nvar)%modelname =      'eta_n+depth(depth<0)' ! name of variable in model
      varlist(nvar)%standard_name =  'sea surface elevation'  ! full variable name in the NetCDF file
      varlist(nvar)%units =          'meter'  ! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        0! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      1 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =         nvar !1 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'u' ! name
      varlist(nvar)%modelname =      'U_n' ! name of variable in model
      varlist(nvar)%standard_name =  'velocity in x direction'! full variable name in the NetCDF file
      varlist(nvar)%units =          'm/s'! variable units in the NetCDF file
      varlist(nvar)%onnode =         0 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        2 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      2 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !2 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'v' ! name
      varlist(nvar)%modelname =      'V_n' ! name of variable in model
      varlist(nvar)%standard_name =  'velocity in y direction'! full variable name in the NetCDF file
      varlist(nvar)%units =          'm/s'! variable units in the NetCDF file
      varlist(nvar)%onnode =         0 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        2 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      2 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !3 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'ssu' ! name
      varlist(nvar)%modelname =      'U_n(1,:)' ! name of variable in model
      varlist(nvar)%standard_name =  'surface velocity in x direction'! full variable name in the NetCDF file
      varlist(nvar)%units =          'm/s'! variable units in the NetCDF file
      varlist(nvar)%onnode =         0 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      2 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !4 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'ssv' ! name
      varlist(nvar)%modelname =      'V_n(1,:)' ! name of variable in model
      varlist(nvar)%standard_name =  'surface velocity in y direction'! full variable name in the NetCDF file
      varlist(nvar)%units =          'm/s'! variable units in the NetCDF file
      varlist(nvar)%onnode =         0 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      2 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !5 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'u2d' ! name
      varlist(nvar)%modelname =      'U_n_2D(1,:)' ! name of variable in model
      varlist(nvar)%standard_name =  '2D velocity in x direction'! full variable name in the NetCDF file
      varlist(nvar)%units =          'm/s'! variable units in the NetCDF file
      varlist(nvar)%onnode =         0 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      1 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !6 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'v2d' ! name
      varlist(nvar)%modelname =      'U_n_2D(2,:)' ! name of variable in model
      varlist(nvar)%standard_name =  '2D velocity in y direction'! full variable name in the NetCDF file
      varlist(nvar)%units =          'm/s'! variable units in the NetCDF file
      varlist(nvar)%onnode =         0 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      1 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !7 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'temperature' ! name
      varlist(nvar)%modelname =      'TF' ! name of variable in model
      varlist(nvar)%standard_name =  'temperature'! full variable name in the NetCDF file
      varlist(nvar)%units =          'grad C'! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        2 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !8 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'salinity' ! name
      varlist(nvar)%modelname =      'SF' ! name of variable in model
      varlist(nvar)%standard_name =  'salinity'! full variable name in the NetCDF file
      varlist(nvar)%units =          'psu'! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        2 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !9 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'sst' ! name
      varlist(nvar)%modelname =      'TF(1,:)' ! name of variable in model
      varlist(nvar)%standard_name =  'temperature'! full variable name in the NetCDF file
      varlist(nvar)%units =          'grad C'! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !10 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'sss' ! name
      varlist(nvar)%modelname =      'SF(1,:)' ! name of variable in model
      varlist(nvar)%standard_name =  'salinity'! full variable name in the NetCDF file
      varlist(nvar)%units =          'psu'! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !11 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'taux' ! name
      varlist(nvar)%modelname =      'taux' ! name of variable in model
      varlist(nvar)%standard_name =  'wind stress along x direction'! full variable name in the NetCDF file
      varlist(nvar)%units =          'N/m2'! variable units in the NetCDF file
      varlist(nvar)%onnode =         0 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      1 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !12 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'tauy' ! name
      varlist(nvar)%modelname =      'tauy' ! name of variable in model
      varlist(nvar)%standard_name =  'wind stress along y direction'! full variable name in the NetCDF file
      varlist(nvar)%units =          'N/m2'! variable units in the NetCDF file
      varlist(nvar)%onnode =         0 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      1 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !13 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'taux_node' ! name
      varlist(nvar)%modelname =      'taux_node' ! name of variable in model
      varlist(nvar)%standard_name =  'wind stress along x direction on node'! full variable name in the NetCDF file
      varlist(nvar)%units =          'N/m2'! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      1 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !14 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'tauy_node' ! name
      varlist(nvar)%modelname =      'tauy_node' ! name of variable in model
      varlist(nvar)%standard_name =  'wind stress along y direction on node'! full variable name in the NetCDF file
      varlist(nvar)%units =          'N/m2'! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      1 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !15 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'mslp' ! name
      varlist(nvar)%modelname =      'mslp' ! name of variable in model
      varlist(nvar)%standard_name =  'mean sea level pressure'! full variable name in the NetCDF file
      varlist(nvar)%units =          'Pascals'! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      1 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !16 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'emp' ! name
      varlist(nvar)%modelname =      'emp' ! name of variable in model
      varlist(nvar)%standard_name =  'evaporation minus precipitation'! full variable name in the NetCDF file
      varlist(nvar)%units =          'kg/m2/s'! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      1 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !17 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'qns' ! name
      varlist(nvar)%modelname =      'qns' ! name of variable in model
      varlist(nvar)%standard_name =  'downward non solar heat over the ocean'! full variable name in the NetCDF file
      varlist(nvar)%units =          'W/m2'! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      1 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !18 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'qsr' ! name
      varlist(nvar)%modelname =      'qsr' ! name of variable in model
      varlist(nvar)%standard_name =  'downward solar heat over the ocean'! full variable name in the NetCDF file
      varlist(nvar)%units =          'W/m2'! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      1 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !19 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'zbar' ! name
      varlist(nvar)%modelname =      'zbar' ! name of variable in model
      varlist(nvar)%standard_name =  'node depth'! full variable name in the NetCDF file
      varlist(nvar)%units =          'm'! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        1 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      2 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !20 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'Z' ! name
      varlist(nvar)%modelname =      'Z' ! name of variable in model
      varlist(nvar)%standard_name =  'cell depth'! full variable name in the NetCDF file
      varlist(nvar)%units =          'm'! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        2 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      2 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !21 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'mask_wd_node' ! name
      varlist(nvar)%modelname =      'mask_wd_node' ! name of variable in model
      varlist(nvar)%standard_name =  'mask (wetting_drying) nodes'! full variable name in the NetCDF file
      varlist(nvar)%units =          ''! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      1 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !22 ! ID of variable
       nvar=nvar+1
      varlist(nvar)%name =           'mask_wd' ! name
      varlist(nvar)%modelname =      'mask_wd' ! name of variable in model
      varlist(nvar)%standard_name =  'mask (wetting_drying) elements'! full variable name in the NetCDF file
      varlist(nvar)%units =          ''! variable units in the NetCDF file
      varlist(nvar)%onnode =         0 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      1 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !23 ! ID of variable
       nvar=nvar+1
      varlist(nvar)%name =           'Cd' ! name
      varlist(nvar)%modelname =      'C_d_el' ! name of variable in model
      varlist(nvar)%standard_name =  'bottom friction elements'! full variable name in the NetCDF file
      varlist(nvar)%units =          ''! variable units in the NetCDF file
      varlist(nvar)%onnode =         0 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      1 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !24 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'windx' ! name
      varlist(nvar)%modelname =      'winds' ! name of variable in model
      varlist(nvar)%standard_name =  'wind along x direction on node'!
      varlist(nvar)%units =          'm/s'! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly
      varlist(nvar)%type_task =      1 ! T,S,... =3, U_n,V_n,... =2, eta_n,...
      varlist(nvar)%varid =          nvar !25 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'windy' ! name
      varlist(nvar)%modelname =      'windy' ! name of variable in model
      varlist(nvar)%standard_name =  'wind along y direction on node'!
      varlist(nvar)%units =          'm/s'! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly
      varlist(nvar)%type_task =      1 ! T,S,... =3, U_n,V_n,... =2, eta_n,...
      varlist(nvar)%varid =          nvar !26 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'w' ! name
      varlist(nvar)%modelname =      'Wvel' ! name of variable in model
      varlist(nvar)%standard_name =  'vertical velocity'! full variable name in the NetCDF file
      varlist(nvar)%units =          'm/s'! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        1 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      2 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !27 ! ID of variable

      nvar=nvar+1
      varlist(nvar)%name =           'snu' ! name
      varlist(nvar)%modelname =      'snu' ! name of variable in model
      varlist(nvar)%standard_name =  'coefficient of turbulent mixing'! full variable name in the NetCDF file
      varlist(nvar)%units =          ''! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        1 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !28 ! ID of variable

      nvar=nvar+1
      varlist(nvar)%name =           'bt' ! name
      varlist(nvar)%modelname =      'bt' ! name of variable in model
      varlist(nvar)%standard_name =  'turbulent kinetic energy'! full variable name in the NetCDF file
      varlist(nvar)%units =          ''! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        1 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !29 ! ID of variable

      nvar=nvar+1
      varlist(nvar)%name =           'ri' ! name
      varlist(nvar)%modelname =      'Ri' ! name of variable in model
      varlist(nvar)%standard_name =  'Richardson number'! full variable name in the NetCDF file
      varlist(nvar)%units =          ''! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        1 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !30 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'm_ice' ! name
      varlist(nvar)%modelname =      'm_ice' ! name of variable in model
      varlist(nvar)%standard_name =  'effective ice thickness'
      varlist(nvar)%units =          'm'! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,...
      varlist(nvar)%varid =          nvar !31 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'm_snow' ! name
      varlist(nvar)%modelname =      'm_snow' ! name of variable in model
      varlist(nvar)%standard_name =  'effective snow thickness'!
      varlist(nvar)%units =          'm'! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,...
      varlist(nvar)%varid =          nvar !32 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'a_ice' ! name
      varlist(nvar)%modelname =      'a_ice' ! name of variable in model
      varlist(nvar)%standard_name =  'ice concentation'!
      varlist(nvar)%units =          'percent'! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,...
      varlist(nvar)%varid =          nvar !33 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           't_skin' ! name
      varlist(nvar)%modelname =      't_skin' ! name of variable in model
      varlist(nvar)%standard_name =  'temperature of snow/ice top surface'!
      varlist(nvar)%units =          'grad C'! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,...
      varlist(nvar)%varid =          nvar !34 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'runoff_rivers' ! name
      varlist(nvar)%modelname =      'runoff_rivers' ! name of variable in model
      varlist(nvar)%standard_name =  'runoff rivers'!
      varlist(nvar)%units =          'm3/s'! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly
      varlist(nvar)%type_task =      1 ! T,S,... =3, U_n,V_n,... =2, eta_n,...
      varlist(nvar)%varid =          nvar !35 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'uice' ! name
      varlist(nvar)%modelname =      'U_n_ice' ! name of variable in model
      varlist(nvar)%standard_name =  'ice velocity in x direction'! full variable name in the NetCDF file
      varlist(nvar)%units =          'm/s'! variable units in the NetCDF file
      varlist(nvar)%onnode =         0 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !36 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'vice' ! name
      varlist(nvar)%modelname =      'U_n_ice' ! name of variable in model
      varlist(nvar)%standard_name =  'ice velocity in y direction'! full variable name in the NetCDF file
      varlist(nvar)%units =          'm/s'! variable units in the NetCDF file
      varlist(nvar)%onnode =         0 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !37 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'real_salt_flux' ! name
      varlist(nvar)%modelname =      'real_salt_flux' ! name of variable in model
      varlist(nvar)%standard_name =  'salinity flux from ice'!
      varlist(nvar)%units =          ''! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,...
      varlist(nvar)%varid =          nvar !38 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'tmp1n' ! name
      varlist(nvar)%modelname =      'anyvar' ! name of variable in model
      varlist(nvar)%standard_name =  'temporary output 1, nodes,2D '!
      varlist(nvar)%units =          ''! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,...
      varlist(nvar)%varid =          nvar !39 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'tmp2n' ! name
      varlist(nvar)%modelname =      'anyvar' ! name of variable in model
      varlist(nvar)%standard_name =  'temporary output 2, nodes,2D '!
      varlist(nvar)%units =          ''! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,...
      varlist(nvar)%varid =          nvar !40 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'tmp1e' ! name
      varlist(nvar)%modelname =      'anyvar' ! name of variable in model
      varlist(nvar)%standard_name =  'temporary output 1, elements,2D '!
      varlist(nvar)%units =          ''! variable units in the NetCDF file
      varlist(nvar)%onnode =         0 ! if on nod =1 , if on elements =0, other
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,...
      varlist(nvar)%varid =          nvar !41 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'tmp2e' ! name
      varlist(nvar)%modelname =      'anyvar' ! name of variable in model
      varlist(nvar)%standard_name =  'temporary output 2, elements,2D '!
      varlist(nvar)%units =          ''! variable units in the NetCDF file
      varlist(nvar)%onnode =         0 ! if on nod =1 , if on elements =0, other
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,...
      varlist(nvar)%varid =          nvar !42 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'xwind' ! name
      varlist(nvar)%modelname =      'atmdata' ! name of variable in model
      varlist(nvar)%standard_name =  'wind X direction'!
      varlist(nvar)%units =          ''! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,...
      varlist(nvar)%varid =          nvar !43 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'ywind' ! name
      varlist(nvar)%modelname =      'atmdata' ! name of variable in model
      varlist(nvar)%standard_name =  'wind Y direction'!
      varlist(nvar)%units =          ''! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,...
      varlist(nvar)%varid =          nvar !44 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'humi' ! name
      varlist(nvar)%modelname =      'atmdata' ! name of variable in model
      varlist(nvar)%standard_name =  'humidity'!
      varlist(nvar)%units =          ''! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,...
      varlist(nvar)%varid =          nvar !45 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'qsr_raw' ! name
      varlist(nvar)%modelname =      'atmdata' ! name of variable in model
      varlist(nvar)%standard_name =  'short wave radiation'!
      varlist(nvar)%units =          ''! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,...
      varlist(nvar)%varid =          nvar !46 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'qlw_raw' ! name
      varlist(nvar)%modelname =      'atmdata' ! name of variable in model
      varlist(nvar)%standard_name =  'long wave radiation'!
      varlist(nvar)%units =          ''! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,...
      varlist(nvar)%varid =          nvar !47 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'tair' ! name
      varlist(nvar)%modelname =      'atmdata' ! name of variable in model
      varlist(nvar)%standard_name =  '2m air temperature'!
      varlist(nvar)%units =          ''! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,...
      varlist(nvar)%varid =          nvar !48 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'prec' ! name
      varlist(nvar)%modelname =      'atmdata' ! name of variable in model
      varlist(nvar)%standard_name =  'total precipitation'!
      varlist(nvar)%units =          ''! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,...
      varlist(nvar)%varid =          nvar !49 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'mslp' ! name
      varlist(nvar)%modelname =      'atmdata' ! name of variable in model
      varlist(nvar)%standard_name =  'mean sea level pressure'!
      varlist(nvar)%units =          ''! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,...
      varlist(nvar)%varid =          nvar !50 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'ice_delta' ! name
      varlist(nvar)%modelname =      'ice_delta' ! name of variable in model
      varlist(nvar)%standard_name =  'ice delta'!
      varlist(nvar)%units =          ''! variable units in the NetCDF file
      varlist(nvar)%onnode =         0 ! if on nod =1 , if on elements =0, other
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,...
      varlist(nvar)%varid =          nvar !51 ! ID of variable
      nvar=nvar+1
      varlist(nvar)%name =           'ice_div' ! name
      varlist(nvar)%modelname =      'ice_div' ! name of variable in model
      varlist(nvar)%standard_name =  'ice divergence'!
      varlist(nvar)%units =          ''! variable units in the NetCDF file
      varlist(nvar)%onnode =         0 ! if on nod =1 , if on elements =0, other
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,...
      varlist(nvar)%varid =          nvar !52 ! ID of variable   
      nvar=nvar+1
      varlist(nvar)%name =           'ice_share' ! name
      varlist(nvar)%modelname =      'ice_share' ! name of variable in model
      varlist(nvar)%standard_name =  'ice share'!
      varlist(nvar)%units =          ''! variable units in the NetCDF file
      varlist(nvar)%onnode =         0 ! if on nod =1 , if on elements =0, other
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in
      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,...
      varlist(nvar)%varid =          nvar !53 ! ID of variable            

      nvar=nvar+1
      varlist(nvar)%name =           'tke_dissip' ! name
      varlist(nvar)%modelname =      'tke_dissip' ! name of variable in model
      varlist(nvar)%standard_name =  'turbulent kinetic energy dissipation'! full variable name in the NetCDF file
      varlist(nvar)%units =          ''! variable units in the NetCDF file
      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other =-1
      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
      varlist(nvar)%vertdim =        1 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
      varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
      varlist(nvar)%varid =          nvar !54 ! ID of variable
!      if (key_passive_tracers) then
!        do i = 1, i_pastr
!          nvar=nvar+1
!          !write(tmp,*) i
!          !write(vname,*) 'tr',trim(ADJUSTL(trim(tmp)))
!          write(vnamef,"(A2,I1)") "tr",i
!          write(*,*) trim(vnamef)
!          varlist(nvar)%name          =  trim(vnamef)
!          varlist(nvar)%modelname =      'Tpass' ! name of variable in model
!          varlist(nvar)%standard_name =  'Tpass' ! full variable name in the NetCDF file
!          varlist(nvar)%units =          ''! variable units in the NetCDF file
!          varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other =-1
!          varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
!          varlist(nvar)%vertdim =        2 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
!          varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
!          varlist(nvar)%varid =          nvar !25..size(biomodel%state_variables) ! ID of variable
!        enddo
!        write(*,*) "total passive tracers: ", nvar-24
!       endif
!      if (key_fabm) then
!        do i = 1, size(biomodel%state_variables)
!          nvar=nvar+1
!          varlist(nvar)%name = biomodel%state_variables(i)%name        ! name
!          varlist(nvar)%modelname =      'Tpass(:,:,1)' ! name of variable in model
!          varlist(nvar)%standard_name =  biomodel%state_variables(i)%long_name ! full variable name in the NetCDF file
!          varlist(nvar)%units =          biomodel%state_variables(i)%units! variable units in the NetCDF file
!          varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other =-1
!          varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
!          varlist(nvar)%vertdim =        2 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
!          varlist(nvar)%type_task =      3 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
!          varlist(nvar)%varid =          nvar !25..size(biomodel%state_variables) ! ID of variable
!        enddo
!        write(*,*) "total state_variables: ", nvar-24
!       endif
!      nvar=nvar+1
!      varlist(nvar)%name =           'emp_2' ! name
!      varlist(nvar)%modelname =      'emp_2' ! name of variable in model
!      varlist(nvar)%standard_name =  'evaporation minus precipitation 2'! full variable name in the NetCDF file
!      varlist(nvar)%units =          'kg/m2/s'! variable units in the NetCDF file
!      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other =-1
!      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
!      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
!      varlist(nvar)%type_task =      1 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
!      varlist(nvar)%varid =          nvar !25 ! ID of variable
!      nvar=nvar+1
!      varlist(nvar)%name =           'qns_2' ! name
!      varlist(nvar)%modelname =      'qns' ! name of variable in model
!      varlist(nvar)%standard_name =  'downward non solar heat over the ocean'! full variable name in the NetCDF file
!      varlist(nvar)%units =          'W/m2'! variable units in the NetCDF file
!      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other =-1
!      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
!      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
!      varlist(nvar)%type_task =      1 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
!      varlist(nvar)%varid =          nvar !26 ! ID of variable
!      nvar=nvar+1
!      varlist(nvar)%name =           'qsr_2' ! name
!      varlist(nvar)%modelname =      'qsr_2' ! name of variable in model
!      varlist(nvar)%standard_name =  'downward solar heat over the ocean'! full variable name in the NetCDF file
!      varlist(nvar)%units =          'W/m2'! variable units in the NetCDF file
!      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other =-1
!      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
!      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
!      varlist(nvar)%type_task =      1 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
!      varlist(nvar)%varid =          nvar !27 ! ID of variable
!      nvar=nvar+1
!      varlist(nvar)%name =           'taux_node_2' ! name
!      varlist(nvar)%modelname =      'taux_node_2' ! name of variable in model
!      varlist(nvar)%standard_name =  'wind stress along x direction on node'! full variable name in the NetCDF file
!      varlist(nvar)%units =          'N/m2'! variable units in the NetCDF file
!      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other =-1
!      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
!      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
!      varlist(nvar)%type_task =      1 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
!      varlist(nvar)%varid =          nvar !28 ! ID of variable
!      nvar=nvar+1
!      varlist(nvar)%name =           'tauy_node_2' ! name
!      varlist(nvar)%modelname =      'tauy_node_2' ! name of variable in model
!      varlist(nvar)%standard_name =  'wind stress along y direction on node'! full variable name in the NetCDF file
!      varlist(nvar)%units =          'N/m2'! variable units in the NetCDF file
!      varlist(nvar)%onnode =         1 ! if on nod =1 , if on elements =0, other =-1
!      varlist(nvar)%timedim =        1 ! time dimension 0 (not in time), 1 (in time)
!      varlist(nvar)%vertdim =        0 ! vertical dimension 0 , 1 (verticaly resolved) sigma, 2 sigmam1
!      varlist(nvar)%type_task =      1 ! T,S,... =3, U_n,V_n,... =2, eta_n,... =1
!      varlist(nvar)%varid =          nvar !29 ! ID of variable


      if (nvar > maxnumberallvars) then
         write(*,*) "ERROR: too much declared variables; change maxnumberallvars in fv_ncoutput"
         STOP "ERROR: check maxnumberallvars"
         !well, it does not have a sense, if so, memory error will be.
      end if

      do i = 1, nvar
      ! count number of dimensions for each variable
         ndim = 0
         select case (varlist(i)%vertdim)
            case (0)
            case (1)
               ndim = ndim+1
            case (2)
               ndim = ndim+1
         end select
         select case (varlist(i)%onnode)
            case (-1)
            case (0)
               ndim = ndim+1
            case (1)
               ndim = ndim+1
         end select
         select case (varlist(i)%timedim)
            case (0)
            case (1)
               ndim = ndim+1
         end select
         varlist(i)%ndim = ndim

         ! get size of variable
         call getvarsize(i,varsize)
         varlist(i)%varsize(1) = varsize(1)
         varlist(i)%varsize(2) = varsize(2)
      end do
   END SUBROUTINE define_varlist
!======================   Put your 0D variables here: ========================================
   SUBROUTINE getvar0d(i)
      IMPLICIT NONE
      integer, intent(in) :: i !

      select case (varbas0d(i)%varid)
!         case(?)
!            varbas0d(i)%snapshot = ???
      end select
   END SUBROUTINE getvar0d
!======================   Put your 1D variables here: ========================================
   SUBROUTINE getvar1d(i)
#ifdef USE_MPI
      USE fv_sbcmpi      ! MPI version of surface boundary cond.
#else
      USE fv_sbc
#endif

      IMPLICIT NONE
      integer, intent(in) :: i !
      integer             :: n !

      select case (varbas1d(i)%varid)
         case(1)
            varbas1d(i)%snapshot = eta_n
            ! if we have wettingdrying and depth<0 then add depth to eta
            do n=1,myDim_nod2D
               if (depth(n) < 0.0_WP) then
                  varbas1d(i)%snapshot(n) = eta_n(n) + depth(n)
               endif
            enddo
         case(4)
            varbas1d(i)%snapshot =  U_n(1,:)
         case(5)
            varbas1d(i)%snapshot =  V_n(1,:)
         case(6)
            varbas1d(i)%snapshot =  U_n_2D(1,:)
         case(7)
            varbas1d(i)%snapshot =  U_n_2D(2,:)
         case(10)
            varbas1d(i)%snapshot =  TF(1,:)
         case(11)
            varbas1d(i)%snapshot =  SF(1,:)
         case(12)
            varbas1d(i)%snapshot =  taux
         case(13)
            varbas1d(i)%snapshot =  tauy
         case(14)
            varbas1d(i)%snapshot =  taux_node
         case(15)
            varbas1d(i)%snapshot =  tauy_node
         case(16)
            varbas1d(i)%snapshot =  mslp
         case(17)
            varbas1d(i)%snapshot =  emp
         case(18)
            varbas1d(i)%snapshot =  qns
         case(19)
            varbas1d(i)%snapshot =  qsr
         case(22)
            varbas1d(i)%snapshot =  mask_wd_node
         case(23)
            varbas1d(i)%snapshot =  mask_wd
         case(24)
            varbas1d(i)%snapshot =  C_d_el
         case(25)
            varbas1d(i)%snapshot =  windx
         case(26)
            varbas1d(i)%snapshot =  windy
         case(31)
            if (use_ice) then
                varbas1d(i)%snapshot =  m_ice
            else
                varbas1d(i)%snapshot =  0.0_WP
            endif
         case(32)
             if (use_ice) then
                varbas1d(i)%snapshot =  m_snow
            else
                varbas1d(i)%snapshot =  0.0_WP
            endif
         case(33)
            if (use_ice) then
                varbas1d(i)%snapshot =  a_ice
            else
                varbas1d(i)%snapshot =  0.0_WP
            endif
         case(34)
            if (use_ice) then
                varbas1d(i)%snapshot =  t_skin
            else
                varbas1d(i)%snapshot =  0.0_WP
            endif
         case(35)
            if (key_rivers) then
                varbas1d(i)%snapshot =  0.0_WP !runoff_rivers - is allocated only on (my_nod_rivers)
            else
                varbas1d(i)%snapshot =  0.0_WP
            endif
         case(36)
            if (use_ice) then
                varbas1d(i)%snapshot =  U_n_ice(1,:)
            else
                varbas1d(i)%snapshot =  0.0_WP
            endif
         case(37)
            if (use_ice) then
                varbas1d(i)%snapshot =  U_n_ice(2,:)
            else
                varbas1d(i)%snapshot =  0.0_WP
            endif
         case(38)
            if (use_ice) then
                varbas1d(i)%snapshot =  real_salt_flux
            else
                varbas1d(i)%snapshot =  0.0_WP
            endif
         case(39)
            if (use_ice) then
                varbas1d(i)%snapshot =  cHrhsi(1,:)
            else
                varbas1d(i)%snapshot =  0.0_WP
            endif
         case(40)
            if (use_ice) then
                varbas1d(i)%snapshot =  tmparr(1,:)
            else
                varbas1d(i)%snapshot =  0.0_WP
            endif
         case(41)
            if (use_ice) then
                varbas1d(i)%snapshot =  tmparr_e(1,:)
            else
                varbas1d(i)%snapshot =  0.0_WP
            endif
         case(42)
            if (use_ice) then
                varbas1d(i)%snapshot =  0.0_WP
            else
                varbas1d(i)%snapshot =  0.0_WP
            endif
         case(43)
            if (key_atm) then
                varbas1d(i)%snapshot =  atmdata(1,:)
            else
                varbas1d(i)%snapshot =  0.0_WP
            endif
         case(44)
            if (key_atm) then
                varbas1d(i)%snapshot =  atmdata(2,:)
            else
                varbas1d(i)%snapshot =  0.0_WP
            endif
         case(45)
            if (key_atm) then
                varbas1d(i)%snapshot =  atmdata(3,:)
            else
                varbas1d(i)%snapshot =  0.0_WP
            endif
         case(46)
            if (key_atm) then
                varbas1d(i)%snapshot =  atmdata(4,:)
            else
                varbas1d(i)%snapshot =  0.0_WP
            endif
         case(47)
            if (key_atm) then
                varbas1d(i)%snapshot =  atmdata(5,:)
            else
                varbas1d(i)%snapshot =  0.0_WP
            endif
         case(48)
            if (key_atm) then
                varbas1d(i)%snapshot =  atmdata(6,:)
            else
                varbas1d(i)%snapshot =  0.0_WP
            endif
         case(49)
            if (key_atm) then
                varbas1d(i)%snapshot =  atmdata(7,:)
            else
                varbas1d(i)%snapshot =  0.0_WP
            endif
         case(50)
            if (key_atm) then
                varbas1d(i)%snapshot =  atmdata(8,:)
            else
                varbas1d(i)%snapshot =  0.0_WP
            endif
         case(51)
            if (use_ice) then
                varbas1d(i)%snapshot =  ice_delta(:)
            else
                varbas1d(i)%snapshot =  0.0_WP
            endif
         case(52)
            if (use_ice) then
                varbas1d(i)%snapshot =  ice_div(:)
            else
                varbas1d(i)%snapshot =  0.0_WP
            endif
         case(53)
            if (use_ice) then
                varbas1d(i)%snapshot =  ice_share(:)
            else
                varbas1d(i)%snapshot =  0.0_WP
            endif                        
      end select
      ! multiply by mask for wetting drying
      if (varlist(varbas1d(i)%varid)%onnode == 1) then
!         write(*,*) "i=",i," varid=",varbas1d(i)%varid," size=",size(varbas1d(i)%snapshot)," size2=",size(mask_wd_node)
         varbas1d(i)%snapshot = varbas1d(i)%snapshot * mask_wd_node
!         write(*,*) "i=",i," varid=",varbas1d(i)%varid," DONE"
      else
         varbas1d(i)%snapshot = varbas1d(i)%snapshot * mask_wd
      endif
   END SUBROUTINE getvar1d
!======================   Put your 2D variables here: ========================================
   SUBROUTINE getvar2d(i)
      USE o_ARRAYS

      IMPLICIT NONE
      integer, intent(in) :: i !
      integer  :: nz, nlev

      select case (varbas2d(i)%varid)
         case(2)
            varbas2d(i)%snapshot =  U_n
         case(3)
            varbas2d(i)%snapshot =  V_n
         case(8)
            varbas2d(i)%snapshot =  TF
         case(9)
            varbas2d(i)%snapshot =  SF
         case(20)
            varbas2d(i)%snapshot =  zbar
         case(21)
            varbas2d(i)%snapshot =  z
         case(27)
            varbas2d(i)%snapshot =  Wvel
         case(28)
            varbas2d(i)%snapshot =  snu
         case(29)
            varbas2d(i)%snapshot =  bt
         case(54)
            varbas2d(i)%snapshot =  tke_dissip

!         case(25:24+9)!size(biomodel%state_variables))
!            varbas2d(i)%snapshot =  Tpass(:,:,varbas2d(i)%varid-24)
      end select
      ! multiply by mask for wetting drying
      if (varlist(varbas2d(i)%varid)%vertdim == 1) then
         nlev = nsigma
      else
         nlev = nsigma - 1
      endif
      if (varlist(varbas2d(i)%varid)%onnode == 1) then
         do nz = 1, nlev
            varbas2d(i)%snapshot(nz,:) = varbas2d(i)%snapshot(nz,:) * mask_wd_node
         enddo
      else
         do nz = 1, nlev
            varbas2d(i)%snapshot(nz,:) = varbas2d(i)%snapshot(nz,:) * mask_wd
         enddo
      endif
   END SUBROUTINE getvar2d
!=============================================================================================
   SUBROUTINE ncoutputmpi_do
      IMPLICIT NONE

      integer :: i

!      call getvar
!      call addvar2basket !add variables to varbas?d lists

      do i = 1, numfiles
         if (filelist(i)%ftype==3 .and. filelist(i)%mystation==0) then
            cycle
         end if
         if ( mod(n_dt,filelist(i)%freq)==0 ) then
            ! if it is time to make output
            call addvar2basket !add variables to varbas?d lists
               ! if we need to calculate average values
!            if (filelist(i)%snapshot == 2) call calculateavg(i)

            call write_data(i) ! write values of snapshots/avg to nc file
             ! put average to 0
!            if (filelist(i)%snapshot == 2) call zeroavg(i)

         end if
      end do

   END SUBROUTINE ncoutputmpi_do

   SUBROUTINE zeroavg(findx)
      IMPLICIT NONE
      integer, intent(in)  :: findx  ! index of file to do averaging
      integer              :: i, fid, avgid
      ! variables one by one
      do i = 1, filelist(findx)%numvars
         fid = filelist(findx)%fid(i)
         avgid = filelist(findx)%avgid(i)
         select case (varlist(filelist(findx)%vars(i))%ndim) ! select proper dimension
            case(1)
               varbas0d(fid)%avg(avgid)%v = 0.0_WP
            case(2)
               varbas1d(fid)%avg(avgid)%v = 0.0_WP
            case(3)
               varbas2d(fid)%avg(avgid)%v = 0.0_WP
            end select
      end do

   END SUBROUTINE zeroavg
   SUBROUTINE calculateavg(findx)
      IMPLICIT NONE
      integer, intent(in)  :: findx  ! index of file to do averaging
      integer              :: i, fid, avgid


      ! variables one by one
      do i = 1, filelist(findx)%numvars
         fid = filelist(findx)%fid(i)
         avgid = filelist(findx)%avgid(i)
         select case (varlist(filelist(findx)%vars(i))%ndim) ! select proper dimension
            case(1)
               varbas0d(fid)%avg(avgid)%v = varbas0d(fid)%avg(avgid)%v / real(filelist(findx)%freq)
            case(2)
               varbas1d(fid)%avg(avgid)%v = varbas1d(fid)%avg(avgid)%v / real(filelist(findx)%freq)
            case(3)
               varbas2d(fid)%avg(avgid)%v = varbas2d(fid)%avg(avgid)%v / real(filelist(findx)%freq)
            end select
      end do

   END SUBROUTINE calculateavg

   SUBROUTINE ncoutputmpi_ini
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE ncoutputmpi_ini ***
      !!
      !! ** Purpose : inizialization of pnetcdf output
      !! ** Method  :
      !! ** Action  :
      !!----------------------------------------------------------------------
      USE pnetcdf
      IMPLICIT NONE


      real(wp), allocatable, dimension(:)  :: x, y, depth_elem
      real(wp), allocatable, dimension(:)  :: nod_x, nod_y
      real(wp)     :: tmpx(1),tmpy(1)
      integer      :: i,j, k, n, elnodes(4), output_alloc

      integer      :: filelist_unit
      integer      :: cmode, errn
      character(len=30)   :: str
      integer            :: yyyy,mm,dd, ihh,imm,isec
      character(4) :: yystr
      character(2) :: mmstr,ddstr,hr,mn,sc
      character(len=130)   :: fpath

      integer              :: ndim, varid,  loc  ! temporal variables
      integer,dimension(maxnumberallvars)  :: var0d,var1d,var2d,numfid0d,numfid1d,numfid2d !temporal lists of variables
      integer, allocatable ,dimension(:,:)  :: fid1d_size,fid2d_size !tmp to save size of avg%v

      integer,dimension(maxnumbervars):: varssh1    ! temporal array

      integer :: mynodstart,myelemstart,felem2D, tmpi
      integer, allocatable ,dimension(:,:)  :: elem2D_nodes_inglob
      integer(kind=MPI_OFFSET_KIND) :: max_nod_in_elem2D_num
      integer(kind=MPI_OFFSET_KIND) :: mpikind !temporal integer type of mpi offset
      var0d = 0
      var1d = 0
      var2d = 0
      numfid0d = 0
      numfid1d = 0
      numfid2d = 0
      if (mype==0) print *, 'inizialization of PnetCDF MPI output: start'
!print *,'nod2D,elem2D'
!print *, nod2D,elem2D,myDim_nod2D,myDim_elem2D
!STOP

!      print *, mype, remPtr_nod2D(mype+1)

      call define_varlist

      call get_fullstartpostitions(mynodstart,myelemstart,felem2D)

      ALLOCATE( x(myDim_elem2D), y(myDim_elem2D), depth_elem(myDim_elem2D), &
                   &      STAT=output_alloc )
      if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays x,y,depth_elem'
      allocate( nod_x(myDim_nod2D), nod_y(myDim_nod2D), &
                   &      STAT=output_alloc )
      if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays nod_x, nod_y'

      ! get coordinates of elements
      do i = 1, myDim_elem2D
         elnodes = elem2D_nodes(:,i) !! 4 nodes in element
         y(i) = w_cv(1,i)*coord_nod2D(2,elnodes(1)) &
                + w_cv(2,i)*coord_nod2D(2,elnodes(2)) &
                + w_cv(3,i)*coord_nod2D(2,elnodes(3)) &
                + w_cv(4,i)*coord_nod2D(2,elnodes(4))
         x(i) = w_cv(1,i)*coord_nod2D(1,elnodes(1)) &
                + w_cv(2,i)*coord_nod2D(1,elnodes(2)) &
                + w_cv(3,i)*coord_nod2D(1,elnodes(3)) &
                + w_cv(4,i)*coord_nod2D(1,elnodes(4))
         depth_elem(i) = w_cv(1,i)*depth(elem2D_nodes(1,i)) &
                       + w_cv(2,i)*depth(elem2D_nodes(2,i)) &
                       + w_cv(3,i)*depth(elem2D_nodes(3,i)) &
                       + w_cv(4,i)*depth(elem2D_nodes(4,i))
      end do
      if (cartesian) then
        y = y*r_earth
        x = x*r_earth
      else
        y = y/rad
        x = x/rad
      end if

      filelist_unit= 15
      ! define files for output
      ! Read structure from outputfilelist.dat
      open( unit=filelist_unit, FILE='outputfilelist.dat', status='old', iostat=iost )
      if( iost == 0 ) then
         if(mype==0) WRITE(*,*) '     file   : outputfilelist.dat open ok; mype=',mype
      else
         WRITE(*,*) 'ERROR: --> bad opening file   outputfilelist.dat ; iostat=',iost
         STOP 'ERROR: --> output_ini'
      endif
      ! reading number of files for output
      READ( unit=filelist_unit, FMT=*,iostat=iost)
      READ( unit=filelist_unit, FMT=*,iostat=iost) numfiles
      ALLOCATE( filelist(numfiles), &
                     &      STAT=output_alloc )
         if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays filelist'
      ! read files options from outputfilelist.dat
      if(mype==0) then
        write(*,*) "start: read/check variable and file lists for output"
        write(*,*) "files for output: "
      endif
      do i = 1, numfiles
         READ( unit=filelist_unit, FMT=*,iostat=iost)
         READ( unit=filelist_unit, FMT=*,iostat=iost) filelist(i)%name
         write(fpath,*) './nc_output/',trim(ADJUSTL(trim(filelist(i)%name)))
         filelist(i)%name = trim(ADJUSTL(trim(fpath)))
         READ( unit=filelist_unit, FMT=*,iostat=iost)
         READ( unit=filelist_unit, FMT=*,iostat=iost) filelist(i)%snapshot
         READ( unit=filelist_unit, FMT=*,iostat=iost)
         READ( unit=filelist_unit, FMT=*,iostat=iost) filelist(i)%ftype
         READ( unit=filelist_unit, FMT=*,iostat=iost)
         READ( unit=filelist_unit, FMT=*,iostat=iost) filelist(i)%posx_ll, &
                    & filelist(i)%posy_ll, filelist(i)%posx_ur, filelist(i)%posy_ur
         READ( unit=filelist_unit, FMT=*,iostat=iost)
         READ( unit=filelist_unit, FMT=*,iostat=iost) filelist(i)%freq
         READ( unit=filelist_unit, FMT=*,iostat=iost)
!         READ( unit=filelist_unit, FMT=*,iostat=iost) filelist(i)%frequnits
!         READ( unit=filelist_unit, FMT=*,iostat=iost)
         READ( unit=filelist_unit, FMT=*,iostat=iost) filelist(i)%specfreq
         READ( unit=filelist_unit, FMT=*,iostat=iost)
         READ( unit=filelist_unit, FMT=*,iostat=iost) filelist(i)%numvars
         if ( filelist(i)%numvars > maxnumbervars ) then
            write(*,*) "ERROR: number of variables",filelist(i)%numvars, "in output file:",trim(filelist(i)%name)
            write(*,*) "          more then maximum allowed (",maxnumbervars,"). hint: change maxnumbervars in fv_ncoutputmpi.f90"
            STOP
         end if
         READ( unit=filelist_unit, FMT=*,iostat=iost)
         filelist(i)%vars=0
         filelist(i)%varsid=0
         READ( unit=filelist_unit, FMT=*,iostat=iost) filelist(i)%vars

         do j =1, filelist(i)%numvars
            if ( filelist(i)%vars(j) < 1 ) then
               write(*,*) 'ERROR: inconsistent number of variables with provided index,',trim(filelist(i)%name)
               stop
            end if
         end do
         j=0
         do while ( j < filelist(i)%numvars )
            j=j+1
            if ( varlist(filelist(i)%vars(j))%type_task > type_task ) then
               if(mype==0) write(*,*) 'WARNING: variable "',trim(varlist(filelist(i)%vars(j))%name), &
                                             & '" removed from file:',trim(filelist(i)%name)
               if(mype==0) write(*,*) '            while type_task inconsistent with current model type_task'
               filelist(i)%numvars = filelist(i)%numvars - 1
               varssh1 = cshift(filelist(i)%vars,1)
               filelist(i)%vars(j:maxnumbervars-1) = varssh1(j:maxnumbervars-1)
               j=j-1
            end if
         end do
          if(mype==0) write(*,*) "          ",trim(filelist(i)%name)
      end do
      if(mype==0) write(*,*) "end: read/check variable and file lists for output"
   !==================== define indexes of nodes and elements where to write
      do i = 1, numfiles
         if ( filelist(i)%ftype==0 ) then
         ! full setup
            filelist(i)%nnods  = myDim_nod2D
            filelist(i)%nelems = myDim_elem2D
            filelist(i)%nfilenods  = nod2D
            filelist(i)%nfileelems = felem2D
            filelist(i)%mynodstart = mynodstart
            filelist(i)%myelemstart = myelemstart
            filelist(i)%mystation = 0

         elseif ( filelist(i)%ftype==3 ) then
         ! station
            !STOP 'output_ini: stations in mpi ;not supported output file type'
            filelist(i)%nnods  = 1
            filelist(i)%nelems = 1
            filelist(i)%nfilenods  = 1
            filelist(i)%nfileelems = 1
            filelist(i)%mynodstart = 1
            filelist(i)%myelemstart = 1
            filelist(i)%mystation = 0
         else
            STOP 'output_ini: wrong/not supported output file type'
         endif

         ALLOCATE(filelist(i)%nods(filelist(i)%nnods) , STAT=output_alloc )
         if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays filelist(i)%nods'
         ALLOCATE(filelist(i)%elems(filelist(i)%nelems) , STAT=output_alloc )
         if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays filelist(i)%elems'
         ALLOCATE(filelist(i)%mynods(filelist(i)%nnods) , STAT=output_alloc )
         if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays filelist(i)%mynods'
         ALLOCATE(filelist(i)%myelems(filelist(i)%nelems) , STAT=output_alloc )
         if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays filelist(i)%myelems'

         if ( filelist(i)%ftype==0 ) then
         ! full setup
            do n = 1, myDim_nod2D
               filelist(i)%nods(n) = myList_nod2D(n)
               filelist(i)%mynods(n) = n
            enddo
            do n = 1, myDim_elem2D
               filelist(i)%elems(n) = myList_elem2D(n)
               filelist(i)%myelems(n) = n
            enddo
         elseif ( filelist(i)%ftype==2 ) then
            !! Add here subroutine to find nodes,elements between points ?
            write(*,*) "write on the path is not implemented yet, stop"
            stop
         elseif ( filelist(i)%ftype==3 ) then
         ! station
            !STOP 'output_ini: stations in mpi ;not supported output file type'
            call ne_nodeelem(filelist(i)%posx_ll, filelist(i)%posy_ll, &
                     & filelist(i)%mynods(1),filelist(i)%myelems(1),filelist(i)%mystation)
            filelist(i)%nods(1) = myList_nod2D(filelist(i)%mynods(1))
            filelist(i)%elems(1) = myList_elem2D(filelist(i)%myelems(1))
            !call ne_element(filelist(i)%posx_ll, filelist(i)%posy_ll, filelist(i)%elems(1))
         endif
      enddo
    !=end================== define indexes of nodes and elements where to write


      !create files
      call MPI_AllREDUCE(maxval(nod_in_elem2D_num),tmpi, 1, MPI_INTEGER, MPI_MAX, &
          MPI_COMM_FESOM_C, MPIerr)
      max_nod_in_elem2D_num = tmpi
      do i = 1, numfiles
         !create new file
         cmode = IOR(NF90_CLOBBER, NF90_64BIT_DATA)
         if (filelist(i)%ftype==0) then
            ! if common file with data at all PE
            iost = nf90mpi_create(MPI_COMM_WORLD, filelist(i)%name, cmode, &
                         MPI_INFO_NULL, filelist(i)%fileid)
         elseif (filelist(i)%ftype==3 .and. filelist(i)%mystation==1) then
            ! if file is a station at current PE
            iost = nf90mpi_create(MPI_COMM_SELF, filelist(i)%name, cmode, &
                         MPI_INFO_NULL, filelist(i)%fileid)
         else
            cycle
         endif

         call check_nferr(iost, 'Error at nf90mpi_create ')
         ! define dimensions
         iost = nf90mpi_def_dim(filelist(i)%fileid, 'node', filelist(i)%nfilenods, filelist(i)%dimid%node)
         call check_nferr(iost, 'Error at nf90mpi_def_dim ')
         errn = nf90mpi_def_dim(filelist(i)%fileid, 'nele', filelist(i)%nfileelems, filelist(i)%dimid%nelem)
         call check_nferr(iost, 'Error at nf90mpi_def_dim ')
         iost = nf90mpi_def_dim(filelist(i)%fileid, 'maxnod', 4_MPI_OFFSET_KIND, filelist(i)%dimid%maxnod)
         call check_nferr(iost, 'Error at nf90mpi_def_dim ')
         iost = nf90mpi_def_dim(filelist(i)%fileid, 'maxnodelem', max_nod_in_elem2D_num, filelist(i)%dimid%maxnodelem)
         call check_nferr(iost, 'Error at nf90mpi_def_dim ')
         iost = nf90mpi_def_dim(filelist(i)%fileid, 'time', NF90MPI_UNLIMITED, filelist(i)%dimid%time)
         call check_nferr(iost, 'Error at nf90mpi_def_dim ')
         if (type_task>1) then
            mpikind = nsigma
            iost = nf90mpi_def_dim(filelist(i)%fileid, 'sigma', mpikind, filelist(i)%dimid%sigma)
            call check_nferr(iost, 'Error at nf90mpi_def_dim ')
            iost = nf90mpi_def_dim(filelist(i)%fileid, 'sigmam1', mpikind-1, filelist(i)%dimid%sigmam1)
            call check_nferr(iost, 'Error at nf90mpi_def_dim ')
         endif
!         !define sigma var
         if (type_task>1) then
            nf_dims(1)=filelist(i)%dimid%sigma
            iost = nfmpi_def_var(filelist(i)%fileid, 'sigma_lev', NF90_DOUBLE,1,nf_dims,filelist(i)%varid%sigma)
            call check_nferr(iost, 'Error at nf90mpi_def_var ')
         endif

!         !define nv
         nf_dims(1) = filelist(i)%dimid%maxnod
         nf_dims(2) = filelist(i)%dimid%nelem
         iost = nfmpi_def_var(filelist(i)%fileid,'nv', NF_INT, 2,nf_dims,filelist(i)%varid%nv)
         call check_nferr(iost, 'Error at nf90mpi_def_var ')
         ! add attributes to face_node_connectivity
         mpikind = len_trim('face_node_connectivity')
         iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%nv,'standard_name', &
                                  & mpikind,'face_node_connectivity')
         call check_nferr(iost, 'Error at nf90mpi_put_att_text ')

         iost = nfmpi_def_var(filelist(i)%fileid,'fesom_c_mesh',NF_INT,0,nf_dims,filelist(i)%varid%mesh)
         call check_nferr(iost, 'Error at nf90mpi_def_var ')
!         ! add attributes to mesh
!         n=2
!!IK it is needed, should be switch on after check         !iost = nf_put_att(filelist(i)%fileid,filelist(i)%varid%mesh,'topology_dimension',NF_INT,1,n)
!         call check_nferr(iost,filelist(i)%name)
         mpikind = len_trim('mesh_topology')
         iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%mesh,'cf_role',mpikind,'mesh_topology')
         call check_nferr(iost, 'Error at nf90mpi_put_att_text ')
         mpikind = len_trim('lon lat')
         iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%mesh,'node_coordinates',mpikind,'lon lat')
         call check_nferr(iost, 'Error at nf90mpi_put_att_text ')
         mpikind = len_trim('nv')
         iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%mesh,'face_node_connectivity',mpikind,'nv')
         call check_nferr(iost, 'Error at nf90mpi_put_att_text ')
!         !  define coordinates on nodes
         nf_dims(1) = filelist(i)%dimid%node
         iost = nfmpi_def_var(filelist(i)%fileid,'lon',NF_DOUBLE,1,nf_dims,filelist(i)%varid%lon)
         call check_nferr(iost, 'Error at nf90mpi_def_var ')
         iost = nfmpi_def_var(filelist(i)%fileid,'lat',NF_DOUBLE,1,nf_dims,filelist(i)%varid%lat)
         call check_nferr(iost, filelist(i)%name)
!         ! define depth
         iost = nfmpi_def_var(filelist(i)%fileid,'depth',NF_DOUBLE,1,nf_dims,filelist(i)%varid%depth)
         call check_nferr(iost,filelist(i)%name)

         nf_dims(1) = filelist(i)%dimid%nelem
         iost = nfmpi_def_var(filelist(i)%fileid,'depth_elem',NF_DOUBLE,1,nf_dims,filelist(i)%varid%depth_elem)
         call check_nferr(iost,filelist(i)%name)

         ! add attributes
         if (cartesian) then
            mpikind = len_trim('X')
            iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%lon,'standard_name',mpikind,'X')
            call check_nferr(iost,filelist(i)%name)
            mpikind = len_trim('meters')
            iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%lon,'units',mpikind,'meters')
            call check_nferr(iost,filelist(i)%name)
            mpikind = len_trim('Y')
            iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%lat,'standard_name',mpikind,'Y')
            call check_nferr(iost,filelist(i)%name)
            mpikind = len_trim('meters')
            iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%lat,'units',mpikind,'meters')
            call check_nferr(iost,filelist(i)%name)
         else
            mpikind = len_trim('longitude')
            iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%lon,'standard_name',mpikind,'longitude')
            call check_nferr(iost,filelist(i)%name)
            mpikind = len_trim('degrees_east')
            iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%lon,'units',mpikind,'degrees_east')
            call check_nferr(iost,filelist(i)%name)
            mpikind = len_trim('latitude')
            iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%lat,'standard_name',mpikind,'latitude')
            call check_nferr(iost,filelist(i)%name)
            mpikind = len_trim('degrees_east')
            iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%lat,'units',mpikind,'degrees_north')
            call check_nferr(iost,filelist(i)%name)
         end if
!         !depth
         mpikind = len_trim('depth')
         iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%depth,'standard_name',mpikind,'depth')
         call check_nferr(iost,filelist(i)%name)
         mpikind = len_trim('m')
         iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%depth,'units',mpikind,'m')
         call check_nferr(iost,filelist(i)%name)
         mpikind = len_trim('depth on elem')
         iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%depth_elem,'standard_name',&
                                                                   mpikind,'depth on elem')
         call check_nferr(iost,filelist(i)%name)
         mpikind = len_trim('m')
         iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%depth_elem,'units',mpikind,'m')
         call check_nferr(iost,filelist(i)%name)
         ! define coordinates on elements
         nf_dims(1) = filelist(i)%dimid%nelem
         iost = nfmpi_def_var(filelist(i)%fileid,'lon_elem',NF_DOUBLE,1,nf_dims,filelist(i)%varid%lon_elem)
         call check_nferr(iost,filelist(i)%name)
         iost = nfmpi_def_var(filelist(i)%fileid,'lat_elem',NF_DOUBLE,1,nf_dims,filelist(i)%varid%lat_elem)
         call check_nferr(iost,filelist(i)%name)
!         ! add attributes
         if (cartesian) then
            mpikind = len_trim('X_elements')
            iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%lon_elem,'standard_name', &
                                      & mpikind,'X_elements')
            call check_nferr(iost,filelist(i)%name)
            mpikind = len_trim('meters')
            iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%lon_elem,'units', &
                                      & mpikind,'meters')
            call check_nferr(iost,filelist(i)%name)
            mpikind = len_trim('Y_elements')
            iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%lat_elem,'standard_name', &
                                      & mpikind,'Y_elements')
            call check_nferr(iost,filelist(i)%name)
            mpikind = len_trim('meters')
            iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%lat_elem,'units', &
                                      & mpikind,'meters')
            call check_nferr(iost,filelist(i)%name)
         else
            mpikind = len_trim('longitude_elements')
            iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%lon_elem,'standard_name', &
                                      & mpikind,'longitude_elements')
            call check_nferr(iost,filelist(i)%name)
            mpikind = len_trim('degrees_east')
            iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%lon_elem,'units', &
                                  & mpikind,'degrees_east')
            call check_nferr(iost,filelist(i)%name)
            mpikind = len_trim('latitude_elements')
            iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%lat_elem,'standard_name', &
                                  & mpikind,'latitude_elements')
            call check_nferr(iost,filelist(i)%name)
            mpikind = len_trim('degrees_east')
            iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%lat_elem,'units', &
                                  & mpikind,'degrees_north')
            call check_nferr(iost,filelist(i)%name)
         end if
         ! define time
         nf_dims(1)=filelist(i)%dimid%time
         iost = nfmpi_def_var(filelist(i)%fileid,'time',NF_DOUBLE,1,nf_dims,filelist(i)%varid%time)
         call check_nferr(iost,filelist(i)%name)
         mpikind = len_trim('time')
         iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%time,'standard_name',mpikind,'time')
         call check_nferr(iost,filelist(i)%name)

         !  define elements area, nods area and weights
         nf_dims(1) = filelist(i)%dimid%node
         iost = nfmpi_def_var(filelist(i)%fileid,'area',NF_DOUBLE,1,nf_dims,filelist(i)%varid%area)
         call check_nferr(iost,filelist(i)%name)

         nf_dims(1) = filelist(i)%dimid%nelem
         iost = nfmpi_def_var(filelist(i)%fileid,'elem_area',NF_DOUBLE,1,nf_dims,filelist(i)%varid%elem_area)
         call check_nferr(iost,filelist(i)%name)

         nf_dims(1) = filelist(i)%dimid%maxnod
         nf_dims(2) = filelist(i)%dimid%nelem
         iost = nfmpi_def_var(filelist(i)%fileid,'w_cv', NF_DOUBLE, 2,nf_dims,filelist(i)%varid%w_cv)
         call check_nferr(iost,filelist(i)%name)
         ! define nod_in_elem2D and nod_in_elem2D_num
         nf_dims(1) = filelist(i)%dimid%node
         iost = nfmpi_def_var(filelist(i)%fileid,'nod_in_elem2d_num', NF_INT, 1,nf_dims,filelist(i)%varid%nod_in_elem2d_num)
         call check_nferr(iost,filelist(i)%name)
         nf_dims(1) = filelist(i)%dimid%maxnodelem
         nf_dims(2) = filelist(i)%dimid%node
         iost = nfmpi_def_var(filelist(i)%fileid,'nod_in_elem2d', NF_INT, 2,nf_dims,filelist(i)%varid%nod_in_elem2d)
         call check_nferr(iost,filelist(i)%name)
         ! add attributes
         mpikind = len_trim('cell area')
         iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%area,&
                                'standard_name',mpikind,'cell area')
         call check_nferr(iost,filelist(i)%name)
         mpikind = len_trim('m**2')
         iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%area,&
                                'units',mpikind,'m**2')
         call check_nferr(iost,filelist(i)%name)
         mpikind = len_trim('elemens area')
         iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%elem_area,&
                                'standard_name',mpikind,'elemens area')
         call check_nferr(iost,filelist(i)%name)
         mpikind = len_trim('m**2')
         iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%elem_area,&
                                'units',mpikind,'m**2')
         call check_nferr(iost,filelist(i)%name)
         mpikind = len_trim(' Scv/Sc')
         iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%w_cv,&
                                'standard_name',mpikind,' Scv/Sc')
         call check_nferr(iost,filelist(i)%name)
         mpikind = len_trim('node is in N elemets')
         iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%nod_in_elem2D_num,&
                                'standard_name',mpikind,'node is in N elemets')
         call check_nferr(iost,filelist(i)%name)
         mpikind = len_trim('node is in elemets')
         iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%nod_in_elem2D,&
                                'standard_name',mpikind,'node is in elemets')
         call check_nferr(iost,filelist(i)%name)
         str='days since 1900-01-01 00:00:00'
         mpikind = len_trim(str)
         iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%time,'units', &
                        & mpikind ,str)
         call check_nferr(iost,filelist(i)%name)
         mpikind = len_trim('T')
         iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%time,'axis',mpikind,'T')
         call check_nferr(iost,filelist(i)%name)
         mpikind = len_trim('Time')
         iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%time,'title',mpikind,'Time')
         call check_nferr(iost,filelist(i)%name)
!      iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%time,'calendar',len_trim('julian'),'julian')
!      call check_nferr(iost,filelist(i)%name)

         !  define nods and elements indexes
         nf_dims(1) = filelist(i)%dimid%node
         iost = nfmpi_def_var(filelist(i)%fileid,'node_idx',NF_INT,1,nf_dims,filelist(i)%varid%node)
         call check_nferr(iost,filelist(i)%name)
         mpikind = len_trim('nodes indexes')
         iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%node,'standard_name',mpikind,'nodes indexes')
         call check_nferr(iost,filelist(i)%name)

         nf_dims(1) = filelist(i)%dimid%nelem
         iost = nfmpi_def_var(filelist(i)%fileid,'elem_idx',NF_INT,1,nf_dims,filelist(i)%varid%elem)
         call check_nferr(iost,filelist(i)%name)
         mpikind = len_trim('elements indexes')
         iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varid%elem,'standard_name',mpikind,'elements indexes')
         call check_nferr(iost,filelist(i)%name)

         do j = 1, filelist(i)%numvars
            varid = filelist(i)%vars(j)
            ndim=0
            select case (varlist(varid)%vertdim)
               case (0)
               case (1)
                  ndim = ndim+1
                  nf_dims(ndim)=filelist(i)%dimid%sigma
               case (2)
                  ndim = ndim+1
                  nf_dims(ndim)=filelist(i)%dimid%sigmam1
            end select
            select case (varlist(varid)%onnode)
               case (-1)
               case (0)
                  ndim = ndim+1
                  nf_dims(ndim)=filelist(i)%dimid%nelem
               case (1)
                  ndim = ndim+1
                  nf_dims(ndim)=filelist(i)%dimid%node
            end select
            select case (varlist(varid)%timedim)
               case (0)
               case (1)
                  ndim = ndim+1
                  nf_dims(ndim)=filelist(i)%dimid%time
            end select

            iost = nfmpi_def_var(filelist(i)%fileid,varlist(varid)%name,&
                                  NF_DOUBLE,ndim,nf_dims,filelist(i)%varsid(j))
            call check_nferr(iost,filelist(i)%name)
            mpikind = len_trim(varlist(varid)%units)
            iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varsid(j),&
                   'units',mpikind,varlist(varid)%units)
            call check_nferr(iost,filelist(i)%name)
            mpikind = len_trim(varlist(varid)%standard_name)
            iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varsid(j),&
                   'standard_name',mpikind,varlist(varid)%standard_name)
            call check_nferr(iost,filelist(i)%name)
            mpikind = len_trim(varlist(varid)%modelname)
            iost = nfmpi_put_att_text(filelist(i)%fileid,filelist(i)%varsid(j),&
                   'modelname',mpikind,varlist(varid)%modelname)
            call check_nferr(iost,filelist(i)%name)


         end do
         ! leaving define mode
         iost = nf90mpi_enddef(filelist(i)%fileid)
         call check_nferr(iost, 'Error at nf90mpi_enddef ')
!         ! write down coordinates (nodes)
         nf_start(1)=filelist(i)%mynodstart
         nf_edges(1)=filelist(i)%nnods
         if (cartesian) then
            nod_x = coord_nod2d(1,1:myDim_nod2D)*r_earth
            nod_y = coord_nod2d(2,1:myDim_nod2D)*r_earth
         else
            nod_x = coord_nod2d(1,1:myDim_nod2D)/rad
            nod_y = coord_nod2d(2,1:myDim_nod2D)/rad
         end if
         if ((filelist(i)%ftype==3 .and. filelist(i)%mystation==1)) then
            !write coordinates for station
            tmpx = nod_x(filelist(i)%mynods(1))
            iost = nfmpi_put_vara_double_all(filelist(i)%fileid,filelist(i)%varid%lon,&
                            nf_start,nf_edges,tmpx(1))
            call check_nferr(iost,filelist(i)%name)
            tmpy = nod_y(filelist(i)%mynods(1))
            iost = nfmpi_put_vara_double_all(filelist(i)%fileid,filelist(i)%varid%lat,&
                            nf_start,nf_edges,tmpy(1))
            call check_nferr(iost,filelist(i)%name)
         else
            !for all other files
            iost = nfmpi_put_vara_double_all(filelist(i)%fileid,filelist(i)%varid%lon,nf_start,nf_edges,&
                                          nod_x(1:filelist(i)%nnods))
            call check_nferr(iost,filelist(i)%name)
            iost = nfmpi_put_vara_double_all(filelist(i)%fileid,filelist(i)%varid%lat,nf_start,nf_edges,&
                                          nod_y(1:filelist(i)%nnods))
            call check_nferr(iost,filelist(i)%name)
         end if

         ! write depth
         iost = nfmpi_put_vara_double_all(filelist(i)%fileid,filelist(i)%varid%depth,nf_start,nf_edges,&
                                          depth(1:filelist(i)%nnods))
         call check_nferr(iost,filelist(i)%name)
!print *, mype,filelist(i)%nelems,filelist(i)%myelemstart,felem2d, myDim_elem2D
         nf_start(1)=filelist(i)%myelemstart
         nf_edges(1)=filelist(i)%nelems
         iost = nfmpi_put_vara_double_all(filelist(i)%fileid,filelist(i)%varid%depth_elem,nf_start,nf_edges,&
                                          depth_elem(1:filelist(i)%nelems))
         call check_nferr(iost,filelist(i)%name)
!         ! write down coordinates (elements)
         if ((filelist(i)%ftype==3 .and. filelist(i)%mystation==1)) then
            !write coordinates for station
            tmpx = x(filelist(i)%myelems(1))
            tmpy = y(filelist(i)%myelems(1))
            iost = nfmpi_put_vara_double_all(filelist(i)%fileid,filelist(i)%varid%lon_elem,nf_start,nf_edges,&
                                            tmpx(1:filelist(i)%nelems))
            call check_nferr(iost,filelist(i)%name)
            iost = nfmpi_put_vara_double_all(filelist(i)%fileid,filelist(i)%varid%lat_elem,nf_start,nf_edges,&
                                            tmpy(1:filelist(i)%nelems))
            call check_nferr(iost,filelist(i)%name)

         else
            iost = nfmpi_put_vara_double_all(filelist(i)%fileid,filelist(i)%varid%lon_elem,nf_start,nf_edges,&
                                            x(1:filelist(i)%nelems))
            call check_nferr(iost,filelist(i)%name)
            iost = nfmpi_put_vara_double_all(filelist(i)%fileid,filelist(i)%varid%lat_elem,nf_start,nf_edges,&
                                            y(1:filelist(i)%nelems))
            call check_nferr(iost,filelist(i)%name)
         endif
!         ! write down NV (
!         !define nv elem2D_nodes nod_in_elem2D_inglob
         ALLOCATE(elem2D_nodes_inglob(4,myDim_elem2D), &
                   &      STAT=output_alloc )
         if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays elem2D_nodes_inglob(1:4,myDim_elem2D)'
         do j=1,myDim_elem2D
            elem2D_nodes_inglob(:,j) = myList_nod2D(elem2D_nodes(:,j))
         end do
         nf_start(2)=filelist(i)%myelemstart
         nf_edges(2)=filelist(i)%nelems
         nf_start(1)=1
         nf_edges(1)=4
         iost = nfmpi_put_vara_int_all(filelist(i)%fileid,filelist(i)%varid%nv,nf_start,nf_edges,&
                                        elem2D_nodes_inglob(1:4,1:filelist(i)%nelems))
         call check_nferr(iost,filelist(i)%name)
         deallocate(elem2D_nodes_inglob)

         iost = nfmpi_put_vara_double_all(filelist(i)%fileid,filelist(i)%varid%w_cv,nf_start,nf_edges,&
                                        w_cv(1:4,1:filelist(i)%nelems))
         call check_nferr(iost,filelist(i)%name)

!         !  write elements area, nods area and weights
         nf_start(1)=filelist(i)%mynodstart
         nf_edges(1)=filelist(i)%nnods
         iost = nfmpi_put_vara_double_all(filelist(i)%fileid,filelist(i)%varid%area,nf_start,nf_edges,&
                                                         area(1:filelist(i)%nnods))
         call check_nferr(iost,filelist(i)%name)

         iost = nfmpi_put_vara_int_all(filelist(i)%fileid,filelist(i)%varid%nod_in_elem2D_num,nf_start,nf_edges,&
                                    nod_in_elem2D_num(1:filelist(i)%nnods))
         call check_nferr(iost,filelist(i)%name)

         nf_start(1)=filelist(i)%myelemstart
         nf_edges(1)=filelist(i)%nelems
         iost = nfmpi_put_vara_double_all(filelist(i)%fileid,filelist(i)%varid%elem_area,nf_start,nf_edges,&
                                        elem_area(1:filelist(i)%nelems))
         call check_nferr(iost,filelist(i)%name)
         nf_start(2)=filelist(i)%mynodstart
         nf_edges(2)=filelist(i)%nnods
         nf_start(1)=1
         nf_edges(1)=maxval(nod_in_elem2D_num)
         iost = nfmpi_put_vara_int_all(filelist(i)%fileid,filelist(i)%varid%nod_in_elem2D,nf_start,nf_edges,&
                                    nod_in_elem2D(1:nf_edges(1),1:filelist(i)%nnods))
         call check_nferr(iost,filelist(i)%name)

         nf_start(1)=filelist(i)%mynodstart
         nf_edges(1)=filelist(i)%nnods
         iost = nfmpi_put_vara_int_all(filelist(i)%fileid,filelist(i)%varid%node,nf_start,nf_edges,&
                                                         filelist(i)%nods(1:filelist(i)%nnods))
         call check_nferr(iost,filelist(i)%name)
         nf_start(1)=filelist(i)%myelemstart
         nf_edges(1)=filelist(i)%nelems
         iost = nfmpi_put_vara_int_all(filelist(i)%fileid,filelist(i)%varid%elem,nf_start,nf_edges,&
                                        filelist(i)%elems(1:filelist(i)%nelems))
         call check_nferr(iost,filelist(i)%name)

         ! close file
         iost = nf90mpi_close(filelist(i)%fileid)
         call check_nferr(iost,filelist(i)%name)
         ! define 0 records (time) in nc outfile
         filelist(i)%nc_rec = 0
         ! write sigma to file if PE=0 or if file is a station at current PE
         if ((mype==0 .and. filelist(i)%ftype==0) &
              .or. (filelist(i)%ftype==3 .and. filelist(i)%mystation==1)) then
             if (type_task>1) then
                iost = nf90mpi_open(MPI_COMM_SELF, filelist(i)%name, NF90_WRITE, &
                                      MPI_INFO_NULL, filelist(i)%fileid)
                call check_nferr(iost,filelist(i)%name)
                ! sigma
                nf_start(1)=1
                nf_edges(1)=nsigma
                iost = nfmpi_put_vara_double_all(filelist(i)%fileid,filelist(i)%varid%sigma,nf_start,nf_edges,sigma)
                call check_nferr(iost,filelist(i)%name)
                iost = nf90mpi_close(filelist(i)%fileid)
                call check_nferr(iost,filelist(i)%name)
             end if
          end if

      end do

      DEALLOCATE( x, y, nod_x, nod_y)

      ALLOCATE( fid1d_size(maxnumberallvars,numfiles), fid2d_size(maxnumberallvars,numfiles), &
                   &      STAT=output_alloc )
      if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays fid1d_size,fid2d_size'

      ! count number of variables of each type
      nvar0d = 0
      nvar1d = 0
      nvar2d = 0

      do i = 1, numfiles
         if (filelist(i)%ftype==3 .and. filelist(i)%mystation==0) then
            cycle
         end if
         do j = 1, filelist(i)%numvars
            filelist(i)%fid(j) = 0
            filelist(i)%avgid(j) = 0
            select case (varlist(filelist(i)%vars(j))%ndim) ! select proper dimension
               case(1)
                  ! find if variable is in the list already
                  ! if not then add
                  ! link var in file to var in a list
                  loc = 0
                  do k = 1, nvar0d
                     if ( filelist(i)%vars(j) == var0d(k) ) then
                        loc = k
                        exit
                     end if
                  end do
                  if ( loc >0 ) then ! if variable is in the list
                     filelist(i)%fid(j) = loc
                     filelist(i)%avgid(j) = numfid0d(loc) + 1
                     numfid0d(loc) = numfid0d(loc) + 1
                  else   ! add one more variable in the list for averaging
                     nvar0d = nvar0d + 1
                     filelist(i)%fid(j) = nvar0d
                     var0d(nvar0d) = filelist(i)%varsid(j)
                     numfid0d(nvar0d) = 1
                     filelist(i)%avgid(j) = 1
                  endif
               case(2)
                  loc = 0
                  do k = 1, nvar1d
                     if ( filelist(i)%vars(j) == var1d(k) ) then
                        loc = k
                        exit
                     end if
                  end do
                  if ( loc >0 ) then ! if variable is in the list
                     filelist(i)%fid(j) = loc
                     filelist(i)%avgid(j) = numfid1d(loc) + 1
                     numfid1d(loc) = numfid1d(loc) + 1
                  else   ! add one more variable in the list for averaging
                     nvar1d = nvar1d + 1
                     filelist(i)%fid(j) = nvar1d
                     var1d(nvar1d) = filelist(i)%vars(j)
                     numfid1d(nvar1d) = 1
                     filelist(i)%avgid(j) = 1
                     loc = nvar1d
                  endif

                  if (varlist(var1d(loc))%onnode == 1 ) then
                     fid1d_size(loc,numfid1d(loc)) = filelist(i)%nnods
                  else
                     fid1d_size(loc,numfid1d(loc)) = filelist(i)%nelems
                  endif
               case(3)
                  loc = 0
                  do k = 1, nvar2d
                     if ( filelist(i)%vars(j) == var2d(k) ) then
                        loc = k
                        exit
                     end if
                  end do
                  if ( loc >0 ) then ! if variable is in the list
                     filelist(i)%fid(j) = loc
                     filelist(i)%avgid(j) = numfid2d(loc) + 1
                     numfid2d(loc) = numfid2d(loc) + 1
                  else   ! add one more variable in the list for averaging
                     nvar2d = nvar2d + 1
                     filelist(i)%fid(j) = nvar2d
                     var2d(nvar2d) = filelist(i)%vars(j)
                     numfid2d(nvar2d) = 1
                     filelist(i)%avgid(j) = 1
                     loc = nvar2d
                  endif
                  if (varlist(var2d(loc))%onnode == 1 ) then
                     fid2d_size(loc,numfid2d(loc)) = filelist(i)%nnods
                  else
                     fid2d_size(loc,numfid2d(loc)) = filelist(i)%nelems
                  endif
            end select
         end do
      end do
      if (mype==0) then
         write(*,*) "ncoutputmpi_ini: myPE=",mype
         write(*,*) "            Number of 0D variables for output:", nvar0d
         write(*,*) "            Number of 1D variables for output:", nvar1d
         write(*,*) "            Number of 2D variables for output:", nvar2d
      end if

      ! allocate memory for list of variables for output
      ! 3 lists for each dimension (time is a dimension and x-y plane is one dimension)
      if (nvar0d > 0) then
         ALLOCATE( varbas0d(nvar0d), &
                   &      STAT=output_alloc )
         if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays varbas0d'
         do i = 1, nvar0d
            varbas0d(i)%varid = var0d(i)
            ALLOCATE( varbas0d(i)%avg(numfid0d(i)), &
                      &      STAT=output_alloc )
            if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays varbas0d(i)%avg'
            do j = 1, numfid0d(i)
               varbas0d(i)%avg(j)%v = 0
            enddo
            varbas0d(i)%nfid = numfid0d(i)
         end do
      end if

      if (nvar1d > 0) then
         ALLOCATE( varbas1d(nvar1d), &
                   &      STAT=output_alloc )
         if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays varbas1d '
         do i = 1, nvar1d
            !
            varbas1d(i)%varid = var1d(i)
            !allocate memory for snapshot of every variable and
            ALLOCATE(varbas1d(i)%snapshot(varlist(var1d(i))%varsize(1)), &
                      & varbas1d(i)%avg(numfid1d(i) ), &
                      &      STAT=output_alloc )
            if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays varbas1d(i)%snapshot,varbas1d(i)%avg'

            do j = 1, numfid1d(i)
               ALLOCATE( varbas1d(i)%avg(j)%v(fid1d_size(i,j)), &
                         varbas1d(i)%avg(j)%sn(fid1d_size(i,j)), &
                      &      STAT=output_alloc )
               if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays varbas1d(i)%avg(j)%v,varbas1d(i)%avg(j)%sn'
               varbas1d(i)%avg(j)%sz = fid1d_size(i,j)
               varbas1d(i)%avg(j)%v = 0
            enddo

            varbas1d(i)%nfid = numfid1d(i)
         end do
      end if

      if (nvar2d > 0) then
         ALLOCATE( varbas2d(nvar2d), &
                   &      STAT=output_alloc )
         if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays varbas2d'
         do i = 1, nvar2d
            varbas2d(i)%varid = var2d(i)
            ALLOCATE( varbas2d(i)%snapshot(varlist(var2d(i))%varsize(1), varlist(var2d(i))%varsize(2) ),  &
                      varbas2d(i)%avg(numfid2d(i)), &
                        &      STAT=output_alloc )
            if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays varbas2d(i)%snapshot, varbas2d(i)%avg'
            do j = 1, numfid2d(i)
               ALLOCATE( varbas2d(i)%avg(j)%v(varlist(var2d(i))%varsize(1),   &
                               fid2d_size(i,j) ),  &
                         varbas2d(i)%avg(j)%sn(varlist(var2d(i))%varsize(1),   &
                               fid2d_size(i,j) ),  &
                        &      STAT=output_alloc )
               if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays varbas2d(i)%avg(j)%v, varbas2d(i)%avg(j)%sn'
               varbas2d(i)%avg(j)%sz(1) = varlist(var2d(i))%varsize(1)
               varbas2d(i)%avg(j)%sz(2) = fid2d_size(i,j)
               varbas2d(i)%avg(j)%v = 0
            enddo
            varbas2d(i)%nfid = numfid2d(i)
         end do
      end if

      DEALLOCATE( fid1d_size, fid2d_size)

   END SUBROUTINE ncoutputmpi_ini

   SUBROUTINE getvarsize(varid,varsize)

      IMPLICIT NONE
      integer, intent(in)  :: varid  ! index of variable in varlist
      integer, intent(out) :: varsize(2)   ! size of variable
      integer              :: ndim
      ndim=0
      varsize(1) = 0
      varsize(2) = 0
      select case (varlist(varid)%vertdim)
         case (0)
         case (1)
            ndim = ndim+1
            varsize(ndim) = nsigma
         case (2)
            ndim = ndim+1
            varsize(ndim) = nsigma-1
      end select
      select case (varlist(varid)%onnode)
         case (-1)
         case (0)
            ndim = ndim+1
            varsize(ndim) = myDim_elem2D+eDim_elem2D+eXDim_elem2D
         case (1)
            ndim = ndim+1
            varsize(ndim) = myDim_nod2D+eDim_nod2D
      end select

   END SUBROUTINE getvarsize


   SUBROUTINE write_data(findx)
      IMPLICIT NONE

      integer, intent(in)  :: findx  ! index of file to write in filelist
      integer              :: i
      real(kind=WP)        :: tshift
      real(kind=WP)        :: newtime(1)
      integer              :: requests(maxnumbervars), statuses(maxnumbervars)
      !increase time record for ncout
      filelist(findx)%nc_rec = filelist(findx)%nc_rec + 1

      !open for writing netCDF file with 0 Pe or if station at current PE
      if (((mype==0) .and. filelist(findx)%ftype==0) &
             .or. (filelist(findx)%ftype==3 .and. filelist(findx)%mystation==1)) then
          iost = nf90mpi_open(MPI_COMM_SELF, filelist(findx)%name, NF90_WRITE, &
                                          MPI_INFO_NULL, filelist(findx)%fileid)
          call check_nferr(iost,filelist(findx)%name)

          ! write time
          nf_start(1)=filelist(findx)%nc_rec
          nf_edges(1)=1



          ! time_jd - time_jd0 - is days since begining of simulations
          ! changed to 2415020.5 (1900-01-01 00:00)
          tshift = 0.0_WP
          if (filelist(findx)%snapshot == 2) then
            tshift = filelist(findx)%freq*dt*0.5_WP/86400.0_WP
          endif
          newtime(1) = time_jd - 2415020.5 - tshift
          iost = nfmpi_put_vara_double_all(filelist(findx)%fileid,&
                  & filelist(findx)%varid%time,nf_start,nf_edges, newtime) !time_jd0)
          call check_nferr(iost,filelist(findx)%name)
          iost = nf90mpi_close(filelist(findx)%fileid)
          call check_nferr(iost,filelist(findx)%name)
      end if

      ! write variables one by one
      !OPEN
      if (filelist(findx)%ftype==0) then
        iost = nf90mpi_open(MPI_COMM_WORLD, filelist(findx)%name, NF90_WRITE, &
                          MPI_INFO_NULL, filelist(findx)%fileid)
        call check_nferr(iost,filelist(findx)%name)
      else
         iost = nf90mpi_open(MPI_COMM_SELF, filelist(findx)%name, NF90_WRITE, &
                                          MPI_INFO_NULL, filelist(findx)%fileid)
         call check_nferr(iost,filelist(findx)%name)
      end if
      do i = 1, filelist(findx)%numvars
         call write_var(findx,i,requests(i))
      end do

      iost = nf90mpi_wait_all(filelist(findx)%fileid, filelist(findx)%numvars, requests, statuses)

      call check_nferr(iost,filelist(findx)%name)
      do i = 1, filelist(findx)%numvars
         call check_nferr(statuses(i),filelist(findx)%name)
      end do

!      ! close file
      iost = nf90mpi_close(filelist(findx)%fileid)
      call check_nferr(iost,filelist(findx)%name)

   END SUBROUTINE write_data

   SUBROUTINE write_var(findx,vindx,req)
    ! write each variable to file
    ! here is DANGER, we use netcdf IDs declared during creation of file,
    ! it could be, by opening files in different order we could mix it somehow
    ! I have no idea why it worked before with serial netcdf (but it worked quiet well)
      IMPLICIT NONE

      integer, intent(in)  :: findx  ! index of file to write in filelist
      integer, intent(in)  :: vindx  ! index of variable to write
      integer, intent(out) :: req    ! PnetCDF request ID, (for non-blocking write interface )
      integer              :: i, ndim, varid
      real(kind=WP)        :: tmpreal(1)


      varid = filelist(findx)%vars(vindx)
      ndim=0
      select case (varlist(varid)%vertdim)
         case (0)
         case (1)
            ndim = ndim+1
            nf_start(ndim)=1
            nf_edges(ndim)=nsigma
         case (2)
            ndim = ndim+1
            nf_start(ndim)=1
            nf_edges(ndim)=nsigma-1
      end select
      select case (varlist(varid)%onnode)
         case (-1)
         case (0)
            ndim = ndim+1
            nf_start(ndim)=filelist(findx)%myelemstart
            nf_edges(ndim)=filelist(findx)%nelems
         case (1)
            ndim = ndim+1
            nf_start(ndim)=filelist(findx)%mynodstart
            nf_edges(ndim)=filelist(findx)%nnods
      end select
      select case (varlist(varid)%timedim)
         case (0)
         case (1)
            ndim = ndim+1
            nf_start(ndim)=filelist(findx)%nc_rec
            nf_edges(ndim)=1
      end select

      if (filelist(findx)%snapshot == 1) then
        ! if we write snapshots
         select case (varlist(filelist(findx)%vars(vindx))%ndim) ! select proper dimension
            case(1)
               tmpreal(1) = varbas0d(filelist(findx)%fid(vindx))%avg(filelist(findx)%avgid(vindx))%sn
               iost = nfmpi_iput_vara_double(filelist(findx)%fileid,filelist(findx)%varsid(vindx),nf_start,nf_edges,&
                   & tmpreal(1), req )
            case(2)
               iost = nfmpi_iput_vara_double(filelist(findx)%fileid,filelist(findx)%varsid(vindx),nf_start,nf_edges,&
                   & varbas1d(filelist(findx)%fid(vindx))%avg(filelist(findx)%avgid(vindx))%sn , req)
            case(3)
               iost = nfmpi_iput_vara_double(filelist(findx)%fileid,filelist(findx)%varsid(vindx),nf_start,nf_edges,&
                   & varbas2d(filelist(findx)%fid(vindx))%avg(filelist(findx)%avgid(vindx))%sn, req )
         end select
      else
        ! if we write averaged values
         select case (varlist(filelist(findx)%vars(vindx))%ndim) ! select proper dimension
            case(1)
               tmpreal(1) = varbas0d(filelist(findx)%fid(vindx))%avg(filelist(findx)%avgid(vindx))%v
               iost = nfmpi_iput_vara_double(filelist(findx)%fileid,filelist(findx)%varsid(vindx),nf_start,nf_edges,&
                   &  tmpreal(1), req)
            case(2)
               iost = nfmpi_iput_vara_double(filelist(findx)%fileid,filelist(findx)%varsid(vindx),nf_start,nf_edges,&
                   & varbas1d(filelist(findx)%fid(vindx))%avg(filelist(findx)%avgid(vindx))%v, req )
            case(3)
               iost = nfmpi_iput_vara_double(filelist(findx)%fileid,filelist(findx)%varsid(vindx),nf_start,nf_edges,&
                   & varbas2d(filelist(findx)%fid(vindx))%avg(filelist(findx)%avgid(vindx))%v, req )
         end select
      end if
      call check_nferr(iost,filelist(findx)%name)

   END SUBROUTINE write_var

   SUBROUTINE getvar
      IMPLICIT NONE

      integer   :: i

      do i = 1, nvar0d
         call getvar0d(i)
      end do
      do i = 1, nvar1d
         call getvar1d(i)
      end do
      do i = 1, nvar2d
         call getvar2d(i)
      end do

   END SUBROUTINE getvar

   SUBROUTINE addvar2basket
      IMPLICIT NONE

      integer   :: i, j
      integer   :: fid, avgid, vid

      do i = 1, nvar0d
         call getvar0d(i)
      end do
      do i = 1, nvar1d
         call getvar1d(i)
      end do
      do i = 1, nvar2d
         call getvar2d(i)
      end do
! this is very complicated way to do simple things
      do i = 1, numfiles
         if (filelist(i)%ftype==3 .and. filelist(i)%mystation==0) then
            cycle
         end if
         do j = 1, filelist(i)%numvars
            fid   = filelist(i)%fid(j)
            avgid = filelist(i)%avgid(j)
            vid   = filelist(i)%vars(j)
            select case (varlist(vid)%ndim) ! select proper dimension
               case(1)
                     varbas0d(fid)%avg(avgid)%sn = varbas0d(fid)%snapshot
               case(2)
                  if (varlist(vid)%onnode == 1 ) then
!                     varbas1d(fid)%avg(avgid)%sn = varbas1d(fid)%snapshot(1:filelist(i)%nnods)
                     varbas1d(fid)%avg(avgid)%sn = varbas1d(fid)%snapshot(filelist(i)%mynods)
                  else
!                     varbas1d(fid)%avg(avgid)%sn = varbas1d(fid)%snapshot(1:filelist(i)%nelems)
                     varbas1d(fid)%avg(avgid)%sn = varbas1d(fid)%snapshot(filelist(i)%myelems)
                  end if
               case(3)
                  if (varlist(vid)%onnode == 1 ) then
!                     varbas2d(fid)%avg(avgid)%sn = varbas2d(fid)%snapshot(:,1:filelist(i)%nnods)
                     varbas2d(fid)%avg(avgid)%sn = varbas2d(fid)%snapshot(:,filelist(i)%mynods)
                  else
!                     varbas2d(fid)%avg(avgid)%sn = varbas2d(fid)%snapshot(:,1:filelist(i)%nelems)
                     varbas2d(fid)%avg(avgid)%sn = varbas2d(fid)%snapshot(:,filelist(i)%myelems)
                  end if
            end select

         end do
      end do

      do i = 1, nvar0d
         do j = 1, varbas0d(i)%nfid
            varbas0d(i)%avg(j)%v = varbas0d(i)%avg(j)%v + varbas0d(i)%avg(j)%sn
         end do
      end do
      do i = 1, nvar1d
         do j = 1, varbas1d(i)%nfid
            varbas1d(i)%avg(j)%v = varbas1d(i)%avg(j)%v + varbas1d(i)%avg(j)%sn
         end do
      end do
      do i = 1, nvar2d
         do j = 1, varbas2d(i)%nfid
            varbas2d(i)%avg(j)%v = varbas2d(i)%avg(j)%v + varbas2d(i)%avg(j)%sn
         end do
      end do

   END SUBROUTINE addvar2basket




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
   SUBROUTINE ne_nodeelem(x,y,ind,ind_elem,mystation)
!===================================================================================|
! Looking for a closest node to x, y and return node and element index.             |
! if station in PE, than mystation=1.                       After looking for       |
!                         first element in current PE in nod_in_elem2D              |
!   if the coordinates are given in Cartesian coord., then usual Eulerian norm is   |
!   used to calculate the distance.                                                 |
!   If the coordinates are given in spherical coord., than Haversine Formula is used|
!===================================================================================|
   IMPLICIT NONE
   real(kind=WP),intent(in) :: x, y
   integer,intent(out)      :: ind, mystation, ind_elem
   real(kind=WP)            :: mindist, d, mindist_glob
   integer                  :: n

   ind = 1
   ind_elem = 1
   mystation = 0
   mindist = 999999999.9
   if (cartesian) then
      do n = 1, myDim_nod2D
         d = (y/r_earth - coord_nod2D(2,n))**2 + (x/r_earth - coord_nod2D(1,n))**2
         d = sqrt(d)
         if ( mindist > d ) then
            ind = n
            mindist = d
         endif
      enddo
   else
      do n = 1, myDim_nod2D
         d = sin((y*rad - coord_nod2D(2,n))*0.5_WP)**2 + &
             cos(coord_nod2D(2,n)) * cos(y*rad) * sin((x*rad - coord_nod2D(1,n))*0.5_WP)**2
         d = asin(min(1.0_WP,sqrt(d)))
         if ( mindist > d ) then
            ind = n
            mindist = d
         endif
      enddo
   endif
   call MPI_AllREDUCE(mindist,mindist_glob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
          MPI_COMM_FESOM_C, MPIerr)
   if (mindist <= mindist_glob) then
      mystation = 1
      ind_elem = nod_in_elem2D(1,ind)
   end if
   END SUBROUTINE ne_nodeelem

   SUBROUTINE ne_element(x,y,ind)
!===================================================================================|
! Looking for a closest element to x, y and return element index                             |
!   if the coordinates are given in Cartesian coord., then usual Eulerian norm is   |
!   used to calculate the distance.                                                 |
!   If the coordinates are given in spherical coord., than Haversine Formula is used|
!===================================================================================|
   IMPLICIT NONE
   real(kind=WP),intent(in) :: x, y
   integer,intent(out)      :: ind
   real(kind=WP)            :: mindist, d
   integer                  :: n, elnodes(4)
   real(wp),dimension(elem2D)  :: xe, ye

      ! get coordinates of elements
   do n = 1, elem2D
      elnodes = elem2D_nodes(:,n) !! 4 nodes in element
      ye(n) = w_cv(1,n)*coord_nod2D(2,elnodes(1)) &
             + w_cv(2,n)*coord_nod2D(2,elnodes(2)) &
             + w_cv(3,n)*coord_nod2D(2,elnodes(3)) &
             + w_cv(4,n)*coord_nod2D(2,elnodes(4))
      xe(n) = w_cv(1,n)*coord_nod2D(1,elnodes(1)) &
             + w_cv(2,n)*coord_nod2D(1,elnodes(2)) &
             + w_cv(3,n)*coord_nod2D(1,elnodes(3)) &
             + w_cv(4,n)*coord_nod2D(1,elnodes(4))
   end do

   mindist = 999999999.9

   if (cartesian) then
      do n = 1, elem2D
         d = (y/r_earth - ye(n))**2 + (x/r_earth - xe(n))**2
         d = sqrt(d)
         if ( mindist > d ) then
            ind = n
            mindist = d
         endif
      enddo
   else
      do n = 1, elem2D
         d = sin((y*rad - ye(n))*0.5_WP)**2 + &
             cos(ye(n)) * cos(y*rad) * sin((x*rad - xe(n))*0.5_WP)**2
         d = asin(min(1.0_WP,sqrt(d)))
         if ( mindist > d ) then
            ind = n
            mindist = d
         endif
      enddo
   endif
   END SUBROUTINE ne_element

   subroutine check_nferr(errn, message)
   !use mpi
   use pnetcdf
   implicit none
   integer errn
   character(len=*) message
   ! It is a good idea to check returned value for possible error
   if (errn .NE. NF90_NOERR) then
      write(6,*) trim(message), trim(nf90mpi_strerror(errn))
       call MPI_Abort(MPI_COMM_WORLD, -1, errn)
   end if
   end subroutine check_nferr

subroutine get_fullstartpostitions(mynodstart,myelemstart,felem2D)

  implicit none
  integer,intent(out) :: mynodstart, myelemstart, felem2D
  integer, allocatable ::  pos_nod(:)
  integer, allocatable ::  pos_elem(:)

  integer :: n2D, e2D
  integer :: n, output_alloc

  if (npes > 1) then
     allocate(pos_nod(npes+1), &
                   &      STAT=output_alloc )
     if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays pos_nod'
     allocate(pos_elem(npes+1), &
                   &      STAT=output_alloc )
     if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays pos_elem'

     if (mype==0) then
        pos_nod(1) = 1
        pos_elem(1) = 1
        pos_nod(2) = myDim_nod2D+1
        pos_elem(2) = myDim_elem2D+1

        do n=1, npes-1
           call MPI_RECV(n2D, 1, MPI_INTEGER, n, 0, MPI_COMM_FESOM_C, MPI_STATUS_IGNORE, MPIerr )
           call MPI_RECV(e2D, 1, MPI_INTEGER, n, 1, MPI_COMM_FESOM_C, MPI_STATUS_IGNORE, MPIerr )

           pos_nod(n+2)  = pos_nod(n+1)  + n2D
           pos_elem(n+2) = pos_elem(n+1) + e2D
        enddo

     else
        call MPI_SEND(myDim_nod2D,   1,            MPI_INTEGER, 0, 0, MPI_COMM_FESOM_C, MPIerr )
        call MPI_SEND(myDim_elem2D,  1,            MPI_INTEGER, 0, 1, MPI_COMM_FESOM_C, MPIerr )
     end if

     call MPI_BCast(pos_nod, npes+1, MPI_INTEGER, 0, MPI_COMM_FESOM_C, MPIerr)
     call MPI_BCast(pos_elem, npes+1, MPI_INTEGER, 0, MPI_COMM_FESOM_C, MPIerr)

     mynodstart = pos_nod(mype+1)
     myelemstart = pos_elem(mype+1)
     felem2D = pos_elem(npes+1)-1

     deallocate(pos_nod)
     deallocate(pos_elem)
  endif

end subroutine get_fullstartpostitions

!   SUBROUTINE write_mesh2nc
!
!      IMPLICIT NONE
!

!   END SUBROUTINE write_mesh2nc

#endif
END MODULE fv_ncoutputmpi

