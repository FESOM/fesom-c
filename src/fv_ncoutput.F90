MODULE fv_ncoutput
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
   USE i_ARRAYS
   USE i_PARAM
   IMPLICIT NONE

   include 'netcdf.inc'

   public  output_ini    ! routine called before 1st time step (open files, read namelist,...)
   public  output_do     ! routine called each time step

   private
      integer      :: iost !I/O status
      integer      :: ncid ! netcdf file ID for temporal
      integer,save :: ncid_m ! netcdf file ID for output (Main (snapshots))
      integer,save :: nc_rec ! record in NC file (=0 at inizialization)

      ! temporal arrays
      integer      :: nf_dims(4) ! dimensions (temporal)
      integer      :: nf_start(4)
      integer      :: nf_edges(4)
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
         integer              :: nv, mesh, lon, lat, lon_elem, lat_elem, time, &
                                 sigma, sigmam1, elem_area, w_cv, area, &
                                 nod_in_elem2D, nod_in_elem2D_num, depth, depth_elem
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
         integer,allocatable,dimension(:) :: nods ! list of nodes were to put output (one for station,...)
         integer,allocatable,dimension(:) :: elems ! list of elements were to put output (one for station,...)
         integer                  :: nnods ! number of nodes were to put output (one for station,...)
         integer                  :: nelems ! number of elements were to put output (one for station,...)
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
      varlist(nvar)%standard_name =  'snu'! full variable name in the NetCDF file
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
      IMPLICIT NONE
      integer, intent(in) :: i !
      integer             :: n !

      select case (varbas1d(i)%varid)
         case(1)
            varbas1d(i)%snapshot = eta_n
            ! if we have wettingdrying and depth<0 then add depth to eta
            do n=1,nod2d
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
                varbas1d(i)%snapshot =  runoff_rivers
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

      end select
      ! multiply by mask for wetting drying
      if (varlist(varbas1d(i)%varid)%onnode == 1) then
         varbas1d(i)%snapshot = varbas1d(i)%snapshot * mask_wd_node
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
   SUBROUTINE output_do
      IMPLICIT NONE

      integer :: i


      call addvar2basket !add variables to varbas?d lists

      do i = 1, numfiles
         if ( mod(n_dt,filelist(i)%freq)==0 .OR. n_dt == 1) then
            ! if it is time to make output
               ! if we need to calculate average values
            if (filelist(i)%snapshot == 2) call calculateavg(i)

            call write_data(i) ! write values of snapshots/avg to nc file
             ! put average to 0
            if (filelist(i)%snapshot == 2) call zeroavg(i)

         end if
      end do

   END SUBROUTINE output_do

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

   SUBROUTINE output_ini
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE output_ini ***
      !!
      !! ** Purpose : inizialization of netcdf output
      !! ** Method  :
      !! ** Action  :
      !!----------------------------------------------------------------------
      IMPLICIT NONE


      real(wp), allocatable, dimension(:)  :: x, y, depth_elem
      real(wp)     :: mindist
      integer      :: i,j, k, n, numnodes, elnodes(4), output_alloc

      integer      :: filelist_unit

      character(len=30)   :: str
      integer            :: yyyy,mm,dd, ihh,imm,isec
      character(4) :: yystr
      character(2) :: mmstr,ddstr,hr,mn,sc
      character(len=130)   :: fpath

      integer              :: ndim, varid,  loc  ! temporal variables
      integer,dimension(maxnumberallvars)  :: var0d,var1d,var2d,numfid0d,numfid1d,numfid2d !temporal lists of variables
      integer, allocatable ,dimension(:,:)  :: fid1d_size,fid2d_size !tmp to save size of avg%v

      integer,dimension(maxnumbervars):: varssh1    ! temporal array
      var0d = 0
      var1d = 0
      var2d = 0
      numfid0d = 0
      numfid1d = 0
      numfid2d = 0
      write(*,*) 'inizialization of netCDF output: start'

      ALLOCATE( x(elem2D), y(elem2D), depth_elem(elem2D), &
                   &      STAT=output_alloc )
      if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays'

      call define_varlist

      ! get coordinates of elements
      do i = 1, elem2D
         elnodes = elem2D_nodes(:,i) !! 4 nodes in element
         numnodes = 4
         if(elnodes(1)==elnodes(4)) numnodes = 3  !! set to 3 if we have triangle
         y(i) = w_cv(1,i)*coord_nod2D(2,elnodes(1)) &
                + w_cv(2,i)*coord_nod2D(2,elnodes(2)) &
                + w_cv(3,i)*coord_nod2D(2,elnodes(3)) &
                + w_cv(4,i)*coord_nod2D(2,elnodes(4))
         x(i) = w_cv(1,i)*coord_nod2D(1,elnodes(1)) &
                + w_cv(2,i)*coord_nod2D(1,elnodes(2)) &
                + w_cv(3,i)*coord_nod2D(1,elnodes(3)) &
                + w_cv(4,i)*coord_nod2D(1,elnodes(4))
         !sum(coord_nod2d(1,elnodes(1:numnodes)))/dble(numnodes) /rad
         depth_elem(i) = w_cv(1,i)*depth(elem2D_nodes(1,i)) &
                       + w_cv(2,i)*depth(elem2D_nodes(2,i)) &
                       + w_cv(3,i)*depth(elem2D_nodes(3,i)) &
                       + w_cv(4,i)*depth(elem2D_nodes(4,i))
      end do
      filelist_unit= 15
      ! define files for output
      ! Read structure from outputfilelist.dat
      open( unit=filelist_unit, FILE='outputfilelist.dat', status='old', iostat=iost )
      if( iost == 0 ) then
         WRITE(*,*) '     file   : outputfilelist.dat open ok'
      else
         WRITE(*,*) 'ERROR: --> bad opening file   outputfilelist.dat ; iostat=',iost
         STOP 'ERROR: --> output_ini'
      endif
      ! reading number of files for output
      READ( unit=filelist_unit, FMT=*,iostat=iost)
      READ( unit=filelist_unit, FMT=*,iostat=iost) numfiles
      ALLOCATE( filelist(numfiles), &
                     &      STAT=output_alloc )
         if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays'
      ! read files options from outputfilelist.dat
      write(*,*) "start: read/check variable and file lists for output"
      write(*,*) "files for output: "
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
            write(*,*) "          more then maximum allowed (",maxnumbervars,"). hint: change maxnumbervars in fv_ncoutput.f90"
            STOP
         end if
         READ( unit=filelist_unit, FMT=*,iostat=iost)
         filelist(i)%vars=0
         filelist(i)%varsid=0
         READ( unit=filelist_unit, FMT=*,iostat=iost) filelist(i)%vars

         do j =1, filelist(i)%numvars
            if ( filelist(i)%vars(j) < 1 ) then
               write(*,*) 'ERROR: inconsist number of variables with provided index,',trim(filelist(i)%name)
               stop
            end if
         end do
         j=0
         do while ( j < filelist(i)%numvars )
            j=j+1
            if ( varlist(filelist(i)%vars(j))%type_task > type_task ) then
               write(*,*) 'WARNING: variable "',trim(varlist(filelist(i)%vars(j))%name), &
                                             & '" removed from file:',trim(filelist(i)%name)
               write(*,*) '            while type_task inconsist with curent model type_task'
               filelist(i)%numvars = filelist(i)%numvars - 1
               varssh1 = cshift(filelist(i)%vars,1)
               filelist(i)%vars(j:maxnumbervars-1) = varssh1(j:maxnumbervars-1)
               j=j-1
            end if
         end do
         write(*,*) "          ",trim(filelist(i)%name)
      end do
      write(*,*) "end: read/check variable and file lists for output"
   !==================== define indexes of nodes and elements where to write
      do i = 1, numfiles
         if ( filelist(i)%ftype==0 ) then
         ! full setup
            filelist(i)%nnods  = nod2D
            filelist(i)%nelems = elem2D
         elseif ( filelist(i)%ftype==3 ) then
         ! station
            filelist(i)%nnods  = 1
            filelist(i)%nelems = 1
         else
            STOP 'output_ini: wrong/not supported output file type'
         endif

         ALLOCATE(filelist(i)%nods(filelist(i)%nnods) , STAT=output_alloc )
         if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays'
         ALLOCATE(filelist(i)%elems(filelist(i)%nelems) , STAT=output_alloc )
         if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays'

         if ( filelist(i)%ftype==0 ) then
         ! full setup
            do n = 1, nod2D
               filelist(i)%nods(n) = n
            enddo
            do n = 1, elem2D
               filelist(i)%elems(n) = n
            enddo
         elseif ( filelist(i)%ftype==3 ) then
         ! station
            call ne_node(filelist(i)%posx_ll*rad, filelist(i)%posy_ll*rad, filelist(i)%nods(1))
            call ne_element(filelist(i)%posx_ll*rad, filelist(i)%posy_ll*rad, filelist(i)%elems(1))
         endif
      enddo
    !=end================== define indexes of nodes and elements where to write

      !create files
      do i = 1, numfiles
         !create new file
         iost = nf_create(filelist(i)%name,NF_CLOBBER,filelist(i)%fileid)
         call check_nferr(iost,filelist(i)%name)
         ! define dimensions
         iost = nf_def_dim(filelist(i)%fileid, 'node', filelist(i)%nnods, filelist(i)%dimid%node)
         call check_nferr(iost,filelist(i)%name)
         iost = nf_def_dim(filelist(i)%fileid, 'nele', filelist(i)%nelems, filelist(i)%dimid%nelem)
         call check_nferr(iost,filelist(i)%name)
         iost = nf_def_dim(filelist(i)%fileid, 'maxnod', 4, filelist(i)%dimid%maxnod)
         call check_nferr(iost,filelist(i)%name)
         iost = nf_def_dim(filelist(i)%fileid, 'maxnodelem', maxval(nod_in_elem2D_num), filelist(i)%dimid%maxnodelem)
         call check_nferr(iost,filelist(i)%name)
         iost = nf_def_dim(filelist(i)%fileid, 'time', NF_UNLIMITED, filelist(i)%dimid%time)
         call check_nferr(iost,filelist(i)%name)
         if (type_task>1) then
            iost = nf_def_dim(filelist(i)%fileid, 'sigma', nsigma, filelist(i)%dimid%sigma)
            call check_nferr(iost,filelist(i)%name)
            iost = nf_def_dim(filelist(i)%fileid, 'sigmam1', nsigma-1, filelist(i)%dimid%sigmam1)
            call check_nferr(iost,filelist(i)%name)
         endif
         !define sigma var
         if (type_task>1) then
            nf_dims(1)=filelist(i)%dimid%sigma
!            nf_dims(2)=filelist(i)%dimid%node
            iost = nf_def_var(filelist(i)%fileid, 'sigma_lev', NF_DOUBLE,1,nf_dims,filelist(i)%varid%sigma)
            call check_nferr(iost,filelist(i)%name)
         endif
         !define nv
         nf_dims(1) = filelist(i)%dimid%maxnod
         nf_dims(2) = filelist(i)%dimid%nelem
         iost = nf_def_var(filelist(i)%fileid,'nv', NF_INT, 2,nf_dims,filelist(i)%varid%nv)
         call check_nferr(iost,filelist(i)%name)
         ! add attributes to face_node_connectivity
         iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%nv,'standard_name', &
                                  & len_trim('face_node_connectivity'),'face_node_connectivity')
         call check_nferr(iost,filelist(i)%name)
         n=1
         ! iost = nf_put_att(filelist(i)%fileid,id_m_nv,'start_index',NF_INT,1,n)
         ! call check_nferr(iost,filelist(i)%name)
         iost = nf_def_var(filelist(i)%fileid,'fvcoastal_mesh',NF_INT,0,nf_dims,filelist(i)%varid%mesh)
         call check_nferr(iost,filelist(i)%name)
         ! add attributes to mesh
         n=2
!IK it is needed, should be switch on after check         !iost = nf_put_att(filelist(i)%fileid,filelist(i)%varid%mesh,'topology_dimension',NF_INT,1,n)
         call check_nferr(iost,filelist(i)%name)
         iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%mesh,'cf_role',len_trim('mesh_topology'),'mesh_topology')
         call check_nferr(iost,filelist(i)%name)
         iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%mesh,'node_coordinates',len_trim('lon lat'),'lon lat')
         call check_nferr(iost,filelist(i)%name)
         iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%mesh,'face_node_connectivity',len_trim('nv'),'nv')
         call check_nferr(iost,filelist(i)%name)
         !  define coordinates on nodes
         nf_dims(1) = filelist(i)%dimid%node
         iost = nf_def_var(filelist(i)%fileid,'lon',NF_DOUBLE,1,nf_dims,filelist(i)%varid%lon)
         call check_nferr(iost,filelist(i)%name)
         iost = nf_def_var(filelist(i)%fileid,'lat',NF_DOUBLE,1,nf_dims,filelist(i)%varid%lat)
         call check_nferr(iost,filelist(i)%name)
         ! define depth
         iost = nf_def_var(filelist(i)%fileid,'depth',NF_DOUBLE,1,nf_dims,filelist(i)%varid%depth)
         call check_nferr(iost,filelist(i)%name)

         nf_dims(1) = filelist(i)%dimid%nelem
         iost = nf_def_var(filelist(i)%fileid,'depth_elem',NF_DOUBLE,1,nf_dims,filelist(i)%varid%depth_elem)
         call check_nferr(iost,filelist(i)%name)

         ! add attributes
         if (cartesian) then
             iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%lon,'standard_name',len_trim('X'),'X')
             call check_nferr(iost,filelist(i)%name)
             iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%lon,'units',len_trim('meters'),'meters')
             call check_nferr(iost,filelist(i)%name)
             iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%lat,'standard_name',len_trim('Y'),'Y')
             call check_nferr(iost,filelist(i)%name)
             iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%lat,'units',len_trim('meters'),'meters')
             call check_nferr(iost,filelist(i)%name)
         else
         iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%lon,'standard_name',len_trim('longitude'),'longitude')
         call check_nferr(iost,filelist(i)%name)
         iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%lon,'units',len_trim('degrees_east'),'degrees_east')
         call check_nferr(iost,filelist(i)%name)
         iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%lat,'standard_name',len_trim('latitude'),'latitude')
         call check_nferr(iost,filelist(i)%name)
         iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%lat,'units',len_trim('degrees_east'),'degrees_north')
         call check_nferr(iost,filelist(i)%name)
         end if
         !depth
         iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%depth,'standard_name',len_trim('depth'),'depth')
         call check_nferr(iost,filelist(i)%name)
         iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%depth,'units',len_trim('m'),'m')
         call check_nferr(iost,filelist(i)%name)
         iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%depth_elem,'standard_name',&
                                                                   len_trim('depth on elem'),'depth on elem')
         call check_nferr(iost,filelist(i)%name)
         iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%depth_elem,'units',len_trim('m'),'m')
         call check_nferr(iost,filelist(i)%name)
         ! define coordinates on elements
         nf_dims(1) = filelist(i)%dimid%nelem
         iost = nf_def_var(filelist(i)%fileid,'lon_elem',NF_DOUBLE,1,nf_dims,filelist(i)%varid%lon_elem)
         call check_nferr(iost,filelist(i)%name)
         iost = nf_def_var(filelist(i)%fileid,'lat_elem',NF_DOUBLE,1,nf_dims,filelist(i)%varid%lat_elem)
         call check_nferr(iost,filelist(i)%name)
         ! add attributes

         if (cartesian) then
         iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%lon_elem,'standard_name', &
                                      & len_trim('X_elements'),'X_elements')
             call check_nferr(iost,filelist(i)%name)
             iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%lon_elem,'units', &
                                      & len_trim('meters'),'meters')
             call check_nferr(iost,filelist(i)%name)
             iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%lat_elem,'standard_name', &
                                      & len_trim('Y_elements'),'Y_elements')
             call check_nferr(iost,filelist(i)%name)
             iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%lat_elem,'units', &
                                      & len_trim('meters'),'meters')
             call check_nferr(iost,filelist(i)%name)
         else
             iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%lon_elem,'standard_name', &
                                      & len_trim('longitude_elements'),'longitude_elements')
         call check_nferr(iost,filelist(i)%name)
         iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%lon_elem,'units', &
                                  & len_trim('degrees_east'),'degrees_east')
         call check_nferr(iost,filelist(i)%name)
         iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%lat_elem,'standard_name', &
                                  & len_trim('latitude_elements'),'latitude_elements')
         call check_nferr(iost,filelist(i)%name)
         iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%lat_elem,'units', &
                                  & len_trim('degrees_east'),'degrees_north')
         call check_nferr(iost,filelist(i)%name)
         end if
         ! define time
         nf_dims(1)=filelist(i)%dimid%time
         iost = nf_def_var(filelist(i)%fileid,'time',NF_DOUBLE,1,nf_dims,filelist(i)%varid%time)
         call check_nferr(iost,filelist(i)%name)
         iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%time,'standard_name',len_trim('time'),'time')
         call check_nferr(iost,filelist(i)%name)
         !  define elements area, nods area and weights
         nf_dims(1) = filelist(i)%dimid%node
         iost = nf_def_var(filelist(i)%fileid,'area',NF_DOUBLE,1,nf_dims,filelist(i)%varid%area)
         call check_nferr(iost,filelist(i)%name)

         nf_dims(1) = filelist(i)%dimid%nelem
         iost = nf_def_var(filelist(i)%fileid,'elem_area',NF_DOUBLE,1,nf_dims,filelist(i)%varid%elem_area)
         call check_nferr(iost,filelist(i)%name)
         nf_dims(1) = filelist(i)%dimid%maxnod
         nf_dims(2) = filelist(i)%dimid%nelem
         iost = nf_def_var(filelist(i)%fileid,'w_cv', NF_DOUBLE, 2,nf_dims,filelist(i)%varid%w_cv)
         call check_nferr(iost,filelist(i)%name)
         ! define nod_in_elem2D and nod_in_elem2D_num
         nf_dims(1) = filelist(i)%dimid%node
         iost = nf_def_var(filelist(i)%fileid,'nod_in_elem2d_num', NF_INT, 1,nf_dims,filelist(i)%varid%nod_in_elem2d_num)

         call check_nferr(iost,filelist(i)%name)
         nf_dims(1) = filelist(i)%dimid%maxnodelem
         nf_dims(2) = filelist(i)%dimid%node
         iost = nf_def_var(filelist(i)%fileid,'nod_in_elem2d', NF_INT, 2,nf_dims,filelist(i)%varid%nod_in_elem2d)
         call check_nferr(iost,filelist(i)%name)

         ! add attributes
         iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%area,&
                                'standard_name',len_trim('cell area'),'cell area')
         call check_nferr(iost,filelist(i)%name)
         iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%area,&
                                'units',len_trim('m**2'),'m**2')
         call check_nferr(iost,filelist(i)%name)
         iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%elem_area,&
                                'standard_name',len_trim('elemens area'),'elemens area')
         call check_nferr(iost,filelist(i)%name)
         iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%elem_area,&
                                'units',len_trim('m**2'),'m**2')
         call check_nferr(iost,filelist(i)%name)
         iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%w_cv,&
                                'standard_name',len_trim(' Scv/Sc'),' Scv/Sc')
         call check_nferr(iost,filelist(i)%name)
         iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%nod_in_elem2D_num,&
                                'standard_name',len_trim('node is in N elemets'),'node is in N elemets')
         call check_nferr(iost,filelist(i)%name)
         iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%nod_in_elem2D,&
                                'standard_name',len_trim('node is in elemets'),'node is in elemets')
         call check_nferr(iost,filelist(i)%name)
!         call calendar_date(floor(time_jd0),yyyy,mm,dd)
!         ihh=floor((time_jd0-floor(time_jd0))*86400.0/60/60)
!         imm=floor((86400.0*(time_jd0-floor(time_jd0))-ihh*3600.0)/60)
!         isec=floor(86400.0*(time_jd0-floor(time_jd0))-ihh*3600.0-imm*60.0)
!         write(yystr,"(I4.4)") yyyy
!         write(mmstr,"(I2.2)") mm
!         write(ddstr,"(I2.2)") dd
!         write(hr,"(I2.2)") ihh
!         write(mn,"(I2.2)") imm
!         write(sc,"(I2.2)") isec
!         str='days since ' // yystr // '-' // mmstr // '-' // ddstr &
!                   &// ' ' // hr    // ':' // mn    // ':' // sc
         str='days since 1900-01-01 00:00:00'
         iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%time,'units', &
                        &   len_trim(str),str)
         call check_nferr(iost,filelist(i)%name)
         iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%time,'axis',len_trim('T'),'T')
         call check_nferr(iost,filelist(i)%name)
         iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%time,'title',len_trim('Time'),'Time')
         call check_nferr(iost,filelist(i)%name)
!      iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varid%time,'calendar',len_trim('julian'),'julian')
!      call check_nferr(iost,filelist(i)%name)

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

            iost = nf_def_var(filelist(i)%fileid,varlist(varid)%name,&
                                  NF_DOUBLE,ndim,nf_dims,filelist(i)%varsid(j))
            call check_nferr(iost,filelist(i)%name)
            iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varsid(j),&
                   'units',len_trim(varlist(varid)%units),varlist(varid)%units)
            call check_nferr(iost,filelist(i)%name)
            iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varsid(j),&
                   'standard_name',len_trim(varlist(varid)%standard_name),varlist(varid)%standard_name)
            call check_nferr(iost,filelist(i)%name)
            iost = nf_put_att_text(filelist(i)%fileid,filelist(i)%varsid(j),&
                   'modelname',len_trim(varlist(varid)%modelname),varlist(varid)%modelname)
            call check_nferr(iost,filelist(i)%name)


         end do

         ! leaving define mode
         iost = nf_enddef(filelist(i)%fileid)
         call check_nferr(iost,filelist(i)%name)

         if (type_task>1) then
            ! sigma
            nf_start(1)=1
            nf_edges(1)=nsigma
            iost = nf_put_vara_double(filelist(i)%fileid,filelist(i)%varid%sigma,nf_start,nf_edges,sigma)
            call check_nferr(iost,filelist(i)%name)
         end if

         ! write down coordinates (nodes)
         nf_start(1)=1
         nf_edges(1)=filelist(i)%nnods
         if (cartesian) then
            iost = nf_put_vara_double(filelist(i)%fileid,filelist(i)%varid%lon,nf_start,nf_edges,&
                                          coord_nod2d(1,filelist(i)%nods(1:filelist(i)%nnods))*r_earth)
            call check_nferr(iost,filelist(i)%name)
            iost = nf_put_vara_double(filelist(i)%fileid,filelist(i)%varid%lat,nf_start,nf_edges,&
                                          coord_nod2d(2,filelist(i)%nods(1:filelist(i)%nnods))*r_earth)
            call check_nferr(iost,filelist(i)%name)
         else
         iost = nf_put_vara_double(filelist(i)%fileid,filelist(i)%varid%lon,nf_start,nf_edges,&
                                          coord_nod2d(1,filelist(i)%nods(1:filelist(i)%nnods))/rad )
         call check_nferr(iost,filelist(i)%name)
         iost = nf_put_vara_double(filelist(i)%fileid,filelist(i)%varid%lat,nf_start,nf_edges,&
                                          coord_nod2d(2,filelist(i)%nods(1:filelist(i)%nnods))/rad)
         call check_nferr(iost,filelist(i)%name)
         end if
         ! write depth
         iost = nf_put_vara_double(filelist(i)%fileid,filelist(i)%varid%depth,nf_start,nf_edges,&
                                          depth(filelist(i)%nods(1:filelist(i)%nnods)))
         call check_nferr(iost,filelist(i)%name)
         nf_start(1)=1
         nf_edges(1)=filelist(i)%nelems
         iost = nf_put_vara_double(filelist(i)%fileid,filelist(i)%varid%depth_elem,nf_start,nf_edges,&
                                          depth_elem(filelist(i)%elems(1:filelist(i)%nelems)))
         call check_nferr(iost,filelist(i)%name)
         ! write down coordinates (elements)
         nf_start(1)=1
         nf_edges(1)=filelist(i)%nelems
         if (cartesian) then
         iost = nf_put_vara_double(filelist(i)%fileid,filelist(i)%varid%lon_elem,nf_start,nf_edges,&
                                         x(filelist(i)%elems(1:filelist(i)%nelems))*r_earth)
         call check_nferr(iost,filelist(i)%name)
         iost = nf_put_vara_double(filelist(i)%fileid,filelist(i)%varid%lat_elem,nf_start,nf_edges,&
                                         y(filelist(i)%elems(1:filelist(i)%nelems))*r_earth)
         call check_nferr(iost,filelist(i)%name)
         else
            iost = nf_put_vara_double(filelist(i)%fileid,filelist(i)%varid%lon_elem,nf_start,nf_edges,&
                                         x(filelist(i)%elems(1:filelist(i)%nelems))/rad)
         call check_nferr(iost,filelist(i)%name)
            iost = nf_put_vara_double(filelist(i)%fileid,filelist(i)%varid%lat_elem,nf_start,nf_edges,&
                                         y(filelist(i)%elems(1:filelist(i)%nelems))/rad)
            call check_nferr(iost,filelist(i)%name)
         end if
         ! write down NV (
         !define nv elem2D_nodes
         nf_start(2)=1
         nf_edges(2)=filelist(i)%nelems
         nf_start(1)=1
         nf_edges(1)=4
         iost = nf_put_vara_int(filelist(i)%fileid,filelist(i)%varid%nv,nf_start,nf_edges,&
                                        elem2D_nodes(1:4,filelist(i)%elems(1:filelist(i)%nelems)))
         call check_nferr(iost,filelist(i)%name)

         !  write elements area, nods area and weights
         nf_start(1)=1
         nf_edges(1)=filelist(i)%nnods
         iost = nf_put_vara_double(filelist(i)%fileid,filelist(i)%varid%area,nf_start,nf_edges,&
                                                         area(filelist(i)%nods(1:filelist(i)%nnods)))
         call check_nferr(iost,filelist(i)%name)
         nf_start(1)=1
         nf_edges(1)=filelist(i)%nelems
         iost = nf_put_vara_double(filelist(i)%fileid,filelist(i)%varid%elem_area,nf_start,nf_edges,&
                                        elem_area(filelist(i)%elems(1:filelist(i)%nelems)))
         call check_nferr(iost,filelist(i)%name)
         nf_start(2)=1
         nf_edges(2)=filelist(i)%nelems
         nf_start(1)=1
         nf_edges(1)=4
         iost = nf_put_vara_double(filelist(i)%fileid,filelist(i)%varid%w_cv,nf_start,nf_edges,&
                                        w_cv(1:4,filelist(i)%elems(1:filelist(i)%nelems)))
         call check_nferr(iost,filelist(i)%name)
         nf_start(1)=1
         nf_edges(1)=filelist(i)%nnods
         iost = nf_put_vara_int(filelist(i)%fileid,filelist(i)%varid%nod_in_elem2D_num,nf_start,nf_edges,&
                                    nod_in_elem2D_num(filelist(i)%nods(1:filelist(i)%nnods)))
         call check_nferr(iost,filelist(i)%name)
         nf_start(2)=1
         nf_edges(2)=filelist(i)%nnods
         nf_start(1)=1
         nf_edges(1)=maxval(nod_in_elem2D_num)
         iost = nf_put_vara_int(filelist(i)%fileid,filelist(i)%varid%nod_in_elem2D,nf_start,nf_edges,&
                                    nod_in_elem2D(1:nf_edges(1),filelist(i)%nods(1:filelist(i)%nnods)))
         call check_nferr(iost,filelist(i)%name)

         ! close file
         iost = nf_close(filelist(i)%fileid)
         call check_nferr(iost,filelist(i)%name)
         ! define 0 records (time) in nc outfile
         filelist(i)%nc_rec = 0
      end do

      DEALLOCATE( x, y)

      ALLOCATE( fid1d_size(maxnumberallvars,numfiles), fid2d_size(maxnumberallvars,numfiles), &
                   &      STAT=output_alloc )
      if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays'

      ! count number of variables of each type
      nvar0d = 0
      nvar1d = 0
      nvar2d = 0
      do i = 1, numfiles
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
      write(*,*) "output_ini:"
      write(*,*) "            Number of 0D variables for output:", nvar0d
      write(*,*) "            Number of 1D variables for output:", nvar1d
      write(*,*) "            Number of 2D variables for output:", nvar2d

      ! allocate memory for list of variables for output
      ! 3 lists for each dimension (time is a dimension and x-y plane is one dimension)
      if (nvar0d > 0) then
         ALLOCATE( varbas0d(nvar0d), &
                   &      STAT=output_alloc )
         if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays'
         do i = 1, nvar0d
            varbas0d(i)%varid = var0d(i)
            ALLOCATE( varbas0d(i)%avg(numfid0d(i)), &
                      &      STAT=output_alloc )
            if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays'
            do j = 1, numfid0d(i)
               varbas0d(i)%avg(j)%v = 0
            enddo
            varbas0d(i)%nfid = numfid0d(i)
         end do
      end if

      if (nvar1d > 0) then
         ALLOCATE( varbas1d(nvar1d), &
                   &      STAT=output_alloc )
         if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays'
         do i = 1, nvar1d
            !
            varbas1d(i)%varid = var1d(i)
            !allocate memory for snapshot of every variable and
            ALLOCATE(varbas1d(i)%snapshot(varlist(var1d(i))%varsize(1)), &
                      & varbas1d(i)%avg(numfid1d(i) ), &
                      &      STAT=output_alloc )
            if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays'

            do j = 1, numfid1d(i)
               ALLOCATE( varbas1d(i)%avg(j)%v(fid1d_size(i,j)), &
                         varbas1d(i)%avg(j)%sn(fid1d_size(i,j)), &
                      &      STAT=output_alloc )
               if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays'
               varbas1d(i)%avg(j)%sz = fid1d_size(i,j)
               varbas1d(i)%avg(j)%v = 0
            enddo

            varbas1d(i)%nfid = numfid1d(i)
         end do
      end if

      if (nvar2d > 0) then
         ALLOCATE( varbas2d(nvar2d), &
                   &      STAT=output_alloc )
         if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays'
         do i = 1, nvar2d
            varbas2d(i)%varid = var2d(i)
            ALLOCATE( varbas2d(i)%snapshot(varlist(var2d(i))%varsize(1), varlist(var2d(i))%varsize(2) ),  &
                      varbas2d(i)%avg(numfid2d(i)), &
                        &      STAT=output_alloc )
            if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays'
            do j = 1, numfid2d(i)
               ALLOCATE( varbas2d(i)%avg(j)%v(varlist(var2d(i))%varsize(1),   &
                               fid2d_size(i,j) ),  &
                         varbas2d(i)%avg(j)%sn(varlist(var2d(i))%varsize(1),   &
                               fid2d_size(i,j) ),  &
                        &      STAT=output_alloc )
               if( output_alloc /= 0 )   STOP 'output_ini: failed to allocate arrays'
               varbas2d(i)%avg(j)%sz(1) = varlist(var2d(i))%varsize(1)
               varbas2d(i)%avg(j)%sz(2) = fid2d_size(i,j)
               varbas2d(i)%avg(j)%v = 0
            enddo
            varbas2d(i)%nfid = numfid2d(i)
         end do
      end if

      DEALLOCATE( fid1d_size, fid2d_size)

   END SUBROUTINE output_ini

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
            varsize(ndim) = elem2D
         case (1)
            ndim = ndim+1
            varsize(ndim) = nod2D
      end select

   END SUBROUTINE getvarsize


   SUBROUTINE write_data(findx)
      IMPLICIT NONE

      integer, intent(in)  :: findx  ! index of file to write in filelist
      integer              :: i
      real(kind=WP)        :: tshift
      !open for writing netCDF file
      iost = nf_open(filelist(findx)%name,NF_WRITE,filelist(findx)%fileid)
      call check_nferr(iost,filelist(findx)%name)

      !increase time record for ncout
      filelist(findx)%nc_rec = filelist(findx)%nc_rec + 1
      ! write time
      nf_start(1)=filelist(findx)%nc_rec
      nf_edges(1)=1
      ! time_jd - time_jd0 - is days since begining of simulations
      ! changed to 2415020.5 (1900-01-01 00:00)
      tshift = 0.0_WP
      if (filelist(findx)%snapshot == 2) then
        tshift = filelist(findx)%freq*dt*0.5_WP/86400.0_WP
      endif
      iost = nf_put_vara_double(filelist(findx)%fileid,&
              & filelist(findx)%varid%time,nf_start,nf_edges, time_jd - 2415020.5 - tshift) !time_jd0)
      call check_nferr(iost,filelist(findx)%name)
      ! write variables one by one
      do i = 1, filelist(findx)%numvars
         call write_var(findx,i)
      end do

      ! close file
      iost = nf_close(filelist(findx)%fileid)
      call check_nferr(iost,filelist(findx)%name)

   END SUBROUTINE write_data

   SUBROUTINE write_var(findx,vindx)
      IMPLICIT NONE

      integer, intent(in)  :: findx  ! index of file to write in filelist
      integer, intent(in)  :: vindx  ! index of variable to write
      integer              :: i, ndim, varid


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
            nf_start(ndim)=1
            nf_edges(ndim)=filelist(findx)%nelems
         case (1)
            ndim = ndim+1
            nf_start(ndim)=1
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
               iost = nf_put_vara_double(filelist(findx)%fileid,filelist(findx)%varsid(vindx),nf_start,nf_edges,&
                   & varbas0d(filelist(findx)%fid(vindx))%avg(filelist(findx)%avgid(vindx))%sn )
            case(2)
               iost = nf_put_vara_double(filelist(findx)%fileid,filelist(findx)%varsid(vindx),nf_start,nf_edges,&
                   & varbas1d(filelist(findx)%fid(vindx))%avg(filelist(findx)%avgid(vindx))%sn )
            case(3)
               iost = nf_put_vara_double(filelist(findx)%fileid,filelist(findx)%varsid(vindx),nf_start,nf_edges,&
                   & varbas2d(filelist(findx)%fid(vindx))%avg(filelist(findx)%avgid(vindx))%sn )
         end select
      else
        ! if we write averaged values
         select case (varlist(filelist(findx)%vars(vindx))%ndim) ! select proper dimension
            case(1)
               iost = nf_put_vara_double(filelist(findx)%fileid,filelist(findx)%varsid(vindx),nf_start,nf_edges,&
                   & varbas0d(filelist(findx)%fid(vindx))%avg(filelist(findx)%avgid(vindx))%v )
            case(2)
               iost = nf_put_vara_double(filelist(findx)%fileid,filelist(findx)%varsid(vindx),nf_start,nf_edges,&
                   & varbas1d(filelist(findx)%fid(vindx))%avg(filelist(findx)%avgid(vindx))%v )
            case(3)
               iost = nf_put_vara_double(filelist(findx)%fileid,filelist(findx)%varsid(vindx),nf_start,nf_edges,&
                   & varbas2d(filelist(findx)%fid(vindx))%avg(filelist(findx)%avgid(vindx))%v )
         end select
      end if
      call check_nferr(iost,filelist(findx)%name)

   END SUBROUTINE write_var



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
         do j = 1, filelist(i)%numvars
            fid   = filelist(i)%fid(j)
            avgid = filelist(i)%avgid(j)
            vid   = filelist(i)%vars(j)
            select case (varlist(vid)%ndim) ! select proper dimension
               case(1)
                     varbas0d(fid)%avg(avgid)%sn = varbas0d(fid)%snapshot
               case(2)
                  if (varlist(vid)%onnode == 1 ) then
                     varbas1d(fid)%avg(avgid)%sn = varbas1d(fid)%snapshot(filelist(i)%nods(1:filelist(i)%nnods))
                  else
                     varbas1d(fid)%avg(avgid)%sn = varbas1d(fid)%snapshot(filelist(i)%elems(1:filelist(i)%nelems))
                  end if
               case(3)
                  if (varlist(vid)%onnode == 1 ) then
                     varbas2d(fid)%avg(avgid)%sn = varbas2d(fid)%snapshot(:,filelist(i)%nods(1:filelist(i)%nnods))
                  else
                     varbas2d(fid)%avg(avgid)%sn = varbas2d(fid)%snapshot(:,filelist(i)%elems(1:filelist(i)%nelems))
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



   SUBROUTINE write_mesh2nc
   !!!! NOT in USE anymore
      IMPLICIT NONE

      integer      :: ncid ! netcdf status
      ! ID dimensions and variables:
      integer      :: id_node
      integer      :: id_lon
      integer      :: id_lat
      integer      :: id_nele
      integer      :: id_maxnod
      integer      :: id_maxnodelem
      integer      :: id_nv
      integer      :: id_mesh
      ! temporal arrays
      integer      :: nf_dims(4) ! dimensions (temporal)
      integer      :: nf_start(4)
      integer      :: nf_edges(4)
      !   real(wp)  :: elnodes(4), x, y
      integer      :: n, elem
      character(len=256)  :: file_name


      file_name='mesh.nc'

      write(*,*) 'write mesh to ',trim(file_name),' start'

      !create new file
      iost = nf_create(file_name,NF_CLOBBER,ncid)
      call check_nferr(iost,file_name)

      ! define dimensions
      iost = nf_def_dim(ncid, 'node', nod2D, id_node)
      call check_nferr(iost,file_name)
      iost = nf_def_dim(ncid, 'nele', elem2D, id_nele)
      call check_nferr(iost,file_name)
      iost = nf_def_dim(ncid, 'maxnod', 4, id_maxnod)
      call check_nferr(iost,file_name)
      iost = nf_def_dim(ncid, 'maxnodelem', maxval(nod_in_elem2D_num), id_maxnodelem)
      call check_nferr(iost,file_name)

      !define nv
      nf_dims(1) = id_maxnod
      nf_dims(2) = id_nele
      iost = nf_def_var(ncid,'nv', NF_INT, 2,nf_dims,id_nv)
      call check_nferr(iost,file_name)
      ! add attributes to face_node_connectivity
      iost = nf_put_att_text(ncid,id_nv,'standard_name',len_trim('face_node_connectivity'),'face_node_connectivity')
      call check_nferr(iost,file_name)
      n=1
   !   iost = nf_put_att(ncid,id_nv,'start_index',NF_INT,1,n)
      call check_nferr(iost,file_name)
      iost = nf_def_var(ncid,'fvcoastal_mesh',NF_INT,0,nf_dims,id_mesh)
      call check_nferr(iost,file_name)

      ! add attributes to mesh
      n=2
   !   iost = nf_put_att(ncid,id_mesh,'topology_dimension',NF_INT,1,n)
      call check_nferr(iost,file_name)
      iost = nf_put_att_text(ncid,id_mesh,'cf_role',len_trim('mesh_topology'),'mesh_topology')
      call check_nferr(iost,file_name)
      iost = nf_put_att_text(ncid,id_mesh,'node_coordinates',len_trim('lon lat'),'lon lat')
      call check_nferr(iost,file_name)
      iost = nf_put_att_text(ncid,id_mesh,'face_node_connectivity',len_trim('nv'),'nv')
      call check_nferr(iost,file_name)

      !  define coordinates
      nf_dims(1) = id_node
      iost = nf_def_var(ncid,'lon',NF_DOUBLE,1,nf_dims,id_lon)
      call check_nferr(iost,file_name)
      iost = nf_def_var(ncid,'lat',NF_DOUBLE,1,nf_dims,id_lat)
      call check_nferr(iost,file_name)
      ! add attributes
      iost = nf_put_att_text(ncid,id_lon,'standard_name',len_trim('longitude'),'longitude')
      call check_nferr(iost,file_name)
      iost = nf_put_att_text(ncid,id_lon,'units',len_trim('degrees_east'),'degrees_east')
      call check_nferr(iost,file_name)
      iost = nf_put_att_text(ncid,id_lat,'standard_name',len_trim('latitude'),'latitude')
      call check_nferr(iost,file_name)
      iost = nf_put_att_text(ncid,id_lat,'units',len_trim('degrees_east'),'degrees_east')
      call check_nferr(iost,file_name)

      ! leaving define mode
      iost = nf_enddef(ncid)
      call check_nferr(iost,file_name)
 !     write(*,*) coord_nod2d(1,1:10)/rad
      nf_start(1)=1
      nf_edges(1)=nod2D
      iost = nf_put_vara_double(ncid,id_lon,nf_start,nf_edges,coord_nod2d(1,:)/rad)
!      write(*,*) iost, trim(nf_strerror(iost))
      call check_nferr(iost,file_name)
      iost = nf_put_vara_double(ncid,id_lat,nf_start,nf_edges,coord_nod2d(2,:)/rad)
      call check_nferr(iost,file_name)
      !close netcdf file

      iost = nf_close(ncid)
      call check_nferr(iost,file_name)
      write(*,*) 'write mesh to ',len_trim(file_name),' done'
   !!!! NOT in USE anymore
   END SUBROUTINE write_mesh2nc

   SUBROUTINE check_nferr(iost,fname)
      IMPLICIT NONE

      character(len=256), intent(in) :: fname
      integer, intent(in) :: iost

      if (iost .ne. NF_NOERR) then
         write(*,*) 'ERROR: I/O status= "',trim(nf_strerror(iost)),'";',iost,' file= ',fname
         STOP 'ERROR: stop'
      endif
   END SUBROUTINE

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
   SUBROUTINE ne_node(x,y,ind)
!===================================================================================|
! Looking for a closest node to x, y and return node index                             |
!   if the coordinates are given in Cartesian coord., then usual Eulerian norm is   |
!   used to calculate the distance.                                                 |
!   If the coordinates are given in spherical coord., than Haversine Formula is used|
!===================================================================================|
   IMPLICIT NONE
   real(kind=WP),intent(in) :: x, y
   integer,intent(out)      :: ind
   real(kind=WP)            :: mindist, d
   integer                  :: n

   mindist = 999999999.9
   if (cartesian) then
      do n = 1, nod2D
         d = (y/r_earth - coord_nod2D(2,n))**2 + (x/r_earth - coord_nod2D(1,n))**2
         d = sqrt(d)
         if ( mindist > d ) then
            ind = n
            mindist = d
         endif
      enddo
   else
      do n = 1, nod2D
         d = sin((y*rad - coord_nod2D(2,n))*0.5_WP)**2 + &
             cos(coord_nod2D(2,n)) * cos(y*rad) * sin((x*rad - coord_nod2D(1,n))*0.5_WP)**2
         d = asin(min(1.0_WP,sqrt(d)))
         if ( mindist > d ) then
            ind = n
            mindist = d
         endif
      enddo
   endif

   END SUBROUTINE ne_node

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

END MODULE fv_ncoutput

