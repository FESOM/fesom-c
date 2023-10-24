#ifndef NO_NETCDF
SUBROUTINE nc_init

!----- uses

  USE o_PARAM
  USE o_MESH
  USE o_ARRAYS

!-----local declarations


  IMPLICIT NONE

  include 'netcdf.inc'

  CHARACTER(len=100)          :: ncfile, file_name
  CHARACTER(len=200)          :: attstr
  INTEGER                      :: attValInt
  REAL(KIND=8)                 :: attVal, notyetknown
  CHARACTER(len=2)            :: par_num
  CHARACTER(len=3)            :: StrMeshID
  INTEGER                     :: i
  INTEGER                     :: ndummy
  INTEGER                     :: fileid 			     ! IDs netCDF file
  INTEGER                     :: DimId_iter			     ! dimension: iteration
  INTEGER                     :: DimId_n2D			     ! dimension: nodes
  INTEGER                     :: DimId_el2D			     ! dimension: elements
  INTEGER                     :: DimId_lyrs			     ! dimension: vertical layers
  INTEGER                     :: Dim1, Dim2, Dim3, Dim4	     ! dimensions: 
  INTEGER                     :: VarId_time, VarId_iter	     ! variables: time, iteration 
  INTEGER                     :: VarId_trian			     ! variable: triangle
  INTEGER                     :: VarId_loc2D, VarId_lon, VarId_lat ! variables: 2D location, longitude, latitude
  INTEGER                     :: VarId_ssh			     ! variable: sea surface height
  INTEGER                     :: VarId_temp			     ! variable: temperature
  INTEGER                     :: VarId_salt			     ! variable: salinity
!!$    INTEGER                     :: VarId_tri_area		     ! variable: area
!!$    INTEGER                     :: VarId_topo, VarId_index	     ! variables: topography, index
!!$    INTEGER                     :: VarId_arrival                     ! variable: arrival time
!!$    INTEGER                     :: VarId_ttt                         ! variable: estimated arrival time ttt
!!$    INTEGER                     :: VarId_mwh  			     ! variable: maximum wave height
!!$    INTEGER                     :: VarId_flux, VarId_vel             ! variables: flux and mean absolute velocity
!!$    INTEGER                     :: VarId_nodal_vel                   ! variable: velocity in nodes
!!$
  INTEGER                     :: s       			     ! auxiliary: status counter
  INTEGER, DIMENSION(3)       :: dimarray			     ! auxiliary: array dimension
  INTEGER, DIMENSION(500)     :: stat    			     ! auxiliary: status array
!!$        
  REAL, DIMENSION(2,nod2D)    :: AUX_loc
!!$    REAL(KIND=8)  ::   X0, Y0, D, L, W, TH, DL, HH, RD
!!$
!!$    REAL(KIND=8)  :: att_mw
!!$    integer ::   att_epi_i,  att_epi_j
!!$
!!$    ncfile = ncfile_all

  ncfile='./fesom_c_out.nc'

!----- open file and write global attributes

  s = 1
  stat(s) = NF_CREATE(TRIM(ncfile), 0, fileid) 
  s = s + 1
  attstr  = 'FESOM_C Simulation'
  stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'title',       LEN_TRIM(attstr),    &
                              TRIM(attstr)) 
    s = s + 1
!!$    attstr  = 'AWI'
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'institution', LEN_TRIM(attstr),    &
!!$                              TRIM(attstr)) 
!!$    s = s + 1
!!$    IF (enable_ruptgen_scenario) THEN
!!$       stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'TSID',     LEN_TRIM(TSID), TSID) 
!!$       s = s + 1
!!$    END IF
    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'Conventions', LEN_TRIM('CF-1.0'),  &
                              'CF-1.0') 

!!$
!!$    WRITE(StrMeshID, '(I3.3)') MeshID
!!$    attstr  = 'TsunAWI Revision '//revision//' - Mesh: '//StrMeshID
!!$    s = s + 1
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'source', LEN_TRIM(attstr),     &
!!$                              TRIM(attstr)) 
!!$    s = s + 1
!!$    attstr  = 'no entry'
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'history', LEN_TRIM(attstr),    &
!!$                              TRIM(attstr)) 
!!$    s = s + 1
!!$    attstr  = 'none'
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'references', LEN_TRIM(attstr), &
!!$                              TRIM(attstr)) 
!!$    s = s + 1
!!$    attstr  = 'none'
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'comment', LEN_TRIM(attstr),    &
!!$                              TRIM(attstr)) 
!!$    s = s + 1
!!$
!!$    attstr  = 'AWI'
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'provider', LEN_TRIM(attstr),    &
!!$                              TRIM(attstr)) 
!!$    s = s + 1
!!$
!!$    attstr  = 'TsunAWI'
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'model', LEN_TRIM(attstr),    &
!!$                              TRIM(attstr)) 
!!$    s = s + 1
!!$
!!$    attstr  = revision
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'modelVersion', LEN_TRIM(attstr),    &
!!$                              TRIM(attstr)) 
!!$    s = s + 1
!!$
!!$
!!$    attstr  = 'regional'
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'scope', LEN_TRIM(attstr),    &
!!$                              TRIM(attstr)) 
!!$    s = s + 1
!!$
!!$    notyetknown=-99999.
!!$
!!$    stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'globalMinSSH',         NF_DOUBLE, 1, notyetknown) 
!!$    s = s + 1
!!$
!!$    stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'globalMaxSSH',         NF_DOUBLE, 1, notyetknown) 
!!$    s = s + 1
!!$
!!$    attval=maxval(coord_nod2D(1,:))/rad
!!$    stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'boundingBox_UpperRight_Lon',         NF_DOUBLE, 1, attval) 
!!$    s = s + 1
!!$
!!$    attval=maxval(coord_nod2D(2,:))/rad
!!$    stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'boundingBox_UpperRight_Lat',         NF_DOUBLE, 1, attval) 
!!$    s = s + 1
!!$
!!$    attval=minval(coord_nod2D(1,:))/rad
!!$    stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'boundingBox_LowerLeft_Lon',         NF_DOUBLE, 1, attval) 
!!$    s = s + 1
!!$
!!$    attval=minval(coord_nod2D(2,:))/rad
!!$    stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'boundingBox_LowerLeft_Lat',         NF_DOUBLE, 1, attval) 
!!$
!!$    
!!$    s = s + 1
!!$    attval=T_out
!!$    stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'outputTimestepSize',         NF_DOUBLE, 1, attval) 
!!$    s = s + 1
!!$    
!!$    attval=T_end
!!$    stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'integrationTime',         NF_DOUBLE, 1, attval)
!!$    s = s + 1
!!$    
!!$    attval=longestEdge
!!$    stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'spatialResolution_Coarse',         NF_DOUBLE, 1, attval)
!!$    s = s + 1
!!$    
!!$    attval=shortestEdge
!!$    stat(s) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'spatialResolution_Fine',         NF_DOUBLE, 1, attval)
!!$    s = s + 1
!!$
!!$    attValInt=nod2D
!!$    stat(s) = NF_PUT_ATT_INT(fileid,    NF_GLOBAL,  'degreeOfFreedom',          NF_INT,    1, attValInt)
!!$
!!$    if (enable_benchmark .and. benchmark_ident<=2) then
!!$       s = s + 1
!!$       if  (benchmark_ident==1) attstr  = 'Bathymetry / topography data: Sloping beach experiment'
!!$       if  (benchmark_ident==1) attstr  = 'Bathymetry / topography data: Monai beach channel experiment'
!!$
!!$       stat(s) = NF_PUT_ATT_TEXT(fileid, NF_GLOBAL, 'meshDataDescription', LEN_TRIM(attstr),    &
!!$            TRIM(attstr)) 
!!$    endif
!!$    
!!$
!!$
!!$
!!$!----- DEFINE DIMENSIONS ------------------------------
!!$
    s = s + 1
    stat(s) = NF_DEF_DIM(fileid, 'nodes_2D',     nod2D, DimId_n2D)             
    s = s + 1
    stat(s) = NF_DEF_DIM(fileid, 'elements_2D', elem2D, DimId_el2D)        
    s = s + 1
    stat(s) = NF_DEF_DIM(fileid, 'vert_layers', nsigma-1 , DimId_lyrs)        
    s = s + 1
    stat(s) = NF_DEF_DIM(fileid, 'one',   1, dim1)                           
    s = s + 1
    stat(s) = NF_DEF_DIM(fileid, 'two',   2, dim2)                           
    s = s + 1
    stat(s) = NF_DEF_DIM(fileid, 'three', 3, dim3)                         
    s = s + 1
    stat(s) = NF_DEF_DIM(fileid, 'four',  4, dim4)                          
    s = s + 1
    stat(s) = NF_DEF_DIM(fileid, 'iteration', NF_UNLIMITED, DimId_iter)

    DO i = 1,  s
       IF (stat(i) /= NF_NOERR) &
          WRITE(*, *) 'NetCDF error',stat(i),' in dimension definitions, no.', i
    END DO


!----- DEFINE VARIABLES ---------------------------------

    s = 1

!----- coordinates

    dimarray(1) = dim2
    dimarray(2) = DimId_n2D
    
    stat(s) = NF_DEF_VAR(fileid, 'surface_locations', NF_FLOAT, 2, dimarray(1:2), VarId_loc2D) 

    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'longitude',         NF_FLOAT, 1, DimId_n2D,     VarId_lon) 

    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'latitude',          NF_FLOAT, 1, DimId_n2D,     VarId_lat) 


!----- grid

    dimarray(1) = dim4
    dimarray(2) = DimId_el2D

    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'triangles', NF_INT,   2, dimarray(1:2), VarId_trian) 

!!$    stat(s) = NF_DEF_VAR(fileid, 'tri_area',  NF_FLOAT, 1, DimId_el2D,    VarId_tri_area) 
!!$    s = s + 1
!!$
!!$!----- node indices
!!$    stat(s) = NF_DEF_VAR(fileid, 'node_index', NF_INT,  1, DimId_n2D,     VarId_index) 

!!$
!!$!----- numbers
!!$
    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'iteration', NF_INT,   1, DimId_iter,    VarId_iter) 

    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'time',   NF_DOUBLE,   1, DimId_iter,    VarId_time) 

!!$    s = s + 1
!!$
!!$
!!$!----- topography
!!$
!!$    stat(s) = NF_DEF_VAR(fileid, 'topography', NF_FLOAT, 1, DimId_n2D,    VarId_topo) 
!!$    s = s + 1

!----- F I E L D S 

!-----  scalar variables

    dimarray(1) = DimId_n2D
    dimarray(2) = DimId_iter

    stat(s) = NF_DEF_VAR(fileid, 'ssh', NF_FLOAT, 2, dimarray(1:2), VarId_ssh); 
    s = s + 1

    dimarray(1) = DimId_lyrs
    dimarray(2) = DimId_n2D
    dimarray(3) = DimId_iter

    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'temperature', NF_FLOAT, 3, dimarray(1:3), VarId_temp); 

    s = s + 1
    stat(s) = NF_DEF_VAR(fileid, 'salinity', NF_FLOAT, 3, dimarray(1:3), VarId_salt); 

!!$
!!$!----- arrival times, max. wave height, max. abs. velocity, max. flux, initial uplift
!!$
!!$    stat(s) = NF_DEF_VAR(fileid, 'first_arrival',     NF_FLOAT, 1, DimId_n2D, VarId_arrival)
!!$ 
!!$    if (write_ttt) then
!!$       s = s + 1
!!$       stat(s) = NF_DEF_VAR(fileid, 'first_arrival_ttt', NF_FLOAT, 1, DimId_n2D, VarId_ttt) 
!!$    endif
!!$    s = s + 1
!!$    stat(s) = NF_DEF_VAR(fileid, 'max_wave_height',   NF_FLOAT, 1, DimId_n2D, VarId_mwh) 
!!$    s = s + 1
!!$    stat(s) = NF_DEF_VAR (fileid, 'max_abs_vel',       NF_FLOAT, 1, DimId_n2D, VarId_vel) 
!!$    s = s + 1
!!$    stat(s) = NF_DEF_VAR(fileid, 'max_flux',          NF_FLOAT, 1, DimId_n2D, VarId_flux) 
!!$    s = s + 1
!!$    stat(s) = NF_DEF_VAR(fileid, 'initial_uplift',    NF_FLOAT, 1, DimId_n2D, VarId_uplift) 
!!$
!!$!----- vector variables
!!$
!!$!-----  velocity in nodes
!!$    IF (write_velocity) THEN
!!$       dimarray(1) = dim2
!!$       dimarray(2) = DimId_n2D
!!$       dimarray(3) = DimId_iter
!!$
!!$       s = s + 1
!!$       stat(s) = NF_DEF_VAR(fileid, 'velocity_in_nodes', NF_FLOAT, 3, dimarray(1:3), VarId_nodal_vel)
!!$    END IF
!!$
!!$ 
    DO i = 1,  s
        IF (stat(i) /= NF_NOERR)                                          &
            WRITE(*, *) 'NetCDF error',stat(i),' in variable definition, no.', i
    END DO

!----- DEFINE ATTRIBUTES -----------------------------------------

    s = 1

    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_loc2D,    'long_name',   19, 'surface_coordinates') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_loc2D,    'units',       18, 'degrees east/north') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_lon,      'long_name',    9, 'longitude') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_lon,      'units',       12, 'degrees_east') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_lat,      'long_name',    8, 'latitude') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_lat,      'units',       13, 'degrees_north') 
    s = s + 1
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_index,    'long_name',   16, 'nodal attributes') 
!!$    s = s + 1
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_tri_area, 'long_name',   13, 'triangle area') 
!!$    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ssh,      'long_name',   21, 'sea surface elevation') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ssh,      'units',        5, 'meter') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ssh,      'field',       19, 'ssh, scalar, series') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ssh,      'connections', 20, 'triangles, triangles') 
    s = s + 1
    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ssh,      'positions',   17, 'surface_locations') 
    s = s + 1
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_topo,     'long_name',   10, 'topography') 
!!$    s = s + 1
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_topo,     'field',       10, 'topography') 
!!$    s = s + 1
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_topo,     'connections', 20, 'triangles, triangles') 
!!$    s = s + 1
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_topo,     'positions',   17, 'surface_locations') 
!!$    s = s + 1
!!$    
!!$    attstr  = 'meters (after shock bathymetry and topography)'
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_topo,     'units',       LEN_TRIM(attstr), attstr) 
!!$    s = s + 1
!!$    
!!$    attstr='initial uplift due to earthquake'
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_uplift,   'long_name',   LEN_TRIM(attstr), attstr) 
!!$    s = s + 1
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_uplift,   'field',       14, 'initial_uplift') 
!!$    s = s + 1
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_uplift,   'connections', 20, 'triangles, triangles') 
!!$    s = s + 1
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_uplift,   'positions',   17, 'surface_locations') 
!!$    s = s + 1
!!$    
!!$    attstr  = 'meters'
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_time,     'units',       LEN_TRIM(attstr), attstr) 
!!$    s = s + 1
!!$    attstr  = 'first wave arrival (sea surface elevation larger than +0.01m)'
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_arrival,   'long_name',   LEN_TRIM(attstr), attstr) 
!!$    s = s + 1    
!!$    attstr = 'eta_thd_orig'
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_arrival,   'eta_type', LEN_TRIM(attstr),attstr)
!!$    s = s + 1
!!$    attstr  = 'seconds after rupture'
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_arrival,   'units',       LEN_TRIM(attstr), attstr) 
!!$    s = s + 1
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_arrival,   'field',       13, 'first_arrival') 
!!$    s = s + 1
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_arrival,   'connections', 20, 'triangles, triangles') 
!!$    s = s + 1
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_arrival,   'positions',   17, 'surface_locations') 
!!$
!!$    if (write_ttt) then
!!$       s = s + 1
!!$       attstr  = 'first wave arrival (estimated arrival time - ttt)'
!!$       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ttt      , 'long_name',   LEN_TRIM(attstr), attstr) 
!!$       s = s + 1
!!$       attstr = 'ttt_est_orig'
!!$       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ttt,   'eta_type', LEN_TRIM(attstr),attstr)
!!$       s = s + 1
!!$ 
!!$       attstr  = 'seconds after rupture'
!!$       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ttt,       'units',       LEN_TRIM(attstr), attstr) 
!!$       s = s + 1
!!$       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ttt,       'field',       17, 'first_arrival_ttt') 
!!$       s = s + 1
!!$       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ttt,       'connections', 20, 'triangles, triangles') 
!!$       s = s + 1
!!$       stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_ttt,         'positions',   17, 'surface_locations') 
!!$    endif
!!$
!!$    s = s + 1
!!$    attstr  = 'maximum wave height\n (above mean sea level on sea \n and above after shock topography on land)'
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_mwh,       'long_name',   LEN_TRIM(attstr), attstr) 
!!$    s = s + 1
!!$    
!!$    attstr  = 'meters'
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_mwh,       'units',       LEN_TRIM(attstr), attstr) 
!!$    s = s + 1
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_mwh,       'field',       15, 'max_wave_height') 
!!$    s = s + 1
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_mwh,       'connections', 20, 'triangles, triangles') 
!!$    s = s + 1
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_mwh,       'positions',   17, 'surface_locations') 
!!$    s = s + 1
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_vel,       'long_name',   36, 'maximum nodal velocity approximation') 
!!$    s = s + 1
!!$    
!!$    attstr  = 'm/sec'
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_vel,       'units',       LEN_TRIM(attstr), attstr) 
!!$    s = s + 1
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_vel,       'field',       11, 'max_abs_vel'); s=s+1
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_vel,       'connections', 20, 'triangles, triangles') 
!!$    s = s + 1
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_vel,       'positions',   17, 'surface_locations') 
!!$    s = s + 1
!!$
!!$    attstr  = 'maximum flux in nodes\n (maximum value of the product of velocity and wave height)'
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_flux,      'long_name',   LEN_TRIM(attstr), attstr) 
!!$    s = s + 1
!!$    
!!$    attstr  = 'm^2/sec'
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_flux,     'units',        LEN_TRIM(attstr), attstr) 
!!$    s = s + 1
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_flux,     'field',         8, 'max_flux') 
!!$    s = s + 1
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_flux,     'connections',  20, 'triangles, triangles') 
!!$    s = s + 1
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_flux,     'positions',    17, 'surface_locations') 
!!$    s = s + 1
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_time,     'long_name',     4, 'time'); 
!!$    s = s + 1
!!$    
!!$    attstr  = 'seconds since rupture'
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_time,     'units',        LEN_TRIM(attstr), attstr) 
!!$    s = s + 1
!!$    
!!$    attstr  = '0'
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_time,     'time_origin',  LEN_TRIM(attstr), attstr) 
!!$    s = s + 1
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_time,     'field',        12, 'time, series') 
!!$    s = s + 1
!!$    stat(s) = NF_PUT_ATT_TEXT(fileid, VarId_iter,     'long_name',     9, 'iteration') 
!!$
!!$    DO i = 1, s
!!$        IF (stat(i) /= NF_NOERR)                                          &
!!$            WRITE(*, *) 'NetCDF error',stat(i),' in attribute assignments, no.', i
!!$    END DO
!!$
    s = 1

    stat(s) = NF_ENDDEF(fileid) 
    

!----- write constant variables

!----- coordinates

!!$    IF (coordinate_type == 1) THEN
!!$        AUX_loc(1, :) = coord_nod2D(1, :) / rad
!!$        AUX_loc(2, :) = coord_nod2D(2, :) / rad
!!$    ELSE
        AUX_loc(1, :) = coord_nod2D_glob(1, :)
        AUX_loc(2, :) = coord_nod2D_glob(2, :)
!!$    END IF
!!$
    s = s + 1
    !stat(s) = NF_PUT_VAR_REAL(fileid, VarId_loc2D,   AUX_loc(1:2, 1:nod2D)) 
    stat(s) = NF_PUT_VAR_REAL(fileid, VarId_loc2D,   AUX_loc(1:2, 1:nod2D)) 
    s = s + 1

    stat(s) = NF_PUT_VAR_REAL(fileid, VarId_lon,     AUX_loc(1, 1:nod2D) ) 
    s = s + 1
    stat(s) = NF_PUT_VAR_REAL(fileid, VarId_lat,     AUX_loc(2, 1:nod2D) ) 
    s = s + 1

!----- triangles

    !SH Drastic workaround to get undistributed elements
    file_name=trim(meshpath)//'elem2d.out'
    open(101, file=file_name)
    read(101,*) ndummy
    do i=1,elem2D
      read(101,*) elem2D_nodes_glob(:,i)
    end do
    close(101)


    stat(s) = NF_PUT_VAR_INT(fileid,  VarId_trian,    elem2D_nodes_glob(:, :) - 1) 
    s = s + 1
!!$    stat(s) = NF_PUT_VAR_REAL(fileid, VarId_tri_area, REAL(voltriangle(:),        4)) 
!!$    s = s + 1
!!$    stat(s) = NF_PUT_VAR_REAL(fileid, VarId_topo,     REAL(nodhn(:),              4)) 
!!$    s = s + 1
!!$    stat(s) = NF_PUT_VAR_REAL(fileid, VarId_uplift,   REAL(ssh_init(:),           4)) 
!!$    s = s + 1
!!$    stat(s) =NF_PUT_VAR_INT(fileid,   VarId_index,    index_nod2D(:)) 
!!$
!!$    if (write_ttt) then
!!$       s = s + 1
!!$       stat(s) = NF_PUT_VAR_REAL(fileid, VarId_ttt,   ttt(:)) 
!!$    endif
!!$
!!$    DO i = 1, s
!!$        IF (stat(i) /= NF_NOERR) &
!!$            WRITE(*, *)'NetCDF error',stat(i),' in writing variables, no.', i
!!$    END DO

    stat(1) = NF_CLOSE(fileid)

    IF (stat(1) /= NF_NOERR) THEN
        WRITE(*, *) 'NetCDF error',stat(1),' in closing NetCDF file'
    END IF

END SUBROUTINE nc_init
#endif


#ifndef NO_NETCDF
SUBROUTINE nc_out(iteration)

!----- uses

  USE o_PARAM
  USE o_MESH
  USE o_ARRAYS


!----- local declarations

    IMPLICIT NONE

    include 'netcdf.inc'
    
    CHARACTER(len=100) 	         :: ncfile
    INTEGER, INTENT(in) 	 :: iteration
    INTEGER, SAVE 	         :: iter=1
    INTEGER 		         :: FileId
    INTEGER 		         :: VarId_iter, VarId_time   ! numbers
    INTEGER 		         :: VarId_topo               ! topography
    INTEGER 		         :: VarId_ssh                ! fields: ssh
    INTEGER 		         :: VarId_temp               ! fields: temperature
    INTEGER 		         :: VarId_salt               ! fields: salinity
!!$    INTEGER 		         :: VarId_arrival            ! fields: Arrival times
!!$    INTEGER 		         :: VarId_mwh                ! fields: maximum wave height
!!$    INTEGER 		         :: VarId_flux, VarId_vel    !fields: flux and mean absolute velocity
!!$    INTEGER		         :: VarId_nodal_vel 	     !fields: vector
!!$
    INTEGER 		         :: i, n
    INTEGER 		         :: pos1, nmb
    INTEGER, DIMENSION(2)      :: pos1vec, nmbvec
    INTEGER, DIMENSION(3)      :: pos1vec3, nmbvec3
 
    INTEGER,    DIMENSION(50)    :: stat                     ! auxiliaries
    INTEGER                      :: s                        ! auxiliaries
!!$    INTEGER                      :: restart_iter, dimid_iter
    REAL(kind=WP),       DIMENSION(nod2D) :: AUX_n2D
    REAL(kind=WP),       DIMENSION(nsigma-1,nod2D) :: AUX_n3D
!!$    REAL,     DIMENSION(2,nod2D) :: AUX_n2D_2
!!$    REAL                         :: acc, inv_acc
!!$   
!!$    ncfile = ncfile_all

  ncfile='./fesom_c_out.nc'


  stat(1) = NF_OPEN(TRIM(ncfile), NF_WRITE, FileId)

  IF (stat(1) /= NF_NOERR) STOP 'nc-file error'
!!$
!!$    IF (read_restart .AND. iter == 1) THEN
!!$        stat(2) = NF_INQ_DIMID(fileid,  'iteration', dimid_iter)
!!$        stat(2) = NF_INQ_DIMLEN(fileid, dimid_iter,  restart_iter)
!!$        iter    = restart_iter + 1
!!$    END IF
!!$
!----- INQUIRE VARIABLE IDs

  s = 1
  stat(s) = NF_INQ_VARID(fileid, "iteration",         VarId_iter) 
  s = s + 1
  stat(s) = NF_INQ_VARID(fileid, "time",              VarId_time) 
  s = s + 1
  stat(s) = NF_INQ_VARID(fileid, "ssh",               VarId_ssh) 
  s = s + 1
  stat(s) = NF_INQ_VARID(fileid, "temperature",       VarId_temp) 
  s = s + 1
  stat(s) = NF_INQ_VARID(fileid, "salinity",          VarId_salt) 
    !s = s + 1
    !stat(s) = NF_INQ_VARID(fileid, "ssh",               VarId_salt) 
!!$    if (write_diag_in_nc_snapshots) then
!!$       s = s + 1
!!$       stat(s) = NF_INQ_VARID(fileid, "first_arrival",     VarId_arrival) 
!!$       s = s + 1
!!$       stat(s) = NF_INQ_VARID(fileid, "max_wave_height",   VarId_mwh) 
!!$       s = s + 1
!!$       stat(s) = NF_INQ_VARID(fileid, "max_abs_vel",       VarId_vel) 
!!$       s = s + 1
!!$       stat(s) = NF_INQ_VARID(fileid, "max_flux",          VarId_flux) 
!!$    endif
!!$
!!$    IF (write_velocity) THEN
!!$       s = s + 1
!!$       stat(s) = NF_INQ_VARID(fileid, "velocity_in_nodes",     VarId_nodal_vel) 
!!$    END IF
!!$
!!$    DO i = 1, s 
!!$        IF (stat(i) /= NF_NOERR) &
!!$            WRITE(*, *) 'NetCDF error',stat(i),' inquiring variable IDs, no.', i
!!$    END DO
!!$
!----- WRITE VARIABLES
    s = 1
    stat(s) = NF_PUT_VARA_INT(fileid,  VarId_iter, iter, 1, iteration) 

    s = s + 1
    stat(s) = NF_PUT_VARA_REAL(fileid, VarId_time, iter, 1, REAL(time, 4)) 


!----- 2D FIELDS

    pos1vec = (/ 1, iter  /)
    nmbvec  = (/ nod2D, 1 /)
!!$
!!$!----- sea surface height
!!$!----- subtract terrain height on land
!!$

    AUX_n2D = 0.
    AUX_n2D(:) = eta_n_2_glob(:)

!!$    if (nc_snapshot_accuracy_ssh <= 0.) then
!!$!$OMP PARALLEL DO
!!$       DO n = 1, nod2D
!!$          IF (nodhn(n) < 0.) THEN
!!$             AUX_n2D(n) = max(0., ssh1(n) + nodhn(n))
!!$          ELSE
!!$             AUX_n2D(n) = ssh1(n)
!!$          END IF
!!$       END DO
!!$!$OMP END PARALLEL DO
!!$    else
!!$       acc = real(nc_snapshot_accuracy_ssh ,4)
!!$       inv_acc = 1._4/nc_snapshot_accuracy_ssh 
!!$!$OMP PARALLEL DO
!!$       DO n = 1, nod2D
!!$          IF (nodhn(n) < 0.) THEN
!!$             AUX_n2D(n) =  acc* INT(inv_acc*max(0., real(ssh1(n) + nodhn(n),4)))
!!$          ELSE
!!$             AUX_n2D(n) =  acc* INT(inv_acc*real(ssh1(n),4))
!!$          END IF
!!$       END DO
!!$!$OMP END PARALLEL DO
!!$    end if
!!$       
!!$

    s = s + 1
    stat(s) = NF_PUT_VARA_REAL(fileid, VarId_ssh, pos1vec, nmbvec, REAL(AUX_n2D(:),4))
!!$
!!$    if (write_diag_in_nc_snapshots) then
!!$       s = s + 1
!!$       stat(s) = NF_PUT_VARA_REAL(fileid, VarId_arrival, 1, nod2D, arrival_time(:)) 
!!$       s = s + 1
!!$       stat(s) = NF_PUT_VARA_REAL(fileid, VarId_mwh,     1, nod2D, mwh(:))
!!$       s = s + 1
!!$       stat(s) = NF_PUT_VARA_REAL(fileid, VarId_vel,     1, nod2D, max_abs_vel(:))
!!$       s = s + 1
!!$       stat(s) = NF_PUT_VARA_REAL(fileid, VarId_flux,    1, nod2D, max_flux(:))
!!$    endif




    pos1vec3 = (/ 1, 1, iter  /)
    nmbvec3  = (/ nsigma-1, nod2D, 1 /)

    AUX_n3D = 0.
    AUX_n3D   = TF_glob

    s = s + 1
    stat(s) = NF_PUT_VARA_REAL(fileid, VarId_temp, pos1vec3, nmbvec3, REAL(AUX_n3D(:,:),4))

!!$
!!$!----- VECTOR FIELDS
!!$
!!$!----- nodal velocity
!!$
!!$    IF (write_velocity) THEN
!!$       AUX_n2D_2(1,:)=u_node(1:nod2D)
!!$       AUX_n2D_2(2,:)=v_node(1:nod2D)
!!$
!!$       pos1vec3 = (/ 1, 1, iter /)
!!$       nmbvec3  = (/ 2, nod2D, 1 /)
!!$       s = s + 1
!!$       stat(s) = NF_PUT_VARA_REAL(fileid, VarId_nodal_vel, pos1vec3, nmbvec3, AUX_n2D_2)
!!$    END IF
!!$
    DO i = 1, s 
        IF (stat(i) /= NF_NOERR) &
            WRITE(*, *) 'NetCDF error',stat(i),' in writing variables, no.', i
    END DO


!----- CLOSE THE FILE

    stat(1) = NF_CLOSE(fileid)

    IF (stat(1) /= NF_NOERR) THEN
        WRITE(*, *) 'NetCDF error',stat(1),' in closing NetCDF file'
    END IF

    iter = iter + 1

END SUBROUTINE nc_out
#endif


#ifndef NO_NETCDF
subroutine nc_finalize

!----- uses

  USE o_PARAM
  USE o_MESH
  USE o_ARRAYS


!----- local declarations

    IMPLICIT NONE

    include 'netcdf.inc'
    
!!$    CHARACTER(len=100) 	    :: ncfile
!!$    INTEGER 		    :: FileId
!!$    INTEGER                 :: stat(32), s, i           ! auxiliaries
!!$    REAL(KIND=8)            :: attVal
!!$    INTEGER 		    :: VarId_arrival            ! fields: Arrival times
!!$    INTEGER 		    :: VarId_mwh                ! fields: maximum wave height
!!$    INTEGER 		    :: VarId_flux, VarId_vel    !fields: flux and mean absolute velocity
!!$
!!$    ncfile=ncfile_all
!!$
!!$    stat(1) = NF_OPEN(TRIM(ncfile), NF_WRITE, FileId)
!!$
!!$    if (stat(1) /= NF_NOERR) stop 'nc-file error'
!!$
!!$    attVal=-99999.
!!$    stat(1) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'globalMinSSH',         NF_DOUBLE, 1, attVal) 
!!$    if (stat(1) /= NF_NOERR) WRITE(*, *) 'NetCDF error',stat(1),' in rewriting global attribute "globalMinSSH"'
!!$
!!$    attVal=maxval(mwh(:))
!!$    stat(1) = NF_PUT_ATT_DOUBLE(fileid, NF_GLOBAL, 'globalMaxSSH',         NF_DOUBLE, 1, attVal) 
!!$    if (stat(1) /= NF_NOERR) WRITE(*, *) 'NetCDF error',stat(1),' in rewriting global attribute "globalMaxSSH"'
!!$    
!!$    if (.not. write_diag_in_nc_snapshots) then
!!$       
!!$       if (verbosity >= 1) print *,'Writing MWH, ETA, Max Vel, Max Flux'
!!$
!!$       s = 1
!!$       stat(s) = NF_INQ_VARID(fileid, "first_arrival",     VarId_arrival) 
!!$       s = s + 1
!!$       stat(s) = NF_INQ_VARID(fileid, "max_wave_height",   VarId_mwh) 
!!$       s = s + 1
!!$       stat(s) = NF_INQ_VARID(fileid, "max_abs_vel",       VarId_vel) 
!!$       s = s + 1
!!$       stat(s) = NF_INQ_VARID(fileid, "max_flux",          VarId_flux) 
!!$
!!$       DO i = 1, s 
!!$          IF (stat(i) /= NF_NOERR) &
!!$               WRITE(*, *) 'NetCDF error',stat(i),', inquiring variable IDs, no.', i
!!$       END DO
!!$
!!$       s = 1
!!$       stat(s) = NF_PUT_VARA_REAL(fileid, VarId_arrival, 1, nod2D, arrival_time(:)) 
!!$       s = s + 1
!!$       stat(s) = NF_PUT_VARA_REAL(fileid, VarId_mwh,     1, nod2D, mwh(:))
!!$       s = s + 1
!!$       stat(s) = NF_PUT_VARA_REAL(fileid, VarId_vel,     1, nod2D, max_abs_vel(:))
!!$       s = s + 1
!!$       stat(s) = NF_PUT_VARA_REAL(fileid, VarId_flux,    1, nod2D, max_flux(:))
!!$      
!!$       DO i = 1, s 
!!$          IF (stat(i) /= NF_NOERR) &
!!$               WRITE(*, *) 'NetCDF error,',stat(i),' writing variables, no.', i
!!$       END DO
!!$    endif
!!$
!!$
!!$    stat(1) = NF_CLOSE(fileid)
!!$    if (stat(1) /= NF_NOERR) write(*, *) 'NetCDF error',stat(1),' in closing NetCDF file'
    

end subroutine nc_finalize
#endif
