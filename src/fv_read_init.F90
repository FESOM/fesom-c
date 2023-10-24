!==============================================================================|
!   Input parameters and flags, which control the simulation                               |
!==============================================================================|
!==============================================================================|
!VF, updated 02.09
SUBROUTINE READ_DATA_INIT()           

USE o_PARAM
USE o_MESH
USE o_UTILIT

   IMPLICIT NONE
   INTEGER ::  ISCAN 
   CHARACTER(LEN=120) :: FNAME

   FNAME = "./fv_run.dat"


!------------------------------------------------------------------------------|
!  Type of coordinates   
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"cartesian",LVAL = cartesian)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING cartesian: ',ISCAN
   END IF   

!------------------------------------------------------------------------------|
!     "MESHPATH"             
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"MESHPATH",CVAL = MESHPATH)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING MESHPATH: ',ISCAN
     STOP 
   END IF

!==============================================================================|
!   CASENAME                                         |
!==============================================================================|

   ISCAN = SCAN_FILE(TRIM(FNAME),"TITLE",CVAL = TITLE)
   IF(ISCAN /= 0)THEN
     WRITE(*,*)'ERROR READING TITLE: ',ISCAN
     STOP 
   END IF 
!------------------------------------------------------------------------------|
!     "cyclic_length"     
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(FNAME,"cyclic_length",FSCAL = cyclic_length)
   IF((ISCAN /= 0).or. (cyclic_length>360).or. (cyclic_length<=0)) THEN
     WRITE(*,*)'ERROR READING cyclic_length: ',ISCAN
     STOP
   END IF
   cyclic_length=cyclic_length*rad
   
!==============================================================================|
!            SCREEN REPORT OF SET VARIABlES                                    !
!==============================================================================|

   if (cartesian) then 
   WRITE(*,*)'!  # Type of the coordinates             : cartesian '
   else
   WRITE(*,*)'!  # Type of the coordinates             : spherical '
   endif
   WRITE(*,*)'!  # MESHPATH                            : ',trim(MESHPATH)

   RETURN
 
END SUBROUTINE READ_DATA_INIT 

