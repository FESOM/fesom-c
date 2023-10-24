!VF, updated 19.10
!===================================================================================|
MODULE o_UTILIT

contains
   SUBROUTINE binarysearch(length, array, value, ind)!, delta)
      ! Given an array and a value, returns the index of the element that
      ! is closest to, but less than, the given value.
      ! Uses a binary search algorithm.
      ! "delta" is the tolerance used to determine if two values are equal
      ! if ( abs(x1 - x2) <= delta) then
      !    assume x1 = x2
      ! endif
      !org. source from: https://github.com/cfinch/Shocksolution_Examples/blob/master/FORTRAN/BilinearInterpolation/interpolation.f90
      USE o_PARAM
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

 FUNCTION Dist2in1(X1, X2)

!===================================================================================|
!   This function calculates the distance from point(X1,X2) to the open boundary    |
!   or any other curve represented by set of points.                                |
!   if the coordinates are given in Cartesian coord., then usual Eulerian norm is   |
!   used to calculate the distance.                                                 |
!   If the coordinates are given in spherical coord., than Haversine Formula is used|
!===================================================================================|

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  IMPLICIT NONE

  REAL(kind=WP)            :: Dist2in1
  real(kind=WP),INTENT(IN) :: X1, X2
  real(kind=WP)            :: a, d, cos_X2
  integer                  :: n,j

   if (cartesian) then
      d = minval((X1obn(:) - X1)**2 + (X2obn(:) - X2)**2)
      Dist2in1 = sqrt(d) *r_earth
   else

      cos_X2 = cos(X2)

      a = minval( sin((X2obn(:) - X2)*0.5_WP)**2 &
                + cos_X2 * cos(X2obn(:)) * sin((X1obn(:) - X1)*0.5_WP)**2 )

      Dist2in1 = 2._WP *  r_earth * asin(min(1.0_WP,sqrt(a)))
   endif
   RETURN

   END FUNCTION Dist2in1

function scal_diff_grad(var1,elnodes)
!=====================================================================================|
USE o_MESH
USE o_ARRAYS
USE o_PARAM
   IMPLICIT NONE
   REAL(kind=WP)            :: scal_diff_grad(8)
   REAL(kind=WP)            :: varg(8)
   real(kind=WP),INTENT(IN) :: var1(4)
   integer,INTENT(IN)       ::elnodes(4)
   real(kind=WP)            :: delta31, delta21, x1,x2,xd
   integer                  :: enodes(6),j

   if(elnodes(1)==elnodes(4)) then

   delta31=var1(3)-var1(1)
   delta21=var1(2)-var1(1)

   varg(1)=(delta31-delta21)
   varg(2)=-delta31
   varg(3)= delta21
   varg(4)=0.0_WP
   varg(5:8)=-varg(1:4)
   else
   enodes(2:5)=(/1,2,3,4/)
   enodes(1)=enodes(5)
   enodes(6)=enodes(2)

   varg(:)=0.0_WP
   DO j=2,5    ! Nodes are listed clockwise n=(-dy, dx)
   varg(j-1)=(var1(enodes(j+1))-var1(enodes(j-1)))
   varg(j-1+4)=(var1(enodes(j+1))-var1(enodes(j-1)))
   end do
   end if
scal_diff_grad=varg
return
 end function scal_diff_grad

!=====================================================================================|
subroutine update_info_sed(sk)

USE o_MESH
USE o_ARRAYS
USE o_PARAM
   IMPLICIT NONE
   integer,INTENT(INOUT)    ::sk

if ((sk<size(sed_vt2)).and.(sed_vt2(sk)<dt*dble(n_dt)/3600.0_WP)) then
write(*,*) 'sed. conc. at the ob'
        sk=sk+1
	cr_distr2=cr_distr(:,:,sk)
else
if (sed_vt2(sk)<dt*dble(n_dt)/3600.0_WP) then
sed_boun_flux=.FALSE.
write(*,*) 'No information available any more for the sed. conc. at the ob'
endif
endif

end subroutine update_info_sed


subroutine ac_create_Xobn

!=====================================================================================|
!   This function calculates the coefficient ac (0<=ac<=1), which will be multiplied  |
! on advection and diffusion terms within the sponge layer at the open boundary.      |
! 1 - keep full terms, 0 - kill full terms                                            |
!=====================================================================================|

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  use g_parsup

  IMPLICIT NONE

  integer                    :: n,j

  j=0
  do n=1,nod2D
     if (index_nod2D_glob(n) == 2) then
        j= j + 1
        X1obn(j)=coord_nod2D_glob(1,n)
        X2obn(j)=coord_nod2D_glob(2,n)
        in_obn(j)=n
     endif
  enddo
  if (j/=nobn) then
     print *,'WARNING AC node discrepancy:',mype,j,nobn
  end if

END subroutine ac_create_Xobn

subroutine ac_create(gc)

!=====================================================================================|
!   This function calculates the coefficient ac (0<=ac<=1), which will be multiplied  |
! on advection and diffusion terms within the sponge layer at the open boundary.      |
! 1 - keep full terms, 0 - kill full terms                                            |
!=====================================================================================|

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  use g_parsup

  IMPLICIT NONE

  REAL(kind=WP),INTENT(INOUT):: gc(myDim_nod2D+eDim_nod2D)
  real(kind=WP)              :: d2obn
  integer                    :: n,j, fID

!SHDB   print *,'Entering ac_create',mynobn,nobn
  do n=1,myDim_nod2D+eDim_nod2D

     d2obn = Dist2in1( coord_nod2D(1,n),  coord_nod2D(2,n))

     if (d2obn < SL_radius) gc(n)=d2obn/SL_radius

  enddo

!SHDB   print *,'AC_CREATE: ',mype,minval(ac),maxval(ac)


!do n=myDim_nod2D+1,myDim_nod2D+eDim_nod2D
!print *,mype,myList_nod2D(n),gc(n)
!end do

END subroutine ac_create


 FUNCTION SCAN_FILE(FNAME,VNAME,ISCAL,FSCAL,IVEC,FVEC,CVEC,NSZE,CVAL,LVAL)

!==============================================================================|
!   This function is taken from FVCOM source code                              |
!   Scan an Input File for a Variable                                          |
!   RETURN VALUE:                                                              |
!        0 = FILE FOUND, VARIABLE VALUE FOUND                                  |
!       -1 = FILE DOES NOT EXIST OR PERMISSIONS ARE INCORRECT                  |
!       -2 = VARIABLE NOT FOUND OR IMPROPERLY SET                              |
!       -3 = VARIABLE IS OF DIFFERENT TYPE, CHECK INPUT FILE                   |
!       -4 = VECTOR PROVIDED BUT DATA IS SCALAR TYPE                           |
!       -5 = NO DATATYPE DESIRED, EXITING                                      |
!							                       |
!   REQUIRED INPUT:		        				       |
!        FNAME = File Name					               |
!        FSIZE = Length of Filename					       |
!                                                                              |
!   OPTIONAL (MUST PROVIDE ONE)        					       |
!        ISCAL = INTEGER SCALAR					               |
!        FSCAL = FLOAT SCALAR  						       |
!        CVAL = CHARACTER VARIABLE                                             |
!        LVAL = LOGICAL VARIABLE                                               |
!        IVEC = INTEGER VECTOR **                                              |
!        FVEC = FLOAT VECTOR **                                                |
!        CVEC = STRING VECTOR **                         |
!      **NSZE = ARRAY SIZE (MUST BE PROVIDED WITH IVEC/FVEC)                   |
!                                                                              |
!==============================================================================|

USE o_PARAM

   IMPLICIT NONE
   CHARACTER(LEN=*) :: FNAME,VNAME
   INTEGER, INTENT(INOUT), OPTIONAL :: ISCAL,IVEC(*)
   REAL(kind=WP),INTENT(INOUT), OPTIONAL :: FSCAL,FVEC(*)
   CHARACTER(LEN=3), OPTIONAL      :: CVEC(*)
   CHARACTER(LEN=80), OPTIONAL      :: CVAL
   LOGICAL, INTENT(INOUT), OPTIONAL :: LVAL
   INTEGER, INTENT(INOUT), OPTIONAL :: NSZE

!------------------------------------------------------------------------------|

   INTEGER :: SCAN_FILE
   REAL(kind=WP) REALVAL(150)
   INTEGER  INTVAL(150)
   CHARACTER(LEN=20 ) :: VARNAME
   CHARACTER(LEN=80 ) :: STRINGVAL(150)
   CHARACTER(LEN=80 ) :: INPLINE
   CHARACTER(LEN=400) :: TLINE
   CHARACTER(LEN=7  ) :: VARTYPE
   CHARACTER(LEN=20 ), DIMENSION(200)  :: SET
   INTEGER I,NVAL,J,NSET,NLINE,NREP
   LOGICAL SETYES,ALLSET,CHECK,LOGVAL


   SCAN_FILE = 0
!==============================================================================|
!            OPEN THE INPUT FILE                                               |
!==============================================================================|
   INQUIRE(FILE=TRIM(FNAME),EXIST=CHECK)
   IF(.NOT.CHECK)THEN
     SCAN_FILE = -1
     RETURN
   END IF

   OPEN(10,FILE=TRIM(FNAME)) ; REWIND(10)

!==============================================================================|
!            SCAN THE FILE FOR THE VARIABLE NAME                               |
!==============================================================================|

   NSET = 0
   NLINE = 0
   DO WHILE(.TRUE.)
     TLINE(1:LEN(TLINE)) = ' '
     NREP  = 0
     NLINE = NLINE + 1
     READ(10,'(a)',END=20) INPLINE
     TLINE(1:80) = INPLINE(1:80)

!----PROCESS LINE CONTINUATIONS------------------------------------------------!
 110 CONTINUE
     I = LEN_TRIM(INPLINE)
     IF(I > 1)THEN
     IF( INPLINE(I-1:I) == '\\')THEN
       NREP = NREP + 1
       READ(10,'(a)',END=20) INPLINE
       NLINE = NLINE + 1
       TLINE( NREP*80 + 1 : NREP*80 +80) = INPLINE(1:80)
       GOTO 110
     END IF
     END IF
     IF(NREP > 4) WRITE(*,*) 'CANNOT HAVE > 4 LINE CONTINUATIONS'

!----REMOVE LINE CONTINUATION CHARACTER \\-------------------------------------!
     IF(NREP > 0)THEN
       DO I=2,LEN_TRIM(TLINE)
         IF( TLINE(I-1:I) == '\\') TLINE(I-1:I) = '  '
       END DO
     END IF

!----PROCESS THE LINE----------------------------------------------------------!
     CALL GET_VAL(NLINE,LEN_TRIM(TLINE),ADJUSTL(TLINE),VARNAME,VARTYPE,LOGVAL,&
                 STRINGVAL,REALVAL,INTVAL,NVAL)

!----IF VARNAME MATCHES, PROCESS VARIABLE AND ERROR-CHECK----------------------!

     IF(TRIM(VARNAME) == TRIM(VNAME))THEN

       IF(PRESENT(ISCAL))THEN
         IF(VARTYPE == 'integer')THEN
           ISCAL = INTVAL(1)
		   CLOSE(10)
           RETURN
         ELSE
           SCAN_FILE = -3
         END IF
       ELSE IF(PRESENT(FSCAL))THEN
         IF(VARTYPE == 'float')THEN
           FSCAL = REALVAL(1)
		   CLOSE(10)
           RETURN
         ELSE
           SCAN_FILE = -3
         END IF
       ELSE IF(PRESENT(CVAL))THEN
         IF(VARTYPE == 'string')THEN
           CVAL = STRINGVAL(1)
		!   write (*,*), cval
           CLOSE(10)
           RETURN
         ELSE
           SCAN_FILE = -3
         END IF
       ELSE IF(PRESENT(LVAL))THEN
         IF(VARTYPE == 'logical')THEN
           LVAL = LOGVAL
		   CLOSE(10)
           RETURN
         ELSE
           SCAN_FILE = -3
         END IF
       ELSE IF(PRESENT(IVEC))THEN
         IF(NVAL > 1)THEN
           IF(VARTYPE == 'integer')THEN
             IVEC(1:NVAL) = INTVAL(1:NVAL)
             NSZE = NVAL
             CLOSE(10)
             RETURN
           ELSE
             SCAN_FILE = -3
           END IF
           ELSE
           SCAN_FILE = -4
         END IF
       ELSE IF(PRESENT(FVEC))THEN
         IF(NVAL > 1)THEN
           IF(VARTYPE == 'float')THEN
             FVEC(1:NVAL) = REALVAL(1:NVAL)
             NSZE = NVAL
             CLOSE(10)
             RETURN
           ELSE
             SCAN_FILE = -3
           END IF
         ELSE
           SCAN_FILE = -4
         END IF
       ELSE IF(PRESENT(CVEC))THEN
         IF(NVAL > 0)THEN
           IF(VARTYPE == 'string')THEN
             CVEC(1:NVAL) = STRINGVAL(2:NVAL+1)
             NSZE = NVAL
			 CLOSE(10)
             RETURN
           ELSE
             SCAN_FILE = -3
           END IF
         ELSE
           SCAN_FILE = -4
         END IF
       ELSE
         SCAN_FILE = -5
       END IF
     END IF  !!VARIABLE IS CORRECT

   END DO !!LOOP OVER INPUT FILE
 20 CLOSE(10)
   SCAN_FILE = -2
   RETURN

   END FUNCTION SCAN_FILE

   !==============================================================================!
!  DECOMPOSE INPUT LINE INTO VARIABLE NAME AND VARIABLE VALUE(S)               !
!==============================================================================!

SUBROUTINE GET_VAL(LNUM,NUMCHAR,TEXT_LINE,VARNAME,VARTYPE,LOGVAL,STRINGVAL,&
                    REALVAL,INTVAL,NVAL)

!==============================================================================!
USE o_PARAM
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: LNUM,NUMCHAR
  CHARACTER(LEN=NUMCHAR) :: TEXT_LINE
  CHARACTER(LEN=20), INTENT(OUT) :: VARNAME
  CHARACTER(LEN=7), INTENT(OUT) :: VARTYPE
  LOGICAL, INTENT(OUT) :: LOGVAL
  CHARACTER(LEN=80), INTENT(OUT) :: STRINGVAL(150)
  REAL(kind=WP), INTENT(INOUT) :: REALVAL(150)
  INTEGER, INTENT(INOUT) :: INTVAL(150)
  INTEGER, INTENT(OUT) :: NVAL
!------------------------------------------------------------------------------!
  CHARACTER(LEN=NUMCHAR) :: VARVAL,TEMP,FRAG(200)
  CHARACTER(LEN=80) :: TSTRING
  CHARACTER(LEN=6) :: ERRSTRING
  CHARACTER(LEN=16) :: NUMCHARS
  INTEGER LENGTH,EQLOC,LVARVAL,DOTLOC
  INTEGER I,J,LOCEX,NP
  LOGICAL ONFRAG

!==============================================================================!
  FRAG = " "
  NUMCHARS = "0123456789+-Ee. "
  VARTYPE = "error"
  LOGVAL = .FALSE.
  LENGTH = LEN_TRIM(TEXT_LINE)
  WRITE(ERRSTRING,"(I6)") LNUM
  LOCEX = INDEX(TEXT_LINE,"!")

!
!-----------------------CHECK FOR BLANK LINE OR COMMENT------------------------!
!
  IF(LENGTH == 0 .OR. LOCEX==1)THEN
    VARTYPE = "no data"
    VARNAME = "no data"
    RETURN
  END IF

!
!-----------------------CHANGE COMMAS TO BLANKS--------------------------------!
!
  DO I=1,LENGTH
    IF(TEXT_LINE(I:I) == ",") TEXT_LINE(I:I) = " "
  END DO
!
!-----------------------REMOVING TRAILING COMMENTS-----------------------------!
!
  IF(LOCEX /= 0)THEN
    TEMP = TEXT_LINE(1:LOCEX-1)
    TEXT_LINE = TEMP
   END IF
!
!--------------------ENSURE "=" EXISTS AND DETERMINE LOCATION------------------!
!
   EQLOC = INDEX(TEXT_LINE,"=")
  ! IF((EQLOC == 0).and.(length >0)) WRITE(*,*) 'DATA LINE MUST CONTAIN "=" ', TEXT_LINE

!
!--------------------SPLIT OFF VARNAME AND VARVAL STRINGS----------------------!
!
   VARNAME = TEXT_LINE(1:EQLOC-1)
   VARVAL  = ADJUSTL(TEXT_LINE(EQLOC+1:LENGTH))
   LVARVAL = LEN_TRIM(VARVAL)
   IF(LVARVAL == 0) WRITE(*,*) 'IN DATA PARAMETER FILE VARIABLE LINE HAS NO ASSOCIATED VALUE'
!
!-----------------DETERMINE TYPE OF VARVAL-------------------------------------!
!

!
!  CHECK FOR LOGICAL
!
   IF((VARVAL(1:1) == "T" .OR. VARVAL(1:1) == "F") .AND. LVARVAL == 1)THEN
     VARTYPE = "logical"
     IF(VARVAL(1:1) == "T") LOGVAL = .TRUE.
     RETURN
   END IF

!
!  CHECK IF IT IS A STRING  (CONTAINS CHARACTERS OTHER THAN 0-9,+,-,e,E,.)
!
   DO I=1,LVARVAL
     IF(INDEX(NUMCHARS,VARVAL(I:I)) == 0) VARTYPE = "string"
   END DO

!
!  PROCESS STRING (MAY BE MULTIPLE)
!
   IF(VARTYPE == "string") THEN
     TSTRING = VARVAL
     STRINGVAL(1) = TSTRING
     NVAL = 1
     ONFRAG = .TRUE.
     DO I=1,LVARVAL
       IF(VARVAL(I:I) /= " ")THEN
         FRAG(NVAL) = TRIM(FRAG(NVAL))//VARVAL(I:I)
         ONFRAG = .TRUE.
       ELSE
         IF(ONFRAG) NVAL = NVAL + 1
         ONFRAG = .FALSE.
       END IF
     END DO
     DO I=1,NVAL
       STRINGVAL(I+1) = TRIM(FRAG(I))
     END DO
     RETURN
   END IF

!
!  CHECK IF IT IS A FLOAT
!

   DOTLOC = INDEX(VARVAL,".")
   IF(DOTLOC /= 0) THEN
     VARTYPE = "float"
   ELSE
     VARTYPE = "integer"
   END IF
!
!-----------------FRAGMENT INTO STRINGS FOR MULTIPLE VALUES---------------------!
!
   NP = 1
   ONFRAG = .TRUE.
   DO I=1,LVARVAL
     IF(VARVAL(I:I) /= " ")THEN
       FRAG(NP) = TRIM(FRAG(NP))//VARVAL(I:I)
       ONFRAG = .TRUE.
     ELSE
       IF(ONFRAG) NP = NP + 1
       ONFRAG = .FALSE.
     END IF
   END DO
!
!-----------------EXTRACT NUMBER(S) FROM CHARACTER STRINGS----------------------!
!

   NVAL = NP
   DO I=1,NP
     TEMP = TRIM(FRAG(I))
     IF(VARTYPE == "float") THEN
       READ(TEMP,*)REALVAL(I)
     ELSE
       READ(TEMP,*)INTVAL(I)
     END IF
   END DO

END SUBROUTINE GET_VAL

END MODULE o_UTILIT

 !============================================================

subroutine compute_el2nodes_2D(var_e,var_n)

  USE o_MESH
  USE o_PARAM
  USE o_ARRAYS

  USE g_PARSUP
  use g_comm_auto

  IMPLICIT NONE

  integer            :: n, num

  REAL(kind=WP),INTENT(IN) :: var_e(myDim_elem2D)
  REAL(kind=WP),INTENT(OUT):: var_n(myDim_nod2D)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n,num)
  DO n=1, myDim_nod2D

     num = nod_in_elem2D_num(n)

     var_n(n) = sum(var_e(nod_in_elem2D(1:num,n))*elem_area(nod_in_elem2D(1:num,n))) &
                                          / sum(elem_area(nod_in_elem2D(1:num,n)))
  END DO
!$OMP END PARALLEL DO

#ifdef USE_MPI
  call exchange_nod(var_n)
#endif
!SHDB print *,'var_n:',mype,minval(var_n),maxval(var_n)

end subroutine compute_el2nodes_2D

!============================================================

subroutine compute_el2nodes_2D_2fields(var1_e,var2_e,var1_n,var2_n)
USE o_MESH
USE o_PARAM
USE o_ARRAYS
!USE g_PARSUP
IMPLICIT NONE

REAL(kind=WP),INTENT(IN) :: var1_e(elem2D), var2_e(elem2D)
REAL(kind=WP),INTENT(OUT):: var1_n(nod2D),  var2_n(nod2D)

integer                  :: n, num
real(kind=WP)            :: area_inv

!$OMP DO
DO n=1, nod2D

   num = nod_in_elem2D_num(n)
   area_inv = 1._WP / sum(elem_area(nod_in_elem2D(1:num,n)))

   var1_n(n) = sum(var1_e(nod_in_elem2D(1:num,n))*elem_area(nod_in_elem2D(1:num,n))) &
              * area_inv

   var2_n(n) = sum(var2_e(nod_in_elem2D(1:num,n))*elem_area(nod_in_elem2D(1:num,n))) &
              * area_inv
END DO
!$OMP END DO
end subroutine compute_el2nodes_2D_2fields

!============================================================

  ! modification by Alexey Androsov
  ! (for sigma coordinates)
  ! 08.10.14

subroutine compute_el2nodes_3D(Ue,Ve,Un,Vn)

  USE o_MESH
  USE o_PARAM
  USE o_ARRAYS

  USE g_PARSUP
  use g_comm_auto

  IMPLICIT NONE

  real(kind=WP), intent(in)  :: Ue(nsigma-1,mydim_elem2D+eDim_elem2D+eXDim_elem2D)
  real(kind=WP), intent(in)  :: Ve(nsigma-1,mydim_elem2D+eDim_elem2D+eXDim_elem2D)
  real(kind=WP), intent(out) :: Un(nsigma-1,myDim_nod2D+eDim_nod2D)
  real(kind=WP), intent(out) :: Vn(nsigma-1,myDim_nod2D+eDim_nod2D)

  integer            :: n, nz, k, el,num
  real(kind=WP)      :: tvol_inv

!SHDB print *,'Ue:',mype,minval(Ue),maxval(Ue)
!SHDB print *,'Ve:',mype,minval(Ve),maxval(Ve)



!$OMP DO
  DO n=1, myDim_nod2D

     num = nod_in_elem2D_num(n)

     tvol_inv = 1._WP/sum(elem_area(nod_in_elem2D(1:num,n)))

     Un(1:nsigma-1,n) = 0._WP
     Vn(1:nsigma-1,n) = 0._WP

     DO k=1,num
        el=nod_in_elem2D(k,n)

        Un(1:nsigma-1,n) = Un(1:nsigma-1,n) + Ue(1:nsigma-1,el)*elem_area(el)
        Vn(1:nsigma-1,n) = Vn(1:nsigma-1,n) + Ve(1:nsigma-1,el)*elem_area(el)

     END DO

     Un(1:nsigma-1,n) = Un(1:nsigma-1,n)*tvol_inv
     Vn(1:nsigma-1,n) = Vn(1:nsigma-1,n)*tvol_inv

  END DO
!$OMP END DO

#ifdef USE_MPI
  call exchange_nod(Un)
  call exchange_nod(Vn)
#endif

!SHDB print *,'UnA:',mype,minval(Un),maxval(Un)
!SHDB print *,'VnA:',mype,minval(Vn),maxval(Vn)

end subroutine compute_el2nodes_3D

!====================================================================
