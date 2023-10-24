SUBROUTINE potential
!AA
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
use g_comm_auto
  IMPLICIT NONE
!шведерский (schvederskiy) приливы
!kowalik proshutinskiy Barents sea tides

  REAL(kind=WP) :: DY, TY, TY2, TY3
  REAL(kind=WP) :: h0, s0, p0, ta
  REAL(kind=WP) :: Apt(11), Fpt(11), sxi(11)
  REAL(kind=WP) :: year, ssh_sd, ssh_lp, ssh_d, day_number, c,time_cor
  INTEGER :: iz, n, j, i

!     amplitude

  DATA Apt/0.242334_WP,0.112841_WP,0.046398_WP,0.030704_WP & ! M2, S2, N2, K2
       ,0.141565_WP,0.100514_WP,0.046843_WP,0.019256_WP &          ! K1, O1, P1, Q1
       ,0.041742_WP,0.022026_WP,0.019446_WP/                              ! Mf, Mm, Ssa

!     phase

  DATA Fpt/1.40519d-4,1.45444d-4,1.37880d-4,1.45842d-4  &
       ,0.72921d-4,0.67598d-4,0.72523d-4,0.64959d-4          &
       ,0.53234d-5,0.26392d-5,0.03982d-5/

  !     day_number ---> is the day number of the year (day_number = 1 for January 1)
  !     year >= 1975 ---> is the year number

  !year = 1975.0_WP
  !day_number=FLOOR(time_jd/(24._WP*3600._WP))+1.0_WP


  !DY = day_number + 365._WP*(year - 1975._WP) + floor((year - 1975._WP)/4._WP)

  DY =  time_jd - 2442413.5_WP  ! days from 01.01.1975 year

  TY = (27392.500528_WP + 1.0000000356_WP*DY)/36525.0_WP
  TY2 = TY*TY
  TY3 = TY2*TY

  h0 = 279.69668_WP + 36000.768930485_WP*TY + 3.03d-4*TY2
  CALL control_360(iz,h0)

  s0 = 270.434358_WP + 481267.88314137_WP*TY - 0.001133_WP*TY2 + &
       1.9d-6*TY3
  CALL control_360(iz,s0)

  p0 = 334.329653_WP + 4069.0340329575_WP*TY - 0.010325_WP*TY2 -&
       1.2d-5*TY3
  CALL control_360(iz,p0)

  sxi(1) = (h0 - s0)*2.0_WP*rad                        ! M2
  sxi(2) = 0.0_WP                                            ! S2
  sxi(3) = (2.0_WP*h0 - 3.0_WP*s0 + p0)*rad     ! N2
  sxi(4) = 2.0_WP*h0*rad                                 ! K2

  sxi(5) = (h0 + 90.0_WP)*rad                          ! K1
  sxi(6) = (h0 - 2.0_WP*s0 - 90.0_WP)*rad         ! O1
  sxi(7) = (h0 - 90.0_WP)*rad                           ! P1
  sxi(8) = (h0 - 3.0_WP*s0 + p0 - 90.0_WP)*rad  ! Q1

  sxi(9) = 2.0_WP*s0*rad                                  ! Mf
  sxi(10)= (s0 - p0)*rad                                    ! Mm
  sxi(11)= 2.0_WP*h0*rad                                 ! Ssa

  ! calculation equilibrium tides
  ! ssh_sd ----> semidiurnal (M2, S2, N2, K2)
  ! ssh_d  ----> diurnal     (K1, O1, P1, Q1)
  ! ssh_lp ----> long-period (Mf, Mm, Ssa)


  DO n=1,myDim_nod2D+eDim_nod2D

     ssh_sd = 0.0_WP
     DO j=1,4
        time_cor = DY*86400.0_WP  - sxi(j)/Fpt(j)
        ssh_sd = ssh_sd + a_th(j)*Apt(j)*SIN(0.5_WP*pi - coord_nod2d(2,n))*SIN(0.5_WP*pi - coord_nod2d(2,n))&
             *COS(Fpt(j)*time_cor + 2.0_WP*coord_nod2d(1,n) + sxi(j))
     ENDDO
     ssh_d = 0.0_WP
     DO j=5,8
            time_cor = DY*86400.0_WP - sxi(j)/Fpt(j)
            ssh_d = ssh_d + a_th(j)*Apt(j)*SIN(2.0_WP*(0.5_WP*pi - coord_nod2d(2,n)))&
              *COS(Fpt(j)*time_cor + coord_nod2d(1,n) + sxi(j))
     ENDDO
     ssh_lp = 0.0_WP
     DO j=9,11
            time_cor = DY*86400.0_WP - sxi(j)/Fpt(j)
                 ssh_lp = ssh_lp + a_th(j)*Apt(j)*(3.0_WP*SIN(0.5_WP*pi - coord_nod2d(2,n))**2.0_WP - 2.0_WP)&
                      *COS(Fpt(j)*time_cor + sxi(j))
     ENDDO
!++++++++++++++++++++++++++++++++++++++++
! equilibrium tide for 11 harmonics
!++++++++++++++++++++++++++++++++++++++++

     ssh_gp(n) = ssh_sd + ssh_d + ssh_lp    !!! now, for pressure gradient (ssh - ssh_gp)
     if (time_jd-time_jd0<5.0_WP) then
            ssh_gp(n)=ssh_gp(n)*(time_jd-time_jd0)/5.0_WP
     end if

  ENDDO

END SUBROUTINE potential
!===========================================================================
SUBROUTINE control_360(iz,fr)

USE o_PARAM
  implicit none

  REAL(KIND=WP) :: fr
  integer :: iz

  iz = INT(fr/360.0_WP)
  IF (iz >= 1) fr = fr - float(iz)*360.0_WP

END SUBROUTINE control_360
