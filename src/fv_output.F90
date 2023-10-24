subroutine out
USE o_MESH
USE o_ARRAYS
USE o_PARAM
IMPLICIT NONE

REAL(KIND=8)  :: elnodes(4), x, y

INTEGER :: n, elem

! ascii output for visualization
open(10, file='ssh.dat')
open(11, file='vel.dat')
do n=1,nod2D
write(10,'(3f12.4)') coord_nod2d(1,n)/rad,coord_nod2d(2,n)/rad,eta_n(n)
enddo
 DO elem=1, elem2D
  elnodes=elem2D_nodes(:,elem)
  x = .25_WP*sum(coord_nod2d(1,elem2D_nodes(:,elem)))/rad
  y = .25_WP*sum(coord_nod2d(2,elem2D_nodes(:,elem)))/rad
write(11,'(2f12.4,1x,f12.4, 1x, f12.4)') x,y,U_n_2D(1,elem),U_n_2D(2,elem)
enddo

close(10)
close(11)
end subroutine out
!========================================================================================
subroutine cross_sec
  use o_MESH
  use o_ARRAYS
  use o_PARAM
  !
  implicit none
  integer                            :: i, n, nz, elem, elnodes(4)
  real(kind=WP)                       :: x,y,x1,x2,y1,y2,xy
  character*23                       :: cross_file='output_TF\cross0000.dat'
  character*25                       :: cross_file20='output_TF20\cross0000.dat'
  character*23                       :: cross_fileU='output_Un\cross0000.dat'
  character*23                       :: cross_vort2D='output_V2\cross0000.dat'
  character*23                       :: cross_vort3D='output_V3\cross0000.dat'
  character*23                       :: cross_SSH='output_SH\cross0000.dat'
  character*23                       :: cross_fileBar='output_Br\cross0000.dat'
!
  write(cross_file(16:19),'(I4.4)') idnint(time/3600.0_WP/24.0_WP/1.0_WP)
  write(cross_file20(18:21),'(I4.4)') idnint(time/3600.0_WP/24.0_WP/1.0_WP)
  write(cross_fileU(16:19),'(I4.4)') idnint(time/3600.0_WP/24.0_WP/1.0_WP)
  write(cross_vort2D(16:19),'(I4.4)') idnint(time/3600.0_WP/24.0_WP/1.0_WP)
  write(cross_vort3D(16:19),'(I4.4)') idnint(time/3600.0_WP/24.0_WP/1.0_WP)
  write(cross_SSH(16:19),'(I4.4)') idnint(time/3600.0_WP/24.0_WP/1.0_WP)
  write(cross_fileBar(16:19),'(I4.4)') idnint(time/3600.0_WP/24.0_WP/1.0_WP)

  open(1,file=cross_file)
  open(2,file=cross_fileU)
  open(3,file=cross_vort2D)
  open(4,file=cross_vort3D)
  open(5,file=cross_SSH)
 ! open(9,file=cross_file20)
!  open(8,file=cross_fileBar)
  do n = 1,nod2D
    y=coord_nod2D(2, n)/rad
    x=coord_nod2D(1, n)/rad
    write(1,'(2f9.4,5f11.6)') x, y, TF(1,n), TF(5,n), TF(10,n), TF(17,n), TF(22,n)
    write(3,'(2f9.4,e16.6)') x, y, vorticity_2D(n)
    write(4,'(2f9.4,5e16.6)') x, y, vorticity_3D(1,n), vorticity_3D(5,n), vorticity_3D(10,n), &
    vorticity_3D(17,n), vorticity_3D(22,n)
    write(5,'(2f9.4,f12.5)') x, y, eta_n(n)
  enddo
 DO elem=1, elem2D
  elnodes=elem2D_nodes(:,elem)
  x = sum(w_cv(1:4,elem)*coord_nod2d(1,elnodes))/rad
  y = sum(w_cv(1:4,elem)*coord_nod2d(2,elnodes))/rad
write(2,'(2f9.4,4f10.4)') x,y,U_n(1,elem), V_n(1,elem), U_n(10,elem), V_n(10,elem)
enddo

  do n = 1,nod2D
    x=coord_nod2D(1, n)/rad
    y=coord_nod2D(2, n)/rad
    do nz=1,nsigma-1
  !  if (x == 20.0) write(9,'(2f11.4,f11.6)') y, -Z(nz,n)/50.0_WP,TF(nz,n)
   enddo

  enddo

!  DO elem=1, elem2D
!  elnodes=elem2D_nodes(:,elem)
!  x = sum(w_cv(1:4,elem)*coord_nod2d(1,elnodes))/rad
!  y = sum(w_cv(1:4,elem)*coord_nod2d(2,elnodes))/rad
!  do nz=1,nsigma-1
!  xy = sum(w_cv(1:4,elem)*Z(nz,elnodes))
!if (x>20. .and. x<20.06) write(8,'(2f11.4,4e16.6)') y,-xy/50.0_WP,Bar_pru_3D(nz,elem)/elem_area(elem), &
!       Bar_prv_3D(nz,elem)/elem_area(elem),x
!       enddo
!enddo

  DO n=1, nod2D
  x = coord_nod2D(1,n)/rad
  y = coord_nod2D(2,n)/rad
! write(8,'(2f11.4,e16.6)') x,y,hpre_2D(n)/density_0
enddo

  close(1)
  close(2)
  close(3)
  close(4)
  close(5)
  close(9)
  close(8)

end subroutine cross_sec
!===========================================================================================
subroutine output_vert_vel
  use o_MESH
  use o_ARRAYS
  use o_PARAM
  !
  implicit none
  integer                                :: n
  real(kind=WP)                      :: x, y
  character*22                      :: cross_Wvel='output_W\cross0000.dat'

  call vert_vel_cart

 write(cross_Wvel(15:18),'(I4.4)')  idnint(time/3600.0_WP/24.0_WP/1.0_WP)
 open(3,file=cross_Wvel)
 DO n=1, nod2D
    y=coord_nod2D(2, n)/rad
    x=coord_nod2D(1, n)/rad
write(3,'(2f12.4,1x,2e16.8)') x, y, W_n(3,n), W_n(12,n)
enddo

close(3)

end subroutine output_vert_vel
!========================================================================================
subroutine cross_sec_LE

  use o_MESH
  use o_ARRAYS
  use o_PARAM

  implicit none

  integer                            :: i, n, nz
  real(kind=WP)                       :: x,y,x1,x2,y1,y2
  character*23                       :: cross_file='output_TF/cross0000.dat'

  write(cross_file(16:19),'(I4.4)') idnint(time/60.0_WP/10.0_WP)
  y1=-1.E-4*rad
  y2=1.E-4*rad
  open(1,file=cross_file)
  do n = 1,nod2D
     y=coord_nod2D(2, n)
     x=coord_nod2D(1, n)
     if( (y > y1).AND.(y < y2) ) then
        do nz=1,nsigma-1
           write(1,'(f15.3,f10.4, f11.6)') x*r_earth/1000.0_WP, -Z(nz,n), TF(nz,n)
        enddo
     endif
  enddo

  close(1)


end subroutine cross_sec_LE
!========================================================================

! ====================================================================
!subroutine write_restart(n)
!  use o_ARRAYS
!  implicit none
!
!  integer  :: n
!  character(len = 256) :: fname,tmp ! file name
!
!  write(*,*) 'Write restart', n
!  !
!  ! Write restart
!  write(tmp,*) n
!  write(fname,*) 'restart_',trim(ADJUSTL(trim(tmp))),'.out'
!  fname=ADJUSTL(trim(fname))
!  write(*,*) 'Write restart file: ', fname
!  open(25,file=fname, form='unformatted')
!!
!! barotropic part
!!
!     write(25) eta_n, eta_n_1, eta_n_2, U_n_2D, U_n_1, U_n_2
!!
!! 3D velocity
!!
!     write(25) U_n, V_n, U_rhsAB, V_rhsAB, Wvel
!!
!! tracer
!!
!     write(25) TF, SF, T_old, S_old
!!
!! vertical mixing
!!
!     if (ver_mix == 2) write(25) tke, Av_node, teps, Kv
!     if (ver_mix == 3) write(25) bt, snu, KV
!!
!! time interval (+1)
!!
!     write(25) n
!     write(25) time_jd, time_jd0, dt, dt_2D
!!
!! T_counter (for tidal output)
!!
!     write(25) T_counter
!  close(25)
!end subroutine write_restart
