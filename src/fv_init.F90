!!$! Shallow water code based cell-vertex finite volume discretization 
!!$! AB3-AM4 time stepping as suggested in Shchepetkin and McWilliams
!!$! (2005)
!!$! Serial version
!!$! July 2012
!!$! sergey.danilov@awi.de 
!!$
!!$
!!$
!===========================================================================


!===========================================================================
subroutine read_mesh2D
USE o_MESH
USE o_PARAM
USE o_UTILIT
!
IMPLICIT NONE
!
INTEGER               :: nq, nt, n_quad, n_tri
INTEGER               :: n1,n2,n3, n4, nod(4)
INTEGER               :: n,ind, idx, index
REAL(kind=WP)         :: x1, x2
INTEGER, allocatable  :: elem_data(:)
INTEGER               :: i_error

! Requires new mesh format: elem2D lists 4 nodes in each row.
! For triangles the forth node coincides with the first one or just list the numbers of 3 nodes

  open(20,file=trim(MeshPath)//trim(TITLE)//'_nod2d.out', status='old')
  open(21,file=trim(MeshPath)//trim(TITLE)//'_elem2d.out', status='old')
  
  READ(20,*) nod2D
  ALLOCATE(coord_nod2D(2,nod2D),index_nod2D(nod2D))
  nobn=0
  
  if (cartesian) then 
    do n=1,nod2D
     read(20,*) nq, x1, x2, ind
     index_nod2D(nq) = ind
	 coord_nod2D(1,nq)=x1/r_earth
     coord_nod2D(2,nq)=x2/r_earth
	 if (ind==2) then
	 nobn=nobn+1
	 endif
  end do
  else
     do n=1,nod2D
     read(20,*) nq, x1, x2, ind
     index_nod2D(nq) = ind
     coord_nod2D(1,nq)=x1*rad
     coord_nod2D(2,nq)=x2*rad
	 if (ind==2) then
	 nobn=nobn+1
	 endif
  end do
  endif
  
  close(20)

  read(21,*)  elem2D
  ALLOCATE(elem2D_nodes(4,elem2D))
  ALLOCATE(elem_data(4*elem2D))
  elem_data(:)=-1

  ! meshes with quads have 4 columns, but TsunAWI grids may be
  ! purely triangular, with 3 columns each. Test, how many
  ! columns there are!
  
  read(21,*,iostat=i_error) elem_data(1:4*elem2D)
  write(*,*), i_error
  if (i_error == 0) then
     ! There is a fourth column => quad or mixed mesh
!$OMP PARALLEL DO PRIVATE(n)     
     do n=1,elem2D
        elem2D_nodes(1:4,n) = elem_data((n-1)*4+1:n*4)
     enddo
!$OMP END PARALLEL DO
     
  else
     ! No third column => triangles only

!$OMP PARALLEL DO PRIVATE(n)     
     do n=1,elem2D
        elem2D_nodes(1:3,n) = elem_data((n-1)*3+1:n*3)
        elem2D_nodes(4,n) = elem_data((n-1)*3+1)
     enddo
!$OMP END PARALLEL DO
  end if
     
  deallocate(elem_data)


  write(*,*) 'Mesh is read     ', 'nod2D=', nod2D,' elem2D=', elem2D
  write(*,*) 'Amount of open boundary nodes ', nobn


END SUBROUTINE read_mesh2D


!==========================================================================
program main

USE o_MESH
USE g_PARSUP

IMPLICIT NONE

integer           :: n, i, i_ok
character(len=10) :: arg1, arg2, string_npes=""
character(len=200) :: command
logical           :: l_ok

!$ real(kind=WP) :: t0, t1, t2, t3, t4, t5


! Determine number of partitions
! (argument to this program)
npes = -1
l_ok = (iargc() == 2)

if (l_ok) then

   call getarg(1,arg1)
   l_ok = (trim(arg1)=="-np")

   if (l_ok) then
      call getarg(2,arg2)
      read(arg2,*,iostat=i_ok) npes
      l_ok = (i_ok==0)
      l_ok = l_ok .and. (npes>1)
   endif
endif

if (.not. l_ok) then
   print *,'Usage:'
      print *,' ./fv_init.x -np <number of MPI partitions>'
   stop
end if


! We need mesh path, title, cyclic length, grid type (cartesian true/false)
 call READ_DATA_INIT

print *,'Distributing the mesh in ',trim(MeshPath),' into ',npes,' partitions'

!==================================================
! mkdir: directory for partition information as subdirectory of MeshPath

write(string_npes,'(i0)') npes
DistPath = trim(MeshPath)//"dist_"//trim(adjustl(string_npes))//"/"

command="mkdir "//trim(DistPath)
call system(trim(command))
! Fortran 2008 standard. Does not work with ifort 14, should with ifort 15:
! call execute_command_line(trim(command), exitstat=i_ok)
! if (i_ok/=0) then
!   print *,'The directory ',trim(DistPath),' could not be created sucessfully'
!   stop
! endif

call read_mesh2D 

call find_edges

call find_elem_neighbors

call find_mesh_partition
call save_partitioning

allocate(myList_nod2D(nod2D))
allocate(myList_elem2D(elem2D))

do mype = 0, npes-1
   call communication_nodn
   call communication_elemn
   
   call mymesh
   call save_dist_mesh
! call communication_edgen
end do

deallocate(myList_nod2D)
deallocate(myList_elem2D)

deallocate(part)
end program main

