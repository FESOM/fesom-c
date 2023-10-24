subroutine find_mesh_partition

  USE iso_c_binding, only: C_INT, idx_t=>C_INT32_T  ! Check in metis.h if metis is compiled 
                                                    ! with 32bit (default) or 64bit integer
                                                    ! Don't forget the interface below.
  USE o_MESH
  USE g_PARSUP
  
  implicit none
  
  integer, parameter :: MAX_ADJACENT=32
  integer            :: check_max_adjacent

  integer(idx_t) :: rowptr(0:nod2D)
!!$  integer(idx_t) :: colind(0:2*edge2D-1)
  integer(idx_t), allocatable :: colind(:)
  integer(idx_t) :: nod2D_c, npes_c, part_c(nod2D), ec_c
  integer(C_INT) :: ierr_c
  integer        :: num_ne(0:nod2D-1), ne(MAX_ADJACENT,0:nod2D-1)
  integer        :: nodes_per_partition(0:npes-1)
  real           ::  nod2D_inv

  integer :: ed, n, is, ie, a,b,c,d, q, el, nod(4), i, j

  interface
     subroutine  metis_wrapper(n, ptr, adj, np, part, ec, ierr) bind(C)

       USE iso_c_binding, only: C_INT, idx_t=>C_INT32_T  ! Check in metis.h if metis is compiled 
                                                         ! with 32bit (default) or 64bit integer
       USE o_MESH, only: nod2D, edge2D

       integer(idx_t), intent(in)  :: n, ptr(0:nod2d), adj(0:2*edge2D-1), np
       integer(idx_t), intent(out) :: part(nod2D), ec
       integer(C_INT), intent(out) :: ierr

     end subroutine metis_wrapper
  end interface

! build adjacency matrix 

  num_ne(:) = 0
  ne(:,:)   = -1
  check_max_adjacent = 0

  do el=1,elem2D

! all nodes in an element are adjacent in the sense of being halo nodes
! (also the opposite nodes of a quad: no edge, but the indirect connection
!  should be taken into account by metis domain decomposition)

     nod(1:4) = elem2D_nodes(1:4,el)-1  ! C-numbering
     q=4
     if (nod(1) == nod(4)) q=3  ! triangle

     do i=2,q
        do j=1,i-1
           if (all(ne(:,nod(i)) /= nod(j))) then

              num_ne(nod(i)) = num_ne(nod(i)) + 1
              num_ne(nod(j)) = num_ne(nod(j)) + 1
              
              check_max_adjacent = max(check_max_adjacent, num_ne(nod(i)), num_ne(nod(j)))

              if (check_max_adjacent <= MAX_ADJACENT ) then
                 ne(num_ne(nod(i)), nod(i)) = nod(j)
                 ne(num_ne(nod(j)), nod(j)) = nod(i)
              endif
           endif
                      
        end do
     end do
     
     if (check_max_adjacent > MAX_ADJACENT ) then
        print *,'Parameter in fv_distribute_mesh.f90,  find_mesh_partition, too small.'
        print *,'Recompile with larger value for MAX_ADJACENT.'
        stop
     endif
  end do

  print *,'Maximum number of adjacent (halo) nodes for one node: ',check_max_adjacent

! copy adjacency matrix to CSR-format

  rowptr(0) = 0   ! C-numbering
  do n=0,nod2D-1
     rowptr(n+1) = rowptr(n) + num_ne(n)
  enddo

  allocate(colind(0:rowptr(nod2D)-1))
  colind(:) = -1
  do n=0,nod2D-1
     is = rowptr(n)
     ie = rowptr(n+1) - 1
     colind(is:ie) = ne(1:num_ne(n),n)
  enddo
  

! Make sure we have the right type when calling Metis
  npes_c  = int(npes, idx_t)
  nod2D_c = int(nod2D, idx_t)


 call metis_wrapper(nod2D_c, rowptr, colind, npes_c, part_c, ec_c, ierr_c)
  
  if (ierr_c == 1) then
     print *,'Metis finished sucessfully with edge cut ',ec_c
  else
     print *,'Metis finished with error code ',ierr_c
     stop
  end if
     
  allocate(part(nod2D))
  nodes_per_partition(:)=0
  do n=1,nod2D
     part(n) = int(part_c(n))
     nodes_per_partition(part(n)) = nodes_per_partition(part(n)) +1
  enddo

  print *,'==='
  print *,'Load balancing: nodes per partition'
  nod2D_inv = 1./real(nod2D)
  do n=0,npes-1
     print *,'partition ',n,':',nodes_per_partition(n),',',100.*real(nodes_per_partition(n))*nod2D_inv,'%'
  enddo
  print *,'Mininum is ',minval(nodes_per_partition),', Maximum ',maxval(nodes_per_partition),&
          ', Inbalance ', 100.*real(maxval(nodes_per_partition)-minval(nodes_per_partition)) &
                             / real(minval(nodes_per_partition)),'%'

  open(22,file=trim(DistPath)//'rpart.out')
  do n=1,nod2D
     write(22,*) part(n)
  enddo
  close(22)

end subroutine find_mesh_partition
