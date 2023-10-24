SUBROUTINE test_elem

  USE o_MESH
  use g_parsup

  IMPLICIT NONE

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Check the order of nodes in elements; correct it if necessary to make
  ! it same sense (clockwise)
  ! Requirement: vertices of an element should be listed so that they
  ! form a full circle (trivial for triangles)
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(kind=WP)   ::  a(2), b(2), c(2), d(2),  r, r1
  integer             ::  n, nx, elnodes(4)
  real(kind=WP)   :: t0, t1 ! to measure time if desired in MPI

  DO n=1, myDim_elem2D

     elnodes=elem2D_nodes(:,n)

     if(elnodes(1)==elnodes(4)) then	! the case of triangle

        a=coord_nod2D(:,elnodes(1))
        b=coord_nod2D(:,elnodes(2))-a
        c=coord_nod2D(:,elnodes(3))-a

        if(b(1)>cyclic_length/2.0_WP) b(1)=b(1)-cyclic_length
        if(b(1)<-cyclic_length/2.0_WP) b(1)=b(1)+cyclic_length
        if(c(1)>cyclic_length/2.0_WP) c(1)=c(1)-cyclic_length
        if(c(1)<-cyclic_length/2.0_WP) c(1)=c(1)+cyclic_length

        r=b(1)*c(2)-b(2)*c(1)
        if (r>0.0_WP) then
           ! Vector b is to right of c
           ! Exchange second and third nodes:
           write(*,*) 'Vertices exchanged for triangle ', n
           nx=elnodes(2)
           elnodes(2)=elnodes(3)
           elnodes(3)=nx
           elem2D_nodes(:,n)=elnodes
        end if

     else                             ! the case of a quadrilateral

        a=coord_nod2D(:,elnodes(1))
        b=coord_nod2D(:,elnodes(2))-a
        c=coord_nod2D(:,elnodes(3))-a
        d=coord_nod2D(:,elnodes(4))-a


        if(b(1)>cyclic_length/2.0_WP) b(1)=b(1)-cyclic_length
        if(b(1)<-cyclic_length/2.0_WP) b(1)=b(1)+cyclic_length
        if(c(1)>cyclic_length/2.0_WP) c(1)=c(1)-cyclic_length
        if(c(1)<-cyclic_length/2.0_WP) c(1)=c(1)+cyclic_length
        if(d(1)>cyclic_length/2.0_WP) d(1)=d(1)-cyclic_length
        if(d(1)<-cyclic_length/2.0_WP) d(1)=d(1)+cyclic_length

        r=b(1)*c(2)-b(2)*c(1)
        r1=b(1)*d(2)-b(2)*d(1)
        if(r1*r<0.0_WP) then            ! b is to right of one, and to left of
                                        ! the other one.

           write(*,*) 'Node list problem ', n
           stop
        else
           if(r>0) then ! vector b is to right of c (and of d)
                        ! wrong rotation sense. We need to replace
                        ! b and d:
              write(*,*) 'Vertices exchanged for quad ', n
              nx=elnodes(2)
              elnodes(2)=elnodes(4)
              elnodes(4)=nx
              elem2D_nodes(:,n)=elnodes
           end if
        end if
     end if
  END DO


END SUBROUTINE  test_elem

!=========================================================================

SUBROUTINE load_edges

  USE o_MESH
  USE o_PARAM
  USE g_PARSUP

  IMPLICIT NONE

  character*1000                        :: file_name
  integer                               :: counter, n, m, nn, k, q, fileID, q1
  integer                               :: elems(2), elem
  integer                               :: elnodes(4), ed(2), eledges(4), iaux4(4)
  integer, allocatable                  :: aux(:)
  real(kind=WP)                         :: t0, t1
  integer                               :: nchunk, chunk_size, ipos, iofs, mesh_check
  integer, allocatable, dimension(:)    :: mapping
  integer, allocatable, dimension(:,:)  :: ibuff
  integer                               :: ierror              ! return error code

  t0=MPI_Wtime()

  !==============================
  ! Edge array is already available (we computed it in the init phase)
  ! (a) Read list of edges and tri containing them from file
  !==============================
  !mesh related files will be read in chunks of chunk_size
  chunk_size=100000
  !==============================
  ! Allocate mapping array (chunk_size)
  ! It will be used for several purposes
  !==============================
  allocate(mapping(chunk_size))
  allocate(ibuff(4,chunk_size))

  ! 0 proc reads the data and sends it to other procs
  fileID=10
  if (mype==0) then
     file_name=trim(meshpath)//'edgenum.out'
     write(*,*) 'reading '// trim(file_name)
     open(fileID, file=trim(file_name))
     read(fileID,*) edge2D
     read(fileID,*) edge2D_in
     write(*,*) '2D mesh info : edge2D=', edge2D
     close(fileID)
  end if

  call MPI_BCast(edge2D,    1, MPI_INTEGER, 0, MPI_COMM_FESOM_C, ierror)
  call MPI_BCast(edge2D_in, 1, MPI_INTEGER, 0, MPI_COMM_FESOM_C, ierror)

  allocate(edge_nodes(2,myDim_edge2D+eDim_edge2D))
  allocate(edge_tri(2,myDim_edge2D+eDim_edge2D))

  ! 0 proc reads the data in chunks and distributes it between other procs
  if (mype==0) then
     file_name=trim(meshpath)//'edges.out'
     open(fileID,   file=trim(file_name))
     write(*,*) 'reading '// trim(file_name)

     file_name=trim(meshpath)//'edge_tri.out'
     open(fileID+1, file=trim(file_name))
     write(*,*) 'reading '// trim(file_name)
  end if

  mesh_check=0
  do nchunk=0, (edge2D-1)/chunk_size
     !create the mapping for the current chunk
     mapping(1:chunk_size)=0
     do n=1, myDim_edge2D+eDim_edge2D
        ipos=(myList_edge2D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_edge2D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do
     !read the chunk into the buffers
     k=min(chunk_size, edge2D-nchunk*chunk_size)
     if (mype==0) then
        do n=1, k
           read(fileID  ,*) ibuff(1:2,n) !edge nodes
           read(fileID+1,*) ibuff(3:4,n) !edge elements
        end do
     end if
     call MPI_BCast(ibuff(1:4,1:k), 4*k, MPI_INTEGER, 0, MPI_COMM_FESOM_C, ierror)
     !    call MPI_BCast(ibuff(1:k,2), k, MPI_INTEGER, 0, MPI_COMM_FESOM_C, ierror)
     !    call MPI_BCast(ibuff(1:k,3), k, MPI_INTEGER, 0, MPI_COMM_FESOM_C, ierror)
     !    call MPI_BCast(ibuff(1:k,4), k, MPI_INTEGER, 0, MPI_COMM_FESOM_C, ierror)
     ! fill the local arrays
     do n=1, k
        if (mapping(n)>0) then
           mesh_check=mesh_check+1
           edge_nodes   (:, mapping(n))=ibuff(1:2,n)
           edge_tri(:, mapping(n))=ibuff(3:4,n)
        end if
     end do

  end do

  if (mesh_check/=myDim_edge2D+eDim_edge2D) then
     write(*,*) 'ERROR while reading edges.out/edge_tri.out on mype=', mype
     write(*,*) mesh_check, ' values have been read in according to partitioning'
     write(*,*) 'it does not equal to myDim_edge2D+eDim_edge2D = ', myDim_edge2D+eDim_edge2D
     stop
  end if

  if (mype==0) then
     close(fileID)
     close(fileID+1)
  end if

  ! =========
  ! edges are in global numbering, transform it to local indexing
  ! =========
  mesh_check=0
  do nchunk=0, (nod2D-1)/chunk_size
     mapping(1:chunk_size)=0
     do n=1, myDim_nod2D+eDim_nod2D
        ipos=(myList_nod2D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_nod2D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do
     do n=1, myDim_edge2D+eDim_edge2D
        do m=1, 2
           nn=edge_nodes(m, n)
           ipos=(nn-1)/chunk_size
           if (ipos==nchunk) then
              mesh_check=mesh_check+1
              iofs=nn-nchunk*chunk_size
              ! minus sign is required to avoid modified entry being modified in another chunk
              ! will be changed to plus at the end
              edge_nodes(m,n)=-mapping(iofs)
           end if
        end do
     end do
  end do
  edge_nodes=-edge_nodes
  mesh_check=mesh_check/2
  if (mesh_check/=myDim_edge2D+eDim_edge2D) then
     write(*,*) 'ERROR while transforming edge nodes to local indexing on mype=', mype
     write(*,*) mesh_check, ' edges have been transformed!'
     write(*,*) 'It does not equal to myDim_edge2D+eDim_edge2D = ', myDim_edge2D+eDim_edge2D
     stop
  end if
  ! =========
  ! edge_tri are in global numbering, transform it to local indexing
  ! =========

  !*** Vadym pointed to the bug with bilding edge_tri
  ! if not doing so the land boundary elements will be set to 999 instead of being zeros or negative!
  where (edge_tri<0)
     edge_tri=0
  end where
  !***
  mesh_check=0
  do nchunk=0, (elem2D-1)/chunk_size
     mapping(1:chunk_size)=0
     do n=1, myDim_elem2D+eDim_elem2D+eXDim_elem2D
        ipos=(myList_elem2D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_elem2D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do
     do n=1, myDim_edge2D+eDim_edge2D
        do m=1, 2
           nn=edge_tri(m, n)
           ipos=(nn-1)/chunk_size
           if (ipos==nchunk .and. nn > 0) then
              mesh_check=mesh_check+abs(m-2) !only first triangle will contribute to statistic
              iofs=nn-nchunk*chunk_size
              ! minus sign is required to avoid modified entry being modified in another chunk
              ! will be changed to plus at the end
              edge_tri(m,n)=-mapping(iofs)
           end if
        end do
     end do
  end do
  edge_tri=-edge_tri

  if (mesh_check/=myDim_edge2D+eDim_edge2D) then
     write(*,*) 'ERROR while transforming edge elements to local indexing on mype=', mype
     write(*,*) mesh_check, ' edges have been transformed!'
     write(*,*) 'It does not equal to myDim_edge2D+eDim_edge2D = ', myDim_edge2D+eDim_edge2D
     stop
  end if

  !SH Re-establish negative indices for boundary edges
  where (edge_tri<=0)
     edge_tri=-999
  end where

  ! =========


  ! Now the only way to check whether an edge is on boundary is
  ! through myList_edge2D(n):  myList_edge2D(n)>edge2D_in == boundary edge

  ! (b) We need an array inverse to edge_tri listing edges
  ! of a given element (quad/triangle)
  allocate(elem_edges(4,myDim_elem2D))
  elem_edges(:,:)=0

  allocate(aux(myDim_elem2D))
  aux=0
  DO n=1, myDim_edge2D+eDim_edge2D
     DO k=1,2
        q=edge_tri(k,n)   ! triangle number
        if((q>0).and.(q<=myDim_elem2D)) then
           aux(q)=aux(q)+1
           elem_edges(aux(q),q)=n
        end if
     END DO
  END DO
  deallocate(aux)

if (mype==0) then
! print *,'PE=0 - myDim_elem2D,eDim_elem2D',myDim_elem2D,eDim_elem2D
 !do n=1,myDim_elem2D
 ! print *,n,elem_edges(:,n)
 !end do
! print *,'========================'
! print *,elem_edges(:,16)
! print *,elem2D_nodes(:,16)
! print *,edge_nodes(:,elem_edges(1,16))
! print *,edge_nodes(:,elem_edges(2,16))
! print *,edge_nodes(:,elem_edges(3,16))
! print *,edge_nodes(:,elem_edges(4,16))

end if


!SH Skipped the reordering for now since with quads it is
!SH not unique!!!

!SH Caution!!! took the part from find_edges
  DO elem=1,myDim_elem2D
     elnodes=elem2D_nodes(:,elem)
     q1=4
     if(elnodes(1)==elnodes(4)) q1=3
     eledges=elem_edges(:,elem)
     DO q=1,q1-1
        DO k=1,q1
           if(((edge_nodes(1,eledges(k))==elnodes(q)).and. &
                (edge_nodes(2,eledges(k))==elnodes(q+1))).or. &
                ((edge_nodes(1,eledges(k))==elnodes(q+1)).and. &
                (edge_nodes(2,eledges(k))==elnodes(q)))) then
              elem_edges(q,elem)=eledges(k)
              exit
           end if
        END DO
     END DO
     DO k=1,q1
        if(((edge_nodes(1,eledges(k))==elnodes(q1)).and. &
             (edge_nodes(2,eledges(k))==elnodes(1))).or. &
             ((edge_nodes(1,eledges(k))==elnodes(1)).and. &
             (edge_nodes(2,eledges(k))==elnodes(q1)))) then
           elem_edges(q1,elem)=eledges(k)
           exit
        end if
     END DO
     if(q1==3) elem_edges(4,elem)=elem_edges(1,elem)
  END DO


!!$  ! The edges in this list should be ordered so that they
!!$  ! are listed in the same rotation sense as nodes.
!!$  DO elem=1,myDim_elem2D
!!$     iaux4(:)=0
!!$     elnodes=elem2D_nodes(:,elem)
!!$     eledges=elem_edges(:,elem)
!!$     DO q=1,4
!!$        DO k=1,4
!!$           if((edge_nodes(1,eledges(k)).ne.elnodes(q)).and. &
!!$              (edge_nodes(2,eledges(k)).ne.elnodes(q)) .and. &
!!$               ALL(iaux4/=eledges(k)) ) then
!!$              !elem_edges(q,elem)=eledges(k)
!!$              iaux4(q)=eledges(k)
!!$              exit
!!$           end if
!!$        END DO
!!$     END DO
!!$     elem_edges(:,elem)=iaux4(:)
!!$  END DO
!!$  ! The edge and elem lists agree in the sense that edge1 does not
!!$  ! contain node 1 and so on

  deallocate(ibuff)
  deallocate(mapping)

if (mype==0) then
 !print *,myDim_elem2D,eDim_elem2D
 !do n=1,myDim_elem2D
 ! print *,n,elem_edges(:,n)
 !end do

end if


  t1=MPI_Wtime()
  if (mype==0) then
     write(*,*) 'load_edges finished in ', t1-t0, ' seconds'
     write(*,*) '========================='
  endif

END SUBROUTINE load_edges

!===========================================================================

SUBROUTINE find_edges

  USE o_MESH
  USE o_PARAM
  IMPLICIT NONE

  integer, allocatable, dimension(:)    :: aux1
  integer                               :: counter, counter_in, n, k, q
  integer                               :: elem, elem1, elems(2), q1, q2
  integer                               :: elnodes(4), ed(2),flag, eledges(4)
  integer                               :: temp(100), node
  real(kind=WP)                    :: xc(2), xe(2), ax(3), amin

  ! ====================
  ! (a) find edges. To make the procedure fast
  ! one needs neighbourhood arrays
  ! ====================

  allocate(ne_num(nod2d))
  ne_num=0
  DO n=1,elem2D
     elnodes=elem2D_nodes(:,n)
     if(elnodes(1)==elnodes(4)) then
        ne_num(elnodes(1:3))=ne_num(elnodes(1:3))+1
     else
        ne_num(elnodes)=ne_num(elnodes)+1
     end if
  END DO
  k=maxval(ne_num)               ! maximum number of neighbour elements

  allocate(ne_pos(k, nod2D),nn_num(nod2D))
  ne_num=0
  DO n=1,elem2D
     elnodes=elem2D_nodes(:,n)
     q1=4
     if(elnodes(1)==elnodes(4)) q1=3
     DO q=1,q1
        ne_num(elnodes(q))=ne_num(elnodes(q))+1
        ne_pos(ne_num(elnodes(q)),elnodes(q))=n
     END DO
  END DO			        ! neighbor elements are found

  ! count neighbour nodes
  ! In quads we should count the nodes that are
  ! connected by edges!
  allocate(aux1(nod2D))
  aux1=0

  DO n=1, nod2D
     counter=0
     DO k=1, ne_num(n)
        elem=ne_pos(k,n)
        elnodes=elem2D_nodes(:,elem)

        if(elnodes(1)==elnodes(4)) then
           DO q=1,3
              if(elnodes(q)==n) CYCLE
              if(aux1(elnodes(q)).ne.1) then
                 counter=counter+1
                 aux1(elnodes(q))=1
                 temp(counter)=elnodes(q)
              end if
           END DO
        else
           ! Find the position of n in elnodes:
           if(elnodes(1)==n) then
              ed(1)=elnodes(2)
              ed(2)=elnodes(4)
           end if
           if(elnodes(2)==n) then
              ed(1)=elnodes(1)
              ed(2)=elnodes(3)
           end if
           if(elnodes(3)==n) then
              ed(1)=elnodes(2)
              ed(2)=elnodes(4)
           end if
           if(elnodes(4)==n) then
              ed(1)=elnodes(1)
              ed(2)=elnodes(3)
           end if
           DO q=1,2
              if(aux1(ed(q)).ne.1) then
                 counter=counter+1
                 aux1(ed(q))=1
                 temp(counter)=ed(q)
              end if
           END DO
        end if
     END DO
     nn_num(n)=counter
     aux1(temp(1:counter))=0
  END DO
  deallocate(aux1)

  allocate(nn_pos(maxval(nn_num)+1,nod2D))
  allocate(aux1(nod2D))
  aux1=0

  DO n=1, nod2D
     counter=0
     DO k=1, ne_num(n)
        elem=ne_pos(k,n)
        elnodes=elem2D_nodes(:,elem)

        if(elnodes(1)==elnodes(4)) then
           DO q=1,3
              if(elnodes(q)==n) CYCLE
              if(aux1(elnodes(q)).ne.1) then
                 counter=counter+1
                 aux1(elnodes(q))=1
                 temp(counter)=elnodes(q)
              end if
           END DO
        else
           ! Find the position of n in elnodes:
           if(elnodes(1)==n) then
              ed(1)=elnodes(2)
              ed(2)=elnodes(4)
           end if
           if(elnodes(2)==n) then
              ed(1)=elnodes(1)
              ed(2)=elnodes(3)
           end if
           if(elnodes(3)==n) then
              ed(1)=elnodes(2)
              ed(2)=elnodes(4)
           end if
           if(elnodes(4)==n) then
              ed(1)=elnodes(1)
              ed(2)=elnodes(3)
           end if
           DO q=1,2
              if(aux1(ed(q)).ne.1) then
                 counter=counter+1
                 aux1(ed(q))=1
                 temp(counter)=ed(q)
              end if
           END DO
        end if
     END DO
     nn_num(n)=counter+1
     aux1(temp(1:counter))=0
     nn_pos(2:counter+1,n)=temp(1:counter)
     nn_pos(1,n)=n

  END DO
  deallocate(aux1)
  ! neighboring nodes are found. First in the list is the node itself

  ! ====================
  ! (b) Find edges and elements containing them.
  ! ====================
  counter=0
  ! Count edges:
  DO n=1,nod2D
     ! ====================
     ! form edges with n by cycling over neighboring
     ! nodes (if edges are not accounted yet).
     ! New edges are added only if neighbor>n
     ! ====================
     DO q=2,nn_num(n)
        node=nn_pos(q,n)
        if(node<n) CYCLE
        counter=counter+1   ! new edge (n,node)
     END DO
  END DO
  edge2D=counter

  allocate(edge_nodes(2,counter), edge_tri(2, counter))
  counter=0
  counter_in=0
  DO n=1,nod2D
     DO q=2,nn_num(n)
        node=nn_pos(q,n)
        if(node<n) CYCLE
        counter=counter+1   ! new edge (n,node)
        ! find elements containing n and node
        flag=0
        DO k=1, ne_num(n)
           elem=ne_pos(k,n)
           elnodes=elem2D_nodes(:,elem)
           q2=4
           if(elnodes(1)==elnodes(4)) q2=3
           DO q1=1,q2
              if (elnodes(q1)==node) then
                 flag=flag+1
                 elems(flag)=elem
                 EXIT
              end if
           END DO
        END DO
        if(flag==2) then
           counter_in=counter_in+1
           edge_nodes(1,counter_in)=n
           edge_nodes(2,counter_in)=node
           edge_tri(:,counter_in)=elems
        end if
     END DO
  END DO
  edge2D_in=counter_in

  ! Repeate to collect boundary edges:
  counter=0
  DO n=1,nod2D
     DO q=2,nn_num(n)
        node=nn_pos(q,n)
        if(node<n) CYCLE
        counter=counter+1   ! new edge (n,node)
        ! find triangles containing n and node
        flag=0
        DO k=1, ne_num(n)
           elem=ne_pos(k,n)
           elnodes=elem2D_nodes(:,elem)
           q2=4
           if(elnodes(1)==elnodes(4)) q2=3
           DO q1=1,q2
	      if (elnodes(q1)==node) then
                 flag=flag+1
                 elems(flag)=elem
                 EXIT
	      end if
           END DO
        END DO
        if(flag==1) then
           counter_in=counter_in+1
           edge_nodes(1,counter_in)=n
           edge_nodes(2,counter_in)=node
           elems(2)=-999
           edge_tri(:,counter_in)=elems
        end if
     END DO
  END DO

  ! Edges from edge2D_in+1 to edge2D lie on the lateral boundary
  ! The rest (1:edge2D_in) are internal edges

  ! ====================
  ! (d) the list of elements on both sides of edge e, edge_tri(:,e)
  ! should be ordered so that the first is the element which is to the left of
  ! the edge (vector pointing from the first to the second node of the edge.
  ! If the edge is on the boundary, there is only the first element.
  ! ====================

  DO n=1, edge2D
     ed=edge_nodes(:,n)
     elem=edge_tri(1,n)
     call elem_center(elem, xc(1), xc(2), amin)
     xc=xc-coord_nod2D(:,ed(1))
     xe=coord_nod2D(:,ed(2))-coord_nod2D(:,ed(1))
     if(xe(1)>cyclic_length/2.0_WP) xe(1)=xe(1)-cyclic_length
     if(xe(1)<-cyclic_length/2.0_WP) xe(1)=xe(1)+cyclic_length
     if(xc(1)>cyclic_length/2.0_WP) xc(1)=xc(1)-cyclic_length
     if(xc(1)<-cyclic_length/2.0_WP) xc(1)=xc(1)+cyclic_length
     if(xc(1)*xe(2)-xc(2)*xe(1)>0.0_WP) then
        ! Vector drawn to the center of the first triangle is to the right
        ! of the edge vector. Triangles have to be exchanged:

        elem=edge_tri(1,n)
        elem1=edge_tri(2,n)
        if(elem1>0) then
	   edge_tri(1,n)=elem1
	   edge_tri(2,n)=elem
        else
	   elem=edge_nodes(2,n)         ! change the order of nodes
	   edge_nodes(2,n)=edge_nodes(1,n)
	   edge_nodes(1,n)=elem
        end if
     end if
  END DO


  ! ====================
  ! (e) We need an array inverse to edge_tri listing edges
  ! of a given element
  ! ====================
  allocate(elem_edges(4,elem2D))
  elem_edges=0
  allocate(aux1(elem2D))
  aux1=0

  DO n=1, edge2D
     DO k=1,2
        q=edge_tri(k,n)   ! element number
        if (q>0) then
           aux1(q)=aux1(q)+1
           elem_edges(aux1(q),q)=n
        end if
     END DO
  END DO
  deallocate(aux1)

  ! We order the edges in this list so that they
  ! are listed in the same rotation sense as nodes.
  ! First is the edge formed by elnodes(1:2), and so on
  DO elem=1,elem2D
     elnodes=elem2D_nodes(:,elem)
     q1=4
     if(elnodes(1)==elnodes(4)) q1=3
     eledges=elem_edges(:,elem)
     DO q=1,q1-1
        DO k=1,q1
           if(((edge_nodes(1,eledges(k))==elnodes(q)).and. &
                (edge_nodes(2,eledges(k))==elnodes(q+1))).or. &
                ((edge_nodes(1,eledges(k))==elnodes(q+1)).and. &
                (edge_nodes(2,eledges(k))==elnodes(q)))) then
              elem_edges(q,elem)=eledges(k)
              exit
           end if
        END DO
     END DO
     DO k=1,q1
        if(((edge_nodes(1,eledges(k))==elnodes(q1)).and. &
             (edge_nodes(2,eledges(k))==elnodes(1))).or. &
             ((edge_nodes(1,eledges(k))==elnodes(1)).and. &
             (edge_nodes(2,eledges(k))==elnodes(q1)))) then
           elem_edges(q1,elem)=eledges(k)
	   exit
        end if
     END DO
     if(q1==3) elem_edges(4,elem)=elem_edges(1,elem)
  END DO

END SUBROUTINE find_edges
!===================================================================
subroutine edge_center(n1, n2, x, y)

  USE o_MESH
  USE o_PARAM

  ! Returns coordinates of edge center in x and y

  implicit none

  integer       :: n1, n2   ! nodes of the edge
  real(kind=WP) :: x, y, a(2), b(2)

  a=coord_nod2D(:,n1)
  b=coord_nod2D(:,n2)

  if(a(1)-b(1)>cyclic_length/2.0_WP) a(1)=a(1)-cyclic_length
  if(a(1)-b(1)<-cyclic_length/2.0_WP) b(1)=b(1)-cyclic_length
  x=0.5_WP*(a(1)+b(1))
  y=0.5_WP*(a(2)+b(2))

end subroutine edge_center
!====================================================================
subroutine elem_center(elem, x, y, s)

  ! Returns coordinates of elem center in x and y

  USE o_MESH
  USE o_PARAM

  implicit none

  integer       :: elem, elnodes(4), k
  real(kind=WP)  :: x, y, x1, y1, s1, s2, ax(4), ay(4),  amin, s

  elnodes=elem2D_nodes(:,elem)

  if (any(elnodes==0)) then
!     print *,elem,elnodes(:)
  end if

  if(elnodes(1)==elnodes(4)) then
     ax=coord_nod2D(1, elnodes)
     ay=coord_nod2D(2, elnodes)
     amin=minval(ax(1:3))
     DO k=1,3
        if(ax(k)-amin>cyclic_length/2.0) ax(k)=ax(k)-cyclic_length
     END DO
     x=sum(ax(1:3))/3.0_WP
     y=sum(ay(1:3))/3.0_WP
     s=0.5_WP*abs((ax(2)-ax(1))*(ay(3)-ay(1))-(ay(2)-ay(1))*(ax(3)-ax(1)))
  else
     ax=coord_nod2D(1, elnodes)
     ay=coord_nod2D(2, elnodes)
     amin=minval(ax)
     DO k=1,4
        if(ax(k)-amin>cyclic_length/2.0_WP) ax(k)=ax(k)-cyclic_length
     END DO
     x=ax(1)+ax(3)+ax(4)
     y=ay(1)+ay(3)+ay(4)
     x1=ax(1)+ax(2)+ax(3)
     y1=ay(1)+ay(2)+ay(3)
     ! areas
     s1=abs((ax(3)-ax(4))*(ay(1)-ay(4))-(ay(3)-ay(4))*(ax(1)-ax(4)))
     s2=abs((ax(3)-ax(2))*(ay(1)-ay(2))-(ay(3)-ay(2))*(ax(1)-ax(2)))
     x=(x*s1+x1*s2)/3.0_WP/(s1+s2)
     y=(y*s1+y1*s2)/3.0_WP/(s1+s2)
     s=0.5_WP*(s1+s2)
  end if

end subroutine elem_center
!=======================================================================
SUBROUTINE find_elem_neighbors

  ! For each element three or four its element neighbors are found
  USE o_PARAM
  USE o_MESH

  use g_parsup
  use g_comm_auto

  implicit none
  integer               :: elem, eledges(4), elem1, j, n, nd, q, elnodes(4), node, el
  integer, allocatable  :: temp_i(:)
  integer               :: mymax(npes), rmax(npes)


#ifdef USE_MPI
  CALL MPI_BARRIER(MPI_COMM_FESOM_C, MPIerr)
#endif

  allocate(elem_neighbors(4,myDim_elem2D))
  elem_neighbors=0

  DO elem=1,myDim_elem2D
     eledges=elem_edges(:,elem)
     q=4

     !if(eledges(1)==eledges(4)) q=3 !SH is it really true for edges or only for nodes???
     if(eledges(1)==eledges(4) .or. eledges(4)==0) q=3
     DO j=1,q
        elem1=edge_tri(1,eledges(j))
        if(elem1==elem) elem1=edge_tri(2,eledges(j))
        elem_neighbors(j,elem)=elem1
!!$#ifdef USE_MPI
!!$if (elem_neighbors(j,elem)<=0) print *,'LOOK: ',mylist_elem2D(elem),j,elem_neighbors(j,elem)
!!$#else
!!$if (elem_neighbors(j,elem)<=0) print *,'LOOK: ',elem,j,elem_neighbors(j,elem)
!!$#endif
     END DO
     if(q==3) then
        elem_neighbors(4,elem)=elem_neighbors(1,elem)
     end if

  END DO

  ! =============
  ! Node neighbourhood
  ! == elements that contain node n
  ! We need eDim neighborhood too for MUSCL advection.
  ! And we already have the place allocated for all
  ! these neighbor elements: it is eDim_elem2D+eXDim_elem2D
  ! =============

  allocate(nod_in_elem2D_num(myDim_nod2D+eDim_nod2D))
  nod_in_elem2D_num=0
  do el=1,myDim_elem2D
     q=4
     if ( elem2D_nodes(1,el) == elem2D_nodes(4,el)) q=3  ! triangle
     do j=1,q
        nd = elem2D_nodes(j,el)
        if (nd>myDim_nod2D) cycle
        nod_in_elem2D_num(nd)=nod_in_elem2D_num(nd)+1
     end do
  end do

#ifdef USE_MPI

  CALL MPI_BARRIER(MPI_COMM_FESOM_C, MPIerr)

  mymax=0
  rmax=0
  mymax(mype+1)=maxval(nod_in_elem2D_num(1:myDim_nod2D))
  call MPI_AllREDUCE( mymax, rmax, &
       npes, MPI_INTEGER,MPI_SUM, &
       MPI_COMM_FESOM_C, MPIerr)

  allocate(nod_in_elem2D(maxval(rmax),myDim_nod2D+eDim_nod2D))

#else

  allocate(nod_in_elem2D(maxval(nod_in_elem2D_num),nod2D))

#endif

  nod_in_elem2D=0
  nod_in_elem2D_num=0
  do el=1,myDim_elem2D
     q=4
     if ( elem2D_nodes(1,el) == elem2D_nodes(4,el)) q=3  ! triangle
     do j=1,q
        nd=elem2D_nodes(j,el)
        nod_in_elem2D_num(nd)=nod_in_elem2D_num(nd)+1
        nod_in_elem2D(nod_in_elem2D_num(nd),nd)=el
     end do
  end do

#ifdef USE_MPI

  call exchange_nod(nod_in_elem2D_num)

  allocate (temp_i(myDim_nod2D+eDim_nod2D))
  temp_i=0
  DO n=1, maxval(rmax)
     ! Exchange global element numbers
     do j=1,myDim_nod2D
        if (nod_in_elem2D(n,j)>0) temp_i(j)=myList_elem2D(nod_in_elem2D(n,j))
     enddo
     call exchange_nod(temp_i)
     nod_in_elem2D(n,:)=temp_i
  END DO
  deallocate(temp_i)
  ! Substitute back local element numbers
  allocate(temp_i(elem2D))
  temp_i=0
  Do n=1, myDim_elem2D+eDim_elem2D+eXDim_elem2D
     temp_i(myList_elem2D(n))=n
  END DO
  DO n=1, myDim_nod2D+eDim_nod2D
     DO j=1, nod_in_elem2D_num(n)
        nod_in_elem2D(j,n)=temp_i(nod_in_elem2D(j,n))
     END DO
  END DO

  deallocate(temp_i)

  ! Among elem_neighbors there can be negative numbers. These correspond to
  ! boundary elements for which neighbours are absent. However, an element
  ! should have at least two valid neighbors, otherwise least square
  ! interpolation for velocities will not work properly.

  ! The rotation sense: corresponds to edges, and edges correspond
  ! to nodes

  ! ============
  ! Test that there are at least two neighbors at the surface:
  ! ============

  DO elem=1,myDim_elem2D
     elem1=0
     DO j=1,4
        if(elem_neighbors(j,elem)>0) elem1=elem1+1
     END DO
     if (elem1<2) then
        write(*,*) 'Only one neighbor in elem (loc/glob) ', elem,myList_elem2D(elem)
        !call par_ex(1)
        !STOP
     end if
  END DO

#endif

  if (mype==0) then
     write(*,*) 'find_neighbors finished'
     write(*,*) '======================='
  endif

END SUBROUTINE find_elem_neighbors
!===================================================================
SUBROUTINE mesh_arrays1

  USE o_MESH
  USE o_PARAM

  use g_parsup
  use g_comm_auto

  IMPLICIT NONE

  ! Collects auxilliary information on the mesh
  ! Allocated and filled in are:
  ! elem_area(elem2D)
  ! area(nod2D)
  ! metric_factor(elem2D)
  ! elem_cos(elem2D)
  ! coriolis(elem2D)

  integer              :: j, elem, elnodes(6), n
  real(kind=WP)	     :: ax, ay, aa, rc(2), r1(2), r0(2), rr(2), rl(2)
  real(kind=WP)  :: min_area, max_area, min_scarea, max_scarea, min_w, max_w

  allocate(elem_area(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
  allocate(area(myDim_nod2D+eDim_nod2D))
  allocate(metric_factor(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
  allocate(elem_cos(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
  allocate(coriolis(myDim_elem2D))
  allocate(w_cv(4,myDim_elem2D+eDim_elem2D+eXDim_elem2D))

  elem_area = 0.0_WP
  elem_cos = 0.0_WP
  area = 0.0_WP
  metric_factor = 0.0_WP
  coriolis = 0.0_WP
  w_cv = 0.0_WP
  ! ============
  ! The areas of elements:
  ! Coriolis
  ! cosines and metric factor
  ! ============

  DO elem=1,myDim_elem2D
     call elem_center(elem,ax,ay,aa)
     !call r2g(lon, lat, ax, ay)   ! reserved for rotated meshes, see 3D model
     coriolis(elem)=2.0_WP*omega*sin(ay)
     elem_cos(elem)=cos(ay)
     metric_factor(elem)=tan(ay)/r_earth
     elem_area(elem)=aa
  END DO

  if (cartesian) then
     elem_cos(:)=1.0_WP
     metric_factor(:)=0.0_WP
  endif

  !VF The possibility is to turn on/off Coriolis, if it is nesessary.
  if (Coriolis_TF) then
     if (cartesian) then
        coriolis(:)=2.0_WP*omega*sin(lat_cartesian)
        !coriolis=2.0_WP*omega*0.71_WP
     endif
!     write(*,*) 'Coriolis max, min:',maxval(coriolis), minval(coriolis)
  else
     coriolis=0.0_WP
  endif


  elem_area=elem_area*elem_cos

#ifdef USE_MPI
  call exchange_elem(elem_area)
  call exchange_elem(metric_factor)
  call exchange_elem(elem_cos)
#endif


  ! =============
  ! Scalar control element and
  ! compute the area weight
  ! w_cv = 1/3 on triangles
  ! w_cv = Scv/Sc  where
  ! Sc - cell area
  ! Scv - the part of it in the scalar control volume
  !         around vertex
  ! 12.11.2014 Androsov Alexey
  ! =============

  area(:)   = 0.0_WP
  w_cv(:,:) = 0.0_WP

  DO elem=1, myDim_elem2D !+eDim_elem2D+eXDim_elem2D
     elnodes(2:5)=elem2D_nodes(:,elem)
     if(elnodes(2)==elnodes(5)) then
        ! Triangle
        DO j=2,4
           area(elnodes(j))=area(elnodes(j))+elem_area(elem)/3.0_WP
           w_cv(j-1,elem) = 1.0_WP/3.0_WP
        END DO
     else
        ! We have a quadrilateral
        call elem_center(elem, rc(1), rc(2), aa)
        elnodes(1)=elnodes(5)
        elnodes(6)=elnodes(2)

        DO j=2,5
           r1=coord_nod2D(:,elnodes(j-1))
           r0=coord_nod2D(:,elnodes(j))
           call cdiff(r1,r0,rr)
           r1=coord_nod2D(:,elnodes(j+1))
           call cdiff(r1,r0,rl)
           call cdiff(rc,r0,r0)
           area(elnodes(j))=area(elnodes(j))+0.25_WP*elem_cos(elem)*( &
                abs(rl(1)*r0(2)-rl(2)*r0(1))+abs(rr(1)*r0(2)-rr(2)*r0(1)))
           w_cv(j-1,elem) = 0.25_WP*elem_cos(elem)*( &
                abs(rl(1)*r0(2)-rl(2)*r0(1))+abs(rr(1)*r0(2)-rr(2)*r0(1)))/elem_area(elem)
!           if(mype==0 .AND. elem==1) write(*,*) 'Elem 1 vectors x, node ',j-1
!           if(mype==0 .AND. elem==1) write(*,*) rr(1), rl(1), rc(1)
!           if(mype==0 .AND. elem==1) write(*,*) rr(2), rl(2), rc(2)
        END DO
     end if
  END DO
  ! ===========
  ! Update to proper dimension
  ! ===========

  elem_area=elem_area*r_earth*r_earth
  area=area*r_earth*r_earth



#ifdef USE_MPI


  call exchange_nod(area)

  call exchange_elem(w_cv)

  call MPI_REDUCE(maxval(w_cv(1,:)), max_w, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
       0, MPI_COMM_FESOM_C, MPIerr)
  call MPI_REDUCE(minval(w_cv(1,:)), min_w, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
       0, MPI_COMM_FESOM_C, MPIerr)

  if (mype==0) write(*,*) 'max_ min_ area_weight',max_w, min_w

#else

  write(*,*) 'max_ min_ area_weight',maxval(w_cv(1,:)),minval(w_cv(1,:))

#endif





!!$  open(10,file='area_dump.out')
!!$
!!$  do n=1,myDim_nod2d
!!$     write(10,*) mype,n,area(n)
!!$  end do
!!$
!!$  close(10)

#ifdef USE_MPI

  call MPI_REDUCE(maxval(elem_area), max_area, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
       0, MPI_COMM_FESOM_C, MPIerr)
  call MPI_REDUCE(minval(elem_area), min_area, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
       0, MPI_COMM_FESOM_C, MPIerr)
  call MPI_REDUCE(maxval(area), max_scarea, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
       0, MPI_COMM_FESOM_C, MPIerr)
  call MPI_REDUCE(minval(area), min_scarea, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
       0, MPI_COMM_FESOM_C, MPIerr)

  if (mype==0) then
     write(*,*) 'Mesh statistics:'
     write(*,*) 'maxArea ',max_area, '   MinArea ', min_area
     write(*,*) 'maxScArea ',max_scarea, '   MinScArea ', min_scarea
     write(*,*) 'Edges:    ', edge2D, ' internal ', edge2D_in
  end if

#else

  !write(*,*) area(1), area(2), area(250), area(50000)
  write(*,*) 'Mesh statistics:'
  write(*,*) 'maxArea ',maxval(elem_area), '   MinArea ', minval(elem_area)
  write(*,*) 'maxScArea ',maxval(area(:)), '   MinScArea ', minval(area(:))
  write(*,*)   'Edges:    ', edge2D, ' internal ', edge2D_in

#endif


END SUBROUTINE mesh_arrays1
!==================================================================
subroutine cdiff(x2,x1,xdiff)

  use o_PARAM
  implicit none

  real(kind=WP) :: x1(2), x2(2), xdiff(2)

  xdiff=x2-x1
  if(xdiff(1)>cyclic_length/2.0_WP) xdiff(1)=xdiff(1)-cyclic_length
  if(xdiff(1)<-cyclic_length/2.0_WP) xdiff(1)=xdiff(1)+cyclic_length

end subroutine cdiff
!===================================================================

SUBROUTINE mesh_arrays2

  ! Collects auxiliary information needed to speed up computations
  ! of gradients, div. This also makes implementation of cyclicity
  ! much more straightforward
  ! Allocated and filled in are:
  ! edge_dxdy(2,edge2D)
  ! edge_cross_dxdy(4,edge2D)
  ! gradient_sca(8,elem2D)
  ! gradient_vec(8,elem2d)

  USE o_MESH
  USE o_PARAM
  use g_comm_auto
  !USE g_ROTATE_grid

  use g_parsup

  IMPLICIT NONE

  integer              :: n,j,q, elnodes(4), ed(2), edg, elem, el(2)
  integer              :: enodes(6)
  real(kind=WP)         :: x1(2), x2(2), xd(2), aa
  real(kind=WP)	     :: a(2), b(2), ax, ay, dfactor, lon, lat
  real(kind=WP)	     :: deltaX31, deltaX21, deltaY31, deltaY21
  real(kind=WP)         :: x(4), y(4), cxx, cxy, cyy, d
  real(kind=WP), allocatable :: center_x(:), center_y(:), temp(:)


  allocate(edge_dxdy(2,myDim_edge2D+eDim_edge2D))        !SH If fields must have this size, check below
  allocate(edge_cross_dxdy(4,myDim_edge2D+eDim_edge2D))  !SH elem_center only valid for myDim_elem2d (NOT +eDim_elem2D)
  !allocate(edge_dxdy(2,myDim_edge2D))
  !allocate(edge_cross_dxdy(4,myDim_edge2D))
  allocate(gradient_sca(8,myDim_elem2D))
  allocate(gradient_vec(8,myDim_elem2D))

  allocate(center_x(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
  allocate(center_y(myDim_elem2D+eDim_elem2D+eXDim_elem2D))


  ! Fill auxiliary arrays center_x and center_y

  do n=1,myDim_elem2D
     call elem_center(n, ax, ay, aa)
     center_x(n)=ax
     center_y(n)=ay
  end do

#ifdef USE_MPI
  call exchange_elem(center_x)
  call exchange_elem(center_y)
#endif

!if (mype==1) then
!  do n=1,myDim_elem2D+eDim_elem2D+eXDim_elem2D
!     print *,n,center_x(n)
!  end do
!end if


!print *,'CENTERX',mype,minval(center_x(:)),maxval(center_x(:))
!print *,'CENTERY',mype,minval(center_y(:)),maxval(center_y(:))

  ! ===========
  ! Distances along the edge
  ! We need them in radian measure!
  ! ===========

  DO n=1, myDim_edge2D+eDim_edge2D
     ed=edge_nodes(:,n)
     a=coord_nod2D(:,ed(2))-coord_nod2D(:, ed(1))
     if(a(1)>cyclic_length/2.0_WP) a(1)=a(1)-cyclic_length
     if(a(1)<-cyclic_length/2.0_WP) a(1)=a(1)+cyclic_length
     edge_dxdy(:,n)=a
  END DO

  ! ===========
  ! Cross-distances for the edge
  ! They are in physical measure!
  ! ===========

  DO edg=1, myDim_edge2D+eDim_edge2D

     ed=edge_nodes(:,edg)
     el=edge_tri(:,edg)

     !call elem_center(el(1), b(1), b(2), aa)
     call edge_center(ed(1), ed(2), a(1), a(2))
     b(1)=center_x(el(1))
     b(2)=center_y(el(1))
     b=b-a

     if(b(1)>cyclic_length/2.0_WP)  b(1)=b(1)-cyclic_length
     if(b(1)<-cyclic_length/2.0_WP) b(1)=b(1)+cyclic_length

     b(1)=b(1)*elem_cos(el(1))
     b=b*r_earth
     edge_cross_dxdy(1:2,edg)=b(1:2)

     if(el(2)>0) then
        b(1)=center_x(el(2))
        b(2)=center_y(el(2))
        b=b-a

        if(b(1)>cyclic_length/2.0_WP) b(1)=b(1)-cyclic_length
        if(b(1)<-cyclic_length/2.0_WP) b(1)=b(1)+cyclic_length

        b(1)=b(1)*elem_cos(el(2))
        b=b*r_earth
        edge_cross_dxdy(3:4,edg)=b(1:2)
     else
        edge_cross_dxdy(3:4,edg)=0.0_WP
     end if

  END DO

!print *,'EDGECR1',mype,minval(edge_cross_dxdy(1,:)),maxval(edge_cross_dxdy(1,:))
!print *,'EDGECR2',mype,minval(edge_cross_dxdy(2,:)),maxval(edge_cross_dxdy(2,:))
!print *,'EDGECR3',mype,minval(edge_cross_dxdy(3,:)),maxval(edge_cross_dxdy(3,:))
!print *,'EDGECR4',mype,minval(edge_cross_dxdy(4,:)),maxval(edge_cross_dxdy(4,:))

  ! ==========================
  ! Derivatives of scalar quantities
  ! ==========================

  DO elem=1, myDim_elem2D

     elnodes=elem2D_nodes(:,elem)
     if(elnodes(1)==elnodes(4)) then
        deltaX31=coord_nod2D(1,elnodes(3))-coord_nod2D(1,elnodes(1))
        if(deltaX31>cyclic_length/2.0_WP) deltaX31=deltaX31-cyclic_length
        if(deltaX31<-cyclic_length/2.0_WP) deltaX31=deltaX31+cyclic_length
        deltaX31=elem_cos(elem)*deltaX31

        deltaX21=coord_nod2D(1,elnodes(2))-coord_nod2D(1,elnodes(1))
        if(deltaX21>cyclic_length/2) deltaX21=deltaX21-cyclic_length
        if(deltaX21<-cyclic_length/2) deltaX21=deltaX21+cyclic_length
        deltaX21=elem_cos(elem)*deltaX21

        deltaY31=coord_nod2D(2,elnodes(3))-coord_nod2D(2,elnodes(1))
        deltaY21=coord_nod2D(2,elnodes(2))-coord_nod2D(2,elnodes(1))

        dfactor=-0.5_WP*r_earth/elem_area(elem)
        gradient_sca(1,elem)=(-deltaY31+deltaY21)*dfactor
        gradient_sca(2,elem)=deltaY31*dfactor
        gradient_sca(3,elem)=-deltaY21*dfactor
        gradient_sca(4,elem)=0.0_WP

        gradient_sca(5,elem)=(deltaX31-deltaX21)*dfactor
        gradient_sca(6,elem)=-deltaX31*dfactor
        gradient_sca(7,elem)=deltaX21*dfactor
        gradient_sca(8,elem)=0.0_WP
     else
        enodes(2:5)=elnodes
        enodes(1)=enodes(5)
        enodes(6)=enodes(2)
        dfactor=0.5_WP*r_earth/elem_area(elem)
        gradient_sca(:,elem)=0.0_WP
        DO j=2,5    ! Nodes are listed clockwise n=(-dy, dx)
           x1=coord_nod2D(:,enodes(j))
           x2=coord_nod2D(:,enodes(j-1))
           call cdiff(x1,x2,xd)           ! xd=(j)-(j-1)
           x2=coord_nod2D(:,enodes(j+1))
           call cdiff(x2,x1,x2)           ! x2=(j+1)-(j)
           xd=(xd+x2)
           xd(1)=xd(1)*elem_cos(elem)
           gradient_sca(j-1,elem)=-xd(2)*dfactor
           gradient_sca(j-1+4,elem)=xd(1)*dfactor
        end do
     end if
  END DO


  ! ==========================
  ! Derivatives of vector quantities
  ! Least squares interpolation is used
  ! ==========================

  DO elem=1,myDim_elem2D
     !call elem_center(elem,a(1),a(2),aa)
     a(1)=center_x(elem)
     a(2)=center_y(elem)
     q=4
     if(elem2D_nodes(1,elem)==elem2D_nodes(4,elem)) q=3
     DO j=1,q
        el(1)=elem_neighbors(j,elem)
        if (el(1)>0) then
           !call elem_center(el(1), b(1), b(2),aa)
           b(1)=center_x(el(1))
           b(2)=center_y(el(1))
           x(j)=b(1)-a(1)
           if(x(j)>cyclic_length/2.0_WP) x(j)=x(j)-cyclic_length
           if(x(j)<-cyclic_length/2.0_WP) x(j)=x(j)+cyclic_length
           y(j)=b(2)-a(2)
        else
           ! Virtual element center is taken
           ed=edge_nodes(:,elem_edges(j,elem))
           call edge_center(ed(1), ed(2), b(1), b(2))
           x(j)=(b(1)-a(1))
           if(x(j)>cyclic_length/2.0_WP)   x(j)=x(j)-cyclic_length
           if(x(j)<-cyclic_length/2.0_WP)  x(j)=x(j)+cyclic_length
           x(j)=2.0_WP*x(j)
           y(j)=2.0_WP*(b(2)-a(2))
        end if
     END DO
     if(q==3) then
        x(4)=0.0_WP
        y(4)=0.0_WP
     end if
     x=x*elem_cos(elem)*r_earth
     y=y*r_earth
     cxx=sum(x**2)
     cxy=sum(x*y)
     cyy=sum(y**2)
     d=cxy*cxy-cxx*cyy
     ! coefficients to compute gradients of velocity components
     ! Metric terms have to be accounted separately if needed

     gradient_vec(1:4,elem)=(cxy*y-cyy*x)/d
     gradient_vec(5:8,elem)=(cxy*x-cxx*y)/d

  END DO

  deallocate(center_y, center_x)

  if (mype==0) print *,'SUBROUTINE mesh_arrays2: COMPLETED'


!print *,'YIPPEE',mype,minval(gradient_vec(1,:)),maxval(gradient_vec(1,:))
!print *,'YIPPEEsca',mype,minval(gradient_sca(1,:)),maxval(gradient_sca(1,:))

END SUBROUTINE mesh_arrays2

! ==================================================================

SUBROUTINE find_up_downwind_triangles

  USE o_MESH
  USE o_ARRAYS
  USE o_PARAM

  USE g_PARSUP
  use g_comm_auto

  IMPLICIT NONE

  integer                    :: n, k, ednodes(2), elem, counter, elnodes(4)
  real(kind=WP)              :: x(2),b(2), c(2), cr, bx, by, xx, xy, ab, ax
  real(kind=WP), allocatable :: coord_elem(:,:,:), temp(:)
  integer, allocatable       :: temp_i(:), e_nodes(:,:)
  integer, allocatable       :: eltp(:) ! Type of element (tri or quad 3 or 4 neighbors)

  edge_up_dn_tri=0
  ! edge_up_dn_grad not used here, allocated only when type_task>1 and initialized with 0,
  ! in case of type_task=1 this will give error!!!!
  ! edge_up_dn_grad=0.0_WP

  ! =====
  ! In order that this procedure works, we need to know nodes and their coordinates
  ! on the extended set of elements (not only my, but myDim+eDim+eXDim)
  ! =====

  allocate(coord_elem(2,4,myDim_elem2D+eDim_elem2D+eXDim_elem2D))
  allocate(temp(myDim_elem2D+eDim_elem2D+eXDim_elem2D))

  DO n=1,4
     DO k=1,2
        temp(1:myDim_elem2D)=coord_nod2D(k,elem2D_nodes(n,:))
#ifdef USE_MPI
        call exchange_elem(temp)
#endif
        coord_elem(k,n,:)=temp(:)
     END DO
  END DO

  deallocate(temp)

#ifndef USE_MPI
  allocate(myList_nod2D(nod2D))
  counter=0
  do n=1, nod2D
     counter=counter+1
     myList_nod2D(counter)=n
  end do
#endif

  allocate(e_nodes(4,myDim_elem2D+eDim_elem2D+eXDim_elem2D))
  allocate(temp_i(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
  allocate(eltp(myDim_elem2D+eDim_elem2D+eXDim_elem2D))

  do elem=1,myDim_elem2D
     elnodes=elem2D_nodes(:,elem)
     if (elnodes(4)==elnodes(1)) then
        eltp(elem)=3
     else
        eltp(elem)=4
     end if
  end do

#ifdef USE_MPI
  call exchange_elem(eltp)
#endif


  DO n=1,4
     temp_i(1:myDim_elem2D) = myList_nod2D(elem2D_nodes(n,:))
#ifdef USE_MPI
     call exchange_elem(temp_i)
#endif
     e_nodes(n,:)=temp_i(:)
  END DO

  deallocate(temp_i)
  !Write(*,*) mype, 'XX', maxval(e_nodes), minval(e_nodes)

  DO n=1,myDim_edge2D
     ednodes=edge_nodes(:,n)
     x=coord_nod2D(:,ednodes(2))-coord_nod2D(:,ednodes(1))

     if(x(1)>cyclic_length/2.0_WP)  x(1)=x(1)-cyclic_length
     if(x(1)<-cyclic_length/2.0_WP) x(1)=x(1)+cyclic_length

     ! Find upwind (in the sense of x) triangle, i. e.
     ! find which triangle contains -x:

     x=-x
     DO k=1,nod_in_elem2D_num(ednodes(1))
        elem=nod_in_elem2D(k,ednodes(1))
        !elnodes=elem2D_nodes(:,elem)
        if (eltp(elem)==3) then      ! triangles
           if(e_nodes(1,elem)==myList_nod2D(ednodes(1))) then
              b=coord_elem(:,2,elem)-coord_elem(:,1,elem)
              c=coord_elem(:,3,elem)-coord_elem(:,1,elem)
           elseif(e_nodes(2,elem)==myList_nod2D(ednodes(1))) then
              b=coord_elem(:,1,elem)-coord_elem(:,2,elem)
              c=coord_elem(:,3,elem)-coord_elem(:,2,elem)
           else
              b=coord_elem(:,1,elem)-coord_elem(:,3,elem)
              c=coord_elem(:,2,elem)-coord_elem(:,3,elem)
           end if
        else                                                ! quads
           if(e_nodes(1,elem)==myList_nod2D(ednodes(1))) then
              b=coord_elem(:,2,elem)-coord_elem(:,1,elem)
              c=coord_elem(:,4,elem)-coord_elem(:,1,elem)
           elseif(e_nodes(2,elem)==myList_nod2D(ednodes(1))) then
              b=coord_elem(:,3,elem)-coord_elem(:,2,elem)
              c=coord_elem(:,1,elem)-coord_elem(:,2,elem)
           elseif(e_nodes(3,elem)==myList_nod2D(ednodes(1))) then
              b=coord_elem(:,4,elem)-coord_elem(:,3,elem)
              c=coord_elem(:,2,elem)-coord_elem(:,3,elem)
           else
              b=coord_elem(:,1,elem)-coord_elem(:,4,elem)
              c=coord_elem(:,3,elem)-coord_elem(:,4,elem)
           end if
        end if
        if(b(1)>cyclic_length/2.) b(1)=b(1)-cyclic_length
        if(b(1)<-cyclic_length/2.) b(1)=b(1)+cyclic_length
        if(c(1)>cyclic_length/2.) c(1)=c(1)-cyclic_length
        if(c(1)<-cyclic_length/2.) c(1)=c(1)+cyclic_length
        ! the vector x has to be between b and c
        ! Decompose b and x into parts along c and along (-cy,cx), i.e.
        ! 90 degree counterclockwise
        cr=sum(c*c)
        if (cr==0.0_WP) then
            write(*,*) "mype=",mype," fv_mesh_array, find_up_downwind_triangles: cr=0 at elem",elem," coord.:",coord_elem(:,:,elem)
            write(*,*) "mype=",mype,"myDim_elem2D = ",myDim_elem2D,eDim_elem2D,eXDim_elem2D
        end if
        bx=sum(b*c)/cr
        by=(-b(1)*c(2)+b(2)*c(1))/cr
        xx=sum(x*c)/cr
        xy=(-x(1)*c(2)+x(2)*c(1))/cr
        ab=atan2(by,bx)
        ax=atan2(xy,xx)
        ! Since b and c are the sides of triangle, |ab|<pi, and atan2 should
        ! be what is needed
        if((ab>0).and.(ax>0).and.(ax<ab)) then
           edge_up_dn_tri(1,n)=elem
           cycle
        endif
        if((ab<0).and.(ax<0).and.(ax>ab)) then
           edge_up_dn_tri(1,n)=elem
           cycle
        endif
        if((ab==ax).or.(ax==0.0)) then
           edge_up_dn_tri(1,n)=elem
           cycle
        endif
     END DO

     ! Find downwind element
     x=-x
     DO k=1,nod_in_elem2D_num(ednodes(2))
        elem=nod_in_elem2D(k,ednodes(2))
        !elnodes=elem2D_nodes(:,elem)
        if(eltp(elem)==3) then      ! triangles
           if(e_nodes(1,elem)==myList_nod2D(ednodes(2))) then
              b=coord_elem(:,2,elem)-coord_elem(:,1,elem)
              c=coord_elem(:,3,elem)-coord_elem(:,1,elem)
           elseif(e_nodes(2, elem)==myList_nod2D(ednodes(2))) then
              b=coord_elem(:,1,elem)-coord_elem(:,2,elem)
              c=coord_elem(:,3,elem)-coord_elem(:,2,elem)
           else
              b=coord_elem(:,1,elem)-coord_elem(:,3,elem)
              c=coord_elem(:,2,elem)-coord_elem(:,3,elem)
           end if
        else                                 ! quads
           if(e_nodes(1,elem)==myList_nod2D(ednodes(2))) then
              b=coord_elem(:,2,elem)-coord_elem(:,1,elem)
              c=coord_elem(:,4,elem)-coord_elem(:,1,elem)
           elseif(e_nodes(2,elem)==myList_nod2D(ednodes(2))) then
              b=coord_elem(:,3,elem)-coord_elem(:,2,elem)
              c=coord_elem(:,1,elem)-coord_elem(:,2,elem)
           elseif(e_nodes(3,elem)==myList_nod2D(ednodes(2))) then
              b=coord_elem(:,4,elem)-coord_elem(:,3,elem)
              c=coord_elem(:,2,elem)-coord_elem(:,3,elem)
           else
              b=coord_elem(:,1,elem)-coord_elem(:,4,elem)
              c=coord_elem(:,3,elem)-coord_elem(:,4,elem)
           end if
        end if
        if(b(1)>cyclic_length/2.0_WP) b(1)=b(1)-cyclic_length
        if(b(1)<-cyclic_length/2.0_WP) b(1)=b(1)+cyclic_length
        if(c(1)>cyclic_length/2.0_WP) c(1)=c(1)-cyclic_length
        if(c(1)<-cyclic_length/2.0_WP) c(1)=c(1)+cyclic_length
        ! the vector x has to be between b and c
        ! Decompose b and x into parts along c and along (-cy,cx), i.e.
        ! 90 degree counterclockwise
        cr=sum(c*c)
        bx=sum(b*c)/cr
        by=(-b(1)*c(2)+b(2)*c(1))/cr
        xx=sum(x*c)/cr
        xy=(-x(1)*c(2)+x(2)*c(1))/cr
        ab=atan2(by,bx)
        ax=atan2(xy,xx)
        ! Since b and c are the sides of triangle, |ab|<pi, and atan2 should
        ! be what is needed
        if((ab>0).and.(ax>0).and.(ax<ab)) then
           edge_up_dn_tri(2,n)=elem
           cycle
        endif
        if((ab<0).and.(ax<0).and.(ax>ab)) then
           edge_up_dn_tri(2,n)=elem
           cycle
        endif
        if((ab==ax).or.(ax==0.0)) then
           edge_up_dn_tri(2,n)=elem
           cycle
        endif
     END DO
  END DO

  ! There is problem with edges close to boundary --- they may be lacking
  ! up or downwind elements. We have to return to the standard Miura at nodes that
  ! belong to such edges. Same issue is occurring at the depth.

  deallocate(e_nodes, coord_elem, eltp)

#ifndef USE_MPI
  deallocate(myList_nod2D)
#endif

  if (mype==0) print *,'SUBROUTINE find_up_downwind_triangles COMPLETED'

end SUBROUTINE find_up_downwind_triangles
!=======================================================================

