! Driving routine. The distributed mesh information and mesh proper
! are read from files.
! Auxiliary arrays with mesh information are assembled.
! At the beginning of each routine I list arrays it initializes.
! Array sizes vary (sometimes we need only myDim, yet sometimes more)!
! S. Danilov, 2012

SUBROUTINE mesh_setup

  USE g_parsup
  USE g_ROTATE_grid

  IMPLICIT NONE

  call set_mesh_transform_matrix  !(rotated grid)
  call read_mesh_par
  call set_par_support
!!$  call find_levels
!!$  call test_tri
!!$  call load_edges
!!$  call find_neighbors
!!$  call mesh_areas
!!$  call mesh_auxiliary_arrays

END SUBROUTINE mesh_setup
!======================================================================
! Reads distributed mesh
! The mesh will be read only by 0 proc and broadcasted to the others.
SUBROUTINE read_mesh_par

  USE o_PARAM
  !!SH USE g_CONFIG
  USE o_MESH
  USE o_ARRAYS
  USE g_PARSUP
  USE g_rotate_grid
  use o_utilit

  use g_comm_auto

  IMPLICIT NONE

  integer        :: n, nn, k, m, fileID
  integer        :: error_status !0/1=no error/error
  integer        :: vert_nodes(1000)
  integer        :: nchunk, chunk_size, ipos, iofs, mesh_check
  real(kind=WP)  :: x, y, rx, ry
  real(kind=WP)  :: t0, t1
  character*10   :: mype_string,npes_string
  character*500  :: file_name
  character*500  :: dist_mesh_dir
  integer       :: ierror              ! return error code
  integer, allocatable, dimension(:)        :: mapping
  integer, allocatable, dimension(:,:)      :: ibuff
  real(kind=WP), allocatable, dimension(:,:) :: rbuff
  integer, allocatable, dimension(:,:)      :: auxbuff ! will be used for reading aux3d.out

  integer :: ioerror
  integer :: n_quad, ntri, nmb_nds_el
  INTEGER, allocatable  :: elem_data(:)
  integer :: ind

  !mesh related files will be read in chunks of chunk_size
  chunk_size=100000
  !==============================
  ! Allocate mapping array (chunk_size)
  ! It will be used for several purposes
  !==============================
  allocate(mapping(chunk_size))
  allocate(ibuff(chunk_size,4), rbuff(chunk_size,3))

  mapping=0
  !==============================
  t0=MPI_Wtime()
  write(mype_string,'(i5.5)') mype
  write(npes_string,"(I10)") npes
  dist_mesh_dir=trim(meshpath)//'dist_'//trim(ADJUSTL(npes_string))//'/'

  !=======================
  ! rank partitioning vector
  ! will be read by 0 proc
  !=======================
  if (mype==0) then
     file_name=trim(dist_mesh_dir)//'rpart.out'
     fileID=10
     open(fileID, file=trim(file_name), iostat=ioerror)
     if (ioerror/=0) then
        error_status=2
     else
        allocate(part(npes+1))
        read(fileID,*) n
        error_status=0
        if (n/=npes) error_status=1 !set the error status for consistency in rpart
        part(1)=1
        read(fileID,*) part(2:npes+1)
        DO n=2, npes+1
           part(n)=part(n-1)+part(n)
        END DO
        close(fileID)
     end if
  end if

  ! communicate and check the error status
  call MPI_BCast(error_status, 1, MPI_INTEGER, 0, MPI_COMM_FESOM_C, ierror)
  if (error_status/=0) then
     if (error_status==1) then
        if (mype==0) write(*,*) 'Number of partitions: n = ',n
        if (mype==0) write(*,*) 'error: NPES does not coincide with that of the mesh (n)'
     elseif (error_status==2) then
        if (mype==0) write(*,*) 'File rpart.out does not exist for specified number of cores NPES=', npes
     end if
     call par_ex(1)
     STOP
  end if

  ! broadcasting partitioning vector to the other procs
  if (mype/=0) then
     allocate(part(npes+1))
  end if
  call MPI_BCast(part, npes+1, MPI_INTEGER, 0, MPI_COMM_FESOM_C, ierror)
  if (mype==0) write(*,*) mype,'rpart is read'

  !===========================
  ! Lists of nodes and elements in global indexing.
  ! every proc reads its file
  !===========================

  file_name=trim(dist_mesh_dir)//'my_list'//trim(mype_string)//'.out'
  fileID=10+mype

  open(fileID, file=trim(file_name))
  read(fileID,*) n

  read(fileID,*) myDim_nod2D
  read(fileID,*) eDim_nod2D
  allocate(myList_nod2D(myDim_nod2D+eDim_nod2D))
  read(fileID,*) myList_nod2D

  read(fileID,*) myDim_elem2D
  read(fileID,*) eDim_elem2D
  read(fileID,*) eXDim_elem2D
  allocate(myList_elem2D(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
  read(fileID,*) myList_elem2D

  read(fileID,*) myDim_edge2D
  read(fileID,*) eDim_edge2D
  allocate(myList_edge2D(myDim_edge2D+eDim_edge2D))
  read(fileID,*) myList_edge2D ! m

  close(fileID)
  if (mype==0) write(*,*) 'myLists are read'

  !==============================
  ! read 2d node data
  !==============================
  ! read the nod2D from nod2d.out and check whether it is equal to part(npes+1)-1
  nod2D=part(npes+1)-1
  allocate(coord_nod2D(2,myDim_nod2D+eDim_nod2D))
  allocate(index_nod2D(myDim_nod2D+eDim_nod2D))
  index_nod2D(:)=0
  if (mype==0) then
     file_name=trim(meshpath)//'nod2d.out'
     open(fileID, file=file_name)
     read(fileID,*) n      ! nod2D, we know it already

     error_status=0
     nobn=0

     if (n/=nod2D) error_status=1 !set the error status for consistency between rpart and nod2D
     write(*,*) 'reading '// trim(file_name)
  end if
  ! check the error status
  call MPI_BCast(error_status, 1, MPI_INTEGER, 0, MPI_COMM_FESOM_C, ierror)
  if (error_status/=0) then
     write(*,*) n
     write(*,*) 'error: nod2D/=part(npes+1)-1'
     call par_ex(1)
     STOP
  end if

  ! 0 proc reads the data in chunks and distributes it between other procs
  mesh_check=0
  do nchunk=0, (nod2D-1)/chunk_size
     !create the mapping for the current chunk
     mapping(1:chunk_size)=0
     do n=1, myDim_nod2D+eDim_nod2D
        ipos=(myList_nod2D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_nod2D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do

     !read the chunk into the buffers
     k=min(chunk_size, nod2D-nchunk*chunk_size)
     if (mype==0) then
        do n=1, k
           read(fileID,*) ibuff(n,1), rbuff(n,1), rbuff(n,2), ibuff(n,2)
            if (ibuff(n,2)==2) nobn=nobn+1
        end do
     end if
     call MPI_BCast(rbuff(1:k,1), k, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
     call MPI_BCast(rbuff(1:k,2), k, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)
     call MPI_BCast(ibuff(1:k,2), k, MPI_INTEGER, 0, MPI_COMM_FESOM_C, ierror)
     ! fill the local arrays
     do n=1, k

        if (cartesian) then
           x = rbuff(n,1)/r_earth
           y = rbuff(n,2)/r_earth
        else
           x = rbuff(n,1)*rad
           y = rbuff(n,2)*rad
        end if

        if (force_rotation) then
           rx=x
           ry=y
           call g2r(rx, ry, x, y)
        end if

        if (mapping(n)>0) then
           mesh_check=mesh_check+1
           coord_nod2D(1,mapping(n))=x
           coord_nod2D(2,mapping(n))=y
           index_nod2D(mapping(n))=ibuff(n,2)
        end if

     end do
  end do
  if (mype==0) close(fileID)


  call MPI_BCast(nobn, 1, MPI_INTEGER, 0, MPI_COMM_FESOM_C, ierror)

  if (mype==0) print *,'Total number of open boundary nodes : ',nobn

  if (mesh_check/=myDim_nod2D+eDim_nod2D) then
     write(*,*) 'ERROR while reading nod2d.out on mype=', mype
     write(*,*) mesh_check, ' values have been read in according to partitioning'
     write(*,*) 'it does not equal to myDim_nod2D+eDim_nod2D = ', myDim_nod2D+eDim_nod2D
  end if

  !local number of boundary nodes
  mynobn=0
  do n=1,myDim_nod2D
     if (index_nod2D(n)==2) mynobn=mynobn+1
  end do
!  print *,'PE: ',mype,' Local number of open boundary nodes: ',mynobn

  if (mynobn>0) then
     allocate(my_in_obn(mynobn), my_in_obn_idx(mynobn))
     my_in_obn(:)=0; my_in_obn_idx(:)=0
     ind=0
     do n=1,myDim_nod2D
        if (index_nod2D(n)==2) then
           ind=ind+1
           my_in_obn(ind)=n
        endif
     end do
  end if

  if (mynobn>0 .AND. ind/=mynobn) print *,'ATTENTION: boundary node discrepancy'
!  if (mynobn>10) print *,'global nodes MY_IN_OBN(1:10):',mype,myList_nod2D(my_in_obn(1:10))

  !==============================
  ! read 2d elem data
  !==============================
  ! read the elem2D from elem2d.out
  if (mype==0)  then
     file_name=trim(meshpath)//'elem2d.out'
     open(fileID, file=file_name)
     read(fileID,*) elem2d
     write(*,*) 'reading '// trim(file_name)

     ! Check if there are 4 columns in the file
     allocate(elem_data(4*elem2D))
     read(fileID,*,iostat=ioerror) elem_data(1:4*elem2D)
     write(*,*) ioerror
     if (ioerror == 0) then
        nmb_nds_el=4
     else
        nmb_nds_el=3
     end if
     write(*,*) 'Number of nodes per element : ',nmb_nds_el
     deallocate(elem_data)

     close(fileID)

     !Re-open the file for further processing
     open(fileID, file=file_name)
     read(fileID,*) elem2d

  end if

  call MPI_BCast(nmb_nds_el, 1, MPI_INTEGER, 0, MPI_COMM_FESOM_C, ierror)
  call MPI_BCast(elem2d, 1, MPI_INTEGER, 0, MPI_COMM_FESOM_C, ierror)
  !SH Adjusted the size for fesom_c ...
  !allocate(elem2D_nodes(4, myDim_elem2D+eDim_elem2D+eXDim_elem2D))
  !...and back
  allocate(elem2D_nodes(4, myDim_elem2D))

  ! 0 proc reads the data in chunks and distributes it between other procs
  do nchunk=0, (elem2D-1)/chunk_size
     mapping(1:chunk_size)=0
     !do n=1, myDim_elem2D
     !SH Adjusted the size for fesom_c
     do n=1, myDim_elem2D !and back +eDim_elem2D+eXDim_elem2D
        ipos=(myList_elem2D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_elem2D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do

     k=min(chunk_size, elem2D-nchunk*chunk_size)
     if (mype==0) then
        do n=1, k
           if (nmb_nds_el==4) then
              read(fileID,*) ibuff(n,1), ibuff(n,2), ibuff(n,3), ibuff(n,4)
           else
              read(fileID,*) ibuff(n,1), ibuff(n,2), ibuff(n,3)
           end if
        end do
     end if

     call MPI_BCast(ibuff(1:k,1), k, MPI_INTEGER, 0, MPI_COMM_FESOM_C, ierror)
     call MPI_BCast(ibuff(1:k,2), k, MPI_INTEGER, 0, MPI_COMM_FESOM_C, ierror)
     call MPI_BCast(ibuff(1:k,3), k, MPI_INTEGER, 0, MPI_COMM_FESOM_C, ierror)
     if (nmb_nds_el==4) then
        call MPI_BCast(ibuff(1:k,4), k, MPI_INTEGER, 0, MPI_COMM_FESOM_C, ierror)
     end if

     do n=1, k
        if (mapping(n)>0) then
           elem2D_nodes(1,mapping(n))=ibuff(n,1)
           elem2D_nodes(2,mapping(n))=ibuff(n,2)
           elem2D_nodes(3,mapping(n))=ibuff(n,3)
           if (nmb_nds_el==4) then
              elem2D_nodes(4,mapping(n))=ibuff(n,4)
           end if
        end if
     end do
  end do

  if (mype==0) close(fileID)

!print *,'READING',mype,myDim_elem2D,eDim_elem2D,eXDim_elem2D

  ! nodes in elem2d are in global numbering. convert to local:
  do nchunk=0, (nod2D-1)/chunk_size
     mapping(1:chunk_size)=0
     do n=1, myDim_nod2D+eDim_nod2D
        ipos=(myList_nod2D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_nod2D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do
     do n=1, myDim_elem2D !and back +eDim_elem2D+eXDim_elem2D
        do m=1,4  !SH ????
           nn=elem2D_nodes(m, n)
           ipos=(nn-1)/chunk_size
           if (ipos==nchunk) then
              iofs=nn-nchunk*chunk_size
              ! minus sign is required to avoid modified entry being modified in another chunk
              ! will be changed to plus at the end
              elem2D_nodes(m,n)=-mapping(iofs)
           end if
        end do
     end do
  end do
  elem2D_nodes=-elem2D_nodes
  if (mype==0) write(*,*) 'elements are read'

  ! Set vertical coordinates
  ! FESOM:   number of levels:  nl
  !          array with levels: zbar
  !          read depth levels from file aux3d.out
  ! FESOM_C: Currently call SET_SIGMA:
  !          number of levels: nsigma
  !          array with levels: sigma


  call SET_SIGMA

  ! Replace lines below by proper rules depending on experiment
  allocate(depth(myDim_nod2D+eDim_nod2D))
  depth(:) = 20._WP


  if (len(trim(title))==2 .AND. title=='LE') then
     depth=20.0_WP
  else

     if (mype==0) open(22,file=trim(meshpath)//'depth.out', status='old')

     ! 0 proc reads the data in chunks and distributes it between other procs
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

        k=min(chunk_size, nod2D-nchunk*chunk_size)
        if (mype==0) then
           do n=1, k
              read(22,*) rbuff(n,1)
           end do
!SHTEST TOPOGRAPHY
!do n=1,k
!  if (rbuff(n,1)<10.0) rbuff(n,1)=10.0
!end do
        end if
        call MPI_BCast(rbuff(1:k,1), k, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM_C, ierror)

        do n=1, k
           x=rbuff(n,1)
           !if (x>zbar(5)) x=zbar(5) !threshold for depth
           if (mapping(n)>0) then
              mesh_check=mesh_check+1
              depth(mapping(n))=x
           end if
        end do
     end do

     if (mype==0) close(22)

     if (mesh_check/=myDim_nod2D+eDim_nod2D) then
        write(*,*) 'ERROR while reading depth.out on mype=', mype
        write(*,*) mesh_check, ' values have been read in according to partitioning'
        write(*,*) 'it does not equal to myDim_nod2D+eDim_nod2D = ', myDim_nod2D+eDim_nod2D
     end if

  end if

!do n=1,mydim_nod2D
! if (index_nod2D(n)==2) print *,'DEPP',mype,myList_nod2D(n),depth(n)
!end do
!do n=1,mydim_nod2D
! if (myList_nod2D(n)==3345)  print *,'KONTROLL',mype,n,myList_nod2D(n),depth(n)
!end do

!  if (mype==0) print *,'DEPTH values are set.'


  ! ==============================
  ! Communication information
  ! every proc reads its file
  ! ==============================
  file_name=trim(dist_mesh_dir)//'com_info'//trim(mype_string)//'.out'
  fileID=10+mype
  open(fileID, file=file_name)
  read(fileID,*)  n
  read(fileID,*) com_nod2D%rPEnum
  if (com_nod2D%rPEnum > MAX_NEIGHBOR_PARTITIONS) then
     print *,'Increase MAX_NEIGHBOR_PARTITIONS in gen_modules_partitioning.F90 and recompile'
     stop
  endif
!!$ ALLOCATE(com_nod2D%rPE(com_nod2D%rPEnum))
  read(fileID,*) com_nod2D%rPE(1:com_nod2D%rPEnum)
!!$  ALLOCATE(com_nod2D%rptr(com_nod2D%rPEnum+1))
  read(fileID,*) com_nod2D%rptr(1:com_nod2D%rPEnum+1)
  ALLOCATE(com_nod2D%rlist(eDim_nod2D))
  read(fileID,*) com_nod2D%rlist

  read(fileID,*) com_nod2D%sPEnum
  if (com_nod2D%sPEnum > MAX_NEIGHBOR_PARTITIONS) then
     print *,'Increase MAX_NEIGHBOR_PARTITIONS in gen_modules_partitioning.F90 and recompile'
     stop
  endif
!!$  ALLOCATE(com_nod2D%sPE(com_nod2D%sPEnum))
  read(fileID,*) com_nod2D%sPE(1:com_nod2D%sPEnum)
!!$  ALLOCATE(com_nod2D%sptr(com_nod2D%sPEnum+1))
  read(fileID,*) com_nod2D%sptr(1:com_nod2D%sPEnum+1)
  n=com_nod2D%sptr(com_nod2D%sPEnum+1)-1
  ALLOCATE(com_nod2D%slist(n))
  read(fileID,*) com_nod2D%slist

  read(fileID,*) com_elem2D%rPEnum
  if (com_elem2D%rPEnum > MAX_NEIGHBOR_PARTITIONS) then
     print *,'Increase MAX_NEIGHBOR_PARTITIONS in gen_modules_partitioning.F90 and recompile'
     stop
  endif
!!$  ALLOCATE(com_elem2D%rPE(com_elem2D%rPEnum))
  read(fileID,*) com_elem2D%rPE(1:com_elem2D%rPEnum)
!!$  ALLOCATE(com_elem2D%rptr(com_elem2D%rPEnum+1))
  read(fileID,*) com_elem2D%rptr(1:com_elem2D%rPEnum+1)
  ALLOCATE(com_elem2D%rlist(eDim_elem2D))
  read(fileID,*) com_elem2D%rlist

  read(fileID,*) com_elem2D%sPEnum
  if (com_elem2D%sPEnum > MAX_NEIGHBOR_PARTITIONS) then
     print *,'Increase MAX_NEIGHBOR_PARTITIONS in gen_modules_partitioning.F90 and recompile'
     stop
  endif
!!$  ALLOCATE(com_elem2D%sPE(com_elem2D%sPEnum))
  read(fileID,*) com_elem2D%sPE(1:com_elem2D%sPEnum)
!!$  ALLOCATE(com_elem2D%sptr(com_elem2D%sPEnum+1))
  read(fileID,*) com_elem2D%sptr(1:com_elem2D%sPEnum+1)
  n=com_elem2D%sptr(com_elem2D%sPEnum+1)-1
  ALLOCATE(com_elem2D%slist(n))
  read(fileID,*) com_elem2D%slist

  read(fileID,*) com_elem2D_full%rPEnum
  if (com_elem2D_full%rPEnum > MAX_NEIGHBOR_PARTITIONS) then
     print *,'Increase MAX_NEIGHBOR_PARTITIONS in gen_modules_partitioning.F90 and recompile'
     stop
  endif
!!$  ALLOCATE(com_elem2D_full%rPE(com_elem2D_full%rPEnum))
  read(fileID,*) com_elem2D_full%rPE(1:com_elem2D_full%rPEnum)
!!$  ALLOCATE(com_elem2D_full%rptr(com_elem2D_full%rPEnum+1))
  read(fileID,*) com_elem2D_full%rptr(1:com_elem2D_full%rPEnum+1)
  ALLOCATE(com_elem2D_full%rlist(eDim_elem2D+eXDim_elem2D))
  read(fileID,*) com_elem2D_full%rlist

  read(fileID,*) com_elem2D_full%sPEnum
  if (com_elem2D_full%sPEnum > MAX_NEIGHBOR_PARTITIONS) then
     print *,'Increase MAX_NEIGHBOR_PARTITIONS in gen_modules_partitioning.F90 and recompile'
     stop
  endif
!!$  ALLOCATE(com_elem2D_full%sPE(com_elem2D_full%sPEnum))
  read(fileID,*) com_elem2D_full%sPE(1:com_elem2D_full%sPEnum)
!!$  ALLOCATE(com_elem2D_full%sptr(com_elem2D_full%sPEnum+1))
  read(fileID,*) com_elem2D_full%sptr(1:com_elem2D_full%sPEnum+1)
  n=com_elem2D_full%sptr(com_elem2D_full%sPEnum+1)-1
  ALLOCATE(com_elem2D_full%slist(n))
  read(fileID,*) com_elem2D_full%slist

!!$ read(fileID,*) com_edge2D%rPEnum
!!$ ALLOCATE(com_edge2D%rPE(com_edge2D%rPEnum))
!!$ read(fileID,*) com_edge2D%rPE
!!$ ALLOCATE(com_edge2D%rptr(com_edge2D%rPEnum+1))
!!$ read(fileID,*) com_edge2D%rptr
!!$ ALLOCATE(com_edge2D%rlist(eDim_edge2D))
!!$ read(fileID,*) com_edge2D%rlist
!!$
!!$ read(fileID,*) com_edge2D%sPEnum
!!$ ALLOCATE(com_edge2D%sPE(com_edge2D%sPEnum))
!!$ read(fileID,*) com_edge2D%sPE
!!$ ALLOCATE(com_edge2D%sptr(com_edge2D%sPEnum+1))
!!$ read(fileID,*) com_edge2D%sptr
!!$ n=com_edge2D%sptr(com_edge2D%sPEnum+1)-1
!!$ ALLOCATE(com_edge2D%slist(n))
!!$ read(fileID,*) com_edge2D%slist
  close(fileID)
  if (mype==0) write(*,*) 'communication arrays are read'

  CALL MPI_BARRIER(MPI_COMM_FESOM_C, MPIerr)

  t1=MPI_Wtime()
  if (mype==0) then
     write(*,*) '========================='
     write(*,*) '2D mesh was read in ', t1-t0, ' seconds'
     write(*,*) '2D mesh info : ', 'nod2D=', nod2D,' elem2D=', elem2D
     write(*,*) '========================='
  endif

  deallocate(rbuff, ibuff)
  deallocate(mapping)

END subroutine  read_mesh_par
!============================================================
