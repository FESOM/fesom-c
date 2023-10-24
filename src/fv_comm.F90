!NR Adapted from FESOM for mixed meshes (quads + triangles).
!NR !!!!!! Not tested !!!!!

! ===============================================================
!=======================================================================
module local_parameters_fv_comm
  integer, parameter   :: MAX_LAENDERECK=8
end module local_parameters_fv_comm

subroutine communication_nodn
  use o_MESH
  use g_PARSUP
  use local_parameters_fv_comm
  implicit none

  integer              :: n,np, prank, el, r_count, s_count, q, i, j, nod, k
  integer              :: num_send(0:npes-1), num_recv(0:npes-1)
  integer              :: recv_from_pe(nod2D)
  integer              :: send_to_pes(MAX_LAENDERECK,nod2D)

  logical              :: max_laendereck_too_small=.false.
  
  ! Assume we have 2D partitioning vector in part. Find communication
  ! rules
  ! allocate(aux(npes,nod2D))
  ! Reduce allocation: find all neighboring PE


  num_send(0:npes-1) = 0
  num_recv(0:npes-1) = 0
  recv_from_pe(1:nod2D) = -1
  send_to_pes(1:MAX_LAENDERECK,1:nod2D) = -1
  
  ! For the local nodes, run through the patch and collect all nodes in the patch
  ! (It would be simpler to re-use the adjacency matrix, but it is not a global
  !  variable... and I am lazy and want to keep the data structure simple)
  do n=1,nod2D
     if (part(n)==mype) then
        do i = 1, nod_in_elem2D_num(n)
           el = nod_in_elem2D(i,n)

           q = 4
           if (elem2D_nodes(1,el) == elem2D_nodes(4,el)) q=3  ! triangle

           do j = 1,q
              nod = elem2D_nodes(j,el)

              if (part(nod) /= mype) then
                 
                 if (recv_from_pe(nod)==-1) then  ! nod already collected to be received?
                    ! We have to receive this node from part(nod)
                    num_recv(part(nod)) = num_recv(part(nod)) + 1
                    recv_from_pe(nod) = part(nod)  ! no new information, just handy
                 endif

                 ! And we have to send n to part(nod). Do we know this already?
                 do k=1,MAX_LAENDERECK
                    if (send_to_pes(k,n) == part(nod)) then
                       exit  ! already collected
                    elseif (send_to_pes(k,n) == -1) then
                       send_to_pes(k,n) = part(nod)
                       num_send(part(nod)) = num_send(part(nod)) + 1
                       exit
                    elseif (k== MAX_LAENDERECK) then
                       max_laendereck_too_small = .true.  ! Problem
                    endif
                 enddo
              endif
           enddo
           
        enddo
     endif
  enddo

  if (max_laendereck_too_small) then
     print *,'Set MAX_LAENDERECK in fv_comm.f90 to a larger value and recompile'
     stop
  endif
    
! Now, build the send and recv communication data structure
  com_nod2D%rPEnum = count(num_recv(0:npes-1) > 0)
  com_nod2D%sPEnum = count(num_send(0:npes-1) > 0)


  allocate(com_nod2D%rPE(com_nod2D%rPEnum))
  allocate(com_nod2D%sPE(com_nod2D%sPEnum))

  r_count = 0
  s_count = 0
  do np = 0, npes-1
     if(num_recv(np) /= 0) then
        r_count = r_count+1
        com_nod2D%rPE(r_count) = np
     end if
     if(num_send(np) /= 0) then
        s_count = s_count+1
        com_nod2D%sPE(s_count) = np
     end if
  enddo


  r_count = 0
  s_count = 0
  allocate(com_nod2D%rptr(com_nod2D%rPEnum+1)) 
  allocate(com_nod2D%sptr(com_nod2D%sPEnum+1))
  
  com_nod2D%rptr(1) = 1
  com_nod2D%sptr(1) = 1
  
  do r_count = 1, com_nod2D%rPEnum
     np = com_nod2D%rPE(r_count)
     com_nod2D%rptr(r_count+1) =  com_nod2D%rptr(r_count)+ num_recv(np)
  enddo
  do s_count = 1, com_nod2D%sPEnum
     np = com_nod2D%sPE(s_count)
     com_nod2D%sptr(s_count+1) =  com_nod2D%sptr(s_count)+ num_send(np)
  enddo

  ! Lists itself

  r_count = 0
  allocate(com_nod2D%rlist(com_nod2D%rptr(com_nod2D%rPEnum+1)-1)) 
  do np = 1,com_nod2D%rPEnum
     prank = com_nod2D%rPE(np)
     do n = 1, nod2D
        if (recv_from_pe(n) == prank) then
           r_count = r_count+1
           com_nod2D%rlist(r_count) = n
        end if
     end do
  end do
  
  s_count = 0
  allocate(com_nod2D%slist(com_nod2D%sptr(com_nod2D%sPEnum+1)-1)) 
  do np = 1,com_nod2D%sPEnum
     prank = com_nod2D%sPE(np)
     do n = 1, nod2D
        if( any(send_to_pes(:,n) == prank)) then 
           s_count = s_count+1
           com_nod2D%slist(s_count) = n
        end if
     end do
  end do

  ! Summary of this piece: mype receives
  ! information on external 2D nodes from
  ! comm_nod2D%rPEnum external PEs
  ! Their ranks (numbers) are in array
  ! comm_nod2D%rPE(:)
  ! Pointers to external node numbers are in
  ! comm_nod2D%rptr(:)
  ! The node numbers are in 
  ! comm_nod2D%list(:)
  ! Putting everything into structure takes many operations, but
  ! without the structure we will need to many names and arrays
  ! Do not forget that we need also send part.

  ! mype sends its data to
  ! comm_nod2D%sPEnum external PEs
  ! Their ranks (numbers) are in array
  ! comm_nod2D%sPE(:)
  ! Pointers to external node numbers are in
  ! comm_nod2D%sptr(:)
  ! The node numbers are in 
  ! comm_nod2D%list(:)
  

end subroutine communication_nodn

!==========================================================================
subroutine communication_elemn
 
  use o_MESH
  use g_PARSUP
  use local_parameters_fv_comm
  implicit none

  integer              :: recv_from_pe(elem2D)
  integer              :: send_to_pes(MAX_LAENDERECK,elem2D)

  logical              :: max_laendereck_too_small=.false.

  integer              :: n, k, ep, np, prank, el, nod
  integer              :: p, q, j, elem, i, r_count, s_count
  
  integer              :: num_send(0:npes-1), num_recv(0:npes-1)

  
  ! Assume we have 2D partitioning vector in part. Find communication
  ! rules. An elem is external to element n if neither of its nodes 
  ! belongs to PE, but it is among the neighbors. Element n belongs to PE if 
  ! any of its nodes does. 
  
  ! This routine takes into account 
  ! com_elem2D_full: all  neighbors  
  ! com_elem2D:      only those sharing an edge 


  !===========================================
  !  com_elem2D
  !===========================================
  
  num_send(0:npes-1) = 0
  num_recv(0:npes-1) = 0
  recv_from_pe(1:elem2D) = -1
  send_to_pes(1:MAX_LAENDERECK,1:elem2D) = -1
  
  ! For the local elements, collect all adjacent (sharing an edge) elements
  ! that belong to other PEs
  
  do el=1,elem2D
     ! Only elements that belong to mype _and_ some other PE are part of the inner halo
     if (any(part(elem2D_nodes(1:4,el)) == mype) .and. &
         any(part(elem2D_nodes(1:4,el)) /= mype) )  then

        q = 4
        if (elem2D_nodes(1,el) == elem2D_nodes(4,el)) q=3  ! triangle

        do j = 1,q

           elem = elem_neighbors(j,el)

           if (elem < 1) cycle  ! boundary, "ghost element"
           
           if (all(part(elem2D_nodes(:,elem)) /= mype)) then  ! This element does not belong to mype
       
              if (recv_from_pe(elem)==-1) then  ! elem to be received already collected?
                 ! We have to receive elem from PE ep:
                 ep = part(elem2D_nodes(1,elem))  ! PE of first node is "main" owner
                 num_recv(ep) = num_recv(ep) + 1
                 recv_from_pe(elem) = ep  
              endif

              ! And maybe, we have to send el to the owners of elem
              ! This gets more complicated:
              ! 1. Is mype the main owner of el?
              if (part(elem2D_nodes(1,el)) == mype) then
                    
                 ! 2. who owns elem? We must check all nodes!
                 ! For triangles, we check elem2D_nodes(1,el)=elem2D_nodes(4,el) twice, this does not matter
                 do i=1,4
                    ep = part(elem2D_nodes(i,elem))

                    ! 3. Is ep also an owner of el and no send is needed? This excludes also mype==ep.
                    if (any(part(elem2D_nodes(1:4,el)) == ep)) cycle
                    
                    ! 4. Ok, for the owner ep, check if sending el is already collected
                    do k=1,MAX_LAENDERECK
                       if (send_to_pes(k,el) == ep) then
                          exit  ! already collected
                       elseif (send_to_pes(k,el) == -1) then
                          send_to_pes(k,el) = ep
                          num_send(ep) = num_send(ep) + 1
                          exit
                       elseif (k== MAX_LAENDERECK) then
                          max_laendereck_too_small = .true.  ! Problem
                       endif
                    enddo
                    
                 end do
                 
              endif
           end if
           
        end do

     elseif (all(part(elem2D_nodes(1:4,el)) == mype)) then
        ! To these elements, all neighbours are also on PE mype.
        ! But maybe, they have to be sent! Less checks are necessary, as
        ! mype is the only owner of mype.
        
        q = 4
        if (elem2D_nodes(1,el) == elem2D_nodes(4,el)) q=3  ! triangle

        do j = 1,q

           elem = elem_neighbors(j,el)

           if (elem < 1) cycle  ! boundary, "ghost element"
           
           ! 1. who owns elem? We must check all nodes!
           ! For triangles, we check elem2D_nodes(1,el)=elem2D_nodes(4,el) twice, this does not matter
           do i=1,4
              ep = part(elem2D_nodes(i,elem))
              if (ep==mype) cycle
              
              ! 2. Ok, for the owner ep, check if sending el is already collected
              do k=1,MAX_LAENDERECK
                 if (send_to_pes(k,el) == ep) then
                    exit  ! already collected
                 elseif (send_to_pes(k,el) == -1) then
                    send_to_pes(k,el) = ep
                    num_send(ep) = num_send(ep) + 1
                    exit
                 elseif (k== MAX_LAENDERECK) then
                    max_laendereck_too_small = .true.  ! Problem
                 endif
              enddo
              
           end do
                 
        end do
        
     end if
  enddo

  if (max_laendereck_too_small) then
     print *,'Set MAX_LAENDERECK in fv_comm.f90 to a larger value and recompile'
     stop
  endif
    
! Now, build the send and recv communication data structure
  com_elem2D%rPEnum = count(num_recv(0:npes-1) > 0)
  com_elem2D%sPEnum = count(num_send(0:npes-1) > 0)

  allocate(com_elem2D%rPE(com_elem2D%rPEnum))
  allocate(com_elem2D%sPE(com_elem2D%sPEnum))

  r_count = 0
  s_count = 0
  do np = 0, npes-1
     if(num_recv(np) /= 0) then
        r_count = r_count+1
        com_elem2D%rPE(r_count) = np
     end if
     if(num_send(np) /= 0) then
        s_count = s_count+1
        com_elem2D%sPE(s_count) = np
     end if
  enddo

  r_count = 0
  s_count = 0
  allocate(com_elem2D%rptr(com_elem2D%rPEnum+1)) 
  allocate(com_elem2D%sptr(com_elem2D%sPEnum+1))
  
  com_elem2D%rptr(1) = 1
  com_elem2D%sptr(1) = 1
  
  do r_count = 1, com_elem2D%rPEnum
     np = com_elem2D%rPE(r_count)
     com_elem2D%rptr(r_count+1) =  com_elem2D%rptr(r_count)+ num_recv(np)
  enddo
  do s_count = 1, com_elem2D%sPEnum
     np = com_elem2D%sPE(s_count)
     com_elem2D%sptr(s_count+1) =  com_elem2D%sptr(s_count)+ num_send(np)
  enddo

  ! Lists itself

  r_count = 0
  allocate(com_elem2D%rlist(com_elem2D%rptr(com_elem2D%rPEnum+1)-1)) 
  do np = 1,com_elem2D%rPEnum
     prank = com_elem2D%rPE(np)
     do el = 1, elem2D
        if (recv_from_pe(el) == prank) then
           r_count = r_count+1
           com_elem2D%rlist(r_count) = el
        end if
     end do
  end do
  
  s_count = 0
  allocate(com_elem2D%slist(com_elem2D%sptr(com_elem2D%sPEnum+1)-1)) 
  do np = 1,com_elem2D%sPEnum
     prank = com_elem2D%sPE(np)
     do el = 1, elem2D
        if( any(send_to_pes(:,el) == prank)) then 
           s_count = s_count+1
           com_elem2D%slist(s_count) = el
        end if
     end do
  end do



  
  !===========================================
  !  com_elem2D_full
  !===========================================

  ! The halo relations that are already determined can be kept.
  ! Just add the elements connected only via nodes.
!  num_send(0:npes-1) = 0
!  num_recv(0:npes-1) = 0
!  recv_from_pe(1:elem2D) = -1
!  send_to_pes(1:MAX_LAENDERECK,1:elem2D) = -1
  
  ! For the local elements, collect all adjacent (sharing a node) elements
  ! that belong to other PEs
  
  do el=1,elem2D
     if (any(part(elem2D_nodes(1:4,el)) == mype) .and. &
         any(part(elem2D_nodes(1:4,el)) /= mype))  then

        q = 4
        if (elem2D_nodes(1,el) == elem2D_nodes(4,el)) q=3  ! triangle

        do n = 1,q                       ! cycle through all nodes
           
           nod = elem2D_nodes(n,el)

           ! my node? Then, the whole patch is already mine.
           if (part(nod)==mype) cycle
           
           do j = 1, nod_in_elem2D_num(nod)  ! and for each node, through its patch
              elem = nod_in_elem2D(j,nod)
           
              if (all(part(elem2D_nodes(1:4,elem)) /= mype)) then
                 
                 if (recv_from_pe(elem)==-1) then  ! elem to be received already collected?
                    ! We have to receive elem from PE ep:
                    ep = part(elem2D_nodes(1,elem))  ! PE of first node is "main" owner
                    num_recv(ep) = num_recv(ep) + 1
                    recv_from_pe(elem) = ep  
                 endif
              endif
              
                 ! And maybe, we have to send el to the owners of elem
                 ! This gets more complicated:
                 ! 1. Is mype the main owner of el?
                 if (part(elem2D_nodes(1,el)) == mype) then
                    
                    ! 2. who owns elem and needs to get el? We must check all nodes!
                    p=4
                    if (elem2D_nodes(1,elem) == elem2D_nodes(4,elem)) p=3  ! triangle
                    
                    do i=1,p
                       ep = part(elem2D_nodes(i,elem))
                       
                       ! 3. Is ep also an owner of el and no send is needed? This excludes also mype==ep.
                       if (any(part(elem2D_nodes(1:4,el)) == ep)) cycle

                       ! 4. Ok, for the owner ep, check if sending el is already collected
                       do k=1,MAX_LAENDERECK
                          if (send_to_pes(k,el) == ep) then
                             exit  ! already collected
                          elseif (send_to_pes(k,el) == -1) then
                             send_to_pes(k,el) = ep
                             num_send(ep) = num_send(ep) + 1
                             exit
                          elseif (k== MAX_LAENDERECK) then
                             max_laendereck_too_small = .true.  ! Problem
                          endif
                       enddo
                       
                    end do
                 
                 endif
              end do
        end do

     elseif (all(part(elem2D_nodes(1:4,el)) == mype)) then
        ! To these elements, all neighbours are also on PE mype.
        ! But maybe, they have to be sent! Less checks are necessary, as
        ! mype is the only owner of mype.

        q = 4
        if (elem2D_nodes(1,el) == elem2D_nodes(4,el)) q=3  ! triangle

        do n = 1,q                       ! cycle through all nodes
           
           nod = elem2D_nodes(n,el)
           
           do j = 1, nod_in_elem2D_num(nod)  ! and for each node, through its patch
              elem = nod_in_elem2D(j,nod)
                                                                 
              ! 2. who owns elem and needs to get el? We must check all nodes!
              p=4
              if (elem2D_nodes(1,elem) == elem2D_nodes(4,elem)) p=3  ! triangle
                    
              do i=1,p
                 ep = part(elem2D_nodes(i,elem))
                       
                 ! 3. Is ep also an owner of el and no send is needed? This excludes also mype==ep.
                 if (any(part(elem2D_nodes(1:4,el)) == ep)) cycle

                 ! 4. Ok, for the owner ep, check if sending el is already collected
                 do k=1,MAX_LAENDERECK
                    if (send_to_pes(k,el) == ep) then
                       exit  ! already collected
                    elseif (send_to_pes(k,el) == -1) then
                       send_to_pes(k,el) = ep
                       num_send(ep) = num_send(ep) + 1
                       exit
                    elseif (k== MAX_LAENDERECK) then
                       max_laendereck_too_small = .true.  ! Problem
                    endif
                 enddo
                       
              end do
              
           end do
        end do

           
     end if
  enddo

  
  if (max_laendereck_too_small) then
     print *,'Set MAX_LAENDERECK in fv_comm.f90 to a larger value and recompile'
     stop
  endif
    
! Now, build the send and recv communication data structure
  com_elem2D_full%rPEnum = count(num_recv(0:npes-1) > 0)
  com_elem2D_full%sPEnum = count(num_send(0:npes-1) > 0)

  allocate(com_elem2D_full%rPE(com_elem2D_full%rPEnum))
  allocate(com_elem2D_full%sPE(com_elem2D_full%sPEnum))

  r_count = 0
  s_count = 0
  do np = 0, npes-1
     if(num_recv(np) /= 0) then
        r_count = r_count+1
        com_elem2D_full%rPE(r_count) = np
     end if
     if(num_send(np) /= 0) then
        s_count = s_count+1
        com_elem2D_full%sPE(s_count) = np
     end if
  enddo

  r_count = 0
  s_count = 0
  allocate(com_elem2D_full%rptr(com_elem2D_full%rPEnum+1)) 
  allocate(com_elem2D_full%sptr(com_elem2D_full%sPEnum+1))
  
  com_elem2D_full%rptr(1) = 1
  com_elem2D_full%sptr(1) = 1
  
  do r_count = 1, com_elem2D_full%rPEnum
     np = com_elem2D_full%rPE(r_count)
     com_elem2D_full%rptr(r_count+1) =  com_elem2D_full%rptr(r_count)+ num_recv(np)
  enddo
  do s_count = 1, com_elem2D_full%sPEnum
     np = com_elem2D_full%sPE(s_count)
     com_elem2D_full%sptr(s_count+1) =  com_elem2D_full%sptr(s_count)+ num_send(np)
  enddo

  ! Lists itself

  r_count = 0
  allocate(com_elem2D_full%rlist(com_elem2D_full%rptr(com_elem2D_full%rPEnum+1)-1)) 
  do np = 1,com_elem2D_full%rPEnum
     prank = com_elem2D_full%rPE(np)
     do el = 1, elem2D
        if (recv_from_pe(el) == prank) then
           r_count = r_count+1
           com_elem2D_full%rlist(r_count) = el
        end if
     end do
  end do
  
  s_count = 0
  allocate(com_elem2D_full%slist(com_elem2D_full%sptr(com_elem2D_full%sPEnum+1)-1)) 
  do np = 1,com_elem2D_full%sPEnum
     prank = com_elem2D_full%sPE(np)
     do el = 1, elem2D
        if( any(send_to_pes(:,el) == prank)) then 
           s_count = s_count+1
           com_elem2D_full%slist(s_count) = el
        end if
     end do
  end do


  ! mype sends its data to
  ! comm_elem2D_full%sPEnum external PEs
  ! Their ranks (numbers) are in array
  ! comm_elem2D_full%sPE(:)
  ! Pointers to external node numbers are in
  ! comm_elem2D_full%sptr(:)
  ! The node numbers are in 
  ! comm_elem2D_full%list(:)

end subroutine communication_elemn

!========================================================================
subroutine mymesh
  use o_MESH
  use g_PARSUP 
  implicit none
  integer                :: n, count, q, k, j, el, ed
  integer :: aux(com_elem2D_full%rptr(com_elem2D_full%rPEnum+1)-1)
  
  !======= NODES 

  ! Owned nodes + external nodes which I need:
  
  count=0   
  do n=1, nod2D
     if (part(n)==mype) then
        count=count+1
        myList_nod2D(count)=n
     end if
  end do
  myDim_nod2D = count
   eDim_nod2D = com_nod2D%rptr(com_nod2D%rPEnum+1)-1   
  myList_nod2D(myDim_nod2D+1:myDim_nod2D+eDim_nod2D) = com_nod2D%rlist(1:eDim_nod2D)

  ! Summary:  	     
  ! myList_nod2D(myDim_nod2D+1:myDim_nod2D+eDim_nod2D)
  ! contains external nodes which mype needs;    
  ! myList_nod2D(1:myDim_nod2D) contains owned nodes

  !======= ELEMENTS
  ! 2D elements 
  ! Element belongs to PE if any of its nodes is owned by PE
  ! Element is external if it is a neighbor and does not contain owned nodes
  ! The external part is needed for FVCOM type of discretization.
  
  count=0
  do el = 1, elem2D
     if( any(part(elem2D_nodes(1:4,el)) == mype) ) then
        count=count+1
        myList_elem2D(count) = el
     endif
  end do
  myDim_elem2D = count
  eDim_elem2D  = com_elem2D%rptr(com_elem2D%rPEnum+1)-1   

  ! =======
  ! full element neighbourhood requires 
  ! a longer list     
  ! ======= 

  aux=0
  do n = 1,com_elem2D%rptr(com_elem2D%rPEnum+1)-1         
     k = com_elem2D%rlist(n)
     do q = 1,com_elem2D_full%rptr(com_elem2D_full%rPEnum+1)-1   

       if (com_elem2D_full%rlist(q)==k) aux(q) = 1

     end do
  end do
  ! Test: 
  if(sum(aux) /= eDim_elem2D) write(*,*) 'mymesh problem, sum(aux)/=eDim_elem2D',sum(aux),eDim_elem2D
!  print *,'plain',com_elem2D%rptr(com_elem2D%rPEnum+1)-1 , com_elem2D%rlist(:)
!  print *,'full',com_elem2D_full%rptr(com_elem2D_full%rPEnum+1)-1,  com_elem2D_full%rlist(:)
  
  
 
  myList_elem2D(myDim_elem2D+1 : myDim_elem2D+eDim_elem2D) = com_elem2D%rlist(:)
  
  count = 0
  
  do q=1,com_elem2D_full%rptr(com_elem2D_full%rPEnum+1)-1   
     if (aux(q)==0) then
        count=count+1
        myList_elem2D( myDim_elem2D+eDim_elem2D + count)=com_elem2D_full%rlist(q)
     end if
  end do
  eXDim_elem2D = count 

              
  ! Summary: 
  ! myList_elem2D(1:myDim_elem2D) contains elements with at least single owned node.
  ! myList_elem2D(myDim_elem2D+1:myDim_elem2D+eDim_elem2D) contains elements-neighbours.
  ! They are needed when interpolation is done in FV type of code.
  ! The final piece from myDim_elem2D+eDim_elem2D to 
  ! myDim_elem2D+eDim_elem2D+eXDim_elem2D is needed for MUSCL-type advection
end subroutine mymesh
!=============================================================================
subroutine save_partitioning

USE o_MESH
USE g_PARSUP 
IMPLICIT NONE

integer      :: n, ncount(0:npes), temp(nod2D)
 
 open(21, file=trim(DistPath)//'rpart.out')

 
  ! ================================
 ! mapping ( PE contiguous 2D numbering) 	 
 ! ================================  
 
 ncount=0
 DO n=1, nod2D
    ncount(part(n))=ncount(part(n))+1
    temp(n) = ncount(part(n))
 END DO
 write(21,*) npes
 write(21,*) ncount(0:npes-1)
    
!!$ ncount=0
!!$ DO n=1, nod2D
!!$    ncount(part(n)) = ncount(part(n))+1
!!$    temp(n) = ncount(part(n))
!!$ END DO

 ncount(0) = 0
 Do n=1,npes	  
    ncount(n) = ncount(n)+ncount(n-1)
 end do
 ! Now count == part in range partitioning   
    
 do n=1,nod2D
    temp(n) = temp(n)+ncount(part(n))
    write(21,*) temp(n)  
 end do
 close(21)

 ! ==============================
 ! The domain decomposition itself
 ! ==============================
 open(22,file=trim(DistPath)//'partitioning.out')
 do n=1,nod2D
    write(22,*) part(n)
 enddo
 close(22)
end subroutine save_partitioning
!=============================================================================
SUBROUTINE save_dist_mesh

USE o_MESH
USE g_PARSUP 
IMPLICIT NONE

 Integer      :: fileID
 character*10 :: mype_string
 character*80 :: file_name

 write(mype_string,'(i5.5)') mype      
 file_name=trim(DistPath)//'my_list'//trim(mype_string)//'.out'  
 fileID = 10+mype  
 open(fileID, file=file_name)

 ! =============================   
 ! lists of owned nodes, elements 
 ! =============================
 write(fileID,*) mype
 write(fileID,*) myDim_nod2D
 write(fileID,*) eDim_nod2D 	 
 write(fileID,*) myList_nod2D(1:myDim_nod2D+eDim_nod2D)

 write(fileID,*) myDim_elem2D
 write(fileID,*) eDim_elem2D
 write(fileID,*) eXDim_elem2D	 	 
 write(fileID,*) myList_elem2D(1:myDim_elem2D +eDim_elem2D +eXDim_elem2D)
 close(fileID)

 close(fileID)
   
      ! =========================  
      ! communication information for nodes
      ! ========================= 
 call com_global2local_nodes   

 file_name=trim(DistPath)//'com_info'//trim(mype_string)//'.out' 
 fileID = npes+10+mype  
 open(fileID, file=file_name)
 write(fileID,*) mype
 write(fileID,*) com_nod2D%rPEnum
 write(fileID,*) com_nod2D%rPE
 write(fileID,*) com_nod2D%rptr
 write(fileID,*) com_nod2D%rlist
 write(fileID,*) com_nod2D%sPEnum
 write(fileID,*) com_nod2D%sPE
 write(fileID,*) com_nod2D%sptr
 write(fileID,*) com_nod2D%slist



 deallocate(com_nod2D%rPE)
 deallocate(com_nod2D%rptr)
 deallocate(com_nod2D%rlist)
 deallocate(com_nod2D%sPE)
 deallocate(com_nod2D%sptr)
 deallocate(com_nod2D%slist)

   
      ! =========================  
      ! communication information for elements
      ! ========================= 
 call com_global2local_elements

 write(fileID,*) com_elem2D%rPEnum
 write(fileID,*) com_elem2D%rPE
 write(fileID,*) com_elem2D%rptr
 write(fileID,*) com_elem2D%rlist
 write(fileID,*) com_elem2D%sPEnum
 write(fileID,*) com_elem2D%sPE
 write(fileID,*) com_elem2D%sptr
 write(fileID,*) com_elem2D%slist

 write(fileID,*) com_elem2D_full%rPEnum
 write(fileID,*) com_elem2D_full%rPE
 write(fileID,*) com_elem2D_full%rptr
 write(fileID,*) com_elem2D_full%rlist
 write(fileID,*) com_elem2D_full%sPEnum
 write(fileID,*) com_elem2D_full%sPE
 write(fileID,*) com_elem2D_full%sptr
 write(fileID,*) com_elem2D_full%slist
 close(fileID)

 deallocate(com_elem2D%rPE)
 deallocate(com_elem2D%rptr)
 deallocate(com_elem2D%rlist)
 deallocate(com_elem2D%sPE)
 deallocate(com_elem2D%sptr)
 deallocate(com_elem2D%slist)
 
 deallocate(com_elem2D_full%rPE)
 deallocate(com_elem2D_full%rptr)
 deallocate(com_elem2D_full%rlist)
 deallocate(com_elem2D_full%sPE)
 deallocate(com_elem2D_full%sptr)
 deallocate(com_elem2D_full%slist)


  write(*,*) 'Distributed mesh is saved for mype=',mype
 
END subroutine  save_dist_mesh

!===========================================================================
SUBROUTINE com_global2local_nodes
USE g_PARSUP
USE o_MESH
IMPLICIT NONE
INTEGER   :: n, m
INTEGER   :: temp(nod2D) 
! =========
! nodes
! =========
  ! Replace global numbering with a local one
  temp(1:nod2D)=0
  DO n=1, myDim_nod2D+eDim_nod2D
     temp(myList_nod2D(n))=n
  END DO
  DO n=1, com_nod2D%sptr(com_nod2D%sPEnum+1)-1

	 com_nod2D%slist(n) = temp(com_nod2D%slist(n))
  END DO

  DO n=1, com_nod2D%rptr(com_nod2D%rPEnum+1)-1
	 com_nod2D%rlist(n) = temp(com_nod2D%rlist(n))
  END DO
	 ! com_nod2D%rlist should be  
	 ! myDim_nod2D+1:myDim_nod2D+eDim_nod2D
end SUBROUTINE com_global2local_nodes

!===========================================================================
SUBROUTINE com_global2local_elements
USE g_PARSUP
USE o_MESH
IMPLICIT NONE
INTEGER   :: n
INTEGER   :: temp(elem2D) 
! =========
! elements
! =========
  ! Replace global numbering with a local one
  temp(1:elem2D)=0
  DO n=1, myDim_elem2D+eDim_elem2D+eXDim_elem2D
     temp(myList_elem2D(n))=n
  END DO
  DO n=1, com_elem2D%sptr(com_elem2D%sPEnum+1)-1

	 com_elem2D%slist(n) = temp(com_elem2D%slist(n))
  END DO

  DO n=1, com_elem2D%rptr(com_elem2D%rPEnum+1)-1
         
	 com_elem2D%rlist(n) = temp(com_elem2D%rlist(n))
  END DO

! =========
! elements (extra needed for MUSCL advection)
! =========

  DO n=1, com_elem2D_full%sptr(com_elem2D_full%sPEnum+1)-1

	 com_elem2D_full%slist(n) = temp(com_elem2D_full%slist(n))
  END DO

  DO n=1, com_elem2D_full%rptr(com_elem2D_full%rPEnum+1)-1
        
	 com_elem2D_full%rlist(n) = temp(com_elem2D_full%rlist(n))
  END DO
  
	 ! com_elem2%rlist should be 
	 ! myDim_elem2D+1:myDim_elem2D+eDim_elem2D
end SUBROUTINE com_global2local_elements
 

!========================================================================


!========================================================================
!NR Edges - maybe, the halos are not needed at all. Skip any work here
!NR until we are sure that it is needed.
!========================================================================

! subroutine communication_edgen

!   use o_MESH
!   use g_PARSUP 
!   implicit none

  
!   integer, parameter   :: MAX_NEIGH_PART=64
!   integer              :: j,n,np, nz, nu,prank, el, eledges(4), count
!   integer              :: q, ed
!   integer              :: aux(MAX_NEIGH_PART, edge2D)

!   integer :: pmap(0:npes-1), pattern(0:npes-1,0:npes-1), pnum(0:npes-1,0:npes-1)

!   ! Assume we have 2D partitioning vector in part. Find communication
!   ! rules. An edge is external to a node n if the triple
!   ! n, edge_nodes(:,edge) forms an element while both PEs of edge 
!   ! are different from PE of n

!   ! an edge is owned by PE if any its nodes belongs to PE 

!   allocate(myList_edge2D(edge2D))

!   pattern(:,:)= 0
!   pmap(:)  = 0

!   do el=1,elem2D

!      if ( all(part(elem2D_nodes(1:3,el)) == part(elem2D_nodes(4,el))) ) &
!           cycle    ! all nodes of this element on the same partition

!      ! if both PEs of an edge opposing a node
!      ! are different from node's PE then the edge should be communicated
!      q=4 ! quad
!      if (elem2D_nodes(1,el) == elem2D_nodes(4,el)) q=3 !triangle

!      eledges(1:q) = elem_edges(1:q,el)
!      do j=1,q
!         n=elem2D_nodes(j,el)
!         ! find edges that do not contain n

!         do nz=1,q
! 	   if((edge_nodes(1,eledges(nz)) /= n) .and. &
!               (edge_nodes(2,eledges(nz)) /= n)) then

!               if((part(n) /= part(edge_nodes(1,eledges(nz)))) .and.&   
!                  (part(n) /= part(edge_nodes(2,eledges(nz))))) then 

!                  pattern(part(n), part(edge_nodes(1,eledges(nz)))) = 1
!                  pattern(part(n), part(edge_nodes(2,eledges(nz)))) = 1
!               end if

! 	   end if
!         end do
!      end do
!   end do

!   do mype=0,npes-1
!   count=0
!   DO np=0,npes-1
!      if(pattern(mype,np) /= 0  .or.  np == mype ) then
!         count    = count+1
!         pmap(np) = count
!      end if
!   END DO

!   if (count > MAX_NEIGH_PART) then
!      print *,'Set MAX_NEIGH_PART to at least',count,' and recompile'
!      stop
!   endif
!   ! This has a much smaller size if npes is large
!   ! but we need pmap to address it
!   aux(:,:)  = 0
  
!   do el=1,elem2D

!      if ( all(part(elem2D_nodes(1:3,el)) == part(elem2D_nodes(4,el))) ) &
!           cycle    ! all nodes of this element on the same partition

!      ! if both PEs of an edge opposing a node
!      ! are different from node's PE then the edge should be communicated
!      q=4 ! quad
!      if (elem2D_nodes(1,el) == elem2D_nodes(4,el)) q=3 !triangle

!      eledges(1:q) = elem_edges(1:q,el)

!      do j=1,q
!         n = elem2D_nodes(j,el)
!         ! find edges that do not contain n
!         if (pmap(part(n))/=0) then
!            do nz=1,q
!               if((edge_nodes(1,eledges(nz)) /= n) .and. &
!                  (edge_nodes(2,eledges(nz)) /= n)) then

!                  if((part(n) /= part(edge_nodes(1,eledges(nz)))) .and.&   
!                     (part(n) /= part(edge_nodes(2,eledges(nz))))) then

!                     aux(pmap(part(n)),eledges(nz))=1
!                  end if

!               end if
!            end do
!         end if
!      end do
!   end do

!   pnum(:,:) = 0
!   do ed=1, edge2D
!      do np = 0, npes-1
!         if(pmap(np) /= 0) then
!            if(aux(pmap(np),ed) /= 0) pnum(np,part(edge_nodes(1,ed))) = pnum(np,part(edge_nodes(1,ed)))+1
! 	end if
!      end do
!   end do

!   ! We know how many external edges each PE needs
!   ! This is the 'receive' list   
!   ! com_edge2D for 2D nodes
!   ! 

!   ! The number of external PEs I receive information from
!   com_edge2D%rPEnum=0
!   ! SENDING PART
!   com_edge2D%sPEnum=0

!   do np=0, npes-1
!      if(pnum(mype,np) /= 0) com_edge2D%rPEnum = com_edge2D%rPEnum+1
!      if(pnum(np,mype) /= 0) com_edge2D%sPEnum = com_edge2D%sPEnum+1
!   end do

!   ! Their ranks (PE numbers)

!   count=0
!   allocate(com_edge2D%rPE(com_edge2D%rPEnum))
!   do np = 0, npes-1
!      if(pnum(mype,np) /= 0) then
!         count = count+1
!         com_edge2D%rPE(count) = np
!      end if
!   end do


!   ! Ptr to list of external nodes ordered by external PE ranks

!   count=0
!   allocate(com_edge2D%rptr(com_edge2D%rPEnum+1)) 
!   com_edge2D%rptr(1) = 1
!   do np=0, npes-1
!      if(pnum(mype,np) /= 0) then
!         count = count+1
!         com_edge2D%rptr(count+1) = com_edge2D%rptr(count)+ pnum(mype,np)
!      end if
!   end do

!   ! List itself

!   count=0
!   allocate(com_edge2D%rlist(com_edge2D%rptr(com_edge2D%rPEnum+1)-1)) 
!   do np=1,com_edge2D%rPEnum
!      prank=com_edge2D%rPE(np)
!      do ed=1, edge2D
!         if((aux(pmap(mype),ed)==1).and.(part(edge_nodes(1,ed))==prank)) then
!            count=count+1
!            com_edge2D%rlist(count)=ed
!         end if
!      end do
!   end do

!   ! Summary of this piece: mype receives
!   ! information on external 2D edges from
!   ! comm_edge2D%rPEnum external PEs
!   ! Their ranks (numbers) are in array
!   ! comm_edge2D%rPE(:)
!   ! Pointers to external node numbers are in
!   ! comm_edge2D%rptr(:)
!   ! The node numbers are in 
!   ! comm_edge2D%list(:)
!   ! Putting everything into structure takes many operations, but
!   ! without the structure we will need to many names and arrays
!   ! Do not forget that we need also send part, and we need analogous
!   ! information for 3D nodes and edges.

!   ! Their ranks (PE numbers)

!   count=0
!   allocate(com_edge2D%sPE(com_edge2D%sPEnum))
!   do np = 0, npes-1
!      if(pnum(np,mype) /= 0) then
!         count = count+1
!         com_edge2D%sPE(count) = np
!      end if
!   end do

!   ! Ptr to list of external nodes ordered by external PE ranks
!   count=0
!   allocate(com_edge2D%sptr(com_edge2D%sPEnum+1)) 
!   com_edge2D%sptr=1
!   do np = 0, npes-1
!      if(pnum(np,mype) /= 0) then
!         count=count+1
!         com_edge2D%sptr(count+1) = com_edge2D%sptr(count)+ pnum(np,mype)
!      end if
!   end do

!   ! List itself

!   count=0
!   allocate(com_edge2D%slist(com_edge2D%sptr(com_edge2D%sPEnum+1)-1)) 
!   do np=1,com_edge2D%sPEnum
!      prank=com_edge2D%sPE(np)
!      if(pmap(prank) /= 0) then
!         do ed=1, edge2D
!            if((aux(pmap(prank),ed)==1) .and. (part(edge_nodes(1,ed))==mype)) then
!               count=count+1
!               com_edge2D%slist(count)=ed
!            end if
!         end do
!      end if
!   end do

!   ! mype sends its data to
!   ! comm_edge2D%sPEnum external PEs
!   ! Their ranks (numbers) are in array
!   ! comm_edge2D%sPE(:)
!   ! Pointers to external node numbers are in
!   ! comm_edge2D%sptr(:)
!   ! The node numbers are in 
!   ! comm_edge2D%list(:)

!   call mymesh_edges
!   call save_dist_mesh_edges
! enddo ! mype

!   deallocate(myList_edge2D)

! end subroutine communication_edgen

! !========================================================================
! subroutine mymesh_edges

!   use o_MESH
!   use g_PARSUP 
!   implicit none
!   integer                :: counter, ed
!   ! ======== EDGES 
!   ! Owned edges (both nodes are mine)+ shared edges I do computations 
!   ! at (only one node is mine; some other PE updates them
!   ! simultaneously with me but I do not care ) + 
!   ! external edges which I need (neither of nodes is mine, but they
!   ! belong to elements in myList_elem:
!   counter=0
!   do ed = 1, edge2D 
!      if (part(edge_nodes(1,ed))==mype .or. part(edge_nodes(1,ed))==mype) then
        
!         counter=counter+1
!         myList_edge2D(counter) = ed
!      endif
!   end do
!   myDim_edge2D = counter
!   eDim_edge2D  = com_edge2D%rptr(com_edge2D%rPEnum+1)-1   
  
!   myList_edge2D(myDim_edge2D+1 : myDim_edge2D+eDim_edge2D) = com_edge2D%rlist(1:eDim_edge2D)

!   ! Summary:  	     
!   ! myList_edge2D(myDim_edge2D+1:myDim_edge2D+eDim_edge2D)
!   ! contains external edges which mype needs;    
!   ! myList_edge2D(1:myDim_edge2D) contains owned edges +
!   ! shared edges which mype updates

! end subroutine mymesh_edges

! !=============================================================================
! SUBROUTINE save_dist_mesh_edges

! USE o_MESH
! USE g_PARSUP 
! IMPLICIT NONE
!  Integer      :: fileID
!  character*10 :: mype_string
!  character*80 :: file_name

!  write(mype_string,'(i5.5)') mype      
!  file_name=trim(DistPath)//'my_list'//trim(mype_string)//'.out'  
!  fileID=10+mype  
!  open(fileID, file=file_name,position='append')
!  ! =============================   
!  ! lists of owned nodes, elements, edges 
!  ! =============================
 
!  write(fileID,*) myDim_edge2D
!  write(fileID,*) eDim_edge2D 	 
!  write(fileID,*) myList_edge2D(1:myDim_edge2D +eDim_edge2D)
!  close(fileID)       
 
   
!       ! =========================  
!       ! communication information
!       ! ========================= 
!  call com_global2local_edges   
 
!  file_name=trim(DistPath)//'com_info'//trim(mype_string)//'.out' 
!  fileID = npes+10+mype  
!  open(fileID, file=file_name,position='append')
 
!  write(fileID,*) com_edge2D%rPEnum
!  write(fileID,*) com_edge2D%rPE
!  write(fileID,*) com_edge2D%rptr
!  write(fileID,*) com_edge2D%rlist
!  write(fileID,*) com_edge2D%sPEnum
!  write(fileID,*) com_edge2D%sPE
!  write(fileID,*) com_edge2D%sptr
!  write(fileID,*) com_edge2D%slist
!  close(fileID)

!  deallocate(com_edge2D%rPE)
!  deallocate(com_edge2D%rptr)
!  deallocate(com_edge2D%rlist)
!  deallocate(com_edge2D%sPE)
!  deallocate(com_edge2D%sptr)
!  deallocate(com_edge2D%slist)
 	 

!  write(*,*) 'Distributed mesh is saved for mype=',mype
! END subroutine  save_dist_mesh_edges
!===========================================================================
! SUBROUTINE com_global2local_edges
! USE g_PARSUP
! USE o_MESH
! IMPLICIT NONE
! INTEGER   :: n
! INTEGER   :: temp(edge2D) 

! ! =========
! ! edges
! ! =========
!   ! Replace global numbering with a local one
!   temp(1:edge2D)=0
!   DO n = 1, myDim_edge2D+eDim_edge2D
!      temp(myList_edge2D(n)) = n
!   END DO
!   DO n=1, com_edge2D%sptr(com_edge2D%sPEnum+1)-1

! 	 com_edge2D%slist(n) = temp(com_edge2D%slist(n))
!   END DO

!   DO n=1, com_edge2D%rptr(com_edge2D%rPEnum+1)-1

! 	 com_edge2D%rlist(n) = temp(com_edge2D%rlist(n))
!   END DO
! 	 ! com_edge2%rlist should be 
! 	 ! myDim_edge2D+1:myDim_edge2D+eDim_edge2D
! end SUBROUTINE com_global2local_edges
