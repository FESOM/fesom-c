module o_PARAM ! from fv_var.F90

    integer, parameter            :: WP=8   ! Working precision
    real(kind=WP), parameter      :: g=9.8_WP   

    integer                       :: nsigma ! number of layers (vertical sigma layers)
    real(kind=WP)                 :: dt     ! time step for baroclinic mode
    real(kind=WP)                 :: cka = 0.4_WP !Karman constant

    real(kind=WP)                 :: Dmin  ! Minimal depth

    real(kind=WP)                 :: z0b_min      ! bottom roughness 
    real(kind=WP)                 :: z0s_min      ! surface roughness 
    real(kind=WP)                 :: beta_scale=3.9_WP      ! 0<=beta_scale<=4  Cut off function for scale
    real(kind=WP)                 :: snul ! minimum diffusion
    real(kind=WP)                 :: PR_num=0.0001_WP          ! Prandtl number ratio of momentum diffusivity

    real(kind=WP)                 :: density_0 ! zero density 

    logical                       :: comp_sediment
    real(kind=WP)                 :: plop_s     ! plop_s - density of individual grains [kg/m3]

end module o_PARAM

module o_ARRAYS ! from fv_var.F90

    use o_PARAM

    real(kind=WP), allocatable    :: bt(:,:), snu(:,:), Ri(:,:), tke_dissip(:,:)  !turbulence characteristics (energy, Av, Richardson, dissipation)
    real(kind=WP), allocatable    :: Kv(:,:), Av(:,:)!, Av_node(:,:)

    real(kind=WP), allocatable    :: eta_n(:) ! sea level elevation (SSH)

    real(kind=WP), allocatable    :: Unode(:,:), Vnode(:,:) ! velocity on nodes
    real(kind=WP), allocatable    :: TF(:,:)  ! temperature
    real(kind=WP), Allocatable    :: SF(:,:)  ! salinity

    real(kind=WP), allocatable    :: zbar(:,:)! node depth

    real(KIND=WP), allocatable    :: windx(:), windy(:)!, wind(:) ! wind at 10 m in X and Y directions and abs

    real(kind=WP), allocatable    :: Jc(:,:)

    !Arrays to check
    real(kind=WP), allocatable    :: bt_out(:,:), snu_out(:,:)
    real(kind=WP), allocatable    :: Kv_out(:,:), Av_out(:,:)

end module o_ARRAYS

module o_MESH
    use o_PARAM

    real(kind=WP),allocatable,dimension(:)      :: depth   !depth (defined on every node)
    real(kind=WP), Allocatable, dimension(:)    :: sigma   ! vertical levels (from 0 to 1) 

    integer                                     :: myDim_elem2D, eDim_elem2D, eXDim_elem2D
    integer, allocatable, dimension(:,:)        :: elem2D_nodes

    real(kind=WP), allocatable, dimension(:,:)  :: w_cv


    

end module o_MESH    

module g_PARSUP !from gen_modules_partitioning_c

    use o_PARAM

    integer                       :: myDim_nod2D, eDim_nod2D




end module g_PARSUP    

module kernel_setup

    use o_PARAM
    use o_MESH
    use o_ARRAYS
    use g_PARSUP

    implicit none

    public      :: initialization_d3

    contains


        subroutine initialization_d3()
            
            implicit none

            !character(len=*), intent(in) :: filename
            character(len=100) :: filename
            integer :: io_status

            write(*,*) "initialization_d3: Start inizialization ..."
            filename = 'd3_end_dump/d3_param.dat'
            open(unit=20, file=filename, form='unformatted', status='old', action='read', iostat=io_status)
            if (io_status /= 0) then
              print *, 'Error opening file: ',filename
              STOP
            end if
          
            read(20) nsigma, dt, cka, Dmin, z0b_min, z0s_min, beta_scale, snul, PR_num, density_0, comp_sediment, plop_s
            read(20) myDim_nod2D, eDim_nod2D, myDim_elem2D, eDim_elem2D, eXDim_elem2D
            close(20)

            write(*,*) "initialization_d3: Arrays allocation ..."
            call allocate_arrays()

            write(*,*) "initialization_d3: Read general arrays ..."
            filename = 'd3_end_dump/d3_param_arrays.dat'
            open(unit=20, file=filename, form='unformatted', status='old', action='read', iostat=io_status)
            if (io_status /= 0) then
                print *, 'Error opening file: ',filename
                STOP
            end if
          
            read(20) depth, sigma, elem2D_nodes, w_cv
            close(20)
    

         
            ! write arrays
!            filename='d3_param_arrays.dat'
!            open(unit=77, file=filename, form='unformatted', status='replace', action='write', iostat=io_status)
!            if (io_status /= 0) then
!              print *, 'Error opening file!'
!              return
!            end if
!            write(77) depth, sigma, elem2D_nodes, w_cv
!            close(77) 

        end subroutine initialization_d3

        subroutine allocate_arrays()

            implicit none

            integer :: alloc_status
            integer :: node_size, elem_size


            node_size=myDim_nod2D+eDim_nod2D
            elem_size=myDim_elem2D+eDim_elem2D+eXDim_elem2D

            allocate(depth(node_size), &
                    sigma(nsigma), &
                    elem2D_nodes(4, myDim_elem2D),  &
                    w_cv(4,elem_size),  &
                           stat=alloc_status)
            if( alloc_status /= 0 )   STOP 'allocate_arrays: allocation error'                           

            allocate(bt(nsigma,node_size), &
                    snu(nsigma,node_size), &
                    Ri(nsigma,node_size), &
                    tke_dissip(nsigma,node_size), &
                    Kv(nsigma,node_size), &
                    Av(nsigma,elem_size), &
                    eta_n(node_size), &
                    Unode(nsigma-1,node_size), &
                    Vnode(nsigma-1,node_size), &
                    TF(nsigma-1,node_size), &
                    SF(nsigma-1,node_size), &
                    zbar(nsigma,node_size), &
                    windx(node_size), &
                    windy(node_size), &
                    Jc(nsigma-1,node_size), &
                            stat=alloc_status)
            if( alloc_status /= 0 )   STOP 'allocate_arrays: allocation error'                           

            allocate(bt_out(nsigma,node_size), &
                    snu_out(nsigma,node_size), &
                    Kv_out(nsigma,node_size), &
                    Av_out(nsigma,elem_size), &
                            stat=alloc_status)
            if( alloc_status /= 0 )   STOP 'allocate_arrays: allocation error'                           





        end subroutine allocate_arrays

 
 


end module kernel_setup

module arrays_io

    use o_PARAM
    use o_MESH
    use o_ARRAYS
    use g_PARSUP

    implicit none

    public      :: read_arrays_in, read_arrays_out 

    contains

        subroutine read_arrays_in(stepindx)

            implicit none

            integer, intent(in) :: stepindx
            character(len=100)  :: filename            
            character(len=100)  :: stepstr
            integer             :: io_status

            write(stepstr, '(I10.10)')  stepindx
            write(filename,*) 'd3_end_dump/',trim(ADJUSTL(trim(stepstr))),'_arrays_in.dat'
            filename=ADJUSTL(trim(filename))
            open(unit=20, file=filename, form='unformatted', status='old', action='read', iostat=io_status)
            if (io_status /= 0) then
              print *, 'Error opening file: ',filename
              STOP
            end if
          
            read(20) bt, snu, Kv, Av, eta_n, Unode, Vnode, TF, SF, zbar, windx, windy, Jc
            close(20)

        end subroutine read_arrays_in

        subroutine read_arrays_out(stepindx)

            implicit none

            integer, intent(in) :: stepindx
            character(len=100)  :: filename          
            character(len=100)  :: stepstr              
            integer             :: io_status

            write(stepstr, '(I10.10)')  stepindx
            write(filename,*) 'd3_end_dump/',trim(ADJUSTL(trim(stepstr))),'_arrays_in.dat'
            filename=ADJUSTL(trim(filename))
            open(unit=20, file=filename, form='unformatted', status='old', action='read', iostat=io_status)
            if (io_status /= 0) then
                print *, 'Error opening file: ',filename
                STOP
            end if
          
            read(20) bt_out, snu_out, Kv_out, Av_out
            close(20)

        end subroutine read_arrays_out     
        
        subroutine compare_arrays()

            implicit none

            real(kind=WP), dimension(nsigma,myDim_nod2D+eDim_nod2D) :: array_diff
            real(kind=WP), dimension(nsigma,myDim_elem2D+eDim_elem2D+eXDim_elem2D) :: array_diff_elem

            array_diff = 0.0_WP
            array_diff_elem = 0.0_WP
            
            write(*,*) "Max/Min values for differences between in and out arrays"            
            ! (:,1:myDim_nod2D) we compare only values relevant to current CPU
            array_diff = bt - bt_out
            write(*,*) "bt:", maxval(array_diff(:,1:myDim_nod2D)), minval(array_diff(:,1:myDim_nod2D))
            array_diff = snu - snu_out
            write(*,*) "snu:", maxval(array_diff(:,1:myDim_nod2D)), minval(array_diff(:,1:myDim_nod2D))
            array_diff = Kv - Kv_out
            write(*,*) "Kv:", maxval(array_diff(:,1:myDim_nod2D)), minval(array_diff(:,1:myDim_nod2D))
            array_diff_elem = Av - Av_out
            write(*,*) "Av:", maxval(array_diff(:,1:myDim_elem2D)), minval(array_diff(:,1:myDim_elem2D))

        end subroutine compare_arrays        


end module arrays_io    
