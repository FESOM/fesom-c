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

end module o_ARRAYS

module o_MESH
    use o_PARAM

    real(kind=WP),allocatable,dimension(:)      :: depth   !depth (defined on every node)
    real(kind=WP), Allocatable, dimension(:)    :: sigma   ! vertical levels (from 0 to 1) 

    integer                                     :: myDim_elem2D!, eDim_elem2D, eXDim_elem2D
    integer, allocatable, dimension(:,:)        :: elem2D_nodes

    real(kind=WP), allocatable, dimension(:,:)  :: w_cv


    

end module o_MESH    

module g_PARSUP !from gen_modules_partitioning_c

    use o_PARAM

    integer                       :: myDim_nod2D, eDim_nod2D




end module g_PARSUP    