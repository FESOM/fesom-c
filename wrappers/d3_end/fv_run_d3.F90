program run_d3

    use kernel_setup
    use d3_module
    use arrays_io
    
    implicit none

    !character(len=100) :: filename
    integer :: stepindx

    call initialization_d3()


    write(*,*) "Have fun ... "

    ! index of file from d3_end_dump/
    ! stepindx is a model time step
    stepindx = 263000 
    ! reading arrays that in use by d3_end (model output)
    call read_arrays_in(stepindx)
    ! reading arrays calculated by d3_end (model output)
    call read_arrays_out(stepindx)

    ! run d3_end
    call d3_end

    ! compare local d3_end with model output
    call compare_arrays    




end program run_d3