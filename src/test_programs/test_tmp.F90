program test_tmp
    use mod_mmpol

    implicit none
    
    type(ommp_system), allocatable :: mysys

    allocate(mysys)
    call mmpol_init(mysys, 1, 10, 10)
    call mmpol_terminate(mysys)
    deallocate(mysys)
end program

