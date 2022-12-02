program test_tmp
    use mod_mmpol
    use mod_inputloader

    implicit none
    
    type(ommp_system), allocatable :: mysys

    allocate(mysys)
    call mmpol_init_from_mmp("/home/mattia/LibEnv/open-mmpol/tests/N-methylacetamide/input_AMOEBA.mmp", mysys)
    deallocate(mysys)
end program

