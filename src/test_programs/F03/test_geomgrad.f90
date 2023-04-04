program test_geomgrad
    use iso_c_binding, only: c_char
    use iso_fortran_env, only: error_unit
    use ommp_interface

    implicit none

    integer :: argc, i, retcode=0, io_file=101, mm_atoms
    character(kind=c_char, len=120), dimension(3) :: argv
    real(ommp_real), allocatable :: grad_ana(:,:)
    type(ommp_system), pointer :: my_system

    argc = command_argument_count()
    if(argc /= 2) then
        write(*, *) "Syntax expected"
        write(*, *) "    $ test_init.exe <INPUT FILE> <OUTPUT FILE>"
        stop 1
    end if

    do i=1, argc
        call get_command_argument(i, argv(i))
    end do

    call ommp_init_mmp(my_system, argv(1))
    mm_atoms = my_system%top%mm_atoms

    allocate(grad_ana(3, my_system%top%mm_atoms))
    open(io_file, file=argv(2))

    grad_ana = 0.0
    call ommp_fixedelec_geomgrad(my_system, grad_ana)

    write(*, *) "Grad EM"
    do i=1, mm_atoms
        write(io_file, "(3f12.8)")grad_ana(:,i)
    end do
    
    grad_ana = 0.0
    call ommp_polelec_geomgrad(my_system, grad_ana)

    write(*, *) "Grad EP"
    do i=1, mm_atoms
        write(io_file, "(3f12.8)")grad_ana(:,i)
    end do
    
    deallocate(grad_ana)
    close(io_file)

    call ommp_terminate(my_system)
    
    if(retcode > 0) then
        stop 1
    else
        stop 0
    end if
end program
