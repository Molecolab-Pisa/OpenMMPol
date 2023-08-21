program test_energy
    use iso_c_binding, only: c_char
    use ommp_interface

    implicit none

    integer :: argc, i
    character(kind=c_char, len=120), dimension(3) :: argv
    real(ommp_real), allocatable :: ef(:,:), ef_pol(:,:), E_MMMM, E_MMPOL
    type(ommp_system), pointer :: my_system

    argc = command_argument_count()
    if(argc /= 2 .and. argc /= 3 ) then
        write(*, *) "Syntax expected"
        write(*, *) "    $ test_init.exe <INPUT FILE> <OUTPUT FILE>&
                    & [<ELECTRIC FIELD FILE>]"
        stop 1
    end if

    do i=1, argc
        call get_command_argument(i, argv(i))
    end do

    call ommp_init_mmp(my_system, argv(1))

    allocate(ef(3, my_system%top%mm_atoms))
    allocate(ef_pol(3, my_system%eel%pol_atoms))
    
    if(argc == 3) then
        open(unit=101, file=argv(3)) 
        do i=1, my_system%top%mm_atoms
            read(101, *) ef(:, i)
        end do
        close(101)
    else
        ef = 0.0
    end if

    do i=1, my_system%eel%pol_atoms
        ef_pol(:,i) = ef(:, my_system%eel%polar_mm(i))
    end do
    
    deallocate(ef)
    
    E_MMMM = ommp_get_fixedelec_energy(my_system)
    call ommp_set_external_field(my_system, ef_pol, OMMP_SOLVER_NONE, OMMP_MATV_NONE)
    E_MMPOL = ommp_get_polelec_energy(my_system)

    ! Write ipd
    open(unit=101, file=argv(2))
    write(101, "(F20.12)") E_MMMM
    write(101, "(F20.12)") E_MMPOL
    close(101)

    deallocate(ef_pol)
    call ommp_terminate(my_system)

end program test_energy
