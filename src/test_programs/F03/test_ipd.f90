program test_ipd
    use iso_c_binding, only: c_char
    use ommp_interface

    implicit none

    integer :: argc, i, j, k, solver
    character(kind=c_char, len=120), dimension(4) :: argv
    real(ommp_real), allocatable :: ef(:,:), ef_pol(:,:)
    type(ommp_system), pointer :: my_system

    argc = command_argument_count()
    if(argc /= 3 .and. argc /= 4 ) then
        write(*, *) "Syntax expected"
        write(*, *) "    $ test_init.exe <INPUT FILE> <OUTPUT FILE> <SOLVER>&
                    & [<ELECTRIC FIELD FILE>]"
        stop 1
    end if

    do i=1, argc
        call get_command_argument(i, argv(i))
    end do

    call ommp_set_verbose(OMMP_VERBOSE_DEBUG)
    call ommp_init_mmp(my_system, argv(1))

    select case(argv(3))
        case("cg")
            solver = OMMP_SOLVER_CG
        case("inversion")
            solver = OMMP_SOLVER_INVERSION
        case("diis")
            solver = OMMP_SOLVER_DIIS
        case default
            write(*, "('Unrecognized solver (',A,'). Exiting.')") argv(3)
            stop 1
    end select

    allocate(ef(3, my_system%top%mm_atoms))
    allocate(ef_pol(3, my_system%eel%pol_atoms))
    
    if(argc == 4) then
        open(unit=101, file=argv(4)) 
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
    call ommp_set_external_field(my_system, ef_pol, solver)

    ! Write ipd
    open(unit=101, file=argv(2))
    do k = 1, my_system%eel%n_ipd
        do j = 1, my_system%eel%pol_atoms
            write(101, "(3F20.12)") my_system%eel%ipd(:, j, k)
        end do
    end do
    close(101)

    deallocate(ef_pol)
    call ommp_terminate(my_system)

end program test_ipd
