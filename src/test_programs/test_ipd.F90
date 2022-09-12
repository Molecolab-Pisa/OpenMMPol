program test_ipd
    use iso_c_binding, only: c_char
    use mod_ommp_interface

    implicit none

    integer :: argc, i, j, k, solver
    character(kind=c_char, len=120), dimension(4) :: argv
    real(ommp_real), allocatable :: ef(:,:)

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

    call ommp_init_mmp(argv(1))
    call ommp_set_verbose(OMMP_VERBOSE_DEBUG)

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

    allocate(ef(3, ommp_pol_atoms))
    
    if(argc == 4) then
        open(unit=101, file=argv(4)) 
        do i=1, ommp_pol_atoms
            read(101, *) ef(:, i)
        end do
        close(101)
    else
        ef = 0.0
    end if
    
    call ommp_set_external_field(ef, solver)

    ! Write ipd
    open(unit=101, file=argv(2))
    do k = 1, ommp_n_ipd
        do j = 1, ommp_pol_atoms
            write(101, "(3F20.12)") ommp_ipd(:, j, k)
        end do
    end do
    close(101)

end program test_ipd
