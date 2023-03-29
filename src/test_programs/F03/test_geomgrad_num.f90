module mod_test_numerical_geomgrad
    use ommp_interface
    
    abstract interface
    function energy_term(s)
            use mod_mmpol, only: ommp_system
            use mod_memory, only: rp
            type(ommp_system), intent(inout), target :: s
            real(rp) :: energy_term
    end function
    end interface

    abstract interface
    subroutine grad_term(s, grad)
        use mod_mmpol, only: ommp_system
        use mod_memory, only: rp
        type(ommp_system), intent(inout), target :: s
        real(rp) :: grad(:,:)
    end subroutine
    end interface

    contains

        subroutine numerical_geomgrad(s, ene_f, grad)
            use mod_mmpol, only: update_coordinates
            implicit none
            
            type(ommp_system), intent(inout) :: s
            !! System data structure
            procedure(energy_term), pointer :: ene_f
            !! The energy function (from interface module) for which
            !! numerical gradients are needed
            real(ommp_real), dimension(3,s%top%mm_atoms), intent(inout) :: grad
            !! Geometrical gradients in output, results will be added

            integer(ommp_integer) :: i, j
            real(ommp_real), allocatable :: new_c(:,:)
            real(ommp_real) :: tmp
            real(ommp_real), parameter :: dd = 1e-5

            allocate(new_c(3, s%top%mm_atoms))
            new_c = s%top%cmm
            
            do i=1, s%top%mm_atoms
                do j=1, 3
                    new_c(j,i) = new_c(j,i) + dd
                    call update_coordinates(s, new_c)
                    tmp = ene_f(s)

                    new_c(j,i) = new_c(j,i) - 2*dd
                    call update_coordinates(s, new_c)
                    tmp = tmp - ene_f(s)
                    
                    grad(j,i) = grad(j,i) + tmp / (2*dd)
                    
                    new_c(j,i) = new_c(j,i) + dd
                    call update_coordinates(s, new_c)
                end do
            end do
            
            deallocate(new_c) 
        end subroutine
        
        function num_ana_compare(s, ene_f, grd_f, io_file, n, del) result(ok)
            use mod_mmpol, only: update_coordinates
            implicit none
            
            type(ommp_system), intent(inout) :: s
            !! System data structure
            procedure(energy_term), pointer :: ene_f
            !! The energy function (from interface module) for which
            !! numerical gradients are needed
            procedure(grad_term), pointer :: grd_f 
            !! The analytical gradients function
            integer :: io_file
            !! File handler for I/O
            character(len=*) :: n
            !! Name of the component for I/O
            real(ommp_real) :: del
            !! Maximum absolute difference between numerical and 
            !! analytical

            real(ommp_real), allocatable :: grad_num(:,:), grad_ana(:,:)
            real(ommp_real) :: MDelta, delta(3)
            integer(ommp_integer) :: nat, i, ok

            nat = s%top%mm_atoms

            allocate(grad_num(3,nat))
            allocate(grad_ana(3,nat))

            call grd_f(s, grad_ana)
            call numerical_geomgrad(s, ene_f, grad_num)
            MDelta = 0.0
            do i=1, nat
                write(io_file, "('[',I5,'] (A) ', 3F12.8)") i, grad_ana(:,i)
                write(io_file, "('        (N) ', 3F12.8)") i, grad_num(:,i)
                delta = grad_num(:,i) - grad_ana(:,i)
                if(maxval(abs(delta)) > MDelta) MDelta = maxval(abs(delta))
                write(io_file, "('        (D) ', 3F12.8)") i, delta
                write(io_file, *)
            end do

            if(Mdelta > del) then
                write(io_file, "('WARNING delta (', F5.3, ') > max_delta (', F5.3, ')')") MDelta, del
                ok = 1
            else
                ok = 0
            end if
        end function
end module

program test_geomgrad_num
    use iso_c_binding, only: c_char
    use iso_fortran_env, only: error_unit
    use ommp_interface
    use mod_test_numerical_geomgrad

    implicit none

    integer :: argc, i, retcode=0
    character(kind=c_char, len=120), dimension(3) :: argv
    real(ommp_real), allocatable :: ef(:,:), ef_pol(:,:), &
                                    grad_num(:,:), grad_ana(:,:), delta(:,:)
    type(ommp_system), pointer :: my_system
    procedure(energy_term), pointer :: ene_f

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

    call ommp_set_verbose(OMMP_VERBOSE_LOW)
    call ommp_init_mmp(my_system, argv(1))

    allocate(ef(3, my_system%top%mm_atoms))
    allocate(ef_pol(3, my_system%eel%pol_atoms))
    
    if(argc == 3) then
        write(*, *) "Currently not supported!"
        write(error_unit, *) retcode
        stop
    else
        ef = 0.0
    end if

    do i=1, my_system%eel%pol_atoms
        ef_pol(:,i) = ef(:, my_system%eel%polar_mm(i))
    end do
    
    deallocate(ef)
    
    allocate(grad_num(3, my_system%top%mm_atoms))
    allocate(grad_ana(3, my_system%top%mm_atoms))
    allocate(delta(3, my_system%top%mm_atoms))

    grad_num = 0.0
    ene_f => ommp_get_fixedelec_energy
    call numerical_geomgrad(my_system, ene_f, grad_num)
    grad_ana = 0.0
    call ommp_fixedelec_geomgrad(my_system, grad_ana)
    delta = grad_num-grad_ana

    write(*, *) "DELTA NUM - ANA FIXEDELEC"
    do i=1, my_system%top%mm_atoms
        write(*, *) delta(:,i)
    end do
    if(maxval(abs(delta)) > 1e-11) then
        write(*, *) "Numerical-Analytical gradients difference is too large (fixed)."
        retcode = retcode + 1
    end if
    
    grad_num = 0.0
    ene_f => ommp_get_polelec_energy
    call numerical_geomgrad(my_system, ene_f, grad_num)
    grad_ana = 0.0
    call ommp_polelec_geomgrad(my_system, grad_ana)
    delta = grad_num-grad_ana

    write(*, *) "DELTA NUM - ANA POLELEC"
    do i=1, my_system%top%mm_atoms
        write(*, *) delta(:,i)
    end do
    if(maxval(abs(delta)) > 1e-8) then
        write(*, *) "Numerical-Analytical gradients difference is too large (pol)."
        retcode = retcode + 2
    end if
    
    deallocate(grad_num, grad_ana)
    deallocate(ef_pol)

    call ommp_terminate(my_system)
    
    if(retcode > 0) then
        stop 1
    else
        stop 0
    end if
end program
