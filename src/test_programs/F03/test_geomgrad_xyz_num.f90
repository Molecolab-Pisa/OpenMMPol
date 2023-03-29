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
        real(rp), intent(out) :: grad(3,s%top%mm_atoms)
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
            
            type(ommp_system), intent(inout), target :: s
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
            write(io_file, "('DELTA NUM - ANA ', A)") n
            do i=1, nat
                write(io_file, "('[',I5,'] (A) ', 3F12.8)") i, grad_ana(:,i)
                write(io_file, "('        (N) ', 3F12.8)") grad_num(:,i)
                delta = grad_num(:,i) - grad_ana(:,i)
                if(maxval(abs(delta)) > MDelta) MDelta = maxval(abs(delta))
                write(io_file, "('        (D) ', 3F12.8)") delta
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

program test_geomgrad_xyz_num
    use iso_c_binding, only: c_char
    use iso_fortran_env, only: error_unit
    use ommp_interface
    use mod_test_numerical_geomgrad

    implicit none

    integer :: argc, i, rc=0, fp=102
    character(kind=c_char, len=120), dimension(3) :: argv
    real(ommp_real), allocatable :: ef(:,:), ef_pol(:,:), &
                                    grad_num(:,:), grad_ana(:,:), delta(:,:)
    type(ommp_system), pointer :: my_system
    procedure(energy_term), pointer :: ene_f
    procedure(grad_term), pointer :: grd_f

    argc = command_argument_count()
    if(argc /= 4 .and. argc /= 3 ) then
        write(*, *) "Syntax expected"
        write(*, *) "    $ test_init.exe <XYZ FILE> <PRM FILE> <OUTPUT FILE>&
                    & [<ELECTRIC FIELD FILE>]"
        stop 1
    end if

    do i=1, argc
        call get_command_argument(i, argv(i))
    end do

    call ommp_set_verbose(OMMP_VERBOSE_LOW)
    call ommp_init_xyz(my_system, argv(1), argv(2))

    allocate(ef(3, my_system%top%mm_atoms))
    allocate(ef_pol(3, my_system%eel%pol_atoms))
    
    if(argc == 4) then
        write(*, *) "Currently not supported!"
        write(error_unit, *) 1
        stop 1
    else
        ef = 0.0
    end if

    do i=1, my_system%eel%pol_atoms
        ef_pol(:,i) = ef(:, my_system%eel%polar_mm(i))
    end do
    
    deallocate(ef)
    
    open(fp, file=argv(3))
    ene_f => ommp_get_full_energy
    grd_f => ommp_full_geomgrad
    rc = num_ana_compare(my_system, ene_f, grd_f, fp, "CompletePotential", 1e-7_ommp_real)
    
    if(rc == 0) stop 0
    
    rc = 0
    ene_f => ommp_get_vdw_energy
    grd_f => ommp_vdw_geomgrad
    rc = rc + num_ana_compare(my_system, ene_f, grd_f, fp, "non-bonded", 1e-7_ommp_real)
    
    ene_f => ommp_get_fixedelec_energy
    grd_f => ommp_fixedelec_geomgrad
    rc = rc + num_ana_compare(my_system, ene_f, grd_f, fp, "fixedelec", 1e-7_ommp_real)
    
    ene_f => ommp_get_polelec_energy
    grd_f => ommp_polelec_geomgrad
    rc = rc + num_ana_compare(my_system, ene_f, grd_f, fp, "polelec", 1e-7_ommp_real)
    
    ene_f => ommp_get_bond_energy
    grd_f => ommp_bond_geomgrad
    rc = rc + num_ana_compare(my_system, ene_f, grd_f, fp, "bond", 1e-7_ommp_real)
    
    ene_f => ommp_get_angle_energy
    grd_f => ommp_angle_geomgrad
    rc = rc + num_ana_compare(my_system, ene_f, grd_f, fp, "angle", 1e-7_ommp_real)
    
    ene_f => ommp_get_urey_energy
    grd_f => ommp_urey_geomgrad
    rc = rc + num_ana_compare(my_system, ene_f, grd_f, fp, "urey", 1e-7_ommp_real)
    
    ene_f => ommp_get_imptorsion_energy
    grd_f => ommp_imptorsion_geomgrad
    rc = rc + num_ana_compare(my_system, ene_f, grd_f, fp, "imptorsion", 1e-7_ommp_real)
    
    ene_f => ommp_get_strbnd_energy
    grd_f => ommp_strbnd_geomgrad
    rc = rc + num_ana_compare(my_system, ene_f, grd_f, fp, "strbnd", 1e-7_ommp_real)
    
    ene_f => ommp_get_angtor_energy
    grd_f => ommp_angtor_geomgrad
    rc = rc + num_ana_compare(my_system, ene_f, grd_f, fp, "angtor", 1e-7_ommp_real)
    
    ene_f => ommp_get_opb_energy
    grd_f => ommp_opb_geomgrad
    rc = rc + num_ana_compare(my_system, ene_f, grd_f, fp, "opb", 1e-7_ommp_real)
    
    ene_f => ommp_get_strtor_energy
    grd_f => ommp_strtor_geomgrad
    rc = rc + num_ana_compare(my_system, ene_f, grd_f, fp, "strtor", 1e-7_ommp_real)
    
    ene_f => ommp_get_tortor_energy
    grd_f => ommp_tortor_geomgrad
    rc = rc + num_ana_compare(my_system, ene_f, grd_f, fp, "tortor", 1e-7_ommp_real)
    
    ene_f => ommp_get_pitors_energy
    grd_f => ommp_pitors_geomgrad
    rc = rc + num_ana_compare(my_system, ene_f, grd_f, fp, "pitors", 1e-7_ommp_real)
    close(fp)
    if(rc > 0) then
        write(error_unit, *) 1
        stop 1
    else
        write(error_unit, *) 0
        stop 0
    end if
end program
