module anagg_print
    abstract interface
    subroutine grad_term(s, grad)
        use mod_mmpol, only: ommp_system
        use mod_memory, only: rp
        type(ommp_system), intent(inout), target :: s
        real(rp), intent(out) :: grad(3,s%top%mm_atoms)
    end subroutine
    end interface

    contains

    subroutine ana_grd_print(s, grd_f, io_file, n)
        use ommp_interface
        implicit none
        
        type(ommp_system), intent(inout), target :: s
        !! System data structure
        procedure(grad_term), pointer :: grd_f 
        !! The analytical gradients function
        integer :: io_file
        !! File handler for I/O
        character(len=*) :: n
        !! Name of the component for I/O

        real(ommp_real), allocatable :: grad_ana(:,:)
        integer(ommp_integer) :: nat, i

        nat = s%top%mm_atoms

        allocate(grad_ana(3,nat))

        call grd_f(s, grad_ana)
        write(io_file, "('Grad ', A, ' ')") n
        do i=1, nat
            write(io_file, "(F20.8, ' ', F20.8, ' ', F20.8)") grad_ana(:,i) * OMMP_AU2KCALMOL * OMMP_ANG2AU
        end do
        write(io_file, *)
    end subroutine
end module

program test_geomgrad_xyz
    use iso_c_binding, only: c_char
    use iso_fortran_env, only: error_unit
    use ommp_interface
    use anagg_print

    implicit none

    integer :: argc, i, rc=0, fp=102
    character(kind=c_char, len=120), dimension(3) :: argv
    real(ommp_real), allocatable :: ef(:,:), ef_pol(:,:), &
                                    grad_num(:,:), grad_ana(:,:), delta(:,:)
    type(ommp_system), pointer :: my_system
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
    grd_f => ommp_full_geomgrad
    call ana_grd_print(my_system, grd_f, fp, "ETOT")
    
    grd_f => ommp_vdw_geomgrad
    call ana_grd_print(my_system, grd_f, fp, "EV")
    
    grd_f => ommp_fixedelec_geomgrad
    call ana_grd_print(my_system, grd_f, fp, "EM")
    
    grd_f => ommp_polelec_geomgrad
    call ana_grd_print(my_system, grd_f, fp, "EP")
    
    grd_f => ommp_bond_geomgrad
    call ana_grd_print(my_system, grd_f, fp, "EB")
    
    grd_f => ommp_torsion_geomgrad
    call ana_grd_print(my_system, grd_f, fp, "ET")
    
    grd_f => ommp_angle_geomgrad
    call ana_grd_print(my_system, grd_f, fp, "EA")
    
    grd_f => ommp_urey_geomgrad
    call ana_grd_print(my_system, grd_f, fp, "EUB")
    
    grd_f => ommp_imptorsion_geomgrad
    call ana_grd_print(my_system, grd_f, fp, "EIT")
    
    grd_f => ommp_strbnd_geomgrad
    call ana_grd_print(my_system, grd_f, fp, "EBA")
    
    grd_f => ommp_angtor_geomgrad
    call ana_grd_print(my_system, grd_f, fp, "EAT")
    
    grd_f => ommp_opb_geomgrad
    call ana_grd_print(my_system, grd_f, fp, "EOPB")
    
    grd_f => ommp_strtor_geomgrad
    call ana_grd_print(my_system, grd_f, fp, "EBT")
    
    grd_f => ommp_tortor_geomgrad
    call ana_grd_print(my_system, grd_f, fp, "ETT")
    
    grd_f => ommp_pitors_geomgrad
    call ana_grd_print(my_system, grd_f, fp, "EPT")
    close(fp)
    if(rc > 0) then
        write(error_unit, *) 1
        stop 1
    else
        write(error_unit, *) 0
        stop 0
    end if
end program
