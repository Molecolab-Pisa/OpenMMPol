module test_geomgrad
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
    function energy_term_qmmm(s, qmh, fakeqm)
        use mod_mmpol, only: ommp_system
        use mod_qm_helper, only: ommp_qm_helper
        use mod_memory, only: rp
        type(ommp_system), intent(inout), target :: s, fakeqm
        type(ommp_qm_helper), intent(inout), target :: qmh
        real(rp) :: energy_term_qmmm
    end function
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
       
        subroutine update_coordinates_qmmm(s, qmh, fakeqm, new_mm_c, new_qm_c)
            
            implicit none
            
            type(ommp_system), intent(inout) :: s, fakeqm
            !! System data structure
            type(ommp_qm_helper), intent(inout) :: qmh
            real(ommp_real), intent(inout) :: new_mm_c(:,:), new_qm_c(:,:)

            call ommp_update_coordinates(s, new_mm_c)
            call ommp_qm_helper_update_coord(qmh, new_qm_c)
            if(s%use_linkatoms) then
                call ommp_update_link_atoms_position(qmh, s)
                call ommp_update_coordinates(fakeqm, qmh%qm_top%cmm)
            end if
        end subroutine

        subroutine numerical_geomgrad_qmmm(s, qmh, fakeqm, ene_f, grad, qmgrad)
            implicit none
            
            type(ommp_system), intent(inout) :: s, fakeqm
            !! System data structure
            type(ommp_qm_helper), intent(inout) :: qmh
            procedure(energy_term_qmmm), pointer :: ene_f
            !! The energy function (from interface module) for which
            !! numerical gradients are needed
            real(ommp_real), dimension(3,s%top%mm_atoms), intent(inout) :: grad
            !! Geometrical gradients in output, results will be added
            real(ommp_real), dimension(3,qmh%qm_top%mm_atoms), intent(inout) :: qmgrad

            integer(ommp_integer) :: i, j
            real(ommp_real), allocatable :: new_c_mm(:,:), new_c_qm(:,:)
            real(ommp_real) :: tmp
            real(ommp_real), parameter :: dd = 1e-5

            allocate(new_c_mm(3, s%top%mm_atoms))
            new_c_mm = s%top%cmm
            
            allocate(new_c_qm(3, qmh%qm_top%mm_atoms))
            new_c_qm = qmh%qm_top%cmm
            
            do i=1, s%top%mm_atoms
                do j=1, 3
                    new_c_mm(j,i) = new_c_mm(j,i) + dd
                    call update_coordinates_qmmm(s, qmh, fakeqm, new_c_mm, new_c_qm)
                    tmp = ene_f(s, qmh, fakeqm)

                    new_c_mm(j,i) = new_c_mm(j,i) - 2*dd
                    call update_coordinates_qmmm(s, qmh, fakeqm, new_c_mm, new_c_qm)
                    tmp = tmp - ene_f(s, qmh, fakeqm)
                    
                    grad(j,i) = grad(j,i) + tmp / (2*dd)
                    
                    new_c_mm(j,i) = new_c_mm(j,i) + dd
                    call update_coordinates_qmmm(s, qmh, fakeqm, new_c_mm, new_c_qm)
                end do
            end do
            
            do i=1, qmh%qm_top%mm_atoms
                do j=1, 3
                    new_c_qm(j,i) = new_c_qm(j,i) + dd
                    call update_coordinates_qmmm(s, qmh, fakeqm, new_c_mm, new_c_qm)
                    tmp = ene_f(s, qmh, fakeqm)

                    new_c_qm(j,i) = new_c_qm(j,i) - 2*dd
                    call update_coordinates_qmmm(s, qmh, fakeqm, new_c_mm, new_c_qm)
                    tmp = tmp - ene_f(s, qmh, fakeqm)
                    
                    qmgrad(j,i) = qmgrad(j,i) + tmp / (2*dd)
                    
                    new_c_qm(j,i) = new_c_qm(j,i) + dd
                    call update_coordinates_qmmm(s, qmh, fakeqm, new_c_mm, new_c_qm)
                end do
            end do
            
            deallocate(new_c_qm, new_c_mm) 
        end subroutine

    subroutine print_qmmm_grad(n, mmat, qmat, mmg, qmg)
        character(len=*) :: n
        integer(ommp_integer) :: mmat, qmat
        real(ommp_real) :: mmg(3,mmat), qmg(3,qmat)

        character(len=OMMP_STR_CHAR_MAX) :: msg
        integer(ommp_integer) :: i
        write(msg, "('Grad ', A)") n
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-GRD")

        do i=1, mmat
            write(msg, "('MM:', I0, E20.12, ' ', E20.12, ' ', E20.12)") &
                i, mmg(:,i)
            call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-GRD")
        end do

        if(qmat > 0) then
            do i=1, qmat
                write(msg, "('QM:', I08, E20.12, ' ', E20.12, ' ', E20.12)") &
                    i, qmg(:,i)
                call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-GRD")
            end do
        end if

        call ommp_message("", OMMP_VERBOSE_NONE, "TEST-GRD")
    end subroutine

    subroutine num_grd_print(sys, ene_f, n)
        type(ommp_system) :: sys
        procedure(energy_term), pointer :: ene_f
        character(len=*) :: n
        
        real(ommp_real), allocatable, dimension(:,:) :: gmm
        real(ommp_real), dimension(3,0) :: gqm

        allocate(gmm(3, sys%top%mm_atoms))
        gmm = 0.0
        call numerical_geomgrad(sys, ene_f, gmm)
        gmm = gmm * OMMP_AU2KCALMOL*OMMP_ANG2AU
        call print_qmmm_grad(n, sys%top%mm_atoms, 0, gmm, gqm)
        deallocate(gmm)

    end subroutine

    !subroutine num_grd_print_qmmm(sys, qmh, fakeqm, grad_f, n)
    !    type(ommp_system) :: sys, fakeqm
    !    type(ommp_qm_helper) :: qmh
    !    procedure(grad_term_qmmm), pointer :: grad_f
    !    character(len=*) :: n
    !    
    !    real(ommp_real), allocatable, dimension(:,:) :: gmm
    !    real(ommp_real), allocatable, dimension(:,:) :: gqm

    !    allocate(gmm(3, sys%top%mm_atoms))
    !    gmm = 0.0
    !    allocate(gqm(3, qmh%qm_top%mm_atoms))
    !    gqm = 0.0
    !    call grad_f(sys, qmh, fakeqm, gmm, gqm)
    !    gmm = gmm * OMMP_AU2KCALMOL*OMMP_ANG2AU
    !    gqm = gqm * OMMP_AU2KCALMOL*OMMP_ANG2AU
    !    call print_qmmm_grad(n, sys%top%mm_atoms, qmh%qm_top%mm_atoms, gmm, gqm)
    !    deallocate(gmm)
    !    deallocate(gqm)

    !end subroutine
end module

program test_SI_geomgrad_num
    use iso_c_binding, only: c_char
    use ommp_interface
    use test_geomgrad

    implicit none
    character(kind=c_char, len=120), dimension(3) :: args
    character(len=OMMP_STR_CHAR_MAX) :: msg, prm_file
    integer :: narg
    type(ommp_system), pointer :: my_system, fake_qm
    type(ommp_qm_helper), pointer :: my_qmh
    logical :: use_qm = .false., use_fake_qm = .false.

    procedure(energy_term), pointer :: gt
    !procedure(energy_term_qmmm), pointer :: gt_qmmm
  
    narg = command_argument_count()
    if (narg /= 2 .and. narg /= 1) then
        write(6, *) "Syntax expected "
        write(6, *) "   $ test_SI_geomgrad_num.exe <JSON FILE> <OUTPUT FILE>"
        call exit(1)
    else 
        call get_command_argument(1, args(1))
        call get_command_argument(2, args(2))
       
        call ommp_smartinput(trim(args(1)), my_system, my_qmh)
        call ommp_set_outputfile(trim(args(2)))

        if(associated(my_qmh)) then
            use_qm = .true.
            if(my_system%use_linkatoms) then
                use_fake_qm = .true.
                call ommp_smartinput_cpstr(trim(args(1)), "qm/prm_file/path", prm_file);
                call ommp_system_from_qm_helper(my_qmh, trim(prm_file), fake_qm);
                call ommp_turn_pol_off(fake_qm, fake_qm%top%mm_atoms, fake_qm%eel%polar_mm)
            else
                allocate(fake_qm)
            end if
        end if
        
        gt => ommp_get_fixedelec_energy
        call num_grd_print(my_system, gt, "EM")
        gt => ommp_get_polelec_energy
        call num_grd_print(my_system, gt, "EP");
        gt => ommp_get_vdw_energy
        call num_grd_print(my_system, gt, "EV");
        gt => ommp_get_bond_energy
        call num_grd_print(my_system, gt, "EB");
        gt => ommp_get_angle_energy
        call num_grd_print(my_system, gt, "EA");
        gt => ommp_get_strbnd_energy
        call num_grd_print(my_system, gt, "EBA");
        gt => ommp_get_urey_energy
        call num_grd_print(my_system, gt, "EUB");
        gt => ommp_get_opb_energy
        call num_grd_print(my_system, gt, "EOPB");
        gt => ommp_get_pitors_energy
        call num_grd_print(my_system, gt, "EPT");
        gt => ommp_get_torsion_energy
        call num_grd_print(my_system, gt, "ET");
        gt => ommp_get_tortor_energy
        call num_grd_print(my_system, gt, "ETT");
        gt => ommp_get_imptorsion_energy
        call num_grd_print(my_system, gt, "EIT");
        gt => ommp_get_strtor_energy
        call num_grd_print(my_system, gt, "EBT");
        gt => ommp_get_angtor_energy
        call num_grd_print(my_system, gt, "EAT");
        !if(use_qm) then
        !    gt_qmmm => ommptest_qm_helper_vdw_geomgrad
        !    call num_grd_print_qmmm(my_system, my_qmh, fake_qm, gt_qmmm, "EVQMMM")
        !    if(use_fake_qm) then
        !        gt_qmmm => ommptest_fakeqm_internal_geomgrad
        !        call num_grd_print_qmmm(my_system, my_qmh, fake_qm, gt_qmmm, "EQ")
        !        gt_qmmm => ommptest_fakeqm_linkatom_geomgrad
        !        call num_grd_print_qmmm(my_system, my_qmh, fake_qm, gt_qmmm, "ELA")
        !    end if
        !    gt_qmmm => ommptest_totalqmmm_geomgrad
        !    call num_grd_print_qmmm(my_system, my_qmh, fake_qm, gt_qmmm, "ETOT")
        !else
        !    gt => ommp_get_full_geomgrad
        !    call num_grd_print(my_system, gt, "ETOT")
        !end if

        if(associated(my_qmh)) call ommp_terminate_qm_helper(my_qmh)
        if(associated(my_system)) call ommp_terminate(my_system)
    end if
end program 
