module test_geomgrad
    use ommp_interface

    abstract interface
    subroutine grad_term(s, grad)
        use mod_mmpol, only: ommp_system
        use mod_memory, only: rp
        type(ommp_system), intent(inout), target :: s
        real(rp), intent(out) :: grad(3,s%top%mm_atoms)
    end subroutine
    end interface
    
    abstract interface
    subroutine grad_term_qmmm(sys, qmh, fakeqm, gmm, gqm)
        use mod_mmpol, only: ommp_system
        use mod_qm_helper, only: ommp_qm_helper
        use mod_memory, only: rp
        type(ommp_system), intent(inout), target :: sys, fakeqm
        type(ommp_qm_helper), intent(inout), target :: qmh
        real(rp), intent(out) :: gmm(3,sys%top%mm_atoms), gqm(3,qmh%qm_top%mm_atoms)
    end subroutine
    end interface
    
    contains

    subroutine ommptest_fakeqm_internal_geomgrad(sys, qmh, fakeqm, gmm, gqm)
        type(ommp_system), intent(inout), target :: sys, fakeqm
        type(ommp_qm_helper), intent(inout), target :: qmh
        real(ommp_real), intent(out) :: gmm(3,sys%top%mm_atoms), gqm(3,qmh%qm_top%mm_atoms)

        gmm = 0.0
        call ommp_full_geomgrad(fakeqm, gqm)
    end subroutine

    subroutine ommptest_fakeqm_linkatom_geomgrad(sys, qmh, fakeqm, gmm, gqm)
        type(ommp_system), intent(inout), target :: sys, fakeqm
        type(ommp_qm_helper), intent(inout), target :: qmh
        real(ommp_real), intent(out) :: gmm(3,sys%top%mm_atoms), gqm(3,qmh%qm_top%mm_atoms)
        
        real(ommp_real), allocatable, dimension(:,:) :: tmp
        allocate(tmp(3, qmh%qm_top%mm_atoms))
        tmp = 0.0
        call ommp_qm_helper_vdw_geomgrad(qmh, sys, gqm, gmm)
        call ommp_full_geomgrad(fakeqm, tmp)
        tmp = tmp + gqm
        call ommp_qm_helper_link_atom_geomgrad(qmh, sys, gqm, gmm, tmp);
        deallocate(tmp)
    end subroutine

    subroutine ommptest_qm_helper_vdw_geomgrad(sys, qmh, fakeqm, gmm, gqm)
        type(ommp_system), intent(inout), target :: sys, fakeqm
        type(ommp_qm_helper), intent(inout), target :: qmh
        real(ommp_real), intent(out) :: gmm(3,sys%top%mm_atoms), gqm(3,qmh%qm_top%mm_atoms)

        call ommp_qm_helper_vdw_geomgrad(qmh, sys, gqm, gmm)
    end subroutine

    subroutine ommptest_totalqmmm_geomgrad(sys, qmh, fakeqm, gmm, gqm)
        type(ommp_system), intent(inout), target :: sys, fakeqm
        type(ommp_qm_helper), intent(inout), target :: qmh
        real(ommp_real), intent(out) :: gmm(3,sys%top%mm_atoms), gqm(3,qmh%qm_top%mm_atoms)

        real(ommp_real), allocatable, dimension(:,:) :: tmpmm, tmpqm

        gmm = 0.0
        gqm = 0.0
        allocate(tmpmm(3, sys%top%mm_atoms))
        tmpmm = 0.0
        allocate(tmpqm(3, qmh%qm_top%mm_atoms))
        tmpqm = 0.0

        call ommp_qm_helper_vdw_geomgrad(qmh, sys, tmpqm, tmpmm)
        gmm = gmm + tmpmm
        gqm = gqm + tmpqm
        if(sys%use_linkatoms) then
            call ommp_full_geomgrad(fakeqm, tmpqm)
            gqm = gqm + tmpqm
        end if
        call ommp_full_geomgrad(sys, tmpmm)
        gmm = gmm + tmpmm
        call ommp_qm_helper_link_atom_geomgrad(qmh, sys, tmpqm, tmpmm, gqm)
        gmm = gmm + tmpmm
        gqm = gqm + tmpqm
        deallocate(tmpmm, tmpqm)
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

    subroutine ana_grd_print(sys, grad_f, n)
        type(ommp_system) :: sys
        procedure(grad_term), pointer :: grad_f
        character(len=*) :: n
        
        real(ommp_real), allocatable, dimension(:,:) :: gmm
        real(ommp_real), dimension(3,0) :: gqm

        allocate(gmm(3, sys%top%mm_atoms))
        gmm = 0.0
        call grad_f(sys, gmm)
        gmm = gmm * OMMP_AU2KCALMOL*OMMP_ANG2AU
        call print_qmmm_grad(n, sys%top%mm_atoms, 0, gmm, gqm)
        deallocate(gmm)

    end subroutine

    subroutine ana_grd_print_qmmm(sys, qmh, fakeqm, grad_f, n)
        type(ommp_system) :: sys, fakeqm
        type(ommp_qm_helper) :: qmh
        procedure(grad_term_qmmm), pointer :: grad_f
        character(len=*) :: n
        
        real(ommp_real), allocatable, dimension(:,:) :: gmm
        real(ommp_real), allocatable, dimension(:,:) :: gqm

        allocate(gmm(3, sys%top%mm_atoms))
        gmm = 0.0
        allocate(gqm(3, qmh%qm_top%mm_atoms))
        gqm = 0.0
        call grad_f(sys, qmh, fakeqm, gmm, gqm)
        gmm = gmm * OMMP_AU2KCALMOL*OMMP_ANG2AU
        gqm = gqm * OMMP_AU2KCALMOL*OMMP_ANG2AU
        call print_qmmm_grad(n, sys%top%mm_atoms, qmh%qm_top%mm_atoms, gmm, gqm)
        deallocate(gmm)
        deallocate(gqm)

    end subroutine
end module

program test_SI_geomgrad
    use iso_c_binding, only: c_char
    use ommp_interface
    use test_geomgrad

    implicit none
    character(kind=c_char, len=120), dimension(3) :: args
    character(len=OMMP_STR_CHAR_MAX) :: prm_file
    integer :: narg
    type(ommp_system), pointer :: my_system, fake_qm
    type(ommp_qm_helper), pointer :: my_qmh
    logical :: use_qm = .false., use_fake_qm = .false.

    procedure(grad_term), pointer :: gt
    procedure(grad_term_qmmm), pointer :: gt_qmmm
  
    narg = command_argument_count()
    if (narg /= 2 .and. narg /= 1) then
        write(6, *) "Syntax expected "
        write(6, *) "   $ test_SI_geomgrad.exe <JSON FILE> <OUTPUT FILE>"
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
        
        gt => ommp_fixedelec_geomgrad
        call ana_grd_print(my_system, gt, "EM")
        gt => ommp_polelec_geomgrad
        call ana_grd_print(my_system, gt, "EP");
        gt => ommp_vdw_geomgrad
        call ana_grd_print(my_system, gt, "EV");
        gt => ommp_bond_geomgrad
        call ana_grd_print(my_system, gt, "EB");
        gt => ommp_angle_geomgrad
        call ana_grd_print(my_system, gt, "EA");
        gt => ommp_strbnd_geomgrad
        call ana_grd_print(my_system, gt, "EBA");
        gt => ommp_urey_geomgrad
        call ana_grd_print(my_system, gt, "EUB");
        gt => ommp_opb_geomgrad
        call ana_grd_print(my_system, gt, "EOPB");
        gt => ommp_pitors_geomgrad
        call ana_grd_print(my_system, gt, "EPT");
        gt => ommp_torsion_geomgrad
        call ana_grd_print(my_system, gt, "ET");
        gt => ommp_tortor_geomgrad
        call ana_grd_print(my_system, gt, "ETT");
        gt => ommp_imptorsion_geomgrad
        call ana_grd_print(my_system, gt, "EIT");
        gt => ommp_strtor_geomgrad
        call ana_grd_print(my_system, gt, "EBT");
        gt => ommp_angtor_geomgrad
        call ana_grd_print(my_system, gt, "EAT");
        if(use_qm) then
            gt_qmmm => ommptest_qm_helper_vdw_geomgrad
            call ana_grd_print_qmmm(my_system, my_qmh, fake_qm, gt_qmmm, "EVQMMM")
            if(use_fake_qm) then
                gt_qmmm => ommptest_fakeqm_internal_geomgrad
                call ana_grd_print_qmmm(my_system, my_qmh, fake_qm, gt_qmmm, "EQ")
                gt_qmmm => ommptest_fakeqm_linkatom_geomgrad
                call ana_grd_print_qmmm(my_system, my_qmh, fake_qm, gt_qmmm, "ELA")
            end if
            gt_qmmm => ommptest_totalqmmm_geomgrad
            call ana_grd_print_qmmm(my_system, my_qmh, fake_qm, gt_qmmm, "ETOT")
        else
            gt => ommp_full_geomgrad
            call ana_grd_print(my_system, gt, "ETOT")
        end if

        if(associated(my_qmh)) call ommp_terminate_qm_helper(my_qmh)
        if(associated(my_system)) call ommp_terminate(my_system)
    end if
end program
