module test_geomhess
    use ommp_interface

    abstract interface
    subroutine hess_term(s, hess)
        use mod_mmpol, only: ommp_system
        use mod_memory, only: rp
        type(ommp_system), intent(inout), target :: s
        real(rp), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)
    end subroutine
    end interface

    contains
        subroutine polelec_geomhess_guarded(s, hess)
            !! polelec_geomhess itself, skipped when there is nothing
            !! polarizable to differentiate -- mirrors the guard that
            !! ommp_polelec_geomgrad applies at the gradient level.
            use mod_geomhess, only: polelec_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            if(s%eel%pol_atoms > 0) call polelec_geomhess(s, hess)
        end subroutine

        subroutine full_elec_geomhess(s, hess)
            !! Fixed-multipole + polarization Hessian, combined. Not
            !! "the full Hessian" in the same sense as ommp_full_geomgrad
            !! (which also includes bonded and vdW terms): kept scoped to
            !! electrostatics only, matching full_elec_geomgrad in
            !! test_geomhess_num, so the two remain directly comparable
            !! (bonded and vdW terms are tested separately, as BNDTOT/VDW).
            use mod_geomhess, only: fixedelec_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            call fixedelec_geomhess(s, hess)
            call polelec_geomhess_guarded(s, hess)
        end subroutine

        subroutine bond_geomhess_guarded(s, hess)
            use mod_bonded, only: bond_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            if(s%use_bonded) call bond_geomhess(s%bds, hess)
        end subroutine

        subroutine urey_geomhess_guarded(s, hess)
            use mod_bonded, only: urey_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            if(s%use_bonded) call urey_geomhess(s%bds, hess)
        end subroutine

        subroutine angle_geomhess_guarded(s, hess)
            use mod_bonded, only: angle_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            if(s%use_bonded) call angle_geomhess(s%bds, hess)
        end subroutine

        subroutine strbnd_geomhess_guarded(s, hess)
            use mod_bonded, only: strbnd_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            if(s%use_bonded) call strbnd_geomhess(s%bds, hess)
        end subroutine

        subroutine opb_geomhess_guarded(s, hess)
            use mod_bonded, only: opb_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            if(s%use_bonded) call opb_geomhess(s%bds, hess)
        end subroutine

        subroutine pitors_geomhess_guarded(s, hess)
            use mod_bonded, only: pitors_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            if(s%use_bonded) call pitors_geomhess(s%bds, hess)
        end subroutine

        subroutine torsion_geomhess_guarded(s, hess)
            use mod_bonded, only: torsion_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            if(s%use_bonded) call torsion_geomhess(s%bds, hess)
        end subroutine

        subroutine imptorsion_geomhess_guarded(s, hess)
            use mod_bonded, only: imptorsion_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            if(s%use_bonded) call imptorsion_geomhess(s%bds, hess)
        end subroutine

        subroutine angtor_geomhess_guarded(s, hess)
            use mod_bonded, only: angtor_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            if(s%use_bonded) call angtor_geomhess(s%bds, hess)
        end subroutine

        subroutine strtor_geomhess_guarded(s, hess)
            use mod_bonded, only: strtor_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            if(s%use_bonded) call strtor_geomhess(s%bds, hess)
        end subroutine

        subroutine tortor_geomhess_guarded(s, hess)
            use mod_bonded, only: tortor_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            if(s%use_bonded) call tortor_geomhess(s%bds, hess)
        end subroutine

        subroutine full_bnd_geomhess(s, hess)
            !! Sum of all 11 analytical bonded-term Hessians, combined --
            !! the analytical counterpart of ommp_full_bnd_geomgrad, kept
            !! in the same scope (bonded only, no link-atom corrections,
            !! since no analytical link-atom Hessian exists) so it is
            !! directly comparable to full_bnd_geomgrad in test_geomhess_num.
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            call bond_geomhess_guarded(s, hess)
            call urey_geomhess_guarded(s, hess)
            call angle_geomhess_guarded(s, hess)
            call strbnd_geomhess_guarded(s, hess)
            call opb_geomhess_guarded(s, hess)
            call pitors_geomhess_guarded(s, hess)
            call torsion_geomhess_guarded(s, hess)
            call imptorsion_geomhess_guarded(s, hess)
            call angtor_geomhess_guarded(s, hess)
            call strtor_geomhess_guarded(s, hess)
            call tortor_geomhess_guarded(s, hess)
        end subroutine

        subroutine vdw_geomhess_guarded(s, hess)
            use mod_nonbonded, only: vdw_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            if(s%use_nonbonded) call vdw_geomhess(s%vdw, hess)
        end subroutine

        subroutine print_hess(n, h, name)
            character(len=*) :: name
            integer(ommp_integer) :: n
            real(ommp_real) :: h(3,3,n,n)

            character(len=OMMP_STR_CHAR_MAX) :: msg
            integer(ommp_integer) :: i, j

            write(msg, "('Hess ', A)") name
            call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-HES")

            do i=1, n
                do j=1, n
                    write(msg, "('IJ:', I0, ':', I0, 9(' ', E20.12))") &
                        i, j, h(:,1,i,j), h(:,2,i,j), h(:,3,i,j)
                    call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-HES")
                end do
            end do

            call ommp_message("", OMMP_VERBOSE_NONE, "TEST-HES")
        end subroutine

        subroutine ana_hess_print(sys, hess_f, n)
            type(ommp_system) :: sys
            procedure(hess_term), pointer :: hess_f
            character(len=*) :: n

            real(ommp_real), allocatable, dimension(:,:,:,:) :: hmm
            integer(ommp_integer) :: natm

            natm = sys%top%mm_atoms
            allocate(hmm(3,3,natm,natm))
            hmm = 0.0
            call hess_f(sys, hmm)
            ! one extra power of ANG2AU vs. the gradient's own scaling,
            ! since the Hessian is a second position-derivative
            hmm = hmm * OMMP_AU2KCALMOL*OMMP_ANG2AU*OMMP_ANG2AU
            call print_hess(natm, hmm, n)
            deallocate(hmm)
        end subroutine
end module

program test_SI_geomhess
    use iso_c_binding, only: c_char
    use ommp_interface
    use test_geomhess
    use mod_geomhess, only: fixedelec_geomhess

    implicit none
    character(kind=c_char, len=120), dimension(2) :: args
    integer :: narg
    type(ommp_system), pointer :: my_system
    type(ommp_qm_helper), pointer :: my_qmh

    procedure(hess_term), pointer :: ht

    narg = command_argument_count()
    if (narg /= 2) then
        write(6, *) "Syntax expected "
        write(6, *) "   $ test_SI_geomhess.exe <JSON FILE> <OUTPUT FILE>"
        call exit(1)
    else
        call get_command_argument(1, args(1))
        call get_command_argument(2, args(2))

        call ommp_smartinput(trim(args(1)), my_system, my_qmh)
        call ommp_set_outputfile(trim(args(2)))

        ht => fixedelec_geomhess
        call ana_hess_print(my_system, ht, "EM")
        ht => polelec_geomhess_guarded
        call ana_hess_print(my_system, ht, "EP")
        ht => full_elec_geomhess
        call ana_hess_print(my_system, ht, "ETOT")

        ht => bond_geomhess_guarded
        call ana_hess_print(my_system, ht, "BOND")
        ht => urey_geomhess_guarded
        call ana_hess_print(my_system, ht, "UREY")
        ht => angle_geomhess_guarded
        call ana_hess_print(my_system, ht, "ANGLE")
        ht => strbnd_geomhess_guarded
        call ana_hess_print(my_system, ht, "STRBND")
        ht => opb_geomhess_guarded
        call ana_hess_print(my_system, ht, "OPB")
        ht => pitors_geomhess_guarded
        call ana_hess_print(my_system, ht, "PITORS")
        ht => torsion_geomhess_guarded
        call ana_hess_print(my_system, ht, "TORSION")
        ht => imptorsion_geomhess_guarded
        call ana_hess_print(my_system, ht, "IMPTORSION")
        ht => angtor_geomhess_guarded
        call ana_hess_print(my_system, ht, "ANGTOR")
        ht => strtor_geomhess_guarded
        call ana_hess_print(my_system, ht, "STRTOR")
        ht => tortor_geomhess_guarded
        call ana_hess_print(my_system, ht, "TORTOR")
        ht => full_bnd_geomhess
        call ana_hess_print(my_system, ht, "BNDTOT")

        ht => vdw_geomhess_guarded
        call ana_hess_print(my_system, ht, "VDW")

        if(associated(my_qmh)) call ommp_terminate_qm_helper(my_qmh)
        if(associated(my_system)) call ommp_terminate(my_system)
    end if
end program
