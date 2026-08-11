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
            !! (which also includes bonded and vdW terms): analytical
            !! bonded/vdW Hessians don't exist yet, so this only covers
            !! electrostatics -- the same scope as full_elec_geomgrad in
            !! test_geomhess_num, so the two remain directly comparable.
            use mod_geomhess, only: fixedelec_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            call fixedelec_geomhess(s, hess)
            call polelec_geomhess_guarded(s, hess)
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

        if(associated(my_qmh)) call ommp_terminate_qm_helper(my_qmh)
        if(associated(my_system)) call ommp_terminate(my_system)
    end if
end program
