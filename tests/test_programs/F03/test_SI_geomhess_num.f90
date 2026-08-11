module test_geomhess_num
    use ommp_interface

    abstract interface
    subroutine grad_term(s, grad)
        use mod_mmpol, only: ommp_system
        use mod_memory, only: rp
        type(ommp_system), intent(inout), target :: s
        real(rp), intent(out) :: grad(3,s%top%mm_atoms)
    end subroutine
    end interface

    contains
        subroutine full_elec_geomgrad(s, grad)
            !! Fixed-multipole + polarization gradient, combined -- the
            !! numerical counterpart of full_elec_geomhess in
            !! test_geomhess (see that module for why this only covers
            !! electrostatics, not the true ommp_full_geomgrad).
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(out) :: grad(3,s%top%mm_atoms)
            real(ommp_real), allocatable :: tmp(:,:)

            allocate(tmp(3,s%top%mm_atoms))
            call ommp_fixedelec_geomgrad(s, grad)
            call ommp_polelec_geomgrad(s, tmp)
            grad = grad + tmp
            deallocate(tmp)
        end subroutine

        subroutine numerical_geomhess(s, grad_f, hess)
            !! Numerical Hessian obtained by central-differencing an
            !! analytical gradient routine (grad_f) wrt every Cartesian
            !! coordinate of every atom -- much better conditioned than
            !! double-differencing the energy, and it's how every
            !! analytical piece of HessPolExpl/Hess1 was validated during
            !! development.
            use mod_mmpol, only: update_coordinates
            implicit none

            type(ommp_system), intent(inout) :: s
            !! System data structure
            procedure(grad_term), pointer :: grad_f
            !! The gradient function (from interface module) for which
            !! a numerical Hessian is needed
            real(ommp_real), dimension(3,3,s%top%mm_atoms,s%top%mm_atoms), intent(inout) :: hess
            !! Geometrical Hessian in output, results will be added

            integer(ommp_integer) :: i, j, k, l
            real(ommp_real), allocatable :: new_c(:,:), g_p(:,:), g_m(:,:)
            real(ommp_real), parameter :: dd = 1.0e-4

            allocate(new_c(3, s%top%mm_atoms))
            new_c = s%top%cmm
            allocate(g_p(3, s%top%mm_atoms), g_m(3, s%top%mm_atoms))

            do i=1, s%top%mm_atoms
                do j=1, 3
                    new_c(j,i) = new_c(j,i) + dd
                    call update_coordinates(s, new_c)
                    g_p = 0.0
                    call grad_f(s, g_p)

                    new_c(j,i) = new_c(j,i) - 2*dd
                    call update_coordinates(s, new_c)
                    g_m = 0.0
                    call grad_f(s, g_m)

                    new_c(j,i) = new_c(j,i) + dd
                    call update_coordinates(s, new_c)

                    ! hess(k-direction; atom l, j-direction; atom i) --
                    ! same (kappa,lambda,atom3,atom4) convention as
                    ! fixedelec_geomhess/polelec_geomhess.
                    do k=1, s%top%mm_atoms
                        do l=1, 3
                            hess(l,j,k,i) = hess(l,j,k,i) + (g_p(l,k) - g_m(l,k)) / (2*dd)
                        end do
                    end do
                end do
            end do

            deallocate(new_c, g_p, g_m)
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

        subroutine num_hess_print(sys, grad_f, n)
            type(ommp_system) :: sys
            procedure(grad_term), pointer :: grad_f
            character(len=*) :: n

            real(ommp_real), allocatable, dimension(:,:,:,:) :: hmm
            integer(ommp_integer) :: natm

            natm = sys%top%mm_atoms
            allocate(hmm(3,3,natm,natm))
            hmm = 0.0
            call numerical_geomhess(sys, grad_f, hmm)
            ! one extra power of ANG2AU vs. the gradient's own scaling,
            ! since the Hessian is a second position-derivative
            hmm = hmm * OMMP_AU2KCALMOL*OMMP_ANG2AU*OMMP_ANG2AU
            call print_hess(natm, hmm, n)
            deallocate(hmm)
        end subroutine
end module

program test_SI_geomhess_num
    use iso_c_binding, only: c_char
    use ommp_interface
    use test_geomhess_num

    implicit none
    character(kind=c_char, len=120), dimension(2) :: args
    integer :: narg
    type(ommp_system), pointer :: my_system
    type(ommp_qm_helper), pointer :: my_qmh

    procedure(grad_term), pointer :: gt

    narg = command_argument_count()
    if (narg /= 2) then
        write(6, *) "Syntax expected "
        write(6, *) "   $ test_SI_geomhess_num.exe <JSON FILE> <OUTPUT FILE>"
        call exit(1)
    else
        call get_command_argument(1, args(1))
        call get_command_argument(2, args(2))

        call ommp_smartinput(trim(args(1)), my_system, my_qmh)
        call ommp_set_outputfile(trim(args(2)))

        gt => ommp_fixedelec_geomgrad
        call num_hess_print(my_system, gt, "EM")
        gt => ommp_polelec_geomgrad
        call num_hess_print(my_system, gt, "EP")
        gt => full_elec_geomgrad
        call num_hess_print(my_system, gt, "ETOT")

        if(associated(my_qmh)) call ommp_terminate_qm_helper(my_qmh)
        if(associated(my_system)) call ommp_terminate(my_system)
    end if
end program
