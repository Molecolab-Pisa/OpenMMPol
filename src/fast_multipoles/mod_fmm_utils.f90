module mod_fmm_utils
    use mod_constants, only : ip
    
    implicit none
    private

    public :: fmm_error, n_sph_harm, ntot_sph_harm

    contains
        pure function n_sph_harm(l) result(n)
            !! Return the number of components for spherical harmonics
            !! of order n
            integer(ip), intent(in) :: l

            integer(ip) :: n

            n = 2*l+1
            return
        end function

        pure function ntot_sph_harm(l) result(n)
            !! Return the number of components for spherical harmonics
            !! up to order n
            integer(ip), intent(in) :: l

            integer(ip) :: n

            n = (l+1) ** 2
            return
        end function


        subroutine fmm_error(s)
            character(len=*), optional :: s
            if(present(s)) write(*, *) s
            stop 1
        end subroutine
end module mod_fmm_utils
