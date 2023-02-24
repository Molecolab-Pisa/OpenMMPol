module mod_jacobian_mat
    use mod_memory, only: rp, ip
    implicit none
    private

    public :: Rij_jacobian

    contains
        
    subroutine Rij_jacobian(ci, cj, Rij, J_i, J_j)
        !! Compute the Jacobian matrix of distance 
        !! Rij = sqrt((ci(_x_)-cj(_x_))**2 + (ci(_y_)-cj(_y_))**2 + (ci(_z_)-cj(_z_))**2)
        !! Derivatives wrt ci(:) are saved in J_i and wrt cj(:) in J_j; the
        !! distance between the two points is also provided in output in Rij.

        implicit none

        real(rp), intent(in) :: ci(3), cj(3)
        real(rp), intent(out) :: Rij, J_i(3), J_j(3)
        integer(ip) :: i

        Rij = norm2(ci-cj)
        do i=1, 3
            J_i(i) = (ci(i) - cj(i))/Rij
            J_j(i) = -J_i(i)
        end do
    end subroutine

end module mod_jacobian_mat


