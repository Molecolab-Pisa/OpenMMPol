module mod_jacobian_mat
    use mod_memory, only: rp, ip
    implicit none
    private

    public :: Rij_jacobian, simple_angle_jacobian

    contains
        
    pure subroutine Rij_jacobian(ci, cj, Rij, J_i, J_j)
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

    subroutine simple_angle_jacobian(ca, cb, cc, J_a, J_b, J_c)
        implicit none

        real(rp), intent(inout), dimension(3) :: ca, cb, cc
        real(rp), intent(out), dimension(3) :: J_a, J_b, J_c

        real(rp), dimension(3) :: dr1, dr2, J_anum
        real(rp) :: l1, l2, dr1_d_dr2, acosd
        integer(ip) :: i

        dr1 = ca - cb
        dr2 = cc - cb
        dr1_d_dr2 = dot_product(dr1, dr2)
        l1 = norm2(dr1)
        l2 = norm2(dr2)
        acosd = 1.0 / sqrt(1.0 - (dr1_d_dr2/(l1*l2))**2)

        do i=1,3
            J_a(i) = ((cc(i) - cb(i)) * l1 - (ca(i)-cb(i)) * dr1_d_dr2 / l1) / (l1**2*l2)
            J_c(i) = ((ca(i) - cb(i)) * l2 - (cc(i)-cb(i)) * dr1_d_dr2 / l2) / (l1*l2**2)
        end do

        J_a = J_a * acosd
        J_b = J_b * acosd
        J_c = J_c * acosd

        J_anum = 0.0
        do i=1,3
            !ca(i) = ca(i) + 1e-5
            cc(i) = cc(i) + 1e-5
            dr1 = ca - cb
            dr2 = cc - cb
            l1 = norm2(dr1)
            l2 = norm2(dr2)

            J_anum(i) = acos(dot_product(dr1, dr2)/(l1*l2)) 
            !ca(i) = ca(i) - 2e-5
            cc(i) = cc(i) - 2e-5
            dr1 = ca - cb
            dr2 = cc - cb
            l1 = norm2(dr1)
            l2 = norm2(dr2)
            J_anum(i) = J_anum(i) - acos(dot_product(dr1, dr2)/(l1*l2)) 
            J_anum(i) = J_anum(i) / 2e-5
            
            !ca(i) = ca(i) + 1e-5
            cc(i) = cc(i) + 1e-5
        end do

        write(*,*) J_c 
        write(*,*) J_anum

    end subroutine

end module mod_jacobian_mat


