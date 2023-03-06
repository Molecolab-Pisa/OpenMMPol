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

    pure subroutine simple_angle_jacobian(ca, cb, cc, thet, J_a, J_b, J_c)
        implicit none

        real(rp), intent(in), dimension(3) :: ca, cb, cc
        !! Coordinates of the atoms defining the angle
        real(rp), intent(out), dimension(3) :: J_a, J_b, J_c
        !! The Jacobian components on atoms a, b and c respectively
        real(rp), intent(out) :: thet
        !! The angle (in rad) defined by ca-cb-cc

        real(rp), dimension(3) :: dr1, dr2
        real(rp) :: l1, l2, dr1_d_dr2, acosd, cost
        integer(ip) :: i

        dr1 = ca - cb
        dr2 = cc - cb
        dr1_d_dr2 = dot_product(dr1, dr2)
        l1 = norm2(dr1)
        l2 = norm2(dr2)
        cost = dr1_d_dr2/(l1*l2)
        thet = acos(cost)
        acosd = 1.0 / sqrt(1.0 - cost**2)

        do i=1,3
            J_a(i) = -acosd * (dr2(i) * l1 - dr1(i) * dr1_d_dr2 / l1) / (l1**2*l2)
            J_c(i) = -acosd * (dr1(i) * l2 - dr2(i) * dr1_d_dr2 / l2) / (l1*l2**2)
            J_b(i) = -(J_a(i) + J_c(i))
        end do
    end subroutine
    
    pure subroutine inplane_angle_jacobian(ca, cb, cc, cx, thet, J_a, J_b, J_c, J_x)
        implicit none

        real(rp), intent(in), dimension(3) :: ca, cb, cc, cx
        !! Coordinates of the atoms defining the angle
        real(rp), intent(out), dimension(3) :: J_a, J_b, J_c, J_x
        !! The Jacobian components on atoms a, b and c respectively
        real(rp), intent(out) :: thet
        !! The angle (in rad) defined by ca-cb-cc

        real(rp), dimension(3) :: dr1, dr2
        real(rp) :: l1, l2, dr1_d_dr2, acosd, cost
        integer(ip) :: i

        dr1 = ca - cb
        dr2 = cc - cb
        dr1_d_dr2 = dot_product(dr1, dr2)
        l1 = norm2(dr1)
        l2 = norm2(dr2)
        cost = dr1_d_dr2/(l1*l2)
        thet = acos(cost)
        acosd = 1.0 / sqrt(1.0 - cost**2)

        do i=1,3
            J_a(i) = -acosd * (dr2(i) * l1 - dr1(i) * dr1_d_dr2 / l1) / (l1**2*l2)
            J_c(i) = -acosd * (dr1(i) * l2 - dr2(i) * dr1_d_dr2 / l2) / (l1*l2**2)
            J_b(i) = -(J_a(i) + J_c(i))
        end do
    end subroutine

end module mod_jacobian_mat


