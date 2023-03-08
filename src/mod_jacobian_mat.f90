module mod_jacobian_mat
    use mod_memory, only: rp, ip
    implicit none
    private

    public :: Rij_jacobian
    public :: simple_angle_jacobian, inplane_angle_jacobian

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
    
    pure subroutine inplane_angle_jacobian(ca, cb, cc, cx, thet, &
                                           J_a, J_b, J_c, J_x)
        !! Computes the Jacobian matrix for the inplane angle definition.
        !! It computes the Jacobian for the normal angle using the projected 
        !! point (R) as central point. Then projects \(J_r\) onto A, B, C, and X 
        !! (auxiliary point). The projection is done computing the 3x3 matrices
        !! of partial derivative of \(\vec{R}\) wrt any actual point and using
        !! them to project \(J_r\).
        !! \[\frac{\partial \vec{R}}{\partial \vec{A}} = 
        !!      \begin{bmatrix} 
        !!          \frac{\partial \vec{R}_x}{\partial \vec{A}_x} & 
        !!          \frac{\partial \vec{R}_y}{\partial \vec{A}_x} & 
        !!          \frac{\partial \vec{R}_z}{\partial \vec{A}_x} \\ 
        !!          \frac{\partial \vec{R}_x}{\partial \vec{A}_y} & 
        !!          \frac{\partial \vec{R}_y}{\partial \vec{A}_y} & 
        !!          \frac{\partial \vec{R}_z}{\partial \vec{A}_y} \\ 
        !!          \frac{\partial \vec{R}_x}{\partial \vec{A}_z} & 
        !!          \frac{\partial \vec{R}_y}{\partial \vec{A}_z} & 
        !!          \frac{\partial \vec{R}_z}{\partial \vec{A}_z} \\ 
        !!      \end{bmatrix} \]
        !! Those matrices are computed using the chain rule.
        !!
        !! Exemple:
        !! \[\vec{V} = \vec{B} - \vec{X}\]
        !! \[\vec{P} = (\vec{A} - \vec{X}) \times (\vec{C} - \vec{X})\]
        !! \[\hat{P} = \frac{\vec{P}}{||\vec{P}||}\]
        !! \[\vec{R} = \vec{B} - (\vec{V}\cdot\hat{P})\hat{P}\]
        !! \[J_a = \frac{\partial \theta}{\partial \vec{A}} 
        !! = \left( \frac{\partial \theta}{\partial \vec{A}} \right)_\vec{R} + 
        !!   \frac{\partial \vec{R}}{\partial \vec{A}} \times
        !!   \frac{\partial \theta}{\partial \vec{R}}
        !! = \left( \frac{\partial \theta}{\partial \vec{A}} \right)_\vec{R} +
        !!   - \frac{\partial \vec{P}}{\partial \vec{A}} \times
        !!     \frac{\partial \hat{P}}{\partial \vec{P}} \times 
        !!     \frac{\partial (\vec{V} \cdot \hat{P})}{\partial \hat{P}}
        !! \]
        use mod_utils, only: cross_product, vec_skw

        implicit none

        real(rp), intent(in), dimension(3) :: ca, cb, cc, cx
        !! Coordinates of the atoms defining the angle
        real(rp), intent(out), dimension(3) :: J_a, J_b, J_c, J_x
        !! The Jacobian components on atoms a, b and c respectively
        real(rp), intent(out) :: thet
        !! The angle (in rad) defined by ca-cb-cc

        real(rp), dimension(3) :: cr, cv, cp, cpp, cq, cs, J_r
        real(rp), dimension(3,3) :: drda, drdb, drdc, drdx, dkpdp, dppda, &
                                    dpdpp, dpda, dpdc, dppdc, dpdx, dppdx, &
                                    dkpdv
        real(rp) :: dd, k

        integer(ip) :: i, j

        cq = ca - cx
        cs = cc - cx
        cv = cb - cx

        cpp = cross_product(cq, cs)
        cp = cpp / norm2(cpp)

        k = dot_product(cv, cp)
        cr = cb - k * cp
        
        call simple_angle_jacobian(ca, cr, cc, thet, J_a, J_r, J_c)

        do i=1,3
            dkpdp(i,:) = cv(i) * cp
            dkpdp(i,i) = dkpdp(i,i) + k
        end do

        do i=1,3
            dpdpp(i,:) = -cpp(i) * cpp
            dpdpp(i,i) = dpdpp(i,i) + norm2(cpp)**2
        end do
        dpdpp = dpdpp / norm2(cpp)**3

        do i=1,3
            dkpdv(i,:) = cp(i) * cp
        end do

        dppda = -vec_skw(cs)
        dppdc = vec_skw(cq)
       
        dppdx = -vec_skw(cs-cq)

        dpda = matmul(dppda, dpdpp)
        drda = matmul(dpda, dkpdp)

        dpdc = matmul(dppdc, dpdpp)
        drdc = matmul(dpdc, dkpdp)

        dpdx = matmul(dppdx, dpdpp)
        drdx = -(matmul(dpdx, dkpdp) - dkpdv)

        J_a = J_a + matmul(drda, J_r)
        J_x = matmul(drdx, J_r)
        J_c = J_c + matmul(drdc, J_r)
        J_b = -(J_a+J_c+J_x)
        
    end subroutine

end module mod_jacobian_mat


