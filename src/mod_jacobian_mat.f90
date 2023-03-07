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
    
    subroutine inplane_angle_jacobian(ca, cb, cc, cx, thet, J_a, J_b, J_c, J_x)
        use mod_utils, only: cross_product

        implicit none

        real(rp), intent(inout), dimension(3) :: ca, cb, cc, cx
        !! Coordinates of the atoms defining the angle
        real(rp), intent(out), dimension(3) :: J_a, J_b, J_c, J_x
        !! The Jacobian components on atoms a, b and c respectively
        real(rp), intent(out) :: thet
        !! The angle (in rad) defined by ca-cb-cc

        real(rp), dimension(3) :: cr, cv, cp, cpp, cq, cs, J_r
        real(rp), dimension(3,3) :: drda, drdb, drdc, drdx
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

        drda = 0.0
        drdb = 0.0
        drdc = 0.0
        drdx = 0.0

        ! Numerical
        dd = 1e-5
        do i=1, 3
            ca(i) = ca(i) + dd
            
            cq = ca - cx
            cs = cc - cx
            cv = cb - cx

            cpp = cross_product(cq, cs)
            cp = cpp / norm2(cpp)

            k = dot_product(cv, cp)
            cr = cb - k * cp
            !drda(i,:) = cr
            call simple_angle_jacobian(ca, cr, cc, thet, J_r, J_r, J_r)
            J_a(i) = thet
            
            ca(i) = ca(i) - 2*dd
            
            cq = ca - cx
            cs = cc - cx
            cv = cb - cx

            cpp = cross_product(cq, cs)
            cp = cpp / norm2(cpp)

            k = dot_product(cv, cp)
            cr = cb - k * cp
            call simple_angle_jacobian(ca, cr, cc, thet, J_r, J_r, J_r)
            !drda(i,:) = drda(i,:) - cr
            !drda(i,:) = drda(i,:) / (2*dd)
            J_a(i) = (J_a(i) - thet)/(2*dd)

            ca(i) = ca(i) + dd
        end do

        do i=1, 3
            cb(i) = cb(i) + dd
            
            cq = ca - cx
            cs = cc - cx
            cv = cb - cx

            cpp = cross_product(cq, cs)
            cp = cpp / norm2(cpp)

            k = dot_product(cv, cp)
            cr = cb - k * cp
            !drdb(i,:) = cr
            call simple_angle_jacobian(ca, cr, cc, thet, J_r, J_r, J_r)
            J_b(i) = thet
            
            cb(i) = cb(i) - 2*dd
            
            cq = ca - cx
            cs = cc - cx
            cv = cb - cx

            cpp = cross_product(cq, cs)
            cp = cpp / norm2(cpp)

            k = dot_product(cv, cp)
            cr = cb - k * cp
            call simple_angle_jacobian(ca, cr, cc, thet, J_r, J_r, J_r)
            !drdb(i,:) = drdb(i,:) - cr
            !drdb(i,:) = drdb(i,:) / (2*dd)
            J_b(i) = (J_b(i) - thet)/(2*dd)

            cb(i) = cb(i) + dd
        end do
        
        do i=1, 3
            cc(i) = cc(i) + dd
            
            cq = ca - cx
            cs = cc - cx
            cv = cb - cx

            cpp = cross_product(cq, cs)
            cp = cpp / norm2(cpp)

            k = dot_product(cv, cp)
            cr = cb - k * cp
            !drdc(i,:) = cr
            call simple_angle_jacobian(ca, cr, cc, thet, J_r, J_r, J_r)
            J_c(i) = thet
            
            cc(i) = cc(i) - 2*dd
            
            cq = ca - cx
            cs = cc - cx
            cv = cb - cx

            cpp = cross_product(cq, cs)
            cp = cpp / norm2(cpp)

            k = dot_product(cv, cp)
            cr = cb - k * cp
            call simple_angle_jacobian(ca, cr, cc, thet, J_r, J_r, J_r)
            !drdc(i,:) = drdc(i,:) - cr
            !drdc(i,:) = drdc(i,:) / (2*dd)
            J_c(i) = (J_c(i) - thet)/(2*dd)

            cc(i) = cc(i) + dd
        end do
        
        do i=1, 3
            cx(i) = cx(i) + dd
            
            cq = ca - cx
            cs = cc - cx
            cv = cb - cx

            cpp = cross_product(cq, cs)
            cp = cpp / norm2(cpp)

            k = dot_product(cv, cp)
            cr = cb - k * cp
            !drdb(i,:) = cr
            call simple_angle_jacobian(ca, cr, cc, thet, J_r, J_r, J_r)
            J_x(i) = thet
            
            cx(i) = cx(i) - 2*dd
            
            cq = ca - cx
            cs = cc - cx
            cv = cb - cx

            cpp = cross_product(cq, cs)
            cp = cpp / norm2(cpp)

            k = dot_product(cv, cp)
            cr = cb - k * cp
            call simple_angle_jacobian(ca, cr, cc, thet, J_r, J_r, J_r)
            !drdb(i,:) = drdb(i,:) - cr
            !drdb(i,:) = drdb(i,:) / (2*dd)
            J_x(i) = (J_x(i) - thet)/(2*dd)

            cx(i) = cx(i) + dd
        end do
        
        !J_a = J_a + matmul(drda, J_r)
        !J_b = matmul(drdb, J_r)
        !J_x = matmul(drdx, J_r)
        !J_c = J_c + matmul(drdc, J_r)
    end subroutine

end module mod_jacobian_mat


