module mod_jacobian_mat
    use mod_memory, only: rp, ip
    implicit none
    private

    public :: Rij_jacobian
    public :: simple_angle_jacobian, inplane_angle_jacobian, &
              torsion_angle_jacobian, opb_angle_jacobian

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
        !! point (R) as central point. Then projects 
        !! \(J_r = \frac{\partial \theta}{\partial \vec{R}} \) onto A, B, C, and X 
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
        !!   J_r
        !! = \left( \frac{\partial \theta}{\partial \vec{A}} \right)_\vec{R} +
        !!   - \frac{\partial \vec{P}}{\partial \vec{A}} \times
        !!     \frac{\partial \hat{P}}{\partial \vec{P}} \times 
        !!     \frac{\partial (\vec{V} \cdot \hat{P})}{\partial \hat{P}} \times
        !!     J_r
        !! \]
        use mod_utils, only: cross_product, vec_skw, versor_der

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
        real(rp) :: k

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

        dpdpp = versor_der(cpp)

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

    subroutine torsion_angle_jacobian(ca, cb, cc, cd, thet, &
                                           J_a, J_b, J_c, J_d)
        !! Computes the Jacobian matrix for torsion angle defined by points 
        !! \(\vec{A}\), \(\vec{B}\), \(\vec{C}\) and \(\vec{D}\) (connected in
        !! this order). The angle is defined as follow:
        !! \[ \vec{U} = (\vec{B} - \vec{C}) \times (\vec{D} - \vec{C}) \\
        !!    \vec{T} = (\vec{B} - \vec{A}) \times (\vec{B} - \vec{C}) \\
        !!    cos(\theta) = \vec{U} \cdot \vec{T}
        !! \]
        !! Applying the chain rule:
        !! \[J_a = \frac{\partial \theta}{\partial \vec{A}} 
        !!       = -\frac{1}{\sqrt{1-\theta^2}} 
        !!          \frac{\partial cos(\theta)}{\partial \vec{A}} 
        !!       = -\frac{1}{\sqrt{1-\theta^2}}
        !!          \frac{\partial \vec{U} \cdot \vec{T}}{\partial \vec{A}} \\
        !!       = -\frac{1}{\sqrt{1-\theta^2}}
        !!          \frac{\partial (\vec{B} - \vec{A})}{\partial \vec{A}} \times
        !!          \frac{\partial \vec{T}}{\partial (\vec{B} - \vec{A})} \times
        !!          \frac{\partial \hat{T}}{\partial \vec{T}} \times \vec{U} \\
        !!       = -\frac{1}{\sqrt{1-\theta^2}}
        !!          (- \mathbb{I}_3) \times
        !!          skw(B-C) \times
        !!          \frac{\partial \hat{T}}{\partial \vec{T}} \times \vec{U}
        !! \]
        !! \[J_d = \frac{\partial \theta}{\partial \vec{D}} =
        !!        = -\frac{1}{\sqrt{1-\theta^2}}
        !!          \mathbb{I}_3 \times (-skw(B-C)) \times
        !!          \frac{\partial \hat{U}}{\partial \vec{U}} \times \vec{T}
        !! \]
        use mod_utils, only: cross_product, vec_skw, versor_der
        use mod_constants, only: pi, eps_rp

        implicit none

        real(rp), intent(inout), dimension(3) :: ca, cb, cc, cd
        !! Coordinates of the atoms defining the angle
        real(rp), intent(out), dimension(3) :: J_a, J_b, J_c, J_d
        !! The Jacobian components on atoms a, b and c and d respectively
        real(rp), intent(out) :: thet
        !! The torsion angle

        real(rp), dimension(3) :: a_b, c_d, c_b, t, u, ht, hu
        real(rp), dimension(3,3) :: dhudu, dhtdt, dtda, dudd, dhudd, dhtda, &
                                    dhudb, dhtdb, dudb, dtdb, &
                                    dhudc, dhtdc, dtdc, dudc
        real(rp) :: costhet, dacost, dd

        integer(ip) :: i, j
        
        a_b = cb - ca
        c_d = cd - cc
        c_b = cb - cc

        t = cross_product(a_b,c_b)
        ht = t / norm2(t)
        
        u = cross_product(c_b,c_d)
        hu = u / norm2(u)

        costhet = dot_product(hu,ht)
        if(costhet + 1.0 <= eps_rp) then
            thet = pi
        else 
            thet = acos(costhet)
        end if
        
        dacost = - 1.0/sqrt(1.0-costhet**2)
        
        dhtdt = versor_der(t)
        dhudu = versor_der(u)

        dtda = vec_skw(c_b)
        dudd = dtda 
        
        dhtda = matmul(dtda, dhtdt)
        J_a = -dacost * matmul(dhtda,hu) 

        dudb = vec_skw(c_d)
        dtdb = vec_skw(ca-cc)

        dhudb = matmul(dudb, dhudu)
        dhtdb = matmul(dtdb, dhtdt)
        J_b = dacost * (matmul(dhudb,ht) + matmul(dhtdb,hu))

        dtdc = vec_skw(a_b)
        dudc = vec_skw(cb - cd)

        dhtdb = matmul(dtdc, dhtdt)
        dhudc = matmul(dudc, dhudu) 
        J_c = dacost * (matmul(dhudc,ht) + matmul(dhtdb,hu))

        dhudd = matmul(dudd, dhudu)
        J_d = -dacost * matmul(dhudd,ht)

    end subroutine
    
    subroutine opb_angle_jacobian(ca, cb, cc, cd, thet, &
                                  J_a, J_b, J_c, J_d)
        use mod_utils, only: cross_product, vec_skw, versor_der
        use mod_constants, only: pi, eps_rp

        implicit none

        real(rp), intent(inout), dimension(3) :: ca, cb, cc, cd
        !! Coordinates of the atoms defining the angle
        real(rp), intent(out), dimension(3) :: J_a, J_b, J_c, J_d
        !! The Jacobian components on atoms a, b and c and d respectively
        real(rp), intent(out) :: thet
        !! The out-of-plane angle

        real(rp), dimension(3) :: a_b, a_c, v, p, hp, hv
        real(rp), dimension(3,3) :: dhpdp, dhvdv, dpda, dhpda, dvda, &
                                    dhvda, dhpdb, dpdb, dhpdc, dpdc, dvdd
        real(rp) :: costhet, dacost, dd, thet0

        integer(ip) :: i, j
        
        a_b = ca - cb
        a_c = ca - cc
        v = ca - cd

        p = cross_product(a_b,a_c)
        
        hp = p / norm2(p)
        hv = v / norm2(v)
        
        costhet = dot_product(hv, hp)
        if(costhet + 1.0 <= eps_rp) then
            thet0 = pi
        else 
            thet0 = acos(costhet)
        end if

        thet = abs(pi/2.0 - thet0)
        dacost = 1.0/sqrt(1.0-costhet**2)
        if(pi/2.0 - thet0 < 0.0) dacost = -dacost
        
        dhpdp = versor_der(p)
        dhvdv = versor_der(v)
        
        dpdb = vec_skw(a_c)
        dhpdb = matmul(dpdb,dhpdp)
        J_b = -dacost * matmul(dhpdb,hv)

        dpdc = -vec_skw(a_b)
        dhpdc = matmul(dpdc, dhpdp)
        J_c = -dacost * matmul(dhpdc,hv) 

        dvdd = -dhvdv
        J_d = dacost * matmul(dvdd,hp)
        
        J_a = -(J_b + J_c +J_d)
    end subroutine
end module mod_jacobian_mat
