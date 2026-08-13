module mod_jacobian_mat
    use mod_memory, only: rp, ip
    implicit none
    private

    public :: Rij_jacobian
    public :: simple_angle_jacobian, inplane_angle_jacobian, &
              torsion_angle_jacobian, opb_angle_jacobian, &
              pitors_angle_jacobian
    public :: Rij_hessian
    public :: simple_angle_hessian
    public :: inplane_angle_hessian
    public :: opb_angle_hessian
    public :: torsion_angle_hessian
    public :: pitors_angle_hessian

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

    pure subroutine Rij_hessian(ci, cj, Rij, J_i, J_j, H_ii, H_ij, H_jj)
        !! Compute the Hessian matrix of the distance Rij (see Rij_jacobian).
        !! With \(\hat u = (\vec{c}_i-\vec{c}_j)/R_{ij}\) and
        !! \(M = \frac{1}{R_{ij}}(\mathbb{1}_3 - \hat u\hat u^\dagger)\)
        !! (i.e. \code{versor\_der(ci-cj)}):
        !! \[ \frac{\partial^2 R_{ij}}{\partial \vec c_i \partial \vec c_i} =
        !!    \frac{\partial^2 R_{ij}}{\partial \vec c_j \partial \vec c_j} = M,
        !!    \qquad
        !!    \frac{\partial^2 R_{ij}}{\partial \vec c_i \partial \vec c_j} = -M \]
        use mod_utils, only: versor_der

        implicit none

        real(rp), intent(in) :: ci(3), cj(3)
        real(rp), intent(out) :: Rij, J_i(3), J_j(3)
        real(rp), intent(out) :: H_ii(3,3), H_ij(3,3), H_jj(3,3)

        call Rij_jacobian(ci, cj, Rij, J_i, J_j)
        H_ii = versor_der(ci-cj)
        H_jj = H_ii
        H_ij = -H_ii
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

    pure subroutine simple_angle_hessian(ca, cb, cc, thet, J_a, J_b, J_c, &
                                         H_aa, H_ab, H_ac, H_bb, H_bc, H_cc)
        !! Compute the Hessian matrix of the valence angle A-B-C (see
        !! simple_angle_jacobian). Returns the six independent 3x3 blocks of
        !! the full (9x9) Hessian; the missing three follow by symmetry of
        !! the Hessian of a scalar, e.g. d^2(thet)/dB dA = transpose(H_ab).
        !!
        !! With \(\bd_1=\bA-\bB\), \(\bd_2=\bC-\bB\), \(\hn{\bd}_1,\hn{\bd}_2\)
        !! their versors, \(c=\cos\theta=\hn{\bd}_1\cdot\hn{\bd}_2\),
        !! \(s=\sin\theta\), \(V_k = \partial\hn{\bd}_k/\partial\bd_k\)
        !! (\code{versor\_der}) and \(\Xi_1=\Xi(\bd_1,\hn{\bd}_2)\),
        !! \(\Xi_2=\Xi(\bd_2,\hn{\bd}_1)\) (\code{versor\_der2}), the Hessian
        !! of \(\theta\) wrt any two atoms X,Y is
        !! \(H_{XY} = -\frac{1}{s}\left(\partial^2 c/\partial X \partial Y +
        !! c\, J_X J_Y^\dagger\right)\), with
        !! \[ \frac{\partial^2 c}{\partial \vec A\partial \vec A} = \Xi_1, \quad
        !!    \frac{\partial^2 c}{\partial \vec A\partial \vec C} = V_1 V_2, \quad
        !!    \frac{\partial^2 c}{\partial \vec A\partial \vec B} = -(\Xi_1+V_1V_2), \]
        !! \[ \frac{\partial^2 c}{\partial \vec C\partial \vec C} = \Xi_2, \quad
        !!    \frac{\partial^2 c}{\partial \vec B\partial \vec C} = -(V_1V_2+\Xi_2), \quad
        !!    \frac{\partial^2 c}{\partial \vec B\partial \vec B} = \Xi_1+\Xi_2+V_1V_2+V_2V_1. \]
        use mod_utils, only: versor_der, versor_der2

        implicit none

        real(rp), intent(in), dimension(3) :: ca, cb, cc
        !! Coordinates of the atoms defining the angle
        real(rp), intent(out), dimension(3) :: J_a, J_b, J_c
        !! The Jacobian components on atoms a, b and c respectively
        real(rp), intent(out) :: thet
        !! The angle (in rad) defined by ca-cb-cc
        real(rp), intent(out), dimension(3,3) :: H_aa, H_ab, H_ac, H_bb, H_bc, H_cc
        !! The independent 3x3 blocks of the Hessian of thet

        real(rp), dimension(3) :: dr1, dr2, d1h, d2h
        real(rp) :: l1, l2, cost, sint
        real(rp), dimension(3,3) :: V1, V2, V1V2, V2V1, Xi1, Xi2, &
                                    c2_aa, c2_ab, c2_ac, c2_bb, c2_bc, c2_cc
        integer(ip) :: k

        call simple_angle_jacobian(ca, cb, cc, thet, J_a, J_b, J_c)

        dr1 = ca - cb
        dr2 = cc - cb
        l1 = norm2(dr1)
        l2 = norm2(dr2)
        d1h = dr1 / l1
        d2h = dr2 / l2
        cost = dot_product(d1h, d2h)
        sint = sqrt(1.0_rp - cost**2)

        V1 = versor_der(dr1)
        V2 = versor_der(dr2)
        V1V2 = matmul(V1, V2)
        V2V1 = matmul(V2, V1)
        Xi1 = versor_der2(dr1, d2h)
        Xi2 = versor_der2(dr2, d1h)

        c2_aa = Xi1
        c2_ac = V1V2
        c2_ab = -(Xi1 + V1V2)
        c2_cc = Xi2
        c2_bc = -(V1V2 + Xi2)
        c2_bb = Xi1 + Xi2 + V1V2 + V2V1

        do k=1,3
            H_aa(k,:) = c2_aa(k,:) + cost*J_a(k)*J_a
            H_ac(k,:) = c2_ac(k,:) + cost*J_a(k)*J_c
            H_ab(k,:) = c2_ab(k,:) + cost*J_a(k)*J_b
            H_cc(k,:) = c2_cc(k,:) + cost*J_c(k)*J_c
            H_bc(k,:) = c2_bc(k,:) + cost*J_b(k)*J_c
            H_bb(k,:) = c2_bb(k,:) + cost*J_b(k)*J_b
        end do
        H_aa = -H_aa/sint
        H_ac = -H_ac/sint
        H_ab = -H_ab/sint
        H_cc = -H_cc/sint
        H_bc = -H_bc/sint
        H_bb = -H_bb/sint
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
        real(rp), dimension(3,3) :: drda, drdc, drdx, dkpdp, dppda, &
                                    dpdpp, dpda, dpdc, dppdc, dpdx, dppdx, &
                                    dkpdv
        real(rp) :: k
        integer(ip) :: i

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

    pure subroutine inplane_angle_hessian(ca, cb, cc, cx, thet, J_a, J_b, J_c, J_x, &
                                     H_aa, H_ab, H_ac, H_ax, H_bb, H_bc, H_bx, &
                                     H_cc, H_cx, H_xx)
        !! Hessian of the in-plane angle of inplane_angle_jacobian. Reuses
        !! the notation there: \(\bq=\bA-\bX\), \(\bs=\bC-\bX\), \(\bv=\bB-\bX\),
        !! \(\bP=\bq\times\bs\), \(\hn\bP=\bP/|\bP|\), \(k=\bv\cdot\hn\bP\),
        !! \(\bR=\bB-k\hn\bP\), \(\theta=\angle(\bA,\bR,\bC)\).
        !!
        !! Writing \(\theta=\theta^{(0)}(\bA,\bR(\bA,\bB,\bC,\bX),\bC)\), the
        !! Hessian follows from the ordinary two-level chain rule, with
        !! \(\theta^{(0)}\)'s own Hessian blocks (\(H^{(0)}\)) coming from
        !! simple_angle_hessian(A,R,C) (with R playing the role of the
        !! vertex), \(D_Y \equiv \partial\bR/\partial Y\) for
        !! \(Y\in\{\bA,\bB,\bC,\bX\}\), and \(H^\phi_{YZ}\equiv
        !! \sum_n (J_R)_n\, \partial^2 R_n/\partial Y\partial Z\) (the part of
        !! the chain rule coming from R's own curvature, i.e. the Hessian of
        !! the scalar \(\phi=\bR\cdot\bw\) at fixed \(\bw=J_R\)):
        !! \[ H_{YZ} = \left(\frac{\partial \bu}{\partial Y}\right)^\dagger
        !!    H^{(0)}\left(\frac{\partial \bu}{\partial Z}\right) +
        !!    H^\phi_{YZ}, \qquad \bu=(\bA,\bR,\bC) \]
        !! (e.g. \(H_{AA}=H^{(0)}_{AA}+H^{(0)}_{AR}D_A+D_A^\dagger(H^{(0)}_{AR})^\dagger
        !! +D_A^\dagger H^{(0)}_{RR}D_A+H^\phi_{AA}\),
        !! \(H_{BX}=D_B^\dagger H^{(0)}_{RR}D_X+H^\phi_{BX}\), etc).
        !!
        !! \(\phi\) itself is \(\bB\cdot\bw-k\,m\) with \(m=\bw\cdot\hn\bP\), so
        !! \(H^\phi\) reduces to derivatives of \(k\) and \(m\). Both, in turn,
        !! reduce to two pieces: \((\partial\bP/\partial Y)^\dagger\Xi(\bP,\bw')
        !! (\partial\bP/\partial Z)\) (\(\Xi\equiv\)versor_der2, \(\bw'=\bw\) for
        !! \(m\) or \(\bw'=\bv\) for \(k\)) plus the Hessian of the scalar
        !! triple product \(\bP\cdot g(\bP)\bw'\) at fixed weight \(g(\bP)\bw'\)
        !! (\(g\equiv\)versor_der), which is exactly \(\pm[g(\bP)\bw']_\times\)
        !! on the A-C/A-X/C-X blocks and exactly 0 on the A-A/C-C/X-X blocks,
        !! since \(\bP\) is a bilinear function of \(\bA,\bC,\bX\) (only the
        !! first, \(\Xi\)-based piece survives on the A-A/C-C/X-X blocks).
        use mod_utils, only: cross_product, vec_skw, versor_der, versor_der2

        implicit none

        real(rp), intent(in), dimension(3) :: ca, cb, cc, cx
        real(rp), intent(out), dimension(3) :: J_a, J_b, J_c, J_x
        real(rp), intent(out) :: thet
        real(rp), intent(out), dimension(3,3) :: H_aa, H_ab, H_ac, H_ax, &
                                                 H_bb, H_bc, H_bx, &
                                                 H_cc, H_cx, H_xx

        real(rp), dimension(3) :: q, s, v, pp, phat, r, J_r, &
                                  dk_da, dk_dc, dk_dx, dk_db, &
                                  dm_da, dm_dc, dm_dx
        real(rp), dimension(3,3) :: Haa0, Har0, Hac0, Hrr0, Hrc0, Hcc0, &
                                    Pda, Pdc, Pdx, Phda, Phdc, Phdx, &
                                    Rda, Rdb, Rdc, Rdx, &
                                    Xiw, Xiv, skw_w, skw_v, &
                                    Hm_aa, Hm_ac, Hm_ax, Hm_cc, Hm_cx, Hm_xx, &
                                    Hk_aa, Hk_ab, Hk_ac, Hk_ax, Hk_cc, Hk_bc, Hk_bx, Hk_cx, Hk_xx, &
                                    Hphi_aa, Hphi_ab, Hphi_ac, Hphi_ax, &
                                    Hphi_bb, Hphi_bc, Hphi_bx, &
                                    Hphi_cc, Hphi_cx, Hphi_xx
        real(rp) :: k, m
        integer(ip) :: i

        q = ca - cx
        s = cc - cx
        v = cb - cx

        pp = cross_product(q, s)
        phat = pp / norm2(pp)
        k = dot_product(v, phat)
        r = cb - k*phat

        call simple_angle_hessian(ca, r, cc, thet, J_a, J_r, J_c, &
                                  Haa0, Har0, Hac0, Hrr0, Hrc0, Hcc0)
        ! J_a, J_c above are (d thet/dA)|_R, (d thet/dC)|_R at this point;
        ! corrected below to the full derivatives.

        ! dP/dY blocks (standard convention: row=P-component, col=Y-component)
        Pda = -vec_skw(s)
        Pdc = vec_skw(q)
        Pdx = vec_skw(s-q)

        ! dPhat/dY = versor_der(P) * dP/dY
        Phda = matmul(versor_der(pp), Pda)
        Phdc = matmul(versor_der(pp), Pdc)
        Phdx = matmul(versor_der(pp), Pdx)

        dk_da = matmul(transpose(Phda), v)
        dk_dc = matmul(transpose(Phdc), v)
        dk_dx = matmul(transpose(Phdx), v) - phat
        dk_db = phat

        dm_da = matmul(transpose(Phda), J_r)
        dm_dc = matmul(transpose(Phdc), J_r)
        dm_dx = matmul(transpose(Phdx), J_r)
        m = dot_product(phat, J_r)

        ! dR/dY blocks
        do i=1,3
            Rda(i,:) = -phat(i)*dk_da - k*Phda(i,:)
            Rdc(i,:) = -phat(i)*dk_dc - k*Phdc(i,:)
            Rdx(i,:) = -phat(i)*dk_dx - k*Phdx(i,:)
            Rdb(i,:) = -phat(i)*phat
        end do
        do i=1,3
            Rdb(i,i) = Rdb(i,i) + 1.0_rp
        end do

        ! Full Jacobian (matches inplane_angle_jacobian)
        J_a = J_a + matmul(transpose(Rda), J_r)
        J_c = J_c + matmul(transpose(Rdc), J_r)
        J_x = matmul(transpose(Rdx), J_r)
        J_b = matmul(transpose(Rdb), J_r)

        ! H^m: Hessian of the scalar m = J_r . Phat(A,C,X), weight J_r fixed.
        ! Any block touching B is 0 (Phat has no B-dependence).
        Xiw = versor_der2(pp, J_r)
        skw_w = vec_skw(matmul(versor_der(pp), J_r))
        Hm_aa = matmul(transpose(Pda), matmul(Xiw, Pda))
        Hm_ac = matmul(transpose(Pda), matmul(Xiw, Pdc)) - skw_w
        Hm_ax = matmul(transpose(Pda), matmul(Xiw, Pdx)) + skw_w
        Hm_cc = matmul(transpose(Pdc), matmul(Xiw, Pdc))
        Hm_cx = matmul(transpose(Pdc), matmul(Xiw, Pdx)) - skw_w
        Hm_xx = matmul(transpose(Pdx), matmul(Xiw, Pdx))

        ! H^mu: same as H^m but with weight v=B-X (current value) instead of
        ! J_r; used only to build H^k below.
        Xiv = versor_der2(pp, v)
        skw_v = vec_skw(matmul(versor_der(pp), v))

        ! H^k: Hessian of k = v.Phat(A,B,C,X) (v=B-X depends on B,X too)
        Hk_aa = matmul(transpose(Pda), matmul(Xiv, Pda))
        Hk_ab = transpose(Phda)
        Hk_ac = matmul(transpose(Pda), matmul(Xiv, Pdc)) - skw_v
        Hk_ax = -transpose(Phda) + matmul(transpose(Pda), matmul(Xiv, Pdx)) + skw_v
        Hk_cc = matmul(transpose(Pdc), matmul(Xiv, Pdc))
        Hk_bc = Phdc
        Hk_bx = Phdx
        Hk_cx = matmul(transpose(Pdc), matmul(Xiv, Pdx)) - skw_v &
                - transpose(Phdc)
        Hk_xx = -Phdx - transpose(Phdx) &
                + matmul(transpose(Pdx), matmul(Xiv, Pdx))

        ! H^phi: Hessian of phi = R.J_r = B.J_r - k*m (J_r fixed weight).
        ! Every block is -(m*Hk_XY + outer(dk_X,dm_Y) + outer(dm_X,dk_Y) +
        ! k*Hm_XY), with dk_B=phat, dm_B=0 (m has no B-dependence) and
        ! Hm_AB=Hm_BB=Hm_BC=Hm_BX=0 (Phat has no B-dependence).
        do i=1,3
            Hphi_aa(i,:) = -(m*Hk_aa(i,:) + dk_da(i)*dm_da + dm_da(i)*dk_da + k*Hm_aa(i,:))
            Hphi_cc(i,:) = -(m*Hk_cc(i,:) + dk_dc(i)*dm_dc + dm_dc(i)*dk_dc + k*Hm_cc(i,:))
            Hphi_bb(i,:) = 0.0_rp
            Hphi_ab(i,:) = -(m*Hk_ab(i,:) + dm_da(i)*phat)
            Hphi_ac(i,:) = -(m*Hk_ac(i,:) + dk_da(i)*dm_dc + dm_da(i)*dk_dc + k*Hm_ac(i,:))
            Hphi_ax(i,:) = -(m*Hk_ax(i,:) + dk_da(i)*dm_dx + dm_da(i)*dk_dx + k*Hm_ax(i,:))
            Hphi_bc(i,:) = -(m*Hk_bc(i,:) + phat(i)*dm_dc)
            Hphi_bx(i,:) = -(m*Hk_bx(i,:) + phat(i)*dm_dx)
            Hphi_cx(i,:) = -(m*Hk_cx(i,:) + dk_dc(i)*dm_dx + dm_dc(i)*dk_dx + k*Hm_cx(i,:))
            Hphi_xx(i,:) = -(m*Hk_xx(i,:) + dk_dx(i)*dm_dx + dm_dx(i)*dk_dx + k*Hm_xx(i,:))
        end do

        ! Final assembly: theta = theta0(A, R(A,B,C,X), C)
        H_aa = Haa0 + matmul(Har0,Rda) + matmul(transpose(Rda),transpose(Har0)) &
               + matmul(transpose(Rda),matmul(Hrr0,Rda)) + Hphi_aa
        H_ac = Hac0 + matmul(Har0,Rdc) + matmul(transpose(Rda),Hrc0) &
               + matmul(transpose(Rda),matmul(Hrr0,Rdc)) + Hphi_ac
        H_ab = matmul(Har0,Rdb) + matmul(transpose(Rda),matmul(Hrr0,Rdb)) + Hphi_ab
        H_ax = matmul(Har0,Rdx) + matmul(transpose(Rda),matmul(Hrr0,Rdx)) + Hphi_ax
        H_cc = Hcc0 + matmul(transpose(Hrc0),Rdc) + matmul(transpose(Rdc),Hrc0) &
               + matmul(transpose(Rdc),matmul(Hrr0,Rdc)) + Hphi_cc
        H_bc = matmul(transpose(Rdb),Hrc0) + matmul(transpose(Rdb),matmul(Hrr0,Rdc)) &
               + Hphi_bc
        H_cx = matmul(transpose(Hrc0),Rdx) + matmul(transpose(Rdc),matmul(Hrr0,Rdx)) &
               + Hphi_cx
        H_bb = matmul(transpose(Rdb),matmul(Hrr0,Rdb)) + Hphi_bb
        H_bx = matmul(transpose(Rdb),matmul(Hrr0,Rdx)) + Hphi_bx
        H_xx = matmul(transpose(Rdx),matmul(Hrr0,Rdx)) + Hphi_xx

    end subroutine inplane_angle_hessian

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
                                    dhudc, dtdc, dudc
        real(rp) :: costhet, dacost, s

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
            dacost = -1.0 / sin(costhet)
        else if(abs(costhet - 1.0) <= eps_rp) then
            thet = 0.0
            dacost = -1.0 / sin(costhet)
        else
            thet = acos(costhet)
            dacost = - 1.0/sqrt(1.0-costhet**2)
        end if

        s = sign(1.0_rp, -dot_product(a_b, u))
        thet = thet * s
        
        dacost = dacost * s
        
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

    pure subroutine torsion_angle_hessian(ca, cb, cc, cd, thet, J_a, J_b, J_c, J_d, &
                                          H_aa, H_ab, H_ac, H_ad, H_bb, H_bc, H_bd, &
                                          H_cc, H_cd, H_dd)
        !! Hessian of the torsion (dihedral) angle of torsion_angle_jacobian.
        !! With \(\ba_b=\bB-\bA\), \(\bc_b=\bB-\bC\), \(\bc_d=\bD-\bC\),
        !! \(\bT=\ba_b\times\bc_b\), \(\bU=\bc_b\times\bc_d\), \(c=\cos\theta_0=
        !! \hn\bU\cdot\hn\bT\), \(s_0=\sin\theta_0\), \(\sigma=\mathrm{sign}
        !! (-\ba_b\cdot\bU)\) (so \(\theta=\sigma\theta_0\) and
        !! \(J_X=-(\sigma/s_0)\,\partial c/\partial X\)):
        !! \[ H_{XY} = -\frac{\sigma c}{s_0}J_XJ_Y^\dagger -
        !!    \frac{\sigma}{s_0}\frac{\partial^2c}{\partial X\partial Y} \]
        !! \(\partial^2c/\partial X\partial Y\) follows the same recipe as
        !! opb_angle_hessian, except now BOTH \(\bT\) and \(\bU\) are curved
        !! (bilinear) functions of the atoms (unlike OPB's \(\bv\), which was
        !! linear), so there are two independent triple-product corrections
        !! \(\psi_T\) (weight \(\bw_T=V(\bT)\hn\bU\)) and \(\psi_U\) (weight
        !! \(\bw_U=V(\bU)\hn\bT\)), each \(\pm[\bw]_\times\) on the three
        !! blocks pairing the two atoms \(\bT\) (resp. \(\bU\)) depends on,
        !! 0 elsewhere:
        !! \[ (\psi_T)_{AB}=[\bw_T]_\times,\ (\psi_T)_{AC}=-[\bw_T]_\times,\
        !!    (\psi_T)_{BC}=[\bw_T]_\times,\ (\psi_T)_{AA}=(\psi_T)_{BB}=
        !!    (\psi_T)_{CC}=0,\ \text{0 on every D-block} \]
        !! \[ (\psi_U)_{BC}=[\bw_U]_\times,\ (\psi_U)_{BD}=-[\bw_U]_\times,\
        !!    (\psi_U)_{CD}=[\bw_U]_\times,\ (\psi_U)_{BB}=(\psi_U)_{CC}=
        !!    (\psi_U)_{DD}=0,\ \text{0 on every A-block} \]
        !! with \(T_A=[\bc_b]_\times\), \(T_B=[\bC-\bA]_\times\),
        !! \(T_C=-[\ba_b]_\times\), \(T_D=0\), \(U_A=0\),
        !! \(U_B=-[\bc_d]_\times\), \(U_C=[\bD-\bB]_\times\),
        !! \(U_D=[\bc_b]_\times\).
        use mod_utils, only: cross_product, vec_skw, versor_der, versor_der2
        use mod_constants, only: pi, eps_rp

        implicit none

        real(rp), intent(in), dimension(3) :: ca, cb, cc, cd
        real(rp), intent(out), dimension(3) :: J_a, J_b, J_c, J_d
        real(rp), intent(out) :: thet
        real(rp), intent(out), dimension(3,3) :: H_aa, H_ab, H_ac, H_ad, &
                                                 H_bb, H_bc, H_bd, &
                                                 H_cc, H_cd, H_dd

        real(rp), dimension(3) :: a_b, c_b, c_d, t, u, ht, hu, gt, gu, &
                                  dc_da, dc_db, dc_dc, dc_dd
        real(rp), dimension(3,3) :: Ta, Tb, Tc, Ub, Uc, Ud, &
                                    Vt, Vu, XiTU, XiUT, VtVu, skw_wt, skw_wu
        real(rp) :: costhet, s0, sigma, sigfac, pref
        integer(ip) :: i

        a_b = cb - ca
        c_d = cd - cc
        c_b = cb - cc

        t = cross_product(a_b, c_b)
        ht = t / norm2(t)
        u = cross_product(c_b, c_d)
        hu = u / norm2(u)

        costhet = dot_product(hu, ht)
        if(costhet + 1.0 <= eps_rp) then
            thet = pi
        else if(abs(costhet - 1.0) <= eps_rp) then
            thet = 0.0_rp
        else
            thet = acos(costhet)
        end if
        sigma = sign(1.0_rp, -dot_product(a_b, u))
        thet = thet * sigma

        s0 = sqrt(max(1.0_rp - costhet**2, eps_rp))
        sigfac = -sigma / s0
        pref = -sigma * costhet / s0

        Vt = versor_der(t)
        Vu = versor_der(u)
        gt = matmul(Vt, hu)
        gu = matmul(Vu, ht)

        Ta = vec_skw(c_b)
        Tb = vec_skw(cc - ca)
        Tc = -vec_skw(a_b)
        Ub = -vec_skw(c_d)
        Uc = vec_skw(cd - cb)
        Ud = vec_skw(c_b)

        dc_da = matmul(transpose(Ta), gt)
        dc_db = matmul(transpose(Tb), gt) + matmul(transpose(Ub), gu)
        dc_dc = matmul(transpose(Tc), gt) + matmul(transpose(Uc), gu)
        dc_dd = matmul(transpose(Ud), gu)

        J_a = sigfac * dc_da
        J_b = sigfac * dc_db
        J_c = sigfac * dc_dc
        J_d = sigfac * dc_dd

        XiTU = versor_der2(t, hu)
        XiUT = versor_der2(u, ht)
        VtVu = matmul(Vt, Vu)
        skw_wt = vec_skw(gt)
        skw_wu = vec_skw(gu)

        do i=1,3
            H_aa(i,:) = pref*J_a(i)*J_a
            H_ab(i,:) = pref*J_a(i)*J_b
            H_ac(i,:) = pref*J_a(i)*J_c
            H_ad(i,:) = pref*J_a(i)*J_d
            H_bb(i,:) = pref*J_b(i)*J_b
            H_bc(i,:) = pref*J_b(i)*J_c
            H_bd(i,:) = pref*J_b(i)*J_d
            H_cc(i,:) = pref*J_c(i)*J_c
            H_cd(i,:) = pref*J_c(i)*J_d
            H_dd(i,:) = pref*J_d(i)*J_d
        end do

        H_aa = H_aa + sigfac*( matmul(transpose(Ta),matmul(XiTU,Ta)) )
        H_ab = H_ab + sigfac*( matmul(transpose(Ta),matmul(XiTU,Tb)) &
                              + matmul(transpose(Ta),matmul(VtVu,Ub)) + skw_wt )
        H_ac = H_ac + sigfac*( matmul(transpose(Ta),matmul(XiTU,Tc)) &
                              + matmul(transpose(Ta),matmul(VtVu,Uc)) - skw_wt )
        H_ad = H_ad + sigfac*( matmul(transpose(Ta),matmul(VtVu,Ud)) )
        H_bb = H_bb + sigfac*( matmul(transpose(Tb),matmul(XiTU,Tb)) &
                              + matmul(transpose(Tb),matmul(VtVu,Ub)) &
                              + matmul(transpose(Ub),matmul(transpose(VtVu),Tb)) &
                              + matmul(transpose(Ub),matmul(XiUT,Ub)) )
        H_bc = H_bc + sigfac*( matmul(transpose(Tb),matmul(XiTU,Tc)) &
                              + matmul(transpose(Tb),matmul(VtVu,Uc)) &
                              + matmul(transpose(Ub),matmul(transpose(VtVu),Tc)) &
                              + matmul(transpose(Ub),matmul(XiUT,Uc)) + skw_wt + skw_wu )
        H_bd = H_bd + sigfac*( matmul(transpose(Tb),matmul(VtVu,Ud)) &
                              + matmul(transpose(Ub),matmul(XiUT,Ud)) - skw_wu )
        H_cc = H_cc + sigfac*( matmul(transpose(Tc),matmul(XiTU,Tc)) &
                              + matmul(transpose(Tc),matmul(VtVu,Uc)) &
                              + matmul(transpose(Uc),matmul(transpose(VtVu),Tc)) &
                              + matmul(transpose(Uc),matmul(XiUT,Uc)) )
        H_cd = H_cd + sigfac*( matmul(transpose(Tc),matmul(VtVu,Ud)) &
                              + matmul(transpose(Uc),matmul(XiUT,Ud)) + skw_wu )
        H_dd = H_dd + sigfac*( matmul(transpose(Ud),matmul(XiUT,Ud)) )

    end subroutine torsion_angle_hessian
    
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
        real(rp), dimension(3,3) :: dhpdp, dhvdv, &
                                    dhpdb, dpdb, dhpdc, dpdc, dvdd
        real(rp) :: costhet, dacost, thet0

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

        dpdc = vec_skw(a_b)
        dpdc = -dpdc
        dhpdc = matmul(dpdc, dhpdp)
        J_c = -dacost * matmul(dhpdc,hv) 

        dvdd = -dhvdv
        J_d = dacost * matmul(dvdd,hp)
        
        J_a = -(J_b + J_c +J_d)
    end subroutine

    pure subroutine opb_angle_hessian(ca, cb, cc, cd, thet, J_a, J_b, J_c, J_d, &
                                      H_aa, H_ab, H_ac, H_ad, H_bb, H_bc, H_bd, &
                                      H_cc, H_cd, H_dd)
        !! Hessian of the out-of-plane (Allinger) angle of
        !! opb_angle_jacobian. With \(\ba_b=\bA-\bB\), \(\ba_c=\bA-\bC\),
        !! \(\bv=\bA-\bD\), \(\bP=\ba_b\times\ba_c\), \(c=\cos\theta_0=
        !! \hn\bv\cdot\hn\bP\), \(s=\sin\theta_0\), \(\sigma=
        !! \mathrm{sign}(\pi/2-\theta_0)\) (so \(\theta=\sigma(\pi/2-\theta_0)\)
        !! and \(J_X=(\sigma/s)\,\partial c/\partial X\)):
        !! \[ H_{XY} = \frac{\sigma c}{s}J_XJ_Y^\dagger +
        !!    \frac{\sigma}{s}\frac{\partial^2 c}{\partial X\partial Y} \]
        !! \(\partial^2c/\partial X\partial Y\) follows the same two-step
        !! recipe as simple_angle_hessian/inplane_angle_hessian: an "internal"
        !! 2x2 (v,P) chain-rule block (\(\Xi\equiv\)versor_der2,
        !! \(V\equiv\)versor_der) plus, since only \(\bP\) (not \(\bv\)) is a
        !! curved (bilinear) function of the atoms, a triple-product
        !! correction \(\chi\) that is \(\mp[\bw']_\times\)
        !! (\(\bw'=V(\bP)\hn\bv\), the current value of \(\partial c/\partial
        !! \bP\)) on the A-B/A-C/B-C blocks and 0 elsewhere (including every
        !! D-block, since \(\bP\) has no D-dependence):
        !! \[ \frac{\partial^2c}{\partial A\partial A} = \Xi(\bv,\hn\bP) +
        !!    V(\bv)V(\bP)P_A + P_A^\dagger V(\bP)V(\bv) +
        !!    P_A^\dagger\Xi(\bP,\hn\bv)P_A \]
        !! \[ \frac{\partial^2c}{\partial A\partial B} =
        !!    V(\bv)V(\bP)P_B + P_A^\dagger\Xi(\bP,\hn\bv)P_B - [\bw']_\times,
        !!    \quad
        !!    \frac{\partial^2c}{\partial A\partial D} =
        !!    -\Xi(\bv,\hn\bP) - P_A^\dagger V(\bP)V(\bv) \]
        !! \[ \frac{\partial^2c}{\partial B\partial B} =
        !!    P_B^\dagger\Xi(\bP,\hn\bv)P_B, \quad
        !!    \frac{\partial^2c}{\partial B\partial C} =
        !!    P_B^\dagger\Xi(\bP,\hn\bv)P_C - [\bw']_\times, \quad
        !!    \frac{\partial^2c}{\partial B\partial D} = -P_B^\dagger V(\bP)V(\bv) \]
        !! \[ \frac{\partial^2c}{\partial D\partial D} = \Xi(\bv,\hn\bP) \]
        !! with \(P_A=\partial\bP/\partial\bA=[\bC-\bB]_\times\),
        !! \(P_B=[\bA-\bC]_\times\), \(P_C=[\bB-\bA]_\times\), and the A-C, C-C,
        !! C-D blocks obtained from the A-B, B-B, B-D ones by \(B\to C\)
        !! (with the \([\bw']_\times\) sign flipped on A-C, matching
        !! \(\chi_{AC}=+[\bw']_\times\) vs \(\chi_{AB}=\chi_{BC}=-[\bw']_\times\)).
        use mod_utils, only: cross_product, vec_skw, versor_der, versor_der2
        use mod_constants, only: pi

        implicit none

        real(rp), intent(in), dimension(3) :: ca, cb, cc, cd
        real(rp), intent(out), dimension(3) :: J_a, J_b, J_c, J_d
        real(rp), intent(out) :: thet
        real(rp), intent(out), dimension(3,3) :: H_aa, H_ab, H_ac, H_ad, &
                                                 H_bb, H_bc, H_bd, &
                                                 H_cc, H_cd, H_dd

        real(rp), dimension(3) :: a_b, a_c, v, p, hp, hv, gv, gp, dc_da, dc_db, dc_dc, dc_dd
        real(rp), dimension(3,3) :: Pda, Pdb, Pdc, Vv, VP, XiVP, XiPV, VvVP, skw_w
        real(rp) :: costhet, thet0, s, sigma, sigfac, pref
        integer(ip) :: i

        a_b = ca - cb
        a_c = ca - cc
        v = ca - cd

        p = cross_product(a_b, a_c)
        hp = p / norm2(p)
        hv = v / norm2(v)
        costhet = dot_product(hv, hp)
        thet0 = acos(costhet)
        thet = abs(pi/2.0_rp - thet0)
        s = sqrt(1.0_rp - costhet**2)
        sigma = sign(1.0_rp, pi/2.0_rp - thet0)
        sigfac = sigma / s

        Vv = versor_der(v)
        VP = versor_der(p)
        gv = matmul(Vv, hp)
        gp = matmul(VP, hv)

        Pda = vec_skw(a_b - a_c)
        Pdb = vec_skw(a_c)
        Pdc = -vec_skw(a_b)

        dc_da = gv + matmul(transpose(Pda), gp)
        dc_db = matmul(transpose(Pdb), gp)
        dc_dc = matmul(transpose(Pdc), gp)
        dc_dd = -gv

        J_a = sigfac * dc_da
        J_b = sigfac * dc_db
        J_c = sigfac * dc_dc
        J_d = sigfac * dc_dd

        XiVP = versor_der2(v, hp)
        XiPV = versor_der2(p, hv)
        VvVP = matmul(Vv, VP)
        skw_w = vec_skw(gp)

        pref = sigma * costhet / s

        do i=1,3
            H_aa(i,:) = pref*J_a(i)*J_a
            H_ab(i,:) = pref*J_a(i)*J_b
            H_ac(i,:) = pref*J_a(i)*J_c
            H_ad(i,:) = pref*J_a(i)*J_d
            H_bb(i,:) = pref*J_b(i)*J_b
            H_bc(i,:) = pref*J_b(i)*J_c
            H_bd(i,:) = pref*J_b(i)*J_d
            H_cc(i,:) = pref*J_c(i)*J_c
            H_cd(i,:) = pref*J_c(i)*J_d
            H_dd(i,:) = pref*J_d(i)*J_d
        end do

        H_aa = H_aa + sigfac*( XiVP + matmul(VvVP,Pda) + matmul(transpose(Pda),transpose(VvVP)) &
                              + matmul(transpose(Pda),matmul(XiPV,Pda)) )
        H_ab = H_ab + sigfac*( matmul(VvVP,Pdb) + matmul(transpose(Pda),matmul(XiPV,Pdb)) - skw_w )
        H_ac = H_ac + sigfac*( matmul(VvVP,Pdc) + matmul(transpose(Pda),matmul(XiPV,Pdc)) + skw_w )
        H_ad = H_ad + sigfac*( -XiVP - matmul(transpose(Pda),transpose(VvVP)) )
        H_bb = H_bb + sigfac*( matmul(transpose(Pdb),matmul(XiPV,Pdb)) )
        H_bc = H_bc + sigfac*( matmul(transpose(Pdb),matmul(XiPV,Pdc)) - skw_w )
        H_bd = H_bd + sigfac*( -matmul(transpose(Pdb),transpose(VvVP)) )
        H_cc = H_cc + sigfac*( matmul(transpose(Pdc),matmul(XiPV,Pdc)) )
        H_cd = H_cd + sigfac*( -matmul(transpose(Pdc),transpose(VvVP)) )
        H_dd = H_dd + sigfac*( XiVP )

    end subroutine opb_angle_hessian


    subroutine pitors_angle_jacobian(ca, cb, cc, cd, ce, cf, thet, &
                                     J_a, J_b, J_c, J_d, J_e, J_f)
        use mod_utils, only: cross_product, vec_skw, versor_der

        implicit none

        real(rp), dimension(3), intent(inout) :: ca, cb, cc, cd, ce, cf
        real(rp), intent(out) :: thet
        real(rp), dimension(3), intent(out) :: J_a, J_b, J_c, J_d, J_e, J_f

        real(rp), dimension(3) :: d_b, c_b, f_a, e_a, s, p, r, t, ht, u, hu
        real(rp) :: costhet, dcosthet

        real(rp), dimension(3,3) :: skw_s, dhtdt, dhudu, &
                                    skw_c_b, skw_d_b, skw_f_a, skw_e_a, &
                                    skw_r, skw_p, dpdb, dtdb
        real(rp), dimension(3) :: dcostdp, dcostdr, dcostdt, dcostdu

        d_b = cd - cb
        c_b = cc - cb
        f_a = cf - ca
        e_a = ce - ca
        s = cb - ca
        p = cross_product(d_b, c_b)
        r = cross_product(f_a, e_a)

        t = cross_product(p, s)
        ht = t / norm2(t)
        u = cross_product(s, r)
        hu = u / norm2(u)

        costhet = dot_product(hu, ht)
        thet = acos(costhet)
        dcosthet = - 1.0 / sqrt(1.0 - costhet**2)
        
        skw_s = vec_skw(s)
        skw_r = vec_skw(r)
        skw_p = vec_skw(p)
        
        dhudu = versor_der(u)
        dhtdt = versor_der(t)

        dcostdt = matmul(dhtdt, hu)
        dcostdu = matmul(dhudu, ht)
        
        dcostdp = matmul(skw_s, dcostdt)
        dcostdr = -matmul(skw_s, dcostdu)

        skw_c_b = vec_skw(c_b)
        skw_d_b = vec_skw(d_b)
        skw_e_a = vec_skw(e_a)
        skw_f_a = vec_skw(f_a)

        dpdb = skw_d_b - skw_c_b
        dtdb = matmul(dpdb, skw_s) - skw_p

        J_b = dcosthet * (matmul(skw_r, dcostdu) + matmul(dtdb, dcostdt))
        J_c = dcosthet * (-matmul(skw_d_b, dcostdp))
        J_d = dcosthet * (matmul(skw_c_b, dcostdp))
        J_e = dcosthet * (-matmul(skw_f_a, dcostdr))
        J_f = dcosthet * (matmul(skw_e_a, dcostdr))

        J_a = -(J_b + J_c + J_d + J_e + J_f)
    end subroutine

    pure subroutine pitors_angle_hessian(ca, cb, cc, cd, ce, cf, thet, &
                                         J_a, J_b, J_c, J_d, J_e, J_f, Hblk)
        !! Hessian of the pi-torsion angle of pitors_angle_jacobian, with
        !! \(\bd_b=\bD-\bB\), \(\bc_b=\bC-\bB\), \(\bf_a=\bF-\bA\),
        !! \(\be_a=\bE-\bA\), \(\bs=\bB-\bA\), \(\bp=\bd_b\times\bc_b\),
        !! \(\br=\bf_a\times\be_a\), \(\bt=\bp\times\bs\), \(\bu=\bs\times\br\),
        !! \(c=\cos\theta=\hn\bu\cdot\hn\bt\) (unsigned, as in
        !! pitors_angle_jacobian). Since \(\theta=\arccos c\) with no extra
        !! sign flip, this uses simple_angle_hessian's exact final formula
        !! \(H_{XY}=-\frac1s(\partial^2c/\partial X\partial Y+cJ_XJ_Y^\dagger)\).
        !!
        !! \(\partial^2c/\partial X\partial Y\) is assembled from
        !! \(T_X\equiv\partial\bt/\partial X=-[\bs]_\times P_X+[\bp]_\times S_X\)
        !! and \(U_X\equiv\partial\bu/\partial X=-[\br]_\times S_X+[\bs]_\times R_X\)
        !! (\(P_X,R_X,S_X\) the Jacobians of \(\bp,\br,\bs\); \(S_A=-\mathbb1,
        !! S_B=\mathbb1\); \(P_B=[\bc_b-\bd_b]_\times,P_C=[\bd_b]_\times,
        !! P_D=-[\bc_b]_\times\); \(R_A=[\be_a-\bf_a]_\times,R_E=[\bf_a]_\times,
        !! R_F=-[\be_a]_\times\); all else 0) via the usual two-level chain
        !! rule (\(\Xi\equiv\)versor_der2, \(V\equiv\)versor_der,
        !! \(g_t=V(\bt)\hn\bu\), \(g_u=V(\bu)\hn\bt\)):
        !! \[ \frac{\partial^2c}{\partial X\partial Y} = T_X^\dagger\Xi(\bt,\hn\bu)T_Y
        !!    +T_X^\dagger V(\bt)V(\bu)U_Y+U_X^\dagger V(\bu)V(\bt)T_Y
        !!    +U_X^\dagger\Xi(\bu,\hn\bt)U_Y + (\chi_t)_{XY}+(\chi_u)_{XY} \]
        !! plus two triple-product corrections (since both \(\bt\) and \(\bu\)
        !! are themselves built from a curved piece crossed with a linear
        !! one): \((\chi_t)_{XY}=P_X^\dagger Z_Y+Z_X^\dagger P_Y+
        !! \text{tripleHess}(B,D,C;\,\bs\times g_t)\) with \(Z_A=[g_t]_\times,
        !! Z_B=-[g_t]_\times\); \((\chi_u)_{XY}=S_X^\dagger Z'_Y+
        !! Z'^\dagger_X S_Y+\text{tripleHess}(A,F,E;\,g_u\times\bs)\) with
        !! \(Z'_A=-[g_u]_\times R_A,\ Z'_E=-[g_u]_\times R_E,\
        !! Z'_F=-[g_u]_\times R_F\); tripleHess(piv,X,Y;w) is 0 on the
        !! piv-piv/X-X/Y-Y blocks and \(\mp[\bw]_\times\) on
        !! piv-X/piv-Y/X-Y (see opb_angle_hessian's \(\chi\)).
        use mod_utils, only: cross_product, vec_skw, versor_der, versor_der2

        implicit none

        real(rp), dimension(3), intent(in) :: ca, cb, cc, cd, ce, cf
        real(rp), intent(out) :: thet
        real(rp), dimension(3), intent(out) :: J_a, J_b, J_c, J_d, J_e, J_f
        real(rp), intent(out) :: Hblk(3,3,6,6)
        !! Hblk(:,:,i,j) = d^2(thet)/d(atom i) d(atom j), atom order
        !! 1=A,2=B,3=C,4=D,5=E,6=F. Hblk(:,:,j,i) = transpose(Hblk(:,:,i,j)).

        real(rp), dimension(3) :: d_b, c_b, f_a, e_a, s, p, r, t, u, ht, hu, &
                                  gt, gu, wp, wr, dcX
        real(rp), dimension(3,3) :: skw_s, skw_p, skw_r, Vt, Vu, Xit, Xiu, VtVu
        real(rp), dimension(3,3,6) :: Pmat, Rmat, Smat, Tmat, Umat, Zt, Zu
        real(rp) :: costhet, sint
        integer(ip) :: i, j, k

        d_b = cd - cb
        c_b = cc - cb
        f_a = cf - ca
        e_a = ce - ca
        s = cb - ca
        p = cross_product(d_b, c_b)
        r = cross_product(f_a, e_a)
        t = cross_product(p, s)
        ht = t / norm2(t)
        u = cross_product(s, r)
        hu = u / norm2(u)
        costhet = dot_product(hu, ht)
        thet = acos(costhet)
        sint = sqrt(1.0_rp - costhet**2)

        skw_s = vec_skw(s)
        skw_p = vec_skw(p)
        skw_r = vec_skw(r)

        ! Jacobians of p, r, s wrt the 6 atoms (1=A,2=B,3=C,4=D,5=E,6=F)
        Pmat = 0.0_rp
        Pmat(:,:,2) = vec_skw(c_b - d_b)
        Pmat(:,:,3) = vec_skw(d_b)
        Pmat(:,:,4) = -vec_skw(c_b)

        Rmat = 0.0_rp
        Rmat(:,:,1) = vec_skw(e_a - f_a)
        Rmat(:,:,5) = vec_skw(f_a)
        Rmat(:,:,6) = -vec_skw(e_a)

        Smat = 0.0_rp
        do i=1,3
            Smat(i,i,1) = -1.0_rp
            Smat(i,i,2) = 1.0_rp
        end do

        do k=1,6
            Tmat(:,:,k) = -matmul(skw_s, Pmat(:,:,k)) + matmul(skw_p, Smat(:,:,k))
            Umat(:,:,k) = -matmul(skw_r, Smat(:,:,k)) + matmul(skw_s, Rmat(:,:,k))
        end do

        Vt = versor_der(t)
        Vu = versor_der(u)
        gt = matmul(Vt, hu)
        gu = matmul(Vu, ht)

        do k=1,6
            dcX = matmul(transpose(Tmat(:,:,k)), gt) + matmul(transpose(Umat(:,:,k)), gu)
            if(k==1) J_a = -dcX/sint
            if(k==2) J_b = -dcX/sint
            if(k==3) J_c = -dcX/sint
            if(k==4) J_d = -dcX/sint
            if(k==5) J_e = -dcX/sint
            if(k==6) J_f = -dcX/sint
        end do

        Xit = versor_der2(t, hu)
        Xiu = versor_der2(u, ht)
        VtVu = matmul(Vt, Vu)

        ! chi_t: Z_A=[gt]x, Z_B=-[gt]x; triple-product correction on p
        ! (pivot B=2, X=D=4, Y=C=3), weight wp = s x gt
        Zt = 0.0_rp
        Zt(:,:,1) = vec_skw(gt)
        Zt(:,:,2) = -vec_skw(gt)
        wp = cross_product(s, gt)

        ! chi_u: Z'_A=-[gu]x*R_A, Z'_E=-[gu]x*R_E, Z'_F=-[gu]x*R_F; triple
        ! product correction on r (pivot A=1, X=F=6, Y=E=5), weight
        ! wr = gu x s
        Zu = 0.0_rp
        Zu(:,:,1) = -matmul(vec_skw(gu), Rmat(:,:,1))
        Zu(:,:,5) = -matmul(vec_skw(gu), Rmat(:,:,5))
        Zu(:,:,6) = -matmul(vec_skw(gu), Rmat(:,:,6))
        wr = cross_product(gu, s)

        do i=1,6
            do j=1,6
                Hblk(:,:,i,j) = matmul(transpose(Tmat(:,:,i)),matmul(Xit,Tmat(:,:,j))) &
                    + matmul(transpose(Tmat(:,:,i)),matmul(VtVu,Umat(:,:,j))) &
                    + matmul(transpose(Umat(:,:,i)),matmul(transpose(VtVu),Tmat(:,:,j))) &
                    + matmul(transpose(Umat(:,:,i)),matmul(Xiu,Umat(:,:,j))) &
                    + matmul(transpose(Pmat(:,:,i)),Zt(:,:,j)) + matmul(transpose(Zt(:,:,i)),Pmat(:,:,j)) &
                    + matmul(transpose(Smat(:,:,i)),Zu(:,:,j)) + matmul(transpose(Zu(:,:,i)),Smat(:,:,j))
            end do
        end do

        ! triple-product correction for p: pivot=2(B), X=4(D), Y=3(C), weight wp
        call add_triple(Hblk, 2, 4, 3, wp)
        ! triple-product correction for r: pivot=1(A), X=6(F), Y=5(E), weight wr
        call add_triple(Hblk, 1, 6, 5, wr)

        do i=1,6
            do j=1,6
                if(i==1) dcX = J_a
                if(i==2) dcX = J_b
                if(i==3) dcX = J_c
                if(i==4) dcX = J_d
                if(i==5) dcX = J_e
                if(i==6) dcX = J_f
                block
                    real(rp) :: dcY(3)
                    integer :: m
                    if(j==1) dcY = J_a
                    if(j==2) dcY = J_b
                    if(j==3) dcY = J_c
                    if(j==4) dcY = J_d
                    if(j==5) dcY = J_e
                    if(j==6) dcY = J_f
                    do m=1,3
                        Hblk(m,:,i,j) = Hblk(m,:,i,j) + costhet*dcX(m)*dcY
                    end do
                end block
                Hblk(:,:,i,j) = -Hblk(:,:,i,j)/sint
            end do
        end do

    contains
        pure subroutine add_triple(H, piv, X, Y, w)
            real(rp), intent(inout) :: H(3,3,6,6)
            integer(ip), intent(in) :: piv, X, Y
            real(rp), intent(in) :: w(3)
            real(rp) :: skw_w(3,3)
            skw_w = vec_skw(w)
            H(:,:,piv,X) = H(:,:,piv,X) - skw_w
            H(:,:,X,piv) = H(:,:,X,piv) + skw_w
            H(:,:,piv,Y) = H(:,:,piv,Y) + skw_w
            H(:,:,Y,piv) = H(:,:,Y,piv) - skw_w
            H(:,:,X,Y) = H(:,:,X,Y) - skw_w
            H(:,:,Y,X) = H(:,:,Y,X) + skw_w
        end subroutine
    end subroutine pitors_angle_hessian
end module mod_jacobian_mat
