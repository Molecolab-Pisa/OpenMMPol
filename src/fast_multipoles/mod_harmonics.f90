module mod_harmonics
    use mod_constants, only : ip, rp, pi, eps_rp
    use mod_fmm_utils, only: n_sph_harm, ntot_sph_harm
    private

    real(rp), allocatable :: vscales(:), vscales_rel(:), vcnk(:), m2l_ztranslate_coef(:,:,:)
    integer(ip) :: vscales_p = 0, vcnk_dmax = 0, m2l_pm = 0, m2l_pl = 0
    
    public :: fmm_m2m, fmm_m2l, fmm_l2l, fmm_m2p, prepare_fmmm_constants

    contains

    subroutine prepare_fmmm_constants(pm, pl)
        integer(ip), intent(in) :: pm, pl

        
        call make_vscales(max(pm,pl))
        call make_vcnk(pm)
        call make_m2l_ztranslate_coef(pm, pl)
    end subroutine
!> Compute scaling factors of real normalized spherical harmonics
!!
!! Output values of scaling factors of \f$ Y_\ell^m \f$ harmonics are filled
!! only for non-negative \f$ m \f$ since scaling factor of \f$ Y_\ell^{-m} \f$
!! is the same as scaling factor of \f$ Y_\ell^m \f$.
!!
!! @param[in] p: Maximal degree of spherical harmonics. `p` >= 0
!! @param[out] vscales: Array of scaling factors. Dimension is `(p+1)**2`
!! @param[out] vscales: Array of values 4pi/(2l+1). Dimension is `p+1`
!! @param[out] vscales_rel: Array of relative scaling factors.
!!      Dimension is `(p+1)**, vscales2`.
    subroutine make_vscales(p)
        ! Input
        integer, intent(in) :: p
        real(rp) :: tmp, twolp1, tmp2
        integer :: l, ind, m

        if(p <= vscales_p) then
            ! Job is already done
            return
        end if 

        if(allocated(vscales)) then
            ! vscales is allocated but it should be expanded.
            ! Just remove everything and restart from scratch
            deallocate(vscales)
        end if
        
        if(allocated(vscales_rel)) then
            deallocate(vscales_rel)
        end if

        vscales_p = p
        allocate(vscales((vscales_p+1)**2))
        allocate(vscales_rel((vscales_p+1)**2))

        twolp1 = 1.0_rp
        do l = 0, p
            ! m = 0
            ind = l*l + l + 1
            tmp = 4.0*pi / twolp1
            tmp2 = tmp
            tmp = sqrt(tmp)
            vscales_rel(ind) = tmp
            vscales(ind) = 1.0_rp / tmp
            twolp1 = twolp1 + 2.0_rp
            tmp = vscales(ind) * sqrt(2.0)
            ! m != 0
            do m = 1, l
                tmp = -tmp / sqrt(dble((l-m+1)*(l+m)))
                vscales(ind+m) = tmp
                vscales(ind-m) = tmp
                vscales_rel(ind+m) = tmp * tmp2
                vscales_rel(ind-m) = vscales_rel(ind+m)
            end do
        end do
    end subroutine make_vscales


    !> Compute FMM-related constants
!!
!! @param[in] dmax: Maximal degree of spherical harmonics to be evaluated.
!!      `dmax` >= 0
!! @param[in] pm: Maximal degree of the multipole expansion. `pm` >= 0.
!! @param[in] pl: Maximal degree of the local expansion. `pl` >= 0.
!! @param[out] vcnk: Array of squre roots of combinatorial factors C_n^k.
!!      Dimension is `(2*dmax+1)*(dmax+1)`.
!! @param[out] m2l_ztranslate_coef: Constants for M2L translation over OZ axis.
!!      Dimension is `(pm+1, pl+1, pl+1)`.
!! @param[out] m2l_ztranslate_coef: Constants for adjoint M2L translation over
!!      OZ axis. Dimension is `(pl+1, pl+1, pm+1)`.
subroutine make_vcnk(dmax)
    ! Inputs
    integer, intent(in) :: dmax
    ! Local variables

    integer(ip) :: i, j, indi

    if(dmax <= vcnk_dmax) then
        ! Job is already done
        return
    end if 

    if(allocated(vcnk)) then
        ! vscales is allocated but it should be expanded.
        ! Just remove everything and restart from scratch
        deallocate(vcnk)
    end if

    vcnk_dmax = dmax
    allocate(vcnk((2*dmax+1)*(dmax+1)))

    vcnk(1) = 1.0_rp
    do i = 2, 2*dmax+1
        ! Offset to the C_{i-2}^{i-2}, next item to be stored is C_{i-1}^0
        indi = (i-1) * i / 2
        ! C_{i-1}^0 = 1
        vcnk(indi+1) = 1.0_rp
        ! C_{i-1}^{i-1} = 1
        vcnk(indi+i) = 1.0_rp
        ! C_{i-1}^{j-1} = C_{i-2}^{j-1} + C_{i-2}^{j-2}
        ! Offset to C_{i-3}^{i-3} is indi-i+1
        do j = 2, i-1
            vcnk(indi+j) = vcnk(indi-i+j+1) + vcnk(indi-i+j)
        end do
    end do
    ! Get square roots of C_n^k. sqrt(1.0_rp) is 1.0_rp, so no need to update C_n^0
    ! and C_n^n
    do i = 3, 2*dmax+1
        indi = (i-1) * i / 2
        do j = 2, i-1
            vcnk(indi+j) = sqrt(vcnk(indi+j))
        end do
    end do
end subroutine make_vcnk

!> Compute FMM-related constants
!!
!! @param[in] dmax: Maximal degree of spherical harmonics to be evaluated.
!!      `dmax` >= 0
!! @param[in] pm: Maximal degree of the multipole expansion. `pm` >= 0.
!! @param[in] pl: Maximal degree of the local expansion. `pl` >= 0.
!! @param[out] vcnk: Array of squre roots of combinatorial factors C_n^k.
!!      Dimension is `(2*dmax+1)*(dmax+1)`.
!! @param[out] m2l_ztranslate_coef: Constants for M2L translation over OZ axis.
!!      Dimension is `(pm+1, pl+1, pl+1)`.
!! @param[out] m2l_ztranslate_coef: Constants for adjoint M2L translation over
!!      OZ axis. Dimension is `(pl+1, pl+1, pm+1)`.
subroutine make_m2l_ztranslate_coef(pm, pl)
    ! Inputs
    integer, intent(in) :: pm, pl

    integer(ip) :: j, k, n, indjn
    real(rp) :: tmp1
    if(pm <= m2l_pm .and. pl <= m2l_pl) then
        ! Job is already done
        return
    end if 

    if(allocated(m2l_ztranslate_coef)) then
        ! vscales is allocated but it should be expanded.
        ! Just remove everything and restart from scratch
        deallocate(m2l_ztranslate_coef)
    end if
    
    m2l_pm = pm
    m2l_pl = pl
    allocate(m2l_ztranslate_coef(m2l_pm+1, m2l_pl+1, m2l_pl+1))
    ! Check if vcnk are already populated, otherwise prepare them
    call make_vcnk(m2l_pm)

    ! Fill in m2l_ztranslate_coef
    do j = 0, pl
        do k = 0, j
            tmp1 = 1.0_rp
            do n = k, pm
                indjn = (j+n)*(j+n+1)/2 + 1
                m2l_ztranslate_coef(n-k+1, k+1, j-k+1) = &
                    & tmp1 * vcnk(indjn+j-k) * vcnk(indjn+j+k)
                !m2l_ztranslate_adj_coef(j-k+1, k+1, n-k+1) = &
                !    & m2l_ztranslate_coef(n-k+1, k+1, j-k+1)
                tmp1 = -tmp1
            end do
        end do
    end do
end subroutine

subroutine make_vfact(p, vfact)
    integer, intent(in) :: p
    ! Outputs
    real(rp), intent(out) :: vfact(2*p+1)

    integer(ip) :: i

    vfact(1) = 1
    do i=2, 2*p+1
        vfact(i) = vfact(i-1) * sqrt(real(i-1, rp))
    end do
end subroutine
    

!> Compute arrays of \f$ \cos(m \phi) \f$ and \f$ \sin(m \phi) \f$
!!
!! All values are computed recurrently from input \f$ \cos(\phi) \f$ and \f$
!! \sin(\phi) \f$ without accessing arccos or arcsin functions.
!!
!! @param[in] cphi: \f$ \cos(\phi) \f$. -1 <= `cphi` <= 1
!! @param[in] sphi: \f$ \sin(\phi) \f$. -1 <= `sphi` <= 1
!! @param[in] p: Maximal value of \f$ m \f$, for which to compute \f$ \cos(m
!!      \phi) \f$ and \f$ \sin(m\phi) \f$. `p` >= 0
!! @param[out] vcos: Array of \f$ \cos(m\phi) \f$ for \f$ m=0..p \f$. Dimension
!!      is `(p+1)`
!! @param[out] vsin: Array of \f$ \sin(m\phi) \f$ for \f$ m=0..p \f$. Dimension
!!      is `(p+1)`
    subroutine trgev(cphi, sphi, p, vcos, vsin)
        ! Inputs
        integer, intent(in) :: p
        real(rp), intent(in) :: cphi, sphi
        ! Output
        real(rp), intent(out) :: vcos(p+1), vsin(p+1)
        ! Local variables
        integer :: m
        real(rp) :: c4phi, s4phi
        !! Treat values of p from 0 to 3 differently
        select case(p)
            ! Do nothing if p < 0
            case (:-1)
                return
            ! p = 0
            case (0)
                vcos(1) = 1.0
                vsin(1) = 0.0
                return
            ! p = 1
            case (1)
                vcos(1) = 1.0
                vsin(1) = 0.0
                vcos(2) = cphi
                vsin(2) = sphi
                return
            ! p = 2
            case (2)
                vcos(1) = 1.0
                vsin(1) = 0.0
                vcos(2) = cphi
                vsin(2) = sphi
                vcos(3) = cphi**2 - sphi**2
                vsin(3) = 2 * cphi * sphi
                return
            ! p = 3
            case (3)
                vcos(1) = 1.0
                vsin(1) = 0.0
                vcos(2) = cphi
                vsin(2) = sphi
                vcos(3) = cphi**2 - sphi**2
                vsin(3) = 2 * cphi * sphi
                vcos(4) = vcos(3)*cphi - vsin(3)*sphi
                vsin(4) = vcos(3)*sphi + vsin(3)*cphi
                return
            ! p >= 4
            case default
                vcos(1) = 1.0
                vsin(1) = 0.0
                vcos(2) = cphi
                vsin(2) = sphi
                vcos(3) = cphi**2 - sphi**2
                vsin(3) = 2 * cphi * sphi
                vcos(4) = vcos(3)*cphi - vsin(3)*sphi
                vsin(4) = vcos(3)*sphi + vsin(3)*cphi
                vcos(5) = vcos(3)**2 - vsin(3)**2
                vsin(5) = 2 * vcos(3) * vsin(3)
                c4phi = vcos(5)
                s4phi = vsin(5)
        end select
        ! Define cos(m*phi) and sin(m*phi) recurrently 4 values at a time
        do m = 6, p-2, 4
            vcos(m:m+3) = vcos(m-4:m-1)*c4phi - vsin(m-4:m-1)*s4phi
            vsin(m:m+3) = vcos(m-4:m-1)*s4phi + vsin(m-4:m-1)*c4phi
        end do
        ! Work with leftover
        select case(m-p)
            case (-1)
                vcos(p-1) = vcos(p-2)*cphi - vsin(p-2)*sphi
                vsin(p-1) = vcos(p-2)*sphi + vsin(p-2)*cphi
                vcos(p) = vcos(p-2)*vcos(3) - vsin(p-2)*vsin(3)
                vsin(p) = vcos(p-2)*vsin(3) + vsin(p-2)*vcos(3)
                vcos(p+1) = vcos(p-2)*vcos(4) - vsin(p-2)*vsin(4)
                vsin(p+1) = vcos(p-2)*vsin(4) + vsin(p-2)*vcos(4)
            case (0)
                vcos(p) = vcos(p-1)*cphi - vsin(p-1)*sphi
                vsin(p) = vcos(p-1)*sphi + vsin(p-1)*cphi
                vcos(p+1) = vcos(p-1)*vcos(3) - vsin(p-1)*vsin(3)
                vsin(p+1) = vcos(p-1)*vsin(3) + vsin(p-1)*vcos(3)
            case (1)
                vcos(p+1) = vcos(p)*cphi - vsin(p)*sphi
                vsin(p+1) = vcos(p)*sphi + vsin(p)*cphi
        end select
    end subroutine trgev

!> Direct M2M translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of multipole-to-multipole
!! translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] alpha: Scalar multipler for `alpha`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size (2*(p+1))
    subroutine fmm_m2m_ztranslate_work(z, p, vscales, vcnk, alpha, &
        & src_m, beta, dst_m, work)
    implicit none
    ! Inputs
    integer, intent(in) :: p
    real(rp), intent(in) :: z, vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_m((p+1)*(p+1)), beta
    ! Output
    real(rp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Temporary workspace
    real(rp), intent(out), target :: work(2*(p+1))
    ! Local variables
    real(rp) :: r2, tmp1, tmp2, tmp3, res1, res2
    integer :: j, k, n, indj, indjn, indjk1, indjk2
    ! Pointers for temporary values of powers
    real(rp), pointer :: pow_r2(:)
    ! In case alpha is 0.0 just do a proper scaling of output
    if (abs(alpha) < eps_rp) then
        if (abs(beta) < eps_rp) then
            dst_m = 0.0
        else
            dst_m = beta * dst_m
        end if
        return
    end if
    ! Now alpha is non-0.0
    ! If harmonics have different centers
    if (abs(z) > eps_rp) then
        ! Prepare pointers
        pow_r2(1:p+1) => work(1:p+1)
        ! Get ratios r1 and r2
        r2 = z 
        pow_r2(1) = 1.0
        do j = 2, p+1
            pow_r2(j) = pow_r2(j-1) * r2
        end do
        ! Do actual M2M
        ! Overwrite output if beta is 0.0
        if (abs(beta) < eps_rp) then
            do j = 0, p
                ! Offset for dst_m
                indj = j*j + j + 1
                ! k = 0
                tmp1 = alpha * vscales(indj)
                tmp2 = tmp1
                res1 = 0.0
                ! Offset for vcnk
                indjk1 = j*(j+1)/2 + 1
                do n = 0, j
                    ! Offset for src_m
                    indjn = (j-n)**2 + (j-n) + 1
                    tmp3 = pow_r2(n+1) / &
                        & vscales(indjn) * vcnk(indjk1+n)**2
                    res1 = res1 + tmp3*src_m(indjn)
                end do
                dst_m(indj) = tmp2 * res1
                ! k != 0
                do k = 1, j
                    tmp2 = tmp1
                    res1 = 0.0
                    res2 = 0.0
                    ! Offsets for vcnk
                    indjk1 = (j-k)*(j-k+1)/2 + 1
                    indjk2 = (j+k)*(j+k+1)/2 + 1
                    do n = 0, j-k
                        ! Offset for src_m
                        indjn = (j-n)**2 + (j-n) + 1
                        tmp3 = pow_r2(n+1) / &
                            & vscales(indjn) * vcnk(indjk1+n) * &
                            & vcnk(indjk2+n)
                        res1 = res1 + tmp3*src_m(indjn+k)
                        res2 = res2 + tmp3*src_m(indjn-k)
                    end do
                    dst_m(indj+k) = tmp2 * res1
                    dst_m(indj-k) = tmp2 * res2
                end do
            end do
        ! Update output if beta is non-0.0
        else
            do j = 0, p
                ! Offset for dst_m
                indj = j*j + j + 1
                ! k = 0
                tmp1 = alpha * vscales(indj)
                tmp2 = tmp1
                res1 = 0.0
                ! Offset for vcnk
                indjk1 = j * (j+1) /2 + 1
                do n = 0, j
                    ! Offset for src_m
                    indjn = (j-n)**2 + (j-n) + 1
                    tmp3 = pow_r2(n+1) / &
                        & vscales(indjn) * vcnk(indjk1+n)**2
                    res1 = res1 + tmp3*src_m(indjn)
                end do
                dst_m(indj) = beta*dst_m(indj) + tmp2*res1
                ! k != 0
                do k = 1, j
                    tmp2 = tmp1
                    res1 = 0.0
                    res2 = 0.0
                    ! Offsets for vcnk
                    indjk1 = (j-k)*(j-k+1)/2 + 1
                    indjk2 = (j+k)*(j+k+1)/2 + 1
                    do n = 0, j-k
                        ! Offset for src_m
                        indjn = (j-n)**2 + (j-n) + 1
                        tmp3 = pow_r2(n+1) / &
                            & vscales(indjn) * vcnk(indjk1+n) * &
                            & vcnk(indjk2+n)
                        res1 = res1 + tmp3*src_m(indjn+k)
                        res2 = res2 + tmp3*src_m(indjn-k)
                    end do
                    dst_m(indj+k) = beta*dst_m(indj+k) + tmp2*res1
                    dst_m(indj-k) = beta*dst_m(indj-k) + tmp2*res2
                end do
            end do
        end if
    ! If harmonics are located at the same point
    else
        ! Overwrite output if beta is 0.0
        if (abs(beta) < eps_rp) then
            tmp1 = alpha 
            do j = 0, p
                indj = j*j + j + 1
                do k = indj-j, indj+j
                    dst_m(k) = tmp1 * src_m(k)
                end do
            end do
        ! Update output if beta is non-0.0
        else
            tmp1 = alpha
            do j = 0, p
                indj = j*j + j + 1
                do k = indj-j, indj+j
                    dst_m(k) = beta*dst_m(k) + tmp1*src_m(k)
                end do
                tmp1 = tmp1
            end do
        end if
    end if
end subroutine fmm_m2m_ztranslate_work 

!> Rotate spherical harmonics around OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha R \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics corresponding to a new cartesion system of coordinates, \f$
!! \mathrm{src} \f$ is a vector of coefficients of input spherical harmonics
!! corresponding to the standard cartesian system of coordinates and \f$
!! R \f$ is a matrix of rotation of coordinates around OZ axis on angle \f$
!! \phi \f$, presented by \f$ \cos(m \phi) \f$ and \f$ \sin(m \phi) \f$.
!!
!!
!! @param[in] p: Maximal order of spherical harmonics
!! @param[in] vcos: Vector \f$ \{ \cos(m \phi) \}_{m=0}^p \f$
!! @param[in] vsin: Vector \f$ \{ \sin(m \phi) \}_{m=0}^p \f$
!! @param[in] alpha: Scalar multiplier for `src`
!! @param[in] src: Coefficients of initial spherical harmonics
!! @param[in] beta: Scalar multipler for `dst`
!! @param[inout] dst: Coefficients of rotated spherical harmonics
subroutine fmm_sph_rotate_oz_work(p, vcos, vsin, alpha, src, beta, dst)
    ! Inputs
    integer, intent(in) :: p
    real(rp), intent(in) :: vcos(p+1), vsin(p+1), alpha, src((p+1)*(p+1)), beta
    ! Output
    real(rp), intent(inout) :: dst((p+1)*(p+1))
    ! Local variables
    integer :: l, m, ind
    real(rp) :: v1, v2, v3, v4
    ! In case alpha is 0.0 just scale output
    if (abs(alpha) < eps_rp) then
        ! Set output to 0.0 if beta is also 0.0
        if (abs(beta) < eps_rp) then
            dst = 0.0
        else
            dst = beta * dst
        end if
        ! Exit subroutine
        return
    end if
    ! Now alpha is non-0.0
    ! In case beta is 0.0 output is just overwritten without being read
    if (abs(beta) < eps_rp) then
        ! l = 0
        dst(1) = alpha*src(1)
        ! l > 0
        !!GCC$ unroll 4
        do l = 1, p
            ind = l*l + l + 1
            ! m = 0
            dst(ind) = alpha*src(ind)
            ! m != 0
            !!GCC$ unroll 4
            do m = 1, l
                v1 = src(ind+m)
                v2 = src(ind-m)
                v3 = vcos(1+m)
                v4 = vsin(1+m)
                ! m > 0
                dst(ind+m) = alpha * (v1*v3-v2*v4)
                ! m < 0
                dst(ind-m) = alpha * (v1*v4+v2*v3)
            end do
        end do
    else
        ! l = 0
        dst(1) = beta*dst(1) + alpha*src(1)
        ! l > 0
        !!GCC$ unroll 4
        do l = 1, p
            ind = l*l + l + 1
            ! m = 0
            dst(ind) = beta*dst(ind) + alpha*src(ind)
            ! m != 0
            !!GCC$ unroll 4
            do m = 1, l
                v1 = src(ind+m)
                v2 = src(ind-m)
                v3 = vcos(1+m)
                v4 = vsin(1+m)
                ! m > 0
                dst(ind+m) = beta*dst(ind+m) + alpha*(v1*v3-v2*v4)
                ! m < 0
                dst(ind-m) = beta*dst(ind-m) + alpha*(v1*v4+v2*v3)
            end do
        end do
    end if
end subroutine fmm_sph_rotate_oz_work

!> Convert input cartesian coordinate into spherical coordinate
!!
!! Output coordinate \f$ (\rho, \theta, \phi) \f$ is presented by \f$ (\rho,
!! \cos \theta, \sin \theta, \cos \phi, \sin\phi) \f$.
!!
!! @param[in] x: Cartesian coordinate
!! @param[out] rho: \f$ \rho \f$
!! @param[out] ctheta: \f$ \cos \theta \f$
!! @param[out] stheta: \f$ \sin \theta \f$
!! @param[out] cphi: \f$ \cos \phi \f$
!! @param[out] sphi: \f$ \sin \phi \f$
    subroutine carttosph(x, rho, ctheta, stheta, cphi, sphi)
        ! Input
        real(rp), intent(in) :: x(3)
        ! Output
        real(rp), intent(out) :: rho, ctheta, stheta, cphi, sphi
        ! Local variables
        real(rp) :: max12, ssq12
        ! Check x(1:2) = 0
        if ((abs(x(1)) < eps_rp) .and. (abs(x(2)) < eps_rp)) then
            rho = abs(x(3))
            ctheta = sign(1.0_rp, x(3))
            stheta = 0.0
            cphi = 1.0
            sphi = 0.0
            return
        end if
        ! In other situations use sum-of-scaled-squares technique
        ! Get norm of x(1:2) and cphi with sphi outputs
        if (abs(x(2)) .gt. abs(x(1))) then
            max12 = abs(x(2))
            ssq12 = 1.0 + (x(1)/x(2))**2
        else
            max12 = abs(x(1))
            ssq12 = 1.0 + (x(2)/x(1))**2
        end if
        stheta = max12 * sqrt(ssq12)
        cphi = x(1) / stheta
        sphi = x(2) / stheta
        ! Then compute rho, ctheta and stheta outputs
        if (abs(x(3)) .gt. max12) then
            rho = 1.0 + ssq12*(max12/x(3))**2
            rho = abs(x(3)) * sqrt(rho)
            stheta = stheta / rho
            ctheta = x(3) / rho
        else
            rho = ssq12 + (x(3)/max12)**2
            rho = max12 * sqrt(rho)
            stheta = stheta / rho
            ctheta = x(3) / rho
        end if
    end subroutine carttosph

!> Transform spherical harmonics in the OXZ plane
!!
!! This function implements @ref fmm_sph_rotate_oxz with predefined values
!! of parameters \p alpha=1.0 and \p beta=0.0.
!! 
!! @param[in] p: maximum order of spherical harmonics
!! @param[in] r1xz: 2D transformation matrix in the OXZ plane
!! @param[in] alpha: Scalar multiplier for `src`
!! @param[in] src: Coefficients of initial spherical harmonics
!! @param[in] beta: Scalar multipler for `dst`
!! @param[out] dst: coefficients of rotated spherical harmonics
!! @param[out] work: Temporary workspace of a size (2*(2*p+1)*(2*p+3))
    subroutine fmm_sph_rotate_oxz_work(p, ctheta, stheta, alpha, src, beta, dst, &
        & work)
    ! Inputs
    integer, intent(in) :: p
    real(rp), intent(in) :: ctheta, stheta, alpha, src((p+1)**2), beta
    ! Output
    real(rp), intent(out) :: dst((p+1)*(p+1))
    ! Temporary workspace
    real(rp), intent(out), target :: work(4*p*p+13*p+4)
    ! Local variables
    real(rp) :: u, v, w, fl, fl2, tmp1, tmp2, vu(2), vv(2), vw(2), &
        ctheta2, stheta2, cstheta
    integer :: l, m, n, ind
    ! Pointers for a workspace
    real(rp), pointer :: r(:, :, :), r_prev(:, :, :), scal_uvw_m(:), &
        & scal_u_n(:), scal_v_n(:), scal_w_n(:), r_swap(:, :, :), vsqr(:)
    !! Spherical harmonics Y_l^m with negative m transform into harmonics Y_l^m
    !! with the same l and negative m, while harmonics Y_l^m with non-negative
    !! m transform into harmonics Y_l^m with the same l and non-negative m.
    !! Transformations for non-negative m will be stored in r(1, :, :) and for
    !! negative m will be in r(2, :, :)
    ! In case alpha is 0.0 just scale output
    if (abs(alpha) < eps_rp) then
        ! Set output to 0.0 if beta is also 0.0
        if (abs(beta) < eps_rp) then
            dst = 0.0
        else
            dst = beta * dst
        end if
        ! Exit subroutine
        return
    end if
    ! Now alpha is non-0.0
    ! In case beta is 0.0 output is just overwritten without being read
    if (abs(beta) < eps_rp) then
        ! Compute rotations/reflections
        ! l = 0
        dst(1) = alpha * src(1)
        if (p .eq. 0) then
            return
        end if
        ! l = 1
        dst(2) = alpha * src(2)
        dst(3) = alpha * (src(3)*ctheta - src(4)*stheta)
        dst(4) = alpha * (src(3)*stheta + src(4)*ctheta)
        if (p .eq. 1) then
            return
        end if
        ! Set pointers
        l = 2 * (p+1) * (p+1)
        r(1:2, 0:p, 0:p) => work(1:l)
        m = 2 * l
        r_prev(1:2, 0:p, 0:p) => work(l+1:m)
        l = m + p + 1
        scal_uvw_m(0:p) => work(m+1:l)
        m = l + p
        scal_u_n(0:p-1) => work(l+1:m)
        l = m + p + 1
        scal_v_n(0:p) => work(m+1:l)
        m = l + p - 2
        scal_w_n(1:p-2) => work(l+1:m)
        l = m + p
        vsqr(1:p) => work(m+1:l)
        ! l = 2, m >= 0
        ctheta2 = ctheta * ctheta
        cstheta = ctheta * stheta
        stheta2 = stheta * stheta
        r(1, 2, 2) = (ctheta2 + 1.0) / 2.0
        r(1, 1, 2) = cstheta
        r(1, 0, 2) = sqrt(3.0) / 2.0 * stheta2
        dst(9) = alpha * (src(9)*r(1, 2, 2) + src(8)*r(1, 1, 2) + &
            & src(7)*r(1, 0, 2))
        r(1, 2, 1) = -cstheta
        r(1, 1, 1) = ctheta2 - stheta2
        r(1, 0, 1) = sqrt(3.0) * cstheta
        dst(8) = alpha * (src(9)*r(1, 2, 1) + src(8)*r(1, 1, 1) + &
            & src(7)*r(1, 0, 1))
        r(1, 2, 0) = sqrt(3.0) / 2.0 * stheta2
        r(1, 1, 0) = -sqrt(3.0) * cstheta
        r(1, 0, 0) = (3.0*ctheta2-1.0) / 2.0
        dst(7) = alpha * (src(9)*r(1, 2, 0) + src(8)*r(1, 1, 0) + &
            & src(7)*r(1, 0, 0))
        ! l = 2,  m < 0
        r(2, 1, 1) = ctheta
        r(2, 2, 1) = -stheta
        dst(6) = alpha * (src(6)*r(2, 1, 1) + src(5)*r(2, 2, 1))
        r(2, 1, 2) = stheta
        r(2, 2, 2) = ctheta
        dst(5) = alpha * (src(6)*r(2, 1, 2) + src(5)*r(2, 2, 2))
        ! l > 2
        vsqr(1) = 1.0
        vsqr(2) = 4.0
        do l = 3, p
            ! Swap previous and current rotation matrices
            r_swap => r_prev
            r_prev => r
            r => r_swap
            ! Prepare scalar factors
            fl = dble(l)
            fl2 = fl * fl
            vsqr(l) = fl2
            scal_uvw_m(0) = 1.0 / fl
            !!GCC$ unroll 4
            do m = 1, l-1
                u = sqrt(fl2 - vsqr(m))
                scal_uvw_m(m) = 1.0 / u
            end do
            u = 2.0 * dble(l)
            u = sqrt(dble(u*(u-1.0)))
            scal_uvw_m(l) = 1.0 / u
            scal_u_n(0) = dble(l)
            scal_v_n(0) = -sqrt(dble(l*(l-1))) / sqrt(2.0)
            !!GCC$ unroll 4
            do n = 1, l-2
                u = sqrt(fl2-vsqr(n))
                scal_u_n(n) = u
            end do
            !!GCC$ unroll 4
            do n = 1, l-2
                v = dble(l+n)
                v = sqrt(v*v-v) / 2.0
                scal_v_n(n) = v
                w = dble(l-n)
                w = -sqrt(w*w-w) / 2.0
                scal_w_n(n) = w
            end do
            u = sqrt(dble(2*l-1))
            scal_u_n(l-1) = u
            v = sqrt(dble((2*l-1)*(2*l-2))) / 2.0
            scal_v_n(l-1) = v
            v = sqrt(dble(2*l*(2*l-1))) / 2.0
            scal_v_n(l) = v
            ind = l*l + l + 1
            ! m = l, n = l and m = -l, n = - l
            vv = ctheta*r_prev(:, l-1, l-1) + r_prev(2:1:-1, l-1, l-1)
            r(:, l, l) = vv * scal_v_n(l) * scal_uvw_m(l)
            tmp1 = src(ind+l) * r(1, l, l)
            tmp2 = src(ind-l) * r(2, l, l)
            ! m = l, n = l-1 and m = -l, n = 1-l
            vu = stheta * r_prev(:, l-1, l-1)
            vv = ctheta*r_prev(:, l-2, l-1) + r_prev(2:1:-1, l-2, l-1)
            r(:, l-1, l) = (vu*scal_u_n(l-1)+vv*scal_v_n(l-1)) * scal_uvw_m(l)
            tmp1 = tmp1 + src(ind+l-1)*r(1, l-1, l)
            tmp2 = tmp2 + src(ind-l+1)*r(2, l-1, l)
            ! m = l, n = 1 and m = -l, n = -1
            vu = stheta * r_prev(:, 1, l-1)
            vv(1) = ctheta * r_prev(1, 0, l-1)
            vv(2) = r_prev(1, 0, l-1)
            vw = ctheta*r_prev(:, 2, l-1) - r_prev(2:1:-1, 2, l-1)
            r(:, 1, l) = vu*scal_u_n(1) + vw*scal_w_n(1) + sqrt(2.0)*scal_v_n(1)*vv
            r(:, 1, l) = r(:, 1, l) * scal_uvw_m(l)
            tmp1 = tmp1 + src(ind+1)*r(1, 1, l)
            tmp2 = tmp2 + src(ind-1)*r(2, 1, l)
            ! m = l, n = 0
            u = stheta * r_prev(1, 0, l-1)
            v = ctheta*r_prev(1, 1, l-1) - r_prev(2, 1, l-1)
            r(1, 0, l) = (u*scal_u_n(0) + v*scal_v_n(0)) * scal_uvw_m(l)
            tmp1 = tmp1 + src(ind)*r(1, 0, l)
            ! m = l, n = 2..l-2 and m = -l, n = 2-l..-2
            !!GCC$ unroll 4
            do n = 2, l-2
                vu = stheta * r_prev(:, n, l-1)
                vv = ctheta*r_prev(:, n-1, l-1) + r_prev(2:1:-1, n-1, l-1)
                vw = ctheta*r_prev(:, n+1, l-1) - r_prev(2:1:-1, n+1, l-1)
                vu = vu*scal_u_n(n) + vv*scal_v_n(n) + vw*scal_w_n(n)
                r(:, n, l) = vu * scal_uvw_m(l)
                tmp1 = tmp1 + src(ind+n)*r(1, n, l)
                tmp2 = tmp2 + src(ind-n)*r(2, n, l)
            end do
            dst(ind+l) = alpha * tmp1
            dst(ind-l) = alpha * tmp2
            ! Now deal with m = 0
            ! n = l and n = -l
            v = -stheta * r_prev(1, l-1, 0)
            u = scal_v_n(l) * scal_uvw_m(0)
            r(1, l, 0) = v * u
            tmp1 = src(ind+l) * r(1, l, 0)
            ! n = l-1
            u = ctheta * r_prev(1, l-1, 0)
            v = -stheta * r_prev(1, l-2, 0)
            w = u*scal_u_n(l-1) + v*scal_v_n(l-1)
            r(1, l-1, 0) = w * scal_uvw_m(0)
            tmp1 = tmp1 + src(ind+l-1)*r(1, l-1, 0)
            ! n = 0
            u = ctheta * r_prev(1, 0, 0)
            v = -stheta * r_prev(1, 1, 0)
            w = u*scal_u_n(0) + v*scal_v_n(0)
            r(1, 0, 0) = w * scal_uvw_m(0)
            tmp1 = tmp1 + src(ind)*r(1, 0, 0)
            ! n = 1
            v = sqrt(2.0)*scal_v_n(1)*r_prev(1, 0, 0) + &
                & scal_w_n(1)*r_prev(1, 2, 0)
            u = ctheta * r_prev(1, 1, 0)
            w = scal_u_n(1)*u - stheta*v
            r(1, 1, 0) = w * scal_uvw_m(0)
            tmp1 = tmp1 + src(ind+1)*r(1, 1, 0)
            ! n = 2..l-2
            !!GCC$ unroll 4
            do n = 2, l-2
                v = scal_v_n(n)*r_prev(1, n-1, 0) + &
                    & scal_w_n(n)*r_prev(1, n+1, 0)
                u = ctheta * r_prev(1, n, 0)
                w = scal_u_n(n)*u - stheta*v
                r(1, n, 0) = w * scal_uvw_m(0)
                tmp1 = tmp1 + src(ind+n)*r(1, n, 0)
            end do
            dst(ind) = alpha * tmp1
            ! Now deal with m=1..l-1 and m=1-l..-1
            !!GCC$ unroll 4
            do m = 1, l-1
                ! n = l and n = -l
                vv = -stheta * r_prev(:, l-1, m)
                u = scal_v_n(l) * scal_uvw_m(m)
                r(:, l, m) = vv * u
                tmp1 = src(ind+l) * r(1, l, m)
                tmp2 = src(ind-l) * r(2, l, m)
                ! n = l-1 and n = 1-l
                vu = ctheta * r_prev(:, l-1, m)
                vv = -stheta * r_prev(:, l-2, m)
                vw = vu*scal_u_n(l-1) + vv*scal_v_n(l-1)
                r(:, l-1, m) = vw * scal_uvw_m(m)
                tmp1 = tmp1 + src(ind+l-1)*r(1, l-1, m)
                tmp2 = tmp2 + src(ind-l+1)*r(2, l-1, m)
                ! n = 0
                u = ctheta * r_prev(1, 0, m)
                v = -stheta * r_prev(1, 1, m)
                w = u*scal_u_n(0) + v*scal_v_n(0)
                r(1, 0, m) = w * scal_uvw_m(m)
                tmp1 = tmp1 + src(ind)*r(1, 0, m)
                ! n = 1
                v = sqrt(2.0)*scal_v_n(1)*r_prev(1, 0, m) + &
                    & scal_w_n(1)*r_prev(1, 2, m)
                u = ctheta * r_prev(1, 1, m)
                w = scal_u_n(1)*u - stheta*v
                r(1, 1, m) = w * scal_uvw_m(m)
                tmp1 = tmp1 + src(ind+1)*r(1, 1, m)
                ! n = -1
                u = ctheta * r_prev(2, 1, m)
                w = -stheta * r_prev(2, 2, m)
                v = u*scal_u_n(1) + w*scal_w_n(1)
                r(2, 1, m) = v * scal_uvw_m(m)
                tmp2 = tmp2 + src(ind-1)*r(2, 1, m)
                ! n = 2..l-2 and n = 2-l..-2
                !!GCC$ unroll 4
                do n = 2, l-2
                    vv = scal_v_n(n)*r_prev(:, n-1, m) + &
                        & scal_w_n(n)*r_prev(:, n+1, m)
                    vu = ctheta * r_prev(:, n, m)
                    vw = scal_u_n(n)*vu - stheta*vv
                    r(:, n, m) = vw * scal_uvw_m(m)
                    tmp1 = tmp1 + src(ind+n)*r(1, n, m)
                    tmp2 = tmp2 + src(ind-n)*r(2, n, m)
                end do
                dst(ind+m) = alpha * tmp1
                dst(ind-m) = alpha * tmp2
            end do
        end do
    else
        stop "Not Implemented"
    end if
end subroutine fmm_sph_rotate_oxz_work

!> Rotate spherical harmonics around OZ axis in an opposite direction
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha R \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics corresponding to a new cartesion system of coordinates, \f$
!! \mathrm{src} \f$ is a vector of coefficients of input spherical harmonics
!! corresponding to the standard cartesian system of coordinates and \f$
!! R \f$ is a matrix of rotation of coordinates around OZ axis on angle \f$
!! \phi \f$, presented by \f$ \cos(m \phi) \f$ and \f$ \sin(m \phi) \f$.
!!
!!
!! @param[in] p: Maximal order of spherical harmonics
!! @param[in] vcos: Vector \f$ \{ \cos(m \phi) \}_{m=0}^p \f$
!! @param[in] vsin: Vector \f$ \{ \sin(m \phi) \}_{m=0}^p \f$
!! @param[in] alpha: Scalar multiplier for `src`
!! @param[in] src: Coefficients of initial spherical harmonics
!! @param[in] beta: Scalar multipler for `dst`
!! @param[inout] dst: Coefficients of rotated spherical harmonics
    subroutine fmm_sph_rotate_oz_adj_work(p, vcos, vsin, alpha, src, beta, dst)
        ! Inputs
        integer, intent(in) :: p
        real(rp), intent(in) :: vcos(p+1), vsin(p+1), alpha, src((p+1)*(p+1)), beta
        ! Output
        real(rp), intent(inout) :: dst((p+1)*(p+1))
        ! Local variables
        integer :: l, m, ind
        real(rp) :: v1, v2, v3, v4
        ! In case alpha is 0.0 just scale output
        if (abs(alpha) < eps_rp) then
            ! Set output to 0.0 if beta is also 0.0
            if (abs(beta) < eps_rp) then
                dst = 0.0
            else
                dst = beta * dst
            end if
            ! Exit subroutine
            return
        end if
        ! Now alpha is non-0.0
        ! In case beta is 0.0 output is just overwritten without being read
        if (abs(beta) < eps_rp) then
            ! l = 0
            dst(1) = alpha*src(1)
            ! l > 0
            do l = 1, p
                ind = l*l + l + 1
                ! m = 0
                dst(ind) = alpha*src(ind)
                ! m != 0
                do m = 1, l
                    v1 = src(ind+m)
                    v2 = src(ind-m)
                    v3 = vcos(1+m)
                    v4 = vsin(1+m)
                    ! m > 0
                    dst(ind+m) = alpha * (v1*v3+v2*v4)
                    ! m < 0
                    dst(ind-m) = alpha * (v2*v3-v1*v4)
                end do
            end do
        else
            ! l = 0
            dst(1) = beta*dst(1) + alpha*src(1)
            ! l > 0
            do l = 1, p
                ind = l*l + l + 1
                ! m = 0
                dst(ind) = beta*dst(ind) + alpha*src(ind)
                ! m != 0
                do m = 1, l
                    v1 = src(ind+m)
                    v2 = src(ind-m)
                    v3 = vcos(1+m)
                    v4 = vsin(1+m)
                    ! m > 0
                    dst(ind+m) = beta*dst(ind+m) + alpha*(v1*v3+v2*v4)
                    ! m < 0
                    dst(ind-m) = beta*dst(ind-m) + alpha*(v2*v3-v1*v4)
                end do
            end do
        end if
    end subroutine fmm_sph_rotate_oz_adj_work

    !> Direct M2L translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_M \f$ is a matrix of multipole-to-local
!! translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old (multipole) harmonics
!! @param[in] dst_r: Radius of new (local) harmonics
!! @parma[in] pm: Maximal degree of multipole spherical harmonics
!! @parma[in] pl: Maximal degree of local spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] m2l_ztranslate_coef:
!! @param[in] alpha: Scalar multipler for `src_m`
!! @param[in] src_m: Expansion in old (multipole) harmonics
!! @param[in] beta: Scalar multipler for `dst_l`
!! @param[inout] dst_l: Expansion in new (local) harmonics
!! @param[out] work: Temporary workspace of a size (pm+2)*(pm+1)
subroutine fmm_m2l_ztranslate_work(z, pm, pl, vscales, &
    & m2l_ztranslate_coef, alpha, src_m, beta, dst_l, work)
! Inputs
integer, intent(in) :: pm, pl
real(rp), intent(in) :: z, vscales((pm+pl+1)*(pm+pl+1)), &
    & m2l_ztranslate_coef(pm+1, pl+1, pl+1), alpha, src_m((pm+1)*(pm+1)), &
    & beta
! Output
real(rp), intent(inout) :: dst_l((pl+1)*(pl+1))
! Temporary workspace
real(rp), intent(out), target :: work((pm+2)*(pm+1))
! Local variables
real(rp) :: tmp1, r1, r2, res1, res2, pow_r2
integer :: j, k, n, indj, indk1, indk2
! Pointers for temporary values of powers
real(rp), pointer :: src_m2(:), pow_r1(:)
! In case alpha is 0.0_rp just do a proper scaling of output
if (abs(alpha) < eps_rp) then
    if (abs(beta) < eps_rp) then
        dst_l = 0.0_rp
    else
        dst_l = beta * dst_l
    end if
    return
end if
! Now alpha is non-0.0_rp
! z cannot be 0.0_rp, as input sphere (multipole) must not intersect with
! output sphere (local)
if (abs(z) < eps_rp) then
    return
end if
! Prepare pointers
n = (pm+1) ** 2
src_m2(1:n) => work(1:n)
pow_r1(1:pm+1) => work(n+1:n+pm+1)
! Get powers of r1 and r2
r1 = 1.0_rp / z
r2 = 1.0_rp / z
! This abs(r1) makes it possible to work with negative z to avoid
! unnecessary rotation to positive z
tmp1 = abs(r1)
do j = 0, pm
    indj = j*j + j + 1
    pow_r1(j+1) = tmp1 / vscales(indj)
    tmp1 = tmp1 * r1
end do
pow_r2 = 1.0_rp
! Reorder source harmonics from (degree, order) to (order, degree)
! 0.0_rp order k=0 at first
do j = 0, pm
    indj = j*j + j + 1
    src_m2(j+1) = pow_r1(j+1) * src_m(indj)
end do
! Non-0.0_rp orders next, a positive k followed by a negative -k
indk1 = pm + 2
do k = 1, pm
    n = pm - k + 1
    indk2 = indk1 + n
    !!GCC$ unroll 4
    do j = k, pm
        indj = j*j + j + 1
        src_m2(indk1+j-k) = pow_r1(j+1) * src_m(indj+k)
        src_m2(indk2+j-k) = pow_r1(j+1) * src_m(indj-k)
    end do
    indk1 = indk2 + n
end do
! Do actual M2L
! Overwrite output if beta is 0.0_rp
if (abs(beta) < eps_rp) then
    do j = 0, pl
        ! Offset for dst_l
        indj = j*j + j + 1
        ! k = 0
        tmp1 = alpha * vscales(indj) * pow_r2
        pow_r2 = pow_r2 * r2
        res1 = 0.0_rp
        !!GCC$ unroll 4
        do n = 0, pm
            res1 = res1 + m2l_ztranslate_coef(n+1, 1, j+1)*src_m2(n+1)
        end do
        dst_l(indj) = tmp1 * res1
        ! k != 0
        !!GCC$ unroll 4
        do k = 1, j
            ! Offsets for src_m2
            indk1 = pm + 2 + (2*pm-k+2)*(k-1)
            indk2 = indk1 + pm - k + 1
            res1 = 0.0_rp
            res2 = 0.0_rp
            !!GCC$ unroll 4
            do n = k, pm
                res1 = res1 + &
                    & m2l_ztranslate_coef(n-k+1, k+1, j-k+1)* &
                    & src_m2(indk1+n-k)
                res2 = res2 + &
                    & m2l_ztranslate_coef(n-k+1, k+1, j-k+1)* &
                    & src_m2(indk2+n-k)
            end do
            dst_l(indj+k) = tmp1 * res1
            dst_l(indj-k) = tmp1 * res2
        end do
    end do
else
    do j = 0, pl
        ! Offset for dst_l
        indj = j*j + j + 1
        ! k = 0
        tmp1 = alpha * vscales(indj) * pow_r2
        pow_r2 = pow_r2 * r2
        res1 = 0.0_rp
        do n = 0, pm
            res1 = res1 + m2l_ztranslate_coef(n+1, 1, j+1)*src_m2(n+1)
        end do
        dst_l(indj) = beta*dst_l(indj) + tmp1*res1
        ! k != 0
        do k = 1, j
            ! Offsets for src_m2
            indk1 = pm + 2 + (2*pm-k+2)*(k-1)
            indk2 = indk1 + pm - k + 1
            res1 = 0.0_rp
            res2 = 0.0_rp
            do n = k, pm
                res1 = res1 + &
                    & m2l_ztranslate_coef(n-k+1, k+1, j-k+1)* &
                    & src_m2(indk1+n-k)
                res2 = res2 + &
                    & m2l_ztranslate_coef(n-k+1, k+1, j-k+1)* &
                    & src_m2(indk2+n-k)
            end do
            dst_l(indj+k) = beta*dst_l(indj+k) + tmp1*res1
            dst_l(indj-k) = beta*dst_l(indj-k) + tmp1*res2
        end do
    end do
end if
end subroutine fmm_m2l_ztranslate_work

!> Direct L2L translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L \f$ is a matrix of local-to-local
!! translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multipler for `src_l`
!! @param[in] src_l: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size (2*(p+1))
subroutine fmm_l2l_ztranslate_work(z, p, vscales, vfact, alpha, &
    & src_l, beta, dst_l, work)
    ! Inputs
    integer, intent(in) :: p
    real(rp), intent(in) :: z, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_l((p+1)*(p+1)), beta
    ! Output
    real(rp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Temporary workspace
    real(rp), intent(out), target :: work(2*(p+1))
    ! Local variables
    real(rp) :: r1, r2, tmp1, tmp2
    integer :: j, k, n, indj, indn
    ! Pointers for temporary values of powers
    real(rp), pointer :: pow_r1(:), pow_r2(:)
    ! Init output
    if (abs(beta) < eps_rp) then
        dst_l = 0.0_rp
    else
        dst_l = beta * dst_l
    end if
    ! If harmonics have different centers
    if (abs(z) > eps_rp) then
        ! Prepare pointers
        pow_r1(1:p+1) => work(1:p+1)
        pow_r2(1:p+1) => work(p+2:2*(p+1))
        ! Get powers of r1 and r2
        r1 = z / 1.0_rp
        r2 = 1.0_rp / z
        pow_r1(1) = 1
        pow_r2(1) = 1
        do j = 2, p+1
            pow_r1(j) = pow_r1(j-1) * r1
            pow_r2(j) = pow_r2(j-1) * r2
        end do
        ! Do actual L2L
        do j = 0, p
            indj = j*j + j + 1
            do k = 0, j
                tmp1 = alpha * pow_r2(j+1) / vfact(j-k+1) / vfact(j+k+1) * &
                    & vscales(indj)
                do n = j, p
                    indn = n*n + n + 1
                    tmp2 = tmp1 * pow_r1(n+1) / vscales(indn) * &
                        & vfact(n-k+1) * vfact(n+k+1) / vfact(n-j+1) / &
                        & vfact(n-j+1)
                    if (mod(n+j, 2) .eq. 1) then
                        tmp2 = -tmp2
                    end if
                    if (k .eq. 0) then
                        dst_l(indj) = dst_l(indj) + tmp2*src_l(indn)
                    else
                        dst_l(indj+k) = dst_l(indj+k) + tmp2*src_l(indn+k)
                        dst_l(indj-k) = dst_l(indj-k) + tmp2*src_l(indn-k)
                    end if
                end do
            end do
        end do
    ! If harmonics are located at the same point
    else
        tmp1 = alpha
        do j = 0, p
            indj = j*j + j + 1
            do k = indj-j, indj+j
                dst_l(k) = dst_l(k) + src_l(k)*tmp1
            end do
            tmp1 = tmp1
        end do
    end if
end subroutine fmm_l2l_ztranslate_work

!> Direct M2M translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of a multipole-to-multipole
!! translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vcnk: Square roots of combinatorial numbers C_n^k
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size 6*p*p+19*p+8
subroutine fmm_m2m_rotation_work(c, p, vscales, vcnk, alpha, &
    & src_m, beta, dst_m, work)
    ! Inputs
    integer, intent(in) :: p
    real(rp), intent(in) :: c(3), vscales((p+1)*(p+1)), &
        & vcnk((2*p+1)*(p+1)), alpha, src_m((p+1)*(p+1)), beta
    ! Output
    real(rp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Temporary workspace
    real(rp), intent(out), target :: work(6*p*p + 19*p + 8)
    ! Local variables
    real(rp) :: rho, ctheta, stheta, cphi, sphi
    integer :: m, n
    ! Pointers for temporary values of harmonics
    real(rp), pointer :: tmp_m(:), tmp_m2(:), vcos(:), vsin(:)
    ! Convert Cartesian coordinates into spherical
    call carttosph(c, rho, ctheta, stheta, cphi, sphi)
    ! If no need for rotations, just do translation along z
    if (abs(stheta) < eps_rp) then
        ! Workspace here is 2*(p+1)
        call fmm_m2m_ztranslate_work(c(3), p, vscales, vcnk, &
            & alpha, src_m, beta, dst_m, work)
        return
    end if
    ! Prepare pointers
    m = (p+1)**2
    n = 4*m + 5*p ! 4*p*p + 13*p + 4
    tmp_m(1:m) => work(n+1:n+m) ! 5*p*p + 15*p + 5
    n = n + m
    tmp_m2(1:m) => work(n+1:n+m) ! 6*p*p + 17*p + 6
    n = n + m
    m = p + 1
    vcos => work(n+1:n+m) ! 6*p*p + 18*p + 7
    n = n + m
    vsin => work(n+1:n+m) ! 6*p*p + 19*p + 8
    ! Compute arrays of cos and sin that are needed for rotations of harmonics
    call trgev(cphi, sphi, p, vcos, vsin)
    ! Rotate around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_adj_work(p, vcos, vsin, alpha, src_m, 0.0_rp, tmp_m)
    ! Perform rotation in the OXZ plane, work size is 4*p*p+13*p+4
    call fmm_sph_rotate_oxz_work(p, ctheta, -stheta, 1.0_rp, tmp_m, 0.0_rp, &
        & tmp_m2, work)
    ! OZ translation, workspace here is 2*(p+1)
    call fmm_m2m_ztranslate_work(rho, p, vscales, vcnk, 1.0_rp, &
        & tmp_m2, 0.0_rp, tmp_m, work)
    ! Backward rotation in the OXZ plane, work size is 4*p*p+13*p+4
    call fmm_sph_rotate_oxz_work(p, ctheta, stheta, 1.0_rp, tmp_m, 0.0_rp, tmp_m2, &
        & work)
    ! Backward rotation around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_work(p, vcos, vsin, 1.0_rp, tmp_m2, beta, dst_m)
end subroutine fmm_m2m_rotation_work

!> Direct M2L translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_M \f$ is a matrix of a multipole-to-local
!! translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] pl: Maximal degree of local spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] m2l_ztranslate_coef:
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size 6*p*p+19*p+8 where p is a
!!      maximum of pm and pl
subroutine fmm_m2l_rotation_work(c, pm, pl, vscales, &
    & m2l_ztranslate_coef, alpha, src_m, beta, dst_l, work)
    ! Inputs
    integer, intent(in) :: pm, pl
    real(rp), intent(in) :: c(3), vscales((pm+pl+1)**2), &
        & m2l_ztranslate_coef(pm+1, pl+1, pl+1), alpha, src_m((pm+1)*(pm+1)), &
        & beta
    ! Output
    real(rp), intent(inout) :: dst_l((pl+1)*(pl+1))
    ! Temporary workspace
    real(rp), intent(out), target :: &
        & work(6*max(pm, pl)**2 + 19*max(pm, pl) + 8)
    ! Local variables
    real(rp) :: rho, ctheta, stheta, cphi, sphi
    integer :: m, n, p
    ! Pointers for temporary values of harmonics
    real(rp), pointer :: tmp_ml(:), tmp_ml2(:), vcos(:), vsin(:)
    !call time_push
    !call time_push
    ! Covert Cartesian coordinates into spherical
    call carttosph(c, rho, ctheta, stheta, cphi, sphi)
    !call time_pull('carttosph')
    ! If no need for rotations, just do translation along z
    if (abs(stheta) < eps_rp) then
        ! Workspace here is (pm+2)*(pm+1)
        !call time_push
        call fmm_m2l_ztranslate_work(c(3), pm, pl, vscales, &
            & m2l_ztranslate_coef, alpha, src_m, beta, dst_l, work)
        !call time_pull("fmm_m2l_ztranslate_work")
        return
    end if
    !call time_push
    ! Prepare pointers
    p = max(pm, pl)
    m = (p+1)**2
    n = 4*m + 5*p ! 4*p*p + 13*p + 4
    tmp_ml(1:m) => work(n+1:n+m) ! 5*p*p + 15*p + 5
    n = n + m
    tmp_ml2(1:m) => work(n+1:n+m) ! 6*p*p + 17*p + 6
    n = n + m
    m = p + 1
    vcos => work(n+1:n+m) ! 6*p*p + 18*p + 7
    n = n + m
    vsin => work(n+1:n+m) ! 6*p*p + 19*p + 8
    !call time_pull('Pointer preparation')
    !call time_push
    ! Compute arrays of cos and sin that are needed for rotations of harmonics
    call trgev(cphi, sphi, p, vcos, vsin)
    !call time_pull('trgev')
    ! Rotate around OZ axis (work array might appear in the future)
    !call time_push
    call fmm_sph_rotate_oz_adj_work(pm, vcos, vsin, alpha, src_m, 0.0_rp, tmp_ml)
    !call time_pull('oza')
    ! Perform rotation in the OXZ plane, work size is 4*pm*pm+13*pm+4
    !call time_push
    call fmm_sph_rotate_oxz_work(pm, ctheta, -stheta, 1.0_rp, tmp_ml, 0.0_rp, &
        & tmp_ml2, work)
    !call time_pull('oxz')
    ! OZ translation, workspace here is (pm+2)*(pm+1)
    !call time_push
    call fmm_m2l_ztranslate_work(rho, pm, pl, vscales, &
        & m2l_ztranslate_coef, 1.0_rp, tmp_ml2, 0.0_rp, tmp_ml, work)
    !call time_pull('zt')
    ! Backward rotation in the OXZ plane, work size is 4*pl*pl+13*pl+4
    !call time_push
    call fmm_sph_rotate_oxz_work(pl, ctheta, stheta, 1.0_rp, tmp_ml, 0.0_rp, &
        & tmp_ml2, work)
    !call time_pull('oxz')
    ! Backward rotation around OZ axis (work array might appear in the future)
    !call time_push
    call fmm_sph_rotate_oz_work(pl, vcos, vsin, 1.0_rp, tmp_ml2, beta, dst_l)
    !call time_pull('oz')
    !call time_pull("fmm_m2l_work")
end subroutine fmm_m2l_rotation_work

!> Direct L2L translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L \f$ is a matrix of a local-to-local
!! translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
!! @param[out] work: Temporary workspace of a size 6*p*p+19*p+8
subroutine fmm_l2l_rotation_work(c, p, vscales, vfact, alpha, &
    & src_l, beta, dst_l, work)
    ! Inputs
    integer, intent(in) :: p
    real(rp), intent(in) :: c(3), vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_l((p+1)*(p+1)), beta
    ! Output
    real(rp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Temporary workspace
    real(rp), intent(out), target :: work(6*p*p + 19*p + 8)
    ! Local variables
    real(rp) :: rho, ctheta, stheta, cphi, sphi
    integer :: m, n
    ! Pointers for temporary values of harmonics
    real(rp), pointer :: tmp_l(:), tmp_l2(:), vcos(:), vsin(:)
    ! Covert Cartesian coordinates into spherical
    call carttosph(c, rho, ctheta, stheta, cphi, sphi)
    ! If no need for rotations, just do translation along z
    if (abs(stheta) < eps_rp) then
        ! Workspace here is 2*(p+1)
        call fmm_l2l_ztranslate_work(c(3), p, vscales, vfact, &
            & alpha, src_l, beta, dst_l, work)
        return
    end if
    ! Prepare pointers
    m = (p+1)**2
    n = 4*m + 5*p ! 4*p*p + 13*p + 4
    tmp_l(1:m) => work(n+1:n+m) ! 5*p*p + 15*p + 5
    n = n + m
    tmp_l2(1:m) => work(n+1:n+m) ! 6*p*p + 17*p + 6
    n = n + m
    m = p + 1
    vcos => work(n+1:n+m) ! 6*p*p + 18*p + 7
    n = n + m
    vsin => work(n+1:n+m) ! 6*p*p + 19*p + 8
    ! Compute arrays of cos and sin that are needed for rotations of harmonics
    call trgev(cphi, sphi, p, vcos, vsin)
    ! Rotate around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_adj_work(p, vcos, vsin, alpha, src_l, 0.0_rp, tmp_l)
    ! Perform rotation in the OXZ plane, work size is 4*p*p+13*p+4
    call fmm_sph_rotate_oxz_work(p, ctheta, -stheta, 1.0_rp, tmp_l, 0.0_rp, &
        & tmp_l2, work)
    ! OZ translation, workspace here is 2*(p+1)
    call fmm_l2l_ztranslate_work(rho, p, vscales, vfact, 1.0_rp, &
        & tmp_l2, 0.0_rp, tmp_l, work)
    ! Backward rotation in the OXZ plane, work size is 4*p*p+13*p+4
    call fmm_sph_rotate_oxz_work(p, ctheta, stheta, 1.0_rp, tmp_l, 0.0_rp, tmp_l2, &
        & work)
    ! Backward rotation around OZ axis (work array might appear in the future)
    call fmm_sph_rotate_oz_work(p, vcos, vsin, 1.0_rp, tmp_l2, beta, dst_l)
end subroutine fmm_l2l_rotation_work


subroutine fmm_m2m(c_st, pm, s, t)
    implicit none

    real(rp) :: c_st(3)
    !! Distance vector from source to target
    integer(ip) :: pm
    !! Maximum level of spherical harmonics expansion for multipoles
    real(rp) :: s(:)
    !! Source distribution expansion coefficients
    real(rp) :: t(:)
    !! Target distribution expansion coefficients

    real(rp), allocatable :: work(:)

    ! Allocate local variables
    allocate(work(6*pm**2 + 19*pm + 8))

    call fmm_m2m_rotation_work(c_st, pm, vscales, vcnk, 1.0_rp, s, 0.0_rp, t, work)

    deallocate(work)
end subroutine

subroutine fmm_m2l(c_st, pm, pl, s, t)
    implicit none
    
    real(rp) :: c_st(3)
    !! Distance vector from source to target
    integer(ip) :: pl
    !! Maximum level of spherical harmonics expansion for local exp. 
    integer(ip) :: pm
    !! Maximum level of spherical harmonics expansion for multipoles
    real(rp) :: s(:)
    !! Source distribution expansion coefficients
    real(rp) :: t(:)
    !! Target distribution expansion coefficients

    real(rp), allocatable :: work(:)

    ! Allocate local variables
    !call time_push
    allocate(work(6*max(pm, pl)**2 + 19*max(pm, pl) + 8))
    !call time_pull('Memory allocation')
    !call time_push
    !call time_pull('Const gen')

    !call time_push
    call fmm_m2l_rotation_work(c_st, pm, pl, vscales, m2l_ztranslate_coef, 1.0_rp, s, 0.0_rp, t, work)
    !call time_pull("M2L rotation")

    deallocate(work)
end subroutine

subroutine fmm_l2l(c_st, r_s, r_t, pl, s, t)
    implicit none
    
    real(rp) :: c_st(3)
    !! Distance vector from source to target
    real(rp) :: r_s
    !! Size of source node
    real(rp) :: r_t
    !! Size of target node
    integer(ip) :: pl
    !! Maximum level of spherical harmonics expansion for local exp. 
    real(rp) :: s(:)
    !! Source distribution expansion coefficients
    real(rp) :: t(:)
    !! Target distribution expansion coefficients

    real(rp), allocatable :: vfact(:), work(:)

    ! Allocate local variables
    allocate(vfact(2*pl+1))
    allocate(work(6*pl**2 + 19*pl + 8))

    call make_vfact(pl, vfact)

    call fmm_l2l_rotation_work(c_st, pl, vscales, vfact, 1.0_rp, s, 0.0_rp, t, work)

    deallocate(vfact, work)
end subroutine
!> Accumulate potential, induced by multipole spherical harmonics
!!
!! This function relies on a user-provided temporary workspace
!!
!! Computes the following sum:
!! \f[
!!      v = \beta v + \alpha \sum_{\ell=0}^p \frac{4\pi}{\sqrt{2\ell+1}}
!!      \left( \frac{r}{\|c\|} \right)^{\ell+1} \sum_{m=-\ell}^\ell
!!      M_\ell^m Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ M \f$ is a vector of coefficients of input harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \alpha \f$ and \f$ \beta \f$ are scaling factors
!! and \f$ c \f$ is a location of a particle.
!!
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{-\ell} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^\ell \f$.
!!
!! @param[in] c: Coordinates of a particle (relative to center of harmonics)
!! @param[in] r: Radius of spherical harmonics
!! @param[in] p: Maximal degree of multipole basis functions
!! @param[in] vscales_rel: Relative normalization constants for
!!      \f$ Y_\ell^m \f$. Dimension is `(p+1)**2`
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole coefficients. Dimension is `(p+1)**2`
!! @param[in] beta: Scalar multiplier for `v`
!! @param[inout] v: Value of induced potential
!! @param[out] work: Temporary workspace of size (p+1)
subroutine fmm_m2p_work(c, src_r, p, vscales_rel, alpha, src_m, &
    & beta, dst_v, work)
    ! Inputs
    integer, intent(in) :: p
    real(rp), intent(in) :: c(3), src_r, vscales_rel((p+1)*(p+1)), alpha, &
        & src_m((p+1)*(p+1)), beta
    ! Output
    real(rp), intent(inout) :: dst_v
    ! Workspace
    real(rp), intent(out), target :: work(p+1)
    ! Local variables
    real(rp) :: rho, ctheta, stheta, cphi, sphi, rcoef, t, tmp, tmp1, tmp2, &
        & tmp3, max12, ssq12, pl2m, pl1m, plm, pmm, cmphi, smphi, ylm
    integer :: l, m, indl
    ! Scale output
    if (abs(beta) < eps_rp) then
        dst_v = 0.0_rp
    else
        dst_v = beta * dst_v
    end if
    ! In case of 0.0_rp alpha nothing else is required no matter what is the
    ! value of the induced potential
    if (abs(alpha) < eps_rp) then
        return
    end if
    ! Get spherical coordinates
    if (abs(c(1)) <  eps_rp) then
        max12 = abs(c(2))
        ssq12 = 1.0_rp
    else if (abs(c(2)) .gt. abs(c(1))) then
        max12 = abs(c(2))
        ssq12 = 1.0_rp + (c(1)/c(2))**2
    else
        max12 = abs(c(1))
        ssq12 = 1.0_rp + (c(2)/c(1))**2
    end if
    ! Then we compute rho
    if (abs(c(3)) < eps_rp) then
        rho = max12 * sqrt(ssq12)
    else if (abs(c(3)) .gt. max12) then
        rho = 1.0_rp + ssq12 *(max12/c(3))**2
        rho = abs(c(3)) * sqrt(rho)
    else
        rho = ssq12 + (c(3)/max12)**2
        rho = max12 * sqrt(rho)
    end if
    ! In case of a singularity (rho=0.0_rp) induced potential is infinite and is
    ! not taken into account.
    if (abs(rho) < eps_rp) then
        return
    end if
    ! Compute the actual induced potential
    rcoef = src_r / rho
    ! Length of a vector x(1:2)
    stheta = max12 * sqrt(ssq12)
    ! Case x(1:2) != 0
    if (abs(stheta) > eps_rp) then
        ! Normalize cphi and sphi
        cphi = c(1) / stheta
        sphi = c(2) / stheta
        ! Normalize ctheta and stheta
        ctheta = c(3) / rho
        stheta = stheta / rho
        ! Treat easy cases
        select case(p)
            case (0)
                ! l = 0
                dst_v = dst_v + alpha*rcoef*vscales_rel(1)*src_m(1)
                return
            case (1)
                ! l = 0
                tmp = src_m(1) * vscales_rel(1)
                ! l = 1
                tmp2 = ctheta * src_m(3) * vscales_rel(3)
                tmp3 = vscales_rel(4) * stheta
                tmp2 = tmp2 - tmp3*sphi*src_m(2)
                tmp2 = tmp2 - tmp3*cphi*src_m(4)
                dst_v = dst_v + alpha*rcoef*(tmp+rcoef*tmp2)
                return
        end select
        ! Now p>1
        ! Precompute alpha*rcoef^{l+1}
        work(1) = alpha * rcoef
        do l = 1, p
            work(l+1) = rcoef * work(l)
        end do
        ! Case m = 0
        ! P_{l-2}^m which is P_0^0 now
        pl2m = 1.0_rp
        dst_v = dst_v + work(1)*src_m(1)*vscales_rel(1)
        ! Update P_m^m for the next iteration
        pmm = -stheta
        ! P_{l-1}^m which is P_{m+1}^m now
        pl1m = ctheta
        ylm = pl1m * vscales_rel(3)
        dst_v = dst_v + work(2)*src_m(3)*ylm
        ! P_l^m for l>m+1
        do l = 2, p
            plm = dble(2*l-1)*ctheta*pl1m - dble(l-1)*pl2m
            plm = plm / dble(l)
            ylm = plm * vscales_rel(l*l+l+1)
            dst_v = dst_v + work(l+1)*src_m(l*l+l+1)*ylm
            ! Update P_{l-2}^m and P_{l-1}^m for the next iteration
            pl2m = pl1m
            pl1m = plm
        end do
        ! Prepare cos(m*phi) and sin(m*phi) for m=1
        cmphi = cphi
        smphi = sphi
        ! Case 0<m<p
        do m = 1, p-1
            ! P_{l-2}^m which is P_m^m now
            pl2m = pmm
            ylm = pmm * vscales_rel((m+1)**2)
            tmp1 = cmphi*src_m((m+1)**2) + smphi*src_m(m*m+1)
            dst_v = dst_v + work(m+1)*ylm*tmp1
            ! Temporary to reduce number of operations
            tmp1 = dble(2*m+1) * pmm
            ! Update P_m^m for the next iteration
            pmm = -stheta * tmp1
            ! P_{l-1}^m which is P_{m+1}^m now
            pl1m = ctheta * tmp1
            ylm = pl1m * vscales_rel((m+1)*(m+3))
            tmp1 = cmphi*src_m((m+1)*(m+3)) + smphi*src_m((m+1)*(m+2)+1-m)
            dst_v = dst_v + work(m+2)*ylm*tmp1
            ! P_l^m for l>m+1
            do l = m+2, p
                plm = dble(2*l-1)*ctheta*pl1m - dble(l+m-1)*pl2m
                plm = plm / dble(l-m)
                ylm = plm * vscales_rel(l*l+l+1+m)
                tmp1 = cmphi*src_m(l*l+l+1+m) + smphi*src_m(l*l+l+1-m)
                dst_v = dst_v + work(l+1)*ylm*tmp1
                ! Update P_{l-2}^m and P_{l-1}^m for the next iteration
                pl2m = pl1m
                pl1m = plm
            end do
            ! Update cos(m*phi) and sin(m*phi) for the next iteration
            tmp1 = cmphi
            cmphi = cmphi*cphi - smphi*sphi
            smphi = tmp1*sphi + smphi*cphi
        end do
        ! Case m=p requires only to use P_m^m
        ylm = pmm * vscales_rel((p+1)**2)
        tmp1 = cmphi*src_m((p+1)**2) + smphi*src_m(p*p+1)
        dst_v = dst_v + work(p+1)*ylm*tmp1
    ! Case of x(1:2) = 0 and x(3) != 0
    else
        ! In this case Y_l^m = 0 for m != 0, so only case m = 0 is taken into
        ! account. Y_l^0 = ctheta^l in this case where ctheta is either +1 or
        ! -1. So, we copy sign(ctheta) into rcoef. But before that we
        ! initialize alpha/r factor for l=0 case without taking possible sign
        ! into account.
        t = alpha * rcoef
        if (c(3) .lt. 0.0_rp) then
            rcoef = -rcoef
        end if
        ! Proceed with accumulation of a potential
        indl = 1
        do l = 0, p
            ! Index of Y_l^0
            indl = indl + 2*l
            ! Add 4*pi/(2*l+1)*rcoef^{l+1}*Y_l^0 contribution
            dst_v = dst_v + t*src_m(indl)*vscales_rel(indl)
            ! Update t
            t = t * rcoef
        end do
    end if
end subroutine fmm_m2p_work

subroutine fmm_m2p(c_st, r_s,  pm, s, potential)
    implicit none
    
    real(rp) :: c_st(3)
    !! Distance vector from source to target
    real(rp) :: r_s
    !! Size of source node
    integer(ip) :: pm
    !! Maximum level of spherical harmonics expansion for multipoles
    real(rp) :: s(:)
    !! Target distribution expansion coefficients
    real(rp), intent(out) :: potential

    real(rp), allocatable :: work(:)

    ! Allocate local variables
    allocate(work(pm+1))


    call fmm_m2p_work(c_st, r_s, pm, vscales_rel, 1.0_rp, s, 0.0_rp, potential, work)

    deallocate(work)
end subroutine

end module