subroutine new_coulomb_kernel(dr, maxder, res)
    !! This function compute the coulomb kernel for the distance vector dr and 
    !! its derivative up to the value required by maxder.
    use mod_memory, only: ip, rp

    implicit none

    real(rp) :: dr(3)
    !! Distance vector from atom i to atom j
    integer(ip), intent(in) :: maxder
    !! Maximum derivative to be computed
    real(rp), intent(out), dimension(maxder+1) :: res
    !! Results vector
    
    integer(ip) :: ii
    real(rp) :: norm2_r, inv_norm_sq 

    if(maxder < 0) then
        ! Strange request, maybe a warning should be printed
        return 
    end if
    
    norm2_r = sqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))
    res(1) = 1.0_rp / norm2_r
    inv_norm_sq = res(1) * res(1)

    do ii = 1, maxder
        res(ii+1) = res(ii) * inv_norm_sq
    end do
    
end subroutine new_coulomb_kernel

subroutine new_damped_coulomb_kernel(i, j, maxder, res, dr)
    !! This subroutine computes the damped coulomb kernel between two atoms.
    !! Note that this only makes sense between two MM atoms, as it is only used
    !! to compute the field that induces the point dipoles!

    use mod_mmpol, only: amoeba, thole, fatal_error, cmm
    use mod_memory, only: ip, rp
    use mod_constants, only: eps_rp
    
    implicit none

    integer(ip), intent(in) :: i, j
    !! Index of atoms (in MM atom list) for which the kernel should be
    !! computed
    integer(ip), intent(in) :: maxder
    !! Maximum derivative to be computed
    real(rp), intent(out), dimension(maxder+1) :: res
    !! Results vector
    real(rp), intent(out), dimension(3) :: dr
    !! Distance vector between i and j
    
    real(rp) :: s, u, u3, u4, fexp, eexp
    
    s = thole(i) * thole(j)
    
    ! Compute undamped kernels
    dr = cmm(:,j) - cmm(:,i)
    call new_coulomb_kernel(dr, maxder, res)

    if(abs(s) < eps_rp) then
        ! either thole(i) or thole(j) are zero, so the damped kernel
        ! is just the same as the regular one. Job is done.
        return
    end if

    u = 1.0_rp / (res(1) * s)
    u3 = u**3

    if(amoeba .and. u3 * 0.39_rp < 50.0_rp) then
        ! 1.0 / (res(1) * s)**3 * 0.39_rp < 50_rp TODO Remove 0.39, this is a
        ! FF constants, ask someone why this condition is here.
        ! Here basically we multiply the standard interactions kernel for
        ! dumping coefficients. The equations implemented here correspond to 
        ! eq. (5) of 10.1021/jp027815
        fexp = -0.39_rp * u3
        eexp = exp(fexp)
        if(maxder >= 1) res(2) = res(2) * (1.0_rp - eexp)
        if(maxder >= 2) res(3) = res(3) * (1.0_rp - (1.0_rp - fexp) * eexp)
        if(maxder >= 3) res(4) = res(4) * &
                        (1.0_rp - (1.0_rp - fexp + 0.6_rp * fexp * fexp) * eexp)
        if(maxder >= 4) res(5) = res(5) * &
                        (1.0_rp - (1.0_rp - fexp + 18.0_rp/35.0_rp * fexp *  fexp - &
                        9.0_rp/35.0_rp * fexp * fexp * fexp) * eexp)
        if(maxder >= 5) then
            call fatal_error("Damped Coulomb kernel (AMOEBA) only supports up to the 5th derivative")
        end if
    else if(.not. amoeba .and. res(1) > 1_rp/s) then
        ! TODO Again it is not clear to me why condition res(1) > 1_rp/s is here.
        u4 = u3*u
        if(maxder >= 1) res(2) = res(2) * (4.0_rp * u3 - 3.0_rp * u4)
        if(maxder >= 2) res(3) = res(3) * u4
        if(maxder >= 3) res(4) = res(4) * 0.2_rp * u4
        if(maxder >= 4) then
            call fatal_error("Damped Coulomb kernel (WANG) only supports up to the 4th derivative")
        end if
    end if
end subroutine new_damped_coulomb_kernel

subroutine coulomb_kernel(damp,maxder,x,y,z,thole_i,thole_j,rm1,rm3,rm5,rm7,rm9,rm11)
  use mod_mmpol, only: ff_type
  use mod_constants, only: zero, one, three, four, nine
  use mod_memory, only: ip, rp
  implicit none
!
! given in input the cartesian coordinates of a distance vector x, y, z,
! and, if required, the thole factors of the two atoms, computes f(r)/r 
! and its derivatives up to the fifth, where f(r) is a damping function.
!
! the damping functions currently implemented are
!
! - no damping, which is used for the interactions between static sources
!   in both MMPol and AMOEBA;, for derivatives up to the fifth
!
! - Wang's polynomial damping for the field and field gradient (required
!   for MMPol energy and forces) produced by static charges and induced
!   dipoles, corresponding to derivatives up to the third.
!
!   this function is defined starting from the 1/r^3 term as follows:
!
!           / 4*(r/s)**3 - 3*(r/s)**4, if 0 < r < s
!     fl3 = 
!           \ 1,                       if r >= s
!
!   where s = thole_i*thole_j. 
!   higher order terms are computed by differentiating fl3/r**3 and collecting the
!   terms that multiply the same inverse power of r. 
!   we get, for 0 < r < s (otherwise all terms are one):
!
!     fl5 = (r/s)**4
!     fl7 = 0.2*(r/s)**4.
!
! - AMOEBA's exponential damping for the field and field gradient produced
!   by (static) multipoles up to the quadrupole and induced dipoles, 
!   corresponding to derivatives up to the fourth.
!
!   this function is defined starting from the 1/r^3 term as follows:
!
!     fl3 = 1 - exp(-a * (r/s)**3)
!
!   where s = thole_i*thole_j and a = 0.39.
!   higher order terms are computed by differentiating fl3/r**3 and collecting the
!   terms that multiply the same inverse power of r. 
!   we get:
!
!     fl5 = 1 - (1 + a * (r/s)**3                                                   ) * exp(-a * (r/s)**3)
!     fl7 = 1 - (1 + a * (r/s)**3 +   3/5 * a**2 * (r/s)**6                         ) * exp(-a * (r/s)**3)
!     fl9 = 1 - (1 + a * (r/s)**3 + 18/35 * a**2 * (r/s)**6 + 9/35 * a**3 * (r/s)**9) * exp(-a * (r/s)**3)
!
  logical,     intent(in)    :: damp
  integer(ip), intent(in)    :: maxder
  real(rp),    intent(in)    :: x, y, z, thole_i, thole_j
  real(rp),    intent(inout) :: rm1, rm3, rm5, rm7, rm9, rm11
!
  logical  :: w_damp, exp_damp
  real(rp) :: r, r2, rm2, fl3, fl5, fl7, fl9, s, u, v3, v4, fexp, eexp, um3, um5, um7, um9
  
  real(rp), parameter :: pt2 = 0.2_rp, pt39 = 0.39_rp, pt6 = 0.6_rp, f18 = 18.0_rp, &
                         f35  = 35.0_rp, f50 = 50.0_rp
    !TODO Miche: PT39 è un parametro di FF in realtà

!
  r2   = x*x + y*y + z*z
  r    = sqrt(r2)
  rm1  = one/r
!
  fl3  = one
  fl5  = one
  fl7  = one
  fl9  = one
  w_damp   = .false.
  exp_damp = .false.
!
  if (damp) then
!
!   check that the damping is actually active:
!
    s = thole_i*thole_j
    if (s.gt.zero) then
      u  = r/s
      v3 = u**3
      w_damp   = ff_type.eq.0 .and. u.lt.one 
      exp_damp = ff_type.eq.1 .and. pt39*v3.lt.f50
    end if
  end if
  if (maxder.ge.1) then
!
! compute rm3, possibly including the appropriate damping function:
!
    rm2 = rm1*rm1
    um3 = rm1*rm2
    if (w_damp) then
      v4  = v3*u 
      fl3 = four*v3 - three*v4
    else if(exp_damp) then
      fexp = -pt39*v3
      eexp = exp(fexp)
      fl3  = one - eexp
    end if
    rm3 = fl3*um3
  end if
!
  if (maxder.ge.2) then
!
! compute rm5, possibly including the appropriate damping function:
!
    um5 = um3*rm2
    if (w_damp) then
      fl5 = v4
    else if(exp_damp) then
      fl5  = one - (one - fexp)*eexp
    end if
    rm5 = fl5*um5
  end if
!
  if (maxder.ge.3) then
!
! compute rm5, possibly including the appropriate damping function:
!
    um7 = um5*rm2
    if (w_damp) then
      fl7 = pt2*v4
    else if(exp_damp) then
      fl7  = one - (one - fexp + pt6*fexp*fexp)*eexp
    end if
    rm7 = fl7*um7
  end if
!
  if (maxder.ge.4) then
!
! compute rm9, possibly including the appropriate damping function:
!
    um9 = um7*rm2
    if(exp_damp) then
      fl9  = one - (one - fexp + f18*fexp*fexp/f35 - nine*fexp*fexp*fexp/f35)*eexp
    end if
    rm9 = fl9*um9
  end if
!
  if (maxder.ge.5) then
!
! compute rm11. no damping currently used.
!
    rm11 = um9*rm2
  end if
!
  
end subroutine coulomb_kernel
