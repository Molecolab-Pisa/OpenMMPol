subroutine coulomb_kernel(damp,maxder,x,y,z,thole_i,thole_j,rm1,rm3,rm5,rm7,rm9,rm11)
  use mmpol
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
  
  real(rp), parameter :: one = 1.0_rp, three = 3.0_rp, five = 5.0_rp, seven = 7.0_rp, pt2 = 0.2_rp, &
                         pt39 = 0.39_rp, pt6 = 0.6_rp, four = 4.0_rp, nine = 9.0_rp, f18 = 18.0_rp, &
                         f35  = 35.0_rp, f50 = 50.0_rp

!
  r2   = x*x + y*y + z*z
  r    = sqrt(r2)
  rm1  = one/r
!
  fl3  = one
  fl5  = one
  fl7  = one
  fl9  = one
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
      fl9  = one - (one - fexp + f18*fexp*fexp/f35 + f9*fexp*fexp*fexp/35)*eexp
    end if
    rm9 = fl9*um9
  end if
!
  if (maxder.ge.5) then
!
! compute rm11. no damping currently used.
!
    rm9 = um9*rm2
  end if
!
  
end subroutine coulomb_kernel
