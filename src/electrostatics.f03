subroutine electrostatics(isrc,itrg,ider,v,e,dv,de)
  use mmpol
  use elstat
  use mod_memory, only: ip, rp
  implicit none
!
! main driver for the calculation of electrostatic properties.
! ============================================================
!
!   isrc ... 0: the sources are the charges (multipoles), contained in the array
!               q and with coordinates cmm (see the mmpol module)
!            1: the sources are the ipds, contained in the array ipd and with
!               coordinates cpol (see the mmpol module)
!
!   itrg ... 0: compute the desired electrostatic property at cmm
!            1: compute the desired electrostatic property at cpol
!
!   ider ... 0: do not differentiate. 
!               - if itrg is 0, compute the electrostatic properties of the same
!                 order than the sources in q, i.e., the potential if q are charges,
!                 the potential, field and field gradient if q are multipoles up to
!                 the quadrupole (AMOEBA). 
!                 the output is given in v in the following order:
!
!                    v, ex, ey, ez, gxx, gxy, gyy, gxz, gyz, gzz
!
!               - if itrg is 0, compute the electric field of the source at the 
!                 polarizable sites. note that for AMOEBA, this requires computing
!                 two fields: the direct and polarization field if the sources are
!                 the multipoles, or the field of both sets of dipoles.
!                 the output is given in e in the following order:
!
!                   ex, ey, ez
!
!                 for amoeba, e has dimensions (3,pol_atoms,2). 
!                 the direct field is in e(:,:,1), the polarization field in e(:,:,2).
!
!        ... 1: differentiate: compute the derivatives of the quantities described before
!               with respect to the position of the target.
!
!               - if itrg is 0, the output is given in dv in the following order:
!
!                   ex, ey, ez, gxx, gxy, gyy, gxz, gyz, gzz, hxxx, hxxy, hxxz, hxyy, hxyz,
!                   hxzz, hyyy, hyyz, hyzz, hzzz.
!               
!               - if itrg is 1, the output is given in de in the following order:
!
!                   gxx, gxy, gyy, gxz, gyz, gzz
!
!                 note again that for AMOEBA, de has dimensions (6,pol_atoms,2).
!                 the direct field gradient is in de(:,:,1) and the polarization field
!                 gradient in de(6,pol_atoms,2). 
!
! - For AMOEBA force-field, the scaling of the nearest neighbor 
!   interactions for p and d dipoles is reversed dscale is used for 
!   scaling potential of p dipoles and pscale is used for scaling 
!   potential of d dipoles.
!
  integer(ip) :: isrc,itrg,ider
  real(rp),   dimension(ld_cart,mm_atoms),   intent(inout) :: v
  real(rp),   dimension(ld_cder,mm_atoms),   intent(inout) :: dv
  real(rp),   dimension(3,pol_atoms,n_ipd), intent(inout) :: e
  real(rp),   dimension(6,pol_atoms,n_ipd), intent(inout) :: de
!
! internal variables:
!

!
! figure how what to compute
!
  if (ider.eq.0) then         ! Do not differentiate
    if (itrg.eq.0) then       ! Calculate potential at MM sites
      call multipoles_potential(isrc,v)
    elseif (itrg.eq.1) then   ! Calculate electric field at MMpol sites
      call multipoles_field(isrc,e)
    else
      write(*,*) "Unknown value of itrg in electrostatics."
    end if 
  elseif (ider.eq.1) then     ! Calculate derivatives
    if (itrg.eq.0) then       ! Calculate derivative of the potential at MM sites
      call multipoles_potential_deriv(isrc,dv)
    elseif (itrg.eq.1) then   ! Calculate electric field derivatives at MMpol sites\
      call multipoles_field_deriv(isrc,de)
    else
      write(*,*) "Unknown value of itrg in electrostatics."
    end if
  else
    write(*,*) "Unknown value of ider in electrostatics."
  end if

end subroutine electrostatics
