subroutine polarization(iscr,e,ipds)
  use mmpol
  use polar
  use solvers
  implicit none
!
! main driver for the calculation of induced dipoles.
! ============================================================
!
!   iscr ... 1: Induced dipoles are computed by conjugate gradient
!               solver 
!            2: Induced dipoles are computed by jacobi solver
!            3: Induced dipoles are computed by matrix inversion.  
!  
!   e    ... Total electric field at MMpol sites.
!               - for amoeba, e has dimensions (3,pol_atoms,2). 
!                 the direct field is in e(:,:,1), the polarization
!                 field in e(:,:,2).
!               
!               - for amber, e has dimension (3,pol_atoms,1).
! 
! The dipoles are induced by the electric field of the static multipoles.
! For the Amoeba force-field only polarization field and p-dipoles are 
! considered.
!   - For AMBER force-field the induced dipoles will be stored in:
!
!       ipd(:,:,1)
!
!   - For Amoeba force-field the direct and polarization dipoles are
!
!       ipd(:,:,1)      ! d dipoles induced by direct field
!       ipd(:,:,2)      ! p dipoles induced by polarization field 
!                         and other p dipoles.
!

integer(ip),intent(in) :: iscr
real(rp), dimension(3,pol_atoms,n_ipd), intent(in) :: e
real(rp), dimension(3,pol_atoms,n_ipd), intent(out) :: ipds
!
! internal variables:
!
real(rp),dimension(3*pol_atoms) :: ep_vec, ed_vec, ipd0_p, ipd0_d
integer(ip) :: I,n
logical :: status
!
integer(ip), parameter :: diis_max = 20, norm_jacobi = 2

!
! Dimension of the system
!
n = 3*pol_atoms

!
! Allocate dipole polarization matrix (deallocated at the end of subroutine)
!
call r_alloc2('polarization [TMat]',n,n,TMat)
    
!
! Compute dipole polarization matrix
!
call create_TMat(TMat)
    
!
! Reshape electric field matrix into a vector
!
ed_vec = RESHAPE(e(:,:,1), (/3*pol_atoms/))    ! the direct field for Amoeba FF or AMBER
if (Amoeba) then
    ep_vec = RESHAPE(e(:,:,2), (/3*pol_atoms/))    ! the polarization field for Amoeba FF
end if 
    
!
! Initial guess for iterative solvers (zero for matrix inversion)
!
ipd0_p = Zero   ! Initial guess for amoeba p dipoles
ipd0_d = Zero   ! Initial guess for amber or amoeba d dipoles
if (iscr.le.2) then
    call PolVec(3*pol_atoms,ed_vec,ipd0_d)
    if (Amoeba) call PolVec(3*pol_atoms,ep_vec,ipd0_p)
end if
!
! Decide which solver to use for computing induced point dipoles 
!
if  (iscr.eq.1) then            ! preconditioned conjugate gradient
    call conjugate_gradient(n,verbose,convergence,ed_vec,ipd0_d,nmax,status,TMatVec,PolVec)
    if (Amoeba) call conjugate_gradient(n,verbose,convergence,ep_vec,ipd0_p,nmax,status,TMatVec,PolVec)
elseif (iscr.eq.2) then         ! jacobi iterations with DIIS extrapolation
    call jacobi_diis(n,verbose,diis_max,norm_jacobi,convergence,ed_vec,ipd0_d,nmax,status,TMatVec_offdiag,PolVec)
    if (Amoeba) call jacobi_diis(n,verbose,diis_max,norm_jacobi,convergence,ep_vec,ipd0_p,nmax,status,TMatVec_offdiag,PolVec)
elseif (iscr.eq.3) then         !  matrix inversion
    call induce_dipoles_MatInv(ed_vec,ipd0_d)  !TODO: read polarization matrix TMat and not calculate it separately
    if (Amoeba) call induce_dipoles_MatInv(ep_vec,ipd0_p)
else
    write(*,*) "Unknown solver for calculation of the induced point dipoles in polarization"
end if

    
!
! Reshape dipole vector into the matrix 
!
if (Amoeba) then
    ipds( :, :, 2) = RESHAPE(ipd0_p, (/3,pol_atoms/)) 
    ipds( :, :, 1) = RESHAPE(ipd0_d, (/3,pol_atoms/)) 
else
    ipds( :, :, 1) = RESHAPE(ipd0_d, (/3,pol_atoms/)) 
end if

!
!   Deallocate dipole polarization matrix when is not needed
!
call r_free2('polarization [TMat]',TMat)

end subroutine polarization
