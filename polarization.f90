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
real(rp),dimension(3*pol_atoms) :: e_vec, ipd0
integer(ip) :: I,n
logical :: status

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
if (Amoeba) then
    e_vec = RESHAPE(e(:,:,2), (/3*pol_atoms/))    ! the polarization field for Amoeba FF
else
    e_vec = RESHAPE(e(:,:,1), (/3*pol_atoms/))
end if 
    
!
! Initial guess for iterative solvers 
!
ipd0 = Zero
if (iscr.le.2) then
    call PolVec(3*pol_atoms,ipd0,ipd0)
!    do I = 1, pol_atoms
!        ipd0(3*I-2) = pol(I)*e(1,I,1)   ! ipd0x
!        ipd0(3*I-1) = pol(I)*e(2,I,1)   ! ipd0y
!        ipd0(3*I)   = pol(I)*e(3,I,1)   ! ipd0z
!    enddo
end if
        
!
! Decide which solver to use for computing induced point dipoles 
!
if  (iscr.eq.1) then            ! preconditioned conjugate gradient
    call conjugate_gradient(n,verbose,convergence,e_vec,ipd0,nmax,status,TMatVec,PolVec)
elseif (iscr.eq.2) then         ! jacobi iterations with DIIS extrapolation
    call jacobi_diis(n,verbose,10,2,convergence,e_vec,ipd0,nmax,status,TMatVec_offdiag,PolVec)
elseif (iscr.eq.3) then         !  matrix inversion
    call induce_dipoles_MatInv(e_vec,ipd0)  !TODO: read polarization matrix TMat and not calculate it separately
else
    write(*,*) "Unknown solver for calculation of the induced point dipoles in polarization"
end if
        
!        
! For Amoeba FF compute induced dipoles by the direct field
!
if (Amoeba) then
    do I=1,pol_atoms
        ipd( 1, I, 1) = pol(I)*e(1,I,1)
        ipd( 2, I, 1) = pol(I)*e(2,I,1)
        ipd( 3, I, 1) = pol(I)*e(3,I,1)
    enddo
end if
    
!
! Reshape dipole vector into the matrix 
!
if (Amoeba) then
    ipds( :, :, 2) = RESHAPE(ipd0, (/3,pol_atoms/)) 
else
    ipds( :, :, 1) = RESHAPE(ipd0, (/3,pol_atoms/)) 
end if

!
!   Deallocate dipole polarization matrix when is not needed
!
call r_free2('polarization [TMat]',TMat)

end subroutine polarization
