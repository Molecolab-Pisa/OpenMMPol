module mod_polarization
    use mod_memory, only: ip, rp
    
    implicit none 
    private
    
    real(rp), allocatable :: TMat(:,:)

    public :: polarization
    
    contains
    
    subroutine polarization(e, ipds, arg_solver, arg_mvmethod)
        !! Main driver for the calculation of induced dipoles. 
        !! Takes electric field at induced dipole sites as input and -- if
        !! solver converges -- provides induced dipoles as output.
        !! Since AMOEBA requires the calculations of two sets of induced dipoles
        !! generated from two different electric fields (normally called direct (D) 
        !! and polarization (P)) both electric field and induced dipoles are shaped
        !! with an extra dimension and this routine calls the solver twice to 
        !! solve the two linear systems in the case of AMOEBA FF. Direct electric
        !! field and induced dipoles are stored in e(:,:,1)/ipds(:,:,1) while
        !! polarization field/dipole are stored in e(:,:,2)/ipds(:,:,2).

        use mod_mmpol, only: amoeba, pol_atoms, n_ipd, fatal_error
        use mod_solvers, only: jacobi_diis_solver, conjugate_gradient_solver, &
                               inversion_solver
        use mod_memory, only: ip, rp, mallocate, mfree
        use mod_constants, only: OMMP_MATV_DEFAULT, &
                                 OMMP_SOLVER_DEFAULT, &
                                 OMMP_SOLVER_CG, &
                                 OMMP_SOLVER_DIIS, &
                                 OMMP_SOLVER_INVERSION
      
        implicit none

        real(rp), dimension(3, pol_atoms, n_ipd), intent(in) :: e
        !! Total electric field that induces the dipoles
        real(rp), dimension(3, pol_atoms, n_ipd), intent(out) :: ipds
        !! Induced dipoles

        integer(ip), intent(in), optional :: arg_solver
        !! Flag for the solver to be used; optional, should be one OMMP_SOLVER_
        !! if not provided [[OMMP_SOLVER_DEFAULT]] is used.
        integer(ip), intent(in), optional :: arg_mvmethod
        !! Flag for the matrix-vector method to be used; optional, should be one of
        !! OMMP_MATV_ if not provided [[OMMP_MATV_DEFAULT]] is used.
        
        real(rp), dimension(3*pol_atoms) :: ep_vec, ed_vec, ipd0_p, ipd0_d
        integer(ip) :: n, solver, mvmethod
        
        ! Handling of optional arguments
        if(present(arg_solver)) then
            solver = arg_solver
        else
            solver = OMMP_SOLVER_DEFAULT
        end if

        if(present(arg_mvmethod)) then
            mvmethod = arg_mvmethod
        else
            mvmethod = OMMP_MATV_DEFAULT
        end if

        ! Dimension of the system
        n = 3*pol_atoms

        ! Allocate and compute dipole polarization matrix
        call mallocate('polarization [TMat]',n,n,TMat)
        call create_TMat(TMat)

        ! Reshape electric field matrix into a vector
        ! direct field for Wang and Amoeba
        ed_vec = reshape(e(:,:,1), (/ n /))
        ! polarization field just for Amoeba
        if(amoeba) ep_vec = reshape(e(:,:,2), (/ n /))
        
        ! Initial guess for iterative solvers (zero for matrix inversion)

        ipd0_p = 0.0_rp   ! Initial guess for amoeba p dipoles
        ipd0_d = 0.0_rp   ! Initial guess for wang or amoeba d dipoles
        
        if(solver /= OMMP_SOLVER_INVERSION) then
            call PolVec(3*pol_atoms, ed_vec, ipd0_d)
            if(amoeba) call PolVec(3*pol_atoms, ep_vec, ipd0_p)
        end if

        select case (solver)
            case(OMMP_SOLVER_CG)
                call conjugate_gradient_solver(n, ed_vec, ipd0_d, TMatVec, PolVec)
                if(amoeba) call conjugate_gradient_solver(n, ep_vec, ipd0_p, TMatVec, PolVec)

            case(OMMP_SOLVER_DIIS)
                call jacobi_diis_solver(n, ed_vec, ipd0_d, TMatVec_offdiag, PolVec)
                if(amoeba) call jacobi_diis_solver(n, ep_vec, ipd0_p, TMatVec_offdiag, PolVec)

            case(OMMP_SOLVER_INVERSION)
                call inversion_solver(n, ed_vec, ipd0_d, TMat)  
                if(amoeba) call inversion_solver(n, ep_vec, ipd0_p, TMat)
                
            case default
                call fatal_error("Unknown solver for calculation of the induced point dipoles") 
        end select
        
        ! Reshape dipole vector into the matrix 
        ipds(:,:,1) = reshape(ipd0_d, (/3_ip, pol_atoms/)) 
        if(amoeba) ipds(:,:,2) = reshape(ipd0_p, (/3_ip, pol_atoms/)) 

        call mfree('polarization [TMat]',TMat)

    end subroutine polarization
    
    subroutine dipole_T(scalef,I,J,TTens)
        use mod_mmpol, only : thole, cpol, pol_atoms, pol, polar_mm
        implicit none
        !                      
        ! Compute element of the polarization tensor TTens between
        ! polarizable cpol atom I and polarizable cpol atom J. On the
        ! TTens diagonal (I=J) are inverse polarizabilities and on the 
        ! off-diagonal dipole field.
        !
        ! Polarizabilities pol are defined for polarizable atoms only while 
        ! Thole factors are defined for all of them
        !
        real(rp), dimension(3,pol_atoms,3,pol_atoms), intent(inout) :: TTens
        real(rp), intent(in) :: scalef
        integer(ip), intent(in) :: I,J
        !
        real(rp)    :: dx, dy, dz, rm1, rm3, rm5, rm7, rm9, rm11
        !
        real(rp), parameter   :: Zero = 0.0_rp, One = 1.0_rp, Three = 3.0_rp
        !
         
        ! Initialize variables
        rm1 = Zero
        rm3 = Zero
        rm5 = Zero
        rm7 = Zero
        rm9 = Zero
        rm11 = Zero
        
        if (I.eq.J) then
            TTens(1,I,1,I) = One/pol(I)
            TTens(2,I,2,I) = One/pol(I)
            TTens(3,I,3,I) = One/pol(I)
        else
            
            ! intermediates
            dx = cpol(1,J) - cpol(1,I)
            dy = cpol(2,J) - cpol(2,I)
            dz = cpol(3,J) - cpol(3,I)

            call coulomb_kernel(.true.,2,dx,dy,dz,thole(polar_mm(I)),thole(polar_mm(J)),rm1,rm3,rm5,rm7,rm9,rm11)

            ! Apply scaling
            rm3 = rm3*scalef
            rm5 = rm5*scalef

            ! Fill the matrix elemets
            TTens(1,J,1,I) = TTens(1,J,1,I) + rm3 - Three*rm5*dx*dx
            TTens(2,J,1,I) = TTens(2,J,1,I) - Three*rm5*dx*dy
            TTens(3,J,1,I) = TTens(3,J,1,I) - Three*rm5*dx*dz

            TTens(1,J,2,I) = TTens(1,J,2,I) - Three*rm5*dy*dx
            TTens(2,J,2,I) = TTens(2,J,2,I) + rm3 - Three*rm5*dy*dy
            TTens(3,J,2,I) = TTens(3,J,2,I) - Three*rm5*dy*dz

            TTens(1,J,3,I) = TTens(1,J,3,I) - Three*rm5*dz*dx
            TTens(2,J,3,I) = TTens(2,J,3,I) - Three*rm5*dz*dy
            TTens(3,J,3,I) = TTens(3,J,3,I) + rm3 - Three*rm5*dz*dz
        end if
    end subroutine dipole_T
        
        

    subroutine create_TMat(TMat)
        use mod_mmpol, only : pol_atoms, uscale, mm_polar, amoeba, &
                              verbose, polar_mm, conn, pg_conn, mmat_polgrp, &
                              polgrp_mmat

        use mod_io, only: print_matrix
        use mod_constants, only: eps_rp

        implicit none
        !                      
        ! Construct polarization tensor for matrix inversion solution
        ! of the induced dipoles by matrix inversion. On the diagonal
        ! of the TMat there are inverse polarizabilities and on the 
        ! off-diagonal dipole field.
        !
        ! The tensor is constructed with dimension (3,pol_atoms,3,pol_atoms)
        ! then we subtract the excluded interaction and reshape the tensor
        ! into matrix with dimension (3*pol_atoms,3*pol_atoms) which is 
        ! the one on the output.
        !
        ! Polarizabilities pol are defined for polarizable atoms only while 
        ! n12, i12... factors are defined for all of them (even mm atoms)
        !
        real(rp), dimension(3*pol_atoms,3*pol_atoms), intent(out) :: TMat
        !
        real(rp), dimension(3,pol_atoms,3,pol_atoms) :: TTens
        integer(ip) :: I, J, K, L, IJ, ineigh, ipg
        !
        real(rp), parameter   :: Zero = 0.0_rp, One = 1.0_rp
        !
        
        ! Initialize the tensor with zeros
        TTens = Zero
        
        ! Assemble the diagonal elements of the MMpol matrix
    !    do I = 1,pol_atoms
    !        do J = 1,3
    !            TTens(J,I,J,I) = One/pol(I)
    !        enddo
    !    enddo
        
        ! Assemble the upper triangular MMPol matrix and diagonal
        do I = 1,pol_atoms
            do J = 1,I !-1
                
                call dipole_T(One,I,J,TTens)
                
            enddo
        enddo
        
        ! Fill the lower diagonal
        do I = 1, pol_atoms
            do J = I+1, pol_atoms
              do K = 1, 3
                do L = 1, 3
                    TTens(L,J,K,I) = TTens(K,I,L,J)
                enddo
              enddo
            enddo
        enddo
        
        ! Remove excluded interactions
        if  (Amoeba) then      ! AMOEBA and sources are induced dipoles (only from p dipoles)
            do I = 1,pol_atoms
                ! For polarization group pg_conn(1) is the identity matrix
                do ineigh=1, 4
                    if (abs(uscale(ineigh) - 1.0_rp) > eps_rp) then
                        do ipg=pg_conn(ineigh)%ri(mmat_polgrp(i)), &
                               pg_conn(ineigh)%ri(mmat_polgrp(i)+1)-1

                            do ij = polgrp_mmat%ri(ipg), &
                                    polgrp_mmat%ri(ipg+1)-1
                                j = polgrp_mmat%ci(ij)

                                if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J

                                call dipole_T(uscale(ineigh)-One,I, mm_polar(J) ,TTens)
                            end do
                        enddo
                    end if 
                end do
            enddo
            
        else     
            ! AMBER and sources are induced dipoles
            do I = 1,pol_atoms
                do ineigh=1, 4
                    ! Field from dipoles is scaled by uscale
                    if(abs(uscale(ineigh) - one) > eps_rp) then
                        do ij = conn(ineigh)%ri(polar_mm(i)), &
                                conn(ineigh)%ri(polar_mm(i)+1)-1
                            j = conn(ineigh)%ci(ij)
                            
                            if (mm_polar(J) /= 0) then
                                ! if J is not polarizable atom, skip J
                                ! mm_polar changes mm index index into pol index
                            
                                call dipole_T(uscale(ineigh)-One,I,mm_polar(J),TTens)
                            end if
                        enddo
                    end if
                end do
            enddo 
        end if

        ! Convert polarization tensor to the matrix
        TMat = RESHAPE(TTens, (/3*pol_atoms, 3*pol_atoms/))  
        
        ! Print the matrix if verbose output is requested
        if (verbose >= 2) then
            call print_matrix(.true.,'K:',3*pol_atoms,3*pol_atoms,3*pol_atoms,3*pol_atoms,TMat)
        end if
        
    end subroutine create_TMat

    subroutine TMatVec(n,x,y)
        use mod_constants, only : One, Zero
        implicit none
        ! Perform matrix vector multiplication y = TMat*x,
        ! where TMat is polarization matrix and x,y are 
        ! vectors
        !
        !   variables:
        !  
        !     n        : integer, input, size of the matrix
        ! 
        !     x        : real, vector which is multiplied by the matrix
        !
        !     y        : real, solution vector
        !
        integer(ip), intent(in) :: n
        real(rp),dimension(n),intent(in)   :: x
        real(rp),dimension(n), intent(out) :: y
        !
        ! initialize solution vector
        !
        y = 0
        !
        ! Compute the matrix vector product
        !
        call DGEMM('N','N',n,1,n,One,TMat,n,x,n,Zero,y,n) !XGEMM
        !
    end subroutine TMatVec
        
    subroutine TMatVec_offdiag(n,x,y)
        use mod_constants, only : One, Zero
        implicit none
        ! Perform matrix vector multiplication y = [TMat-diag(TMat)]*x,
        ! where TMat is polarization matrix and x,y are 
        ! vectors
        !
        !   variables:
        !  
        !     n        : integer, input, size of the matrix
        ! 
        !     x        : real, vector which is multiplied by the matrix
        !
        !     y        : real, solution vector
        !
        integer(ip), intent(in) :: n
        real(rp),dimension(n),intent(in)   :: x
        real(rp),dimension(n), intent(out) :: y
        !
        real(rp),dimension(n,n) :: O
        integer(ip) :: I
        !
        ! initialize solution vector and matrix
        !
        y = 0
        O = TMat
        !
        ! Subtract the diagonal from the O matrix
        !
        do I = 1, n
            O(I,I) = Zero
        enddo
        !
        ! Compute the matrix vector product
        !
        call DGEMM('N','N',n,1,n,One,O,n,x,n,Zero,y,n) !XGEMM
        !
    end subroutine TMatVec_offdiag

    subroutine PolVec(n,x,y)
        use mod_mmpol, only : pol
        implicit none
        ! Perform matrix vector multiplication y = pol*x,
        ! where pol polarizability vector and x,y are 
        ! vectors
        !
        !   variables:
        !  
        !     n        : integer, input, size of the matrix
        ! 
        !     x        : real, vector which is multiplied by the 
        !                inverse diagonal matrix
        !
        !     y        : real, solution vector
        !
        integer(ip), intent(in) :: n
        real(rp),dimension(n),intent(in)   :: x
        real(rp),dimension(n), intent(out) :: y
        !
        integer(ip) :: I,indx
        !
        ! initialize solution vector
        !
        y = 0
        !
        ! Compute the matrix vector product
        !
        do I = 1, n
            indx = (I+2)/3
            y(I) = pol(indx)*x(I)   
        enddo
        !
    end subroutine PolVec

end module mod_polarization
