module mod_polarization
    implicit none 
    
    private

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

        use mod_mmpol, only: amoeba, pol_atoms, n_ipd, verbose, fatal_error
        use polar  
        use mod_solvers, only: jacobi_diis, conjugate_gradient
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
        logical :: status
        
        integer(ip), parameter :: diis_max = 20, norm_jacobi = 2 ! TODO
        real(rp) :: convergence = 1e-8
        integer(ip) :: nmax = 200

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
                call conjugate_gradient(n, ed_vec, ipd0_d, TMatVec, PolVec)
                if(amoeba) call conjugate_gradient(n, ep_vec, ipd0_p, TMatVec, PolVec)

            case(OMMP_SOLVER_DIIS)
                call jacobi_diis(n,verbose,diis_max,norm_jacobi,convergence,ed_vec,ipd0_d,nmax,status,TMatVec_offdiag,PolVec)
                if(amoeba) call jacobi_diis(n,verbose,diis_max, &
                norm_jacobi,convergence,ep_vec,ipd0_p,nmax,status,TMatVec_offdiag,PolVec)

            case(OMMP_SOLVER_INVERSION)
                call induce_dipoles_MatInv(ed_vec,ipd0_d)  !TODO: read polarization matrix TMat and not calculate it separately
                if(amoeba) call induce_dipoles_MatInv(ep_vec,ipd0_p)
                
            case default
                call fatal_error("Unknown solver for calculation of the induced point dipoles") 
        end select
        
        ! Reshape dipole vector into the matrix 
        ipds(:,:,1) = reshape(ipd0_d, (/3_ip, pol_atoms/)) 
        if(amoeba) ipds(:,:,2) = reshape(ipd0_p, (/3_ip, pol_atoms/)) 

        call mfree('polarization [TMat]',TMat)

    end subroutine polarization
end module mod_polarization
