module mod_mmpol
    !! Main module for the control of openMMPol library. It contains
    !! all the scalar and vector (allocatable) quantities needed to
    !! build up the atomistic polarizable embedding model and perform
    !! the calculation required from the quantum chemical software.
    
    use mod_memory

    implicit none 
    !private
    
    ! scalar control variables

    integer(ip), parameter :: maxn12 = 8, maxn13 = 24, maxn14 = 72, maxn15 = 216
    !! maximum allowed number of 1-2, ..., 1-5 neighbors.

    integer(ip), parameter :: maxpgp = 120
    !! maximum number of members for the same polarization group
    
    integer(ip), parameter :: ws_max = 100
    !! maximum number of neighbors to consider for building the Wang-Skeel preconditioner
    
    integer(ip) :: verbose
    !! verbosity flag, allowed range 0 (no printing at all) -- 3 (debug printing)
    ! TODO adapt to this convention!!
    
    integer(ip) :: ff_type
    !! Force field type selection flag (0 for AMBER, 1 for AMOEBA)
    
    integer(ip) :: ff_rules
    !! Force field exclusion rules (0 for Wang AL, 1 for Wang DL)

    integer(ip) :: solver
    !! Parameter that control how the polarization equations are solved.
    !! 1 = preconditioned conjugate gradient (default),
    !! 2 = jacobi iterations with DIIS extrapolation,
    !! 3 = matrix inversion
    
    integer(ip) :: matrix_vector
    !! Selection flag for the matrix-vectror product routine
    !! 1 = assemble the matrix using O(n^2) storage and use dgemv
    !! 2 = compute the matrix vector products in a direct fashion (default)
    !! 3 = use a fast multiplication technique (TODO)
  
    integer(ip) :: nmax
    !! maximum number of steps for iterative solvers
    
    real(rp) :: convergence
    !! convergence threshold (rms norm of the residual/increment) for iterative solvers
  
    logical :: gradient
    !! ??
    
    real(rp) :: ws_shift
    !! diagonal shift parameter for the Wang-Skeel preconditioner (default: 1.0 for AMBER, 2.0 for AMOEBA).
    
    integer(ip) :: len_inname
    !! length of input file name
    
    character(len=120) :: input_file
    !! name of the input file
    
    logical :: amoeba
    
    integer(ip) :: mm_atoms !! number of MM atoms
    integer(ip) :: pol_atoms !! number of polarizable atoms
    integer(ip) :: ld_cart, ld_cder
!!     size of the cartesian multipolar distribution (i.e., (l+1)*(l+2)*(l+3)/6)
!!     this is 1 for AMBER (charges only), 10 for AMOEBA (up to quadrupoles). 
!!     this is also the size of the array that contains the electrostatic properties
!!     of the sources at the sources. ld_cder is the leading size of the derivative of
!!     such a distribution, which is 3 for AMBER and 19 for AMOEBA.
    integer(ip) :: n_ipd 
    !! number of induced point dipoles distributions 
    !! this is 1 for AMBER and 2 for AMOEBA
    
    !   arrays for the force field dependent exclusion factors. 
    
    ! TODO MB22 why those are vectors?
    real(rp) :: mscale(4)
    !! factors for charge-charge (or multipole-multipole) interactions

    real(rp) :: pscale(5)
    !! factors for chrage-ipd (or multipole-ipd) interactions.
    !! in AMOEBA, this is used to define the polarization field, i.e., the right-hand
    !! side to the polarization equations, and depends on the connectivity.
    
    real(rp) :: dscale(4)
    !! factors for multipoles-ipd interactions used to compute the direct field,
    !! which is used to define the polarization energy. these factors depend on 
    !! the polarization group "connectivity" (AMOEBA only)

    real(rp) :: uscale(4)
    !! factor for ipd-ipd interactions. these depend on the connectivity (AMBER)
    !! or on the polarization group " connectivity (AMOEBA)

    ! allocatable arrays which describe the polarizable system
    
    real(rp), allocatable, target :: cmm(:,:)
    !! Coordinates of MM atoms (3:mm_atoms)
    
    real(rp), allocatable, target :: cpol(:,:)
    !! Coordinates of polarizable atoms (3:pol_atoms)
    
    real(rp), allocatable, target :: q(:,:)
    !! Mutlipolar distribution (ld_cart:mm_atoms)
    !! For AMOEBA this is the rotated distribution.
    !! The order for the stored multipoles is
    !! q, px, py, pz, Qxx, Qxy, Qyy, Qxz, Qyx, Qzz.

    real(rp), allocatable, target :: q0(:,:)
    !! Unrotated utlipolar distribution (ld_cart:mm_atoms)
    !! (AMOEBA only)
    
    real(rp), allocatable, target :: ipd(:,:,:)
    !! induced point dipoles (3:pol_atoms:ipd) 
    
    real(rp), allocatable :: pol(:)
    !! Polarizabilities for each polarizable atom
    
    integer(ip), allocatable :: mm_polar(:)
    !! indices of the MM atoms that are polarizable

    integer(ip), allocatable :: polar_mm(:)
    !! positions of a polarizable atom in the mm atoms list
    
    integer(ip), allocatable :: n12(:), i12(:,:), n13(:), i13(:,:), n14(:), i14(:,:), n15(:), i15(:,:)
    !! connectivity (number of neighbors and list of neighbors)
    
    integer(ip), allocatable :: group(:)
    !! polarization group or fragment
    
    integer(ip), allocatable :: np11(:), ip11(:,:), np12(:), ip12(:,:), np13(:), ip13(:,:), np14(:), ip14(:,:)
    !! polarization group "connectivity":

    ! parameters for the definition of the rotation matrices for the multipoles:
    
    integer(ip), allocatable :: mol_frame(:)
    !! definition of the molecular frame
    !! convention: 0 ... do not rotate
    !!             1 ... z-then-x
    !!             2 ... bisector
    !!             3 ... z-only
    !!             4 ... z-bisector
    !!             5 ... 3-fold

    integer(ip), allocatable :: ix(:), iy(:), iz(:)
    !! neighboring atoms used to define the axes of the molecular frame

    ! scalars and arrays for various useful intermediates and results
    
    real(rp) :: e_ele, e_pol, e_qd, e_dd
    !! electrostatic and polarization energies, including their breakdown into contributoins

!
!     potential (and higher order terms) of the multipoles 
!      at the charges (multipoles) and its derivatives 
!
    real(rp), allocatable :: v_qq(:,:), dv_qq(:,:)
    !! potential (and higher order terms) of the multipoles 
    !! at the charges (multipoles) and its derivatives
  
    real(rp), allocatable :: ef_qd(:,:,:), def_qd(:,:,:)
    !! field of the charges (multipoles) at the ipd and its derivatives
  
    real(rp), allocatable :: v_dq(:,:), dv_dq(:,:)
    !! potential (and higher order terms) of the induced point 
    !! dipoles at the charges (multipoles) and its derivatives
    
    real(rp), allocatable :: ef_dd(:,:,:), def_dd(:,:,:)
    !! field of the ipd at the ipd and its derivatives
    
    integer(ip), allocatable :: n_ws(:), list_ws(:)
    !! things related to the wang-skeel preconditioner
    
    real(rp),    allocatable :: block_ws(:,:)
    !! things related to the wang-skeel preconditioner
    
    real(rp),    allocatable :: thole(:)
    !! array to store the thole factors for computing damping functions
    
    contains
    
    subroutine fatal_error(message)
        !! Print a message and exit from the program. This
        !! function should be used in all the conditions 
        !! where the program cannot proceed.

        implicit none
      
        character (len=*) message
        !! Message to print before the program termination

        write(6, '(t3,a)') message
        stop '   error termination for open_mmpol.'
    end subroutine fatal_error

    subroutine set_screening_parameters()
        !! Subroutine to initialize the screening parameters
        use mod_constants, only: one, zero, pt5
        
        implicit none
        real(rp), parameter :: pt4 = 0.40_rp, pt8 = 0.80_rp
        if (ff_type.eq.0 .and. ff_rules.eq.0) then
            ! WangAL
            mscale(1) = zero
            mscale(2) = zero
            mscale(3) = one
            mscale(4) = one
            pscale(1) = zero
            pscale(2) = zero
            pscale(3) = one
            pscale(4) = one
            pscale(5) = one
            dscale(1) = zero
            dscale(2) = zero
            dscale(3) = one
            dscale(4) = one
            uscale(1) = zero
            uscale(2) = zero
            uscale(3) = one
            uscale(4) = one
        else if (ff_type.eq.0 .and. ff_rules.eq.1) then
            ! WangDL
            mscale(1) = one
            mscale(2) = one
            mscale(3) = one
            mscale(4) = one
            pscale(1) = one
            pscale(2) = one
            pscale(3) = one
            pscale(4) = one
            pscale(5) = one
            dscale(1) = one
            dscale(2) = one
            dscale(3) = one
            dscale(4) = one
            uscale(1) = one
            uscale(2) = one
            uscale(3) = one
            uscale(4) = one
        else if (ff_type.eq.1) then
            ! AMOEBA
            mscale(1) = zero
            mscale(2) = zero
            mscale(3) = pt4
            mscale(4) = pt8
            pscale(1) = zero
            pscale(2) = zero
            pscale(3) = one
            pscale(4) = one
            pscale(5) = pt5
            dscale(1) = zero
            dscale(2) = one
            dscale(3) = one
            dscale(4) = one
            uscale(1) = one
            uscale(2) = one
            uscale(3) = one
            uscale(4) = one
        else
            call fatal_error('the required force field is not implemented.')
        end if
    end subroutine set_screening_parameters

end module mod_mmpol
