module mod_mmpol
    !! Main module for the control of openMMPol library. It contains
    !! all the scalar and vector (allocatable) quantities needed to
    !! build up the atomistic polarizable embedding model and perform
    !! the calculation required from the quantum chemical software.
    
    use mod_memory, only: ip, rp

    implicit none 
    !private
    
    ! scalar control variables

    integer(ip), parameter :: maxn12 = 8, maxn13 = 24, maxn14 = 72, &
                              maxn15 = 216
    !! maximum allowed number of 1-2, ..., 1-5 neighbors.

    integer(ip), parameter :: maxpgp = 120
    !! maximum number of members for the same polarization group
    
    integer(ip), parameter :: ws_max = 100
    !! maximum number of neighbors to consider for building 
    !! the Wang-Skeel preconditioner
    
    integer(ip), protected :: verbose = 0_ip
    !! verbosity flag, allowed range 0 (no printing at all) -- 
    !! 3 (debug printing)
    ! TODO adapt to this convention!!
    
    integer(ip), protected :: ff_type
    !! Force field type selection flag (0 for AMBER, 1 for AMOEBA)
    
    integer(ip), protected :: ff_rules
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
    !! convergence threshold (rms norm of the residual/increment) for 
    !! iterative solvers
  
    logical :: gradient
    !! ??
    
    real(rp) :: ws_shift
    !! diagonal shift parameter for the Wang-Skeel preconditioner 
    !! (default: 1.0 for AMBER, 2.0 for AMOEBA).
    
    logical, protected :: amoeba
    
    integer(ip), protected :: mm_atoms !! number of MM atoms
    integer(ip), protected :: pol_atoms !! number of polarizable atoms
    integer(ip), protected :: ld_cart, ld_cder
!!     size of the cartesian multipolar distribution (i.e., (l+1)*(l+2)*(l+3)/6)
!!     this is 1 for AMBER (charges only), 10 for AMOEBA (up to quadrupoles). 
!!     this is also the size of the array that contains the electrostatic properties
!!     of the sources at the sources. ld_cder is the leading size of the derivative of
!!     such a distribution, which is 3 for AMBER and 19 for AMOEBA.
    integer(ip), protected :: n_ipd 
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
    ! TODO those quantities could probably be used only at need and
    ! then removed.

    ! scalars and arrays for various useful intermediates and results
    
    real(rp) :: e_ele, e_pol, e_qd, e_dd
    !! electrostatic and polarization energies, including 
    !! their breakdown into contributoins

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

    subroutine set_verbosity(v)
        integer(ip), intent(in) :: v

        if( v < 0 ) then
            verbose = 0_ip
        else if( v > 3 ) then 
            verbose = 3_ip
        else
            verbose = v
        end if

    end subroutine set_verbosity

    subroutine mmpol_init(l_ff_type, l_ff_rules, l_mm_atoms, l_pol_atoms)
        !! Performs all the memory allocation and vector initialization
        !! needed to run the openMMPol library
        
        use mod_memory, only: ip, rp, mallocate

        implicit none

        integer(ip), intent(in) :: l_ff_type
        !! Force field type used in initialization
        
        integer(ip), intent(in) :: l_ff_rules
        !! Exclusion rules type used in initialization
        
        integer(ip), intent(in) :: l_mm_atoms
        !! Number of MM atoms used in initialization
        
        integer(ip), intent(in) :: l_pol_atoms
        !! Number of polarizable atoms used in initialization
        
        ! FF related settings
        ff_type = l_ff_type
        ff_rules = l_ff_rules
        mm_atoms = l_mm_atoms
        pol_atoms = l_pol_atoms

        if(ff_type == 1) then
            amoeba = .true.
            ld_cart = 10_ip
            ld_cder = 19_ip
            n_ipd = 2_ip
        else if(ff_type == 0) then
            amoeba = .false.
            ld_cart = 1_ip
            ld_cder = 3_ip
            n_ipd = 1_ip
        end if
  
        call set_screening_parameters
        
        ! Memory allocation
        call mallocate('mmpol_init [cmm]', 3_ip, mm_atoms, cmm)
        call mallocate('mmpol_init [q]', ld_cart, mm_atoms, q)
        call mallocate('mmpol_init [n12]', mm_atoms, n12)
        call mallocate('mmpol_init [i12]', maxn12, mm_atoms, i12)
        call mallocate('mmpol_init [n12]', mm_atoms, n13)
        call mallocate('mmpol_init [i13]', maxn13, mm_atoms, i13)
        call mallocate('mmpol_init [n14]', mm_atoms, n14)
        call mallocate('mmpol_init [i14]', maxn14, mm_atoms, i14)
        call mallocate('mmpol_init [group]', mm_atoms, group)
        call mallocate('mmpol_init [pol]', pol_atoms, pol)
        call mallocate('mmpol_init [cpol]', 3_ip, pol_atoms, cpol)
        call mallocate('mmpol_init [polar_mm]', pol_atoms, polar_mm)
        call mallocate('mmpol_init [mm_polar]', mm_atoms, mm_polar)
        call mallocate('mmpol_init [thole]', mm_atoms, thole)
        call mallocate('mmpol_init [idp]', 3_ip, pol_atoms, n_ipd, ipd) 
        ipd = 0.0_rp

        if (amoeba) then
            ! Extra quantities that should be allocated only
            ! for AMOEBA
            call mallocate('mmpol_init [q0]', ld_cart, mm_atoms, q0)
            call mallocate('mmpol_init [n15]', mm_atoms, n15)
            call mallocate('mmpol_init [i15]', maxn15, mm_atoms, i15)
            call mallocate('mmpol_init [np11]', mm_atoms, np11)
            call mallocate('mmpol_init [ip11]', maxpgp, mm_atoms, ip11)
            call mallocate('mmpol_init [np12]', mm_atoms, np12)
            call mallocate('mmpol_init [ip12]', maxpgp, mm_atoms, ip12)
            call mallocate('mmpol_init [np13]', mm_atoms, np13)
            call mallocate('mmpol_init [ip13]', maxpgp, mm_atoms, ip13)
            call mallocate('mmpol_init [np14]', mm_atoms, np14)
            call mallocate('mmpol_init [ip14]', maxpgp, mm_atoms, ip14)
            call mallocate('mmpol_init [mol_frame]', mm_atoms, mol_frame)
            call mallocate('mmpol_init [ix]', mm_atoms, ix)
            call mallocate('mmpol_init [iy]', mm_atoms, iy)
            call mallocate('mmpol_init [iz]', mm_atoms, iz)
        end if
  
        call mallocate('mmpol_init [v_qq]', ld_cart, mm_atoms, v_qq)
        v_qq = 0.0_rp
        call mallocate('mmpol_init [ef_qd]', 3_ip, pol_atoms, n_ipd, ef_qd)
        ef_qd = 0.0_rp
        call mallocate('mmpol_init [dv_qq]', ld_cder, mm_atoms, dv_qq)
        dv_qq = 0.0_rp
        call mallocate('mmpol_init [def_qd]', 3_ip, pol_atoms, n_ipd, def_qd)
        def_qd = 0.0_rp
        call mallocate('mmpol_init [v_dq]', ld_cart, mm_atoms, v_dq)
        v_dq = 0.0_rp
        call mallocate('mmpol_init [ef_dd]', 3_ip, pol_atoms, n_ipd, ef_dd)
        ef_dd = 0.0_rp
        call mallocate('mmpol_init [dv_dq]', ld_cder, mm_atoms, dv_dq)
        dv_qq = 0.0_rp
        call mallocate('mmpol_init [def_dd]', 6_ip, pol_atoms, n_ipd, def_dd)
        def_dd = 0.0_rp
      
    end subroutine mmpol_init

    subroutine mmpol_prepare()
        !! Compute some derived quantities from the input that 
        !! are used during the calculation. The upstream code have
        !! to provide cmm, ... TODO.
        !! This routine 
        !!   * compute connectivity lists from connected atoms
        !!   * invert polar_mm list creating mm_polar
        !!   * populate cpol list of coordinates
        !!   * compute factors for thole damping 
        !!   * scales by 1/3 AMOEBA quadrupoles (?)
        !!   * TODO pol groups?
        !!   * performs multipoles rotation
        
        implicit none

        integer(ip) :: i
        real(rp) :: xx(3) ! TODO remove this variable

        ! compute connectivity lists from connected atoms
        call connectivity_list_init()
        
        ! invert mm_polar list creating mm_polar
        mm_polar(:) = 0
        do i = 1, pol_atoms
            mm_polar(polar_mm(i)) = i
        end do

        ! populate cpol list of coordinates
        do i = 1, pol_atoms
            cpol(:,i) = cmm(:, polar_mm(i))
        end do

        ! compute factors for thole damping
        call thole_init()

        if(amoeba) then
            ! Copy multipoles from q to q0
            q0 = q

            ! scales by 1/3 AMOEBA quadrupoles (?)
            ! Mysterious division of multipoles by three
            ! FL told me that it was done like that in
            ! Tinker
            q0(5:10,:) = q0(5:10,:) / 3.0_rp

            ! pol groups?

            ! performs multipoles rotation
            call rotate_multipoles(.false.,xx,xx)
        end if

    end subroutine mmpol_prepare

    subroutine mmpol_terminate()
        !! Performs all the deallocation needed at the end of the 
        !! calculation
        use mod_memory, only: mfree

        call mfree('mmpol_terminate [cmm]', cmm)
        call mfree('mmpol_terminate [q]', q)
        call mfree('mmpol_terminate [n12]', n12)
        call mfree('mmpol_terminate [i12]', i12)
        call mfree('mmpol_terminate [n12]', n13)
        call mfree('mmpol_terminate [i13]', i13)
        call mfree('mmpol_terminate [n14]', n14)
        call mfree('mmpol_terminate [i14]', i14)
        call mfree('mmpol_terminate [group]', group)
        call mfree('mmpol_terminate [pol]', pol)
        call mfree('mmpol_terminate [cpol]', cpol)
        call mfree('mmpol_terminate [polar_mm]', polar_mm)
        call mfree('mmpol_terminate [mm_polar]', mm_polar)
        call mfree('mmpol_terminate [thole]', thole)
        call mfree('mmpol_terminate [idp]', ipd) 
        call mfree('mmpol_terminate [v_qq]', v_qq)
        call mfree('mmpol_terminate [ef_qd]', ef_qd)
        call mfree('mmpol_terminate [dv_qq]',  dv_qq)
        call mfree('mmpol_terminate [def_qd]', def_qd)
        call mfree('mmpol_terminate [v_dq]', v_dq)
        call mfree('mmpol_terminate [ef_dd]', ef_dd)
        call mfree('mmpol_terminate [dv_dq]', dv_dq)
        call mfree('mmpol_terminate [def_dd]', def_dd)

        if (amoeba) then
            ! Extra quantities that should be deallocated only
            ! for AMOEBA
            call mfree('mmpol_terminate [q0]', q0)
            call mfree('mmpol_terminate [n15]', n15)
            call mfree('mmpol_terminate [i15]', i15)
            call mfree('mmpol_terminate [np11]', np11)
            call mfree('mmpol_terminate [ip11]', ip11)
            call mfree('mmpol_terminate [np12]', np12)
            call mfree('mmpol_terminate [ip12]', ip12)
            call mfree('mmpol_terminate [np13]', np13)
            call mfree('mmpol_terminate [ip13]', ip13)
            call mfree('mmpol_terminate [np14]', np14)
            call mfree('mmpol_terminate [ip14]', ip14)
            call mfree('mmpol_terminate [mol_frame]', mol_frame)
            call mfree('mmpol_terminate [ix]', ix)
            call mfree('mmpol_terminate [iy]', iy)
            call mfree('mmpol_terminate [iz]', iz)
        end if

    end subroutine mmpol_terminate
    
    subroutine fatal_error(message)
        !! Prints a message and exit from the program. This
        !! function should be used in all the conditions 
        !! where the program cannot proceed.

        implicit none
      
        character (len=*) message
        !! Message to print before the program termination

        write(6, '(t3,a)') message
        stop '   error termination for open_mmpol.'
    end subroutine fatal_error

    subroutine thole_init()
        ! This routine compute the thole factors and stores
        ! them in a vector. TODO add reference
        use mod_constants, only: a_wal, a_wdl
        
        implicit none
        
        integer(ip) :: i, j
        
        thole = 0.0_rp
        
        do i = 1, pol_atoms
            j = polar_mm(i)
            thole(j) = pol(i) ** (1.0_rp/6.0_rp)
        end do
        
        if(.not. amoeba) then
            if(ff_rules == 0) &
                thole = thole * sqrt(a_wal)
            if(ff_rules == 1) &
                thole = thole * sqrt(a_wdl)
        end if
    end subroutine thole_init

    subroutine connectivity_list_init()
        !! Subroutine to initialize the connectivity lists.
        !! Basically it starts from a list of atoms that are
        !! connected () and it computes list of atoms that are
        !! n-bonds away. n=2,3,4
        
        implicit none

        integer(ip) :: i, j, k, l, ij, neigh
        logical :: bad_neigh

        n12 = 0
        n13 = 0
        i13 = 0
        n14 = 0
        i14 = 0
        
        ! n12:
        do i = 1, mm_atoms
            do j = 1, maxn12
                if(i12(j,i).ne.0) n12(i) = n12(i) + 1
            end do
            ! sort i12:
            call ihsort(.true.,n12(i),i12(:,i))
        end do

        ! n13, i13:
        do i = 1, mm_atoms
            do j = 1, n12(i)
                ij = i12(j,i)
                do k = 1, n12(ij)
                    neigh = i12(k,ij)
                    ! check that neigh is actually a 1-3 neighbor:
                    bad_neigh = .false.
                    if (neigh.eq.i) bad_neigh = .true.
                    do l = 1, n12(i)
                        if (neigh.eq.i12(l,i)) bad_neigh = .true.
                    end do
                    if (.not. bad_neigh) then
                        n13(i) = n13(i) + 1
                        i13(n13(i),i) = neigh
                    end if
                end do
            end do
            ! sort i13 and remove duplicates:
            call ihsort(.true.,n13(i),i13(:,i))
        end do
    
        ! n14, i14:
        do i = 1, mm_atoms
            do j = 1, n13(i)
                ij = i13(j,i)
                do k = 1, n12(ij)
                    neigh = i12(k,ij)
                    ! check that neigh is actually a 1-4 neighbor

                    bad_neigh = .false.
                    if(neigh.eq.i) bad_neigh = .true.
                    do l = 1, n12(i)
                        if(neigh.eq.i12(l,i)) bad_neigh = .true.
                    end do
                    do l = 1, n13(i)
                        if (neigh.eq.i13(l,i)) bad_neigh = .true.
                    end do
                    if (.not. bad_neigh) then
                        n14(i) = n14(i) + 1
                        i14(n14(i),i) = neigh
                    end if
                end do
            end do
            ! sort i14 and remove duplicates:
            call ihsort(.true.,n14(i),i14(:,i))
        end do
    end subroutine connectivity_list_init

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
