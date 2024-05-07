#include "f_cart_components.h"
module mod_electrostatics
    use mod_io, only: fatal_error, ommp_message
    use mod_constants, only: OMMP_VERBOSE_DEBUG
    use mod_profiling, only: time_push, time_pull
    use mod_memory, only: ip, rp, lp
    use mod_adjacency_mat, only: yale_sparse
    use mod_topology, only: ommp_topology_type
    use mod_profiling
    use fmmlib_interface

    !! TODO Check the signs in electrostatic elemental functions
    !! TODO [OPT] Use Laplace equation to simplify the calculations:
    !! TODO  1. Egrd(_zz_) = -(Egrd(_xx_) + Egrd(_yy_))
    !! TODO  2. EHes(_zzz_) = -(EHes(_xxz_) + EHes(_yyz_))
    !! TODO  3. EHes(_zzx_) = -(EHes(_xxx_) + EHes(_yyx_))
    !! TODO  4. EHes(_zzy_) = -(EHes(_xxy_) + EHes(_yyy_))
    !! TODO [OPT] Fundamental electrostatic functions should be pure/elemental
    !! TODO [BUG] Handling of flags gg
    implicit none 
    private

    type ommp_electrostatics_type
        integer(ip) :: def_solver
        !! Solver to be used by default for this eel object.
        integer(ip) :: def_matv
        !! Matrix-vector method to be used by default for this eel object.

        type(ommp_topology_type), pointer :: top
        !! Data structure containing all the topological informations
        integer(ip) :: pol_atoms 
        !! number of polarizable atoms
        logical(lp) :: amoeba
        !! True if AMOEBA FF is used
        integer(ip) :: ld_cart, ld_cder
        !! size of the cartesian multipolar distribution (i.e., (l+1)*(l+2)*(l+3)/6)
        !! this is 1 for AMBER (charges only), 10 for AMOEBA (up to quadrupoles). 
        !! this is also the size of the array that contains the electrostatic properties
        !! of the sources at the sources. ld_cder is the leading size of the derivative of
        !! such a distribution, which is 3 for AMBER and 19 for AMOEBA.
        integer(ip) :: n_ipd 
        !! number of induced point dipoles distributions 
        !! this is 1 for AMBER and 2 for AMOEBA
        
        real(rp) :: mscale(4)
        !! factors for charge-charge (or multipole-multipole) interactions

        real(rp) :: pscale(4)
        !! factors for chrage-ipd (or multipole-ipd) interactions.
        !! in AMOEBA, this is used to define the polarization field, i.e., the right-hand
        !! side to the polarization equations, and depends on the connectivity.

        real(rp) :: pscale_intra(4)
        !! Only used for AMOEBA, same as pscale but for atoms that belong to the 
        !! same polarization group
        
        real(rp) :: dscale(4)
        !! factors for multipoles-ipd interactions used to compute the direct field,
        !! which is used to define the polarization energy. these factors depend on 
        !! the polarization group "connectivity" (AMOEBA only)

        real(rp) :: uscale(4)
        !! factor for ipd-ipd interactions. these depend on the connectivity (AMBER)
        !! or on the polarization group " connectivity (AMOEBA)

        ! allocatable arrays which describe the polarizable system
        
        real(rp), allocatable :: thole(:)
        !! array to store the thole factors for computing damping functions
        
        real(rp) :: thole_scale
        !! Scale factor for thole damping (only used by non-AMOEBA FF); all
        !! the element of thole(:) are multiplied by thole_scale ** 0.5

        real(rp), allocatable :: cpol(:,:)
        !! Coordinates of polarizable atoms (3:pol_atoms)
        
        real(rp), allocatable :: q(:,:)
        !! Mutlipolar distribution (ld_cart:mm_atoms)
        !! For AMOEBA this is the rotated distribution.
        !! The order for the stored multipoles is
        !! q, px, py, pz, Qxx, Qxy, Qyy, Qxz, Qyx, Qzz.

        real(rp), allocatable :: q0(:,:)
        !! Unrotated utlipolar distribution (ld_cart:mm_atoms)
        !! (AMOEBA only)
        
        real(rp), allocatable :: pol(:)
        !! Polarizabilities for each polarizable atom
        
        integer(ip), allocatable :: mm_polar(:)
        !! for each mm atom: 0 if is not polarizable, index in 
        !! polarizable atom list otherwise

        integer(ip), allocatable :: polar_mm(:)
        !! positions of a polarizable atom in the mm atoms list

        integer(ip), allocatable :: C_polar_mm(:)
        !! [[polar_mm]] with 0-based C-indexing, only allocated at need.
        
        integer(ip), allocatable :: mmat_polgrp(:)
        !! Polarizability group index for each MM site

        type(yale_sparse) :: polgrp_mmat
        !! For each polarization group index, list all the MM atoms included.
        !! It basically is a sparse boolean matrix of dimension 
        !! N_polgroups x N_mmatoms

        type(yale_sparse), allocatable :: pg_conn(:)
        !! Adjacency and connectivity matytrices between polarizability groups.
        !! Two groups are said to be adjacent if they are connected by a chemical 
        !! bond. The 1st element is the identity matrix for code simplicity.
        
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

        !- FMM quantities here
        logical(lp) :: use_fmm = .false.
        !! flag to use fast multipoles
        integer(ip) :: fmm_maxl = 0
        !! Maximum angular moment used in fast multipoles
        type(fmm_type), allocatable :: fmm_static
        !! Fast multipoles object for static multipoles sources
        logical(lp) :: fmm_static_done = .false.
        !! Flag for a fresh solution of ipd fmm
        type(fmm_type), allocatable :: fmm_ipd(:)
        !! Fast multipoles object for static multipoles sources
        logical(lp), allocatable :: fmm_ipd_done(:)
        !! Flag for a fresh solution of ipd fmm
        type(fmm_tree_type), allocatable :: tree
        !! Tree object
        type(yale_sparse) :: fmm_near_field_list
        !! For each particle, all the particles that should be included in near field
    
        !- Intermediate data allocate here -!
        logical(lp) :: M2M_done = .false.
        !! flag to set when M2M electrostatic quantities are computed.
        logical(lp) :: M2Mgg_done = .false.
        !! flag to set when M2M electrostatic quantities for geometrical
        !! gradients are computed.
        real(rp), allocatable :: V_M2M(:)
        !! potential of MM permanent multipoles at MM sites; 
        real(rp), allocatable :: E_M2M(:,:)
        !! electric_field of MM permanent multipoles at MM sites; 
        real(rp), allocatable :: Egrd_M2M(:,:)
        !! electric_field gradient of MM permanent multipoles at MM sites;
        real(rp), allocatable :: EHes_M2M(:,:) 
        !! electric field Hessian of MM permanent multipoles at MM sites;

        logical(lp) :: M2D_done = .false.
        !! Flag to set when M2D electrostatics have been computed.
        logical(lp) :: M2Dgg_done = .false.
        !! Flag to set when M2D electrostatics for geometrical gradients 
        !! have been computed.
        real(rp), allocatable :: V_M2D(:,:)
        ! electrostatic potential of MM permanent multipoles at POL sites; unused.
        real(rp), allocatable :: E_M2D(:,:,:) ! TODO the third dimension is used?
        !! electric field of MM permanent multipoles at POL sites; 
        real(rp), allocatable :: Egrd_M2D(:,:,:)
        !! electric field of MM permanent multipoles at POL sites;
        real(rp), allocatable :: EHes_M2D(:,:,:)
        ! electric field Hessian of MM permanent multipoles at POL sites; unused.

        logical(lp) :: D2Mgg_done = .false.
        real(rp), allocatable :: V_D2M(:)
        real(rp), allocatable :: E_D2M(:,:)
        real(rp), allocatable :: Egrd_D2M(:,:)
        real(rp), allocatable :: EHes_D2M(:,:)

        logical(lp) :: D2Dgg_done = .false.
        real(rp), allocatable :: V_D2D(:,:)
        real(rp), allocatable :: E_D2D(:,:,:)
        real(rp), allocatable :: Egrd_D2D(:,:,:)
        real(rp), allocatable :: EHes_D2D(:,:,:)
        
        logical(lp) :: ipd_done = .false.
        !! Flag to set when IPD have been computed.
        logical(lp) :: ipd_use_guess = .false.
        !! Flag to set when current value of IPD can be
        !! used as guess for next solution of LS.
        real(rp), allocatable :: ipd(:,:,:)
        !! induced point dipoles (3:pol_atoms:ipd) 
    
        real(rp), allocatable :: TMat(:,:)
        !! Interaction tensor, only allocated for the methods that explicitly 
        !! requires it.
        
        logical(lp) :: screening_list_done = .false.
        !! Flag to check if screening list have already been prepared
        type(yale_sparse), allocatable :: list_S_S, list_P_P, list_S_P_P, list_S_P_D                                                             
        !! Sparse matrix containg the scale factors for the scaled elements                                                                      
        !! of electrostatic interactions (all the element that are not                                                                           
        !! present in the sparse matrix have a scaling factor 1.0).                                                                              
        !! When FMM are enabled, those lists are only for near-field                                                                             
        !! interactions                                                                                                                          
        type(yale_sparse), allocatable :: list_S_S_fmm_far, &                                                                                    
                                          list_P_P_fmm_far, &                                                                                    
                                          list_S_P_P_fmm_far, &                                                                                  
                                          list_S_P_D_fmm_far                                                                                     
        !! Sparse matrices containing the scale factors for the scaled elements                                                                  
        !! of electrostatic far-field interactions (only when FMM are                                                                            
        !! enabled). If a reasonable FMM distance is chosen, those lists                                                                         
        !! should be almost empty.     
        logical(lp), dimension(:), allocatable :: todo_S_S, todo_P_P, todo_S_P_P, todo_S_P_D
        !! Logical array of the same dimension of column-index vector; true if 
        !! the scaling factor is zero, false otherwise
        logical(lp), dimension(:), allocatable :: todo_S_S_fmm_far, todo_P_P_fmm_far, &
                                                  todo_S_P_P_fmm_far, todo_S_P_D_fmm_far
        !! Logical array of the same dimension of column-index vector; true if 
        !! the scaling factor is zero, false otherwise
        real(rp), dimension(:), allocatable :: scalef_S_S_fmm_far, scalef_P_P_fmm_far, &
                                               scalef_S_P_P_fmm_far, scalef_S_P_D_fmm_far
        !! Array of the same dimension of column-index vector; contains the 
        !! value of the scaling factors different from 1.0
        real(rp), dimension(:), allocatable :: scalef_S_S, scalef_P_P, &
                                               scalef_S_P_P, scalef_S_P_D
        !! Array of the same dimension of column-index vector; contains the 
        !! value of the scaling factors different from 1.0
    end type ommp_electrostatics_type

    public :: ommp_electrostatics_type
    public :: electrostatics_init, electrostatics_terminate
    public :: set_def_solver, set_def_matv
    public :: thole_init, remove_null_pol, set_screening_parameters
    public :: screening_rules, make_screening_lists
    public :: damped_coulomb_kernel, field_extD2D
    public :: energy_MM_MM, energy_MM_pol
    public :: prepare_fixedelec, prepare_polelec
    public :: q_elec_prop, coulomb_kernel
    public :: potential_M2E, potential_D2E
    public :: field_M2E, field_D2E
    public :: fmm_coordinates_update

    contains

    subroutine electrostatics_init(eel_obj, amoeba, pol_atoms, top_obj)
        use mod_memory, only: mallocate
        use mod_constants, only: OMMP_MATV_DEFAULT, OMMP_SOLVER_DEFAULT, &
                                 OMMP_FMM_DEFAULT_MAXL, OMMP_FMM_ENABLE_THR

        implicit none 

        logical(lp), intent(in) :: amoeba
        integer(ip), intent(in) :: pol_atoms
        type(ommp_topology_type), intent(in), target :: top_obj
        type(ommp_electrostatics_type), intent(inout) :: eel_obj

        integer(ip) :: mm_atoms

        mm_atoms = top_obj%mm_atoms
        eel_obj%amoeba = amoeba
        eel_obj%pol_atoms = pol_atoms
        eel_obj%top => top_obj
        eel_obj%def_solver = OMMP_SOLVER_DEFAULT
        eel_obj%def_matv = OMMP_MATV_DEFAULT

        if(amoeba) then
            eel_obj%ld_cart = 10_ip
            eel_obj%ld_cder = 19_ip
            eel_obj%n_ipd = 2_ip
        else
            eel_obj%ld_cart = 1_ip
            eel_obj%ld_cder = 3_ip
            eel_obj%n_ipd = 1_ip
        endif

        if(mm_atoms > OMMP_FMM_ENABLE_THR) then
            eel_obj%use_fmm = .true.
            eel_obj%fmm_maxl = OMMP_FMM_DEFAULT_MAXL
            allocate(eel_obj%tree)
            allocate(eel_obj%fmm_static)
            eel_obj%fmm_static_done = .false.
            allocate(eel_obj%fmm_ipd(eel_obj%n_ipd))
            allocate(eel_obj%fmm_ipd_done(eel_obj%n_ipd))
            eel_obj%fmm_ipd_done = .false.
        end if

        call mallocate('electrostatics_init [q]', eel_obj%ld_cart, &
                       mm_atoms, eel_obj%q)
        call mallocate('electrostatics_init [pol]', eel_obj%pol_atoms, &
                       eel_obj%pol)
        call mallocate('electrostatics_init [cpol]', 3_ip, eel_obj%pol_atoms, &
                       eel_obj%cpol)
        call mallocate('electrostatics_init [polar_mm]', eel_obj%pol_atoms, &
                       eel_obj%polar_mm)
        call mallocate('electrostatics_init [mm_polar]', mm_atoms, &
                       eel_obj%mm_polar)
        call mallocate('electrostatics_init [thole]', mm_atoms, &
                       eel_obj%thole)
        
        call mallocate('electrostatics_init [idp]', 3_ip, eel_obj%pol_atoms, &
                       eel_obj%n_ipd, eel_obj%ipd) 
        eel_obj%ipd_done = .false.
        eel_obj%ipd = 0.0_rp

        if (eel_obj%amoeba) then
            ! Extra quantities that should be allocated only
            ! for AMOEBA
            call mallocate('electrostatics_init [q0]', eel_obj%ld_cart, &
                           mm_atoms, eel_obj%q0)
            
            call mallocate('electrostatics_init [mmat_polgrp]', & 
                           mm_atoms, eel_obj%mmat_polgrp)

            call mallocate('electrostatics_init [mol_frame]', &
                           mm_atoms, eel_obj%mol_frame)
            call mallocate('electrostatics_init [ix]', mm_atoms, eel_obj%ix)
            call mallocate('electrostatics_init [iy]', mm_atoms, eel_obj%iy)
            call mallocate('electrostatics_init [iz]', mm_atoms, eel_obj%iz)
        end if

    end subroutine electrostatics_init

    subroutine electrostatics_terminate(eel_obj)
        use mod_memory, only: mfree
        use mod_adjacency_mat, only: matfree

        implicit none

        type(ommp_electrostatics_type), intent(inout) :: eel_obj
        integer(ip) :: i

        call mfree('electrostatics_terminate [q]', eel_obj%q)
        call mfree('electrostatics_terminate [pol]', eel_obj%pol)
        call mfree('electrostatics_terminate [cpol]', eel_obj%cpol)
        call mfree('electrostatics_terminate [polar_mm]', eel_obj%polar_mm)
        call mfree('electrostatics_terminate [mm_polar]', eel_obj%mm_polar)
        call mfree('electrostatics_terminate [thole]', eel_obj%thole)
        call mfree('electrostatics_terminate [idp]', eel_obj%ipd) 

        if (eel_obj%amoeba) then
            call mfree('electrostatics_terminate [q0]', eel_obj%q0)
            call mfree('electrostatics_terminate [mmat_polgrp]', &
                       eel_obj%mmat_polgrp)
            if(allocated(eel_obj%pg_conn)) then
                do i=1, size(eel_obj%pg_conn)
                    call matfree(eel_obj%pg_conn(i))
                end do
                deallocate(eel_obj%pg_conn)
            end if
            
            call matfree(eel_obj%polgrp_mmat)
            call mfree('electrostatics_terminate [mol_frame]', &
                       eel_obj%mol_frame)
            call mfree('electrostatics_terminate [ix]', eel_obj%ix)
            call mfree('electrostatics_terminate [iy]', eel_obj%iy)
            call mfree('electrostatics_terminate [iz]', eel_obj%iz)
        end if
        
        call mfree('electrostatics_terminate [E_M2D]', eel_obj%E_M2D)
        call mfree('electrostatics_terminate [V_M2M]', eel_obj%V_M2M)
        call mfree('electrostatics_terminate [E_M2M]', eel_obj%E_M2M)
        call mfree('electrostatics_terminate [Egrd_M2M]', eel_obj%Egrd_M2M)

        if(allocated(eel_obj%todo_S_S)) deallocate(eel_obj%todo_S_S)
        if(allocated(eel_obj%todo_P_P)) deallocate(eel_obj%todo_P_P)
        call mfree('electrostatics_terminate [scalef_S_S]', eel_obj%scalef_S_S)
        call mfree('electrostatics_terminate [scalef_P_P]', eel_obj%scalef_P_P)
        if(allocated(eel_obj%list_S_S)) then
            call matfree(eel_obj%list_S_S)
            deallocate(eel_obj%list_S_S)
        end if
        if(allocated(eel_obj%list_P_P)) then
            call matfree(eel_obj%list_P_P)
            deallocate(eel_obj%list_P_P)
        end if

        if(eel_obj%use_fmm) then
            call free_fmm(eel_obj%fmm_static)
            do i=1, eel_obj%n_ipd
                call free_fmm(eel_obj%fmm_ipd(i))
            end do
            call free_tree(eel_obj%tree)
            deallocate(eel_obj%fmm_static, eel_obj%fmm_ipd, eel_obj%tree)
            call matfree(eel_obj%fmm_near_field_list)
        end if

    end subroutine electrostatics_terminate

    subroutine set_def_solver(eel_obj, solver)
        use mod_constants, only: OMMP_SOLVER_CG, OMMP_SOLVER_INVERSION, OMMP_SOLVER_DIIS
        implicit none
        
        type(ommp_electrostatics_type), intent(inout) :: eel_obj
        integer(ip), intent(in) :: solver

        if(solver /= OMMP_SOLVER_CG .and. &
           solver /= OMMP_SOLVER_INVERSION .and. &
           solver /= OMMP_SOLVER_DIIS) &
            call fatal_error("Unrecognized setting for default solver method")
        eel_obj%def_solver = solver
    end subroutine

    subroutine set_def_matv(eel_obj, matv)
        use mod_constants, only: OMMP_MATV_INCORE, OMMP_MATV_DIRECT
        implicit none
        
        type(ommp_electrostatics_type), intent(inout) :: eel_obj
        integer(ip), intent(in) :: matv

        if(matv /= OMMP_MATV_INCORE .and. &
           matv /= OMMP_MATV_DIRECT) &
            call fatal_error("Unrecognized setting for default matrix-vector method")
        eel_obj%def_matv = matv
    end subroutine
    
    subroutine set_screening_parameters(eel_obj, m, p, d, u, i)
        !! Subroutine to initialize the screening parameters
       
        implicit none

        type(ommp_electrostatics_type), intent(inout) :: eel_obj
        real(rp), intent(in) :: m(4), p(4), d(4), u(4)
        real(rp), optional, intent(in) :: i(4)
        
        eel_obj%mscale = m
        eel_obj%pscale = p
        eel_obj%dscale = d
        eel_obj%uscale = u
        
        if(present(i)) then
            if(eel_obj%amoeba) then
                eel_obj%pscale_intra = i
            else
                call fatal_error("Scale factors for atoms of the same group &
                                 &cannot be set outside AMOEBA FF")
            end if
        else
            if(eel_obj%amoeba) &
                call fatal_error("Scale factors for atoms of the same group &
                                 &should be defined in AMOEBA FF")
        end if
        
    end subroutine set_screening_parameters

    subroutine remove_null_pol(eel)
        !! Check which polarizabilities are close enough to 0 to be 
        !! just excluded from the calculation, and remove them.

        use mod_constants, only: OMMP_STR_CHAR_MAX, OMMP_VERBOSE_LOW, eps_rp
        use mod_memory, only: mallocate, mfree
        implicit none

        type(ommp_electrostatics_type), intent(inout) :: eel
        integer(ip), allocatable :: idx(:), polar_mm(:)
        integer(ip) :: i, nidx
        real(rp), allocatable :: tmp(:)
        character(len=OMMP_STR_CHAR_MAX) :: msg

        if(eel%pol_atoms > 0 .and. allocated(eel%pol)) then
            call mallocate('remove_null_pol [idx]', eel%pol_atoms, idx)
            nidx = 0
            do i=1, eel%pol_atoms
                if(abs(eel%pol(i)) > eps_rp) then
                    nidx = nidx + 1
                    idx(nidx) = i
                end if
            end do

            if(nidx < eel%pol_atoms) then
                write(msg, '(A,I0,A,I0)') "Removing ", eel%pol_atoms - nidx, &
                                          " polarizabilities out of ", eel%pol_atoms
                call ommp_message(msg, OMMP_VERBOSE_LOW)
                call mallocate('remove_null_pol [tmp]', eel%pol_atoms, tmp)
                ! Polarizabilities
                tmp = eel%pol
                call mfree('remove_null_pol [eel%pol]', eel%pol)
                call mallocate('remove_null_pol [eel%pol]', nidx, eel%pol)
                do i=1, nidx
                    eel%pol(i) = tmp(idx(i))
                end do
                call mfree('remove_null_pol [tmp]', tmp)
                    
                call mallocate('remove_null_pol [polar_mm]', eel%pol_atoms, polar_mm)
                polar_mm = eel%polar_mm
                call mfree('remove_null_pol [eel%polar_mm]', eel%polar_mm)
                call mallocate('remove_null_pol [eel%polar_mm]', nidx, eel%polar_mm)
                call mfree('remove_null_pol [eel%cpol]', eel%cpol)
                call mallocate('remove_null_pol [eel%cpol]', 3, nidx, eel%cpol)
                eel%mm_polar = 0
                do i=1, nidx
                    eel%polar_mm(i) = polar_mm(idx(i))
                    eel%mm_polar(eel%polar_mm(i)) = i
                    eel%cpol(:,i) = eel%top%cmm(:,eel%polar_mm(i))
                end do

                call mfree('remove_null_pol [polar_mm]', polar_mm)
                
                eel%pol_atoms = nidx
                call mfree('remove_null_pol [eel%ipd]', eel%ipd)
                call mallocate('electrostatics_init [idp]', 3_ip, eel%pol_atoms, &
                               eel%n_ipd, eel%ipd) 
                eel%ipd_done = .false.
                eel%ipd_use_guess = .false.
                eel%ipd = 0.0_rp
                
                call ommp_message("Removing null polarizable sites done", OMMP_VERBOSE_LOW)
            end if
            
            call mfree('remove_null_pol [idx]', idx)
        end if

    end subroutine

    subroutine make_screening_lists(eel)
        use mod_memory, only: mallocate, mfree
        use mod_adjacency_mat, only: compress_list, compress_data, matfree
        use mod_constants, only: eps_rp
        
        implicit none

        type(ommp_electrostatics_type), intent(inout) :: eel
        
        integer(ip) :: i, ineigh, ij, j, ns, ns_guess, n, npol, &
                       ipp, jp, ns_guess_grp, pg_i, igrp, grp, ib, ie
        logical :: to_do, to_scale, usenl
        real(rp) :: scalf
        integer(ip), allocatable :: itmp(:,:), ntmp(:), itmp_far(:,:), ntmp_far(:)
        real(rp), allocatable :: rtmp(:,:), rtmp_far(:,:)

        if(eel%screening_list_done) return

        n = eel%top%mm_atoms
        npol = eel%pol_atoms

        ! First estimate the maximum number of elements to scale
        ns_guess = 0
        !$omp parallel do default(shared) private(ns, ineigh, i)
        do i=1, n
            ns = 0
            do ineigh=1, 4
                ns = ns + eel%top%conn(ineigh)%ri(i+1) - eel%top%conn(ineigh)%ri(i)
            end do
            !$omp atomic
            ns_guess = max(ns, ns_guess)
        end do
        
        if(eel%amoeba) then
            ! If amoeba is used, some scaling (D-field) are performed using 
            ! polarization group, so a second guess should be done and the 
            ! maximum is used in memory allocation.
            ns_guess_grp = 0
            
            !$omp parallel do default(shared) & 
            !$omp private(i, pg_i, ineigh, igrp, grp, ns)
            do i=1, n
                pg_i = eel%mmat_polgrp(i)
                ns = 0
                do ineigh=1, 4
                    do igrp=eel%pg_conn(ineigh)%ri(pg_i), &
                            eel%pg_conn(ineigh)%ri(pg_i+1)-1

                        grp = eel%pg_conn(ineigh)%ci(igrp)
                        ns = ns + eel%polgrp_mmat%ri(grp+1)-eel%polgrp_mmat%ri(grp)
                    end do
                end do
                !$omp atomic
                ns_guess_grp = max(ns, ns_guess_grp)
            end do
            
            ns_guess = max(ns_guess_grp, ns_guess)
        end if

        allocate(eel%list_S_S)
        allocate(eel%list_P_P)
        allocate(eel%list_S_P_P)
        if(eel%amoeba) allocate(eel%list_S_P_D)

        call mallocate('make_screening_list [rtmp]', ns_guess, n, rtmp)
        call mallocate('make_screening_list [itmp]', ns_guess, n, itmp)
        call mallocate('make_screening_list [ntmp]', n, ntmp)

        if(eel%use_fmm) then
            allocate(eel%list_S_S_fmm_far)
            allocate(eel%list_P_P_fmm_far)
            allocate(eel%list_S_P_P_fmm_far)
            if(eel%amoeba) allocate(eel%list_S_P_D_fmm_far)
        
            call mallocate('make_screening_list [rtmp_far]', ns_guess, n, rtmp_far)
            call mallocate('make_screening_list [itmp_far]', ns_guess, n, itmp_far)
            call mallocate('make_screening_list [ntmp_far]', n, ntmp_far)

            if(.not. allocated(eel%fmm_near_field_list%ri)) call fmm_make_neigh_list(eel) 
        end if


        ! Build S S list
        !$omp parallel do schedule(dynamic) default(shared) &
        !$omp private(i,ineigh,ij,j, scalf,usenl)
        do i=1, n
            ntmp(i) = 0
            if(eel%use_fmm) ntmp_far(i) = 0
            do ineigh=1, 4
                do ij=eel%top%conn(ineigh)%ri(i), eel%top%conn(ineigh)%ri(i+1)-1
                    j = eel%top%conn(ineigh)%ci(ij)
                    scalf = screening_rules(eel, i, 'S', j, 'S', '-')
                    if(abs(scalf-1.0) > eps_rp) then
                        usenl = .true.
                        if(eel%use_fmm) usenl = fmm_list_are_near(eel,i,j)
                        if(usenl) then
                            ntmp(i) = ntmp(i) + 1
                            itmp(ntmp(i),i) = j
                            rtmp(ntmp(i),i) = scalf
                        else
                            ntmp_far(i) = ntmp_far(i) + 1
                            itmp_far(ntmp_far(i),i) = j
                            rtmp_far(ntmp_far(i),i) = scalf
                        end if
                    end if
                end do
            end do
        end do

        ! Compress the list
        call compress_list(n, itmp, ntmp, eel%list_S_S)
        call compress_data(eel%list_S_S, rtmp, eel%scalef_S_S)
        call mallocate('make_screening_list [todo_S_S]', &
                       size(eel%scalef_S_S), eel%todo_S_S)
        !$omp parallel do
        do i=1, size(eel%scalef_S_S)
            eel%todo_S_S(i) = (abs(eel%scalef_S_S(i)) > eps_rp)
        end do

        if(eel%use_fmm) then
            if(sum(ntmp_far) > 0) then
                ! There are elements to consider
                call compress_list(n, itmp_far, ntmp_far, eel%list_S_S_fmm_far)
                call compress_data(eel%list_S_S_fmm_far, rtmp_far, eel%scalef_S_S_fmm_far)
                call mallocate('make_screening_list [todo_S_S_fmm_far]', &
                               size(eel%scalef_S_S_fmm_far), eel%todo_S_S_fmm_far)
                !$omp parallel do
                do i=1, size(eel%scalef_S_S_fmm_far)
                    eel%todo_S_S_fmm_far(i) = (abs(eel%scalef_S_S_fmm_far(i)) > eps_rp)
                end do
            else
                ! The list is empty
                call matfree(eel%list_S_S_fmm_far)
                deallocate(eel%list_S_S_fmm_far)
            end if
        end if

        ! Build P P list
        ntmp = 0
        if(eel%use_fmm) ntmp_far = 0
        !$omp parallel do schedule(dynamic) default(shared) & 
        !$omp private(ipp,i,ineigh,ij,j,jp,scalf,usenl)
        do ipp=1, npol
            i = eel%polar_mm(ipp)
            do ineigh=1, 4
                do ij=eel%top%conn(ineigh)%ri(i), eel%top%conn(ineigh)%ri(i+1)-1
                    j = eel%top%conn(ineigh)%ci(ij)
                    jp = eel%mm_polar(j)
                    if(jp > 0) then
                        scalf = screening_rules(eel, ipp, 'P', jp, 'P', '-')
                        if(abs(scalf-1.0) > eps_rp) then
                            usenl = .true.
                            if(eel%use_fmm) usenl = fmm_list_are_near(eel,i,j)
                            if(usenl) then
                                ntmp(ipp) = ntmp(ipp) + 1
                                itmp(ntmp(ipp),ipp) = jp
                                rtmp(ntmp(ipp),ipp) = scalf
                            else
                                ntmp_far(ipp) = ntmp_far(ipp) + 1
                                itmp_far(ntmp_far(ipp),ipp) = jp
                                rtmp_far(ntmp_far(ipp),ipp) = scalf
                            end if
                        end if
                    end if
                end do
            end do
        end do

        ! Compress the list
        call compress_list(npol, itmp, ntmp(1:npol), eel%list_P_P)
        call compress_data(eel%list_P_P, rtmp, eel%scalef_P_P)
        call mallocate('make_screening_list [todo_P_P]', &
                       size(eel%scalef_P_P), eel%todo_P_P)
        !$omp parallel do
        do i=1, size(eel%scalef_P_P)
            eel%todo_P_P(i) = (abs(eel%scalef_P_P(i)) > eps_rp)
        end do
        
        if(eel%use_fmm) then
            if(sum(ntmp_far) > 0) then
                call compress_list(npol, itmp_far, ntmp_far(1:npol), eel%list_P_P_fmm_far)
                call compress_data(eel%list_P_P_fmm_far, rtmp, eel%scalef_P_P_fmm_far)
                call mallocate('make_screening_list [todo_P_P]', &
                            size(eel%scalef_P_P_fmm_far), eel%todo_P_P_fmm_far)
                !$omp parallel do
                do i=1, size(eel%scalef_P_P_fmm_far)
                    eel%todo_P_P_fmm_far(i) = (abs(eel%scalef_P_P_fmm_far(i)) > eps_rp)
                end do
            else
                ! The list is empty
                call matfree(eel%list_P_P_fmm_far)
                deallocate(eel%list_P_P_fmm_far)
            end if
        end if
        
        ! Build S P lists
        !$omp parallel do schedule(dynamic) default(shared) &
        !$omp private(i,ineigh,ij,j,jp,scalf,usenl)
        do i=1, n
            ntmp(i) = 0
            if(eel%use_fmm) ntmp_far(i) = 0
            do ineigh=1, 4
                do ij=eel%top%conn(ineigh)%ri(i), eel%top%conn(ineigh)%ri(i+1)-1
                    j = eel%top%conn(ineigh)%ci(ij)
                    jp = eel%mm_polar(j)
                    if(jp > 0) then
                        ! Needed for either amoeba or amber
                        scalf = screening_rules(eel, i, 'S', jp, 'P', 'P')
                        if(abs(scalf-1.0) > eps_rp) then
                          usenl = .true.
                            if(eel%use_fmm) usenl = fmm_list_are_near(eel,i,j)
                            if(usenl) then
                                ntmp(i) = ntmp(i) + 1
                                itmp(ntmp(i),i) = jp
                                rtmp(ntmp(i),i) = scalf
                            else
                                ntmp_far(i) = ntmp_far(i) + 1
                                itmp_far(ntmp_far(i),i) = jp
                                rtmp_far(ntmp_far(i),i) = scalf
                            end if
                        end if
                    end if
                end do
            end do
        end do

        ! Compress the list
        call compress_list(n, itmp, ntmp, eel%list_S_P_P)
        call compress_data(eel%list_S_P_P, rtmp, eel%scalef_S_P_P)
        call mallocate('make_screening_list [todo_S_P_P]', &
                       size(eel%scalef_S_P_P), eel%todo_S_P_P)
        !$omp parallel do
        do i=1, size(eel%scalef_S_P_P)
            eel%todo_S_P_P(i) = (abs(eel%scalef_S_P_P(i)) > eps_rp)
        end do
        
        if(eel%use_fmm) then
            if(sum(ntmp_far) > 0) then
                call compress_list(n, itmp_far, ntmp_far, eel%list_S_P_P_fmm_far)
                call compress_data(eel%list_S_P_P_fmm_far, rtmp_far, eel%scalef_S_P_P_fmm_far)
                call mallocate('make_screening_list [todo_S_P_P_fmm_far]', &
                            size(eel%scalef_S_P_P_fmm_far), eel%todo_S_P_P_fmm_far)
                !$omp parallel do
                do i=1, size(eel%scalef_S_P_P_fmm_far)
                    eel%todo_S_P_P_fmm_far(i) = (abs(eel%scalef_S_P_P_fmm_far(i)) > eps_rp)
                end do
            else
                ! The list is empty
                call matfree(eel%list_S_P_P_fmm_far)
                deallocate(eel%list_S_P_P_fmm_far)
            end if
        end if
        
        if(eel%amoeba) then
            !$omp parallel do schedule(dynamic) default(shared) &
            !$omp private(i,ineigh,ij,j,jp,scalf,pg_i,igrp,grp,usenl)
            do i=1, n
                ntmp(i) = 0
                if(eel%use_fmm) ntmp_far(i) = 0
                pg_i = eel%mmat_polgrp(i)
                
                do ineigh=1, 4
                    do igrp=eel%pg_conn(ineigh)%ri(pg_i), &
                            eel%pg_conn(ineigh)%ri(pg_i+1)-1
                        
                        grp = eel%pg_conn(ineigh)%ci(igrp)

                        do j=eel%polgrp_mmat%ri(grp), &
                             eel%polgrp_mmat%ri(grp+1)-1
                            jp = eel%mm_polar(j)
                            if(jp > 0) then
                                scalf = screening_rules(eel, i, 'S', jp, 'P', 'D')
                                if(abs(scalf-1.0) > eps_rp) then
                                    usenl = .true.
                                    if(eel%use_fmm) usenl = fmm_list_are_near(eel,i,j)
                                    
                                    if(usenl) then
                                        ntmp(i) = ntmp(i) + 1
                                        itmp(ntmp(i),i) = jp
                                        rtmp(ntmp(i),i) = scalf
                                    else
                                        ntmp_far(i) = ntmp_far(i) + 1
                                        itmp_far(ntmp_far(i),i) = jp
                                        rtmp_far(ntmp_far(i),i) = scalf
                                    end if
                                end if
                            end if
                        end do
                    end do
                end do
            end do

            ! Compress the list
            call compress_list(n, itmp, ntmp, eel%list_S_P_D)
            call compress_data(eel%list_S_P_D, rtmp, eel%scalef_S_P_D)
            call mallocate('make_screening_list [todo_S_P_D]', &
                           size(eel%scalef_S_P_D), eel%todo_S_P_D)
            !$omp parallel do
            do i=1, size(eel%scalef_S_P_D)
                eel%todo_S_P_D(i) = (abs(eel%scalef_S_P_D(i)) > eps_rp)
            end do
            
            if(eel%use_fmm) then
                if(sum(ntmp_far) > 0) then
                    ! Compress the list
                    call compress_list(n, itmp_far, ntmp_far, eel%list_S_P_D_fmm_far)
                    call compress_data(eel%list_S_P_D_fmm_far, rtmp_far, eel%scalef_S_P_D_fmm_far)
                    call mallocate('make_screening_list [todo_S_P_D_fmm_far]', &
                                size(eel%scalef_S_P_D_fmm_far), eel%todo_S_P_D_fmm_far)
                    !$omp parallel do
                    do i=1, size(eel%scalef_S_P_D_fmm_far)
                        eel%todo_S_P_D_fmm_far(i) = (abs(eel%scalef_S_P_D_fmm_far(i)) > eps_rp)
                    end do
                else
                    ! The list is empty
                    call matfree(eel%list_S_P_D_fmm_far)
                    deallocate(eel%list_S_P_D_fmm_far)
                end if
            end if
        end if
        
        eel%screening_list_done = .true.

        call mfree('make_screening_list [rtmp]', rtmp)
        call mfree('make_screening_list [itmp]', itmp)
        call mfree('make_screening_list [ntmp]', ntmp)
    end subroutine

    subroutine thole_init(eel)
        ! This routine compute the thole factors and stores
        ! them in a vector. TODO add reference
        ! TODO in AMOEBA should be read from FF
       
        use mod_constants, only: OMMP_VERBOSE_LOW, eps_rp
        implicit none

        type(ommp_electrostatics_type), intent(inout) :: eel
        
        integer(ip) :: i, j
        
        eel%thole = 0.0_rp
        
        do i = 1, eel%pol_atoms
            j = eel%polar_mm(i)
            eel%thole(j) = eel%pol(i) ** (1.0_rp/6.0_rp)
        end do
        
        if(.not. eel%amoeba) then
            if(eel%thole_scale > eps_rp) then
                eel%thole = eel%thole * sqrt(eel%thole_scale)
            else
                call fatal_error("Scale factor for Thole damping should be &
                                 &greater than 0.0 when non-AMOEBA FF are used")
            end if
        else
            if(eel%thole_scale > eps_rp) then
                call ommp_message("Scale factor set for thole damping is &
                    &ignored because AMOEBA FF is used", OMMP_VERBOSE_LOW)
            end if
        end if
    end subroutine thole_init

    subroutine energy_MM_MM(eel, ene)
        !! This function computes the interaction energy of 
        !! static electric multipoles
        use mod_memory, only: mallocate, mfree

        implicit none

        type(ommp_electrostatics_type), intent(inout) :: eel
        !! Electrostatics data structure
        real(rp), intent(inout) :: ene
        !! Energy (results will be added)
        real(rp) :: eMM

        integer(ip) :: i

        call prepare_fixedelec(eel)
        eMM = 0.0

        if(eel%amoeba) then
            eMM = eMM + dot_product(eel%q(1,:), eel%V_M2M)
            do i=1, 3
                eMM = eMM - dot_product(eel%q(i+1,:), eel%E_M2M(i,:))
            end do

            do i=1,6
                if(i == _xx_ .or. i == _yy_ .or. i == _zz_) then
                    ! diagonal elements
                    eMM = eMM + dot_product(eel%q(i+4,:), eel%Egrd_M2M(i,:))
                else
                    ! off-diagonal elements (are stored once, but 
                    ! should be counted twice
                    eMM = eMM + 2.0 * dot_product(eel%q(i+4,:), eel%Egrd_M2M(i,:))
                end if
            end do
        else
            eMM = eMM + dot_product(eel%q(1,:), eel%V_M2M)
        end if

        ! Since potential is computed using all the sites each 
        ! interaction is counted twice
        ene = ene + 0.5_rp * eMM
    
    end subroutine energy_MM_MM
    
    subroutine energy_MM_pol(eel, ene)
        !! This function computes the interaction energy of 
        !! static electric multipoles
        use mod_memory, only: mallocate, mfree

        implicit none

        type(ommp_electrostatics_type), intent(inout) :: eel
        !! Electrostatics data structure
        real(rp), intent(inout) :: ene
        !! Energy (results will be added)
        real(rp) :: eMM

        integer(ip) :: i

        call prepare_polelec(eel)
        if(.not. eel%ipd_done) call fatal_error("IPD should be computed before&
                                                & computing MM-Pol energy.")
        eMM = 0.0
        
        if(eel%amoeba) then
            do i=1, 3
                eMM = eMM - dot_product(eel%ipd(i,:,_amoeba_D_), &
                                        eel%E_M2D(i,:,_amoeba_P_))
            end do
        else
            do i=1, 3
                eMM = eMM - dot_product(eel%ipd(i,:,1), eel%E_M2D(i,:,1))
            end do
        end if

        ene = ene + 0.5_rp * eMM

    end subroutine energy_MM_pol

    subroutine coulomb_kernel(dr, maxder, res)
        !! This function compute the coulomb kernel for the distance vector dr and 
        !! its derivative up to the value required by maxder.
        use mod_memory, only: ip, rp
        use mod_io, only: fatal_error
        use mod_constants, only: eps_rp

        implicit none

        real(rp), intent(in) :: dr(3)
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
        if(norm2_r < eps_rp) then
            call fatal_error("Requesting Coulomb kernel for two atoms &
                             &placed in the same point, this could be &
                             &an internal bug or a problem in your input &
                             &file, please check.")
        end if
        res(1) = 1.0_rp / norm2_r
        inv_norm_sq = res(1) * res(1)

        do ii = 1, maxder
            res(ii+1) = res(ii) * inv_norm_sq
        end do
        
    end subroutine coulomb_kernel

    subroutine damped_coulomb_kernel(eel, i, j, maxder, res, dr)
        !! This subroutine computes the damped coulomb kernel between two atoms.
        !! Note that this only makes sense between two MM atoms, as it is only used
        !! to compute the field that induces the point dipoles!

        use mod_memory, only: ip, rp
        use mod_constants, only: eps_rp
        
        implicit none

        type(ommp_electrostatics_type), intent(in) :: eel
        !! Electrostatics data structure
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
        
        s = eel%thole(i) * eel%thole(j)
        
        ! Compute undamped kernels
        dr = eel%top%cmm(:,j) - eel%top%cmm(:,i)
        call coulomb_kernel(dr, maxder, res)

        if(abs(s) < eps_rp) then
            ! either thole(i) or thole(j) are zero, so the damped kernel
            ! is just the same as the regular one. Job is done.
            return
        end if

        u = 1.0_rp / (res(1) * s)
        u3 = u**3

        if(eel%amoeba .and. u3 * 0.39_rp < 50.0_rp) then
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
        else if(.not. eel%amoeba .and. res(1) > 1_rp/s) then
            ! TODO Again it is not clear to me why condition res(1) > 1_rp/s is here.
            u4 = u3*u
            if(maxder >= 1) res(2) = res(2) * (4.0_rp * u3 - 3.0_rp * u4)
            if(maxder >= 2) res(3) = res(3) * u4
            if(maxder >= 3) res(4) = res(4) * 0.2_rp * u4
            if(maxder >= 4) then
                call fatal_error("Damped Coulomb kernel (WANG) only supports up to the 4th derivative")
            end if
        end if
    end subroutine damped_coulomb_kernel

    subroutine q_elec_prop(q, dr, kernel, &
                           do_V, V, do_E, E, do_grdE, grdE, do_HE, HE)
        !! TODO
        !! Computes the electric potential of a charge \(q\) at position
        !! \(\mathbf{dr}\) from the charge itself. Pre-computed kernel should
        !! be provided as input. 
        !! The result is added to \(V\).
        !! $$ V_q = q \frac{f(r)}{||\mathbf{dr}||} $$

        implicit none

        real(rp), intent(in) :: q
        !! Charge
        real(rp), intent(in) :: dr(3)
        !! Distance vector
        real(rp), intent(in) :: kernel(:)
        !! Array of coulomb kernel (either damped or undamped)
        logical, intent(in) :: do_V, do_E, do_grdE, do_HE
        !! Flags to enable/disable calculation of different electrostatic 
        !! properties
        real(rp), intent(inout) :: V, E(3), grdE(6), HE(10)
        !! Electric potential
        
        if(do_V) then
            V =  V + kernel(1) * q
        end if
        
        if(do_E) then
            E = E + q * kernel(2) * dr
        end if

        if(do_grdE) then
            ! xx
            grdE(1) = grdE(1) + q * 3.0_rp * kernel(3) * dr(1) * dr(1) - q * kernel(2)
            ! xy
            grdE(2) = grdE(2) + q * 3.0_rp * kernel(3) * dr(1) * dr(2)
            ! yy
            grdE(3) = grdE(3) + q * 3.0_rp * kernel(3) * dr(2) * dr(2) - q * kernel(2)
            ! xz
            grdE(4) = grdE(4) + q * 3.0_rp * kernel(3) * dr(1) * dr(3)
            ! yz
            grdE(5) = grdE(5) + q * 3.0_rp * kernel(3) * dr(2) * dr(3)
            ! zz
            grdE(6) = grdE(6) + q * 3.0_rp * kernel(3) * dr(3) * dr(3) - q * kernel(2)
        end if

        if(do_HE) then
            ! xxx
            HE(1) = HE(1) + (15.0_rp * kernel(4) * dr(1) * dr(1) * dr(1) - &
                             9.0_rp * kernel(3) * dr(1)) * q
            ! xxy
            HE(2) = HE(2) + (15.0_rp * kernel(4) * dr(1) * dr(1) * dr(2) - &
                             3.0_rp * kernel(3) * dr(2)) * q
            ! xxz
            HE(3) = HE(3) + (15.0_rp * kernel(4) * dr(1) * dr(1) * dr(3) - &
                             3.0_rp * kernel(3) * dr(3)) * q
            ! xyy
            HE(4) = HE(4) + (15.0_rp * kernel(4) * dr(2) * dr(2) * dr(1) - & 
                             3.0_rp * kernel(3) * dr(1)) * q 
            ! xyz
            HE(5) = HE(5) + 15.0_rp * kernel(4) * dr(1) * dr(2) * dr(3) * q
            
            ! xzz
            HE(6) = HE(6) + (15.0_rp * kernel(4) * dr(3) * dr(3) * dr(1) - &
                             3.0_rp * kernel(3) * dr(1)) * q
            ! yyy
            HE(7) = HE(7) + (15.0_rp * kernel(4) * dr(2) * dr(2) * dr(2) - & 
                             9.0_rp * kernel(3) * dr(2)) * q 
            ! yyz
            HE(8) = HE(8) + (15.0_rp * kernel(4) * dr(2) * dr(2) * dr(3) - &
                             3.0_rp * kernel(3) * dr(3)) * q
            ! yzz
            HE(9) = HE(9) + (15.0_rp * kernel(4) * dr(3) * dr(3) * dr(2) - &
                             3.0_rp * kernel(3) * dr(2)) * q
            ! zzz
            HE(10) = HE(10) + (15.0_rp * kernel(4) * dr(3) * dr(3) * dr(3) - &
                               9.0_rp * kernel(3) * dr(3)) * q
        end if

    end subroutine q_elec_prop
    
    subroutine mu_elec_prop(mu, dr, kernel, &
                            do_V, V, do_E, E, do_grdE, grdE, do_HE, HE)
        
        implicit none

        real(rp), intent(in) :: mu(3)
        !! point dipole
        real(rp), intent(in) :: dr(3)
        !! Distance vector
        real(rp), intent(in) :: kernel(:)
        !! Array of coulomb kernel (either damped or undamped)
        logical, intent(in) :: do_V, do_E, do_grdE, do_HE
        !! Flags to enable/disable calculation of different electrostatic 
        !! properties
        real(rp), intent(inout) :: V, E(3), grdE(6), HE(10)
        !! Electric potential
        
        real(rp) :: mu_dot_dr
        
        mu_dot_dr = dot_product(mu, dr)

        if(do_V) then
            V = V + mu_dot_dr * kernel(2)
        end if
        
        if(do_E) then
            E = E + 3.0_rp * mu_dot_dr * dr * kernel(3) - mu * kernel(2)
        end if

        if(do_grdE) then
            ! xx
            grdE(1) = grdE(1) + mu_dot_dr * kernel(4) * 15.0 * dr(1) * dr(1) - &
                              (mu_dot_dr + 2.0 * mu(1) * dr(1)) * 3.0 * kernel(3)
            ! xy
            grdE(2) = grdE(2) + mu_dot_dr * kernel(4) * 15.0 * dr(1) * dr(2) - &
                              (mu(1)*dr(2) + mu(2)*dr(1)) * 3.0 * kernel(3)
            ! yy
            grdE(3) = grdE(3) + mu_dot_dr * kernel(4) * 15.0 * dr(2) * dr(2) - &
                              (mu_dot_dr + 2.0 * mu(2) * dr(2)) * 3.0 * kernel(3)
            ! xz
            grdE(4) = grdE(4) + mu_dot_dr * kernel(4) * 15.0 * dr(1) * dr(3) - &
                              (mu(1)*dr(3) + mu(3)*dr(1)) * 3.0 * kernel(3)
            ! yz
            grdE(5) = grdE(5) + mu_dot_dr * kernel(4) * 15.0 * dr(2) * dr(3) - &
                              (mu(2)*dr(3) + mu(3)*dr(2)) * 3.0 * kernel(3)
            ! zz
            grdE(6) = grdE(6) + mu_dot_dr * kernel(4) * 15.0 * dr(3) * dr(3) - &
                              (mu_dot_dr + 2.0 * mu(3) * dr(3)) * 3.0 * kernel(3)
        end if

        if(do_HE) then
            ! xxx
            HE(1) = HE(1) + 105.0_rp * mu_dot_dr * kernel(5) * dr(1)*dr(1)*dr(1) &
                          - 45.0_rp * kernel(4) * dr(1) * (mu(1)*dr(1) + mu_dot_dr) &
                          + 9.0_rp * kernel(3) * mu(1)
            ! xxy
            HE(2) = HE(2) + 105.0_rp * kernel(5) * mu_dot_dr * dr(1)*dr(1)*dr(2) &
                          - 15.0_rp * kernel(4) * (mu(2)*dr(1)*dr(1) + &
                                                   2.0_rp *mu(1)*dr(1)*dr(2) + &
                                                   mu_dot_dr*dr(2)) &
                          + 3.0_rp * kernel(3) * mu(2)

            ! xxz
            HE(3) = HE(3) + 105.0_rp * kernel(5) * mu_dot_dr * dr(1)*dr(1)*dr(3) &
                          - 15.0_rp * kernel(4) * (mu(3)*dr(1)*dr(1) + &
                                                   2.0_rp *mu(1)*dr(1)*dr(3) + &
                                                   mu_dot_dr*dr(3)) &
                          + 3.0_rp * kernel(3) * mu(3)
            ! xyy
            HE(4) = HE(4) + 105.0_rp * kernel(5) * mu_dot_dr * dr(2)*dr(2)*dr(1) &
                          - 15.0_rp * kernel(4) * (mu(1)*dr(2)*dr(2) + &
                                                   2.0_rp *mu(2)*dr(2)*dr(1) + &
                                                   mu_dot_dr*dr(1)) &
                          + 3.0_rp * kernel(3) * mu(1)
            ! xyz
            HE(5) = HE(5) + 105.0_rp * mu_dot_dr * kernel(5) * dr(1)*dr(2)*dr(3) &
                          - 15.0_rp * kernel(4) * (mu(1)*dr(2)*dr(3) + &
                                                   dr(1)*mu(2)*dr(3) + &
                                                   dr(1)*dr(2)*mu(3))
            ! xzz
            HE(6) = HE(6) + 105.0_rp * kernel(5) * mu_dot_dr * dr(3)*dr(3)*dr(1) &
                          - 15.0_rp * kernel(4) * (mu(1)*dr(3)*dr(3) + &
                                                   2.0_rp *mu(3)*dr(3)*dr(1) + &
                                                   mu_dot_dr*dr(1)) &
                          + 3.0_rp * kernel(3) * mu(1)
            ! yyy
            HE(7) = HE(7) + 105.0_rp * mu_dot_dr * kernel(5) * dr(2)*dr(2)*dr(2) &
                          - 45.0_rp * kernel(4) * dr(2) * (mu(2)*dr(2) + mu_dot_dr) &
                          + 9.0_rp * kernel(3) * mu(2)
            ! yyz
            HE(8) = HE(8) + 105.0_rp * kernel(5) * mu_dot_dr * dr(2)*dr(2)*dr(3) &
                          - 15.0_rp * kernel(4) * (mu(3)*dr(2)*dr(2) + &
                                                   2.0_rp *mu(2)*dr(2)*dr(3) + &
                                                   mu_dot_dr*dr(3)) &
                          + 3.0_rp * kernel(3) * mu(3)
            ! yzz
            HE(9) = HE(9) + 105.0_rp * kernel(5) * mu_dot_dr * dr(3)*dr(3)*dr(2) &
                          - 15.0_rp * kernel(4) * (mu(2)*dr(3)*dr(3) + &
                                                   2.0_rp *mu(3)*dr(3)*dr(2) + &
                                                   mu_dot_dr*dr(2)) &
                          + 3.0_rp * kernel(3) * mu(2)
            ! zzz
            HE(10) = HE(10) + 105.0_rp * mu_dot_dr * kernel(5) * dr(3)*dr(3)*dr(3) &
                            - 45.0_rp * kernel(4) * dr(3) * (mu(3)*dr(3) + mu_dot_dr) &
                            + 9.0_rp * kernel(3) * mu(3)
        end if
    end subroutine mu_elec_prop
    
    subroutine quad_elec_prop(quad, dr, kernel, &
                              do_V, V, do_E, E, do_grdE, grdE, do_HE, HE)
        
        implicit none

        real(rp), intent(in) :: quad(6)
        !! point quadrupole stored as (xx, xy, yy, xz, yz, zz)
        real(rp), intent(in) :: dr(3)
        !! Distance vector
        real(rp), intent(in) :: kernel(:)
        !! Array of coulomb kernel (either damped or undamped)
        logical, intent(in) :: do_V, do_E, do_grdE, do_HE
        !! Flags to enable/disable calculation of different electrostatic 
        !! properties
        real(rp), intent(inout) :: V, E(3), grdE(6), HE(10)
        !! Electric potential
        
        real(rp) :: quadxr(3), quadxr_dot_r

        quadxr(1) = quad(1)*dr(1) + quad(2)*dr(2) + quad(4)*dr(3)
        quadxr(2) = quad(2)*dr(1) + quad(3)*dr(2) + quad(5)*dr(3)
        quadxr(3) = quad(4)*dr(1) + quad(5)*dr(2) + quad(6)*dr(3)

        quadxr_dot_r = dot_product(quadxr, dr)

        if(do_V) then
            V = V + 3.0_rp * quadxr_dot_r * kernel(3)
        end if
        
        if(do_E) then
            E = E + 15.0_rp * quadxr_dot_r * dr * kernel(4) &
                - 6.0_rp * quadxr * kernel(3)
        end if

        if(do_grdE) then
            ! xx
            grdE(1) = grdE(1) + 105.0 * kernel(5) * quadxr_dot_r * dr(1)*dr(1) + &
                              6.0 * kernel(3) * quad(1) - &
                              15.0*kernel(4)*(quadxr_dot_r+4.0*quadxr(1)*dr(1))
            ! xy
            grdE(2) = grdE(2) + 105.0 * kernel(5) * quadxr_dot_r * dr(1)*dr(2) + &
                              6.0 * kernel(3) * quad(2) - &
                              30.0*kernel(4)*(quadxr(1)*dr(2)+quadxr(2)*dr(1))
            ! yy
            grdE(3) = grdE(3) + 105.0 * kernel(5) * quadxr_dot_r * dr(2)*dr(2) + &
                              6.0 * kernel(3) * quad(3) - &
                              15.0*kernel(4)*(quadxr_dot_r+4.0*quadxr(2)*dr(2))
            ! xz
            grdE(4) = grdE(4) + 105.0 * kernel(5) * quadxr_dot_r * dr(1)*dr(3) + &
                              6.0 * kernel(3) * quad(4) - &
                              30.0*kernel(4)*(quadxr(1)*dr(3)+quadxr(3)*dr(1))
            ! yz
            grdE(5) = grdE(5) + 105.0 * kernel(5) * quadxr_dot_r * dr(2)*dr(3) + &
                              6.0 * kernel(3) * quad(5) - &
                              30.0*kernel(4)*(quadxr(2)*dr(3)+quadxr(3)*dr(2))
            ! zz
            grdE(6) = grdE(6) + 105.0 * kernel(5) * quadxr_dot_r * dr(3)*dr(3) + &
                              6.0 * kernel(3) * quad(6) - &
                              15.0*kernel(4)*(quadxr_dot_r+4.0*quadxr(3)*dr(3))
        end if

        if(do_HE) then
            ! 3
            HE(_xxx_) = HE(_xxx_) + 945*kernel(6)*dr(_x_)**3 * quadxr_dot_r &
                          - 315*kernel(5)*dr(_x_)*(2*quadxr(_x_)*dr(_x_) + quadxr_dot_r) &
                          + 90*kernel(4)*(quad(_xx_)*dr(_x_) + quadxr(_x_))   
            
            HE(_yyy_) = HE(_yyy_) + 945*kernel(6)*dr(_y_)**3 * quadxr_dot_r &
                          - 315*kernel(5)*dr(_y_)*(2*quadxr(_y_)*dr(_y_) + quadxr_dot_r) &
                          + 90*kernel(4)*(quad(_yy_)*dr(_y_) + quadxr(_y_))   
            
            HE(_zzz_) = HE(_zzz_) + 945*kernel(6)*dr(_z_)**3 * quadxr_dot_r &
                          - 315*kernel(5)*dr(_z_)*(2*quadxr(_z_)*dr(_z_) + quadxr_dot_r) &
                          + 90*kernel(4)*(quad(_zz_)*dr(_z_) + quadxr(_z_))   
           
            ! 2 + 1
            HE(_xxy_) = HE(_xxy_) + 945*kernel(6)*dr(_x_)**2*dr(_y_)*quadxr_dot_r &
                          - 105*kernel(5)*(4*quadxr(_x_)*dr(_x_)*dr(_y_) + &
                                           2*quadxr(_y_)*dr(_x_)*dr(_x_) + &
                                           dr(_y_)*quadxr_dot_r) &
                          + 30*kernel(4)*(quad(_xx_)*dr(_y_) + 2*quad(_xy_)*dr(_x_) + quadxr(_y_)) 
            
            HE(_xxz_) = HE(_xxz_) + 945*kernel(6)*dr(_x_)**2*dr(_z_)*quadxr_dot_r &
                          - 105*kernel(5)*(4*quadxr(_x_)*dr(_x_)*dr(_z_) + &
                                           2*quadxr(_z_)*dr(_x_)*dr(_x_) + &
                                           dr(_z_)*quadxr_dot_r) &
                          + 30*kernel(4)*(quad(_xx_)*dr(_z_) + 2*quad(_xz_)*dr(_x_) + quadxr(_z_)) 

            HE(_yyx_) = HE(_yyx_) + 945*kernel(6)*dr(_y_)**2*dr(_x_)*quadxr_dot_r &
                          - 105*kernel(5)*(4*quadxr(_y_)*dr(_y_)*dr(_x_) + &
                                           2*quadxr(_x_)*dr(_y_)*dr(_y_) + &
                                           dr(_x_)*quadxr_dot_r) &
                          + 30*kernel(4)*(quad(_yy_)*dr(_x_) + 2*quad(_xy_)*dr(_y_) + quadxr(_x_)) 
            
            HE(_yyz_) = HE(_yyz_) + 945*kernel(6)*dr(_y_)**2*dr(_z_)*quadxr_dot_r &
                          - 105*kernel(5)*(4*quadxr(_y_)*dr(_y_)*dr(_z_) + &
                                           2*quadxr(_z_)*dr(_y_)*dr(_y_) + &
                                           dr(_z_)*quadxr_dot_r) &
                          + 30*kernel(4)*(quad(_yy_)*dr(_z_) + 2*quad(_zy_)*dr(_y_) + quadxr(_z_)) 

            HE(_zzx_) = HE(_zzx_) + 945*kernel(6)*dr(_z_)**2*dr(_x_)*quadxr_dot_r &
                          - 105*kernel(5)*(4*quadxr(_z_)*dr(_z_)*dr(_x_) + &
                                           2*quadxr(_x_)*dr(_z_)*dr(_z_) + &
                                           dr(_x_)*quadxr_dot_r) &
                          + 30*kernel(4)*(quad(_zz_)*dr(_x_) + 2*quad(_zx_)*dr(_z_) + quadxr(_x_)) 
            
            HE(_zzy_) = HE(_zzy_) + 945*kernel(6)*dr(_z_)**2*dr(_y_)*quadxr_dot_r &
                          - 105*kernel(5)*(4*quadxr(_z_)*dr(_z_)*dr(_y_) + &
                                           2*quadxr(_y_)*dr(_z_)*dr(_z_) + &
                                           dr(_y_)*quadxr_dot_r) &
                          + 30*kernel(4)*(quad(_zz_)*dr(_y_) + 2*quad(_zy_)*dr(_z_) + quadxr(_y_)) 
            
            ! 1 + 1 + 1
            HE(_xyz_) = HE(_xyz_) + 945*kernel(6)*dr(_x_)*dr(_y_)*dr(_z_)*quadxr_dot_r &
                          - 210*kernel(5)*(quadxr(_x_)*dr(_y_)*dr(_z_) + &
                                           quadxr(_y_)*dr(_x_)*dr(_z_) + &
                                           quadxr(_z_)*dr(_x_)*dr(_y_)) &
                          + 30*kernel(4)*(quad(_xy_)*dr(_z_) + quad(_xz_)*dr(_y_) + quad(_yz_)*dr(_x_))  
        end if
    end subroutine quad_elec_prop

    subroutine prepare_fixedelec(eel, arg_dogg) 
        !! This function allocate and populate array of electrostatic 
        !! properties of static multipoles at static multipoles sites.
        !! It should be called blindly before any calculation that requires
        !! V_M2M etc.
   
        use mod_memory, only: mallocate

        implicit none

        type(ommp_electrostatics_type), intent(inout) :: eel
        logical, optional, intent(in) :: arg_dogg

        integer(ip) :: mm_atoms
        logical :: do_gg

        mm_atoms = eel%top%mm_atoms

        !! TODO Improve logic
        do_gg = .false.
        if(present(arg_dogg)) then
            if(arg_dogg) do_gg = .true.
        end if

        if(.not. do_gg .and. eel%M2M_done) return
        if(do_gg .and. eel%M2M_done .and. eel%M2Mgg_done) return

        if(eel%amoeba) then
            if(.not. allocated(eel%V_M2M)) then
                call mallocate('prepare_fixedelec [V_M2M]', mm_atoms, eel%V_M2M)
            end if

            if(.not. allocated(eel%E_M2M)) then
                call mallocate('prepare_fixedelec [E_M2M]', 3, mm_atoms, eel%E_M2M)
            end if

            if(.not. allocated(eel%Egrd_M2M)) then
                call mallocate('prepare_fixedelec [Egrd_M2M]', 6, mm_atoms, eel%Egrd_M2M)
            end if
            
            if(do_gg .and. .not. allocated(eel%EHes_M2M)) then
                call mallocate('prepare_fixedelec [EHes_M2M]', 10, mm_atoms, eel%EHes_M2M)
            end if
            
            eel%V_M2M = 0.0_rp
            eel%E_M2M = 0.0_rp
            eel%Egrd_M2M = 0.0_rp
            if(do_gg) eel%EHes_M2M = 0.0_rp
            
            if(do_gg) then
                call elec_prop_M2M(eel, .true., .true., .true., .true.)
            else
                call elec_prop_M2M(eel, .true., .true., .true., .false.)
            end if
        else
            if(.not. allocated(eel%V_M2M)) then
                call mallocate('prepare_fixedelec [V_M2M]', mm_atoms, eel%V_M2M)
            end if

            if(do_gg .and. .not. allocated(eel%E_M2M)) then
                call mallocate('prepare_fixedelec [E_M2M]', 3, mm_atoms, eel%E_M2M)
            end if

            eel%V_M2M = 0.0_rp
            if(do_gg) eel%E_M2M = 0.0_rp

            if(do_gg) then
                call elec_prop_M2M(eel, .true., .true., .false., .false.)
            else
                call elec_prop_M2M(eel, .true., .false., .false., .false.)
            end if
        end if
        
        eel%M2M_done = .true.
        if(do_gg) eel%M2Mgg_done = .true.

    end subroutine prepare_fixedelec

    subroutine prepare_polelec(eel, arg_dogg)
        use mod_memory, only: mallocate
        implicit none

        type(ommp_electrostatics_type), intent(inout) :: eel
        logical, optional, intent(in) :: arg_dogg

        logical :: do_gg

        do_gg = .false.
        if(present(arg_dogg)) then
            if(arg_dogg) do_gg = .true.
        end if
        
        if(.not. do_gg .and. eel%M2D_done) return
        if(do_gg .and. eel%M2D_done .and. eel%M2Dgg_done) return
        
        if(.not. allocated(eel%V_M2D) .and. eel%use_fmm) then
            ! This is needed just as a placeholder in fmm calls
            call mallocate('prepare_polelec [V_M2D]', eel%pol_atoms, &
                            eel%n_ipd, eel%V_M2D)
        end if
        
        if(.not. allocated(eel%E_M2D)) then
            call mallocate('prepare_polelec [E_M2D]', 3, eel%pol_atoms, &
                            eel%n_ipd, eel%E_M2D)
        end if
        
        if(do_gg .or. eel%use_fmm) then
            if(.not. eel%ipd_done .and. do_gg) call fatal_error("IPD should be computed &
                &before computing analytical geometrical gradients of &
                &polarization energy.")
            
            if(.not. allocated(eel%Egrd_M2D)) then
                call mallocate('prepare_polelec [Egrd_M2D]', 6, eel%pol_atoms, &
                                eel%n_ipd, eel%Egrd_M2D)
            end if
            
            if(.not. allocated(eel%V_D2D) .and. eel%use_fmm) then
                call mallocate('prepare_polelec [V_D2D]', eel%pol_atoms, &
                               eel%n_ipd, eel%V_D2D)
            end if
            
            if(.not. allocated(eel%E_D2D) .and. eel%use_fmm) then
                call mallocate('prepare_polelec [E_D2D]', 3, eel%pol_atoms, &
                               eel%n_ipd, eel%E_D2D)
            end if
            
            if(.not. allocated(eel%Egrd_D2D)) then
                call mallocate('prepare_polelec [Egrd_D2D]', 6, eel%pol_atoms, &
                               eel%n_ipd, eel%Egrd_D2D)
            end if
            
            if(.not. allocated(eel%EHes_D2D) .and. eel%use_fmm) then
                call mallocate('prepare_polelec [EHes_D2D]', 10, eel%pol_atoms, &
                               eel%n_ipd, eel%EHes_D2D)
            end if
            
            if(.not. allocated(eel%V_D2M) .and. eel%use_fmm) then
                call mallocate('prepare_polelec [V_D2M]', eel%top%mm_atoms, &
                               eel%V_D2M)
            end if
            
            if(.not. allocated(eel%E_D2M)) then
                call mallocate('prepare_polelec [E_D2M]', 3, eel%top%mm_atoms, &
                               eel%E_D2M)
            end if
            if(.not. allocated(eel%Egrd_D2M) .and. eel%amoeba) then
                call mallocate('prepare_polelec [Egrd_D2M]', 6, eel%top%mm_atoms, &
                               eel%Egrd_D2M)
            end if
            if(.not. allocated(eel%EHes_D2M) .and. eel%amoeba) then
                call mallocate('prepare_polelec [E_D2M]', 10, eel%top%mm_atoms, &
                               eel%EHes_D2M)
            end if

        end if
        
        if(.not. do_gg) then
            eel%E_M2D = 0.0_rp
            call elec_prop_M2D(eel, .false., .true., .false., .false.)
        else
            eel%E_M2D = 0.0_rp
            eel%Egrd_M2D = 0.0_rp
            call elec_prop_M2D(eel, .false., .true., .true., .false.)

            eel%E_D2M = 0.0_rp
            eel%Egrd_D2D = 0.0_rp
            
            if(eel%amoeba) then
                eel%Egrd_D2M = 0.0_rp
                eel%EHes_D2M = 0.0_rp
                call elec_prop_D2M(eel, 'P', .false., .true., .true., .true.)
                call elec_prop_D2M(eel, 'D', .false., .true., .true., .true.)
        
                eel%E_D2M = eel%E_D2M * 0.5
                eel%Egrd_D2M = eel%Egrd_D2M * 0.5
                eel%EHes_D2M = eel%EHes_D2M * 0.5

                call elec_prop_D2D(eel, 'P', .false., .false., .true., .false.)
                call elec_prop_D2D(eel, 'D', .false., .false., .true., .false.)
            else
                call elec_prop_D2M(eel, '-', .false., .true., .false., .false.)
                call elec_prop_D2D(eel, '-', .false., .false., .true., .false.)
            end if
        end if
        
        if(do_gg) eel%M2Dgg_done = .true.
        eel%M2D_done = .true.
    end subroutine prepare_polelec

    subroutine preapare_fmm_static(eel)
        use mod_memory, only: mallocate, mfree
        implicit none
        
        type(ommp_electrostatics_type), intent(inout) :: eel
        !! Electrostatics data structure
        real(rp), allocatable :: tmp_q(:), tmp_mu(:,:), tmp_quad(:,:)
        
        if(.not. eel%fmm_static_done) then
            call mallocate('prepare_fmm_static [tmp_q]', eel%top%mm_atoms, tmp_q)
            tmp_q(:) = eel%q(1,:)
            if(eel%amoeba) then
                call mallocate('prepare_fmm_static [tmp_mu]', 3, eel%top%mm_atoms, tmp_mu)
                tmp_mu(:,:) = eel%q(2:4,:)
                call mallocate('prepare_fmm_static [tmp_quad]', 6, eel%top%mm_atoms, tmp_quad)
                tmp_quad(:,:) = eel%q(5:10,:)
            else
                call mallocate('prepare_fmm_static [tmp_mu]', 3, 1, tmp_mu)
                call mallocate('prepare_fmm_static [tmp_quad]', 6, 1, tmp_quad)
            end if

            call fmm_solve_for_multipoles(eel%fmm_static, &
                                        tmp_q, logical(.true., lp), &
                                        tmp_mu, eel%amoeba, &
                                        tmp_quad, eel%amoeba)
            eel%fmm_static_done = .true.
            
            call mfree('prepare_fmm_static [tmp_q]', tmp_q)
            call mfree('prepare_fmm_static [tmp_mu]', tmp_mu)
            call mfree('prepare_fmm_static [tmp_quad]', tmp_quad)
        end if
    end subroutine
    
    subroutine prepare_fmm_ext_ipd(eel, fmm, ipd)
        use mod_memory, only: mallocate, mfree
        implicit none
        
        type(ommp_electrostatics_type), intent(in) :: eel
        !! Electrostatics data structure
        type(fmm_type), intent(inout) :: fmm
        !! fmm object used to run the calculation
        real(rp), intent(in) :: ipd(3, eel%pol_atoms)
        !! External induced point dipoles at polarizable sites

        real(rp), allocatable :: tmp_mu(:,:)
        real(rp) :: fake_q(1), fake_quad(6,1)
        integer(ip) :: i
        
        call mallocate('prepare_fmm_ext_ipd [tmp_mu]', 3, eel%top%mm_atoms, tmp_mu)
        tmp_mu = 0.0
        do i=1, eel%pol_atoms
            tmp_mu(:,eel%polar_mm(i)) = ipd(:,i)
        end do
        
        call time_push
        call fmm_solve_for_multipoles(fmm, &
                                      fake_q, logical(.false., lp), &
                                      tmp_mu, logical(.true., lp), &
                                      fake_quad, logical(.false., lp))
        call time_pull('fmm_solve_for_multipoles')
        call mfree('prepare_fmm_ext_ipd [tmp_mu]', tmp_mu)
    end subroutine

    subroutine prepare_fmm_ipd(eel, knd)
        use mod_memory, only: mallocate, mfree
        implicit none
        
        type(ommp_electrostatics_type), intent(inout) :: eel
        !! Electrostatics data structure
        integer(ip), intent(in) :: knd
        !! Index for the set of dipoles to use.

        if(eel%ipd_done) then
            if(.not. eel%fmm_ipd_done(knd)) then
                call prepare_fmm_ext_ipd(eel, eel%fmm_ipd(knd), eel%ipd(:, :, knd))
                eel%fmm_ipd_done(knd) = .true.
            end if
        else
            call fatal_error('Converged IPD are needed to call prepare_fmm_ipd.')
        end if
    end subroutine

    subroutine elec_prop_M2M(eel, do_V, do_E, do_Egrd, do_EHes)
        !! Computes the electric potential, field and field gradients of 
        !! static multipoles at all sites (polarizable sites are a 
        !! subset of static ones)
        implicit none
        
        type(ommp_electrostatics_type), intent(inout) :: eel
        !! Electrostatics data structure
        logical, intent(in) :: do_V, do_E, do_Egrd, do_EHes
        !! Flags to enable/disable the calculation of different components


        real(rp) :: kernel(6), dr(3), tmpV, tmpE(3), tmpEgr(6), tmpHE(10), scalf
        integer(ip) :: i, j, idx, ikernel
        logical :: to_do, to_scale
        type(ommp_topology_type), pointer :: top

        top => eel%top

        if(do_EHes) then
            ikernel = 3 
        elseif(do_Egrd) then
            ikernel = 2
        elseif(do_E) then
            ikernel = 1
        elseif(do_V) then
            ikernel = 0
        else
            return
        end if
        if(eel%amoeba) ikernel = ikernel + 2

        if(eel%use_fmm) then
            call preapare_fmm_static(eel)

            ! TODO !
            if(allocated(eel%list_S_S_fmm_far)) &
                call fatal_error('Currently not supported')
            !$omp parallel do default(shared) schedule(dynamic) &
            !$omp private(i,j,idx,scalf,dr,kernel,tmpV,tmpE,tmpEgr,tmpHE)
            do i=1, top%mm_atoms 
                call cart_prop_at_ipart(eel%fmm_static, i, &
                                        do_V, eel%V_M2M(i), &
                                        do_E, eel%E_M2M(:,i), &
                                        do_Egrd, eel%Egrd_M2M(:,i), &
                                        do_EHes, eel%EHes_M2M(:,i))
                
                ! Now removed the scaled parts only
                do idx=eel%list_S_S%ri(i), eel%list_S_S%ri(i+1)-1
                    j = eel%list_S_S%ci(idx)
                    scalf = 1.0 - eel%scalef_S_S(idx)
                        
                    dr = top%cmm(:,i) - top%cmm(:, j)
                    call coulomb_kernel(dr, ikernel, kernel)
                    
                    if(do_V) tmpV = 0.0_rp
                    if(do_E) tmpE = 0.0_rp
                    if(do_Egrd) tmpEgr = 0.0_rp
                    if(do_EHes) tmpHE = 0.0_rp

                    call q_elec_prop(eel%q(1,j), dr, kernel, &
                                        do_V, tmpV, & 
                                        do_E, tmpE, &
                                        do_Egrd, tmpEgr, &
                                        do_EHes, tmpHE)
                    if(eel%amoeba) then
                        call mu_elec_prop(eel%q(2:4,j), dr, kernel, &
                                            do_V, tmpV, & 
                                            do_E, tmpE, &
                                            do_Egrd, tmpEgr, &
                                            do_EHes, tmpHE)

                        call quad_elec_prop(eel%q(5:10,j), dr, kernel, &
                                            do_V, tmpV, & 
                                            do_E, tmpE, &
                                            do_Egrd, tmpEgr, &
                                            do_EHes, tmpHE)
                    end if

                    if(do_V) eel%V_M2M(i) = eel%V_M2M(i) - tmpV * scalf
                    if(do_E) eel%E_M2M(:,i) = eel%E_M2M(:,i) - tmpE * scalf
                    if(do_Egrd) eel%Egrd_M2M(:,i) = eel%Egrd_M2M(:,i) - tmpEgr * scalf
                    if(do_EHes) eel%EHes_M2M(:,i) = eel%EHes_M2M(:,i) - tmpHE * scalf

                end do
            end do
        else
        if(eel%amoeba) then
            !$omp parallel do default(shared) schedule(dynamic) &
            !$omp private(i,j,idx,to_do,to_scale,scalf,dr,kernel,tmpV,tmpE,tmpEgr,tmpHE)
            do j=1, top%mm_atoms
                ! loop on sources
                do i=1, top%mm_atoms
                    if(j == i) cycle
                    !loop on target

                    to_do = .true.
                    to_scale = .false.
                    scalf = 1.0

                    ! Check if the element should be scaled
                    do idx=eel%list_S_S%ri(i), eel%list_S_S%ri(i+1)-1
                        if(eel%list_S_S%ci(idx) == j) then
                            to_scale = .true.
                            exit
                        end if
                    end do

                    !If it should set the correct variables
                    if(to_scale) then
                        to_do = eel%todo_S_S(idx)
                        scalf = eel%scalef_S_S(idx)
                    end if

                    if(to_do) then
                        dr = top%cmm(:,j) - top%cmm(:, i)
                        call coulomb_kernel(dr, ikernel, kernel)
                        
                        if(do_V) tmpV = 0.0_rp
                        if(do_E) tmpE = 0.0_rp
                        if(do_Egrd) tmpEgr = 0.0_rp
                        if(do_EHes) tmpHE = 0.0_rp

                        call q_elec_prop(eel%q(1,i), dr, kernel, &
                                         do_V, tmpV, & 
                                         do_E, tmpE, &
                                         do_Egrd, tmpEgr, &
                                         do_EHes, tmpHE)

                        call mu_elec_prop(eel%q(2:4,i), dr, kernel, &
                                          do_V, tmpV, & 
                                          do_E, tmpE, &
                                          do_Egrd, tmpEgr, &
                                          do_EHes, tmpHE)

                        call quad_elec_prop(eel%q(5:10,i), dr, kernel, &
                                            do_V, tmpV, & 
                                            do_E, tmpE, &
                                            do_Egrd, tmpEgr, &
                                            do_EHes, tmpHE)

                        if(to_scale) then
                            if(do_V) eel%V_M2M(j) = eel%V_M2M(j) + tmpV * scalf
                            if(do_E) eel%E_M2M(:,j) = eel%E_M2M(:,j) + tmpE * scalf
                            if(do_Egrd) eel%Egrd_M2M(:,j) = eel%Egrd_M2M(:,j) + tmpEgr * scalf
                            if(do_EHes) eel%EHes_M2M(:,j) = eel%EHes_M2M(:,j) + tmpHE * scalf
                        else
                            if(do_V) eel%V_M2M(j) = eel%V_M2M(j) + tmpV
                            if(do_E) eel%E_M2M(:,j) = eel%E_M2M(:,j) + tmpE
                            if(do_Egrd) eel%Egrd_M2M(:,j) = eel%Egrd_M2M(:,j) + tmpEgr
                            if(do_EHes) eel%EHes_M2M(:,j) = eel%EHes_M2M(:,j) + tmpHE
                        end if
                    end if
                end do
            end do
        else
            !$omp parallel do default(shared) schedule(dynamic) &
            !$omp private(i,j,idx,to_do,to_scale,scalf,dr,kernel,tmpV,tmpE,tmpEgr,tmpHE) 
            do j=1, top%mm_atoms
                ! loop on sources
                do i=1, top%mm_atoms
                    if(j == i) cycle
                    !loop on target
                    to_do = .true.
                    to_scale = .false.
                    scalf = 1.0

                    ! Check if the element should be scaled
                    do idx=eel%list_S_S%ri(i), eel%list_S_S%ri(i+1)-1
                        if(eel%list_S_S%ci(idx) == j) then
                            to_scale = .true.
                            exit
                        end if
                    end do

                    !If it should set the correct variables
                    if(to_scale) then
                        to_do = eel%todo_S_S(idx)
                        scalf = eel%scalef_S_S(idx)
                    end if
                    
                    if(to_do) then
                        dr = top%cmm(:,j) - top%cmm(:, i)
                        call coulomb_kernel(dr, ikernel, kernel)
                        
                        if(do_V) tmpV = 0.0_rp
                        if(do_E) tmpE = 0.0_rp
                        if(do_Egrd) tmpEgr = 0.0_rp
                        if(do_EHes) tmpHE = 0.0_rp

                        call q_elec_prop(eel%q(1,i), dr, kernel, &
                                         do_V, tmpV, & 
                                         do_E, tmpE, &
                                         do_Egrd, tmpEgr, &
                                         do_EHes, tmpHE)

                        if(to_scale) then
                            if(do_V) eel%V_M2M(j) = eel%V_M2M(j) + tmpV * scalf
                            if(do_E) eel%E_M2M(:,j) = eel%E_M2M(:,j) + tmpE * scalf
                            if(do_Egrd) eel%Egrd_M2M(:,j) = eel%Egrd_M2M(:,j) + tmpEgr * scalf
                            if(do_EHes) eel%EHes_M2M(:,j) = eel%EHes_M2M(:,j) + tmpHE * scalf
                        else
                            if(do_V) eel%V_M2M(j) = eel%V_M2M(j) + tmpV
                            if(do_E) eel%E_M2M(:,j) = eel%E_M2M(:,j) + tmpE
                            if(do_Egrd) eel%Egrd_M2M(:,j) = eel%Egrd_M2M(:,j) + tmpEgr
                            if(do_EHes) eel%EHes_M2M(:,j) = eel%EHes_M2M(:,j) + tmpHE
                        end if
                    end if
                end do
            end do
        end if
        end if
    end subroutine elec_prop_M2M
    
    subroutine field_extD2D(eel, ext_ipd, E)
        !! Computes the electric field of a trial set of induced point dipoles
        !! at polarizable sites. This is intended to be used as matrix-vector
        !! routine in the solution of the linear system.
        
        implicit none

        type(ommp_electrostatics_type), intent(in) :: eel
        !! Data structure for electrostatic part of the system
        real(rp), intent(in) :: ext_ipd(3, eel%pol_atoms)
        !! External induced point dipoles at polarizable sites
        real(rp), intent(inout) :: E(3, eel%pol_atoms)
        !! Electric field (results will be added)

        integer(ip) :: i, j, ipol, jpol, ij, idx
        logical :: to_scale, to_do
        real(rp) :: kernel(5), dr(3), tmpV, tmpE(3), tmpEgr(6), tmpHE(10), scalf
        type(fmm_type), allocatable :: fmm_ipd

        if(eel%use_fmm) then
            call time_push()
            call time_push()
            allocate(fmm_ipd)
            call time_pull("Allocation")
            call time_push()
            call fmm_init(fmm_ipd, eel%fmm_maxl, eel%tree)
            call time_pull("FMM initialization")
            call time_push()
            call prepare_fmm_ext_ipd(eel, fmm_ipd, ext_ipd)
            call time_pull("FMM solution")

            call time_push()
            !$omp parallel do default(shared) schedule(dynamic) &
            !$omp private(i,j,ij,ipol,jpol,idx,dr,kernel,to_do,to_scale,scalf,tmpV,tmpE,tmpEgr,tmpHE) 
            do ipol=1, eel%pol_atoms 
                i = eel%polar_mm(ipol)
                tmpE = 0.0
                call cart_propfar_at_ipart(fmm_ipd, i, &
                                           .false., tmpV, &
                                           .true. , tmpE, &
                                           .false., tmpEgr, &
                                           .false., tmpHE)
                E(:, ipol) = tmpE
            end do
            call time_pull("Far field")
            
            call time_push()
            !$omp parallel do default(shared) schedule(dynamic) &
            !$omp private(i,j,ij,ipol,jpol,idx,dr,kernel,to_do,to_scale,scalf,tmpV,tmpE,tmpEgr,tmpHE) 
            do ipol=1, eel%pol_atoms 
                i = eel%polar_mm(ipol)
                 
                ! Near field is computed internally because dumped kernel is required
                do ij=eel%fmm_near_field_list%ri(i), eel%fmm_near_field_list%ri(i+1)-1
                    j = eel%fmm_near_field_list%ci(ij)
                    jpol = eel%polar_mm(j)
                    ! If the atom is not polarizable, skip
                    if(jpol < 1) cycle 

                    !loop on target
                    to_do = .true.
                    to_scale = .false.
                    scalf = 1.0

                    ! Check if the element should be scaled
                    do idx=eel%list_P_P%ri(ipol), eel%list_P_P%ri(ipol+1)-1
                        if(eel%list_P_P%ci(idx) == jpol) then
                            to_scale = .true.
                            exit
                        end if
                    end do

                    !If it should set the correct variables
                    if(to_scale) then
                        to_do = eel%todo_P_P(idx)
                        scalf = eel%scalef_P_P(idx)
                    end if
                    
                    if(to_do) then
                        call damped_coulomb_kernel(eel, j, i,& 
                                                   2, kernel(1:3), dr)
                        
                        tmpE = 0.0_rp

                        call mu_elec_prop(ext_ipd(:,jpol), dr, kernel, .false., tmpV, &
                                        .true., tmpE, .false., tmpEgr, & 
                                        .false., tmpHE)
                        if(to_scale) then
                            E(:, ipol) = E(:, ipol) + tmpE * scalf
                        else
                            E(:, ipol) = E(:, ipol) + tmpE
                        end if
                    end if
                end do
            end do
            call time_pull("Near field")
            call time_push
            
            if(allocated(eel%list_P_P_fmm_far)) then
                ! Now remove screened interactions from far-field
                do ipol=1, eel%pol_atoms
                    i = eel%polar_mm(ipol)
                    
                    do idx=eel%list_P_P_fmm_far%ri(ipol), eel%list_P_P_fmm_far%ri(ipol+1)-1
                        jpol = eel%list_P_P_fmm_far%ci(idx)
                        j = eel%polar_mm(jpol)

                        scalf = 1.0 - eel%scalef_P_P_fmm_far(idx)
                        
                        call damped_coulomb_kernel(eel, j, i,& 
                                                    2, kernel(1:3), dr)
                        
                        tmpE = 0.0_rp
                        call mu_elec_prop(ext_ipd(:,jpol), dr, kernel, .false., tmpV, &
                                        .true., tmpE, .false., tmpEgr, & 
                                        .false., tmpHE)
                        
                        E(:, ipol) = E(:, ipol) - tmpE * scalf
                    end do
                end do
            end if
            call time_pull("Removing screened ff")
            call time_push()
            deallocate(fmm_ipd%multipoles)
            deallocate(fmm_ipd%local_expansion)
            deallocate(fmm_ipd)
            call time_pull("Deallocation")
            call time_pull("Total")
        else
        
        !$omp parallel do default(shared) schedule(dynamic) &
        !$omp private(i,j,to_do,to_scale,scalf,idx,tmpV,tmpE,tmpEgr,tmpHE,kernel,dr)
        do j=1, eel%pol_atoms
            do i=1, eel%pol_atoms
                if(j == i) cycle
                !loop on target
                to_do = .true.
                to_scale = .false.
                scalf = 1.0

                ! Check if the element should be scaled
                do idx=eel%list_P_P%ri(i), eel%list_P_P%ri(i+1)-1
                    if(eel%list_P_P%ci(idx) == j) then
                        to_scale = .true.
                        exit
                    end if
                end do

                !If it should set the correct variables
                if(to_scale) then
                    to_do = eel%todo_P_P(idx)
                    scalf = eel%scalef_P_P(idx)
                end if
                
                if(to_do) then
                    call damped_coulomb_kernel(eel, eel%polar_mm(i), &
                                               eel%polar_mm(j),& 
                                               2, kernel(1:3), dr)
                    
                    tmpE = 0.0_rp

                    call mu_elec_prop(ext_ipd(:,i), dr, kernel, .false., tmpV, &
                                      .true., tmpE, .false., tmpEgr, & 
                                      .false., tmpHE)
                    if(to_scale) then
                        E(:, j) = E(:, j) + tmpE * scalf
                    else
                        E(:, j) = E(:, j) + tmpE
                    end if
                end if
            end do
        end do
        end if
    end subroutine field_extD2D
    
    subroutine elec_prop_D2D(eel, in_kind, do_V, do_E, do_Egrd, do_EHes)
        !! Computes the electric field of a trial set of induced point dipoles
        !! at polarizable sites. This is intended to be used as matrix-vector
        !! routine in the solution of the linear system.
        implicit none

        type(ommp_electrostatics_type), intent(inout) :: eel
        !! Data structure for electrostatic part of the system
        logical, intent(in) :: do_V, do_E, do_Egrd, do_EHes
        !! Flag to control which properties have to be computed.
        character, intent(in) :: in_kind

        integer(ip) :: i, j, jpol, ipol, ij, idx, ikernel, knd
        logical :: to_scale, to_do
        real(rp) :: kernel(5), dr(3), tmpV, tmpE(3), tmpEgr(6), tmpHE(10), scalf

        knd = 1 ! Default
        if(in_kind == 'P') then
            knd = _amoeba_P_
        elseif(in_kind == 'D') then 
            knd = _amoeba_D_
        elseif(eel%amoeba) then
            call fatal_error("Unrecognized interaction '"//in_kind//"' in elec&
                             &_prop_D2D.")
        else
            knd = 1
        end if

        if(do_EHes) then
            ikernel = 5 
        elseif(do_Egrd) then
            ikernel = 4
        elseif(do_E) then
            ikernel = 3
        elseif(do_V) then
            ikernel = 2
        else
            return
        end if

        if(eel%use_fmm) then
            call prepare_fmm_ipd(eel, knd)


            !$omp parallel do default(shared) schedule(dynamic) &
            !$omp private(i,j,ij,ipol,jpol,idx,dr,kernel,to_do,to_scale,scalf,tmpV,tmpE,tmpEgr,tmpHE) 
            do ipol=1, eel%pol_atoms 
                i = eel%polar_mm(ipol)
                
                call cart_propfar_at_ipart(eel%fmm_ipd(knd), i, &
                                           do_V, eel%V_D2D(ipol,knd), &
                                           do_E, eel%E_D2D(:,ipol,knd), &
                                           do_Egrd, eel%Egrd_D2D(:,ipol,knd), &
                                           do_EHes, eel%EHes_D2D(:,ipol,knd))

                ! Near field is computed internally because dumped kernel is required
                do ij=eel%fmm_near_field_list%ri(i), eel%fmm_near_field_list%ri(i+1)-1
                    j = eel%fmm_near_field_list%ci(ij)
                    jpol = eel%polar_mm(j)
                    ! If the atom is not polarizable, skip
                    if(jpol < 1) cycle 

                    !loop on target
                    to_do = .true.
                    to_scale = .false.
                    scalf = 1.0

                    ! Check if the element should be scaled
                    do idx=eel%list_P_P%ri(ipol), eel%list_P_P%ri(ipol+1)-1
                        if(eel%list_P_P%ci(idx) == jpol) then
                            to_scale = .true.
                            exit
                        end if
                    end do

                    !If it should set the correct variables
                    if(to_scale) then
                        to_do = eel%todo_P_P(idx)
                        scalf = eel%scalef_P_P(idx)
                    end if
                    
                    if(to_do) then
                        call damped_coulomb_kernel(eel, j, i,& 
                                                   ikernel, kernel, dr)
                        

                        if(do_V) tmpV = 0.0_rp
                        if(do_E) tmpE = 0.0_rp
                        if(do_Egrd) tmpEgr = 0.0_rp
                        if(do_EHes) tmpHE = 0.0_rp

                        call mu_elec_prop(eel%ipd(:,jpol,knd), dr, kernel, &
                                        do_V, tmpV, &
                                        do_E, tmpE, &
                                        do_Egrd, tmpEgr, & 
                                        do_EHes, tmpHE)
                        if(to_scale) then
                            if(do_V) eel%V_D2D(ipol,knd) = eel%V_D2D(ipol,knd) + tmpV * scalf
                            if(do_E) eel%E_D2D(:, ipol,knd) = eel%E_D2D(:, ipol,knd) + tmpE * scalf
                            if(do_Egrd) eel%Egrd_D2D(:, ipol,knd) = eel%Egrd_D2D(:, ipol,knd) + tmpEgr * scalf
                            if(do_EHes) eel%EHes_D2D(:, ipol,knd) = eel%EHes_D2D(:, ipol,knd) + tmpHE * scalf
                        else
                            if(do_V) eel%V_D2D(ipol,knd) = eel%V_D2D(ipol,knd) + tmpV
                            if(do_E) eel%E_D2D(:, ipol,knd) = eel%E_D2D(:, ipol,knd) + tmpE
                            if(do_Egrd) eel%Egrd_D2D(:, ipol,knd) = eel%Egrd_D2D(:, ipol,knd) + tmpEgr
                            if(do_EHes) eel%EHes_D2D(:, ipol,knd) = eel%EHes_D2D(:, ipol,knd) + tmpHE
                        end if
                    end if
                end do
            end do

            if(allocated(eel%list_P_P_fmm_far)) then
                ! Now remove screened interactions from far-field
                do ipol=1, eel%pol_atoms
                    i = eel%polar_mm(ipol)
                    
                    do idx=eel%list_P_P_fmm_far%ri(ipol), eel%list_P_P_fmm_far%ri(ipol+1)-1
                        jpol = eel%list_P_P_fmm_far%ci(idx)
                        j = eel%polar_mm(jpol)

                        scalf = 1.0 - eel%scalef_P_P_fmm_far(idx)
                        
                        call damped_coulomb_kernel(eel, j, i,& 
                                                    2, kernel(1:3), dr)
                        
                        if(do_V) tmpV = 0.0_rp
                        if(do_E) tmpE = 0.0_rp
                        if(do_Egrd) tmpEgr = 0.0_rp
                        if(do_EHes) tmpHE = 0.0_rp

                        call mu_elec_prop(eel%ipd(:,jpol,knd), dr, kernel, &
                                        do_V, tmpV, &
                                        do_E, tmpE, &
                                        do_Egrd, tmpEgr, & 
                                        do_EHes, tmpHE)

                        if(do_V) eel%V_D2D(ipol,knd) = eel%V_D2D(ipol,knd) + tmpV * scalf
                        if(do_E) eel%E_D2D(:, ipol,knd) = eel%E_D2D(:, ipol,knd) + tmpE * scalf
                        if(do_Egrd) eel%Egrd_D2D(:, ipol,knd) = eel%Egrd_D2D(:, ipol,knd) + tmpEgr * scalf
                        if(do_EHes) eel%EHes_D2D(:, ipol,knd) = eel%EHes_D2D(:, ipol,knd) + tmpHE * scalf
                    end do
                end do
            end if
        else
        !$omp parallel do default(shared) schedule(dynamic) &
        !$omp private(i,j,idx,dr,kernel,to_do,to_scale,scalf,tmpV,tmpE,tmpEgr,tmpHE) 
        do j=1, eel%pol_atoms
            do i=1, eel%pol_atoms
                if(j == i) cycle
                !loop on target
                to_do = .true.
                to_scale = .false.
                scalf = 1.0

                ! Check if the element should be scaled
                do idx=eel%list_P_P%ri(i), eel%list_P_P%ri(i+1)-1
                    if(eel%list_P_P%ci(idx) == j) then
                        to_scale = .true.
                        exit
                    end if
                end do

                !If it should set the correct variables
                if(to_scale) then
                    to_do = eel%todo_P_P(idx)
                    scalf = eel%scalef_P_P(idx)
                end if
                
                if(to_do) then
                    call damped_coulomb_kernel(eel, eel%polar_mm(i), &
                                               eel%polar_mm(j),& 
                                               ikernel, kernel, dr)
                    
                    if(do_V) tmpV = 0.0_rp
                    if(do_E) tmpE = 0.0_rp
                    if(do_Egrd) tmpEgr = 0.0_rp
                    if(do_EHes) tmpHE = 0.0_rp

                    call mu_elec_prop(eel%ipd(:,i,knd), dr, kernel, &
                                      do_V, tmpV, &
                                      do_E, tmpE, &
                                      do_Egrd, tmpEgr, & 
                                      do_EHes, tmpHE)
                    if(to_scale) then
                        if(do_V) eel%V_D2D(j,knd) = eel%V_D2D(j,knd) + tmpV * scalf
                        if(do_E) eel%E_D2D(:, j,knd) = eel%E_D2D(:, j,knd) + tmpE * scalf
                        if(do_Egrd) eel%Egrd_D2D(:, j,knd) = eel%Egrd_D2D(:, j,knd) + tmpEgr * scalf
                        if(do_EHes) eel%EHes_D2D(:, j,knd) = eel%EHes_D2D(:, j,knd) + tmpHE * scalf
                    else
                        if(do_V) eel%V_D2D(j,knd) = eel%V_D2D(j,knd) + tmpV
                        if(do_E) eel%E_D2D(:, j,knd) = eel%E_D2D(:, j,knd) + tmpE
                        if(do_Egrd) eel%Egrd_D2D(:, j,knd) = eel%Egrd_D2D(:, j,knd) + tmpEgr
                        if(do_EHes) eel%EHes_D2D(:, j,knd) = eel%EHes_D2D(:, j,knd) + tmpHE
                    end if
                end if
            end do
        end do
        end if
    end subroutine elec_prop_D2D
    
    subroutine elec_prop_M2D(eel, do_V, do_E, do_Egrd, do_EHes)
        !! Computes the electric field of static multipoles at induced dipoles
        !! sites. This is only intended to be used to build the RHS of the 
        !! linear system. This field is modified by the indroduction of the 
        !! damped kernels and by the connectivity-based screening rules.

        implicit none

        type(ommp_electrostatics_type), intent(inout) :: eel
        !! Electrostatics data structure
        logical, intent(in) :: do_V, do_E, do_Egrd, do_EHes
        !! Flag to control which properties have to be computed.

        integer(ip) :: i, ipol, j, jnode, ij, idx, ikernel
        logical :: to_do_p, to_scale_p, to_do_d, to_scale_d, to_do, to_scale, &
                   amoeba
        real(rp) :: kernel(5), dr(3), tmpV, tmpE(3), tmpEgr(6), tmpHE(10), &
                    scalf_p, scalf_d, scalf
        type(ommp_topology_type), pointer :: top
      
        write(*, *) 'M2D'
        call time_push
        ! Shortcuts
        top => eel%top
        amoeba = eel%amoeba

        if(do_EHes) then
            ikernel = 3 
        elseif(do_Egrd) then
            ikernel = 2
        elseif(do_E) then
            ikernel = 1
        elseif(do_V) then
            ikernel = 0
        else
            return
        end if
        if(eel%amoeba) ikernel = ikernel + 2_ip
        
        if(eel%use_fmm) then
            call preapare_fmm_static(eel)

            !$omp parallel do default(shared) schedule(dynamic) &
            !$omp private(i,j,ij,ipol,idx,dr,kernel,to_do_p,to_do_d,to_scale_p,to_scale_d,scalf_p,scalf_d,tmpV,tmpE,tmpEgr,tmpHE) 
            do ipol=1, eel%pol_atoms 
                i = eel%polar_mm(ipol)
                call cart_propfar_at_ipart(eel%fmm_static, i, &
                                           do_V, eel%V_M2D(ipol, _amoeba_D_), &
                                           do_E, eel%E_M2D(:, ipol, _amoeba_D_), &
                                           do_Egrd, eel%Egrd_M2D(:, ipol, _amoeba_D_), &
                                           do_EHes, eel%EHes_M2D(:, ipol, _amoeba_D_))
                if(eel%amoeba) then
                    if(do_V) eel%V_M2D(ipol, _amoeba_P_) = eel%V_M2D(ipol, _amoeba_D_)
                    if(do_E) eel%E_M2D(:, ipol, _amoeba_P_) = eel%E_M2D(:, ipol, _amoeba_D_)
                    if(do_Egrd) eel%Egrd_M2D(:, ipol, _amoeba_P_) = eel%Egrd_M2D(:, ipol, _amoeba_D_)
                    if(do_EHes) eel%EHes_M2D(:, ipol, _amoeba_P_) = eel%EHes_M2D(:, ipol, _amoeba_D_)
                end if
               
                ! Near field is computed internally because dumped kernel is required
                do ij=eel%fmm_near_field_list%ri(i), eel%fmm_near_field_list%ri(i+1)-1
                    j = eel%fmm_near_field_list%ci(ij)

                    !loop on target
                    to_do_p = .true.
                    to_scale_p = .false.
                    scalf_p = 1.0

                    ! Check if the element should be scaled
                    do idx=eel%list_S_P_P%ri(j), eel%list_S_P_P%ri(j+1)-1
                        if(eel%list_S_P_P%ci(idx) == ipol) then
                            to_scale_p = .true.
                            exit
                        end if
                    end do

                    !If it should set the correct variables
                    if(to_scale_p) then
                        to_do_p = eel%todo_S_P_P(idx)
                        scalf_p = eel%scalef_S_P_P(idx)
                    end if
                   
                    if(amoeba) then
                        to_do_d = .true.
                        to_scale_d = .false.
                        scalf_d = 1.0

                        ! Check if the element should be scaled
                        do idx=eel%list_S_P_D%ri(j), eel%list_S_P_D%ri(j+1)-1
                            if(eel%list_S_P_D%ci(idx) == ipol) then
                                to_scale_d = .true.
                                exit
                            end if
                        end do

                        !If it should set the correct variables
                        if(to_scale_d) then
                            to_do_d = eel%todo_S_P_D(idx)
                            scalf_d = eel%scalef_S_P_D(idx)
                        end if
                    else
                        to_do_d = .false.
                        to_scale_d = .true.
                        scalf_d = 0.0
                    end if

                    if(to_do_p .or. to_do_d) then
                    
                        call damped_coulomb_kernel(eel, j, i, &
                                                ikernel, kernel, dr)

                        if(do_V) tmpV = 0.0_rp
                        if(do_E) tmpE = 0.0_rp
                        if(do_Egrd) tmpEgr = 0.0_rp
                        if(do_EHes) tmpHE = 0.0_rp

                        call q_elec_prop(eel%q(1,j), dr, kernel, &
                                            do_V, tmpV, & 
                                            do_E, tmpE, &
                                            do_Egrd, tmpEgr, &
                                            do_EHes, tmpHE)
                        if(eel%amoeba) then
                            call mu_elec_prop(eel%q(2:4,j), dr, kernel, &
                                                do_V, tmpV, & 
                                                do_E, tmpE, &
                                                do_Egrd, tmpEgr, &
                                                do_EHes, tmpHE)

                            call quad_elec_prop(eel%q(5:10,j), dr, kernel, &
                                                do_V, tmpV, & 
                                                do_E, tmpE, &
                                                do_Egrd, tmpEgr, &
                                                do_EHes, tmpHE)
                        end if

                        if(eel%amoeba) then
                            if(to_scale_p) then
                                if(do_V) eel%V_M2D(ipol, _amoeba_P_) = eel%V_M2D(ipol, _amoeba_P_) + tmpV * scalf_p 
                                if(do_E) eel%E_M2D(:, ipol, _amoeba_P_) = eel%E_M2D(:, ipol, _amoeba_P_) + tmpE * scalf_p
                                if(do_Egrd) eel%Egrd_M2D(:, ipol, _amoeba_P_) = eel%Egrd_M2D(:, ipol, _amoeba_P_) + tmpEgr * scalf_p
                                if(do_EHes) eel%EHes_M2D(:, ipol, _amoeba_P_) = eel%EHes_M2D(:, ipol, _amoeba_P_) + tmpHE * scalf_p
                            else if(to_do_p) then
                                if(do_V) eel%V_M2D(ipol, _amoeba_P_) = eel%V_M2D(ipol, _amoeba_P_) + tmpV 
                                if(do_E) eel%E_M2D(:, ipol, _amoeba_P_) = eel%E_M2D(:, ipol, _amoeba_P_) + tmpE
                                if(do_Egrd) eel%Egrd_M2D(:, ipol, _amoeba_P_) = eel%Egrd_M2D(:, ipol, _amoeba_P_) + tmpEgr
                                if(do_EHes) eel%EHes_M2D(:, ipol, _amoeba_P_) = eel%EHes_M2D(:, ipol, _amoeba_P_) + tmpHE
                            end if
                        
                            if(to_scale_d) then
                                if(do_V) eel%V_M2D(ipol, _amoeba_D_) = eel%V_M2D(ipol, _amoeba_D_) + tmpV * scalf_d 
                                if(do_E) eel%E_M2D(:, ipol, _amoeba_D_) = eel%E_M2D(:, ipol, _amoeba_D_) + tmpE * scalf_d
                                if(do_Egrd) eel%Egrd_M2D(:, ipol, _amoeba_D_) = eel%Egrd_M2D(:, ipol, _amoeba_D_) + tmpEgr * scalf_d
                                if(do_EHes) eel%EHes_M2D(:, ipol, _amoeba_D_) = eel%EHes_M2D(:, ipol, _amoeba_D_) + tmpHE * scalf_d
                            else if(to_do_d) then
                                if(do_V) eel%V_M2D(ipol, _amoeba_D_) = eel%V_M2D(ipol, _amoeba_D_) + tmpV 
                                if(do_E) eel%E_M2D(:, ipol, _amoeba_D_) = eel%E_M2D(:, ipol, _amoeba_D_) + tmpE
                                if(do_Egrd) eel%Egrd_M2D(:, ipol, _amoeba_D_) = eel%Egrd_M2D(:, ipol, _amoeba_D_) + tmpEgr
                                if(do_EHes) eel%EHes_M2D(:, ipol, _amoeba_D_) = eel%EHes_M2D(:, ipol, _amoeba_D_) + tmpHE
                            end if
                        else
                            if(to_scale_p) then
                                if(do_V) eel%V_M2D(ipol, 1) = eel%V_M2D(ipol, 1) + tmpV * scalf_p 
                                if(do_E) eel%E_M2D(:, ipol, 1) = eel%E_M2D(:, ipol, 1) + tmpE * scalf_p
                                if(do_Egrd) eel%Egrd_M2D(:, ipol, 1) = eel%Egrd_M2D(:, ipol, 1) + tmpEgr * scalf_p
                                if(do_EHes) eel%EHes_M2D(:, ipol, 1) = eel%EHes_M2D(:, ipol, 1) + tmpHE * scalf_p
                            else if(to_do_p) then
                                if(do_V) eel%V_M2D(ipol, 1) = eel%V_M2D(ipol, 1) + tmpV 
                                if(do_E) eel%E_M2D(:, ipol, 1) = eel%E_M2D(:, ipol, 1) + tmpE
                                if(do_Egrd) eel%Egrd_M2D(:, ipol, 1) = eel%Egrd_M2D(:, ipol, 1) + tmpEgr
                                if(do_EHes) eel%EHes_M2D(:, ipol, 1) = eel%EHes_M2D(:, ipol, 1) + tmpHE
                            end if
                        end if

                    end if        
                end do
            end do
            
            if(allocated(eel%list_S_P_P_fmm_far)) then
                ! Now remove screened interactions from far-field
                do j=1, top%mm_atoms
                    to_scale_p = .false.

                    do idx=eel%list_S_P_P_fmm_far%ri(j), eel%list_S_P_P_fmm_far%ri(j+1)-1
                        ipol = eel%list_S_P_P_fmm_far%ci(idx)
                        i = eel%polar_mm(ipol)
                        to_scale_p = .true.
                    
                        ! The interaction should be corrected
                        call damped_coulomb_kernel(eel, j, i, &
                                                ikernel, kernel, dr)

                        if(do_V) tmpV = 0.0_rp
                        if(do_E) tmpE = 0.0_rp
                        if(do_Egrd) tmpEgr = 0.0_rp
                        if(do_EHes) tmpHE = 0.0_rp

                        call q_elec_prop(eel%q(1,j), dr, kernel, &
                                            do_V, tmpV, & 
                                            do_E, tmpE, &
                                            do_Egrd, tmpEgr, &
                                            do_EHes, tmpHE)
                        if(eel%amoeba) then
                            call mu_elec_prop(eel%q(2:4,j), dr, kernel, &
                                                do_V, tmpV, & 
                                                do_E, tmpE, &
                                                do_Egrd, tmpEgr, &
                                                do_EHes, tmpHE)

                            call quad_elec_prop(eel%q(5:10,j), dr, kernel, &
                                                do_V, tmpV, & 
                                                do_E, tmpE, &
                                                do_Egrd, tmpEgr, &
                                                do_EHes, tmpHE)
                        end if

                        if(do_V) eel%V_M2D(ipol, _amoeba_P_) = eel%V_M2D(ipol, _amoeba_P_) - tmpV * scalf_p 
                        if(do_E) eel%E_M2D(:, ipol, _amoeba_P_) = eel%E_M2D(:, ipol, _amoeba_P_) - tmpE * scalf_p
                        if(do_Egrd) eel%Egrd_M2D(:, ipol, _amoeba_P_) = eel%Egrd_M2D(:, ipol, _amoeba_P_) - tmpEgr * scalf_p
                        if(do_EHes) eel%EHes_M2D(:, ipol, _amoeba_P_) = eel%EHes_M2D(:, ipol, _amoeba_P_) - tmpHE * scalf_p
                    end do
                end do
            end if
           
            if(allocated(eel%list_S_P_D_fmm_far)) then
                ! Now remove screened interactions from far-field
                do j=1, top%mm_atoms
                    to_scale_d = .false.

                    do idx=eel%list_S_P_D_fmm_far%ri(j), eel%list_S_P_D_fmm_far%ri(j+1)-1
                        ipol = eel%list_S_P_D_fmm_far%ci(idx)
                        i = eel%polar_mm(ipol)
                        to_scale_d = .true.
                    
                        ! The interaction should be corrected
                        call damped_coulomb_kernel(eel, j, i, &
                                                ikernel, kernel, dr)

                        if(do_V) tmpV = 0.0_rp
                        if(do_E) tmpE = 0.0_rp
                        if(do_Egrd) tmpEgr = 0.0_rp
                        if(do_EHes) tmpHE = 0.0_rp

                        call q_elec_prop(eel%q(1,j), dr, kernel, &
                                            do_V, tmpV, & 
                                            do_E, tmpE, &
                                            do_Egrd, tmpEgr, &
                                            do_EHes, tmpHE)
                        if(eel%amoeba) then
                            call mu_elec_prop(eel%q(2:4,j), dr, kernel, &
                                                do_V, tmpV, & 
                                                do_E, tmpE, &
                                                do_Egrd, tmpEgr, &
                                                do_EHes, tmpHE)

                            call quad_elec_prop(eel%q(5:10,j), dr, kernel, &
                                                do_V, tmpV, & 
                                                do_E, tmpE, &
                                                do_Egrd, tmpEgr, &
                                                do_EHes, tmpHE)
                        end if

                        if(do_V) eel%V_M2D(ipol, _amoeba_D_) = eel%V_M2D(ipol, _amoeba_D_) - tmpV * scalf_d
                        if(do_E) eel%E_M2D(:, ipol, _amoeba_D_) = eel%E_M2D(:, ipol, _amoeba_D_) - tmpE * scalf_d
                        if(do_Egrd) eel%Egrd_M2D(:, ipol, _amoeba_D_) = eel%Egrd_M2D(:, ipol, _amoeba_D_) - tmpEgr * scalf_d
                        if(do_EHes) eel%EHes_M2D(:, ipol, _amoeba_D_) = eel%EHes_M2D(:, ipol, _amoeba_D_) - tmpHE * scalf_d
                    end do
                end do
            end if
        else
        if(amoeba) then
            !$omp parallel do default(shared) schedule(dynamic) &
            !$omp private(i,j,idx,dr,kernel,to_do_p,to_do_d,to_scale_p,to_scale_d,scalf_p,scalf_d,tmpV,tmpE,tmpEgr,tmpHE) 
            do j=1, eel%pol_atoms
                ! loop on sources
                do i=1, top%mm_atoms
                    if(eel%polar_mm(j) == i) cycle
                    !loop on target
                    to_do_p = .true.
                    to_scale_p = .false.
                    scalf_p = 1.0

                    ! Check if the element should be scaled
                    do idx=eel%list_S_P_P%ri(i), eel%list_S_P_P%ri(i+1)-1
                        if(eel%list_S_P_P%ci(idx) == j) then
                            to_scale_p = .true.
                            exit
                        end if
                    end do

                    !If it should set the correct variables
                    if(to_scale_p) then
                        to_do_p = eel%todo_S_P_P(idx)
                        scalf_p = eel%scalef_S_P_P(idx)
                    end if
                    
                    to_do_d = .true.
                    to_scale_d = .false.
                    scalf_d = 1.0

                    ! Check if the element should be scaled
                    do idx=eel%list_S_P_D%ri(i), eel%list_S_P_D%ri(i+1)-1
                        if(eel%list_S_P_D%ci(idx) == j) then
                            to_scale_d = .true.
                            exit
                        end if
                    end do

                    !!If it should set the correct variables
                    if(to_scale_d) then
                        to_do_d = eel%todo_S_P_D(idx)
                        scalf_d = eel%scalef_S_P_D(idx)
                    end if
                    
                    if(to_do_p .or. to_do_d) then
                        call damped_coulomb_kernel(eel, i, eel%polar_mm(j), &
                                                   ikernel, kernel, dr)
                        
                        if(do_V) tmpV = 0.0_rp
                        if(do_E) tmpE = 0.0_rp
                        if(do_Egrd) tmpEgr = 0.0_rp
                        if(do_EHes) tmpHE = 0.0_rp

                        call q_elec_prop(eel%q(1,i), dr, kernel, &
                                         do_V, tmpV, &
                                         do_E, tmpE, &
                                         do_Egrd, tmpEgr, & 
                                         do_EHes, tmpHE)
                        call mu_elec_prop(eel%q(2:4,i), dr, kernel, &
                                          do_V, tmpV, &
                                          do_E, tmpE, &
                                          do_Egrd, tmpEgr, & 
                                          do_EHes, tmpHE)
                        call quad_elec_prop(eel%q(5:10,i), dr, kernel, &
                                            do_V, tmpV, &
                                            do_E, tmpE, &
                                            do_Egrd, tmpEgr, & 
                                            do_EHes, tmpHE)

                        if(to_do_p) then
                            if(to_scale_p) then
                                if(do_V) eel%V_M2D(j, _amoeba_P_) = eel%V_M2D(j, _amoeba_P_) + tmpV * scalf_p
                                if(do_E) eel%E_M2D(:, j, _amoeba_P_) = eel%E_M2D(:, j, _amoeba_P_) + tmpE * scalf_p
                                if(do_Egrd) eel%Egrd_M2D(:, j, _amoeba_P_) = eel%Egrd_M2D(:, j, _amoeba_P_) + tmpEgr * scalf_p
                                if(do_EHes) eel%EHes_M2D(:, j, _amoeba_P_) = eel%EHes_M2D(:, j, _amoeba_P_) + tmpHE * scalf_p
                            else
                                if(do_V) eel%V_M2D(j, _amoeba_P_) = eel%V_M2D(j, _amoeba_P_) + tmpV 
                                if(do_E) eel%E_M2D(:, j, _amoeba_P_) = eel%E_M2D(:, j, _amoeba_P_) + tmpE
                                if(do_Egrd) eel%Egrd_M2D(:, j, _amoeba_P_) = eel%Egrd_M2D(:, j, _amoeba_P_) + tmpEgr
                                if(do_EHes) eel%EHes_M2D(:, j, _amoeba_P_) = eel%EHes_M2D(:, j, _amoeba_P_) + tmpHE
                            end if
                        end if

                        if(to_do_d) then
                            if(to_scale_d) then
                                if(do_V) eel%V_M2D(j, _amoeba_D_) = eel%V_M2D(j, _amoeba_D_) + tmpV * scalf_d
                                if(do_E) eel%E_M2D(:, j, _amoeba_D_) = eel%E_M2D(:, j, _amoeba_D_) + tmpE * scalf_d
                                if(do_Egrd) eel%Egrd_M2D(:, j, _amoeba_D_) = eel%Egrd_M2D(:, j, _amoeba_D_) + tmpEgr * scalf_d
                                if(do_EHes) eel%EHes_M2D(:, j, _amoeba_D_) = eel%EHes_M2D(:, j, _amoeba_D_) + tmpHE * scalf_d
                            else
                                if(do_V) eel%V_M2D(j, _amoeba_D_) = eel%V_M2D(j, _amoeba_D_) + tmpV 
                                if(do_E) eel%E_M2D(:, j, _amoeba_D_) = eel%E_M2D(:, j, _amoeba_D_) + tmpE
                                if(do_Egrd) eel%Egrd_M2D(:, j, _amoeba_D_) = eel%Egrd_M2D(:, j, _amoeba_D_) + tmpEgr
                                if(do_EHes) eel%EHes_M2D(:, j, _amoeba_D_) = eel%EHes_M2D(:, j, _amoeba_D_) + tmpHE
                            end if
                        end if
                    end if
                end do
            end do
        else
            !$omp parallel do default(shared) schedule(dynamic) &
            !$omp private(i,j,idx,dr,kernel,to_do,to_scale,scalf,tmpV,tmpE,tmpEgr,tmpHE) 
            do j=1, eel%pol_atoms
                ! loop on sources
                do i=1, top%mm_atoms
                    if(eel%polar_mm(j) == i) cycle
                    !loop on target
                    to_do = .true.
                    to_scale = .false.
                    scalf = 1.0

                    ! Check if the element should be scaled
                    do idx=eel%list_S_P_P%ri(i), eel%list_S_P_P%ri(i+1)-1
                        if(eel%list_S_P_P%ci(idx) == j) then
                            to_scale = .true.
                            exit
                        end if
                    end do

                    !If it should set the correct variables
                    if(to_scale) then
                        to_do = eel%todo_S_P_P(idx)
                        scalf = eel%scalef_S_P_P(idx)
                    end if

                    if(to_do) then
                        call damped_coulomb_kernel(eel, i, eel%polar_mm(j), & 
                                                   ikernel, kernel, dr)
                       
                        if(do_V) tmpV = 0.0_rp
                        if(do_E) tmpE = 0.0_rp
                        if(do_Egrd) tmpEgr = 0.0_rp
                        if(do_EHes) tmpHE = 0.0_rp
                        
                        call q_elec_prop(eel%q(1,i), dr, kernel, & 
                                         do_V, tmpV, &
                                         do_E, tmpE, &
                                         do_Egrd, tmpEgr, &
                                         do_EHes, tmpHE)
                        if(to_scale) then
                            if(do_V) eel%V_M2D(j, 1) = eel%V_M2D(j, 1) + tmpV * scalf
                            if(do_E) eel%E_M2D(:, j, 1) = eel%E_M2D(:, j, 1) + tmpE * scalf
                            if(do_Egrd) eel%Egrd_M2D(:, j, 1) = eel%Egrd_M2D(:, j, 1) + tmpEgr * scalf
                            if(do_EHes) eel%EHes_M2D(:, j, 1) = eel%EHes_M2D(:, j, 1) + tmpHE * scalf
                        else
                            if(do_V) eel%V_M2D(j, 1) = eel%V_M2D(j, 1) + tmpV
                            if(do_E) eel%E_M2D(:, j, 1) = eel%E_M2D(:, j, 1) + tmpE
                            if(do_Egrd) eel%Egrd_M2D(:, j, 1) = eel%Egrd_M2D(:, j, 1) + tmpEgr
                            if(do_EHes) eel%EHes_M2D(:, j, 1) = eel%EHes_M2D(:, j, 1) + tmpHE
                        end if 
                    end if
                end do
            end do
        end if
        end if
        call time_pull('M2D')
        write(*, *) 'M2D'
    end subroutine
    
    subroutine elec_prop_D2M(eel, in_kind, do_V, do_E, do_Egrd, do_EHes)

        implicit none

        type(ommp_electrostatics_type), intent(inout) :: eel
        !! Electrostatics data structure
        character, intent(in) :: in_kind

        logical, intent(in) :: do_V, do_E, do_Egrd, do_EHes
        !! Flag to control which properties have to be computed.

        integer(ip) :: i, j, ij, ipol, jpol, idx, ikernel, knd
        logical :: to_do, to_scale, amoeba
        real(rp) :: kernel(5), dr(3), tmpV, tmpE(3), tmpEgr(6), tmpHE(10), &
                    scalf
        type(ommp_topology_type), pointer :: top
        character :: screening_type
        
        ! Shortcuts
        top => eel%top
        amoeba = eel%amoeba
       
        knd =  1 ! Default
        if(eel%amoeba .and. in_kind /= 'P' .and. in_kind /= 'D') then
            call fatal_error("Unrecognized field '"//in_kind//"' for AMOEBA &
                             &force-field.")
        elseif(eel%amoeba .and. in_kind == 'P') then
            knd = _amoeba_P_
            screening_type = 'D'
        elseif(eel%amoeba .and. in_kind == 'D') then
            knd = _amoeba_D_
            screening_type = 'P'
        elseif(.not. eel%amoeba) then
            screening_type = '-'
        else
            call fatal_error("Unexpected error in elec_prop_D2M.")
        end if

        if(do_EHes) then
            ikernel = 4 
        elseif(do_Egrd) then
            ikernel = 3
        elseif(do_E) then
            ikernel = 2
        elseif(do_V) then 
            ikernel = 1
        else
            return
        end if
        
        if(eel%use_fmm) then
            
            call prepare_fmm_ipd(eel, knd)

            !$omp parallel do default(shared) schedule(dynamic) &
            !$omp private(i,j,ij,jpol,idx,dr,kernel,to_do,to_scale,scalf,tmpV,tmpE,tmpEgr,tmpHE) 
            do i=1, top%mm_atoms
                
                call cart_propfar_at_ipart(eel%fmm_ipd(knd), i, &
                                           do_V, eel%V_D2M(i), &
                                           do_E, eel%E_D2M(:,i), &
                                           do_Egrd, eel%Egrd_D2M(:,i), &
                                           do_EHes, eel%EHes_D2M(:,i))

                ! Near field is computed internally because dumped kernel is required
                do ij=eel%fmm_near_field_list%ri(i), eel%fmm_near_field_list%ri(i+1)-1
                    j = eel%fmm_near_field_list%ci(ij)
                    jpol = eel%polar_mm(j)
                    ! If the atom is not polarizable, skip
                    if(jpol < 1) cycle 

                    if(screening_type == 'P') then
                        to_do = .true.
                        to_scale = .false.
                        scalf = 1.0

                        ! Check if the element should be scaled
                        do idx=eel%list_S_P_P%ri(i), eel%list_S_P_P%ri(i+1)-1
                            if(eel%list_S_P_P%ci(idx) == jpol) then
                                to_scale = .true.
                                exit
                            end if
                        end do

                        !If it should set the correct variables
                        if(to_scale) then
                            to_do = eel%todo_S_P_P(idx)
                            scalf = eel%scalef_S_P_P(idx)
                        end if
                    else
                        to_do = .true.
                        to_scale = .false.
                        scalf = 1.0

                        ! Check if the element should be scaled
                        do idx=eel%list_S_P_D%ri(i), eel%list_S_P_D%ri(i+1)-1
                            if(eel%list_S_P_D%ci(idx) == jpol) then
                                to_scale = .true.
                                exit
                            end if
                        end do

                        !If it should set the correct variables
                        if(to_scale) then
                            to_do = eel%todo_S_P_D(idx)
                            scalf = eel%scalef_S_P_D(idx)
                        end if
                    end if

                    
                    if(to_do) then
                        call damped_coulomb_kernel(eel, j, i, & 
                                                   ikernel, kernel, dr)
                       
                        if(do_V) tmpV = 0.0_rp
                        if(do_E) tmpE = 0.0_rp
                        if(do_Egrd) tmpEgr = 0.0_rp
                        if(do_EHes) tmpHE = 0.0_rp
                        
                        call mu_elec_prop(eel%ipd(:,jpol, knd), dr, kernel, & 
                                         do_V, tmpV, &
                                         do_E, tmpE, &
                                         do_Egrd, tmpEgr, &
                                         do_EHes, tmpHE)
                        if(to_scale) then
                            if(do_V) eel%V_D2M(i) = eel%V_D2M(i) + tmpV * scalf
                            if(do_E) eel%E_D2M(:, i) = eel%E_D2M(:, i) + tmpE * scalf
                            if(do_Egrd) eel%Egrd_D2M(:, i) = eel%Egrd_D2M(:, i) + tmpEgr * scalf
                            if(do_EHes) eel%EHes_D2M(:, i) = eel%EHes_D2M(:, i) + tmpHE * scalf
                        else
                            if(do_V) eel%V_D2M(i) = eel%V_D2M(i) + tmpV
                            if(do_E) eel%E_D2M(:, i) = eel%E_D2M(:, i) + tmpE
                            if(do_Egrd) eel%Egrd_D2M(:, i) = eel%Egrd_D2M(:, i) + tmpEgr
                            if(do_EHes) eel%EHes_D2M(:, i) = eel%EHes_D2M(:, i) + tmpHE
                        end if 
                    end if
                end do
            end do

            ! Now remove screened interactions from far-field, hopefully they should be almost absent
            if(screening_type == 'P' .and. allocated(eel%list_S_P_P_fmm_far)) then
                do i=1, top%mm_atoms
                    do idx=eel%list_S_P_P_fmm_far%ri(i), eel%list_S_P_P_fmm_far%ri(i+1)-1
                        jpol = eel%list_S_P_P_fmm_far%ci(idx)
                        j = eel%polar_mm(jpol)

                        scalf = 1.0 - eel%scalef_S_P_P_fmm_far(idx)
                        
                        call damped_coulomb_kernel(eel, j, i,& 
                                                    ikernel, kernel, dr)
                        
                        if(do_V) tmpV = 0.0_rp
                        if(do_E) tmpE = 0.0_rp
                        if(do_Egrd) tmpEgr = 0.0_rp
                        if(do_EHes) tmpHE = 0.0_rp

                        call mu_elec_prop(eel%ipd(:,jpol,knd), dr, kernel, &
                                        do_V, tmpV, &
                                        do_E, tmpE, &
                                        do_Egrd, tmpEgr, & 
                                        do_EHes, tmpHE)

                        if(do_V) eel%V_D2M(i) = eel%V_D2M(i) - tmpV * scalf
                        if(do_E) eel%E_D2M(:, i) = eel%E_D2M(:, i) - tmpE * scalf
                        if(do_Egrd) eel%Egrd_D2M(:, i) = eel%Egrd_D2M(:, i) - tmpEgr * scalf
                        if(do_EHes) eel%EHes_D2M(:, i) = eel%EHes_D2M(:, i) - tmpHE * scalf
                    end do
                end do
            else if(screening_type == 'D'.and. allocated(eel%scalef_S_P_D_fmm_far)) then
                do i=1, top%mm_atoms
                    do idx=eel%list_S_P_D_fmm_far%ri(i), eel%list_S_P_D_fmm_far%ri(i+1)-1
                        jpol = eel%list_S_P_D_fmm_far%ci(idx)
                        j = eel%polar_mm(jpol)

                        scalf = 1.0 - eel%scalef_S_P_D_fmm_far(idx)
                        
                        call damped_coulomb_kernel(eel, j, i,& 
                                                    ikernel, kernel, dr)
                        
                        if(do_V) tmpV = 0.0_rp
                        if(do_E) tmpE = 0.0_rp
                        if(do_Egrd) tmpEgr = 0.0_rp
                        if(do_EHes) tmpHE = 0.0_rp

                        call mu_elec_prop(eel%ipd(:,jpol,knd), dr, kernel, &
                                        do_V, tmpV, &
                                        do_E, tmpE, &
                                        do_Egrd, tmpEgr, & 
                                        do_EHes, tmpHE)

                        if(do_V) eel%V_D2M(i) = eel%V_D2M(i) - tmpV * scalf
                        if(do_E) eel%E_D2M(:, i) = eel%E_D2M(:, i) - tmpE * scalf
                        if(do_Egrd) eel%Egrd_D2M(:, i) = eel%Egrd_D2M(:, i) - tmpEgr * scalf
                        if(do_EHes) eel%EHes_D2M(:, i) = eel%EHes_D2M(:, i) - tmpHE * scalf
                    end do
                end do
            end if
        else
        if(amoeba) then
            !$omp parallel do default(shared) schedule(dynamic) &
            !$omp private(i,j,idx,dr,kernel,to_do,to_scale,scalf,tmpV,tmpE,tmpEgr,tmpHE) 
            do j=1, top%mm_atoms
                ! loop on sources
                do i=1, eel%pol_atoms
                    if(eel%polar_mm(i) == j) cycle
                    !loop on target
                    if(screening_type == 'P') then
                        to_do = .true.
                        to_scale = .false.
                        scalf = 1.0

                        ! Check if the element should be scaled
                        do idx=eel%list_S_P_P%ri(j), eel%list_S_P_P%ri(j+1)-1
                            if(eel%list_S_P_P%ci(idx) == i) then
                                to_scale = .true.
                                exit
                            end if
                        end do

                        !If it should set the correct variables
                        if(to_scale) then
                            to_do = eel%todo_S_P_P(idx)
                            scalf = eel%scalef_S_P_P(idx)
                        end if
                    else
                        to_do = .true.
                        to_scale = .false.
                        scalf = 1.0

                        ! Check if the element should be scaled
                        do idx=eel%list_S_P_D%ri(j), eel%list_S_P_D%ri(j+1)-1
                            if(eel%list_S_P_D%ci(idx) == i) then
                                to_scale = .true.
                                exit
                            end if
                        end do

                        !If it should set the correct variables
                        if(to_scale) then
                            to_do = eel%todo_S_P_D(idx)
                            scalf = eel%scalef_S_P_D(idx)
                        end if
                    end if

                    
                    if(to_do) then
                        call damped_coulomb_kernel(eel, eel%polar_mm(i), j, & 
                                                   ikernel, kernel, dr)
                       
                        if(do_V) tmpV = 0.0_rp
                        if(do_E) tmpE = 0.0_rp
                        if(do_Egrd) tmpEgr = 0.0_rp
                        if(do_EHes) tmpHE = 0.0_rp
                        
                        call mu_elec_prop(eel%ipd(:,i, knd), dr, kernel, & 
                                         do_V, tmpV, &
                                         do_E, tmpE, &
                                         do_Egrd, tmpEgr, &
                                         do_EHes, tmpHE)
                        if(to_scale) then
                            if(do_V) eel%V_D2M(j) = eel%V_D2M(j) + tmpV * scalf
                            if(do_E) eel%E_D2M(:, j) = eel%E_D2M(:, j) + tmpE * scalf
                            if(do_Egrd) eel%Egrd_D2M(:, j) = eel%Egrd_D2M(:, j) + tmpEgr * scalf
                            if(do_EHes) eel%EHes_D2M(:, j) = eel%EHes_D2M(:, j) + tmpHE * scalf
                        else
                            if(do_V) eel%V_D2M(j) = eel%V_D2M(j) + tmpV
                            if(do_E) eel%E_D2M(:, j) = eel%E_D2M(:, j) + tmpE
                            if(do_Egrd) eel%Egrd_D2M(:, j) = eel%Egrd_D2M(:, j) + tmpEgr
                            if(do_EHes) eel%EHes_D2M(:, j) = eel%EHes_D2M(:, j) + tmpHE
                        end if 
                    end if
                end do
            end do
        else
            !$omp parallel do default(shared) schedule(dynamic) &
            !$omp private(i,j,idx,dr,kernel,to_do,to_scale,scalf,tmpV,tmpE,tmpEgr,tmpHE) 
            do j=1, top%mm_atoms
                ! loop on sources
                do i=1, eel%pol_atoms
                    if(eel%polar_mm(i) == j) cycle
                    !loop on target
                    to_do = .true.
                    to_scale = .false.
                    scalf = 1.0

                    ! Check if the element should be scaled
                    do idx=eel%list_S_P_P%ri(j), eel%list_S_P_P%ri(j+1)-1
                        if(eel%list_S_P_P%ci(idx) == i) then
                            to_scale = .true.
                            exit
                        end if
                    end do

                    !If it should set the correct variables
                    if(to_scale) then
                        to_do = eel%todo_S_P_P(idx)
                        scalf = eel%scalef_S_P_P(idx)
                    end if
                    
                    if(to_do) then
                        call damped_coulomb_kernel(eel, eel%polar_mm(i), j, & 
                                                   ikernel, kernel, dr)
                       
                        if(do_V) tmpV = 0.0_rp
                        if(do_E) tmpE = 0.0_rp
                        if(do_Egrd) tmpEgr = 0.0_rp
                        if(do_EHes) tmpHE = 0.0_rp
                        
                        call mu_elec_prop(eel%ipd(:,i, knd), dr, kernel, & 
                                         do_V, tmpV, &
                                         do_E, tmpE, &
                                         do_Egrd, tmpEgr, &
                                         do_EHes, tmpHE)
                        if(to_scale) then
                            if(do_V) eel%V_D2M(j) = eel%V_D2M(j) + tmpV * scalf
                            if(do_E) eel%E_D2M(:, j) = eel%E_D2M(:, j) + tmpE * scalf
                            if(do_Egrd) eel%Egrd_D2M(:, j) = eel%Egrd_D2M(:, j) + tmpEgr * scalf
                            if(do_EHes) eel%EHes_D2M(:, j) = eel%EHes_D2M(:, j) + tmpHE * scalf
                        else
                            if(do_V) eel%V_D2M(j) = eel%V_D2M(j) + tmpV
                            if(do_E) eel%E_D2M(:, j) = eel%E_D2M(:, j) + tmpE
                            if(do_Egrd) eel%Egrd_D2M(:, j) = eel%Egrd_D2M(:, j) + tmpEgr
                            if(do_EHes) eel%EHes_D2M(:, j) = eel%EHes_D2M(:, j) + tmpHE
                        end if 
                    end if
                end do
            end do
        end if
        end if
    end subroutine

    subroutine potential_D2E(eel, cpt, V, amoeba_P_insted_of_D_)
        !! This subroutine computes the potential generated by the induced 
        !! point dipoles to a set of arbitrary coordinates, without applying
        !! any screening rules. Note: for AMOEBA D dipoles should be used. 
        
        implicit none

        type(ommp_electrostatics_type), intent(in) :: eel
        !! Electrostatics data structure
        real(rp), intent(inout) :: V(:)
        !! Electric field (results will be added)
        real(rp), intent(in) :: cpt(:,:)
        !! Coordinates at which the electric field is requested
        logical, optional, intent(in) :: amoeba_P_insted_of_D_
        !! For AMOEBA FF, if true the potential of P dipoles
        !! is computed, otherwise potential of D dipoles is computed

        integer(ip) :: i, j, n_cpt
        logical :: amoeba_P_insted_of_D
        real(rp) :: kernel(5), dr(3), tmpV, tmpE(3), tmpEgr(6), tmpHE(10)

        if(eel%pol_atoms < 1) return

        if(present(amoeba_P_insted_of_D_)) then
            amoeba_P_insted_of_D = amoeba_P_insted_of_D_
        else
            amoeba_P_insted_of_D = .false.
        end if

        if(.not. eel%ipd_done) call fatal_error("IPD should be computed before&
                                                & D2E potential.")
        n_cpt = size(cpt, 2)

        if(eel%amoeba) then
            if(.not. amoeba_P_insted_of_D) then
                !$omp parallel do default(shared) schedule(dynamic) collapse(2) &
                !$omp private(i,j,dr,kernel,tmpV,tmpE,tmpEgr,tmpHE) reduction(+:V)
                do i=1, eel%pol_atoms
                    do j=1, n_cpt
                        dr = cpt(:,j) - eel%cpol(:,i)
                        call coulomb_kernel(dr, 2, kernel(1:2))
                        tmpV = 0.0_rp
                        
                        call mu_elec_prop(eel%ipd(:,i,_amoeba_D_), &
                                        dr, kernel, .true., tmpV, &
                                        .false., tmpE, .false., tmpEgr, & 
                                        .false., tmpHE)

                        V(j) = V(j) + tmpV
                    end do
                end do
            else
                !$omp parallel do default(shared) schedule(dynamic) collapse(2) &
                !$omp private(i,j,dr,kernel,tmpV,tmpE,tmpEgr,tmpHE) reduction(+:V)
                do i=1, eel%pol_atoms
                    do j=1, n_cpt
                        dr = cpt(:,j) - eel%cpol(:,i)
                        call coulomb_kernel(dr, 2, kernel(1:2))
                        tmpV = 0.0_rp
                        
                        call mu_elec_prop(eel%ipd(:,i,_amoeba_P_), &
                                        dr, kernel, .true., tmpV, &
                                        .false., tmpE, .false., tmpEgr, & 
                                        .false., tmpHE)

                        V(j) = V(j) + tmpV
                    end do
                end do
            end if
        else
            !$omp parallel do default(shared) schedule(dynamic) collapse(2) &
            !$omp private(i,j,dr,kernel,tmpV,tmpE,tmpEgr,tmpHE) reduction(+:V)
            do i=1, eel%pol_atoms
                ! loop on sources
                do j=1, n_cpt
                    dr = cpt(:,j) - eel%cpol(:,i)
                    call coulomb_kernel(dr, 2, kernel(1:2))
                    tmpV = 0.0_rp
                    
                    call mu_elec_prop(eel%ipd(:,i,1), &
                                      dr, kernel, .true., tmpV, &
                                      .false., tmpE, .false., tmpEgr, & 
                                      .false., tmpHE)
                    
                    V(j) = V(j) + tmpV
                end do
            end do
        end if
    end subroutine potential_D2E

    subroutine potential_M2E(eel, cpt, V)
        !! This subroutine computes the potential generated by the static
        !! multipoles to a set of arbitrary coordinates, without applying
        !! any screening rules.
        
        implicit none

        type(ommp_electrostatics_type), intent(in) :: eel
        !! Electrostatics data structure
        real(rp), intent(inout) :: V(:)
        !! Electric field (results will be added)
        real(rp), intent(in) :: cpt(:,:)
        !! Coordinates at which the electric field is requested

        integer(ip) :: i, j, n_cpt
        real(rp) :: kernel(5), dr(3), tmpV, tmpE(3), tmpEgr(6), tmpHE(10)

        n_cpt = size(cpt, 2)

        if(eel%amoeba) then
            !$omp parallel do default(shared) schedule(dynamic) collapse(2) &
            !$omp private(i,j,dr,kernel,tmpV,tmpE,tmpEgr,tmpHE) reduction(+:V)
            do i=1, eel%top%mm_atoms
                do j=1, n_cpt
                    dr = cpt(:,j) - eel%top%cmm(:,i)
                    call coulomb_kernel(dr, 2, kernel(1:3))
                    tmpV = 0.0_rp
                    
                    call q_elec_prop(eel%q(1,i), dr, kernel, .true., tmpV, &
                                     .false., tmpE, .false., tmpEgr, & 
                                     .false., tmpHE)
                    call mu_elec_prop(eel%q(2:4,i), dr, kernel, .true., tmpV, &
                                      .false., tmpE, .false., tmpEgr, & 
                                      .false., tmpHE)
                    call quad_elec_prop(eel%q(5:10,i), dr, kernel, .true., tmpV, &
                                        .false., tmpE, .false., tmpEgr, & 
                                        .false., tmpHE)

                    V(j) = V(j) + tmpV
                end do
            end do
        else
            !$omp parallel do default(shared) schedule(dynamic) collapse(2) &
            !$omp private(i,j,dr,kernel,tmpV,tmpE,tmpEgr,tmpHE) reduction(+:V)
            do i=1, eel%top%mm_atoms
                ! loop on sources
                do j=1, n_cpt
                    dr = cpt(:,j) - eel%top%cmm(:,i)
                    call coulomb_kernel(dr, 1, kernel(1:2))
                    tmpV = 0.0_rp
                    
                    call q_elec_prop(eel%q(1,i), dr, kernel, .true., tmpV, &
                                     .false., tmpE, .false., tmpEgr, & 
                                     .false., tmpHE)
                    
                    V(j) = V(j) + tmpV
                end do
            end do
        end if
    end subroutine potential_M2E
    
    subroutine field_D2E(eel, cpt, E)
        !! This subroutine computes the potential generated by the static
        !! multipoles to a set of arbitrary coordinates, without applying
        !! any screening rules.
        
        implicit none

        type(ommp_electrostatics_type), intent(in) :: eel
        !! Electrostatics data structure
        real(rp), intent(inout) :: E(:,:)
        !! Electric field (results will be added)
        real(rp), intent(in) :: cpt(:,:)
        !! Coordinates at which the electric field is requested

        integer(ip) :: i, j, n_cpt
        real(rp) :: kernel(5), dr(3), tmpV, tmpE(3), tmpEgr(6), tmpHE(10)

        if(eel%pol_atoms < 1) return

        if(.not. eel%ipd_done) call fatal_error("IPD should be computed before&
                                                & computing D2E field.")
        n_cpt = size(cpt, 2)

        if(eel%amoeba) then
            !$omp parallel do default(shared) schedule(dynamic) &
            !$omp private(i,j,dr,kernel,tmpV,tmpE,tmpEgr,tmpHE) reduction(+:E)
            do i=1, eel%pol_atoms
                do j=1, n_cpt
                    dr = cpt(:,j) - eel%cpol(:,i)
                    call coulomb_kernel(dr, 3, kernel)
                    tmpE = 0.0_rp
                    !TODO
                    call mu_elec_prop(0.5*(eel%ipd(:,i, _amoeba_P_) + eel%ipd(:,i, _amoeba_D_)), dr, kernel, .false., tmpV, &
                                      .true., tmpE, .false., tmpEgr, & 
                                      .false., tmpHE)

                    E(:,j) = E(:,j) + tmpE 
                end do
            end do
        else
            !$omp parallel do default(shared) schedule(dynamic) &
            !$omp private(i,j,dr,kernel,tmpV,tmpE,tmpEgr,tmpHE) reduction(+:E)
            do i=1, eel%pol_atoms
                ! loop on sources
                do j=1, n_cpt
                    dr = cpt(:,j) - eel%cpol(:,i)
                    call coulomb_kernel(dr, 3, kernel)
                    tmpE = 0.0_rp
                    
                    call mu_elec_prop(eel%ipd(:,i,1), dr, kernel, .false., tmpV, &
                                     .true., tmpE, .false., tmpEgr, & 
                                     .false., tmpHE)
                    
                    E(:,j) = E(:,j) + tmpE
                end do
            end do
        end if
    end subroutine field_D2E

    subroutine field_M2E(eel, cpt, E)
        !! This subroutine computes the potential generated by the static
        !! multipoles to a set of arbitrary coordinates, without applying
        !! any screening rules.
        
        implicit none

        type(ommp_electrostatics_type), intent(in) :: eel
        !! Electrostatics data structure
        real(rp), intent(inout) :: E(:,:)
        !! Electric field (results will be added)
        real(rp), intent(in) :: cpt(:,:)
        !! Coordinates at which the electric field is requested

        integer(ip) :: i, j, n_cpt
        real(rp) :: kernel(5), dr(3), tmpV, tmpE(3), tmpEgr(6), tmpHE(10)

        n_cpt = size(cpt, 2)

        if(eel%amoeba) then
            !$omp parallel do default(shared) schedule(dynamic) &
            !$omp private(i,j,dr,kernel,tmpE,tmpV,tmpEgr,tmpHE) reduction(+:E)
            do i=1, eel%top%mm_atoms
                do j=1, n_cpt
                    dr = cpt(:,j) - eel%top%cmm(:,i)
                    call coulomb_kernel(dr, 4, kernel)
                    tmpE = 0.0_rp
                    
                    call q_elec_prop(eel%q(1,i), dr, kernel, .false., tmpV, &
                                     .true., tmpE, .false., tmpEgr, & 
                                     .false., tmpHE)
                    call mu_elec_prop(eel%q(2:4,i), dr, kernel, .false., tmpV, &
                                      .true., tmpE, .false., tmpEgr, & 
                                      .false., tmpHE)
                    call quad_elec_prop(eel%q(5:10,i), dr, kernel, .false., tmpV, &
                                        .true., tmpE, .false., tmpEgr, & 
                                        .false., tmpHE)

                    E(:,j) = E(:,j) + tmpE
                end do
            end do
        else
            !$omp parallel do default(shared) schedule(dynamic) &
            !$omp private(i,j,dr,kernel,tmpE,tmpV,tmpEgr,tmpHE) reduction(+:E)
            do i=1, eel%top%mm_atoms
                ! loop on sources
                do j=1, n_cpt
                    dr = cpt(:,j) - eel%top%cmm(:,i)
                    call coulomb_kernel(dr, 2, kernel)
                    tmpE = 0.0_rp
                    
                    call q_elec_prop(eel%q(1,i), dr, kernel, .false., tmpV, &
                                     .true., tmpE, .false., tmpEgr, & 
                                     .false., tmpHE)
                    
                    E(:,j) = E(:,j) + tmpE
                end do
            end do
        end if
    end subroutine field_M2E
    
    function screening_rules(eel, i, kind_i, j, kind_j, in_field) result(scalf)
        !! Utility function used to decide if an interaction between sites i and j
        !! should be computed and eventually scaled by a factor.
        !! This function is intended to be used in \(\mathcalO(n^2)\) code, for
        !! linear scaling code lists should be built.
        !! This is written to minimize code repetitions, all the screening rules
        !! are handled in two possible cases: 
        !! 1. rules based on adjacency matrix
        !! 2. rules based on AMOEBA polarization groups

        use mod_constants, only: eps_rp

        type(ommp_electrostatics_type), intent(in) :: eel
        !! Electrostatics object
        integer(ip), intent(in) :: i
        !! Index of source site (MM index is used for static sites, while
        !! Pol index is used for polarizable sites)
        integer(ip), intent(in) :: j
        !! Index of target site (MM index is used for static sites, while
        !! Pol index is used for polarizable sites)
        character, intent(in) :: in_field
        !! Which screening rules have to be applied? 'D' = screening rules
        !! for direct field; 'P' = screening rules for polarization field
        character, intent(in) :: kind_i, kind_j
        !! Type of sites i and j in the interaction for which the screening
        !! rules are required; possible choices are 'P' (polarizable site) or 
        !! 'S' (static site). Any other option will cause a fatal error.
        real(rp) :: scalf
        !! Scale factor for the interaction

        integer(ip) :: ineigh, grp, igrp, pg_i
        integer(ip) :: j_mm, i_mm
        character :: field, interaction
        logical :: amoeba
        type(ommp_topology_type), pointer :: top 
        real(rp) :: myscale(4), myscale_intra(4)

        ! Some shortcut
        top => eel%top
        amoeba = eel%amoeba

        ! Decide which kind of rule should be used
        if(kind_i == 'P') then
            i_mm = eel%polar_mm(i)
        else if(kind_i == 'S') then
            i_mm = i
        else
            call fatal_error('Unexpected value of kind_i in screening_rules')
            i_mm = 0
        end if

        if(kind_j == 'P') then
            j_mm = eel%polar_mm(j)
        else if(kind_j == 'S') then
            j_mm = j
        else
            call fatal_error('Unexpected value of kind_j in screening_rules')
            j_mm = 0
        end if

        myscale = 0.0_rp
        if(kind_j == 'P' .and. kind_i == 'P') then
            ! Use IPD-IPD screening rules
            myscale = eel%uscale
            interaction = 'P' !Pol-pol interaction is named P
        else if(kind_j == 'S' .and. kind_i == 'S') then
            ! Use static multipoles-static multipoles screening rules
            myscale = eel%mscale
            field = '-'
            interaction = 'S' !Static-static interaction is named S
        else
            ! Use static multipoles-IPD screening rules
            if(in_field == 'P' .and. amoeba) then
                myscale = eel%pscale
                myscale_intra = eel%pscale_intra
                field = 'P'
            else if(in_field == 'D' .and. amoeba) then
                myscale = eel%dscale
                field = 'D'
            else if(.not. amoeba) then
                myscale = eel%pscale
                field = 'P'
            else
                call fatal_error('Unexpected value of field in screening rules')
            end if
            interaction = 'M' !Mixed interaction is named M
        end if

        ! Default return value
        scalf = 1.0_rp
        
        if((.not. amoeba) .or. &
           (amoeba .and. interaction == 'M' .and. field == 'P') .or. &
           (amoeba .and. interaction == 'S')) then
            ! Screening rules based on connectivity matrix: 
            ! 1. all the ones of non-amoeba FF, 
            ! 2. the static to IPD polarization field
            ! 3. the static-static interaction of AMOEBA
            do ineigh=1,4
                ! Look if j is at distance ineigh from i
                if(any(top%conn(ineigh)%ci(top%conn(ineigh)%ri(i_mm): &
                                           top%conn(ineigh)%ri(i_mm+1)-1) == j_mm)) then
                   
                    if(interaction == 'M' .and. amoeba) then
                        ! Amoeba uses two different screening rules for atoms 
                        ! of the same group and for atoms of different group
                        if(eel%mmat_polgrp(i_mm) == eel%mmat_polgrp(j_mm)) then
                            scalf = myscale_intra(ineigh)
                        else
                            scalf = myscale(ineigh)
                        end if
                    else
                        ! In any other case just use the scale factor for this
                        ! distance
                        scalf = myscale(ineigh)
                    end if
                    
                    ! Exit the loop
                    exit 

                end if
            end do
        else if((amoeba .and. interaction == 'M' .and. field == 'D') .or. &
                (amoeba .and. interaction == 'P')) then
            ! Screening rules based on polarization groups: 
            ! 1. the static to IPD direct field
            ! 2. the IPD to IPD interaction of AMOEBA
            pg_i = eel%mmat_polgrp(i_mm)

            outer: do ineigh=1,4
                do igrp=eel%pg_conn(ineigh)%ri(pg_i), &
                        eel%pg_conn(ineigh)%ri(pg_i+1)-1
                    ! Indexes of groups at distance ineigh from group pg_i
                    grp = eel%pg_conn(ineigh)%ci(igrp)
                    
                    if(any(eel%polgrp_mmat%ci(eel%polgrp_mmat%ri(grp): &
                                              eel%polgrp_mmat%ri(grp+1)-1) == j_mm)) then
                        ! If atom j is in a group at distance ineigh from the 
                        ! one of atom i, the field is scaled according to dscale
                        scalf = myscale(ineigh)
                        exit outer
                    end if
                end do
            end do outer
        else 
            ! Unexpected error
            call fatal_error("Unexpected combination of parameter for screening_rules")
        end if

    end function screening_rules

    subroutine fmm_coordinates_update(eel)
        use mod_constants, only: angstrom2au
        implicit none

        type(ommp_electrostatics_type), intent(inout) :: eel
        integer(ip) :: i
       
        call time_push()
        call free_tree(eel%tree)
        call init_as_octatree(eel%tree, eel%top%cmm, 12.0_rp)
        !call init_as_ribtree(eel%tree, eel%top%cmm, 12.0_rp)
        call fmm_make_neigh_list(eel)
        call time_pull("Tree initialization")
        
        call time_push()
        call fmm_init(eel%fmm_static, eel%fmm_maxl, eel%tree)
        do i=1, eel%n_ipd
            call fmm_init(eel%fmm_ipd(i), eel%fmm_maxl, eel%tree)
        end do
        call time_pull("FMM initialization")
    end subroutine

    pure function fmm_list_are_near(eel, i, j) result(near)
        ! Check if i and j (mm atoms) are marked as near in fmm lists
        type(ommp_electrostatics_type), intent(in) :: eel
        integer(ip), intent(in) :: i, j
        logical :: near
        integer(ip) :: idx
        
        near = (i == j)
        if(near) return
        
        do idx=eel%fmm_near_field_list%ri(i), eel%fmm_near_field_list%ri(i+1)-1
            if(eel%fmm_near_field_list%ci(idx) == j) then
                near = .true.
                exit
            end if
        end do
        return
    end function

    subroutine fmm_make_neigh_list(eel)
        use mod_memory, only: mallocate
        use mod_adjacency_mat, only: reallocate_mat
        implicit none

        type(ommp_electrostatics_type), intent(inout) :: eel
        integer(ip), parameter :: n_guess_neigh = 30
        integer(ip) :: i, j, inode, part_idx, jnode, node_idx

        eel%fmm_near_field_list%n = eel%top%mm_atoms
        call mallocate('fmm_make_neigh_list [ri]', eel%top%mm_atoms+1, eel%fmm_near_field_list%ri)
        call mallocate('fmm_make_neigh_list [ci]', eel%top%mm_atoms * n_guess_neigh, eel%fmm_near_field_list%ci)
        eel%fmm_near_field_list%ri(1) = 1

        do i=1, eel%top%mm_atoms
            eel%fmm_near_field_list%ri(i+1) = eel%fmm_near_field_list%ri(i)
            inode = eel%tree%particle_to_node(i)

            do node_idx=eel%tree%near_nl%ri(inode), eel%tree%near_nl%ri(inode+1)-1
                jnode = eel%tree%near_nl%ci(node_idx)

                do part_idx=eel%tree%particle_list%ri(jnode), eel%tree%particle_list%ri(jnode+1)-1
                    j = eel%tree%particle_list%ci(part_idx)
                    if(i == j) cycle
                    eel%fmm_near_field_list%ci(eel%fmm_near_field_list%ri(i+1)) = j
                    eel%fmm_near_field_list%ri(i+1) = eel%fmm_near_field_list%ri(i+1) + 1
                    if(eel%fmm_near_field_list%ri(i+1) > size(eel%fmm_near_field_list%ci)) then
                        ! Next iteration would cause an overflow
                        call reallocate_mat(eel%fmm_near_field_list, size(eel%fmm_near_field_list%ci) &
                                                                     + eel%top%mm_atoms * n_guess_neigh)
                    end if
                end do
            end do
        end do
        ! In the end shrink the row list
        call reallocate_mat(eel%fmm_near_field_list, eel%fmm_near_field_list%ri(eel%top%mm_atoms+1)-1)
    end subroutine

    subroutine fmm_solve_for_multipoles(fmm_obj, q, use_q, mu, use_mu, quad, use_quad)
        use mod_constants, only: pi
        
        implicit none

        type(fmm_type), intent(inout) :: fmm_obj
        !! FMM object, it should be already initialized
        real(rp), intent(in) :: q(:)
        logical(lp), intent(in) :: use_q
        
        real(rp), intent(in) :: mu(:, :)
        logical(lp), intent(in) :: use_mu
        
        real(rp), intent(in) :: quad(:, :)
        logical(lp), intent(in) :: use_quad

        real(rp), allocatable :: multipoles_sphe(:, :)
        integer(ip) :: n_comp, i

        call time_push
        if(use_quad) then
            n_comp = 9
        else if(use_mu) then
            n_comp = 4
        else if(use_q) then 
            n_comp = 1 
        else
            n_comp = 0
        end if

        if(n_comp == 0) then
            call fatal_error("fmm_solve_for_multipoles called without any input source.")
        end if
        
        allocate(multipoles_sphe(n_comp, fmm_obj%tree%n_particles))
        
        multipoles_sphe = 0.0
        
        if(use_q) then
            if(size(q,1) /= fmm_obj%tree%n_particles) then
                call fatal_error("charges array has wrong size in fmm_solve_for_multipoles")
            end if
            multipoles_sphe(1,:) = q / sqrt(4.0 * pi)
        end if
        
        if(use_mu) then
            if(size(mu,1) /= 3 .or. size(mu,2) /= fmm_obj%tree%n_particles) then
                call fatal_error("dipoles array has wrong size in fmm_solve_for_multipoles")
            end if
            multipoles_sphe(2,:) = mu(2,:) / sqrt(4.0 * pi / 3.0)
            multipoles_sphe(3,:) = mu(3,:) / sqrt(4.0 * pi / 3.0)
            multipoles_sphe(4,:) = mu(1,:) / sqrt(4.0 * pi / 3.0)
        end if
        
        if(use_quad) then
            if(size(quad,1) /= 6 .or. size(quad,2) /= fmm_obj%tree%n_particles) then
                call fatal_error("dipoles array has wrong size in fmm_solve_for_multipoles")
            end if
            multipoles_sphe(5,:) = 2.0 * quad(_xy_,:) * sqrt(15.0 / (4.0 * pi))
            multipoles_sphe(6,:) = 2.0 *quad(_yz_,:) * sqrt(15.0 / (4.0 * pi))
            multipoles_sphe(7,:) = 6.0 / 2.0 * quad(_zz_,:) * sqrt(5.0/(4.0*pi))
            multipoles_sphe(8,:) = 2.0 * quad(_xz_,:) * sqrt(15.0 / (4.0 * pi))
            multipoles_sphe(9,:) = (quad(_xx_,:) - quad(_yy_,:)) * sqrt(15.0/(4.0*pi))
        end if
        call time_pull('cart2sphe')
        
        ! Load FMM
        call time_push
        call tree_p2m(fmm_obj, multipoles_sphe, 2)
        call time_pull('tree P2M')
        call time_push()
        call time_push
        call tree_m2m(fmm_obj)
        call time_pull("tree M2M")
        call time_push
        call tree_m2l(fmm_obj)
        call time_pull("tree M2L")
        call time_push
        call tree_l2l(fmm_obj)
        call time_pull("tree L2L")
        call time_pull('fmm_solve')

        deallocate(multipoles_sphe)
        
    end subroutine

end module mod_electrostatics
