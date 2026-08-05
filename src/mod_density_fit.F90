#include "f_cart_components.h"
module mod_density_fit
!! This module implements density fitting of QM charge distributions
!! onto a set of fitting points (typically MM atom positions).

    use mod_memory, only: ip, rp, lp, mallocate, mfree
    use mod_constants, only: ommp_df_solver_svd, &
                             ommp_df_svd_rcond_default, &
                             ommp_df_atoms, &
                             ommp_df_fibonacci, &
                             ommp_df_cubic, &
                             ommp_df_mm_top, ommp_df_qm_top, &
                             OMMP_VERBOSE_DEBUG
    use mod_io, only: fatal_error, ommp_message
    use mod_topology, only: ommp_topology_type
    use mod_adjacency_mat, only: yale_sparse, allocate_yale_sparse, free_yale_sparse
    use mod_profiling, only: time_pull, time_push

    implicit none
    private

    type ommp_density_fit_type
        integer(ip) :: n_pts = 0
        !! Number of fitting points

        integer(ip) :: n_charges = 0
        !! Number of target charges

        integer(ip) :: df_method = ommp_df_solver_svd ! NOT USED
        !! Method used for solving the linear system (PINV, SVD, or NORMAL)

        real(rp) :: svd_rcond = ommp_df_svd_rcond_default ! NOT USED
        !! SVD truncation threshold (singular values < rcond * max are zeroed)

        !! Point generation strategies
        integer(ip) :: charge_point_type = 0
        !! Grid type for charge points: ommp_df_atoms, ommp_df_fibonacci, ommp_df_cubic

        integer(ip) :: charge_n_pts_per_atom = 0
        !! Number of charge points per source atom (for fibonacci grids)

        real(rp) :: charge_radius = 0.0_rp
        !! Radius parameter for charge point generation

        integer(ip) :: fit_point_type = 0
        !! Grid type for fit points: ommp_df_atoms, ommp_df_fibonacci, ommp_df_cubic

        integer(ip) :: fit_n_pts_per_atom = 0
        !! Number of fit points per source atom (for grid-based strategies)

        real(rp) :: fit_radius = 0.0_rp
        !! Radius parameter for fit point generation

        !! Source topologies for grid generation (can be qm or mm)
        type(ommp_topology_type), pointer :: charge_top
        !! Pointer to topology providing charge point coordinates
        integer(ip) :: charge_top_type = 0
        !! Type of charge point topology: ommp_df_mm_top (1) or ommp_df_qm_top (2)

        type(ommp_topology_type), pointer :: fit_top
        !! Pointer to topology providing fit point coordinates
        integer(ip) :: fit_top_type = 0
        !! Type of fit point topology: ommp_df_mm_top (1) or ommp_df_qm_top (2)

        !! Full QM and MM topology pointers (always set at init, used by gradient routines)
        type(ommp_topology_type), pointer :: qm_top
        !! Pointer to the QM topology (used in gradient computations)
        type(ommp_topology_type), pointer :: mm_top
        !! Pointer to the MM topology (used in gradient computations)

        real(rp), allocatable :: charge_coord(:,:)
        !! Coordinates of the charge positions (3 x n_charges)

        real(rp), allocatable :: fit_point_coord(:,:)
        !! Coordinates of the fitting points (3 x n_pts)

        real(rp), allocatable :: target_charges(:)
        !! Target charges to be fitted

        real(rp), allocatable :: fit_potential(:)
        !! Potential values at fitting points

        real(rp), allocatable :: X(:,:)
        !! Design matrix, dimensions (n_pts x n_charges)

        real(rp), allocatable :: Xinv(:,:)
        !! Pseudoinverse design matrix, dimensions (n_charges x n_pts).
        !! Satisfies fitted_charges = Xinv @ fit_potential

        logical(lp) :: initialized = .false.
        !! Flag indicating whether the object is initialized

        logical(lp) :: xinv_done = .false.
        !! Flag indicating whether Xinv has been computed

        logical(lp) :: x_done = .false.
        !! Flag indicating whether X has been computed

        logical(lp) :: fit_done = .false.
        !! Flag indicating whether charges have been computed from fit_potential

        real(rp), allocatable :: V_m2q(:)
        !! Electrostatic potential from MM static multipoles at charge coordinates

        real(rp), allocatable :: V_p2q(:,:)
        !! Electrostatic potential from MM induced dipoles at charge coordinates

        logical(lp) :: V_m2q_done = .false.
        !! Flag indicating whether V_mm2df has been computed

        logical(lp) :: V_p2q_done = .false.
        !! Flag indicating whether V_pd2df has been computed

        real(rp), allocatable :: VXI_m(:)
        !! Projected static quantity: V_m2q @ Xinv, for Fock matrix

        real(rp), allocatable :: VXI_p(:,:)
        !! Projected dipole quantity: V_p2q @ Xinv, for Fock matrix

        logical(lp) :: VXI_m_done = .false.
        !! Flag indicating whether VXI_m has been computed

        logical(lp) :: VXI_p_done = .false.
        !! Flag indicating whether VXI_p has been computed

        real(rp), allocatable :: E_q2p(:,:)
        !! Electric field from fitted charges at polarizable sites

        logical(lp) :: E_q2p_done = .false.
        !! Flag indicating whether E_q2p has been computed

        real(rp), allocatable :: E_m2q(:,:,:)
        !! Electric field generated by MM static multipoles at charge sites.
        !! Shape: (3, n_mm, n_charges) - E_m2q(:,i,j) is the field from
        !! MM atom i at charge site j.

        logical(lp) :: E_m2q_done = .false.
        !! Flag indicating whether E_m2q has been computed

        real(rp), allocatable :: E_p2q(:,:,:,:)
        !! Electric field generated by IPDs at charge sites.
        !! Shape: (3, n_mm, n_charges, n_ipd)
        !!   - only populated for polarizable atoms (where mm_polar(i) > 0)
        !!   - E_p2q(:,i,j,k) is the field from IPD component k of MM atom i
        !!     evaluated at charge site j

        logical(lp) :: E_p2q_done = .false.
        !! Flag indicating whether E_p2q has been computed

        real(rp), allocatable :: E_q2m(:,:)
        !! Electric field from fitted charges at MM atom positions.
        !! Shape: (3, n_mm_atoms) - E_q2M(:,k) is the electric field at MM atom k
        !! generated by all fitted charges.

        logical(lp) :: E_q2m_done = .false.
        !! Flag indicating whether E_q2M has been computed

        real(rp), allocatable :: GEF_q2m(:,:)
        !! Gradient of the electric field from fitted charges at MM atom positions.
        !! Shape: (6, n_mm_atoms) - 6-component field gradient tensor (symmetric)
        !! at each MM atom: (xx, xy, yy, xz, yz, zz).

        logical(lp) :: GEF_q2m_done = .false.
        !! Flag indicating whether GEF_q2M has been computed

        real(rp), allocatable :: dX_dr(:,:,:)
        !! Derivative of the design matrix X w.r.t. fit-point coordinates.
        !! Shape: (3,nfit,ncharges)
        !! It is compressed with a lacking dimension because
        !! dX_dr_i = \delta_i,k kernel(2) (r_k, r_q)
        !! so there is no need to store it bigger !

        logical(lp) :: dX_dr_done = .false.
        !! Flag indicating whether dX_dr has been computed

        real(rp) :: E_pol_ene
        !! Polarization energy from fitted-charge electric field:
        !! E = -0.5 * sum_j ipd(:,j) .dot. E_q2p(:,j)

        logical(lp) :: E_pol_ene_done = .false.
        !! Flag indicating whether E_pol_ene has been computed

        !! Gradient matrices all shaped 3 ngrid (or ncharges) x 3 nmm (or nqm)
        logical(lp) :: nabla_done = .false.
        !! Label for indicating that nabla have been computed.
        real(rp), allocatable :: nabla_g_mm(:,:)
        !! Gradient matrix: grid wrt MM coordinates
        logical(lp) :: nabla_g_mm_is_null = .false.
        !! .true. if nabla_g_mm is the null matrix
        logical(lp) :: nabla_g_mm_is_identity = .false.
        !! .true. if nabla_g_mm is the identity matrix
        logical(lp) :: nabla_g_mm_is_sparse = .false.
        !! .true. if nabla_g_mm is stored in Yale sparse format
        type(yale_sparse) :: nabla_g_mm_sparse

        real(rp), allocatable :: nabla_g_qm(:,:)
        !! Gradient matrix: grid wrt QM coordinates
        logical(lp) :: nabla_g_qm_is_null = .false.
        !! .true. if nabla_g_qm is the null matrix
        logical(lp) :: nabla_g_qm_is_identity = .false.
        !! .true. if nabla_g_qm is the identity matrix
        logical(lp) :: nabla_g_qm_is_sparse = .false.
        !! .true. if nabla_g_qm is stored in Yale sparse format
        type(yale_sparse) :: nabla_g_qm_sparse

        real(rp), allocatable :: nabla_q_qm(:,:)
        !! Gradient matrix: q wrt QM coordinates
        logical(lp) :: nabla_q_qm_is_null = .false.
        !! .true. if nabla_q_qm is the null matrix
        logical(lp) :: nabla_q_qm_is_identity = .false.
        !! .true. if nabla_q_qm is the identity matrix
        logical(lp) :: nabla_q_qm_is_sparse = .false.
        !! .true. if nabla_q_qm is stored in Yale sparse format
        type(yale_sparse) :: nabla_q_qm_sparse

        real(rp), allocatable :: nabla_q_mm(:,:)
        !! Gradient matrix: q wrt MM coordinates
        logical(lp) :: nabla_q_mm_is_null = .false.
        !! .true. if nabla_q_mm is the null matrix
        logical(lp) :: nabla_q_mm_is_identity = .false.
        !! .true. if nabla_q_mm is the identity matrix
        logical(lp) :: nabla_q_mm_is_sparse = .false.
        !! .true. if nabla_q_mm is stored in Yale sparse format
        type(yale_sparse) :: nabla_q_mm_sparse

        !! Snapshot of the last MM and QM coordinates.
        !! Used by df_update to detect coordinate changes and
        !! regenerate grids only when necessary.
        real(rp), allocatable :: last_mm_coords(:,:)
        !! Snapshot of mm_top%cmm

        real(rp), allocatable :: last_qm_coords(:,:)
        !! Snapshot of qm_top%cmm

    end type ommp_density_fit_type

    public :: ommp_density_fit_type
    public :: df_init, df_terminate
    public :: df_solve
    public :: df_electrostatic_static, df_electrostatic_dipoles
    public :: df_project_static, df_project_dipoles
    public :: df_e_field_to_pol, df_e_field_pol_ene
    public :: df_geomgrad
    public :: df_update
    !! Accessors for gradient-related quantities (exposed outside)

contains

    subroutine generate_grid_from_topo(top, out_coord, n_tot, n_pts_per_atom, radius, grid_type, label)
        !! Generic grid generation from a topology.
        !!
        !! Supported grid_type values (from mod_constants):
        !!   1 (ommp_df_atoms) - atom positions
        !!   2 (ommp_df_fibonacci) - fibonacci sphere per atom
        !!   3 (ommp_df_cubic) - 7-point cubic grid (1 center + 6 face centers)

        implicit none

        type(ommp_topology_type), intent(in), pointer :: top
        !! Topology providing the base atom coordinates

        real(rp), allocatable, intent(inout) :: out_coord(:,:)
        !! Output coordinate array (3 x n_tot), allocated by this subroutine

        integer(ip), intent(out) :: n_tot
        !! Total number of points generated

        integer(ip), intent(in) :: n_pts_per_atom
        !! Points per atom (for fibonacci only)

        real(rp), intent(in) :: radius
        !! Radius parameter (for fibonacci and cubic)

        integer(ip), intent(in) :: grid_type
        !! Grid type selector

        character(len=*), intent(in) :: label
        !! Label for this grid (used in messages)

        integer(ip) :: n_atoms, ii, i, idx
        real(rp) :: x0, y0, z0, theta, phi, y
        character(len=256) :: msg

        real(rp), parameter :: golden_angle = 2.0_rp * atan(1.0_rp) * (3.0_rp - sqrt(5.0_rp))
        !! Golden angle in radians (~2.39996 rad) for fibonacci sphere distribution

        call time_push

        if(.not. associated(top)) then
            call fatal_error('generate_grid_from_topo: topology is not associated.')
        end if
        n_atoms = top%mm_atoms
        if(n_atoms == 0) then
            call fatal_error('generate_grid_from_topo: topology has no atoms.')
        end if

        call mfree('[generate_grid_from_topo] out_coord', out_coord)

        select case(grid_type)
        case(ommp_df_atoms)
            !! Atom positions
            call mallocate('[generate_grid_from_topo] out_coord', 3_ip, n_atoms, out_coord)
            out_coord(1:3, 1:n_atoms) = top%cmm(1:3, 1:n_atoms)
            n_tot = n_atoms
            write(msg, '(A,I0,A,A,A)') 'Generated ', n_atoms, ' ', label, ' points from atom positions'
            call ommp_message(trim(msg), 1, 'df')

        case(ommp_df_fibonacci)
            !! Fibonacci sphere per atom
            if(n_pts_per_atom <= 0) then
                call fatal_error('generate_grid_from_topo: fibonacci grid requires n_pts_per_atom > 0.')
            end if
            if(radius <= 0.0_rp) then
                call fatal_error('generate_grid_from_topo: fibonacci grid requires radius > 0.')
            end if
            call mallocate('[generate_grid_from_topo] out_coord', 3_ip, &
                           n_atoms * n_pts_per_atom, out_coord)
            do ii = 1, n_atoms
                x0 = top%cmm(1, ii)
                y0 = top%cmm(2, ii)
                z0 = top%cmm(3, ii)
                do i = 1, n_pts_per_atom
                    y = (2.0_rp * real(i, rp) - 1.0_rp) / real(n_pts_per_atom, rp) - 1.0_rp
                    theta = real(golden_angle, rp) * real(i, rp)
                    phi = sqrt(1.0_rp - y**2)
                    idx = (ii - 1) * n_pts_per_atom + i
                    out_coord(1, idx) = x0 + cos(theta) * phi * radius
                    out_coord(2, idx) = y0 + y * radius
                    out_coord(3, idx) = z0 + sin(theta) * phi * radius
                end do
            end do
            n_tot = n_atoms * n_pts_per_atom
            write(msg, '(A,I0,A,A,A,I0,A)') 'Generated ', n_pts_per_atom, ' ', &
                   label, ' fibonacci points per atom (', n_atoms, ' atoms)'
            call ommp_message(trim(msg), 1, 'df')

        case(ommp_df_cubic)
            !! Cubic grid: 7 points per atom (1 center + 6 face centers)
            if(radius <= 0.0_rp) then
                call fatal_error('generate_grid_from_topo: cubic grid requires radius > 0.')
            end if
            call mallocate('[generate_grid_from_topo] out_coord', 3_ip, n_atoms * 7_ip, out_coord)
            do ii = 1, n_atoms
                i = (ii - 1) * 7 + 1
                x0 = top%cmm(1, ii)
                y0 = top%cmm(2, ii)
                z0 = top%cmm(3, ii)
                out_coord(1, i:i+6) = x0
                out_coord(2, i:i+6) = y0
                out_coord(3, i:i+6) = z0
                out_coord(1, i+1) = x0 + radius
                out_coord(1, i+2) = x0 - radius
                out_coord(2, i+3) = y0 + radius
                out_coord(2, i+4) = y0 - radius
                out_coord(3, i+5) = z0 + radius
                out_coord(3, i+6) = z0 - radius
            end do
            n_tot = n_atoms * 7_ip
            write(msg, '(A,A,A,I0,A)') 'Generated 7 ', label, ' cubic points per atom (', &
                   n_atoms, ' atoms)'
            call ommp_message(trim(msg), 1, 'df')

        case default
            call fatal_error('generate_grid_from_topo: unknown grid_type.')
        end select
        
        call time_pull("Grid generation")

    end subroutine generate_grid_from_topo

    subroutine compute_nabla_matrices(df)

        implicit none

        type(ommp_density_fit_type), intent(inout) :: df
        integer(ip) :: i, j

        call df_update(df)

        if(df%nabla_done) return

        call time_push

        select case(df%charge_top_type)
        case(ommp_df_qm_top)
            ! Charge positions depend only on QM position

            df%nabla_q_mm_is_null = .true.
            df%nabla_q_mm_is_identity = .false.
            df%nabla_q_mm_is_sparse = .false.

            select case(df%charge_point_type)
            case(ommp_df_atoms)
                df%nabla_q_qm_is_null = .false.
                df%nabla_q_qm_is_identity = .true.

            case(ommp_df_fibonacci, ommp_df_cubic)

                df%nabla_q_qm_is_null = .false.
                df%nabla_q_qm_is_identity = .false.
                df%nabla_q_qm_is_sparse = .true.

                ! In this cases each point belong to a specific QM atom
                call allocate_yale_sparse(df%nabla_q_qm_sparse, df%qm_top%mm_atoms, df%n_charges)
                do i=1, df%qm_top%mm_atoms
                    df%nabla_q_qm_sparse%ri(i) = (i-1) * df%charge_n_pts_per_atom + 1
                    df%nabla_q_qm_sparse%ri(i+1) = i * df%charge_n_pts_per_atom + 1
                    do j=df%nabla_q_qm_sparse%ri(i), df%nabla_q_qm_sparse%ri(i+1)-1
                        df%nabla_q_qm_sparse%ci(j) = j
                    end do
                end do

            case default
                call fatal_error('compute_nabla_matrices: unknown charge_point_type.')
            end select
        case(ommp_df_mm_top)
            call fatal_error("Unsupported case: charge points generated from MM topology only.")
        case default
            call fatal_error('compute_nabla_matrices: unknown charge_top_type.')
        end select

        select case(df%fit_top_type)
        case(ommp_df_mm_top)
            ! Obvius case fit grid is built on mm atoms

            df%nabla_g_qm_is_null = .true.
            df%nabla_g_qm_is_identity = .false.
            df%nabla_g_qm_is_sparse = .false.

            select case(df%fit_point_type)
            case(ommp_df_atoms)
                df%nabla_g_mm_is_null = .false.
                df%nabla_g_mm_is_identity = .true.
                df%nabla_g_mm_is_sparse = .false.

            case(ommp_df_fibonacci, ommp_df_cubic)

                df%nabla_g_mm_is_null = .false.
                df%nabla_g_mm_is_identity = .false.
                df%nabla_g_mm_is_sparse = .true.

                ! In this cases each point belong to a specific MM atom
                call allocate_yale_sparse(df%nabla_g_mm_sparse, df%mm_top%mm_atoms, df%n_pts)
                do i=1, df%mm_top%mm_atoms
                    df%nabla_g_mm_sparse%ri(i) = (i-1) * df%fit_n_pts_per_atom + 1
                    df%nabla_g_mm_sparse%ri(i+1) = i * df%fit_n_pts_per_atom + 1
                    do j=df%nabla_g_mm_sparse%ri(i), df%nabla_g_mm_sparse%ri(i+1)-1
                        df%nabla_g_mm_sparse%ci(j) = j
                    end do
                end do

            case default
                call fatal_error('compute_nabla_matrices: unknown fit_point_type.')
            end select
        case(ommp_df_qm_top)
            call fatal_error("Unsupported case: grid points generated from QM topology only.")
        case default
            call fatal_error('compute_nabla_matrices: unknown fit_top_type.')
        end select

        df%nabla_done = .true.

        call time_pull("DF - Nabla Matrices")

    end subroutine compute_nabla_matrices

    subroutine df_update(df)
        !! Update density fitting grids when coordinates change.
        !!
        !! Detects whether coordinates have changed by comparing against
        !! snapshot arrays. On first call (snapshots not yet allocated),
        !! always regenerates. Grids are deterministic functions of the
        !! topology, so a simple coordinate comparison suffices.
        !!
        !! This is the single entry point for coordinate updates.
        !! Called from update_coordinates (MM changes) and from
        !! wherever QM coordinates are updated. Also called once from
        !! df_init for the initial grid generation.

        use mod_constants, only: eps_rp

        implicit none

        type(ommp_density_fit_type), intent(inout) :: df

        logical(lp) :: need_update

        !! First call: snapshots not allocated → always update
        if(.not. allocated(df%last_mm_coords) .or. .not. allocated(df%last_qm_coords)) then
            need_update = .true.
            if(.not. allocated(df%last_mm_coords)) &
                call mallocate('[df_update] last_mm_coords', 3_ip, df%mm_top%mm_atoms, df%last_mm_coords)
            if(.not. allocated(df%last_qm_coords)) &
                call mallocate('[df_update] last_qm_coords', 3_ip, df%qm_top%mm_atoms, df%last_qm_coords)
        else
            !! Compare against stored snapshots
            need_update = any(abs(df%mm_top%cmm - df%last_mm_coords) > eps_rp) .or. &
                          any(abs(df%qm_top%cmm - df%last_qm_coords) > eps_rp)
        end if

        if(.not. need_update) return

        !! Regenerate both grids from current topology
        call generate_grid_from_topo(df%charge_top, df%charge_coord, df%n_charges, &
                                     df%charge_n_pts_per_atom, df%charge_radius, &
                                     df%charge_point_type, 'charge')
        call generate_grid_from_topo(df%fit_top, df%fit_point_coord, df%n_pts, &
                                     df%fit_n_pts_per_atom, df%fit_radius, &
                                     df%fit_point_type, 'fit')

        !! Invalidate all derived quantities
        df%xinv_done = .false.
        df%x_done = .false.
        df%fit_done = .false.
        df%V_m2q_done = .false.
        df%V_p2q_done = .false.
        df%VXI_m_done = .false.
        df%VXI_p_done = .false.
        df%E_q2p_done = .false.
        df%E_q2M_done = .false.
        df%GEF_q2M_done = .false.
        df%E_pol_ene_done = .false.
        df%E_m2q_done = .false.
        df%E_p2q_done = .false.
        df%dX_dr_done = .false.

        df%last_mm_coords(1:3, 1:df%mm_top%mm_atoms) = df%mm_top%cmm(1:3, 1:df%mm_top%mm_atoms)
        df%last_qm_coords(1:3, 1:df%qm_top%mm_atoms) = df%qm_top%cmm(1:3, 1:df%qm_top%mm_atoms)

    end subroutine df_update

    subroutine df_init(df, mm_top, qm_top, &
                       charge_point_type, charge_n_pts_per_atom, charge_radius, &
                       fit_point_type, fit_n_pts_per_atom, fit_radius, &
                       charge_top_type, fit_top_type)
        !! Initialize the density fit object.
        !! Coordinates are generated from the provided topology pointers
        !! according to the point generation strategy.

        use mod_memory, only: mallocate
        use mod_constants, only: ommp_df_mm_top, ommp_df_qm_top

        implicit none

        type(ommp_density_fit_type), intent(inout) :: df
        type(ommp_topology_type), intent(in), target :: mm_top
        !! MM topology (always provided, used when charge/fit type = ommp_df_mm_top)
        type(ommp_topology_type), intent(in), target :: qm_top
        !! QM topology (always provided, used when charge/fit type = ommp_df_qm_top)
        integer(ip), intent(in) :: charge_point_type
        !! Type of charge point source
        integer(ip), intent(in) :: charge_n_pts_per_atom
        !! Number of charge points per source atom
        real(rp), intent(in) :: charge_radius
        !! Radius parameter for charge point generation
        integer(ip), intent(in) :: fit_point_type
        !! Type of fit point source
        integer(ip), intent(in) :: fit_n_pts_per_atom
        !! Number of fit points per source atom
        real(rp), intent(in) :: fit_radius
        !! Radius parameter for fit point generation
        integer(ip), intent(in) :: charge_top_type
        !! Type of charge topology: ommp_df_mm_top (1) or ommp_df_qm_top (2)
        integer(ip), intent(in) :: fit_top_type
        !! Type of fit topology: ommp_df_mm_top (1) or ommp_df_qm_top (2)

        integer(ip) :: i, j, idx
        real(rp) :: x0, y0, z0, theta, phi, y

        if(df%initialized) then
            call fatal_error("Density fit object already initialized!")
        end if

        df%mm_top => mm_top
        df%qm_top => qm_top

        select case(charge_top_type)
        case(ommp_df_mm_top)
            df%charge_top => mm_top
        case(ommp_df_qm_top)
            df%charge_top => qm_top
        case default
            call fatal_error('df_init: unknown charge_top_type.')
        end select

        select case(fit_top_type)
        case(ommp_df_mm_top)
            df%fit_top => mm_top
        case(ommp_df_qm_top)
            df%fit_top => qm_top
        case default
            call fatal_error('df_init: unknown fit_top_type.')
        end select

        if(df%charge_top%mm_atoms == 0) then
            call fatal_error('df_init: charge topology has no atoms.')
        end if
        if(df%fit_top%mm_atoms == 0) then
            call fatal_error('df_init: fit topology has no atoms.')
        end if

        !! Store point generation configuration
        df%charge_point_type = charge_point_type
        df%charge_n_pts_per_atom = charge_n_pts_per_atom
        df%charge_radius = charge_radius
        df%fit_point_type = fit_point_type
        df%fit_n_pts_per_atom = fit_n_pts_per_atom
        df%fit_radius = fit_radius
        df%charge_top_type = charge_top_type
        df%fit_top_type = fit_top_type

        call df_update(df)

        call mallocate('df_init [target_charges]', df%n_charges, df%target_charges)
        df%target_charges = 0.0_rp

        call mallocate('df_init [fit_potential]', df%n_pts, df%fit_potential)
        df%fit_potential = 0.0_rp

        call mallocate('df_init [X]', df%n_pts, df%n_charges, df%X)
        call mallocate('df_init [Xinv]', df%n_charges, df%n_pts, df%Xinv)

        df%E_pol_ene = 0.0_rp

        df%initialized = .true.
    end subroutine df_init

    subroutine df_terminate(df)
        !! Terminate and free memory in the density fit object

        use mod_memory, only: mfree

        implicit none

        type(ommp_density_fit_type), intent(inout) :: df

        if(.not. df%initialized) return

        call mfree('df_terminate [charge_coord]', df%charge_coord)
        call mfree('df_terminate [fit_point_coord]', df%fit_point_coord)
        call mfree('df_terminate [target_charges]', df%target_charges)
        call mfree('df_terminate [fit_potential]', df%fit_potential)
        call mfree('df_terminate [X]', df%X)
        call mfree('df_terminate [Xinv]', df%Xinv)
        call mfree('df_terminate [V_mm2df]', df%V_m2q)
        call mfree('df_terminate [V_pd2df]', df%V_p2q)
        call mfree('df_terminate [VXI_m]', df%VXI_m)
        call mfree('df_terminate [VXI_p]', df%VXI_p)
        call mfree('df_terminate [E_q2p]', df%E_q2p)
        call mfree('df_terminate [E_MM2Q]', df%E_m2q)
        call mfree('df_terminate [E_p2q]', df%E_p2q)
        call mfree('df_terminate [E_q2M]', df%E_q2M)
        call mfree('df_terminate [GEF_q2M]', df%GEF_q2M)
        call mfree('df_terminate [dX_dr]', df%dX_dr)
        call mfree('df_terminate [last_mm_coords]', df%last_mm_coords)
        call mfree('df_terminate [last_qm_coords]', df%last_qm_coords)
        call mfree('df_terminate [nabla_g_mm]', df%nabla_g_mm)
        call mfree('df_terminate [nabla_g_qm]', df%nabla_g_qm)
        call mfree('df_terminate [nabla_q_qm]', df%nabla_q_qm)
        call mfree('df_terminate [nabla_q_mm]', df%nabla_q_mm)

        if(df%nabla_g_mm_is_sparse) call free_yale_sparse(df%nabla_g_mm_sparse)
        if(df%nabla_g_qm_is_sparse) call free_yale_sparse(df%nabla_g_qm_sparse)
        if(df%nabla_q_qm_is_sparse) call free_yale_sparse(df%nabla_q_qm_sparse)
        if(df%nabla_q_mm_is_sparse) call free_yale_sparse(df%nabla_q_mm_sparse)

    end subroutine df_terminate

    subroutine df_build_X(df)
        !! Build the design matrix X where X(i,j) = 1 / |r_charge_i - r_fit_j|.
        !! The matrix has dimensions (n_charges x n_pts) in Fortran layout.

        use mod_constants, only: eps_rp

        implicit none

        type(ommp_density_fit_type), intent(inout) :: df
        integer(ip) :: i, j
        real(rp) :: dr(3), dist

        call df_update(df)

        call time_push

        if(.not. df%initialized) then
            call fatal_error("Density fit object not initialized!")
        end if

        !$omp parallel do default(shared) schedule(static) &
        !$omp private(i,j,dr,dist)
        do i = 1, df%n_charges
            do j = 1, df%n_pts
                dr = df%charge_coord(:,i) - df%fit_point_coord(:,j)
                dist = norm2(dr)
                if(dist < eps_rp) then
                    !$omp critical
                    call fatal_error('Charge and fitting point coincide!')
                    !$omp end critical
                end if
                df%X(j,i) = 1.0_rp / dist
            end do
        end do

        call time_pull("DF - Computing X")
    end subroutine df_build_X

    subroutine df_solve(df)
        !! Solve the density fitting problem by matrix multiplication:
        !!   fitted_charges = Xinv @ fit_potential
        !!
        !! Lazy evaluation:
        !!   - X is built on first request (if not already done)
        !!   - Xinv is computed from X on first request (if not already done)
        !!   - Matrix-vector multiply is performed every time (fit_done is reset
        !!     when fit_potential changes, so this runs after every update)
        !!
        !! Residual metric (local, not stored):
        !!   residual_norm = || fit_potential - X @ fitted_charges || / || fit_potential ||
        !!   This is computed from local temporaries and discarded after df_solve returns.

        use mod_memory, only: mallocate, mfree

        implicit none

        type(ommp_density_fit_type), intent(inout) :: df
        real(rp), allocatable :: predicted_pot(:)
        real(rp) :: fit_norm, residual_norm
        integer(ip) :: i
        character(len=256) :: msg

        if(.not. df%initialized) then
            call fatal_error("Density fit object not initialized!")
        end if

        call df_update(df)

        !! Lazy: compute Xinv if not already done
        if(.not. df%xinv_done) then
            call df_compute_Xinv_svd(df)
        end if

        df%fit_done = .false.
        df%V_p2q_done = .false.
        df%VXI_p_done = .false.
        df%E_q2p_done = .false.
        df%E_q2M_done = .false.
        df%GEF_q2M_done = .false.
        df%E_pol_ene_done = .false.

        call time_push
        !! Solve: fitted_charges = Xinv @ fit_potential
        call dgemv('N', df%n_charges, df%n_pts, 1.0_rp, df%Xinv, df%n_charges, df%fit_potential, 1, 0.0_rp, df%target_charges, 1)
        call time_pull("DF - Computing fit charges")

        !! Forward: predicted_potential = X @ fitted_charges
        call time_push
        call mallocate('df_solve [predicted_pot]', df%n_pts, predicted_pot)
        call dgemv('N', df%n_pts, df%n_charges, 1.0_rp, df%X, df%n_pts, df%target_charges, 1, 0.0_rp, predicted_pot, 1)

        !! Residual L2 norm: || fit_potential - predicted_potential || / || fit_potential ||
        fit_norm = 0.0_rp
        residual_norm = 0.0_rp
        do i = 1, df%n_pts
            residual_norm = residual_norm + (df%fit_potential(i) - predicted_pot(i))**2
            fit_norm = fit_norm + df%fit_potential(i)**2
        end do
        residual_norm = sqrt(residual_norm / fit_norm)

        call mfree('df_solve [predicted_pot]', predicted_pot)
        
        call time_pull("DF - Computing residues")
        df%fit_done = .true.

        !! Report residual norm as percentage
        write(msg, '(A,F10.6,A)') 'Density fit residual norm = ', residual_norm * 100.0_rp, ' %'
        call ommp_message(trim(msg), 1, 'df')

        !! Report total fitted charge
        write(msg, '(A,F12.2)') 'Total fitted charge = ', sum(df%target_charges)
        call ommp_message(trim(msg), 1, 'df')
    end subroutine df_solve

subroutine df_compute_Xinv_svd(df)

    implicit none

    type(ommp_density_fit_type), intent(inout) :: df
    integer(ip) :: min_dim, lwork, info, i
    integer(ip), dimension(:), allocatable :: iwork
    real(rp), dimension(:), allocatable :: s, work
    real(rp), dimension(:,:), allocatable :: u, vt, tmp_X

    if(.not. df%initialized) then
        call fatal_error("Density fit object not initialized!")
    end if

    call df_update(df)

    if(.not. df%x_done) then
        call df_build_X(df)
    end if

    call time_push
    if(df%n_charges == 0 .or. df%n_pts == 0) &
        call fatal_error("Either target or fit grids in density-fitting have no points.")

    min_dim = min(df%n_charges, df%n_pts)

    !! Allocate SVD workspace — U is now n_pts x min_dim (economy SVD)
    allocate(tmp_X(df%n_pts, df%n_charges))
    allocate(s(min_dim))
    allocate(u(df%n_pts, min_dim))      ! was (n_pts, n_pts)
    allocate(vt(min_dim, df%n_charges)) ! was (n_charges, n_charges)
    allocate(iwork(8 * min_dim))

    call time_push
    tmp_X = df%X

    !! Query optimal workspace size using economy SVD ('S')
    lwork = -1
    allocate(work(1))
    call dgesdd('S', df%n_pts, df%n_charges, tmp_X, df%n_pts, s, &
                u, df%n_pts, vt, min_dim, &   ! ldvt is now min_dim, not n_charges
                work, lwork, iwork, info)
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))

    !! Compute economy SVD
    call dgesdd('S', df%n_pts, df%n_charges, tmp_X, df%n_pts, s, &
                u, df%n_pts, vt, min_dim, &
                work, lwork, iwork, info)
    if(info /= 0) call fatal_error('dgesdd SVD failed')
    call time_pull("DF - SVD dgesdd")

    call time_push
    !! Scale columns of U by 1/s(i): u(:,i) <- u(:,i) / s(i)
    do i = 1, min_dim
        u(:,i) = u(:,i) / s(i)
    end do

    !! Pseudoinverse: Xinv = Vt^T @ (U @ S^{-1})^T
    !! = Vt^T (min_dim x n_charges)^T  @  U_scaled^T (n_pts x min_dim)^T
    !! dgemm: Xinv(n_charges, n_pts) = Vt^T(n_charges, min_dim) @ U_scaled^T(min_dim, n_pts)
    call dgemm('T', 'T', df%n_charges, df%n_pts, min_dim, &
               1.0_rp, vt, min_dim, &       ! lda is now min_dim, not n_charges
               u, df%n_pts, &
               0.0_rp, df%Xinv, df%n_charges)

    call time_pull("Pseudoinverse")

    deallocate(s, u, vt, work, iwork, tmp_X)
    df%xinv_done = .true.
    call ommp_message('SVD-based Xinv computed', 2, 'df')
    call time_pull("DF - Computing X+")
end subroutine df_compute_Xinv_svd

    subroutine df_electrostatic_static(df, eel)
        !! Compute the electrostatic potential generated by MM static multipoles
        !! at the charge coordinates. Follows the same pattern as V_m2n in
        !! electrostatic_for_ene from mod_qm_helper.F90.

        use mod_electrostatics, only: potential_M2E, ommp_electrostatics_type
        use mod_memory, only: mallocate

        implicit none

        type(ommp_density_fit_type), intent(inout) :: df
        type(ommp_electrostatics_type), intent(in) :: eel

        call df_update(df)

        if(.not. df%V_m2q_done) then
            call time_push
            if(.not. allocated(df%V_m2q)) then
                call mallocate('df_electrostatic_static [V_mm2df]', &
                               df%n_charges, df%V_m2q)
            end if

            df%V_m2q = 0.0_rp
            call potential_M2E(eel, df%charge_coord, df%V_m2q)
            call time_pull("DF - Computing V M2Q")
            df%V_m2q_done = .true.
        end if
    end subroutine df_electrostatic_static

    subroutine df_electrostatic_dipoles(df, eel)
        !! Compute the electrostatic potential generated by MM induced dipoles
        !! at the charge coordinates. Follows the same pattern as V_p2n in
        !! electrostatic_for_ene from mod_qm_helper.F90.

        use mod_electrostatics, only: potential_D2E, ommp_electrostatics_type
        use mod_memory, only: mallocate

        implicit none

        type(ommp_density_fit_type), intent(inout) :: df
        type(ommp_electrostatics_type), intent(in) :: eel

        call df_update(df)

        if(eel%pol_atoms == 0) return

        if(.not. df%V_p2q_done .and. eel%ipd_done) then
            call time_push
            if(.not. allocated(df%V_p2q)) then
                call mallocate('df_electrostatic_dipoles [V_pd2df]', &
                               df%n_charges, eel%n_ipd, df%V_p2q)
            end if

            df%V_p2q = 0.0_rp
            if(eel%amoeba) then
                call potential_D2E(eel, df%charge_coord, df%V_p2q(:,_amoeba_D_))
                call potential_D2E(eel, df%charge_coord, df%V_p2q(:,_amoeba_P_), .true.)
            else
                call potential_D2E(eel, df%charge_coord, df%V_p2q(:,1))
            end if
            df%V_p2q_done = .true.
            call time_pull("DF - Computing V P2Q")
        end if

    end subroutine df_electrostatic_dipoles

    subroutine df_project_static(df, eel)
        !! Compute the projected static quantity: VXI_m = V_m2q @ Xinv.
        !! This is used to construct the Fock matrix elements.

        use mod_memory, only: mallocate
        use mod_electrostatics, only: ommp_electrostatics_type

        implicit none

        type(ommp_density_fit_type), intent(inout) :: df
        type(ommp_electrostatics_type), intent(in) :: eel

        call df_update(df)

        if(.not. df%VXI_m_done) then
            call time_push
            if(.not. df%xinv_done) then
                call df_compute_Xinv_svd(df)
            end if

            if(.not. df%V_m2q_done) then
                call df_electrostatic_static(df, eel)
            end if

            if(.not. allocated(df%VXI_m)) then
                call mallocate('df_project_static [VXI_m]', df%n_pts, df%VXI_m)
            end if

            df%VXI_m = 0.0_rp
            call dgemv('T', df%n_charges, df%n_pts, 1.0_rp, &
                       df%Xinv, df%n_charges, &
                       df%V_m2q, 1, 0.0_rp, df%VXI_m, 1)
            df%VXI_m_done = .true.
            call time_pull("DF - Computing VX+ (MM)")
        end if
    end subroutine df_project_static

    subroutine df_project_dipoles(df, eel)
        !! Compute the projected dipole quantity: VXI_p = V_p2q @ Xinv.
        !! This is used to construct the Fock matrix elements.

        use mod_memory, only: mallocate
        use mod_electrostatics, only: ommp_electrostatics_type

        implicit none

        type(ommp_density_fit_type), intent(inout) :: df
        type(ommp_electrostatics_type), intent(in) :: eel

        call df_update(df)

        if(eel%pol_atoms < 1) return

        if(.not. df%VXI_p_done .and. eel%ipd_done) then
            if(.not. df%Xinv_done) then
                call df_compute_Xinv_svd(df)
            end if

            if(.not. df%V_p2q_done) then
                call df_electrostatic_dipoles(df, eel)
            end if

            if(.not. allocated(df%VXI_p)) then
                call mallocate('df_project_dipoles [VXI_p]', df%n_pts, eel%n_ipd, df%VXI_p)
            end if
            call time_push

            df%VXI_p = 0.0_rp
            if(eel%amoeba) then
                call dgemv('T', df%n_charges, df%n_pts, 1.0_rp, &
                        df%Xinv, df%n_charges, &
                        df%V_p2q(:,_amoeba_D_), 1, 0.0_rp, df%VXI_p(:,_amoeba_D_), 1)
                call dgemv('T', df%n_charges, df%n_pts, 1.0_rp, &
                        df%Xinv, df%n_charges, &
                        df%V_p2q(:,_amoeba_P_), 1, 0.0_rp, df%VXI_p(:,_amoeba_P_), 1)
            else
                call dgemv('T', df%n_charges, df%n_pts, 1.0_rp, &
                        df%Xinv, df%n_charges, &
                        df%V_p2q(:,1), 1, 0.0_rp, df%VXI_p(:,1), 1)
            endif
            df%VXI_p_done = .true.

            call time_pull("DF - Computing VX+ (pol)")
        end if
    end subroutine df_project_dipoles

    subroutine df_e_field_to_pol(df, eel)
        !! Compute the electric field generated by fitted charges at
        !! polarizable sites. Analogous to E_n2p in mod_qm_helper.F90.

        use mod_memory, only: mallocate
        use mod_electrostatics, only: q_elec_prop, coulomb_kernel, ommp_electrostatics_type

        implicit none

        type(ommp_density_fit_type), intent(inout) :: df
        type(ommp_electrostatics_type), intent(in) :: eel

        real(rp) :: kernel(5), dr(3), tmpV, tmpE(3), tmpEgr(6), tmpHE(10)
        integer(ip) :: i, j

        call df_update(df)

        if(.not. df%E_q2p_done) then
            if(.not. allocated(df%E_q2p)) then
                call mallocate('df_e_field_to_pol [E_q2p]', &
                               3_ip, eel%pol_atoms, df%E_q2p)
            end if

            df%E_q2p = 0.0_rp
            !$omp parallel do private(i, j, dr, kernel, tmpE) &
            !$omp schedule(static)
            do j = 1, eel%pol_atoms
                do i = 1, df%n_charges
                    dr = eel%cpol(:,j) - df%charge_coord(:,i)
                    call coulomb_kernel(dr, 1, kernel)

                    tmpE = 0.0
                    call q_elec_prop(df%target_charges(i), dr, kernel, &
                                     .false., tmpV, &
                                     .true., tmpE, &
                                     .false., tmpEgr, &
                                     .false., tmpHE)

                    df%E_q2p(:,j) = df%E_q2p(:,j) + tmpE
                end do
            end do
            df%E_q2p_done = .true.
        end if
    end subroutine df_e_field_to_pol

    subroutine df_rotation_geomgrad(df, eel, mmg)

        use mod_memory, only: mallocate
        use mod_electrostatics, only: q_elec_prop, coulomb_kernel, ommp_electrostatics_type

        implicit none

        type(ommp_density_fit_type), intent(inout) :: df
        type(ommp_electrostatics_type), intent(in) :: eel
        real(rp), intent(inout) :: mmg(:,:)

        real(rp) :: kernel(5), dr(3), tmpV, tmpE(3), tmpEgr(6), tmpHE(10)
        integer(ip) :: i, j, n_mm

        if(.not. eel%amoeba) return
        call df_update(df)

        if(.not. df%E_q2M_done) then
            n_mm = eel%top%mm_atoms

            if(.not. allocated(df%E_q2M)) then
                call mallocate('df_e_field_q2M [E_q2M]', 3_ip, n_mm, df%E_q2M)
            end if

            if(.not. allocated(df%GEF_q2M)) then
                call mallocate('df_e_field_q2M [GEF_q2M]', 6_ip, n_mm, df%GEF_q2M)
            end if

            df%E_q2M = 0.0_rp
            df%GEF_q2M = 0.0_rp


            !$omp parallel do default(shared) schedule(static) &
            !$omp private(i,j,dr,kernel,tmpE,tmpV,tmpEgr,tmpHE)
            do j = 1, n_mm
                do i = 1, df%n_charges
                    dr = eel%top%cmm(:,j) - df%charge_coord(:,i)
                    call coulomb_kernel(dr, 2, kernel)

                    tmpV = 0.0_rp
                    tmpE = 0.0_rp
                    tmpEgr = 0.0_rp
                    tmpHE = 0.0_rp

                    call q_elec_prop(df%target_charges(i), dr, kernel, &
                                     .false., tmpV, &
                                     .true., tmpE, &
                                     .true., tmpEgr, &
                                     .false., tmpHE)

                    df%E_q2M(:,j) = df%E_q2M(:,j) + tmpE
                    df%GEF_q2M(:,j) = df%GEF_q2M(:,j) + tmpEgr
                end do
            end do

            df%E_q2M_done = .true.
            df%GEF_q2M_done = .true.
        end if

        call rotation_geomgrad(eel, df%E_q2m, df%GEF_q2M, mmg)

    end subroutine df_rotation_geomgrad

    subroutine df_e_field_pol_ene(df, eel)
        !! Compute the polarization energy from the electric field generated
        !! by fitted charges: E = -0.5 * sum_j ipd(:,j) .dot. E_q2p(:,j).
        !! Follows the same pattern as energy_MM_pol in mod_electrostatics.F90.

        use mod_memory, only: mallocate
        use mod_electrostatics, only: ommp_electrostatics_type

        implicit none

        type(ommp_density_fit_type), intent(inout) :: df
        type(ommp_electrostatics_type), intent(in) :: eel

        real(rp) :: eMM
        integer(ip) :: i, j

        call df_update(df)

        if(.not. df%E_pol_ene_done) then
            !! Ensure E_q2p is computed first
            if(.not. df%E_q2p_done) then
                call df_e_field_to_pol(df, eel)
            end if

            eMM = 0.0_rp
            ! if(eel%amoeba) then
                !! Use only _amoeba_D_ dipoles, contracted with E_q2p
            !    do i = 1, 3
            !        do j = 1, eel%pol_atoms
            !           eMM = eMM - eel%ipd(i,j,_amoeba_D_) * df%E_q2p(i,j)
            !        end do
            !    end do
            !else
            !    do i = 1, 3
            !        do j = 1, eel%pol_atoms
            !            eMM = eMM - eel%ipd(i,j,1) * df%E_q2p(i,j)
            !        end do
            !    end do
            !end if

            if(eel%amoeba) then
                !$omp parallel do default(shared) schedule(static) &
                !$omp private(i) reduction(+:eMM)
                do i=1, eel%pol_atoms
                    eMM = eMM - dot_product(eel%ipd(:,i,_amoeba_D_), df%E_q2p(:,i))
                end do
            else
                !$omp parallel do default(shared) schedule(static) &
                !$omp private(i) reduction(+:eMM)
                do i=1, eel%pol_atoms
                    eMM = eMM - dot_product(eel%ipd(:,i,1), df%E_q2p(:,i))
                end do
            end if

            df%E_pol_ene = 0.5_rp * eMM
            df%E_pol_ene_done = .false.
        end if
    end subroutine df_e_field_pol_ene



    subroutine df_electrostatics_for_geomgrad(df, eel)
        !! Compute electric field quantities needed for geometry gradients.
        !! 1. E_m2q : electric field of MM static multipoles at charge sites
        !! 2. E_p2q : electric field of IPDs at charge sites (per atom, per component)

        use mod_memory, only: mallocate
        use mod_electrostatics, only: q_elec_prop, mu_elec_prop, quad_elec_prop, &
                                      coulomb_kernel, ommp_electrostatics_type

        implicit none

        type(ommp_density_fit_type), intent(inout) :: df
        type(ommp_electrostatics_type), intent(in) :: eel

        integer(ip) :: i, j, a, k, n_cpt, n_mm, n_pts, n_charges, n_pol, ipol
        real(rp) :: kernel(5), dr(3), tmpV, tmpE(3), tmpEgr(6), tmpHE(10)

        call df_update(df)

        !! Part 1a: E_m2q - static multipole field from all MM atoms

        if(.not. df%E_m2q_done) then
            n_mm = eel%top%mm_atoms
            n_cpt = df%n_charges

            if(.not. allocated(df%E_m2q)) then
                call mallocate('[df_electrostatics_for_geomgrad] E_MM2Q', &
                               3, n_mm, df%n_charges, df%E_m2q)
            end if

            df%E_m2q = 0.0_rp

            if(eel%amoeba) then
                !$omp parallel do default(shared) schedule(dynamic) &
                !$omp private(i,j,dr,kernel,tmpV,tmpE,tmpEgr,tmpHE)
                do i = 1, n_mm
                    do j = 1, n_cpt
                        dr = df%charge_coord(:,j) - eel%top%cmm(:,i)
                        call coulomb_kernel(dr, 4, kernel)
                        tmpE = 0.0_rp

                        call q_elec_prop(eel%q(1,i), dr, kernel, &
                                         .false., tmpV, &
                                         .true., tmpE, &
                                         .false., tmpEgr, &
                                         .false., tmpHE)
                        call mu_elec_prop(eel%q(2:4,i), dr, kernel, &
                                          .false., tmpV, &
                                          .true., tmpE, &
                                          .false., tmpEgr, &
                                          .false., tmpHE)
                        call quad_elec_prop(eel%q(5:10,i), dr, kernel, &
                                            .false., tmpV, &
                                            .true., tmpE, &
                                            .false., tmpEgr, &
                                            .false., tmpHE)

                        df%E_m2q(:,i,j) = df%E_m2q(:,i,j) + tmpE
                    end do
                end do
            else
                !$omp parallel do default(shared) schedule(dynamic) &
                !$omp private(i,j,dr,kernel,tmpV,tmpE,tmpEgr,tmpHE)
                do i = 1, n_mm
                    do j = 1, n_cpt
                        dr = df%charge_coord(:,j) - eel%top%cmm(:,i)
                        call coulomb_kernel(dr, 2, kernel)
                        tmpE = 0.0_rp

                        call q_elec_prop(eel%q(1,i), dr, kernel, &
                                         .false., tmpV, &
                                         .true., tmpE, &
                                         .false., tmpEgr, &
                                         .false., tmpHE)

                        df%E_m2q(:,i,j) = df%E_m2q(:,i,j) + tmpE
                    end do
                end do
            end if

            df%E_m2q_done = .true.
        end if

        !! Part 1b: E_p2q - IPD field from polarizable atoms at charge sites
        !! Shape: (3, n_mm, n_charges, n_ipd), only populated where mm_polar(i) > 0

        if(.not. df%E_p2q_done .and. eel%ipd_done) then
            n_mm = eel%top%mm_atoms
            n_cpt = df%n_charges
            n_pol = eel%pol_atoms

            if(n_pol > 0) then
                if(.not. allocated(df%E_p2q)) then
                    call mallocate('[df_electrostatics_for_geomgrad] E_p2q', &
                                   3, n_mm, df%n_charges, eel%n_ipd, df%E_p2q)
                end if

                df%E_p2q = 0.0_rp

                do k = 1, eel%n_ipd
                    !$omp parallel do collapse(2) default(shared) schedule(static) &
                    !$omp private(i,ipol,j,dr,kernel,tmpV,tmpE,tmpEgr,tmpHE)
                    do i = 1, n_pol
                        ipol = eel%polar_mm(i)
                        do j = 1, n_cpt
                            dr = df%charge_coord(:,j) - eel%cpol(:,i)
                            call coulomb_kernel(dr, 2, kernel)
                            tmpE = 0.0_rp
                            call mu_elec_prop(eel%ipd(:,i,k), dr, kernel, &
                                                .false., tmpV, &
                                                .true., tmpE, &
                                                .false., tmpEgr, &
                                                .false., tmpHE)
                            df%E_p2q(:,ipol,j,k) = df%E_p2q(:,ipol,j,k) + tmpE
                        end do
                    end do
                end do
                df%E_p2q_done = .true.
            end if
        end if

        !! Part 2: dX/drfit and dXinv/drfit
        if(.not. df%dX_dr_done) then
            n_charges = df%n_charges
            n_pts     = df%n_pts

            if(.not. allocated(df%dX_dr)) then
                call mallocate('[df_electrostatics_for_geomgrad] dX_dr', &
                               3_ip, n_pts, n_charges, df%dX_dr)
            end if

            df%dX_dr = 0.0_rp

            !$omp parallel do default(shared) schedule(dynamic) &
            !$omp private(i,j,dr,kernel)
            do i = 1, n_charges
                do j = 1, n_pts
                    dr = df%charge_coord(:,i) - df%fit_point_coord(:,j)
                    call coulomb_kernel(dr, 2, kernel)

                    df%dX_dr(:,j,i) = kernel(2) * dr
                end do
            end do

            df%dX_dr_done = .true.
        end if

        if(.not. df%V_m2q_done) call df_project_static(df, eel)
        if(.not. df%V_p2q_done) call df_project_dipoles(df, eel)

    end subroutine df_electrostatics_for_geomgrad

    subroutine df_geomgrad(df, eel, qmg, mmg, doqm, domm, ef_qm2fit)
        !! Compute the gradient (force) contribution from density fitting
        !! with respect to nuclear coordinates.
        !!
        !! The forces are accumulated into qmg and mmg arrays.
        !! The electric field at grid points and the Lagrange multipliers
        !! are computed internally.
        !!
        !! doqm / domm control which sub-blocks are updated; they allow
        !! the caller to compose this gradient with other contributors
        !! (e.g. QM forces) without double-counting.
        !!
        !! ef(3, n_pts) : electric field at fitting points (input)

        use mod_io, only: fatal_error, print_matrix
        use mod_electrostatics, only: ommp_electrostatics_type

        implicit none

#define USE_OPTIMIZED

        type(ommp_density_fit_type), intent(inout) :: df
        type(ommp_electrostatics_type), intent(in) :: eel
        real(rp), intent(inout) :: qmg(:,:)   ! (3, n_qm_atoms)
        real(rp), intent(inout) :: mmg(:,:)   ! (3, n_mm_atoms)
        logical, intent(in) :: doqm
        logical, intent(in) :: domm
        real(rp), intent(in) :: ef_qm2fit(3,df%n_pts)

        integer(ip) :: i, a, j, k, kk
        integer(ip) :: n_qm, n_mm, n_q, n_fit
#ifndef USE_OPTIMIZED
        real(rp), allocatable :: XtX_inv(:,:)
        integer(ip) :: b, l, s
#endif
        real(rp), allocatable :: resids(:)
        real(rp), allocatable :: V_mmpol2q(:)
        real(rp), allocatable :: field_combined(:,:,:)
        real(rp), allocatable :: w1(:), w4(:), vi(:)
        real(rp), allocatable :: tmp(:,:,:), tmp2(:,:)
        real(rp), allocatable :: w2(:,:), inner(:,:)

        if(.not. df%initialized) then
            call fatal_error("df_geomgrad: density fit object not initialized.")
        end if

        call df_update(df)

        call df_electrostatics_for_geomgrad(df, eel)

        call compute_nabla_matrices(df)

        if(doqm) then
            qmg= 0.0_rp
        end if
        if(domm) then
            mmg = 0.0_rp
        end if

        call df_rotation_geomgrad(df, eel, mmg)

        n_mm = df%mm_top%mm_atoms
        n_qm = df%qm_top%mm_atoms
        n_q = df%n_charges
        n_fit = df%n_pts

#ifndef USE_OPTIMIZED
        call ommp_message("Using multiple loop DF gradients", OMMP_VERBOSE_DEBUG, 'df')
#define USE_LOOPS
#else
        call ommp_message("Using optimized DF gradients", OMMP_VERBOSE_DEBUG, 'df')
#endif

        !! Allocate local intermediate quantities
#ifdef USE_LOOPS
        call mallocate('[df_geomgrad] XtX_inv', n_q, n_q, XtX_inv)
#endif
        if(n_q < n_fit) then
            call mallocate('[df_geomgrad] resids', n_fit, resids)
            call mallocate('[df_geomgrad] inner', 3, n_fit, inner)
            call mallocate('[df_geomgrad] w4', n_q, w4)
        else if(n_q > n_fit) then
            call mallocate('[df_geomgrad] resids', n_q, resids)
            call mallocate('[df_geomgrad] w4', n_fit, w4)
        end if
        call mallocate('[df_geomgrad] V_mmpol2q', n_q, V_mmpol2q)
        call mallocate('[df_geomgrad] field_combined', 3, n_mm, n_q, field_combined)
        
        call mallocate('[df_geomgrad] w1', n_fit, w1)
        call mallocate('[df_geomgrad] w2', 3, n_fit, w2)

        call mallocate('[df_geomgrad] vi', n_fit, vi)
        call mallocate('[df_geomgrad] tmp', 3_ip, n_q, n_fit, tmp)
        call mallocate('[df_geomgrad] tmp2', 3_ip, n_q, tmp2)

        if(.not. (df%nabla_g_qm_is_null .and. &
                  df%nabla_q_mm_is_null)) then
            call fatal_error("Only grid dependent on MM atoms only and &
                             &charges dependent on QM atoms only is currently supported")
        end if

        if(eel%pol_atoms > 0) then
            if(eel%amoeba) then
                field_combined(:,:,:) = df%E_m2q(:,:,:) + &
                                        (df%E_p2q(:,:,:,_amoeba_D_) + &
                                        df%E_p2q(:,:,:,_amoeba_P_)) * 0.5
                V_mmpol2q = df%V_m2q + &
                            (df%V_p2q(:,_amoeba_D_) + &
                            df%V_p2q(:,_amoeba_P_)) * 0.5
                vi = df%VXI_m + 0.5_rp * df%VXI_p(:,_amoeba_D_) + &
                     0.5_rp * df%VXI_p(:,_amoeba_P_)
            else
                field_combined(:,:,:) = df%E_m2q(:,:,:) + df%E_p2q(:,:,:,1)
                V_mmpol2q = df%V_m2q + df%V_p2q(:,1)
                vi = df%VXI_m + df%VXI_p(:,1)
            end if
        else
            field_combined(:,:,:) = df%E_m2q(:,:,:)
            V_mmpol2q(:) = df%V_m2q(:)
            vi = df%VXI_m
        end if

        ! On MM atoms
        ! 1st term: static field dVmm(rq)/drmm @ qfit
#ifdef USE_LOOPS
        do i=1, n_mm
            do j=1, n_q
                do a=1, 3
                    ! TODO AMOEBA!!!
                    if(eel%pol_atoms>0) then
                        if(eel%amoeba) then
                            mmg(a,i) = mmg(a,i) + (df%E_m2q(a,i,j)  + 0.5 * df%E_p2q(a,i,j,1) + 0.5 * df%E_p2q(a,i,j,2) ) &
                                        * df%target_charges(j)
                        else
                            call fatal_error("Polarizable AMBER NOT IMPLEMENTED")
                        end if
                    else
                        mmg(a,i) = mmg(a,i) + df%E_m2q(a,i,j) * df%target_charges(j) 
                    end if
                end do
            end do
        end do
#endif

#ifdef USE_OPTIMIZED
        call dgemv('N', 3*n_mm, n_q, 1.0_rp, &
                    field_combined, 3*n_mm, &
                    df%target_charges, 1, 1.0_rp, &
                    mmg, 1)
        if(.not. df%nabla_q_mm_is_null) then
            call fatal_error("Depndency between MM coordinates and charges points is not supported.")
            ! TODO If r_charges depend on MM atoms, another therm should be added here!
        end if
#endif
        ! 2nd term Vmm(rq) dA+/drmm Vqm(rfit)

#ifdef USE_LOOPS
        XtX_inv = 0.0
        do i=1, n_q
            do j=1, n_fit
                do k=1, n_q
                    XtX_inv(i,k) = XtX_inv(i,k) + df%Xinv(i,j) * df%Xinv(k,j)
                end do
            end do
        end do

        resids = 0.0
        do i=1, n_fit
            do j=1, n_q
                resids(i) = resids(i) - df%X(i,j) * df%target_charges(j)
            end do
            resids(i) = resids(i) + df%fit_potential(i) 
        end do
#endif

#ifdef USE_OPTIMIZED
        if(n_q < n_fit) then
            resids = df%fit_potential
            call dgemv('N', n_fit, n_q, -1.0_rp, df%X, n_fit, df%target_charges, 1, 1.0_rp, resids, 1)
        else if(n_q > n_fit) then
            resids = V_mmpol2q
            call dgemv('T', n_fit, n_q, -1.0_rp, df%X, n_fit, vi, 1, 1.0_rp, resids, 1)
        end if
#endif

#ifdef USE_LOOPS
        ! dq_dm_term2 = np.einsum('li,ija,j->ial', A_inv, dA_dr, qfit)
        do i=1, n_mm
            do a=1, 3
                do l=1, n_q
                    do s=1, n_q
                        mmg(a,i) = mmg(a,i) - df%Xinv(l,i) * df%dX_dr(a,i,s) * df%target_charges(s) * V_mmpol2q(l)
                    end do
                end do
            end do
        end do
        
        ! dq_dm_term3 = -np.einsum('lk,ika,i->ial', AtA_inv, dA_dr, resids)
        do i=1, n_mm
            do a=1, 3
                do l=1, n_q
                    do s=1, n_q
                        mmg(a,i) = mmg(a,i) + XtX_inv(l,s) * df%dX_dr(a,i,s) * resids(i) * V_mmpol2q(l)
                    end do
                end do
            end do
        end do 
#endif

#ifdef USE_OPTIMIZED
        ! w1 = Xinv^T * V_mmpol2q
        call dgemv('T', n_q, n_fit, 1.0_rp, df%Xinv, n_q, V_mmpol2q, 1, 0.0_rp, w1, 1)

        ! w2 = dX_dr * target_charges
        call dgemv('N', 3*n_fit, n_q, 1.0_rp, df%dX_dr, 3*n_fit, df%target_charges, 1, 0.0_rp, w2, 1)

        if(df%nabla_g_mm_is_identity .and. df%nabla_q_mm_is_null) then
            !$omp parallel do collapse(2) default(shared)
            do i = 1, n_mm
                do a = 1, 3
                    mmg(a,i) = mmg(a,i) - w1(i) * w2(a,i)
                end do
            end do

            if(n_q < n_fit) then

                ! w4 = Xinv * w1  (same as XtX_inv^T*V_mmpol2q, but O(n_q*n_fit))
                call dgemv('N', n_q, n_fit, 1.0_rp, df%Xinv, n_q, w1, 1, 0.0_rp, w4, 1)

                ! inner = dX_dr * w4
                call dgemv('N', 3*n_fit, n_q, 1.0_rp, df%dX_dr, 3*n_fit, w4, 1, 0.0_rp, inner, 1)

                !$omp parallel do collapse(2) default(shared)
                do i = 1, n_mm
                    do a = 1, 3
                        mmg(a,i) = mmg(a,i) + resids(i) * inner(a,i)
                    end do
                end do
            else if(n_q > n_fit) then
                ! w1 = X_inv^T q
                call dgemv('T', n_q, n_fit, 1.0_rp, df%Xinv, n_q, df%target_charges, 1, 0.0_rp, w4, 1)

                do i=1, n_mm
                    do a=1, 3
                        do j=1, n_q
                            mmg(a, i) = mmg(a,i) + df%dX_dr(a,i,j) * w4(i) * resids(j)
                        end do
                    end do
                end do
            end if
        else
            if(.not. df%nabla_g_mm_is_identity) then
                call fatal_error("Unsupported dependency between grid points and MM atoms")
            end if

            if(.not. df%nabla_q_mm_is_null) then
                call fatal_error("Unsupported dependency between fit points and MM atoms")
            end if
        end if
#endif 

        !! 3rd term Vmm(rq) A+ Eqm(rfit)

#ifdef USE_LOOPS
        if(eel%pol_atoms > 0) then
            if(eel%amoeba) then
                do i=1, n_mm
                    do a=1, 3
                        mmg(a,i) = mmg(a,i) - (df%VXI_m(i) + 0.5 * df%VXI_p(i,_amoeba_D_) &
                        + 0.5 * df%VXI_p(i,_amoeba_P_)) * ef_qm2fit(a,i)
                    end do
                end do
            else
                do i=1, n_mm
                    do a=1, 3
                        mmg(a,i) = mmg(a,i) - (df%VXI_m(i) + df%VXI_p(i,1)) * ef_qm2fit(a,i)
                    end do
                end do
            end if
        else
            do i=1, n_mm
                do a=1, 3
                    mmg(a,i) = mmg(a,i) - df%VXI_m(i) * ef_qm2fit(a,i)
                end do
            end do
        end if
#endif

#ifdef USE_OPTIMIZED
        ! TODO
        !$omp parallel do collapse(2) default(shared) schedule(static)
        do a = 1, 3
            do i = 1, n_fit
                mmg(a,i) = mmg(a,i) - vi(i) * ef_qm2fit(a,i)
            end do
        end do     
#endif

        ! On QM atoms
        ! 1st term: static field dVmm(rq)/drqm @ qfit
#ifdef USE_LOOPS
        do i=1, n_mm
            do j=1, n_q
                do a=1,3
                    if(eel%pol_atoms>0) then
                        if(eel%amoeba) then
                            qmg(a,j) = qmg(a,j) - (df%E_m2q(a,i,j)  + 0.5* df%E_p2q(a,i,j,1) &
                                        + 0.5 * df%E_p2q(a,i,j,2) )  * df%target_charges(j)
                        else
                            call fatal_error("Polarizable Amber not implemented")
                        end if
                    else
                        qmg(a,j) = qmg(a,j) - df%E_m2q(a,i,j)  * df%target_charges(j)
                    end if
                end do
            end do
        end do
#endif

#ifdef USE_OPTIMIZED
        if(df%nabla_q_qm_is_identity) then
            !$omp parallel do collapse(2) default(shared) &
            !$omp private(i)
            do j = 1, n_qm
                do a = 1, 3
                    do i = 1, n_mm
                        qmg(a,j) = qmg(a,j) - field_combined(a,i,j) * df%target_charges(j)
                    end do
                end do
            end do
        else if(df%nabla_q_qm_is_sparse) then
            !$omp parallel do collapse(2) default(shared) &
            !$omp private(i, k)
            do j = 1, n_qm
                do a = 1, 3
                    do k=df%nabla_q_qm_sparse%ci(df%nabla_q_qm_sparse%ri(j)), &
                         df%nabla_q_qm_sparse%ci(df%nabla_q_qm_sparse%ri(j+1)-1)
                        do i = 1, n_mm
                            qmg(a,j) = qmg(a,j) - field_combined(a,i,k) * df%target_charges(k)
                        end do
                    end do
                end do
            end do
        else if(.not. df%nabla_q_qm_is_null) then
            call fatal_error("Dense nabala Q-QM is still not implemented.")
        end if
#endif

        ! 2nd term on QM atoms
        ! Part 1: sum_j V_mmpol2q(j) * Xinv(j,k) * dX_dr(a,k,i) * target_charges(i)
        ! Part 2: -sum_j V_mmpol2q(j) * XtX_inv(j,i) * dX_dr(a,k,i) * resids(k)

#ifdef USE_LOOPS
        ! Part 1
         do i=1, n_q
             do a=1,3
                 do j=1, n_q
                     do k=1, n_fit
                         qmg(a,i) = qmg(a,i) + df%Xinv(j,k) * df%dX_dr(a,k,i) * df%target_charges(i) * V_mmpol2q(j)
                     end do
                 end do
             end do
         end do

        ! Part 2
         do i=1, n_q
             do a=1,3
                 do j=1, n_q
                     do k=1, n_fit
                         qmg(a,i) = qmg(a,i) - XtX_inv(j,i) * df%dX_dr(a,k,i) * resids(k) * V_mmpol2q(j)
                     end do
                 end do
             end do
         end do
#endif

#ifdef USE_OPTIMIZED
        call dgemv('T', n_q, n_fit, 1.0_rp, &
                    df%Xinv, n_q, V_mmpol2q, 1, 0.0_rp, w1, 1)

        !$omp parallel do collapse(3) schedule(static)
        do a = 1, 3
            do i = 1, n_q
                do k = 1, n_fit
                    tmp(a, i, k) = df%dX_dr(a, k, i)
                end do
            end do
        end do

        call dgemv('N', 3*n_q, n_fit, 1.0_rp, &
                    tmp, 3*n_q, w1, 1, 0.0_rp, tmp2, 1)


        if(df%nabla_q_qm_is_identity) then
            !$omp parallel do collapse(2) schedule(static)
            do i = 1, n_qm
                do a = 1, 3
                    qmg(a, i) = qmg(a, i) + df%target_charges(i) * tmp2(a, i)
                end do
            end do
        else if(df%nabla_q_qm_is_sparse) then
            !$omp parallel do collapse(2) schedule(static) private(k)
            do i = 1, n_qm
                do a = 1, 3
                    do k=df%nabla_q_qm_sparse%ci(df%nabla_q_qm_sparse%ri(i)), &
                         df%nabla_q_qm_sparse%ci(df%nabla_q_qm_sparse%ri(i+1)-1)
                        qmg(a, i) = qmg(a, i) + df%target_charges(k) * tmp2(a, k)
                    end do
                end do
            end do
        else if(.not. df%nabla_q_qm_is_null) then
            call fatal_error("Dense nabala Q-QM is still not implemented.")
        end if

        if(n_q < n_fit) then

            call dgemv('N', 3*n_q, n_fit, 1.0_rp, &
                        tmp, 3*n_q, resids, 1, 0.0_rp, tmp2, 1)
            

            if(df%nabla_q_qm_is_identity) then
                !$omp parallel do collapse(2) schedule(static)
                do i = 1, n_qm
                    do a = 1, 3
                        qmg(a, i) = qmg(a, i) - w4(i) * tmp2(a, i)
                    end do
                end do
            else if(df%nabla_q_qm_is_sparse) then
                !$omp parallel do collapse(2) schedule(static) private(k)
                do i = 1, n_qm
                    do a = 1, 3
                        do k=df%nabla_q_qm_sparse%ci(df%nabla_q_qm_sparse%ri(i)), &
                            df%nabla_q_qm_sparse%ci(df%nabla_q_qm_sparse%ri(i+1)-1)
                            qmg(a, i) = qmg(a, i) - w4(k) * tmp2(a, k)
                        end do
                    end do
                end do
            else if(.not. df%nabla_q_qm_is_null) then
                call fatal_error("Dense nabala Q-QM is still not implemented.")
            end if
        else if(n_q > n_fit) then
            if(df%nabla_q_qm_is_identity) then
                !$omp parallel do collapse(2) schedule(static)
                do j=1, n_q
                    do a=1, 3
                        do i=1, n_mm
                            qmg(a, j) = qmg(a,j) - df%dX_dr(a,i,j) * w4(i) * resids(j)
                        end do
                    end do
                end do
            else if(df%nabla_q_qm_is_sparse) then
                !$omp parallel do collapse(2) schedule(static) private(k)
                do j=1, n_qm
                    do a=1, 3
                        do k=df%nabla_q_qm_sparse%ci(df%nabla_q_qm_sparse%ri(j)), &
                            df%nabla_q_qm_sparse%ci(df%nabla_q_qm_sparse%ri(j+1)-1)
                            do i=1, n_mm
                                qmg(a, j) = qmg(a,j) - df%dX_dr(a,i,k) * w4(i) * resids(k)
                            end do
                        end do
                    end do
                end do
            else if(.not. df%nabla_q_qm_is_null) then
                call fatal_error("Dense nabala Q-QM is still not implemented.")
            end if
        end if
#endif

        ! Third term is computed outside!

        !! Deallocate local intermediate quantities
#ifdef USE_LOOPS
        call mfree('[df_geomgrad] XtX_inv', XtX_inv)
#endif
        call mfree('[df_geomgrad] resids', resids)
        call mfree('[df_geomgrad] V_mmpol2q', V_mmpol2q)
        call mfree('[df_geomgrad] field_combined', field_combined)
        call mfree('[df_geomgrad] w1', w1)
        call mfree('[df_geomgrad] w4', w4)
        call mfree('[df_geomgrad] w2', w2)
        call mfree('[df_geomgrad] inner', inner)
        call mfree('[df_geomgrad] vi', vi)
        call mfree('[df_geomgrad] tmp', tmp)
        call mfree('[df_geomgrad] tmp2', tmp2)
    end subroutine df_geomgrad


end module
