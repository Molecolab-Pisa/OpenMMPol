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
                             ommp_df_mm_top, ommp_df_qm_top
    use mod_io, only: fatal_error, ommp_message
    use mod_topology, only: ommp_topology_type

    implicit none
    private

    type ommp_density_fit_type
        integer(ip) :: n_pts = 0
        !! Number of fitting points

        integer(ip) :: n_charges = 0
        !! Number of target charges

        integer(ip) :: df_method = ommp_df_solver_svd
        !! Method used for solving the linear system (PINV, SVD, or NORMAL)

        real(rp) :: svd_rcond = ommp_df_svd_rcond_default
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

        real(rp), allocatable :: E_MM2Q(:,:,:)
        !! Electric field generated by MM static multipoles at charge sites.
        !! Shape: (n_charges, 3, n_mm_atoms) — E_MM2Q(:,k,i) is the k-th
        !! Cartesian component of the field at charge site i due to MM atom i.

        logical(lp) :: E_MM2Q_done = .false.
        !! Flag indicating whether E_MM2Q has been computed

        real(rp), allocatable :: dX_dr(:,:,:)
        !! Derivative of the design matrix X w.r.t. fit-point coordinates.
        !! Shape: (n_charges, 3, n_pts) — dX_dr(i,k,j) is the k-th
        !! Cartesian component of dX/d(r_grid)_j for charge i.
        !! For Coulomb kernel X(i,j) = 1/|r_charge_i - r_grid_j|, we have:
        !!   dX/dr_grid_j = kernel(2) * dr
        !! where dr = r_grid_j - r_charge_i.

        logical(lp) :: dX_dr_done = .false.
        !! Flag indicating whether dX_dr has been computed

        real(rp), allocatable :: lambda(:)
        !! Lagrange multipliers: lambda(s) = sum_l Xinv(l,s) * V_m2q(l)
        !! Used in the gradient computation [ d/dr A ] . lambda . q

        logical(lp) :: lambda_done = .false.
        !! Flag indicating whether lambda has been computed

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

        real(rp), allocatable :: nabla_g_qm(:,:)
        !! Gradient matrix: grid wrt QM coordinates
        logical(lp) :: nabla_g_qm_is_null = .false.
        !! .true. if nabla_g_qm is the null matrix
        logical(lp) :: nabla_g_qm_is_identity = .false.
        !! .true. if nabla_g_qm is the identity matrix

        real(rp), allocatable :: nabla_q_qm(:,:)
        !! Gradient matrix: q wrt QM coordinates
        logical(lp) :: nabla_q_qm_is_null = .false.
        !! .true. if nabla_q_qm is the null matrix
        logical(lp) :: nabla_q_qm_is_identity = .false.
        !! .true. if nabla_q_qm is the identity matrix

        real(rp), allocatable :: nabla_q_mm(:,:)
        !! Gradient matrix: q wrt MM coordinates
        logical(lp) :: nabla_q_mm_is_null = .false.
        !! .true. if nabla_q_mm is the null matrix
        logical(lp) :: nabla_q_mm_is_identity = .false.
        !! .true. if nabla_q_mm is the identity matrix

    end type ommp_density_fit_type

    public :: ommp_density_fit_type
    public :: df_init, df_terminate
    public :: df_solve, df_generate_grid
    public :: df_electrostatic_static, df_electrostatic_dipoles
    public :: df_project_static, df_project_dipoles
    public :: df_e_field_to_pol, df_e_field_pol_ene
    public :: df_compute_lambda
    public :: compute_nabla_matrices, df_geomgrad
    public :: df_electrostatics_for_geomgrad
    !! Accessors for gradient-related quantities (exposed outside)

contains

    subroutine generate_grid_from_topo(top, out_coord, n_tot, n_pts_per_atom, radius, grid_type, label)
        !! Generic grid generation from a topology.
        !!
        !! Supported grid_type values (from mod_constants):
        !!   1 (ommp_df_atoms) — atom positions
        !!   2 (ommp_df_fibonacci) — fibonacci sphere per atom
        !!   3 (ommp_df_cubic) — 7-point cubic grid (1 center + 6 face centers)

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

    end subroutine generate_grid_from_topo

    subroutine compute_nabla_matrices(df)

        implicit none

        type(ommp_density_fit_type), intent(inout) :: df
        integer(ip) :: i

        select case(df%charge_top_type)
        case(ommp_df_qm_top)
            
            df%nabla_q_mm_is_null = .true.
            df%nabla_q_mm_is_identity = .false.

            select case(df%charge_point_type)
            case(ommp_df_atoms)
                df%nabla_q_qm_is_null = .false.
                df%nabla_q_qm_is_identity = .true.

            case(ommp_df_fibonacci)
                
                df%nabla_q_qm_is_null = .false.
                df%nabla_q_qm_is_identity = .false.

                if(.not. allocated(df%nabla_q_qm)) then
                    call mallocate('compute_nabla_matrices [nabla_q_qm]', 3 * df%n_charges, 3 * df%qm_top%mm_atoms, df%nabla_q_qm)
                end if
                df%nabla_q_qm = 0.0

                do i=1, df%charge_top%mm_atoms
                    df%nabla_q_qm((i-1)*3*df%charge_n_pts_per_atom+1:(i)*3*df%charge_n_pts_per_atom, (i-1)*3+1:i*3) = 1.0_rp
                end do

            case(ommp_df_cubic)
                
                df%nabla_q_qm_is_null = .false.
                df%nabla_q_qm_is_identity = .false.

                if(.not. allocated(df%nabla_q_qm)) then
                    call mallocate('compute_nabla_matrices [nabla_q_qm]', 3 * df%n_charges, 3 * df%qm_top%mm_atoms, df%nabla_q_qm)
                end if
                df%nabla_q_qm = 0.0

                do i=1, df%charge_top%mm_atoms
                    df%nabla_q_qm((i-1)*3*7+1:(i)*3*7, (i-1)*3+1:i*3) = 1.0_rp
                end do

            case default
                call fatal_error('compute_nabla_matrices: unknown charge_point_type.')
            end select
        case(ommp_df_mm_top)
            call fatal_error("Unsupported case: charge points generated from MM topology.")
        case default
            call fatal_error('compute_nabla_matrices: unknown charge_top_type.')
        end select

        select case(df%fit_top_type)
        case(ommp_df_mm_top)
            ! Obvius case fit grid is built on mm atoms
            
            df%nabla_g_qm_is_null = .true.
            df%nabla_g_qm_is_identity = .false.

            select case(df%fit_point_type)
            case(ommp_df_atoms)
                df%nabla_g_mm_is_null = .false.
                df%nabla_g_mm_is_identity = .true.

            case(ommp_df_fibonacci)
                
                df%nabla_g_mm_is_null = .false.
                df%nabla_g_mm_is_identity = .false.

                if(.not. allocated(df%nabla_g_mm)) then
                    call mallocate('compute_nabla_matrices [nabla_g_mm]', 3 * df%n_pts, 3 * df%mm_top%mm_atoms, df%nabla_g_mm)
                end if
                df%nabla_g_mm = 0.0

                do i=1, df%fit_top%mm_atoms
                    df%nabla_g_mm((i-1)*3*df%fit_n_pts_per_atom:(i)*3*df%fit_n_pts_per_atom, (i-1)*3:i*3) = 1.0_rp
                end do

            case(ommp_df_cubic)
                
                df%nabla_g_mm_is_null = .false.
                df%nabla_g_mm_is_identity = .false.

                if(.not. allocated(df%nabla_g_mm)) then
                    call mallocate('compute_nabla_matrices [nabla_g_mm]', 3 * df%n_pts, 3 * df%mm_top%mm_atoms, df%nabla_g_mm)
                end if
                df%nabla_g_mm = 0.0

                do i=1, df%fit_top%mm_atoms
                    df%nabla_g_mm((i-1)*3*7+1:(i)*3*7, (i-1)*3+1:i*3) = 1.0_rp
                end do

            case default
                call fatal_error('compute_nabla_matrices: unknown fit_point_type.')
            end select
        case(ommp_df_qm_top)
            call fatal_error("Unsupported case: grid points generated from QM topology.")
        case default
            call fatal_error('compute_nabla_matrices: unknown fit_top_type.')
        end select

        df%nabla_done = .true.

    end subroutine compute_nabla_matrices

    subroutine df_generate_grid(df)
        !! Generate charge and fit point coordinates from the stored
        !! topology pointers and generation parameters.
        !! This is called internally by df_init and can be called again
        !! if topologies have changed (e.g. after coordinate update).
        !! Delegates to generate_grid_from_topo for each grid.

        implicit none

        type(ommp_density_fit_type), intent(inout) :: df

        !! Delegate to generic grid generator
        call generate_grid_from_topo(df%charge_top, df%charge_coord, df%n_charges, &
                                     df%charge_n_pts_per_atom, df%charge_radius, &
                                     df%charge_point_type, 'charge')
        call generate_grid_from_topo(df%fit_top, df%fit_point_coord, df%n_pts, &
                                     df%fit_n_pts_per_atom, df%fit_radius, &
                                     df%fit_point_type, 'fit')

        !! Invalidate dependent computed quantities
        df%xinv_done = .false.
        df%x_done = .false.
        df%fit_done = .false.
        df%V_m2q_done = .false.
        df%V_p2q_done = .false.
        df%VXI_m_done = .false.
        df%VXI_p_done = .false.
        df%E_q2p_done = .false.
        df%E_pol_ene_done = .false.
        df%E_MM2Q_done = .false.
        df%dX_dr_done = .false.
        df%lambda_done = .false.

        !! Invalidate gradient matrices
        df%nabla_g_mm_is_null = .false.
        df%nabla_g_mm_is_identity = .false.
        df%nabla_g_qm_is_null = .false.
        df%nabla_g_qm_is_identity = .false.
        df%nabla_q_qm_is_null = .false.
        df%nabla_q_qm_is_identity = .false.
        df%nabla_q_mm_is_null = .false.
        df%nabla_q_mm_is_identity = .false.
        if(allocated(df%nabla_g_mm)) call mfree('[df_generate_grid] nabla_g_mm', df%nabla_g_mm)
        if(allocated(df%nabla_g_qm)) call mfree('[df_generate_grid] nabla_g_qm', df%nabla_g_qm)
        if(allocated(df%nabla_q_qm)) call mfree('[df_generate_grid] nabla_q_qm', df%nabla_q_qm)
        if(allocated(df%nabla_q_mm)) call mfree('[df_generate_grid] nabla_q_mm', df%nabla_q_mm)
    end subroutine df_generate_grid

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

        call df_generate_grid(df)

        call mallocate('df_init [target_charges]', df%n_charges, df%target_charges)
        df%target_charges = 0.0_rp

        call mallocate('df_init [fit_potential]', df%n_pts, df%fit_potential)
        df%fit_potential = 0.0_rp

        call mallocate('df_init [X]', df%n_pts, df%n_charges, df%X)
        call mallocate('df_init [Xinv]', df%n_charges, df%n_pts, df%Xinv)

        df%xinv_done = .false.
        df%x_done = .false.
        df%fit_done = .false.
        df%V_m2q_done = .false.
        df%V_p2q_done = .false.
        df%VXI_m_done = .false.
        df%VXI_p_done = .false.
        df%E_q2p_done = .false.
        df%E_pol_ene_done = .false.
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
        call mfree('df_terminate [E_MM2Q]', df%E_MM2Q)
        call mfree('df_terminate [dX_dr]', df%dX_dr)
        call mfree('df_terminate [lambda]', df%lambda)
        call mfree('df_terminate [nabla_g_mm]', df%nabla_g_mm)
        call mfree('df_terminate [nabla_g_qm]', df%nabla_g_qm)
        call mfree('df_terminate [nabla_q_qm]', df%nabla_q_qm)
        call mfree('df_terminate [nabla_q_mm]', df%nabla_q_mm)

    end subroutine df_terminate

    subroutine df_build_X(df)
        !! Build the design matrix X where X(i,j) = 1 / |r_charge_i - r_fit_j|.
        !! The matrix has dimensions (n_charges x n_pts) in Fortran layout.

        use mod_constants, only: eps_rp

        implicit none

        type(ommp_density_fit_type), intent(inout) :: df
        integer(ip) :: i, j
        real(rp) :: dr_x, dr_y, dr_z, dist2, dist

        if(.not. df%initialized) then
            call fatal_error("Density fit object not initialized!")
        end if

        !! TODO parallelize
        do i = 1, df%n_charges
            do j = 1, df%n_pts
                dr_x = df%charge_coord(1,i) - df%fit_point_coord(1,j)
                dr_y = df%charge_coord(2,i) - df%fit_point_coord(2,j)
                dr_z = df%charge_coord(3,i) - df%fit_point_coord(3,j)
                dist2 = dr_x*dr_x + dr_y*dr_y + dr_z*dr_z
                if(dist2 < eps_rp) then
                    call fatal_error('Charge and fitting point coincide!')
                end if
                dist = sqrt(dist2)
                df%X(j,i) = 1.0_rp / dist
            end do
        end do
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

        !! Lazy: compute Xinv if not already done
        if(.not. df%xinv_done) then
            call df_compute_Xinv_svd(df)
        end if

        df%fit_done = .false.
        df%V_p2q_done = .false.
        df%VXI_p_done = .false.
        df%E_q2p_done = .false.
        df%E_pol_ene_done = .false.

        !! Solve: fitted_charges = Xinv @ fit_potential
        call dgemv('N', df%n_charges, df%n_pts, 1.0_rp, df%Xinv, df%n_charges, df%fit_potential, 1, 0.0_rp, df%target_charges, 1)

        !! Forward: predicted_potential = X @ fitted_charges
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

        df%fit_done = .true.

        !! Report residual norm as percentage
        write(msg, '(A,F10.6,A)') 'Density fit residual norm = ', residual_norm * 100.0_rp, ' %'
        call ommp_message(trim(msg), 1, 'df')

        !! Report total fitted charge
        write(msg, '(A,F10.6)') 'Total fitted charge = ', sum(df%target_charges)
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

        if(.not. df%x_done) then
            call df_build_X(df)
        end if

        !! Edge case: nothing to do
        if(df%n_charges == 0 .or. df%n_pts == 0) call fatal_error("Either target or fit grids in density-fitting have no points.")

        min_dim = min(df%n_charges, df%n_pts)

        !! Allocate SVD workspace
        allocate(tmp_X(df%n_pts, df%n_charges))
        allocate(s(min_dim))
        allocate(u(df%n_pts, df%n_pts))
        allocate(vt(df%n_charges, df%n_charges))
        allocate(iwork(8 * min_dim))

        tmp_X = df%X
        !! Query optimal workspace size
        lwork = -1
        allocate(work(1))
        call dgesdd('A', df%n_pts, df%n_charges, tmp_X, df%n_pts, s, u, df%n_pts, vt, df%n_charges, work, lwork, iwork, info)
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))

        !! Compute full SVD (U is computed but only VT is needed for the pseudoinverse)
        call dgesdd('A', df%n_pts, df%n_charges, tmp_X, df%n_pts, s, u, df%n_pts, vt, df%n_charges, work, lwork, iwork, info)
        if(info /= 0) then
            call fatal_error('dgesdd SVD failed')
        end if
        
        
        !! Compute the pseudoinverse matrix
        !! 1. Compute (S_inv @ U^T)^T inplace starting from U
        do i = 1, min_dim
            u(:,i) = u(:,i) / s(i)
        end do
        !! Actually the desired matrix is in u(:,:n_charges)
        
        !! Now compute the pseudoinverse as Vt^T @ (U @ S-1) ^ T
        call dgemm('T', 'T', df%n_charges, df%n_pts, min_dim, 1.0_rp, vt, df%n_charges, u, df%n_pts, &
                   0.0_rp, df%Xinv, df%n_charges)

        !! Cleanup
        deallocate(s, u, vt, work, iwork, tmp_X)
        df%xinv_done = .true.
        call ommp_message('SVD-based Xinv computed', 2, 'df')
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

        if(.not. df%V_m2q_done) then
            if(.not. allocated(df%V_m2q)) then
                call mallocate('df_electrostatic_static [V_mm2df]', &
                               df%n_charges, df%V_m2q)
            end if

            df%V_m2q = 0.0_rp
            call potential_M2E(eel, df%charge_coord, df%V_m2q)
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

        if(.not. df%V_p2q_done .and. (eel%ipd_done .or. eel%pol_atoms == 0)) then
            if(.not. allocated(df%V_p2q)) then
                call mallocate('df_electrostatic_dipoles [V_pd2df]', &
                               df%n_charges, eel%n_ipd, df%V_p2q)
            end if

            df%V_p2q = 0.0_rp
            if(eel%pol_atoms > 0) then
                if(eel%amoeba) then
                    call potential_D2E(eel, df%charge_coord, df%V_p2q(:,_amoeba_D_))
                    call potential_D2E(eel, df%charge_coord, df%V_p2q(:,_amoeba_P_), .true.)
                else
                    call potential_D2E(eel, df%charge_coord, df%V_p2q(:,1))
                end if
            end if
            df%V_p2q_done = .true.
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

        if(.not. df%VXI_m_done) then
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

        if(.not. df%E_q2p_done) then
            if(.not. allocated(df%E_q2p)) then
                call mallocate('df_e_field_to_pol [E_q2p]', &
                               3_ip, eel%pol_atoms, df%E_q2p)
            end if

            df%E_q2p = 0.0_rp
            do i = 1, df%n_charges
                do j = 1, eel%pol_atoms
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

        if(.not. df%E_pol_ene_done) then
            !! Ensure E_q2p is computed first
            if(.not. df%E_q2p_done) then
                call df_e_field_to_pol(df, eel)
            end if

            ! TODO improve this double loop.
            eMM = 0.0_rp
            if(eel%amoeba) then
                !! Use only _amoeba_D_ dipoles, contracted with E_q2p
                do i = 1, 3
                    do j = 1, eel%pol_atoms
                        eMM = eMM - eel%ipd(i,j,_amoeba_D_) * df%E_q2p(i,j)
                    end do
                end do
            else
                do i = 1, 3
                    do j = 1, eel%pol_atoms
                        eMM = eMM - eel%ipd(i,j,1) * df%E_q2p(i,j)
                    end do
                end do
            end if

            df%E_pol_ene = 0.5_rp * eMM
            df%E_pol_ene_done = .true.
        end if
    end subroutine df_e_field_pol_ene

    subroutine df_electrostatics_for_geomgrad(df, eel)
        !! Compute electric field quantities needed for geometry gradients.
        !!
        !! 1. E_MM2Q(n_charges, 3, n_mm_atoms): electric field generated by
        !!    MM static multipoles (charges, dipoles, quadrupoles) at the
        !!    charge fitting sites.  This follows the same pattern as
        !!    field_M2E in mod_electrostatics.F90, but instead of accumulating
        !!    into a 2D array (3, n_sites), we store the contribution from each
        !!    MM atom separately so that the gradient w.r.t. MM coordinates can
        !!    be assembled by matrix-vector multiplication.
        !!
        !! 2. dX_dr(n_charges, 3, n_pts): derivative of the design matrix X
        !!    w.r.t. fit-point coordinates.  For each charge site i and grid
        !!    point j, dX_dr(i,:,j) = kernel(2) * dr where dr =
        !!    r_grid_j - r_charge_i.  This stores dX/dr_grid which is needed
        !!    for gradient computations where we need to know how the
        !!    potential at grid points changes with nuclear coordinates.

        use mod_memory, only: mallocate
        use mod_electrostatics, only: q_elec_prop, mu_elec_prop, quad_elec_prop, &
                                      coulomb_kernel, ommp_electrostatics_type

        implicit none

        type(ommp_density_fit_type), intent(inout) :: df
        type(ommp_electrostatics_type), intent(in) :: eel

        integer(ip) :: i, j, n_cpt, n_mm, n_pts, n_charges
        real(rp) :: kernel(5), dr(3), tmpV, tmpE(3), tmpEgr(6), tmpHE(10)

        !! ---------------------------------------------------------------
        !! Part 1: E_MM2Q — MM multipoles → charge sites
        !! ---------------------------------------------------------------
        if(.not. df%E_MM2Q_done) then
            n_mm = eel%top%mm_atoms
            n_cpt = df%n_charges

            if(.not. allocated(df%E_MM2Q)) then
                call mallocate('[df_electrostatics_for_geomgrad] E_MM2Q', &
                               df%n_charges, 3_ip, n_mm, df%E_MM2Q)
            end if

            df%E_MM2Q = 0.0_rp

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

                        df%E_MM2Q(j,:,i) = df%E_MM2Q(j,:,i) + tmpE
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

                        df%E_MM2Q(j,:,i) = df%E_MM2Q(j,:,i) + tmpE
                    end do
                end do
            end if

            df%E_MM2Q_done = .true.
        end if

        !! ---------------------------------------------------------------
        !! Part 2: dX_dr — derivative of design matrix X w.r.t. grid coords
        !! ---------------------------------------------------------------
        if(.not. df%dX_dr_done) then
            n_charges = df%n_charges
            n_pts     = df%n_pts

            if(.not. allocated(df%dX_dr)) then
                call mallocate('[df_electrostatics_for_geomgrad] dX_dr', &
                               3_ip, n_charges, n_pts, df%dX_dr)
            end if

            df%dX_dr = 0.0_rp

            !$omp parallel do default(shared) schedule(dynamic) &
            !$omp private(i,j,dr,kernel)
            do i = 1, n_charges
                do j = 1, n_pts
                    dr = df%charge_coord(:,i) - df%fit_point_coord(:,j)
                    call coulomb_kernel(dr, 2, kernel)

                    df%dX_dr(:,i,j) = kernel(2) * dr
                end do
            end do

            df%dX_dr_done = .true.
        end if

        if(.not. df%V_m2q_done) call df_electrostatic_static(df, eel)
        if(.not. df%V_p2q_done) call df_electrostatic_dipoles(df, eel)
    end subroutine df_electrostatics_for_geomgrad

    subroutine df_geomgrad(df, qmg, mmg, doqm, domm)
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

        use mod_io, only: fatal_error, print_matrix

        implicit none

        type(ommp_density_fit_type), intent(inout) :: df
        real(rp), intent(inout) :: qmg(:,:)   ! (3, n_qm_atoms)
        real(rp), intent(inout) :: mmg(:,:)   ! (3, n_mm_atoms)
        logical, intent(in) :: doqm
        logical, intent(in) :: domm

        integer(ip) :: a, b, i, j, k, l, s, n_pts
        real(rp) :: ef_grid(3,df%n_pts)

        if(.not. df%initialized) then
            call fatal_error("df_geomgrad: density fit object not initialized.")
        end if

        if(doqm) then
            qmg= 0.0_rp
        end if
        if(domm) then
            mmg = 0.0_rp
        end if

        !! Compute electric field at grid points from fitted charges:
        !! ef_grid(a,j) = sum_i target_charges(i) * dX_dr(i,a,j)
        ef_grid = 0.0_rp
        n_pts = df%n_pts
        do a = 1, 3
            do j = 1, n_pts
                do l = 1, df%n_charges
                    ef_grid(a,j) = ef_grid(a,j) + df%target_charges(l) * df%dX_dr(a,l,j)
                end do
            end do
        end do

        !! Ensure lambda is computed
        call df_compute_lambda(df)

        ! There are three components:
        ! 1. [ d/dr V(MM+POL) at q] @ q

        do i=1, df%mm_top%mm_atoms
            do a=1, 3
                do l=1, df%n_charges
                    mmg(a,i) = mmg(a,i) + df%target_charges(l) * df%E_MM2Q(l,a,i)
                end do
            end do
        end do

        if(.not. df%nabla_q_mm_is_null) then
            ! It should be something really similar to the derivatives on QM atoms
            call fatal_error("Dependency of fit charges from MM atoms' coordinates is not implemented.")
        end if

        if(df%nabla_q_qm_is_identity) then
            do i=1, df%mm_top%mm_atoms
                do l=1, df%n_charges ! Same of QM atoms
                    do a=1, 3
                        qmg(a,l) = qmg(a,l) - df%target_charges(l) * df%E_MM2Q(l,a,i)
                    end do
                end do
            end do    
        else if(df%nabla_q_qm_is_null) then
            ! nothing to do
            continue
        else
            do j=1, df%qm_top%mm_atoms
                do b=1, 3
                    do i=1, df%mm_top%mm_atoms
                            do l=1, df%n_charges ! Same of QM atoms
                                do a=1, 3
                                    qmg(b,j) = qmg(b,j) - &
                                    df%target_charges(l) * df%E_MM2Q(l,a,i) * df%nabla_q_qm(3*(l-1)+a, 3*(j-1)+b)
                                end do
                        end do
                    end do
                end do
            end do    
        end if

        ! 2. [ d/dr A] . lambda . q
        !! lambda is already computed and stored in df%lambda
        
        if(df%nabla_g_mm_is_identity) then
            do i=1, df%mm_top%mm_atoms
                do a=1, 3
                    do l=1, df%n_charges
                        mmg(a,i) = mmg(a,i) - df%dX_dr(a,l,i) * df%lambda(i) * df%target_charges(l)
                    end do
                end do
            end do

        else if(df%nabla_g_mm_is_null) then
            continue
        else
            do i=1, df%mm_top%mm_atoms
                do a=1, 3
                    do s=1, df%n_pts
                        do l=1, df%n_charges
                            do b=1, 3
                                mmg(a,i) = mmg(a,i) + &
                                df%dX_dr(b,l,s) * df%nabla_g_mm(b+(s-1)*3, a+(i-1)*3) * df%lambda(s) * df%target_charges(l)
                            end do
                        end do
                    end do
                end do
            end do
        end if

        if(df%nabla_q_qm_is_identity) then
            do i=1, df%qm_top%mm_atoms
                do a=1, 3
                    do s=1, df%n_pts
                        qmg(a,i) = qmg(a,i) + df%dX_dr(a,i,s) * df%lambda(s) * df%target_charges(i)
                    end do
                end do
            end do

        else if(df%nabla_g_mm_is_null) then
            continue
        else
            do i=1, df%qm_top%mm_atoms
                do a=1, 3
                    do s=1, df%n_pts
                        do l=1, df%n_charges
                            do b=1, 3
                                qmg(a,i) = qmg(a,i) + &
                                df%dX_dr(b,l,s) * df%nabla_q_qm(b+(l-1)*3, a+(i-1)*3) * df%lambda(s) * df%target_charges(l)
                            end do
                        end do
                    end do
                end do
            end do
        end if


        if(.not. df%nabla_q_mm_is_null) then
            ! It should be something really similar to the derivatives on QM atoms
            call fatal_error("Dependency of fit charges from MM atoms' coordinates is not implemented.")
        end if

        if(df%nabla_g_qm_is_identity) then
            call fatal_error("Dependency of grid points from QM atoms' coordinates is not implemented.")
            ! Wired dependency from QM atoms positions and grid points, this should not happen for now
            do i=1, df%qm_top%mm_atoms
                do a=1, 3
                    do l=1, df%n_charges
                        qmg(a,i) = qmg(a,i) + df%dX_dr(a,l,i) * df%lambda(i) * df%target_charges(l)
                    end do
                end do
            end do

        else if(df%nabla_g_qm_is_null) then
            continue
        else
            do i=1, df%mm_top%mm_atoms
                do a=1, 3
                    do s=1, df%n_pts
                        do l=1, df%n_charges
                            do b=1, 3
                                qmg(a,i) = qmg(a,i) + &
                                df%dX_dr(b,l,s) * df%nabla_g_qm(b+(s-1)*3, a+(i-1)*3) * df%lambda(s) * df%target_charges(l)
                            end do
                        end do
                    end do
                end do
            end do
        end if
        

    end subroutine df_geomgrad

    subroutine df_compute_lambda(df)
        !! Compute the Lagrange multiplier vector.
        !!
        !!   lambda(s) = sum_l Xinv(l,s) * V_m2q(l)  (= VXI_m)
        !!
        !! This is a lazy-evaluated operation: it triggers computation of
        !! Xinv if not already available, then performs a single dgemv
        !! to compute lambda = Xinv^T @ V_m2q.
        !!
        !! NOTE: V_m2q must already be populated (e.g. via
        !! df_electrostatic_static or df_project_static) before calling
        !! this function.  The higher-level ommp_df_compute_lambda in
        !! mod_interface handles that dependency.

        use mod_memory, only: mallocate

        implicit none

        type(ommp_density_fit_type), intent(inout) :: df
        integer(ip) :: i, j

        if(.not. df%xinv_done) then
            call df_compute_Xinv_svd(df)
        end if

        if(.not. df%V_m2q_done .or. .not. df%V_p2q_done) then
            call fatal_error("Call df_electrostatics_for_geomgrad before df_compute_lambda")
        end if

        if(.not. df%lambda_done) then

            if(.not. allocated(df%lambda)) then
                call mallocate('[df_compute_lambda] lambda', df%n_pts, df%lambda)
            end if
            df%lambda = 0.0_rp

            do i=1, df%n_pts
                do j=1, df%n_charges
                    df%lambda(i) = df%lambda(i) + df%V_m2q(j) * df%Xinv(j,i)
                end do
            end do

            df%lambda_done = .true.
        end if
    end subroutine df_compute_lambda

end module
