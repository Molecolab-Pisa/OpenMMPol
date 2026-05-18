#include "f_cart_components.h"
module mod_density_fit
!! This module implements density fitting of QM charge distributions
!! onto a set of fitting points (typically MM atom positions).

    use mod_memory, only: ip, rp, lp, mallocate, mfree
    use mod_constants, only: ommp_df_solver_svd, &
                             ommp_df_svd_rcond_default, &
                             ommp_df_charge_qm_atoms, &
                             ommp_df_charge_fibonacci, &
                             ommp_df_charge_cubic, &
                             ommp_df_fit_mm_atoms, &
                             ommp_df_fit_cubic
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
        !! Type of charge points: OMMP_DF_CHARGE_QM_ATOMS = 1, etc.

        integer(ip) :: charge_n_pts_per_atom = 0
        !! Number of charge points per source atom (for grid-based strategies)

        real(rp) :: charge_radius = 0.0_rp
        !! Radius parameter for charge point generation

        integer(ip) :: fit_point_type = 0
        !! Type of fit points: OMMP_DF_FIT_MM_ATOMS = 1, etc.

        integer(ip) :: fit_n_pts_per_atom = 0
        !! Number of fit points per source atom (for grid-based strategies)

        real(rp) :: fit_radius = 0.0_rp
        !! Radius parameter for fit point generation

        !! Source topologies for grid generation
        type(ommp_topology_type), pointer :: qm_top
        !! Pointer to QM topology (source of charge points)

        type(ommp_topology_type), pointer :: mm_top
        !! Pointer to MM topology (source of fit points)

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

        real(rp) :: E_pol_ene
        !! Polarization energy from fitted-charge electric field:
        !! E = -0.5 * sum_j ipd(:,j) .dot. E_q2p(:,j)

        logical(lp) :: E_pol_ene_done = .false.
        !! Flag indicating whether E_pol_ene has been computed

    end type ommp_density_fit_type

    public :: ommp_density_fit_type
    public :: OMMP_DF_CHARGE_QM_ATOMS, OMMP_DF_CHARGE_FIBONACCI, OMMP_DF_CHARGE_CUBIC, OMMP_DF_FIT_MM_ATOMS, OMMP_DF_FIT_CUBIC
    public :: df_init, df_terminate
    public :: df_solve, df_generate_grid
    public :: df_electrostatic_static, df_electrostatic_dipoles
    public :: df_project_static, df_project_dipoles
    public :: df_e_field_to_pol, df_e_field_pol_ene

contains

    subroutine df_generate_grid(df)
        !! Generate charge and fit point coordinates from the stored
        !! topology pointers and generation parameters.
        !! This is called internally by df_init and can be called again
        !! if topologies have changed (e.g. after coordinate update).

        implicit none

        type(ommp_density_fit_type), intent(inout) :: df
        integer(ip) :: n_qm, n_mm, i, j, idx, ii
        real(rp) :: r(3), theta, phi, y, x0, y0, z0
        character(len=256) :: msg

        real(rp), parameter :: golden_angle = 2.0_rp * atan(1.0_rp) * (3.0_rp - sqrt(5.0_rp))

        !! Free previously generated coordinates
        if(allocated(df%charge_coord)) call mfree('[df_generate_grid] charge_coord', df%charge_coord)
        if(allocated(df%fit_point_coord)) call mfree('[df_generate_grid] fit_point_coord', df%fit_point_coord)

        !! Generate charge points
        select case(df%charge_point_type)
        case(OMMP_DF_CHARGE_QM_ATOMS)
            
            if(.not. associated(df%qm_top)) then
                call fatal_error('df_generate_grid: OMMP_DF_CHARGE_QM_ATOMS selected but qm_top is not associated.')
            end if
            n_qm = df%qm_top%mm_atoms
            if(n_qm == 0) then
                call fatal_error('df_generate_grid: QM topology has no atoms.')
            end if
            call mallocate('[df_generate_grid] charge_coord', 3_ip, n_qm, df%charge_coord)
            df%charge_coord = df%qm_top%cmm(:,1:n_qm)
            df%n_charges = n_qm
            write(msg, '(A,I0,A)') 'Generated ', n_qm, ' charge points from QM atoms'
            call ommp_message(trim(msg), 1, 'df')

        case(OMMP_DF_CHARGE_FIBONACCI)
            if(.not. associated(df%qm_top)) then
                call fatal_error('df_generate_grid: OMMP_DF_CHARGE_FIBONACCI selected but qm_top is not associated.')
            end if
            n_qm = df%qm_top%mm_atoms
            if(n_qm == 0) then
                call fatal_error('df_generate_grid: QM topology has no atoms.')
            end if
            if(df%charge_n_pts_per_atom <= 0) then
                call fatal_error('df_generate_grid: OMMP_DF_CHARGE_FIBONACCI requires n_pts_per_atom > 0.')
            end if
            if(df%charge_radius <= 0.0_rp) then
                call fatal_error('df_generate_grid: OMMP_DF_CHARGE_FIBONACCI requires radius > 0.')
            end if
            call mallocate('[df_generate_grid] charge_coord', 3_ip, &
                           n_qm * df%charge_n_pts_per_atom, df%charge_coord)
            do i = 1, n_qm
                x0 = df%qm_top%cmm(1, i)
                y0 = df%qm_top%cmm(2, i)
                z0 = df%qm_top%cmm(3, i)
                do j = 1, df%charge_n_pts_per_atom
                    y = (2.0_rp * real(j, rp) - 1.0_rp) / real(df%charge_n_pts_per_atom, rp) - 1.0_rp
                    theta = real(golden_angle, rp) * real(j, rp)
                    phi = sqrt(1.0_rp - y**2)
                    idx = (i - 1) * df%charge_n_pts_per_atom + j
                    df%charge_coord(1, idx) = x0 + cos(theta) * phi * df%charge_radius
                    df%charge_coord(2, idx) = y0 + y * df%charge_radius
                    df%charge_coord(3, idx) = z0 + sin(theta) * phi * df%charge_radius
                end do
            end do
            df%n_charges = n_qm * df%charge_n_pts_per_atom
            write(msg, '(A,I0,A,I0,A)') 'Generated ', df%charge_n_pts_per_atom, ' fibonacci points per atom (', &
                   n_qm, ' atoms)'
            call ommp_message(trim(msg), 1, 'df')

        case(OMMP_DF_CHARGE_CUBIC)
            if(.not. associated(df%qm_top)) then
                call fatal_error('df_generate_grid: OMMP_DF_CHARGE_CUBIC selected but qm_top is not associated.')
            end if
            n_qm = df%qm_top%mm_atoms
            if(n_qm == 0) then
                call fatal_error('df_generate_grid: QM topology has no atoms.')
            end if
            if(df%charge_radius <= 0.0_rp) then
                call fatal_error('df_generate_grid: OMMP_DF_CHARGE_CUBIC requires radius > 0.')
            end if
            !! 7 points per atom: 1 center + 6 face centers of a cube
            call mallocate('[df_generate_grid] charge_coord', 3_ip, n_qm * 7_ip, df%charge_coord)
            do ii = 1, n_qm
                i = (ii-1) * 7 + 1
                x0 = df%qm_top%cmm(1, ii)
                y0 = df%qm_top%cmm(2, ii)
                z0 = df%qm_top%cmm(3, ii)

                df%charge_coord(1, i:i+6) = x0
                df%charge_coord(2, i:i+6) = y0
                df%charge_coord(3, i:i+6) = z0
                df%charge_coord(1, i+1) = x0 + df%charge_radius
                df%charge_coord(1, i+2) = x0 - df%charge_radius
                df%charge_coord(2, i+3) = y0 + df%charge_radius
                df%charge_coord(2, i+4) = y0 - df%charge_radius
                df%charge_coord(3, i+5) = z0 + df%charge_radius
                df%charge_coord(3, i+6) = z0 - df%charge_radius

            end do
            df%n_charges = n_qm * 7_ip
            write(msg, '(A,I0,A)') 'Generated 7 cubic points per atom (', &
                   n_qm, ' atoms)'
            call ommp_message(trim(msg), 1, 'df')
        case default
            call fatal_error('df_generate_grid: unknown charge_point_type.')
        end select

        !! Generate fit points
        select case(df%fit_point_type)
        case(OMMP_DF_FIT_MM_ATOMS)
            if(.not. associated(df%mm_top)) then
                call fatal_error('df_generate_grid: OMMP_DF_FIT_MM_ATOMS selected but mm_top is not associated.')
            end if
            n_mm = df%mm_top%mm_atoms
            if(n_mm == 0) then
                call fatal_error('df_generate_grid: MM topology has no atoms.')
            end if
            call mallocate('[df_generate_grid] fit_point_coord', 3_ip, n_mm, df%fit_point_coord)
            df%fit_point_coord = df%mm_top%cmm(:,1:n_mm)
            df%n_pts = n_mm
            write(msg, '(A,I0,A)') 'Generated ', n_mm, ' fit points from MM atoms'
            call ommp_message(trim(msg), 1, 'df')
            
        case(OMMP_DF_FIT_CUBIC)
            if(.not. associated(df%mm_top)) then
                call fatal_error('df_generate_grid: OMMP_DF_FIT_CUBIC selected but mm_top is not associated.')
            end if
            n_mm = df%mm_top%mm_atoms
            if(n_mm == 0) then
                call fatal_error('df_generate_grid: MM topology has no atoms.')
            end if
            if(df%fit_radius <= 0.0_rp) then
                call fatal_error('df_generate_grid: OMMP_DF_FIT_CUBIC requires radius > 0.')
            end if
            !! 7 points per atom: 1 center + 6 face centers of a cube
            call mallocate('[df_generate_grid] fit_point_coord', 3_ip, n_mm * 7_ip, df%fit_point_coord)
            do ii = 1, n_mm
                x0 = df%mm_top%cmm(1, ii)
                y0 = df%mm_top%cmm(2, ii)
                z0 = df%mm_top%cmm(3, ii)
                i = (ii-1) * 7 + 1
                df%fit_point_coord(1, i:i+6) = x0
                df%fit_point_coord(2, i:i+6) = y0
                df%fit_point_coord(3, i:i+6) = z0
                df%fit_point_coord(1, i+1) = x0 + df%fit_radius
                df%fit_point_coord(1, i+2) = x0 - df%fit_radius
                df%fit_point_coord(2, i+3) = y0 + df%fit_radius
                df%fit_point_coord(2, i+4) = y0 - df%fit_radius
                df%fit_point_coord(3, i+5) = z0 + df%fit_radius
                df%fit_point_coord(3, i+6) = z0 - df%fit_radius
            end do
            df%n_pts = n_mm * 7_ip
            write(msg, '(A,I0,A)') 'Generated 7 cubic fit points per atom (', &
                   n_mm, ' atoms)'
            call ommp_message(trim(msg), 1, 'df')

        case default
            call fatal_error('df_generate_grid: unknown fit_point_type.')
        end select

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
    end subroutine df_generate_grid

    subroutine df_init(df, qm_top, mm_top, &
                       charge_point_type, charge_n_pts_per_atom, charge_radius, &
                       fit_point_type, fit_n_pts_per_atom, fit_radius)
        !! Initialize the density fit object.
        !! Coordinates are generated from the provided topology pointers
        !! according to the point generation strategy.

        use mod_memory, only: mallocate

        implicit none

        type(ommp_density_fit_type), intent(inout) :: df
        type(ommp_topology_type), intent(in), target :: qm_top
        !! QM topology providing charge point coordinates
        type(ommp_topology_type), intent(in), target :: mm_top
        !! MM topology providing fit point coordinates
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

        integer(ip) :: i, j, idx
        real(rp) :: x0, y0, z0, theta, phi, y

        if(df%initialized) then
            call fatal_error("Density fit object already initialized!")
        end if

        if(qm_top%mm_atoms == 0) then
            call fatal_error('df_init: QM topology has no atoms.')
        end if
        if(mm_top%mm_atoms == 0) then
            call fatal_error('df_init: MM topology has no atoms.')
        end if            

        !! Store point generation configuration
        df%charge_point_type = charge_point_type
        df%charge_n_pts_per_atom = charge_n_pts_per_atom
        df%charge_radius = charge_radius
        df%fit_point_type = fit_point_type
        df%fit_n_pts_per_atom = fit_n_pts_per_atom
        df%fit_radius = fit_radius
        df%qm_top => qm_top
        df%mm_top => mm_top

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

        if(.not. df%V_p2q_done .and. eel%ipd_done) then
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

end module
