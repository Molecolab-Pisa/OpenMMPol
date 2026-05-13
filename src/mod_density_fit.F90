module mod_density_fit
!! This module implements density fitting of QM charge distributions
!! onto a set of fitting points (typically MM atom positions).

    use mod_memory, only: ip, rp, lp
    use mod_constants, only: ommp_df_solver_svd, ommp_df_svd_rcond_default
    use mod_io, only: fatal_error, ommp_message

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

        real(rp), allocatable :: V_p2q(:)
        !! Electrostatic potential from MM induced dipoles at charge coordinates

        logical(lp) :: V_m2q_done = .false.
        !! Flag indicating whether V_mm2df has been computed

        logical(lp) :: V_p2q_done = .false.
        !! Flag indicating whether V_pd2df has been computed

        real(rp), allocatable :: VXI_m(:)
        !! Projected static quantity: V_m2q @ Xinv, for Fock matrix

        real(rp), allocatable :: VXI_p(:)
        !! Projected dipole quantity: V_p2q @ Xinv, for Fock matrix

        logical(lp) :: VXI_m_done = .false.
        !! Flag indicating whether VXI_m has been computed

        logical(lp) :: VXI_p_done = .false.
        !! Flag indicating whether VXI_p has been computed

    end type ommp_density_fit_type

    public :: ommp_density_fit_type
    public :: df_init, df_terminate
    public :: df_solve
    public :: df_electrostatic_static, df_electrostatic_dipoles
    public :: df_project_static, df_project_dipoles

contains

    subroutine df_init(df, charge_coord, fit_point_coord)
        !! Initialize the density fit object

        use mod_memory, only: mallocate

        implicit none

        type(ommp_density_fit_type), intent(inout) :: df
        real(rp), intent(in) :: charge_coord(:,:)
        !! Coordinates of the charge positions (3 x n_charges)
        real(rp), intent(in) :: fit_point_coord(:,:)
        !! Coordinates of the fitting points (3 x n_pts)

        integer(ip) :: n_pts
        integer(ip) :: n_charges

        if(df%initialized) then
            call fatal_error("Density fit object already initialized!")
        end if

        n_charges = size(charge_coord, 2)
        n_pts = size(fit_point_coord, 2)

        df%n_pts = n_pts
        df%n_charges = n_charges

        call mallocate('df_init [charge_coord]', 3_ip, n_charges, df%charge_coord)
        df%charge_coord = charge_coord

        call mallocate('df_init [fit_point_coord]', 3_ip, n_pts, df%fit_point_coord)
        df%fit_point_coord = fit_point_coord

        call mallocate('df_init [target_charges]', n_charges, df%target_charges)
        df%target_charges = 0.0_rp

        call mallocate('df_init [fit_potential]', n_pts, df%fit_potential)
        df%fit_potential = 0.0_rp

        call mallocate('df_init [X]', n_pts, n_charges, df%X)
        call mallocate('df_init [Xinv]', n_charges, n_pts, df%Xinv)

        df%xinv_done = .false.
        df%x_done = .false.
        df%fit_done = .false.
        df%V_m2q_done = .false.
        df%V_p2q_done = .false.
        df%VXI_m_done = .false.
        df%VXI_p_done = .false.

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
                               df%n_charges, df%V_p2q)
            end if

            df%V_p2q = 0.0_rp
            call potential_D2E(eel, df%charge_coord, df%V_p2q)
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
                call mallocate('df_project_dipoles [VXI_p]', df%n_pts, df%VXI_p)
            end if

            df%VXI_p = 0.0_rp
            call dgemv('T', df%n_charges, df%n_pts, 1.0_rp, &
                       df%Xinv, df%n_charges, &
                       df%V_p2q, 1, 0.0_rp, df%VXI_p, 1)
            df%VXI_p_done = .true.
        end if
    end subroutine df_project_dipoles

end module
