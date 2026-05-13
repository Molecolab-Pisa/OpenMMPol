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
        !! Design matrix, dimensions (n_charges x n_pts)

        real(rp), allocatable :: Xinv(:,:)
        !! Pseudoinverse design matrix, dimensions (n_charges x n_pts).
        !! Satisfies fitted_charges = Xinv @ fit_potential

        real(rp), allocatable :: fitted_charges(:)
        !! Fitted charges

        logical(lp) :: initialized = .false.
        !! Flag indicating whether the object is initialized

        logical(lp) :: xinv_done = .false.
        !! Flag indicating whether Xinv has been computed

        logical(lp) :: x_done = .false.
        !! Flag indicating whether X has been computed

        logical(lp) :: fit_done = .false.
        !! Flag indicating whether charges have been computed from fit_potential

    end type ommp_density_fit_type

    public :: ommp_density_fit_type
    public :: df_init, df_terminate
    public :: df_solve

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

        implicit none

        type(ommp_density_fit_type), intent(inout) :: df
        real(rp) :: alpha

        if(.not. df%initialized) then
            call fatal_error("Density fit object not initialized!")
        end if

        !! Lazy: build X if not already done
        if(.not. df%x_done) then
            call df_build_X(df)
            df%x_done = .true.
        end if

        !! Lazy: compute Xinv if not already done
        if(.not. df%xinv_done) then
            call df_compute_Xinv_svd(df)
            df%xinv_done = .true.
        end if

        call dgemv('N', df%n_charges, df%n_pts, 1.0_rp, df%Xinv, df%n_charges, df%fit_potential, 1, 0.0_rp, df%target_charges, 1)
        df%fit_done = .true.
    end subroutine df_solve

    subroutine df_compute_Xinv_svd(df)

        implicit none

        type(ommp_density_fit_type), intent(inout) :: df
        real(rp) :: svd_rcond, max_sv
        integer(ip) :: min_dim, lwork, info, i, j, k
        integer(ip), dimension(:), allocatable :: iwork
        real(rp), dimension(:), allocatable :: s, work
        real(rp), dimension(:,:), allocatable :: u, vt, tmp_X

        if(.not. df%initialized) then
            call fatal_error("Density fit object not initialized!")
        end if

        if(.not. df%x_done) then
            call df_build_X(df)
        end if

        svd_rcond = df%svd_rcond

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

        call ommp_message('SVD-based Xinv computed', 2, 'df')
    end subroutine df_compute_Xinv_svd

end module
