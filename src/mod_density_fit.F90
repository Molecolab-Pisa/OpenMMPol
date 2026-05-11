module mod_density_fit
!! This module implements ...

    use mod_memory, only: ip, rp, lp

    implicit none
    private

    type ommp_density_fit_type
        integer(ip) :: n_pts = 0
        !! Number of fitting points

        integer(ip) :: n_charges = 0
        !! Number of target charges

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
        !! Inverse (pseudoinverse) design matrix, dimensions (n_pts x n_charges)

        logical(lp) :: initialized = .false.
        !! Flag indicating whether the object is initialized

    end type ommp_density_fit_type

    public :: ommp_density_fit_type
    public :: df_init, df_terminate

contains

    subroutine df_init(df, charge_coord, fit_point_coord)
        !! Initialize the density fit object

        use mod_memory, only: mallocate
        use mod_io, only: fatal_error

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

        call mallocate('df_init [X]', n_charges, n_pts, df%X)
        df%X = 0.0_rp

        call mallocate('df_init [Xinv]', n_pts, n_charges, df%Xinv)
        df%Xinv = 0.0_rp

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

        df%n_pts = 0
        df%n_charges = 0
        df%initialized = .false.
    end subroutine df_terminate

end module
