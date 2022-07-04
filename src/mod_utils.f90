module mod_utils
    use mod_memory, only: ip

    implicit none
    private

    public :: sort_ivec, sort_ivec_inplace

    contains

    subroutine sort_ivec(iv, ov)
        !! This is a simple -- and unefficient -- routine to sort a vector of int
        !! it is just used during some output to simplify comparisons with older 
        !! version of the code
        use mod_memory, only: mallocate

        implicit none

        !! Input vector
        integer(ip), dimension(:), intent(in) :: iv
        !! Output, sorted vector
        integer(ip), allocatable, intent(out) :: ov(:)
        integer(ip) :: i, imin(1)
        logical, allocatable :: mask(:)

        call mallocate('sort_ivec', size(iv), ov)
        allocate(mask(size(iv)))
        mask = .true.

        do i = 1, size(iv)
            imin = minloc(iv, mask)
            ov(i) = iv(imin(1))
            mask(imin(1)) = .false.
        end do

        deallocate(mask)

    end subroutine sort_ivec
    
    subroutine sort_ivec_inplace(iv)
        !! Inplace equivalent of [[sort_ivec]] routine.

        implicit none

        !! Vector to be reordered in place
        integer(ip), dimension(:), intent(inout) :: iv
        integer(ip), allocatable :: ov(:)

        call sort_ivec(iv, ov)
        iv = ov

    end subroutine sort_ivec_inplace

end module mod_utils
