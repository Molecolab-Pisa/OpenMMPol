module mod_utils
    use mod_memory, only: ip

    implicit none
    private

    public :: sort_ivec, sort_ivec_inplace, skip_lines
    public :: starts_with_alpha, isreal, isint, tokenize

    contains
   
    function starts_with_alpha(s)
        ! Decide if a string starts with a letter
        implicit none

        character(len=*), intent(in) :: s
        logical :: starts_with_alpha

        starts_with_alpha = (scan(s(1:1), &
        'qwertyuiopasdfghjklzxcvbnmQWERTYUIOPASDFGHJKLZXCVBNM') /= 0)
        return
    end function

    function isint(s)
        ! Decide if a string contains an integer
        implicit none

        character(len=*), intent(in) :: s
        logical :: isint
        
        isint = (verify(s, '+-1234567890') == 0)
        return
    end function

    function isreal(s)
        ! Decide if a string contains a real
        implicit none

        character(len=*), intent(in) :: s
        logical :: isreal

        isreal = (verify(s, '+-1234567890.') == 0)
        isreal = isreal .and. (scan(s, '.') /= 0)
        return
    end function

    function tokenize(s, ib, ntok)
        ! This function, given a string returns the first printable character
        ! if ib is absent or the last printable character after s(ib) if ib is 
        ! present. It's used to subdivide a string in tokens.

        use mod_memory, only: ip
        implicit none

        character(len=120), intent(in) :: s
        integer(ip), intent(inout), optional :: ib
        integer(ip), intent(in), optional :: ntok
        integer(ip) :: tokenize

        integer(ip) :: i, slen, ich, itok

        slen = len(s)
        ! Default return for end of string
        tokenize = -1
        if(present(ib)) then

            ! This is a very unreasonable case
            if(ib > slen) return
        
            do i=ib, slen
                ! Search the first valid char and save it in ib
                ich = iachar(s(i:i))
                if(ich > 32 .and. ich /= 127) exit
            end do
            ib = i
            if(ib >= slen) return
            
            if(present(ntok)) then
                itok = ntok
            else
                itok = 1
            end if

            do i=ib+1, slen
                ich = iachar(s(i:i))
                if(ich <= 32 .or. ich == 127) then 
                    ich = iachar(s(i-1:i-1))
                    if(ich > 32 .and. ich /= 127) itok = itok - 1
                end if
                if(itok == 0) exit
            end do
            tokenize = i-1
        else 
            do i=1, slen
                ich = iachar(s(i:i))
                if(ich > 32 .and. ich /= 127) exit
            end do
            tokenize = i
        end if
    end function
    
    subroutine skip_lines(f, n)
        !! This subroutine just skip lines while reading an input
        !! file
        use mod_memory, only: ip
        
        implicit none

        !! unit file
        integer(ip), intent(in) :: f
        !! number of line to be skipped
        integer(ip), intent(in) :: n

        integer(ip) :: i
        ! character(len=512) :: s

        do i=1, n
            read(f, *)
            ! read(f, *) s
            ! write(6, *) "SKIPPING ", i, "/", n, "  '", trim(s),"'"
        end do

    end subroutine skip_lines


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
