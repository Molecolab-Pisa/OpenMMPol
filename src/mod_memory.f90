module mod_memory
    !! This module is used to handle the memory, the variable
    !! kinds, the dynamic allocation and the optional soft
    !! memory limit of the openMMPol library.

    use iso_c_binding
    use mod_constants, only: ip, rp, lp, OMMP_STR_CHAR_MAX
    use mod_io, only: fatal_error, ommp_message

    implicit none
    private 

    integer(ip), parameter :: iof_memory = 6
    !! Unit file for memory errors, warning and debug
    integer(ip) :: maxmem !! Max memory that can be allocated in bytes
    integer(ip) :: usedmem !! Memory that is currently used by the code
    integer(ip) :: size_of_int !! Number of bytes for an integer
    integer(ip) :: size_of_real !! Number of bytes for a real
    integer(ip) :: size_of_logical !! Number of bytes for a logical (?)
    logical(lp) :: is_init = .false.
    logical :: do_chk_limit !! Decide if the soft memory limit is on

    public :: rp, ip, lp
    public :: mallocate, mfree, memory_init
    public :: use_8bytes_int 
    
    interface mallocate
        !! Interface to perform memory allocation within the
        !! openMMPol library, it can be called for 1,2 and 
        !! 3-dimensional arrays of either integer or real
        module procedure r_alloc1
        module procedure r_alloc2
        module procedure r_alloc3
        module procedure i_alloc1
        module procedure i_alloc2
        module procedure i_alloc3
        module procedure l_alloc1
        module procedure l_alloc2
    end interface mallocate

    interface mfree
        !! Interface to perform memory deallocation within the
        !! openMMPol library, it can be called for 1,2 and 
        !! 3-dimensional arrays of either integer or real
        module procedure r_free1
        module procedure r_free2
        module procedure r_free3
        module procedure i_free1
        module procedure i_free2
        module procedure i_free3
        module procedure l_free1
        module procedure l_free2
    end interface mfree

    contains

    function use_8bytes_int() bind(c, name='__use_8bytes_int')
        !! This function is used to know if the library is
        !! compiled using integer 8 or 4 bytes long.
        logical(kind=c_bool) :: use_8bytes_int

#ifdef USE_I8
        use_8bytes_int = .true.
#else
        use_8bytes_int = .false.
#endif
    end function use_8bytes_int

    subroutine memory_init(do_chk, max_bytes)
        !! Routine used to initialize the memory module. It should
        !! be called during the module initialization.
        implicit none

        logical :: do_chk !! Switch for memory soft limit 
        integer(ip), intent(in) :: max_bytes !! Amount of memory available in bytes
        integer(ip) :: my_int !! Integer used only as target for sizeof
        real(rp) :: my_real !! Real used only as target for sizeof
        logical(lp) :: my_bool
        intrinsic :: sizeof

        if(.not. is_init) then
            do_chk_limit = do_chk
            maxmem = max_bytes
            usedmem = 0
            size_of_real = sizeof(my_real)
            size_of_int = sizeof(my_int)
            size_of_logical = sizeof(my_bool)
            is_init = .true.
        end if
    end subroutine memory_init
    
    subroutine r_alloc1(string, len1, v)
        !! Allocate a 1-dimensional array of reals
        implicit none

        character(len=*), intent(in) :: string
        !! Human-readable description string of the allocation
        !! operation, just for output purpose.
        integer(ip), intent(in) :: len1
        !! Dimension of the vector
        real(rp), allocatable, intent(inout) :: v(:)
        !! Vector to allocate

        integer(ip) :: istat
 
        if(.not. is_init) call memory_init(.false., 0)
        allocate(v(len1), stat=istat)
        call chk_alloc(string, len1*size_of_real, istat)
    end subroutine r_alloc1
  
    subroutine r_alloc2(string, len1, len2, v)
        !! Allocate a 2-dimensional array of reals
        implicit none

        character(len=*), intent(in) :: string
        !! Human-readable description string of the allocation
        !! operation, just for output purpose.
        integer(ip), intent(in) :: len1, len2
        !! Dimensions of the vector
        real(rp), allocatable, intent(inout) :: v(:,:)
        !! Vector to allocate

        integer(ip) :: istat

        if(.not. is_init) call memory_init(.false., 0)
        allocate(v(len1, len2), stat=istat)
        call chk_alloc(string, len1*len2*size_of_real, istat)
    end subroutine r_alloc2
  
    subroutine r_alloc3(string, len1, len2, len3, v)
        !! Allocate a 3-dimensional array of reals
        implicit none

        character(len=*), intent(in) :: string
        !! Human-readable description string of the allocation
        !! operation, just for output purpose.
        integer(ip), intent(in) :: len1, len2, len3
        !! Dimensions of the vector
        real(rp), allocatable, intent(inout) :: v(:,:,:)
        !! Vector to allocate
        
        integer(ip) :: istat

        if(.not. is_init) call memory_init(.false., 0)
        allocate(v(len1, len2, len3), stat=istat)
        call chk_alloc(string, len1*len2*len3*size_of_real, istat)
    end subroutine r_alloc3

    subroutine i_alloc1(string, len1, v)
        !! Allocate a 1-dimensional array of integers
        implicit none

        character (len=*), intent(in) :: string
        !! Human-readable description string of the allocation
        !! operation, just for output purpose.
        integer(ip), intent(in) :: len1
        !! Dimension of the vector
        integer(ip), allocatable, intent(inout) :: v(:)
        !! Vector to allocate

        integer(ip) :: istat 

        if(.not. is_init) call memory_init(.false., 0)
        allocate(v(len1), stat=istat)
        call chk_alloc(string, len1*size_of_int, istat)
    end subroutine i_alloc1

    subroutine i_alloc2(string, len1, len2, v)
        !! Allocate a 2-dimensional array of integers
        implicit none

        character (len=*), intent(in) :: string
        !! Human-readable description string of the allocation
        !! operation, just for output purpose.
        integer(ip), intent(in) :: len1, len2
        !! Dimensions of the vector
        integer(ip), allocatable, intent(inout) :: v(:,:)
        !! Vector to allocate
 
        integer(ip) :: istat 

        if(.not. is_init) call memory_init(.false., 0)
        allocate(v(len1, len2), stat=istat)
        call chk_alloc(string, len1*len2*size_of_int, istat)
    end subroutine i_alloc2
 
    subroutine i_alloc3(string, len1, len2, len3, v)
        !! Allocate a 3-dimensional array of integers
        implicit none

        character (len=*), intent(in) :: string
        !! Human-readable description string of the allocation
        !! operation, just for output purpose.
        integer(ip), intent(in) :: len1, len2, len3
        !! Dimensions of the vector
        integer(ip), allocatable, intent(inout) :: v(:,:,:)
        !! Vector to allocate

        integer(ip) :: istat

        if(.not. is_init) call memory_init(.false., 0)
        allocate(v(len1, len2, len3), stat=istat)
        call chk_alloc(string, len1*len2*len3*size_of_int, istat)
    end subroutine i_alloc3
    
    subroutine l_alloc1(string, len1, v)
        !! Allocate a 1-dimensional array of reals
        implicit none

        character(len=*), intent(in) :: string
        !! Human-readable description string of the allocation
        !! operation, just for output purpose.
        integer(ip), intent(in) :: len1
        !! Dimension of the vector
        logical(lp), allocatable, intent(inout) :: v(:)
        !! Vector to allocate

        integer(ip) :: istat
 
        if(.not. is_init) call memory_init(.false., 0)
        allocate(v(len1), stat=istat)
        call chk_alloc(string, len1*size_of_logical, istat)
    end subroutine l_alloc1
  
    subroutine l_alloc2(string, len1, len2, v)
        !! Allocate a 2-dimensional array of reals
        implicit none

        character(len=*), intent(in) :: string
        !! Human-readable description string of the allocation
        !! operation, just for output purpose.
        integer(ip), intent(in) :: len1, len2
        !! Dimensions of the vector
        logical(lp), allocatable, intent(inout) :: v(:,:)
        !! Vector to allocate

        integer(ip) :: istat

        if(.not. is_init) call memory_init(.false., 0)
        allocate(v(len1, len2), stat=istat)
        call chk_alloc(string, len1*len2*size_of_logical, istat)
    end subroutine l_alloc2

    subroutine chk_alloc(string, lall, istat)
        !! Handles the memory errors (including soft limit)
        !! during memory allocation 
        implicit none

        integer(ip), intent(in) :: lall 
        !! Amount of allocated memory in bytes
        integer(ip), intent(in) :: istat
        !! Status flag from allocate()
        character(len=*), intent(in) :: string
        !! Human-readable description string of the allocation
        !! operation, just for output purpose.

        character(len=OMMP_STR_CHAR_MAX) :: msg
        !! Message string for errors
        
        if(istat /= 0) then
            write(msg, "('Allocation error in subroutine ', a ,'. stat= ', i5)") string, istat
            call fatal_error(msg)
        else if(do_chk_limit .and. usedmem+lall > maxmem) then
            write(msg, "('Allocation error in subroutine ', a ,'. Not enough memory (internal limit ', i8, ' W).')") string, maxmem
            call fatal_error(msg)
        else
            usedmem = usedmem + lall
        end if
    end subroutine chk_alloc

    subroutine r_free1(string, v)
        !! Free a 1-dimensional array of reals

        character(len=*), intent(in) :: string
        !! Human-readable description string of the deallocation
        !! operation, just for output purpose.
        real(rp), allocatable, intent(inout) :: v(:)
        !! Array to free
        
        integer(ip) :: istat, ltot

        if(allocated(v)) then
            ltot = size(v) * size_of_real
            deallocate(v, stat=istat)
            call chk_free(string, ltot, istat)
        end if
    end subroutine r_free1

    subroutine r_free2(string, v)
        !! Free a 2-dimensional array of reals
        
        character(len=*), intent(in) :: string
        !! Human-readable description string of the deallocation
        !! operation, just for output purpose.
        real(rp),  allocatable, intent(inout) :: v(:,:)
        !! Array to free
        
        integer(ip) :: istat, ltot

        if(allocated(v)) then
            ltot = size(v) * size_of_real
            deallocate(v, stat=istat)
            call chk_free(string, ltot, istat)
        end if
    end subroutine r_free2

    subroutine r_free3(string, v)
        !! Free a 3-dimensional array of reals
        
        character (len=*), intent(in) :: string
        !! Human-readable description string of the deallocation
        !! operation, just for output purpose.
        real(rp), allocatable, intent(inout) :: v(:,:,:)
        !! Array to free
        
        integer(ip) :: istat, ltot

        if(allocated(v)) then
            ltot = size(v) * size_of_real
            deallocate(v, stat=istat)
            call chk_free(string, ltot, istat)
        end if
    end subroutine r_free3

    subroutine i_free1(string, v)
        !! Free a 1-dimensional array of integers
        
        character (len=*), intent(in) :: string
        !! Human-readable description string of the deallocation
        !! operation, just for output purpose.
        integer(ip), allocatable, intent(inout) :: v(:)
        !! Array to free
        
        integer(ip) :: istat, ltot
  
        if(allocated(v)) then
            ltot = size(v) * size_of_int
            deallocate(v, stat=istat)
            call chk_free(string, ltot, istat)
        end if
    end subroutine i_free1
  
    subroutine i_free2(string, v)
        !! Free a 2-dimensional array of integers
        
        character (len=*), intent(in) :: string
        !! Human-readable description string of the deallocation
        !! operation, just for output purpose.
        integer(ip), allocatable, intent(inout) :: v(:,:)
        !! Array to free
        
        integer(ip) :: istat, ltot
  
        if(allocated(v)) then
            ltot = size(v) * size_of_int
            deallocate(v, stat=istat)
            call chk_free(string, ltot, istat)
        end if
    end subroutine i_free2
  
    subroutine i_free3(string, v)
        !! Free a 3-dimensional array of integers
        
        character (len=*), intent(in) :: string
        !! Human-readable description string of the deallocation
        !! operation, just for output purpose.
        integer(ip), allocatable, intent(inout) :: v(:,:,:)
        !! Array to free
        
        integer(ip) :: istat, ltot

        if(allocated(v)) then
            ltot = size(v) * size_of_int
            deallocate (v, stat=istat)
            call chk_free(string, ltot, istat)
        end if
    end subroutine i_free3
    
    subroutine l_free1(string, v)
        !! Free a 1-dimensional array of integers
        
        character (len=*), intent(in) :: string
        !! Human-readable description string of the deallocation
        !! operation, just for output purpose.
        logical(lp), allocatable, intent(inout) :: v(:)
        !! Array to free
        
        integer(ip) :: istat, ltot
  
        if(allocated(v)) then
            ltot = size(v) * size_of_logical
            deallocate(v, stat=istat)
            call chk_free(string, ltot, istat)
        end if
    end subroutine l_free1
  
    subroutine l_free2(string, v)
        !! Free a 2-dimensional array of integers
        
        character (len=*), intent(in) :: string
        !! Human-readable description string of the deallocation
        !! operation, just for output purpose.
        logical(lp), allocatable, intent(inout) :: v(:,:)
        !! Array to free
        
        integer(ip) :: istat, ltot
  
        if(allocated(v)) then
            ltot = size(v) * size_of_logical
            deallocate(v, stat=istat)
            call chk_free(string, ltot, istat)
        end if
    end subroutine l_free2
 
    subroutine chk_free(string, lfree, istat)
        !! Handles the memory errors (including soft limits)
        !! during the deallocation
        implicit none

        character(len=*), intent(in) :: string
        !! Human-readable description string of the deallocation
        !! operation, just for output purpose.
        integer(ip), intent(in) :: lfree
        !! amount of memory (in bytes) to free
        integer(ip), intent(in) :: istat
        !! return flag of deallocate
        
        character(len=OMMP_STR_CHAR_MAX) :: msg
        !! Message string for errors

        if(istat /= 0)then
            write(msg, "('Deallocation error in subroutine ', a ,'. stat= ', i5)") string, istat
            call fatal_error(msg)
        else
            usedmem = usedmem - lfree
        end if

    end subroutine chk_free

end module mod_memory
