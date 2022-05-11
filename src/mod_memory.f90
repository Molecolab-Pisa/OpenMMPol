module memory
    use iso_c_binding

    implicit none
    private 

!    #ifdef USE_I8
!        integer(kind=c_int64_t), parameter :: ip = c_int64_t
!    #else
        integer(kind=c_int32_t), parameter :: ip = c_int32_t
!    #endif
    integer(ip), parameter :: rp = c_double
    
    integer(ip) :: maxmem ! Max memory that can be allocated in bytes
    integer(ip) :: size_of_int ! Number of bytes for an integer
    integer(ip) :: size_of_real ! Number of bytes for a real

    public :: rp, ip ! precision for real and integers
    public :: mallocate, mfree, print_memory_info, \
              memory_init
    
    integer(ip), parameter   :: memout = 6
    integer(ip), private     :: istat, ltot

    interface mallocate
        module procedure r_alloc1
        module procedure r_alloc2
        module procedure r_alloc3
        module procedure i_alloc1
        module procedure i_alloc2
        module procedure i_alloc3
    end interface mallocate

    interface mfree
        module procedure r_free1
        module procedure r_free2
        module procedure r_free3
        module procedure i_free1
        module procedure i_free2
        module procedure i_free3
    end interface mfree

    contains

    subroutine print_memory_info()
        implicit none 

        write(6, *) "Precision for integer numbers: ", ip
        write(6, *) "An integer occupies ", size_of_int, " bytes."
        write(6, *) "Precision for real numbers: ", rp
        write(6, *) "An integer occupies ", size_of_real, " bytes."
        write(6, *) "The total memory available is ", maxmem, " bytes."
    end subroutine print_memory_info

    subroutine memory_init(max_bytes)
        implicit none

        integer(ip), intent(in) :: max_bytes
        integer(ip) :: my_int
        real(rp) :: my_real
        intrinsic :: sizeof

        maxmem = max_bytes
        size_of_real = sizeof(my_real)
        size_of_int = sizeof(my_int)
    end subroutine memory_init
    
    subroutine r_alloc1(string, len1, v)
        implicit none

        character(len=*), intent(in) :: string
        integer(ip), intent(in) :: len1
        real(rp), allocatable, intent(inout) :: v(:)

        integer(ip) :: istat
 
        allocate(v(len1), stat=istat)
        call chk_all(string, len1*size_of_real, istat)
    end subroutine r_alloc1
  
    subroutine r_alloc2(string, len1, len2, v)
        implicit none

        character(len=*), intent(in) :: string
        integer(ip), intent(in) :: len1, len2
        real(rp), allocatable, intent(inout) :: v(:,:)

        integer(ip) :: istat

        allocate(v(len1, len2), stat=istat)
        call chk_all(string, len1*len2*size_of_real, istat)
    end subroutine r_alloc2
  
    subroutine r_alloc3(string, len1, len2, len3, v)
        implicit none

        character(len=*), intent(in) :: string
        integer(ip), intent(in) :: len1, len2, len3
        real(rp), allocatable, intent(inout) :: v(:,:,:)
        
        integer(ip) :: istat

        allocate(v(len1, len2, len3), stat=istat)
        call chk_all(string, len1*len2*len3*size_of_real, istat)
    end subroutine r_alloc3

    subroutine i_alloc1(string, len1, v)
        implicit none

        character (len=*), intent(in) :: string
        integer(ip), intent(in) :: len1
        integer(ip), allocatable, intent(inout) :: v(:)

        integer(ip) :: istat 

        allocate(v(len1), stat=istat)
        call chk_all(string, len1*size_of_int, istat)
    end subroutine i_alloc1

    subroutine i_alloc2(string, len1, len2, v)
      implicit none

      character (len=*), intent(in) :: string
      integer(ip), intent(in) :: len1, len2
      integer(ip), allocatable, intent(inout) :: v(:,:)
 
      integer(ip) :: istat 

      allocate(v(len1, len2), stat=istat)
      call chk_all(string, len1*len2*size_of_int, istat)
    end subroutine i_alloc2
 
    subroutine i_alloc3(string, len1, len2, len3, v)
      implicit none

      character (len=*), intent(in) :: string
      integer(ip), intent(in) :: len1, len2, len3
      integer(ip), allocatable, intent(inout) :: v(:,:,:)

      integer(ip) :: istat

      allocate(v(len1, len2, len3), stat=istat)
      call chk_all(string, len1*len2*len3*size_of_int, istat)
    end subroutine i_alloc3

    subroutine chk_all(string,lall,istat)
      implicit none
      integer(ip),       intent(in) :: lall, istat
      character (len=*), intent(in) :: string
!  
!     memory error format:
!  
   9000 format(t3,'allocation error in subroutine ',a,'. stat= ',i5)
   9010 format(t3,'allocation error in subroutine ',a,'.',/,    &
               t3,'not enough memory. ',i8,' words required',/, &
               t3,'                   ',i8,' words available.')
!  
      if (istat.ne.0) then
        write(memout,9000) string, istat
        stop
      else if (lall.gt.maxmem) then
        write(memout,9010) string, lall, maxmem
        stop
      else
        maxmem = maxmem - lall
      end if
      return
    end subroutine chk_all
!  
    subroutine r_free1(string,v)
      character (len=*),              intent(in)    :: string
      real(rp),          allocatable, intent(inout) :: v(:)
!  
      if (.not. allocated(v)) return
      ltot = size(v)
      deallocate (v, stat=istat)
      call chk_free(string,ltot,istat)
      return
    end subroutine r_free1
!  
    subroutine r_free2(string,v)
      character (len=*),              intent(in)    :: string
      real(rp),          allocatable, intent(inout) :: v(:,:)
!  
      if (.not. allocated(v)) return
      ltot = size(v)
      deallocate (v, stat=istat)
      call chk_free(string,ltot,istat)
      return
    end subroutine r_free2
!  
    subroutine r_free3(string,v)
      character (len=*),              intent(in)    :: string
      real(rp),          allocatable, intent(inout) :: v(:,:,:)
!  
      if (.not. allocated(v)) return
      ltot = size(v)
      deallocate (v, stat=istat)
      call chk_free(string,ltot,istat)
      return
    end subroutine r_free3
!  
    subroutine i_free1(string,v)
      character (len=*),              intent(in)    :: string
      integer(ip),       allocatable, intent(inout) :: v(:)
!  
      if (.not. allocated(v)) return
      ltot = size(v)
      deallocate (v, stat=istat)
      call chk_free(string,ltot,istat)
      return
    end subroutine i_free1
!  
    subroutine i_free2(string,v)
      character (len=*),              intent(in)    :: string
      integer(ip),       allocatable, intent(inout) :: v(:,:)
!  
      if (.not. allocated(v)) return
      ltot = size(v)
      deallocate (v, stat=istat)
      call chk_free(string,ltot,istat)
      return
    end subroutine i_free2
!  
    subroutine i_free3(string,v)
      character (len=*),              intent(in)    :: string
      integer(ip),       allocatable, intent(inout) :: v(:,:,:)
!  
      if (.not. allocated(v)) return
      ltot = size(v)
      deallocate (v, stat=istat)
      call chk_free(string,ltot,istat)
      return
    end subroutine i_free3
!  
    subroutine chk_free(string,lfree,istat)
      implicit none
      integer(ip),       intent(in) :: lfree, istat
      character (len=*), intent(in) :: string
!  
!     memory error format:
!  
   9000 format(t3,'deallocation error in subroutine ',a,'. stat= ',i5)
!  
      if (istat.ne.0) then
        write(memout,9000) string, istat
        stop
      else
        maxmem = maxmem + lfree
      end if
      return
    end subroutine chk_free

end module memory
