module precision
  implicit none
  integer(kind=4), parameter :: rp = 8
  integer(kind=4), parameter :: ip = 8
!#ifdef USE_I8
!  integer(kind=8), parameter :: ip = 8
!#else
!  integer(kind=4), parameter :: ip = 4
!#endif
end module precision
