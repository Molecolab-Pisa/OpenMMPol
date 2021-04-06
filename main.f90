program main
  use mmpol
  implicit none
  character (len=120), dimension(2) :: args
!
! the name of the mmpol input file and, optionally, the name of the
! mmpol binary file are provided in input as arguments to the program:
!
  if (iargc() .eq. 0) then
!
!   no argument provided
!
  else if (iargc() .eq. 1) then
    call getarg(1, args(1))
    input_file = trim(args(1))
    scratch_file = input_file(1:len_inname-4)//'.rwf' 
  else if (iargc() .eq. 2) then
    call getarg(1, args(1))
    call getarg(2, args(2))
    input_file = trim(args(1))
    scratch_file = trim(args(2))
  end if
!
  call mmpol_init
  
end program main
