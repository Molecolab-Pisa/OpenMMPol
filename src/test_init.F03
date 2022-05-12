program test_init
  use mmpol
  use mod_memory, only: print_memory_info

  implicit none
  character (len=120), dimension(2) :: args
!
! the name of the mmpol input file and, optionally, the name of the
! mmpol binary file are provided in input as arguments to the program:
!
  call memory_init(.false., 10000000_ip)
  call print_memory_info()

  if (command_argument_count() .eq. 0) then
     write(6, *) "Syntax expected "
     write(6, *) "   $ test_init.exe input_file.mmp [ scratch_file.rwf ]"
!
!   no argument provided
!
  else if (command_argument_count() .eq. 1) then
    call get_command_argument(1, args(1))
    input_file = trim(args(1))
    scratch_file = input_file(1:len_inname-4)//'.rwf' 
    call mmpol_init()
  else if (command_argument_count() .eq. 2) then
    call get_command_argument(1, args(1))
    call get_command_argument(2, args(2))
    input_file = trim(args(1))
    scratch_file = trim(args(2))
    call mmpol_init()
  end if
!
  
end program test_init
