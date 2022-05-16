program test_init
  use mmpol
  use mod_memory, only: print_memory_info

  implicit none
  character (len=120), dimension(2) :: args
  
  call memory_init(.false., 10000000_ip)
  call print_memory_info()

  if (command_argument_count() /= 1) then
     write(6, *) "Syntax expected "
     write(6, *) "   $ test_init.exe input_file.mmp"
  else
    call get_command_argument(1, args(1))
    input_file = trim(args(1))
    call mmpol_init()
  end if
end program test_init
