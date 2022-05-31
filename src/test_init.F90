program test_init
  use mod_inputloader, only: mmpol_init_from_mmp
  use mod_memory, only: print_memory_info, memory_init, ip

  implicit none
  character (len=120), dimension(2) :: args
  
  call memory_init(.false., 10000000_ip)
  call print_memory_info()

  if (command_argument_count() /= 1) then
     write(6, *) "Syntax expected "
     write(6, *) "   $ test_init.exe input_file.mmp"
  else
    call get_command_argument(1, args(1))
    call mmpol_init_from_mmp(trim(args(1)))
  end if
end program test_init
