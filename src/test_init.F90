program test_init
    use iso_c_binding, only: c_char
    use mod_interface, only: w_mmpol_init, write_hdf5

    implicit none
    character(kind=c_char, len=120), dimension(2) :: args
    integer :: hdf5_werr
  
    if (command_argument_count() /= 2) then
        write(6, *) "Syntax expected "
        write(6, *) "   $ test_init.exe <INPUT FILE> <OUTPUT FILE>"
    else
        call get_command_argument(1, args(1))
        call get_command_argument(2, args(2))
        
        write(6, *) "Input file: ", args(1)
        write(6, *) "Output file: ", args(2)
        
        call w_mmpol_init(trim(args(1)))
        hdf5_werr = write_hdf5(trim(args(2)))
    end if
end program test_init
