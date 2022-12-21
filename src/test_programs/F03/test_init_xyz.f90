program test_init
    use iso_c_binding, only: c_char
    use ommp_interface

    implicit none
    character(kind=c_char, len=120), dimension(3) :: args
    integer :: narg
    type(ommp_system), pointer :: my_system
  
    narg = command_argument_count()
    if (narg /= 2 .and. narg /= 3) then
        write(6, *) "Syntax expected "
        write(6, *) "   $ test_init.exe <XYZ FILE> <PRM FILE> [<OUTPUT FILE>]"
    else 
        call get_command_argument(1, args(1))
        call get_command_argument(2, args(2))
        if(narg == 3) &
            call get_command_argument(3, args(3))
       
        call ommp_set_verbose(OMMP_VERBOSE_DEBUG)
        call ommp_init_xyz(my_system, trim(args(1)), trim(args(2)))
        if(narg == 3) then
             call ommp_print_summary_to_file(my_system, args(3))
        else
            call ommp_print_summary(my_system)
        end if
        
        call ommp_terminate(my_system)
    end if
end program test_init
