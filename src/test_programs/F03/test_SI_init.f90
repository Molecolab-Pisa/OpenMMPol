program test_SI_init
    use iso_c_binding, only: c_char
    use ommp_interface

    implicit none
    character(kind=c_char, len=120), dimension(2) :: args
    integer :: narg
    type(ommp_system), pointer :: my_system
    type(ommp_qm_helper), pointer :: my_qmh
  
    narg = command_argument_count()
    if (narg /= 2 .and. narg /= 1) then
        write(6, *) "Syntax expected "
        write(6, *) "   $ test_init.exe <JSON FILE> [<OUTPUT FILE>]"
    else 
        call get_command_argument(1, args(1))
        if(narg == 2) &
            call get_command_argument(2, args(2))
       
        call ommp_smartinput(trim(args(1)), my_system, my_qmh)
        if(narg == 2) then
             call ommp_print_summary_to_file(my_system, &
                                             trim(args(2)))
        else
            call ommp_print_summary(my_system)
        end if
        
        if(associated(my_qmh)) call ommp_terminate_qm_helper(my_qmh)
        if(associated(my_system)) call ommp_terminate(my_system)
    end if
end program test_SI_init
