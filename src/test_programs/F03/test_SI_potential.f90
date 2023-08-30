program test_SI_potential
    use iso_c_binding, only: c_char
    use ommp_interface

    implicit none
    character(kind=c_char, len=120), dimension(3) :: args
    character(len=OMMP_STR_CHAR_MAX) :: msg
    integer :: narg
    type(ommp_system), pointer :: my_system
    type(ommp_qm_helper), pointer :: my_qmh
    real(ommp_real) :: eb, ea, eba, eub, eaa, eopb, eopd, eid, eit, et, &
                       ept, ebt, eat, ett, ev, er, edsp, ec, ecd, ed, em, ep, &
                       ect, erxf, es, elf, eg, ex, evqmmm
    real(ommp_real), allocatable :: ef(:,:)
    logical :: use_qm = .false., use_fake_qm = .false., use_ext_ef = .false.
    integer(ommp_integer) :: i, k
  
    narg = command_argument_count()
    if (narg /= 2 .and. narg /= 3) then
        write(6, *) "Syntax expected "
        write(6, *) "   $ test_SI_potential.exe <JSON FILE> <OUTPUT FILE> [<EF_FILE>]"
        call exit(1)
    else 
        call get_command_argument(1, args(1))
        call get_command_argument(2, args(2))
        if(narg == 3) then
            call get_command_argument(3, args(3))
            use_ext_ef = .true.
        end if
       
        call ommp_smartinput(trim(args(1)), my_system, my_qmh)
        call ommp_set_outputfile(trim(args(2)))

        if(associated(my_qmh)) then
            use_qm = .true.
            if(my_system%use_linkatoms) then
                use_fake_qm = .true.
                call ommp_fatal("Smart Input Error!!!")
                call exit(1)
            end if
        end if
        
        ! Read the electric field from file if required
        allocate(ef(3, my_system%eel%pol_atoms))
        ef = 0.0
        if(use_ext_ef) then
            open(unit=101, file=args(3)) 
            do i=1, my_system%top%mm_atoms
                if(my_system%eel%mm_polar(i) > 0) then
                    read(101, *) ef(:, my_system%eel%mm_polar(i))
                else
                    read(101, *)
                end if
            end do
            close(101)
        end if

        em   = ommp_get_fixedelec_energy(my_system)
        call ommp_set_external_field(my_system, ef, OMMP_SOLVER_NONE, OMMP_MATV_NONE)
        
        ep   = ommp_get_polelec_energy(my_system)
        ev   = ommp_get_vdw_energy(my_system)
        if(use_qm) then
            evqmmm = ommp_qm_helper_vdw_energy(my_qmh, my_system)
        else
            evqmmm = 0.0
        end if
        eb   = ommp_get_bond_energy(my_system)
        ea   = ommp_get_angle_energy(my_system)
        eba  = ommp_get_strbnd_energy(my_system)
        eub  = ommp_get_urey_energy(my_system)
        eopb = ommp_get_opb_energy(my_system)
        ept  = ommp_get_pitors_energy(my_system)
        et   = ommp_get_torsion_energy(my_system)
        ett  = ommp_get_tortor_energy(my_system)
        eat  = ommp_get_angtor_energy(my_system)
        ebt  = ommp_get_strtor_energy(my_system)
        eaa  = 0.0
        eopd = 0.0
        eid  = 0.0  
        eit  = ommp_get_imptorsion_energy(my_system) 
        er   = 0.0
        edsp = 0.0
        ec   = 0.0
        ecd  = 0.0
        ed   = 0.0
        ect  = 0.0
        erxf = 0.0
        es   = 0.0
        elf  = 0.0
        eg   = 0.0
        ex   = 0.0
        
        write(msg, '(A4, E20.12)') "EM  ", em   
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        write(msg, '(A4, E20.12)') "EP  ", ep   
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        write(msg, '(A4, E20.12)') "EV  ", ev   
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        write(msg, '(A4, E20.12)') "EVQMMM  ", evqmmm 
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        write(msg, '(A4, E20.12)') "EB  ", eb   
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        write(msg, '(A4, E20.12)') "EA  ", ea   
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        write(msg, '(A4, E20.12)') "EBA ", eba  
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        write(msg, '(A4, E20.12)') "EUB ", eub  
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        write(msg, '(A4, E20.12)') "EOPB", eopb 
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        write(msg, '(A4, E20.12)') "EPT ", ept  
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        write(msg, '(A4, E20.12)') "ET  ", et   
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        write(msg, '(A4, E20.12)') "ETT ", ett  
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        write(msg, '(A4, E20.12)') "EAT ", eat  
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        write(msg, '(A4, E20.12)') "EBT ", ebt  
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        write(msg, '(A4, E20.12)') "EAA ", eaa  
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        write(msg, '(A4, E20.12)') "EOPD", eopd 
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        write(msg, '(A4, E20.12)') "EID ", eid  
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        write(msg, '(A4, E20.12)') "EIT ", eit  
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        write(msg, '(A4, E20.12)') "ER  ", er   
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        write(msg, '(A4, E20.12)') "EDSP", edsp 
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        write(msg, '(A4, E20.12)') "EC  ", ec   
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        write(msg, '(A4, E20.12)') "ECD ", ecd  
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        write(msg, '(A4, E20.12)') "ED  ", ed   
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        write(msg, '(A4, E20.12)') "ECT ", ect  
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        write(msg, '(A4, E20.12)') "ERXF", erxf 
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        write(msg, '(A4, E20.12)') "ES  ", es   
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        write(msg, '(A4, E20.12)') "ELF ", elf  
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        write(msg, '(A4, E20.12)') "EG  ", eg   
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        write(msg, '(A4, E20.12)') "EX  ", ex   
        call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE")
        
        do k=1, my_system%eel%n_ipd
            do i=1, my_system%eel%pol_atoms
                write(msg, "(E20.12, ' ', E20.12, ' ', E20.12)") my_system%eel%ipd(:,i,k)
                call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-IPD")
            end do
        end do
        if(associated(my_qmh)) call ommp_terminate_qm_helper(my_qmh)
        if(associated(my_system)) call ommp_terminate(my_system)
    end if
end program test_SI_potential
