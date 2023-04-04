program test_potential_xyz
    use iso_c_binding, only: c_char
    use ommp_interface

    implicit none
    character(kind=c_char, len=120), dimension(3) :: args
    integer :: narg
    type(ommp_system), pointer :: my_system
    real(ommp_real) :: eb, ea, eba, eub, eaa, eopb, eopd, eid, eit, et, &
                       ept, ebt, eat, ett, ev, er, edsp, ec, ecd, ed, em, ep, &
                       ect, erxf, es, elf, eg, ex
    real(ommp_real), allocatable :: ef(:,:)
  
    narg = command_argument_count()
    if (narg /= 3) then
        write(6, *) "Syntax expected "
        write(6, *) "   $ test_init.exe <XYZ FILE> <PRM FILE> <OUTPUT FILE>"
    else 
        call get_command_argument(1, args(1))
        call get_command_argument(2, args(2))
        call get_command_argument(3, args(3))
       
        call ommp_init_xyz(my_system, trim(args(1)), trim(args(2)))
        
        allocate(ef(3, my_system%eel%pol_atoms))
        ef = 0.0
        
        call ommp_set_external_field(my_system, ef, OMMP_SOLVER_DEFAULT)
        
        em   = ommp_get_fixedelec_energy(my_system)
        ep   = ommp_get_polelec_energy(my_system)
        ev   = ommp_get_vdw_energy(my_system)
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

        open(unit=101, file=args(3))
        write(101, '(A4, E20.12)') "EM  ", em   * OMMP_AU2KCALMOL 
        write(101, '(A4, E20.12)') "EP  ", ep   * OMMP_AU2KCALMOL 
        write(101, '(A4, E20.12)') "EV  ", ev   * OMMP_AU2KCALMOL 
        write(101, '(A4, E20.12)') "EB  ", eb   * OMMP_AU2KCALMOL 
        write(101, '(A4, E20.12)') "EA  ", ea   * OMMP_AU2KCALMOL 
        write(101, '(A4, E20.12)') "EBA ", eba  * OMMP_AU2KCALMOL 
        write(101, '(A4, E20.12)') "EUB ", eub  * OMMP_AU2KCALMOL 
        write(101, '(A4, E20.12)') "EOPB", eopb * OMMP_AU2KCALMOL 
        write(101, '(A4, E20.12)') "EPT ", ept  * OMMP_AU2KCALMOL 
        write(101, '(A4, E20.12)') "ET  ", et   * OMMP_AU2KCALMOL 
        write(101, '(A4, E20.12)') "ETT ", ett  * OMMP_AU2KCALMOL 
        write(101, '(A4, E20.12)') "EAT ", eat  * OMMP_AU2KCALMOL 
        write(101, '(A4, E20.12)') "EBT ", ebt  * OMMP_AU2KCALMOL 
        write(101, '(A4, E20.12)') "EAA ", eaa  * OMMP_AU2KCALMOL 
        write(101, '(A4, E20.12)') "EOPD", eopd * OMMP_AU2KCALMOL 
        write(101, '(A4, E20.12)') "EID ", eid  * OMMP_AU2KCALMOL 
        write(101, '(A4, E20.12)') "EIT ", eit  * OMMP_AU2KCALMOL 
        write(101, '(A4, E20.12)') "ER  ", er   * OMMP_AU2KCALMOL 
        write(101, '(A4, E20.12)') "EDSP", edsp * OMMP_AU2KCALMOL 
        write(101, '(A4, E20.12)') "EC  ", ec   * OMMP_AU2KCALMOL 
        write(101, '(A4, E20.12)') "ECD ", ecd  * OMMP_AU2KCALMOL 
        write(101, '(A4, E20.12)') "ED  ", ed   * OMMP_AU2KCALMOL 
        write(101, '(A4, E20.12)') "ECT ", ect  * OMMP_AU2KCALMOL 
        write(101, '(A4, E20.12)') "ERXF", erxf * OMMP_AU2KCALMOL 
        write(101, '(A4, E20.12)') "ES  ", es   * OMMP_AU2KCALMOL 
        write(101, '(A4, E20.12)') "ELF ", elf  * OMMP_AU2KCALMOL 
        write(101, '(A4, E20.12)') "EG  ", eg   * OMMP_AU2KCALMOL 
        write(101, '(A4, E20.12)') "EX  ", ex   * OMMP_AU2KCALMOL 
        close(101)
       
        call ommp_terminate(my_system)
    end if
end program test_potential_xyz
