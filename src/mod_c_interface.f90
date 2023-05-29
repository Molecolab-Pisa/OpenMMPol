! Wrapper function for open-mmpol library
module mod_ommp_C_interface
    !! The interface of the library, basically all the operation performed
    !! by an external code should be done through the routines of this
    !! module. 
    !! The interface is conceived to work naturally with C and Fortran; the C
    !! interface is also used to build the interface for Python.
    !! In a fortran code, this module can be imported and it should expose 
    !! directly all the vector and scalar quantities needed.
    !! In a C code, routines are provided to get the pointer or the values of 
    !! vector and scalar quantites respectively.

    use iso_c_binding
    use ommp_interface
    use mod_constants, only: OMMP_STR_CHAR_MAX

    implicit none
    
    private :: c2f_string, OMMP_STR_CHAR_MAX

    contains
        !! Internal utilities for Fortran -> C interface
        pure subroutine c2f_string(c_str, f_str)
            !! Convert a string coming from C into a Fortran string
            implicit none
            
            character(kind=c_char), intent(in) :: c_str(:)
            !! Input string to be converted
            character(len=*), intent(out) :: f_str

            integer :: i 

            i = 1
            do while(c_str(i) /= c_null_char)
                f_str(i:i) = c_str(i)
                i = i + 1
            end do

            do i = i, len(f_str)
                f_str(i:i) = ' '
            end do

            f_str = trim(f_str)
        end subroutine c2f_string

        ! Functions directly mapped on OMMP internal functions (which are 
        ! exposed on fortran side by 
        !     use mod_xxx, only a => b
        subroutine C_ommp_set_verbose(verb) bind(c, name='ommp_set_verbose')
            !! Set the verbosity level of the library to verb
            implicit none 

            integer(ommp_integer), intent(in), value :: verb
            
            !! Requested verbosityi of library
            call ommp_set_verbose(verb)
        end subroutine C_ommp_set_verbose

        subroutine C_ommp_print_summary(s_prt) bind(c, name='ommp_print_summary')
            !! Print a summary of the system input on standard output.
            implicit none
            
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
           
            call c_f_pointer(s_prt, s)
            call ommp_print_summary(s)

        end subroutine C_ommp_print_summary
        
        subroutine C_ommp_print_summary_to_file(s_prt, filename) &
                bind(c, name='ommp_print_summary_to_file')
            !! Print a summary of the system input on file.
            implicit none
            
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            character(kind=c_char), intent(in) :: filename(OMMP_STR_CHAR_MAX)
            !! File where the summary will be printed
            character(len=OMMP_STR_CHAR_MAX) :: output_file
            
            call c_f_pointer(s_prt, s)
            call c2f_string(filename, output_file)
            call ommp_print_summary_to_file(s, output_file)

        end subroutine C_ommp_print_summary_to_file
        
        subroutine C_ommp_save_mmp(s_prt, filename, version) &
                   bind(c, name='ommp_save_mmp')
            implicit none
            type(c_ptr), value :: s_prt
            character(kind=c_char), intent(in) :: filename(OMMP_STR_CHAR_MAX)
            integer(ommp_integer), value :: version

            character(len=OMMP_STR_CHAR_MAX) :: output_file
            type(ommp_system), pointer :: s

            call c_f_pointer(s_prt, s)

            call c2f_string(filename, output_file)
            call ommp_save_mmp(s, output_file, version)
        end subroutine 
        
        subroutine C_ommp_update_coordinates(s_prt, new_c_p) &
                bind(C, name='ommp_update_coordinates')
            use mod_mmpol, only: update_coordinates
            
            implicit none
            
            type(c_ptr), value :: s_prt
            type(c_ptr), value :: new_c_p
            
            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: new_c(:,:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(new_c_p, new_c, [3, s%top%mm_atoms])
            
            call update_coordinates(s, new_c)
        end subroutine
        
        ! Functions mapped on actual Fortran interface functions
        ! OMMP System Object housekeeping
        function C_ommp_init_mmp(filename) &
                result(c_prt) bind(c, name='ommp_init_mmp')
           !! Initalize OMMP System Object from .mmp file  
            implicit none

            type(ommp_system), pointer :: s
            character(kind=c_char), intent(in) :: filename(OMMP_STR_CHAR_MAX)
            character(len=OMMP_STR_CHAR_MAX) :: input_file
            type(c_ptr) :: c_prt
            
            allocate(s)
            
            call c2f_string(filename, input_file)
            call ommp_init_mmp(s, input_file)
            c_prt = c_loc(s)
        end function
        
        function C_ommp_init_xyz(xyzfile, prmfile) &
                result(c_prt) bind(c, name='ommp_init_xyz')
            !! Initialize the library using a Tinker xyz and a Tinker prm
            implicit none
            
            type(ommp_system), pointer :: s
            character(kind=c_char), intent(in) :: xyzfile(OMMP_STR_CHAR_MAX), & 
                                                  prmfile(OMMP_STR_CHAR_MAX)
            character(len=OMMP_STR_CHAR_MAX) :: xyz_file, prm_file
            type(c_ptr) :: c_prt

            allocate(s)
            
            call c2f_string(prmfile, prm_file)
            call c2f_string(xyzfile, xyz_file)
            call ommp_init_xyz(s, xyz_file, prm_file)
            c_prt = c_loc(s)
        end function
        
        subroutine C_ommp_set_frozen_atoms(s_prt, n, frozen) &
                bind(c, name='ommp_set_frozen_atoms')
            implicit none
            type(c_ptr), value :: s_prt, frozen
            integer(ommp_integer), value :: n
            type(ommp_system), pointer :: s
            integer(ommp_integer), pointer :: f(:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(frozen, f, [n])
            f = f + 1
            call ommp_set_frozen_atoms(s, n, f)
            f = f - 1
        end subroutine
        
        subroutine C_ommp_turn_pol_off(s_prt, n, nopol) &
                bind(c, name='ommp_turn_pol_off')
            implicit none
            type(c_ptr), value :: s_prt, nopol
            integer(ommp_integer), value :: n
            type(ommp_system), pointer :: s
            integer(ommp_integer), pointer :: f(:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(nopol, f, [n])
            f = f + 1
            call ommp_turn_pol_off(s, n, f)
            f = f - 1
        end subroutine
        
        subroutine C_ommp_terminate(s_prt) bind(c, name='ommp_terminate')
            !! Terminate a OMMP System Object
            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s

            call c_f_pointer(s_prt, s)
            call ommp_terminate(s)
        end subroutine
       
        ! Interface for normal operation of OMMP System Object
        subroutine C_ommp_set_external_field(s_prt, ext_field_prt, solver) &
                bind(c, name='ommp_set_external_field')
            implicit none
            
            type(c_ptr), value :: s_prt
            type(c_ptr), value :: ext_field_prt
            integer(ommp_integer), intent(in), value :: solver
            
            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: ext_field(:,:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(ext_field_prt, ext_field, [3, s%eel%pol_atoms])

            call ommp_set_external_field(s, ext_field, solver, .true.)
        end subroutine C_ommp_set_external_field
        
        subroutine C_ommp_set_external_field_nomm(s_prt, ext_field_prt, solver) &
                bind(c, name='ommp_set_external_field_nomm')
            !!use mod_mmpol, only: pol_atoms
            
            implicit none
            
            type(c_ptr), value :: s_prt
            type(c_ptr), value :: ext_field_prt
            integer(ommp_integer), intent(in), value :: solver
            
            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: ext_field(:,:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(ext_field_prt, ext_field, [3, s%eel%pol_atoms])
            
            call ommp_set_external_field(s, ext_field, solver, .false.)
        end subroutine C_ommp_set_external_field_nomm

        subroutine C_ommp_potential_mmpol2ext(s_prt, n, cext, v) &
                bind(c, name='ommp_potential_mmpol2ext')
            ! Compute the electric potential of static sites at
            ! arbitrary coordinates
            implicit none
            
            integer(ommp_integer), intent(in), value :: n
            type(c_ptr), value :: s_prt, cext, v
            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: fcext(:,:), fv(:)
           
            call c_f_pointer(s_prt, s)
            call c_f_pointer(cext, fcext, [3,n])
            call c_f_pointer(v, fv, [n])
            call ommp_potential_mmpol2ext(s, n, fcext, fv)
        end subroutine

        subroutine C_ommp_potential_pol2ext(s_prt, n, cext, v) &
                bind(c, name='ommp_potential_pol2ext')
            ! Compute the electric potential of static sites at
            ! arbitrary coordinates
            implicit none
            
            integer(ommp_integer), intent(in), value :: n
            type(c_ptr), value :: s_prt, cext, v
            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: fcext(:,:), fv(:)
           
            call c_f_pointer(s_prt, s)
            call c_f_pointer(cext, fcext, [3,n])
            call c_f_pointer(v, fv, [n])
            call ommp_potential_pol2ext(s, n, fcext, fv)
        end subroutine
        
        subroutine C_ommp_potential_mm2ext(s_prt, n, cext, v) &
                bind(c, name='ommp_potential_mm2ext')
            ! Compute the electric potential of static sites at
            ! arbitrary coordinates
            implicit none
            
            integer(ommp_integer), intent(in), value :: n
            type(c_ptr), value :: s_prt, cext, v
            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: fcext(:,:), fv(:)
           
            call c_f_pointer(s_prt, s)
            call c_f_pointer(cext, fcext, [3,n])
            call c_f_pointer(v, fv, [n])
            call ommp_potential_mm2ext(s, n, fcext, fv)
        end subroutine
        
        function C_ommp_get_polelec_energy(s_prt) &
                result(epol) bind(c, name='ommp_get_polelec_energy')
            implicit none
            
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            real(ommp_real) :: epol

            call c_f_pointer(s_prt, s)
            epol = ommp_get_polelec_energy(s)
        end function

        function C_ommp_get_fixedelec_energy(s_prt) &
                result(emm) bind(c, name='ommp_get_fixedelec_energy')
            ! Get the interaction energy of fixed multipoles
            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            real(ommp_real) :: emm

            call c_f_pointer(s_prt, s)

            emm = ommp_get_fixedelec_energy(s)
        end function
        
        function C_ommp_get_full_ele_energy(s_prt) &
                result(ene) bind(c, name='ommp_get_full_ele_energy')
            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: sys_obj
            real(ommp_real) :: ene

            call c_f_pointer(s_prt, sys_obj)

            ene = ommp_get_full_ele_energy(sys_obj) 
        end function
        
        function C_ommp_get_vdw_energy(s_prt) &
                result(evdw) bind(c, name='ommp_get_vdw_energy')
            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            real(ommp_real) :: evdw

            call c_f_pointer(s_prt, s)

            evdw = ommp_get_vdw_energy(s)
        
        end function
        
        function C_ommp_get_bond_energy(s_prt) &
                result(ebnd) bind(c, name='ommp_get_bond_energy')
            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            real(ommp_real) :: ebnd

            call c_f_pointer(s_prt, s)

            ebnd = ommp_get_bond_energy(s)
        end function
        
        function C_ommp_get_angle_energy(s_prt) &
                result(eang) bind(c, name='ommp_get_angle_energy')
            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            real(ommp_real) :: eang

            call c_f_pointer(s_prt, s)

            eang = ommp_get_angle_energy(s)
        end function
        
        function C_ommp_get_strbnd_energy(s_prt) &
                result(eba) bind(c, name='ommp_get_strbnd_energy')
            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            real(ommp_real) :: eba

            call c_f_pointer(s_prt, s)

            eba = ommp_get_strbnd_energy(s)
        end function
        
        function C_ommp_get_urey_energy(s_prt) &
                result(eub) bind(c, name='ommp_get_urey_energy')
            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            real(ommp_real) :: eub

            call c_f_pointer(s_prt, s)

            eub = ommp_get_urey_energy(s)
        end function
        
        function C_ommp_get_opb_energy(s_prt) &
            result(eopb) bind(c, name='ommp_get_opb_energy')
            
            use mod_bonded, only: opb_potential

            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            real(ommp_real) :: eopb

            call c_f_pointer(s_prt, s)

            eopb = ommp_get_opb_energy(s)
        end function
        
        function C_ommp_get_imptorsion_energy(s_prt) &
                result(et) bind(c, name='ommp_get_imptorsion_energy')
            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            real(ommp_real) :: et

            call c_f_pointer(s_prt, s)

            et = ommp_get_imptorsion_energy(s)
        end function
        
        function C_ommp_get_torsion_energy(s_prt) &
                result(et) bind(c, name='ommp_get_torsion_energy')
            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            real(ommp_real) :: et

            call c_f_pointer(s_prt, s)

            et = ommp_get_torsion_energy(s)
        end function
        
        function C_ommp_get_pitors_energy(s_prt) &
                result(epitors) bind(c, name='ommp_get_pitors_energy')
            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            real(ommp_real) :: epitors

            call c_f_pointer(s_prt, s)

            epitors = ommp_get_pitors_energy(s)
        end function
        
        function C_ommp_get_strtor_energy(s_prt) &
                result(ebt) bind(c, name='ommp_get_strtor_energy')
            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            real(ommp_real) :: ebt

            call c_f_pointer(s_prt, s)

            ebt = ommp_get_strtor_energy(s)
        end function
        
        function C_ommp_get_angtor_energy(s_prt) &
                result(eat) bind(c, name='ommp_get_angtor_energy')
            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            real(ommp_real) :: eat

            call c_f_pointer(s_prt, s)

            eat = ommp_get_angtor_energy(s)
        end function
        
        function C_ommp_get_tortor_energy(s_prt) &
                result(ett) bind(c, name='ommp_get_tortor_energy')
            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            real(ommp_real) :: ett

            call c_f_pointer(s_prt, s)

            ett = ommp_get_tortor_energy(s)
        end function
        
        function C_ommp_get_full_bnd_energy(s_prt) &
                result(ene) bind(c, name='ommp_get_full_bnd_energy')
            
            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: sys_obj
            real(ommp_real) :: ene

            call c_f_pointer(s_prt, sys_obj)

            ene = ommp_get_full_bnd_energy(sys_obj) 
        end function
        
        function C_ommp_get_full_energy(s_prt) &
                result(ene) bind(c, name='ommp_get_full_energy')
            
            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: sys_obj
            real(ommp_real) :: ene

            call c_f_pointer(s_prt, sys_obj)

            ene = ommp_get_full_energy(sys_obj)
        end function

        ! Functions for advanced operation and gradients
        subroutine C_ommp_vdw_geomgrad(s_prt, grd_prt) &
                bind(C, name='ommp_vdw_geomgrad')
            implicit none
            
            type(c_ptr), value :: s_prt
            type(c_ptr), value :: grd_prt
            
            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: grd(:,:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(grd_prt, grd, [3, s%top%mm_atoms])
            
            call ommp_vdw_geomgrad(s, grd)    
        end subroutine
        
        subroutine C_ommp_rotation_geomgrad(s_prt, pE, pE_grd, grd_prt) &
                bind(C, name="ommp_rotation_geomgrad")
            implicit none

            type(c_ptr), value :: s_prt
            type(c_ptr), value :: pE
            type(c_ptr), value :: pE_grd
            type(c_ptr), value :: grd_prt
            
            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: grd(:,:), E(:,:), Egrd(:,:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(grd_prt, grd, [3, s%top%mm_atoms])
            call c_f_pointer(pE, E, [3, s%top%mm_atoms])
            call c_f_pointer(pE_grd, Egrd, [6, s%top%mm_atoms])
            
            call ommp_rotation_geomgrad(s, E, Egrd, grd)
        end subroutine
        
        subroutine C_ommp_bond_geomgrad(s_prt, grd_prt) &
                bind(C, name='ommp_bond_geomgrad')
            implicit none
            
            type(c_ptr), value :: s_prt
            type(c_ptr), value :: grd_prt
            
            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: grd(:,:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(grd_prt, grd, [3, s%top%mm_atoms])
            
            call ommp_bond_geomgrad(s, grd)    
        end subroutine
        
        subroutine C_ommp_angle_geomgrad(s_prt, grd_prt) &
                bind(C, name='ommp_angle_geomgrad')
            implicit none
            
            type(c_ptr), value :: s_prt
            type(c_ptr), value :: grd_prt
            
            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: grd(:,:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(grd_prt, grd, [3, s%top%mm_atoms])
            
            call ommp_angle_geomgrad(s, grd)    
        end subroutine
        
        subroutine C_ommp_strbnd_geomgrad(s_prt, grd_prt) &
                bind(C, name='ommp_strbnd_geomgrad')
            implicit none
            
            type(c_ptr), value :: s_prt
            type(c_ptr), value :: grd_prt
            
            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: grd(:,:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(grd_prt, grd, [3, s%top%mm_atoms])
            
            call ommp_strbnd_geomgrad(s, grd)    
        end subroutine
        
        subroutine C_ommp_urey_geomgrad(s_prt, grd_prt) &
                bind(C, name='ommp_urey_geomgrad')
            implicit none
            
            type(c_ptr), value :: s_prt
            type(c_ptr), value :: grd_prt
            
            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: grd(:,:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(grd_prt, grd, [3, s%top%mm_atoms])
            
            call ommp_urey_geomgrad(s, grd)    
        end subroutine
        
        subroutine C_ommp_torsion_geomgrad(s_prt, grd_prt) &
                bind(C, name='ommp_torsion_geomgrad')
            implicit none
            
            type(c_ptr), value :: s_prt
            type(c_ptr), value :: grd_prt
            
            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: grd(:,:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(grd_prt, grd, [3, s%top%mm_atoms])
            
            call ommp_torsion_geomgrad(s, grd)    
        end subroutine
        
        subroutine C_ommp_imptorsion_geomgrad(s_prt, grd_prt) &
                bind(C, name='ommp_imptorsion_geomgrad')
            implicit none
            
            type(c_ptr), value :: s_prt
            type(c_ptr), value :: grd_prt
            
            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: grd(:,:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(grd_prt, grd, [3, s%top%mm_atoms])
            
            call ommp_imptorsion_geomgrad(s, grd)    
        end subroutine
        
        subroutine C_ommp_angtor_geomgrad(s_prt, grd_prt) &
                bind(C, name='ommp_angtor_geomgrad')
            implicit none
            
            type(c_ptr), value :: s_prt
            type(c_ptr), value :: grd_prt
            
            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: grd(:,:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(grd_prt, grd, [3, s%top%mm_atoms])
            
            call ommp_angtor_geomgrad(s, grd)    
        end subroutine
        
        subroutine C_ommp_opb_geomgrad(s_prt, grd_prt) &
                bind(C, name='ommp_opb_geomgrad')
            implicit none
            
            type(c_ptr), value :: s_prt
            type(c_ptr), value :: grd_prt
            
            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: grd(:,:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(grd_prt, grd, [3, s%top%mm_atoms])
            
            call ommp_opb_geomgrad(s, grd)    
        end subroutine
        
        subroutine C_ommp_strtor_geomgrad(s_prt, grd_prt) &
                bind(C, name='ommp_strtor_geomgrad')
            implicit none
            
            type(c_ptr), value :: s_prt
            type(c_ptr), value :: grd_prt
            
            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: grd(:,:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(grd_prt, grd, [3, s%top%mm_atoms])
            
            call ommp_strtor_geomgrad(s, grd)    
        end subroutine
        
        subroutine C_ommp_tortor_geomgrad(s_prt, grd_prt) &
                bind(C, name='ommp_tortor_geomgrad')
            implicit none
            
            type(c_ptr), value :: s_prt
            type(c_ptr), value :: grd_prt
            
            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: grd(:,:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(grd_prt, grd, [3, s%top%mm_atoms])
            
            call ommp_tortor_geomgrad(s, grd)    
        end subroutine
         
        subroutine C_ommp_pitors_geomgrad(s_prt, grd_prt) &
                bind(C, name='ommp_pitors_geomgrad')
            implicit none
            
            type(c_ptr), value :: s_prt
            type(c_ptr), value :: grd_prt
            
            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: grd(:,:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(grd_prt, grd, [3, s%top%mm_atoms])
            
            call ommp_pitors_geomgrad(s, grd)    
        end subroutine
         
        subroutine C_ommp_full_bnd_geomgrad(s_prt, grd_prt) &
                bind(C, name='ommp_full_bnd_geomgrad')
            implicit none
            
            type(c_ptr), value :: s_prt
            type(c_ptr), value :: grd_prt
            
            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: grd(:,:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(grd_prt, grd, [3, s%top%mm_atoms])
            
            call ommp_full_bnd_geomgrad(s, grd)    
        end subroutine

        subroutine C_ommp_fixedelec_geomgrad(s_prt, grd_prt) &
                bind(C, name='ommp_fixedelec_geomgrad')
            implicit none
            
            type(c_ptr), value :: s_prt
            type(c_ptr), value :: grd_prt
            
            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: grd(:,:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(grd_prt, grd, [3, s%top%mm_atoms])
            call ommp_fixedelec_geomgrad(s, grd)
        end subroutine
        
        subroutine C_ommp_polelec_geomgrad(s_prt, grd_prt) &
                bind(C, name='ommp_polelec_geomgrad')
            implicit none
            
            type(c_ptr), value :: s_prt
            type(c_ptr), value :: grd_prt
            
            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: grd(:,:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(grd_prt, grd, [3, s%top%mm_atoms])
            call ommp_polelec_geomgrad(s, grd)
        end subroutine
        
        subroutine C_ommp_full_geomgrad(s_prt, grd_prt) &
                bind(C, name='ommp_full_geomgrad')
            use ommp_interface, only: ommp_full_geomgrad

            implicit none
            
            type(c_ptr), value :: s_prt
            type(c_ptr), value :: grd_prt
            
            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: grd(:,:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(grd_prt, grd, [3, s%top%mm_atoms])
            
            call ommp_full_geomgrad(s, grd)
        end subroutine

#ifdef USE_HDF5
        function C_ommp_init_hdf5(filename, namespace) &
                result(c_prt) bind(c, name='ommp_init_hdf5')
            !! This function is an interface for saving an HDF5 file 
            !! with all the data contained in mmpol module using
            !! [[mod_io:mmpol_save_as_hdf5]]
            use mod_iohdf5, only: mmpol_init_from_hdf5
            
            implicit none
            
            type(ommp_system), pointer :: s
            character(kind=c_char), intent(in) :: filename(OMMP_STR_CHAR_MAX), &
                                                  namespace(OMMP_STR_CHAR_MAX)
            character(len=OMMP_STR_CHAR_MAX) :: hdf5in, nms
            integer(ommp_integer) :: ok
            type(c_ptr) :: c_prt

            allocate(s)

            call c2f_string(filename, hdf5in)
            call c2f_string(namespace, nms)
            call mmpol_init_from_hdf5(hdf5in, trim(nms), s, ok)
            
            c_prt = c_loc(s)
            
        end function C_ommp_init_hdf5
        
        subroutine C_ommp_save_as_hdf5(s_prt, filename, namespace) &
                bind(c, name='ommp_save_as_hdf5')
            
            use mod_iohdf5, only: save_system_as_hdf5 

            implicit none
            
            character(kind=c_char), intent(in) :: filename(OMMP_STR_CHAR_MAX), &
                                                  namespace(OMMP_STR_CHAR_MAX)
            character(len=OMMP_STR_CHAR_MAX) :: hdf5out, nms
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            integer(kind=4) :: err

            call c_f_pointer(s_prt, s)

            call c2f_string(filename, hdf5out)
            call c2f_string(namespace, nms)
            call save_system_as_hdf5(hdf5out, s, err, trim(nms), .false.)
            
        end subroutine C_ommp_save_as_hdf5
        
        subroutine C_ommp_checkpoint(s_prt, filename, namespace) &
                bind(c, name='ommp_checkpoint')
            
            use mod_iohdf5, only: save_system_as_hdf5 

            implicit none
            
            character(kind=c_char), intent(in) :: filename(OMMP_STR_CHAR_MAX), &
                                                  namespace(OMMP_STR_CHAR_MAX)
            character(len=OMMP_STR_CHAR_MAX) :: hdf5out, nms
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            integer(kind=4) :: err

            call c_f_pointer(s_prt, s)

            call c2f_string(filename, hdf5out)
            call c2f_string(namespace, nms)
            call save_system_as_hdf5(hdf5out, s, err, trim(nms), .true.)
            
        end subroutine C_ommp_checkpoint
#endif
        ! Functions to provide direct access to Fortran objects/memory from
        ! C and derived codes.
        function C_ommp_get_cmm(s_prt) bind(c, name='ommp_get_cmm')
            !! Return the c-pointer to the array containing the coordinates of
            !! MM atoms.
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            type(c_ptr) :: C_ommp_get_cmm

            call c_f_pointer(s_prt, s)
            C_ommp_get_cmm = c_loc(s%top%cmm)
        end function C_ommp_get_cmm

        function C_ommp_get_zmm(s_prt) bind(c, name='ommp_get_zmm')
            !! Return the c-pointer to the array containing the coordinates of
            !! MM atoms.
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            type(c_ptr) :: C_ommp_get_zmm

            call c_f_pointer(s_prt, s)
            C_ommp_get_zmm = c_loc(s%top%atz)
        end function C_ommp_get_zmm

        function C_ommp_get_cpol(s_prt) bind(c, name='ommp_get_cpol')
            !! Return the c-pointer to the array containing the coordinates of
            !! polarizable atoms.
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            type(c_ptr) :: C_ommp_get_cpol

            call c_f_pointer(s_prt, s)
            C_ommp_get_cpol = c_loc(s%eel%cpol)
        end function C_ommp_get_cpol

        function C_ommp_get_q(s_prt) bind(c, name='ommp_get_q')
            !! Return the c-pointer to the array containing the static source of 
            !! the electrostatic field.
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            type(c_ptr) :: C_ommp_get_q

            call c_f_pointer(s_prt, s)
            C_ommp_get_q = c_loc(s%eel%q)
        end function C_ommp_get_q

        function C_ommp_get_ipd(s_prt) bind(c, name='ommp_get_ipd')
            !! Return the c-pointer to the array containing the induced dipoles 
            !! on polarizable sites.
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            type(c_ptr) :: C_ommp_get_ipd

            call c_f_pointer(s_prt, s)
            C_ommp_get_ipd = c_loc(s%eel%ipd)
        end function C_ommp_get_ipd
        
        function C_ommp_get_polar_mm(s_prt) bind(c, name='ommp_get_polar_mm')
            !! Return the c-pointer to the array containing the map from 
            !! polarizable to MM atoms.
            use mod_memory, only: mallocate
            implicit none

            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            type(c_ptr) :: C_ommp_get_polar_mm
            
            call c_f_pointer(s_prt, s)
            if(.not. allocated(s%eel%C_polar_mm)) then 
                call mallocate('C_ommp_get_polar_mm [C_polar_mm]', &
                               size(s%eel%polar_mm), s%eel%C_polar_mm)
                s%eel%C_polar_mm = s%eel%polar_mm - 1
            end if

            C_ommp_get_polar_mm = c_loc(s%eel%C_polar_mm)
        end function C_ommp_get_polar_mm
        
        function C_ommp_use_frozen(s_prt) bind(c, name='ommp_use_frozen')
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            logical(c_bool) :: C_ommp_use_frozen

            call c_f_pointer(s_prt, s)
            C_ommp_use_frozen = s%top%use_frozen
        end function C_ommp_use_frozen

        function C_ommp_get_frozen(s_prt) bind(c, name='ommp_get_frozen')
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            type(c_ptr) :: C_ommp_get_frozen

            call c_f_pointer(s_prt, s)
            C_ommp_get_frozen = c_loc(s%top%frozen)
        end function C_ommp_get_frozen

        function C_ommp_get_mm_atoms(s_prt) bind(c, name='ommp_get_mm_atoms')
            !! Return the number of MM atoms in the system.
            implicit none

            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            integer(ommp_integer) :: C_ommp_get_mm_atoms

            call c_f_pointer(s_prt, s)
            C_ommp_get_mm_atoms = s%top%mm_atoms
        end function C_ommp_get_mm_atoms
        
        function C_ommp_get_pol_atoms(s_prt) bind(c, name='ommp_get_pol_atoms')
            !! Return the number of polarizable atoms in the system.
            implicit none

            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            integer(ommp_integer) :: C_ommp_get_pol_atoms
            
            call c_f_pointer(s_prt, s)

            C_ommp_get_pol_atoms = s%eel%pol_atoms
        end function C_ommp_get_pol_atoms

        function C_ommp_get_n_ipd(s_prt) bind(c, name='ommp_get_n_ipd')
            !! Return the number of dipole's set for the current Force-Field.
            implicit none

            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            integer(ommp_integer) :: C_ommp_get_n_ipd

            call c_f_pointer(s_prt, s)
            C_ommp_get_n_ipd = s%eel%n_ipd
        end function C_ommp_get_n_ipd

        function C_ommp_get_ld_cart(s_prt) bind(c, name='ommp_get_ld_cart')
            implicit none

            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            integer(ommp_integer) :: C_ommp_get_ld_cart

            call c_f_pointer(s_prt, s)
            C_ommp_get_ld_cart = s%eel%ld_cart
        end function C_ommp_get_ld_cart

        function C_ommp_ff_is_amoeba(s_prt) bind(c, name='ommp_ff_is_amoeba')
            !! Return true if the current forcefield is AMOEBA, and false in
            !! all other cases.
            implicit none

            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            logical(c_bool) :: C_ommp_ff_is_amoeba

            call c_f_pointer(s_prt, s)
            C_ommp_ff_is_amoeba = s%amoeba
        end function C_ommp_ff_is_amoeba
        
        function C_ommp_use_linkatoms(s_ptr) &
                result(u) bind(C, name='ommp_use_linkatoms')
            implicit none

            type(c_ptr), value :: s_ptr
            !! C pointer to system object
            type(ommp_system), pointer :: s
            logical(c_bool) :: u

            call c_f_pointer(s_ptr, s)
            u = s%use_linkatoms
        end function
        
        !??
        subroutine C_ommp_field_mmpol2ext(s_prt, n, cext, E) &
                bind(c, name='ommp_field_mmpol2ext')
            ! Compute the electric potential of static sites at
            ! arbitrary coordinates
            implicit none
            
            integer(ommp_integer), intent(in), value :: n
            type(c_ptr), value :: s_prt, cext, E
            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: fcext(:,:), fE(:,:)
           
            call c_f_pointer(s_prt, s)
            call c_f_pointer(cext, fcext, [3,n])
            call c_f_pointer(E, fE, [3,n])
            call ommp_field_mmpol2ext(s, n, fcext, fE)
        end subroutine
        
        subroutine C_ommp_field_mm2ext(s_prt, n, cext, E) &
                bind(c, name='ommp_field_mm2ext')
            ! Compute the electric potential of static sites at
            ! arbitrary coordinates
            use mod_electrostatics, only: field_M2E

            implicit none
            
            integer(ommp_integer), intent(in), value :: n
            type(c_ptr), value :: s_prt, cext, E
            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: fcext(:,:), fE(:,:)
           
            call c_f_pointer(s_prt, s)
            call c_f_pointer(cext, fcext, [3,n])
            call c_f_pointer(E, fE, [3,n])
            call field_M2E(s%eel, fcext, fE)
        end subroutine

        subroutine C_ommp_field_pol2ext(s_prt, n, cext, E) &
                bind(c, name='ommp_field_pol2ext')
            ! Compute the electric potential of static sites at
            ! arbitrary coordinates
            use mod_electrostatics, only: field_D2E

            implicit none
            
            integer(ommp_integer), intent(in), value :: n
            type(c_ptr), value :: s_prt, cext, E
            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: fcext(:,:), fE(:,:)
           
            call c_f_pointer(s_prt, s)
            call c_f_pointer(cext, fcext, [3,n])
            call c_f_pointer(E, fE, [3,n])
            call field_D2E(s%eel, fcext, fE)
        end subroutine

        ! Interface for QM Helper module
        function C_ommp_init_qm_helper(n, cqm, qqm, zqm) &
                result(c_prt) bind(c, name='ommp_init_qm_helper')
            implicit none

            type(ommp_qm_helper), pointer :: s
            integer(ommp_integer), value, intent(in) :: n
            type(c_ptr), value, intent(in) :: cqm, qqm, zqm
            
            real(ommp_real), pointer :: fcqm(:,:), fqqm(:)
            integer(ommp_integer), pointer :: fzqm(:)
            type(c_ptr) :: c_prt
            
            call c_f_pointer(cqm, fcqm, [3,n])
            call c_f_pointer(qqm, fqqm, [n])
            call c_f_pointer(zqm, fzqm, [n])
            call ommp_init_qm_helper(s, n, fcqm, fqqm, fzqm)
            c_prt = c_loc(s)
        end function
        
        subroutine C_ommp_qm_helper_set_frozen_atoms(s_prt, n, frozen) &
                bind(c, name='ommp_qm_helper_set_frozen_atoms')
            implicit none
            type(c_ptr), value :: s_prt, frozen
            integer(ommp_integer), value :: n
            type(ommp_qm_helper), pointer :: s
            integer(ommp_integer), pointer :: f(:)
            call c_f_pointer(s_prt, s)
            call c_f_pointer(frozen, f, [n])

            f = f + 1
            call ommp_qm_helper_set_frozen_atoms(s, n, f)
            f = f - 1
        end subroutine

        subroutine C_ommp_terminate_qm_helper(s_ptr) &
                bind(c, name='ommp_terminate_qm_helper')
            
            use mod_qm_helper, only: qm_helper_terminate, ommp_qm_helper
            
            implicit none

            type(c_ptr), value :: s_ptr
            type(ommp_qm_helper), pointer :: s
            
            call c_f_pointer(s_ptr, s)
            call qm_helper_terminate(s)
            deallocate(s)
        end subroutine
        
        subroutine C_ommp_qm_helper_update_coord(s_ptr, cqm) &
                bind(c, name='ommp_qm_helper_update_coord')
            implicit none

            type(c_ptr), value :: s_ptr
            type(c_ptr), value, intent(in) :: cqm
            
            real(ommp_real), pointer :: fcqm(:,:)
            type(ommp_qm_helper), pointer :: s
            
            call c_f_pointer(s_ptr, s)
            call c_f_pointer(cqm, fcqm, [3,s%qm_top%mm_atoms])
            call ommp_qm_helper_update_coord(s, fcqm)
        end subroutine

        subroutine C_ommp_qm_helper_set_attype(pqm, pattype) &
                 bind(c, name='ommp_qm_helper_set_attype')
            implicit none

            type(c_ptr), value, intent(in) :: pqm, pattype
            
            type(ommp_qm_helper), pointer :: qm
            integer(ommp_integer), pointer :: attype(:)

            call c_f_pointer(pqm, qm)
            call c_f_pointer(pattype, attype, [qm%qm_top%mm_atoms])

            call ommp_qm_helper_set_attype(qm, attype)
        end subroutine
        
        subroutine C_ommp_qm_helper_init_vdw_prm(pqm, cprmfile) &
                 bind(c, name='ommp_qm_helper_init_vdw_prm')
            implicit none

            type(c_ptr), value, intent(in) :: pqm
            character(kind=c_char), intent(in) :: cprmfile(OMMP_STR_CHAR_MAX)
            
            type(ommp_qm_helper), pointer :: qm
            character(len=OMMP_STR_CHAR_MAX) :: prmfile

            call c_f_pointer(pqm, qm)
            call c2f_string(cprmfile, prmfile)

            call ommp_qm_helper_init_vdw_prm(qm, prmfile)
        end subroutine
        
        subroutine C_ommp_qm_helper_init_vdw(pqm, peps, prad, pfac, &
                                             cvdw_type, cradius_rule, &
                                             cradius_size, cradius_type, &
                                             ceps_rule) &
                 bind(c, name='ommp_qm_helper_init_vdw')
            implicit none

            type(c_ptr), value, intent(in) :: pqm, peps, prad, pfac
            character(kind=c_char), intent(in) :: cvdw_type(OMMP_STR_CHAR_MAX), &
                                                  cradius_rule(OMMP_STR_CHAR_MAX), &
                                                  cradius_size(OMMP_STR_CHAR_MAX), &
                                                  cradius_type(OMMP_STR_CHAR_MAX), &
                                                  ceps_rule(OMMP_STR_CHAR_MAX)
            
            type(ommp_qm_helper), pointer :: qm
            character(len=OMMP_STR_CHAR_MAX) :: vdw_type, radius_rule, &
                                                radius_size, radius_type, &
                                                eps_rule
            real(ommp_real), pointer :: eps(:), rad(:), fac(:)

            call c_f_pointer(pqm, qm)
            call c_f_pointer(peps, eps, [qm%qm_top%mm_atoms])
            call c_f_pointer(prad, rad, [qm%qm_top%mm_atoms])
            call c_f_pointer(pfac, fac, [qm%qm_top%mm_atoms])
            call c2f_string(cvdw_type, vdw_type)
            call c2f_string(cradius_rule, radius_rule)
            call c2f_string(cradius_size, radius_size)
            call c2f_string(cradius_type, radius_type)
            call c2f_string(ceps_rule, eps_rule)

            call ommp_qm_helper_init_vdw(qm, eps, rad, fac, vdw_type, radius_rule, &
                                         radius_size, radius_type, eps_rule)
        end subroutine
        
        function C_ommp_qm_helper_vdw_energy(qm_prt, s_prt) &
                result(evdw) bind(c, name='ommp_qm_helper_vdw_energy')
            implicit none

            type(c_ptr), value :: qm_prt, s_prt
            type(ommp_system), pointer :: s
            type(ommp_qm_helper), pointer :: qm
            real(ommp_real) :: evdw

            call c_f_pointer(s_prt, s)
            call c_f_pointer(qm_prt, qm)

            evdw = ommp_qm_helper_vdw_energy(qm, s)
        end function
        
        subroutine C_ommp_qm_helper_vdw_geomgrad(qm_prt, s_prt, qmg_prt, mmg_prt) &
                bind(c, name='ommp_qm_helper_vdw_geomgrad')
            implicit none

            type(c_ptr), value :: qm_prt, s_prt, qmg_prt, mmg_prt
            type(ommp_system), pointer :: s
            type(ommp_qm_helper), pointer :: qm
            real(ommp_real), pointer :: qmg(:,:), mmg(:,:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(qm_prt, qm)
            call c_f_pointer(qmg_prt, qmg, [3,qm%qm_top%mm_atoms])
            call c_f_pointer(mmg_prt, mmg, [3,s%top%mm_atoms])

            call ommp_qm_helper_vdw_geomgrad(qm, s, qmg, mmg)
        end subroutine

        subroutine C_ommp_qm_helper_linkatom_geomgrad(qm_prt, s_prt, qmg_prt, mmg_prt, old_qmg_ptr) &
                bind(c, name='ommp_qm_helper_linkatom_geomgrad')
            implicit none

            type(c_ptr), value :: qm_prt, s_prt, qmg_prt, mmg_prt, old_qmg_ptr
            type(ommp_system), pointer :: s
            type(ommp_qm_helper), pointer :: qm
            real(ommp_real), pointer :: qmg(:,:), mmg(:,:), old_qmg(:,:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(qm_prt, qm)
            call c_f_pointer(qmg_prt, qmg, [3,qm%qm_top%mm_atoms])
            call c_f_pointer(old_qmg_ptr, old_qmg, [3,qm%qm_top%mm_atoms])
            call c_f_pointer(mmg_prt, mmg, [3,s%top%mm_atoms])

            call ommp_qm_helper_linkatom_geomgrad(qm, s, qmg, mmg, old_qmg)
        end subroutine

        subroutine C_ommp_prepare_qm_ele_ene(s_ptr, qm_ptr) &
                bind(c, name='ommp_prepare_qm_ele_ene')
            implicit none

            type(c_ptr), value :: s_ptr
            !! C pointer to system object
            type(c_ptr), value :: qm_ptr
            !! C pointer to qm_helper object

            type(ommp_qm_helper), pointer :: qm_help
            type(ommp_system), pointer :: s

            call c_f_pointer(qm_ptr, qm_help)
            call c_f_pointer(s_ptr, s)

            call ommp_prepare_qm_ele_ene(s, qm_help)
        end subroutine 
        
        subroutine C_ommp_prepare_qm_ele_grd(s_ptr, qm_ptr) &
                bind(c, name='ommp_prepare_qm_ele_grd')
            implicit none

            type(c_ptr), value :: s_ptr
            !! C pointer to system object
            type(c_ptr), value :: qm_ptr
            !! C pointer to qm_helper object

            type(ommp_qm_helper), pointer :: qm_help
            type(ommp_system), pointer :: s

            call c_f_pointer(qm_ptr, qm_help)
            call c_f_pointer(s_ptr, s)

            call ommp_prepare_qm_ele_grd(s, qm_help)
        end subroutine 

        function C_ommp_qm_helper_get_npol(qm_ptr) &
                result(npol) bind(C, name='ommp_qm_helper_get_npol')
            use mod_qm_helper, only: ommp_qm_helper
            implicit none

            type(c_ptr), value :: qm_ptr
            !! C pointer to qm_helper object

            type(ommp_qm_helper), pointer :: qm_help
            integer(ommp_integer) :: npol

            call c_f_pointer(qm_ptr, qm_help)
            if(qm_help%E_n2p_done) then
                npol = size(qm_help%E_n2p, 2, ommp_integer)
            else
                npol = 0
            end if
        end function
        
        function C_ommp_qm_helper_get_nmm(qm_ptr) &
                result(nmm) bind(C, name='ommp_qm_helper_get_nmm')
            use mod_qm_helper, only: ommp_qm_helper
            implicit none

            type(c_ptr), value :: qm_ptr
            !! C pointer to qm_helper object

            type(ommp_qm_helper), pointer :: qm_help
            integer(ommp_integer) :: nmm

            call c_f_pointer(qm_ptr, qm_help)
            if(qm_help%E_n2m_done) then
                nmm = size(qm_help%E_n2m, 2, ommp_integer)
            else
                nmm = 0
            end if
        end function

        function C_ommp_qm_helper_get_cqm(qm_ptr) &
                result(ptr) bind(C, name='ommp_qm_helper_get_cqm')
            use mod_qm_helper, only: ommp_qm_helper
            implicit none

            type(c_ptr), value :: qm_ptr
            !! C pointer to qm_helper object

            type(ommp_qm_helper), pointer :: qm_help
            type(c_ptr) :: ptr
           
            call c_f_pointer(qm_ptr, qm_help)
            ptr = c_loc(qm_help%qm_top%cmm)
        end function
        
        function C_ommp_qm_helper_get_E_n2p(qm_ptr) &
                result(ptr) bind(C, name='ommp_qm_helper_get_E_n2p')
            use mod_qm_helper, only: ommp_qm_helper
            implicit none

            type(c_ptr), value :: qm_ptr
            !! C pointer to qm_helper object

            type(ommp_qm_helper), pointer :: qm_help
            type(c_ptr) :: ptr
            
            call c_f_pointer(qm_ptr, qm_help)
            if(qm_help%E_n2p_done) then
                ptr = c_loc(qm_help%E_n2p)
            else
                ptr = c_null_ptr
            end if
        end function
        
        function C_ommp_qm_helper_get_G_n2p(qm_ptr) &
                result(ptr) bind(C, name='ommp_qm_helper_get_G_n2p')
            use mod_qm_helper, only: ommp_qm_helper
            implicit none

            type(c_ptr), value :: qm_ptr
            !! C pointer to qm_helper object

            type(ommp_qm_helper), pointer :: qm_help
            type(c_ptr) :: ptr
            
            call c_f_pointer(qm_ptr, qm_help)
            if(qm_help%G_n2p_done) then
                ptr = c_loc(qm_help%G_n2p)
            else
                ptr = c_null_ptr
            end if
        end function

        function C_ommp_qm_helper_get_E_n2m(qm_ptr) &
                result(ptr) bind(C, name='ommp_qm_helper_get_E_n2m')
            use mod_qm_helper, only: ommp_qm_helper
            implicit none

            type(c_ptr), value :: qm_ptr
            !! C pointer to qm_helper object

            type(ommp_qm_helper), pointer :: qm_help
            type(c_ptr) :: ptr
            
            call c_f_pointer(qm_ptr, qm_help)
            if(qm_help%E_n2m_done) then
                ptr = c_loc(qm_help%E_n2m)
            else
                ptr = c_null_ptr
            end if
        end function
        
        function C_ommp_qm_helper_get_G_n2m(qm_ptr) &
                result(ptr) bind(C, name='ommp_qm_helper_get_G_n2m')
            use mod_qm_helper, only: ommp_qm_helper
            implicit none

            type(c_ptr), value :: qm_ptr
            !! C pointer to qm_helper object

            type(ommp_qm_helper), pointer :: qm_help
            type(c_ptr) :: ptr
            
            call c_f_pointer(qm_ptr, qm_help)
            if(qm_help%G_n2m_done) then
                ptr = c_loc(qm_help%G_n2m)
            else
                ptr = c_null_ptr
            end if
        end function

        function C_ommp_qm_helper_get_H_n2m(qm_ptr) &
                result(ptr) bind(C, name='ommp_qm_helper_get_H_n2m')
            use mod_qm_helper, only: ommp_qm_helper
            implicit none

            type(c_ptr), value :: qm_ptr
            !! C pointer to qm_helper object

            type(ommp_qm_helper), pointer :: qm_help
            type(c_ptr) :: ptr
            
            call c_f_pointer(qm_ptr, qm_help)
            if(qm_help%H_n2m_done) then
                ptr = c_loc(qm_help%H_n2m)
            else
                ptr = c_null_ptr
            end if
        end function
        
        function C_ommp_qm_helper_get_E_m2n(qm_ptr) &
                result(ptr) bind(C, name='ommp_qm_helper_get_E_m2n')
            use mod_qm_helper, only: ommp_qm_helper
            implicit none

            type(c_ptr), value :: qm_ptr
            !! C pointer to qm_helper object

            type(ommp_qm_helper), pointer :: qm_help
            type(c_ptr) :: ptr
            
            call c_f_pointer(qm_ptr, qm_help)
            if(qm_help%E_m2n_done) then
                ptr = c_loc(qm_help%E_m2n)
            else
                ptr = c_null_ptr
            end if
        end function
        
        function C_ommp_qm_helper_get_E_p2n(qm_ptr) &
                result(ptr) bind(C, name='ommp_qm_helper_get_E_p2n')
            use mod_qm_helper, only: ommp_qm_helper
            implicit none

            type(c_ptr), value :: qm_ptr
            !! C pointer to qm_helper object

            type(ommp_qm_helper), pointer :: qm_help
            type(c_ptr) :: ptr
            
            call c_f_pointer(qm_ptr, qm_help)
            if(qm_help%E_p2n_done) then
                ptr = c_loc(qm_help%E_p2n)
            else
                ptr = c_null_ptr
            end if
        end function
        
        function C_ommp_qm_helper_get_V_m2n(qm_ptr) &
                result(ptr) bind(C, name='ommp_qm_helper_get_V_m2n')
            use mod_qm_helper, only: ommp_qm_helper
            implicit none

            type(c_ptr), value :: qm_ptr
            !! C pointer to qm_helper object

            type(ommp_qm_helper), pointer :: qm_help
            type(c_ptr) :: ptr
            
            call c_f_pointer(qm_ptr, qm_help)
            if(qm_help%V_m2n_done) then
                ptr = c_loc(qm_help%V_m2n)
            else
                ptr = c_null_ptr
            end if
        end function
       
        function C_ommp_qm_helper_get_V_p2n(qm_ptr) &
                result(ptr) bind(C, name='ommp_qm_helper_get_V_p2n')
            use mod_qm_helper, only: ommp_qm_helper
            implicit none

            type(c_ptr), value :: qm_ptr
            !! C pointer to qm_helper object

            type(ommp_qm_helper), pointer :: qm_help
            type(c_ptr) :: ptr
            
            call c_f_pointer(qm_ptr, qm_help)
            if(qm_help%V_p2n_done) then
                ptr = c_loc(qm_help%V_p2n)
            else
                ptr = c_null_ptr
            end if
        end function
       
        function C_ommp_qm_helper_get_qm_atoms(qm_ptr) &
                result(n) bind(C, name='ommp_qm_helper_get_qm_atoms')
            use mod_qm_helper, only: ommp_qm_helper
            implicit none

            type(c_ptr), value :: qm_ptr
            !! C pointer to qm_helper object

            type(ommp_qm_helper), pointer :: qm_help
            integer(ommp_integer) :: n
            
            call c_f_pointer(qm_ptr, qm_help)
            n = qm_help%qm_top%mm_atoms
        end function

        function C_ommp_qm_helper_use_nonbonded(qm_ptr) &
                result(u) bind(C, name='ommp_qm_helper_use_nonbonded')
            use mod_qm_helper, only: ommp_qm_helper
            implicit none

            type(c_ptr), value :: qm_ptr
            !! C pointer to qm_helper object
            type(ommp_qm_helper), pointer :: qm_help
            logical(c_bool) :: u

            call c_f_pointer(qm_ptr, qm_help)
            u = qm_help%use_nonbonded
        end function
        
        function C_ommp_qm_helper_use_frozen(s_prt) bind(c, name='ommp_qm_helper_use_frozen')
            use mod_qm_helper, only: ommp_qm_helper
            implicit none

            type(c_ptr), value :: s_prt
            type(ommp_qm_helper), pointer :: s
            logical(c_bool) :: C_ommp_qm_helper_use_frozen

            call c_f_pointer(s_prt, s)
            C_ommp_qm_helper_use_frozen = s%qm_top%use_frozen
        end function C_ommp_qm_helper_use_frozen

        function C_ommp_qm_helper_get_frozen(s_prt) bind(c, name='ommp_qm_helper_get_frozen')
            use mod_qm_helper, only: ommp_qm_helper
            implicit none
            
            type(c_ptr), value :: s_prt
            type(ommp_qm_helper), pointer :: s
            type(c_ptr) :: C_ommp_qm_helper_get_frozen

            call c_f_pointer(s_prt, s)
            C_ommp_qm_helper_get_frozen = c_loc(s%qm_top%frozen)
        end function C_ommp_qm_helper_get_frozen

        function C_ommp_create_link_atom(qm_prt, s_prt, imm, iqm, ila, prmfile, &
                                         ladist, neel_remove) &
                result(la_idx) &
                bind(c, name='ommp_create_link_atom')
            implicit none
            
            type(c_ptr), value :: s_prt, qm_prt
            integer(ommp_integer), value :: iqm, imm, ila, neel_remove
            real(ommp_real), value :: ladist
            character(kind=c_char), intent(in) :: prmfile(OMMP_STR_CHAR_MAX)
            

            type(ommp_system), pointer :: s
            type(ommp_qm_helper), pointer :: qm

            integer(ommp_integer) :: la_idx
            character(len=OMMP_STR_CHAR_MAX) :: prm_file
            
            call c2f_string(prmfile, prm_file)
            call c_f_pointer(s_prt, s)
            call c_f_pointer(qm_prt, qm)

            la_idx = ommp_create_link_atom(qm, s, imm, iqm, ila, prm_file, &
                                           ladist, neel_remove)
        end function

        subroutine C_ommp_get_link_atom_coordinates(s_p, la_idx, crd_p) &
                bind(c, name="ommp_get_link_atom_coordinates")
            implicit none

            type(c_ptr), value :: s_p, crd_p
            integer(ommp_integer), value :: la_idx

            type(ommp_system), pointer :: s
            real(ommp_real), dimension(:), pointer :: crd

            call c_f_pointer(s_p, s)
            call c_f_pointer(crd_p, crd, [3])

            call ommp_get_link_atom_coordinates(s, la_idx, crd)

        end subroutine 
        
        subroutine C_ommp_update_link_atoms_position(qm_prt, s_prt) &
                bind(c, name='ommp_update_link_atoms_position')
            implicit none
            
            type(c_ptr), value :: s_prt, qm_prt

            type(ommp_system), pointer :: s
            type(ommp_qm_helper), pointer :: qm

            call c_f_pointer(s_prt, s)
            call c_f_pointer(qm_prt, qm)
            
            call ommp_update_link_atoms_position(s, qm)
        end subroutine
end module mod_ommp_C_interface
