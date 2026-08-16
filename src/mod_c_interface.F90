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

        function C_ommp_yst_get_n(yst_p) result(r) bind(c, name='ommp_yst_get_n')
            use mod_adjacency_mat, only : yale_sparse

            implicit none

            type(c_ptr), value :: yst_p
            type(yale_sparse), pointer :: yst
            integer(ommp_integer) :: r

            call c_f_pointer(yst_p, yst)
            r = yst%n
        end function

        function C_ommp_yst_get_ri(yst_p) result(r) bind(c, name='ommp_yst_get_ri')
            use mod_adjacency_mat, only : yale_sparse

            implicit none

            type(c_ptr), value :: yst_p
            type(yale_sparse), pointer :: yst
            type(c_ptr) :: r

            call c_f_pointer(yst_p, yst)
            r = c_loc(yst%ri)
        end function

        function C_ommp_yst_get_ci(yst_p) result(r) bind(c, name='ommp_yst_get_ci')
            use mod_adjacency_mat, only : yale_sparse

            implicit none

            type(c_ptr), value :: yst_p
            type(yale_sparse), pointer :: yst
            type(c_ptr) :: r

            call c_f_pointer(yst_p, yst)
            r = c_loc(yst%ci)
        end function

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
        
        subroutine C_ommp_set_default_solver(s_prt, solver) bind(c, name='ommp_set_default_solver')
            implicit none 

            integer(ommp_integer), intent(in), value :: solver
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
           
            call c_f_pointer(s_prt, s)
            
            call ommp_set_default_solver(s, solver)
        end subroutine C_ommp_set_default_solver
        
        subroutine C_ommp_set_default_matv(s_prt, matv) bind(c, name='ommp_set_default_matv')
            implicit none 

            integer(ommp_integer), intent(in), value :: matv
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
           
            call c_f_pointer(s_prt, s)
            
            call ommp_set_default_matv(s, matv)
        end subroutine C_ommp_set_default_matv

        subroutine C_ommp_set_polarization_conv_thr(s_prt, conv_thr) &
                bind(c, name='ommp_set_polarization_conv_thr')
            implicit none

            real(c_double), value :: conv_thr
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
           
            call c_f_pointer(s_prt, s)
            
            call ommp_set_polarization_conv_thr(s, conv_thr)
        end subroutine C_ommp_set_polarization_conv_thr

        subroutine C_ommp_set_polarization_use_guess(s_prt, use_guess) &
                bind(c, name='ommp_set_polarization_use_guess')
            implicit none

            logical(c_bool), value :: use_guess
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            
            call c_f_pointer(s_prt, s)
            
            call ommp_set_polarization_use_guess(s, use_guess)
        end subroutine C_ommp_set_polarization_use_guess

        subroutine C_ommp_fatal(c_msg) &
                bind(c, name='ommp_fatal')
            implicit none
            
            character(kind=c_char), intent(in) :: c_msg(OMMP_STR_CHAR_MAX)
            character(len=OMMP_STR_CHAR_MAX) :: msg
            
            call c2f_string(c_msg, msg)
            call ommp_fatal(msg)

        end subroutine C_ommp_fatal
        
        subroutine C_ommp_message(c_msg, level, c_pre) &
                bind(c, name='ommp_message')
            implicit none
            
            character(kind=c_char), intent(in) :: c_msg(OMMP_STR_CHAR_MAX)
            integer(kind=ommp_integer), value :: level
            character(kind=c_char), intent(in) :: c_pre(OMMP_STR_CHAR_MAX)
            character(len=OMMP_STR_CHAR_MAX) :: msg
            character(len=OMMP_STR_CHAR_MAX) :: pre
            
            call c2f_string(c_msg, msg)
            call c2f_string(c_pre, pre)
            call ommp_message(msg, level, pre)

        end subroutine C_ommp_message
        
        subroutine C_ommp_time_pull(c_msg) &
                bind(c, name='ommp_time_pull')
            implicit none
            
            character(kind=c_char), intent(in) :: c_msg(OMMP_STR_CHAR_MAX)
            character(len=OMMP_STR_CHAR_MAX) :: msg
            
            call c2f_string(c_msg, msg)
            call ommp_time_pull(trim(msg))

        end subroutine C_ommp_time_pull
        
        subroutine C_ommp_time_push() &
                bind(c, name='ommp_time_push')
            implicit none
            
            call ommp_time_push

        end subroutine C_ommp_time_push
        
        subroutine C_ommp_set_outputfile(fname) &
                bind(c, name='ommp_set_outputfile')
            implicit none
            
            character(kind=c_char), intent(in) :: fname(OMMP_STR_CHAR_MAX)
            character(len=OMMP_STR_CHAR_MAX) :: ffname
            
            call c2f_string(fname, ffname)
            call ommp_set_outputfile(ffname)

        end subroutine C_ommp_set_outputfile
        
        subroutine C_ommp_close_outputfile() &
                bind(c, name='ommp_close_outputfile')
            implicit none
            
            call ommp_close_outputfile()

        end subroutine C_ommp_close_outputfile
        
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
            call c_f_pointer(new_c_p, new_c, [3_ommp_integer, s%top%mm_atoms])
            
            call update_coordinates(s, new_c)
        end subroutine
        
        ! Functions mapped on actual Fortran interface functions
        ! OMMP System Object housekeeping
        function C_ommp_init_mmp(filename) &
                result(c_prt) bind(c, name='ommp_init_mmp')
           !! Initalize OMMP System Object from .mmp file  
            implicit none

            type(ommp_system), pointer, save :: s
            character(kind=c_char), intent(in) :: filename(OMMP_STR_CHAR_MAX)
            character(len=OMMP_STR_CHAR_MAX) :: input_file
            type(c_ptr) :: c_prt
            
            !allocate(s)
            
            call c2f_string(filename, input_file)
            call ommp_init_mmp(s, input_file)
            c_prt = c_loc(s)
        end function
        
        function C_ommp_init_xyz(xyzfile, prmfile) &
                result(c_prt) bind(c, name='ommp_init_xyz')
            !! Initialize the library using a Tinker xyz and a Tinker prm
            implicit none
            
            type(ommp_system), pointer, save :: s
            character(kind=c_char), intent(in) :: xyzfile(OMMP_STR_CHAR_MAX), & 
                                                  prmfile(OMMP_STR_CHAR_MAX)
            character(len=OMMP_STR_CHAR_MAX) :: xyz_file, prm_file
            type(c_ptr) :: c_prt

            !allocate(s)
            
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
            call ommp_set_frozen_atoms(s, n, f)
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
            call ommp_turn_pol_off(s, n, f)
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
        subroutine C_ommp_set_external_field(s_prt, ext_field_prt, solver, matv) &
                bind(c, name='ommp_set_external_field')
            implicit none
            
            type(c_ptr), value :: s_prt
            type(c_ptr), value :: ext_field_prt
            integer(ommp_integer), intent(in), value :: solver
            integer(ommp_integer), intent(in), value :: matv
            
            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: ext_field(:,:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(ext_field_prt, ext_field, [3_ommp_integer, s%eel%pol_atoms])

            call ommp_set_external_field(s, ext_field, solver, matv, .true.)
        end subroutine C_ommp_set_external_field
        
        subroutine C_ommp_set_external_field_nomm(s_prt, ext_field_prt, solver, matv) &
                bind(c, name='ommp_set_external_field_nomm')
            !!use mod_mmpol, only: pol_atoms
            
            implicit none
            
            type(c_ptr), value :: s_prt
            type(c_ptr), value :: ext_field_prt
            integer(ommp_integer), intent(in), value :: solver
            integer(ommp_integer), intent(in), value :: matv
            
            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: ext_field(:,:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(ext_field_prt, ext_field, [3_ommp_integer, s%eel%pol_atoms])
            
            call ommp_set_external_field(s, ext_field, solver, matv, .false.)
        end subroutine C_ommp_set_external_field_nomm

        subroutine C_ommp_set_fit_potential(s_prt, fit_pot_prt, n_pts) &
                bind(c, name='ommp_set_fit_potential')
            !! Set the electrostatic potential at the density fitting points.
            implicit none
            
            type(c_ptr), value :: s_prt
            type(c_ptr), value :: fit_pot_prt
            integer(ommp_integer), intent(in), value :: n_pts

            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: fit_pot(:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(fit_pot_prt, fit_pot, [n_pts])

            call ommp_set_fit_potential(s, fit_pot)
        end subroutine C_ommp_set_fit_potential

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
            call c_f_pointer(cext, fcext, [3_ommp_integer,n])
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
            call c_f_pointer(cext, fcext, [3_ommp_integer,n])
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
            call c_f_pointer(cext, fcext, [3_ommp_integer,n])
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
            call c_f_pointer(grd_prt, grd, [3_ommp_integer, s%top%mm_atoms])
            
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
            call c_f_pointer(grd_prt, grd, [3_ommp_integer, s%top%mm_atoms])
            call c_f_pointer(pE, E, [3_ommp_integer, s%top%mm_atoms])
            call c_f_pointer(pE_grd, Egrd, [6_ommp_integer, s%top%mm_atoms])
            
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
            call c_f_pointer(grd_prt, grd, [3_ommp_integer, s%top%mm_atoms])
            
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
            call c_f_pointer(grd_prt, grd, [3_ommp_integer, s%top%mm_atoms])
            
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
            call c_f_pointer(grd_prt, grd, [3_ommp_integer, s%top%mm_atoms])
            
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
            call c_f_pointer(grd_prt, grd, [3_ommp_integer, s%top%mm_atoms])
            
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
            call c_f_pointer(grd_prt, grd, [3_ommp_integer, s%top%mm_atoms])
            
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
            call c_f_pointer(grd_prt, grd, [3_ommp_integer, s%top%mm_atoms])
            
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
            call c_f_pointer(grd_prt, grd, [3_ommp_integer, s%top%mm_atoms])
            
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
            call c_f_pointer(grd_prt, grd, [3_ommp_integer, s%top%mm_atoms])
            
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
            call c_f_pointer(grd_prt, grd, [3_ommp_integer, s%top%mm_atoms])
            
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
            call c_f_pointer(grd_prt, grd, [3_ommp_integer, s%top%mm_atoms])
            
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
            call c_f_pointer(grd_prt, grd, [3_ommp_integer, s%top%mm_atoms])
            
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
            call c_f_pointer(grd_prt, grd, [3_ommp_integer, s%top%mm_atoms])
            
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
            call c_f_pointer(grd_prt, grd, [3_ommp_integer, s%top%mm_atoms])
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
            call c_f_pointer(grd_prt, grd, [3_ommp_integer, s%top%mm_atoms])
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
            call c_f_pointer(grd_prt, grd, [3_ommp_integer, s%top%mm_atoms])
            
            call ommp_full_geomgrad(s, grd)
        end subroutine

#ifdef WITH_HDF5
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
            call save_system_as_hdf5(hdf5out, s, err, trim(nms), logical(.false., kind=ommp_logical))
            
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
            call save_system_as_hdf5(hdf5out, s, err, trim(nms), logical(.true., kind=ommp_logical))
            
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
        
        function C_ommp_get_attypemm(s_prt) bind(c, name='ommp_get_attypemm')
            !! Return the c-pointer to the array containing the coordinates of
            !! MM atoms.
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            type(c_ptr) :: C_ommp_get_attypemm

            call c_f_pointer(s_prt, s)
            if(s%top%attype_initialized) then
                C_ommp_get_attypemm = c_loc(s%top%attype)
            else
                C_ommp_get_attypemm = c_null_ptr
            end if
        end function C_ommp_get_attypemm


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
                               int(size(s%eel%polar_mm), ommp_integer), s%eel%C_polar_mm)
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
            call c_f_pointer(cext, fcext, [3_ommp_integer,n])
            call c_f_pointer(E, fE, [3_ommp_integer,n])
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
            call c_f_pointer(cext, fcext, [3_ommp_integer,n])
            call c_f_pointer(E, fE, [3_ommp_integer,n])
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
            call c_f_pointer(cext, fcext, [3_ommp_integer,n])
            call c_f_pointer(E, fE, [3_ommp_integer,n])
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
            
            call c_f_pointer(cqm, fcqm, [3_ommp_integer,n])
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

            call ommp_qm_helper_set_frozen_atoms(s, n, f)
        end subroutine

        subroutine C_ommp_terminate_qm_helper(s_ptr) &
                bind(c, name='ommp_terminate_qm_helper')
            
            use mod_qm_helper, only: qm_helper_terminate, ommp_qm_helper
            
            implicit none

            type(c_ptr), value :: s_ptr
            type(ommp_qm_helper), pointer :: s
            
            if(c_associated(s_ptr)) then
                call c_f_pointer(s_ptr, s)
                call qm_helper_terminate(s)
                deallocate(s)
            end if
        end subroutine
        
        subroutine C_ommp_qm_helper_update_coord(s_ptr, cqm) &
                bind(c, name='ommp_qm_helper_update_coord')
            implicit none

            type(c_ptr), value :: s_ptr
            type(c_ptr), value, intent(in) :: cqm
            
            real(ommp_real), pointer :: fcqm(:,:)
            type(ommp_qm_helper), pointer :: s
            
            call c_f_pointer(s_ptr, s)
            call c_f_pointer(cqm, fcqm, [3_ommp_integer,s%qm_top%mm_atoms])
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

        subroutine C_ommp_qm_helper_vdw_energy_by_atoms(qm_prt, s_prt, evdw_ba_ptr) &
                bind(c, name='ommp_qm_helper_vdw_energy_by_atoms')
            implicit none

            type(c_ptr), value :: qm_prt, s_prt, evdw_ba_ptr
            type(ommp_system), pointer :: s
            type(ommp_qm_helper), pointer :: qm
            real(ommp_real), pointer :: evdw_ba(:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(qm_prt, qm)
            call c_f_pointer(evdw_ba_ptr, evdw_ba, [qm%qm_top%mm_atoms])

            call ommp_qm_helper_vdw_energy_by_atom(qm, s, evdw_ba)
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
            call c_f_pointer(qmg_prt, qmg, [3_ommp_integer,qm%qm_top%mm_atoms])
            call c_f_pointer(mmg_prt, mmg, [3_ommp_integer,s%top%mm_atoms])

            call ommp_qm_helper_vdw_geomgrad(qm, s, qmg, mmg)
        end subroutine
        
        subroutine C_ommp_qm_helper_link_atom_geomgrad(qm_prt, s_prt, qmg_prt, mmg_prt, old_qmg_ptr) &
                bind(c, name='ommp_qm_helper_link_atom_geomgrad')
            implicit none

            type(c_ptr), value :: qm_prt, s_prt, qmg_prt, mmg_prt, old_qmg_ptr
            type(ommp_system), pointer :: s
            type(ommp_qm_helper), pointer :: qm
            real(ommp_real), pointer :: qmg(:,:), mmg(:,:), old_qmg(:,:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(qm_prt, qm)
            call c_f_pointer(qmg_prt, qmg, [3_ommp_integer,qm%qm_top%mm_atoms])
            call c_f_pointer(old_qmg_ptr, old_qmg, [3_ommp_integer,qm%qm_top%mm_atoms])
            call c_f_pointer(mmg_prt, mmg, [3_ommp_integer,s%top%mm_atoms])

            call ommp_qm_helper_link_atom_geomgrad(qm, s, qmg, mmg, old_qmg)
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
            
            call ommp_update_link_atoms_position(qm, s)
        end subroutine

        function C_ommp_system_from_qm_helper(cqmh, cprm_file) &
                result(csys) bind(c, name='ommp_system_from_qm_helper')
            implicit none
            
            type(c_ptr), value, intent(in) :: cqmh

            type(ommp_system), pointer, save :: s
            type(ommp_qm_helper), pointer :: qm
            type(c_ptr) :: csys

            character(kind=c_char), intent(in) :: cprm_file(OMMP_STR_CHAR_MAX)
            character(len=OMMP_STR_CHAR_MAX) :: prm_file
            
            call c2f_string(cprm_file, prm_file)
            call c_f_pointer(cqmh, qm)
            
            !allocate(s)

            call ommp_system_from_qm_helper(qm, prm_file, s)
            csys = c_loc(s)
        end function
    
        subroutine C_ommp_set_vdw_cutoff(sp, cutoff) &
                bind(c, name='ommp_set_vdw_cutoff')

            implicit none

             type(c_ptr), value, intent(in) :: sp
            real(ommp_real), intent(in), value :: cutoff
           
            type(ommp_system), pointer :: s
            
            call c_f_pointer(sp, s)
            call ommp_set_vdw_cutoff(s, cutoff)
        end subroutine
        
        subroutine C_ommp_set_fmm_lmax_pol(sp, l) &
                bind(c, name='ommp_set_fmm_lmax_pol')

            implicit none

            type(c_ptr), value, intent(in) :: sp
            integer(ommp_integer), intent(in), value :: l
           
            type(ommp_system), pointer :: s
            
            call c_f_pointer(sp, s)
            s%eel%fmm_maxl_pol = l
        end subroutine
        
        subroutine C_ommp_set_fmm_lmax(sp, l) &
                bind(c, name='ommp_set_fmm_lmax')

            implicit none

            type(c_ptr), value, intent(in) :: sp
            integer(ommp_integer), intent(in), value :: l
           
            type(ommp_system), pointer :: s
            
            call c_f_pointer(sp, s)
            s%eel%fmm_maxl_static = l
        end subroutine

        subroutine C_ommp_set_fmm_cache_mode(sp, m) &
                bind(c, name='ommp_set_fmm_cache_mode')
            !! Controls the M2L/L2L/M2M/P2M/L2P rotation-matrix caching
            !! used across CG iterations (see fmm_cg_init in
            !! mod_electrostatics, fmm_rotcache_try_enable in mod_fmm):
            !! 0 (default) = cache M2L+L2P if memory allows, silently fall
            !! back to no caching otherwise; 1 = caching disabled; 2 =
            !! force-cache all five terms, erroring out if memory doesn't
            !! allow it. Just sets the field -- consulted lazily the next
            !! time a CG/DIIS polarization solve starts, no immediate
            !! side effect.

            implicit none

            type(c_ptr), value, intent(in) :: sp
            integer(ommp_integer), intent(in), value :: m

            type(ommp_system), pointer :: s

            call c_f_pointer(sp, s)
            s%eel%fmm_cache_mode = m
        end subroutine

        subroutine C_ommp_set_fmm_distance(sp, d) &
                bind(c, name='ommp_set_fmm_distance')

            use mod_electrostatics, only: fmm_coordinates_update

            implicit none

            type(c_ptr), value, intent(in) :: sp
            real(ommp_real), intent(in), value :: d
           
            type(ommp_system), pointer :: s
            
            call c_f_pointer(sp, s)
            s%eel%fmm_distance = d
            call fmm_coordinates_update(s%eel)
        end subroutine
        
        subroutine C_ommp_set_fmm_min_cell_size(sp, d) &
                bind(c, name='ommp_set_fmm_min_cell_size')

            use mod_electrostatics, only: fmm_coordinates_update
            
            implicit none

            type(c_ptr), value, intent(in) :: sp
            real(ommp_real), intent(in), value :: d
           
            type(ommp_system), pointer :: s
            
            call c_f_pointer(sp, s)
            s%eel%fmm_min_cell_size = d
            call fmm_coordinates_update(s%eel)
        end subroutine

        subroutine C_ommp_set_fmm_params(sp, distance, min_cell_size) &
                bind(c, name='ommp_set_fmm_params')
            !! Sets fmm_distance AND fmm_min_cell_size together, triggering
            !! fmm_coordinates_update only ONCE after both are consistent.
            !! Setting them one at a time via ommp_set_fmm_distance/
            !! ommp_set_fmm_min_cell_size (each of which independently
            !! triggers its own tree rebuild) exposes an intermediate state
            !! where one of the two is still at its Fortran type-default of
            !! 0.0 -- a degenerate value that has been observed to crash
            !! free_tree for small systems (root-caused in the `hessian`
            !! branch FMM-caching work, 2026-08-15). Use this instead of
            !! the two individual setters when both values are being
            !! configured together (e.g. from a fresh smart-input read).

            use mod_electrostatics, only: fmm_coordinates_update

            implicit none

            type(c_ptr), value, intent(in) :: sp
            real(ommp_real), intent(in), value :: distance, min_cell_size

            type(ommp_system), pointer :: s

            call c_f_pointer(sp, s)
            s%eel%fmm_distance = distance
            s%eel%fmm_min_cell_size = min_cell_size
            call fmm_coordinates_update(s%eel)
        end subroutine
        
        function C_ommp_use_fmm(s_prt) bind(c, name='ommp_use_fmm')
            !! Return true if the current forcefield is AMOEBA, and false in
            !! all other cases.
            implicit none

            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            logical(c_bool) :: C_ommp_use_fmm

            call c_f_pointer(s_prt, s)
            C_ommp_use_fmm = s%eel%use_fmm
        end function C_ommp_use_fmm
        
        subroutine C_ommp_enable_fmm(sp) &
                bind(c, name='ommp_enable_fmm')

            implicit none

            type(c_ptr), value, intent(in) :: sp
           
            type(ommp_system), pointer :: s
            
            call c_f_pointer(sp, s)
            s%eel%use_fmm = .true.
        end subroutine
        
        subroutine C_ommp_disable_fmm(sp) &
                bind(c, name='ommp_disable_fmm')

            implicit none

            type(c_ptr), value, intent(in) :: sp
           
            type(ommp_system), pointer :: s
            
            call c_f_pointer(sp, s)
            s%eel%use_fmm = .false.
        end subroutine

        function C_ommp_use_density_fit(s_prt) bind(c, name='ommp_use_density_fit')
            !! Return true if density fitting is enabled, false otherwise.
            implicit none

            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            logical(c_bool) :: C_ommp_use_density_fit

            call c_f_pointer(s_prt, s)
            C_ommp_use_density_fit = s%use_density_fit
        end function C_ommp_use_density_fit

        subroutine C_ommp_ignore_duplicated_angle_prm() &
                bind(c, name='ommp_ignore_duplicated_angle_prm')

            implicit none

            call ommp_ignore_duplicated_angle_prm
        end subroutine
        
        subroutine C_ommp_ignore_duplicated_opb_prm() &
                bind(c, name='ommp_ignore_duplicated_opb_prm')

            implicit none

            call ommp_ignore_duplicated_opb_prm
        end subroutine

        subroutine C_ommp_init_density_fit(s_prt, qmh_prt, &
                                           charge_point_type, charge_n_pts_per_atom, charge_radius, &
                                           fit_point_type, fit_n_pts_per_atom, fit_radius, &
                                           charge_top_source, fit_top_source) &
                bind(c, name='ommp_init_density_fit')

            implicit none

            type(c_ptr), value, intent(in) :: s_prt
            type(c_ptr), value, intent(in) :: qmh_prt
            integer(ommp_integer), intent(in), value :: charge_point_type
            integer(ommp_integer), intent(in), value :: charge_n_pts_per_atom
            real(ommp_real), intent(in), value :: charge_radius
            integer(ommp_integer), intent(in), value :: fit_point_type
            integer(ommp_integer), intent(in), value :: fit_n_pts_per_atom
            real(ommp_real), intent(in), value :: fit_radius
            character(kind=c_char), intent(in) :: charge_top_source(OMMP_STR_CHAR_MAX)
            character(kind=c_char), intent(in) :: fit_top_source(OMMP_STR_CHAR_MAX)

            type(ommp_system), pointer :: s
            type(ommp_qm_helper), pointer :: qmh
            character(len=OMMP_STR_CHAR_MAX) :: charge_src, fit_src

            call c_f_pointer(s_prt, s)
            call c_f_pointer(qmh_prt, qmh)

            !! Convert C strings to Fortran strings
            call c2f_string(charge_top_source, charge_src)
            call c2f_string(fit_top_source, fit_src)

            !! Delegate topology selection and initialization to ommp_init_density_fit
            call ommp_init_density_fit(s, qmh, &
                                       charge_point_type, charge_n_pts_per_atom, charge_radius, &
                                       fit_point_type, fit_n_pts_per_atom, fit_radius, &
                                       charge_src, fit_src)

        end subroutine

        function C_ommp_get_df_n_pts(s_prt) bind(c, name='ommp_get_df_n_pts')
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            integer(ommp_integer) :: C_ommp_get_df_n_pts

            call c_f_pointer(s_prt, s)
            if(allocated(s%df)) then
                C_ommp_get_df_n_pts = s%df%n_pts
            else
                C_ommp_get_df_n_pts = 0
            end if
        end function C_ommp_get_df_n_pts

        function C_ommp_get_df_n_charges(s_prt) bind(c, name='ommp_get_df_n_charges')
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            integer(ommp_integer) :: C_ommp_get_df_n_charges

            call c_f_pointer(s_prt, s)
            if(allocated(s%df)) then
                C_ommp_get_df_n_charges = s%df%n_charges
            else
                C_ommp_get_df_n_charges = 0
            end if
        end function C_ommp_get_df_n_charges

        function C_ommp_get_df_n_qm_atoms(s_prt) bind(c, name='ommp_get_df_n_qm_atoms')
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            integer(ommp_integer) :: C_ommp_get_df_n_qm_atoms

            call c_f_pointer(s_prt, s)
            if(allocated(s%df)) then
                C_ommp_get_df_n_qm_atoms = s%df%qm_top%mm_atoms
            else
                C_ommp_get_df_n_qm_atoms = 0
            end if
        end function C_ommp_get_df_n_qm_atoms

        function C_ommp_get_df_initialized(s_prt) bind(c, name='ommp_get_df_initialized')
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            logical(c_bool) :: C_ommp_get_df_initialized

            call c_f_pointer(s_prt, s)
            if(allocated(s%df)) then
                C_ommp_get_df_initialized = s%df%initialized
            else
                C_ommp_get_df_initialized = .false.
            end if
        end function C_ommp_get_df_initialized

        function C_ommp_get_df_charge_coord(s_prt) bind(c, name='ommp_get_df_charge_coord')
            !! Return the c-pointer to the charge coordinates array
            !! (3 x n_charges). Null if not available.
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            type(c_ptr) :: C_ommp_get_df_charge_coord

            call c_f_pointer(s_prt, s)
            if(allocated(s%df)) then
                C_ommp_get_df_charge_coord = c_loc(s%df%charge_coord)
            else
                C_ommp_get_df_charge_coord = c_null_ptr
            end if
        end function C_ommp_get_df_charge_coord

        function C_ommp_get_df_fit_point_coord(s_prt) bind(c, name='ommp_get_df_fit_point_coord')
            !! Return the c-pointer to the fitting point coordinates
            !! array (3 x n_pts). Null if not available.
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            type(c_ptr) :: C_ommp_get_df_fit_point_coord

            call c_f_pointer(s_prt, s)
            if(allocated(s%df)) then
                C_ommp_get_df_fit_point_coord = c_loc(s%df%fit_point_coord)
            else
                C_ommp_get_df_fit_point_coord = c_null_ptr
            end if
        end function C_ommp_get_df_fit_point_coord

        function C_ommp_get_df_target_charges(s_prt) bind(c, name='ommp_get_df_target_charges')
            !! Return the c-pointer to the target charges array
            !! (n_charges). Null if not available.
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            type(c_ptr) :: C_ommp_get_df_target_charges

            call c_f_pointer(s_prt, s)
            if(allocated(s%df)) then
                C_ommp_get_df_target_charges = c_loc(s%df%target_charges)
            else
                C_ommp_get_df_target_charges = c_null_ptr
            end if
        end function C_ommp_get_df_target_charges

        function C_ommp_get_df_X(s_prt) bind(c, name='ommp_get_df_X')
            !! Return the c-pointer to the design matrix X
            !! (n_charges x n_pts). Null if not available.
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            type(c_ptr) :: C_ommp_get_df_X
            call c_f_pointer(s_prt, s)
            if(allocated(s%df)) then
                C_ommp_get_df_X = c_loc(s%df%X)
            else
                C_ommp_get_df_X = c_null_ptr
            end if
        end function C_ommp_get_df_X

        function C_ommp_get_df_Xinv(s_prt) bind(c, name='ommp_get_df_Xinv')
            !! Return the c-pointer to the inverse/pseudoinverse
            !! design matrix Xinv (n_pts x n_charges). Null if not
            !! available.
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            type(c_ptr) :: C_ommp_get_df_Xinv

            call c_f_pointer(s_prt, s)
            if(allocated(s%df)) then
                C_ommp_get_df_Xinv = c_loc(s%df%Xinv)
            else
                C_ommp_get_df_Xinv = c_null_ptr
            end if
        end function C_ommp_get_df_Xinv

        function C_ommp_get_df_VXI_m(s_prt) bind(c, name='ommp_get_df_VXI_m')
            !! Return the c-pointer to the projected static quantity
            !! VXI_m = V_m2q @ Xinv. Null if not available.
            !! Triggers computation via the Fortran interface.
            use mod_density_fit, only: df_project_static

            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            type(c_ptr) :: C_ommp_get_df_VXI_m

            call c_f_pointer(s_prt, s)
            call df_project_static(s%df, s%eel)

            if(s%df%VXI_m_done .and. allocated(s%df%VXI_m)) then
                C_ommp_get_df_VXI_m = c_loc(s%df%VXI_m)
            else
                C_ommp_get_df_VXI_m = c_null_ptr
            end if
        end function C_ommp_get_df_VXI_m

        function C_ommp_get_df_VXI_p(s_prt) bind(c, name='ommp_get_df_VXI_p')
            !! Return the c-pointer to the projected dipole quantity
            !! VXI_p = V_p2q @ Xinv. Null if not available.
            !! Triggers computation via the Fortran interface.
            use mod_density_fit, only: df_project_dipoles

            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            type(c_ptr) :: C_ommp_get_df_VXI_p

            call c_f_pointer(s_prt, s)
            call df_project_dipoles(s%df, s%eel)

            if(s%df%VXI_p_done .and. allocated(s%df%VXI_p)) then
                C_ommp_get_df_VXI_p = c_loc(s%df%VXI_p)
            else
                C_ommp_get_df_VXI_p = c_null_ptr
            end if
        end function C_ommp_get_df_VXI_p

        subroutine C_ommp_df_compute_induced_dipoles(s_prt, solver, matv, add_mm_field, &
                                                      add_nuclei_field, qm_helper_prt, &
                                                      exclude_df_field) &
                bind(c, name='ommp_df_compute_induced_dipoles')
            !! Compute induced dipoles from fitted charges electric field.
            !! Wrapper for ommp_df_compute_induced_dipoles.
            type(c_ptr), value :: s_prt
            integer(ommp_integer), value :: solver
            integer(ommp_integer), value :: matv
            integer(ommp_integer), value :: add_mm_field
            integer(ommp_integer), value :: add_nuclei_field
            type(c_ptr), value :: qm_helper_prt
            integer(ommp_integer), value :: exclude_df_field
            type(ommp_system), pointer :: s
            type(ommp_qm_helper), pointer :: qm_help
            logical :: do_mm_f
            logical :: do_nuc_f
            logical :: do_exc_df_f
            logical :: has_qm

            call c_f_pointer(s_prt, s)
            do_mm_f = (add_mm_field /= 0)
            do_nuc_f = (add_nuclei_field /= 0)
            do_exc_df_f = (exclude_df_field /= 0)

            has_qm = (c_associated(qm_helper_prt))

            if(has_qm) then
                call c_f_pointer(qm_helper_prt, qm_help)
            end if

            if(has_qm) then
                call ommp_df_compute_induced_dipoles(s, solver, matv, do_mm_f, &
                                                    do_nuc_f, qm_help, &
                                                    do_exc_df_f)
            else
                call ommp_df_compute_induced_dipoles(s, solver, matv, do_mm_f, &
                                                    do_nuc_f, exclude_df_field=do_exc_df_f)
            end if
        end subroutine C_ommp_df_compute_induced_dipoles

        function C_ommp_get_df_e_field_pol_ene(s_prt) bind(c, name='ommp_get_df_e_field_pol_ene')
            !! Return the polarization energy from fitted-charge electric field.
            !! Wrapper for ommp_df_get_e_field_pol_ene.
            type(c_ptr), value :: s_prt
            real(c_double) :: C_ommp_get_df_e_field_pol_ene

            type(ommp_system), pointer :: s
            real(ommp_real) :: ene

            call c_f_pointer(s_prt, s)
            call ommp_df_get_e_field_pol_ene(s, ene)
            C_ommp_get_df_e_field_pol_ene = real(ene, c_double)
        end function C_ommp_get_df_e_field_pol_ene

        function C_ommp_get_df_E_q2p(s_prt) bind(c, name='ommp_get_df_E_q2p')
            !! Return the c-pointer to the electric field from fitted charges
            !! at polarizable sites array (3 x n_polarizable_atoms).
            !! Null if not available or not computed.
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            type(c_ptr) :: C_ommp_get_df_E_q2p

            call c_f_pointer(s_prt, s)
            if(allocated(s%df) .and. s%df%E_q2p_done) then
                C_ommp_get_df_E_q2p = c_loc(s%df%E_q2p)
            else
                C_ommp_get_df_E_q2p = c_null_ptr
            end if
        end function C_ommp_get_df_E_q2p

        subroutine C_ommp_df_geomgrad(s_prt, qmg_prt, mmg_prt, ef_prt) &
                bind(c, name='ommp_df_geomgrad')
            !! Compute the gradient (force) contribution from density fitting.
            !! Wrapper for ommp_df_geomgrad.
            !!
            !! qmg (3, n_qm)   - QM atom gradient/output array
            !! mmg (3, n_mm)   - MM atom gradient/output array
            !! ef  (3, n_pts)  - electric field at fitting points (input)
            !!
            !! The electric field at grid points and Lagrange multipliers
            !! are computed internally.

            implicit none

            type(c_ptr), value :: s_prt
            type(c_ptr), value :: qmg_prt
            type(c_ptr), value :: mmg_prt
            type(c_ptr), value :: ef_prt

            type(ommp_system), pointer :: s
            real(ommp_real), pointer :: qmg(:,:)
            real(ommp_real), pointer :: mmg(:,:)
            real(ommp_real), pointer :: ef(:,:)

            call c_f_pointer(s_prt, s)
            call c_f_pointer(qmg_prt, qmg, [3_ommp_integer, s%df%qm_top%mm_atoms])
            call c_f_pointer(mmg_prt, mmg, [3_ommp_integer, s%df%mm_top%mm_atoms])
            call c_f_pointer(ef_prt, ef, [3_ommp_integer, s%df%n_pts])

            call ommp_df_geomgrad(s, qmg, mmg, .true., .true., ef)
        end subroutine C_ommp_df_geomgrad

        function C_ommp_get_df_dX_dr(s_prt) result(ptr) bind(c, name='ommp_get_df_dX_dr')
            !! Get the dX_dr matrix (n_charges x 3 x n_pts).
            !! Returns c_null_ptr if dX_dr is not yet computed.
            type(c_ptr), value :: s_prt
            type(c_ptr) :: ptr
            type(ommp_system), pointer :: s
            call c_f_pointer(s_prt, s)
            if(.not. s%df%dX_dr_done) then
                ptr = c_null_ptr
            else
                ptr = c_loc(s%df%dX_dr)
            end if
        end function C_ommp_get_df_dX_dr

        function C_ommp_get_df_nabla_g_mm(s_prt, is_null, is_identity) result(ptr) &
                bind(c, name='ommp_get_df_nabla_g_mm')

            !! Get the gradient-grid w.r.t. MM coordinates nabla matrix.
            type(c_ptr), value :: s_prt
            logical(c_bool), intent(out) :: is_null, is_identity
            type(c_ptr) :: ptr
            type(ommp_system), pointer :: s

            call c_f_pointer(s_prt, s)

            if(.not. s%df%nabla_done) call ommp_fatal("Nabla Matrices are not available call ommp_df_geomgrad first.")

            is_null = s%df%nabla_g_mm_is_null
            is_identity = s%df%nabla_g_mm_is_identity

            if(is_null .or. is_identity .or. &
               .not. allocated(s%df%nabla_g_mm)) then
                ptr = c_null_ptr
            else
                ptr = c_loc(s%df%nabla_g_mm)
            end if
        end function C_ommp_get_df_nabla_g_mm

        function C_ommp_get_df_nabla_g_qm(s_prt, is_null, is_identity) result(ptr) &
                bind(c, name='ommp_get_df_nabla_g_qm')

            !! Get the gradient-grid w.r.t. QM coordinates nabla matrix.
            type(c_ptr), value :: s_prt
            logical(c_bool), intent(out) :: is_null, is_identity
            type(c_ptr) :: ptr
            type(ommp_system), pointer :: s

            call c_f_pointer(s_prt, s)

            if(.not. s%df%nabla_done) call ommp_fatal("Nabla Matrices are not available call ommp_df_geomgrad first.")

            is_null = s%df%nabla_g_qm_is_null
            is_identity = s%df%nabla_g_qm_is_identity

            if(is_null .or. is_identity .or. &
               .not. allocated(s%df%nabla_g_qm)) then
                ptr = c_null_ptr
            else
                ptr = c_loc(s%df%nabla_g_qm)
            end if
        end function C_ommp_get_df_nabla_g_qm

        function C_ommp_get_df_nabla_q_qm(s_prt, is_null, is_identity, is_sparse) result(ptr) &
                bind(c, name='ommp_get_df_nabla_q_qm')

            !! Get the fit-charge w.r.t. QM coordinates nabla matrix.
            type(c_ptr), value :: s_prt
            logical(c_bool), intent(out) :: is_null, is_identity, is_sparse
            type(c_ptr) :: ptr
            type(ommp_system), pointer :: s

            call c_f_pointer(s_prt, s)

            if(.not. s%df%nabla_done) call ommp_fatal("Nabla Matrices are not available call ommp_df_geomgrad first.")

            is_null = s%df%nabla_q_qm_is_null
            is_identity = s%df%nabla_q_qm_is_identity
            is_sparse = s%df%nabla_q_qm_is_sparse

            if(is_sparse) then
                ptr = c_loc(s%df%nabla_q_qm_sparse)
            else if(allocated(s%df%nabla_q_qm)) then
                ptr = c_loc(s%df%nabla_q_qm)
            else
                ptr = c_null_ptr
            end if
        end function C_ommp_get_df_nabla_q_qm

        function C_ommp_get_df_nabla_q_mm(s_prt, is_null, is_identity) result(ptr) &
                bind(c, name='ommp_get_df_nabla_q_mm')

            !! Get the fit-charge w.r.t. MM coordinates nabla matrix.
            type(c_ptr), value :: s_prt
            logical(c_bool), intent(out) :: is_null, is_identity
            type(c_ptr) :: ptr
            type(ommp_system), pointer :: s

            call c_f_pointer(s_prt, s)

            if(.not. s%df%nabla_done) call ommp_fatal("Nabla Matrices are not available call ommp_df_geomgrad first.")

            is_null = s%df%nabla_q_mm_is_null
            is_identity = s%df%nabla_q_mm_is_identity

            if(is_null .or. is_identity .or. &
               .not. allocated(s%df%nabla_q_mm)) then
                ptr = c_null_ptr
            else
                ptr = c_loc(s%df%nabla_q_mm)
            end if
        end function C_ommp_get_df_nabla_q_mm

end module mod_ommp_C_interface