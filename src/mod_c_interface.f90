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

        !!TODO function c_ommp_is_initialized() bind(c, name='ommp_is_initialized')
        !!TODO     implicit none 
        !!TODO     logical(c_bool) :: c_ommp_is_initialized

        !!TODO     c_ommp_is_initialized = ommp_is_initialized
        !!TODO     
        !!TODO     return
        !!TODO end function c_ommp_is_initialized
        !!TODO 
        !!TODO function ommp_get_cmm() bind(c, name='ommp_get_cmm')
        !!TODO     !! Return the c-pointer to the array containing the coordinates of
        !!TODO     !! MM atoms.
        !!TODO     use mod_mmpol, only: cmm
        !!TODO     type(c_ptr) :: ommp_get_cmm

        !!TODO     ommp_get_cmm = c_loc(cmm)
        !!TODO end function ommp_get_cmm

        !!TODO function ommp_get_cpol() bind(c, name='ommp_get_cpol')
        !!TODO     !! Return the c-pointer to the array containing the coordinates of
        !!TODO     !! polarizable atoms.
        !!TODO     use mod_mmpol, only: cpol
        !!TODO     type(c_ptr) :: ommp_get_cpol

        !!TODO     ommp_get_cpol = c_loc(cpol)
        !!TODO end function ommp_get_cpol

        !!TODO function ommp_get_q() bind(c, name='ommp_get_q')
        !!TODO     !! Return the c-pointer to the array containing the static source of 
        !!TODO     !! the electrostatic field.
        !!TODO     use mod_mmpol, only: q
        !!TODO     type(c_ptr) :: ommp_get_q

        !!TODO     ommp_get_q = c_loc(q)
        !!TODO end function ommp_get_q

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
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            type(c_ptr) :: C_ommp_get_polar_mm

            call c_f_pointer(s_prt, s)
            C_ommp_get_polar_mm = c_loc(s%eel%polar_mm)
        end function C_ommp_get_polar_mm

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

        function ommp_get_n_ipd(s_prt) bind(c, name='ommp_get_n_ipd')
            !! Return the number of dipole's set for the current Force-Field.
            implicit none

            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            integer(ommp_integer) :: ommp_get_n_ipd

            call c_f_pointer(s_prt, s)
            ommp_get_n_ipd = s%eel%n_ipd
        end function ommp_get_n_ipd

        !!TODO function ommp_get_ld_cart() bind(c, name='ommp_get_ld_cart')
        !!TODO     use mod_mmpol, only: ld_cart
        !!TODO     implicit none

        !!TODO     integer(ommp_integer) :: ommp_get_ld_cart

        !!TODO     ommp_get_ld_cart = ld_cart
        !!TODO end function ommp_get_ld_cart

        !!TODO function ommp_ff_is_amoeba() bind(c, name='ommp_ff_is_amoeba')
        !!TODO     !! Return true if the current forcefield is AMOEBA, and false in
        !!TODO     !! all other cases.
        !!TODO     use mod_mmpol, only: amoeba
        !!TODO     implicit none

        !!TODO     logical(c_bool) :: ommp_ff_is_amoeba

        !!TODO     ommp_ff_is_amoeba = amoeba
        !!TODO end function ommp_ff_is_amoeba
        
        subroutine C_ommp_set_verbose(verb) bind(c, name='ommp_set_verbose')
            !! Set the verbosity level of the library to verb
            use mod_io, only: set_verbosity
            implicit none 

            integer(ommp_integer), intent(in), value :: verb
            
            !! Requested verbosityi of library
            call set_verbosity(verb)
        end subroutine C_ommp_set_verbose

        subroutine C_ommp_print_summary(s_prt) bind(c, name='ommp_print_summary')
            !! Print a summary of the system input on standard output.
            use mod_mmpol, only: mmpol_ommp_print_summary

            implicit none
            
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
           
            call c_f_pointer(s_prt, s)
            call mmpol_ommp_print_summary(s)

        end subroutine C_ommp_print_summary
        
        subroutine C_ommp_print_summary_to_file(s_prt, filename) &
                bind(c, name='ommp_print_summary_to_file')
            !! Print a summary of the system input on file.
            use mod_mmpol, only: mmpol_ommp_print_summary

            implicit none
            
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            character(kind=c_char), intent(in) :: filename(OMMP_STR_CHAR_MAX)
            !! File where the summary will be printed
            character(len=OMMP_STR_CHAR_MAX) :: output_file
            
            call c_f_pointer(s_prt, s)
            call c2f_string(filename, output_file)
            call mmpol_ommp_print_summary(s, output_file)

        end subroutine C_ommp_print_summary_to_file

        subroutine c2f_string(c_str, f_str)
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
        
        function C_ommp_init_xyz(xyzfile, prmfile) &
                result(c_prt) bind(c, name='ommp_init_xyz')
            !! Initialize the library using a Tinker xyz and a Tinker prm
            use mod_inputloader, only : mmpol_init_from_xyz
            
            implicit none
            
            type(ommp_system), pointer :: s
            character(kind=c_char), intent(in) :: xyzfile(OMMP_STR_CHAR_MAX), & 
                                                  prmfile(OMMP_STR_CHAR_MAX)
            character(len=OMMP_STR_CHAR_MAX) :: xyz_file, prm_file
            type(c_ptr) :: c_prt

            allocate(s)
            
            call c2f_string(prmfile, prm_file)
            call c2f_string(xyzfile, xyz_file)
            call mmpol_init_from_xyz(s, xyz_file, prm_file)
            c_prt = c_loc(s)
        end function
        
        function C_ommp_init_mmp(filename) &
                result(c_prt) bind(c, name='ommp_init_mmp')
            
            use mod_inputloader, only : mmpol_init_from_mmp
            
            implicit none

            type(ommp_system), pointer :: s
            character(kind=c_char), intent(in) :: filename(OMMP_STR_CHAR_MAX)
            character(len=OMMP_STR_CHAR_MAX) :: input_file
            type(c_ptr) :: c_prt
            
            allocate(s)
            
            call c2f_string(filename, input_file)
            call mmpol_init_from_mmp(input_file, s)
            c_prt = c_loc(s)
        end function
        
        subroutine C_ommp_save_mmp(s_prt, filename, version) &
                   bind(c, name='ommp_save_mmp')
            use mod_mmpol, only : mmpol_save_as_mmp
            
            implicit none
            type(c_ptr), value :: s_prt
            character(kind=c_char), intent(in) :: filename(OMMP_STR_CHAR_MAX)
            integer(ommp_integer), value :: version

            character(len=OMMP_STR_CHAR_MAX) :: output_file
            type(ommp_system), pointer :: s

            call c_f_pointer(s_prt, s)

            call c2f_string(filename, output_file)
            call mmpol_save_as_mmp(s, output_file, version)
        end subroutine 
        
        subroutine C_ommp_set_external_field(s_prt, ext_field_prt, solver) &
                bind(c, name='ommp_set_external_field')
            !!use mod_mmpol, only: pol_atoms
            
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

        function C_ommp_get_polelec_energy(s_prt) &
                result(epol) bind(c, name='ommp_get_polelec_energy')
            !! Solve the polarization equation for a certain external field
            !! and compute the interaction energy of the induced dipoles with
            !! themselves and fixed multipoles.

            use mod_polarization, only: polarization
            use mod_electrostatics, only: prepare_M2D, energy_MM_pol, &
                                          ommp_electrostatics_type
            use mod_constants, only: OMMP_SOLVER_DEFAULT

            implicit none
            
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            type(ommp_electrostatics_type), pointer :: eel
            real(ommp_real) :: epol

            call c_f_pointer(s_prt, s)
            eel => s%eel
            
            if(.not. eel%ipd_done) then
                !! Solve the polarization system without external field
                call prepare_M2D(eel)
                call polarization(s, eel%e_m2d, OMMP_SOLVER_DEFAULT)
            end if
            epol = 0.0
            call energy_MM_pol(eel, epol) 
        end function

        function C_ommp_get_fixedelec_energy(s_prt) &
                result(emm) bind(c, name='ommp_get_fixedelec_energy')
            ! Get the interaction energy of fixed multipoles
            use mod_electrostatics, only: energy_MM_MM

            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            real(ommp_real) :: emm

            call c_f_pointer(s_prt, s)

            emm = 0.0
            call energy_MM_MM(s%eel, emm)
        end function
!TODO
!TODO        subroutine ommp_potential_mmpol2ext(n, cext, v) &
!TODO                bind(c, name='ommp_potential_mmpol2ext')
!TODO            ! Compute the electric potential of static sites at
!TODO            ! arbitrary coordinates
!TODO            use mod_electrostatics, only: potential_D2E, &
!TODO                                          potential_M2E
!TODO
!TODO            implicit none
!TODO            
!TODO            integer(ommp_integer), intent(in), value :: n
!TODO            real(ommp_real), intent(in) :: cext(3,n)
!TODO            real(ommp_real), intent(inout) :: v(n)
!TODO            
!TODO            call potential_M2E(cext, v)
!TODO            call potential_D2E(cext, v)
!TODO        end subroutine
!TODO        
!TODO        subroutine ommp_potential_pol2ext(n, cext, v) &
!TODO                bind(c, name='ommp_potential_pol2ext')
!TODO            ! Compute the electric potential of static sites at
!TODO            ! arbitrary coordinates
!TODO            use mod_electrostatics, only: potential_D2E
!TODO
!TODO            implicit none
!TODO            
!TODO            integer(ommp_integer), intent(in), value :: n
!TODO            real(ommp_real), intent(in) :: cext(3,n)
!TODO            real(ommp_real), intent(inout) :: v(n)
!TODO            
!TODO            call potential_D2E(cext, v)
!TODO        end subroutine
!TODO        
!TODO        subroutine ommp_potential_mm2ext(n, cext, v) &
!TODO                bind(c, name='ommp_potential_mm2ext')
!TODO            ! Compute the electric potential of static sites at
!TODO            ! arbitrary coordinates
!TODO            use mod_electrostatics, only: potential_M2E
!TODO
!TODO            implicit none
!TODO            
!TODO            integer(ommp_integer), intent(in), value :: n
!TODO            real(ommp_real), intent(in) :: cext(3,n)
!TODO            real(ommp_real), intent(inout) :: v(n)
!TODO            
!TODO            call potential_M2E(cext, v)
!TODO        end subroutine
!TODO
        function C_ommp_get_urey_energy(s_prt) &
                result(eub) bind(c, name='ommp_get_urey_energy')
            
            !! use mod_bonded, only: urey_potential

            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            real(ommp_real) :: eub

            call c_f_pointer(s_prt, s)

            eub = 0.0
            !! call urey_potential(eub)
        end function
        
        function C_ommp_get_strbnd_energy(s_prt) &
                result(eba) bind(c, name='ommp_get_strbnd_energy')
            
            !! use mod_bonded, only: strbnd_potential

            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            real(ommp_real) :: eba

            call c_f_pointer(s_prt, s)

            eba = 0.0
            !! call strbnd_potential(eba)
        end function
        
        function C_ommp_get_angle_energy(s_prt) &
                result(eang) bind(c, name='ommp_get_angle_energy')
            
            !!use mod_bonded, only: angle_potential

            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            real(ommp_real) :: eang

            call c_f_pointer(s_prt, s)

            eang = 0.0
            !! call angle_potential(eang)
        end function
        
        function C_ommp_get_angtor_energy(s_prt) &
                result(eat) bind(c, name='ommp_get_angtor_energy')
            
            !! use mod_bonded, only: angtor_potential

            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            real(ommp_real) :: eat

            call c_f_pointer(s_prt, s)

            eat = 0.0
            !! call angtor_potential(eat)
        end function
        
        function C_ommp_get_strtor_energy(s_prt) &
                result(ebt) bind(c, name='ommp_get_strtor_energy')
            
            !! use mod_bonded, only: strtor_potential

            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            real(ommp_real) :: ebt

            call c_f_pointer(s_prt, s)

            ebt = 0.0
            !! call strtor_potential(ebt)
        end function
        
        function C_ommp_get_bond_energy(s_prt) &
                result(ebnd) bind(c, name='ommp_get_bond_energy')
            
            !! use mod_bonded, only: bond_potential

            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            real(ommp_real) :: ebnd

            call c_f_pointer(s_prt, s)

            ebnd = 0.0
            !! call bond_potential(ebnd)
        end function
        
        function C_ommp_get_opb_energy(s_prt) &
            result(eopb) bind(c, name='ommp_get_opb_energy')
            
            !!use mod_bonded, only: opb_potential

            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            real(ommp_real) :: eopb

            call c_f_pointer(s_prt, s)

            eopb = 0.0
            !! call opb_potential(eopb)
        end function
        
        function C_ommp_get_pitors_energy(s_prt) &
                result(epitors) bind(c, name='ommp_get_pitors_energy')
            
            !! use mod_bonded, only: pitors_potential

            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            real(ommp_real) :: epitors

            call c_f_pointer(s_prt, s)

            epitors = 0.0
            !! call pitors_potential(epitors)
        end function
        
        function C_ommp_get_torsion_energy(s_prt) &
                result(et) bind(c, name='ommp_get_torsion_energy')
            
            !! use mod_bonded, only: torsion_potential

            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            real(ommp_real) :: et

            call c_f_pointer(s_prt, s)

            et = 0.0
            !! call torsion_potential(et)
        end function
        
        function C_ommp_get_tortor_energy(s_prt) &
                result(ett) bind(c, name='ommp_get_tortor_energy')
            
            !! use mod_bonded, only: tortor_potential

            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            real(ommp_real) :: ett

            call c_f_pointer(s_prt, s)

            ett = 0.0
            !! call tortor_potential(ett)
        end function
        
        function C_ommp_get_vdw_energy(s_prt) &
                result(evdw) bind(c, name='ommp_get_vdw_energy')
            
            use mod_nonbonded, only: vdw_potential
            
            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            real(ommp_real) :: evdw

            call c_f_pointer(s_prt, s)

            evdw = 0.0
            call vdw_potential(s%vdw, evdw)

        end function

        subroutine C_ommp_terminate(s_prt) bind(c, name='ommp_terminate')
            use mod_mmpol, only: mmpol_terminate

            implicit none
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s

            call c_f_pointer(s_prt, s)
            call mmpol_terminate(s)
            
            deallocate(s)

        end subroutine
!TODO
!TODO#ifdef USE_HDF5
!TODO        subroutine ommp_init_hdf5(filename) bind(c, name='ommp_init_hdf5')
!TODO            !! This function is an interface for saving an HDF5 file 
!TODO            !! with all the data contained in mmpol module using
!TODO            !! [[mod_io:mmpol_save_as_hdf5]]
!TODO            use mod_iohdf5, only: mmpol_init_from_hdf5
!TODO
!TODO            implicit none
!TODO            
!TODO            character(kind=c_char), intent(in) :: filename(OMMP_STR_CHAR_MAX)
!TODO            character(len=OMMP_STR_CHAR_MAX) :: hdf5in
!TODO            integer(ommp_integer) :: ok
!TODO
!TODO            call c2f_string(filename, hdf5in)
!TODO            call mmpol_init_from_hdf5(hdf5in, ok)
!TODO            
!TODO        end subroutine ommp_init_hdf5
!TODO        
!TODO        function ommp_write_hdf5(filename) bind(c, name='ommp_write_hdf5')
!TODO            !! This function is an interface for saving an HDF5 file 
!TODO            !! with all the data contained in mmpol module using
!TODO            !! [[mod_io:mmpol_save_as_hdf5]]
!TODO            use mod_iohdf5, only: mmpol_save_as_hdf5
!TODO
!TODO            implicit none
!TODO            
!TODO            character(kind=c_char), intent(in) :: filename(OMMP_STR_CHAR_MAX)
!TODO            character(len=OMMP_STR_CHAR_MAX) :: hdf5out
!TODO            integer(ommp_integer) :: ommp_write_hdf5
!TODO
!TODO            call c2f_string(filename, hdf5out)
!TODO            call mmpol_save_as_hdf5(hdf5out, ommp_write_hdf5)
!TODO            
!TODO        end function ommp_write_hdf5
!TODO#endif

end module mod_ommp_C_interface

