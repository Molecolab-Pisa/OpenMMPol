! Wrapper function for open-mmpol library
module mod_ommp_interface
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
    
    ! Renamed import of several global variable that should be available
    ! in the interface
    use mod_memory, only: ommp_integer => ip, &
                          ommp_real => rp

    use mod_constants, only: OMMP_FF_AMOEBA, OMMP_FF_WANG_AL, OMMP_FF_WANG_DL, &
                             OMMP_SOLVER_CG, OMMP_SOLVER_DIIS, &
                             OMMP_SOLVER_INVERSION, OMMP_SOLVER_DEFAULT, &
                             OMMP_MATV_INCORE, OMMP_MATV_DIRECT, &
                             OMMP_MATV_DEFAULT, &
                             OMMP_AMOEBA_D, OMMP_AMOEBA_P, &
                             OMMP_VERBOSE_DEBUG, OMMP_VERBOSE_HIGH, &
                             OMMP_VERBOSE_LOW, OMMP_VERBOSE_NONE
    use mod_constants, only: OMMP_STR_CHAR_MAX
    use mod_mmpol, only: ommp_system

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

        !!TODO function ommp_get_ipd() bind(c, name='ommp_get_ipd')
        !!TODO     !! Return the c-pointer to the array containing the induced dipoles 
        !!TODO     !! on polarizable sites.
        !!TODO     use mod_mmpol, only: ipd
        !!TODO     type(c_ptr) :: ommp_get_ipd

        !!TODO     ommp_get_ipd = c_loc(ipd)
        !!TODO end function ommp_get_ipd
        !!TODO 
        !!TODO function ommp_get_polar_mm() bind(c, name='ommp_get_polar_mm')
        !!TODO     !! Return the c-pointer to the array containing the map from 
        !!TODO     !! polarizable to MM atoms.
        !!TODO     use mod_mmpol, only: polar_mm
        !!TODO     type(c_ptr) :: ommp_get_polar_mm

        !!TODO     ommp_get_polar_mm = c_loc(polar_mm)
        !!TODO end function ommp_get_polar_mm

        !!TODO function ommp_get_mm_atoms() bind(c, name='ommp_get_mm_atoms')
        !!TODO     !! Return the number of MM atoms in the system.
        !!TODO     use mod_mmpol, only: mm_atoms
        !!TODO     implicit none

        !!TODO     integer(ommp_integer) :: ommp_get_mm_atoms

        !!TODO     ommp_get_mm_atoms = mm_atoms
        !!TODO end function ommp_get_mm_atoms
        !!TODO 
        function ommp_get_pol_atoms(s_prt) bind(c, name='ommp_get_pol_atoms')
            !! Return the number of polarizable atoms in the system.
            implicit none

            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            integer(ommp_integer) :: ommp_get_pol_atoms
            
            call c_f_pointer(s_prt, s)

            ommp_get_pol_atoms = s%eel%pol_atoms
        end function ommp_get_pol_atoms

        !!TODO function get_n_ipd() bind(c, name='get_n_ipd')
        !!TODO     !! Return the number of dipole's set for the current Force-Field.
        !!TODO     use mod_mmpol, only: n_ipd
        !!TODO     implicit none

        !!TODO     integer(ommp_integer) :: get_n_ipd

        !!TODO     get_n_ipd = n_ipd
        !!TODO end function get_n_ipd

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
        
        subroutine ommp_set_verbose(s_prt, verb) bind(c, name='ommp_set_verbose')
            !! Set the verbosity level of the library to verb
            use mod_io, only: set_verbosity
            implicit none 

            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
            integer(ommp_integer), intent(in), value :: verb
            !! Requested verbosity library
            call c_f_pointer(s_prt, s)
            call set_verbosity(verb)
        end subroutine ommp_set_verbose

        subroutine ommp_print_summary(s_prt) bind(c, name='ommp_print_summary')
            !! Print a summary of the system input on standard output.
            use mod_mmpol, only: mmpol_ommp_print_summary

            implicit none
            
            type(c_ptr), value :: s_prt
            type(ommp_system), pointer :: s
           
            call c_f_pointer(s_prt, s)
            call mmpol_ommp_print_summary(s)

        end subroutine ommp_print_summary
        
        subroutine ommp_print_summary_to_file(s_prt, filename) &
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

        end subroutine ommp_print_summary_to_file

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
            
            type(ommp_system), save, allocatable, target :: s
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
            
            type(ommp_system), save, allocatable, target :: s
            character(kind=c_char), intent(in) :: filename(OMMP_STR_CHAR_MAX)
            character(len=OMMP_STR_CHAR_MAX) :: input_file
            type(c_ptr) :: c_prt
            
            allocate(s)
            
            call c2f_string(filename, input_file)
            call mmpol_init_from_mmp(input_file, s)
            c_prt = c_loc(s)
        end function
        
!TODO        subroutine ommp_init_mmp(filename)
!TODO            use mod_inputloader, only : mmpol_init_from_mmp
!TODO            
!TODO            implicit none
!TODO            
!TODO            character(len=*) :: filename
!TODO
!TODO            call mmpol_init_from_mmp(trim(filename))
!TODO        end subroutine 
!TODO        
!TODO        subroutine C_ommp_save_mmp(filename) bind(c, name='ommp_save_mmp')
!TODO            use mod_mmpol, only : mmpol_save_as_mmp
!TODO            
!TODO            implicit none
!TODO            
!TODO            character(kind=c_char), intent(in) :: filename(OMMP_STR_CHAR_MAX)
!TODO            character(len=OMMP_STR_CHAR_MAX) :: output_file
!TODO
!TODO            call c2f_string(filename, output_file)
!TODO            call mmpol_save_as_mmp(output_file)
!TODO        end subroutine 
!TODO
!TODO        subroutine ommp_save_mmp(filename, version)
!TODO            use mod_mmpol, only : mmpol_save_as_mmp
!TODO
!TODO            implicit none
!TODO
!TODO            character(len=*) :: filename
!TODO            integer :: version
!TODO
!TODO            call mmpol_save_as_mmp(trim(filename), version)
!TODO        end subroutine
!TODO
!TODO        subroutine C_ommp_set_external_field(ext_field, solver) &
!TODO                bind(c, name='ommp_set_external_field')
!TODO            use mod_mmpol, only: pol_atoms
!TODO            
!TODO            implicit none
!TODO            
!TODO            real(ommp_real), intent(in) :: ext_field(3, pol_atoms)
!TODO            integer(ommp_integer), intent(in), value :: solver
!TODO
!TODO            call ommp_set_external_field(ext_field, solver, .true.)
!TODO        end subroutine C_ommp_set_external_field
!TODO        
!TODO        subroutine C_ommp_set_external_field_nomm(ext_field, solver) &
!TODO                bind(c, name='ommp_set_external_field_nomm')
!TODO            use mod_mmpol, only: pol_atoms
!TODO            
!TODO            implicit none
!TODO            
!TODO            real(ommp_real), intent(in) :: ext_field(3, pol_atoms)
!TODO            integer(ommp_integer), intent(in), value :: solver
!TODO
!TODO            call ommp_set_external_field(ext_field, solver, .false.)
!TODO        end subroutine C_ommp_set_external_field_nomm
!TODO
!TODO        subroutine ommp_set_external_field(ext_field, solver, add_mm_field)
!TODO            !! This function get an external field and solve the polarization
!TODO            !! system in the presence of the provided external field.
!TODO            use mod_polarization, only: polarization, ipd_done
!TODO            use mod_mmpol, only: pol_atoms, n_ipd, ipd
!TODO            use mod_electrostatics, only: e_m2d, prepare_M2D
!TODO            use mod_memory, only: mallocate, mfree
!TODO
!TODO            implicit none
!TODO            
!TODO            real(ommp_real), intent(in) :: ext_field(3, pol_atoms)
!TODO            integer(ommp_integer), intent(in), value :: solver
!TODO            logical, intent(in), value, optional :: add_mm_field
!TODO
!TODO            real(ommp_real), allocatable :: ef(:,:,:)
!TODO            integer :: i
!TODO            logical :: do_mm_f
!TODO
!TODO            if(present(add_mm_field)) then
!TODO                do_mm_f = add_mm_field
!TODO            else
!TODO                do_mm_f = .true.
!TODO            end if
!TODO
!TODO            ipd_done = .false.
!TODO
!TODO            if(do_mm_f) then
!TODO                call mallocate('ommp_get_polelec_energy [ef]', &
!TODO                               3, pol_atoms, n_ipd, ef)
!TODO                call prepare_M2D()
!TODO                do i=1, n_ipd
!TODO                    ef(:,:,i) = e_m2d(:,:,i) + ext_field
!TODO                end do
!TODO                call polarization(ef, ipd, solver)
!TODO                call mfree('ommp_get_polelec_energy [ef]', ef)
!TODO            else
!TODO                call mallocate('ommp_get_polelec_energy [ef]', &
!TODO                               3, pol_atoms, n_ipd, ef)
!TODO                
!TODO                ef(:,:,1) = ext_field
!TODO                call polarization(ef, ipd, solver, &
!TODO                                  OMMP_MATV_DEFAULT, [.true., .false.] )
!TODO                
!TODO                call mfree('ommp_get_polelec_energy [ef]', ef)
!TODO            end if
!TODO        end subroutine ommp_set_external_field
!TODO
!TODO        subroutine ommp_get_polelec_energy(epol) &
!TODO                bind(c, name='ommp_get_polelec_energy')
!TODO            !! Solve the polarization equation for a certain external field
!TODO            !! and compute the interaction energy of the induced dipoles with
!TODO            !! themselves and fixed multipoles.
!TODO
!TODO            use mod_polarization, only: polarization, ipd_done
!TODO            use mod_mmpol, only: ipd
!TODO            use mod_electrostatics, only: e_m2d, prepare_M2D, energy_MM_pol
!TODO            use mod_constants, only: OMMP_SOLVER_DEFAULT
!TODO
!TODO            implicit none
!TODO            
!TODO            real(ommp_real), intent(out) :: epol
!TODO
!TODO            if(.not. ipd_done) then
!TODO                !! Solve the polarization system without external field
!TODO                call prepare_M2D()
!TODO                call polarization(e_m2d, ipd, OMMP_SOLVER_DEFAULT)
!TODO            end if
!TODO            call energy_MM_pol(epol) 
!TODO        end subroutine
!TODO
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
!TODO        subroutine ommp_get_urey_energy(eub) bind(c, name='ommp_get_urey_energy')
!TODO            
!TODO            use mod_bonded, only: urey_potential
!TODO
!TODO            implicit none
!TODO            real(ommp_real), intent(inout) :: eub
!TODO
!TODO            eub = 0.0
!TODO            call urey_potential(eub)
!TODO        end subroutine
!TODO        
!TODO        subroutine ommp_get_strbnd_energy(eba) bind(c, name='ommp_get_strbnd_energy')
!TODO            
!TODO            use mod_bonded, only: strbnd_potential
!TODO
!TODO            implicit none
!TODO            real(ommp_real), intent(inout) :: eba
!TODO
!TODO            eba = 0.0
!TODO            call strbnd_potential(eba)
!TODO        end subroutine
!TODO        
!TODO        subroutine ommp_get_angle_energy(eang) bind(c, name='ommp_get_angle_energy')
!TODO            
!TODO            use mod_bonded, only: angle_potential
!TODO
!TODO            implicit none
!TODO            real(ommp_real), intent(inout) :: eang
!TODO
!TODO            eang = 0.0
!TODO            call angle_potential(eang)
!TODO        end subroutine
!TODO        
!TODO        subroutine ommp_get_angtor_energy(eat) bind(c, name='ommp_get_angtor_energy')
!TODO            
!TODO            use mod_bonded, only: angtor_potential
!TODO
!TODO            implicit none
!TODO            real(ommp_real), intent(inout) :: eat
!TODO
!TODO            eat = 0.0
!TODO            call angtor_potential(eat)
!TODO        end subroutine
!TODO        
!TODO        subroutine ommp_get_strtor_energy(ebt) bind(c, name='ommp_get_strtor_energy')
!TODO            
!TODO            use mod_bonded, only: strtor_potential
!TODO
!TODO            implicit none
!TODO            real(ommp_real), intent(inout) :: ebt
!TODO
!TODO            ebt = 0.0
!TODO            call strtor_potential(ebt)
!TODO        end subroutine
!TODO        
!TODO        subroutine ommp_get_bond_energy(ebnd) bind(c, name='ommp_get_bond_energy')
!TODO            
!TODO            use mod_bonded, only: bond_potential
!TODO
!TODO            implicit none
!TODO            real(ommp_real), intent(inout) :: ebnd
!TODO
!TODO            ebnd = 0.0
!TODO            call bond_potential(ebnd)
!TODO        end subroutine
!TODO        
!TODO        subroutine ommp_get_opb_energy(eopb) bind(c, name='ommp_get_opb_energy')
!TODO            
!TODO            use mod_bonded, only: opb_potential
!TODO
!TODO            implicit none
!TODO            real(ommp_real), intent(inout) :: eopb
!TODO
!TODO            eopb = 0.0
!TODO            call opb_potential(eopb)
!TODO        end subroutine
!TODO        
!TODO        subroutine ommp_get_pitors_energy(epitors) bind(c, name='ommp_get_pitors_energy')
!TODO            
!TODO            use mod_bonded, only: pitors_potential
!TODO
!TODO            implicit none
!TODO            real(ommp_real), intent(inout) :: epitors
!TODO
!TODO            epitors = 0.0
!TODO            call pitors_potential(epitors)
!TODO        end subroutine
!TODO        
!TODO        subroutine ommp_get_torsion_energy(et) bind(c, name='ommp_get_torsion_energy')
!TODO            
!TODO            use mod_bonded, only: torsion_potential
!TODO
!TODO            implicit none
!TODO            real(ommp_real), intent(inout) :: et
!TODO
!TODO            et = 0.0
!TODO            call torsion_potential(et)
!TODO        end subroutine
!TODO        
!TODO        subroutine ommp_get_tortor_energy(ett) bind(c, name='ommp_get_tortor_energy')
!TODO            
!TODO            use mod_bonded, only: tortor_potential
!TODO
!TODO            implicit none
!TODO            real(ommp_real), intent(inout) :: ett
!TODO
!TODO            ett = 0.0
!TODO            call tortor_potential(ett)
!TODO        end subroutine
!TODO        
!TODO        subroutine ommp_get_vdw_energy(evdw) bind(c, name='ommp_get_vdw_energy')
!TODO            
!TODO            use mod_nonbonded, only: vdw_potential
!TODO
!TODO            implicit none
!TODO            real(ommp_real), intent(inout) :: evdw
!TODO
!TODO            evdw = 0.0
!TODO            call vdw_potential(evdw)
!TODO        end subroutine
!TODO
        subroutine ommp_terminate(s_prt) bind(c, name='ommp_terminate')
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

end module mod_ommp_interface

