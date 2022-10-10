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

    use mod_mmpol, only: ommp_cmm => cmm, &
                         ommp_cpol => cpol, &
                         ommp_q => q, &
                         ommp_ipd => ipd, &
                         ommp_polar_mm => polar_mm, &
                         ommp_mm_atoms => mm_atoms, &
                         ommp_n_ipd => n_ipd, &
                         ommp_pol_atoms => pol_atoms, &
                         ommp_ld_cart => ld_cart, &
                         ommp_is_amoeba => amoeba, &
                         ommp_is_initialized => mmpol_is_init

    use mod_constants, only: OMMP_FF_AMOEBA, OMMP_FF_WANG_AL, OMMP_FF_WANG_DL, &
                             OMMP_SOLVER_CG, OMMP_SOLVER_DIIS, &
                             OMMP_SOLVER_INVERSION, OMMP_SOLVER_DEFAULT, &
                             OMMP_MATV_INCORE, OMMP_MATV_DIRECT, &
                             OMMP_MATV_DEFAULT, &
                             OMMP_AMOEBA_D, OMMP_AMOEBA_P, &
                             OMMP_VERBOSE_DEBUG, OMMP_VERBOSE_HIGH, &
                             OMMP_VERBOSE_LOW, OMMP_VERBOSE_NONE
    use mod_constants, only: OMMP_STR_CHAR_MAX

    implicit none
    
    private :: c2f_string, OMMP_STR_CHAR_MAX

    contains

        function c_ommp_is_initialized() bind(c, name='ommp_is_initialized')
            implicit none 
            logical(c_bool) :: c_ommp_is_initialized

            c_ommp_is_initialized = ommp_is_initialized
            
            return
        end function c_ommp_is_initialized
        
        function ommp_get_cmm() bind(c, name='ommp_get_cmm')
            !! Return the c-pointer to the array containing the coordinates of
            !! MM atoms.
            use mod_mmpol, only: cmm
            type(c_ptr) :: ommp_get_cmm

            ommp_get_cmm = c_loc(cmm)
        end function ommp_get_cmm

        function ommp_get_cpol() bind(c, name='ommp_get_cpol')
            !! Return the c-pointer to the array containing the coordinates of
            !! polarizable atoms.
            use mod_mmpol, only: cpol
            type(c_ptr) :: ommp_get_cpol

            ommp_get_cpol = c_loc(cpol)
        end function ommp_get_cpol

        function ommp_get_q() bind(c, name='ommp_get_q')
            !! Return the c-pointer to the array containing the static source of 
            !! the electrostatic field.
            use mod_mmpol, only: q
            type(c_ptr) :: ommp_get_q

            ommp_get_q = c_loc(q)
        end function ommp_get_q

        function ommp_get_ipd() bind(c, name='ommp_get_ipd')
            !! Return the c-pointer to the array containing the induced dipoles 
            !! on polarizable sites.
            use mod_mmpol, only: ipd
            type(c_ptr) :: ommp_get_ipd

            ommp_get_ipd = c_loc(ipd)
        end function ommp_get_ipd
        
        function ommp_get_polar_mm() bind(c, name='ommp_get_polar_mm')
            !! Return the c-pointer to the array containing the map from 
            !! polarizable to MM atoms.
            use mod_mmpol, only: polar_mm
            type(c_ptr) :: ommp_get_polar_mm

            ommp_get_polar_mm = c_loc(polar_mm)
        end function ommp_get_polar_mm

        function ommp_get_mm_atoms() bind(c, name='ommp_get_mm_atoms')
            !! Return the number of MM atoms in the system.
            use mod_mmpol, only: mm_atoms
            implicit none

            integer(ommp_integer) :: ommp_get_mm_atoms

            ommp_get_mm_atoms = mm_atoms
        end function ommp_get_mm_atoms
        
        function ommp_get_pol_atoms() bind(c, name='ommp_get_pol_atoms')
            !! Return the number of polarizable atoms in the system.
            use mod_mmpol, only: pol_atoms
            implicit none

            integer(ommp_integer) :: ommp_get_pol_atoms

            ommp_get_pol_atoms = pol_atoms
        end function ommp_get_pol_atoms

        function get_n_ipd() bind(c, name='get_n_ipd')
            !! Return the number of dipole's set for the current Force-Field.
            use mod_mmpol, only: n_ipd
            implicit none

            integer(ommp_integer) :: get_n_ipd

            get_n_ipd = n_ipd
        end function get_n_ipd

        function ommp_get_ld_cart() bind(c, name='ommp_get_ld_cart')
            use mod_mmpol, only: ld_cart
            implicit none

            integer(ommp_integer) :: ommp_get_ld_cart

            ommp_get_ld_cart = ld_cart
        end function ommp_get_ld_cart

        function ommp_ff_is_amoeba() bind(c, name='ommp_ff_is_amoeba')
            !! Return true if the current forcefield is AMOEBA, and false in
            !! all other cases.
            use mod_mmpol, only: amoeba
            implicit none

            logical(c_bool) :: ommp_ff_is_amoeba

            ommp_ff_is_amoeba = amoeba
        end function ommp_ff_is_amoeba
        
        subroutine ommp_set_verbose(verb) bind(c, name='ommp_set_verbose')
            !! Set the verbosity level of the library to verb
            use mod_io, only: set_verbosity
            implicit none 

            integer(ommp_integer), intent(in), value :: verb
            !! Requested verbosity library

            call set_verbosity(verb)
        end subroutine ommp_set_verbose

        subroutine ommp_print_summary() bind(c, name='ommp_print_summary')
            !! Print a summary of the system input on standard output.
            use mod_mmpol, only: mmpol_ommp_print_summary

            implicit none
            
            call mmpol_ommp_print_summary()

        end subroutine ommp_print_summary
        
        subroutine ommp_print_summary_to_file(filename) &
                bind(c, name='ommp_print_summary_to_file')
            !! Print a summary of the system input on file.
            use mod_mmpol, only: mmpol_ommp_print_summary

            implicit none
            
            character(kind=c_char), intent(in) :: filename(OMMP_STR_CHAR_MAX)
            !! File where the summary will be printed
            character(len=OMMP_STR_CHAR_MAX) :: output_file
            
            call c2f_string(filename, output_file)
            call mmpol_ommp_print_summary(output_file)

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
        
        subroutine ommp_init_xyz(xyzfile, prmfile) &
                bind(c, name='ommp_init_xyz')
            !! Initialize the library using a Tinker xyz and a Tinker prm
            use mod_inputloader, only : mmpol_init_from_xyz
            
            implicit none
            
            character(kind=c_char), intent(in) :: xyzfile(OMMP_STR_CHAR_MAX), & 
                                                  prmfile(OMMP_STR_CHAR_MAX)
            character(len=OMMP_STR_CHAR_MAX) :: xyz_file, prm_file

            call c2f_string(prmfile, prm_file)
            call c2f_string(xyzfile, xyz_file)
            call mmpol_init_from_xyz(xyz_file, prm_file)
        end subroutine 
        
        subroutine C_ommp_init_mmp(filename) bind(c, name='ommp_init_mmp')
            use mod_inputloader, only : mmpol_init_from_mmp
            
            implicit none
            
            character(kind=c_char), intent(in) :: filename(OMMP_STR_CHAR_MAX)
            character(len=OMMP_STR_CHAR_MAX) :: input_file

            call c2f_string(filename, input_file)
            call mmpol_init_from_mmp(input_file)
        end subroutine 
        
        subroutine ommp_init_mmp(filename)
            use mod_inputloader, only : mmpol_init_from_mmp
            
            implicit none
            
            character(len=*) :: filename

            call mmpol_init_from_mmp(trim(filename))
        end subroutine 
        
        subroutine C_ommp_save_mmp(filename) bind(c, name='ommp_save_mmp')
            use mod_mmpol, only : mmpol_save_as_mmp
            
            implicit none
            
            character(kind=c_char), intent(in) :: filename(OMMP_STR_CHAR_MAX)
            character(len=OMMP_STR_CHAR_MAX) :: output_file

            call c2f_string(filename, output_file)
            call mmpol_save_as_mmp(output_file)
        end subroutine 

        subroutine ommp_save_mmp(filename, version)
            use mod_mmpol, only : mmpol_save_as_mmp

            implicit none

            character(len=*) :: filename
            integer :: version

            call mmpol_save_as_mmp(trim(filename), version)
        end subroutine

        subroutine C_ommp_set_external_field(ext_field, solver) &
                bind(c, name='ommp_set_external_field')
            use mod_mmpol, only: pol_atoms
            
            implicit none
            
            real(ommp_real), intent(in) :: ext_field(3, pol_atoms)
            integer(ommp_integer), intent(in), value :: solver

            call ommp_set_external_field(ext_field, solver, .true.)
        end subroutine C_ommp_set_external_field
        
        subroutine C_ommp_set_external_field_nomm(ext_field, solver) &
                bind(c, name='ommp_set_external_field_nomm')
            use mod_mmpol, only: pol_atoms
            
            implicit none
            
            real(ommp_real), intent(in) :: ext_field(3, pol_atoms)
            integer(ommp_integer), intent(in), value :: solver

            call ommp_set_external_field(ext_field, solver, .false.)
        end subroutine C_ommp_set_external_field_nomm

        subroutine ommp_set_external_field(ext_field, solver, add_mm_field)
            !! This function get an external field and solve the polarization
            !! system in the presence of the provided external field.
            use mod_polarization, only: polarization, ipd_done
            use mod_mmpol, only: pol_atoms, n_ipd, ipd
            use mod_electrostatics, only: e_m2d, prepare_M2D
            use mod_memory, only: mallocate, mfree

            implicit none
            
            real(ommp_real), intent(in) :: ext_field(3, pol_atoms)
            integer(ommp_integer), intent(in), value :: solver
            logical, intent(in), value, optional :: add_mm_field

            real(ommp_real), allocatable :: ef(:,:,:)
            integer :: i
            logical :: do_mm_f

            if(present(add_mm_field)) then
                do_mm_f = add_mm_field
            else
                do_mm_f = .true.
            end if

            ipd_done = .false.

            if(do_mm_f) then
                call mallocate('ommp_get_polelec_energy [ef]', &
                               3, pol_atoms, n_ipd, ef)
                call prepare_M2D()
                do i=1, n_ipd
                    ef(:,:,i) = e_m2d(:,:,i) + ext_field
                end do
                call polarization(ef, ipd, solver)
                call mfree('ommp_get_polelec_energy [ef]', ef)
            else
                call mallocate('ommp_get_polelec_energy [ef]', &
                               3, pol_atoms, n_ipd, ef)
                
                ef(:,:,1) = ext_field
                call polarization(ef, ipd, solver, &
                                  OMMP_MATV_DEFAULT, [.true., .false.] )
                
                call mfree('ommp_get_polelec_energy [ef]', ef)
            end if
        end subroutine ommp_set_external_field

        subroutine ommp_get_polelec_energy(epol) &
                bind(c, name='ommp_get_polelec_energy')
            !! Solve the polarization equation for a certain external field
            !! and compute the interaction energy of the induced dipoles with
            !! themselves and fixed multipoles.

            use mod_polarization, only: polarization, ipd_done
            use mod_mmpol, only: ipd
            use mod_electrostatics, only: e_m2d, prepare_M2D, energy_MM_pol
            use mod_constants, only: OMMP_SOLVER_DEFAULT

            implicit none
            
            real(ommp_real), intent(out) :: epol

            if(.not. ipd_done) then
                !! Solve the polarization system without external field
                call prepare_M2D()
                call polarization(e_m2d, ipd, OMMP_SOLVER_DEFAULT)
            end if
            call energy_MM_pol(epol) 
        end subroutine

        subroutine ommp_get_fixedelec_energy(emm) &
                bind(c, name='ommp_get_fixedelec_energy')
            ! Get the interaction energy of fixed multipoles
            use mod_electrostatics, only: energy_MM_MM, energy_MM_pol

            implicit none
            real(ommp_real), intent(out) :: emm

            emm = 0.0
            
            call energy_MM_MM(emm)
        end subroutine

        subroutine ommp_potential_mm2ext(n, cext, v) &
                bind(c, name='ommp_potential_mm2ext')
            ! Compute the electric potential of static sites at
            ! arbitrary coordinates
            use mod_electrostatics, only: potential_M2E

            implicit none
            
            integer(ommp_integer), intent(in), value :: n
            real(ommp_real), intent(in) :: cext(3,n)
            real(ommp_real), intent(inout) :: v(n)
            
            call potential_M2E(cext, v)
        end subroutine

        subroutine ommp_get_urey_energy(eub) bind(c, name='ommp_get_urey_energy')
            
            use mod_bonded, only: urey_potential

            implicit none
            real(ommp_real), intent(inout) :: eub

            eub = 0.0
            call urey_potential(eub)
        end subroutine
        
        subroutine ommp_get_strbnd_energy(eba) bind(c, name='ommp_get_strbnd_energy')
            
            use mod_bonded, only: strbnd_potential

            implicit none
            real(ommp_real), intent(inout) :: eba

            eba = 0.0
            call strbnd_potential(eba)
        end subroutine
        
        subroutine ommp_get_angle_energy(eang) bind(c, name='ommp_get_angle_energy')
            
            use mod_bonded, only: angle_potential

            implicit none
            real(ommp_real), intent(inout) :: eang

            eang = 0.0
            call angle_potential(eang)
        end subroutine
        
        subroutine ommp_get_angtor_energy(eat) bind(c, name='ommp_get_angtor_energy')
            
            use mod_bonded, only: angtor_potential

            implicit none
            real(ommp_real), intent(inout) :: eat

            eat = 0.0
            call angtor_potential(eat)
        end subroutine
        
        subroutine ommp_get_strtor_energy(ebt) bind(c, name='ommp_get_strtor_energy')
            
            use mod_bonded, only: strtor_potential

            implicit none
            real(ommp_real), intent(inout) :: ebt

            ebt = 0.0
            call strtor_potential(ebt)
        end subroutine
        
        subroutine ommp_get_bond_energy(ebnd) bind(c, name='ommp_get_bond_energy')
            
            use mod_bonded, only: bond_potential

            implicit none
            real(ommp_real), intent(inout) :: ebnd

            ebnd = 0.0
            call bond_potential(ebnd)
        end subroutine
        
        subroutine ommp_get_opb_energy(eopb) bind(c, name='ommp_get_opb_energy')
            
            use mod_bonded, only: opb_potential

            implicit none
            real(ommp_real), intent(inout) :: eopb

            eopb = 0.0
            call opb_potential(eopb)
        end subroutine
        
        subroutine ommp_get_pitors_energy(epitors) bind(c, name='ommp_get_pitors_energy')
            
            use mod_bonded, only: pitors_potential

            implicit none
            real(ommp_real), intent(inout) :: epitors

            epitors = 0.0
            call pitors_potential(epitors)
        end subroutine
        
        subroutine ommp_get_torsion_energy(et) bind(c, name='ommp_get_torsion_energy')
            
            use mod_bonded, only: torsion_potential

            implicit none
            real(ommp_real), intent(inout) :: et

            et = 0.0
            call torsion_potential(et)
        end subroutine
        
        subroutine ommp_get_tortor_energy(ett) bind(c, name='ommp_get_tortor_energy')
            
            use mod_bonded, only: tortor_potential

            implicit none
            real(ommp_real), intent(inout) :: ett

            ett = 0.0
            call tortor_potential(ett)
        end subroutine
        
        subroutine ommp_get_vdw_energy(evdw) bind(c, name='ommp_get_vdw_energy')
            
            use mod_nonbonded, only: vdw_potential

            implicit none
            real(ommp_real), intent(inout) :: evdw

            evdw = 0.0
            call vdw_potential(evdw)
        end subroutine

        subroutine ommp_terminate() bind(c, name='ommp_terminate')
            use mod_mmpol, only: mmpol_terminate
            use mod_bonded, only: terminate_bonded
            use mod_nonbonded, only: vdw_terminate
            use mod_electrostatics, only: electrostatics_terminate
            use mod_polarization, only: polarization_terminate

            implicit none

            call mmpol_terminate()
            call terminate_bonded()
            call vdw_terminate()
            call electrostatics_terminate()
            call polarization_terminate()

        end subroutine

#ifdef USE_HDF5
        subroutine ommp_init_hdf5(filename) bind(c, name='ommp_init_hdf5')
            !! This function is an interface for saving an HDF5 file 
            !! with all the data contained in mmpol module using
            !! [[mod_io:mmpol_save_as_hdf5]]
            use mod_iohdf5, only: mmpol_init_from_hdf5

            implicit none
            
            character(kind=c_char), intent(in) :: filename(OMMP_STR_CHAR_MAX)
            character(len=OMMP_STR_CHAR_MAX) :: hdf5in
            integer(ommp_integer) :: ok

            call c2f_string(filename, hdf5in)
            call mmpol_init_from_hdf5(hdf5in, ok)
            
        end subroutine ommp_init_hdf5
        
        function ommp_write_hdf5(filename) bind(c, name='ommp_write_hdf5')
            !! This function is an interface for saving an HDF5 file 
            !! with all the data contained in mmpol module using
            !! [[mod_io:mmpol_save_as_hdf5]]
            use mod_iohdf5, only: mmpol_save_as_hdf5

            implicit none
            
            character(kind=c_char), intent(in) :: filename(OMMP_STR_CHAR_MAX)
            character(len=OMMP_STR_CHAR_MAX) :: hdf5out
            integer(ommp_integer) :: ommp_write_hdf5

            call c2f_string(filename, hdf5out)
            call mmpol_save_as_hdf5(hdf5out, ommp_write_hdf5)
            
        end function ommp_write_hdf5
#endif

end module mod_ommp_interface

