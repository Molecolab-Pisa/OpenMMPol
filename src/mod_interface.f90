! Wrapper function for open-mmpol library
module ommp_interface
    !! The interface of the library, basically all the operation performed
    !! by an external code should be done through the routines of this
    !! module. 
    !! The interface is conceived to work naturally with C and Fortran; the C
    !! interface is also used to build the interface for Python.
    !! In a fortran code, this module can be imported and it should expose 
    !! directly all the vector and scalar quantities needed.
    !! In a C code, routines are provided to get the pointer or the values of 
    !! vector and scalar quantites respectively.

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
    
    use mod_mmpol, only: ommp_system
    use mod_electrostatics, only: ommp_electrostatics_type
    use mod_topology, only: ommp_topology_type

    use mod_mmpol, only: ommp_save_mmp => mmpol_save_as_mmp, &
                         ommp_print_summary => mmpol_ommp_print_summary, &
                         ommp_print_summary_to_file => mmpol_ommp_print_summary

    use mod_io, only: ommp_set_verbose => set_verbosity

    implicit none
    
    contains
        
        subroutine ommp_init_mmp(s, filename)
            use mod_inputloader, only : mmpol_init_from_mmp
            
            implicit none
            
            type(ommp_system), pointer, intent(inout) :: s
            character(len=*) :: filename

            allocate(s)
            call mmpol_init_from_mmp(trim(filename), s)
        end subroutine
        
        subroutine ommp_init_xyz(s, xyzfile, prmfile)
            use mod_inputloader, only : mmpol_init_from_xyz
            
            implicit none
            
            type(ommp_system), pointer, intent(inout) :: s
            character(len=*) :: xyzfile, prmfile

            allocate(s)
            call mmpol_init_from_xyz(s, trim(xyzfile), trim(prmfile))
        end subroutine
        
        subroutine ommp_terminate(s)
            use mod_mmpol, only: mmpol_terminate

            implicit none
            
            type(ommp_system), pointer, intent(inout) :: s

            call mmpol_terminate(s)
            
            deallocate(s)

        end subroutine

        subroutine ommp_set_external_field(sys_obj, ext_field, solver, &
                                           add_mm_field)
            !! This function get an external field and solve the polarization
            !! system in the presence of the provided external field.
            use mod_polarization, only: polarization
            use mod_electrostatics, only: prepare_M2D
            use mod_memory, only: mallocate, mfree

            implicit none
            
            type(ommp_system), intent(inout), target :: sys_obj
            real(ommp_real), intent(in) :: ext_field(3,sys_obj%eel%pol_atoms)
            integer(ommp_integer), intent(in), value :: solver
            logical, intent(in), value, optional :: add_mm_field
            
            type(ommp_electrostatics_type), pointer :: eel
            real(ommp_real), allocatable :: ef(:,:,:)
            integer :: i
            logical :: do_mm_f

            eel => sys_obj%eel

            if(present(add_mm_field)) then
                do_mm_f = add_mm_field
            else
                do_mm_f = .true.
            end if

            eel%ipd_done = .false.

            if(do_mm_f) then
                call mallocate('ommp_get_polelec_energy [ef]', &
                               3, eel%pol_atoms, eel%n_ipd, ef)
                call prepare_M2D(eel)
                do i=1, eel%n_ipd
                    ef(:,:,i) = eel%e_m2d(:,:,i) + ext_field
                end do
                call polarization(sys_obj, ef, solver)
                call mfree('ommp_get_polelec_energy [ef]', ef)
            else
                call mallocate('ommp_get_polelec_energy [ef]', &
                               3, eel%pol_atoms, eel%n_ipd, ef)
                
                ef(:,:,1) = ext_field
                call polarization(sys_obj, ef, solver, &
                                  OMMP_MATV_DEFAULT, [.true., .false.] )
                
                call mfree('ommp_get_polelec_energy [ef]', ef)
            end if
        end subroutine ommp_set_external_field

        function ommp_get_polelec_energy(sys_obj) result(ene)
            
            use mod_electrostatics, only: energy_MM_pol, prepare_M2D
            use mod_polarization, only: polarization
            use mod_constants, only: OMMP_SOLVER_DEFAULT

            implicit none
            
            type(ommp_system), intent(inout), target :: sys_obj
            real(ommp_real) :: ene

            if(.not. sys_obj%eel%ipd_done) then
                !! Solve the polarization system without external field
                call prepare_M2D(sys_obj%eel)
                call polarization(sys_obj, sys_obj%eel%e_m2d, OMMP_SOLVER_DEFAULT)
            end if

            ene = 0.0
            call energy_MM_pol(sys_obj%eel, ene)
        end function
        
        function ommp_get_fixedelec_energy(sys_obj) result(ene)
            
            use mod_electrostatics, only: energy_MM_MM

            implicit none
            
            type(ommp_system), intent(inout), target :: sys_obj
            real(ommp_real) :: ene

            ene = 0.0
            call energy_MM_MM(sys_obj%eel, ene)

        end function
        
        function ommp_get_vdw_energy(sys_obj) result(evdw)
            
            use mod_nonbonded, only: vdw_potential
            
            implicit none
            type(ommp_system), intent(inout), target :: sys_obj
            real(ommp_real) :: evdw

            evdw = 0.0
            call vdw_potential(sys_obj%vdw, evdw)
        
        end function
        
end module ommp_interface

