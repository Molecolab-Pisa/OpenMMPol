#include "f_cart_components.h"
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

    ! Renamed import of several global variables that should be available
    ! in the interface
    use mod_constants, only: OMMP_FF_AMOEBA, OMMP_FF_WANG_AL, OMMP_FF_WANG_DL, &
                             OMMP_SOLVER_CG, OMMP_SOLVER_DIIS, &
                             OMMP_SOLVER_INVERSION, OMMP_SOLVER_DEFAULT, &
                             OMMP_MATV_INCORE, OMMP_MATV_DIRECT, &
                             OMMP_MATV_DEFAULT, &
                             OMMP_VERBOSE_DEBUG, OMMP_VERBOSE_HIGH, &
                             OMMP_VERBOSE_LOW, OMMP_VERBOSE_NONE, &
                             OMMP_AU2KCALMOL => au2kcalmol, &
                             OMMP_ANG2AU => angstrom2au
    
    ! Internal types
    use mod_memory, only: ommp_integer => ip, &
                          ommp_real => rp
    use mod_mmpol, only: ommp_system
    use mod_electrostatics, only: ommp_electrostatics_type
    use mod_topology, only: ommp_topology_type
    use mod_qm_helper, only: ommp_qm_helper

    use mod_mmpol, only: ommp_save_mmp => mmpol_save_as_mmp, &
                         ommp_print_summary => mmpol_ommp_print_summary, &
                         ommp_update_coordinates => update_coordinates, &
                         ommp_print_summary_to_file => mmpol_ommp_print_summary

    use mod_io, only: ommp_set_verbose => set_verbosity, ommp_version
    
    use mod_qm_helper, only: ommp_qm_helper_init_vdw_prm => qm_helper_init_vdw_prm, &
                             ommp_qm_helper_init_vdw => qm_helper_init_vdw, &
                             ommp_prepare_qm_ele_ene => electrostatic_for_ene, &
                             ommp_prepare_qm_ele_grd => electrostatic_for_grad

    implicit none
    
    contains
        
        subroutine ommp_init_mmp(s, filename)
            use mod_inputloader, only : mmpol_init_from_mmp
            
            implicit none
            
            type(ommp_system), pointer, intent(inout) :: s
            character(len=*) :: filename

            call ommp_version(OMMP_VERBOSE_LOW)
            allocate(s)
            call mmpol_init_from_mmp(trim(filename), s)
        end subroutine
        
        subroutine ommp_init_xyz(s, xyzfile, prmfile)
            use mod_inputloader, only : mmpol_init_from_xyz
            
            implicit none
            
            type(ommp_system), pointer, intent(inout) :: s
            character(len=*) :: xyzfile, prmfile

            call ommp_version(OMMP_VERBOSE_LOW)
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
            use mod_electrostatics, only: prepare_polelec
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
                call prepare_polelec(eel)
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

        subroutine ommp_set_external_field_nomm(sys_obj, ext_field, solver)
            !! This is just the same as [[ommp_set_external_field]] but 
            !! implicitly assuming [[ommp_set_external_field:add_mm_field]] as 
            !! false, mainly here for interface consistency with C

            implicit none
            
            type(ommp_system), intent(inout), target :: sys_obj
            real(ommp_real), intent(in) :: ext_field(3,sys_obj%eel%pol_atoms)
            integer(ommp_integer), intent(in), value :: solver

            call ommp_set_external_field(sys_obj, ext_field, solver, .false.)
        end subroutine
        
        subroutine ommp_potential_mmpol2ext(s, n, cext, v)
            ! Compute the electric potential of static sites at
            ! arbitrary coordinates
            use mod_electrostatics, only: potential_D2E, &
                                          potential_M2E

            implicit none
            
            type(ommp_system), intent(inout), target :: s
            integer(ommp_integer), intent(in) :: n
            real(ommp_real), intent(in) :: cext(3,n)
            real(ommp_real), intent(inout) :: v(n)
            
            call potential_M2E(s%eel, cext, v)
            call potential_D2E(s%eel, cext, v)
        end subroutine
        
        subroutine ommp_potential_pol2ext(s, n, cext, v) 
            ! Compute the electric potential of static sites at
            ! arbitrary coordinates
            use mod_electrostatics, only: potential_D2E

            implicit none
            
            type(ommp_system), intent(inout), target :: s
            integer(ommp_integer), intent(in) :: n
            real(ommp_real), intent(in) :: cext(3,n)
            real(ommp_real), intent(inout) :: v(n)
            
            call potential_D2E(s%eel, cext, v)
        end subroutine
        
        subroutine ommp_potential_mm2ext(s, n, cext, v)
            ! Compute the electric potential of static sites at
            ! arbitrary coordinates
            use mod_electrostatics, only: potential_M2E

            implicit none
            
            type(ommp_system), intent(inout), target :: s
            integer(ommp_integer), intent(in) :: n
            real(ommp_real), intent(in) :: cext(3,n)
            real(ommp_real), intent(inout) :: v(n)
            
            call potential_M2E(s%eel, cext, v)
        end subroutine

        function ommp_get_polelec_energy(sys_obj) result(ene)
            !! Solve the polarization equation for a certain external field
            !! and compute the interaction energy of the induced dipoles with
            !! themselves and fixed multipoles.

            use mod_electrostatics, only: energy_MM_pol, prepare_polelec
            use mod_polarization, only: polarization
            use mod_constants, only: OMMP_SOLVER_DEFAULT

            implicit none
            
            type(ommp_system), intent(inout), target :: sys_obj
            real(ommp_real) :: ene

            if(.not. sys_obj%eel%ipd_done) then
                !! Solve the polarization system without external field
                call prepare_polelec(sys_obj%eel)
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
        
        function ommp_get_full_ele_energy(sys_obj) result(ene)

            implicit none
            
            type(ommp_system), intent(inout), target :: sys_obj
            real(ommp_real) :: ene

            ene = ommp_get_fixedelec_energy(sys_obj)
            ene = ene + ommp_get_polelec_energy(sys_obj)

        end function
        
        function ommp_get_vdw_energy(sys_obj) result(evdw)
            
            use mod_nonbonded, only: vdw_potential
            
            implicit none
            type(ommp_system), intent(inout), target :: sys_obj
            real(ommp_real) :: evdw

            evdw = 0.0
            if(sys_obj%use_nonbonded) call vdw_potential(sys_obj%vdw, evdw)
        
        end function
        
        function ommp_get_bond_energy(sys_obj) result(eb)
            
            use mod_bonded, only: bond_potential
            
            implicit none
            type(ommp_system), intent(inout), target :: sys_obj
            real(ommp_real) :: eb

            eb = 0.0
            if(sys_obj%use_bonded) call bond_potential(sys_obj%bds, eb)
        
        end function
        
        function ommp_get_angle_energy(sys_obj) result(ea)
            
            use mod_bonded, only: angle_potential
            
            implicit none
            type(ommp_system), intent(inout), target :: sys_obj
            real(ommp_real) :: ea

            ea = 0.0
            if(sys_obj%use_bonded) call angle_potential(sys_obj%bds, ea)
        
        end function
        
        function ommp_get_strbnd_energy(sys_obj) result(eba)
            
            use mod_bonded, only: strbnd_potential
            
            implicit none
            type(ommp_system), intent(inout), target :: sys_obj
            real(ommp_real) :: eba

            eba = 0.0
            if(sys_obj%use_bonded) call strbnd_potential(sys_obj%bds, eba)
        
        end function
        
        function ommp_get_urey_energy(sys_obj) result(eub)
            
            use mod_bonded, only: urey_potential
            
            implicit none
            type(ommp_system), intent(inout), target :: sys_obj
            real(ommp_real) :: eub

            eub = 0.0
            if(sys_obj%use_bonded) call urey_potential(sys_obj%bds, eub)
        
        end function
        
        function ommp_get_opb_energy(sys_obj) result(eopb)
            
            use mod_bonded, only: opb_potential
            
            implicit none
            type(ommp_system), intent(inout), target :: sys_obj
            real(ommp_real) :: eopb

            eopb = 0.0
            if(sys_obj%use_bonded) call opb_potential(sys_obj%bds, eopb)
        
        end function
        
        function ommp_get_imptorsion_energy(sys_obj) result(et)
            
            use mod_bonded, only: imptorsion_potential
            
            implicit none
            type(ommp_system), intent(inout), target :: sys_obj
            real(ommp_real) :: et

            et = 0.0
            if(sys_obj%use_bonded) call imptorsion_potential(sys_obj%bds, et)
        
        end function
        
        function ommp_get_torsion_energy(sys_obj) result(et)
            
            use mod_bonded, only: torsion_potential
            
            implicit none
            type(ommp_system), intent(inout), target :: sys_obj
            real(ommp_real) :: et

            et = 0.0
            if(sys_obj%use_bonded) call torsion_potential(sys_obj%bds, et)
        
        end function
        
        function ommp_get_pitors_energy(sys_obj) result(ept)
            
            use mod_bonded, only: pitors_potential
            
            implicit none
            type(ommp_system), intent(inout), target :: sys_obj
            real(ommp_real) :: ept

            ept = 0.0
            if(sys_obj%use_bonded) call pitors_potential(sys_obj%bds, ept)
        
        end function
        
        function ommp_get_strtor_energy(sys_obj) result(est)
            
            use mod_bonded, only: strtor_potential
            
            implicit none
            type(ommp_system), intent(inout), target :: sys_obj
            real(ommp_real) :: est

            est = 0.0
            if(sys_obj%use_bonded) call strtor_potential(sys_obj%bds, est)
        
        end function

        function ommp_get_angtor_energy(sys_obj) result(eat)
         
            use mod_bonded, only: angtor_potential
             
            implicit none
            type(ommp_system), intent(inout), target :: sys_obj
            real(ommp_real) :: eat

            eat = 0.0
            if(sys_obj%use_bonded) call angtor_potential(sys_obj%bds, eat)

        end function
        
        function ommp_get_tortor_energy(sys_obj) result(ett)
            
            use mod_bonded, only: tortor_potential
            
            implicit none
            type(ommp_system), intent(inout), target :: sys_obj
            real(ommp_real) :: ett

            ett = 0.0
            if(sys_obj%use_bonded) call tortor_potential(sys_obj%bds, ett)
        
        end function
        
        function ommp_get_full_bnd_energy(sys_obj) result(ene)
            
            use mod_bonded, only: bond_potential, angtor_potential, &
                                  strbnd_potential, urey_potential, &
                                  opb_potential, pitors_potential, &
                                  torsion_potential, tortor_potential, &
                                  strtor_potential, angle_potential, &
                                  imptorsion_potential
            
            implicit none
            type(ommp_system), intent(inout), target :: sys_obj
            real(ommp_real) :: ene

            ene = 0.0
            
            if(sys_obj%use_bonded) then
                call bond_potential(sys_obj%bds, ene)
                call angle_potential(sys_obj%bds, ene)
                call strbnd_potential(sys_obj%bds, ene)
                call urey_potential(sys_obj%bds, ene)
                call opb_potential(sys_obj%bds, ene)
                call imptorsion_potential(sys_obj%bds, ene) 
                call torsion_potential(sys_obj%bds, ene)
                call pitors_potential(sys_obj%bds, ene)
                call strtor_potential(sys_obj%bds, ene)
                call angtor_potential(sys_obj%bds, ene)
                call tortor_potential(sys_obj%bds, ene)
            end if
        end function
        
        function ommp_get_full_energy(sys_obj) result(ene)
            
            implicit none
            type(ommp_system), intent(inout), target :: sys_obj
            real(ommp_real) :: ene

            ene = ommp_get_vdw_energy(sys_obj)
            ene = ene + ommp_get_full_ele_energy(sys_obj)
            ene = ene + ommp_get_full_bnd_energy(sys_obj)
        end function

        ! Functions for advanced operation and gradients
        subroutine ommp_fixedelec_geomgrad(s, grd)
            use mod_geomgrad, only: fixedelec_geomgrad
            
            implicit none 
            
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(out) :: grd(3,s%top%mm_atoms)

            grd = 0.0
            call fixedelec_geomgrad(s, grd)
        end subroutine
        
        subroutine ommp_polelec_geomgrad(s, grd)
            use mod_geomgrad, only: polelec_geomgrad
            
            implicit none 
            
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(out) :: grd(3,s%top%mm_atoms)

            grd = 0.0
            call polelec_geomgrad(s, grd)
        end subroutine

        subroutine ommp_vdw_geomgrad(s, grd)
            use mod_nonbonded, only: vdw_geomgrad
            
            implicit none 
            
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(out) :: grd(3,s%top%mm_atoms)

            grd = 0.0
            if(s%use_nonbonded) call vdw_geomgrad(s%vdw, grd)
        end subroutine
        
        subroutine ommp_rotation_geomgrad(s, E, Egrd, grd )
            implicit none

            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(in) :: E(:,:), Egrd(:,:)
            real(ommp_real), intent(out) :: grd(:,:)
            
            grd = 0.0
            call rotation_geomgrad(s%eel, E, Egrd, grd)
        end subroutine

        subroutine ommp_bond_geomgrad(s, grd)
            use mod_bonded, only: bond_geomgrad 
            
            implicit none 
        
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(out) :: grd(3,s%top%mm_atoms)

            grd = 0.0
            if(s%use_bonded) call bond_geomgrad(s%bds, grd)
        end subroutine
        
        subroutine ommp_angle_geomgrad(s, grd)
            use mod_bonded, only: angle_geomgrad 
            
            implicit none 
            
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(out) :: grd(3,s%top%mm_atoms)

            grd = 0.0
            if(s%use_bonded) call angle_geomgrad(s%bds, grd)
        end subroutine
        
        subroutine ommp_strbnd_geomgrad(s, grd)
            use mod_bonded, only: strbnd_geomgrad 
            
            implicit none 
            
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(out) :: grd(3,s%top%mm_atoms)

            grd = 0.0
            if(s%use_bonded) call strbnd_geomgrad(s%bds, grd)
        end subroutine
        
        subroutine ommp_urey_geomgrad(s, grd)
            use mod_bonded, only: urey_geomgrad 
            
            implicit none 
            
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(out) :: grd(3,s%top%mm_atoms)

            grd = 0.0
            if(s%use_bonded) call urey_geomgrad(s%bds, grd)
        end subroutine
        
        subroutine ommp_torsion_geomgrad(s, grd)
            use mod_bonded, only: torsion_geomgrad 
            
            implicit none 
            
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(out) :: grd(3,s%top%mm_atoms)

            grd = 0.0
            if(s%use_bonded) call torsion_geomgrad(s%bds, grd)
        end subroutine
        
        subroutine ommp_imptorsion_geomgrad(s, grd)
            use mod_bonded, only: imptorsion_geomgrad 
            
            implicit none 
            
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(out) :: grd(3,s%top%mm_atoms)

            grd = 0.0
            if(s%use_bonded) call imptorsion_geomgrad(s%bds, grd)
        end subroutine
        
        subroutine ommp_angtor_geomgrad(s, grd)
            use mod_bonded, only: angtor_geomgrad 
            
            implicit none 
            
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(out) :: grd(3,s%top%mm_atoms)

            grd = 0.0
            if(s%use_bonded) call angtor_geomgrad(s%bds, grd)
        end subroutine
        
        subroutine ommp_opb_geomgrad(s, grd)
            use mod_bonded, only: opb_geomgrad 
            
            implicit none 
            
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(out) :: grd(3,s%top%mm_atoms)

            grd = 0.0
            if(s%use_bonded) call opb_geomgrad(s%bds, grd)
        end subroutine
        
        subroutine ommp_strtor_geomgrad(s, grd)
            use mod_bonded, only: strtor_geomgrad 
            
            implicit none 
            
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(out) :: grd(3,s%top%mm_atoms)

            grd = 0.0
            if(s%use_bonded) call strtor_geomgrad(s%bds, grd)
        end subroutine
        
        subroutine ommp_tortor_geomgrad(s, grd)
            use mod_bonded, only: tortor_geomgrad 
            
            implicit none 
            
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(out) :: grd(3,s%top%mm_atoms)

            grd = 0.0
            if(s%use_bonded) call tortor_geomgrad(s%bds, grd)
        end subroutine
        
        subroutine ommp_pitors_geomgrad(s, grd)
            use mod_bonded, only: pitors_geomgrad 
            
            implicit none 
            
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(out) :: grd(3,s%top%mm_atoms)

            grd = 0.0
            if(s%use_bonded) call pitors_geomgrad(s%bds, grd)
        end subroutine
        
        subroutine ommp_full_bnd_geomgrad(s, grd)
            use mod_bonded, only: bond_geomgrad, &
                                  angle_geomgrad, &
                                  strbnd_geomgrad, &
                                  urey_geomgrad, &
                                  opb_geomgrad, &
                                  imptorsion_geomgrad, &
                                  torsion_geomgrad, &
                                  pitors_geomgrad, &
                                  strtor_geomgrad, &
                                  angtor_geomgrad, &
                                  tortor_geomgrad
            
            implicit none 
            
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(out) :: grd(3,s%top%mm_atoms)

            grd = 0.0
            if(s%use_bonded) then
                call bond_geomgrad(s%bds, grd)
                call angle_geomgrad(s%bds, grd)
                call strbnd_geomgrad(s%bds, grd)
                call urey_geomgrad(s%bds, grd)
                call opb_geomgrad(s%bds, grd)
                call imptorsion_geomgrad(s%bds, grd) 
                call torsion_geomgrad(s%bds, grd)
                call pitors_geomgrad(s%bds, grd)
                call strtor_geomgrad(s%bds, grd)
                call angtor_geomgrad(s%bds, grd)
                call tortor_geomgrad(s%bds, grd)
            end if
        end subroutine
        
        subroutine ommp_full_geomgrad(s, grd)
            use mod_memory, only: mallocate
            use mod_nonbonded, only: vdw_geomgrad
            use mod_geomgrad, only: polelec_geomgrad, fixedelec_geomgrad

            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(out) :: grd(3,s%top%mm_atoms)

            grd = 0.0
            call ommp_full_bnd_geomgrad(s, grd)
            call fixedelec_geomgrad(s, grd)
            call polelec_geomgrad(s, grd)
            if(s%use_nonbonded) call vdw_geomgrad(s%vdw, grd)

        end subroutine

#ifdef USE_HDF5
        subroutine ommp_init_hdf5(s, filename, namespace)
            !! This function is an interface for saving an HDF5 file 
            !! with all the data contained in mmpol module using
            !! [[mod_io:mmpol_save_as_hdf5]]
            use mod_iohdf5, only: mmpol_init_from_hdf5
            
            implicit none
            
            type(ommp_system), pointer :: s
            character(len=*) :: filename, namespace
            integer(ommp_integer) :: ok

            call ommp_version(OMMP_VERBOSE_LOW)
            allocate(s)
            call mmpol_init_from_hdf5(filename, namespace, s, ok)
            
        end subroutine ommp_init_hdf5
        
        subroutine ommp_save_as_hdf5(s, filename, namespace) 
            
            use mod_iohdf5, only: save_system_as_hdf5 

            implicit none
            
            character(len=*) :: filename, namespace
            type(ommp_system), pointer :: s
            integer(kind=4) :: err

            call save_system_as_hdf5(filename, s, err, namespace, .false.)
            
        end subroutine ommp_save_as_hdf5
        
        subroutine ommp_checkpoint(s, filename, namespace)
            
            use mod_iohdf5, only: save_system_as_hdf5 

            implicit none
            
            character(len=*) :: filename, namespace
            type(ommp_system), pointer :: s
            integer(kind=4) :: err

            call save_system_as_hdf5(filename, s, err, namespace, .true.)
            
        end subroutine ommp_checkpoint
#endif

    ! QM Helper Object housekeeping
    subroutine ommp_init_qm_helper(s, n, cqm, qqm, zqm)
        
        use mod_qm_helper, only: qm_helper_init
        
        implicit none

        type(ommp_qm_helper), pointer, intent(inout) :: s
        integer(ommp_integer) :: n
        real(ommp_real), intent(in) :: cqm(:,:), qqm(:)
        integer(ommp_integer), intent(in) :: zqm(:)

        allocate(s)
        call qm_helper_init(s, n, cqm, qqm, zqm)
    end subroutine
    
    subroutine ommp_terminate_qm_helper(s) 
        
        use mod_qm_helper, only: qm_helper_terminate
        
        implicit none

        type(ommp_qm_helper), pointer, intent(inout) :: s
        
        call qm_helper_terminate(s)
        deallocate(s)
    end subroutine
    
    function ommp_qm_helper_vdw_energy(qm, s) result(evdw)
        use mod_qm_helper, only: qm_helper_vdw_energy

        implicit none

        type(ommp_system), intent(in) :: s
        type(ommp_qm_helper), intent(in) :: qm
        real(ommp_real) :: evdw

        evdw = 0.0
        call qm_helper_vdw_energy(qm, s, evdw)
    end function
    
    subroutine ommp_qm_helper_vdw_geomgrad(qm, s, qmg, mmg)
        
        use mod_qm_helper, only: qm_helper_vdw_geomgrad

        implicit none

        type(ommp_system), intent(in) :: s
        type(ommp_qm_helper), intent(in) :: qm
        real(ommp_real), intent(out) :: qmg(:,:), mmg(:,:)

        mmg = 0.0
        qmg = 0.0
        call qm_helper_vdw_geomgrad(qm, s, qmg, mmg)
    end subroutine

end module ommp_interface

