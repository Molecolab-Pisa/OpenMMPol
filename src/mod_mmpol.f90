module mod_mmpol
    !! Main module for the control of openMMPol library. It contains
    !! all the scalar and vector (allocatable) quantities needed to
    !! build up the atomistic polarizable embedding model and perform
    !! the calculation required from the quantum chemical software.
    
    use mod_memory, only: ip, rp, lp
    use mod_adjacency_mat, only: yale_sparse
    use mod_topology, only: ommp_topology_type, topology_init, &
                            topology_terminate
    use mod_electrostatics, only: ommp_electrostatics_type
    use mod_nonbonded, only: ommp_nonbonded_type
    use mod_bonded, only: ommp_bonded_type
    use mod_link_atom, only: ommp_link_atom_type
    use mod_io, only: ommp_message, fatal_error
    use mod_constants, only: OMMP_STR_CHAR_MAX

    implicit none 

    type ommp_system
        logical :: mmpol_is_init = .false.
        !! Initialization flag
        integer(ip) :: ff_type
        !! Force field type selection flag (0 for AMBER, 1 for AMOEBA)
        logical(lp) :: amoeba
        !! AMOEBA FF = True; WANG-AMBER = False

        type(ommp_topology_type), allocatable :: top
        !! Data structure containing the topology of the system
        type(ommp_electrostatics_type), allocatable :: eel
        !! Data structure containing all the information needed to run the
        !! elctrostatics related calculations
        logical(lp) :: use_bonded = .false.
        type(ommp_bonded_type), allocatable :: bds
        !! Data structure containing all the information needed to run the
        !! bonded terms calculations
        logical(lp) :: use_nonbonded = .false.
        type(ommp_nonbonded_type), allocatable :: vdw
        !! Data structure containing all the information needed to run the
        !! non-bonded terms calculations
        logical(lp) :: use_linkatoms = .false.
        type(ommp_link_atom_type), allocatable :: la
        !! Data structure containing all the information needed to handle 
        !! link atoms with a certain QM part described by a QM Helper object
    end type ommp_system
    
    contains

    subroutine mmpol_init(sys_obj, l_ff_type, l_mm_atoms, l_pol_atoms)
        !! Performs all the memory allocation and vector initialization
        !! needed to run the openMMPol library
        
        use mod_constants, only: OMMP_FF_AMBER, OMMP_FF_AMOEBA, OMMP_FF_UNKNOWN
        use mod_electrostatics, only: electrostatics_init
        use mod_io, only: print_matrix, fatal_error
        use mod_profiling, only: time_push, time_pull

        implicit none

        type(ommp_system), intent(inout) :: sys_obj
        !! The object to be initialized
        integer(ip), intent(in) :: l_ff_type
        !! Force field type used in initialization
        integer(ip), intent(in) :: l_mm_atoms
        !! Number of MM atoms used in initialization
        integer(ip), intent(in) :: l_pol_atoms
        !! Number of polarizable atoms used in initialization
        
        call time_push()

        !! Allocation topology...
        allocate(sys_obj%top)
        !! ... and electrostatics
        allocate(sys_obj%eel)
        
        ! FF related settings
        if(l_ff_type == OMMP_FF_UNKNOWN) then
            call fatal_error("Cannot initialize an UNKNOWN forcefield!")
        end if
        sys_obj%ff_type = l_ff_type
        
        if(sys_obj%ff_type == OMMP_FF_AMOEBA) then
            sys_obj%amoeba = .true.
        else if(sys_obj%ff_type == OMMP_FF_AMBER) then
            sys_obj%amoeba = .false.
        end if
        
        ! Initialization of sub-modules:
        !   a. topology
        call topology_init(sys_obj%top, l_mm_atoms)
        !   b. electrostatics
        call electrostatics_init(sys_obj%eel, sys_obj%amoeba, l_pol_atoms, &
                                 sys_obj%top)

        ! Everything is done
        sys_obj%mmpol_is_init = .true.
        call time_pull('MMPol object creation (mmpol_init)')
    end subroutine mmpol_init

    subroutine mmpol_init_nonbonded(sys_obj)
        !! Enable nonbonded part of pontential
        implicit none

        type(ommp_system), intent(inout) :: sys_obj
        !! The object to be initialized

        allocate(sys_obj%vdw)
        sys_obj%use_nonbonded = .true.

    end subroutine mmpol_init_nonbonded
    
    subroutine mmpol_init_bonded(sys_obj)
        !! Enable nonbonded part of pontential
        implicit none

        type(ommp_system), intent(inout), target :: sys_obj
        !! The object to be initialized

        allocate(sys_obj%bds)
        sys_obj%use_bonded = .true.
        sys_obj%bds%top => sys_obj%top

    end subroutine mmpol_init_bonded

    subroutine mmpol_init_link_atom(sys_obj)
        !! Enable link atom
        implicit none

        type(ommp_system), intent(inout), target :: sys_obj
        !! The object to be initialized
        
        if(sys_obj%use_linkatoms .and. allocated(sys_obj%la)) return
        if(allocated(sys_obj%la)) deallocate(sys_obj%la)

        allocate(sys_obj%la)
        sys_obj%use_linkatoms = .true.
    end subroutine
        
    subroutine mmpol_prepare(sys_obj)
        !! Compute some derived quantities from the input that 
        !! are used during the calculation. The upstream code have
        !! to provide cmm, q, pol, adjacency matrix and in
        !! the case of AMOEBA also multipoles rotation information, and 
        !! polarization group information.   
        !! This routine    
        !!   * compute connectivity lists from connected atoms    
        !!   * invert polar_mm list creating mm_polar   
        !!   * populate cpol list of coordinates   
        !!   * compute factors for thole damping    
        !!   * scales by 1/3 AMOEBA quadrupoles (?)    
        !!   * Build list for polarization groups, compute groups connectivity   
        !!   * performs multipoles rotation   

        use mod_adjacency_mat, only: build_conn_upto_n, matcpy, reverse_grp_tab
        use mod_io, only: ommp_message
        use mod_profiling, only: time_push, time_pull
        use mod_constants, only: OMMP_VERBOSE_DEBUG
        use mod_electrostatics, only: thole_init, remove_null_pol, &
                                      make_screening_lists, fmm_coordinates_update

        implicit none
        
        type(ommp_system), intent(inout) :: sys_obj
        !! The system object to bi initialized

        integer(ip) :: i
        
        type(yale_sparse) :: adj, pg_adj

        call time_push()
        call ommp_message("Building connectivity lists", OMMP_VERBOSE_DEBUG)
        
        ! compute connectivity lists from connected atoms
        if(size(sys_obj%top%conn) < 4) then
            call matcpy(sys_obj%top%conn(1), adj)
            deallocate(sys_obj%top%conn)
            call build_conn_upto_n(adj, 4, sys_obj%top%conn, .false.)
        end if

        call remove_null_pol(sys_obj%eel)
       
        ! invert mm_polar list creating mm_polar
        sys_obj%eel%mm_polar(:) = 0
        if(sys_obj%eel%pol_atoms > 0) then
            call ommp_message("Creating MM->polar and polar->MM lists", &
                              OMMP_VERBOSE_DEBUG)
            !$omp parallel do default(shared) private(i)
            do i = 1, sys_obj%eel%pol_atoms
                sys_obj%eel%mm_polar(sys_obj%eel%polar_mm(i)) = i
            end do

            call ommp_message("Populating coordinates of polarizable atoms", &
                              OMMP_VERBOSE_DEBUG)
            ! populate cpol list of coordinates
            !$omp parallel do default(shared) private(i)
            do i = 1, sys_obj%eel%pol_atoms
                sys_obj%eel%cpol(:,i) = sys_obj%top%cmm(:, sys_obj%eel%polar_mm(i))
            end do

            call ommp_message("Setting Thole factors", OMMP_VERBOSE_DEBUG)
            ! compute factors for thole damping
            call thole_init(sys_obj%eel)
        else
            sys_obj%eel%thole = 0.0
        end if

        if(sys_obj%amoeba) then
            ! Copy multipoles from q to q0
            sys_obj%eel%q0 = sys_obj%eel%q

            ! scales by 1/3 AMOEBA quadrupoles (?)
            ! Mysterious division of multipoles by three
            ! FL told me that it was done like that in
            ! Tinker
            sys_obj%eel%q0(5:10,:) = sys_obj%eel%q0(5:10,:) / 3.0_rp

            ! polarization groups connectivity list
            call reverse_grp_tab(sys_obj%eel%mmat_polgrp, &
                                 sys_obj%eel%polgrp_mmat)
            call build_pg_adjacency_matrix(sys_obj%eel, pg_adj)
            call build_conn_upto_n(pg_adj, 3, sys_obj%eel%pg_conn, .true.)

            ! performs multipoles rotation
            call rotate_multipoles(sys_obj%eel)
        end if

        if(sys_obj%eel%use_fmm) then
            call ommp_message("Building FMM near lists", OMMP_VERBOSE_DEBUG)
            call fmm_coordinates_update(sys_obj%eel)
        end if
        call ommp_message("Building screening lists", OMMP_VERBOSE_DEBUG)
        call time_push()
        call make_screening_lists(sys_obj%eel)
        call time_pull("Preparing screening lists")
        call ommp_message("MMPol initialization (mmpol_prepare) completed.", OMMP_VERBOSE_DEBUG)
        call time_pull('MMPol object initialization (mmpol_prepare)')
    end subroutine mmpol_prepare

    subroutine mmpol_terminate(sys_obj)
        !! Performs all the deallocation needed at the end of the 
        !! calculation
        use mod_memory, only: mfree
        use mod_electrostatics, only: electrostatics_terminate
        use mod_nonbonded, only: vdw_terminate
        use mod_bonded, only: bonded_terminate 

        implicit none 

        type(ommp_system), intent(inout) :: sys_obj

        call electrostatics_terminate(sys_obj%eel)
        deallocate(sys_obj%eel)

        call topology_terminate(sys_obj%top)
        deallocate(sys_obj%top)

        if(sys_obj%use_nonbonded) then
            call vdw_terminate(sys_obj%vdw)
            sys_obj%use_nonbonded = .false.
        end if
        
        if(sys_obj%use_bonded) then
            call bonded_terminate(sys_obj%bds)
            sys_obj%use_bonded = .false.
        end if

        sys_obj%mmpol_is_init = .false.

    end subroutine mmpol_terminate
    
    !TODO move to eel module
    subroutine build_pg_adjacency_matrix(eel, adj)
        !! Builds the adjacency matrix of polarization groups starting from
        !! atomic adjacency matrix and list of polarization groups indices.

        use mod_adjacency_mat, only: reallocate_mat 

        implicit none

        type(ommp_electrostatics_type), intent(in) :: eel
        type(yale_sparse), intent(out) :: adj
        !! The group adjacency matrix to be saved.

        integer(ip) :: npg, pg1, atm1, atm2, i, j

        npg = eel%polgrp_mmat%n

        adj%n = npg
        allocate(adj%ri(adj%n+1))
        allocate(adj%ci(adj%n*2))
        adj%ri(1) = 1

        do pg1=1, npg
            ! For each polarization group
            adj%ri(pg1+1) = adj%ri(pg1)

            do i=eel%polgrp_mmat%ri(pg1), eel%polgrp_mmat%ri(pg1+1)-1
                ! Loop on every atom of the group
                atm1 = eel%polgrp_mmat%ci(i)
                do j=eel%top%conn(1)%ri(atm1), eel%top%conn(1)%ri(atm1+1)-1
                    ! Loop on each connected atom...
                    atm2 = eel%top%conn(1)%ci(j)

                    ! If the two atoms are in different PG, then the two
                    ! polarization groups are connected. 
                    if(eel%mmat_polgrp(atm1) /= eel%mmat_polgrp(atm2) .and. &
                       ! if the group is not already present in the matrix
                       all(adj%ci(adj%ri(pg1):adj%ri(pg1+1)-1) /= eel%mmat_polgrp(atm2))) then
                        adj%ci(adj%ri(pg1+1)) = eel%mmat_polgrp(atm2)
                        adj%ri(pg1+1) = adj%ri(pg1+1) + 1
                        if(adj%ri(pg1+1) > size(adj%ci)) then
                            ! If matrix is too small, it could be enlarged...
                            call reallocate_mat(adj, size(adj%ci)+adj%n)
                        end if
                    end if
                end do
            end do
        end do
        
        ! Finally trim the output matrix
        call reallocate_mat(adj, adj%ri(adj%n+1)-1)

    end subroutine build_pg_adjacency_matrix
    
    subroutine update_coordinates(sys_obj, new_c) 
        !! Interface to change the coordinates of the system (eg. during a 
        !! MD or a geometry optimization). This function clears all the 
        !! relevant, flags and update the needed quantities. All those 
        !! operations are needed for a correct functionality of the program 
        !! therefore coordinates should never be updated without passing from
        !! this interface.
       
        use mod_memory, only: mfree
        use mod_link_atom, only: link_atom_update_merged_topology
        use mod_electrostatics, only: fmm_coordinates_update
        implicit none

        type(ommp_system), intent(inout), target :: sys_obj
        !! System data structure
        real(rp), dimension(3,sys_obj%top%mm_atoms), intent(in) :: new_c
        !! New coordinates to be updated

        type(ommp_topology_type), pointer :: top
        type(ommp_electrostatics_type), pointer :: eel
        integer(ip) :: i

        top => sys_obj%top
        eel => sys_obj%eel

        ! 1. Copy coordinates
        top%cmm = new_c

        ! 2. Update electrostatics module
        ! 2.1 Coordinates
        do i=1, eel%pol_atoms
            eel%cpol(:,i) = top%cmm(:,eel%polar_mm(i))
        end do
        ! 2.2 Flags and allocated quantities
        eel%M2M_done = .false.
        eel%M2Mgg_done = .false.
        eel%M2D_done = .false.
        eel%M2Dgg_done = .false.
        eel%ipd_done = .false.
        eel%ipd_use_guess = .false.
        if(allocated(eel%TMat)) call mfree('update_coordinates [TMat]',eel%TMat)
        ! 2.3 Multipoles rotation
        if(sys_obj%amoeba) call rotate_multipoles(sys_obj%eel)
        ! 2.3 Update coordinates inside link atom object
        if(sys_obj%use_linkatoms) call link_atom_update_merged_topology(sys_obj%la)
        ! 2.4 Update fast-multipoles tree if needed
        if(eel%use_fmm) call fmm_coordinates_update(eel)
    end subroutine
        
    
    subroutine mmpol_ommp_print_summary(sys_obj, of_name)
        !! Prints a complete summary of all the quantities stored 
        !! in the MMPol module
        use mod_memory, only: mallocate, mfree
        use mod_io, only: iof_mmpol, print_matrix, print_int_vec
        use mod_utils, only: sort_ivec
        
        implicit none

        type(ommp_system), intent(in) :: sys_obj
        character(len=*), intent(in), optional :: of_name
        
        integer(ip) :: of_unit

        integer(ip) :: i, j, grp, igrp, lst(1000), ilst
        real(rp), allocatable :: polar(:) ! Polarizabilities of all atoms
        integer(ip), allocatable :: tmp(:)
        character(len=OMMP_STR_CHAR_MAX) :: str

        if(present(of_name)) then
            of_unit = 101
            open(unit=of_unit, &
                 file=of_name(1:len(trim(of_name))), &
                 action='write')
        else
            of_unit = iof_mmpol
        end if

        call mallocate('mmpol_ommp_print_summary [polar]', &
                       sys_obj%top%mm_atoms, polar)
        polar = 0.0_rp
        do i=1, sys_obj%eel%pol_atoms
            polar(sys_obj%eel%polar_mm(i)) = sys_obj%eel%pol(i)
        end do

        write(of_unit, '(A, 4F8.4)') 'mscale: ', sys_obj%eel%mscale
        if(sys_obj%eel%pol_atoms > 0) then
            write(of_unit, '(A, 4F8.4)') 'pscale: ', sys_obj%eel%pscale
            if(sys_obj%amoeba) write(of_unit, '(A, 4F8.4)') 'pscale (intra): ', &
                                                          sys_obj%eel%pscale_intra
            write(of_unit, '(A, 4F8.4)') 'dscale: ', sys_obj%eel%dscale
            write(of_unit, '(A, 4F8.4)') 'uscale: ', sys_obj%eel%uscale
        end if

        call print_matrix(.true., 'coordinates:', sys_obj%top%cmm, of_unit)
        if(sys_obj%amoeba) then
            call print_matrix(.true., 'multipoles - non rotated:', &
                              sys_obj%eel%q0, of_unit)
        end if
        call print_matrix(.true., 'multipoles :', sys_obj%eel%q, of_unit)
        call print_matrix(.true., 'coordinates of polarizable atoms:', &
                          sys_obj%eel%cpol, of_unit)
        call print_matrix(.false., 'polarizabilities:', polar, of_unit)
        call print_matrix(.false., 'thole factors:', sys_obj%eel%thole, of_unit)
        call print_int_vec('mm_polar list:', sys_obj%eel%mm_polar, &
                           of_unit)
        call print_int_vec('polar_mm list:', sys_obj%eel%polar_mm, &
                           of_unit)

        ! write the connectivity information for each atom:
1000    format(t3,'connectivity information for the ',i8,'-th atom:')
    
        do i = 1, sys_obj%top%mm_atoms
            write(of_unit, 1000) i
            
            do j=1, 4
                if(j == 4 .and. .not. sys_obj%amoeba) cycle
                
                write(str, "('1-', I1, ' neighors:')") j+1
                call sort_ivec(sys_obj%top%conn(j)%ci(&
                                sys_obj%top%conn(j)%ri(i):&
                                sys_obj%top%conn(j)%ri(i+1)-1), tmp)
                call print_int_vec(trim(str), tmp, of_unit) 
            end do
            
            if(sys_obj%amoeba) then 
                do j=1, 4
                    ilst = 1
                    do igrp = &
                        sys_obj%eel%pg_conn(j)%ri(sys_obj%eel%mmat_polgrp(i)), &
                        sys_obj%eel%pg_conn(j)%ri(sys_obj%eel%mmat_polgrp(i)+1)-1
                    
                        grp = sys_obj%eel%pg_conn(j)%ci(igrp)
                        lst(ilst:ilst+sys_obj%eel%polgrp_mmat%ri(grp+1) - &
                                 sys_obj%eel%polgrp_mmat%ri(grp)-1) = &
                        sys_obj%eel%polgrp_mmat%ci(&
                            sys_obj%eel%polgrp_mmat%ri(grp):&
                            sys_obj%eel%polgrp_mmat%ri(grp+1)-1)

                        ilst = ilst+&
                               sys_obj%eel%polgrp_mmat%ri(grp+1)- &
                               sys_obj%eel%polgrp_mmat%ri(grp)
                    end do
                   
                    write(str, "('1-', I1, ' polarization neighors:')") j
                    if(ilst == 1) ilst = 0
                    call sort_ivec(lst(1:ilst-1), tmp)
                    ! needed to addres the empty array case
                    call print_int_vec(trim(str), tmp, of_unit)
                end do
            end if
        end do
        
        if(present(of_name)) close(of_unit)
        
        if(allocated(tmp)) &
            call mfree('mmpol_ommp_print_summary [tmp]', tmp)
        call mfree('mmpol_ommp_print_summary [polar]', polar)

    end subroutine mmpol_ommp_print_summary
    
    
    subroutine mmpol_save_as_mmp(sys_obj, of_name, r_version)
        !! Save the loaded system in mmpol format. Only the electrostatic part
        !! is saved, everything else is just ignored. Version 2 and 3 of the 
        !! mmp format are supported, 3 is used as default.
        
        use mod_io, only: ommp_message
        use mod_constants, only: mscale_wang_al, &
                                 pscale_wang_al, &
                                 dscale_wang_al, &
                                 uscale_wang_al, &
                                 mscale_wang_dl, &
                                 pscale_wang_dl, &
                                 dscale_wang_dl, &
                                 uscale_wang_dl, &
                                 eps_rp, angstrom2au

        implicit none

        type(ommp_system), target, intent(in) :: sys_obj
        !! System data structure to be saved
        character(len=*), intent(in) :: of_name
        !! Name of the output file
        integer(ip), intent(in), optional :: r_version
        !! Revision version requested for .mmp

        integer(ip) :: of_unit, version, i, jb, je, igrp, inta(120)
        type(ommp_electrostatics_type), pointer :: eel
        type(ommp_topology_type), pointer :: top

        eel => sys_obj%eel
        top => sys_obj%top

        if(.not. sys_obj%mmpol_is_init) then
            call fatal_error('OpenMMPol is not initialized, cannot save .mmp &
                             &file.')
        end if
        
        if(present(r_version)) then
            version = r_version
        else
            version = 3
        end if

        of_unit = 101
        open(unit=of_unit, &
             file=of_name(1:len(trim(of_name))), &
             action='write')
        ! 1.  Integer, revision number used to check if the MMPol.mmp file is 
        !     consistent with the Gaussian version being used
        if(version == 2 .or. version == 3) then
            write(of_unit, '(I0,T2)') version
        else
            call fatal_error("Requested version of .mmp file is not supported &
                             &by openMMPol")
        end if
        
        ! 2.  Integer, MMPol job type:
        !     
        !     0: for a QM/MM calculation (
        !     1: for a Tinker QM/MM calculation
        !     2: for a QM/MM electronic energy transfer (EET) calculation
        write(of_unit, '(I0,T2)') 0
    
        ! 3.  Integer, verbosity flag
        write(of_unit, '(I0,T2)') 3

        ! 4.  Integer, force field type:
        !     
        !     0: Amber-like
        !     1: AMOEBA-like
        if(sys_obj%amoeba) then
            write(of_unit, '(I0,T2)') 1
        else
            write(of_unit, '(I0,T2)') 0
        end if
        
        ! 5.  Integer, force field sub-type. For AMBER:
        !     
        !     0: AMBER Wang-AL
        !     1: AMBER Wang-DL
        !     2: AMBER Thole
        !     0: AMOEBA
        if(sys_obj%amoeba) then
            write(of_unit, '(I0,T2)') 0
        else
            if(all(abs(eel%mscale-mscale_wang_al) < eps_rp) .and. &
               all(abs(eel%pscale-pscale_wang_al) < eps_rp) .and. &
               all(abs(eel%dscale-dscale_wang_al) < eps_rp) .and. &
               all(abs(eel%uscale-uscale_wang_al) < eps_rp)) then
                write(of_unit, '(I0,T2)') 0
            else if(all(abs(eel%mscale-mscale_wang_dl) < eps_rp) .and. &
                    all(abs(eel%pscale-pscale_wang_dl) < eps_rp) .and. &
                    all(abs(eel%dscale-dscale_wang_dl) < eps_rp) .and. &
                    all(abs(eel%uscale-uscale_wang_dl) < eps_rp)) then
                write(of_unit, '(I0,T2)') 1
            else
                call fatal_error("The scaling scheme used is non-standard and &
                                 &therefore cannot be saved in .mmp format")
            end if
        end if
        
        ! 6.  Integer, disabled flag
        write(of_unit, '(I0,T2)') 0

        ! 7.  Real, damping parameter for the dipole-dipole interactions
        !     (Angstrom)
        write(of_unit, '(F16.8)') 0.0

        ! 8.  Integer, solution method for the polarization equations 
        !     (suggested 3):
        !     
        !     0: default, same as 3
        !     1: matrix inversion
        !     2: iterative using Jacobi/DIIS
        !     3: iterative using preconditioned CG solver
        write(of_unit, '(I0,T2)') 0
        
        ! 9.  Integer, how to compute the matrix/vector products for 
        !     MMPol (suggested 3):
        !     
        !     0: default (1 for matrix inversion, 3 for iterative)
        !     1: assemble the matrix incore
        !     2: on-the-fly products with O(N^2) scaling
        !     3: on-the-fly products with the use of FMM if the system is large enough
        !     4: on-the-fly products with the use of FMM in all cases
        write(of_unit, '(I0,T2)') 0
        
        ! 10. Integer, convergence threshold (10^-N) for the iterative 
        !     solvers (suggested 8)
        write(of_unit, '(I0,T2)') 8

        ! 11. Integer, accuracy of the FMM calculations (suggested 0):
        !     
        !     <0: 10^-5 RMS 10^-4 Max Error
        !     0: 10^-6 RMS 10^-5 Max Error
        !     1: 10^-7 RMS 10^-6 Max Error
        !     >2: maximum accuracy
        write(of_unit, '(I0,T2)') 0     
        
        ! 12. Real, FMM box size for MM interactions (Angstrom) (suggested 12.0)
        write(of_unit, '(F16.8)') 12.0

        ! 13. Real, FMM box size for MM/PCM interactions (Angstrom) 
        !     (suggested 6.0)
        if(version == 3) then
            write(of_unit, '(F16.8)') 6.0
        end if

        ! 14. Integer, continuum solvation model:
        !     
        !     0: no continuum solvation model
        !     1: ddCOSMO
        !     2: ddPCM
        write(of_unit, '(I0,T2)') 0     
        
        ! 15. Integer, maximum angular momentum for PCM (suggested 10)
        write(of_unit, '(I0,T2)') 10

        ! 16. Integer, number of Lebedev integration points for PCM 
        !     (suggested 302)
        write(of_unit, '(I0,T2)') 302

        ! 17. Integer, convergence threshold (10^-N) for the PCM 
        !     iterative solvers (suggested 8)
        write(of_unit, '(I0,T2)') 8

        ! 18. Real, dielectric constant of the solvent for PCM (78.3 for water)
        write(of_unit, '(F16.8)') 78.3

        ! 19. Real, optical dielectric constant of the solvent for PCM 
        !     (1.7 for water)
        if(version == 3) then
            write(of_unit, '(F16.8)') 0.00
        end if

        ! 20. Real, relative size of the switching region for PCM 
        !     (suggested 0.1)
        write(of_unit, '(F16.8)') 0.1

        ! 21. Integer, cavity type for PCM (suggested 3):
        !     
        !     0: default, same as 3
        !     1: SAS cavity using UFF radii plus the probe radius
        !     2: Read the cavity from the input file 
        !        (the input must contain at least one sphere per atom)
        !     3: VdW cavity using Bondi's radii scaled by 1.1
        !     4: SAS cavity using Bondi's radii plus the probe radius
        write(of_unit, '(I0,T2)') 0
        
        ! 22. Real, probe radius for the SAS cavity (Angstrom) (1.4 for water)
        write(of_unit, '(F16.8)') 1.4

        ! 23. Integer, number of MM atoms
        write(of_unit, '(I0,T2)') top%mm_atoms

        ! 24. Integer, number of spheres composing the cavity. 
        !     Mandatory if the cavity type is 2, additional spheres 
        !     if the cavity is automatically built.
        write(of_unit, '(I0,T2)') 0
       
        ! 25. Integer array, size (N): atomic numbers
        if(top%atz_initialized) then
            do i=1, top%mm_atoms
                write(of_unit, '(I0,T2)') top%atz(i)
            end do
        else
            do i=1, top%mm_atoms
                write(of_unit, '(I0,T2)') 1
            end do
        end if


        ! 26. Real array, size (N, 3): coordinates
        do i=1, top%mm_atoms
            write(of_unit, '(3F16.8)') top%cmm(:,i) / angstrom2au
        end do

        ! 27. Integer array, size (N): residue numbers
        do i=1, top%mm_atoms
            write(of_unit, '(I0,T2)') 1
        end do

        ! 28. Real array, size (N): charges or fiixed multipoles
        do i=1, top%mm_atoms
            if(sys_obj%amoeba) then
                write(of_unit, '(10F16.8)') eel%q0(:,i) * [1.0, & !q
                                                           1.0, 1.0, 1.0, & !mu
                                                           3.0, 3.0, 3.0, & ! Q
                                                           3.0, 3.0, 3.0]   ! Q

            else
                write(of_unit, '(F16.8)') eel%q(:,i)
            end if
        end do

        ! 29. Real array, size (N): polarizabilities
        do i=1, top%mm_atoms
            if(eel%mm_polar(i) > 0) then
                write(of_unit, '(F16.8)') eel%pol(i) / (angstrom2au**3)
            else
                write(of_unit, '(F16.8)') 0.0
            end if
        end do

        ! 30. Integer array, size (N, 8): connectivity matrix
        do i=1, top%mm_atoms
            jb = top%conn(1)%ri(i)
            je = top%conn(1)%ri(i+1)-1
            inta(1:8) = 0
            inta(1:je-jb+1) = top%conn(1)%ci(jb:je)
            write(of_unit, '(8I8)') inta(1:8)
        end do

        ! 31. Integer array, size (N, 120): polarization group members
        do i=1, top%mm_atoms
            igrp = eel%mmat_polgrp(i)
            jb=eel%polgrp_mmat%ri(igrp)
            je=eel%polgrp_mmat%ri(igrp+1)-1
            inta = 0
            inta(1:je-jb+1) = eel%polgrp_mmat%ci(jb:je)
            write(of_unit, '(120I8)') inta
        end do
        ! 32. Integer array, size (N, 4): rotation conventions and 
        !     reference atoms for the AMOEBA multipoles
        do i=1, top%mm_atoms
            write(of_unit, '(4I8)') eel%mol_frame(i), eel%iz(i), eel%ix(i), &
                                    eel%iy(i)
        end do

    end subroutine mmpol_save_as_mmp

end module mod_mmpol
