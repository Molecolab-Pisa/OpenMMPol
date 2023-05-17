module mod_qm_helper
!! This is an utility module, that is not actually used my openMMpol itself,
!! but can be initialized and used by a QM program interfaced with openMMPol
!! to simplify certain steps of the interface using already well tested code.
    
    use mod_memory, only: ip, rp, lp
    use mod_mmpol, only: ommp_system
    use mod_topology, only: ommp_topology_type
    use mod_nonbonded, only: ommp_nonbonded_type

    implicit none
    private

    type ommp_qm_helper
        type(ommp_topology_type), allocatable :: qm_top
        !! Topology of the QM system
        logical(lp) :: reconnect = .false.
        !! Flag to decide if the topology should be rebuily
        !! at each change of geometry.
        real(rp), allocatable :: qqm(:)
        !! Charges of QM nuclei
        logical :: E_n2p_done = .false.
        !! Flag for [[E_n2p]]
        real(rp), allocatable :: E_n2p(:,:)
        !! Electric field of nuclei at polarizable sites
        logical :: G_n2p_done = .false.
        !! Flag for [[G_n2p]]
        real(rp), allocatable :: G_n2p(:,:)
        !! Electric field gradients of nuclei at polarizable sites
        logical :: E_n2m_done = .false.
        !! Flag for [[E_n2m]]
        real(rp), allocatable :: E_n2m(:,:)
        !! Electric field of nuclei at static sites
        logical :: G_n2m_done = .false.
        !! Flag for [[G_n2m]]
        real(rp), allocatable :: G_n2m(:,:)
        !! Electric field gradients of nuclei at static sites
        logical :: H_n2m_done = .false.
        !! Flag for [[H_n2m]]
        real(rp), allocatable :: H_n2m(:,:)
        !! Electric field Hessian of nuclei at static sites
        logical :: V_m2n_done = .false.
        !! Flag for [[V_m2n]]
        real(rp), allocatable :: V_m2n(:)
        !! Electrostatic potential of MMPol atoms (static) at QM nuclei
        logical :: E_m2n_done = .false.
        !! Flag for [[E_m2n]]
        real(rp), allocatable :: E_m2n(:,:)
        !! Electrostatic potential of MMPol atoms (static) at QM nuclei
        logical :: V_p2n_done = .false.
        !! Flag for [[V_p2n]]
        real(rp), allocatable :: V_p2n(:)
        !! Electrostatic potential of MMPol atoms (polarizable) at QM nuclei
        logical :: E_p2n_done = .false.
        !! Flag for [[E_p2n]]
        real(rp), allocatable :: E_p2n(:,:)
        !! Electrostatic potential of MMPol atoms (polarizable) at QM nuclei
        logical :: use_nonbonded = .false.
        !! Flag for using QM nonbonded terms
        type(ommp_nonbonded_type), allocatable :: qm_vdw
        !! Structure to store VdW parameter for QM atoms
    end type ommp_qm_helper

    public :: ommp_qm_helper
    public :: qm_helper_init, qm_helper_terminate
    public :: qm_helper_init_vdw, qm_helper_init_vdw_prm, &
              qm_helper_vdw_energy, qm_helper_vdw_geomgrad, &
              qm_helper_update_coord, qm_helper_set_attype
    public :: electrostatic_for_ene, electrostatic_for_grad

    contains
        subroutine qm_helper_init(qm, qmat, cqm, qqm, zqm)
            use mod_memory, only: mallocate
            use mod_topology, only: topology_init, guess_connectivity

            implicit none
            type(ommp_qm_helper), intent(inout) :: qm
            !! [[ommp_qm_helper]] object to be initialized
            real(rp), intent(in) :: cqm(:,:)
            !! Coordinates of QM atoms
            real(rp), intent(in) :: qqm(:)
            !! Nuclear charges of QM atoms
            integer(ip), intent(in) :: zqm(:)
            !! Atomic number of QM atoms
            integer(ip), intent(in) :: qmat
            !! Number of QM atoms

            allocate(qm%qm_top)
            call mallocate('qm_helper_init [qqm]', qmat, qm%qqm)
            
            call topology_init(qm%qm_top, qmat)
            qm%qm_top%cmm = cqm
            qm%qm_top%atz = zqm
            qm%qm_top%atz_initialized = .true.
            qm%qqm = qqm

            call guess_connectivity(qm%qm_top)
        end subroutine
        
        subroutine qm_helper_set_attype(qm, attype) 
            implicit none

            type(ommp_qm_helper), intent(inout) :: qm
            integer(ip), intent(in) :: attype(qm%qm_top%mm_atoms)
            
            qm%qm_top%attype = attype
            qm%qm_top%attype_initialized = .true.
        end subroutine

        subroutine qm_helper_update_coord(qm, cnew, reconnect_in)
            use mod_adjacency_mat, only: matfree
            use mod_topology, only: guess_connectivity

            implicit none
            type(ommp_qm_helper), intent(inout) :: qm
            !! [[ommp_qm_helper]] object to be initialized
            real(rp), intent(in) :: cnew(3,qm%qm_top%mm_atoms)
            !! Coordinates of QM atoms
            logical(lp), intent(in), optional :: reconnect_in

            integer(ip) :: i
            logical(lp) :: rc

            if(present(reconnect_in)) then
                rc = reconnect_in
            else
                rc = qm%reconnect
            end if

            qm%qm_top%cmm = cnew
            qm%E_n2p_done = .false.
            qm%G_n2p_done = .false.
            qm%E_n2m_done = .false.
            qm%G_n2m_done = .false.
            qm%H_n2m_done = .false.
            qm%V_m2n_done = .false.
            qm%E_m2n_done = .false.
            if(rc) then
                do i=1, size(qm%qm_top%conn)
                    call matfree(qm%qm_top%conn(i))
                end do
                deallocate(qm%qm_top%conn)
                allocate(qm%qm_top%conn(1))

                call guess_connectivity(qm%qm_top)
            end if
        end subroutine

        subroutine qm_helper_init_vdw(qm, eps, rad, fac, &
                                      vdw_type, radius_rule, radius_size, &
                                      radius_type, epsrule)
            use mod_memory, only: mallocate
            use mod_nonbonded, only: vdw_init
            use mod_io, only: fatal_error

            implicit none

            type(ommp_qm_helper), intent(inout) :: qm
            real(rp), intent(in) :: eps(qm%qm_top%mm_atoms)
            real(rp), intent(in) :: rad(qm%qm_top%mm_atoms)
            real(rp), intent(in) :: fac(qm%qm_top%mm_atoms)
            character(len=*) :: vdw_type, radius_rule, radius_size, &
                                radius_type, epsrule
    
            if(qm%use_nonbonded) then
                call fatal_error("VdW is already initialized!")
            end if
            
            allocate(qm%qm_vdw)

            call vdw_init(qm%qm_vdw, qm%qm_top, vdw_type, radius_rule, &
                          radius_size, radius_type, epsrule)

            qm%qm_vdw%vdw_e = eps
            qm%qm_vdw%vdw_r = rad
            qm%qm_vdw%vdw_f = fac
            qm%use_nonbonded = .true.

        end subroutine
       
        subroutine qm_helper_init_vdw_prm(qm, prmfile)
            !! Assign vdw parameters of the QM part from attype and prm file
            use mod_prm, only: assign_vdw
            use mod_io, only: fatal_error

            implicit none

            type(ommp_qm_helper), intent(inout) :: qm
            character(len=*) :: prmfile
            
            if(qm%use_nonbonded) then
                call fatal_error("VdW is already initialized!")
            end if

            if(.not. qm%qm_top%attype_initialized) then
                call fatal_error("Set the types for QM helper atoms before &
                                 &requesting creation of VdW.")
            end if
            
            allocate(qm%qm_vdw)
            call assign_vdw(qm%qm_vdw, qm%qm_top, prmfile)
            qm%use_nonbonded = .true.
        end subroutine

        subroutine qm_helper_vdw_energy(qm, mm, V)
            use mod_nonbonded, only: vdw_potential_inter, vdw_potential_inter_restricted
            use mod_mmpol, only: ommp_system

            implicit none

            type(ommp_system), intent(in) :: mm
            type(ommp_qm_helper), intent(in) :: qm
            real(rp), intent(inout) :: V

        
            if(mm%use_nonbonded .and. qm%use_nonbonded) then
                call vdw_potential_inter(mm%vdw, qm%qm_vdw, V)
                if(mm%use_linkatoms) then
                    ! Screening due to the presence of link atom
                    call vdw_potential_inter_restricted(mm%vdw, qm%qm_vdw, &
                                                        mm%la%vdw_screening_pairs,&
                                                        mm%la%vdw_screening_f, &
                                                        mm%la%vdw_n_screening, V)
                end if
            end if

        end subroutine
        
        subroutine qm_helper_vdw_geomgrad(qm, mm, qmg, mmg)
            use mod_nonbonded, only: vdw_geomgrad_inter
            use mod_mmpol, only: ommp_system

            implicit none

            type(ommp_system), intent(in) :: mm
            type(ommp_qm_helper), intent(in) :: qm
            real(rp), intent(inout) :: qmg(3,qm%qm_top%mm_atoms), &
                                       mmg(3,mm%top%mm_atoms)
        
            if(mm%use_nonbonded .and. qm%use_nonbonded) then
                call vdw_geomgrad_inter(mm%vdw, qm%qm_vdw, mmg, qmg)
            end if
        end subroutine
    
        subroutine qm_helper_terminate(qm)
            use mod_memory, only: mfree
            use mod_topology, only: topology_terminate
            use mod_nonbonded, only: vdw_terminate

            implicit none
            type(ommp_qm_helper), intent(inout) :: qm

            call mfree('qm_helper_terminate [qqm]', qm%qqm)
            if(allocated(qm%qm_vdw)) then
                call vdw_terminate(qm%qm_vdw)
                deallocate(qm%qm_vdw)
                qm%use_nonbonded = .false.
            end if
            if(allocated(qm%qm_top)) then
                call topology_terminate(qm%qm_top)
                deallocate(qm%qm_top)
            end if
        end subroutine

        subroutine electrostatic_for_ene(system, qm)
            !! Computes the electrostatic quantities (that means nuclear-MM    
            !! interaction terms) needed to perform an energy calculation.
            !! Computes:
            !!    (1) EF of nuclei at polarizable sites
            !!    (2) V of the whole MM system at QM sites
            use mod_memory, only: mallocate
            use mod_electrostatics, only: potential_M2E, potential_D2E, &
                                          q_elec_prop, coulomb_kernel

            implicit none 

            type(ommp_system), intent(in) :: system
            type(ommp_qm_helper), intent(inout) :: qm
            
            real(rp) :: kernel(5), dr(3), tmpV, tmpE(3), tmpEgr(6), &
                        tmpHE(10)
            integer(ip) :: i, j
           
            if(.not. qm%E_n2p_done) then
                if(.not. allocated(qm%E_n2p)) then
                    call mallocate('electrostatic_for_ene [E_n2p]', &
                                3, system%eel%pol_atoms, qm%E_n2p)
                end if

                qm%E_n2p = 0.0
                do i=1, qm%qm_top%mm_atoms
                    do j=1, system%eel%pol_atoms
                        dr = system%eel%cpol(:,j) - qm%qm_top%cmm(:,i)
                        call coulomb_kernel(dr, 1, kernel)

                        tmpE = 0.0
                        call q_elec_prop(qm%qqm(i), dr, kernel, &
                                        .false., tmpV, &
                                        .true., tmpE, &
                                        .false., tmpEgr, & 
                                        .false., tmpHE)
                        
                        qm%E_n2p(:,j) = qm%E_n2p(:,j) + tmpE
                    end do
                end do
                qm%E_n2p_done = .true.
            end if
            
            if(.not. qm%V_m2n_done) then
                if(.not. allocated(qm%V_m2n)) then
                    call mallocate('electrostatic_for_ene [V_m2n]', &
                                qm%qm_top%mm_atoms, qm%V_m2n)
                end if
            
                qm%V_m2n = 0.0
                call potential_M2E(system%eel, qm%qm_top%cmm, qm%V_m2n)
                qm%V_m2n_done = .true.
            end if

            if(.not. qm%V_p2n_done .and. system%eel%ipd_done) then
                if(.not. allocated(qm%V_p2n)) then
                    call mallocate('electrostatic_for_ene [V_p2n]', &
                                qm%qm_top%mm_atoms, qm%V_p2n)
                end if

                qm%V_p2n = 0.0
                call potential_D2E(system%eel, qm%qm_top%cmm, qm%V_p2n)
                qm%V_p2n_done = .true.
            endif

        end subroutine

        subroutine electrostatic_for_grad(system, qm)
            !! Computes the electrostatic quantities (that means nuclear-MM    
            !! interaction terms) needed to perform a gradient calculation.
            !! Computes:
            !!    (1) GEF on nuclei at polarizable sites
            !!    (2) EF, GEF, HEF of nuclei at static sites
            !!    (3) EF of whole MM system at QM sites
            use mod_memory, only: mallocate
            use mod_electrostatics, only: field_M2E, field_D2E, &
                                          q_elec_prop, coulomb_kernel

            implicit none 

            type(ommp_system), intent(in) :: system
            type(ommp_qm_helper), intent(inout) :: qm
            
            real(rp) :: kernel(5), dr(3), tmpV, tmpE(3), tmpEgr(6), &
                        tmpHE(10)
            integer(ip) :: i, j
           
            if(.not. allocated(qm%G_n2p)) then
                call mallocate('electrostatic_for_ene [G_n2p]', &
                               6, system%eel%pol_atoms, qm%G_n2p)
            end if
            if(.not. allocated(qm%E_n2m)) then
                call mallocate('electrostatic_for_ene [E_n2m]', &
                               3, system%top%mm_atoms, qm%E_n2m)
            end if
            if(.not. allocated(qm%G_n2m)) then
                call mallocate('electrostatic_for_ene [G_n2m]', &
                               6, system%top%mm_atoms, qm%G_n2m)
            end if
            if(.not. allocated(qm%H_n2m)) then
                call mallocate('electrostatic_for_ene [H_n2m]', &
                               10, system%top%mm_atoms, qm%H_n2m)
            end if

            qm%E_n2m = 0.0
            qm%G_n2m = 0.0
            qm%H_n2m = 0.0
            
            do i=1, qm%qm_top%mm_atoms
                do j=1, system%top%mm_atoms
                    dr = system%top%cmm(:,j) - qm%qm_top%cmm(:,i)
                    call coulomb_kernel(dr, 3, kernel)

                    tmpE = 0.0
                    tmpEgr = 0.0
                    tmpHE = 0.0
                    call q_elec_prop(qm%qqm(i), dr, kernel, &
                                     .false., tmpV, &
                                     .true., tmpE, &
                                     .true., tmpEgr, & 
                                     .true., tmpHE)
                    
                    qm%E_n2m(:,j) = qm%E_n2m(:,j) + tmpE
                    qm%G_n2m(:,j) = qm%G_n2m(:,j) + tmpEgr
                    qm%H_n2m(:,j) = qm%H_n2m(:,j) + tmpHE
                end do
            end do
            
            qm%E_n2m_done = .true.
            qm%G_n2m_done = .true.
            qm%H_n2m_done = .true.

            do j=1, system%eel%pol_atoms
                qm%G_n2p(:,j) = qm%G_n2m(:,system%eel%polar_mm(j))
            end do
            qm%G_n2p_done = .true.
            
            if(.not. allocated(qm%E_m2n)) then
                call mallocate('electrostatic_for_ene [E_m2n]', &
                               3, qm%qm_top%mm_atoms, qm%E_m2n)
            end if
            qm%E_m2n = 0.0
            call field_M2E(system%eel, qm%qm_top%cmm, qm%E_m2n)
            qm%E_m2n_done = .true.
            
            if(.not. allocated(qm%E_p2n)) then
                call mallocate('electrostatic_for_ene [E_p2n]', &
                               3, qm%qm_top%mm_atoms, qm%E_p2n)
            end if
            qm%E_p2n = 0.0
            call field_D2E(system%eel, qm%qm_top%cmm, qm%E_p2n)
            qm%E_p2n_done = .true.
        end subroutine
end module
