module mod_qm_helper
!! This is an utility module, that is not actually used my openMMpol itself,
!! but can be initialized and used by a QM program interfaced with openMMPol
!! to simplify certain steps of the interface using already well tested code.
    
    use mod_memory, only: ip, rp
    use mod_mmpol, only: ommp_system
    use mod_topology, only: ommp_topology_type
    use mod_nonbonded, only: ommp_nonbonded_type

    implicit none
    private

    type ommp_qm_helper
        integer(ip) :: qm_atoms
        !! Number of QM atoms
        type(ommp_topology_type), allocatable :: qm_top
        !! Topology of the QM system (not always needed, it could be 
        !! uninitialized.
        real(rp), allocatable :: cqm(:,:)
        !! Coordinates of QM atoms, shape (3,qm_atoms)
        real(rp), allocatable :: qqm(:)
        !! Charges of QM nuclei, shape(qm_atoms)
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
        !! Electrostatic potential of MMPol atoms at QM nuclei
        logical :: E_m2n_done = .false.
        !! Flag for [[E_m2n]]
        real(rp), allocatable :: E_m2n(:,:)
        !! Electrostatic potential of MMPol atoms at QM nuclei
        type(ommp_nonbonded_type), allocatable :: qm_vdw
        !! Structure to store VdW parameter for QM atoms
    end type ommp_qm_helper

    public :: ommp_qm_helper
    public :: qm_helper_init, qm_helper_terminate
    public :: electrostatic_for_ene, electrostatic_for_grad

    contains
        subroutine qm_helper_init(qm, qmat, cqm, qqm)
            use mod_memory, only: mallocate

            implicit none
            type(ommp_qm_helper), intent(inout) :: qm
            real(rp), intent(in) :: cqm(:,:), qqm(:)
            integer(ip), intent(in) :: qmat

            call mallocate('qm_helper_init [cqm]', 3, qmat, qm%cqm)
            call mallocate('qm_helper_init [qqm]', qmat, qm%qqm)
            
            qm%qm_atoms = qmat
            qm%cqm = cqm
            qm%qqm = qqm

        end subroutine

        subroutine qm_helper_terminate(qm)
            use mod_memory, only: mfree
            use mod_topology, only: topology_terminate
            use mod_nonbonded, only: vdw_terminate

            implicit none
            type(ommp_qm_helper), intent(inout) :: qm

            call mfree('qm_helper_terminate [cqm]', qm%cqm)
            call mfree('qm_helper_terminate [qqm]', qm%qqm)
            if(allocated(qm%qm_vdw)) then
                call vdw_terminate(qm%qm_vdw)
                deallocate(qm%qm_vdw)
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
            
            if(.not. allocated(qm%E_n2p)) then
                call mallocate('electrostatic_for_ene [E_n2p]', &
                               3, system%eel%pol_atoms, qm%E_n2p)
            end if

            qm%E_n2p = 0.0
            do i=1, qm%qm_atoms
                do j=1, system%eel%pol_atoms
                    dr = system%eel%cpol(:,j) - qm%cqm(:,i)
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
            
            if(.not. allocated(qm%V_m2n)) then
                call mallocate('electrostatic_for_ene [V_m2n]', &
                               qm%qm_atoms, qm%V_m2n)
            end if
            qm%V_m2n = 0.0
            call potential_M2E(system%eel, qm%cqm, qm%V_m2n)
            call potential_D2E(system%eel, qm%cqm, qm%V_m2n)
            qm%V_m2n_done = .true.

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
            
            do i=1, qm%qm_atoms
                do j=1, system%top%mm_atoms
                    dr = system%top%cmm(:,j) - qm%cqm(:,i)
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
                               3, qm%qm_atoms, qm%E_m2n)
            end if
            qm%E_m2n = 0.0
            call field_M2E(system%eel, qm%cqm, qm%E_m2n)
            call field_D2E(system%eel, qm%cqm, qm%E_m2n)
            qm%E_m2n_done = .true.
        end subroutine
end module
