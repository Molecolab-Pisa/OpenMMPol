#define _MM_ 1
#define _QM_ 2
#define _LA_ 3

module mod_link_atom

    use mod_memory, only: rp, ip, lp, mallocate, mfree
    use mod_topology, only: ommp_topology_type
    use mod_constants, only: angstrom2au, default_link_atom_dist, &
                             default_link_atom_n_eel_remove
    use mod_nonbonded, only: vdw_geomgrad_inter
    use mod_bonded, only: ommp_bonded_type
    use mod_io, only: large_file_read

    implicit none
    private
    
    integer(ip), parameter :: link_atom_allocation_chunk = 20

    type ommp_link_atom_type
        type(ommp_topology_type), pointer :: mmtop
        !! Topology of MM part of the system
        type(ommp_topology_type), pointer :: qmtop
        !! Topology of QM part of the system
        type(ommp_topology_type), allocatable :: qmmmtop
        !! Linked QM/MM topology
        integer(ip), allocatable :: mm2full(:), qm2full(:)
        !! Maps from qm and mm systems to full topology indexes
        integer(ip) :: nla
        !! Number of link atoms
        integer(ip), allocatable :: links(:,:)
        !! Indexes of link-involved atoms: 
        !!    links(_MM_,:) contains the index of MM atom in mmtop
        !!    links(_QM_,:) contains the index of QM atom in qmtop
        !!    links(_LA_,:) contains the index of link atom in qmtop
        real(rp), allocatable :: la_distance(:)
        !! Distance of link atom from QM atom
        integer(ip), allocatable :: vdw_screening_pairs(:,:)
        !! Pairs of vdw interactions that should be screened
        real(rp), allocatable :: vdw_screening_f(:)
        !! Screening factors for each vdw pair
        integer(ip) :: vdw_n_screening
        !! Number of screening vdw to be used
        type(ommp_bonded_type), allocatable :: bds
    end type

    public :: ommp_link_atom_type, init_link_atom, add_link_atom
    public :: link_atom_position, init_vdw_for_link_atom, init_bonded_for_link_atom
    public :: init_eel_for_link_atom
    public :: default_link_atom_dist, default_link_atom_n_eel_remove, link_atom_update_merged_topology
    public :: link_atom_bond_geomgrad, link_atom_angle_geomgrad, link_atom_torsion_geomgrad, &
              link_atom_project_grd

    contains
        subroutine init_link_atom(la, qmtop, mmtop)
            use mod_topology, only: merge_top
            use mod_constants, only: OMMP_VERBOSE_LOW
            use mod_io, only: ommp_message

            implicit none
            
            type(ommp_link_atom_type), intent(inout) :: la
            type(ommp_topology_type), intent(in), target :: qmtop, mmtop

            call mallocate('init_linkatoms [links]', &
                           3, link_atom_allocation_chunk, &
                           la%links)
            call mallocate('init_linkatoms [la_distance]', &
                           link_atom_allocation_chunk, la%la_distance)
            call mallocate('init_linkatoms [vdw_screening_f]', &
                           link_atom_allocation_chunk, la%vdw_screening_f)
            call mallocate('init_linkatoms [vdw_screening_pairs]', 2, &
                           link_atom_allocation_chunk, la%vdw_screening_pairs)

            la%nla = 0
            la%vdw_n_screening = 0
            la%qmtop => qmtop
            la%mmtop => mmtop
            
            call ommp_message("Creating merged topology", OMMP_VERBOSE_LOW, 'linkatom')
            allocate(la%qmmmtop)
            call merge_top(la%mmtop, la%qmtop, la%qmmmtop, la%mm2full, la%qm2full)
            allocate(la%bds)
        end subroutine

        subroutine link_atom_update_merged_topology(la)
            !! Update merged topology in linkatom object so that its coordinates
            !! are the same of mmtop and qmtop.
            implicit none

            type(ommp_link_atom_type), intent(inout) :: la

            integer(ip) :: i

            do i=1, la%mmtop%mm_atoms
                la%qmmmtop%cmm(:,la%mm2full(i)) = la%mmtop%cmm(:,i)
            end do

            do i=1, la%qmtop%mm_atoms
                la%qmmmtop%cmm(:,la%qm2full(i)) = la%qmtop%cmm(:,i)
            end do
        end subroutine

        subroutine add_link_atom(la, imm, iqm, ila, la_dist)
            use mod_constants, only: OMMP_STR_CHAR_MAX, OMMP_VERBOSE_LOW
            use mod_io, only: ommp_message
            use mod_topology, only: create_new_bond
            
            implicit none

            type(ommp_link_atom_type), intent(inout) :: la
            integer(ip), intent(in) :: imm
            integer(ip), intent(in) :: iqm
            integer(ip), intent(in) :: ila
            real(rp), intent(in) :: la_dist
        
            integer(ip), allocatable :: tmp(:,:)
            real(rp), allocatable :: rtmp(:)
            integer(ip) :: nmax
            character(len=OMMP_STR_CHAR_MAX) :: message
            
            nmax = size(la%links, 2)
            if(la%nla + 1 > nmax) then
                ! Reallocate links to accomodate new link atoms
                call mallocate('create_link_atom [tmp]', 3, nmax, tmp)
                call mallocate('create_link_atom [rtmp]', nmax, rtmp)
                tmp = la%links
                rtmp = la%la_distance
                call mfree('create_link_atom [la%links]', la%links)
                call mallocate('create_link_atom [la%links]', 3, nmax+link_atom_allocation_chunk, la%links)
                call mfree('create_link_atom [la%la_distance]', la%la_distance)
                call mallocate('create_link_atom [la%la_distance]', nmax+link_atom_allocation_chunk, la%la_distance)
                la%links = tmp(:,0:nmax)
                la%la_distance = rtmp(0:nmax)
                call mfree('create_link_atom [tmp]', tmp)
                call mfree('create_link_atom [rtmp]', rtmp)
            end if
            ! 0.2. Link atom creation and positioning
            la%nla = la%nla + 1
            la%links(_MM_, la%nla) = imm
            la%links(_QM_, la%nla) = iqm
            la%links(_LA_, la%nla) = ila
            la%la_distance(la%nla) = la_dist
           
            if(la%qmmmtop%frozen(la%qm2full(iqm)) .and. la%qmmmtop%frozen(la%mm2full(imm))) then
                ! if QM atom and MM atoms are frozen, also LA is frozen
                la%qmmmtop%frozen(la%qm2full(ila)) = .true.
            end if
            call create_new_bond(la%qmmmtop, la%mm2full(imm), la%qm2full(iqm))

            write(message, "(A, I0, A, I0, A, I0, A)") &
                  "Created link atom MM [", imm, &
                  "] - (LA) [", ila, "] - QM [", iqm, "]"
            
            call ommp_message(message, OMMP_VERBOSE_LOW, 'linkatom')
        end subroutine

        subroutine init_eel_for_link_atom(la, imm, ila, eel, prmfile)
            use mod_memory, only: mallocate, mfree
            use mod_prm, only: assign_mpoles
            use mod_electrostatics, only: ommp_electrostatics_type, &
                                          electrostatics_init, &
                                          remove_null_pol
            use mod_io, only: fatal_error, ommp_message
            use mod_constants, only: eps_rp, OMMP_VERBOSE_LOW, &
                                     OMMP_STR_CHAR_MAX, &
                                     OMMP_VERBOSE_DEBUG

            implicit none

            integer(ip), intent(in) :: imm
            integer(ip), intent(in) :: ila
            type(ommp_link_atom_type), intent(in) :: la
            type(ommp_electrostatics_type), intent(inout) :: eel
            character(len=*), intent(in) :: prmfile

            real(rp) :: removed_charge, qred, old_q
            integer(ip) :: n_eel_remove = default_link_atom_n_eel_remove
            integer(ip) :: ist, i, j, idx, ii
            type(ommp_electrostatics_type) :: tmp_eel
            character(len=OMMP_STR_CHAR_MAX) :: msg
            character(len=OMMP_STR_CHAR_MAX), allocatable :: prm_buf(:)
            integer(ip), allocatable :: attocheck(:)
            
            ! Check if in the complete topology some of the atoms connected
            ! to imm or imm itself do change their parameters. This is a 
            ! problem that affects force-field in which parameters are 
            ! assigned based on connectivity (AMOEBA). In particular we
            ! check only 1,2 and 1,3 neighbours of QM atom that are the ones 
            ! that could be influenced by the new bond. Moreover we only
            ! take care of charges, as multipoles will be removed on those
            ! atoms in a while.

            call mallocate('init_eel_for_link_atom [attocheck]', &
                           eel%top%conn(1)%ri(imm+1) - eel%top%conn(1)%ri(imm)+1, &
                           attocheck)
            attocheck(1) = imm
            ii = 2
            do i=eel%top%conn(1)%ri(imm), eel%top%conn(1)%ri(imm+1)-1
                j = eel%top%conn(1)%ci(i)
                attocheck(ii) = j
                ii = ii + 1
            end do
            write(msg, "(A, I0, A)") "Assigning electrostatic parameter to merged&
                                     & topology. Ignore all the warnings on &
                                     &link atom (", la%qm2full(ila), ")." 
            call ommp_message(msg, OMMP_VERBOSE_LOW, 'linkatom')

            call electrostatics_init(tmp_eel, eel%amoeba, la%qmmmtop%mm_atoms, &
                                     la%qmmmtop)
            
            call large_file_read(prmfile, prm_buf)
            call assign_mpoles(tmp_eel, prm_buf)
            deallocate(prm_buf) 

            old_q = sum(eel%q(1,:))
            do i=1, size(attocheck)
                j = attocheck(i)
                    if(abs(tmp_eel%q(1,la%mm2full(j)) - eel%q(1,j)) > eps_rp) then
                        write(msg, "(A, I0)") "Reassigning electrostatic parameter to&
                                            & atom ", j
                        call ommp_message(msg, OMMP_VERBOSE_LOW, 'linkatom')
                        write(msg, '("  Old parameters: q=", F6.2)') eel%q(1,j)
                        call ommp_message(msg, OMMP_VERBOSE_DEBUG, 'linkatom')
                        write(msg, '("  New parameters: q=", F6.2)') tmp_eel%q(1,la%mm2full(j))
                        call ommp_message(msg, OMMP_VERBOSE_DEBUG, 'linkatom')
                        if(eel%amoeba) then
                            eel%q0(:,j) = tmp_eel%q(:,la%mm2full(j))
                        else
                            eel%q(:,j) = tmp_eel%q(:,la%mm2full(j))
                        end if
                    end if
            end do
            call mfree('init_eel_for_link_atom [attocheck]', attocheck)

            ! Remove dipoles, multipoles, charges and polarizabilities
            !    on all the atoms that have a distance from (QM) atom less or equal to
            !    n_eel_remove. If n_eel_remove is 0, the MM electrostatic is not changed;
            !    if n_eel_remove is 1, only the connected atom is removed and so on.
            !    Removed charges is distributed to all the atoms 1 bond further away.

            if(n_eel_remove > 0) then
                if(n_eel_remove < 2) then
                    call fatal_error("Electrostatic interaction should be removed at least on 1,2 and &
                                    &1,3 neighbour of linked QM atom")
                end if

                removed_charge = 0.0
                
                if(n_eel_remove > size(la%mmtop%conn)) then
                    call fatal_error("Connectivity rebuild is not implemented in link atom")
                end if
                
                ! Remove charges, multipoles and polarizabilities
                if(eel%amoeba) then
                    removed_charge = removed_charge + eel%q0(1,imm)
                    eel%q0(:,imm) = 0.0
                    if(eel%mm_polar(imm) > 0) &
                        eel%pol(eel%mm_polar(imm)) = 0.0
                else
                    removed_charge = removed_charge + eel%q(1,imm)
                    eel%q(:,imm) = 0.0

                    if(eel%mm_polar(imm) > 0) &
                        eel%pol(eel%mm_polar(imm)) = 0.0
                end if
                
                do i=1, n_eel_remove-1
                    do j=eel%top%conn(i)%ri(imm), eel%top%conn(i)%ri(imm+1)-1
                        idx = eel%top%conn(i)%ci(j)
                        if(eel%amoeba) then
                            removed_charge = removed_charge + eel%q0(1,idx)
                            eel%q0(:,idx) = 0.0
                            if(eel%mm_polar(idx) > 0) &
                                eel%pol(eel%mm_polar(idx)) = 0.0
                        else
                            removed_charge = removed_charge + eel%q(1,idx)
                            eel%q(:,idx) = 0.0
                            if(eel%mm_polar(idx) > 0) &
                                eel%pol(eel%mm_polar(idx)) = 0.0
                        end if
                        write(msg, '("Removed charge, multipoles and polarizabilities on atom ", I0)') idx
                        call ommp_message(msg, OMMP_VERBOSE_LOW, 'linkatom')
                    end do
                end do

                ! Redistribute removed charge to preserve neutrality
                qred = removed_charge / (eel%top%conn(n_eel_remove)%ri(imm+1) - eel%top%conn(i)%ri(imm))

                do j=eel%top%conn(n_eel_remove)%ri(imm), &
                    eel%top%conn(n_eel_remove)%ri(imm+1)-1
                    idx = eel%top%conn(n_eel_remove)%ci(j)
                    if(eel%amoeba) then
                        eel%q0(1,idx) = eel%q0(1,idx) + qred
                    else
                        eel%q(1,idx) = eel%q(1,idx) + qred
                    end if
                    write(msg, '("Redistributing charge (", F6.3, " A.U.) on atom ", I0)') qred, idx
                    call ommp_message(msg, OMMP_VERBOSE_LOW, 'linkatom')
                end do
            end if
            
            call remove_null_pol(eel)

            eel%M2M_done = .false.
            eel%M2Mgg_done = .false.
            eel%M2D_done = .false.
            eel%M2Dgg_done = .false.
            eel%ipd_done = .false.
            if(allocated(eel%TMat)) call mfree('update_coordinates [TMat]',eel%TMat)
            if(eel%amoeba) call rotate_multipoles(eel)
            write(msg, '("Charge of the systems passed from ", F6.3, " to ", F6.3, "A.U.")') &
                old_q, sum(eel%q(1,:))
            call ommp_message(msg, OMMP_VERBOSE_LOW, 'linkatom')

            
        end subroutine

        subroutine link_atom_position(la, idx, crd)
            !! Compute the cartesian coordinates of link atom idx at
            !! the current geometry.
            implicit none

            type(ommp_link_atom_type), intent(in) :: la
            integer(ip), intent(in) :: idx
            real(rp), intent(out) :: crd(3)

            real(rp), dimension(3) :: rmm, rqm, dmmqm

            rmm = la%mmtop%cmm(:,la%links(_MM_,idx))
            rqm = la%qmtop%cmm(:,la%links(_QM_,idx))
            dmmqm = rmm-rqm
            dmmqm = dmmqm / norm2(dmmqm)

            crd = rqm + la%la_distance(idx) * dmmqm
        end subroutine

        subroutine check_vdw_pairs(la, n)
            !! Check if n new screening pairs could be allocated
            !! in la structure. If the allocated arrays are too 
            !! small, they are reallocated on-the-fly.
            implicit none

            type(ommp_link_atom_type), intent(inout) :: la
            integer(ip), intent(in) :: n

            integer(ip), allocatable :: itmp(:,:)
            real(rp), allocatable :: rtmp(:)
            integer(ip) :: nnew, nold

            nnew = n + la%vdw_n_screening
            nold = size(la%vdw_screening_f)
            if(nnew > nold) then
                call mallocate('check_vdw_pairs [itmp]', 2, nold, itmp)
                call mallocate('check_vdw_pairs [rtmp]', nold, rtmp)
                itmp = la%vdw_screening_pairs
                rtmp = la%vdw_screening_f
                call mfree('check_vdw_pairs [vdw_screening_pairs]', la%vdw_screening_pairs)
                call mfree('check_vdw_pairs [vdw_screening_f]', la%vdw_screening_f)
                
                call mallocate('check_vdw_pairs [vdw_screening_pairs]', 2, nnew, la%vdw_screening_pairs)
                call mallocate('check_vdw_pairs [vdw_screening_f]', nnew, la%vdw_screening_f)

                la%vdw_screening_pairs = itmp
                la%vdw_screening_f = rtmp

                call mfree('check_vdw_pairs [itmp]', itmp)
                call mfree('check_vdw_pairs [rtmp]', rtmp)
            end if
        end subroutine

        subroutine add_screening_pair(la, iqm, imm, s)
            use mod_io, only: ommp_message
            use mod_constants, only: OMMP_STR_CHAR_MAX, OMMP_VERBOSE_DEBUG
            !! Insert a VdW screening pair in the link atom structure
            implicit none

            type(ommp_link_atom_type), intent(inout) :: la
            integer(ip), intent(in) :: iqm, imm
            real(rp), intent(in) :: s
            character(len=OMMP_STR_CHAR_MAX) :: message

            ! Just for safety, it should be done elsewhere for all the pairs
            ! ones want to insert in the structure.
            call check_vdw_pairs(la, 1)

            la%vdw_n_screening = la%vdw_n_screening + 1
            la%vdw_screening_pairs(_QM_, la%vdw_n_screening) = iqm
            la%vdw_screening_pairs(_MM_, la%vdw_n_screening) = imm
            la%vdw_screening_f(la%vdw_n_screening) = s
            write(message, "(A, I0, A, I0, A, F3.2)") &
                  "Screening VdW interactions between atoms ", imm, " (MM) &
                  &and ", iqm, "(QM) by a factor ", 1.0 + s
            call ommp_message(message, OMMP_VERBOSE_DEBUG, 'linkatom')

        end subroutine

        subroutine init_vdw_for_link_atom(la, iqm, imm, vdw_screening)
            !! Initialize the quantities needed for vdw screening due to the
            !! presence of a link atom between iqm and imm
            use mod_topology, only: check_conn_matrix
            use mod_io, only: fatal_error, ommp_message
            use mod_constants, only: eps_rp, OMMP_STR_CHAR_MAX, OMMP_VERBOSE_DEBUG


            implicit none

            type(ommp_link_atom_type), intent(inout) :: la
            integer(ip), intent(in) :: iqm, imm
            real(rp), intent(in) :: vdw_screening(:)

            integer(ip) :: i, j, ineigh_qm, ineigh_mm, idist, iscr, &
                           qmneigh(4), mmneigh(4)
            real(rp) :: screen
            character(len=OMMP_STR_CHAR_MAX) :: message

            ! Check that the connectivity matrices have been built up to the
            ! required order.
            call check_conn_matrix(la%qmtop, size(vdw_screening)-1)
            call check_conn_matrix(la%mmtop, size(vdw_screening)-1)
            
            ! Count how many interactions should be screened
            qmneigh(1) = 1
            mmneigh(1) = 1
            do i=2, size(vdw_screening)
                qmneigh(i) = la%qmtop%conn(i-1)%ri(iqm+1) - &
                             la%qmtop%conn(i-1)%ri(iqm)
                mmneigh(i) = la%mmtop%conn(i-1)%ri(imm+1) - &
                             la%mmtop%conn(i-1)%ri(imm)
            end do

            iscr = 0
            do idist=1, size(vdw_screening)
                if(abs(vdw_screening(idist) - 1.0) > eps_rp) then
                    do i=1, idist
                        iscr = iscr + qmneigh(i) * mmneigh(idist-i+1)
                    end do
                end if
            end do
            write(message, "(I0, A)") &
                  iscr, " VdW interactions will be screened due to link atom"
            call ommp_message(message, OMMP_VERBOSE_DEBUG, 'linkatom')

            ! Check the allocation of vectors inside link atom structure
            call check_vdw_pairs(la, iscr)

            ! Insert the new screened interactions inside link atom structure
            do idist=1, size(vdw_screening)
                screen = vdw_screening(idist) - 1.0
                if(abs(screen) > eps_rp) then
                    do ineigh_qm=1, idist
                        ineigh_mm = idist - ineigh_qm + 1
                        if(ineigh_qm == 1) then
                            if(ineigh_mm == 1) then
                                call add_screening_pair(la, iqm, imm, screen)
                            else
                                do i=la%mmtop%conn(ineigh_mm-1)%ri(imm), &
                                     la%mmtop%conn(ineigh_mm-1)%ri(imm+1) - 1
                                    call add_screening_pair(la, &
                                                            iqm, &
                                                            la%mmtop%conn(ineigh_mm-1)%ci(i), &
                                                            screen)
                                end do
                            end if
                        else
                            if(ineigh_mm == 1) then
                                do i=la%qmtop%conn(ineigh_qm-1)%ri(iqm), &
                                     la%qmtop%conn(ineigh_qm-1)%ri(iqm+1) -1
                                    call add_screening_pair(la, &
                                                            la%qmtop%conn(ineigh_qm-1)%ci(i), &
                                                            imm, &
                                                            screen)
                                end do
                            else
                                do i=la%mmtop%conn(ineigh_mm-1)%ri(imm), &
                                     la%mmtop%conn(ineigh_mm-1)%ri(imm+1) - 1
                                    do j=la%qmtop%conn(ineigh_qm-1)%ri(iqm), &
                                         la%qmtop%conn(ineigh_qm-1)%ri(iqm+1) -1
                                        call add_screening_pair(la, &
                                                                la%qmtop%conn(ineigh_qm-1)%ci(j), &
                                                                la%mmtop%conn(ineigh_mm-1)%ci(i), &
                                                                screen)
                                    end do
                                end do
                            end if
                        end if
                    end do
                end if
            end do

        end subroutine

        subroutine init_bonded_for_link_atom(la, prmfile)
            !! Insert in the bonded parameter required for the link atom between iqm and imm
            use mod_prm, only: assign_bond, assign_angle, assign_torsion
            use mod_bonded, only: bonded_terminate, &
                                  bond_init, angle_init, torsion_init, &
                                  bond_terminate, angle_terminate, torsion_terminate
            use mod_bonded, only: OMMP_ANG_SIMPLE, OMMP_ANG_H0, OMMP_ANG_H1, &
                                  OMMP_ANG_H2, OMMP_ANG_INPLANE, &
                                  OMMP_ANG_INPLANE_H0, OMMP_ANG_INPLANE_H1
            use mod_io, only: ommp_message, fatal_error
            use mod_constants, only: OMMP_STR_CHAR_MAX, OMMP_VERBOSE_LOW
            use mod_topology, only: check_conn_matrix

            implicit none

            type(ommp_link_atom_type), intent(inout), target :: la
            character(len=*), intent(in) :: prmfile
            
            type(ommp_bonded_type) :: tmp_bnd
            integer(ip), parameter :: maxt = 1024
            integer(ip) :: ist, i, j, nqm, nterms, iterms(maxt)
            character(len=OMMP_STR_CHAR_MAX) :: message
            character(len=OMMP_STR_CHAR_MAX), allocatable :: prm_buf(:)

            call check_conn_matrix(la%qmmmtop, 4)
            tmp_bnd%top => la%qmmmtop
            la%bds%top => la%qmmmtop

            call large_file_read(prmfile, prm_buf)
            ! Bonded terms
            call assign_bond(tmp_bnd, prm_buf, la%qm2full, 2)
            if(tmp_bnd%use_bond) then
                nterms = 0
                do i=1, tmp_bnd%nbond
                    nqm = 0
                    if(any(la%qm2full == tmp_bnd%bondat(1,i))) nqm = nqm+1
                    if(any(la%qm2full == tmp_bnd%bondat(2,i))) nqm = nqm+1
                    if(nqm == 1) then
                        nterms = nterms + 1
                        if(nterms > maxt) then
                            call fatal_error("Maximum number of terms to be added due to link atom &
                                            &exceeded. This is probably a bug.")
                        end if
                        iterms(nterms) = i
                    end if
                end do
                
                if(la%bds%use_bond) then
                    ! Bond terms are already initializad by a previous
                    ! link atom
                    call bond_terminate(la%bds) 
                end if

                call bond_init(la%bds, nterms)

                do i=1, nterms
                    la%bds%bondat(:,i) = tmp_bnd%bondat(:,iterms(i))
                    la%bds%kbond(i) = tmp_bnd%kbond(iterms(i))
                    la%bds%l0bond(i) = tmp_bnd%l0bond(iterms(i))
                end do
                la%bds%bond_cubic = tmp_bnd%bond_cubic
                la%bds%bond_quartic = tmp_bnd%bond_quartic
                
                write(message, "(I0, A)") nterms, " bond terms added due to link atoms."
                call ommp_message(message, OMMP_VERBOSE_LOW, "linkatom")
            end if
            
            call assign_angle(tmp_bnd, prm_buf, la%qm2full, 2)
            if(tmp_bnd%use_angle) then
                nterms = 0
                do i=1, tmp_bnd%nangle
                    if(tmp_bnd%anglety(i) == OMMP_ANG_SIMPLE .or. &
                       tmp_bnd%anglety(i) == OMMP_ANG_H0 .or. &
                       tmp_bnd%anglety(i) == OMMP_ANG_H1 .or. &
                       tmp_bnd%anglety(i) == OMMP_ANG_H2) then
                        nqm = 0
                        if(any(la%qm2full == tmp_bnd%angleat(1,i))) nqm = nqm+1
                        if(any(la%qm2full == tmp_bnd%angleat(2,i))) nqm = nqm+1
                        if(any(la%qm2full == tmp_bnd%angleat(3,i))) nqm = nqm+1

                        if(nqm == 1) then
                            nterms = nterms + 1
                            if(nterms > maxt) then
                                call fatal_error("Maximum number of terms to be added due to link atom &
                                                &exceeded. This is probably a bug.")
                            end if
                            iterms(nterms) = i
                        end if
                    else if(tmp_bnd%anglety(i) == OMMP_ANG_INPLANE .or. &
                            tmp_bnd%anglety(i) == OMMP_ANG_INPLANE_H0 .or. &
                            tmp_bnd%anglety(i) == OMMP_ANG_INPLANE_H1) then
                        nqm = 0
                        if(any(la%qm2full == tmp_bnd%angleat(1,i))) nqm = nqm+1
                        if(any(la%qm2full == tmp_bnd%angleat(2,i))) nqm = nqm+1
                        if(any(la%qm2full == tmp_bnd%angleat(3,i))) nqm = nqm+1
                        if(nqm == 1 .and. .not. any(la%qm2full == tmp_bnd%angauxat(i))) then
                            nterms = nterms + 1
                            if(nterms > maxt) then
                                call fatal_error("Maximum number of terms to be added due to link atom &
                                                &exceeded. This is probably a bug.")
                            end if
                            iterms(nterms) = i
                        end if
                    end if
                end do
                
                if(la%bds%use_angle) then
                    call angle_terminate(la%bds)
                end if
                call angle_init(la%bds, nterms)

                do i=1, nterms
                    la%bds%angleat(:,i) = tmp_bnd%angleat(:,iterms(i))
                    la%bds%anglety(i) = tmp_bnd%anglety(iterms(i))
                    la%bds%angauxat(i) = tmp_bnd%angauxat(iterms(i))
                    la%bds%kangle(i) = tmp_bnd%kangle(iterms(i))
                    la%bds%eqangle(i) = tmp_bnd%eqangle(iterms(i))
                end do
                la%bds%angle_cubic = tmp_bnd%angle_cubic
                la%bds%angle_quartic = tmp_bnd%angle_quartic
                la%bds%angle_pentic = tmp_bnd%angle_pentic
                la%bds%angle_sextic = tmp_bnd%angle_sextic
                
                write(message, "(I0, A)") nterms, " angle terms added due to link atoms."
                call ommp_message(message, OMMP_VERBOSE_LOW, "linkatom")
            end if
            
            call assign_torsion(tmp_bnd, prm_buf)
            if(tmp_bnd%use_torsion) then
                nterms = 0
                do i=1, tmp_bnd%ntorsion
                    nqm = 0
                    if(any(la%qm2full == tmp_bnd%torsionat(1,i))) nqm = nqm+1
                    if(any(la%qm2full == tmp_bnd%torsionat(2,i))) nqm = nqm+1
                    if(any(la%qm2full == tmp_bnd%torsionat(3,i))) nqm = nqm+1
                    if(any(la%qm2full == tmp_bnd%torsionat(4,i))) nqm = nqm+1
                    !! Link atoms should never appear in bonded interactions!
                    do j=1, 4
                        if(any(la%qm2full(la%links(_LA_,1:la%nla)) == tmp_bnd%torsionat(j,i))) nqm = -1
                    end do

                    if(nqm == 1 .or. nqm == 2) then
                        nterms = nterms + 1
                        if(nterms > maxt) then
                            call fatal_error("Maximum number of terms to be added due to link atom &
                                            &exceeded. This is probably a bug.")
                        end if
                        iterms(nterms) = i
                    end if
                end do

                if(la%bds%use_torsion) then
                    call torsion_terminate(la%bds)
                end if
                call torsion_init(la%bds, nterms)
                
                do i=1, nterms
                    la%bds%torsionat(:,i) = tmp_bnd%torsionat(:,iterms(i))
                    la%bds%torsn(:,i) = tmp_bnd%torsn(:,iterms(i))
                    la%bds%torsamp(:,i) = tmp_bnd%torsamp(:,iterms(i))
                    la%bds%torsphase(:,i) = tmp_bnd%torsphase(:,iterms(i))
                end do
                
                write(message, "(I0, A)") nterms, " torsion terms added due to link atoms."
                call ommp_message(message, OMMP_VERBOSE_LOW, "linkatom")
            end if
            
            deallocate(prm_buf)
            call bonded_terminate(tmp_bnd)
        end subroutine

        subroutine link_atom_bond_geomgrad(la, qmg, mmg, doqm, domm) 
            use mod_bonded, only: bond_geomgrad

            implicit none

            type(ommp_link_atom_type), intent(in) :: la
            real(rp), intent(inout) :: qmg(:,:), &
                                       mmg(:,:)
            logical, intent(in) :: doqm, domm

            real(rp), allocatable :: grd(:,:)
            integer(ip) :: i

            call mallocate('link_atom_bond_geomgrad [grd]', &
                           3, la%qmmmtop%mm_atoms, grd)
            grd = 0.0
            call bond_geomgrad(la%bds, grd)
            
            if(doqm) then
                do i=1, la%qmtop%mm_atoms
                    qmg(:,i) = qmg(:,i) + grd(:,la%qm2full(i))
                end do
            end if

            if(domm) then
                do i=1, la%mmtop%mm_atoms
                    mmg(:,i) = mmg(:,i) + grd(:,la%mm2full(i))
                end do
            end if

            call mfree('link_atom_bond_geomgrad [grd]', grd)
        end subroutine
        
        subroutine link_atom_angle_geomgrad(la, qmg, mmg, doqm, domm) 
            use mod_bonded, only: angle_geomgrad

            implicit none

            type(ommp_link_atom_type), intent(in) :: la
            real(rp), intent(inout) :: qmg(:,:), &
                                       mmg(:,:)
            logical, intent(in) :: doqm, domm

            real(rp), allocatable :: grd(:,:)
            integer(ip) :: i

            call mallocate('link_atom_bond_geomgrad [grd]', &
                           3, la%qmmmtop%mm_atoms, grd)
            grd = 0.0
            call angle_geomgrad(la%bds, grd)
            
            if(doqm) then
                do i=1, la%qmtop%mm_atoms
                    qmg(:,i) = qmg(:,i) + grd(:,la%qm2full(i))
                end do
            end if

            if(domm) then
                do i=1, la%mmtop%mm_atoms
                    mmg(:,i) = mmg(:,i) + grd(:,la%mm2full(i))
                end do
            end if

            call mfree('link_atom_bond_geomgrad [grd]', grd)
        end subroutine

        subroutine link_atom_torsion_geomgrad(la, qmg, mmg, doqm, domm)
            use mod_bonded, only: torsion_geomgrad

            implicit none

            type(ommp_link_atom_type), intent(in) :: la
            real(rp), intent(inout) :: qmg(:,:), &
                                       mmg(:,:)
            logical, intent(in) :: doqm, domm

            real(rp), allocatable :: grd(:,:)
            integer(ip) :: i

            call mallocate('link_atom_bond_geomgrad [grd]', &
                           3, la%qmmmtop%mm_atoms, grd)
            grd = 0.0
            call torsion_geomgrad(la%bds, grd)
            
            if(doqm) then
                do i=1, la%qmtop%mm_atoms
                    qmg(:,i) = qmg(:,i) + grd(:,la%qm2full(i))
                end do
            end if

            if(domm) then
                do i=1, la%mmtop%mm_atoms
                    mmg(:,i) = mmg(:,i) + grd(:,la%mm2full(i))
                end do
            end if

            call mfree('link_atom_bond_geomgrad [grd]', grd)
        end subroutine
        
        subroutine link_atom_project_grd(la, laforces, qmg, mmg)
            use mod_utils, only: versor_der

            implicit none

            type(ommp_link_atom_type), intent(in) :: la
            real(rp), intent(in) :: laforces(3, la%nla)
            real(rp), intent(inout) :: qmg(3,la%qmtop%mm_atoms), &
                                       mmg(3,la%mmtop%mm_atoms)

            integer(ip) :: i, iqm, imm, ila
            real(rp) :: delta(3), rmm(3), rqm(3), dedqm(3,3), dedmm(3,3)

            do i=1, la%nla
                iqm = la%links(_QM_,i)
                imm = la%links(_MM_,i)
                ila = la%links(_LA_,i)

                rmm = la%mmtop%cmm(:,imm)
                rqm = la%qmtop%cmm(:,iqm)
                delta = rmm-rqm

                dedmm = la%la_distance(i) * versor_der(delta)
               
                dedqm = 0.0
                dedqm(1,1) = 1.0
                dedqm(2,2) = 1.0
                dedqm(3,3) = 1.0
                dedqm = dedqm - dedmm
                
                if(.not. la%qmmmtop%frozen(la%qm2full(iqm))) &
                    qmg(:, iqm) = qmg(:, iqm) + matmul(dedqm, laforces(:,i))
                if(.not. la%qmmmtop%frozen(la%mm2full(imm))) &
                    mmg(:, imm) = mmg(:, imm) + matmul(dedmm, laforces(:,i))
                qmg(:, ila) = -laforces(:,i)
            end do
        end subroutine
end module
