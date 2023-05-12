#define _MM_ 1
#define _QM_ 2
#define _LA_ 3

module mod_link_atom

    use mod_memory, only: rp, ip, lp, mallocate, mfree
    use mod_topology, only: ommp_topology_type
    use mod_constants, only: angstrom2au

    implicit none
    private

    integer(ip), parameter :: la_allocation_chunk = 20, &
                              default_la_n_eel_remove = 2
    real(rp), parameter :: default_la_dist = 0.91 * angstrom2au 

    type ommp_link_atom_type
        type(ommp_topology_type), pointer :: mmtop
        !! Topology of MM part of the system
        type(ommp_topology_type), pointer :: qmtop
        !! Topology of QM part of the system
        integer(ip), allocatable :: links(:,:)
        !! Indexes of link-involved atoms: 
        !!    links(_MM_,:) contains the index of MM atom in mmtop
        !!    links(_QM_,:) contains the index of QM atom in qmtop
        !!    links(_LA_,:) contains the index of link atom in qmtop
        real(rp), allocatable :: la_distance(:)
        integer(ip) :: nla
        !! Number of link atoms
    end type

    public :: ommp_link_atom_type, init_link_atom, add_link_atom
    public :: link_atom_position
    public :: default_la_dist, default_la_n_eel_remove

    contains
        subroutine init_link_atom(la, qmtop, mmtop)
            implicit none
            
            type(ommp_link_atom_type), intent(inout) :: la
            type(ommp_topology_type), intent(in), target :: qmtop, mmtop

            call mallocate('init_linkatoms [links]', &
                           3, la_allocation_chunk, &
                           la%links)
            call mallocate('init_linkatoms [la_distance]', &
                           la_allocation_chunk, la%la_distance)
            la%nla = 0
            la%qmtop => qmtop
            la%mmtop => mmtop
        end subroutine

        subroutine add_link_atom(la, imm, iqm, ila, la_dist)
            implicit none

            type(ommp_link_atom_type), intent(inout) :: la
            integer(ip), intent(in) :: imm
            integer(ip), intent(in) :: iqm
            integer(ip), intent(in) :: ila
            real(rp), intent(in) :: la_dist
        
            integer(ip), allocatable :: tmp(:,:)
            real(rp), allocatable :: rtmp(:)
            integer(ip) :: nmax
            
            nmax = size(la%links, 2)
            if(la%nla + 1 > nmax) then
                ! Reallocate links to accomodate new link atoms
                call mallocate('create_link_atom [tmp]', 3, nmax, tmp)
                call mallocate('create_link_atom [rtmp]', nmax, rtmp)
                tmp = la%links
                rtmp = la%la_distance
                call mfree('create_link_atom [la%links]', la%links)
                call mallocate('create_link_atom [la%links]', 3, nmax+la_allocation_chunk, la%links)
                call mfree('create_link_atom [la%la_distance]', la%la_distance)
                call mallocate('create_link_atom [la%la_distance]', nmax+la_allocation_chunk, la%la_distance)
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

end module
