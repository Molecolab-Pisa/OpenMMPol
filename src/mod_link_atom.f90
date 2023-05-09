#define _MM_ 1
#define _QM_ 2

module mod_link_atom

    use mod_memory, only: rp, ip, lp, mallocate, mfree
    use mod_topology, only: ommp_topology_type

    implicit none
    private

    integer(ip), parameter :: la_allocation_chunk = 20

    type ommp_link_atom_type
        type(ommp_topology_type), pointer :: mmtop
        !! Topology of MM part of the system
        type(ommp_topology_type), pointer :: qmtop
        !! Topology of QM part of the system
        integer(ip), allocatable :: links(:,:)
        !! Atoms that are linked
        integer(ip) :: nla
        !! Number of link atoms
    end type

    public :: ommp_link_atom_type, init_link_atom, create_link_atom

    contains
        subroutine init_link_atom(la, qmtop, mmtop)
            implicit none
            
            type(ommp_link_atom_type), intent(inout), allocatable :: la
            type(ommp_topology_type), intent(in), target :: qmtop, mmtop

            if(.not. allocated(la)) allocate(la)
            call mallocate('init_linkatoms [links]', &
                           2, la_allocation_chunk, &
                           la%links)
            la%nla = 0
            la%qmtop => qmtop
            la%mmtop => mmtop

        end subroutine

        subroutine create_link_atom(la, imm, iqm)
            !! Create a bond between atoms imm of s and iqm of qmh.
            use mod_io, only: fatal_error, ommp_message
            use mod_constants, only: OMMP_STR_CHAR_MAX, OMMP_VERBOSE_LOW

            implicit none

            type(ommp_link_atom_type), intent(inout), target :: la
            integer(ip), intent(in) :: imm
            integer(ip), intent(in) :: iqm

            character(len=OMMP_STR_CHAR_MAX) :: message
            
            ! 0. Initialization
            !if(.not. associated(s%la%qmh, qmh)) then
            !    call fatal_error("Link Atoms can be created with a single QM &
            !                     &helper object.")
            !end if

            write(message, "(A, I5, A, I5, A)") "Creating link atom between atom ", &
                                                imm, "(MM) and ", iqm, "(QM)."
            call ommp_message(message, OMMP_VERBOSE_LOW) 
            ! 1. Electrostatic
            ! 2. VdW
            ! 3. Bonded

        end subroutine

end module
