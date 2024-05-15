module mod_topology
    use mod_memory, only: ip, rp, lp
    use mod_adjacency_mat, only: yale_sparse

    implicit none 
    private

    type ommp_topology_type
        integer(ip) :: mm_atoms 
        !! number of MM atoms
        real(rp), allocatable :: cmm(:,:)
        !! Coordinates of MM atoms (3:mm_atoms)
        type(yale_sparse), allocatable :: conn(:)
        !! connectivity matrices listing atoms separetad by 1, 2, 3 
        !! (and 4 -- only for AMOEBA) bonds. 1st element is the adjacency 
        !! matrix.
        logical(lp) :: use_frozen = .false.
        !! Flag to use the frozen atom feature, if it is set to true frozen
        !! array is read/used otherwise all atoms are active
        logical(lp), allocatable :: frozen(:)
        !! For each atom, if set to true, it will contribute to the total
        !! energy but its coordinates are locked to the initial value.
        integer(ip), allocatable :: atz(:)
        !! The atomic number for each atom of the system, when it is not 
        !! initialized, it contains only zeros.
        logical(lp) :: atz_initialized = .false.
        !! Initialization flag for atz, when it is filled with actual values
        !! it should be set to true
        real(rp), allocatable :: atmass(:)
        !! The atomic mass for each atom in the system, when it is not initialized
        !! it contains only zeros
        logical(lp) :: atmass_initialized = .false.
        !! Initialization flag for atmass, when it is filled with actual values
        !! it should be set to true
        integer(ip), allocatable :: attype(:)
        !! Atom class for each atom in the system.
        !! It is only used during certain parameters assignaments; when it
        !! is not initialized, it contains only zeros.
        logical(lp) :: attype_initialized = .false.
        !! Initialization flag for attype, when it is filled with actual values
        !! it should be set to true
        integer(ip), allocatable :: atclass(:)
        !! Atom class for each atom in the system.
        !! It is only used during certain parameters assignaments; when it
        !! is not initialized, it contains only zeros.
        logical(lp) :: atclass_initialized = .false.
        !! Initialization flag for atclass, when it is filled with actual values
        !! it should be set to true
    end type ommp_topology_type

    public :: ommp_topology_type
    public :: topology_init, topology_terminate, guess_connectivity
    public :: set_frozen, check_conn_matrix, merge_top, create_new_bond

    contains

        subroutine topology_init(top_obj, mm_atoms)
            use mod_memory, only: mallocate

            implicit none

            type(ommp_topology_type), intent(out) :: top_obj
            integer(ip), intent(in) :: mm_atoms

            top_obj%mm_atoms = mm_atoms
            ! Memory allocation
            call mallocate('topology_init [cmm]', 3_ip, top_obj%mm_atoms, &
                           top_obj%cmm)
            ! Temporary allocation, it should be allocated of the proper
            ! size when all the connectivity matricies are built, now
            ! it should only contain adjacency matrix.
            allocate(top_obj%conn(1))
            top_obj%conn(1)%n = mm_atoms
            
            call mallocate('topology_init [atz]', top_obj%mm_atoms, top_obj%atz)
            top_obj%atz = 0
            top_obj%atz_initialized = .false.
            call mallocate('topology_init [atmass]', top_obj%mm_atoms, top_obj%atmass)
            top_obj%atz = 0.0
            top_obj%atz_initialized = .false.
            call mallocate('topology_init [atclass]', top_obj%mm_atoms, &
                           top_obj%atclass)
            top_obj%atclass = 0
            top_obj%atclass_initialized = .false.
            call mallocate('topology_init [attype]', top_obj%mm_atoms, &
                           top_obj%attype)
            top_obj%attype = 0
            top_obj%attype_initialized = .false.
            
            top_obj%use_frozen = .false.
            if(.not. allocated(top_obj%frozen)) &
                allocate(top_obj%frozen(top_obj%mm_atoms))
            top_obj%frozen = .false.
        end subroutine

        subroutine guess_connectivity(top, exclude_list)
            !! This subroutine guess the connectivity of the system from the 
            !! coordinates and the atomic number of its atoms. It is based
            !! on distances and atomic radii, so it can easily fail on 
            !! distorted geometries. It should be used only when the 
            !!  the bonds of the molecule are not availble in any 
            !! other way; it is often used to assign connectivity to a QM part
            !! that does not have any.
            use mod_constants, only: angstrom2au, &
                                     OMMP_VERBOSE_DEBUG, OMMP_STR_CHAR_MAX
            use mod_io, only: fatal_error, ommp_message
            use mod_memory, only: mallocate, mfree
            use mod_adjacency_mat, only: adj_mat_from_conn
            
            implicit none

            type(ommp_topology_type), intent(inout) :: top
            !! Topology for which connectivity should be guessed.
            integer(ip), optional, intent(in) :: exclude_list(:)
            !! List of atoms that should not be connected when guessing the 
            !! topology, this is used for LA that could have an ambiguous or
            !! not defined position.

            real(rp) :: atomic_radii(118) = 0.0
            !! Covalent atomic radii; 0.0 is used for every unknown element
            !! (that is every element different from H, B, C, N, O, F, Si, P, 
            !! S, Cl, As, Se, Br, Te, I) and means that this is an unexpected 
            !! element and will likely left alone.
            !! The other values are taken from table I of  
            !! J. Chem. Inf. Comput. Sci. 1992, 32, 401-406
            !! "Automatic Assignment of Chemical Connectivity to Organic 
            !!  Molecules in the Cambridge Structural Database"
            !! Values are converted to A.U.
            real(rp), parameter :: tolerance = 0.45 * angstrom2au
            !! Tolerance value for connected atoms, 0.45 Angstrom is used in
            !! agreement with J. Chem. Inf. Comput. Sci. 1992, 32, 401-406
            integer(ip), parameter :: maxbond = 10
            !! Maximum number of atomic of bond for each atom
            integer(ip), allocatable :: i12(:,:)
            !! Temporary bond list that is passed to 
            !! [[mod_adjacency_mat::adj_mat_from_conn]]
            integer(ip), allocatable :: n12(:)
            !! Number of connected atoms already assigned to the i-th atom. 
            integer(ip) :: i, j
            real(rp) :: l, l0
            character(len=OMMP_STR_CHAR_MAX) :: msg
            !! Actual and expected bond length
            
            atomic_radii(1)  = 0.23 * angstrom2au !! H
            atomic_radii(5)  = 0.83 * angstrom2au !! B
            atomic_radii(6)  = 0.68 * angstrom2au !! C
            atomic_radii(7)  = 0.68 * angstrom2au !! N
            atomic_radii(8)  = 0.68 * angstrom2au !! O
            atomic_radii(9)  = 0.64 * angstrom2au !! F
            atomic_radii(14) = 1.20 * angstrom2au !! Si
            atomic_radii(15) = 1.05 * angstrom2au !! P
            atomic_radii(16) = 1.02 * angstrom2au !! S
            atomic_radii(17) = 0.99 * angstrom2au !! Cl
            atomic_radii(33) = 1.21 * angstrom2au !! As
            atomic_radii(34) = 1.22 * angstrom2au !! Se
            atomic_radii(35) = 1.21 * angstrom2au !! Br
            atomic_radii(52) = 1.47 * angstrom2au !! Te
            atomic_radii(53) = 1.40 * angstrom2au !! I
            
            if(.not. top%atz_initialized) then
                call fatal_error('Connectivity cannot be guessed without atomic&
                                 & nubers.')
            end if

            if(allocated(top%conn(1)%ri) .or. allocated(top%conn(1)%ci)) then
                call fatal_error("Connectivity is already present in this &
                                 &topology.")
            end if

            call mallocate('guess_connectivity [i12]', &
                           maxbond, top%mm_atoms, i12)
            call mallocate('guess_connectivity [n12]', &
                           top%mm_atoms, n12)
            i12 = 0
            n12 = 1
            
            do i=1, top%mm_atoms
                if(present(exclude_list)) then
                    if(any(exclude_list == i)) cycle
                end if
                do j=i+1, top%mm_atoms
                    if(present(exclude_list)) then
                        if(any(exclude_list == j)) cycle
                    end if
                    l = norm2(top%cmm(:,i) - top%cmm(:,j))
                    l0 = atomic_radii(top%atz(i)) + atomic_radii(top%atz(j))
                    if(l < l0 + tolerance) then
                        ! There is a bond between i and j
                        i12(n12(i),i) = j
                        i12(n12(j),j) = i
                        n12(i) = n12(i) + 1
                        n12(j) = n12(j) + 1
                        write(msg, '(A, I0, A, I0)') &
                            "Bond found between atoms ", j, &
                            " and ", i
                        call ommp_message(msg, OMMP_VERBOSE_DEBUG)
                    end if
                end do
            end do

            call adj_mat_from_conn(i12, top%conn(1))

            call mfree('guess_connectivity [i12]', i12)
            call mfree('guess_connectivity [n12]', n12)
        end subroutine

        subroutine set_frozen(top_obj, frozen_atoms)
            !! Set the frozen atoms in the current topology, if the
            !! the frozen atoms has already been set, it reinitialize
            !! the whole list, without taking into account the content
            !! of [[top_obj%frozen]]
            use mod_io, only: fatal_error
            implicit none

            type(ommp_topology_type), intent(inout) :: top_obj
            !! Topology object to use
            integer(ip), intent(in) :: frozen_atoms(:)
            !! Indexes of atoms to be frozen

            integer(ip) :: n, i

            top_obj%use_frozen = .true.
            if(.not. allocated(top_obj%frozen)) &
                allocate(top_obj%frozen(top_obj%mm_atoms))

            top_obj%frozen = .false.
            n = size(frozen_atoms)
            do i=1, n
                top_obj%frozen(frozen_atoms(i)) = .true.
            end do
        end subroutine

        subroutine check_conn_matrix(top_obj, n)
            !! Check if adjacency matrix up to nth order is present in
            !! topology object. If it is not present, update the topology
            !! accordingly.
            use mod_adjacency_mat, only: matfree, matcpy, build_conn_upto_n
            use mod_io, only: fatal_error

            implicit none
            type(ommp_topology_type), intent(inout) :: top_obj
            integer(ip), intent(in) :: n

            type(yale_sparse) :: adj
            integer(ip) :: i
            
            if(.not. allocated(top_obj%conn) .or. size(top_obj%conn) < 1) then
                call fatal_error("The topology does not contain a connectivity.")
            end if 

            if(size(top_obj%conn) < n) then
                call matcpy(top_obj%conn(1), adj)
                do i=1, size(top_obj%conn)
                    call matfree(top_obj%conn(i))
                end do

                call build_conn_upto_n(adj, n, top_obj%conn, .false.)

            end if
        end subroutine
        
        subroutine topology_terminate(top_obj)
            use mod_memory, only: mfree
            use mod_adjacency_mat, only: matfree

            implicit none

            type(ommp_topology_type), intent(out) :: top_obj
            integer(ip) :: i

            call mfree('topology_terminate [cmm]', top_obj%cmm)
            call mfree('topology_terminate [atz]', top_obj%atz)
            call mfree('topology_terminate [atmass]', top_obj%atmass)
            call mfree('topology_terminate [atclass]', top_obj%atclass)
            call mfree('topology_terminate [attype]', top_obj%attype)
            
            if(allocated(top_obj%frozen)) &
                deallocate(top_obj%frozen)
            
            if(allocated(top_obj%conn)) then
                do i=1, size(top_obj%conn)
                    call matfree(top_obj%conn(i))
                end do
                deallocate(top_obj%conn)
            endif

        end subroutine

        subroutine create_new_bond(top, i, j)
            !! Create a new bond between atoms i and j in the 
            !! topology.
            use mod_adjacency_mat, only: reallocate_mat
            use mod_io, only: fatal_error
            
            implicit none
            
            type(ommp_topology_type), intent(inout) :: top
            integer(ip), intent(in) :: i, j

            integer(ip) :: l, h, n
            ! Sanity check i and j should be in the topology
            ! and not bonded
            if(i > top%mm_atoms .or. j > top%mm_atoms) then
                call fatal_error("Requested creation of bond outside topology")
            end if
           
            if(any(top%conn(1)%ci(top%conn(1)%ri(i):top%conn(1)%ri(i+1)-1) == j) .or. &
               any(top%conn(1)%ci(top%conn(1)%ri(j):top%conn(1)%ri(j+1)-1) == i)) then
               call fatal_error("Requested bond to create is already present in topology")
            end if

            l = min(i,j)
            h = max(i,j)
            n = size(top%conn(1)%ci) + 2
            call reallocate_mat(top%conn(1), n)
            top%conn(1)%ci(top%conn(1)%ri(l+1)+1:n-1) = &
                top%conn(1)%ci(top%conn(1)%ri(l+1):n-2)
            top%conn(1)%ci(top%conn(1)%ri(l+1)) = h
            top%conn(1)%ri(l+1:) = top%conn(1)%ri(l+1:) + 1
            
            top%conn(1)%ci(top%conn(1)%ri(h+1)+1:n) = &
                top%conn(1)%ci(top%conn(1)%ri(h+1):n-1)
            top%conn(1)%ci(top%conn(1)%ri(h+1)) = l
            top%conn(1)%ri(h+1:) = top%conn(1)%ri(h+1:) + 1

        end subroutine

        subroutine merge_top(top1, top2, top3, map13, map23)
            !! Merge topologies top1 and top2 to create top3
            !! no link between the two topologies are created.
            !! map13 and map23 are arrays mapping the atoms of
            !! top1 to top3 and the ones of top2 to top3 
            !! respectively

            implicit none

            type(ommp_topology_type), intent(in) :: top1, top2
            type(ommp_topology_type), intent(out) :: top3
            integer(ip), allocatable, intent(out) :: map13(:), map23(:)

            integer(ip) :: n, i

            ! Create the new topology
            n = top1%mm_atoms + top2%mm_atoms
            call topology_init(top3, n)
            
            ! Merge coordinates
            top3%cmm(:,1:top1%mm_atoms) = top1%cmm
            top3%cmm(:,top1%mm_atoms+1:n) = top2%cmm

            ! Merge connectivities
            ! TODO check that they are present !
            allocate(top3%conn(1)%ri(n+1))
            allocate(top3%conn(1)%ci(size(top1%conn(1)%ci) + size(top2%conn(1)%ci)))
            allocate(map13(top1%mm_atoms))
            allocate(map23(top2%mm_atoms))

            top3%conn(1)%ri(1) = top1%conn(1)%ri(1)
            do i=1, top1%mm_atoms
                map13(i) = i
                top3%conn(1)%ri(i+1) = top1%conn(1)%ri(i+1)
                top3%conn(1)%ci(top3%conn(1)%ri(i):top3%conn(1)%ri(i+1)-1) = &
                    top1%conn(1)%ci(top1%conn(1)%ri(i):top1%conn(1)%ri(i+1)-1)
            end do
            
            do i=1, top2%mm_atoms
                map23(i) = top1%mm_atoms+i
                top3%conn(1)%ri(top1%mm_atoms+i+1) = &
                    top2%conn(1)%ri(i+1) + top1%conn(1)%ri(top1%mm_atoms+1)-1
                top3%conn(1)%ci(top3%conn(1)%ri(top1%mm_atoms+i):top3%conn(1)%ri(top1%mm_atoms+i+1)-1) = &
                    top2%conn(1)%ci(top2%conn(1)%ri(i):top2%conn(1)%ri(i+1)-1) + top1%mm_atoms
            end do

            ! Merge other properties
            top3%use_frozen = top1%use_frozen .or. top2%use_frozen
            if(top3%use_frozen) top3%frozen = .false.
            if(top1%use_frozen) top3%frozen(1:top1%mm_atoms) = top1%frozen
            if(top2%use_frozen) top3%frozen(top1%mm_atoms+1:n) = top2%frozen
            
            top3%atz_initialized = top1%atz_initialized .and. top2%atz_initialized
            if(top1%atz_initialized) top3%atz(1:top1%mm_atoms) = top1%atz
            if(top2%atz_initialized) top3%atz(top1%mm_atoms+1:n) = top2%atz
            
            top3%atmass_initialized = top1%atmass_initialized .and. top2%atmass_initialized
            if(top1%atmass_initialized) top3%atmass(1:top1%mm_atoms) = top1%atmass(:)
            if(top2%atmass_initialized) top3%atmass(top1%mm_atoms+1:n) = top2%atmass(:)
            
            top3%attype_initialized = top1%attype_initialized .and. top2%attype_initialized
            if(top1%attype_initialized) top3%attype(1:top1%mm_atoms) = top1%attype
            if(top2%attype_initialized) top3%attype(top1%mm_atoms+1:n) = top2%attype

            top3%atclass_initialized = top1%atclass_initialized .and. top2%atclass_initialized
            if(top1%atclass_initialized) top3%atclass(1:top1%mm_atoms) = top1%atclass
            if(top2%atclass_initialized) top3%atclass(top1%mm_atoms+1:n) = top2%atclass
        end subroutine

end module mod_topology
