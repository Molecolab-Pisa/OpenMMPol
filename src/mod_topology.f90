module mod_topology
    use mod_memory, only: ip, rp
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
    end type ommp_topology_type

    public :: ommp_topology_type
    public :: topology_init, topology_terminate

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
            allocate(top_obj%conn(1)) 
            ! Temporary allocation, it should be allocated of the proper
            ! size when all the connectivity matricies are built, now
            ! it should only contain adjacency matrix.
        end subroutine
        
        subroutine topology_terminate(top_obj)
            use mod_memory, only: mfree
            use mod_adjacency_mat, only: matfree

            implicit none

            type(ommp_topology_type), intent(out) :: top_obj
            integer(ip) :: i

            call mfree('topology_terminate [cmm]', top_obj%cmm)
            
            if(allocated(top_obj%conn)) then
                do i=1, size(top_obj%conn)
                    call matfree(top_obj%conn(i))
                end do
                deallocate(top_obj%conn)
            endif

        end subroutine

end module mod_topology
