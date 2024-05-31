module mod_tree
    use mod_constants, only: ip, rp, lp
    use mod_adjacency_mat, only: yale_sparse
    use mod_fmm_utils, only: fmm_error
    use mod_profiling, only: time_push, time_pull

    implicit none
    private

    type fmm_tree_type
        !! Data structore to represent the tree for storing
        !! sources and targets of the fmm problem
        integer(ip) :: tree_degree = 0
        !! Maximum number of children for each node
        integer(ip) :: n_nodes = 0
        !! Number of nodes in the tree
        integer(ip) :: breadth = 0
        !! Number of levels in the tree
        integer(ip), allocatable :: children(:,:)
        !! For each node in the tree, its children (0 is the null element)
        integer(ip), allocatable :: parent(:)
        !! For each node in the tree, its parent (0 is the null element)
        integer(ip), allocatable :: node_level(:)
        !! For each node the level of the tree in which it is located (1-based)
        logical(lp), allocatable :: is_leaf(:)
        !! For each node true if is a leaf and false otherwise
        type(yale_sparse) :: level_list
        !! For each level, the list of the nodes contained
        integer(ip) :: n_particles = 0
        !! Number of particles stored in the tree
        real(rp), pointer :: particles_coords(:, :)
        !! Pointer to particle coordinates
        type(yale_sparse) :: particle_list
        !! List of particles contained in each node, for all non-leaf nodes
        !! this list should be empty
        integer(ip), allocatable :: particle_to_node(:)
        ! For each particle the node in which is placed (basically [[particle_list]] inverted)
        real(rp), allocatable :: node_centroid(:,:)
        !! For each node its centroid (in cartesian coordinates)
        real(rp), allocatable :: node_dimension(:)
        !! For each node, its dimension (or radius)
        type(yale_sparse) :: near_nl
        !! List of nodes pair eligible for near-field
        type(yale_sparse) :: far_nl
        !! List of nodes pair eligible for far-field
    end type

    public :: fmm_tree_type, free_tree, print_tree, allocate_tree, tree_populate_farnear_lists, &
              populate_level_list, populate_leaf_list
    public :: tree_populate_farnear_lists_safe

    contains

    subroutine allocate_tree(t)
        !! Given a tree with all the dimension set, it allocates all the required 
        !! arrays and data stractures for populating the tree
        use mod_adjacency_mat, only: allocate_yale_sparse

        implicit none
        
        type(fmm_tree_type), intent(inout) :: t

        if(t%tree_degree == 0 .or. t%n_nodes == 0 .or. &
           t%breadth == 0 .or. t%n_particles == 0) then
            call fmm_error()
        end if
        
        allocate(t%children(t%tree_degree, t%n_nodes))
        allocate(t%parent(t%n_nodes))
        allocate(t%node_level(t%n_nodes))
        allocate(t%node_centroid(3,t%n_nodes))
        allocate(t%node_dimension(t%n_nodes))
        allocate(t%is_leaf(t%n_nodes))
        allocate(t%particle_to_node(t%n_particles))
        call allocate_yale_sparse(t%level_list, t%breadth, t%n_nodes)
        call allocate_yale_sparse(t%particle_list, t%n_nodes, t%n_particles)
        ! Guess on required size for total number of neigbour, if needed
        ! it will be shrinked or enlarged in [[tree_populate_farnear_lists]]
        call allocate_yale_sparse(t%near_nl, t%n_nodes, t%n_nodes)
        call allocate_yale_sparse(t%far_nl, t%n_nodes, t%n_nodes)
    end subroutine
    
    subroutine free_tree(t)
        !! Frees all allocatable quantities contained inside a tree

        use mod_adjacency_mat, only: free_yale_sparse

        implicit none
        
        type(fmm_tree_type), intent(inout) :: t

        t%tree_degree = 0
        t%n_nodes = 0
        t%breadth = 0 
        t%n_particles = 0
        if(allocated(t%children)) deallocate(t%children)
        if(allocated(t%parent)) deallocate(t%parent)
        if(allocated(t%node_level)) deallocate(t%node_level)
        if(allocated(t%node_centroid)) deallocate(t%node_centroid)
        if(allocated(t%node_dimension)) deallocate(t%node_dimension)
        if(allocated(t%is_leaf)) deallocate(t%is_leaf)
        if(allocated(t%particle_to_node)) deallocate(t%particle_to_node)
        call free_yale_sparse(t%level_list)
        call free_yale_sparse(t%particle_list)
        call free_yale_sparse(t%near_nl)
        call free_yale_sparse(t%far_nl)
    end subroutine

    subroutine print_tree(t)
        implicit none
        
        type(fmm_tree_type), intent(in) :: t
        !! Tree data structure
        integer(ip) :: i, j

        do i=1, t%n_nodes
            write(*, *) "NODE", i
            write(*, *) "** PARENT", t%parent(i)
            write(*, *) "** CHILDREN", t%children(:,i)
            write(*, *) "** LEVEL", t%node_level(i)
            write(*, *) "** CENTROID", t%node_centroid(:,i)
            write(*, *) "** PARTICLES"
            do j=t%particle_list%ri(i), t%particle_list%ri(i+1)-1
                write(*, *) "**** PT", t%particle_list%ci(j)
            end do
        end do

        do i=1, t%breadth
            write(*, *) "LEVEL", i
            do j=t%level_list%ri(i), t%level_list%ri(i+1)-1
                write(*, *) "** ", t%level_list%ci(j)
            end do 
        end do
    end subroutine

    subroutine tree_populate_farnear_lists_safe(t, min_dist_thr)
        !! Just for testing, creates far and near list using a double loop
        !! algorithm, it is basically just the application of the following
        !! definition:
        !!     1. Two nodes are near IF they are both leaves and if the 
        !!        distance is below [min_dist_thr]
        !!     2. Two nodes are far IF none of their discendent are near
        !!     3. Descendent of two far nodes are not present in any list

        use mod_adjacency_mat, only: compress_list

        implicit none 

        type(fmm_tree_type), intent(inout) :: t
        !! Tree data structure to populate
        real(rp), intent(in), optional :: min_dist_thr
        !! Minimum threshold for two nodes to be near, every nodes within 
        !! this threshold are guaranteed to be near
        integer(ip) :: i, j
        real(rp) :: d

        integer(ip), allocatable :: uncompressed_far(:,:), uncompressed_n_far(:)
        integer(ip), allocatable :: uncompressed_near(:,:), uncompressed_n_near(:)
        integer(ip) :: chunk_size, size_near

        call time_push
        ! Just a guess, if needed the vectors are enlarged, a better guess is a good TODO
        chunk_size = min(2000, t%n_nodes)
        
        size_near = chunk_size
        allocate(uncompressed_far(size_near, t%n_nodes))
        allocate(uncompressed_n_far(t%n_nodes))
        uncompressed_n_far = 0
        allocate(uncompressed_near(size_near, t%n_nodes))
        allocate(uncompressed_n_near(t%n_nodes))
        uncompressed_n_near = 0

        do i=1, t%n_nodes
            do j=i, t%n_nodes
                if(t%is_leaf(i) .and. t%is_leaf(j)) then
                    ! Both of them are leavs
                    d = norm2(t%node_centroid(:,i) - t%node_centroid(:,j))
                    
                    if(d - t%node_dimension(i) - t%node_dimension(j) < min_dist_thr) then
                        uncompressed_n_near(i) = uncompressed_n_near(i) + 1
                        uncompressed_near(uncompressed_n_near(i),i) = j
                       
                        if(i /= j) then
                            uncompressed_n_near(j) = uncompressed_n_near(j) + 1
                            uncompressed_near(uncompressed_n_near(j),j) = i
                        end if
                    else
                        uncompressed_n_far(i) = uncompressed_n_far(i) + 1
                        uncompressed_far(uncompressed_n_far(i),i) = j
                        
                        if(i /= j) then
                            uncompressed_n_far(j) = uncompressed_n_far(j) + 1
                            uncompressed_far(uncompressed_n_far(j),j) = i
                        end if
                    end if
                end if
            end do
        end do
        
        call time_push
        call aggregative_pass(t, uncompressed_n_far, uncompressed_far)
        call time_pull("Aggregative pass")
        
        call time_push
        !! Now compress in yale sparse format and delete the uncompressed lists
        call compress_list(t%n_nodes, uncompressed_far, uncompressed_n_far, t%far_nl)
        call compress_list(t%n_nodes, uncompressed_near, uncompressed_n_near, t%near_nl)
        call time_pull("List compression")

        ! In the end remove the pair list that is still allocated
        if(allocated(uncompressed_near)) deallocate(uncompressed_near)
        if(allocated(uncompressed_n_near)) deallocate(uncompressed_n_near)
        if(allocated(uncompressed_far)) deallocate(uncompressed_far)
        if(allocated(uncompressed_n_far)) deallocate(uncompressed_n_far)
        call time_pull("Far near list")

    end subroutine
    
    subroutine aggregative_pass(t, uncompressed_n_far, uncompressed_far)
        !! Aggregates far nodes to upper levels of the tree. Operation are 
        !! performed on uncompressed list

        implicit none

        type(fmm_tree_type), intent(in) :: t
        integer(ip), intent(inout) :: uncompressed_n_far(t%n_nodes)
        integer(ip), intent(inout) :: uncompressed_far(:,:)

        integer(ip) :: l, j, i, k, kk, jj, parent, &
                       max_pass, pass, nb, na
        logical :: all_children_far
        integer(ip), parameter :: max_total_pass = 100

        max_pass = max(max_total_pass, t%breadth)

        do pass=1, max_pass
            nb = sum(uncompressed_n_far)
            do l=t%breadth, 2, -1
                !$omp parallel do default(shared) private(i,jj, parent, all_children_far, k, kk)
                do i=1, t%n_nodes
                    do jj=1, uncompressed_n_far(i)
                        j = uncompressed_far(jj,i)
                        if(t%node_level(j) /= l) cycle

                        parent = t%parent(j)
                       
                        all_children_far = .true.
                        do k=1, t%tree_degree
                            if(t%children(k,parent) == 0) cycle
                            all_children_far = all_children_far .and. &
                                               any(uncompressed_far(:uncompressed_n_far(i),i) == t%children(k,parent))
                        end do

                        if(all_children_far) then
                            ! Replace every children with parent in i-node list
                            kk = 1
                            do k=1, uncompressed_n_far(i)
                                if(all(t%children(:,parent) /= uncompressed_far(k,i))) then
                                    uncompressed_far(kk,i) = uncompressed_far(k,i)
                                    kk = kk + 1
                                end if
                            end do
                            uncompressed_n_far(i) = kk - 1
                            if(all(uncompressed_far(1:uncompressed_n_far(i),i) /= parent)) then
                                uncompressed_n_far(i) = uncompressed_n_far(i) + 1
                                uncompressed_far(uncompressed_n_far(i),i) = parent
                            end if
                        end if
                    end do
                end do
            end do
            na = sum(uncompressed_n_far)
            if(na >= nb) exit
        end do
    end subroutine

    subroutine tree_populate_farnear_lists(t, min_dist_thr)
        use mod_adjacency_mat, only: compress_list

        implicit none 

        type(fmm_tree_type), intent(inout) :: t
        !! Tree data structure to populate
        real(rp), intent(in), optional :: min_dist_thr
        !! Minimum threshold for two nodes to be near, every nodes within 
        !! this threshold are guaranteed to be near
        integer(ip) :: i, j, k, ii, jj, ilev
        real(rp) :: d

        integer(ip), allocatable :: uncompressed_far(:,:), uncompressed_n_far(:)
        integer(ip), allocatable :: uncompressed_near(:,:), uncompressed_n_near(:)
        integer(ip), allocatable :: uncompressed_nnear(:,:), uncompressed_n_nnear(:)
        integer(ip), allocatable :: i_checklist(:), j_checklist(:)
        integer(ip) :: chunk_size, size_near, kk, m, n, m_node, n_node

        call time_push
        ! Just a guess, if needed the vectors are enlarged, a better guess is a good TODO
        chunk_size = min(10000, t%n_nodes)
        
        size_near = chunk_size
        allocate(uncompressed_far(size_near, t%n_nodes))
        allocate(uncompressed_n_far(t%n_nodes))
        uncompressed_n_far(:) = 0
        allocate(uncompressed_near(size_near, t%n_nodes))
        allocate(uncompressed_n_near(t%n_nodes))
        uncompressed_n_near(:) = 0
        
        ! Used to keep track of non-leaf nodes that are near
        allocate(uncompressed_nnear(size_near, t%n_nodes))
        allocate(uncompressed_n_nnear(t%n_nodes))
        uncompressed_n_nnear(:) = 0

        allocate(i_checklist(t%tree_degree+1))
        allocate(j_checklist(t%tree_degree+1))

        ! Root node is near to itself
        uncompressed_n_nnear(1) = 1
        uncompressed_nnear(1,1) = 1
        
        call time_push
        do ilev=1, t%breadth-1
            ! For each level starting from the first one
            do ii=t%level_list%ri(ilev), t%level_list%ri(ilev+1)-1
                ! For each node in this level
                i = t%level_list%ci(ii)

                i_checklist(:) = 0
                if(t%is_leaf(i)) then
                    ! This node is a leaf
                    i_checklist(1) = i
                else
                    kk = 1
                    do k=1, t%tree_degree
                        if(t%children(k,i) /= 0) then
                            i_checklist(kk) = t%children(k,i)
                            kk = kk + 1
                        end if
                    end do
                end if

                do jj=1, uncompressed_n_nnear(i)
                    j = uncompressed_nnear(jj,i)
                
                    j_checklist(:) = 0
                    if(t%is_leaf(j)) then
                        ! This node is a leaf
                        j_checklist(1) = j
                    else
                        kk = 1
                        do k=1, t%tree_degree
                            if(t%children(k,j) /= 0) then
                                j_checklist(kk) = t%children(k,j)
                                kk = kk + 1
                            end if
                        end do
                    end if

                    do m=1, t%tree_degree + 1
                        m_node = i_checklist(m)
                        if(m_node == 0) cycle
                        do n=1, t%tree_degree + 1
                            n_node = j_checklist(n)
                            if(n_node == 0) cycle

                            d = norm2(t%node_centroid(:,m_node) - t%node_centroid(:,n_node))
                            if(d - t%node_dimension(m_node) - t%node_dimension(n_node) < min_dist_thr) then
                                if(t%is_leaf(m_node) .and. t%is_leaf(n_node)) then
                                    if(all(uncompressed_near(1:uncompressed_n_near(m_node), m_node) /= n_node)) then
                                        uncompressed_n_near(m_node) = uncompressed_n_near(m_node) + 1
                                        uncompressed_near(uncompressed_n_near(m_node), m_node) = n_node
                                    end if
                                    
                                    if(all(uncompressed_near(1:uncompressed_n_near(n_node), n_node) /= m_node)) then
                                        uncompressed_n_near(n_node) = uncompressed_n_near(n_node) + 1
                                        uncompressed_near(uncompressed_n_near(n_node), n_node) = m_node
                                    end if
                                else
                                    uncompressed_n_nnear(m_node) = uncompressed_n_nnear(m_node) + 1
                                    uncompressed_nnear(uncompressed_n_nnear(m_node), m_node) = n_node
                                end if
                            else
                                if(.not. any(uncompressed_far(1:uncompressed_n_far(m_node), m_node) == n_node)) then
                                    uncompressed_n_far(m_node) = uncompressed_n_far(m_node) + 1
                                    uncompressed_far(uncompressed_n_far(m_node), m_node) = n_node
                                end if
                                if(.not. any(uncompressed_far(1:uncompressed_n_far(n_node), n_node) == m_node)) then
                                    uncompressed_n_far(n_node) = uncompressed_n_far(n_node) + 1
                                    uncompressed_far(uncompressed_n_far(n_node), n_node) = m_node
                                end if
                            end if
                        end do
                    end do
                end do
            end do
        end do
        call time_pull("Near search")

        call time_push
        call aggregative_pass(t, uncompressed_n_far, uncompressed_far)
        call time_pull("Aggregative pass")
        
        call time_push
        !! Now compress in yale sparse format and delete the uncompressed lists
        call compress_list(t%n_nodes, uncompressed_far, uncompressed_n_far, t%far_nl)
        call compress_list(t%n_nodes, uncompressed_near, uncompressed_n_near, t%near_nl)
        call time_pull("List compression")

        deallocate(uncompressed_far, uncompressed_n_far, &
                   uncompressed_near, uncompressed_n_near, &
                   uncompressed_nnear, uncompressed_n_nnear)
        call time_pull("Far near list")
    end subroutine

    subroutine populate_leaf_list(t)
        implicit none 

        type(fmm_tree_type), intent(inout) :: t
        !! Tree data structure to populate

        integer(ip) :: i
       
        t%is_leaf(:) = .false.
        do i=1, t%n_nodes
            if(all(t%children(:,i) == 0)) t%is_leaf(i) = .true.
        end do

    end subroutine
    
    subroutine populate_level_list(t)
        use mod_adjacency_mat, only: compress_list

        implicit none 

        type(fmm_tree_type), intent(inout) :: t
        !! Tree data structure to populate

        integer(ip) :: i, lev
        integer(ip), allocatable :: tmp_levels(:,:), idx_level(:)
        
        ! There shouldn't be any level with more than t%n_particles nodes.
        allocate(tmp_levels(t%n_particles, t%breadth))
        allocate(idx_level(t%breadth))
        idx_level = 0
        ! TODO this can probably be improved for efficiency
        do i=1, t%n_nodes
            lev = t%node_level(i)
            idx_level(lev) = idx_level(lev) + 1
            tmp_levels(idx_level(lev), lev) = i
        end do

        call compress_list(t%breadth, tmp_levels, idx_level, t%level_list)

        deallocate(tmp_levels)
        deallocate(idx_level)
    end subroutine
    
end module
