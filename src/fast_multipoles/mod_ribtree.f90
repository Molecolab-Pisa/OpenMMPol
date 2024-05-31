module mod_ribtree
    use mod_tree
    use mod_constants, only: ip, rp
    use mod_fmm_utils, only: fmm_error
    use mod_profiling

    implicit none
    private

    public :: init_as_ribtree
    
contains

    subroutine init_as_ribtree(t, c_particle, dfar)
        !! Build a recursive inertial binary tree
        !!
        !! Uses inertial bisection in a recursive manner until each leaf node has only
        !! one particle inside. Number of tree nodes is always 2*n_particle-1.
        use mod_adjacency_mat, only: compress_list, free_yale_sparse

        implicit none 

        type(fmm_tree_type), intent(inout) :: t
        !! Tree data structure to populate
        real(rp), target, intent(in) :: c_particle(:,:)
        !! Coordinates of the particles to insert in the tree
        real(rp), intent(in) :: dfar
        !! Minimum size of a node

        integer(ip) :: i_node, j_node, i_particle, i, j, s, e, n, div
        integer(ip), allocatable :: order(:), cluster(:,:)
        real(rp) :: r1, r2, c(3), c1(3), c2(3), d
        
        call time_push

        t%tree_degree = 2
        t%n_particles = size(c_particle, 2)
        t%particles_coords => c_particle
        t%n_nodes = 2 * t%n_particles  - 1
        ! The number of levels is the log2(n_particles) approximated by excess;
        ! this is an educated guess, and it is going to be likely higher, so at each
        ! iteration the level should be checked and  in case level list should be 
        ! enlarged accordingly
        t%breadth = ceiling(log(real(t%n_particles))/log(2.0)) 

        call allocate_tree(t)
        
        ! Init the root node
        t%parent(1) = 0
        t%node_level(1) = 1
        t%particle_list%ri(1) = 1

        ! Index of the first unassigned node
        i_node = 2

        ! Init particle ordering
        allocate(order(t%n_particles))
        do i = 1, t%n_particles
            order(i) = i
        end do
        ! Init cluster bipartition
        allocate(cluster(2,t%n_nodes))
        cluster(1,1) = 1
        cluster(2,1) = t%n_particles

        t%particle_to_node = 0

        ! Divide nodes until there is just a single particle per (leaf) node
        do i = 1, t%n_nodes
            s = cluster(1, i)
            e = cluster(2, i)
            n = e - s + 1
            div = 0
            ! Divide only if there are 2 or more particles
            if (n > 1) then
                ! Use inertial bisection to reorder particles and cut into the 
                ! particles below this node into two halves
                call tree_rib_node_bisect(c_particle, order(s:e), div)
                
                ! Since this is not a leaf node, it have no particles
                t%particle_list%ri(i+1) = t%particle_list%ri(i)

                ! Assign the first half
                cluster(1, i_node) = s
                cluster(2, i_node) = s + div - 1
                ! Assign the second half to the (j+1)-th node
                cluster(1, i_node+1) = s + div
                cluster(2, i_node+1) = e
                ! Update list of children of i-th node
                t%children(1, i) = i_node
                t%children(2, i) = i_node + 1
                ! Set parents of new nodes
                t%parent(i_node) = i
                t%parent(i_node+1) = i
                ! Set the level for the newly created nodes
                t%node_level(i_node) = t%node_level(i) + 1
                t%node_level(i_node+1) = t%node_level(i) + 1
                
                ! Shift index of the first unassigned node
                i_node = i_node + 2

            ! Set information for a leaf node
            else
                ! This is a leaf node that means:
                ! a) no children nodes
                t%children(:, i) = 0
                ! b) it contains a single particle
                t%particle_list%ri(i+1) = t%particle_list%ri(i) + 1
                t%particle_list%ci(t%particle_list%ri(i)) = order(s)
                t%particle_to_node(order(s)) = i
            end if
        end do

        if(any(t%particle_to_node == 0)) then
            ! Sanity check
            call fmm_error("All particles should be assigned to a node, this is a bug.")
        end if

        ! Update levels
        t%breadth = maxval(t%node_level)
        call populate_level_list(t)
        call populate_leaf_list(t)

        !! Compute centroids
        ! Bottom-to-top pass over all nodes of the tree, level by level
        do i=t%breadth, 1, -1
            ! For each level in the tree
            do j=t%level_list%ri(i), t%level_list%ri(i+1)-1
                ! For each node in the level (TODO parallelize here)
                i_node = t%level_list%ci(j)
                if(t%particle_list%ri(i_node) < t%particle_list%ri(i_node+1)) then
                    ! This is a leaf node
                    ! We know that each leaf node have one particle only.
                    i_particle = t%particle_list%ci(t%particle_list%ri(i_node))
                    t%node_centroid(:,i_node) = c_particle(:,i_particle)
                    ! Node dimension
                    t%node_dimension(i_node) = 0.0
                else
                    ! The centroid is the average of the children's centroids
                    j_node = t%children(1, i_node)
                    r1 = t%node_dimension(j_node)
                    c1 = t%node_centroid(:,j_node)

                    j_node = t%children(2, i_node)
                    r2 = t%node_dimension(j_node)
                    c2 = t%node_centroid(:,j_node)

                    c = c1 - c2
                    d = norm2(c)
                    ! If sphere #2 is completely inside sphere #1 use the first sphere
                    if(r1 >= r2 + d) then
                        t%node_centroid(:,i_node) = c1
                        t%node_dimension(i_node) = r1
                    ! If sphere #1 is completely inside sphere #2 use the second sphere
                    else if(r2 >= r1+d) then
                        t%node_centroid(:,i_node) = c2
                        t%node_dimension(i_node) = r2
                    ! Otherwise use special formula to find a minimal sphere
                    else
                        t%node_dimension(i_node) = (r1+r2+d) * 0.5
                        t%node_centroid(:,i_node) = c2 + c * (t%node_dimension(i_node)-r2) / d
                    end if
                    
                end if
            end do    
        end do

        deallocate(order)
        deallocate(cluster)

        call tree_populate_farnear_lists(t, dfar)

        call time_pull('RIB tree initialization')
    end subroutine
    
    subroutine tree_rib_node_bisect(c_particle, order, div)
        !> Divide given cluster of spheres into two subclusters by inertial bisection
        !!
        !! @param[in] nsph: Number of all input spheres
        !! @param[in] csph: Centers of all input spheres
        !! @param[in] n: Number of spheres in a given cluster
        !! @param[inout] order: Indexes of spheres in a given cluster. On exit, indexes
        !!      `order(1:div)` correspond to the first subcluster and indexes
        !!      `order(div+1:n)` correspond to the second subcluster.
        !! @param[out] div: Break point of `order` array between two clusters.
        !! @param[inout] ddx_error: ddX error
        implicit none

        real(rp), intent(in) :: c_particle(:, :)
        integer(ip), intent(inout) :: order(:)
        integer(ip), intent(out) :: div
        
        real(rp) :: c(3),  s(3)
        real(rp), allocatable :: tmp_c(:, :), work(:)
        external :: dgesvd
        integer(ip) :: i, l, r, lwork, info, n, n_particle
        integer(ip), allocatable :: tmp_order(:)


        n = size(order)
        n_particle = size(c_particle, 2)

        allocate(tmp_c(3, n))
        allocate(tmp_order(n))

        ! Get coordinates in a contiguous array
        do i = 1, n
            tmp_c(:, i) = c_particle(:, order(i))
        end do
        ! Remove the geometrical center for the coordinates
        c = sum(tmp_c, 2) / n
        do i=1, n
            tmp_c(:,i) = tmp_c(:,i) - c
        end do
        
        !! Find right singular vectors
        ! Get proper size of temporary workspace
        lwork = -1
        call dgesvd('N', 'O', 3, n, tmp_c, 3, s, tmp_c, 3, tmp_c, 3, &
                    s, lwork, info)
        lwork = int(s(1))
        allocate(work(lwork))
        ! Get right singular vectors
        call dgesvd('N', 'O', 3, n, tmp_c, 3, s, tmp_c, 3, tmp_c, 3, &
                    work, lwork, info)
        if (info.ne.0) then
            call fmm_error('DGESVD failed in tree_node_bisect')
            return
        end if
        deallocate(work)
        !! Sort spheres by sign of the above scalar product, which is equal to
        !! the leading right singular vector scaled by the leading singular value.
        !! However, we only care about sign, so we take into account only the
        !! leading right singular vector.
        ! First empty index from the left
        l = 1
        ! First empty index from the right
        r = n
        ! Cycle over values of the singular vector
        do i = 1, n
            ! Positive scalar products are moved to the beginning of `order`
            if (tmp_c(1, i) > 0.0) then
                tmp_order(l) = order(i)
                l = l + 1
            ! Negative scalar products are moved to the end of `order`
            else
                tmp_order(r) = order(i)
                r = r - 1
            end if
        end do
        ! Set divider and update order
        div = r
        order = tmp_order
        deallocate(tmp_c)
        deallocate(tmp_order)
    end subroutine tree_rib_node_bisect
    
end module