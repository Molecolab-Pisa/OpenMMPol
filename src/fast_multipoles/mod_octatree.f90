module mod_octatree
    use mod_tree
    use mod_constants, only: ip, rp
    use mod_fmm_utils, only: fmm_error
    use mod_profiling

    implicit none
    private

    public :: init_as_octatree
    
contains

    subroutine init_as_octatree(t, c_particle, dfar, min_cell_size_in)
        !! Build an adaptive octatree
        use iso_c_binding, only: c_int32_t
        use mod_adjacency_mat, only: compress_list, free_yale_sparse, yale_sparse

        implicit none 

        type(fmm_tree_type), intent(inout) :: t
        !! Tree data structure to populate
        real(rp), target, intent(in) :: c_particle(:,:)
        !! Coordinates of the particles to insert in the tree
        real(rp), intent(in) :: dfar
        !! Threshold distance for near to far field
        real(rp), intent(in), optional :: min_cell_size_in
        !! Minimum radius for a node, if a node is below this threshold it 
        !! won't be split

        integer(ip) :: d, i_node, nnode_guess, max_lev, i_part, idx, pl_idx(8), iparent, ilev, j, &
                       b_level, e_level, i, npart, ii, jj, l
        integer(ip), allocatable :: t_children(:,:), particle_lst(:,:)
        real(rp), allocatable :: n_dimension(:, :), n_centers(:, :)
        real(rp) :: minc(3), maxc(3), root_dim(3), root_center(3), dist, min_cell_size
        type(yale_sparse) :: node2particle        

        call time_push

        if(present(min_cell_size_in)) then
            min_cell_size = max(0.0_rp, min_cell_size_in)
        else
            min_cell_size = 0.0_rp
        end if

        t%tree_degree = 8
        t%n_particles = size(c_particle, 2)
        t%particles_coords => c_particle
        
        ! Compute the center of the space and its dimension
        do d=1, 3
            minc = minval(c_particle(d,:))
            maxc = maxval(c_particle(d,:))
            root_dim(d) = maxc(d) - minc(d)
            root_center(d) = minc(d) + root_dim(d) * .5
        end do

        max_lev = ceiling(log(maxval(root_dim)/2.0)/log(2.0)) + 1
        
        if(ip == c_int32_t .and. max_lev >= 10) then
          call fmm_error("The code is compiled with int32 but for the &
                         &size of this system int64 are required!")
        end if

        nnode_guess = (8**(max_lev) - 1) / 7 ! Sum of a geometric seies

        allocate(t_children(8, nnode_guess))
        allocate(n_dimension(3, nnode_guess))
        allocate(n_centers(3, nnode_guess))
        allocate(particle_lst(8, t%n_particles))

        node2particle%n = nnode_guess
        allocate(node2particle%ri(node2particle%n+1))
        allocate(node2particle%ci(max_lev*t%n_particles))
      
        t_children = 0
        ! Initialize root node
        i_node = 1
        b_level = 1
        n_centers(:, i_node) = root_center(:)
        n_dimension(:, i_node) = root_dim(:)
        e_level = 1

        node2particle%ri(1) = 1
        node2particle%ri(2) = t%n_particles + 1
        do i=1, t%n_particles
            node2particle%ci(i) = i
        end do
        
        do ilev=1, max_lev
            do iparent=b_level, e_level
                ! We should stop here?
                npart = node2particle%ri(iparent+1) - node2particle%ri(iparent)
                if(npart < 2) then
                    cycle
                end if

                ! TODO: breaks far field
                if(norm2(n_dimension(:,iparent)) / 2.0 < min_cell_size) then
                    ! Splitting this cell is going to create all-near children for sure,
                    ! so we can stop here
                    cycle
                end if
             
                ! If all the tests above are passed, the node should be split
                ! Divide the particles in parent node into children nodes
                pl_idx = 0
                do i=node2particle%ri(iparent), node2particle%ri(iparent+1)-1
                    i_part = node2particle%ci(i)

                    ! Compute the index of the cell in which the particle will be assigned
                    idx = 1
                    if(c_particle(1,i_part) > n_centers(1,iparent)) idx = idx + 1
                    if(c_particle(2,i_part) > n_centers(2,iparent)) idx = idx + 2
                    if(c_particle(3,i_part) > n_centers(3,iparent)) idx = idx + 4
                    
                    pl_idx(idx) = pl_idx(idx) + 1
                    particle_lst(idx,pl_idx(idx)) = i_part
                end do

                do i=1, 8
                    if(pl_idx(i) /= 0) then
                        ! The node is not empty so it should be created
                        i_node = i_node + 1
                        t_children(i,iparent) = i_node
                        node2particle%ri(i_node+1) = node2particle%ri(i_node)
                        ! The node is not empty
                        if(pl_idx(i) == 1) then
                            ! Only one particle in this node
                            n_dimension(:,i_node) = 1.0
                            n_centers(:,i_node) = c_particle(:,particle_lst(i,1))
                            node2particle%ci(node2particle%ri(i_node+1)) = particle_lst(i,1)
                            node2particle%ri(i_node+1) = node2particle%ri(i_node+1) + 1
                        else
                            n_centers(:,i_node) = 0.0
                            minc = c_particle(:,particle_lst(i,1))
                            maxc = minc
                            do j=1, pl_idx(i)
                                node2particle%ci(node2particle%ri(i_node+1)) = particle_lst(i,j)
                                node2particle%ri(i_node+1) = node2particle%ri(i_node+1) + 1

                                do d=1, 3
                                    minc(d) = min(minc(d), c_particle(d,particle_lst(i,j)))
                                    maxc(d) = max(maxc(d), c_particle(d,particle_lst(i,j)))
                                end do
                            end do
                            n_dimension(:,i_node) = (maxc - minc)
                            n_centers(:,i_node) = minc + n_dimension(:,i_node) * .5
                        end if
                    end if
                end do
            end do

            b_level = e_level + 1
            e_level = i_node
                    
            if(e_level - b_level < 0) exit
        end do
       
        ! Populate the octatree
        t%breadth = ilev
        t%n_nodes = i_node 
        call allocate_tree(t)
        
        t%parent(1) = 0
        t%node_level(1) = 1
        t%node_centroid(:,:) = n_centers(:,:t%n_nodes)
        
        t%children(:,:) = 0
        t%particle_list%ri(1) = 1
        do i=1, t%n_nodes
            if(i > 1) t%node_level(i) = t%node_level(t%parent(i)) + 1
            !t%node_dimension(i) = maxval(n_dimension(:,i)) * 0.7072 !sqrt(2.0) / 2.0
            
            ii = 1
            do j=1, t%tree_degree
                if(t_children(j,i) > 0) then
                    t%children(ii,i) = t_children(j,i)
                    t%parent(t%children(ii,i)) = i
                    ii = ii + 1
                end if
            end do

            t%particle_list%ri(i+1) = t%particle_list%ri(i)
            if(all(t%children(:,i) == 0)) then
                ! Node is a leaf, list all particles in this node
                do j=node2particle%ri(i), node2particle%ri(i+1)-1
                    t%particle_list%ci(t%particle_list%ri(i+1)) = node2particle%ci(j)
                    t%particle_list%ri(i+1) = t%particle_list%ri(i+1) + 1
                end do
            end if
        end do

        call populate_level_list(t)
        call populate_leaf_list(t)
        
        do i=1, t%n_nodes
            if(all(t%children(:,i) == 0)) then
                t%node_dimension(i) = norm2(n_dimension(:,i)) / 2.0
            end if 
        end do

        do i=t%breadth, 1, -1
            do jj=t%level_list%ri(i), t%level_list%ri(i+1)-1
                j = t%level_list%ci(jj)
                if((any(t%children(:,j) /= 0))) then
                    do ii=1, t%tree_degree
                        l = t%children(ii,j)
                        if(l == 0)  cycle
                        dist = norm2(t%node_centroid(:,j)- t%node_centroid(:,l))
                        t%node_dimension(j) = max(t%node_dimension(j), dist + t%node_dimension(l))
                    end do 
                end if
            end do
        end do

        call tree_populate_farnear_lists(t, dfar)
        
        t%particle_to_node = 0
        do i=1, t%n_nodes
            if(t%children(1,i) /= 0) cycle
            ! The node is a leaf
            do j=t%particle_list%ri(i), t%particle_list%ri(i+1) - 1
                t%particle_to_node(t%particle_list%ci(j)) = i
            end do
        end do

        do i=1, t%n_particles
            if(t%particle_to_node(i) == 0) then
                call fmm_error("Some particles aren't in any node")
            end if
        end do
        call time_pull('Octatree initialization')
    end subroutine
    
end module
