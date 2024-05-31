#define _x_ 1
#define _y_ 2
#define _z_ 3

#define _xx_ 1
#define _xy_ 2
#define _yy_ 3
#define _xz_ 4
#define _yz_ 5
#define _zz_ 6
#define _yx_ _xy_
#define _zx_ _xz_
#define _zy_ _yz_

#define _xxx_ 1
#define _xxy_ 2
#define _xxz_ 3
#define _xyy_ 4
#define _xyz_ 5
#define _xzz_ 6
#define _yyy_ 7
#define _yyz_ 8
#define _yzz_ 9
#define _zzz_ 10
#define _xyx_ _xxy_
#define _xzx_ _xxz_
#define _xzy_ _xyz_
#define _yxx_ _xxy_
#define _yxy_ _xyy_
#define _yxz_ _xyz_
#define _yyx_ _xyy_
#define _yzx_ _xyz_
#define _yzy_ _yyz_
#define _zxx_ _xxz_
#define _zxy_ _xyz_
#define _zxz_ _xzz_
#define _zyx_ _xyz_
#define _zyy_ _yyz_
#define _zyz_ _yzz_
#define _zzx_ _xzz_
#define _zzy_ _zyz_

module mod_fmm

    use mod_constants, only : ip, rp
    use mod_tree, only: fmm_tree_type
    use mod_fmm_utils, only: fmm_error
    use mod_profiling, only: time_pull, time_push
    use mod_harmonics, only: prepare_fmmm_constants

    implicit none

    type fmm_type
        type(fmm_tree_type), pointer :: tree
        !! Tree data structure to store the particles
        integer(ip) :: pmax_mm
        !! Maximum order of spherical harmonics used in multipolar expansion
        integer(ip) :: pmax_le
        !! Maximum order of spherical harmonics used in local expansion
        real(rp), allocatable :: multipoles_p(:,:)
        !! Multipole expansion for each particle
        real(rp), allocatable :: multipoles(:,:)
        !! Multipole expansion for each node of the tree
        real(rp), allocatable :: local_expansion(:, :)
        !! Local expansion for each node of the tree
    end type

    contains

    subroutine fmm_init(fmm_obj, pmax, tree)
        use mod_fmm_utils, only: ntot_sph_harm

        implicit none

        integer(ip), intent(in) :: pmax
        type(fmm_type), intent(inout) :: fmm_obj
        type(fmm_tree_type), intent(in), target :: tree

        call time_push()
        fmm_obj%tree => tree
        fmm_obj%pmax_mm = pmax
        fmm_obj%pmax_le = pmax

        allocate(fmm_obj%multipoles_p(ntot_sph_harm(fmm_obj%pmax_mm), tree%n_particles))
        allocate(fmm_obj%multipoles(ntot_sph_harm(fmm_obj%pmax_mm), tree%n_nodes))
        allocate(fmm_obj%local_expansion(ntot_sph_harm(fmm_obj%pmax_le), tree%n_nodes))
        
        fmm_obj%multipoles = 0.0
        fmm_obj%local_expansion = 0.0
        call prepare_fmmm_constants(fmm_obj%pmax_mm, fmm_obj%pmax_le)
        call time_pull("FMM Init")
    end subroutine

    subroutine free_fmm(fmm_obj)
        use mod_tree, only: free_tree

        implicit none

        type(fmm_type), intent(inout) :: fmm_obj

        ! call free_tree(fmm_obj%tree)
        if(allocated(fmm_obj%multipoles_p)) deallocate(fmm_obj%multipoles_p)
        if(allocated(fmm_obj%multipoles)) deallocate(fmm_obj%multipoles)
        if(allocated(fmm_obj%local_expansion)) deallocate(fmm_obj%local_expansion)
    end subroutine

    subroutine fmm_solve(fmm_obj)
        implicit none

        type(fmm_type), intent(inout) :: fmm_obj

        call tree_m2m(fmm_obj)
        call tree_m2l(fmm_obj)
        call tree_l2l(fmm_obj)

    end subroutine

    subroutine cart_prop_at_ipart(fmm_obj, i_part, do_V, V, do_E, E, do_grdE, grdE, do_HE, HE)
        implicit none

        type(fmm_type), intent(in) :: fmm_obj
        integer(ip) :: i_part
        logical, intent(in) :: do_V, do_E, do_grdE, do_HE
        real(rp), intent(inout) :: V, E(3), grdE(6), HE(10)

        call cart_propfar_at_ipart(fmm_obj, i_part, do_V, V, do_E, E, do_grdE, grdE, do_HE, HE)
        call cart_propnear_at_ipart(fmm_obj, i_part, do_V, V, do_E, E, do_grdE, grdE, do_HE, HE)
    end subroutine
    
    subroutine cart_propfar_at_ipart(fmm_obj, i_part, do_V, V, do_E, E, do_grdE, grdE, do_HE, HE)
        use mod_constants, only: pi
        use mod_fmm_utils, only: ntot_sph_harm
        use mod_harmonics, only: fmm_l2l
        implicit none

        type(fmm_type), intent(in) :: fmm_obj
        integer(ip) :: i_part
        logical, intent(in) :: do_V, do_E, do_grdE, do_HE
        real(rp), intent(inout) :: V, E(3), grdE(6), HE(10)
        
        type(fmm_tree_type), pointer :: t
        integer(ip) :: i_node
        real(rp) :: x2_y2, z2, x2z_y2z, z3, xz2, yz2, x3_3xy2, y3_3x2y, dr(3)
        real(rp), allocatable :: tmp_local(:)
        
        t => fmm_obj%tree

        i_node = t%particle_to_node(i_part)
        allocate(tmp_local(ntot_sph_harm(fmm_obj%pmax_le)))
        tmp_local = 0.0

        dr = t%node_centroid(:,i_node) - t%particles_coords(:,i_part)
        ! Local expansion needs a further translation
        call fmm_l2l(dr, &
                     1.0_rp, 1.0_rp, &
                     fmm_obj%pmax_le, fmm_obj%local_expansion(:,i_node), &
                     tmp_local)

        if(do_V) then
            v = v + sqrt(4.0*pi) * tmp_local(1)
        end if

        if(do_E) then
            E(3) = E(3) - sqrt(4.0/3.0*pi) * tmp_local(3) 
            E(1) = E(1) - sqrt(4.0/3.0*pi) * tmp_local(4) 
            E(2) = E(2) - sqrt(4.0/3.0*pi) * tmp_local(2)
        end if

        if(do_grdE) then
            x2_y2 = sqrt(16.0*pi/15.0) * tmp_local(9) * 3.0
            z2 = (sqrt(16.0*pi/5.0) * tmp_local(7)) 
            grdE(6) = grdE(6) + z2 ! zz
            grdE(1) = grdE(1) + (x2_y2 - z2) / 2.0
            grdE(3) = grdE(3) - (x2_y2 + z2) / 2.0
            grdE(2) = grdE(2) + 3.0 * sqrt(4.0*pi/15.0) * tmp_local(5) !xy
            grdE(4) = grdE(4) + 3.0 * sqrt(4.0*pi/15.0) * tmp_local(8) !xz
            grdE(5) = grdE(5) + 3.0 * sqrt(4.0*pi/15.0) * tmp_local(6) !yz
        end if

        if(do_HE) then
            z3 = 15.0 * 4.0 / 5.0 * sqrt(pi / 7.0) *            tmp_local(13)
            x2z_y2z = 15.0 * 4.0 * sqrt(pi / 105.0) *           tmp_local(15)
            xz2 = 15.0 * 4.0 / 5.0 * sqrt(2.0 * pi / 21.0) *    tmp_local(14)
            yz2 = 15.0 * 4.0 / 5.0 * sqrt(2.0 * pi / 21.0) *    tmp_local(12)
            x3_3xy2 = 15.0 * 4.0 * sqrt(2.0 * pi / 35.0) *      tmp_local(16)
            y3_3x2y = - 15.0 * 4.0 * sqrt(2.0 * pi / 35.0) *    tmp_local(10)
            HE(_xyz_) = HE(_xyz_) - 15.0 * sqrt(4.0*pi/105.0) * tmp_local(11)
            HE(_yzz_) = HE(_yzz_) - yz2 
            HE(_xzz_) = HE(_xzz_) - xz2 
            HE(_zzz_) = HE(_zzz_) - z3
            HE(_xxz_) = HE(_xxz_) - (x2z_y2z - z3) / 2.0
            HE(_yyz_) = HE(_yyz_) + (x2z_y2z + z3) / 2.0
            HE(_xxx_) = HE(_xxx_) - (x3_3xy2 - 3.0 * xz2) / 4.0
            HE(_yyy_) = HE(_yyy_) - (y3_3x2y - 3.0 * yz2) / 4.0
            HE(_xyy_) = HE(_xyy_) + (x3_3xy2 + xz2) / 4.0
            HE(_xxy_) = HE(_xxy_) + (y3_3x2y + yz2) / 4.0
        end if
    deallocate(tmp_local)
    end subroutine
    
    subroutine cart_propnear_at_ipart(fmm_obj, i_part, do_V, V, do_E, E, do_grdE, grdE, do_HE, HE)
        use mod_constants, only: pi
        use mod_fmm_utils, only: ntot_sph_harm
        use mod_harmonics, only: fmm_m2l
        implicit none

        type(fmm_type), intent(in) :: fmm_obj
        integer(ip) :: i_part
        logical, intent(in) :: do_V, do_E, do_grdE, do_HE
        real(rp), intent(inout) :: V, E(3), grdE(6), HE(10)

        real(rp), allocatable :: local_tmp(:), local(:)
        type(fmm_tree_type), pointer :: t
        integer(ip) :: i_node, j, j_node, j_particle, jj
        real(rp) :: c_st(3), x2_y2, z2, x2z_y2z, z3, xz2, yz2, x3_3xy2, y3_3x2y
        t => fmm_obj%tree
        
        i_node = t%particle_to_node(i_part)
        
        allocate(local_tmp(ntot_sph_harm(fmm_obj%pmax_le)))
        allocate(local(ntot_sph_harm(fmm_obj%pmax_le)))
       
        local = 0.0

        do j=t%near_nl%ri(i_node), t%near_nl%ri(i_node+1)-1
            j_node = t%near_nl%ci(j)
            do jj=t%particle_list%ri(j_node), t%particle_list%ri(j_node+1)-1
                j_particle = t%particle_list%ci(jj)
                if(i_part == j_particle) cycle
                c_st = t%particles_coords(:,j_particle) - t%particles_coords(:,i_part)
                
                call fmm_m2l(c_st, &
                                fmm_obj%pmax_mm, &
                                fmm_obj%pmax_le, &
                                fmm_obj%multipoles_p(:,j_particle), & 
                                local_tmp)
                    local = local + local_tmp
            end do
        end do
        
        deallocate(local_tmp)

        if(do_V) then
            v = v + sqrt(4.0*pi) * local(1)
        end if

        if(do_E) then
            E(3) = E(3) - sqrt(4.0/3.0*pi) * local(3) 
            E(1) = E(1) - sqrt(4.0/3.0*pi) * local(4) 
            E(2) = E(2) - sqrt(4.0/3.0*pi) * local(2)
        end if
    
        if(do_grdE) then
            x2_y2 = sqrt(16.0*pi/15.0) * local(9) * 3.0
            z2 = sqrt(16.0*pi/5.0) * local(7)
            grdE(6) = grdE(6) + z2 ! zz
            grdE(1) = grdE(1) + (x2_y2 - z2) / 2.0
            grdE(3) = grdE(3) - (x2_y2 + z2) / 2.0
            grdE(2) = grdE(2) + 3.0 * sqrt(4.0*pi/15.0) * local(5) !xy
            grdE(4) = grdE(4) + 3.0 * sqrt(4.0*pi/15.0) * local(8) !xz
            grdE(5) = grdE(5) + 3.0 * sqrt(4.0*pi/15.0) * local(6) !yz
        end if
        
        if(do_HE) then
            z3 = 15.0 * 4.0 / 5.0 * sqrt(pi / 7.0) * local(13)
            x2z_y2z = 15.0 * 4.0 * sqrt(pi / 105.0) * local(15)
            xz2 = 15.0 * 4.0 / 5.0 * sqrt(2.0 * pi / 21.0) * local(14)
            yz2 = 15.0 * 4.0 / 5.0 * sqrt(2.0 * pi / 21.0) * local(12)
            x3_3xy2 = 15.0 * 4.0 * sqrt(2.0 * pi / 35.0) * local(16)
            y3_3x2y = - 15.0 * 4.0 * sqrt(2.0 * pi / 35.0) * local(10)
            HE(_xyz_) = HE(_xyz_) - 15.0 * sqrt(4.0*pi/105.0) * local(11)
            HE(_yzz_) = HE(_yzz_) - yz2 
            HE(_xzz_) = HE(_xzz_) - xz2 
            HE(_zzz_) = HE(_zzz_) - z3
            HE(_xxz_) = HE(_xxz_) - (x2z_y2z - z3) / 2.0
            HE(_yyz_) = HE(_yyz_) + (x2z_y2z + z3) / 2.0
            HE(_xxx_) = HE(_xxx_) - (x3_3xy2 - 3.0 * xz2) / 4.0
            HE(_yyy_) = HE(_yyy_) - (y3_3x2y - 3.0 * yz2) / 4.0
            HE(_xyy_) = HE(_xyy_) + (x3_3xy2 + xz2) / 4.0
            HE(_xxy_) = HE(_xxy_) + (y3_3x2y + yz2) / 4.0
        end if
        deallocate(local)

    end subroutine
    
    subroutine tree_p2m(fmm_obj, particle_multipoles, pmax_particles)
        use mod_fmm_utils, only: ntot_sph_harm
        use mod_harmonics, only: fmm_m2m

        implicit none

        type(fmm_type), intent(inout) :: fmm_obj
        real(rp), intent(in) :: particle_multipoles(:,:)
        integer(ip), intent(in) :: pmax_particles

        type(fmm_tree_type), pointer :: t
        integer(ip) :: j, i_node, j_particle, n_expansion, p_expansion
        real(rp), allocatable :: expansion_s(:), expansion_t(:)
        real(rp) :: c_st(3), r_s, r_t

        call time_push

        t => fmm_obj%tree

        ! Sources are particles and target are nodes, in princile they could have different
        ! expansion levels (in general we should expect pmax_mm >> pmax_particles, but 
        ! people are strange).
        
        if(fmm_obj%pmax_mm < pmax_particles) then
            call fmm_error("Multipoles expansion should be at least as high as particles expansion")
        end if
        
        fmm_obj%multipoles_p(:,:) = 0.0
        fmm_obj%multipoles_p(1_ip:ntot_sph_harm(pmax_particles),:) = particle_multipoles(:,:)

        p_expansion = max(fmm_obj%pmax_mm, pmax_particles)
        n_expansion = ntot_sph_harm(p_expansion)
        allocate(expansion_s(n_expansion), expansion_t(n_expansion))
      
        !$omp parallel do default(shared) private(i_node, expansion_s, expansion_t, j, j_particle, c_st, r_t, r_s) schedule(dynamic)
        do i_node=1, t%n_nodes
            ! For each node
            fmm_obj%multipoles(:,i_node) = 0.0
            expansion_s = 0.0 ! Needed if pmax_mm > particles_pmax
            do j=t%particle_list%ri(i_node), t%particle_list%ri(i_node+1)-1
                j_particle = t%particle_list%ci(j)
                expansion_s(1:ntot_sph_harm(pmax_particles)) = particle_multipoles(:, j_particle)
                c_st = t%particles_coords(:,j_particle) - t%node_centroid(:,i_node)
                call fmm_m2m(c_st, p_expansion, expansion_s, expansion_t)
                fmm_obj%multipoles(:,i_node) = fmm_obj%multipoles(:,i_node) + expansion_t
            end do
        end do

        call time_pull('Tree P2M')
    end subroutine
    
    subroutine tree_m2m(fmm_obj)
        use mod_fmm_utils, only: ntot_sph_harm
        use mod_harmonics, only: fmm_m2m

        implicit none

        type(fmm_type), intent(inout) :: fmm_obj

        type(fmm_tree_type), pointer :: t
        integer(ip) :: l, i, j, i_node, j_node
        real(rp), allocatable :: expansion_s(:), expansion_t(:)
        real(rp) :: c_st(3), r_s, r_t

         call time_push
        t => fmm_obj%tree
        allocate(expansion_s(ntot_sph_harm(fmm_obj%pmax_mm)), &
                 expansion_t(ntot_sph_harm(fmm_obj%pmax_mm)))


        do l=t%breadth, 1, -1
            ! For each level, leaves to root
            !$omp parallel do default(shared) private(i, j, i_node, expansion_s, j_node, c_st, r_s, r_t, expansion_t)
            do i=t%level_list%ri(l), t%level_list%ri(l+1)-1
                ! For each node in the level
                i_node = t%level_list%ci(i)

                if(t%children(1,i_node) /= 0) then
                    ! P2M should already have populated leaves nodes
                    ! This node is not a leaf
                    fmm_obj%multipoles(:,i_node) = 0.0 ! Initialization
                    expansion_s = 0.0 ! Needed if pmax_mm < particles_pmax (unreasonable)
                    do j=1, t%tree_degree
                        j_node = t%children(j,i_node)
                        if(j_node == 0) then
                            ! This child is not present
                            cycle
                        end if

                        expansion_s = fmm_obj%multipoles(:, j_node)
                        c_st = t%node_centroid(:,j_node) - t%node_centroid(:,i_node)
                        call fmm_m2m(c_st, fmm_obj%pmax_mm, expansion_s, expansion_t)
                        fmm_obj%multipoles(:,i_node) = fmm_obj%multipoles(:,i_node) + expansion_t
                    end do
                end if
            end do
        end do

        deallocate(expansion_s, expansion_t)
        call time_pull('Tree M2M')
    end subroutine

    subroutine tree_m2l(fmm_obj)
        use mod_fmm_utils, only: ntot_sph_harm
        use mod_harmonics, only: fmm_m2l

        implicit none

        type(fmm_type), intent(inout) :: fmm_obj

        type(fmm_tree_type), pointer :: t
        real(rp) :: c_st(3), r_s, r_t
        real(rp), allocatable :: mme_s(:), le_t(:)
        integer(ip) :: i_node, j_node, j
       
        call time_push
        t => fmm_obj%tree
        allocate(mme_s(ntot_sph_harm(fmm_obj%pmax_mm)))
        allocate(le_t(ntot_sph_harm(fmm_obj%pmax_le)))
       
        !$omp parallel do default(shared) &
        !$omp private(i_node, j_node, c_st, r_s, r_t, mme_s, le_t)
        do i_node = 1, t%n_nodes
            fmm_obj%local_expansion(:,i_node) = 0.0
            do j=t%far_nl%ri(i_node), t%far_nl%ri(i_node+1)-1
                j_node = t%far_nl%ci(j)
                
                c_st = t%node_centroid(:,j_node) - t%node_centroid(:,i_node)
                mme_s = fmm_obj%multipoles(:,j_node)
                
                call fmm_m2l(c_st, fmm_obj%pmax_mm, fmm_obj%pmax_le, mme_s, le_t)
                fmm_obj%local_expansion(:,i_node) = fmm_obj%local_expansion(:,i_node) + le_t
            end do
        end do

        deallocate(mme_s, le_t)
        call time_pull('Tree M2L')
    end subroutine

    subroutine tree_l2l(fmm_obj)
        use mod_fmm_utils, only: ntot_sph_harm
        use mod_harmonics, only: fmm_l2l

        implicit none

        type(fmm_type), intent(inout) :: fmm_obj

        type(fmm_tree_type), pointer :: t
        integer(ip) :: l, i, j, i_node, j_node
        real(rp), allocatable :: le_s(:), le_t(:)
        real(rp) :: c_st(3), r_s, r_t

        call time_push
        t => fmm_obj%tree
        
        allocate(le_s(ntot_sph_harm(fmm_obj%pmax_le)), &
                 le_t(ntot_sph_harm(fmm_obj%pmax_le)))

        do l=1, t%breadth-1
            ! For each level, root to leaves
            !$omp parallel do default(shared) private(i_node, le_s, j, i, r_s, j_node, c_st, r_t, le_t)
            do i=t%level_list%ri(l), t%level_list%ri(l+1)-1
                ! For each node in the level

                ! Propagate local expansion on each node on its children
                i_node = t%level_list%ci(i)
                le_s = fmm_obj%local_expansion(:, i_node)

                ! If node is not a leaf
                do j=1, t%tree_degree
                    j_node = t%children(j,i_node)
                    if(j_node == 0) then
                        ! This child is not present
                        cycle
                    end if

                    c_st = t%node_centroid(:,i_node) - t%node_centroid(:,j_node)
                    !call time_push()
                    call fmm_l2l(c_st, 1.0_rp, 1.0_rp, fmm_obj%pmax_le, le_s, le_t)
                    !call time_pull("Single L2L")
                    fmm_obj%local_expansion(:,j_node) = fmm_obj%local_expansion(:,j_node) + le_t
                end do
            end do
        end do
        
        deallocate(le_s, le_t)
        call time_pull('Tree L2L')
    end subroutine

    
end module
