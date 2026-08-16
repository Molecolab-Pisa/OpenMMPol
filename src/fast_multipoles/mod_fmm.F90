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

#define _xxxx_ 1
#define _xxxy_ 2
#define _xxxz_ 3
#define _xxyy_ 4
#define _xxyz_ 5
#define _xxzz_ 6
#define _xyyy_ 7
#define _xyyz_ 8
#define _xyzz_ 9
#define _xzzz_ 10
#define _yyyy_ 11
#define _yyyz_ 12
#define _yyzz_ 13
#define _yzzz_ 14
#define _zzzz_ 15

module mod_fmm

    use mod_constants, only : ip, rp
    use mod_tree, only: fmm_tree_type
    use mod_fmm_utils, only: fmm_error
    use mod_harmonics, only: prepare_fmmm_constants

    implicit none

    type fmm_type
        type(fmm_tree_type), pointer :: tree
        !! Tree data structure to store the particles
        integer(ip) :: pmax_mm
        !! Maximum order of spherical harmonics used in multipolar expansion
        integer(ip) :: pmax_le
        !! Maximum order of spherical harmonics used in local expansion
        integer(ip) :: nrhs
        !! Number of independent right-hand-sides (source sets) batched
        !! together through the SAME tree/geometry in one P2M->M2M->M2L->L2L
        !! pass. nrhs=1 (the default) is the single-source case used
        !! everywhere today; cart_prop_at_ipart and friends only support
        !! nrhs=1 objects for now (see their own guards).
        real(rp), allocatable :: multipoles_p(:,:,:)
        !! Multipole expansion for each particle, shape (ncoef, nrhs, n_particles)
        real(rp), allocatable :: multipoles(:,:,:)
        !! Multipole expansion for each node of the tree, shape (ncoef, nrhs, n_nodes)
        real(rp), allocatable :: local_expansion(:,:,:)
        !! Local expansion for each node of the tree, shape (ncoef, nrhs, n_nodes)

        logical :: m2m_rotcache_active = .false.
        !! Whether the M2M OXZ rotation-matrix cache below is in use. Only
        !! worth enabling when tree_m2m is going to be called many times in
        !! a row on this SAME geometry (e.g. once per CG iteration of a
        !! polarization solve) -- positions (and hence rotations) must not
        !! change while the cache is active, see fmm_m2m_rotcache_try_enable.
        logical :: m2m_rotcache_built = .false.
        !! Whether m2m_rotcache_fwd/bwd already hold a valid build (set by
        !! tree_m2m after its first pass with the cache active)
        real(rp), allocatable :: m2m_rotcache_fwd(:,:), m2m_rotcache_bwd(:,:)
        !! Per-node OXZ rotation-matrix stacks (see oxz_rot_cache_size),
        !! shape (oxz_rot_cache_size(pmax_mm), n_nodes). Slot j_node holds
        !! the rotation data for the M2M edge translating node j_node's
        !! multipole into its parent -- every non-root node has exactly one
        !! parent edge, so node index alone is a safe, race-free cache slot.

        logical :: m2l_rotcache_active = .false.
        !! Same idea as m2m_rotcache_active, for M2L far-list pairs.
        logical :: m2l_rotcache_built = .false.
        real(rp), allocatable :: m2l_rotcache_fwd(:,:), m2l_rotcache_bwd(:,:)
        !! Per-far-pair OXZ rotation-matrix stacks, shape
        !! (oxz_rot_cache_size(pmax_mm or pmax_le), n_far_pairs). Slot j is
        !! the position of a given far-pair within tree%far_nl's flat CSR
        !! arrays (tree%far_nl%ci/ri) -- every far pair has exactly one such
        !! position, and tree_m2l's parallel loop only ever touches disjoint
        !! j-ranges per i_node (the row structure of a CSR matrix), so this
        !! is a safe, race-free cache slot too.

        logical :: l2l_rotcache_active = .false.
        !! Same idea as m2m_rotcache_active, for L2L parent->child edges.
        logical :: l2l_rotcache_built = .false.
        real(rp), allocatable :: l2l_rotcache_fwd(:,:), l2l_rotcache_bwd(:,:)
        !! Per-node OXZ rotation-matrix stacks, shape
        !! (oxz_rot_cache_size(pmax_le), n_nodes). Slot j_node holds the
        !! rotation data for the L2L edge shifting the parent's local
        !! expansion down into child j_node -- same reasoning as
        !! m2m_rotcache_fwd/bwd, just the opposite direction along the same
        !! parent/child edges.

        logical :: p2m_rotcache_active = .false.
        !! Same idea as m2m_rotcache_active, for P2M per-particle shifts
        !! (tree_p2m calls fmm_m2m once per particle, shifting it into its
        !! leaf node -- same routine, same rotation cost, as M2M).
        logical :: p2m_rotcache_built = .false.
        real(rp), allocatable :: p2m_rotcache_fwd(:,:), p2m_rotcache_bwd(:,:)
        !! Per-particle OXZ rotation-matrix stacks, shape
        !! (oxz_rot_cache_size(pmax_mm), n_particles). Slot j_particle is
        !! the particle's own global index -- every particle belongs to
        !! exactly one leaf node's particle_list, so this is a safe,
        !! race-free cache slot under tree_p2m's parallel loop (parallelized
        !! over i_node, with disjoint particle_list ranges per node).
        !! Only used when tree_p2m is called with pmax_particles such that
        !! max(pmax_mm, pmax_particles) == pmax_mm (always true in current
        !! usage, see tree_p2m) -- falls back to uncached otherwise.

        logical :: l2p_rotcache_active = .false.
        !! Same idea as p2m_rotcache_active, for the L2P shift
        !! (cart_propfar_at_ipart calls fmm_l2l once per particle, shifting
        !! the local expansion from its leaf node's centroid down onto the
        !! exact particle position -- same routine, same rotation cost, as
        !! L2L). Unlike the tree_* routines, cart_propfar_at_ipart is
        !! called once PER PARTICLE from an external loop (field_extD2D and
        !! friends), so it cannot itself decide when "the whole build pass"
        !! is done -- l2p_rotcache_built is instead flipped by the caller
        !! after its own per-particle loop completes (see field_extD2D).
        logical :: l2p_rotcache_built = .false.
        real(rp), allocatable :: l2p_rotcache_fwd(:,:), l2p_rotcache_bwd(:,:)
        !! Per-particle OXZ rotation-matrix stacks, shape
        !! (oxz_rot_cache_size(pmax_le), n_particles), keyed by the
        !! particle's own global index -- same reasoning as p2m's cache
        !! slot.
    end type

    contains

    subroutine fmm_init(fmm_obj, pmax, tree, nrhs)
        use mod_fmm_utils, only: ntot_sph_harm

        implicit none

        integer(ip), intent(in) :: pmax
        type(fmm_type), intent(inout) :: fmm_obj
        type(fmm_tree_type), intent(in), target :: tree
        integer(ip), intent(in), optional :: nrhs

        fmm_obj%tree => tree
        fmm_obj%pmax_mm = pmax
        fmm_obj%pmax_le = pmax
        fmm_obj%nrhs = 1
        if(present(nrhs)) fmm_obj%nrhs = nrhs

        allocate(fmm_obj%multipoles_p(ntot_sph_harm(fmm_obj%pmax_mm), fmm_obj%nrhs, tree%n_particles))
        allocate(fmm_obj%multipoles(ntot_sph_harm(fmm_obj%pmax_mm), fmm_obj%nrhs, tree%n_nodes))
        allocate(fmm_obj%local_expansion(ntot_sph_harm(fmm_obj%pmax_le), fmm_obj%nrhs, tree%n_nodes))

        fmm_obj%multipoles = 0.0
        fmm_obj%local_expansion = 0.0
        call prepare_fmmm_constants(fmm_obj%pmax_mm, fmm_obj%pmax_le)
    end subroutine

    subroutine free_fmm(fmm_obj)
        use mod_tree, only: free_tree

        implicit none

        type(fmm_type), intent(inout) :: fmm_obj

        ! call free_tree(fmm_obj%tree)
        if(allocated(fmm_obj%multipoles_p)) deallocate(fmm_obj%multipoles_p)
        if(allocated(fmm_obj%multipoles)) deallocate(fmm_obj%multipoles)
        if(allocated(fmm_obj%local_expansion)) deallocate(fmm_obj%local_expansion)
        if(allocated(fmm_obj%m2m_rotcache_fwd)) deallocate(fmm_obj%m2m_rotcache_fwd)
        if(allocated(fmm_obj%m2m_rotcache_bwd)) deallocate(fmm_obj%m2m_rotcache_bwd)
        fmm_obj%m2m_rotcache_active = .false.
        fmm_obj%m2m_rotcache_built = .false.
        if(allocated(fmm_obj%m2l_rotcache_fwd)) deallocate(fmm_obj%m2l_rotcache_fwd)
        if(allocated(fmm_obj%m2l_rotcache_bwd)) deallocate(fmm_obj%m2l_rotcache_bwd)
        fmm_obj%m2l_rotcache_active = .false.
        fmm_obj%m2l_rotcache_built = .false.
        if(allocated(fmm_obj%l2l_rotcache_fwd)) deallocate(fmm_obj%l2l_rotcache_fwd)
        if(allocated(fmm_obj%l2l_rotcache_bwd)) deallocate(fmm_obj%l2l_rotcache_bwd)
        fmm_obj%l2l_rotcache_active = .false.
        fmm_obj%l2l_rotcache_built = .false.
        if(allocated(fmm_obj%p2m_rotcache_fwd)) deallocate(fmm_obj%p2m_rotcache_fwd)
        if(allocated(fmm_obj%p2m_rotcache_bwd)) deallocate(fmm_obj%p2m_rotcache_bwd)
        fmm_obj%p2m_rotcache_active = .false.
        fmm_obj%p2m_rotcache_built = .false.
        if(allocated(fmm_obj%l2p_rotcache_fwd)) deallocate(fmm_obj%l2p_rotcache_fwd)
        if(allocated(fmm_obj%l2p_rotcache_bwd)) deallocate(fmm_obj%l2p_rotcache_bwd)
        fmm_obj%l2p_rotcache_active = .false.
        fmm_obj%l2p_rotcache_built = .false.
    end subroutine

    function fmm_available_memory_bytes() result(avail)
        !! Best-effort estimate of currently available system memory, in
        !! bytes. Falls back to a conservative fixed default if the query
        !! is unavailable (e.g. non-Linux platforms without /proc/meminfo).
        implicit none
        real(rp) :: avail
        integer :: iunit, ios
        character(len=256) :: line
        integer(8) :: kb
        real(rp), parameter :: fallback_bytes = 2.0e9_rp ! 2 GB

        avail = fallback_bytes
        open(newunit=iunit, file='/proc/meminfo', status='old', action='read', iostat=ios)
        if(ios /= 0) return
        do
            read(iunit, '(a)', iostat=ios) line
            if(ios /= 0) exit
            if(line(1:13) == 'MemAvailable:') then
                read(line(14:), *, iostat=ios) kb
                if(ios == 0) avail = real(kb, rp) * 1024.0_rp
                exit
            end if
        end do
        close(iunit)
    end function fmm_available_memory_bytes

    function fmm_m2m_rotcache_bytes(fmm_obj) result(need_bytes)
        !! Bytes needed for the M2M rotation cache (see fmm_type), or 0 if
        !! there is nothing to cache (degrees 0,1 have no rotation matrix).
        use mod_harmonics, only: oxz_rot_cache_size
        implicit none
        type(fmm_type), intent(in) :: fmm_obj
        real(rp) :: need_bytes

        need_bytes = 0.0_rp
        if(fmm_obj%pmax_mm < 2) return
        need_bytes = 2.0_rp * real(oxz_rot_cache_size(fmm_obj%pmax_mm), rp) * &
            real(fmm_obj%tree%n_nodes, rp) * 8.0_rp
    end function fmm_m2m_rotcache_bytes

    subroutine fmm_m2m_rotcache_allocate(fmm_obj)
        !! Unconditionally allocate the M2M rotation cache -- no budget
        !! check, see fmm_rotcache_try_enable for the gated entry point.
        use mod_harmonics, only: oxz_rot_cache_size
        implicit none
        type(fmm_type), intent(inout) :: fmm_obj
        integer(ip) :: csize

        if(fmm_obj%m2m_rotcache_active) return
        if(fmm_obj%pmax_mm < 2) return
        csize = oxz_rot_cache_size(fmm_obj%pmax_mm)
        allocate(fmm_obj%m2m_rotcache_fwd(csize, fmm_obj%tree%n_nodes))
        allocate(fmm_obj%m2m_rotcache_bwd(csize, fmm_obj%tree%n_nodes))
        fmm_obj%m2m_rotcache_active = .true.
        fmm_obj%m2m_rotcache_built = .false.
    end subroutine fmm_m2m_rotcache_allocate

    function fmm_m2l_rotcache_bytes(fmm_obj) result(need_bytes)
        !! Bytes needed for the M2L rotation cache (see fmm_type), or 0 if
        !! there is nothing to cache.
        use mod_harmonics, only: oxz_rot_cache_size
        implicit none
        type(fmm_type), intent(in) :: fmm_obj
        real(rp) :: need_bytes
        integer(ip) :: n_far

        need_bytes = 0.0_rp
        if(fmm_obj%pmax_mm < 2 .and. fmm_obj%pmax_le < 2) return
        n_far = fmm_obj%tree%far_nl%ri(fmm_obj%tree%n_nodes+1) - 1
        need_bytes = real(oxz_rot_cache_size(fmm_obj%pmax_mm) + &
            oxz_rot_cache_size(fmm_obj%pmax_le), rp) * real(n_far, rp) * 8.0_rp
    end function fmm_m2l_rotcache_bytes

    subroutine fmm_m2l_rotcache_allocate(fmm_obj)
        !! Unconditionally allocate the M2L rotation cache -- no budget
        !! check, see fmm_rotcache_try_enable for the gated entry point.
        use mod_harmonics, only: oxz_rot_cache_size
        implicit none
        type(fmm_type), intent(inout) :: fmm_obj
        integer(ip) :: n_far

        if(fmm_obj%m2l_rotcache_active) return
        if(fmm_obj%pmax_mm < 2 .and. fmm_obj%pmax_le < 2) return
        n_far = fmm_obj%tree%far_nl%ri(fmm_obj%tree%n_nodes+1) - 1
        allocate(fmm_obj%m2l_rotcache_fwd(oxz_rot_cache_size(fmm_obj%pmax_mm), n_far))
        allocate(fmm_obj%m2l_rotcache_bwd(oxz_rot_cache_size(fmm_obj%pmax_le), n_far))
        fmm_obj%m2l_rotcache_active = .true.
        fmm_obj%m2l_rotcache_built = .false.
    end subroutine fmm_m2l_rotcache_allocate

    function fmm_l2l_rotcache_bytes(fmm_obj) result(need_bytes)
        !! Bytes needed for the L2L rotation cache (see fmm_type), or 0 if
        !! there is nothing to cache.
        use mod_harmonics, only: oxz_rot_cache_size
        implicit none
        type(fmm_type), intent(in) :: fmm_obj
        real(rp) :: need_bytes

        need_bytes = 0.0_rp
        if(fmm_obj%pmax_le < 2) return
        need_bytes = 2.0_rp * real(oxz_rot_cache_size(fmm_obj%pmax_le), rp) * &
            real(fmm_obj%tree%n_nodes, rp) * 8.0_rp
    end function fmm_l2l_rotcache_bytes

    subroutine fmm_l2l_rotcache_allocate(fmm_obj)
        !! Unconditionally allocate the L2L rotation cache -- no budget
        !! check, see fmm_rotcache_try_enable for the gated entry point.
        use mod_harmonics, only: oxz_rot_cache_size
        implicit none
        type(fmm_type), intent(inout) :: fmm_obj
        integer(ip) :: csize

        if(fmm_obj%l2l_rotcache_active) return
        if(fmm_obj%pmax_le < 2) return
        csize = oxz_rot_cache_size(fmm_obj%pmax_le)
        allocate(fmm_obj%l2l_rotcache_fwd(csize, fmm_obj%tree%n_nodes))
        allocate(fmm_obj%l2l_rotcache_bwd(csize, fmm_obj%tree%n_nodes))
        fmm_obj%l2l_rotcache_active = .true.
        fmm_obj%l2l_rotcache_built = .false.
    end subroutine fmm_l2l_rotcache_allocate

    function fmm_p2m_rotcache_bytes(fmm_obj) result(need_bytes)
        !! Bytes needed for the P2M rotation cache (see fmm_type), or 0 if
        !! there is nothing to cache.
        use mod_harmonics, only: oxz_rot_cache_size
        implicit none
        type(fmm_type), intent(in) :: fmm_obj
        real(rp) :: need_bytes

        need_bytes = 0.0_rp
        if(fmm_obj%pmax_mm < 2) return
        need_bytes = 2.0_rp * real(oxz_rot_cache_size(fmm_obj%pmax_mm), rp) * &
            real(fmm_obj%tree%n_particles, rp) * 8.0_rp
    end function fmm_p2m_rotcache_bytes

    subroutine fmm_p2m_rotcache_allocate(fmm_obj)
        !! Unconditionally allocate the P2M rotation cache -- no budget
        !! check, see fmm_rotcache_try_enable for the gated entry point.
        use mod_harmonics, only: oxz_rot_cache_size
        implicit none
        type(fmm_type), intent(inout) :: fmm_obj
        integer(ip) :: csize

        if(fmm_obj%p2m_rotcache_active) return
        if(fmm_obj%pmax_mm < 2) return
        csize = oxz_rot_cache_size(fmm_obj%pmax_mm)
        allocate(fmm_obj%p2m_rotcache_fwd(csize, fmm_obj%tree%n_particles))
        allocate(fmm_obj%p2m_rotcache_bwd(csize, fmm_obj%tree%n_particles))
        fmm_obj%p2m_rotcache_active = .true.
        fmm_obj%p2m_rotcache_built = .false.
    end subroutine fmm_p2m_rotcache_allocate

    function fmm_l2p_rotcache_bytes(fmm_obj) result(need_bytes)
        !! Bytes needed for the L2P rotation cache (see fmm_type), or 0 if
        !! there is nothing to cache.
        use mod_harmonics, only: oxz_rot_cache_size
        implicit none
        type(fmm_type), intent(in) :: fmm_obj
        real(rp) :: need_bytes

        need_bytes = 0.0_rp
        if(fmm_obj%pmax_le < 2) return
        need_bytes = 2.0_rp * real(oxz_rot_cache_size(fmm_obj%pmax_le), rp) * &
            real(fmm_obj%tree%n_particles, rp) * 8.0_rp
    end function fmm_l2p_rotcache_bytes

    subroutine fmm_l2p_rotcache_allocate(fmm_obj)
        !! Unconditionally allocate the L2P rotation cache -- no budget
        !! check, see fmm_rotcache_try_enable for the gated entry point.
        use mod_harmonics, only: oxz_rot_cache_size
        implicit none
        type(fmm_type), intent(inout) :: fmm_obj
        integer(ip) :: csize

        if(fmm_obj%l2p_rotcache_active) return
        if(fmm_obj%pmax_le < 2) return
        csize = oxz_rot_cache_size(fmm_obj%pmax_le)
        allocate(fmm_obj%l2p_rotcache_fwd(csize, fmm_obj%tree%n_particles))
        allocate(fmm_obj%l2p_rotcache_bwd(csize, fmm_obj%tree%n_particles))
        fmm_obj%l2p_rotcache_active = .true.
        fmm_obj%l2p_rotcache_built = .false.
    end subroutine fmm_l2p_rotcache_allocate

    subroutine fmm_rotcache_try_enable(fmm_obj, mode)
        !! Joint gate for every per-node-pair OXZ rotation-matrix cache
        !! (M2M, M2L, L2L, P2M, L2P). `mode` selects which terms are even
        !! considered, matching eel%fmm_cache_mode (see mod_electrostatics):
        !!   0 (default, or argument absent) = M2L+L2P only -- these two
        !!     alone capture ~97% of the achievable speedup on a real
        !!     18k-atom system (see the `hessian` branch FMM-caching work,
        !!     2026-08-15), while M2M/L2L/P2M together buy <1% more for
        !!     real memory/complexity cost. Falls back to no caching at
        !!     all (today's on-the-fly, uncached rotations) if even this
        !!     smaller pair doesn't fit in budget.
        !!   1 = caching disabled entirely, no budget check performed.
        !!   2 = force-cache all five terms; if they don't fit in budget,
        !!     this is a hard fmm_error (NOT a silent fallback) -- the
        !!     caller explicitly asked to force it.
        !! Whichever subset is considered, the byte estimate is summed
        !! JOINTLY and checked against ONE budget before allocating
        !! anything -- independent per-type checks could each pass
        !! individually and compound past 100% of available memory.
        use mod_io, only: ommp_message
        use mod_constants, only: OMMP_VERBOSE_HIGH, OMMP_STR_CHAR_MAX

        implicit none

        type(fmm_type), intent(inout) :: fmm_obj
        integer(ip), intent(in), optional :: mode
        integer(ip) :: m
        real(rp) :: need_bytes, budget_bytes
        character(len=OMMP_STR_CHAR_MAX) :: msg
        real(rp), parameter :: budget_fraction = 0.25_rp
        !! Never use more than this fraction of available memory for the
        !! caches combined -- keeps headroom for the rest of the calculation.

        m = 0
        if(present(mode)) m = mode
        if(m == 1) return

        budget_bytes = budget_fraction * fmm_available_memory_bytes()

        if(m == 2) then
            need_bytes = fmm_m2m_rotcache_bytes(fmm_obj) + fmm_m2l_rotcache_bytes(fmm_obj) + &
                fmm_l2l_rotcache_bytes(fmm_obj) + fmm_p2m_rotcache_bytes(fmm_obj) + &
                fmm_l2p_rotcache_bytes(fmm_obj)
            if(need_bytes > budget_bytes) then
                write(msg, "(a,f0.1,a,f0.1,a)") "fmm_cache=2 (force all) requested but FMM &
                    &rotation-matrix caches need ", need_bytes/1.0e6_rp, " MB, over budget (", &
                    budget_bytes/1.0e6_rp, " MB)"
                call fmm_error(msg)
            end if
            call fmm_m2m_rotcache_allocate(fmm_obj)
            call fmm_m2l_rotcache_allocate(fmm_obj)
            call fmm_l2l_rotcache_allocate(fmm_obj)
            call fmm_p2m_rotcache_allocate(fmm_obj)
            call fmm_l2p_rotcache_allocate(fmm_obj)
            write(msg, "(a,f0.1,a)") "FMM rotation-matrix caches enabled, all five terms (", &
                need_bytes/1.0e6_rp, " MB)"
            call ommp_message(msg, OMMP_VERBOSE_HIGH)
            return
        end if

        ! mode 0 (default): M2L+L2P only
        need_bytes = fmm_m2l_rotcache_bytes(fmm_obj) + fmm_l2p_rotcache_bytes(fmm_obj)

        if(need_bytes <= budget_bytes) then
            call fmm_m2l_rotcache_allocate(fmm_obj)
            call fmm_l2p_rotcache_allocate(fmm_obj)
            write(msg, "(a,f0.1,a)") "FMM rotation-matrix caches enabled, M2L+L2P (", &
                need_bytes/1.0e6_rp, " MB)"
            call ommp_message(msg, OMMP_VERBOSE_HIGH)
        else
            write(msg, "(a,f0.1,a,f0.1,a)") "FMM rotation-matrix caches (M2L+L2P) would need ", &
                need_bytes/1.0e6_rp, " MB, over budget (", budget_bytes/1.0e6_rp, &
                " MB) -- falling back to on-the-fly rotations"
            call ommp_message(msg, OMMP_VERBOSE_HIGH)
        end if
    end subroutine fmm_rotcache_try_enable

    subroutine fmm_solve(fmm_obj)
        implicit none

        type(fmm_type), intent(inout) :: fmm_obj

        call tree_m2m(fmm_obj)
        call tree_m2l(fmm_obj)
        call tree_l2l(fmm_obj)

    end subroutine

    subroutine cart_prop_at_ipart(fmm_obj, i_part, do_V, V, do_E, E, do_grdE, grdE, do_HE, HE, do_D3E, D3E)
        implicit none

        type(fmm_type), intent(inout) :: fmm_obj
        integer(ip) :: i_part
        logical, intent(in) :: do_V, do_E, do_grdE, do_HE, do_D3E
        real(rp), intent(inout) :: V(fmm_obj%nrhs), E(3,fmm_obj%nrhs), grdE(6,fmm_obj%nrhs), &
                                    HE(10,fmm_obj%nrhs), D3E(15,fmm_obj%nrhs)

        call cart_propfar_at_ipart(fmm_obj, i_part, do_V, V, do_E, E, do_grdE, grdE, do_HE, HE, do_D3E, D3E)
        call cart_propnear_at_ipart(fmm_obj, i_part, do_V, V, do_E, E, do_grdE, grdE, do_HE, HE, do_D3E, D3E)
    end subroutine

    subroutine cart_propfar_at_ipart(fmm_obj, i_part, do_V, V, do_E, E, do_grdE, grdE, do_HE, HE, do_D3E, D3E)
        use mod_constants, only: pi
        use mod_fmm_utils, only: ntot_sph_harm
        use mod_harmonics, only: fmm_l2l
        implicit none

        ! intent(inout), not (in): the L2P rotation cache below (when
        ! active) is stashed into fmm_obj%l2p_rotcache_fwd/bwd via an
        ! array-section actual argument, which needs fmm_obj itself to be
        ! definable.
        type(fmm_type), intent(inout) :: fmm_obj
        integer(ip) :: i_part
        logical, intent(in) :: do_V, do_E, do_grdE, do_HE, do_D3E
        real(rp), intent(inout) :: V(fmm_obj%nrhs), E(3,fmm_obj%nrhs), grdE(6,fmm_obj%nrhs), &
                                    HE(10,fmm_obj%nrhs), D3E(15,fmm_obj%nrhs)

        type(fmm_tree_type), pointer :: t
        integer(ip) :: i_node, nrhs
        real(rp) :: dr(3)
        real(rp), allocatable :: x2_y2(:), z2(:), x2z_y2z(:), z3(:), xz2(:), yz2(:), x3_3xy2(:), y3_3x2y(:)
        real(rp), allocatable :: u_m4(:), u_m3(:), u_m2(:), u_m1(:), u_0(:), u_p1(:), u_p2(:), u_p3(:), u_p4(:)
        real(rp), allocatable :: tmp_local(:,:)

        t => fmm_obj%tree
        nrhs = fmm_obj%nrhs

        if(do_D3E .and. fmm_obj%pmax_le < 4) &
            call fmm_error("D3E (4th order) properties require a local expansion of order >= 4.")

        i_node = t%particle_to_node(i_part)
        allocate(tmp_local(ntot_sph_harm(fmm_obj%pmax_le), nrhs))
        allocate(x2_y2(nrhs), z2(nrhs), x2z_y2z(nrhs), z3(nrhs), xz2(nrhs), yz2(nrhs), x3_3xy2(nrhs), y3_3x2y(nrhs))
        allocate(u_m4(nrhs), u_m3(nrhs), u_m2(nrhs), u_m1(nrhs), u_0(nrhs), u_p1(nrhs), u_p2(nrhs), u_p3(nrhs), u_p4(nrhs))

        dr = t%node_centroid(:,i_node) - t%particles_coords(:,i_part)
        ! Local expansion needs a further translation
        if(fmm_obj%l2p_rotcache_active) then
            ! i_part uniquely identifies this particle's L2P shift (every
            ! particle belongs to exactly one leaf node), same reasoning as
            ! tree_p2m's cache slot. Unlike the tree_* routines,
            ! l2p_rotcache_built is NOT flipped here -- this routine is
            ! called once per particle from an external loop, so it cannot
            ! itself know when every particle has been visited; the caller
            ! (field_extD2D) flips it after its own per-particle loop ends.
            call fmm_l2l(dr, 1.0_rp, 1.0_rp, fmm_obj%pmax_le, nrhs, &
                         fmm_obj%local_expansion(:,:,i_node), tmp_local, &
                         rstack_fwd=fmm_obj%l2p_rotcache_fwd(:,i_part), &
                         rstack_bwd=fmm_obj%l2p_rotcache_bwd(:,i_part), &
                         cache_built=fmm_obj%l2p_rotcache_built)
        else
            call fmm_l2l(dr, &
                         1.0_rp, 1.0_rp, &
                         fmm_obj%pmax_le, nrhs, fmm_obj%local_expansion(:,:,i_node), &
                         tmp_local)
        end if

        if(do_V) then
            V = V + sqrt(4.0*pi) * tmp_local(1,:)
        end if

        if(do_E) then
            E(3,:) = E(3,:) - sqrt(4.0/3.0*pi) * tmp_local(3,:)
            E(1,:) = E(1,:) - sqrt(4.0/3.0*pi) * tmp_local(4,:)
            E(2,:) = E(2,:) - sqrt(4.0/3.0*pi) * tmp_local(2,:)
        end if

        if(do_grdE) then
            x2_y2 = sqrt(16.0*pi/15.0) * tmp_local(9,:) * 3.0
            z2 = (sqrt(16.0*pi/5.0) * tmp_local(7,:))
            grdE(6,:) = grdE(6,:) + z2 ! zz
            grdE(1,:) = grdE(1,:) + (x2_y2 - z2) / 2.0
            grdE(3,:) = grdE(3,:) - (x2_y2 + z2) / 2.0
            grdE(2,:) = grdE(2,:) + 3.0 * sqrt(4.0*pi/15.0) * tmp_local(5,:) !xy
            grdE(4,:) = grdE(4,:) + 3.0 * sqrt(4.0*pi/15.0) * tmp_local(8,:) !xz
            grdE(5,:) = grdE(5,:) + 3.0 * sqrt(4.0*pi/15.0) * tmp_local(6,:) !yz
        end if

        if(do_HE) then
            z3 = 15.0 * 4.0 / 5.0 * sqrt(pi / 7.0) *            tmp_local(13,:)
            x2z_y2z = 15.0 * 4.0 * sqrt(pi / 105.0) *           tmp_local(15,:)
            xz2 = 15.0 * 4.0 / 5.0 * sqrt(2.0 * pi / 21.0) *    tmp_local(14,:)
            yz2 = 15.0 * 4.0 / 5.0 * sqrt(2.0 * pi / 21.0) *    tmp_local(12,:)
            x3_3xy2 = 15.0 * 4.0 * sqrt(2.0 * pi / 35.0) *      tmp_local(16,:)
            y3_3x2y = - 15.0 * 4.0 * sqrt(2.0 * pi / 35.0) *    tmp_local(10,:)
            HE(_xyz_,:) = HE(_xyz_,:) - 15.0 * sqrt(4.0*pi/105.0) * tmp_local(11,:)
            HE(_yzz_,:) = HE(_yzz_,:) - yz2
            HE(_xzz_,:) = HE(_xzz_,:) - xz2
            HE(_zzz_,:) = HE(_zzz_,:) - z3
            HE(_xxz_,:) = HE(_xxz_,:) - (x2z_y2z - z3) / 2.0
            HE(_yyz_,:) = HE(_yyz_,:) + (x2z_y2z + z3) / 2.0
            HE(_xxx_,:) = HE(_xxx_,:) - (x3_3xy2 - 3.0 * xz2) / 4.0
            HE(_yyy_,:) = HE(_yyy_,:) - (y3_3x2y - 3.0 * yz2) / 4.0
            HE(_xyy_,:) = HE(_xyy_,:) + (x3_3xy2 + xz2) / 4.0
            HE(_xxy_,:) = HE(_xxy_,:) + (y3_3x2y + yz2) / 4.0
        end if

        if(do_D3E) then
            ! l=4 block of tmp_local is tmp_local(17:25,:), m = -4..4.
            ! D3E_abcd = d^3 E_a / dr_b dr_c dr_d = -d^4 V / dr_a dr_b dr_c dr_d,
            ! derived and validated (against a symbolic reference and against
            ! the l<=3 formulas above, using the same procedure) by fitting
            ! the local-expansion -> Cartesian map for a point charge; see
            ! scratchpad/fit_d3e2.py in the session that introduced this.
            u_m4 = sqrt(140.0*pi) * tmp_local(17,:)
            u_m3 = sqrt(70.0*pi)  * tmp_local(18,:)
            u_m2 = sqrt(20.0*pi)  * tmp_local(19,:)
            u_m1 = sqrt(10.0*pi)  * tmp_local(20,:)
            u_0  = sqrt(4.0*pi)   * tmp_local(21,:)
            u_p1 = sqrt(10.0*pi)  * tmp_local(22,:)
            u_p2 = sqrt(80.0*pi)  * tmp_local(23,:)
            u_p3 = sqrt(70.0*pi)  * tmp_local(24,:)
            u_p4 = sqrt(140.0*pi) * tmp_local(25,:)

            D3E(_xxxx_,:) = D3E(_xxxx_,:) - 3.0*u_0 + u_p2 - u_p4
            D3E(_xxxy_,:) = D3E(_xxxy_,:) - u_m4 + u_m2
            D3E(_xxxz_,:) = D3E(_xxxz_,:) + 3.0*u_p1 - u_p3
            D3E(_xxyy_,:) = D3E(_xxyy_,:) - u_0 + u_p4
            D3E(_xxyz_,:) = D3E(_xxyz_,:) - u_m3 + u_m1
            D3E(_xxzz_,:) = D3E(_xxzz_,:) + 4.0*u_0 - u_p2
            D3E(_xyyy_,:) = D3E(_xyyy_,:) + u_m4 + u_m2
            D3E(_xyyz_,:) = D3E(_xyyz_,:) + u_p1 + u_p3
            D3E(_xyzz_,:) = D3E(_xyzz_,:) - 2.0*u_m2
            D3E(_xzzz_,:) = D3E(_xzzz_,:) - 4.0*u_p1
            D3E(_yyyy_,:) = D3E(_yyyy_,:) - 3.0*u_0 - u_p2 - u_p4
            D3E(_yyyz_,:) = D3E(_yyyz_,:) + u_m3 + 3.0*u_m1
            D3E(_yyzz_,:) = D3E(_yyzz_,:) + 4.0*u_0 + u_p2
            D3E(_yzzz_,:) = D3E(_yzzz_,:) - 4.0*u_m1
            D3E(_zzzz_,:) = D3E(_zzzz_,:) - 8.0*u_0
        end if
    deallocate(tmp_local)
    deallocate(x2_y2, z2, x2z_y2z, z3, xz2, yz2, x3_3xy2, y3_3x2y)
    deallocate(u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4)
    end subroutine

    subroutine cart_propnear_at_ipart(fmm_obj, i_part, do_V, V, do_E, E, do_grdE, grdE, do_HE, HE, do_D3E, D3E)
        use mod_constants, only: pi
        use mod_fmm_utils, only: ntot_sph_harm
        use mod_harmonics, only: fmm_m2l
        implicit none

        type(fmm_type), intent(in) :: fmm_obj
        integer(ip) :: i_part
        logical, intent(in) :: do_V, do_E, do_grdE, do_HE, do_D3E
        real(rp), intent(inout) :: V(fmm_obj%nrhs), E(3,fmm_obj%nrhs), grdE(6,fmm_obj%nrhs), &
                                    HE(10,fmm_obj%nrhs), D3E(15,fmm_obj%nrhs)

        real(rp), allocatable :: local_tmp(:,:), local(:,:)
        type(fmm_tree_type), pointer :: t
        integer(ip) :: i_node, j, j_node, j_particle, jj, nrhs
        real(rp) :: c_st(3)
        real(rp), allocatable :: x2_y2(:), z2(:), x2z_y2z(:), z3(:), xz2(:), yz2(:), x3_3xy2(:), y3_3x2y(:)
        real(rp), allocatable :: u_m4(:), u_m3(:), u_m2(:), u_m1(:), u_0(:), u_p1(:), u_p2(:), u_p3(:), u_p4(:)
        t => fmm_obj%tree
        nrhs = fmm_obj%nrhs

        if(do_D3E .and. fmm_obj%pmax_le < 4) &
            call fmm_error("D3E (4th order) properties require a local expansion of order >= 4.")

        i_node = t%particle_to_node(i_part)

        allocate(local_tmp(ntot_sph_harm(fmm_obj%pmax_le), nrhs))
        allocate(local(ntot_sph_harm(fmm_obj%pmax_le), nrhs))
        allocate(x2_y2(nrhs), z2(nrhs), x2z_y2z(nrhs), z3(nrhs), xz2(nrhs), yz2(nrhs), x3_3xy2(nrhs), y3_3x2y(nrhs))
        allocate(u_m4(nrhs), u_m3(nrhs), u_m2(nrhs), u_m1(nrhs), u_0(nrhs), u_p1(nrhs), u_p2(nrhs), u_p3(nrhs), u_p4(nrhs))

        local = 0.0

        do j=t%near_nl%ri(i_node), t%near_nl%ri(i_node+1)-1
            j_node = t%near_nl%ci(j)
            do jj=t%particle_list%ri(j_node), t%particle_list%ri(j_node+1)-1
                j_particle = t%particle_list%ci(jj)
                if(i_part == j_particle) cycle
                c_st = t%particles_coords(:,j_particle) - t%particles_coords(:,i_part)

                call fmm_m2l(c_st, &
                                fmm_obj%pmax_mm, &
                                fmm_obj%pmax_le, nrhs, &
                                fmm_obj%multipoles_p(:,:,j_particle), &
                                local_tmp)
                    local = local + local_tmp
            end do
        end do

        deallocate(local_tmp)

        if(do_V) then
            V = V + sqrt(4.0*pi) * local(1,:)
        end if

        if(do_E) then
            E(3,:) = E(3,:) - sqrt(4.0/3.0*pi) * local(3,:)
            E(1,:) = E(1,:) - sqrt(4.0/3.0*pi) * local(4,:)
            E(2,:) = E(2,:) - sqrt(4.0/3.0*pi) * local(2,:)
        end if

        if(do_grdE) then
            x2_y2 = sqrt(16.0*pi/15.0) * local(9,:) * 3.0
            z2 = sqrt(16.0*pi/5.0) * local(7,:)
            grdE(6,:) = grdE(6,:) + z2 ! zz
            grdE(1,:) = grdE(1,:) + (x2_y2 - z2) / 2.0
            grdE(3,:) = grdE(3,:) - (x2_y2 + z2) / 2.0
            grdE(2,:) = grdE(2,:) + 3.0 * sqrt(4.0*pi/15.0) * local(5,:) !xy
            grdE(4,:) = grdE(4,:) + 3.0 * sqrt(4.0*pi/15.0) * local(8,:) !xz
            grdE(5,:) = grdE(5,:) + 3.0 * sqrt(4.0*pi/15.0) * local(6,:) !yz
        end if

        if(do_HE) then
            z3 = 15.0 * 4.0 / 5.0 * sqrt(pi / 7.0) * local(13,:)
            x2z_y2z = 15.0 * 4.0 * sqrt(pi / 105.0) * local(15,:)
            xz2 = 15.0 * 4.0 / 5.0 * sqrt(2.0 * pi / 21.0) * local(14,:)
            yz2 = 15.0 * 4.0 / 5.0 * sqrt(2.0 * pi / 21.0) * local(12,:)
            x3_3xy2 = 15.0 * 4.0 * sqrt(2.0 * pi / 35.0) * local(16,:)
            y3_3x2y = - 15.0 * 4.0 * sqrt(2.0 * pi / 35.0) * local(10,:)
            HE(_xyz_,:) = HE(_xyz_,:) - 15.0 * sqrt(4.0*pi/105.0) * local(11,:)
            HE(_yzz_,:) = HE(_yzz_,:) - yz2
            HE(_xzz_,:) = HE(_xzz_,:) - xz2
            HE(_zzz_,:) = HE(_zzz_,:) - z3
            HE(_xxz_,:) = HE(_xxz_,:) - (x2z_y2z - z3) / 2.0
            HE(_yyz_,:) = HE(_yyz_,:) + (x2z_y2z + z3) / 2.0
            HE(_xxx_,:) = HE(_xxx_,:) - (x3_3xy2 - 3.0 * xz2) / 4.0
            HE(_yyy_,:) = HE(_yyy_,:) - (y3_3x2y - 3.0 * yz2) / 4.0
            HE(_xyy_,:) = HE(_xyy_,:) + (x3_3xy2 + xz2) / 4.0
            HE(_xxy_,:) = HE(_xxy_,:) + (y3_3x2y + yz2) / 4.0
        end if

        if(do_D3E) then
            u_m4 = sqrt(140.0*pi) * local(17,:)
            u_m3 = sqrt(70.0*pi)  * local(18,:)
            u_m2 = sqrt(20.0*pi)  * local(19,:)
            u_m1 = sqrt(10.0*pi)  * local(20,:)
            u_0  = sqrt(4.0*pi)   * local(21,:)
            u_p1 = sqrt(10.0*pi)  * local(22,:)
            u_p2 = sqrt(80.0*pi)  * local(23,:)
            u_p3 = sqrt(70.0*pi)  * local(24,:)
            u_p4 = sqrt(140.0*pi) * local(25,:)

            D3E(_xxxx_,:) = D3E(_xxxx_,:) - 3.0*u_0 + u_p2 - u_p4
            D3E(_xxxy_,:) = D3E(_xxxy_,:) - u_m4 + u_m2
            D3E(_xxxz_,:) = D3E(_xxxz_,:) + 3.0*u_p1 - u_p3
            D3E(_xxyy_,:) = D3E(_xxyy_,:) - u_0 + u_p4
            D3E(_xxyz_,:) = D3E(_xxyz_,:) - u_m3 + u_m1
            D3E(_xxzz_,:) = D3E(_xxzz_,:) + 4.0*u_0 - u_p2
            D3E(_xyyy_,:) = D3E(_xyyy_,:) + u_m4 + u_m2
            D3E(_xyyz_,:) = D3E(_xyyz_,:) + u_p1 + u_p3
            D3E(_xyzz_,:) = D3E(_xyzz_,:) - 2.0*u_m2
            D3E(_xzzz_,:) = D3E(_xzzz_,:) - 4.0*u_p1
            D3E(_yyyy_,:) = D3E(_yyyy_,:) - 3.0*u_0 - u_p2 - u_p4
            D3E(_yyyz_,:) = D3E(_yyyz_,:) + u_m3 + 3.0*u_m1
            D3E(_yyzz_,:) = D3E(_yyzz_,:) + 4.0*u_0 + u_p2
            D3E(_yzzz_,:) = D3E(_yzzz_,:) - 4.0*u_m1
            D3E(_zzzz_,:) = D3E(_zzzz_,:) - 8.0*u_0
        end if
        deallocate(local)
        deallocate(x2_y2, z2, x2z_y2z, z3, xz2, yz2, x3_3xy2, y3_3x2y)
        deallocate(u_m4, u_m3, u_m2, u_m1, u_0, u_p1, u_p2, u_p3, u_p4)

    end subroutine
    
    subroutine tree_p2m(fmm_obj, particle_multipoles, pmax_particles)
        use mod_fmm_utils, only: ntot_sph_harm
        use mod_harmonics, only: fmm_m2m

        implicit none

        type(fmm_type), intent(inout) :: fmm_obj
        real(rp), intent(in) :: particle_multipoles(:,:,:)
        !! Shape (ncoef_particles, nrhs, n_particles)
        integer(ip), intent(in) :: pmax_particles

        type(fmm_tree_type), pointer :: t
        integer(ip) :: j, i_node, j_particle, n_expansion, p_expansion, nrhs
        real(rp), allocatable :: expansion_s(:,:), expansion_t(:,:)
        real(rp) :: c_st(3)
        logical :: use_cache, build_now

        t => fmm_obj%tree
        nrhs = fmm_obj%nrhs

        ! Sources are particles and target are nodes, in princile they could have different
        ! expansion levels (in general we should expect pmax_mm >> pmax_particles, but
        ! people are strange).

        if(fmm_obj%pmax_mm < pmax_particles) then
            call fmm_error("Multipoles expansion should be at least as high as particles expansion")
        end if
        if(size(particle_multipoles,2) /= nrhs) &
            call fmm_error("tree_p2m: particle_multipoles' nrhs dimension does not match fmm_obj%nrhs.")

        fmm_obj%multipoles_p(:,:,:) = 0.0
        fmm_obj%multipoles_p(1_ip:ntot_sph_harm(pmax_particles),:,:) = particle_multipoles(:,:,:)

        p_expansion = max(fmm_obj%pmax_mm, pmax_particles)
        n_expansion = ntot_sph_harm(p_expansion)
        allocate(expansion_s(n_expansion,nrhs), expansion_t(n_expansion,nrhs))

        ! The P2M rotation cache is sized for degree pmax_mm (see
        ! fmm_p2m_rotcache_allocate); only usable when that matches the
        ! degree this call actually translates at (always true in current
        ! usage, where pmax_particles=2 << pmax_mm, but guard it rather
        ! than assume).
        use_cache = fmm_obj%p2m_rotcache_active .and. (p_expansion == fmm_obj%pmax_mm)
        build_now = use_cache .and. .not. fmm_obj%p2m_rotcache_built

        !$omp parallel do default(shared) private(i_node, expansion_s, expansion_t, j, j_particle, c_st) schedule(dynamic)
        do i_node=1, t%n_nodes
            ! For each node
            fmm_obj%multipoles(:,:,i_node) = 0.0
            expansion_s = 0.0 ! Needed if pmax_mm > particles_pmax
            do j=t%particle_list%ri(i_node), t%particle_list%ri(i_node+1)-1
                j_particle = t%particle_list%ci(j)
                expansion_s(1:ntot_sph_harm(pmax_particles),:) = particle_multipoles(:,:, j_particle)
                c_st = t%particles_coords(:,j_particle) - t%node_centroid(:,i_node)
                ! all nrhs sources translated together in ONE call, sharing
                ! the per-node-pair rotation/translation setup
                if(use_cache) then
                    ! j_particle uniquely identifies this P2M shift (every
                    ! particle belongs to exactly one leaf node), same
                    ! reasoning as tree_m2m's cache slot
                    call fmm_m2m(c_st, p_expansion, nrhs, expansion_s, expansion_t, &
                                 rstack_fwd=fmm_obj%p2m_rotcache_fwd(:,j_particle), &
                                 rstack_bwd=fmm_obj%p2m_rotcache_bwd(:,j_particle), &
                                 cache_built=(.not. build_now))
                else
                    call fmm_m2m(c_st, p_expansion, nrhs, expansion_s, expansion_t)
                end if
                fmm_obj%multipoles(:,:,i_node) = fmm_obj%multipoles(:,:,i_node) + expansion_t
            end do
        end do

        if(build_now) fmm_obj%p2m_rotcache_built = .true.
    end subroutine

    subroutine tree_m2m(fmm_obj)
        use mod_fmm_utils, only: ntot_sph_harm
        use mod_harmonics, only: fmm_m2m

        implicit none

        type(fmm_type), intent(inout) :: fmm_obj

        type(fmm_tree_type), pointer :: t
        integer(ip) :: l, i, j, i_node, j_node, nrhs
        real(rp), allocatable :: expansion_s(:,:), expansion_t(:,:)
        real(rp) :: c_st(3)
        logical :: use_cache, build_now

        t => fmm_obj%tree
        nrhs = fmm_obj%nrhs
        allocate(expansion_s(ntot_sph_harm(fmm_obj%pmax_mm),nrhs), &
                 expansion_t(ntot_sph_harm(fmm_obj%pmax_mm),nrhs))

        use_cache = fmm_obj%m2m_rotcache_active
        build_now = use_cache .and. .not. fmm_obj%m2m_rotcache_built

        do l=t%breadth, 1, -1
            ! For each level, leaves to root
            !$omp parallel do default(shared) private(i, j, i_node, expansion_s, j_node, c_st, expansion_t)
            do i=t%level_list%ri(l), t%level_list%ri(l+1)-1
                ! For each node in the level
                i_node = t%level_list%ci(i)

                if(t%children(1,i_node) /= 0) then
                    ! P2M should already have populated leaves nodes
                    ! This node is not a leaf
                    fmm_obj%multipoles(:,:,i_node) = 0.0 ! Initialization
                    expansion_s = 0.0 ! Needed if pmax_mm < particles_pmax (unreasonable)
                    do j=1, t%tree_degree
                        j_node = t%children(j,i_node)
                        if(j_node == 0) then
                            ! This child is not present
                            cycle
                        end if

                        expansion_s = fmm_obj%multipoles(:,:, j_node)
                        c_st = t%node_centroid(:,j_node) - t%node_centroid(:,i_node)
                        ! all nrhs sources translated together in ONE call
                        if(use_cache) then
                            ! j_node uniquely identifies this M2M edge (every
                            ! non-root node has exactly one parent), so it is
                            ! a safe, race-free cache slot across the omp
                            ! parallel loop above
                            call fmm_m2m(c_st, fmm_obj%pmax_mm, nrhs, expansion_s, expansion_t, &
                                         rstack_fwd=fmm_obj%m2m_rotcache_fwd(:,j_node), &
                                         rstack_bwd=fmm_obj%m2m_rotcache_bwd(:,j_node), &
                                         cache_built=(.not. build_now))
                        else
                            call fmm_m2m(c_st, fmm_obj%pmax_mm, nrhs, expansion_s, expansion_t)
                        end if
                        fmm_obj%multipoles(:,:,i_node) = fmm_obj%multipoles(:,:,i_node) + expansion_t
                    end do
                end if
            end do
        end do

        if(build_now) fmm_obj%m2m_rotcache_built = .true.

        deallocate(expansion_s, expansion_t)
    end subroutine

    subroutine tree_m2l(fmm_obj)
        use mod_fmm_utils, only: ntot_sph_harm
        use mod_harmonics, only: fmm_m2l

        implicit none

        type(fmm_type), intent(inout) :: fmm_obj

        type(fmm_tree_type), pointer :: t
        real(rp) :: c_st(3)
        real(rp), allocatable :: mme_s(:,:), le_t(:,:)
        integer(ip) :: i_node, j_node, j, nrhs
        logical :: use_cache, build_now

        t => fmm_obj%tree
        nrhs = fmm_obj%nrhs
        allocate(mme_s(ntot_sph_harm(fmm_obj%pmax_mm),nrhs))
        allocate(le_t(ntot_sph_harm(fmm_obj%pmax_le),nrhs))

        use_cache = fmm_obj%m2l_rotcache_active
        build_now = use_cache .and. .not. fmm_obj%m2l_rotcache_built

        !$omp parallel do default(shared) &
        !$omp private(i_node, j, j_node, c_st, mme_s, le_t)
        do i_node = 1, t%n_nodes
            fmm_obj%local_expansion(:,:,i_node) = 0.0
            do j=t%far_nl%ri(i_node), t%far_nl%ri(i_node+1)-1
                j_node = t%far_nl%ci(j)

                c_st = t%node_centroid(:,j_node) - t%node_centroid(:,i_node)
                mme_s = fmm_obj%multipoles(:,:,j_node)

                ! all nrhs sources translated together in ONE call
                if(use_cache) then
                    ! j (this far-pair's position in the far_nl CSR arrays)
                    ! is a safe, race-free cache slot: the outer i_node loop
                    ! is parallelized, but each i_node's j-range is disjoint
                    ! from every other's (CSR row structure), so no two
                    ! threads ever touch the same slot.
                    call fmm_m2l(c_st, fmm_obj%pmax_mm, fmm_obj%pmax_le, nrhs, mme_s, le_t, &
                                 rstack_fwd=fmm_obj%m2l_rotcache_fwd(:,j), &
                                 rstack_bwd=fmm_obj%m2l_rotcache_bwd(:,j), &
                                 cache_built=(.not. build_now))
                else
                    call fmm_m2l(c_st, fmm_obj%pmax_mm, fmm_obj%pmax_le, nrhs, mme_s, le_t)
                end if
                fmm_obj%local_expansion(:,:,i_node) = fmm_obj%local_expansion(:,:,i_node) + le_t
            end do
        end do

        if(build_now) fmm_obj%m2l_rotcache_built = .true.

        deallocate(mme_s, le_t)
    end subroutine

    subroutine tree_l2l(fmm_obj)
        use mod_fmm_utils, only: ntot_sph_harm
        use mod_harmonics, only: fmm_l2l

        implicit none

        type(fmm_type), intent(inout) :: fmm_obj

        type(fmm_tree_type), pointer :: t
        integer(ip) :: l, i, j, i_node, j_node, nrhs
        real(rp), allocatable :: le_s(:,:), le_t(:,:)
        real(rp) :: c_st(3)
        logical :: use_cache, build_now

        t => fmm_obj%tree
        nrhs = fmm_obj%nrhs

        allocate(le_s(ntot_sph_harm(fmm_obj%pmax_le),nrhs), &
                 le_t(ntot_sph_harm(fmm_obj%pmax_le),nrhs))

        use_cache = fmm_obj%l2l_rotcache_active
        build_now = use_cache .and. .not. fmm_obj%l2l_rotcache_built

        do l=1, t%breadth-1
            ! For each level, root to leaves
            !$omp parallel do default(shared) private(i_node, le_s, j, i, j_node, c_st, le_t)
            do i=t%level_list%ri(l), t%level_list%ri(l+1)-1
                ! For each node in the level

                ! Propagate local expansion on each node on its children
                i_node = t%level_list%ci(i)
                le_s = fmm_obj%local_expansion(:,:, i_node)

                ! If node is not a leaf
                do j=1, t%tree_degree
                    j_node = t%children(j,i_node)
                    if(j_node == 0) then
                        ! This child is not present
                        cycle
                    end if

                    c_st = t%node_centroid(:,i_node) - t%node_centroid(:,j_node)
                    ! all nrhs sources translated together in ONE call
                    if(use_cache) then
                        ! j_node uniquely identifies this L2L edge (every
                        ! non-root node has exactly one parent), same
                        ! reasoning as tree_m2m's cache slot
                        call fmm_l2l(c_st, 1.0_rp, 1.0_rp, fmm_obj%pmax_le, nrhs, le_s, le_t, &
                                     rstack_fwd=fmm_obj%l2l_rotcache_fwd(:,j_node), &
                                     rstack_bwd=fmm_obj%l2l_rotcache_bwd(:,j_node), &
                                     cache_built=(.not. build_now))
                    else
                        call fmm_l2l(c_st, 1.0_rp, 1.0_rp, fmm_obj%pmax_le, nrhs, le_s, le_t)
                    end if
                    fmm_obj%local_expansion(:,:,j_node) = fmm_obj%local_expansion(:,:,j_node) + le_t
                end do
            end do
        end do

        if(build_now) fmm_obj%l2l_rotcache_built = .true.

        deallocate(le_s, le_t)
    end subroutine

    
end module
