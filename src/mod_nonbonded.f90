module mod_nonbonded

    use mod_memory, only: rp, ip
    use mod_constants, only: OMMP_STR_CHAR_MAX
    use mod_adjacency_mat, only: yale_sparse
    use mod_topology, only: ommp_topology_type
    use mod_constants, only: OMMP_VDWTYPE_LJ, & 
                             OMMP_VDWTYPE_BUF714, &
                             OMMP_RADRULE_ARITHMETIC, &
                             OMMP_RADRULE_CUBIC, &
                             OMMP_RADTYPE_RMIN, & 
                             OMMP_EPSRULE_GEOMETRIC, &
                             OMMP_EPSRULE_HHG

    implicit none
    private

    integer(ip), parameter :: pair_allocation_chunk = 20

    type ommp_nonbonded_type
        !! Derived type for storing the information relative to the calculation
        !! of non-bonding interactions
        type(ommp_topology_type), pointer :: top
        !! Data structure for topology
        real(rp), allocatable, dimension(:) :: vdw_r
        !! VdW radii for the atoms of the system
        real(rp), allocatable, dimension(:) :: vdw_e
        !! Vdw energies for the atoms of fthe system
        real(rp), allocatable, dimension(:) :: vdw_f
        !! Scale factor for displacing the interaction center
    
        type(yale_sparse) :: vdw_pair
        !! If a pair is present in this sparse matrix, its VdW interaction
        !! should be calculated using the corresponding radius and energy
        !! not the standard ones, derived from the single-atom values.
        real(rp), allocatable :: vdw_pair_r(:)
        !! Radii for the VdW atom pairs
        real(rp), allocatable :: vdw_pair_e(:)
        !! VdW energies for atom pairs
    
        real(rp), dimension(4) :: vdw_screening = [0.0, 0.0, 1.0, 1.0]
        !! Screening factors for 1,2 1,3 and 1,4 neighbours.
        !! Default vaules from tinker manual.
        real(rp) :: radf = 1.0
        !! Scal factor for atomic radii/diameters (1.0 is used for diameters,
        !! 2.0 for radii)
        integer(ip) :: radrule
        !! Flag for the radius rule to be used
        integer(ip) :: radtype
        !! Flag for the radius type
        integer(ip) :: vdwtype
        !! Flag for the vdW type
        integer(ip) :: epsrule
        !! Flag for eps rule !!TODO
    end type ommp_nonbonded_type
    
    abstract interface
    subroutine vdw_term(Rij, Rij0, Eij, V)
        use mod_memory, only: rp
        real(rp), intent(in) :: Rij
        real(rp), intent(in) :: Rij0
        real(rp), intent(in) :: Eij
        real(rp), intent(inout) :: V
    end subroutine
    end interface
   
    public :: ommp_nonbonded_type
    public :: vdw_init, vdw_terminate, vdw_set_pair, vdw_potential, vdw_geomgrad
    public :: vdw_potential_inter, vdw_geomgrad_inter

    contains

    subroutine vdw_init(vdw, top, vdw_type, radius_rule, radius_size, &
                        radius_type, epsrule)
        !! Initialize the non-bonded object allocating the parameters vectors
        
        use mod_memory, only: mallocate
        use mod_io, only: fatal_error

        implicit none

        type(ommp_nonbonded_type), intent(inout) :: vdw
        type(ommp_topology_type), intent(in), target :: top
        character(len=*) :: vdw_type, radius_rule, radius_size, radius_type, &
                            epsrule

        select case(trim(vdw_type))
            case("lennard-jones")
                vdw%vdwtype = OMMP_VDWTYPE_LJ
                continue
            case("buckingham")
                call fatal_error("VdW type buckingham is not implemented")
            case("buffered-14-7")
                vdw%vdwtype = OMMP_VDWTYPE_BUF714
                continue
            case("mm3-hbond")
                call fatal_error("VdW type mm3-hbond is not implemented")
            case("gaussian")
                call fatal_error("VdW type gaussian is not implemented")
            case default
                call fatal_error("VdW type specified is not understood")
        end select
        
        select case(trim(radius_rule))
            case("arithmetic")
                vdw%radrule = OMMP_RADRULE_ARITHMETIC
                continue
            case("geometric")
                call fatal_error("radiusrule geometric is not implemented")
            case("cubic-mean")
                vdw%radrule = OMMP_RADRULE_CUBIC
                continue
            case default
                call fatal_error("radiusrule specified is not understood")
        end select
        
        select case(trim(radius_size))
            case("radius")
                vdw%radf = 2.0
                continue
            case("diameter")
                vdw%radf = 1.0
                continue
            case default
                call fatal_error("radiussize specified is not understood")
        end select
        
        select case(trim(radius_type))
            case("sigma")
                call fatal_error("radiustype sigma is not implemented")
            case("r-min")
                vdw%radtype = OMMP_RADTYPE_RMIN
                continue
            case default
                call fatal_error("radiustype specified is not understood")
        end select
        
        select case(trim(epsrule))
            case("geometric")
                vdw%epsrule = OMMP_EPSRULE_GEOMETRIC
                continue
            case("arithmetic")
                call fatal_error("epsilonrule arithmetic is not implemented")
            case("harmonic")
                call fatal_error("epsilonrule harmonic is not implemented")
            case("w-h")
                call fatal_error("epsilonrule w-h is not implemented")
            case("hhg")
                vdw%epsrule = OMMP_EPSRULE_HHG
                continue
            case default
                call fatal_error("epsilonrule specified is not understood")
        end select

        vdw%top => top
        call mallocate('vdw_init [vdw_r]', top%mm_atoms, vdw%vdw_r)
        call mallocate('vdw_init [vdw_e]', top%mm_atoms, vdw%vdw_e)
        call mallocate('vdw_init [vdw_f]', top%mm_atoms, vdw%vdw_f)
        call mallocate('vdw_init [vdw_pair%ri]', top%mm_atoms+1, &
                       vdw%vdw_pair%ri)
        vdw%vdw_pair%ri = 1 ! The matrix is empty for now
        call mallocate('vdw_init [vdw_pair%ci]', pair_allocation_chunk, &
                       vdw%vdw_pair%ci)
        call mallocate('vdw_init [vdw_pair_r]', pair_allocation_chunk, &
                       vdw%vdw_pair_r)
        call mallocate('vdw_init [vdw_pair_e]', pair_allocation_chunk, &
                       vdw%vdw_pair_e)

        vdw%vdw_f = 1.0_rp
    end subroutine vdw_init

    subroutine vdw_terminate(vdw)
        use mod_memory, only: mfree
        use mod_adjacency_mat, only: matfree

        implicit none
        
        type(ommp_nonbonded_type), intent(inout) :: vdw
        
        call mfree('vdw_terminate [vdw_r]', vdw%vdw_r)
        call mfree('vdw_terminate [vdw_e]', vdw%vdw_e)
        call mfree('vdw_terminate [vdw_f]', vdw%vdw_f)
        call matfree(vdw%vdw_pair)
        call mfree('vdw_terminate [vdw_pair_r]', vdw%vdw_pair_r)
        call mfree('vdw_terminate [vdw_pair_e]', vdw%vdw_pair_e)

    end subroutine

    subroutine vdw_set_pair(vdw, ia, ib, r, e)
        !! Set VdW interaction parameters for a specific atom pair, those
        !! parameters overwrite the one obtained combining the mono-atomic
        !! ones. If a specific interaction is already set for this atom pair,
        !! it is overwritten with a warning print
        
        use mod_io, only: ommp_message, fatal_error
        use mod_constants, only: OMMP_VERBOSE_LOW
        use mod_adjacency_mat, only: reallocate_mat
        use mod_memory, only: mallocate, mfree

        implicit none

        type(ommp_nonbonded_type), intent(inout) :: vdw
        !! Nonbonded data structure
        integer(ip), intent(in) :: ia
        !! First atom for which the interaction is defined
        integer(ip), intent(in) :: ib
        !! Second atom for which the interaction is defined
        real(rp), intent(in) :: r
        !! Equilibrium distance for the pair
        real(rp), intent(in) :: e 
        !! Depth of the potential 

        integer(ip) :: i, jc, jr, newsize, oldsize
        real(rp), allocatable :: tmp(:)
        character(len=OMMP_STR_CHAR_MAX) :: msg

        if(ia == ib) then 
            call fatal_error("VdW parameters could not be set for a self-interaction")
        end if

        ! To avoid saving the same interaction two times, we always use 
        ! min(ia, ib) as row index and max(ia, ib) as column index
        jr = min(ia, ib)
        jc = max(ia, ib)
        if(any(vdw%vdw_pair%ci(vdw%vdw_pair%ri(jr):&
                               vdw%vdw_pair%ri(jr+1)-1) == jc)) then
            ! The pair is already present in the matrix
            write(msg, "(A, I5, I5, A)") "VdW parameter for pair ", jr, jc, "will be overwritten"
            call ommp_message(msg, OMMP_VERBOSE_LOW)

            do i=vdw%vdw_pair%ri(jr), vdw%vdw_pair%ri(jr+1)-1
                if(vdw%vdw_pair%ci(i) == jc) then
                    vdw%vdw_pair_r = r
                    vdw%vdw_pair_e = e
                end if
            end do
        else
            ! The pair is not present and should be created
            ! 1. check if there is space in the vectors
            if(size(vdw%vdw_pair%ci) < vdw%vdw_pair%ri(vdw%top%mm_atoms+1) + 1) then
                ! 1b. if there is no space, allocate a new chunk
                oldsize = size(vdw%vdw_pair%ci)
                newsize = oldsize + pair_allocation_chunk
                call reallocate_mat(vdw%vdw_pair, newsize)
                call mallocate('vdw_set_pair [tmp]', oldsize, tmp)
                tmp = vdw%vdw_pair_r
                call mfree('vdw_set_pair [vdw_pair_r]', vdw%vdw_pair_r)
                call mallocate('vdw_set_pair [vdw_pair_r]', newsize, vdw%vdw_pair_r)
                vdw%vdw_pair_r(1:oldsize) = tmp
                tmp = vdw%vdw_pair_e
                call mfree('vdw_set_pair [vdw_pair_e]', vdw%vdw_pair_e)
                call mallocate('vdw_set_pair [vdw_pair_e]', newsize, vdw%vdw_pair_e)
                vdw%vdw_pair_e(1:oldsize) = tmp
                call mfree('vdw_set_pair [tmp]', tmp)
            end if
            ! 2. rewrite the r and e vectors
            do i=vdw%vdw_pair%ri(vdw%top%mm_atoms+1)-1, vdw%vdw_pair%ri(jr+1), -1
                vdw%vdw_pair%ci(i+1) = vdw%vdw_pair%ci(i)
                vdw%vdw_pair_r(i+1) = vdw%vdw_pair_r(i)
                vdw%vdw_pair_e(i+1) = vdw%vdw_pair_e(i)
            end do
            
            ! 3. rewrite the index vectors
            vdw%vdw_pair%ri(jr+1:) = vdw%vdw_pair%ri(jr+1:) + 1

            ! 4. write the new parameters
            vdw%vdw_pair_r(vdw%vdw_pair%ri(jr+1)-1) = r
            vdw%vdw_pair_e(vdw%vdw_pair%ri(jr+1)-1) = e
            vdw%vdw_pair%ci(vdw%vdw_pair%ri(jr+1)-1) = jc
        end if

    end subroutine vdw_set_pair

    subroutine vdw_lennard_jones(Rij, Rij0, Eij, V)
        implicit none

        real(rp), intent(in) :: Rij
        real(rp), intent(in) :: Rij0
        real(rp), intent(in) :: Eij
        real(rp), intent(inout) :: V

        real(rp) :: sigma_ov_r

        sigma_ov_r = Rij0 / Rij
        V = V + Eij*(sigma_ov_r**12 - 2*sigma_ov_r**6)

    end subroutine vdw_lennard_jones
    
    subroutine vdw_buffered_7_14(Rij, Rij0, Eij, V)
        !! Compute the dispersion-repulsion energy using the buffered 7-14 
        !! potential. Details can be found in ref: 10.1021/jp027815
        
        implicit none

        real(rp), intent(in) :: Rij
        real(rp), intent(in) :: Rij0
        real(rp), intent(in) :: Eij
        real(rp), intent(inout) :: V

        real(rp) :: rho
        real(rp), parameter :: delta = 0.07, gam = 0.12
        integer(ip), parameter :: n = 14, m = 7

        rho = Rij / Rij0
        V = V + Eij*((1+delta)/(rho+delta))**(n-m) * ((1+gam)/(rho**m+gam)-2)

    end subroutine vdw_buffered_7_14
    
    subroutine vdw_buffered_7_14_Rijgrad(Rij, Rij0, Eij, Rijgrad)
        !! Compute the gradient of vdw energy (using the buffered 7-14 
        !! potential ref: 10.1021/jp027815) with respect to the distance 
        !! between the two atoms (Rij).
        
        implicit none

        real(rp), intent(in) :: Rij
        real(rp), intent(in) :: Rij0
        real(rp), intent(in) :: Eij
        real(rp), intent(out) :: Rijgrad

        real(rp) :: rho, rhom
        real(rp), parameter :: delta = 0.07, gam = 0.12
        integer(ip), parameter :: n = 14, m = 7

        rho = Rij / Rij0
        rhom = rho**m
        Rijgrad = Eij * ((delta+1)/(delta+rho))**(n-m) * &
                  (rho*(gam+rhom)*(m-n)*(1-2*rhom-gam) - &
                  rhom*m*(delta+rho)*(gam+1)) / &
                  (rho*(delta+rho)*(gam+rhom)**2)
        Rijgrad = Rijgrad / Rij0

    end subroutine vdw_buffered_7_14_Rijgrad

    pure function get_Rij0(vdw, i, j) result(Rij0)
        use mod_io, only : fatal_error
        
        implicit none

        type(ommp_nonbonded_type), intent(in), target :: vdw
        !! Nonbonded data structure
        integer(ip), intent(in) :: i, j
        !! Indices of interacting atoms
        real(rp) :: Rij0

        select case(vdw%radrule) 
            case(OMMP_RADRULE_ARITHMETIC)
                Rij0 = vdw%radf*(vdw%vdw_r(i) + vdw%vdw_r(j))/2
            case(OMMP_RADRULE_CUBIC)
                Rij0 = (vdw%vdw_r(i)**3 + vdw%vdw_r(j)**3) / &
                       (vdw%vdw_r(i)**2 + vdw%vdw_r(j)**2)
            case default
                Rij0 = 0.0
                continue
                !call fatal_error("Unexpected error in get_Rij0")
        end select
    end function
    
    pure function get_eij(vdw, i, j) result(eij)
        use mod_io, only : fatal_error
        
        implicit none

        type(ommp_nonbonded_type), intent(in) :: vdw
        !! Nonbonded data structure
        integer(ip), intent(in) :: i, j
        !! Indices of interacting atoms
        real(rp) :: eij

        select case(vdw%epsrule) 
            case(OMMP_EPSRULE_GEOMETRIC)
                eij = sqrt(vdw%vdw_e(i)*vdw%vdw_e(j))
            case(OMMP_EPSRULE_HHG)
                eij = (4*vdw%vdw_e(i)*vdw%vdw_e(j)) / &
                      (vdw%vdw_e(i)**0.5 + vdw%vdw_e(j)**0.5)**2
            case default
                eij = 0.0
                continue
                !call fatal_error("Unexpected error in get_eij")
        end select
    end function

    subroutine vdw_potential(vdw, V)
        !! Compute the dispersion repulsion energy for the whole system
        !! using a double loop algorithm

        use mod_io, only : fatal_error
        use mod_constants, only: eps_rp
        implicit none

        type(ommp_nonbonded_type), intent(in), target :: vdw
        !! Nonbonded data structure
        real(rp), intent(inout) :: V
        !! Potential, result will be added

        integer(ip) :: i, j, l, ipair, ineigh
        real(rp) :: eij, rij0, rij, ci(3), cj(3), s, vtmp
        type(ommp_topology_type), pointer :: top
        procedure(vdw_term), pointer :: vdw_func

        top => vdw%top
        select case(vdw%vdwtype)
            case(OMMP_VDWTYPE_LJ)
                vdw_func => vdw_lennard_jones
            case(OMMP_VDWTYPE_BUF714)
                vdw_func => vdw_buffered_7_14
            case default
                call fatal_error("Unexpected error in vdw_potential")
        end select

        do i=1, top%mm_atoms
            if(abs(vdw%vdw_f(i) - 1.0) < eps_rp) then
                ci = top%cmm(:,i)
            else
                ! Scale factors are used only for monovalent atoms, in that
                ! case the vdw center is displaced along the axis connecting
                ! the atom to its neighbour
                if(top%conn(1)%ri(i+1) - top%conn(1)%ri(i) /= 1) then
                    call fatal_error("Scale factors are only expected for &
                                     &monovalent atoms")
                end if
                ineigh = top%conn(1)%ci(top%conn(1)%ri(i))

                ci = top%cmm(:,ineigh) + (top%cmm(:,i) - top%cmm(:,ineigh)) &
                     * vdw%vdw_f(i)
            endif
                
            do j=i+1, top%mm_atoms
                ! Compute the screening factor for this pair
                s = 1.0
                do ineigh=1,4
                    ! Look if j is at distance ineigh from i
                    if(any(top%conn(ineigh)%ci(top%conn(ineigh)%ri(i): &
                                               top%conn(ineigh)%ri(i+1)-1) == j)) then
                       
                        s = vdw%vdw_screening(ineigh)
                        ! Exit the loop
                        exit 
                    end if
                end do
                
                if(s > eps_rp) then
                    ipair = -1
                    do l=vdw%vdw_pair%ri(i), vdw%vdw_pair%ri(i+1)-1
                        if(vdw%vdw_pair%ci(l) == j) then
                            ipair = l
                            exit
                        end if
                    end do

                    if(ipair > 0) then
                        Rij0 = vdw%vdw_pair_r(ipair)
                        eij = vdw%vdw_pair_e(ipair)
                    else
                        Rij0 = get_Rij0(vdw, i, j)
                        eij = get_eij(vdw, i, j)
                    end if
                
                    if(abs(vdw%vdw_f(j) - 1.0) < eps_rp) then
                        cj = top%cmm(:,j)
                    else
                        ! Scale factors are used only for monovalent atoms, in that
                        ! case the vdw center is displaced along the axis connecting
                        ! the atom to its neighbour
                        if(top%conn(1)%ri(j+1) - top%conn(1)%ri(j) /= 1) then
                            call fatal_error("Scale factors are only expected &
                                             & for monovalent atoms")
                        end if
                        ineigh = top%conn(1)%ci(top%conn(1)%ri(j))

                        cj = top%cmm(:,ineigh) + &
                             (top%cmm(:,j) - top%cmm(:,ineigh)) * vdw%vdw_f(j)
                    endif
                    Rij = norm2(ci-cj)
                    
                    vtmp = 0.0
                    call vdw_func(Rij, Rij0, Eij, vtmp)
                    v = v + vtmp*s
                end if
            end do
        end do
    end subroutine vdw_potential
    
    subroutine vdw_geomgrad(vdw, grad)
        !! Compute the dispersion repulsion energy for the whole system
        !! using a double loop algorithm

        use mod_io, only : fatal_error
        use mod_constants, only: eps_rp
        use mod_jacobian_mat, only: Rij_jacobian
        implicit none

        type(ommp_nonbonded_type), intent(in), target :: vdw
        !! Nonbonded data structure
        real(rp), intent(inout) :: grad(3,vdw%top%mm_atoms)
        !! Gradients, result will be added

        integer(ip) :: i, j, l, ipair, ineigh, ineigh_i, ineigh_j
        real(rp) :: eij, rij0, rij, ci(3), cj(3), s, J_i(3), J_j(3), Rijg, &
                    f_i, f_j
        type(ommp_topology_type), pointer :: top

        top => vdw%top

        do i=1, top%mm_atoms
            if(abs(vdw%vdw_f(i) - 1.0) < eps_rp) then
                ci = top%cmm(:,i)
                ineigh_i = 0 ! This is needed later for force projection
            else
                ! Scale factors are used only for monovalent atoms, in that
                ! case the vdw center is displaced along the axis connecting
                ! the atom to its neighbour
                if(top%conn(1)%ri(i+1) - top%conn(1)%ri(i) /= 1) then
                    call fatal_error("Scale factors are only expected for &
                                     &monovalent atoms")
                end if
                ineigh_i = top%conn(1)%ci(top%conn(1)%ri(i))
                f_i = vdw%vdw_f(i)

                ci = top%cmm(:,ineigh_i) + (top%cmm(:,i) - &
                                            top%cmm(:,ineigh_i)) * f_i
            endif
                
            do j=i+1, top%mm_atoms
                ! Compute the screening factor for this pair
                s = 1.0
                do ineigh=1,4
                    ! Look if j is at distance ineigh from i
                    if(any(top%conn(ineigh)%ci(top%conn(ineigh)%ri(i): &
                                               top%conn(ineigh)%ri(i+1)-1) == j)) then
                       
                        s = vdw%vdw_screening(ineigh)
                        ! Exit the loop
                        exit 
                    end if
                end do
                
                if(s > eps_rp) then
                    ipair = -1
                    do l=vdw%vdw_pair%ri(i), vdw%vdw_pair%ri(i+1)-1
                        if(vdw%vdw_pair%ci(l) == j) then
                            ipair = l
                            exit
                        end if
                    end do

                    if(ipair > 0) then
                        Rij0 = vdw%vdw_pair_r(ipair)
                        eij = vdw%vdw_pair_e(ipair)
                    else 
                        Rij0 = (vdw%vdw_r(i)**3 + vdw%vdw_r(j)**3) / &
                               (vdw%vdw_r(i)**2 + vdw%vdw_r(j)**2)
                        eij = (4*vdw%vdw_e(i)*vdw%vdw_e(j)) / &
                              (vdw%vdw_e(i)**0.5 + vdw%vdw_e(j)**0.5)**2
                    end if
                
                    if(abs(vdw%vdw_f(j) - 1.0) < eps_rp) then
                        cj = top%cmm(:,j)
                        ineigh_j = 0 ! This is needed later for force projection
                    else
                        ! Scale factors are used only for monovalent atoms, in that
                        ! case the vdw center is displaced along the axis connecting
                        ! the atom to its neighbour
                        if(top%conn(1)%ri(j+1) - top%conn(1)%ri(j) /= 1) then
                            call fatal_error("Scale factors are only expected &
                                             & for monovalent atoms")
                        end if
                        ineigh_j = top%conn(1)%ci(top%conn(1)%ri(j))
                        f_j = vdw%vdw_f(j)

                        cj = top%cmm(:,ineigh_j) + &
                             (top%cmm(:,j) - top%cmm(:,ineigh_j)) * f_j
                    endif

                    call Rij_jacobian(ci, cj, Rij, J_i, J_j)
                    call vdw_buffered_7_14_Rijgrad(Rij, Rij0, Eij, Rijg)

                    if(ineigh_i == 0) then
                        grad(:,i) =  grad(:,i) + J_i * Rijg
                    else
                        ! If the center is displaced, the forces should be 
                        ! projected onto the two atoms that determine the
                        ! position of the center
                        grad(:,i) = grad(:,i) + J_i * Rijg * f_i
                        grad(:,ineigh_i) = grad(:,ineigh_i) + J_i * Rijg * (1-f_i)
                    end if

                    if(ineigh_j == 0) then
                        grad(:,j) =  grad(:,j) + J_j * Rijg
                    else
                        ! If the center is displaced, the forces should be 
                        ! projected onto the two atoms that determine the
                        ! position of the center
                        grad(:,j) = grad(:,j) + J_j * Rijg * f_j
                        grad(:,ineigh_j) = grad(:,ineigh_j) + J_j * Rijg * (1-f_j)
                    endif
                end if
            end do
        end do
    end subroutine
    
    subroutine vdw_potential_inter(vdw1, vdw2, V)
        !! Compute the dispersion repulsion energy for the whole system
        !! using a double loop algorithm

        use mod_io, only : fatal_error
        use mod_constants, only: eps_rp
        implicit none

        type(ommp_nonbonded_type), intent(in), target :: vdw1, vdw2
        !! Nonbonded data structure
        real(rp), intent(inout) :: V
        !! Potential, result will be added

        integer(ip) :: i, j, l, ipair, ineigh
        real(rp) :: eij, rij0, rij, ci(3), cj(3), s, vtmp
        type(ommp_topology_type), pointer :: top1, top2

        top1 => vdw1%top
        top2 => vdw2%top
        
        do i=1, top1%mm_atoms
            if(abs(vdw1%vdw_f(i) - 1.0) < eps_rp) then
                ci = top1%cmm(:,i)
            else
                ! Scale factors are used only for monovalent atoms, in that
                ! case the vdw center is displaced along the axis connecting
                ! the atom to its neighbour
                if(top1%conn(1)%ri(i+1) - top1%conn(1)%ri(i) /= 1) then
                    call fatal_error("Scale factors are only expected for &
                                     &monovalent atoms")
                end if
                ineigh = top1%conn(1)%ci(top1%conn(1)%ri(i))

                ci = top1%cmm(:,ineigh) + (top1%cmm(:,i) - top1%cmm(:,ineigh)) &
                     * vdw1%vdw_f(i)
            endif
                
            do j=1, top2%mm_atoms
                Rij0 = (vdw1%vdw_r(i)**3 + vdw2%vdw_r(j)**3) / &
                       (vdw1%vdw_r(i)**2 + vdw2%vdw_r(j)**2)
                eij = (4*vdw1%vdw_e(i)*vdw2%vdw_e(j)) / &
                      (vdw1%vdw_e(i)**0.5 + vdw2%vdw_e(j)**0.5)**2
            
                if(abs(vdw2%vdw_f(j) - 1.0) < eps_rp) then
                    cj = top2%cmm(:,j)
                else
                    ! Scale factors are used only for monovalent atoms, in that
                    ! case the vdw center is displaced along the axis connecting
                    ! the atom to its neighbour
                    if(top2%conn(1)%ri(j+1) - top2%conn(1)%ri(j) /= 1) then
                        call fatal_error("Scale factors are only expected &
                                         & for monovalent atoms")
                    end if
                    ineigh = top2%conn(1)%ci(top2%conn(1)%ri(j))

                    cj = top2%cmm(:,ineigh) + &
                         (top2%cmm(:,j) - top2%cmm(:,ineigh)) * vdw2%vdw_f(j)
                endif
                Rij = norm2(ci-cj)
                
                vtmp = 0.0
                call vdw_buffered_7_14(Rij, Rij0, Eij, vtmp)
                v = v + vtmp
            end do
        end do
    end subroutine vdw_potential_inter
    
    subroutine vdw_geomgrad_inter(vdw1, vdw2, grad1, grad2)
        !! Compute the dispersion repulsion energy for the whole system
        !! using a double loop algorithm

        use mod_io, only : fatal_error
        use mod_constants, only: eps_rp
        use mod_jacobian_mat, only: Rij_jacobian
        implicit none

        type(ommp_nonbonded_type), intent(in), target :: vdw1, vdw2
        !! Nonbonded data structure
        real(rp), intent(inout) :: grad1(3,vdw1%top%mm_atoms), &
                                   grad2(3,vdw1%top%mm_atoms)
        !! Potential, result will be added

        integer(ip) :: i, j, l, ipair, ineigh_i, ineigh_j
        real(rp) :: eij, rij0, rij, ci(3), cj(3), s, vtmp, Rijg, f_i, f_j, &
                    J_i(3), J_j(3)
        type(ommp_topology_type), pointer :: top1, top2

        top1 => vdw1%top
        top2 => vdw2%top
        
        do i=1, top1%mm_atoms
            if(abs(vdw1%vdw_f(i) - 1.0) < eps_rp) then
                ci = top1%cmm(:,i)
                ineigh_i = 0
            else
                ! Scale factors are used only for monovalent atoms, in that
                ! case the vdw center is displaced along the axis connecting
                ! the atom to its neighbour
                if(top1%conn(1)%ri(i+1) - top1%conn(1)%ri(i) /= 1) then
                    call fatal_error("Scale factors are only expected for &
                                     &monovalent atoms")
                end if
                ineigh_i = top1%conn(1)%ci(top1%conn(1)%ri(i))
                f_i = vdw1%vdw_f(i)
                ci = top1%cmm(:,ineigh_i) + (top1%cmm(:,i) - top1%cmm(:,ineigh_i)) &
                     * f_i
            endif
                
            do j=1, top2%mm_atoms
                Rij0 = (vdw1%vdw_r(i)**3 + vdw2%vdw_r(j)**3) / &
                       (vdw1%vdw_r(i)**2 + vdw2%vdw_r(j)**2)
                eij = (4*vdw1%vdw_e(i)*vdw2%vdw_e(j)) / &
                      (vdw1%vdw_e(i)**0.5 + vdw2%vdw_e(j)**0.5)**2
            
                if(abs(vdw2%vdw_f(j) - 1.0) < eps_rp) then
                    cj = top2%cmm(:,j)
                    ineigh_j = 0
                else
                    ! Scale factors are used only for monovalent atoms, in that
                    ! case the vdw center is displaced along the axis connecting
                    ! the atom to its neighbour
                    if(top2%conn(1)%ri(j+1) - top2%conn(1)%ri(j) /= 1) then
                        call fatal_error("Scale factors are only expected &
                                         & for monovalent atoms")
                    end if
                    ineigh_j = top2%conn(1)%ci(top2%conn(1)%ri(j))
                    f_j = vdw2%vdw_f(j)

                    cj = top2%cmm(:,ineigh_j) + &
                         (top2%cmm(:,j) - top2%cmm(:,ineigh_j)) * f_j
                endif
                
                call Rij_jacobian(ci, cj, Rij, J_i, J_j)
                call vdw_buffered_7_14_Rijgrad(Rij, Rij0, Eij, Rijg)

                if(ineigh_i == 0) then
                    grad1(:,i) =  grad1(:,i) + J_i * Rijg
                else
                    ! If the center is displaced, the forces should be 
                    ! projected onto the two atoms that determine the
                    ! position of the center
                    grad1(:,i) = grad1(:,i) + J_i * Rijg * f_i
                    grad1(:,ineigh_i) = grad1(:,ineigh_i) + J_i * Rijg * (1-f_i)
                end if

                if(ineigh_j == 0) then
                    grad2(:,j) =  grad2(:,j) + J_j * Rijg
                else
                    ! If the center is displaced, the forces should be 
                    ! projected onto the two atoms that determine the
                    ! position of the center
                    grad2(:,j) = grad2(:,j) + J_j * Rijg * f_j
                    grad2(:,ineigh_j) = grad2(:,ineigh_j) + J_j * Rijg * (1-f_j)
                endif
            end do
        end do
    end subroutine

end module mod_nonbonded
