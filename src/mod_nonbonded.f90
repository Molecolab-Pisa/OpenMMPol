module mod_nonbonded

    use mod_memory, only: rp, ip
    use mod_adjacency_mat, only: yale_sparse

    implicit none
    private

    real(rp), allocatable, dimension(:) :: vdw_r, vdw_e, vdw_f
    
    type(yale_sparse) :: vdw_pair
    real(rp), allocatable :: vdw_pair_r(:), vdw_pair_e(:)
    integer(ip), parameter :: pair_allocation_chunk = 20
    
    real(rp), dimension(4) :: vdw_screening

    logical :: use_nonbonded = .false.
    
    public :: vdw_init, vdw_terminate, vdw_potential, vdw_set_pair
    public :: vdw_r, vdw_e, vdw_f, vdw_screening, vdw_pair, vdw_pair_r, &
              vdw_pair_e, use_nonbonded

    contains

    subroutine vdw_init(vdw_type, radius_rule, radius_size, radius_type, epsrule)
        !! Initialize the non-bonded module allocating the parameters vectors
        
        use mod_memory, only: mallocate
        use mod_mmpol, only: mm_atoms, fatal_error

        implicit none

        character(len=*) :: vdw_type, radius_rule, radius_size, radius_type, &
                            epsrule

        select case(trim(vdw_type))
            case("lennard-jones")
                call fatal_error("VdW type lennard-jones is not implemented")
            case("buckingham")
                call fatal_error("VdW type buckingham is not implemented")
            case("buffered-14-7")
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
                call fatal_error("radiusrule arithmetic is not implemented")
            case("geometric")
                call fatal_error("radiusrule geometric is not implemented")
            case("cubic-mean")
                continue
            case default
                call fatal_error("radiusrule specified is not understood")
        end select
        
        select case(trim(radius_size))
            case("radius")
                call fatal_error("radiussize radius is not implemented")
            case("diameter")
                continue
            case default
                call fatal_error("radiussize specified is not understood")
        end select
        
        select case(trim(radius_type))
            case("sigma")
                call fatal_error("radiustype sigma is not implemented")
            case("r-min")
                continue
            case default
                call fatal_error("radiustype specified is not understood")
        end select
        
        select case(trim(epsrule))
            case("geometric")
                call fatal_error("epsilonrule geometric is not implemented")
            case("arithmetic")
                call fatal_error("epsilonrule arithmetic is not implemented")
            case("harmonic")
                call fatal_error("epsilonrule harmonic is not implemented")
            case("w-h")
                call fatal_error("epsilonrule w-h is not implemented")
            case("hhg")
                continue
            case default
                call fatal_error("epsilonrule specified is not understood")
        end select

        if(.not. use_nonbonded) then
            call mallocate('vdw_init [vdw_r]', mm_atoms, vdw_r)
            call mallocate('vdw_init [vdw_e]', mm_atoms, vdw_e)
            call mallocate('vdw_init [vdw_f]', mm_atoms, vdw_f)
            call mallocate('vdw_init [vdw_pair%ri]', mm_atoms+1, vdw_pair%ri)
            vdw_pair%ri = 1 ! The matrix is empty for now
            call mallocate('vdw_init [vdw_pair%ci]', pair_allocation_chunk, vdw_pair%ci)
            call mallocate('vdw_init [vdw_pair_r]', pair_allocation_chunk, vdw_pair_r)
            call mallocate('vdw_init [vdw_pair_e]', pair_allocation_chunk, vdw_pair_e)

            vdw_f = 1.0_rp
            use_nonbonded = .true.
        end if

    end subroutine vdw_init

    subroutine vdw_terminate()
        use mod_memory, only: mfree
        use mod_adjacency_mat, only: matfree

        implicit none
        
        if(.not. use_nonbonded) return

        call mfree('vdw_terminate [vdw_r]', vdw_r)
        call mfree('vdw_terminate [vdw_e]', vdw_e)
        call mfree('vdw_terminate [vdw_f]', vdw_f)
        call matfree(vdw_pair)
        call mfree('vdw_terminate [vdw_pair_r]', vdw_pair_r)
        call mfree('vdw_terminate [vdw_pair_e]', vdw_pair_e)

    end subroutine

    subroutine vdw_set_pair(ia, ib, r, e)
        !! Set VdW interaction parameters for a specific atom pair, those
        !! parameters overwrite the one obtained combining the mono-atomic
        !! ones. If a specific interaction is already set for this atom pair,
        !! it is overwritten with a warning print
        
        use mod_mmpol, only: fatal_error, mm_atoms
        use mod_io, only: ommp_message
        use mod_constants, only: OMMP_VERBOSE_LOW
        use mod_adjacency_mat, only: reallocate_mat
        use mod_memory, only: mallocate, mfree

        implicit none

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
        character(len=120) :: msg

        if(ia == ib) then 
            call fatal_error("VdW parameters could not be set for a self-interaction")
        end if

        ! To avoid saving the same interaction two times, we always use 
        ! min(ia, ib) as row index and max(ia, ib) as column index
        jr = min(ia, ib)
        jc = max(ia, ib)
        if(any(vdw_pair%ci(vdw_pair%ri(jr):vdw_pair%ri(jr+1)-1) == jc)) then
            ! The pair is already present in the matrix
            write(msg, "(A, I5, I5, A)") "VdW parameter for pair ", jr, jc, "will be overwritten"
            call ommp_message(msg, OMMP_VERBOSE_LOW)

            do i=vdw_pair%ri(jr), vdw_pair%ri(jr+1)-1
                if(vdw_pair%ci(i) == jc) then
                    vdw_pair_r = r
                    vdw_pair_e = e
                end if
            end do
        else
            ! The pair is not present and should be created
            ! 1. check if there is space in the vectors
            if(size(vdw_pair%ci) < vdw_pair%ri(mm_atoms+1) + 1) then
                ! 1b. if there is no space, allocate a new chunk
                oldsize = size(vdw_pair%ci)
                newsize = oldsize + pair_allocation_chunk
                call reallocate_mat(vdw_pair, newsize)
                call mallocate('vdw_set_pair [tmp]', oldsize, tmp)
                tmp = vdw_pair_r
                call mfree('vdw_set_pair [vdw_pair_r]', vdw_pair_r)
                call mallocate('vdw_set_pair [vdw_pair_r]', newsize, vdw_pair_r)
                vdw_pair_r(1:oldsize) = tmp
                tmp = vdw_pair_e
                call mfree('vdw_set_pair [vdw_pair_e]', vdw_pair_e)
                call mallocate('vdw_set_pair [vdw_pair_e]', newsize, vdw_pair_e)
                vdw_pair_e(1:oldsize) = tmp
                call mfree('vdw_set_pair [tmp]', tmp)
            end if
            ! 2. rewrite the r and e vectors
            do i=vdw_pair%ri(mm_atoms+1)-1, vdw_pair%ri(jr+1), -1
                vdw_pair%ci(i+1) = vdw_pair%ci(i)
                vdw_pair_r(i+1) = vdw_pair_r(i)
                vdw_pair_e(i+1) = vdw_pair_e(i)
            end do
            
            ! 3. rewrite the index vectors
            vdw_pair%ri(jr+1:) = vdw_pair%ri(jr+1:) + 1

            ! 4. write the new parameters
            vdw_pair_r(vdw_pair%ri(jr+1)-1) = r
            vdw_pair_e(vdw_pair%ri(jr+1)-1) = e
            vdw_pair%ci(vdw_pair%ri(jr+1)-1) = jc
        end if

    end subroutine vdw_set_pair

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

    subroutine vdw_potential(V)
        !! Compute the dispersion repulsion energy for the whole system
        !! using a double loop algorithm

        use mod_mmpol, only : mm_atoms, cmm, conn, fatal_error
        use mod_constants, only: eps_rp
        implicit none

        real(rp), intent(inout) :: V

        integer(ip) :: i, j, l, ipair, ineigh
        real(rp) :: eij, rij0, rij, ci(3), cj(3), s, vtmp

        do i=1, mm_atoms
            if(abs(vdw_f(i) - 1.0) < eps_rp) then
                ci = cmm(:,i)
            else
                ! Scale factors are used only for monovalent atoms, in that
                ! case the vdw center is displaced along the axis connecting
                ! the atom to its neighbour
                if(conn(1)%ri(i+1) - conn(1)%ri(i) /= 1) then
                    call fatal_error("Scale factors are only expected for monovalent atoms")
                end if
                ineigh = conn(1)%ci(conn(1)%ri(i))

                ci = cmm(:,ineigh) + (cmm(:,i) - cmm(:,ineigh)) * vdw_f(i)
            endif
                
            do j=i+1, mm_atoms
                ! Compute the screening factor for this pair
                s = 1.0
                do ineigh=1,4
                    ! Look if j is at distance ineigh from i
                    if(any(conn(ineigh)%ci(conn(ineigh)%ri(i): &
                                           conn(ineigh)%ri(i+1)-1) == j)) then
                       
                        s = vdw_screening(ineigh)
                        ! Exit the loop
                        exit 
                    end if
                end do
                
                if(s > eps_rp) then
                    ipair = -1
                    do l=vdw_pair%ri(i), vdw_pair%ri(i+1)-1
                        if(vdw_pair%ci(l) == j) then
                            ipair = l
                            exit
                        end if
                    end do

                    if(ipair > 0) then
                        Rij0 = vdw_pair_r(ipair)
                        eij = vdw_pair_e(ipair)
                    else 
                        Rij0 = (vdw_r(i)**3 + vdw_r(j)**3)/(vdw_r(i)**2 + vdw_r(j)**2)
                        eij = (4*vdw_e(i)*vdw_e(j))/(vdw_e(i)**0.5 + vdw_e(j)**0.5)**2
                    end if
                
                    if(abs(vdw_f(j) - 1.0) < eps_rp) then
                        cj = cmm(:,j)
                    else
                        ! Scale factors are used only for monovalent atoms, in that
                        ! case the vdw center is displaced along the axis connecting
                        ! the atom to its neighbour
                        if(conn(1)%ri(j+1) - conn(1)%ri(j) /= 1) then
                            write(*,*) "ATOM i J ", i, j
                            call fatal_error("Scale factors are only expected for monovalent atoms")
                        end if
                        ineigh = conn(1)%ci(conn(1)%ri(j))

                        cj = cmm(:,ineigh) + (cmm(:,j) - cmm(:,ineigh)) * vdw_f(j)
                    endif
                    Rij = ((ci(1)-cj(1))**2 + (ci(2)-cj(2))**2 + (ci(3)-cj(3))**2)**0.5
                    vtmp = 0.0
                    
                    call vdw_buffered_7_14(Rij, Rij0, Eij, vtmp)
                    v = v + vtmp*s
                end if
            end do
        end do
    end subroutine vdw_potential

end module mod_nonbonded
