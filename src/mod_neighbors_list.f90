#include "f_cart_components.h"
module mod_neighbor_list
    use mod_memory, only: ip, rp, lp, mallocate, mfree
    use mod_adjacency_mat, only: yale_sparse
    use mod_io, only: ommp_message, fatal_error

    implicit none
    private

    type ommp_neigh_list
        integer(ip) :: n
        !! Number of particles in the system
        real(rp) :: cutoff
        !! Cutoff distance for the neighbor list
        integer(ip) :: cellf = 1
        !! Number of subdivision of 'standard' cell. If you create smaller
        !! cells ([[cutoff]]/[[cellf]]) you can sensibly increase the performanc
        !! of the algorithm. See 10.1016/S0010-4655(98)00203-3
        real(rp) :: celld
        !! Dimension of each cell
        integer(ip) :: ncell(3)
        !! Number of cell along each dimension
        integer(ip) :: ncells
        !! Total number of cells
        real(rp) :: offset(3)
        !! Coordinates offset along each dimension (that is the minimum value of
        !! each coordinate.
        integer(ip) :: nneigh
        !! Number of neighbor cells
        integer(ip), allocatable :: neigh_offset(:)
        !! Offsets of neighbor cells
        integer(ip), allocatable :: p2c(:)
        !! Cell of each particle
        type(yale_sparse) :: c2p
        !! Particles contained in each cell, sparse matrix format
    end type ommp_neigh_list

    public :: ommp_neigh_list, nl_init, nl_terminate, nl_update, get_ith_nl

    contains

        subroutine nl_init(nl, c, cutoff, f)
            implicit none
            
            real(rp), intent(in) :: c(:,:)
            !! Coordinates in input
            real(rp), intent(in) :: cutoff
            !! Cut off distance
            integer(ip), intent(in) :: f
            !! Subdivision required for each cell

            type(ommp_neigh_list), intent(inout) :: nl
            !! Neigh list object to initialize

            if(size(c,1) /= 3) then
                call fatal_error("In nl_init, coordinates should be shaped 3xn")
            end if
            nl%n = size(c,2)
            nl%cutoff = cutoff
            nl%cellf = f
            if(nl%cellf > 2 .or. nl%cellf < 0) then
                call fatal_error("Subdivision required is not implemented")
            end if
            nl%nneigh = (nl%cellf*2+1)**3
            call mallocate('nl_init [neigh_offset]', nl%nneigh, nl%neigh_offset)

            nl%celld = nl%cutoff / nl%cellf
            call mallocate('nl_init [p2c]', nl%n, nl%p2c)
            call nl_update(nl, c)
        end subroutine

        subroutine nl_terminate(nl)
            use mod_adjacency_mat, only : matfree
            
            implicit none

            type(ommp_neigh_list), intent(inout) :: nl

            call matfree(nl%c2p)
            call mfree('nl_terminate [p2c]', nl%p2c)
            call mfree('nl_terminate [neigh_offset]', nl%neigh_offset)
        end subroutine

        subroutine nl_update(nl, c)
            use mod_io, only: ommp_message, time_push, time_pull
            use mod_adjacency_mat, only: reverse_grp_tab
            use mod_constants, only: OMMP_VERBOSE_LOW
            implicit none
            
            type(ommp_neigh_list), intent(inout) :: nl
            !! Neigh list object to initialize
            real(rp), intent(in) :: c(3,nl%n)
            !! Coordinates in input

            integer(ip) :: i, j, k, l, cc(3), ccmap(3)

            call time_push()
            call ommp_message('Updating neighbor lists', OMMP_VERBOSE_LOW)
            do i=1, 3
                !! TODO this should be improved 
                nl%offset(i) = minval(c(i,:))
                nl%ncell(i) = ceiling((maxval(c(i,:)) - nl%offset(i)) / nl%celld)
            end do
            nl%ncells = product(nl%ncell)
            
            ccmap(_x_) = nl%ncell(_y_) * nl%ncell(_z_)
            ccmap(_y_) = nl%ncell(_z_)
            ccmap(_z_) = 1

            l = 1
            do i=-nl%cellf, nl%cellf
                do j=-nl%cellf, nl%cellf
                    do k=-nl%cellf, nl%cellf
                        nl%neigh_offset(l) = i * ccmap(1) + j * ccmap(2) + k * ccmap(3)
                        l = l + 1
                    end do
                end do
            end do

            ! Each particle is assigned to a cell
            do i=1, nl%n
                do j=1, 3
                    cc(j) = floor((c(j,i)-nl%offset(j)) / nl%celld)
                end do
                nl%p2c(i) = dot_product(cc, ccmap) + 1
            end do
            
            ! Revert assignation to get neighbor list!
            ! The number of cell could be different...
            if(allocated(nl%c2p%ri)) then
              if(size(nl%c2p%ri) /= nl%ncells+1) then
                ! This automatically calls for the reallocation in 
                ! reverse_grp_tab.
                call mfree('nl_update [ri]', nl%c2p%ri)
              end if
            end if

            call reverse_grp_tab(nl%p2c, nl%c2p, nl%ncells)
            call time_pull("Neighbor list update")
        end subroutine

        subroutine get_ith_nl(nl, i, c, neigh, dist)
            !! Once that the neighbor list have been initialized and
            !! updated, this function provide a logical array for atom
            !! [[i]] with all interactions that should be computed and
            !! corresponding distances.
            implicit none

            type(ommp_neigh_list), intent(in) :: nl
            !! Neigh list object
            integer(ip), intent(in) :: i
            !! Index of atom for which the neigbor list is required
            real(rp), intent(in) :: c(3,nl%n)
            !! Coordinates in input
            logical(lp), intent(out) :: neigh(nl%n)
            !! Logical array for marking neighbors
            real(rp), intent(out) :: dist(nl%n)
            !! Array for returning distances

            integer(ip) :: icell, j, jid, jp, jjp
            real(rp) :: vdist(3), d2, thr2

            icell = nl%p2c(i)
            thr2 = nl%cutoff * nl%cutoff
            neigh = .false.
            dist = 0.0

            do j=1, nl%nneigh
                jid = icell + nl%neigh_offset(j)
                if(jid > 0 .and. jid <= nl%ncells) then
                    do jp=nl%c2p%ri(jid), nl%c2p%ri(jid+1)-1
                        jjp = nl%c2p%ci(jp)
                        vdist = c(:,i)-c(:,jjp)
                        d2 = dot_product(vdist, vdist)
                        if(d2 < thr2) then
                            dist(jjp) = sqrt(d2)
                            neigh(jjp) = .true.
                        end if
                    end do
                end if
            end do
        end subroutine
end module
