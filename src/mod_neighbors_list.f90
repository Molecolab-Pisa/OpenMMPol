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

    public :: ommp_neigh_list, nl_init, nl_terminate, nl_update

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
            implicit none
            
            type(ommp_neigh_list), intent(inout) :: nl
            !! Neigh list object to initialize
            real(rp), intent(in) :: c(3,nl%n)
            !! Coordinates in input

            integer(ip) :: i, j, k, l, cc(3), ccmap(3)

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

            !write(*, *) "[MB23] DIMS", nl%ncell,"TOT", nl%ncells

            ! Each particle is assigned to a cell
            do i=1, nl%n
                do j=1, 3
                    cc(j) = floor((c(j,i)-nl%offset(j)) / nl%celld)
                end do
                nl%p2c(i) = dot_product(cc, ccmap) + 1
                !write(*, *) "[MB23] ", i, nl%p2c(i)
            end do
            
            ! Revert assignation to get neighbor list!
            nl%c2p%n = nl%ncells
            if(allocated(nl%c2p%ri)) call mfree('nl_update [ri]', nl%c2p%ri)
            call mallocate('nl_update [ri]', nl%c2p%n+1, nl%c2p%ri)
            if(.not. allocated(nl%c2p%ci)) call mallocate('nl_update [ci]', nl%n, nl%c2p%ci)
            nl%c2p%ri(1) = 1
            do i=1, nl%ncells
                nl%c2p%ri(i+1) = nl%c2p%ri(i)
                do j=1, nl%n
                    if(nl%p2c(j) == i) then
                        nl%c2p%ci(nl%c2p%ri(i+1)) = j
                        nl%c2p%ri(i+1) = nl%c2p%ri(i+1) + 1
                    end if
                end do
            end do
            call nl_test(nl, c)
        end subroutine
        
        subroutine nl_test(nl, c)
            use mod_constants, only: eps_rp
            implicit none
            
            type(ommp_neigh_list), intent(inout) :: nl
            !! Neigh list object to initialize
            real(rp), intent(in) :: c(3,nl%n)
            !! Coordinates in input

            integer(ip) :: i, j, ip, jp, iip, jjp, jid, pi, pj, cnt
            real(rp), allocatable :: dn2(:,:), dnl(:,:)

            allocate(dn2(nl%n,nl%n))
            allocate(dnl(nl%n,nl%n))

            do i=1, nl%n
                do j=i, nl%n
                    dn2(i,j) = norm2(c(:,i)-c(:,j))
                    if(dn2(i,j) > nl%cutoff) dn2(i,j) = -1.0
                    dn2(j,i) = dn2(i,j)
                end do
            end do

            dnl = -1.0
            cnt = 0
            do i=1, nl%ncells
                do j=1, nl%nneigh
                    jid = i + nl%neigh_offset(j)
                    if(jid > 0 .and. jid <= nl%ncells .and. i >= jid) then
                        do ip=nl%c2p%ri(i), nl%c2p%ri(i+1)-1
                            iip = nl%c2p%ci(ip)
                            do jp=nl%c2p%ri(jid), nl%c2p%ri(jid+1)-1
                                jjp = nl%c2p%ci(jp)
                                dnl(iip,jjp) = norm2(c(:,iip)-c(:,jjp))
                                if(dnl(iip,jjp) > nl%cutoff) dnl(iip,jjp) = -1.0
                                dnl(jjp,iip) = dnl(iip,jjp)
                                cnt = cnt +1
                            end do
                        end do
                    end if
                end do
            end do
            
            write(*, *) "CALC ", cnt, nl%n*(nl%n-1)/2, real(cnt)/real(nl%n*(nl%n-1)/2)
            do i=1, nl%n
                do j=i, nl%n
                    if(abs(dnl(i,j) - dn2(i,j)) > eps_rp) &
                        write(*, *) i, j, dnl(i,j), dn2(i,j)
                end do
            end do
            stop 1
        end subroutine
end module
