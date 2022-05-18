subroutine mmpol_process()
  use mod_mmpol
  use mod_memory, only: mallocate
  use mod_io, only: iof_mmpol
  use mod_constants, only: zero, one, three, six, thres
  use mod_constants, only: angstrom2au
  implicit none
!
! this routine processes the input for mmpol and amoeba calculations.
! it performs the following tasks:
!
! 1) fix the units: transform the coordinates into atomic units, as they are read
!                   in Angstrom from the input file.
!                   for amoeba, transform the polarizabilities in atomic units. 
!                   also for amoeba, scale the quadrupoles by one third.
!
! 2) create the polarizabilities array compressing the array read from input
!    and create the correspondence lists between mm atoms and polarizable atoms
!
! 3) build the connectivity lists
!
! 4) for amoeba, build the polarization group connectivity lists
!
! 5) for amoeba, rotate the multipoles to the lab frame.
!
  logical                 :: bad_neigh
  integer(ip)             :: i, ii, ij, j, k, l, neigh
  integer(ip)             :: nkeep, nlist
  real(rp)                :: third, tobohr, tobohr3, fa, fexp, xx(3)
  integer(ip),    allocatable :: mask(:), keep(:), list(:)

!
! 3) build the connectivity lists:
!
!
! amoeba block:
!
  if (amoeba) then
!
!   get memory to build the group connectivity lists:
!
    call mallocate('mmpol_process [mask]',mm_atoms,mask)
    call mallocate('mmpol_process [keep]',int(120,ip),keep)
    call mallocate('mmpol_process [list]',int(1000,ip),list)
!
!   clear the relevant arrays:
!
    n15  = 0
    i15  = 0
    np11 = 0
    np12 = 0
    ip12 = 0
    np13 = 0
    ip13 = 0
    np14 = 0
    ip14 = 0
!
!   n15, i15:
!
    do i = 1, mm_atoms
      do j = 1, n14(i)
        ij = i14(j,i)
        do k = 1, n12(ij)
          neigh = i12(k,ij)
!  
!         check that neigh is actually a 1-4 neighbor:
!  
          bad_neigh = .false.
          if (neigh.eq.i) bad_neigh = .true.
          do l = 1, n12(i)
            if (neigh.eq.i12(l,i)) bad_neigh = .true.
          end do
          do l = 1, n13(i)
            if (neigh.eq.i13(l,i)) bad_neigh = .true.
          end do
          do l = 1, n14(i)
            if (neigh.eq.i14(l,i)) bad_neigh = .true.
          end do
          if (.not. bad_neigh) then
            n15(i) = n15(i) + 1
            i15(n15(i),i) = neigh
          end if
        end do
      end do
!  
!     sort i15 and remove duplicates:
!  
      call ihsort(.true.,n15(i),i15(:,i))
    end do
!
!   np11:
!
    do i = 1, mm_atoms
      do j = 1, maxpgp
        if (ip11(j,i).ne.0) np11(i) = np11(i) + 1
      end do
!
!   sort ip11 and remove duplicates:
!
      call ihsort(.true.,np11(i),ip11(:,i))
    end do
!
!   np12, ip12:
!
    mask = 0
    keep = 0
    list = 0
!
    do i = 1, mm_atoms
!
!     create a mask to record what atoms belong to the same polarization group
!     as the i-th:
!
      do j = 1, np11(i)
        ij       = ip11(j,i)
        mask(ij) = i
      end do
!
!     now scan the 1-2 neighbors of the atoms in the same group and check whether 
!     they also belong to the same group.
!     if not, keep them:
!
      nkeep = 0
      nlist = 0
!
      do j = 1, np11(i)
        ij = ip11(j,i)
        do k = 1, n12(ij)
          l = i12(k,ij)
          if (mask(l).ne.i) then
            nkeep       = nkeep + 1
            keep(nkeep) = l 
          end if
        end do
      end do
!
!     now, for each kept element, add to a temporary list all the atoms belonging to the
!     same polarization group:
!
      do j = 1, nkeep
        ij = keep(j)
        do k = 1, np11(ij)
          l           = ip11(k,ij)
          nlist       = nlist + 1
          list(nlist) = l
        end do
      end do
!
!     sort the list and get rid of the duplicates:
!
      call ihsort(.true.,nlist,list)
!
!     if there are too many atoms in the list, abort the calculation:
!
      if (nlist.gt.maxpgp) then
        call fatal_error('too many atoms in 1-2 polarization group.')
      else
        np12(i)         = nlist
        ip12(1:nlist,i) = list(1:nlist)
      end if
    end do
!
!   np13, ip13:
!
    mask = 0
    list = 0
!
    do i = 1, mm_atoms
!
!     create a mask to record what atoms belong to the same polarization group
!     as the i-th or to a 1-2 polarization neighbor:
!
      do j = 1, np11(i)
        ij       = ip11(j,i)
        mask(ij) = i
      end do
      do j = 1, np12(i)
        ij       = ip12(j,i)
        mask(ij) = i
      end do
!
!     proceed as for standard connectivity lists. 
!     look at 1-2 polarization neighbors of 1-2 polarization neighbors:
!
      nlist = 0
!
      do j = 1, np12(i)
        ij = ip12(j,i)
        do k = 1, np12(ij)
          l = ip12(k,ij)
          if (mask(l).ne.i) then
            nlist       = nlist + 1
            list(nlist) = l 
          end if
        end do
      end do
!
!     sort the list and get rid of the duplicates:
!
      call ihsort(.true.,nlist,list)
!
!     if there are too many atoms in the list, abort the calculation:
!
      if (nlist.gt.maxpgp) then
        call fatal_error('too many atoms in 1-3 polarization group.')
      else
        np13(i)         = nlist
        ip13(1:nlist,i) = list(1:nlist)
      end if
    end do
!
!   np14, ip14:
!
    mask = 0
    list = 0
!
    do i = 1, mm_atoms
!
!     create a mask to record what atoms belong to the same polarization group
!     as the i-th or to a 1-2 or 1-3 polarization neighbor:
!
      do j = 1, np11(i)
        ij       = ip11(j,i)
        mask(ij) = i
      end do
      do j = 1, np12(i)
        ij       = ip12(j,i)
        mask(ij) = i
      end do
      do j = 1, np13(i)
        ij       = ip13(j,i)
        mask(ij) = i
      end do
!
!     proceed as for standard connectivity lists. 
!     look at 1-2 polarization neighbors of 1-3 polarization neighbors:
!
      nlist = 0
!
      do j = 1, np13(i)
        ij = ip13(j,i)
        do k = 1, np12(ij)
          l = ip12(k,ij)
          if (mask(l).ne.i) then
            nlist       = nlist + 1
            list(nlist) = l 
          end if
        end do
      end do
!
!     sort the list and get rid of the duplicates:
!
      call ihsort(.true.,nlist,list)
!
!     if there are too many atoms in the list, abort the calculation:
!
      if (nlist.gt.maxpgp) then
        call fatal_error('too many atoms in 1-4 polarization group.')
      else
        np14(i)         = nlist
        ip14(1:nlist,i) = list(1:nlist)
      end if
    end do
  end if
!
! if required, print all the relevant information:
!
  if (verbose >= 2) then
    call print_matrix(.true.,'coordinates:',3,mm_atoms,3,mm_atoms,cmm)
    if (amoeba) call print_matrix(.true.,'multipoles - non rotated:',ld_cart,mm_atoms,ld_cart,mm_atoms,q0)
    call print_matrix(.true.,'multipoles :',ld_cart,mm_atoms,ld_cart,mm_atoms,q)
    call print_matrix(.true.,'coordinates of polarizable atoms:',3,pol_atoms,3,pol_atoms,cpol)
    call print_matrix(.false.,'polarizabilities:',mm_atoms,1,mm_atoms,1,pol)
    call print_matrix(.false.,'thole factors:',mm_atoms,1,mm_atoms,1,thole)
    call print_int_vec('mm_polar list:',mm_atoms,0,0,mm_polar)
    call print_int_vec('polar_mm list:',pol_atoms,0,0,polar_mm)
!
!   write the connectivity information for each atom:
!
  1000 format(t3,'connectivity information for the ',i8,'-th atom:')
    do i = 1, mm_atoms
      write(iof_mmpol, 1000) i
      call print_int_vec('1-2 neighors:',n12(i),0,0,i12(:,i))
      call print_int_vec('1-3 neighors:',n13(i),0,0,i13(:,i))
      call print_int_vec('1-4 neighors:',n14(i),0,0,i14(:,i))
      if (amoeba) then 
        call print_int_vec('1-5 neighors:',n15(i),0,0,i15(:,i))
        call print_int_vec('1-1 polarization neighors:',np11(i),0,0,ip11(:,i))
        call print_int_vec('1-2 polarization neighors:',np12(i),0,0,ip12(:,i))
        call print_int_vec('1-3 polarization neighors:',np13(i),0,0,ip13(:,i))
        call print_int_vec('1-4 polarization neighors:',np14(i),0,0,ip14(:,i))
      end if
    end do
  end if
!
  return
end subroutine mmpol_process
