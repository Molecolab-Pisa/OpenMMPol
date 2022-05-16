!
! various utilities for mmpol.
!
subroutine print_header
  use mmpol
  use mod_io, only : iof_mmpol
  implicit none
!
  character (len=20) ffprt, lssolv, mvalg
!
  9000 format(t3,' .d88888b.                             888b     d888 888b     d888 8888888b.          888 ',/,&
              t3,'d88P" "Y88b                            8888b   d8888 8888b   d8888 888   Y88b         888 ',/,&
              t3,'888     888                            88888b.d88888 88888b.d88888 888    888         888 ',/,&
              t3,'888     888 88888b.   .d88b.  88888b.  888Y88888P888 888Y88888P888 888   d88P .d88b.  888 ',/,&
              t3,'888     888 888 "88b d8P  Y8b 888 "88b 888 Y888P 888 888 Y888P 888 8888888P" d88""88b 888 ',/,&
              t3,'888     888 888  888 88888888 888  888 888  Y8P  888 888  Y8P  888 888       888  888 888 ',/,&
              t3,'Y88b. .d88P 888 d88P Y8b.     888  888 888   "   888 888   "   888 888       Y88..88P 888 ',/,&
              t3,' "Y88888P"  88888P"   "Y8888  888  888 888       888 888       888 888        "Y88P"  888 ',/,&
              t3,'            888                                                                           ',/,&
              t3,'            888                                                                           ',/,&
              t3,'            888                                                                           ')
  9100 format(t3,'an open-source implementation of MMPol and AMOEBA embedding for polarizable QM/MM',/, &
              t5,'by Vladislav Slama, Lorenzo Cupellini, Benedetta Mennucci, ..., Filippo Lipparini',/, &
              t5,'MoLECoLab Pisa')

  1000 format(t3,'parameters:     ',/, &
              t3,'=============================================',/, &
              t5,'force field:           ',a,/, &
              t5,'solution method:       ',a,/, &
              t5,'algorithm:             ',a,/, &
              t5,'convergence threshold: ',d18.2,/, &
              t5,'# mm atoms:            ',i18,/, &
              t5,'# polarizable atoms:   ',i18,/, &
              t3,'=============================================')
!
  write(iof_mmpol,9000)
  write(iof_mmpol,9100)
  write(iof_mmpol,*)
  if (ff_type.eq.0) then
    if (ff_rules.eq.0) ffprt = '    AMBER (WangAL)'
    if (ff_rules.eq.1) ffprt = '    AMBER (WangDL)'
  else
    ffprt = '            AMOEBA'
  end if
! 
  if (solver.eq.0 .or. solver.eq.1) then
    lssolv = 'conjugate gradient'
  else if (solver.eq.2) then
    lssolv = '       jacobi/DIIS'
  else
    lssolv = '  matrix inversion'
  end if
! 
  if (matrix_vector.eq.0 .or. matrix_vector.eq.2) then
    mvalg = '        on-the-fly'
  else if (matrix_vector.eq.1) then
    mvalg = '            incore'
  end if
  write(iof_mmpol,1000) ffprt, lssolv, mvalg, convergence, mm_atoms, pol_atoms
  write(iof_mmpol,*)
!
  return
end subroutine print_header

subroutine print_matrix(trans,label,lda,ldb,n,m,matrix)
  use mmpol
  use mod_io, only : iof_mmpol
  implicit none
  logical,                         intent(in) :: trans
  character    (len=*),            intent(in) :: label
  integer(ip),                     intent(in) :: lda, ldb, n, m
  real(rp),    dimension(lda,ldb), intent(in) :: matrix
!
  integer(ip)        :: i, j, nbatch, nres, icol(5)
  character (len=24) :: iform, rform
!
  1000 format(t3,a)
  1010 format(t3,5i16)
  1020 format(t3,5f16.8)
  write(iof_mmpol,1000) label
  if (trans) then
!
    nbatch = n/5
    nres   = n - 5*nbatch
    write(iform,'("(t3,",i1,"i16)")') nres
    write(rform,'("(t3,",i1,"f16.8)")') nres
!
    do i = 1, nbatch
      icol(1) = (i-1)*5 + 1
      icol(2) = (i-1)*5 + 2
      icol(3) = (i-1)*5 + 3
      icol(4) = (i-1)*5 + 4
      icol(5) = (i-1)*5 + 5
      write(iof_mmpol,1010) icol(1:5)
      do j = 1, m
        write(iof_mmpol,1020) matrix(icol(1):icol(5),j)
      end do
      write(iof_mmpol,*)
    end do
!
    if (nres.ne.0) then
      do i = 1, nres
        icol(i) = 5*nbatch + i
      end do
      write(iof_mmpol,iform) icol(1:nres)
      do j = 1, m
        write(iof_mmpol,rform) matrix(icol(1):icol(nres),j)
      end do
    end if
!
  else
!
    nbatch = m/5
    nres   = m - 5*nbatch
    write(iform,'("(t3,",i1,"i16)")') nres
    write(rform,'("(t3,",i1,"f16.8)")') nres
!
    do i = 1, nbatch
      icol(1) = (i-1)*5 + 1
      icol(2) = (i-1)*5 + 2
      icol(3) = (i-1)*5 + 3
      icol(4) = (i-1)*5 + 4
      icol(5) = (i-1)*5 + 5
      write(iof_mmpol,1010) icol(1:5)
      do j = 1, n
        write(iof_mmpol,1020) matrix(j,icol(1):icol(5))
      end do
      write(iof_mmpol,*)
    end do
!
    if (nres.ne.0) then
      do i = 1, nres
        icol(i) = 5*nbatch + i
      end do
      write(iof_mmpol,iform) icol(1:nres)
      do j = 1, n
        write(iof_mmpol,rform) matrix(j,icol(1):icol(nres))
       end do
    end if
  end if
end subroutine print_matrix
!
subroutine print_int_vec(label,n,ibeg,iend,vec)
  use mmpol
  use mod_io, only : iof_mmpol
  implicit none
  character    (len=*),      intent(in) :: label
  integer(ip),               intent(in) :: n, ibeg, iend
  integer(ip), dimension(n), intent(in) :: vec
!
  integer(ip) :: ib, ie
  1000 format(t3,a)
  1010 format(t5,10i8)
!
  ib = ibeg
  ie = iend
  if (ib.eq.0) ib = 1
  if (ie.eq.0) ie = n
!
  write(iof_mmpol,1000) label
  write(iof_mmpol,1010) vec(ib:ie)
  return
!
end subroutine print_int_vec
!
subroutine ihsort(remove_duplicates,n,a)
  use mmpol
!
! sort an integer array, possibly removing duplicate entries.
! mostly taken from numerical recipes in fortran 77.
!
  implicit none
  logical,                   intent(in)    :: remove_duplicates
  integer(ip),               intent(inout) :: n
  integer(ip), dimension(n), intent(inout) :: a
!
  integer(ip) :: i, j, k, ind, asave
!
  k   = n/2 + 1
  ind = n
  do while (n.gt.1)
    if (k.gt.1) then
      k = k - 1
      asave  = a(k)
    else
      asave   = a(ind)
      a(ind)  = a(1)
      ind     = ind - 1
      if (ind.le.1) then
        a(1)  = asave
!
!       if required, get rid of duplicates:
!
        if (remove_duplicates) then
          k = 1
          do j = 2, n
            if (a(j-1) .ne. a(j)) then
              k = k + 1
              a(k) = a(j)
            end if
          end do
!
!         update the length of the array:
!
          n = k
        end if
!
!       all done.
!
        return
      end if
    end if
!
    i = k
    j = k + k
    do while (j.le.ind) 
      if (j.lt.ind) then
        if (a(j).lt.a(j+1)) j = j + 1
      end if
      if (asave .lt. a(j)) then
        a(i)  = a(j)
        i     = j
        j     = j + j
      else
        j     = ind + 1
      end if
    end do
  a(i)  = asave
  end do
!
  return
end subroutine ihsort
    

