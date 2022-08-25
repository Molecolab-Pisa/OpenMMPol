module mod_io
    use mod_memory, only: ip, rp

    implicit none
    private

    integer, parameter :: iof_memory = 6
    integer, parameter :: iof_mmpol = 6
    integer, parameter :: iof_mmpinp = 100
    
    integer(ip), protected :: verbose = 0_ip
    !! verbosity flag, allowed range 0 (no printing at all) -- 
    !! 3 (debug printing)
    
    public :: iof_memory, iof_mmpol, iof_mmpinp
    public :: set_verbosity, ommp_message
    public :: print_header
    public :: print_matrix, print_int_vec

    contains
    
    subroutine set_verbosity(v)
        integer(ip), intent(in) :: v

        if( v < 0 ) then
            verbose = 0_ip
        else if( v > 3 ) then 
            verbose = 3_ip
        else
            verbose = v
        end if
    end subroutine set_verbosity

    subroutine ommp_message(s, level, logpre, u)
        use mod_constants, only: OMMP_VERBOSE_DEBUG, &
                                 OMMP_VERBOSE_HIGH, &
                                 OMMP_VERBOSE_LOW, &
                                 OMMP_VERBOSE_NONE

        implicit none

        character(len=*), intent(in) :: s
        character(len=*), intent(in), optional :: logpre
        integer(ip), intent(in) :: level
        integer(ip), intent(in), optional :: u

        integer(ip) :: outunit
        character(len=12) :: pre

        if(level == OMMP_VERBOSE_NONE) return

        if(present(u)) then
            outunit = u
        else
            outunit = 6
        end if

        if(present(logpre)) then
            write(pre, '(A12)') "["//logpre//"]"
        else
            select case(level)
                case(OMMP_VERBOSE_LOW)
                    write(pre, '(A12)') '[warning]'
                case(OMMP_VERBOSE_HIGH)
                    write(pre, '(A12)') '[info]'
                case(OMMP_VERBOSE_DEBUG)
                    write(pre, '(A12)') '[debug]'
            end select
        end if

        write(outunit, '(A12, " ", A)') pre, s
    end subroutine ommp_message
    

    subroutine print_header
      implicit none
      
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

      write(iof_mmpol,9000)
      write(iof_mmpol,9100)
      write(iof_mmpol,*)
      return
    end subroutine print_header

    subroutine print_matrix(trans,label,lda,ldb,n,m,matrix, ofunit)
      implicit none
      logical,                         intent(in) :: trans
      character    (len=*),            intent(in) :: label
      integer(ip),                     intent(in) :: lda, ldb, n, m
      real(rp),    dimension(lda,ldb), intent(in) :: matrix
      integer(ip), intent(in), optional :: ofunit
    !
      integer(ip)        :: i, j, nbatch, nres, icol(5), out_unit
      character (len=24) :: iform, rform

      if(present(ofunit)) then
          out_unit = ofunit
      else
          out_unit = iof_mmpol
      end if
    !
      1000 format(t3,a)
      1010 format(t3,5i16)
      1020 format(t3,5f16.8)
      write(out_unit,1000) label
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
          write(out_unit,1010) icol(1:5)
          do j = 1, m
            write(out_unit,1020) matrix(icol(1):icol(5),j)
          end do
          write(out_unit,'("")')
        end do
    !
        if (nres.ne.0) then
          do i = 1, nres
            icol(i) = 5*nbatch + i
          end do
          write(out_unit,iform) icol(1:nres)
          do j = 1, m
            write(out_unit,rform) matrix(icol(1):icol(nres),j)
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
          write(out_unit,1010) icol(1:5)
          do j = 1, n
            write(out_unit,1020) matrix(j,icol(1):icol(5))
          end do
          write(out_unit,*)
        end do
    !
        if (nres.ne.0) then
          do i = 1, nres
            icol(i) = 5*nbatch + i
          end do
          write(out_unit,iform) icol(1:nres)
          do j = 1, n
            write(out_unit,rform) matrix(j,icol(1):icol(nres))
           end do
        end if
      end if
    end subroutine print_matrix
    !
    subroutine print_int_vec(label, n, ibeg, iend, vec, dosort, ofunit)
        use mod_utils, only: sort_ivec
        
        implicit none
        
        character(len=*), intent(in) :: label
        integer(ip), intent(in) :: n, ibeg, iend
        logical, optional :: dosort
        integer(ip), intent(in), optional :: ofunit
        logical :: sort
        integer(ip), dimension(n), intent(in) :: vec
        
        integer(ip) :: ib, ie, out_unit

        integer(ip), allocatable, dimension(:) :: sorted_vec

        ib = ibeg
        ie = iend
        if(ib == 0) ib = 1
        if(ie == 0) ie = n
        if(.not. present(dosort)) then
            sort = .false.
        else
            sort =  dosort
        end if
        
        if(present(ofunit)) then
            out_unit = ofunit
        else
            out_unit = iof_mmpol
        end if
        
        write(out_unit, '(t3, a)') label
        
        if(ib > ie) then
            write(out_unit,'(t5)')
            return
        end if 

        if(sort) then
            call sort_ivec(vec(ib:ie), sorted_vec)
            write(out_unit,'(t5, 10i8)') sorted_vec
        else
            write(out_unit,'(t5, 10i8)') vec(ib:ie)
        end if

        return

    end subroutine print_int_vec

end module mod_io
