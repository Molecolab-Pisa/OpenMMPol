#include "version.h"

module mod_io
   !! Unified Input/Output handling across the code.

   use mod_constants, only: OMMP_VERBOSE_DEBUG, &
      OMMP_VERBOSE_HIGH, &
      OMMP_VERBOSE_LOW, &
      OMMP_VERBOSE_NONE, &
      OMMP_VERBOSE_DEFAULT, &
      ip, rp

   implicit none
   private

   integer :: iof_mmpol = 6
   integer, parameter :: iof_mmpinp = 100

   integer(ip), protected :: verbose = OMMP_VERBOSE_DEFAULT
   !! verbosity flag, allowed range 0 (no printing at all) --
   !! 3 (debug printing)

   public :: iof_mmpol, iof_mmpinp
   public :: set_iof_mmpol, close_output
   public :: set_verbosity, ommp_message, fatal_error, ommp_version
   public :: print_matrix, print_int_vec
   public :: large_file_read

   interface print_matrix
      !! Interface for matrix printing function
      module procedure d1_print_matrix
      module procedure d2_print_matrix
   end interface print_matrix

contains

   subroutine set_iof_mmpol(filename)
      !! This subroutine changes the output file for mmpol to a file defined by filename.
      use mod_constants, only: OMMP_STR_CHAR_MAX

      implicit none

      character(len=*), intent(in) :: filename
      !! File name for the new output stream
      character(len=OMMP_STR_CHAR_MAX) :: msg, oldfname
      integer(ip) :: ist


      if(iof_mmpol /= 6) then
         !! A file has already been set, close it before proceed.
         inquire(unit=iof_mmpol, name=oldfname)
         write(msg, '("Switching output from ", a, " to ", a,".")') trim(oldfname), trim(filename)
         call ommp_message(msg, OMMP_VERBOSE_LOW)
         close(iof_mmpol)
      else
         write(msg, '("Switching output from stdout to ", a,".")') trim(filename)
         call ommp_message(msg, OMMP_VERBOSE_LOW)
      end if

      open(unit=110, &
         file=filename, &
         form='formatted', &
         iostat=ist, &
         action='write')

      if(ist /= 0) then
         call fatal_error('Error while opening output input file')
      end if

      iof_mmpol = 110
   end subroutine

   subroutine close_output()
      !! This subroutine changes the output file for mmpol to a file defined by filename.
      use mod_constants, only: OMMP_STR_CHAR_MAX

      implicit none

      character(len=OMMP_STR_CHAR_MAX) :: msg, oldfname


      if(iof_mmpol /= 6) then
         !! A file has already been set, close it before proceed.
         inquire(unit=iof_mmpol, name=oldfname)
         write(msg, '("Closing output file ", a,".")') trim(oldfname)
         call ommp_message(msg, OMMP_VERBOSE_LOW)
         close(iof_mmpol)
      end if

   end subroutine

   subroutine set_verbosity(v)
      !! Set the verbosity level for the output, this is a library-level
      !! function, that changes the behaviour of several I/O functions. It
      !! also enforces min/max verbosity levels (currently no output is 0,
      !! while debug output is 3).

      integer(ip), intent(in) :: v
      !! Requested level of verbosity

      if( v < 0 ) then
         verbose = 0_ip
      else if( v > 3 ) then
         verbose = 3_ip
      else
         verbose = v
      end if
   end subroutine set_verbosity

   subroutine ommp_message(s, level, logpre, u)
      !! Output a message according to the verbosity level.
      !! @note All the output messages should pass for this function, if it
      !! is necessary, first use write to format the message on a string and
      !! then pass the string to [[ommp_message]]

      implicit none

      character(len=*), intent(in) :: s
      !! Message to be printed
      character(len=*), intent(in), optional :: logpre
      !! String that explains the type of message, if missing a default type
      !! is assigned based on requested verbosity level
      integer(ip), intent(in) :: level
      !! Requested verbosity level
      integer(ip), intent(in), optional :: u
      !! Output unit for the message, if missing, [[iof_mmpol]] is used.

      integer(ip) :: outunit
      character(len=12) :: pre

      if(level > verbose) return

      if(present(u)) then
         outunit = u
      else
         outunit = iof_mmpol
      end if

      if(present(logpre)) then
         write(pre, '(A12)') "["//trim(logpre)//"]"
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

      write(outunit, '(A6, A12, " ", A)') '[OMMP]', pre, trim(s)
   end subroutine ommp_message

   subroutine fatal_error(message)
      !! Prints a message and exit from the program. This
      !! function should be used in all the conditions
      !! where the program cannot proceed.

      implicit none

      character (len=*), intent(in) :: message
      !! Message to print before the program termination
      call ommp_message(message, OMMP_VERBOSE_LOW, 'stop')
      call ommp_message("Unrecoverable error in openMMPol &
      &library. Exiting.", OMMP_VERBOSE_LOW, &
         'stop')
      !! Close output file
      call close_output()

      !TODO call mmpol_terminate()

      !stop 1
      call exit(1)
   end subroutine fatal_error

   subroutine ommp_version(v)
      integer(ip), intent(in) :: v

      call ommp_message("OpenMMPol version: "//_OMMP_VERSION, &
         v, "version")

   end subroutine

   subroutine d1_print_matrix(trans, label, matrix, ofunit)
      !! Output a 1D-matrix of real in a well formatted way
      implicit none

      logical, intent(in) :: trans
      !! It the matrix is transposed or not
      character(len=*), intent(in) :: label
      !! Label to be printed before the matrix
      real(rp), dimension(:), intent(in) :: matrix
      !! Matrix to be printed
      integer(ip), intent(in), optional :: ofunit
      !! Unit where the matrix should be printed, if not present [[iof_mmpol]]
      !! is used.

      call d2_print_matrix(trans, label, &
         reshape(matrix, [size(matrix), 1]), ofunit)

   end subroutine d1_print_matrix

   subroutine d2_print_matrix(trans, label, matrix, ofunit)
      !! Output a 2D-matrix of real in a well formatted way
      implicit none

      logical, intent(in) :: trans
      !! It the matrix is transposed or not
      character(len=*), intent(in) :: label
      !! Label to be printed before the matrix
      real(rp), dimension(:, :), intent(in) :: matrix
      !! Matrix to be printed
      integer(ip), intent(in), optional :: ofunit
      !! Unit where the matrix should be printed, if not present [[iof_mmpol]]
      !! is used.

      integer(ip) :: i, j, nbatch, nres, out_unit, nrow, ncol
      integer(ip), parameter :: colperrow = 5
      integer(ip), dimension(colperrow) :: icol
      character (len=24) :: iform, rform

      if(present(ofunit)) then
         out_unit = ofunit
      else
         out_unit = iof_mmpol
      end if

      nrow = size(matrix, 1)
      ncol = size(matrix, 2)

      icol = [(i, i=1, colperrow)]

1000  format(t3,a)
1010  format(t3,5i16)
1020  format(t3,5f16.8)

      write(out_unit,1000) label

      if(trans) then
         nbatch = nrow/colperrow
         nres   = mod(nrow, colperrow)

         write(iform,'("(t3,",i1,"i16)")') nres
         write(rform,'("(t3,",i1,"f16.8)")') nres

         do i = 1, nbatch
            write(out_unit, 1010) icol(1:colperrow)

            do j = 1, ncol
               write(out_unit, 1020) matrix(icol(1):icol(colperrow), j)
            end do

            write(out_unit,'("")')

            icol = icol + colperrow
         end do

         if (nres > 0) then
            write(out_unit,iform) icol(1:nres)

            do j = 1, ncol
               write(out_unit,rform) matrix(icol(1):icol(nres),j)
            end do
         end if
      else
         nbatch = ncol/colperrow
         nres   = mod(ncol,colperrow)
         write(iform,'("(t3,",i1,"i16)")') nres
         write(rform,'("(t3,",i1,"f16.8)")') nres

         do i = 1, nbatch
            write(out_unit,1010) icol(1:colperrow)

            do j = 1, nrow
               write(out_unit,1020) matrix(j,icol(1):icol(colperrow))
            end do

            write(out_unit,*)
            icol = icol + colperrow
         end do

         if(nres.ne.0) then
            write(out_unit,iform) icol(1:nres)

            do j = 1, nrow
               write(out_unit,rform) matrix(j,icol(1):icol(nres))
            end do
         end if
      end if
   end subroutine d2_print_matrix

   subroutine print_int_vec(label, vec, ofunit, ibeg, iend)
      !! Print an array of integers in a well formatted way.

      implicit none

      character(len=*), intent(in) :: label
      !! Label to print before the array
      integer(ip), dimension(:), intent(in) :: vec
      !! Integer vector to be printed
      integer(ip), intent(in), optional :: ofunit
      !! If present specify the unit where the array should be printed,
      !! otherwise [[iof_mmpol]] is used.
      integer(ip), intent(in), optional :: ibeg
      !! Index of the first element to be printed, if 0 or absent the first
      !! element is used
      integer(ip), intent(in), optional :: iend
      !! Index of the last element to be printed, if 0 or absent the last

      integer(ip) :: ib, ie, out_unit

      if(present(ibeg)) then
         ib = ibeg
      else
         ib = 0
      end if
      if(present(ibeg)) then
         ie = iend
      else
         ie = 0
      end if
      if(ib == 0) ib = 1
      if(ie == 0) ie = size(vec)

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

      write(out_unit,'(t5, 10i8)') vec(ib:ie)

   end subroutine print_int_vec

   subroutine  large_file_read(fname, outstr)
      !! This function attempt some magic to speed up the read of a
      !! large text file in order to transfer it to memory.
      use mod_constants, only: OMMP_STR_CHAR_MAX
      implicit none

      character(len=*), intent(in) :: fname
      !! File to be read
      character(len=OMMP_STR_CHAR_MAX), allocatable, intent(out) :: outstr(:)
      !! Data structure to be filled with the data from file

      character(len=:), allocatable :: buf
      character :: nlc

      integer(ip) :: inu, i, nline
      integer(ip), allocatable :: lineidx(:)
      integer(8) :: fs
      integer(ip) :: err_r, ierr

      err_r = 0

      inquire(file = fname, size = fs, iostat = ierr)
      if(fs < 0 .or. ierr > 0) then
         call fatal_error("Error while checking size of file '"//fname//"'. Cannot continue.")
      end if

      if(allocated(outstr)) then
         deallocate(outstr)
      end if

      allocate(character(len=fs) :: buf)

      open(newunit=inu, &
         file=fname, &
         form='unformatted', &
         action='read', &
         access='stream', &
         status='old')

      read(inu,pos=1,iostat=err_r) buf
      close(inu)

      if(err_r /= 0) call fatal_error("Error while reading file '"//fname//"'. Cannot continue.")
      nlc = new_line(buf(1:1))

      !nline = 1
      !do i=1, fs
      !    if(buf(i:i) == nlc) then
      !      nline = nline + 1
      !    end if
      !end do

      ! I think it's just faster to allocate the maximum number of lines
      ! that could possibly be there...
      allocate(lineidx(fs))

      nline = 2
      lineidx(1) = 0
      do i=1, fs
         if(buf(i:i) == nlc) then
            lineidx(nline) = i
            nline = nline + 1
         end if
      end do
      lineidx(nline) = fs

      allocate(outstr(nline-1))

      do i=1, nline-1
         outstr(i) = buf(lineidx(i)+1:lineidx(i+1)-1)
      end do

      deallocate(lineidx)
      deallocate(buf)
   end subroutine

end module mod_io
