#include "version.h"

module mod_io
    !! Unified Input/Output handling across the code.
    
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
    public :: set_verbosity, ommp_message, fatal_error, ommp_version
    public :: print_header
    public :: print_matrix, print_int_vec

    interface print_matrix
        !! Interface for matrix printing function
        module procedure d1_print_matrix
        module procedure d2_print_matrix
    end interface print_matrix

    contains
    
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
        use mod_constants, only: OMMP_VERBOSE_DEBUG, &
                                 OMMP_VERBOSE_HIGH, &
                                 OMMP_VERBOSE_LOW

        implicit none

        character(len=*), intent(in) :: s
        !! Message to be printed
        character(len=*), intent(in), optional :: logpre
        !! String that explains the type of message, if missing a default type
        !! is assigned based on requested verbosity level
        integer(ip), intent(in) :: level
        !! Requested verbosity level
        integer(ip), intent(in), optional :: u
        !! Output unit for the message, if missing, unit 6/stdout is used.

        integer(ip) :: outunit
        character(len=12) :: pre
        
        if(level > verbose) return

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

        write(outunit, '(A6, A12, " ", A)') '[OMMP]', pre, trim(s)
    end subroutine ommp_message
    
    subroutine fatal_error(message)
        !! Prints a message and exit from the program. This
        !! function should be used in all the conditions 
        !! where the program cannot proceed.

        use mod_constants, only: OMMP_VERBOSE_LOW
        implicit none
      
        character (len=*), intent(in) :: message
        !! Message to print before the program termination

        !write(6, '(t3,a)') message
        call ommp_message(message, OMMP_VERBOSE_LOW, 'stop')
        call ommp_message("Unrecoverable error in openMMPol &
                          &library. Exiting.", OMMP_VERBOSE_LOW, &
                          'stop')
        !TODO call mmpol_terminate()

        stop
    end subroutine fatal_error

    subroutine ommp_version(v)
        integer(ip), intent(in) :: v

        call ommp_message("OpenMMPol version: "//OMMP_VERSION, &
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

   1000 format(t3,a)
   1010 format(t3,5i16)
   1020 format(t3,5f16.8)

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

    subroutine print_int_vec(label, vec, dosort, ofunit, ibeg, iend)
        !! Print an array of integers in a well formatted way.
        use mod_utils, only: sort_ivec
        
        implicit none
        
        character(len=*), intent(in) :: label
        !! Label to print before the array
        integer(ip), dimension(:), intent(in) :: vec
        !! Integer vector to be printed
        logical, optional :: dosort
        !! If present and .true., print the sorted array
        !! otherwise [[iof_mmpol]] is used.
        integer(ip), intent(in), optional :: ofunit
        !! If present specify the unit where the array should be printed, 
        integer(ip), intent(in), optional :: ibeg
        !! Index of the first element to be printed, if 0 or absent the first 
        !! element is used
        integer(ip), intent(in), optional :: iend
        !! Index of the last element to be printed, if 0 or absent the last 
        !! element of the array is used
        logical :: sort

        
        integer(ip) :: ib, ie, out_unit

        integer(ip), allocatable, dimension(:) :: sorted_vec

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

        !! TODO extend the scope of this function...
        if(sort) then
            call sort_ivec(vec(ib:ie), sorted_vec)
            write(out_unit,'(t5, 10i8)') sorted_vec
        else
            write(out_unit,'(t5, 10i8)') vec(ib:ie)
        end if

        return

    end subroutine print_int_vec

end module mod_io
