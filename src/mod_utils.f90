module mod_utils
    !! This module contains some very generic utils for string manipulation,
    !! or very basic computational/mathematic operation. It should not depend
    !! on any module except from [[mod_memory]].
    use mod_memory, only: ip
    use mod_constants, only: OMMP_STR_CHAR_MAX

    implicit none
    private

    public :: sort_ivec, sort_ivec_inplace, skip_lines
    public :: starts_with_alpha, isreal, isint, tokenize, &
              count_substr_occurence, str_to_lower
    public :: cyclic_spline, compute_bicubic_interp

    contains

    function str_to_lower(s)
        !! Convert string in input from upper case to lower case and return
        !! the lower case string as output.

        implicit none

        character(len=*), intent(in) :: s
        !! String to be converted in lowercase
        character(len(s)) :: str_to_lower
        !! String converted to lowercase

        integer :: ic, i

        character(26), parameter :: capital = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        character(26), parameter :: lower = 'abcdefghijklmnopqrstuvwxyz'

        str_to_lower = s
        do i=1, len_trim(s)
            ic = index(capital, s(i:i))
            if(ic > 0) str_to_lower(i:i) = lower(ic:ic)
        end do
        
        return

    end function str_to_lower
   
    function count_substr_occurence(s, c)
        !! Count the number of occurence of substring c in string s, and return
        !! the number of occurence, if c is not contained in s, zero is returned.
        implicit none

        character(len=*), intent(in) :: s
        !! String where to search the substring
        character(len=*), intent(in) :: c
        !! Substring to search
        integer(ip) :: count_substr_occurence, i, lens, lenc
        
        count_substr_occurence = 0
        lens = len(s)
        lenc = len(c)

        do i=1, lens-lenc+1
            if(s(i:i+lenc-1) == c) then
                count_substr_occurence = count_substr_occurence+1
            end if
        end do
        
        return
    end function

    function starts_with_alpha(s)
        !! Decide if a string starts with a letter or not.
        implicit none

        character(len=*), intent(in) :: s
        !! String to analyze

        logical :: starts_with_alpha

        starts_with_alpha = (scan(s(1:1), &
        'qwertyuiopasdfghjklzxcvbnmQWERTYUIOPASDFGHJKLZXCVBNM') /= 0)
        return
    end function

    function isint(s)
        !! Decide if a string can be interpreted as an integer or not

        implicit none

        character(len=*), intent(in) :: s
        !! String to analyze
        logical :: isint
        
        isint = (verify(s, '+-1234567890') == 0)
        return
    end function

    function isreal(s)
        !! Decide if a string can be interpreted as a real
        implicit none

        character(len=*), intent(in) :: s
        !! String to analyze
        logical :: isreal

        isreal = (verify(s, '+-1234567890.') == 0)
        isreal = isreal .and. (scan(s, '.') /= 0)
        return
    end function

    function tokenize(s, ib, ntok)
        !! This function is used to subsequently break a string into tokens.
        !! Tokens separators are any number of spaces.   
        !! If just the string is provided, the function returns the position of
        !! the first printable character;   
        !! If also ib is provided it saves the position of the first printable 
        !! character after position ib (or ib itself) in ib and return the 
        !! position of the last printable character before the first space after
        !! ib. If ntok is specified instead of a single token, ntok are 
        !! returned. In case of last token hitten -1 is returned.   
        !! To divide a string follow the following scheme:   
        !! 1. ib = tokenize(s)   
        !! 2. ie = tokenize(s, ib)   
        !! 3. tok1 = s(ib:ie)   
        !! 4a. ib = ib+1   
        !! 4b. ie = tokenize(s, ib)   
        !! 5. tok2 = s(ib:ie)   

        use mod_memory, only: ip
        implicit none

        character(len=OMMP_STR_CHAR_MAX), intent(in) :: s
        !! String to subdivide in token
        integer(ip), intent(inout), optional :: ib
        !! Index where to start token research (input)/Index where token 
        !! begins (output)
        integer(ip), intent(in), optional :: ntok
        !! Number of token to be extracted
        integer(ip) :: tokenize
        !! Index where token ends.

        integer(ip) :: i, slen, ich, itok

        slen = len(s)
        ! Default return for end of string
        tokenize = -1
        if(present(ib)) then

            ! This is a very unreasonable case
            if(ib > slen) return
        
            do i=ib, slen
                ! Search the first valid char and save it in ib
                ich = iachar(s(i:i))
                if(ich > 32 .and. ich /= 127) exit
            end do
            ib = i
            if(ib >= slen) return
            
            if(present(ntok)) then
                itok = ntok
            else
                itok = 1
            end if

            do i=ib+1, slen
                ich = iachar(s(i:i))
                if(ich <= 32 .or. ich == 127) then 
                    ich = iachar(s(i-1:i-1))
                    if(ich > 32 .and. ich /= 127) itok = itok - 1
                end if
                if(itok == 0) exit
            end do
            tokenize = i-1
        else 
            do i=1, slen
                ich = iachar(s(i:i))
                if(ich > 32 .and. ich /= 127) exit
            end do
            tokenize = i
        end if
    end function
    
    subroutine skip_lines(f, n)
        !! Skips n lines while reading an input file
        use mod_memory, only: ip
        
        implicit none

        integer(ip), intent(in) :: f
        !! unit file
        integer(ip), intent(in) :: n
        !! number of line to be skipped

        integer(ip) :: i

        do i=1, n
            read(f, *)
        end do
    end subroutine skip_lines


    subroutine sort_ivec(iv, ov)
        !! This is a simple -- and unefficient -- routine to sort a vector of 
        !! integers.
        !! It is just used during some output to simplify comparisons with older 
        !! version of the code.   
        !! @warning This function should not be used in efficiency-critical
        !! part of the code!
        use mod_memory, only: mallocate

        implicit none

        integer(ip), dimension(:), intent(in) :: iv
        !! Input vector
        integer(ip), allocatable, intent(out) :: ov(:)
        !! Output, sorted vector
        
        integer(ip) :: i, imin(1)
        logical, allocatable :: mask(:)

        call mallocate('sort_ivec', size(iv), ov)
        allocate(mask(size(iv)))
        mask = .true.

        do i = 1, size(iv)
            imin = minloc(iv, mask)
            ov(i) = iv(imin(1))
            mask(imin(1)) = .false.
        end do

        deallocate(mask)

    end subroutine sort_ivec
    
    subroutine sort_ivec_inplace(iv)
        !! Inplace equivalent of [[sort_ivec]] routine.

        implicit none

        integer(ip), dimension(:), intent(inout) :: iv
        !! Vector to be reordered in place
        integer(ip), allocatable :: ov(:)

        call sort_ivec(iv, ov)
        iv = ov

    end subroutine sort_ivec_inplace

    subroutine compute_bicubic_interp(x, y, z, nx, ny, xgrd, ygrd, &
                                      v, vx, vy, vxy)
        !! Evaluate the z value at position (x, y) of a surface built as a 
        !! bicubic spline interpolating the points     
        !! (xgrd(\(x_i\), \(y_i\)), ygrd(\(x_i\), \(y_i\)),
        !! vxgrd(\(x_i\), \(y_i\))).     
        !! In order to do so also the derivatives of the surface at the points
        !! of the grid arre needed and in particular if we consider the surface
        !! points as
        !! \[(x_i, y_i, V(x_i, y_i))\]
        !! also value of \(\frac{\partial V}{\partial x}\),
        !! \(\frac{\partial V}{\partial y}\) and 
        !! \(\frac{\partial^2 V}{\partial x \partial y}\) at grid points
        !! are needed.
        
        use mod_memory, only: ip, rp

        implicit none

        real(rp), intent(in) :: x, y
        !! Coordinates at which the z value of the surface should be computed
        real(rp), intent(out) :: z
        !! The z value of the surface
        integer(ip), intent(in) :: nx, ny
        !! Number of grid points along x and y direction
        real(rp), dimension(ny,nx), intent(in) :: xgrd
        !! X-coordinate value of points on the grid
        real(rp), dimension(ny,nx), intent(in) :: ygrd 
        !! Y-coordinate value of points on the grid
        real(rp), dimension(ny,nx), intent(in) :: v
        !! V(x, y) of points on the grid
        real(rp), dimension(ny,nx), intent(in) :: vx 
        !! \(\frac{\partial V}{\partial x}\) of points on the grid
        real(rp), dimension(ny,nx), intent(in) :: vy 
        !! \(\frac{\partial V}{\partial y}\) of points on the grid
        real(rp), dimension(ny,nx), intent(in) :: vxy
        !! \(\frac{\partial^2 V}{\partial x \partial y}\) of points on the grid

        integer(ip) :: ix, iy, ii, jj
        logical :: done

        real(rp), parameter :: A(4,4) = reshape([ 1.0,  0.0, -3.0,  2.0, &
                                                  0.0,  0.0,  3.0, -2.0, &
                                                  0.0,  1.0, -2.0,  1.0, &
                                                  0.0,  0.0, -1.0,  1.0],&
                                                 shape(A))

        real(rp) :: alpha(4,4), f(4,4), xx(4), yy(4), &
                    deltax, deltay, dx, dy

        done = .false.
        ix = 0
        iy = 0
        do ix=1, nx-1
            do iy=1, ny-1
                if(x < xgrd(ix+1,iy) .and. x > xgrd(ix,iy) .and. &
                   y > ygrd(ix,iy) .and. y < ygrd(ix,iy+1)) then 
                    done = .true.
                end if
                
                if(done) exit
            end do
            if(done) exit
        end do
        
        if(.not. done) then
            write(*, *) "Cannot find the required point on the grid"
            stop
        end if
        
        dx = xgrd(ix+1,iy)-xgrd(ix,iy)
        dy = ygrd(ix,iy+1)-ygrd(ix,iy)
        do ii=0, 1
            do jj=0, 1
                f(jj+1,ii+1) = v(ix+jj,iy+ii)
                f(jj+1,ii+3) = vy(ix+jj,iy+ii) * dy
                f(jj+3,ii+1) = vx(ix+jj,iy+ii) * dx
                f(jj+3,ii+3) = vxy(ix+jj,iy+ii) * dx*dy
            end do
        end do

        alpha = matmul(A, matmul(f, transpose(A)))

        xx(1) = 1.0
        yy(1) = 1.0
        deltax = (x-xgrd(ix,iy)) / dx
        deltay = (y-ygrd(ix,iy)) / dy
        do ii=2, 4
            xx(ii) = xx(ii-1) * deltax
            yy(ii) = yy(ii-1) * deltay
        end do

        z = 0.0
        do ii=1,4
            do jj=1, 4
                z = z + alpha(ii, jj) * yy(jj) * xx(ii)
            end do
        end do
    end subroutine compute_bicubic_interp

    subroutine cyclic_spline(n, x, y, a, b, c, d)
        !! Compute the cyclic interpolating cubic spline (2D) that passes for  
        !! points \((x_i, y_i)\). Each segment is described by the curve:
        !! \[S_i(x) = a_i + b_i(x-x_i) + c_i (x-x_i)^2 + d_i(x-x_i)^3\]
        !! The algorithm used to compute the coefficients is taken from 
        !! "Numerical Algorithm with C" - Gisela ENGELN-MULLGES Frank UHLIG
        !! (10.1007/978-3-642-61074-5) 

        use mod_memory, only: ip, rp, mallocate, mfree
        implicit none

        integer(ip), intent(in) :: n 
        !! Dimension of the data series
        real(rp), intent(in) :: x(n)
        !! X-values of input data
        real(rp), intent(in) :: y(n) 
        !! Y-values of input data

        real(rp), intent(out) :: a(n), b(n), c(n), d(n)
        !! Coefficients of the cubic spline in each segment of the spline

        real(rp), allocatable :: h(:), g(:), am(:,:), work(:)
        integer(ip), allocatable :: ipiv(:)
        integer(ip) :: i, info

        ! check that x_i are in increasing order
        do i=1, n-1
            if(x(i) >= x(i+1)) then
                write(*, *) "Data provided to cyclic_spline function should be &
                            &in strictly increasing order."
                stop
            end if
        end do

        ! Assemble the matrix and the rhs of the linear system
        call mallocate('cyclic_spline [h]', n, h)
        call mallocate('cyclic_spline [g]', n-1, g)
        call mallocate('cyclic_spline [am]', n-1, n-1, am)
        call mallocate('cyclic_spline [work]', n-1, work)
        call mallocate('cyclic_spline [ipiv]', n-1, ipiv)

        do i=1, n-1
            h(i) = x(i+1) - x(i)
            a(i) = y(i)            
        end do
        h(n) = h(1)
        a(n) = a(1)

        do i=2, n-1
            g(i-1) = (a(i+1)-a(i))/h(i) - (a(i)-a(i-1))/h(i-1)
        end do
        g(n-1) = (a(2)-a(n))/h(n) - (a(n)-a(n-1))/h(n-1)
        g = 3.0 * g
        
        am = 0.0
        do i=1, n-1
            am(i, i) = 2 * (h(i) + h(i+1))
            if(i > 1) then
                am(i-1, i) = h(i)
            else
                am(n-1, i) = h(i)
            end if
            if(i < n-1) then
                am(i+1, i) = h(i+1)
            else
                am(1, i) = h(i+1)
            end if
        end do

        ! Now solve the linear system am @ c = g that is c = am^-1 g
        call dgetrf(n-1, n-1, am, n-1, ipiv, info)
        call dgetri(n-1, am, n-1, ipiv, work, n-1, info)
        call dgemm('N', 'N', n-1, 1, n-1, 1.0_rp, am, n-1, &
                    g, n-1, 0.0_rp, c(2:n), n-1)
        c(1) = c(n)

        do i=1, n-1
            b(i) = (a(i+1)-a(i))/h(i) - h(i)*(c(i+1)+2.0*c(i))/3.0
            d(i) = (c(i+1)-c(i))/(3.0*h(i))
        end do
        b(n) = (a(2)-a(n))/h(n) - h(n)*(c(2)+2.0*c(n))/3.0
        d(n) = (c(2)-c(n))/(3.0*h(n))

        call mfree('cyclic_spline [g]', g)
        call mfree('cyclic_spline [h]', h)
        call mfree('cyclic_spline [am]', am)
        call mfree('cyclic_spline [work]', work)
        call mfree('cyclic_spline [ipiv]', ipiv)

    end subroutine

end module mod_utils
