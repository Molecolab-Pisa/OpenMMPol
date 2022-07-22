module mod_bicubic_interp
    use mod_memory, only: ip, rp

    implicit none
    
    contains

    subroutine control_points_sphere(indata, cpt, n, m)
        use mod_memory, only: mallocate, mfree

        implicit none

        real(rp), intent(in) :: indata(:)
        real(rp), intent(out) :: cpt(:)
        integer(ip), intent(in) :: n, m
        
        real(rp), allocatable :: phi(:,:), work(:), pad_data(:)
        integer(ip) :: i, j, k, jj, ld, info
        integer(ip), allocatable :: ipiv(:)

        ld = (n+2)*(m+2)
        call mallocate('control_points_sphere [phi]', ld, ld, phi)
        call mallocate('control_points_sphere [work]', ld*ld, work)
        call mallocate('control_points_sphere [pad_data]', ld, pad_data)
        call mallocate('control_points_sphere [ipiv]', ld, ipiv)

        pad_data = 0.0
        pad_data(1:n*m) = indata

        ! Populate Phi matrix

        ! (interpolation)
        do j=0, m-1
            do i=0, n-1
                phi(i+j*n+1, i+j*(n+2)+1) = 1.0
                phi(i+j*n+1, i+j*(n+2)+2) = 4.0
                phi(i+j*n+1, i+j*(n+2)+3) = 1.0

                phi(i+j*n+1, i+(j+1)*(n+2)+1) = 4.0
                phi(i+j*n+1, i+(j+1)*(n+2)+2) = 16.0
                phi(i+j*n+1, i+(j+1)*(n+2)+3) = 4.0
                
                phi(i+j*n+1, i+(j+2)*(n+2)+1) = 1.0
                phi(i+j*n+1, i+(j+2)*(n+2)+2) = 4.0
                phi(i+j*n+1, i+(j+2)*(n+2)+3) = 1.0
            end do
        end do

        ! x borders
        do i=0, m+1
            phi(n*m+i+1, modulo(i*(n+2), ld)+1) = &
            phi(n*m+i+1, modulo(i*(n+2), ld)+1) - 3.0
            phi(n*m+i+1, modulo(i*(n+2)+2, ld)+1) = &
            phi(n*m+i+1, modulo(i*(n+2)+2, ld)+1) + 3.0
            
            phi(n*m+i+1, modulo((i+1)*(n+2), ld)+1) = &
            phi(n*m+i+1, modulo((i+1)*(n+2), ld)+1) - 12.0
            phi(n*m+i+1, modulo((i+1)*(n+2)+2, ld)+1) = &
            phi(n*m+i+1, modulo((i+1)*(n+2)+2, ld)+1) + 12.0

            phi(n*m+i+1, modulo((i+2)*(n+2), ld)+1) = &
            phi(n*m+i+1, modulo((i+2)*(n+2), ld)+1) - 3.0
            phi(n*m+i+1, modulo((i+2)*(n+2)+2, ld)+1) = &
            phi(n*m+i+1, modulo((i+2)*(n+2)+2, ld)+1) + 3.0
        end do

        do i=0, m+1
            phi(n*m+i+m+3, modulo((i+1)*(n+2)-1, ld)+1) = &
            phi(n*m+i+m+3, modulo((i+1)*(n+2)-1, ld)+1) - 3.0
            phi(n*m+i+m+3, modulo((i+1)*(n+2)-3, ld)+1) = &
            phi(n*m+i+m+3, modulo((i+1)*(n+2)-3, ld)+1) + 3.0
            
            phi(n*m+i+m+3, modulo((i+2)*(n+2)-1, ld)+1) = &
            phi(n*m+i+m+3, modulo((i+2)*(n+2)-1, ld)+1) - 12.0
            phi(n*m+i+m+3, modulo((i+2)*(n+2)-3, ld)+1) = &
            phi(n*m+i+m+3, modulo((i+2)*(n+2)-3, ld)+1) + 12.0
            
            phi(n*m+i+m+3, modulo((i+3)*(n+2)-1, ld)+1) = &
            phi(n*m+i+m+3, modulo((i+3)*(n+2)-1, ld)+1) - 3.0
            phi(n*m+i+m+3, modulo((i+3)*(n+2)-3, ld)+1) = &
            phi(n*m+i+m+3, modulo((i+3)*(n+2)-3, ld)+1) + 3.0
        end do

        do i=0, n-1
            phi(n*m+2*(m+2)+i+1,i+1)=-3.0
            phi(n*m+2*(m+2)+i+1,i+2)=-12.0
            phi(n*m+2*(m+2)+i+1,i+3)=-3.0

            phi(n*m+2*(m+2)+i+1,2*(n+2)+i+1)=3.0
            phi(n*m+2*(m+2)+i+1,2*(n+2)+i+2)=12.0
            phi(n*m+2*(m+2)+i+1,2*(n+2)+i+3)=3.0

            phi(n*m+2*(m+2)+i+1,(m-1)*(n+2)+i+1)=3.0
            phi(n*m+2*(m+2)+i+1,(m-1)*(n+2)+i+2)=12.0
            phi(n*m+2*(m+2)+i+1,(m-1)*(n+2)+i+3)=3.0

            phi(n*m+2*(m+2)+i+1,(m+1)*(n+2)+i+1)=-3.0
            phi(n*m+2*(m+2)+i+1,(m+1)*(n+2)+i+2)=-12.0
            phi(n*m+2*(m+2)+i+1,(m+1)*(n+2)+i+3)=-3.0
        end do

        do i=0, n-1
            phi(n*m+2*(m+2)+n+i+1,i+1)=6.0
            phi(n*m+2*(m+2)+n+i+1,i+2)=24.0
            phi(n*m+2*(m+2)+n+i+1,i+3)=6.0

            phi(n*m+2*(m+2)+n+i+1,(n+2)+i+1)=-12.0
            phi(n*m+2*(m+2)+n+i+1,(n+2)+i+2)=-48.0
            phi(n*m+2*(m+2)+n+i+1,(n+2)+i+3)=-12.0

            phi(n*m+2*(m+2)+n+i+1,2*(n+2)+i+1)=6.0
            phi(n*m+2*(m+2)+n+i+1,2*(n+2)+i+2)=24.0
            phi(n*m+2*(m+2)+n+i+1,2*(n+2)+i+3)=6.0

            phi(n*m+2*(m+2)+n+i+1,(m-1)*(n+2)+i+1)=-6.0
            phi(n*m+2*(m+2)+n+i+1,(m-1)*(n+2)+i+2)=-24.0
            phi(n*m+2*(m+2)+n+i+1,(m-1)*(n+2)+i+3)=-6.0

            phi(n*m+2*(m+2)+n+i+1,(m)*(n+2)+i+1)=12.0
            phi(n*m+2*(m+2)+n+i+1,(m)*(n+2)+i+2)=48.0
            phi(n*m+2*(m+2)+n+i+1,(m)*(n+2)+i+3)=12.0

            phi(n*m+2*(m+2)+n+i+1,(m+1)*(n+2)+i+1)=-6.0
            phi(n*m+2*(m+2)+n+i+1,(m+1)*(n+2)+i+2)=-24.0
            phi(n*m+2*(m+2)+n+i+1,(m+1)*(n+2)+i+3)=-6.0
        end do
            
        ! Invert phi matrix
        call dgetrf(ld, ld, phi, ld, ipiv, info)
        call dgetri(ld, phi, ld, ipiv, work, ld, info)

        ! Now compute the values of Qz
        cpt = 0.0

        do i=1, ld
            do k=1, n
                do j=1, m
                    jj = (j-1)*n+k
                    cpt(i) = cpt(i) + 36.0*phi(i,(k-1)*m+j)*pad_data(jj)
                end do
            end do
        end do
                
        call mfree('control_points_sphere [work]', work)
        call mfree('control_points_sphere [phi]', phi)
        call mfree('control_points_sphere [pad_data]', pad_data)
        call mfree('control_points_sphere [ipiv]', ipiv)

    end subroutine control_points_sphere

    subroutine evaluate_bicubic(x, y, z, xmin, xmax, nx, ymin, ymax, ny, qz)

        implicit none

        real(rp), intent(in) :: x, y, xmin, xmax, ymin, ymax, qz(:)
        real(rp), intent(out) :: z
        integer(ip) :: nx, ny

        integer(ip) :: pu, pv, i, j, ii, jj
        real(rp) :: u, v, deltax, deltay, q(4,4), vv(1,4), uu(4, 1), res(1,1)

        real(rp), parameter, dimension(4,4) :: M = reshape([-1,  3, -3,  1, &
                                                             3, -6,  0,  4, &
                                                            -3,  3,  3,  1, &
                                                             1,  0,  0,  0],&
                                                           shape(M))

        deltax = (xmax-xmin)/(nx-1)
        deltay = (ymax-ymin)/(ny-1)

        pu = floor((x-xmin)/deltax)
        pv = floor((y-ymin)/deltay)
        
        u = (x-(pu*deltax+xmin))/deltax
        v = (y-(pv*deltay+ymin))/deltay

        vv = reshape([v*v*v, v*v, v, 1.0_rp], shape(vv))
        uu = reshape([u*u*u, u*u, u, 1.0_rp], shape(uu))
        
        do i=1, 4
            ii = (pv+i-1) 
            if(ii >= ny+2) ii = ii - ny+2
            ii = ii * (nx+2) 
            do j=1, 4
                jj = (pu+j) 
                if(jj >= nx+2) jj = jj - nx+2
                q(j, i) = qz(ii+jj)
            end do
        end do

        q=transpose(q)
        res = matmul(matmul(matmul(matmul(vv, M), q), transpose(M)), uu)
        z = res(1,1) / 36.0

    end subroutine
end module
