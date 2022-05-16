module solvers
  use mod_memory, only: ip, rp
  implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                  !
! this fortran module contains a two of routines to iteratively solve a linear     !
! system of equations, together with the required ancillary routines.              !
! at the moment, the algorithms implemented are:                                   !
!                                                                                  !
!   - (preconditioned) conjugate gradient: as the polarization equations are       !
!     symmetric and positive definite, this is the optimal choice.                 !
!   - jacobi iterations accelerated with Pulay's direct inversion in the iterative !
!     subspace: this is a pretty robust solver that can be use for general systems !
!     and that is less sensitive to small errors in the symmetry of the matrix     !
!                                                                                  !
! both solvers require two user-provided routines with a fixed interface:          !
!                                                                                  !
!   - matvec(n,x,y): computes y = Ax, where A is the linear system matrix          !
!   - precnd(n,x,y): computes y = Mx, where M is a preconditioner.                 !
!                                                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  contains
!
  subroutine conjugate_gradient(n,lprint,tol,rhs,x,n_iter,ok,matvec,precnd)
    implicit none
!
!   preconditioned conjugated gradient solver.
!  
!   variables:
!  
!     n        : integer, input, size of the matrix
!  
!     lprint   : integer, input, printing flag.
!  
!     tol      : real, input, convergence criterion. if norm = 3, convergence is 
!                achieved when rms (grad) < tol and max |grad| < 10*tol.
!  
!     rhs      : real, dimension(n), input, right-hand side of the linear system
!  
!     x        : real, dimension(n). In input, a guess of the solution (can be zero).
!                In output, the solution
!  
!     n_iter   : integer, in input, the maximum number of iterations. In output,
!                the number of iterations needed to converge.
!  
!     ok       : logical, output, T if the solver converged, false otherwise.
!  
!     matvec   : external, subroutine to compute the required matrix-vector multiplication
!                format: subroutine matvec(n,x,y)
!  
!     precnd   : external, subroutine to apply the inverse diagonal matrix to a vector.
!                format: subroutine precnd(n,x,y)
!  
    integer(ip),               intent(in)    :: n, lprint
    real(rp),                  intent(in)    :: tol
    real(rp),    dimension(n), intent(in)    :: rhs
    real(rp),    dimension(n), intent(inout) :: x
    integer(ip),               intent(inout) :: n_iter
    logical,                   intent(inout) :: ok
    external                                 :: matvec, precnd
!  
    integer(ip)              :: it, istatus
    real(rp)                 :: rms_norm, alpha, gnew, gold, gama
    real(rp),    allocatable :: r(:), p(:), h(:), z(:)
!  
    real(rp),    parameter   :: zero = 0.0_rp, one = 1.0_rp
!
!   initialize and allocate memory:
!
    ok = .false.
    allocate (r(n), p(n), h(n), z(n), stat=istatus)
    if (istatus.ne.0) then
      write(6,*) ' conjugate_gradient: failed allocation'
      stop
    end if
!
!   compute a guess, if required:
!
    rms_norm = dot_product(x,x)
    if (rms_norm.eq.zero) call precnd(n,rhs,x)
!
!   compute the residual:
!
    call matvec(n,x,z)
    r = rhs - z
!
!   apply the preconditioner and get the first direction:
!
    call precnd(n,r,z)
    p = z
    gold = dot_product(r,z)
!
!   main loop:
!   ==========
!
    do it = 1, n_iter
!
!     compute the step:
!
      call matvec(n,p,h)
      gama = dot_product(h,p)
!
!     unlikely quick return:
!
      if (gama.eq.zero) then
        ok = .true.
        n_iter = it
        return
      end if
!
      alpha = gold/gama
      x = x + alpha*p
      r = r - alpha*h
!
!     apply the preconditioner:
!
      call precnd(n,r,z)
      gnew = dot_product(r,z)
      rms_norm = sqrt(gnew/dble(n))
!  
!     printing
!
      if ( lprint.gt.0 )  write(6,100) it, rms_norm
    100   format(t3,'iter=',i4,' residual rms norm: ', d14.4 )
!
!     check for convergence:
!
      ok = rms_norm.lt.tol
!
!     if done, exit:
!
      if (ok) exit
!
!     compute the next direction:
!
      gama = gnew/gold
      p    = gama*p + z
      gold = gnew
!
    end do
!
    n_iter = it
       
    return
  end subroutine conjugate_gradient

  subroutine jacobi_diis(n,lprint,diis_max,norm,tol,rhs,x,n_iter,ok,matvec,precnd)
    implicit none
!
!   jacobi/diis solver.
!  
!   variables:
!  
!     n        : integer, input, size of the matrix
!  
!     lprint   : integer, input, printing flag.
!  
!     diis_max : integer, input, number of points to be used for diis extrapolation
!  
!                if diis_max = 0, this is just a Jacobi solver.
!  
!     norm     : integer, input, norm to be used to evaluate convergence
!                1: max |x_new - x|
!                2: rms (x_new - x)
!                3: rms (x_new - x) and max |x_new - x|
!  
!     tol      : real, input, convergence criterion. if norm = 3, convergence is 
!                achieved when rms (x_new - x) < tol and max |x_new - x| < 10*tol.
!  
!     rhs      : real, dimension(n), input, right-hand side of the linear system
!  
!     x        : real, dimension(n). In input, a guess of the solution (can be zero).
!                In output, the solution
!  
!     n_iter   : integer, in input, the maximum number of iterations. In output,
!                the number of iterations needed to converge.
!  
!     ok       : logical, output, T if the solver converged, false otherwise.
!  
!     matvec   : external, subroutine to compute the required matrix-vector multiplication
!                format: subroutine matvec(n,x,y)
!  
!     precnd   : external, subroutine to apply the inverse diagonal matrix to a vector.
!                format: subroutine precnd(n,x,y)
!  
    integer(ip),               intent(in)    :: n, diis_max, norm, lprint
    real(rp),                  intent(in)    :: tol
    real(rp),    dimension(n), intent(in)    :: rhs
    real(rp),    dimension(n), intent(inout) :: x
    integer(ip),               intent(inout) :: n_iter
    logical,                   intent(inout) :: ok
    external                                 :: matvec, precnd
!  
    integer(ip)           :: it, nmat, istatus, lenb
    real(rp)              :: rms_norm, max_norm, tol_max
    logical               :: dodiis
!  
    real(rp), allocatable :: x_new(:), y(:), x_diis(:,:), e_diis(:,:), bmat(:,:)
    real(rp), parameter   :: zero = 0.0_rp
!  
!   DIIS extrapolation flag
!  
    dodiis =  (diis_max.ne.0)
!  
!   set tolerance
!  
    tol_max = 10.0d0 * tol
!  
!   extrapolation required
!  
    if (dodiis) then
!  
!     allocate workspaces
!  
      lenb = diis_max + 1
      allocate( x_diis(n,diis_max), e_diis(n,diis_max), bmat(lenb,lenb) , stat=istatus )
      if (istatus .ne. 0) then
        write(*,*) ' jacobi_diis: [1] failed allocation (diis)'
        stop
      endif
!      
      nmat = 1
!      
    endif
!  
!   allocate workspaces
!  
    allocate( x_new(n), y(n) , stat=istatus )
    if (istatus .ne. 0) then
      write(*,*) ' jacobi_diis: [2] failed allocation (scratch)' 
      stop
    endif
!
!   if required, compute a guess:
!
    rms_norm = dot_product(x,x)
    if (rms_norm.eq.zero) call precnd(n,rhs,x)
!  
!   Jacobi iterations
!   =================
    do it = 1, n_iter
!  
!     y = rhs - O x
      call matvec( n, x, y )
      y = rhs - y
!  
!     x_new = D^-1 y
      call precnd(n,y,x_new)
!  
!     DIIS extrapolation
!     ==================
      if (dodiis) then
!  
        x_diis(:,nmat) = x_new
        e_diis(:,nmat) = x_new - x
!  
        call diis(n,nmat,diis_max,x_diis,e_diis,bmat,x_new)
!  
      endif
!  
!     increment
      x = x_new - x
!  
!     rms/max norm of increment
      if ( norm.le.3 ) then
!  
!       compute norm
        call rmsvec( n, x, rms_norm, max_norm )
!  
!       check norm
        if ( norm.eq.1 ) then
!                
          ok = (rms_norm.lt.tol)
          
        elseif ( norm.eq.2 ) then
!                
          ok = (max_norm.lt.tol)
!          
        else 

          ok = (rms_norm.lt.tol) .and. (max_norm.lt.tol_max)
!          
        endif
!  
      endif
!  
!     printing
!
      if ( lprint.gt.0 )  write(*,100) it, rms_norm, max_norm
    100   format(t3,'iter=',i4,' residual norm (rms,max): ', 2d14.4 )
!  
!     update
      x = x_new
!  
!     EXIT Jacobi loop here
!     =====================
      if (ok) exit
!  
    enddo
!
! record number of Jacobi iterations
  n_iter = it
!
  return
!
!
  end subroutine jacobi_diis
!
  subroutine diis(n,nmat,ndiis,x,e,b,xnew)
    implicit none
!  
!   perform Pulay's direct inversion in the iterative subspace extrapolation:
!  
    integer(ip),                             intent(in)    :: n, ndiis
    integer(ip),                             intent(inout) :: nmat
    real(rp),    dimension(n,ndiis),         intent(inout) :: x, e
    real(rp),    dimension(ndiis+1,ndiis+1), intent(inout) :: b
    real(rp),    dimension(n),               intent(inout) :: xnew
!  
    integer(ip) :: nmat1, i, istatus, info
    integer(ip) :: j, k
!  
    real(rp),    allocatable :: bloc(:,:), cex(:)
    integer(ip), allocatable :: ipiv(:)
!  
    real(rp), parameter   :: zero = 0.0_rp, one = 1.0_rp
!  
!  ------------------------------------------------------------------------------
!  
    if (nmat.ge.ndiis) then
      do j = 2, nmat - 10
        do k = 2, nmat - 10
          b(j,k) = b(j+10,k+10)
        end do
      end do
      do j = 1, nmat - 10
        x(:,j) = x(:,j+10)
        e(:,j) = e(:,j+10)
      end do
      nmat = nmat - 10
    end if
    nmat1 = nmat + 1
    allocate (bloc(nmat1,nmat1),cex(nmat1), ipiv(nmat1), stat=istatus)
    if ( istatus.ne.0 ) then 
      write(*,*) 'diis: allocation failed!'
      stop
    endif
    
    call makeb(n,nmat,ndiis,e,b)
    bloc   = b(1:nmat1,1:nmat1)
    cex    = zero
    cex(1) = one
!  
    call dgesv(nmat1,1,bloc,nmat1,ipiv,cex,nmat1,info)
    if (info .ne. 0) then
!  
!     inversion failed. discard the previous points and restart.
!  
      nmat = 1
      deallocate (bloc,cex,ipiv, stat=istatus)
      if ( istatus.ne.0 ) then 
        write(*,*) 'diis: deallocation failed!'
        stop
      endif
      return
    end if
!  
    xnew = zero
    do i = 1, nmat
      xnew = xnew + cex(i+1)*x(:,i)
    end do
    nmat = nmat + 1
!  
    deallocate (bloc,cex,ipiv, stat=istatus)
    if ( istatus.ne.0 ) then 
      write(*,*) 'diis: deallocation failed!'
      stop
    endif
!  
    return
  end subroutine diis
!
  subroutine makeb(n,nmat,ndiis,e,b)
!
!   assemble the DIIS B matrix:
!  
    implicit none
    integer(ip),                             intent(in)    :: n, nmat, ndiis
    real(rp),    dimension(n,ndiis),         intent(in)    :: e
    real(rp),    dimension(ndiis+1,ndiis+1), intent(inout) :: b
!  
    integer(ip)         :: i
    real(rp)            :: bij
    real(rp), parameter :: zero = 0.0_rp, one = 1.0_rp
      
!   1st built
!  
    if (nmat.eq.1) then
!  
!         [ 0 |  1  ]
!     b = [ --+---- ]
!         [ 1 | e*e ]
!  
      b(1,1) = zero
      b(1,2) = one
      b(2,1) = one
      b(2,2) = dot_product(e(:,1),e(:,1))
!  
!   subsequent builts
!  
    else
!  
!     first, update the lagrangian line:
!  
      b(nmat+1,1) = one
      b(1,nmat+1) = one
!  
!     now, compute the new matrix elements:
!  
      do i = 1, nmat - 1
        bij = dot_product(e(:,i),e(:,nmat))
        b(nmat+1,i+1) = bij
        b(i+1,nmat+1) = bij
      end do
      b(nmat+1,nmat+1) = dot_product(e(:,nmat),e(:,nmat))
    end if
    !
    return
  end subroutine makeb
!
  subroutine rmsvec( n, v, vrms, vmax )
    implicit none
!
!   compute root-mean-square and max norms of a vector.
!
    integer(ip),               intent(in)    :: n
    real(rp),    dimension(n), intent(in)    :: v
    real(rp),                  intent(inout) :: vrms, vmax
!
    integer(ip)           :: i
    real(rp),   parameter :: zero=0.0d0
!      
!   initialize
    vrms = zero
    vmax = zero
!
!   loop over entries
    do i = 1,n
!
!     max norm
      vmax = max(vmax,abs(v(i)))
!
!     rms norm
      vrms = vrms + v(i)*v(i)
!      
    enddo
!
!   the much neglected square root
    vrms = sqrt(vrms/dble(n))
!    
    return
!      
!      
  end subroutine rmsvec
!
end module solvers