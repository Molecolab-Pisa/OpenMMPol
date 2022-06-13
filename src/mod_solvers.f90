module mod_solvers
  use mod_memory, only: ip, rp
  use mod_constants, only: OMMP_VERBOSE_HIGH, &
                           OMMP_VERBOSE_LOW, &
                           OMMP_VERBOSE_DEBUG
  use mod_mmpol, only: verbose, fatal_error

  implicit none
  private
  
  real(rp), parameter :: OMMP_DEFAULT_SOLVER_TOL = 1e-8_rp
  integer(ip), parameter :: OMMP_DEFAULT_SOLVER_ITER = 200

  public :: inversion_solver, conjugate_gradient_solver, jacobi_diis_solver

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
  subroutine inversion_solver(n, rhs, x, tmat)
      use mod_mmpol, only : pol_atoms
      use mod_memory, only: mallocate, mfree
      !
      ! Induce dipoles on the polarizable atoms. The dipoles are induced
      ! by the electric field of the static multipoles (Later change to 
      ! the total electric field including the QM system). For the Amoeba
      ! force-field only polarization field and p-dipoles are considered.
      !
      
      implicit none
      
      integer(ip), intent(in) :: n
      !! Size of the matrix
      real(rp), dimension(n), intent(in) :: rhs
      !! Right hand side of the linear system
      real(rp), dimension(n), intent(out) :: x
      !! In output the solution of the linear system
      real(rp), dimension(n, n), intent(in) :: tmat
      !! Polarization matrix TODO


      integer(ip) :: info
      integer(ip), dimension(:), allocatable :: ipiv
      real(rp), dimension(:), allocatable :: work
      real(rp), dimension(:,:), allocatable :: TMatI
      
      call mallocate('inversion_solver [TMatI]', n, n, TMatI)
      call mallocate('inversion_solver [work]', n, work)
      call mallocate('inversion_solver [ipiv]', n, ipiv)
      
      ! Initialize inverse polarization matrix
      TMatI = TMat
      
      !Compute the inverse of TMat
      call dgetrf(n, n, TMatI, n, iPiv, info)
      call dgetri(n, TMatI, n, iPiv, Work, n, info)
      
      ! Calculate dipoles with matrix inversion
      call dgemm('N', 'N', n, 1, n, 1.0_rp, TMatI, n, rhs, n, 0.0_rp, x, n)
      
      call mfree('inversion_solver [TMatI]', TMatI)
      call mfree('inversion_solver [work]', work)
      call mfree('inversion_solver [ipiv]', ipiv)
    
  end subroutine inversion_solver

  subroutine conjugate_gradient_solver(n, rhs, x, matvec, precnd, arg_tol, arg_n_iter)
    !! Conjugate gradient solver (TODO)
    ! TODO add more printing
    use mod_constants, only: eps_rp
    use mod_memory, only: mallocate, mfree
    
    implicit none
    
    integer(ip), intent(in) :: n
    !! Size of the matrix
    real(rp), intent(in), optional :: arg_tol
    !! Optional convergence criterion in input, if not present
    !! OMMP_DEFAULT_SOLVER_TOL is used.
    real(rp) :: tol
    !! Convergence criterion, it is required that RMS norm < tol
    
    integer(ip), intent(in), optional :: arg_n_iter
    !! Optional maximum number of iterations for the solver, if not present
    !! OMMP_DEFAULT_SOLVER_ITER is used.
    integer(ip) :: n_iter
    !! Maximum number of iterations for the solver 

    real(rp), dimension(n), intent(in) :: rhs
    !! Right hand side of the linear system
    real(rp), dimension(n), intent(inout) :: x
    !! In input, initial guess for the solver, in output the solution
    external :: matvec
    !! Routine to perform matrix-vector product
    external :: precnd
    !! Preconditioner routine

    integer(ip) :: it
    real(rp) :: rms_norm, alpha, gnew, gold, gama
    real(rp), allocatable :: r(:), p(:), h(:), z(:)

    ! Optional arguments handling
    if(present(arg_tol)) then
        tol = arg_tol
    else
        tol = OMMP_DEFAULT_SOLVER_TOL
    end if
    
    if(present(arg_n_iter)) then
        n_iter = arg_n_iter
    else
        n_iter = OMMP_DEFAULT_SOLVER_ITER
    end if

    if(verbose >= OMMP_VERBOSE_LOW) then
        write(6, *) "Solving linear system with CG solver"
        write(6, *) "Max iter:", n_iter
        write(6, *) "Tolerance: ", tol
    end if

    call mallocate('conjugate_gradient_solver [r]', n, r)
    call mallocate('conjugate_gradient_solver [p]', n, p)
    call mallocate('conjugate_gradient_solver [h]', n, h)
    call mallocate('conjugate_gradient_solver [z]', n, z)
    
    ! compute a guess, if required:
    rms_norm = dot_product(x,x)
    if(rms_norm < eps_rp) call precnd(n, rhs, x)

    ! compute the residual:
    call matvec(n, x, z)
    r = rhs - z

    ! apply the preconditioner and get the first direction:
    call precnd(n, r, z)
    p = z
    gold = dot_product(r, z)
    
    do it = 1, n_iter
      ! compute the step:
      call matvec(n, p, h)
      gama = dot_product(h, p)

      ! unlikely quick return:
      if(abs(gama) < eps_rp) then
        write(6, *) "Direction vector with zero norm, exiting iterative solver."
        exit
      end if

      alpha = gold / gama
      x = x + alpha * p
      r = r - alpha * h

      ! apply the preconditioner:
      call precnd(n, r, z)
      gnew = dot_product(r, z)
      rms_norm = sqrt(gnew/dble(n))
      
      if(verbose >= OMMP_VERBOSE_HIGH) then
          write(6, "(t3,'iter=',i4,' residual rms norm: ', d14.4)") it, rms_norm
      end if
      
      ! Check convergence
      if(rms_norm < tol) then
          write(6, *) "Required convergence threshold reached, exiting iterative solver."
          exit
      end if

      ! compute the next direction:
      gama = gnew/gold
      p    = gama*p + z
      gold = gnew
    end do
    
    if(rms_norm > tol .and. abs(gama) > eps_rp) then
        call fatal_error("Iterative solver did not converged")
    end if

    call mfree('conjugate_gradient_solver [r]', r)
    call mfree('conjugate_gradient_solver [p]', p)
    call mfree('conjugate_gradient_solver [h]', h)
    call mfree('conjugate_gradient_solver [z]', z)
  
  end subroutine conjugate_gradient_solver

  subroutine jacobi_diis_solver(n,lprint,diis_max,norm,tol,rhs,x,n_iter,ok,matvec,precnd)
    use mod_constants, only: eps_rp
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
    if (rms_norm < eps_rp ) call precnd(n,rhs,x)
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
  end subroutine jacobi_diis_solver
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
end module mod_solvers
