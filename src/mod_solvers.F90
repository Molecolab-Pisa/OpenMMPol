module mod_solvers
    !! Module that contains the routines used to solve the polarization linear 
    !! system \(\mathbf A \mathbf x = \mathbf B\).
    !! Currently three methods are implemented:    
    !! 1.  __matrix inversion__;    
    !! 2.  __(preconditioned) conjugate gradients__ - since polarization equations 
    !!     are symmetric and positive definite, this is the optimal choice;    
    !! 3.  __jacobi iterations__ accelerated with Pulay's direct inversion 
    !!     in the iterative subspace (__DIIS__): this is a pretty robust solver that 
    !!     can be use for general systems and that is less sensitive to small 
    !!     errors in the symmetry of the matrix.    
    !!       
    !! Iterative solvers need two additional routines to be passed as arguments,
    !! namely matvec that computes a generic product 
    !! \(\mathbf y = \mathbf A \mathbf v\)
    !! and precond that computes \(\mathbf y = \mathbf M \mathbf v\), where 
    !! \(M\) is a precontioner

    use mod_memory, only: ip, rp
    use mod_constants, only: OMMP_VERBOSE_HIGH, &
                             OMMP_VERBOSE_LOW, &
                             OMMP_VERBOSE_DEBUG, &
                             OMMP_STR_CHAR_MAX
    use mod_io, only: ommp_message, fatal_error
    use mod_electrostatics, only: ommp_electrostatics_type

    implicit none
    private
  
    real(rp), parameter :: OMMP_DEFAULT_SOLVER_TOL = 1e-8_rp
    !! Default tolerance for iterative solvers
    integer(ip), parameter :: OMMP_DEFAULT_SOLVER_ITER = 200
    !! Default maximum number of iteration for iterative solvers
    integer(ip), parameter :: OMMP_DEFAULT_DIIS_MAX_POINTS = 20
    !! Default maximum number of points in DIIS extrapolation

    public :: inversion_solver, conjugate_gradient_solver, jacobi_diis_solver, &
              OMMP_DEFAULT_SOLVER_TOL, cp_inversion_solver, block_conjugate_gradient_solver, &
              batched_block_cg_solver, deflated_subspace_cg_solver

    contains
    
    subroutine inversion_solver(n, rhs, x, tmat)
        !! Solve the linear system directly inverting the matrix:
        !! $$\mathbf A \mathbf x = \mathbf B $$
        !! $$ \mathbf x = \mathbf A ^-1 \mathbf B $$
        !! This is highly unefficient and should only be used for testing 
        !! other methods of solution.

        use mod_memory, only: mallocate, mfree
        
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

    subroutine cp_inversion_solver(n, nrhs, rhs, x, tmat)
        !! Solve the linear system A X = B for many right-hand-sides B at
        !! once (X, B are n x nrhs), factoring A only once (LU, dgetrf)
        !! and reusing the factorization for all columns (dgetrs).
        !!
        !! Meant for the coupled-perturbed (CP) equations of the
        !! polarization Hessian: the same polarization matrix tmat is
        !! solved against one right-hand-side per Cartesian perturbation
        !! (3*mm_atoms of them), so factoring once and solving all of them
        !! together is the whole point -- unlike inversion_solver (which
        !! is only ever used for a single right-hand-side), redoing the
        !! factorization per column here would turn an O(n^3) solve into
        !! an O(n^4) one.

        use mod_memory, only: mallocate, mfree

        implicit none

        integer(ip), intent(in) :: n
        !! Size of the matrix
        integer(ip), intent(in) :: nrhs
        !! Number of right-hand-sides (columns of rhs/x)
        real(rp), dimension(n, nrhs), intent(in) :: rhs
        !! Right hand sides of the linear system
        real(rp), dimension(n, nrhs), intent(out) :: x
        !! In output the solutions of the linear system, one per column
        real(rp), dimension(n, n), intent(in) :: tmat
        !! Polarization matrix

        integer(ip) :: info
        integer(ip), dimension(:), allocatable :: ipiv
        real(rp), dimension(:,:), allocatable :: TMatLU

        call mallocate('cp_inversion_solver [TMatLU]', n, n, TMatLU)
        call mallocate('cp_inversion_solver [ipiv]', n, ipiv)

        TMatLU = tmat
        call dgetrf(n, n, TMatLU, n, ipiv, info)

        x = rhs
        call dgetrs('N', n, nrhs, TMatLU, n, ipiv, x, n, info)

        call mfree('cp_inversion_solver [TMatLU]', TMatLU)
        call mfree('cp_inversion_solver [ipiv]', ipiv)

    end subroutine cp_inversion_solver

    subroutine conjugate_gradient_solver(n, rhs, x, eel, matvec, precnd, &
                                         arg_tol, arg_n_iter)
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
        type(ommp_electrostatics_type), intent(in) :: eel
        !! Electrostatics data structure
        external :: matvec
        !! Routine to perform matrix-vector product
        external :: precnd
        !! Preconditioner routine

        integer(ip) :: it
        real(rp) :: rms_norm, alpha, gnew, gold, gama
        real(rp), allocatable :: r(:), p(:), h(:), z(:)
        character(len=OMMP_STR_CHAR_MAX) :: msg

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

        call ommp_message("Solving linear system with CG solver", OMMP_VERBOSE_LOW)
        write(msg, "(A, I4)") "Max iter:", n_iter
        call ommp_message(msg, OMMP_VERBOSE_LOW)
        write(msg, "(A, E8.1)") "Tolerance: ", tol
        call ommp_message(msg, OMMP_VERBOSE_LOW)

        call mallocate('conjugate_gradient_solver [r]', n, r)
        call mallocate('conjugate_gradient_solver [p]', n, p)
        call mallocate('conjugate_gradient_solver [h]', n, h)
        call mallocate('conjugate_gradient_solver [z]', n, z)

        ! compute a guess, if required:
        rms_norm = dot_product(x,x)
        if(rms_norm < eps_rp) then
            call ommp_message("Input guess has zero norm, generating a guess&
                              & from preconditioner.", OMMP_VERBOSE_HIGH)
            call precnd(eel, x, x)
        else
            call ommp_message("Using input guess as a starting point for&
                              & iterative solver.", OMMP_VERBOSE_HIGH)
        end if

        ! compute the residual:
        call matvec(eel, x, z, .true.)
        r = rhs - z
        ! apply the preconditioner and get the first direction:
        call precnd(eel, r, z)
        p = z
        gold = dot_product(r, z)
        gama = 0.0_rp

        do it = 1, n_iter
            ! compute the step:
            call matvec(eel, p, h, .true.)
            gama = dot_product(h, p)

            ! unlikely quick return:
            if(abs(gama) < eps_rp) then
                call ommp_message("Direction vector with zero norm, exiting &
                                  &iterative solver.", OMMP_VERBOSE_HIGH)
                exit
            end if

            alpha = gold / gama
            x = x + alpha * p
            r = r - alpha * h

            ! apply the preconditioner:
            call precnd(eel, r, z)
            gnew = dot_product(r, z)
            rms_norm = sqrt(gnew/dble(n))

            write(msg, "('iter=',i4,' residual rms norm: ', d14.4)") it, rms_norm
            call ommp_message(msg, OMMP_VERBOSE_HIGH)

            ! Check convergence
            if(rms_norm < tol) then
                call ommp_message("Required convergence threshold reached, &
                                  &exiting iterative solver.", OMMP_VERBOSE_HIGH)
                exit
            end if

            ! compute the next direction:
            gama = gnew/gold
            p    = gama*p + z
            gold = gnew
        end do

        call mfree('conjugate_gradient_solver [r]', r)
        call mfree('conjugate_gradient_solver [p]', p)
        call mfree('conjugate_gradient_solver [h]', h)
        call mfree('conjugate_gradient_solver [z]', z)

        if(rms_norm > tol .and. abs(gama) > eps_rp) then
            call fatal_error("Iterative solver did not converged")
        end if

    end subroutine conjugate_gradient_solver

    subroutine block_conjugate_gradient_solver(n, nrhs, rhs, x, eel, matvec, precnd, &
                                               arg_tol, arg_n_iter, arg_matvec_cols)
        !! Block conjugate gradient solver: solves A X = RHS for X, where X
        !! and RHS are (n, nrhs) matrices, simultaneously for all nrhs
        !! columns. Standard block-CG recursion (e.g. O'Leary 1980): the
        !! nrhs x nrhs "step size" and "conjugation" coefficients (alpha,
        !! beta below) are found by solving small nrhs x nrhs linear systems
        !! each iteration (via LAPACK dgesv -- nrhs is expected to stay
        !! small, e.g. 2 for AMOEBA's D/P dipoles, so this is not a
        !! bottleneck), rather than being scalars as in plain CG.
        !!
        !! The actual performance case for using this over nrhs independent
        !! calls to conjugate_gradient_solver is that matvec/precnd are
        !! themselves block routines that can share work (e.g. an FMM
        !! far-field tree pass) across all nrhs columns of one matrix-vector
        !! product -- if matvec/precnd are just naive per-column loops
        !! around a single-vector routine, block CG mainly buys whatever
        !! convergence benefit comes from the shared/larger Krylov subspace,
        !! and pays a bit extra for the small linear solves each iteration.

        use mod_constants, only: eps_rp
        use mod_memory, only: mallocate, mfree

        implicit none

        integer(ip), intent(in) :: n
        !! Size of each column of the (block) linear system
        integer(ip), intent(in) :: nrhs
        !! Number of right-hand-sides solved simultaneously
        real(rp), intent(in), optional :: arg_tol
        real(rp) :: tol
        !! Convergence criterion: max over columns of RMS norm of residual < tol

        integer(ip), intent(in), optional :: arg_n_iter
        integer(ip) :: n_iter

        real(rp), dimension(n, nrhs), intent(in) :: rhs
        !! Right hand sides of the linear system, one per column
        real(rp), dimension(n, nrhs), intent(inout) :: x
        !! In input, initial guess for the solver (one column per rhs), in
        !! output the solutions
        type(ommp_electrostatics_type), intent(in) :: eel
        !! Electrostatics data structure
        external :: matvec
        !! Block matrix-vector routine: matvec(eel, nrhs, x, y, dodiag),
        !! x/y shape (n, nrhs)
        external :: precnd
        !! Block preconditioner routine: precnd(eel, nrhs, x, y), x/y shape
        !! (n, nrhs)
        integer(ip), intent(inout), optional :: arg_matvec_cols
        !! If present, incremented by nrhs for every matvec call issued
        !! (i.e. accumulates total column-matvecs) -- caller must
        !! initialize to 0 before the first call to get a meaningful total
        !! across repeated/batched invocations. Diagnostic only, used to
        !! compare solver strategies.

        integer(ip) :: it, i, info
        real(rp) :: rms_norm
        real(rp), allocatable :: r(:,:), p(:,:), q(:,:), z(:,:)
        real(rp), allocatable :: gold(:,:), gnew(:,:), h(:,:), alpha(:,:), beta(:,:)
        integer(ip), allocatable :: ipiv(:)
        character(len=OMMP_STR_CHAR_MAX) :: msg

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

        write(msg, "(A, I4, A, I4)") "Solving block linear system with CG solver, nrhs=", &
                                     nrhs, ", max iter:", n_iter
        call ommp_message(msg, OMMP_VERBOSE_LOW)
        write(msg, "(A, E8.1)") "Tolerance: ", tol
        call ommp_message(msg, OMMP_VERBOSE_LOW)

        call mallocate('block_conjugate_gradient_solver [r]', n, nrhs, r)
        call mallocate('block_conjugate_gradient_solver [p]', n, nrhs, p)
        call mallocate('block_conjugate_gradient_solver [q]', n, nrhs, q)
        call mallocate('block_conjugate_gradient_solver [z]', n, nrhs, z)
        call mallocate('block_conjugate_gradient_solver [gold]', nrhs, nrhs, gold)
        call mallocate('block_conjugate_gradient_solver [gnew]', nrhs, nrhs, gnew)
        call mallocate('block_conjugate_gradient_solver [h]', nrhs, nrhs, h)
        call mallocate('block_conjugate_gradient_solver [alpha]', nrhs, nrhs, alpha)
        call mallocate('block_conjugate_gradient_solver [beta]', nrhs, nrhs, beta)
        call mallocate('block_conjugate_gradient_solver [ipiv]', nrhs, ipiv)

        ! compute the residual for the (given) initial guess:
        call matvec(eel, nrhs, x, z, .true.)
        if(present(arg_matvec_cols)) arg_matvec_cols = arg_matvec_cols + nrhs
        r = rhs - z
        ! apply the preconditioner and get the first direction:
        call precnd(eel, nrhs, r, z)
        p = z
        gold = matmul(transpose(r), z)

        rms_norm = huge(1.0_rp)
        do it = 1, n_iter
            ! compute the step:
            call matvec(eel, nrhs, p, q, .true.)
            if(present(arg_matvec_cols)) arg_matvec_cols = arg_matvec_cols + nrhs
            h = matmul(transpose(p), q)

            ! alpha solves h * alpha = gold
            alpha = gold
            call dgesv(nrhs, nrhs, h, nrhs, ipiv, alpha, nrhs, info)
            if(info /= 0) then
                call ommp_message("Block Gram matrix is singular, exiting &
                                  &iterative solver.", OMMP_VERBOSE_HIGH)
                exit
            end if

            x = x + matmul(p, alpha)
            r = r - matmul(q, alpha)

            ! apply the preconditioner:
            call precnd(eel, nrhs, r, z)
            gnew = matmul(transpose(r), z)

            ! Convergence check: worst column RMS norm
            rms_norm = 0.0_rp
            do i = 1, nrhs
                rms_norm = max(rms_norm, sqrt(sum(r(:,i)**2) / dble(n)))
            end do

            write(msg, "('iter=',i4,' worst-column residual rms norm: ', d14.4)") it, rms_norm
            call ommp_message(msg, OMMP_VERBOSE_HIGH)

            if(rms_norm < tol) then
                call ommp_message("Required convergence threshold reached, &
                                  &exiting iterative solver.", OMMP_VERBOSE_HIGH)
                exit
            end if

            ! beta solves gold * beta = gnew (gold is the OLD gram matrix,
            ! used here as the matrix to invert -- NOT h)
            beta = gnew
            call dgesv(nrhs, nrhs, gold, nrhs, ipiv, beta, nrhs, info)
            if(info /= 0) then
                call ommp_message("Block Gram matrix is singular, exiting &
                                  &iterative solver.", OMMP_VERBOSE_HIGH)
                exit
            end if

            p = z + matmul(p, beta)
            gold = gnew
        end do

        call mfree('block_conjugate_gradient_solver [r]', r)
        call mfree('block_conjugate_gradient_solver [p]', p)
        call mfree('block_conjugate_gradient_solver [q]', q)
        call mfree('block_conjugate_gradient_solver [z]', z)
        call mfree('block_conjugate_gradient_solver [gold]', gold)
        call mfree('block_conjugate_gradient_solver [gnew]', gnew)
        call mfree('block_conjugate_gradient_solver [h]', h)
        call mfree('block_conjugate_gradient_solver [alpha]', alpha)
        call mfree('block_conjugate_gradient_solver [beta]', beta)
        call mfree('block_conjugate_gradient_solver [ipiv]', ipiv)

        if(rms_norm > tol) then
            call fatal_error("Block iterative solver did not converge")
        end if

    end subroutine block_conjugate_gradient_solver

    subroutine batched_block_cg_solver(n, nrhs, rhs, x, eel, matvec, precnd, &
                                       arg_tol, arg_n_iter, arg_batch_size, arg_matvec_cols)
        !! Thin wrapper around block_conjugate_gradient_solver that splits a
        !! large RHS matrix into fixed-size column batches and solves each
        !! batch independently. block_conjugate_gradient_solver's block Gram
        !! matrix (nrhs x nrhs, see there) becomes numerically unsafe as
        !! nrhs approaches n (near/exactly rank-deficient, undetected by
        !! dgesv's exact-singularity check -- observed to blow up to NaN in
        !! as few as ~85 iterations on a case with nrhs==n); batching keeps
        !! each block width far below n regardless of how large the caller's
        !! total nrhs is (e.g. CPID's nrhs=3*mm_atoms, which is always >=
        !! n=3*pol_atoms).

        implicit none

        integer(ip), intent(in) :: n, nrhs
        real(rp), intent(in), optional :: arg_tol
        integer(ip), intent(in), optional :: arg_n_iter, arg_batch_size
        real(rp), dimension(n, nrhs), intent(in) :: rhs
        real(rp), dimension(n, nrhs), intent(inout) :: x
        type(ommp_electrostatics_type), intent(in) :: eel
        external :: matvec, precnd
        integer(ip), intent(inout), optional :: arg_matvec_cols

        integer(ip) :: batch_size, i0, i1

        if(present(arg_batch_size)) then
            batch_size = arg_batch_size
        else
            batch_size = max(1_ip, nint(0.05_rp*real(nrhs, rp)))
        end if

        i0 = 1
        do while(i0 <= nrhs)
            i1 = min(i0 + batch_size - 1, nrhs)
            call block_conjugate_gradient_solver(n, i1-i0+1, rhs(:,i0:i1), x(:,i0:i1), &
                                                 eel, matvec, precnd, arg_tol, arg_n_iter, &
                                                 arg_matvec_cols)
            i0 = i1 + 1
        end do
    end subroutine batched_block_cg_solver

    subroutine deflated_subspace_cg_solver(n, nrhs, rhs, x, eel, matvec, precnd, &
                                           arg_tol, arg_n_iter, arg_batch_size, &
                                           arg_defl_tol, arg_matvec_cols, &
                                           arg_final_subspace_size, &
                                           arg_max_subspace_size)
        !! Deflated/recycled subspace solver for many RHS against one shared
        !! SPD operator A, for the regime where nrhs is comparable to (or
        !! larger than) n -- exactly CPID's case (nrhs=3*mm_atoms >=
        !! n=3*pol_atoms always). Maintains ONE growing orthonormal basis V
        !! (with cached AV and the small projected Gram matrix
        !! G=(AV)^T(AV)) across the whole RHS sweep, instead of independent
        !! Krylov spaces per batch. RHS columns are processed in batches;
        !! each batch is first solved by projecting onto the CURRENT V (a
        !! small dense solve, zero new matvecs); only the resulting
        !! preconditioned residual directions that are NOT already
        !! representable in V (checked by Gram-Schmidt against V, dropping
        !! near-zero results -- "deflation") turn into new matvecs and grow
        !! V. If many RHS columns' true response lives in a shared
        !! low-dimensional subspace (plausible here: nearby geometric
        !! perturbations should excite very similar induced-dipole response
        !! patterns), V plateaus well below n and later batches converge
        !! for free.
        !!
        !! The projection minimizes the residual directly (solve
        !! (AV)^T(AV) y = (AV)^T b, i.e. the normal equations of
        !! min_y ||b - AVy||_2, GMRES/MINRES-style) rather than using a
        !! Galerkin/FOM-style V^T A V y = V^T b projection: the latter only
        !! guarantees a monotonic decrease of the A-norm error, NOT of the
        !! plain residual norm this solver's convergence check (and
        !! block_conjugate_gradient_solver's, for a fair comparison) both
        !! use -- confirmed by hitting exactly this pathology during
        !! development (residual would shrink for a while then drift back
        !! up over hundreds of iterations even as V kept growing with
        !! seemingly-new orthogonal directions). Minimizing the residual
        !! directly guarantees it can only shrink (in exact arithmetic) as
        !! V grows, regardless of how many columns are being solved at once.
        !!
        !! Unlike block_conjugate_gradient_solver, this also stays
        !! numerically safe even if the deflation assumption above fails
        !! and V has to grow all the way to n: V is orthonormal by
        !! construction, so G's conditioning tracks A's own squared
        !! conditioning (in an orthonormal basis) instead of the raw
        !! block-CG Gram matrix, which was observed to become
        !! catastrophically (and silently, past dgesv's exact-singularity
        !! check) ill-conditioned exactly in the nrhs>=n regime CPID always
        !! operates in.
        !!
        !! The projected Gram matrix G=(AV)^T(AV) is maintained as a
        !! Cholesky factor L (G=L L^T) that is extended INCREMENTALLY
        !! (bordering: solve L11 y = G12, then Cholesky-factor the small
        !! Schur complement G22 - y^T y) every time the subspace grows,
        !! instead of being refactorized from scratch every iteration.
        !! Refactorizing the whole (growing) G on every iteration was the
        !! actual, avoidable, dominant cost during development (an
        !! O(mcur^3) LU factorization repeated at every single iteration,
        !! not anything intrinsic to the problem's size) -- bordering turns
        !! the total Cholesky cost across the whole sweep into
        !! O(mcur_final^3), i.e. one dense-equivalent factorization, same
        !! order as the cost Filippo expected from dense linear algebra at
        !! this problem size.
        !!
        !! The subspace is capped at arg_max_subspace_size (default 10% of
        !! min(n,nrhs)): if the deflation assumption fails badly enough
        !! that V would need to grow past that, this method's overhead no
        !! longer pays for itself, so the remainder of the sweep (this
        !! batch's still-unconverged columns, plus every later batch) is
        !! handed off to batched_block_cg_solver instead -- warm-started by
        !! projecting onto the frozen V first, so whatever this method
        !! already learned is not thrown away.

        use mod_memory, only: mallocate, mfree

        implicit none

        integer(ip), intent(in) :: n, nrhs
        real(rp), dimension(n, nrhs), intent(in) :: rhs
        real(rp), dimension(n, nrhs), intent(inout) :: x
        type(ommp_electrostatics_type), intent(in) :: eel
        external :: matvec
        !! Block matrix-vector routine: matvec(eel, nrhs, x, y, dodiag)
        external :: precnd
        !! Block preconditioner routine: precnd(eel, nrhs, x, y)
        real(rp), intent(in), optional :: arg_tol, arg_defl_tol
        integer(ip), intent(in), optional :: arg_n_iter, arg_batch_size, arg_max_subspace_size
        integer(ip), intent(inout), optional :: arg_matvec_cols
        integer(ip), intent(out), optional :: arg_final_subspace_size

        real(rp) :: tol, defl_tol
        integer(ip) :: n_iter, batch_size, max_subspace_size

        real(rp), allocatable :: V(:,:), AV(:,:), G(:,:), L(:,:)
        integer(ip) :: mcur, mcap

        real(rp), allocatable :: B(:,:), Xb(:,:), R(:,:), Z(:,:)
        real(rp), allocatable :: y(:,:), rhs_proj(:,:), Ycol(:,:)
        real(rp), allocatable :: xu(:,:), Zproj(:,:)
        integer(ip), allocatable :: keep_idx(:)
        real(rp) :: rms, nrm, nrm0
        integer(ip) :: i0, i1, bs, it, i, j, info, mnew
        logical :: accept_stagnation
        character(len=OMMP_STR_CHAR_MAX) :: msg

        integer(ip) :: clock_rate, clock_t0, clock_t1
        real(rp) :: t_matvec, t_precnd, t_solve, t_cholesky, t_screen, t_ortho, t_fallback, t_total0, t_total1

        call system_clock(count_rate=clock_rate)
        call system_clock(clock_t0)
        t_matvec = 0.0_rp; t_precnd = 0.0_rp; t_solve = 0.0_rp
        t_cholesky = 0.0_rp; t_screen = 0.0_rp; t_ortho = 0.0_rp; t_fallback = 0.0_rp
        t_total0 = real(clock_t0, rp) / real(clock_rate, rp)

        if(present(arg_tol)) then; tol = arg_tol; else; tol = OMMP_DEFAULT_SOLVER_TOL; end if
        if(present(arg_n_iter)) then; n_iter = arg_n_iter; else; n_iter = OMMP_DEFAULT_SOLVER_ITER; end if
        if(present(arg_batch_size)) then
            batch_size = arg_batch_size
        else
            batch_size = max(1_ip, nint(0.05_rp*real(nrhs, rp)))
        end if
        ! Deflation tolerance: ABSOLUTE bound on the leftover (not-yet-in-V)
        ! norm of a candidate direction, defaulting to 1e-4 * tol. Deriving
        ! it from tol (rather than a fixed machine-precision-scale
        ! constant, or comparing to the candidate's OWN norm) matters once
        ! the batch's residual has shrunk close to tol: the earlier
        ! nrm0-relative version compared two quantities that both shrink
        ! together as convergence is approached, which is exactly where
        ! Z - V*(V^T*Z)-style cancellation is least trustworthy -- an
        ! absolute floor tied to tol instead stays comfortably (4 orders
        ! of magnitude) above double-precision roundoff while still being
        ! far enough below tol to keep any direction that could plausibly
        ! still matter for reaching it.
        if(present(arg_defl_tol)) then; defl_tol = arg_defl_tol; else; defl_tol = tol * 1.0e-4_rp; end if
        if(present(arg_max_subspace_size)) then
            max_subspace_size = arg_max_subspace_size
        else
            max_subspace_size = max(batch_size, nint(0.10_rp*real(min(n, nrhs), rp)))
        end if
        ! Bound against n only (the subspace lives in R^n, so it can never
        ! usefully exceed that) -- NOT against nrhs: nrhs is just how many
        ! columns THIS call happens to solve, which matters for the
        ! DEFAULT heuristic above (no way to know the "true" full-scale
        ! problem size otherwise) but must not silently override an
        ! explicit arg_max_subspace_size when testing on a smaller subset
        ! of a larger production nrhs.
        max_subspace_size = min(max_subspace_size, n)

        mcur = 0
        mcap = min(max_subspace_size, max(4_ip, 4*batch_size))
        call mallocate('deflated_subspace_cg_solver [V]', n, mcap, V)
        call mallocate('deflated_subspace_cg_solver [AV]', n, mcap, AV)
        call mallocate('deflated_subspace_cg_solver [G]', mcap, mcap, G)
        call mallocate('deflated_subspace_cg_solver [L]', mcap, mcap, L)

        i0 = 1
        do while(i0 <= nrhs)
            i1 = min(i0 + batch_size - 1, nrhs)
            bs = i1 - i0 + 1

            call mallocate('deflated_subspace_cg_solver [B]', n, bs, B)
            call mallocate('deflated_subspace_cg_solver [Xb]', n, bs, Xb)
            call mallocate('deflated_subspace_cg_solver [R]', n, bs, R)
            B = rhs(:, i0:i1)
            Xb = x(:, i0:i1)

            rms = huge(1.0_rp)
            accept_stagnation = .false.
            do it = 1, n_iter
                if(mcur > 0) then
                    call system_clock(clock_t0)
                    call mallocate('deflated_subspace_cg_solver [rhs_proj]', mcur, bs, rhs_proj)
                    call mallocate('deflated_subspace_cg_solver [y]', mcur, bs, y)

                    ! Minimum-residual projection: y solves the normal
                    ! equations of min_y ||B - AV y||_2, NOT the Galerkin
                    ! system V^T A V y = V^T B (see header note -- Galerkin
                    ! does not guarantee a monotonic plain residual norm).
                    ! Solved via the incrementally-maintained Cholesky
                    ! factor L of G (back-substitution only, O(mcur^2*bs) --
                    ! see the subroutine header for why this replaced a
                    ! from-scratch dgesv factorization here).
                    rhs_proj = matmul(transpose(AV(:,1:mcur)), B)
                    y = rhs_proj
                    call dpotrs('L', mcur, bs, L, mcap, y, mcur, info)
                    if(info /= 0) call fatal_error("deflated_subspace_cg_solver: back-&
                                                    &substitution against the projected &
                                                    &Gram matrix's Cholesky factor failed &
                                                    &(should not happen for an orthonormal &
                                                    &basis).")
                    Xb = matmul(V(:,1:mcur), y)
                    R = B - matmul(AV(:,1:mcur), y)

                    call mfree('deflated_subspace_cg_solver [rhs_proj]', rhs_proj)
                    call mfree('deflated_subspace_cg_solver [y]', y)
                    call system_clock(clock_t1)
                    t_solve = t_solve + real(clock_t1-clock_t0, rp)/real(clock_rate, rp)
                else
                    Xb = 0.0_rp
                    R = B
                end if

                rms = 0.0_rp
                do i = 1, bs
                    rms = max(rms, sqrt(sum(R(:,i)**2) / dble(n)))
                end do
                write(msg, "('deflated subspace solver: batch [',i6,':',i6,'] it=',i4,&
                             &' subspace size=',i6,' worst-column residual rms norm: ',d14.4)") &
                             i0, i1, it, mcur, rms
                call ommp_message(msg, OMMP_VERBOSE_HIGH)
                if(rms < tol) exit

                ! Precondition the residual. Then, decide which columns are
                ! worth turning into new basis vectors via a cheap BLAS3
                ! screen (project against the CURRENT V, compare the
                ! resulting residual norm to the original column norm) --
                ! this is only a SCREEN for the deflation decision. The
                ! actual orthogonalization/orthonormalization of whichever
                ! columns survive is done afterwards by ortho_vs_x/ortho_cd
                ! (Cholesky-based, BLAS3, iteratively refined until tight
                ! orthogonality is reached -- see their header comments),
                ! which replaces the previous hand-rolled two-pass
                ! Gram-Schmidt: that scheme did exactly two passes
                ! regardless of whether that achieved tight orthogonality,
                ! which matters most exactly in the regime we operate in
                ! here -- candidates near the deflation threshold are, by
                ! definition, close to parallel to existing directions, the
                ! textbook case where naive Gram-Schmidt is least reliable.
                call system_clock(clock_t0)
                call mallocate('deflated_subspace_cg_solver [Z]', n, bs, Z)
                call precnd(eel, bs, R, Z)
                call system_clock(clock_t1)
                t_precnd = t_precnd + real(clock_t1-clock_t0, rp)/real(clock_rate, rp)

                call system_clock(clock_t0)
                call mallocate('deflated_subspace_cg_solver [Zproj]', n, bs, Zproj)
                if(mcur > 0) then
                    call mallocate('deflated_subspace_cg_solver [xu]', mcur, bs, xu)
                    xu = matmul(transpose(V(:,1:mcur)), Z)
                    Zproj = Z - matmul(V(:,1:mcur), xu)
                    ! second pass (classical Gram-Schmidt reorthogonalization):
                    ! a single projection pass can leave a spuriously large
                    ! "leftover" norm for a vector that is actually already
                    ! well represented in V (Kahan's classical cancellation
                    ! issue with one-pass GS), which would fool the deflation
                    ! threshold below into keeping near-duplicate directions.
                    ! A second pass against the same V recovers the true
                    ! residual-vs-V norm to good accuracy -- this is exactly
                    ! why the routine this replaced always did two passes for
                    ! its norm estimate too.
                    xu = matmul(transpose(V(:,1:mcur)), Zproj)
                    Zproj = Zproj - matmul(V(:,1:mcur), xu)
                    call mfree('deflated_subspace_cg_solver [xu]', xu)
                else
                    Zproj = Z
                end if
                call system_clock(clock_t1)
                t_screen = t_screen + real(clock_t1-clock_t0, rp)/real(clock_rate, rp)

                call mallocate('deflated_subspace_cg_solver [keep_idx]', bs, keep_idx)
                mnew = 0
                do i = 1, bs
                    nrm0 = norm2(Z(:,i))
                    if(nrm0 <= 0.0_rp) cycle
                    nrm = norm2(Zproj(:,i))

                    if(nrm > defl_tol) then
                        mnew = mnew + 1
                        keep_idx(mnew) = i
                    end if
                end do
                call mfree('deflated_subspace_cg_solver [Zproj]', Zproj)

                ! Deflation-exhaustion acceptance: if NO candidate's
                ! leftover-vs-V component exceeds defl_tol (default
                ! tol*1e-4, i.e. 4 orders of magnitude below the target
                ! residual tolerance), there is genuinely nothing left, in
                ! a double-precision-safe sense, that this batch's current
                ! best projection onto V is missing -- accept it as
                ! converged even though the raw residual may sit just
                ! above tol, rather than grinding out the last bit with a
                ! block CG fallback. This is a deliberate accuracy/cost
                ! trade-off (Filippo's proposal): forcing exact convergence
                ! past this point would mean resolving directions the
                ! deflation screen cannot distinguish from noise at
                ! defl_tol precision.
                accept_stagnation = (mnew == 0)

                ! Respect the subspace cap: never grow V past
                ! max_subspace_size (see header). Unlike genuine
                ! deflation-exhaustion above, running out of ROOM for
                ! otherwise-genuine new directions is a real resource
                ! limit, not evidence there is nothing left to find -- it
                ! still falls back to block CG below (accept_stagnation
                ! stays false in that case).
                mnew = min(mnew, max(0_ip, max_subspace_size - mcur))

                if(mnew == 0) then
                    call mfree('deflated_subspace_cg_solver [Z]', Z)
                    call mfree('deflated_subspace_cg_solver [keep_idx]', keep_idx)
                    exit
                end if

                if(mcur+mnew > mcap) call grow_subspace_(n, mcap, mcur+mnew, V, AV, G, L)
                do j = 1, mnew
                    V(:,mcur+j) = Z(:,keep_idx(j))
                end do

                call mfree('deflated_subspace_cg_solver [Z]', Z)
                call mfree('deflated_subspace_cg_solver [keep_idx]', keep_idx)

                ! Robust orthogonalization of the new block against V, and
                ! orthonormalization among the new columns themselves
                ! (handles near-duplicate directions arising within a
                ! single batch too, not just against pre-existing V).
                call system_clock(clock_t0)
                call ortho_vs_x(n, mcur, mnew, V(:,1:mcur), V(:,mcur+1:mcur+mnew))
                call system_clock(clock_t1)
                t_ortho = t_ortho + real(clock_t1-clock_t0, rp)/real(clock_rate, rp)

                call system_clock(clock_t0)
                call matvec(eel, mnew, V(:,mcur+1:mcur+mnew), AV(:,mcur+1:mcur+mnew), .true.)
                if(present(arg_matvec_cols)) arg_matvec_cols = arg_matvec_cols + mnew
                call system_clock(clock_t1)
                t_matvec = t_matvec + real(clock_t1-clock_t0, rp)/real(clock_rate, rp)

                call system_clock(clock_t0)
                if(mcur > 0) then
                    G(1:mcur, mcur+1:mcur+mnew) = matmul(transpose(AV(:,1:mcur)), AV(:,mcur+1:mcur+mnew))
                    G(mcur+1:mcur+mnew, 1:mcur) = transpose(G(1:mcur, mcur+1:mcur+mnew))
                end if
                G(mcur+1:mcur+mnew, mcur+1:mcur+mnew) = &
                    matmul(transpose(AV(:,mcur+1:mcur+mnew)), AV(:,mcur+1:mcur+mnew))

                ! Extend the Cholesky factor L incrementally (bordering)
                ! instead of refactorizing G(1:mcur+mnew,1:mcur+mnew) from
                ! scratch (see subroutine header).
                if(mcur > 0) then
                    call mallocate('deflated_subspace_cg_solver [Ycol]', mcur, mnew, Ycol)
                    Ycol = G(1:mcur, mcur+1:mcur+mnew)
                    ! Y = L11^-1 G12 (forward substitution, L11 lower triangular)
                    call dtrsm('L', 'L', 'N', 'N', mcur, mnew, 1.0_rp, L, mcap, Ycol, mcur)
                    L(mcur+1:mcur+mnew, 1:mcur) = transpose(Ycol)
                    ! Schur complement: S = G22 - L21 L21^T = G22 - Y^T Y
                    L(mcur+1:mcur+mnew, mcur+1:mcur+mnew) = &
                        G(mcur+1:mcur+mnew, mcur+1:mcur+mnew) - matmul(transpose(Ycol), Ycol)
                    call mfree('deflated_subspace_cg_solver [Ycol]', Ycol)
                else
                    L(1:mnew, 1:mnew) = G(1:mnew, 1:mnew)
                end if
                call dpotrf('L', mnew, L(mcur+1,mcur+1), mcap, info)
                if(info /= 0) then
                    ! Extremely rare fallback: the Schur complement was not
                    ! numerically PD (would indicate a near-duplicate
                    ! direction slipped past deflation). Recompute L from
                    ! scratch for this one growth event instead of
                    ! aborting -- a safety net, not the common path, so its
                    ! O(mcur^3) cost here is acceptable.
                    L(1:mcur+mnew,1:mcur+mnew) = G(1:mcur+mnew,1:mcur+mnew)
                    call dpotrf('L', mcur+mnew, L, mcap, info)
                    if(info /= 0) call fatal_error("deflated_subspace_cg_solver: Cholesky &
                                                    &factorization of the projected Gram &
                                                    &matrix failed (should not happen for &
                                                    &an orthonormal basis).")
                end if
                call system_clock(clock_t1)
                t_cholesky = t_cholesky + real(clock_t1-clock_t0, rp)/real(clock_rate, rp)

                mcur = mcur + mnew
            end do

            x(:, i0:i1) = Xb

            if(rms > tol .and. accept_stagnation) then
                ! Deflation-exhaustion acceptance (see above): no candidate
                ! carried a leftover component above defl_tol, so the
                ! current Xb is accepted as converged even though the raw
                ! residual sits (slightly) above tol -- report it, but do
                ! not fall back.
                write(msg, "('deflated subspace solver: batch [',i6,':',i6,'] &
                             &accepted via deflation exhaustion (defl_tol=',&
                             &d10.2,') at residual rms norm=',d14.4,' (tol=',&
                             &d10.2,').')") i0, i1, defl_tol, rms, tol
                call ommp_message(msg, OMMP_VERBOSE_LOW)
            else if(rms > tol) then
                ! This batch did not converge via the subspace method alone
                ! and was NOT accepted via deflation exhaustion (i.e. the
                ! max_subspace_size cap was reached while genuine new
                ! information was still available -- see above). Finish
                ! just THIS batch with a plain, already-numerically-safe
                ! block CG call (warm-started from the best Xb found so
                ! far), then resume attempting deflation on subsequent
                ! batches: one batch running out of subspace room does not
                ! mean later, potentially quite different, RHS columns
                ! won't deflate against V just fine (though if the cap is
                ! genuinely saturated, every later batch will fall back the
                ! same way, one at a time, which is the expected behavior).
                write(msg, "('deflated subspace solver: subspace cap on batch &
                             &[',i6,':',i6,'], finishing it with block CG.')") i0, i1
                call ommp_message(msg, OMMP_VERBOSE_LOW)

                call system_clock(clock_t0)
                call block_conjugate_gradient_solver(n, bs, rhs(:,i0:i1), x(:,i0:i1), &
                                                     eel, matvec, precnd, tol, n_iter, &
                                                     arg_matvec_cols)
                call system_clock(clock_t1)
                t_fallback = t_fallback + real(clock_t1-clock_t0, rp)/real(clock_rate, rp)
            end if

            call mfree('deflated_subspace_cg_solver [B]', B)
            call mfree('deflated_subspace_cg_solver [Xb]', Xb)
            call mfree('deflated_subspace_cg_solver [R]', R)

            i0 = i1 + 1
        end do

        if(present(arg_final_subspace_size)) arg_final_subspace_size = mcur

        call mfree('deflated_subspace_cg_solver [V]', V)
        call mfree('deflated_subspace_cg_solver [AV]', AV)
        call mfree('deflated_subspace_cg_solver [G]', G)
        call mfree('deflated_subspace_cg_solver [L]', L)

        call system_clock(clock_t1)
        t_total1 = real(clock_t1, rp) / real(clock_rate, rp)
        write(msg, "('deflated subspace solver phase timings (s): matvec=',f10.3,&
                     &' precnd=',f10.3,' solve(dpotrs)=',f10.3,' cholesky-update=',&
                     &f10.3,' defl-screen=',f10.3,' ortho=',f10.3,' fallback(block CG)=',&
                     &f10.3,' total=',f10.3)") &
                     t_matvec, t_precnd, t_solve, t_cholesky, t_screen, t_ortho, t_fallback, &
                     t_total1-t_total0
        call ommp_message(msg, OMMP_VERBOSE_LOW)

    contains

        subroutine grow_subspace_(n, mcap, need, V, AV, G, L)
            integer(ip), intent(in) :: n, need
            integer(ip), intent(inout) :: mcap
            real(rp), allocatable, intent(inout) :: V(:,:), AV(:,:), G(:,:), L(:,:)
            real(rp), allocatable :: Vnew(:,:), AVnew(:,:), Gnew(:,:), Lnew(:,:)
            integer(ip) :: newcap

            newcap = min(max(2*mcap, need), max_subspace_size)
            call mallocate('grow_subspace_ [Vnew]', n, newcap, Vnew)
            call mallocate('grow_subspace_ [AVnew]', n, newcap, AVnew)
            call mallocate('grow_subspace_ [Gnew]', newcap, newcap, Gnew)
            call mallocate('grow_subspace_ [Lnew]', newcap, newcap, Lnew)
            Vnew = 0.0_rp; AVnew = 0.0_rp; Gnew = 0.0_rp; Lnew = 0.0_rp
            Vnew(:,1:mcap) = V
            AVnew(:,1:mcap) = AV
            Gnew(1:mcap,1:mcap) = G
            Lnew(1:mcap,1:mcap) = L
            call mfree('grow_subspace_ [V]', V)
            call mfree('grow_subspace_ [AV]', AV)
            call mfree('grow_subspace_ [G]', G)
            call mfree('grow_subspace_ [L]', L)
            call move_alloc(Vnew, V)
            call move_alloc(AVnew, AV)
            call move_alloc(Gnew, G)
            call move_alloc(Lnew, L)
            mcap = newcap
        end subroutine grow_subspace_

    end subroutine deflated_subspace_cg_solver

    subroutine ortho_vs_x(n, m, k, x, u)
        !! Given two sets x(n,m) and u(n,k) of vectors, where x is assumed
        !! to be orthogonal (to machine precision), orthogonalize u against
        !! x, and orthonormalize u among its own columns. The u vs x
        !! orthogonalization and the subsequent orthonormalization of u are
        !! iterated until the overlap between x and the orthogonalized u is
        !! smaller than a (tight) threshold -- unlike a fixed two-pass
        !! Gram-Schmidt, this keeps refining for as long as needed (up to
        !! maxit), which matters most for candidate columns that are close
        !! to parallel with x, the case where a fixed number of passes is
        !! least reliable. Adapted from Filippo's diaglib (Davidson
        !! diagonalization toolkit): originally free-standing, ip/rp are
        !! here host-associated from mod_solvers instead of a separate
        !! `utils` module, and hard `stop` aborts are replaced with
        !! fatal_error for consistency with the rest of this module.

        implicit none

        integer(ip), intent(in) :: n, m, k
        real(rp), dimension(n,m), intent(in) :: x
        real(rp), dimension(n,k), intent(inout) :: u

        logical :: done, ok
        integer(ip) :: it
        real(rp) :: xu_norm, growth
        real(rp), allocatable :: xu(:,:)

        integer(ip), parameter :: maxit = 10
        real(rp), parameter :: tol_ortho = 2.0_rp * epsilon(1.0_rp)
        real(rp), parameter :: zero = 0.0_rp, one = 1.0_rp

        ! quick return: nothing to orthonormalize.
        if(k == 0) return

        ! start with an initial orthogonalization to improve conditioning.
        call ortho_cd(n, k, u, growth, ok)

        ! quick return: x is empty (e.g. the very first subspace-growth
        ! event, projecting against an as-yet-empty V) -- there is nothing
        ! to project u against, u is already orthonormalized by the call
        ! above, and calling dgemm below with a leading dimension of m=0
        ! is rejected by some BLAS implementations (MKL in particular
        ! validates LD>=1 unconditionally, even for operations that are a
        ! mathematical no-op because a dimension is 0).
        if(m == 0) return

        ok = .false.
        allocate(xu(m,k))
        done = .false.
        it = 0

        ! iteratively orthogonalize u against x, and then orthonormalize u.
        do while(.not. done)
            it = it + 1

            ! u = u - x (x^t u)
            call dgemm('t', 'n', m, k, n, one, x, n, u, n, zero, xu, m)
            call dgemm('n', 'n', n, k, m, -one, x, n, xu, m, one, u, n)

            ! now, orthonormalize u.
            call ortho_cd(n, k, u, growth, ok)

            ! the orthogonalization has introduced an error that makes the
            ! new vectors no longer fully orthogonal to x. assuming that u
            ! was orthogonal to x to machine precision before, we estimate
            ! the error with growth * eps, where growth is the product of
            ! the norms of all the linear transformations applied to u.
            xu_norm = growth * epsilon(one)
            done = xu_norm < tol_ortho

            if(it > maxit) call fatal_error("ortho_vs_x: catastrophic failure &
                                             &(orthogonality vs previous subspace &
                                             &did not converge in maxit iterations).")
        end do

        deallocate(xu)

    end subroutine ortho_vs_x

    subroutine ortho_cd(n, m, u, growth, ok)
        !! Orthogonalize m vectors of length n using the Cholesky
        !! factorization of their overlap: metric = U^t U = L L^t, and the
        !! orthogonal vectors are obtained by solving the triangular linear
        !! system U(ortho) L^t = U. As Cholesky is not the most stable way
        !! of orthogonalizing a set of vectors, the orthogonalization is
        !! refined iteratively, using a conservative estimate of the
        !! orthogonalization error to assess convergence. If dpotrf fails
        !! (near-singular overlap, i.e. near-linearly-dependent input
        !! vectors), increasingly large shifts are added to the metric's
        !! diagonal until it factorizes -- this is what makes the routine
        !! safe to call on candidate blocks that may be close to rank
        !! deficient, exactly the regime the deflation screen in
        !! deflated_subspace_cg_solver selects for.
        !!
        !! Returns a growth factor, used by ortho_vs_x to estimate the
        !! orthogonality error this routine introduces w.r.t. a separate
        !! previously-orthogonal set x.

        implicit none

        integer(ip), intent(in) :: n, m
        real(rp), dimension(n,m), intent(inout) :: u
        real(rp), intent(inout) :: growth
        logical, intent(inout) :: ok

        real(rp), parameter :: tol_ortho = 2.0_rp * epsilon(1.0_rp)
        real(rp), parameter :: zero = 0.0_rp, one = 1.0_rp
        real(rp), parameter :: tol_ortho_cd = tol_ortho
        integer(ip), parameter :: maxit = 10

        integer(ip) :: it, it_micro, info
        real(rp) :: error, alpha, unorm, shift
        real(rp) :: rcond, l_norm, linv_norm
        logical :: macro_done, micro_done

        real(rp), allocatable :: metric(:,:), msave(:,:)
        real(rp), external :: dnrm2

        allocate(metric(m,m), msave(m,m))
        metric = zero
        macro_done = .false.

        it = 0
        growth = one
        do while(.not. macro_done)
            it = it + 1
            if(it > maxit) then
                ! ortho_cd failed to converge in maxit macro-iterations.
                ok = .false.
                call fatal_error("ortho_cd: maximum number of iterations reached &
                                  &while orthogonalizing a candidate block.")
            end if

            call dgemm('t', 'n', m, m, n, one, u, n, u, n, zero, metric, m)
            msave = metric

            ! compute the cholesky factorization of the metric.
            call dpotrf('l', m, metric, m, info)

            if(info /= 0) then
                ! dpotrf failed: try again after level-shifting the diagonal
                ! of the metric, with larger and larger shifts, until it
                ! manages to factorize it.
                alpha = 100.0_rp
                unorm = dnrm2(n*m, u, 1)
                it_micro = 0
                micro_done = .false.

                do while(.not. micro_done)
                    it_micro = it_micro + 1
                    if(it_micro > maxit) then
                        ok = .false.
                        call fatal_error("ortho_cd: maximum number of iterations &
                                         &reached while factorizing the (shifted) &
                                         &overlap metric of a candidate block.")
                    end if

                    shift = max(epsilon(one)*alpha*unorm, tol_ortho)
                    metric = msave
                    call diag_shift(m, shift, metric)
                    call dpotrf('l', m, metric, m, info)
                    alpha = alpha * 10.0_rp
                    micro_done = (info == 0)
                end do
            end if

            ! we assume that the error on the orthogonality is of order
            ! k(l)^2 * eps, where eps is the machine precision, and the
            ! condition number k(l) is estimated as ||l|| ||l^-1|| (see
            ! norm_est for the norm used).
            !
            ! compute l^-1, using msave to store the inverse cholesky factor.
            msave = metric
            call dtrtri('l', 'n', m, msave, m, info)

            l_norm = norm_est(m, metric)
            linv_norm = norm_est(m, msave)
            rcond = l_norm * linv_norm

            ! in each iteration of ortho_cd, we apply l^-t to u, which
            ! introduces a numerical error of order ||l^-1||. this error is
            ! saved in growth and used in ortho_vs_x to check how much
            ! ortho_cd spoiled the previously computed orthogonality to x.
            growth = growth * linv_norm

            ! orthogonalize u by applying l^(-t)
            call dtrmm('r', 'l', 't', 'n', n, m, one, msave, m, u, n)

            error = epsilon(one) * rcond*rcond
            macro_done = error < tol_ortho_cd
        end do

        ok = .true.

        deallocate(metric, msave)

    end subroutine ortho_cd

    function norm_est(m, a) result(res)
        !! Cheap estimate of the norm of a lower triangular matrix a = d +
        !! o, where d = diag(a): since ||a|| <= ||d|| + ||o||, compute
        !! ||d|| as max_i |d(i)| and ||o|| as its Frobenius norm. Tight
        !! enough for the condition-number estimate in ortho_cd, and goes
        !! to 1 when a approaches the identity.

        implicit none

        integer(ip), intent(in) :: m
        real(rp), dimension(m,m), intent(in) :: a
        real(rp) :: res

        integer(ip) :: i, j
        real(rp) :: diag_norm, od_norm

        diag_norm = 0.0_rp
        do i = 1, m
            diag_norm = max(diag_norm, abs(a(i,i)))
        end do

        od_norm = 0.0_rp
        do i = 1, m
            do j = 1, i - 1
                od_norm = od_norm + a(i,j)**2
            end do
        end do
        od_norm = sqrt(od_norm)

        res = diag_norm + od_norm

    end function norm_est

    subroutine diag_shift(n, shift, a)
        !! Add shift to the diagonal elements of the matrix a.

        implicit none

        integer(ip), intent(in) :: n
        real(rp), intent(in) :: shift
        real(rp), dimension(n,n), intent(inout) :: a

        integer(ip) :: i

        do i = 1, n
            a(i,i) = a(i,i) + shift
        end do

    end subroutine diag_shift

    subroutine jacobi_diis_solver(n, rhs, x, eel, matvec, inv_diag, arg_tol, &
                                  arg_n_iter, arg_diis_max)
    
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
        
        integer(ip), intent(in), optional :: arg_diis_max
        !! Optional maximum number of points for diis extrapolation, if not present
        !! OMMP_DEFAULT_DIIS_MAX_POINTS is used.
        integer(ip) :: diis_max
        !! Maximum number of points for diis extrapolation, if zero or negative,
        !! diis extrapolation is not used.

        real(rp), dimension(n), intent(in) :: rhs
        !! Right hand side of the linear system
        real(rp), dimension(n), intent(inout) :: x
        !! In input, initial guess for the solver, in output the solution
        type(ommp_electrostatics_type), intent(in) :: eel
        !! Electrostatics data structure
        real(rp), dimension(n), intent(in) :: inv_diag
        !! Element-wise inverse of diagonal of LHS matrix
        external :: matvec
        !! Routine to perform matrix-vector product
        
        integer(ip) :: it, nmat
        real(rp) :: rms_norm, max_norm
        logical :: do_diis
        real(rp), allocatable :: x_new(:), y(:), x_diis(:,:), e_diis(:,:), bmat(:,:)
        character(len=OMMP_STR_CHAR_MAX) :: msg
        
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
        
        if(present(arg_diis_max)) then
            diis_max = arg_diis_max
        else
            diis_max = OMMP_DEFAULT_DIIS_MAX_POINTS
        end if

        do_diis =  (diis_max > 0)
        
        call ommp_message("Solving linear system with jacobi solver", OMMP_VERBOSE_LOW)
        write(msg, "(A, I4)") "Max iter:", n_iter
        call ommp_message(msg, OMMP_VERBOSE_LOW)
        write(msg, "(A, E8.1)") "Tolerance: ", tol
        call ommp_message(msg, OMMP_VERBOSE_LOW)
        if(do_diis) then
            write(msg, "(A, I4)") "DIIS is enabled with n = ", diis_max
        else
            write(msg, "(A)") "DIIS is disabled"
        endif
        call ommp_message(msg, OMMP_VERBOSE_LOW)
        
        ! Memory allocation
        call mallocate('jacobi_diis_solver [x_new]', n, x_new)
        call mallocate('jacobi_diis_solver [y]', n, y)
        if(do_diis) then
            call mallocate('jacobi_diis_solver [x_diis]', n, diis_max, x_diis)
            call mallocate('jacobi_diis_solver [e_diis]', n, diis_max, e_diis)
            call mallocate('jacobi_diis_solver [bmat]', diis_max+1, diis_max+1, bmat)
            nmat = 1
        endif
        
        ! if required, compute a guess
        rms_norm = dot_product(x, x)
        if(rms_norm < eps_rp) then
            call ommp_message("Input guess has zero norm, generating a guess&
                              & from preconditioner.", OMMP_VERBOSE_HIGH)
            x = inv_diag * rhs
        else
            call ommp_message("Using input guess as a starting point for&
                              & iterative solver.", OMMP_VERBOSE_HIGH)
        end if
        
        ! Jacobi iterations
        do it = 1, n_iter
            ! y = rhs - O x
            call matvec(eel, x, y, .false.)
            y = rhs - y

            ! x_new = D^-1 y
            x_new = inv_diag * y
            !call precnd(y, x_new)
            
            ! DIIS extrapolation
            if(do_diis) then
                x_diis(:,nmat) = x_new
                e_diis(:,nmat) = x_new - x
                call diis(n, nmat, diis_max, x_diis, e_diis, bmat, x_new)
            endif

            ! increment
            x = x_new - x
            ! compute norm
            call rmsvec(n, x, rms_norm, max_norm)
            ! update
            x = x_new
            
            write(msg, "('iter=',i4,' residual norm (rms, max): ', 2d14.4)") it, rms_norm, max_norm
            call ommp_message(msg, OMMP_VERBOSE_HIGH)

            ! Check convergence
            if(max_norm < tol) then
                call ommp_message("Required convergence threshold reached, &
                                  &exiting iterative solver.", OMMP_VERBOSE_HIGH)
                exit
            end if
        enddo
        
        call mfree('jacobi_diis_solver [x_new]', x_new)
        call mfree('jacobi_diis_solver [y]', y)
        if(do_diis) then
            call mfree('jacobi_diis_solver [x_diis]', x_diis)
            call mfree('jacobi_diis_solver [e_diis]', e_diis)
            call mfree('jacobi_diis_solver [bmat]', bmat)
        endif
      
        if(max_norm > tol) then
            call fatal_error("Iterative solver did not converged")
        end if

    end subroutine jacobi_diis_solver
  
    subroutine diis(n,nmat,ndiis,x,e,b,xnew)
        !! perform Pulay's direct inversion in the iterative subspace extrapolation:
        use mod_memory, only: mallocate, mfree

        implicit none
        ! TODO doc
        integer(ip), intent(in) :: n, ndiis
        integer(ip), intent(inout) :: nmat
        real(rp), dimension(n, ndiis), intent(inout) :: x, e
        real(rp), dimension(ndiis+1, ndiis+1), intent(inout) :: b
        real(rp), dimension(n), intent(inout) :: xnew
        
        integer(ip) :: nmat1, i, info
        integer(ip) :: j, k

        real(rp),    allocatable :: bloc(:,:), cex(:)
        integer(ip), allocatable :: ipiv(:)

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
    
        call mallocate('diis [bloc]', nmat1, nmat1, bloc)
        call mallocate('diis [cex]', nmat1, cex)
        call mallocate('diis [ipiv]', nmat1, ipiv)
    
        call makeb(n, nmat, ndiis, e, b)
        bloc = b(1:nmat1,1:nmat1)
        cex = 0.0_rp
        cex(1) = 1.0_rp

        call dgesv(nmat1, 1, bloc, nmat1, ipiv, cex, nmat1, info)
        
        if(info /= 0) then
            ! inversion failed. discard the previous points and restart.
            nmat = 1
            call mfree('diis [bloc]', bloc)
            call mfree('diis [cex]', cex)
            call mfree('diis [ipiv]', ipiv)
            return
        end if

        xnew = 0.0_rp
        do i = 1, nmat
            xnew = xnew + cex(i+1)*x(:,i)
        end do
        nmat = nmat + 1
        
        call mfree('diis [bloc]', bloc)
        call mfree('diis [cex]', cex)
        call mfree('diis [ipiv]', ipiv)
  
    end subroutine diis

    subroutine makeb(n,nmat,ndiis,e,b)
        !! assemble the DIIS B matrix:
        implicit none
        
        integer(ip), intent(in) :: n, nmat, ndiis
        real(rp), dimension(n, ndiis), intent(in) :: e
        real(rp), dimension(ndiis+1, ndiis+1), intent(inout) :: b

        integer(ip) :: i
        real(rp) :: bij
      
        if(nmat == 1) then
            ! 1st built:
            !         [ 0 |  1  ]
            !     b = [ --+---- ]
            !         [ 1 | e*e ]
            b(1,1) = 0.0_rp
            b(1,2) = 1.0_rp
            b(2,1) = 1.0_rp
            b(2,2) = dot_product(e(:,1),e(:,1))
        else
            ! subsequent builts
            ! first, update the lagrangian line:
            b(nmat+1,1) = 1.0_rp
            b(1,nmat+1) = 1.0_rp

            ! now, compute the new matrix elements:
            do i = 1, nmat - 1
                bij = dot_product(e(:,i),e(:,nmat))
                b(nmat+1,i+1) = bij
                b(i+1,nmat+1) = bij
            end do
            
            b(nmat+1,nmat+1) = dot_product(e(:,nmat),e(:,nmat))
        end if
    end subroutine makeb

    subroutine rmsvec( n, v, vrms, vmax )
        !! compute root-mean-square and max norms of a vector.
        implicit none
    
        integer(ip), intent(in) :: n
        real(rp), dimension(n), intent(in) :: v
        real(rp), intent(inout) :: vrms, vmax
    
        integer(ip) :: i
    
        ! initialize
        vrms = 0.0_rp
        vmax = 0.0_rp
        
        ! loop over entries
        do i = 1, n
            ! max norm
            vmax = max(vmax,abs(v(i)))
            ! rms norm
            vrms = vrms + v(i)*v(i)
        enddo
        
        vrms = sqrt(vrms/dble(n))
    end subroutine rmsvec

end module mod_solvers
