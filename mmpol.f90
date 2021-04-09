module mmpol
  use precision
!
! mmpol - an opensource library for polarizable molecular mechanics based embedding
!
!   MoLECoLab Pisa - Modelling Light and Environment effects on Complex Systems Lab
!   Department of Chemistry and Industrial Chemistry
!   University of Pisa, Italy
!
! 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                  !
!   this fortran module contains all the shared scalar and array quantities that   !
!   are used by the mmpol library.                                                 !
!   it also contains wrapper routines to allocate memory, keeping track of the     !
!   amount used.                                                                   !
!                                                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   scalar control variables:
!   =========================
!
!     revision control: parameter used to check the consistency between the code
!     and the version of qmip used to generate the input file.
!
  integer(ip), parameter   :: revision = 1
!
!     units used for i/o:
!
  integer(ip), parameter   :: mmpinp = 100, mmpfile = 101, mmpout=6
!
!     maximum allowed number of 1-2, ..., 1-5 neighbors.
!
  integer(ip), parameter   :: maxn12 = 8, maxn13 = 24, maxn14 = 72, maxn15 = 216
!
!     maximum number of members for the same polarization group
!
  integer(ip), parameter   :: maxpgp = 120
!
!     maximum number of neighbors to consider for building the Wang-Skeel preconditioner
!
  integer(ip), parameter   :: ws_max = 100
!
!     maximum number of integer and real parameters (for i/o purposes)
!
  integer(ip), parameter   :: nipars = 32, nrpars = 32
!
!     maximum amount of memory (in gb) and number of cores that can be used by mmpol:
!
  integer(ip)              :: maxcor, nproc
!
!     verbosity flag:
!
  integer(ip)              :: verbose
!
!     force field selection. ff_type is 0 for AMBER and 1 for AMOEBA. for AMBER
!     ff_rules is 0 for Wang AL and 1 for Wang DL exclusion rules, respectively.
!
  integer(ip)              :: ff_type, ff_rules
!
!     parameters that control how the polarization equations are solved.
!     solver:            1: preconditioned conjugate gradient (default)
!                        2: jacobi iterations with DIIS extrapolation
!                        3: matrix inversion
!
  integer(ip)              :: solver         
!
!     matrix_vector:     1: assemble the matrix using O(n^2) storage and use dgemv
!                        2: compute the matrix vector products in a direct fashion (default)
!                        3: use a fast multiplication technique (NYI.)
!
  integer(ip)              :: matrix_vector
!
!     convergence:       convergence threshold (rms norm of the residual/increment) for iterative solvers
!
  real(rp)                 :: convergence
!
!
!
  logical                  :: gradient
!     
!     diagonal shift parameter for the Wang-Skeel preconditioner.
!     default: one for AMBER, two for AMOEBA
!
  real(rp)                 :: ws_shift
!
!     name of the input and scratch files and associated lenght
!
  integer(ip)             :: len_inname, len_scname
  character(len=120)      :: input_file, scratch_file
!
  logical                 :: amoeba
!
!   variables that control the size of various arrays:  
!   ==================================================
!
!     number of MM atoms
!
  integer(ip)              :: mm_atoms
!
!     number of atoms bearing a polarizability:
!
  integer(ip)              :: pol_atoms
!
!     size of the cartesian multipolar distribution (i.e., (l+1)*(l+2)*(l+3)/6)
!     this is 1 for AMBER (charges only), 10 for AMOEBA (up to quadrupoles). 
!     this is also the size of the array that contains the electrostatic properties
!     of the sources at the sources. ld_cder is the leading size of the derivative of
!     such a distribution, which is 3 for AMBER and 19 for AMOEBA.
!
  integer(ip)              :: ld_cart, ld_cder
!
!     number of induced point dipoles distributions
!     this is 1 for AMBER, 2 for AMOEBA. 
!
  integer(ip)              :: n_ipd
!
!   arrays for the force field dependent exclusion factors. 
!   =======================================================
!
!     factors for charge-charge (or multipole-multipole) interactions
!
  real(rp)                 :: mscale(4)
!
!     factors for chrage-ipd (or multipole-ipd) interactions. 
!     in AMOEBA, this is used to define the polarization field, i.e., the right-hand
!     side to the polarization equations, and depends on the connectivity.
!
  real(rp)                 :: pscale(5)
!
!     factors for multipoles-ipd interactions used to compute the direct field,
!     which is used to define the polarization energy. these factors depend on 
!     the polarization group "connectivity" (AMOEBA only)
!
  real(rp)                 :: dscale(4)
!
!     factor for ipd-ipd interactions. these depend on the connectivity (AMBER)
!     or on the polarization group " connectivity (AMOEBA)
!
  real(rp)                 :: uscale(4)
!
!     arrays to pack the integer and real parameters for i/o purposed:
!
  integer(ip)              :: int_pars(nipars)
  real(rp)                 :: real_pars(nrpars)
!
!   arrays:
!   =======
!
!     coordinates of the mm atoms:
!
  real(rp),    allocatable :: cmm(:,:)
!
!     coordinates of the polarizable mm atoms:
!
  real(rp),    allocatable :: cpol(:,:)
!
!     mutlipolar distribution. note that there are two arrays, one for the
!     distribution rotated to the lab frame and one for the force field 
!     parameters. this is only relevant for AMOEBA.
!
  real(rp),    allocatable :: q(:,:), q0(:,:)
!
!     induced point dipoles. note the third dimension to account for more than
!     one ipd distribution in AMOEBA.
!
  real(rp),    allocatable :: ipd(:,:,:)
!
!     polarizabilities:
!
  real(rp),    allocatable :: pol(:)
!
!     indices pointing to the polarizable atoms and viceversa. 
!     these indices identify the subset of the mm atoms that is polarizable (mm_polar)
!     and the position of a polarizable atom in the mm atoms list (polar_mm).
!
  integer(ip), allocatable :: mm_polar(:), polar_mm(:)
!
!     connectivity (number of neighbors and list of neighbors):
!
  integer(ip), allocatable :: n12(:), i12(:,:), n13(:), i13(:,:), n14(:), i14(:,:), n15(:), i15(:,:)
!
!     polarization group or fragment
!
  integer(ip), allocatable :: group(:)
!
!     polarization group "connectivity":
!
  integer(ip), allocatable :: np11(:), ip11(:,:), np12(:), ip12(:,:), np13(:), ip13(:,:), np14(:), ip14(:,:)
!
!   parameters for the definition of the rotation matrices for the multipoles:
!   ==========================================================================
!
!     definition of the molecular frame
!     convention: 0 ... do not rotate
!                 1 ... z-then-x
!                 2 ... bisector
!                 3 ... z-only
!                 4 ... z-bisector
!                 5 ... 3-fold
!
  integer(ip), allocatable :: mol_frame(:)
!
!     neighboring atoms used to define the axes of the molecular frame:
!
  integer(ip), allocatable :: ix(:), iy(:), iz(:)
!
!   scalars and arrays for various useful intermediates and results:
!   ================================================================
!
!     electrostatic and polarization energies, including their breakdown into contributoins:
!
  real(rp)                 :: e_ele, e_pol, e_qd, e_dd
!
!     potential (and higher order terms) at the charges (multipoles)
!     and its derivatives
!
  real(rp),    allocatable :: v_qq(:,:), dv_qq(:,:)
!
!     field of the charges (multipoles) at the ipd and its derivatives
!
  real(rp),    allocatable :: ef_qd(:,:,:), def_qd(:,:,:)
!
!     polarization matrix (only allocated if explicitly requested)
!
  real(rp),    allocatable :: t_matrix(:,:)
!
!     arrays for the wang-skeel preconditioner:
!
  integer(ip), allocatable :: n_ws(:), list_ws(:)
  real(rp),    allocatable :: block_ws(:,:)
!
!     array to store the thole factors for computing damping functions
!
  real(rp),    allocatable :: thole(:)
!
!   constants:
!   ==========
!
  real(rp),    parameter   :: thres = 1.0e-8_rp
  real(rp),    parameter   :: zero = 0.0_rp, pt5 = 0.5_rp, one = 1.0_rp, two = 2.0_rp, three = 3.0_rp, &
                              four = 4.0_rp, five = 5.0_rp, six = 6.0_rp, seven = 7.0_rp, nine = 9.0_rp, &
                              ten = 10.0_rp
!
! -------------------------------------------------------------------------------------------------------
!
!   internal variables for memory allocation and definition of the interface for
!   the wrappers:
!
  integer(ip)              :: maxmem
  integer(ip), parameter   :: memout = 6
  integer(ip), private     :: istat, ltot
!
  interface mallocate
    module procedure r_alloc1
    module procedure r_alloc2
    module procedure r_alloc3
    module procedure i_alloc1
    module procedure i_alloc2
    module procedure i_alloc3
  end interface mallocate

  interface mfree
    module procedure r_free1
    module procedure r_free2
    module procedure r_free3
    module procedure i_free1
    module procedure i_free2
    module procedure i_free3
  end interface mfree
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                  !
!   set of memory allocation routines that keep track of how much memory is used   !
!                                                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  contains
    subroutine r_alloc1(string,len1,v)
      character (len=*),              intent(in)    :: string
      integer(ip),                    intent(in)    :: len1
      real(rp),          allocatable, intent(inout) :: v(:)
!  
      allocate (v(len1), stat=istat)
      call chk_all(string,len1,istat)
      return
    end subroutine r_alloc1
!  
    subroutine r_alloc2(string,len1,len2,v)
      character (len=*),              intent(in)    :: string
      integer(ip),                    intent(in)    :: len1, len2
      real(rp),          allocatable, intent(inout) :: v(:,:)
!  
      allocate (v(len1,len2), stat=istat)
      call chk_all(string,len1*len2,istat)
      return
    end subroutine r_alloc2
!  
    subroutine r_alloc3(string,len1,len2,len3,v)
      character (len=*),              intent(in)    :: string
      integer(ip),                    intent(in)    :: len1, len2, len3
      real(rp),          allocatable, intent(inout) :: v(:,:,:)
!  
      allocate (v(len1,len2,len3), stat=istat)
      call chk_all(string,len1*len2*len3,istat)
      return
    end subroutine r_alloc3
!  
    subroutine i_alloc1(string,len1,v)
      character (len=*),              intent(in)    :: string
      integer(ip),                    intent(in)    :: len1
      integer(ip),       allocatable, intent(inout) :: v(:)
!  
      allocate (v(len1), stat=istat)
      call chk_all(string,len1,istat)
      return
    end subroutine i_alloc1
!  
    subroutine i_alloc2(string,len1,len2,v)
      character (len=*),              intent(in)    :: string
      integer(ip),                    intent(in)    :: len1, len2
      integer(ip),       allocatable, intent(inout) :: v(:,:)
!  
      allocate (v(len1,len2), stat=istat)
      call chk_all(string,len1*len2,istat)
      return
    end subroutine i_alloc2
!  
    subroutine i_alloc3(string,len1,len2,len3,v)
      character (len=*),              intent(in)    :: string
      integer(ip),                    intent(in)    :: len1, len2, len3
      integer(ip),       allocatable, intent(inout) :: v(:,:,:)
!  
      allocate (v(len1,len2,len3), stat=istat)
      call chk_all(string,len1*len2*len3,istat)
      return
    end subroutine i_alloc3
!  
    subroutine chk_all(string,lall,istat)
      implicit none
      integer(ip),       intent(in) :: lall, istat
      character (len=*), intent(in) :: string
!  
!     memory error format:
!  
   9000 format(t3,'allocation error in subroutine ',a,'. stat= ',i5)
   9010 format(t3,'allocation error in subroutine ',a,'.',/,    &
               t3,'not enough memory. ',i8,' words required',/, &
               t3,'                   ',i8,' words available.')
!  
      if (istat.ne.0) then
        write(memout,9000) string, istat
        stop
      else if (lall.gt.maxmem) then
        write(memout,9010) string, lall, maxmem
        stop
      else
        maxmem = maxmem - lall
      end if
      return
    end subroutine chk_all
!  
    subroutine r_free1(string,v)
      character (len=*),              intent(in)    :: string
      real(rp),          allocatable, intent(inout) :: v(:)
!  
      if (.not. allocated(v)) return
      ltot = size(v)
      deallocate (v, stat=istat)
      call chk_free(string,ltot,istat)
      return
    end subroutine r_free1
!  
    subroutine r_free2(string,v)
      character (len=*),              intent(in)    :: string
      real(rp),          allocatable, intent(inout) :: v(:,:)
!  
      if (.not. allocated(v)) return
      ltot = size(v)
      deallocate (v, stat=istat)
      call chk_free(string,ltot,istat)
      return
    end subroutine r_free2
!  
    subroutine r_free3(string,v)
      character (len=*),              intent(in)    :: string
      real(rp),          allocatable, intent(inout) :: v(:,:,:)
!  
      if (.not. allocated(v)) return
      ltot = size(v)
      deallocate (v, stat=istat)
      call chk_free(string,ltot,istat)
      return
    end subroutine r_free3
!  
    subroutine i_free1(string,v)
      character (len=*),              intent(in)    :: string
      integer(ip),       allocatable, intent(inout) :: v(:)
!  
      if (.not. allocated(v)) return
      ltot = size(v)
      deallocate (v, stat=istat)
      call chk_free(string,ltot,istat)
      return
    end subroutine i_free1
!  
    subroutine i_free2(string,v)
      character (len=*),              intent(in)    :: string
      integer(ip),       allocatable, intent(inout) :: v(:,:)
!  
      if (.not. allocated(v)) return
      ltot = size(v)
      deallocate (v, stat=istat)
      call chk_free(string,ltot,istat)
      return
    end subroutine i_free2
!  
    subroutine i_free3(string,v)
      character (len=*),              intent(in)    :: string
      integer(ip),       allocatable, intent(inout) :: v(:,:,:)
!  
      if (.not. allocated(v)) return
      ltot = size(v)
      deallocate (v, stat=istat)
      call chk_free(string,ltot,istat)
      return
    end subroutine i_free3
!  
    subroutine chk_free(string,lfree,istat)
      implicit none
      integer,           intent(in) :: lfree, istat
      character (len=*), intent(in) :: string
!  
!     memory error format:
!  
   9000 format(t3,'deallocation error in subroutine ',a,'. stat= ',i5)
!  
      if (istat.ne.0) then
        write(memout,9000) string, istat
        stop
      else
        maxmem = maxmem + lfree
      end if
      return
    end subroutine chk_free
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                  !
!   silly routine to abort the calculation with an error message:                  !
!                                                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    subroutine fatal_error(message)
      implicit none
      character (len=*) message
   1000 format(t3,a)
      write(memout,1000) message
      stop '   error termination for open_mmpol.'
      return
    end subroutine fatal_error
!
    subroutine set_screening_parameters
      implicit none
      real(rp), parameter :: pt4 = 0.40_rp, pt8 = 0.80_rp
      if (ff_type.eq.0 .and. ff_rules.eq.0) then
!
!       exclusion rules for WangAL
!
        mscale(1) = zero
        mscale(2) = zero
        mscale(3) = one
        mscale(4) = one
        pscale(1) = zero
        pscale(2) = zero
        pscale(3) = one
        pscale(4) = one
        pscale(5) = one
        dscale(1) = zero
        dscale(2) = zero
        dscale(3) = one
        dscale(4) = one
        uscale(1) = zero
        uscale(2) = zero
        uscale(3) = one
        uscale(4) = one
      else if (ff_type.eq.0 .and. ff_rules.eq.1) then
!
!       exclusion rules for WangDL
!
        mscale(1) = one
        mscale(2) = one
        mscale(3) = one
        mscale(4) = one
        pscale(1) = one
        pscale(2) = one
        pscale(3) = one
        pscale(4) = one
        pscale(5) = one
        dscale(1) = one
        dscale(2) = one
        dscale(3) = one
        dscale(4) = one
        uscale(1) = one
        uscale(2) = one
        uscale(3) = one
        uscale(4) = one
      else if (ff_type.eq.1) then
!
!       exclusion rules for AMOEBA
!
        mscale(1) = zero
        mscale(2) = zero
        mscale(3) = pt4
        mscale(4) = pt8
        pscale(1) = zero
        pscale(2) = zero
        pscale(3) = one
        pscale(4) = one
        pscale(5) = pt5
        dscale(1) = zero
        dscale(2) = one
        dscale(3) = one
        dscale(4) = one
        uscale(1) = one
        uscale(2) = one
        uscale(3) = one
        uscale(4) = one
      else
        call fatal_error('the required force field is not implemented.')
      end if
    end subroutine set_screening_parameters

end module mmpol
