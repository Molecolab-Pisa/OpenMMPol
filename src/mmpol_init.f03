subroutine mmpol_init
  use mmpol
  use mod_memory, only: ip, rp
  use mod_io, only: iof_mmpinp
  implicit none
!
! read the input for the mmpol calculation and process it.
!
  integer(ip)              :: input_revision, iconv, i, j
!
! integer(ip), allocatable :: 
  real(rp),    allocatable :: polar(:)
!
!
! open the (formatted) input file and the (binary) internal scratch file:
!
  len_inname = len(trim(input_file))
  open (unit=iof_mmpinp, file=input_file(1:len_inname),form='formatted',access='sequential')
!
! start reading the integer control parameters:
!
  read(iof_mmpinp,*) input_revision
  if (input_revision .ne. revision) call fatal_error('input and internal revision conflict.')
  read(iof_mmpinp,*) maxcor, nproc
  call memory_init(.true., maxcor*1024*1024*1024/8)
  read(iof_mmpinp,*) verbose
  read(iof_mmpinp,*) ff_type
  read(iof_mmpinp,*) ff_rules
  read(iof_mmpinp,*) solver
  read(iof_mmpinp,*) matrix_vector
  read(iof_mmpinp,*) iconv
  read(iof_mmpinp,*) mm_atoms
!
! decode a few scalar parameters:
!
  convergence = ten**(-iconv)
!
  amoeba    = ff_type .eq. 1
  ld_cart   = 1
  ld_cder   = 3
  n_ipd     = 1
  if (amoeba) then
    ld_cart = 10
    ld_cder = 19
    n_ipd   = 2
  end if
!
! allocate memory for the mmpol parameters:
!
  call mallocate('mmpol_init [cmm]',int(3,ip),mm_atoms,cmm)
  call mallocate('mmpol_init [q]',ld_cart,mm_atoms,q)
  call mallocate('mmpol_init [polar]',mm_atoms,polar)
  call mallocate('mmpol_init [n12]',mm_atoms,n12)
  call mallocate('mmpol_init [i12]',maxn12,mm_atoms,i12)
  call mallocate('mmpol_init [n12]',mm_atoms,n13)
  call mallocate('mmpol_init [i13]',maxn13,mm_atoms,i13)
  call mallocate('mmpol_init [n14]',mm_atoms,n14)
  call mallocate('mmpol_init [i14]',maxn14,mm_atoms,i14)
  call mallocate('mmpol_init [group]',mm_atoms,group)
  if (amoeba) then
    call mallocate('mmpol_init [q0]',ld_cart,mm_atoms,q0)
    call mallocate('mmpol_init [n15',mm_atoms,n15)
    call mallocate('mmpol_init [i15',maxn15,mm_atoms,i15)
    call mallocate('mmpol_init [np11',mm_atoms,np11)
    call mallocate('mmpol_init [ip11',maxpgp,mm_atoms,ip11)
    call mallocate('mmpol_init [np12',mm_atoms,np12)
    call mallocate('mmpol_init [ip12',maxpgp,mm_atoms,ip12)
    call mallocate('mmpol_init [np13',mm_atoms,np13)
    call mallocate('mmpol_init [ip13',maxpgp,mm_atoms,ip13)
    call mallocate('mmpol_init [np14',mm_atoms,np14)
    call mallocate('mmpol_init [ip14',maxpgp,mm_atoms,ip14)
    call mallocate('mmpol_init [mol_frame]',mm_atoms,mol_frame)
    call mallocate('mmpol_init [ix]',mm_atoms,ix)
    call mallocate('mmpol_init [iy]',mm_atoms,iy)
    call mallocate('mmpol_init [iz]',mm_atoms,iz)
  end if
!
! read the input file:
!
! coordinates:
!
  do i = 1, mm_atoms
    read(iof_mmpinp,*) cmm(1:3,i)
  end do
!
! group/fragment/residue:
!
  do i = 1, mm_atoms
    read(iof_mmpinp,*) group(i)
  end do
!
! charges/multipoles:
!
  do i = 1, mm_atoms
    read(iof_mmpinp,*) q(1:ld_cart,i)
  end do
!
! polarizabilities:
!
  do i = 1, mm_atoms
    read(iof_mmpinp,*) polar(i)
  end do
!
! count how many atoms are polarizable:
!
  pol_atoms = 0
  do i = 1, mm_atoms
    if (polar(i).gt.thres) pol_atoms = pol_atoms + 1
  end do
  call print_header
!
! allocate memory for the polarizabilities array and lists
!
  call mallocate('mmpol_init [pol]',pol_atoms,pol)
  call mallocate('mmpol_init [cpol]',int(3,ip),pol_atoms,cpol)
  call mallocate('mmpol_init [polar_mm]',pol_atoms,polar_mm)
  call mallocate('mmpol_init [mm_polar]',mm_atoms,mm_polar)
!
! 1-2 connectivity:
!
  do i = 1, mm_atoms
    read(iof_mmpinp,*) i12(1:maxn12,i)
  end do
!
! the following input is only relevant for amoeba:
!
  if (amoeba) then
!
!   group 11 connectivity:
!   (to be replaced with polarization group)
!
    do i = 1, mm_atoms
      read(iof_mmpinp,*) ip11(1:maxpgp,i)
    end do
!
!   information to rotate the multipoles to the lab frame.
!   mol_frame, iz, ix, iy:
!
    do i = 1, mm_atoms
      read(iof_mmpinp,*) mol_frame(i), iz(i), ix(i), iy(i)
    end do
  end if
!
! now, process the input, create all the required arrays and the correspondence lists:
!
  
  call mmpol_process(polar)
  call mfree('mmpol_init [pol]',polar)
!
!
! allocate memory for the electrostatic properties of static multipoles
!
  call mallocate('mmpol_init [v_qq]',ld_cart,mm_atoms,v_qq)
  call mallocate('mmpol_init [ef_qd]',int(3,ip),pol_atoms,n_ipd,ef_qd)
  call mallocate('mmpol_init [dv_qq]',ld_cder,mm_atoms,dv_qq)
  call mallocate('mmpol_init [def_qd]',int(6,ip),pol_atoms,n_ipd,def_qd)
  call mallocate('mmpol_init [v_dq]',ld_cart,mm_atoms,v_dq)
  call mallocate('mmpol_init [ef_dd]',int(3,ip),pol_atoms,n_ipd,ef_dd)
  call mallocate('mmpol_init [dv_dq]',ld_cder,mm_atoms,dv_dq)
  call mallocate('mmpol_init [def_dd]',int(6,ip),pol_atoms,n_ipd,def_dd)
  v_qq   = Zero
  ef_qd  = Zero
  dv_qq  = Zero
  def_qd = Zero
  v_dq   = Zero
  ef_dd  = Zero
  dv_dq  = Zero
  def_dd = Zero
!
! allocate memory for the induced point dipoles
!
  call mallocate('mmpol_init [idp]',int(3,ip),pol_atoms,n_ipd,ipd) 
  ipd    = Zero

end subroutine mmpol_init
