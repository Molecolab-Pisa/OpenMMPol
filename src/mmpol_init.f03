subroutine mmpol_init_from_mmp(input_file)
    use mod_mmpol
    use mod_memory, only: ip, rp, mfree, mallocate, memory_init
    use mod_io, only: iof_mmpinp
    use mod_constants, only: zero, ten, thres

    implicit none
    character(len=*), intent(in) :: input_file
    !! name of the input file

    integer(ip), parameter :: revision = 1
    !! mmpol file revision expected from the code
    ! TODO
        
    ! read the input for the mmpol calculation and process it.
    integer(ip) :: input_revision, iconv
    integer(ip) :: maxcor, nproc
    integer(ip) :: my_mm_atoms, my_pol_atoms, my_ff_type, my_ff_rules
    integer(ip) :: my_ld_cart, verbosity
    
    integer(ip) :: i
    
    real(rp), allocatable :: my_pol(:), my_cmm(:,:), my_q(:,:)
    integer(ip), allocatable :: my_group(:)
    
    ! open the (formatted) input file
    open (unit=iof_mmpinp, &
          file=input_file(1:len(trim(input_file))), &
          form='formatted', &
          access='sequential')

    ! start reading the integer control parameters:

    read(iof_mmpinp,*) input_revision
    if (input_revision .ne. revision) &
        call fatal_error('input and internal revision conflict.')
    read(iof_mmpinp,*) maxcor, nproc
    
    call memory_init(.true., maxcor*1024*1024*1024/8)
    
    read(iof_mmpinp,*) verbosity
    read(iof_mmpinp,*) my_ff_type
    read(iof_mmpinp,*) my_ff_rules
    read(iof_mmpinp,*) solver
    read(iof_mmpinp,*) matrix_vector
    read(iof_mmpinp,*) iconv
    read(iof_mmpinp,*) my_mm_atoms
    
    ! decode a few scalar parameters:
    convergence = ten**(-iconv)

    ! TODO improve this part
    if(my_ff_type == 1) then
        my_ld_cart = 10
    else
        my_ld_cart = 1
    end if
    
    call mallocate('mmpol_init_from_mmp [my_cmm]', 3_ip, my_mm_atoms, my_cmm)
    call mallocate('mmpol_init_from_mmp [my_group]', my_mm_atoms, my_group)
    call mallocate('mmpol_init_from_mmp [my_q]', my_ld_cart, my_mm_atoms, my_q)
    call mallocate('mmpol_init_from_mmp [my_pol]', my_mm_atoms, my_pol)
   
    ! coordinates:
    do i = 1, my_mm_atoms
        read(iof_mmpinp,*) my_cmm(1:3,i)
    end do

    ! group/fragment/residue:
    do i = 1, my_mm_atoms
        read(iof_mmpinp,*) my_group(i)
    end do
    
    ! charges/multipoles:
    do i = 1, my_mm_atoms
        read(iof_mmpinp,*) my_q(1:my_ld_cart,i)
    end do

    ! polarizabilities:
    my_pol = 0.0_rp
    do i = 1, my_mm_atoms
        read(iof_mmpinp,*) my_pol(i)
    end do

    ! count how many atoms are polarizable:
    my_pol_atoms = 0
    do i = 1, my_mm_atoms
        if (my_pol(i).gt.thres) my_pol_atoms = my_pol_atoms + 1
    end do
    
    call mmpol_init(my_ff_type, my_ff_rules, my_pol_atoms, my_mm_atoms)
    call set_verbosity(verbosity)
    cmm = my_cmm
    group = my_group
    q = my_q
    pol = my_pol(:pol_atoms)
    
    call mfree('mmpol_init_from_mmp [my_cmm]', my_cmm)
    call mfree('mmpol_init_from_mmp [my_group]', my_group)
    call mfree('mmpol_init_from_mmp [my_q]', my_q)
    
    call print_header()

    ! 1-2 connectivity:
    do i = 1, mm_atoms
        read(iof_mmpinp,*) i12(1:maxn12,i)
    end do
    
    if(amoeba) then
        ! group 11 connectivity:
        ! (to be replaced with polarization group)
        do i = 1, mm_atoms
            read(iof_mmpinp,*) ip11(1:maxpgp,i)
        end do

        ! information to rotate the multipoles to the lab frame.
        ! mol_frame, iz, ix, iy:
        do i = 1, mm_atoms
            read(iof_mmpinp,*) mol_frame(i), iz(i), ix(i), iy(i)
        end do
    end if
    ! now, process the input, create all the required arrays 
    ! and the correspondence lists:
  
    call mmpol_process(my_pol)
    call mfree('mmpol_init_from_mmp [pol]', my_pol)

end subroutine mmpol_init_from_mmp
