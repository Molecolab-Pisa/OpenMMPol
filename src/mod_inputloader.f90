module mod_inputloader
    implicit none
    private

    public :: mmpol_init_from_mmp

    contains
    
    subroutine skip_lines(f, n)
        !! This subroutine just skip lines while reading an input
        !! file
        use mod_memory, only: ip
        
        implicit none

        !! unit file
        integer(ip), intent(in) :: f
        !! number of line to be skipped
        integer(ip), intent(in) :: n

        integer(ip) :: i
        ! character(len=512) :: s

        do i=1, n
            read(f, *)
            ! read(f, *) s
            ! write(6, *) "SKIPPING ", i, "/", n, "  '", trim(s),"'"
        end do

    end subroutine skip_lines

    subroutine mmpol_init_from_mmp(input_file)
        !! This function read a .mmp file (revision 2 and 3) are supported
        !! and initialize all the quantities need to describe the environment
        !! within this library.

        use mod_mmpol, only: cmm, q, q0, pol
        use mod_mmpol, only: mm_atoms, pol_atoms, ff_rules, &
                             ff_type, amoeba, &
                             polar_mm,  conn, mmat_polgrp
        use mod_mmpol, only: mol_frame, iz, ix, iy
        use mod_mmpol, only: fatal_error, mmpol_prepare, mmpol_init
        
        use mod_memory, only: ip, rp, mfree, mallocate, memory_init
        use mod_io, only: mmpol_print_summary, iof_mmpinp
        use mod_constants, only: zero, ten, thres, angstrom2au
        use mod_adjacency_mat, only: adj_mat_from_conn

        implicit none

        character(len=*), intent(in) :: input_file
        !! name of the input file

        ! read the input for the mmpol calculation and process it.
        integer(ip) :: input_revision, iconv
        integer(ip) :: maxcor, nproc
        integer(ip) :: my_mm_atoms, my_pol_atoms, my_ff_type, my_ff_rules
        integer(ip) :: my_ld_cart, verbosity
        
        integer(ip) :: i, j
        
        real(rp), allocatable :: my_pol(:), my_cmm(:,:), my_q(:,:)
        integer(ip), allocatable :: pol_atoms_list(:), my_ip11(:,:)
        
        integer(ip), parameter :: maxn12 = 8
        integer(ip), parameter :: maxpgp = 120
        !! maximum number of members for the same polarization group
        integer(ip), allocatable :: i12(:,:)
        
        ! open the (formatted) input file
        open (unit=iof_mmpinp, &
              file=input_file(1:len(trim(input_file))), &
              form='formatted', &
              access='sequential')

        ! Read input revision, supported revisions are 2 and 3.
        read(iof_mmpinp,*) input_revision
        
        if (input_revision /= 3 .and. input_revision /= 2) &
            call fatal_error('input and internal revision conflict.')

        call memory_init(.false., 0)
        
        call skip_lines(iof_mmpinp, 2)
        read(iof_mmpinp,*) my_ff_type
        read(iof_mmpinp,*) my_ff_rules
        
        if(input_revision == 3) then
            call skip_lines(iof_mmpinp, 17)
        else if(input_revision == 2) then
            call skip_lines(iof_mmpinp, 15)
        end if
        
        read(iof_mmpinp,*) my_mm_atoms
        
        if(my_ff_type == 1) then
            my_ld_cart = 10
        else
            my_ld_cart = 1
        end if
        
        call mallocate('mmpol_init_from_mmp [my_cmm]', 3_ip, my_mm_atoms, my_cmm)
        call mallocate('mmpol_init_from_mmp [my_q]', my_ld_cart, my_mm_atoms, my_q)
        call mallocate('mmpol_init_from_mmp [my_pol]', my_mm_atoms, my_pol)
        call mallocate('mmpol_init_from_mmp [my_pol_atoms]', my_mm_atoms, pol_atoms_list)
        call mallocate('mmpol_init_from_mmp [my_ip11]', maxpgp, my_mm_atoms, my_ip11)
       
        call skip_lines(iof_mmpinp, my_mm_atoms+1) ! Skip a zero and the section of atomic numbers

        ! coordinates:
        do i = 1, my_mm_atoms
            read(iof_mmpinp,*) my_cmm(1:3,i)
        end do

        call skip_lines(iof_mmpinp, my_mm_atoms) ! Skip section of residues number

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
        ! TODO this is more efficiently done with pack and count
        my_pol_atoms = 0
        pol_atoms_list(:) = 0
        do i = 1, my_mm_atoms
            if (my_pol(i) > thres) then
                my_pol_atoms = my_pol_atoms + 1
                pol_atoms_list(my_pol_atoms) = i
            end if
        end do
        
        ! remove null polarizabilities from the list
        do i = 1, my_pol_atoms
            my_pol(i) = my_pol(pol_atoms_list(i))
        end do
        
        ! mmpol module initialization
        call mmpol_init(my_ff_type, my_ff_rules, my_pol_atoms, my_mm_atoms)
        
        ! Copy data in the correct units (this means AU)
        cmm = my_cmm * angstrom2au
        q = my_q
        if(amoeba) then
            my_pol = my_pol*angstrom2au**3
        end if
        pol = my_pol
        polar_mm = pol_atoms_list
        
        call mfree('mmpol_init_from_mmp [my_cmm]', my_cmm)
        call mfree('mmpol_init_from_mmp [my_q]', my_q)
        call mfree('mmpol_init_from_mmp [pol]', my_pol)

        ! 1-2 connectivity:
        call mallocate('mmpol_init_from_mmp [i12]', maxn12, mm_atoms, i12)
        do i = 1, mm_atoms
            read(iof_mmpinp,*) i12(1:maxn12,i)
        end do
        
        ! Writes the adjacency matrix in Yale sparse format in conn(1)
        call adj_mat_from_conn(i12, conn(1)) 
        call mfree('mmpol_init_from_mmp [i12]', i12)
         
        if(amoeba) then
            ! group 11 connectivity:
            ! (to be replaced with polarization group)
            do i = 1, mm_atoms
                read(iof_mmpinp,*) my_ip11(1:maxpgp,i)
            end do

            call polgroup11_to_mm2pg(my_ip11, mmat_polgrp)
            call mfree('mmpol_init_from_mmp [my_ip11]', my_ip11)

            ! information to rotate the multipoles to the lab frame.
            ! mol_frame, iz, ix, iy:
            do i = 1, mm_atoms
                read(iof_mmpinp,*) mol_frame(i), iz(i), ix(i), iy(i)
            end do
        end if
        ! now, process the input, create all the required arrays 
        ! and the correspondence lists:
     
        call mmpol_prepare()

    end subroutine mmpol_init_from_mmp

    subroutine polgroup11_to_mm2pg(polgroup_neigh, mm2pol)
        !! Take as input a matrix in which the n-th row stores the atoms
        !! that are in the same polarization group as the n-th atom (rows are 
        !! padded with zeros) and assign to each atom a group index, according to
        !! the input information. This is a way of compressing and making more 
        !! clear the handling of polarization groups.

        use mod_memory, only : ip
        
        implicit none

        integer(ip), intent(in) :: polgroup_neigh(:,:)
        integer(ip), intent(out) :: mm2pol(:)

        integer(ip) :: i, j, maxn, ipg, mm_atoms

        mm2pol = 0
        ipg = 1
        maxn = size(polgroup_neigh, 1)
        mm_atoms = size(polgroup_neigh, 2)

        do i=1, mm_atoms
            if(mm2pol(i) == 0) then
                ! This is a new group
                do j=1, maxn
                    if(polgroup_neigh(j,i) == 0) exit
                    if(mm2pol(polgroup_neigh(j,i)) /= 0) write(*, *) "MB22 Unexpected error!" !TODO
                    mm2pol(polgroup_neigh(j,i)) = ipg
                end do

                ipg = ipg + 1
            end if
        end do

        if( any(mm2pol == 0) ) write(*, *) "MB22 2 unexpected error" !TODO

    end subroutine polgroup11_to_mm2pg

end module mod_inputloader
