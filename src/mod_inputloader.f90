module mod_inputloader
    use mod_io, only: ommp_message

    implicit none
    private

    public :: mmpol_init_from_mmp, mmpol_init_from_xyz

    contains
    
    subroutine mmpol_init_from_mmp(input_file)
        !! This function read a .mmp file (revision 2 and 3) are supported
        !! and initialize all the quantities need to describe the environment
        !! within this library.

        use mod_mmpol, only: cmm, q, pol
        use mod_mmpol, only: mm_atoms, amoeba, &
                             polar_mm, conn, mmat_polgrp, thole_scale
        use mod_mmpol, only: mol_frame, iz, ix, iy
        use mod_mmpol, only: fatal_error, mmpol_prepare, mmpol_init, &
                             set_screening_parameters
        
        use mod_memory, only: ip, rp, mfree, mallocate, memory_init
        use mod_io, only: iof_mmpinp
        use mod_adjacency_mat, only: adj_mat_from_conn
        use mod_utils, only: skip_lines
        use mod_constants, only: angstrom2au, OMMP_VERBOSE_DEBUG
        use mod_constants, only: mscale_wang_al, &
                                 pscale_wang_al, &
                                 dscale_wang_al, &
                                 uscale_wang_al, & 
                                 mscale_wang_dl, &
                                 pscale_wang_dl, &
                                 dscale_wang_dl, &
                                 uscale_wang_dl, & 
                                 mscale_amoeba, &
                                 pscale_amoeba, &
                                 dscale_amoeba, &
                                 uscale_amoeba, &
                                 pscale_intra_amoeba, &
                                 thole_scale_wang_dl, &
                                 thole_scale_wang_al

        implicit none

        character(len=*), intent(in) :: input_file
        !! name of the input file

        ! read the input for the mmpol calculation and process it.
        integer(ip) :: input_revision
        integer(ip) :: my_mm_atoms, my_pol_atoms, my_ff_type, my_ff_rules
        integer(ip) :: my_ld_cart
        
        integer(ip) :: i, ist
        
        real(rp), allocatable :: my_pol(:), my_cmm(:,:), my_q(:,:)
        integer(ip), allocatable :: pol_atoms_list(:), my_ip11(:,:)
        
        integer(ip), parameter :: maxn12 = 8
        !! maximum number of neighbour atoms
        integer(ip), parameter :: maxpgp = 120
        !! maximum number of members for the same polarization group
        real(rp), parameter :: thres = 1e-8
        !! Threshold used to decide if a polarizability is zero
        ! TODO why we do not use eps_rp directly?

        integer(ip), allocatable :: i12(:,:)
        character(len=120) :: msg
       
        write(msg, "(A)") "Reading MMP file: "//input_file(1:len(trim(input_file)))
        call ommp_message(msg, OMMP_VERBOSE_DEBUG)

        ! open the (formatted) input file
        open (unit=iof_mmpinp, &
              file=input_file(1:len(trim(input_file))), &
              form='formatted', &
              access='sequential', &
              iostat=ist)

        if(ist /= 0) then
            call fatal_error('Error while opening MMP input file')
        end if

        call ommp_message("Reading input parameters", OMMP_VERBOSE_DEBUG)

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
        
        call ommp_message("Allocating memory", OMMP_VERBOSE_DEBUG)

        call mallocate('mmpol_init_from_mmp [my_cmm]', 3_ip, my_mm_atoms, my_cmm)
        call mallocate('mmpol_init_from_mmp [my_q]', my_ld_cart, my_mm_atoms, my_q)
        call mallocate('mmpol_init_from_mmp [my_pol]', my_mm_atoms, my_pol)
        call mallocate('mmpol_init_from_mmp [my_pol_atoms]', my_mm_atoms, pol_atoms_list)
        call mallocate('mmpol_init_from_mmp [my_ip11]', maxpgp, my_mm_atoms, my_ip11)
       
        call skip_lines(iof_mmpinp, my_mm_atoms+1) ! Skip a zero and the section of atomic numbers
        
        call ommp_message("Reading coordinates", OMMP_VERBOSE_DEBUG)
        ! coordinates:
        do i = 1, my_mm_atoms
            read(iof_mmpinp,*) my_cmm(1:3,i)
        end do

        call skip_lines(iof_mmpinp, my_mm_atoms) ! Skip section of residues number

        call ommp_message("Reading fixed multipoles", OMMP_VERBOSE_DEBUG)
        ! charges/multipoles:
        do i = 1, my_mm_atoms
            read(iof_mmpinp,*) my_q(1:my_ld_cart,i)
        end do

        call ommp_message("Reading polarizabilities", OMMP_VERBOSE_DEBUG)
        ! polarizabilities:
        my_pol = 0.0_rp
        do i = 1, my_mm_atoms
            read(iof_mmpinp,*) my_pol(i)
        end do

        call ommp_message("Processing polarizabilities", OMMP_VERBOSE_DEBUG)
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
        my_pol(my_pol_atoms+1:my_mm_atoms) = 0.0_rp
        
        call ommp_message("Initializing open mmpol module", OMMP_VERBOSE_DEBUG)
        ! mmpol module initialization
        call mmpol_init(my_ff_type, my_mm_atoms, my_pol_atoms)
        if(amoeba) then
            call set_screening_parameters(mscale_amoeba, &
                                          pscale_amoeba, &
                                          dscale_amoeba, &
                                          uscale_amoeba, &
                                          pscale_intra_amoeba)
        else
            if(my_ff_rules == 1) then
                call set_screening_parameters(mscale_wang_dl, &
                                              pscale_wang_dl, &
                                              dscale_wang_dl, &
                                              uscale_wang_dl)
                thole_scale = thole_scale_wang_dl
            else if(my_ff_rules == 0) then
                call set_screening_parameters(mscale_wang_al, &
                                              pscale_wang_al, &
                                              dscale_wang_al, &
                                              uscale_wang_al)
                thole_scale = thole_scale_wang_al
            end if
        end if 
        
        call ommp_message("Converting input units to A.U.", OMMP_VERBOSE_DEBUG)
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

        call ommp_message("Processing connectivity informations", OMMP_VERBOSE_DEBUG)

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
     
        call ommp_message("Populating utility arrays", OMMP_VERBOSE_DEBUG)
        call mmpol_prepare()
        
        call ommp_message("Initialization from MMP file done.", OMMP_VERBOSE_DEBUG)

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

    subroutine mmpol_init_from_xyz(xyz_file, prm_file)
        !! This function read a .mmp file (revision 2 and 3) are supported
        !! and initialize all the quantities need to describe the environment
        !! within this library.
        use mod_mmpol, only:cmm, polar_mm, conn
        use mod_mmpol, only: fatal_error, mmpol_prepare, mmpol_init
        
        use mod_memory, only: ip, mfree, mallocate, memory_init
        use mod_constants, only: angstrom2au, OMMP_VERBOSE_DEBUG, OMMP_FF_AMOEBA
        use mod_adjacency_mat, only: adj_mat_from_conn, yale_sparse, &
                                     build_conn_upto_n
        use mod_prm, only: assign_vdw, assign_pol, assign_mpoles, assign_bond, &
                           assign_angle, assign_urey, assign_strbnd, &
                           assign_opb, assign_pitors, assign_torsion, &
                           assign_tortors, assign_angtor, assign_strtor, &
                           check_keyword, terminate_prm
        use mod_utils, only: starts_with_alpha, isreal, isint, tokenize

        implicit none

        character(len=*), intent(in) :: xyz_file
        !! name of the input XYZ file
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file

        integer(ip), parameter :: iof_xyzinp = 200, &
                                  maxn12 = 8
        integer(ip) :: my_mm_atoms, ist, i, j, atom_id, tokb, toke
        integer(ip), allocatable :: i12(:,:), attype(:)
        character(len=120) :: line, msg
        type(yale_sparse) :: adj

        
        write(msg, "(A)") "Reading XYZ file: "//xyz_file(1:len(trim(xyz_file)))
        call ommp_message(msg, OMMP_VERBOSE_DEBUG)

        ! open tinker xyz file
        open(unit=iof_xyzinp, &
             file=xyz_file(1:len(trim(xyz_file))), &
             form='formatted', &
             access='sequential', &
             iostat=ist)

        if(ist /= 0) then
            call fatal_error('Error while opening XYZ input file')
        end if

        ! First line contains as first word the number of atoms and
        ! then a comment that could be ignored.
        read(iof_xyzinp, *) my_mm_atoms
        
        ! Initialize the mmpol module
        ! TODO I'm assuming that it is AMOEBA and fully polarizable
        call mmpol_init(OMMP_FF_AMOEBA, my_mm_atoms, my_mm_atoms)        
        do i=1, my_mm_atoms
           polar_mm(i) = i 
        end do

        ! Temporary quantities that are only used during the initialization
        call mallocate('mmpol_init_from_xyz [attype]', my_mm_atoms, attype)
        call mallocate('mmpol_init_from_xyz [i12]', maxn12, my_mm_atoms, i12)
        attype = 0_ip
        i12 = 0_ip

        do i=1, my_mm_atoms
            read(iof_xyzinp, '(A)') line

            ! First token contains an atom ID. Only sequential numbering is
            ! currently supported.
            tokb = tokenize(line)
            toke = tokenize(line, tokb)

            read(line(tokb:toke), *) atom_id
            if(atom_id /= i) then
                call fatal_error('Non-sequential atom ids in xyz cannot be handled')
            end if
            
            ! This token should contain an atom name, so it should start
            ! with a letter. If this is not true, an unexpected error should
            ! be raised.
            tokb=toke+1
            toke = tokenize(line, tokb)
            if(.not. starts_with_alpha(line(tokb:toke))) then
                call fatal_error('Atom symbol missing or PBC string present in XYZ')
            end if

            tokb=toke+1
            toke = tokenize(line, tokb, 3)
            ! The remaining part contains cartesian coordinates, atom type
            ! and connectivity for the current atom.
            read(line(tokb:toke), *) cmm(1,i), cmm(2,i), cmm(3,i)
            
            tokb=toke+1
            toke = tokenize(line, tokb)
            read(line(tokb:toke), *) attype(i)

            do j=1, maxn12
                tokb=toke+1
                toke = tokenize(line, tokb)
                if(toke < 0) exit
                read(line(tokb:toke), *) i12(j,i)
            end do
            
        end do
        ! Close XYZ file after reading
        close(iof_xyzinp)

        cmm = cmm * angstrom2au
        
        ! Writes the adjacency matrix in Yale sparse format in adj and then
        ! build the connectivity up to 4th order. This is needed here to be
        ! able to assign the parameters
        call adj_mat_from_conn(i12, adj)
        call build_conn_upto_n(adj, 4, conn, .false.)

        call mfree('mmpol_init_from_xyz [i12]', i12)
        
        if( .not. check_keyword(prm_file)) then
            call fatal_error("PRM file cannot be completely understood")
        end if
        
        call ommp_message("Assigning electrostatic parameters", OMMP_VERBOSE_DEBUG)
        call assign_pol(prm_file, attype)
        call assign_mpoles(prm_file, attype)
        
        call ommp_message("Assigning non-bonded parameters", OMMP_VERBOSE_DEBUG)
        call assign_vdw(prm_file, attype)
        
        call ommp_message("Assigning bonded parameters", OMMP_VERBOSE_DEBUG)
        call assign_bond(prm_file, attype)
        call assign_angle(prm_file, attype)
        call assign_urey(prm_file, attype)
        call assign_strbnd(prm_file, attype)
        call assign_opb(prm_file, attype)
        call assign_pitors(prm_file, attype)
        call assign_torsion(prm_file, attype)
        call assign_tortors(prm_file, attype)
        call assign_angtor(prm_file, attype)
        call assign_strtor(prm_file, attype)

        call terminate_prm()
        call mfree('mmpol_init_from_xyz [attype]', attype)
        
        call mmpol_prepare()

    end subroutine mmpol_init_from_xyz

end module mod_inputloader
