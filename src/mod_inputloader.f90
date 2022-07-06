module mod_inputloader
    implicit none
    private

    public :: mmpol_init_from_mmp, mmpol_init_from_xyz

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

        use mod_mmpol, only: verbose, cmm, q, pol
        use mod_mmpol, only: mm_atoms, amoeba, &
                             polar_mm, conn, mmat_polgrp
        use mod_mmpol, only: mol_frame, iz, ix, iy
        use mod_mmpol, only: fatal_error, mmpol_prepare, mmpol_init
        
        use mod_memory, only: ip, rp, mfree, mallocate, memory_init
        use mod_io, only: mmpol_print_summary, iof_mmpinp
        use mod_constants, only: angstrom2au, OMMP_VERBOSE_DEBUG
        use mod_adjacency_mat, only: adj_mat_from_conn

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
       
        if(verbose == OMMP_VERBOSE_DEBUG) then
            write(6, *) "Reading MMP file: ", input_file(1:len(trim(input_file)))
        end if

        ! open the (formatted) input file
        open (unit=iof_mmpinp, &
              file=input_file(1:len(trim(input_file))), &
              form='formatted', &
              access='sequential', &
              iostat=ist)

        if(ist /= 0) then
            call fatal_error('Error while opening MMP input file')
        end if

        if(verbose == OMMP_VERBOSE_DEBUG) then
            write(6, *) "Reading input parameters"
        end if

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
        
        if(verbose == OMMP_VERBOSE_DEBUG) then
            write(6, *) "Allocating memory"
        end if

        call mallocate('mmpol_init_from_mmp [my_cmm]', 3_ip, my_mm_atoms, my_cmm)
        call mallocate('mmpol_init_from_mmp [my_q]', my_ld_cart, my_mm_atoms, my_q)
        call mallocate('mmpol_init_from_mmp [my_pol]', my_mm_atoms, my_pol)
        call mallocate('mmpol_init_from_mmp [my_pol_atoms]', my_mm_atoms, pol_atoms_list)
        call mallocate('mmpol_init_from_mmp [my_ip11]', maxpgp, my_mm_atoms, my_ip11)
       
        call skip_lines(iof_mmpinp, my_mm_atoms+1) ! Skip a zero and the section of atomic numbers
        
        if(verbose == OMMP_VERBOSE_DEBUG) then
            write(6, *) "Reading coordinates"
        end if
        ! coordinates:
        do i = 1, my_mm_atoms
            read(iof_mmpinp,*) my_cmm(1:3,i)
        end do

        call skip_lines(iof_mmpinp, my_mm_atoms) ! Skip section of residues number

        if(verbose == OMMP_VERBOSE_DEBUG) then
            write(6, *) "Reading fixed multipoles"
        end if
        ! charges/multipoles:
        do i = 1, my_mm_atoms
            read(iof_mmpinp,*) my_q(1:my_ld_cart,i)
        end do

        if(verbose == OMMP_VERBOSE_DEBUG) then
            write(6, *) "Reading polarizabilities"
        end if
        ! polarizabilities:
        my_pol = 0.0_rp
        do i = 1, my_mm_atoms
            read(iof_mmpinp,*) my_pol(i)
        end do

        if(verbose == OMMP_VERBOSE_DEBUG) then
            write(6, *) "Processing polarizabilities"
        end if
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
        
        if(verbose == OMMP_VERBOSE_DEBUG) then
            write(6, *) "Initializing open mmpol module"
        end if
        ! mmpol module initialization
        call mmpol_init(my_ff_type, my_ff_rules, my_mm_atoms, my_pol_atoms)
        
        if(verbose == OMMP_VERBOSE_DEBUG) then
            write(6, *) "Converting input units to A.U."
        end if
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

        if(verbose == OMMP_VERBOSE_DEBUG) then
            write(6, *) "Processing connectivity informations"
        end if

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
     
        if(verbose == OMMP_VERBOSE_DEBUG) then
            write(6, *) "Populating utility arrays"
        end if
        
        call mmpol_prepare()
        
        if(verbose == OMMP_VERBOSE_DEBUG) then
            write(6, *) "Initialization from MMP file done."
        end if

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
        
        use mod_mmpol, only: verbose, cmm, q, pol
        use mod_mmpol, only: mm_atoms, amoeba, &
                             polar_mm, conn, mmat_polgrp
        use mod_mmpol, only: mol_frame, iz, ix, iy
        use mod_mmpol, only: fatal_error, mmpol_prepare, mmpol_init
        
        use mod_nonbonded, only: vdw_init
        use mod_memory, only: ip, rp, mfree, mallocate, memory_init
        use mod_io, only: mmpol_print_summary, iof_mmpinp
        use mod_constants, only: angstrom2au, OMMP_VERBOSE_DEBUG, OMMP_FF_AMOEBA, angstrom2au
        use mod_adjacency_mat, only: adj_mat_from_conn

        implicit none

        character(len=*), intent(in) :: xyz_file
        !! name of the input XYZ file
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file

        integer(ip), parameter :: iof_xyzinp = 200, &
                                  iof_prminp = 201, &
                                  maxn12 = 8
        integer(ip) :: my_mm_atoms, ist, i, j, l, atom_id, tokb, toke
        integer(ip), allocatable :: i12(:,:), attype(:)
        character(len=120) :: line, errstring

        
        if(verbose == OMMP_VERBOSE_DEBUG) then
            write(6, *) "Reading XYZ file: ", xyz_file(1:len(trim(xyz_file)))
        end if

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
        call mmpol_init(OMMP_FF_AMOEBA, 0_ip, my_mm_atoms, my_mm_atoms)        
        call vdw_init()
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
            if(scan(line(tokb:tokb), 'qwertyuiopasdfghjklzxcvbnmQWERTYUIOPASDFGHJKLZXCVBNM') == 0) then
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
        
        ! Writes the adjacency matrix in Yale sparse format in conn(1)
        call adj_mat_from_conn(i12, conn(1))
        call mfree('mmpol_init_from_xyz [i12]', i12)
        
        call read_prm(prm_file, attype)
        call mfree('mmpol_init_from_xyz [attype]', attype)
        
        call mmpol_prepare()

    end subroutine mmpol_init_from_xyz

    subroutine read_prm(prm_file, my_attype)

        use mod_memory, only: mallocate, mfree, ip, rp
        use mod_mmpol, only: fatal_error, pol, q, mol_frame, iz, ix, iy, &
                             mmat_polgrp, conn, mm_atoms
        use mod_nonbonded, only: vdw_e, vdw_f, vdw_r, vdw_set_pair
        use mod_constants, only: AMOEBA_ROT_NONE, &
                                 AMOEBA_ROT_Z_THEN_X, &
                                 AMOEBA_ROT_BISECTOR, &
                                 AMOEBA_ROT_Z_ONLY, &
                                 AMOEBA_ROT_Z_BISECT, &
                                 AMOEBA_ROT_3_FOLD, &
                                 angstrom2au, kcalmol2au
        
        implicit none
        
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file
        integer(ip), intent(in) :: my_attype(:)
        !! List of atom types that shoukd be used to populate parameter
        !! vectors

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: ist, i, j, k, l, iat, tokb, toke
        character(len=120) :: line, errstring
        integer(ip), allocatable :: polat(:), pgspec(:,:), multat(:), &
                                    multax(:,:), multframe(:), atclass(:), &
                                    vdwat(:), vdwpr_a(:), vdwpr_b(:)
        real(rp), allocatable :: thf(:), isopol(:), cmult(:,:), &
                                 vdw_e_prm(:), vdw_r_prm(:), vdw_f_prm(:), &
                                 vdwpr_r(:), vdwpr_e(:)
        integer(ip) :: natype, npolarize, ipolarize, nmult, imult, iax(3), &
                       nvdw, ivdw, atc, nvdwpr, ivdwpr
        logical :: ax_found(3), done


        ! open tinker xyz file
        open(unit=iof_prminp, &
             file=prm_file(1:len(trim(prm_file))), &
             form='formatted', &
             access='sequential', &
             iostat=ist)
        
        if(ist /= 0) then
           call fatal_error('Error while opening PRM input file')
        end if

        ! Read all the lines of file just to count how large vector should be 
        ! allocated 
        ist = 0
        natype = 0
        npolarize = 0
        nmult = 0
        nvdw = 0
        nvdwpr = 0
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            if(line(:5) == 'atom ') natype = natype + 1
            if(line(:9) == 'polarize ') npolarize = npolarize + 1
            if(line(:11) == 'multipole ') nmult = nmult + 1
            if(line(:4) == 'vdw ') nvdw = nvdw + 1
            if(line(:6) == 'vdwpr ') nvdwpr = nvdwpr + 1
        end do
        ! ATOM 
        call mallocate('read_prm [atclass]', natype, atclass)
        atclass = 0

        ! POLARIZE
        call mallocate('read_prm [polat]', npolarize, polat)
        call mallocate('read_prm [isopol]', npolarize, isopol)
        call mallocate('read_prm [thf]', npolarize, thf)
        call mallocate('read_prm [pgspec]', 8, npolarize, pgspec)
        pgspec = 0
        ipolarize = 1
        
        ! MULTIPOLE
        call mallocate('read_prm [multat]', nmult, multat)
        call mallocate('read_prm [multframe]', nmult, multframe)
        call mallocate('read_prm [multax]', 3, nmult, multax)
        call mallocate('read_prm [cmult]', 10, nmult, cmult)
        multax = 0
        imult = 1

        ! VDW
        call mallocate('read_prm [vdwat]', nvdw, vdwat)
        call mallocate('read_prm [vdw_r_prm]', nvdw, vdw_r_prm)
        call mallocate('read_prm [vdw_e_prm]', nvdw, vdw_e_prm)
        call mallocate('read_prm [vdw_f_prm]', nvdw, vdw_f_prm)
        vdw_f_prm = 1.0_rp
        ivdw = 1

        ! VDWPR
        call mallocate('read_prm [vdwpr_a]', nvdwpr, vdwpr_a)
        call mallocate('read_prm [vdwpr_b]', nvdwpr, vdwpr_b)
        call mallocate('read_prm [vdwpr_e]', nvdwpr, vdwpr_e)
        call mallocate('read_prm [vdwpr_r]', nvdwpr, vdwpr_r)
        ivdwpr = 1

        ! Restart the reading from the beginning to actually save the parameters
        rewind(iof_prminp)
        ist = 0
        i=1
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            
            if(line(:5) == 'atom ') then
                tokb = 6
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong ATOM card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) iat

                tokb = toke+1
                toke = tokenize(line, tokb)
                read(line(tokb:toke), *) atclass(iat)
                ! Only partial reading of ATOM card is needed for now.

            else if(line(:11) == 'multipole ') then
                tokb = 12 ! len of keyword + 1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong MULTIPOLE card"
                    call fatal_error(errstring)
                endif
                read(line(tokb:toke), *) multat(imult)
                if(multat(ipolarize) < 0) then
                    write(errstring, *) "Wrong MULTIPOLE card (specific atom not supported)"
                    call fatal_error(errstring)
                end if
                
                ! Rotation axis information
                tokb = toke+1
                toke = tokenize(line, tokb, 2)
                read(line(tokb:toke), *) multax(1:2,imult)
                ! TODO maybe some check here is not bad.

                ! For some centers also y axis is specified (integer) otherwise
                ! this is the zeroth-order multipole (charge)
                tokb = toke+1
                toke = tokenize(line, tokb)
                if(isint(line(tokb:toke))) then
                    ! Read the y axis
                    read(line(tokb:toke), *) multax(3,imult)
                    ! Get another token, this should be the charge.
                    tokb = toke+1 
                    toke = tokenize(line, tokb)
                end if

                ! The type of rotation is encoded in the sign/nullness of 
                ! of the axis specification
                if(multax(1,imult) == 0) then
                    multframe(imult) = AMOEBA_ROT_NONE
                else if(multax(2, imult) == 0) then
                    multframe(imult) = AMOEBA_ROT_Z_ONLY
                else if(multax(1,imult) < 0 .and. multax(2,imult) < 0 &
                        .and. multax(3,imult) < 0) then
                    multframe(imult) = AMOEBA_ROT_3_FOLD
                else if(multax(1,imult) < 0 .or. multax(2,imult) < 0) then
                    multframe(imult) = AMOEBA_ROT_BISECTOR
                else if(multax(2,imult) < 0 .and. multax(3,imult) < 0) then
                    multframe(imult) = AMOEBA_ROT_Z_BISECT
                else
                    multframe(imult) = AMOEBA_ROT_Z_THEN_X
                end if
                ! Remove the encoded information after saving it in a reasonable
                ! way
                multax(:,imult) = abs(multax(:,imult))

                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong MULTIPOLE card"
                    call fatal_error(errstring)
                end if

                read(line(tokb:toke), *) cmult(1, imult)

                read(iof_prminp, '(A)', iostat=ist) line
                read(line, *) cmult(2:4, imult)
                
                read(iof_prminp, '(A)', iostat=ist) line
                read(line, *) cmult(5, imult)
                
                read(iof_prminp, '(A)', iostat=ist) line
                read(line, *) cmult(6:7, imult)

                read(iof_prminp, '(A)', iostat=ist) line
                read(line, *) cmult(8:10, imult)

                imult = imult + 1
            else if(line(:9) == 'polarize ') then
                tokb = 10 ! len of keyword + 1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong POLARIZE card"
                    call fatal_error(errstring)
                endif
                read(line(tokb:toke), *) polat(ipolarize)
                if(polat(ipolarize) < 0) then
                    write(errstring, *) "Wrong POLARIZE card (specific atom not supported)"
                    call fatal_error(errstring)
                end if
                
                ! Isotropic polarizability
                tokb = toke+1
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong POLARIZE card"
                    call fatal_error(errstring)
                endif
                read(line(tokb:toke), *) isopol(ipolarize)

                ! Thole dumping factor
                tokb = toke+1
                toke = tokenize(line, tokb)
                if(isint(line(tokb:toke))) then
                    ! If this card is skipped then it is HIPPO
                    write(errstring, *) "HIPPO FF is not currently supported"
                    call fatal_error(errstring)
                else if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong POLARIZE card"
                    call fatal_error(errstring)
                endif
                read(line(tokb:toke), *) thf(ipolarize)
                
                ! Polarization group information
                tokb = toke+1
                toke = tokenize(line, tokb)
                if(isreal(line(tokb:toke))) then
                    ! If there is an additional real modifier then it is AMOEBA+
                    write(errstring, *) "AMOEBA+ FF is not currently supported"
                    call fatal_error(errstring)
                end if
                j = 1
                do while(toke > 0) 
                    if(.not. isint(line(tokb:toke))) then
                        write(errstring, *) "Wrong POLARIZE card"
                        call fatal_error(errstring)
                    endif
                    read(line(tokb:toke), *) pgspec(j,ipolarize)

                    tokb = toke+1
                    toke = tokenize(line, tokb)
                    j = j + 1
                end do
                ipolarize = ipolarize + 1
            else if(line(:4) == 'vdw ') then
                tokb = 5
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong VDW card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) vdwat(ivdw)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong VDW card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) vdw_r_prm(ivdw)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong VDW card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) vdw_e_prm(ivdw)
                
                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(toke > 0) then
                    if(.not. isreal(line(tokb:toke))) then
                        write(errstring, *) "Wrong VDW card"
                        call fatal_error(errstring)
                    end if
                    read(line(tokb:toke), *) vdw_f_prm(ivdw)
                endif

                ivdw = ivdw + 1
            else if(line(:6) == 'vdwpr ') then
                tokb = 7
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong VDWPR card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) vdwpr_a(ivdwpr)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong VDWPR card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) vdwpr_b(ivdwpr)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong VDWPR card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) vdwpr_r(ivdwpr)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong VDWPR card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) vdwpr_e(ivdwpr)

                ivdwpr = ivdwpr + 1
            end if
            i = i+1
        end do
        close(iof_prminp)
        
        mmat_polgrp = 0
        mol_frame = 0
        ix = 0
        iy = 0
        iz = 0

        do i=1, size(my_attype)
            ! Atom class for current atom
            atc = atclass(my_attype(i))
            
            ! Polarization
            do j=1, npolarize
                if(polat(j) == my_attype(i)) then
                    pol(i) = isopol(j) * angstrom2au**3
                    !TODO Thole factors.
                    ! Assign a polgroup label to each atom
                    if(mmat_polgrp(i) == 0) &
                        mmat_polgrp(i) = maxval(mmat_polgrp) + 1
                    
                    ! loop over the atoms connected to ith atom
                    do k=conn(1)%ri(i), conn(1)%ri(i+1)-1
                        iat = conn(1)%ci(k)
                        if(any(my_attype(iat) == pgspec(:,j))) then
                            ! The two atoms are in the same group
                            if(mmat_polgrp(iat) == 0) then
                                mmat_polgrp(iat) = mmat_polgrp(i)
                            else if(mmat_polgrp(iat) /= mmat_polgrp(i)) then
                                ! TODO This code have never been tested, as no
                                ! suitable case have been found
                                do l=1, mm_atoms
                                    if(mmat_polgrp(l) == 0) then
                                        continue
                                    else if(mmat_polgrp(l) == mmat_polgrp(iat) &
                                            .or. mmat_polgrp(l) == mmat_polgrp(i)) then
                                        mmat_polgrp(l) = min(mmat_polgrp(iat), mmat_polgrp(i))
                                    else if(mmat_polgrp(l) > max(mmat_polgrp(iat),mmat_polgrp(i))) then
                                        mmat_polgrp(l) = mmat_polgrp(l) - 1
                                    else
                                        continue
                                    end if
                                end do
                            end if
                        end if
                    end do
                end if
            end do

            ! Multipoles
            do j=1, nmult
                if(multat(j) == my_attype(i)) then
                    ! For each center different multipoles are defined for 
                    ! different environment. So first check if the environment
                    ! is the correct one
                    ! write(*, *) "ATOM", i, "MULTIPOLE", j
                    
                    ! Assignament with only 1,2-neighbours.
                    ax_found = .false.
                    iax = 0_ip

                    if(multframe(j) == AMOEBA_ROT_NONE) then
                        ! No axis needed
                        ax_found = .true.
                        ! if(all(ax_found)) write(*, *) ">>> NOAX"
                    else if(multframe(j) == AMOEBA_ROT_Z_ONLY) then
                        ! Assignament with only-z
                        ax_found(2:3) = .true.
                        do k=conn(1)%ri(i), conn(1)%ri(i+1)-1
                            iat = conn(1)%ci(k)
                            if(my_attype(iat) == multax(1,j) &
                               .and. .not. ax_found(1)) then
                                ax_found(1) = .true.
                                iax(1) = iat
                            end if
                        end do
                        ! if(all(ax_found)) write(*, *) ">>> ONLYZ"
                    else
                        ! 2 or 3 axis needed
                        if(multax(3,j) == 0) ax_found(3) = .true.
                        
                        ! Using only 1,2 connectivity
                        do k=conn(1)%ri(i), conn(1)%ri(i+1)-1
                            iat = conn(1)%ci(k)
                            if(my_attype(iat) == multax(1,j) &
                               .and. .not. ax_found(1)) then
                                ax_found(1) = .true.
                                iax(1) = iat
                            else if(my_attype(iat) == multax(2,j) &
                                    .and. .not. ax_found(2)) then
                                ax_found(2) = .true.
                                iax(2) = iat
                            else if(my_attype(iat) == multax(3,j) &
                                    .and. .not. ax_found(3)) then
                                ax_found(3) = .true.
                                iax(3) = iat
                            end if
                        end do
                        !if(all(ax_found)) write(*, *) ">>> ONLY12"

                        ! Using also 1,3 connectivity
                        if(ax_found(1) .and. .not. ax_found(2)) then
                            do k=conn(1)%ri(iax(1)), conn(1)%ri(iax(1)+1)-1
                                iat = conn(1)%ci(k)
                                if(iat == i .or. iat == iax(1)) cycle
                                if(my_attype(iat) == multax(2,j) &
                                   .and. .not. ax_found(2) &
                                   .and. iat /= iax(1)) then
                                    ax_found(2) = .true.
                                    iax(2) = iat
                                else if(my_attype(iat) == multax(3,j) &
                                        .and. .not. ax_found(3) & 
                                        .and. iat /= iax(1) &
                                        .and. iat /= iax(2)) then
                                    ax_found(3) = .true.
                                    iax(3) = iat
                                end if
                            end do
                            !if(all(ax_found)) write(*, *) ">>> ONLY13"
                        end if
                    end if

                    ! Everything is done, do the assignament
                    if(all(ax_found)) then
                        ix(i) = iax(2)
                        iy(i) = iax(3)
                        iz(i) = iax(1)
                        mol_frame(i) = multframe(j)
                        q(:,i) = cmult(:,j) 
                        exit
                    end if
                end if
            end do

            ! VdW parameters
            done = .false.
            do j=1, nvdw
                if(vdwat(j) == atc) then
                    done = .true.
                    vdw_e(i) = vdw_e_prm(j) * kcalmol2au
                    vdw_r(i) = vdw_r_prm(j) * angstrom2au
                    vdw_f(i) = vdw_f_prm(j)
                    exit
                end if
            end do
            if(.not. done) then
                call fatal_error("VdW parameter not found!")
            end if

            ! VdW pair parameters
            do l=1, nvdwpr
                if(vdwpr_a(l) == atc) then
                    do j=i+1, mm_atoms
                        if(atclass(my_attype(j)) == vdwpr_b(l)) then
                            call vdw_set_pair(i, j, &
                                              vdwpr_r(l) * angstrom2au, &
                                              vdwpr_e(l) * kcalmol2au)
                        end if
                    end do
                else if(vdwpr_b(l) == atc) then
                    do j=i+1, mm_atoms
                        if(atclass(my_attype(j)) == vdwpr_a(l)) then
                            call vdw_set_pair(i, j, &
                                              vdwpr_r(l) * angstrom2au, &
                                              vdwpr_e(l) * kcalmol2au)
                        end if
                    end do
                end if
            end do
        end do
    end subroutine read_prm

    function isint(s)
        implicit none

        character(len=*), intent(in) :: s
        logical :: isint
        
        isint = (verify(s, '+-1234567890') == 0)
        return
    end function

    function isreal(s)
        implicit none

        character(len=*), intent(in) :: s
        logical :: isreal

        isreal = (verify(s, '+-1234567890.') == 0)
        isreal = isreal .and. (scan(s, '.') /= 0)
        return
    end function

    function tokenize(s, ib, ntok)
        ! This function, given a string returns the first printable character
        ! if ib is absent or the last printable character after s(ib) if ib is 
        ! present. It's used to subdivide a string in tokens.

        use mod_memory, only: ip
        implicit none

        character(len=120), intent(in) :: s
        integer(ip), intent(inout), optional :: ib
        integer(ip), intent(in), optional :: ntok
        integer(ip) :: tokenize

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

end module mod_inputloader
