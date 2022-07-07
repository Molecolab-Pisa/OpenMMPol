module mod_prm
    !! This module handles the reading of a parameter file in .prm format and
    !! the asignament of parameters based on atom type and connectivity.
    use mod_memory, only: ip, rp
    use mod_utils, only: isreal, isint, tokenize

    implicit none
    private

    integer(ip), allocatable :: atclass(:)

    public :: assign_vdw, assign_pol, assign_mpoles

    contains 

    subroutine read_atom_cards(prm_file)
        use mod_memory, only: mallocate
        use mod_mmpol, only: fatal_error
        
        implicit none
        
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: ist, iat, toke, tokb 
        character(len=120) :: line, errstring
        integer(ip) :: natype


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
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            if(line(:5) == 'atom ') natype = natype + 1
        end do
        ! ATOM 
        call mallocate('read_prm [atclass]', natype, atclass)
        atclass = 0

        ! Restart the reading from the beginning to actually save the parameters
        rewind(iof_prminp)
        ist = 0
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
            end if
        end do
        close(iof_prminp)
    end subroutine read_atom_cards

    subroutine assign_vdw(prm_file, my_attype)
        use mod_memory, only: mallocate, mfree
        use mod_mmpol, only: fatal_error, mm_atoms
        use mod_nonbonded, only: vdw_e, vdw_f, vdw_r, vdw_set_pair, &
                                 vdw_screening, vdw_init
        use mod_constants, only: angstrom2au, kcalmol2au
        
        implicit none
        
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file
        integer(ip), intent(in) :: my_attype(:)
        !! List of atom types that shoukd be used to populate parameter
        !! vectors

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: ist, i, j, l, tokb, toke
        character(len=120) :: line, errstring
        integer(ip), allocatable :: vdwat(:), vdwpr_a(:), vdwpr_b(:)
        real(rp), allocatable :: vdw_e_prm(:), vdw_r_prm(:), vdw_f_prm(:), &
                                 vdwpr_r(:), vdwpr_e(:)
        integer(ip) :: nvdw, ivdw, atc, nvdwpr, ivdwpr
        logical :: done

        if(.not. allocated(atclass)) call read_atom_cards(prm_file)
        
        call vdw_init()

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
        nvdw = 0
        nvdwpr = 0
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            if(line(:4) == 'vdw ') nvdw = nvdw + 1
            if(line(:6) == 'vdwpr ') nvdwpr = nvdwpr + 1
        end do

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
           
            if(line(:13) == 'vdw-12-scale ') then
                tokb = 14
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong VDW-12-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) vdw_screening(1)
            
            else if(line(:13) == 'vdw-13-scale ') then
                tokb = 14
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong VDW-12-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) vdw_screening(2)
            
            else if(line(:13) == 'vdw-14-scale ') then
                tokb = 14
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong VDW-12-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) vdw_screening(3)
            
            else if(line(:13) == 'vdw-15-scale ') then
                tokb = 14
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong VDW-12-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) vdw_screening(4)

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
        
        do i=1, size(my_attype)
            ! Atom class for current atom
            atc = atclass(my_attype(i))
            
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
        
        call mfree('read_prm [vdwat]', vdwat)
        call mfree('read_prm [vdw_r_prm]', vdw_r_prm)
        call mfree('read_prm [vdw_e_prm]', vdw_e_prm)
        call mfree('read_prm [vdw_f_prm]', vdw_f_prm)
        call mfree('read_prm [vdwpr_a]', vdwpr_a)
        call mfree('read_prm [vdwpr_b]', vdwpr_b)
        call mfree('read_prm [vdwpr_e]', vdwpr_e)
        call mfree('read_prm [vdwpr_r]', vdwpr_r)
    
    end subroutine assign_vdw
    
    subroutine assign_pol(prm_file, my_attype)
        use mod_memory, only: mallocate, mfree, ip, rp
        use mod_mmpol, only: fatal_error, pol, mmat_polgrp, conn, mm_atoms
        use mod_mmpol, only: mscale, set_screening_parameters
        use mod_constants, only: angstrom2au
        
        implicit none
        
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file
        integer(ip), intent(in) :: my_attype(:)
        !! List of atom types that shoukd be used to populate parameter
        !! vectors

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: ist, i, j, k, l, iat, tokb, toke
        character(len=120) :: line, errstring
        
        integer(ip), allocatable :: polat(:), pgspec(:,:) 
        real(rp), allocatable :: thf(:), isopol(:)
        real(rp) :: usc(4), psc(4), pisc(4), dsc(4)

        integer(ip) :: npolarize, ipolarize


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
        npolarize = 0
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            if(line(:9) == 'polarize ') npolarize = npolarize + 1
        end do
        
        call mallocate('read_prm [polat]', npolarize, polat)
        call mallocate('read_prm [isopol]', npolarize, isopol)
        call mallocate('read_prm [thf]', npolarize, thf)
        call mallocate('read_prm [pgspec]', 8, npolarize, pgspec)
        pgspec = 0
        ipolarize = 1
        
        ! Restart the reading from the beginning to actually save the parameters
        rewind(iof_prminp)
        ist = 0
        i=1
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
           
            if(line(:15) == 'polar-12-intra ') then
                tokb = 16
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong POLAR-12-INTRA card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) pisc(1)

            else if(line(:15) == 'polar-13-intra ') then
                tokb = 16
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong POLAR-13-INTRA card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) pisc(2)

            else if(line(:15) == 'polar-14-intra ') then
                tokb = 16
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong POLAR-14-INTRA card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) pisc(3)

            else if(line(:15) == 'polar-15-intra ') then
                tokb = 16
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong POLAR-15-INTRA card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) pisc(4)

            else if(line(:15) == 'polar-12-scale ') then
                tokb = 16
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong POLAR-12-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) psc(1)

            else if(line(:15) == 'polar-13-scale ') then
                tokb = 16
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong POLAR-13-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) psc(2)

            else if(line(:15) == 'polar-14-scale ') then
                tokb = 16
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong POLAR-14-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) psc(3)

            else if(line(:15) == 'polar-15-scale ') then
                tokb = 16
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong POLAR-15-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) psc(4)

            else if(line(:16) == 'direct-11-scale ') then
                tokb = 17
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong DIRECT-11-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) dsc(1)

            else if(line(:16) == 'direct-12-scale ') then
                tokb = 17
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong DIRECT-12-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) dsc(2)

            else if(line(:16) == 'direct-13-scale ') then
                tokb = 17
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong DIRECT-13-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) dsc(3)

            else if(line(:16) == 'direct-14-scale ') then
                tokb = 17
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong DIRECT-14-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) dsc(4)
            
            else if(line(:16) == 'mutual-11-scale ') then
                tokb = 17
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong MUTUAL-11-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) usc(1)
            else if(line(:16) == 'mutual-12-scale ') then
                tokb = 17
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong MUTUAL-12-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) usc(2)

            else if(line(:16) == 'mutual-13-scale ') then
                tokb = 17
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong MUTUAL-13-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) usc(3)

            else if(line(:16) == 'mutual-14-scale ') then
                tokb = 17
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong MUTUAL-14-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) usc(4)
            
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
            end if
            i = i+1
        end do
        close(iof_prminp)
        
        call set_screening_parameters(mscale, psc, dsc, usc, pisc)

        ! Now assign the parameters to the atoms
        mmat_polgrp = 0
        do i=1, size(my_attype)
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
        end do
        
        call mfree('read_prm [polat]', polat)
        call mfree('read_prm [isopol]', isopol)
        call mfree('read_prm [thf]', thf)
        call mfree('read_prm [pgspec]', pgspec)
    
    end subroutine assign_pol
    
    subroutine assign_mpoles(prm_file, my_attype)
        use mod_memory, only: mallocate, mfree
        use mod_mmpol, only: fatal_error, q, mol_frame, iz, ix, iy, &
                             conn
        use mod_mmpol, only: pscale, dscale, uscale, pscale_intra, &
                             set_screening_parameters
        use mod_constants, only: AMOEBA_ROT_NONE, &
                                 AMOEBA_ROT_Z_THEN_X, &
                                 AMOEBA_ROT_BISECTOR, &
                                 AMOEBA_ROT_Z_ONLY, &
                                 AMOEBA_ROT_Z_BISECT, &
                                 AMOEBA_ROT_3_FOLD
        
        implicit none
        
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file
        integer(ip), intent(in) :: my_attype(:)
        !! List of atom types that shoukd be used to populate parameter
        !! vectors

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: ist, i, j, k, iat, tokb, toke
        character(len=120) :: line, errstring
        integer(ip), allocatable :: multat(:), multax(:,:), multframe(:)
        real(rp), allocatable :: cmult(:,:)
        real(rp) :: msc(4)
        integer(ip) :: nmult, imult, iax(3)
        logical :: ax_found(3), found13, only12


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
        nmult = 0
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            if(line(:11) == 'multipole ') nmult = nmult + 1
        end do
        
        ! MULTIPOLE
        call mallocate('read_prm [multat]', nmult, multat)
        call mallocate('read_prm [multframe]', nmult, multframe)
        call mallocate('read_prm [multax]', 3, nmult, multax)
        call mallocate('read_prm [cmult]', 10, nmult, cmult)
        multax = 0
        imult = 1

        ! Restart the reading from the beginning to actually save the parameters
        rewind(iof_prminp)
        ist = 0
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            
            if(line(:15) == 'mpole-12-scale ') then
                tokb = 16
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong MPOLE-12-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) msc(1)
            
            else if(line(:15) == 'mpole-13-scale ') then
                tokb = 16
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong MPOLE-13-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) msc(2)
            
            else if(line(:15) == 'mpole-14-scale ') then
                tokb = 16
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong MPOLE-14-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) msc(3)
            
            else if(line(:15) == 'mpole-15-scale ') then
                tokb = 16
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong MPOLE-15-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) msc(4)

            else if(line(:11) == 'multipole ') then
                tokb = 12 ! len of keyword + 1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong MULTIPOLE card"
                    call fatal_error(errstring)
                endif
                read(line(tokb:toke), *) multat(imult)
                if(multat(imult) < 0) then
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
            end if
        end do
        close(iof_prminp)
        
        call set_screening_parameters(msc, pscale, dscale, uscale, pscale_intra)

        mol_frame = 0
        ix = 0
        iy = 0
        iz = 0

        do i=1, size(my_attype)
            ! Multipoles
            only12 = .false. ! Only search for params based on 12 connectivity

            do j=1, nmult
                found13 = .false. ! Parameter found is based on 13 connectivity
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
                            if(all(ax_found)) found13 = .true.
                        end if
                    end if

                    ! Everything is done, no further improvement is possible
                    if(all(ax_found) .and. .not. (only12 .and. found13)) then
                        ix(i) = iax(2)
                        iy(i) = iax(3)
                        iz(i) = iax(1)
                        mol_frame(i) = multframe(j)
                        q(:,i) = cmult(:,j) 
                        if(.not. found13) then
                            exit ! No further improvement is possible
                        else
                            only12 = .true.
                        end if
                    end if
                end if
            end do
        end do
        
        call mfree('read_prm [multat]', multat)
        call mfree('read_prm [multframe]', multframe)
        call mfree('read_prm [multax]', multax)
        call mfree('read_prm [cmult]', cmult)
    
    end subroutine assign_mpoles

end module mod_prm
