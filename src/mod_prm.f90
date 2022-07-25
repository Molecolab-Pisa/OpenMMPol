module mod_prm
    !! This module handles the reading of a parameter file in .prm format and
    !! the asignament of parameters based on atom type and connectivity.

    use mod_memory, only: ip, rp
    use mod_utils, only: isreal, isint, tokenize, count_substr_occurence, &
                         str_to_lower

    implicit none
    private

    integer(ip), allocatable :: atclass(:), atz(:)

    public :: assign_vdw, assign_pol, assign_mpoles, assign_bond, &
              assign_angle, assign_urey, assign_strbnd, assign_opb, &
              assign_pitors, assign_torsion, assign_tortors, terminate_prm
    public :: check_keyword

    contains 
    
#include "prm_keywords.f90"

    subroutine terminate_prm()
        use mod_memory, only: mfree
        
        implicit none
        
        if(allocated(atclass)) call mfree('terminate_prm [atclass]', atclass)
        if(allocated(atz)) call mfree('terminate_prm [atz]', atz)

    end subroutine terminate_prm
    
    subroutine read_atom_cards(prm_file)
        use mod_memory, only: mallocate
        use mod_mmpol, only: fatal_error
        
        implicit none
        
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: ist, iat, toke, tokb, tokb1, nquote
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
            line = str_to_lower(line)
            if(line(:5) == 'atom ') natype = natype + 1
        end do
        ! ATOM 
        call mallocate('read_prm [atclass]', natype, atclass)
        call mallocate('read_prm [atz]', natype, atz)
        atclass = 0
        atz = 0

        ! Restart the reading from the beginning to actually save the parameters
        rewind(iof_prminp)
        ist = 0
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            line = str_to_lower(line)
           
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

                tokb = toke+1
                toke = tokenize(line, tokb)
                ! This token contain the atom name

                tokb = toke+1
                toke = tokenize(line, tokb)
                nquote = count_substr_occurence(line(tokb:toke), '"')
                do while(nquote < 2)
                    tokb1 = toke+1
                    toke = tokenize(line, tokb1)
                    nquote = nquote + count_substr_occurence(line(tokb1:toke), '"')
                end do
                ! This token contains the description string
                tokb = toke+1
                toke = tokenize(line, tokb)
                read(line(tokb:toke), *) atz(iat)

                ! Only partial reading of ATOM card is needed for now.
            end if
        end do
        close(iof_prminp)
    end subroutine read_atom_cards

    subroutine assign_bond(prm_file, my_attype)
        use mod_memory, only: mallocate, mfree
        use mod_mmpol, only: fatal_error, mm_atoms, conn
        use mod_bonded, only: bond_init, bond_potential, bondat, &
                              kbond, l0bond, bond_cubic, bond_quartic
        use mod_constants, only: angstrom2au, kcalmol2au
        
        implicit none
        
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file
        integer(ip), intent(in) :: my_attype(:)
        !! List of atom types that shoukd be used to populate parameter
        !! vectors

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: ist, i, j, l, jat, tokb, toke, ibnd, nbnd, &
                       cla, clb
        character(len=120) :: line, errstring
        integer(ip), allocatable :: classa(:), classb(:)
        real(rp), allocatable :: kbnd(:), l0bnd(:)
        logical :: done

        if(.not. allocated(atclass)) call read_atom_cards(prm_file)
        
        ! We assume that all pair of bonded atoms have a bonded 
        ! parameter
        call bond_init((conn(1)%ri(mm_atoms+1)-1) / 2)
        kbond = 0
        l0bond = 0
        l=1
        do i=1, mm_atoms
            do j=conn(1)%ri(i), conn(1)%ri(i+1)-1
                jat = conn(1)%ci(j)
                if(i < jat) then
                    bondat(1,l) = i
                    bondat(2,l) = jat
                    l = l+1
                end if
            end do
        end do

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
        nbnd = 0
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            line = str_to_lower(line)
            if(line(:5) == 'bond ') nbnd = nbnd + 1
        end do

        call mallocate('assign_bond [classa]', nbnd, classa)
        call mallocate('assign_bond [classb]', nbnd, classb)
        call mallocate('assign_bond [l0bnd]', nbnd, l0bnd)
        call mallocate('assign_bond [kbnd]', nbnd, kbnd)

        ! Restart the reading from the beginning to actually save the parameters
        rewind(iof_prminp)
        ist = 0
        ibnd = 1
        i=1
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            line = str_to_lower(line)
           
            if(line(:11) == 'bond-cubic ') then
                tokb = 12
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong BOND-CUBIC card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) bond_cubic
                ! This parameter is 1/Angstrom
                bond_cubic = bond_cubic / angstrom2au
            
            else if(line(:13) == 'bond-quartic ') then
                tokb = 14
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong BOND-QUARTIC card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) bond_quartic
                bond_quartic = bond_quartic / (angstrom2au**2)
            
            else if(line(:5) == 'bond ') then
                tokb = 6
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong BOND card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classa(ibnd)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong BOND card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classb(ibnd)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong BOND card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) kbnd(ibnd)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong BOND card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) l0bnd(ibnd)
                
                ibnd = ibnd + 1
            end if
            i = i+1
        end do
        close(iof_prminp)
        
        do i=1, size(bondat,2)
            ! Atom class for current pair
            cla = atclass(my_attype(bondat(1,i)))
            clb = atclass(my_attype(bondat(2,i)))
            
            done = .false.
            do j=1, nbnd
                if((classa(j)==cla .and. classb(j)==clb) .or. &
                   (classa(j)==clb .and. classb(j)==cla)) then
                    done = .true.
                    kbond(i) = kbnd(j) * kcalmol2au / (angstrom2au**2)
                    l0bond(i) = l0bnd(j) * angstrom2au
                    exit
                end if
            end do
            if(.not. done) then
                call fatal_error("Bond parameter not found!")
            end if
        end do
        
        call mfree('assign_bond [classa]', classa)
        call mfree('assign_bond [classb]', classb)
        call mfree('assign_bond [l0bnd]', l0bnd)
        call mfree('assign_bond [kbnd]', kbnd)
    
    end subroutine assign_bond
    
    subroutine assign_urey(prm_file, my_attype)
        use mod_memory, only: mallocate, mfree
        use mod_mmpol, only: fatal_error, mm_atoms, conn
        use mod_bonded, only: urey_init, urey_potential, ureyat, &
                              kurey, l0urey, urey_cubic, urey_quartic
        use mod_constants, only: angstrom2au, kcalmol2au
        
        implicit none
        
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file
        integer(ip), intent(in) :: my_attype(:)
        !! List of atom types that shoukd be used to populate parameter
        !! vectors

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: ist, i, j, tokb, toke, iub, nub, &
                       cla, clb, clc, maxub, a, b, c, jc, jb 
        character(len=120) :: line, errstring
        integer(ip), allocatable :: classa(:), classb(:), classc(:), ubtmp(:)
        real(rp), allocatable :: kub(:), l0ub(:)
        logical :: done

        if(.not. allocated(atclass)) call read_atom_cards(prm_file)
        
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
        nub = 0
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            line = str_to_lower(line)
            if(line(:9) == 'ureybrad ') nub = nub + 1
        end do

        maxub = conn(2)%ri(mm_atoms+1)-1 
        ! Maximum number of UB terms (each angle have an UB term)
        call mallocate('assign_urey [classa]', nub, classa)
        call mallocate('assign_urey [classb]', nub, classb)
        call mallocate('assign_urey [classc]', nub, classc)
        call mallocate('assign_urey [l0ub]', nub, l0ub)
        call mallocate('assign_urey [kub]', nub, kub)
        call mallocate('assign_urey [ubtmp]', maxub, ubtmp)

        ! Restart the reading from the beginning to actually save the parameters
        rewind(iof_prminp)
        ist = 0
        iub = 1
        i=1
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            line = str_to_lower(line)
           
            if(line(:11) == 'urey-cubic ') then
                tokb = 12
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong UREY-CUBIC card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) urey_cubic
                ! This parameter is 1/Angstrom
                urey_cubic = urey_cubic / angstrom2au
            
            else if(line(:13) == 'urey-quartic ') then
                tokb = 14
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong UREY-QUARTIC card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) urey_quartic
                urey_quartic = urey_quartic / (angstrom2au**2)
            
            else if(line(:9) == 'ureybrad ') then
                tokb = 10
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong UREYBRAD card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classa(iub)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong UREYBRAD card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classb(iub)
                
                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong UREYBRAD card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classc(iub)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong UREYBRAD card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) kub(iub)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong UREYBRAD card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) l0ub(iub)
                
                iub = iub + 1
            end if
            i = i+1
        end do
        close(iof_prminp)
        
        ubtmp = -1
        do a=1, mm_atoms
            cla = atclass(my_attype(a))
            do jb=conn(2)%ri(a), conn(2)%ri(a+1)-1
                b = conn(2)%ci(jb)
                done = .false.
                if(a > b) cycle
                clb = atclass(my_attype(b))
                
                do jc=conn(1)%ri(a), conn(1)%ri(a+1)-1
                    c = conn(1)%ci(jc)
                    if(all(conn(1)%ci(conn(1)%ri(b):conn(1)%ri(b+1)-1) /= c)) cycle
                    ! There is an angle in the form A-C-B
                    clc = atclass(my_attype(c))
                    do j=1, nub
                        if((cla == classa(j) &
                            .and. clb == classc(j) &
                            .and. clc == classb(j)) .or. &
                           (clb == classa(j) &
                            .and. cla == classc(j) &
                            .and. clc == classb(j))) then

                            ubtmp(jb) = j 
                            ! Temporary assignament in a sparse matrix logic
                            done = .true.
                            exit
                        end if
                    end do
                        
                    if(done) exit 
                    ! If we have already found a parameter for A-B pair, stop 
                    ! the research
                end do
            end do
        end do

        call urey_init(count(ubtmp > 0))
        iub = 1
        do a=1, mm_atoms
            do jb=conn(2)%ri(a), conn(2)%ri(a+1)-1
                if(ubtmp(jb) > 0) then
                    ureyat(1,iub) = a
                    ureyat(2,iub) = conn(2)%ci(jb)
                    kurey(iub) = kub(ubtmp(jb)) * kcalmol2au / (angstrom2au**2) 
                    l0urey(iub) = l0ub(ubtmp(jb)) * angstrom2au
                    iub = iub + 1
                end if
            end do
        end do

        call mfree('assign_urey [classa]', classa)
        call mfree('assign_urey [classb]', classb)
        call mfree('assign_urey [classc]', classc)
        call mfree('assign_urey [l0ub]', l0ub)
        call mfree('assign_urey [kub]', kub)
        call mfree('assign_urey [ubtmp]', ubtmp)
        
    end subroutine assign_urey
    
    subroutine assign_strbnd(prm_file, my_attype)
        use mod_memory, only: mallocate, mfree
        use mod_mmpol, only: fatal_error, mm_atoms, conn
        use mod_bonded, only: strbnd_init, strbndat, strbndk1, strbndk2, &
                              strbndl10, strbndl20, strbndthet0, &
                              bondat, l0bond, angleat, eqangle 

        use mod_constants, only: kcalmol2au, angstrom2au
        
        implicit none
        
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file
        integer(ip), intent(in) :: my_attype(:)
        !! List of atom types that shoukd be used to populate parameter
        !! vectors

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: ist, i, j, tokb, toke, isb, nstrbnd, &
                       cla, clb, clc, a, b, c, jc, jb, maxsb, &
                       l1a, l1b, l2a, l2b
        character(len=120) :: line, errstring
        integer(ip), allocatable :: classa(:), classb(:), classc(:), sbtmp(:), &
                                    sbattmp(:, :)
        real(rp), allocatable :: k1(:), k2(:)
        logical :: done, thet_done, l1_done, l2_done

        if(.not. allocated(atclass)) call read_atom_cards(prm_file)
        
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
        nstrbnd = 0
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            line = str_to_lower(line)
            if(line(:7) == 'strbnd ') nstrbnd = nstrbnd + 1
        end do

        maxsb = (conn(2)%ri(mm_atoms+1)-1) / 2
        call mallocate('assign_angle [classa]', nstrbnd, classa)
        call mallocate('assign_angle [classb]', nstrbnd, classb)
        call mallocate('assign_angle [classc]', nstrbnd, classc)
        call mallocate('assign_angle [eqang]', nstrbnd, k1)
        call mallocate('assign_angle [kang]', nstrbnd, k2)
        call mallocate('assign_angle [sbtmp]', maxsb, sbtmp)
        call mallocate('assign_angle [sbattmp]', 3, maxsb, sbattmp)

        ! Restart the reading from the beginning to actually save the parameters
        rewind(iof_prminp)
        ist = 0
        isb = 1
        i=1
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            line = str_to_lower(line)
           
            if(line(:7) == 'strbnd ') then
                tokb = 8
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong STRBND card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classa(isb)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong STRBND card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classb(isb)
                
                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong STRBND card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classc(isb)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong STRBND card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) k1(isb)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong STRBND card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) k2(isb)

                isb = isb + 1
            end if
            i = i+1
        end do
        close(iof_prminp)
        
        isb = 1
        do a=1, mm_atoms
            cla = atclass(my_attype(a))
            do jb=conn(2)%ri(a), conn(2)%ri(a+1)-1
                b = conn(2)%ci(jb)
                if(a > b) cycle
                clb = atclass(my_attype(b))
                
                do jc=conn(1)%ri(a), conn(1)%ri(a+1)-1
                    c = conn(1)%ci(jc)
                    if(all(conn(1)%ci(conn(1)%ri(b):conn(1)%ri(b+1)-1) /= c)) cycle
                    ! There is an angle in the form A-C-B
                    clc = atclass(my_attype(c))
                    done = .false.

                    do j=1, nstrbnd
                        if((cla == classa(j) &
                            .and. clb == classc(j) &
                            .and. clc == classb(j)) .or. &
                           (clb == classa(j) &
                            .and. cla == classc(j) &
                            .and. clc == classb(j))) then
                            sbattmp(1,isb) = a
                            sbattmp(2,isb) = c
                            sbattmp(3,isb) = b
                            if(cla == classa(j)) then
                                ! Assign the correct k to each bond stretching!
                                sbtmp(isb) = j
                            else
                                sbtmp(isb) = -j
                            end if
                            isb = isb + 1
                            exit
                        end if
                    end do
                end do
            end do
        end do

        call strbnd_init(isb-1)

        do i=1, isb-1
            ! First assign the parameters
            strbndat(:,i) = sbattmp(:,i)
            j = abs(sbtmp(i))
            if(sbtmp(i) > 0) then
                strbndk1(i) = k1(j) * kcalmol2au / angstrom2au
                strbndk2(i) = k2(j) * kcalmol2au / angstrom2au
                l1a = min(strbndat(1,i), strbndat(2,i))
                l1b = max(strbndat(1,i), strbndat(2,i))
                l2a = min(strbndat(3,i), strbndat(2,i))
                l2b = max(strbndat(3,i), strbndat(2,i))
            else
                strbndk1(i) = k2(j) * kcalmol2au / angstrom2au
                strbndk2(i) = k1(j) * kcalmol2au / angstrom2au
                l1a = min(strbndat(3,i), strbndat(2,i))
                l1b = max(strbndat(3,i), strbndat(2,i))
                l2a = min(strbndat(1,i), strbndat(2,i))
                l2b = max(strbndat(1,i), strbndat(2,i))
            end if
           
            ! Now search for the corresponding bond and angle parameters to
            ! set the equilibrium distances and angle
            l1_done = .false.
            l2_done = .false.
            thet_done = .false.
            

            do j=1, size(bondat, 2)
                if(l1a == bondat(1,j) .and. l1b == bondat(2,j)) then
                    l1_done = .true.
                    strbndl10(i) = l0bond(j)
                else if(l2a == bondat(1,j) .and. l2b == bondat(2,j)) then
                    l2_done = .true.
                    strbndl20(i) = l0bond(j)
                end if

                if(l1_done .and. l2_done) exit
            end do
            
            do j=1, size(angleat, 2)
                if(all(strbndat(:,i) == angleat(:,j))) then
                    thet_done = .true.
                    strbndthet0(i) = eqangle(j)
                    exit
                end if
            end do
            
            if(.not.(l1_done .and. l2_done .and. thet_done)) then
                call fatal_error("Ill-defined stretch-bending cross term")
            end if
        end do

        call mfree('assign_strbnd [classa]', classa)
        call mfree('assign_strbnd [classb]', classb)
        call mfree('assign_strbnd [classc]', classc)
        call mfree('assign_strbnd [eqang]', k1)
        call mfree('assign_strbnd [kang]', k2)
        call mfree('assign_strbnd [sbtmp]', sbtmp)
        call mfree('assign_strbnd [sbattmp]', sbattmp)
    
    end subroutine assign_strbnd
    
    subroutine assign_opb(prm_file, my_attype)
        use mod_memory, only: mallocate, mfree
        use mod_mmpol, only: fatal_error, mm_atoms, conn
        use mod_bonded, only: opb_cubic, opb_quartic, &
                              opb_pentic, opb_sextic, opbat, kopb, opb_init

        use mod_constants, only: kcalmol2au, rad2deg
        
        implicit none
        
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file
        integer(ip), intent(in) :: my_attype(:)
        !! List of atom types that shoukd be used to populate parameter
        !! vectors

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: ist, i, tokb, toke, iopb, nopb, &
                       cla, clb, clc, cld, maxopb, a, b, c, d, jc, jb, iprm
        character(len=120) :: line, errstring, opb_type
        integer(ip), allocatable :: classa(:), classb(:), classc(:), & 
                                    classd(:), tmpat(:,:)
        real(rp), allocatable :: kopbend(:), tmpk(:)

        if(.not. allocated(atclass)) call read_atom_cards(prm_file)
        
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
        nopb = 0
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            line = str_to_lower(line)
            if(line(:7) == 'opbend ') nopb = nopb + 1
        end do

        maxopb = mm_atoms*3
        call mallocate('assign_opb [classa]', nopb, classa)
        call mallocate('assign_opb [classb]', nopb, classb)
        call mallocate('assign_opb [classc]', nopb, classc)
        call mallocate('assign_opb [classd]', nopb, classd)
        call mallocate('assign_opb [kopbend]', nopb, kopbend)
        call mallocate('assign_opb [tmpat]', 4, maxopb, tmpat)
        call mallocate('assign_opb [tmpk]', maxopb, tmpk)

        ! Restart the reading from the beginning to actually save the parameters
        rewind(iof_prminp)
        ist = 0
        iopb = 1
        i=1
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            line = str_to_lower(line)
           
            if(line(:11) == 'opbendtype ') then
                tokb = 12
                toke = tokenize(line, tokb)
                read(line(tokb:toke), *) opb_type
            
            else if(line(:13) == 'opbend-cubic ') then
                tokb = 14
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong OPBEND-CUBIC card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) opb_cubic
                opb_cubic = opb_cubic * rad2deg
            
            else if(line(:15) == 'opbend-quartic ') then
                tokb = 16
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong OPBEND-QUARTIC card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) opb_quartic
                opb_quartic = opb_quartic * rad2deg**2
            
            else if(line(:14) == 'opbend-pentic ') then
                tokb = 15
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong OPBEND-PENTIC card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) opb_pentic
                opb_pentic = opb_pentic * rad2deg**3
            
            else if(line(:14) == 'opbend-sextic ') then
                tokb = 15
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong OPBEND-SEXTIC card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) opb_sextic
                opb_sextic = opb_sextic * rad2deg**4
            
            else if(line(:7) == 'opbend ') then
                tokb = 8
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong OPBEND card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classa(iopb)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong OPBEND card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classb(iopb)
                
                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong OPBEND card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classc(iopb)
                
                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong OPBEND card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classd(iopb)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong OPBEND card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) kopbend(iopb)
                
                iopb = iopb + 1
            end if
            i = i+1
        end do
        close(iof_prminp)
       
        iopb = 1
        do a=1, mm_atoms
            ! Check if the center is trigonal
            if(conn(1)%ri(a+1) - conn(1)%ri(a) /= 3) cycle
            cla = atclass(my_attype(a))
            ! Loop over the atoms connected to the trignonal center
            do jb=conn(1)%ri(a), conn(1)%ri(a+1)-1
                b = conn(1)%ci(jb)
                clb = atclass(my_attype(b))
              
                c = -1
                d = -1
                clc = 0
                cld = 0
                do jc=conn(1)%ri(a), conn(1)%ri(a+1)-1
                    if(conn(1)%ci(jc) == b) cycle
                    if(c < 0) then
                        c = conn(1)%ci(jc)
                        clc = atclass(my_attype(c))
                    else if(d < 0) then
                        d = conn(1)%ci(jc)
                        cld = atclass(my_attype(d))
                    end if
                end do

                do iprm=1, nopb
                    if((classa(iprm) == clb) .and. &
                       (classb(iprm) == cla) .and. &
                       ((classc(iprm) == clc .and. & 
                         classd(iprm) == cld) .or. &
                        (classd(iprm) == clc .and. &
                         classc(iprm) == cld) .or. &
                        (classd(iprm) == 0 .and. &
                         (classc(iprm) == cld .or. classc(iprm) == clc)) .or. &
                        (classc(iprm) == 0 .and. &
                         (classd(iprm) == cld .or. classd(iprm) == clc)) .or. &
                        (classc(iprm) == 0 .or. classd(iprm) == 0))) then
                        ! The parameter is ok
                        tmpat(1,iopb) = a
                        tmpat(2,iopb) = b
                        tmpat(3,iopb) = c
                        tmpat(4,iopb) = d
                        tmpk(iopb) = kopbend(iprm)
                        iopb = iopb + 1
                        exit
                    endif
                end do
            end do
        end do

        call opb_init(iopb-1, trim(opb_type))
        
        do i=1, iopb-1
            kopb(i) = tmpk(i) * kcalmol2au
            opbat(:,i) = tmpat(:,i)
        end do

        call mfree('assign_opb [classa]', classa)
        call mfree('assign_opb [classb]', classb)
        call mfree('assign_opb [classc]', classc)
        call mfree('assign_opb [classd]', classd)
        call mfree('assign_opb [kopbend]', kopbend)
        call mfree('assign_opb [tmpat]', tmpat)
        call mfree('assign_opb [tmpk]', tmpk)
    
    end subroutine assign_opb
    
    subroutine assign_pitors(prm_file, my_attype)
        use mod_memory, only: mallocate, mfree
        use mod_mmpol, only: fatal_error, mm_atoms, conn
        use mod_bonded, only: kpitors, pitorsat, pitors_init
        use mod_constants, only: kcalmol2au
        
        implicit none
        
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file
        integer(ip), intent(in) :: my_attype(:)
        !! List of atom types that shoukd be used to populate parameter
        !! vectors

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: ist, i, tokb, toke, ipitors, npitors, &
                       cla, clb, maxpi, a, b, c, jb, iprm
        character(len=120) :: line, errstring
        integer(ip), allocatable :: classa(:), classb(:), tmpat(:,:)
        real(rp), allocatable :: kpi(:), tmpk(:)

        if(.not. allocated(atclass)) call read_atom_cards(prm_file)
        
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
        npitors = 1
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            line = str_to_lower(line)
            if(line(:7) == 'pitors ') npitors = npitors + 1
        end do

        maxpi = mm_atoms 
        ! TODO This is maybe excessive, all trivalent atomso should be enough
        call mallocate('assign_pitors [classa]', npitors, classa)
        call mallocate('assign_pitors [classb]', npitors, classb)
        call mallocate('assign_pitors [kpi]', npitors, kpi)
        call mallocate('assign_pitors [tmpat]', 6, maxpi, tmpat)
        call mallocate('assign_pitors [tmpk]', maxpi, tmpk)

        ! Restart the reading from the beginning to actually save the parameters
        rewind(iof_prminp)
        ist = 0
        ipitors = 1
        i=1
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            line = str_to_lower(line)
           
            if(line(:7) == 'pitors ') then
                tokb = 8
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong PITORS card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classa(ipitors)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong PITORS card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classb(ipitors)
                
                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong PITORS card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) kpi(ipitors)
                
                ipitors = ipitors + 1
            end if
            i = i+1
        end do
        close(iof_prminp)
       
        ipitors = 1
        do a=1, mm_atoms
            ! Check if the center is trigonal
            if(conn(1)%ri(a+1) - conn(1)%ri(a) /= 3) cycle
            cla = atclass(my_attype(a))
            ! Loop over the atoms connected to the trignonal center
            do jb=conn(1)%ri(a), conn(1)%ri(a+1)-1
                b = conn(1)%ci(jb)
                ! This avoid to compute the functions two times
                if(a > b) cycle
                
                !Check if the second center is trigonal
                if(conn(1)%ri(b+1) - conn(1)%ri(b) /= 3) cycle
                clb = atclass(my_attype(b))

                do iprm=1, npitors
                    if((cla == classa(iprm) .and. clb == classb(iprm)) .or. &
                       (clb == classa(iprm) .and. cla == classb(iprm))) then
                        ! The parameter is the right one
                        ! Save the atoms in the following way:
                        !
                        !  2        5            a => 1
                        !   \      /             b => 4
                        !    1 -- 4  
                        !   /      \
                        !  3        6
                        
                        tmpat(:,ipitors) = 0
                        tmpat(1,ipitors) = a
                        do i=conn(1)%ri(a), conn(1)%ri(a+1)-1
                            c = conn(1)%ci(i)
                            if(c /= b) then
                                if(tmpat(2,ipitors) == 0) then
                                    tmpat(2,ipitors) = c
                                else
                                    tmpat(3,ipitors) = c
                                end if
                            end if
                        end do

                        tmpat(4,ipitors) = b
                        do i=conn(1)%ri(b), conn(1)%ri(b+1)-1
                            c = conn(1)%ci(i)
                            if(c /= a) then
                                if(tmpat(5,ipitors) == 0) then
                                    tmpat(5,ipitors) = c
                                else
                                    tmpat(6,ipitors) = c
                                end if
                            end if
                        end do
                        tmpk(ipitors) = kpi(iprm)
                        
                        ipitors = ipitors+1
                        exit
                    end if
                end do
            end do
        end do

        call pitors_init(ipitors-1)
        
        do i=1, ipitors-1
            kpitors(i) = tmpk(i) * kcalmol2au
            pitorsat(:,i) = tmpat(:,i)
        end do

        call mfree('assign_pitors [classa]', classa)
        call mfree('assign_pitors [classb]', classb)
        call mfree('assign_pitors [kpi]', kpi)
        call mfree('assign_pitors [tmpat]', tmpat)
        call mfree('assign_pitors [tmpk]', tmpk)
    
    end subroutine assign_pitors
    
    subroutine assign_torsion(prm_file, my_attype)
        use mod_memory, only: mallocate, mfree
        use mod_mmpol, only: fatal_error, mm_atoms, conn
        use mod_bonded, only: torsion_init, torsionat, torsamp, torsphase, torsn
        use mod_constants, only: kcalmol2au, deg2rad, eps_rp
        
        implicit none
        
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file
        integer(ip), intent(in) :: my_attype(:)
        !! List of atom types that shoukd be used to populate parameter
        !! vectors

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: ist, i, j, tokb, toke, it, nt, &
                       cla, clb, clc, cld, maxt, a, b, c, d, jb, jc, jd, iprm, ji, period
        character(len=120) :: line, errstring
        integer(ip), allocatable :: classa(:), classb(:), classc(:), classd(:), &
                                    t_n(:,:), tmpat(:,:), tmpprm(:)
        real(rp), allocatable :: t_amp(:,:), t_pha(:,:)
        real(rp) :: amp, phase, torsion_unit

        if(.not. allocated(atclass)) call read_atom_cards(prm_file)
        
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
        nt = 1
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            line = str_to_lower(line)
            if(line(:8) == 'torsion ') nt = nt + 1
        end do

        maxt = conn(4)%ri(mm_atoms+1)-1 
        call mallocate('assign_torsion [classa]', nt, classa)
        call mallocate('assign_torsion [classb]', nt, classb)
        call mallocate('assign_torsion [classc]', nt, classc)
        call mallocate('assign_torsion [classd]', nt, classd)
        call mallocate('assign_torsion [t_amp]', 6, nt, t_amp)
        call mallocate('assign_torsion [t_pha]', 6, nt, t_pha)
        call mallocate('assign_torsion [t_n]', 6, nt, t_n)
        call mallocate('assign_torsion [tmpat]', 4, maxt, tmpat)
        call mallocate('assign_torsion [tmpprm]', maxt, tmpprm)

        ! Restart the reading from the beginning to actually save the parameters
        rewind(iof_prminp)
        ist = 0
        it = 1
        i=1
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            line = str_to_lower(line)
           
            if(line(:12) == 'torsionunit ') then
                tokb = 13
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong TORSIONUNIT card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) torsion_unit

            else if(line(:8) == 'torsion ') then
                tokb = 9
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong TORSION card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classa(it)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong TORSION card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classb(it)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong TORSION card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classc(it)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong TORSION card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classd(it)
                
                ji = 1
                t_n(:,it) = -1
                do j=1, 6
                    tokb = toke + 1
                    toke = tokenize(line, tokb)
                    if(toke < 0) exit

                    if(.not. isreal(line(tokb:toke))) then
                        write(errstring, *) "Wrong TORSION card"
                        call fatal_error(errstring)
                    end if
                    read(line(tokb:toke), *) amp
                    
                    tokb = toke + 1
                    toke = tokenize(line, tokb)
                    if(.not. isreal(line(tokb:toke))) then
                        write(errstring, *) "Wrong TORSION card"
                        call fatal_error(errstring)
                    end if
                    read(line(tokb:toke), *) phase
                
                    tokb = toke + 1
                    toke = tokenize(line, tokb)
                    if(.not. isint(line(tokb:toke))) then
                        write(errstring, *) "Wrong TORSION card"
                        call fatal_error(errstring)
                    end if
                    read(line(tokb:toke), *) period

                    if(abs(amp) > eps_rp) then
                        t_amp(ji, it) = amp 
                        t_pha(ji, it) = phase
                        t_n(ji, it) = period
                        ji = ji + 1
                    end if
                end do

                if(j == 1) then
                    ! No parameter found
                    write(errstring, *) "Wrong TORSION card"
                    call fatal_error(errstring)
                end if
                
                it = it + 1
            end if
            i = i+1
        end do
        close(iof_prminp)

        it = 1
        do a=1, mm_atoms
            cla = atclass(my_attype(a))
            do jb=conn(1)%ri(a), conn(1)%ri(a+1)-1
                b = conn(1)%ci(jb)
                clb = atclass(my_attype(b))
                do jc=conn(1)%ri(b), conn(1)%ri(b+1)-1
                    c = conn(1)%ci(jc)
                    if(c == a) cycle
                    clc = atclass(my_attype(c))
                    do jd=conn(1)%ri(c), conn(1)%ri(c+1)-1
                        d = conn(1)%ci(jd)
                        if(d == a .or. d == b) cycle
                        if(a > d) cycle
                        cld = atclass(my_attype(d))
                        ! There is a dihedral A-B-C-D
                        do iprm=1, nt
                            if((classa(iprm) == cla .and. &
                                classb(iprm) == clb .and. &
                                classc(iprm) == clc .and. &
                                classd(iprm) == cld) .or. &
                               (classa(iprm) == cld .and. &
                                classb(iprm) == clc .and. &
                                classc(iprm) == clb .and. &
                                 classd(iprm) == cla)) then
                                ! The parameter is ok
                                tmpat(:,it) = [a, b, c, d]
                                tmpprm(it) = iprm
                                it = it+1
                                exit
                            end if
                        end do
                    end do
                end do
            end do
        end do

        call torsion_init(it-1)
        do i=1, it-1
           torsionat(:,i) = tmpat(:,i) 
           torsamp(:,i) = t_amp(:,tmpprm(i)) * kcalmol2au * torsion_unit
           torsphase(:,i) = t_pha(:,tmpprm(i)) * deg2rad
           torsn(:,i) = t_n(:,tmpprm(i))
        end do
        
        call mfree('assign_torsion [classa]', classa)
        call mfree('assign_torsion [classb]', classb)
        call mfree('assign_torsion [classc]', classc)
        call mfree('assign_torsion [classd]', classd)
        call mfree('assign_torsion [t_amp]', t_amp)
        call mfree('assign_torsion [t_pha]', t_pha)
        call mfree('assign_torsion [t_n]', t_n)
        call mfree('assign_torsion [tmpat]', tmpat)
        call mfree('assign_torsion [tmpprm]', tmpprm)
       
    end subroutine assign_torsion
    
    subroutine assign_angle(prm_file, my_attype)
        use mod_memory, only: mallocate, mfree
        use mod_mmpol, only: fatal_error, mm_atoms, conn
        use mod_bonded, only: OMMP_ANG_SIMPLE, &
                              OMMP_ANG_H0, &
                              OMMP_ANG_H1, &
                              OMMP_ANG_H2, &
                              OMMP_ANG_INPLANE, &
                              OMMP_ANG_INPLANE_H0, &
                              OMMP_ANG_INPLANE_H1
        use mod_bonded, only: angle_cubic, angle_quartic, &
                              angle_pentic, angle_sextic, angleat, anglety, &
                              kangle, eqangle, angle_init

        use mod_constants, only: kcalmol2au, rad2deg, deg2rad
        
        implicit none
        
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file
        integer(ip), intent(in) :: my_attype(:)
        !! List of atom types that shoukd be used to populate parameter
        !! vectors

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: ist, i, j, tokb, toke, iang, nang, &
                       cla, clb, clc, maxang, a, b, c, jc, jb, k, nhenv
        character(len=120) :: line, errstring
        integer(ip), allocatable :: classa(:), classb(:), classc(:), angtype(:)
        real(rp), allocatable :: kang(:), th0ang(:)
        logical :: done

        if(.not. allocated(atclass)) call read_atom_cards(prm_file)
        
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
        nang = 1
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            line = str_to_lower(line)
            if(line(:6) == 'angle ') nang = nang + 3 
            ! One angle keyourd could stand for 3 parameters for different H-env
            if(line(:7) == 'anglep ') nang = nang + 2
            ! One angle keyourd could stand for 2 parameters for different H-env
        end do

        maxang = (conn(2)%ri(mm_atoms+1)-1) / 2
        call mallocate('assign_angle [classa]', nang, classa)
        call mallocate('assign_angle [classb]', nang, classb)
        call mallocate('assign_angle [classc]', nang, classc)
        call mallocate('assign_angle [eqang]', nang, th0ang)
        call mallocate('assign_angle [kang]', nang, kang)
        call mallocate('assign_angle [angtype]', nang, angtype)
        call angle_init(maxang)

        ! Restart the reading from the beginning to actually save the parameters
        rewind(iof_prminp)
        ist = 0
        iang = 1
        i=1
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            line = str_to_lower(line)
           
            if(line(:12) == 'angle-cubic ') then
                tokb = 13
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong ANGLE-CUBIC card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) angle_cubic
                angle_cubic = angle_cubic * rad2deg
            
            else if(line(:14) == 'angle-quartic ') then
                tokb = 15
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong ANGLE-QUARTIC card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) angle_quartic
                angle_quartic = angle_quartic * rad2deg**2
            
            else if(line(:13) == 'angle-pentic ') then
                tokb = 13
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong ANGLE-PENTIC card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) angle_pentic
                angle_pentic = angle_pentic * rad2deg**3
            
            else if(line(:13) == 'angle-sextic ') then
                tokb = 14
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong ANGLE-SEXTIC card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) angle_sextic
                angle_sextic = angle_sextic * rad2deg**4
            
            else if(line(:6) == 'angle ') then
                tokb = 7
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong ANGLE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classa(iang)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong ANGLE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classb(iang)
                
                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong ANGLE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classc(iang)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong ANGLE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) kang(iang)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong ANGLE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) th0ang(iang)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(toke < 0) then
                    ! Only one angle parameter is specified so it is good 
                    ! for all the H-envirnoment
                    angtype(iang) = OMMP_ANG_SIMPLE
                    iang = iang + 1
                else
                    ! Three equilibrim angles are specified for three different
                    ! H-environment
                    angtype(iang) = OMMP_ANG_H0
                    iang = iang + 1
                    
                    classa(iang) = classa(iang-1)
                    classb(iang) = classb(iang-1)
                    classc(iang) = classc(iang-1)
                    kang(iang) = kang(iang-1)
                    if(.not. isreal(line(tokb:toke))) then
                        write(errstring, *) "Wrong ANGLE card"
                        call fatal_error(errstring)
                    end if
                    read(line(tokb:toke), *) th0ang(iang)
                    angtype(iang) = OMMP_ANG_H1
                    
                    iang = iang + 1
                    
                    classa(iang) = classa(iang-1)
                    classb(iang) = classb(iang-1)
                    classc(iang) = classc(iang-1)
                    kang(iang) = kang(iang-1)
                    angtype(iang) = OMMP_ANG_H2
                    tokb = toke + 1
                    toke = tokenize(line, tokb)
                    if(.not. isreal(line(tokb:toke))) then
                        write(errstring, *) "Wrong ANGLE card"
                        call fatal_error(errstring)
                    end if
                    read(line(tokb:toke), *) th0ang(iang)
                    iang = iang + 1
                end if
            
            else if(line(:7) == 'anglep ') then
                tokb = 8
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong ANGLEP card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classa(iang)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong ANGLEP card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classb(iang)
                
                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong ANGLEP card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classc(iang)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong ANGLEP card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) kang(iang)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong ANGLEP card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) th0ang(iang)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(toke < 0) then
                    ! Only one angle parameter is specified so it is good 
                    ! for all the H-envirnoment
                    angtype(iang) = OMMP_ANG_INPLANE
                    iang = iang + 1
                else
                    ! Three equilibrim angles are specified for three different
                    ! H-environment
                    angtype(iang) = OMMP_ANG_INPLANE_H0
                    iang = iang + 1
                    
                    classa(iang) = classa(iang-1)
                    classb(iang) = classb(iang-1)
                    classc(iang) = classc(iang-1)
                    kang(iang) = kang(iang-1)
                    if(.not. isreal(line(tokb:toke))) then
                        write(errstring, *) "Wrong ANGLEP card"
                        call fatal_error(errstring)
                    end if
                    read(line(tokb:toke), *) th0ang(iang)
                    angtype(iang) = OMMP_ANG_INPLANE_H1
                    
                    iang = iang + 1
                end if
            end if
            i = i+1
        end do
        close(iof_prminp)
        nang = iang
        
        iang = 1
        do a=1, mm_atoms
            cla = atclass(my_attype(a))
            do jb=conn(2)%ri(a), conn(2)%ri(a+1)-1
                b = conn(2)%ci(jb)
                if(a > b) cycle
                clb = atclass(my_attype(b))
                
                do jc=conn(1)%ri(a), conn(1)%ri(a+1)-1
                    c = conn(1)%ci(jc)
                    if(all(conn(1)%ci(conn(1)%ri(b):conn(1)%ri(b+1)-1) /= c)) cycle
                    ! There is an angle in the form A-C-B
                    clc = atclass(my_attype(c))
                    done = .false.

                    do j=1, nang
                        if((cla == classa(j) &
                            .and. clb == classc(j) &
                            .and. clc == classb(j)) .or. &
                           (clb == classa(j) &
                            .and. cla == classc(j) &
                            .and. clc == classb(j))) then
                            
                            if(angtype(j) == OMMP_ANG_SIMPLE .or. &
                               angtype(j) == OMMP_ANG_INPLANE) then
                                ! For those types no check of the H 
                                ! environment is required
                                done = .true.
                                exit
                            else
                                ! Check the H-environment
                                nhenv = 0
                                do k=conn(1)%ri(c), conn(1)%ri(c+1)-1
                                    if(atz(my_attype(conn(1)%ci(k))) == 1) &
                                        nhenv = nhenv + 1
                                end do
                                if(atz(my_attype(a)) == 1) nhenv = nhenv-1 
                                if(atz(my_attype(b)) == 1) nhenv = nhenv-1 
                                
                                if(nhenv == 0 .and. ( &
                                   angtype(j) == OMMP_ANG_H0 .or. &
                                   angtype(j) == OMMP_ANG_INPLANE_H0)) then
                                    done = .true.
                                    exit
                                else if(nhenv == 1 .and. ( &
                                        angtype(j) == OMMP_ANG_H1 .or. & 
                                        angtype(j) == OMMP_ANG_INPLANE_H1)) then
                                    done = .true.
                                    exit
                                else if(nhenv == 2 .and. (&
                                        angtype(j) == OMMP_ANG_H2)) then
                                    done = .true.
                                    exit
                                end if
                            end if
                        end if
                    end do

                    if(done) then
                        angleat(1,iang) = a
                        angleat(2,iang) = c
                        angleat(3,iang) = b
                        anglety(iang) = angtype(j)
                        kangle(iang) = kang(j) * kcalmol2au
                        eqangle(iang) = th0ang(j) * deg2rad
                        iang = iang + 1
                    else
                        write(*, *) "MB22 NO ANGLE PARAM FOUND FOR", a,b, c
                    end if
                end do
            end do
        end do

        call mfree('assign_angle [classa]', classa)
        call mfree('assign_angle [classb]', classb)
        call mfree('assign_angle [classc]', classc)
        call mfree('assign_angle [eqang]', th0ang)
        call mfree('assign_angle [kang]', kang)
        call mfree('assign_angle [angtype]', angtype)
    
    end subroutine assign_angle

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
        character(len=20) :: radrule, radsize, radtype, vdwtype, epsrule
        integer(ip), allocatable :: vdwat(:), vdwpr_a(:), vdwpr_b(:)
        real(rp), allocatable :: vdw_e_prm(:), vdw_r_prm(:), vdw_f_prm(:), &
                                 vdwpr_r(:), vdwpr_e(:)
        integer(ip) :: nvdw, ivdw, atc, nvdwpr, ivdwpr
        logical :: done

        if(.not. allocated(atclass)) call read_atom_cards(prm_file)
        
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
            line = str_to_lower(line)
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

        ! Default rules for VDW (from Tinker Manual)
        radrule = "arithmetic"
        radsize = "radius"
        radtype = "r-min"
        vdwtype = "lennard-jones"
        epsrule = "geometric"

        ! Restart the reading from the beginning to actually save the parameters
        rewind(iof_prminp)
        ist = 0
        i=1
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            line = str_to_lower(line)
           
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

            else if(line(:12) == 'epsilonrule ') then
                tokb = 12
                toke = tokenize(line, tokb)
                read(line(tokb:toke), *) epsrule
            
            else if(line(:8) == 'vdwtype ') then
                tokb = 9
                toke = tokenize(line, tokb)
                read(line(tokb:toke), *) vdwtype
            
            else if(line(:11) == 'radiusrule ') then
                tokb = 12
                toke = tokenize(line, tokb)
                read(line(tokb:toke), *) radrule

            else if(line(:11) == 'radiussize ') then
                tokb = 12
                toke = tokenize(line, tokb)
                read(line(tokb:toke), *) radsize

            else if(line(:11) == 'radiustype ') then
                tokb = 12
                toke = tokenize(line, tokb)
                read(line(tokb:toke), *) radtype

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
        
        call vdw_init(vdwtype, radrule, radsize, radtype, epsrule)
        
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
            line = str_to_lower(line)
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
            line = str_to_lower(line)
          
            if(line(:13) == 'polarization ') then
                tokb = 14
                toke = tokenize(line, tokb)
                select case(line(tokb:toke))
                    case('mutual')
                        continue
                    case('direct')
                        call fatal_error("Polarization DIRECT is not supported")
                    case default
                        call fatal_error("Wrong POLARIZATION card")
                end select

            else if(line(:15) == 'polar-12-intra ') then
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
            line = str_to_lower(line)
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
            line = str_to_lower(line)
            
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
    
    subroutine assign_tortors(prm_file, my_attype)
        use mod_memory, only: mallocate, mfree
        use mod_mmpol, only: fatal_error, mm_atoms, conn
        use mod_bonded, only: tortorat, tortorprm, tortor_newmap, tortor_init 
        
        implicit none
        
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file
        integer(ip), intent(in) :: my_attype(:)
        !! List of atom types that shoukd be used to populate parameter
        !! vectors

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: ist, i, j, tokb, toke, iprm, jd, je, e, d, cle,it,cld,&
                       cla, clb, clc, a, b, c, jc, jb, itt, ndata, ntt, ibeg, iend, maxtt
        character(len=120) :: line, errstring
        integer(ip), allocatable :: classx(:,:), map_dimension(:,:), tmpat(:,:), tmpprm(:), savedmap(:)
        real(rp), allocatable :: data_map(:), ang_map(:,:)

        if(.not. allocated(atclass)) call read_atom_cards(prm_file)
        
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
        ntt = 0
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            line = str_to_lower(line)
            if(line(:8) == 'tortors ') ntt = ntt + 1
        end do

        maxtt = conn(4)%ri(mm_atoms+1)-1 
        call mallocate('assign_tortors [classx]', 5, ntt, classx)
        call mallocate('assign_tortors [map_dimension]', 2, ntt, map_dimension)
        call mallocate('assign_tortors [savedmap]', ntt, savedmap)
        call mallocate('assign_tortors [tmpat]', 5, maxtt, tmpat)
        call mallocate('assign_tortors [tmpprm]', maxtt, tmpprm)

        ! Restart the reading from the beginning to actually save the parameters
        rewind(iof_prminp)
        ist = 0
        itt = 1
        i=1
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            line = str_to_lower(line)
           
            if(line(:8) == 'tortors ') then
                tokb = 9
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong TORTORS card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classx(1,itt)
                
                tokb = toke+1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong TORTORS card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classx(2,itt)
        
                tokb = toke+1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong TORTORS card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classx(3,itt)
        
                tokb = toke+1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong TORTORS card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classx(4,itt)
        
                tokb = toke+1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong TORTORS card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classx(5,itt)

                tokb = toke+1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong TORTORS card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) map_dimension(1,itt)
                
                tokb = toke+1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong TORTORS card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) map_dimension(2,itt)
                
                itt = itt + 1
            end if
            i = i+1
        end do

        ! Allocate data space and finally read the map
        ndata = dot_product(map_dimension(1,:), map_dimension(2,:))
        call mallocate('assign_tortors [data_map]', ndata, data_map)
        call mallocate('assign_tortors [ang_map]', 2, ndata, ang_map)
        
        rewind(iof_prminp)
        ist = 0
        itt = 1
        i=1
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            line = str_to_lower(line)
           
            if(line(:8) == 'tortors ') then
                ndata = map_dimension(1,itt)*map_dimension(2,itt)
                do j=1, ndata
                    read(iof_prminp, '(A)', iostat=ist) line
                    
                    tokb = tokenize(line)
                    toke = tokenize(line, tokb)
                    if(.not. isreal(line(tokb:toke))) then
                        write(errstring, *) "Wrong TORTORS data card"
                        call fatal_error(errstring)
                    end if
                    read(line(tokb:toke), *) ang_map(1, i)
                    
                    tokb = toke+1
                    toke = tokenize(line, tokb)
                    if(.not. isreal(line(tokb:toke))) then
                        write(errstring, *) "Wrong TORTORS data card"
                        call fatal_error(errstring)
                    end if
                    read(line(tokb:toke), *) ang_map(2, i)
                    
                    tokb = toke+1
                    toke = tokenize(line, tokb)
                    if(.not. isreal(line(tokb:toke))) then
                        write(errstring, *) "Wrong TORTORS data card"
                        call fatal_error(errstring)
                    end if
                    read(line(tokb:toke), *) data_map(i)
                    i = i + 1
                end do
                itt = itt + 1
            end if
        end do
        close(iof_prminp)
        
        it = 1
        do a=1, mm_atoms
            cla = atclass(my_attype(a))
            do jb=conn(1)%ri(a), conn(1)%ri(a+1)-1
                b = conn(1)%ci(jb)
                clb = atclass(my_attype(b))
                do jc=conn(1)%ri(b), conn(1)%ri(b+1)-1
                    c = conn(1)%ci(jc)
                    if(c == a) cycle
                    clc = atclass(my_attype(c))
                    do jd=conn(1)%ri(c), conn(1)%ri(c+1)-1
                        d = conn(1)%ci(jd)
                        if(d == a .or. d == b) cycle
                        cld = atclass(my_attype(d))
                        do je=conn(1)%ri(d), conn(1)%ri(d+1)-1
                            e = conn(1)%ci(je)
                            if(e == a .or. e == b .or. e == c) cycle
                            if(a > e) cycle
                            cle = atclass(my_attype(e))
                            ! There is a dihedral pair A-B-C-D-E
                            do iprm=1, ntt
                                if((classx(1,iprm) == cla .and. &
                                    classx(2,iprm) == clb .and. &
                                    classx(3,iprm) == clc .and. &
                                    classx(4,iprm) == cld .and. &
                                    classx(5,iprm) == cle) .or. &
                                   (classx(1,iprm) == cle .and. &
                                    classx(2,iprm) == cld .and. &
                                    classx(3,iprm) == clc .and. &
                                    classx(4,iprm) == clb .and. &
                                    classx(5,iprm) == cla)) then
                                    ! The parameter is ok
                                    tmpat(:,it) = [a, b, c, d, e]
                                    tmpprm(it) = iprm
                                    it = it+1
                                    exit
                                end if
                            end do
                        end do
                    end do
                end do
            end do
        end do
        
        call tortor_init(it-1)
        savedmap = -1
        iprm = 1
        do i=1, it-1
            if(savedmap(tmpprm(i)) < 0) then
                ! If needed, save the map in the module
                ibeg = 1 
                do j=1, tmpprm(i)-1
                    ibeg = ibeg + map_dimension(1,j)*map_dimension(2,j)
                end do
                iend = ibeg + map_dimension(1,tmpprm(i))*map_dimension(2,tmpprm(i)) - 1
                call tortor_newmap(map_dimension(1,tmpprm(i)), &
                                   map_dimension(2,tmpprm(i)), &
                                   ang_map(1,ibeg:iend), &
                                   ang_map(2,ibeg:iend), &
                                   data_map(ibeg:iend))
                savedmap(tmpprm(i)) = iprm
                iprm = iprm + 1
            end if

            tortorat(:,i) = tmpat(:,i)
            tortorprm(i) = savedmap(tmpprm(i))
        end do

        call mfree('assign_tortors [classx]', classx)
        call mfree('assign_tortors [map_dimension]', map_dimension)
        call mfree('assign_tortors [savedmap]', savedmap)
        call mfree('assign_tortors [data_map]', data_map)
        call mfree('assign_tortors [ang_map]', ang_map)
        call mfree('assign_tortors [tmpat]', tmpat)
        call mfree('assign_tortors [tmpprm]', tmpprm)
    
    end subroutine assign_tortors

end module mod_prm
