module mod_prm
    !! This module handles the reading of a parameter file in .prm format and
    !! the asignament of parameters based on atom type and connectivity.

    use mod_memory, only: ip, rp
    use mod_io, only: fatal_error, ommp_message
    use mod_topology, only: ommp_topology_type
    use mod_bonded, only: ommp_bonded_type
    use mod_electrostatics, only: ommp_electrostatics_type
    use mod_constants, only: OMMP_STR_CHAR_MAX, OMMP_VERBOSE_LOW
    use mod_utils, only: isreal, isint, tokenize, count_substr_occurence, &
                         str_to_lower

    implicit none
    private


    !!public :: assign_vdw, assign_pol, assign_mpoles, assign_bond, &
    !!          assign_angle, assign_urey, assign_strbnd, assign_opb, &
    !!          assign_pitors, assign_torsion, assign_tortors, &
    !!          assign_angtor, assign_strtor, terminate_prm
    public :: assign_pol, assign_mpoles
    public :: assign_vdw
    public :: assign_bond, assign_angle, assign_urey, assign_strbnd, assign_opb
    public :: assign_pitors, assign_torsion, assign_tortors, assign_angtor
    public :: assign_strtor
    public :: check_keyword, get_prm_ff_type

    contains 
    
#include "prm_keywords.f90"

    function get_prm_ff_type(prm_file) result(ff_type)
        !! This function is intended to check if the ff described by prm_type
        !! is AMOEBA (or amoeba-like) or AMBER or FF of another kind.
        !! A FF is considered to be AMOEBA if: it contains multipole keywords 
        !! (and no charge keywords) and polarization MUTUAL.
        !! A FF is considered do be AMBER if contains charge keyword and no
        !! multipole keyword.
        use mod_constants, only: OMMP_FF_AMBER, &
                                 OMMP_FF_AMOEBA, &
                                 OMMP_FF_UNKNOWN
        use mod_io, only: fatal_error
        
        implicit none
        
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: i, ist, nmultipole, ncharge, tokb, toke, ff_type
        character(len=OMMP_STR_CHAR_MAX) :: line, polarization
        
        ! open tinker prm file
        open(unit=iof_prminp, &
             file=prm_file(1:len(trim(prm_file))), &
             form='formatted', &
             access='sequential', &
             iostat=ist, &
             action='read')
        
        if(ist /= 0) then
           call fatal_error('Error while opening PRM input file')
        end if

        ! Read all the lines of file just to count how large vector should be 
        ! allocated 
        ist = 0
        nmultipole = 0
        ncharge = 0
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            line = str_to_lower(line)
            if(line(:7) == 'charge ') ncharge = ncharge + 1
            if(line(:10) == 'multipole ') nmultipole = nmultipole + 1
            if(line(:13) == 'polarization ') then
                tokb = 13
                toke = tokenize(line, tokb)
                read(line(tokb:toke), '(A)') polarization
            end if
        end do

        close(iof_prminp)

        if(nmultipole > 0 .and. &
           ncharge == 0 .and. &
           polarization(:6) == 'mutual') then
           ff_type = OMMP_FF_AMOEBA
        else if(nmultipole == 0 .and. &
                ncharge > 0) then
            ff_type = OMMP_FF_AMBER
        else
            ff_type = OMMP_FF_UNKNOWN
        end if
    end function

    subroutine read_atom_cards(top, prm_file)
        use mod_memory, only: mallocate, mfree
        use mod_io, only: fatal_error
        
        implicit none
        
        type(ommp_topology_type), intent(inout) :: top
        !! Topology object
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: i, ist, iat, toke, tokb, tokb1, nquote
        character(len=OMMP_STR_CHAR_MAX) :: line, errstring
        integer(ip) :: natype
        integer(ip), allocatable, dimension(:) :: typez, typeclass


        if(.not. top%attype_initialized) then
            call fatal_error("Atom type array in topology should be initialized&
                            & before performing atomclass asignament.")
        end if
        
        ! open tinker xyz file
        open(unit=iof_prminp, &
             file=prm_file(1:len(trim(prm_file))), &
             form='formatted', &
             access='sequential', &
             iostat=ist, &
             action='read')
        
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
            if(line(:5) == 'atom ') then
                tokb = 6
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong ATOM card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) iat
                natype = max(natype, iat)
            end if
        end do
        
        call mallocate('read_prm [typeclass]', natype, typeclass)
        call mallocate('read_prm [atz]', natype, typez)
        typeclass = 0
        typez = 0

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
                read(line(tokb:toke), *) typeclass(iat)

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
                read(line(tokb:toke), *) typez(iat)

                ! Only partial reading of ATOM card is needed for now.
            end if
        end do
        close(iof_prminp)

        do i = 1, top%mm_atoms
            top%atclass(i) = typeclass(top%attype(i)) 
            top%atz(i) = typez(top%attype(i)) 
        end do
        
        top%atclass_initialized = .true.
        top%atz_initialized = .true.
        
        call mfree('read_prm [typeclass]', typeclass)
        call mfree('read_prm [atz]', typez)

    end subroutine read_atom_cards

    subroutine assign_bond(bds, prm_file)
        use mod_memory, only: mallocate, mfree
        use mod_io, only: fatal_error
        use mod_bonded, only: bond_init, ommp_bonded_type
        use mod_constants, only: angstrom2au, kcalmol2au
        
        implicit none

        type(ommp_bonded_type), intent(inout) :: bds
        !! Bonded potential data structure
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: ist, i, j, l, jat, tokb, toke, ibnd, nbnd, &
                       cla, clb
        character(len=OMMP_STR_CHAR_MAX) :: line, errstring
        integer(ip), allocatable :: classa(:), classb(:)
        real(rp), allocatable :: kbnd(:), l0bnd(:)
        logical :: done
        type(ommp_topology_type), pointer :: top

        top => bds%top

        if(.not. top%atclass_initialized .or. .not. top%atz_initialized) then
            call read_atom_cards(top, prm_file)
        end if
        
        ! We assume that all pair of bonded atoms have a bonded 
        ! parameter
        call bond_init(bds, (top%conn(1)%ri(top%mm_atoms+1)-1) / 2)
        bds%kbond = 0
        bds%l0bond = 0

        l=1
        do i=1, top%mm_atoms
            do j=top%conn(1)%ri(i), top%conn(1)%ri(i+1)-1
                jat = top%conn(1)%ci(j)
                if(i < jat) then
                    bds%bondat(1,l) = i
                    bds%bondat(2,l) = jat
                    l = l+1
                end if
            end do
        end do

        ! open tinker xyz file
        open(unit=iof_prminp, &
             file=prm_file(1:len(trim(prm_file))), &
             form='formatted', &
             access='sequential', &
             iostat=ist, &
             action='read')
        
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
                read(line(tokb:toke), *) bds%bond_cubic
                ! This parameter is 1/Angstrom
                bds%bond_cubic = bds%bond_cubic / angstrom2au
            
            else if(line(:13) == 'bond-quartic ') then
                tokb = 14
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong BOND-QUARTIC card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) bds%bond_quartic
                bds%bond_quartic = bds%bond_quartic / (angstrom2au**2)
            
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
        
        do i=1, size(bds%bondat,2)
            ! Atom class for current pair
            cla = top%atclass(bds%bondat(1,i))
            clb = top%atclass(bds%bondat(2,i))
            
            done = .false.
            do j=1, nbnd
                if((classa(j)==cla .and. classb(j)==clb) .or. &
                   (classa(j)==clb .and. classb(j)==cla)) then
                    done = .true.
                    bds%kbond(i) = kbnd(j) * kcalmol2au / (angstrom2au**2)
                    bds%l0bond(i) = l0bnd(j) * angstrom2au
                    exit
                end if
            end do
            if(.not. done) then
                write(*, *) cla, clb, bds%bondat(1,i), bds%bondat(2,i)
                call fatal_error("Bond parameter not found!")
            end if
        end do
        
        call mfree('assign_bond [classa]', classa)
        call mfree('assign_bond [classb]', classb)
        call mfree('assign_bond [l0bnd]', l0bnd)
        call mfree('assign_bond [kbnd]', kbnd)
    
    end subroutine assign_bond
    
    subroutine assign_urey(bds, prm_file)
        use mod_memory, only: mallocate, mfree
        use mod_bonded, only: urey_init
        use mod_constants, only: angstrom2au, kcalmol2au
        
        implicit none
        

        type(ommp_bonded_type), intent(inout) :: bds
        !! Bonded potential data structure
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: ist, i, j, tokb, toke, iub, nub, &
                       cla, clb, clc, maxub, a, b, c, jc, jb 
        character(len=OMMP_STR_CHAR_MAX) :: line, errstring
        integer(ip), allocatable :: classa(:), classb(:), classc(:), ubtmp(:)
        real(rp), allocatable :: kub(:), l0ub(:)
        logical :: done
        type(ommp_topology_type), pointer :: top

        top => bds%top

        if(.not. top%atclass_initialized .or. .not. top%atz_initialized) then
            call read_atom_cards(top, prm_file)
        end if
        
        ! open tinker xyz file
        open(unit=iof_prminp, &
             file=prm_file(1:len(trim(prm_file))), &
             form='formatted', &
             access='sequential', &
             iostat=ist, &
             action='read')
        
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

        maxub = top%conn(2)%ri(top%mm_atoms+1)-1 
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
                read(line(tokb:toke), *) bds%urey_cubic
                ! This parameter is 1/Angstrom
                bds%urey_cubic = bds%urey_cubic / angstrom2au
            
            else if(line(:13) == 'urey-quartic ') then
                tokb = 14
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong UREY-QUARTIC card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) bds%urey_quartic
                bds%urey_quartic = bds%urey_quartic / (angstrom2au**2)
            
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
        do a=1, top%mm_atoms
            cla = top%atclass(a)
            do jb=top%conn(2)%ri(a), top%conn(2)%ri(a+1)-1
                b = top%conn(2)%ci(jb)
                done = .false.
                if(a > b) cycle
                clb = top%atclass(b)
                
                do jc=top%conn(1)%ri(a), top%conn(1)%ri(a+1)-1
                    c = top%conn(1)%ci(jc)
                    if(all(top%conn(1)%ci(top%conn(1)%ri(b): &
                                          top%conn(1)%ri(b+1)-1) /= c)) cycle
                    ! There is an angle in the form A-C-B
                    clc = top%atclass(c)
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

        call urey_init(bds, count(ubtmp > 0))
        iub = 1
        do a=1, top%mm_atoms
            do jb=top%conn(2)%ri(a), top%conn(2)%ri(a+1)-1
                if(ubtmp(jb) > 0) then
                    bds%ureyat(1,iub) = a
                    bds%ureyat(2,iub) = top%conn(2)%ci(jb)
                    bds%kurey(iub) = kub(ubtmp(jb)) * kcalmol2au / (angstrom2au**2) 
                    bds%l0urey(iub) = l0ub(ubtmp(jb)) * angstrom2au
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
    
    subroutine assign_strbnd(bds, prm_file)
        use mod_memory, only: mallocate, mfree
        use mod_bonded, only: strbnd_init 
        use mod_constants, only: kcalmol2au, angstrom2au
        
        implicit none
       
        type(ommp_bonded_type), intent(inout) :: bds
        !! Bonded potential data structure
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: ist, i, j, tokb, toke, isb, nstrbnd, &
                       cla, clb, clc, a, b, c, jc, jb, maxsb, &
                       l1a, l1b, l2a, l2b
        character(len=OMMP_STR_CHAR_MAX) :: line, errstring
        integer(ip), allocatable :: classa(:), classb(:), classc(:), sbtmp(:), &
                                    sbattmp(:, :)
        real(rp), allocatable :: k1(:), k2(:)
        logical :: done, thet_done, l1_done, l2_done
        type(ommp_topology_type), pointer :: top

        top => bds%top

        if(.not. top%atclass_initialized .or. .not. top%atz_initialized) then
            call read_atom_cards(top, prm_file)
        end if
        
        ! open tinker xyz file
        open(unit=iof_prminp, &
             file=prm_file(1:len(trim(prm_file))), &
             form='formatted', &
             access='sequential', &
             iostat=ist, &
             action='read')
        
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

        maxsb = (top%conn(2)%ri(top%mm_atoms+1)-1) / 2
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
        do a=1, top%mm_atoms
            cla = top%atclass(a)
            do jb=top%conn(2)%ri(a), top%conn(2)%ri(a+1)-1
                b = top%conn(2)%ci(jb)
                if(a > b) cycle
                clb = top%atclass(b)
                
                do jc=top%conn(1)%ri(a), top%conn(1)%ri(a+1)-1
                    c = top%conn(1)%ci(jc)
                    if(all(top%conn(1)%ci(top%conn(1)%ri(b):&
                                          top%conn(1)%ri(b+1)-1) /= c)) cycle
                    ! There is an angle in the form A-C-B
                    clc = top%atclass(c)
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

        call strbnd_init(bds, isb-1)

        do i=1, isb-1
            ! First assign the parameters
            bds%strbndat(:,i) = sbattmp(:,i)
            j = abs(sbtmp(i))
            if(sbtmp(i) > 0) then
                bds%strbndk1(i) = k1(j) * kcalmol2au / angstrom2au
                bds%strbndk2(i) = k2(j) * kcalmol2au / angstrom2au
            else
                bds%strbndk1(i) = k2(j) * kcalmol2au / angstrom2au
                bds%strbndk2(i) = k1(j) * kcalmol2au / angstrom2au
            end if
            
            l1a = min(bds%strbndat(1,i), bds%strbndat(2,i))
            l1b = max(bds%strbndat(1,i), bds%strbndat(2,i))
            l2a = min(bds%strbndat(3,i), bds%strbndat(2,i))
            l2b = max(bds%strbndat(3,i), bds%strbndat(2,i))
           
            ! Now search for the corresponding bond and angle parameters to
            ! set the equilibrium distances and angle
            l1_done = .false.
            l2_done = .false.
            thet_done = .false.

            do j=1, size(bds%bondat, 2)
                if(l1a == bds%bondat(1,j) .and. l1b == bds%bondat(2,j)) then
                    l1_done = .true.
                    bds%strbndl10(i) = bds%l0bond(j)
                else if(l2a == bds%bondat(1,j) .and. l2b == bds%bondat(2,j)) then
                    l2_done = .true.
                    bds%strbndl20(i) = bds%l0bond(j)
                end if

                if(l1_done .and. l2_done) exit
            end do
            
            do j=1, size(bds%angleat, 2)
                if(all(bds%strbndat(:,i) == bds%angleat(:,j))) then
                    thet_done = .true.
                    bds%strbndthet0(i) = bds%eqangle(j)
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
    
    subroutine assign_opb(bds, prm_file)
        use mod_memory, only: mallocate, mfree
        use mod_bonded, only: opb_init

        use mod_constants, only: kcalmol2au, rad2deg
        
        implicit none
        
        type(ommp_bonded_type), intent(inout) :: bds
        !! Bonded potential data structure
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: ist, i, tokb, toke, iopb, nopb, &
                       cla, clb, clc, cld, maxopb, a, b, c, d, jc, jb, iprm
        character(len=OMMP_STR_CHAR_MAX) :: line, errstring, opb_type
        integer(ip), allocatable :: classa(:), classb(:), classc(:), & 
                                    classd(:), tmpat(:,:)
        real(rp), allocatable :: kopbend(:), tmpk(:)
        type(ommp_topology_type), pointer :: top

        top => bds%top

        if(.not. top%atclass_initialized .or. .not. top%atz_initialized) then
            call read_atom_cards(top, prm_file)
        end if
        
        ! open tinker xyz file
        open(unit=iof_prminp, &
             file=prm_file(1:len(trim(prm_file))), &
             form='formatted', &
             access='sequential', &
             iostat=ist, &
             action='read')
        
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


        if(nopb == 0) then
           ! If there are no OPB terms, just stop here.
           return     
        end if
        maxopb = top%mm_atoms*3
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
                read(line(tokb:toke), *) bds%opb_cubic
                bds%opb_cubic = bds%opb_cubic * rad2deg
            
            else if(line(:15) == 'opbend-quartic ') then
                tokb = 16
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong OPBEND-QUARTIC card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) bds%opb_quartic
                bds%opb_quartic = bds%opb_quartic * rad2deg**2
            
            else if(line(:14) == 'opbend-pentic ') then
                tokb = 15
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong OPBEND-PENTIC card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) bds%opb_pentic
                bds%opb_pentic = bds%opb_pentic * rad2deg**3
            
            else if(line(:14) == 'opbend-sextic ') then
                tokb = 15
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong OPBEND-SEXTIC card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) bds%opb_sextic
                bds%opb_sextic = bds%opb_sextic * rad2deg**4
            
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
        do a=1, top%mm_atoms
            ! Check if the center is trigonal
            if(top%conn(1)%ri(a+1) - top%conn(1)%ri(a) /= 3) cycle
            cla = top%atclass(a)
            ! Loop over the atoms connected to the trignonal center
            do jb=top%conn(1)%ri(a), top%conn(1)%ri(a+1)-1
                b = top%conn(1)%ci(jb)
                clb = top%atclass(b)
              
                c = -1
                d = -1
                clc = 0
                cld = 0
                do jc=top%conn(1)%ri(a), top%conn(1)%ri(a+1)-1
                    if(top%conn(1)%ci(jc) == b) cycle
                    if(c < 0) then
                        c = top%conn(1)%ci(jc)
                        clc = top%atclass(c)
                    else if(d < 0) then
                        d = top%conn(1)%ci(jc)
                        cld = top%atclass(d)
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

        call opb_init(bds, iopb-1, trim(opb_type))
        
        do i=1, iopb-1
            bds%kopb(i) = tmpk(i) * kcalmol2au
            bds%opbat(:,i) = tmpat(:,i)
        end do

        call mfree('assign_opb [classa]', classa)
        call mfree('assign_opb [classb]', classb)
        call mfree('assign_opb [classc]', classc)
        call mfree('assign_opb [classd]', classd)
        call mfree('assign_opb [kopbend]', kopbend)
        call mfree('assign_opb [tmpat]', tmpat)
        call mfree('assign_opb [tmpk]', tmpk)
    
    end subroutine assign_opb
    
    subroutine assign_pitors(bds, prm_file)
        use mod_memory, only: mallocate, mfree
        use mod_bonded, only: pitors_init
        use mod_constants, only: kcalmol2au
        
        implicit none
        
        type(ommp_bonded_type), intent(inout) :: bds
        !! Bonded potential data structure
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: ist, i, tokb, toke, ipitors, npitors, &
                       cla, clb, maxpi, a, b, c, jb, iprm
        character(len=OMMP_STR_CHAR_MAX) :: line, errstring
        integer(ip), allocatable :: classa(:), classb(:), tmpat(:,:)
        real(rp), allocatable :: kpi(:), tmpk(:)
        type(ommp_topology_type), pointer :: top

        top => bds%top

        if(.not. top%atclass_initialized .or. .not. top%atz_initialized) then
            call read_atom_cards(top, prm_file)
        end if
        
        ! open tinker xyz file
        open(unit=iof_prminp, &
             file=prm_file(1:len(trim(prm_file))), &
             form='formatted', &
             access='sequential', &
             iostat=ist, &
             action='read')
        
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

        maxpi = top%mm_atoms 
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
        do a=1, top%mm_atoms
            ! Check if the center is trigonal
            if(top%conn(1)%ri(a+1) - top%conn(1)%ri(a) /= 3) cycle
            cla = top%atclass(a)
            ! Loop over the atoms connected to the trignonal center
            do jb=top%conn(1)%ri(a), top%conn(1)%ri(a+1)-1
                b = top%conn(1)%ci(jb)
                ! This avoid to compute the functions two times
                if(a > b) cycle
                
                !Check if the second center is trigonal
                if(top%conn(1)%ri(b+1) - top%conn(1)%ri(b) /= 3) cycle
                clb = top%atclass(b)

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
                        do i=top%conn(1)%ri(a), top%conn(1)%ri(a+1)-1
                            c = top%conn(1)%ci(i)
                            if(c /= b) then
                                if(tmpat(2,ipitors) == 0) then
                                    tmpat(2,ipitors) = c
                                else
                                    tmpat(3,ipitors) = c
                                end if
                            end if
                        end do

                        tmpat(4,ipitors) = b
                        do i=top%conn(1)%ri(b), top%conn(1)%ri(b+1)-1
                            c = top%conn(1)%ci(i)
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
        
        call pitors_init(bds, ipitors-1)
        
        do i=1, ipitors-1
            bds%kpitors(i) = tmpk(i) * kcalmol2au
            bds%pitorsat(:,i) = tmpat(:,i)
        end do

        call mfree('assign_pitors [classa]', classa)
        call mfree('assign_pitors [classb]', classb)
        call mfree('assign_pitors [kpi]', kpi)
        call mfree('assign_pitors [tmpat]', tmpat)
        call mfree('assign_pitors [tmpk]', tmpk)
    
    end subroutine assign_pitors
    
    subroutine assign_torsion(bds, prm_file)
        use mod_memory, only: mallocate, mfree
        use mod_bonded, only: torsion_init
        use mod_constants, only: kcalmol2au, deg2rad, eps_rp
        
        implicit none
        
        type(ommp_bonded_type), intent(inout) :: bds
        !! Bonded potential data structure
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: ist, i, j, tokb, toke, it, nt, &
                       cla, clb, clc, cld, maxt, a, b, c, d, jb, jc, jd, iprm, ji, period
        character(len=OMMP_STR_CHAR_MAX) :: line, errstring
        integer(ip), allocatable :: classa(:), classb(:), classc(:), classd(:), &
                                    t_n(:,:), tmpat(:,:), tmpprm(:)
        real(rp), allocatable :: t_amp(:,:), t_pha(:,:)
        real(rp) :: amp, phase, torsion_unit = 1.0
        type(ommp_topology_type), pointer :: top

        top => bds%top

        if(.not. top%atclass_initialized .or. .not. top%atz_initialized) then
            call read_atom_cards(top, prm_file)
        end if
        
        ! open tinker xyz file
        open(unit=iof_prminp, &
             file=prm_file(1:len(trim(prm_file))), &
             form='formatted', &
             access='sequential', &
             iostat=ist, &
             action='read')
        
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

        maxt = top%conn(4)%ri(top%mm_atoms+1)-1 
        call mallocate('assign_torsion [classa]', nt, classa)
        call mallocate('assign_torsion [classb]', nt, classb)
        call mallocate('assign_torsion [classc]', nt, classc)
        call mallocate('assign_torsion [classd]', nt, classd)
        call mallocate('assign_torsion [t_amp]', 6, nt, t_amp)
        call mallocate('assign_torsion [t_pha]', 6, nt, t_pha)
        call mallocate('assign_torsion [t_n]', 6, nt, t_n)
        call mallocate('assign_torsion [tmpat]', 4, maxt, tmpat)
        call mallocate('assign_torsion [tmpprm]', maxt, tmpprm)
        t_amp = 0.0
        t_pha = 0.0
        t_n = 1
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
        do a=1, top%mm_atoms
            cla = top%atclass(a)
            do jb=top%conn(1)%ri(a), top%conn(1)%ri(a+1)-1
                b = top%conn(1)%ci(jb)
                clb = top%atclass(b)
                do jc=top%conn(1)%ri(b), top%conn(1)%ri(b+1)-1
                    c = top%conn(1)%ci(jc)
                    if(c == a) cycle
                    clc = top%atclass(c)
                    do jd=top%conn(1)%ri(c), top%conn(1)%ri(c+1)-1
                        d = top%conn(1)%ci(jd)
                        if(d == a .or. d == b) cycle
                        if(a > d) cycle
                        cld = top%atclass(d)
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

        call torsion_init(bds, it-1)
        do i=1, it-1
           bds%torsionat(:,i) = tmpat(:,i) 
           bds%torsamp(:,i) = t_amp(:,tmpprm(i)) * kcalmol2au * torsion_unit
           bds%torsphase(:,i) = t_pha(:,tmpprm(i)) * deg2rad
           bds%torsn(:,i) = t_n(:,tmpprm(i))
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
    
    subroutine assign_strtor(bds, prm_file)
        use mod_memory, only: mallocate, mfree
        use mod_bonded, only: strtor_init
        use mod_constants, only: kcalmol2au, deg2rad, eps_rp, angstrom2au
        
        implicit none
        
        type(ommp_bonded_type), intent(inout) :: bds
        !! Bonded potential data structure
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: ist, i, j, tokb, toke, it, nt, &
                       cla, clb, clc, cld, maxt, a, b, c, d, jb, jc, jd, iprm, ji, period
        character(len=OMMP_STR_CHAR_MAX) :: line, errstring
        integer(ip), allocatable :: classa(:), classb(:), classc(:), classd(:), &
                                    tmpat(:,:), tmpprm(:)
        real(rp), allocatable :: kat(:,:)
        real(rp) :: phase, torsion_unit = 1.0
        logical :: tor_done, bnd1_done, bnd2_done, bnd3_done
        type(ommp_topology_type), pointer :: top

        top => bds%top

        if(.not. top%atclass_initialized .or. .not. top%atz_initialized) then
            call read_atom_cards(top, prm_file)
        end if
        
        
        ! open tinker xyz file
        open(unit=iof_prminp, &
             file=prm_file(1:len(trim(prm_file))), &
             form='formatted', &
             access='sequential', &
             iostat=ist, &
             action='read')
        
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
            if(line(:8) == 'strtors ') nt = nt + 1
        end do

        maxt = top%conn(4)%ri(top%mm_atoms+1)-1 
        call mallocate('assign_strtor [classa]', nt, classa)
        call mallocate('assign_strtor [classb]', nt, classb)
        call mallocate('assign_strtor [classc]', nt, classc)
        call mallocate('assign_strtor [classd]', nt, classd)
        call mallocate('assign_strtor [kat]', 9, nt, kat)
        call mallocate('assign_strtor [tmpat]', 4, maxt, tmpat)
        call mallocate('assign_strtor [tmpprm]', maxt, tmpprm)

        ! Restart the reading from the beginning to actually save the parameters
        rewind(iof_prminp)
        ist = 0
        it = 1
        i=1
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            line = str_to_lower(line)
           
            if(line(:8) == 'strtors ') then
                tokb = 9
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong STRTORS card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classa(it)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong STRTORS card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classb(it)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong STRTORS card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classc(it)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong STRTORS card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classd(it)
                
                do j=1, 9
                    tokb = toke + 1
                    toke = tokenize(line, tokb)

                    if(.not. isreal(line(tokb:toke))) then
                        write(errstring, *) "Wrong STRTORS card"
                        call fatal_error(errstring)
                    end if
                    read(line(tokb:toke), *) kat(j,it)
                end do

                it = it + 1
            end if
            i = i+1
        end do
        close(iof_prminp)

        it = 1
        do a=1, top%mm_atoms
            cla = top%atclass(a)
            do jb=top%conn(1)%ri(a), top%conn(1)%ri(a+1)-1
                b = top%conn(1)%ci(jb)
                clb = top%atclass(b)
                do jc=top%conn(1)%ri(b), top%conn(1)%ri(b+1)-1
                    c = top%conn(1)%ci(jc)
                    if(c == a) cycle
                    clc = top%atclass(c)
                    do jd=top%conn(1)%ri(c), top%conn(1)%ri(c+1)-1
                        d = top%conn(1)%ci(jd)
                        if(d == a .or. d == b) cycle
                        if(a > d) cycle
                        cld = top%atclass(d)
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

        call strtor_init(bds, it-1)
        do i=1, it-1
            bds%strtorat(:,i) = tmpat(:,i) 
            if(classa(tmpprm(i)) == top%atclass(bds%strtorat(1,i))) then
                bds%strtork(:,i) = kat(:,tmpprm(i))
            else
                bds%strtork(1:3,i) = kat(7:9,tmpprm(i))
                bds%strtork(4:6,i) = kat(4:6,tmpprm(i))
                bds%strtork(7:9,i) = kat(1:3,tmpprm(i))
            end if
            bds%strtork(:,i) = bds%strtork(:,i) * kcalmol2au / angstrom2au

            tor_done = .false.
            do j=1, size(bds%torsionat, 2)
                if(all(bds%strtorat(:,i) == bds%torsionat(:,j))) then
                    tor_done = .true.
                    bds%strtor_t(i) = j
                    exit
                end if
            end do
            
            bnd1_done = .false.
            bnd2_done = .false.
            bnd3_done = .false.
            do j=1, size(bds%bondat, 2)
                if(all(bds%strtorat(1:2,i) == bds%bondat(:,j)) .or. &
                   all(bds%strtorat(2:1:-1,i) == bds%bondat(:,j))) then
                    bnd1_done = .true.
                    bds%strtor_b(1,i) = j
                else if(all(bds%strtorat(2:3,i) == bds%bondat(:,j)) .or. &
                        all(bds%strtorat(3:2:-1,i) == bds%bondat(:,j))) then
                    bnd2_done = .true.
                    bds%strtor_b(2,i) = j
                else if(all(bds%strtorat(3:4,i) == bds%bondat(:,j)) .or. &
                        all(bds%strtorat(4:3:-1,i) == bds%bondat(:,j))) then
                    bnd3_done = .true.
                    bds%strtor_b(3,i) = j
                end if
                if(bnd1_done .and. bnd2_done .and. bnd3_done) exit
            end do

            if(.not. (tor_done .and. bnd1_done .and. bnd2_done .and. bnd3_done)) then
                call fatal_error('Ill defined stretching-torsion coupling parameter')
            end if
        end do
        
        call mfree('assign_strtor [classa]', classa)
        call mfree('assign_strtor [classb]', classb)
        call mfree('assign_strtor [classc]', classc)
        call mfree('assign_strtor [classd]', classd)
        call mfree('assign_strtor [kat]', kat)
        call mfree('assign_strtor [tmpat]', tmpat)
        call mfree('assign_strtor [tmpprm]', tmpprm)
       
    end subroutine assign_strtor

    subroutine assign_angtor(bds, prm_file)
        use mod_memory, only: mallocate, mfree
        use mod_bonded, only: angtor_init
        use mod_constants, only: kcalmol2au
        
        implicit none
        
        type(ommp_bonded_type), intent(inout) :: bds
        !! Bonded potential data structure
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: ist, i, j, tokb, toke, it, nt, &
                       cla, clb, clc, cld, maxt, a, b, c, d, jb, jc, jd, iprm
        character(len=OMMP_STR_CHAR_MAX) :: line, errstring
        integer(ip), allocatable :: classa(:), classb(:), classc(:), classd(:), &
                                    tmpat(:,:), tmpprm(:)
        real(rp), allocatable :: kat(:,:)
        logical :: tor_done, ang1_done, ang2_done
        type(ommp_topology_type), pointer :: top

        top => bds%top

        if(.not. top%atclass_initialized .or. .not. top%atz_initialized) then
            call read_atom_cards(top, prm_file)
        end if
        
        ! open tinker xyz file
        open(unit=iof_prminp, &
             file=prm_file(1:len(trim(prm_file))), &
             form='formatted', &
             access='sequential', &
             iostat=ist, & 
             action='read')
        
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
            if(line(:8) == 'angtors ') nt = nt + 1
        end do

        maxt = top%conn(4)%ri(top%mm_atoms+1)-1 
        call mallocate('assign_angtor [classa]', nt, classa)
        call mallocate('assign_angtor [classb]', nt, classb)
        call mallocate('assign_angtor [classc]', nt, classc)
        call mallocate('assign_angtor [classd]', nt, classd)
        call mallocate('assign_angtor [kat]', 6, nt, kat)
        call mallocate('assign_angtor [tmpat]', 4, maxt, tmpat)
        call mallocate('assign_angtor [tmpprm]', maxt, tmpprm)

        ! Restart the reading from the beginning to actually save the parameters
        rewind(iof_prminp)
        ist = 0
        it = 1
        i=1
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            line = str_to_lower(line)
           
            if(line(:8) == 'angtors ') then
                tokb = 9
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong ANGTORS card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classa(it)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong ANGTORS card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classb(it)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong ANGOTORS card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classc(it)

                tokb = toke + 1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong ANGTORS card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) classd(it)
                
                do j=1, 6
                    tokb = toke + 1
                    toke = tokenize(line, tokb)

                    if(.not. isreal(line(tokb:toke))) then
                        write(errstring, *) "Wrong ANGTORS card"
                        call fatal_error(errstring)
                    end if
                    read(line(tokb:toke), *) kat(j,it)
                end do

                it = it + 1
            end if
            i = i+1
        end do
        close(iof_prminp)

        it = 1
        do a=1, top%mm_atoms
            cla = top%atclass(a)
            do jb=top%conn(1)%ri(a), top%conn(1)%ri(a+1)-1
                b = top%conn(1)%ci(jb)
                clb = top%atclass(b)
                do jc=top%conn(1)%ri(b), top%conn(1)%ri(b+1)-1
                    c = top%conn(1)%ci(jc)
                    if(c == a) cycle
                    clc = top%atclass(c)
                    do jd=top%conn(1)%ri(c), top%conn(1)%ri(c+1)-1
                        d = top%conn(1)%ci(jd)
                        if(d == a .or. d == b) cycle
                        if(a > d) cycle
                        cld = top%atclass(d)
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
        
        call angtor_init(bds, it-1)
        do i=1, it-1
            bds%angtorat(:,i) = tmpat(:,i) 
            if(classa(tmpprm(i)) == top%atclass(bds%angtorat(1,i))) then
                bds%angtork(:,i) = kat(:,tmpprm(i)) * kcalmol2au
            else
                bds%angtork(1:3,i) = kat(4:6,tmpprm(i)) * kcalmol2au
                bds%angtork(4:6,i) = kat(1:3,tmpprm(i)) * kcalmol2au
            end if

            tor_done = .false.
            do j=1, size(bds%torsionat, 2)
                if(all(bds%angtorat(:,i) == bds%torsionat(:,j))) then
                    tor_done = .true.
                    bds%angtor_t(i) = j
                    exit
                end if
            end do
            
            ang1_done = .false.
            ang2_done = .false.
            do j=1, size(bds%angleat, 2)
                if(all(bds%angtorat(1:3,i) == bds%angleat(:,j)) .or. &
                   all(bds%angtorat(1:3,i) == bds%angleat(3:1:-1,j))) then
                    ang1_done = .true.
                    bds%angtor_a(1,i) = j
                else if(all(bds%angtorat(2:4,i) == bds%angleat(:,j)) .or. &
                        all(bds%angtorat(2:4,i) == bds%angleat(3:1:-1,j))) then
                    ang2_done = .true.
                    bds%angtor_a(2,i) = j
                end if
                if(ang1_done .and. ang2_done) exit
            end do

            if(.not. (tor_done .and. ang1_done .and. ang2_done)) then
                call fatal_error('Ill defined angle-torsion coupling parameter')
            end if
            
        end do
        
        call mfree('assign_angtor [classa]', classa)
        call mfree('assign_angtor [classb]', classb)
        call mfree('assign_angtor [classc]', classc)
        call mfree('assign_angtor [classd]', classd)
        call mfree('assign_angtor [kat]', kat)
        call mfree('assign_angtor [tmpat]', tmpat)
        call mfree('assign_angtor [tmpprm]', tmpprm)
       
    end subroutine assign_angtor
    
    subroutine assign_angle(bds, prm_file)
        use mod_memory, only: mallocate, mfree
        use mod_bonded, only: OMMP_ANG_SIMPLE, &
                              OMMP_ANG_H0, &
                              OMMP_ANG_H1, &
                              OMMP_ANG_H2, &
                              OMMP_ANG_INPLANE, &
                              OMMP_ANG_INPLANE_H0, &
                              OMMP_ANG_INPLANE_H1
        use mod_bonded, only: angle_init

        use mod_constants, only: kcalmol2au, rad2deg, deg2rad
        
        implicit none
        
        type(ommp_bonded_type), intent(inout) :: bds
        !! Bonded potential data structure
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: ist, i, j, tokb, toke, iang, nang, &
                       cla, clb, clc, maxang, a, b, c, jc, jb, k, nhenv
        character(len=OMMP_STR_CHAR_MAX) :: line, errstring
        integer(ip), allocatable :: classa(:), classb(:), classc(:), angtype(:)
        real(rp), allocatable :: kang(:), th0ang(:)
        logical :: done
        type(ommp_topology_type), pointer :: top

        top => bds%top

        if(.not. top%atclass_initialized .or. .not. top%atz_initialized) then
            call read_atom_cards(top, prm_file)
        end if
        
        ! open tinker xyz file
        open(unit=iof_prminp, &
             file=prm_file(1:len(trim(prm_file))), &
             form='formatted', &
             access='sequential', &
             iostat=ist, &
             action='read')
        
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

        maxang = (top%conn(2)%ri(top%mm_atoms+1)-1) / 2
        call mallocate('assign_angle [classa]', nang, classa)
        call mallocate('assign_angle [classb]', nang, classb)
        call mallocate('assign_angle [classc]', nang, classc)
        call mallocate('assign_angle [eqang]', nang, th0ang)
        call mallocate('assign_angle [kang]', nang, kang)
        call mallocate('assign_angle [angtype]', nang, angtype)
        call angle_init(bds, maxang)

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
                read(line(tokb:toke), *) bds%angle_cubic
                bds%angle_cubic = bds%angle_cubic * rad2deg
            
            else if(line(:14) == 'angle-quartic ') then
                tokb = 15
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong ANGLE-QUARTIC card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) bds%angle_quartic
                bds%angle_quartic = bds%angle_quartic * rad2deg**2
            
            else if(line(:13) == 'angle-pentic ') then
                tokb = 13
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong ANGLE-PENTIC card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) bds%angle_pentic
                bds%angle_pentic = bds%angle_pentic * rad2deg**3
            
            else if(line(:13) == 'angle-sextic ') then
                tokb = 14
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong ANGLE-SEXTIC card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) bds%angle_sextic
                bds%angle_sextic = bds%angle_sextic * rad2deg**4
            
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
        do a=1, top%mm_atoms
            cla = top%atclass(a)
            do jb=top%conn(2)%ri(a), top%conn(2)%ri(a+1)-1
                b = top%conn(2)%ci(jb)
                if(a > b) cycle
                clb = top%atclass(b)
                
                do jc=top%conn(1)%ri(a), top%conn(1)%ri(a+1)-1
                    c = top%conn(1)%ci(jc)
                    if(all(top%conn(1)%ci(top%conn(1)%ri(b):top%conn(1)%ri(b+1)-1) /= c)) cycle
                    ! There is an angle in the form A-C-B
                    clc = top%atclass(c)
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
                                do k=top%conn(1)%ri(c), top%conn(1)%ri(c+1)-1
                                    if(top%atz(top%conn(1)%ci(k)) == 1) &
                                        nhenv = nhenv + 1
                                end do
                                if(top%atz(a) == 1) nhenv = nhenv-1 
                                if(top%atz(b) == 1) nhenv = nhenv-1 
                                
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
                        bds%angleat(1,iang) = a
                        bds%angleat(2,iang) = c
                        bds%angleat(3,iang) = b
                        bds%anglety(iang) = angtype(j)
                        bds%kangle(iang) = kang(j) * kcalmol2au
                        bds%eqangle(iang) = th0ang(j) * deg2rad
                        iang = iang + 1
                    else
                        write(errstring, *) "No angle parameter found for &
                            &atoms ", a, b, c
                        call fatal_error(errstring)
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

    subroutine assign_vdw(vdw, top, prm_file)
        use mod_memory, only: mallocate, mfree
        use mod_io, only: fatal_error
        use mod_nonbonded, only: ommp_nonbonded_type, vdw_init, vdw_set_pair
        use mod_constants, only: angstrom2au, kcalmol2au
        
        implicit none
       
        type(ommp_nonbonded_type), intent(inout) :: vdw
        !! Non-bonded structure to be initialized
        type(ommp_topology_type), intent(inout) :: top
        !! Topology structure
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: ist, i, j, l, tokb, toke
        character(len=OMMP_STR_CHAR_MAX) :: line, errstring
        character(len=20) :: radrule, radsize, radtype, vdwtype, epsrule
        integer(ip), allocatable :: vdwat(:), vdwpr_a(:), vdwpr_b(:)
        real(rp), allocatable :: vdw_e_prm(:), vdw_r_prm(:), vdw_f_prm(:), &
                                 vdwpr_r(:), vdwpr_e(:)
        integer(ip) :: nvdw, ivdw, atc, nvdwpr, ivdwpr
        logical :: done

        if(.not. top%atclass_initialized .or. .not. top%atz_initialized) then
            call read_atom_cards(top, prm_file)
        end if
        
        ! open tinker xyz file
        open(unit=iof_prminp, &
             file=prm_file(1:len(trim(prm_file))), &
             form='formatted', &
             access='sequential', &
             iostat=ist, &
             action='read')
        
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
            if(line(:6) == 'vdwpr ' .or. line(:8) == 'vdwpair ') &
                nvdwpr = nvdwpr + 1
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
                read(line(tokb:toke), *) vdw%vdw_screening(1)
            
            else if(line(:13) == 'vdw-13-scale ') then
                tokb = 14
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong VDW-12-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) vdw%vdw_screening(2)
            
            else if(line(:13) == 'vdw-14-scale ') then
                tokb = 14
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong VDW-12-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) vdw%vdw_screening(3)
            
            else if(line(:13) == 'vdw-15-scale ') then
                tokb = 14
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong VDW-12-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) vdw%vdw_screening(4)

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
            else if(line(:6) == 'vdwpr ' .or. line(:8) == 'vdwpair ') then
                if(line(:6) == 'vdwpr ') then
                    tokb = 7
                else
                    tokb = 9
                end if
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
        
        call vdw_init(vdw, top, vdwtype, radrule, radsize, radtype, epsrule)
        
        do i=1, top%mm_atoms
            ! Atom class for current atom
            atc = top%atclass(i)
            
            ! VdW parameters
            done = .false.
            do j=1, nvdw
                if(vdwat(j) == atc) then
                    done = .true.
                    vdw%vdw_e(i) = vdw_e_prm(j) * kcalmol2au
                    vdw%vdw_r(i) = vdw_r_prm(j) * angstrom2au
                    vdw%vdw_f(i) = vdw_f_prm(j)
                    exit
                end if
            end do
            if(.not. done) then
                call fatal_error("VdW parameter not found!")
            end if

            ! VdW pair parameters
            do l=1, nvdwpr
                if(vdwpr_a(l) == atc) then
                    do j=i+1, top%mm_atoms
                        if(top%atclass(j) == vdwpr_b(l)) then
                            call vdw_set_pair(vdw, i, j, &
                                              vdwpr_r(l) * angstrom2au, &
                                              vdwpr_e(l) * kcalmol2au)
                        end if
                    end do
                else if(vdwpr_b(l) == atc) then
                    do j=i+1, top%mm_atoms
                        if(top%atclass(j) == vdwpr_a(l)) then
                            call vdw_set_pair(vdw, i, j, &
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
    
    subroutine assign_pol(eel, prm_file)
        use mod_memory, only: mallocate, mfree, ip, rp
        use mod_mmpol, only: set_screening_parameters
        use mod_constants, only: angstrom2au
        
        implicit none
        
        type(ommp_electrostatics_type), intent(inout), target :: eel
        !! Electrostatics data structure to be initialized
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: ist, i, j, k, l, iat, tokb, toke
        character(len=OMMP_STR_CHAR_MAX) :: line, errstring
        
        integer(ip), allocatable :: polat(:), pgspec(:,:) 
        real(rp), allocatable :: thf(:), isopol(:)
        real(rp) :: usc(4), psc(4), pisc(4), dsc(4)

        integer(ip) :: npolarize, ipolarize
        
        type(ommp_topology_type), pointer :: top

        top => eel%top

        if(.not. top%attype_initialized) then
            call fatal_error("Atom type array in topology should be initialized&
                            & before performing polarization asignament.")
        end if

        ! open tinker xyz file
        open(unit=iof_prminp, &
             file=prm_file(1:len(trim(prm_file))), &
             form='formatted', &
             access='sequential', &
             iostat=ist, &
             action='read')
        
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
       
        if(eel%amoeba) then
            call set_screening_parameters(eel, eel%mscale, psc, dsc, usc, pisc)
            eel%mmat_polgrp = 0
        else
            call set_screening_parameters(eel, eel%mscale, psc, dsc, usc)
        end if


        ! Now assign the parameters to the atoms
        do i=1, size(top%attype)
            ! Polarization
            do j=1, npolarize
                if(polat(j) == top%attype(i)) then
                    eel%pol(i) = isopol(j) * angstrom2au**3
                    !TODO Thole factors.
                    ! Assign a polgroup label to each atom
                    if(eel%mmat_polgrp(i) == 0) &
                        eel%mmat_polgrp(i) = maxval(eel%mmat_polgrp) + 1
                    
                    ! loop over the atoms connected to ith atom
                    do k=top%conn(1)%ri(i), top%conn(1)%ri(i+1)-1
                        iat = top%conn(1)%ci(k)
                        if(any(top%attype(iat) == pgspec(:,j))) then
                            ! The two atoms are in the same group
                            if(eel%mmat_polgrp(iat) == 0) then
                                eel%mmat_polgrp(iat) = eel%mmat_polgrp(i)
                            else if(eel%mmat_polgrp(iat) /= eel%mmat_polgrp(i)) then
                                ! TODO This code have never been tested, as no
                                ! suitable case have been found
                                do l=1, top%mm_atoms
                                    if(eel%mmat_polgrp(l) == 0) then
                                        continue
                                    else if(eel%mmat_polgrp(l) == eel%mmat_polgrp(iat) &
                                            .or. eel%mmat_polgrp(l) == eel%mmat_polgrp(i)) then
                                        eel%mmat_polgrp(l) = min(eel%mmat_polgrp(iat), eel%mmat_polgrp(i))
                                    else if(eel%mmat_polgrp(l) > max(eel%mmat_polgrp(iat),eel%mmat_polgrp(i))) then
                                        eel%mmat_polgrp(l) = eel%mmat_polgrp(l) - 1
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
    
    subroutine assign_mpoles(eel, prm_file)
        use mod_memory, only: mallocate, mfree
        use mod_mmpol, only: set_screening_parameters
        use mod_constants, only: AMOEBA_ROT_NONE, &
                                 AMOEBA_ROT_Z_THEN_X, &
                                 AMOEBA_ROT_BISECTOR, &
                                 AMOEBA_ROT_Z_ONLY, &
                                 AMOEBA_ROT_Z_BISECT, &
                                 AMOEBA_ROT_3_FOLD, &
                                 eps_rp
        
        implicit none
        
        type(ommp_electrostatics_type), intent(inout) :: eel
        !! The electrostatic object to be initialized
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: ist, i, j, k, iat, tokb, toke
        character(len=OMMP_STR_CHAR_MAX) :: line, errstring
        integer(ip), allocatable :: multat(:), multax(:,:), multframe(:)
        real(rp), allocatable :: cmult(:,:)
        real(rp) :: msc(4), csc(4) 
        real(rp) :: eel_au2kcalmol, eel_scale
        real(rp) :: default_eel_au2kcalmol = 332.063713
        ! Default conversion from A.U. to kcal/mol used in electrostatics of
        ! Tinker, only used to handle electric keyword
        integer(ip) :: nmult, nchg, imult, iax(3)
        logical :: ax_found(3), found13, only12
        type(ommp_topology_type), pointer :: top

        top => eel%top
        
        if(.not. top%attype_initialized) then
            call fatal_error("Atom type array in topology should be initialized&
                            & before performing multipoles asignament.")
        end if

        ! open tinker xyz file
        open(unit=iof_prminp, &
             file=prm_file(1:len(trim(prm_file))), &
             form='formatted', &
             access='sequential', &
             iostat=ist, &
             action='read')
        
        if(ist /= 0) then
           call fatal_error('Error while opening PRM input file')
        end if

        ! Read all the lines of file just to count how large vector should be 
        ! allocated 
        ist = 0
        nmult = 0
        nchg = 0
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            line = str_to_lower(line)
            if(line(:11) == 'multipole ') nmult = nmult + 1
            if(line(:7) == 'charge ') nchg = nchg + 1
        end do
        
        ! MULTIPOLE
        call mallocate('read_prm [multat]', nmult+nchg, multat)
        call mallocate('read_prm [multframe]', nmult+nchg, multframe)
        call mallocate('read_prm [multax]', 3, nmult+nchg, multax)
        call mallocate('read_prm [cmult]', 10, nmult+nchg, cmult)
        multax = AMOEBA_ROT_NONE
        imult = 1
        
        ! Default values from Tinker manual
        msc = [0.0, 0.0, 1.0, 1.0]
        csc = [0.0, 0.0, 1.0, 1.0] 
        eel_scale = 1.0

        ! Restart the reading from the beginning to actually save the parameters
        rewind(iof_prminp)
        ist = 0
        do while(ist == 0) 
            read(iof_prminp, '(A)', iostat=ist) line
            line = str_to_lower(line)
            
            if(line(:13) == 'chg-12-scale ') then
                tokb = 14
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong CHG-12-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) csc(1)
                if(csc(1) > 1.0) csc(1) = 1/csc(1)
            
            else if(line(:13) == 'chg-13-scale ') then
                tokb = 14
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong CHG-13-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) csc(2)
                if(csc(2) > 1.0) csc(2) = 1/csc(2)
            
            else if(line(:13) == 'chg-14-scale ') then
                tokb = 14
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong CHG-14-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) csc(3)
                if(csc(3) > 1.0) csc(3) = 1/csc(3)
            
            else if(line(:13) == 'chg-15-scale ') then
                tokb = 14
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong CHG-15-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) csc(4)
                if(csc(4) > 1.0) csc(4) = 1/csc(4)
            
            else if(line(:15) == 'mpole-12-scale ') then
                tokb = 16
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong MPOLE-12-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) msc(1)
                if(msc(1) > 1.0) msc(1) = 1/msc(1)
            
            else if(line(:15) == 'mpole-13-scale ') then
                tokb = 16
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong MPOLE-13-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) msc(2)
                if(msc(2) > 1.0) msc(2) = 1/msc(2)
            
            else if(line(:15) == 'mpole-14-scale ') then
                tokb = 16
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong MPOLE-14-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) msc(3)
                if(msc(3) > 1.0) msc(3) = 1/msc(3)
            
            else if(line(:15) == 'mpole-15-scale ') then
                tokb = 16
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong MPOLE-15-SCALE card"
                    call fatal_error(errstring)
                end if
                read(line(tokb:toke), *) msc(4)
                if(msc(4) > 1.0) msc(4) = 1/msc(4)
            
            else if(line(:9) == 'electric ') then
                tokb = 10
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong ELECTRIC card"
                    call fatal_error(errstring)
                end if
                ! This keyword is used to change the conversion from A.U. to
                ! kcal/mol for the electrostatic interaction.
                ! It is a bit creepy and unclear what it's correct to do here,
                ! I think that best thing is to scale the electrostatic itself
                ! by a factor (electric/default_electric) ** 0.5
                read(line(tokb:toke), *) eel_au2kcalmol
                eel_scale = (eel_au2kcalmol/default_eel_au2kcalmol) ** 0.5

            else if(line(:7) == 'charge ') then
                tokb = 8 ! len of keyword + 1
                toke = tokenize(line, tokb)
                if(.not. isint(line(tokb:toke))) then
                    write(errstring, *) "Wrong CHARGE card"
                    call fatal_error(errstring)
                endif
                read(line(tokb:toke), *) multat(imult)
                if(multat(imult) < 0) then
                    write(errstring, *) "Wrong CHARGE card (specific atom not supported)"
                    call fatal_error(errstring)
                end if
                ! Rotation axis information, charge does not contain any
                multax(:,imult) = 0
                multframe(imult) = AMOEBA_ROT_NONE

                tokb = toke+1
                toke = tokenize(line, tokb)
                if(.not. isreal(line(tokb:toke))) then
                    write(errstring, *) "Wrong CHARGE card"
                    call fatal_error(errstring)
                end if

                read(line(tokb:toke), *) cmult(1, imult)
                cmult(2:10, imult) = 0.0 ! Fixed dipole and quadrupole are not present
                imult = imult + 1

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
        
        if(nmult > 0 .and. nchg == 0) then
            call set_screening_parameters(eel, msc, eel%pscale, eel%dscale, &
                                          eel%uscale, eel%pscale_intra)
        else if(nmult == 0 .and. nchg > 0) then
            call set_screening_parameters(eel, csc, eel%pscale, eel%dscale, &
                                          eel%uscale)
        else if(nmult > 0 .and. nchg > 0) then
            write(errstring, *) "Unexpected FF with both multipoles and charges"
            call fatal_error(errstring)
        end if
        
        if(eel%amoeba) then
            eel%mol_frame = 0
            eel%ix = 0
            eel%iy = 0
            eel%iz = 0
        end if

        do i=1, size(top%attype)
            ! Multipoles
            only12 = .false. ! Only search for params based on 12 connectivity
            do j=1, max(nmult, nchg)
                found13 = .false. ! Parameter found is based on 13 connectivity
                if(multat(j) == top%attype(i)) then
                    ! For each center different multipoles are defined for 
                    ! different environment. So first check if the environment
                    ! is the correct one
                    
                    ! Assignament with only 1,2-neighbours.
                    ax_found = .false.
                    iax = 0_ip

                    if(multframe(j) == AMOEBA_ROT_NONE) then
                        ! No axis needed
                        ax_found = .true.
                    else if(multframe(j) == AMOEBA_ROT_Z_ONLY) then
                        ! Assignament with only-z
                        ax_found(2:3) = .true.
                        do k=top%conn(1)%ri(i), top%conn(1)%ri(i+1)-1
                            iat = top%conn(1)%ci(k)
                            if(top%attype(iat) == multax(1,j) &
                               .and. .not. ax_found(1)) then
                                ax_found(1) = .true.
                                iax(1) = iat
                            end if
                        end do
                    else
                        ! 2 or 3 axis needed
                        if(multax(3,j) == 0) ax_found(3) = .true.
                        
                        ! Using only 1,2 connectivity
                        do k=top%conn(1)%ri(i), top%conn(1)%ri(i+1)-1
                            iat = top%conn(1)%ci(k)
                            if(top%attype(iat) == multax(1,j) &
                               .and. .not. ax_found(1)) then
                                ax_found(1) = .true.
                                iax(1) = iat
                            else if(top%attype(iat) == multax(2,j) &
                                    .and. .not. ax_found(2)) then
                                ax_found(2) = .true.
                                iax(2) = iat
                            else if(top%attype(iat) == multax(3,j) &
                                    .and. .not. ax_found(3)) then
                                ax_found(3) = .true.
                                iax(3) = iat
                            end if
                        end do

                        ! Using also 1,3 connectivity
                        if(ax_found(1) .and. .not. ax_found(2)) then
                            do k=top%conn(1)%ri(iax(1)), top%conn(1)%ri(iax(1)+1)-1
                                iat = top%conn(1)%ci(k)
                                if(iat == i .or. iat == iax(1)) cycle
                                if(top%attype(iat) == multax(2,j) &
                                   .and. .not. ax_found(2) &
                                   .and. iat /= iax(1)) then
                                    ax_found(2) = .true.
                                    iax(2) = iat
                                else if(top%attype(iat) == multax(3,j) &
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
                        if(eel%amoeba) then
                            eel%ix(i) = iax(2)
                            eel%iy(i) = iax(3)
                            eel%iz(i) = iax(1)
                            eel%mol_frame(i) = multframe(j)
                            eel%q(:,i) = cmult(:,j) 
                        else
                            eel%q(1,i) = cmult(1,j) 
                        end if

                        
                        if(.not. found13) then
                            exit ! No further improvement is possible
                        else
                            only12 = .true.
                        end if
                    end if
                end if
            end do
        end do

        if(abs(eel_scale - 1.0) > eps_rp) then
            write(errstring, '(A, F10.6)') "Scaling charges by", eel_scale
            call ommp_message(errstring, OMMP_VERBOSE_LOW)
            eel%q = eel%q * eel_scale
        end if
        
        call mfree('read_prm [multat]', multat)
        call mfree('read_prm [multframe]', multframe)
        call mfree('read_prm [multax]', multax)
        call mfree('read_prm [cmult]', cmult)
    
    end subroutine assign_mpoles
    
    subroutine assign_tortors(bds, prm_file)
        use mod_memory, only: mallocate, mfree
        use mod_bonded, only: tortor_newmap, tortor_init
        use mod_constants, only: deg2rad, kcalmol2au
        
        implicit none
        
        type(ommp_bonded_type), intent(inout) :: bds
        !! Bonded potential data structure
        character(len=*), intent(in) :: prm_file
        !! name of the input PRM file

        integer(ip), parameter :: iof_prminp = 201
        integer(ip) :: ist, i, j, tokb, toke, iprm, jd, je, e, d, cle,it,cld,&
                       cla, clb, clc, a, b, c, jc, jb, itt, ndata, ntt, ibeg, iend, maxtt
        character(len=OMMP_STR_CHAR_MAX) :: line, errstring
        integer(ip), allocatable :: classx(:,:), map_dimension(:,:), tmpat(:,:), tmpprm(:), savedmap(:)
        real(rp), allocatable :: data_map(:), ang_map(:,:)
        type(ommp_topology_type), pointer :: top

        top => bds%top

        if(.not. top%atclass_initialized .or. .not. top%atz_initialized) then
            call read_atom_cards(top, prm_file)
        end if
        
        ! open tinker xyz file
        open(unit=iof_prminp, &
             file=prm_file(1:len(trim(prm_file))), &
             form='formatted', &
             access='sequential', &
             iostat=ist, &
             action='read')
        
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

        maxtt = top%conn(4)%ri(top%mm_atoms+1)-1 
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
        do a=1, top%mm_atoms
            cla = top%atclass(a)
            do jb=top%conn(1)%ri(a), top%conn(1)%ri(a+1)-1
                b = top%conn(1)%ci(jb)
                clb = top%atclass(b)
                do jc=top%conn(1)%ri(b), top%conn(1)%ri(b+1)-1
                    c = top%conn(1)%ci(jc)
                    if(c == a) cycle
                    clc = top%atclass(c)
                    do jd=top%conn(1)%ri(c), top%conn(1)%ri(c+1)-1
                        d = top%conn(1)%ci(jd)
                        if(d == a .or. d == b) cycle
                        cld = top%atclass(d)
                        do je=top%conn(1)%ri(d), top%conn(1)%ri(d+1)-1
                            e = top%conn(1)%ci(je)
                            if(e == a .or. e == b .or. e == c) cycle
                            if(a > e) cycle
                            cle = top%atclass(e)
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
        
        call tortor_init(bds, it-1)
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
                call tortor_newmap(bds, map_dimension(1,tmpprm(i)), &
                                   map_dimension(2,tmpprm(i)), &
                                   ang_map(1,ibeg:iend) * deg2rad, &
                                   ang_map(2,ibeg:iend) * deg2rad, &
                                   data_map(ibeg:iend) * kcalmol2au)
                savedmap(tmpprm(i)) = iprm
                iprm = iprm + 1
            end if

            bds%tortorat(:,i) = tmpat(:,i)
            bds%tortorprm(i) = savedmap(tmpprm(i))
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
