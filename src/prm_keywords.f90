function keyword_is_implemented(kw)
    implicit none 

    character(len=*), parameter :: imp_kwd(62) = ["angle               ", &
                                                 "angle-cubic         ", &
                                                 "anglep              ", &
                                                 "angle-pentic        ", &
                                                 "angle-quartic       ", &
                                                 "angle-sextic        ", &
                                                 "atom                ", &
                                                 "bond                ", &
                                                 "bond-cubic          ", &
                                                 "bond-quartic        ", &
                                                 "direct-11-scale     ", &
                                                 "direct-12-scale     ", &
                                                 "direct-13-scale     ", &
                                                 "direct-14-scale     ", &
                                                 "epsilonrule         ", &
                                                 "induce-12-scale     ", &
                                                 "induce-13-scale     ", &
                                                 "induce-14-scale     ", &
                                                 "induce-15-scale     ", &
                                                 "mpole-12-scale      ", &
                                                 "mpole-13-scale      ", &
                                                 "mpole-14-scale      ", &
                                                 "mpole-15-scale      ", &
                                                 "multipole           ", &
                                                 "mutual-11-scale     ", &
                                                 "mutual-12-scale     ", &
                                                 "mutual-13-scale     ", &
                                                 "mutual-14-scale     ", &
                                                 "opbend              ", &
                                                 "opbend-cubic        ", &
                                                 "opbend-pentic       ", &
                                                 "opbend-quartic      ", &
                                                 "opbend-sextic       ", &
                                                 "opbendtype          ", & 
                                                 "pitors              ", &
                                                 "polar-12-intra      ", &
                                                 "polar-12-scale      ", &
                                                 "polar-13-intra      ", &
                                                 "polar-13-scale      ", &
                                                 "polar-14-intra      ", &
                                                 "polar-14-scale      ", &
                                                 "polar-15-intra      ", &
                                                 "polar-15-scale      ", &
                                                 "polarization        ", & 
                                                 "polarize            ", &
                                                 "radiusrule          ", &
                                                 "radiustype          ", &
                                                 "radiussize          ", &
                                                 "strbnd              ", &
                                                 "torsion             ", &
                                                 "torsionunit         ", &
                                                 "tortors             ", &
                                                 "ureybrad            ", &
                                                 "urey-cubic          ", &
                                                 "urey-quartic        ", &
                                                 "vdw                 ", &
                                                 "vdw-12-scale        ", &
                                                 "vdw-13-scale        ", &
                                                 "vdw-14-scale        ", &
                                                 "vdw-15-scale        ", &
                                                 "vdwtype             ", &
                                                 "vdwpr               "]
    character(len=*) :: kw
    integer(ip) :: i, l 
    logical :: keyword_is_implemented

    keyword_is_implemented = .false.

    do i = 1, size(imp_kwd)
        l = len(trim(imp_kwd(i)))
        if(imp_kwd(i)(1:l) == kw) then
            keyword_is_implemented = .true.
            exit
        end if
    end do
end function


function keyword_is_ignored(kw)
    implicit none

    character(len=*), parameter :: ign_kwd(7) = [&
                                                 "biotype             ", &
                                                 "forcefield          ", &
                                                 "solute              ", &
                                                 "dielectric          ", & ! TODO
                                                 "strtors             ", & ! TODO
                                                 "angtors             ", & ! TODO
                                                 "anglef              "  & ! TODO
                                                 ]
    character(len=*) :: kw
    integer(ip) :: i, l 
    logical :: keyword_is_ignored

    keyword_is_ignored = .false.

    do i = 1, size(ign_kwd)
        l = len(trim(ign_kwd(i)))
        if(ign_kwd(i)(1:l) == kw) then
            keyword_is_ignored = .true.
            exit
        end if
    end do
end function

function keyword_is_recognized(kw)
    implicit none

    character(len=*), parameter :: rec_kwd(170) = ["a-expterm           ", &
                                              "angang              ", &
                                              "angangunit          ", &
                                              "angcflux            ", &
                                              "angle               ", &
                                              "angle3              ", &
                                              "angle4              ", &
                                              "angle5              ", &
                                              "angle-cubic         ", &
                                              "anglef              ", &
                                              "anglep              ", &
                                              "angle-pentic        ", &
                                              "angle-quartic       ", &
                                              "angle-sextic        ", &
                                              "angleunit           ", &
                                              "angtors             ", &
                                              "angtorunit          ", &
                                              "atom                ", &
                                              "b-expterm           ", &
                                              "biotype             ", &
                                              "bond                ", &
                                              "bond3               ", &
                                              "bond4               ", &
                                              "bond5               ", &
                                              "bond-cubic          ", &
                                              "bond-quartic        ", &
                                              "bondtype            ", &
                                              "bondunit            ", &
                                              "born-radius         ", &
                                              "c-expterm           ", &
                                              "charge              ", &
                                              "chargetransfer      ", &
                                              "chg-11-scale        ", &
                                              "chg-12-scale        ", &
                                              "chg-13-scale        ", &
                                              "chg-14-scale        ", &
                                              "chg-15-scale        ", &
                                              "chg-buffer          ", &
                                              "chgpen              ", &
                                              "chgtrn              ", &
                                              "deform              ", &
                                              "delta-halgren       ", &
                                              "d-equals-p          ", &
                                              "dielectric          ", &
                                              "dielectric-offset   ", &
                                              "diffuse-charge      ", &
                                              "diffuse-torsion     ", &
                                              "diffuse-vdw         ", &
                                              "dipole              ", &
                                              "dipole3             ", &
                                              "dipole4             ", &
                                              "dipole5             ", &
                                              "direct-11-scale     ", &
                                              "direct-12-scale     ", &
                                              "direct-13-scale     ", &
                                              "direct-14-scale     ", &
                                              "disp-12-scale       ", &
                                              "disp-13-scale       ", &
                                              "disp-14-scale       ", &
                                              "disp-15-scale       ", &
                                              "disp-correction     ", &
                                              "dispersion          ", &
                                              "electneg            ", &
                                              "electric            ", &
                                              "epsilonrule         ", &
                                              "forcefield          ", &
                                              "gamma-halgren       ", &
                                              "gausstype           ", &
                                              "hbond               ", &
                                              "improper            ", &
                                              "impropterm          ", &
                                              "impropunit          ", &
                                              "imptors             ", &
                                              "imptorterm          ", &
                                              "imptorunit          ", &
                                              "induce-12-scale     ", &
                                              "induce-13-scale     ", &
                                              "induce-14-scale     ", &
                                              "induce-15-scale     ", &
                                              "mmffangle           ", &
                                              "mmffarom            ", &
                                              "mmffbci             ", &
                                              "mmffbond            ", &
                                              "mmffbonder          ", &
                                              "mmffcovrad          ", &
                                              "mmffdefstbn         ", &
                                              "mmffequiv           ", &
                                              "mmffopbend          ", &
                                              "mmffpbci            ", &
                                              "mmff-pibond         ", &
                                              "mmffprop            ", &
                                              "mmffstrbnd          ", &
                                              "mmfftorsion         ", &
                                              "mmffvdw             ", &
                                              "mpole-12-scale      ", &
                                              "mpole-13-scale      ", &
                                              "mpole-14-scale      ", &
                                              "mpole-15-scale      ", &
                                              "multipole           ", &
                                              "mutual-11-scale     ", &
                                              "mutual-12-scale     ", &
                                              "mutual-13-scale     ", &
                                              "mutual-14-scale     ", &
                                              "opbend              ", &
                                              "opbend-cubic        ", &
                                              "opbend-pentic       ", &
                                              "opbend-quartic      ", &
                                              "opbend-sextic       ", &
                                              "opbendtype          ", &
                                              "opbendunit          ", &
                                              "opdist              ", &
                                              "opdist-cubic        ", &
                                              "opdist-pentic       ", &
                                              "opdist-quartic      ", &
                                              "opdist-sextic       ", &
                                              "opdistunit          ", &
                                              "penetration         ", &
                                              "piatom              ", &
                                              "pibond              ", &
                                              "pibond4             ", &
                                              "pibond5             ", &
                                              "pisystem            ", &
                                              "pitors              ", &
                                              "pitorsunit          ", &
                                              "polar-12-intra      ", &
                                              "polar-12-scale      ", &
                                              "polar-13-intra      ", &
                                              "polar-13-scale      ", &
                                              "polar-14-intra      ", &
                                              "polar-14-scale      ", &
                                              "polar-15-intra      ", &
                                              "polar-15-scale      ", &
                                              "polar-eps           ", &
                                              "polarization        ", &
                                              "polarize            ", &
                                              "radiusrule          ", &
                                              "radiussize          ", &
                                              "radiustype          ", &
                                              "rep-12-scale        ", &
                                              "rep-13-scale        ", &
                                              "rep-14-scale        ", &
                                              "rep-15-scale        ", &
                                              "repulsion           ", &
                                              "smoothing           ", &
                                              "solute              ", &
                                              "solvate             ", &
                                              "strbnd              ", &
                                              "strbndunit          ", &
                                              "strtors             ", &
                                              "strtorunit          ", &
                                              "torsion             ", &
                                              "torsion4            ", &
                                              "torsion5            ", &
                                              "torsionunit         ", &
                                              "tortors             ", &
                                              "tortorunit          ", &
                                              "ureybrad            ", &
                                              "urey-cubic          ", &
                                              "urey-quartic        ", &
                                              "ureyunit            ", &
                                              "vdw                 ", &
                                              "vdw-12-scale        ", &
                                              "vdw-13-scale        ", &
                                              "vdw14               ", &
                                              "vdw-14-scale        ", &
                                              "vdw-15-scale        ", &
                                              "vdwindex            ", &
                                              "vdwpr               ", &
                                              "vdwterm             ", &
                                              "vdwtype             "]
    character(len=*) :: kw
    integer(ip) :: i, l 
    logical :: keyword_is_recognized

    keyword_is_recognized = .false.

    do i = 1, size(rec_kwd)
        l = len(trim(rec_kwd(i)))
        if(rec_kwd(i)(1:l) == kw) then
            keyword_is_recognized = .true.
            exit
        end if
    end do
end function

function check_keyword(prm_file)
    use mod_memory, only : ip
    use mod_mmpol, only : fatal_error
    use mod_utils, only : starts_with_alpha, str_to_lower, tokenize
    
    implicit none

    character(len=*), intent(in) :: prm_file
    !! name of the input PRM file
    logical :: check_keyword
    
    integer(ip), parameter :: iof_prminp = 201
    integer(ip) :: ist, ibeg, iend
    character(len=120) :: line, kw
    
    ! open tinker xyz file
    open(unit=iof_prminp, &
         file=prm_file(1:len(trim(prm_file))), &
         form='formatted', &
         access='sequential', &
         iostat=ist)
    
    if(ist /= 0) then
       call fatal_error('Error while opening PRM input file')
    end if
    
    check_keyword = .true.

    do while(ist == 0) 
        read(iof_prminp, '(A)', iostat=ist) line
        line = str_to_lower(line)

        ! Only lines that start with a char do contain keyword
        if(starts_with_alpha(line)) then
            ibeg = 1
            iend = tokenize(line, ibeg)
            kw = line(ibeg:iend)
            if(keyword_is_recognized(kw)) then
                if(.not. keyword_is_implemented(kw)) then
                    if(keyword_is_ignored(kw)) then
                        write(*, *) "'", trim(kw), "' - keyword&
                                   & ignored"
                    else
                        write(*, *) "'", trim(kw), "' - keyword&
                                   & is not implemented and &
                                   &cannot be ignored."
                        check_keyword = .false.
                    end if
                end if
            else
                ! TODO only with some debug value
                write(*, *) "'", trim(kw), "' - keyword&
                           & is not recognized."
            end if
        end if
    end do

    close(iof_prminp)
end function
