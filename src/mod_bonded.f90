module mod_bonded
    !! Module to handle the bonded part of the FF, it closely follows the 
    !! AMOEBA functional form.
    
    use mod_memory, only: ip, rp
    use mod_topology, only: ommp_topology_type
    
    implicit none
    private
    
    ! Those constants are used as shorthand for the type of angle parameter
    ! that is used for a certain term. They consider two aspects: the functional
    ! form that could be a simple armonic constaint on the angle or something
    ! more involved (as the \testit{in-plane}) angle); the second aspects is 
    ! the hydrogen-environment that is the introduction of different force 
    ! constatns when the central atom is connected to a different number of
    ! hydrogen atoms.
    integer(ip), parameter, public :: OMMP_ANG_SIMPLE = 0
    !! Simple angle with no difference for hydrogen environments
    integer(ip), parameter, public :: OMMP_ANG_H0 = 1
    !! Simple angle with two different hydrogen environments
    integer(ip), parameter, public :: OMMP_ANG_H1 = 2
    !! Simple angle with three different hydrogen environments
    integer(ip), parameter, public :: OMMP_ANG_H2 = 3
    !! Simple angle with four different hydrogen environments
    integer(ip), parameter, public :: OMMP_ANG_INPLANE = 4
    !! In-plane angle with no difference for hydrogen environments
    integer(ip), parameter, public :: OMMP_ANG_INPLANE_H0 = 5
    !! In-plane angle with two different hydrogen environments
    integer(ip), parameter, public :: OMMP_ANG_INPLANE_H1 = 6
    !! In-plane angle with three different hydrogen environments

    type ommp_bonded_type
        type(ommp_topology_type), pointer :: top
        !! Data structure for topology
        
        ! Bond
        integer(ip) :: nbond
        !! Number of bond terms in the potential energy function
        integer(ip), allocatable :: bondat(:,:)
        !! Atoms involved in the ith bond term
        real(rp) :: bond_cubic, bond_quartic
        !! 3rd and 4th order terms coefficients, corresponding to 
        !! \(k^{(2)}\) and \(k^{(3)}\) 
        real(rp), allocatable :: kbond(:)    !! Force constants for bond terms
        real(rp), allocatable :: l0bond(:)   !! Equilibrium lengths for bonds
        logical :: use_bond = .false.
        !! Flag to enable the calculation of bond terms in potential 
        !! energy function

        ! Angle
        integer(ip) :: nangle
        !! Number of angle terms in the potential energy function
        integer(ip), allocatable :: angleat(:,:)
        !! Atoms involved in the ith angle term
        integer(ip), allocatable :: anglety(:)
        !! Type of function to be used for ith angle term
        integer(ip), allocatable :: angauxat(:)
        !! Auxiliary atom to be used in calculaton of ith angle term
        real(rp) :: angle_cubic, angle_quartic, angle_pentic, angle_sextic
        !! Coefficients for 3rd to 6th order terms corresponding to 
        !! \(k^{(3)}\) ... \(k^{(6)}\). 
        real(rp), allocatable :: kangle(:)
        !! Force constants for the ith angle term
        real(rp), allocatable :: eqangle(:)
        !! Equilibrium angle for the ith angle term
        logical :: use_angle = .false.
        !! Flag to enable the calculation of angle terms in potential energy 
        !! function

        ! Stretch-Bend
        integer(ip) :: nstrbnd
        !! Number of stretching-bending coupling terms in potential energy function
        integer(ip), allocatable :: strbndat(:,:)
        !! Atoms involved in the ith stretching-bending term
        real(rp), allocatable :: strbndk1(:), strbndk2(:)
        !! Force constants for the ith stretching-bending term (\(k_1\) and \(k_2\))
        real(rp), allocatable :: strbndthet0(:)
        !! Equilibrium angle for the ith stretching-bending term
        real(rp), allocatable :: strbndl10(:), strbndl20(:)
        !! Equilibrium distances for the ith stretching-bending term
        logical :: use_strbnd = .false.
        !! Flag to enable calculation of stretching-bending coupling terms in 
        !! potential energy function
        
        ! Angle-Torsion coupling
        integer(ip) :: nangtor
        integer(ip), allocatable :: angtorat(:,:), angtor_t(:), angtor_a(:,:)
        real(rp), allocatable :: angtork(:,:)
        logical :: use_angtor = .false.
        
        ! Bond-Torsion coupling
        integer(ip) :: nstrtor
        integer(ip), allocatable :: strtorat(:,:), strtor_t(:), strtor_b(:,:)
        real(rp), allocatable :: strtork(:,:)
        logical :: use_strtor = .false.
        
        ! Urey-Bradley
        integer(ip) :: nurey
        !! Number of Urey-Bradley terms in potential energy function
        integer(ip), allocatable :: ureyat(:,:)
        !! Atoms involved in ith Urey-Bradley term
        real(rp) :: urey_cubic, urey_quartic
        !! 3rd and 4th order constants for U-B potential (
        !! \(k^{(3)}\) and \(k^{(4)}\))
        real(rp), allocatable :: kurey(:)
        !! Force constants for U-B terms
        real(rp), allocatable :: l0urey(:)
        !! Equilibrium distance for U-B potentials
        logical :: use_urey = .false.
        !! Flag to enable calculation of U-B terms in the potential energy function

        ! Out-of-Plane Bending
        integer(ip) :: nopb
        !! Number of out-of-plane bending function in potential energy func.
        integer(ip), allocatable :: opbat(:,:)
        !! Atoms involved in ith oop bending function
        real(rp) :: opb_cubic=0.0, opb_quartic=0.0, opb_pentic=0.0, opb_sextic=0.0
        !! Coefficients for 3rd to 6th order terms corresponding to 
        !! \(k^{(3)}\) ... \(k^{(6)}\) for out-of-plane bending. 
        real(rp), allocatable :: kopb(:)
        !! Force constants for ith out-of plane bending
        logical :: use_opb = .false.
        !! Flag to enable out-of-plane bending calculation
        
        ! Pi-torsion 
        integer(ip) :: npitors
        integer(ip), allocatable :: pitorsat(:,:)
        real(rp), allocatable :: kpitors(:)
        logical :: use_pitors = .false.

        ! Torsion
        integer(ip) :: ntorsion
        integer(ip), allocatable :: torsionat(:,:), torsn(:,:)
        real(rp), allocatable :: torsamp(:,:), torsphase(:,:)
        logical :: use_torsion = .false.

        ! Torsion-torsion coupling (cmap)
        integer(ip) :: ntortor
        integer(ip), allocatable :: tortorat(:,:), tortorprm(:), ttmap_shape(:,:)
        real(rp), allocatable :: ttmap_ang1(:), ttmap_ang2(:), ttmap_v(:), &
                                 ttmap_vx(:), ttmap_vy(:), ttmap_vxy(:)
        logical :: use_tortor = .false.
    end type ommp_bonded_type

    public :: ommp_bonded_type
    public :: bond_init, bond_potential, bond_terminate
    public :: bonded_terminate
    
    contains

    subroutine bond_init(bds, n) 
        !! Initialize array used in calculation of bond stratching terms of
        !! potential energy

        use mod_memory, only: mallocate

        implicit none
        
        type(ommp_bonded_type) :: bds
        ! Bonded potential data structure
        integer(ip) :: n
        !! Number of bond stretching functions in the potential
        !! energy of the system
        
        if( n < 1 ) return
        bds%use_bond = .true.

        call mallocate('bond_init [bondat]', 2, n, bds%bondat)
        call mallocate('bond_init [kbond]', n, bds%kbond)
        call mallocate('bond_init [l0bond]', n, bds%l0bond)
        
        bds%nbond = n
        bds%bond_cubic = 0.0_rp
        bds%bond_quartic = 0.0_rp

    end subroutine bond_init

    subroutine bond_potential(bds, V)
        !! Compute the bond-stretching terms of the potential energy.  
        !! They are computed according to the general formula adopted in AMOEBA
        !! Force Field:
        !! \[U_{bond} = \sum_i k_i \Delta l_i^2 \large(1 + k^{(3)}\Delta l_i + 
        !! k^{(4)}\Delta l_i^2 \large)\]
        !! \[\Delta l_i = l_i - l^{(eq)}_i\]

        use mod_constants, only : eps_rp

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        ! Bonded potential data structure
        real(rp), intent(inout) :: V
        !! Bond potential, result will be added to V

        integer :: i
        logical :: use_cubic, use_quartic
        real(rp) :: dr(3), l, dl, dl2

        use_cubic = (abs(bds%bond_cubic) > eps_rp)
        use_quartic = (abs(bds%bond_quartic) > eps_rp)
        
        if(.not. bds%use_bond) return

        if(.not. use_cubic .and. .not. use_quartic) then
            ! This is just a regular harmonic potential
            do i=1, bds%nbond
                dr = bds%top%cmm(:,bds%bondat(1,i)) - &
                     bds%top%cmm(:,bds%bondat(2,i))
                l = sqrt(dot_product(dr, dr))
                dl = l - bds%l0bond(i)
                
                V = V + bds%kbond(i) * dl * dl
            end do
        else
            do i=1, bds%nbond
                dr = bds%top%cmm(:,bds%bondat(1,i)) - &
                     bds%top%cmm(:,bds%bondat(2,i))
                l = sqrt(dot_product(dr, dr))
                dl = l - bds%l0bond(i)
                dl2 = dl * dl

                V = V + bds%kbond(i)*dl2 * &
                    (1.0_rp + bds%bond_cubic*dl + bds%bond_quartic*dl2)
            end do
        end if
        
    end subroutine bond_potential
!!
!!    subroutine angle_init(n)
!!        !! Initialize arrays used in calculation of angle bending functions
!!
!!        use mod_memory, only: mallocate
!!
!!        implicit none
!!
!!        type(ommp_bonded_type) :: bds
!!        ! Bonded potential data structure
!!        integer(ip) :: n
!!        !! Number of angle bending functions in the potential
!!        !! energy of the system
!!
!!        if( n < 1 ) return
!!        use_angle = .true.
!!
!!        call mallocate('angle_init [angleat]', 3, n, angleat)
!!        call mallocate('angle_init [anglety]', n, anglety)
!!        call mallocate('angle_init [angauxat]', n, angauxat)
!!        call mallocate('angle_init [kangle]', n, kangle)
!!        call mallocate('angle_init [eqangle]', n, eqangle)
!!        
!!        nangle = n
!!        angauxat = 0
!!        angle_cubic = 0.0_rp
!!        angle_quartic = 0.0_rp
!!        angle_pentic = 0.0_rp
!!        angle_sextic = 0.0_rp
!!
!!    end subroutine angle_init
!!
!!    subroutine angle_potential(V)
!!        !! Compute angle-bending terms of the potential energy function.   
!!        !! Simple angle terms are computed according to the formula:
!!        !! \[U_{angle} = \sum_i k_i \Delta \theta_i^2 \large(1 +  
!!        !!  \sum_{j=1}^4 k^{(j+2)} \Delta \theta_i^j \large)\]
!!        !! \[\Delta \theta_i = \theta_i - \theta^{(eq)}_i\]    
!!        !! Out-of plane angle are more complex. First, central atom has to be
!!        !! a trigonal center, the other two atoms together with the auxliary 
!!        !! atom (that is the remaining one connected to the trigonal center) 
!!        !! define the projection plane. During the first run the auxiliary atom
!!        !! is found and saved.
!!        !! Then, the trigonal center is projected on the plane defined by the 
!!        !! other three atoms, and the angle is the one defined by the projection
!!        !! (which is the vertex, and the other two atoms -- the auxiliary is
!!        !! excluded). Then the same formula used for simple angle terms is used.
!!        use mod_mmpol, only: cmm, conn, fatal_error
!!
!!        implicit none
!!
!!        type(ommp_bonded_type) :: bds
!!        ! Bonded potential data structure
!!        real(rp), intent(inout) :: V
!!        !! Bond potential, result will be added to V
!!        
!!        integer(ip) :: i, j
!!        real(rp) :: l1, l2, dr1(3), dr2(3), thet, d_theta
!!        real(rp), dimension(3) :: v_dist, plv1, plv2, pln, a, b, c, prj_b, aux
!!
!!        if(.not. use_angle) return
!!        
!!        do i=1, nangle
!!            if(anglety(i) == OMMP_ANG_SIMPLE .or. &
!!               anglety(i) == OMMP_ANG_H0 .or. &
!!               anglety(i) == OMMP_ANG_H1 .or. &
!!               anglety(i) == OMMP_ANG_H2) then
!!                dr1 = cmm(:, angleat(1,i)) - cmm(:, angleat(2,i))
!!                dr2 = cmm(:, angleat(3,i)) - cmm(:, angleat(2,i))
!!                l1 = sqrt(dot_product(dr1, dr1))
!!                l2 = sqrt(dot_product(dr2, dr2))
!!
!!                thet = acos(dot_product(dr1, dr2)/(l1*l2))
!!                
!!                d_theta = thet-eqangle(i) 
!!                
!!                V = V + kangle(i) * d_theta**2 * (1.0 + angle_cubic*d_theta &
!!                    + angle_quartic*d_theta**2 + angle_pentic*d_theta**3 &
!!                    + angle_sextic*d_theta**4)
!!
!!            else if(anglety(i) == OMMP_ANG_INPLANE .or. &
!!                    anglety(i) == OMMP_ANG_INPLANE_H0 .or. &
!!                    anglety(i) == OMMP_ANG_INPLANE_H1) then
!!                if(angauxat(i) < 1) then
!!                    ! Find the auxiliary atom used to define the projection
!!                    ! plane
!!                    if(conn(1)%ri(angleat(2,i)+1) - conn(1)%ri(angleat(2,i)) /= 3) then
!!                        call fatal_error("Angle IN-PLANE defined for a non-&
!!                                         &trigonal center")
!!                    end if 
!!                    do j=conn(1)%ri(angleat(2,i)), conn(1)%ri(angleat(2,i)+1)-1
!!                        if(conn(1)%ci(j) /= angleat(1,i) .and. &
!!                           conn(1)%ci(j) /= angleat(3,i)) then
!!                            angauxat(i) = conn(1)%ci(j)
!!                        endif
!!                    end do
!!                end if
!!                a = cmm(:, angleat(1,i))
!!                b = cmm(:, angleat(2,i)) !! Trigonal center
!!                c = cmm(:, angleat(3,i))
!!
!!                aux = cmm(:, angauxat(i))
!!                plv1 = a - aux
!!                plv2 = c - aux
!!                pln(1) = plv1(2)*plv2(3) - plv1(3)*plv2(2)
!!                pln(2) = plv1(3)*plv2(1) - plv1(1)*plv2(3)
!!                pln(3) = plv1(1)*plv2(2) - plv1(2)*plv2(1)
!!                !! Normal vector of the projection plane
!!                pln = pln / sqrt(dot_product(pln, pln))
!!
!!                v_dist = b - aux
!!                prj_b = b - dot_product(v_dist, pln) * pln 
!!
!!                dr1 = cmm(:, angleat(1,i)) - prj_b
!!                dr2 = cmm(:, angleat(3,i)) - prj_b
!!                l1 = sqrt(dot_product(dr1, dr1))
!!                l2 = sqrt(dot_product(dr2, dr2))
!!
!!                thet = acos(dot_product(dr1, dr2)/(l1*l2))
!!                
!!                d_theta = thet-eqangle(i) 
!!                
!!                V = V + kangle(i) * d_theta**2 * (1.0 + angle_cubic*d_theta &
!!                    + angle_quartic*d_theta**2 + angle_pentic*d_theta**3 &
!!                    + angle_sextic*d_theta**4)
!!            end if
!!        end do
!!    end subroutine angle_potential
!!    
!!    subroutine strbnd_init(n)
!!        !! Initialize arrays for calculation of stretch-bend cross term 
!!        !! potential
!!
!!        use mod_memory, only: mallocate
!!
!!        implicit none
!!
!!        type(ommp_bonded_type) :: bds
!!        ! Bonded potential data structure
!!        integer(ip) :: n
!!        !! Number of stretch-bend functions in the potential
!!        !! energy of the system
!!
!!        if( n < 1 ) return
!!        use_strbnd = .true.
!!
!!        call mallocate('strbnd_init [strbndat]', 3, n, strbndat)
!!        call mallocate('strbnd_init [strbndl10]', n, strbndl10)
!!        call mallocate('strbnd_init [strbndl20]', n, strbndl20)
!!        call mallocate('strbnd_init [strbndthet0]', n, strbndthet0)
!!        call mallocate('strbnd_init [strbndk1]', n, strbndk1)
!!        call mallocate('strbnd_init [strbndk2]', n, strbndk2)
!!        nstrbnd = n
!!
!!    end subroutine strbnd_init
!!
!!    subroutine strbnd_potential(V)
!!        !! Compute the stretch-bend cross term potential.   
!!        !! Those terms are computed according the following formula:
!!        !! \[U_{bond/angle} = (k_i \Delta l_i + k_j \Delta l_j) 
!!        !! \Delta \theta_{ij} \]
!!        !! where \(\theta_{ij}\) is the angle delimited by the bond \(i\) and 
!!        !! \(j\).   
!!        !! The force constants \(k_i\) and \(k_j\) are explicitely defined in
!!        !! the FF, while the equilibrium values are the same as for stretching
!!        !! and bending terms.
!!
!!        use mod_mmpol, only: cmm, fatal_error
!!
!!        implicit none
!!
!!        type(ommp_bonded_type) :: bds
!!        ! Bonded potential data structure
!!        real(rp), intent(inout) :: V
!!        !! Stretch-bend cross term potential, result will be added to V
!!
!!        integer(ip) :: i
!!        real(rp) :: d_l1, d_l2, d_thet, dr1(3), dr2(3), l1, l2, thet
!!        
!!        if(.not. use_strbnd) return
!!
!!        do i=1, nstrbnd
!!            dr1 = cmm(:, strbndat(2,i)) - cmm(:, strbndat(1,i))
!!            l1 = norm2(dr1)
!!            d_l1 = l1 - strbndl10(i)
!!            
!!            dr2 = cmm(:,strbndat(2,i)) - cmm(:, strbndat(3,i))
!!            l2 = norm2(dr2)
!!            d_l2 = l2 - strbndl20(i)
!!
!!            thet = acos(dot_product(dr1, dr2)/(l1*l2))
!!            d_thet = thet - strbndthet0(i) 
!!            
!!            V = V + (d_l1*strbndk1(i) + d_l2*strbndk2(i)) * d_thet
!!        end do
!!    end subroutine strbnd_potential
!!    
!!    subroutine urey_init(n) 
!!        !! Initialize Urey-Bradley potential arrays
!!
!!        use mod_memory, only: mallocate
!!
!!        implicit none
!!
!!        type(ommp_bonded_type) :: bds
!!        ! Bonded potential data structure
!!        integer(ip) :: n
!!        !! Number of Urey-Bradley functions in the potential
!!        !! energy of the system
!!        
!!        if( n < 1 ) return
!!        use_urey = .true.
!!
!!        call mallocate('urey_init [ureya]', 2, n, ureyat)
!!        call mallocate('urey_init [kurey]', n, kurey)
!!        call mallocate('urey_init [l0urey]', n, l0urey)
!!        nurey = n
!!        urey_cubic = 0.0_rp
!!        urey_quartic = 0.0_rp
!!
!!    end subroutine urey_init
!!
!!    subroutine urey_potential(V)
!!        !! Compute the Urey-Bradley potential.  
!!        !! This is basically a virtual bond, with its stretching harminic 
!!        !! potential that connect two otherwise un-connected bonds. The same
!!        !! potential formula used for normal stretching is used.
!!
!!        use mod_constants, only : eps_rp
!!        use mod_mmpol, only: cmm
!!
!!        implicit none
!!
!!        type(ommp_bonded_type) :: bds
!!        ! Bonded potential data structure
!!        real(rp), intent(inout) :: V
!!        !! Urey-Bradley potential, result will be added to V
!!
!!        integer :: i
!!        logical :: use_cubic, use_quartic
!!        real(rp) :: dr(3), l, dl, dl2
!!        
!!        if(.not. use_urey) return
!!
!!        use_cubic = (abs(urey_cubic) > eps_rp)
!!        use_quartic = (abs(urey_quartic) > eps_rp)
!!
!!        if(.not. use_cubic .and. .not. use_quartic) then
!!            ! This is just a regular harmonic potential
!!            do i=1, nurey
!!                dr = cmm(:,ureyat(1,i)) - cmm(:,ureyat(2,i))
!!                l = sqrt(dot_product(dr, dr))
!!                dl = l - l0urey(i)
!!                V = V + kurey(i) * dl * dl
!!            end do
!!        else
!!            do i=1, nurey
!!                dr = cmm(:,ureyat(1,i)) - cmm(:,ureyat(2,i))
!!                l = sqrt(dot_product(dr, dr))
!!                dl = l - l0urey(i)
!!                dl2 = dl * dl
!!
!!                V = V + kurey(i)*dl2 * (1.0_rp + urey_cubic*dl + urey_quartic*dl2)
!!            end do
!!        end if
!!    end subroutine urey_potential
!!
!!    subroutine opb_init(n, opbtype)
!!        !! Initialize arrays for out-of-plane bending potential calculation.   
!!        !! @todo Currently only Allinger functional form is supported 
!!        use mod_io, only: ommp_message
!!        use mod_constants, only: OMMP_VERBOSE_LOW
!!        use mod_mmpol, only: fatal_error
!!        use mod_memory, only: mallocate
!!
!!        implicit none
!!
!!        type(ommp_bonded_type) :: bds
!!        ! Bonded potential data structure
!!        integer(ip) :: n
!!        !! Number of out of plane Bending functions in the potential
!!        !! energy of the system
!!        character(len=*) :: opbtype
!!
!!        select case(opbtype)
!!            case('allinger')
!!                continue
!!            case('w-d-c')
!!                call fatal_error('Out-of-plane bend W-D-C is not implemented')
!!            case default
!!                call ommp_message("Found OPB type: '"//opbtype//"'", OMMP_VERBOSE_LOW)
!!                call fatal_error('Out-of-plane type specified is not understood')
!!        end select
!!
!!        if( n < 1 ) return
!!        use_opb = .true.
!!
!!        call mallocate('opb_init [opbat]', 4, n, opbat)
!!        call mallocate('opb_init [kopb]', n, kopb)
!!        nopb = n
!!
!!    end subroutine opb_init
!!
!!    subroutine opb_potential(V)
!!        !! Computes the out-of-plane bending potential.  
!!        !! With Allinger formula: similarly to in plane angles, here we are 
!!        !! considering a trigonal center, where D is the central atom and 
!!        !! A, B, C are connected to D. Allinger formula consider the angle 
!!        !! between vector \(\vec{AD}\) and the normal vector of plane ABC, 
!!        !! using \(\frac{\pi}{2}\) as implicit equilibrium value. The formula
!!        !! for this potential term is:
!!        !! \[U_{out-of-plane} = \sum_i k_i \chi_i^2 \large(1 + 
!!        !! \sum_{j=1}^4 k^{(j+2)} \chi_i^j \large) \]
!!
!!        use mod_constants, only : pi
!!        use mod_mmpol, only: cmm
!!
!!        implicit none
!!
!!        type(ommp_bonded_type) :: bds
!!        ! Bonded potential data structure
!!        real(rp), intent(inout) :: V
!!        !! out-of-plane potential, result will be added to V
!!        real(rp), dimension(3) :: a, b, c, d, plv1, plv2, pln, vad
!!        real(rp) :: lpln, lvad, thet, thet2, thet3, thet4
!!        integer(ip) :: i
!!
!!        if(.not. use_opb) return
!!        
!!        do i=1, nopb
!!            ! A* -- D -- C
!!            !       |
!!            !       B 
!!            a = cmm(:,opbat(2,i))
!!            d = cmm(:,opbat(1,i))
!!            c = cmm(:,opbat(3,i))
!!            b = cmm(:,opbat(4,i))
!!
!!            ! Compute the normal vector of the plane
!!            plv1 = a - b
!!            plv2 = a - c
!!            pln(1) = plv1(2)*plv2(3) - plv1(3)*plv2(2)
!!            pln(2) = plv1(3)*plv2(1) - plv1(1)*plv2(3)
!!            pln(3) = plv1(1)*plv2(2) - plv1(2)*plv2(1)
!!            lpln = norm2(pln)
!!
!!            ! Vector from A to D
!!            vad = a - d
!!            lvad = norm2(vad)
!!
!!            thet = abs(pi/2.0 - acos(dot_product(vad, pln)/(lvad*lpln)))
!!            thet2 = thet*thet
!!            thet3 = thet2*thet
!!            thet4 = thet3*thet
!!            V = V +  kopb(i) * thet2 * (1 + opb_cubic*thet + opb_quartic*thet2 &
!!                + opb_pentic*thet3 + opb_sextic*thet4)
!!        end do
!!    end subroutine opb_potential
!!    
!!    subroutine pitors_init(n)
!!        !! Initialize arrays needed to compute pi-torsion potential
!!
!!        use mod_memory, only: mallocate
!!
!!        implicit none
!!
!!        type(ommp_bonded_type) :: bds
!!        ! Bonded potential data structure
!!        integer(ip) :: n
!!        !! Number of out of plane pi-torsion functions in the potential
!!        !! enerpgy of the system
!!        
!!        if( n < 1 ) return
!!        use_pitors = .true.
!!
!!        call mallocate('pitors_init [pitorsat]', 6, n, pitorsat)
!!        call mallocate('pitors_init [kpitors]', n, kpitors)
!!        npitors = n
!!
!!    end subroutine pitors_init
!!
!!    subroutine pitors_potential(V)
!!        !! Compute pi-torsion terms of the potential.  
!!        !! This potential is defined on a \(\pi\)-system, and uses the 
!!        !! coordinates of six atoms A...F the central "double" bond is A-B, then
!!        !! C and D are connected to A while E and F are connected to B. So two
!!        !! plane ACD and BEF are defined. The potential is computed using the 
!!        !! dihedral angle of the normal vector of those two planes, connected 
!!        !! by segment A-B (\(\theta\)).  
!!        !! The formula used is:
!!        !! \[U_{\pi-torsion} = \sum_i k_i \large(1 + cos(2\theta-\pi) \large)\]
!!
!!        use mod_constants, only : pi
!!        use mod_mmpol, only: cmm
!!
!!        implicit none
!!
!!        type(ommp_bonded_type) :: bds
!!        ! Bonded potential data structure
!!        real(rp), intent(inout) :: V
!!        !! pi-torsion potential, result will be added to V
!!        real(rp), dimension(3) :: a, b, c, d, e, f, u, t, cd, plv1, plv2, pln1, pln2
!!        real(rp) :: thet, costhet
!!        integer(ip) :: i
!!
!!        if(.not. use_pitors) return
!!        
!!        do i=1, npitors
!!            !
!!            !  2(c)        5(e)         a => 1
!!            !   \         /             b => 4
!!            !    1(a) -- 4(b)  
!!            !   /         \
!!            !  3(d)        6(f)
!!            
!!            ! Atoms that defines the two planes
!!            a = cmm(:,pitorsat(1,i))
!!            c = cmm(:,pitorsat(2,i))
!!            d = cmm(:,pitorsat(3,i))
!!
!!            b = cmm(:,pitorsat(4,i))
!!            e = cmm(:,pitorsat(5,i))
!!            f = cmm(:,pitorsat(6,i))
!!
!!
!!            ! Compute the normal vector of the first plane
!!            plv1 = d - b
!!            plv2 = c - b
!!            pln1(1) = plv1(2)*plv2(3) - plv1(3)*plv2(2)
!!            pln1(2) = plv1(3)*plv2(1) - plv1(1)*plv2(3)
!!            pln1(3) = plv1(1)*plv2(2) - plv1(2)*plv2(1)
!!
!!            ! Compute the normal vector of the second plane
!!            plv1 = f - a
!!            plv2 = e - a
!!            pln2(1) = plv1(2)*plv2(3) - plv1(3)*plv2(2)
!!            pln2(2) = plv1(3)*plv2(1) - plv1(1)*plv2(3)
!!            pln2(3) = plv1(1)*plv2(2) - plv1(2)*plv2(1)
!!
!!            cd = b - a
!!
!!            t(1) = pln1(2)*cd(3) - pln1(3)*cd(2)
!!            t(2) = pln1(3)*cd(1) - pln1(1)*cd(3)
!!            t(3) = pln1(1)*cd(2) - pln1(2)*cd(1)
!!            t = t / norm2(t)
!!            
!!            u(1) = cd(2)*pln2(3) - cd(3)*pln2(2)
!!            u(2) = cd(3)*pln2(1) - cd(1)*pln2(3)
!!            u(3) = cd(1)*pln2(2) - cd(2)*pln2(1)
!!            u = u / norm2(u)
!!            
!!            costhet = dot_product(u,t)
!!            
!!            thet = acos(costhet)
!!
!!            V = V +  kpitors(i) * (1 + cos(2.0*thet-pi))
!!        end do
!!
!!    end subroutine pitors_potential
!!    
!!    subroutine torsion_init(n)
!!        !! Initialize torsion potential arrays
!!
!!        use mod_memory, only: mallocate
!!
!!        implicit none
!!
!!        type(ommp_bonded_type) :: bds
!!        ! Bonded potential data structure
!!        integer(ip) :: n
!!        !! Number of torsion functions in the potential
!!        !! energy of the system
!!        
!!        if( n < 1 ) return
!!        use_torsion = .true.
!!
!!        call mallocate('torsion_init [torsionat]', 4, n, torsionat)
!!        call mallocate('torsion_init [torsamp]', 6, n, torsamp)
!!        call mallocate('torsion_init [torsphase]', 6, n, torsphase)
!!        call mallocate('torsion_init [torsn]', 6, n, torsn)
!!
!!        ntorsion = n
!!
!!    end subroutine torsion_init
!!
!!    subroutine torsion_potential(V)
!!        !! Compute torsion potential
!!        use mod_constants, only: pi, eps_rp
!!
!!        implicit none
!!
!!        type(ommp_bonded_type) :: bds
!!        ! Bonded potential data structure
!!        real(rp), intent(inout) :: V
!!        !! torsion potential, result will be added to V
!!        real(rp) :: thet, costhet
!!        integer(ip) :: i, j
!!        
!!        if(.not. use_torsion) return
!!
!!        do i=1, ntorsion
!!            ! Atoms that defines the dihedral angle
!!            costhet = cos_torsion(torsionat(:,i))
!!            
!!            if(costhet + 1.0 <= eps_rp) then
!!                thet = pi
!!            else
!!                thet = acos(costhet)
!!            end if
!!
!!            do j=1, 6
!!                if(torsn(j,i) < 1) exit
!!                V = V + torsamp(j,i) * (1+cos(real(torsn(j,i))*thet - torsphase(j,i)))
!!            end do
!!        end do
!!
!!    end subroutine torsion_potential
!!    
!!    subroutine angtor_init(n)
!!        !! Initialize angle-torsion coupling potential arrays
!!
!!        use mod_memory, only: mallocate
!!
!!        implicit none
!!
!!        type(ommp_bonded_type) :: bds
!!        ! Bonded potential data structure
!!        integer(ip) :: n
!!        !! Number of angle torsion coupling functions in the potential
!!        !! energy of the system
!!        
!!        if( n < 1 ) return
!!        use_angtor = .true.
!!
!!        call mallocate('angtor_init [angtorat]', 4, n, angtorat)
!!        call mallocate('angtor_init [angtork]', 6, n, angtork)
!!        call mallocate('angtor_init [angtor_t]', n, angtor_t)
!!        call mallocate('angtor_init [angtor_a]', 2, n, angtor_a)
!!
!!        nangtor = n
!!
!!    end subroutine angtor_init
!!    
!!    subroutine strtor_init(n)
!!        
!!        use mod_memory, only: mallocate
!!
!!        implicit none
!!
!!        type(ommp_bonded_type) :: bds
!!        ! Bonded potential data structure
!!        integer(ip) :: n
!!        
!!        if( n < 1 ) return
!!        use_strtor = .true.
!!
!!        call mallocate('strtor_init [strtorat]', 4, n, strtorat)
!!        call mallocate('strtor_init [strtork]', 9, n, strtork)
!!        call mallocate('strtor_init [strtor_t]', n, strtor_t)
!!        call mallocate('strtor_init [strtor_a]', 3, n, strtor_b)
!!
!!        nstrtor = n
!!
!!    end subroutine strtor_init
!!    
!!    subroutine angtor_potential(V)
!!        use mod_mmpol, only: cmm
!!
!!        implicit none
!!
!!        type(ommp_bonded_type) :: bds
!!        ! Bonded potential data structure
!!        real(rp), intent(inout) :: V
!!        real(rp) :: thet, costhet, dihef(3), delta_a(2), vat, l1, l2, &
!!                    dr1(3), dr2(3), angle1, angle2
!!        integer(ip) :: i, j, k, ia1, ia2
!!        
!!        if(.not. use_torsion) return
!!
!!        do i=1, nangtor
!!            ! Atoms that defines the dihedral angle
!!            costhet = cos_torsion(angtorat(:,i))
!!            thet = acos(costhet)
!!            do j=1, 3
!!                dihef(j) = 1.0 + cos(j*thet+torsphase(j,angtor_t(i)))
!!            end do
!!
!!            ia1 = angtor_a(1,i)
!!            ia2 = angtor_a(2,i)
!!            
!!            dr1 = cmm(:, angleat(1,ia1)) - cmm(:, angleat(2,ia1))
!!            dr2 = cmm(:, angleat(3,ia1)) - cmm(:, angleat(2,ia1))
!!            l1 = norm2(dr1)
!!            l2 = norm2(dr2)
!!            angle1 = acos(dot_product(dr1, dr2)/(l1*l2))
!!
!!            dr1 = cmm(:, angleat(1,ia2)) - cmm(:, angleat(2,ia2))
!!            dr2 = cmm(:, angleat(3,ia2)) - cmm(:, angleat(2,ia2))
!!            l1 = norm2(dr1)
!!            l2 = norm2(dr2)
!!            angle2 = acos(dot_product(dr1, dr2)/(l1*l2))
!!           
!!            delta_a(1) = angle1 - eqangle(angtor_a(1,i))
!!            delta_a(2) = angle2 - eqangle(angtor_a(2,i))
!!
!!            do j=1,2
!!                vat = 0.0
!!                do k=1, 3
!!                    vat = vat + angtork((j-1)*3+k,i) * dihef(k)
!!                end do
!!                V = V + vat * delta_a(j)
!!            end do
!!        end do
!!
!!    end subroutine angtor_potential
!!    
!!    subroutine strtor_potential(V)
!!        use mod_mmpol, only: cmm
!!        use mod_constants
!!
!!        implicit none
!!
!!        type(ommp_bonded_type) :: bds
!!        ! Bonded potential data structure
!!        real(rp), intent(inout) :: V
!!        real(rp) :: thet, costhet, dihef(3), dr(3), r(3), vst
!!        integer(ip) :: i, j, k, ib1, ib2, ib3
!!        
!!        if(.not. use_torsion) return
!!
!!        do i=1, nstrtor
!!            ! Atoms that defines the dihedral angle
!!            costhet = cos_torsion(strtorat(:,i))
!!            thet = acos(costhet)
!!            do j=1, 3
!!                dihef(j) = 1.0 + cos(j*thet+torsphase(j,strtor_t(i)))
!!            end do
!!
!!            ib1 = strtor_b(1,i) 
!!            ib2 = strtor_b(2,i)
!!            ib3 = strtor_b(3,i)
!!            r(1) = norm2(cmm(:, bondat(1,ib1)) - cmm(:, bondat(2,ib1)))
!!            r(2) = norm2(cmm(:, bondat(1,ib2)) - cmm(:, bondat(2,ib2)))
!!            r(3) = norm2(cmm(:, bondat(1,ib3)) - cmm(:, bondat(2,ib3)))
!!            dr(1) = r(1) - l0bond(ib1)  
!!            dr(2) = r(2) - l0bond(ib2)  
!!            dr(3) = r(3) - l0bond(ib3)  
!!            
!!            do j=1,3
!!                vst = 0.0
!!                do k=1, 3
!!                    vst = vst + strtork((j-1)*3+k,i) * dihef(k)
!!                end do
!!                V = V + vst * dr(j)
!!            end do
!!        end do
!!
!!    end subroutine strtor_potential
!!    
!!    subroutine tortor_init(n)
!!        !! Initialize torsion-torsion correction potential arrays
!!
!!        use mod_memory, only: mallocate
!!
!!        implicit none
!!
!!        type(ommp_bonded_type) :: bds
!!        ! Bonded potential data structure
!!        integer(ip) :: n
!!        !! Number of torsion-torsion 'map' functions in the potential
!!        !! energy of the system
!!
!!        if( n < 1 ) return
!!        use_tortor = .true.
!!        
!!        call mallocate('torsion_init [tortorprm]', n, tortorprm )
!!        call mallocate('torsion_init [tortorat]', 5, n, tortorat)
!!
!!        ntortor = n
!!
!!    end subroutine tortor_init
!!
!!    subroutine tortor_newmap(d1, d2, ang1, ang2, v)
!!        !! Store in module memory the data describing a new torsion-torsion 
!!        !! map
!!        use mod_memory, only: mallocate, mfree
!!        use mod_utils, only: cyclic_spline
!!
!!        implicit none
!!
!!        type(ommp_bonded_type) :: bds
!!        ! Bonded potential data structure
!!        integer(ip), intent(in) :: d1, d2
!!        !! Dimensions of the new map to be saved
!!        real(rp), intent(in) :: ang1(:)
!!        !! Value of torsion1 for the new map 
!!        real(rp), intent(in) :: ang2(:)
!!        !! Value of torsion2 for the new map 
!!        real(rp), intent(in) :: v(:)
!!        !! Value of potential for the new map 
!!
!!        integer :: i, j, ii
!!        real(rp), allocatable, dimension(:) :: a, b, c, d, dx, dy, dxy, tmpx, tmpy
!!        
!!        real(rp), allocatable :: rtmp(:)
!!        integer(ip), allocatable :: itmp(:,:)
!!        integer(ip) :: n_data, n_map
!!
!!        if(allocated(ttmap_ang1)) then
!!            ! Reallocate the arrays to make space for the new data
!!            n_data = size(ttmap_ang1)
!!            call mallocate('torstors_newmap [rtmp]', n_data, rtmp)
!!            
!!            rtmp = ttmap_ang1
!!            call mfree('torstors_newmap [ttmap_ang1]', ttmap_ang1)
!!            call mallocate('torstors_newmap [ttmap_ang1]', &
!!                           n_data+d1*d2,  ttmap_ang1)
!!            ttmap_ang1(:n_data) = rtmp
!!            
!!            rtmp = ttmap_ang2
!!            call mfree('torstors_newmap [ttmap_ang2]', ttmap_ang2)
!!            call mallocate('torstors_newmap [ttmap_ang2]', &
!!                           n_data+d1*d2,  ttmap_ang2)
!!            ttmap_ang2(:n_data) = rtmp
!!            
!!            
!!            call mfree('torstors_newmap [rtmp]', rtmp)
!!            n_data = size(ttmap_v)
!!            call mallocate('torstors_newmap [rtmp]', n_data, rtmp)
!!            
!!            rtmp = ttmap_v
!!            call mfree('torstors_newmap [ttmap_v]', ttmap_v)
!!            call mallocate('torstors_newmap [ttmap_v]', &
!!                           n_data+d1*d2,  ttmap_v)
!!            ttmap_v(:n_data) = rtmp
!!            
!!            rtmp = ttmap_vx
!!            call mfree('torstors_newmap [ttmap_vx]', ttmap_vx)
!!            call mallocate('torstors_newmap [ttmap_vx]', &
!!                           n_data+d1*d2,  ttmap_vx)
!!            ttmap_vx(:n_data) = rtmp
!!
!!            rtmp = ttmap_vy
!!            call mfree('torstors_newmap [ttmap_vy]', ttmap_vy)
!!            call mallocate('torstors_newmap [ttmap_vy]', &
!!                           n_data+d1*d2,  ttmap_vy)
!!            ttmap_vy(:n_data) = rtmp
!!
!!            rtmp = ttmap_vxy
!!            call mfree('torstors_newmap [ttmap_vxy]', ttmap_vxy)
!!            call mallocate('torstors_newmap [ttmap_vxy]', &
!!                           n_data+d1*d2,  ttmap_vxy)
!!            ttmap_vxy(:n_data) = rtmp
!!            call mfree('torstors_newmap [rtmp]', rtmp)
!!
!!            n_map = size(ttmap_shape, 2)
!!            call mallocate('torstors_newmap [itmp]', 2, n_map, itmp)
!!            itmp = ttmap_shape
!!            call mfree('torstors_newmap [ttmap_shape]', ttmap_shape)
!!            call mallocate('torstors_newmap [ttmap_shape]', &
!!                           2, n_map+1, ttmap_shape)
!!            ttmap_shape(:,:n_map) = itmp
!!
!!            call mfree('torstors_newmap [itmp]', itmp)
!!        else 
!!            ! First allocation, n_data and n_map are just set for consistency
!!            n_data = 0
!!            n_map = 0
!!            call mallocate('torstors_newmap [ttmap_ang1]', d1*d2,  ttmap_ang1)
!!            call mallocate('torstors_newmap [ttmap_ang2]', d1*d2,  ttmap_ang2)
!!            call mallocate('torstors_newmap [ttmap_v]', d1*d2,  ttmap_v)
!!            call mallocate('torstors_newmap [ttmap_vx]', d1*d2,  ttmap_vx)
!!            call mallocate('torstors_newmap [ttmap_vy]', d1*d2,  ttmap_vy)
!!            call mallocate('torstors_newmap [ttmap_vxy]', d1*d2,  ttmap_vxy)
!!            call mallocate('torstors_newmap [ttmap_shape]', 2, 1, ttmap_shape)
!!        end if
!!
!!        call mallocate('tortor_newmap [a]', max(d1,d2), a)
!!        call mallocate('tortor_newmap [b]', max(d1,d2), b)
!!        call mallocate('tortor_newmap [c]', max(d1,d2), c)
!!        call mallocate('tortor_newmap [d]', max(d1,d2), d)
!!        call mallocate('tortor_newmap [dx]', d1*d2, dx)
!!        call mallocate('tortor_newmap [dy]', d1*d2, dy)
!!        call mallocate('tortor_newmap [dxy]', d1*d2, dxy)
!!
!!        ! This part of the code computes df/dx, df/dy and d^2f/dxdy on the grid.
!!        ! Since we are basically interpolating on a sphere, we extract the 
!!        ! coordinate on a meridian, we interpolate it with a cubic spline, and
!!        ! finally we compute the derivative of this curve at the grid intersection
!!        ! The same is done in the second direction.
!!        ! To compute the mixed derivative we apply the same procedure but using
!!        ! the derivative data (basically we apply the procedure used to compute
!!        ! df/dx but using  df/dy data instead of actual f values.
!!        do i=1, d2
!!            call cyclic_spline(d1, ang1((i-1)*d1+1:i*d1), v((i-1)*d1+1:i*d1), &
!!                               a(1:d1), b(1:d1), c(1:d1), d(1:d1))
!!            dx((i-1)*d1+1:i*d1) = b(1:d1)
!!        end do
!!        
!!        ! df/dy since in this direction data are not contiguous, wa allocate 
!!        ! temporary arrays
!!        call mallocate('tortor_newmap [tmpx]', d2, tmpx)
!!        call mallocate('tortor_newmap [tmpy]', d2, tmpy)
!!        do i=1, d1
!!            ii = 1
!!            do j=i, (d2-1)*d1+i, d2
!!                tmpx(ii) = ang2(j)
!!                tmpy(ii) = v(j)
!!                ii = ii + 1
!!            end do
!!            call cyclic_spline(d2, tmpx, tmpy, &
!!                               a(1:d2), b(1:d2), c(1:d2), d(1:d2))
!!            
!!            ii = 1
!!            do j=i, (d2-1)*d1+i, d2
!!                dy(j) = b(ii)
!!                ii = ii + 1
!!            end do
!!        end do
!!        
!!        ! d^2f/dxdy in this case we use df/dx procedure to exploit data contiguity.
!!        do i=1, d2
!!            call cyclic_spline(d1, ang1((i-1)*d1+1:i*d1), dy((i-1)*d1+1:i*d1), &
!!                               a(1:d1), b(1:d1), c(1:d1), d(1:d1))
!!            dxy((i-1)*d1+1:i*d1) = b(1:d1)
!!        end do
!!        call mfree('tortor_newmap [tmpx]', tmpx)
!!        call mfree('tortor_newmap [tmpy]', tmpy)
!!
!!        ttmap_ang1(n_data+1:) = ang1
!!        ttmap_ang2(n_data+1:) = ang2
!!        ttmap_shape(1,n_map+1) = d1
!!        ttmap_shape(2,n_map+1) = d2
!!        ttmap_v(n_data+1:) = v
!!        ttmap_vx(n_data+1:) = dx
!!        ttmap_vy(n_data+1:) = dy
!!        ttmap_vxy(n_data+1:) = dxy
!!        
!!        call mfree('tortor_newmap [a]', a)
!!        call mfree('tortor_newmap [b]', b)
!!        call mfree('tortor_newmap [c]', c)
!!        call mfree('tortor_newmap [d]', d)
!!        call mfree('tortor_newmap [dx]', dx)
!!        call mfree('tortor_newmap [dy]', dy)
!!        call mfree('tortor_newmap [dxy]', dxy)
!!
!!    end subroutine tortor_newmap
!!
!!    subroutine tortor_potential(V)
!!        !! Compute torsion potential
!!
!!        use mod_utils, only: compute_bicubic_interp
!!
!!        implicit none
!!
!!        type(ommp_bonded_type) :: bds
!!        ! Bonded potential data structure
!!        real(rp), intent(inout) :: V
!!        !! torsion potential, result will be added to V
!!        real(rp) :: thetx, thety, vtt
!!
!!        integer(ip) :: i, j, iprm, ibeg, iend
!!
!!        if(.not. use_tortor) return
!!        
!!        do i=1, ntortor
!!            ! Atoms that defines the two angles
!!            iprm = tortorprm(i)
!!            ibeg = 1
!!            do j=1, iprm-1
!!                ibeg = ibeg + ttmap_shape(1,j)*ttmap_shape(2,j)
!!            end do
!!            iend = ibeg + ttmap_shape(1,iprm)*ttmap_shape(2,iprm) - 1
!!           
!!            thetx = ang_torsion(tortorat(1:4,i))
!!            thety = ang_torsion(tortorat(2:5,i))
!!           
!!            call compute_bicubic_interp(thetx, thety, vtt, &
!!                                        ttmap_shape(1,iprm), &
!!                                        ttmap_shape(2,iprm), &
!!                                        ttmap_ang1(ibeg:iend), &
!!                                        ttmap_ang2(ibeg:iend), &
!!                                        ttmap_v(ibeg:iend), &
!!                                        ttmap_vx(ibeg:iend), &
!!                                        ttmap_vy(ibeg:iend), &
!!                                        ttmap_vxy(ibeg:iend))
!!            V = V + vtt
!!        end do
!!
!!    end subroutine tortor_potential
!!
!!    pure function cos_torsion(idx)
!!        !! Compute the cosine of torsional angle between four atoms specified
!!        !! with indices idx
!!        use mod_mmpol, only: cmm
!!        
!!        implicit none
!!
!!        integer(ip), intent(in) :: idx(4)
!!        real(rp) :: cos_torsion
!!
!!        real(rp), dimension(3) :: a, b, c, d, ab, cd, cb, t, u
!!            
!!        a = cmm(:,idx(1))
!!        b = cmm(:,idx(2))
!!        c = cmm(:,idx(3))
!!        d = cmm(:,idx(4))
!!
!!        ab = b - a
!!        cd = d - c
!!        cb = b - c
!!
!!        t(1) = ab(2)*cb(3) - ab(3)*cb(2)
!!        t(2) = ab(3)*cb(1) - ab(1)*cb(3)
!!        t(3) = ab(1)*cb(2) - ab(2)*cb(1)
!!        t = t / norm2(t)
!!            
!!        u(1) = cb(2)*cd(3) - cb(3)*cd(2)
!!        u(2) = cb(3)*cd(1) - cb(1)*cd(3)
!!        u(3) = cb(1)*cd(2) - cb(2)*cd(1)
!!        u = u / norm2(u)
!!            
!!        cos_torsion = dot_product(u,t)
!!        return 
!!
!!    end function
!!    
!!    pure function ang_torsion(idx)
!!        !! Compute the torsional angle between four atoms specified
!!        !! with indices idx; results are in range [-pi;pi]
!!        use mod_mmpol, only: cmm
!!        
!!        implicit none
!!
!!        integer(ip), intent(in) :: idx(4)
!!        real(rp) :: cos_torsion, ang_torsion
!!
!!        real(rp), dimension(3) :: a, b, c, d, ab, cd, cb, t, u
!!            
!!        a = cmm(:,idx(1))
!!        b = cmm(:,idx(2))
!!        c = cmm(:,idx(3))
!!        d = cmm(:,idx(4))
!!
!!        ab = b - a
!!        cd = d - c
!!        cb = b - c
!!
!!        t(1) = ab(2)*cb(3) - ab(3)*cb(2)
!!        t(2) = ab(3)*cb(1) - ab(1)*cb(3)
!!        t(3) = ab(1)*cb(2) - ab(2)*cb(1)
!!        t = t / norm2(t)
!!            
!!        u(1) = cb(2)*cd(3) - cb(3)*cd(2)
!!        u(2) = cb(3)*cd(1) - cb(1)*cd(3)
!!        u(3) = cb(1)*cd(2) - cb(2)*cd(1)
!!        u = u / norm2(u)
!!            
!!        cos_torsion = dot_product(u,t)
!!        ang_torsion = acos(cos_torsion)
!!        if(dot_product(ab, u) > 0) ang_torsion = - ang_torsion
!!        return 
!!
!!    end function

    subroutine bonded_terminate(bds)
        !! Just terminate every "submodule" in bonded, 
        !! deallocating arrays and disabling the potential terms
        implicit none
    
        type(ommp_bonded_type), intent(inout) :: bds
        ! Bonded potential data structure

        call bond_terminate(bds)
        call angle_terminate(bds)
        call strbnd_terminate(bds)
        call urey_terminate(bds)
        call opb_terminate(bds)
        call pitors_terminate(bds)
        call torsion_terminate(bds)
        call tortor_terminate(bds)
        call angtor_terminate(bds)
        call strtor_terminate(bds)

    end subroutine bonded_terminate
    
    subroutine bond_terminate(bds)
        use mod_memory, only: mfree

        implicit none

        type(ommp_bonded_type), intent(inout) :: bds
        ! Bonded potential data structure
        if( .not. bds%use_bond ) return

        bds%use_bond = .false.
        call mfree('bond_terminate [bondat]', bds%bondat)
        call mfree('bond_terminate [kbond]', bds%kbond)
        call mfree('bond_terminate [l0bond]', bds%l0bond)

    end subroutine bond_terminate
    
    subroutine angle_terminate(bds)
        use mod_memory, only: mfree

        implicit none

        type(ommp_bonded_type), intent(inout) :: bds
        ! Bonded potential data structure
        if( .not. bds%use_angle ) return

        bds%use_angle = .false.
        call mfree('angle_terminate [angleat]', bds%angleat)
        call mfree('angle_terminate [anglety]', bds%anglety)
        call mfree('angle_terminate [angauxat]', bds%angauxat)
        call mfree('angle_terminate [kangle]', bds%kangle)
        call mfree('angle_terminate [eqangle]', bds%eqangle)

    end subroutine angle_terminate
    
    subroutine strbnd_terminate(bds)
        use mod_memory, only: mfree

        implicit none

        type(ommp_bonded_type), intent(inout) :: bds
        ! Bonded potential data structure
        if( .not. bds%use_strbnd ) return

        bds%use_strbnd = .false.
        call mfree('strbnd_terminate [strbndat]', bds%strbndat)
        call mfree('strbnd_terminate [strbndl10]', bds%strbndl10)
        call mfree('strbnd_terminate [strbndl20]', bds%strbndl20)
        call mfree('strbnd_terminate [strbndthet0]', bds%strbndthet0)
        call mfree('strbnd_terminate [strbndk1]', bds%strbndk1)
        call mfree('strbnd_terminate [strbndk2]', bds%strbndk2)

    end subroutine strbnd_terminate
    
    subroutine urey_terminate(bds) 
        use mod_memory, only: mfree

        implicit none

        type(ommp_bonded_type), intent(inout) :: bds
        ! Bonded potential data structure
        if( .not. bds%use_urey ) return
        
        bds%use_urey = .false.
        call mfree('urey_terminate [ureya]',  bds%ureyat)
        call mfree('urey_terminate [kurey]',  bds%kurey)
        call mfree('urey_terminate [l0urey]', bds%l0urey)

    end subroutine urey_terminate
    
    subroutine opb_terminate(bds)
        use mod_memory, only: mfree

        implicit none

        type(ommp_bonded_type), intent(inout) :: bds
        ! Bonded potential data structure
        if( .not. bds%use_opb ) return
        
        bds%use_opb = .false.
        call mfree('opb_terminate [opbat]', bds%opbat)
        call mfree('opb_terminate [kopb]', bds%kopb)

    end subroutine opb_terminate

    subroutine pitors_terminate(bds)
        use mod_memory, only: mfree

        implicit none

        type(ommp_bonded_type), intent(inout) :: bds
        ! Bonded potential data structure
        if( .not. bds%use_pitors ) return
        
        bds%use_pitors = .false.
        call mfree('pitors_terminate [pitorsat]', bds%pitorsat)
        call mfree('p_terminate [kpitors]', bds%kpitors)

    end subroutine pitors_terminate
    
    subroutine torsion_terminate(bds)
        use mod_memory, only: mfree

        implicit none

        type(ommp_bonded_type), intent(inout) :: bds
        ! Bonded potential data structure
        if( .not. bds%use_torsion ) return
        
        bds%use_torsion = .false.
        call mfree('torsion_terminate [torsionat]', bds%torsionat)
        call mfree('torsion_terminate [torsamp]', bds%torsamp)
        call mfree('torsion_terminate [torsphase]', bds%torsphase)
        call mfree('torsion_terminate [torsn]', bds%torsn)

    end subroutine torsion_terminate
    
    subroutine tortor_terminate(bds)
        use mod_memory, only: mfree

        implicit none

        type(ommp_bonded_type), intent(inout) :: bds
        ! Bonded potential data structure
        if( .not. bds%use_tortor ) return
        
        bds%use_tortor = .false.
        call mfree('torsion_terminate [tortorprm]', bds%tortorprm )
        call mfree('torsion_terminate [tortorat]', bds%tortorat)
        call mfree('torsion_terminate [ttmap_shape]', bds%ttmap_shape)
        call mfree('torsion_terminate [ttmap_ang1]', bds%ttmap_ang1)
        call mfree('torsion_terminate [ttmap_ang2]', bds%ttmap_ang2)
        call mfree('torsion_terminate [ttmap_v]', bds%ttmap_v)
        call mfree('torsion_terminate [ttmap_vx]', bds%ttmap_vx)
        call mfree('torsion_terminate [ttmap_vy]', bds%ttmap_vy)
        call mfree('torsion_terminate [ttmap_vxy]', bds%ttmap_vxy)

    end subroutine tortor_terminate
    
    subroutine angtor_terminate(bds)
        use mod_memory, only: mfree

        implicit none

        type(ommp_bonded_type), intent(inout) :: bds
        ! Bonded potential data structure
        if( .not. bds%use_angtor ) return
        
        bds%use_angtor = .false.
        call mfree('angtor_terminate [angtorat]', bds%angtorat)
        call mfree('angtor_terminate [angtork]', bds%angtork)
        call mfree('angtor_terminate [angtor_t]', bds%angtor_t)
        call mfree('angtor_terminate [angtor_a]', bds%angtor_a)

    end subroutine angtor_terminate

    subroutine strtor_terminate(bds)
        use mod_memory, only: mfree

        implicit none

        type(ommp_bonded_type), intent(inout) :: bds
        ! Bonded potential data structure
        if( .not. bds%use_strtor ) return
        
        bds%use_strtor = .false.
        call mfree('strtor_terminate [strtorat]', bds%strtorat)
        call mfree('strtor_terminate [strtork]', bds%strtork)
        call mfree('strtor_terminate [strtor_t]', bds%strtor_t)
        call mfree('strtor_terminate [strtor_b]', bds%strtor_b)

    end subroutine strtor_terminate

end module mod_bonded
