module mod_bonded
    !! Module to handle the bonded part of the FF, it closely follows the 
    !! AMOEBA functional form.
    
    use mod_memory, only: ip, rp
    use mod_topology, only: ommp_topology_type
    use mod_io, only: fatal_error
    
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
        real(rp), allocatable :: kbond(:)    
        !! Force constants for bond terms
        real(rp), allocatable :: l0bond(:)   
        !! Equilibrium lengths for bonds
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
        
        ! Imporoper Torsion
        integer(ip) :: nimptorsion
        integer(ip), allocatable :: imptorsionat(:,:), imptorsn(:,:)
        real(rp), allocatable :: imptorsamp(:,:), imptorsphase(:,:)
        logical :: use_imptorsion = .false.

        ! Torsion-torsion coupling (cmap)
        integer(ip) :: ntortor
        integer(ip), allocatable :: tortorat(:,:), tortorprm(:), ttmap_shape(:,:)
        real(rp), allocatable :: ttmap_ang1(:), ttmap_ang2(:), ttmap_v(:), &
                                 ttmap_vx(:), ttmap_vy(:), ttmap_vxy(:)
        logical :: use_tortor = .false.
    end type ommp_bonded_type

    public :: ommp_bonded_type
    public :: bond_init, bond_potential, bond_geomgrad, bond_terminate
    public :: angle_init, angle_potential, angle_geomgrad, angle_terminate
    public :: urey_init, urey_potential, urey_geomgrad, urey_terminate
    public :: strbnd_init, strbnd_potential, strbnd_geomgrad, strbnd_terminate
    public :: opb_init, opb_potential, opb_geomgrad, opb_terminate
    public :: pitors_init, pitors_potential, pitors_geomgrad, pitors_terminate
    public :: torsion_init, torsion_potential, torsion_geomgrad, &
              torsion_terminate
    public :: imptorsion_init, imptorsion_potential, imptorsion_geomgrad, &
              imptorsion_terminate
    public :: tortor_init, tortor_potential, tortor_geomgrad, &
              tortor_terminate, tortor_newmap
    public :: strtor_init, strtor_potential, strtor_geomgrad, strtor_terminate
    public :: angtor_init, angtor_potential, angtor_geomgrad, angtor_terminate
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
    
    subroutine bond_geomgrad(bds, grad)
        use mod_constants, only : eps_rp
        use mod_jacobian_mat, only: Rij_jacobian

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        !! Bonded potential data structure
        real(rp), intent(inout) :: grad(3,bds%top%mm_atoms)
        !! Gradients of bond stretching terms of potential energy

        integer :: i, ia, ib
        logical :: use_cubic, use_quartic, sk_a, sk_b
        real(rp) :: ca(3), cb(3), J_a(3), J_b(3), l, dl, g

        use_cubic = (abs(bds%bond_cubic) > eps_rp)
        use_quartic = (abs(bds%bond_quartic) > eps_rp)
        
        if(.not. bds%use_bond) return

        if(.not. use_cubic .and. .not. use_quartic) then
            ! This is just a regular harmonic potential
            do i=1, bds%nbond
                ia = bds%bondat(1,i)
                ib = bds%bondat(2,i)
                
                if(bds%top%use_frozen) then
                    sk_a = bds%top%frozen(ia)
                    sk_b = bds%top%frozen(ib)
                    if(sk_a .and. sk_b) cycle
                else
                    sk_a = .false.
                    sk_b = .false.
                end if
                
                ca = bds%top%cmm(:,ia)
                cb = bds%top%cmm(:,ib)
                
                call Rij_jacobian(ca, cb, l, J_a, J_b)
                dl = l - bds%l0bond(i)
                
                g = 2 * bds%kbond(i) * dl
                if(.not. sk_a) grad(:,ia) = grad(:,ia) + J_a * g
                if(.not. sk_b) grad(:,ib) = grad(:,ib) + J_b * g
            end do
        else
            do i=1, bds%nbond
                ia = bds%bondat(1,i)
                ib = bds%bondat(2,i)
                
                if(bds%top%use_frozen) then
                    sk_a = bds%top%frozen(ia)
                    sk_b = bds%top%frozen(ib)
                    if(sk_a .and. sk_b) cycle
                else
                    sk_a = .false.
                    sk_b = .false.
                end if
                
                ca = bds%top%cmm(:,ia)
                cb = bds%top%cmm(:,ib)
                
                call Rij_jacobian(ca, cb, l, J_a, J_b)
                dl = l - bds%l0bond(i)
                
                g = 2 * bds%kbond(i) * dl * (1.0_rp + 3.0/2.0*bds%bond_cubic*dl&
                                             + 2.0*bds%bond_quartic*dl**2)

                if(.not. sk_a) grad(:,ia) = grad(:,ia) + J_a * g
                if(.not. sk_b) grad(:,ib) = grad(:,ib) + J_b * g
            end do
        end if
        
    end subroutine bond_geomgrad

    subroutine angle_init(bds, n)
        !! Initialize arrays used in calculation of angle bending functions

        use mod_memory, only: mallocate

        implicit none

        type(ommp_bonded_type), intent(inout) :: bds
        ! Bonded potential data structure
        integer(ip) :: n
        !! Number of angle bending functions in the potential
        !! energy of the system

        if( n < 1 ) return
        bds%use_angle = .true.

        call mallocate('angle_init [angleat]', 3, n, bds%angleat)
        call mallocate('angle_init [anglety]', n, bds%anglety)
        call mallocate('angle_init [angauxat]', n, bds%angauxat)
        call mallocate('angle_init [kangle]', n, bds%kangle)
        call mallocate('angle_init [eqangle]', n, bds%eqangle)
        
        bds%nangle = n
        bds%angauxat = 0
        bds%angle_cubic = 0.0_rp
        bds%angle_quartic = 0.0_rp
        bds%angle_pentic = 0.0_rp
        bds%angle_sextic = 0.0_rp

    end subroutine angle_init

    subroutine angle_potential(bds, V)
        !! Compute angle-bending terms of the potential energy function.   
        !! Simple angle terms are computed according to the formula:
        !! \[U_{angle} = \sum_i k_i \Delta \theta_i^2 \large(1 +  
        !!  \sum_{j=1}^4 k^{(j+2)} \Delta \theta_i^j \large)\]
        !! \[\Delta \theta_i = \theta_i - \theta^{(eq)}_i\]    
        !! Out-of plane angle are more complex. First, central atom has to be
        !! a trigonal center, the other two atoms together with the auxliary 
        !! atom (that is the remaining one connected to the trigonal center) 
        !! define the projection plane. During the first run the auxiliary atom
        !! is found and saved.
        !! Then, the trigonal center is projected on the plane defined by the 
        !! other three atoms, and the angle is the one defined by the projection
        !! (which is the vertex, and the other two atoms -- the auxiliary is
        !! excluded). Then the same formula used for simple angle terms is used.

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        ! Bonded potential data structure
        real(rp), intent(inout) :: V
        !! Bond potential, result will be added to V
        
        integer(ip) :: i
        real(rp) :: l1, l2, dr1(3), dr2(3), thet, d_theta
        real(rp), dimension(3) :: v_dist, plv1, plv2, pln, a, b, c, prj_b, aux

        if(.not. bds%use_angle) return
        
        do i=1, bds%nangle
            if(bds%anglety(i) == OMMP_ANG_SIMPLE .or. &
               bds%anglety(i) == OMMP_ANG_H0 .or. &
               bds%anglety(i) == OMMP_ANG_H1 .or. &
               bds%anglety(i) == OMMP_ANG_H2) then
                dr1 = bds%top%cmm(:, bds%angleat(1,i)) - bds%top%cmm(:, bds%angleat(2,i))
                dr2 = bds%top%cmm(:, bds%angleat(3,i)) - bds%top%cmm(:, bds%angleat(2,i))
                l1 = sqrt(dot_product(dr1, dr1))
                l2 = sqrt(dot_product(dr2, dr2))

                thet = acos(dot_product(dr1, dr2)/(l1*l2))
                
                d_theta = thet-bds%eqangle(i) 
                
                V = V + bds%kangle(i) * d_theta**2 * (1.0 + bds%angle_cubic*d_theta &
                    + bds%angle_quartic*d_theta**2 + bds%angle_pentic*d_theta**3 &
                    + bds%angle_sextic*d_theta**4)

            else if(bds%anglety(i) == OMMP_ANG_INPLANE .or. &
                    bds%anglety(i) == OMMP_ANG_INPLANE_H0 .or. &
                    bds%anglety(i) == OMMP_ANG_INPLANE_H1) then
                
                a = bds%top%cmm(:, bds%angleat(1,i))
                b = bds%top%cmm(:, bds%angleat(2,i)) !! Trigonal center
                c = bds%top%cmm(:, bds%angleat(3,i))

                aux = bds%top%cmm(:, bds%angauxat(i))
                plv1 = a - aux
                plv2 = c - aux
                pln(1) = plv1(2)*plv2(3) - plv1(3)*plv2(2)
                pln(2) = plv1(3)*plv2(1) - plv1(1)*plv2(3)
                pln(3) = plv1(1)*plv2(2) - plv1(2)*plv2(1)
                !! Normal vector of the projection plane
                pln = pln / sqrt(dot_product(pln, pln))

                v_dist = b - aux
                prj_b = b - dot_product(v_dist, pln) * pln 

                dr1 = bds%top%cmm(:, bds%angleat(1,i)) - prj_b
                dr2 = bds%top%cmm(:, bds%angleat(3,i)) - prj_b
                l1 = sqrt(dot_product(dr1, dr1))
                l2 = sqrt(dot_product(dr2, dr2))

                thet = acos(dot_product(dr1, dr2)/(l1*l2))
                
                d_theta = thet-bds%eqangle(i) 
                
                V = V + bds%kangle(i) * d_theta**2 * (1.0 + bds%angle_cubic*d_theta &
                    + bds%angle_quartic*d_theta**2 + bds%angle_pentic*d_theta**3 &
                    + bds%angle_sextic*d_theta**4)
            end if
        end do
    end subroutine angle_potential
    
    subroutine angle_geomgrad(bds, grad)
        use mod_jacobian_mat, only: simple_angle_jacobian, &
                                    inplane_angle_jacobian

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        !! Bonded potential data structure
        real(rp), intent(inout) :: grad(3,bds%top%mm_atoms)
        !! Gradients of bond stretching terms of potential energy
        
        real(rp) :: a(3), b(3), c(3), Ja(3), Jb(3), Jc(3), Jx(3), g, thet, &
                    d_theta, aux(3)
        integer(ip) :: i
        logical :: sk_a, sk_b, sk_c, sk_x

        if(.not. bds%use_angle) return
        
        do i=1, bds%nangle
            if(bds%anglety(i) == OMMP_ANG_SIMPLE .or. &
               bds%anglety(i) == OMMP_ANG_H0 .or. &
               bds%anglety(i) == OMMP_ANG_H1 .or. &
               bds%anglety(i) == OMMP_ANG_H2) then
                if(bds%top%use_frozen) then
                    sk_a = bds%top%frozen(bds%angleat(1,i))
                    sk_b = bds%top%frozen(bds%angleat(2,i))
                    sk_c = bds%top%frozen(bds%angleat(3,i))
                    if(sk_a .and. sk_b .and. sk_c) cycle
                else
                    sk_a = .false.
                    sk_b = .false.
                    sk_c = .false.
                end if

                a = bds%top%cmm(:, bds%angleat(1,i)) 
                b = bds%top%cmm(:, bds%angleat(2,i))
                c = bds%top%cmm(:, bds%angleat(3,i))
                call simple_angle_jacobian(a, b, c, thet, Ja, Jb, Jc)
                d_theta = thet - bds%eqangle(i) 
           
                g = bds%kangle(i) * d_theta * (2.0 &
                                               + 3.0 * bds%angle_cubic * d_theta &
                                               + 4.0 * bds%angle_quartic * d_theta**2 &
                                               + 5.0 * bds%angle_pentic * d_theta**3 &
                                               + 6.0 * bds%angle_sextic * d_theta**4)

                if(.not. sk_a) grad(:,bds%angleat(1,i)) = grad(:,bds%angleat(1,i)) + g * Ja
                if(.not. sk_b) grad(:,bds%angleat(2,i)) = grad(:,bds%angleat(2,i)) + g * Jb
                if(.not. sk_c) grad(:,bds%angleat(3,i)) = grad(:,bds%angleat(3,i)) + g * Jc
            
            else if(bds%anglety(i) == OMMP_ANG_INPLANE .or. &
                    bds%anglety(i) == OMMP_ANG_INPLANE_H0 .or. &
                    bds%anglety(i) == OMMP_ANG_INPLANE_H1) then
                
                if(bds%top%use_frozen) then
                    sk_a = bds%top%frozen(bds%angleat(1,i))
                    sk_b = bds%top%frozen(bds%angleat(2,i))
                    sk_c = bds%top%frozen(bds%angleat(3,i))
                    sk_x = bds%top%frozen(bds%angauxat(i))
                    if(sk_a .and. sk_b .and. sk_c .and. sk_x) cycle
                else
                    sk_a = .false.
                    sk_b = .false.
                    sk_c = .false.
                    sk_x = .false.
                end if
                
                a = bds%top%cmm(:, bds%angleat(1,i))
                b = bds%top%cmm(:, bds%angleat(2,i)) !! Trigonal center
                c = bds%top%cmm(:, bds%angleat(3,i))

                aux = bds%top%cmm(:, bds%angauxat(i))
                
                call inplane_angle_jacobian(a, b, c, aux, thet, Ja, Jb, Jc, Jx)
                d_theta = thet - bds%eqangle(i) 
                g = bds%kangle(i) * d_theta * (2.0 &
                                               + 3.0 * bds%angle_cubic * d_theta &
                                               + 4.0 * bds%angle_quartic * d_theta**2 &
                                               + 5.0 * bds%angle_pentic * d_theta**3 &
                                               + 6.0 * bds%angle_sextic * d_theta**4)
               
                if(.not. sk_a) grad(:,bds%angleat(1,i)) = grad(:,bds%angleat(1,i)) + g * Ja
                if(.not. sk_b) grad(:,bds%angleat(2,i)) = grad(:,bds%angleat(2,i)) + g * Jb
                if(.not. sk_c) grad(:,bds%angleat(3,i)) = grad(:,bds%angleat(3,i)) + g * Jc
                if(.not. sk_x) grad(:,bds%angauxat(i)) = grad(:,bds%angauxat(i)) + g * Jx

            end if
        end do
    end subroutine angle_geomgrad
   
    subroutine strbnd_init(bds, n)
        !! Initialize arrays for calculation of stretch-bend cross term 
        !! potential

        use mod_memory, only: mallocate

        implicit none

        type(ommp_bonded_type), intent(inout) :: bds
        ! Bonded potential data structure
        integer(ip) :: n
        !! Number of stretch-bend functions in the potential
        !! energy of the system

        if( n < 1 ) return
        bds%use_strbnd = .true.

        call mallocate('strbnd_init [strbndat]', 3, n, bds%strbndat)
        call mallocate('strbnd_init [strbndl10]', n, bds%strbndl10)
        call mallocate('strbnd_init [strbndl20]', n, bds%strbndl20)
        call mallocate('strbnd_init [strbndthet0]', n, bds%strbndthet0)
        call mallocate('strbnd_init [strbndk1]', n, bds%strbndk1)
        call mallocate('strbnd_init [strbndk2]', n, bds%strbndk2)
        bds%nstrbnd = n

    end subroutine strbnd_init

    subroutine strbnd_potential(bds, V)
        !! Compute the stretch-bend cross term potential.   
        !! Those terms are computed according the following formula:
        !! \[U_{bond/angle} = (k_i \Delta l_i + k_j \Delta l_j) 
        !! \Delta \theta_{ij} \]
        !! where \(\theta_{ij}\) is the angle delimited by the bond \(i\) and 
        !! \(j\).   
        !! The force constants \(k_i\) and \(k_j\) are explicitely defined in
        !! the FF, while the equilibrium values are the same as for stretching
        !! and bending terms.

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        ! Bonded potential data structure
        real(rp), intent(inout) :: V
        !! Stretch-bend cross term potential, result will be added to V

        integer(ip) :: i
        real(rp) :: d_l1, d_l2, d_thet, dr1(3), dr2(3), l1, l2, thet
        
        if(.not. bds%use_strbnd) return

        do i=1, bds%nstrbnd
            dr1 = bds%top%cmm(:, bds%strbndat(2,i)) - &
                  bds%top%cmm(:, bds%strbndat(1,i))
            l1 = norm2(dr1)
            d_l1 = l1 - bds%strbndl10(i)
            
            dr2 = bds%top%cmm(:, bds%strbndat(2,i)) - &
                  bds%top%cmm(:, bds%strbndat(3,i))
            l2 = norm2(dr2)
            d_l2 = l2 - bds%strbndl20(i)

            thet = acos(dot_product(dr1, dr2)/(l1*l2))
            d_thet = thet - bds%strbndthet0(i) 
            
            V = V + (d_l1*bds%strbndk1(i) + d_l2*bds%strbndk2(i)) * d_thet
        end do
    end subroutine strbnd_potential
    
    subroutine strbnd_geomgrad(bds, grad)
        use mod_jacobian_mat, only: Rij_jacobian, simple_angle_jacobian

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        ! Bonded potential data structure
        real(rp), intent(inout) :: grad(3,bds%top%mm_atoms)
        !! Gradients of bond stretching terms of potential energy

        integer(ip) :: i, ia, ib, ic
        real(rp) :: d_l1, d_l2, d_thet, l1, l2, thet, g1, g2, g3
        real(rp), dimension(3) :: a, b, c, &
                                  J1_a, J1_b, &
                                  J2_b, J2_c, &
                                  J3_a, J3_b, J3_c
        logical :: sk_a, sk_b, sk_c
        
        if(.not. bds%use_strbnd) return

        do i=1, bds%nstrbnd
            ia = bds%strbndat(1,i)
            ib = bds%strbndat(2,i)
            ic = bds%strbndat(3,i)
            
            if(bds%top%use_frozen) then
                sk_a = bds%top%frozen(ia)
                sk_b = bds%top%frozen(ib)
                sk_c = bds%top%frozen(ic)
                if(sk_a .and. sk_b .and. sk_c) cycle
            else
                sk_a = .false.
                sk_b = .false.
                sk_c = .false.
            end if

            a = bds%top%cmm(:, ia)
            b = bds%top%cmm(:, ib)
            c = bds%top%cmm(:, ic)

            call Rij_jacobian(a, b, l1, J1_a, J1_b)
            call Rij_jacobian(b, c, l2, J2_b, J2_c)
            call simple_angle_jacobian(a, b, c, thet, J3_a, J3_b, J3_c)
            
            d_l1 = l1 - bds%strbndl10(i)
            d_l2 = l2 - bds%strbndl20(i)
            d_thet = thet - bds%strbndthet0(i) 
           
            g1 = bds%strbndk1(i) * d_thet
            g2 = bds%strbndk2(i) * d_thet
            g3 = bds%strbndk1(i) * d_l1 + bds%strbndk2(i) * d_l2

            if(.not. sk_a) grad(:,ia) = grad(:,ia) + J1_a * g1 + J3_a * g3
            if(.not. sk_b) grad(:,ib) = grad(:,ib) + J1_b * g1 + J2_b * g2 + J3_b * g3 
            if(.not. sk_c) grad(:,ic) = grad(:,ic) + J2_c * g2 + J3_c * g3
        end do

    end subroutine strbnd_geomgrad
    
    subroutine urey_init(bds, n) 
        !! Initialize Urey-Bradley potential arrays

        use mod_memory, only: mallocate

        implicit none

        type(ommp_bonded_type), intent(inout) :: bds
        ! Bonded potential data structure
        integer(ip) :: n
        !! Number of Urey-Bradley functions in the potential
        !! energy of the system
        
        if( n < 1 ) return
        bds%use_urey = .true.

        call mallocate('urey_init [ureya]', 2, n, bds%ureyat)
        call mallocate('urey_init [kurey]', n, bds%kurey)
        call mallocate('urey_init [l0urey]', n, bds%l0urey)
        bds%nurey = n
        bds%urey_cubic = 0.0_rp
        bds%urey_quartic = 0.0_rp

    end subroutine urey_init

    subroutine urey_potential(bds, V)
        !! Compute the Urey-Bradley potential.  
        !! This is basically a virtual bond, with its stretching harminic 
        !! potential that connect two otherwise un-connected bonds. The same
        !! potential formula used for normal stretching is used.

        use mod_constants, only : eps_rp

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        ! Bonded potential data structure
        real(rp), intent(inout) :: V
        !! Urey-Bradley potential, result will be added to V

        integer :: i
        logical :: use_cubic, use_quartic
        real(rp) :: dr(3), l, dl, dl2
        
        if(.not. bds%use_urey) return

        use_cubic = (abs(bds%urey_cubic) > eps_rp)
        use_quartic = (abs(bds%urey_quartic) > eps_rp)

        if(.not. use_cubic .and. .not. use_quartic) then
            ! This is just a regular harmonic potential
            do i=1, bds%nurey
                dr = bds%top%cmm(:,bds%ureyat(1,i)) - &
                     bds%top%cmm(:,bds%ureyat(2,i))
                l = sqrt(dot_product(dr, dr))
                dl = l - bds%l0urey(i)
                V = V + bds%kurey(i) * dl * dl
            end do
        else
            do i=1, bds%nurey
                dr = bds%top%cmm(:,bds%ureyat(1,i)) - &
                     bds%top%cmm(:,bds%ureyat(2,i))
                l = sqrt(dot_product(dr, dr))
                dl = l - bds%l0urey(i)
                dl2 = dl * dl

                V = V + bds%kurey(i)*dl2 * (1.0_rp + bds%urey_cubic*dl + &
                                            bds%urey_quartic*dl2)
            end do
        end if
    end subroutine urey_potential
    
    subroutine urey_geomgrad(bds, grad)
        use mod_constants, only : eps_rp
        use mod_jacobian_mat, only: Rij_jacobian

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        !! Bonded potential data structure
        real(rp), intent(inout) :: grad(3,bds%top%mm_atoms)
        !! Gradients of bond stretching terms of potential energy

        integer :: i, ia, ib
        logical :: use_cubic, use_quartic, sk_a, sk_b
        real(rp) :: l, dl, J_a(3), J_b(3), g
        
        if(.not. bds%use_urey) return

        use_cubic = (abs(bds%urey_cubic) > eps_rp)
        use_quartic = (abs(bds%urey_quartic) > eps_rp)

        if(.not. use_cubic .and. .not. use_quartic) then
            ! This is just a regular harmonic potential
            do i=1, bds%nurey
                ia = bds%ureyat(1,i)
                ib = bds%ureyat(2,i)

                if(bds%top%use_frozen) then
                    sk_a = bds%top%frozen(ia)
                    sk_b = bds%top%frozen(ib)
                    if(sk_a .and. sk_b) cycle
                else
                    sk_a = .false.
                    sk_b = .false.
                end if

                call Rij_jacobian(bds%top%cmm(:,ia), &
                                  bds%top%cmm(:,ib), &
                                  l, J_a, J_b)
                dl = l - bds%l0urey(i)
                g = 2 * bds%kurey(i) * dl
                if(.not. sk_a) grad(:,ia) = grad(:,ia) + J_a * g
                if(.not. sk_b) grad(:,ib) = grad(:,ib) + J_b * g
            end do
        else
            do i=1, bds%nurey
                ia = bds%ureyat(1,i)
                ib = bds%ureyat(2,i)

                if(bds%top%use_frozen) then
                    sk_a = bds%top%frozen(ia)
                    sk_b = bds%top%frozen(ib)
                    if(sk_a .and. sk_b) cycle
                else
                    sk_a = .false.
                    sk_b = .false.
                end if

                call Rij_jacobian(bds%top%cmm(:,ia), &
                                  bds%top%cmm(:,ib), &
                                  l, J_a, J_b)
                dl = l - bds%l0urey(i)
                g = 2 * bds%kurey(i) * dl * (1.0 &
                                             + 3.0/2.0 * bds%urey_cubic*dl &
                                             + 2.0 * bds%urey_quartic*dl**2)

                if(.not. sk_a) grad(:,ia) = grad(:,ia) + J_a * g
                if(.not. sk_b) grad(:,ib) = grad(:,ib) + J_b * g
            end do
        end if
    end subroutine urey_geomgrad

    subroutine opb_init(bds, n, opbtype)
        !! Initialize arrays for out-of-plane bending potential calculation.   
        !! @todo Currently only Allinger functional form is supported 
        use mod_io, only: ommp_message
        use mod_constants, only: OMMP_VERBOSE_LOW
        use mod_memory, only: mallocate

        implicit none

        type(ommp_bonded_type), intent(inout) :: bds
        ! Bonded potential data structure
        integer(ip) :: n
        !! Number of out of plane Bending functions in the potential
        !! energy of the system
        character(len=*) :: opbtype

        select case(opbtype)
            case('allinger')
                continue
            case('w-d-c')
                call fatal_error('Out-of-plane bend W-D-C is not implemented')
            case default
                call ommp_message("Found OPB type: '"//opbtype//"'", OMMP_VERBOSE_LOW)
                call fatal_error('Out-of-plane type specified is not understood')
        end select

        if( n < 1 ) return
        bds%use_opb = .true.

        call mallocate('opb_init [opbat]', 4, n, bds%opbat)
        call mallocate('opb_init [kopb]', n, bds%kopb)
        bds%nopb = n

    end subroutine opb_init

    subroutine opb_potential(bds, V)
        !! Computes the out-of-plane bending potential.  
        !! With Allinger formula: similarly to in plane angles, here we are 
        !! considering a trigonal center, where D is the central atom and 
        !! A, B, C are connected to D. Allinger formula consider the angle 
        !! between vector \(\vec{AD}\) and the normal vector of plane ABC, 
        !! using \(\frac{\pi}{2}\) as implicit equilibrium value. The formula
        !! for this potential term is:
        !! \[U_{out-of-plane} = \sum_i k_i \chi_i^2 \large(1 + 
        !! \sum_{j=1}^4 k^{(j+2)} \chi_i^j \large) \]

        use mod_constants, only : pi

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        ! Bonded potential data structure
        real(rp), intent(inout) :: V
        !! out-of-plane potential, result will be added to V
        real(rp), dimension(3) :: a, b, c, d, plv1, plv2, pln, vad
        real(rp) :: lpln, lvad, thet, thet2, thet3, thet4
        integer(ip) :: i

        if(.not. bds%use_opb) return
        
        do i=1, bds%nopb
            ! A* -- D -- C
            !       |
            !       B 
            a = bds%top%cmm(:,bds%opbat(2,i))
            d = bds%top%cmm(:,bds%opbat(1,i))
            c = bds%top%cmm(:,bds%opbat(3,i))
            b = bds%top%cmm(:,bds%opbat(4,i))

            ! Compute the normal vector of the plane
            plv1 = a - b
            plv2 = a - c
            pln(1) = plv1(2)*plv2(3) - plv1(3)*plv2(2)
            pln(2) = plv1(3)*plv2(1) - plv1(1)*plv2(3)
            pln(3) = plv1(1)*plv2(2) - plv1(2)*plv2(1)
            lpln = norm2(pln)

            ! Vector from A to D
            vad = a - d
            lvad = norm2(vad)

            thet = abs(pi/2.0 - acos(dot_product(vad, pln)/(lvad*lpln)))
            thet2 = thet*thet
            thet3 = thet2*thet
            thet4 = thet3*thet
            V = V +  bds%kopb(i) * thet2 * (1 + bds%opb_cubic*thet &
                + bds%opb_quartic*thet2 + bds%opb_pentic*thet3 &
                + bds%opb_sextic*thet4)
        end do
    end subroutine opb_potential
    
    subroutine opb_geomgrad(bds, grad)
        use mod_jacobian_mat, only: opb_angle_jacobian

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        ! Bonded potential data structure
        real(rp), intent(inout) :: grad(3,bds%top%mm_atoms)
        !! Gradients of bond stretching terms of potential energy
        real(rp) :: thet, g, J_a(3), J_b(3), J_c(3), J_d(3)
        integer(ip) :: i, ia, ib, ic, id
        logical :: sk_a, sk_b, sk_c, sk_d
        
        if(.not. bds%use_opb) return
        
        do i=1, bds%nopb
            ia = bds%opbat(2,i)
            ib = bds%opbat(4,i)
            ic = bds%opbat(3,i)
            id = bds%opbat(1,i) 

            if(bds%top%use_frozen) then
                sk_a = bds%top%frozen(ia)
                sk_b = bds%top%frozen(ib)
                sk_c = bds%top%frozen(ic)
                sk_d = bds%top%frozen(id)
                if(sk_a .and. sk_b .and. sk_c .and. sk_d) cycle
            else
                sk_a = .false.
                sk_b = .false.
                sk_c = .false.
                sk_d = .false.
            end if

            call opb_angle_jacobian(bds%top%cmm(:,ia), & 
                                    bds%top%cmm(:,ib), &
                                    bds%top%cmm(:,ic), &
                                    bds%top%cmm(:,id), &
                                    thet, J_a, J_b, J_c, J_d)
            
            g = bds%kopb(i) * thet * (2.0 + 3.0*bds%opb_cubic*thet &
                + 4.0*bds%opb_quartic*thet**2 + 5.0*bds%opb_pentic*thet**3 &
                + 6.0*bds%opb_sextic*thet**4)

            if(.not. sk_a) grad(:,ia) = grad(:,ia) + g * J_a
            if(.not. sk_b) grad(:,ib) = grad(:,ib) + g * J_b
            if(.not. sk_c) grad(:,ic) = grad(:,ic) + g * J_c
            if(.not. sk_d) grad(:,id) = grad(:,id) + g * J_d
        end do
    end subroutine opb_geomgrad
    
    subroutine pitors_init(bds, n)
        !! Initialize arrays needed to compute pi-torsion potential

        use mod_memory, only: mallocate

        implicit none

        type(ommp_bonded_type), intent(inout) :: bds
        ! Bonded potential data structure
        integer(ip) :: n
        !! Number of out of plane pi-torsion functions in the potential
        !! enerpgy of the system
        
        if( n < 1 ) return

        bds%use_pitors = .true.

        call mallocate('pitors_init [pitorsat]', 6, n, bds%pitorsat)
        call mallocate('pitors_init [kpitors]', n, bds%kpitors)
        bds%npitors = n

    end subroutine pitors_init

    subroutine pitors_potential(bds, V)
        !! Compute pi-torsion terms of the potential.  
        !! This potential is defined on a \(\pi\)-system, and uses the 
        !! coordinates of six atoms A...F the central "double" bond is A-B, then
        !! C and D are connected to A while E and F are connected to B. So two
        !! plane ACD and BEF are defined. The potential is computed using the 
        !! dihedral angle of the normal vector of those two planes, connected 
        !! by segment A-B (\(\theta\)).  
        !! The formula used is:
        !! \[U_{\pi-torsion} = \sum_i k_i \large(1 + cos(2\theta-\pi) \large)\]

        use mod_constants, only : pi

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        ! Bonded potential data structure
        real(rp), intent(inout) :: V
        !! pi-torsion potential, result will be added to V
        real(rp), dimension(3) :: a, b, c, d, e, f, u, t, cd, plv1, plv2, pln1, pln2
        real(rp) :: thet, costhet
        integer(ip) :: i

        if(.not. bds%use_pitors) return
        
        do i=1, bds%npitors
            !
            !  2(c)        5(e)         a => 1
            !   \         /             b => 4
            !    1(a) -- 4(b)  
            !   /         \
            !  3(d)        6(f)
            
            ! Atoms that defines the two planes
            a = bds%top%cmm(:,bds%pitorsat(1,i))
            c = bds%top%cmm(:,bds%pitorsat(2,i))
            d = bds%top%cmm(:,bds%pitorsat(3,i))

            b = bds%top%cmm(:,bds%pitorsat(4,i))
            e = bds%top%cmm(:,bds%pitorsat(5,i))
            f = bds%top%cmm(:,bds%pitorsat(6,i))


            ! Compute the normal vector of the first plane
            plv1 = d - b
            plv2 = c - b
            pln1(1) = plv1(2)*plv2(3) - plv1(3)*plv2(2)
            pln1(2) = plv1(3)*plv2(1) - plv1(1)*plv2(3)
            pln1(3) = plv1(1)*plv2(2) - plv1(2)*plv2(1)

            ! Compute the normal vector of the second plane
            plv1 = f - a
            plv2 = e - a
            pln2(1) = plv1(2)*plv2(3) - plv1(3)*plv2(2)
            pln2(2) = plv1(3)*plv2(1) - plv1(1)*plv2(3)
            pln2(3) = plv1(1)*plv2(2) - plv1(2)*plv2(1)

            cd = b - a

            t(1) = pln1(2)*cd(3) - pln1(3)*cd(2)
            t(2) = pln1(3)*cd(1) - pln1(1)*cd(3)
            t(3) = pln1(1)*cd(2) - pln1(2)*cd(1)
            t = t / norm2(t)
            
            u(1) = cd(2)*pln2(3) - cd(3)*pln2(2)
            u(2) = cd(3)*pln2(1) - cd(1)*pln2(3)
            u(3) = cd(1)*pln2(2) - cd(2)*pln2(1)
            u = u / norm2(u)
            
            costhet = dot_product(u,t)
            
            thet = acos(costhet)
            
            V = V +  bds%kpitors(i) * (1 + cos(2.0*thet-pi))
        end do

    end subroutine pitors_potential
    
    subroutine pitors_geomgrad(bds, grad)
        use mod_jacobian_mat, only: pitors_angle_jacobian
        use mod_constants, only : pi

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        ! Bonded potential data structure
        real(rp), intent(inout) :: grad(3,bds%top%mm_atoms)
        !! improper torsion potential, result will be added to V
        real(rp) :: thet, g, J_a(3), J_b(3), J_c(3), J_d(3), J_e(3), J_f(3)
        integer(ip) :: i, ia, ib, ic, id, ie, if_
        logical :: sk_a, sk_b, sk_c, sk_d, sk_e, sk_f

        if(.not. bds%use_pitors) return
        
        do i=1, bds%npitors
            ia = bds%pitorsat(1,i)
            ic = bds%pitorsat(2,i)
            id = bds%pitorsat(3,i)
            ib = bds%pitorsat(4,i)
            ie = bds%pitorsat(5,i)
            if_ = bds%pitorsat(6,i)
            
            if(bds%top%use_frozen) then
                sk_a = bds%top%frozen(ia)
                sk_b = bds%top%frozen(ib)
                sk_c = bds%top%frozen(ic)
                sk_d = bds%top%frozen(id)
                sk_e = bds%top%frozen(ie)
                sk_f = bds%top%frozen(if_)
                if(sk_a .and. sk_b .and. sk_c .and. sk_d .and. sk_e .and. sk_f) cycle
            else
                sk_a = .false.
                sk_b = .false.
                sk_c = .false.
                sk_d = .false.
                sk_e = .false.
                sk_f = .false.
            end if

            call pitors_angle_jacobian(bds%top%cmm(:,ia), &
                                       bds%top%cmm(:,ib), &
                                       bds%top%cmm(:,ic), &
                                       bds%top%cmm(:,id), &
                                       bds%top%cmm(:,ie), &
                                       bds%top%cmm(:,if_), &
                                       thet, J_a, J_b, J_c, J_d, J_e, J_f)

            g = -2.0 * bds%kpitors(i) * sin(2.0*thet-pi)

            if(.not. sk_a) grad(:,ia) = grad(:,ia) + g * J_a
            if(.not. sk_b) grad(:,ib) = grad(:,ib) + g * J_b
            if(.not. sk_c) grad(:,ic) = grad(:,ic) + g * J_c
            if(.not. sk_d) grad(:,id) = grad(:,id) + g * J_d
            if(.not. sk_e) grad(:,ie) = grad(:,ie) + g * J_e
            if(.not. sk_f) grad(:,if_) = grad(:,if_) + g * J_f
        end do
    end subroutine pitors_geomgrad
    
    subroutine torsion_init(bds, n)
        !! Initialize torsion potential arrays

        use mod_memory, only: mallocate

        implicit none

        type(ommp_bonded_type), intent(inout) :: bds
        ! Bonded potential data structure
        integer(ip) :: n
        !! Number of torsion functions in the potential
        !! energy of the system
        
        if( n < 1 ) return
        bds%use_torsion = .true.

        call mallocate('torsion_init [torsionat]', 4, n, bds%torsionat)
        call mallocate('torsion_init [torsamp]', 6, n, bds%torsamp)
        call mallocate('torsion_init [torsphase]', 6, n, bds%torsphase)
        call mallocate('torsion_init [torsn]', 6, n, bds%torsn)

        bds%ntorsion = n

    end subroutine torsion_init

    subroutine torsion_potential(bds, V)
        !! Compute torsion potential
        use mod_constants, only: pi, eps_rp

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        ! Bonded potential data structure
        real(rp), intent(inout) :: V
        !! torsion potential, result will be added to V
        real(rp) :: thet, costhet
        integer(ip) :: i, j
        
        if(.not. bds%use_torsion) return

        do i=1, bds%ntorsion
            ! Atoms that defines the dihedral angle
            costhet = cos_torsion(bds%top, bds%torsionat(:,i))
            
            if(costhet + 1.0 <= eps_rp) then
                thet = pi
            else
                thet = acos(costhet)
            end if

            do j=1, 6
                if(bds%torsn(j,i) < 1) exit
                V = V + bds%torsamp(j,i) * (1+cos(real(bds%torsn(j,i))*thet &
                                            - bds%torsphase(j,i)))
            end do
        end do

    end subroutine torsion_potential
    
    subroutine torsion_geomgrad(bds, grad)
        !! Compute torsion potential
        use mod_jacobian_mat, only: torsion_angle_jacobian

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        ! Bonded potential data structure
        real(rp), intent(inout) :: grad(3,bds%top%mm_atoms)
        !! Gradients of bond stretching terms of potential energy
        real(rp) :: thet, g, J_a(3), J_b(3), J_c(3), J_d(3)
        integer(ip) :: i, j, ia, ib, ic, id
        logical :: sk_a, sk_b, sk_c, sk_d
        
        if(.not. bds%use_torsion) return

        do i=1, bds%ntorsion
            ia = bds%torsionat(1,i)
            ib = bds%torsionat(2,i)
            ic = bds%torsionat(3,i)
            id = bds%torsionat(4,i) 

            if(bds%top%use_frozen) then
                sk_a = bds%top%frozen(ia)
                sk_b = bds%top%frozen(ib)
                sk_c = bds%top%frozen(ic)
                sk_d = bds%top%frozen(id)
                if(sk_a .and. sk_b .and. sk_c .and. sk_d) cycle
            else
                sk_a = .false.
                sk_b = .false.
                sk_c = .false.
                sk_d = .false.
            end if

            call torsion_angle_jacobian(bds%top%cmm(:,ia), &
                                        bds%top%cmm(:,ib), &
                                        bds%top%cmm(:,ic), &
                                        bds%top%cmm(:,id), &
                                        thet, J_a, J_b, J_c, J_d)
            
            do j=1, 6
                if(bds%torsn(j,i) < 1) exit
                g = -real(bds%torsn(j,i)) * sin(real(bds%torsn(j,i))* thet &
                                                - bds%torsphase(j,i)) &
                    * bds%torsamp(j,i)
                if(.not. sk_a) grad(:, ia) = grad(:, ia) + J_a * g
                if(.not. sk_b) grad(:, ib) = grad(:, ib) + J_b * g
                if(.not. sk_c) grad(:, ic) = grad(:, ic) + J_c * g
                if(.not. sk_d) grad(:, id) = grad(:, id) + J_d * g
            end do
        end do

    end subroutine torsion_geomgrad
    
    subroutine imptorsion_potential(bds, V)
        !! Compute torsion potential
        use mod_constants, only: pi, eps_rp

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        ! Bonded potential data structure
        real(rp), intent(inout) :: V
        !! improper torsion potential, result will be added to V
        real(rp) :: thet, costhet
        integer(ip) :: i, j
        
        if(.not. bds%use_imptorsion) return
        
        do i=1, bds%nimptorsion
            ! Atoms that defines the dihedral angle
            costhet = cos_torsion(bds%top, bds%imptorsionat(:,i))
            
            if(costhet + 1.0 <= eps_rp) then
                thet = pi
            else
                thet = acos(costhet)
            end if
            
            do j=1, 3
                if(bds%imptorsn(j,i) < 1) exit
                V = V + bds%imptorsamp(j,i) * (1+cos(real(bds%imptorsn(j,i))*thet &
                                            - bds%imptorsphase(j,i)))
            end do
        end do

    end subroutine imptorsion_potential
    
    subroutine imptorsion_geomgrad(bds, grad)
        !! Compute torsion potential
        use mod_jacobian_mat, only: torsion_angle_jacobian

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        ! Bonded potential data structure
        real(rp), intent(inout) :: grad(3,bds%top%mm_atoms)
        !! improper torsion potential, result will be added to V
        real(rp) :: thet, g, J_a(3), J_b(3), J_c(3), J_d(3)
        integer(ip) :: i, j, ia, ib, ic, id
        logical :: sk_a, sk_b, sk_c, sk_d
        
        if(.not. bds%use_imptorsion) return
        
        do i=1, bds%nimptorsion
            ! Atoms that defines the dihedral angle
            ia = bds%imptorsionat(1,i)
            ib = bds%imptorsionat(2,i)
            ic = bds%imptorsionat(3,i)
            id = bds%imptorsionat(4,i) 
            
            if(bds%top%use_frozen) then
                sk_a = bds%top%frozen(ia)
                sk_b = bds%top%frozen(ib)
                sk_c = bds%top%frozen(ic)
                sk_d = bds%top%frozen(id)
                if(sk_a .and. sk_b .and. sk_c .and. sk_d) cycle
            else
                sk_a = .false.
                sk_b = .false.
                sk_c = .false.
                sk_d = .false.
            end if

            call torsion_angle_jacobian(bds%top%cmm(:,ia), &
                                        bds%top%cmm(:,ib), &
                                        bds%top%cmm(:,ic), &
                                        bds%top%cmm(:,id), &
                                        thet, J_a, J_b, J_c, J_d)
            
            do j=1, 3
                if(bds%imptorsn(j,i) < 1) exit
                g = -real(bds%imptorsn(j,i)) * sin(real(bds%imptorsn(j,i))* thet &
                                                   - bds%imptorsphase(j,i)) &
                                             * bds%imptorsamp(j,i)
                if(.not. sk_a) grad(:, ia) = grad(:, ia) + J_a * g
                if(.not. sk_b) grad(:, ib) = grad(:, ib) + J_b * g
                if(.not. sk_c) grad(:, ic) = grad(:, ic) + J_c * g
                if(.not. sk_d) grad(:, id) = grad(:, id) + J_d * g
            end do
        end do

    end subroutine imptorsion_geomgrad
    
    subroutine imptorsion_init(bds, n)
        !! Initialize improper torsion potential arrays

        use mod_memory, only: mallocate

        implicit none

        type(ommp_bonded_type), intent(inout) :: bds
        ! Bonded potential data structure
        integer(ip) :: n
        !! Number of improper torsion functions in the potential
        !! energy of the system
        
        if( n < 1 ) return
        bds%use_imptorsion = .true.

        call mallocate('imptorsion_init [imptorsionat]', 4, n, bds%imptorsionat)
        call mallocate('imptorsion_init [imptorsamp]', 3, n, bds%imptorsamp)
        call mallocate('imptorsion_init [imptorsphase]', 3, n, bds%imptorsphase)
        call mallocate('imptorsion_init [imptorsn]', 3, n, bds%imptorsn)

        bds%nimptorsion = n

    end subroutine imptorsion_init
    
    subroutine angtor_init(bds, n)
        !! Initialize angle-torsion coupling potential arrays

        use mod_memory, only: mallocate

        implicit none

        type(ommp_bonded_type), intent(inout) :: bds
        ! Bonded potential data structure
        integer(ip) :: n
        !! Number of angle torsion coupling functions in the potential
        !! energy of the system
        
        if( n < 1 ) return
        bds%use_angtor = .true.

        call mallocate('angtor_init [angtorat]', 4, n, bds%angtorat)
        call mallocate('angtor_init [angtork]', 6, n, bds%angtork)
        call mallocate('angtor_init [angtor_t]', n, bds%angtor_t)
        call mallocate('angtor_init [angtor_a]', 2, n, bds%angtor_a)

        bds%nangtor = n

    end subroutine angtor_init
    
    subroutine strtor_init(bds, n)
        
        use mod_memory, only: mallocate

        implicit none

        type(ommp_bonded_type), intent(inout) :: bds
        ! Bonded potential data structure
        integer(ip) :: n
        
        if( n < 1 ) return
        bds%use_strtor = .true.

        call mallocate('strtor_init [strtorat]', 4, n, bds%strtorat)
        call mallocate('strtor_init [strtork]', 9, n, bds%strtork)
        call mallocate('strtor_init [strtor_t]', n, bds%strtor_t)
        call mallocate('strtor_init [strtor_a]', 3, n, bds%strtor_b)

        bds%nstrtor = n

    end subroutine strtor_init
    
    subroutine angtor_potential(bds, V)

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        ! Bonded potential data structure
        real(rp), intent(inout) :: V
        real(rp) :: thet, costhet, dihef(3), delta_a(2), vat, l1, l2, &
                    dr1(3), dr2(3), angle1, angle2
        integer(ip) :: i, j, k, ia1, ia2
        
        if(.not. bds%use_angtor) return

        do i=1, bds%nangtor
            ! Atoms that defines the dihedral angle
            costhet = cos_torsion(bds%top, bds%angtorat(:,i))
            thet = acos(costhet)
            do j=1, 3
                dihef(j) = 1.0 + cos(j*thet+bds%torsphase(j,bds%angtor_t(i)))
            end do

            ia1 = bds%angtor_a(1,i)
            ia2 = bds%angtor_a(2,i)
            
            dr1 = bds%top%cmm(:, bds%angleat(1,ia1)) - &
                  bds%top%cmm(:, bds%angleat(2,ia1))
            dr2 = bds%top%cmm(:, bds%angleat(3,ia1)) - &
                  bds%top%cmm(:, bds%angleat(2,ia1))
            l1 = norm2(dr1)
            l2 = norm2(dr2)
            angle1 = acos(dot_product(dr1, dr2)/(l1*l2))

            dr1 = bds%top%cmm(:, bds%angleat(1,ia2)) - &
                  bds%top%cmm(:, bds%angleat(2,ia2))
            dr2 = bds%top%cmm(:, bds%angleat(3,ia2)) - &
                  bds%top%cmm(:, bds%angleat(2,ia2))
            l1 = norm2(dr1)
            l2 = norm2(dr2)
            angle2 = acos(dot_product(dr1, dr2)/(l1*l2))
           
            delta_a(1) = angle1 - bds%eqangle(bds%angtor_a(1,i))
            delta_a(2) = angle2 - bds%eqangle(bds%angtor_a(2,i))

            do j=1,2
                vat = 0.0
                do k=1, 3
                    vat = vat + bds%angtork((j-1)*3+k,i) * dihef(k)
                end do
                V = V + vat * delta_a(j)
            end do
        end do

    end subroutine angtor_potential
    
    subroutine angtor_geomgrad(bds, grad)
        use mod_jacobian_mat, only: simple_angle_jacobian, torsion_angle_jacobian

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        ! Bonded potential data structure
        real(rp), intent(inout) :: grad(3,bds%top%mm_atoms)
        !! improper torsion potential, result will be added to V
        real(rp) :: thet, gt(3), dihef(3), da1, da2, angle1, angle2, &
                    Jt_a(3), Jt_b(3), Jt_c(3), Jt_d(3), &
                    Ja1_a(3), Ja1_b(3), Ja1_c(3), &
                    Ja2_a(3), Ja2_b(3), Ja2_c(3)

        integer(ip) :: i, j, k, ia1, ia2, &
                       it_a, it_b, it_c, it_d, &
                       ia1_a, ia1_b, ia1_c, &
                       ia2_a, ia2_b, ia2_c
        logical :: sk_ta, sk_tb, sk_tc, sk_td, &
                   sk_1a, sk_1b, sk_1c, &
                   sk_2a, sk_2b, sk_2c
        
        if(.not. bds%use_angtor) return

        do i=1, bds%nangtor
            ! Atoms that defines the dihedral angle
            it_a = bds%angtorat(1,i)
            it_b = bds%angtorat(2,i)
            it_c = bds%angtorat(3,i)
            it_d = bds%angtorat(4,i) 

            ia1 = bds%angtor_a(1,i)
            ia1_a = bds%angleat(1,ia1)
            ia1_b = bds%angleat(2,ia1)
            ia1_c = bds%angleat(3,ia1)
            
            ia2 = bds%angtor_a(2,i)
            ia2_a = bds%angleat(1,ia2)
            ia2_b = bds%angleat(2,ia2)
            ia2_c = bds%angleat(3,ia2)
            
            if(bds%top%use_frozen) then
                sk_ta = bds%top%frozen(it_a)
                sk_tb = bds%top%frozen(it_b)
                sk_tc = bds%top%frozen(it_c)
                sk_td = bds%top%frozen(it_d)
                
                sk_1a = bds%top%frozen(ia1_a)
                sk_1b = bds%top%frozen(ia1_b)
                sk_1c = bds%top%frozen(ia1_c)
                
                sk_2a = bds%top%frozen(ia2_a)
                sk_2b = bds%top%frozen(ia2_b)
                sk_2c = bds%top%frozen(ia2_c)
                if(sk_ta .and. sk_tb .and. sk_tc .and. sk_td .and. &
                   sk_1a .and. sk_1b .and. sk_1c .and. &
                   sk_2a .and. sk_2b .and. sk_2c) cycle
            else
                sk_ta = .false.
                sk_tb = .false.
                sk_tc = .false.
                sk_td = .false.
                sk_1a = .false.
                sk_1b = .false.
                sk_1c = .false.
                sk_2a = .false.
                sk_2b = .false.
                sk_2c = .false.
            end if

            call torsion_angle_jacobian(bds%top%cmm(:,it_a), &
                                        bds%top%cmm(:,it_b), &
                                        bds%top%cmm(:,it_c), &
                                        bds%top%cmm(:,it_d), &
                                        thet, Jt_a, Jt_b, Jt_c, Jt_d)
            do j=1, 3
                gt(j) = -real(j) * sin(j*thet+bds%torsphase(j,bds%angtor_t(i)))
                dihef(j) = 1.0 + cos(j*thet+bds%torsphase(j,bds%angtor_t(i)))
            end do

            call simple_angle_jacobian(bds%top%cmm(:,ia1_a), &
                                       bds%top%cmm(:,ia1_b), &
                                       bds%top%cmm(:,ia1_c), &
                                       angle1, Ja1_a, Ja1_b, Ja1_c)
            
            call simple_angle_jacobian(bds%top%cmm(:,ia2_a), &
                                       bds%top%cmm(:,ia2_b), &
                                       bds%top%cmm(:,ia2_c), &
                                       angle2, Ja2_a, Ja2_b, Ja2_c)

            da1 = angle1 - bds%eqangle(ia1)
            da2 = angle2 - bds%eqangle(ia2)

            do k=1, 3
                if(.not. sk_ta) grad(:,it_a) = grad(:,it_a) + bds%angtork(k,i) * da1 * gt(k) * Jt_a
                if(.not. sk_tb) grad(:,it_b) = grad(:,it_b) + bds%angtork(k,i) * da1 * gt(k) * Jt_b
                if(.not. sk_tc) grad(:,it_c) = grad(:,it_c) + bds%angtork(k,i) * da1 * gt(k) * Jt_c
                if(.not. sk_td) grad(:,it_d) = grad(:,it_d) + bds%angtork(k,i) * da1 * gt(k) * Jt_d
                
                if(.not. sk_1a) grad(:,ia1_a) = grad(:,ia1_a) + bds%angtork(k,i) * dihef(k) * Ja1_a
                if(.not. sk_1b) grad(:,ia1_b) = grad(:,ia1_b) + bds%angtork(k,i) * dihef(k) * Ja1_b
                if(.not. sk_1c) grad(:,ia1_c) = grad(:,ia1_c) + bds%angtork(k,i) * dihef(k) * Ja1_c
                
                if(.not. sk_ta) grad(:,it_a) = grad(:,it_a) + bds%angtork(3+k,i) * da2 * gt(k) * Jt_a
                if(.not. sk_tb) grad(:,it_b) = grad(:,it_b) + bds%angtork(3+k,i) * da2 * gt(k) * Jt_b
                if(.not. sk_tc) grad(:,it_c) = grad(:,it_c) + bds%angtork(3+k,i) * da2 * gt(k) * Jt_c
                if(.not. sk_td) grad(:,it_d) = grad(:,it_d) + bds%angtork(3+k,i) * da2 * gt(k) * Jt_d
                
                if(.not. sk_2a) grad(:,ia2_a) = grad(:,ia2_a) + bds%angtork(3+k,i) * dihef(k) * Ja2_a
                if(.not. sk_2b) grad(:,ia2_b) = grad(:,ia2_b) + bds%angtork(3+k,i) * dihef(k) * Ja2_b
                if(.not. sk_2c) grad(:,ia2_c) = grad(:,ia2_c) + bds%angtork(3+k,i) * dihef(k) * Ja2_c
            end do
        end do
    end subroutine angtor_geomgrad
    
    subroutine strtor_potential(bds, V)
        use mod_constants

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        ! Bonded potential data structure
        real(rp), intent(inout) :: V
        real(rp) :: thet, costhet, dihef(3), dr(3), r(3), vst
        integer(ip) :: i, j, k, ib1, ib2, ib3
        
        if(.not. bds%use_strtor) return

        do i=1, bds%nstrtor
            ! Atoms that defines the dihedral angle
            costhet = cos_torsion(bds%top, bds%strtorat(:,i))
            thet = acos(costhet)
            do j=1, 3
                dihef(j) = 1.0 + cos(j*thet+bds%torsphase(j,bds%strtor_t(i)))
            end do

            ib1 = bds%strtor_b(1,i) 
            ib2 = bds%strtor_b(2,i)
            ib3 = bds%strtor_b(3,i)
            r(1) = norm2(bds%top%cmm(:, bds%bondat(1,ib1)) - &
                         bds%top%cmm(:, bds%bondat(2,ib1)))
            r(2) = norm2(bds%top%cmm(:, bds%bondat(1,ib2)) - &
                         bds%top%cmm(:, bds%bondat(2,ib2)))
            r(3) = norm2(bds%top%cmm(:, bds%bondat(1,ib3)) - &
                         bds%top%cmm(:, bds%bondat(2,ib3)))
            dr(1) = r(1) - bds%l0bond(ib1)  
            dr(2) = r(2) - bds%l0bond(ib2)  
            dr(3) = r(3) - bds%l0bond(ib3)  
            
            do j=1,3
                vst = 0.0
                do k=1, 3
                    vst = vst + bds%strtork((j-1)*3+k,i) * dihef(k)
                end do
                V = V + vst * dr(j)
            end do
        end do

    end subroutine strtor_potential
    
    subroutine strtor_geomgrad(bds, grad)
        use mod_jacobian_mat, only: Rij_jacobian, torsion_angle_jacobian

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        ! Bonded potential data structure
        real(rp), intent(inout) :: grad(3,bds%top%mm_atoms)
        !! improper torsion potential, result will be added to V
        
        real(rp) :: thet, gt(3), dihef(3), dr1, dr2,  dr3, r1, r2, r3, &
                    Jt_a(3), Jt_b(3), Jt_c(3), Jt_d(3), &
                    Jb1_a(3), Jb1_b(3), &
                    Jb2_a(3), Jb2_b(3), &
                    Jb3_a(3), Jb3_b(3)

        integer(ip) :: i, j, k, ib1, ib2, ib3, &
                       it_a, it_b, it_c, it_d, &
                       ib1_a, ib1_b, &
                       ib2_a, ib2_b, &
                       ib3_a, ib3_b
        logical :: sk_ta, sk_tb, sk_tc, sk_td, &
                   sk_1a, sk_1b, &
                   sk_2a, sk_2b, &
                   sk_3a, sk_3b

        
        if(.not. bds%use_strtor) return

        do i=1, bds%nstrtor
            ! Atoms that defines the dihedral angle
            it_a = bds%strtorat(1,i)
            it_b = bds%strtorat(2,i)
            it_c = bds%strtorat(3,i)
            it_d = bds%strtorat(4,i) 

            ib1 = bds%strtor_b(1,i)
            ib1_a = bds%bondat(1,ib1)
            ib1_b = bds%bondat(2,ib1)
            
            ib2 = bds%strtor_b(2,i)
            ib2_a = bds%bondat(1,ib2)
            ib2_b = bds%bondat(2,ib2)
            
            ib3 = bds%strtor_b(3,i)
            ib3_a = bds%bondat(1,ib3)
            ib3_b = bds%bondat(2,ib3)
            
            if(bds%top%use_frozen) then
                sk_ta = bds%top%frozen(it_a)
                sk_tb = bds%top%frozen(it_b)
                sk_tc = bds%top%frozen(it_c)
                sk_td = bds%top%frozen(it_d)
                
                sk_1a = bds%top%frozen(ib1_a)
                sk_1b = bds%top%frozen(ib1_b)
                
                sk_2a = bds%top%frozen(ib2_a)
                sk_2b = bds%top%frozen(ib2_b)

                sk_3a = bds%top%frozen(ib3_a)
                sk_3b = bds%top%frozen(ib3_b)

                if(sk_ta .and. sk_tb .and. sk_tc .and. sk_td .and. &
                   sk_1a .and. sk_1b .and. &
                   sk_2a .and. sk_2b .and. &
                   sk_3a .and. sk_3b) cycle
            else
                sk_ta = .false.
                sk_tb = .false.
                sk_tc = .false.
                sk_td = .false.
                sk_1a = .false.
                sk_1b = .false.
                sk_2a = .false.
                sk_2b = .false.
                sk_3a = .false.
                sk_3b = .false.
            end if

            call torsion_angle_jacobian(bds%top%cmm(:,it_a), &
                                        bds%top%cmm(:,it_b), &
                                        bds%top%cmm(:,it_c), &
                                        bds%top%cmm(:,it_d), &
                                        thet, Jt_a, Jt_b, Jt_c, Jt_d)
            do j=1, 3
                gt(j) = -real(j) * sin(j*thet+bds%torsphase(j,bds%angtor_t(i)))
                dihef(j) = 1.0 + cos(j*thet+bds%torsphase(j,bds%angtor_t(i)))
            end do

            call Rij_jacobian(bds%top%cmm(:,ib1_a), &
                              bds%top%cmm(:,ib1_b), &
                              r1, Jb1_a, Jb1_b) 
            dr1 = r1 - bds%l0bond(ib1)  
            
            call Rij_jacobian(bds%top%cmm(:,ib2_a), &
                              bds%top%cmm(:,ib2_b), &
                              r2, Jb2_a, Jb2_b) 
            dr2 = r2 - bds%l0bond(ib2)  
            
            call Rij_jacobian(bds%top%cmm(:,ib3_a), &
                              bds%top%cmm(:,ib3_b), &
                              r3, Jb3_a, Jb3_b) 
            dr3 = r3 - bds%l0bond(ib3)  
            
            do k=1, 3
                if(.not. sk_ta) grad(:,it_a) = grad(:,it_a) + bds%strtork(k,i) * dr1 * gt(k) * Jt_a
                if(.not. sk_tb) grad(:,it_b) = grad(:,it_b) + bds%strtork(k,i) * dr1 * gt(k) * Jt_b
                if(.not. sk_tc) grad(:,it_c) = grad(:,it_c) + bds%strtork(k,i) * dr1 * gt(k) * Jt_c
                if(.not. sk_td) grad(:,it_d) = grad(:,it_d) + bds%strtork(k,i) * dr1 * gt(k) * Jt_d
                
                if(.not. sk_1a) grad(:,ib1_a) = grad(:,ib1_a) + bds%strtork(k,i) * dihef(k) * Jb1_a
                if(.not. sk_1b) grad(:,ib1_b) = grad(:,ib1_b) + bds%strtork(k,i) * dihef(k) * Jb1_b
                
                if(.not. sk_ta) grad(:,it_a) = grad(:,it_a) + bds%strtork(3+k,i) * dr2 * gt(k) * Jt_a
                if(.not. sk_tb) grad(:,it_b) = grad(:,it_b) + bds%strtork(3+k,i) * dr2 * gt(k) * Jt_b
                if(.not. sk_tc) grad(:,it_c) = grad(:,it_c) + bds%strtork(3+k,i) * dr2 * gt(k) * Jt_c
                if(.not. sk_td) grad(:,it_d) = grad(:,it_d) + bds%strtork(3+k,i) * dr2 * gt(k) * Jt_d
                
                if(.not. sk_2a) grad(:,ib2_a) = grad(:,ib2_a) + bds%strtork(3+k,i) * dihef(k) * Jb2_a
                if(.not. sk_2b) grad(:,ib2_b) = grad(:,ib2_b) + bds%strtork(3+k,i) * dihef(k) * Jb2_b

                if(.not. sk_ta) grad(:,it_a) = grad(:,it_a) + bds%strtork(6+k,i) * dr3 * gt(k) * Jt_a
                if(.not. sk_tb) grad(:,it_b) = grad(:,it_b) + bds%strtork(6+k,i) * dr3 * gt(k) * Jt_b
                if(.not. sk_tc) grad(:,it_c) = grad(:,it_c) + bds%strtork(6+k,i) * dr3 * gt(k) * Jt_c
                if(.not. sk_td) grad(:,it_d) = grad(:,it_d) + bds%strtork(6+k,i) * dr3 * gt(k) * Jt_d
                
                if(.not. sk_3a) grad(:,ib3_a) = grad(:,ib3_a) + bds%strtork(6+k,i) * dihef(k) * Jb3_a
                if(.not. sk_3b) grad(:,ib3_b) = grad(:,ib3_b) + bds%strtork(6+k,i) * dihef(k) * Jb3_b
            end do
        end do
    end subroutine strtor_geomgrad
    
    subroutine tortor_init(bds, n)
        !! Initialize torsion-torsion correction potential arrays

        use mod_memory, only: mallocate

        implicit none

        type(ommp_bonded_type), intent(inout) :: bds
        ! Bonded potential data structure
        integer(ip) :: n
        !! Number of torsion-torsion 'map' functions in the potential
        !! energy of the system

        if( n < 1 ) return
        bds%use_tortor = .true.
        
        call mallocate('torsion_init [tortorprm]', n, bds%tortorprm )
        call mallocate('torsion_init [tortorat]', 5, n, bds%tortorat)

        bds%ntortor = n

    end subroutine tortor_init

    subroutine tortor_newmap(bds, d1, d2, ang1, ang2, v)
        !! Store in module memory the data describing a new torsion-torsion 
        !! map
        use mod_memory, only: mallocate, mfree
        use mod_utils, only: cyclic_spline

        implicit none

        type(ommp_bonded_type), intent(inout) :: bds
        ! Bonded potential data structure
        integer(ip), intent(in) :: d1, d2
        !! Dimensions of the new map to be saved
        real(rp), intent(in) :: ang1(:)
        !! Value of torsion1 for the new map 
        real(rp), intent(in) :: ang2(:)
        !! Value of torsion2 for the new map 
        real(rp), intent(in) :: v(:)
        !! Value of potential for the new map 

        integer :: i, j, ii
        real(rp), allocatable, dimension(:) :: a, b, c, d, dx, dy, dxy, tmpx, tmpy
        
        real(rp), allocatable :: rtmp(:)
        integer(ip), allocatable :: itmp(:,:)
        integer(ip) :: n_data, n_map

        if(allocated(bds%ttmap_ang1)) then
            ! Reallocate the arrays to make space for the new data
            n_data = size(bds%ttmap_ang1)
            call mallocate('torstors_newmap [rtmp]', n_data, rtmp)
            
            rtmp = bds%ttmap_ang1
            call mfree('torstors_newmap [ttmap_ang1]', bds%ttmap_ang1)
            call mallocate('torstors_newmap [ttmap_ang1]', &
                           n_data+d1*d2,  bds%ttmap_ang1)
            bds%ttmap_ang1(:n_data) = rtmp
            
            rtmp = bds%ttmap_ang2
            call mfree('torstors_newmap [ttmap_ang2]', bds%ttmap_ang2)
            call mallocate('torstors_newmap [ttmap_ang2]', &
                           n_data+d1*d2,  bds%ttmap_ang2)
            bds%ttmap_ang2(:n_data) = rtmp
            
            
            call mfree('torstors_newmap [rtmp]', rtmp)
            n_data = size(bds%ttmap_v)
            call mallocate('torstors_newmap [rtmp]', n_data, rtmp)
            
            rtmp = bds%ttmap_v
            call mfree('torstors_newmap [ttmap_v]', bds%ttmap_v)
            call mallocate('torstors_newmap [ttmap_v]', &
                           n_data+d1*d2,  bds%ttmap_v)
            bds%ttmap_v(:n_data) = rtmp
            
            rtmp = bds%ttmap_vx
            call mfree('torstors_newmap [ttmap_vx]', bds%ttmap_vx)
            call mallocate('torstors_newmap [ttmap_vx]', &
                           n_data+d1*d2,  bds%ttmap_vx)
            bds%ttmap_vx(:n_data) = rtmp

            rtmp = bds%ttmap_vy
            call mfree('torstors_newmap [ttmap_vy]', bds%ttmap_vy)
            call mallocate('torstors_newmap [ttmap_vy]', &
                           n_data+d1*d2,  bds%ttmap_vy)
            bds%ttmap_vy(:n_data) = rtmp

            rtmp = bds%ttmap_vxy
            call mfree('torstors_newmap [ttmap_vxy]', bds%ttmap_vxy)
            call mallocate('torstors_newmap [ttmap_vxy]', &
                           n_data+d1*d2,  bds%ttmap_vxy)
            bds%ttmap_vxy(:n_data) = rtmp
            call mfree('torstors_newmap [rtmp]', rtmp)

            n_map = size(bds%ttmap_shape, 2)
            call mallocate('torstors_newmap [itmp]', 2, n_map, itmp)
            itmp = bds%ttmap_shape
            call mfree('torstors_newmap [ttmap_shape]', bds%ttmap_shape)
            call mallocate('torstors_newmap [ttmap_shape]', &
                           2, n_map+1, bds%ttmap_shape)
            bds%ttmap_shape(:,:n_map) = itmp

            call mfree('torstors_newmap [itmp]', itmp)
        else 
            ! First allocation, n_data and n_map are just set for consistency
            n_data = 0
            n_map = 0
            call mallocate('torstors_newmap [ttmap_ang1]', d1*d2,  bds%ttmap_ang1)
            call mallocate('torstors_newmap [ttmap_ang2]', d1*d2,  bds%ttmap_ang2)
            call mallocate('torstors_newmap [ttmap_v]', d1*d2,  bds%ttmap_v)
            call mallocate('torstors_newmap [ttmap_vx]', d1*d2,  bds%ttmap_vx)
            call mallocate('torstors_newmap [ttmap_vy]', d1*d2,  bds%ttmap_vy)
            call mallocate('torstors_newmap [ttmap_vxy]', d1*d2,  bds%ttmap_vxy)
            call mallocate('torstors_newmap [ttmap_shape]', 2, 1, bds%ttmap_shape)
        end if

        call mallocate('tortor_newmap [a]', max(d1,d2), a)
        call mallocate('tortor_newmap [b]', max(d1,d2), b)
        call mallocate('tortor_newmap [c]', max(d1,d2), c)
        call mallocate('tortor_newmap [d]', max(d1,d2), d)
        call mallocate('tortor_newmap [dx]', d1*d2, dx)
        call mallocate('tortor_newmap [dy]', d1*d2, dy)
        call mallocate('tortor_newmap [dxy]', d1*d2, dxy)

        ! This part of the code computes df/dx, df/dy and d^2f/dxdy on the grid.
        ! Since we are basically interpolating on a sphere, we extract the 
        ! coordinate on a meridian, we interpolate it with a cubic spline, and
        ! finally we compute the derivative of this curve at the grid intersection
        ! The same is done in the second direction.
        ! To compute the mixed derivative we apply the same procedure but using
        ! the derivative data (basically we apply the procedure used to compute
        ! df/dx but using  df/dy data instead of actual f values.
        do i=1, d2
            call cyclic_spline(d1, ang1((i-1)*d1+1:i*d1), v((i-1)*d1+1:i*d1), &
                               a(1:d1), b(1:d1), c(1:d1), d(1:d1))
            dx((i-1)*d1+1:i*d1) = b(1:d1)
        end do
        
        ! df/dy since in this direction data are not contiguous, wa allocate 
        ! temporary arrays
        call mallocate('tortor_newmap [tmpx]', d2, tmpx)
        call mallocate('tortor_newmap [tmpy]', d2, tmpy)
        do i=1, d1
            ii = 1
            do j=i, (d2-1)*d1+i, d2
                tmpx(ii) = ang2(j)
                tmpy(ii) = v(j)
                ii = ii + 1
            end do
            call cyclic_spline(d2, tmpx, tmpy, &
                               a(1:d2), b(1:d2), c(1:d2), d(1:d2))
            
            ii = 1
            do j=i, (d2-1)*d1+i, d2
                dy(j) = b(ii)
                ii = ii + 1
            end do
        end do
        
        ! d^2f/dxdy in this case we use df/dx procedure to exploit data contiguity.
        do i=1, d2
            call cyclic_spline(d1, ang1((i-1)*d1+1:i*d1), dy((i-1)*d1+1:i*d1), &
                               a(1:d1), b(1:d1), c(1:d1), d(1:d1))
            dxy((i-1)*d1+1:i*d1) = b(1:d1)
        end do
        call mfree('tortor_newmap [tmpx]', tmpx)
        call mfree('tortor_newmap [tmpy]', tmpy)

        bds%ttmap_ang1(n_data+1:) = ang1
        bds%ttmap_ang2(n_data+1:) = ang2
        bds%ttmap_shape(1,n_map+1) = d1
        bds%ttmap_shape(2,n_map+1) = d2
        bds%ttmap_v(n_data+1:) = v
        bds%ttmap_vx(n_data+1:) = dx
        bds%ttmap_vy(n_data+1:) = dy
        bds%ttmap_vxy(n_data+1:) = dxy
        
        call mfree('tortor_newmap [a]', a)
        call mfree('tortor_newmap [b]', b)
        call mfree('tortor_newmap [c]', c)
        call mfree('tortor_newmap [d]', d)
        call mfree('tortor_newmap [dx]', dx)
        call mfree('tortor_newmap [dy]', dy)
        call mfree('tortor_newmap [dxy]', dxy)

    end subroutine tortor_newmap

    subroutine tortor_potential(bds, V)
        !! Compute torsion potential

        use mod_utils, only: compute_bicubic_interp

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        ! Bonded potential data structure
        real(rp), intent(inout) :: V
        !! torsion potential, result will be added to V
        real(rp) :: thetx, thety, vtt, dvttdx, dvttdy

        integer(ip) :: i, j, iprm, ibeg, iend

        if(.not. bds%use_tortor) return
        
        do i=1, bds%ntortor
            ! Atoms that defines the two angles
            iprm = bds%tortorprm(i)
            ibeg = 1
            do j=1, iprm-1
                ibeg = ibeg + bds%ttmap_shape(1,j)*bds%ttmap_shape(2,j)
            end do
            iend = ibeg + bds%ttmap_shape(1,iprm)*bds%ttmap_shape(2,iprm) - 1
           
            thetx = ang_torsion(bds%top, bds%tortorat(1:4,i))
            thety = ang_torsion(bds%top, bds%tortorat(2:5,i))
           
            call compute_bicubic_interp(thetx, thety, vtt, &
                                        dvttdx, dvttdy, &
                                        bds%ttmap_shape(1,iprm), &
                                        bds%ttmap_shape(2,iprm), &
                                        bds%ttmap_ang1(ibeg:iend), &
                                        bds%ttmap_ang2(ibeg:iend), &
                                        bds%ttmap_v(ibeg:iend), &
                                        bds%ttmap_vx(ibeg:iend), &
                                        bds%ttmap_vy(ibeg:iend), &
                                        bds%ttmap_vxy(ibeg:iend))

            V = V + vtt
        end do

    end subroutine tortor_potential
    
    subroutine tortor_geomgrad(bds, grad)
        !! Compute torsion potential

        use mod_utils, only: compute_bicubic_interp
        use mod_jacobian_mat, only: torsion_angle_jacobian

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        ! Bonded potential data structure
        real(rp), intent(inout) :: grad(3,bds%top%mm_atoms)
        !! improper torsion potential, result will be added to V
        real(rp) :: thetx, thety, vtt, dvttdx, dvttdy
        real(rp), dimension(3) :: J1_a, J1_b, J2_b, J1_c, &
                                  J2_c, J1_d, J2_d, J2_e

        integer(ip) :: i, j, iprm, ibeg, iend, ia, ib, ic, id, ie
        logical :: sk_a, sk_b, sk_c, sk_d, sk_e

        if(.not. bds%use_tortor) return
        
        do i=1, bds%ntortor
            ! Atoms that defines the two angles
            iprm = bds%tortorprm(i)
            ibeg = 1
            do j=1, iprm-1
                ibeg = ibeg + bds%ttmap_shape(1,j)*bds%ttmap_shape(2,j)
            end do
            iend = ibeg + bds%ttmap_shape(1,iprm)*bds%ttmap_shape(2,iprm) - 1
            
            ia = bds%tortorat(1,i)
            ib = bds%tortorat(2,i)
            ic = bds%tortorat(3,i)
            id = bds%tortorat(4,i)
            ie = bds%tortorat(5,i)
            
            if(bds%top%use_frozen) then
                sk_a = bds%top%frozen(ia)
                sk_b = bds%top%frozen(ib)
                sk_c = bds%top%frozen(ic)
                sk_d = bds%top%frozen(id)
                sk_e = bds%top%frozen(ie)
                if(sk_a .and. sk_b .and. sk_c .and. sk_d .and. sk_e) cycle
            else
                sk_a = .false.
                sk_b = .false.
                sk_c = .false.
                sk_d = .false.
                sk_e = .false.
            end if

            call torsion_angle_jacobian(bds%top%cmm(:,ia), &
                                        bds%top%cmm(:,ib), &
                                        bds%top%cmm(:,ic), &
                                        bds%top%cmm(:,id), &
                                        thetx, &
                                        J1_a, J1_b, J1_c, J1_d)
            thetx = ang_torsion(bds%top, bds%tortorat(1:4,i))
            
            call torsion_angle_jacobian(bds%top%cmm(:,ib), &
                                        bds%top%cmm(:,ic), &
                                        bds%top%cmm(:,id), &
                                        bds%top%cmm(:,ie), &
                                        thety, &
                                        J2_b, J2_c, J2_d, J2_e)
            thety = ang_torsion(bds%top, bds%tortorat(2:5,i))

            call compute_bicubic_interp(thetx, thety, vtt, &
                                        dvttdx, dvttdy, &
                                        bds%ttmap_shape(1,iprm), &
                                        bds%ttmap_shape(2,iprm), &
                                        bds%ttmap_ang1(ibeg:iend), &
                                        bds%ttmap_ang2(ibeg:iend), &
                                        bds%ttmap_v(ibeg:iend), &
                                        bds%ttmap_vx(ibeg:iend), &
                                        bds%ttmap_vy(ibeg:iend), &
                                        bds%ttmap_vxy(ibeg:iend))

            if(.not. sk_a) grad(:,ia) = grad(:,ia) + J1_a * dvttdx
            if(.not. sk_b) grad(:,ib) = grad(:,ib) + J1_b * dvttdx + J2_b * dvttdy
            if(.not. sk_c) grad(:,ic) = grad(:,ic) + J1_c * dvttdx + J2_c * dvttdy
            if(.not. sk_d) grad(:,id) = grad(:,id) + J1_d * dvttdx + J2_d * dvttdy
            if(.not. sk_e) grad(:,ie) = grad(:,ie) + J2_e * dvttdy
        end do

    end subroutine tortor_geomgrad

    pure function cos_torsion(top, idx)
        !! Compute the cosine of torsional angle between four atoms specified
        !! with indices idx
        
        implicit none

        type(ommp_topology_type), intent(in) :: top
        integer(ip), intent(in) :: idx(4)
        real(rp) :: cos_torsion

        real(rp), dimension(3) :: a, b, c, d, ab, cd, cb, t, u
            
        a = top%cmm(:,idx(1))
        b = top%cmm(:,idx(2))
        c = top%cmm(:,idx(3))
        d = top%cmm(:,idx(4))

        ab = b - a
        cd = d - c
        cb = b - c

        t(1) = ab(2)*cb(3) - ab(3)*cb(2)
        t(2) = ab(3)*cb(1) - ab(1)*cb(3)
        t(3) = ab(1)*cb(2) - ab(2)*cb(1)
        t = t / norm2(t)
            
        u(1) = cb(2)*cd(3) - cb(3)*cd(2)
        u(2) = cb(3)*cd(1) - cb(1)*cd(3)
        u(3) = cb(1)*cd(2) - cb(2)*cd(1)
        u = u / norm2(u)
            
        cos_torsion = dot_product(u,t)
        return 

    end function
    
    pure function ang_torsion(top, idx)
        !! Compute the torsional angle between four atoms specified
        !! with indices idx; results are in range [-pi;pi]
        
        implicit none

        type(ommp_topology_type), intent(in) :: top
        integer(ip), intent(in) :: idx(4)
        real(rp) :: cos_torsion, ang_torsion

        real(rp), dimension(3) :: a, b, c, d, ab, cd, cb, t, u
            
        a = top%cmm(:,idx(1))
        b = top%cmm(:,idx(2))
        c = top%cmm(:,idx(3))
        d = top%cmm(:,idx(4))

        ab = b - a
        cd = d - c
        cb = b - c

        t(1) = ab(2)*cb(3) - ab(3)*cb(2)
        t(2) = ab(3)*cb(1) - ab(1)*cb(3)
        t(3) = ab(1)*cb(2) - ab(2)*cb(1)
        t = t / norm2(t)
            
        u(1) = cb(2)*cd(3) - cb(3)*cd(2)
        u(2) = cb(3)*cd(1) - cb(1)*cd(3)
        u(3) = cb(1)*cd(2) - cb(2)*cd(1)
        u = u / norm2(u)
            
        cos_torsion = dot_product(u,t)
        ang_torsion = acos(cos_torsion)
        !if(dot_product(ab, u) > 0) ang_torsion = - ang_torsion
        ang_torsion = ang_torsion * sign(1.0_rp, -dot_product(ab,u))

    end function

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
        call imptorsion_terminate(bds)
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
    
    subroutine imptorsion_terminate(bds)
        use mod_memory, only: mfree

        implicit none

        type(ommp_bonded_type), intent(inout) :: bds
        ! Bonded potential data structure
        if( .not. bds%use_imptorsion ) return
        
        bds%use_imptorsion = .false.
        call mfree('imptorsion_terminate [imptorsionat]', bds%imptorsionat)
        call mfree('imptorsion_terminate [imptorsamp]', bds%imptorsamp)
        call mfree('imptorsion_terminate [imptorsphase]', bds%imptorsphase)
        call mfree('imptorsion_terminate [imptorsn]', bds%imptorsn)

    end subroutine imptorsion_terminate
    
    subroutine tortor_terminate(bds)
        use mod_memory, only: mfree

        implicit none

        type(ommp_bonded_type), intent(inout) :: bds
        ! Bonded potential data structure
        if( .not. bds%use_tortor ) return
        
        bds%use_tortor = .false.
        call mfree('tortor_terminate [tortorprm]', bds%tortorprm )
        call mfree('tortor_terminate [tortorat]', bds%tortorat)
        call mfree('tortor_terminate [ttmap_shape]', bds%ttmap_shape)
        call mfree('tortor_terminate [ttmap_ang1]', bds%ttmap_ang1)
        call mfree('tortor_terminate [ttmap_ang2]', bds%ttmap_ang2)
        call mfree('tortor_terminate [ttmap_v]', bds%ttmap_v)
        call mfree('tortor_terminate [ttmap_vx]', bds%ttmap_vx)
        call mfree('tortor_terminate [ttmap_vy]', bds%ttmap_vy)
        call mfree('tortor_terminate [ttmap_vxy]', bds%ttmap_vxy)

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
