module mod_bonded
    !! Module to handle the bonded part of the FF, it closely follows the 
    !! AMOEBA functional form.
    
    use mod_memory, only: ip, rp, lp
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
        logical(lp) :: use_bond = .false.
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
        logical(lp) :: use_angle = .false.
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
        logical(lp) :: use_strbnd = .false.
        !! Flag to enable calculation of stretching-bending coupling terms in 
        !! potential energy function
        
        ! Angle-Torsion coupling
        integer(ip) :: nangtor
        integer(ip), allocatable :: angtorat(:,:), angtor_t(:), angtor_a(:,:)
        real(rp), allocatable :: angtork(:,:)
        logical(lp) :: use_angtor = .false.
        
        ! Bond-Torsion coupling
        integer(ip) :: nstrtor
        integer(ip), allocatable :: strtorat(:,:), strtor_t(:), strtor_b(:,:)
        real(rp), allocatable :: strtork(:,:)
        logical(lp) :: use_strtor = .false.
        
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
        logical(lp) :: use_urey = .false.
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
        logical(lp) :: use_opb = .false.
        !! Flag to enable out-of-plane bending calculation
        
        ! Pi-torsion 
        integer(ip) :: npitors
        integer(ip), allocatable :: pitorsat(:,:)
        real(rp), allocatable :: kpitors(:)
        logical(lp) :: use_pitors = .false.

        ! Torsion
        integer(ip) :: ntorsion
        integer(ip), allocatable :: torsionat(:,:), torsn(:,:)
        real(rp), allocatable :: torsamp(:,:), torsphase(:,:)
        logical(lp) :: use_torsion = .false.
        
        ! Imporoper Torsion
        integer(ip) :: nimptorsion
        integer(ip), allocatable :: imptorsionat(:,:), imptorsn(:,:)
        real(rp), allocatable :: imptorsamp(:,:), imptorsphase(:,:)
        logical(lp) :: use_imptorsion = .false.

        ! Torsion-torsion coupling (cmap)
        integer(ip) :: ntortor
        integer(ip), allocatable :: tortorat(:,:), tortorprm(:), ttmap_shape(:,:)
        real(rp), allocatable :: ttmap_ang1(:), ttmap_ang2(:), ttmap_v(:), &
                                 ttmap_vx(:), ttmap_vy(:), ttmap_vxy(:)
        logical(lp) :: use_tortor = .false.
    end type ommp_bonded_type

    public :: ommp_bonded_type
    public :: bond_init, bond_potential, bond_geomgrad, bond_geomhess, bond_terminate
    public :: angle_init, angle_potential, angle_geomgrad, angle_geomhess, angle_terminate
    public :: urey_init, urey_potential, urey_geomgrad, urey_geomhess, urey_terminate
    public :: strbnd_init, strbnd_potential, strbnd_geomgrad, strbnd_geomhess, strbnd_terminate
    public :: opb_init, opb_potential, opb_geomgrad, opb_geomhess, opb_terminate
    public :: pitors_init, pitors_potential, pitors_geomgrad, pitors_geomhess, pitors_terminate
    public :: torsion_init, torsion_potential, torsion_geomgrad, torsion_geomhess, &
              torsion_terminate
    public :: imptorsion_init, imptorsion_potential, imptorsion_geomgrad, &
              imptorsion_geomhess, imptorsion_terminate
    public :: tortor_init, tortor_potential, tortor_geomgrad, tortor_geomhess, &
              tortor_terminate, tortor_newmap
    public :: strtor_init, strtor_potential, strtor_geomgrad, strtor_geomhess, strtor_terminate
    public :: angtor_init, angtor_potential, angtor_geomgrad, angtor_geomhess, angtor_terminate
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

        call mallocate('bond_init [bondat]', 2_ip, n, bds%bondat)
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
        logical(lp) :: use_cubic, use_quartic
        real(rp) :: dr(3), l, dl, dl2

        use_cubic = (abs(bds%bond_cubic) > eps_rp)
        use_quartic = (abs(bds%bond_quartic) > eps_rp)
        
        if(.not. bds%use_bond) return

        if(.not. use_cubic .and. .not. use_quartic) then
            ! This is just a regular harmonic potential
            !$omp parallel do default(shared) schedule(static) & 
            !$omp private(i,dr,l,dl) reduction(+:v) 
            do i=1, bds%nbond
                dr = bds%top%cmm(:,bds%bondat(1,i)) - &
                     bds%top%cmm(:,bds%bondat(2,i))
                l = sqrt(dot_product(dr, dr))
                dl = l - bds%l0bond(i)
                
                V = V + bds%kbond(i) * dl * dl
            end do
        else
            !$omp parallel do default(shared) schedule(static) & 
            !$omp private(i,dr,l,dl,dl2) reduction(+:v) 
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
        logical(lp) :: use_cubic, use_quartic
        logical :: sk_a, sk_b
        real(rp) :: ca(3), cb(3), J_a(3), J_b(3), l, dl, g

        use_cubic = (abs(bds%bond_cubic) > eps_rp)
        use_quartic = (abs(bds%bond_quartic) > eps_rp)

        if(.not. bds%use_bond) return

        if(.not. use_cubic .and. .not. use_quartic) then
            ! This is just a regular harmonic potential
            !$omp parallel do default(shared) schedule(dynamic) & 
            !$omp private(i,ia,ib,sk_a,sk_b,ca,cb,dl,l,g,J_a,J_b) 
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
                
                if(.not. sk_a) then
                    !$omp atomic update
                    grad(1,ia) = grad(1,ia) + J_a(1) * g
                    !$omp atomic update
                    grad(2,ia) = grad(2,ia) + J_a(2) * g
                    !$omp atomic update
                    grad(3,ia) = grad(3,ia) + J_a(3) * g
                end if

                if(.not. sk_b) then
                    !$omp atomic update
                    grad(1,ib) = grad(1,ib) + J_b(1) * g
                    !$omp atomic update
                    grad(2,ib) = grad(2,ib) + J_b(2) * g
                    !$omp atomic update
                    grad(3,ib) = grad(3,ib) + J_b(3) * g
                end if
            end do
        else
            !$omp parallel do default(shared) schedule(dynamic) & 
            !$omp private(i,ia,ib,sk_a,sk_b,ca,cb,dl,l,g,J_a,J_b) 
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

                g = 2 * bds%kbond(i) * dl * (1.0_rp + 3.0/2.0*bds%bond_cubic*dl &
                                             + 2.0*bds%bond_quartic*dl**2)
                
                if(.not. sk_a) then
                    !$omp atomic update
                    grad(1,ia) = grad(1,ia) + J_a(1) * g
                    !$omp atomic update
                    grad(2,ia) = grad(2,ia) + J_a(2) * g
                    !$omp atomic update
                    grad(3,ia) = grad(3,ia) + J_a(3) * g
                end if

                if(.not. sk_b) then
                    !$omp atomic update
                    grad(1,ib) = grad(1,ib) + J_b(1) * g
                    !$omp atomic update
                    grad(2,ib) = grad(2,ib) + J_b(2) * g
                    !$omp atomic update
                    grad(3,ib) = grad(3,ib) + J_b(3) * g
                end if
            end do
        end if

    end subroutine bond_geomgrad

    subroutine bond_geomhess(bds, hess)
        !! Compute the Hessian of the bond-stretching terms of the potential
        !! energy. With \(g=\partial U_i/\partial \Delta l_i\) as in
        !! bond_geomgrad and \(h=\partial^2 U_i/\partial \Delta l_i^2 =
        !! 2k_i(1+3k^{(3)}\Delta l_i+6k^{(4)}\Delta l_i^2)\), and \(M\) the
        !! symmetric matrix from Rij_hessian (equal for both diagonal
        !! blocks, with the off-diagonal block equal to \(-M\)):
        !! \[ H_{aa} = H_{bb} = -H_{ab} = -H_{ba} = h\, J_a J_a^\dagger + g M \]
        use mod_constants, only : eps_rp
        use mod_jacobian_mat, only: Rij_hessian

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        !! Bonded potential data structure
        real(rp), intent(inout) :: hess(3,3,bds%top%mm_atoms,bds%top%mm_atoms)
        !! Hessian of bond stretching terms of potential energy

        integer :: i, ia, ib, p
        logical(lp) :: use_cubic, use_quartic
        logical :: sk_a, sk_b
        real(rp) :: ca(3), cb(3), J_a(3), J_b(3), l, dl, g, h
        real(rp) :: M(3,3), H_ij(3,3), H_jj(3,3), Haa(3,3)

        use_cubic = (abs(bds%bond_cubic) > eps_rp)
        use_quartic = (abs(bds%bond_quartic) > eps_rp)

        if(.not. bds%use_bond) return

        !$omp parallel do default(shared) schedule(dynamic) &
        !$omp private(i,ia,ib,sk_a,sk_b,ca,cb,dl,l,g,h,J_a,J_b,M,H_ij,H_jj,Haa,p)
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

            call Rij_hessian(ca, cb, l, J_a, J_b, M, H_ij, H_jj)
            dl = l - bds%l0bond(i)

            if(.not. use_cubic .and. .not. use_quartic) then
                g = 2 * bds%kbond(i) * dl
                h = 2 * bds%kbond(i)
            else
                g = 2 * bds%kbond(i) * dl * (1.0_rp + 3.0/2.0*bds%bond_cubic*dl &
                                             + 2.0*bds%bond_quartic*dl**2)
                h = 2 * bds%kbond(i) * (1.0_rp + 3.0*bds%bond_cubic*dl &
                                        + 6.0*bds%bond_quartic*dl**2)
            end if

            do p=1,3
                Haa(p,:) = h*J_a(p)*J_a + g*M(p,:)
            end do

            if(.not. sk_a) then
                !$omp critical
                hess(:,:,ia,ia) = hess(:,:,ia,ia) + Haa
                !$omp end critical
            end if
            if(.not. sk_b) then
                !$omp critical
                hess(:,:,ib,ib) = hess(:,:,ib,ib) + Haa
                !$omp end critical
            end if
            if(.not. sk_a .and. .not. sk_b) then
                !$omp critical
                hess(:,:,ia,ib) = hess(:,:,ia,ib) - Haa
                hess(:,:,ib,ia) = hess(:,:,ib,ia) - Haa
                !$omp end critical
            end if
        end do

    end subroutine bond_geomhess

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

        call mallocate('angle_init [angleat]', 3_ip, n, bds%angleat)
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
        use mod_constants, only: eps_rp
        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        ! Bonded potential data structure
        real(rp), intent(inout) :: V
        !! Bond potential, result will be added to V
        
        integer(ip) :: i
        real(rp) :: l1, l2, dr1(3), dr2(3), thet, d_theta
        real(rp), dimension(3) :: v_dist, plv1, plv2, pln, a, b, c, prj_b, aux

        if(.not. bds%use_angle) return
        
        !$omp parallel do default(shared) schedule(static) reduction(+:V) &
        !$omp private(i,dr1,dr2,l1,l2,thet,d_theta,a,b,c,aux,plv1,plv2,pln,v_dist,prj_b) 
        do i=1, bds%nangle
            if(abs(bds%kangle(i)) < eps_rp) cycle
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
        use mod_constants, only: eps_rp

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
        
        !$omp parallel do default(shared) schedule(dynamic) &
        !$omp private(i,sk_a,sk_b,sk_c,sk_x,a,b,c,aux,thet,d_theta,g,Ja,Jb,Jc,Jx)
        do i=1, bds%nangle
            if(abs(bds%kangle(i)) < eps_rp) cycle
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

                if(.not. sk_a) then
                    !$omp atomic update
                    grad(1,bds%angleat(1,i)) = grad(1,bds%angleat(1,i)) + g * Ja(1)
                    !$omp atomic update
                    grad(2,bds%angleat(1,i)) = grad(2,bds%angleat(1,i)) + g * Ja(2)
                    !$omp atomic update
                    grad(3,bds%angleat(1,i)) = grad(3,bds%angleat(1,i)) + g * Ja(3)
                end if

                if(.not. sk_b) then
                    !$omp atomic update
                    grad(1,bds%angleat(2,i)) = grad(1,bds%angleat(2,i)) + g * Jb(1)
                    !$omp atomic update
                    grad(2,bds%angleat(2,i)) = grad(2,bds%angleat(2,i)) + g * Jb(2)
                    !$omp atomic update
                    grad(3,bds%angleat(2,i)) = grad(3,bds%angleat(2,i)) + g * Jb(3)
                end if

                if(.not. sk_c) then
                    !$omp atomic update
                    grad(1,bds%angleat(3,i)) = grad(1,bds%angleat(3,i)) + g * Jc(1)
                    !$omp atomic update
                    grad(2,bds%angleat(3,i)) = grad(2,bds%angleat(3,i)) + g * Jc(2)
                    !$omp atomic update
                    grad(3,bds%angleat(3,i)) = grad(3,bds%angleat(3,i)) + g * Jc(3)
                end if
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
                if(.not. sk_a) then
                    !$omp atomic update
                    grad(1,bds%angleat(1,i)) = grad(1,bds%angleat(1,i)) + g * Ja(1)
                    !$omp atomic update
                    grad(2,bds%angleat(1,i)) = grad(2,bds%angleat(1,i)) + g * Ja(2)
                    !$omp atomic update
                    grad(3,bds%angleat(1,i)) = grad(3,bds%angleat(1,i)) + g * Ja(3)
                end if

                if(.not. sk_b) then
                    !$omp atomic update
                    grad(1,bds%angleat(2,i)) = grad(1,bds%angleat(2,i)) + g * Jb(1)
                    !$omp atomic update
                    grad(2,bds%angleat(2,i)) = grad(2,bds%angleat(2,i)) + g * Jb(2)
                    !$omp atomic update
                    grad(3,bds%angleat(2,i)) = grad(3,bds%angleat(2,i)) + g * Jb(3)
                end if

                if(.not. sk_c) then
                    !$omp atomic update
                    grad(1,bds%angleat(3,i)) = grad(1,bds%angleat(3,i)) + g * Jc(1)
                    !$omp atomic update
                    grad(2,bds%angleat(3,i)) = grad(2,bds%angleat(3,i)) + g * Jc(2)
                    !$omp atomic update
                    grad(3,bds%angleat(3,i)) = grad(3,bds%angleat(3,i)) + g * Jc(3)
                end if

                if(.not. sk_x) then
                    !$omp atomic update
                    grad(1,bds%angauxat(i)) = grad(1,bds%angauxat(i)) + g * Jx(1)
                    !$omp atomic update
                    grad(2,bds%angauxat(i)) = grad(2,bds%angauxat(i)) + g * Jx(2)
                    !$omp atomic update
                    grad(3,bds%angauxat(i)) = grad(3,bds%angauxat(i)) + g * Jx(3)
                end if
            end if
        end do
    end subroutine angle_geomgrad

    subroutine angle_geomhess(bds, hess)
        !! Compute the Hessian of the angle-bending terms of the potential
        !! energy. With \(g=\partial U_i/\partial \Delta\theta_i\) as in
        !! angle_geomgrad and
        !! \(h=\partial^2 U_i/\partial\Delta\theta_i^2 = 2k_i(1+3k^{(3)}
        !! \Delta\theta_i+6k^{(4)}\Delta\theta_i^2+10k^{(5)}\Delta\theta_i^3
        !! +15k^{(6)}\Delta\theta_i^4)\), for any two atoms X,Y of the term:
        !! \[ H_{XY} = h\, J_X J_Y^\dagger + g\, H_{XY}(\theta) \]
        !! where \(H_{XY}(\theta)\) is the corresponding block from
        !! simple_angle_hessian (OMMP_ANG_SIMPLE/H0/H1/H2) or
        !! inplane_angle_hessian (OMMP_ANG_INPLANE/H0/H1).
        use mod_jacobian_mat, only: simple_angle_hessian, inplane_angle_hessian
        use mod_constants, only: eps_rp

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        !! Bonded potential data structure
        real(rp), intent(inout) :: hess(3,3,bds%top%mm_atoms,bds%top%mm_atoms)
        !! Hessian of angle bending terms of potential energy

        real(rp) :: a(3), b(3), c(3), x(3), Ja(3), Jb(3), Jc(3), Jx(3), &
                    g, h, thet, d_theta
        real(rp), dimension(3,3) :: Haa, Hab, Hac, Hax, Hbb, Hbc, Hbx, Hcc, Hcx, Hxx
        integer(ip) :: i, ia, ib, ic, ix, p
        logical :: sk_a, sk_b, sk_c, sk_x

        if(.not. bds%use_angle) return

        !$omp parallel do default(shared) schedule(dynamic) &
        !$omp private(i,ia,ib,ic,ix,sk_a,sk_b,sk_c,sk_x,a,b,c,x,thet,d_theta,g,h,p) &
        !$omp private(Ja,Jb,Jc,Jx,Haa,Hab,Hac,Hax,Hbb,Hbc,Hbx,Hcc,Hcx,Hxx)
        do i=1, bds%nangle
            if(abs(bds%kangle(i)) < eps_rp) cycle
            if(bds%anglety(i) == OMMP_ANG_SIMPLE .or. &
               bds%anglety(i) == OMMP_ANG_H0 .or. &
               bds%anglety(i) == OMMP_ANG_H1 .or. &
               bds%anglety(i) == OMMP_ANG_H2) then
                ia = bds%angleat(1,i)
                ib = bds%angleat(2,i)
                ic = bds%angleat(3,i)

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

                a = bds%top%cmm(:,ia)
                b = bds%top%cmm(:,ib)
                c = bds%top%cmm(:,ic)
                call simple_angle_hessian(a, b, c, thet, Ja, Jb, Jc, &
                                          Haa, Hab, Hac, Hbb, Hbc, Hcc)
                d_theta = thet - bds%eqangle(i)

                g = bds%kangle(i) * d_theta * (2.0 &
                                               + 3.0 * bds%angle_cubic * d_theta &
                                               + 4.0 * bds%angle_quartic * d_theta**2 &
                                               + 5.0 * bds%angle_pentic * d_theta**3 &
                                               + 6.0 * bds%angle_sextic * d_theta**4)
                h = 2.0 * bds%kangle(i) * (1.0 &
                                           + 3.0 * bds%angle_cubic * d_theta &
                                           + 6.0 * bds%angle_quartic * d_theta**2 &
                                           + 10.0 * bds%angle_pentic * d_theta**3 &
                                           + 15.0 * bds%angle_sextic * d_theta**4)

                do p=1,3
                    Haa(p,:) = h*Ja(p)*Ja + g*Haa(p,:)
                    Hab(p,:) = h*Ja(p)*Jb + g*Hab(p,:)
                    Hac(p,:) = h*Ja(p)*Jc + g*Hac(p,:)
                    Hbb(p,:) = h*Jb(p)*Jb + g*Hbb(p,:)
                    Hbc(p,:) = h*Jb(p)*Jc + g*Hbc(p,:)
                    Hcc(p,:) = h*Jc(p)*Jc + g*Hcc(p,:)
                end do

                !$omp critical
                if(.not. sk_a) hess(:,:,ia,ia) = hess(:,:,ia,ia) + Haa
                if(.not. sk_b) hess(:,:,ib,ib) = hess(:,:,ib,ib) + Hbb
                if(.not. sk_c) hess(:,:,ic,ic) = hess(:,:,ic,ic) + Hcc
                if(.not. sk_a .and. .not. sk_b) then
                    hess(:,:,ia,ib) = hess(:,:,ia,ib) + Hab
                    hess(:,:,ib,ia) = hess(:,:,ib,ia) + transpose(Hab)
                end if
                if(.not. sk_a .and. .not. sk_c) then
                    hess(:,:,ia,ic) = hess(:,:,ia,ic) + Hac
                    hess(:,:,ic,ia) = hess(:,:,ic,ia) + transpose(Hac)
                end if
                if(.not. sk_b .and. .not. sk_c) then
                    hess(:,:,ib,ic) = hess(:,:,ib,ic) + Hbc
                    hess(:,:,ic,ib) = hess(:,:,ic,ib) + transpose(Hbc)
                end if
                !$omp end critical
            else if(bds%anglety(i) == OMMP_ANG_INPLANE .or. &
                    bds%anglety(i) == OMMP_ANG_INPLANE_H0 .or. &
                    bds%anglety(i) == OMMP_ANG_INPLANE_H1) then

                ia = bds%angleat(1,i)
                ib = bds%angleat(2,i)
                ic = bds%angleat(3,i)
                ix = bds%angauxat(i)

                if(bds%top%use_frozen) then
                    sk_a = bds%top%frozen(ia)
                    sk_b = bds%top%frozen(ib)
                    sk_c = bds%top%frozen(ic)
                    sk_x = bds%top%frozen(ix)
                    if(sk_a .and. sk_b .and. sk_c .and. sk_x) cycle
                else
                    sk_a = .false.
                    sk_b = .false.
                    sk_c = .false.
                    sk_x = .false.
                end if

                a = bds%top%cmm(:,ia)
                b = bds%top%cmm(:,ib)
                c = bds%top%cmm(:,ic)
                x = bds%top%cmm(:,ix)
                call inplane_angle_hessian(a, b, c, x, thet, Ja, Jb, Jc, Jx, &
                                           Haa, Hab, Hac, Hax, Hbb, Hbc, Hbx, &
                                           Hcc, Hcx, Hxx)
                d_theta = thet - bds%eqangle(i)

                g = bds%kangle(i) * d_theta * (2.0 &
                                               + 3.0 * bds%angle_cubic * d_theta &
                                               + 4.0 * bds%angle_quartic * d_theta**2 &
                                               + 5.0 * bds%angle_pentic * d_theta**3 &
                                               + 6.0 * bds%angle_sextic * d_theta**4)
                h = 2.0 * bds%kangle(i) * (1.0 &
                                           + 3.0 * bds%angle_cubic * d_theta &
                                           + 6.0 * bds%angle_quartic * d_theta**2 &
                                           + 10.0 * bds%angle_pentic * d_theta**3 &
                                           + 15.0 * bds%angle_sextic * d_theta**4)

                do p=1,3
                    Haa(p,:) = h*Ja(p)*Ja + g*Haa(p,:)
                    Hab(p,:) = h*Ja(p)*Jb + g*Hab(p,:)
                    Hac(p,:) = h*Ja(p)*Jc + g*Hac(p,:)
                    Hax(p,:) = h*Ja(p)*Jx + g*Hax(p,:)
                    Hbb(p,:) = h*Jb(p)*Jb + g*Hbb(p,:)
                    Hbc(p,:) = h*Jb(p)*Jc + g*Hbc(p,:)
                    Hbx(p,:) = h*Jb(p)*Jx + g*Hbx(p,:)
                    Hcc(p,:) = h*Jc(p)*Jc + g*Hcc(p,:)
                    Hcx(p,:) = h*Jc(p)*Jx + g*Hcx(p,:)
                    Hxx(p,:) = h*Jx(p)*Jx + g*Hxx(p,:)
                end do

                !$omp critical
                if(.not. sk_a) hess(:,:,ia,ia) = hess(:,:,ia,ia) + Haa
                if(.not. sk_b) hess(:,:,ib,ib) = hess(:,:,ib,ib) + Hbb
                if(.not. sk_c) hess(:,:,ic,ic) = hess(:,:,ic,ic) + Hcc
                if(.not. sk_x) hess(:,:,ix,ix) = hess(:,:,ix,ix) + Hxx
                if(.not. sk_a .and. .not. sk_b) then
                    hess(:,:,ia,ib) = hess(:,:,ia,ib) + Hab
                    hess(:,:,ib,ia) = hess(:,:,ib,ia) + transpose(Hab)
                end if
                if(.not. sk_a .and. .not. sk_c) then
                    hess(:,:,ia,ic) = hess(:,:,ia,ic) + Hac
                    hess(:,:,ic,ia) = hess(:,:,ic,ia) + transpose(Hac)
                end if
                if(.not. sk_a .and. .not. sk_x) then
                    hess(:,:,ia,ix) = hess(:,:,ia,ix) + Hax
                    hess(:,:,ix,ia) = hess(:,:,ix,ia) + transpose(Hax)
                end if
                if(.not. sk_b .and. .not. sk_c) then
                    hess(:,:,ib,ic) = hess(:,:,ib,ic) + Hbc
                    hess(:,:,ic,ib) = hess(:,:,ic,ib) + transpose(Hbc)
                end if
                if(.not. sk_b .and. .not. sk_x) then
                    hess(:,:,ib,ix) = hess(:,:,ib,ix) + Hbx
                    hess(:,:,ix,ib) = hess(:,:,ix,ib) + transpose(Hbx)
                end if
                if(.not. sk_c .and. .not. sk_x) then
                    hess(:,:,ic,ix) = hess(:,:,ic,ix) + Hcx
                    hess(:,:,ix,ic) = hess(:,:,ix,ic) + transpose(Hcx)
                end if
                !$omp end critical
            end if
        end do
    end subroutine angle_geomhess

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

        call mallocate('strbnd_init [strbndat]', 3_ip, n, bds%strbndat)
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

        !$omp parallel do default(shared) reduction(+:V) &
        !$omp private(i,dr1,l1,l2,d_l1,d_l2,dr2,thet,d_thet)
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

        !$omp parallel do default(shared) schedule(dynamic) &
        !$omp private(i,ia,ib,ic,sk_a,sk_b,sk_c,a,b,c,l1,l2,d_l1,d_l2,thet,d_thet) &
        !$omp private(J1_a,J1_b,J2_b,J2_c,J3_a,J3_b,J3_c,g1,g2,g3)
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

            if(.not. sk_a) then
                !$omp atomic update
                grad(1,ia) = grad(1,ia) + J1_a(1) * g1 + J3_a(1) * g3
                !$omp atomic update
                grad(2,ia) = grad(2,ia) + J1_a(2) * g1 + J3_a(2) * g3
                !$omp atomic update
                grad(3,ia) = grad(3,ia) + J1_a(3) * g1 + J3_a(3) * g3
            end if

            if(.not. sk_b) then
                !$omp atomic update
                grad(1,ib) = grad(1,ib) + J1_b(1) * g1 + J2_b(1) * g2 + J3_b(1) * g3
                !$omp atomic update
                grad(2,ib) = grad(2,ib) + J1_b(2) * g1 + J2_b(2) * g2 + J3_b(2) * g3
                !$omp atomic update
                grad(3,ib) = grad(3,ib) + J1_b(3) * g1 + J2_b(3) * g2 + J3_b(3) * g3
            end if

            if(.not. sk_c) then
                !$omp atomic update
                grad(1,ic) = grad(1,ic) + J2_c(1) * g2 + J3_c(1) * g3
                !$omp atomic update
                grad(2,ic) = grad(2,ic) + J2_c(2) * g2 + J3_c(2) * g3
                !$omp atomic update
                grad(3,ic) = grad(3,ic) + J2_c(3) * g2 + J3_c(3) * g3
            end if
        end do

    end subroutine strbnd_geomgrad

    subroutine strbnd_geomhess(bds, hess)
        !! Compute the Hessian of the stretch-bend cross term. Treating
        !! \(U=(k_1\Delta l_1+k_2\Delta l_2)\Delta\theta\) as a function of
        !! the three internal coordinates \(l_1=R_{AB}\), \(l_2=R_{BC}\),
        !! \(\theta=\angle(A,B,C)\), the only nonzero internal second
        !! derivatives are \(\partial^2U/\partial l_1\partial\theta=k_1\) and
        !! \(\partial^2U/\partial l_2\partial\theta=k_2\), so for any two
        !! atoms X,Y of the term:
        !! \[ H_{XY} = k_1\left(J^{l_1}_XJ^{\theta\dagger}_Y+J^\theta_XJ^{l_1\dagger}_Y\right)
        !!    + k_2\left(J^{l_2}_XJ^{\theta\dagger}_Y+J^\theta_XJ^{l_2\dagger}_Y\right)
        !!    + g_1H^{l_1}_{XY} + g_2H^{l_2}_{XY} + g_3H^\theta_{XY} \]
        !! with \(g_1,g_2,g_3\) as in strbnd_geomgrad, and
        !! \(J^{l_1},H^{l_1}\) (resp. \(J^{l_2},H^{l_2}\)) from Rij_hessian
        !! on the A-B (resp. B-C) bond (zero on blocks not touching A,B --
        !! resp. B,C), \(J^\theta,H^\theta\) from simple_angle_hessian(A,B,C).
        use mod_jacobian_mat, only: Rij_hessian, simple_angle_hessian

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        !! Bonded potential data structure
        real(rp), intent(inout) :: hess(3,3,bds%top%mm_atoms,bds%top%mm_atoms)
        !! Hessian of stretch-bend terms of potential energy

        integer(ip) :: i, ia, ib, ic, p
        real(rp) :: d_l1, d_l2, d_thet, l1, l2, thet, g1, g2, g3
        real(rp), dimension(3) :: a, b, c, J1_a, J1_b, J2_b, J2_c, J3_a, J3_b, J3_c
        real(rp), dimension(3,3) :: H1_aa, H1_ab, H1_bb, H2_bb, H2_bc, H2_cc, &
                                    H3_aa, H3_ab, H3_ac, H3_bb, H3_bc, H3_cc, &
                                    Haa, Hab, Hac, Hbb, Hbc, Hcc
        logical :: sk_a, sk_b, sk_c

        if(.not. bds%use_strbnd) return

        !$omp parallel do default(shared) schedule(dynamic) &
        !$omp private(i,ia,ib,ic,sk_a,sk_b,sk_c,a,b,c,l1,l2,d_l1,d_l2,thet,d_thet) &
        !$omp private(J1_a,J1_b,J2_b,J2_c,J3_a,J3_b,J3_c,g1,g2,g3,p) &
        !$omp private(H1_aa,H1_ab,H1_bb,H2_bb,H2_bc,H2_cc,H3_aa,H3_ab,H3_ac,H3_bb,H3_bc,H3_cc) &
        !$omp private(Haa,Hab,Hac,Hbb,Hbc,Hcc)
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

            call Rij_hessian(a, b, l1, J1_a, J1_b, H1_aa, H1_ab, H1_bb)
            call Rij_hessian(b, c, l2, J2_b, J2_c, H2_bb, H2_bc, H2_cc)
            call simple_angle_hessian(a, b, c, thet, J3_a, J3_b, J3_c, &
                                      H3_aa, H3_ab, H3_ac, H3_bb, H3_bc, H3_cc)

            d_l1 = l1 - bds%strbndl10(i)
            d_l2 = l2 - bds%strbndl20(i)
            d_thet = thet - bds%strbndthet0(i)

            g1 = bds%strbndk1(i) * d_thet
            g2 = bds%strbndk2(i) * d_thet
            g3 = bds%strbndk1(i) * d_l1 + bds%strbndk2(i) * d_l2

            do p=1,3
                Haa(p,:) = bds%strbndk1(i)*(J1_a(p)*J3_a + J3_a(p)*J1_a) &
                          + g1*H1_aa(p,:) + g3*H3_aa(p,:)
                Hab(p,:) = bds%strbndk1(i)*(J1_a(p)*J3_b + J3_a(p)*J1_b) &
                          + bds%strbndk2(i)*J3_a(p)*J2_b &
                          + g1*H1_ab(p,:) + g3*H3_ab(p,:)
                Hac(p,:) = bds%strbndk1(i)*J1_a(p)*J3_c + bds%strbndk2(i)*J3_a(p)*J2_c &
                          + g3*H3_ac(p,:)
                Hbb(p,:) = bds%strbndk1(i)*(J1_b(p)*J3_b + J3_b(p)*J1_b) &
                          + bds%strbndk2(i)*(J2_b(p)*J3_b + J3_b(p)*J2_b) &
                          + g1*H1_bb(p,:) + g2*H2_bb(p,:) + g3*H3_bb(p,:)
                Hbc(p,:) = bds%strbndk1(i)*J1_b(p)*J3_c &
                          + bds%strbndk2(i)*(J2_b(p)*J3_c + J3_b(p)*J2_c) &
                          + g2*H2_bc(p,:) + g3*H3_bc(p,:)
                Hcc(p,:) = bds%strbndk2(i)*(J2_c(p)*J3_c + J3_c(p)*J2_c) &
                          + g2*H2_cc(p,:) + g3*H3_cc(p,:)
            end do

            !$omp critical
            if(.not. sk_a) hess(:,:,ia,ia) = hess(:,:,ia,ia) + Haa
            if(.not. sk_b) hess(:,:,ib,ib) = hess(:,:,ib,ib) + Hbb
            if(.not. sk_c) hess(:,:,ic,ic) = hess(:,:,ic,ic) + Hcc
            if(.not. sk_a .and. .not. sk_b) then
                hess(:,:,ia,ib) = hess(:,:,ia,ib) + Hab
                hess(:,:,ib,ia) = hess(:,:,ib,ia) + transpose(Hab)
            end if
            if(.not. sk_a .and. .not. sk_c) then
                hess(:,:,ia,ic) = hess(:,:,ia,ic) + Hac
                hess(:,:,ic,ia) = hess(:,:,ic,ia) + transpose(Hac)
            end if
            if(.not. sk_b .and. .not. sk_c) then
                hess(:,:,ib,ic) = hess(:,:,ib,ic) + Hbc
                hess(:,:,ic,ib) = hess(:,:,ic,ib) + transpose(Hbc)
            end if
            !$omp end critical
        end do

    end subroutine strbnd_geomhess

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

        call mallocate('urey_init [ureya]', 2_ip, n, bds%ureyat)
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
        logical(lp) :: use_cubic, use_quartic
        real(rp) :: dr(3), l, dl, dl2
        
        if(.not. bds%use_urey) return

        use_cubic = (abs(bds%urey_cubic) > eps_rp)
        use_quartic = (abs(bds%urey_quartic) > eps_rp)

        if(.not. use_cubic .and. .not. use_quartic) then
            ! This is just a regular harmonic potential
            !$omp parallel do default(shared) reduction(+:V) &
            !$omp private(i,dr,l,dl)
            do i=1, bds%nurey
                dr = bds%top%cmm(:,bds%ureyat(1,i)) - &
                     bds%top%cmm(:,bds%ureyat(2,i))
                l = sqrt(dot_product(dr, dr))
                dl = l - bds%l0urey(i)
                V = V + bds%kurey(i) * dl * dl
            end do
        else
            !$omp parallel do default(shared) reduction(+:V) &
            !$omp private(i,dr,l,dl,dl2)
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
        logical(lp) :: use_cubic, use_quartic
        logical :: sk_a, sk_b
        real(rp) :: l, dl, J_a(3), J_b(3), g
        
        if(.not. bds%use_urey) return

        use_cubic = (abs(bds%urey_cubic) > eps_rp)
        use_quartic = (abs(bds%urey_quartic) > eps_rp)

        if(.not. use_cubic .and. .not. use_quartic) then
            ! This is just a regular harmonic potential
            !$omp parallel do default(shared)  &
            !$omp private(i,ia,ib,sk_a,sk_b,l,dl,g,J_a,J_b)
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

                if(.not. sk_a) then
                    !$omp atomic update
                    grad(1,ia) = grad(1,ia) + J_a(1) * g
                    !$omp atomic update
                    grad(2,ia) = grad(2,ia) + J_a(2) * g
                    !$omp atomic update
                    grad(3,ia) = grad(3,ia) + J_a(3) * g
                end if

                if(.not. sk_b) then
                    !$omp atomic update
                    grad(1,ib) = grad(1,ib) + J_b(1) * g
                    !$omp atomic update
                    grad(2,ib) = grad(2,ib) + J_b(2) * g
                    !$omp atomic update
                    grad(3,ib) = grad(3,ib) + J_b(3) * g
                end if
            end do
        else
            !$omp parallel do default(shared) &
            !$omp private(i,ia,ib,sk_a,sk_b,l,dl,g,J_a,J_b)
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

                if(.not. sk_a) then
                    !$omp atomic update
                    grad(1,ia) = grad(1,ia) + J_a(1) * g
                    !$omp atomic update
                    grad(2,ia) = grad(2,ia) + J_a(2) * g
                    !$omp atomic update
                    grad(3,ia) = grad(3,ia) + J_a(3) * g
                end if

                if(.not. sk_b) then
                    !$omp atomic update
                    grad(1,ib) = grad(1,ib) + J_b(1) * g
                    !$omp atomic update
                    grad(2,ib) = grad(2,ib) + J_b(2) * g
                    !$omp atomic update
                    grad(3,ib) = grad(3,ib) + J_b(3) * g
                end if
            end do
        end if
    end subroutine urey_geomgrad

    subroutine urey_geomhess(bds, hess)
        !! Compute the Hessian of the Urey-Bradley terms of the potential
        !! energy. Formally identical to bond_geomhess (see there), applied
        !! to the ureyat/kurey/l0urey/urey_cubic/urey_quartic parameters.
        use mod_constants, only : eps_rp
        use mod_jacobian_mat, only: Rij_hessian

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        !! Bonded potential data structure
        real(rp), intent(inout) :: hess(3,3,bds%top%mm_atoms,bds%top%mm_atoms)
        !! Hessian of Urey-Bradley terms of potential energy

        integer :: i, ia, ib, p
        logical(lp) :: use_cubic, use_quartic
        logical :: sk_a, sk_b
        real(rp) :: J_a(3), J_b(3), l, dl, g, h
        real(rp) :: M(3,3), H_ij(3,3), H_jj(3,3), Haa(3,3)

        if(.not. bds%use_urey) return

        use_cubic = (abs(bds%urey_cubic) > eps_rp)
        use_quartic = (abs(bds%urey_quartic) > eps_rp)

        !$omp parallel do default(shared) schedule(dynamic) &
        !$omp private(i,ia,ib,sk_a,sk_b,dl,l,g,h,J_a,J_b,M,H_ij,H_jj,Haa,p)
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

            call Rij_hessian(bds%top%cmm(:,ia), bds%top%cmm(:,ib), &
                             l, J_a, J_b, M, H_ij, H_jj)
            dl = l - bds%l0urey(i)

            if(.not. use_cubic .and. .not. use_quartic) then
                g = 2 * bds%kurey(i) * dl
                h = 2 * bds%kurey(i)
            else
                g = 2 * bds%kurey(i) * dl * (1.0_rp + 3.0/2.0*bds%urey_cubic*dl &
                                             + 2.0*bds%urey_quartic*dl**2)
                h = 2 * bds%kurey(i) * (1.0_rp + 3.0*bds%urey_cubic*dl &
                                        + 6.0*bds%urey_quartic*dl**2)
            end if

            do p=1,3
                Haa(p,:) = h*J_a(p)*J_a + g*M(p,:)
            end do

            if(.not. sk_a) then
                !$omp critical
                hess(:,:,ia,ia) = hess(:,:,ia,ia) + Haa
                !$omp end critical
            end if
            if(.not. sk_b) then
                !$omp critical
                hess(:,:,ib,ib) = hess(:,:,ib,ib) + Haa
                !$omp end critical
            end if
            if(.not. sk_a .and. .not. sk_b) then
                !$omp critical
                hess(:,:,ia,ib) = hess(:,:,ia,ib) - Haa
                hess(:,:,ib,ia) = hess(:,:,ib,ia) - Haa
                !$omp end critical
            end if
        end do

    end subroutine urey_geomhess

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

        call mallocate('opb_init [opbat]', 4_ip, n, bds%opbat)
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
        
        !$omp parallel do default(shared) reduction(+:V) &
        !$omp private(i,a,b,c,d,plv1,plv2,pln,lpln,vad,lvad,thet,thet2,thet3,thet4)
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

        !$omp parallel do default(shared) schedule(dynamic)&
        !$omp private(i,ia,ib,ic,id,sk_a,sk_b,sk_c,sk_d,thet,J_a,J_b,J_c,J_d,g)
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

            if(.not. sk_a) then
                !$omp atomic update
                grad(1,ia) = grad(1,ia) + J_a(1) * g
                !$omp atomic update
                grad(2,ia) = grad(2,ia) + J_a(2) * g
                !$omp atomic update
                grad(3,ia) = grad(3,ia) + J_a(3) * g
            end if

            if(.not. sk_b) then
                !$omp atomic update
                grad(1,ib) = grad(1,ib) + J_b(1) * g
                !$omp atomic update
                grad(2,ib) = grad(2,ib) + J_b(2) * g
                !$omp atomic update
                grad(3,ib) = grad(3,ib) + J_b(3) * g
            end if

            if(.not. sk_c) then
                !$omp atomic update
                grad(1,ic) = grad(1,ic) + J_c(1) * g
                !$omp atomic update
                grad(2,ic) = grad(2,ic) + J_c(2) * g
                !$omp atomic update
                grad(3,ic) = grad(3,ic) + J_c(3) * g
            end if

            if(.not. sk_d) then
                !$omp atomic update
                grad(1,id) = grad(1,id) + J_d(1) * g
                !$omp atomic update
                grad(2,id) = grad(2,id) + J_d(2) * g
                !$omp atomic update
                grad(3,id) = grad(3,id) + J_d(3) * g
            end if
        end do
    end subroutine opb_geomgrad

    subroutine opb_geomhess(bds, hess)
        !! Compute the Hessian of the out-of-plane bending terms. Same
        !! g=dU/dtheta as opb_geomgrad, and h=d^2U/dtheta^2 =
        !! 2k(1+3k^(3)*theta+6k^(4)*theta^2+10k^(5)*theta^3+15k^(6)*theta^4)
        !! (identical functional form to angle_geomhess, since theta already
        !! measures the deviation from the implicit equilibrium of zero).
        !! For any two atoms X,Y of the term: H_XY = h*J_X*J_Y^T + g*H_XY(theta),
        !! with H_XY(theta) from opb_angle_hessian.
        use mod_jacobian_mat, only: opb_angle_hessian

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        !! Bonded potential data structure
        real(rp), intent(inout) :: hess(3,3,bds%top%mm_atoms,bds%top%mm_atoms)
        !! Hessian of out-of-plane bending terms of potential energy

        real(rp) :: thet, g, h, J_a(3), J_b(3), J_c(3), J_d(3)
        real(rp), dimension(3,3) :: Haa, Hab, Hac, Had, Hbb, Hbc, Hbd, Hcc, Hcd, Hdd
        integer(ip) :: i, ia, ib, ic, id, p
        logical :: sk_a, sk_b, sk_c, sk_d

        if(.not. bds%use_opb) return

        !$omp parallel do default(shared) schedule(dynamic)&
        !$omp private(i,ia,ib,ic,id,sk_a,sk_b,sk_c,sk_d,thet,J_a,J_b,J_c,J_d,g,h,p) &
        !$omp private(Haa,Hab,Hac,Had,Hbb,Hbc,Hbd,Hcc,Hcd,Hdd)
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

            call opb_angle_hessian(bds%top%cmm(:,ia), bds%top%cmm(:,ib), &
                                   bds%top%cmm(:,ic), bds%top%cmm(:,id), &
                                   thet, J_a, J_b, J_c, J_d, &
                                   Haa, Hab, Hac, Had, Hbb, Hbc, Hbd, Hcc, Hcd, Hdd)

            g = bds%kopb(i) * thet * (2.0 + 3.0*bds%opb_cubic*thet &
                + 4.0*bds%opb_quartic*thet**2 + 5.0*bds%opb_pentic*thet**3 &
                + 6.0*bds%opb_sextic*thet**4)
            h = 2.0 * bds%kopb(i) * (1.0 + 3.0*bds%opb_cubic*thet &
                + 6.0*bds%opb_quartic*thet**2 + 10.0*bds%opb_pentic*thet**3 &
                + 15.0*bds%opb_sextic*thet**4)

            do p=1,3
                Haa(p,:) = h*J_a(p)*J_a + g*Haa(p,:)
                Hab(p,:) = h*J_a(p)*J_b + g*Hab(p,:)
                Hac(p,:) = h*J_a(p)*J_c + g*Hac(p,:)
                Had(p,:) = h*J_a(p)*J_d + g*Had(p,:)
                Hbb(p,:) = h*J_b(p)*J_b + g*Hbb(p,:)
                Hbc(p,:) = h*J_b(p)*J_c + g*Hbc(p,:)
                Hbd(p,:) = h*J_b(p)*J_d + g*Hbd(p,:)
                Hcc(p,:) = h*J_c(p)*J_c + g*Hcc(p,:)
                Hcd(p,:) = h*J_c(p)*J_d + g*Hcd(p,:)
                Hdd(p,:) = h*J_d(p)*J_d + g*Hdd(p,:)
            end do

            !$omp critical
            if(.not. sk_a) hess(:,:,ia,ia) = hess(:,:,ia,ia) + Haa
            if(.not. sk_b) hess(:,:,ib,ib) = hess(:,:,ib,ib) + Hbb
            if(.not. sk_c) hess(:,:,ic,ic) = hess(:,:,ic,ic) + Hcc
            if(.not. sk_d) hess(:,:,id,id) = hess(:,:,id,id) + Hdd
            if(.not. sk_a .and. .not. sk_b) then
                hess(:,:,ia,ib) = hess(:,:,ia,ib) + Hab
                hess(:,:,ib,ia) = hess(:,:,ib,ia) + transpose(Hab)
            end if
            if(.not. sk_a .and. .not. sk_c) then
                hess(:,:,ia,ic) = hess(:,:,ia,ic) + Hac
                hess(:,:,ic,ia) = hess(:,:,ic,ia) + transpose(Hac)
            end if
            if(.not. sk_a .and. .not. sk_d) then
                hess(:,:,ia,id) = hess(:,:,ia,id) + Had
                hess(:,:,id,ia) = hess(:,:,id,ia) + transpose(Had)
            end if
            if(.not. sk_b .and. .not. sk_c) then
                hess(:,:,ib,ic) = hess(:,:,ib,ic) + Hbc
                hess(:,:,ic,ib) = hess(:,:,ic,ib) + transpose(Hbc)
            end if
            if(.not. sk_b .and. .not. sk_d) then
                hess(:,:,ib,id) = hess(:,:,ib,id) + Hbd
                hess(:,:,id,ib) = hess(:,:,id,ib) + transpose(Hbd)
            end if
            if(.not. sk_c .and. .not. sk_d) then
                hess(:,:,ic,id) = hess(:,:,ic,id) + Hcd
                hess(:,:,id,ic) = hess(:,:,id,ic) + transpose(Hcd)
            end if
            !$omp end critical
        end do
    end subroutine opb_geomhess

    
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

        call mallocate('pitors_init [pitorsat]', 6_ip, n, bds%pitorsat)
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
        
        !$omp parallel do default(shared) reduction(+:V) &
        !$omp private(i,a,b,c,d,e,f,plv1,plv2,pln1,pln2,t,u,cd,thet,costhet)
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

        !$omp parallel do default(shared) schedule(dynamic) &
        !$omp private(i,ia,ib,ic,id,ie,if_,sk_a,sk_b,sk_c,sk_d,sk_e,sk_f) &
        !$omp private(J_a,J_b,J_c,J_d,J_e,J_f,g,thet)
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

        if(.not. sk_a) then
            !$omp atomic update
            grad(1,ia) = grad(1,ia) + g * J_a(1)
            !$omp atomic update
            grad(2,ia) = grad(2,ia) + g * J_a(2)
            !$omp atomic update
            grad(3,ia) = grad(3,ia) + g * J_a(3)
        end if

        if(.not. sk_b) then
            !$omp atomic update
            grad(1,ib) = grad(1,ib) + g * J_b(1)
            !$omp atomic update
            grad(2,ib) = grad(2,ib) + g * J_b(2)
            !$omp atomic update
            grad(3,ib) = grad(3,ib) + g * J_b(3)
        end if

        if(.not. sk_c) then
            !$omp atomic update
            grad(1,ic) = grad(1,ic) + g * J_c(1)
            !$omp atomic update
            grad(2,ic) = grad(2,ic) + g * J_c(2)
            !$omp atomic update
            grad(3,ic) = grad(3,ic) + g * J_c(3)
        end if

        if(.not. sk_d) then
            !$omp atomic update
            grad(1,id) = grad(1,id) + g * J_d(1)
            !$omp atomic update
            grad(2,id) = grad(2,id) + g * J_d(2)
            !$omp atomic update
            grad(3,id) = grad(3,id) + g * J_d(3)
        end if

        if(.not. sk_e) then
            !$omp atomic update
            grad(1,ie) = grad(1,ie) + g * J_e(1)
            !$omp atomic update
            grad(2,ie) = grad(2,ie) + g * J_e(2)
            !$omp atomic update
            grad(3,ie) = grad(3,ie) + g * J_e(3)
        end if

        if(.not. sk_f) then
            !$omp atomic update
            grad(1,if_) = grad(1,if_) + g * J_f(1)
            !$omp atomic update
            grad(2,if_) = grad(2,if_) + g * J_f(2)
            !$omp atomic update
            grad(3,if_) = grad(3,if_) + g * J_f(3)
        end if
        end do
    end subroutine pitors_geomgrad

    subroutine pitors_geomhess(bds, hess)
        !! Compute the Hessian of the pi-torsion potential. With g as in
        !! pitors_geomgrad and h = d^2U/dtheta^2 = -4*k*cos(2*theta-pi), for
        !! any two atoms X,Y of the term: H_XY = h*J_X*J_Y^T + g*H_XY(theta),
        !! with H_XY(theta) from pitors_angle_hessian (atom order
        !! 1=A,2=B,3=C,4=D,5=E,6=F matching pitors_geomgrad's ia,ib,ic,id,ie,if_).
        use mod_jacobian_mat, only: pitors_angle_hessian
        use mod_constants, only : pi

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        !! Bonded potential data structure
        real(rp), intent(inout) :: hess(3,3,bds%top%mm_atoms,bds%top%mm_atoms)
        !! Hessian of pi-torsion terms of potential energy

        real(rp) :: thet, g, h, J_a(3), J_b(3), J_c(3), J_d(3), J_e(3), J_f(3)
        real(rp) :: Hblk(3,3,6,6), Jall(3,6), block_(3,3)
        integer(ip) :: i, k, l, iat(6)
        logical :: sk(6)

        if(.not. bds%use_pitors) return

        !$omp parallel do default(shared) schedule(dynamic) &
        !$omp private(i,iat,sk,J_a,J_b,J_c,J_d,J_e,J_f,g,h,thet,Hblk,Jall,k,l,block_)
        do i=1, bds%npitors
            iat(1) = bds%pitorsat(1,i)
            iat(4) = bds%pitorsat(3,i)
            iat(2) = bds%pitorsat(4,i)
            iat(3) = bds%pitorsat(2,i)
            iat(5) = bds%pitorsat(5,i)
            iat(6) = bds%pitorsat(6,i)

            if(bds%top%use_frozen) then
                do k=1,6
                    sk(k) = bds%top%frozen(iat(k))
                end do
                if(all(sk)) cycle
            else
                sk = .false.
            end if

            call pitors_angle_hessian(bds%top%cmm(:,iat(1)), bds%top%cmm(:,iat(2)), &
                                      bds%top%cmm(:,iat(3)), bds%top%cmm(:,iat(4)), &
                                      bds%top%cmm(:,iat(5)), bds%top%cmm(:,iat(6)), &
                                      thet, J_a, J_b, J_c, J_d, J_e, J_f, Hblk)
            Jall(:,1) = J_a; Jall(:,2) = J_b; Jall(:,3) = J_c
            Jall(:,4) = J_d; Jall(:,5) = J_e; Jall(:,6) = J_f

            g = -2.0 * bds%kpitors(i) * sin(2.0*thet-pi)
            h = -4.0 * bds%kpitors(i) * cos(2.0*thet-pi)

            !$omp critical
            do k=1,6
                if(sk(k)) cycle
                do l=1,6
                    if(sk(l)) cycle
                    block_ = h*outer3(Jall(:,k), Jall(:,l)) + g*Hblk(:,:,k,l)
                    hess(:,:,iat(k),iat(l)) = hess(:,:,iat(k),iat(l)) + block_
                end do
            end do
            !$omp end critical
        end do

    contains
        pure function outer3(u, v) result(m)
            real(rp), intent(in) :: u(3), v(3)
            real(rp) :: m(3,3)
            integer :: p
            do p=1,3
                m(p,:) = u(p)*v
            end do
        end function
    end subroutine pitors_geomhess

    
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

        call mallocate('torsion_init [torsionat]', 4_ip, n, bds%torsionat)
        call mallocate('torsion_init [torsamp]', 6_ip, n, bds%torsamp)
        call mallocate('torsion_init [torsphase]', 6_ip, n, bds%torsphase)
        call mallocate('torsion_init [torsn]', 6_ip, n, bds%torsn)

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

        !$omp parallel do default(shared) &
        !$omp private(i,costhet,thet,j) reduction(+:V)
        do i=1, bds%ntorsion
            ! Atoms that defines the dihedral angle
            costhet = cos_torsion(bds%top, bds%torsionat(:,i))
            
            if(costhet + 1.0 <= eps_rp) then
                thet = pi
            else if(abs(costhet - 1.0) <= eps_rp) then
                thet = 0.0
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

        !$omp parallel do default(shared) &
        !$omp private(i,ia,ib,ic,id,sk_a,sk_b,sk_c,sk_d,j,thet,J_a,J_b,J_c,J_d,g)
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
                if(.not. sk_a) then
                    !$omp atomic update
                    grad(1, ia) = grad(1, ia) + J_a(1) * g
                    !$omp atomic update
                    grad(2, ia) = grad(2, ia) + J_a(2) * g
                    !$omp atomic update
                    grad(3, ia) = grad(3, ia) + J_a(3) * g
                end if
                if(.not. sk_b) then
                    !$omp atomic update
                    grad(1, ib) = grad(1, ib) + J_b(1) * g
                    !$omp atomic update
                    grad(2, ib) = grad(2, ib) + J_b(2) * g
                    !$omp atomic update
                    grad(3, ib) = grad(3, ib) + J_b(3) * g
                end if
                if(.not. sk_c) then
                    !$omp atomic update
                    grad(1, ic) = grad(1, ic) + J_c(1) * g
                    !$omp atomic update
                    grad(2, ic) = grad(2, ic) + J_c(2) * g
                    !$omp atomic update
                    grad(3, ic) = grad(3, ic) + J_c(3) * g
                end if
                if(.not. sk_d) then
                    !$omp atomic update
                    grad(1, id) = grad(1, id) + J_d(1) * g
                    !$omp atomic update
                    grad(2, id) = grad(2, id) + J_d(2) * g
                    !$omp atomic update
                    grad(3, id) = grad(3, id) + J_d(3) * g
                end if
            end do
        end do

    end subroutine torsion_geomgrad

    subroutine torsion_geomhess(bds, hess)
        !! Compute the Hessian of the torsion potential. With g as in
        !! torsion_geomgrad (summed over the active Fourier terms) and
        !! h = d^2U/dtheta^2 = sum_j -n_j^2*A_j*cos(n_j*theta-phi_j), for any
        !! two atoms X,Y of the term: H_XY = h*J_X*J_Y^T + g*H_XY(theta),
        !! with H_XY(theta) from torsion_angle_hessian.
        use mod_jacobian_mat, only: torsion_angle_hessian

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        !! Bonded potential data structure
        real(rp), intent(inout) :: hess(3,3,bds%top%mm_atoms,bds%top%mm_atoms)
        !! Hessian of torsion terms of potential energy

        real(rp) :: thet, g, h, n, J_a(3), J_b(3), J_c(3), J_d(3)
        real(rp), dimension(3,3) :: Haa, Hab, Hac, Had, Hbb, Hbc, Hbd, Hcc, Hcd, Hdd
        integer(ip) :: i, j, ia, ib, ic, id, p
        logical :: sk_a, sk_b, sk_c, sk_d

        if(.not. bds%use_torsion) return

        !$omp parallel do default(shared) &
        !$omp private(i,ia,ib,ic,id,sk_a,sk_b,sk_c,sk_d,j,n,thet,J_a,J_b,J_c,J_d,g,h,p) &
        !$omp private(Haa,Hab,Hac,Had,Hbb,Hbc,Hbd,Hcc,Hcd,Hdd)
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

            call torsion_angle_hessian(bds%top%cmm(:,ia), bds%top%cmm(:,ib), &
                                       bds%top%cmm(:,ic), bds%top%cmm(:,id), &
                                       thet, J_a, J_b, J_c, J_d, &
                                       Haa, Hab, Hac, Had, Hbb, Hbc, Hbd, Hcc, Hcd, Hdd)

            g = 0.0_rp
            h = 0.0_rp
            do j=1, 6
                if(bds%torsn(j,i) < 1) exit
                n = real(bds%torsn(j,i), rp)
                g = g - n * sin(n*thet - bds%torsphase(j,i)) * bds%torsamp(j,i)
                h = h - n**2 * cos(n*thet - bds%torsphase(j,i)) * bds%torsamp(j,i)
            end do

            do p=1,3
                Haa(p,:) = h*J_a(p)*J_a + g*Haa(p,:)
                Hab(p,:) = h*J_a(p)*J_b + g*Hab(p,:)
                Hac(p,:) = h*J_a(p)*J_c + g*Hac(p,:)
                Had(p,:) = h*J_a(p)*J_d + g*Had(p,:)
                Hbb(p,:) = h*J_b(p)*J_b + g*Hbb(p,:)
                Hbc(p,:) = h*J_b(p)*J_c + g*Hbc(p,:)
                Hbd(p,:) = h*J_b(p)*J_d + g*Hbd(p,:)
                Hcc(p,:) = h*J_c(p)*J_c + g*Hcc(p,:)
                Hcd(p,:) = h*J_c(p)*J_d + g*Hcd(p,:)
                Hdd(p,:) = h*J_d(p)*J_d + g*Hdd(p,:)
            end do

            !$omp critical
            if(.not. sk_a) hess(:,:,ia,ia) = hess(:,:,ia,ia) + Haa
            if(.not. sk_b) hess(:,:,ib,ib) = hess(:,:,ib,ib) + Hbb
            if(.not. sk_c) hess(:,:,ic,ic) = hess(:,:,ic,ic) + Hcc
            if(.not. sk_d) hess(:,:,id,id) = hess(:,:,id,id) + Hdd
            if(.not. sk_a .and. .not. sk_b) then
                hess(:,:,ia,ib) = hess(:,:,ia,ib) + Hab
                hess(:,:,ib,ia) = hess(:,:,ib,ia) + transpose(Hab)
            end if
            if(.not. sk_a .and. .not. sk_c) then
                hess(:,:,ia,ic) = hess(:,:,ia,ic) + Hac
                hess(:,:,ic,ia) = hess(:,:,ic,ia) + transpose(Hac)
            end if
            if(.not. sk_a .and. .not. sk_d) then
                hess(:,:,ia,id) = hess(:,:,ia,id) + Had
                hess(:,:,id,ia) = hess(:,:,id,ia) + transpose(Had)
            end if
            if(.not. sk_b .and. .not. sk_c) then
                hess(:,:,ib,ic) = hess(:,:,ib,ic) + Hbc
                hess(:,:,ic,ib) = hess(:,:,ic,ib) + transpose(Hbc)
            end if
            if(.not. sk_b .and. .not. sk_d) then
                hess(:,:,ib,id) = hess(:,:,ib,id) + Hbd
                hess(:,:,id,ib) = hess(:,:,id,ib) + transpose(Hbd)
            end if
            if(.not. sk_c .and. .not. sk_d) then
                hess(:,:,ic,id) = hess(:,:,ic,id) + Hcd
                hess(:,:,id,ic) = hess(:,:,id,ic) + transpose(Hcd)
            end if
            !$omp end critical
        end do

    end subroutine torsion_geomhess

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
        real(rp), intent(inout) :: grad(3, bds%top%mm_atoms)
        !! improper torsion potential, result will be added to V
        real(rp) :: thet, g, J_a(3), J_b(3), J_c(3), J_d(3)
        integer(ip) :: i, j, ia, ib, ic, id
        logical :: sk_a, sk_b, sk_c, sk_d

        if (.not. bds%use_imptorsion) return

        !$omp parallel do default(shared) &
        !$omp private(i, ia, ib, ic, id, sk_a, sk_b, sk_c, sk_d, j, thet, J_a, J_b, J_c, J_d, g)
        do i = 1, bds%nimptorsion
            ! Atoms that define the dihedral angle
            ia = bds%imptorsionat(1, i)
            ib = bds%imptorsionat(2, i)
            ic = bds%imptorsionat(3, i)
            id = bds%imptorsionat(4, i)

            if (bds%top%use_frozen) then
                sk_a = bds%top%frozen(ia)
                sk_b = bds%top%frozen(ib)
                sk_c = bds%top%frozen(ic)
                sk_d = bds%top%frozen(id)
                if (sk_a .and. sk_b .and. sk_c .and. sk_d) cycle
            else
                sk_a = .false.
                sk_b = .false.
                sk_c = .false.
                sk_d = .false.
            end if

            call torsion_angle_jacobian(bds%top%cmm(:, ia), &
                                        bds%top%cmm(:, ib), &
                                        bds%top%cmm(:, ic), &
                                        bds%top%cmm(:, id), &
                                        thet, J_a, J_b, J_c, J_d)

            do j = 1, 3
                if (bds%imptorsn(j, i) < 1) exit
                g = -real(bds%imptorsn(j, i)) * sin(real(bds%imptorsn(j, i)) * thet &
                                                    - bds%imptorsphase(j, i)) &
                                                * bds%imptorsamp(j, i)
                if (.not. sk_a) then
                    !$omp atomic update
                    grad(1, ia) = grad(1, ia) + J_a(1) * g
                    !$omp atomic update
                    grad(2, ia) = grad(2, ia) + J_a(2) * g
                    !$omp atomic update
                    grad(3, ia) = grad(3, ia) + J_a(3) * g
                end if
                if (.not. sk_b) then
                    !$omp atomic update
                    grad(1, ib) = grad(1, ib) + J_b(1) * g
                    !$omp atomic update
                    grad(2, ib) = grad(2, ib) + J_b(2) * g
                    !$omp atomic update
                    grad(3, ib) = grad(3, ib) + J_b(3) * g
                end if
                if (.not. sk_c) then
                    !$omp atomic update
                    grad(1, ic) = grad(1, ic) + J_c(1) * g
                    !$omp atomic update
                    grad(2, ic) = grad(2, ic) + J_c(2) * g
                    !$omp atomic update
                    grad(3, ic) = grad(3, ic) + J_c(3) * g
                end if
                if (.not. sk_d) then
                    !$omp atomic update
                    grad(1, id) = grad(1, id) + J_d(1) * g
                    !$omp atomic update
                    grad(2, id) = grad(2, id) + J_d(2) * g
                    !$omp atomic update
                    grad(3, id) = grad(3, id) + J_d(3) * g
                end if
            end do
        end do
    end subroutine imptorsion_geomgrad

    subroutine imptorsion_geomhess(bds, hess)
        !! Compute the Hessian of the improper torsion potential. Identical
        !! in structure to torsion_geomhess (see there), only over the up to
        !! 3 imptorsn/imptorsamp/imptorsphase Fourier terms.
        use mod_jacobian_mat, only: torsion_angle_hessian

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        !! Bonded potential data structure
        real(rp), intent(inout) :: hess(3,3,bds%top%mm_atoms,bds%top%mm_atoms)
        !! Hessian of improper torsion terms of potential energy

        real(rp) :: thet, g, h, n, J_a(3), J_b(3), J_c(3), J_d(3)
        real(rp), dimension(3,3) :: Haa, Hab, Hac, Had, Hbb, Hbc, Hbd, Hcc, Hcd, Hdd
        integer(ip) :: i, j, ia, ib, ic, id, p
        logical :: sk_a, sk_b, sk_c, sk_d

        if(.not. bds%use_imptorsion) return

        !$omp parallel do default(shared) &
        !$omp private(i,ia,ib,ic,id,sk_a,sk_b,sk_c,sk_d,j,n,thet,J_a,J_b,J_c,J_d,g,h,p) &
        !$omp private(Haa,Hab,Hac,Had,Hbb,Hbc,Hbd,Hcc,Hcd,Hdd)
        do i=1, bds%nimptorsion
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

            call torsion_angle_hessian(bds%top%cmm(:,ia), bds%top%cmm(:,ib), &
                                       bds%top%cmm(:,ic), bds%top%cmm(:,id), &
                                       thet, J_a, J_b, J_c, J_d, &
                                       Haa, Hab, Hac, Had, Hbb, Hbc, Hbd, Hcc, Hcd, Hdd)

            g = 0.0_rp
            h = 0.0_rp
            do j=1, 3
                if(bds%imptorsn(j,i) < 1) exit
                n = real(bds%imptorsn(j,i), rp)
                g = g - n * sin(n*thet - bds%imptorsphase(j,i)) * bds%imptorsamp(j,i)
                h = h - n**2 * cos(n*thet - bds%imptorsphase(j,i)) * bds%imptorsamp(j,i)
            end do

            do p=1,3
                Haa(p,:) = h*J_a(p)*J_a + g*Haa(p,:)
                Hab(p,:) = h*J_a(p)*J_b + g*Hab(p,:)
                Hac(p,:) = h*J_a(p)*J_c + g*Hac(p,:)
                Had(p,:) = h*J_a(p)*J_d + g*Had(p,:)
                Hbb(p,:) = h*J_b(p)*J_b + g*Hbb(p,:)
                Hbc(p,:) = h*J_b(p)*J_c + g*Hbc(p,:)
                Hbd(p,:) = h*J_b(p)*J_d + g*Hbd(p,:)
                Hcc(p,:) = h*J_c(p)*J_c + g*Hcc(p,:)
                Hcd(p,:) = h*J_c(p)*J_d + g*Hcd(p,:)
                Hdd(p,:) = h*J_d(p)*J_d + g*Hdd(p,:)
            end do

            !$omp critical
            if(.not. sk_a) hess(:,:,ia,ia) = hess(:,:,ia,ia) + Haa
            if(.not. sk_b) hess(:,:,ib,ib) = hess(:,:,ib,ib) + Hbb
            if(.not. sk_c) hess(:,:,ic,ic) = hess(:,:,ic,ic) + Hcc
            if(.not. sk_d) hess(:,:,id,id) = hess(:,:,id,id) + Hdd
            if(.not. sk_a .and. .not. sk_b) then
                hess(:,:,ia,ib) = hess(:,:,ia,ib) + Hab
                hess(:,:,ib,ia) = hess(:,:,ib,ia) + transpose(Hab)
            end if
            if(.not. sk_a .and. .not. sk_c) then
                hess(:,:,ia,ic) = hess(:,:,ia,ic) + Hac
                hess(:,:,ic,ia) = hess(:,:,ic,ia) + transpose(Hac)
            end if
            if(.not. sk_a .and. .not. sk_d) then
                hess(:,:,ia,id) = hess(:,:,ia,id) + Had
                hess(:,:,id,ia) = hess(:,:,id,ia) + transpose(Had)
            end if
            if(.not. sk_b .and. .not. sk_c) then
                hess(:,:,ib,ic) = hess(:,:,ib,ic) + Hbc
                hess(:,:,ic,ib) = hess(:,:,ic,ib) + transpose(Hbc)
            end if
            if(.not. sk_b .and. .not. sk_d) then
                hess(:,:,ib,id) = hess(:,:,ib,id) + Hbd
                hess(:,:,id,ib) = hess(:,:,id,ib) + transpose(Hbd)
            end if
            if(.not. sk_c .and. .not. sk_d) then
                hess(:,:,ic,id) = hess(:,:,ic,id) + Hcd
                hess(:,:,id,ic) = hess(:,:,id,ic) + transpose(Hcd)
            end if
            !$omp end critical
        end do

    end subroutine imptorsion_geomhess
    
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

        call mallocate('imptorsion_init [imptorsionat]', 4_ip, n, bds%imptorsionat)
        call mallocate('imptorsion_init [imptorsamp]', 3_ip, n, bds%imptorsamp)
        call mallocate('imptorsion_init [imptorsphase]', 3_ip, n, bds%imptorsphase)
        call mallocate('imptorsion_init [imptorsn]', 3_ip, n, bds%imptorsn)

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

        call mallocate('angtor_init [angtorat]', 4_ip, n, bds%angtorat)
        call mallocate('angtor_init [angtork]', 6_ip, n, bds%angtork)
        call mallocate('angtor_init [angtor_t]', n, bds%angtor_t)
        call mallocate('angtor_init [angtor_a]', 2_ip, n, bds%angtor_a)

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

        call mallocate('strtor_init [strtorat]', 4_ip, n, bds%strtorat)
        call mallocate('strtor_init [strtork]', 9_ip, n, bds%strtork)
        call mallocate('strtor_init [strtor_t]', n, bds%strtor_t)
        call mallocate('strtor_init [strtor_a]', 3_ip, n, bds%strtor_b)

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

        !$omp parallel do default(shared) reduction(+:V) &
        !$omp private(i,costhet,thet,j,dihef,ia1,ia2,dr1,dr2,l1,l2,angle1,angle2,delta_a,vat,k)
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
        real(rp) :: thet, gt(3), dihef(3), da1, da2, angle1, angle2, f1, f2, f3, &
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

        !$omp parallel do default(shared) &
        !$omp private(thet, gt, dihef, da1, da2, angle1, angle2, f1, f2, f3, Jt_a, Jt_b) &
        !$omp private(Jt_c, Jt_d, Ja1_a, Ja1_b, Ja1_c) &
        !$omp private(Ja2_a, Ja2_b, Ja2_c, i, j, k, ia1, ia2) &
        !$omp private(it_a, it_b, it_c, it_d, ia1_a, ia1_b, ia1_c, ia2_a, ia2_b, ia2_c) &
        !$omp private(sk_ta, sk_tb, sk_tc, sk_td, sk_1a, sk_1b, sk_1c, sk_2a, sk_2b, sk_2c)
        do i=1, bds%nangtor
            ! Atoms that define the dihedral angle
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
 
            do k = 1, 3
                if(.not.(sk_ta .and. sk_tb .and. sk_tc .and. sk_td)) &
                    f1 = (bds%angtork(k, i) * da1 + bds%angtork(3+k,i) * da2) * gt(k)
                if(.not.(sk_1a .and. sk_1b .and. sk_1c)) &
                    f2 = bds%angtork(k, i) * dihef(k)
                if(.not.(sk_2a .and. sk_2b .and. sk_2c)) &
                    f3 = bds%angtork(3+k, i) * dihef(k)

                if (.not. sk_ta) then
                    !$omp atomic update
                    grad(1, it_a) = grad(1, it_a) + f1 * Jt_a(1)
                    !$omp atomic update
                    grad(2, it_a) = grad(2, it_a) + f1 * Jt_a(2)
                    !$omp atomic update
                    grad(3, it_a) = grad(3, it_a) + f1 * Jt_a(3)
                end if
                
                if (.not. sk_tb) then
                    !$omp atomic update
                    grad(1, it_b) = grad(1, it_b) + f1 * Jt_b(1)
                    !$omp atomic update
                    grad(2, it_b) = grad(2, it_b) + f1 * Jt_b(2)
                    !$omp atomic update
                    grad(3, it_b) = grad(3, it_b) + f1 * Jt_b(3)
                end if
                if (.not. sk_tc) then
                    !$omp atomic update
                    grad(1, it_c) = grad(1, it_c) + f1 * Jt_c(1)
                    !$omp atomic update
                    grad(2, it_c) = grad(2, it_c) + f1 * Jt_c(2)
                    !$omp atomic update
                    grad(3, it_c) = grad(3, it_c) + f1 * Jt_c(3)
                end if
                if (.not. sk_td) then
                    !$omp atomic update
                    grad(1, it_d) = grad(1, it_d) + f1 * Jt_d(1)
                    !$omp atomic update
                    grad(2, it_d) = grad(2, it_d) + f1 * Jt_d(2)
                    !$omp atomic update
                    grad(3, it_d) = grad(3, it_d) + f1 * Jt_d(3)
                end if

                if (.not. sk_1a) then
                    !$omp atomic update
                    grad(1, ia1_a) = grad(1, ia1_a) + f2 * Ja1_a(1)
                    !$omp atomic update
                    grad(2, ia1_a) = grad(2, ia1_a) + f2 * Ja1_a(2)
                    !$omp atomic update
                    grad(3, ia1_a) = grad(3, ia1_a) + f2 * Ja1_a(3)
                end if
                if (.not. sk_1b) then
                    !$omp atomic update
                    grad(1, ia1_b) = grad(1, ia1_b) + f2 * Ja1_b(1)
                    !$omp atomic update
                    grad(2, ia1_b) = grad(2, ia1_b) + f2 * Ja1_b(2)
                    !$omp atomic update
                    grad(3, ia1_b) = grad(3, ia1_b) + f2 * Ja1_b(3)
                end if
                if (.not. sk_1c) then
                    !$omp atomic update
                    grad(1, ia1_c) = grad(1, ia1_c) + f2 * Ja1_c(1)
                    !$omp atomic update
                    grad(2, ia1_c) = grad(2, ia1_c) + f2 * Ja1_c(2)
                    !$omp atomic update
                    grad(3, ia1_c) = grad(3, ia1_c) + f2 * Ja1_c(3)
                end if

                if (.not. sk_2a) then
                    !$omp atomic update
                    grad(1, ia2_a) = grad(1, ia2_a) + f3 * Ja2_a(1)
                    !$omp atomic update
                    grad(2, ia2_a) = grad(2, ia2_a) + f3 * Ja2_a(2)
                    !$omp atomic update
                    grad(3, ia2_a) = grad(3, ia2_a) + f3 * Ja2_a(3)
                end if
                if (.not. sk_2b) then
                    !$omp atomic update
                    grad(1, ia2_b) = grad(1, ia2_b) + f3 * Ja2_b(1)
                    !$omp atomic update
                    grad(2, ia2_b) = grad(2, ia2_b) + f3 * Ja2_b(2)
                    !$omp atomic update
                    grad(3, ia2_b) = grad(3, ia2_b) + f3 * Ja2_b(3)
                end if
                if (.not. sk_2c) then
                    !$omp atomic update
                    grad(1, ia2_c) = grad(1, ia2_c) + f3 * Ja2_c(1)
                    !$omp atomic update
                    grad(2, ia2_c) = grad(2, ia2_c) + f3 * Ja2_c(2)
                    !$omp atomic update
                    grad(3, ia2_c) = grad(3, ia2_c) + f3 * Ja2_c(3)
                end if
            end do
        end do
    end subroutine angtor_geomgrad

    subroutine angtor_geomhess(bds, hess)
        !! Compute the Hessian of the angle-torsion coupling term. Treating
        !! U as a function of the three internal coordinates
        !! (theta=torsion, alpha1, alpha2), the only nonzero internal second
        !! derivatives are d^2U/dtheta^2 = f11, d^2U/dtheta/dalpha1 = f12 and
        !! d^2U/dtheta/dalpha2 = f13 (d^2U/dalpha_m^2 = d^2U/dalpha1/dalpha2
        !! = 0, since U is linear in each alpha_m separately), giving, for
        !! any two atoms X,Y of the term (each belonging to the torsion
        !! and/or one of the two angles):
        !! \[ H_{XY} = f_{11}J^\theta_XJ^{\theta\dagger}_Y
        !!    + f_{12}\left(J^\theta_XJ^{\alpha_1\dagger}_Y+J^{\alpha_1}_XJ^{\theta\dagger}_Y\right)
        !!    + f_{13}\left(J^\theta_XJ^{\alpha_2\dagger}_Y+J^{\alpha_2}_XJ^{\theta\dagger}_Y\right)
        !!    + f_1H^\theta_{XY}+f_2H^{\alpha_1}_{XY}+f_3H^{\alpha_2}_{XY} \]
        !! where the last three terms are only present when X,Y both belong
        !! to the torsion (resp. angle1, angle2), with f1,f2,f3 as in
        !! angtor_geomgrad and H^theta from torsion_angle_hessian,
        !! H^alpha1/2 from simple_angle_hessian.
        use mod_jacobian_mat, only: simple_angle_hessian, torsion_angle_hessian

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        !! Bonded potential data structure
        real(rp), intent(inout) :: hess(3,3,bds%top%mm_atoms,bds%top%mm_atoms)
        !! Hessian of angle-torsion terms of potential energy

        real(rp) :: thet, gt(3), ht_(3), dihef(3), da1, da2, angle1, angle2, &
                    f1, f2, f3, f11, f12, f13
        real(rp), dimension(3) :: Jt_a, Jt_b, Jt_c, Jt_d, Ja1_a, Ja1_b, Ja1_c, &
                                  Ja2_a, Ja2_b, Ja2_c
        real(rp), dimension(3,3) :: Haa, Hab, Hac, Had, Hbb, Hbc, Hbd, Hcc, Hcd, Hdd, &
                                    H1aa, H1ab, H1ac, H1bb, H1bc, H1cc, &
                                    H2aa, H2ab, H2ac, H2bb, H2bc, H2cc, block_
        real(rp) :: Jslot(3,10), Mmat(3,3)
        real(rp) :: Ht(3,3,4,4), Ha1(3,3,3,3), Ha2(3,3,3,3)
        integer(ip) :: i, j, k, p, q, tp, tq, ap, aq, ia1, ia2, iat(10)
        logical :: sk(10)

        if(.not. bds%use_angtor) return

        !$omp parallel do default(shared) schedule(dynamic) &
        !$omp private(i,j,k,ia1,ia2,iat,sk,thet,gt,ht_,dihef,da1,da2,angle1,angle2) &
        !$omp private(f1,f2,f3,f11,f12,f13,Jt_a,Jt_b,Jt_c,Jt_d,Ja1_a,Ja1_b,Ja1_c) &
        !$omp private(Ja2_a,Ja2_b,Ja2_c,Haa,Hab,Hac,Had,Hbb,Hbc,Hbd,Hcc,Hcd,Hdd) &
        !$omp private(H1aa,H1ab,H1ac,H1bb,H1bc,H1cc,H2aa,H2ab,H2ac,H2bb,H2bc,H2cc) &
        !$omp private(Jslot,Mmat,Ht,Ha1,Ha2,p,q,tp,tq,ap,aq,block_)
        do i=1, bds%nangtor
            iat(1) = bds%angtorat(1,i)
            iat(2) = bds%angtorat(2,i)
            iat(3) = bds%angtorat(3,i)
            iat(4) = bds%angtorat(4,i)

            ia1 = bds%angtor_a(1,i)
            iat(5) = bds%angleat(1,ia1)
            iat(6) = bds%angleat(2,ia1)
            iat(7) = bds%angleat(3,ia1)

            ia2 = bds%angtor_a(2,i)
            iat(8) = bds%angleat(1,ia2)
            iat(9) = bds%angleat(2,ia2)
            iat(10) = bds%angleat(3,ia2)

            if(bds%top%use_frozen) then
                do k=1,10
                    sk(k) = bds%top%frozen(iat(k))
                end do
                if(all(sk)) cycle
            else
                sk = .false.
            end if

            call torsion_angle_hessian(bds%top%cmm(:,iat(1)), bds%top%cmm(:,iat(2)), &
                                       bds%top%cmm(:,iat(3)), bds%top%cmm(:,iat(4)), &
                                       thet, Jt_a, Jt_b, Jt_c, Jt_d, &
                                       Haa, Hab, Hac, Had, Hbb, Hbc, Hbd, Hcc, Hcd, Hdd)
            call simple_angle_hessian(bds%top%cmm(:,iat(5)), bds%top%cmm(:,iat(6)), &
                                      bds%top%cmm(:,iat(7)), angle1, Ja1_a, Ja1_b, Ja1_c, &
                                      H1aa, H1ab, H1ac, H1bb, H1bc, H1cc)
            call simple_angle_hessian(bds%top%cmm(:,iat(8)), bds%top%cmm(:,iat(9)), &
                                      bds%top%cmm(:,iat(10)), angle2, Ja2_a, Ja2_b, Ja2_c, &
                                      H2aa, H2ab, H2ac, H2bb, H2bc, H2cc)

            do j=1,3
                gt(j) = -real(j,rp) * sin(j*thet+bds%torsphase(j,bds%angtor_t(i)))
                ht_(j) = -real(j,rp)**2 * cos(j*thet+bds%torsphase(j,bds%angtor_t(i)))
                dihef(j) = 1.0_rp + cos(j*thet+bds%torsphase(j,bds%angtor_t(i)))
            end do
            da1 = angle1 - bds%eqangle(ia1)
            da2 = angle2 - bds%eqangle(ia2)

            f1 = 0.0_rp; f2 = 0.0_rp; f3 = 0.0_rp
            f11 = 0.0_rp; f12 = 0.0_rp; f13 = 0.0_rp
            do j=1,3
                f1 = f1 + (bds%angtork(j,i)*da1 + bds%angtork(3+j,i)*da2) * gt(j)
                f2 = f2 + bds%angtork(j,i) * dihef(j)
                f3 = f3 + bds%angtork(3+j,i) * dihef(j)
                f11 = f11 + (bds%angtork(j,i)*da1 + bds%angtork(3+j,i)*da2) * ht_(j)
                f12 = f12 + bds%angtork(j,i) * gt(j)
                f13 = f13 + bds%angtork(3+j,i) * gt(j)
            end do

            Jslot(:,1)=Jt_a; Jslot(:,2)=Jt_b; Jslot(:,3)=Jt_c; Jslot(:,4)=Jt_d
            Jslot(:,5)=Ja1_a; Jslot(:,6)=Ja1_b; Jslot(:,7)=Ja1_c
            Jslot(:,8)=Ja2_a; Jslot(:,9)=Ja2_b; Jslot(:,10)=Ja2_c

            Mmat(1,:) = [f11, f12, f13]
            Mmat(2,:) = [f12, 0.0_rp, 0.0_rp]
            Mmat(3,:) = [f13, 0.0_rp, 0.0_rp]

            Ht(:,:,1,1)=Haa; Ht(:,:,1,2)=Hab; Ht(:,:,1,3)=Hac; Ht(:,:,1,4)=Had
            Ht(:,:,2,1)=transpose(Hab); Ht(:,:,2,2)=Hbb; Ht(:,:,2,3)=Hbc; Ht(:,:,2,4)=Hbd
            Ht(:,:,3,1)=transpose(Hac); Ht(:,:,3,2)=transpose(Hbc); Ht(:,:,3,3)=Hcc; Ht(:,:,3,4)=Hcd
            Ht(:,:,4,1)=transpose(Had); Ht(:,:,4,2)=transpose(Hbd); Ht(:,:,4,3)=transpose(Hcd); Ht(:,:,4,4)=Hdd

            Ha1(:,:,1,1)=H1aa; Ha1(:,:,1,2)=H1ab; Ha1(:,:,1,3)=H1ac
            Ha1(:,:,2,1)=transpose(H1ab); Ha1(:,:,2,2)=H1bb; Ha1(:,:,2,3)=H1bc
            Ha1(:,:,3,1)=transpose(H1ac); Ha1(:,:,3,2)=transpose(H1bc); Ha1(:,:,3,3)=H1cc

            Ha2(:,:,1,1)=H2aa; Ha2(:,:,1,2)=H2ab; Ha2(:,:,1,3)=H2ac
            Ha2(:,:,2,1)=transpose(H2ab); Ha2(:,:,2,2)=H2bb; Ha2(:,:,2,3)=H2bc
            Ha2(:,:,3,1)=transpose(H2ac); Ha2(:,:,3,2)=transpose(H2bc); Ha2(:,:,3,3)=H2cc

            !$omp critical
            do p=1,10
                if(sk(p)) cycle
                tp = slot_type(p); ap = slot_local(p)
                do q=1,10
                    if(sk(q)) cycle
                    tq = slot_type(q); aq = slot_local(q)
                    block_ = Mmat(tp,tq) * outer10(Jslot(:,p), Jslot(:,q))
                    if(tp == tq) then
                        if(tp == 1) block_ = block_ + f1*Ht(:,:,ap,aq)
                        if(tp == 2) block_ = block_ + f2*Ha1(:,:,ap,aq)
                        if(tp == 3) block_ = block_ + f3*Ha2(:,:,ap,aq)
                    end if
                    hess(:,:,iat(p),iat(q)) = hess(:,:,iat(p),iat(q)) + block_
                end do
            end do
            !$omp end critical
        end do

    contains
        pure function slot_type(k) result(tt)
            integer(ip), intent(in) :: k
            integer(ip) :: tt
            if(k <= 4) then
                tt = 1
            else if(k <= 7) then
                tt = 2
            else
                tt = 3
            end if
        end function
        pure function slot_local(k) result(ll)
            integer(ip), intent(in) :: k
            integer(ip) :: ll
            if(k <= 4) then
                ll = k
            else if(k <= 7) then
                ll = k - 4
            else
                ll = k - 7
            end if
        end function
        pure function outer10(u, v) result(m)
            real(rp), intent(in) :: u(3), v(3)
            real(rp) :: m(3,3)
            integer :: r
            do r=1,3
                m(r,:) = u(r)*v
            end do
        end function
    end subroutine angtor_geomhess
    
    subroutine strtor_potential(bds, V)
        use mod_constants

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        ! Bonded potential data structure
        real(rp), intent(inout) :: V
        real(rp) :: thet, costhet, dihef(3), dr(3), r(3), vst
        integer(ip) :: i, j, k, ib1, ib2, ib3
        
        if(.not. bds%use_strtor) return

        !$omp parallel do default(shared) reduction(+:V) &
        !$omp private(i,costhet,thet,j,dihef,ib1,ib2,ib3,r,dr,vst,k)
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
        real(rp), intent(inout) :: grad(3, bds%top%mm_atoms)
        !! improper torsion potential, result will be added to V

        real(rp) :: thet, gt(3), dihef(3), dr1, dr2, dr3, r1, r2, r3, &
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

        if (.not. bds%use_strtor) return

        !$omp parallel do default(shared) &
        !$omp private(thet, gt, dihef, dr1, dr2, dr3, r1, r2, r3) &
        !$omp private(Jt_a, Jt_b, Jt_c, Jt_d, Jb1_a, Jb1_b, Jb2_a, Jb2_b) &
        !$omp private(Jb3_a, Jb3_b, i, j, k, ib1, ib2, ib3, it_a, it_b, it_c, it_d) &
        !$omp private(ib1_a, ib1_b, ib2_a, ib2_b, ib3_a, ib3_b, sk_ta, sk_tb, sk_tc, sk_td) &
        !$omp private(sk_1a, sk_1b, sk_2a, sk_2b, sk_3a, sk_3b)
        do i = 1, bds%nstrtor
            ! Atoms that define the dihedral angle
            it_a = bds%strtorat(1, i)
            it_b = bds%strtorat(2, i)
            it_c = bds%strtorat(3, i)
            it_d = bds%strtorat(4, i)

            ib1 = bds%strtor_b(1, i)
            ib1_a = bds%bondat(1, ib1)
            ib1_b = bds%bondat(2, ib1)

            ib2 = bds%strtor_b(2, i)
            ib2_a = bds%bondat(1, ib2)
            ib2_b = bds%bondat(2, ib2)

            ib3 = bds%strtor_b(3, i)
            ib3_a = bds%bondat(1, ib3)
            ib3_b = bds%bondat(2, ib3)

            if (bds%top%use_frozen) then
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

                if (sk_ta .and. sk_tb .and. sk_tc .and. sk_td .and. &
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

            call torsion_angle_jacobian(bds%top%cmm(:, it_a), &
                                        bds%top%cmm(:, it_b), &
                                        bds%top%cmm(:, it_c), &
                                        bds%top%cmm(:, it_d), &
                                        thet, Jt_a, Jt_b, Jt_c, Jt_d)
            do j = 1, 3
                gt(j) = -real(j) * sin(j * thet + bds%torsphase(j, bds%strtor_t(i)))
                dihef(j) = 1.0 + cos(j * thet + bds%torsphase(j, bds%strtor_t(i)))
            end do

            call Rij_jacobian(bds%top%cmm(:, ib1_a), &
                              bds%top%cmm(:, ib1_b), &
                              r1, Jb1_a, Jb1_b)
            dr1 = r1 - bds%l0bond(ib1)

            call Rij_jacobian(bds%top%cmm(:, ib2_a), &
                              bds%top%cmm(:, ib2_b), &
                              r2, Jb2_a, Jb2_b)
            dr2 = r2 - bds%l0bond(ib2)

            call Rij_jacobian(bds%top%cmm(:, ib3_a), &
                              bds%top%cmm(:, ib3_b), &
                              r3, Jb3_a, Jb3_b)
            dr3 = r3 - bds%l0bond(ib3)
        
            do k = 1, 3
                if (.not. sk_ta) then
                    !$omp atomic update
                    grad(1, it_a) = grad(1, it_a) + bds%strtork(k, i) * dr1 * gt(k) * Jt_a(1)
                    !$omp atomic update
                    grad(2, it_a) = grad(2, it_a) + bds%strtork(k, i) * dr1 * gt(k) * Jt_a(2)
                    !$omp atomic update
                    grad(3, it_a) = grad(3, it_a) + bds%strtork(k, i) * dr1 * gt(k) * Jt_a(3)
                end if
                if (.not. sk_tb) then
                    !$omp atomic update
                    grad(1, it_b) = grad(1, it_b) + bds%strtork(k, i) * dr1 * gt(k) * Jt_b(1)
                    !$omp atomic update
                    grad(2, it_b) = grad(2, it_b) + bds%strtork(k, i) * dr1 * gt(k) * Jt_b(2)
                    !$omp atomic update
                    grad(3, it_b) = grad(3, it_b) + bds%strtork(k, i) * dr1 * gt(k) * Jt_b(3)
                end if
                if (.not. sk_tc) then
                    !$omp atomic update
                    grad(1, it_c) = grad(1, it_c) + bds%strtork(k, i) * dr1 * gt(k) * Jt_c(1)
                    !$omp atomic update
                    grad(2, it_c) = grad(2, it_c) + bds%strtork(k, i) * dr1 * gt(k) * Jt_c(2)
                    !$omp atomic update
                    grad(3, it_c) = grad(3, it_c) + bds%strtork(k, i) * dr1 * gt(k) * Jt_c(3)
                end if
                if (.not. sk_td) then
                    !$omp atomic update
                    grad(1, it_d) = grad(1, it_d) + bds%strtork(k, i) * dr1 * gt(k) * Jt_d(1)
                    !$omp atomic update
                    grad(2, it_d) = grad(2, it_d) + bds%strtork(k, i) * dr1 * gt(k) * Jt_d(2)
                    !$omp atomic update
                    grad(3, it_d) = grad(3, it_d) + bds%strtork(k, i) * dr1 * gt(k) * Jt_d(3)
                end if
                if (.not. sk_1a) then
                    !$omp atomic update
                    grad(1, ib1_a) = grad(1, ib1_a) + bds%strtork(k, i) * dihef(k) * Jb1_a(1)
                    !$omp atomic update
                    grad(2, ib1_a) = grad(2, ib1_a) + bds%strtork(k, i) * dihef(k) * Jb1_a(2)
                    !$omp atomic update
                    grad(3, ib1_a) = grad(3, ib1_a) + bds%strtork(k, i) * dihef(k) * Jb1_a(3)
                end if
                if (.not. sk_1b) then
                    !$omp atomic update
                    grad(1, ib1_b) = grad(1, ib1_b) + bds%strtork(k, i) * dihef(k) * Jb1_b(1)
                    !$omp atomic update
                    grad(2, ib1_b) = grad(2, ib1_b) + bds%strtork(k, i) * dihef(k) * Jb1_b(2)
                    !$omp atomic update
                    grad(3, ib1_b) = grad(3, ib1_b) + bds%strtork(k, i) * dihef(k) * Jb1_b(3)
                end if
                if (.not. sk_ta) then
                    !$omp atomic update
                    grad(1, it_a) = grad(1, it_a) + bds%strtork(3 + k, i) * dr2 * gt(k) * Jt_a(1)
                    !$omp atomic update
                    grad(2, it_a) = grad(2, it_a) + bds%strtork(3 + k, i) * dr2 * gt(k) * Jt_a(2)
                    !$omp atomic update
                    grad(3, it_a) = grad(3, it_a) + bds%strtork(3 + k, i) * dr2 * gt(k) * Jt_a(3)
                end if
                if (.not. sk_tb) then
                    !$omp atomic update
                    grad(1, it_b) = grad(1, it_b) + bds%strtork(3 + k, i) * dr2 * gt(k) * Jt_b(1)
                    !$omp atomic update
                    grad(2, it_b) = grad(2, it_b) + bds%strtork(3 + k, i) * dr2 * gt(k) * Jt_b(2)
                    !$omp atomic update
                    grad(3, it_b) = grad(3, it_b) + bds%strtork(3 + k, i) * dr2 * gt(k) * Jt_b(3)
                end if
                if (.not. sk_tc) then
                    !$omp atomic update
                    grad(1, it_c) = grad(1, it_c) + bds%strtork(3 + k, i) * dr2 * gt(k) * Jt_c(1)
                    !$omp atomic update
                    grad(2, it_c) = grad(2, it_c) + bds%strtork(3 + k, i) * dr2 * gt(k) * Jt_c(2)
                    !$omp atomic update
                    grad(3, it_c) = grad(3, it_c) + bds%strtork(3 + k, i) * dr2 * gt(k) * Jt_c(3)
                end if
                if (.not. sk_td) then
                    !$omp atomic update
                    grad(1, it_d) = grad(1, it_d) + bds%strtork(3 + k, i) * dr2 * gt(k) * Jt_d(1)
                    !$omp atomic update
                    grad(2, it_d) = grad(2, it_d) + bds%strtork(3 + k, i) * dr2 * gt(k) * Jt_d(2)
                    !$omp atomic update
                    grad(3, it_d) = grad(3, it_d) + bds%strtork(3 + k, i) * dr2 * gt(k) * Jt_d(3)
                end if
                if (.not. sk_2a) then
                    !$omp atomic update
                    grad(1, ib2_a) = grad(1, ib2_a) + bds%strtork(3 + k, i) * dihef(k) * Jb2_a(1)
                    !$omp atomic update
                    grad(2, ib2_a) = grad(2, ib2_a) + bds%strtork(3 + k, i) * dihef(k) * Jb2_a(2)
                    !$omp atomic update
                    grad(3, ib2_a) = grad(3, ib2_a) + bds%strtork(3 + k, i) * dihef(k) * Jb2_a(3)
                end if
                if (.not. sk_2b) then
                    !$omp atomic update
                    grad(1, ib2_b) = grad(1, ib2_b) + bds%strtork(3 + k, i) * dihef(k) * Jb2_b(1)
                    !$omp atomic update
                    grad(2, ib2_b) = grad(2, ib2_b) + bds%strtork(3 + k, i) * dihef(k) * Jb2_b(2)
                    !$omp atomic update
                    grad(3, ib2_b) = grad(3, ib2_b) + bds%strtork(3 + k, i) * dihef(k) * Jb2_b(3)
                end if
                if (.not. sk_ta) then
                    !$omp atomic update
                    grad(1, it_a) = grad(1, it_a) + bds%strtork(6 + k, i) * dr3 * gt(k) * Jt_a(1)
                    !$omp atomic update
                    grad(2, it_a) = grad(2, it_a) + bds%strtork(6 + k, i) * dr3 * gt(k) * Jt_a(2)
                    !$omp atomic update
                    grad(3, it_a) = grad(3, it_a) + bds%strtork(6 + k, i) * dr3 * gt(k) * Jt_a(3)
                end if
                if (.not. sk_tb) then
                    !$omp atomic update
                    grad(1, it_b) = grad(1, it_b) + bds%strtork(6 + k, i) * dr3 * gt(k) * Jt_b(1)
                    !$omp atomic update
                    grad(2, it_b) = grad(2, it_b) + bds%strtork(6 + k, i) * dr3 * gt(k) * Jt_b(2)
                    !$omp atomic update
                    grad(3, it_b) = grad(3, it_b) + bds%strtork(6 + k, i) * dr3 * gt(k) * Jt_b(3)
                end if
                if (.not. sk_tc) then
                    !$omp atomic update
                    grad(1, it_c) = grad(1, it_c) + bds%strtork(6 + k, i) * dr3 * gt(k) * Jt_c(1)
                    !$omp atomic update
                    grad(2, it_c) = grad(2, it_c) + bds%strtork(6 + k, i) * dr3 * gt(k) * Jt_c(2)
                    !$omp atomic update
                    grad(3, it_c) = grad(3, it_c) + bds%strtork(6 + k, i) * dr3 * gt(k) * Jt_c(3)
                end if
                if (.not. sk_td) then
                    !$omp atomic update
                    grad(1, it_d) = grad(1, it_d) + bds%strtork(6 + k, i) * dr3 * gt(k) * Jt_d(1)
                    !$omp atomic update
                    grad(2, it_d) = grad(2, it_d) + bds%strtork(6 + k, i) * dr3 * gt(k) * Jt_d(2)
                    !$omp atomic update
                    grad(3, it_d) = grad(3, it_d) + bds%strtork(6 + k, i) * dr3 * gt(k) * Jt_d(3)
                end if
                if (.not. sk_3a) then
                    !$omp atomic update
                    grad(1, ib3_a) = grad(1, ib3_a) + bds%strtork(6 + k, i) * dihef(k) * Jb3_a(1)
                    !$omp atomic update
                    grad(2, ib3_a) = grad(2, ib3_a) + bds%strtork(6 + k, i) * dihef(k) * Jb3_a(2)
                    !$omp atomic update
                    grad(3, ib3_a) = grad(3, ib3_a) + bds%strtork(6 + k, i) * dihef(k) * Jb3_a(3)
                end if
                if (.not. sk_3b) then
                    !$omp atomic update
                    grad(1, ib3_b) = grad(1, ib3_b) + bds%strtork(6 + k, i) * dihef(k) * Jb3_b(1)
                    !$omp atomic update
                    grad(2, ib3_b) = grad(2, ib3_b) + bds%strtork(6 + k, i) * dihef(k) * Jb3_b(2)
                    !$omp atomic update
                    grad(3, ib3_b) = grad(3, ib3_b) + bds%strtork(6 + k, i) * dihef(k) * Jb3_b(3)
                end if
            end do
        end do
    end subroutine strtor_geomgrad

    subroutine strtor_geomhess(bds, hess)
        !! Compute the Hessian of the stretch-torsion coupling term.
        !! Identical in spirit to angtor_geomhess (see there), but the
        !! internal coordinates are (theta, l1, l2, l3) instead of
        !! (theta, alpha1, alpha2): d^2U/dtheta^2 = f11 and
        !! d^2U/dtheta/dl_m = f1(1+m) are the only nonzero internal second
        !! derivatives (d^2U/dl_m^2 = d^2U/dl_m/dl_n = 0). For any two atoms
        !! X,Y of the term:
        !! \[ H_{XY} = f_{11}J^\theta_XJ^{\theta\dagger}_Y +
        !!    \sum_{m=1}^3 f_{1,1+m}\left(J^\theta_XJ^{l_m\dagger}_Y+J^{l_m}_XJ^{\theta\dagger}_Y\right)
        !!    + f_1H^\theta_{XY} + \sum_{m=1}^3 f_{1+m}H^{l_m}_{XY} \]
        !! with f1..f4 as in strtor_geomgrad, H^theta from
        !! torsion_angle_hessian, H^{l_m} from Rij_hessian on each bond.
        use mod_jacobian_mat, only: Rij_hessian, torsion_angle_hessian

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        !! Bonded potential data structure
        real(rp), intent(inout) :: hess(3,3,bds%top%mm_atoms,bds%top%mm_atoms)
        !! Hessian of stretch-torsion terms of potential energy

        real(rp) :: thet, gt(3), ht_(3), dihef(3), r1, r2, r3, dr(3)
        real(rp) :: f(4), fth(4)
        real(rp), dimension(3) :: Jt_a, Jt_b, Jt_c, Jt_d, &
                                  Jb1_a, Jb1_b, Jb2_a, Jb2_b, Jb3_a, Jb3_b
        real(rp), dimension(3,3) :: Haa, Hab, Hac, Had, Hbb, Hbc, Hbd, Hcc, Hcd, Hdd, &
                                    Hb1ii, Hb1ij, Hb1jj, Hb2ii, Hb2ij, Hb2jj, &
                                    Hb3ii, Hb3ij, Hb3jj, block_
        real(rp) :: Jslot(3,10), Mmat(4,4)
        real(rp) :: Ht(3,3,4,4), Hl(3,3,2,2,3)
        integer(ip) :: i, j, k, m, ib(3), it_a, it_b, it_c, it_d, p, q, tp, tq, ap, aq, &
                       iat(10)
        logical :: sk(10)

        if(.not. bds%use_strtor) return

        !$omp parallel do default(shared) schedule(dynamic) &
        !$omp private(i,j,k,m,ib,it_a,it_b,it_c,it_d,iat,sk,thet,gt,ht_,dihef,r1,r2,r3,dr) &
        !$omp private(f,fth,Jt_a,Jt_b,Jt_c,Jt_d,Jb1_a,Jb1_b,Jb2_a,Jb2_b,Jb3_a,Jb3_b) &
        !$omp private(Haa,Hab,Hac,Had,Hbb,Hbc,Hbd,Hcc,Hcd,Hdd) &
        !$omp private(Hb1ii,Hb1ij,Hb1jj,Hb2ii,Hb2ij,Hb2jj,Hb3ii,Hb3ij,Hb3jj) &
        !$omp private(Jslot,Mmat,Ht,Hl,p,q,tp,tq,ap,aq,block_)
        do i=1, bds%nstrtor
            it_a = bds%strtorat(1,i)
            it_b = bds%strtorat(2,i)
            it_c = bds%strtorat(3,i)
            it_d = bds%strtorat(4,i)
            iat(1)=it_a; iat(2)=it_b; iat(3)=it_c; iat(4)=it_d

            ib(1) = bds%strtor_b(1,i)
            ib(2) = bds%strtor_b(2,i)
            ib(3) = bds%strtor_b(3,i)
            iat(5) = bds%bondat(1,ib(1)); iat(6) = bds%bondat(2,ib(1))
            iat(7) = bds%bondat(1,ib(2)); iat(8) = bds%bondat(2,ib(2))
            iat(9) = bds%bondat(1,ib(3)); iat(10) = bds%bondat(2,ib(3))

            if(bds%top%use_frozen) then
                do k=1,10
                    sk(k) = bds%top%frozen(iat(k))
                end do
                if(all(sk)) cycle
            else
                sk = .false.
            end if

            call torsion_angle_hessian(bds%top%cmm(:,it_a), bds%top%cmm(:,it_b), &
                                       bds%top%cmm(:,it_c), bds%top%cmm(:,it_d), &
                                       thet, Jt_a, Jt_b, Jt_c, Jt_d, &
                                       Haa, Hab, Hac, Had, Hbb, Hbc, Hbd, Hcc, Hcd, Hdd)
            call Rij_hessian(bds%top%cmm(:,iat(5)), bds%top%cmm(:,iat(6)), &
                             r1, Jb1_a, Jb1_b, Hb1ii, Hb1ij, Hb1jj)
            call Rij_hessian(bds%top%cmm(:,iat(7)), bds%top%cmm(:,iat(8)), &
                             r2, Jb2_a, Jb2_b, Hb2ii, Hb2ij, Hb2jj)
            call Rij_hessian(bds%top%cmm(:,iat(9)), bds%top%cmm(:,iat(10)), &
                             r3, Jb3_a, Jb3_b, Hb3ii, Hb3ij, Hb3jj)

            do j=1,3
                gt(j) = -real(j,rp) * sin(j*thet + bds%torsphase(j,bds%strtor_t(i)))
                ht_(j) = -real(j,rp)**2 * cos(j*thet + bds%torsphase(j,bds%strtor_t(i)))
                dihef(j) = 1.0_rp + cos(j*thet + bds%torsphase(j,bds%strtor_t(i)))
            end do
            dr(1) = r1 - bds%l0bond(ib(1))
            dr(2) = r2 - bds%l0bond(ib(2))
            dr(3) = r3 - bds%l0bond(ib(3))

            f = 0.0_rp
            fth = 0.0_rp
            do m=1,3
                do j=1,3
                    f(1) = f(1) + bds%strtork((m-1)*3+j,i) * dr(m) * gt(j)
                    f(1+m) = f(1+m) + bds%strtork((m-1)*3+j,i) * dihef(j)
                    fth(1) = fth(1) + bds%strtork((m-1)*3+j,i) * dr(m) * ht_(j)
                    fth(1+m) = fth(1+m) + bds%strtork((m-1)*3+j,i) * gt(j)
                end do
            end do

            Jslot(:,1)=Jt_a; Jslot(:,2)=Jt_b; Jslot(:,3)=Jt_c; Jslot(:,4)=Jt_d
            Jslot(:,5)=Jb1_a; Jslot(:,6)=Jb1_b
            Jslot(:,7)=Jb2_a; Jslot(:,8)=Jb2_b
            Jslot(:,9)=Jb3_a; Jslot(:,10)=Jb3_b

            Mmat = 0.0_rp
            Mmat(1,1) = fth(1)
            Mmat(1,2) = fth(2); Mmat(2,1) = fth(2)
            Mmat(1,3) = fth(3); Mmat(3,1) = fth(3)
            Mmat(1,4) = fth(4); Mmat(4,1) = fth(4)

            Ht(:,:,1,1)=Haa; Ht(:,:,1,2)=Hab; Ht(:,:,1,3)=Hac; Ht(:,:,1,4)=Had
            Ht(:,:,2,1)=transpose(Hab); Ht(:,:,2,2)=Hbb; Ht(:,:,2,3)=Hbc; Ht(:,:,2,4)=Hbd
            Ht(:,:,3,1)=transpose(Hac); Ht(:,:,3,2)=transpose(Hbc); Ht(:,:,3,3)=Hcc; Ht(:,:,3,4)=Hcd
            Ht(:,:,4,1)=transpose(Had); Ht(:,:,4,2)=transpose(Hbd); Ht(:,:,4,3)=transpose(Hcd); Ht(:,:,4,4)=Hdd

            Hl(:,:,1,1,1)=Hb1ii; Hl(:,:,1,2,1)=Hb1ij; Hl(:,:,2,1,1)=transpose(Hb1ij); Hl(:,:,2,2,1)=Hb1jj
            Hl(:,:,1,1,2)=Hb2ii; Hl(:,:,1,2,2)=Hb2ij; Hl(:,:,2,1,2)=transpose(Hb2ij); Hl(:,:,2,2,2)=Hb2jj
            Hl(:,:,1,1,3)=Hb3ii; Hl(:,:,1,2,3)=Hb3ij; Hl(:,:,2,1,3)=transpose(Hb3ij); Hl(:,:,2,2,3)=Hb3jj

            !$omp critical
            do p=1,10
                if(sk(p)) cycle
                call slot_info(p, tp, ap)
                do q=1,10
                    if(sk(q)) cycle
                    call slot_info(q, tq, aq)
                    block_ = Mmat(tp,tq) * outer10(Jslot(:,p), Jslot(:,q))
                    if(tp == tq) then
                        if(tp == 1) then
                            block_ = block_ + f(1)*Ht(:,:,ap,aq)
                        else
                            block_ = block_ + f(tp)*Hl(:,:,ap,aq,tp-1)
                        end if
                    end if
                    hess(:,:,iat(p),iat(q)) = hess(:,:,iat(p),iat(q)) + block_
                end do
            end do
            !$omp end critical
        end do

    contains
        pure subroutine slot_info(k, tt, ll)
            !! slot->(type,local-index): type 1=theta(atoms1-4), 2=bond1(5-6),
            !! 3=bond2(7-8), 4=bond3(9-10)
            integer(ip), intent(in) :: k
            integer(ip), intent(out) :: tt, ll
            if(k <= 4) then
                tt = 1; ll = k
            else
                tt = 2 + (k-5)/2
                ll = mod(k-5,2) + 1
            end if
        end subroutine
        pure function outer10(u, v) result(m)
            real(rp), intent(in) :: u(3), v(3)
            real(rp) :: m(3,3)
            integer :: r
            do r=1,3
                m(r,:) = u(r)*v
            end do
        end function
    end subroutine strtor_geomhess

    
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
        
        call mallocate('tortor_init [tortorprm]', n, bds%tortorprm )
        call mallocate('tortor_init [tortorat]', 5_ip, n, bds%tortorat)

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
            call mallocate('torstors_newmap [itmp]', 2_ip, n_map, itmp)
            itmp = bds%ttmap_shape
            call mfree('torstors_newmap [ttmap_shape]', bds%ttmap_shape)
            call mallocate('torstors_newmap [ttmap_shape]', &
                           2_ip, n_map+1, bds%ttmap_shape)
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
            call mallocate('torstors_newmap [ttmap_shape]', 2_ip, 1_ip, bds%ttmap_shape)
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
        
        !$omp parallel do default(shared) reduction(+:V) &
        !$omp private(i,iprm,ibeg,j,iend,thetx,thety,vtt,dvttdx,dvttdy)
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

        !$omp parallel do default(shared) schedule(dynamic) &
        !$omp private(i,iprm,ibeg,j,iend,ia,ib,ic,id,ie,sk_a,sk_b,sk_c,sk_d,sk_e) &
        !$omp private(thetx,thety,J1_a,J1_b,J1_c,J1_d,J2_b,J2_c,J2_d,J2_e,vtt,dvttdx,dvttdy)
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


            if(.not. sk_a) then
                !$omp atomic update
                grad(1,ia) = grad(1,ia) + J1_a(1) * dvttdx
                !$omp atomic update
                grad(2,ia) = grad(2,ia) + J1_a(2) * dvttdx
                !$omp atomic update
                grad(3,ia) = grad(3,ia) + J1_a(3) * dvttdx
            end if
            if(.not. sk_b) then
                !$omp atomic update
                grad(1,ib) = grad(1,ib) + J1_b(1) * dvttdx + J2_b(1) * dvttdy
                !$omp atomic update
                grad(2,ib) = grad(2,ib) + J1_b(2) * dvttdx + J2_b(2) * dvttdy
                !$omp atomic update
                grad(3,ib) = grad(3,ib) + J1_b(3) * dvttdx + J2_b(3) * dvttdy
            end if
            if(.not. sk_c) then
                !$omp atomic update
                grad(1,ic) = grad(1,ic) + J1_c(1) * dvttdx + J2_c(1) * dvttdy
                !$omp atomic update
                grad(2,ic) = grad(2,ic) + J1_c(2) * dvttdx + J2_c(2) * dvttdy
                !$omp atomic update
                grad(3,ic) = grad(3,ic) + J1_c(3) * dvttdx + J2_c(3) * dvttdy
            end if
            if(.not. sk_d) then
                !$omp atomic update
                grad(1,id) = grad(1,id) + J1_d(1) * dvttdx + J2_d(1) * dvttdy
                !$omp atomic update
                grad(2,id) = grad(2,id) + J1_d(2) * dvttdx + J2_d(2) * dvttdy
                !$omp atomic update
                grad(3,id) = grad(3,id) + J1_d(3) * dvttdx + J2_d(3) * dvttdy
            end if
            if(.not. sk_e) then
                !$omp atomic update
                grad(1,ie) = grad(1,ie) + J2_e(1) * dvttdy
                !$omp atomic update
                grad(2,ie) = grad(2,ie) + J2_e(2) * dvttdy
                !$omp atomic update
                grad(3,ie) = grad(3,ie) + J2_e(3) * dvttdy
            end if
        end do

    end subroutine tortor_geomgrad

    subroutine tortor_geomhess(bds, hess)
        !! Compute the Hessian of the torsion-torsion (CMAP) potential.
        !! With phi1=torsion(A,B,C,D), phi2=torsion(B,C,D,E) (both signed, as
        !! in tortor_geomgrad) and the CMAP patch's second partials
        !! d2V/dphi1^2, d2V/dphi1/dphi2, d2V/dphi2^2 from
        !! compute_bicubic_interp_hess, for any two atoms X,Y of the term:
        !! \[ H_{XY} = \frac{\partial^2V}{\partial\varphi_1^2}J^{\varphi_1}_XJ^{\varphi_1\dagger}_Y
        !!    + \frac{\partial^2V}{\partial\varphi_1\partial\varphi_2}\left(J^{\varphi_1}_XJ^{\varphi_2\dagger}_Y+J^{\varphi_2}_XJ^{\varphi_1\dagger}_Y\right)
        !!    + \frac{\partial^2V}{\partial\varphi_2^2}J^{\varphi_2}_XJ^{\varphi_2\dagger}_Y
        !!    + \frac{\partial V}{\partial\varphi_1}H^{\varphi_1}_{XY} + \frac{\partial V}{\partial\varphi_2}H^{\varphi_2}_{XY} \]
        !! with H^phi1 (resp. H^phi2, both from torsion_angle_hessian) only
        !! contributing when X,Y are both among A,B,C,D (resp. B,C,D,E), and
        !! J^phi1=0 for atom E, J^phi2=0 for atom A.
        use mod_jacobian_mat, only: torsion_angle_hessian
        use mod_utils, only: compute_bicubic_interp_hess

        implicit none

        type(ommp_bonded_type), intent(in) :: bds
        !! Bonded potential data structure
        real(rp), intent(inout) :: hess(3,3,bds%top%mm_atoms,bds%top%mm_atoms)
        !! Hessian of torsion-torsion terms of potential energy

        real(rp) :: thetx, thety, vtt, dvttdx, dvttdy, d2vdx2, d2vdxdy, d2vdy2
        real(rp), dimension(3) :: J1_a, J1_b, J1_c, J1_d, J2_b, J2_c, J2_d, J2_e
        real(rp), dimension(3,3) :: H1aa, H1ab, H1ac, H1ad, H1bb, H1bc, H1bd, H1cc, H1cd, H1dd
        real(rp), dimension(3,3) :: H2bb, H2bc, H2bd, H2be, H2cc, H2cd, H2ce, H2dd, H2de, H2ee
        real(rp) :: J1f(3,5), J2f(3,5), H1f(3,3,5,5), H2f(3,3,5,5), block_(3,3)
        integer(ip) :: i, j, iprm, ibeg, iend, iat(5), p, q
        logical :: sk(5)

        if(.not. bds%use_tortor) return

        !$omp parallel do default(shared) schedule(dynamic) &
        !$omp private(i,iprm,ibeg,j,iend,iat,sk,thetx,thety,vtt,dvttdx,dvttdy) &
        !$omp private(d2vdx2,d2vdxdy,d2vdy2,J1_a,J1_b,J1_c,J1_d,J2_b,J2_c,J2_d,J2_e) &
        !$omp private(H1aa,H1ab,H1ac,H1ad,H1bb,H1bc,H1bd,H1cc,H1cd,H1dd) &
        !$omp private(H2bb,H2bc,H2bd,H2be,H2cc,H2cd,H2ce,H2dd,H2de,H2ee) &
        !$omp private(J1f,J2f,H1f,H2f,p,q,block_)
        do i=1, bds%ntortor
            iprm = bds%tortorprm(i)
            ibeg = 1
            do j=1, iprm-1
                ibeg = ibeg + bds%ttmap_shape(1,j)*bds%ttmap_shape(2,j)
            end do
            iend = ibeg + bds%ttmap_shape(1,iprm)*bds%ttmap_shape(2,iprm) - 1

            iat(1) = bds%tortorat(1,i)
            iat(2) = bds%tortorat(2,i)
            iat(3) = bds%tortorat(3,i)
            iat(4) = bds%tortorat(4,i)
            iat(5) = bds%tortorat(5,i)

            if(bds%top%use_frozen) then
                do j=1,5
                    sk(j) = bds%top%frozen(iat(j))
                end do
                if(all(sk)) cycle
            else
                sk = .false.
            end if

            call torsion_angle_hessian(bds%top%cmm(:,iat(1)), bds%top%cmm(:,iat(2)), &
                                       bds%top%cmm(:,iat(3)), bds%top%cmm(:,iat(4)), &
                                       thetx, J1_a, J1_b, J1_c, J1_d, &
                                       H1aa, H1ab, H1ac, H1ad, H1bb, H1bc, H1bd, H1cc, H1cd, H1dd)
            call torsion_angle_hessian(bds%top%cmm(:,iat(2)), bds%top%cmm(:,iat(3)), &
                                       bds%top%cmm(:,iat(4)), bds%top%cmm(:,iat(5)), &
                                       thety, J2_b, J2_c, J2_d, J2_e, &
                                       H2bb, H2bc, H2bd, H2be, H2cc, H2cd, H2ce, H2dd, H2de, H2ee)

            call compute_bicubic_interp_hess(thetx, thety, vtt, dvttdx, dvttdy, &
                                             d2vdx2, d2vdxdy, d2vdy2, &
                                             bds%ttmap_shape(1,iprm), bds%ttmap_shape(2,iprm), &
                                             bds%ttmap_ang1(ibeg:iend), bds%ttmap_ang2(ibeg:iend), &
                                             bds%ttmap_v(ibeg:iend), bds%ttmap_vx(ibeg:iend), &
                                             bds%ttmap_vy(ibeg:iend), bds%ttmap_vxy(ibeg:iend))

            J1f = 0.0_rp
            J1f(:,1)=J1_a; J1f(:,2)=J1_b; J1f(:,3)=J1_c; J1f(:,4)=J1_d
            J2f = 0.0_rp
            J2f(:,2)=J2_b; J2f(:,3)=J2_c; J2f(:,4)=J2_d; J2f(:,5)=J2_e

            H1f = 0.0_rp
            H1f(:,:,1,1)=H1aa; H1f(:,:,1,2)=H1ab; H1f(:,:,1,3)=H1ac; H1f(:,:,1,4)=H1ad
            H1f(:,:,2,1)=transpose(H1ab); H1f(:,:,2,2)=H1bb; H1f(:,:,2,3)=H1bc; H1f(:,:,2,4)=H1bd
            H1f(:,:,3,1)=transpose(H1ac); H1f(:,:,3,2)=transpose(H1bc); H1f(:,:,3,3)=H1cc; H1f(:,:,3,4)=H1cd
            H1f(:,:,4,1)=transpose(H1ad); H1f(:,:,4,2)=transpose(H1bd); H1f(:,:,4,3)=transpose(H1cd); H1f(:,:,4,4)=H1dd

            H2f = 0.0_rp
            H2f(:,:,2,2)=H2bb; H2f(:,:,2,3)=H2bc; H2f(:,:,2,4)=H2bd; H2f(:,:,2,5)=H2be
            H2f(:,:,3,2)=transpose(H2bc); H2f(:,:,3,3)=H2cc; H2f(:,:,3,4)=H2cd; H2f(:,:,3,5)=H2ce
            H2f(:,:,4,2)=transpose(H2bd); H2f(:,:,4,3)=transpose(H2cd); H2f(:,:,4,4)=H2dd; H2f(:,:,4,5)=H2de
            H2f(:,:,5,2)=transpose(H2be); H2f(:,:,5,3)=transpose(H2ce); H2f(:,:,5,4)=transpose(H2de); H2f(:,:,5,5)=H2ee

            !$omp critical
            do p=1,5
                if(sk(p)) cycle
                do q=1,5
                    if(sk(q)) cycle
                    block_ = d2vdx2*outer5(J1f(:,p),J1f(:,q)) &
                        + d2vdxdy*(outer5(J1f(:,p),J2f(:,q))+outer5(J2f(:,p),J1f(:,q))) &
                        + d2vdy2*outer5(J2f(:,p),J2f(:,q)) &
                        + dvttdx*H1f(:,:,p,q) + dvttdy*H2f(:,:,p,q)
                    hess(:,:,iat(p),iat(q)) = hess(:,:,iat(p),iat(q)) + block_
                end do
            end do
            !$omp end critical
        end do

    contains
        pure function outer5(u, v) result(m)
            real(rp), intent(in) :: u(3), v(3)
            real(rp) :: m(3,3)
            integer :: r
            do r=1,3
                m(r,:) = u(r)*v
            end do
        end function
    end subroutine tortor_geomhess

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
