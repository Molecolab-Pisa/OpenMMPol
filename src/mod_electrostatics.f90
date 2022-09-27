module mod_electrostatics
    use mod_memory, only: ip, rp

    implicit none 
    private
    
    real(rp), allocatable :: V_M2M(:), E_M2M(:,:), Egrd_M2M(:,:)
    !! potential of MM permanent multipoles at MM sites; 
    !! shaped (ld_cart, mm_atoms).
    logical :: M2M_done
  
    real(rp), allocatable :: E_M2D(:,:,:)
    !! electric field of MM permanent multipoles at POL sites; 
    logical :: M2D_done
    

    public :: screening_rules, damped_coulomb_kernel, field_extD2D
    public :: energy_MM_MM, energy_MM_pol
    public :: electrostatics_terminate, prepare_M2M, prepare_M2D
    public :: e_m2d

    contains

    subroutine electrostatics_terminate()
        use mod_memory, only: mfree

        implicit none

        if(allocated(E_M2D)) &
            call mfree('electrostatics_terminate [E_M2D]', E_M2D)
        if(allocated(V_M2M)) &
            call mfree('electrostatics_terminate [V_M2M]', V_M2M)
        if(allocated(E_M2M)) &
            call mfree('electrostatics_terminate [E_M2M]', E_M2M)
        if(allocated(Egrd_M2M)) &
            call mfree('electrostatics_terminate [Egrd_M2M]', Egrd_M2M)

    end subroutine electrostatics_terminate

    subroutine energy_MM_MM(ene)
        !! This function computes the interaction energy of 
        !! static electric multipoles
        use mod_mmpol, only: amoeba, q
        use mod_memory, only: mallocate, mfree

        implicit none

        real(rp), intent(inout) :: ene
        !! Energy (results will be added)
        real(rp) :: eMM

        integer(ip) :: i

        call prepare_M2M()
        
        eMM = 0.0

        if(amoeba) then
            eMM = eMM + dot_product(q(1,:), V_M2M)
            do i=1, 3
                eMM = eMM - dot_product(q(i+1,:), E_M2M(i,:))
            end do

            do i=1,6
                if(i == 1 .or. i == 3 .or. i == 6) then
                    ! diagonal elements
                    eMM = eMM + dot_product(q(i+4,:), Egrd_M2M(i,:))
                else
                    ! off-diagonal elements (are stored once, but 
                    ! should be counted twice
                    eMM = eMM + 2.0 * dot_product(q(i+4,:), Egrd_M2M(i,:))
                end if
            end do
        else
            eMM = eMM + dot_product(q(1,:), V_M2M)
        end if

        ! Since potential is computed using all the sites each 
        ! interaction is counted twice
        ene = ene + 0.5_rp * eMM
    
    end subroutine energy_MM_MM
    
    subroutine energy_MM_pol(ene)
        !! This function computes the interaction energy of 
        !! static electric multipoles
        use mod_mmpol, only: amoeba, ipd
        use mod_memory, only: mallocate, mfree
        use mod_constants, only: OMMP_AMOEBA_D, OMMP_AMOEBA_P

        implicit none

        real(rp), intent(inout) :: ene
        !! Energy (results will be added)
        real(rp) :: eMM

        integer(ip) :: i

        call prepare_M2D()
        eMM = 0.0
        
        if(amoeba) then
            do i=1, 3
                eMM = eMM - dot_product(ipd(i,:,OMMP_AMOEBA_D), &
                                        E_M2D(i,:,OMMP_AMOEBA_P))
            end do
        else
            do i=1, 3
                eMM = eMM - dot_product(ipd(i,:,1), E_M2D(i,:,1))
            end do
        end if

        ene = ene + 0.5_rp * eMM

    end subroutine energy_MM_pol

    subroutine coulomb_kernel(dr, maxder, res)
        !! This function compute the coulomb kernel for the distance vector dr and 
        !! its derivative up to the value required by maxder.
        use mod_memory, only: ip, rp

        implicit none

        real(rp) :: dr(3)
        !! Distance vector from atom i to atom j
        integer(ip), intent(in) :: maxder
        !! Maximum derivative to be computed
        real(rp), intent(out), dimension(maxder+1) :: res
        !! Results vector
        
        integer(ip) :: ii
        real(rp) :: norm2_r, inv_norm_sq 

        if(maxder < 0) then
            ! Strange request, maybe a warning should be printed
            return 
        end if
        
        norm2_r = sqrt(dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3))
        res(1) = 1.0_rp / norm2_r
        inv_norm_sq = res(1) * res(1)

        do ii = 1, maxder
            res(ii+1) = res(ii) * inv_norm_sq
        end do
        
    end subroutine coulomb_kernel

    subroutine damped_coulomb_kernel(i, j, maxder, res, dr)
        !! This subroutine computes the damped coulomb kernel between two atoms.
        !! Note that this only makes sense between two MM atoms, as it is only used
        !! to compute the field that induces the point dipoles!

        use mod_mmpol, only: amoeba, thole, fatal_error, cmm
        use mod_memory, only: ip, rp
        use mod_constants, only: eps_rp
        
        implicit none

        integer(ip), intent(in) :: i, j
        !! Index of atoms (in MM atom list) for which the kernel should be
        !! computed
        integer(ip), intent(in) :: maxder
        !! Maximum derivative to be computed
        real(rp), intent(out), dimension(maxder+1) :: res
        !! Results vector
        real(rp), intent(out), dimension(3) :: dr
        !! Distance vector between i and j
        
        real(rp) :: s, u, u3, u4, fexp, eexp
        
        s = thole(i) * thole(j)
        
        ! Compute undamped kernels
        dr = cmm(:,j) - cmm(:,i)
        call coulomb_kernel(dr, maxder, res)

        if(abs(s) < eps_rp) then
            ! either thole(i) or thole(j) are zero, so the damped kernel
            ! is just the same as the regular one. Job is done.
            return
        end if

        u = 1.0_rp / (res(1) * s)
        u3 = u**3

        if(amoeba .and. u3 * 0.39_rp < 50.0_rp) then
            ! 1.0 / (res(1) * s)**3 * 0.39_rp < 50_rp TODO Remove 0.39, this is a
            ! FF constants, ask someone why this condition is here.
            ! Here basically we multiply the standard interactions kernel for
            ! dumping coefficients. The equations implemented here correspond to 
            ! eq. (5) of 10.1021/jp027815
            fexp = -0.39_rp * u3
            eexp = exp(fexp)
            if(maxder >= 1) res(2) = res(2) * (1.0_rp - eexp)
            if(maxder >= 2) res(3) = res(3) * (1.0_rp - (1.0_rp - fexp) * eexp)
            if(maxder >= 3) res(4) = res(4) * &
                            (1.0_rp - (1.0_rp - fexp + 0.6_rp * fexp * fexp) * eexp)
            if(maxder >= 4) res(5) = res(5) * &
                            (1.0_rp - (1.0_rp - fexp + 18.0_rp/35.0_rp * fexp *  fexp - &
                            9.0_rp/35.0_rp * fexp * fexp * fexp) * eexp)
            if(maxder >= 5) then
                call fatal_error("Damped Coulomb kernel (AMOEBA) only supports up to the 5th derivative")
            end if
        else if(.not. amoeba .and. res(1) > 1_rp/s) then
            ! TODO Again it is not clear to me why condition res(1) > 1_rp/s is here.
            u4 = u3*u
            if(maxder >= 1) res(2) = res(2) * (4.0_rp * u3 - 3.0_rp * u4)
            if(maxder >= 2) res(3) = res(3) * u4
            if(maxder >= 3) res(4) = res(4) * 0.2_rp * u4
            if(maxder >= 4) then
                call fatal_error("Damped Coulomb kernel (WANG) only supports up to the 4th derivative")
            end if
        end if
    end subroutine damped_coulomb_kernel

    subroutine q_elec_prop(q, dr, kernel, &
                           do_V, V, do_E, E, do_grdE, grdE, do_HE, HE)
        !! TODO
        !! Computes the electric potential of a charge \(q\) at position
        !! \(\mathbf{dr}\) from the charge itself. Pre-computed kernel should
        !! be provided as input. 
        !! The result is added to \(V\).
        !! $$ V_q = q \frac{f(r)}{||\mathbf{dr}||} $$
        use mod_mmpol, only: fatal_error

        implicit none

        real(rp), intent(in) :: q
        !! Charge
        real(rp), intent(in) :: dr(3)
        !! Distance vector
        real(rp), intent(in) :: kernel(:)
        !! Array of coulomb kernel (either damped or undamped)
        logical, intent(in) :: do_V, do_E, do_grdE, do_HE
        !! Flags to enable/disable calculation of different electrostatic 
        !! properties
        real(rp), intent(out) :: V, E(3), grdE(6), HE(9)
        !! Electric potential
        
        if(do_V) then
            V =  V + kernel(1) * q
        end if
        
        if(do_E) then
            E = E + q * kernel(2) * dr
        end if

        if(do_grdE) then
            ! xx
            grdE(1) = grdE(1) + q * 3.0_rp * kernel(3) * dr(1) * dr(1) - q * kernel(2)
            ! xy
            grdE(2) = grdE(2) + q * 3.0_rp * kernel(3) * dr(1) * dr(2)
            ! yy
            grdE(3) = grdE(3) + q * 3.0_rp * kernel(3) * dr(2) * dr(2) - q * kernel(2)
            ! xz
            grdE(4) = grdE(4) + q * 3.0_rp * kernel(3) * dr(1) * dr(3)
            ! yz
            grdE(5) = grdE(5) + q * 3.0_rp * kernel(3) * dr(2) * dr(3)
            ! zz
            grdE(6) = grdE(6) + q * 3.0_rp * kernel(3) * dr(3) * dr(3) - q * kernel(2)
        end if

        if(do_HE) then
            HE = 0.0
            call fatal_error("Field Hessian not implemented")
        end if

    end subroutine q_elec_prop
    
    subroutine mu_elec_prop(mu, dr, kernel, &
                            do_V, V, do_E, E, do_grdE, grdE, do_HE, HE)
        
        use mod_mmpol, only: fatal_error

        implicit none

        real(rp), intent(in) :: mu(3)
        !! point dipole
        real(rp), intent(in) :: dr(3)
        !! Distance vector
        real(rp), intent(in) :: kernel(:)
        !! Array of coulomb kernel (either damped or undamped)
        logical, intent(in) :: do_V, do_E, do_grdE, do_HE
        !! Flags to enable/disable calculation of different electrostatic 
        !! properties
        real(rp), intent(out) :: V, E(3), grdE(6), HE(9)
        !! Electric potential
        
        real(rp) :: mu_dot_dr
        
        mu_dot_dr = dot_product(mu, dr)

        if(do_V) then
            V = V + mu_dot_dr * kernel(2)
        end if
        
        if(do_E) then
            E = E + 3.0_rp * mu_dot_dr * dr * kernel(3) - mu * kernel(2)
        end if

        if(do_grdE) then
            ! xx
            grdE(1) = grdE(1) + mu_dot_dr * kernel(4) * 15.0 * dr(1) * dr(1) - &
                              (mu_dot_dr + 2.0 * mu(1) * dr(1)) * 3.0 * kernel(3)
            ! xy
            grdE(2) = grdE(2) + mu_dot_dr * kernel(4) * 15.0 * dr(1) * dr(2) - &
                              (mu(1)*dr(2) + mu(2)*dr(1)) * 3.0 * kernel(3)
            ! yy
            grdE(3) = grdE(3) + mu_dot_dr * kernel(4) * 15.0 * dr(2) * dr(2) - &
                              (mu_dot_dr + 2.0 * mu(2) * dr(2)) * 3.0 * kernel(3)
            ! xz
            grdE(4) = grdE(4) + mu_dot_dr * kernel(4) * 15.0 * dr(1) * dr(3) - &
                              (mu(1)*dr(3) + mu(3)*dr(1)) * 3.0 * kernel(3)
            ! yz
            grdE(5) = grdE(5) + mu_dot_dr * kernel(4) * 15.0 * dr(2) * dr(3) - &
                              (mu(2)*dr(3) + mu(3)*dr(2)) * 3.0 * kernel(3)
            ! zz
            grdE(6) = grdE(6) + mu_dot_dr * kernel(4) * 15.0 * dr(3) * dr(3) - &
                              (mu_dot_dr + 2.0 * mu(3) * dr(3)) * 3.0 * kernel(3)
        end if

        if(do_HE) then
            HE = 0.0
            call fatal_error("Field Hessian not implemented")
        end if
    end subroutine mu_elec_prop
    
    subroutine quad_elec_prop(quad, dr, kernel, &
                              do_V, V, do_E, E, do_grdE, grdE, do_HE, HE)
        
        use mod_mmpol, only: fatal_error

        implicit none

        real(rp), intent(in) :: quad(6)
        !! point quadrupole stored as (xx, xy, yy, xz, yz, zz)
        real(rp), intent(in) :: dr(3)
        !! Distance vector
        real(rp), intent(in) :: kernel(:)
        !! Array of coulomb kernel (either damped or undamped)
        logical, intent(in) :: do_V, do_E, do_grdE, do_HE
        !! Flags to enable/disable calculation of different electrostatic 
        !! properties
        real(rp), intent(out) :: V, E(3), grdE(6), HE(9)
        !! Electric potential
        
        real(rp) :: quadxr(3), quadxr_dot_r

        quadxr(1) = quad(1)*dr(1) + quad(2)*dr(2) + quad(4)*dr(3)
        quadxr(2) = quad(2)*dr(1) + quad(3)*dr(2) + quad(5)*dr(3)
        quadxr(3) = quad(4)*dr(1) + quad(5)*dr(2) + quad(6)*dr(3)

        quadxr_dot_r = dot_product(quadxr, dr)

        if(do_V) then
            V = V + 3.0_rp * quadxr_dot_r * kernel(3)
        end if
        
        if(do_E) then
            E = E + 15.0_rp * quadxr_dot_r * dr * kernel(4) &
                - 6.0_rp * quadxr * kernel(3)
        end if

        if(do_grdE) then
            ! xx
            grdE(1) = grdE(1) + 105.0 * kernel(5) * quadxr_dot_r * dr(1)*dr(1) + &
                              6.0 * kernel(3) * quad(1) - &
                              15.0*kernel(4)*(quadxr_dot_r+4.0*quadxr(1)*dr(1))
            ! xy
            grdE(2) = grdE(2) + 105.0 * kernel(5) * quadxr_dot_r * dr(1)*dr(2) + &
                              6.0 * kernel(3) * quad(2) - &
                              30.0*kernel(4)*(quadxr(1)*dr(2)+quadxr(2)*dr(1))
            ! yy
            grdE(3) = grdE(3) + 105.0 * kernel(5) * quadxr_dot_r * dr(2)*dr(2) + &
                              6.0 * kernel(3) * quad(3) - &
                              15.0*kernel(4)*(quadxr_dot_r+4.0*quadxr(2)*dr(2))
            ! xz
            grdE(4) = grdE(4) + 105.0 * kernel(5) * quadxr_dot_r * dr(1)*dr(3) + &
                              6.0 * kernel(3) * quad(4) - &
                              30.0*kernel(4)*(quadxr(1)*dr(3)+quadxr(3)*dr(1))
            ! yz
            grdE(5) = grdE(5) + 105.0 * kernel(5) * quadxr_dot_r * dr(2)*dr(3) + &
                              6.0 * kernel(3) * quad(5) - &
                              30.0*kernel(4)*(quadxr(2)*dr(3)+quadxr(3)*dr(2))
            ! zz
            grdE(6) = grdE(6) + 105.0 * kernel(5) * quadxr_dot_r * dr(3)*dr(3) + &
                              6.0 * kernel(3) * quad(6) - &
                              15.0*kernel(4)*(quadxr_dot_r+4.0*quadxr(3)*dr(3))
        end if

        if(do_HE) then
            HE = 0.0
            call fatal_error("Field Hessian not implemented")
        end if
    end subroutine quad_elec_prop

    subroutine prepare_M2M() 
        !! This function allocate and populate array of electrostatic 
        !! properties of static multipoles at static multipoles sites.
        !! It should be called blindly before any calculation that requires
        !! V_M2M etc.

        use mod_mmpol, only: mm_atoms, amoeba
        use mod_memory, only: mallocate

        implicit none

        if(M2M_done) return

        if(amoeba) then
            if(.not. allocated(V_M2M)) then
                call mallocate('prepare_m2m [V_M2M]', mm_atoms, V_M2M)
            end if

            if(.not. allocated(E_M2M)) then
                call mallocate('prepare_m2m [E_M2M]', 3, mm_atoms, E_M2M)
            end if

            if(.not. allocated(Egrd_M2M)) then
                call mallocate('prepare_m2m [Egrd_M2M]', 6, mm_atoms, Egrd_M2M)
            end if
            
            V_M2M = 0.0_rp
            E_M2M = 0.0_rp
            Egrd_M2M = 0.0_rp
            call potential_field_fieldgr_M2M(V_M2M, E_M2M, Egrd_M2M)
        else
            if(.not. allocated(V_M2M)) then
                call mallocate('prepare_m2m [V_M2M]', mm_atoms, V_M2M)
            end if
            
            V_M2M = 0.0_rp
            call potential_M2M(V_M2M)
        end if
        
        M2M_done = .true.

    end subroutine prepare_M2M

    subroutine prepare_M2D
        use mod_mmpol, only: n_ipd, pol_atoms
        use mod_memory, only: mallocate
        implicit none

        if(M2D_done) return

        if(.not. allocated(E_M2D)) then
            call mallocate('prepare_m2d [E_M2D]', 3, pol_atoms, n_ipd, E_M2D)
        end if

        E_M2D = 0.0_rp
        call field_M2D(E_M2D)

        M2D_done = .true.

    end subroutine prepare_M2D

  
    subroutine potential_field_fieldgr_M2M(V, E, Egr)
        !! Computes the electric potential, field and field gradients of 
        !! static multipoles at all sites (polarizable sites are a 
        !! subset of static ones)
        use mod_mmpol, only: amoeba, mm_atoms, cmm, q
        implicit none
        
        real(rp), intent(inout) :: V(mm_atoms)
        !! Potential on MM sites, results will be added
        real(rp), intent(inout) :: E(3, mm_atoms)
        !! Electric field on MM sites, order [x,y,z] results will be added
        real(rp), intent(inout) :: Egr(6, mm_atoms)
        !! Electric field gradients on MM sites, order [xx,xy,yy,xz,yz,zz]
        !! results will be added

        real(rp) :: kernel(5), dr(3), tmpV, tmpE(3), tmpEgr(6), tmpHE(9), scalf
        integer(ip) :: i, j
        logical :: to_do, to_scale

        if(amoeba) then
            !$omp parallel do private(to_do, to_scale, scalf, dr, kernel, tmpV, tmpE, tmpEgr) reduction(+:V, E, Egr)
            do i=1, mm_atoms
                ! loop on sources
                do j=1, mm_atoms
                    if(j == i) cycle
                    !loop on target
                    call screening_rules(i, 'S', j, 'S', '-', &
                                         to_do, to_scale, scalf)
                    
                    if(to_do) then
                        dr = cmm(:,j) - cmm(:, i)
                        call coulomb_kernel(dr, 4, kernel)
                        
                        tmpV = 0.0_rp
                        tmpE = 0.0_rp
                        tmpEgr = 0.0_rp

                        call q_elec_prop(q(1,i), dr, kernel, &
                                         .true., tmpV, & 
                                         .true., tmpE, &
                                         .true., tmpEgr, &
                                         .false., tmpHE)

                        call mu_elec_prop(q(2:4,i), dr, kernel, &
                                          .true., tmpV, & 
                                          .true., tmpE, &
                                          .true., tmpEgr, &
                                          .false., tmpHE)

                        call quad_elec_prop(q(5:10,i), dr, kernel, &
                                          .true., tmpV, & 
                                          .true., tmpE, &
                                          .true., tmpEgr, &
                                          .false., tmpHE)

                        if(to_scale) then
                            V(j) = V(j) + tmpV * scalf
                            E(:,j) = E(:,j) + tmpE * scalf
                            Egr(:,j) = Egr(:,j) + tmpEgr * scalf
                        else
                            V(j) = V(j) + tmpV
                            E(:,j) = E(:,j) + tmpE 
                            Egr(:,j) = Egr(:,j) + tmpEgr
                        end if

                    end if
                end do
            end do
            !$omp end parallel do
        else
            !$omp parallel do private(to_do, to_scale, scalf, dr, kernel, tmpV, tmpE, tmpEgr) reduction(+:V, E, Egr)
            do i=1, mm_atoms
                ! loop on sources
                do j=1, mm_atoms
                    if(j == i) cycle
                    !loop on target
                    call screening_rules(i, 'S', j, 'S', '-', &
                                         to_do, to_scale, scalf)
                    
                    if(to_do) then
                        dr = cmm(:,j) - cmm(:, i)
                        call coulomb_kernel(dr, 3, kernel)
                        
                        tmpV = 0.0_rp
                        tmpE = 0.0_rp
                        tmpEgr = 0.0_rp

                        call q_elec_prop(q(1,i), dr, kernel, &
                                         .true., tmpV, &
                                         .true., tmpE, &
                                         .true., tmpEgr, &
                                         .false., tmpHE)

                        if(to_scale) then
                            V(j) = V(j) + tmpV * scalf
                            E(:,j) = E(:,j) + tmpE * scalf
                            Egr(:,j) = Egr(:,j) + tmpEgr * scalf
                        else
                            V(j) = V(j) + tmpV
                            E(:,j) = E(:,j) + tmpE 
                            Egr(:,j) = Egr(:,j) + tmpEgr
                        end if

                    end if
                end do
            end do
            !$omp end parallel do
        end if
    end subroutine potential_field_fieldgr_M2M
    
    subroutine potential_M2M(V)
        !! Computes the electric potential of static multipoles at all sites
        !! (polarizable sites are a subset of static ones)
        use mod_mmpol, only: amoeba, mm_atoms, cmm, q
        implicit none
        
        real(rp), intent(inout) :: V(mm_atoms)
        !! Potential on MM sites, results will be added

        real(rp) :: kernel(3), dr(3), tmpV, tmpE(3), tmpEgr(6), tmpHE(9), scalf
        integer(ip) :: i, j
        logical :: to_do, to_scale

        if(amoeba) then
            do i=1, mm_atoms
                ! loop on sources
                do j=1, mm_atoms
                    if(j == i) cycle
                    !loop on target
                    call screening_rules(i, 'S', j, 'S', '-', &
                                         to_do, to_scale, scalf)
                    
                    if(to_do) then
                        dr = cmm(:,j) - cmm(:, i)
                        call coulomb_kernel(dr, 2, kernel)
                        
                        tmpV = 0.0_rp
                        call q_elec_prop(q(1,i), dr, kernel, .true., tmpV, &
                                         .false., tmpE, .false., tmpEgr, &
                                         .false., tmpHE)
                        
                        call mu_elec_prop(q(2:4,i), dr, kernel, .true., tmpV, &
                                          .false., tmpE, .false., tmpEgr, &
                                          .false., tmpHE)
                        
                        call quad_elec_prop(q(5:10,i), dr, kernel, .true., tmpV, &
                                          .false., tmpE, .false., tmpEgr, &
                                          .false., tmpHE)

                        if(to_scale) then
                            V(j) = V(j) + tmpV * scalf
                        else
                            V(j) = V(j) + tmpV
                        end if

                    end if
                end do
            end do
        else
            do i=1, mm_atoms
                ! loop on sources
                do j=1, mm_atoms
                    if(j == i) cycle
                    !loop on target
                    call screening_rules(i, 'S', j, 'S', '-', &
                                         to_do, to_scale, scalf)
                    
                    if(to_do) then
                        dr = cmm(:,j) - cmm(:, i)
                        call coulomb_kernel(dr, 0, kernel)
                        
                        tmpV = 0.0_rp
                        call q_elec_prop(q(1,i), dr, kernel, .true., tmpV, &
                                         .false., tmpE, .false., tmpEgr, &
                                         .false., tmpHE)

                        if(to_scale) then
                            V(j) = V(j) + tmpV * scalf
                        else
                            V(j) = V(j) + tmpV
                        end if

                    end if
                end do
            end do
        end if
        
    end subroutine potential_M2M

    subroutine field_extD2D(E, ext_ipd)
        !! Computes the electric field of a trial set of induced point dipoles
        !! at polarizable sites. This is intended to be used as matrix-vector
        !! routine in the solution of the linear system.
        
        use mod_mmpol, only: pol_atoms, polar_mm

        implicit none

        real(rp), intent(in) :: ext_ipd(3, pol_atoms)
        !! External induced point dipoles at polarizable sites
        real(rp), intent(inout) :: E(3, pol_atoms)
        !! Electric field (results will be added)

        integer(ip) :: i, j
        logical :: to_scale, to_do
        real(rp) :: kernel(5), dr(3), tmpV, tmpE(3), tmpEgr(6), tmpHE(9), scalf

        do i=1, pol_atoms
            do j=1, pol_atoms
                if(j == i) cycle
                !loop on target
                call screening_rules(i, 'P', j, 'P', '-', &
                                     to_do, to_scale, scalf)
                
                if(to_do) then
                    call damped_coulomb_kernel(polar_mm(i), polar_mm(j),& 
                                               2, kernel(1:3), dr)
                    
                    tmpE = 0.0_rp

                    call mu_elec_prop(ext_ipd(:,i), dr, kernel, .false., tmpV, &
                                      .true., tmpE, .false., tmpEgr, & 
                                      .false., tmpHE)
                    if(to_scale) then
                        E(:, j) = E(:, j) + tmpE * scalf
                    else
                        E(:, j) = E(:, j) + tmpE
                    end if
                end if
            end do
        end do
    end subroutine field_extD2D

    subroutine field_D2D(E)
        !! Computes the electric field of indued dipoles at induced dipoles
        !! sites.
        
        use mod_mmpol, only: amoeba, ipd, pol_atoms, n_ipd
        use mod_constants, only : OMMP_AMOEBA_P, OMMP_AMOEBA_D

        implicit none

        real(rp), intent(out) :: E(3, pol_atoms, n_ipd)
        !! Electric field (results will be added)

        if(amoeba) then
            call field_extD2D(E(:,:,OMMP_AMOEBA_P), ipd(:,:,OMMP_AMOEBA_P))
            call field_extD2D(E(:,:,OMMP_AMOEBA_D), ipd(:,:,OMMP_AMOEBA_D))
        else
            call field_extD2D(E(:,:,1), ipd(:,:,1))
        end if
    end subroutine field_D2D
    
    subroutine field_M2D(E)
        !! Computes the electric field of static multipoles at induced dipoles
        !! sites. This is only intended to be used to build the RHS of the 
        !! linear system. This field is modified by the indroduction of the 
        !! damped kernels and by the connectivity-based screening rules.

        use mod_mmpol, only: amoeba, polar_mm, mm_atoms, pol_atoms, n_ipd, q
        use mod_constants, only : OMMP_AMOEBA_P, OMMP_AMOEBA_D

        implicit none

        real(rp), intent(out) :: E(3, pol_atoms, n_ipd)
        !! Electric field (results will be added)

        integer(ip) :: i, j
        logical :: to_do_p, to_scale_p, to_do_d, to_scale_d, to_do, to_scale
        real(rp) :: kernel(5), dr(3), tmpV, tmpE(3), tmpEgr(6), tmpHE(9), &
                    scalf_p, scalf_d, scalf

        if(amoeba) then
            do i=1, mm_atoms
                ! loop on sources
                do j=1, pol_atoms
                    if(polar_mm(j) == i) cycle
                    !loop on target
                    call screening_rules(i, 'S', j, 'P', 'P', &
                                         to_do_p, to_scale_p, scalf_p)
                    call screening_rules(i, 'S', j, 'P', 'D', &
                                         to_do_d, to_scale_d, scalf_d)
                    
                    if(to_do_p .or. to_do_d) then
                        call damped_coulomb_kernel(i, polar_mm(j), 3, kernel(1:4), dr)
                        
                        tmpE = 0.0_rp
                        call q_elec_prop(q(1,i), dr, kernel, .false., tmpV, &
                                         .true., tmpE, .false., tmpEgr, & 
                                         .false., tmpHE)
                        call mu_elec_prop(q(2:4,i), dr, kernel, .false., tmpV, &
                                          .true., tmpE, .false., tmpEgr, & 
                                          .false., tmpHE)
                        call quad_elec_prop(q(5:10,i), dr, kernel, .false., tmpV, &
                                            .true., tmpE, .false., tmpEgr, & 
                                            .false., tmpHE)

                        if(to_do_p) then
                            if(to_scale_p) then
                                E(:, j, OMMP_AMOEBA_P) = E(:, j, OMMP_AMOEBA_P) + tmpE * scalf_p
                            else
                                E(:, j, OMMP_AMOEBA_P) = E(:, j, OMMP_AMOEBA_P) + tmpE
                            end if
                        end if

                        if(to_do_d) then
                            if(to_scale_d) then
                                E(:, j, OMMP_AMOEBA_D) = E(:, j, OMMP_AMOEBA_D) + tmpE * scalf_d
                            else
                                E(:, j, OMMP_AMOEBA_D) = E(:, j, OMMP_AMOEBA_D) + tmpE
                            end if
                        end if
                    end if
                end do
            end do
        else
            do i=1, mm_atoms
                ! loop on sources
                do j=1, pol_atoms
                    if(polar_mm(j) == i) cycle
                    !loop on target
                    call screening_rules(i, 'S', j, 'P', 'P', &
                                         to_do, to_scale, scalf)
                    if(to_do) then
                        call damped_coulomb_kernel(i, polar_mm(j), 1, kernel(1:2), dr)
                        
                        tmpE = 0.0_rp
                        call q_elec_prop(q(1,i), dr, kernel, .false., tmpV, &
                                         .true., tmpE, .false., tmpEgr, &
                                         .false., tmpHE)
                        if(to_scale) then
                            E(:, j, 1) = E(:, j, 1) + tmpE * scalf
                        else
                            E(:, j, 1) = E(:, j, 1) + tmpE
                        end if 
                    end if
                end do
            end do
        end if

    end subroutine field_M2D

    !! subroutine field_M2E(cpt, E)
    !!     !! This subroutine computes the electric field generated by the static
    !!     !! multipoles to a set of arbitrary coordinates, without applying
    !!     !! any screening rules.
    !!     
    !!     ! Test this function TODO
    !!     use mod_mmpol, only: amoeba, mm_atoms, pol_atoms, q, cmm
    !!     ! use mod_constants, only : OMMP_AMOEBA_P, OMMP_AMOEBA_D

    !!     implicit none

    !!     real(rp), intent(out) :: E(:,:)
    !!     !! Electric field (results will be added)
    !!     real(rp), intent(in) :: cpt(:,:)
    !!     !! Coordinates at which the electric field is requested

    !!     integer(ip) :: i, j, n_cpt
    !!     real(rp) :: kernel(5), dr(3), tmpV, tmpE(3), tmpEgr(6), tmpHE(9)

    !!     n_cpt = size(cpt, 2)

    !!     if(amoeba) then
    !!         do i=1, mm_atoms
    !!             do j=1, n_cpt
    !!                 dr = cpt(:,j) - cmm(:,i)
    !!                 call coulomb_kernel(dr, 3, kernel(1:4))
    !!                 tmpE = 0.0_rp
    !!                 
    !!                 call q_elec_prop(q(1,i), dr, kernel, .false., tmpV, &
    !!                                  .true., tmpE, .false., tmpEgr, & 
    !!                                  .false., tmpHE)
    !!                 call mu_elec_prop(q(2:4,i), dr, kernel, .false., tmpV, &
    !!                                   .true., tmpE, .false., tmpEgr, & 
    !!                                   .false., tmpHE)
    !!                 call quad_elec_prop(q(5:10,i), dr, kernel, .false., tmpV, &
    !!                                     .true., tmpE, .false., tmpEgr, & 
    !!                                     .false., tmpHE)

    !!                 E(:, j) = E(:, j) + tmpE
    !!             end do
    !!         end do
    !!     else
    !!         do i=1, mm_atoms
    !!             ! loop on sources
    !!             do j=1, pol_atoms
    !!                 dr = cpt(:,j) - cmm(:,i)
    !!                 call coulomb_kernel(dr, 3, kernel(1:4))
    !!                 tmpE = 0.0_rp
    !!                 
    !!                 call q_elec_prop(q(1,i), dr, kernel, .false., tmpV, &
    !!                                  .true., tmpE, .false., tmpEgr, & 
    !!                                  .false., tmpHE)

    !!                 E(:, j) = E(:, j) + tmpE
    !!             end do
    !!         end do
    !!     end if
    !! end subroutine field_M2E

    subroutine screening_rules(i, kind_i, j, kind_j, in_field, &
                               tocompute, toscale, scalf)
        !! Utility function used to decide if an interaction between sites i and j
        !! should be computed and eventually scaled by a factor.
        !! This function is intended to be used in \(\mathcalO(n^2)\) code, for
        !! linear scaling code lists should be built.
        !! This is written to minimize code repetitions, all the screening rules
        !! are handled in two possible cases: 
        !! 1. rules based on adjacency matrix
        !! 2. rules based on AMOEBA polarization groups

        use mod_mmpol, only: amoeba, conn, pg_conn, polgrp_mmat, polar_mm, &
                             mmat_polgrp, mscale, pscale, dscale, uscale, &
                             pscale_intra, fatal_error
        use mod_constants, only: eps_rp

        integer(ip), intent(in) :: i
        !! Index of source site (MM index is used for static sites, while
        !! Pol index is used for polarizable sites)
        integer(ip), intent(in) :: j
        !! Index of target site (MM index is used for static sites, while
        !! Pol index is used for polarizable sites)
        character, intent(in) :: in_field
        !! Which screening rules have to be applied? 'D' = screening rules
        !! for direct field; 'P' = screening rules for polarization field
        character, intent(in) :: kind_i, kind_j
        !! Type of sites i and j in the interaction for which the screening
        !! rules are required; possible choices are 'P' (polarizable site) or 
        !! 'S' (static site). Any other option will cause a fatal error.
        logical, intent(out) :: tocompute
        !! Interaction should be computed?
        logical, intent(out) :: toscale
        !! Interaction should be scaled?
        real(rp), intent(out) :: scalf
        !! Scale factor for the interaction

        integer(ip) :: ineigh, grp, igrp, pg_i
        integer(ip) :: j_mm, i_mm
        character :: field, interaction

        real(rp) :: myscale(4), myscale_intra(4)

        ! Decide which kind of rule should be used
        if(kind_i == 'P') then
            i_mm = polar_mm(i)
        else if(kind_i == 'S') then
            i_mm = i
        else
            call fatal_error('Unexpected value of kind_i in screening_rules')
            i_mm = 0
        end if

        if(kind_j == 'P') then
            j_mm = polar_mm(j)
        else if(kind_j == 'S') then
            j_mm = j
        else
            call fatal_error('Unexpected value of kind_j in screening_rules')
            j_mm = 0
        end if

        myscale = 0.0_rp
        if(kind_j == 'P' .and. kind_i == 'P') then
            ! Use IPD-IPD screening rules
            myscale = uscale
            interaction = 'P' !Pol-pol interaction is named P
        else if(kind_j == 'S' .and. kind_i == 'S') then
            ! Use static multipoles-static multipoles screening rules
            myscale = mscale
            field = '-'
            interaction = 'S' !Static-static interaction is named S
        else
            ! Use static multipoles-IPD screening rules
            if(in_field == 'P') then
                myscale = pscale
                if(amoeba) myscale_intra = pscale_intra
                field = 'P'
            else if(in_field == 'D' .and. amoeba) then
                myscale = dscale
                field = 'D'
            else
                call fatal_error('Unexpected value of field in screening rules')
            end if
            interaction = 'M' !Mixed interaction is named M
        end if

        ! Default return values
        tocompute = .true.
        toscale = .false.
        scalf = 1.0_rp
        
        if((.not. amoeba) .or. &
           (amoeba .and. interaction == 'M' .and. field == 'P') .or. &
           (amoeba .and. interaction == 'S')) then
            ! Screening rules based on connectivity matrix: 
            ! 1. all the ones of non-amoeba FF, 
            ! 2. the static to IPD polarization field
            ! 3. the static-static interaction of AMOEBA
            do ineigh=1,4
                ! Look if j is at distance ineigh from i
                if(any(conn(ineigh)%ci(conn(ineigh)%ri(i_mm): &
                                       conn(ineigh)%ri(i_mm+1)-1) == j_mm)) then
                   
                    if(interaction == 'M' .and. amoeba) then
                        ! Amoeba uses two different screening rules for atoms 
                        ! of the same group and for atoms of different group
                        if(mmat_polgrp(i_mm) == mmat_polgrp(j_mm)) then
                            scalf = myscale_intra(ineigh)
                        else
                            scalf = myscale(ineigh)
                        end if
                    else
                        ! In any other case just use the scale factor for this
                        ! distance
                        scalf = myscale(ineigh)
                    end if
                    
                    ! Exit the loop
                    exit 

                end if
            end do
        else if((amoeba .and. interaction == 'M' .and. field == 'D') .or. &
                (amoeba .and. interaction == 'P')) then
            ! Screening rules based on polarization groups: 
            ! 1. the static to IPD direct field
            ! 2. the IPD to IPD interaction of AMOEBA
            pg_i = mmat_polgrp(i_mm)

            outer: do ineigh=1,4
                do igrp=pg_conn(ineigh)%ri(pg_i), &
                        pg_conn(ineigh)%ri(pg_i+1)-1
                    ! Indexes of groups at distance ineigh from group pg_i
                    grp = pg_conn(ineigh)%ci(igrp)
                    
                    if(any(polgrp_mmat%ci(polgrp_mmat%ri(grp): &
                                          polgrp_mmat%ri(grp+1)-1) == j_mm)) then
                        ! If atom j is in a group at distance ineigh from the 
                        ! one of atom i, the field is scaled according to dscale
                        scalf = myscale(ineigh)
                        exit outer
                    end if
                end do
            end do outer
        else 
            ! Unexpected error
            call fatal_error("Unexpected combination of parameter for screening_rules")
        end if
                    
        if(abs(scalf) < eps_rp) then 
            ! Scaling factor is 0.0 so interaction should not be
            ! computed
            tocompute = .false.
        else if(abs(scalf - 1.0_rp) < eps_rp) then
            ! Scaling factor is 1.0 so it's just a normal 
            ! interaction
        else
            toscale = .true.
        end if

    end subroutine screening_rules

end module mod_electrostatics
