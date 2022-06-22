module mod_electrostatics
    use mod_memory, only: ip, rp

    implicit none 
    private

    public :: new_field_M2D, new_field_extD2D

    contains

    subroutine q_V(q, dr, kernel, V)
        !! Computes the electric potential of a charge \(q\) at position
        !! \(\mathbf{dr}\) from the charge itself. Pre-computed kernel should
        !! be provided as input. TODO
        !! The result is added to \(V\).
        !! $$ V_q = q \frac{f(r)}{||\mathbf{dr}||} $$
        implicit none

        real(rp), intent(in) :: q
        !! Charge
        real(rp), intent(in) :: dr(3)
        !! Distance vector
        real(rp), intent(in) :: kernel(:)
        !! Array of coulomb kernel (either damped or undamped)
        real(rp), intent(out) :: V
        !! Electric potential

        V =  V + kernel(1) * q
    end subroutine q_V

    subroutine q_E(q, dr, kernel, E)
        implicit none

        real(rp), intent(in) :: q
        !! Charge
        real(rp), intent(in) :: dr(3)
        !! Distance vector
        real(rp), intent(in) :: kernel(:)
        !! Array of coulomb kernel (either damped or undamped)
        real(rp), intent(out) :: E(3)
        !! Electric field

        E = E + q * kernel(2) * dr
    end subroutine q_E

    !subroutine q_grdE(...)
    !end subroutine q_grdE

    subroutine mu_V(mu, dr, kernel, V)
        implicit none

        real(rp), intent(in) :: mu(3)
        !! point dipole
        real(rp), intent(in) :: dr(3)
        !! Distance vector
        real(rp), intent(in) :: kernel(:)
        !! Array of coulomb kernel (either damped or undamped)
        real(rp), intent(out) :: V
        !! Electric potential (result will be added)

        V = V + dot_product(mu, dr) * kernel(2)
    end subroutine mu_V

    subroutine mu_E(mu, dr, kernel, E)
        implicit none

        real(rp), intent(in) :: mu(3)
        !! Point dipole
        real(rp), intent(in) :: dr(3)
        !! Distance vector
        real(rp), intent(in) :: kernel(:)
        !! Array of coulomb kernel (either damped or undamped)
        real(rp), intent(out) :: E(3)
        !! Electric potential (result will be added)

        E = E + 3.0_rp * dot_product(mu, dr) * dr * kernel(3) - mu * kernel(2)
    end subroutine mu_E

    !subroutine mu_grdE(...)
    !end subroutine mu_grdE

    subroutine quad_V(quad, dr, kernel, V)
        implicit none

        real(rp), intent(in) :: quad(6)
        !! point quadrupole stored as (xx, xy, yy, xz, yz, zz)
        real(rp), intent(in) :: dr(3)
        !! Distance vector
        real(rp), intent(in) :: kernel(:)
        !! Array of coulomb kernel (either damped or undamped)
        real(rp), intent(out) :: V
        !! Eectric potential

        real(rp) :: quadxr(3)

        quadxr(1) = quad(1)*dr(1) + quad(2)*dr(2) + quad(4)*dr(3)
        quadxr(2) = quad(2)*dr(1) + quad(3)*dr(2) + quad(5)*dr(3)
        quadxr(3) = quad(4)*dr(1) + quad(5)*dr(2) + quad(6)*dr(3)

        V = V + 3.0_rp * dot_product(quadxr, dr) * kernel(3)
    end subroutine quad_V

    subroutine quad_E(quad, dr, kernel, E)
        implicit none

        real(rp), intent(in) :: quad(6)
        !! point quadrupole stored as (xx, xy, yy, xz, yz, zz)
        real(rp), intent(in) :: dr(3)
        !! Distance vector
        real(rp), intent(in) :: kernel(:)
        !! Array of coulomb kernel (either damped or undamped)
        real(rp), intent(out) :: E(3)
        !! Electric potential

        real(rp) :: quadxr(3)

        quadxr(1) = quad(1)*dr(1) + quad(2)*dr(2) + quad(4)*dr(3)
        quadxr(2) = quad(2)*dr(1) + quad(3)*dr(2) + quad(5)*dr(3)
        quadxr(3) = quad(4)*dr(1) + quad(5)*dr(2) + quad(6)*dr(3)

        E = E + 15.0_rp * dot_product(quadxr, dr) * dr * kernel(4) &
            - 6.0_rp * quadxr * kernel(3)
    end subroutine quad_E

    !subroutine quad_grdE(...)
    !end subroutine quad_grdE
    
    subroutine new_field_extD2D(E, ext_ipd)
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
        real(rp) :: kernel(5), dr(3), tmpE(3), scalf

        do i=1, pol_atoms
            do j=1, pol_atoms
                if(j == i) cycle
                !loop on target
                call screening_rules(i, 'P', j, 'P', '-', &
                                     to_do, to_scale, scalf)
                
                if(to_do) then
                    call new_damped_coulomb_kernel(polar_mm(i), polar_mm(j),& 
                                                   2, kernel(1:3), dr)
                    
                    tmpE = 0.0_rp

                    call mu_E(ext_ipd(:,i), dr, kernel, tmpE)
                    if(to_scale) then
                        E(:, j) = E(:, j) + tmpE * scalf
                    else
                        E(:, j) = E(:, j) + tmpE
                    end if
                end if
            end do
        end do
    end subroutine new_field_extD2D

    subroutine new_field_D2D(E)
        !! Computes the electric field of indued dipoles at induced dipoles
        !! sites.
        
        use mod_mmpol, only: amoeba, ipd, pol_atoms, n_ipd
        use mod_constants, only : OMMP_AMOEBA_P, OMMP_AMOEBA_D

        implicit none

        real(rp), intent(out) :: E(3, pol_atoms, n_ipd)
        !! Electric field (results will be added)

        if(amoeba) then
            call new_field_extD2D(E(:,:,OMMP_AMOEBA_P), ipd(:,:,OMMP_AMOEBA_P))
            call new_field_extD2D(E(:,:,OMMP_AMOEBA_D), ipd(:,:,OMMP_AMOEBA_D))
        else
            call new_field_extD2D(E(:,:,1), ipd(:,:,1))
        end if
    end subroutine new_field_D2D
    
    subroutine new_field_M2D(E)
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
        real(rp) :: kernel(5), dr(3), tmpE(3), scalf_p, scalf_d, scalf

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
                        call new_damped_coulomb_kernel(i, polar_mm(j), 3, kernel(1:4), dr)
                        
                        tmpE = 0.0_rp
                        call q_E(q(1,i), dr, kernel, tmpE)
                        call mu_E(q(2:4,i), dr, kernel, tmpE)
                        call quad_E(q(5:10,i), dr, kernel, tmpE)

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
                        call new_damped_coulomb_kernel(i, polar_mm(j), 1, kernel(1:2), dr)
                        
                        tmpE = 0.0_rp
                        call q_E(q(1,i), dr, kernel, tmpE)
                        if(to_scale) then
                            E(:, j, 1) = E(:, j, 1) + tmpE * scalf
                        else
                            E(:, j, 1) = E(:, j, 1) + tmpE
                        end if 
                    end if
                end do
            end do
        end if

    end subroutine new_field_M2D

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
                             fatal_error
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

        real(rp) :: myscale(5)

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
            myscale(1:4) = uscale
            interaction = 'P' !Pol-pol interaction is named P
        else if(kind_j == 'S' .and. kind_i == 'S') then
            ! Use static multipoles-static multipoles screening rules
            myscale(1:4) = mscale
            field = '-'
            interaction = 'S' !Static-static interaction is named S
        else
            ! Use static multipoles-IPD screening rules
            if(in_field == 'P') then
                myscale(1:5) = pscale
                field = 'P'
            else if(in_field == 'D' .and. amoeba) then
                myscale(1:4) = dscale
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
                    
                    if(ineigh == 3 .and. amoeba .and. interaction == 'M') then
                        ! This is a bit a wired way to compute the actual scale
                        ! factor for the interaction in the case of AMOEBA 
                        ! polarization field that also takes into account
                        if(mmat_polgrp(i_mm) == mmat_polgrp(j_mm)) then
                            scalf = myscale(3) * myscale(5)
                        else
                            scalf = myscale(3)
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
