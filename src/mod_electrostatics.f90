module mod_electrostatics
    use mod_memory, only: ip, rp

    implicit none 
    private

    public :: new_field_M2D

    contains

    subroutine q_V(q, dr, kernel, V)
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

    subroutine new_field_M2D(E)
        !! Compute the electric field of static multipoles at induced dipoles
        !! sites. This is only intended to be used to build the RHS of the 
        !! linear system. This field is modified by the indroduction of the 
        !! damped kernels and by the connectivity-based screening rules.

        use mod_mmpol, only: amoeba, polar_mm, mm_atoms, pol_atoms, n_ipd, q

        implicit none

        real(rp), intent(out) :: E(3, pol_atoms, n_ipd)
        !! Electric field

        integer(ip) :: i, j
        logical :: to_do_p, to_scale_p, to_do_d, to_scale_d, to_do, to_scale
        real(rp) :: kernel(5), dr(3), tmpE(3), scalf_p, scalf_d, scalf

        E = 0.0_rp

        if(amoeba) then
            do i=1, mm_atoms
                ! loop on sources
                do j=1, pol_atoms
                    if(polar_mm(j) == i) cycle
                    !loop on target
                    call screening_rules_M2D(i, j, 'P', to_do_p, to_scale_p, scalf_p)
                    call screening_rules_M2D(i, j, 'D', to_do_d, to_scale_d, scalf_d)
                    
                    if(to_do_p .or. to_do_d) then
                        call new_damped_coulomb_kernel(i, polar_mm(j), 3, kernel(1:4), dr)
                        
                        tmpE = 0.0_rp
                        call q_E(q(1,i), dr, kernel, tmpE)
                        call mu_E(q(2:4,i), dr, kernel, tmpE)
                        call quad_E(q(5:10,i), dr, kernel, tmpE)

                        if(to_do_p) then
                            if(to_scale_p) then
                                E(:, j, 1) = E(:, j, 1) + tmpE * scalf_p
                            else
                                E(:, j, 1) = E(:, j, 1) + tmpE
                            end if
                        end if

                        if(to_do_d) then
                            if(to_scale_d) then
                                E(:, j, 2) = E(:, j, 2) + tmpE * scalf_d
                            else
                                E(:, j, 2) = E(:, j, 2) + tmpE
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
                    call screening_rules_M2D(i, j, 'P', to_do, to_scale, scalf)
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

    subroutine screening_rules_M2D(i, j, field, tocompute, toscale, scalf)
        !! Utility function used to decide if the electric field of MM site i 
        !! should be computed for polarizable site j when creating the RHS of
        !! the linear system
        use mod_mmpol, only: amoeba, conn, pg_conn, polgrp_mmat, polar_mm, &
                             mmat_polgrp, pscale, dscale, fatal_error
        use mod_constants, only: eps_rp

        integer(ip), intent(in) :: i
        !! Index (in MM atoms list) of source site
        integer(ip), intent(in) :: j
        !! Index (in polarizable atoms list) of target site
        character, intent(in) :: field
        !! Which screening rules have to be applied? 'D' = screening rules
        !! for direct field; 'P' = screening rules for polarization field
        logical, intent(out) :: tocompute
        !! Interaction should be computed?
        logical, intent(out) :: toscale
        !! Interaction should be scaled?
        real(rp), intent(out) :: scalf
        !! Scale factor for the interaction

        integer(ip) :: ineigh, grp, igrp, pg_i
        integer(ip) :: j_mm

        ! Default return values
        tocompute = .true.
        toscale = .false.
        scalf = 1.0_rp
        
        j_mm = polar_mm(j)
        
        if((amoeba .and. field == 'P') .or. (.not. amoeba)) then
            ! P field, screening rules only based on connectivity
            ! this applies with very tiny differences on the calculation of
            ! scale factor to both WANG and AMOEBA FF
            do ineigh=1,4
                if(any(conn(ineigh)%ci(conn(ineigh)%ri(i): &
                                       conn(ineigh)%ri(i+1)-1) == j_mm)) then
                    ! This is a bit a wired way to compute the actual scale
                    ! factor for the interaction
                    if(ineigh == 3 .and. amoeba) then
                        if(mmat_polgrp(i) == mmat_polgrp(j_mm)) then
                            scalf = pscale(3) * pscale(5)
                        else
                            scalf = pscale(3)
                        end if
                    else
                        scalf = pscale(ineigh)
                    end if
                    
                    if(abs(scalf) < eps_rp) then 
                        ! Scaling factor is 0.0 so interaction should not be
                        ! computed
                        tocompute = .false.
                        return
                    else if(abs(scalf - 1.0_rp) < eps_rp) then
                        ! Scaling factor is 1.0 so it's just a normal 
                        ! interaction
                        return
                    else
                        toscale = .true.
                        return
                    end if
                end if
            end do
        else if(amoeba .and. field == 'D') then
            ! D field, this only applies to AMOEBA
            pg_i = mmat_polgrp(i)

            do ineigh=1,4
                do igrp=pg_conn(ineigh)%ri(pg_i), &
                        pg_conn(ineigh)%ri(pg_i+1)-1
                    ! Indexes of groups at distance ineigh from group pg_i
                    grp = pg_conn(ineigh)%ci(igrp)
                    
                    if(any(polgrp_mmat%ci(polgrp_mmat%ri(grp): &
                                          polgrp_mmat%ri(grp+1)-1) == j_mm)) then
                        ! If atom j is in a group at distance ineigh from the 
                        ! one of atom i, the field is scaled according to dscale
                        scalf = dscale(ineigh)
                        if(abs(scalf) < eps_rp) then 
                            ! Scaling factor is 0.0 so interaction should not be
                            ! computed
                            tocompute = .false.
                            return
                        else if(abs(scalf - 1.0_rp) < eps_rp) then
                            ! Scaling factor is 1.0 so it's just a normal 
                            ! interaction
                            return
                        else
                            toscale = .true.
                            return
                        end if
                    end if
                end do
            end do
        else 
            ! Unexpected error
            call fatal_error("Unexpected combination of parameter for screening_rules_M2D")
        end if
    end subroutine screening_rules_M2D

end module mod_electrostatics
