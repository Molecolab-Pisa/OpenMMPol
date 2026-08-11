#include "f_cart_components.h"

module mod_geomhess
    use mod_io, only: fatal_error, ommp_message
    use mod_memory, only: ip, rp
    use mod_mmpol, only: ommp_system
    use mod_topology, only: ommp_topology_type
    use mod_profiling, only: time_push, time_pull

    implicit none
    private

    public :: fixedelec_geomhess, polelec_geomhess, build_cpid_rhs, solve_cpid, &
              hess_pol_selfterm, hess_pol_pairterm, hess_pol_torque_pair

    contains

        subroutine solve_cpid(s, RHSd, RHSp, dmud, dmup)
            !! Solves the coupled-perturbed induced dipole (CPID) equations
            !! (eq. CPDD/CPDP): T*(dmu_d/dr) = RHSd, T*(dmu_p/dr) = RHSp,
            !! for every Cartesian perturbation of every mm_atom at once
            !! (nrhs = 3*mm_atoms), factoring T only once. RHSd/RHSp must
            !! have been computed beforehand by build_cpid_rhs (kept as a
            !! separate input here, rather than built internally, since
            !! polelec_geomhess also needs RHSd/RHSp themselves, not just
            !! the solution, to assemble the mixed term of eq. HessPol).

            use mod_polarization, only: create_TMat
            use mod_solvers, only: cp_inversion_solver
            use mod_memory, only: mallocate

            implicit none

            type(ommp_system), intent(inout), target :: s
            real(rp), dimension(3*s%eel%pol_atoms, 3*s%top%mm_atoms), intent(in) :: RHSd, RHSp
            real(rp), dimension(3*s%eel%pol_atoms, 3*s%top%mm_atoms), intent(out) :: dmud, dmup

            integer(ip) :: n, nrhs

            n = 3*s%eel%pol_atoms
            nrhs = 3*s%top%mm_atoms

            if(.not. allocated(s%eel%TMat)) call mallocate('solve_cpid [TMat]', n, n, s%eel%TMat)
            call create_TMat(s%eel)

            call cp_inversion_solver(n, nrhs, RHSd, dmud, s%eel%TMat)
            call cp_inversion_solver(n, nrhs, RHSp, dmup, s%eel%TMat)

        end subroutine solve_cpid

        subroutine add_dipole_ehes(hblock, sgn, mu, ehes)
            !! Adds sgn * Sum_gamma mu(gamma)*ehes(alpha,beta,gamma) into
            !! hblock(alpha,beta), upper triangle only (symmetric, fill
            !! the rest once after all contributions are accumulated).
            !! Same contraction pattern as hiiterm's dipole-EHes term.
            implicit none
            real(rp), intent(inout) :: hblock(3,3)
            real(rp), intent(in) :: sgn, mu(3), ehes(10)

            hblock(_x_,_x_) = hblock(_x_,_x_) + sgn*(mu(_x_)*ehes(_xxx_) + mu(_y_)*ehes(_xxy_) + mu(_z_)*ehes(_xxz_))
            hblock(_x_,_y_) = hblock(_x_,_y_) + sgn*(mu(_x_)*ehes(_xxy_) + mu(_y_)*ehes(_xyy_) + mu(_z_)*ehes(_xyz_))
            hblock(_x_,_z_) = hblock(_x_,_z_) + sgn*(mu(_x_)*ehes(_xxz_) + mu(_y_)*ehes(_xyz_) + mu(_z_)*ehes(_xzz_))
            hblock(_y_,_y_) = hblock(_y_,_y_) + sgn*(mu(_x_)*ehes(_xyy_) + mu(_y_)*ehes(_yyy_) + mu(_z_)*ehes(_yyz_))
            hblock(_y_,_z_) = hblock(_y_,_z_) + sgn*(mu(_x_)*ehes(_xyz_) + mu(_y_)*ehes(_yyz_) + mu(_z_)*ehes(_yzz_))
            hblock(_z_,_z_) = hblock(_z_,_z_) + sgn*(mu(_x_)*ehes(_xzz_) + mu(_y_)*ehes(_yzz_) + mu(_z_)*ehes(_zzz_))
        end subroutine add_dipole_ehes

        subroutine add_quad_e3d(hblock, sgn, qq, f3d)
            !! Adds sgn * Sum_{gamma,delta} qq(gamma,delta)*f3d(alpha,beta,
            !! gamma,delta) into hblock(alpha,beta), upper triangle only.
            !! Same contraction pattern as hiiterm's quadrupole-E3D term.
            implicit none
            real(rp), intent(inout) :: hblock(3,3)
            real(rp), intent(in) :: sgn, qq(6), f3d(15)

            hblock(_x_,_x_) = hblock(_x_,_x_) + sgn*( &
                f3d(_xxxx_)*qq(_xx_) + 2.0_rp*f3d(_xxxy_)*qq(_xy_) + 2.0_rp*f3d(_xxxz_)*qq(_xz_) &
              + f3d(_xxyy_)*qq(_yy_) + 2.0_rp*f3d(_xxyz_)*qq(_yz_) + f3d(_xxzz_)*qq(_zz_))
            hblock(_x_,_y_) = hblock(_x_,_y_) + sgn*( &
                f3d(_xyxx_)*qq(_xx_) + 2.0_rp*f3d(_xyxy_)*qq(_xy_) + 2.0_rp*f3d(_xyxz_)*qq(_xz_) &
              + f3d(_xyyy_)*qq(_yy_) + 2.0_rp*f3d(_xyyz_)*qq(_yz_) + f3d(_xyzz_)*qq(_zz_))
            hblock(_x_,_z_) = hblock(_x_,_z_) + sgn*( &
                f3d(_xzxx_)*qq(_xx_) + 2.0_rp*f3d(_xzxy_)*qq(_xy_) + 2.0_rp*f3d(_xzxz_)*qq(_xz_) &
              + f3d(_xzyy_)*qq(_yy_) + 2.0_rp*f3d(_xzyz_)*qq(_yz_) + f3d(_xzzz_)*qq(_zz_))
            hblock(_y_,_y_) = hblock(_y_,_y_) + sgn*( &
                f3d(_yyxx_)*qq(_xx_) + 2.0_rp*f3d(_yyxy_)*qq(_xy_) + 2.0_rp*f3d(_yyxz_)*qq(_xz_) &
              + f3d(_yyyy_)*qq(_yy_) + 2.0_rp*f3d(_yyyz_)*qq(_yz_) + f3d(_yyzz_)*qq(_zz_))
            hblock(_y_,_z_) = hblock(_y_,_z_) + sgn*( &
                f3d(_yzxx_)*qq(_xx_) + 2.0_rp*f3d(_yzxy_)*qq(_xy_) + 2.0_rp*f3d(_yzxz_)*qq(_xz_) &
              + f3d(_yzyy_)*qq(_yy_) + 2.0_rp*f3d(_yzyz_)*qq(_yz_) + f3d(_yzzz_)*qq(_zz_))
            hblock(_z_,_z_) = hblock(_z_,_z_) + sgn*( &
                f3d(_zzxx_)*qq(_xx_) + 2.0_rp*f3d(_zzxy_)*qq(_xy_) + 2.0_rp*f3d(_zzxz_)*qq(_xz_) &
              + f3d(_zzyy_)*qq(_yy_) + 2.0_rp*f3d(_zzyz_)*qq(_yz_) + f3d(_zzzz_)*qq(_zz_))
        end subroutine add_quad_e3d

        subroutine hess_pol_selfterm(s, k, hkk)
            !! Group 1 of HessPolExpl: the k=l "self" diagonal block,
            !! obtained by differentiating explicit_only_geomgrad(k) (the
            !! non-torque part of polelec_geomgrad) wrt r_k at fixed
            !! induced dipoles. explicit_only_geomgrad(k) has 4 pieces:
            !!   T1 = -q_k*E_D2M(k)
            !!   T2 = +dip_k.Egrd_D2M(k)
            !!   T3 = -quad_k:EHes_D2M(k)
            !!   T4 = 0.5*mud_k.[Egrd_M2D(k,P)+Egrd_D2D(k,P)]
            !!      + 0.5*mup_k.[Egrd_M2D(k,D)+Egrd_D2D(k,D)]
            !! using raw_n = -d(raw_{n-1})/dr_target (target=k throughout,
            !! since only k's own position varies here):
            !!   d(T1)/dr_k = +q_k*Egrd_D2M(k)                       [charge]
            !!   d(T2)/dr_k = -dip_k.EHes_D2M(k)                     [[c]+[d]]
            !!   d(T3)/dr_k = +quad_k:E3D_D2M(k)                     [[c]+[d]]
            !!   d(T4)/dr_k = -0.5*mud_k.EHes_M2D(k,P)                    [a]
            !!              - 0.5*mud_k.EHes_D2D(k,P)                    [new]
            !!              - 0.5*mup_k.EHes_M2D(k,D)                    [b]
            !!              - 0.5*mup_k.EHes_D2D(k,D)                    [e, was wrongly -1.0]
            !! EHes_D2M/E3D_D2M are pre-averaged over D/P sources inside
            !! prepare_polelec, so [c]+[d] are recomputed here from
            !! scratch (un-averaged, D-sourced/P-screened and
            !! P-sourced/D-screened separately) mirroring elec_prop_D2M's
            !! own pairwise construction, rather than reusing the shared
            !! (pre-averaged) eel%EHes_D2M/E3D_D2M state.
            !! EHes_M2D/EHes_D2D are already properly D/P-indexed (no
            !! averaging), so [a]/[b]/[new]/[e] reuse them directly --
            !! note EHes_D2D's D/P index labels the SOURCE dipole type
            !! (elec_prop_D2D uses eel%ipd(:,:,knd) as source, single
            !! list_P_P screening, no crossing), unlike M2D/D2M's
            !! source-independent D/P-as-screening-type convention.
            !!
            !! Sign convention throughout (established and validated
            !! against hiiterm/rotation_geomhess/build_cpid_rhs): a field
            !! variable X of rank n (E=1, Egrd=2, EHes=3, E3D=4), built by
            !! plain +accumulation of q/mu/quad_elec_prop outputs with
            !! dr=target-source, equals Phi^n = (-1)^n * X.

            use mod_electrostatics, only: ommp_electrostatics_type, damped_coulomb_kernel, &
                                          screening_rules, mu_elec_prop
            use mod_constants, only: eps_rp

            implicit none

            type(ommp_system), intent(inout), target :: s
            integer(ip), intent(in) :: k
            real(rp), intent(inout) :: hkk(3,3)

            type(ommp_electrostatics_type), pointer :: eel
            integer(ip) :: kpol, jpol, j
            real(rp) :: mud_k(3), mup_k(3), qqk(6)
            real(rp) :: EHesP(10), EHesD(10), E3DP(15), E3DD(15)
            real(rp) :: dr(3), kernel(6), scr_p, scr_d
            real(rp) :: tmpV, tmpE(3), tmpEgrd(6), tmpHE(10), tmpD3E(15)

            eel => s%eel
            kpol = eel%mm_polar(k)
            qqk(_xx_) = eel%q(4+_xx_,k); qqk(_xy_) = eel%q(4+_xy_,k); qqk(_yy_) = eel%q(4+_yy_,k)
            qqk(_xz_) = eel%q(4+_xz_,k); qqk(_yz_) = eel%q(4+_yz_,k); qqk(_zz_) = eel%q(4+_zz_,k)

            ! --- charge term: d(T1)/dr_k = +q_k*Egrd_D2M(k) ---
            hkk(_x_,_x_) = hkk(_x_,_x_) + eel%q(1,k) * eel%Egrd_D2M(_xx_,k)
            hkk(_x_,_y_) = hkk(_x_,_y_) + eel%q(1,k) * eel%Egrd_D2M(_xy_,k)
            hkk(_x_,_z_) = hkk(_x_,_z_) + eel%q(1,k) * eel%Egrd_D2M(_xz_,k)
            hkk(_y_,_y_) = hkk(_y_,_y_) + eel%q(1,k) * eel%Egrd_D2M(_yy_,k)
            hkk(_y_,_z_) = hkk(_y_,_z_) + eel%q(1,k) * eel%Egrd_D2M(_yz_,k)
            hkk(_z_,_z_) = hkk(_z_,_z_) + eel%q(1,k) * eel%Egrd_D2M(_zz_,k)

            ! --- (a),(new),(b),(e): reuse EHes_M2D / EHes_D2D ---
            if(kpol > 0) then
                mud_k = eel%ipd(:,kpol,_amoeba_D_)
                mup_k = eel%ipd(:,kpol,_amoeba_P_)
                call add_dipole_ehes(hkk, -0.5_rp, mud_k, eel%EHes_M2D(:,kpol,_amoeba_P_))
                call add_dipole_ehes(hkk, -0.5_rp, mud_k, eel%EHes_D2D(:,kpol,_amoeba_P_))
                call add_dipole_ehes(hkk, -0.5_rp, mup_k, eel%EHes_M2D(:,kpol,_amoeba_D_))
                call add_dipole_ehes(hkk, -0.5_rp, mup_k, eel%EHes_D2D(:,kpol,_amoeba_D_))
            end if

            ! --- (c),(d): fresh pairwise sums (un-averaged D2M direction) ---
            EHesP = 0.0_rp; EHesD = 0.0_rp; E3DP = 0.0_rp; E3DD = 0.0_rp
            do jpol = 1, eel%pol_atoms
                j = eel%polar_mm(jpol)
                if(j == k) cycle

                scr_p = screening_rules(eel, jpol, 'P', k, 'S', 'P')
                scr_d = screening_rules(eel, jpol, 'P', k, 'S', 'D')
                if(abs(scr_p) < eps_rp .and. abs(scr_d) < eps_rp) cycle

                call damped_coulomb_kernel(eel, j, k, 5_ip, kernel, dr)

                if(abs(scr_p) > eps_rp) then
                    ! mu_d source, P-screening -> feeds piece (c)
                    tmpV=0.0_rp; tmpE=0.0_rp; tmpEgrd=0.0_rp; tmpHE=0.0_rp; tmpD3E=0.0_rp
                    call mu_elec_prop(eel%ipd(:,jpol,_amoeba_D_), dr, kernel, &
                                      .false.,tmpV,.false.,tmpE,.false.,tmpEgrd,.true.,tmpHE,.true.,tmpD3E)
                    EHesP = EHesP + scr_p*tmpHE
                    E3DP = E3DP + scr_p*tmpD3E
                end if
                if(abs(scr_d) > eps_rp) then
                    ! mu_p source, D-screening -> feeds piece (d)
                    tmpV=0.0_rp; tmpE=0.0_rp; tmpEgrd=0.0_rp; tmpHE=0.0_rp; tmpD3E=0.0_rp
                    call mu_elec_prop(eel%ipd(:,jpol,_amoeba_P_), dr, kernel, &
                                      .false.,tmpV,.false.,tmpE,.false.,tmpEgrd,.true.,tmpHE,.true.,tmpD3E)
                    EHesD = EHesD + scr_d*tmpHE
                    E3DD = E3DD + scr_d*tmpD3E
                end if
            end do

            call add_dipole_ehes(hkk, -0.5_rp, eel%q(2:4,k), EHesP)
            call add_quad_e3d(hkk, 0.5_rp, qqk, E3DP)
            call add_dipole_ehes(hkk, -0.5_rp, eel%q(2:4,k), EHesD)
            call add_quad_e3d(hkk, 0.5_rp, qqk, E3DD)

            hkk(_y_,_x_) = hkk(_x_,_y_)
            hkk(_z_,_x_) = hkk(_x_,_z_)
            hkk(_z_,_y_) = hkk(_y_,_z_)

        end subroutine hess_pol_selfterm

        subroutine hess_pol_pairterm(s, k, l, hkl)
            !! Group 2 of HessPolExpl: the k/=l off-diagonal block,
            !! obtained by differentiating explicit_only_geomgrad(k) (see
            !! hess_pol_selfterm's docstring for its T1-T4 decomposition)
            !! wrt r_l (l/=k), at fixed induced dipoles. Only the j=l term
            !! of each sum survives; differentiating wrt the SOURCE side
            !! (r_l, not the target r_k) flips the sign of the raw_n/
            !! raw_{n+1} relation used in hess_pol_selfterm:
            !!   raw_{n+1} = +d(raw_n)/dr_source   (vs. -d/dr_target there)
            !! giving:
            !!   d(T1)/dr_l = -0.5*q_k*[s^p_kl*Egrd(mu_d_l) + s^d_kl*Egrd(mu_p_l)]
            !!   d(T2)/dr_l = +0.5*dip_k.[s^p_kl*EHes(mu_d_l) + s^d_kl*EHes(mu_p_l)]
            !!   d(T3)/dr_l = -0.5*quad_k:[s^p_kl*E3D(mu_d_l) + s^d_kl*E3D(mu_p_l)]
            !!   d(T4)/dr_l = +0.5*s^p_kl*mu_d_k.EHes(Theta_l)
            !!              + 0.5*s^u_kl*mu_d_k.EHes(mu_p_l)
            !!              + 0.5*s^d_kl*mu_p_k.EHes(Theta_l)
            !!              + 0.5*s^u_kl*mu_p_k.EHes(mu_d_l)
            !! where Egrd/EHes/E3D(mu_d_l) etc. are fresh pairwise sums
            !! (kernel(l,k), dr=r_k-r_l) exactly mirroring the source used
            !! for hess_pol_selfterm's (c)/(d) pieces, and EHes(Theta_l)
            !! is the same rank-3 field built from l's FULL permanent
            !! multipole (charge+dipole+quadrupole), mirroring
            !! elec_prop_M2D's own q_elec_prop+mu_elec_prop+quad_elec_prop
            !! construction (needed since G^p_k(Theta)/G^d_k(Theta), which
            !! T4's M2D piece reuses, has a generic multipole as source).
            !!
            !! d, u screening note (see hess_pol_selfterm and
            !! build_cpid_rhs for the general pattern): T1/T2/T3 and T4's
            !! M2D piece pair mu_d with s^p and mu_p with s^d throughout,
            !! exactly like hess_pol_selfterm's (a)-(d). T4's D2D piece is
            !! the one exception, screened with s^u (uscale/list_P_P) and
            !! pairing mu_d_k with mu_p_l and mu_p_k with mu_d_l (source
            !! and screening kind are decoupled there, since elec_prop_D2D
            !! always screens with list_P_P regardless of which dipole
            !! type is the source -- see hess_pol_selfterm's comment on
            !! EHes_D2D's D/P label meaning "source type", not "screening
            !! type").
            !!
            !! This block turns out to be symmetric in (kappa,lambda) at
            !! fixed (k,l): every piece above is built by contracting a
            !! fixed vector/tensor (dip_k, quad_k, mu_d_k, mu_p_k) into
            !! ONE slot of a fully index-symmetric rank-3/4 kernel tensor,
            !! leaving the other two slots (one from r_k, one from r_l)
            !! interchangeable -- a standard feature of two-body
            !! (pairwise, function-of-dr-only) interaction Hessians.

            use mod_electrostatics, only: ommp_electrostatics_type, damped_coulomb_kernel, &
                                          screening_rules, q_elec_prop, mu_elec_prop, quad_elec_prop
            use mod_constants, only: eps_rp

            implicit none

            type(ommp_system), intent(inout), target :: s
            integer(ip), intent(in) :: k, l
            real(rp), intent(inout) :: hkl(3,3)

            type(ommp_electrostatics_type), pointer :: eel
            integer(ip) :: kpol, lpol
            real(rp) :: dip_k(3), qqk(6), mud_k(3), mup_k(3)
            real(rp) :: dr(3), kernel(6)
            real(rp) :: scr_p_l, scr_d_l, scr_p_k, scr_d_k, scr_u
            real(rp) :: Egrd_d(6), EHes_d(10), E3D_d(15)
            real(rp) :: Egrd_p(6), EHes_p(10), E3D_p(15)
            real(rp) :: EHes_thl(10)
            real(rp) :: tmpV, tmpE(3), tmpEgrd(6), tmpD3E(15)

            eel => s%eel
            kpol = eel%mm_polar(k)
            lpol = eel%mm_polar(l)

            hkl = 0.0_rp
            if(kpol <= 0 .and. lpol <= 0) return

            scr_p_l = 0.0_rp; scr_d_l = 0.0_rp
            if(lpol > 0) then
                scr_p_l = screening_rules(eel, lpol, 'P', k, 'S', 'P')
                scr_d_l = screening_rules(eel, lpol, 'P', k, 'S', 'D')
            end if
            scr_p_k = 0.0_rp; scr_d_k = 0.0_rp
            if(kpol > 0) then
                scr_p_k = screening_rules(eel, kpol, 'P', l, 'S', 'P')
                scr_d_k = screening_rules(eel, kpol, 'P', l, 'S', 'D')
            end if
            scr_u = 0.0_rp
            if(kpol > 0 .and. lpol > 0) then
                scr_u = screening_rules(eel, kpol, 'P', lpol, 'P', '-')
            end if

            if(abs(scr_p_l) < eps_rp .and. abs(scr_d_l) < eps_rp .and. &
               abs(scr_p_k) < eps_rp .and. abs(scr_d_k) < eps_rp .and. &
               abs(scr_u) < eps_rp) return

            call damped_coulomb_kernel(eel, l, k, 5_ip, kernel, dr)

            if(lpol > 0) then
                tmpV=0.0_rp; tmpE=0.0_rp
                Egrd_d=0.0_rp; EHes_d=0.0_rp; E3D_d=0.0_rp
                call mu_elec_prop(eel%ipd(:,lpol,_amoeba_D_), dr, kernel, &
                                  .false.,tmpV,.false.,tmpE,.true.,Egrd_d,.true.,EHes_d,.true.,E3D_d)
                Egrd_p=0.0_rp; EHes_p=0.0_rp; E3D_p=0.0_rp
                call mu_elec_prop(eel%ipd(:,lpol,_amoeba_P_), dr, kernel, &
                                  .false.,tmpV,.false.,tmpE,.true.,Egrd_p,.true.,EHes_p,.true.,E3D_p)
            end if

            if(kpol > 0) then
                EHes_thl = 0.0_rp
                tmpV=0.0_rp; tmpE=0.0_rp; tmpEgrd=0.0_rp; tmpD3E=0.0_rp
                call q_elec_prop(eel%q(1,l), dr, kernel, &
                                 .false.,tmpV,.false.,tmpE,.false.,tmpEgrd,.true.,EHes_thl,.false.,tmpD3E)
                call mu_elec_prop(eel%q(2:4,l), dr, kernel, &
                                  .false.,tmpV,.false.,tmpE,.false.,tmpEgrd,.true.,EHes_thl,.false.,tmpD3E)
                call quad_elec_prop(eel%q(5:10,l), dr, kernel, &
                                    .false.,tmpV,.false.,tmpE,.false.,tmpEgrd,.true.,EHes_thl,.false.,tmpD3E)
            end if

            if(lpol > 0) then
                dip_k = eel%q(2:4,k)
                qqk(_xx_) = eel%q(4+_xx_,k); qqk(_xy_) = eel%q(4+_xy_,k); qqk(_yy_) = eel%q(4+_yy_,k)
                qqk(_xz_) = eel%q(4+_xz_,k); qqk(_yz_) = eel%q(4+_yz_,k); qqk(_zz_) = eel%q(4+_zz_,k)

                ! T1 (charge): direct rank-2 add
                hkl(_x_,_x_) = hkl(_x_,_x_) - 0.5_rp*eel%q(1,k)*(scr_p_l*Egrd_d(_xx_) + scr_d_l*Egrd_p(_xx_))
                hkl(_x_,_y_) = hkl(_x_,_y_) - 0.5_rp*eel%q(1,k)*(scr_p_l*Egrd_d(_xy_) + scr_d_l*Egrd_p(_xy_))
                hkl(_x_,_z_) = hkl(_x_,_z_) - 0.5_rp*eel%q(1,k)*(scr_p_l*Egrd_d(_xz_) + scr_d_l*Egrd_p(_xz_))
                hkl(_y_,_y_) = hkl(_y_,_y_) - 0.5_rp*eel%q(1,k)*(scr_p_l*Egrd_d(_yy_) + scr_d_l*Egrd_p(_yy_))
                hkl(_y_,_z_) = hkl(_y_,_z_) - 0.5_rp*eel%q(1,k)*(scr_p_l*Egrd_d(_yz_) + scr_d_l*Egrd_p(_yz_))
                hkl(_z_,_z_) = hkl(_z_,_z_) - 0.5_rp*eel%q(1,k)*(scr_p_l*Egrd_d(_zz_) + scr_d_l*Egrd_p(_zz_))

                ! T2 (dipole)
                call add_dipole_ehes(hkl, 0.5_rp*scr_p_l, dip_k, EHes_d)
                call add_dipole_ehes(hkl, 0.5_rp*scr_d_l, dip_k, EHes_p)

                ! T3 (quadrupole)
                call add_quad_e3d(hkl, -0.5_rp*scr_p_l, qqk, E3D_d)
                call add_quad_e3d(hkl, -0.5_rp*scr_d_l, qqk, E3D_p)
            end if

            if(kpol > 0) then
                mud_k = eel%ipd(:,kpol,_amoeba_D_)
                mup_k = eel%ipd(:,kpol,_amoeba_P_)

                ! T4, M2D piece
                call add_dipole_ehes(hkl, 0.5_rp*scr_p_k, mud_k, EHes_thl)
                call add_dipole_ehes(hkl, 0.5_rp*scr_d_k, mup_k, EHes_thl)

                ! T4, D2D piece (u-screened, source/screening decoupled)
                if(lpol > 0) then
                    call add_dipole_ehes(hkl, 0.5_rp*scr_u, mud_k, EHes_p)
                    call add_dipole_ehes(hkl, 0.5_rp*scr_u, mup_k, EHes_d)
                end if
            end if

            hkl(_y_,_x_) = hkl(_x_,_y_)
            hkl(_z_,_x_) = hkl(_x_,_z_)
            hkl(_z_,_y_) = hkl(_y_,_z_)

        end subroutine hess_pol_pairterm

        subroutine hess_pol_torque_pair(s, j, l, ddip, dqua, hess)
            !! Groups 3/4 of HessPolExpl, "cross" half: the piece of the
            !! torque contribution to grad(p) (p a frame atom of j, from
            !! rotation_geomgrad's -matmul(ddip_j,E_D2M(j)) +
            !! dqua_j:Egrd_D2M(j)) that comes from l/=j being one of the
            !! polarizable atoms feeding E_D2M(j)/Egrd_D2M(j)'s sum,
            !! rather than from j itself moving (that l=j case is already
            !! covered by rotation_geomhess's on-site call, reusing
            !! E_D2M/Egrd_D2M/EHes_D2M as-is -- see that subroutine's
            !! call in polelec_geomhess for why its sign matches T2/T3's
            !! own explicit-term sign exactly, one rank up).
            !!
            !! Differentiating the l=j case used a target-side derivative
            !! (raw_{n+1} = -d(raw_n)/dr_target, j being E_D2M(j)'s own
            !! target); here l is a SOURCE, so it's a source-side
            !! derivative instead (raw_{n+1} = +d(raw_n)/dr_source, same
            !! rule as hess_pol_pairterm), flipping the sign relative to
            !! the on-site piece:
            !!   d(-matmul(ddip_j,E_D2M(j)))/dr_l = -matmul(ddip_j, +Egrd(l))
            !!   d(+dqua_j:Egrd_D2M(j))/dr_l      = +dqua_j:(+EHes(l))
            !! landing in hess(:,:,p,l) and its transpose hess(:,:,l,p),
            !! for every frame atom p of j. mu_d/mu_p at l pair with s^p/
            !! s^d exactly as in hess_pol_selfterm/hess_pol_pairterm's
            !! (a)-(d) pieces (no s^u exception here -- this is a D2M-type
            !! field, not D2D).
            !!
            !! ddip/dqua must have been computed beforehand by a call to
            !! rotate_multipoles(eel, >=1, ddip, dqua, ...). Meant to be
            !! called once per ordered pair (j,l), j having an amoeba
            !! frame and l polarizable (screened, nonzero s^p or s^d).

            use mod_electrostatics, only: ommp_electrostatics_type, damped_coulomb_kernel, &
                                          screening_rules, mu_elec_prop
            use mod_rotate_multipoles, only: nref_atoms
            use mod_constants, only: eps_rp

            implicit none

            type(ommp_system), intent(inout), target :: s
            integer(ip), intent(in) :: j, l
            real(rp), intent(in) :: ddip(3,3,4,s%top%mm_atoms)
            real(rp), intent(in) :: dqua(3,3,3,4,s%top%mm_atoms)
            real(rp), dimension(3,3,s%top%mm_atoms,s%top%mm_atoms), intent(inout) :: hess

            type(ommp_electrostatics_type), pointer :: eel
            integer(ip) :: lpol, jat, dir, nact, atom(4)
            real(rp) :: scr_p, scr_d
            real(rp) :: dr(3), kernel(6)
            real(rp) :: Egrd_d(6), EHes_d(10), Egrd_p(6), EHes_p(10)
            real(rp) :: tmpV, tmpE(3), tmpD3E(15)
            real(rp) :: gmat_d(3,3), gmat_p(3,3), y(3,3)

            eel => s%eel
            lpol = eel%mm_polar(l)
            if(lpol <= 0) return

            nact = nref_atoms(eel%mol_frame(j))
            if(nact == 0) return

            ! the 0.5 here matches E_D2M/Egrd_D2M/EHes_D2M's own internal
            ! 0.5*(P-call+D-call) pre-averaging (see prepare_polelec),
            ! which rotation_geomgrad/rotation_geomhess's on-site reuse
            ! inherits automatically but this pairwise reconstruction of
            ! l's own contribution must apply explicitly (same pattern as
            ! hess_pol_selfterm/hess_pol_pairterm's own 0.5 factors).
            scr_p = 0.5_rp * screening_rules(eel, lpol, 'P', j, 'S', 'P')
            scr_d = 0.5_rp * screening_rules(eel, lpol, 'P', j, 'S', 'D')
            if(abs(scr_p) < eps_rp .and. abs(scr_d) < eps_rp) return

            atom(_self_) = j
            atom(_iz_) = eel%iz(j); if(atom(_iz_) == 0) atom(_iz_) = j
            atom(_ix_) = eel%ix(j); if(atom(_ix_) == 0) atom(_ix_) = j
            atom(_iy_) = eel%iy(j); if(atom(_iy_) == 0) atom(_iy_) = j

            call damped_coulomb_kernel(eel, l, j, 5_ip, kernel, dr)

            tmpV=0.0_rp; tmpE=0.0_rp; Egrd_d=0.0_rp; EHes_d=0.0_rp; tmpD3E=0.0_rp
            call mu_elec_prop(eel%ipd(:,lpol,_amoeba_D_), dr, kernel, &
                              .false.,tmpV,.false.,tmpE,.true.,Egrd_d,.true.,EHes_d,.false.,tmpD3E)
            tmpV=0.0_rp; tmpE=0.0_rp; Egrd_p=0.0_rp; EHes_p=0.0_rp; tmpD3E=0.0_rp
            call mu_elec_prop(eel%ipd(:,lpol,_amoeba_P_), dr, kernel, &
                              .false.,tmpV,.false.,tmpE,.true.,Egrd_p,.true.,EHes_p,.false.,tmpD3E)

            gmat_d(_x_,_x_)=Egrd_d(_xx_); gmat_d(_x_,_y_)=Egrd_d(_xy_); gmat_d(_x_,_z_)=Egrd_d(_xz_)
            gmat_d(_y_,_x_)=Egrd_d(_xy_); gmat_d(_y_,_y_)=Egrd_d(_yy_); gmat_d(_y_,_z_)=Egrd_d(_yz_)
            gmat_d(_z_,_x_)=Egrd_d(_xz_); gmat_d(_z_,_y_)=Egrd_d(_yz_); gmat_d(_z_,_z_)=Egrd_d(_zz_)
            gmat_p(_x_,_x_)=Egrd_p(_xx_); gmat_p(_x_,_y_)=Egrd_p(_xy_); gmat_p(_x_,_z_)=Egrd_p(_xz_)
            gmat_p(_y_,_x_)=Egrd_p(_xy_); gmat_p(_y_,_y_)=Egrd_p(_yy_); gmat_p(_y_,_z_)=Egrd_p(_yz_)
            gmat_p(_z_,_x_)=Egrd_p(_xz_); gmat_p(_z_,_y_)=Egrd_p(_yz_); gmat_p(_z_,_z_)=Egrd_p(_zz_)

            do jat = 1, nact
                if(eel%top%use_frozen) then
                    if(eel%top%frozen(atom(jat))) cycle
                end if

                ! dipole part: -matmul(ddip_j, Egrd(l))
                y = -scr_p*matmul(ddip(:,:,jat,j), gmat_d) - scr_d*matmul(ddip(:,:,jat,j), gmat_p)

                ! quadrupole part: +dqua_j:EHes(l)
                do dir = 1, 3
                    y(dir,_x_) = y(dir,_x_) &
                               + scr_p*(dqua(dir,_x_,_x_,jat,j) * EHes_d(_xxx_) + dqua(dir,_y_,_y_,jat,j) * EHes_d(_yyx_) &
                               +        dqua(dir,_z_,_z_,jat,j) * EHes_d(_zzx_) &
                               +        2.0_rp*(dqua(dir,_x_,_y_,jat,j) * EHes_d(_xyx_) + dqua(dir,_x_,_z_,jat,j) * EHes_d(_xzx_) &
                               +                dqua(dir,_y_,_z_,jat,j) * EHes_d(_yzx_))) &
                               + scr_d*(dqua(dir,_x_,_x_,jat,j) * EHes_p(_xxx_) + dqua(dir,_y_,_y_,jat,j) * EHes_p(_yyx_) &
                               +        dqua(dir,_z_,_z_,jat,j) * EHes_p(_zzx_) &
                               +        2.0_rp*(dqua(dir,_x_,_y_,jat,j) * EHes_p(_xyx_) + dqua(dir,_x_,_z_,jat,j) * EHes_p(_xzx_) &
                               +                dqua(dir,_y_,_z_,jat,j) * EHes_p(_yzx_)))
                    y(dir,_y_) = y(dir,_y_) &
                               + scr_p*(dqua(dir,_x_,_x_,jat,j) * EHes_d(_xxy_) + dqua(dir,_y_,_y_,jat,j) * EHes_d(_yyy_) &
                               +        dqua(dir,_z_,_z_,jat,j) * EHes_d(_zzy_) &
                               +        2.0_rp*(dqua(dir,_x_,_y_,jat,j) * EHes_d(_xyy_) + dqua(dir,_x_,_z_,jat,j) * EHes_d(_xzy_) &
                               +                dqua(dir,_y_,_z_,jat,j) * EHes_d(_yzy_))) &
                               + scr_d*(dqua(dir,_x_,_x_,jat,j) * EHes_p(_xxy_) + dqua(dir,_y_,_y_,jat,j) * EHes_p(_yyy_) &
                               +        dqua(dir,_z_,_z_,jat,j) * EHes_p(_zzy_) &
                               +        2.0_rp*(dqua(dir,_x_,_y_,jat,j) * EHes_p(_xyy_) + dqua(dir,_x_,_z_,jat,j) * EHes_p(_xzy_) &
                               +                dqua(dir,_y_,_z_,jat,j) * EHes_p(_yzy_)))
                    y(dir,_z_) = y(dir,_z_) &
                               + scr_p*(dqua(dir,_x_,_x_,jat,j) * EHes_d(_xxz_) + dqua(dir,_y_,_y_,jat,j) * EHes_d(_yyz_) &
                               +        dqua(dir,_z_,_z_,jat,j) * EHes_d(_zzz_) &
                               +        2.0_rp*(dqua(dir,_x_,_y_,jat,j) * EHes_d(_xyz_) + dqua(dir,_x_,_z_,jat,j) * EHes_d(_xzz_) &
                               +                dqua(dir,_y_,_z_,jat,j) * EHes_d(_yzz_))) &
                               + scr_d*(dqua(dir,_x_,_x_,jat,j) * EHes_p(_xxz_) + dqua(dir,_y_,_y_,jat,j) * EHes_p(_yyz_) &
                               +        dqua(dir,_z_,_z_,jat,j) * EHes_p(_zzz_) &
                               +        2.0_rp*(dqua(dir,_x_,_y_,jat,j) * EHes_p(_xyz_) + dqua(dir,_x_,_z_,jat,j) * EHes_p(_xzz_) &
                               +                dqua(dir,_y_,_z_,jat,j) * EHes_p(_yzz_)))
                end do

                hess(:,:,atom(jat),l) = hess(:,:,atom(jat),l) + y
                hess(:,:,l,atom(jat)) = hess(:,:,l,atom(jat)) + transpose(y)
            end do

        end subroutine hess_pol_torque_pair

        subroutine build_cpid_rhs(s, RHSd, RHSp)
            !! Builds the right-hand-sides of the coupled-perturbed induced
            !! dipole (CPID) equations (eq. CPDD/CPDP in main.tex), one
            !! column per Cartesian perturbation of every mm_atom:
            !!   RHSd(:, 3*(k-1)+beta) is RHS_d(k) with k perturbed along beta
            !!   RHSp(:, 3*(k-1)+beta) is RHS_p(k) with k perturbed along beta
            !! Each column has 3*pol_atoms rows, matching mu_d/mu_p and T.
            !!
            !! Requires prepare_polelec(eel, .true.) and up-to-date induced
            !! dipoles (eel%ipd) to have been computed beforehand (same
            !! preconditions as polelec_geomgrad's torque term).
            !!
            !! Three pieces per eq. CPDD (d shown; p is identical with
            !! d<->p and amoeba_D_<->amoeba_P_ swapped throughout):
            !!  A (k=i, "self"): reuses Egrd_M2D/Egrd_D2D as-is (these are
            !!     already exactly Gd_i(Theta)/Gu_i(mu_d) from the gradient
            !!     derivation, eq. dpol3).
            !!  B (k/=i, single pair): d_r_k of Ed_i's own pairwise term
            !!     with k, i.e. -Egrd(Theta_k, damped kernel) for the static
            !!     part and +Egrd(mu_d_k, damped kernel) for the "u" part
            !!     (sign flip because dr = r_i - r_k so d/dr_k = -d/dr).
            !!  C (torque): field at i generated by the frame-differentiated
            !!     multipoles of every j (d_k Theta_j), same rank as Ed_i
            !!     itself (E, not Egrd) since the frame-differentiated
            !!     multipole simply replaces Theta_j as a source.

            use mod_electrostatics, only: ommp_electrostatics_type, &
                                          damped_coulomb_kernel, screening_rules, &
                                          q_elec_prop, mu_elec_prop, quad_elec_prop
            use mod_rotate_multipoles, only: rotate_multipoles, nref_atoms
            use mod_constants, only: eps_rp

            implicit none

            type(ommp_system), intent(inout), target :: s
            real(rp), dimension(3*s%eel%pol_atoms, 3*s%top%mm_atoms), intent(out) :: RHSd, RHSp

            type(ommp_electrostatics_type), pointer :: eel
            real(rp), allocatable :: ddip(:,:,:,:), dqua(:,:,:,:,:)
            integer(ip) :: ipol, i, j, k, kpol, jat, nact, dir
            integer(ip) :: atomj(4)
            real(rp) :: gmat(3,3)
            real(rp) :: dr(3), kernel(5)
            real(rp) :: scalf_d, scalf_p, scalf_u
            real(rp) :: tmpV, tmpE(3), tmpEgrd(6), tmpHE(10), tmpD3E(15)

            eel => s%eel
            RHSd = 0.0_rp
            RHSp = 0.0_rp

            allocate(ddip(3,3,4,s%top%mm_atoms))
            allocate(dqua(3,3,3,4,s%top%mm_atoms))
            call rotate_multipoles(eel, 1_ip, ddip, dqua)

            do ipol = 1, eel%pol_atoms
                i = eel%polar_mm(ipol)

                ! --- Term A: k = i (self term), pure reuse ---
                gmat(_x_,_x_) = eel%Egrd_M2D(_xx_,ipol,_amoeba_D_) + eel%Egrd_D2D(_xx_,ipol,_amoeba_D_)
                gmat(_x_,_y_) = eel%Egrd_M2D(_xy_,ipol,_amoeba_D_) + eel%Egrd_D2D(_xy_,ipol,_amoeba_D_)
                gmat(_x_,_z_) = eel%Egrd_M2D(_xz_,ipol,_amoeba_D_) + eel%Egrd_D2D(_xz_,ipol,_amoeba_D_)
                gmat(_y_,_x_) = gmat(_x_,_y_)
                gmat(_y_,_y_) = eel%Egrd_M2D(_yy_,ipol,_amoeba_D_) + eel%Egrd_D2D(_yy_,ipol,_amoeba_D_)
                gmat(_y_,_z_) = eel%Egrd_M2D(_yz_,ipol,_amoeba_D_) + eel%Egrd_D2D(_yz_,ipol,_amoeba_D_)
                gmat(_z_,_x_) = gmat(_x_,_z_)
                gmat(_z_,_y_) = gmat(_y_,_z_)
                gmat(_z_,_z_) = eel%Egrd_M2D(_zz_,ipol,_amoeba_D_) + eel%Egrd_D2D(_zz_,ipol,_amoeba_D_)
                RHSd(3*(ipol-1)+1:3*ipol, 3*(i-1)+1:3*i) = RHSd(3*(ipol-1)+1:3*ipol, 3*(i-1)+1:3*i) - gmat

                gmat(_x_,_x_) = eel%Egrd_M2D(_xx_,ipol,_amoeba_P_) + eel%Egrd_D2D(_xx_,ipol,_amoeba_P_)
                gmat(_x_,_y_) = eel%Egrd_M2D(_xy_,ipol,_amoeba_P_) + eel%Egrd_D2D(_xy_,ipol,_amoeba_P_)
                gmat(_x_,_z_) = eel%Egrd_M2D(_xz_,ipol,_amoeba_P_) + eel%Egrd_D2D(_xz_,ipol,_amoeba_P_)
                gmat(_y_,_x_) = gmat(_x_,_y_)
                gmat(_y_,_y_) = eel%Egrd_M2D(_yy_,ipol,_amoeba_P_) + eel%Egrd_D2D(_yy_,ipol,_amoeba_P_)
                gmat(_y_,_z_) = eel%Egrd_M2D(_yz_,ipol,_amoeba_P_) + eel%Egrd_D2D(_yz_,ipol,_amoeba_P_)
                gmat(_z_,_x_) = gmat(_x_,_z_)
                gmat(_z_,_y_) = gmat(_y_,_z_)
                gmat(_z_,_z_) = eel%Egrd_M2D(_zz_,ipol,_amoeba_P_) + eel%Egrd_D2D(_zz_,ipol,_amoeba_P_)
                RHSp(3*(ipol-1)+1:3*ipol, 3*(i-1)+1:3*i) = RHSp(3*(ipol-1)+1:3*ipol, 3*(i-1)+1:3*i) - gmat

                ! --- Term B: k /= i, single-pair term ---
                do k = 1, s%top%mm_atoms
                    if(k == i) cycle

                    scalf_d = screening_rules(eel, k, 'S', ipol, 'P', 'D')
                    scalf_p = screening_rules(eel, k, 'S', ipol, 'P', 'P')
                    kpol = eel%mm_polar(k)
                    scalf_u = 0.0_rp
                    if(kpol > 0) scalf_u = screening_rules(eel, kpol, 'P', ipol, 'P', '-')

                    if(abs(scalf_d) < eps_rp .and. abs(scalf_p) < eps_rp .and. abs(scalf_u) < eps_rp) cycle

                    call damped_coulomb_kernel(eel, k, i, 4_ip, kernel, dr)

                    ! static (d/p) part: -d_r_k[Ed_i's pairwise term with k]
                    tmpV = 0.0_rp; tmpE = 0.0_rp; tmpEgrd = 0.0_rp; tmpHE = 0.0_rp; tmpD3E = 0.0_rp
                    call q_elec_prop(eel%q(1,k), dr, kernel, &
                                     .false., tmpV, .false., tmpE, .true., tmpEgrd, &
                                     .false., tmpHE, .false., tmpD3E)
                    call mu_elec_prop(eel%q(2:4,k), dr, kernel, &
                                      .false., tmpV, .false., tmpE, .true., tmpEgrd, &
                                      .false., tmpHE, .false., tmpD3E)
                    call quad_elec_prop(eel%q(5:10,k), dr, kernel, &
                                        .false., tmpV, .false., tmpE, .true., tmpEgrd, &
                                        .false., tmpHE, .false., tmpD3E)

                    gmat(_x_,_x_) = tmpEgrd(_xx_); gmat(_x_,_y_) = tmpEgrd(_xy_); gmat(_x_,_z_) = tmpEgrd(_xz_)
                    gmat(_y_,_x_) = tmpEgrd(_xy_); gmat(_y_,_y_) = tmpEgrd(_yy_); gmat(_y_,_z_) = tmpEgrd(_yz_)
                    gmat(_z_,_x_) = tmpEgrd(_xz_); gmat(_z_,_y_) = tmpEgrd(_yz_); gmat(_z_,_z_) = tmpEgrd(_zz_)

                    RHSd(3*(ipol-1)+1:3*ipol, 3*(k-1)+1:3*k) = RHSd(3*(ipol-1)+1:3*ipol, 3*(k-1)+1:3*k) &
                                                              + scalf_d * gmat
                    RHSp(3*(ipol-1)+1:3*ipol, 3*(k-1)+1:3*k) = RHSp(3*(ipol-1)+1:3*ipol, 3*(k-1)+1:3*k) &
                                                              + scalf_p * gmat

                    ! "u" part: +d_r_k[T_ik . mu_k]
                    if(kpol > 0 .and. abs(scalf_u) > eps_rp) then
                        tmpV = 0.0_rp; tmpE = 0.0_rp; tmpEgrd = 0.0_rp; tmpHE = 0.0_rp; tmpD3E = 0.0_rp
                        call mu_elec_prop(eel%ipd(:,kpol,_amoeba_D_), dr, kernel, &
                                          .false., tmpV, .false., tmpE, .true., tmpEgrd, &
                                          .false., tmpHE, .false., tmpD3E)
                        gmat(_x_,_x_) = tmpEgrd(_xx_); gmat(_x_,_y_) = tmpEgrd(_xy_); gmat(_x_,_z_) = tmpEgrd(_xz_)
                        gmat(_y_,_x_) = tmpEgrd(_xy_); gmat(_y_,_y_) = tmpEgrd(_yy_); gmat(_y_,_z_) = tmpEgrd(_yz_)
                        gmat(_z_,_x_) = tmpEgrd(_xz_); gmat(_z_,_y_) = tmpEgrd(_yz_); gmat(_z_,_z_) = tmpEgrd(_zz_)
                        RHSd(3*(ipol-1)+1:3*ipol, 3*(k-1)+1:3*k) = RHSd(3*(ipol-1)+1:3*ipol, 3*(k-1)+1:3*k) &
                                                                  + scalf_u * gmat

                        tmpV = 0.0_rp; tmpE = 0.0_rp; tmpEgrd = 0.0_rp; tmpHE = 0.0_rp; tmpD3E = 0.0_rp
                        call mu_elec_prop(eel%ipd(:,kpol,_amoeba_P_), dr, kernel, &
                                          .false., tmpV, .false., tmpE, .true., tmpEgrd, &
                                          .false., tmpHE, .false., tmpD3E)
                        gmat(_x_,_x_) = tmpEgrd(_xx_); gmat(_x_,_y_) = tmpEgrd(_xy_); gmat(_x_,_z_) = tmpEgrd(_xz_)
                        gmat(_y_,_x_) = tmpEgrd(_xy_); gmat(_y_,_y_) = tmpEgrd(_yy_); gmat(_y_,_z_) = tmpEgrd(_yz_)
                        gmat(_z_,_x_) = tmpEgrd(_xz_); gmat(_z_,_y_) = tmpEgrd(_yz_); gmat(_z_,_z_) = tmpEgrd(_zz_)
                        RHSp(3*(ipol-1)+1:3*ipol, 3*(k-1)+1:3*k) = RHSp(3*(ipol-1)+1:3*ipol, 3*(k-1)+1:3*k) &
                                                                  + scalf_u * gmat
                    end if
                end do

                ! --- Term C: torque, field of frame-differentiated multipoles ---
                do j = 1, s%top%mm_atoms
                    if(j == i) cycle

                    scalf_d = screening_rules(eel, j, 'S', ipol, 'P', 'D')
                    scalf_p = screening_rules(eel, j, 'S', ipol, 'P', 'P')
                    if(abs(scalf_d) < eps_rp .and. abs(scalf_p) < eps_rp) cycle

                    nact = nref_atoms(eel%mol_frame(j))
                    if(nact == 0) cycle

                    atomj(_self_) = j
                    atomj(_iz_) = eel%iz(j); if(atomj(_iz_) == 0) atomj(_iz_) = j
                    atomj(_ix_) = eel%ix(j); if(atomj(_ix_) == 0) atomj(_ix_) = j
                    atomj(_iy_) = eel%iy(j); if(atomj(_iy_) == 0) atomj(_iy_) = j

                    call damped_coulomb_kernel(eel, j, i, 4_ip, kernel, dr)

                    do jat = 1, nact
                        do dir = 1, 3
                            tmpV = 0.0_rp; tmpE = 0.0_rp; tmpEgrd = 0.0_rp; tmpHE = 0.0_rp; tmpD3E = 0.0_rp
                            call mu_elec_prop(ddip(dir,:,jat,j), dr, kernel, &
                                              .false., tmpV, .true., tmpE, .false., tmpEgrd, &
                                              .false., tmpHE, .false., tmpD3E)
                            call quad_elec_prop( &
                                (/dqua(dir,_x_,_x_,jat,j), dqua(dir,_x_,_y_,jat,j), dqua(dir,_y_,_y_,jat,j), &
                                  dqua(dir,_x_,_z_,jat,j), dqua(dir,_y_,_z_,jat,j), dqua(dir,_z_,_z_,jat,j)/), &
                                dr, kernel, &
                                .false., tmpV, .true., tmpE, .false., tmpEgrd, &
                                .false., tmpHE, .false., tmpD3E)

                            RHSd(3*(ipol-1)+1:3*ipol, 3*(atomj(jat)-1)+dir) = &
                                RHSd(3*(ipol-1)+1:3*ipol, 3*(atomj(jat)-1)+dir) + scalf_d * tmpE
                            RHSp(3*(ipol-1)+1:3*ipol, 3*(atomj(jat)-1)+dir) = &
                                RHSp(3*(ipol-1)+1:3*ipol, 3*(atomj(jat)-1)+dir) + scalf_p * tmpE
                        end do
                    end do
                end do
            end do

            deallocate(ddip)
            deallocate(dqua)

        end subroutine build_cpid_rhs

        subroutine fixedelec_geomhess(s, hess)
            use mod_electrostatics, only: prepare_fixedelec, &
                                          ommp_electrostatics_type
            use mod_rotate_multipoles, only: rotate_multipoles, rotation_geomhess, &
                                              rotation_geomhess_pair, rotation_geomhess_pair8

            implicit none

            type(ommp_system), intent(inout), target :: s
            !! System data structure
            real(rp), dimension(3,3,s%top%mm_atoms,s%top%mm_atoms), intent(inout) :: hess
            !! Geometrical Hessian in output, results will be added

            integer(ip) :: i, j, idx
            logical     :: to_do, to_scale
            real(rp)    :: scalf
            real(rp), allocatable :: ddip(:,:,:,:), dqua(:,:,:,:,:), &
                                      d2dip(:,:,:,:,:,:), d2qua(:,:,:,:,:,:,:)

            type(ommp_electrostatics_type), pointer :: eel
            eel => s%eel
            
            call time_push
            call prepare_fixedelec(eel, .false., .true.)
            call time_pull("Prepare fixedelec")

            call time_push
            if(eel%amoeba) then
                allocate(ddip(3,3,4,s%top%mm_atoms))
                allocate(dqua(3,3,3,4,s%top%mm_atoms))
                allocate(d2dip(3,3,3,4,4,s%top%mm_atoms))
                allocate(d2qua(3,3,3,3,4,4,s%top%mm_atoms))
                call rotate_multipoles(eel, 2_ip, ddip, dqua, d2dip, d2qua)

                !$omp parallel do
                do j=1, s%top%mm_atoms
                    ! If the atom is frozen, there are no contribution to compute
                    if(s%top%use_frozen) then
                        if(s%top%frozen(j)) cycle
                    end if
                    do i=1, s%top%mm_atoms
                        if(s%top%use_frozen) then
                            if(s%top%frozen(i)) cycle
                        end if
                        if (i.ne.j) then
                            to_do = .true.
                            to_scale = .false.
                            scalf = 1.0_rp

                            ! Check if the element should be scaled
                            do idx=eel%list_S_S%ri(i), eel%list_S_S%ri(i+1)-1
                                if(eel%list_S_S%ci(idx) == j) then
                                    to_scale = .true.
                                    exit
                                end if
                            end do

                            !If it should set the correct variables
                            if(to_scale) then
                                to_do = eel%todo_S_S(idx)
                                scalf = eel%scalef_S_S(idx)
                            end if

                            if(to_do) call hijterm(s,scalf,i,j,hess(:,:,i,j))
                        else
                            call hiiterm(s,i,hess(:,:,i,i))
                        end if
                    end do
                end do
                call time_push
                ! Torque contributions from multipoles rotation (on-site
                ! terms 3,4,7 of eq. Hess1: field of the plain multipoles
                ! contracted with the first/second derivative of the
                ! rotated multipoles themselves)
                call rotation_geomhess(eel, eel%E_M2M, eel%Egrd_M2M, eel%EHes_M2M, &
                                       ddip, dqua, d2dip, d2qua, hess)
                call time_pull("Rotation hess")

                call time_push
                ! Cross contributions (terms 5,6 of eq. Hess1: field
                ! generated at i by j's differentiated multipoles,
                ! contracted with i's plain multipole). Kept in its own
                ! serial loop: unlike the parallel loop above, the writes
                ! here land on hess(:,:,frame_atom_of_j,i) and its
                ! transpose, and frame_atom_of_j need not be j itself, so
                ! two different j's (handled by different threads) could
                ! alias the same hess column.
                do j=1, s%top%mm_atoms
                    if(s%top%use_frozen) then
                        if(s%top%frozen(j)) cycle
                    end if
                    do i=1, s%top%mm_atoms
                        if(s%top%use_frozen) then
                            if(s%top%frozen(i)) cycle
                        end if
                        if (i.ne.j) then
                            to_do = .true.
                            to_scale = .false.
                            scalf = 1.0_rp

                            do idx=eel%list_S_S%ri(i), eel%list_S_S%ri(i+1)-1
                                if(eel%list_S_S%ci(idx) == j) then
                                    to_scale = .true.
                                    exit
                                end if
                            end do

                            if(to_scale) then
                                to_do = eel%todo_S_S(idx)
                                scalf = eel%scalef_S_S(idx)
                            end if

                            if(to_do) then
                                call rotation_geomhess_pair(eel, scalf, i, j, ddip, dqua, hess)
                                call rotation_geomhess_pair8(eel, scalf, i, j, ddip, dqua, hess)
                            end if
                        end if
                    end do
                end do
                call time_pull("Rotation hess pair")

                deallocate(ddip)
                deallocate(dqua)
                deallocate(d2dip)
                deallocate(d2qua)
            else
                call fatal_error("WangAL switching has no second derivatives.")
            end if
            call time_pull("Hess sum")
        end subroutine

        subroutine polelec_geomhess(s, hess)
            !! Polarization contribution to the geometrical Hessian
            !! (eq. HessPol): the explicit second derivative of the
            !! Lagrangian (HessPolExpl, eq. 761-766) plus the mixed term
            !! that comes from the coupled-perturbed dipoles/Lagrange
            !! multipliers.
            !!
            !! The mixed term is [d2L/drk dmu_d] . [dmu_d/drl]
            !! + [d2L/drk dmu_p] . [dmu_p/drl] (eq. 743); since
            !! RHSd(i,k) = -2*(d2L/drk dmu_p)_i and RHSp(i,k) =
            !! -2*(d2L/drk dmu_d)_i (that's exactly how RHSd/RHSp were
            !! built -- they are minus twice the mixed derivative of the
            !! Lagrangian, by construction of the CPID equations), the
            !! mixed term is just -1/2*(RHSp^T . dmu_d + RHSd^T . dmu_p),
            !! reusing the already-validated RHSd/RHSp/dmud/dmup with no
            !! new kernel derivatives.

            use mod_electrostatics, only: prepare_polelec, ommp_electrostatics_type
            use mod_polarization, only: polarization
            use mod_rotate_multipoles, only: rotate_multipoles, rotation_geomhess

            implicit none

            type(ommp_system), intent(inout), target :: s
            !! System data structure
            real(rp), dimension(3,3,s%top%mm_atoms,s%top%mm_atoms), intent(inout) :: hess
            !! Geometrical Hessian in output, results will be added

            type(ommp_electrostatics_type), pointer :: eel
            real(rp), allocatable :: RHSd(:,:), RHSp(:,:), dmud(:,:), dmup(:,:), mixed(:,:)
            real(rp), allocatable :: ddip(:,:,:,:), dqua(:,:,:,:,:), d2dip(:,:,:,:,:,:), d2qua(:,:,:,:,:,:,:)
            integer(ip) :: n, npol, k, l, kappa, lambda

            eel => s%eel

            if(.not. eel%ipd_done) then
                call prepare_polelec(eel, .false.)
                call polarization(s, eel%e_M2D)
            end if
            call prepare_polelec(eel, .true., .true.)

            n = 3*s%top%mm_atoms
            npol = 3*eel%pol_atoms

            allocate(RHSd(npol,n), RHSp(npol,n))
            call build_cpid_rhs(s, RHSd, RHSp)

            allocate(dmud(npol,n), dmup(npol,n))
            call solve_cpid(s, RHSd, RHSp, dmud, dmup)

            allocate(mixed(n,n))
            mixed = -0.5_rp * (matmul(transpose(RHSp), dmud) + matmul(transpose(RHSd), dmup))

            do l = 1, s%top%mm_atoms
                do k = 1, s%top%mm_atoms
                    do lambda = 1, 3
                        do kappa = 1, 3
                            hess(kappa,lambda,k,l) = hess(kappa,lambda,k,l) &
                                                    + mixed(3*(k-1)+kappa, 3*(l-1)+lambda)
                        end do
                    end do
                end do
            end do

            ! Group 1 of HessPolExpl: k=l self term
            block
                real(rp) :: hkk(3,3)
                do k = 1, s%top%mm_atoms
                    hkk = 0.0_rp
                    call hess_pol_selfterm(s, k, hkk)
                    hess(:,:,k,k) = hess(:,:,k,k) + hkk
                end do
            end block

            ! Group 2 of HessPolExpl: k/=l off-diagonal pairwise term
            block
                real(rp) :: hkl(3,3)
                do l = 1, s%top%mm_atoms
                    do k = 1, s%top%mm_atoms
                        if(k == l) cycle
                        hkl = 0.0_rp
                        call hess_pol_pairterm(s, k, l, hkl)
                        hess(:,:,k,l) = hess(:,:,k,l) + hkl
                    end do
                end do
            end block

            ! Groups 3/4/5 of HessPolExpl (on-site torque piece): the
            ! same rotation_geomhess machinery used by fixedelec_geomhess
            ! for the M2M field, fed with the D2M field (from induced
            ! dipoles) instead -- (d_p Theta_j).Egrd_D2M(j) -
            ! (d_p qua_j):EHes_D2M(j), landing in hess(:,:,p,j)/(j,p) for
            ! p in j's frame (single torque), and (d_p d_q Theta_j) field
            ! contractions landing in hess(:,:,p,q) (double torque).
            if(eel%amoeba) then
                allocate(ddip(3,3,4,s%top%mm_atoms))
                allocate(dqua(3,3,3,4,s%top%mm_atoms))
                allocate(d2dip(3,3,3,4,4,s%top%mm_atoms))
                allocate(d2qua(3,3,3,3,4,4,s%top%mm_atoms))
                call rotate_multipoles(eel, 2_ip, ddip, dqua, d2dip, d2qua)
                call rotation_geomhess(eel, eel%E_D2M, eel%Egrd_D2M, eel%EHes_D2M, &
                                       ddip, dqua, d2dip, d2qua, hess)

                ! Groups 3/4 "cross" half: l/=j feeding E_D2M(j)/Egrd_D2M(j)
                do l = 1, s%top%mm_atoms
                    do k = 1, s%top%mm_atoms
                        if(k == l) cycle
                        call hess_pol_torque_pair(s, k, l, ddip, dqua, hess)
                    end do
                end do

                deallocate(ddip)
                deallocate(dqua)
                deallocate(d2dip)
                deallocate(d2qua)
            end if

            deallocate(RHSd)
            deallocate(RHSp)
            deallocate(dmud)
            deallocate(dmup)
            deallocate(mixed)

        end subroutine polelec_geomhess
!
        subroutine hiiterm(s,i,hii)
            use mod_electrostatics
            implicit none
            type(ommp_system), intent(inout), target         :: s
            integer(ip),       intent(in)                    :: i
            real(rp),          intent(inout), dimension(3,3) :: hii
!
            real(rp)     :: qi, di(3), qqi(6), fg(6), fh(10), f3d(15)
!
            hii = 0.0_rp
            qi  = s%eel%q(1,i)
            di  = s%eel%q(2:4,i)
            qqi = s%eel%q(5:10,i)
            fg  = s%eel%Egrd_M2M(1:6,i)
            fh  = s%eel%EHes_M2M(1:10,i)
            f3d = s%eel%E3D_M2M(1:15,i)

            hii(_x_,_x_) = hii(_x_,_x_) + fg(_xx_) * qi &
                - fh(_xxx_)*di(_x_) - fh(_xxy_)*di(_y_) - fh(_xxz_)*di(_z_) &
                + f3d(_xxxx_)*qqi(_xx_) + 2.0_rp*f3d(_xxxy_)*qqi(_xy_) &
                + 2.0_rp*f3d(_xxxz_)*qqi(_xz_) + f3d(_xxyy_)*qqi(_yy_) &
                + 2.0_rp*f3d(_xxyz_)*qqi(_yz_) + f3d(_xxzz_)*qqi(_zz_)

            hii(_x_,_y_) = hii(_x_,_y_) + fg(_xy_) * qi &
                - fh(_xyx_)*di(_x_) - fh(_xyy_)*di(_y_) - fh(_xyz_)*di(_z_) &
                + f3d(_xyxx_)*qqi(_xx_) + 2.0_rp*f3d(_xyxy_)*qqi(_xy_) &
                + 2.0_rp*f3d(_xyxz_)*qqi(_xz_) + f3d(_xyyy_)*qqi(_yy_) &
                + 2.0_rp*f3d(_xyyz_)*qqi(_yz_) + f3d(_xyzz_)*qqi(_zz_)

            hii(_y_,_x_) = hii(_x_,_y_)

            hii(_y_,_y_) = hii(_y_,_y_) + fg(_yy_) * qi &
                - fh(_yyx_)*di(_x_) - fh(_yyy_)*di(_y_) - fh(_yyz_)*di(_z_) &
                + f3d(_yyxx_)*qqi(_xx_) + 2.0_rp*f3d(_yyxy_)*qqi(_xy_) &
                + 2.0_rp*f3d(_yyxz_)*qqi(_xz_) + f3d(_yyyy_)*qqi(_yy_) &
                + 2.0_rp*f3d(_yyyz_)*qqi(_yz_) + f3d(_yyzz_)*qqi(_zz_)

            hii(_x_,_z_) = hii(_x_,_z_) + fg(_xz_) * qi &
                - fh(_xzx_)*di(_x_) - fh(_xzy_)*di(_y_) - fh(_xzz_)*di(_z_) &
                + f3d(_xzxx_)*qqi(_xx_) + 2.0_rp*f3d(_xzxy_)*qqi(_xy_) &
                + 2.0_rp*f3d(_xzxz_)*qqi(_xz_) + f3d(_xzyy_)*qqi(_yy_) &
                + 2.0_rp*f3d(_xzyz_)*qqi(_yz_) + f3d(_xzzz_)*qqi(_zz_)

            hii(_z_,_x_) = hii(_x_,_z_)

            hii(_y_,_z_) = hii(_y_,_z_) + fg(_yz_) * qi &
                - fh(_yzx_)*di(_x_) - fh(_yzy_)*di(_y_) - fh(_yzz_)*di(_z_) &
                + f3d(_yzxx_)*qqi(_xx_) + 2.0_rp*f3d(_yzxy_)*qqi(_xy_) &
                + 2.0_rp*f3d(_yzxz_)*qqi(_xz_) + f3d(_yzyy_)*qqi(_yy_) &
                + 2.0_rp*f3d(_yzyz_)*qqi(_yz_) + f3d(_yzzz_)*qqi(_zz_)

            hii(_z_,_y_) = hii(_y_,_z_)

            hii(_z_,_z_) = hii(_z_,_z_) + fg(_zz_) * qi &
                - fh(_zzx_)*di(_x_) - fh(_zzy_)*di(_y_) - fh(_zzz_)*di(_z_) &
                + f3d(_zzxx_)*qqi(_xx_) + 2.0_rp*f3d(_zzxy_)*qqi(_xy_) &
                + 2.0_rp*f3d(_zzxz_)*qqi(_xz_) + f3d(_zzyy_)*qqi(_yy_) &
                + 2.0_rp*f3d(_zzyz_)*qqi(_yz_) + f3d(_zzzz_)*qqi(_zz_)
        end subroutine hiiterm
!
        subroutine hijterm(s,scalf,i,j,hij)
            use mod_electrostatics
            implicit none
            type(ommp_system), intent(inout), target         :: s
            real(rp),          intent(in)                    :: scalf
            integer(ip),       intent(in)                    :: i, j
            real(rp),          intent(inout), dimension(3,3) :: hij
!
            integer   :: ix, jx
            real(rp)  :: kernel(7), ri(3), rj(3), dr(3), si(3), sj(3), &
                         pij(3), pji(3), nij(3), nji(3), mij(3,3), mji(3,3)
            real(rp)  :: qi, qj, di(3), dj(3), qqi(3,3), qqj(3,3)
            real(rp)  :: mi, mj, rqri, rqrj, cij, cji, dij, didj
!
            ri  = s%top%cmm(:,i)
            rj  = s%top%cmm(:,j)
            dr  = ri - rj
            call coulomb_kernel(dr, 6, kernel)
            kernel = kernel * scalf
            qi  = s%eel%q(1,i)
            di  = s%eel%q(2:4,i)
            qj  = s%eel%q(1,j)
            dj  = s%eel%q(2:4,j)
!
            qqi(1,1) = s%eel%q(4+_xx_,i)
            qqi(2,1) = s%eel%q(4+_yx_,i)
            qqi(3,1) = s%eel%q(4+_zx_,i)
            qqi(1,2) = s%eel%q(4+_xy_,i)
            qqi(2,2) = s%eel%q(4+_yy_,i)
            qqi(3,2) = s%eel%q(4+_zy_,i)
            qqi(1,3) = s%eel%q(4+_xz_,i)
            qqi(2,3) = s%eel%q(4+_yz_,i)
            qqi(3,3) = s%eel%q(4+_zz_,i)
            qqj(1,1) = s%eel%q(4+_xx_,j)
            qqj(2,1) = s%eel%q(4+_yx_,j)
            qqj(3,1) = s%eel%q(4+_zx_,j)
            qqj(1,2) = s%eel%q(4+_xy_,j)
            qqj(2,2) = s%eel%q(4+_yy_,j)
            qqj(3,2) = s%eel%q(4+_zy_,j)
            qqj(1,3) = s%eel%q(4+_xz_,j)
            qqj(2,3) = s%eel%q(4+_yz_,j)
            qqj(3,3) = s%eel%q(4+_zz_,j)
!
!           assemble the various intermediates:
!
            mi    = dot_product(di,dr)
            si    = matmul(qqi,dr)
            rqri  = dot_product(si,dr)
            pij   = matmul(qqi,dj)

            mj    = dot_product(dj,dr)
            sj    = matmul(qqj,dr)
            rqrj  = dot_product(sj,dr)
            pji   = matmul(qqj,di)

            mij   = matmul(qqi,qqj)
            mji   = matmul(qqj,qqi)
            nij   = matmul(mij,dr)
            nji   = matmul(mji,dr)

            cij   = dot_product(pji,dr)  ! di . (Qj r)
            cji   = dot_product(pij,dr)  ! dj . (Qi r)
            dij   = mij(1,1) + mij(2,2) + mij(3,3)
            didj  = dot_product(di,dj)
!
!           charge-charge, charge-dipole, charge-quadrupole,
!           dipole-charge, dipole-dipole, dipole-quadrupole,
!           quadrupole-charge, quadrupole-dipole, quadrupole-quadrupole
!
            hij  = 0.0_rp

            do ix = 1, 3
              do jx = 1, 3

!
!               charge-charge: qi T11 qj
!
                if (jx.eq.ix) hij(ix,jx) = hij(ix,jx) + qi*qj*kernel(2)
                hij(ix,jx) = hij(ix,jx) - 3.0_rp*kernel(3)*qi*qj*dr(ix)*dr(jx)

!
!               charge-dipole: qi T12 dj
!
                hij(ix,jx) = hij(ix,jx) &
                  + qi*( &
                    -15.0_rp*kernel(4)*mj*dr(ix)*dr(jx) &
                    + 3.0_rp*kernel(3)*( &
                        merge(mj,0.0_rp,ix.eq.jx) &
                      + dj(ix)*dr(jx) + dj(jx)*dr(ix) ) )

!
!               dipole-charge: di T21 qj
!
                hij(ix,jx) = hij(ix,jx) &
                  + qj*( &
                     15.0_rp*kernel(4)*mi*dr(ix)*dr(jx) &
                    - 3.0_rp*kernel(3)*( &
                        merge(mi,0.0_rp,ix.eq.jx) &
                      + di(ix)*dr(jx) + di(jx)*dr(ix) ) )

!
!               charge-quadrupole: qi T13 Qj
!
                hij(ix,jx) = hij(ix,jx) &
                  + qi*( &
                    -105.0_rp*kernel(5)*rqrj*dr(ix)*dr(jx) &
                    + 15.0_rp*kernel(4)*( &
                        merge(rqrj,0.0_rp,ix.eq.jx) &
                      + 2.0_rp*sj(ix)*dr(jx) &
                      + 2.0_rp*sj(jx)*dr(ix) ) &
                    - 6.0_rp*kernel(3)*qqj(ix,jx) )

!
!               quadrupole-charge: Qi T31 qj
!
                hij(ix,jx) = hij(ix,jx) &
                  + qj*( &
                    -105.0_rp*kernel(5)*rqri*dr(ix)*dr(jx) &
                    + 15.0_rp*kernel(4)*( &
                        merge(rqri,0.0_rp,ix.eq.jx) &
                      + 2.0_rp*si(ix)*dr(jx) &
                      + 2.0_rp*si(jx)*dr(ix) ) &
                    - 6.0_rp*kernel(3)*qqi(ix,jx) )

!
!               dipole-dipole: di T22 dj
!
                hij(ix,jx) = hij(ix,jx) &
                  + 105.0_rp*kernel(5)*mi*mj*dr(ix)*dr(jx) &
                  - 15.0_rp*kernel(4)*( &
                        merge(mi*mj,0.0_rp,ix.eq.jx) &
                      + di(ix)*dr(jx)*mj &
                      + di(jx)*dr(ix)*mj &
                      + dj(ix)*dr(jx)*mi &
                      + dj(jx)*dr(ix)*mi &
                      + didj*dr(ix)*dr(jx) ) &
                  + 3.0_rp*kernel(3)*( &
                        merge(didj,0.0_rp,ix.eq.jx) &
                      + di(ix)*dj(jx) + dj(ix)*di(jx) )

!
!               dipole-quadrupole: di T23 Qj
!
                hij(ix,jx) = hij(ix,jx) &
                  + 945.0_rp*kernel(6)*mi*rqrj*dr(ix)*dr(jx) &
                  - 105.0_rp*kernel(5)*( &
                        merge(mi*rqrj,0.0_rp,ix.eq.jx) &
                      + di(ix)*dr(jx)*rqrj &
                      + di(jx)*dr(ix)*rqrj &
                      + 2.0_rp*mi*(sj(ix)*dr(jx) + sj(jx)*dr(ix)) &
                      + 2.0_rp*cij*dr(ix)*dr(jx) ) &
                  + 15.0_rp*kernel(4)*( &
                        2.0_rp*merge(cij,0.0_rp,ix.eq.jx) &
                      + 2.0_rp*di(ix)*sj(jx) &
                      + 2.0_rp*di(jx)*sj(ix) &
                      + 2.0_rp*mi*qqj(ix,jx) &
                      + 2.0_rp*pji(ix)*dr(jx) &
                      + 2.0_rp*pji(jx)*dr(ix) )

!
!               quadrupole-dipole: Qi T32 dj = - perm(i,j) of di T23 Qj
!
                hij(ix,jx) = hij(ix,jx) &
                  - 945.0_rp*kernel(6)*mj*rqri*dr(ix)*dr(jx) &
                  + 105.0_rp*kernel(5)*( &
                        merge(mj*rqri,0.0_rp,ix.eq.jx) &
                      + dj(ix)*dr(jx)*rqri &
                      + dj(jx)*dr(ix)*rqri &
                      + 2.0_rp*mj*(si(ix)*dr(jx) + si(jx)*dr(ix)) &
                      + 2.0_rp*cji*dr(ix)*dr(jx) ) &
                  - 15.0_rp*kernel(4)*( &
                        2.0_rp*merge(cji,0.0_rp,ix.eq.jx) &
                      + 2.0_rp*dj(ix)*si(jx) &
                      + 2.0_rp*dj(jx)*si(ix) &
                      + 2.0_rp*mj*qqi(ix,jx) &
                      + 2.0_rp*pij(ix)*dr(jx) &
                      + 2.0_rp*pij(jx)*dr(ix) )

!
!               quadrupole-quadrupole: Qi T33 Qj
!
                hij(ix,jx) = hij(ix,jx) &
                  - 10395.0_rp*kernel(7)*rqri*rqrj*dr(ix)*dr(jx) &
                  + 945.0_rp*kernel(6)*( &
                        merge(rqri*rqrj,0.0_rp,ix.eq.jx) &
                      + 2.0_rp*rqrj*(si(ix)*dr(jx) + si(jx)*dr(ix)) &
                      + 2.0_rp*rqri*(sj(ix)*dr(jx) + sj(jx)*dr(ix)) &
                      + 4.0_rp*dot_product(si,sj)*dr(ix)*dr(jx) ) &
                  - 105.0_rp*kernel(5)*( &
                        4.0_rp*merge(dot_product(si,sj),0.0_rp,ix.eq.jx) &
                      + 2.0_rp*qqi(ix,jx)*rqrj &
                      + 2.0_rp*qqj(ix,jx)*rqri &
                      + 2.0_rp*dij*dr(ix)*dr(jx) &
                      + 4.0_rp*(si(ix)*sj(jx) + sj(ix)*si(jx)) &
                      + 4.0_rp*(nij(ix)*dr(jx) + nij(jx)*dr(ix)) &
                      + 4.0_rp*(nji(ix)*dr(jx) + nji(jx)*dr(ix)) ) &
                  + 15.0_rp*kernel(4)*( &
                        2.0_rp*merge(dij,0.0_rp,ix.eq.jx) &
                      + 4.0_rp*(mij(ix,jx) + mji(ix,jx)) )

              end do
            end do
        end subroutine hijterm

end module
