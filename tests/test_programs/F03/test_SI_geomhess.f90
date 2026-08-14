module test_geomhess
    use ommp_interface

    abstract interface
    subroutine hess_term(s, hess)
        use mod_mmpol, only: ommp_system
        use mod_memory, only: rp
        type(ommp_system), intent(inout), target :: s
        real(rp), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)
    end subroutine
    subroutine dmu_term(s, dmu)
        use mod_mmpol, only: ommp_system
        use mod_memory, only: rp
        type(ommp_system), intent(inout), target :: s
        real(rp), intent(out) :: dmu(3,3,s%top%mm_atoms)
    end subroutine
    end interface

    contains
        subroutine polelec_geomhess_guarded(s, hess)
            !! polelec_geomhess itself, skipped when there is nothing
            !! polarizable to differentiate -- mirrors the guard that
            !! ommp_polelec_geomgrad applies at the gradient level.
            use mod_geomhess, only: polelec_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            if(s%eel%pol_atoms > 0) call polelec_geomhess(s, hess)
        end subroutine

        subroutine full_elec_geomhess(s, hess)
            !! Fixed-multipole + polarization Hessian, combined. Not
            !! "the full Hessian" in the same sense as ommp_full_geomgrad
            !! (which also includes bonded and vdW terms): kept scoped to
            !! electrostatics only, matching full_elec_geomgrad in
            !! test_geomhess_num, so the two remain directly comparable
            !! (bonded and vdW terms are tested separately, as BNDTOT/VDW).
            use mod_geomhess, only: fixedelec_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            call fixedelec_geomhess(s, hess)
            call polelec_geomhess_guarded(s, hess)
        end subroutine

        subroutine bond_geomhess_guarded(s, hess)
            use mod_bonded, only: bond_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            if(s%use_bonded) call bond_geomhess(s%bds, hess)
        end subroutine

        subroutine urey_geomhess_guarded(s, hess)
            use mod_bonded, only: urey_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            if(s%use_bonded) call urey_geomhess(s%bds, hess)
        end subroutine

        subroutine angle_geomhess_guarded(s, hess)
            use mod_bonded, only: angle_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            if(s%use_bonded) call angle_geomhess(s%bds, hess)
        end subroutine

        subroutine strbnd_geomhess_guarded(s, hess)
            use mod_bonded, only: strbnd_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            if(s%use_bonded) call strbnd_geomhess(s%bds, hess)
        end subroutine

        subroutine opb_geomhess_guarded(s, hess)
            use mod_bonded, only: opb_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            if(s%use_bonded) call opb_geomhess(s%bds, hess)
        end subroutine

        subroutine pitors_geomhess_guarded(s, hess)
            use mod_bonded, only: pitors_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            if(s%use_bonded) call pitors_geomhess(s%bds, hess)
        end subroutine

        subroutine torsion_geomhess_guarded(s, hess)
            use mod_bonded, only: torsion_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            if(s%use_bonded) call torsion_geomhess(s%bds, hess)
        end subroutine

        subroutine imptorsion_geomhess_guarded(s, hess)
            use mod_bonded, only: imptorsion_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            if(s%use_bonded) call imptorsion_geomhess(s%bds, hess)
        end subroutine

        subroutine angtor_geomhess_guarded(s, hess)
            use mod_bonded, only: angtor_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            if(s%use_bonded) call angtor_geomhess(s%bds, hess)
        end subroutine

        subroutine strtor_geomhess_guarded(s, hess)
            use mod_bonded, only: strtor_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            if(s%use_bonded) call strtor_geomhess(s%bds, hess)
        end subroutine

        subroutine tortor_geomhess_guarded(s, hess)
            use mod_bonded, only: tortor_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            if(s%use_bonded) call tortor_geomhess(s%bds, hess)
        end subroutine

        subroutine full_bnd_geomhess(s, hess)
            !! Sum of all 11 analytical bonded-term Hessians, combined --
            !! the analytical counterpart of ommp_full_bnd_geomgrad, kept
            !! in the same scope (bonded only, no link-atom corrections,
            !! since no analytical link-atom Hessian exists) so it is
            !! directly comparable to full_bnd_geomgrad in test_geomhess_num.
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            call bond_geomhess_guarded(s, hess)
            call urey_geomhess_guarded(s, hess)
            call angle_geomhess_guarded(s, hess)
            call strbnd_geomhess_guarded(s, hess)
            call opb_geomhess_guarded(s, hess)
            call pitors_geomhess_guarded(s, hess)
            call torsion_geomhess_guarded(s, hess)
            call imptorsion_geomhess_guarded(s, hess)
            call angtor_geomhess_guarded(s, hess)
            call strtor_geomhess_guarded(s, hess)
            call tortor_geomhess_guarded(s, hess)
        end subroutine

        subroutine vdw_geomhess_guarded(s, hess)
            use mod_nonbonded, only: vdw_geomhess
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            if(s%use_nonbonded) call vdw_geomhess(s%vdw, hess)
        end subroutine

        subroutine full_geomhess(s, hess)
            !! Sum of every analytical Hessian contribution implemented so
            !! far (bonded + fixed-multipole electrostatics + polarization
            !! + vdW) -- the Hessian counterpart of ommp_full_geomgrad,
            !! kept in the same scope (no analytical link-atom Hessian
            !! exists, see full_bnd_geomhess) and used as the input for
            !! the harmonic-frequency analysis below.
            !!
            !! Order matters: fixedelec_geomhess (called first, inside
            !! full_elec_geomhess) zeroes every hess(:,:,i,j) block it
            !! touches before accumulating its own contribution into it
            !! (see hiiterm/hijterm in mod_geomhess.F90), instead of
            !! adding to whatever is already there -- so it must run on
            !! a still-zero hess. Bonded and vdW terms use true += ,
            !! hence are safe to accumulate afterwards.
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(inout) :: hess(3,3,s%top%mm_atoms,s%top%mm_atoms)

            call full_elec_geomhess(s, hess)
            call full_bnd_geomhess(s, hess)
            call vdw_geomhess_guarded(s, hess)
        end subroutine

        function total_dipole(s) result(mu)
            !! Total molecular dipole moment (a.u., e*bohr), the same
            !! observable derived in vibrational_analysis's IR-intensity
            !! docstring: mu = sum_i(q_i*r_i + d_i) + 1/2*sum_i(mu^d_i +
            !! mu^p_i), i.e. the field-derivative of the energy at F=0.
            !! Quadrupoles (and higher) don't contribute to the net
            !! dipole of a distributed-multipole expansion.
            use mod_electrostatics, only: prepare_polelec
            use mod_polarization, only: polarization
            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real) :: mu(3)

            integer(ommp_integer) :: i, ipol
            ! ipd's 3rd index: 1 = "D" set, 2 = "P" set (matching
            ! _amoeba_D_/_amoeba_P_ in f_cart_components.h -- this file
            ! is plain .f90, not preprocessed, so the macros aren't
            ! available and the literals are used directly instead).

            if(s%eel%pol_atoms > 0 .and. .not. s%eel%ipd_done) then
                call prepare_polelec(s%eel)
                call polarization(s, s%eel%e_m2d)
            end if

            mu = 0.0
            do i=1, s%top%mm_atoms
                mu = mu + s%eel%q(1,i)*s%top%cmm(:,i) + s%eel%q(2:4,i)
            end do
            do ipol=1, s%eel%pol_atoms
                mu = mu + 0.5_ommp_real*(s%eel%ipd(:,ipol,1) + s%eel%ipd(:,ipol,2))
            end do
        end function

        subroutine full_dipole_geomderiv(s, dmu)
            !! Geometrical derivative of the total dipole moment (see
            !! total_dipole), dmu(alpha,beta,l) = d(mu_alpha)/d(r_l,beta),
            !! rederived (per the user's correction) as minus the mixed
            !! second derivative of the energy wrt an external field F and
            !! r_l, using the field-augmented polarization Lagrangian:
            !!   L_full(mu^d,mu^p;r,F) = L_pol(mu^d,mu^p;r; E^d+F,E^p+F)
            !!                          - F.sum_i(q_i*r_i + d_i)
            !! Stationarity in mu^d/mu^p gives mu = -dL_full/dF =
            !! sum_i(q_i*r_i+d_i) + 1/2*sum_i(mu^d_i+mu^p_i) (total_dipole
            !! above), and since that expression is just an ordinary
            !! (non-stationary) function of r once mu^d(r)/mu^p(r) are
            !! taken as the CPID solutions, its r-derivative needs BOTH
            !! dmu^d/dr and dmu^p/dr (not just one of the two dipole
            !! sets -- mu^p enters the observable with the same 1/2
            !! weight as mu^d, even though it started out as a Lagrange
            !! multiplier) -- exactly the dmud/dmup pair already produced
            !! by build_cpid_rhs+solve_cpid for the polarization Hessian.
            !!
            !! Three pieces, matching total_dipole term by term:
            !!  - monopole: d(q_l*r_l)/d(r_l,beta) = q_l*delta_ab
            !!  - permanent dipole rotation: ddip(beta,alpha,jat,j) is
            !!    d(d_j)_alpha / d(r_{atom(jat)})_beta (see
            !!    rotation_geomgrad's matmul(ddip,E) contraction), so it
            !!    lands transposed into dmu(:,:,atom(jat)).
            !!  - induced dipoles: 0.5*(dmud+dmup), summed over every
            !!    polarizable atom (total_dipole sums mu^d/mu^p over all
            !!    polarizable sites), reusing the exact RHSd/RHSp/dmud/
            !!    dmup already validated for polelec_geomhess.
            use mod_electrostatics, only: prepare_polelec, ommp_electrostatics_type
            use mod_polarization, only: polarization
            use mod_geomhess, only: build_cpid_rhs, solve_cpid
            use mod_rotate_multipoles, only: rotate_multipoles, nref_atoms

            implicit none
            type(ommp_system), intent(inout), target :: s
            real(ommp_real), intent(out) :: dmu(3,3,s%top%mm_atoms)

            type(ommp_electrostatics_type), pointer :: eel
            real(ommp_real), allocatable :: ddip(:,:,:,:), dqua(:,:,:,:,:)
            real(ommp_real), allocatable :: RHSd(:,:), RHSp(:,:), dmud(:,:), dmup(:,:)
            integer(ommp_integer) :: natm, n, npol, i, j, l, jat, nact, atom(4), alpha, beta, ipol

            eel => s%eel
            natm = s%top%mm_atoms
            dmu = 0.0

            ! --- monopole term ---
            do i=1, natm
                dmu(1,1,i) = dmu(1,1,i) + eel%q(1,i)
                dmu(2,2,i) = dmu(2,2,i) + eel%q(1,i)
                dmu(3,3,i) = dmu(3,3,i) + eel%q(1,i)
            end do

            ! --- permanent dipole rotation term ---
            if(eel%amoeba) then
                allocate(ddip(3,3,4,natm))
                allocate(dqua(3,3,3,4,natm))
                call rotate_multipoles(eel, 1_ommp_integer, ddip, dqua)

                do j=1, natm
                    nact = nref_atoms(eel%mol_frame(j))
                    if(nact == 0) cycle
                    atom(1) = j
                    atom(2) = eel%iz(j); if(atom(2) == 0) atom(2) = j
                    atom(3) = eel%ix(j); if(atom(3) == 0) atom(3) = j
                    atom(4) = eel%iy(j); if(atom(4) == 0) atom(4) = j
                    do jat=1, nact
                        dmu(:,:,atom(jat)) = dmu(:,:,atom(jat)) + transpose(ddip(:,:,jat,j))
                    end do
                end do

                deallocate(ddip, dqua)
            end if

            ! --- induced dipole (CPID) term ---
            if(eel%pol_atoms > 0) then
                if(.not. eel%ipd_done) then
                    call prepare_polelec(eel, .false.)
                    call polarization(s, eel%e_M2D)
                end if
                call prepare_polelec(eel, .true., .true.)

                n = 3*natm
                npol = 3*eel%pol_atoms
                allocate(RHSd(npol,n), RHSp(npol,n))
                call build_cpid_rhs(s, RHSd, RHSp)
                allocate(dmud(npol,n), dmup(npol,n))
                call solve_cpid(s, RHSd, RHSp, dmud, dmup)

                do l=1, natm
                    do beta=1, 3
                        do ipol=1, eel%pol_atoms
                            do alpha=1, 3
                                dmu(alpha,beta,l) = dmu(alpha,beta,l) + 0.5_ommp_real * &
                                    (dmud(3*(ipol-1)+alpha, 3*(l-1)+beta) + dmup(3*(ipol-1)+alpha, 3*(l-1)+beta))
                            end do
                        end do
                    end do
                end do

                deallocate(RHSd, RHSp, dmud, dmup)
            end if
        end subroutine

        subroutine print_dmu(n, dmu, name)
            character(len=*) :: name
            integer(ommp_integer) :: n
            real(ommp_real) :: dmu(3,3,n)

            character(len=OMMP_STR_CHAR_MAX) :: msg
            integer(ommp_integer) :: i

            write(msg, "('DMU ', A)") name
            call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-HES")

            do i=1, n
                write(msg, "('I:', I0, 9(' ', E20.12))") i, dmu(:,1,i), dmu(:,2,i), dmu(:,3,i)
                call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-HES")
            end do

            call ommp_message("", OMMP_VERBOSE_NONE, "TEST-HES")
        end subroutine

        subroutine print_hess(n, h, name)
            character(len=*) :: name
            integer(ommp_integer) :: n
            real(ommp_real) :: h(3,3,n,n)

            character(len=OMMP_STR_CHAR_MAX) :: msg
            integer(ommp_integer) :: i, j

            write(msg, "('Hess ', A)") name
            call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-HES")

            do i=1, n
                do j=1, n
                    write(msg, "('IJ:', I0, ':', I0, 9(' ', E20.12))") &
                        i, j, h(:,1,i,j), h(:,2,i,j), h(:,3,i,j)
                    call ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-HES")
                end do
            end do

            call ommp_message("", OMMP_VERBOSE_NONE, "TEST-HES")
        end subroutine

        subroutine ana_hess_print(sys, hess_f, n)
            type(ommp_system) :: sys
            procedure(hess_term), pointer :: hess_f
            character(len=*) :: n

            real(ommp_real), allocatable, dimension(:,:,:,:) :: hmm
            integer(ommp_integer) :: natm

            natm = sys%top%mm_atoms
            allocate(hmm(3,3,natm,natm))
            hmm = 0.0
            call hess_f(sys, hmm)
            ! one extra power of ANG2AU vs. the gradient's own scaling,
            ! since the Hessian is a second position-derivative
            hmm = hmm * OMMP_AU2KCALMOL*OMMP_ANG2AU*OMMP_ANG2AU
            call print_hess(natm, hmm, n)
            deallocate(hmm)
        end subroutine

        subroutine vibrational_analysis(sys, hess_f, dmu_f)
            !! Builds the mass-weighted Hessian K_ij = H_ij / sqrt(m_i m_j)
            !! (the 3x3 Cartesian block of atoms i,j implicit in that
            !! notation) out of the full analytical Hessian, diagonalizes
            !! it and reports the harmonic frequencies in cm^-1, ascending,
            !! together with the double-harmonic IR intensity of each mode
            !! (in km/mol). Negative eigenvalues (saddle-point directions)
            !! are reported as imaginary frequencies: the modulus,
            !! followed by 'i'.
            !!
            !! The IR intensity of mode k is the standard double-harmonic
            !! result I_k = (N_A*pi)/(3c^2) * |d(mu)/dQ_k|^2. The
            !! eigenvectors K(:,k) (dsyev job 'V') give x_{i,a}[Ang] =
            !! sum_k K(3*(i-1)+a,k)*Q_k/sqrt(m_i[amu]), so
            !! d(mu_alpha)/dQ_k = sum_{i,a} dmu(alpha,a,i)*K(3(i-1)+a,k)
            !! /sqrt(m_i). dmu_f (full_dipole_geomderiv) returns raw a.u.
            !! d(mu)/d(r); as derived in that routine's docstring, a
            !! mixed derivative like d(mu)/dr is numerically the same
            !! whether mu/r are expressed in a consistent bohr- or
            !! Angstrom-based pair of units, so it can be contracted
            !! directly against K (defined for x in Angstrom, since hmm
            !! below is pre-scaled to kcal/mol/Ang^2) with no extra
            !! Ang2au factor.
            implicit none
            type(ommp_system), intent(inout), target :: sys
            procedure(hess_term), pointer :: hess_f
            procedure(dmu_term), pointer :: dmu_f

            ! cm^-1 per sqrt(kcal/mol/Ang^2/amu), i.e. the unit system
            ! that hess_f is reported in once scaled by
            ! OMMP_AU2KCALMOL*OMMP_ANG2AU**2 (see ana_hess_print) with
            ! masses (sys%top%atmass) in amu.
            real(ommp_real), parameter :: cm1_per_sqrt_unit = 108.5913586111_ommp_real
            ! km/mol per (e/amu^0.5)^2: I_k[cm/mol] = (N_A*pi)/(3c^2) *
            ! (dmu/dQ_k [esu/g^0.5])^2, with dmu/dQ_k [e/amu^0.5]
            ! converted to [esu/g^0.5] via e_esu/sqrt(amu_in_g)
            ! (CODATA: e=1.602176634e-19 C, 1 esu=3.33564095198e-10 C,
            ! amu=1.66053906660e-24 g, N_A=6.02214076e23, c=2.99792458e10
            ! cm/s), then cm/mol -> km/mol (x1e-5). Cross-checked against
            ! the textbook ~42.256 km/mol constant (which uses D instead
            ! of e*Ang for the dipole: 1 e*Ang = 4.8032 D, and
            ! 4.8032^2*42.256 = 974.88, matching).
            real(ommp_real), parameter :: kmmol_per_sqrt_unit2 = 974.8801113_ommp_real

            real(ommp_real), allocatable :: hmm(:,:,:,:), K(:,:), w(:), work(:), dmu(:,:,:)
            real(ommp_real) :: lwork_query(1), freq, dmudq(3), intens
            character(len=OMMP_STR_CHAR_MAX) :: msg
            integer(ommp_integer) :: natm, n, i, j, a, b, info, lwork, imode, idx

            if(.not. sys%top%atmass_initialized) then
                ! Atomic masses are only ever set while parsing a .prm
                ! file (mod_prm.F90); systems loaded from a bare .mmp
                ! file never populate them, so the mass-weighted analysis
                ! is skipped rather than dividing by an undefined mass.
                call ommp_message("Atomic masses not available (no .prm "// &
                    "file was used), skipping vibrational analysis", &
                    OMMP_VERBOSE_NONE, "FREQ")
                return
            end if

            natm = sys%top%mm_atoms
            n = 3*natm

            allocate(hmm(3,3,natm,natm))
            hmm = 0.0
            call hess_f(sys, hmm)
            hmm = hmm * OMMP_AU2KCALMOL*OMMP_ANG2AU*OMMP_ANG2AU

            allocate(K(n,n))
            do i=1, natm
                do j=1, natm
                    do b=1, 3
                        do a=1, 3
                            K(3*(i-1)+a,3*(j-1)+b) = hmm(a,b,i,j) / &
                                sqrt(sys%top%atmass(i)*sys%top%atmass(j))
                        end do
                    end do
                end do
            end do
            deallocate(hmm)

            allocate(dmu(3,3,natm))
            call dmu_f(sys, dmu)

            allocate(w(n))
            call dsyev('V', 'U', n, K, n, w, lwork_query, -1_ommp_integer, info)
            lwork = int(lwork_query(1), ommp_integer)
            allocate(work(lwork))
            call dsyev('V', 'U', n, K, n, w, work, lwork, info)
            deallocate(work)

            if(info /= 0) then
                write(6, *) "Diagonalization of the mass-weighted Hessian failed, info=", info
                call exit(1)
            end if

            call ommp_message("Harmonic frequencies (cm^-1) and IR intensities (km/mol)", &
                OMMP_VERBOSE_NONE, "FREQ")
            do imode=1, n
                dmudq = 0.0_ommp_real
                do i=1, natm
                    do a=1, 3
                        idx = 3*(i-1)+a
                        dmudq = dmudq + dmu(:,a,i) * K(idx,imode) / sqrt(sys%top%atmass(i))
                    end do
                end do
                intens = kmmol_per_sqrt_unit2 * sum(dmudq**2)

                if(w(imode) >= 0.0) then
                    freq = cm1_per_sqrt_unit * sqrt(w(imode))
                    write(msg, "(I0, ': ', F0.2, ' cm^-1   I = ', F0.2, ' km/mol')") imode, freq, intens
                else
                    freq = cm1_per_sqrt_unit * sqrt(-w(imode))
                    write(msg, "(I0, ': ', F0.2, 'i cm^-1   I = ', F0.2, ' km/mol')") imode, freq, intens
                end if
                call ommp_message(msg, OMMP_VERBOSE_NONE, "FREQ")
            end do
            call ommp_message("", OMMP_VERBOSE_NONE, "FREQ")

            deallocate(w, K, dmu)
        end subroutine
end module

program test_SI_geomhess
    use iso_c_binding, only: c_char
    use ommp_interface
    use test_geomhess
    use mod_geomhess, only: fixedelec_geomhess

    implicit none
    character(kind=c_char, len=120), dimension(2) :: args
    integer :: narg
    type(ommp_system), pointer :: my_system
    type(ommp_qm_helper), pointer :: my_qmh

    procedure(hess_term), pointer :: ht
    procedure(dmu_term), pointer :: dmt

    narg = command_argument_count()
    if (narg /= 2) then
        write(6, *) "Syntax expected "
        write(6, *) "   $ test_SI_geomhess.exe <JSON FILE> <OUTPUT FILE>"
        call exit(1)
    else
        call get_command_argument(1, args(1))
        call get_command_argument(2, args(2))

        call ommp_smartinput(trim(args(1)), my_system, my_qmh)
        call ommp_set_outputfile(trim(args(2)))

        ht => fixedelec_geomhess
        call ana_hess_print(my_system, ht, "EM")
        ht => polelec_geomhess_guarded
        call ana_hess_print(my_system, ht, "EP")
        ht => full_elec_geomhess
        call ana_hess_print(my_system, ht, "ETOT")

        ht => bond_geomhess_guarded
        call ana_hess_print(my_system, ht, "BOND")
        ht => urey_geomhess_guarded
        call ana_hess_print(my_system, ht, "UREY")
        ht => angle_geomhess_guarded
        call ana_hess_print(my_system, ht, "ANGLE")
        ht => strbnd_geomhess_guarded
        call ana_hess_print(my_system, ht, "STRBND")
        ht => opb_geomhess_guarded
        call ana_hess_print(my_system, ht, "OPB")
        ht => pitors_geomhess_guarded
        call ana_hess_print(my_system, ht, "PITORS")
        ht => torsion_geomhess_guarded
        call ana_hess_print(my_system, ht, "TORSION")
        ht => imptorsion_geomhess_guarded
        call ana_hess_print(my_system, ht, "IMPTORSION")
        ht => angtor_geomhess_guarded
        call ana_hess_print(my_system, ht, "ANGTOR")
        ht => strtor_geomhess_guarded
        call ana_hess_print(my_system, ht, "STRTOR")
        ht => tortor_geomhess_guarded
        call ana_hess_print(my_system, ht, "TORTOR")
        ht => full_bnd_geomhess
        call ana_hess_print(my_system, ht, "BNDTOT")

        ht => vdw_geomhess_guarded
        call ana_hess_print(my_system, ht, "VDW")

        ht => full_geomhess
        call ana_hess_print(my_system, ht, "FULLHESS")

        block
            real(ommp_real), allocatable :: dmu(:,:,:)
            allocate(dmu(3,3,my_system%top%mm_atoms))
            call full_dipole_geomderiv(my_system, dmu)
            call print_dmu(my_system%top%mm_atoms, dmu, "DIPDERIV")
            deallocate(dmu)
        end block

        dmt => full_dipole_geomderiv
        call vibrational_analysis(my_system, ht, dmt)

        if(associated(my_qmh)) call ommp_terminate_qm_helper(my_qmh)
        if(associated(my_system)) call ommp_terminate(my_system)
    end if
end program
