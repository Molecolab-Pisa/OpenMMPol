#include "f_cart_components.h"

module mod_rotate_multipoles
    !! Rotation of AMOEBA distributed multipoles from the molecular frame
    !! to the lab frame, and the geometrical gradient ("torque") terms
    !! that stem from the dependence of that rotation on atomic positions.
    implicit none
    private

    public :: rotate_multipoles, rotation_geomgrad, rotation_geomhess, &
              rotation_geomhess_pair, rotation_geomhess_pair8

    contains

pure function nref_atoms(mol_frame) result(nact)
    !! Number of frame-defining atoms actually used by a rotation
    !! convention. As can be seen in rotation_matrix, the active
    !! reference-atom slots for any convention are always a prefix
    !! self[,iz[,ix[,iy]]] of the full list (this is why iy is kept last:
    !! not every convention uses it, and ix/iy are used progressively
    !! less often). So looping 1,nact instead of 1,4 in rotate_multipoles/
    !! rotation_geomgrad/rotation_geomhess skips exactly the slots whose
    !! ddip/dqua/d2dip/d2qua are provably zero, without changing the
    !! result.

    use mod_memory, only: ip

    implicit none

    integer(ip), intent(in) :: mol_frame
    integer(ip) :: nact

    select case(mol_frame)
        case(0)
            nact = 0 ! no molecular frame, identity rotation
        case(3)
            nact = 2 ! z-only: self, iz
        case(1, 2)
            nact = 3 ! z-then-x, bisector: self, iz, ix
        case(4, 5)
            nact = 4 ! z-bisector, 3-fold: self, iz, ix, iy
        case default
            nact = 4 ! rotation_matrix will fatal_error on this anyway
    end select
end function nref_atoms

subroutine rotation_geomgrad(eel, E, Egrd, ddip, dqua, grad)
    !! Adds to grad the "torque" contribution to the geometrical gradient
    !! that stems from the dependence of the rotated dipoles and
    !! quadrupoles on the positions of the atoms that define each site's
    !! molecular frame (self, iz, ix, iy -- see [[rotate_multipoles]]).
    !!
    !! ddip/dqua must have been computed beforehand by a call to
    !! rotate_multipoles(eel, 1, ddip, dqua). They are taken as arguments,
    !! rather than (re-)computed in here, so that eel can stay intent(in):
    !! E/Egrd are typically actual arguments aliased to eel%E_M2M/
    !! eel%Egrd_M2M, and calling rotate_multipoles from in here would
    !! require eel intent(inout), aliasing it with E/Egrd -- undefined
    !! behaviour in Fortran.

    use mod_memory, only: ip, rp
    use mod_electrostatics, only: ommp_electrostatics_type

    implicit none

    type(ommp_electrostatics_type), intent(in) :: eel
    real(rp), intent(in) :: E(3, eel%top%mm_atoms), Egrd(6, eel%top%mm_atoms)
    real(rp), intent(in) :: ddip(3,3,4,eel%top%mm_atoms)
    real(rp), intent(in) :: dqua(3,3,3,4,eel%top%mm_atoms)
    real(rp), dimension(3, eel%top%mm_atoms), intent(inout) :: grad

    integer(ip) :: j, jat, nact, atom(4)
    logical :: frozen(4)

    do j = 1, eel%top%mm_atoms
        nact = nref_atoms(eel%mol_frame(j))
        if(nact == 0) cycle

        atom(_self_) = j
        atom(_iz_) = eel%iz(j); if(atom(_iz_) == 0) atom(_iz_) = j
        atom(_ix_) = eel%ix(j); if(atom(_ix_) == 0) atom(_ix_) = j
        atom(_iy_) = eel%iy(j); if(atom(_iy_) == 0) atom(_iy_) = j

        frozen = .false.
        if(eel%top%use_frozen) then
            do jat = 1, nact
                frozen(jat) = eel%top%frozen(atom(jat))
            end do
            if(all(frozen(1:nact))) cycle
        end if

        do jat = 1, nact
            if(frozen(jat)) cycle

            ! torque from the dipole: grad -= (d_beta dip) . E
            grad(:,atom(jat)) = grad(:,atom(jat)) - matmul(ddip(:,:,jat,j), E(:,j))

            ! torque from the quadrupole: grad += (d_beta qua) : Egrd
            grad(:,atom(jat)) = grad(:,atom(jat)) &
                              + dqua(:,_x_,_x_,jat,j) * Egrd(_xx_,j) &
                              + dqua(:,_y_,_y_,jat,j) * Egrd(_yy_,j) &
                              + dqua(:,_z_,_z_,jat,j) * Egrd(_zz_,j) &
                              + 2.0_rp*(dqua(:,_x_,_y_,jat,j) * Egrd(_xy_,j) &
                              +         dqua(:,_x_,_z_,jat,j) * Egrd(_xz_,j) &
                              +         dqua(:,_y_,_z_,jat,j) * Egrd(_yz_,j))
        end do
    end do

end subroutine rotation_geomgrad

subroutine rotation_geomhess(eel, E, Egrd, EHes, ddip, dqua, d2dip, d2qua, hess)
    !! Adds to hess the "on-site" contributions to the geometrical Hessian
    !! that stem from the dependence of the rotated dipoles/quadrupoles on
    !! the positions of their frame-defining atoms (self, iz, ix, iy):
    !!   - (d_p Theta_j) . Phi_j^{L+1}(Theta), landing in hess(:,:,p,j),
    !!     and its transpose (d_p Theta_j) . Phi_j^{L+1}(Theta) read the
    !!     other way round, landing in hess(:,:,j,p) -- this is one rank
    !!     higher than rotation_geomgrad's torque term (Egrd instead of E
    !!     for the dipole, EHes instead of Egrd for the quadrupole);
    !!   - (d_p d_q Theta_j) . Phi_j^L(Theta), landing in hess(:,:,p,q):
    !!     the direct second-derivative analogue of rotation_geomgrad's
    !!     torque term, using d2dip/d2qua against the same E/Egrd fields.
    !! ddip/dqua/d2dip/d2qua must have been computed beforehand by a call
    !! to rotate_multipoles(eel, 2, ddip, dqua, d2dip, d2qua).

    use mod_memory, only: ip, rp
    use mod_electrostatics, only: ommp_electrostatics_type

    implicit none

    type(ommp_electrostatics_type), intent(in) :: eel
    real(rp), intent(in) :: E(3, eel%top%mm_atoms), Egrd(6, eel%top%mm_atoms)
    real(rp), intent(in) :: EHes(10, eel%top%mm_atoms)
    real(rp), intent(in) :: ddip(3,3,4,eel%top%mm_atoms)
    real(rp), intent(in) :: dqua(3,3,3,4,eel%top%mm_atoms)
    real(rp), intent(in) :: d2dip(3,3,3,4,4,eel%top%mm_atoms)
    real(rp), intent(in) :: d2qua(3,3,3,3,4,4,eel%top%mm_atoms)
    real(rp), dimension(3,3,eel%top%mm_atoms,eel%top%mm_atoms), intent(inout) :: hess

    integer(ip) :: j, jat, jat2, delta, p, q, atom(4)
    logical :: frozen(4)
    real(rp) :: gmat(3,3), x(3,3)
    integer(ip) :: nact

    do j = 1, eel%top%mm_atoms
        nact = nref_atoms(eel%mol_frame(j))
        if(nact == 0) cycle

        atom(_self_) = j
        atom(_iz_) = eel%iz(j); if(atom(_iz_) == 0) atom(_iz_) = j
        atom(_ix_) = eel%ix(j); if(atom(_ix_) == 0) atom(_ix_) = j
        atom(_iy_) = eel%iy(j); if(atom(_iy_) == 0) atom(_iy_) = j

        frozen = .false.
        if(eel%top%use_frozen) then
            do jat = 1, nact
                frozen(jat) = eel%top%frozen(atom(jat))
            end do
            if(all(frozen(1:nact))) cycle
        end if

        ! unpack Egrd(j) into a full 3x3 matrix for the dipole part of the
        ! one-rank-up torque term below.
        gmat(_x_,_x_) = Egrd(_xx_,j); gmat(_x_,_y_) = Egrd(_xy_,j); gmat(_x_,_z_) = Egrd(_xz_,j)
        gmat(_y_,_x_) = Egrd(_xy_,j); gmat(_y_,_y_) = Egrd(_yy_,j); gmat(_y_,_z_) = Egrd(_yz_,j)
        gmat(_z_,_x_) = Egrd(_xz_,j); gmat(_z_,_y_) = Egrd(_yz_,j); gmat(_z_,_z_) = Egrd(_zz_,j)

        do jat = 1, nact
            if(frozen(jat)) cycle
            p = atom(jat)

            ! (d_p dip_j) . Egrd_j
            x = matmul(ddip(:,:,jat,j), gmat)

            ! - (d_p qua_j) : EHes_j
            x(:,_x_) = x(:,_x_) &
                     - (dqua(:,_x_,_x_,jat,j) * EHes(_xxx_,j) + dqua(:,_y_,_y_,jat,j) * EHes(_yyx_,j) &
                     +  dqua(:,_z_,_z_,jat,j) * EHes(_zzx_,j) &
                     +  2.0_rp*(dqua(:,_x_,_y_,jat,j) * EHes(_xyx_,j) + dqua(:,_x_,_z_,jat,j) * EHes(_xzx_,j) &
                     +          dqua(:,_y_,_z_,jat,j) * EHes(_yzx_,j)))
            x(:,_y_) = x(:,_y_) &
                     - (dqua(:,_x_,_x_,jat,j) * EHes(_xxy_,j) + dqua(:,_y_,_y_,jat,j) * EHes(_yyy_,j) &
                     +  dqua(:,_z_,_z_,jat,j) * EHes(_zzy_,j) &
                     +  2.0_rp*(dqua(:,_x_,_y_,jat,j) * EHes(_xyy_,j) + dqua(:,_x_,_z_,jat,j) * EHes(_xzy_,j) &
                     +          dqua(:,_y_,_z_,jat,j) * EHes(_yzy_,j)))
            x(:,_z_) = x(:,_z_) &
                     - (dqua(:,_x_,_x_,jat,j) * EHes(_xxz_,j) + dqua(:,_y_,_y_,jat,j) * EHes(_yyz_,j) &
                     +  dqua(:,_z_,_z_,jat,j) * EHes(_zzz_,j) &
                     +  2.0_rp*(dqua(:,_x_,_y_,jat,j) * EHes(_xyz_,j) + dqua(:,_x_,_z_,jat,j) * EHes(_xzz_,j) &
                     +          dqua(:,_y_,_z_,jat,j) * EHes(_yzz_,j)))

            hess(:,:,p,j) = hess(:,:,p,j) + x
            hess(:,:,j,p) = hess(:,:,j,p) + transpose(x)

            ! (d_p d_q dip_j) . E_j - (d_p d_q qua_j) : Egrd_j
            do jat2 = 1, nact
                if(frozen(jat2)) cycle
                q = atom(jat2)

                do delta = 1, 3
                    hess(:,delta,p,q) = hess(:,delta,p,q) &
                                       - matmul(d2dip(:,delta,:,jat,jat2,j), E(:,j)) &
                                       + d2qua(:,delta,_x_,_x_,jat,jat2,j) * Egrd(_xx_,j) &
                                       + d2qua(:,delta,_y_,_y_,jat,jat2,j) * Egrd(_yy_,j) &
                                       + d2qua(:,delta,_z_,_z_,jat,jat2,j) * Egrd(_zz_,j) &
                                       + 2.0_rp*(d2qua(:,delta,_x_,_y_,jat,jat2,j) * Egrd(_xy_,j) &
                                       +         d2qua(:,delta,_x_,_z_,jat,jat2,j) * Egrd(_xz_,j) &
                                       +         d2qua(:,delta,_y_,_z_,jat,jat2,j) * Egrd(_yz_,j))
                end do
            end do
        end do
    end do

end subroutine rotation_geomhess

subroutine rotation_geomhess_pair(eel, scalf, i, j, ddip, dqua, hess)
    !! Adds to hess the "cross" contributions to the geometrical Hessian
    !! (terms 5 and 6 of eq. Hess1) that come from the field generated at
    !! atom i by the frame-differentiated multipoles of atom j:
    !!   T5: Theta_i . Phi_i^{L+1}(d_p Theta_j), landing in hess(:,:,p,i);
    !!   T6: Theta_i . Phi_i^{L+1}(d_p Theta_j) read the other way round
    !!       (i.e. with i and p swapped), landing in hess(:,:,i,p) -- the
    !!       transpose of T5's block.
    !! for each frame atom p of j (this is the same relationship between
    !! T5/T6 as between T3/T4 in rotation_geomhess, just with the field
    !! recomputed pairwise at i instead of reused from eel%*_M2M, since
    !! here it is generated by j's differentiated multipoles rather than
    !! by the plain ones).
    !!
    !! Meant to be called once per ordered pair (i,j) exactly where
    !! hijterm is (same screening scalf) -- the two orderings (i,j) and
    !! (j,i), visited separately by that loop, together cover all four
    !! (target,source) combinations that make up T5+T6.
    !!
    !! ddip/dqua must have been computed beforehand by a call to
    !! rotate_multipoles(eel, >=1, ddip, dqua, ...).

    use mod_memory, only: ip, rp
    use mod_electrostatics, only: ommp_electrostatics_type, coulomb_kernel, &
                                   mu_elec_prop, quad_elec_prop

    implicit none

    type(ommp_electrostatics_type), intent(in) :: eel
    real(rp), intent(in) :: scalf
    integer(ip), intent(in) :: i, j
    real(rp), intent(in) :: ddip(3,3,4,eel%top%mm_atoms)
    real(rp), intent(in) :: dqua(3,3,3,4,eel%top%mm_atoms)
    real(rp), dimension(3,3,eel%top%mm_atoms,eel%top%mm_atoms), intent(inout) :: hess

    integer(ip) :: jat, dir, nact, atom(4)
    real(rp) :: dr(3), kernel(7), qpk(6)
    real(rp) :: qi, di(3), qqi(6)
    real(rp) :: tmpV, tmpE(3), tmpEgrd(6), tmpHE(10), tmpD3E(15)
    real(rp) :: gmat(3,3), y(3,3)

    nact = nref_atoms(eel%mol_frame(j))
    if(nact == 0) return

    atom(_self_) = j
    atom(_iz_) = eel%iz(j); if(atom(_iz_) == 0) atom(_iz_) = j
    atom(_ix_) = eel%ix(j); if(atom(_ix_) == 0) atom(_ix_) = j
    atom(_iy_) = eel%iy(j); if(atom(_iy_) == 0) atom(_iy_) = j

    dr = eel%top%cmm(:,i) - eel%top%cmm(:,j)
    call coulomb_kernel(dr, 6, kernel)
    kernel = kernel * scalf

    qi = eel%q(1,i)
    di = eel%q(2:4,i)
    qqi(_xx_) = eel%q(4+_xx_,i); qqi(_xy_) = eel%q(4+_xy_,i); qqi(_yy_) = eel%q(4+_yy_,i)
    qqi(_xz_) = eel%q(4+_xz_,i); qqi(_yz_) = eel%q(4+_yz_,i); qqi(_zz_) = eel%q(4+_zz_,i)

    do jat = 1, nact
        if(eel%top%use_frozen) then
            if(eel%top%frozen(atom(jat))) cycle
        end if

        do dir = 1, 3
            tmpV = 0.0_rp; tmpE = 0.0_rp; tmpEgrd = 0.0_rp; tmpHE = 0.0_rp; tmpD3E = 0.0_rp

            ! pack (d_dir qua_j)(:,:,jat) as (xx,xy,yy,xz,yz,zz), the
            ! convention expected by quad_elec_prop.
            qpk(_xx_) = dqua(dir,_x_,_x_,jat,j); qpk(_xy_) = dqua(dir,_x_,_y_,jat,j)
            qpk(_yy_) = dqua(dir,_y_,_y_,jat,j); qpk(_xz_) = dqua(dir,_x_,_z_,jat,j)
            qpk(_yz_) = dqua(dir,_y_,_z_,jat,j); qpk(_zz_) = dqua(dir,_z_,_z_,jat,j)

            ! field at i generated by the pseudo dipole/quadrupole
            ! (d_dir dip_j)/(d_dir qua_j) sitting at j.
            call mu_elec_prop(ddip(dir,:,jat,j), dr, kernel, &
                              .false., tmpV, .true., tmpE, .true., tmpEgrd, &
                              .true., tmpHE, .false., tmpD3E)
            call quad_elec_prop(qpk, dr, kernel, &
                                .false., tmpV, .true., tmpE, .true., tmpEgrd, &
                                .true., tmpHE, .false., tmpD3E)

            ! contract with Theta_i: Phi^1(...) = -E, Phi^2(...) = +Egrd,
            ! Phi^3(...) = -EHes (same sign convention as rotation_geomhess).
            gmat(_x_,_x_) = tmpEgrd(_xx_); gmat(_x_,_y_) = tmpEgrd(_xy_); gmat(_x_,_z_) = tmpEgrd(_xz_)
            gmat(_y_,_x_) = tmpEgrd(_xy_); gmat(_y_,_y_) = tmpEgrd(_yy_); gmat(_y_,_z_) = tmpEgrd(_yz_)
            gmat(_z_,_x_) = tmpEgrd(_xz_); gmat(_z_,_y_) = tmpEgrd(_yz_); gmat(_z_,_z_) = tmpEgrd(_zz_)

            y(dir,:) = -qi * tmpE + matmul(di, gmat)

            y(dir,_x_) = y(dir,_x_) &
                       - (qqi(_xx_) * tmpHE(_xxx_) + qqi(_yy_) * tmpHE(_yyx_) + qqi(_zz_) * tmpHE(_zzx_) &
                       +  2.0_rp*(qqi(_xy_) * tmpHE(_xyx_) + qqi(_xz_) * tmpHE(_xzx_) + qqi(_yz_) * tmpHE(_yzx_)))
            y(dir,_y_) = y(dir,_y_) &
                       - (qqi(_xx_) * tmpHE(_xxy_) + qqi(_yy_) * tmpHE(_yyy_) + qqi(_zz_) * tmpHE(_zzy_) &
                       +  2.0_rp*(qqi(_xy_) * tmpHE(_xyy_) + qqi(_xz_) * tmpHE(_xzy_) + qqi(_yz_) * tmpHE(_yzy_)))
            y(dir,_z_) = y(dir,_z_) &
                       - (qqi(_xx_) * tmpHE(_xxz_) + qqi(_yy_) * tmpHE(_yyz_) + qqi(_zz_) * tmpHE(_zzz_) &
                       +  2.0_rp*(qqi(_xy_) * tmpHE(_xyz_) + qqi(_xz_) * tmpHE(_xzz_) + qqi(_yz_) * tmpHE(_yzz_)))
        end do

        hess(:,:,atom(jat),i) = hess(:,:,atom(jat),i) + y
        hess(:,:,i,atom(jat)) = hess(:,:,i,atom(jat)) + transpose(y)
    end do

end subroutine rotation_geomhess_pair

subroutine rotation_geomhess_pair8(eel, scalf, i, j, ddip, dqua, hess)
    !! Adds to hess term 8 of eq. Hess1: the field generated at i by j's
    !! frame-differentiated multipoles, contracted with i's OWN
    !! frame-differentiated multipoles (rather than i's plain multipole,
    !! as in rotation_geomhess_pair for T5/T6) -- landing in
    !! hess(:,:,frame_atom_of_i,frame_atom_of_j).
    !!
    !! Unlike T5/T6, this pairs same-rank quantities (Phi_i^L, not
    !! Phi_i^{L+1}), since both sides are already differentiated once, so
    !! it only needs E/Egrd (not Egrd/EHes) from mu_elec_prop/
    !! quad_elec_prop.
    !!
    !! Called once per ordered pair (i,j), same as hijterm/
    !! rotation_geomhess_pair (same screening scalf); unlike those, no
    !! transpose-fill trick is used here, since the two orderings (i,j)
    !! and (j,i), visited separately by that loop, are already the two
    !! distinct terms of the underlying Sum_i Sum_{m!=i} double sum.
    !!
    !! ddip/dqua must have been computed beforehand by a call to
    !! rotate_multipoles(eel, >=1, ddip, dqua, ...).

    use mod_memory, only: ip, rp
    use mod_electrostatics, only: ommp_electrostatics_type, coulomb_kernel, &
                                   mu_elec_prop, quad_elec_prop

    implicit none

    type(ommp_electrostatics_type), intent(in) :: eel
    real(rp), intent(in) :: scalf
    integer(ip), intent(in) :: i, j
    real(rp), intent(in) :: ddip(3,3,4,eel%top%mm_atoms)
    real(rp), intent(in) :: dqua(3,3,3,4,eel%top%mm_atoms)
    real(rp), dimension(3,3,eel%top%mm_atoms,eel%top%mm_atoms), intent(inout) :: hess

    integer(ip) :: jat, jat2, dir, dir2, nacti, nactj, atomi(4), atomj(4)
    real(rp) :: dr(3), kernel(7), qpk(6)
    real(rp) :: tmpV, tmpE(3), tmpEgrd(6), tmpHE(10), tmpD3E(15)
    real(rp) :: gmat(3,3), y

    nacti = nref_atoms(eel%mol_frame(i))
    nactj = nref_atoms(eel%mol_frame(j))
    if(nacti == 0 .or. nactj == 0) return

    atomi(_self_) = i
    atomi(_iz_) = eel%iz(i); if(atomi(_iz_) == 0) atomi(_iz_) = i
    atomi(_ix_) = eel%ix(i); if(atomi(_ix_) == 0) atomi(_ix_) = i
    atomi(_iy_) = eel%iy(i); if(atomi(_iy_) == 0) atomi(_iy_) = i

    atomj(_self_) = j
    atomj(_iz_) = eel%iz(j); if(atomj(_iz_) == 0) atomj(_iz_) = j
    atomj(_ix_) = eel%ix(j); if(atomj(_ix_) == 0) atomj(_ix_) = j
    atomj(_iy_) = eel%iy(j); if(atomj(_iy_) == 0) atomj(_iy_) = j

    dr = eel%top%cmm(:,i) - eel%top%cmm(:,j)
    call coulomb_kernel(dr, 6, kernel)
    kernel = kernel * scalf

    do jat2 = 1, nactj
        if(eel%top%use_frozen) then
            if(eel%top%frozen(atomj(jat2))) cycle
        end if

        do dir2 = 1, 3
            ! field at i generated by the pseudo dipole/quadrupole
            ! (d_dir2 dip_j)/(d_dir2 qua_j) sitting at j -- Phi_i^L only
            ! needs up to rank L=2 (do_grdE) here, unlike T5/T6's
            ! Phi^{L+1} (do_HE).
            tmpV = 0.0_rp; tmpE = 0.0_rp; tmpEgrd = 0.0_rp; tmpHE = 0.0_rp; tmpD3E = 0.0_rp

            qpk(_xx_) = dqua(dir2,_x_,_x_,jat2,j); qpk(_xy_) = dqua(dir2,_x_,_y_,jat2,j)
            qpk(_yy_) = dqua(dir2,_y_,_y_,jat2,j); qpk(_xz_) = dqua(dir2,_x_,_z_,jat2,j)
            qpk(_yz_) = dqua(dir2,_y_,_z_,jat2,j); qpk(_zz_) = dqua(dir2,_z_,_z_,jat2,j)

            call mu_elec_prop(ddip(dir2,:,jat2,j), dr, kernel, &
                              .false., tmpV, .true., tmpE, .true., tmpEgrd, &
                              .false., tmpHE, .false., tmpD3E)
            call quad_elec_prop(qpk, dr, kernel, &
                                .false., tmpV, .true., tmpE, .true., tmpEgrd, &
                                .false., tmpHE, .false., tmpD3E)

            gmat(_x_,_x_) = tmpEgrd(_xx_); gmat(_x_,_y_) = tmpEgrd(_xy_); gmat(_x_,_z_) = tmpEgrd(_xz_)
            gmat(_y_,_x_) = tmpEgrd(_xy_); gmat(_y_,_y_) = tmpEgrd(_yy_); gmat(_y_,_z_) = tmpEgrd(_yz_)
            gmat(_z_,_x_) = tmpEgrd(_xz_); gmat(_z_,_y_) = tmpEgrd(_yz_); gmat(_z_,_z_) = tmpEgrd(_zz_)

            do jat = 1, nacti
                if(eel%top%use_frozen) then
                    if(eel%top%frozen(atomi(jat))) cycle
                end if

                do dir = 1, 3
                    ! contract with (d_dir dip_i)/(d_dir qua_i):
                    ! Phi_i^1(...) = -E, Phi_i^2(...) = +Egrd
                    y = -dot_product(ddip(dir,:,jat,i), tmpE) &
                      + (dqua(dir,_x_,_x_,jat,i) * gmat(_x_,_x_) + dqua(dir,_y_,_y_,jat,i) * gmat(_y_,_y_) &
                       + dqua(dir,_z_,_z_,jat,i) * gmat(_z_,_z_) &
                       + 2.0_rp*(dqua(dir,_x_,_y_,jat,i) * gmat(_x_,_y_) + dqua(dir,_x_,_z_,jat,i) * gmat(_x_,_z_) &
                       +         dqua(dir,_y_,_z_,jat,i) * gmat(_y_,_z_)))

                    hess(dir,dir2,atomi(jat),atomj(jat2)) = hess(dir,dir2,atomi(jat),atomj(jat2)) + y
                end do
            end do
        end do
    end do

end subroutine rotation_geomhess_pair8

subroutine rotate_multipoles(eel, ider, ddip, dqua, d2dip, d2qua)
    !! Rotates the atomic multipoles from the molecular frame, where they
    !! are defined as force field parameters, to the lab frame.
    !!
    !! If ider >= 1, also computes the first derivatives of the rotated
    !! dipole and quadrupole with respect to the positions of the (up to
    !! four) atoms that define each site's molecular frame: self, iz, ix,
    !! iy, in this order (iy is last as not every convention uses it; see
    !! [[rotation_matrix]]). If ider = 2, the second derivatives are
    !! computed too.
    !!
    !! Storage convention: the differentiation direction(s) always come
    !! first, so these arrays contract directly against a field/field
    !! gradient tensor, e.g. grad -= matmul(ddip(:,:,jat,j), E):
    !!   ddip(beta,alpha,jat,j)       = d [dip_j]_alpha / d r_{atom(jat)}^beta
    !!   dqua(beta,alpha,gamma,jat,j) = d [qua_j]_{alpha,gamma} / d r_{atom(jat)}^beta
    !!   d2dip(beta,delta,alpha,jat,jat2,j)
    !!       = d^2 [dip_j]_alpha / d r_{atom(jat)}^beta d r_{atom(jat2)}^delta
    !!   d2qua(beta,delta,alpha,gamma,jat,jat2,j)
    !!       = d^2 [qua_j]_{alpha,gamma} / d r_{atom(jat)}^beta d r_{atom(jat2)}^delta
    !! where jat, jat2 in [1,4] select self/iz/ix/iy (see _self_/_iz_/_ix_/
    !! _iy_ in f_cart_components.h).

    use mod_memory, only: ip, rp
    use mod_electrostatics, only: ommp_electrostatics_type
    use mod_io, only: fatal_error

    implicit none

    type(ommp_electrostatics_type), intent(inout) :: eel
    integer(ip), intent(in), optional :: ider
    real(rp), intent(out), optional :: ddip(3,3,4,eel%top%mm_atoms)
    real(rp), intent(out), optional :: dqua(3,3,3,4,eel%top%mm_atoms)
    real(rp), intent(out), optional :: d2dip(3,3,3,4,4,eel%top%mm_atoms)
    real(rp), intent(out), optional :: d2qua(3,3,3,3,4,4,eel%top%mm_atoms)

    integer(ip) :: j, jx, jy, jz, myider, nact
    integer(ip) :: jat, jat2, beta, delta
    real(rp), dimension(3,3) :: r, rt, qua, rqua, tmp
    real(rp), dimension(3,3,3,4) :: drmat
    real(rp), dimension(3,3,3,3,4,4) :: d2rmat
    real(rp), dimension(3) :: dip0
    real(rp), dimension(3,3) :: amat, bmat, cmat

    myider = 0
    if(present(ider)) myider = ider

    if(myider < 0 .or. myider > 2) &
        call fatal_error("rotate_multipoles: ider must be 0, 1 or 2")
    if(myider >= 1 .and. .not. (present(ddip) .and. present(dqua))) &
        call fatal_error("rotate_multipoles: ddip and dqua are required when ider >= 1")
    if(myider >= 2 .and. .not. (present(d2dip) .and. present(d2qua))) &
        call fatal_error("rotate_multipoles: d2dip and d2qua are required when ider >= 2")

    ! loop over the mm sites and build the rotation matrices.
    do j = 1, eel%top%mm_atoms
        nact = nref_atoms(eel%mol_frame(j))

        jz = eel%iz(j)
        if(jz == 0) jz = j
        jx = eel%ix(j)
        if(jx == 0) jx = j
        jy = eel%iy(j)
        if(jy == 0) jy = j

        call rotation_matrix(myider, &
                             eel%top%cmm(:,j), eel%top%cmm(:,jx), &
                             eel%top%cmm(:,jy), eel%top%cmm(:,jz), &
                             eel%mol_frame(j), &
                             r, drmat, d2rmat)

        rt = transpose(r)
        ! copy the monopole:
        eel%q(1,j) = eel%q0(1,j)

        ! rotate the dipole
        dip0 = eel%q0(2:4,j)
        eel%q(2:4,j) = matmul(r,dip0)

        ! extract, rotate and put back the quadrupole:
        qua(_x_,_x_)  = eel%q0(4+_xx_,j)
        qua(_x_,_y_)  = eel%q0(4+_xy_,j)
        qua(_x_,_z_)  = eel%q0(4+_xz_,j)
        qua(_y_,_x_)  = eel%q0(4+_yx_,j)
        qua(_y_,_y_)  = eel%q0(4+_yy_,j)
        qua(_y_,_z_)  = eel%q0(4+_yz_,j)
        qua(_z_,_x_)  = eel%q0(4+_zx_,j)
        qua(_z_,_y_)  = eel%q0(4+_zy_,j)
        qua(_z_,_z_)  = eel%q0(4+_zz_,j)

        tmp  = matmul(r,qua)
        rqua = matmul(tmp,rt)
        eel%q(4+_xx_,j)  = rqua(_x_,_x_)
        eel%q(4+_yy_,j)  = rqua(_y_,_y_)
        eel%q(4+_zz_,j)  = rqua(_z_,_z_)
        eel%q(4+_xy_,j)  = rqua(_x_,_y_)
        eel%q(4+_xz_,j)  = rqua(_x_,_z_)
        eel%q(4+_yz_,j)  = rqua(_y_,_z_)

        if(myider < 1) cycle

        ! first derivatives:
        !   d_beta dip = (d_beta R) dip0
        !   d_beta qua = (d_beta R) Q0 R^T - Q (d_beta R) R^T
        do jat = 1, nact
            do beta = 1, 3
                amat = drmat(:,:,beta,jat)
                ddip(beta,:,jat,j) = matmul(amat,dip0)
                dqua(beta,:,:,jat,j) = matmul(matmul(amat,qua),rt) &
                                     - matmul(rqua,matmul(amat,rt))
            end do
        end do

        if(myider < 2) cycle

        ! second derivatives:
        !   d2 dip = (d2 R) dip0
        !   d2 qua = (d2 R) Q0 R^T - Q (d2 R) R^T
        !          - (d_beta R) Q0 R^T (d_delta R) R^T
        !          - (d_delta R) Q0 R^T (d_beta R) R^T
        !          + Q (d_beta R) R^T (d_delta R) R^T
        !          + Q (d_delta R) R^T (d_beta R) R^T
        do jat2 = 1, nact
            do delta = 1, 3
                bmat = drmat(:,:,delta,jat2)
                do jat = 1, nact
                    do beta = 1, 3
                        amat = drmat(:,:,beta,jat)
                        cmat = d2rmat(:,:,beta,delta,jat,jat2)

                        d2dip(beta,delta,:,jat,jat2,j) = matmul(cmat,dip0)

                        d2qua(beta,delta,:,:,jat,jat2,j) = &
                              matmul(matmul(cmat,qua),rt) &
                            - matmul(rqua,matmul(cmat,rt)) &
                            - matmul(matmul(matmul(amat,qua),rt), matmul(bmat,rt)) &
                            - matmul(matmul(matmul(bmat,qua),rt), matmul(amat,rt)) &
                            + matmul(rqua, matmul(matmul(amat,rt), matmul(bmat,rt))) &
                            + matmul(rqua, matmul(matmul(bmat,rt), matmul(amat,rt)))
                    end do
                end do
            end do
        end do
    end do

end subroutine rotate_multipoles

subroutine rotation_matrix(ider,c,cx,cy,cz,mol_frame,r,dr,d2r)
    use mod_io,        only: fatal_error
    use mod_memory,    only: ip, rp
    use mod_constants, only: eps_rp

    !! given an atom j and the reference atoms jx, jy, and jz, this routine
    !! computes the rotation matrix needed to rotate the multipoles on the
    !! i-th atom from the molecular frame to the lab frame.
    !! if required, it also return the first and second derivatives of the
    !! rotation matrices with respect to the coordinates of all the atoms
    !! involved in its definition.
    !!
    !! this routine is completely general and can be easily augmented with
    !! new rotation conventions.
    !!
    !! given the atoms used to define the molecolar frame, identified by the
    !! indices jx, jy, jz, this routine builds the vectors
    !!
    !!   xi   = cmm(:,jz) - cmm(:,i)
    !!   eta  = cmm(:,jx) - cmm(:,i)
    !!   zeta = cmm(:,jy) - cmm(:,i)
    !!
    !! it then decodes the rotatoin conventions by introducing two vectors,
    !! u and v, that span the xz plane.
    !! this is the only convention-dependent part: everything else works
    !! automatically in a general way.
    !!
    !! for the definition of u and v, the unit vectors that identify the
    !! orthogonal molecular systems are built as follows:
    !!
    !! ez = u/|u|
    !! ex = gram-schmidt (v,ez)
    !! ey = ez x ex
    !!
    !! output:
    !! =======
    !!
    !! r(i,j) is the rotation matrix, whose columns are (ex,ey,ez)
    !!
    !! dr(i,j,k,jat) contains the derivative of the i-th component of e_j
    !! (i.e. of r(i,j)) with respect to the k-th cartesian component of
    !! atom jat = self, iz, ix, iy (see _self_/_iz_/_ix_/_iy_ in
    !! f_cart_components.h). Note the column-select index j comes before
    !! the derivative direction k, matching how it is populated below
    !! (e.g. dr(:,1,k,jtype) = ex_r(:,k,jtype)).
    !!
    !! d2r(i,j,k,l,jat,jat2) contains the second derivative of the i-th
    !! component of e_j with respect to the k-th component of atom jat
    !! and the l-th component of atom jat2, with the same jat/jat2
    !! convention as above.
    !!

    implicit none

    integer(ip),                         intent(in)    :: ider
    real(rp),    dimension(3),           intent(in)    :: c, cx, cy, cz
    integer(ip),                         intent(in)    :: mol_frame
    real(rp),    dimension(3,3),         intent(inout) :: r
    real(rp),    dimension(3,3,3,4),     intent(inout) :: dr
    real(rp),    dimension(3,3,3,3,4,4), intent(inout) :: d2r

    integer(ip) :: k, a, b, jtype, jtype2
    real(rp) :: dot
    real(rp) :: xi_norm, eta_norm, zeta_norm
    real(rp), parameter :: zofac = 0.866_rp

    real(rp), dimension(3)         :: xi, eta, zeta
    real(rp), dimension(3)         :: xihat, etahat, zetahat
    real(rp), dimension(3)         :: u, v, w, ex, ey, ez
    real(rp), dimension(3,3)       :: xi_xi, eta_eta, zeta_zeta
    real(rp), dimension(3,3,3)     :: xi_xixi, eta_etaeta, zeta_zetazeta
    real(rp), dimension(3,3)       :: ez_u, ex_w, w_v, w_ez
    real(rp), dimension(3,3,3)     :: ez_uu, ex_ww, w_ezez, w_ezv
    real(rp), dimension(3,3,4)     :: u_r, v_r, ez_r, w_r, ex_r, ey_r
    real(rp), dimension(3,3,4,3,4) :: u_rr, v_rr, ez_rr, w_rr, ex_rr, ey_rr

    if (ider.lt.0 .or. ider.gt.2) then
        call fatal_error("wrong ider in rotation_matrix2")
    end if

    r = 0.0_rp
    if (ider.ge.1) dr  = 0.0_rp
    if (ider.ge.2) d2r = 0.0_rp

    xi = cz - c
    xi_norm = sqrt(dot_product(xi,xi))

    eta = cx - c
    eta_norm = sqrt(dot_product(eta,eta))

    zeta = cy - c
    zeta_norm = sqrt(dot_product(zeta,zeta))

    u = 0.0_rp
    v = 0.0_rp

!
!   initialize the unit vectors derivatives to zero.
!
    eta_eta       = 0.0_rp
    eta_etaeta    = 0.0_rp
    zeta_zeta     = 0.0_rp
    zeta_zetazeta = 0.0_rp

    if (mol_frame.eq.0) then

        r(1,1) = 1.0_rp
        r(2,2) = 1.0_rp
        r(3,3) = 1.0_rp
        return

    else if (mol_frame.eq.1) then

        ! z-then-x convention
        call unit_vector_derivatives(ider,xi,xihat,xi_xi,xi_xixi)
        call unit_vector_derivatives(ider,eta,etahat,eta_eta,eta_etaeta)
        u = xi
        v = eta

    else if (mol_frame.eq.2) then

        ! bisector convention, normalized form
        call unit_vector_derivatives(ider,xi,xihat,xi_xi,xi_xixi)
        call unit_vector_derivatives(ider,eta,etahat,eta_eta,eta_etaeta)
        u = xihat + etahat
        v = eta

    else if (mol_frame.eq.3) then

        ! z-only convention
        u = xi
        dot = u(3)/xi_norm
        if (dot.le.zofac .and. abs(u(2)) > eps_rp) then
            v(1) = 1.0_rp
        else
            v(2) = 1.0_rp
        end if
        call unit_vector_derivatives(ider,xi,xihat,xi_xi,xi_xixi)

    else if (mol_frame.eq.4) then

        ! z-bisector convention, normalized form
        call unit_vector_derivatives(ider,xi,xihat,xi_xi,xi_xixi)
        call unit_vector_derivatives(ider,eta,etahat,eta_eta,eta_etaeta)
        call unit_vector_derivatives(ider,zeta,zetahat,zeta_zeta,zeta_zetazeta)
        u = xi
        v = etahat + zetahat

    else if (mol_frame.eq.5) then

        ! 3-fold convention, normalized form
        call unit_vector_derivatives(ider,xi,xihat,xi_xi,xi_xixi)
        call unit_vector_derivatives(ider,zeta,zetahat,zeta_zeta,zeta_zetazeta)
        call unit_vector_derivatives(ider,eta,etahat,eta_eta,eta_etaeta)
        u = xihat + etahat + zetahat
        v = eta

    else

        call fatal_error('the required rotation convention is not implemented.')

    end if

!
!   Build ez = u/|u|, w = v - (v.ez) ez, ex = w/|w|, ey = ez x ex.
!
    call unit_vector_derivatives(ider,u,ez,ez_u,ez_uu)
    call w_derivatives(ider,ez,v,w,w_ez,w_v,w_ezez,w_ezv)
    call unit_vector_derivatives(ider,w,ex,ex_w,ex_ww)

    ey = cross_product(ez,ex)

    r(:,1) = ex
    r(:,2) = ey
    r(:,3) = ez

    if (ider.eq.0) return
!
!   Convention-specific first and second derivatives of u and v.
!
    u_r  = 0.0_rp
    v_r  = 0.0_rp
    u_rr = 0.0_rp
    v_rr = 0.0_rp

    if (mol_frame.eq.1) then

        ! z-then-x: u = xi, v = eta
        do k = 1, 3
            u_r(k,k,1) = -1.0_rp
            u_r(k,k,2) =  1.0_rp

            v_r(k,k,1) = -1.0_rp
            v_r(k,k,3) =  1.0_rp
        end do

    else if (mol_frame.eq.2) then

        ! bisector: u = xihat + etahat, v = eta
        call add_unit_vector_chain(2, xi_xi,  xi_xixi,  u_r, u_rr)
        call add_unit_vector_chain(3, eta_eta,eta_etaeta,u_r, u_rr)

        do k = 1, 3
            v_r(k,k,1) = -1.0_rp
            v_r(k,k,3) =  1.0_rp
        end do

    else if (mol_frame.eq.3) then

        ! z-only: u = xi, v = fixed Cartesian axis
        do k = 1, 3
            u_r(k,k,1) = -1.0_rp
            u_r(k,k,2) =  1.0_rp
        end do

    else if (mol_frame.eq.4) then

        ! z-bisector: u = xi, v = etahat + zetahat
        do k = 1, 3
            u_r(k,k,1) = -1.0_rp
            u_r(k,k,2) =  1.0_rp
        end do

        call add_unit_vector_chain(3, eta_eta,  eta_etaeta,  v_r, v_rr)
        call add_unit_vector_chain(4, zeta_zeta,zeta_zetazeta,v_r, v_rr)

    else if (mol_frame.eq.5) then

        ! 3-fold: u = xihat + etahat + zetahat, v = eta
        call add_unit_vector_chain(2, xi_xi,    xi_xixi,       u_r, u_rr)
        call add_unit_vector_chain(3, eta_eta,  eta_etaeta,    u_r, u_rr)
        call add_unit_vector_chain(4, zeta_zeta,zeta_zetazeta, u_r, u_rr)

        do k = 1, 3
            v_r(k,k,1) = -1.0_rp
            v_r(k,k,3) =  1.0_rp
        end do

    end if

!
!   First derivatives through the chain:
!   u -> ez, (v,ez) -> w, w -> ex, (ez,ex) -> ey.
!
    do jtype = 1, 4

        ez_r(:,:,jtype) = matmul(ez_u,u_r(:,:,jtype))

        w_r(:,:,jtype) = matmul(w_v, v_r(:,:,jtype)) &
                       + matmul(w_ez,ez_r(:,:,jtype))

        ex_r(:,:,jtype) = matmul(ex_w,w_r(:,:,jtype))

        do k = 1, 3
            ey_r(:,k,jtype) = cross_product(ez_r(:,k,jtype),ex) &
                            + cross_product(ez,ex_r(:,k,jtype))
        end do

    end do

!
!   Store first derivatives.
!
    do jtype = 1, 4
        do k = 1, 3
            dr(:,1,k,jtype) = ex_r(:,k,jtype)
            dr(:,2,k,jtype) = ey_r(:,k,jtype)
            dr(:,3,k,jtype) = ez_r(:,k,jtype)
        end do
!       dr(:,:,1,jtype) = transpose(ex_r(:,:,jtype))
!       dr(:,:,2,jtype) = transpose(ey_r(:,:,jtype))
!       dr(:,:,3,jtype) = transpose(ez_r(:,:,jtype))
    end do

    if (ider.eq.1) return

!
!   Second derivatives through the same chain.
!
    do jtype = 1, 4
        do a = 1, 3
            do jtype2 = 1, 4
                do b = 1, 3

!
!                   ez = ez(u)
!
                    ez_rr(:,a,jtype,b,jtype2) = &
                          matmul(ez_u, u_rr(:,a,jtype,b,jtype2)) &
                        + hess_vec(ez_uu, u_r(:,a,jtype), u_r(:,b,jtype2))

!
!                   w = w(v,ez)
!
                    w_rr(:,a,jtype,b,jtype2) = &
                          matmul(w_v,  v_rr(:,a,jtype,b,jtype2)) &
                        + matmul(w_ez, ez_rr(:,a,jtype,b,jtype2)) &
                        + hess_vec(w_ezez, ez_r(:,a,jtype), ez_r(:,b,jtype2)) &
                        + hess_vec(w_ezv,  ez_r(:,a,jtype), v_r(:,b,jtype2)) &
                        + hess_vec_ezv_swap(w_ezv, v_r(:,a,jtype), ez_r(:,b,jtype2))

!
!                   ex = ex(w)
!
                    ex_rr(:,a,jtype,b,jtype2) = &
                          matmul(ex_w, w_rr(:,a,jtype,b,jtype2)) &
                        + hess_vec(ex_ww, w_r(:,a,jtype), w_r(:,b,jtype2))

!
!                   ey = ez x ex
!
                    ey_rr(:,a,jtype,b,jtype2) = &
                          cross_product(ez_rr(:,a,jtype,b,jtype2), ex) &
                        + cross_product(ez_r(:,a,jtype), ex_r(:,b,jtype2)) &
                        + cross_product(ez_r(:,b,jtype2), ex_r(:,a,jtype)) &
                        + cross_product(ez, ex_rr(:,a,jtype,b,jtype2))

                end do
            end do
        end do
    end do

!
!   Store second derivatives.
!
    do jtype = 1, 4
        do jtype2 = 1, 4
            do a = 1, 3
                do b = 1, 3
                    d2r(:,1,a,b,jtype,jtype2) = ex_rr(:,a,jtype,b,jtype2)
                    d2r(:,2,a,b,jtype,jtype2) = ey_rr(:,a,jtype,b,jtype2)
                    d2r(:,3,a,b,jtype,jtype2) = ez_rr(:,a,jtype,b,jtype2)
                end do
            end do
        end do
    end do

contains

    pure function cross_product(a,b) result(c)
        implicit none
        real(rp), intent(in) :: a(3), b(3)
        real(rp)             :: c(3)

        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function cross_product

    pure function hess_vec(h,x,y) result(z)
        implicit none
        real(rp), intent(in) :: h(3,3,3), x(3), y(3)
        real(rp) :: z(3)
        integer(ip) :: aa, bb, cc

        z = 0.0_rp
        do aa = 1, 3
            do bb = 1, 3
                do cc = 1, 3
                    z(aa) = z(aa) + h(aa,bb,cc)*x(bb)*y(cc)
                end do
            end do
        end do
    end function hess_vec

    pure function hess_vec_ezv_swap(h,x,y) result(z)
        implicit none
        real(rp), intent(in) :: h(3,3,3), x(3), y(3)
        real(rp) :: z(3)
        integer(ip) :: aa, bb, cc

!
!       h stores d2w_a / d ez_b d v_c.
!       This helper contracts d2w_a / d v_b d ez_c
!       using h(a,c,b).
!
        z = 0.0_rp
        do aa = 1, 3
            do bb = 1, 3
                do cc = 1, 3
                    z(aa) = z(aa) + h(aa,cc,bb)*x(bb)*y(cc)
                end do
            end do
        end do
    end function hess_vec_ezv_swap

end subroutine rotation_matrix
!
subroutine unit_vector_derivatives(ider, v, e, e_v, e_vv)
    use mod_memory,    only: ip, rp
    use mod_io,        only: fatal_error
    use mod_constants, only: eps_rp

    implicit none

    integer(ip), intent(in)  :: ider
    real(rp),    intent(in)  :: v(3)
    real(rp),    intent(out) :: e(3)
    real(rp),    intent(out) :: e_v(3,3)
    real(rp),    intent(out) :: e_vv(3,3,3)

    integer(ip) :: a, b, c
    real(rp) :: v_norm, v_norm2

    if (ider.lt.0 .or. ider.gt.2) then
        call fatal_error("wrong ider in unit_vector_derivatives")
    end if

    e    = 0.0_rp
    e_v  = 0.0_rp
    e_vv = 0.0_rp

    v_norm = sqrt(dot_product(v,v))

    if (v_norm < eps_rp) then
        call fatal_error("unit_vector_derivatives: zero norm vector.")
    end if

    v_norm2 = v_norm*v_norm
    e = v / v_norm
!
!   done for the matrix only.
!
    if (ider.eq.0) return

    do a = 1, 3
        e_v(a,a) = 1.0_rp
        do b = 1, 3
            e_v(a,b) = e_v(a,b) - e(a)*e(b)
        end do
    end do
    e_v = e_v/v_norm
!
!   done for first derivatives only.
!
    if (ider.eq.1) return

    do a = 1, 3
        do b = 1, 3
            e_vv(a,a,b) = e_vv(a,a,b) - e(b)
            e_vv(a,b,a) = e_vv(a,b,a) - e(b)
            e_vv(b,a,a) = e_vv(b,a,a) - e(b)
            do c = 1, 3
                e_vv(a,b,c) = e_vv(a,b,c) + 3.0_rp*e(a)*e(b)*e(c)
            end do
        end do
    end do
    e_vv = e_vv / v_norm2

end subroutine unit_vector_derivatives
!
subroutine w_derivatives(ider, ez, v, w, w_ez, w_v, w_ezez, w_ezv)
    use mod_memory, only: ip, rp
    implicit none

    integer(ip), intent(in)  :: ider
    real(rp),    intent(in)  :: ez(3), v(3)
    real(rp),    intent(out) :: w(3)
    real(rp),    intent(out) :: w_ez(3,3)
    real(rp),    intent(out) :: w_v(3,3)
    real(rp),    intent(out) :: w_ezez(3,3,3)
    real(rp),    intent(out) :: w_ezv(3,3,3)

    integer(ip) :: a, b, c
    real(rp) :: vez

    vez = dot_product(v, ez)

    w = v - vez*ez

    w_ez   = 0.0_rp
    w_v    = 0.0_rp
    w_ezez = 0.0_rp
    w_ezv  = 0.0_rp
!
    if (ider.eq.0) return
!
!   First derivatives:
!
!   w_a = v_a - (v.ez) ez_a
!
!   d w_a / d ez_b = - ez_a v_b - (v.ez) delta_ab
!   d w_a / d v_b  =   delta_ab - ez_a ez_b
!
    do a = 1, 3
        w_ez(a,a) = - vez
        w_v(a,a) = 1.0_rp
        do b = 1, 3
            w_ez(a,b) = w_ez(a,b) - ez(a)*v(b)
            w_v(a,b)  = w_v(a,b)  - ez(a)*ez(b)
        end do
    end do
!
    if (ider.eq.1) return
!
!   Second derivatives:
!
!   d2 w_a / d ez_b d ez_c = -v_b delta_ac - v_c delta_ab
!
!   d2 w_a / d ez_b d v_c  = -delta_bc ez_a - ez_c delta_ab
!
!   d2 w_a / d v_b d v_c = 0, not stored.
!
    do a = 1, 3
        do b = 1, 3
            w_ezez(a,b,a) = w_ezez(a,b,a) - v(b)
            w_ezez(a,a,b) = w_ezez(a,a,b) - v(b)
            w_ezv(a,b,b)  = w_ezv(a,b,b)  - ez(a)
            w_ezv(a,a,b)  = w_ezv(a,a,b)  - ez(b)
        end do
    end do

end subroutine w_derivatives

subroutine add_unit_vector_chain(jref, e_x, e_xx, f_r, f_rr)
    use mod_memory, only: ip, rp
    implicit none

    integer(ip), intent(in) :: jref
    real(rp), intent(in)    :: e_x(3,3)
    real(rp), intent(in)    :: e_xx(3,3,3)
    real(rp), intent(inout) :: f_r(3,3,4)
    real(rp), intent(inout) :: f_rr(3,3,4,3,4)

    integer(ip) :: aa, bb, p, q
    real(rp) :: sp, sq

    do p = 1, 4

        sp = 0.0_rp
        if (p.eq.1)    sp = -1.0_rp
        if (p.eq.jref) sp =  1.0_rp

        if (sp.ne.0.0_rp) then
            do aa = 1, 3
                f_r(:,aa,p) = f_r(:,aa,p) + sp*e_x(:,aa)
            end do
        end if

        do q = 1, 4

            sq = 0.0_rp
            if (q.eq.1)    sq = -1.0_rp
            if (q.eq.jref) sq =  1.0_rp

            if (sp.ne.0.0_rp .and. sq.ne.0.0_rp) then
                do aa = 1, 3
                    do bb = 1, 3
                        f_rr(:,aa,p,bb,q) = f_rr(:,aa,p,bb,q) &
                                          + sp*sq*e_xx(:,aa,bb)
                    end do
                end do
            end if

        end do
    end do

end subroutine add_unit_vector_chain

end module mod_rotate_multipoles
