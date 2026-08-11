#include "f_cart_components.h"

module mod_geomhess
    use mod_io, only: fatal_error, ommp_message
    use mod_memory, only: ip, rp
    use mod_mmpol, only: ommp_system
    use mod_topology, only: ommp_topology_type
    use mod_profiling, only: time_push, time_pull

    implicit none
    private
    
    public :: fixedelec_geomhess, polelec_geomhess

    contains

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
!           !use mod_electrostatics, only: prepare_M2D, ommp_electrostatics_type
            use mod_polarization, only: polarization
            use mod_electrostatics

            implicit none
!           
            type(ommp_system), intent(inout), target :: s
            !! System data structure
            real(rp), dimension(3,s%top%mm_atoms,3,s%top%mm_atoms), intent(inout) :: hess
            !! Geometrical Hessian in output, results will be added
!           
!           integer(ip) :: i
!           type(ommp_electrostatics_type), pointer :: eel 
!           eel => s%eel

        end subroutine
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
