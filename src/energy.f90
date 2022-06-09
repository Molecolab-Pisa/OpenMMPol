subroutine energy(scr,EMM)
    use mod_mmpol, only: q, v_qq, ef_qd, amoeba, ipd, mm_atoms, pol_atoms
    use mod_memory, only : ip, rp
    use mod_constants, only: zero, pt5, two
    implicit none
    !
    ! Compute the total MM interaction energy between static multipole
    ! distributions (MM-MM) and between static multipole distribution
    ! and induced dipoles (MM-MMpol). Which one of these two energies is 
    ! calculated is controlled by scr:
    ! 
    ! scr  ...  0:  The electrostatic interaction energy between static
    !               multipole distributions on MM atoms is calculated (MM-MM)  
    !           1:  The electrostatic interaction energy between static 
    !               multipole distribution of MM atoms and induced dipoles
    !               on MMpol atoms (MM-MMpol)
    !
    real(rp), intent(out) :: EMM
    integer(ip), intent(in) :: scr
    !
    integer(ip) :: I
    real(rp)    :: ddot
    real(rp)    :: qq, px, py, pz, qxx, qxy, qxz, qyy, qyz, qzz
    real(rp)    :: fv, ex, ey, ez, gxx, gxy, gxz, gyy, gyz, gzz
    !
    ! real(rp), dimension(10), parameter :: fac = (/ 1.0_rp, -1.0_rp, -1.0_rp, -1.0_rp, 1.0_rp, &
    !                                                2.0_rp, 1.0_rp, 2.0_rp, 2.0_rp, 1.0_rp /)
    !
    
    ! Initialize the interaction energy
    EMM = Zero
    
    ! Compute intercation energy
    if (scr.eq.0) then          ! MM-MM interaction energy
        !
        ! Faster and more compact approach
        !
    
        !Nqq =  SIZE(q, 1)
        !do I= 1, Nqq
        !    EMM = EMM + fac(I)*ddot(mm_atoms,v_qq(I,:),1,q( I,:),1)*pt5
        !enddo
    
        !
        ! More human readable approach
        !
        if (amoeba) then        ! AMOEBA force-field
            do I = 1, mm_atoms
                !
                ! initialize multipoles
                !
                qq  = q( 1,I)
                px  = q( 2,I)
                py  = q( 3,I)
                pz  = q( 4,I)
                qxx = q( 5,I)
                qxy = q( 6,I)
                qyy = q( 7,I)
                qxz = q( 8,I)
                qyz = q( 9,I)
                qzz = q(10,I)
                !
                ! Initialize potential, field and field gradients
                !
                fv  = v_qq(1,I)
                ex  = v_qq(2,I)
                ey  = v_qq(3,I)
                ez  = v_qq(4,I)
                gxx = v_qq(5,I)
                gxy = v_qq(6,I)
                gyy = v_qq(7,I)
                gxz = v_qq(8,I)
                gyz = v_qq(9,I)
                gzz = v_qq(10,I)
                !
                ! Compute interaction energy
                !
                EMM = EMM + qq*fv - (px*ex + py*ey + pz*ez) + qxx*gxx + qyy*gyy + qzz*gzz  + two*(qxy*gxy + qxz*gxz + qyz*gyz)
            enddo
        else                     ! ABMER force-field
            do I = int(1,ip), mm_atoms
                !
                ! initialize multipoles
                !
                qq  = q( int(1,ip),I)
                !
                ! Initialize potential, field and field gradients
                !
                fv  = v_qq(int(1,ip),I)
                !
                ! Compute interaction energy
                !
                EMM = EMM + qq*fv 
            enddo
        end if
        EMM = EMM*pt5        
        
        
    elseif (scr.eq.1) then      ! MM-MMpol interaction energy
        
        if (Amoeba) then        ! For Amoeba interaction energy is d-ipd*ep
            do I=1,3
                !EMM = EMM - ddot(pol_atoms,ef_qd(I,:,1),1,ipd(I,:,2),1)*pt5   ! factor 1/2 because we have interaction with
                                                                               ! induced dipoles
                EMM = EMM - ddot(pol_atoms,ef_qd(I,:,2),1,ipd(I,:,1),1)*pt5    ! factor 1/2 because we have interaction with
                                                                               ! induced dipoles
            enddo
        else                    ! AMBER polarizable forcefield
            do I=1,3
                EMM = EMM - ddot(pol_atoms,ef_qd(I,:,1),1,ipd(I,:,1),1)*pt5   ! factor 1/2 because we have interaction with
                                                                             ! induced dipoles
            enddo
        end if
    end if
    
end subroutine
