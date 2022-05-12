!-----------------------------------------------------------
!   Not finished
!-----------------------------------------------------------
module force
contains


subroutine forces_MM(dv,ff_qq)
    use mmpol
    implicit none
    !
    ! Computes forces of the static multipole J at the position  
    ! of MM atom I. The computed forces are added to the forces ff_qq.
    ! The derivatives of the potential for the calculation of the forces
    ! are specified in dv, which is organized as follows:
    !
    ! - For AMBER force-field only potential derivative at mm sites are needed
    !
    !       ex, ey, ez
    !
    ! - For AMOEBA force-field derivatives of the potential, electric 
    !   field and field gradients at MM sites are needed
    !   
    !       ex, ey, ez, gxx, gxy, gyy, gxz, gyz, gzz, hxxx, hxxy, hxxz,  
    !       hxyy, hxyz, hxzz, hyyy, hyyz, hyzz, hzzz
    !
    !
    real(rp), dimension(ld_cder,mm_atoms), intent(in) :: dv
    real(rp), dimension(3,mm_atoms), intent(inout) :: ff_qq
    !
    integer(ip) :: I
    real(rp)    :: qq, px, py, pz, qxx, qxy, qyy, qxz, qyz, qzz
    real(rp)    :: ex, ey, ez, gxx, gxy, gyy, gxz, gyz, gzz, hxxx, &
                   hxxy, hxxz, hxyy, hxyz, hxzz, hyyy, hyyz, hyzz, hzzz
    !
    If(.not.Amoeba) then  ! AMBER force-field
        do I = 1,mm_atoms
            !
            ! Intermediates
            qq   = q( 1,I)
            ex   = dv( 1,I)
            ey   = dv( 2,I)
            ez   = dv( 3,I)
            !
            ! Compute forces
            ff_qq(1,I) = -qq*ex
            ff_qq(2,I) = -qq*ey
            ff_qq(3,I) = -qq*ez
        enddo
    else                ! AMOEBA force-field
        do I = 1,mm_atoms
            !
            ! Intermediates - multipoles
            qq   = q( 1,I)
            px   = q( 2,I)
            py   = q( 3,I)
            pz   = q( 4,I)
            qxx  = q( 5,I)
            qxy  = q( 6,I)
            qyy  = q( 7,I)
            qxz  = q( 8,I)
            qyz  = q( 9,I)
            qzz  = q(10,I)
            !
            ! Intermediates - potential derivative
            ex   = dv( 1,I)
            ey   = dv( 2,I)
            ez   = dv( 3,I)
            gxx  = dv( 4,I)
            gxy  = dv( 5,I)
            gyy  = dv( 6,I)
            gxz  = dv( 7,I)
            gyz  = dv( 8,I)
            gzz  = dv( 9,I)
            hxxx = dv(10,I)
            hxxy = dv(11,I)
            hxxz = dv(12,I)
            hxyy = dv(13,I)
            hxyz = dv(14,I)
            hxzz = dv(15,I)
            hyyy = dv(16,I)
            hyyz = dv(17,I)
            hyzz = dv(18,I)
            hzzz = dv(19,I)
            !
            ! Compute forces
            ff_qq(1,I) = -(qq*ex - gxx*px - gxy*py - gxz*pz + hxxx*qxx + hxyy*qyy + hxzz*qzz + &
                      two*(hxxy*qxy + hxyz*qyz + hxxz*qxz))
            ff_qq(2,I) = -(qq*ey - gxy*px - gyy*py - gyz*pz + hxxy*qxx + hyyy*qyy + hyzz*qzz + &
                      two*(hxyy*qxy + hxyz*qxz + hyyz*qyz))
            ff_qq(3,I) = -(qq*ez - gxz*px - gyz*py - gzz*pz + hxxz*qxx + hyyz*qyy + hzzz*qzz + &
                      two*(hxyz*qxy + hyzz*qyz + hxzz*qxz))
        enddo
    end if
end subroutine forces_MM

subroutine forces_DD(def,ff_dd) 
    use mmpol
    implicit none
    !
    ! Computes forces of the static multipole J at the position  
    ! of MM atom I. The computed forces are added to the forces ff_qq.
    ! The derivatives of the electric field for the calculation of the 
    ! forces are specified in def, which are organized as follows:
    !
    ! - For AMBER force-field electric field gradients at polarizable
    !   site are needed in following order
    !
    !       gxx, gxy, gyy, gxz, gyz, gzz
    !
    ! - For AMOEBA force-field electric field gradients for the direct
    !   and polarization field are needed in following order
    !
    !       gxx, gxy, gyy, gxz, gyz, gzz
    !
    !   the direct field gradients (for d dipoles) are in de(:,:,1)
    !   the polarization field gradients (for p dipoles) are in de(:,:,2)
    !
    real(rp), dimension(3,pol_atoms), intent(inout) :: ff_dd
    real(rp), dimension(6,pol_atoms,n_ipd), intent(in) :: def
    !
    integer(ip) :: I
    real(rp)    :: pdx, pdy, pdz, ppx, ppy, ppz
    real(rp)    :: gdxx, gdxy, gdyy, gdxz, gdyz, gdzz, &
                   gpxx, gpxy, gpyy, gpxz, gpyz, gpzz
    !
    if (.not.Amoeba) then
        do I=1, pol_atoms
            !
            ! Intermediates - field gradients
            gdxx = def(1,I,1)
            gdxy = def(2,I,1)
            gdyy = def(3,I,1)
            gdxz = def(4,I,1)
            gdyz = def(5,I,1)
            gdzz = def(6,I,1)
            !
            ! Intermediates - dipoles
            pdx = ipd(1,I,1)
            pdy = ipd(2,I,1)
            pdz = ipd(3,I,1)
            !
            ! Compute forces
            ff_dd(1,I) = (pdx*gdxx + pdy*gdxy + pdz*gdxz)
            ff_dd(2,I) = (pdx*gdxy + pdy*gdyy + pdz*gdyz)
            ff_dd(3,I) = (pdx*gdxz + pdy*gdyz + pdz*gdzz)
        enddo
    else
        do I=1, pol_atoms
            !
            ! Intermediates - field gradients
            gdxx = def(1,I,1)
            gdxy = def(2,I,1)
            gdyy = def(3,I,1)
            gdxz = def(4,I,1)
            gdyz = def(5,I,1)
            gdzz = def(6,I,1)
            gpxx = def(1,I,2)
            gpxy = def(2,I,2)
            gpyy = def(3,I,2)
            gpxz = def(4,I,2)
            gpyz = def(5,I,2)
            gpzz = def(6,I,2)
            !
            ! Intermediates - dipoles
            pdx  = ipd(1,I,1)
            pdy  = ipd(2,I,1)
            pdz  = ipd(3,I,1)
            ppx  = ipd(1,I,2)
            ppy  = ipd(2,I,2)
            ppz  = ipd(3,I,2)
            !
            ! Compute forces
            ff_dd(1,I) = pt5*(gdxx*ppx + gdxy*ppy + gdxz*ppz + gpxx*pdx + gpxy*pdy + gpxz*pdz)
            ff_dd(2,I) = pt5*(gdxy*ppx + gdyy*ppy + gdyz*ppz + gpxy*pdx + gpyy*pdy + gpyz*pdz)
            ff_dd(3,I) = pt5*(gdxz*ppx + gdyz*ppy + gdzz*ppz + gpxz*pdx + gpyz*pdy + gpzz*pdz)
        enddo
    end if
end subroutine forces_DD
    
end module force

