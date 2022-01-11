!-----------------------------------------------------------
!   Not finished
!-----------------------------------------------------------
module elstat
contains
subroutine multipoles_potential(scr,v)
    use mmpol
    implicit none
    !
    ! Computes potential, electric field and field gradients at position 
    ! of MM atoms. The computed potential is added to the potential v.
    !
    ! scr  ...  0:  The sources are static multipoles of mm atoms  
    !               contained in the array q and with coordinates cmm
    !           1:  The sources are induced multipoles ipd, contained
    !               in the array ipd with coordinates cpol
    ! 
    ! - For AMBER force-field only potential at MM sites is computed.
    !   The output is given in following order:
    !
    !       v
    !
    ! - For AMOEBA force-field potential, electric field and field gradients
    !   are computed. The output is given in following order:
    !
    !       v, ex, ey, ez, Gxx, Gxy, Gyy, Gxz, Gyz, Gzz 
    !
    
    real(rp), dimension(ld_cart,mm_atoms), intent(inout) :: v
    integer(ip), intent(in) :: scr
    !
    !logical     :: Amoeba
    integer(ip) :: I, J
    !
    !real(rp), parameter :: One = 1.0_rp
    
    
    !Amoeba = (ld_cart.gt.1)

    if (scr.eq.0) then
        do I = 1,mm_atoms
            do J = 1,mm_atoms
                if (I.eq.J) cycle

                call potential_M2M(One,I,J,v) 

            enddo
        enddo
    elseif (scr.eq.1) then
        do I = 1,mm_atoms
            do J = 1,pol_atoms
                if (polar_mm(J).eq.I) cycle

                call potential_D2M(Zero,One,I,J,v) 

            enddo
        enddo
    end if
    
    ! Remove unwanted contributions from the bonded atoms
    call multipoles_potential_remove(scr,v)
        
end subroutine multipoles_potential

subroutine multipoles_field(scr,e)
    use mmpol
    implicit none
    !
    ! Computes electric field at polarizable sites cpol. The computed
    ! field is added to the field e.
    !
    ! scr  ...  0:  The sources are static multipoles of mm atoms  
    !               contained in the array q and with coordinates cmm
    !           1:  The sources are induced multipoles ipd, contained
    !               in the array ipd with coordinates cpol.
    ! 
    ! - For AMBER force-field only electric field at polarizable site
    !   is computed. The output is given in following order:
    !
    !       ex, ey,ez
    !
    ! - For AMOEBA force-field electric fields for p and d dipoles
    !   are computed. The output is given in following order:
    !
    !       ex, ey,ez
    !
    !   the direct field (for d dipoles) is in e(:,:,1)
    !   the polarization field (for p dipoles) is in e(:,:,2)
    !
    ! - Wang's polynomial damping for of the coulomb potential is used for
    !   the electric field generated at polarizable site for both AMOEBA 
    !   and AMBER polarizable force-fields 
    !

    
    real(rp), dimension(3,pol_atoms,n_ipd), intent(inout) :: e
    integer(ip), intent(in) :: scr
    !
    !logical     :: Amoeba
    integer(ip) :: I, J
    !real(rp), parameter :: one = 1.0_rp
    
    !Amoeba = n_ipd.eq.2
    
    if (scr.eq.0) then
        do I = 1,pol_atoms
            do J = 1,mm_atoms
                if (polar_mm(I).eq.J) cycle

                call field_M2D(one,one,I,J,e)

            enddo
        enddo
    elseif (scr.eq.1) then
        do I = 1,pol_atoms
            do J = 1,pol_atoms
                if (I.eq.J) cycle
                
                call field_D2D(one,one,I,J,e)
            enddo
        enddo
    end if
    
    ! Remove unwanted contributions from the bonded atoms
    call multipoles_field_remove(scr,e)
    
end subroutine multipoles_field

subroutine multipoles_potential_deriv(scr,dv)
    use mmpol
    implicit none
    !
    ! Computes derivatives of the potential, electric field and field 
    ! gradients at position of MM atoms. The computed derivatives of
    ! the potential are added to the dv. 
    !
    ! scr  ...  0:  The sources are static multipoles of mm atoms  
    !               contained in the array q and with coordinates cmm
    !           1:  The sources are induced multipoles ipd, contained
    !               in the array ipd with coordinates cpol
    ! 
    ! - For AMBER force-field only potential derivative at mm sites is 
    !   computed. The output is given in following order:
    !
    !       ex, ey, ez
    !
    ! - For AMOEBA force-field derivatives of the potential, electric 
    !   field and field gradients are computed. The output is given in
    !   following order:
    !   
    !       ex, ey, ez, gxx, gxy, gyy, gxz, gyz, gzz, hxxx, hxxy, hxxz,  
    !       hxyy, hxyz, hxzz, hyyy, hyyz, hyzz, hzzz
    !
    ! - For AMOEBA force-field, the scaling of the nearest neighbor 
    !   interactions for p and d dipoles is reversed dscale is used for 
    !   scaling potential of p dipoles and pscale is used for scaling 
    !   potential of d dipoles.
    !
    
    real(rp), dimension(ld_cder,mm_atoms), intent(inout) :: dv
    integer(ip), intent(in) :: scr
    !
    !logical     :: Amoeba
    integer(ip) :: I, J
    !
    !real(rp), parameter :: One = 1.0_rp
    
    
    !Amoeba = (ld_cder.gt.3)

    if (scr.eq.0) then
        do I = 1,mm_atoms
            do J = 1,mm_atoms
                if (I.eq.J) cycle

                call potential_deriv_M2M(One,I,J,dv) 

            enddo
        enddo
    elseif (scr.eq.1) then
        do I = 1,mm_atoms
            do J = 1,pol_atoms
                if (polar_mm(J).eq.I) cycle

                call potential_deriv_D2M(One,One,I,J,dv) 
                
            enddo
        enddo
    end if
    
    ! Remove unwanted contributions from the bonded atoms
    call multipoles_potential_deriv_remove(scr,dv)
    
end subroutine multipoles_potential_deriv

subroutine multipoles_field_deriv(scr,de)
    use mmpol
    implicit none
    !
    ! Computes derivatives of electric field at polarizable sites cpol.
    !  The computed field derivatives are added to the de.
    !
    ! scr  ...  0:  The sources are static multipoles of mm atoms  
    !               contained in the array q and with coordinates cmm
    !           1:  The sources are induced multipoles ipd, contained
    !               in the array ipd with coordinates cpol.
    ! 
    ! - For AMBER force-field electric field gradients at polarizable
    !   site are computed. The output is given in following order:
    !
    !       gxx, gxy, gyy, gxz, gyz, gzz
    !
    ! - For AMOEBA force-field electric field gradients for the direct
    !   and polarization field are computed. The output is given in 
    !   following order:
    !
    !       gxx, gxy, gyy, gxz, gyz, gzz
    !
    !   the direct field gradients (for d dipoles) are in de(:,:,1)
    !   the polarization field gradients (for p dipoles) are in de(:,:,2)
    !
    ! - Wang's polynomial damping for of the coulomb potential is used for
    !   the electric field generated at polarizable site for both AMOEBA 
    !   and AMBER polarizable force-fields 
    !

    
    real(rp), dimension(6,pol_atoms,n_ipd), intent(inout) :: de
    integer(ip), intent(in) :: scr
    !
    !logical     :: Amoeba
    integer(ip) :: I, J
    !real(rp), parameter :: one = 1.0_rp
    
    !Amoeba = n_ipd.eq.2
    
    if (scr.eq.0) then
        do I = 1,pol_atoms
            do J = 1,mm_atoms
                if (polar_mm(I).eq.J) cycle

                call field_deriv_M2D(one,one,I,J,de)

            enddo
        enddo
    elseif (scr.eq.1) then
        do I = 1,pol_atoms
            do J = 1,pol_atoms
                if (I.eq.J) cycle
                
                call field_deriv_D2D(one,one,I,J,de)
            enddo
        enddo
    end if
    
    ! Remove unwanted contributions from the bonded atoms
    call multipoles_field_deriv_remove(scr,de)
    
end subroutine multipoles_field_deriv


subroutine multipoles_potentialQM(scr,v)
    use mmpol
    implicit none
    !
    ! Computes potential, electric field and field gradients at position 
    ! of MM atoms. The computed potential is added to the potential v.
    !
    ! scr  ...  0:  The sources are static multipoles of mm atoms  
    !               contained in the array q and with coordinates cmm
    !           1:  The sources are induced multipoles ipd, contained
    !               in the array ipd with coordinates cpol
    ! 
    ! - The electrostatic potential at QM sites is computed.
    !   The output is given in following order:
    !
    !       v
    !
    ! - For the AMOEBA forcefield only potential from the p dipoles is
    !   computed
    !
    
    real(rp), dimension(qm_atoms), intent(inout) :: v
    integer(ip), intent(in) :: scr
    !
    integer(ip) :: I, J
    !

    if (scr.eq.0) then
        do I = 1,qm_atoms
            do J = 1,mm_atoms

                call potential_M2Q(One,I,J,v) 

            enddo
        enddo
    elseif (scr.eq.1) then
        do I = 1,qm_atoms
            do J = 1,pol_atoms

                call potential_D2Q(Zero,One,I,J,v) 

            enddo
        enddo
    end if

end subroutine multipoles_potentialQM

subroutine multipoles_fieldQM(scr,e)
    use mmpol
    implicit none
    !
    ! Computes electric field at polarizable sites cpol. The computed
    ! field is added to the field e.
    !
    ! scr  ...  0:  The sources are static multipoles of mm atoms  
    !               contained in the array q and with coordinates cmm
    !           1:  The sources are induced multipoles ipd, contained
    !               in the array ipd with coordinates cpol.
    ! 
    ! - For AMBER force-field only electric field at QM site
    !   is computed. The output is given in following order:
    !
    !       ex, ey,ez
    !
    ! - For AMOEBA force-field electric fields only for  p dipole
    !   is computed. The output is given in following order:
    !
    !       ex, ey,ez
    !
    ! - No damping of the coulomb potential is used
    !

    
    real(rp), dimension(3,pol_atoms,n_ipd), intent(inout) :: e
    integer(ip), intent(in) :: scr
    !
    integer(ip) :: I, J
    
    !Amoeba = n_ipd.eq.2
    
    if (scr.eq.0) then
        do I = 1,qm_atoms
            do J = 1,mm_atoms
                if (polar_mm(I).eq.J) cycle

                call field_M2Q(one,I,J,e)

            enddo
        enddo
    elseif (scr.eq.1) then
        do I = 1,qm_atoms
            do J = 1,pol_atoms
                if (I.eq.J) cycle
                
                call field_D2Q(Zero,One,I,J,e)
            enddo
        enddo
    end if
    
end subroutine multipoles_fieldQM

subroutine multipoles_potential_derivQM(scr,dv)
    use mmpol
    implicit none
    !
    ! Computes derivatives of the potential, electric field and field 
    ! gradients at position of MM atoms. The computed derivatives of
    ! the potential are added to the dv. 
    !
    ! scr  ...  0:  The sources are static multipoles of mm atoms  
    !               contained in the array q and with coordinates cmm
    !           1:  The sources are induced multipoles ipd, contained
    !               in the array ipd with coordinates cpol
    ! 
    ! - For AMBER force-field only potential derivative at mm sites is 
    !   computed. The output is given in following order:
    !
    !       ex, ey, ez
    !
    ! - For AMOEBA force-field derivatives of the potential, electric 
    !   field and field gradients are computed. The output is given in
    !   following order:
    !   
    !       ex, ey, ez, gxx, gxy, gyy, gxz, gyz, gzz, hxxx, hxxy, hxxz,  
    !       hxyy, hxyz, hxzz, hyyy, hyyz, hyzz, hzzz
    !
    ! - For AMOEBA force-field, the scaling of the nearest neighbor 
    !   interactions for p and d dipoles is reversed dscale is used for 
    !   scaling potential of p dipoles and pscale is used for scaling 
    !   potential of d dipoles.
    !
    
    real(rp), dimension(ld_cder,mm_atoms), intent(inout) :: dv
    integer(ip), intent(in) :: scr
    !
    !logical     :: Amoeba
    integer(ip) :: I, J
    !
    !real(rp), parameter :: One = 1.0_rp
    
    
    !Amoeba = (ld_cder.gt.3)

    if (scr.eq.0) then
        do I = 1,mm_atoms
            do J = 1,mm_atoms
                if (I.eq.J) cycle

                call potential_deriv_M2M(One,I,J,dv) 

            enddo
        enddo
    elseif (scr.eq.1) then
        do I = 1,mm_atoms
            do J = 1,pol_atoms
                if (polar_mm(J).eq.I) cycle

                call potential_deriv_D2M(One,One,I,J,dv) 
                
            enddo
        enddo
    end if
    
    ! Remove unwanted contributions from the bonded atoms
    call multipoles_potential_deriv_remove(scr,dv)
    
end subroutine multipoles_potential_derivQM

subroutine multipoles_field_derivQM(scr,de)
    use mmpol
    implicit none
    !
    ! Computes derivatives of electric field at polarizable sites cpol.
    !  The computed field derivatives are added to the de.
    !
    ! scr  ...  0:  The sources are static multipoles of mm atoms  
    !               contained in the array q and with coordinates cmm
    !           1:  The sources are induced multipoles ipd, contained
    !               in the array ipd with coordinates cpol.
    ! 
    ! - For AMBER force-field electric field gradients at polarizable
    !   site are computed. The output is given in following order:
    !
    !       gxx, gxy, gyy, gxz, gyz, gzz
    !
    ! - For AMOEBA force-field electric field gradients for the direct
    !   and polarization field are computed. The output is given in 
    !   following order:
    !
    !       gxx, gxy, gyy, gxz, gyz, gzz
    !
    !   the direct field gradients (for d dipoles) are in de(:,:,1)
    !   the polarization field gradients (for p dipoles) are in de(:,:,2)
    !
    ! - Wang's polynomial damping for of the coulomb potential is used for
    !   the electric field generated at polarizable site for both AMOEBA 
    !   and AMBER polarizable force-fields 
    !

    
    real(rp), dimension(6,pol_atoms,n_ipd), intent(inout) :: de
    integer(ip), intent(in) :: scr
    !
    !logical     :: Amoeba
    integer(ip) :: I, J
    !real(rp), parameter :: one = 1.0_rp
    
    !Amoeba = n_ipd.eq.2
    
    if (scr.eq.0) then
        do I = 1,pol_atoms
            do J = 1,mm_atoms
                if (polar_mm(I).eq.J) cycle

                call field_deriv_M2D(one,one,I,J,de)

            enddo
        enddo
    elseif (scr.eq.1) then
        do I = 1,pol_atoms
            do J = 1,pol_atoms
                if (I.eq.J) cycle
                
                call field_deriv_D2D(one,one,I,J,de)
            enddo
        enddo
    end if
    
    ! Remove unwanted contributions from the bonded atoms
    call multipoles_field_deriv_remove(scr,de)
    
end subroutine multipoles_field_derivQM


subroutine potential_M2M(scalef,I,J,v)
    use mmpol
    implicit none
    !
    ! Computes potential, electric field and field gradients at position 
    ! of mm atom I generated by multipoles of mm atom J. The computed
    ! potential is added to the potential v
    ! 
    ! - For AMBER force-field only potential at mm sites is computed.
    !
    ! - For AMOEBA force-field potential, electric field and field gradients
    !   are computed. 
    !
    ! - no coulomb potential damping is which used for the interactions 
    !   between static sources in both MMPol and AMOEBA
    !
    ! - The total potential is scaled by factor scalef. For interaction
    !   without shielding scalef = 1
    !

    
    real(rp), dimension(ld_cart,mm_atoms), intent(inout) :: v
    !logical,intent(in)      :: Amoeba
    integer(ip), intent(in) :: I, J
    real(rp), intent(in)    :: scalef
    !
    real(rp)    :: x, y, z, dx, dy, dz, rm1, rm3, rm5, rm7, rm9, rm11
    real(rp)    :: qq, px, py, pz, qxx, qxy, qxz, qyy, qyz, qzz, DdR, QRx, QRy, QRz, QRR, fac
    !
    
    dx   = cmm(1,I) - cmm(1,J)
    dy   = cmm(2,I) - cmm(2,J)
    dz   = cmm(3,I) - cmm(3,J)
    rm1  = 0.0_rp
    rm3  = 0.0_rp
    rm5  = 0.0_rp
    rm7  = 0.0_rp
    rm9  = 0.0_rp
    rm11 = 0.0_rp

    ! If multipoles q are only charges calculate only potential
    if (.not. Amoeba) then

        call coulomb_kernel(.false.,0,dx,dy,dz,0.0_rp,0.0_rp,rm1,rm3,rm5,rm7,rm9,rm11)
        rm1 = scalef*rm1
        
        v(1,I) = v(1,I) + q( 1,J)*rm1

    ! Else the potential electric field and field gradients are calculated
    elseif (Amoeba) then
        call coulomb_kernel(.false.,4,dx,dy,dz,0.0_rp,0.0_rp,rm1,rm3,rm5,rm7,rm9,rm11)
        !
        ! Scale distances
        !
        rm1 = scalef*rm1
        rm3 = scalef*rm3
        rm5 = scalef*Three*rm5
        rm7 = scalef*f15*rm7
        rm9 = scalef*f105*rm9
        !
        ! Set multipole moments 
        !
        qq  = q( 1,J)
        px  = q( 2,J)
        py  = q( 3,J)
        pz  = q( 4,J)
        qxx = q( 5,J)
        qxy = q( 6,J)
        qyy = q( 7,J)
        qxz = q( 8,J)
        qyz = q( 9,J)
        qzz = q(10,J)

        ! intermediates
        DdR  = px*dx  + py*dy  + pz*dz
        QRx  = qxx*dx + qxy*dy + qxz*dz
        QRy  = qxy*dx + qyy*dy + qyz*dz
        QRz  = qxz*dx + qyz*dy + qzz*dz
        QRR  = QRx*dx + QRy*dy + QRz*dz

        ! Potential
        v(1,I) = v(1,I) + (qq*rm1 + DdR*rm3 + QRR*rm5)     ! q

        ! Field
        fac    = qq*rm3 + DdR*rm5 + QRR*rm7
        v(2,I) = v(2,I) + (fac*dx - px*rm3 - Two*QRx*rm5)    ! Ex
        v(3,I) = v(3,I) + (fac*dy - py*rm3 - Two*QRy*rm5)    ! Ey
        v(4,I) = v(4,I) + (fac*dz - pz*rm3 - Two*QRz*rm5)    ! Ez

        ! Field gradients
        fac  = qq*rm5 + DdR*rm7 + QRR*rm9
        v(5,I)  = v(5,I)  + (fac*dx*dx - qq*rm3 + (- DdR - Two*px*dx + Two*qxx)*rm5 - (QRR + Four*QRx*dx)*rm7)  ! Gxx
        v(6,I)  = v(6,I)  + (fac*dx*dy + (- px*dy - py*dx + Two*qxy)*rm5 - Two*(QRx*dy + QRy*dx)*rm7)           ! Gxy
        v(7,I)  = v(7,I)  + (fac*dy*dy - qq*rm3 + (- DdR - Two*py*dy + Two*qyy)*rm5 - (QRR + Four*QRy*dy)*rm7)  ! Gyy
        v(8,I)  = v(8,I)  + (fac*dx*dz + (- px*dz - pz*dx + Two*qxz)*rm5 - Two*(QRx*dz + QRz*dx)*rm7)           ! Gxz
        v(9,I)  = v(9,I)  + (fac*dy*dz + (- py*dz - pz*dy + Two*qyz)*rm5 - Two*(QRy*dz + QRz*dy)*rm7)           ! Gyz
        v(10,I) = v(10,I) + (fac*dz*dz - qq*rm3 + (- DdR - Two*pz*dz + Two*qzz)*rm5 - (QRR + Four*QRz*dz)*rm7)  ! Gzz
    end if
end subroutine potential_M2M

subroutine field_M2D(scalefd,scalefp,I,J,e)
    use mmpol
    implicit none
    !
    ! Computes electric field at position of polarizable cpol atom I
    ! generated by static multipoles of mm atom J. The computed 
    ! field is added to the field e
    ! 
    ! - For AMBER force-field only electric field at polarizable site
    !   is computed.
    !
    ! - For AMOEBA force-field electric fields for p and d dipoles
    !   are computed. 
    !
    ! - Wang's polynomial damping for of the coulomb potential is used for
    !   the electric field generated by static multipoles at polarizable 
    !   site for both AMOEBA and AMBER polarizable force-fields.
    !
    ! - The total electric field for AMBER dipoles and AMOEBA p dipoles is 
    !   scaled by factor scalefp. For the AMOEBA ddipoles the electric 
    !   field is scaled by factor scalefd. For the interaction without
    !   shielding scalefd = 1 and scalefp = 1
    !

    
    real(rp), dimension(3,pol_atoms,n_ipd), intent(inout) :: e
    real(rp), intent(in)    :: scalefd, scalefp
    !logical,intent(in)      :: Amoeba
    integer(ip),intent(in)  :: I, J
    !
    real(rp)    :: x, y, z, dx, dy, dz, rm1, rm3, rm5, rm7, rm9, rm11
    real(rp)    :: qq, px, py, pz, qxx, qxy, qxz, qyy, qyz, qzz, DdR, QRx, QRy, QRz, QRR, fac
    !
    !real(rp), parameter :: Two = 2.0_rp
            
    dx  = cpol(1,I) - cmm(1,J)
    dy  = cpol(2,I) - cmm(2,J)
    dz  = cpol(3,I) - cmm(3,J)
    rm1  = 0.0_rp
    rm3  = 0.0_rp
    rm5  = 0.0_rp
    rm7  = 0.0_rp
    rm9  = 0.0_rp
    rm11 = 0.0_rp

    ! Check if AMBER or AMOEBA polarizable force-field is used
    if (.not. Amoeba) then        ! AMBER

        call coulomb_kernel(.true.,1,dx,dy,dz,thole( polar_mm(I) ),thole(J),rm1,rm3,rm5,rm7,rm9,rm11)

        ! Field
        fac  = q( 1,J)*rm3
        e(1,I,1) = e(1,I,1) + scalefp*fac*dx
        e(2,I,1) = e(2,I,1) + scalefp*fac*dy
        e(3,I,1) = e(3,I,1) + scalefp*fac*dz

    elseif (Amoeba) then   ! AMOEBA

        call coulomb_kernel(.true.,3,dx,dy,dz,thole( polar_mm(I) ),thole(J),rm1,rm3,rm5,rm7,rm9,rm11)
        !
        ! Scale distances
        !
        rm5 = Three*rm5
        rm7 = f15*rm7
        !
        ! Set multipole moments 
        !
        qq  = q( 1,J)
        px  = q( 2,J)
        py  = q( 3,J)
        pz  = q( 4,J)
        qxx = q( 5,J)
        qxy = q( 6,J)
        qyy = q( 7,J)
        qxz = q( 8,J)
        qyz = q( 9,J)
        qzz = q(10,J)

        ! intermediates
        DdR  = px*dx  + py*dy  + pz*dz
        QRx  = qxx*dx + qxy*dy + qxz*dz
        QRy  = qxy*dx + qyy*dy + qyz*dz
        QRz  = qxz*dx + qyz*dy + qzz*dz
        QRR  = QRx*dx + QRy*dy + QRz*dz

        ! Field
        fac  = qq*rm3 + DdR*rm5 + QRR*rm7
        e(1,I,1) = e(1,I,1) + scalefd*(fac*dx - px*rm3 - Two*QRx*rm5)
        e(2,I,1) = e(2,I,1) + scalefd*(fac*dy - py*rm3 - Two*QRy*rm5)
        e(3,I,1) = e(3,I,1) + scalefd*(fac*dz - pz*rm3 - Two*QRz*rm5)
        e(1,I,2) = e(1,I,2) + scalefp*(fac*dx - px*rm3 - Two*QRx*rm5)
        e(2,I,2) = e(2,I,2) + scalefp*(fac*dy - py*rm3 - Two*QRy*rm5)
        e(3,I,2) = e(3,I,2) + scalefp*(fac*dz - pz*rm3 - Two*QRz*rm5)
    end if
end subroutine field_M2D
                
subroutine potential_deriv_M2M(scalef,I,J,dv)
    use mmpol
    implicit none
    !
    ! Computes derivatives potential, electric field and field gradients
    ! at position of mm atom I generated by multipoles of mm atom J. The 
    ! computed potential is added to the potential derivative dv
    ! 
    ! - For AMBER force-field only potential derivative at mm sites is 
    !   computed.
    !
    ! - For AMOEBA force-field derivatives of the potential, electric 
    !   field and field gradients are computed. 
    !
    ! - no coulomb potential damping is used for the interactions 
    !   between static sources in both MMPol and AMOEBA
    !
    ! - The total potential is scaled by factor scalef. For interaction
    !   without shielding scalef = 1
    !

    
    real(rp), dimension(ld_cder,mm_atoms), intent(inout) :: dv
    !logical,intent(in)      :: Amoeba
    integer(ip), intent(in) :: I, J
    real(rp), intent(in)    :: scalef
    !
    real(rp)    :: x, y, z, dx, dy, dz, rm1, rm3, rm5, rm7, rm9, rm11
    real(rp)    :: qq, px, py, pz, qxx, qxy, qxz, qyy, qyz, qzz, DdR, QRx, QRy, QRz, QRR, fac
    !
    real(rp), parameter :: f945 = 945.0_rp
    
    dx   = cmm(1,I) - cmm(1,J)
    dy   = cmm(2,I) - cmm(2,J)
    dz   = cmm(3,I) - cmm(3,J)
    rm1  = 0.0_rp
    rm3  = 0.0_rp
    rm5  = 0.0_rp
    rm7  = 0.0_rp
    rm9  = 0.0_rp
    rm11 = 0.0_rp

    ! If multipoles q are only charges calculate only potential
    if (.not. Amoeba) then
        call coulomb_kernel(.false.,1,dx,dy,dz,0.0_rp,0.0_rp,rm1,rm3,rm5,rm7,rm9,rm11)

        ! Scaling of the interactions
        rm3 = scalef*rm3
        
        ! Field
        fac = q( 1,J)*rm3
        dv(1,I) = dv(1,I) + fac*dx  ! ex
        dv(2,I) = dv(2,I) + fac*dy  ! ey
        dv(3,I) = dv(3,I) + fac*dz  ! ez
            
    ! Else the potential electric field and field gradients are calculated
    elseif (Amoeba) then
        call coulomb_kernel(.false.,5,dx,dy,dz,0.0_rp,0.0_rp,rm1,rm3,rm5,rm7,rm9,rm11)
 
        ! Scaling of the interactions
        rm1  = scalef*rm1
        rm3  = scalef*rm3
        rm5  = scalef*Three*rm5
        rm7  = scalef*f15*rm7
        rm9  = scalef*f105*rm9
        rm11 = scalef*f945*rm11
        
        ! Multipoles
        qq  = q( 1,J)
        px  = q( 2,J)
        py  = q( 3,J)
        pz  = q( 4,J)
        qxx = q( 5,J)
        qxy = q( 6,J)
        qyy = q( 7,J)
        qxz = q( 8,J)
        qyz = q( 9,J)
        qzz = q(10,J)
        
        ! Intermediates
        DdR  = px*dx  + py*dy  + pz*dz
        QRx  = qxx*dx + qxy*dy + qxz*dz
        QRy  = qxy*dx + qyy*dy + qyz*dz
        QRz  = qxz*dx + qyz*dy + qzz*dz
        QRR  = QRx*dx + QRy*dy + QRz*dz

        ! Field
        fac     = qq*rm3 + DdR*rm5 + QRR*rm7
        dv(1,I) = dv(1,I) + fac*dx - px*rm3 - Two*QRx*rm5  ! ex
        dv(2,I) = dv(2,I) + fac*dy - py*rm3 - Two*QRy*rm5  ! ey
        dv(3,I) = dv(3,I) + fac*dz - pz*rm3 - Two*QRz*rm5  ! ez

        ! Field gradients
        fac     = qq*rm5 + DdR*rm7 + QRR*rm9
        dv(4,I) = dv(4,I) + fac*dx*dx - qq*rm3 + (- DdR - Two*px*dx + Two*qxx)*rm5 - (QRR + Four*QRx*dx)*rm7  ! gxx
        dv(5,I) = dv(5,I) + fac*dx*dy + (- px*dy - py*dx + Two*qxy)*rm5 - Two*(QRx*dy + QRy*dx)*rm7           ! gxy
        dv(6,I) = dv(6,I) + fac*dy*dy - qq*rm3 + (- DdR - Two*py*dy + Two*qyy)*rm5 - (QRR + Four*QRy*dy)*rm7  ! gyy
        dv(7,I) = dv(7,I) + fac*dx*dz + (- px*dz - pz*dx + Two*qxz)*rm5 - Two*(QRx*dz + QRz*dx)*rm7           ! gxz
        dv(8,I) = dv(8,I) + fac*dy*dz + (- py*dz - pz*dy + Two*qyz)*rm5 - Two*(QRy*dz + QRz*dy)*rm7           ! gyz
        dv(9,I) = dv(9,I) + fac*dz*dz - qq*rm3 + (- DdR - Two*pz*dz + Two*qzz)*rm5 - (QRR + Four*QRz*dz)*rm7  ! gzz
        
        ! Field hessians
        fac = qq*rm7 + DdR*rm9 + QRR*rm11
        dv(10,I) = dv(10,I) + fac*dx*dx*dx + rm5*Three*(px - qq*dx) + rm7*Three*(Two*dx*Qxx + Two*QRx - px*dx*dx - DdR*dx) &
                   - rm9*Three*(Two*QRx*dx*dx + QRR*dx)                                                 ! hxxx
        dv(11,I) = dv(11,I) + fac*dx*dx*dy + rm5*(py - qq*dy) + rm7*(Four*dx*Qxy + Two*dy*Qxx + Two*QRy    &
                   - Two*dx*dy*px - py*dx*dx - DdR*dy) - rm9*(Four*QRx*dx*dy + Two*QRy*dx*dx + QRR*dy)  ! hxxy
        dv(12,I) = dv(12,I) + fac*dx*dx*dz + rm5*(pz - qq*dz) + rm7*(Four*dx*Qxz + Two*dz*Qxx + Two*QRz    &
                   - Two*dx*dz*px - pz*dx*dx - DdR*dz) - rm9*(Four*QRx*dx*dz + Two*QRz*dx*dx + QRR*dz)  ! hxxz
        dv(13,I) = dv(13,I) + fac*dx*dy*dy + rm5*(px - qq*dx) + rm7*(Four*dy*Qxy + Two*dx*Qyy + Two*QRx    &
                   - Two*dx*dy*py - px*dy*dy - DdR*dx) - rm9*(Four*QRy*dx*dy + Two*QRx*dy*dy + QRR*dx)  ! hxyy
        dv(14,I) = dv(14,I) + fac*dx*dy*dz + rm7*(Two*dx*Qyz + Two*dy*Qxz + Two*dz*Qxy - dx*dy*pz - dy*dz*px - dz*dx*py)  &
                   - rm9*Two*(QRx*dy*dz + QRy*dx*dz + QRz*dx*dy)                                        ! hxyz
        dv(15,I) = dv(15,I) + fac*dx*dz*dz + rm5*(px - qq*dx) + rm7*(Four*dz*Qxz + Two*dx*Qzz + Two*QRx    &
                   - Two*dx*dz*pz - px*dz*dz - DdR*dx)- rm9*(Four*QRz*dx*dz + Two*QRx*dz*dz + QRR*dx)   ! hxzz
        dv(16,I) = dv(16,I) + fac*dy*dy*dy + rm5*Three*(py - qq*dy) + rm7*Three*(Two*dy*Qyy + Two*QRy - py*dy*dy - DdR*dy) &
                   - rm9*Three*(Two*QRy*dy*dy + QRR*dy)                                                 ! hyyy
        dv(17,I) = dv(17,I) + fac*dy*dy*dz + rm5*(pz - qq*dz) + rm7*(Four*dy*Qyz + Two*dz*Qyy + Two*QRz    &
                   - Two*dy*dz*py - pz*dy*dy - DdR*dz) - rm9*(Four*QRy*dy*dz + Two*QRz*dy*dy + QRR*dz)  ! hyyz
        dv(18,I) = dv(18,I) + fac*dy*dz*dz + rm5*(py - qq*dy) + rm7*(Four*dz*Qyz + Two*dy*Qzz + Two*QRy    &
                   - Two*dy*dz*pz - py*dz*dz - DdR*dy) - rm9*(Four*QRz*dy*dz + Two*QRy*dz*dz + QRR*dy)  ! hyzz
        dv(19,I) = dv(19,I) + fac*dz*dz*dz + rm5*Three*(pz - qq*dz) + rm7*Three*(Two*dz*Qzz + Two*QRz - pz*dz*dz - DdR*dz) &
                   - rm9*Three*(Two*QRz*dz*dz + QRR*dz)                                                 ! hzzz  
    end if
end subroutine potential_deriv_M2M


subroutine field_deriv_M2D(scalefd,scalefp,I,J,de)
    use mmpol
    implicit none
    !
    ! Computes electric field derivatives at position of polarizable cpol
    ! atom I generated by static multipoles of mm atom J. The computed 
    ! field is added to the field de
    ! 
    ! - For AMBER force-field only electric field derivatives at polarizable
    !   site is computed.
    !
    ! - For AMOEBA force-field electric field derivatives at p and d dipoles
    !   is computed.
    !
    ! - Wang's polynomial damping for of the coulomb potential is used for
    !   the electric field generated by the static multipoles at polarizable 
    !   site for both AMOEBA and AMBER polarizable force-fields 
    !
    ! - The total electric field for AMBER dipoles and AMOEBA p dipoles are 
    !   scaled by factor scalefp. For the AMOEBA d dipoles the electric 
    !   field are scaled by the factor scalefd. For the interaction without
    !   shielding scalefd = 1 and scalefp = 1
    
    real(rp), dimension(6,pol_atoms,n_ipd), intent(inout) :: de
    real(rp), intent(in)    :: scalefd, scalefp
    !logical,intent(in)      :: Amoeba
    integer(ip),intent(in)  :: I, J
    !
    real(rp)    :: x, y, z, dx, dy, dz, rm1, rm3, rm5, rm7, rm9, rm11
    real(rp)    :: qq, px, py, pz, qxx, qxy, qxz, qyy, qyz, qzz, DdR, QRx, QRy, QRz, QRR, fac
    real(rp)    :: fac1, fac2, fGxx, fGxy, fGxz, fGyy, fGyz, fGzz
    !
    !real(rp), parameter :: Two = 2.0_rp, Three = 3.0_rp, Four = 4.0_rp
    
    dx   = cpol(1,I) - cmm(1,J)
    dy   = cpol(2,I) - cmm(2,J)
    dz   = cpol(3,I) - cmm(3,J)
    rm1  = 0.0_rp
    rm3  = 0.0_rp
    rm5  = 0.0_rp
    rm7  = 0.0_rp
    rm9  = 0.0_rp
    rm11 = 0.0_rp
    
    ! Check if AMBER or AMOEBA polarizable force-field is used
    if (.not. Amoeba) then        ! AMBER

        call coulomb_kernel(.true.,2,dx,dy,dz,thole( polar_mm(I) ),thole(J),rm1,rm3,rm5,rm7,rm9,rm11)

        ! Field derivatives
        fac1 = scalefp*q(1,J)*rm3
        fac2 = scalefp*Three*q(1,J)*rm5
        de(1,I,1) = de(1,I,1) - fac1 + fac2*dx*dx  ! gxx
        de(2,I,1) = de(2,I,1)  + fac2*dx*dy        ! gxy
        de(3,I,1) = de(3,I,1) - fac1 + fac2*dy*dy  ! gyy
        de(4,I,1) = de(4,I,1) + fac2*dx*dz         ! gxz
        de(5,I,1) = de(5,I,1) + fac2*dy*dz         ! gyz
        de(6,I,1) = de(6,I,1) - fac1 + fac2*dz*dz  ! gzz
        
    elseif (Amoeba) then   ! AMOEBA

        call coulomb_kernel(.true.,4,dx,dy,dz,thole( polar_mm(I) ),thole(J),rm1,rm3,rm5,rm7,rm9,rm11)
        
        ! Scaling of the distances
        rm5  = Three*rm5
        rm7  = f15*rm7
        rm9  = f105*rm9
        
        ! Multipoles
        qq  = q( 1,J)
        px  = q( 2,J)
        py  = q( 3,J)
        pz  = q( 4,J)
        qxx = q( 5,J)
        qxy = q( 6,J)
        qyy = q( 7,J)
        qxz = q( 8,J)
        qyz = q( 9,J)
        qzz = q(10,J)

        ! Intermediates
        DdR  = px*dx  + py*dy  + pz*dz
        QRx  = qxx*dx + qxy*dy + qxz*dz
        QRy  = qxy*dx + qyy*dy + qyz*dz
        QRz  = qxz*dx + qyz*dy + qzz*dz
        QRR  = QRx*dx + QRy*dy + QRz*dz
        
        ! Field deriatives
        fac  = qq*rm5 + DdR*rm7 + QRR*rm9
        fGxx = fac*dx*dx - qq*rm3 + (- DdR - Two*px*dx + Two*qxx)*rm5 - (QRR + Four*QRx*dx)*rm7
        fGxy = fac*dx*dy + (- px*dy - py*dx + Two*qxy)*rm5 - Two*(QRx*dy + QRy*dx)*rm7
        fGyy = fac*dy*dy - qq*rm3 + (- DdR - Two*py*dy + Two*qyy)*rm5 - (QRR + Four*QRy*dy)*rm7
        fGxz = fac*dx*dz + (- px*dz - pz*dx + Two*qxz)*rm5 - Two*(QRx*dz + QRz*dx)*rm7
        fGyz = fac*dy*dz + (- py*dz - pz*dy + Two*qyz)*rm5 - Two*(QRy*dz + QRz*dy)*rm7
        fGzz = fac*dz*dz - qq*rm3 + (- DdR - Two*pz*dz + Two*qzz)*rm5 - (QRR + Four*QRz*dz)*rm7
        
        ! d dipoles
        de(1,I,1) = de(1,I,1) + scalefd*fGxx  ! gxx
        de(2,I,1) = de(2,I,1) + scalefd*fGxy  ! gxy
        de(3,I,1) = de(3,I,1) + scalefd*fGyy  ! gyy
        de(4,I,1) = de(4,I,1) + scalefd*fGxz  ! gxz
        de(5,I,1) = de(5,I,1) + scalefd*fGyz  ! gyz
        de(6,I,1) = de(6,I,1) + scalefd*fGzz  ! gzz
        
        ! p dipoles
        de(1,I,2) = de(1,I,2) + scalefp*fGxx  ! gxx
        de(2,I,2) = de(2,I,2) + scalefp*fGxy  ! gxy
        de(3,I,2) = de(3,I,2) + scalefp*fGyy  ! gyy
        de(4,I,2) = de(4,I,2) + scalefp*fGxz  ! gxz
        de(5,I,2) = de(5,I,2) + scalefp*fGyz  ! gyz
        de(6,I,2) = de(6,I,2) + scalefp*fGzz  ! gzz        
    end if

end subroutine field_deriv_M2D


subroutine potential_D2M(scalefd,scalefp,I,J,v)
    use mmpol
    implicit none
    !
    ! Computes potential, electric field and field gradients at position 
    ! of mm atom I generated by induced dipoles of cpol atom J. The computed
    ! potential is added to the potential v.
    ! 
    ! - For AMBER force-field only potential at mm sites is computed.
    !
    ! - For AMOEBA force-field potential, electric field and field gradients
    !   are computed. 
    !
    ! - Wang's polynomial damping for of the coulomb potential is used for
    !   the electric field generated by the induced dipoles at static multipoles 
    !   for both AMOEBA and AMBER polarizable force-fields 
    !
    ! - The total potential from AMBER dipoles and AMOEBA p dipoles is 
    !   scaled by factor scalefp. For the AMOEBA d dipoles the potential 
    !   is scaled by the factor scalefd. For the interaction without
    !   shielding scalefd = 1 and scalefp = 1
    !

    
    real(rp), dimension(ld_cart,mm_atoms), intent(inout) :: v
    !logical,intent(in)      :: Amoeba
    integer(ip), intent(in) :: I, J
    real(rp), intent(in)    :: scalefd,scalefp
    !
    real(rp)    :: x, y, z, dx, dy, dz, rm1, rm3, rm5, rm7, rm9, rm11
    real(rp)    :: dpx, dpy, dpz, ppx, ppy, ppz, px, py, pz, DdR, DpdR, DddR
    real(rp)    :: facd, facp, fac
    !
    !real(rp), parameter :: Two = 2.0_rp, Four = 4.0_rp
    
    dx   = cmm(1,I) - cpol(1,J)
    dy   = cmm(2,I) - cpol(2,J)
    dz   = cmm(3,I) - cpol(3,J)
    rm1  = 0.0_rp
    rm3  = 0.0_rp
    rm5  = 0.0_rp
    rm7  = 0.0_rp
    rm9  = 0.0_rp
    rm11 = 0.0_rp

    ! If multipoles q are only charges calculate only potential
    if (.not. Amoeba) then

        call coulomb_kernel(.true.,1,dx,dy,dz,thole(I),thole( polar_mm(J) ),rm1,rm3,rm5,rm7,rm9,rm11)
        
        ! Scaling of the interaction
        rm3 = scalefp*rm3
        
        ! Intermediates
        px  = ipd( 1, J, 1)
        py  = ipd( 2, J, 1)
        pz  = ipd( 3, J, 1)
        DdR  = px*dx  + py*dy  + pz*dz
        
        ! Potential of the induced dipoles
        v(1,I) = v(1,I) + DdR*rm3
        
    ! Else the potential electric field and field gradients are calculated
    elseif (Amoeba) then
        call coulomb_kernel(.true.,3,dx,dy,dz,thole(I),thole( polar_mm(J) ),rm1,rm3,rm5,rm7,rm9,rm11)
        
        rm1 = rm1
        rm3 = rm3
        rm5 = Three*rm5
        rm7 = f15*rm7
        rm9 = f105*rm9
        
        
        ! Intermediates p and d dipoles
        dpx  = scalefd*ipd( 1, J, 1)
        dpy  = scalefd*ipd( 2, J, 1)
        dpz  = scalefd*ipd( 3, J, 1)
        ppx  = scalefp*ipd( 1, J, 2)
        ppy  = scalefp*ipd( 2, J, 2)
        ppz  = scalefp*ipd( 3, J, 2)
        DddR = dpx*dx  + dpy*dy  + dpz*dz   ! DdR for d induced dipoles 
        DpdR = ppx*dx  + ppy*dy  + ppz*dz   ! DdR for p induced dipoles 

        ! Potential
        v(1,I) = v(1,I) + DddR*rm3  ! Vd
        v(1,I) = v(1,I) + DpdR*rm3  ! Vp

        
        ! Field
        facd   = DddR*rm5
        facp   = DpdR*rm5
        v(2,I) = v(2,I) + (facd*dx - dpx*rm3)   ! Ex from d dipoles
        v(3,I) = v(3,I) + (facd*dy - dpy*rm3)   ! Ey from d dipoles
        v(4,I) = v(4,I) + (facd*dz - dpz*rm3)   ! Ez from d dipoles
        
        v(2,I) = v(2,I) + (facp*dx - ppx*rm3)   ! Ex from p dipoles
        v(3,I) = v(3,I) + (facp*dy - ppy*rm3)   ! Ey from p dipoles
        v(4,I) = v(4,I) + (facp*dz - ppz*rm3)   ! Ez from p dipoles

       
        ! Field gradients
        facd  = DddR*rm7
        facp  = DpdR*rm7
        v(5,I)  = v(5,I)  + (facd*dx*dx + (- DddR - Two*dpx*dx)*rm5)  ! Gxx from d dipoles
        v(6,I)  = v(6,I)  + (facd*dx*dy + (- dpx*dy - dpy*dx)*rm5)    ! Gxy from d dipoles
        v(7,I)  = v(7,I)  + (facd*dy*dy + (- DddR - Two*dpy*dy)*rm5)  ! Gyy from d dipoles
        v(8,I)  = v(8,I)  + (facd*dx*dz + (- dpx*dz - dpz*dx)*rm5)    ! Gxz from d dipoles
        v(9,I)  = v(9,I)  + (facd*dy*dz + (- dpy*dz - dpz*dy)*rm5)    ! Gyz from d dipoles
        v(10,I) = v(10,I) + (facd*dz*dz + (- DddR - Two*dpz*dz)*rm5)  ! Gzz from d dipoles
        
        !write(*,*) I,J,v(8,I),rm5,dpz,dpx,dz,dx,facd
        
        v(5,I)  = v(5,I)  + (facp*dx*dx + (- DpdR - Two*ppx*dx)*rm5)  ! Gxx from p dipoles
        v(6,I)  = v(6,I)  + (facp*dx*dy + (- ppx*dy - ppy*dx)*rm5)    ! Gxy from p dipoles
        v(7,I)  = v(7,I)  + (facp*dy*dy + (- DpdR - Two*ppy*dy)*rm5)  ! Gyy from p dipoles
        v(8,I)  = v(8,I)  + (facp*dx*dz + (- ppx*dz - ppz*dx)*rm5)    ! Gxz from p dipoles
        v(9,I)  = v(9,I)  + (facp*dy*dz + (- ppy*dz - ppz*dy)*rm5)    ! Gyz from p dipoles
        v(10,I) = v(10,I) + (facp*dz*dz + (- DpdR - Two*ppz*dz)*rm5)  ! Gzz from p dipoles
        
        !write(*,*) I,J,v(8,I),rm5,ppz,ppx,dz,dx,facp
        
!        ! Alternative formulation (more compact and efficient)
!        
!        ! Intermediates
!        px  = (scalefd*ipd( 1, J, 1) + scalefp*ipd( 1, J, 2))
!        py  = (scalefd*ipd( 2, J, 1) + scalefp*ipd( 2, J, 2))
!        pz  = (scalefd*ipd( 3, J, 1) + scalefp*ipd( 3, J, 2))
!        DdR = px*dx  + py*dy  + pz*dz
!
!        ! Potential
!        v(1,I) = v(1,I) + DdR*rm3  ! V
!        
!        ! Field
!        fac   = DdR*rm5
!        v(2,I) = v(2,I) + (fac*dx - px*rm3)   ! Ex
!        v(3,I) = v(3,I) + (fac*dy - py*rm3)   ! Ey
!        v(4,I) = v(4,I) + (fac*dz - pz*rm3)   ! Ez
!        
!        ! Field gradients
!        fac  = DdR*rm7
!        v(5,I)  = v(5,I)  + (fac*dx*dx + (- DdR - Two*px*dx)*rm5)  ! Gxx
!        v(6,I)  = v(6,I)  + (fac*dx*dy + (- px*dy - py*dx)*rm5)    ! Gxy
!        v(7,I)  = v(7,I)  + (fac*dx*dz + (- px*dz - pz*dx)*rm5)    ! Gxz
!        v(8,I)  = v(8,I)  + (fac*dy*dy + (- DdR - Two*py*dy)*rm5)  ! Gyy
!        v(9,I)  = v(9,I)  + (fac*dy*dz + (- py*dz - pz*dy)*rm5)    ! Gyz
!        v(10,I) = v(10,I) + (fac*dz*dz + (- DdR - Two*pz*dz)*rm5)  ! Gzz
        
    end if
end subroutine potential_D2M


subroutine field_D2D(scalefd,scalefp,I,J,e)
    use mmpol
    implicit none
    !
    ! Computes electric field at position of polarizable cpol atom I
    ! generated by polarizable cpol atom J. The computed field is added
    ! to the field e
    ! 
    ! - For AMBER force-field only electric field at polarizable site
    !   is computed.
    !
    ! - For AMOEBA only p dipoles contribute to the electric field at 
    ! ` polarizable sites
    !
    ! - Wang's polynomial damping for of the coulomb potential is used for
    !   the electric field generated by induced dipoles at polarizable 
    !   site for both AMOEBA and AMBER polarizable force-fields 
    !
    ! - The total electric field for AMBER dipoles and AMOEBA p dipoles is 
    !   scaled by factor scalefp. For the AMOEBA d dipoles the electric 
    !   field is scaled by factor scalefd. For the interaction without
    !   shielding scalefd = 1 and scalefp = 1
    !

    real(rp), dimension(3,pol_atoms,n_ipd), intent(inout) :: e
    real(rp), intent(in)    :: scalefd, scalefp
    !logical,intent(in)      :: Amoeba
    integer(ip),intent(in)  :: I, J
    !
    real(rp)    :: x, y, z, dx, dy, dz, rm1, rm3, rm5, rm7, rm9, rm11
    real(rp)    :: px, py, pz, DdR
    
    dx   = cpol(1,I) - cpol(1,J)
    dy   = cpol(2,I) - cpol(2,J)
    dz   = cpol(3,I) - cpol(3,J)
    rm1  = 0.0_rp
    rm3  = 0.0_rp
    rm5  = 0.0_rp
    rm7  = 0.0_rp
    rm9  = 0.0_rp
    rm11 = 0.0_rp

    ! If multipoles q are only charges calculate only potential
    if (.not. Amoeba) then

        call coulomb_kernel(.true.,2,dx,dy,dz,thole( polar_mm(I) ),thole( polar_mm(J) ),rm1,rm3,rm5,rm7,rm9,rm11)
        rm1 = rm1
        rm3 = rm3
        rm5 = Three*rm5
        
        ! Intermediates
        px  = ipd( 1, J, 1)
        py  = ipd( 2, J, 1)
        pz  = ipd( 3, J, 1)
        DdR  = px*dx  + py*dy  + pz*dz
        
        ! Field from the induced dipoles
        e(1,I,1) = e(1,I,1) + scalefp*(DdR*rm5*dx - px*rm3)
        e(2,I,1) = e(2,I,1) + scalefp*(DdR*rm5*dy - py*rm3)
        e(3,I,1) = e(3,I,1) + scalefp*(DdR*rm5*dz - pz*rm3)
        
    elseif (Amoeba) then
        call coulomb_kernel(.true.,2,dx,dy,dz,thole( polar_mm(I) ),thole( polar_mm(J) ),rm1,rm3,rm5,rm7,rm9,rm11)
        rm1 = rm1
        rm3 = rm3
        rm5 = Three*rm5
        
        ! Intermediates p dipoles
        px  = ipd( 1, J, 2)
        py  = ipd( 2, J, 2)
        pz  = ipd( 3, J, 2)
        
        ! Field from the induced p dipoles
        DdR = px*dx + py*dy + pz*dz
        e(1,I,2) = e(1,I,2) + scalefp*(DdR*rm5*dx - px*rm3)
        e(2,I,2) = e(2,I,2) + scalefp*(DdR*rm5*dy - py*rm3)
        e(3,I,2) = e(3,I,2) + scalefp*(DdR*rm5*dz - pz*rm3)
        
        ! Intermediates d dipoles
        px  = ipd( 1, J, 1)
        py  = ipd( 2, J, 1)
        pz  = ipd( 3, J, 1)
        
        ! Field from the induced d dipoles
        DdR = px*dx + py*dy + pz*dz
        e(1,I,1) = e(1,I,1) + scalefd*(DdR*rm5*dx - px*rm3)
        e(2,I,1) = e(2,I,1) + scalefd*(DdR*rm5*dy - py*rm3)
        e(3,I,1) = e(3,I,1) + scalefd*(DdR*rm5*dz - pz*rm3)
        
!          DdRD = DD(1,J)*dx + DD(2,J)*dy + DD(3,J)*dz
!          rm3  = fl3/r3
!          rm5  = Three*fl5/r5
!          fEDx = fEDx + Scr(J)*(DdRD*rm5*dx - DD(1,J)*rm3)
!          fEDy = fEDy + Scr(J)*(DdRD*rm5*dy - DD(2,J)*rm3)
!          fEDz = fEDz + Scr(J)*(DdRD*rm5*dz - DD(3,J)*rm3)
        
    end if
end subroutine field_D2D


subroutine potential_deriv_D2M(scalefd,scalefp,I,J,dv)
    use mmpol
    implicit none
    !
    ! Computes derivatives potential, electric field and field gradients at position 
    ! of mm atom I generated by induced dipoles of cpol atom J. The computed
    ! potential is added to the potential derivative dv
    ! 
    ! - For AMBER force-field only potential at mm sites is computed.
    !
    ! - For AMOEBA force-field potential, electric field and field gradients
    !   are computed. 
    !
    ! - The total potential from AMBER dipoles and AMOEBA p dipoles is 
    !   scaled by factor scalefp. For the AMOEBA d dipoles the potential 
    !   is scaled by the factor scalefd. For the interaction without
    !   shielding scalefd = 1 and scalefp = 1
    !

    
    real(rp), dimension(ld_cder,mm_atoms), intent(inout) :: dv
    !logical,intent(in)      :: Amoeba
    integer(ip), intent(in) :: I, J
    real(rp), intent(in)    :: scalefd, scalefp
    !
    real(rp)    :: x, y, z, dx, dy, dz, rm1, rm3, rm5, rm7, rm9, rm11
    real(rp)    :: px, py,pz, DdR, dpx, dpy, dpz, ppx, ppy, ppz, DpdR, DddR, facd, facp !,fac,DdR
    !
    real(rp), parameter :: f945 = 945.0_rp
    
    dx   = cmm(1,I) - cpol(1,J)
    dy   = cmm(2,I) - cpol(2,J)
    dz   = cmm(3,I) - cpol(3,J)
    rm1  = 0.0_rp
    rm3  = 0.0_rp
    rm5  = 0.0_rp
    rm7  = 0.0_rp
    rm9  = 0.0_rp
    rm11 = 0.0_rp

    ! If multipoles q are only charges calculate only potential
    if (.not. Amoeba) then
        call coulomb_kernel(.true.,2,dx,dy,dz,thole(I),thole( polar_mm(J) ),rm1,rm3,rm5,rm7,rm9,rm11)
        
        rm1 = rm1
        rm3 = rm3
        rm5 = Three*rm5
        
        ! Intermediates
        px  = ipd( 1, J, 1)
        py  = ipd( 2, J, 1)
        pz  = ipd( 3, J, 1)
        DdR = px*dx + py*dy + pz*dz

        ! Field
        dv(1,I) = dv(1,I) + scalefp*(DdR*rm5*dx - px*rm3)  ! ex
        dv(2,I) = dv(2,I) + scalefp*(DdR*rm5*dy - py*rm3)  ! ey
        dv(3,I) = dv(3,I) + scalefp*(DdR*rm5*dz - pz*rm3)  ! ez
            
    ! Else the potential electric field and field gradients are calculated
    elseif (Amoeba) then
        call coulomb_kernel(.true.,4,dx,dy,dz,thole(I),thole( polar_mm(J) ),rm1,rm3,rm5,rm7,rm9,rm11)
        rm1 = rm1
        rm3 = rm3
        rm5 = Three*rm5
        rm7 = f15*rm7
        rm9 = f105*rm9
        !rm11 = scalef*f945*rm11
        !write(*,*) I,J,rm1, rm3, rm5, rm7, rm9
        
        
        ! Intermediates p and d dipoles
        dpx  = scalefd*ipd( 1, J, 1)
        dpy  = scalefd*ipd( 2, J, 1)
        dpz  = scalefd*ipd( 3, J, 1)
        ppx  = scalefp*ipd( 1, J, 2)
        ppy  = scalefp*ipd( 2, J, 2)
        ppz  = scalefp*ipd( 3, J, 2)
        
        ! All the quantities are linear in the dipole, therefore we can just sum the p and d dipoles and then use these
        ! dipoles for calculation of the potential derivatives (we dont need to do it separately for the p and d dipoles)
        
        DddR = dpx*dx + dpy*dy + dpz*dz
        DpdR = ppx*dx + ppy*dy + ppz*dz
        facd = DddR*rm5
        facp = DpdR*rm5
        
        ! Field
        dv(1,I) = dv(1,I) + pt5*(facp*dx - ppx*rm3)   ! ex from p dipoles
        dv(2,I) = dv(2,I) + pt5*(facp*dy - ppy*rm3)   ! ey from p dipoles
        dv(3,I) = dv(3,I) + pt5*(facp*dz - ppz*rm3)   ! ez from p dipoles
        dv(1,I) = dv(1,I) + pt5*(facd*dx - dpx*rm3)   ! ex from d dipoles
        dv(2,I) = dv(2,I) + pt5*(facd*dy - dpy*rm3)   ! ey from d dipoles
        dv(3,I) = dv(3,I) + pt5*(facd*dz - dpz*rm3)   ! ez from d dipoles
        
        ! Field gradients
        facd = DddR*rm7
        facp = DpdR*rm7
        dv(4,I) = dv(4,I) + pt5*(facp*dx*dx-(DpdR+Two*ppx*dx)*rm5)    ! gxx from p dipoles
        dv(5,I) = dv(5,I) + pt5*(facp*dx*dy-(ppx*dy+ppy*dx)*rm5)      ! gxy from p dipoles
        dv(6,I) = dv(6,I) + pt5*(facp*dy*dy-(DpdR+Two*ppy*dy)*rm5)    ! gyy from p dipoles
        dv(7,I) = dv(7,I) + pt5*(facp*dx*dz-(ppx*dz+ppz*dx)*rm5)      ! gxz from p dipoles
        dv(8,I) = dv(8,I) + pt5*(facp*dy*dz-(ppy*dz+ppz*dy)*rm5)      ! gyz from p dipoles
        dv(9,I) = dv(9,I) + pt5*(facp*dz*dz-(DpdR+two*ppz*dz)*rm5)    ! gzz from p dipoles

        dv(4,I) = dv(4,I) + pt5*(facd*dx*dx-(DddR+Two*dpx*dx)*rm5)    ! gxx from d dipoles
        dv(5,I) = dv(5,I) + pt5*(facd*dx*dy-(dpx*dy+dpy*dx)*rm5)      ! gxy from d dipoles
        dv(6,I) = dv(6,I) + pt5*(facd*dy*dy-(DddR+Two*dpy*dy)*rm5)    ! gyy from d dipoles
        dv(7,I) = dv(7,I) + pt5*(facd*dx*dz-(dpx*dz+dpz*dx)*rm5)      ! gxz from d dipoles
        dv(8,I) = dv(8,I) + pt5*(facd*dy*dz-(dpy*dz+dpz*dy)*rm5)      ! gyz from d dipoles
        dv(9,I) = dv(9,I) + pt5*(facd*dz*dz-(DddR+two*dpz*dz)*rm5)    ! gzz from d dipoles

        
        ! Field hessians
        facd = DddR*rm9
        facp = DpdR*rm9
        dv(10,I) = dv(10,I) + pt5*(facp*dx*dx*dx + rm5*Three*ppx - rm7*Three*(ppx*dx*dx + DpdR*dx))       ! hxxx from p dipoles
        dv(11,I) = dv(11,I) + pt5*(facp*dx*dx*dy + rm5*ppy - rm7*(Two*dx*dy*ppx + ppy*dx*dx + DpdR*dy))   ! hxxy from p dipoles
        dv(12,I) = dv(12,I) + pt5*(facp*dx*dx*dz + rm5*ppz - rm7*(Two*dx*dz*ppx + ppz*dx*dx + DpdR*dz))   ! hxxz from p dipoles
        dv(13,I) = dv(13,I) + pt5*(facp*dx*dy*dy + rm5*ppx - rm7*(Two*dx*dy*ppy + ppx*dy*dy + DpdR*dx))   ! hxyy from p dipoles
        dv(14,I) = dv(14,I) + pt5*(facp*dx*dy*dz - rm7*(dx*dy*ppz + dy*dz*ppx + dz*dx*ppy))               ! hxyz from p dipoles
        dv(15,I) = dv(15,I) + pt5*(facp*dx*dz*dz + rm5*ppx - rm7*(Two*dx*dz*ppz + ppx*dz*dz + DpdR*dx))   ! hxzz from p dipoles
        dv(16,I) = dv(16,I) + pt5*(facp*dy*dy*dy + rm5*Three*ppy - rm7*Three*(ppy*dy*dy + DpdR*dy))       ! hyyy from p dipoles
        dv(17,I) = dv(17,I) + pt5*(facp*dy*dy*dz + rm5*ppz - rm7*(Two*dy*dz*ppy + ppz*dy*dy + DpdR*dz))   ! hyyz from p dipoles
        dv(18,I) = dv(18,I) + pt5*(facp*dy*dz*dz + rm5*ppy - rm7*(Two*dy*dz*ppz + ppy*dz*dz + DpdR*dy))   ! hyzz from p dipoles
        dv(19,I) = dv(19,I) + pt5*(facp*dz*dz*dz + rm5*Three*ppz - rm7*Three*(ppz*dz*dz + DpdR*dz))       ! hzzz from p dipoles
        
        dv(10,I) = dv(10,I) + pt5*(facd*dx*dx*dx + rm5*Three*dpx - rm7*Three*(dpx*dx*dx + DddR*dx))       ! hxxx from d dipoles
        dv(11,I) = dv(11,I) + pt5*(facd*dx*dx*dy + rm5*dpy - rm7*(Two*dx*dy*dpx + dpy*dx*dx + DddR*dy))   ! hxxy from d dipoles
        dv(12,I) = dv(12,I) + pt5*(facd*dx*dx*dz + rm5*dpz - rm7*(Two*dx*dz*dpx + dpz*dx*dx + DddR*dz))   ! hxxz from d dipoles
        dv(13,I) = dv(13,I) + pt5*(facd*dx*dy*dy + rm5*dpx - rm7*(Two*dx*dy*dpy + dpx*dy*dy + DddR*dx))   ! hxyy from d dipoles
        dv(14,I) = dv(14,I) + pt5*(facd*dx*dy*dz - rm7*(dx*dy*dpz + dy*dz*dpx + dz*dx*dpy))               ! hxyz from d dipoles
        dv(15,I) = dv(15,I) + pt5*(facd*dx*dz*dz + rm5*dpx - rm7*(Two*dx*dz*dpz + dpx*dz*dz + DddR*dx))   ! hxzz from d dipoles
        dv(16,I) = dv(16,I) + pt5*(facd*dy*dy*dy + rm5*Three*dpy - rm7*Three*(dpy*dy*dy + DddR*dy))       ! hyyy from d dipoles
        dv(17,I) = dv(17,I) + pt5*(facd*dy*dy*dz + rm5*dpz - rm7*(Two*dy*dz*dpy + dpz*dy*dy + DddR*dz))   ! hyyz from d dipoles
        dv(18,I) = dv(18,I) + pt5*(facd*dy*dz*dz + rm5*dpy - rm7*(Two*dy*dz*dpz + dpy*dz*dz + DddR*dy))   ! hyzz from d dipoles
        dv(19,I) = dv(19,I) + pt5*(facd*dz*dz*dz + rm5*Three*dpz - rm7*Three*(dpz*dz*dz + DddR*dz))       ! hzzz from d dipoles
        
        
        ! An alternative formulation:
!        ! Intermediates
!        px  = scalefd*ipd( 1, J, 1) + scalefp*ipd( 1, J, 2)
!        py  = scalefd*ipd( 2, J, 1) + scalefp*ipd( 2, J, 2)
!        pz  = scalefd*ipd( 3, J, 1) + scalefp*ipd( 3, J, 2)
!        DdR = px*dx + py*dy + pz*dz
!        fac = DdR*rm5
!        
!        ! Field
!        dv(1,I) = dv(1,I) + fac*dx - px*rm3   ! ex
!        dv(2,I) = dv(2,I) + fac*dy - py*rm3   ! ey
!        dv(3,I) = dv(3,I) + fac*dz - pz*rm3   ! ez
!        
!        ! Field gradients
!        fac = DdR*rm7
!        dv(4,I) = dv(4,I) + fac*dx*dx-(DdR+Two*px*dx)*rm5    ! gxx
!        dv(5,I) = dv(5,I) + fac*dx*dy-(px*dy+py*dx)*rm5      ! gxy
!        dv(6,I) = dv(6,I) + fac*dx*dz-(px*dz+pz*dx)*rm5      ! gxz
!        dv(7,I) = dv(7,I) + fac*dy*dy-(DdR+Two*py*dy)*rm5    ! gyy
!        dv(8,I) = dv(8,I) + fac*dy*dz-(py*dz+pz*dy)*rm5      ! gyz
!        dv(9,I) = dv(9,I) + fac*dz*dz-(DdR+two*pz*dz)*rm5    ! gzz
!        
!        ! Field hessians
!        fac = DdR*rm9
!        dv(10,I) = dv(10,I) + fac*dx*dx*dx + rm5*Three*px - rm7*Three*(px*dx*dx + DdR*dx)      ! hxxx
!        dv(11,I) = dv(11,I) + fac*dx*dx*dy + rm5*py - rm7*(Two*dx*dy*px + py*dx*dx + DdR*dy)   ! hxxy
!        dv(11,I) = dv(11,I) + fac*dx*dx*dz + rm5*pz - rm7*(Two*dx*dz*px + pz*dx*dx + DdR*dz)   ! hxxz
!        dv(13,I) = dv(13,I) + fac*dx*dy*dy + rm5*px - rm7*(Two*dx*dy*py + px*dy*dy + DdR*dx)   ! hxyy
!        dv(14,I) = dv(14,I) + fac*dx*dy*dz - rm7*(dx*dy*pz + dy*dz*px + dz*dx*py)              ! hxyz
!        dv(15,I) = dv(15,I) + fac*dx*dz*dz + rm5*px - rm7*(Two*dx*dz*pz + px*dz*dz + DdR*dx)   ! hxzz
!        dv(16,I) = dv(16,I) + fac*dy*dy*dy + rm5*Three*py - rm7*Three*(py*dy*dy + DdR*dy)      ! hyyy
!        dv(17,I) = dv(17,I) + fac*dy*dy*dz + rm5*pz - rm7*(Two*dy*dz*py + pz*dy*dy + DdR*dz)   ! hyyz
!        dv(18,I) = dv(18,I) + fac*dy*dz*dz + rm5*py - rm7*(Two*dy*dz*pz + py*dz*dz + DdR*dy)   ! hyzz
!        dv(19,I) = dv(19,I) + fac*dz*dz*dz + rm5*Three*pz - rm7*Three*(pz*dz*dz + DdR*dz)      ! hzzz
!        

    end if
end subroutine potential_deriv_D2M


subroutine field_deriv_D2D(scalefd,scalefp,I,J,de) !Amoeba,
    use mmpol
    implicit none
    !
    ! Computes electric field derivatives at position of polarizable cpol
    ! atom I generated by induced dipoles of cpol atom J. The computed 
    ! field is added to the field de
    ! 
    ! - For AMBER force-field only electric field derivatives at polarizable
    !   site is computed.
    !
    ! - For AMOEBA force-field electric field derivatives at p dipoles are
    !   calculated because only these dipoles can be polarized by induced 
    !   dipoles. Also only p dipoles contribute to this field
    !
    ! - Wang's polynomial damping for of the coulomb potential is used for
    !   the electric field generated by induced dipoles at polarizable 
    !   site for both AMOEBA and AMBER polarizable force-fields 
    !
    ! - The total electric field for AMBER dipoles and AMOEBA p dipoles is 
    !   scaled by factor scalefp. For the AMOEBA d dipoles the electric 
    !   field are scaled by factor scalefd. For the interaction without
    !   shielding scalefd = 1 and scalefp = 1
    !
    
    real(rp), dimension(6,pol_atoms,n_ipd), intent(inout) :: de
    real(rp), intent(in)    :: scalefd, scalefp
    !logical, intent(in)      :: Amoeba
    integer(ip), intent(in)  :: I, J
    !
    real(rp)    :: x, y, z, dx, dy, dz, rm1, rm3, rm5, rm7, rm9, rm11
    real(rp)    :: px, py, pz, DdR
    !
    !real(rp), parameter :: Two = 2.0_rp, Three = 3.0_rp, Four = 4.0_rp
    
    dx   = cpol(1,I) - cpol(1,J)
    dy   = cpol(2,I) - cpol(2,J)
    dz   = cpol(3,I) - cpol(3,J)
    rm1  = 0.0_rp
    rm3  = 0.0_rp
    rm5  = 0.0_rp
    rm7  = 0.0_rp
    rm9  = 0.0_rp
    rm11 = 0.0_rp
    
    ! Check if AMBER or AMOEBA polarizable force-field is used
    if (.not. Amoeba) then        ! AMBER

        call coulomb_kernel(.true.,3,dx,dy,dz,thole( polar_mm(I) ),thole( polar_mm(J) ),rm1,rm3,rm5,rm7,rm9,rm11)
        rm1 = rm1
        rm3 = rm3
        rm5 = Three*rm5
        rm7 = f15*rm7
        
        
!        rm5 = Three*Scr(IndJ)*fl5/r5
!        rm7 = F15*Scr(IndJ)*DdR*fl7/r7
!        fGxx = fGxx - rm5*(DdR + Two*dx*Dip(1,J))     + dx*dx*rm7
!        fGxy = fGxy - rm5*(dx*Dip(2,J) + dy*Dip(1,J)) + dx*dy*rm7
!        fGyy = fGyy - rm5*(DdR + Two*dy*Dip(2,J))     + dy*dy*rm7
!        fGxz = fGxz - rm5*(dx*Dip(3,J) + dz*Dip(1,J)) + dx*dz*rm7
!        fGyz = fGyz - rm5*(dy*Dip(3,J) + dz*Dip(2,J)) + dy*dz*rm7
!        fGzz = fGzz - rm5*(DdR + Two*dz*Dip(3,J))     + dz*dz*rm7
        
        
        ! intermediates
        px  = ipd( 1, J, 1)
        py  = ipd( 2, J, 1)
        pz  = ipd( 3, J, 1)
        DdR = px*dx + py*dy + pz*dz
          
        ! Field derivatives
        de(1,I,1) = de(1,I,1) + scalefp*(-rm5*(DdR + Two*dx*px) + dx*dx*DdR*rm7)   ! gxx
        de(2,I,1) = de(2,I,1) + scalefp*(-rm5*(dx*py + dy*px) + dx*dy*DdR*rm7)     ! gxy
        de(3,I,1) = de(3,I,1) + scalefp*(-rm5*(DdR + Two*dy*py) + dy*dy*DdR*rm7)   ! gyy
        de(4,I,1) = de(4,I,1) + scalefp*(-rm5*(dx*pz + dz*px) + dx*dz*DdR*rm7)     ! gxz
        de(5,I,1) = de(5,I,1) + scalefp*(-rm5*(dy*pz + dz*py) + dy*dz*DdR*rm7)     ! gyz
        de(6,I,1) = de(6,I,1) + scalefp*(-rm5*(DdR + Two*dz*pz) + dz*dz*DdR*rm7)   ! gzz
        
    elseif (Amoeba) then   ! AMOEBA

        call coulomb_kernel(.true.,3,dx,dy,dz,thole( polar_mm(I) ),thole( polar_mm(J) ),rm1,rm3,rm5,rm7,rm9,rm11)
        rm1 = rm1
        rm3 = rm3
        rm5 = Three*rm5
        rm7 = f15*rm7
        
        ! intermediates only p dipoles
        px  = ipd( 1, J, 2)
        py  = ipd( 2, J, 2)
        pz  = ipd( 3, J, 2)
        DdR = px*dx + py*dy + pz*dz
        
        ! Field derivatives p dipoles
        de(1,I,2) = de(1,I,2) + scalefp*(-rm5*(DdR + Two*dx*px) + dx*dx*DdR*rm7)   ! gxx
        de(2,I,2) = de(2,I,2) + scalefp*(-rm5*(dx*py + dy*px) + dx*dy*DdR*rm7)     ! gxy
        de(3,I,2) = de(3,I,2) + scalefp*(-rm5*(DdR + Two*dy*py) + dy*dy*DdR*rm7)   ! gyy
        de(4,I,2) = de(4,I,2) + scalefp*(-rm5*(dx*pz + dz*px) + dx*dz*DdR*rm7)     ! gxz
        de(5,I,2) = de(5,I,2) + scalefp*(-rm5*(dy*pz + dz*py) + dy*dz*DdR*rm7)     ! gyz
        de(6,I,2) = de(6,I,2) + scalefp*(-rm5*(DdR + Two*dz*pz) + dz*dz*DdR*rm7)   ! gzz       
        
        ! intermediates only d dipoles
        px  = ipd( 1, J, 1)
        py  = ipd( 2, J, 1)
        pz  = ipd( 3, J, 1)
        DdR = px*dx + py*dy + pz*dz

        ! Field derivatives d dipoles
        de(1,I,1) = de(1,I,1) + scalefd*(-rm5*(DdR + Two*dx*px) + dx*dx*DdR*rm7)   ! gxx
        de(2,I,1) = de(2,I,1) + scalefd*(-rm5*(dx*py + dy*px) + dx*dy*DdR*rm7)     ! gxy
        de(3,I,1) = de(3,I,1) + scalefd*(-rm5*(DdR + Two*dy*py) + dy*dy*DdR*rm7)   ! gyy
        de(4,I,1) = de(4,I,1) + scalefd*(-rm5*(dx*pz + dz*px) + dx*dz*DdR*rm7)     ! gxz
        de(5,I,1) = de(5,I,1) + scalefd*(-rm5*(dy*pz + dz*py) + dy*dz*DdR*rm7)     ! gyz
        de(6,I,1) = de(6,I,1) + scalefd*(-rm5*(DdR + Two*dz*pz) + dz*dz*DdR*rm7)   ! gzz
        
        
    end if
end subroutine field_deriv_D2D


subroutine potential_M2Q(scalef,I,J,v)
    use mmpol
    implicit none
    !
    ! Computes potential, electric field and field gradients at position 
    ! of qm atom I generated by multipoles of mm atom J. The computed
    ! potential is added to the potential v
    ! 

    real(rp), dimension(qm_atoms), intent(inout) :: v
    !logical,intent(in)      :: Amoeba
    integer(ip), intent(in) :: I, J
    real(rp), intent(in)    :: scalef
    !
    real(rp)    :: x, y, z, dx, dy, dz, rm1, rm3, rm5, rm7, rm9, rm11
    real(rp)    :: qq, px, py, pz, qxx, qxy, qxz, qyy, qyz, qzz, DdR, QRx, QRy, QRz, QRR, fac
    !
    
    dx   = cqm(1,I) - cmm(1,J)
    dy   = cqm(2,I) - cmm(2,J)
    dz   = cqm(3,I) - cmm(3,J)
    rm1  = 0.0_rp
    rm3  = 0.0_rp
    rm5  = 0.0_rp
    rm7  = 0.0_rp
    rm9  = 0.0_rp
    rm11 = 0.0_rp

    ! If multipoles q are only charges calculate only potential
    if (.not. Amoeba) then

        call coulomb_kernel(.false.,0,dx,dy,dz,0.0_rp,0.0_rp,rm1,rm3,rm5,rm7,rm9,rm11)
        rm1 = scalef*rm1
        
        v(I) = v(I) + q( 1,J)*rm1

    ! Else the potential electric field and field gradients are calculated
    elseif (Amoeba) then
        call coulomb_kernel(.false.,4,dx,dy,dz,0.0_rp,0.0_rp,rm1,rm3,rm5,rm7,rm9,rm11)
        !
        ! Scale distances
        !
        rm1 = scalef*rm1
        rm3 = scalef*rm3
        rm5 = scalef*Three*rm5
        !
        ! Set multipole moments 
        !
        qq  = q( 1,J)
        px  = q( 2,J)
        py  = q( 3,J)
        pz  = q( 4,J)
        qxx = q( 5,J)
        qxy = q( 6,J)
        qyy = q( 7,J)
        qxz = q( 8,J)
        qyz = q( 9,J)
        qzz = q(10,J)

        ! intermediates
        DdR  = px*dx  + py*dy  + pz*dz
        QRx  = qxx*dx + qxy*dy + qxz*dz
        QRy  = qxy*dx + qyy*dy + qyz*dz
        QRz  = qxz*dx + qyz*dy + qzz*dz
        QRR  = QRx*dx + QRy*dy + QRz*dz

        ! Potential
        v(I) = v(I) + (qq*rm1 + DdR*rm3 + QRR*rm5)     ! q

    end if
end subroutine potential_M2Q

subroutine field_M2Q(scalef,I,J,e)
    use mmpol
    implicit none
    !
    ! Computes electric field at position of polarizable cqm atom I
    ! generated by static multipoles of mm atom J. The computed 
    ! field is added to the field e
    !
    
    real(rp), dimension(3,qm_atoms), intent(inout) :: e
    real(rp), intent(in)    :: scalef
    !logical,intent(in)      :: Amoeba
    integer(ip),intent(in)  :: I, J
    !
    real(rp)    :: x, y, z, dx, dy, dz, rm1, rm3, rm5, rm7, rm9, rm11
    real(rp)    :: qq, px, py, pz, qxx, qxy, qxz, qyy, qyz, qzz, DdR, QRx, QRy, QRz, QRR, fac
    !
    !real(rp), parameter :: Two = 2.0_rp
            
    dx  = cqm(1,I) - cmm(1,J)
    dy  = cqm(2,I) - cmm(2,J)
    dz  = cqm(3,I) - cmm(3,J)
    rm1  = 0.0_rp
    rm3  = 0.0_rp
    rm5  = 0.0_rp
    rm7  = 0.0_rp
    rm9  = 0.0_rp
    rm11 = 0.0_rp

    ! Check if AMBER or AMOEBA polarizable force-field is used
    if (.not. Amoeba) then        ! AMBER

        call coulomb_kernel(.false.,1,dx,dy,dz,0.0,0.0,rm1,rm3,rm5,rm7,rm9,rm11)

        ! Field
        fac  = q( 1,J)*rm3
        e(1,I) = e(1,I) + scalef*fac*dx
        e(2,I) = e(2,I) + scalef*fac*dy
        e(3,I) = e(3,I) + scalef*fac*dz

    elseif (Amoeba) then   ! AMOEBA

        call coulomb_kernel(.false.,3,dx,dy,dz,0.0,0.0,rm1,rm3,rm5,rm7,rm9,rm11)
        !
        ! Scale distances
        !
        rm5 = Three*rm5
        rm7 = f15*rm7
        !
        ! Set multipole moments 
        !
        qq  = q( 1,J)
        px  = q( 2,J)
        py  = q( 3,J)
        pz  = q( 4,J)
        qxx = q( 5,J)
        qxy = q( 6,J)
        qyy = q( 7,J)
        qxz = q( 8,J)
        qyz = q( 9,J)
        qzz = q(10,J)

        ! intermediates
        DdR  = px*dx  + py*dy  + pz*dz
        QRx  = qxx*dx + qxy*dy + qxz*dz
        QRy  = qxy*dx + qyy*dy + qyz*dz
        QRz  = qxz*dx + qyz*dy + qzz*dz
        QRR  = QRx*dx + QRy*dy + QRz*dz

        ! Field
        fac  = qq*rm3 + DdR*rm5 + QRR*rm7
        e(1,I) = e(1,I) + scalef*(fac*dx - px*rm3 - Two*QRx*rm5)
        e(2,I) = e(2,I) + scalef*(fac*dy - py*rm3 - Two*QRy*rm5)
        e(3,I) = e(3,I) + scalef*(fac*dz - pz*rm3 - Two*QRz*rm5)
    end if
end subroutine field_M2Q


subroutine potential_D2Q(scalefd,scalefp,I,J,v)
    use mmpol
    implicit none
    !
    ! Computes potential at position of qm atom I generated by induced dipoles
    ! of cpol atom J. The computed potential is added to the potential v.
    ! 

    
    real(rp), dimension(qm_atoms), intent(inout) :: v
    !logical,intent(in)      :: Amoeba
    integer(ip), intent(in) :: I, J
    real(rp), intent(in)    :: scalefd,scalefp
    !
    real(rp)    :: x, y, z, dx, dy, dz, rm1, rm3, rm5, rm7, rm9, rm11
    real(rp)    :: dpx, dpy, dpz, ppx, ppy, ppz, px, py, pz, DdR, DpdR, DddR
    real(rp)    :: facd, facp, fac
    !
    !real(rp), parameter :: Two = 2.0_rp, Four = 4.0_rp
    
    dx   = cqm(1,I) - cpol(1,J)
    dy   = cqm(2,I) - cpol(2,J)
    dz   = cqm(3,I) - cpol(3,J)
    rm1  = 0.0_rp
    rm3  = 0.0_rp
    rm5  = 0.0_rp
    rm7  = 0.0_rp
    rm9  = 0.0_rp
    rm11 = 0.0_rp

    ! If multipoles q are only charges calculate only potential
    if (.not. Amoeba) then

        call coulomb_kernel(.false.,1,dx,dy,dz,0.0,0.0,rm1,rm3,rm5,rm7,rm9,rm11)
        
        ! Scaling of the interaction
        rm3 = scalefp*rm3
        
        ! Intermediates
        px  = ipd( 1, J, 1)
        py  = ipd( 2, J, 1)
        pz  = ipd( 3, J, 1)
        DdR  = px*dx  + py*dy  + pz*dz
        
        ! Potential of the induced dipoles
        v(I) = v(I) + DdR*rm3
        
    ! Else the potential electric field and field gradients are calculated
    elseif (Amoeba) then
        call coulomb_kernel(.false.,3,dx,dy,dz,0.0,0.0,rm1,rm3,rm5,rm7,rm9,rm11)
        
        rm1 = rm1
        rm3 = rm3
        rm5 = Three*rm5
        rm7 = f15*rm7
        rm9 = f105*rm9
        
        
        ! Intermediates p and d dipoles
        dpx  = scalefd*ipd( 1, J, 1)
        dpy  = scalefd*ipd( 2, J, 1)
        dpz  = scalefd*ipd( 3, J, 1)
        ppx  = scalefp*ipd( 1, J, 2)
        ppy  = scalefp*ipd( 2, J, 2)
        ppz  = scalefp*ipd( 3, J, 2)
        DddR = dpx*dx  + dpy*dy  + dpz*dz   ! DdR for d induced dipoles 
        DpdR = ppx*dx  + ppy*dy  + ppz*dz   ! DdR for p induced dipoles 

        ! Potential
        v(I) = v(I) + DddR*rm3  ! Vd
        v(I) = v(I) + DpdR*rm3  ! Vp

        
    end if
end subroutine potential_D2Q


subroutine field_D2Q(scalefd,scalefp,I,J,e)
    use mmpol
    implicit none
    !
    ! Computes electric field at position of polarizable cqm atom I
    ! generated by polarizable cpol atom J. The computed field is added
    ! to the field e
    !
    ! For the AMOEBA force-field the electric field of the p and d dipoles
    ! is summed into single electric field
    ! 

    real(rp), dimension(3,pol_atoms), intent(inout) :: e
    real(rp), intent(in)    :: scalefd, scalefp
    !logical,intent(in)      :: Amoeba
    integer(ip),intent(in)  :: I, J
    !
    real(rp)    :: x, y, z, dx, dy, dz, rm1, rm3, rm5, rm7, rm9, rm11
    real(rp)    :: px, py, pz, DdR
    
    dx   = cqm(1,I) - cpol(1,J)
    dy   = cqm(2,I) - cpol(2,J)
    dz   = cqm(3,I) - cpol(3,J)
    rm1  = 0.0_rp
    rm3  = 0.0_rp
    rm5  = 0.0_rp
    rm7  = 0.0_rp
    rm9  = 0.0_rp
    rm11 = 0.0_rp

    ! If multipoles q are only charges calculate only potential
    if (.not. Amoeba) then

        call coulomb_kernel(.false.,2,dx,dy,dz,0.0,0.0,rm1,rm3,rm5,rm7,rm9,rm11)
        rm1 = rm1
        rm3 = rm3
        rm5 = Three*rm5
        
        ! Intermediates
        px  = ipd( 1, J, 1)
        py  = ipd( 2, J, 1)
        pz  = ipd( 3, J, 1)
        DdR  = px*dx  + py*dy  + pz*dz
        
        ! Field from the induced dipoles
        e(1,I) = e(1,I) + scalefp*(DdR*rm5*dx - px*rm3)
        e(2,I) = e(2,I) + scalefp*(DdR*rm5*dy - py*rm3)
        e(3,I) = e(3,I) + scalefp*(DdR*rm5*dz - pz*rm3)
        
    elseif (Amoeba) then
        call coulomb_kernel(.false.,2,dx,dy,dz,0.0,0.0,rm1,rm3,rm5,rm7,rm9,rm11)
        rm1 = rm1
        rm3 = rm3
        rm5 = Three*rm5
        
        ! Intermediates p dipoles
        px  = ipd( 1, J, 2)
        py  = ipd( 2, J, 2)
        pz  = ipd( 3, J, 2)
        
        ! Field from the induced p dipoles
        DdR = px*dx + py*dy + pz*dz
        e(1,I) = e(1,I) + scalefp*(DdR*rm5*dx - px*rm3)
        e(2,I) = e(2,I) + scalefp*(DdR*rm5*dy - py*rm3)
        e(3,I) = e(3,I) + scalefp*(DdR*rm5*dz - pz*rm3)
        
        ! Intermediates d dipoles
        px  = ipd( 1, J, 1)
        py  = ipd( 2, J, 1)
        pz  = ipd( 3, J, 1)
        
        ! Field from the induced d dipoles
        DdR = px*dx + py*dy + pz*dz
        e(1,I) = e(1,I) + scalefd*(DdR*rm5*dx - px*rm3)
        e(2,I) = e(2,I) + scalefd*(DdR*rm5*dy - py*rm3)
        e(3,I) = e(3,I) + scalefd*(DdR*rm5*dz - pz*rm3)
        
        
    end if
end subroutine field_D2Q


subroutine multipoles_potential_remove(scr,v)
    use mmpol
    implicit none
    !
    ! Remove unwanted 1-5 contributions to the potential at MM sites 
    ! from the bonded atoms. 
    !
    ! scr  ...  0:  The sources are static multipoles of mm atoms  
    !               contained in the array q and with coordinates cmm
    !           1:  The sources are induced multipoles ipd, contained
    !               in the array ipd with coordinates cpol.
    ! 
    ! - For AMBER force-field only potential at MM sites is computed.
    !
    ! - For AMOEBA force-field potential, electric field and field gradients
    !   are computed. 
    !
    
    real(rp), dimension(ld_cart,mm_atoms), intent(inout) :: v
    integer(ip), intent(in) :: scr
    !
    !logical     :: Amoeba
    integer(ip) :: I, J, K, IJ
    real(rp)    :: scale
    !
    !real(rp), parameter :: Zero = 0.0_rp, One = 1.0_rp
    
    !Amoeba = ld_cart.gt.1
    
    if (scr.eq.0) then
        do I = 1,mm_atoms
            ! For every mm atom subtract interaction with it's neighbors 
            if(mscale(1).ne.one) then
                do IJ = 1, n12(I)
                    J = i12(IJ,I)
                    call potential_M2M(mscale(1)-One,I,J,v)
                enddo
            end if

            if(mscale(2).ne.one) then
                do IJ = 1, n13(I)
                    J = i13(IJ,I)
                    call potential_M2M(mscale(2)-One,I,J,v)
                enddo
            end if 

            if(mscale(3).ne.one) then
                do IJ = 1, n14(I)
                    J = i14(IJ,I)
                    call potential_M2M(mscale(3)-One,I,J,v)
                enddo
            end if 

            if(mscale(4).ne.one) then
                do IJ = 1, n15(I)
                    J = i15(IJ,I)
                    call potential_M2M(mscale(4)-One,I,J,v)
                enddo
            end if 
        enddo
    
! TODO: Check the potential of the induced dipoles if is correctly calculated (factors and if they are the same for p and d dipoles)
    elseif ((scr.eq.1).and.(Amoeba)) then
        
        ! For AMOEBA different scaling for the p and d dipoles
        do J = 1,pol_atoms
!            ! For every pol atom subtract interaction with it's neighbors for polarization field (p dipoles) 
!            if(pscale(1).ne.one) then
!                do IJ = 1, n12(polar_mm(J))
!                    I = i12(IJ,polar_mm(J))
!                    call potential_D2M(Zero,pscale(1)-One,I,J,v)
!                enddo
!            endif 
!            
!            if(pscale(2).ne.one) then
!                do IJ = 1, n13(polar_mm(J))
!                    I = i13(IJ,polar_mm(J))
!                    call potential_D2M(Zero,pscale(2)-One,I,J,v)
!                enddo
!            endif 
!            
!            if(pscale(3).ne.one) then
!                do IJ = 1, n14(polar_mm(J))
!                    I = i14(IJ,polar_mm(J))
!                    
!                    scale = pscale(3)
!                    do K = 1,np11(J)
!                        if (I.eq.ip11(K,J)) scale = pscale(3)*pscale(5)
!                    enddo
!                    
!                    call potential_D2M(Zero,scale-One,I,J,v)
!                enddo
!            endif 
!            
!            if(pscale(4).ne.one) then
!                do IJ = 1, n15(polar_mm(J))
!                    I = i15(IJ,polar_mm(J))
!                    call potential_D2M(Zero,pscale(4)-One,I,J,v)
!                enddo
!            endif 
            
            ! For every pol atom subtract interaction with it's neighbors for direct field (d dipoles) 
            if(dscale(1).ne.one) then
                do IJ = 1, np11(polar_mm(J))
                    I = ip11(IJ,polar_mm(J))
                    if (I.eq.polar_mm(J)) cycle
                    !call potential_D2M(dscale(1)-One,Zero,I,J,v)
                    call potential_D2M(Zero,dscale(1)-One,I,J,v)
                enddo
            endif 
            
            if(dscale(2).ne.one) then
                do IJ = 1, np12(polar_mm(J))
                    I = ip12(IJ,polar_mm(J))
                    if (I.eq.polar_mm(J)) cycle
                    !call potential_D2M(dscale(2)-One,Zero,I,J,v)
                    call potential_D2M(Zero,dscale(2)-One,I,J,v)
                enddo
            endif 
            
            if(dscale(3).ne.one) then
                do IJ = 1, np13(polar_mm(J))
                    I = ip13(IJ,polar_mm(J))
                    if (I.eq.polar_mm(J)) cycle
                    !call potential_D2M(dscale(3)-One,Zero,I,J,v)
                    call potential_D2M(Zero,dscale(3)-One,I,J,v)
                enddo
            endif 
            
            if(dscale(4).ne.one) then
                do IJ = 1, np14(polar_mm(J))
                    I = ip14(IJ,polar_mm(J))
                    if (I.eq.polar_mm(J)) cycle
                    !call potential_D2M(dscale(4)-One,Zero,I,J,v)
                    call potential_D2M(Zero,dscale(4)-One,I,J,v)
                enddo
            endif 
        enddo
    elseif ((scr.eq.1).and.(.not. Amoeba)) then
        
        ! For AMBER the potential from the dipoles is scaled by pscale (so far)
        do J = 1,pol_atoms
            if(pscale(1).ne.one) then
                do IJ = 1, n12(polar_mm(J))
                    I = i12(IJ,polar_mm(J))
                    call potential_D2M(Zero,pscale(1)-One,I,J,v)
                enddo
            endif 
            
            if(pscale(2).ne.one) then
                do IJ = 1, n13(polar_mm(J))
                    I = i13(IJ,polar_mm(J))
                    call potential_D2M(Zero,pscale(2)-One,I,J,v)
                enddo
            endif 
            
            if(pscale(3).ne.one) then
                do IJ = 1, n14(polar_mm(J))
                    I = i14(IJ,polar_mm(J))
                    call potential_D2M(Zero,pscale(3)-One,I,J,v)
                enddo
            endif 
            
            if(pscale(4).ne.one) then
                do IJ = 1, n15(polar_mm(J))
                    I = i15(IJ,polar_mm(J))
                    call potential_D2M(Zero,pscale(4)-One,I,J,v)
                enddo
            endif 
        enddo 
    end if 
end subroutine multipoles_potential_remove

subroutine multipoles_field_remove(scr,e)
    use mmpol
    implicit none
    !
    ! Remove unwanted 1-4 contributions to the electric field from
    ! the bonded atoms generated by the static multipoles of the MM
    ! atoms at polarizable atom sites. 
    ! 
    ! scr  ...  0:  The sources are static multipoles of mm atoms  
    !               contained in the array q and with coordinates cmm
    !           1:  The sources are induced multipoles ipd, contained
    !               in the array ipd with coordinates cpol.
    !
    ! - For AMBER force-field only electric field at polarizable site
    !   is computed.
    !
    ! - For AMOEBA force-field direct electric field from d dipoles and 
    !   polarization electric field from p dipoles are calculated
    !
    ! - Wang's polynomial damping for of the coulomb potential is used for
    !   the electric field at the polarizable site for both AMOEBA and 
    !   AMBER polarizable force-fields 
    !
    ! The function can be simplified by adding p scaling for both amoeba
    ! and AMBER dipoles and d scaling only to the AMOEBA (shorter form)
    !
    
    real(rp), dimension(3,pol_atoms,n_ipd), intent(inout) :: e
    integer(ip), intent(in) :: scr
    !
    !logical     :: Amoeba
    integer(ip) :: I, J, K, IJ
    real(rp)    :: scale
    !
    !real(rp), parameter :: Zero = 0.0_rp ,One = 1.0_rp
    
    
    !Amoeba = n_ipd.ge.2
    
    if ((scr.eq.0).and.(Amoeba)) then      ! AMOEBA and sources are MM atoms
        do I = 1,pol_atoms

            ! p field dipoles
            if (pscale(1).ne.one) then
                do IJ = 1, n12(polar_mm(I))
                    J = i12(IJ,polar_mm(I))
                    call field_M2D(Zero,pscale(1)-One,I,J,e)
                enddo
            end if

            if (pscale(2).ne.one) then
                do IJ = 1, n13(polar_mm(I))
                    J = i13(IJ,polar_mm(I))
                    call field_M2D(Zero,pscale(2)-One,I,J,e)
                enddo
            end if

            ! Strange way of rescaling 1-4 interactions - if 1-4 atom also in ip11 then scale by pscale(3)*pscale(5) else pscale(3)
            do IJ = 1, n14(polar_mm(I))
                J = i14(IJ,polar_mm(I))
                scale = pscale(3)
                do K = 1,np11(I)
                    if (J.eq.ip11(K,I)) scale = pscale(3)*pscale(5)
                enddo
                if (scale.ne.one) then
                    call field_M2D(Zero,scale-One,I,J,e)
                end if
            enddo
            

            if (pscale(4).ne.one) then
                do IJ = 1, n15(polar_mm(I))
                    J = i15(IJ,polar_mm(I))
                    call field_M2D(Zero,pscale(4)-One,I,J,e)
                enddo
            end if

            ! d field dipoles
            if (dscale(1).ne.One) then
                do IJ = 1, np11(polar_mm(I))
                    J = ip11(IJ,polar_mm(I))
                    if (polar_mm(I).eq.J) cycle
                    call field_M2D(dscale(1)-One,Zero,I,J,e)
                enddo
            end if 

            if (dscale(2).ne.one) then
                do IJ = 1, np12(polar_mm(I))
                    J = ip12(IJ,polar_mm(I))
                    if (polar_mm(I).eq.J) cycle
                    call field_M2D(dscale(2)-One,Zero,I,J,e)
                enddo
            end if 

            if (dscale(3).ne.one) then
                do IJ = 1, np13(polar_mm(I))
                    J = ip13(IJ,polar_mm(I))
                    if (polar_mm(I).eq.J) cycle
                    call field_M2D(dscale(3)-One,Zero,I,J,e)
                enddo
            end if 

            if (dscale(4).ne.one) then
                do IJ = 1, np14(polar_mm(I))
                    J = ip14(IJ,polar_mm(I))
                    if (polar_mm(I).eq.J) cycle
                    call field_M2D(dscale(4)-One,Zero,I,J,e)
                enddo
            end if 
        enddo
        
    elseif ((scr.eq.0).and.(.not. Amoeba)) then     ! AMBER and sources are MM atoms
        do I = 1,pol_atoms

            ! Field from dipoles is scaled by pscale
            if (pscale(1).ne.one) then
                do IJ = 1, n12(polar_mm(I))
                    J = i12(IJ,polar_mm(I))
                    call field_M2D(Zero,pscale(1)-One,I,J,e)
                enddo
            end if

            if (pscale(2).ne.one) then
                do IJ = 1, n13(polar_mm(I))
                    J = i13(IJ,polar_mm(I))
                    call field_M2D(Zero,pscale(2)-One,I,J,e)
                enddo
            end if

            if (pscale(3).ne.one) then
                do IJ = 1, n14(polar_mm(I))
                    J = i14(IJ,polar_mm(I))
                    call field_M2D(Zero,pscale(3)-One,I,J,e)
                enddo
            end if

            if (pscale(4).ne.one) then
                do IJ = 1, n15(polar_mm(I))
                    J = i15(IJ,polar_mm(I))
                    call field_M2D(Zero,pscale(4)-One,I,J,e)
                enddo
            end if
        enddo 
    
    elseif ((scr.eq.1).and.(Amoeba)) then      ! AMOEBA and sources are induced dipoles (only from p dipoles)
        do I = 1,pol_atoms
            if (uscale(1).ne.One) then
                do IJ = 1, np11( polar_mm(I) )
                    J = ip11(IJ, polar_mm(I) )

                    if ((mm_polar(J).eq.0).or.(polar_mm(I).eq.J)) cycle ! if J is not polarizable atom, skip J

                    call field_D2D(uscale(1)-One,uscale(1)-One,I,mm_polar(J),e)
                enddo
            end if 

            if (uscale(2).ne.one) then
                do IJ = 1, np12(polar_mm(I))
                    J = ip12(IJ,polar_mm(I))
                    
                    if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J
                    
                    call field_D2D(uscale(2)-One,uscale(2)-One,I,mm_polar(J),e)
                enddo
            end if 

            if (uscale(3).ne.one) then
                do IJ = 1, np13(polar_mm(I))
                    J = ip13(IJ,polar_mm(I))
                    
                    if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J
                    
                    call field_D2D(uscale(3)-One,uscale(3)-One,I,mm_polar(J),e)
                enddo
            end if 

            if (uscale(4).ne.one) then
                do IJ = 1, np14(polar_mm(I))
                    J = ip14(IJ,polar_mm(I))
                    
                    if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J
                    
                    call field_D2D(uscale(4)-One,uscale(4)-One,I,mm_polar(J),e)
                enddo
            end if 
        enddo
        
    elseif ((scr.eq.1).and.(.not. Amoeba)) then     ! AMBER and sources are induced dipoles
        do I = 1,pol_atoms

            ! Field from dipoles is scaled by uscale
            if (uscale(1).ne.one) then
                do IJ = 1, n12(polar_mm(I))
                    J = i12(IJ,polar_mm(I))
                    
                    if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J
                    
                    call field_D2D(Zero,uscale(1)-One,I, mm_polar(J) ,e)
                enddo
            end if

            if (uscale(2).ne.one) then
                do IJ = 1, n13(polar_mm(I))
                    J = i13(IJ,polar_mm(I))
                    
                    if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J
                    
                    call field_D2D(Zero,uscale(2)-One,I, mm_polar(J) ,e)
                enddo
            end if

            if (uscale(3).ne.one) then
                do IJ = 1, n14(polar_mm(I))
                    J = i14(IJ,polar_mm(I))
                    
                    if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J
                    
                    call field_D2D(Zero,uscale(3)-One,I, mm_polar(J) ,e)
                enddo
            end if

            if (uscale(4).ne.one) then
                do IJ = 1, n15(polar_mm(I))
                    J = i15(IJ,polar_mm(I))
                    
                    if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J
                    
                    call field_D2D(Zero,uscale(4)-One,I, mm_polar(J) ,e)
                enddo
            end if
        enddo 
        
    end if

end subroutine multipoles_field_remove

subroutine multipoles_potential_deriv_remove(scr,dv)
    use mmpol
    implicit none
    !
    ! Remove unwanted 1-5 contributions to the potential derivatives at 
    ! MM sites from the bonded atoms. 
    !
    ! scr  ...  0:  The sources are static multipoles of mm atoms  
    !               contained in the array q and with coordinates cmm
    !           1:  The sources are induced multipoles ipd, contained
    !               in the array ipd with coordinates cpol.
    ! 
    ! - For AMBER force-field only potential derivatives at MM sites 
    !   are computed.
    !
    ! - For AMOEBA force-field derivatives of the potential, electric 
    !   field and field gradients are computed. 
    !
    
    real(rp), dimension(ld_cder,mm_atoms), intent(inout) :: dv
    integer(ip), intent(in) :: scr
    !
    !logical     :: Amoeba
    integer(ip) :: I, J, K, IJ
    real(rp)    :: scale
    !
    !real(rp), parameter :: Zero = 0.0_rp, One = 1.0_rp
    
    !Amoeba = ld_cder.gt.3
    
    if (scr.eq.0) then
        do I = 1,mm_atoms
            ! For every mm atom subtract interaction with it's neighbors 
            if(mscale(1).ne.one) then
                do IJ = 1, n12(I)
                    J = i12(IJ,I)
                    call potential_deriv_M2M(mscale(1)-One,I,J,dv)
                enddo
            end if

            if(mscale(2).ne.one) then
                do IJ = 1, n13(I)
                    J = i13(IJ,I)
                    call potential_deriv_M2M(mscale(2)-One,I,J,dv)
                enddo
            end if 

            if(mscale(3).ne.one) then
                do IJ = 1, n14(I)
                    J = i14(IJ,I)
                    call potential_deriv_M2M(mscale(3)-One,I,J,dv)
                enddo
            end if 

            if(mscale(4).ne.one) then
                do IJ = 1, n15(I)
                    J = i15(IJ,I)
                    call potential_deriv_M2M(mscale(4)-One,I,J,dv)
                enddo
            end if 
        enddo
    
! TODO: Check the potential of the induced dipoles if is correctly calculated (factors and if they are the same for p and d dipoles)
    elseif ((scr.eq.1).and.(Amoeba)) then
        
        ! For AMOEBA different scaling for the p and d dipoles
        do J = 1,pol_atoms
            ! For every pol atom subtract interaction with it's neighbors for polarization field (p dipoles) 
            if(pscale(1).ne.one) then
                do IJ = 1, n12(polar_mm(J))
                    I = i12(IJ,polar_mm(J))
                    if (I.eq.polar_mm(J)) cycle
                    !call potential_deriv_D2M(Zero,pscale(1)-One,I,J,dv)
                    call potential_deriv_D2M(pscale(1)-One,Zero,I,J,dv)
                enddo
            endif 
            
            if(pscale(2).ne.one) then
                do IJ = 1, n13(polar_mm(J))
                    I = i13(IJ,polar_mm(J))
                    if (I.eq.polar_mm(J)) cycle
                    !call potential_deriv_D2M(Zero,pscale(2)-One,I,J,dv)
                    call potential_deriv_D2M(pscale(2)-One,Zero,I,J,dv)
                enddo
            endif 
            
            ! Strange way of rescaling 1-4 interactions - if 1-4 atom also in ip11 then scale by pscale(3)*pscale(5) else pscale(3)
            do IJ = 1, n14(polar_mm(J))
                I = i14(IJ,polar_mm(J))
                if (I.eq.polar_mm(J)) cycle
                scale = pscale(3)
                do K = 1,np11(polar_mm(J))
                    if (I.eq.ip11(K,polar_mm(J))) scale = pscale(3)*pscale(5)
                enddo
                if (scale.ne.one) then
                    !call potential_deriv_D2M(Zero,scale-One,I,J,dv)
                    call potential_deriv_D2M(scale-One,Zero,I,J,dv)
                end if
            enddo
!            if(pscale(3).ne.one) then
!                do IJ = 1, n14(polar_mm(J))
!                    I = i14(IJ,polar_mm(J))
!                    if (I.eq.polar_mm(J)) cycle
!                    call potential_deriv_D2M(Zero,pscale(3)-One,I,J,dv)
!                enddo
!            endif 
            
            if(pscale(4).ne.one) then
                do IJ = 1, n15(polar_mm(J))
                    I = i15(IJ,polar_mm(J))
                    if (I.eq.polar_mm(J)) cycle
                    !call potential_deriv_D2M(Zero,pscale(4)-One,I,J,dv
                    call potential_deriv_D2M(pscale(4)-One,Zero,I,J,dv)
                enddo
            endif 
            
            ! For every pol atom subtract interaction with it's neighbors for direct field (d dipoles) 
            if(dscale(1).ne.one) then
                do IJ = 1, np11(polar_mm(J))
                    I = ip11(IJ,polar_mm(J))
                    if (I.eq.polar_mm(J)) cycle
                    !call potential_deriv_D2M(dscale(1)-One,Zero,I,J,dv)
                    call potential_deriv_D2M(Zero,dscale(1)-One,I,J,dv)
                    
                    
                enddo
            endif 
            
            if(dscale(2).ne.one) then
                do IJ = 1, np12(polar_mm(J))
                    I = ip12(IJ,polar_mm(J))
                    if (I.eq.polar_mm(J)) cycle
                    !call potential_deriv_D2M(dscale(2)-One,Zero,I,J,dv)
                    call potential_deriv_D2M(Zero,dscale(2)-One,I,J,dv)
                enddo
            endif 
            
            if(dscale(3).ne.one) then
                do IJ = 1, np13(polar_mm(J))
                    I = ip13(IJ,polar_mm(J))
                    if (I.eq.polar_mm(J)) cycle
                    !call potential_deriv_D2M(dscale(3)-One,Zero,I,J,dv)
                    call potential_deriv_D2M(Zero,dscale(3)-One,I,J,dv)
                enddo
            endif 
            
            if(dscale(4).ne.one) then
                do IJ = 1, np14(polar_mm(J))
                    I = ip14(IJ,polar_mm(J))
                    if (I.eq.polar_mm(J)) cycle
                    !call potential_deriv_D2M(dscale(4)-One,Zero,I,J,dv)
                    call potential_deriv_D2M(Zero,dscale(4)-One,I,J,dv)
                enddo
            endif 
        enddo
    elseif ((scr.eq.1).and.(.not. Amoeba)) then
        
        ! For AMBER the potential from the dipoles is scaled by pscale (so far)
        do J = 1,pol_atoms
            if(pscale(1).ne.one) then
                do IJ = 1, n12(polar_mm(J))
                    I = i12(IJ,polar_mm(J))
                    call potential_deriv_D2M(Zero,pscale(1)-One,I,J,dv)
                enddo
            endif 
            
            if(pscale(2).ne.one) then
                do IJ = 1, n13(polar_mm(J))
                    I = i13(IJ,polar_mm(J))
                    call potential_deriv_D2M(Zero,pscale(2)-One,I,J,dv)
                enddo
            endif 
            
            if(pscale(3).ne.one) then
                do IJ = 1, n14(polar_mm(J))
                    I = i14(IJ,polar_mm(J))
                    call potential_deriv_D2M(Zero,pscale(3)-One,I,J,dv)
                enddo
            endif 
            
            if(pscale(4).ne.one) then
                do IJ = 1, n15(polar_mm(J))
                    I = i15(IJ,polar_mm(J))
                    call potential_deriv_D2M(Zero,pscale(4)-One,I,J,dv)
                enddo
            endif 
        enddo 
    end if 
end subroutine multipoles_potential_deriv_remove

subroutine multipoles_field_deriv_remove(scr,de)
    use mmpol
    implicit none
    !
    ! Remove unwanted 1-4 contributions to the electric field derivatives
    ! from the bonded atoms generated by the static multipoles of the MM
    ! atoms at polarizable atom sites. 
    ! 
    ! scr  ...  0:  The sources are static multipoles of mm atoms  
    !               contained in the array q and with coordinates cmm
    !           1:  The sources are induced multipoles ipd, contained
    !               in the array ipd with coordinates cpol.
    !
    ! - For AMBER force-field only electric derivatives field at polarizable
    !   site is computed.
    !
    ! - For AMOEBA force-field electric field derivatives for p and d dipoles
    !   are computed. 
    !
    ! - Wang's polynomial damping for of the coulomb potential is used for
    !   the electric field at the polarizable site for both AMOEBA and 
    !   AMBER polarizable force-fields 
    !
    ! The function can be simplified by adding p scaling for both amoeba
    ! and AMBER dipoles and d scaling only to the AMOEBA (shorter form)
    !
    ! Scaling of the nearest neighbor interactions for p and d dipoles is
    ! reversed dscale is used for scaling potential of p dipoles and pscale
    ! is used for scaling potential of d dipoles
    !
    
    real(rp), dimension(6,pol_atoms,n_ipd), intent(inout) :: de
    integer(ip), intent(in) :: scr
    !
    !logical     :: Amoeba
    integer(ip) :: I, J, K, IJ
    real(rp)    :: scale
    !
    !real(rp), parameter :: Zero = 0.0_rp ,One = 1.0_rp
    
    
    !Amoeba = n_ipd.ge.2
    
    if ((scr.eq.0).and.(Amoeba)) then      ! AMOEBA and sources are MM atoms
        do I = 1,pol_atoms 
            ! p field dipoles
            if (pscale(1).ne.one) then
                do IJ = 1, n12(polar_mm(I))
                    J = i12(IJ,polar_mm(I))
                    call field_deriv_M2D(Zero,pscale(1)-One,I,J,de)
                enddo
            end if

            if (pscale(2).ne.one) then
                do IJ = 1, n13(polar_mm(I))
                    J = i13(IJ,polar_mm(I))
                    call field_deriv_M2D(Zero,pscale(2)-One,I,J,de)
                enddo
            end if

            do IJ = 1, n14(polar_mm(I))
                J = i14(IJ,polar_mm(I))
                scale = pscale(3)
                do K = 1,np11(polar_mm(I))
                    if (J.eq.ip11(K,polar_mm(I))) scale = pscale(3)*pscale(5)
                enddo
                if (scale.ne.one) then
                    call field_deriv_M2D(Zero,scale-One,I,J,de)
                end if
            enddo
!            if (pscale(3).ne.one) then
!                do IJ = 1, n14(polar_mm(I))
!                    J = i14(IJ,polar_mm(I))
!                    call field_deriv_M2D(Zero,pscale(3)-One,I,J,de)
!                enddo
!            end if

            if (pscale(4).ne.one) then
                do IJ = 1, n15(polar_mm(I))
                    J = i15(IJ,polar_mm(I))
                    call field_deriv_M2D(Zero,pscale(4)-One,I,J,de)
                enddo
            end if


            ! d field dipoles
            if (dscale(1).ne.One) then
                do IJ = 1, np11(polar_mm(I))
                    J = ip11(IJ,polar_mm(I))
                    if (polar_mm(I).eq.J) cycle
                    call field_deriv_M2D(dscale(1)-One,Zero,I,J,de)
                enddo
            end if 

            if (dscale(2).ne.one) then
                do IJ = 1, np12(polar_mm(I))
                    J = ip12(IJ,polar_mm(I))
                    if (polar_mm(I).eq.J) cycle
                    call field_deriv_M2D(dscale(2)-One,Zero,I,J,de)
                enddo
            end if 

            if (dscale(3).ne.one) then
                do IJ = 1, np13(polar_mm(I))
                    J = ip13(IJ,polar_mm(I))
                    if (polar_mm(I).eq.J) cycle
                    call field_deriv_M2D(dscale(3)-One,Zero,I,J,de)
                enddo
            end if 

            if (dscale(4).ne.one) then
                do IJ = 1, np14(polar_mm(I))
                    J = ip14(IJ,polar_mm(I))
                    if (polar_mm(I).eq.J) cycle
                    call field_deriv_M2D(dscale(4)-One,Zero,I,J,de)
                enddo
            end if 
        enddo
        
    elseif ((scr.eq.0).and.(.not. Amoeba)) then     ! AMBER and sources are MM atoms
        do I = 1,pol_atoms

            ! Field from dipoles is scaled by pscale
            if (pscale(1).ne.one) then
                do IJ = 1, n12(polar_mm(I))
                    J = i12(IJ,polar_mm(I))
                    call field_deriv_M2D(Zero,pscale(1)-One,I,J,de)
                enddo
            end if

            if (pscale(2).ne.one) then
                do IJ = 1, n13(polar_mm(I))
                    J = i13(IJ,polar_mm(I))
                    call field_deriv_M2D(Zero,pscale(2)-One,I,J,de)
                enddo
            end if

            if (pscale(3).ne.one) then
                do IJ = 1, n14(polar_mm(I))
                    J = i14(IJ,polar_mm(I))
                    call field_deriv_M2D(Zero,pscale(3)-One,I,J,de)
                enddo
            end if

            if (pscale(4).ne.one) then
                do IJ = 1, n15(polar_mm(I))
                    J = i15(IJ,polar_mm(I))
                    call field_deriv_M2D(Zero,pscale(4)-One,I,J,de)
                enddo
            end if
        enddo 
    
    elseif ((scr.eq.1).and.(Amoeba)) then      ! AMOEBA and sources are induced dipoles (only from p dipoles)
        do I = 1,pol_atoms
            if (uscale(1).ne.One) then
                do IJ = 1, np11(polar_mm(I))
                    J = ip11(IJ,polar_mm(I))

                    if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J

                    call field_deriv_D2D(uscale(1)-One,uscale(1)-One,I,mm_polar(J),de) !Amoeba,
                enddo
            end if 

            if (uscale(2).ne.one) then
                do IJ = 1, np12(polar_mm(I))
                    J = ip12(IJ,polar_mm(I))
                    
                    if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J
                    
                    call field_deriv_D2D(uscale(2)-One,uscale(2)-One,I,mm_polar(J),de) !Amoeba,
                enddo
            end if 

            if (uscale(3).ne.one) then
                do IJ = 1, np13(polar_mm(I))
                    J = ip13(IJ,polar_mm(I))
                    
                    if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J
                    
                    call field_deriv_D2D(uscale(3)-One,uscale(3)-One,I,mm_polar(J),de) !Amoeba,
                enddo
            end if 

            if (uscale(4).ne.one) then
                do IJ = 1, np14(polar_mm(I))
                    J = ip14(IJ,polar_mm(I))
                    
                    if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J
                    
                    call field_deriv_D2D(uscale(4)-One,uscale(4)-One,I,mm_polar(J),de) !Amoeba,
                enddo
            end if 
        enddo
        
    elseif ((scr.eq.1).and.(.not. Amoeba)) then     ! AMBER and sources are induced dipoles
        do I = 1,pol_atoms

            ! Field from dipoles is scaled by uscale
            if (uscale(1).ne.one) then
                do IJ = 1, n12(polar_mm(I))
                    J = i12(IJ,polar_mm(I))
                    
                    if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J
                    
                    call field_deriv_D2D(Zero,uscale(1)-One,I,mm_polar(J),de) !Amoeba,
                enddo
            end if

            if (uscale(2).ne.one) then
                do IJ = 1, n13(polar_mm(I))
                    J = i13(IJ,polar_mm(I))
                    
                    if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J
                    
                    call field_deriv_D2D(Zero,uscale(2)-One,I,mm_polar(J),de) !Amoeba,
                enddo
            end if

            if (uscale(3).ne.one) then
                do IJ = 1, n14(polar_mm(I))
                    J = i14(IJ,polar_mm(I))
                    
                    if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J
                    
                    call field_deriv_D2D(Zero,uscale(3)-One,I,mm_polar(J),de)
                enddo
            end if

            if (uscale(4).ne.one) then
                do IJ = 1, n15(polar_mm(I))
                    J = i15(IJ,polar_mm(I))
                    
                    if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J
                    
                    call field_deriv_D2D(Zero,uscale(4)-One,I,mm_polar(J),de)
                enddo
            end if
        enddo 
        
    end if

end subroutine multipoles_field_deriv_remove
    

subroutine multipoles_potential_FMM(scr,v)
!
end subroutine multipoles_potential_FMM
    
end module elstat


! TODO:
! Check if np11, np12,... is defined with polarizable atom indexes or with mm atom indexes
! check for _D2D quatities that we always include porazizable atom indexes and not mm atom ones 
!
! QUESTIONS:
! What should be outut of potential calculation from induced dipoles? 
!    - for interaction energy it should be potential of p dipoles ipd(:,:,2) with scaling scaled
!      and exclusion neighbors np1..
!
! Scaling of the nearest neighbor interactions for p and d dipoles is
! reversed dscale is used for scaling potential of p dipoles and pscale
! is used for scaling potential of d dipoles
!   - Why?
!
! For dE D@D both d and p dipoles are scaled by uscale and excluded parameters are both in np1..
!
    !
    ! Scaling of the nearest neighbor interactions for p and d dipoles is
    ! reversed dscale is used for scaling potential of p dipoles and pscale
    ! is used for scaling potential of d dipoles
    !