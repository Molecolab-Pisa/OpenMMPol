!-----------------------------------------------------------
!   Not finished
!-----------------------------------------------------------
module polar

contains

subroutine dipole_T(scalef,I,J,TTens)
    use precision
    use mmpol, only : thole, cpol, pol_atoms, pol
    implicit none
    !                      
    ! Compute element of the polarization tensor TTens between
    ! polarizable cpol atom I and polarizable cpol atom J. On the
    ! TTens diagonal (I=J) are inverse polarizabilities and on the 
    ! off-diagonal dipole field.
    !
    real(rp), dimension(3,pol_atoms,3,pol_atoms), intent(inout) :: TTens
    real(rp), intent(in) :: scalef
    integer(ip), intent(in) :: I,J
    !
    real(rp)    :: dx, dy, dz, rm1, rm3, rm5, rm7, rm9, rm11
    !
    real(rp), parameter   :: Zero = 0.0_rp, One = 1.0_rp, Three = 3.0_rp
    !
     
    ! Initialize variables
    rm1 = Zero
    rm3 = Zero
    rm5 = Zero
    rm7 = Zero
    rm9 = Zero
    rm11 = Zero
    
    if (I.eq.J) then
        TTens(1,I,1,I) = One/pol(I)
        TTens(2,I,2,I) = One/pol(I)
        TTens(3,I,3,I) = One/pol(I)
    else
        
        ! intermediates
        dx = cpol(1,J) - cpol(1,I)
        dy = cpol(2,J) - cpol(2,I)
        dz = cpol(3,J) - cpol(3,I)

        call coulomb_kernel(.true.,2,dx,dy,dz,thole(I),thole(J),rm1,rm3,rm5,rm7,rm9,rm11)

        ! Apply scaling
        rm3 = rm3*scalef
        rm5 = rm5*scalef

        ! Fill the matrix elemets
        TTens(1,J,1,I) = TTens(1,J,1,I) + rm3 - Three*rm5*dx*dx
        TTens(2,J,1,I) = TTens(2,J,1,I) - Three*rm5*dx*dy
        TTens(3,J,1,I) = TTens(3,J,1,I) - Three*rm5*dx*dz

        TTens(1,J,2,I) = TTens(1,J,2,I) - Three*rm5*dy*dx
        TTens(2,J,2,I) = TTens(2,J,2,I) + rm3 - Three*rm5*dy*dy
        TTens(3,J,2,I) = TTens(3,J,2,I) - Three*rm5*dy*dz

        TTens(1,J,3,I) = TTens(1,J,3,I) - Three*rm5*dz*dx
        TTens(2,J,3,I) = TTens(2,J,3,I) - Three*rm5*dz*dy
        TTens(3,J,3,I) = TTens(3,J,3,I) + rm3 - Three*rm5*dz*dz
    end if
end subroutine dipole_T
    
    

subroutine create_TMat(TMat)
    use precision
    use mmpol, only : pol_atoms, cpol, uscale, np11, np12, np13, np14, &
                           ip11, ip12, ip13, ip14, n12, n13, n14, n15, &
                           i12, i13, i14, i15, pol, mm_polar, Amoeba
    implicit none
    !                      
    ! Construct polarization tensor for matrix inversion solution
    ! of the induced dipoles by matrix inversion. On the diagonal
    ! of the TMat there are inverse polarizabilities and on the 
    ! off-diagonal dipole field.
    !
    ! The tensor is constructed with dimension (3,pol_atoms,3,pol_atoms)
    ! then we subtract the excluded interaction and reshape the tensor
    ! into matrix with dimension (3*pol_atoms,3*pol_atoms) which is 
    ! the one on the output.
    !
    real(rp), dimension(3*pol_atoms,3*pol_atoms), intent(out) :: TMat
    !
    real(rp), dimension(3,pol_atoms,3,pol_atoms) :: TTens
    integer(ip) :: I, J, K, L, IJ
    !
    real(rp), parameter   :: Zero = 0.0_rp, One = 1.0_rp
    !
    
    ! Initialize the tensor with zeros
    TTens = Zero
    
    ! Assemble the diagonal elements of the MMpol matrix
!    do I = 1,pol_atoms
!        do J = 1,3
!            TTens(J,I,J,I) = One/pol(I)
!        enddo
!    enddo
    
    ! Assemble the upper triangular MMPol matrix and diagonal
    do I = 1,pol_atoms
        do J = 1,I !-1
            
            call dipole_T(One,I,J,TTens)
            
        enddo
    enddo
    
    ! Fill the lower diagonal
    do I = 1, pol_atoms
        do J = I+1, pol_atoms
          do K = 1, 3
            do L = 1, 3
                TTens(L,J,K,I) = TTens(K,I,L,J)
            enddo
          enddo
        enddo
    enddo
    
    ! Remove excluded interactions
    if  (Amoeba) then      ! AMOEBA and sources are induced dipoles (only from p dipoles)
        do I = 1,pol_atoms
            if (uscale(1).ne.One) then
                do IJ = 1, np11(I)
                    J = ip11(IJ,I)

                    if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J

                    call dipole_T(uscale(1)-One,I,J,TTens)

                enddo
            end if 

            if (uscale(2).ne.one) then
                do IJ = 1, np12(I)
                    J = ip12(IJ,I)
                    
                    if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J
                    
                    call dipole_T(uscale(2)-One,I,J,TTens)
                enddo
            end if 

            if (uscale(3).ne.one) then
                do IJ = 1, np13(I)
                    J = ip13(IJ,I)
                    
                    if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J
                    
                    call dipole_T(uscale(3)-One,I,J,TTens)
                    
                enddo
            end if 

            if (uscale(4).ne.one) then
                do IJ = 1, np14(I)
                    J = ip14(IJ,I)
                    
                    if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J
                    
                    call dipole_T(uscale(4)-One,I,J,TTens)

                enddo
            end if 
        enddo
        
    else     ! AMBER and sources are induced dipoles
        do I = 1,pol_atoms

            ! Field from dipoles is scaled by uscale
            if (uscale(1).ne.one) then
                do IJ = 1, n12(I)
                    J = i12(IJ,I)
                    
                    if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J
                    
                    call dipole_T(uscale(1)-One,I,J,TTens)

                enddo
            end if

            if (uscale(2).ne.one) then
                do IJ = 1, n13(I)
                    J = i13(IJ,I)
                    
                    if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J
                    
                    call dipole_T(uscale(2)-One,I,J,TTens)
                    
                enddo
            end if

            if (uscale(3).ne.one) then
                do IJ = 1, n14(I)
                    J = i14(IJ,I)
                    
                    if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J
                    
                    call dipole_T(uscale(3)-One,I,J,TTens)

                enddo
            end if

            if (uscale(4).ne.one) then
                do IJ = 1, n15(I)
                    J = i15(IJ,I)
                    
                    if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J
                    
                    call dipole_T(uscale(4)-One,I,J,TTens)
                    
                enddo
            end if
        enddo 
    end if

    ! Convert polarization tensor to the matrix
    TMat = RESHAPE(TTens, (/3*pol_atoms, 3*pol_atoms/))  
    
end subroutine create_TMat

!subroutine inverse_TMat(TMat)
!    use precision
!    use mmpol, only : pol_atoms
!    
!    implicit none
!    !
!    real(rp), dimension(3*pol_atoms,3*pol_atoms), intent(inout) :: TMat
!    !
!    integer(ip) :: LDT
!    integer(ip) :: Info
!    integer(ip),dimension(3*pol_atoms) :: Piv
!    real(rp),dimension(3*pol_atoms) :: Work
!    !
!
!    !  Compute the LU decomposition of T
!    LDT   = 3*pol_atoms
!    call DGETRF(LDT,LDT,TMat,LDT,Piv,Info)
!
!    !Compute the inverse of T
!    call DGETRI(LDT,TMat,LDT,Piv,Work,LDT,Info)
!    
!end subroutine inverse_TMat

subroutine induce_dipoles(e,ipd) !,ipd)
    use precision
    use mmpol, only : pol_atoms, n_ipd, Amoeba, pol
    !
    ! Induce dipoles on the polarizable atoms. The dipoles are induced
    ! by the electric field of the static multipoles (Later change to 
    ! the total electric field including the QM system). For the Amoeba
    ! force-field only polarization field and p-dipoles are considered.
    !
    ! - For AMBER force-field the induced dipoles will be stored in:
    !
    !       ipd(:,:,1)
    !
    ! - For Amoeba force-field the direct and polarization dipoles are
    !
    !       ipd(:,:,1)      ! d dipoles induced by direct field
    !       ipd(:,:,2)      ! p dipoles induced by polarization field 
    !                         and other p dipoles.
    !
    
    implicit none
    !
    real(rp), dimension(3,pol_atoms,n_ipd), intent(in) :: e
    real(rp), dimension(3,pol_atoms,n_ipd), intent(out) :: ipd
    !
    integer(ip) :: LDT,I
    integer(ip) :: Info
    real(rp), dimension(3*pol_atoms,3*pol_atoms) :: TMat
    real(rp),dimension(3*pol_atoms) :: d_vec
    real(rp),dimension(3*pol_atoms) :: e_vec
    integer(ip),dimension(3*pol_atoms) :: iPiv
    real(rp),dimension(3*pol_atoms) :: Work
    !
    real(rp), parameter   :: Zero = 0.0_rp, One = 1.0_rp
    !
    
    ! intermediates
    LDT = 3*pol_atoms
    
    ! Create dipole polarization matrix
    call create_TMat(TMat)
    
    
    !Compute the inverse of TMat
    call DGETRF(LDT,LDT,TMat,LDT,iPiv,Info)
    call DGETRI(LDT,TMat,LDT,iPiv,Work,LDT,Info)
    
    
    ! Reshape electric field matrix into a vector
    if (Amoeba) then
        e_vec = RESHAPE(e(:,:,2), (/3*pol_atoms/))    ! the polarization field for Amoeba FF
    else
        e_vec = RESHAPE(e(:,:,1), (/3*pol_atoms/))
    end if 
    
    
    ! Calculate dipoles with matrix inversion (for Amoeba only p dipoles)
    call DGEMM('N','N',LDT,1,LDT,One,TMat,LDT,e_vec,LDT,Zero,d_vec,LDT) !XGEMM
    
    ! For Amoeba FF compute induced dipoles by the direct field
    if (Amoeba) then
        do I=1,pol_atoms
            ipd( 1, I, 1) = pol(I)*e(1,I,1)
            ipd( 2, I, 1) = pol(I)*e(2,I,1)
            ipd( 3, I, 1) = pol(I)*e(3,I,1)
        enddo
    end if
    
    
    ! Reshape dipole vector into the matrix 
    if (Amoeba) then
        ipd( :, :, 2) = RESHAPE(d_vec, (/3,pol_atoms/)) 
    else
        ipd( :, :, 1) = RESHAPE(d_vec, (/3,pol_atoms/)) 
    end if
    
end subroutine induce_dipoles
    
end module polar


! Initialize everything by zeros