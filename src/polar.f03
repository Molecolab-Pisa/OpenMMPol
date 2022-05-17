module polar
    use mod_memory, only: ip, rp
    
    real(rp), allocatable :: TMat(:,:)
    real(rp), allocatable :: TMatI(:,:)

    contains

    subroutine dipole_T(scalef,I,J,TTens)
        use mod_mmpol, only : thole, cpol, pol_atoms, pol, polar_mm
        implicit none
        !                      
        ! Compute element of the polarization tensor TTens between
        ! polarizable cpol atom I and polarizable cpol atom J. On the
        ! TTens diagonal (I=J) are inverse polarizabilities and on the 
        ! off-diagonal dipole field.
        !
        ! Polarizabilities pol are defined for polarizable atoms only while 
        ! Thole factors are defined for all of them
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

            call coulomb_kernel(.true.,2,dx,dy,dz,thole(polar_mm(I)),thole(polar_mm(J)),rm1,rm3,rm5,rm7,rm9,rm11)

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
        use mod_mmpol, only : pol_atoms, cpol, uscale, np11, np12, np13, np14, &
                               ip11, ip12, ip13, ip14, n12, n13, n14, n15, &
                               i12, i13, i14, i15, pol, mm_polar, Amoeba,  &
                               verbose, polar_mm
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
        ! Polarizabilities pol are defined for polarizable atoms only while 
        ! n12, i12... factors are defined for all of them (even mm atoms)
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
                    do IJ = 1, np11( polar_mm(I) )
                        J = ip11(IJ, polar_mm(I) )

                        if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J

                        call dipole_T(uscale(1)-One,I, mm_polar(J) ,TTens)

                    enddo
                end if 

                if (uscale(2).ne.one) then
                    do IJ = 1, np12( polar_mm(I) )
                        J = ip12(IJ, polar_mm(I) )
                        
                        if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J
                        
                        call dipole_T(uscale(2)-One,I, mm_polar(J) ,TTens)
                    enddo
                end if 

                if (uscale(3).ne.one) then
                    do IJ = 1, np13( polar_mm(I) )
                        J = ip13(IJ, polar_mm(I) )
                        
                        if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J
                        
                        call dipole_T(uscale(3)-One,I, mm_polar(J) ,TTens)
                        
                    enddo
                end if 

                if (uscale(4).ne.one) then
                    do IJ = 1, np14( polar_mm(I) )
                        J = ip14(IJ, polar_mm(I) )
                        
                        if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J
                        
                        call dipole_T(uscale(4)-One,I, mm_polar(J) ,TTens)

                    enddo
                end if 
            enddo
            
        else     ! AMBER and sources are induced dipoles
            do I = 1,pol_atoms

                ! Field from dipoles is scaled by uscale
                if (uscale(1).ne.one) then
                    do IJ = 1, n12( polar_mm(I) )   ! polar_mm changes pol index into mm index
                        J = i12(IJ, polar_mm(I) )
                        
                        if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J
                                                    ! mm_polar changes mm index index into pol index
                        
                        call dipole_T(uscale(1)-One,I,mm_polar(J),TTens)

                    enddo
                end if

                if (uscale(2).ne.one) then
                    do IJ = 1, n13( polar_mm(I) )   ! polar_mm changes pol index into mm index
                        J = i13(IJ, polar_mm(I) )
                        
                        if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J
                                                    ! mm_polar changes mm index index into pol index
                        
                        call dipole_T(uscale(2)-One,I,mm_polar(J),TTens)
                        
                    enddo
                end if

                if (uscale(3).ne.one) then
                    do IJ = 1, n14( polar_mm(I) )   ! polar_mm changes pol index into mm index
                        J = i14(IJ, polar_mm(I) )
                        
                        if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J
                                                    ! mm_polar changes mm index index into pol index
                        
                        call dipole_T(uscale(3)-One,I,mm_polar(J),TTens)

                    enddo
                end if

                if (uscale(4).ne.one) then
                    do IJ = 1, n15( polar_mm(I) )   ! polar_mm changes pol index into mm index
                        J = i15(IJ, polar_mm(I) )
                        
                        if (mm_polar(J).eq.0) cycle ! if J is not polarizable atom, skip J
                                                    ! mm_polar changes mm index index into pol index

                        call dipole_T(uscale(4)-One,I,mm_polar(J),TTens)
                        
                    enddo
                end if
            enddo 
        end if

        ! Convert polarization tensor to the matrix
        TMat = RESHAPE(TTens, (/3*pol_atoms, 3*pol_atoms/))  
        
        ! Print the matrix if verbose output is requested
        if (verbose.ge.1) then
            call print_matrix(.true.,'K:',3*pol_atoms,3*pol_atoms,3*pol_atoms,3*pol_atoms,TMat)
        end if
        
    end subroutine create_TMat

    subroutine TMatVec(n,x,y)
        use mod_constants, only : One, Zero
        implicit none
        ! Perform matrix vector multiplication y = TMat*x,
        ! where TMat is polarization matrix and x,y are 
        ! vectors
        !
        !   variables:
        !  
        !     n        : integer, input, size of the matrix
        ! 
        !     x        : real, vector which is multiplied by the matrix
        !
        !     y        : real, solution vector
        !
        integer(ip), intent(in) :: n
        real(rp),dimension(n),intent(in)   :: x
        real(rp),dimension(n), intent(out) :: y
        !
        ! initialize solution vector
        !
        y = 0
        !
        ! Compute the matrix vector product
        !
        call DGEMM('N','N',n,1,n,One,TMat,n,x,n,Zero,y,n) !XGEMM
        !
    end subroutine TMatVec
        
    subroutine TMatVec_offdiag(n,x,y)
        use mod_constants, only : One, Zero
        implicit none
        ! Perform matrix vector multiplication y = [TMat-diag(TMat)]*x,
        ! where TMat is polarization matrix and x,y are 
        ! vectors
        !
        !   variables:
        !  
        !     n        : integer, input, size of the matrix
        ! 
        !     x        : real, vector which is multiplied by the matrix
        !
        !     y        : real, solution vector
        !
        integer(ip), intent(in) :: n
        real(rp),dimension(n),intent(in)   :: x
        real(rp),dimension(n), intent(out) :: y
        !
        real(rp),dimension(n,n) :: O
        integer(ip) :: I
        !
        ! initialize solution vector and matrix
        !
        y = 0
        O = TMat
        !
        ! Subtract the diagonal from the O matrix
        !
        do I = 1, n
            O(I,I) = Zero
        enddo
        !
        ! Compute the matrix vector product
        !
        call DGEMM('N','N',n,1,n,One,O,n,x,n,Zero,y,n) !XGEMM
        !
    end subroutine TMatVec_offdiag

    subroutine PolVec(n,x,y)
        use mod_mmpol, only : pol
        implicit none
        ! Perform matrix vector multiplication y = pol*x,
        ! where pol polarizability vector and x,y are 
        ! vectors
        !
        !   variables:
        !  
        !     n        : integer, input, size of the matrix
        ! 
        !     x        : real, vector which is multiplied by the 
        !                inverse diagonal matrix
        !
        !     y        : real, solution vector
        !
        integer(ip), intent(in) :: n
        real(rp),dimension(n),intent(in)   :: x
        real(rp),dimension(n), intent(out) :: y
        !
        integer(ip) :: I,indx
        !
        ! initialize solution vector
        !
        y = 0
        !
        ! Compute the matrix vector product
        !
        do I = 1, n
            indx = (I+2)/3
            y(I) = pol(indx)*x(I)   
        enddo
        !
    end subroutine PolVec

    subroutine induce_dipoles_MatInv(elf,dip) !,ipd)
        use mod_mmpol, only : pol_atoms, n_ipd, Amoeba, pol
        use mod_memory, only: mallocate
        !
        ! Induce dipoles on the polarizable atoms. The dipoles are induced
        ! by the electric field of the static multipoles (Later change to 
        ! the total electric field including the QM system). For the Amoeba
        ! force-field only polarization field and p-dipoles are considered.
        !
        
        implicit none
        !
        real(rp), dimension(3*pol_atoms), intent(in) :: elf
        real(rp), dimension(3*pol_atoms), intent(out) :: dip
        !
        integer(ip) :: LDT,I
        integer(ip) :: Info
        !real(rp), dimension(3*pol_atoms,3*pol_atoms) :: TMat
        integer(ip),dimension(3*pol_atoms) :: iPiv
        real(rp),dimension(3*pol_atoms) :: Work
        !
        real(rp), parameter   :: Zero = 0.0_rp, One = 1.0_rp
        !
        
        ! intermediates
        LDT = 3*pol_atoms
        
        ! Create dipole polarization matrix
        call create_TMat(TMat)
        
        if (.not.ALLOCATED(TMatI)) then
            call mallocate('polar [TMatI]',3*pol_atoms,3*pol_atoms,TMatI)
        end if
        
        ! Initialize inverse polarization matrix
        TMatI = TMat
        
        !Compute the inverse of TMat
        call DGETRF(LDT,LDT,TMatI,LDT,iPiv,Info)
        call DGETRI(LDT,TMatI,LDT,iPiv,Work,LDT,Info)
        
        
        ! Calculate dipoles with matrix inversion (for Amoeba only p dipoles)
        call DGEMM('N','N',LDT,1,LDT,One,TMatI,LDT,elf,LDT,Zero,dip,LDT) !XGEMM
        !call TMatVec(LDT,elf,dip)
       
    end subroutine induce_dipoles_MatInv

end module polar
