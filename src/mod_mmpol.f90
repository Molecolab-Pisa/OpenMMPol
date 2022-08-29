module mod_mmpol
    !! Main module for the control of openMMPol library. It contains
    !! all the scalar and vector (allocatable) quantities needed to
    !! build up the atomistic polarizable embedding model and perform
    !! the calculation required from the quantum chemical software.
    
    use mod_memory, only: ip, rp
    use mod_adjacency_mat, only: yale_sparse
    use mod_io, only: ommp_message

    implicit none 
    !private TODO
    
    integer(ip), protected :: ff_type
    !! Force field type selection flag (0 for AMBER, 1 for AMOEBA)

    logical, protected :: amoeba
    !! AMOEBA FF = True; WANG-AMBER = False
    
    integer(ip), protected :: mm_atoms !! number of MM atoms
    integer(ip), protected :: pol_atoms !! number of polarizable atoms
    integer(ip), protected :: ld_cart, ld_cder
!!     size of the cartesian multipolar distribution (i.e., (l+1)*(l+2)*(l+3)/6)
!!     this is 1 for AMBER (charges only), 10 for AMOEBA (up to quadrupoles). 
!!     this is also the size of the array that contains the electrostatic properties
!!     of the sources at the sources. ld_cder is the leading size of the derivative of
!!     such a distribution, which is 3 for AMBER and 19 for AMOEBA.
    integer(ip), protected :: n_ipd 
    !! number of induced point dipoles distributions 
    !! this is 1 for AMBER and 2 for AMOEBA
    
    ! arrays for the force field dependent exclusion factors. 
    
    real(rp), protected :: mscale(4)
    !! factors for charge-charge (or multipole-multipole) interactions

    real(rp), protected :: pscale(4)
    !! factors for chrage-ipd (or multipole-ipd) interactions.
    !! in AMOEBA, this is used to define the polarization field, i.e., the right-hand
    !! side to the polarization equations, and depends on the connectivity.

    real(rp), protected :: pscale_intra(4)
    !! Only used for AMOEBA, same as pscale but for atoms that belong to the 
    !! same polarization group
    
    real(rp), protected :: dscale(4)
    !! factors for multipoles-ipd interactions used to compute the direct field,
    !! which is used to define the polarization energy. these factors depend on 
    !! the polarization group "connectivity" (AMOEBA only)

    real(rp), protected :: uscale(4)
    !! factor for ipd-ipd interactions. these depend on the connectivity (AMBER)
    !! or on the polarization group " connectivity (AMOEBA)

    ! allocatable arrays which describe the polarizable system
    
    real(rp), allocatable :: thole(:)
    !! array to store the thole factors for computing damping functions
    
    real(rp) :: thole_scale
    !! Scale factor for thole damping (only used by non-AMOEBA FF); all
    !! the element of thole(:) are multiplied by thole_scale ** 0.5
    
    real(rp), allocatable, target :: cmm(:,:)
    !! Coordinates of MM atoms (3:mm_atoms)
    
    real(rp), allocatable, target :: cpol(:,:)
    !! Coordinates of polarizable atoms (3:pol_atoms)
    
    real(rp), allocatable, target :: q(:,:)
    !! Mutlipolar distribution (ld_cart:mm_atoms)
    !! For AMOEBA this is the rotated distribution.
    !! The order for the stored multipoles is
    !! q, px, py, pz, Qxx, Qxy, Qyy, Qxz, Qyx, Qzz.

    real(rp), allocatable, target :: q0(:,:)
    !! Unrotated utlipolar distribution (ld_cart:mm_atoms)
    !! (AMOEBA only)
    
    real(rp), allocatable, target :: ipd(:,:,:)
    !! induced point dipoles (3:pol_atoms:ipd) 
    logical :: ipd_done
    
    real(rp), allocatable :: pol(:)
    !! Polarizabilities for each polarizable atom
    
    integer(ip), allocatable :: mm_polar(:)
    !! indices of the MM atoms that are polarizable

    integer(ip), allocatable, target :: polar_mm(:)
    !! positions of a polarizable atom in the mm atoms list
    
    type(yale_sparse), allocatable :: conn(:)
    !! connectivity matrices listing atoms separetad by 1, 2, 3 (and 4 -- only 
    !! for AMOEBA) bonds. 1st element is the adjacency matrix.

    integer(ip), allocatable :: mmat_polgrp(:)
    !! Polarizability group index for each MM site

    type(yale_sparse) :: polgrp_mmat
    !! For each polarization group index, list all the MM atoms included.
    !! It basically is a sparse boolean matrix of dimension 
    !! N_polgroups x N_mmatoms

    type(yale_sparse), allocatable :: pg_conn(:)
    !! Adjacency and connectivity matytrices between polarizability groups.
    !! Two groups are said to be adjacent if they are connected by a chemical 
    !! bond. The 1st element is the identity matrix for code simplicity.
    
    ! parameters for the definition of the rotation matrices for the multipoles:
    integer(ip), allocatable :: mol_frame(:)
    !! definition of the molecular frame
    !! convention: 0 ... do not rotate
    !!             1 ... z-then-x
    !!             2 ... bisector
    !!             3 ... z-only
    !!             4 ... z-bisector
    !!             5 ... 3-fold

    integer(ip), allocatable :: ix(:), iy(:), iz(:)
    !! neighboring atoms used to define the axes of the molecular frame
    
    contains

    subroutine mmpol_init(l_ff_type, l_mm_atoms, l_pol_atoms)
        !! Performs all the memory allocation and vector initialization
        !! needed to run the openMMPol library
        
        use mod_memory, only: ip, rp, mallocate

        implicit none

        integer(ip), intent(in) :: l_ff_type
        !! Force field type used in initialization
        
        integer(ip), intent(in) :: l_mm_atoms
        !! Number of MM atoms used in initialization
        
        integer(ip), intent(in) :: l_pol_atoms
        !! Number of polarizable atoms used in initialization
        
        ! FF related settings
        ff_type = l_ff_type
        mm_atoms = l_mm_atoms
        pol_atoms = l_pol_atoms

        if(ff_type == 1) then
            amoeba = .true.
            ld_cart = 10_ip
            ld_cder = 19_ip
            n_ipd = 2_ip
        else if(ff_type == 0) then
            amoeba = .false.
            ld_cart = 1_ip
            ld_cder = 3_ip
            n_ipd = 1_ip
        end if
  
        ! Memory allocation
        call mallocate('mmpol_init [cmm]', 3_ip, mm_atoms, cmm)
        call mallocate('mmpol_init [q]', ld_cart, mm_atoms, q)
        call mallocate('mmpol_init [pol]', pol_atoms, pol)
        call mallocate('mmpol_init [cpol]', 3_ip, pol_atoms, cpol)
        call mallocate('mmpol_init [polar_mm]', pol_atoms, polar_mm)
        call mallocate('mmpol_init [mm_polar]', mm_atoms, mm_polar)
        call mallocate('mmpol_init [thole]', mm_atoms, thole)
        
        call mallocate('mmpol_init [idp]', 3_ip, pol_atoms, n_ipd, ipd) 
        ipd_done = .false.
        ipd = 0.0_rp

        allocate(conn(1)) 
        ! Temporary allocation, it should be allocated of the proper
        ! size when all the connectivity matricies are built, now
        ! it should only contain adjacency matrix.

        if (amoeba) then
            ! Extra quantities that should be allocated only
            ! for AMOEBA
            call mallocate('mmpol_init [q0]', ld_cart, mm_atoms, q0)
            
            call mallocate('mmpol_init [mmat_polgrp]', mm_atoms, mmat_polgrp)

            call mallocate('mmpol_init [mol_frame]', mm_atoms, mol_frame)
            call mallocate('mmpol_init [ix]', mm_atoms, ix)
            call mallocate('mmpol_init [iy]', mm_atoms, iy)
            call mallocate('mmpol_init [iz]', mm_atoms, iz)
        end if
  
    end subroutine mmpol_init

    subroutine mmpol_prepare()
        !! Compute some derived quantities from the input that 
        !! are used during the calculation. The upstream code have
        !! to provide cmm, q, pol, adjacency matrix and in
        !! the case of AMOEBA also multipoles rotation information, and 
        !! polarization group information.   
        !! This routine    
        !!   * compute connectivity lists from connected atoms    
        !!   * invert polar_mm list creating mm_polar   
        !!   * populate cpol list of coordinates   
        !!   * compute factors for thole damping    
        !!   * scales by 1/3 AMOEBA quadrupoles (?)    
        !!   * Build list for polarization groups, compute groups connectivity   
        !!   * performs multipoles rotation   

        use mod_adjacency_mat, only: build_conn_upto_n, matcpy
        use mod_io, only: ommp_message
        use mod_constants, only: OMMP_VERBOSE_DEBUG

        implicit none

        integer(ip) :: i
        real(rp) :: xx(3) ! TODO remove this variable
        
        type(yale_sparse) :: adj, pg_adj

        call ommp_message("Building connectivity lists", OMMP_VERBOSE_DEBUG)
        
        ! compute connectivity lists from connected atoms
        if(size(conn) < 4) then
            call matcpy(conn(1), adj)
            deallocate(conn)
            call build_conn_upto_n(adj, 4, conn, .false.)
        end if
        
        call ommp_message("Creating MM->polar and polar->MM lists", OMMP_VERBOSE_DEBUG)

        ! invert mm_polar list creating mm_polar
        mm_polar(:) = 0
        do i = 1, pol_atoms
            mm_polar(polar_mm(i)) = i
        end do

        call ommp_message("Populating coordinates of polarizable atoms", OMMP_VERBOSE_DEBUG)
        ! populate cpol list of coordinates
        do i = 1, pol_atoms
            cpol(:,i) = cmm(:, polar_mm(i))
        end do

        call ommp_message("Setting Thole factors", OMMP_VERBOSE_DEBUG)
        ! compute factors for thole damping
        
        if(amoeba) then
            call thole_init()
        else
            call thole_init(thole_scale)
        end if

        if(amoeba) then
            ! Copy multipoles from q to q0
            q0 = q

            ! scales by 1/3 AMOEBA quadrupoles (?)
            ! Mysterious division of multipoles by three
            ! FL told me that it was done like that in
            ! Tinker
            q0(5:10,:) = q0(5:10,:) / 3.0_rp

            ! polarization groups connectivity list
            call reverse_polgrp_tab(mmat_polgrp, polgrp_mmat)
            call build_pg_adjacency_matrix(pg_adj)
            call build_conn_upto_n(pg_adj, 3, pg_conn, .true.)

            ! performs multipoles rotation
            call rotate_multipoles(.false.,xx,xx)
        end if

    end subroutine mmpol_prepare

    subroutine mmpol_terminate()
        !! Performs all the deallocation needed at the end of the 
        !! calculation
        use mod_memory, only: mfree
        use mod_adjacency_mat, only: matfree

        implicit none 

        integer(ip) :: i

        call mfree('mmpol_terminate [cmm]', cmm)
        call mfree('mmpol_terminate [q]', q)
        call mfree('mmpol_terminate [pol]', pol)
        call mfree('mmpol_terminate [cpol]', cpol)
        call mfree('mmpol_terminate [polar_mm]', polar_mm)
        call mfree('mmpol_terminate [mm_polar]', mm_polar)
        call mfree('mmpol_terminate [thole]', thole)
        call mfree('mmpol_terminate [idp]', ipd) 
        
        do i=1, size(conn)
            call matfree(conn(i))
        end do
        deallocate(conn)

        if (amoeba) then
            ! Extra quantities that should be deallocated only
            ! for AMOEBA
            
            ! Second set of multipoles (q0 = unrotated, q=rotated)
            call mfree('mmpol_terminate [q0]', q0)
            
            ! Polarization groups
            call mfree('mmpol_terminate [mmat_polgrp]', mmat_polgrp)
            do i=1, size(pg_conn)
                call matfree(pg_conn(i))
            end do
            deallocate(pg_conn)
            call matfree(polgrp_mmat)

            ! Information for multipoles rotation
            call mfree('mmpol_terminate [mol_frame]', mol_frame)
            call mfree('mmpol_terminate [ix]', ix)
            call mfree('mmpol_terminate [iy]', iy)
            call mfree('mmpol_terminate [iz]', iz)
        end if

    end subroutine mmpol_terminate
    
    subroutine fatal_error(message)
        !! Prints a message and exit from the program. This
        !! function should be used in all the conditions 
        !! where the program cannot proceed.

        implicit none
      
        character (len=*) message
        !! Message to print before the program termination

        write(6, '(t3,a)') message
        stop '   error termination for open_mmpol.'
    end subroutine fatal_error

    subroutine thole_init(asc)
        ! This routine compute the thole factors and stores
        ! them in a vector. TODO add reference
        
        use mod_constants, only: OMMP_VERBOSE_LOW
        implicit none

        real(rp), optional, intent(in) :: asc
        
        integer(ip) :: i, j
        
        thole = 0.0_rp
        
        do i = 1, pol_atoms
            j = polar_mm(i)
            thole(j) = pol(i) ** (1.0_rp/6.0_rp)
        end do
        
        if(.not. amoeba) then
            if(present(asc)) then
                thole = thole * sqrt(asc)
            else
                call fatal_error("Scale factor for Thole damping should be passed &
                                 &to thole_init() when non-AMOEBA FF are used")
            end if
        else
            if(present(asc)) then
                call ommp_message("Scale factor passed to thole_init is &
                    &ignored because AMOEBA FF is used", OMMP_VERBOSE_LOW)
            end if
        end if
    end subroutine thole_init

    subroutine reverse_polgrp_tab(mm2pg, pg2mm)
        !! Takes as argument an array of polarization group index for each
        !! atom, and create a list of atms in each group using the boolean
        !! sparse matrix format (saved as Yale format).
        
        implicit none

        integer(ip), intent(in) :: mm2pg(mm_atoms)
        !! Index of polarization group for each MM atom
        type(yale_sparse), intent(out) :: pg2mm
        !! Indices of atoms included in each polarization group;
        !! Atom indeces for the n-th group are found at 
        !! pg2mm%ci(pg2mm%ri(n):pg2mm%ri(n+1)-1)

        integer(ip) :: i, j

        ! Allocation of Yale fmt sparse matrix
        pg2mm%n = maxval(mm2pg)
        allocate(pg2mm%ri(pg2mm%n+1))
        allocate(pg2mm%ci(mm_atoms))
        pg2mm%ri(1) = 1

        do i=1, pg2mm%n
            pg2mm%ri(i+1) = pg2mm%ri(i)
            
            do j=1, mm_atoms
                if(mm2pg(j) /= i) cycle
                
                pg2mm%ci(pg2mm%ri(i+1)) = j
                pg2mm%ri(i+1) = pg2mm%ri(i+1) + 1
            end do
        end do
    end subroutine reverse_polgrp_tab

    subroutine build_pg_adjacency_matrix(adj)
        !! Builds the adjacency matrix of polarization groups starting from
        !! atomic adjacency matrix and list of polarization groups indices.

        use mod_adjacency_mat, only: reallocate_mat 

        implicit none

        type(yale_sparse), intent(out) :: adj
        !! The group adjacency matrix to be saved.

        integer(ip) :: npg, pg1, atm1, atm2, i, j

        npg = polgrp_mmat%n

        adj%n = npg
        allocate(adj%ri(adj%n+1))
        allocate(adj%ci(adj%n*2))
        adj%ri(1) = 1

        do pg1=1, npg
            ! For each polarization group
            adj%ri(pg1+1) = adj%ri(pg1)

            do i=polgrp_mmat%ri(pg1), polgrp_mmat%ri(pg1+1)-1
                ! Loop on every atom of the group
                atm1 = polgrp_mmat%ci(i)
                do j=conn(1)%ri(atm1), conn(1)%ri(atm1+1)-1
                    ! Loop on each connected atom...
                    atm2 = conn(1)%ci(j)

                    ! If the two atoms are in different PG, then the two
                    ! polarization groups are connected. 
                    if(mmat_polgrp(atm1) /= mmat_polgrp(atm2) .and. &
                       ! if the group is not already present in the matrix
                       all(adj%ci(adj%ri(pg1):adj%ri(pg1+1)-1) /= mmat_polgrp(atm2))) then
                        adj%ci(adj%ri(pg1+1)) = mmat_polgrp(atm2)
                        adj%ri(pg1+1) = adj%ri(pg1+1) + 1
                        if(adj%ri(pg1+1) > size(adj%ci)) then
                            ! If matrix is too small, it could be enlarged...
                            call reallocate_mat(adj, size(adj%ci)+adj%n)
                        end if
                    end if
                end do
            end do
        end do
        
        ! Finally trim the output matrix
        call reallocate_mat(adj, adj%ri(adj%n+1)-1)

    end subroutine build_pg_adjacency_matrix
    
    subroutine set_screening_parameters(m, p, d, u, i)
        !! Subroutine to initialize the screening parameters
       
        implicit none

        real(rp), intent(in) :: m(4), p(4), d(4), u(4)
        real(rp), optional, intent(in) :: i(4)
        
        mscale = m
        pscale = p
        dscale = d
        uscale = u
        
        if(present(i)) then
            if(amoeba) then
                pscale_intra = i
            else
                call fatal_error("Scale factors for atoms of the same group &
                                 &cannot be set outside AMOEBA FF")
            end if
        else
            if(amoeba) &
                call fatal_error("Scale factors for atoms of the same group &
                                 &should be defined in AMOEBA FF")
        end if
        
    end subroutine set_screening_parameters
    
    subroutine mmpol_print_summary(of_name)
        !! Prints a complete summary of all the quantities stored 
        !! in the MMPol module

        use mod_io, only: iof_mmpol, print_matrix, print_int_vec
        
        implicit none

        character(len=*), intent(in), optional :: of_name
        
        integer(ip) :: of_unit

        integer(ip) :: i, j, grp, igrp, lst(1000), ilst
        real(rp), dimension(mm_atoms) :: polar ! Polarizabilities of all atoms
        character(len=120) :: str

        if(present(of_name)) then
            of_unit = 101
            open(unit=of_unit, &
                 file=of_name(1:len(trim(of_name))), &
                 action='write')
        else
            of_unit = iof_mmpol
        end if

        polar = 0.0_rp
        do i=1, pol_atoms
            polar(polar_mm(i)) = pol(i)
        end do

        write(of_unit, '(A, 4F8.4)') 'mscale: ', mscale
        write(of_unit, '(A, 4F8.4)') 'pscale: ', pscale
        if(amoeba) write(of_unit, '(A, 4F8.4)') 'pscale (intra): ', pscale_intra
        write(of_unit, '(A, 4F8.4)') 'dscale: ', dscale
        write(of_unit, '(A, 4F8.4)') 'uscale: ', uscale

        call print_matrix(.true., 'coordinates:', cmm, of_unit)
        if (amoeba) then
            call print_matrix(.true., 'multipoles - non rotated:', q0, of_unit)
        end if
        call print_matrix(.true., 'multipoles :', q, of_unit)
        call print_matrix(.true., 'coordinates of polarizable atoms:', cpol, of_unit)
        call print_matrix(.false., 'polarizabilities:', polar, of_unit)
        call print_matrix(.false., 'thole factors:',thole, of_unit)
        call print_int_vec('mm_polar list:', mm_polar, .false., of_unit)
        call print_int_vec('polar_mm list:', polar_mm, .false., of_unit)

        ! write the connectivity information for each atom:
  1000  format(t3,'connectivity information for the ',i8,'-th atom:')
    
        do i = 1, mm_atoms
            write(of_unit, 1000) i
            
            do j=1, 4
                if(j == 4 .and. .not. amoeba) cycle
                
                write(str, "('1-', I1, ' neighors:')") j+1
                call print_int_vec(trim(str), &
                                   conn(j)%ci, &
                                   .true., of_unit, &
                                   conn(j)%ri(i), &
                                   conn(j)%ri(i+1)-1) 
            end do
            
            if(amoeba) then 
                do j=1, 4
                    ilst = 1
                    do igrp=pg_conn(j)%ri(mmat_polgrp(i)), &
                            pg_conn(j)%ri(mmat_polgrp(i)+1)-1
                        grp = pg_conn(j)%ci(igrp)
                        lst(ilst:ilst+polgrp_mmat%ri(grp+1)-polgrp_mmat%ri(grp)-1) = &
                        polgrp_mmat%ci(polgrp_mmat%ri(grp):polgrp_mmat%ri(grp+1)-1)
                        ilst = ilst+polgrp_mmat%ri(grp+1)-polgrp_mmat%ri(grp)
                    end do
                   
                    write(str, "('1-', I1, ' polarization neighors:')") j
                    if(ilst == 1) ilst = 0
                    ! needed to addres the empty array case
                    call print_int_vec(trim(str), &
                                       lst, &
                                       .true., of_unit, &
                                       1, ilst-1)
                end do
            end if
        end do
        
        if(present(of_name)) close(of_unit)

    end subroutine mmpol_print_summary

end module mod_mmpol
