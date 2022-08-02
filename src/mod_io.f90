module mod_io
#ifdef USE_HDF5
    use hdf5, only: hid_t
#endif
    
    implicit none
    private

    integer, parameter :: iof_memory = 6
    integer, parameter :: iof_mmpol = 6
    integer, parameter :: iof_mmpinp = 100
    public :: iof_memory, iof_mmpol, iof_mmpinp
    public :: mmpol_print_summary, print_header
    public :: print_matrix, print_int_vec

#ifdef USE_HDF5
#define H5T_RP H5T_NATIVE_DOUBLE
#define H5T_LOGICAL H5T_BITFIELD_F
#ifdef USE_I8
#define H5T_IP H5T_STD_I64LE
#else
#define H5T_IP H5T_STD_I32LE
#endif
    integer(hid_t) :: iof_hdf5_out = 101
    public :: mmpol_save_as_hdf5
    
    interface hdf5_add_scalar
        ! Write a scalar as an attribute of the group
        module procedure r_hdf5_add_scalar
        module procedure i_hdf5_add_scalar
        module procedure l_hdf5_add_scalar
    end interface hdf5_add_scalar
    
    interface hdf5_add_array
        ! Write a scalar as an attribute of the group
        module procedure r1_hdf5_add_array
        module procedure r2_hdf5_add_array
        module procedure r3_hdf5_add_array

        module procedure i1_hdf5_add_array
        module procedure i2_hdf5_add_array
        module procedure i3_hdf5_add_array
    end interface hdf5_add_array
#endif
    
    contains
    
#ifdef USE_HDF5
    ! Subroutines dedicated to HDF5 I/O
    !subroutine hdf5_add_array(hid, label, array)

    !end subroutine

    subroutine r_hdf5_add_scalar(hid, label, scalar)
        use hdf5
        use mod_memory, only: ip, rp
        
        implicit none
        
        integer(hid_t), intent(in) :: hid
        character(len=*), intent(in) :: label
        real(rp), intent(in) :: scalar

        integer(hsize_t), dimension(1), parameter :: dims = [1]
        integer(hid_t) :: cur_dst, cur_dsp
        integer(kind=4) :: eflag
        
        call H5Screate_f(H5S_SCALAR_F, cur_dsp, eflag)
        
        call H5Acreate_f(hid, &
                         label, &
                         H5T_RP, &
                         cur_dsp, cur_dst, eflag)
        call H5Awrite_f(cur_dst, H5T_RP, scalar, dims, eflag)
    end subroutine
    
    subroutine i_hdf5_add_scalar(hid, label, scalar)
        use hdf5
        use mod_memory, only: ip
        
        implicit none
        
        integer(hid_t), intent(in) :: hid
        character(len=*), intent(in) :: label
        integer(ip), intent(in) :: scalar

        integer(hsize_t), dimension(1), parameter :: dims = [1]
        integer(hid_t) :: cur_dst, cur_dsp
        integer(kind=4) :: eflag
        
        call H5Screate_f(H5S_SCALAR_F, cur_dsp, eflag)
        call H5Acreate_f(hid, &
                         label, &
                         H5T_IP, &
                         cur_dsp, cur_dst, eflag)
        call H5Awrite_f(cur_dst, H5T_IP, scalar, dims, eflag)
    end subroutine
    
    subroutine l_hdf5_add_scalar(hid, label, scalar)
        use hdf5
        use mod_memory, only: ip
        
        implicit none
        
        integer(hid_t), intent(in) :: hid
        character(len=*), intent(in) :: label
        logical, intent(in) :: scalar

        integer(hsize_t), dimension(1), parameter :: dims = [1]
        integer(hid_t) :: cur_dst, cur_dsp
        integer(kind=4) :: eflag
        
        call H5Screate_f(H5S_SCALAR_F, cur_dsp, eflag)
        call H5Acreate_f(hid, &
                         label, &
                         H5T_IP, &
                         cur_dsp, cur_dst, eflag)
        if(scalar) then
            call H5Awrite_f(cur_dst, H5T_IP, 1, dims, eflag)
        else
            call H5Awrite_f(cur_dst, H5T_IP, 0, dims, eflag)
        end if

    end subroutine
    
    subroutine r1_hdf5_add_array(hid, label, v)
        use hdf5
        use mod_memory, only: ip, rp
        
        implicit none
        
        integer(hid_t), intent(in) :: hid
        character(len=*), intent(in) :: label
        real(rp), intent(in), dimension(:) :: v

        integer(hsize_t), dimension(1) :: dims
        integer(hid_t) :: cur_dst, cur_dsp
        integer(kind=4) :: eflag
        
        dims = shape(v)
        call h5screate_simple_f(1, dims, cur_dsp, eflag)
        call h5dcreate_f(hid, &
                         label, &
                         H5T_RP, &
                         cur_dsp, cur_dst, eflag)
        call h5dwrite_f(cur_dst, H5T_RP, v, dims, eflag)
    end subroutine
    
    subroutine r2_hdf5_add_array(hid, label, v)
        use hdf5
        use mod_memory, only: ip, rp
        
        implicit none
        
        integer(hid_t), intent(in) :: hid
        character(len=*), intent(in) :: label
        real(rp), intent(in), dimension(:,:) :: v

        integer(hsize_t), dimension(2) :: dims
        integer(hid_t) :: cur_dst, cur_dsp
        integer(kind=4) :: eflag
        
        dims = shape(v)
        call h5screate_simple_f(2, dims, cur_dsp, eflag)
        call h5dcreate_f(hid, &
                         label, &
                         H5T_RP, &
                         cur_dsp, cur_dst, eflag)
        call h5dwrite_f(cur_dst, H5T_RP, v, dims, eflag)
    end subroutine
    
    subroutine r3_hdf5_add_array(hid, label, v)
        use hdf5
        use mod_memory, only: ip, rp
        
        implicit none
        
        integer(hid_t), intent(in) :: hid
        character(len=*), intent(in) :: label
        real(rp), intent(in), dimension(:,:,:) :: v

        integer(hsize_t), dimension(3) :: dims
        integer(hid_t) :: cur_dst, cur_dsp
        integer(kind=4) :: eflag
        
        dims = shape(v)
        call h5screate_simple_f(3, dims, cur_dsp, eflag)
        call h5dcreate_f(hid, &
                         label, &
                         H5T_RP, &
                         cur_dsp, cur_dst, eflag)
        call h5dwrite_f(cur_dst, H5T_RP, v, dims, eflag)
    end subroutine
    
    subroutine i1_hdf5_add_array(hid, label, v)
        use hdf5
        use mod_memory, only: ip, rp
        
        implicit none
        
        integer(hid_t), intent(in) :: hid
        character(len=*), intent(in) :: label
        integer(ip), intent(in), dimension(:) :: v

        integer(hsize_t), dimension(1) :: dims
        integer(hid_t) :: cur_dst, cur_dsp
        integer(kind=4) :: eflag
        
        dims = shape(v)
        call h5screate_simple_f(1, dims, cur_dsp, eflag)
        call h5dcreate_f(hid, &
                         label, &
                         H5T_IP, &
                         cur_dsp, cur_dst, eflag)
        call h5dwrite_f(cur_dst, H5T_IP, v, dims, eflag)
    end subroutine
    
    subroutine i2_hdf5_add_array(hid, label, v)
        use hdf5
        use mod_memory, only: ip, rp
        
        implicit none
        
        integer(hid_t), intent(in) :: hid
        character(len=*), intent(in) :: label
        integer(ip), intent(in), dimension(:,:) :: v

        integer(hsize_t), dimension(2) :: dims
        integer(hid_t) :: cur_dst, cur_dsp
        integer(kind=4) :: eflag
        
        dims = shape(v)
        call h5screate_simple_f(2, dims, cur_dsp, eflag)
        call h5dcreate_f(hid, &
                         label, &
                         H5T_IP, &
                         cur_dsp, cur_dst, eflag)
        call h5dwrite_f(cur_dst, H5T_IP, v, dims, eflag)
    end subroutine
    
    subroutine i3_hdf5_add_array(hid, label, v)
        use hdf5
        use mod_memory, only: ip, rp
        
        implicit none
        
        integer(hid_t), intent(in) :: hid
        character(len=*), intent(in) :: label
        integer(ip), intent(in), dimension(:,:,:) :: v

        integer(hsize_t), dimension(3) :: dims 
        integer(hid_t) :: cur_dst, cur_dsp
        integer(kind=4) :: eflag
       
        dims = shape(v)
        call h5screate_simple_f(3, dims, cur_dsp, eflag)
        call h5dcreate_f(hid, &
                         label, &
                         H5T_IP, &
                         cur_dsp, cur_dst, eflag)
        call h5dwrite_f(cur_dst, H5T_IP, v, dims, eflag)
    end subroutine

    subroutine mmpol_save_as_hdf5(filename, out_fail)
        use hdf5
        use mod_memory, only: ip
        use mod_mmpol, only: mm_atoms, pol_atoms, cmm, polar_mm, ld_cart, q, &
                             q0, amoeba, pol, conn, ff_type, &
                             ix, iy, iz, mol_frame, mmat_polgrp
        use mod_bonded

        implicit none

        character(len=*), intent(in) :: filename
        integer(ip), intent(out) :: out_fail
        
        integer(hid_t) :: hg_sysmodel, hg_res, hg_cur, hg_amoeba, &
                          hg_top, hg_cur_param, hg_cur_bp, cur_dst, cur_dsp
        integer(hsize_t), dimension(4) :: dims
        integer(kind=4) :: eflag

        ! Initialize interface
        call h5open_f(eflag)
        if(eflag /= 0) then
            write(iof_mmpol, *) "Unable to initialize HDF5 module. Failure in &
                               &h5open_f subroutine."
            out_fail = -1_ip
            return
        end if
        
        call h5fcreate_f(filename, H5F_ACC_EXCL_F, iof_hdf5_out, eflag)
        if( eflag /= 0) then 
            write(iof_mmpol, *) "Unable to create HDF5 file. Failure in &
                               &h5fcreate_f subroutine."
            out_fail = -1_ip
            return
        end if
        
        call h5gcreate_f(iof_hdf5_out, "system model", hg_sysmodel, eflag)
        if( eflag /= 0) then 
            write(iof_mmpol, *) "Error while creating group 'system model.'&
                                &Failure in h5gcreate_f subroutine."
            out_fail = -1_ip
            return
        end if
        
        ! Write H5G:system model
        ! Those are the quantities needed to describe the model potential of
        ! the system.
        ! It should not contain any information on/derived from atoms positions

        call hdf5_add_scalar(hg_sysmodel, "N-atoms", mm_atoms)
        call hdf5_add_scalar(hg_sysmodel, "N-pol-atoms", pol_atoms)
        call hdf5_add_scalar(hg_sysmodel, "FF type", ff_type)
        
        call h5gcreate_f(hg_sysmodel, "topology", hg_top, eflag)
        
        call h5gcreate_f(hg_top, "connectivity", hg_cur, eflag)
        call hdf5_add_array(hg_cur, "ADJ1-RowIdx", conn(1)%ri) 
        call hdf5_add_array(hg_cur, "ADJ1-ColIdx", conn(1)%ci) 
        call h5gclose_f(hg_cur, eflag)
        
        if(amoeba) then
            call h5gcreate_f(hg_top, "AMOEBA", hg_amoeba, eflag)
            call hdf5_add_array(hg_amoeba, "polarization_group_id", mmat_polgrp)
            call h5gclose_f(hg_amoeba, eflag)
        end if
        
        call h5gclose_f(hg_top, eflag)

        call h5gcreate_f(hg_sysmodel, "parameters", hg_cur, eflag)
        
        call h5gcreate_f(hg_cur, "bonded", hg_cur_param, eflag)
        
        ! Bond stretching
        call h5gcreate_f(hg_cur_param, "stretching", hg_cur_bp, eflag)
        call hdf5_add_scalar(hg_cur_bp, "enabled", use_bond)
        if(use_bond) then
            call hdf5_add_scalar(hg_cur_bp, "cubic", bond_cubic)
            call hdf5_add_scalar(hg_cur_bp, "quartic", bond_quartic)
            call hdf5_add_array(hg_cur_bp, "k", kbond)
            call hdf5_add_array(hg_cur_bp, "l0", l0bond)
            call hdf5_add_array(hg_cur_bp, "atoms", bondat)
        end if
        call h5gclose_f(hg_cur_bp, eflag)
        
        ! Angle bending
        call h5gcreate_f(hg_cur_param, "bending", hg_cur_bp, eflag)
        call hdf5_add_scalar(hg_cur_bp, "enabled", use_angle)
        if(use_angle) then
            call hdf5_add_scalar(hg_cur_bp, "cubic", angle_cubic)
            call hdf5_add_scalar(hg_cur_bp, "quartic", angle_quartic)
            call hdf5_add_scalar(hg_cur_bp, "pentic", angle_pentic)
            call hdf5_add_scalar(hg_cur_bp, "sextic", angle_sextic)
            call hdf5_add_array(hg_cur_bp, "k", kangle)
            call hdf5_add_array(hg_cur_bp, "ang0", eqangle)
            call hdf5_add_array(hg_cur_bp, "atoms", angleat)
            call hdf5_add_array(hg_cur_bp, "type", anglety)
        end if
        call h5gclose_f(hg_cur_bp, eflag)
        
        ! Dihedral torsion
        call h5gcreate_f(hg_cur_param, "torsion", hg_cur_bp, eflag)
        call hdf5_add_scalar(hg_cur_bp, "enabled", use_torsion)
        if(use_torsion) then
            call hdf5_add_array(hg_cur_bp, "amplitudes", torsamp)
            call hdf5_add_array(hg_cur_bp, "phase", torsphase)
            call hdf5_add_array(hg_cur_bp, "atoms", torsionat)
            call hdf5_add_array(hg_cur_bp, "period", torsn)
        end if
        call h5gclose_f(hg_cur_bp, eflag)
        
        ! Stretching-bending coupling
        call h5gcreate_f(hg_cur_param, "stretching-bending", hg_cur_bp, eflag)
        call hdf5_add_scalar(hg_cur_bp, "enabled", use_strbnd)
        if(use_strbnd) then
            call hdf5_add_array(hg_cur_bp, "k1", strbndk1)
            call hdf5_add_array(hg_cur_bp, "k2", strbndk2)
            call hdf5_add_array(hg_cur_bp, "l1_0", strbndl10)
            call hdf5_add_array(hg_cur_bp, "l2_0", strbndl20)
            call hdf5_add_array(hg_cur_bp, "ang0", strbndthet0)
            call hdf5_add_array(hg_cur_bp, "atoms", strbndat)
        end if
        call h5gclose_f(hg_cur_bp, eflag)
        
        ! Stretching-torsion coupling
        call h5gcreate_f(hg_cur_param, "stretching-torsion", hg_cur_bp, eflag)
        call hdf5_add_scalar(hg_cur_bp, "enabled", use_strtor)
        if(use_strtor) then
            call hdf5_add_array(hg_cur_bp, "k", strtork)
            call hdf5_add_array(hg_cur_bp, "bonds_idx", strtor_b)
            call hdf5_add_array(hg_cur_bp, "torsion_idx", strtor_t)
            call hdf5_add_array(hg_cur_bp, "atoms", strtorat)
        end if
        call h5gclose_f(hg_cur_bp, eflag)

        ! Bending-torsion coupling 
        call h5gcreate_f(hg_cur_param, "bending-torsion", hg_cur_bp, eflag)
        call h5gclose_f(hg_cur_bp, eflag)
        ! Torsion-torsion coupling
        call h5gcreate_f(hg_cur_param, "torsion-torsion", hg_cur_bp, eflag)
        call h5gclose_f(hg_cur_bp, eflag)
        ! Pi-torsion
        call h5gcreate_f(hg_cur_param, "pi-torsion", hg_cur_bp, eflag)
        call h5gclose_f(hg_cur_bp, eflag)
        ! Out-of-plane bending
        call h5gcreate_f(hg_cur_param, "out-of-plane-bending", hg_cur_bp, eflag)
        call h5gclose_f(hg_cur_bp, eflag)
        ! Urey-Bradley stretching
        call h5gcreate_f(hg_cur_param, "urey-bradley", hg_cur_bp, eflag)
        call h5gclose_f(hg_cur_bp, eflag)

        call h5gclose_f(hg_cur_param, eflag)
        
        call h5gcreate_f(hg_cur, "non-bonded", hg_cur_param, eflag)
        call h5gclose_f(hg_cur_param, eflag)
        
        call h5gcreate_f(hg_cur, "electrostatics", hg_cur_param, eflag)
        if(amoeba) then
            ! Write the unrotated multipoles, that are coordinates independent
            call hdf5_add_array(hg_cur_param, "fixed_multipoles", q0)
            ! Write all the information needed to perform the rotation of the 
            ! multipoles
            call h5gcreate_f(hg_cur_param, "AMOEBA", hg_amoeba, eflag)
            call hdf5_add_array(hg_amoeba, "fixed_mmpoles_rot_Z", iz)
            call hdf5_add_array(hg_amoeba, "fixed_mmpoles_rot_X", ix)
            call hdf5_add_array(hg_amoeba, "fixed_mmpoles_rot_Y", iy)
            call hdf5_add_array(hg_amoeba, "fixed_mmpoles_rot_CONV", mol_frame)
            call h5gclose_f(hg_amoeba, eflag)
        else
            call hdf5_add_array(hg_cur_param, "fixed_multipoles", q)
        end if
        call hdf5_add_array(hg_cur_param, "polarizable_atoms_idx", polar_mm)
        call hdf5_add_array(hg_cur_param, "polarizabilities", pol) 

        call h5gclose_f(hg_cur_param, eflag)
        call h5gclose_f(hg_cur, eflag)
        
        call h5gclose_f(hg_sysmodel, eflag)
        if( eflag /= 0) then 
            write(iof_mmpol, *) "Error while closing group 'system fundamental.'&
                                &Failure in h5gclose_f subroutine."
            out_fail = -1_ip
            return
        end if

        
        !call h5gcreate_f(iof_hdf5_out, "system derived", hg_sysder, eflag)
        !if( eflag /= 0) then 
        !    write(iof_mmpol, *) "Error while creating group 'system derived'.&
        !                        &Failure in h5gcreate_f subroutine."
        !    out_fail = -1_ip
        !    return
        !end if
        
        call h5gcreate_f(iof_hdf5_out, "results", hg_res, eflag)
        if( eflag /= 0) then 
            write(iof_mmpol, *) "Error while creating group 'results'.&
                                &Failure in h5gcreate_f subroutine."
            out_fail = -1_ip
            return
        end if

        call h5fclose_f(iof_hdf5_out, eflag)
        if( eflag /= 0) then 
            write(iof_mmpol, *) "Error while closing HDF5 file. Failure in &
                               &h5fclose_f subroutine."
            out_fail = -1_ip
            return
        end if

        out_fail = 0_ip
    
    end subroutine mmpol_save_as_hdf5
    
#endif
    
    subroutine mmpol_print_summary(of_name)
        !! Prints a complete summary of all the quantities stored 
        !! in the MMPol module

        use mod_mmpol ! TODO only 

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

        call print_matrix(.true.,'coordinates:', &
                          3,mm_atoms,3,mm_atoms,cmm,of_unit)
        if (amoeba) then
            call print_matrix(.true.,'multipoles - non rotated:', &
                              ld_cart,mm_atoms,ld_cart,mm_atoms,q0,of_unit)
        end if
        call print_matrix(.true.,'multipoles :', &
                          ld_cart,mm_atoms,ld_cart,mm_atoms,q,of_unit)
        call print_matrix(.true.,'coordinates of polarizable atoms:', &
                          3,pol_atoms,3,pol_atoms,cpol,of_unit)
        call print_matrix(.false.,'polarizabilities:', &
                          mm_atoms,1,mm_atoms,1,polar,of_unit)
        call print_matrix(.false.,'thole factors:', &
                          mm_atoms,1,mm_atoms,1,thole,of_unit)
        call print_int_vec('mm_polar list:', &
                           mm_atoms,0,0,mm_polar,.false., of_unit)
        call print_int_vec('polar_mm list:', &
                           pol_atoms,0,0,polar_mm, .false.,of_unit)

        ! write the connectivity information for each atom:
  1000  format(t3,'connectivity information for the ',i8,'-th atom:')
    
        do i = 1, mm_atoms
            write(of_unit, 1000) i
            
            do j=1, 4
                if(j == 4 .and. .not. amoeba) cycle
                
                write(str, "('1-', I1, ' neighors:')") j+1
                call print_int_vec(trim(str), &
                                   size(conn(j)%ci), &
                                   conn(j)%ri(i), &
                                   conn(j)%ri(i+1)-1, & 
                                   conn(j)%ci, &
                                   .true.,of_unit)
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
                    call print_int_vec(trim(str), & 
                                       ilst-1,0,0,lst,.true., of_unit)
                end do
            end if
        end do
        
        if(present(of_name)) close(of_unit)

    end subroutine mmpol_print_summary

    subroutine print_header
      implicit none
      
      9000 format(t3,' .d88888b.                             888b     d888 888b     d888 8888888b.          888 ',/,&
                  t3,'d88P" "Y88b                            8888b   d8888 8888b   d8888 888   Y88b         888 ',/,&
                  t3,'888     888                            88888b.d88888 88888b.d88888 888    888         888 ',/,&
                  t3,'888     888 88888b.   .d88b.  88888b.  888Y88888P888 888Y88888P888 888   d88P .d88b.  888 ',/,&
                  t3,'888     888 888 "88b d8P  Y8b 888 "88b 888 Y888P 888 888 Y888P 888 8888888P" d88""88b 888 ',/,&
                  t3,'888     888 888  888 88888888 888  888 888  Y8P  888 888  Y8P  888 888       888  888 888 ',/,&
                  t3,'Y88b. .d88P 888 d88P Y8b.     888  888 888   "   888 888   "   888 888       Y88..88P 888 ',/,&
                  t3,' "Y88888P"  88888P"   "Y8888  888  888 888       888 888       888 888        "Y88P"  888 ',/,&
                  t3,'            888                                                                           ',/,&
                  t3,'            888                                                                           ',/,&
                  t3,'            888                                                                           ')
      9100 format(t3,'an open-source implementation of MMPol and AMOEBA embedding for polarizable QM/MM',/, &
                  t5,'by Vladislav Slama, Lorenzo Cupellini, Benedetta Mennucci, ..., Filippo Lipparini',/, &
                  t5,'MoLECoLab Pisa')

      write(iof_mmpol,9000)
      write(iof_mmpol,9100)
      write(iof_mmpol,*)
      return
    end subroutine print_header

    subroutine print_matrix(trans,label,lda,ldb,n,m,matrix, ofunit)
      !use mmpol
      use mod_memory, only: ip, rp
      implicit none
      logical,                         intent(in) :: trans
      character    (len=*),            intent(in) :: label
      integer(ip),                     intent(in) :: lda, ldb, n, m
      real(rp),    dimension(lda,ldb), intent(in) :: matrix
      integer(ip), intent(in), optional :: ofunit
    !
      integer(ip)        :: i, j, nbatch, nres, icol(5), out_unit
      character (len=24) :: iform, rform

      if(present(ofunit)) then
          out_unit = ofunit
      else
          out_unit = iof_mmpol
      end if
    !
      1000 format(t3,a)
      1010 format(t3,5i16)
      1020 format(t3,5f16.8)
      write(out_unit,1000) label
      if (trans) then
    !
        nbatch = n/5
        nres   = n - 5*nbatch
        write(iform,'("(t3,",i1,"i16)")') nres
        write(rform,'("(t3,",i1,"f16.8)")') nres
    !
        do i = 1, nbatch
          icol(1) = (i-1)*5 + 1
          icol(2) = (i-1)*5 + 2
          icol(3) = (i-1)*5 + 3
          icol(4) = (i-1)*5 + 4
          icol(5) = (i-1)*5 + 5
          write(out_unit,1010) icol(1:5)
          do j = 1, m
            write(out_unit,1020) matrix(icol(1):icol(5),j)
          end do
          write(out_unit,'("")')
        end do
    !
        if (nres.ne.0) then
          do i = 1, nres
            icol(i) = 5*nbatch + i
          end do
          write(out_unit,iform) icol(1:nres)
          do j = 1, m
            write(out_unit,rform) matrix(icol(1):icol(nres),j)
          end do
        end if
    !
      else
    !
        nbatch = m/5
        nres   = m - 5*nbatch
        write(iform,'("(t3,",i1,"i16)")') nres
        write(rform,'("(t3,",i1,"f16.8)")') nres
    !
        do i = 1, nbatch
          icol(1) = (i-1)*5 + 1
          icol(2) = (i-1)*5 + 2
          icol(3) = (i-1)*5 + 3
          icol(4) = (i-1)*5 + 4
          icol(5) = (i-1)*5 + 5
          write(out_unit,1010) icol(1:5)
          do j = 1, n
            write(out_unit,1020) matrix(j,icol(1):icol(5))
          end do
          write(out_unit,*)
        end do
    !
        if (nres.ne.0) then
          do i = 1, nres
            icol(i) = 5*nbatch + i
          end do
          write(out_unit,iform) icol(1:nres)
          do j = 1, n
            write(out_unit,rform) matrix(j,icol(1):icol(nres))
           end do
        end if
      end if
    end subroutine print_matrix
    !
    subroutine print_int_vec(label, n, ibeg, iend, vec, dosort, ofunit)
        use mod_memory, only: ip
        use mod_utils, only: sort_ivec
        
        implicit none
        
        character(len=*), intent(in) :: label
        integer(ip), intent(in) :: n, ibeg, iend
        logical, optional :: dosort
        integer(ip), intent(in), optional :: ofunit
        logical :: sort
        integer(ip), dimension(n), intent(in) :: vec
        
        integer(ip) :: ib, ie, out_unit

        integer(ip), allocatable, dimension(:) :: sorted_vec

        ib = ibeg
        ie = iend
        if(ib == 0) ib = 1
        if(ie == 0) ie = n
        if(.not. present(dosort)) then
            sort = .false.
        else
            sort =  dosort
        end if
        
        if(present(ofunit)) then
            out_unit = ofunit
        else
            out_unit = iof_mmpol
        end if
        
        write(out_unit, '(t3, a)') label
        
        if(ib > ie) then
            write(out_unit,'(t5)')
            return
        end if 

        if(sort) then
            call sort_ivec(vec(ib:ie), sorted_vec)
            write(out_unit,'(t5, 10i8)') sorted_vec
        else
            write(out_unit,'(t5, 10i8)') vec(ib:ie)
        end if

        return

    end subroutine print_int_vec

end module mod_io
