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
    integer(hid_t) :: iof_hdf5_out = 301, iof_hdf5_io = 302
    public :: mmpol_save_as_hdf5, mmpol_init_from_hdf5
    
    interface hdf5_add_scalar
        ! Write a scalar as an attribute of the group
        module procedure r_hdf5_add_scalar
        module procedure i_hdf5_add_scalar
        module procedure l_hdf5_add_scalar
    end interface hdf5_add_scalar

    interface hdf5_read_scalar
        module procedure r_hdf5_read_scalar
        module procedure i_hdf5_read_scalar
        module procedure l_hdf5_read_scalar
    end interface hdf5_read_scalar
    
    interface hdf5_add_array
        ! Write a scalar as an attribute of the group
        module procedure r1_hdf5_add_array
        module procedure r2_hdf5_add_array
        module procedure r3_hdf5_add_array

        module procedure i1_hdf5_add_array
        module procedure i2_hdf5_add_array
        module procedure i3_hdf5_add_array
    end interface hdf5_add_array

    interface hdf5_read_array
        module procedure r1_hdf5_read_array
        module procedure r2_hdf5_read_array
        module procedure r3_hdf5_read_array
        
        module procedure i1_hdf5_read_array
        module procedure i2_hdf5_read_array
        module procedure i3_hdf5_read_array
    end interface hdf5_read_array
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
    
    function hdf5_array_len(hid, dataset_name)
        use hdf5
        use mod_memory, only: ip

        implicit none

        integer(hid_t), intent(in) :: hid
        character(len=*), intent(in) :: dataset_name
        integer(ip) :: hdf5_array_len

        integer(hsize_t), dimension(4) :: dims, maxdims
        integer(ip) :: rank
        integer(hid_t) :: dataset, dataspace
        integer(kind=4) :: eflag
        
        call h5dopen_f(hid, dataset_name, dataset, eflag)
        call h5dget_space_f(dataset, dataspace, eflag)
        call h5sget_simple_extent_ndims_f(dataspace, rank, eflag)
        call h5sget_simple_extent_dims_f(dataspace, dims, maxdims, eflag)
        hdf5_array_len = dims(rank)
    end function
    
    subroutine r1_hdf5_read_array(hid, dataset_name, v)
        use hdf5
        use mod_mmpol, only: fatal_error
        use mod_memory, only: ip, rp, mallocate

        implicit none

        integer(hid_t), intent(in) :: hid
        character(len=*), intent(in) :: dataset_name
        real(rp), allocatable, dimension(:) :: v

        integer(hsize_t), dimension(1) :: dims, maxdims
        integer(hid_t) :: dataset, dataspace
        integer(kind=4) :: eflag
        
        write(*, *) dataset_name
        call h5dopen_f(hid, dataset_name, dataset, eflag)
        call h5dget_space_f(dataset, dataspace, eflag)
        call h5sget_simple_extent_dims_f(dataspace, dims, maxdims, eflag)
        if(.not. allocated(v)) then
            call mallocate('r1_hdf5_read_array [v]', int(dims(1), ip), v)
        else 
            if(size(v, 1) /= int(dims(1), ip)) then
                call fatal_error("Reading HDF5 array on a buffer of wrong size")
            end if
        end if
        call h5dread_f(dataset, H5T_RP, v, dims, eflag)
        write(*, *) "OK"
    end subroutine

    subroutine r2_hdf5_read_array(hid, dataset_name, v)
        use hdf5
        use mod_mmpol, only: fatal_error
        use mod_memory, only: ip, rp, mallocate

        implicit none

        integer(hid_t), intent(in) :: hid
        character(len=*), intent(in) :: dataset_name
        real(rp), allocatable, dimension(:,:) :: v

        integer(hsize_t), dimension(2) :: dims, maxdims
        integer(hid_t) :: dataset, dataspace
        integer(kind=4) :: eflag
        
        write(*, *) dataset_name
        call h5dopen_f(hid, dataset_name, dataset, eflag)
        call h5dget_space_f(dataset, dataspace, eflag)
        call h5sget_simple_extent_dims_f(dataspace, dims, maxdims, eflag)
        if(.not. allocated(v)) then
            call mallocate('r2_hdf5_read_array [v]', int(dims(1), ip), &
                           int(dims(2), ip), v)
        else 
            if(size(v, 1) /= int(dims(1), ip) .or. &
               size(v, 2) /= int(dims(2), ip)) then
                call fatal_error("Reading HDF5 array on a buffer of wrong size")
            end if
        end if
        call h5dread_f(dataset, H5T_RP, v, dims, eflag)
    end subroutine

    subroutine r3_hdf5_read_array(hid, dataset_name, v)
        use hdf5
        use mod_mmpol, only: fatal_error
        use mod_memory, only: ip, rp, mallocate

        implicit none

        integer(hid_t), intent(in) :: hid
        character(len=*), intent(in) :: dataset_name
        real(rp), allocatable, dimension(:,:,:) :: v

        integer(hsize_t), dimension(3) :: dims, maxdims
        integer(hid_t) :: dataset, dataspace
        integer(kind=4) :: eflag
        
        write(*, *) dataset_name
        call h5dopen_f(hid, dataset_name, dataset, eflag)
        call h5dget_space_f(dataset, dataspace, eflag)
        call h5sget_simple_extent_dims_f(dataspace, dims, maxdims, eflag)
        if(.not. allocated(v)) then
            call mallocate('r3_hdf5_read_array [v]', int(dims(1), ip), &
                           int(dims(2), ip), int(dims(3), ip), v)
        else 
            if(size(v, 1) /= int(dims(1), ip) .or. &
               size(v, 2) /= int(dims(2), ip) .or. &
               size(v, 3) /= int(dims(3), ip)) then
                call fatal_error("Reading HDF5 array on a buffer of wrong size")
            end if
        end if
        call h5dread_f(dataset, H5T_RP, v, dims, eflag)
    end subroutine
    
    subroutine i1_hdf5_read_array(hid, dataset_name, v)
        use hdf5
        use mod_mmpol, only: fatal_error
        use mod_memory, only: ip, rp, mallocate

        implicit none

        integer(hid_t), intent(in) :: hid
        character(len=*), intent(in) :: dataset_name
        integer(ip), allocatable, dimension(:) :: v

        integer(hsize_t), dimension(1) :: dims, maxdims
        integer(hid_t) :: dataset, dataspace
        integer(kind=4) :: eflag
        
        write(*, *) dataset_name
        call h5dopen_f(hid, dataset_name, dataset, eflag)
        call h5dget_space_f(dataset, dataspace, eflag)
        call h5sget_simple_extent_dims_f(dataspace, dims, maxdims, eflag)
        if(.not. allocated(v)) then
            call mallocate('i1_hdf5_read_array [v]', int(dims(1), ip), v)
        else 
            if(size(v, 1) /= int(dims(1), ip)) then
                call fatal_error("Reading HDF5 array on a buffer of wrong size")
            end if
        end if
        call h5dread_f(dataset, H5T_IP, v, dims, eflag)
    end subroutine

    subroutine i2_hdf5_read_array(hid, dataset_name, v)
        use hdf5
        use mod_mmpol, only: fatal_error
        use mod_memory, only: ip, rp, mallocate

        implicit none

        integer(hid_t), intent(in) :: hid
        character(len=*), intent(in) :: dataset_name
        integer(ip), allocatable, dimension(:,:) :: v

        integer(hsize_t), dimension(2) :: dims, maxdims
        integer(hid_t) :: dataset, dataspace
        integer(kind=4) :: eflag
        
        write(*, *) dataset_name
        call h5dopen_f(hid, dataset_name, dataset, eflag)
        call h5dget_space_f(dataset, dataspace, eflag)
        call h5sget_simple_extent_dims_f(dataspace, dims, maxdims, eflag)
        if(.not. allocated(v)) then
            call mallocate('i2_hdf5_read_array [v]', int(dims(1), ip), &
                           int(dims(2), ip), v)
        else 
            if(size(v, 1) /= int(dims(1), ip) .or. &
               size(v, 2) /= int(dims(2), ip)) then
                call fatal_error("Reading HDF5 array on a buffer of wrong size")
            end if
        end if
        call h5dread_f(dataset, H5T_IP, v, dims, eflag)
    end subroutine

    subroutine i3_hdf5_read_array(hid, dataset_name, v)
        use hdf5
        use mod_mmpol, only: fatal_error
        use mod_memory, only: ip, rp, mallocate

        implicit none

        integer(hid_t), intent(in) :: hid
        character(len=*), intent(in) :: dataset_name
        integer(ip), allocatable, dimension(:,:,:) :: v

        integer(hsize_t), dimension(3) :: dims, maxdims
        integer(hid_t) :: dataset, dataspace
        integer(kind=4) :: eflag
        
        write(*, *) dataset_name
        call h5dopen_f(hid, dataset_name, dataset, eflag)
        call h5dget_space_f(dataset, dataspace, eflag)
        call h5sget_simple_extent_dims_f(dataspace, dims, maxdims, eflag)
        if(.not. allocated(v)) then
            call mallocate('i3_hdf5_read_array [v]', int(dims(1), ip), &
                           int(dims(2), ip), int(dims(3), ip), v)
        else 
            if(size(v, 1) /= int(dims(1), ip) .or. &
               size(v, 2) /= int(dims(2), ip) .or. &
               size(v, 3) /= int(dims(3), ip)) then
                call fatal_error("Reading HDF5 array on a buffer of wrong size")
            end if
        end if
        call h5dread_f(dataset, H5T_IP, v, dims, eflag)
    end subroutine

    subroutine r_hdf5_read_scalar(hid, location, attname, s)
        use hdf5
        use mod_memory, only: ip, rp

        implicit none

        integer(hid_t), intent(in) :: hid
        character(len=*), intent(in) :: location, attname
        real(rp) :: s
        
        integer(hsize_t), dimension(1), parameter :: dims = [1]
        integer(kind=4) :: eflag
        integer(hid_t) :: att_id, dataset
        
        call h5gopen_f(hid, location, dataset, eflag)
        call H5Aopen_name_f(dataset, attname, att_id, eflag)
        call H5Aread_f(att_id, H5T_RP, s, dims, eflag)
        call h5gclose_f(dataset, eflag)
    end subroutine
    
    subroutine i_hdf5_read_scalar(hid, location, attname, s)
        use hdf5
        use mod_memory, only: ip, rp

        implicit none

        integer(hid_t), intent(in) :: hid
        character(len=*), intent(in) :: location, attname
        integer(ip) :: s
        
        integer(hsize_t), dimension(1), parameter :: dims = [1]
        integer(kind=4) :: eflag
        integer(hid_t) :: att_id, dataset
        
        call h5gopen_f(hid, location, dataset, eflag)
        call H5Aopen_name_f(dataset, attname, att_id, eflag)
        call H5Aread_f(att_id, H5T_IP, s, dims, eflag)
        call h5gclose_f(dataset, eflag)
    end subroutine
    
    subroutine l_hdf5_read_scalar(hid, location, attname, s)
        use hdf5
        use mod_memory, only: ip, rp

        implicit none

        integer(hid_t), intent(in) :: hid
        character(len=*), intent(in) :: location, attname
        integer(ip) :: is
        logical :: s
        
        integer(hsize_t), dimension(1), parameter :: dims = [1]
        integer(kind=4) :: eflag
        integer(hid_t) :: att_id, dataset
       
        write(*, *) location, "--", attname
        call h5gopen_f(hid, location, dataset, eflag)
        call H5Aopen_name_f(dataset, attname, att_id, eflag)
        call H5Aread_f(att_id, H5T_IP, is, dims, eflag)
        call h5gclose_f(dataset, eflag)
        if(is == 0) then
            s = .false.
        else
            s = .true.
        end if
        write(*, *) "OK"
    end subroutine

    subroutine mmpol_save_as_hdf5(filename, out_fail)
        use hdf5
        use mod_memory, only: ip
        use mod_mmpol, only: mm_atoms, pol_atoms, cmm, polar_mm, &
                             q0, q, amoeba, pol, conn, ff_type, &
                             ix, iy, iz, mol_frame, mmat_polgrp, &
                             mscale, pscale, dscale, uscale, pscale_intra
        use mod_bonded
        use mod_nonbonded

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
        
        call h5gcreate_f(iof_hdf5_out, "system_model", hg_sysmodel, eflag)
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
        call hdf5_add_scalar(hg_sysmodel, "FF_type", ff_type)

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
        call hdf5_add_scalar(hg_cur_bp, "enabled", use_angtor)
        if(use_angtor) then
            call hdf5_add_array(hg_cur_bp, "k", angtork)
            call hdf5_add_array(hg_cur_bp, "angles_idx", angtor_a)
            call hdf5_add_array(hg_cur_bp, "torsion_idx", angtor_t)
            call hdf5_add_array(hg_cur_bp, "atoms", angtorat)
        end if
        call h5gclose_f(hg_cur_bp, eflag)
        
        ! Torsion-torsion coupling
        call h5gcreate_f(hg_cur_param, "torsion-torsion", hg_cur_bp, eflag)
        call hdf5_add_scalar(hg_cur_bp, "enabled", use_tortor)
        if(use_tortor) then
            call hdf5_add_array(hg_cur_bp, "atoms", tortorat)
            call hdf5_add_array(hg_cur_bp, "map_id", tortorprm)
            call hdf5_add_array(hg_cur_bp, "maps_ang1", ttmap_ang1)
            call hdf5_add_array(hg_cur_bp, "maps_ang2", ttmap_ang2)
            call hdf5_add_array(hg_cur_bp, "maps_pot", ttmap_v)
            call hdf5_add_array(hg_cur_bp, "maps_shapes", ttmap_shape)
        end if
        call h5gclose_f(hg_cur_bp, eflag)
        
        ! Pi-torsion
        call h5gcreate_f(hg_cur_param, "pi-torsion", hg_cur_bp, eflag)
        call hdf5_add_scalar(hg_cur_bp, "enabled", use_pitors)
        if(use_pitors) then
            call hdf5_add_array(hg_cur_bp, "atoms", pitorsat)
            call hdf5_add_array(hg_cur_bp, "k", kpitors)
        end if
        call h5gclose_f(hg_cur_bp, eflag)

        ! Out-of-plane bending
        call h5gcreate_f(hg_cur_param, "out-of-plane-bending", hg_cur_bp, eflag)
        call hdf5_add_scalar(hg_cur_bp, "enabled", use_opb)
        if(use_opb) then
            call hdf5_add_scalar(hg_cur_bp, "cubic", opb_cubic)
            call hdf5_add_scalar(hg_cur_bp, "quartic", opb_quartic)
            call hdf5_add_scalar(hg_cur_bp, "pentic", opb_pentic)
            call hdf5_add_scalar(hg_cur_bp, "sextic", opb_sextic)
            call hdf5_add_array(hg_cur_bp, "k", kopb)
            call hdf5_add_array(hg_cur_bp, "atoms", opbat)
        endif
        call h5gclose_f(hg_cur_bp, eflag)
        
        ! Urey-Bradley stretching
        call h5gcreate_f(hg_cur_param, "urey-bradley", hg_cur_bp, eflag)
        call hdf5_add_scalar(hg_cur_bp, "enabled", use_urey)
        if(use_urey) then
            call hdf5_add_scalar(hg_cur_bp, "cubic", urey_cubic)
            call hdf5_add_scalar(hg_cur_bp, "quartic", urey_quartic)
            call hdf5_add_array(hg_cur_bp, "k", kurey)
            call hdf5_add_array(hg_cur_bp, "l0", l0urey)
            call hdf5_add_array(hg_cur_bp, "atoms", ureyat)
        end if
        call h5gclose_f(hg_cur_bp, eflag)
        call h5gclose_f(hg_cur_param, eflag)
        
        ! Non-bonded (VdW and similar, all except electrostatics)
        call h5gcreate_f(hg_cur, "non-bonded", hg_cur_param, eflag)
        call hdf5_add_scalar(hg_cur_param, "enabled", use_nonbonded)
        if(use_nonbonded) then
            call hdf5_add_array(hg_cur_param, "screening", vdw_screening)
            call hdf5_add_array(hg_cur_param, "radius", vdw_r)
            call hdf5_add_array(hg_cur_param, "energy", vdw_e)
            call hdf5_add_array(hg_cur_param, "scale_factor", vdw_f)
            call hdf5_add_array(hg_cur_param, "pair_row_idx", vdw_pair%ri)
            call hdf5_add_array(hg_cur_param, "pair_col_idx", vdw_pair%ci)
            call hdf5_add_array(hg_cur_param, "pair_radius", vdw_pair_r)
            call hdf5_add_array(hg_cur_param, "pair_energy", vdw_pair_e)
        end if
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
        call hdf5_add_array(hg_cur_param, "fixed_fixed_scale_f", mscale)
        call hdf5_add_array(hg_cur_param, "fixed_ipd_scale_f", pscale)
        call hdf5_add_array(hg_cur_param, "ipd_ipd_scale_f", uscale)
        if(amoeba) then
            call hdf5_add_array(hg_cur_param, "fixed_direct_ipd_scale_f", dscale)
            call hdf5_add_array(hg_cur_param, "fixed_intragroup_ipd_scale_f", pscale_intra)
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
        call hdf5_add_array(iof_hdf5_out, "coordinates", cmm)
        
        if( eflag /= 0) then 
            write(iof_mmpol, *) "Error while creating group 'system model.'&
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

    subroutine mmpol_init_from_hdf5(filename, out_fail)
        use hdf5
        use mod_adjacency_mat, only: build_conn_upto_n, yale_sparse
        use mod_memory, only: ip, rp
        use mod_mmpol, only: mm_atoms, pol_atoms, cmm, polar_mm, &
                             q0, q, amoeba, pol, conn, ff_type, &
                             ix, iy, iz, mol_frame, mmat_polgrp, &
                             mmpol_init, set_screening_parameters, &
                             mmpol_prepare, pg_conn
        use mod_bonded
        use mod_nonbonded
        use mod_constants, only: OMMP_FF_AMOEBA

        implicit none

        character(len=*), intent(in) :: filename
        integer(ip), intent(out) :: out_fail
        
        integer(hid_t) :: hg_sysmodel, hg_res, hg_cur, hg_amoeba, &
                          hg_top, hg_cur_param, hg_cur_bp, cur_dst, cur_dsp
        integer(hsize_t), dimension(4) :: dims
        integer(kind=4) :: eflag
        real(rp), dimension(:), allocatable :: l_mscale, l_pscale, l_dscale, &
                                               l_uscale, l_ipscale, l_vdwscale
        type(yale_sparse) :: conn_1

        ! For handling torsion maps
        integer(ip) :: i, ibeg, iend
        integer(ip), allocatable, dimension(:,:) :: tmp_shape
        real(rp), allocatable, dimension(:) :: tmp_ang1, tmp_ang2, tmp_v

        ! Initialize interface
        call h5open_f(eflag)
        if(eflag /= 0) then
            write(iof_mmpol, *) "Unable to initialize HDF5 module. Failure in &
                               &h5open_f subroutine."
            out_fail = -1_ip
            return
        end if
        
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, iof_hdf5_io, eflag)
        if( eflag /= 0) then 
            write(iof_mmpol, *) "Unable to open HDF5 file. Failure in &
                               &h5fopen_f subroutine."
            out_fail = -1_ip
            return
        end if
    
        call hdf5_read_scalar(iof_hdf5_io, &
                              'system_model', 'N-atoms', &
                              mm_atoms)
        call hdf5_read_scalar(iof_hdf5_io, &
                              'system_model', 'N-pol-atoms', &
                              pol_atoms)
        call hdf5_read_scalar(iof_hdf5_io, &
                              'system_model', 'FF_type', &
                              ff_type)
        
        ! Initialize mmpol module
        call mmpol_init(ff_type, mm_atoms, pol_atoms)
        
        ! Connectivity 
        call hdf5_read_array(iof_hdf5_io, &
                             'system_model/topology/connectivity/ADJ1-RowIdx', &
                             conn_1%ri)
        call hdf5_read_array(iof_hdf5_io, &
                             'system_model/topology/connectivity/ADJ1-ColIdx', &
                             conn_1%ci)
        conn_1%n = size(conn_1%ri) - 1
        call build_conn_upto_n(conn_1, 4, conn, .false.)

        ! AMOEBA
        if(amoeba) then
            call hdf5_read_array(iof_hdf5_io, &
                                 'system_model/topology/AMOEBA/polarization_group_id', &
                                  mmat_polgrp)
        end if
        
        ! Bonded Parameters
        ! Bond stretching
        call hdf5_read_scalar(iof_hdf5_io, &
                              'system_model/parameters/bonded/stretching', &
                              'enabled', &
                              use_bond)
        if(use_bond) then
            call bond_init(hdf5_array_len(iof_hdf5_io, &
                                          'system_model/parameters/bonded/stretching/atoms'))

            call hdf5_read_scalar(iof_hdf5_io, &
                                  'system_model/parameters/bonded/stretching', &
                                  'quartic', &
                                  bond_quartic)
            call hdf5_read_scalar(iof_hdf5_io, &
                                  'system_model/parameters/bonded/stretching', &
                                  'cubic', &
                                  bond_cubic)
            call hdf5_read_array(iof_hdf5_io, &
                                 'system_model/parameters/bonded/stretching/k', &
                                 kbond)
            call hdf5_read_array(iof_hdf5_io, &
                                 'system_model/parameters/bonded/stretching/l0', &
                                 l0bond)
            call hdf5_read_array(iof_hdf5_io, &
                                 'system_model/parameters/bonded/stretching/atoms', &
                                 bondat)
        end if
        
        ! Angle bending
        call hdf5_read_scalar(iof_hdf5_io, &
                              'system_model/parameters/bonded/bending', &
                              'enabled', &
                              use_angle)
        if(use_angle) then
            call angle_init(hdf5_array_len(iof_hdf5_io, &
                                          'system_model/parameters/bonded/bending/atoms'))
            call hdf5_read_scalar(iof_hdf5_io, &
                                  'system_model/parameters/bonded/bending', &
                                  "cubic", angle_cubic)
            call hdf5_read_scalar(iof_hdf5_io, &
                                  'system_model/parameters/bonded/bending', &
                                  "quartic", angle_quartic)
            call hdf5_read_scalar(iof_hdf5_io, &
                                  'system_model/parameters/bonded/bending', &
                                  "pentic", angle_pentic)
            call hdf5_read_scalar(iof_hdf5_io, &
                                  'system_model/parameters/bonded/bending', &
                                  "sextic", angle_sextic)
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/bending/k", &
                                 kangle)
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/bending/ang0", &
                                 eqangle)
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/bending/atoms", &
                                 angleat)
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/bending/type", &
                                 anglety)
        end if
        
        ! Dihedral torsion
        call hdf5_read_scalar(iof_hdf5_io, &
                              'system_model/parameters/bonded/torsion', &
                              'enabled', &
                              use_torsion)
        if(use_torsion) then
            call torsion_init(hdf5_array_len(iof_hdf5_io, &
                                          'system_model/parameters/bonded/torsion/atoms'))
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/torsion/amplitudes", &
                                 torsamp)
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/torsion/phase", &
                                 torsphase)
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/torsion/atoms", &
                                 torsionat)
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/torsion/period", &
                                 torsn)
        end if
        
        ! Stretching-bending coupling
        call hdf5_read_scalar(iof_hdf5_io, &
                              'system_model/parameters/bonded/stretching-bending', &
                              'enabled', &
                              use_strbnd)
        if(use_strbnd) then
            call strbnd_init(hdf5_array_len(iof_hdf5_io, &
                                          'system_model/parameters/bonded/stretching-bending/atoms'))
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/stretching-bending/k1", &
                                 strbndk1)
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/stretching-bending/k2", &
                                 strbndk2)
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/stretching-bending/l1_0", &
                                 strbndl10)
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/stretching-bending/l2_0", &
                                 strbndl20)
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/stretching-bending/ang0", &
                                 strbndthet0)
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/stretching-bending/atoms", &
                                 strbndat)
        end if
        
        ! Stretching-torsion coupling
        call hdf5_read_scalar(iof_hdf5_io, &
                              'system_model/parameters/bonded/stretching-torsion', &
                              'enabled', &
                              use_strtor)
        if(use_strtor) then
            call strtor_init(hdf5_array_len(iof_hdf5_io, &
                                          'system_model/parameters/bonded/stretching-torsion/atoms'))
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/stretching-torsion/k", &
                                 strtork)
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/stretching-torsion/bonds_idx", &
                                 strtor_b)
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/stretching-torsion/torsion_idx", &
                                 strtor_t)
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/stretching-torsion/atoms", &
                                 strtorat)
        end if

        ! Bending-torsion coupling 
        call hdf5_read_scalar(iof_hdf5_io, &
                              'system_model/parameters/bonded/bending-torsion', &
                              'enabled', &
                              use_angtor)
        if(use_angtor) then
            call angtor_init(hdf5_array_len(iof_hdf5_io, &
                                          'system_model/parameters/bonded/bending-torsion/atoms'))
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/bending-torsion/k", &
                                 angtork)
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/bending-torsion/angles_idx", &
                                 angtor_a)
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/bending-torsion/torsion_idx", &
                                 angtor_t)
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/bending-torsion/atoms", &
                                 angtorat)
        end if
        
        ! Torsion-torsion coupling
        call hdf5_read_scalar(iof_hdf5_io, &
                              'system_model/parameters/bonded/torsion-torsion', &
                              'enabled', &
                              use_tortor)
        if(use_tortor) then
            call tortor_init(hdf5_array_len(iof_hdf5_io, &
                                          'system_model/parameters/bonded/torsion-torsion/atoms'))
            call hdf5_read_array(iof_hdf5_io, &
                                "system_model/parameters/bonded/torsion-torsion/atoms", &
                                tortorat)
            call hdf5_read_array(iof_hdf5_io, &
                                "system_model/parameters/bonded/torsion-torsion/map_id", &
                                tortorprm)
            call hdf5_read_array(iof_hdf5_io, &
                                "system_model/parameters/bonded/torsion-torsion/maps_ang1", &
                                tmp_ang1)
            call hdf5_read_array(iof_hdf5_io, &
                                "system_model/parameters/bonded/torsion-torsion/maps_ang2", &
                                tmp_ang2)
            call hdf5_read_array(iof_hdf5_io, &
                                "system_model/parameters/bonded/torsion-torsion/maps_pot", &
                                tmp_v)
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/torsion-torsion/maps_shapes", &
                                 tmp_shape)
            
            iend = 0
            do i=1, size(tmp_shape, 2)
                ibeg = iend + 1 
                iend = ibeg + tmp_shape(1,i) * tmp_shape(2,i) - 1
                write(*, *) ibeg, iend, shape(tmp_ang1)

                call tortor_newmap(tmp_shape(1,i), &
                                   tmp_shape(2,i), &
                                   tmp_ang1(ibeg:iend), &
                                   tmp_ang2(ibeg:iend), &
                                   tmp_v(ibeg:iend))
            end do
        end if
        
        ! Pi-torsion
        call hdf5_read_scalar(iof_hdf5_io, &
                              'system_model/parameters/bonded/pi-torsion', &
                              'enabled', &
                              use_pitors)
        if(use_pitors) then
            call pitors_init(hdf5_array_len(iof_hdf5_io, &
                                          'system_model/parameters/bonded/pi-torsion/atoms'))
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/pi-torsion/atoms", &
                                 pitorsat)
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/pi-torsion/k", &
                                 kpitors)
        end if

        ! Out-of-plane bending
        call hdf5_read_scalar(iof_hdf5_io, &
                              'system_model/parameters/bonded/out-of-plane-bending', &
                              'enabled', &
                              use_opb)
        if(use_opb) then
            call opb_init(hdf5_array_len(iof_hdf5_io, &
                                          'system_model/parameters/bonded/out-of-plane-bending/atoms'), 'allinger')
            call hdf5_read_scalar(iof_hdf5_io, &
                                  'system_model/parameters/bonded/out-of-plane-bending', &
                                  "cubic", opb_cubic)
            call hdf5_read_scalar(iof_hdf5_io, &
                                  'system_model/parameters/bonded/out-of-plane-bending', &
                                  "quartic", opb_quartic)
            call hdf5_read_scalar(iof_hdf5_io, &
                                  'system_model/parameters/bonded/out-of-plane-bending', &
                                  "pentic", opb_pentic)
            call hdf5_read_scalar(iof_hdf5_io, &
                                  'system_model/parameters/bonded/out-of-plane-bending', &
                                  "sextic", opb_sextic)
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/out-of-plane-bending/k", &
                                 kopb)
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/out-of-plane-bending/atoms", &
                                 opbat)
        endif
        
        ! Urey-Bradley stretching
        call hdf5_read_scalar(iof_hdf5_io, &
                              'system_model/parameters/bonded/urey-bradley', &
                              'enabled', &
                              use_urey)
        if(use_urey) then
            call urey_init(hdf5_array_len(iof_hdf5_io, &
                                          'system_model/parameters/bonded/urey-bradley/atoms'))
            call hdf5_read_scalar(iof_hdf5_io, &
                                  'system_model/parameters/bonded/urey-bradley', &
                                  "cubic", urey_cubic)
            call hdf5_read_scalar(iof_hdf5_io, &
                                  'system_model/parameters/bonded/urey-bradley', &
                                  "quartic", urey_quartic)
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/urey-bradley/k", &
                                 kurey)
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/urey-bradley/l0", &
                                 l0urey)
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/bonded/urey-bradley/atoms", &
                                 ureyat)
        end if
        
        ! Non-bonded (VdW and similar, all except electrostatics)
        call hdf5_read_scalar(iof_hdf5_io, &
                              'system_model/parameters/non-bonded', &
                              'enabled', &
                              use_nonbonded)
        if(use_nonbonded) then
            call vdw_init("buffered-14-7", "cubic-mean", "diameter", "r-min", &
                          "hhg")
            call hdf5_read_array(iof_hdf5_io, & 
                                 "system_model/parameters/non-bonded/screening", &
                                 l_vdwscale)
            vdw_screening = l_vdwscale

            call hdf5_read_array(iof_hdf5_io, & 
                                 "system_model/parameters/non-bonded/radius", &
                                 vdw_r)
            call hdf5_read_array(iof_hdf5_io, & 
                                 "system_model/parameters/non-bonded/energy", &
                                 vdw_e)
            call hdf5_read_array(iof_hdf5_io, & 
                                 "system_model/parameters/non-bonded/scale_factor", &
                                 vdw_f)
            call hdf5_read_array(iof_hdf5_io, & 
                                 "system_model/parameters/non-bonded/pair_row_idx", &
                                 vdw_pair%ri)
            vdw_pair%n = size(vdw_pair%ri)

            call hdf5_read_array(iof_hdf5_io, & 
                                 "system_model/parameters/non-bonded/pair_col_idx", &
                                 vdw_pair%ci)
            call hdf5_read_array(iof_hdf5_io, & 
                                 "system_model/parameters/non-bonded/pair_radius", &
                                 vdw_pair_r)
            call hdf5_read_array(iof_hdf5_io, & 
                                 "system_model/parameters/non-bonded/pair_energy", &
                                 vdw_pair_e)
        end if
        
        if(amoeba) then
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/electrostatics/fixed_multipoles", q0)
            call hdf5_read_array(iof_hdf5_io, &
                                "system_model/parameters/electrostatics/AMOEBA/fixed_mmpoles_rot_Z", &
                                iz)
            call hdf5_read_array(iof_hdf5_io, &
                                "system_model/parameters/electrostatics/AMOEBA/fixed_mmpoles_rot_X", &
                                ix)
            call hdf5_read_array(iof_hdf5_io, &
                                "system_model/parameters/electrostatics/AMOEBA/fixed_mmpoles_rot_Y", &
                                iy)
            call hdf5_read_array(iof_hdf5_io, &
                                "system_model/parameters/electrostatics/AMOEBA/fixed_mmpoles_rot_CONV", &
                                mol_frame)
        else
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/electrostatics/fixed_multipoles", q)
        end if
        call hdf5_read_array(iof_hdf5_io, &
                             "system_model/parameters/electrostatics/fixed_fixed_scale_f", l_mscale)
        call hdf5_read_array(iof_hdf5_io, &
                             "system_model/parameters/electrostatics/fixed_ipd_scale_f", l_pscale)
        call hdf5_read_array(iof_hdf5_io, &
                             "system_model/parameters/electrostatics/ipd_ipd_scale_f", l_uscale)
        if(amoeba) then
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/electrostatics/fixed_direct_ipd_scale_f", &
                                 l_dscale)
            call hdf5_read_array(iof_hdf5_io, &
                                 "system_model/parameters/electrostatics/fixed_intragroup_ipd_scale_f", &
                                 l_ipscale)
        end if
        call set_screening_parameters(l_mscale, l_pscale, l_dscale, l_uscale, &
                                      l_ipscale)
        
        call hdf5_read_array(iof_hdf5_io, &
                             "system_model/parameters/electrostatics/polarizable_atoms_idx", polar_mm)
        call hdf5_read_array(iof_hdf5_io, &
                             "system_model/parameters/electrostatics/polarizabilities", pol)
        
        call hdf5_read_array(iof_hdf5_io, "coordinates", cmm)

        call h5fclose_f(iof_hdf5_io, eflag)
        if( eflag /= 0) then 
            write(iof_mmpol, *) "Error while closing HDF5 file. Failure in &
                               &h5fclose_f subroutine."
            out_fail = -1_ip
            return
        end if
        
        call mmpol_prepare()

        out_fail = 0_ip
    end subroutine mmpol_init_from_hdf5
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
