#ifdef USE_HDF5

#define H5T_RP H5T_NATIVE_DOUBLE
#define H5T_LOGICAL H5T_BITFIELD_F
#ifdef USE_I8
#define H5T_IP H5T_STD_I64LE
#else
#define H5T_IP H5T_STD_I32LE
#endif

module mod_iohdf5
    use hdf5
    use mod_memory, only: ip, rp
    use mod_mmpol, only: ommp_system
    use mod_topology, only: ommp_topology_type
    use mod_electrostatics, only: ommp_electrostatics_type
    use mod_nonbonded, only: ommp_nonbonded_type, vdw_init
    use mod_bonded, only: ommp_bonded_type
    use mod_constants, only: OMMP_VERBOSE_LOW, OMMP_VERBOSE_HIGH
    use mod_io, only: ommp_message, fatal_error
    
    implicit none
    private

    public :: save_system_as_hdf5, mmpol_init_from_hdf5
    
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
    
    contains
    
    ! Subroutines dedicated to HDF5 I/O

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
        
        call h5dopen_f(hid, dataset_name, dataset, eflag)
        call h5dget_space_f(dataset, dataspace, eflag)
        call h5sget_simple_extent_dims_f(dataspace, dims, maxdims, eflag)
        if(.not. allocated(v)) then
            call mallocate('r1_hdf5_read_array [v]', int(dims(1), ip), v)
        else 
            if(size(v, 1) /= int(dims(1), ip)) then
                call fatal_error("Reading HDF5 array on a buffer of wrong size ["//dataset_name//"]")
            end if
        end if
        call h5dread_f(dataset, H5T_RP, v, dims, eflag)
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
        
        call h5dopen_f(hid, dataset_name, dataset, eflag)
        call h5dget_space_f(dataset, dataspace, eflag)
        call h5sget_simple_extent_dims_f(dataspace, dims, maxdims, eflag)
        if(.not. allocated(v)) then
            call mallocate('r2_hdf5_read_array [v]', int(dims(1), ip), &
                           int(dims(2), ip), v)
        else 
            if(size(v, 1) /= int(dims(1), ip) .or. &
               size(v, 2) /= int(dims(2), ip)) then
                call fatal_error("Reading HDF5 array on a buffer of wrong size ["//dataset_name//"]")
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
                call fatal_error("Reading HDF5 array on a buffer of wrong size ["//dataset_name//"]")
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
        
        call h5dopen_f(hid, dataset_name, dataset, eflag)
        call h5dget_space_f(dataset, dataspace, eflag)
        call h5sget_simple_extent_dims_f(dataspace, dims, maxdims, eflag)
        if(.not. allocated(v)) then
            call mallocate('i1_hdf5_read_array [v]', int(dims(1), ip), v)
        else 
            if(size(v, 1) /= int(dims(1), ip)) then
                call fatal_error("Reading HDF5 array on a buffer of wrong size ["//dataset_name//"]")
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
        
        call h5dopen_f(hid, dataset_name, dataset, eflag)
        call h5dget_space_f(dataset, dataspace, eflag)
        call h5sget_simple_extent_dims_f(dataspace, dims, maxdims, eflag)
        if(.not. allocated(v)) then
            call mallocate('i2_hdf5_read_array [v]', int(dims(1), ip), &
                           int(dims(2), ip), v)
        else 
            if(size(v, 1) /= int(dims(1), ip) .or. &
               size(v, 2) /= int(dims(2), ip)) then
                call fatal_error("Reading HDF5 array on a buffer of wrong size ["//dataset_name//"]")
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
                call fatal_error("Reading HDF5 array on a buffer of wrong size ["//dataset_name//"]")
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
       
        call h5gopen_f(hid, location, dataset, eflag)
        call H5Aopen_name_f(dataset, attname, att_id, eflag)
        call H5Aread_f(att_id, H5T_IP, is, dims, eflag)
        call h5gclose_f(dataset, eflag)
        if(is == 0) then
            s = .false.
        else
            s = .true.
        end if
    end subroutine

    subroutine hdf5_group_exists(hid, location, exists)
        implicit none 
        
        integer(hid_t), intent(in) :: hid
        character(len=*), intent(in) :: location
        logical, intent(out) :: exists
        
        integer(kind=4) :: eflag
        integer(hid_t) :: grp_id
       
        call h5gopen_f(hid, location, grp_id, eflag)
        if(eflag == 0) then
            exists = .true.
            call h5gclose_f(grp_id, eflag)
        else
            exists = .false.
        end if
    end subroutine


    subroutine save_system_as_hdf5(filename, s, out_fail, & 
                                   namespace, &
                                   mutable_only)
        implicit none

        character(len=*), intent(in) :: filename, namespace
        type(ommp_system), intent(in) :: s
        integer(ip), intent(out) :: out_fail
        logical, intent(in) :: mutable_only
        
        integer(hid_t) :: hg
        integer(hsize_t), dimension(4) :: dims
        integer(kind=4) :: eflag
        integer(hid_t) :: iof_hdf5 = 301
        logical :: append

        ! Initialize interface
        call h5open_f(eflag)
        if(eflag /= 0) then
            call ommp_message("Unable to initialize HDF5 module. Failure in &
                               &h5open_f subroutine.", OMMP_VERBOSE_LOW)
            out_fail = -1_ip
            return
        end if

        inquire(file=filename, exist=append)
        if(append) then
            call ommp_message("HDF5 file exists, appending.", OMMP_VERBOSE_HIGH)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, iof_hdf5, eflag)
            if( eflag /= 0) then 
                call ommp_message("Unable to open HDF5 file. Failure in &
                                  &h5opene_f subroutine.", OMMP_VERBOSE_LOW)
                out_fail = -1_ip
                return
            end if
        else
            call h5fcreate_f(filename, H5F_ACC_EXCL_F, iof_hdf5, eflag)
            if( eflag /= 0) then 
                call ommp_message("Unable to create HDF5 file. Failure in &
                                  &h5fcreate_f subroutine.", OMMP_VERBOSE_LOW)
                out_fail = -1_ip
                return
            end if
        end if
        
        ! TODO Handle more complex cases like a/b/c in namespace
        call h5gcreate_f(iof_hdf5, namespace, hg, eflag)
        if( eflag /= 0) then 
            call ommp_message("Error while creating group.&
                              &Failure in h5gcreate_f subroutine.", OMMP_VERBOSE_LOW)
            out_fail = -1_ip
            return
        end if
        call hdf5_add_scalar(hg, 'mutable_only', mutable_only)
        call h5gclose_f(hg, eflag)
       
        call save_topology_as_hdf5(iof_hdf5, s%top, eflag, &
                                   namespace//'/topology', mutable_only)
        call save_electrostatics_as_hdf5(iof_hdf5, s%eel, eflag, &
                                         namespace//'/electrostatics', & 
                                         mutable_only)

        if(s%use_nonbonded) &
            call save_nonbonded_as_hdf5(iof_hdf5, s%vdw, eflag, & 
                                        namespace//'/nonbonded', &
                                        mutable_only)
        if(s%use_bonded) &
            call save_bonded_as_hdf5(iof_hdf5, s%bds, eflag, & 
                                     namespace//'/bonded', &
                                     mutable_only)

        call h5fclose_f(iof_hdf5, eflag)
        if( eflag /= 0) then 
            call ommp_message("Error while closing HDF5 file. Failure in &
                               &h5fclose_f subroutine.", OMMP_VERBOSE_LOW)
            out_fail = -1_ip
            return
        end if

        out_fail = 0_ip
    end subroutine

    subroutine save_topology_as_hdf5(iof_hdf5, top, out_fail, &
                                     namespace, &
                                     mutable_only)
        implicit none

        integer(hid_t) :: iof_hdf5
        character(len=*), intent(in) :: namespace
        type(ommp_topology_type), intent(in) :: top
        integer(ip), intent(out) :: out_fail
        logical, intent(in) :: mutable_only
        
        integer(hid_t) :: hg, hg_cur
        integer(hsize_t), dimension(4) :: dims
        integer(kind=4) :: eflag

        call h5gcreate_f(iof_hdf5, namespace, hg, eflag)
        if( eflag /= 0) then 
            call ommp_message("Error while creating group 'system model.'&
                              &Failure in h5gcreate_f subroutine.", OMMP_VERBOSE_LOW)
            out_fail = -1_ip
            return
        end if
        
        if(.not. mutable_only) then
            call hdf5_add_scalar(hg, "N-atoms", top%mm_atoms)
            call h5gcreate_f(hg, "connectivity", hg_cur, eflag)
            call hdf5_add_array(hg_cur, "ADJ1-RowIdx", top%conn(1)%ri) 
            call hdf5_add_array(hg_cur, "ADJ1-ColIdx", top%conn(1)%ci) 
            call h5gclose_f(hg_cur, eflag)

            if(top%atz_initialized) call hdf5_add_array(hg, "Atoms-Z", top%atz)
            if(top%attype_initialized) &
                 call hdf5_add_array(hg, "Atoms-Type", top%attype)
            if(top%atclass_initialized) &
                call hdf5_add_array(hg, "Atoms-Class", top%atclass)
        end if
        
        call hdf5_add_array(hg, "Atoms-Coordinates", top%cmm)
        call h5gclose_f(hg, eflag)

        out_fail = 0_ip

    end subroutine
    
    subroutine save_electrostatics_as_hdf5(iof_hdf5, eel, out_fail, &
                                           namespace, &
                                           mutable_only)

        use mod_memory, only: mallocate, mfree
        implicit none

        integer(hid_t) :: iof_hdf5
        character(len=*), intent(in) :: namespace
        type(ommp_electrostatics_type), intent(in) :: eel
        integer(ip), intent(out) :: out_fail
        logical, intent(in) :: mutable_only
        
        integer(hid_t) :: hg, hg_cur
        integer(hsize_t), dimension(4) :: dims
        integer(kind=4) :: eflag
        real(rp), allocatable :: tmp_q0(:, :)

        call h5gcreate_f(iof_hdf5, namespace, hg, eflag)
        if( eflag /= 0) then 
            call ommp_message("Error while creating group. &
                              &Failure in h5gcreate_f subroutine.", OMMP_VERBOSE_LOW)
            out_fail = -1_ip
            return
        end if
        
        if(.not. mutable_only) then
            call hdf5_add_scalar(hg, "N-pol-atoms", eel%pol_atoms)
            call hdf5_add_scalar(hg, "amoeba", eel%amoeba)
            
            if(eel%amoeba) then
                ! Write the unrotated multipoles, that are coordinates independent
                ! Since quadrupoles have been multiplied by 1/3 
                ! (see [[mod_mmpol::mmpol_prepare]]) we should now save them 
                ! correctly so we multiply quadrupoles terms by 3
                call mallocate('save_as_hdf5 [tmp_q0]', &
                               size(eel%q0, 1), size(eel%q0, 2), tmp_q0)
                tmp_q0 = eel%q0
                tmp_q0(5:10,:) = tmp_q0(5:10,:) * 3.0
                call hdf5_add_array(hg, "fixed_multipoles_unrotated", tmp_q0)
                call mfree('save_as_hdf5 [tmp_q0]', tmp_q0)
                
                ! Write all the information needed to perform the rotation of the 
                ! multipoles
                call hdf5_add_array(hg, "fixed_mmpoles_rot_Z", eel%iz)
                call hdf5_add_array(hg, "fixed_mmpoles_rot_X", eel%ix)
                call hdf5_add_array(hg, "fixed_mmpoles_rot_Y", eel%iy)
                call hdf5_add_array(hg, "fixed_mmpoles_rot_CONV", eel%mol_frame)

                call hdf5_add_array(hg, "polarization_group_id", eel%mmat_polgrp)
            else
                call hdf5_add_array(hg, "fixed_multipoles", eel%q)
            end if
            
            call hdf5_add_array(hg, "fixed_fixed_scale_f", eel%mscale)
            call hdf5_add_array(hg, "fixed_ipd_scale_f", eel%pscale)
            call hdf5_add_array(hg, "ipd_ipd_scale_f", eel%uscale)
            if(eel%amoeba) then
                call hdf5_add_array(hg, "fixed_direct_ipd_scale_f", eel%dscale)
                call hdf5_add_array(hg, "fixed_intragroup_ipd_scale_f", eel%pscale_intra)
            end if
            
            call hdf5_add_array(hg, "polarizable_atoms_idx", eel%polar_mm)
            call hdf5_add_array(hg, "polarizabilities", eel%pol) 
        end if
       
        if(eel%amoeba) then
             call hdf5_add_array(hg, "fixed_multipoles_rotated", eel%q)
        end if

        if(eel%M2M_done) then
            call hdf5_add_array(hg, "potential_M2M", eel%V_M2M)
            call hdf5_add_array(hg, "field_M2M", eel%E_M2M)
            call hdf5_add_array(hg, "field_grd_M2M", eel%Egrd_M2M)
        end if
        
        if(eel%M2D_done) then
            call hdf5_add_array(hg, "field_M2M", eel%E_M2D)
        end if
        
        if(eel%ipd_done) then
            call hdf5_add_array(hg, "induced_point_dipoles", eel%ipd)
        end if

        call h5gclose_f(hg, eflag)

        out_fail = 0_ip

    end subroutine

    subroutine save_nonbonded_as_hdf5(iof_hdf5, vdw, out_fail, &
                                      namespace, &
                                      mutable_only)
        implicit none

        integer(hid_t) :: iof_hdf5
        character(len=*), intent(in) :: namespace
        type(ommp_nonbonded_type), intent(in) :: vdw
        integer(ip), intent(out) :: out_fail
        logical, intent(in) :: mutable_only
        
        integer(hid_t) :: hg, hg_cur
        integer(hsize_t), dimension(4) :: dims
        integer(kind=4) :: eflag

        call h5gcreate_f(iof_hdf5, namespace, hg, eflag)
        if( eflag /= 0) then 
            call ommp_message("Error while creating group. &
                              &Failure in h5gcreate_f subroutine.", OMMP_VERBOSE_LOW)
            out_fail = -1_ip
            return
        end if
        
        if(.not. mutable_only) then
            call hdf5_add_array(hg, "screening", vdw%vdw_screening)
            call hdf5_add_array(hg, "radius", vdw%vdw_r)
            call hdf5_add_array(hg, "energy", vdw%vdw_e)
            call hdf5_add_array(hg, "scale_factor", vdw%vdw_f)
            call hdf5_add_array(hg, "pair_row_idx", vdw%vdw_pair%ri)
            call hdf5_add_array(hg, "pair_col_idx", vdw%vdw_pair%ci)
            call hdf5_add_array(hg, "pair_radius", vdw%vdw_pair_r)
            call hdf5_add_array(hg, "pair_energy", vdw%vdw_pair_e)
        end if
       
        call h5gclose_f(hg, eflag)

        out_fail = 0_ip
    end subroutine
    
    subroutine save_bonded_as_hdf5(iof_hdf5, bds, out_fail, &
                                   namespace, &
                                   mutable_only)
        implicit none

        integer(hid_t) :: iof_hdf5
        character(len=*), intent(in) :: namespace
        type(ommp_bonded_type), intent(in) :: bds
        integer(ip), intent(out) :: out_fail
        logical, intent(in) :: mutable_only
        
        integer(hid_t) :: hg, hg_cur, hg_cur_bp
        integer(hsize_t), dimension(4) :: dims
        integer(kind=4) :: eflag

        call h5gcreate_f(iof_hdf5, namespace, hg, eflag)
        if( eflag /= 0) then 
            call ommp_message("Error while creating group. &
                              &Failure in h5gcreate_f subroutine.", OMMP_VERBOSE_LOW)
            out_fail = -1_ip
            return
        end if
        
        if(.not. mutable_only) then
            call h5gcreate_f(hg, "stretching", hg_cur_bp, eflag)
            call hdf5_add_scalar(hg_cur_bp, "enabled", bds%use_bond)
            if(bds%use_bond) then
                call hdf5_add_scalar(hg_cur_bp, "cubic", bds%bond_cubic)
                call hdf5_add_scalar(hg_cur_bp, "quartic", bds%bond_quartic)
                call hdf5_add_array(hg_cur_bp, "k", bds%kbond)
                call hdf5_add_array(hg_cur_bp, "l0", bds%l0bond)
                call hdf5_add_array(hg_cur_bp, "atoms", bds%bondat)
            end if
            call h5gclose_f(hg_cur_bp, eflag)
            
            ! Angle bending
            call h5gcreate_f(hg, "bending", hg_cur_bp, eflag)
            call hdf5_add_scalar(hg_cur_bp, "enabled", bds%use_angle)
            if(bds%use_angle) then
                call hdf5_add_scalar(hg_cur_bp, "cubic", bds%angle_cubic)
                call hdf5_add_scalar(hg_cur_bp, "quartic", bds%angle_quartic)
                call hdf5_add_scalar(hg_cur_bp, "pentic", bds%angle_pentic)
                call hdf5_add_scalar(hg_cur_bp, "sextic", bds%angle_sextic)
                call hdf5_add_array(hg_cur_bp, "k", bds%kangle)
                call hdf5_add_array(hg_cur_bp, "ang0", bds%eqangle)
                call hdf5_add_array(hg_cur_bp, "atoms", bds%angleat)
                call hdf5_add_array(hg_cur_bp, "type", bds%anglety)
            end if
            call h5gclose_f(hg_cur_bp, eflag)
            
            ! Dihedral torsion
            call h5gcreate_f(hg, "torsion", hg_cur_bp, eflag)
            call hdf5_add_scalar(hg_cur_bp, "enabled", bds%use_torsion)
            if(bds%use_torsion) then
                call hdf5_add_array(hg_cur_bp, "amplitudes", bds%torsamp)
                call hdf5_add_array(hg_cur_bp, "phase", bds%torsphase)
                call hdf5_add_array(hg_cur_bp, "atoms", bds%torsionat)
                call hdf5_add_array(hg_cur_bp, "period", bds%torsn)
            end if
            call h5gclose_f(hg_cur_bp, eflag)
            
            ! Stretching-bending coupling
            call h5gcreate_f(hg, "stretching-bending", hg_cur_bp, eflag)
            call hdf5_add_scalar(hg_cur_bp, "enabled", bds%use_strbnd)
            if(bds%use_strbnd) then
                call hdf5_add_array(hg_cur_bp, "k1", bds%strbndk1)
                call hdf5_add_array(hg_cur_bp, "k2", bds%strbndk2)
                call hdf5_add_array(hg_cur_bp, "l1_0", bds%strbndl10)
                call hdf5_add_array(hg_cur_bp, "l2_0", bds%strbndl20)
                call hdf5_add_array(hg_cur_bp, "ang0", bds%strbndthet0)
                call hdf5_add_array(hg_cur_bp, "atoms", bds%strbndat)
            end if
            call h5gclose_f(hg_cur_bp, eflag)
            
            ! Stretching-torsion coupling
            call h5gcreate_f(hg, "stretching-torsion", hg_cur_bp, eflag)
            call hdf5_add_scalar(hg_cur_bp, "enabled", bds%use_strtor)
            if(bds%use_strtor) then
                call hdf5_add_array(hg_cur_bp, "k", bds%strtork)
                call hdf5_add_array(hg_cur_bp, "bonds_idx", bds%strtor_b)
                call hdf5_add_array(hg_cur_bp, "torsion_idx", bds%strtor_t)
                call hdf5_add_array(hg_cur_bp, "atoms", bds%strtorat)
            end if
            call h5gclose_f(hg_cur_bp, eflag)

            ! Bending-torsion coupling 
            call h5gcreate_f(hg, "bending-torsion", hg_cur_bp, eflag)
            call hdf5_add_scalar(hg_cur_bp, "enabled", bds%use_angtor)
            if(bds%use_angtor) then
                call hdf5_add_array(hg_cur_bp, "k", bds%angtork)
                call hdf5_add_array(hg_cur_bp, "angles_idx", bds%angtor_a)
                call hdf5_add_array(hg_cur_bp, "torsion_idx", bds%angtor_t)
                call hdf5_add_array(hg_cur_bp, "atoms", bds%angtorat)
            end if
            call h5gclose_f(hg_cur_bp, eflag)
            
            ! Torsion-torsion coupling
            call h5gcreate_f(hg, "torsion-torsion", hg_cur_bp, eflag)
            call hdf5_add_scalar(hg_cur_bp, "enabled", bds%use_tortor)
            if(bds%use_tortor) then
                call hdf5_add_array(hg_cur_bp, "atoms", bds%tortorat)
                call hdf5_add_array(hg_cur_bp, "map_id", bds%tortorprm)
                call hdf5_add_array(hg_cur_bp, "maps_ang1", bds%ttmap_ang1)
                call hdf5_add_array(hg_cur_bp, "maps_ang2", bds%ttmap_ang2)
                call hdf5_add_array(hg_cur_bp, "maps_pot", bds%ttmap_v)
                call hdf5_add_array(hg_cur_bp, "maps_shapes", bds%ttmap_shape)
            end if
            call h5gclose_f(hg_cur_bp, eflag)
            
            ! Pi-torsion
            call h5gcreate_f(hg, "pi-torsion", hg_cur_bp, eflag)
            call hdf5_add_scalar(hg_cur_bp, "enabled", bds%use_pitors)
            if(bds%use_pitors) then
                call hdf5_add_array(hg_cur_bp, "atoms", bds%pitorsat)
                call hdf5_add_array(hg_cur_bp, "k", bds%kpitors)
            end if
            call h5gclose_f(hg_cur_bp, eflag)

            ! Out-of-plane bending
            call h5gcreate_f(hg, "out-of-plane-bending", hg_cur_bp, eflag)
            call hdf5_add_scalar(hg_cur_bp, "enabled", bds%use_opb)
            if(bds%use_opb) then
                call hdf5_add_scalar(hg_cur_bp, "cubic", bds%opb_cubic)
                call hdf5_add_scalar(hg_cur_bp, "quartic", bds%opb_quartic)
                call hdf5_add_scalar(hg_cur_bp, "pentic", bds%opb_pentic)
                call hdf5_add_scalar(hg_cur_bp, "sextic", bds%opb_sextic)
                call hdf5_add_array(hg_cur_bp, "k", bds%kopb)
                call hdf5_add_array(hg_cur_bp, "atoms", bds%opbat)
            endif
            call h5gclose_f(hg_cur_bp, eflag)
            
            ! Urey-Bradley stretching
            call h5gcreate_f(hg, "urey-bradley", hg_cur_bp, eflag)
            call hdf5_add_scalar(hg_cur_bp, "enabled", bds%use_urey)
            if(bds%use_urey) then
                call hdf5_add_scalar(hg_cur_bp, "cubic", bds%urey_cubic)
                call hdf5_add_scalar(hg_cur_bp, "quartic", bds%urey_quartic)
                call hdf5_add_array(hg_cur_bp, "k", bds%kurey)
                call hdf5_add_array(hg_cur_bp, "l0", bds%l0urey)
                call hdf5_add_array(hg_cur_bp, "atoms", bds%ureyat)
            end if
            call h5gclose_f(hg_cur_bp, eflag)
        end if
       
        call h5gclose_f(hg, eflag)

        out_fail = 0_ip
    end subroutine

    subroutine mmpol_init_from_hdf5(filename, namespace, s, out_fail)
        use hdf5
        use mod_adjacency_mat, only: build_conn_upto_n, yale_sparse
        use mod_io, only: ommp_message
        use mod_memory, only: ip, rp, mfree
        use mod_mmpol, only: mmpol_init, set_screening_parameters, &
                             mmpol_prepare, mmpol_init_nonbonded, mmpol_init_bonded
        use mod_constants, only: OMMP_FF_AMOEBA, OMMP_VERBOSE_LOW
        use mod_bonded, only: bond_init, angle_init, urey_init, strbnd_init, &
                              opb_init, pitors_init, torsion_init, tortor_init, &
                              strtor_init, angtor_init, tortor_newmap
        use mod_nonbonded, only: vdw_set_pair

        implicit none

        character(len=*), intent(in) :: namespace
        type(ommp_system), intent(inout) :: s
        character(len=*), intent(in) :: filename
        integer(ip), intent(out) :: out_fail
        
        integer(hid_t) :: hg_sysmodel, hg_res, hg_cur, hg_amoeba, &
                          hg_top, hg_cur_param, hg_cur_bp, cur_dst, cur_dsp
        integer(hsize_t), dimension(4) :: dims
        integer(hid_t) :: iof_hdf5 = 301
        integer(kind=4) :: eflag
        real(rp), dimension(:), allocatable :: l_mscale, l_pscale, l_dscale, &
                                               l_uscale, l_ipscale, l_vdwscale
        type(yale_sparse) :: conn_1, tmp_vdw_pair
        integer(ip) :: mm_atoms, pol_atoms
        logical :: amoeba, mutable_only, use_nonbonded, bp_exist

        ! For handling torsion maps
        integer(ip) :: i, j, ibeg, iend
        integer(ip), allocatable, dimension(:,:) :: tmp_shape
        real(rp), allocatable, dimension(:) :: tmp_ang1, tmp_ang2, tmp_v, &
                                               tmp_vdw_pair_e, tmp_vdw_pair_r

        ! Initialize interface
        call h5open_f(eflag)
        if(eflag /= 0) then
            call ommp_message("Unable to initialize HDF5 module. Failure in &
                              &h5open_f subroutine.", OMMP_VERBOSE_LOW)
            out_fail = -1_ip
            return
        end if
        
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, iof_hdf5, eflag)
        if( eflag /= 0) then 
            call ommp_message("Unable to open HDF5 file. Failure in &
                               &h5fopen_f subroutine.", OMMP_VERBOSE_LOW)
            out_fail = -1_ip
            return
        end if

        call hdf5_read_scalar(iof_hdf5, & 
                              namespace, 'mutable_only', &
                              mutable_only)
        if(mutable_only) then
            call ommp_message("Unable to intialize from the selected file/&
                              &namespace as it is marked mutable only", &
                              OMMP_VERBOSE_LOW)
            out_fail = -1_ip
            return
        end if
    
        call hdf5_read_scalar(iof_hdf5, &
                              namespace//'/topology', 'N-atoms', &
                              mm_atoms)
        call hdf5_read_scalar(iof_hdf5, &
                              namespace//'/electrostatics', 'N-pol-atoms', &
                              pol_atoms)
        call hdf5_read_scalar(iof_hdf5, &
                              namespace//'/electrostatics', 'amoeba', &
                              amoeba)
        
        if(amoeba) then
            call mmpol_init(s, 1_ip, mm_atoms, pol_atoms)
        else
            call mmpol_init(s, 0_ip, mm_atoms, pol_atoms)
        end if
        
        ! Connectivity 
        call hdf5_read_array(iof_hdf5, &
                             namespace//'/topology/connectivity/ADJ1-RowIdx', &
                             conn_1%ri)
        call hdf5_read_array(iof_hdf5, &
                             namespace//'/topology/connectivity/ADJ1-ColIdx', &
                             conn_1%ci)
        conn_1%n = size(conn_1%ri) - 1
        call build_conn_upto_n(conn_1, 4, s%top%conn, .false.)

        ! AMOEBA
        if(amoeba) then
            call hdf5_read_array(iof_hdf5, &
                                  namespace//'/electrostatics/polarization_group_id', &
                                  s%eel%mmat_polgrp)
        end if
        
        ! Bonded Parameters
        call hdf5_group_exists(iof_hdf5, namespace//'/bonded', bp_exist)
        if(bp_exist) then
            call mmpol_init_bonded(s)
            ! Bond stretching
            call hdf5_read_scalar(iof_hdf5, &
                                  namespace//'/bonded/stretching', &
                                  'enabled', &
                                  s%bds%use_bond)
            if(s%bds%use_bond) then
                call bond_init(s%bds, hdf5_array_len(iof_hdf5, &
                                              namespace//'/bonded/stretching/atoms'))

                call hdf5_read_scalar(iof_hdf5, &
                                      namespace//'/bonded/stretching', &
                                      'quartic', &
                                      s%bds%bond_quartic)
                call hdf5_read_scalar(iof_hdf5, &
                                      namespace//'/bonded/stretching', &
                                      'cubic', &
                                      s%bds%bond_cubic)
                call hdf5_read_array(iof_hdf5, &
                                     namespace//'/bonded/stretching/k', &
                                     s%bds%kbond)
                call hdf5_read_array(iof_hdf5, &
                                     namespace//'/bonded/stretching/l0', &
                                     s%bds%l0bond)
                call hdf5_read_array(iof_hdf5, &
                                     namespace//'/bonded/stretching/atoms', &
                                     s%bds%bondat)
            end if
            
            ! Angle bending
            call hdf5_read_scalar(iof_hdf5, &
                                  namespace//'/bonded/bending', &
                                  'enabled', &
                                  s%bds%use_angle)
            if(s%bds%use_angle) then
                call angle_init(s%bds, hdf5_array_len(iof_hdf5, &
                                              namespace//'/bonded/bending/atoms'))
                call hdf5_read_scalar(iof_hdf5, &
                                      namespace//'/bonded/bending', &
                                      "cubic", s%bds%angle_cubic)
                call hdf5_read_scalar(iof_hdf5, &
                                      namespace//'/bonded/bending', &
                                      "quartic", s%bds%angle_quartic)
                call hdf5_read_scalar(iof_hdf5, &
                                      namespace//'/bonded/bending', &
                                      "pentic", s%bds%angle_pentic)
                call hdf5_read_scalar(iof_hdf5, &
                                      namespace//'/bonded/bending', &
                                      "sextic", s%bds%angle_sextic)
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/bending/k", &
                                     s%bds%kangle)
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/bending/ang0", &
                                     s%bds%eqangle)
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/bending/atoms", &
                                     s%bds%angleat)
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/bending/type", &
                                     s%bds%anglety)
            end if
            
            ! Dihedral torsion
            call hdf5_read_scalar(iof_hdf5, &
                                  namespace//'/bonded/torsion', &
                                  'enabled', &
                                  s%bds%use_torsion)
            if(s%bds%use_torsion) then
                call torsion_init(s%bds, hdf5_array_len(iof_hdf5, &
                                              namespace//'/bonded/torsion/atoms'))
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/torsion/amplitudes", &
                                     s%bds%torsamp)
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/torsion/phase", &
                                     s%bds%torsphase)
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/torsion/atoms", &
                                     s%bds%torsionat)
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/torsion/period", &
                                     s%bds%torsn)
            end if
            
            ! Stretching-bending coupling
            call hdf5_read_scalar(iof_hdf5, &
                                  namespace//'/bonded/stretching-bending', &
                                  'enabled', &
                                  s%bds%use_strbnd)
            if(s%bds%use_strbnd) then
                call strbnd_init(s%bds, hdf5_array_len(iof_hdf5, &
                                              namespace//'/bonded/stretching-bending/atoms'))
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/stretching-bending/k1", &
                                     s%bds%strbndk1)
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/stretching-bending/k2", &
                                     s%bds%strbndk2)
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/stretching-bending/l1_0", &
                                     s%bds%strbndl10)
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/stretching-bending/l2_0", &
                                     s%bds%strbndl20)
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/stretching-bending/ang0", &
                                     s%bds%strbndthet0)
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/stretching-bending/atoms", &
                                     s%bds%strbndat)
            end if
            
            ! Stretching-torsion coupling
            call hdf5_read_scalar(iof_hdf5, &
                                  namespace//'/bonded/stretching-torsion', &
                                  'enabled', &
                                  s%bds%use_strtor)
            if(s%bds%use_strtor) then
                call strtor_init(s%bds, hdf5_array_len(iof_hdf5, &
                                              namespace//'/bonded/stretching-torsion/atoms'))
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/stretching-torsion/k", &
                                     s%bds%strtork)
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/stretching-torsion/bonds_idx", &
                                     s%bds%strtor_b)
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/stretching-torsion/torsion_idx", &
                                     s%bds%strtor_t)
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/stretching-torsion/atoms", &
                                     s%bds%strtorat)
            end if

            ! Bending-torsion coupling 
            call hdf5_read_scalar(iof_hdf5, &
                                  namespace//'/bonded/bending-torsion', &
                                  'enabled', &
                                  s%bds%use_angtor)
            if(s%bds%use_angtor) then
                call angtor_init(s%bds, hdf5_array_len(iof_hdf5, &
                                              namespace//'/bonded/bending-torsion/atoms'))
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/bending-torsion/k", &
                                     s%bds%angtork)
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/bending-torsion/angles_idx", &
                                     s%bds%angtor_a)
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/bending-torsion/torsion_idx", &
                                     s%bds%angtor_t)
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/bending-torsion/atoms", &
                                     s%bds%angtorat)
            end if
            
            ! Torsion-torsion coupling
            call hdf5_read_scalar(iof_hdf5, &
                                  namespace//'/bonded/torsion-torsion', &
                                  'enabled', &
                                  s%bds%use_tortor)
            if(s%bds%use_tortor) then
                call tortor_init(s%bds, hdf5_array_len(iof_hdf5, &
                                              namespace//'/bonded/torsion-torsion/atoms'))
                call hdf5_read_array(iof_hdf5, &
                                    namespace//"/bonded/torsion-torsion/atoms", &
                                    s%bds%tortorat)
                call hdf5_read_array(iof_hdf5, &
                                    namespace//"/bonded/torsion-torsion/map_id", &
                                    s%bds%tortorprm)
                call hdf5_read_array(iof_hdf5, &
                                    namespace//"/bonded/torsion-torsion/maps_ang1", &
                                    tmp_ang1)
                call hdf5_read_array(iof_hdf5, &
                                    namespace//"/bonded/torsion-torsion/maps_ang2", &
                                    tmp_ang2)
                call hdf5_read_array(iof_hdf5, &
                                    namespace//"/bonded/torsion-torsion/maps_pot", &
                                    tmp_v)
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/torsion-torsion/maps_shapes", &
                                     tmp_shape)
                
                iend = 0
                do i=1, size(tmp_shape, 2)
                    ibeg = iend + 1 
                    iend = ibeg + tmp_shape(1,i) * tmp_shape(2,i) - 1

                    call tortor_newmap(s%bds, tmp_shape(1,i), &
                                       tmp_shape(2,i), &
                                       tmp_ang1(ibeg:iend), &
                                       tmp_ang2(ibeg:iend), &
                                       tmp_v(ibeg:iend))
                end do

                call mfree('mmpol_init_from_hdf5 [tmp_ang1]', tmp_ang1)
                call mfree('mmpol_init_from_hdf5 [tmp_ang2]', tmp_ang2)
                call mfree('mmpol_init_from_hdf5 [tmp_v]', tmp_v)
                call mfree('mmpol_init_from_hdf5 [tmp_shape]', tmp_shape)
            end if
            
            ! Pi-torsion
            call hdf5_read_scalar(iof_hdf5, &
                                  namespace//'/bonded/pi-torsion', &
                                  'enabled', &
                                  s%bds%use_pitors)
            if(s%bds%use_pitors) then
                call pitors_init(s%bds, hdf5_array_len(iof_hdf5, &
                                              namespace//'/bonded/pi-torsion/atoms'))
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/pi-torsion/atoms", &
                                     s%bds%pitorsat)
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/pi-torsion/k", &
                                     s%bds%kpitors)
            end if

            ! Out-of-plane bending
            call hdf5_read_scalar(iof_hdf5, &
                                  namespace//'/bonded/out-of-plane-bending', &
                                  'enabled', &
                                  s%bds%use_opb)
            if(s%bds%use_opb) then
                call opb_init(s%bds, hdf5_array_len(iof_hdf5, &
                                              namespace//'/bonded/out-of-plane-bending/atoms'), 'allinger')
                call hdf5_read_scalar(iof_hdf5, &
                                      namespace//'/bonded/out-of-plane-bending', &
                                      "cubic", s%bds%opb_cubic)
                call hdf5_read_scalar(iof_hdf5, &
                                      namespace//'/bonded/out-of-plane-bending', &
                                      "quartic", s%bds%opb_quartic)
                call hdf5_read_scalar(iof_hdf5, &
                                      namespace//'/bonded/out-of-plane-bending', &
                                      "pentic", s%bds%opb_pentic)
                call hdf5_read_scalar(iof_hdf5, &
                                      namespace//'/bonded/out-of-plane-bending', &
                                      "sextic", s%bds%opb_sextic)
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/out-of-plane-bending/k", &
                                     s%bds%kopb)
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/out-of-plane-bending/atoms", &
                                     s%bds%opbat)
            endif
            
            ! Urey-Bradley stretching
            call hdf5_read_scalar(iof_hdf5, &
                                  namespace//'/bonded/urey-bradley', &
                                  'enabled', &
                                  s%bds%use_urey)
            if(s%bds%use_urey) then
                call urey_init(s%bds, hdf5_array_len(iof_hdf5, &
                                              namespace//'/bonded/urey-bradley/atoms'))
                call hdf5_read_scalar(iof_hdf5, &
                                      namespace//'/bonded/urey-bradley', &
                                      "cubic", s%bds%urey_cubic)
                call hdf5_read_scalar(iof_hdf5, &
                                      namespace//'/bonded/urey-bradley', &
                                      "quartic", s%bds%urey_quartic)
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/urey-bradley/k", &
                                     s%bds%kurey)
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/urey-bradley/l0", &
                                     s%bds%l0urey)
                call hdf5_read_array(iof_hdf5, &
                                     namespace//"/bonded/urey-bradley/atoms", &
                                     s%bds%ureyat)
            end if
        end if
        
        call hdf5_group_exists(iof_hdf5, namespace//"/nonbonded", use_nonbonded)
        if(use_nonbonded) then
            call mmpol_init_nonbonded(s)
            call vdw_init(s%vdw, s%top, "buffered-14-7", "cubic-mean", "diameter", "r-min", &
                          "hhg")
            call hdf5_read_array(iof_hdf5, & 
                                 namespace//"/nonbonded/screening", &
                                 l_vdwscale)
            
            s%vdw%vdw_screening = l_vdwscale
            call mfree('mmpol_init_from_hdf5 [l_vdwscale]', l_vdwscale)

            call hdf5_read_array(iof_hdf5, & 
                                 namespace//"/nonbonded/radius", &
                                 s%vdw%vdw_r)
            call hdf5_read_array(iof_hdf5, & 
                                 namespace//"/nonbonded/energy", &
                                 s%vdw%vdw_e)
            call hdf5_read_array(iof_hdf5, & 
                                 namespace//"/nonbonded/scale_factor", &
                                 s%vdw%vdw_f)
            call hdf5_read_array(iof_hdf5, & 
                                 namespace//"/nonbonded/pair_row_idx", &
                                 tmp_vdw_pair%ri)
            call hdf5_read_array(iof_hdf5, & 
                                 namespace//"/nonbonded/pair_col_idx", &
                                 tmp_vdw_pair%ci)
            call hdf5_read_array(iof_hdf5, & 
                                 namespace//"/nonbonded/pair_radius", &
                                 tmp_vdw_pair_r)
            call hdf5_read_array(iof_hdf5, & 
                                 namespace//"/nonbonded/pair_energy", &
                                 tmp_vdw_pair_e)
        end if

        do i=1, s%top%mm_atoms
            do j=tmp_vdw_pair%ri(i), tmp_vdw_pair%ri(i+1)-1
               call vdw_set_pair(s%vdw, i, tmp_vdw_pair%ci(j), tmp_vdw_pair_r(j), &
                                 tmp_vdw_pair_e(j)) 
            end do
        end do
        call mfree('mmpol_init_hdf5', tmp_vdw_pair_e)
        call mfree('mmpol_init_hdf5', tmp_vdw_pair_r)
        
        call hdf5_read_array(iof_hdf5, &
                             namespace//"/electrostatics/fixed_multipoles_unrotated", s%eel%q)
        if(amoeba) then
            call hdf5_read_array(iof_hdf5, &
                                namespace//"/electrostatics/fixed_mmpoles_rot_Z", &
                                s%eel%iz)
            call hdf5_read_array(iof_hdf5, &
                                namespace//"/electrostatics/fixed_mmpoles_rot_X", &
                                s%eel%ix)
            call hdf5_read_array(iof_hdf5, &
                                namespace//"/electrostatics/fixed_mmpoles_rot_Y", &
                                s%eel%iy)
            call hdf5_read_array(iof_hdf5, &
                                namespace//"/electrostatics/fixed_mmpoles_rot_CONV", &
                                s%eel%mol_frame)
        end if
        call hdf5_read_array(iof_hdf5, &
                             namespace//"/electrostatics/fixed_fixed_scale_f", l_mscale)
        call hdf5_read_array(iof_hdf5, &
                             namespace//"/electrostatics/fixed_ipd_scale_f", l_pscale)
        call hdf5_read_array(iof_hdf5, &
                             namespace//"/electrostatics/ipd_ipd_scale_f", l_uscale)
        if(amoeba) then
            call hdf5_read_array(iof_hdf5, &
                                 namespace//"/electrostatics/fixed_direct_ipd_scale_f", &
                                 l_dscale)
            call hdf5_read_array(iof_hdf5, &
                                 namespace//"/electrostatics/fixed_intragroup_ipd_scale_f", &
                                 l_ipscale)
        end if
        call set_screening_parameters(s%eel, l_mscale, l_pscale, l_dscale, l_uscale, &
                                      l_ipscale)
        call mfree('mmpol_init_from_hdf5 [l_mscale]', l_mscale)
        call mfree('mmpol_init_from_hdf5 [l_pscale]', l_pscale)
        call mfree('mmpol_init_from_hdf5 [l_dscale]', l_dscale)
        call mfree('mmpol_init_from_hdf5 [l_uscale]', l_uscale)
        call mfree('mmpol_init_from_hdf5 [l_ipscale]', l_ipscale)
        
        call hdf5_read_array(iof_hdf5, &
                             namespace//"/electrostatics/polarizable_atoms_idx", s%eel%polar_mm)
        call hdf5_read_array(iof_hdf5, &
                             namespace//"/electrostatics/polarizabilities", s%eel%pol)
        
        call hdf5_read_array(iof_hdf5, namespace//"/topology/Atoms-Coordinates", s%top%cmm)

        call h5fclose_f(iof_hdf5, eflag)
        if( eflag /= 0) then 
            call ommp_message("Error while closing HDF5 file. Failure in &
                               &h5fclose_f subroutine.", OMMP_VERBOSE_LOW)
            out_fail = -1_ip
            return
        end if
        
        call mmpol_prepare(s)

        out_fail = 0_ip
    end subroutine mmpol_init_from_hdf5

end module mod_iohdf5
#endif
