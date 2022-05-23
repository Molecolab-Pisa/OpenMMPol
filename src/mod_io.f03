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
    public :: mmpol_print_summary

#ifdef USE_HDF5
    integer(hid_t) :: iof_hdf5_out = 101
    public :: mmpol_save_as_hdf5
#endif
    
    contains
    
#ifdef USE_HDF5
#define H5T_RP H5T_NATIVE_DOUBLE
#ifdef USE_I8
#define H5T_IP H5T_STD_I64LE
#else
#define H5T_IP H5T_STD_I32LE
#endif
    
    ! Subroutines dedicated to HDF5 I/O
    subroutine mmpol_save_as_hdf5(filename, out_fail)
        use hdf5
        use mod_memory, only: ip
        use mod_mmpol, only: mm_atoms, pol_atoms, cmm, polar_mm, ld_cart, q, &
                             amoeba, pol

        implicit none

        character(len=*), intent(in) :: filename
        integer(ip), intent(out) :: out_fail
        
        integer(hid_t) :: hg_sysfund, hg_sysder, hg_res, cur_dst, cur_dsp
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
        
        call h5gcreate_f(iof_hdf5_out, "system fundamental", hg_sysfund, eflag)
        if( eflag /= 0) then 
            write(iof_mmpol, *) "Error while creating group 'system fundamental.'&
                                &Failure in h5gcreate_f subroutine."
            out_fail = -1_ip
            return
        end if
        
        ! Write H5G:system fundamental
        ! Those are the minimal quantities needed to start a calculation.
        ! It should contain only the information needed by the program to 
        ! correctly describe the system.
       
        ! coordinates
        dims = (/3, mm_atoms, 0, 0/)
        call h5screate_simple_f(2, dims(:2), cur_dsp, eflag)
        call h5dcreate_f(hg_sysfund, &
                         "Coordinates of MM sites", &
                         H5T_RP, &
                         cur_dsp, cur_dst, eflag)
        call h5dwrite_f(cur_dst, H5T_RP, cmm, dims(:2), eflag)

        ! list of polarizable atoms
        dims = (/pol_atoms, 0, 0, 0/)
        call h5screate_simple_f(1, dims(:1), cur_dsp, eflag)
        call h5dcreate_f(hg_sysfund, &
                         "Index of polarizable MM sites (1-based)", &
                         H5T_IP, &
                         cur_dsp, cur_dst, eflag)
        call h5dwrite_f(cur_dst, H5T_IP, polar_mm, dims(:1), eflag)
        
        ! q
        dims = (/ld_cart, mm_atoms, 0, 0/)
        call h5screate_simple_f(2, dims(:2), cur_dsp, eflag)
        call h5dcreate_f(hg_sysfund, &
                         "Multiopolar distributions on MM atoms", &
                         H5T_RP, &
                         cur_dsp, cur_dst, eflag)
        call h5dwrite_f(cur_dst, H5T_RP, q, dims(:2), eflag)
        
        ! polarizabilities 
        dims = (/pol_atoms, 0, 0, 0/)
        call h5screate_simple_f(1, dims(:1), cur_dsp, eflag)
        call h5dcreate_f(hg_sysfund, &
                         "Polarizabilities of POL atoms", &
                         H5T_RP, &
                         cur_dsp, cur_dst, eflag)
        call h5dwrite_f(cur_dst, H5T_RP, pol, dims(:1), eflag)

        ! Which is the best way to save the connectivity??

        if(amoeba) then
            ! Rotation convenction

            ! Group connectivity 
        endif

        call h5gclose_f(hg_sysfund, eflag)
        if( eflag /= 0) then 
            write(iof_mmpol, *) "Error while closing group 'system fundamental.'&
                                &Failure in h5gclose_f subroutine."
            out_fail = -1_ip
            return
        end if

        
        call h5gcreate_f(iof_hdf5_out, "system derived", hg_sysder, eflag)
        if( eflag /= 0) then 
            write(iof_mmpol, *) "Error while creating group 'system derived'.&
                                &Failure in h5gcreate_f subroutine."
            out_fail = -1_ip
            return
        end if
        
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
    
subroutine mmpol_print_summary()
        !! Prints a complete summary of all the quantities stored 
        !! in the MMPol module

        use mod_mmpol ! TODO only 

        implicit none

        integer(ip) :: i

        call print_matrix(.true.,'coordinates:', &
                          3,mm_atoms,3,mm_atoms,cmm)
        if (amoeba) then
            call print_matrix(.true.,'multipoles - non rotated:', &
                              ld_cart,mm_atoms,ld_cart,mm_atoms,q0)
        end if
        call print_matrix(.true.,'multipoles :', &
                          ld_cart,mm_atoms,ld_cart,mm_atoms,q)
        call print_matrix(.true.,'coordinates of polarizable atoms:', &
                          3,pol_atoms,3,pol_atoms,cpol)
        call print_matrix(.false.,'polarizabilities:', &
                          mm_atoms,1,mm_atoms,1,pol)
        call print_matrix(.false.,'thole factors:', &
                          mm_atoms,1,mm_atoms,1,thole)
        call print_int_vec('mm_polar list:', &
                           mm_atoms,0,0,mm_polar)
        call print_int_vec('polar_mm list:', &
                           pol_atoms,0,0,polar_mm)

        ! write the connectivity information for each atom:
  1000  format(t3,'connectivity information for the ',i8,'-th atom:')
    
        do i = 1, mm_atoms
            write(iof_mmpol, 1000) i
            
            call print_int_vec('1-2 neighors:',n12(i),0,0,i12(:,i))
            call print_int_vec('1-3 neighors:',n13(i),0,0,i13(:,i))
            call print_int_vec('1-4 neighors:',n14(i),0,0,i14(:,i))

            if(amoeba) then 
                call print_int_vec('1-5 neighors:',n15(i),0,0,i15(:,i))
                call print_int_vec('1-1 polarization neighors:', &
                                   np11(i),0,0,ip11(:,i))
                call print_int_vec('1-2 polarization neighors:', &
                                   np12(i),0,0,ip12(:,i))
                call print_int_vec('1-3 polarization neighors:', &
                                   np13(i),0,0,ip13(:,i))
                call print_int_vec('1-4 polarization neighors:', &
                                   np14(i),0,0,ip14(:,i))
            end if
        end do

    end subroutine mmpol_print_summary

end module mod_io
