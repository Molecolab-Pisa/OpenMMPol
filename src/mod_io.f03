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
                             amoeba, pol, conn, ff_rules, ff_type, &
                             ix, iy, iz, mol_frame, ip11, maxpgp

        implicit none

        character(len=*), intent(in) :: filename
        integer(ip), intent(out) :: out_fail
        
        integer(hid_t) :: hg_sysfund, hg_sysder, hg_res, hg_cur, hg_amoeba, &
                          cur_dst, cur_dsp
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

        ! Attributes
        dims = (/1, 0, 0, 0/)
        call H5Screate_f(H5S_SCALAR_F, cur_dsp, eflag)
        
        call H5Acreate_f(hg_sysfund, &
                         "numer of MM atoms", &
                         H5T_IP, &
                         cur_dsp, cur_dst, eflag)
        call H5Awrite_f(cur_dst, H5T_IP, mm_atoms, dims(:1), eflag)
        
        call H5Acreate_f(hg_sysfund, &
                         "numer of POL atoms", &
                         H5T_IP, &
                         cur_dsp, cur_dst, eflag)
        call H5Awrite_f(cur_dst, H5T_IP, pol_atoms, dims(:1), eflag)

        call H5Acreate_f(hg_sysfund, &
                         "FF type", &
                         H5T_IP, &
                         cur_dsp, cur_dst, eflag)
        call H5Awrite_f(cur_dst, H5T_IP, ff_type, dims(:1), eflag)
        
        call H5Acreate_f(hg_sysfund, &
                         "FF rules", &
                         H5T_IP, &
                         cur_dsp, cur_dst, eflag)
        call H5Awrite_f(cur_dst, H5T_IP, ff_rules, dims(:1), eflag)

        ! Dataset
        ! coordinates
        dims = (/3, mm_atoms, 0, 0/)
        call H5Screate_simple_f(2, dims(:2), cur_dsp, eflag)
        call H5Dcreate_f(hg_sysfund, &
                         "Coordinates of MM sites", &
                         H5T_RP, &
                         cur_dsp, cur_dst, eflag)
        call H5Dwrite_f(cur_dst, H5T_RP, cmm, dims(:2), eflag)

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

        ! TODO
        ! Connectivity of the environment is saved as an adjacency matrix; 
        ! Since such a matrix is sparse and boolean it can be represented in
        ! a very efficient way using Yale format (sometimes called compressed
        ! sparse row) omitting the V vector. In this way the adj. matrix is 
        ! represented by adj_m_ci(n_bond) and adj_m_ri(mm_atoms+1).
        
        ! For now, just save n1m and i1m, all closed in a subfolder
        call h5gcreate_f(hg_sysfund, "topology", hg_cur, eflag)
        
        dims = (/mm_atoms+1, 0, 0, 0/)
        call h5screate_simple_f(1, dims(:1), cur_dsp, eflag)
        call h5dcreate_f(hg_cur, &
                         "Adjacency matrix (Yale format) RowIdx", &
                         H5T_IP, &
                         cur_dsp, cur_dst, eflag)
        call h5dwrite_f(cur_dst, H5T_IP, conn(1)%ri, dims(:1), eflag)
        
        dims = (/size(conn(1)%ci), 0, 0, 0/)
        call h5screate_simple_f(1, dims(:1), cur_dsp, eflag)
        call h5dcreate_f(hg_cur, &
                         "Adjacency matrix (Yale format) ColumnIdx", &
                         H5T_IP, &
                         cur_dsp, cur_dst, eflag)
        call h5dwrite_f(cur_dst, H5T_IP, conn(1)%ci, dims(:1), eflag)
        
        call h5gclose_f(hg_cur, eflag)

        if(amoeba) then
            call h5gcreate_f(hg_sysfund, "amoeba", hg_amoeba, eflag)
            ! Rotation convenction
            call h5gcreate_f(hg_amoeba, "rotation", hg_cur, eflag)
            
            dims = (/mm_atoms, 0, 0, 0/)
            call h5screate_simple_f(1, dims(:1), cur_dsp, eflag)
            call h5dcreate_f(hg_cur, &
                             "X-axys index", &
                             H5T_IP, &
                             cur_dsp, cur_dst, eflag)
            call h5dwrite_f(cur_dst, H5T_IP, ix, dims(:1), eflag)
            
            call h5dcreate_f(hg_cur, &
                             "Y-axys index", &
                             H5T_IP, &
                             cur_dsp, cur_dst, eflag)
            call h5dwrite_f(cur_dst, H5T_IP, iy, dims(:1), eflag)
            
            call h5dcreate_f(hg_cur, &
                             "Z-axys index", &
                             H5T_IP, &
                             cur_dsp, cur_dst, eflag)
            call h5dwrite_f(cur_dst, H5T_IP, iz, dims(:1), eflag)
            
            call h5dcreate_f(hg_cur, &
                             "Molecular frame def", &
                             H5T_IP, &
                             cur_dsp, cur_dst, eflag)
            call h5dwrite_f(cur_dst, H5T_IP, mol_frame, dims(:1), eflag)
            
            ! Group connectivity 
            call h5gcreate_f(hg_amoeba, "polarization groups connectivity", &
                             hg_cur, eflag)
            
            dims = (/maxpgp, mm_atoms, 0, 0/)
            call h5screate_simple_f(2, dims(:2), cur_dsp, eflag)
            call h5dcreate_f(hg_cur, &
                             "Adjacency indices", &
                             H5T_IP, &
                             cur_dsp, cur_dst, eflag)
            call h5dwrite_f(cur_dst, H5T_IP, ip11, dims(:2), eflag)
            
            call h5gclose_f(hg_cur, eflag)
            
            call h5gclose_f(hg_amoeba, eflag)
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

        integer(ip) :: i, j, grp, igrp, lst(1000), ilst
        character(len=120) :: str

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
            
            do j=1, 4
                if(j == 4 .and. .not. amoeba) cycle
                
                write(str, "('1-', I1, ' neighors:')") j+1
                call print_int_vec(str, &
                                   size(conn(j)%ci), &
                                   conn(j)%ri(i), &
                                   conn(j)%ri(i+1)-1, & 
                                   conn(j)%ci)
            end do
            
            if(amoeba) then 
                
                call print_int_vec('1-1 polarization neighors:', &
                                   mm_atoms, &
                                   polgrp_mmat%ri(mmat_polgrp(i)), &
                                   polgrp_mmat%ri(mmat_polgrp(i)+1)-1, &
                                   polgrp_mmat%ci)
               
                do j=1, 3
                    ilst = 1
                    do igrp=pg_conn(j)%ri(mmat_polgrp(i)), &
                            pg_conn(j)%ri(mmat_polgrp(i)+1)-1
                        grp = pg_conn(j)%ci(igrp)
                        lst(ilst:ilst+polgrp_mmat%ri(grp+1)-polgrp_mmat%ri(grp)-1) = &
                        polgrp_mmat%ci(polgrp_mmat%ri(grp):polgrp_mmat%ri(grp+1)-1)
                        ilst = ilst+polgrp_mmat%ri(grp+1)-polgrp_mmat%ri(grp)
                    end do
                    
                    write(str, "('1-', I1, ' polarization neighors:')") j+1
                    call print_int_vec(str, & 
                                       ilst-1,0,0,lst)
                end do
            end if
        end do

    end subroutine mmpol_print_summary

end module mod_io
