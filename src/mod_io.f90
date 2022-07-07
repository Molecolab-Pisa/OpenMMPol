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
                             amoeba, pol, conn, ff_type, &
                             ix, iy, iz, mol_frame, mmat_polgrp

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

        ! Connectivity of the environment is saved as an adjacency matrix; 
        ! Since such a matrix is sparse and boolean it can be represented in
        ! a very efficient way using Yale format (sometimes called compressed
        ! sparse row) omitting the V vector. In this way the adj. matrix is 
        ! represented by adj_m_ci(n_bond) and adj_m_ri(mm_atoms+1).
        
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
            
            dims = (/mm_atoms, 0, 0, 0/)
            call h5screate_simple_f(1, dims(:1), cur_dsp, eflag)
            call h5dcreate_f(hg_amoeba, &
                             "Polarization groups index", &
                             H5T_IP, &
                             cur_dsp, cur_dst, eflag)
            call h5dwrite_f(cur_dst, H5T_IP, mmat_polgrp, dims(:1), eflag)
            
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
          write(out_unit,*)
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
        
        if(sort) then
            call sort_ivec(vec(ib:ie), sorted_vec)
            write(out_unit,'(t5, 10i8)') sorted_vec
        else
            write(out_unit,'(t5, 10i8)') vec(ib:ie)
        end if

        return

    end subroutine print_int_vec

end module mod_io
