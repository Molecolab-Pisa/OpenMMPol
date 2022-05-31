! Wrapper function for open-mmpol library
module mod_interface
    use iso_c_binding
    use mod_mmpol
    use mod_memory, only: ip, rp

    implicit none
    private

    public :: get_mm_atoms, get_pol_atoms, get_n_ipd, get_ld_cart
    public :: is_amoeba, get_cmm, get_cpol, get_q, get_ipd
    public :: w_mmpol_init, do_mm, do_qmmm, restart, get_energy

    contains
        
        function get_cmm() bind(c, name='get_cmm')
            type(c_ptr) :: get_cmm

            get_cmm = c_loc(cmm)
        end function get_cmm

        function get_cpol() bind(c, name='get_cpol')
            type(c_ptr) :: get_cpol

            get_cpol = c_loc(cpol)
        end function get_cpol

        function get_q() bind(c, name='get_q')
            type(c_ptr) :: get_q

            get_q = c_loc(q)
        end function get_q

        function get_ipd() bind(c, name='get_ipd')
            type(c_ptr) :: get_ipd

            get_ipd = c_loc(ipd)
        end function get_ipd

        function get_mm_atoms() bind(c, name='get_mm_atoms')
            implicit none

            integer(ip) :: get_mm_atoms

            get_mm_atoms = mm_atoms
        end function get_mm_atoms
        
        function get_pol_atoms() bind(c, name='get_pol_atoms')
            implicit none

            integer(ip) :: get_pol_atoms

            get_pol_atoms = pol_atoms
        end function get_pol_atoms

        function get_n_ipd() bind(c, name='get_n_ipd')
            implicit none

            integer(ip) :: get_n_ipd

            get_n_ipd = n_ipd
        end function get_n_ipd

        function get_ld_cart() bind(c, name='get_ld_cart')
            implicit none

            integer(ip) :: get_ld_cart

            get_ld_cart = ld_cart
        end function get_ld_cart

        function is_amoeba() bind(c, name='is_amoeba')
            implicit none

            logical(c_bool) :: is_amoeba

            is_amoeba = amoeba
        end function is_amoeba
        
        subroutine c2f_string(c_str, f_str)
            implicit none
            
            character(kind=c_char), intent(in) :: c_str(:)
            character(len=*), intent(out) :: f_str

            integer :: i 

            i = 1
            do while(c_str(i) /= c_null_char)
                f_str(i:i) = c_str(i)
                i = i + 1
            end do

            do i = i, len(f_str)
                f_str(i:i) = ' '
            end do

            f_str = trim(f_str)
        end subroutine c2f_string
        
        subroutine w_mmpol_init(filename) bind(c, name='w_mmpol_init')
            use mod_inputloader, only : mmpol_init_from_mmp
            implicit none
            character(kind=c_char), intent(in) :: filename(120)
            
            character(len=120) :: input_file

            call c2f_string(filename, input_file)
            call mmpol_init_from_mmp(input_file)
        end subroutine 

        subroutine do_mm() bind(c, name='do_mm') 
            ! Perform the MM-only calculations of pot and field
            ! and return EMM

            implicit none
  
            ! Initialize variables
            v_qq   = 0.0_rp
            ef_qd  = 0.0_rp
            dv_qq  = 0.0_rp
            def_qd = 0.0_rp
    
            call electrostatics(0_ip,0_ip,0_ip,v_qq,ef_qd,dv_qq,def_qd)
            call electrostatics(0_ip,1_ip,0_ip,v_qq,ef_qd,dv_qq,def_qd)

            !  call print_matrix(.true.,'VMM:',mm_atoms,ld_cart,mm_atoms,ld_cart,transpose(v_qq))               
            !  call print_matrix(.true.,'EMM 1:',pol_atoms,3,pol_atoms,3,transpose(ef_qd(:,:,1)))               
            !  call print_matrix(.true.,'EMM 2:',pol_atoms,3,pol_atoms,3,transpose(ef_qd(:,:,2)))               

        end subroutine 


        subroutine do_qmmm(v_qmm,ef_qmd,n1,n2,n3,n4) bind(c, name='do_qmmm')
            ! Compute polarization and QM-MM energy

            implicit none
            ! real(kind=rp), intent(in) :: v_qmm(ld_cart,mm_atoms), ef_qmd(3,pol_atoms,n_ipd)
            integer(kind=ip), intent(in), value :: n1,n2,n3,n4
            real(kind=rp), intent(in) :: v_qmm(n1,n2), ef_qmd(3,n3,n4)

            !call print_matrix(.true.,'VQM:',ld_cart,mm_atoms,ld_cart,mm_atoms,v_qmm)
            !call print_matrix(.true.,'EQM 1:',pol_atoms,3,pol_atoms,3,transpose(ef_qmd(:,:,1))) 
            !call print_matrix(.true.,'EQM 2:',pol_atoms,3,pol_atoms,3,transpose(ef_qmd(:,:,2)))
            ! Initialize variables
            ipd = 0.0_rp 

            ! Compute polarization
            nmax = 200
            !convergence = 5.0E-10
            call polarization(solver,ef_qd+ef_qmd,ipd) ! - jacobi interation solver
            !call polarization(1,ef_qd+ef_qmd,ipd) 

        end subroutine

        subroutine restart() bind(c, name='restart')
            implicit none

            v_qq   = 0.0_rp 
            ef_qd  = 0.0_rp
            dv_qq  = 0.0_rp
            def_qd = 0.0_rp
            ipd    = 0.0_rp

        end subroutine

        subroutine get_energy(EMM,EPol) bind(c, name='get_energy')
            ! Get the energy

            implicit none
            real(kind=rp), intent(out) :: EMM, EPol

            call energy(0_ip,EMM)
            call energy(1_ip,EPol)

        end subroutine

#ifdef USE_HDF5
        function write_hdf5(filename) bind(c, name='write_hdf5')
            !! This function is an interface for saving an HDF5 file 
            !! with all the data contained in mmpol module using
            !! [[mod_io:mmpol_save_as_hdf5]]
            use mod_io, only: mmpol_save_as_hdf5

            implicit none
            
            character(kind=c_char), intent(in) :: filename(120)
            
            character(len=120) :: hdf5out
            integer(ip) :: write_hdf5

            call c2f_string(filename, hdf5out)
            call mmpol_save_as_hdf5(hdf5out, write_hdf5)
            
        end function write_hdf5
#endif

end module mod_interface

