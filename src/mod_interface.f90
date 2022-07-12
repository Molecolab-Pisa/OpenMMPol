! Wrapper function for open-mmpol library
module mod_interface
    use iso_c_binding
    use mod_mmpol
    use mod_memory, only: ip, rp

    implicit none
    private

    public :: get_mm_atoms, get_pol_atoms, get_n_ipd, get_ld_cart
    public :: is_amoeba, get_cmm, get_cpol, get_q, get_ipd
    public :: get_polar_mm
    public :: mmpol_init_mmp, do_mm, do_qmmm, get_energy, &
              write_hdf5

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
        
        function get_polar_mm() bind(c, name='get_polar_mm')
            type(c_ptr) :: get_polar_mm

            get_polar_mm = c_loc(polar_mm)
        end function get_polar_mm

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
        
        subroutine set_verbose(verb) bind(c, name='set_verbose')
            implicit none 

            integer(ip), intent(in), value :: verb

            call set_verbosity(verb)
        end subroutine set_verbose

        subroutine print_summary() bind(c, name='print_summary')
            use mod_io, only: mmpol_print_summary

            implicit none
            
            call mmpol_print_summary()

        end subroutine print_summary
        
        subroutine print_summary_to_file(filename) bind(c, name='print_summary_to_file')
            use mod_io, only: mmpol_print_summary

            implicit none
            
            character(kind=c_char), intent(in) :: filename(120)
            character(len=120) :: output_file
            
            call c2f_string(filename, output_file)
            call mmpol_print_summary(output_file)

        end subroutine print_summary_to_file

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
        
        subroutine mmpol_init_xyz(xyzfile, prmfile) bind(c, name='mmpol_init_xyz')
            use mod_inputloader, only : mmpol_init_from_xyz
            
            implicit none
            
            character(kind=c_char), intent(in) :: xyzfile(120), prmfile(120)
            character(len=120) :: xyz_file, prm_file

            call c2f_string(prmfile, prm_file)
            call c2f_string(xyzfile, xyz_file)
            call mmpol_init_from_xyz(xyz_file, prm_file)
        end subroutine 
        
        subroutine mmpol_init_mmp(filename) bind(c, name='mmpol_init_mmp')
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
            use mod_io, only: print_matrix
            use mod_electrostatics, only: prepare_M2M, prepare_M2D
            
            implicit none
  
            call prepare_M2M()
            call prepare_M2D()

        end subroutine 

        subroutine do_qmmm(ef_qmd, solver) bind(c, name='do_qmmm')
            ! Compute polarization and QM-MM energy
            use mod_polarization, only: polarization
            use mod_mmpol, only: pol_atoms, n_ipd
            use mod_electrostatics, only: e_m2d, prepare_M2D

            implicit none
            
            real(kind=rp), intent(in) :: ef_qmd(3, pol_atoms, n_ipd)
            integer(ip), intent(in), value :: solver

            call prepare_M2D()
            call polarization(e_m2d+ef_qmd, ipd, solver)
        end subroutine

        subroutine get_energy(EMM, EPol) bind(c, name='get_energy')
            ! Get the energy
            use mod_electrostatics, only: energy_MM_MM, energy_MM_pol

            implicit none
            real(kind=rp), intent(out) :: EMM, EPol

            emm = 0.0_rp
            epol = 0.0_rp
            
            call energy_MM_MM(emm)
            call energy_MM_pol(epol)

        end subroutine

        subroutine get_urey_energy(eub) bind(c, name='get_urey_energy')
            
            use mod_bonded, only: urey_potential

            implicit none
            real(rp), intent(inout) :: eub

            eub = 0.0_rp
            call urey_potential(eub)
        end subroutine
        
        subroutine get_angle_energy(eang) bind(c, name='get_angle_energy')
            
            use mod_bonded, only: angle_potential

            implicit none
            real(rp), intent(inout) :: eang

            eang = 0.0_rp
            call angle_potential(eang)
        end subroutine
        
        subroutine get_bond_energy(ebnd) bind(c, name='get_bond_energy')
            
            use mod_bonded, only: bond_potential

            implicit none
            real(rp), intent(inout) :: ebnd

            ebnd = 0.0_rp
            call bond_potential(ebnd)
        end subroutine
        
        subroutine get_vdw_energy(evdw) bind(c, name='get_vdw_energy')
            
            use mod_nonbonded, only: vdw_potential

            implicit none
            real(rp), intent(inout) :: evdw

            evdw = 0.0_rp
            call vdw_potential(evdw)
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

