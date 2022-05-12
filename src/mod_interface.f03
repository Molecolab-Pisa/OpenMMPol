! Wrapper function for open-mmpol library
module mod_interface
    use iso_c_binding
    use mmpol
    use mod_memory, only: ip, rp

    implicit none
    private

    public :: w_mmpol_init, do_mm, do_qmmm, restart, get_energy

    contains
        subroutine c2f_string(c_str, f_str)
            implicit none
            
            character(kind=c_char), intent(in) :: c_str(:)
            character(len=*), intent(out) :: f_str

            integer :: i = 1

            do while(c_str(i) /= c_null_char)
                f_str(i:i) = c_str(i)
            end do

            f_str = trim(f_str)
        end subroutine c2f_string
        
        subroutine w_mmpol_init(filename) bind(c, name='w_mmpol_init')
            implicit none
            character(kind=c_char), intent(in) :: filename(120)
            character(len=120) :: BsName

            call c2f_string(filename, input_file)
            BsName = input_file(1:scan(input_file, '.', back=.true.))
            scratch_file = trim(BsName)//'rwf'

            call mmpol_init()
        end subroutine 

        subroutine do_mm() bind(c, name='do_mm') 
            ! Perform the MM-only calculations of pot and field
            ! and return EMM

            implicit none
  
            ! Initialize variables
            v_qq   = Zero
            ef_qd  = Zero
            dv_qq  = Zero
            def_qd = Zero
    
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
            integer(kind=ip), intent(in) :: n1,n2,n3,n4
            real(kind=rp), intent(in) :: v_qmm(n1,n2), ef_qmd(3,n3,n4)

            !call print_matrix(.true.,'VQM:',mm_atoms,ld_cart,mm_atoms,ld_cart,transpose(v_qmm))
            !call print_matrix(.true.,'EQM 1:',pol_atoms,3,pol_atoms,3,transpose(ef_qmd(:,:,1))) 
            !call print_matrix(.true.,'EQM 2:',pol_atoms,3,pol_atoms,3,transpose(ef_qmd(:,:,2))) 

            ! Initialize variables
            ipd = Zero

            ! Compute polarization
            nmax = 200
            !convergence = 5.0E-10
            call polarization(solver,ef_qd+ef_qmd,ipd) ! - jacobi interation solver
            !call polarization(1,ef_qd+ef_qmd,ipd) 

        end subroutine

        subroutine restart() bind(c, name='restart')
            implicit none

            v_qq   = Zero
            ef_qd  = Zero
            dv_qq  = Zero
            def_qd = Zero
            ipd    = Zero

        end subroutine

        subroutine get_energy(EMM,EPol) bind(c, name='get_energy')
            ! Get the energy

            implicit none
            real(kind=rp), intent(out) :: EMM, EPol

            call energy(int(0,ip),EMM)
            call energy(int(1,ip),EPol)

        end subroutine
end module mod_interface

