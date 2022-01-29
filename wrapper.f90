! Wrapper function for open-mmpol library
subroutine w_mmpol_init(filename)
  use mmpol
  ! 
  implicit none
  character(len=120), intent(in) :: filename
  character(len=120) :: BsName

  input_file = filename
  BsName = filename(1:Scan(filename,'.',BACK=.true.))
  scratch_file = Trim(BsName)//'rwf'

  call mmpol_init
  
end subroutine 

subroutine do_mm 
  use precision
  use mmpol
  ! Perform the MM-only calculations of pot and field
  ! and return EMM

  implicit none
  
  call electrostatics(0,0,0,v_qq,ef_qd,dv_qq,def_qd)
  call electrostatics(0,1,0,v_qq,ef_qd,dv_qq,def_qd)

end subroutine 


subroutine do_qmmm(v_qmm,ef_qmd,n1,n2,n3,n4)
  use precision
  use mmpol
  ! Compute polarization and QM-MM energy

  implicit none
  ! real(kind=rp), intent(in) :: v_qmm(ld_cart,mm_atoms), ef_qmd(3,pol_atoms,n_ipd)
  real(kind=rp), intent(in) :: v_qmm(n1,n2), ef_qmd(3,n3,n4)
  integer(kind=ip), intent(in) :: n1,n2,n3,n4

  ! Compute polarization
  call polarization(1,ef_qd+ef_qmd,ipd) 

end subroutine


subroutine get_energy(EMM,EPol)
  use precision
  use mmpol
  ! Get the energy

  implicit none
  real(kind=rp), intent(out) :: EMM, EPol

  call energy(0,EMM)
  call energy(1,EPol)

end subroutine



