subroutine forces(isrc,itrg,ff_qq,ff_qd,ff_dq,ff_dd)
  use mmpol
  use force
  use precision
  implicit none
!
! main driver for the calculation of electrostatic properties.
! ============================================================
!
!   isrc ... 0: the sources are the charges (multipoles), contained in the array
!               q and with coordinates cmm (see the mmpol module)
!            1: the sources are the ipds, contained in the array ipd and with
!               coordinates cpol (see the mmpol module)
!
!   itrg ... 0: compute the forces at cmm
!            1: compute the forces at cpol
!
 
!
  integer(ip) :: isrc,itrg
  real(rp),   dimension(3,mm_atoms),    intent(inout) :: ff_qq
  real(rp),   dimension(3,mm_atoms),    intent(inout) :: ff_dq
  real(rp),   dimension(3,pol_atoms),   intent(inout) :: ff_qd
  real(rp),   dimension(3,pol_atoms),   intent(inout) :: ff_dd
!
! internal variables:
!
  real(rp),   dimension(3,mm_atoms) ::  ff_torq_qq, ff_torq_dq
!
! Initialize variables
!
ff_torq_qq = Zero
ff_torq_dq = Zero
!
! figure how what to compute
!
if (itrg.eq.0) then       ! Calculate forces at MM sites
    if (isrc.eq.0) then   ! Sources are MM sites
        !
        call forces_MM(dv_qq,ff_qq)
        !
        ! For Amoeba force-field calculate torque forces
        !
        if (Amoeba) then
            call rotate_multipoles(.true.,dv_qq,ff_torq_qq)
        end if
        !
        if (verbose.ge.4) then
            call print_matrix(.true.,'q E(MM) (force from MM at MM) ',mm_atoms,3,mm_atoms,3,transpose(ff_qq))
            call print_matrix(.true.,'qE(MM) Trq (Torque forces MM multipoles): ',mm_atoms,3,mm_atoms,3,transpose(ff_torq_qq))
        end if
        !
        ! Add torque forces to the multipole forces
        !
        ff_qq = ff_qq + ff_torq_qq
        !
    elseif (isrc.eq.1) then   ! Sources are idp
        !
        call forces_MM(dv_dq,ff_dq)
        !
        ! For Amoeba force-field calculate torque forces
        !
        if (Amoeba) then
            call rotate_multipoles(.true.,dv_dq,ff_torq_dq)
        end if
        !
        ! Add torque forces to the multipole forces
        !
        if (verbose.ge.4) then
            call print_matrix(.true.,'q E(mu) (force from ipd at MM) ',mm_atoms,3,mm_atoms,3,transpose(ff_dq))
            call print_matrix(.true.,'q E(mu) Trq (Torque forces from ipd): ',mm_atoms,3,mm_atoms,3,transpose(ff_torq_dq))
        end if
        !
        ff_dq = ff_dq + ff_torq_dq
        !
    else
        write(*,*) "Unknown value of isrc in forces."
    end if
elseif (itrg.eq.1) then       ! Calculate forces at ipd
    if (isrc.eq.0) then       ! Sources are MM sites
        !
        call forces_DD(def_qd,ff_qd) 
        !
        if (verbose.ge.4) then
            call print_matrix(.true.,'mu G(MM) (force from MM at ipd) ',pol_atoms,3,pol_atoms,3,transpose(ff_qd))
        end if
        !
    elseif (isrc.eq.1) then   ! Sources are idp
        !
        call forces_DD(def_dd,ff_dd)
        !
        if (verbose.ge.4) then
            call print_matrix(.true.,'mu G(MM) (force from ipd at ipd) ',pol_atoms,3,pol_atoms,3,transpose(ff_dd))
        end if
        !
    else
        write(*,*) "Unknown value of isrc in forces."
    end if
else
    write(*,*) "Unknown value of itrg in forces."
end if 

end subroutine forces