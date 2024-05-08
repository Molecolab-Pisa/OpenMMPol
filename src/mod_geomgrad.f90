#include "f_cart_components.h"

module mod_geomgrad
    use mod_io, only: fatal_error, ommp_message
    use mod_memory, only: ip, rp
    use mod_mmpol, only: ommp_system
    use mod_topology, only: ommp_topology_type
    use mod_profiling, only: time_push, time_pull

    implicit none
    private
    
    public :: fixedelec_geomgrad, polelec_geomgrad

    contains

        subroutine fixedelec_geomgrad(s, grad)
            use mod_electrostatics, only: prepare_fixedelec, &
                                          ommp_electrostatics_type

            implicit none
            
            type(ommp_system), intent(inout), target :: s
            !! System data structure
            real(rp), dimension(3,s%top%mm_atoms), intent(inout) :: grad
            !! Geometrical gradients in output, results will be added
            
            integer(ip) :: i
            type(ommp_electrostatics_type), pointer :: eel 
            eel => s%eel
            
            call time_push
            call prepare_fixedelec(eel, .true.)
            call time_pull("Prepare fixedelec")

            call time_push
            if(eel%amoeba) then
                !$omp parallel do 
                do i=1, s%top%mm_atoms
                    ! If the atom is frozen, there are no contribution to compute
                    if(s%top%use_frozen) then
                        if(s%top%frozen(i)) cycle
                    end if

                    ! Charges -qE
                    grad(:,i) = grad(:,i) - eel%q(1,i) * eel%E_M2M(:,i)
                    
                    ! Dipoles mu \nablaE
                    grad(_x_,i) = grad(_x_,i) &
                                  + eel%q(1+_x_,i) * eel%Egrd_M2M(_xx_,i) &
                                  + eel%q(1+_y_,i) * eel%Egrd_M2M(_xy_,i) &
                                  + eel%q(1+_z_,i) * eel%Egrd_M2M(_xz_,i)
                    grad(_y_,i) = grad(_y_,i) &
                                  + eel%q(1+_x_,i) * eel%Egrd_M2M(_yx_,i) &
                                  + eel%q(1+_y_,i) * eel%Egrd_M2M(_yy_,i) &
                                  + eel%q(1+_z_,i) * eel%Egrd_M2M(_yz_,i)
                    grad(_z_,i) = grad(_z_,i) &
                                  + eel%q(1+_x_,i) * eel%Egrd_M2M(_zx_,i) &
                                  + eel%q(1+_y_,i) * eel%Egrd_M2M(_zy_,i) &
                                  + eel%q(1+_z_,i) * eel%Egrd_M2M(_zz_,i)
                    
                    ! Quadrupoles Q \nabla^2E
                    grad(_x_,i) = grad(_x_,i) &
                                  - eel%q(4+_xx_,i) * eel%EHes_M2M(_xxx_,i) &
                                  - eel%q(4+_yy_,i) * eel%EHes_M2M(_yyx_,i) &
                                  - eel%q(4+_zz_,i) * eel%EHes_M2M(_zzx_,i) &
                                  - 2*(eel%q(4+_xy_,i) * eel%EHes_M2M(_xyx_,i) &
                                  +    eel%q(4+_xz_,i) * eel%EHes_M2M(_xzx_,i) &
                                  +    eel%q(4+_yz_,i) * eel%EHes_M2M(_yzx_,i))
                    grad(_y_,i) = grad(_y_,i) &
                                  - eel%q(4+_xx_,i) * eel%EHes_M2M(_xxy_,i) &
                                  - eel%q(4+_yy_,i) * eel%EHes_M2M(_yyy_,i) &
                                  - eel%q(4+_zz_,i) * eel%EHes_M2M(_zzy_,i) &
                                  - 2*(eel%q(4+_xy_,i) * eel%EHes_M2M(_xyy_,i) &
                                  +    eel%q(4+_xz_,i) * eel%EHes_M2M(_xzy_,i) &
                                  +    eel%q(4+_yz_,i) * eel%EHes_M2M(_yzy_,i))
                    grad(_z_,i) = grad(_z_,i) &
                                  - eel%q(4+_xx_,i) * eel%EHes_M2M(_xxz_,i) &
                                  - eel%q(4+_yy_,i) * eel%EHes_M2M(_yyz_,i) &
                                  - eel%q(4+_zz_,i) * eel%EHes_M2M(_zzz_,i) &
                                  - 2*(eel%q(4+_xy_,i) * eel%EHes_M2M(_xyz_,i) &
                                  +    eel%q(4+_xz_,i) * eel%EHes_M2M(_xzz_,i) &
                                  +    eel%q(4+_yz_,i) * eel%EHes_M2M(_yzz_,i))
                    
                end do
                call time_push
                ! Torque forces from multipoles rotation
                call rotation_geomgrad(eel, eel%E_M2M, eel%Egrd_M2M, grad)
                call time_pull("Rotation grad")
            else
                !$omp parallel do 
                do i=1, s%top%mm_atoms
                    ! Skip frozen atoms contributions
                    if(s%top%use_frozen) then
                        if(s%top%frozen(i)) cycle
                    end if
                    
                    ! Here the minus sign is due to the definition of Electric
                    ! field that is E = -\nabla V.
                    ! The factor 1/2 that is present in the definition of the
                    ! energy disappear because in the derivative there is a 
                    ! factor 2 that takes in account that when an atom is 
                    ! "displaced", this affect both the point at which the 
                    ! potential of all the other charges is computed and
                    ! the potential computed at the -fixed- position of all
                    ! the other charges (due to the displaced source).
                    grad(:,i) = grad(:,i) - eel%q(1,i) * eel%E_M2M(:,i)
                end do
            end if
            call time_pull("Grad sum")
        end subroutine

        subroutine polelec_geomgrad(s, grad)
            !use mod_electrostatics, only: prepare_M2D, ommp_electrostatics_type
            use mod_polarization, only: polarization
            use mod_electrostatics

            implicit none
            
            type(ommp_system), intent(inout), target :: s
            !! System data structure
            real(rp), dimension(3,s%top%mm_atoms), intent(inout) :: grad
            !! Geometrical gradients in output, results will be added
            
            integer(ip) :: i
            type(ommp_electrostatics_type), pointer :: eel 
            eel => s%eel

            if(.not. eel%ipd_done) then
                call prepare_polelec(eel, .false.)
                call polarization(s, eel%e_M2D)
            end if
            call prepare_polelec(eel, .true.)

            if(eel%amoeba) then
                !$omp parallel do 
                do i=1, eel%top%mm_atoms
                    ! Skip frozen atoms contributions
                    if(s%top%use_frozen) then
                        if(s%top%frozen(i)) cycle
                    end if
                    
                    ! Charges q E
                    grad(:,i) = grad(:,i) - eel%q(1,i) * eel%E_D2M(:,i)
                    
                    ! Dipoles mu \nablaE
                    grad(_x_,i) = grad(_x_,i) &
                                  + eel%q(1+_x_,i) * eel%Egrd_D2M(_xx_,i) &
                                  + eel%q(1+_y_,i) * eel%Egrd_D2M(_xy_,i) &
                                  + eel%q(1+_z_,i) * eel%Egrd_D2M(_xz_,i)
                    grad(_y_,i) = grad(_y_,i) &
                                  + eel%q(1+_x_,i) * eel%Egrd_D2M(_yx_,i) &
                                  + eel%q(1+_y_,i) * eel%Egrd_D2M(_yy_,i) &
                                  + eel%q(1+_z_,i) * eel%Egrd_D2M(_yz_,i)
                    grad(_z_,i) = grad(_z_,i) &
                                  + eel%q(1+_x_,i) * eel%Egrd_D2M(_zx_,i) &
                                  + eel%q(1+_y_,i) * eel%Egrd_D2M(_zy_,i) &
                                  + eel%q(1+_z_,i) * eel%Egrd_D2M(_zz_,i)
                    
                    ! Quadrupoles Q \nabla^2E
                    grad(_x_,i) = grad(_x_,i) &
                                  - eel%q(4+_xx_,i) * eel%EHes_D2M(_xxx_,i) &
                                  - eel%q(4+_yy_,i) * eel%EHes_D2M(_yyx_,i) &
                                  - eel%q(4+_zz_,i) * eel%EHes_D2M(_zzx_,i) &
                                  - 2*(eel%q(4+_xy_,i) * eel%EHes_D2M(_xyx_,i) &
                                  +    eel%q(4+_xz_,i) * eel%EHes_D2M(_xzx_,i) &
                                  +    eel%q(4+_yz_,i) * eel%EHes_D2M(_yzx_,i))
                    grad(_y_,i) = grad(_y_,i) &
                                  - eel%q(4+_xx_,i) * eel%EHes_D2M(_xxy_,i) &
                                  - eel%q(4+_yy_,i) * eel%EHes_D2M(_yyy_,i) &
                                  - eel%q(4+_zz_,i) * eel%EHes_D2M(_zzy_,i) &
                                  - 2*(eel%q(4+_xy_,i) * eel%EHes_D2M(_xyy_,i) &
                                  +    eel%q(4+_xz_,i) * eel%EHes_D2M(_xzy_,i) &
                                  +    eel%q(4+_yz_,i) * eel%EHes_D2M(_yzy_,i))
                    grad(_z_,i) = grad(_z_,i) &
                                  - eel%q(4+_xx_,i) * eel%EHes_D2M(_xxz_,i) &
                                  - eel%q(4+_yy_,i) * eel%EHes_D2M(_yyz_,i) &
                                  - eel%q(4+_zz_,i) * eel%EHes_D2M(_zzz_,i) &
                                  - 2*(eel%q(4+_xy_,i) * eel%EHes_D2M(_xyz_,i) &
                                  +    eel%q(4+_xz_,i) * eel%EHes_D2M(_xzz_,i) &
                                  +    eel%q(4+_yz_,i) * eel%EHes_D2M(_yzz_,i))
                end do
              
                !$omp parallel do 
                do i=1, eel%pol_atoms
                    ! Skip frozen atoms contributions
                    if(s%top%use_frozen) then
                        if(s%top%frozen(eel%polar_mm(i))) cycle
                    end if
                    
                    ! \mu_D Egrd_P
                    grad(_x_,eel%polar_mm(i)) = grad(_x_,eel%polar_mm(i)) &
                                  + 0.5*eel%ipd(_x_,i,_amoeba_D_) * (eel%Egrd_M2D(_xx_,i,_amoeba_P_) &
                                                  + eel%Egrd_D2D(_xx_,i,_amoeba_P_)) &
                                  + 0.5*eel%ipd(_y_,i,_amoeba_D_) * (eel%Egrd_M2D(_xy_,i,_amoeba_P_) &
                                                  + eel%Egrd_D2D(_xy_,i,_amoeba_P_)) &
                                  + 0.5*eel%ipd(_z_,i,_amoeba_D_) * (eel%Egrd_M2D(_xz_,i,_amoeba_P_) & 
                                                  + eel%Egrd_D2D(_xz_,i,_amoeba_P_)) 
                    grad(_y_,eel%polar_mm(i)) = grad(_y_,eel%polar_mm(i)) &
                                  + 0.5*eel%ipd(_x_,i,_amoeba_D_) * (eel%Egrd_M2D(_yx_,i,_amoeba_P_) &
                                                  + eel%Egrd_D2D(_yx_,i,_amoeba_P_)) &
                                  + 0.5*eel%ipd(_y_,i,_amoeba_D_) * (eel%Egrd_M2D(_yy_,i,_amoeba_P_) &
                                                  + eel%Egrd_D2D(_yy_,i,_amoeba_P_)) &
                                  + 0.5*eel%ipd(_z_,i,_amoeba_D_) * (eel%Egrd_M2D(_yz_,i,_amoeba_P_) &
                                                  + eel%Egrd_D2D(_yz_,i,_amoeba_P_)) 
                    grad(_z_,eel%polar_mm(i)) = grad(_z_,eel%polar_mm(i)) &
                                  + 0.5*eel%ipd(_x_,i,_amoeba_D_) * (eel%Egrd_M2D(_zx_,i,_amoeba_P_) &
                                                  + eel%Egrd_D2D(_zx_,i,_amoeba_P_)) &
                                  + 0.5*eel%ipd(_y_,i,_amoeba_D_) * (eel%Egrd_M2D(_zy_,i,_amoeba_P_) &
                                                  + eel%Egrd_D2D(_zy_,i,_amoeba_P_)) &
                                  + 0.5*eel%ipd(_z_,i,_amoeba_D_) * (eel%Egrd_M2D(_zz_,i,_amoeba_P_) &
                                                  + eel%Egrd_D2D(_zz_,i,_amoeba_P_))
                    ! \mu_P Egrd_D
                    grad(_x_,eel%polar_mm(i)) = grad(_x_,eel%polar_mm(i)) &
                                  + 0.5*eel%ipd(_x_,i,_amoeba_P_) * (eel%Egrd_M2D(_xx_,i,_amoeba_D_) &
                                                  + eel%Egrd_D2D(_xx_,i,_amoeba_D_)) &
                                  + 0.5*eel%ipd(_y_,i,_amoeba_P_) * (eel%Egrd_M2D(_xy_,i,_amoeba_D_) &
                                                  + eel%Egrd_D2D(_xy_,i,_amoeba_D_)) &
                                  + 0.5*eel%ipd(_z_,i,_amoeba_P_) * (eel%Egrd_M2D(_xz_,i,_amoeba_D_) & 
                                                  + eel%Egrd_D2D(_xz_,i,_amoeba_D_)) 
                    grad(_y_,eel%polar_mm(i)) = grad(_y_,eel%polar_mm(i)) &
                                  + 0.5*eel%ipd(_x_,i,_amoeba_P_) * (eel%Egrd_M2D(_yx_,i,_amoeba_D_) &
                                                  + eel%Egrd_D2D(_yx_,i,_amoeba_D_)) &
                                  + 0.5*eel%ipd(_y_,i,_amoeba_P_) * (eel%Egrd_M2D(_yy_,i,_amoeba_D_) &
                                                  + eel%Egrd_D2D(_yy_,i,_amoeba_D_)) &
                                  + 0.5*eel%ipd(_z_,i,_amoeba_P_) * (eel%Egrd_M2D(_yz_,i,_amoeba_D_) &
                                                  + eel%Egrd_D2D(_yz_,i,_amoeba_D_)) 
                    grad(_z_,eel%polar_mm(i)) = grad(_z_,eel%polar_mm(i)) &
                                  + 0.5*eel%ipd(_x_,i,_amoeba_P_) * (eel%Egrd_M2D(_zx_,i,_amoeba_D_) &
                                                  + eel%Egrd_D2D(_zx_,i,_amoeba_D_)) &
                                  + 0.5*eel%ipd(_y_,i,_amoeba_P_) * (eel%Egrd_M2D(_zy_,i,_amoeba_D_) &
                                                  + eel%Egrd_D2D(_zy_,i,_amoeba_D_)) &
                                  + 0.5*eel%ipd(_z_,i,_amoeba_P_) * (eel%Egrd_M2D(_zz_,i,_amoeba_D_) &
                                                  + eel%Egrd_D2D(_zz_,i,_amoeba_D_))
                end do
            else
                do i=1, eel%top%mm_atoms
                    ! Skip frozen atoms contributions
                    if(s%top%use_frozen) then
                        if(s%top%frozen(i)) cycle
                    end if

                    grad(:,i) = grad(:,i) - eel%q(1,i) * eel%E_D2M(:,i)
                end do
                
                do i=1, eel%pol_atoms
                    ! Skip frozen atoms contributions
                    if(s%top%use_frozen) then
                        if(s%top%frozen(eel%polar_mm(i))) cycle
                    end if
                    
                    grad(_x_,eel%polar_mm(i)) = grad(_x_,eel%polar_mm(i)) &
                                  + eel%ipd(_x_,i,1) * (eel%Egrd_M2D(_xx_,i,1) &
                                                      + eel%Egrd_D2D(_xx_,i,1)) &
                                  + eel%ipd(_y_,i,1) * (eel%Egrd_M2D(_xy_,i,1) &
                                                      + eel%Egrd_D2D(_xy_,i,1)) &
                                  + eel%ipd(_z_,i,1) * (eel%Egrd_M2D(_xz_,i,1) & 
                                                      + eel%Egrd_D2D(_xz_,i,1)) 
                    grad(_y_,eel%polar_mm(i)) = grad(_y_,eel%polar_mm(i)) &
                                  + eel%ipd(_x_,i,1) * (eel%Egrd_M2D(_yx_,i,1) &
                                                      + eel%Egrd_D2D(_yx_,i,1)) &
                                  + eel%ipd(_y_,i,1) * (eel%Egrd_M2D(_yy_,i,1) &
                                                      + eel%Egrd_D2D(_yy_,i,1)) &
                                  + eel%ipd(_z_,i,1) * (eel%Egrd_M2D(_yz_,i,1) &
                                                      + eel%Egrd_D2D(_yz_,i,1)) 
                    grad(_z_,eel%polar_mm(i)) = grad(_z_,eel%polar_mm(i)) &
                                  + eel%ipd(_x_,i,1) * (eel%Egrd_M2D(_zx_,i,1) &
                                                      + eel%Egrd_D2D(_zx_,i,1)) &
                                  + eel%ipd(_y_,i,1) * (eel%Egrd_M2D(_zy_,i,1) &
                                                      + eel%Egrd_D2D(_zy_,i,1)) &
                                  + eel%ipd(_z_,i,1) * (eel%Egrd_M2D(_zz_,i,1) &
                                                      + eel%Egrd_D2D(_zz_,i,1)) 
                end do
            end if

            if(eel%amoeba) call rotation_geomgrad(eel, eel%E_D2M, eel%Egrd_D2M, grad)
        end subroutine
end module
