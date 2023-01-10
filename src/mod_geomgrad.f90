#include "f_cart_components.h"

module mod_geomgrad
    use mod_io, only: fatal_error, ommp_message
    use mod_memory, only: ip, rp
    use mod_mmpol, only: ommp_system
    use mod_topology, only: ommp_topology_type

    implicit none
    private
    
    abstract interface
    function energy_term(s)
            use mod_mmpol, only: ommp_system
            use mod_memory, only: rp
            type(ommp_system), intent(inout), target :: s
            real(rp) :: energy_term
    end function
    end interface

    public :: numerical_fixedelec_geomgrad, analytical_fixedelec_geomgrad
    public :: numerical_polelec_geomgrad, analytical_polelec_geomgrad

    contains

        subroutine analytical_fixedelec_geomgrad(s, grad)
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

            call prepare_fixedelec(eel, .true.)

            if(eel%amoeba) then
                do i=1, s%top%mm_atoms
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
                ! Torque forces from multipoles rotation
                call fixedelec_rotation_geomgrad(eel, grad)
            else
                do i=1, s%top%mm_atoms
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
        end subroutine

        subroutine numerical_fixedelec_geomgrad(s, grad)
            use ommp_interface, only: ommp_get_fixedelec_energy

            implicit none
            
            type(ommp_system), intent(inout) :: s
            !! System data structure
            real(rp), dimension(3,s%top%mm_atoms), intent(inout) :: grad
            !! Geometrical gradients in output, results will be added

            procedure(energy_term), pointer :: ene_f => ommp_get_fixedelec_energy

            call numerical_geomgrad(s, ene_f, grad)
        end subroutine
        
        subroutine analytical_polelec_geomgrad(s, grad)
            !use mod_electrostatics, only: prepare_M2D, ommp_electrostatics_type
            use mod_polarization, only: polarization
            use mod_constants, only: OMMP_AMOEBA_D, OMMP_AMOEBA_P
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
                do i=1, eel%top%mm_atoms
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
               
                do i=1, eel%pol_atoms
                    ! \mu_D Egrd_P
                    grad(_x_,i) = grad(_x_,i) &
                                  + 0.5*eel%ipd(_x_,i,OMMP_AMOEBA_D) * (eel%Egrd_M2D(_xx_,i,OMMP_AMOEBA_P) &
                                                  + eel%Egrd_D2D(_xx_,i,OMMP_AMOEBA_P)) &
                                  + 0.5*eel%ipd(_y_,i,OMMP_AMOEBA_D) * (eel%Egrd_M2D(_xy_,i,OMMP_AMOEBA_P) &
                                                  + eel%Egrd_D2D(_xy_,i,OMMP_AMOEBA_P)) &
                                  + 0.5*eel%ipd(_z_,i,OMMP_AMOEBA_D) * (eel%Egrd_M2D(_xz_,i,OMMP_AMOEBA_P) & 
                                                  + eel%Egrd_D2D(_xz_,i,OMMP_AMOEBA_P)) 
                    grad(_y_,i) = grad(_y_,i) &
                                  + 0.5*eel%ipd(_x_,i,OMMP_AMOEBA_D) * (eel%Egrd_M2D(_yx_,i,OMMP_AMOEBA_P) &
                                                  + eel%Egrd_D2D(_yx_,i,OMMP_AMOEBA_P)) &
                                  + 0.5*eel%ipd(_y_,i,OMMP_AMOEBA_D) * (eel%Egrd_M2D(_yy_,i,OMMP_AMOEBA_P) &
                                                  + eel%Egrd_D2D(_yy_,i,OMMP_AMOEBA_P)) &
                                  + 0.5*eel%ipd(_z_,i,OMMP_AMOEBA_D) * (eel%Egrd_M2D(_yz_,i,OMMP_AMOEBA_P) &
                                                  + eel%Egrd_D2D(_yz_,i,OMMP_AMOEBA_P)) 
                    grad(_z_,i) = grad(_z_,i) &
                                  + 0.5*eel%ipd(_x_,i,OMMP_AMOEBA_D) * (eel%Egrd_M2D(_zx_,i,OMMP_AMOEBA_P) &
                                                  + eel%Egrd_D2D(_zx_,i,OMMP_AMOEBA_P)) &
                                  + 0.5*eel%ipd(_y_,i,OMMP_AMOEBA_D) * (eel%Egrd_M2D(_zy_,i,OMMP_AMOEBA_P) &
                                                  + eel%Egrd_D2D(_zy_,i,OMMP_AMOEBA_P)) &
                                  + 0.5*eel%ipd(_z_,i,OMMP_AMOEBA_D) * (eel%Egrd_M2D(_zz_,i,OMMP_AMOEBA_P) &
                                                  + eel%Egrd_D2D(_zz_,i,OMMP_AMOEBA_P))
                    ! \mu_P Egrd_D
                    grad(_x_,i) = grad(_x_,i) &
                                  + 0.5*eel%ipd(_x_,i,OMMP_AMOEBA_P) * (eel%Egrd_M2D(_xx_,i,OMMP_AMOEBA_D) &
                                                  + eel%Egrd_D2D(_xx_,i,OMMP_AMOEBA_D)) &
                                  + 0.5*eel%ipd(_y_,i,OMMP_AMOEBA_P) * (eel%Egrd_M2D(_xy_,i,OMMP_AMOEBA_D) &
                                                  + eel%Egrd_D2D(_xy_,i,OMMP_AMOEBA_D)) &
                                  + 0.5*eel%ipd(_z_,i,OMMP_AMOEBA_P) * (eel%Egrd_M2D(_xz_,i,OMMP_AMOEBA_D) & 
                                                  + eel%Egrd_D2D(_xz_,i,OMMP_AMOEBA_D)) 
                    grad(_y_,i) = grad(_y_,i) &
                                  + 0.5*eel%ipd(_x_,i,OMMP_AMOEBA_P) * (eel%Egrd_M2D(_yx_,i,OMMP_AMOEBA_D) &
                                                  + eel%Egrd_D2D(_yx_,i,OMMP_AMOEBA_D)) &
                                  + 0.5*eel%ipd(_y_,i,OMMP_AMOEBA_P) * (eel%Egrd_M2D(_yy_,i,OMMP_AMOEBA_D) &
                                                  + eel%Egrd_D2D(_yy_,i,OMMP_AMOEBA_D)) &
                                  + 0.5*eel%ipd(_z_,i,OMMP_AMOEBA_P) * (eel%Egrd_M2D(_yz_,i,OMMP_AMOEBA_D) &
                                                  + eel%Egrd_D2D(_yz_,i,OMMP_AMOEBA_D)) 
                    grad(_z_,i) = grad(_z_,i) &
                                  + 0.5*eel%ipd(_x_,i,OMMP_AMOEBA_P) * (eel%Egrd_M2D(_zx_,i,OMMP_AMOEBA_D) &
                                                  + eel%Egrd_D2D(_zx_,i,OMMP_AMOEBA_D)) &
                                  + 0.5*eel%ipd(_y_,i,OMMP_AMOEBA_P) * (eel%Egrd_M2D(_zy_,i,OMMP_AMOEBA_D) &
                                                  + eel%Egrd_D2D(_zy_,i,OMMP_AMOEBA_D)) &
                                  + 0.5*eel%ipd(_z_,i,OMMP_AMOEBA_P) * (eel%Egrd_M2D(_zz_,i,OMMP_AMOEBA_D) &
                                                  + eel%Egrd_D2D(_zz_,i,OMMP_AMOEBA_D))
                end do
            else
                do i=1, eel%top%mm_atoms
                    grad(:,i) = grad(:,i) - eel%q(1,i) * eel%E_D2M(:,i)
                end do
                do i=1, eel%pol_atoms
                    grad(_x_,i) = grad(_x_,i) &
                                  + eel%ipd(_x_,i,1) * (eel%Egrd_M2D(_xx_,i,1) &
                                                      + eel%Egrd_D2D(_xx_,i,1)) &
                                  + eel%ipd(_y_,i,1) * (eel%Egrd_M2D(_xy_,i,1) &
                                                      + eel%Egrd_D2D(_xy_,i,1)) &
                                  + eel%ipd(_z_,i,1) * (eel%Egrd_M2D(_xz_,i,1) & 
                                                      + eel%Egrd_D2D(_xz_,i,1)) 
                    grad(_y_,i) = grad(_y_,i) &
                                  + eel%ipd(_x_,i,1) * (eel%Egrd_M2D(_yx_,i,1) &
                                                      + eel%Egrd_D2D(_yx_,i,1)) &
                                  + eel%ipd(_y_,i,1) * (eel%Egrd_M2D(_yy_,i,1) &
                                                      + eel%Egrd_D2D(_yy_,i,1)) &
                                  + eel%ipd(_z_,i,1) * (eel%Egrd_M2D(_yz_,i,1) &
                                                      + eel%Egrd_D2D(_yz_,i,1)) 
                    grad(_z_,i) = grad(_z_,i) &
                                  + eel%ipd(_x_,i,1) * (eel%Egrd_M2D(_zx_,i,1) &
                                                      + eel%Egrd_D2D(_zx_,i,1)) &
                                  + eel%ipd(_y_,i,1) * (eel%Egrd_M2D(_zy_,i,1) &
                                                      + eel%Egrd_D2D(_zy_,i,1)) &
                                  + eel%ipd(_z_,i,1) * (eel%Egrd_M2D(_zz_,i,1) &
                                                      + eel%Egrd_D2D(_zz_,i,1)) 
                end do
            end if

            call polelec_rotation_geomgrad(eel, grad)
        end subroutine
        
        subroutine numerical_polelec_geomgrad(s, grad)
            use ommp_interface, only: ommp_get_polelec_energy

            implicit none
            
            type(ommp_system), intent(inout) :: s
            !! System data structure
            real(rp), dimension(3,s%top%mm_atoms), intent(inout) :: grad
            !! Geometrical gradients in output, results will be added

            procedure(energy_term), pointer :: ene_f => ommp_get_polelec_energy

            call numerical_geomgrad(s, ene_f, grad)
        end subroutine
        
        subroutine numerical_geomgrad(s, ene_f, grad)
            use mod_mmpol, only: update_coordinates
            use mod_memory, only: mallocate, mfree
            implicit none
            
            type(ommp_system), intent(inout) :: s
            !! System data structure
            procedure(energy_term), pointer :: ene_f
            !! The energy function (from interface module) for which
            !! numerical gradients are needed
            real(rp), dimension(3,s%top%mm_atoms), intent(inout) :: grad
            !! Geometrical gradients in output, results will be added

            integer(ip) :: i, j
            real(rp), allocatable :: new_c(:,:)
            real(rp) :: tmp
            real(rp), parameter :: dd = 1e-5

            call mallocate('numerical_geomgrad [new_c]', &
                           3, s%top%mm_atoms, new_c)
            new_c = s%top%cmm
            
            do i=1, s%top%mm_atoms
                do j=1, 3
                    new_c(j,i) = new_c(j,i) + dd
                    call update_coordinates(s, new_c)
                    tmp = ene_f(s)

                    new_c(j,i) = new_c(j,i) - 2*dd
                    call update_coordinates(s, new_c)
                    tmp = tmp - ene_f(s)
                    
                    grad(j,i) = grad(j,i) + tmp / (2*dd)
                    
                    new_c(j,i) = new_c(j,i) + dd
                    call update_coordinates(s, new_c)
                end do
            end do
            
            call mfree('numerical_geomgrad [new_c]', new_c)
        end subroutine
end module
