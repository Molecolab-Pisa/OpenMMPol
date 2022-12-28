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

    contains

        subroutine analytical_fixedelec_geomgrad(s, grad)
            use mod_electrostatics, only: prepare_M2M, ommp_electrostatics_type

            implicit none
            
            type(ommp_system), intent(inout), target :: s
            !! System data structure
            real(rp), dimension(3,s%top%mm_atoms), intent(inout) :: grad
            !! Geometrical gradients in output, results will be added
            
            integer(ip) :: i
            type(ommp_electrostatics_type), pointer :: eel 
            eel => s%eel

            call prepare_M2M(eel, .true.)

            if(eel%amoeba) then
                call fatal_error("Not Implemented")
            else
                do i=1, s%top%mm_atoms
                    grad(:,i) = grad(:,i) + 0.5 * eel%q(1,i) * eel%E_M2M(:,i)
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
