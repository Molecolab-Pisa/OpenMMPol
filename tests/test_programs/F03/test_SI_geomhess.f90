program test_SI_geomhess
    use iso_c_binding, only: c_char
    use ommp_interface
    use mod_geomgrad
    use mod_geomhess
    use mod_memory, only : ip, rp
    use mod_mmpol

    implicit none
    character(kind=c_char, len=120), dimension(3) :: args
    character(len=OMMP_STR_CHAR_MAX) :: prm_file
    integer :: narg
    type(ommp_system), pointer :: my_system, fake_qm
    type(ommp_qm_helper), pointer :: my_qmh
    logical :: use_qm = .false., use_fake_qm = .false.
    integer(ip) :: i, j, ix, jx, n
!
    real(rp), parameter   :: delta = 1.0e-5_rp
    real(rp), allocatable :: hess(:,:,:,:), gph(:,:), gmh(:,:), hessnum(:,:,:,:), cmm(:,:), gg(:,:)
  
    narg = command_argument_count()
    if (narg /= 1) then
        write(6, *) "Syntax expected "
        write(6, *) "   $ test_SI_geomhess.exe <JSON FILE>"
        call exit(1)
    else 
        call get_command_argument(1, args(1))
        call ommp_smartinput(trim(args(1)), my_system, my_qmh)
!
        n = my_system%top%mm_atoms
        write(6,*) '# atoms = ', n
        allocate(hess(3,3,n,n),hessnum(3,3,n,n))
        hess = 0.0_rp
        hessnum = 0.0_rp
!fl!!
!       do i = 1, n
!         my_system%eel%q(1,i) = 0.0_rp
!         my_system%eel%q(2:4,i) = 0.0_rp
!         my_system%eel%q(5:10,i) = 0.0_rp
!       end do
!!!!!
        call fixedelec_geomhess(my_system, hess)
!       do i = 1, n
!         do j = 1, n
!           write(6,*) 'hess for atoms ',j, i
!           write(6,'(10f14.8)') hess(:,:,j,i)
!         end do
!       end do
!       write(6,*) 'charges:'
!       do i = 1, n
!         write(6,*) my_system%eel%q(1,i)
!       end do
!
!       write(6,*) 'field:'
!       do i = 1, n
!         write(6,'(3f12.6)') my_system%eel%e_M2M(:,i)
!       end do
!       write(6,*) 'dipoles:'
!       do i = 1, n
!         write(6,'(3f12.8)') my_system%eel%q(2:4,i)
!       end do
!       write(6,*) 'field gradient:'
!       do i = 1, n
!         write(6,'(10f12.8)') my_system%eel%Egrd_M2M(:,i)
!       end do
!       write(6,*) 'field hessian:'
!       do i = 1, n
!         write(6,'(10f12.8)') my_system%eel%Ehes_M2M(:,i)
!       end do
!
!fl
!       gg = 0.0_rp
!       do i = 1, n
!         gg(1,i) = my_system%eel%q(2,i) * my_system%eel%Egrd_M2M(1,i) + &
!                   my_system%eel%q(3,i) * my_system%eel%Egrd_M2M(2,i) +  &
!                   my_system%eel%q(4,i) * my_system%eel%Egrd_M2M(4,i)
!         gg(2,i) = my_system%eel%q(2,i) * my_system%eel%Egrd_M2M(2,i) + &
!                   my_system%eel%q(3,i) * my_system%eel%Egrd_M2M(3,i) +  &
!                   my_system%eel%q(4,i) * my_system%eel%Egrd_M2M(5,i)
!         gg(3,i) = my_system%eel%q(2,i) * my_system%eel%Egrd_M2M(4,i) + &
!                   my_system%eel%q(3,i) * my_system%eel%Egrd_M2M(5,i) +  &
!                   my_system%eel%q(4,i) * my_system%eel%Egrd_M2M(6,i)
!       end do
!
!       write(6,*) 'gradient a manina:'
!       write(6,'(10f12.8)') gg
!       gg = 0.0_rp
!       call fixedelec_geomgrad(my_system,gg)
!       write(6,*) 'field:'
!       do i = 1, n
!         write(6,'(3f12.8)') my_system%eel%e_M2M(:,i)
!       end do
!       write(6,*) 'field gradient:'
!       do i = 1, n
!         write(6,'(10f12.8)') my_system%eel%Egrd_M2M(:,i)
!       end do
!       write(6,*) 'field hessian:'
!       do i = 1, n
!         write(6,'(10f12.8)') my_system%eel%Ehes_M2M(:,i)
!       end do
!       write(6,*) 'gradient '
!       write(6,'(10f12.8)') gg
!
!fl
        allocate (cmm(3,n), gph(3,n), gmh(3,n), gg(3,n))
        cmm = my_system%top%cmm
        do i = 1, n
          do ix = 1, 3
            cmm(ix,i) = cmm(ix,i) + delta
            call update_coordinates(my_system,cmm)
            gph = 0.0_rp
            call fixedelec_geomgrad(my_system,gph)
            write(6,'(10f12.6)') gph
            cmm(ix,i) = cmm(ix,i) - 2.0_rp * delta
            call update_coordinates(my_system,cmm)
            gmh = 0.0_rp
            call fixedelec_geomgrad(my_system,gmh)
            write(6,'(10f12.6)') gmh
            cmm(ix,i) = cmm(ix,i) + delta
            call update_coordinates(my_system,cmm)
!
            do j = 1, n
              do jx = 1, 3
                hessnum(jx,ix,j,i) = (gph(jx,j) - gmh(jx,j)) / (2.0_rp * delta)
              end do
            end do
!
          end do
        end do
!
        do i = 1, n
          do j = 1, n
            write(6,*) 'diff hess for atoms (analytical/numerical) ',j, i
            write(6,'(3f14.8)') hess(1,:,j,i) - hessnum(1,:,j,i)
            write(6,'(3f14.8)') hess(2,:,j,i) - hessnum(2,:,j,i)
            write(6,'(3f14.8)') hess(3,:,j,i) - hessnum(3,:,j,i)
          end do
        end do
            
        

    end if
end program
