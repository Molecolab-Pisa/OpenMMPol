#include "f_cart_components.h"

subroutine rotation_geomgrad(eel, E, Egrd, grad)
    
    use mod_memory, only: ip, rp
    use mod_electrostatics, only: ommp_electrostatics_type
    
    implicit none
  
    type(ommp_electrostatics_type), intent(in) :: eel
    real(rp), intent(in) :: E(3, eel%top%mm_atoms), Egrd(6, eel%top%mm_atoms)
    real(rp), dimension(3, eel%top%mm_atoms), intent(inout) :: grad

    integer(ip) :: j, jx, jy, jz, k, l, m, n
    real(rp), dimension(3) :: dip
    real(rp), dimension(3,3) :: r, rt, qua, rqua, ddip
    real(rp), dimension(3,3,3) :: dri, driz, drix, driy, dqua, dtmp

    ! TODO prepare_fixedelec

    ! loop over the mm sites and build the derivatives of the rotation
    ! matrices with respect to the positions of all the relevant atoms.
    do j = 1, eel%top%mm_atoms
        jz = eel%iz(j)
        if(jz == 0) jz = j
        jx = eel%ix(j)
        if(jx == 0) jx = j
        jy = eel%iy(j)
        if(jy == 0) jy = j
        
        call rotation_matrix(.true., & 
                             eel%top%cmm(:,j), eel%top%cmm(:,jx), &
                             eel%top%cmm(:,jy), eel%top%cmm(:,jz), &
                             eel%mol_frame(j), &
                             r, dri, driz, drix, driy)
        
        rt = transpose(r)
        
        ! get the multipoles. we also need the rotated quadrupoles:
        dip(_x_) = eel%q0(1+_x_,j)
        dip(_y_) = eel%q0(1+_y_,j)
        dip(_z_) = eel%q0(1+_z_,j)

        qua(_x_,_x_)  = eel%q0(4+_xx_,j)
        qua(_x_,_y_)  = eel%q0(4+_xy_,j)
        qua(_x_,_z_)  = eel%q0(4+_xz_,j)
        qua(_y_,_x_)  = eel%q0(4+_yx_,j)
        qua(_y_,_y_)  = eel%q0(4+_yy_,j)
        qua(_y_,_z_)  = eel%q0(4+_yz_,j)
        qua(_z_,_x_)  = eel%q0(4+_zx_,j)
        qua(_z_,_y_)  = eel%q0(4+_zy_,j)
        qua(_z_,_z_)  = eel%q0(4+_zz_,j)
        
        rqua(_x_,_x_)  = eel%q(4+_xx_,j)
        rqua(_x_,_y_)  = eel%q(4+_xy_,j)
        rqua(_x_,_z_)  = eel%q(4+_xz_,j)
        rqua(_y_,_x_)  = eel%q(4+_yx_,j)
        rqua(_y_,_y_)  = eel%q(4+_yy_,j)
        rqua(_y_,_z_)  = eel%q(4+_yz_,j)
        rqua(_z_,_x_)  = eel%q(4+_zx_,j)
        rqua(_z_,_y_)  = eel%q(4+_zy_,j)
        rqua(_z_,_z_)  = eel%q(4+_zz_,j)
        
        ! contributions to the forces on the j-th atoms:

        ddip = 0.0_rp
        dqua = 0.0_rp
        dtmp = 0.0_rp
        ! compute the differentiated multipoles:
        do k = 1, 3
            do l = 1, 3
                ddip(:,k) = ddip(:,k) + dri(:,k,l)*dip(l)
            end do
        end do
        
        do k = 1, 3
            do l = 1, 3
                do m = 1, 3
                    do n = 1, 3
                        !dtmp(:,l,m)*rt(m,k)
                        dqua(:,k,l) = dqua(:,k,l) + (dri(:,k,m)*qua(m,n) - &
                                      rqua(k,m)*dri(:,m,n))*rt(n,l)  
                    end do
                end do
            end do
        end do

        ! increment the forces for the dipoles...
        grad(:,j) = grad(:,j) - ddip(:,_x_)*E(_x_,j) &
                              - ddip(:,_y_)*E(_y_,j) & 
                              - ddip(:,_z_)*E(_z_,j)
        ! ... and for the quadrupoles:
        grad(:,j) = grad(:,j) + dqua(:,_x_,_x_)*Egrd(_xx_,j) &
                              + dqua(:,_y_,_y_)*Egrd(_yy_,j) &
                              + dqua(:,_z_,_z_)*Egrd(_zz_,j) &
                              + 2*(dqua(:,_x_,_y_)*Egrd(_xy_,j) &
                              +    dqua(:,_x_,_z_)*Egrd(_xz_,j) &
                              +    dqua(:,_y_,_z_)*Egrd(_yz_,j))
     
        ! do jx
        ddip = 0.0_rp
        dqua = 0.0_rp
        dtmp = 0.0_rp
        do k = 1, 3
            do l = 1, 3
                ddip(:,k) = ddip(:,k) + drix(:,k,l)*dip(l)
            end do
        end do
        
        do k = 1, 3
            do l = 1, 3
                do m = 1, 3
                    do n = 1, 3
                        !dtmp(:,l,m)*rt(m,k)
                        dqua(:,k,l) = dqua(:,k,l) + (drix(:,k,m)*qua(m,n) -&
                                      rqua(k,m)*drix(:,m,n))*rt(n,l)
                    end do
                end do
            end do
        end do

        ! increment the forces for the dipoles...
        grad(:,jx) = grad(:,jx) - ddip(:,_x_)*E(_x_,j) &
                                - ddip(:,_y_)*E(_y_,j) & 
                                - ddip(:,_z_)*E(_z_,j)
        ! ... and for the quadrupoles:
        grad(:,jx) = grad(:,jx) + dqua(:,_x_,_x_)*Egrd(_xx_,j) &
                                + dqua(:,_y_,_y_)*Egrd(_yy_,j) &
                                + dqua(:,_z_,_z_)*Egrd(_zz_,j) &
                                + 2*(dqua(:,_x_,_y_)*Egrd(_xy_,j) &
                                +    dqua(:,_x_,_z_)*Egrd(_xz_,j) &
                                +    dqua(:,_y_,_z_)*Egrd(_yz_,j))
                
        ! do jy
        ddip = 0.0_rp
        dqua = 0.0_rp
        dtmp = 0.0_rp
        do k = 1, 3
            do l = 1, 3
                ddip(:,k) = ddip(:,k) + driy(:,k,l)*dip(l)
            end do
        end do
        
        do k = 1, 3
            do l = 1, 3
                do m = 1, 3
                    do n = 1, 3
                        ! dtmp(:,l,m)*rt(m,k)
                        dqua(:,k,l) = dqua(:,k,l) + (driy(:,k,m)*qua(m,n) -&
                                      rqua(k,m)*driy(:,m,n))*rt(n,l)
                    end do
                end do
            end do
        end do
        
        ! increment the forces for the dipoles...
        grad(:,jy) = grad(:,jy) - ddip(:,_x_)*E(_x_,j) &
                                - ddip(:,_y_)*E(_y_,j) & 
                                - ddip(:,_z_)*E(_z_,j)
        ! ... and for the quadrupoles:
        grad(:,jy) = grad(:,jy) + dqua(:,_x_,_x_)*Egrd(_xx_,j) &
                                + dqua(:,_y_,_y_)*Egrd(_yy_,j) &
                                + dqua(:,_z_,_z_)*Egrd(_zz_,j) &
                                + 2*(dqua(:,_x_,_y_)*Egrd(_xy_,j) &
                                +    dqua(:,_x_,_z_)*Egrd(_xz_,j) &
                                +    dqua(:,_y_,_z_)*Egrd(_yz_,j))
         
        ! Do jz
        ddip = 0.0_rp
        dqua = 0.0_rp
        dtmp = 0.0_rp
        do k = 1, 3
            do l = 1, 3
                ddip(:,k) = ddip(:,k) + driz(:,k,l)*dip(l)
            end do
        end do
        
        do k = 1, 3
            do l = 1, 3
                do m = 1, 3
                    do n = 1, 3
                        ! tmp(:,l,m)*rt(m,k)
                        dqua(:,k,l) = dqua(:,k,l) + (driz(:,k,m)*qua(m,n) -&
                                      rqua(k,m)*driz(:,m,n))*rt(n,l)
                    end do
                end do
            end do
        end do
        
        ! increment the forces for the dipoles...
        grad(:,jz) = grad(:,jz) - ddip(:,_x_)*E(_x_,j) &
                                - ddip(:,_y_)*E(_y_,j) & 
                                - ddip(:,_z_)*E(_z_,j)
        ! ... and for the quadrupoles:
        grad(:,jz) = grad(:,jz) + dqua(:,_x_,_x_)*Egrd(_xx_,j) &
                                + dqua(:,_y_,_y_)*Egrd(_yy_,j) &
                                + dqua(:,_z_,_z_)*Egrd(_zz_,j) &
                                + 2*(dqua(:,_x_,_y_)*Egrd(_xy_,j) &
                                +    dqua(:,_x_,_z_)*Egrd(_xz_,j) &
                                +    dqua(:,_y_,_z_)*Egrd(_yz_,j))

    end do 
end subroutine rotation_geomgrad

subroutine rotate_multipoles(eel)
    !! this routine rotates the atomic multipoles from the molecular frame
    !! where they are defined as force field parameters to the lab frame.
    !! if required, it also computes the contribution to the forces that
    !! stems from the derivatives of the rotation matrices, sometimes 
    !! referred to as "torques".
    !! for the latter task, it uses the field and field gradient from 
    !! at the multipoles, which is passed in def. 
    !! for consistency reasons, def is dimensioned (ld_cder,mm_atoms). 
    
    use mod_memory, only: ip, rp
    use mod_electrostatics, only: ommp_electrostatics_type
    
    implicit none
  
    type(ommp_electrostatics_type), intent(inout) :: eel

    integer(ip) :: j, jx, jy, jz 
    real(rp), dimension(3,3) :: r, rt, qua, rqua, tmp
    real(rp), dimension(3,3,3) :: dri, driz, drix, driy

    ! loop over the mm sites and build the rotation matrices.
    do j = 1, eel%top%mm_atoms
        jz = eel%iz(j)
        if(jz == 0) jz = j
        jx = eel%ix(j)
        if(jx == 0) jx = j
        jy = eel%iy(j)
        if(jy == 0) jy = j
        
        call rotation_matrix(.false., & 
                             eel%top%cmm(:,j), eel%top%cmm(:,jx), &
                             eel%top%cmm(:,jy), eel%top%cmm(:,jz), &
                             eel%mol_frame(j), &
                             r, dri, driz, drix, driy)
        
        rt = transpose(r)
        ! copy the monopole:
        eel%q(1,j) = eel%q0(1,j)

        ! rotate the dipole 
        eel%q(2:4,j) = matmul(r,eel%q0(2:4,j))

        ! exctract, rotate and put back the quadrupole:
        qua(_x_,_x_)  = eel%q0(4+_xx_,j)
        qua(_x_,_y_)  = eel%q0(4+_xy_,j)
        qua(_x_,_z_)  = eel%q0(4+_xz_,j)
        qua(_y_,_x_)  = eel%q0(4+_yx_,j)
        qua(_y_,_y_)  = eel%q0(4+_yy_,j)
        qua(_y_,_z_)  = eel%q0(4+_yz_,j)
        qua(_z_,_x_)  = eel%q0(4+_zx_,j)
        qua(_z_,_y_)  = eel%q0(4+_zy_,j)
        qua(_z_,_z_)  = eel%q0(4+_zz_,j)
        
        tmp  = matmul(r,qua)
        rqua = matmul(tmp,rt)
        eel%q(4+_xx_,j)  = rqua(_x_,_x_)
        eel%q(4+_yy_,j)  = rqua(_y_,_y_)
        eel%q(4+_zz_,j)  = rqua(_z_,_z_)
        eel%q(4+_xy_,j)  = rqua(_x_,_y_)
        eel%q(4+_xz_,j)  = rqua(_x_,_z_)
        eel%q(4+_yz_,j)  = rqua(_y_,_z_)
    end do
end subroutine rotate_multipoles


subroutine rotation_matrix(doder,c,cx,cy,cz,mol_frame,r,dri,driz,drix,driy)
    use mod_io, only: fatal_error
    use mod_memory, only: ip, rp
    use mod_constants, only: eps_rp
    
    implicit none
    !! given an atom j and the reference atoms jx, jy, and jz, this routine
    !! computes the rotation matrix needed to rotate the multipoles on the
    !! i-th atom from the molecular frame to the lab frame.
    !! if required, it also return the derivative of the rotation matrices
    !! with respect to the coordinates of all the atoms involved in its
    !! definition. 
    !!
    !! this routine is completely general and can be easily augmented with 
    !! new rotation conventions. 
    !!
    !! given the atoms used to define the molecolar frame, identified by the
    !! indices jx, jy, jz, this routine builds the vectors 
    !!
    !!   xi   = cmm(:,jz) - cmm(:,i) 
    !!   eta  = cmm(:,jx) - cmm(:,i) 
    !!   zeta = cmm(:,jy) - cmm(:,i) 
    !!
    !! it then decodes the rotatoin conventions by introducing two vectors, 
    !! u and v, that span the xz plane. 
    !! this is the only convention-dependent part: everything else works
    !! automatically in a general way.
    !!
    !! for the definition of u and v, the unit vectors that identify the
    !! orthogonal molecular systems are built as follows:
    !!
    !! ez = u/|u|
    !! ex = gram-schmidt (v,ez)
    !! ey = ez x ex
    !!
    !! output:
    !! =======
    !!
    !! r(i,j) is the rotation matrix, whose columns are (ex,ey,ez)
    !!
    !! dri(i,j,k) contains the derivative of the i-th component of e_k 
    !! with respect to ri_j
  
    logical, intent(in)    :: doder
    real(rp), intent(in), dimension(3) :: c, cx, cy, cz
    integer(ip), intent(in) :: mol_frame

    real(rp), dimension(3,3), intent(out) :: r
    real(rp), dimension(3,3,3), intent(inout) :: dri, driz, drix, driy

    integer(ip) :: k
    real(rp) :: dot
    real(rp) :: xi(3), eta(3), zeta(3), xi_norm, eta_norm, zeta_norm
    real(rp) :: u(3), v(3), u_norm, ex(3), ey(3), ez(3), ex_nrm, ex_sq,      &
                u_sq, u_cube, vez
    real(rp), dimension(3,3) :: ez_u, ex_v, ex_ez, ey_ez, ey_ex, u_ri, u_riz, &
                                u_rix, u_riy, v_ri, v_rix, v_riy, ez_ri, &
                                ez_riz, ez_rix, ez_riy, ex_ri, ex_riz, &
                                ex_rix, ex_riy, ey_ri, ey_riz, ey_rix, ey_riy
    
    real(rp), parameter :: zofac = 0.866_rp
    
    r = 0.0_rp
    
    ! get the xi, eta and zeta vectors and their norms:
    xi = cz - c
    xi_norm  = sqrt(dot_product(xi,xi))

    eta = cx - c
    eta_norm = sqrt(dot_product(eta,eta))

    zeta = cy-c
    zeta_norm = sqrt(dot_product(zeta,zeta))
!
! build the u and v vectors, that span the xz plane according to the specific
! convention:
!
  u = 0.0_rp
  v = 0.0_rp
  if (mol_frame.eq.0) then
!
!   do nothing.
!
    r(1,1) = 1.0_rp
    r(2,2) = 1.0_rp
    r(3,3) = 1.0_rp
    if (doder) then
      dri  = 0.0_rp
      driz = 0.0_rp
      drix = 0.0_rp
      driy = 0.0_rp
    end if
    return
!
  else if (mol_frame.eq.1) then
!
!   z-then-x convention:
!
    u = xi
    v = eta
!
  else if (mol_frame.eq.2) then
!
!   bisector convention:
!
    u = eta_norm*xi + xi_norm*eta
    v = eta
!
  else if (mol_frame.eq.3) then
!
!   z-only convention:
!
    u = xi
    dot = u(3)/xi_norm
    if (dot.le.zofac) then
      v(1) = 1.0_rp
    else
      v(2) = 1.0_rp
    end if
!
  else if (mol_frame.eq.4) then
!
!   z-bisector convention:
!
    u = xi
    v = zeta_norm*eta + eta_norm*zeta
!
  else if (mol_frame.eq.5) then
!
!   3-fold convention:
!
    u = eta_norm*zeta_norm*xi + xi_norm*zeta_norm*eta + xi_norm*eta_norm*zeta
    v = eta
!
  else
!
!   convention not yet implemented.
!
    !!TODO call fatal_error('the required rotation convention is not implemented.')
  end if
!
! build the unit vectors:
!
  if(norm2(u) < eps_rp .or. norm2(v) < eps_rp) then
      call fatal_error("Ill defined rotation axis (either u or v vectors have &
                       &0 norm.")
  end if

  u_sq   = dot_product(u,u)
  u_norm = sqrt(u_sq)
  ez     = u/u_norm
!
  vez    = dot_product(v,ez)
  ex     = v - vez*ez
  ex_sq  = dot_product(ex,ex)
  ex_nrm = sqrt(ex_sq)
  ex     = ex/ex_nrm
!
  ey(1)  = ez(2)*ex(3) - ez(3)*ex(2)
  ey(2)  = ez(3)*ex(1) - ez(1)*ex(3)
  ey(3)  = ez(1)*ex(2) - ez(2)*ex(1)
!
  r(:,1) = ex
  r(:,2) = ey
  r(:,3) = ez
!
! return if derivatives are not requested:
!
  if (.not. doder) return
!
! derivatives code
! ================
!
! the derivatives of the rotation matrices can be computed using chain-rule
! differentiation. 
! please, see JCTC 2014, 10, 1638−1651 for details.
!
! clean interemediates:
!
  ez_u  = 0.0_rp
  ex_v  = 0.0_rp
  ex_ez = 0.0_rp
  ey_ez = 0.0_rp
  ey_ex = 0.0_rp
  u_ri  = 0.0_rp
  u_riz = 0.0_rp
  u_rix = 0.0_rp
  u_riy = 0.0_rp
  v_ri  = 0.0_rp
  v_rix = 0.0_rp
  v_riy = 0.0_rp
!
! start by assembling the intermediates that do not depend on the sepcific convention:
!
  u_cube = u_sq*u_norm
  do k = 1, 3
    ez_u(k,k)  = 1.0_rp/u_norm
    ex_v(k,k)  = 1.0_rp/ex_nrm
    ex_ez(k,k) = - vez/ex_nrm
  end do
  do k = 1, 3
    ez_u(:,k)  = ez_u(:,k)  - u(:)*u(k)/u_cube
    ex_v(:,k)  = ex_v(:,k)  - ez(:)*ez(k)/ex_nrm - ex(:)*ex(k)/ex_nrm
    ex_ez(:,k) = ex_ez(:,k) + ex(:)*vez*v(k)/ex_sq - ez(:)*v(k)/ex_nrm
  end do
!
  ey_ez(1,1) =   0.0_rp
  ey_ez(1,2) =   ex(3)
  ey_ez(1,3) = - ex(2)
  ey_ez(2,1) = - ex(3)
  ey_ez(2,2) =   0.0_rp
  ey_ez(2,3) =   ex(1)
  ey_ez(3,1) =   ex(2)
  ey_ez(3,2) = - ex(1)
  ey_ez(3,3) =   0.0_rp
!
  ey_ex(1,1) =   0.0_rp
  ey_ex(1,2) = - ez(3)
  ey_ex(1,3) =   ez(2)
  ey_ex(2,1) =   ez(3)
  ey_ex(2,2) =   0.0_rp
  ey_ex(2,3) = - ez(1)
  ey_ex(3,1) = - ez(2)
  ey_ex(3,2) =   ez(1)
  ey_ex(3,3) =   0.0_rp
!
! now compute the convention-specific terms:
!
  if (mol_frame.eq.1) then
!
!   z-then-x:
!
    do k = 1, 3
      u_ri(k,k)  = - 1.0_rp
      u_riz(k,k) =   1.0_rp
      v_ri(k,k)  = - 1.0_rp
      v_rix(k,k) =   1.0_rp
    end do
!
  else if (mol_frame.eq.2) then
!
!   bisector:
!
    do k = 1, 3
      u_ri(k,k)  = - xi_norm - eta_norm
      u_riz(k,k) =   eta_norm
      u_rix(k,k) =   xi_norm
      v_ri(k,k)  = - 1.0_rp
      v_rix(k,k) =   1.0_rp
    end do
!
    do k = 1, 3
      u_ri(:,k)  = u_ri(:,k)  - eta(:)*xi(k)/xi_norm - xi(:)*eta(k)/eta_norm
      u_riz(:,k) = u_riz(:,k) + eta(:)*xi(k)/xi_norm 
      u_rix(:,k) = u_rix(:,k) + xi(:)*eta(k)/eta_norm
    end do
!
  else if (mol_frame.eq.3) then
!
!   z-only:
!
    do k = 1, 3
      u_ri(k,k)  = - 1.0_rp
      u_riz(k,k) =   1.0_rp
    end do
!
  else if (mol_frame.eq.4) then
!
!   z-bisector:
!
    do k = 1, 3
      u_ri(k,k)  = - 1.0_rp
      u_riz(k,k) =   1.0_rp
      v_ri(k,k)  = - eta_norm - zeta_norm
      v_rix(k,k) =   zeta_norm
      u_riy(k,k) =   eta_norm
    end do
!
    do k = 1, 3
      v_ri(:,k)  = v_ri(:,k)  - eta(:)*zeta(k)/zeta_norm - zeta(:)*eta(k)/eta_norm
      v_rix(:,k) = v_rix(:,k) + zeta(:)*eta(k)/eta_norm
      v_riy(:,k) = v_riy(:,k) + eta(:)*zeta(k)/zeta_norm
    end do
!
  else if (mol_frame.eq.5) then
!
!   3-fold:
!
    do k = 1, 3
      u_ri(k,k)  = - xi_norm*eta_norm - xi_norm*zeta_norm - eta_norm*zeta_norm
      u_riz(k,k) =   eta_norm*zeta_norm
      u_rix(k,k) =   xi_norm*zeta_norm
      u_riy(k,k) =   xi_norm*eta_norm
      v_ri(k,k)  = - 1.0_rp
      v_rix(k,k) =   1.0_rp
    end do
!
    do k = 1, 3
      u_ri(:,k)  = u_ri(:,k)  - (eta*zeta_norm + zeta*eta_norm)*xi(k)/xi_norm &
                              - (xi*zeta_norm  + zeta*xi_norm )*eta(k)/eta_norm &
                              - (xi*eta_norm   + eta*xi_norm  )*zeta(k)/zeta_norm
      u_riz(:,k) = u_riz(:,k) + (eta*zeta_norm + zeta*eta_norm)*xi(k)/xi_norm 
      u_rix(:,k) = u_rix(:,k) + (xi*zeta_norm  + zeta*xi_norm )*eta(k)/eta_norm 
      u_riy(:,k) = u_riy(:,k) + (xi*eta_norm   + eta*xi_norm  )*zeta(k)/zeta_norm
    end do
!
  end if
!
! proceed assembling the derivatives of the unit vectors:
!
! ez:
!
  ez_ri  = matmul(ez_u,u_ri)
  ez_riz = matmul(ez_u,u_riz)
  ez_rix = matmul(ez_u,u_rix)
  ez_riy = matmul(ez_u,u_riy)
!
! ex:
!
  ex_ri  = matmul(ex_v,v_ri)  + matmul(ex_ez,ez_ri)
  ex_riz =                      matmul(ex_ez,ez_riz)
  ex_rix = matmul(ex_v,v_rix) + matmul(ex_ez,ez_rix)
  ex_riy = matmul(ex_v,v_riy) + matmul(ex_ez,ez_riy)
!
! ey:
!
  ey_ri  = matmul(ey_ex,ex_ri)  + matmul(ey_ez,ez_ri)
  ey_riz = matmul(ey_ex,ex_riz) + matmul(ey_ez,ez_riz)
  ey_rix = matmul(ey_ex,ex_rix) + matmul(ey_ez,ez_rix)
  ey_riy = matmul(ey_ex,ex_riy) + matmul(ey_ez,ez_riy)
!
! finally, assemble the derivatives of the rotation matrices:
!
  dri(:,:,1)  = transpose(ex_ri)
  dri(:,:,2)  = transpose(ey_ri)
  dri(:,:,3)  = transpose(ez_ri)
  driz(:,:,1) = transpose(ex_riz)
  driz(:,:,2) = transpose(ey_riz)
  driz(:,:,3) = transpose(ez_riz)
  drix(:,:,1) = transpose(ex_rix)
  drix(:,:,2) = transpose(ey_rix)
  drix(:,:,3) = transpose(ez_rix)
  driy(:,:,1) = transpose(ex_riy)
  driy(:,:,2) = transpose(ey_riy)
  driy(:,:,3) = transpose(ez_riy)
end subroutine rotation_matrix
