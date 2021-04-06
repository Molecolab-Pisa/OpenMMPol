subroutine rotate_multipoles
  use mmpol
  implicit none
end subroutine rotate_multipoles
!
subroutine rotation_matrix(doder,j,jz,jx,jy,r,dri,driz,drix,driy)
  use mmpol
  implicit none
!
! given an atom j and the reference atoms jx, jy, and jz, this routine
! computes the rotation matrix needed to rotate the multipoles on the
! i-th atom from the molecular frame to the lab frame.
! if required, it also return the derivative of the rotation matrices
! with respect to the coordinates of all the atoms involved in its
! definition. 
!
! this routine is completely general and can be easily augmented with 
! new rotation conventions. 
!
! given the atoms used to define the molecolar frame, identified by the
! indices jx, jy, jz, this routine builds the vectors 
!
!   xi   = cmm(:,jz) - cmm(:,i) 
!   eta  = cmm(:,jx) - cmm(:,i) 
!   zeta = cmm(:,jy) - cmm(:,i) 
!
! it then decodes the rotatoin conventions by introducing two vectors, 
! u and v, that span the xz plane. 
! this is the only convention-dependent part: everything else works
! automatically in a general way.
!
! for the definition of u and v, the unit vectors that identify the
! orthogonal molecular systems are built as follows:
!
! ez = u/|u|
! ex = gram-schmidt (v,ez)
! ey = ez x ex
!
  logical,                       intent(in)    :: doder
  integer(ip),                   intent(in)    :: j, jz, jx, jy
  real(rp),    dimension(3,3),   intent(inout) :: r
  real(rp),    dimension(3,3,3), intent(inout) :: dri, driz, drix, driy
!
  integer(ip)              :: k, l, m, n
  real(rp)                 :: xj, yj, zj, dx, dy, dz, dot
  real(rp)                 :: xi(3), eta(3), zeta(3), xi_norm, eta_norm, zeta_norm
  real(rp)                 :: u(3), v(3), u_norm, ex(3), ey(3), ez(3), ex_nrm, ex_sq,      &
                              u_sq, u_cube, e_sq, vez
  real(rp), dimension(3,3) :: ez_u, ex_v, ex_ez, ey_ez, ey_ex, u_ri, u_riz, u_rix,         &
                              u_riy, v_ri, v_rix, v_riy, ez_ri, ez_riz, ez_rix, ez_riy,    &
                              ex_ri, ex_riz, ex_rix, ex_riy, ey_ri, ey_riz, ey_rix, ey_riy
!
  real(rp),    parameter :: zofac = 0.866_rp
!
  r  = zero
!
! get the xi, eta and zeta vectors and their norms:
!
  xi          = cmm(:,jz) - cmm(:,j)
  xi_norm     = sqrt(dot_product(xi,xi))
  if (jx.ne.0) then
    eta       = cmm(:,jx) - cmm(:,j)
    eta_norm  = sqrt(dot_product(eta,eta))
  end if
  if (jy.ne.0) then
    zeta      = cmm(:,jy) - cmm(:,j)
    zeta_norm = sqrt(dot_product(zeta,zeta))
  end if
!
! build the u and v vectors, that span the xz plane according to the specific
! convention:
!
  u = zero
  v = zero
  if (mol_frame(j).eq.0) then
!
!   do nothing.
!
    r(1,1) = one
    r(2,2) = one
    r(3,3) = one
    if (doder) then
      dri  = zero
      driz = zero
      drix = zero
      driy = zero
    end if
    return
!
  else if (mol_frame(j).eq.1) then
!
!   z-then-x convention:
!
    u = xi
    v = eta
!
  else if (mol_frame(j).eq.2) then
!
!   bisector convention:
!
    u = eta_norm*xi + xi_norm*eta
    v = eta
!
  else if (mol_frame(j).eq.3) then
!
!   z-only convention:
!
    u = xi
    dot = u(3)/xi_norm
    if (dot.le.zofac) then
      v(1) = one
    else
      v(2) = one
    end if
!
  else if (mol_frame(j).eq.4) then
!
!   z-bisector convention:
!
    u = xi
    v = zeta_norm*eta + eta_norm*zeta
!
  else if (mol_frame(j).eq.5) then
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
    call fatal_error('the required rotation convention is not implemented.')
  end if
!
! build the unit vectors:
!
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
  ez_u  = zero
  ex_v  = zero
  ex_ez = zero
  ey_ez = zero
  ey_ex = zero
  u_ri  = zero
  u_riz = zero
  u_rix = zero
  u_riy = zero
  v_ri  = zero
  v_rix = zero
  v_riy = zero
!
! start by assembling the intermediates that do not depend on the sepcific convention:
!
  u_cube = u_sq*u_norm
  do k = 1, 3
    ez_u(k,k)  = one/u_norm
    ex_v(k,k)  = one/ex_nrm
    ex_ez(k,k) = - vez/ex_nrm
  end do
  do k = 1, 3
    ez_u(:,k)  = ez_u(:,k)  - u(:)*u(k)/u_cube
    ex_v(:,k)  = ex_v(:,k)  - ez(:)*ez(k)/ex_nrm - ex(:)*ex(k)/ex_nrm
    ex_ez(:,k) = ex_ez(:,k) - ex(:)*vez*v(k)/ex_sq - ez(:)*v(k)/ex_nrm
  end do
!
  ey_ez(1,1) =   zero
  ey_ez(1,2) =   ex(3)
  ey_ez(1,3) = - ex(2)
  ey_ez(2,1) = - ex(3)
  ey_ez(2,2) =   zero
  ey_ez(2,3) =   ex(1)
  ey_ez(3,1) =   ex(2)
  ey_ez(3,2) = - ex(1)
  ey_ez(3,3) =   zero
!
  ey_ez(1,1) =   zero
  ey_ez(1,2) = - ez(3)
  ey_ez(1,3) =   ez(2)
  ey_ez(2,1) =   ez(3)
  ey_ez(2,2) =   zero
  ey_ez(2,3) = - ez(1)
  ey_ez(3,1) = - ez(2)
  ey_ez(3,2) =   ez(1)
  ey_ez(3,3) =   zero
!
! now compute the convention-specific terms:
!
  if (mol_frame(j).eq.1) then
!
!   z-then-x:
!
    do k = 1, 3
      u_ri(k,k)  = - one
      u_riz(k,k) =   one
      v_ri(k,k)  = - one
      v_rix(k,k) =   one
    end do
!
  else if (mol_frame(j).eq.2) then
!
!   bisector:
!
    do k = 1, 3
      u_ri(k,k)  = - xi_norm - eta_norm
      u_riz(k,k) =   eta_norm
      u_rix(k,k) =   xi_norm
      v_ri(k,k)  = - one
      v_rix(k,k) =   one
    end do
!
    do k = 1, 3
      u_ri(:,k)  = u_ri(:,k)  - eta(:)*xi(k)/xi_norm - xi(:)*eta(k)/eta_norm
      u_riz(:,k) = u_riz(:,k) + eta(:)*xi(k)/xi_norm 
      u_rix(:,k) = u_rix(:,k) + xi(:)*eta(k)/eta_norm
    end do
!
  else if (mol_frame(j).eq.3) then
!
!   z-only:
!
    do k = 1, 3
      u_ri(k,k)  = - one
      u_riz(k,k) =   one
    end do
!
  else if (mol_frame(j).eq.4) then
!
!   z-bisector:
!
    do k = 1, 3
      u_ri(k,k)  = - one
      u_riz(k,k) =   one
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
  else if (mol_frame(j).eq.5) then
!
!   3-fold:
!
    do k = 1, 3
      u_ri(k,k)  = - xi_norm*eta_norm - xi_norm*zeta_norm - eta_norm*zeta_norm
      u_riz(k,k) =   eta_norm*zeta_norm
      u_rix(k,k) =   xi_norm*zeta_norm
      u_riy(k,k) =   xi_norm*eta_norm
      v_ri(k,k)  = - one
      v_rix(k,k) =   one
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
  ey_riz = matmul(ey_ex,ex_riy) + matmul(ey_ez,ez_riy)
!
! finally, assemble the derivatives of the rotation matrices:
!
  dri(:,:,1)  = ex_ri
  dri(:,:,2)  = ey_ri
  dri(:,:,3)  = ez_ri
  driz(:,:,1) = ex_riz
  driz(:,:,2) = ey_riz
  driz(:,:,3) = ez_riz
  drix(:,:,1) = ex_rix
  drix(:,:,2) = ey_rix
  drix(:,:,3) = ez_rix
  driy(:,:,1) = ex_riy
  driy(:,:,2) = ey_riy
  driy(:,:,3) = ez_riy
!
  return
end subroutine rotation_matrix
