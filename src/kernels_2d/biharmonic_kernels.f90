! Kernels for solving Biharmonic equation in 2D
!
! the following routines rely on the srcinfo and targinfo arrays
! containing the following information, standardized in the following
! way:
!
!     *info(1:3) = xyz
!     *info(4:6) = tanget vector 1
!     *info(7:9) = tangent vector 2
!     *info(10:12) = normal vector
!     *info(13) = curvature
!
!

subroutine bh2d_g(src, ndt,targ, ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: src(*), targ(ndt),dpars(ndd)
  integer ipars(ndi)
  real *8 :: dx, dy, r2
  real *8 :: over4pi
  real *8 :: val
  data over4pi/0.07957747154594767d0/
  !
  ! returns the biharmonic volumetric potential kernel
  !

  dx=targ(1)-src(1)
  dy=targ(2)-src(2)

  r2=dx**2+dy**2

  !val = -over4pi*log(r2)
  val = r2 * log(r2) * over4pi / 4  

  return
end subroutine bh2d_g 
!
!
!
!
!
subroutine bh2d_gdn(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk, &
   ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(*),targinfo(ndt),dpars(ndd)
  integer ipars(ndi)
  real *8 :: dx, dy, r2, rdotn
  real *8 :: val
  real *8 :: over4pi
  data over4pi/0.07957747154594767d0/
  !
  ! returns the normal derivative of the biharmonic volumetric kernel
  !

  dx = targinfo(1) - srcinfo(1)
  dy = targinfo(2) - srcinfo(2)

  rdotn = dx*targinfo(10) + dy*targinfo(11)
  r2 = dx**2 + dy**2

  !val =  -2*rdotn/r2*over4pi
  val = rdotn*(log(r2)+1)*over4pi/2 

  return
end subroutine bh2d_gdn
!
!
!
!
!
subroutine bh2d_gsupp2(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk, &
   ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(*),targinfo(ndt),dpars(ndd)
  integer ipars(ndi)
  real *8 :: dx, dy, r2, rdotn
  real *8 :: val, nu, gsxx, gsxy, gsyy
  real *8 :: over4pi
  real *8 :: nx, ny
  data over4pi/0.07957747154594767d0/

  dx = targinfo(1) - srcinfo(1)
  dy = targinfo(2) - srcinfo(2)


  rdotn = dx*targinfo(10) + dy*targinfo(11)
  r2 = dx**2 + dy**2

  call bh2d_green_der2(dx,dy,gsxx,gsxy,gsyy)

  nx = targinfo(10)
  ny = targinfo(11)

  nu = dpars(1)  


  val = nu*(gsxx + gsyy) + &
    (1.0d0 - nu)*(nx*nx*gsxx + 2*nx*ny*gsxy + ny*ny*gsyy)

            
  return
end subroutine bh2d_gsupp2
!
!
!
!
!
subroutine bh2d_gfree2(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk, &
   ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(*),targinfo(ndt),dpars(ndd)
  integer ipars(ndi)
  real *8 :: dx, dy, r2, rdotn
  real *8 :: val, nu, gsxx, gsxy ,gsyy
  real *8 :: over4pi
  real *8 :: taux, tauy, taux2, tauy2
  real *8 :: nx, ny, nx2, ny2, nx3, ny3 
  real *8 :: kappa 
  data over4pi/0.07957747154594767d0/
  !
  ! returns the second free plate condition of the 
  ! biharmonic volumetric kernel
  !

  
  dx = targinfo(1) - srcinfo(1)
  dy = targinfo(2) - srcinfo(2)


  nx = targinfo(10)
  ny = targinfo(11)
  kappa = targinfo(13)

  nx2 = nx*nx 
  nx3 = nx*nx2 

  ny2 = ny*ny 
  ny3 = ny*ny2 

  taux = targinfo(4)
  tauy = targinfo(5)

  taux2 = taux*taux 
  tauy2 = tauy*tauy

  nu = dpars(1)

  call bh2d_green_der23(dx,dy,gsxx,gsxy,gsyy,gsxxx,gsxxy,gsxyy,&
  gsyyy)

  val = gsxxx*nx*nx*nx + 3*gsxxy*nx*nx*ny + 3*gsxyy*nx*ny*ny + & 
      gsyyy*ny*ny*ny + (2-nu)*(gsxxx*nx*taux*taux + & 
      gsxxy*(taux*taux*ny + 2*taux*tauy*nx) + & 
      gsxyy*(2*taux*tauy*ny + tauy*tauy*nx) + & 
      gsyyy*tauy*tauy*ny)+kappa*(1-nu)*(gsxx*taux*taux + & 
      2*gsxy*taux*tauy+gsyy*tauy*tauy - gsxx*nx*nx - &
      2*gsxy*nx*ny - gsyy*ny*ny)


  return
end subroutine bh2d_gfree2
!
!
!
!
!
subroutine bh2d_green_der2(dx,dy,gsxx,gsxy,gsyy)
  implicit real *8 (a-h,o-z)
  real *8 :: gsxx, gsxy, gsyy
  real *8 :: dx, dy
  real *8 :: dx2, dy2, r2, r 
  real *8 :: rm1, rm2 
  real *8 :: g0, g1, g21, g321, g4321, g54321

  dx2 = dx*dx 
  dy2 = dy*dy 
  r2 = dx2+dy2 
  r = dsqrt(r2)
  call r2logr_rders(r,g0,g1,g21,g321,g4321,g54321)

  rm1 = 1.0d0/r 
  rm2 = rm1*rm1
  
! second-order derivative 
  gsxx = dx2*g21*rm2+g1*rm1
  gsxy = dx*dy*g21*rm2
  gsyy = dy2*g21*rm2+g1*rm1


  return
end subroutine bh2d_green_der2
!
!
!
!
!
subroutine bh2d_green_der23(dx,dy,gsxx,gsxy,gsyy,gsxxx,gsxxy,gsxyy,&
  gsyyy)
  implicit real *8 (a-h,o-z)
  real *8 :: gsxx, gsxy, gsyy, gsxxx, gsxxy, gsxyy, gsyyy
  real *8 :: dx, dy, dx2, dy2, r2, r
  real *8 :: rm1, rm2, rm3
  real *8 :: g0, g1, g21, g321, g4321, g54321

  ! integer, save :: unit_r = -1
  ! logical, save :: first = .true.

  dx2 = dx*dx  
  dy2 = dy*dy 
  dx3 = dx2*dx 
  dy3 = dy2*dy  
  r2 = dx2+dy2 
  r = dsqrt(r2)

  
  ! if (first) then ! open file once
  !   open(newunit=unit_r, file='r_values.dat', &
  !        status='replace', action='write')
  !   first = .false.
  ! end if
  ! write(unit_r,'(ES25.16)') r   ! save r


  call r2logr_rders(r,g0,g1,g21,g321,g4321,g54321)

  rm1 = 1.0d0/r 
  rm2 = rm1*rm1 
  rm3 = rm1*rm2


! second-order derivative 
  gsxx = dx2*g21*rm2+g1*rm1
  gsxy = dx*dy*g21*rm2
  gsyy = dy2*g21*rm2+g1*rm1


! third-order derivative 
  gsxxx = dx3*g321*rm3 + 3*dx*g21*rm2
  gsxxy = dx2*dy*g321*rm3 + dy*g21*rm2
  gsxyy = dx*dy2*g321*rm3 + dx*g21*rm2
  gsyyy = dy3*g321*rm3 + 3*dy*g21*rm2


  return
end subroutine bh2d_green_der23
!
!
!
!
!
subroutine r2logr_rders(r,g0,g1,g21,g321,g4321,g54321)
  implicit real *8 (a-h,o-z)
  real *8 :: r, g0, g1, g21, g321, g4321, g54321
  real *8 :: o8p 
  real *8 :: r2, r2d1, rm1, rm2, rm3, rm4, rm5 
  real *8 :: logr 
  data over4pi/0.07957747154594767d0/

  
  o8p = over4pi/2.0d0
  ! call prin2_long('o8p=*',o8p,1)
  r2 = r*r;
  r2d1 = 2*r;

  rm1 = 1.0d0/r;
  rm2 = rm1*rm1;
  rm3 = rm2*rm1;
  rm4 = rm3*rm1;
  rm5 = rm4*rm1;

  logr = log(r);

  g0 = o8p*(logr*r2);
  g1 = o8p*(logr*r2d1 + r2*rm1);
  g21 = o8p*(2*r2d1*rm1 - 2*r2*rm2);
  g321 = o8p*(-6*r2d1*rm2 + 8*r2*rm3);
  g4321 = o8p*(32*r2d1*rm3 - 48*r2*rm4);
  g54321 = o8p*(- 240*r2d1*rm4 + 384*r2*rm5);

  ! call prin2_long('g0=*',g0,1)
  ! call prin2_long('g1=*',g1,1)
  ! call prin2_long('g21=*',g21,1)
  ! call prin2_long('g321=*',g321,1)
  ! call prin2_long('g4321=*',g4321,1)
  ! call prin2_long('g54321=*',g54321,1)

  return 
end subroutine r2logr_rders