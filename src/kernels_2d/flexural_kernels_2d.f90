subroutine flex2d_g(src,ndt,targ,ndd,dpars,ndz,zpars,ndi,ipars,val)
  real *8 :: src(*), targ(ndt), dpars(ndd)
  real *8 :: dr 
  integer *8 ipars(ndi)
  real *8 :: dx, dy, r2
  complex *16 :: val
  complex *16 :: zpars(2), zk1, zk2 
  complex *16 :: zt1, zt2 
  complex *16 :: g1, g2, h1 
  complex *16 :: h01,h0x1,h0y1
  complex *16 :: h02,h0x2,h0y2
  complex *16 :: h0xx,h0xy,h0yy,h0xxx,h0xxy,h0xyy,h0yyy
  complex *16 :: gradx, grady


  zk1 = zpars(1)
  zk2 = zpars(2)

  dx = targ(1)-src(1)
  dy = targ(2)-src(2)

  call helmdiffgreen(zk1,dx,dy,g1,h0x1,h0y1,h0xx,h0xy,h0yy,h0xxx,h0xxy,h0xyy,h0yyy)

  if (abs(zk2).gt.1E-6) then 
    call helmdiffgreen(zk2,dx,dy,g2,h0x2,h0y2,h0xx,h0xy,h0yy,h0xxx,h0xxy,h0xyy,h0yyy)
  else 
    g2 = 0
  endif 

  val = (g1-g2)/(zk1*zk1-zk2*zk2)
!  val = g1

end subroutine flex2d_g
!
!
!
!
!
subroutine flex2d_gdn(src,ndt,targ,ndd,dpars,ndz,zpars,ndi,ipars,val)
  real *8 :: src(*), targ(ndt), dpars(ndd)
  real *8 :: dr 
  integer *8 ipars(ndi)
  real *8 :: dx, dy, r2
  real *8 :: nx, ny 
  complex *16 :: val
  complex *16 :: zpars(2), zk1, zk2 
  complex *16 :: zt1, zt2 
  complex *16 :: h01,h0x1,h0y1
  complex *16 :: h02,h0x2,h0y2
  complex *16 :: h0xx,h0xy,h0yy,h0xxx,h0xxy,h0xyy,h0yyy
  complex *16 :: gradx, grady
  

  zk1 = zpars(1)
  zk2 = zpars(2)

  dx = targ(1)-src(1)
  dy = targ(2)-src(2)

  nx = targ(10)
  ny = targ(11)

  call helmdiffgreen(zk1,dx,dy,h01,h0x1,h0y1,h0xx,h0xy,h0yy,h0xxx,h0xxy,h0xyy,h0yyy)

  if (abs(zk2).gt.1E-6) then 
    call helmdiffgreen(zk2,dx,dy,h02,h0x2,h0y2,h0xx,h0xy,h0yy,h0xxx,h0xxy,h0xyy,h0yyy)
  else 
    h02 = 0
    h0x2 = 0 
    h0y2 = 0 
  endif 

  gradx = h0x1-h0x2
  grady = h0y1-h0y2

  val = (nx*gradx + ny*grady)/(zk1*zk1-zk2*zk2)

end subroutine flex2d_gdn
!
!
!
!
!
subroutine flex2d_gsupp2(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk, &
   ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: srcinfo(*),targinfo(ndt),dpars(1)
  integer *8 ipars(ndi)
  real *8 :: dx, dy, r2, rdotn
  real *8 :: nu, nx, ny
  complex *16 :: zk(2), zk1, zk2 
  complex *16 :: val, gsxx, gsxy, gsyy
  complex *16 :: h01,h01x,h01y,h02,h02x,h02y
  complex *16 :: h01xx,h01xy,h01yy,h01xxx,h01xxy,h01xyy,h01yyy
  complex *16 :: h02xx,h02xy,h02yy,h02xxx,h02xxy,h02xyy,h02yyy


  dx = targinfo(1) - srcinfo(1)
  dy = targinfo(2) - srcinfo(2)
  r2 = dx**2 + dy**2

  zk1 = zk(1)
  zk2 = zk(2)

!  rdotn = dx*targinfo(10) + dy*targinfo(11)
  
  call helmdiffgreen(zk1,dx,dy,h01,h01x,h01y,h01xx,h01xy,h01yy,&
    h01xxx,h01xxy,h01xyy,h01yyy)
  call helmdiffgreen(zk2,dx,dy,h02,h02x,h02y,h02xx,h02xy,h02yy,&
    h02xxx,h02xxy,h02xyy,h02yyy)

  gsxx = (h01xx-h02xx)/(zk1*zk1-zk2*zk2)
  gsxy = (h01xy-h02xy)/(zk1*zk1-zk2*zk2)
  gsyy = (h01yy-h02yy)/(zk1*zk1-zk2*zk2)


  nx = targinfo(10)
  ny = targinfo(11)

  nu = dpars(1)  


  val = nu*(gsxx + gsyy) + &
    (1.0d0 - nu)*(nx*nx*gsxx + 2*nx*ny*gsxy + ny*ny*gsyy)

            
  return
end subroutine flex2d_gsupp2
!
!
!
!
!
!
!
!
subroutine flex2d_gfree2(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk, &
   ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: srcinfo(*),targinfo(ndt),dpars(1)
  integer *8 ipars(ndi)
  real *8 :: dx, dy, r2, rdotn
  real *8 :: nu, nx, ny, kappa
  real *8 :: nx2, nx3, ny2, ny3
  real *8 :: taux, tauy, taux2, tauy2
  complex *16 :: zk(2), zk1, zk2
  complex *16 :: val, gsxx, gsxy, gsyy
  complex *16 :: gsxxx, gsxxy, gsxyy, gsyyy
  complex *16 :: h01,h01x,h01y,h02,h02x,h02y
  complex *16 :: h01xx,h01xy,h01yy,h01xxx,h01xxy,h01xyy,h01yyy
  complex *16 :: h02xx,h02xy,h02yy,h02xxx,h02xxy,h02xyy,h02yyy  
  !
  ! returns the second free plate condition of the
  ! flexural volumetric kernel
  !


  dx = targinfo(1) - srcinfo(1)
  dy = targinfo(2) - srcinfo(2)

  zk1 = zk(1)
  zk2 = zk(2)

  nx = targinfo(10)
  ny = targinfo(11)
  kappa = targinfo(13)

  taux = targinfo(4)
  tauy = targinfo(5)

  ds = sqrt(taux*taux + tauy*tauy)
  taux = taux/ds
  tauy = tauy/ds

  nu = dpars(1)

  call helmdiffgreen(zk1,dx,dy,h01,h01x,h01y,h01xx,h01xy,h01yy,&
    h01xxx,h01xxy,h01xyy,h01yyy)
  call helmdiffgreen(zk2,dx,dy,h02,h02x,h02y,h02xx,h02xy,h02yy,&
    h02xxx,h02xxy,h02xyy,h02yyy)

  gsxx = (h01xx-h02xx)/(zk1*zk1-zk2*zk2)
  gsxy = (h01xy-h02xy)/(zk1*zk1-zk2*zk2)
  gsyy = (h01yy-h02yy)/(zk1*zk1-zk2*zk2)

  gsxxx = (h01xxx-h02xxx)/(zk1*zk1-zk2*zk2)
  gsxxy = (h01xxy-h02xxy)/(zk1*zk1-zk2*zk2)
  gsxyy = (h01xyy-h02xyy)/(zk1*zk1-zk2*zk2)
  gsyyy = (h01yyy-h02yyy)/(zk1*zk1-zk2*zk2)


  val = gsxxx*nx*nx*nx + 3*gsxxy*nx*nx*ny + 3*gsxyy*nx*ny*ny + &
      gsyyy*ny*ny*ny + (2-nu)*(gsxxx*nx*taux*taux + &
      gsxxy*(taux*taux*ny + 2*taux*tauy*nx) + &
      gsxyy*(2*taux*tauy*ny + tauy*tauy*nx) + &
      gsyyy*tauy*tauy*ny)+kappa*(1-nu)*(gsxx*taux*taux + &
      2*gsxy*taux*tauy+gsyy*tauy*tauy - gsxx*nx*nx - &
      2*gsxy*nx*ny - gsyy*ny*ny)

!  val = gsxx*taux*taux +2*gsxy*taux*tauy+gsyy*tauy*tauy
!  val = gsxx*nx*nx +2*gsxy*nx*ny+gsyy*ny*ny


  return
end subroutine flex2d_gfree2
!
!
subroutine flex2d_gfree2_var(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk, &
   ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: srcinfo(*),targinfo(ndt),dpars(1)
  integer *8 ipars(ndi)
  real *8 :: dx, dy, r2, rdotn
  real *8 :: nu, nx, ny, kappa, dadn, dadt, a
  real *8 :: nx2, nx3, ny2, ny3
  real *8 :: taux, tauy, taux2, tauy2
  complex *16 :: zk(2), zk1, zk2
  complex *16 :: val, gsxx, gsxy, gsyy
  complex *16 :: gsxxx, gsxxy, gsxyy, gsyyy
  complex *16 :: h01,h01x,h01y,h02,h02x,h02y
  complex *16 :: h01xx,h01xy,h01yy,h01xxx,h01xxy,h01xyy,h01yyy
  complex *16 :: h02xx,h02xy,h02yy,h02xxx,h02xxy,h02xyy,h02yyy  
  !
  ! returns the second free plate condition of the
  ! flexural volumetric kernel in the case of variable thickness
  !


  dx = targinfo(1) - srcinfo(1)
  dy = targinfo(2) - srcinfo(2)

  zk1 = zk(1)
  zk2 = zk(2)

  nx = targinfo(10)
  ny = targinfo(11)
  kappa = targinfo(13)
  a = targinfo(14)
  dadn = targinfo(15)
  dadt = targinfo(16)

  taux = targinfo(4)
  tauy = targinfo(5)

  ds = sqrt(taux*taux + tauy*tauy)
  taux = taux/ds
  tauy = tauy/ds

  nu = dpars(1)

  call helmdiffgreen(zk1,dx,dy,h01,h01x,h01y,h01xx,h01xy,h01yy,&
    h01xxx,h01xxy,h01xyy,h01yyy)
  call helmdiffgreen(zk2,dx,dy,h02,h02x,h02y,h02xx,h02xy,h02yy,&
    h02xxx,h02xxy,h02xyy,h02yyy)

  gsxx = (h01xx-h02xx)/(zk1*zk1-zk2*zk2)
  gsxy = (h01xy-h02xy)/(zk1*zk1-zk2*zk2)
  gsyy = (h01yy-h02yy)/(zk1*zk1-zk2*zk2)

  gsxxx = (h01xxx-h02xxx)/(zk1*zk1-zk2*zk2)
  gsxxy = (h01xxy-h02xxy)/(zk1*zk1-zk2*zk2)
  gsxyy = (h01xyy-h02xyy)/(zk1*zk1-zk2*zk2)
  gsyyy = (h01yyy-h02yyy)/(zk1*zk1-zk2*zk2)


  val = gsxxx*nx*nx*nx + 3*gsxxy*nx*nx*ny + 3*gsxyy*nx*ny*ny + &
      gsyyy*ny*ny*ny + (2-nu)*(gsxxx*nx*taux*taux + &
      gsxxy*(taux*taux*ny + 2*taux*tauy*nx) + &
      gsxyy*(2*taux*tauy*ny + tauy*tauy*nx) + &
      gsyyy*tauy*tauy*ny)+kappa*(1-nu)*(gsxx*taux*taux + &
      2*gsxy*taux*tauy+gsyy*tauy*tauy - gsxx*nx*nx - &
      2*gsxy*nx*ny - gsyy*ny*ny) + &
      dadn/a*(nu*(gsxx + gsyy) + &
      (1.0d0 - nu)*(nx*nx*gsxx + 2*nx*ny*gsxy + ny*ny*gsyy)) + &
      2*(1.0d0 - nu)*dadt/a*(nx*taux*gsxx + & 
      (nx*tauy+ny*taux)*gsxy + ny*tauy*gsyy)

  return
end subroutine flex2d_gfree2_var
!
!
!
!
!
subroutine flex2d_lapgx(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk, &
   ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: srcinfo(*),targinfo(ndt),dpars(1)
  integer *8 ipars(ndi)
  real *8 :: dx, dy, r2, rdotn
  real *8 :: nu, nx, ny
  complex *16 :: zk(2), zk1, zk2
  complex *16 :: val, gsxxx, gsxyy
  complex *16 :: h01,h01x,h01y,h02,h02x,h02y
  complex *16 :: h01xx,h01xy,h01yy,h01xxx,h01xxy,h01xyy,h01yyy
  complex *16 :: h02xx,h02xy,h02yy,h02xxx,h02xxy,h02xyy,h02yyy

  dx = targinfo(1) - srcinfo(1)
  dy = targinfo(2) - srcinfo(2)
  r2 = dx**2 + dy**2

  zk1 = zk(1)
  zk2 = zk(2)

  call helmdiffgreen(zk1,dx,dy,h01,h01x,h01y,h01xx,h01xy,h01yy,&
    h01xxx,h01xxy,h01xyy,h01yyy)
  call helmdiffgreen(zk2,dx,dy,h02,h02x,h02y,h02xx,h02xy,h02yy,&
    h02xxx,h02xxy,h02xyy,h02yyy)

  gsxxx = (h01xxx-h02xxx)/(zk1*zk1-zk2*zk2)
  gsxyy = (h01xyy-h02xyy)/(zk1*zk1-zk2*zk2)

  val = gsxxx + gsxyy

  return
end subroutine flex2d_lapgx
!
!
subroutine flex2d_lapgy(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk, &
   ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: srcinfo(*),targinfo(ndt),dpars(1)
  integer *8 ipars(ndi)
  real *8 :: dx, dy, r2, rdotn
  real *8 :: nu, nx, ny
  complex *16 :: zk(2), zk1, zk2
  complex *16 :: val, gsxxy, gsyyy
  complex *16 :: h01,h01x,h01y,h02,h02x,h02y
  complex *16 :: h01xx,h01xy,h01yy,h01xxx,h01xxy,h01xyy,h01yyy
  complex *16 :: h02xx,h02xy,h02yy,h02xxx,h02xxy,h02xyy,h02yyy

  dx = targinfo(1) - srcinfo(1)
  dy = targinfo(2) - srcinfo(2)
  r2 = dx**2 + dy**2

  zk1 = zk(1)
  zk2 = zk(2)

  call helmdiffgreen(zk1,dx,dy,h01,h01x,h01y,h01xx,h01xy,h01yy,&
    h01xxx,h01xxy,h01xyy,h01yyy)
  call helmdiffgreen(zk2,dx,dy,h02,h02x,h02y,h02xx,h02xy,h02yy,&
    h02xxx,h02xxy,h02xyy,h02yyy)

  gsxxy = (h01xxy-h02xxy)/(zk1*zk1-zk2*zk2)
  gsyyy = (h01yyy-h02yyy)/(zk1*zk1-zk2*zk2)

  val = gsxxy + gsyyy

  return
end subroutine flex2d_lapgy
!
!
subroutine flex2d_lapg(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk, &
   ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: srcinfo(*),targinfo(ndt),dpars(1)
  integer *8 ipars(ndi)
  real *8 :: dx, dy, r2, rdotn
  real *8 :: nu, nx, ny
  complex *16 :: zk(2), zk1, zk2
  complex *16 :: val, gsxx, gsyy
  complex *16 :: h01,h01x,h01y,h02,h02x,h02y
  complex *16 :: h01xx,h01xy,h01yy,h01xxx,h01xxy,h01xyy,h01yyy
  complex *16 :: h02xx,h02xy,h02yy,h02xxx,h02xxy,h02xyy,h02yyy

  dx = targinfo(1) - srcinfo(1)
  dy = targinfo(2) - srcinfo(2)
  r2 = dx**2 + dy**2

  zk1 = zk(1)
  zk2 = zk(2)

  call helmdiffgreen(zk1,dx,dy,h01,h01x,h01y,h01xx,h01xy,h01yy,&
    h01xxx,h01xxy,h01xyy,h01yyy)
  call helmdiffgreen(zk2,dx,dy,h02,h02x,h02y,h02xx,h02xy,h02yy,&
    h02xxx,h02xxy,h02xyy,h02yyy)

  gsxx = (h01xx-h02xx)/(zk1*zk1-zk2*zk2)
  gsyy = (h01yy-h02yy)/(zk1*zk1-zk2*zk2)

  val = gsxx + gsyy

  return
end subroutine flex2d_lapg
!
!
subroutine flex2d_gxx(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk, &
   ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: srcinfo(*),targinfo(ndt),dpars(1)
  integer *8 ipars(ndi)
  real *8 :: dx, dy, r2, rdotn
  real *8 :: nu, nx, ny
  complex *16 :: zk(2), zk1, zk2
  complex *16 :: val, gsxx
  complex *16 :: h01,h01x,h01y,h02,h02x,h02y
  complex *16 :: h01xx,h01xy,h01yy,h01xxx,h01xxy,h01xyy,h01yyy
  complex *16 :: h02xx,h02xy,h02yy,h02xxx,h02xxy,h02xyy,h02yyy

  dx = targinfo(1) - srcinfo(1)
  dy = targinfo(2) - srcinfo(2)
  r2 = dx**2 + dy**2

  zk1 = zk(1)
  zk2 = zk(2)

  call helmdiffgreen(zk1,dx,dy,h01,h01x,h01y,h01xx,h01xy,h01yy,&
    h01xxx,h01xxy,h01xyy,h01yyy)
  call helmdiffgreen(zk2,dx,dy,h02,h02x,h02y,h02xx,h02xy,h02yy,&
    h02xxx,h02xxy,h02xyy,h02yyy)

  gsxx = (h01xx-h02xx)/(zk1*zk1-zk2*zk2)

  val = gsxx

  return
end subroutine flex2d_gxx
!
!
subroutine flex2d_gyy(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk, &
   ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: srcinfo(*),targinfo(ndt),dpars(1)
  integer *8 ipars(ndi)
  real *8 :: dx, dy, r2, rdotn
  real *8 :: nu, nx, ny
  complex *16 :: zk(2), zk1, zk2
  complex *16 :: val, gsyy
  complex *16 :: h01,h01x,h01y,h02,h02x,h02y
  complex *16 :: h01xx,h01xy,h01yy,h01xxx,h01xxy,h01xyy,h01yyy
  complex *16 :: h02xx,h02xy,h02yy,h02xxx,h02xxy,h02xyy,h02yyy

  dx = targinfo(1) - srcinfo(1)
  dy = targinfo(2) - srcinfo(2)
  r2 = dx**2 + dy**2

  zk1 = zk(1)
  zk2 = zk(2)

  call helmdiffgreen(zk1,dx,dy,h01,h01x,h01y,h01xx,h01xy,h01yy,&
    h01xxx,h01xxy,h01xyy,h01yyy)
  call helmdiffgreen(zk2,dx,dy,h02,h02x,h02y,h02xx,h02xy,h02yy,&
    h02xxx,h02xxy,h02xyy,h02yyy)

  gsyy = (h01yy-h02yy)/(zk1*zk1-zk2*zk2)

  val = gsyy

  return
end subroutine flex2d_gyy
!
!
subroutine flex2d_gxy(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk, &
   ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: srcinfo(*),targinfo(ndt),dpars(1)
  integer *8 ipars(ndi)
  real *8 :: dx, dy, r2, rdotn
  real *8 :: nu, nx, ny
  complex *16 :: zk(2), zk1, zk2
  complex *16 :: val, gsxy
  complex *16 :: h01,h01x,h01y,h02,h02x,h02y
  complex *16 :: h01xx,h01xy,h01yy,h01xxx,h01xxy,h01xyy,h01yyy
  complex *16 :: h02xx,h02xy,h02yy,h02xxx,h02xxy,h02xyy,h02yyy

  dx = targinfo(1) - srcinfo(1)
  dy = targinfo(2) - srcinfo(2)
  r2 = dx**2 + dy**2

  zk1 = zk(1)
  zk2 = zk(2)

  call helmdiffgreen(zk1,dx,dy,h01,h01x,h01y,h01xx,h01xy,h01yy,&
    h01xxx,h01xxy,h01xyy,h01yyy)
  call helmdiffgreen(zk2,dx,dy,h02,h02x,h02y,h02xx,h02xy,h02yy,&
    h02xxx,h02xxy,h02xyy,h02yyy)

  gsxy = (h01xy-h02xy)/(zk1*zk1-zk2*zk2)

  val = gsxy

  return
end subroutine flex2d_gxy
!
!


