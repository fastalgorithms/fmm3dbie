subroutine flex2d_gders(zpars,dx,dy,g,gx,gy,gxx,gxy,gyy, &
   gxxx,gxxy,gxyy,gyyy)
  !
  ! returns the flexural Green's function and its target derivatives up
  ! to third order, for the pair of wavenumbers zpars = (zk1, zk2):
  !
  !   both nonzero:  G = (G_{zk1} - G_{zk2})/(zk1^2 - zk2^2)
  !   one zero:      G = (G_{zk} - G_0)/zk^2, zk the nonzero wavenumber
  !   both zero:     G = r^2 log r/(8 pi), the biharmonic Green's function
  !
  ! where G_k = i/4 H_0^(1)(k r) and G_0 = -log r/(2 pi). A wavenumber is
  ! treated as zero if its absolute value is below 1e-6.
  !
  implicit none
  complex *16 :: zpars(2), zk1, zk2, zfac
  real *8 :: dx, dy, s, s2, rl
  complex *16 :: g, gx, gy, gxx, gxy, gyy, gxxx, gxxy, gxyy, gyyy
  complex *16 :: h, hx, hy, hxx, hxy, hyy, hxxx, hxxy, hxyy, hyyy
  real *8 :: over4pi, over8pi
  data over4pi/0.07957747154594767d0/
  data over8pi/0.039788735772973836d0/

  zk1 = zpars(1)
  zk2 = zpars(2)

  if (abs(zk1).le.1d-6) then
    zk1 = zpars(2)
    zk2 = 0
  endif

  if (abs(zk1).le.1d-6) then
    ! biharmonic
    s = dx*dx + dy*dy
    s2 = s*s
    rl = log(s)
    g = s*rl*over8pi/2
    gx = dx*(rl + 1)*over8pi
    gy = dy*(rl + 1)*over8pi
    gxx = (rl + 1 + 2*dx*dx/s)*over8pi
    gxy = dx*dy/s*over4pi
    gyy = (rl + 1 + 2*dy*dy/s)*over8pi
    gxxx = dx*(dx*dx + 3*dy*dy)/s2*over4pi
    gxxy = dy*(dy*dy - dx*dx)/s2*over4pi
    gxyy = dx*(dx*dx - dy*dy)/s2*over4pi
    gyyy = dy*(dy*dy + 3*dx*dx)/s2*over4pi
    return
  endif

  call helmdiffgreen(zk1,dx,dy,g,gx,gy,gxx,gxy,gyy,gxxx,gxxy,gxyy,gyyy)

  if (abs(zk2).gt.1d-6) then
    call helmdiffgreen(zk2,dx,dy,h,hx,hy,hxx,hxy,hyy,hxxx,hxxy,hxyy,hyyy)
    g = g - h
    gx = gx - hx
    gy = gy - hy
    gxx = gxx - hxx
    gxy = gxy - hxy
    gyy = gyy - hyy
    gxxx = gxxx - hxxx
    gxxy = gxxy - hxxy
    gxyy = gxyy - hxyy
    gyyy = gyyy - hyyy
  else
    zk2 = 0
  endif

  zfac = 1/(zk1*zk1 - zk2*zk2)
  g = g*zfac
  gx = gx*zfac
  gy = gy*zfac
  gxx = gxx*zfac
  gxy = gxy*zfac
  gyy = gyy*zfac
  gxxx = gxxx*zfac
  gxxy = gxxy*zfac
  gxyy = gxyy*zfac
  gyyy = gyyy*zfac

end subroutine flex2d_gders
!
!
!
!
!
subroutine flex2d_g(src,ndt,targ,ndd,dpars,ndz,zpars,ndi,ipars,val)
  implicit none
  integer *8 ndt, ndd, ndz, ndi
  real *8 :: src(*), targ(ndt), dpars(ndd)
  integer *8 ipars(ndi)
  real *8 :: dx, dy
  complex *16 :: val
  complex *16 :: zpars(2)
  complex *16 :: gx, gy, gxx, gxy, gyy, gxxx, gxxy, gxyy, gyyy

  dx = targ(1)-src(1)
  dy = targ(2)-src(2)

  call flex2d_gders(zpars,dx,dy,val,gx,gy,gxx,gxy,gyy,gxxx,gxxy,gxyy,gyyy)

end subroutine flex2d_g
!
!
!
!
!
subroutine flex2d_gdn(src,ndt,targ,ndd,dpars,ndz,zpars,ndi,ipars,val)
  implicit none
  integer *8 ndt, ndd, ndz, ndi
  real *8 :: src(*), targ(ndt), dpars(ndd)
  integer *8 ipars(ndi)
  real *8 :: dx, dy
  real *8 :: nx, ny
  complex *16 :: val
  complex *16 :: zpars(2)
  complex *16 :: g, gx, gy, gxx, gxy, gyy, gxxx, gxxy, gxyy, gyyy

  dx = targ(1)-src(1)
  dy = targ(2)-src(2)

  nx = targ(10)
  ny = targ(11)

  call flex2d_gders(zpars,dx,dy,g,gx,gy,gxx,gxy,gyy,gxxx,gxxy,gxyy,gyyy)

  val = nx*gx + ny*gy

end subroutine flex2d_gdn
!
!
!
!
!
subroutine flex2d_gsupp2(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk, &
   ndi,ipars,val)
  implicit none
  integer *8 ndt, ndd, ndz, ndi
  real *8 :: srcinfo(*),targinfo(ndt),dpars(1)
  integer *8 ipars(ndi)
  real *8 :: dx, dy
  real *8 :: nu, nx, ny
  complex *16 :: zk(2)
  complex *16 :: val, g, gx, gy, gsxx, gsxy, gsyy
  complex *16 :: gsxxx, gsxxy, gsxyy, gsyyy
  !
  ! returns the second supported plate condition of the
  ! flexural volumetric kernel
  !

  dx = targinfo(1) - srcinfo(1)
  dy = targinfo(2) - srcinfo(2)

  call flex2d_gders(zk,dx,dy,g,gx,gy,gsxx,gsxy,gsyy, &
    gsxxx,gsxxy,gsxyy,gsyyy)

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
subroutine flex2d_gfree2(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk, &
   ndi,ipars,val)
  implicit none
  integer *8 ndt, ndd, ndz, ndi
  real *8 :: srcinfo(*),targinfo(ndt),dpars(1)
  integer *8 ipars(ndi)
  real *8 :: dx, dy, ds
  real *8 :: nu, nx, ny, kappa
  real *8 :: taux, tauy
  complex *16 :: zk(2)
  complex *16 :: val, g, gx, gy, gsxx, gsxy, gsyy
  complex *16 :: gsxxx, gsxxy, gsxyy, gsyyy
  !
  ! returns the second free plate condition of the
  ! flexural volumetric kernel
  !

  dx = targinfo(1) - srcinfo(1)
  dy = targinfo(2) - srcinfo(2)

  nx = targinfo(10)
  ny = targinfo(11)
  kappa = targinfo(13)

  taux = targinfo(4)
  tauy = targinfo(5)

  ds = sqrt(taux*taux + tauy*tauy)
  taux = taux/ds
  tauy = tauy/ds

  nu = dpars(1)

  call flex2d_gders(zk,dx,dy,g,gx,gy,gsxx,gsxy,gsyy, &
    gsxxx,gsxxy,gsxyy,gsyyy)

  val = gsxxx*nx*nx*nx + 3*gsxxy*nx*nx*ny + 3*gsxyy*nx*ny*ny + &
      gsyyy*ny*ny*ny + (2-nu)*(gsxxx*nx*taux*taux + &
      gsxxy*(taux*taux*ny + 2*taux*tauy*nx) + &
      gsxyy*(2*taux*tauy*ny + tauy*tauy*nx) + &
      gsyyy*tauy*tauy*ny)+kappa*(1-nu)*(gsxx*taux*taux + &
      2*gsxy*taux*tauy+gsyy*tauy*tauy - gsxx*nx*nx - &
      2*gsxy*nx*ny - gsyy*ny*ny)

  return
end subroutine flex2d_gfree2
!
!
!
!
!
subroutine flex2d_gvar(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk, &
   ndi,ipars,val)
  implicit none
  integer *8 ndt, ndd, ndz, ndi
  real *8 :: srcinfo(*),targinfo(ndt),dpars(ndd)
  integer *8 ipars(ndi)
  real *8 :: dx, dy
  complex *16 :: zk(2), c(7)
  complex *16 :: val, g, gx, gy, gxx, gxy, gyy
  complex *16 :: gxxx, gxxy, gxyy, gyyy
  integer *8 j
  !
  ! returns the kernel of the variable coefficient plate operator,
  !
  !   c1 d_x Lap G + c2 d_y Lap G + c3 Lap G + c4 G_yy + c5 G_xx
  !      + c6 G_xy + c7 G
  !
  ! with G the flexural Green's function of flex2d_gders, derivatives
  ! taken in the target, and the target coefficients
  ! c_j = targinfo(13+j) + i targinfo(20+j), j = 1,...,7
  !

  dx = targinfo(1) - srcinfo(1)
  dy = targinfo(2) - srcinfo(2)

  do j = 1,7
    c(j) = dcmplx(targinfo(13+j), targinfo(20+j))
  enddo

  call flex2d_gders(zk,dx,dy,g,gx,gy,gxx,gxy,gyy,gxxx,gxxy,gxyy,gyyy)

  val = c(1)*(gxxx + gxyy) + c(2)*(gxxy + gyyy) + c(3)*(gxx + gyy) + &
      c(4)*gyy + c(5)*gxx + c(6)*gxy + c(7)*g

  return
end subroutine flex2d_gvar
