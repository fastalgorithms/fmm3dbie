


!
! the following routines rely on the srcinfo and targinfo arrays
! containing the following information, standardized in the following
! way:
!
!     *info(1:3) = xyz
!     *info(4:6) = tanget vector 1
!     *info(7:9) = tangent vector 2
!     *info(10:12) = normal vector
!
!



subroutine h3d_slp(src, ndt,targ, ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: src(*), targ(ndt),dpars(ndd)
  integer ipars(ndi)

  complex *16 :: zk, val

  complex *16 :: ima

  data ima/(0.0d0,1.0d0)/
  !
  ! returns the helmholtz single layer potential kernel
  !

  dx=targ(1)-src(1)
  dy=targ(2)-src(2)
  dz=targ(3)-src(3)

  r=sqrt(dx**2+dy**2+dz**2)

  val =  exp(ima*zk*r)/r

  return
end subroutine h3d_slp





subroutine h3d_dlp(srcinfo,ndt,targ,ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(12), targ(ndt),dpars(ndd)
  integer ipars(ndi)
  complex *16 :: zk, val

  real *8 :: src(3), srcnorm(3)
  complex *16 :: ima
  data ima/(0.0d0,1.0d0)/
  !
  ! returns the helmholtz double layer kernel
  !

  src(1)=srcinfo(1)
  src(2)=srcinfo(2)
  src(3)=srcinfo(3)
  srcnorm(1)=srcinfo(10)
  srcnorm(2)=srcinfo(11)
  srcnorm(3)=srcinfo(12)

!  call prin2('srcinfo=*',srcinfo,12)

  dx=targ(1)-src(1)
  dy=targ(2)-src(2)
  dz=targ(3)-src(3)

  d = dx*srcnorm(1) + dy*srcnorm(2) + dz*srcnorm(3)
  r=sqrt(dx**2+dy**2+dz**2)

  val =  d*(1.0d0 - ima*zk*r)*exp(ima*zk*r)/(r**3)

  return
end subroutine h3d_dlp






subroutine h3d_sprime(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(*), targinfo(12),dpars(ndd)
  integer ipars(ndi)
  complex *16 :: zk, val

  complex *16 :: ima
  data ima/(0.0d0,1.0d0)/
  !
  ! returns the normal derivative of the single layer kernel
  !

  dx=targinfo(1)-srcinfo(1)
  dy=targinfo(2)-srcinfo(2)
  dz=targinfo(3)-srcinfo(3)

  d = dx*targinfo(10) + dy*targinfo(11) + dz*targinfo(12)
  r=sqrt(dx**2+dy**2+dz**2)

  val =  -d*(1.0d0 - ima*zk*r)*exp(ima*zk*r)/(r**3)

  return
end subroutine h3d_sprime




subroutine h3d_comb(srcinfo, ndt,targ, ndd,dpars,ndz,zpars,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(12), targ(ndt),dpars(ndd)
  integer ipars(ndi)
  complex *16 :: zk, val, zpars(ndz), alpha,beta

  real *8 :: src(3), srcnorm(3)
  complex *16 :: ima
  data ima/(0.0d0,1.0d0)/
  !
  ! returns the helmholtz double layer kernel
  !

  zk = zpars(1)
  alpha = zpars(2)
  beta = zpars(3)
  


  src(1)=srcinfo(1)
  src(2)=srcinfo(2)
  src(3)=srcinfo(3)
  srcnorm(1)=srcinfo(10)
  srcnorm(2)=srcinfo(11)
  srcnorm(3)=srcinfo(12)

  dx=targ(1)-src(1)
  dy=targ(2)-src(2)
  dz=targ(3)-src(3)

  d = dx*srcnorm(1) + dy*srcnorm(2) + dz*srcnorm(3)
  r=sqrt(dx**2+dy**2+dz**2)


  val =  (beta*d*(1.0d0 - ima*zk*r)/(r**3)+alpha/r)*exp(ima*zk*r)

  return
end subroutine h3d_comb




subroutine h3d_qlp(srcinfo, ndt,targ,ndd, dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(12), targ(ndt),dpars(ndd)
  integer ipars(ndi)
  complex *16 :: zk, val

  real *8 :: src(3), srcnorm(3)
  real *8 dotprod
  complex *16 :: ima
  data ima/(0.0d0,1.0d0)/
  !
  ! returns the helmholtz quadruple layer kernel
  !

  src(1)=srcinfo(1)
  src(2)=srcinfo(2)
  src(3)=srcinfo(3)
  srcnorm(1)=srcinfo(10)
  srcnorm(2)=srcinfo(11)
  srcnorm(3)=srcinfo(12)

!  call prin2('srcinfo=*',srcinfo,12)

  dx=targ(1)-src(1)
  dy=targ(2)-src(2)
  dz=targ(3)-src(3)

  d = dx*srcnorm(1) + dy*srcnorm(2) + dz*srcnorm(3)
  rsq = dx**2 + dy**2 + dz**2
  r=sqrt(rsq)

  val =  (d**2*(-(ima*zk)**2 -3/rsq + 3*ima*zk/r) +  &
     1-ima*zk*r)*exp(ima*zk*r)/(r**3)

  return
end subroutine h3d_qlp




