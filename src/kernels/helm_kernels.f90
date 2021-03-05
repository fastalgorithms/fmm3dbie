


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
  real *8 over4pi
  integer ipars(ndi)
  complex *16 :: zk, val
  complex *16 :: ima

  data ima/(0.0d0,1.0d0)/
  data over4pi/0.07957747154594767d0/
  !
  ! returns the helmholtz single layer potential kernel
  !

  dx=targ(1)-src(1)
  dy=targ(2)-src(2)
  dz=targ(3)-src(3)

  r=sqrt(dx**2+dy**2+dz**2)

  val =  exp(ima*zk*r)/r*over4pi

  return
end subroutine h3d_slp





subroutine h3d_dlp(srcinfo,ndt,targ,ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(12), targ(ndt),dpars(ndd)
  integer ipars(ndi)
  complex *16 :: zk, val

  real *8 :: src(3), srcnorm(3)
  real *8 over4pi
  complex *16 :: ima
  data ima/(0.0d0,1.0d0)/
  data over4pi/0.07957747154594767d0/
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

  val =  d*(1.0d0 - ima*zk*r)*exp(ima*zk*r)/(r**3)*over4pi

  return
end subroutine h3d_dlp






subroutine h3d_sprime(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(*), targinfo(12),dpars(ndd)
  integer ipars(ndi)
  complex *16 :: zk, val

  complex *16 :: ima
  real *8 over4pi
  data ima/(0.0d0,1.0d0)/
  data over4pi/0.07957747154594767d0/
  !
  ! returns the normal derivative of the single layer kernel
  !

  dx=targinfo(1)-srcinfo(1)
  dy=targinfo(2)-srcinfo(2)
  dz=targinfo(3)-srcinfo(3)

  d = dx*targinfo(10) + dy*targinfo(11) + dz*targinfo(12)
  r=sqrt(dx**2+dy**2+dz**2)

  val =  -d*(1.0d0 - ima*zk*r)*exp(ima*zk*r)/(r**3)*over4pi

  return
end subroutine h3d_sprime




subroutine h3d_comb(srcinfo, ndt,targ, ndd,dpars,ndz,zpars,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(12), targ(ndt),dpars(ndd)
  integer ipars(ndi)
  complex *16 :: zk, val, zpars(ndz), alpha,beta

  real *8 :: src(3), srcnorm(3)
  real *8 over4pi
  complex *16 :: ima
  data ima/(0.0d0,1.0d0)/
  data over4pi/0.07957747154594767d0/
  !
  ! returns the helmholtz combined field kernel 
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


  val =  (beta*d*(1.0d0 - ima*zk*r)/(r**3)+alpha/r)*exp(ima*zk*r)* &
    over4pi

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
  real *8 :: over4pi
  data ima/(0.0d0,1.0d0)/
  data over4pi/0.07957747154594767d0/
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
     1-ima*zk*r)*exp(ima*zk*r)/(r**3)*over4pi

  return
end subroutine h3d_qlp






subroutine h3d_dprime_diff(srcinfo, ndt,targ, ndd,dpars,ndz,zpars,ndi, &
  ipars,val)
!
!  This subroutine returns D_{k1}' - D_{k2}'
!  
!  zpars(1) = k1
!  zpars(2) = k2
!
!  This routine does so by computing D_{k1}'-D_{0}' and
!  D_{k2}'-D_{0}' in a stabilized manner and then return
!  the difference
!

  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(12), targ(ndt),dpars(ndd)
  integer ipars(ndi)
  complex *16 :: zk, val, zpars(ndz), alpha,beta

  real *8 :: src(3), rns(3),rnt(3),rnsdot,rntdot,rnstdot 
  complex *16 :: ima,z1,z2,ztmp,zexp,sexp,ztmp2
  real *8 :: over4pi
  data ima/(0.0d0,1.0d0)/
  data over4pi/0.07957747154594767d0/

  done = 1
  pi = atan(done)*4



  src(1)=srcinfo(1)
  src(2)=srcinfo(2)
  src(3)=srcinfo(3)

  rns(1)=srcinfo(10)
  rns(2)=srcinfo(11)
  rns(3)=srcinfo(12)

  rnt(1) = targ(10)
  rnt(2) = targ(11)
  rnt(3) = targ(12)

  dx=targ(1)-src(1)
  dy=targ(2)-src(2)
  dz=targ(3)-src(3)

  rnsdot = dx*rns(1) + dy*rns(2) + dz*rns(3)
  rntdot = dx*rnt(1) + dy*rnt(2) + dz*rnt(3)
  rnstdot = rns(1)*rnt(1) + rns(2)*rnt(2) + rns(3)*rnt(3)

  r=sqrt(dx**2+dy**2+dz**2)
  rinv = 1.0d0/r
  rinv3 = rinv**3
  rinv5 = rinv**5

  ztmp = ima*zpars(1)*r
  ztmp2 = ztmp**2
  zexp = exp(ztmp/2)
  sexp = sin(zpars(1)*r/2)*2*ima*zexp

  z1 = rinv3*((ztmp-1)*sexp + ztmp)
  z2 = rinv5*((ztmp2-3*ztmp+3)*sexp+ztmp2-3*ztmp)

  val = -rntdot*rnsdot*z2 - rnstdot*z1



  ztmp = ima*zpars(2)*r
  ztmp2 = ztmp**2
  zexp = exp(ztmp/2)
  sexp = sin(zpars(2)*r/2)*2*ima*zexp

  z1 = rinv3*((ztmp-1)*sexp + ztmp)
  z2 = rinv5*((ztmp2-3*ztmp+3)*sexp+ztmp2-3*ztmp)
  val = val + rntdot*rnsdot*z2 + rnstdot*z1
  val = val*over4pi



  return
end subroutine h3d_dprime_diff




subroutine h3d_slp_diff(src, ndt,targ, ndd,dpars,ndz,zpars,ndi,ipars,val)
!
!  This subroutine evaluates the kernel
!    a0 S_{k0} - a1 S_{k1}
!
!  where zpars(1) = k0,
!        zpars(2) = k1,
!        zpars(3) = a0
!        zpars(4) = a1

  implicit real *8 (a-h,o-z)
  real *8 :: src(*), targ(ndt),dpars(ndd)
  real *8 over4pi
  integer ipars(ndi)
  complex *16 :: zpars(ndz), val
  complex *16 :: ima

  data ima/(0.0d0,1.0d0)/
  data over4pi/0.07957747154594767d0/

  dx=targ(1)-src(1)
  dy=targ(2)-src(2)
  dz=targ(3)-src(3)

  r=sqrt(dx**2+dy**2+dz**2)
  rinv = over4pi/r

  val=(zpars(3)*exp(ima*zpars(1)*r) -zpars(4)*exp(ima*zpars(2)*r))*rinv

  return
end subroutine h3d_slp_diff





subroutine h3d_dlp_diff(srcinfo,ndt,targ,ndd,dpars,ndz,zpars,ndi, &
  ipars,val)
!
!  This subroutine evaluates the kernel
!    a0 D_{k0} - a_{1} D_{k1}
!
!  where zpars(1) = k0,
!        zpars(2) = k1,
!        zpars(3) = a0
!        zpars(4) = a1

  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(12), targ(ndt),dpars(ndd)
  integer ipars(ndi)
  complex *16 :: zpars(ndz), val

  real *8 :: src(3), srcnorm(3)
  real *8 over4pi
  complex *16 :: ima,zk0,zk1
  data ima/(0.0d0,1.0d0)/
  data over4pi/0.07957747154594767d0/

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
  rtmp = over4pi*d/r**3

  zk0 = zpars(1)
  zk1 = zpars(2)

  val =  (zpars(3)*(1.0d0 - ima*zk0*r)*exp(ima*zk0*r) - &
       zpars(4)*(1.0d0 - ima*zk1*r)*exp(ima*zk1*r))*rtmp

  return
end subroutine h3d_dlp_diff






subroutine h3d_sprime_diff(srcinfo,ndt,targinfo,ndd,dpars,ndz,zpars,&
  ndi,ipars,val)
!
!  This subroutine evaluates the kernel
!    a0 S_{k0}' - a1 S_{k1}'
!
!  where zpars(1) = k0,
!        zpars(2) = k1,
!        zpars(3) = a0
!        zpars(4) = a1

  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(*), targinfo(12),dpars(ndd)
  integer ipars(ndi)
  complex *16 :: zpars(ndz), val,zk0,zk1

  complex *16 :: ima
  real *8 over4pi
  data ima/(0.0d0,1.0d0)/
  data over4pi/0.07957747154594767d0/

  dx=targinfo(1)-srcinfo(1)
  dy=targinfo(2)-srcinfo(2)
  dz=targinfo(3)-srcinfo(3)

  d = dx*targinfo(10) + dy*targinfo(11) + dz*targinfo(12)
  r=sqrt(dx**2+dy**2+dz**2)
  rtmp = d*over4pi/r**3

  zk0 = zpars(1)
  zk1 = zpars(2)

  val =  -(zpars(3)*(1.0d0 - ima*zk0*r)*exp(ima*zk0*r) - &
    zpars(4)*(1.0d0-ima*zk1*r)*exp(ima*zk1*r))*rtmp

  return
end subroutine h3d_sprime_diff









subroutine h3d_sgradx(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(*), targinfo(*),dpars(ndd)
  complex *16 val
  integer ipars(ndi)
  complex *16 :: zk,ima

  real *8 over4pi
  data over4pi/0.07957747154594767d0/
  data ima/(0.0d0,1.0d0)/
  !
  ! returns the normal derivative of the single layer kernel
  !

  dx=targinfo(1)-srcinfo(1)
  dy=targinfo(2)-srcinfo(2)
  dz=targinfo(3)-srcinfo(3)

  r=sqrt(dx**2+dy**2+dz**2)

  val =  -dx*(1.0d0 - ima*zk*r)*exp(ima*zk*r)/(r**3)*over4pi

  return
end subroutine h3d_sgradx








subroutine h3d_sgrady(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(*), targinfo(*),dpars(ndd)
  complex *16 :: val
  integer ipars(ndi)
  complex *16 :: zk,ima

  real *8 over4pi
  data over4pi/0.07957747154594767d0/
  data ima/(0.0d0,1.0d0)/
  !
  ! returns the normal derivative of the single layer kernel
  !

  dx=targinfo(1)-srcinfo(1)
  dy=targinfo(2)-srcinfo(2)
  dz=targinfo(3)-srcinfo(3)

  r=sqrt(dx**2+dy**2+dz**2)

  val =  -dy*(1.0d0 - ima*zk*r)*exp(ima*zk*r)/(r**3)*over4pi

  return
end subroutine h3d_sgrady











subroutine h3d_sgradz(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(*), targinfo(*),dpars(ndd)
  complex *16 :: val
  integer ipars(ndi)
  complex *16 :: zk,ima

  real *8 over4pi
  data over4pi/0.07957747154594767d0/
  data ima/(0.0d0,1.0d0)/
  !
  ! returns the normal derivative of the single layer kernel
  !

  dx=targinfo(1)-srcinfo(1)
  dy=targinfo(2)-srcinfo(2)
  dz=targinfo(3)-srcinfo(3)

  r=sqrt(dx**2+dy**2+dz**2)

  val =  -dz*(1.0d0 - ima*zk*r)*exp(ima*zk*r)/(r**3)*over4pi

  return
end subroutine h3d_sgradz




