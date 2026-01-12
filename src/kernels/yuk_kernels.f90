


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



subroutine y3d_slp(src, ndt,targ, ndd,dzk,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: src(*), targ(ndt),dzk
  real *8 over4pi,val
  integer *8 ipars(ndi)
  complex *16 :: zk
  data over4pi/0.07957747154594767d0/
  !
  ! returns the yukawa single layer potential kernel
  !

  dx=targ(1)-src(1)
  dy=targ(2)-src(2)
  dz=targ(3)-src(3)

  r=sqrt(dx**2+dy**2+dz**2)

  val =  exp(-dzk*r)/r*over4pi

  return
end subroutine y3d_slp





subroutine y3d_dlp(srcinfo,ndt,targ,ndd,dzk,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: srcinfo(12), targ(ndt)
  real *8 :: dzk,val
  integer *8 ipars(ndi)
  complex *16 :: zk

  real *8 :: src(3), srcnorm(3)
  real *8 over4pi
  data over4pi/0.07957747154594767d0/
  !
  ! returns the yukawa double layer kernel
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

  val =  d*(1.0d0 + dzk*r)*exp(-dzk*r)/(r**3)*over4pi

  return
end subroutine y3d_dlp






subroutine y3d_sprime(srcinfo,ndt,targinfo,ndd,dzk,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: srcinfo(*), targinfo(12),dzk,val
  integer *8 ipars(ndi)
  complex *16 :: zk

  real *8 over4pi
  data over4pi/0.07957747154594767d0/
  !
  ! returns the normal derivative of the single layer kernel
  !

  dx=targinfo(1)-srcinfo(1)
  dy=targinfo(2)-srcinfo(2)
  dz=targinfo(3)-srcinfo(3)

  d = dx*targinfo(10) + dy*targinfo(11) + dz*targinfo(12)
  r=sqrt(dx**2+dy**2+dz**2)

  val =  -d*(1.0d0 + dzk*r)*exp(-dzk*r)/(r**3)*over4pi

  return
end subroutine y3d_sprime




subroutine y3d_comb(srcinfo, ndt,targ, ndd,dpars,ndz,zpars,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: srcinfo(12), targ(ndt),dpars(ndd),dzk,val
  real *8 :: alpha,beta
  integer *8 ipars(ndi)
  complex *16 :: zk, zpars(ndz)

  real *8 :: src(3), srcnorm(3)
  real *8 over4pi
  data over4pi/0.07957747154594767d0/
  !
  ! returns the yukawa combined field kernel 
  !

  dzk = dpars(1)
  alpha = dpars(2)
  beta = dpars(3)
  
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


  val =  (beta*d*(1.0d0 + dzk*r)/(r**3)+alpha/r)*exp(-dzk*r)* &
    over4pi

  return
end subroutine y3d_comb





subroutine y3d_sgradx(srcinfo,ndt,targinfo,ndd,dzk,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: srcinfo(*), targinfo(*),dzk,val
  integer *8 ipars(ndi)
  complex *16 :: zk

  real *8 over4pi
  data over4pi/0.07957747154594767d0/
  !
  ! returns the normal derivative of the single layer kernel
  !

  dx=targinfo(1)-srcinfo(1)
  dy=targinfo(2)-srcinfo(2)
  dz=targinfo(3)-srcinfo(3)

  r=sqrt(dx**2+dy**2+dz**2)

  val =  -dx*(1.0d0 + dzk*r)*exp(-dzk*r)/(r**3)*over4pi

  return
end subroutine y3d_sgradx








subroutine y3d_sgrady(srcinfo,ndt,targinfo,ndd,dzk,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: srcinfo(*), targinfo(*),dzk,val
  integer *8 ipars(ndi)
  complex *16 :: zk

  real *8 over4pi
  data over4pi/0.07957747154594767d0/
  !
  ! returns the normal derivative of the single layer kernel
  !

  dx=targinfo(1)-srcinfo(1)
  dy=targinfo(2)-srcinfo(2)
  dz=targinfo(3)-srcinfo(3)

  r=sqrt(dx**2+dy**2+dz**2)

  val =  -dy*(1.0d0 + dzk*r)*exp(-dzk*r)/(r**3)*over4pi

  return
end subroutine y3d_sgrady











subroutine y3d_sgradz(srcinfo,ndt,targinfo,ndd,dzk,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: srcinfo(*), targinfo(*),dzk,val
  integer *8 ipars(ndi)
  complex *16 :: zk

  real *8 over4pi
  data over4pi/0.07957747154594767d0/
  !
  ! returns the normal derivative of the single layer kernel
  !

  dx=targinfo(1)-srcinfo(1)
  dy=targinfo(2)-srcinfo(2)
  dz=targinfo(3)-srcinfo(3)

  r=sqrt(dx**2+dy**2+dz**2)

  val =  -dz*(1.0d0 + dzk*r)*exp(-dzk*r)/(r**3)*over4pi

  return
end subroutine y3d_sgradz




