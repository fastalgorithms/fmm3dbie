! Kernels for solving Laplace's equation in 2D
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



subroutine l2d_g(src, ndt,targ, ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: src(*), targ(ndt),dpars(ndd)
  integer ipars(ndi)
  real *8 :: dx, dy, r2
  real *8 :: over4pi
  real *8 :: val
  data over4pi/0.07957747154594767d0/
  !
  ! returns the Laplace volumetric potential
  !

  dx=targ(1)-src(1)
  dy=targ(2)-src(2)

  r2=dx**2+dy**2

  val = -over4pi*log(r2)

  return
end subroutine l2d_g

subroutine l2d_gdn(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk, &
   ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(*),targinfo(ndt),dpars(ndd)
  integer ipars(ndi)
  real *8 :: dx, dy, r2, rdotn
  real *8 :: val
  real *8 :: over4pi
  data over4pi/0.07957747154594767d0/
  !
  ! returns derivatives of the Laplace volumetric potential
  !     evaluated at the target locations and dotted with 
  !     the target normals defined by targinfo(10:11)
  !

  dx = targinfo(1) - srcinfo(1)
  dy = targinfo(2) - srcinfo(2)

  rdotn = dx*targinfo(10) + dy*targinfo(11)
  r2 = dx**2 + dy**2

  val =  -2*rdotn/r2*over4pi

  return
end subroutine l2d_gdn



subroutine l2d_slp(src, ndt,targ, ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: src(*), targ(ndt),dpars(ndd)
  integer ipars(ndi)
  real *8 :: dx, dy, r2
  real *8 :: over4pi
  real *8 :: val
  data over4pi/0.07957747154594767d0/
  !
  ! returns the Laplace single layer potential kernel
  !

  dx=targ(1)-src(1)
  dy=targ(2)-src(2)

  r2=dx**2+dy**2

  val = -over4pi*log(r2)

  return
end subroutine l2d_slp

subroutine l2d_sprime(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk, &
   ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(*),targinfo(ndt),dpars(ndd)
  integer ipars(ndi)
  real *8 :: dx, dy, r2, rdotn
  real *8 :: val
  real *8 :: over4pi
  data over4pi/0.07957747154594767d0/
  !
  ! returns the normal derivative of the single layer kernel
  !

  dx = targinfo(1) - srcinfo(1)
  dy = targinfo(2) - srcinfo(2)

  rdotn = dx*targinfo(10) + dy*targinfo(11)
  r2 = dx**2 + dy**2

  val =  -2*rdotn/r2*over4pi

  return
end subroutine l2d_sprime
