! Kernels for the Laplace-Beltrami and Helmholtz-Beltrami operators
!
! the following routines rely on the srcinfo and targinfo arrays
! containing the following information, standardized in the following
! way:
!
!     *info(1:3) = xyz
!     *info(4:6) = tangent vector 1
!     *info(7:9) = tangent vector 2
!     *info(10:12) = normal vector
!     *info(13) = mean curvature H
!
!


subroutine lap_bel_res(src, ndt,targ, ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: src(*), targ(ndt),dpars(ndd)
  integer *8 ipars(ndi)
  real *8 over4pi

  complex *16 :: zk
  real *8 :: val

  complex *16 :: ima

  data ima/(0.0d0,1.0d0)/
  data over4pi/0.07957747154594767d0/
  !
  ! returns the Laplace-Beltrami residual kernel applied to log|x-y|
  !
  !   K(x,y) = (n(x).(x-y)/|x-y|^2)^2 - 2H(x)*(n(x).(x-y)/|x-y|^2)
  !
  ! scaled by 1/4pi, where H(x) = targ(13) is the mean curvature at x
  !

  dx=targ(1)-src(1)
  dy=targ(2)-src(2)
  dz=targ(3)-src(3)

  r2=dx**2+dy**2+dz**2
  rn = dx*targ(10) + dy*targ(11) + dz*targ(12)

  rnoverr2 = rn/r2

  val = 4*rnoverr2*(rnoverr2 - targ(13))*over4pi

  return
end subroutine lap_bel_res






subroutine lap_bel_log(src, ndt,targ, ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: src(*), targ(ndt),dpars(ndd)
  integer *8 ipars(ndi)
  real *8 over4pi

  complex *16 :: zk
  real *8 :: val

  complex *16 :: ima

  data ima/(0.0d0,1.0d0)/
  data over4pi/0.07957747154594767d0/
  !
  ! returns the logarithmic kernel log|x-y| / (4pi)
  !   (single layer potential kernel for the Laplace-Beltrami operator)
  !

  dx=targ(1)-src(1)
  dy=targ(2)-src(2)
  dz=targ(3)-src(3)

  r2=dx**2+dy**2+dz**2
  val = log(r2)*over4pi

  return
end subroutine lap_bel_log

subroutine helm_bel_res(src, ndt,targ, ndd,dpars,ndz,zk,ndi,ipars,zval)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: src(*), targ(ndt),dpars(ndd)
  integer *8 ipars(ndi)
  real *8 over4pi

  complex *16 :: zk
  complex *16 :: zval, h0, h1, h2, tmp

  complex *16 :: ima

  data ima/(0.0d0,1.0d0)/
  data over4pi/0.07957747154594767d0/
  !
  ! returns the Helmholtz-Beltrami residual kernel applied to H_0(zk|x-y|)
  !
  !   K(x,y) = (-zk^2*H_2(zk*r)*(n(x).(x-y)/r)^2
  !              + 2H(x)*zk*H_1(zk*r)*(n(x).(x-y)/r)) / (4i)
  !
  ! where H(x) = targ(13) is the mean curvature at x
  !

  dx=targ(1)-src(1)
  dy=targ(2)-src(2)
  dz=targ(3)-src(3)

  r2=dx**2+dy**2+dz**2
  r=sqrt(r2)

  call hank101(zk*r,h0,h1)
  h2 = 2*h1/zk/r - h0

  rn = dx*targ(10) + dy*targ(11) + dz*targ(12)

  rnoverr = rn/r
  tmp = -zk*zk*h2*rnoverr*rnoverr
  zval = (tmp + zk*h1*2*targ(13)*rnoverr)/4/ima

  return
end subroutine helm_bel_res

subroutine helm_bel_hank(src, ndt,targ, ndd,dpars,ndz,zk,ndi,ipars,zval)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: src(*), targ(ndt),dpars(ndd)
  integer *8 ipars(ndi)
  real *8 over4pi

  complex *16 :: zk
  complex *16 :: zval, h0, h1

  complex *16 :: ima

  data ima/(0.0d0,1.0d0)/
  !
  ! returns the Hankel function kernel H_0^{(1)}(zk*|x-y|) / (4i)
  !   (single layer potential kernel for the Helmholtz-Beltrami operator)
  !

  dx=targ(1)-src(1)
  dy=targ(2)-src(2)
  dz=targ(3)-src(3)

  r=sqrt(dx**2+dy**2+dz**2)

  call hank101(zk*r,h0,h1)

  zval = h0/4/ima

  return
end subroutine helm_bel_hank


