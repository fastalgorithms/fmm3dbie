

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
!     *info(13) = meancrvs if it exists (THIS might need to be fixed later!)
!
!



subroutine lap_bel_log(src, ndt,targ, ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: src(*), targ(ndt),dpars(ndd)
  integer ipars(ndi)
  real *8 over4pi

  complex *16 :: zk
  real *8 :: val

  complex *16 :: ima

  data ima/(0.0d0,1.0d0)/
  data over4pi/0.07957747154594767d0/
  !
  ! returns the laplace beltrami operator applied to log, i.e.
  !
  ! K(x,y) = (n(x).(x-y)/|x-y|^2)^2 -2H(x)*(n(x).(x-y)/|x-y|^2)
  !

  dx=targ(1)-src(1)
  dy=targ(2)-src(2)
  dz=targ(3)-src(3)

  

  r2=dx**2+dy**2+dz**2
  rn = dx*targ(10) + dy*targ(11) + dz*targ(12)

  rnoverr2 = rn/r2

  val = rnoverr2*(rnoverr2 - 2*targ(13))*over4pi

  return
end subroutine lap_bel_log



