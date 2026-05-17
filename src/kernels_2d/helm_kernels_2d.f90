subroutine h2d_slp(src, ndt,targ, ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: src(*), targ(ndt),dpars(ndd)
  integer *8 ipars(ndi)
  real *8 :: dx, dy, r2
  real *8 :: over4pi
  complex *16 :: zk,ima,z,val,h0,h1
  data over4pi/0.07957747154594767d0/
  data ima/(0.0d0,1.0d0)/
  !
  ! returns the 2d helmholtz single layer potential kernel
  !

  ifexpon = 1
  
  dx=targ(1)-src(1)
  dy=targ(2)-src(2)
  r2=dx**2+dy**2
  z = zk*sqrt(r2)

  call hank103(z,h0,h1,ifexpon)

  val = 0.25d0*ima*h0

  return
end subroutine h2d_slp


subroutine h2d_sprime(src, ndt,targ, ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: src(*), targ(ndt),dpars(ndd)
  integer *8 ipars(ndi)
  real *8 :: dx, dy
  real *8 :: over4pi
  complex *16 :: zk,ima,z,val,h0,h1
  data over4pi/0.07957747154594767d0/
  data ima/(0.0d0,1.0d0)/
  !
  ! returns the 2d helmholtz single layer potential kernel
  !

  ifexpon = 1
  
  dx=targ(1)-src(1)
  dy=targ(2)-src(2)
  r=sqrt(dx**2+dy**2)
  z = zk*r

  rdotn = dx*targ(10) + dy*targ(11)

  call hank103(z,h0,h1,ifexpon)

  val = -0.25d0*ima*rdotn/r*zk*h1

  return
end subroutine h2d_sprime

