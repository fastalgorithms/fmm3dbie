


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



subroutine l3d_slp(src, ndt,targ, ndd,dpars,ndz,zk,ndi,ipars,val)
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
  ! returns the Laplace single layer potential kernel
  !

  dx=targ(1)-src(1)
  dy=targ(2)-src(2)
  dz=targ(3)-src(3)

  r=sqrt(dx**2+dy**2+dz**2)

  val =  over4pi/r

  return
end subroutine l3d_slp





subroutine l3d_dlp(srcinfo,ndt,targ,ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: srcinfo(12), targ(ndt),dpars(ndd)
  integer *8 ipars(ndi)
  complex *16 :: zk
  real *8 :: val
  real *8 :: over4pi

  real *8 :: src(3), srcnorm(3)
  complex *16 :: ima
  data ima/(0.0d0,1.0d0)/
  data over4pi/0.07957747154594767d0/
  !
  ! returns the Laplace double layer kernel
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
  r = sqrt(dx**2+dy**2+dz**2)

  val =  d/(r**3)*over4pi

  return
end subroutine l3d_dlp






subroutine l3d_sprime(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: srcinfo(*), targinfo(12),dpars(ndd)
  integer *8 ipars(ndi)
  real *8 :: val
  real *8 :: over4pi
  complex *16 :: zk

  data over4pi/0.07957747154594767d0/
  !
  ! returns the normal derivative of the single layer kernel
  !

  dx = targinfo(1) - srcinfo(1)
  dy = targinfo(2) - srcinfo(2)
  dz = targinfo(3) - srcinfo(3)

  d = dx*targinfo(10) + dy*targinfo(11) + dz*targinfo(12)
  r = sqrt(dx**2 + dy**2 + dz**2)

  val =  -d/(r**3)*over4pi

  return
end subroutine l3d_sprime




subroutine l3d_comb(srcinfo, ndt,targ, ndd,dpars,ndz,zpars,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: srcinfo(12), targ(ndt),dpars(ndd)
  integer *8 ipars(ndi)
  complex *16 :: zk, zpars(ndz)
  real *8 :: alpha,beta
  real *8 :: val

  real *8 :: src(3), srcnorm(3)
  real *8 :: over4pi
  complex *16 :: ima
  data ima/(0.0d0,1.0d0)/
  data over4pi/0.07957747154594767d0/
  !
  ! returns the Laplace double layer kernel
  !

  alpha = dpars(1)
  beta = dpars(2)
  


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


  val =  (beta*d/(r**3)+alpha/r)*over4pi

  return
end subroutine l3d_comb
!
!
!
!
!
!

subroutine l3d_qlp(srcinfo, ndt,targ,ndd, dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: srcinfo(12), targ(ndt),dpars(ndd)
  integer *8 ipars(ndi)
  complex *16 :: zk
  real *8 val

  real *8 :: src(3), srcnorm(3)
  real *8 dotprod
  real *8 over4pi
  complex *16 :: ima
  data ima/(0.0d0,1.0d0)/
  data over4pi/0.07957747154594767d0/
  !
  ! returns the Laplace quadruple layer kernel
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

  val =  (d**2*(-3/rsq) + 1)/(r**3)*over4pi

  return
end subroutine l3d_qlp










subroutine l3d_sgradx(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: srcinfo(*), targinfo(ndt),dpars(ndd)
  integer *8 ipars(ndi)
  real *8 :: val
  real *8 :: over4pi
  complex *16 :: zk

  complex *16 :: ima
  data ima/(0.0d0,1.0d0)/
  data over4pi/0.07957747154594767d0/
  !
  ! returns the normal derivative of the single layer kernel
  !

  dx=targinfo(1)-srcinfo(1)
  dy=targinfo(2)-srcinfo(2)
  dz=targinfo(3)-srcinfo(3)

  r=sqrt(dx**2+dy**2+dz**2)

  val =  -dx/(r**3)*over4pi

  return
end subroutine l3d_sgradx








subroutine l3d_sgrady(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: srcinfo(*), targinfo(ndt),dpars(ndd)
  integer *8 ipars(ndi)
  real *8 :: val
  real *8 :: over4pi
  complex *16 :: zk

  complex *16 :: ima
  data ima/(0.0d0,1.0d0)/
  data over4pi/0.07957747154594767d0/
  !
  ! returns the normal derivative of the single layer kernel
  !

  dx=targinfo(1)-srcinfo(1)
  dy=targinfo(2)-srcinfo(2)
  dz=targinfo(3)-srcinfo(3)

  r=sqrt(dx**2+dy**2+dz**2)

  val =  -dy/(r**3)*over4pi

  return
end subroutine l3d_sgrady












subroutine l3d_sgradz(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: srcinfo(*), targinfo(ndt),dpars(ndd)
  integer *8 ipars(ndi)
  real *8 :: val
  real *8 :: over4pi
  complex *16 :: zk

  data over4pi/0.07957747154594767d0/
  !
  ! returns the normal derivative of the single layer kernel
  !

  dx=targinfo(1)-srcinfo(1)
  dy=targinfo(2)-srcinfo(2)
  dz=targinfo(3)-srcinfo(3)

  r=sqrt(dx**2+dy**2+dz**2)

  val =  -dz/(r**3)*over4pi

  return
end subroutine l3d_sgradz






subroutine l3d_spp_sum_dp(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk, &
   ndi,ipars,val) 

!
!---------------------
!  This subroutine evaluates the difference kernel
!   (S'' + D')*4*pi (to be consistent with fmm scaling) 
!
!  Input arguments:
!  
!    - srcinfo: real *8 (12)
!        Source information
!    - ndt: integer *8
!        must be at least 12, with first twelve components 
!        being the standard targ info, xyz,dxyz/du,dxyz/dv,normal
!    - targinfo: real *8 (12)
!        target information
!    - ndd: integer *8
!        dpars not used
!    - dpars: real *8
!        dpars not used
!    - ndz: integer *8
!        zpars not used
!    - zpars: double complex 
!        zpars not used
!    - ndi: integer *8
!        ipars not used
!    - ipars: integer *8
!        ipars not used
!
!  Output arugments:
!    - val: double precision
!        s'' + d'
!    
!------------------

  implicit none
  integer *8 ndd,ndi,ndz,ndt
  real *8 :: srcinfo(*), targinfo(ndt),dpars(ndd)
  integer *8 ipars(ndi)
  real *8, intent(out) :: val
  complex *16 :: zk

  real *8 dx(3),dns(3),dnt(3),r,d,drns,drnt,r3inv,r5inv
  real *8 dnsnt,over4pi,rinv
  data over4pi/0.07957747154594767d0/

  dx(1) = targinfo(1) - srcinfo(1)
  dx(2) = targinfo(2) - srcinfo(2)
  dx(3) = targinfo(3) - srcinfo(3)

  dns(1) = srcinfo(10)
  dns(2) = srcinfo(11)
  dns(3) = srcinfo(12)

  dnt(1) = targinfo(10)
  dnt(2) = targinfo(11)
  dnt(3) = targinfo(12)

  r = sqrt(dx(1)**2 + dx(2)**2 + dx(3)**2)
  rinv = 1.0d0/r
  r3inv = rinv**3
  r5inv = rinv**5
  drns = dx(1)*dns(1) + dx(2)*dns(2) + dx(3)*dns(3)
  drnt = dx(1)*dnt(1) + dx(2)*dnt(2) + dx(3)*dnt(3)
  dnsnt = dnt(1)*dns(1) + dnt(2)*dns(2) + dnt(3)*dns(3)
  
!
!  second derivative of single layer
!
  val = 3*drnt**2*r5inv - r3inv

!
! add in derivative of double layer
!
  val = val - 3*drns*drnt*r5inv + dnsnt*r3inv
  val = val*over4pi

end subroutine l3d_spp_sum_dp
!
!
!
!
!
subroutine l3d_sintprime(srcinfo, ndt, targinfo, ndd, dpars, ndz, & 
   zk, ndi, ipars, val)
!
!  
!  This subroutine returns the Neumann data corresponding
!  to using S + (kap1 + kap2)/2 \phi_a + (kap1 - kap2) \psi_a
!
!  where \phi is obtained by the indefinite integral of the
!  Green's function, and \psi is the harmonic function
!  obtained by integrating the greens function twice
!  and taking two derivatives in the tangent plane
!
!  The functions are defined as follows:
!    \phi = log(\sqrt(x^2 + y^2 + z^2) + z) 
!    \psi = (x^2 - y^2)/(\sqrt(x^2 + y^2 + z^2) + z)^2
!
!    \phi_t = \phi(x,y,z) - \sum_{j=0}^{2} c_{j} \phi(x,y,z-z_{j})
!    \psi_t = \psi(x,y,z) - \sum_{j=0}^{2} c_{j} \psi(x,y,z-z_{j})
!
!  where c_{j} = -\prod_{\ell \neq j} z_{\ell}/(\prod_{\ell\neq j} (z_{j} - z_{\ell}))
!  and z_{j}'s are provided through the target info
!
!  Finally, 
!  \phi_a(x,y,z) = \phi_t(u,v,q)
!  \psi_a(x,y,z) = \psi_t(u,v,q) where
!
!  u,v,q are the coordinates relative with the origin at the target
!  with the frame formed by the tangent vectors and the normal
!  at the target.
!
!  In particular, if u_{r}, v_{r}, are orthogonal unit vectors in the 
!  tangent plane (obtained by srcvals(4:6,:), 
!  srcvals(4:6,:) \times srcvals(10:12,:))
!
!  Q: Are the principal curvatures kap1 and kap2 in the source
!  variable or the target variable? From a representation perspective
!  they have to be in the source variable because the definition 
!  does not extend to the target variable right? 
!  
!
!

  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: srcinfo(*), targinfo(12), dpars(ndd)
  integer *8 :: ipars(ndi)
  real *8 :: val
  real *8 :: over4pi
  complex *16 :: zk

  data over4pi/0.07957747154594767d0/
!
!
! TODO: FIX THIS ROUTINE
!

  !
  ! returns the normal derivative of the single layer kernel
  !

  dx = targinfo(1) - srcinfo(1)
  dy = targinfo(2) - srcinfo(2)
  dz = targinfo(3) - srcinfo(3)

  d = dx*targinfo(10) + dy*targinfo(11) + dz*targinfo(12)
  r = sqrt(dx**2 + dy**2 + dz**2)

  val =  -d/(r**3)*over4pi

  return
end subroutine l3d_sintprime



