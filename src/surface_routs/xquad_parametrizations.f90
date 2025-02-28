!
!
!


subroutine xquad_wtorus_eval(iquad, u, v, xyz, dxyzduv, quadinfo, &
    radii, scales, p4)
  implicit real *8 (a-h,o-z)
  real *8 :: xyz(3), dxyzduv(3,2), quadinfo(3,3,*), scales(3)
  real *8 :: radii(3), dxyzdst(3,2)

  !
  ! project the quad iquad in quadinfo onto a torus
  !
  !    Input:
  ! iquad - quad number to map
  ! u,v - local uv coordinates on quad iquad
  ! quadinfo - flat skeleton quad info
  ! radii - the two radii defining the torus, the third
  !     radius is the radius of osciallation
  ! scales - scaling for x,y,z components from the standard torus
  ! p4 - number of oscillations (must be an integer currently recast
  !   as a double precision number)
  !
  !    Output:
  ! xyz - point on the sphere
  ! dxyzduv - first and second derivative information
  !
  !

  x0=quadinfo(1,1,iquad)
  y0=quadinfo(2,1,iquad)
  z0=quadinfo(3,1,iquad)

  x1=quadinfo(1,2,iquad)
  y1=quadinfo(2,2,iquad)
  z1=quadinfo(3,2,iquad)

  x2=quadinfo(1,3,iquad)
  y2=quadinfo(2,3,iquad)
  z2=quadinfo(3,3,iquad)


  nosc = p4


  s = x0+(1.0d0+u)/2*(x1-x0)+(1.0d0+v)/2*(x2-x0)
  t = y0+(1.0d0+u)/2*(y1-y0)+(1.0d0+v)/2*(y2-y0)

  call wtorus_eval(s, t, radii, scales, nosc, xyz, dxyzdst)

  dsdu = (x1-x0)/2
  dsdv = (x2-x0)/2
  dtdu = (y1-y0)/2
  dtdv = (y2-y0)/2

  dxyzduv(1,1) = dxyzdst(1,1)*dsdu + dxyzdst(1,2)*dtdu
  dxyzduv(2,1) = dxyzdst(2,1)*dsdu + dxyzdst(2,2)*dtdu
  dxyzduv(3,1) = dxyzdst(3,1)*dsdu + dxyzdst(3,2)*dtdu

  dxyzduv(1,2) = dxyzdst(1,1)*dsdv + dxyzdst(1,2)*dtdv
  dxyzduv(2,2) = dxyzdst(2,1)*dsdv + dxyzdst(2,2)*dtdv
  dxyzduv(3,2) = dxyzdst(3,1)*dsdv + dxyzdst(3,2)*dtdv

  return
end subroutine xquad_wtorus_eval
!
!
!
!
!
!
!


subroutine xquad_stell_eval(iquad, u, v, xyz, dxyzduv, quadinfo, &
    deltas, m, n)
  implicit real *8 (a-h,o-z)
  real *8 :: xyz(3), dxyzduv(3,2), quadinfo(3,3,*), deltas(-1:m,-1:n)
  real *8 :: dxyzds(3),dxyzdt(3)

  !
  ! project the quad iquad in quadinfo onto a stellarator
  !
  !    Input:
  ! iquad - quad number to map
  ! u,v - local uv coordinates on quad iquad
  ! quadinfo - flat skeleton quad info
  !
  !    Output:
  ! xyz - point on the sphere
  ! dxyzduv - first and second derivative information
  !
  !


  x0=quadinfo(1,1,iquad)
  y0=quadinfo(2,1,iquad)
  z0=quadinfo(3,1,iquad)

  x1=quadinfo(1,2,iquad)
  y1=quadinfo(2,2,iquad)
  z1=quadinfo(3,2,iquad)

  x2=quadinfo(1,3,iquad)
  y2=quadinfo(2,3,iquad)
  z2=quadinfo(3,3,iquad)


  !
  ! ... process the geometry, return the point location on the almond
  ! and the derivatives with respect to u and v
  !
  s = x0+(1.0d0+u)/2*(x1-x0)+(1.0d0+v)/2*(x2-x0)
  t = y0+(1.0d0+u)/2*(y1-y0)+(1.0d0+v)/2*(y2-y0)

  xyz(1) = 0
  xyz(2) = 0
  xyz(3) = 0

  dxyzds(1) = 0
  dxyzds(2) = 0
  dxyzds(3) = 0

  dxyzdt(1) = 0
  dxyzdt(2) = 0
  dxyzdt(3) = 0
  !call prin2('deltas = *', deltas, (m+2)*(n+2))
  !stop

  ct = cos(t)
  st = sin(t) 
  do i = -1,m
    do j = -1,n

      cst = cos((1.0d0-i)*s + j*t)
      sst = sin((1.0d0-i)*s + j*t)
      xyz(1) = xyz(1) + ct*deltas(i,j)*cst
      xyz(2) = xyz(2) + st*deltas(i,j)*cst
      xyz(3) = xyz(3) + deltas(i,j)*sst


      dxyzds(1) = dxyzds(1) - (1.0d0-i)*ct*deltas(i,j)*sst
      dxyzds(2) = dxyzds(2) - (1.0d0-i)*st*deltas(i,j)*sst
      dxyzds(3) = dxyzds(3) + (1.0d0-i)*deltas(i,j)*cst

      dxyzdt(1) = dxyzdt(1) + deltas(i,j)*(-st*cst -sst*ct*j)
      dxyzdt(2) = dxyzdt(2) + deltas(i,j)*(ct*cst - sst*st*j)
      dxyzdt(3) = dxyzdt(3) + deltas(i,j)*cst*j

    end do
  end do


  dsdu = (x1-x0)/2
  dsdv = (x2-x0)/2
  dtdu = (y1-y0)/2
  dtdv = (y2-y0)/2

  dxyzduv(1,1) = dxyzds(1)*dsdu + dxyzdt(1)*dtdu
  dxyzduv(2,1) = dxyzds(2)*dsdu + dxyzdt(2)*dtdu
  dxyzduv(3,1) = dxyzds(3)*dsdu + dxyzdt(3)*dtdu


  dxyzduv(1,2) = (dxyzds(1)*dsdv + dxyzdt(1)*dtdv)
  dxyzduv(2,2) = (dxyzds(2)*dsdv + dxyzdt(2)*dtdv)
  dxyzduv(3,2) = (dxyzds(3)*dsdv + dxyzdt(3)*dtdv)

  return
end subroutine xquad_stell_eval
!
!
!
!
!
!
!

subroutine xquad_xyz_tensor_fourier_eval(iquad, u, v, xyz, dxyzduv, &
    quadinfo, coefs, ipars, scales)
  implicit real *8 (a-h,o-z)
  integer ipars(2)
  real *8 :: xyz(3), dxyzduv(3,2), quadinfo(3,3,*), coefs(*)
  real *8 :: dxyzdst(3,2), scales(3)

  !
  ! project the quadrilateral iquad in quadinfo onto a double fourier toriodal
  ! surface
  !
  !    Input:
  ! iquad - quadrangle number to map
  ! u,v - local uv coordinates on quadrangle iquad
  ! quadinfo - flat skeleton quadrangle info
  !
  ! surface is given by
  !
  ! \hat(x) = \sum_{i=1}^{2m+1} \sum_{j=1} x_{ij} b_{i}(s) b_{j}(t)
  ! \hat(y) = \sum_{i=1}^{2m+1} \sum_{j=1} y_{ij} b_{i}(s) b_{j}(t)
  ! \hat(z) = \sum_{i=1}^{2m+1} \sum_{j=1} z_{ij} b_{i}(s) b_{j}(t)
  !
  ! x(s,t) = (\hat(x) \cos(s) - \hat(y) \sin(s))*scales(1)
  ! y(s,t) = (\hat(x) \sin(s) + \hat(y) \cos(s))*scales(2)
  ! z(s,t) = \hat(z)*scales(3)
  !
  !    Output:
  ! xyz - point on the 
  ! dxyzduv - first derivative information
  !
  !

  x0=quadinfo(1,1,iquad)
  y0=quadinfo(2,1,iquad)
  z0=quadinfo(3,1,iquad)

  x1=quadinfo(1,2,iquad)
  y1=quadinfo(2,2,iquad)
  z1=quadinfo(3,2,iquad)

  x2=quadinfo(1,3,iquad)
  y2=quadinfo(2,3,iquad)
  z2=quadinfo(3,3,iquad)


  !
  ! ... process the geometry, return the point location on the almond
  ! and the derivatives with respect to u and v
  !
  s = x0+(1.0d0+u)*(x1-x0)/2+(1.0d0+v)*(x2-x0)/2
  t = y0+(1.0d0+u)*(y1-y0)/2+(1.0d0+v)*(y2-y0)/2

  call xyz_tensor_fourier_eval(s, t, coefs, ipars(1), ipars(2), &
    scales, xyz, dxyzdst)

  dsdu = (x1-x0)/2
  dsdv = (x2-x0)/2
  dtdu = (y1-y0)/2
  dtdv = (y2-y0)/2


  dxyzduv(1,1) = dxyzdst(1,1)*dsdu + dxyzdst(1,2)*dtdu
  dxyzduv(2,1) = dxyzdst(2,1)*dsdu + dxyzdst(2,2)*dtdu
  dxyzduv(3,1) = dxyzdst(3,1)*dsdu + dxyzdst(3,2)*dtdu

  dxyzduv(1,2) = dxyzdst(1,1)*dsdv + dxyzdst(1,2)*dtdv
  dxyzduv(2,2) = dxyzdst(2,1)*dsdv + dxyzdst(2,2)*dtdv
  dxyzduv(3,2) = dxyzdst(3,1)*dsdv + dxyzdst(3,2)*dtdv

  return
end subroutine xquad_xyz_tensor_fourier_eval





subroutine xquad_ellipsoid_eval(iquad, u, v, xyz, dxyzduv, & 
  quadinfo,p2, p3, p4)

!
! project the quadrangle iquad in quadinfo onto the sphere
!
!    Input:
! iquad - quadrangle number to map
! u,v - local uv coordinates on quadrangle iquad
! quadinfo - flat skeleton quadrangle info
! p2(3) - semi major axes
! p3(3) - center of ellipsoid
! p4 - dummy parameters
!
!    Output:
! xyz - point on the sphere
! dxyzduv - first derivative information
!
!

  implicit real *8 (a-h,o-z)
  real *8 :: xyz(3), dxyzduv(3,2), quadinfo(3,3,*), p2(3), p3(3)

  x0=quadinfo(1,1,iquad)
  y0=quadinfo(2,1,iquad)
  z0=quadinfo(3,1,iquad)

  x1=quadinfo(1,2,iquad)
  y1=quadinfo(2,2,iquad)
  z1=quadinfo(3,2,iquad)

  x2=quadinfo(1,3,iquad)
  y2=quadinfo(2,3,iquad)
  z2=quadinfo(3,3,iquad)

!
! ... process the geometry, return the point location on the sphere
! and the derivatives with respect to u and v
!
  x=x0+(1.0d0+u)*(x1-x0)/2+(1.0d0+v)*(x2-x0)/2
  y=y0+(1.0d0+u)*(y1-y0)/2+(1.0d0+v)*(y2-y0)/2
  z=z0+(1.0d0+u)*(z1-z0)/2+(1.0d0+v)*(z2-z0)/2

  dxdu = (x1-x0)/2
  dydu = (y1-y0)/2
  dzdu = (z1-z0)/2
    
  dxdv = (x2-x0)/2
  dydv = (y2-y0)/2
  dzdv = (z2-z0)/2

!
! project onto the ellipsoid
!
  r=sqrt(x**2 + y**2 + z**2)
  xyz(1)=p2(1)*x/r + p3(1)
  xyz(2)=p2(2)*y/r + p3(2)
  xyz(3)=p2(3)*z/r + p3(3)

  drdu = (x*dxdu + y*dydu + z*dzdu)/r
  drdv = (x*dxdv + y*dydv + z*dzdv)/r


! du
  dxyzduv(1,1) = p2(1)*(r*dxdu-x*drdu)/r/r
  dxyzduv(2,1) = p2(2)*(r*dydu-y*drdu)/r/r
  dxyzduv(3,1) = p2(3)*(r*dzdu-z*drdu)/r/r

! dv
  dxyzduv(1,2) = p2(1)*(r*dxdv-x*drdv)/r/r
  dxyzduv(2,2) = p2(2)*(r*dydv-y*drdv)/r/r
  dxyzduv(3,2) = p2(3)*(r*dzdv-z*drdv)/r/r

  return
end subroutine xquad_ellipsoid_eval
!
!
!
!
!
subroutine xquad_sphere_eval(iquad, u, v, xyz, dxyzduv, quadinfo, &
    p2, p3, p4)
  implicit real *8 (a-h,o-z)
  real *8 :: xyz(3), dxyzduv(3,2), quadinfo(3,3,*)

  !
  ! project the quad iquad in quadinfo onto the sphere
  !
  !    Input:
  ! iquad - quad number to map
  ! u,v - local uv coordinates on quad iquad
  ! quadinfo - flat skeleton quad info
  ! p2,p3,p4 - dummy parameters
  !
  !    Output:
  ! xyz - point on the sphere
  ! dxyzduv - first derivative information
  !
  !

  x0=quadinfo(1,1,iquad)
  y0=quadinfo(2,1,iquad)
  z0=quadinfo(3,1,iquad)

  x1=quadinfo(1,2,iquad)
  y1=quadinfo(2,2,iquad)
  z1=quadinfo(3,2,iquad)

  x2=quadinfo(1,3,iquad)
  y2=quadinfo(2,3,iquad)
  z2=quadinfo(3,3,iquad)

  !
  ! ... process the geometry, return the point location on the sphere
  ! and the derivatives with respect to u and v
  !
  x=x0+(1.0d0+u)/2*(x1-x0)+(1.0d0+v)/2*(x2-x0)
  y=y0+(1.0d0+u)/2*(y1-y0)+(1.0d0+v)/2*(y2-y0)
  z=z0+(1.0d0+u)/2*(z1-z0)+(1.0d0+v)/2*(z2-z0)

  dxdu = (x1-x0)/2
  dydu = (y1-y0)/2
  dzdu = (z1-z0)/2

  dxdv = (x2-x0)/2
  dydv = (y2-y0)/2
  dzdv = (z2-z0)/2

  !
  ! second derivatives are zero...
  !

  !
  ! project onto the sphere
  !
  r=sqrt(x**2+y**2+z**2)
  xyz(1)=x/r
  xyz(2)=y/r
  xyz(3)=z/r

  a = x0*(x1-x0) + y0*(y1-y0) + z0*(z1-z0)
  b = (x1-x0)*(x2-x0) + (y1-y0)*(y2-y0) + (z1-z0)*(z2-z0)
  c = (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0)

  drdu = (x*dxdu+y*dydu+z*dzdu)/r


  drdv = (x*dxdv+y*dydv+z*dzdv)/r

  ! du
  dxyzduv(1,1) = (r*dxdu-x*drdu)/r/r
  dxyzduv(2,1) = (r*dydu-y*drdu)/r/r
  dxyzduv(3,1) = (r*dzdu-z*drdu)/r/r

  ! dv
  dxyzduv(1,2) = (r*dxdv-x*drdv)/r/r
  dxyzduv(2,2) = (r*dydv-y*drdv)/r/r
  dxyzduv(3,2) = (r*dzdv-z*drdv)/r/r

  return
end subroutine xquad_sphere_eval


!
!
!
subroutine xquad_axissym_eval(iquad, u, v, xyz, dxyzduv, & 
  quadinfo, srccoefs2d, ichuse, ixys2d)


!
!  Project iquad onto axissymmetric domain from a cylindrical
!  chart.
!
!    Input:
! iquad - triangle number to map
! u,v - local uv coordinates on triangle iquad
! quadinfo - flat skeleton triangle info
! srccoefs2d - coefs of points on chunks on generating curve 
! ichuse  - which chunk does the current triangle belong on 
! ixys2d - starting info for chunks 
!
!    Output:
! xyz - point on the sphere
! dxyzduv - first derivative information
!
!
  implicit real *8 (a-h,o-z)
  real *8 :: xyz(3), dxyzduv(3,2), quadinfo(3,3,*)
  real *8 :: pols(100), srccoefs2d(6,*)
  integer :: ichuse(*), ixys2d(*) 
      

  x0=quadinfo(1,1,iquad)
  y0=quadinfo(2,1,iquad)
  z0=quadinfo(3,1,iquad)

  x1=quadinfo(1,2,iquad)
  y1=quadinfo(2,2,iquad)
  z1=quadinfo(3,2,iquad)

  x2=quadinfo(1,3,iquad)
  y2=quadinfo(2,3,iquad)
  z2=quadinfo(3,3,iquad)

  s = x0 + (1.0d0+u)/2*(x1-x0) + (1.0d0+v)/2*(x2-x0)
  t = y0 + u*(y1-y0) + v*(y2-y0)

  dsdu = (x1-x0)/2
  dtdu = (y1-y0)/2
    
  dsdv = (x2-x0)/2
  dtdv = (y2-y0)/2
!
! s is coordinate in parametrization of r,z space
! t is coorindate in theta space
!


!
!  Compute r,z,drds,dzds
! 
!
  r = 0
  z = 0
  drds = 0
  dzds = 0
  pols = 0
  ich = ichuse(iquad)
  k = ixys2d(ich+1) - ixys2d(ich)
  call legepols(s, k-1, pols)
  do j=1,k
    r = r + srccoefs2d(1,ixys2d(ich)+j-1)*pols(j)
    z = z + srccoefs2d(2,ixys2d(ich)+j-1)*pols(j)
    drds = drds + srccoefs2d(3,ixys2d(ich)+j-1)*pols(j)
    dzds = dzds + srccoefs2d(4,ixys2d(ich)+j-1)*pols(j)
  enddo

  xyz(1)= r*cos(t)
  xyz(2)= r*sin(t)
  xyz(3)= z

! du
  dxyzduv(1,1) = drds*cos(t)*dsdu - r*sin(t)*dtdu
  dxyzduv(2,1) = drds*sin(t)*dsdu + r*cos(t)*dtdu 
  dxyzduv(3,1) = dzds*dsdu

! dv
  dxyzduv(1,2) = drds*cos(t)*dsdv - r*sin(t)*dtdv
  dxyzduv(2,2) = drds*sin(t)*dsdv + r*cos(t)*dtdv 
  dxyzduv(3,2) = dzds*dsdv

  return
end subroutine xquad_axissym_eval
!
!
!
!
!
subroutine xquad_axissym_circ_eval(iquad, u, v, xyz, dxyzduv, & 
  ptcoefs, srccoefs2d, ipars, pars)
!
!
!  Project iquad onto axissymmetric domain from a polar cap
!
!    Input:
! iquad - triangle number to map
! u,v - local uv coordinates on triangle iquad
! ptcoefs - flat skeleton triangle info
! srccoefs - basis function expansions for generating curve
! ipars - ipars(1) is starting location in srccoefs array
!         ipars(2) is number of points in expansion
!         ipars(3) is norder for ptcoefs evaluator
!         ipars(4) is iptype for ptcoefs evaluator
!         the rest is the ixys array corresponding to ptinfo
! pars(2)  - pars(1), pars(2) are rescaling of s
!         variable
!
!    Output:
! xyz - point on the sphere
! dxyzduv - first derivative information
!
!
  implicit real *8 (a-h,o-z)
  real *8 :: xyz(3), dxyzduv(3,2), ptcoefs(6,*), srccoefs2d(6,*)
  real *8 :: pols(500), pars(*)
  integer :: ipars(*)
  integer :: norder, iptype, istart
  real *8 :: uv(2)
  real *8 a0, b0

!
!  get values of s,t,dsdu,dtdu,dsdv,dtdv for
!  the given triangle on the circular patch
!
!
  uv(1) = u
  uv(2) = v

  istart2d = ipars(1)
  k = ipars(2)
  norder = ipars(3)
  iptype = ipars(4)
  
  call get_npols(iptype, norder, npols)
  call get_basis_pols(uv, norder, npols, iptype, pols)


  a0 = pars(1)
  b0 = pars(2)
  
  istart = ipars(4+iquad)
  s0 = 0
  t0 = 0
  ds0du = 0
  dt0du = 0
  ds0dv = 0
  dt0dv = 0
  do i = 1,npols
    s0 = s0 + ptcoefs(1,istart+i-1)*pols(i)
    t0 = t0 + ptcoefs(2,istart+i-1)*pols(i)
    ds0du = ds0du + ptcoefs(3,istart+i-1)*pols(i)
    dt0du = dt0du + ptcoefs(4,istart+i-1)*pols(i)
    ds0dv = ds0dv + ptcoefs(5,istart+i-1)*pols(i)
    dt0dv = dt0dv + ptcoefs(6,istart+i-1)*pols(i)
  enddo
  
  rr = sqrt(s0**2 + t0**2)
  s = a0 + (b0-a0)*rr
  t = atan2(t0, s0)


  drrdu = 1.0d0/rr*(s0*ds0du + t0*dt0du)
  drrdv = 1.0d0/rr*(s0*ds0dv + t0*dt0dv)
  
  dsdu = (b0-a0)*drrdu
  dsdv = (b0-a0)*drrdv
  
  dtdu = -t0/rr**2*ds0du + s0/rr**2*dt0du
  dtdv = -t0/rr**2*ds0dv + s0/rr**2*dt0dv

!
! s is coordinate in parametrization of r,z space
! t is coorindate in theta space
!

!
!  Compute r,z,drds,dzds
! 
!
  r = 0
  z = 0
  drds = 0
  dzds = 0

  call legepols(s, k-1, pols)
  do j=1,k
    r = r + srccoefs2d(1,istart2d+j-1)*pols(j)
    z = z + srccoefs2d(2,istart2d+j-1)*pols(j)
    drds = drds + srccoefs2d(3,istart2d+j-1)*pols(j)
    dzds = dzds + srccoefs2d(4,istart2d+j-1)*pols(j)
  enddo


  xyz(1)= r*cos(t)
  xyz(2)= r*sin(t)
  xyz(3)= z

! du
  dxyzduv(1,1) = drds*cos(t)*dsdu - r*sin(t)*dtdu
  dxyzduv(2,1) = drds*sin(t)*dsdu + r*cos(t)*dtdu 
  dxyzduv(3,1) = dzds*dsdu

! dv
  dxyzduv(1,2) = drds*cos(t)*dsdv - r*sin(t)*dtdv
  dxyzduv(2,2) = drds*sin(t)*dsdv + r*cos(t)*dtdv 
  dxyzduv(3,2) = dzds*dsdv

  return
end subroutine xquad_axissym_circ_eval
!
!
!
!
!
!
!
!
!
subroutine xquad_axissym_fun_eval(iquad, u, v, xyz, dxyzduv, & 
  quadinfo, np, pars, fcurve)
!
!  Project iquad onto axissymmetric domain from a cylindrical
!  chart.
!
!    Input:
! iquad - quadrilateral number to map
! u,v - local uv coordinates on quadrilateral iquad
! quadinfo - flat skeleton quadrilateral info
! np - number of parameters for function handle 
! pars(np)  - parameters for function handle 
! fcurve - function handle describing r(s), z(s) 
!
!    Output:
! xyz - point on the sphere
! dxyzduv - first derivative information
!
!
  implicit real *8 (a-h,o-z)
  real *8 :: xyz(3), dxyzduv(3,2), quadinfo(3,3,*)
  real *8 :: pols(100),pars(np)
  integer :: np
  external fcurve
      

  x0=quadinfo(1,1,iquad)
  y0=quadinfo(2,1,iquad)
  z0=quadinfo(3,1,iquad)

  x1=quadinfo(1,2,iquad)
  y1=quadinfo(2,2,iquad)
  z1=quadinfo(3,2,iquad)

  x2=quadinfo(1,3,iquad)
  y2=quadinfo(2,3,iquad)
  z2=quadinfo(3,3,iquad)

  s = x0 + (1.0d0+u)/2*(x1-x0) + (1.0d0+v)/2*(x2-x0)
  t = y0 + (1.0d0+u)/2*(y1-y0) + (1.0d0+v)/2*(y2-y0)

  dsdu = (x1-x0)/2
  dtdu = (y1-y0)/2
    
  dsdv = (x2-x0)/2
  dtdv = (y2-y0)/2
!
! s is coordinate in parametrization of r,z space
! t is coorindate in theta space
!


!
!  Compute r,z,drds,dzds
! 
!
      
  call fcurve(s, np, pars, r, z, drds, dzds, tmp1, tmp1)

  xyz(1)= r*cos(t)
  xyz(2)= r*sin(t)
  xyz(3)= z

! du
  dxyzduv(1,1) = drds*cos(t)*dsdu - r*sin(t)*dtdu
  dxyzduv(2,1) = drds*sin(t)*dsdu + r*cos(t)*dtdu 
  dxyzduv(3,1) = dzds*dsdu

! dv
  dxyzduv(1,2) = drds*cos(t)*dsdv - r*sin(t)*dtdv
  dxyzduv(2,2) = drds*sin(t)*dsdv + r*cos(t)*dtdv 
  dxyzduv(3,2) = dzds*dsdv

  return
end subroutine xquad_axissym_fun_eval
!
!
!
!
!
subroutine xquad_axissym_fun_circ_eval(iquad, u, v, xyz, dxyzduv, & 
  ptcoefs, ipars, pars, fcurve)
!
!  Project iquad onto axissymmetric domain from a polar cap
!
!    Input:
! iquad - quadrilateral number to map
! u,v - local uv coordinates on quadrilateral iquad
! ptcoefs - flat skeleton quadrilateral info
! ipars - ipars(1) is number of parameters for function handle,
!         ipars(2) is norder for ptcoefs evaluator
!         ipars(3) is iptype for ptcoefs evaluator
!         the rest is the ixys array corresponding to ptinfo
! pars(np+2)  - pars(1), pars(2) are rescaling of s
!         variable, and the rest are parameters for function handle 
! fcurve - function handle describing r(s), z(s) 
!
!    Output:
! xyz - point on the sphere
! dxyzduv - first derivative information
!
!
  implicit real *8 (a-h,o-z)
  real *8 :: xyz(3), dxyzduv(3,2), ptcoefs(6,*)
  real *8 :: pols(500), pars(*)
  integer :: ipars(*)
  integer :: np, norder, iptype, istart
  real *8 :: uv(2)
  external fcurve

!
!  get values of s,t,dsdu,dtdu,dsdv,dtdv for
!  the given quadrilateral on the circular patch
!
!

!
!  get values of s,t,dsdu,dtdu,dsdv,dtdv for
!  the given triangle on the circular patch
!
!
  uv(1) = u
  uv(2) = v
  np = ipars(1)
  norder = ipars(2)
  iptype = ipars(3)
  
  call get_npols(iptype, norder, npols)
  call get_basis_pols(uv, norder, npols, iptype, pols)

  a0 = pars(1)
  b0 = pars(2)
  
  istart = ipars(3+iquad)
  s0 = 0
  t0 = 0
  ds0du = 0
  dt0du = 0
  ds0dv = 0
  dt0dv = 0
  do i = 1,npols
    s0 = s0 + ptcoefs(1,istart+i-1)*pols(i)
    t0 = t0 + ptcoefs(2,istart+i-1)*pols(i)
    ds0du = ds0du + ptcoefs(3,istart+i-1)*pols(i)
    dt0du = dt0du + ptcoefs(4,istart+i-1)*pols(i)
    ds0dv = ds0dv + ptcoefs(5,istart+i-1)*pols(i)
    dt0dv = dt0dv + ptcoefs(6,istart+i-1)*pols(i)
  enddo

  rr = sqrt(s0**2 + t0**2)
  s = a0 + (b0-a0)*rr
  t = atan2(t0, s0)

  drrdu = 1.0d0/rr*(s0*ds0du + t0*dt0du)
  drrdv = 1.0d0/rr*(s0*ds0dv + t0*dt0dv)
  
  dsdu = (b0-a0)*drrdu
  dsdv = (b0-a0)*drrdv
  
  dtdu = -t0/rr**2*ds0du + s0/rr**2*dt0du
  dtdv = -t0/rr**2*ds0dv + s0/rr**2*dt0dv

!
! s is coordinate in parametrization of r,z space
! t is coorindate in theta space
!

!
!  Compute r,z,drds,dzds
! 
!
      
  call fcurve(s, np, pars(3), r, z, drds, dzds, tmp1, tmp1)

  xyz(1)= r*cos(t)
  xyz(2)= r*sin(t)
  xyz(3)= z

! du
  dxyzduv(1,1) = drds*cos(t)*dsdu - r*sin(t)*dtdu
  dxyzduv(2,1) = drds*sin(t)*dsdu + r*cos(t)*dtdu 
  dxyzduv(3,1) = dzds*dsdu

! dv
  dxyzduv(1,2) = drds*cos(t)*dsdv - r*sin(t)*dtdv
  dxyzduv(2,2) = drds*sin(t)*dsdv + r*cos(t)*dtdv 
  dxyzduv(3,2) = dzds*dsdv

  return
end subroutine xquad_axissym_fun_circ_eval
!
!
!
!

!
!
!
subroutine xquad_rectmesh_ani(umin, umax, vmin, vmax, nu, nv, &
    nover, maxquad, nquad, quadskel)
  implicit real *8 (a-h,o-z)
  real *8 :: quadskel(3,3,maxquad)

  real *8 :: verts(3,100), ifaces(3,100), verts1(3,3)
  real *8 :: verts2(3,3), verts3(3,3), verts4(3,3)

  !
  ! this routine returns a (possibly oversampled) quadngular mesh of
  ! the rectangle [umin,umax] \times [umin,vmax] beginning with nu
  ! divisions in u and nv divisions in v -- keep in mind the
  ! quads are stored with u,v,w components, with w=0 (in order to
  ! be consistent with other skeleton meshes)
  !
  !      Input:
  ! umin,umax, vmin,vmax - sets the dimensions of the rectangle
  ! nu, nv - number of original quads in each direction
  ! nover - determines how many times to oversample the skeleton
  !    quadngulation (every oversampling generates 4 time as many
  !    quads)
  ! maxquad - maximum number of allowable quads
  !
  !      Output:
  ! nquad - number of quads created
  ! quadskel - skeleton mesh info, basically just vertices of each
  ! quad
  !
  !

  width = umax-umin
  hu = width/nu
  height = vmax-vmin
  hv = height/nv

  !
  ! compute original quads on this grid
  !

  nquad = 0
  do i = 1,nu
    u0 = umin + (i-1)*hu
    u1 = u0 + hu
    do j = 1,nv
      v0 = vmin + (j-1)*hv
      v1 = v0 + hv
      call xquad_rectmesh0(u0, u1, v0, v1, quadskel(1,1,nquad+1))
      nquad = nquad + 1
    end do
  end do



  !
  ! go about oversampling the quads
  !
  do ijk = 1,nover

    nnn = nquad
    if (nnn*4 .gt. maxquad) then
      call prinf('maxquad exceeded, maxquad = *', maxquad, 1)
      stop
    end if

    do i = 1,nnn

      call xquad_refine4_flat(quadskel(1,1,i), verts1, &
          verts2, verts3, verts4)
      do j = 1,3
        do k = 1,3
          quadskel(k,j,i) = verts1(k,j)
          quadskel(k,j,nquad+1) = verts2(k,j)
          quadskel(k,j,nquad+2) = verts3(k,j)
          quadskel(k,j,nquad+3) = verts4(k,j)
        end do
      end do

      nquad = nquad +3
    end do
  end do

  return
end subroutine xquad_rectmesh_ani







subroutine xquad_rectmesh(umin, umax, vmin, vmax, nover, maxquad, &
    nquad, quadskel)
  implicit real *8 (a-h,o-z)
  real *8 :: quadskel(3,3,maxquad)

  real *8 :: verts(3,100), ifaces(3,100), verts1(3,3)
  real *8 :: verts2(3,3), verts3(3,3), verts4(3,3)

  !
  ! this routine returns a (possibly oversampled) quadngular mesh of
  ! the rectangle [0,umax] \times [0,vmax] -- keep in mind the
  ! quads are stored with u,v,w components, with w=0 (in order to
  ! be consistent with other skeleton meshes)
  !
  !      Input:
  ! umax, vmax - sets the dimensions of the rectangle
  ! nover - determines how many times to oversample the skeleton
  !    quadngulation (every oversampling generates 4 time as many
  !    quads)
  ! maxquad - maximum number of allowable quads
  !
  !      Output:
  ! nquad - number of quads created
  ! quadskel - skeleton mesh info, basically just vertices of each
  ! quad
  !
  !

  width = umax-umin
  height = vmax-vmin

  !
  ! split in v if taller than wide
  !
  if (height .ge. width) then

    irat = height/width
    hv = height/irat
    !call prin2('hv = *', hv, 1)

    nquad = 0
    do i = 1,irat
      v0 = vmin + (i-1)*hv
      v1 = v0 + hv
      call xquad_rectmesh0(umin, umax, v0, v1, quadskel(1,1,nquad+1))
      nquad = nquad + 1
    end do

  end if


  !
  ! split in u if wider than tall
  !
  if (height .lt. width) then

    irat = width/height
    hu = width/irat
    !call prin2('hu = *', hv, 1)

    nquad = 0
    do i = 1,irat
      u0 = umin + (i-1)*hu
      u1 = u0 + hu
      call xquad_rectmesh0(u0, u1, vmin, vmax, quadskel(1,1,nquad+1))
      nquad = nquad + 1
    end do

  end if


  !
  ! go about oversampling the quads
  !
  do ijk = 1,nover

    nnn = nquad
    if (nnn*4 .gt. maxquad) then
      call prinf('maxquad exceeded, maxquad = *', maxquad, 1)
      stop
    end if

    do i = 1,nnn

      call xquad_refine4_flat(quadskel(1,1,i), verts1, &
          verts2, verts3, verts4)
      do j = 1,3
        do k = 1,3
          quadskel(k,j,i) = verts1(k,j)
          quadskel(k,j,nquad+1) = verts2(k,j)
          quadskel(k,j,nquad+2) = verts3(k,j)
          quadskel(k,j,nquad+3) = verts4(k,j)
        end do
      end do

      nquad = nquad +3
    end do
  end do

  return
end subroutine xquad_rectmesh
!
!
!
!
!
subroutine xquad_rectmesh0(umin, umax, vmin, vmax, quadskel)
  implicit real *8 (a-h,o-z)
  real *8 :: quadskel(3,3)

  !
  ! this routine returns a skeleton mesh with two quads in
  ! the rectangle [umin,umax] \times [umin,vmax] -- keep in mind the
  ! quads are stored with u,v,w components, with w=0 (in order to
  ! be consistent with other skeleton meshes)
  !
  !      Input:
  ! umin, umax, vmin, vmax - sets the dimensions of the rectangle
  !
  !      Output:
  ! quadskel - skeleton mesh info, basically just vertices of each
  !     quad
  !
  !

  quadskel(1,1) = umin
  quadskel(2,1) = vmin
  quadskel(3,1) = 0
  


  quadskel(1,2) = umax
  quadskel(2,2) = vmin
  quadskel(3,2) = 0

  quadskel(1,3) = umin
  quadskel(2,3) = vmax
  quadskel(3,3) = 0


  return
end subroutine xquad_rectmesh0
!
!
!
!
!
subroutine xquad_rectmesh_3d(v1, v2, v3, v4, nu, nv, npatches, quadskel)
  implicit real *8 (a-h,o-z)
  real *8 quadskel(3,3,npatches), v1(3), v2(3), v3(3), v4(3)
  real *8 vl(3), vr(3), vb(3), vt(3)
  real *8 uvw1(3), uvw2(3), uvw3(3), uvw4(3)

  vl(1:3) = v4(1:3) - v1(1:3)
  vr(1:3) = v3(1:3) - v2(1:3)

  nquad = 0

  do i = 1, nv
    uvw1(1:3) = v1(1:3) + (i-1)*vl(1:3)/(nv+ 0.0d0)
    uvw4(1:3) = uvw1(1:3) + vl(1:3)/(nv+0.0d0)

    vb(1:3) =  v2(1:3) + (i-1)*vr(1:3)/(nv+0.0d0) - uvw1(1:3)
    vt(1:3) = uvw1(1:3) + vb(1:3) + vr(1:3)/(nv+0.0d0) - uvw4(1:3)
        
    uvw2(1:3) = uvw1(1:3) + vb(1:3)/(nu+0.0d0)
    uvw3(1:3) = uvw4(1:3) + vt(1:3)/(nu+0.0d0)

    do j = 1, nu
      call xquad_rectmesh0_3d(uvw1, uvw2, uvw3, uvw4, &
        quadskel(1,1,nquad+1))
      nquad = nquad + 1
      uvw1(1:3) = uvw2(1:3)
      uvw2(1:3) = uvw2(1:3) + vb(1:3)/(nu+0.0d0)
      uvw4(1:3) = uvw3(1:3)
      uvw3(1:3) = uvw3(1:3) + vt(1:3)/(nu+0.0d0)
    enddo
  enddo

end subroutine xquad_rectmesh_3d
!
!
!
!
!
subroutine xquad_rectmesh0_3d(v1, v2, v3, v4, quadskel)
  implicit real *8 (a-h,o-z)
  real *8 v1(3), v2(3), v3(3), v4(3), quadskel(3,3)

  do i=1,3
    quadskel(i,1) = v1(i)
    quadskel(i,2) = v2(i)
    quadskel(i,3) = v4(i)
  enddo

  return
end subroutine xquad_rectmesh0_3d
!
!
!
!
!
subroutine xquad_get_rectparapiped(a, b, c, na, nb, nc, &
  npatches, quadskel)

  implicit real *8 (a-h,o-z)
  real *8 quadskel(3,3,npatches),vs(3,4)
  real *8 vcube(3,8),xnorm(3)

      
      
  vcube(1,1) = -a
  vcube(2,1) = -b
  vcube(3,1) = -c

  vcube(1,2) = a
  vcube(2,2) = -b
  vcube(3,2) = -c

  vcube(1,3) = a
  vcube(2,3) = b
  vcube(3,3) = -c

  vcube(1,4) = -a
  vcube(2,4) = b
  vcube(3,4) = -c

  vcube(1,5) = -a
  vcube(2,5) = -b
  vcube(3,5) = c

  vcube(1,6) = a
  vcube(2,6) = -b
  vcube(3,6) = c

  vcube(1,7) = a
  vcube(2,7) = b
  vcube(3,7) = c

  vcube(1,8) = -a
  vcube(2,8) = b
  vcube(3,8) = c




!       z = -c face      
  vs(1:3,1) = vcube(1:3,1)
  vs(1:3,2) = vcube(1:3,4)
  vs(1:3,3) = vcube(1:3,3)
  vs(1:3,4) = vcube(1:3,2)
  nquad = 0
  call xquad_rectmesh_3d(vs(1,1),vs(1,2),vs(1,3),vs(1,4),nb,na, &
    npatches,quadskel(1,1,nquad+1))

  nquad = nquad + na*nb
      

!       z = c face      
  vs(1:3,1:4) = vcube(1:3,5:8)
  call xquad_rectmesh_3d(vs(1,1),vs(1,2),vs(1,3),vs(1,4), na, nb, &
    npatches,quadskel(1,1,nquad+1))

  nquad = nquad + na*nb

!      y = -b face
!
  vs(1:3,1) = vcube(1:3,1)
  vs(1:3,2) = vcube(1:3,2)
  vs(1:3,3) = vcube(1:3,6)
  vs(1:3,4) = vcube(1:3,5)

  call xquad_rectmesh_3d(vs(1,1),vs(1,2),vs(1,3),vs(1,4),na,nc, &
    npatches,quadskel(1,1,nquad+1))


  nquad = nquad + na*nc

!      y = b face
!
  vs(1:3,1) = vcube(1:3,4)
  vs(1:3,2) = vcube(1:3,8)
  vs(1:3,3) = vcube(1:3,7)
  vs(1:3,4) = vcube(1:3,3)

  call xquad_rectmesh_3d(vs(1,1),vs(1,2),vs(1,3),vs(1,4),nc,na, &
    npatches,quadskel(1,1,nquad+1))

  nquad = nquad + na*nc



!      x = -a face
!
  vs(1:3,1) = vcube(1:3,1)
  vs(1:3,2) = vcube(1:3,5)
  vs(1:3,3) = vcube(1:3,8)
  vs(1:3,4) = vcube(1:3,4)

  call xquad_rectmesh_3d(vs(1,1),vs(1,2),vs(1,3),vs(1,4),nc,nb, &
    npatches,quadskel(1,1,nquad+1))


  nquad = nquad + nb*nc

!      x = a face
!
  vs(1:3,1) = vcube(1:3,2)
  vs(1:3,2) = vcube(1:3,3)
  vs(1:3,3) = vcube(1:3,7)
  vs(1:3,4) = vcube(1:3,6)

  call xquad_rectmesh_3d(vs(1,1),vs(1,2),vs(1,3),vs(1,4),nb,nc, &
    npatches,quadskel(1,1,nquad+1))

  nquad = nquad + nb*nc


  return
end
!
!
!
!
!
!

subroutine get_norm_quadskel(quad, xnorm)
  implicit real *8 (a-h,o-z)
  real *8 quad(3,3), xnorm(3), xu(3), xv(3)
      

  xu(1:3) = quad(1:3,2) - quad(1:3,1)
  xv(1:3) = quad(1:3,3) - quad(1:3,1)

  xnorm(1:3) = 0
  call cross_prod3d(xu, xv, xnorm)


  return
end
!
!
!
!








subroutine xquad_platonic(itype, nover, maxquad, &
    nquads, quadinfo, isides)
  implicit real *8 (a-h,o-z)
  integer :: isides(*)
  real *8 :: quadinfo(3,3,maxquad)

  real *8 :: verts(3,100), ifaces(3,100), verts1(3,3)
  real *8 :: verts2(3,3), verts3(3,3), verts4(3,3)

  !
  ! this routine returns minimal (but possibly oversampled)
  ! quadngulations of the platonic solids (i.e. a cube that isn't
  ! oversampled has 12 quads)
  !
  !      Input:
  ! itype - determine the platonic solid to return
  !      1 - tetrahedron
  !      2 - cube
  !      3 - octahedron
  !      4 - icosahedron
  ! nover - determines how many times to oversample the skeleton
  !    quadngulation (every oversampling generates 4 time as many
  !    quads)
  ! maxquad - maximum number of allowable quads
  !
  !      Output:
  ! nquads - number of quads
  ! quadinfo - quad info
  ! isides - the side that the (possibly oversampled) quad
  !     belongs to (relative to original quadngulation), arbitrarily
  !     ordered.
  !
  !
  call xquad_rsolid(itype, verts, nverts, ifaces, nfaces)
  call xquad_genquadinfo(verts, nverts, ifaces, nfaces, quadinfo)
  nquads = nfaces

  do i = 1,nfaces
    isides(i) = i
  end do

  if (nover .eq. 0) then
    return
  end if

  !
  ! go about oversampling the quads
  !
  do ijk = 1,nover

    nnn = nquads
    if (nnn*4 .gt. maxquad) then
      call prinf('maxquad exceeded, maxquad = *', maxquad, 1)
      stop
    end if

    do i = 1,nnn
      call xquad_refine4_flat(quadinfo(1,1,i), verts1, &
          verts2, verts3, verts4)
      do j = 1,3
        do k = 1,3
          quadinfo(k,j,i) = verts1(k,j)
          quadinfo(k,j,nquads+1) = verts2(k,j)
          quadinfo(k,j,nquads+2) = verts3(k,j)
          quadinfo(k,j,nquads+3) = verts4(k,j)
        end do
        isides(nquads+1) = isides(i)
        isides(nquads+2) = isides(i)
        isides(nquads+3) = isides(i)
      end do

      nquads = nquads +3
    end do
  end do

  return
end subroutine xquad_platonic





subroutine xquad_rsolid(itype, verts, nverts, ifaces, nfaces)
  implicit real *8 (a-h,o-z)
  dimension verts(3,*),ifaces(3,*)
  !
  ! This subroutine returns the vertices and faces of regular (and
  ! not so regular) polyhedra.
  !
  !         Input:
  ! itype - determine the platonic solid to return
  !      2 - cube
  !

  if (itype .eq. 2) then

    !
    ! a cube
    !
    nverts=8
    nfaces=6

    ! ... vertices
    kk=1
    verts(1,kk)=-1
    verts(2,kk)=-1
    verts(3,kk)=-1

    kk=kk+1
    verts(1,kk)=+1
    verts(2,kk)=-1
    verts(3,kk)=-1

    kk=kk+1
    verts(1,kk)=+1
    verts(2,kk)=+1
    verts(3,kk)=-1

    kk=kk+1
    verts(1,kk)=-1
    verts(2,kk)=+1
    verts(3,kk)=-1

    kk=kk+1
    verts(1,kk)=-1
    verts(2,kk)=-1
    verts(3,kk)=+1

    kk=kk+1
    verts(1,kk)=+1
    verts(2,kk)=-1
    verts(3,kk)=+1

    kk=kk+1
    verts(1,kk)=+1
    verts(2,kk)=+1
    verts(3,kk)=+1

    kk=kk+1
    verts(1,kk)=-1
    verts(2,kk)=+1
    verts(3,kk)=+1

    ! ... faces
    kk=1
    ifaces(1,kk)=2
    ifaces(2,kk)=1
    ifaces(3,kk)=3


    kk=kk+1
    ifaces(1,kk)=1
    ifaces(2,kk)=2
    ifaces(3,kk)=5

    kk=kk+1
    ifaces(1,kk)=2
    ifaces(2,kk)=3
    ifaces(3,kk)=6

    kk=kk+1
    ifaces(1,kk)=4
    ifaces(2,kk)=8
    ifaces(3,kk)=3


    kk=kk+1
    ifaces(1,kk)=1
    ifaces(2,kk)=5
    ifaces(3,kk)=4

    kk=kk+1
    ifaces(1,kk)=5
    ifaces(2,kk)=6
    ifaces(3,kk)=8

  endif

  !
  ! scale the thing so that each edge has length 1
  !
  r=sqrt(verts(1,1)**2+verts(2,1)**2+verts(3,1)**2)

  do i=1,nverts
    verts(1,i)=verts(1,i)/r
    verts(2,i)=verts(2,i)/r
    verts(3,i)=verts(3,i)/r
  end do

  return
end subroutine xquad_rsolid





subroutine xquad_genquadinfo(verts,nverts,ifaces,nfaces,quadinfo)
  implicit real *8 (a-h,o-z)
  dimension verts(3,1),ifaces(3,1),quadinfo(3,3,1)

  call prinf('ifaces=*',ifaces,3*nfaces)

  do i=1,nfaces

    quadinfo(1,1,i)=verts(1,ifaces(1,i))
    quadinfo(2,1,i)=verts(2,ifaces(1,i))
    quadinfo(3,1,i)=verts(3,ifaces(1,i))

    quadinfo(1,2,i)=verts(1,ifaces(2,i))
    quadinfo(2,2,i)=verts(2,ifaces(2,i))
    quadinfo(3,2,i)=verts(3,ifaces(2,i))

    quadinfo(1,3,i)=verts(1,ifaces(3,i))
    quadinfo(2,3,i)=verts(2,ifaces(3,i))
    quadinfo(3,3,i)=verts(3,ifaces(3,i))

  end do

  return
end subroutine xquad_genquadinfo





subroutine xquad_refine4_flat(verts, verts1, verts2, verts3, verts4)
  implicit real *8 (a-h,o-z)
  real *8 :: verts(3,3), verts1(3,3), verts2(3,3), verts3(3,3)
  real *8 :: verts4(3,3)

  real *8 :: xyz12(3), xyz13(3), xyz34(3), xyz24(3), xyzm(3),v4(3)

  v4(1) = verts(1,2) + verts(1,3) - verts(1,1)
  v4(2) = verts(2,2) + verts(2,3) - verts(2,1)
  v4(3) = verts(3,2) + verts(3,3) - verts(3,1)


  !
  ! perform a refinement of a flat quad into four other flat
  ! quads
  !
  xyz12(1) = (verts(1,1) + verts(1,2))/2
  xyz12(2) = (verts(2,1) + verts(2,2))/2
  xyz12(3) = (verts(3,1) + verts(3,2))/2

  xyz13(1) = (verts(1,1) + verts(1,3))/2
  xyz13(2) = (verts(2,1) + verts(2,3))/2
  xyz13(3) = (verts(3,1) + verts(3,3))/2

  xyz34(1) = (v4(1) + verts(1,3))/2
  xyz34(2) = (v4(2) + verts(2,3))/2
  xyz34(3) = (v4(3) + verts(3,3))/2

  xyz24(1) = (verts(1,2) + v4(1))/2
  xyz24(2) = (verts(2,2) + v4(2))/2
  xyz24(3) = (verts(3,2) + v4(3))/2

  xyzm(1) = (verts(1,2) + verts(1,3))/2
  xyzm(2) = (verts(2,2) + verts(2,3))/2
  xyzm(3) = (verts(3,2) + verts(3,3))/2

  

  !
  ! first subdivision
  !
  verts1(1,1) = verts(1,1)
  verts1(2,1) = verts(2,1)
  verts1(3,1) = verts(3,1)

  verts1(1,2) = xyz12(1)
  verts1(2,2) = xyz12(2)
  verts1(3,2) = xyz12(3)

  verts1(1,3) = xyz13(1)
  verts1(2,3) = xyz13(2)
  verts1(3,3) = xyz13(3)

  !
  ! second subdivision
  !
  verts2(1,1) = xyz12(1)
  verts2(2,1) = xyz12(2)
  verts2(3,1) = xyz12(3)

  verts2(1,2) = verts(1,2)
  verts2(2,2) = verts(2,2)
  verts2(3,2) = verts(3,2)

  verts2(1,3) = xyzm(1)
  verts2(2,3) = xyzm(2)
  verts2(3,3) = xyzm(3)

  !
  ! third subdivision
  !
  verts3(1,1) = xyz13(1)
  verts3(2,1) = xyz13(2)
  verts3(3,1) = xyz13(3)

  verts3(1,2) = xyzm(1)
  verts3(2,2) = xyzm(2)
  verts3(3,2) = xyzm(3)

  verts3(1,3) = verts(1,3)
  verts3(2,3) = verts(2,3)
  verts3(3,3) = verts(3,3)

  !
  ! fourth subdivision
  !
  verts4(1,1) = xyzm(1)
  verts4(2,1) = xyzm(2)
  verts4(3,1) = xyzm(3)

  verts4(1,2) = xyz24(1)
  verts4(2,2) = xyz24(2)
  verts4(3,2) = xyz24(3)

  verts4(1,3) = xyz34(1)
  verts4(2,3) = xyz34(2)
  verts4(3,3) = xyz34(3)

  return
end subroutine xquad_refine4_flat
