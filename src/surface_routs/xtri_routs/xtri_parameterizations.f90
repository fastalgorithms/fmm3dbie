


subroutine xtri_wtorus_eval(itri, u, v, xyz, dxyzduv, triainfo, &
    radii, scales, p4)
  implicit real *8 (a-h,o-z)
  real *8 :: xyz(3), dxyzduv(3,2), triainfo(3,3,*), scales(3)
  real *8 :: radii(3)

  !
  ! project the triangle itri in triainfo onto a torus
  !
  !    Input:
  ! itri - triangle number to map
  ! u,v - local uv coordinates on triangle itri
  ! triainfo - flat skeleton triangle info
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

  x0=triainfo(1,1,itri)
  y0=triainfo(2,1,itri)
  z0=triainfo(3,1,itri)

  x1=triainfo(1,2,itri)
  y1=triainfo(2,2,itri)
  z1=triainfo(3,2,itri)

  x2=triainfo(1,3,itri)
  y2=triainfo(2,3,itri)
  z2=triainfo(3,3,itri)

  rminor = radii(1)
  rmajor = radii(2)
  rwave = radii(3)

  a = scales(1)
  b = scales(2)
  c = scales(3)

  nosc = p4


  !
  ! ... process the geometry, return the point location on the almond
  ! and the derivatives with respect to u and v
  !
  s = x0+u*(x1-x0)+v*(x2-x0)
  t = y0+u*(y1-y0)+v*(y2-y0)

  rr = rmajor+rminor*cos(t)+rwave*cos(nosc*s)


  xyz(1) = a*rr*cos(s)
  xyz(2) = b*rr*sin(s)
  xyz(3) = c*rminor*sin(t)

  dsdu = (x1-x0)
  dsdv = (x2-x0)
  dtdu = (y1-y0)
  dtdv = (y2-y0)

  drrds = -nosc*rwave*sin(nosc*s)
  drrdt = -rminor*sin(t)

  dxds = a*drrds*cos(s) - a*rr*sin(s)
  dyds = b*drrds*sin(s) + b*rr*cos(s)
  dzds = 0

  dxdt = a*drrdt*cos(s)
  dydt = b*drrdt*sin(s)
  dzdt = c*rminor*cos(t)
  
  dxdu = dxds*dsdu + dxdt*dtdu
  dydu = dyds*dsdu + dydt*dtdu
  dzdu = dzds*dsdu + dzdt*dtdu

  dxdv = dxds*dsdv + dxdt*dtdv
  dydv = dyds*dsdv + dydt*dtdv
  dzdv = dzds*dsdv + dzdt*dtdv

  dxyzduv(1,1) = dxdu
  dxyzduv(2,1) = dydu
  dxyzduv(3,1) = dzdu

  dxyzduv(1,2) = dxdv
  dxyzduv(2,2) = dydv
  dxyzduv(3,2) = dzdv

  return
end subroutine xtri_wtorus_eval
!
!
!
!
!
subroutine xtri_startorus_eval(itri, u, v, xyz, dxyzduv, triainfo, &
    radii, scales, p4)
  implicit real *8 (a-h,o-z)
  real *8 :: xyz(3), dxyzduv(3,2), triainfo(3,3,*), scales(3)
  real *8 :: radii(3)

  !
  ! project the triangle itri in triainfo onto a star-shaped torus
  !
  !    Input:
  ! itri - triangle number to map
  ! u,v - local uv coordinates on triangle itri
  ! triainfo - flat skeleton triangle info
  ! radii - the two radii defining the torus
  ! scales - scaling for x,y,z components from the standard torus
  ! p4 - number of oscillations (must be an integer currently recast
  !   as a double precision number)
  !
  !    Output:
  ! xyz - point on the sphere
  ! dxyzduv - first and second derivative information
  !
  !

  x0=triainfo(1,1,itri)
  y0=triainfo(2,1,itri)
  z0=triainfo(3,1,itri)

  x1=triainfo(1,2,itri)
  y1=triainfo(2,2,itri)
  z1=triainfo(3,2,itri)

  x2=triainfo(1,3,itri)
  y2=triainfo(2,3,itri)
  z2=triainfo(3,3,itri)

  rminor = radii(1)
  rmajor = radii(2)
  rwave = radii(3)

  a = scales(1)
  b = scales(2)
  c = scales(3)

  nosc = p4


  !
  ! ... process the geometry, return the point location on the almond
  ! and the derivatives with respect to u and v
  !
  s = x0+u*(x1-x0)+v*(x2-x0)
  t = y0+u*(y1-y0)+v*(y2-y0)

  rrs = rminor + rwave*cos(nosc*s)
  drrsds = -nosc*rwave*sin(nosc*s)

  rho = rmajor + rrs*cos(s)
  zz = rrs*sin(s)

  drhods = drrsds*cos(s) - rrs*sin(s)
  dzzds = drrsds*sin(s) + rrs*cos(s)


  xyz(1) = a*rho*cos(t)
  xyz(2) = b*rho*sin(t)
  xyz(3) = c*zz

  dsdu = (x1-x0)
  dsdv = (x2-x0)
  dtdu = (y1-y0)
  dtdv = (y2-y0)


  dxds = a*drhods*cos(t)
  dyds = b*drhods*sin(t)
  dzds = c*dzzds

  dxdt = -a*rho*sin(t)
  dydt = b*rho*cos(t)
  dzdt = 0 
  
  dxdu = dxds*dsdu + dxdt*dtdu
  dydu = dyds*dsdu + dydt*dtdu
  dzdu = dzds*dsdu + dzdt*dtdu

  dxdv = dxds*dsdv + dxdt*dtdv
  dydv = dyds*dsdv + dydt*dtdv
  dzdv = dzds*dsdv + dzdt*dtdv

  dxyzduv(1,1) = dxdu
  dxyzduv(2,1) = dydu
  dxyzduv(3,1) = dzdu

  dxyzduv(1,2) = dxdv
  dxyzduv(2,2) = dydv
  dxyzduv(3,2) = dzdv

  return
end subroutine xtri_startorus_eval








subroutine xtri_stell_eval(itri, u, v, xyz, dxyzduv, triainfo, &
    deltas, m, n)
  implicit real *8 (a-h,o-z)
  real *8 :: xyz(3), dxyzduv(3,2), triainfo(3,3,*), deltas(-1:m,-1:n)
  real *8 :: dxyzds(3),dxyzdt(3)

  !
  ! project the triangle itri in triainfo onto a stellarator
  !
  !    Input:
  ! itri - triangle number to map
  ! u,v - local uv coordinates on triangle itri
  ! triainfo - flat skeleton triangle info
  !
  !    Output:
  ! xyz - point on the stellarator
  ! dxyzduv - first and second derivative information
  !
  !


  x0=triainfo(1,1,itri)
  y0=triainfo(2,1,itri)
  z0=triainfo(3,1,itri)

  x1=triainfo(1,2,itri)
  y1=triainfo(2,2,itri)
  z1=triainfo(3,2,itri)

  x2=triainfo(1,3,itri)
  y2=triainfo(2,3,itri)
  z2=triainfo(3,3,itri)


  !
  ! ... process the geometry, return the point location on the almond
  ! and the derivatives with respect to u and v
  !
  s = x0+u*(x1-x0)+v*(x2-x0)
  t = y0+u*(y1-y0)+v*(y2-y0)

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


  dsdu = (x1-x0)
  dsdv = (x2-x0)
  dtdu = (y1-y0)
  dtdv = (y2-y0)

  dxyzduv(1,1) = dxyzds(1)*dsdu + dxyzdt(1)*dtdu
  dxyzduv(2,1) = dxyzds(2)*dsdu + dxyzdt(2)*dtdu
  dxyzduv(3,1) = dxyzds(3)*dsdu + dxyzdt(3)*dtdu


  dxyzduv(1,2) = (dxyzds(1)*dsdv + dxyzdt(1)*dtdv)
  dxyzduv(2,2) = (dxyzds(2)*dsdv + dxyzdt(2)*dtdv)
  dxyzduv(3,2) = (dxyzds(3)*dsdv + dxyzdt(3)*dtdv)

  return
end subroutine xtri_stell_eval
!
!
!
!
!

subroutine xtri_xyz_tensor_fourier_eval(itri, u, v, xyz, dxyzduv, &
    triainfo, coefs, m, scales)
  implicit real *8 (a-h,o-z)
  real *8 :: xyz(3), dxyzduv(3,2), triainfo(3,3,*), coefs(2*m+1,2*m+1,3)
  real *8 :: dxyzds(3), dxyzdt(3)
  real *8 :: hatxyz(3), dhatxyzds(3), dhatxyzdt(3)
  real *8 :: bis(2*m+1), bjs(2*m+1), bdis(2*m+1), bdjs(2*m+1)
  real *8 :: scales(3)
  complex *16 zmuls, zmult, ima, zfacs, zfact
  data ima/(0.0d0,1.0d0)/

  !
  ! project the triangle itri in triainfo onto a stellarator
  !
  !    Input:
  ! itri - triangle number to map
  ! u,v - local uv coordinates on triangle itri
  ! triainfo - flat skeleton triangle info
  !
  ! surface is given by
  !
  ! \hat(x) = \sum_{i=1}^{2m+1} \sum_{j=1} x_{ij} b_{i}(s) b_{j}(t)
  ! \hat(y) = \sum_{i=1}^{2m+1} \sum_{j=1} y_{ij} b_{i}(s) b_{j}(t)
  ! \hat(z) = \sum_{i=1}^{2m+1} \sum_{j=1} z_{ij} b_{i}(s) b_{j}(t)
  !
  ! x(s,t) = (\hat(x) \cos(t) - \hat(y) \sin(t))*scales(1)
  ! y(s,t) = (\hat(x) \sin(t) + \hat(y) \cos(t))*scales(2)
  ! z(s,t) = \hat(z)*scales(3)
  !
  !    Output:
  ! xyz - point on the 
  ! dxyzduv - first derivative information
  !
  !

  x0=triainfo(1,1,itri)
  y0=triainfo(2,1,itri)
  z0=triainfo(3,1,itri)

  x1=triainfo(1,2,itri)
  y1=triainfo(2,2,itri)
  z1=triainfo(3,2,itri)

  x2=triainfo(1,3,itri)
  y2=triainfo(2,3,itri)
  z2=triainfo(3,3,itri)


  !
  ! ... process the geometry, return the point location on the almond
  ! and the derivatives with respect to u and v
  !
  s = x0+u*(x1-x0)+v*(x2-x0)
  t = y0+u*(y1-y0)+v*(y2-y0)

  xyz(1) = 0
  xyz(2) = 0
  xyz(3) = 0

  hatxyz(1:3) = 0

  dxyzds(1) = 0
  dxyzds(2) = 0
  dxyzds(3) = 0

  dhatxyzds(1:3) = 0

  dxyzdt(1) = 0
  dxyzdt(2) = 0
  dxyzdt(3) = 0

  dhatxyzdt(1:3) = 0
  
  ct = cos(t)
  st = sin(t)
  zmuls = cos(s) + ima*sin(s)
  zmult = ct + ima*st


  zfacs = zmuls
  zfact = zmult
  
  bis(1) = 1
  bdis(1) = 0
  
  bjs(1) = 1
  bdjs(1) = 0

  do i=1,m
     bis(i+1) = real(zmuls)
     bis(i+1+m) = imag(zmuls)

     bdis(i+1) = -i*imag(zmuls)
     bdis(i+m+1) = i*real(zmuls)
     
     bjs(i+1) = real(zmult)
     bjs(i+1+m) = imag(zmult)

     bdjs(i+1) = -i*imag(zmult)
     bdjs(i+1+m) = i*real(zmult)


     zmuls = zmuls*zfacs
     zmult = zmult*zfact
  enddo
!
!
  do j = 1,2*m+1
    do i = 1,2*m+1
      hatxyz(1) = hatxyz(1) + coefs(i,j,1)*bis(i)*bjs(j)
      hatxyz(2) = hatxyz(2) + coefs(i,j,2)*bis(i)*bjs(j)
      hatxyz(3) = hatxyz(3) + coefs(i,j,3)*bis(i)*bjs(j)

      dhatxyzds(1) = dhatxyzds(1) + coefs(i,j,1)*bdis(i)*bjs(j)
      dhatxyzds(2) = dhatxyzds(2) + coefs(i,j,2)*bdis(i)*bjs(j)
      dhatxyzds(3) = dhatxyzds(3) + coefs(i,j,3)*bdis(i)*bjs(j)

      dhatxyzdt(1) = dhatxyzdt(1) + coefs(i,j,1)*bis(i)*bdjs(j)
      dhatxyzdt(2) = dhatxyzdt(2) + coefs(i,j,2)*bis(i)*bdjs(j)
      dhatxyzdt(3) = dhatxyzdt(3) + coefs(i,j,3)*bis(i)*bdjs(j)
    enddo
  enddo

  xyz(1) = (hatxyz(1)*ct - hatxyz(2)*st)*scales(1)
  xyz(2) = (hatxyz(1)*st + hatxyz(2)*ct)*scales(2)
  xyz(3) = hatxyz(3)*scales(3)
  
  dxyzds(1) = (dhatxyzds(1)*ct - dhatxyzds(2)*st)*scales(1) 
  dxyzds(2) = (dhatxyzds(1)*st + dhatxyzds(2)*ct)*scales(2) 
  dxyzds(3) = dhatxyzds(3)*scales(3)

  dxyzdt(1) = -(hatxyz(1)*st + hatxyz(2)*ct)
  dxyzdt(1) = dxyzdt(1) + (dhatxyzdt(1)*ct - dhatxyzdt(2)*st)
  dxyzdt(1) = dxyzdt(1)*scales(1)

  dxyzdt(2) = (hatxyz(1)*ct - hatxyz(2)*st)
  dxyzdt(2) = dxyzdt(2) + (dhatxyzdt(1)*st + dhatxyzdt(2)*ct)
  dxyzdt(2) = dxyzdt(2)*scales(2)

  dxyzdt(3) = dhatxyzdt(3)*scales(3)

  

  dsdu = (x1-x0)
  dsdv = (x2-x0)
  dtdu = (y1-y0)
  dtdv = (y2-y0)

  dxyzduv(1,1) = dxyzds(1)*dsdu + dxyzdt(1)*dtdu
  dxyzduv(2,1) = dxyzds(2)*dsdu + dxyzdt(2)*dtdu
  dxyzduv(3,1) = dxyzds(3)*dsdu + dxyzdt(3)*dtdu


  dxyzduv(1,2) = dxyzds(1)*dsdv + dxyzdt(1)*dtdv
  dxyzduv(2,2) = dxyzds(2)*dsdv + dxyzdt(2)*dtdv
  dxyzduv(3,2) = dxyzds(3)*dsdv + dxyzdt(3)*dtdv

  return
end subroutine xtri_xyz_tensor_fourier_eval





subroutine xtri_ellipsoid_eval(itri, u, v, xyz, dxyzduv, & 
  triainfo,p2, p3, p4)

!
! project the triangle itri in triainfo onto the sphere
!
!    Input:
! itri - triangle number to map
! u,v - local uv coordinates on triangle itri
! triainfo - flat skeleton triangle info
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
  real *8 :: xyz(3), dxyzduv(3,2), triainfo(3,3,*), p2(3), p3(3)

  x0=triainfo(1,1,itri)
  y0=triainfo(2,1,itri)
  z0=triainfo(3,1,itri)

  x1=triainfo(1,2,itri)
  y1=triainfo(2,2,itri)
  z1=triainfo(3,2,itri)

  x2=triainfo(1,3,itri)
  y2=triainfo(2,3,itri)
  z2=triainfo(3,3,itri)

!
! ... process the geometry, return the point location on the sphere
! and the derivatives with respect to u and v
!
  x=x0+u*(x1-x0)+v*(x2-x0)
  y=y0+u*(y1-y0)+v*(y2-y0)
  z=z0+u*(z1-z0)+v*(z2-z0)

  dxdu = x1-x0
  dydu = y1-y0
  dzdu = z1-z0
    
  dxdv = x2-x0
  dydv = y2-y0
  dzdv = z2-z0

!
! project onto the ellipsoid
!
  r=sqrt(x**2 + y**2 + z**2)
  xyz(1)=p2(1)*x/r + p3(1)
  xyz(2)=p2(2)*y/r + p3(2)
  xyz(3)=p2(3)*z/r + p3(3)

  a = x0*(x1-x0) + y0*(y1-y0) + z0*(z1-z0)
  b = (x1-x0)*(x2-x0) + (y1-y0)*(y2-y0) + (z1-z0)*(z2-z0)
  c = (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0)

  drdu = (a + v*b + u*c)/r
  drdu2 = (r*c - r*drdu*drdu)/r/r

  e = x0*(x2-x0) + y0*(y2-y0) + z0*(z2-z0)
  f = b
  g = (x2-x0)*(x2-x0) + (y2-y0)*(y2-y0) + (z2-z0)*(z2-z0)

  drdv = (e + u*f + v*g)/r
  drdv2 = (r*g - r*drdv*drdv)/r/r

  drduv = (r*b - r*drdu*drdv)/r/r

! du
  dxyzduv(1,1) = p2(1)*(r*dxdu-x*drdu)/r/r
  dxyzduv(2,1) = p2(2)*(r*dydu-y*drdu)/r/r
  dxyzduv(3,1) = p2(3)*(r*dzdu-z*drdu)/r/r

! dv
  dxyzduv(1,2) = p2(1)*(r*dxdv-x*drdv)/r/r
  dxyzduv(2,2) = p2(2)*(r*dydv-y*drdv)/r/r
  dxyzduv(3,2) = p2(3)*(r*dzdv-z*drdv)/r/r

  return
end subroutine xtri_ellipsoid_eval
!
!
!
!
!
subroutine xtri_rectmesh_ani(umin, umax, vmin, vmax, nu, nv, &
    nover, maxtri, ntri, triaskel)
  implicit real *8 (a-h,o-z)
  real *8 :: triaskel(3,3,maxtri)

  real *8 :: verts(3,100), ifaces(3,100), verts1(3,3)
  real *8 :: verts2(3,3), verts3(3,3), verts4(3,3)

  !
  ! this routine returns a (possibly oversampled) triangular mesh of
  ! the rectangle [umin,umax] \times [umin,vmax] beginning with nu
  ! divisions in u and nv divisions in v -- keep in mind the
  ! triangles are stored with u,v,w components, with w=0 (in order to
  ! be consistent with other skeleton meshes)
  !
  !      Input:
  ! umin,umax, vmin,vmax - sets the dimensions of the rectangle
  ! nu, nv - number of original triangles in each direction
  ! nover - determines how many times to oversample the skeleton
  !    triangulation (every oversampling generates 4 time as many
  !    triangles)
  ! maxtri - maximum number of allowable triangles
  !
  !      Output:
  ! ntri - number of triangles created
  ! triaskel - skeleton mesh info, basically just vertices of each
  ! triangle
  !
  !

  width = umax-umin
  hu = width/nu
  height = vmax-vmin
  hv = height/nv

  !
  ! compute original triangles on this grid
  !

  ntri = 0
  do i = 1,nu
    u0 = umin + (i-1)*hu
    u1 = u0 + hu
    do j = 1,nv
      v0 = vmin + (j-1)*hv
      v1 = v0 + hv
      call xtri_rectmesh0(u0, u1, v0, v1, triaskel(1,1,ntri+1))
      ntri = ntri + 2
    end do
  end do



  !
  ! go about oversampling the triangles
  !
  do ijk = 1,nover

    nnn = ntri
    if (nnn*4 .gt. maxtri) then
      call prinf('maxtri exceeded, maxtri = *', maxtri, 1)
      stop
    end if

    do i = 1,nnn

      call xtri_refine4_flat(triaskel(1,1,i), verts1, &
          verts2, verts3, verts4)
      do j = 1,3
        do k = 1,3
          triaskel(k,j,i) = verts1(k,j)
          triaskel(k,j,ntri+1) = verts2(k,j)
          triaskel(k,j,ntri+2) = verts3(k,j)
          triaskel(k,j,ntri+3) = verts4(k,j)
        end do
      end do

      ntri = ntri +3
    end do
  end do

  return
end subroutine xtri_rectmesh_ani







subroutine xtri_rectmesh(umin, umax, vmin, vmax, nover, maxtri, &
    ntri, triaskel)
  implicit real *8 (a-h,o-z)
  real *8 :: triaskel(3,3,maxtri)

  real *8 :: verts(3,100), ifaces(3,100), verts1(3,3)
  real *8 :: verts2(3,3), verts3(3,3), verts4(3,3)

  !
  ! this routine returns a (possibly oversampled) triangular mesh of
  ! the rectangle [0,umax] \times [0,vmax] -- keep in mind the
  ! triangles are stored with u,v,w components, with w=0 (in order to
  ! be consistent with other skeleton meshes)
  !
  !      Input:
  ! umax, vmax - sets the dimensions of the rectangle
  ! nover - determines how many times to oversample the skeleton
  !    triangulation (every oversampling generates 4 time as many
  !    triangles)
  ! maxtri - maximum number of allowable triangles
  !
  !      Output:
  ! ntri - number of triangles created
  ! triaskel - skeleton mesh info, basically just vertices of each
  ! triangle
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

    ntri = 0
    do i = 1,irat
      v0 = vmin + (i-1)*hv
      v1 = v0 + hv
      call xtri_rectmesh0(umin, umax, v0, v1, triaskel(1,1,ntri+1))
      ntri = ntri + 2
    end do

  end if


  !
  ! split in u if wider than tall
  !
  if (height .lt. width) then

    irat = width/height
    hu = width/irat
    !call prin2('hu = *', hv, 1)

    ntri = 0
    do i = 1,irat
      u0 = umin + (i-1)*hu
      u1 = u0 + hu
      call xtri_rectmesh0(u0, u1, vmin, vmax, triaskel(1,1,ntri+1))
      ntri = ntri + 2
    end do

  end if


  !
  ! go about oversampling the triangles
  !
  do ijk = 1,nover

    nnn = ntri
    if (nnn*4 .gt. maxtri) then
      call prinf('maxtri exceeded, maxtri = *', maxtri, 1)
      stop
    end if

    do i = 1,nnn

      call xtri_refine4_flat(triaskel(1,1,i), verts1, &
          verts2, verts3, verts4)
      do j = 1,3
        do k = 1,3
          triaskel(k,j,i) = verts1(k,j)
          triaskel(k,j,ntri+1) = verts2(k,j)
          triaskel(k,j,ntri+2) = verts3(k,j)
          triaskel(k,j,ntri+3) = verts4(k,j)
        end do
      end do

      ntri = ntri +3
    end do
  end do

  return
end subroutine xtri_rectmesh
!
!
!
!
!
subroutine xtri_rectmesh0(umin, umax, vmin, vmax, triaskel)
  implicit real *8 (a-h,o-z)
  real *8 :: triaskel(3,3,2)

  !
  ! this routine returns a skeleton mesh with two triangles in
  ! the rectangle [umin,umax] \times [umin,vmax] -- keep in mind the
  ! triangles are stored with u,v,w components, with w=0 (in order to
  ! be consistent with other skeleton meshes)
  !
  !      Input:
  ! umin, umax, vmin, vmax - sets the dimensions of the rectangle
  !
  !      Output:
  ! triaskel - skeleton mesh info, basically just vertices of each
  !     triangle
  !
  !

  triaskel(1,1,1) = umin
  triaskel(2,1,1) = vmin
  triaskel(1,2,1) = umax
  triaskel(2,2,1) = vmin
  triaskel(1,3,1) = umin
  triaskel(2,3,1) = vmax

  triaskel(1,1,2) = umax
  triaskel(2,1,2) = vmax
  triaskel(1,2,2) = umin
  triaskel(2,2,2) = vmax
  triaskel(1,3,2) = umax
  triaskel(2,3,2) = vmin

  ntri = 2
  do n = 1,ntri
    do i = 1,3
      triaskel(3,i,n) = 0
    end do
  end do

  return
end subroutine xtri_rectmesh0
!
!
!
!
!
subroutine xtri_rectmesh_3d(v1, v2, v3, v4, nu, nv, npatches, triaskel)
  implicit real *8 (a-h,o-z)
  real *8 triaskel(3,3,npatches), v1(3), v2(3), v3(3), v4(3)
  real *8 vl(3), vr(3), vb(3), vt(3)
  real *8 uvw1(3), uvw2(3), uvw3(3), uvw4(3)

  vl(1:3) = v4(1:3) - v1(1:3)
  vr(1:3) = v3(1:3) - v2(1:3)

  ntri = 0

  do i=1,nv
    uvw1(1:3) = v1(1:3) + (i-1)*vl(1:3)/(nv+ 0.0d0)
    uvw4(1:3) = uvw1(1:3) + vl(1:3)/(nv+0.0d0)

    vb(1:3) =  v2(1:3) + (i-1)*vr(1:3)/(nv+0.0d0) - uvw1(1:3)
    vt(1:3) = uvw1(1:3) + vb(1:3) + vr(1:3)/(nv+0.0d0) - uvw4(1:3)
        
    uvw2(1:3) = uvw1(1:3) + vb(1:3)/(nu+0.0d0)
    uvw3(1:3) = uvw4(1:3) + vt(1:3)/(nu+0.0d0)

    do j=1,nu
      call xtri_rectmesh0_3d(uvw1, uvw2, uvw3, uvw4, &
        triaskel(1,1,ntri+1))
      ntri = ntri + 2
      uvw1(1:3) = uvw2(1:3)
      uvw2(1:3) = uvw2(1:3) + vb(1:3)/(nu+0.0d0)
      uvw4(1:3) = uvw3(1:3)
      uvw3(1:3) = uvw3(1:3) + vt(1:3)/(nu+0.0d0)
    enddo
  enddo

end subroutine xtri_rectmesh_3d
!
!
!
!
!
subroutine xtri_rectmesh0_3d(v1, v2, v3, v4, triaskel)
  implicit real *8 (a-h,o-z)
  real *8 v1(3), v2(3), v3(3), v4(3), triaskel(3,3,2)

  do i=1,3
    triaskel(i,1,1) = v1(i)
    triaskel(i,2,1) = v2(i)
    triaskel(i,3,1) = v4(i)
    triaskel(i,1,2) = v3(i)
    triaskel(i,2,2) = v4(i)
    triaskel(i,3,2) = v2(i)
  enddo

  return
end subroutine xtri_rectmesh0_3d
!
!
!
!
subroutine xtri_get_rectparapiped(a, b, c, na, nb, nc, &
  npatches, triaskel)

  implicit real *8 (a-h,o-z)
  real *8 triaskel(3,3,npatches),vs(3,4)
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
  ntri = 0
  call xtri_rectmesh_3d(vs(1,1),vs(1,2),vs(1,3),vs(1,4),nb,na, &
    npatches,triaskel(1,1,ntri+1))

  call get_norm_triaskel(triaskel(1,1,ntri+1),xnorm)
  call get_norm_triaskel(triaskel(1,1,ntri+2),xnorm)

  ntri = ntri + 2*na*nb
      

!       z = c face      
  vs(1:3,1:4) = vcube(1:3,5:8)
  call xtri_rectmesh_3d(vs(1,1),vs(1,2),vs(1,3),vs(1,4),na,nb, &
    npatches,triaskel(1,1,ntri+1))

  call get_norm_triaskel(triaskel(1,1,ntri+1),xnorm)
  call get_norm_triaskel(triaskel(1,1,ntri+2),xnorm)

  ntri = ntri + 2*na*nb

!      y = -b face
!
  vs(1:3,1) = vcube(1:3,1)
  vs(1:3,2) = vcube(1:3,2)
  vs(1:3,3) = vcube(1:3,6)
  vs(1:3,4) = vcube(1:3,5)

  call xtri_rectmesh_3d(vs(1,1),vs(1,2),vs(1,3),vs(1,4),na,nc, &
    npatches,triaskel(1,1,ntri+1))
  call get_norm_triaskel(triaskel(1,1,ntri+1),xnorm)
  call get_norm_triaskel(triaskel(1,1,ntri+2),xnorm)


  ntri = ntri + 2*na*nc

!      y = b face
!
  vs(1:3,1) = vcube(1:3,4)
  vs(1:3,2) = vcube(1:3,8)
  vs(1:3,3) = vcube(1:3,7)
  vs(1:3,4) = vcube(1:3,3)

  call xtri_rectmesh_3d(vs(1,1),vs(1,2),vs(1,3),vs(1,4),nc,na, &
    npatches,triaskel(1,1,ntri+1))
  call get_norm_triaskel(triaskel(1,1,ntri+1),xnorm)
  call get_norm_triaskel(triaskel(1,1,ntri+2),xnorm)

  ntri = ntri + 2*na*nc



!      x = -a face
!
  vs(1:3,1) = vcube(1:3,1)
  vs(1:3,2) = vcube(1:3,5)
  vs(1:3,3) = vcube(1:3,8)
  vs(1:3,4) = vcube(1:3,4)

  call xtri_rectmesh_3d(vs(1,1),vs(1,2),vs(1,3),vs(1,4),nc,nb, &
    npatches,triaskel(1,1,ntri+1))
  call get_norm_triaskel(triaskel(1,1,ntri+1),xnorm)
  call get_norm_triaskel(triaskel(1,1,ntri+2),xnorm)


  ntri = ntri + 2*nb*nc

!      x = a face
!
  vs(1:3,1) = vcube(1:3,2)
  vs(1:3,2) = vcube(1:3,3)
  vs(1:3,3) = vcube(1:3,7)
  vs(1:3,4) = vcube(1:3,6)

  call xtri_rectmesh_3d(vs(1,1),vs(1,2),vs(1,3),vs(1,4),nb,nc, &
    npatches,triaskel(1,1,ntri+1))
  call get_norm_triaskel(triaskel(1,1,ntri+1),xnorm)
  call get_norm_triaskel(triaskel(1,1,ntri+2),xnorm)

  ntri = ntri + 2*nb*nc


  return
end
!
!
!
!
!
!

subroutine get_norm_triaskel(tria, xnorm)
  implicit real *8 (a-h,o-z)
  real *8 tria(3,3), xnorm(3), xu(3), xv(3)
      

  xu(1:3) = tria(1:3,2) - tria(1:3,1)
  xv(1:3) = tria(1:3,3) - tria(1:3,1)

  xnorm(1:3) = 0
  call cross_prod3d(xu, xv, xnorm)


  return
end
!
!
!
!
subroutine xtri_sphere_eval(itri, u, v, xyz, dxyzduv, triainfo, &
    p2, p3, p4)
  implicit real *8 (a-h,o-z)
  real *8 :: xyz(3), dxyzduv(3,2), triainfo(3,3,*)

  !
  ! project the triangle itri in triainfo onto the sphere
  !
  !    Input:
  ! itri - triangle number to map
  ! u,v - local uv coordinates on triangle itri
  ! triainfo - flat skeleton triangle info
  ! p2,p3,p4 - dummy parameters
  !
  !    Output:
  ! xyz - point on the sphere
  ! dxyzduv - first derivative information
  !
  !

  x0=triainfo(1,1,itri)
  y0=triainfo(2,1,itri)
  z0=triainfo(3,1,itri)

  x1=triainfo(1,2,itri)
  y1=triainfo(2,2,itri)
  z1=triainfo(3,2,itri)

  x2=triainfo(1,3,itri)
  y2=triainfo(2,3,itri)
  z2=triainfo(3,3,itri)

  !
  ! ... process the geometry, return the point location on the sphere
  ! and the derivatives with respect to u and v
  !
  x=x0+u*(x1-x0)+v*(x2-x0)
  y=y0+u*(y1-y0)+v*(y2-y0)
  z=z0+u*(z1-z0)+v*(z2-z0)

  dxdu = x1-x0
  dydu = y1-y0
  dzdu = z1-z0

  dxdv = x2-x0
  dydv = y2-y0
  dzdv = z2-z0

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

  drdu = (a + v*b + u*c)/r
  drdu2 = (r*c - r*drdu*drdu)/r/r

  e = x0*(x2-x0) + y0*(y2-y0) + z0*(z2-z0)
  f = b
  g = (x2-x0)*(x2-x0) + (y2-y0)*(y2-y0) + (z2-z0)*(z2-z0)

  drdv = (e + u*f + v*g)/r
  drdv2 = (r*g - r*drdv*drdv)/r/r

  drduv = (r*b - r*drdu*drdv)/r/r

  ! du
  dxyzduv(1,1) = (r*dxdu-x*drdu)/r/r
  dxyzduv(2,1) = (r*dydu-y*drdu)/r/r
  dxyzduv(3,1) = (r*dzdu-z*drdu)/r/r

  ! dv
  dxyzduv(1,2) = (r*dxdv-x*drdv)/r/r
  dxyzduv(2,2) = (r*dydv-y*drdv)/r/r
  dxyzduv(3,2) = (r*dzdv-z*drdv)/r/r

  return
end subroutine xtri_sphere_eval








subroutine xtri_platonic(itype, nover, maxtri, &
    npatches, triainfo, isides)
  implicit real *8 (a-h,o-z)
  integer :: isides(*)
  real *8 :: triainfo(3,3,maxtri)

  real *8 :: verts(3,100), ifaces(3,100), verts1(3,3)
  real *8 :: verts2(3,3), verts3(3,3), verts4(3,3)

  !
  ! this routine returns minimal (but possibly oversampled)
  ! triangulations of the platonic solids (i.e. a cube that isn't
  ! oversampled has 12 triangles)
  !
  !      Input:
  ! itype - determine the platonic solid to return
  !      1 - tetrahedron
  !      2 - cube
  !      3 - octahedron
  !      4 - icosahedron
  ! nover - determines how many times to oversample the skeleton
  !    triangulation (every oversampling generates 4 time as many
  !    triangles)
  ! maxtri - maximum number of allowable triangles
  !
  !      Output:
  ! npatches - number of triangles
  ! triainfo - triangle info
  ! isides - the side that the (possibly oversampled) triangle
  !     belongs to (relative to original triangulation), arbitrarily
  !     ordered.
  !
  !
  call xtri_rsolid(itype, verts, nverts, ifaces, nfaces)
  call xtri_gentriainfo(verts, nverts, ifaces, nfaces, triainfo)
  npatches = nfaces

  do i = 1,nfaces
    isides(i) = i
  end do

  if (nover .eq. 0) then
    return
  end if

  !
  ! go about oversampling the triangles
  !
  do ijk = 1,nover

    nnn = npatches
    if (nnn*4 .gt. maxtri) then
      call prinf('maxtri exceeded, maxtri = *', maxtri, 1)
      stop
    end if

    do i = 1,nnn

      call xtri_refine4_flat(triainfo(1,1,i), verts1, &
          verts2, verts3, verts4)
      do j = 1,3
        do k = 1,3
          triainfo(k,j,i) = verts1(k,j)
          triainfo(k,j,npatches+1) = verts2(k,j)
          triainfo(k,j,npatches+2) = verts3(k,j)
          triainfo(k,j,npatches+3) = verts4(k,j)
        end do
        isides(npatches+1) = isides(i)
        isides(npatches+2) = isides(i)
        isides(npatches+3) = isides(i)
      end do

      npatches = npatches +3
    end do
  end do

  return
end subroutine xtri_platonic





subroutine xtri_rsolid(itype, verts, nverts, ifaces, nfaces)
  implicit real *8 (a-h,o-z)
  dimension verts(3,*),ifaces(3,*)
  !
  ! This subroutine returns the vertices and faces of regular (and
  ! not so regular) polyhedra.
  !
  !         Input:
  ! itype - determine the platonic solid to return
  !      1 - tetrahedron
  !      2 - cube
  !      3 - octahedron
  !      4 - icosahedron
  !

  if( itype .eq. 1 ) then
    !c
    !c       ... tetrahedron
    !c
    nverts=4
    nfaces=4

    ! ... vertices
    kk=1
    verts(1,kk)=-1
    verts(2,kk)=-1/sqrt(3.0d0)
    verts(3,kk)=-1/sqrt(6.0d0)

    kk=kk+1
    verts(1,kk)=+1
    verts(2,kk)=-1/sqrt(3.0d0)
    verts(3,kk)=-1/sqrt(6.0d0)

    kk=kk+1
    verts(1,kk)=0
    verts(2,kk)=+2/sqrt(3.0d0)
    verts(3,kk)=-1/sqrt(6.0d0)

    kk=kk+1
    verts(1,kk)=0
    verts(2,kk)=0
    verts(3,kk)=+3/sqrt(6.0d0)

    !... faces
    kk=1
    ifaces(1,kk)=2
    ifaces(2,kk)=1
    ifaces(3,kk)=3

    kk=kk+1
    ifaces(1,kk)=1
    ifaces(2,kk)=2
    ifaces(3,kk)=4

    kk=kk+1
    ifaces(1,kk)=2
    ifaces(2,kk)=3
    ifaces(3,kk)=4

    kk=kk+1
    ifaces(1,kk)=3
    ifaces(2,kk)=1
    ifaces(3,kk)=4

  else if (itype .eq. 2) then

    !
    ! a cube
    !
    nverts=8
    nfaces=12

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
    ifaces(1,kk)=4
    ifaces(2,kk)=3
    ifaces(3,kk)=1

    kk=kk+1
    ifaces(1,kk)=1
    ifaces(2,kk)=2
    ifaces(3,kk)=6

    kk=kk+1
    ifaces(1,kk)=6
    ifaces(2,kk)=5
    ifaces(3,kk)=1

    kk=kk+1
    ifaces(1,kk)=2
    ifaces(2,kk)=3
    ifaces(3,kk)=7

    kk=kk+1
    ifaces(1,kk)=7
    ifaces(2,kk)=6
    ifaces(3,kk)=2

    kk=kk+1
    ifaces(1,kk)=3
    ifaces(2,kk)=4
    ifaces(3,kk)=8

    kk=kk+1
    ifaces(1,kk)=8
    ifaces(2,kk)=7
    ifaces(3,kk)=3

    kk=kk+1
    ifaces(1,kk)=4
    ifaces(2,kk)=1
    ifaces(3,kk)=5

    kk=kk+1
    ifaces(1,kk)=5
    ifaces(2,kk)=8
    ifaces(3,kk)=4

    kk=kk+1
    ifaces(1,kk)=5
    ifaces(2,kk)=6
    ifaces(3,kk)=7

    kk=kk+1
    ifaces(1,kk)=7
    ifaces(2,kk)=8
    ifaces(3,kk)=5

  else if (itype .eq. 3) then

    !
    ! an octahedron
    !
    nverts=6
    nfaces=8

    ! ... vertices
    kk=1
    verts(1,kk)=-1
    verts(2,kk)=0
    verts(3,kk)=0

    kk=kk+1
    verts(1,kk)=0
    verts(2,kk)=-1
    verts(3,kk)=0

    kk=kk+1
    verts(1,kk)=+1
    verts(2,kk)=0
    verts(3,kk)=0

    kk=kk+1
    verts(1,kk)=0
    verts(2,kk)=+1
    verts(3,kk)=0

    kk=kk+1
    verts(1,kk)=0
    verts(2,kk)=0
    verts(3,kk)=+1

    kk=kk+1
    verts(1,kk)=0
    verts(2,kk)=0
    verts(3,kk)=-1

    ! ... faces
    kk=1
    ifaces(1,kk)=1
    ifaces(2,kk)=2
    ifaces(3,kk)=5

    kk=kk+1
    ifaces(1,kk)=2
    ifaces(2,kk)=3
    ifaces(3,kk)=5

    kk=kk+1
    ifaces(1,kk)=3
    ifaces(2,kk)=4
    ifaces(3,kk)=5

    kk=kk+1
    ifaces(1,kk)=4
    ifaces(2,kk)=1
    ifaces(3,kk)=5

    kk=kk+1
    ifaces(1,kk)=2
    ifaces(2,kk)=1
    ifaces(3,kk)=6

    kk=kk+1
    ifaces(1,kk)=3
    ifaces(2,kk)=2
    ifaces(3,kk)=6

    kk=kk+1
    ifaces(1,kk)=4
    ifaces(2,kk)=3
    ifaces(3,kk)=6

    kk=kk+1
    ifaces(1,kk)=1
    ifaces(2,kk)=4
    ifaces(3,kk)=6

  else if (itype .eq. 4) then

    !
    ! an icosahedron
    !
    nverts=12
    nfaces=20

    ! . . . vertices
    !
    done=1
    pi=4*atan(done)

    do i=1,5
      verts(1,i)=cos(2*pi*i/5.0d0)
      verts(2,i)=sin(2*pi*i/5.0d0)
      verts(3,i)=0
    end do

    !
    ! find the top vertex
    !
    s=sqrt((verts(1,2)-verts(1,1))**2+(verts(2,2)-verts(2,1))**2)
    d=sqrt(s**2-1)

    verts(1,6)=0
    verts(2,6)=0
    verts(3,6)=d

    x=(1-d**2)/(2*d)
    do i=1,6
      verts(3,i)=verts(3,i)+x
    end do

    do i=1,6
      verts(1,i+6)=-verts(1,i)
      verts(2,i+6)=-verts(2,i)
      verts(3,i+6)=-verts(3,i)
    end do

    !
    ! ... faces
    kk=1
    ifaces(1,kk)=1
    ifaces(2,kk)=2
    ifaces(3,kk)=6

    kk=kk+1
    ifaces(1,kk)=2
    ifaces(2,kk)=3
    ifaces(3,kk)=6

    kk=kk+1
    ifaces(1,kk)=3
    ifaces(2,kk)=4
    ifaces(3,kk)=6

    kk=kk+1
    ifaces(1,kk)=4
    ifaces(2,kk)=5
    ifaces(3,kk)=6

    kk=kk+1
    ifaces(1,kk)=5
    ifaces(2,kk)=1
    ifaces(3,kk)=6

    kk=kk+1
    ifaces(1,kk)=8
    ifaces(2,kk)=7
    ifaces(3,kk)=12

    kk=kk+1
    ifaces(1,kk)=9
    ifaces(2,kk)=8
    ifaces(3,kk)=12

    kk=kk+1
    ifaces(1,kk)=10
    ifaces(2,kk)=9
    ifaces(3,kk)=12

    kk=kk+1
    ifaces(1,kk)=11
    ifaces(2,kk)=10
    ifaces(3,kk)=12

    kk=kk+1
    ifaces(1,kk)=7
    ifaces(2,kk)=11
    ifaces(3,kk)=12

    kk=kk+1
    ifaces(1,kk)=2
    ifaces(2,kk)=1
    ifaces(3,kk)=10

    kk=kk+1
    ifaces(1,kk)=3
    ifaces(2,kk)=2
    ifaces(3,kk)=11

    kk=kk+1
    ifaces(1,kk)=4
    ifaces(2,kk)=3
    ifaces(3,kk)=7

    kk=kk+1
    ifaces(1,kk)=5
    ifaces(2,kk)=4
    ifaces(3,kk)=8

    kk=kk+1
    ifaces(1,kk)=1
    ifaces(2,kk)=5
    ifaces(3,kk)=9

    kk=kk+1
    ifaces(1,kk)=8
    ifaces(2,kk)=9
    ifaces(3,kk)=5

    kk=kk+1
    ifaces(1,kk)=9
    ifaces(2,kk)=10
    ifaces(3,kk)=1

    kk=kk+1
    ifaces(1,kk)=10
    ifaces(2,kk)=11
    ifaces(3,kk)=2

    kk=kk+1
    ifaces(1,kk)=11
    ifaces(2,kk)=7
    ifaces(3,kk)=3

    kk=kk+1
    ifaces(1,kk)=7
    ifaces(2,kk)=8
    ifaces(3,kk)=4

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
end subroutine xtri_rsolid





subroutine xtri_gentriainfo(verts,nverts,ifaces,nfaces,triainfo)
  implicit real *8 (a-h,o-z)
  dimension verts(3,1),ifaces(3,1),triainfo(3,3,1)

  do i=1,nfaces

    triainfo(1,1,i)=verts(1,ifaces(1,i))
    triainfo(2,1,i)=verts(2,ifaces(1,i))
    triainfo(3,1,i)=verts(3,ifaces(1,i))

    triainfo(1,2,i)=verts(1,ifaces(2,i))
    triainfo(2,2,i)=verts(2,ifaces(2,i))
    triainfo(3,2,i)=verts(3,ifaces(2,i))

    triainfo(1,3,i)=verts(1,ifaces(3,i))
    triainfo(2,3,i)=verts(2,ifaces(3,i))
    triainfo(3,3,i)=verts(3,ifaces(3,i))

  end do

  return
end subroutine xtri_gentriainfo





subroutine xtri_refine4_flat(verts, verts1, verts2, verts3, verts4)
  implicit real *8 (a-h,o-z)
  real *8 :: verts(3,3), verts1(3,3), verts2(3,3), verts3(3,3)
  real *8 :: verts4(3,3)

  real *8 :: xyz12(3), xyz23(3), xyz13(3)

  !
  ! perform a refinement of a flat triangle into four other flat
  ! triangles
  !
  xyz12(1) = (verts(1,1) + verts(1,2))/2
  xyz12(2) = (verts(2,1) + verts(2,2))/2
  xyz12(3) = (verts(3,1) + verts(3,2))/2

  xyz23(1) = (verts(1,2) + verts(1,3))/2
  xyz23(2) = (verts(2,2) + verts(2,3))/2
  xyz23(3) = (verts(3,2) + verts(3,3))/2

  xyz13(1) = (verts(1,1) + verts(1,3))/2
  xyz13(2) = (verts(2,1) + verts(2,3))/2
  xyz13(3) = (verts(3,1) + verts(3,3))/2

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
  verts2(1,1) = xyz23(1)
  verts2(2,1) = xyz23(2)
  verts2(3,1) = xyz23(3)

  verts2(1,2) = xyz13(1)
  verts2(2,2) = xyz13(2)
  verts2(3,2) = xyz13(3)

  verts2(1,3) = xyz12(1)
  verts2(2,3) = xyz12(2)
  verts2(3,3) = xyz12(3)

  !
  ! third subdivision
  !
  verts3(1,1) = xyz12(1)
  verts3(2,1) = xyz12(2)
  verts3(3,1) = xyz12(3)

  verts3(1,2) = verts(1,2)
  verts3(2,2) = verts(2,2)
  verts3(3,2) = verts(3,2)

  verts3(1,3) = xyz23(1)
  verts3(2,3) = xyz23(2)
  verts3(3,3) = xyz23(3)

  !
  ! fourth subdivision
  !
  verts4(1,1) = xyz13(1)
  verts4(2,1) = xyz13(2)
  verts4(3,1) = xyz13(3)

  verts4(1,2) = xyz23(1)
  verts4(2,2) = xyz23(2)
  verts4(3,2) = xyz23(3)

  verts4(1,3) = verts(1,3)
  verts4(2,3) = verts(2,3)
  verts4(3,3) = verts(3,3)

  return
end subroutine xtri_refine4_flat



