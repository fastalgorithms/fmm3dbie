


subroutine xtri_wtorus_eval(itri, u, v, xyz, dxyzduv, triainfo, &
    radii, scales, p4)
  implicit real *8 (a-h,o-z)
  real *8 :: xyz(3), dxyzduv(3,2), triainfo(3,3,*), scales(3)
  real *8 :: radii(3), dxyzdst(3,2)

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

  nosc = p4

  s = x0+u*(x1-x0)+v*(x2-x0)
  t = y0+u*(y1-y0)+v*(y2-y0)
  call wtorus_eval(s, t, radii, scales, nosc, xyz, dxyzdst)

  dsdu = (x1-x0)
  dsdv = (x2-x0)
  dtdu = (y1-y0)
  dtdv = (y2-y0)

  dxyzduv(1,1) = dxyzdst(1,1)*dsdu + dxyzdst(1,2)*dtdu
  dxyzduv(2,1) = dxyzdst(2,1)*dsdu + dxyzdst(2,2)*dtdu
  dxyzduv(3,1) = dxyzdst(3,1)*dsdu + dxyzdst(3,2)*dtdu

  dxyzduv(1,2) = dxyzdst(1,1)*dsdv + dxyzdst(1,2)*dtdv
  dxyzduv(2,2) = dxyzdst(2,1)*dsdv + dxyzdst(2,2)*dtdv
  dxyzduv(3,2) = dxyzdst(3,1)*dsdv + dxyzdst(3,2)*dtdv

  return
end subroutine xtri_wtorus_eval
!
!
!
!
!
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
  real *8 :: dxyzdst(3,2), scales(3)

  !
  ! project the triangle itri in triainfo onto a double fourier toriodal
  ! surface
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

  call xyz_tensor_fourier_eval(s, t, coefs, m, scales, xyz, dxyzdst)


  dsdu = (x1-x0)
  dsdv = (x2-x0)
  dtdu = (y1-y0)
  dtdv = (y2-y0)

  dxyzduv(1,1) = dxyzdst(1,1)*dsdu + dxyzdst(1,2)*dtdu
  dxyzduv(2,1) = dxyzdst(2,1)*dsdu + dxyzdst(2,2)*dtdu
  dxyzduv(3,1) = dxyzdst(3,1)*dsdu + dxyzdst(3,2)*dtdu

  dxyzduv(1,2) = dxyzdst(1,1)*dsdv + dxyzdst(1,2)*dtdv
  dxyzduv(2,2) = dxyzdst(2,1)*dsdv + dxyzdst(2,2)*dtdv
  dxyzduv(3,2) = dxyzdst(3,1)*dsdv + dxyzdst(3,2)*dtdv

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





subroutine xtri_trefoilknot_eval(itri, u, v, xyz, dxyzduv,triainfo, radii, scales, p4)
  implicit real *8 (a-h,o-z)
  real *8 :: xyz(3), dxyzduv(3,2), triainfo(3,3,*)
  real *8 :: radii(3),scales(4)

  !
  ! project the triangle itri in triainfo onto a trefoilknot
  ! Note: incomplete documentation
  !
  !    Input:
  ! itri - triangle number to map
  ! u,v - local uv coordinates on triangle itri
  ! triainfo - flat skeleton triangle info
  ! radii - to be explained
  ! scales - to be explained
  ! p4 - to be explained
  !
  !    Output:
  ! xyz - point on the trefoil knot
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

  r = radii(1)

  a = scales(1)
  b = scales(2)
  c = scales(3)
  d = scales(4)

  done = 1
  pi = atan(done)*4


  !
  ! ... process the geometry, return the point location on the trefoil knot
  ! and the derivatives with respect to u and v
  !
  
  s = x0+u*(x1-x0)+v*(x2-x0)
  t = y0+u*(y1-y0)+v*(y2-y0)


  xyz(1) = (a+b*dcos(c*s))*dcos(d*s) + r*dcos(c*s)*dcos(d*s) &
       *dcos(t)+(dsqrt(2.0d0)*r*(d*(a + b*dcos(c*s))*dcos(d*s) &
       *dsin(c*s)-b*c*dsin(d*s))*dsin(t))/dsqrt(2.0d0*(a**2.0d0) &
       *(d**2.0d0)+(b**2.0d0)*(2.0d0*(c**2.0d0) + d**2.0d0) &
       + b*(d**2.0d0)*(4*a*dcos(c*s) + b*dcos(2.0d0*c*s)))

  xyz(2) = (a+b*dcos(c*s))*dsin(d*s) + r*dcos(c*s)*dcos(t)*dsin(d*s) &
       + (dsqrt(2.0d0)*r*(b*c*dcos(d*s) &
       + d*(a + b*dcos(c*s))*dsin(c*s)*dsin(d*s))*dsin(t)) &
       /dsqrt(2.0d0*a**2.0d0*d**2.0d0 + b**2.0d0*(2.0d0*c**2.0d0 &
       + d**2.0d0) + b*d**2.0d0*(4*a*dcos(c*s) + b*dcos(2.0d0*c*s)))


  xyz(3) = b*dsin(c*s) + r*dcos(t)*dsin(c*s) &
       - (dsqrt(2.0d0)*d*r*dcos(c*s)*(a + b*dcos(c*s))*dsin(t)) &
       /dsqrt(2.0d0*a**2.0d0*d**2.0d0 + b**2.0d0*(2.0d0*c**2.0d0 &
       + d**2.0d0) + b*d**2.0d0*(4*a*dcos(c*s) + b*dcos(2.0d0*c*s)))

  dsdu = (x1-x0)
  dsdv = (x2-x0)
  dtdu = (y1-y0)
  dtdv = (y2-y0)

  
  dxds = -(b*c*dcos(d*s)*dsin(c*s)) - c*r*dcos(d*s)*dcos(t)*dsin(c*s) &
       - d*(a + b*dcos(c*s))*dsin(d*s)- d*r*dcos(c*s)*dcos(t)*dsin(d*s)

  dxds = dxds + (2.0d0*dsqrt(2.0d0)*b*c*d**2.0d0*r*(a + b*dcos(c*s)) &
       *dsin(c*s)*(d*(a + b*dcos(c*s))*dcos(d*s) * dsin(c*s) - b*c*dsin(d*s))*dsin(t)) &
       /(2.0d0*a**2.0d0*d**2.0d0 + b**2.0d0*(2.0d0*c**2.0d0 + d**2.0d0) + b*d**2.0d0 &
       *(4.0d0*a*dcos(c*s)+ b*dcos(2.0d0*c*s)))**1.5d0

  dxds = dxds + (dsqrt(2.0d0)*d*r*(c*dcos(d*s)*(a*dcos(c*s) - 2.0d0*b*dsin(c*s)**2.0d0) &
       - d*(a + b*dcos(c*s))*dsin(c*s)*dsin(d*s))*dsin(t)) &
       /dsqrt(2.0d0*a**2.0d0*d**2.0d0 + b**2.0d0*(2.0d0*c**2.0d0 + d**2.0d0) &
       + b*d**2.0d0*(4.0d0*a*dcos(c*s) + b*dcos(2.0d0*c*s)))

  dyds = d*(a + b*dcos(c*s))*dcos(d*s) + d*r*dcos(c*s)*dcos(d*s)*dcos(t) &
       - b*c*dsin(c*s)*dsin(d*s) - c*r*dcos(t)*dsin(c*s)*dsin(d*s)

  dyds = dyds + (2.0d0*dsqrt(2.0d0)*b*c*d**2.0d0*r*(a + b*dcos(c*s))*dsin(c*s) &
       *(b*c*dcos(d*s) + d*(a + b*dcos(c*s)) * dsin(c*s)*dsin(d*s))*dsin(t)) &
       /(2.0d0*a**2.0d0*d**2.0d0 + b**2.0d0*(2.0d0*c**2.0d0 + d**2.0d0) + b*d**2.0d0 &
       * (4.0d0*a*dcos(c*s) + b*dcos(2.0d0*c*s)))**1.5d0

  dyds = dyds + (dsqrt(2.0d0)*d*r*(d*(a + b*dcos(c*s))*dcos(d*s) * dsin(c*s) + c*(a*dcos(c*s) &
       - 2.0d0*b*dsin(c*s)**2.0d0) * dsin(d*s))*dsin(t))/dsqrt(2.0d0*a**2.0d0*d**2.0d0 &
       + b**2.0d0*(2.0d0*c**2.0d0 + d**2.0d0) + b*d**2.0d0*(4.0d0*a*dcos(c*s) + b*dcos(2.0d0*c*s))) 

  dzds = b*c*dcos(c*s) + c*r*dcos(c*s)*dcos(t)

  dzds = dzds + (dsqrt(2.0d0)*c*d*r*(a + b*dcos(c*s))* dsin(c*s)*dsin(t)) &
       /dsqrt(2.0d0*a**2.0d0*d**2.0d0 + b**2.0d0*(2.0d0*c**2.0d0 + d**2.0d0) &
       + b*d**2.0d0*(4*a*dcos(c*s) + b*dcos(2.0d0*c*s)))

  dzds = dzds - (dsqrt(2.0d0)*b*c*d**3*r * (a + b*dcos(c*s))**2.0d0*dsin(2.0d0*c*s)*dsin(t)) & 
       /(2.0d0*a**2.0d0*d**2.0d0 + b**2.0d0*(2.0d0*c**2.0d0 + d**2.0d0) + b*d**2.0d0*(4*a*dcos(c*s) &
       + b*dcos(2.0d0*c*s)))**1.5
  
  dzds = dzds + (b*c*d*r*dsin(2.0d0*c*s)*dsin(t))/(dsqrt(2.0d0)*dsqrt(2.0d0*a**2.0d0*d**2.0d0 & 
       + b**2.0d0*(2.0d0*c**2.0d0 + d**2.0d0)+ b*d**2.0d0*(4*a*dcos(c*s) + b*dcos(2.0d0*c*s))))

  dxdt = (dsqrt(2.0d0)*r*dcos(t)*(d*(a+ b*dcos(c*s))*dcos(d*s)*dsin(c*s)-b*c*dsin(d*s))) & 
       /dsqrt(2.0d0*a**2.0d0*d**2.0d0+b**2.0d0 *(2.0d0*c**2.0d0+d**2.0d0)+b*d**2.0d0 &
       *(4.0d0*a*dcos(c*s)+b*dcos(2.0d0*c*s)))- r*dcos(c*s)*dcos(d*s)*dsin(t)

  dydt = (dsqrt(2.0d0)*r*dcos(t)*(b*c*dcos(d*s) + d*(a+b*dcos(c*s))*dsin(c*s)*dsin(d*s))) &
       /dsqrt(2.0d0*a**2.0d0*d**2.0d0+b**2.0d0*(2.0d0*c**2.0d0 + d**2.0d0)+b*d**2.0d0 &
       * (4.0d0*a*dcos(c*s)+b*dcos(2.0d0*c*s)))-r*dcos(c*s)*dsin(d*s)*dsin(t)
  
  dzdt = -((dsqrt(2.0d0)*d*r*dcos(c*s)*(a+b*dcos(c*s))*dcos(t)) &
       /dsqrt(2.0d0*a**2.0d0*d**2.0d0+b**2.0d0*(2.0d0*c**2.0d0 &
       +d**2.0d0)+b*d**2.0d0*(4.0d0*a*dcos(c*s)+b*dcos(2.0d0*c*s)))) - r*dsin(c*s)*dsin(t)

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
end subroutine xtri_trefoilknot_eval

subroutine xtri_cruller_eval(itri, u, v, xyz, dxyzduv,triainfo, params, p2,p3, p4)
  implicit real *8 (a-h,o-z)
  real *8 :: xyz(3), dxyzduv(3,2), triainfo(3,3,*)
  real *8 :: params(5)

  !
  ! project the triangle itri in triainfo onto a cruller
  !
  !    Input:
  ! itri - triangle number to map
  ! u,v - local uv coordinates on triangle itri
  ! triainfo - flat skeleton triangle info
  !          params(1): major radius
  !          params(2): minor radius
  !          params(3): maximum deviation from minor radius
  !          params(4): number of pertubations
  !          params(5): number of rotations aroung the poloidal direction
  !
  ! params(4) and params(5) must be integers recast as double precision numbers
  !
  !    Output:
  ! xyz - point on the cruller
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

  r = params(1)
  a = params(2)
  b = params(3)

  c = params(4)
  d = params(5)

  done = 1
  pi = atan(done)*4

  !
  ! ... process the geometry, return the point location on the cruller
  ! and the derivatives with respect to u and v
  !
  s = x0+u*(x1-x0)+v*(x2-x0)
  t = y0+u*(y1-y0)+v*(y2-y0)


  h = a+b*dcos(c*s+d*t)
  xyz(1) =  (r+h*dcos(t))*dcos(s)
  xyz(2) =  (r+h*dcos(t))*dsin(s)
  xyz(3) = h*dsin(t)

  dsdu = (x1-x0)
  dsdv = (x2-x0)
  dtdu = (y1-y0)
  dtdv = (y2-y0)

  dxds = -((r+dcos(t)*(a+b*dcos(c*s+d*t)))*dsin(s)) &
       -b*c*dcos(s)*dcos(t)*dsin(c*s+d*t)
  dyds = dcos(s)*(r+dcos(t)*(a+b*dcos(c*s+d*t))) &
       -b*c*dcos(t)*dsin(s)*dsin(c*s+d*t)
  dzds = -b*c*dsin(t)*dsin(c*s+d*t)


  dxdt = dcos(s)*(-((a+b*dcos(c*s+d*t))*dsin(t)) &
       -b*d*dcos(t)*dsin(c*s+d*t))
  dydt = dsin(s)*(-((a+b*dcos(c*s+d*t))*dsin(t)) &
       -b*d*dcos(t)*dsin(c*s+d*t))
  dzdt = dcos(t)*(a+b*dcos(c*s+d*t))-b*d*dsin(t)*dsin(c*s+d*t)

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
end subroutine xtri_cruller_eval



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




subroutine xtri_qtree_rectmesh(umin,umax,vmin,vmax,norder,nlevels,nboxes,ltree,itree,iptr,&
     centers,boxsize,triaskel,maxtri)
  
  implicit none
  integer ibox,ilev,ntri,maxtri
  integer nlevels,nboxes,ltree,norder
  integer iptr(8)
  integer itree(ltree)
  
  real *8 centers(2,nboxes), boxsize(2,0:nlevels)
  real *8 umin,umax,vmin,vmax,width,height,boxlenx,boxleny
  real *8 triaskel(3,3,maxtri)
  real *8 u0,u1,v0,v1,c0x,c0y
!     
!     NOTE: This code assumes that points in the quad-tree
!     are mapped from [-1/2,1/2]^2.
!     
!     this routine returns a mesh consisting of triangles.
!     The leaf level rectangles in a quadtree are split into two triangels.
!     --- keep in mind the triangles are stored with u,v,w components,
!     with w=0 (in order to be consistent with other skeleton meshes)
!     
!     Input:
!     umin,umax, vmin,vmax - sets the dimensions of the rectangle
!     maxtri - maximum number of allowable triangles
!     
!     Output:
!     ntri - number of triangles created
!     triaskel - skeleton mesh info, basically just vertices of each triangle
!

  width = umax-umin
  height = vmax-vmin
  boxlenx = boxsize(1,1)
  boxleny = boxsize(2,1)

  ntri = 0
  !     loop over all boxes, if the box is a leaf, split it
  do ibox = 1,nboxes         
     if(itree(iptr(4)+ibox-1).eq.0) then
        !     is a leaf box -> split into two triangles
        ilev = itree(iptr(2)+ibox-1)
        u0 = (centers(1,ibox) - boxsize(1,ilev)/2 + 1/2) * width + umin
        u1 = (centers(1,ibox) + boxsize(1,ilev)/2 + 1/2) * width + umin
        v0 = (centers(2,ibox) - boxsize(2,ilev)/2 + 1/2) * height + vmin
        v1 = (centers(2,ibox) + boxsize(2,ilev)/2 + 1/2) * height + vmin
        call xtri_rectmesh0(u0, u1, v0, v1, triaskel(1,1,ntri+1))
        ntri= ntri + 2
     endif
  enddo

  return
end subroutine xtri_qtree_rectmesh
