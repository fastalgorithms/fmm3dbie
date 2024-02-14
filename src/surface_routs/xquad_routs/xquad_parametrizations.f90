!
!
!


subroutine xquad_wtorus_eval(iquad, u, v, xyz, dxyzduv, quadinfo, &
    radii, scales, p4)
  implicit real *8 (a-h,o-z)
  real *8 :: xyz(3), dxyzduv(3,2), quadinfo(3,3,*), scales(3)
  real *8 :: radii(3)

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
  s = x0+(1.0d0+u)/2*(x1-x0)+(1.0d0+v)/2*(x2-x0)
  t = y0+(1.0d0+u)/2*(y1-y0)+(1.0d0+v)/2*(y2-y0)

  rr = rmajor+rminor*cos(t)+rwave*cos(nosc*s)


  xyz(1) = a*rr*cos(s)
  xyz(2) = b*rr*sin(s)
  xyz(3) = c*rminor*sin(t)

  dsdu = (x1-x0)/2
  dsdv = (x2-x0)/2
  dtdu = (y1-y0)/2
  dtdv = (y2-y0)/2

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
end subroutine xquad_wtorus_eval
!
!
!
!
!
!
!


subroutine xquad_startorus_eval(iquad, u, v, xyz, dxyzduv, quadinfo, &
    radii, scales, p4)
  implicit real *8 (a-h,o-z)
  real *8 :: xyz(3), dxyzduv(3,2), quadinfo(3,3,*), scales(3)
  real *8 :: radii(3)

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
  s = x0+(1.0d0+u)/2*(x1-x0)+(1.0d0+v)/2*(x2-x0)
  t = y0+(1.0d0+u)/2*(y1-y0)+(1.0d0+v)/2*(y2-y0)

  rrs = rminor + rwave*cos(nosc*s)
  drrsds = -nosc*rwave*sin(nosc*s)

  rho = rmajor + rrs*cos(s)
  zz = rrs*sin(s)

  drhods = drrsds*cos(s) - rrs*sin(s)
  dzzds = drrsds*sin(s) + rrs*cos(s)

  xyz(1) = a*rho*cos(t)
  xyz(2) = b*rho*sin(t)
  xyz(3) = c*zz

  dsdu = (x1-x0)/2
  dsdv = (x2-x0)/2
  dtdu = (y1-y0)/2
  dtdv = (y2-y0)/2

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
end subroutine xquad_startorus_eval
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
