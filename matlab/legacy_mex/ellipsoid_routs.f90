
      subroutine get_ellipsoid_mem(a,b,c,rmax,ifc,npatches)
!
!  This subroutine estimates the number of patches required
!  discretizing an ellipsoid with half axis lengths (a,b,c).
!
!  The mesh can either be conforming or non-conforming with
!  a specified maximum patch size
!
!
!  Input arguments:
!    - a: real *8
!        half axis length in the x direction
!    - b: real *8
!        half axis length in the y direction
!    - c: real *8
!        half axis length in the z direction
!    - rmax: real *8
!        maximum patch size
!    - ifc: integer
!        flag for determining if the mesh should be conforming
!        or not
!        ifc = 1, then mesh is conforming
!        ifc = 0, then mesh is non-conforming
!
!  Output arguments:
!    - npatches: integer
!        number of patches in the discretization
!
      implicit real *8 (a-h,o-z)
      real *8, intent(in) :: a,b,c,rmax
      integer, intent(in) :: ifc
      integer, intent(out) :: npatches

      done = 1.0d0
      pi = atan(done)*4

      npatches = 0

      nthet = ceiling(2*c/rmax)
      if(ifc.eq.1) then

        alpha = a
        beta = b
        hh = (alpha-beta)**2/(alpha+beta)**2 

        ellip_p = pi*(alpha + beta)* &
           (1.0d0 + 3*hh/(10.0d0 + sqrt(4.0d0-3*hh)))
        
        nphi = ceiling(ellip_p/rmax)
        npatches = 2*nthet*nphi
      else
        hthet = pi/(nthet+0.0d0)
        do ithet=1,nthet
          t0 = (ithet-1)*hthet
          t1 = (ithet)*hthet

          tuse = t0
          if(abs(t0-pi/2).ge.abs(t1-pi/2)) tuse = t1
          if((t0-pi/2)*(t1-pi/2).le.0) tuse = pi/2


          alpha = a*sin(tuse)
          beta = b*sin(tuse)
          hh = (alpha-beta)**2/(alpha+beta)**2 

          ellip_p = pi*(alpha + beta)* &
             (1.0d0 + 3*hh/(10.0d0 + sqrt(4.0d0-3*hh)))
        
          nphi = ceiling(ellip_p/rmax)
          npatches = npatches + 2*nphi
        enddo
      endif

      return
      end subroutine get_ellipsoid_mem




      subroutine get_ellipsoid_geom(a,b,c,rmax,ifc,norder,npatches, & 
        npts,norders,ixyzs,iptype,srcvals,srccoefs)
!
!  This subroutine discretizes 
!  an ellipsoid with half axis lengths (a,b,c).
!
!
!  Input arguments:
!    - a: real *8
!        half axis length in the x direction
!    - b: real *8
!        half axis length in the y direction
!    - c: real *8
!        half axis length in the z direction
!    - rmax: real *8
!        maximum patch size
!    - ifc: integer
!        flag for determining if the mesh should be conforming
!        or not
!        ifc = 1, then mesh is conforming
!        ifc = 0, then mesh is non-conforming
!    - norder:
!        order of discretization
!    - npatches: integer
!        number of patches in the discretization can be computed
!        by a call to get_ellipsoid_mem
!    - npts: integer
!        number of points in the discretization = 
!          npatches*(norder+1)*(norder+2)/2
!
!  Output arguments:
!    - norders: integer(npatches)
!        discretization order of patches
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!
      implicit real *8 (a-h,o-z)
      real *8, intent(in) :: a,b,c,rmax
      integer, intent(in) :: ifc,npatches,norder,npts
      integer, intent(out) :: norders(npatches),ixyzs(npatches+1)
      integer, intent(out) :: iptype(npatches)
      real *8, intent(out) :: srcvals(12,npts),srccoefs(9,npts)


      real *8 v1(3),v2(3),v3(3),v4(3),xyz0(1:3)
      real *8, allocatable, target :: triaskel(:,:,:)
      real *8, pointer :: ptr1,ptr2,ptr3,ptr4
      real *8, allocatable :: uvs(:,:),wts(:),umatr(:,:),vmatr(:,:)
      real *8, target :: p1(10),p2(10),p3(10),p4(10)
      external xtri_ellipsoid_eval

      npols = (norder+1)*(norder+2)/2
      do i=1,npatches
        norders(i) = norder
        ixyzs(i) = (i-1)*npols + 1
        iptype(i) = 1
      enddo
      ixyzs(npatches+1) = npts+1

      allocate(triaskel(3,3,npatches))

      done = 1.0d0
      pi = atan(done)*4


      nthet = ceiling(2*c/rmax)
      vmin = 0
      vmax = 2*pi
      nover = 0

      if(ifc.eq.1) then

        n0 = 0
        alpha = a
        beta = b
        hh = (alpha-beta)**2/(alpha+beta)**2 

        ellip_p = pi*(alpha + beta)* &
           (1.0d0 + 3*hh/(10.0d0 + sqrt(4.0d0-3*hh)))
        
        umin = 0
        umax = pi
        nphi = ceiling(ellip_p/rmax)
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,nthet,nphi,nover, &
          npatches,n0,triaskel)
      else
        hthet = pi/(nthet+0.0d0)
        istart = 1
        nthet0 = 1
        do ithet=1,nthet
          n0 = 0
          umin = (ithet-1.0d0)*hthet
          umax = (ithet+0.0d0)*hthet

          tuse = umin
          if(abs(umin-pi/2).ge.abs(umax-pi/2)) tuse = umax
          if((umin-pi/2)*(umax-pi/2).le.0) tuse = pi/2


          alpha = a*sin(tuse)
          beta = b*sin(tuse)
          hh = (alpha-beta)**2/(alpha+beta)**2 

          ellip_p = pi*(alpha + beta)* &
             (1.0d0 + 3*hh/(10.0d0 + sqrt(4.0d0-3*hh)))
        
          nphi = ceiling(ellip_p/rmax)
          call xtri_rectmesh_ani(umin,umax,vmin,vmax,nthet0,nphi,nover, &
            npatches,n0,triaskel(1,1,istart))
          istart = istart + 2*nphi
        enddo
      endif
      xyz0(1:3) = 0

      p2(1) = a
      p2(2) = b
      p2(3) = c

      p3(1:3) = xyz0(1:3)
      
      ptr1 => triaskel(1,1,1)
      ptr2 => p2(1)
      ptr3 => p3(1)
      ptr4 => p4(1)

      npols = (norder+1)*(norder+2)/2
      allocate(uvs(2,npols),wts(npols),umatr(npols,npols), &
        vmatr(npols,npols))
      call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)
      call getgeominfo(npatches,xtri_ellipsoid_eval,ptr1,ptr2,ptr3, &
        ptr4,npols,uvs,umatr,srcvals,srccoefs)


      return
      end subroutine get_ellipsoid_geom
!
!
!
!
!



      subroutine xtri_ellipsoid_eval(itri, u, v, xyz, dxyzduv, & 
          triainfo,p2, p3, p4)
      implicit real *8 (a-h,o-z)
      real *8 :: xyz(3), dxyzduv(3,2), triainfo(3,3,*),p2(3),p3(3)


      !
      ! Evaluate chart of ellipsoid in polar coordinates, u
      ! is elevation, and phi is azimuth
      !
      !    Input:
      ! itri - triangle number to map
      ! u,v - local uv coordinates on triangle itri
      ! triainfo - flat skeleton triangle info
      ! p2,p3,p4 - dummy parameters
      !
      !    Output:
      ! xyz - point on the ellipsoid
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
      ! evaluate ellipsoid chart 
      !
      xyz(1)=p2(1)*sin(x)*cos(y) + p3(1)
      xyz(2)=p2(2)*sin(x)*sin(y) + p3(2)
      xyz(3)=p2(3)*cos(x) + p3(3)


      ! du
      dxyzduv(1,1) = p2(1)*(cos(x)*cos(y)*dxdu - sin(x)*sin(y)*dydu)
      dxyzduv(2,1) = p2(2)*(cos(x)*sin(y)*dxdu + sin(x)*cos(y)*dydu)
      dxyzduv(3,1) = p2(3)*(-sin(x)*dxdu)

      ! dv
      dxyzduv(1,2) = p2(1)*(cos(x)*cos(y)*dxdv - sin(x)*sin(y)*dydv)
      dxyzduv(2,2) = p2(2)*(cos(x)*sin(y)*dxdv + sin(x)*cos(y)*dydv)
      dxyzduv(3,2) = p2(3)*(-sin(x)*dxdv)

      return
      end subroutine xtri_ellipsoid_eval



      subroutine ellipsoid_interp(nd,a,b,c,rmax,ifc,ntarg,xyztarg,npatches,&
        norders,ixyzs,iptype,npts,u,uinterp)
!
!  This function interpolates a collection of functions defined on a triangulated
!  ellipsoid with half axes a,b,c, maximum patch size rmax, at a collection
!  of user specified targets on the ellipsoid provided in their cartesian 
!  coordinates. 
!
!  Note that this subroutine assumes that the targets are on surface, 
!  if they are not, then the interpolated values will be returned
!  at the point with polar coorindates (\theta,\phi) for the target (x/a,y/b,z/c).
!
!  Input arguments:
!    - nd: integer
!        number of functions to interpolate
!    - a: real *8
!        half axis length in the x direction
!    - b: real *8
!        half axis length in the y direction
!    - c: real *8
!        half axis length in the z direction
!    - rmax: real *8
!        maximum patch size
!    - ifc: integer
!        flag for determining if the mesh should be conforming
!        or not
!        ifc = 1, then mesh is conforming
!        ifc = 0, then mesh is non-conforming
!    - ntarg: integer
!        number of targets
!    - xyztarg: real *8 (3,ntarg)
!        xyz coordinates of the targets
!    - npatches: integer
!        number of patches in the discretization can be computed
!        by a call to get_ellipsoid_mem
!    - norders: integer(npatches)
!        discretization order of patches
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer
!        number of points in the discretization = 
!          npatches*(norder+1)*(norder+2)/2
!    - u: real *8 (nd,npts)
!        function values at the boundary points
!
!  Output arguments:
!    - uinterp: real *8 (nd,ntarg)
!        interpolated function values at the target locations
!

      implicit real *8 (a-h,o-z)
      integer, intent(in) :: nd
      real *8, intent(in) :: a,b,c,rmax
      integer, intent(in) :: ifc,ntarg
      real *8, intent(in) :: xyztarg(3,ntarg)
      integer, intent(in) :: npatches,norders(npatches),ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches),npts
      real *8, intent(in) :: u(nd,npts)
      real *8, intent(out) :: uinterp(nd,ntarg)
!
!  Temporary variables
!
      integer, allocatable :: ithet_start(:),nphis(:)
      real *8, allocatable :: ucoefs(:,:),pols(:)
      real *8 xyz0(3),uvs_targ(2)

      done = 1.0d0
      pi = atan(done)*4
      
      allocate(ucoefs(nd,npts))
      call surf_vals_to_coefs(nd,npatches,norders,ixyzs,iptype,npts,&
        u,ucoefs)
      nthet = ceiling(2*c/rmax)
      hthet = pi/(nthet+0.0d0)

      nordermax = maxval(norders(1:npatches))
      npmax = (nordermax+1)*(nordermax+2)/2
      allocate(pols(npmax))
!
!  Figure out how many patches there in theta and phi directions
!
      
      allocate(ithet_start(nthet+1),nphis(nthet))

      if(ifc.eq.1) then

        alpha = a
        beta = b
        hh = (alpha-beta)**2/(alpha+beta)**2 

        ellip_p = pi*(alpha + beta)* &
           (1.0d0 + 3*hh/(10.0d0 + sqrt(4.0d0-3*hh)))
        
        nphi = ceiling(ellip_p/rmax)
        do i=1,nthet
          ithet_start(i) = (i-1)*nphi*2 + 1
          nphis(i) = nphi
        enddo
        ithet_start(nthet+1) = nthet*nphi*2+1
      else
        ithet_start(1) = 1
        do ithet=1,nthet
          t0 = (ithet-1)*hthet
          t1 = (ithet)*hthet

          tuse = t0
          if(abs(t0-pi/2).ge.abs(t1-pi/2)) tuse = t1
          if((t0-pi/2)*(t1-pi/2).le.0) tuse = pi/2


          alpha = a*sin(tuse)
          beta = b*sin(tuse)
          hh = (alpha-beta)**2/(alpha+beta)**2 

          ellip_p = pi*(alpha + beta)* &
             (1.0d0 + 3*hh/(10.0d0 + sqrt(4.0d0-3*hh)))
        
          nphi = ceiling(ellip_p/rmax)
          ithet_start(ithet+1) = ithet_start(ithet) + 2*nphi
          nphis(ithet) = nphi
        enddo
      endif


      alpha = 1.0d0
      beta = 0.0d0
      incx = 1
      incy = 1


!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,xyz0,r,theta,phi) &
!$OMP PRIVATE(ithetuse,hphi,iphiuse,u0,v0,uuse,vuse,ipatch_id) &
!$OMP PRIVATE(uvs_targ,pols,norder,npols,ii)
      do i=1,ntarg
        xyz0(1) = xyztarg(1,i)/a
        xyz0(2) = xyztarg(2,i)/b
        xyz0(3) = xyztarg(3,i)/c

        call cart2polar(xyz0,r,theta,phi)
        if(phi.lt.0) phi = phi + 2*pi

!
!  figure out which ithet and which iphi
!
!
        ithetuse = ceiling(theta/hthet)
        if(ithetuse.eq.0) ithetuse = 1

        hphi = (2*pi/nphis(ithetuse))
        iphiuse = ceiling(phi/hphi)
        if(iphiuse.eq.0) iphiuse = 1

        
        
        u0 = (ithetuse-1)*hthet
        v0 = (iphiuse-1.0d0)*hphi

        uuse = (theta - u0)/hthet
        vuse = (phi - v0)/hphi


!
!
!     find which patch and local uv coordinates on each patch
! 

        if(uuse+vuse.le.1) then
          ipatch_id = ithet_start(ithetuse) + 2*(iphiuse-1) 
          uvs_targ(1) = uuse
          uvs_targ(2) = vuse

        else
          ipatch_id = ithet_start(ithetuse) + 2*(iphiuse-1)+1
          uvs_targ(1) = 1.0d0-uuse
          uvs_targ(2) = 1.0d0-vuse
        endif



!
!  Interpolate
! 

        norder = norders(ipatch_id)
        npols = (norder+1)*(norder+2)/2
        call koorn_pols(uvs_targ,norder,npols,pols)
        ii = ixyzs(ipatch_id)

        call dgemv_guru('n',nd,npols,alpha,ucoefs(1,ii),nd,pols,incx,beta, &
          uinterp(1,i),incy)
      enddo
!$OMP END PARALLEL DO      



      end subroutine ellipsoid_interp
!
!
!
!
!

      subroutine ellipsoid_local_coord_targ(a,b,c,rmax,ifc,ntarg,xyztarg, &
        npatches,norders,ixyzs,iptype,npts,ipatchtarg,uvs_targ)

!
!  This function extracts the local patch id and uv coordinates
!  for a collection of targets on a triangulated ellipsoid
!  ellipsoid with half axes a,b,c, maximum patch size rmax, at a collection
!  of user specified targets on the ellipsoid provided in their cartesian 
!  coordinates. 
!
!  Note that this subroutine assumes that the targets are on surface, 
!  if they are not, then the interpolated values will be returned
!  at the point with polar coorindates (\theta,\phi) for the target (x/a,y/b,z/c).
!
!  Input arguments:
!    - a: real *8
!        half axis length in the x direction
!    - b: real *8
!        half axis length in the y direction
!    - c: real *8
!        half axis length in the z direction
!    - rmax: real *8
!        maximum patch size
!    - ifc: integer
!        flag for determining if the mesh should be conforming
!        or not
!        ifc = 1, then mesh is conforming
!        ifc = 0, then mesh is non-conforming
!    - ntarg: integer
!        number of targets
!    - xyztarg: real *8 (3,ntarg)
!        xyz coordinates of the targets
!    - npatches: integer
!        number of patches in the discretization can be computed
!        by a call to get_ellipsoid_mem
!    - norders: integer(npatches)
!        discretization order of patches
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer
!        number of points in the discretization = 
!          npatches*(norder+1)*(norder+2)/2
!
!  Output arguments:
!    - ipatchtarg: integer(ntarg)
!        ipatchtarg(i) is the patch on which target i is located
!    - uvs_targ: real *8 (2,ntarg)
!        local uv coordinates on patch
!

      implicit real *8 (a-h,o-z)
      real *8, intent(in) :: a,b,c,rmax
      integer, intent(in) :: ifc,ntarg
      real *8, intent(in) :: xyztarg(3,ntarg)
      integer, intent(in) :: npatches,norders(npatches),ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches),npts
      integer, intent(out) :: ipatchtarg(ntarg)
      real *8, intent(out) :: uvs_targ(2,ntarg)
!
!  Temporary variables
!
      integer, allocatable :: ithet_start(:),nphis(:)
      real *8 xyz0(3)

      done = 1.0d0
      pi = atan(done)*4
      
      nthet = ceiling(2*c/rmax)
      hthet = pi/(nthet+0.0d0)

      nordermax = maxval(norders(1:npatches))
      npmax = (nordermax+1)*(nordermax+2)/2
!
!  Figure out how many patches there in theta and phi directions
!
      
      allocate(ithet_start(nthet+1),nphis(nthet))

      if(ifc.eq.1) then

        alpha = a
        beta = b
        hh = (alpha-beta)**2/(alpha+beta)**2 

        ellip_p = pi*(alpha + beta)* &
           (1.0d0 + 3*hh/(10.0d0 + sqrt(4.0d0-3*hh)))
        
        nphi = ceiling(ellip_p/rmax)
        do i=1,nthet
          ithet_start(i) = (i-1)*nphi*2 + 1
          nphis(i) = nphi
        enddo
        ithet_start(nthet+1) = nthet*nphi*2+1
      else
        ithet_start(1) = 1
        do ithet=1,nthet
          t0 = (ithet-1)*hthet
          t1 = (ithet)*hthet

          tuse = t0
          if(abs(t0-pi/2).ge.abs(t1-pi/2)) tuse = t1
          if((t0-pi/2)*(t1-pi/2).le.0) tuse = pi/2


          alpha = a*sin(tuse)
          beta = b*sin(tuse)
          hh = (alpha-beta)**2/(alpha+beta)**2 

          ellip_p = pi*(alpha + beta)* &
             (1.0d0 + 3*hh/(10.0d0 + sqrt(4.0d0-3*hh)))
        
          nphi = ceiling(ellip_p/rmax)
          ithet_start(ithet+1) = ithet_start(ithet) + 2*nphi
          nphis(ithet) = nphi
        enddo
      endif


      alpha = 1.0d0
      beta = 0.0d0
      incx = 1
      incy = 1


!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,xyz0,r,theta,phi) &
!$OMP PRIVATE(ithetuse,hphi,iphiuse,u0,v0,uuse,vuse) 
      do i=1,ntarg
        xyz0(1) = xyztarg(1,i)/a
        xyz0(2) = xyztarg(2,i)/b
        xyz0(3) = xyztarg(3,i)/c

        call cart2polar(xyz0,r,theta,phi)
        if(phi.lt.0) phi = phi + 2*pi

!
!  figure out which ithet and which iphi
!
!
        ithetuse = ceiling(theta/hthet)
        if(ithetuse.eq.0) ithetuse = 1

        hphi = (2*pi/nphis(ithetuse))
        iphiuse = ceiling(phi/hphi)
        if(iphiuse.eq.0) iphiuse = 1

        
        
        u0 = (ithetuse-1)*hthet
        v0 = (iphiuse-1.0d0)*hphi

        uuse = (theta - u0)/hthet
        vuse = (phi - v0)/hphi


!
!
!     find which patch and local uv coordinates on each patch
! 

        if(uuse+vuse.le.1) then
          ipatchtarg(i) = ithet_start(ithetuse) + 2*(iphiuse-1) 
          uvs_targ(1,i) = uuse
          uvs_targ(2,i) = vuse

        else
          ipatchtarg(i) = ithet_start(ithetuse) + 2*(iphiuse-1)+1
          uvs_targ(1,i) = 1.0d0-uuse
          uvs_targ(2,i) = 1.0d0-vuse
        endif

      enddo
!$OMP END PARALLEL DO      



      end subroutine ellipsoid_local_coord_targ
     





      subroutine check_ellipsoid_interior(a,b,c,xyz,inflag)
!
!  This subroutine tests whether a given point xyz is inside
!  an ellipse or not 
!  
!
!  inflag = -1, inside the domain
!  inflag = 0, on the boundary
!  inflag = 1, outside the domain
!
      implicit real *8 (a-h,o-z)
      real *8, intent(in) :: a,b,c,xyz(3)
      integer, intent(out) :: inflag
      real *8 rr

      inflag = 0
      rr = xyz(1)**2/a**2 + xyz(2)**2/b**2 + xyz(3)**2/c**2
      
      if(rr.gt.1) inflag = 1
      if(rr.lt.1) inflag = -1
      

      return
      end
      
