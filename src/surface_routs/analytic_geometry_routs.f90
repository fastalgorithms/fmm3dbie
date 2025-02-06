!
!
!  This file contains initialization of a few analytic
!  sets of surfaces in the go3 surface format.
!
!  There will be three versions of routines:
!   * _npat: where the user specifies the number of patches
!            used to discretize the geometry
!   * _psize: where the discretization guarantees that the
!            minimal bounding spheres around patches are smaller
!            than the user prescribed value
!   * _adap: where the surface is discretized in an adaptive manner
!            to resolve the geometry to a user defined tolerance
!
!  get_ellipsoid_<>
!  get_sphere_<>
!  get_stellarator_<>
!  get_star_torus_<>
!  get_xyz_tensor_fourier_<>
!  get_xyz_tensor_fourier_sparse_<>
!
!  Currently, only the _npat routines are available.
!  All of the routines have an extra _mem routines to get
!  the number of patches in the final discretization
!  
      subroutine get_ellipsoid_npat_mem(abc, nabc, c0, norder, &
        iptype0, npatches, npts)
!
!  Estimate the number of patches in the discretization of an ellipsoid
!  given by
!
!  ((x - x0)/a)**2 + ((y-y0)/b)**2 + ((z-z0)/c)**2 = 1
!
!  where c0 = (x0, y0, z0), abc = (a,b,c)
!
!  The surface is discretized by projecting to the surface from
!  a cube -> sphere and then rescaling the xyz components.
!
!  nabc determines the number of quadrilateral patches on the faces 
!  of the cube, with na patches
!  in the x direction, nb patches in the y direction, and nz patches
!  in the z direction, with nabc = (na, nb, nc)
!
!  Input arugments:
!    - abc: real *8(3)
!        semi major axes (a,b,c) along the x,y, and z directions
!    - nabc: integer(3)
!        number of patches along the coordinate directions of the
!        cube (na, nb, nc)
!    - c0: real *8(3)
!        center of ellipsoid (x0, y0, z0)
!    - norder: integer
!        order of discretization on each patch
!    - iptype0: integer
!        type of patch to be used in the discretization
!        * iptype = 1, triangular patch discretized using RV nodes
!        * iptype = 11, quadrangular patch discretized with GL nodes
!        * iptype = 12, quadrangular patch discretized with Chebyshev 
!                       nodes
!
!  Output arguments:
!    - npatches: integer
!        number of patches
!    - npts: integer
!        Number of points on the discretized surface
!
      implicit none
      real *8, intent(in) :: abc(3), c0(3)
      integer, intent(in) :: iptype0, nabc(3)
      integer, intent(in) :: norder
      
      integer, intent(out) :: npatches, npts

      npatches = 0
      if (iptype0.eq.1) then
        npatches = 4*(nabc(1)*nabc(2) + nabc(2)*nabc(3) + &
          nabc(3)*nabc(1))
        npts = npatches*(norder+1)*(norder+2)/2
      endif

      if (iptype0.eq.11.or.iptype0.eq.12) then
        npatches = 2*(nabc(1)*nabc(2) + nabc(2)*nabc(3) + &
          nabc(3)*nabc(1))
        npts = npatches*(norder+1)*(norder+1)
      endif

      return
      end
!
!
!
!
!
      subroutine get_ellipsoid_npat(abc, nabc, c0, norder, iptype0, &
        npatches, npts, norders, ixyzs, iptype, srccoefs, srcvals)
!
!  Discretize the surface of an ellipsoid given by
!
!  ((x - x0)/a)**2 + ((y-y0)/b)**2 + ((z-z0)/c)**2 = 1
!
!  where c0 = (x0, y0, z0), abc = (a,b,c)
!
!  The surface is discretized by projecting to the surface from
!  a cube -> sphere and then rescaling the xyz components.
!
!  nabc determines the number of quadrilateral patches on the faces 
!  of the cube, with na patches
!  in the x direction, nb patches in the y direction, and nz patches
!  in the z direction, with nabc = (na, nb, nc)
!
!  Input arugments:
!    - abc: real *8(3)
!        semi major axes (a,b,c) along the x,y, and z directions
!    - nabc: integer(3)
!        number of patches along the coordinate directions of the
!        cube (na, nb, nc)
!    - c0: real *8(3)
!        center of ellipsoid (x0, y0, z0)
!    - norder: integer
!        order of discretization on each patch
!    - iptype0: integer
!        type of patch to be used in the discretization
!        * iptype = 1, triangular patch discretized using RV nodes
!        * iptype = 11, quadrangular patch discretized with GL nodes
!        * iptype = 12, quadrangular patch discretized with Chebyshev 
!                       nodes
!    - npatches: integer
!        number of patches
!    - npts: integer
!        Number of points on the discretized surface
!
!  Output arguments:
!    - norders: integer(npatches)
!        order of discretization on each patch
!        norders(i) = norder, for all i=1,2,...,npatches
!    - ixyzs: integer(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!        iptype = 11, quadrangular patch discretized with GL nodes
!        iptype = 12, quadrangular patch discretized with Chebyshev nodes
!        iptype(i) = iptype0, for all i=1,2,...,npatches
!    - srccoefs: real *8 (9,npts)
!        basis expansion coefficients of xyz, dxyz/du,
!        and dxyz/dv on each patch. 
!        For each point 
!          * srccoefs(1:3,i) is xyz info
!          * srccoefs(4:6,i) is dxyz/du info
!          * srccoefs(7:9,i) is dxyz/dv info
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the 
!        discretization nodes on the surface
!          * srcvals(1:3,i) - xyz info
!          * srcvals(4:6,i) - dxyz/du info
!          * srcvals(7:9,i) - dxyz/dv info
!          * srcvals(10:12,i) - normals info
!
      implicit none
      real *8, intent(in) :: abc(3), c0(3)
      integer, intent(in) :: iptype0, nabc(3), norder, npatches, npts

      integer, intent(out) :: norders(npatches), iptype(npatches)
      integer, intent(out) :: ixyzs(npatches+1)
      real *8, intent(out) :: srccoefs(9,npts), srcvals(12,npts)

      real *8, allocatable, target :: skel(:,:,:)
      real *8, pointer :: ptr1, ptr2, ptr3, ptr4
      real *8, allocatable :: uvs(:,:), umatr(:,:), wts(:), vmatr(:,:)
      real *8 abcuse(3)
      real *8 ause, buse, cuse
      real *8, target :: p2(3), p3(3), p4(1)
      integer na, nb, nc

      procedure (), pointer :: patchpnt
      integer i

      integer npols
      
      external xtri_ellipsoid_eval, xquad_ellipsoid_eval

      call get_npols(iptype0, norder, npols)
      allocate(uvs(2,npols), wts(npols), umatr(npols,npols))
      allocate(vmatr(npols,npols))
      call get_disc_exps(norder, npols, iptype0, uvs, umatr, vmatr, wts)
      allocate(skel(3,3,npatches))

      abcuse(1:3) = abc(1:3)/sqrt(3.0d0)
      abcuse(1:3) = 1.0d0/sqrt(3.0d0)
      ause = abcuse(1)
      buse = abcuse(2)
      cuse = abcuse(3)

      na = nabc(1)
      nb = nabc(2)
      nc = nabc(3)

      if (iptype0.eq.1) then
        call xtri_get_rectparapiped(ause, buse, cuse, na, nb, nc, &
          npatches, skel)
        patchpnt => xtri_ellipsoid_eval
      elseif (iptype0.eq.11.or.iptype0.eq.12) then
        call xquad_get_rectparapiped(ause, buse, cuse, na, nb, nc, &
          npatches, skel)
        patchpnt => xquad_ellipsoid_eval
      endif

      p2(1:3) = abc(1:3)
      p3(1:3) = c0(1:3)
      ptr1 => skel(1,1,1)
      ptr2 => p2(1)
      ptr3 => p3(1)
      ptr4 => p4(1)

      call getgeominfo(npatches, patchpnt, ptr1, ptr2, ptr3, ptr4, &
        npols, uvs, umatr, srcvals, srccoefs)

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npatches
        norders(i) = norder
        iptype(i) = iptype0
        ixyzs(i) = (i-1)*npols + 1
      enddo
!$OMP END PARALLEL DO
      
      ixyzs(npatches+1) = npts + 1

      return
      end
!
!
!
!
!
!
      subroutine get_sphere_npat_mem(a, na, c0, norder, &
        iptype0, npatches, npts)
!
!  Estimate the number of patches in the discretization of a
!  sphere of radius 'a', centered at c0
!
!  The surface is discretized by projecting to the surface from
!  a cube -> sphere 
!
!  na determines the number of quadrilateral patches on the faces 
!  of the cube
!
!  Input arugments:
!    - a: real *8
!        radius of the sphere  
!    - na: integer
!        number of patches along the coordinate directions of the cube
!    - c0: real *8(3)
!        center of ellipsoid (x0, y0, z0)
!    - norder: integer
!        order of discretization on each patch
!    - iptype0: integer
!        type of patch to be used in the discretization
!        * iptype = 1, triangular patch discretized using RV nodes
!        * iptype = 11, quadrangular patch discretized with GL nodes
!        * iptype = 12, quadrangular patch discretized with Chebyshev 
!                       nodes
!
!  Output arguments:
!    - npatches: integer
!        number of patches
!    - npts: integer
!        Number of points on the discretized surface
!
      implicit none
      real *8, intent(in) :: a, c0(3)
      integer, intent(in) :: iptype0, na
      integer, intent(in) :: norder
      
      integer, intent(out) :: npatches, npts

      integer nabc(3)
      real *8 abc(3)

      abc(1) = a
      abc(2) = a
      abc(3) = a
      
      nabc(1) = na
      nabc(2) = na
      nabc(3) = na

      call get_ellipsoid_npat_mem(abc, nabc, c0, norder, &
        iptype0, npatches, npts)

      return
      end
!
!
!
!
!
      subroutine get_sphere_npat(a, na, c0, norder, iptype0, &
        npatches, npts, norders, ixyzs, iptype, srccoefs, srcvals)
!
!  Discretize the surface of a sphere of radius 'a', centered at c0
!
!  The surface is discretized by projecting to the surface from
!  a cube -> sphere 
!
!  na determines the number of quadrilateral patches on the faces 
!  of the cube
!
!  Input arugments:
!    - a: real *8
!        radius of the sphere  
!    - na: integer
!        number of patches along the coordinate directions of the cube
!    - c0: real *8(3)
!        center of ellipsoid (x0, y0, z0)
!    - norder: integer
!        order of discretization on each patch
!    - iptype0: integer
!        type of patch to be used in the discretization
!        * iptype = 1, triangular patch discretized using RV nodes
!        * iptype = 11, quadrangular patch discretized with GL nodes
!        * iptype = 12, quadrangular patch discretized with Chebyshev 
!                       nodes
!    - npatches: integer
!        number of patches
!    - npts: integer
!        Number of points on the discretized surface
!
!  Output arguments:
!    - norders: integer(npatches)
!        order of discretization on each patch
!        norders(i) = norder, for all i=1,2,...,npatches
!    - ixyzs: integer(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!        iptype = 11, quadrangular patch discretized with GL nodes
!        iptype = 12, quadrangular patch discretized with Chebyshev nodes
!        iptype(i) = iptype0, for all i=1,2,...,npatches
!    - srccoefs: real *8 (9,npts)
!        basis expansion coefficients of xyz, dxyz/du,
!        and dxyz/dv on each patch. 
!        For each point 
!          * srccoefs(1:3,i) is xyz info
!          * srccoefs(4:6,i) is dxyz/du info
!          * srccoefs(7:9,i) is dxyz/dv info
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the 
!        discretization nodes on the surface
!          * srcvals(1:3,i) - xyz info
!          * srcvals(4:6,i) - dxyz/du info
!          * srcvals(7:9,i) - dxyz/dv info
!          * srcvals(10:12,i) - normals info
!
      implicit none
      real *8, intent(in) :: a, c0(3)
      integer, intent(in) :: iptype0, na, norder, npatches, npts

      integer, intent(out) :: norders(npatches), iptype(npatches)
      integer, intent(out) :: ixyzs(npatches+1)
      real *8, intent(out) :: srccoefs(9,npts), srcvals(12,npts)
      
      integer nabc(3)
      real *8 abc(3)

      nabc(1) = na
      nabc(2) = na
      nabc(3) = na

      abc(1) = a
      abc(2) = a
      abc(3) = a

      call get_ellipsoid_npat(abc, nabc, c0, norder, iptype0, &
        npatches, npts, norders, ixyzs, iptype, srccoefs, srcvals)

      return
      end
!
!
!
!
!
!
      subroutine get_xyz_tensor_fourier_npat_mem(coefs, m, scales, &
        iort, nuv, norder, iptype0, npatches, npts)
!
!  Estimate the number of patches in the discretization of toroidal
!  double fourier surface given by
!
!  hat(x) = (\sum_{i=1}^{2m+1} \sum_{j=0}^{2m+1} x_{ij} b_{i}(u) b_{j}(v))
!  hat(y) = (\sum_{i=1}^{2m+1} \sum_{j=0}^{2m+1} y_{ij} b_{i}(u) b_{j}(v))
!  hat(z) = (\sum_{i=1}^{2m+1} \sum_{j=0}^{2m+1} z_{ij} b_{i}(u) b_{j}(v))
!
!  x(u,v) = (\hat(x) \cos(v) - \hat(y) \sin(v))*scales(1)
!  y(u,v) = (\hat(x) \sin(v) + \hat(y) \cos(v))*scales(2)
!  z(u,v) = \hat(z)*scales(3)
!  
!  u,v \in [0, 2\pi]^2, and 
!
!  b_{i}(u) = \{1, \cos{u}, \ldots \cos{m \cdot u}\, ...
!    \sin{u}, \ldots \sin{m \cdot u} \}.
!
!  
!  Input arugments:
!    - coefs: real *8(2*m+1, 2*m+1, 3)
!        \hat(x) = coefs(2*m+1,2*m+1, 1)
!        \hat(y) = coefs(2*m+1,2*m+1, 2)
!        \hat(z) = coefs(2*m+1,2*m+1, 3)
!    - m: integer
!        max fourier content of \hat(x), \hat(y), and \hat(z)
!    - scales: real *8(3)
!        scaling parameter for the coordinates of the surface
!    - iort: integer
!        orientiation flag
!        if iort = 1, then parameter space is [0,2\pi)^2
!        if iort = -1, then parameter space is [2\pi,0) \times [0,2\pi)
!    - nuv: integer(2)
!        number of quadrilateral patches in the u, and v direction
!        respectively, if iptype is 1, then each quadrilateral
!        patch is subdivided into two triangular patches
!    - norder: integer
!        order of discretization on each patch
!    - iptype0: integer
!        type of patch to be used in the discretization
!        * iptype = 1, triangular patch discretized using RV nodes
!        * iptype = 11, quadrangular patch discretized with GL nodes
!        * iptype = 12, quadrangular patch discretized with Chebyshev 
!                       nodes
!
!  Output arguments:
!    - npatches: integer
!        number of patches = 2*nuv(1)*nuv(2) if iptype0 = 1
!                          = nuv(1)*nuv(2) if iptype0 = 11, or 12
!    - npts: integer
!        Number of points on the discretized surface
!
!
!
      implicit none
      integer, intent(in) :: m, nuv(2), iptype0, norder, iort
      real *8, intent(in) :: coefs(2*m+1,2*m+1,3), scales(3)
      integer, intent(out) :: npatches, npts
      
      npatches = 0
      if (iptype0.eq.1) then
        npatches = 2*nuv(1)*nuv(2)
        npts = npatches*(norder+1)*(norder+2)/2
      endif

      if (iptype0.eq.11.or.iptype0.eq.12) then
        npatches = nuv(1)*nuv(2) 
        npts = npatches*(norder+1)*(norder+1)
      endif

      return
      end
!      
!
!
!
!
      subroutine get_xyz_tensor_fourier_npat(coefs, m, scales, &
        iort, nuv, norder, iptype0, npatches, npts, norders, ixyzs, &
        iptype, srccoefs, srcvals)
!
!  Discretize the toroidal double fourier surface given by
!
!  hat(x) = (\sum_{i=1}^{2m+1} \sum_{j=0}^{2m+1} x_{ij} b_{i}(u) b_{j}(v))
!  hat(y) = (\sum_{i=1}^{2m+1} \sum_{j=0}^{2m+1} y_{ij} b_{i}(u) b_{j}(v))
!  hat(z) = (\sum_{i=1}^{2m+1} \sum_{j=0}^{2m+1} z_{ij} b_{i}(u) b_{j}(v))
!
!  x(u,v) = (\hat(x) \cos(v) - \hat(y) \sin(v))*scales(1)
!  y(u,v) = (\hat(x) \sin(v) + \hat(y) \cos(v))*scales(2)
!  z(u,v) = \hat(z)*scales(3)
!  
!  u,v \in [0, 2\pi]^2, and 
!
!  b_{i}(u) = \{1, \cos{u}, \ldots \cos{m \cdot u}\, ...
!    \sin{u}, \ldots \sin{m \cdot u} \}.
!
!  
!  Input arugments:
!    - coefs: real *8(2*m+1, 2*m+1, 3)
!        \hat(x) = coefs(2*m+1,2*m+1, 1)
!        \hat(y) = coefs(2*m+1,2*m+1, 2)
!        \hat(z) = coefs(2*m+1,2*m+1, 3)
!    - m: integer
!        max fourier content of \hat(x), \hat(y), and \hat(z)
!    - scales: real *8(3)
!        scaling parameter for the coordinates of the surface
!    - iort: integer
!        orientiation flag
!        if iort = 1, then parameter space is [0,2\pi)^2
!        if iort = -1, then parameter space is [2\pi,0) \times [0,2\pi)
!    - nuv: integer(2)
!        number of quadrilateral patches in the u, and v direction
!        respectively, if iptype is 1, then each quadrilateral
!        patch is subdivided into two triangular patches
!    - norder: integer
!        order of discretization on each patch
!    - iptype0: integer
!        type of patch to be used in the discretization
!        * iptype = 1, triangular patch discretized using RV nodes
!        * iptype = 11, quadrangular patch discretized with GL nodes
!        * iptype = 12, quadrangular patch discretized with Chebyshev 
!                       nodes
!    - npatches: integer
!        number of patches = 2*nuv(1)*nuv(2) if iptype0 = 1
!                          = nuv(1)*nuv(2) if iptype0 = 11, or 12
!    - npts: integer
!        Number of points on the discretized surface
!
!  Output arguments:
!    - norders: integer(npatches)
!        order of discretization on each patch
!        norders(i) = norder, for all i=1,2,...,npatches
!    - ixyzs: integer(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!        iptype = 11, quadrangular patch discretized with GL nodes
!        iptype = 12, quadrangular patch discretized with Chebyshev nodes
!        iptype(i) = iptype0, for all i=1,2,...,npatches
!    - srccoefs: real *8 (9,npts)
!        basis expansion coefficients of xyz, dxyz/du,
!        and dxyz/dv on each patch. 
!        For each point 
!          * srccoefs(1:3,i) is xyz info
!          * srccoefs(4:6,i) is dxyz/du info
!          * srccoefs(7:9,i) is dxyz/dv info
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the 
!        discretization nodes on the surface
!          * srcvals(1:3,i) - xyz info
!          * srcvals(4:6,i) - dxyz/du info
!          * srcvals(7:9,i) - dxyz/dv info
!          * srcvals(10:12,i) - normals info
!
!
      implicit none
      integer, intent(in) :: m, nuv(2), iptype0, norder
      real *8, intent(in) :: coefs(2*m+1,2*m+1,3), scales(3)
      integer, intent(in) :: npatches, npts, iort
      
      integer, intent(out) :: norders(npatches), iptype(npatches)
      integer, intent(out) :: ixyzs(npatches+1)
      real *8, intent(out) :: srccoefs(9,npts), srcvals(12,npts)
      
      integer npols
      real *8, allocatable :: uvs(:,:), wts(:), umatr(:,:), vmatr(:,:)
      real *8, allocatable, target :: skel(:,:,:)
      real *8, allocatable, target :: pcoefs(:,:,:)
      real *8, target :: p3(1), p4(3)
      real *8, pointer :: ptr1, ptr2, ptr3, ptr4
      integer, pointer :: iptr3
      integer, target :: muse


      real *8 done, pi
      procedure (), pointer :: patchpnt
      integer i, j, k
      real *8 umin, umax, vmin, vmax
      integer nover, npatches0

      external xtri_xyz_tensor_fourier_eval
      external xquad_xyz_tensor_fourier_eval


      done = 1.0d0
      pi = atan(done)*4
      
      
      call get_npols(iptype0, norder, npols)
      allocate(uvs(2,npols), wts(npols), umatr(npols,npols))
      allocate(vmatr(npols,npols))
      call get_disc_exps(norder, npols, iptype0, uvs, umatr, vmatr, wts)
      allocate(skel(3,3,npatches))


      umin = 0 
      umax = 2*pi 

      vmin = 0
      vmax = 2*pi

      if (iort.lt.0) then
        umin = 2*pi
        umax = 0 

      endif

      
      nover = 0
      npatches0 = 0
      


      if (iptype0.eq.1) then
        call xtri_rectmesh_ani(umin, umax, vmin, vmax, nuv(1), nuv(2), &
            nover, npatches, npatches0, skel)
        patchpnt => xtri_xyz_tensor_fourier_eval
      elseif (iptype0.eq.11.or.iptype0.eq.12) then
        call xquad_rectmesh_ani(umin, umax, vmin, vmax, nuv(1), nuv(2), &
            nover, npatches, npatches0, skel)
        patchpnt => xquad_xyz_tensor_fourier_eval
      endif




      allocate(pcoefs(2*m+1,2*m+1,3))

      do i=1,3
        do j=1,2*m+1
          do k=1,2*m+1
            pcoefs(k,j,i) = coefs(k,j,i)
          enddo
        enddo
      enddo

      muse = m
      p4(1:3) = scales(1:3)

      ptr1 => skel(1,1,1)
      ptr2 => pcoefs(1,1,1)
      iptr3 => muse 
      ptr4 => p4(1)
      
      call getgeominfo(npatches, patchpnt, ptr1, ptr2, iptr3, ptr4, &
        npols, uvs, umatr, srcvals, srccoefs)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,npatches
        norders(i) = norder
        iptype(i) = iptype0
        ixyzs(i) = (i-1)*npols + 1
      enddo
!$OMP END PARALLEL DO
      
      ixyzs(npatches+1) = npts + 1
      

      return
      end
!
!
!
!
!
      subroutine get_startorus_npat_mem(radii, nosc, scales, &
        nuv, norder, iptype0, npatches, npts)
!
!  Estimate the number of patches in the discretization of toroidal
!  double fourier surface given by
!
!  x(u,v) = \rho(u) \cos(v)*scales(1)
!  y(u,v) = \rho(u) \sin(v)*scales(2)
!  z(u,v) = z(u)*scales(3)
!  
!  u,v \in [0, 2\pi]^2, and
!
!  \rho(u) = Rmajor + r(u)*cos(u)
!  z(u) = r(u)*sin(u)
!
!  r(u) = rminor + rwave*cos(nosc*u)
!
!  Input arugments:
!    - radii: real *8(3)
!         radii of the star shaped torus
!         * radii(1) = rmajor
!         * radii(2) = rminor
!         * radii(3) = rwave
!    - nosc: integer
!         number of oscillations in the star-shaped torus
!    - scales: real *8(3)
!        scaling parameter for the coordinates of the surface
!    - nuv: integer(2)
!        number of quadrilateral patches in the u, and v direction
!        respectively, if iptype is 1, then each quadrilateral
!        patch is subdivided into two triangular patches
!    - norder: integer
!        order of discretization on each patch
!    - iptype0: integer
!        type of patch to be used in the discretization
!        * iptype = 1, triangular patch discretized using RV nodes
!        * iptype = 11, quadrangular patch discretized with GL nodes
!        * iptype = 12, quadrangular patch discretized with Chebyshev 
!                       nodes
!
!  Output arguments:
!    - npatches: integer
!        number of patches = 2*nuv(1)*nuv(2) if iptype0 = 1
!                          = nuv(1)*nuv(2) if iptype0 = 11, or 12
!    - npts: integer
!        Number of points on the discretized surface
!
!
!
      implicit none
      integer, intent(in) :: nosc, nuv(2), iptype0, norder
      real *8, intent(in) :: radii(3), scales(3)
      integer, intent(out) :: npatches, npts
      
      npatches = 0
      if (iptype0.eq.1) then
        npatches = 2*nuv(1)*nuv(2)
        npts = npatches*(norder+1)*(norder+2)/2
      endif

      if (iptype0.eq.11.or.iptype0.eq.12) then
        npatches = nuv(1)*nuv(2) 
        npts = npatches*(norder+1)*(norder+1)
      endif

      return
      end
!      
!
!
!
!
      subroutine get_startorus_npat(radii, nosc, scales, &
        nuv, norder, iptype0, npatches, npts, norders, ixyzs, iptype, &
        srccoefs, srcvals)
!
!  Discretize the star-shaped torus 
!
!  x(u,v) = \rho(u) \cos(v)*scales(1)
!  y(u,v) = \rho(u) \sin(v)*scales(2)
!  z(u,v) = z(u)*scales(3)
!
!  u,v \in [0, 2\pi]^2, and 
!
!  \rho(u) = Rmajor + r(u)*cos(u)
!  z(u) = r(u)*sin(u)
!
!  r(u) = rminor + rwave*cos(nosc*u)
!  
!  Input arugments:
!    - radii: real *8(3)
!         radii of the star shaped torus
!         * radii(1) = rmajor
!         * radii(2) = rminor
!         * radii(3) = rwave
!    - nosc: integer
!         number of oscillations in the star-shaped torus
!    - scales: real *8(3)
!        scaling parameter for the coordinates of the surface
!    - nuv: integer(2)
!        number of quadrilateral patches in the u, and v direction
!        respectively, if iptype is 1, then each quadrilateral
!        patch is subdivided into two triangular patches
!    - norder: integer
!        order of discretization on each patch
!    - iptype0: integer
!        type of patch to be used in the discretization
!        * iptype = 1, triangular patch discretized using RV nodes
!        * iptype = 11, quadrangular patch discretized with GL nodes
!        * iptype = 12, quadrangular patch discretized with Chebyshev 
!                       nodes
!    - npatches: integer
!        number of patches = 2*nuv(1)*nuv(2) if iptype0 = 1
!                          = nuv(1)*nuv(2) if iptype0 = 11, or 12
!    - npts: integer
!        Number of points on the discretized surface
!
!  Output arguments:
!    - norders: integer(npatches)
!        order of discretization on each patch
!        norders(i) = norder, for all i=1,2,...,npatches
!    - ixyzs: integer(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!        iptype = 11, quadrangular patch discretized with GL nodes
!        iptype = 12, quadrangular patch discretized with Chebyshev nodes
!        iptype(i) = iptype0, for all i=1,2,...,npatches
!    - srccoefs: real *8 (9,npts)
!        basis expansion coefficients of xyz, dxyz/du,
!        and dxyz/dv on each patch. 
!        For each point 
!          * srccoefs(1:3,i) is xyz info
!          * srccoefs(4:6,i) is dxyz/du info
!          * srccoefs(7:9,i) is dxyz/dv info
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the 
!        discretization nodes on the surface
!          * srcvals(1:3,i) - xyz info
!          * srcvals(4:6,i) - dxyz/du info
!          * srcvals(7:9,i) - dxyz/dv info
!          * srcvals(10:12,i) - normals info
!
!
      implicit none
      integer, intent(in) :: nosc, nuv(2), iptype0, norder
      real *8, intent(in) :: radii(3), scales(3)
      integer, intent(in) :: npatches, npts
      
      integer, intent(out) :: norders(npatches), iptype(npatches)
      integer, intent(out) :: ixyzs(npatches+1)
      real *8, intent(out) :: srccoefs(9,npts), srcvals(12,npts)

      real *8, allocatable :: coefs(:,:,:)
      integer i, j, k, m, iort

      m = nosc + 1

      allocate(coefs(2*m+1,2*m+1,3))

      do i=1,3
        do j=1,2*m+1
          do k=1,2*m+1
            coefs(k,j,i) = 0
          enddo
        enddo
      enddo

      coefs(1,1,1) = radii(1)
      coefs(2,1,1) = radii(2)
      coefs(1+nosc+1,1,1) = coefs(1+nosc+1,1,1) + radii(3)/2
      coefs(1+abs(nosc-1),1,1) = coefs(1+abs(nosc-1),1,1) + radii(3)/2
      
      coefs(m+2,1,3) = radii(2)
      coefs(2*m+1,1,3) = coefs(2*m+1,1,3) + radii(3)/2
      if (nosc-1.gt.0) coefs(m+1+nosc-1,1,3) = coefs(m+1+nosc-1,1,3)-radii(3)/2
      if (nosc-1.lt.0) coefs(m+1+abs(nosc-1),1,3) = coefs(m+1+abs(nosc-1),1,3) + radii(3)/2

      iort = -1

      call get_xyz_tensor_fourier_npat(coefs, m, scales, &
        iort, nuv, norder, iptype0, npatches, npts, norders, ixyzs, &
        iptype, srccoefs, srcvals)

      return
      end
!
!
!
!
!
      subroutine get_stellarator_npat_mem(nuv, norder, iptype0, &
        npatches, npts)
!
!  Estimate the number of patches in the discretization of 
!  the stellarator example in the fmm3dbie paper
!  given by
!
!  x(u,v) = \hat(x)(u,v) \cos(v)
!  y(u,v) = \hat(x)(u,v) \sin(v)
!  z(u,v) = \hat(z)(u,v)
!  
!  u,v \in [0, 2\pi]^2, and
!
!  hat(x) = (\sum_{i=1}^{2m+1} \sum_{j=0}^{2m+1} x_{ij} b_{i}(u) b_{j}(v))
!  hat(z) = (\sum_{i=1}^{2m+1} \sum_{j=0}^{2m+1} z_{ij} b_{i}(u) b_{j}(v))
!
!  with m = 2
!  and the only non-zero coefs given by
!  x_{3,2} = 0.17 
!  x_{5,4} = 0.17 
!  x_{1,1} = 4.5 
!  x_{2,1} = 0.75 
!  x_{3,1} = 0.11 
!  x_{2,2} = 0.07 + (-0.45) 
!  x_{4,4} = -0.07 + (-0.45) 
!
!  z_{5,2} = 0.17 
!  z_{3,4} = -0.17 
!  z_{4,1} = 1.25 
!  z_{5,1} = 0.11 
!  z_{4,2} = 0.07 - (-0.45) 
!  z_{2,4} = 0.07 + (-0.45) 
!
!  Input arugments:
!    - nuv: integer(2)
!        number of quadrilateral patches in the u, and v direction
!        respectively, if iptype is 1, then each quadrilateral
!        patch is subdivided into two triangular patches
!    - norder: integer
!        order of discretization on each patch
!    - iptype0: integer
!        type of patch to be used in the discretization
!        * iptype = 1, triangular patch discretized using RV nodes
!        * iptype = 11, quadrangular patch discretized with GL nodes
!        * iptype = 12, quadrangular patch discretized with Chebyshev 
!                       nodes
!
!  Output arguments:
!    - npatches: integer
!        number of patches = 2*nuv(1)*nuv(2) if iptype0 = 1
!                          = nuv(1)*nuv(2) if iptype0 = 11, or 12
!    - npts: integer
!        Number of points on the discretized surface
!
!
!
      implicit none
      integer, intent(in) :: nuv(2), iptype0, norder
      integer, intent(out) :: npatches, npts
      
      npatches = 0
      if (iptype0.eq.1) then
        npatches = 2*nuv(1)*nuv(2)
        npts = npatches*(norder+1)*(norder+2)/2
      endif

      if (iptype0.eq.11.or.iptype0.eq.12) then
        npatches = nuv(1)*nuv(2) 
        npts = npatches*(norder+1)*(norder+1)
      endif

      return
      end
!      
!
!
!
!
      subroutine get_stellarator_npat(nuv, norder, iptype0, npatches, &
        npts, norders, ixyzs, iptype, srccoefs, srcvals)
!
!  Discretize the stellarator example in the fmm3dbie paper
!  given by
!
!  x(u,v) = \hat(x)(u,v) \cos(v)
!  y(u,v) = \hat(x)(u,v) \sin(v)
!  z(u,v) = \hat(z)(u,v)
!  
!  u,v \in [0, 2\pi]^2, and
!
!  hat(x) = (\sum_{i=1}^{2m+1} \sum_{j=0}^{2m+1} x_{ij} b_{i}(u) b_{j}(v))
!  hat(z) = (\sum_{i=1}^{2m+1} \sum_{j=0}^{2m+1} z_{ij} b_{i}(u) b_{j}(v))
!
!  with m = 2
!  and the only non-zero coefs given by
!  x_{3,2} = 0.17 
!  x_{5,4} = 0.17 
!  x_{1,1} = 4.5 
!  x_{2,1} = 0.75 
!  x_{3,1} = 0.11 
!  x_{2,2} = 0.07 + (-0.45) 
!  x_{4,4} = -0.07 + (-0.45) 
!
!  z_{5,2} = 0.17 
!  z_{3,4} = -0.17 
!  z_{4,1} = 1.25 
!  z_{5,1} = 0.11 
!  z_{4,2} = 0.07 - (-0.45) 
!  z_{2,4} = 0.07 + (-0.45) 
!
!  Input arugments:
!    - nuv: integer(2)
!        number of quadrilateral patches in the u, and v direction
!        respectively, if iptype is 1, then each quadrilateral
!        patch is subdivided into two triangular patches
!    - norder: integer
!        order of discretization on each patch
!    - iptype0: integer
!        type of patch to be used in the discretization
!        * iptype = 1, triangular patch discretized using RV nodes
!        * iptype = 11, quadrangular patch discretized with GL nodes
!        * iptype = 12, quadrangular patch discretized with Chebyshev 
!                       nodes
!    - npatches: integer
!        number of patches = 2*nuv(1)*nuv(2) if iptype0 = 1
!                          = nuv(1)*nuv(2) if iptype0 = 11, or 12
!    - npts: integer
!        Number of points on the discretized surface
!
!  Output arguments:
!    - norders: integer(npatches)
!        order of discretization on each patch
!        norders(i) = norder, for all i=1,2,...,npatches
!    - ixyzs: integer(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!        iptype = 11, quadrangular patch discretized with GL nodes
!        iptype = 12, quadrangular patch discretized with Chebyshev nodes
!        iptype(i) = iptype0, for all i=1,2,...,npatches
!    - srccoefs: real *8 (9,npts)
!        basis expansion coefficients of xyz, dxyz/du,
!        and dxyz/dv on each patch. 
!        For each point 
!          * srccoefs(1:3,i) is xyz info
!          * srccoefs(4:6,i) is dxyz/du info
!          * srccoefs(7:9,i) is dxyz/dv info
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the 
!        discretization nodes on the surface
!          * srcvals(1:3,i) - xyz info
!          * srcvals(4:6,i) - dxyz/du info
!          * srcvals(7:9,i) - dxyz/dv info
!          * srcvals(10:12,i) - normals info
!
!
      implicit none
      integer, intent(in) :: nuv(2), iptype0, norder
      integer, intent(in) :: npatches, npts
      
      integer, intent(out) :: norders(npatches), iptype(npatches)
      integer, intent(out) :: ixyzs(npatches+1)
      real *8, intent(out) :: srccoefs(9,npts), srcvals(12,npts)

      real *8, allocatable :: coefs(:,:,:)
      real *8 scales(3)
      integer i, j, k, m, iort

      m = 2

      allocate(coefs(2*m+1,2*m+1,3))

      do i=1,3
        do j=1,2*m+1
          do k=1,2*m+1
            coefs(k,j,i) = 0
          enddo
        enddo
      enddo

      coefs(3,2,1) = 0.17d0
      coefs(5,4,1) = 0.17d0
      coefs(1,1,1) = 4.5d0
      coefs(2,1,1) = 0.75d0
      coefs(3,1,1) = 0.11d0
      coefs(2,2,1) = -0.38d0
      coefs(4,4,1) = -0.52d0

      coefs(5,2,3) = 0.17d0
      coefs(3,4,3) = -0.17d0
      coefs(4,1,3) = 1.25d0
      coefs(5,1,3) = 0.11d0
      coefs(4,2,3) = 0.52d0 
      coefs(2,4,3) = -0.38d0

      scales(1) = 1
      scales(2) = 1
      scales(3) = 1

      iort = -1


      call get_xyz_tensor_fourier_npat(coefs, m, scales, &
        iort, nuv, norder, iptype0, npatches, npts, norders, ixyzs, &
        iptype, srccoefs, srcvals)

      return
      end
!
!
!
!
      subroutine get_axissym_fcurve_npat_mem(nch2d, tchse, k, fcurve, &
        np, pars, nrts, rmid, nmid, iort, norder, iptype0, &
        npatches, npts)     
!
!  This subroutine estimates the number of patches required to 
!  discretize an axissymmetric geometry where the generating curve is 
!  provided as a function handle.  
!  The axis of rotation is assumed to be the z axis.
!
!  The geometry is given by
!  (r(s) \cos(t), r(s) \sin(t), z(s)), s \in [a,b], 
!                                      t \in [0,2\pi]
!
!  where r(s), z(s) are accesible via the function handle fcurve
!  In three dimensions, this cylindrical domain for 
!  s\in [tchse2d(i), tchse2d(i+1)] is discretized using
!  nrts(1,i) patches in the s direction, and nrts(2,i) patches in the t 
!  direction.
!
!  The above is true for all chunks in the 2d discretization expect 
!  the first and the last chunk, which are assumed to be the polar caps.
!  A circle mesh is used to discretize the polar caps. See documentation
!  in circle_mesh.f90, but briefly the polar map is a map from the
!  circle to the domain, and the circle will be discretized with nrts(1,1),
!  points in a annular region outside a square of size rmid, and the 
!  square is discretized using nmid \times nmid patches
!  
!  Input arguments:
!    - nch2d: integer
!        number of chunks describing the generating curve
!    - tchse: real *8 (nch2d+1)
!        starting and ending location of in parameter space
!        for each 2d chunk
!    - k: integer
!        number of points per chunk (for the 2d curve)
!    - fcurve: function handle
!        function handle for corresponding to generating curve.
!        Should have calling sequence
!        fcurve(s, np, pars, r, z, drds, dzds, d2rds2, d2zds2)
!    - np: integer
!        number of parameters in fcurve
!    - pars: real *8 (np)
!        parameters of fcurve
!    - nrts: interger(2,nch2d)
!        number of quadrilateral patches in the rz (s) direction,
!        and the periodic (t) direction.
!        if iptype is 1, then each quadrilateral patch
!        is subdivided into two triangular patches 
!    - rmid: real *8
!        radius of the inner square on the unit circle mesh
!        at the polar caps, should be between 0 and 1/sqrt(2)
!    - nmid: integer
!        * number of intervals in the middle square on the
!          unit circl mesh at the polar caps
!    - iort: integer
!        orientiation flag
!        if iort = 1, then parameter space is [a,b] \times [0,2\pi)
!        if iort = -1, then parameter space is [a,b] \times [2\pi,0)
!    - norder: integer
!        order of discretization on each patch
!    - iptype0: integer
!        type of patch to be used in the discretization
!        * iptype = 1, triangular patch discretized using RV nodes
!        * iptype = 11, quadrangular patch discretized with GL nodes
!        * iptype = 12, quadrangular patch discretized with Chebyshev 
!
!  Output arguments:
!    - npatches: integer
!        number of patches
!    - npts: integer
!        Number of points on the discretized surface
!
      implicit none
      integer, intent(in) :: nch2d, k
      real *8, intent(in) :: tchse(nch2d+1)
      integer, intent(in) :: nrts(2,nch2d), nmid
      integer, intent(in) :: iptype0, norder
      integer, intent(in) :: np
      real *8, intent(in) :: pars(np)
      real *8, intent(in) :: rmid
      integer, intent(in) :: iort
      integer, intent(out) :: npatches, npts

      integer i, npat0

      external fcurve
      
      npat0 = 2*nmid*nmid
      do i=2,nch2d-1
        npat0 = npat0 + nrts(1,i)*nrts(2,i)
      enddo
      npat0 = npat0 + 4*nrts(1,1)*nrts(2,1)
      npat0 = npat0 + 4*nrts(1,nch2d)*nrts(2,nch2d)

      npatches = 0
      npts = 0
      if (iptype0.eq.1) then
        npatches = 2*npat0 
        npts = npatches*(norder+1)*(norder+2)/2
      endif

      if (iptype0.eq.11.or.iptype0.eq.12) then
        npatches = npat0 
        npts = npatches*(norder+1)*(norder+1)
      endif

      

      return
      end
!
!
!
!
!
      subroutine get_axissym_fcurve_npat(nch2d, tchse, k, fcurve, &
        np, pars, nrts, rmid, nmid, iort, norder, iptype0, &
        npatches, npts, norders, ixyzs, iptype, srccoefs, srcvals)     
!
!
!  This subroutine discretizes an axissymmetric geometry where the 
!  generating curve is provided as a function handle.  
!  The axis of rotation is assumed to be the z axis.
!
!  The geometry is given by
!  (r(s) \cos(t), r(s) \sin(t), z(s)), s \in [a,b], 
!                                      t \in [0,2\pi]
!
!  where r(s), z(s) are accesible via the function handle fcurve
!  In three dimensions, this cylindrical domain for 
!  s\in [tchse2d(i), tchse2d(i+1)] is discretized using
!  nrts(1,i) patches in the s direction, and nrts(2,i) patches in the t 
!  direction.
!
!  The above is true for all chunks in the 2d discretization expect 
!  the first and the last chunk, which are assumed to be the polar caps.
!  A circle mesh is used to discretize the polar caps. See documentation
!  in circle_mesh.f90, but briefly the polar map is a map from the
!  circle to the domain, and the circle will be discretized with nrts(1,1),
!  points in a annular region outside a square of size rmid, and the 
!  square is discretized using nmid \times nmid patches
!  
!  Input arguments:
!    - nch2d: integer
!        number of chunks describing the generating curve
!    - tchse: real *8 (nch2d+1)
!        starting and ending location of in parameter space
!        for each 2d chunk
!    - k: integer
!        number of points per chunk (for the 2d curve)
!    - fcurve: function handle
!        function handle for corresponding to generating curve.
!        Should have calling sequence
!        fcurve(s, np, pars, r, z, drds, dzds, d2rds2, d2zds2)
!    - np: integer
!        number of parameters in fcurve
!    - pars: real *8 (np)
!        parameters of fcurve
!    - nrts: interger(2,nch2d)
!        number of quadrilater patches in the rz (s) direction,
!        and the periodic (t) direction.
!        if iptype is 1, then each quadrilateral patch
!        is subdivided into two triangular patches 
!    - rmid: real *8
!        radius of the inner square on the unit circle mesh
!        at the polar caps, should be between 0 and 1/sqrt(2)
!    - nmid: integer
!        * number of intervals in the middle square on the
!          unit circle mesh at the polar caps
!    - iort: integer
!        orientiation flag
!        if iort = 1, then parameter space is [a,b] \times [0,2\pi)
!        if iort = -1, then parameter space is [a,b] \times [2\pi,0)
!    - norder: integer
!        order of discretization on each patch
!    - iptype0: integer
!        type of patch to be used in the discretization
!        * iptype = 1, triangular patch discretized using RV nodes
!        * iptype = 11, quadrangular patch discretized with GL nodes
!        * iptype = 12, quadrangular patch discretized with Chebyshev 
!    - npatches: integer
!        number of patches
!    - npts: integer
!        Number of points on the discretized surface
!
!  Output arguments:
!    - norders: integer(npatches)
!        order of discretization on each patch
!        norders(i) = norder, for all i=1,2,...,npatches
!    - ixyzs: integer(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!        iptype = 11, quadrangular patch discretized with GL nodes
!        iptype = 12, quadrangular patch discretized with Chebyshev nodes
!        iptype(i) = iptype0, for all i=1,2,...,npatches
!    - srccoefs: real *8 (9,npts)
!        basis expansion coefficients of xyz, dxyz/du,
!        and dxyz/dv on each patch. 
!        For each point 
!          * srccoefs(1:3,i) is xyz info
!          * srccoefs(4:6,i) is dxyz/du info
!          * srccoefs(7:9,i) is dxyz/dv info
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the 
!        discretization nodes on the surface
!          * srcvals(1:3,i) - xyz info
!          * srcvals(4:6,i) - dxyz/du info
!          * srcvals(7:9,i) - dxyz/dv info
!          * srcvals(10:12,i) - normals info
!
!
      implicit none
      integer, intent(in) :: nch2d, k
      real *8, intent(in) :: tchse(nch2d+1)
      integer, intent(in) :: nrts(2,nch2d), nmid
      integer, intent(in) :: iptype0, norder
      integer, intent(in), target :: np
      integer, intent(in) :: iort
      real *8, intent(in), target :: pars(np)
      real *8, intent(in) :: rmid
      integer, intent(in) :: npatches, npts

      integer, intent(out) :: norders(npatches), iptype(npatches)
      integer, intent(out) :: ixyzs(npatches+1)
      real *8, intent(out) :: srccoefs(9,npts), srcvals(12,npts)

      integer i, npat0

      external fcurve
      external xtri_axissym_fun_eval, xtri_axissym_fun_circ_eval
      external xquad_axissym_fun_eval, xquad_axissym_fun_circ_eval
      
      procedure (), pointer :: patchpnt, patchpnt_circ

      real *8 done, pi
      integer npols
      real *8, allocatable, target :: pars1(:), pars_end(:)
      real *8, allocatable, target :: skel(:,:,:)
      integer, allocatable, target :: ipars1(:), ipars_end(:)
      integer, allocatable :: ixys1(:), ixys_end(:)
      real *8, allocatable :: ptinfo1(:,:), ptinfo_end(:,:)
      real *8, allocatable, target :: ptcoefs1(:,:), ptcoefs_end(:,:)
      integer, allocatable :: iptype1(:), iptype_end(:)
      integer, allocatable :: norders1(:), norders_end(:)


      real *8, allocatable :: uvs(:,:), wts(:), umatr(:,:), vmatr(:,:)
      integer npts1, npts_end
      integer npatches1, npatches_end
      real *8, pointer :: ptr1, ptr2, ptr3, ptr4
      integer, pointer :: iptr1, iptr2, iptr3, iptr4
      real *8 umin, umax, vmin, vmax
      integer npatches0, istart, npars(3)
      integer nover
      integer ich, iort_use
      real *8 dh


      done = 1.0d0
      pi = atan(done)*4
      
      
      call get_npols(iptype0, norder, npols)
      allocate(uvs(2,npols), wts(npols), umatr(npols,npols))
      allocate(vmatr(npols,npols))
      call get_disc_exps(norder, npols, iptype0, uvs, umatr, vmatr, wts)

      allocate(skel(3,3,npatches))


      npatches0 = 0
      istart = 1
      vmin = 0
      vmax = 2*pi

      if (iort.eq.-1) then
        vmin = 2*pi
        vmax = 0
      endif

      nover = 0

      if (iptype0.eq.1) then
        do ich = 2,nch2d-1
          umin = tchse(ich)
          umax = tchse(ich+1)
          call xtri_rectmesh_ani(umin, umax, vmin, vmax, nrts(1,ich), &
            nrts(2,ich), nover, npatches, npatches0, skel(1,1,istart))
          istart = istart + npatches0
        enddo
        npatches0 = istart - 1
        

        patchpnt => xtri_axissym_fun_eval
        patchpnt_circ => xtri_axissym_fun_circ_eval
      elseif (iptype0.eq.11.or.iptype0.eq.12) then
        do ich = 2,nch2d-1  
          umin = tchse(ich)
          umax = tchse(ich+1)
          call xquad_rectmesh_ani(umin, umax, vmin, vmax, nrts(1,ich), &
            nrts(2,ich), nover, npatches, npatches0, skel(1,1,istart))
          istart = istart + npatches0
        enddo
        npatches0 = istart-1
        patchpnt => xquad_axissym_fun_eval
        patchpnt_circ => xquad_axissym_fun_circ_eval
      endif

!
!  Now setup stuff for top and bottom polar caps
!
      npars(1) = nmid
      npars(2) = nrts(1,1)
      npars(3) = nrts(2,1)
      call mesh_circle_pts_mem(rmid, npars, iort, norder, iptype0, &
        npatches1, npts1)
      allocate(ixys1(npatches1+1), ptinfo1(6,npts1))
      allocate(ptcoefs1(6,npts1), iptype1(npatches1))
      allocate(norders1(npatches1))
      
      call mesh_circle_pts(rmid, npars, iort, norder, iptype0, & 
        npatches1, npts1, ixys1, ptinfo1)

      do i=1,npatches1
        iptype1(i) = iptype0
        norders1(i) = norder
      enddo

      allocate(pars1(np+2), pars_end(np+2))
      pars1(1) = tchse(1)
      pars1(2) = tchse(2)

      pars1(3:(np+2)) = pars(1:np)

      pars_end(1) = tchse(nch2d+1)
      pars_end(2) = tchse(nch2d)

      pars_end(3:(np+2)) = pars(1:np)



       
      call surf_vals_to_coefs(6, npatches1, norders1, ixys1, iptype1, &
        npts1, ptinfo1, ptcoefs1)

      npars(1) = nmid
      npars(2) = nrts(1,nch2d)
      npars(3) = nrts(2,nch2d)

      iort_use = iort*-1
      call mesh_circle_pts_mem(rmid, npars, iort_use, norder, iptype0, &
        npatches_end, npts_end)
      allocate(ixys_end(npatches_end+1), ptinfo_end(6,npts_end))
      allocate(ptcoefs_end(6,npts_end), iptype_end(npatches_end))
      allocate(norders_end(npatches_end))

      call mesh_circle_pts(rmid, npars, iort_use, norder, iptype0, & 
        npatches_end, npts_end, ixys_end, ptinfo_end)

      do i=1,npatches_end
        iptype_end(i) = iptype0
        norders_end(i) = norder
      enddo
      call surf_vals_to_coefs(6, npatches_end, norders_end, ixys_end, &
        iptype_end, npts_end, ptinfo_end, ptcoefs_end)

!
!  Let the discretrization begin
!
      allocate(ipars1(npatches1+4), ipars_end(npatches_end+4))
      ipars1(1) = np
      ipars1(2) = norder
      ipars1(3) = iptype0

      ipars_end(1) = np
      ipars_end(2) = norder
      ipars_end(3) = iptype0

      do i = 1,npatches1+1
        ipars1(i+3) = ixys1(i)
      enddo

      do i = 1,npatches_end+1
        ipars_end(i+3) = ixys_end(i)
      enddo
!
! top cap
!
      ptr1 => ptcoefs1(1,1)
      iptr2 => ipars1(1)
      ptr3 => pars1(1)


      call getgeominfo(npatches1, patchpnt_circ, ptr1, iptr2, ptr3, &
        fcurve, npols, uvs, umatr, srcvals, srccoefs)

!
! middle cylinder
!
      ptr1 => skel(1,1,1)
      iptr2 => np 
      ptr3 => pars(1)


      call getgeominfo(npatches0, patchpnt, ptr1, iptr2, ptr3, &
        fcurve, npols, uvs, umatr, srcvals(1,npts1+1), &
        srccoefs(1,npts1+1))
      
      istart = npts1 + npatches0*npols

!
! bottom cap
!
      ptr1 => ptcoefs_end(1,1)
      iptr2 => ipars_end(1)
      ptr3 => pars_end(1)

      call getgeominfo(npatches_end, patchpnt_circ, ptr1, iptr2, ptr3, &
        fcurve, npols, uvs, umatr, srcvals(1,istart+1), &
        srccoefs(1,istart+1))

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,npatches
        norders(i) = norder
        iptype(i) = iptype0
        ixyzs(i) = (i-1)*npols + 1
      enddo
!$OMP END PARALLEL DO
      
      ixyzs(npatches+1) = npts + 1

      return
      end

