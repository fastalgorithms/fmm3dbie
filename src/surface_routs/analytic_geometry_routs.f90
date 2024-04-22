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
!         * radii(1) = rminor
!         * radii(2) = rmajor
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
      coefs(1+nosc+1,1,1) = radii(3)/2
      coefs(1+abs(nosc-1),1,1) = radii(3)/2
      
      coefs(m+2,1,3) = radii(2)
      coefs(2*m+1,1,3) = radii(3)/2
      if (nosc-1.gt.0) coefs(m+1+nosc-1,1,3) = -radii(3)/2
      if (nosc-1.lt.0) coefs(m+1+abs(nosc-1),1,3) = radii(3)/2

      iort = 1

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
