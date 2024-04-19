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
!
!  Currently, only the _npat routines are available.
!  All of the routines have an extra _mem routines to get
!  the number of patches in the final discretization
!  
      subroutine get_ellipsoid_npat_mem(abc, nabc, c0, norder, &
        iptype0, npatches, npts)
!
!  Estimate the number of patches in the discretization of an ellipsoid
!  given  by
!
!  ((x - x0)/a)**2 + ((y-y0)/b)**2 + ((z-z0)/c)**2 = 1
!
!  where c0 = (x0, y0, z0), abc = (a,b,c)
!
!  The surface is discretized by projecting to the surface from
!  a rectangular piped with half-widths a/sqrt(3), b/sqrt(3), c/sqrt(3).
!
!  nabc determines the number of quadrilateral patches on the faces of the
!  rectangles of the rectangular parallelopipeds, with na patches
!  in the x direction, nb patches in the y direction, and nz patches
!  in the z direction, with nabc = (na, nb, nc)
!
!  Input arugments:
!    - abc: real *8(3)
!        semi major axes (a,b,c) along the x,y, and z directions
!    - nabc: integer(3)
!        number of patches along the coordinate directions of the
!        rectangular parellopiped (na, nb, nc)
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
!  Estimate the number of patches in the discretization of an ellipsoid
!  given  by
!
!  ((x - x0)/a)**2 + ((y-y0)/b)**2 + ((z-z0)/c)**2 = 1
!
!  where c0 = (x0, y0, z0), abc = (a,b,c)
!
!  The surface is discretized by projecting to the surface from
!  a rectangular piped with half-widths a/sqrt(3), b/sqrt(3), c/sqrt(3).
!
!  nabc determines the number of quadrilateral patches on the faces of the
!  rectangles of the rectangular parallelopipeds, with na patches
!  in the x direction, nb patches in the y direction, and nz patches
!  in the z direction, with nabc = (na, nb, nc)
!
!  Input arugments:
!    - abc: real *8(3)
!        semi major axes (a,b,c) along the x,y, and z directions
!    - nabc: integer(3)
!        number of patches along the coordinate directions of the
!        rectangular parellopiped (na, nb, nc)
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
      real *8, intent(in), target :: abc(3), c0(3)
      integer, intent(in) :: iptype0, nabc(3), norder, npatches, npts

      integer, intent(out) :: norders(npatches), iptype(npatches)
      integer, intent(out) :: ixyzs(npatches+1)
      real *8, intent(out) :: srccoefs(9,npts), srcvals(12,npts)

      real *8, allocatable, target :: skel(:,:,:)
      real *8, pointer :: ptr1, ptr2, ptr3, ptr4
      real *8, allocatable :: uvs(:,:), umatr(:,:), wts(:), vmatr(:,:)
      real *8 abcuse(3)
      real *8 ause, buse, cuse
      real *8, target :: p4(1)
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
      ptr1 => skel(1,1,1)
      ptr2 => abc(1)
      ptr3 => c0(1)
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

