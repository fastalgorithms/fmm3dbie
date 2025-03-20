      
      subroutine getnearquad_helm_open_surface_disk(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, &
        eps, ifscal, zk, iquadtype, nnz, row_ptr, col_ind, &
        iquad, rfac0, nquad, wnear)
!
!  This subroutine generates the near field quadrature
!  for one of the two representations
!
!  u = S[\sigma]  if ifscal = 0
!  u = S[1/sqrt(1-r^2) \sigma]  if ifscal = 1
!
!  required to simulate the Dirichlet problem on an open surface
!
!  The quadrature is computed by the following strategy
!  targets within a sphere of radius rfac0*rs
!  of a patch centroid is handled using adaptive integration
!  where rs is the radius of the bounding sphere
!  for the patch
!
!  All other targets in the near field are handled via
!  oversampled quadrature
!
!  The recommended parameter for rfac0 is 1.25d0
!
!  Input arguments:
!
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization on each patch
!    - ixyzs: integer(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!    - npts: integer
!        total number of discretization points on the boundary
!    - srccoefs: real *8 (9,npts)
!        koornwinder expansion coefficients of xyz, dxyz/du,
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
!    - eps: real *8
!        precision requested
!    - ifscal: integer
!        whether to scale the representation with sqrt of distance 
!        to boundary
!    - zk: complex *16 
!        kernel parameters (See notes above)
!    - iquadtype: integer
!        quadrature type
!          * iquadtype = 1, use ggq for self + adaptive integration
!            for rest
!          * iquadtype = 2, use adap for all (recommended if ifscal = 1)
!    - nnz: integer
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer(npts+1)
!        row_ptr(i) is the pointer
!        to col_ind array where list of relevant source patches
!        for target i start
!    - col_ind: integer (nnz)
!        list of source patches relevant for all targets, sorted
!        by the target number
!    - iquad: integer(nnz+1)
!        location in wnear_ij array where quadrature for col_ind(i)
!        starts for a single kernel. In this case the different kernels
!        are matrix entries are located at (m-1)*nquad+iquad(i), where
!        m is the kernel number
!    - rfac0: integer
!        radius parameter for near field
!    - nquad: integer
!        number of near field entries corresponding to each source target
!        pair. 
!
!  Output arguments
!    - wnear: complex *16(nquad)
!        The desired near field quadrature
!
!

      implicit none
      integer npatches,norders(npatches),npts,nquad
      integer ixyzs(npatches+1),iptype(npatches)
      integer ifscal
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,rfac0
      integer ndtarg,ntarg
      integer iquadtype
      complex *16 zk 
      integer nnz,ipars(2)
      real *8 dpars(1)
      integer row_ptr(npts+1),col_ind(nnz),iquad(nnz+1)
      complex *16 wnear(nquad)

      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)


      complex *16 ima,z
      integer ndi,ndd,ndz,i

      integer ipv

      procedure (), pointer :: fker
      external h3d_slp, h3d_slp_disk

      data ima/(0.0d0,1.0d0)/

      if(ifscal.eq.0) then
        ipv = 0
        fker => h3d_slp
      else
        ipv = 1
        fker => h3d_slp_disk
      endif

      ndz=1
      ndd=0
      ndi=0
      ndtarg = 12
      ntarg = npts


      allocate(ipatch_id(npts),uvs_targ(2,npts))
!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo
!$OMP END PARALLEL DO

      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,&
        ipatch_id,uvs_targ)

      if (iquadtype.eq.1) then

        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, srcvals, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, zk, &
          ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
          wnear)
      elseif (iquadtype.eq.2) then
        call zgetnearquad_adap_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, srcvals, &
          ipatch_id, uvs_targ, eps, fker, ndd, dpars, ndz, zk, &
          ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
          wnear)
      endif


      return
      end subroutine getnearquad_helm_open_surface_disk
!





      subroutine h3d_slp_disk(src, ndt, targ, ndd, dpars, ndz, zk, &
         ndi, ipars, val)
      implicit real *8 (a-h, o-z)
      real *8 src(*), targ(ndt), dpars(ndd) 
      real *8 over4pi
      integer ipars(ndi)
      complex *16 :: zk, val
      complex *16 :: ima

      data ima/(0.0d0,1.0d0)/
      data over4pi/0.07957747154594767d0/
!
!   returns the helmholtz single layer potential kernel
!   weighted by distance to boundary
!

      dx=targ(1)-src(1)
      dy=targ(2)-src(2)
      dz=targ(3)-src(3)

      r=sqrt(dx**2 + dy**2 + dz**2)
      dd = sqrt(1.0d0 - src(1)**2 - src(2)**2)
      if(dd.le.0) dd = 1e-8
      rbdry = 1.0d0/dd

      val =  exp(ima*zk*r)/r*over4pi*rbdry

      return
      end subroutine h3d_slp_disk

