!
!
!  Make notes of pde being solved, and representation used
!
!

      subroutine getnearquad_gravity_all(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, &
        iquad, rfac0, nquad, iker, wnear)
!
!  This subroutine generates the near field quadrature
!  for the representation:
!
!  u = (Enter representation here) 
!
!  On imposing the boundary condition, we get the following operator
!
!  du/dn + ik \lambda u =  
!    z \sigma + S_{k}'[\sigma] + i\alpha S_{i|k|}'^2 [\sigma] + i \alpha 
!       (D_{k}' - D_{i|k|}') S_{i|k|}[\sigma]  + 
!       ik \lambda (S_{k} + i \alpha D_{k} S_{i|k|} + 
!       i \alpha w S_{i|k|}) = f
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
!        iptype = 11, quadrangular patch discretized with GL nodes
!        iptype = 12, quadrangular patch discretized with Chebyshev 
!                     nodes
!    - npts: integer
!        total number of discretization points on the boundary
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
!    - eps: real *8
!        precision requested
!    - zpars: complex *16(*)
!        kernel parameters (Referring to formula (1))
!        fix documentation here
!    - iquadtype: integer
!        quadrature type
!          * iquadtype = 1, use ggq for self + adaptive integration
!            for rest
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
!    - rfac0: real *8
!        radius parameter for switching to predetermined quadarature
!        rule        
!    - nquad: integer
!        number of near field entries corresponding to each source target
!        pair
!    - iker: integer
!        index of kernel, iker = 0 -> gs, iker = 1 -> gphi, iker = 2 -> lapgs
!
!  Output arguments
!    - wnear: complex *16(nquad)
!        The desired near field quadrature
!        stores the quadrature corrections for <enter kernel here> 
  
      implicit none 
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches), ipatch_id(ntarg)
      real *8, intent(in) ::  srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps, uvs_targ(2,ntarg), targs(ndtarg,ntarg)
      complex *16, intent(in) :: zpars(6)
      integer *8, intent(in) :: iquadtype, nnz
      integer *8, intent(in) :: row_ptr(npts+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: nquad
      complex *16, intent(out) :: wnear(nquad)
      
      complex *16 zpars_tmp(3)
      integer *8 ipars(2)
      real *8 dpars(1)
      
      integer *8 ipv, i, ndi, ndd, ndz
      
      integer *8 ndtarg, ntarg

      integer *8 iker

      procedure (), pointer :: fker
      external gphigravkern, gsgravkern

      ndz=0
      ndd=0
      ndi=0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,nquad
        wnear(i) = 0
      enddo
!$OMP END PARALLEL DO

      if (iquadtype.eq.1) then
        ipv=0

        if (iker.eq.0) then
          fker => gsgravkern 
        elseif (iker.eq.1) then
          fker => gphigravkern
        endif
        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, zpars, &
          ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
          wnear)

      endif

      return
      end subroutine getnearquad_gravity_all
!
!
!
      subroutine getnearquad_gravity_all_vpp(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, zpars0, iquadtype, nnz, row_ptr, col_ind, &
        iquad, rfac0, nquad, iker, maxdist, wnear)
!
!  This subroutine generates the near field quadrature
!  for the representation:
!
!  u = (Enter representation here) 
!
!  On imposing the boundary condition, we get the following operator
!
!  du/dn + ik \lambda u =  
!    z \sigma + S_{k}'[\sigma] + i\alpha S_{i|k|}'^2 [\sigma] + i \alpha 
!       (D_{k}' - D_{i|k|}') S_{i|k|}[\sigma]  + 
!       ik \lambda (S_{k} + i \alpha D_{k} S_{i|k|} + 
!       i \alpha w S_{i|k|}) = f
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
!        iptype = 11, quadrangular patch discretized with GL nodes
!        iptype = 12, quadrangular patch discretized with Chebyshev 
!                     nodes
!    - npts: integer
!        total number of discretization points on the boundary
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
!    - eps: real *8
!        precision requested
!    - zpars: complex *16(*)
!        kernel parameters (Referring to formula (1))
!        fix documentation here
!    - iquadtype: integer
!        quadrature type
!          * iquadtype = 1, use ggq for self + adaptive integration
!            for rest
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
!    - rfac0: real *8
!        radius parameter for switching to predetermined quadarature
!        rule        
!    - nquad: integer
!        number of near field entries corresponding to each source target
!        pair
!    - iker: integer
!        index of kernel, iker = 0 -> gs, iker = 1 -> gphi, iker = 2 -> lapgs
!
!  Output arguments
!    - wnear: complex *16(nquad)
!        The desired near field quadrature
!        stores the quadrature corrections for <enter kernel here> 
  
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches), ipatch_id(ntarg)
      real *8, intent(in) ::  srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps, uvs_targ(2,ntarg), targs(ndtarg,ntarg)
      complex *16, intent(in) :: zpars0(6)
      integer *8, intent(in) :: iquadtype, nnz
      integer *8, intent(in) :: row_ptr(npts+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: rfac0, maxdist
      integer *8, intent(in) :: nquad
      complex *16, intent(out) :: wnear(nquad)
      complex *16 zpars_tmp(3), zpars
      integer *8 ipars(5), ipars0
      integer *8, parameter :: ldpars = 10000000
      real *8 :: dpars0
      real *8, allocatable :: dpars(:)
      real *8 :: tol,pi,done,a,b
      integer *8 ipv, i, ndi, ndd, ndz, ndi0, ndd0, ndz0
      integer *8 ndtarg, ntarg, n, nf, maxsub, maxdepth, ietype
      integer *8 iker, ier
      procedure (), pointer :: fker, fker0
      external gphigravkern, gsgravkern, vpp_kern

      allocate(dpars(ldpars))

!      parameters defining the true kernel
      ndz0=0
      ndd0=0
      ndi0=0
      

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,nquad
        wnear(i) = 0
      enddo
!$OMP END PARALLEL DO

      if (iquadtype.eq.1) then
        ipv=0

        done = 1
        pi = atan(done)*4
        nf = 2
        a = 1d-12
        b = maxdist
        n = 16
        tol = 1d-12
        maxsub = 1000
        maxdepth = 50
        ietype = 1
        ier = 0

        if (iker.eq.0) then
          fker0 => gsgravkern 
        elseif (iker.eq.1) then
          fker0 => gphigravkern
        endif
        call vpp_buildkern(fker0,ndd0,dpars0,ndz0,zpars0,ndi0, &
          ipars0,nf,a,b,n,tol,maxsub,maxdepth,ietype, &
          ipars,ndd,ldpars,dpars,ier)

        fker => vpp_kern

        ndi = 5
        ndz = 0
        zpars = 0

        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, zpars, &
          ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
          wnear)

      endif

      return
      end subroutine getnearquad_gravity_all_vpp
!
!
!
