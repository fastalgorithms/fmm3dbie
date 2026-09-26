!
!  Near field quadrature for the flexural surface wave kernels
!

      subroutine getnearquad_flex_all(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, &
        eps, zpars, iquadtype, nnz, row_ptr, col_ind, &
        iquad, rfac0, nquad, iker, wnear)
!
!  This subroutine generates the near field quadrature
!  for the flexural surface wave kernel selected by iker.
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
!        kernel: 1 = G_S, 2 = G_phi, 3 = bilap G_S, 4 = bilap G_phi,
!        5 = S3d G_phi
!
!  Output arguments
!    - wnear: complex *16(nquad)
!        The desired near field quadrature
!        stores the quadrature corrections for <enter kernel here> 
  
      implicit none 
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8, intent(in) ::  srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps
      complex *16, intent(in) :: zpars(13)
      integer *8, intent(in) :: iquadtype, nnz
      integer *8, intent(in) :: row_ptr(npts+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: nquad
      complex *16, intent(out) :: wnear(nquad)
      complex *16 zpars_tmp(3)
      integer *8 ipars(2)
      real *8 dpars(1)
      real *8, allocatable :: uvs_targ(:,:)
      integer *8, allocatable :: ipatch_id(:)
      integer *8 ipv, i, ndi, ndd, ndz
      integer *8 ndtarg, ntarg
      integer *8 iker
      procedure (), pointer :: fker
      external gphiflexkern, gsflexkern
      external bilapgsflexkern, bilapgphiflexkern
      external s3dgphiflexkern
      
      ndz=13
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

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,nquad
        wnear(i) = 0
      enddo
!$OMP END PARALLEL DO

      call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts, &
        ipatch_id, uvs_targ)
      if (iquadtype.eq.1) then
        ipv=0

        if (iker.eq.1) then
          fker => gsflexkern 
        elseif (iker.eq.2) then
          fker => gphiflexkern
        elseif (iker.eq.3) then
          fker => bilapgsflexkern
        elseif (iker.eq.4) then
          fker => bilapgphiflexkern
        elseif (iker.eq.5) then 
          fker => s3dgphiflexkern
        endif
        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, srcvals, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, zpars, &
          ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
          wnear)

      endif

      return
      end subroutine getnearquad_flex_all
!
!
!
!
!
!  getnearquad_flex_all at arbitrary targets. row_ptr, col_ind and
!  iquad are sized by ntarg.
!
!  Additional input arguments (relative to getnearquad_flex_all):
!    - ndtarg: integer
!        leading dimension of targvals
!    - ntarg: integer
!        number of targets
!    - targvals: real *8 (ndtarg,ntarg)
!        target locations (first 3 components must be xyz)
!    - ipatch_id_targ: integer(ntarg)
!        patch id of target if on-surface, else -1
!    - uvs_targ: real *8(2,ntarg)
!        local uv coordinates of target on patch ipatch_id_targ(i),
!        if on-surface
!
!  row_ptr: integer(ntarg+1), col_ind: integer(nnz), iquad: integer(nnz+1)
!
      subroutine getnearquad_flex_all_targ(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, &
        ndtarg, ntarg, targvals, ipatch_id_targ, uvs_targ, &
        eps, zpars, iquadtype, nnz, row_ptr, col_ind, &
        iquad, rfac0, nquad, iker, wnear)

      implicit none
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8, intent(in) ::  srccoefs(9,npts), srcvals(12,npts)
      integer *8, intent(in) :: ndtarg, ntarg
      real *8, intent(in) :: targvals(ndtarg,ntarg)
      integer *8, intent(in) :: ipatch_id_targ(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      real *8, intent(in) :: eps
      complex *16, intent(in) :: zpars(13)
      integer *8, intent(in) :: iquadtype, nnz
      integer *8, intent(in) :: row_ptr(ntarg+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: nquad
      complex *16, intent(out) :: wnear(nquad)
      integer *8 ipars(2)
      real *8 dpars(1)
      integer *8 ipv, i, ndi, ndd, ndz
      integer *8 iker
      procedure (), pointer :: fker
      external gphiflexkern, gsflexkern
      external bilapgsflexkern, bilapgphiflexkern
      external s3dgphiflexkern

      ndz=13
      ndd=0
      ndi=0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,nquad
        wnear(i) = 0
      enddo
!$OMP END PARALLEL DO

      if (iquadtype.eq.1) then
        ipv=0

        if (iker.eq.1) then
          fker => gsflexkern
        elseif (iker.eq.2) then
          fker => gphiflexkern
        elseif (iker.eq.3) then
          fker => bilapgsflexkern
        elseif (iker.eq.4) then
          fker => bilapgphiflexkern
        elseif (iker.eq.5) then
          fker => s3dgphiflexkern
        endif
        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targvals, &
          ipatch_id_targ, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, &
          zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, &
          nquad, wnear)

      endif

      return
      end subroutine getnearquad_flex_all_targ
