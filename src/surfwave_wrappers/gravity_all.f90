!
!  Near field quadrature for the gravity surface wave kernels
!

      subroutine getnearquad_gravity_all(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, col_ind, &
        iquad, rfac0, nquad, iker, wnear)
!
!  This subroutine generates the near field quadrature
!  for the gravity surface wave kernel selected by iker.
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
!        kernel: 0 = G_S, 1 = G_phi
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
      subroutine getnearquad_gravity_grad_all(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, iker, wnear)
!
!  Near field quadrature for the target gradient (d/dx, d/dy) of the
!  kernels of getnearquad_gravity_all.
!
!  Input arguments: same as getnearquad_gravity_all, plus
!    - iker: integer
!        index of kernel, iker = 0 -> grad gs, iker = 1 -> grad gphi
!
!  Output arguments
!    - wnear: complex *16(2,nquad)
!        wnear(1,:) - d/dx of the kernel, near field correction
!        wnear(2,:) - d/dy of the kernel, near field correction
!
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
      complex *16, intent(out) :: wnear(2,nquad)

      integer *8 ipars(2)
      real *8 dpars(1)

      integer *8 ipv, i, ndi, ndd, ndz
      integer *8 ndtarg, ntarg
      integer *8 iker

      complex *16, allocatable :: wneartmp(:)

      procedure (), pointer :: fker
      external gsgravgradxkern, gsgravgradykern
      external gphigravgradxkern, gphigravgradykern

      ndz=0
      ndd=0
      ndi=0

      allocate(wneartmp(nquad))

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,nquad
        wnear(1,i) = 0
        wnear(2,i) = 0
      enddo
!$OMP END PARALLEL DO

      if (iquadtype.eq.1) then
        ipv=0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wneartmp(i) = 0
        enddo
!$OMP END PARALLEL DO

        if (iker.eq.0) then
          fker => gsgravgradxkern
        elseif (iker.eq.1) then
          fker => gphigravgradxkern
        endif
        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, zpars, &
          ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
          wneartmp)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wnear(1,i) = wneartmp(i)
        enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wneartmp(i) = 0
        enddo
!$OMP END PARALLEL DO

        if (iker.eq.0) then
          fker => gsgravgradykern
        elseif (iker.eq.1) then
          fker => gphigravgradykern
        endif
        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, zpars, &
          ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
          wneartmp)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wnear(2,i) = wneartmp(i)
        enddo
!$OMP END PARALLEL DO

      endif

      return
      end subroutine getnearquad_gravity_grad_all
!
!
!
