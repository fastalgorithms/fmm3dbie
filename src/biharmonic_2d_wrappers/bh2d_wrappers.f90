!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!   Near quadrature wrappers for biharmonic equation in 2d
!
!  PDE:
!    \Delta^2 u = 0
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  Subroutines defined in this file:
!    - getnearquad_bh2d_gv2v: near quadrature for the biharmonic Green's
!        function on-surface (volume-to-volume)
!
!    - getnearquad_bh2d_dir: near quadrature for the biharmonic Green's
!        function, Dirichlet targets
!
!    - getnearquad_bh2d_neu: near quadrature for the normal derivative
!        of the biharmonic Green's function, Neumann targets
!
!    - getnearquad_bh2d_supp2: near quadrature for the simply-supported
!        plate boundary operator (2nd-order term)
!
!    - getnearquad_bh2d_free2: near quadrature for the free plate
!        boundary operator (2nd-order term)
!
!    - getnearquad_bh2d_clamped: near quadrature for the clamped plate
!        representation (both bh2d_g and bh2d_gdn)
!
!    - getnearquad_bh2d_b2v_dir: near quadrature for the b2v Dirichlet
!        representation using bh2d_g
!
!    - getnearquad_bh2d_v2b_neu: near quadrature for the v2b Neumann
!        representation using bh2d_gdn
!
!    - getnearquad_bh2d_gdn: near quadrature for the normal derivative
!        of the biharmonic Green's function
!

      subroutine getnearquad_bh2d_gv2v(npatches, norders, &
       ixyzs, iptype, npts, srccoefs, srcvals, eps, iquadtype, nnz, &
       row_ptr, col_ind, iquad, rfac0, nquad, wnear)
!
!  This subroutine generates the near field quadrature
!  for the representation:
!
!  This subroutine generates the near field quadrature
!  for the representation:
!
!  u = G_{bh}[\sigma]
!
!  on-surface (targets = source points)
!
!  Input arguments:
!    - npatches: integer *8
!        number of patches
!    - norders: integer *8(npatches)
!        order of discretization on each patch
!    - ixyzs: integer *8(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer *8(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!        iptype = 11, quadrangular patch discretized with GL nodes
!        iptype = 12, quadrangular patch discretized with Chebyshev
!                     nodes
!    - npts: integer *8
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
!    - iquadtype: integer *8
!        quadrature type
!          * iquadtype = 1, use ggq for self + adaptive integration
!            for rest
!    - nnz: integer *8
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer *8(npts+1)
!        row_ptr(i) is the pointer
!        to col_ind array where list of relevant source patches
!        for target i start
!    - col_ind: integer *8 (nnz)
!        list of source patches relevant for all targets, sorted
!        by the target number
!    - iquad: integer *8(nnz+1)
!        location in wnear_ij array where quadrature for col_ind(i)
!        starts for a single kernel. In this case the different kernels
!        are matrix entries are located at (m-1)*nquad+iquad(i), where
!        m is the kernel number
!    - rfac0: real *8
!        radius parameter for near field
!    - nquad: integer *8
!        number of near field entries corresponding to each source target
!        pair.
!
!  Output arguments
!    - wnear: real *8(nquad)
!        The desired near field quadrature
!               
!

      implicit none 
      integer *8, intent(in) :: npatches, norders(npatches)
      integer *8, intent(in) :: ixyzs(npatches+1), iptype(npatches)
      integer *8, intent(in) :: npts
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps
      integer *8, intent(in) :: iquadtype, nnz
      integer *8, intent(in) :: row_ptr(npts+1), col_ind(nnz), iquad(nnz+1)
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: nquad
      
      real *8, intent(out) :: wnear(nquad)

      integer *8 ndtarg, ntarg, ndd, ndz, ndi
      real *8 dpars(1)
      complex *16 zpars(1)
      integer *8 ipars(1)

      real *8, allocatable :: uvs_targ(:,:)
      integer *8, allocatable :: ipatch_id(:)

      integer *8 i
      integer *8 ipv

      procedure (), pointer :: fker
      external bh2d_g


      ndz=0
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

      call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts, &
        ipatch_id, uvs_targ)

      ipv=0
      fker => bh2d_g
      call dgetnearquad_ggq_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, srcvals, &
        ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear)

      return
      end subroutine getnearquad_bh2d_gv2v
!
!
!
!
!
      subroutine getnearquad_bh2d_dir(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)
!
!
!  This subroutine generates the near field quadrature
!  for the representation:
!
!  u = G_{bh}[\sigma]
!
!  Input arguments:
!    - npatches: integer *8
!        number of patches
!    - norders: integer *8(npatches)
!        order of discretization on each patch
!    - ixyzs: integer *8(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer *8(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!        iptype = 11, quadrangular patch discretized with GL nodes
!        iptype = 12, quadrangular patch discretized with Chebyshev
!                     nodes
!    - npts: integer *8
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
!    - ndtarg: integer *8
!        leading dimension of target information array
!    - ntarg: integer *8
!        number of targets
!    - targs: real *8(ndtarg,ntarg)
!        target information
!    - ipatch_id: integer *8(ntarg)
!        ipatch_id(i) indicates the patch on which target i
!        is, if it is on surface. ipatch_id(i) should be 0
!        otherwise
!    - uvs_targ: real *8(2,ntarg)
!        if ipatch_id(i) > 0, then uvs_targ(1:2,i) are the
!        local uv coordinates of the target on the patch,
!        uvs_targ(1:2,i) is unused otherwise
!    - eps: real *8
!        precision requested
!    - iquadtype: integer *8
!        quadrature type
!          * iquadtype = 1, use ggq for self + adaptive integration
!            for rest
!    - nnz: integer *8
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer *8(ntarg+1)
!        row_ptr(i) is the pointer
!        to col_ind array where list of relevant source patches
!        for target i start
!    - col_ind: integer *8 (nnz)
!        list of source patches relevant for all targets, sorted
!        by the target number
!    - iquad: integer *8(nnz+1)
!        location in wnear_ij array where quadrature for col_ind(i)
!        starts for a single kernel. In this case the different kernels
!        are matrix entries are located at (m-1)*nquad+iquad(i), where
!        m is the kernel number
!    - rfac0: real *8
!        radius parameter for switching to predetermined quadrature
!        rule
!    - nquad: integer *8
!        number of near field entries corresponding to each source target
!        pair
!
!  Output arguments
!    - wnear: real *8(nquad)
!        The desired near field quadrature
!
!

      implicit none
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8, intent(in) ::  srccoefs(9,npts), srcvals(12,npts)
      integer *8, intent(in) :: ndtarg, ntarg
      real *8, intent(in) :: targs(ndtarg, ntarg)
      integer *8, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      real *8, intent(in) :: eps
      integer *8, intent(in) :: iquadtype, nnz
      integer *8, intent(in) :: row_ptr(ntarg+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: nquad
      real *8, intent(out) :: wnear(nquad)

      integer *8 ipv, ndd, ndi, ndz, i
      real *8 dpars
      integer *8 ipars
      complex *16 zpars

      procedure (), pointer :: fker
      external bh2d_g

      ndd = 0
      ndz = 0
      ndi = 0



    
      if (iquadtype.eq.1) then
        ipv = 0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wnear(i) = 0 
        enddo
!$OMP END PARALLEL DO
        fker => bh2d_g
        call dgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, &
          zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, &
          rfac0, nquad, wnear)
      endif 

      return
      end
!
!
!
!
!
      subroutine getnearquad_bh2d_neu(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)
!
!
!  This subroutine generates the near field quadrature
!  for the representation:
!
!  u = \partial G_{bh}/\partial n_{y} [\sigma]
!
!  Input arguments:
!    - npatches: integer *8
!        number of patches
!    - norders: integer *8(npatches)
!        order of discretization on each patch
!    - ixyzs: integer *8(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer *8(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!        iptype = 11, quadrangular patch discretized with GL nodes
!        iptype = 12, quadrangular patch discretized with Chebyshev
!                     nodes
!    - npts: integer *8
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
!    - ndtarg: integer *8
!        leading dimension of target information array
!    - ntarg: integer *8
!        number of targets
!    - targs: real *8(ndtarg,ntarg)
!        target information
!    - ipatch_id: integer *8(ntarg)
!        ipatch_id(i) indicates the patch on which target i
!        is, if it is on surface. ipatch_id(i) should be 0
!        otherwise
!    - uvs_targ: real *8(2,ntarg)
!        if ipatch_id(i) > 0, then uvs_targ(1:2,i) are the
!        local uv coordinates of the target on the patch,
!        uvs_targ(1:2,i) is unused otherwise
!    - eps: real *8
!        precision requested
!    - iquadtype: integer *8
!        quadrature type
!          * iquadtype = 1, use ggq for self + adaptive integration
!            for rest
!    - nnz: integer *8
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer *8(ntarg+1)
!        row_ptr(i) is the pointer
!        to col_ind array where list of relevant source patches
!        for target i start
!    - col_ind: integer *8 (nnz)
!        list of source patches relevant for all targets, sorted
!        by the target number
!    - iquad: integer *8(nnz+1)
!        location in wnear_ij array where quadrature for col_ind(i)
!        starts for a single kernel. In this case the different kernels
!        are matrix entries are located at (m-1)*nquad+iquad(i), where
!        m is the kernel number
!    - rfac0: real *8
!        radius parameter for switching to predetermined quadrature
!        rule
!    - nquad: integer *8
!        number of near field entries corresponding to each source target
!        pair
!
!  Output arguments
!    - wnear: real *8(nquad)
!        The desired near field quadrature
!
!

      implicit none
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8, intent(in) ::  srccoefs(9,npts), srcvals(12,npts)
      integer *8, intent(in) :: ndtarg, ntarg
      real *8, intent(in) :: targs(ndtarg, ntarg)
      integer *8, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      real *8, intent(in) :: eps
      integer *8, intent(in) :: iquadtype, nnz
      integer *8, intent(in) :: row_ptr(ntarg+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: nquad
      real *8, intent(out) :: wnear(nquad)

      integer *8 ipv, ndd, ndi, ndz, i
      real *8 dpars
      integer *8 ipars
      complex *16 zpars

      procedure (), pointer :: fker
      external bh2d_gdn

      ndd = 0
      ndz = 0
      ndi = 0



    
      if (iquadtype.eq.1) then
        ipv = 0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wnear(i) = 0 
        enddo
!$OMP END PARALLEL DO
        fker => bh2d_gdn
        call dgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, &
          zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, &
          rfac0, nquad, wnear)

      endif 
      return
      end
!
!
!
!
!
      subroutine getnearquad_bh2d_supp2(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, dpars, zpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)
!
!
!  This subroutine generates the near field quadrature
!  for the simply-supported plate boundary operator (2nd-order term).
!
!  Input arguments:
!    - npatches: integer *8
!        number of patches
!    - norders: integer *8(npatches)
!        order of discretization on each patch
!    - ixyzs: integer *8(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer *8(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!        iptype = 11, quadrangular patch discretized with GL nodes
!        iptype = 12, quadrangular patch discretized with Chebyshev
!                     nodes
!    - npts: integer *8
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
!    - ndtarg: integer *8
!        leading dimension of target information array
!    - ntarg: integer *8
!        number of targets
!    - targs: real *8(ndtarg,ntarg)
!        target information
!    - ipatch_id: integer *8(ntarg)
!        ipatch_id(i) indicates the patch on which target i
!        is, if it is on surface. ipatch_id(i) should be 0
!        otherwise
!    - uvs_targ: real *8(2,ntarg)
!        if ipatch_id(i) > 0, then uvs_targ(1:2,i) are the
!        local uv coordinates of the target on the patch,
!        uvs_targ(1:2,i) is unused otherwise
!    - eps: real *8
!        precision requested
!    - dpars: real *8 (1)
!        real kernel parameters
!        dpars(1) = nu, Poisson ratio
!    - iquadtype: integer *8
!        quadrature type
!          * iquadtype = 1, use ggq for self + adaptive integration
!            for rest
!    - nnz: integer *8
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer *8(ntarg+1)
!        row_ptr(i) is the pointer
!        to col_ind array where list of relevant source patches
!        for target i start
!    - col_ind: integer *8 (nnz)
!        list of source patches relevant for all targets, sorted
!        by the target number
!    - iquad: integer *8(nnz+1)
!        location in wnear_ij array where quadrature for col_ind(i)
!        starts for a single kernel. In this case the different kernels
!        are matrix entries are located at (m-1)*nquad+iquad(i), where
!        m is the kernel number
!    - rfac0: real *8
!        radius parameter for switching to predetermined quadrature
!        rule
!    - nquad: integer *8
!        number of near field entries corresponding to each source target
!        pair
!
!  Output arguments
!    - wnear: real *8(nquad)
!        The desired near field quadrature
!
!

      implicit none
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8, intent(in) ::  srccoefs(9,npts), srcvals(12,npts)
      integer *8, intent(in) :: ndtarg, ntarg
      real *8, intent(in) :: targs(ndtarg, ntarg)
      integer *8, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      real *8, intent(in) :: eps
      integer *8, intent(in) :: iquadtype, nnz
      integer *8, intent(in) :: row_ptr(ntarg+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: nquad
      real *8, intent(out) :: wnear(nquad)

      integer *8 ipv, ndd, ndi, ndz, i
      real *8 dpars
      integer *8 ipars
      complex *16 zpars

      procedure (), pointer :: fker
      external bh2d_gsupp2

      ndd = 1
      ndz = 0
      ndi = 0



    
      if (iquadtype.eq.1) then
        ipv = 0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wnear(i) = 0 
        enddo
!$OMP END PARALLEL DO
        fker => bh2d_gsupp2
        call dgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, &
          zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, &
          rfac0, nquad, wnear)

      endif 
      return
      end
!
!
!
!
!
      subroutine getnearquad_bh2d_free2(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, dpars, zpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)
!
!
!  This subroutine generates the near field quadrature
!  for the free plate boundary operator (2nd-order term).
!
!  Input arguments:
!    - npatches: integer *8
!        number of patches
!    - norders: integer *8(npatches)
!        order of discretization on each patch
!    - ixyzs: integer *8(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer *8(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!        iptype = 11, quadrangular patch discretized with GL nodes
!        iptype = 12, quadrangular patch discretized with Chebyshev
!                     nodes
!    - npts: integer *8
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
!    - ndtarg: integer *8
!        leading dimension of target information array
!    - ntarg: integer *8
!        number of targets
!    - targs: real *8(ndtarg,ntarg)
!        target information
!    - ipatch_id: integer *8(ntarg)
!        ipatch_id(i) indicates the patch on which target i
!        is, if it is on surface. ipatch_id(i) should be 0
!        otherwise
!    - uvs_targ: real *8(2,ntarg)
!        if ipatch_id(i) > 0, then uvs_targ(1:2,i) are the
!        local uv coordinates of the target on the patch,
!        uvs_targ(1:2,i) is unused otherwise
!    - eps: real *8
!        precision requested
!    - dpars: real *8 (1)
!        real kernel parameters
!        dpars(1) = nu, Poisson ratio
!    - iquadtype: integer *8
!        quadrature type
!          * iquadtype = 1, use ggq for self + adaptive integration
!            for rest
!    - nnz: integer *8
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer *8(ntarg+1)
!        row_ptr(i) is the pointer
!        to col_ind array where list of relevant source patches
!        for target i start
!    - col_ind: integer *8 (nnz)
!        list of source patches relevant for all targets, sorted
!        by the target number
!    - iquad: integer *8(nnz+1)
!        location in wnear_ij array where quadrature for col_ind(i)
!        starts for a single kernel. In this case the different kernels
!        are matrix entries are located at (m-1)*nquad+iquad(i), where
!        m is the kernel number
!    - rfac0: real *8
!        radius parameter for switching to predetermined quadrature
!        rule
!    - nquad: integer *8
!        number of near field entries corresponding to each source target
!        pair
!
!  Output arguments
!    - wnear: real *8(nquad)
!        The desired near field quadrature
!
!

      implicit none
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8, intent(in) ::  srccoefs(9,npts), srcvals(12,npts)
      integer *8, intent(in) :: ndtarg, ntarg
      real *8, intent(in) :: targs(ndtarg, ntarg)
      integer *8, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      real *8, intent(in) :: eps
      integer *8, intent(in) :: iquadtype, nnz
      integer *8, intent(in) :: row_ptr(ntarg+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: nquad
      real *8, intent(out) :: wnear(nquad)

      integer *8 ipv, ndd, ndi, ndz, i
      real *8 dpars
      integer *8 ipars
      complex *16 zpars

      procedure (), pointer :: fker
      external bh2d_gfree2

      ndd = 1
      ndz = 0
      ndi = 0



    
      if (iquadtype.eq.1) then
        ipv = 0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wnear(i) = 0 
        enddo
!$OMP END PARALLEL DO
        fker => bh2d_gfree2
        call dgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, &
          zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, &
          rfac0, nquad, wnear)

      endif
      return
      end
!
!
!
!
!
      subroutine getnearquad_bh2d_clamped(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)
!
!  This subroutine generates the near field quadrature
!  for the clamped plate representation using both
!  the biharmonic Green's function (bh2d_g) and its normal derivative (bh2d_gdn).
!  wnear(1,:) contains quadrature for bh2d_g, wnear(2,:) for bh2d_gdn.
!
      implicit none
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      integer *8, intent(in) :: ndtarg, ntarg
      real *8, intent(in) :: targs(ndtarg,ntarg)
      integer *8, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      real *8, intent(in) :: eps
      integer *8, intent(in) :: iquadtype, nnz
      integer *8, intent(in) :: row_ptr(ntarg+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: nquad
      real *8, intent(out) :: wnear(2,nquad)

      integer *8 ipv, ndd, ndi, ndz, i
      real *8 dpars
      integer *8 ipars
      complex *16 zpars

      procedure (), pointer :: fker
      external bh2d_g, bh2d_gdn

      ndd = 0
      ndz = 0
      ndi = 0

      if (iquadtype.eq.1) then
        ipv = 0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wnear(1,i) = 0
          wnear(2,i) = 0
        enddo
!$OMP END PARALLEL DO

        fker => bh2d_g
        call dgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, &
          zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, &
          rfac0, nquad, wnear(1,:))

        fker => bh2d_gdn
        call dgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, &
          zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, &
          rfac0, nquad, wnear(2,:))
      endif

      return
      end
!
!
!
!
!
      subroutine getnearquad_bh2d_b2v_dir(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)
!
!  This subroutine generates the near field quadrature
!  for the b2v Dirichlet representation using bh2d_g.
!
      implicit none
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      integer *8, intent(in) :: ndtarg, ntarg
      real *8, intent(in) :: targs(ndtarg,ntarg)
      integer *8, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      real *8, intent(in) :: eps
      integer *8, intent(in) :: iquadtype, nnz
      integer *8, intent(in) :: row_ptr(ntarg+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: nquad
      real *8, intent(out) :: wnear(nquad)

      integer *8 ipv, ndd, ndi, ndz, i
      real *8 dpars
      integer *8 ipars
      complex *16 zpars

      procedure (), pointer :: fker
      external bh2d_g

      ndd = 0
      ndz = 0
      ndi = 0

      if (iquadtype.eq.1) then
        ipv = 0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wnear(i) = 0
        enddo
!$OMP END PARALLEL DO

        fker => bh2d_g
        call dgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, &
          zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, &
          rfac0, nquad, wnear)
      endif

      return
      end
!
!
!
!
!
      subroutine getnearquad_bh2d_v2b_neu(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)
!
!  This subroutine generates the near field quadrature
!  for the v2b Neumann representation using bh2d_gdn.
!
      implicit none
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      integer *8, intent(in) :: ndtarg, ntarg
      real *8, intent(in) :: targs(ndtarg,ntarg)
      integer *8, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      real *8, intent(in) :: eps
      integer *8, intent(in) :: iquadtype, nnz
      integer *8, intent(in) :: row_ptr(ntarg+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: nquad
      real *8, intent(out) :: wnear(nquad)

      integer *8 ipv, ndd, ndi, ndz, i
      real *8 dpars
      integer *8 ipars
      complex *16 zpars

      procedure (), pointer :: fker
      external bh2d_gdn

      ndd = 0
      ndz = 0
      ndi = 0

      if (iquadtype.eq.1) then
        ipv = 0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wnear(i) = 0
        enddo
!$OMP END PARALLEL DO

        fker => bh2d_gdn
        call dgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, &
          zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, &
          rfac0, nquad, wnear)
      endif

      return
      end
!
!
!
!
!
      subroutine getnearquad_bh2d_gdn(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)
!
!  This subroutine generates the near field quadrature
!  for the normal derivative of the biharmonic Green's function (bh2d_gdn).
!
      implicit none
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      integer *8, intent(in) :: ndtarg, ntarg
      real *8, intent(in) :: targs(ndtarg,ntarg)
      integer *8, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      real *8, intent(in) :: eps
      integer *8, intent(in) :: iquadtype, nnz
      integer *8, intent(in) :: row_ptr(ntarg+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: nquad
      real *8, intent(out) :: wnear(nquad)

      integer *8 ipv, ndd, ndi, ndz, i
      real *8 dpars
      integer *8 ipars
      complex *16 zpars

      procedure (), pointer :: fker
      external bh2d_gdn

      ndd = 0
      ndz = 0
      ndi = 0

      if (iquadtype.eq.1) then
        ipv = 0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wnear(i) = 0
        enddo
!$OMP END PARALLEL DO

        fker => bh2d_gdn
        call dgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, &
          zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, &
          rfac0, nquad, wnear)
      endif

      return
      end