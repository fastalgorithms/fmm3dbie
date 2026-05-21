!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!   Near quadrature wrappers for flexural-gravity waves in 2d
!
!  PDE:
!    (\Delta^2 - zk1^2)(\Delta^2 - zk2^2) u = 0
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  Subroutines defined in this file:
!    - getnearquad_flex2d_gv2v: near quadrature for the flexural Green's
!        function on-surface (volume-to-volume)
!
!    - getnearquad_flex2d_dir: near quadrature for the flexural Green's
!        function, Dirichlet targets
!
!    - getnearquad_flex2d_neu: near quadrature for the normal derivative
!        of the flexural Green's function, Neumann targets
!
!    - getnearquad_flex2d_supp2: near quadrature for the simply-supported
!        plate boundary operator (2nd-order term)
!
!    - getnearquad_flex2d_free2: near quadrature for the free plate
!        boundary operator (2nd-order term), constant coefficients
!
!    - getnearquad_flex2d_free2_var: near quadrature for the free plate
!        boundary operator (2nd-order term), variable coefficients
!
!    - getnearquad_flex2d_var: near quadrature for a selected flexural
!        kernel, chosen via kernel_id
!

      subroutine getnearquad_flex2d_gv2v(npatches, norders, &
       ixyzs, iptype, npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, &
       row_ptr, col_ind, iquad, rfac0, nquad, wnear)
!
!  This subroutine generates the near field quadrature
!  for the representation:
!
!  u = G_{flex}[\sigma]
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
!    - zpars: complex *16 (2)
!        kernel parameters
!        zpars(1) = zk1, first flexural wavenumber (propagating)
!        zpars(2) = zk2, second flexural wavenumber (evanescent)
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
!    - wnear: complex *16(nquad)
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
      
      complex *16, intent(out) :: wnear(nquad)

      integer *8 ndtarg, ntarg, ndd, ndz, ndi
      real *8 dpars
      complex *16 zpars(2)
      integer *8 ipars(1)

      real *8, allocatable :: uvs_targ(:,:)
      integer *8, allocatable :: ipatch_id(:)

      integer *8 i
      integer *8 ipv

      procedure (), pointer :: fker
      external flex2d_g


      ndz=2
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
      fker => flex2d_g
      call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, srcvals, &
        ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear)

      return
      end subroutine getnearquad_flex2d_gv2v
!
!
!
!
!
      subroutine getnearquad_flex2d_dir(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)
!
!
!  This subroutine generates the near field quadrature
!  for the representation:
!
!  u = G_{flex}[\sigma]
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
!    - zpars: complex *16 (2)
!        kernel parameters
!        zpars(1) = zk1, first flexural wavenumber (propagating)
!        zpars(2) = zk2, second flexural wavenumber (evanescent)
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
!    - wnear: complex *16(nquad)
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
      complex *16, intent(out) :: wnear(nquad)

      integer *8 ipv, ndd, ndi, ndz, i
      real *8 dpars
      integer *8 ipars
      complex *16 zpars(2)

      procedure (), pointer :: fker
      external flex2d_g,flex2d_gdn

      ndz=2
      ndd=0
      ndi=0

      if (iquadtype.eq.1) then
        ipv = 0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wnear(i) = 0 
        enddo
!$OMP END PARALLEL DO

        fker => flex2d_g
!        fker => flex2d_gdn
        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
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
!
!
!
!
!
      subroutine getnearquad_flex2d_neu(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)
!
!
!  This subroutine generates the near field quadrature
!  for the representation:
!
!  u = \partial G_{flex}/\partial n_{y} [\sigma]
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
!    - zpars: complex *16 (2)
!        kernel parameters
!        zpars(1) = zk1, first flexural wavenumber (propagating)
!        zpars(2) = zk2, second flexural wavenumber (evanescent)
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
!    - wnear: complex *16(nquad)
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
      complex *16, intent(out) :: wnear(nquad)

      integer *8 ipv, ndd, ndi, ndz, i
      real *8 dpars
      integer *8 ipars
      complex *16 zpars(2)

      procedure (), pointer :: fker
      external flex2d_gdn

      ndz=2
      ndd=0
      ndi=0


    
      if (iquadtype.eq.1) then
        ipv = 0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wnear(i) = 0 
        enddo
!$OMP END PARALLEL DO
        fker => flex2d_gdn
        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
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
!
!
!
!
!
      subroutine getnearquad_flex2d_supp2(npatches, norders, &
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
!    - zpars: complex *16 (2)
!        kernel parameters
!        zpars(1) = zk1, first flexural wavenumber (propagating)
!        zpars(2) = zk2, second flexural wavenumber (evanescent)
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
!    - wnear: complex *16(nquad)
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
      complex *16, intent(out) :: wnear(nquad)

      integer *8 ipv, ndd, ndi, ndz, i
      real *8 dpars(1)
      integer *8 ipars
      complex *16 zpars(2)

      procedure (), pointer :: fker
      external flex2d_gsupp2

      ndd = 1
      ndz = 2
      ndi = 0



    
      if (iquadtype.eq.1) then
        ipv = 0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wnear(i) = 0 
        enddo
!$OMP END PARALLEL DO
        fker => flex2d_gsupp2
        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
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
!
!
!
      subroutine getnearquad_flex2d_free2(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, dpars, zpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)
!
!
!  This subroutine generates the near field quadrature
!  for the free plate boundary operator (2nd-order term), constant coefficients.
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
!    - zpars: complex *16 (2)
!        kernel parameters
!        zpars(1) = zk1, first flexural wavenumber (propagating)
!        zpars(2) = zk2, second flexural wavenumber (evanescent)
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
!    - wnear: complex *16(nquad)
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
      complex *16, intent(out) :: wnear(nquad)

      integer *8 ipv, ndd, ndi, ndz, i
      real *8 dpars(1)
      integer *8 ipars
      complex *16 zpars(2)

      procedure (), pointer :: fker
      external flex2d_gfree2

      ndd = 1
      ndz = 2
      ndi = 0



    
      if (iquadtype.eq.1) then
        ipv = 0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wnear(i) = 0 
        enddo
!$OMP END PARALLEL DO
        fker => flex2d_gfree2
        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
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
      subroutine getnearquad_flex2d_free2_var(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, dpars, zpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)
!
!
!  This subroutine generates the near field quadrature
!  for the free plate boundary operator (2nd-order term), variable coefficients.
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
!        dpars(1) = nu, Poisson ratio (spatially varying, evaluated at source)
!    - zpars: complex *16 (2)
!        kernel parameters
!        zpars(1) = zk1, first flexural wavenumber (propagating)
!        zpars(2) = zk2, second flexural wavenumber (evanescent)
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
!    - wnear: complex *16(nquad)
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
      complex *16, intent(out) :: wnear(nquad)

      integer *8 ipv, ndd, ndi, ndz, i
      real *8 dpars(1)
      integer *8 ipars
      complex *16 zpars(2)

      procedure (), pointer :: fker
      external flex2d_gfree2_var

      ndd = 1
      ndz = 2
      ndi = 0



    
      if (iquadtype.eq.1) then
        ipv = 0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wnear(i) = 0 
        enddo
!$OMP END PARALLEL DO
        fker => flex2d_gfree2_var
        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
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
!
!
!
      subroutine getnearquad_flex2d_var(npatches, & 
        norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, &
        ntarg, targs, ipatch_id, uvs_targ, eps, zpars, &
        iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
        kernel_id, wnear)
!
!
!  This subroutine generates the near field quadrature for a single
!  selected flexural kernel, chosen by kernel_id.
!
!  kernel_id = 1: flex2d_lapgx  (\partial_x \Delta G_{flex})
!  kernel_id = 2: flex2d_lapgy  (\partial_y \Delta G_{flex})
!  kernel_id = 3: flex2d_lapg   (\Delta G_{flex})
!  kernel_id = 4: flex2d_gyy    (\partial_{yy} G_{flex})
!  kernel_id = 5: flex2d_gxx    (\partial_{xx} G_{flex})
!  kernel_id = 6: flex2d_gxy    (\partial_{xy} G_{flex})
!  kernel_id = 7: flex2d_g      (G_{flex})
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
!    - zpars: complex *16 (2)
!        kernel parameters
!        zpars(1) = zk1, first flexural wavenumber (propagating)
!        zpars(2) = zk2, second flexural wavenumber (evanescent)
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
!    - kernel_id: integer *8
!        selects which flexural kernel to use (1--7, see above)
!
!  Output arguments
!    - wnear: complex *16(nquad)
!        The desired near field quadrature
!
!
!  The order of kernels below is based on the order appearing in the 
!  integral equation for flexural-gravity waves (see 3.15-3.21 in Askham, 
!  Hoskins, Nekrasov, Rachh, 2025).
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
      integer *8, intent(in) :: nquad, kernel_id
      complex *16, intent(out) :: wnear(nquad)

      integer *8 ipv, ndd, ndi, ndz, i
      real *8 dpars(1)
      integer *8 ipars
      complex *16 zpars(2)
  
      procedure (), pointer :: fker
      external flex2d_lapgx, flex2d_lapgy
      external flex2d_lapg, flex2d_gyy, flex2d_gxx
      external flex2d_gxy, flex2d_g

      ndd = 1
      ndz = 2
      ndi = 0


      if (kernel_id.eq.1) then
          fker => flex2d_lapgx
      else if (kernel_id.eq.2) then
          fker => flex2d_lapgy
      else if (kernel_id.eq.3) then 
          fker => flex2d_lapg 
      else if (kernel_id.eq.4) then
          fker => flex2d_gyy 
      else if (kernel_id.eq.5) then
          fker => flex2d_gxx 
      else if (kernel_id.eq.6) then
          fker => flex2d_gxy
      else if (kernel_id.eq.7) then
          fker => flex2d_g
      end if

    
      if (iquadtype.eq.1) then
        ipv = 0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wnear(i) = 0 
        enddo
!$OMP END PARALLEL DO
        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, &
          zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, &
          rfac0, nquad, wnear)

      endif 
      return
      end

