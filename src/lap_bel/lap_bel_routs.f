c
c     Laplace-Beltrami and Helmholtz-Beltrami parametrix/remainder
c     kernels. This file contains the following user callable
c     routines:
c
c       getnearquad_lap_bel_log - near field quadrature for the
c         Laplace-Beltrami parametrix log|x-y|^2/(4 pi) and for the
c         associated remainder
c
c       getnearquad_helm_bel_hank - near field quadrature for the
c         Helmholtz-Beltrami parametrix H_0(zk|x-y|)/(4i) and for the
c         associated remainder
c
c     Both routines require the mean curvature of the target surface to
c     be supplied in the 13th component of the target array.
c
c
      subroutine getnearquad_lap_bel_log(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2   ipatch_id,uvs_targ,eps,iktype,iquadtype,nnz,row_ptr,col_ind,
     3   iquad,rfac0,nquad,wnear)
c
c
c  This subroutine generates the near field quadrature for the
c  Laplace-Beltrami parametrix
c
c        K(x,y) = log|x-y|^2/(4 pi)
c
c  or for the associated remainder
c
c        R(x,y) = 4 (n(x).(x-y)/|x-y|^2)
c                   (n(x).(x-y)/|x-y|^2 - H(x))/(4 pi)
c
c  where H(x) is the mean curvature at the target x.
c
c  The quadrature is computed by the following strategy
c  targets within a sphere of radius rfac0*rs
c  of a patch centroid is handled using adaptive integration
c  where rs is the radius of the bounding sphere
c  for the patch
c
c  All other targets in the near field are handled via
c  oversampled quadrature
c
c  The recommended parameter for rfac0 is 1.25d0
c
c  Input arguments:
c
c    - npatches: integer
c        number of patches
c    - norders: integer(npatches)
c        order of discretization on each patch
c    - ixyzs: integer(npatches+1)
c        ixyzs(i) denotes the starting location in srccoefs,
c        and srcvals array corresponding to patch i
c    - iptype: integer(npatches)
c        type of patch
c        iptype = 1, triangular patch discretized using RV nodes
c    - npts: integer
c        total number of discretization points on the boundary
c    - srccoefs: real *8 (9,npts)
c        koornwinder expansion coefficients of xyz, dxyz/du,
c        and dxyz/dv on each patch.
c        For each point
c          * srccoefs(1:3,i) is xyz info
c          * srccoefs(4:6,i) is dxyz/du info
c          * srccoefs(7:9,i) is dxyz/dv info
c    - srcvals: real *8 (12,npts)
c        xyz(u,v) and derivative info sampled at the
c        discretization nodes on the surface
c          * srcvals(1:3,i) - xyz info
c          * srcvals(4:6,i) - dxyz/du info
c          * srcvals(7:9,i) - dxyz/dv info
c          * srcvals(10:12,i) - normals info
c    - ndtarg: integer
c        leading dimension of target array, must be at least 13
c    - ntarg: integer
c        number of targets
c    - targs: real *8 (ndtarg,ntarg)
c        target information, must contain the surface geometry in
c        components 1:12 in the same layout as srcvals, and the mean
c        curvature of the target in component 13
c    - ipatch_id: integer(ntarg)
c        id of patch of target i, id = -1, if target is off-surface
c    - uvs_targ: real *8 (2,ntarg)
c        local uv coordinates on patch if target is on surface
c    - eps: real *8
c        precision requested
c    - iktype: integer
c        kernel selector
c          * iktype = 1, remainder kernel R
c          * iktype = 2, parametrix kernel K
c    - iquadtype: integer
c        quadrature type
c          * iquadtype = 1, use ggq for self + adaptive integration
c            for rest
c    - nnz: integer
c        number of source patch-> target interactions in the near field
c    - row_ptr: integer(ntarg+1)
c        row_ptr(i) is the pointer
c        to col_ind array where list of relevant source patches
c        for target i start
c    - col_ind: integer (nnz)
c        list of source patches relevant for all targets, sorted
c        by the target number
c    - iquad: integer(nnz+1)
c        location in wnear array where quadrature for col_ind(i)
c        starts
c    - rfac0: real *8
c        radius parameter for near field
c    - nquad: integer
c        number of near field entries corresponding to each source
c        target pair
c
c  Output arguments
c    - wnear: real *8(nquad)
c        The desired near field quadrature
c
c

      implicit none
      integer *8, intent(in) :: npatches,norders(npatches),npts,nquad
      integer *8, intent(in) :: ixyzs(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      integer *8, intent(in) :: ndtarg,ntarg
      real *8, intent(in) :: targs(ndtarg,ntarg)
      integer *8, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: iktype
      integer *8, intent(in) :: iquadtype
      integer *8, intent(in) :: nnz
      integer *8, intent(in) :: row_ptr(ntarg+1),col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(out) :: wnear(nquad)

      integer *8 ipars(1)
      integer *8 ndd,ndz,ndi
      real *8 dpars
      complex *16 zpars

      integer *8 ipv

      procedure (), pointer :: fker
      external lap_bel_log,lap_bel_res

      if(ndtarg.lt.13) then
        call prinf('ndtarg must be at least 13, ndtarg=*',ndtarg,1)
        stop
      endif

c
c
c        initialize the appropriate kernel function
c

      ndd = 0
      ndi = 0
      ndz = 1
      zpars = 0
      if(iquadtype.eq.1) then
        if (iktype.eq.1) then
          ipv = 0
          fker=>lap_bel_res
        elseif (iktype.eq.2) then
          ipv = 0
          fker=>lap_bel_log
        endif
        call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2     ipatch_id,uvs_targ,
     3     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     4     col_ind,iquad,
     5     rfac0,nquad,wnear)
      endif


      return
      end
c
c
c
c
c

      subroutine getnearquad_helm_bel_hank(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2   ipatch_id,uvs_targ,eps,zpars,iktype,iquadtype,nnz,row_ptr,
     3   col_ind,iquad,rfac0,nquad,wnear)
c
c
c  This subroutine generates the near field quadrature for the
c  Helmholtz-Beltrami parametrix
c
c        K(x,y) = H_0(zk|x-y|)/(4i)
c
c  or for the associated remainder
c
c        R(x,y) = (-zk^2 H_2(zk|x-y|) (n(x).(x-y)/|x-y|)^2
c                   + 2 H(x) zk H_1(zk|x-y|) n(x).(x-y)/|x-y|)/(4i)
c
c  where H(x) is the mean curvature at the target x.
c
c  The quadrature is computed by the following strategy
c  targets within a sphere of radius rfac0*rs
c  of a patch centroid is handled using adaptive integration
c  where rs is the radius of the bounding sphere
c  for the patch
c
c  All other targets in the near field are handled via
c  oversampled quadrature
c
c  The recommended parameter for rfac0 is 1.25d0
c
c  Input arguments:
c
c    - npatches: integer
c        number of patches
c    - norders: integer(npatches)
c        order of discretization on each patch
c    - ixyzs: integer(npatches+1)
c        ixyzs(i) denotes the starting location in srccoefs,
c        and srcvals array corresponding to patch i
c    - iptype: integer(npatches)
c        type of patch
c        iptype = 1, triangular patch discretized using RV nodes
c    - npts: integer
c        total number of discretization points on the boundary
c    - srccoefs: real *8 (9,npts)
c        koornwinder expansion coefficients of xyz, dxyz/du,
c        and dxyz/dv on each patch.
c        For each point
c          * srccoefs(1:3,i) is xyz info
c          * srccoefs(4:6,i) is dxyz/du info
c          * srccoefs(7:9,i) is dxyz/dv info
c    - srcvals: real *8 (12,npts)
c        xyz(u,v) and derivative info sampled at the
c        discretization nodes on the surface
c          * srcvals(1:3,i) - xyz info
c          * srcvals(4:6,i) - dxyz/du info
c          * srcvals(7:9,i) - dxyz/dv info
c          * srcvals(10:12,i) - normals info
c    - ndtarg: integer
c        leading dimension of target array, must be at least 13
c    - ntarg: integer
c        number of targets
c    - targs: real *8 (ndtarg,ntarg)
c        target information, must contain the surface geometry in
c        components 1:12 in the same layout as srcvals, and the mean
c        curvature of the target in component 13
c    - ipatch_id: integer(ntarg)
c        id of patch of target i, id = -1, if target is off-surface
c    - uvs_targ: real *8 (2,ntarg)
c        local uv coordinates on patch if target is on surface
c    - eps: real *8
c        precision requested
c    - zpars: complex *16
c        Helmholtz wave number zk
c    - iktype: integer
c        kernel selector
c          * iktype = 1, remainder kernel R
c          * iktype = 2, parametrix kernel K
c    - iquadtype: integer
c        quadrature type
c          * iquadtype = 1, use ggq for self + adaptive integration
c            for rest
c    - nnz: integer
c        number of source patch-> target interactions in the near field
c    - row_ptr: integer(ntarg+1)
c        row_ptr(i) is the pointer
c        to col_ind array where list of relevant source patches
c        for target i start
c    - col_ind: integer (nnz)
c        list of source patches relevant for all targets, sorted
c        by the target number
c    - iquad: integer(nnz+1)
c        location in wnear array where quadrature for col_ind(i)
c        starts
c    - rfac0: real *8
c        radius parameter for near field
c    - nquad: integer
c        number of near field entries corresponding to each source
c        target pair
c
c  Output arguments
c    - wnear: complex *16(nquad)
c        The desired near field quadrature
c
c

      implicit none
      integer *8, intent(in) :: npatches,norders(npatches),npts,nquad
      integer *8, intent(in) :: ixyzs(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      integer *8, intent(in) :: ndtarg,ntarg
      real *8, intent(in) :: targs(ndtarg,ntarg)
      integer *8, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      real *8, intent(in) :: rfac0
      complex *16, intent(in) :: zpars
      integer *8, intent(in) :: iktype
      integer *8, intent(in) :: iquadtype
      integer *8, intent(in) :: nnz
      integer *8, intent(in) :: row_ptr(ntarg+1),col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      complex *16, intent(out) :: wnear(nquad)

      integer *8 ipars(1)
      integer *8 ndd,ndz,ndi
      real *8 dpars

      integer *8 ipv

      procedure (), pointer :: fker
      external helm_bel_hank,helm_bel_res

      if(ndtarg.lt.13) then
        call prinf('ndtarg must be at least 13, ndtarg=*',ndtarg,1)
        stop
      endif

c
c
c        initialize the appropriate kernel function
c

      ndd = 0
      ndi = 0
      ndz = 1
      if(iquadtype.eq.1) then
        if (iktype.eq.1) then
          ipv = 0
          fker=>helm_bel_res
        elseif (iktype.eq.2) then
          ipv = 0
          fker=>helm_bel_hank
        endif
        call zgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2     ipatch_id,uvs_targ,
     3     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     4     col_ind,iquad,
     5     rfac0,nquad,wnear)
      endif


      return
      end
c
c
c
c
c
