!
!  Near field quadrature for a radial kernel represented by a
!  precomputed piecewise polynomial expansion in the vpp format of
!  vpp.f90.
!
!
      subroutine getnearquad_radcheb(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, &
        uvs_targ, eps, iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, &
        nquad, ndd, dpars, ndi, ipars, wnear)
!
!  This subroutine generates the near field quadrature for the
!  representation
!
!  u = \int_{\Gamma} K(|x-y|) \sigma(y) dS_{y}
!
!  where K is given by the piecewise polynomial expansion stored in
!  (ndd, dpars, ndi, ipars).
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
!    - npts: integer
!        total number of discretization points on the boundary
!    - srccoefs: real *8 (9,npts)
!        basis expansion coefficients of xyz, dxyz/du,
!        and dxyz/dv on each patch
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the
!        discretization nodes on the surface
!    - ndtarg: integer
!        leading dimension of target array
!    - ntarg: integer
!        number of targets
!    - targs: real *8 (ndtarg,ntarg)
!        target information
!    - ipatch_id: integer(ntarg)
!        id of patch of target i, id = -1, if target is off-surface
!    - uvs_targ: real *8 (2,ntarg)
!        local uv coordinates on patch if on surface, otherwise unused
!    - eps: real *8
!        precision requested
!    - iquadtype: integer
!        quadrature type
!          * iquadtype = 1, use ggq for self + adaptive integration
!            for rest
!    - nnz: integer
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer(ntarg+1)
!        row_ptr(i) is the pointer to col_ind array where list of
!        relevant source patches for target i start
!    - col_ind: integer (nnz)
!        list of source patches relevant for all targets, sorted
!        by the target number
!    - iquad: integer(nnz+1)
!        location in wnear array where quadrature for col_ind(i) starts
!    - rfac0: real *8
!        radius parameter for switching to predetermined quadrature rule
!    - nquad: integer
!        number of near field entries corresponding to each source
!        target pair
!    - ndd: integer
!        length of dpars
!    - dpars: real *8 (ndd)
!        panel endpoints and coefficients of the vpp expansion
!    - ndi: integer
!        length of ipars, must be 5
!    - ipars: integer *8 (ndi)
!        vpp pointer array, ipars(5) = 2 for a complex valued kernel
!
!  Output arguments
!    - wnear: complex *16(nquad)
!        The desired near field quadrature
!

!  Integer arguments are integer *8.

      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: ndtarg, ntarg
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches), ipatch_id(ntarg)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps, uvs_targ(2,ntarg), targs(ndtarg,ntarg)
      integer *8, intent(in) :: iquadtype, nnz
      integer *8, intent(in) :: row_ptr(ntarg+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: nquad, ndd, ndi
      real *8, intent(in) :: dpars(ndd)
      integer *8, intent(in) :: ipars(ndi)
      complex *16, intent(out) :: wnear(nquad)

      complex *16 zpars
      integer *8 ipv, i, ndz
      procedure (), pointer :: fker
      external vpp_kern


!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i = 1,nquad
        wnear(i) = 0
      enddo
!$OMP END PARALLEL DO

      if (iquadtype.eq.1) then
        ipv = 0
        ndz = 0
        zpars = 0

        fker => vpp_kern

        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, zpars, &
          ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
          wnear)
      endif

      return
      end subroutine getnearquad_radcheb
!
!
!
!
      subroutine getnearquad_radcheb_real(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, iquadtype, nnz, row_ptr, col_ind, &
        iquad, rfac0, nquad, ndd, dpars, ndi, ipars, wnear)
!
!  Real valued version of getnearquad_radcheb, for expansions with
!  ipars(5) = 1. See getnearquad_radcheb for the argument list.
!
!  Output arguments
!    - wnear: real *8(nquad)
!        The desired near field quadrature
!

!  Integer arguments are integer *8.

      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: ndtarg, ntarg
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches), ipatch_id(ntarg)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps, uvs_targ(2,ntarg), targs(ndtarg,ntarg)
      integer *8, intent(in) :: iquadtype, nnz
      integer *8, intent(in) :: row_ptr(ntarg+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: nquad, ndd, ndi
      real *8, intent(in) :: dpars(ndd)
      integer *8, intent(in) :: ipars(ndi)
      real *8, intent(out) :: wnear(nquad)

      complex *16 zpars
      integer *8 ipv, i, ndz
      procedure (), pointer :: fker
      external vpp_kern



!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i = 1,nquad
        wnear(i) = 0
      enddo
!$OMP END PARALLEL DO

      if (iquadtype.eq.1) then
        ipv = 0
        ndz = 0
        zpars = 0

        fker => vpp_kern

        call dgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, zpars, &
          ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
          wnear)
      endif

      return
      end subroutine getnearquad_radcheb_real
