!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  Evaluator for the traction of the Stokes single layer (S' kernel)
!
!  The S' operator (also written T[S] or st3d_strac) is the target-
!  traction of the single-layer potential:
!
!    (S'[\sigma])_i(x) = n_k(x) \sigma_{ik}(u)(x)
!                      = n_k(x)(-p(x)\delta_{ik}
!                              + \partial_i u_k + \partial_k u_i)(x)
!
!  where  u = S_{stok}[\sigma]  and  x  carries the target normal n(x).
!
!  The kernel is t(S)_{ij}(x,y) = T_{ijk}(x,y) n_k(x):
!
!    t(S)_{ij}(x,y) = -3(r . n(x)) r_i r_j / r^5 / (4 pi)
!
!  with r = x - y.  This is the negative transpose of the Stokes
!  double-layer kernel  D_{ij}(x,y) = T_{ijk}(x,y) n_k(y).
!
!  Because t(S)_{ij} = t(S)_{ji}, only the 6 lower-triangular entries
!  of the 3x3 tensor are stored in wnear.
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  Advanced interfaces:
!    - getnearquad_stok_combtrac: compute the near quadrature correction
!        for evaluating S'[\sigma] at target points which can be either
!        on-surface or off-surface, with user-provided near-information
!        prescribed in row-sparse compressed format.
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getnearquad_stok_combtrac(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, &
        targs, ipatch_id, uvs_targ, eps, dpars, iquadtype, nnz, &
        row_ptr, col_ind, iquad, rfac0, nquad, wnear)
!
!  This subroutine generates the near field quadrature for the
!  traction of the Stokes combined field layer:
!
!    S'[\sigma](x) = dpars(1) * n_k(x) sigma_{ik}(S[\sigma])(x) 
!                      + dpars(2) * n_k(x) sigma_{ik}(T[\sigma])(x)
!
!
!  This routine handles both on-surface and off-surface targets.
!  For on-surface targets, only the principal value part is returned;
!  the identity contribution (jump term +-dpars(1)*sigma/2) must be added
!  separately by the caller.
!
!  The quadrature strategy is:
!    * targets within rfac0 * (patch bounding sphere radius) of the
!      patch centroid => adaptive integration (GGQ)
!    * all other near-field targets => oversampled quadrature
!
!  Recommended value: rfac0 = 1.25d0
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
!        iptype = 12, quadrangular patch discretized with Chebyshev nodes
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
!        leading dimension of the target information array.
!        Must be >= 12 so that target normals (entries 10:12) are
!        accessible for on-surface targets and provided off-surface.
!    - ntarg: integer *8
!        number of targets
!    - targs: real *8(ndtarg,ntarg)
!        target information
!          * targs(1:3,i)   - xyz coordinates of target i
!          * targs(10:12,i) - outward unit normal at target i
!    - ipatch_id: integer *8(ntarg)
!        ipatch_id(i) is the patch index if target i lies on the
!        surface, and -1 (or 0) if the target is off-surface
!    - uvs_targ: real *8(2,ntarg)
!        local uv coordinates on patch ipatch_id(i) if on-surface,
!        otherwise set to 0
!    - eps: real *8
!        precision requested
!    - iquadtype: integer *8
!        quadrature type
!          * iquadtype = 1, use GGQ for self + adaptive integration
!            for near-field neighbors
!    - nnz: integer *8
!        number of source-patch -> target interactions in the near field
!    - row_ptr: integer *8(ntarg+1)
!        row_ptr(i) is the pointer into col_ind where the list of
!        relevant source patches for target i starts
!    - col_ind: integer *8(nnz)
!        list of source patches relevant for all targets, sorted
!        by target number
!    - iquad: integer *8(nnz+1)
!        location in wnear where the quadrature for col_ind(i) starts
!        for a single kernel entry. Entry m of the tensor is stored at
!        offset (m-1)*nquad + iquad(i)
!    - rfac0: real *8
!        radius parameter for switching to oversampled quadrature
!    - nquad: integer *8
!        number of near-field quadrature entries per kernel component
!
!  Output arguments:
!    - wnear: real *8(6,nquad)
!        Near-field quadrature weights for S'. Since the tensor is
!        symmetric only the 6 lower-triangular entries are stored:
!        * wnear(1,:) - (1,1) entry
!        * wnear(2,:) - (1,2) entry
!        * wnear(3,:) - (1,3) entry
!        * wnear(4,:) - (2,2) entry
!        * wnear(5,:) - (2,3) entry
!        * wnear(6,:) - (3,3) entry
!

      implicit none
      integer *8, intent(in) :: npatches, norders(npatches), npts, nquad
      integer *8, intent(in) :: ixyzs(npatches+1), iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts), eps
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: iquadtype
      integer *8, intent(in) :: nnz
      integer *8, intent(in) :: ndtarg, ntarg
      real *8, intent(in) :: targs(ndtarg, ntarg)
      integer *8, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      integer *8, intent(in) :: row_ptr(ntarg+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: dpars(2)
      real *8, intent(out) :: wnear(6,nquad)

      real *8, allocatable :: wneartmp(:)
      integer *8 ndd, ndi, ndz
      
      complex *16 zpars(1)
      integer *8 ipars(2)

      integer *8 i, ipv, ii, ijloc(2,6)

      procedure (), pointer :: fker
      external st3d_combtrac

      ndd = 2
      ndi = 2
      ndz = 0
      fker => st3d_combtrac
      ipv = 2

      ! (i,j) pairs for the 6 independent entries of the symmetric tensor
      ijloc(1,1) = 1 ;  ijloc(2,1) = 1
      ijloc(1,2) = 1 ;  ijloc(2,2) = 2
      ijloc(1,3) = 1 ;  ijloc(2,3) = 3
      ijloc(1,4) = 2 ;  ijloc(2,4) = 2
      ijloc(1,5) = 2 ;  ijloc(2,5) = 3
      ijloc(1,6) = 3 ;  ijloc(2,6) = 3

      allocate(wneartmp(nquad))

      if (iquadtype .eq. 1) then

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii)
        do ii = 1,nquad
          wneartmp(ii) = 0
        enddo
!$OMP END PARALLEL DO

        do i = 1,6
          call dgetnearquad_ggq_guru(npatches, norders, ixyzs, &
            iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
            ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, &
            zpars, ndi, ijloc(1,i), nnz, row_ptr, col_ind, iquad, &
            rfac0, nquad, wneartmp)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii)
          do ii = 1,nquad
            wnear(i,ii) = wneartmp(ii)
          enddo
!$OMP END PARALLEL DO
        enddo

      endif

      return
      end

!
!
!