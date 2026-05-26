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
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  User callable routines:
!    - stok_sprime_eval: Given sigma, evaluates S'[\sigma](x) at a
!        collection of targets x (on-surface or off-surface).
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  Advanced interfaces:
!    - getnearquad_stok_sprime: compute the near quadrature correction
!        for evaluating S'[\sigma] at target points which can be either
!        on-surface or off-surface, with user-provided near-information
!        prescribed in row-sparse compressed format.
!
!    - stok_sprime_eval_addsub: compute S'[\sigma] at a collection of
!        targets (on-surface or off-surface), given sigma.
!        On input, user provides precomputed near quadrature in
!        row-sparse compressed format and oversampling surface
!        information for the far-part.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getnearquad_stok_sprime(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, &
        targs, ipatch_id, uvs_targ, eps, iquadtype, nnz, &
        row_ptr, col_ind, iquad, rfac0, nquad, wnear)
!
!  This subroutine generates the near field quadrature for the
!  traction of the Stokes single layer:
!
!    S'[\sigma](x) = n_k(x) sigma_{ik}(S[\sigma])(x)
!
!  The kernel is t(S)_{ij} = T_{ijk} n_k(x) (target normal), which
!  has the same form as the Stokes double-layer kernel but with the
!  normal at the TARGET instead of the source.
!
!  This routine handles both on-surface and off-surface targets.
!  For on-surface targets, only the principal value part is returned;
!  the identity contribution (jump term +-sigma/2) must be added
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
      real *8, intent(out) :: wnear(6,nquad)

      real *8, allocatable :: wneartmp(:)
      integer *8 ndd, ndi, ndz
      real *8 dpars(2)
      complex *16 zpars(1)
      integer *8 ipars(2)

      integer *8 i, ipv, ii, ijloc(2,6)

      procedure (), pointer :: fker
      external st3d_strac

      ndd = 2
      ndi = 2
      ndz = 0
      fker => st3d_strac
      !
      ! The kernel st3d_strac has a principal-value (hypersingular) part
      ! on the surface; ipv=1 signals this to dgetnearquad_ggq_guru.
      !
      ipv = 1

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
!
!
      subroutine stok_sprime_eval(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, sigma, pot)
!
!f2py intent(in) npatches,norders,ixyzs,iptype,npts,srccoefs,srcvals
!f2py intent(in) ndtarg,ntarg,targs,ipatch_id,uvs_targ,eps
!f2py intent(in) sigma
!f2py intent(out) pot
!
!------------------------------
!  This subroutine evaluates the traction of the Stokes single layer:
!
!    pot = S'_{stok}[\sigma](x)
!        = n_k(x) sigma_{ik}(S[\sigma])(x)
!
!  where sigma_{ij}(u) = -p delta_{ij} + (du_i/dx_j + du_j/dx_i)
!  is the Stokes stress tensor and n(x) is the outward unit normal
!  at the target x.
!
!  For targets on the surface, this routine returns only the
!  principal-value part.  The caller is responsible for adding
!  the jump-term contribution +-sigma/2.
!
!  The computation uses the add-and-subtract FMM strategy:
!    * The far field is computed via stfmm3d (requesting pot+pre+grad).
!    * The traction is assembled from the FMM output as
!        pot_i = n_k(x)(-p delta_ik + du_i/dx_k + du_k/dx_i)
!    * Precomputed near-field corrections are added.
!    * The smooth near-field FMM contribution is subtracted.
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
!        basis expansion coefficients of xyz, dxyz/du, dxyz/dv
!    - srcvals: real *8 (12,npts)
!        xyz, tangent vectors, and normals at the discretization nodes
!    - ndtarg: integer *8
!        leading dimension of targs; must be >= 12
!    - ntarg: integer *8
!        number of targets
!    - targs: real *8(ndtarg,ntarg)
!        target information
!          * targs(1:3,i)   - xyz of target i
!          * targs(10:12,i) - unit outward normal at target i
!    - ipatch_id: integer *8(ntarg)
!        patch index for on-surface targets, -1 for off-surface
!    - uvs_targ: real *8(2,ntarg)
!        local uv coordinates (on-surface targets only)
!    - eps: real *8
!        precision requested
!    - sigma: real *8(3,npts)
!        layer density
!
!  Output arguments:
!    - pot: real *8(3,ntarg)
!        S'[\sigma] evaluated at the target points
!        (principal-value part only for on-surface targets)
!
!-----------------------------------

      implicit none
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: ndtarg, ntarg
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts), eps
      real *8, intent(in) :: targs(ndtarg,ntarg)
      real *8, intent(in) :: sigma(3,npts)
      integer *8, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)

      real *8, intent(out) :: pot(3,ntarg)

      ! local arrays
      real *8, allocatable :: wnear(:,:)
      integer *8, allocatable :: row_ptr(:), col_ind(:), iquad(:)
      integer *8, allocatable :: ixyzso(:), novers(:)
      real *8, allocatable :: srcover(:,:), wover(:)
      real *8, allocatable :: cms(:,:), rads(:), rad_near(:)

      integer *8 nnz, nquad, npts_over
      integer *8 iquadtype, nker, ikerorder
      real *8 rfac0, rfac

      integer *8 ndd, ndz, ndi, lwork, ndim
      real *8 dpars(1)
      complex *16 zpars
      integer *8 ipars(1)
      real *8 work(1)

      integer *8 i, iptype_avg, norder_avg

      !
      ! average patch type and order (same convention as stok_comb_vel_eval)
      !
      iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
      norder_avg = floor(sum(norders)/(npatches+0.0d0))

      call get_rfacs(norder_avg, iptype_avg, rfac, rfac0)

      !
      ! Step 1: centroid-based near-field detection
      !
      allocate(cms(3,npatches), rads(npatches), rad_near(npatches))

      call get_centroid_rads(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, cms, rads)

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i = 1,npatches
        rad_near(i) = rads(i)*rfac
      enddo
!$OMP END PARALLEL DO

      call findnearmem(cms, npatches, rad_near, ndtarg, targs, ntarg, &
        nnz)

      allocate(row_ptr(ntarg+1), col_ind(nnz))

      call findnear(cms, npatches, rad_near, ndtarg, targs, ntarg, &
        row_ptr, col_ind)

      allocate(iquad(nnz+1))
      call get_iquad_rsc(npatches, ixyzs, ntarg, nnz, row_ptr, &
        col_ind, iquad)

      !
      ! S' has the same singular order as the DLP (ikerorder = 0)
      !
      ikerorder = 0

      !
      ! Step 2: oversampled geometry for far-field FMM
      !
      allocate(novers(npatches), ixyzso(npatches+1))

      zpars = 0
      call get_far_order(eps, npatches, norders, ixyzs, iptype, cms, &
        rads, npts, srccoefs, ndtarg, ntarg, targs, ikerorder, zpars, &
        nnz, row_ptr, col_ind, rfac, novers, ixyzso)

      npts_over = ixyzso(npatches+1) - 1

      allocate(srcover(12,npts_over), wover(npts_over))

      call oversample_geom(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, srcvals, novers, ixyzso, npts_over, srcover)

      call get_qwts(npatches, novers, ixyzso, iptype, npts_over, &
        srcover, wover)

      !
      ! Step 3: near-field quadrature correction
      !
      nquad = iquad(nnz+1) - 1
      nker  = 6
      allocate(wnear(nker, nquad))

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i = 1,nquad
        wnear(1,i) = 0
        wnear(2,i) = 0
        wnear(3,i) = 0
        wnear(4,i) = 0
        wnear(5,i) = 0
        wnear(6,i) = 0
      enddo
!$OMP END PARALLEL DO

      iquadtype = 1

      call getnearquad_stok_sprime(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, iquadtype, nnz, &
        row_ptr, col_ind, iquad, rfac0, nquad, wnear)

      !
      ! Step 4: call the add-and-subtract evaluator
      !
      ndd   = 1
      ndz   = 0
      ndi   = 1
      lwork = 0
      ndim  = 3

      call stok_sprime_eval_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        eps, ndd, dpars, ndz, zpars, ndi, ipars, &
        nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, &
        novers, npts_over, ixyzso, srcover, wover, lwork, work, &
        ndim, sigma, pot)

      return
      end

!
!
!
!
!
      subroutine stok_sprime_eval_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        eps, ndd, dpars, ndz, zpars, ndi, ipars, &
        nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, &
        novers, nptso, ixyzso, srcover, whtsover, lwork, work, &
        ndim, sigma, pot)
!
!  This subroutine evaluates the traction of the Stokes single layer:
!
!    pot = S'_{stok}[\sigma](x)
!
!  where the near field is precomputed and stored in the row-sparse
!  compressed format, and the far field is handled by stfmm3d.
!
!  The FMM is called requesting velocity + pressure + gradient
!  (ifppregtarg = 3), and the traction at each target x is assembled as
!
!    pot_i(x) = n_k(x) * (-p(x) delta_ik + du_i/dx_k + du_k/dx_i)
!
!  The smooth near-field contribution is subtracted by calling
!  st3ddirectstokstrsg with istress = 0 (stoklet sources only),
!  assembling the traction the same way, then subtracting the result.
!
!  Input arguments:
!    - npatches: integer *8
!        number of patches
!    - norders: integer *8(npatches)
!        order of discretization on each patch
!    - ixyzs: integer *8(npatches+1)
!        starting location in srccoefs/srcvals for patch i
!    - iptype: integer *8(npatches)
!        type of patch
!    - npts: integer *8
!        total number of discretization points on the boundary
!    - srccoefs: real *8(9,npts)
!        basis expansion coefficients of xyz, dxyz/du, dxyz/dv
!    - srcvals: real *8(12,npts)
!        xyz, tangents, and normals at discretization nodes
!    - ndtarg: integer *8
!        leading dimension of targs; must be >= 12
!    - ntarg: integer *8
!        number of targets
!    - targs: real *8(ndtarg,ntarg)
!        target information; targs(1:3,i) = xyz, targs(10:12,i) = n(x)
!    - eps: real *8
!        precision requested
!    - ndd: integer *8
!        number of real kernel parameters (unused, pass 1)
!    - dpars: real *8(ndd)
!        real kernel parameters (unused)
!    - ndz: integer *8
!        number of complex kernel parameters (unused, pass 0)
!    - zpars: complex *16(ndz)
!        complex kernel parameters (unused)
!    - ndi: integer *8
!        number of integer kernel parameters (unused, pass 1)
!    - ipars: integer *8(ndi)
!        integer kernel parameters (unused)
!    - nnz: integer *8
!        number of source-patch -> target interactions in near field
!    - row_ptr: integer *8(ntarg+1)
!        row_ptr(i) points into col_ind for target i's patches
!    - col_ind: integer *8(nnz)
!        source patches for each target
!    - iquad: integer *8(nnz+1)
!        offset into wnear for each source-patch/target pair
!    - nquad: integer *8
!        total number of near-field quadrature entries
!    - nker: integer *8
!        number of kernel components in wnear; must be 6
!    - wnear: real *8(nker,nquad)
!        precomputed near-field quadrature for S'.
!        wnear(1,:)-(1,1), wnear(2,:)-(1,2), wnear(3,:)-(1,3),
!        wnear(4,:)-(2,2), wnear(5,:)-(2,3), wnear(6,:)-(3,3)
!    - novers: integer *8(npatches)
!        oversampling orders for each patch
!    - nptso: integer *8
!        total number of oversampled points
!    - ixyzso: integer *8(npatches+1)
!        starting location in srcover for patch i
!    - srcover: real *8(12,nptso)
!        oversampled surface information
!    - whtsover: real *8(nptso)
!        smooth quadrature weights at oversampled nodes
!    - lwork: integer *8
!        size of work array (unused; pass 0)
!    - work: real *8(lwork)
!        work array (unused)
!    - ndim: integer *8
!        number of densities per point; must be 3
!    - sigma: real *8(3,npts)
!        layer density
!
!  Output arguments:
!    - pot: real *8(3,ntarg)
!        S'[\sigma] evaluated at the targets
!        (principal-value part only for on-surface targets)
!

      implicit none
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)

      integer *8, intent(in) :: ndtarg, ntarg
      real *8, intent(in) :: targs(ndtarg,ntarg)

      real *8, intent(in) :: eps

      integer *8, intent(in) :: ndd, ndz, ndi
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      integer *8, intent(in) :: ipars(ndi)

      integer *8, intent(in) :: nnz, nquad
      integer *8, intent(in) :: row_ptr(ntarg+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)

      integer *8, intent(in) :: nker
      real *8, intent(in) :: wnear(nker,nquad)

      integer *8, intent(in) :: nptso
      integer *8, intent(in) :: ixyzso(npatches+1), novers(npatches)
      real *8, intent(in) :: srcover(12,nptso), whtsover(nptso)

      integer *8, intent(in) :: lwork
      real *8, intent(in) :: work(max(1,lwork))

      integer *8, intent(in) :: ndim
      real *8, intent(in) :: sigma(ndim,npts)

      real *8, intent(out) :: pot(ndim,ntarg)

      ! Local variables
      integer *8 ns, nt, ndsigma
      real *8, allocatable :: sources(:,:), targvals(:,:)
      real *8, allocatable :: sigmaover(:,:)
      real *8, allocatable :: stoklet(:,:)
      real *8, allocatable :: strslet(:,:), strsvec(:,:)
      real *8, allocatable :: rotlet(:,:), rotvec(:,:)
      real *8, allocatable :: doublet(:,:), doubvec(:,:)
      real *8, allocatable :: pottarg(:,:), pretarg(:), gradtarg(:,:,:)

      real *8, allocatable :: srctmp2(:,:), sttmp2(:,:)
      real *8 strstmp2(3), strsvec2(3)

      real *8 potv(3), pv, gradv(3,3)
      real *8 tmp(10)
      real *8 thresh, ra
      integer *8 nmax, nss

      real *8 w11, w12, w13, w21, w22, w23, w31, w32, w33
      real *8 sig1, sig2, sig3
      real *8 tnorm(3)

      integer *8 i, j, l, jpatch, jquadstart, jstart, npols, npover
      integer *8 nd, ntarg0
      integer *8 ier, iper
      integer *8 ifstoklet, ifstrslet, ifrotlet, ifdoublet
      integer *8 ifppreg, ifppregtarg
      integer *8 istress

      real *8 t1, t2
      real *8 done, pi

      parameter (nd=1, ntarg0=1)

      done = 1
      pi = atan(done)*4

      ns = nptso
      nt = ntarg

      !
      ! Estimate the max number of oversampled sources near any target
      !
      nmax = 0
      call get_near_corr_max(nt, row_ptr, nnz, col_ind, npatches, &
        ixyzso, nmax)
      allocate(srctmp2(3,nmax), sttmp2(3,nmax))

      allocate(sources(3,ns), targvals(3,nt))
      allocate(stoklet(3,ns))
      allocate(strslet(3,ns), strsvec(3,ns))
      allocate(rotlet(3,ns), rotvec(3,ns))
      allocate(doublet(3,ns), doubvec(3,ns))
      allocate(sigmaover(3,ns))
      allocate(pottarg(3,nt), pretarg(nt), gradtarg(3,3,nt))

      !
      ! Oversample the density
      !
      ndsigma = 3
      call oversample_fun_surf(ndsigma, npatches, norders, ixyzs, &
        iptype, npts, sigma, novers, ixyzso, ns, sigmaover)

      !
      ! Pack sources and stoklet strengths (weighted density)
      !
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i = 1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)

        stoklet(1,i) = sigmaover(1,i)*whtsover(i)
        stoklet(2,i) = sigmaover(2,i)*whtsover(i)
        stoklet(3,i) = sigmaover(3,i)*whtsover(i)

        strslet(1:3,i) = 0
        strsvec(1:3,i) = 0
        rotlet(1:3,i)  = 0
        rotvec(1:3,i)  = 0
        doublet(1:3,i) = 0
        doubvec(1:3,i) = 0
      enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i = 1,nt
        targvals(1,i) = targs(1,i)
        targvals(2,i) = targs(2,i)
        targvals(3,i) = targs(3,i)
      enddo
!$OMP END PARALLEL DO

      ifstoklet  = 1
      ifstrslet  = 0
      ifrotlet   = 0
      ifdoublet  = 0
      ifppreg    = 0
      ! Request velocity + pressure + gradient at targets
      ifppregtarg = 3

      iper = 0
      ier  = 0

      !
      ! Call the Stokes FMM for velocity, pressure, and velocity gradient
      !
      call cpu_time(t1)
      call stfmm3d(nd, eps, ns, sources, ifstoklet, stoklet, &
        ifstrslet, strslet, strsvec, ifrotlet, rotlet, rotvec, &
        ifdoublet, doublet, doubvec, ifppreg, tmp, tmp, tmp, nt, &
        targvals, ifppregtarg, pottarg, pretarg, gradtarg, ier)
      call cpu_time(t2)

      !
      ! Assemble traction from FMM output:
      !   pot_i = n_k * (-p delta_ik + du_i/dx_k + du_k/dx_i)
      !
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, tnorm)
      do i = 1,nt
        tnorm(1) = targs(10,i)
        tnorm(2) = targs(11,i)
        tnorm(3) = targs(12,i)

        pot(1,i) = tnorm(1)*(-pretarg(i) + 2*gradtarg(1,1,i)) + &
                   tnorm(2)*(gradtarg(1,2,i) + gradtarg(2,1,i)) + &
                   tnorm(3)*(gradtarg(1,3,i) + gradtarg(3,1,i))
        pot(2,i) = tnorm(1)*(gradtarg(1,2,i) + gradtarg(2,1,i)) + &
                   tnorm(2)*(-pretarg(i) + 2*gradtarg(2,2,i)) + &
                   tnorm(3)*(gradtarg(2,3,i) + gradtarg(3,2,i))
        pot(3,i) = tnorm(1)*(gradtarg(1,3,i) + gradtarg(3,1,i)) + &
                   tnorm(2)*(gradtarg(2,3,i) + gradtarg(3,2,i)) + &
                   tnorm(3)*(-pretarg(i) + 2*gradtarg(3,3,i))
      enddo
!$OMP END PARALLEL DO

      !
      ! Compute near-field FMM exclusion threshold
      !
      call get_fmm_thresh(3, ns, sources, 3, nt, targvals, thresh)

      !
      ! Add precomputed near-field quadrature correction
      !
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart) &
!$OMP PRIVATE(jstart,npols,l,sig1,sig2,sig3,w11,w12,w13,w21,w22,w23,w31,w32,w33)
      do i = 1,nt
        do j = row_ptr(i), row_ptr(i+1)-1
          jpatch    = col_ind(j)
          npols     = ixyzs(jpatch+1) - ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart    = ixyzs(jpatch)
          do l = 1,npols
            sig1 = sigma(1, jstart+l-1)
            sig2 = sigma(2, jstart+l-1)
            sig3 = sigma(3, jstart+l-1)

            w11 = wnear(1, jquadstart+l-1)
            w12 = wnear(2, jquadstart+l-1)
            w13 = wnear(3, jquadstart+l-1)
            w21 = w12
            w22 = wnear(4, jquadstart+l-1)
            w23 = wnear(5, jquadstart+l-1)
            w31 = w13
            w32 = w23
            w33 = wnear(6, jquadstart+l-1)

            pot(1,i) = pot(1,i) + w11*sig1 + w12*sig2 + w13*sig3
            pot(2,i) = pot(2,i) + w21*sig1 + w22*sig2 + w23*sig3
            pot(3,i) = pot(3,i) + w31*sig1 + w32*sig2 + w33*sig3
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

      !
      ! Subtract the smooth (oversampled) near-field FMM contribution
      !
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2,sttmp2) &
!$OMP PRIVATE(strstmp2,strsvec2,nss,l,jstart,npover) &
!$OMP PRIVATE(potv,pv,gradv,istress,tnorm)
      do i = 1,nt
        nss = 0
        do j = row_ptr(i), row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l = ixyzso(jpatch), ixyzso(jpatch+1)-1
            nss = nss + 1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)

            sttmp2(1,nss) = stoklet(1,l)
            sttmp2(2,nss) = stoklet(2,l)
            sttmp2(3,nss) = stoklet(3,l)
          enddo
        enddo

        potv(1:3)    = 0
        pv           = 0
        gradv(1:3,1:3) = 0

        istress = 0
        call st3ddirectstokstrsg(nd, srctmp2, sttmp2, istress, &
          strstmp2, strsvec2, nss, targvals(1,i), ntarg0, potv, pv, &
          gradv, thresh)

        ! Assemble traction from the direct evaluation output
        tnorm(1) = targs(10,i)
        tnorm(2) = targs(11,i)
        tnorm(3) = targs(12,i)

        potv(1) = tnorm(1)*(-pv + 2*gradv(1,1)) + &
                  tnorm(2)*(gradv(1,2) + gradv(2,1)) + &
                  tnorm(3)*(gradv(1,3) + gradv(3,1))
        potv(2) = tnorm(1)*(gradv(1,2) + gradv(2,1)) + &
                  tnorm(2)*(-pv + 2*gradv(2,2)) + &
                  tnorm(3)*(gradv(2,3) + gradv(3,2))
        potv(3) = tnorm(1)*(gradv(1,3) + gradv(3,1)) + &
                  tnorm(2)*(gradv(2,3) + gradv(3,2)) + &
                  tnorm(3)*(-pv + 2*gradv(3,3))

        pot(1,i) = pot(1,i) - potv(1)
        pot(2,i) = pot(2,i) - potv(2)
        pot(3,i) = pot(3,i) - potv(3)

      enddo
!$OMP END PARALLEL DO

      !
      ! Negate to match the double-negative sign convention used by the
      ! Stokes DLP evaluator: kern_SP = -t(S), so that SP_mat = -D^T
      ! (with swapped src/targ roles) in the same way that kern_D = -D.
      !
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i = 1,nt
        pot(1,i) = -pot(1,i)
        pot(2,i) = -pot(2,i)
        pot(3,i) = -pot(3,i)
      enddo
!$OMP END PARALLEL DO

      return
      end
