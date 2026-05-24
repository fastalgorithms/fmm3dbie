!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!   Combined prime representation for Laplace problems
!
!  PDE:
!    \Delta u = 0
!
!  This file provides quadrature and evaluation for the
!  combined prime kernel:
!
!    K[\sigma] = (\alpha S_{0}' + \beta D_{0}')[\sigma]
!
!  where
!    S_{0}' = d/dn_x G_0   (normal derivative at target)
!    D_{0}' = d/dn_x d/dn_y G_0  (hypersingular)
!
!  dpars(1) = alpha
!  dpars(2) = beta
!
!  User callable routines:
!    - getnearquad_lap_comb_cprime: compute near quadrature for
!        on-surface targets
!
!    - getnearquad_lap_comb_cprime_eval: compute near quadrature for
!        general (on or off surface) targets
!
!    - lap_comb_cprime_eval: evaluate the layer potential at a
!        collection of targets using the combined prime kernel
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      subroutine getnearquad_lap_comb_cprime(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, &
        eps, dpars, iquadtype, nnz, row_ptr, col_ind, &
        iquad, rfac0, nquad, wnear)
!
!  This subroutine generates the near field quadrature
!  for the combined prime representation
!    K = (\alpha S_{0}' + \beta D_{0}')[\sigma]
!  where the near field is specified by the user
!  in row sparse compressed format, and targets are on-surface.
!
!  The quadrature is computed by the following strategy:
!   targets within a sphere of radius rfac0*rs of a patch centroid
!   are handled using adaptive integration, where rs is the radius
!   of the bounding sphere for the patch.
!   All other targets in the near field are handled via
!   oversampled quadrature.
!
!  The recommended parameter for rfac0 is 1.25d0
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
!    - npts: integer *8
!        total number of discretization points
!    - srccoefs: real *8 (9,npts)
!        koornwinder expansion coefficients of xyz, dxyz/du, dxyz/dv
!    - srcvals: real *8 (12,npts)
!        xyz(u,v), dxyz/du, dxyz/dv, normals
!    - eps: real *8
!        precision requested
!    - dpars: real *8(2)
!        dpars(1) = alpha, dpars(2) = beta
!    - iquadtype: integer *8
!        quadrature type (1 = ggq)
!    - nnz: integer *8
!        number of source patch -> target interactions in near field
!    - row_ptr: integer *8(npts+1)
!    - col_ind: integer *8(nnz)
!    - iquad: integer *8(nnz+1)
!    - rfac0: real *8
!        radius parameter for near field
!    - nquad: integer *8
!        number of entries in wnear
!
!  Output arguments:
!    - wnear: real *8(nquad)
!

      implicit none
      integer *8, intent(in) :: npatches, norders(npatches), npts, nquad
      integer *8, intent(in) :: ixyzs(npatches+1), iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts), eps
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: iquadtype
      integer *8, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      real *8, intent(in) :: dpars(2)
      integer *8, intent(in) :: nnz
      integer *8, intent(in) :: row_ptr(npts+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(out) :: wnear(nquad)
      integer *8 i
      integer *8 ndtarg, ntarg

      ndtarg = 12
      ntarg = npts
      allocate(ipatch_id(npts), uvs_targ(2,npts))

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i = 1,npts
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo
!$OMP END PARALLEL DO

      call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts, &
        ipatch_id, uvs_targ)

      call getnearquad_lap_comb_cprime_eval(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, srcvals, &
        ipatch_id, uvs_targ, eps, dpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)

      return
      end subroutine getnearquad_lap_comb_cprime
!
!
!


      subroutine getnearquad_lap_comb_cprime_eval(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, dpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)
!
!  This subroutine generates the near field quadrature
!  for the combined prime representation
!    K = (\alpha S_{0}' + \beta D_{0}')[\sigma]
!  where the near field is specified by the user
!  in row sparse compressed format.
!  Targets can be on-surface or off-surface.
!
!  ipv = 2 is used for hypersingular quadrature.
!
!  Input arguments:
!    - npatches: integer *8
!    - norders: integer *8(npatches)
!    - ixyzs: integer *8(npatches+1)
!    - iptype: integer *8(npatches)
!    - npts: integer *8
!    - srccoefs: real *8 (9,npts)
!    - srcvals: real *8 (12,npts)
!    - ndtarg: integer *8
!        leading dimension of target array (must be >= 12)
!    - ntarg: integer *8
!        number of targets
!    - targs: real *8 (ndtarg,ntarg)
!        target information (must include normals in rows 10:12)
!    - ipatch_id: integer *8(ntarg)
!        patch id of target (-1 if off-surface)
!    - uvs_targ: real *8 (2,ntarg)
!    - eps: real *8
!    - dpars: real *8(2)
!        dpars(1) = alpha, dpars(2) = beta
!    - iquadtype: integer *8 (1 = ggq)
!    - nnz: integer *8
!    - row_ptr: integer *8(ntarg+1)
!    - col_ind: integer *8(nnz)
!    - iquad: integer *8(nnz+1)
!    - rfac0: real *8
!    - nquad: integer *8
!
!  Output arguments:
!    - wnear: real *8(nquad)
!

      implicit none
      integer *8, intent(in) :: npatches, norders(npatches), npts, nquad
      integer *8, intent(in) :: ixyzs(npatches+1), iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts), eps
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: ndtarg, ntarg
      integer *8, intent(in) :: iquadtype
      real *8, intent(in) :: targs(ndtarg,ntarg)
      integer *8, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      real *8, intent(in) :: dpars(2)
      integer *8, intent(in) :: nnz
      integer *8, intent(in) :: row_ptr(ntarg+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(out) :: wnear(nquad)

      integer *8 ipars(1)
      integer *8 ndd, ndz, ndi
      complex *16 zpars(1)
      integer *8 i, ipv

      procedure (), pointer :: fker
      external l3d_combprime

      ndd = 2
      ndz = 0
      ndi = 0
      ipv = 2
      fker => l3d_combprime

      if (iquadtype .eq. 1) then
        call dgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, &
          zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, &
          nquad, wnear)
      endif

      return
      end subroutine getnearquad_lap_comb_cprime_eval
!
!
!


      subroutine lap_comb_cprime_eval(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, dpars, sigma, pot)
!
!f2py intent(in) npatches,norders,ixyzs,iptype,npts,srccoefs,srcvals
!f2py intent(in) ndtarg,ntarg,targs,ipatch_id,uvs_targ,eps,dpars
!f2py intent(in) sigma
!f2py intent(out) pot
!
!------------------------------
!  This subroutine evaluates the layer potential for the representation
!
!    K[\sigma] = (\alpha S_{0}' + \beta D_{0}')[\sigma]
!
!  where
!    S_{0}' = d/dn_x G_0   (normal derivative at target x)
!    D_{0}' = d/dn_x d/dn_y G_0  (hypersingular)
!
!  The far field is computed via FMM using gradient output at targets,
!  then dotted with the target normal.
!  The near field uses the precomputed hypersingular quadrature (ipv=2).
!
!  Note: For targets on the boundary, this routine computes only the
!  principal value part.
!
!  Input arguments:
!    - npatches: integer *8
!    - norders: integer *8(npatches)
!    - ixyzs: integer *8(npatches+1)
!    - iptype: integer *8(npatches)
!    - npts: integer *8
!    - srccoefs: real *8(9,npts)
!    - srcvals: real *8(12,npts)
!    - ndtarg: integer *8
!        leading dimension of target array (must be >= 12)
!    - ntarg: integer *8
!    - targs: real *8(ndtarg,ntarg)
!        target info; rows 10:12 must contain target normals
!    - ipatch_id: integer *8(ntarg)
!        -1 if off-surface
!    - uvs_targ: real *8(2,ntarg)
!    - eps: real *8
!    - dpars: real *8(2)
!        dpars(1) = alpha, dpars(2) = beta
!    - sigma: real *8(npts)
!
!  Output arguments:
!    - pot: real *8(ntarg)
!
!-----------------------------------

      implicit none
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: ndtarg, ntarg
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts), eps
      real *8, intent(in) :: targs(ndtarg,ntarg)
      real *8, intent(in) :: dpars(2)
      real *8, intent(in) :: sigma(npts)
      integer *8, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      real *8, intent(out) :: pot(ntarg)

      integer *8 nptso, nnz, nquad

      integer *8 nover, npolso
      integer *8 norder, npols
      integer *8, allocatable :: row_ptr(:), col_ind(:), iquad(:)
      real *8, allocatable :: wnear(:)

      real *8, allocatable :: srcover(:,:), wover(:)
      integer *8, allocatable :: ixyzso(:), novers(:)

      real *8, allocatable :: cms(:,:), rads(:), rad_near(:)

      integer *8 i, j, jpatch, jquadstart, jstart

      integer *8 ipars(1)
      complex *16 zpars(1)
      real *8 timeinfo(10), t1, t2

      real *8 ttot, done
      real *8 rfac, rfac0
      integer *8 iptype_avg, norder_avg
      integer *8 ikerorder, iquadtype, npts_over

      real *8 alpha, beta

      integer *8 ns
      real *8, allocatable :: sigmaover(:)
      real *8, allocatable :: sources(:,:), srctmp(:,:), srctmp2(:,:)
      real *8, allocatable :: charges(:), dipvec(:,:)
      real *8, allocatable :: ctmp2(:), dtmp2(:,:)
      real *8, allocatable :: pot_aux(:), grad_aux(:,:)

      integer *8 ifcharge, ifdipole, ifpgh, ifpghtarg
      integer *8 ier, iper
      real *8 pottmp, gradtmp(3)
      integer *8 nss, l
      real *8 thresh, ra
      integer *8 nmax
      integer *8 nd, ntarg0
      integer *8 int8_3, int8_12

      parameter (nd=1, ntarg0=1)
      int8_3  = 3
      int8_12 = 12
      done = 1

      iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
      norder_avg = floor(sum(norders)/(npatches+0.0d0))

      call get_rfacs(norder_avg, iptype_avg, rfac, rfac0)

      allocate(cms(3,npatches), rads(npatches), rad_near(npatches))

      call get_centroid_rads(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, cms, rads)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
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
      call get_iquad_rsc(npatches, ixyzs, ntarg, nnz, row_ptr, col_ind, &
        iquad)

!     ikerorder: 0 for S' only (alpha only), 1 if D' is present
      ikerorder = 0
      if (abs(dpars(2)) .gt. 1.0d-16) ikerorder = 1

      allocate(novers(npatches), ixyzso(npatches+1))

!     Laplace: pass zpars = 0
      zpars(1) = 0.0d0
      call get_far_order(eps, npatches, norders, ixyzs, iptype, cms, &
        rads, npts, srccoefs, ndtarg, ntarg, targs, ikerorder, zpars(1), &
        nnz, row_ptr, col_ind, rfac, novers, ixyzso)

      npts_over = ixyzso(npatches+1)-1

      allocate(srcover(12,npts_over), wover(npts_over))

      call oversample_geom(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, srcvals, novers, ixyzso, npts_over, srcover)

      call get_qwts(npatches, novers, ixyzso, iptype, npts_over, &
        srcover, wover)

      nquad = iquad(nnz+1)-1
      allocate(wnear(nquad))

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i = 1,nquad
        wnear(i) = 0
      enddo
!$OMP END PARALLEL DO

      iquadtype = 1

      call getnearquad_lap_comb_cprime_eval(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, dpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)

!
!   compute layer potential via FMM
!   S'[sigma](x) = grad_x S[sigma](x) . n(x)
!   D'[sigma](x) = grad_x D[sigma](x) . n(x)
!
!   Far field: call lfmm3d_t_cd_g with
!     charges = alpha * sigma_over * wts  (for S')
!     dipoles = beta  * sigma_over * wts * n_src  (for D')
!   then dot grad_aux with target normal
!
      ns = npts_over
      alpha = dpars(1)
      beta  = dpars(2)

      allocate(sources(3,ns), srctmp(3,ntarg))
      allocate(charges(ns), dipvec(3,ns))
      allocate(sigmaover(ns))
      allocate(pot_aux(ntarg), grad_aux(3,ntarg))

      call oversample_fun_surf(nd, npatches, norders, ixyzs, iptype, &
        npts, sigma, novers, ixyzso, ns, sigmaover)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i = 1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)
        charges(i)   = alpha*sigmaover(i)*wover(i)
        dipvec(1,i)  = beta*sigmaover(i)*wover(i)*srcover(10,i)
        dipvec(2,i)  = beta*sigmaover(i)*wover(i)*srcover(11,i)
        dipvec(3,i)  = beta*sigmaover(i)*wover(i)*srcover(12,i)
      enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i = 1,ntarg
        srctmp(1,i)   = targs(1,i)
        srctmp(2,i)   = targs(2,i)
        srctmp(3,i)   = targs(3,i)
        pot(i)        = 0
        pot_aux(i)    = 0
        grad_aux(1,i) = 0
        grad_aux(2,i) = 0
        grad_aux(3,i) = 0
      enddo
!$OMP END PARALLEL DO

      ifcharge = 1
      ifdipole = 1
      if (abs(alpha) .lt. 1.0d-16) ifcharge = 0
      if (abs(beta)  .lt. 1.0d-16) ifdipole = 0

      ier = 0
      if (ifcharge .eq. 1 .and. ifdipole .eq. 1) then
        call lfmm3d_t_cd_g(eps, ns, sources, charges, dipvec, &
          ntarg, srctmp, pot_aux, grad_aux, ier)
      elseif (ifcharge .eq. 1 .and. ifdipole .eq. 0) then
        call lfmm3d_t_c_g(eps, ns, sources, charges, &
          ntarg, srctmp, pot_aux, grad_aux, ier)
      elseif (ifcharge .eq. 0 .and. ifdipole .eq. 1) then
        call lfmm3d_t_d_g(eps, ns, sources, dipvec, &
          ntarg, srctmp, pot_aux, grad_aux, ier)
      endif

!     dot gradient with target normal
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i = 1,ntarg
        pot(i) = grad_aux(1,i)*targs(10,i) + &
                 grad_aux(2,i)*targs(11,i) + &
                 grad_aux(3,i)*targs(12,i)
      enddo
!$OMP END PARALLEL DO

      call get_fmm_thresh(int8_12, ns, srcover, int8_12, ntarg, srctmp, &
        thresh)

!     Add precomputed near-field quadrature
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart) &
!$OMP PRIVATE(jstart,pottmp,npols,l)
      do i = 1,ntarg
        do j = row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch)
          do l = 1,npols
            pot(i) = pot(i) + wnear(jquadstart+l-1)*sigma(jstart+l-1)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

!     Remove near contribution from FMM (add-subtract)
      call get_near_corr_max(ntarg, row_ptr, nnz, col_ind, npatches, &
        ixyzso, nmax)
      allocate(srctmp2(3,nmax), ctmp2(nmax), dtmp2(3,nmax))

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2) &
!$OMP PRIVATE(ctmp2,dtmp2,nss,l,jstart,pottmp,gradtmp)
      do i = 1,ntarg
        nss = 0
        do j = row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l = ixyzso(jpatch),ixyzso(jpatch+1)-1
            nss = nss+1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)
            if (ifcharge .eq. 1) ctmp2(nss)   = charges(l)
            if (ifdipole .eq. 1) then
              dtmp2(1,nss) = dipvec(1,l)
              dtmp2(2,nss) = dipvec(2,l)
              dtmp2(3,nss) = dipvec(3,l)
            endif
          enddo
        enddo

        pottmp     = 0
        gradtmp(1) = 0
        gradtmp(2) = 0
        gradtmp(3) = 0

        if (ifcharge .eq. 1 .and. ifdipole .eq. 1) then
          call l3ddirectcdg(nd, srctmp2, ctmp2, dtmp2, &
            nss, srctmp(1,i), ntarg0, pottmp, gradtmp, thresh)
        elseif (ifcharge .eq. 1 .and. ifdipole .eq. 0) then
          call l3ddirectcg(nd, srctmp2, ctmp2, &
            nss, srctmp(1,i), ntarg0, pottmp, gradtmp, thresh)
        elseif (ifcharge .eq. 0 .and. ifdipole .eq. 1) then
          call l3ddirectdg(nd, srctmp2, dtmp2, &
            nss, srctmp(1,i), ntarg0, pottmp, gradtmp, thresh)
        endif

        pot(i) = pot(i) - (gradtmp(1)*targs(10,i) + &
                            gradtmp(2)*targs(11,i) + &
                            gradtmp(3)*targs(12,i))
      enddo
!$OMP END PARALLEL DO

      return
      end subroutine lap_comb_cprime_eval
!
!
!


      subroutine lap_comb_cprime_eval_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        eps, ndd, dpars, ndz, zpars, ndi, ipars, &
        nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
        nptso, ixyzso, srcover, whtsover, lwork, work, idensflag, &
        ndim_s, sigma, ipotflag, ndim_p, pot)
!
!  This subroutine evaluates the combined prime layer potential
!
!    K[\sigma] = (\alpha S_{0}' + \beta D_{0}')[\sigma]
!
!  where the near field quadrature is precomputed and stored
!  in row sparse compressed format, and the oversampled geometry
!  is provided by the caller.
!
!  For targets on the boundary, this routine returns only the
!  principal value part.
!
!  The FMM computes grad_x (alpha*S_0 + beta*D_0)[sigma] at targets,
!  dotted with the target normal n_x, using add-subtract to remove
!  the near contribution and replace it with the precomputed wnear.
!
!  Input arguments:
!    - npatches: integer *8
!    - norders: integer *8(npatches)
!    - ixyzs: integer *8(npatches+1)
!    - iptype: integer *8(npatches)
!    - npts: integer *8
!    - srccoefs: real *8(9,npts)
!    - srcvals: real *8(12,npts)
!    - ndtarg: integer *8
!        leading dimension of target array (must be >= 12)
!    - ntarg: integer *8
!    - targs: real *8(ndtarg,ntarg)
!        rows 10:12 must contain target normals
!    - eps: real *8
!    - ndd: integer *8
!        must be 2
!    - dpars: real *8(ndd)
!        dpars(1) = alpha, dpars(2) = beta
!    - ndz: integer *8 (unused, pass 0)
!    - zpars: complex *16(ndz) (unused)
!    - ndi: integer *8 (unused, pass 0)
!    - ipars: integer *8(ndi) (unused)
!    - nnz: integer *8
!    - row_ptr: integer *8(ntarg+1)
!    - col_ind: integer *8(nnz)
!    - iquad: integer *8(nnz+1)
!    - nquad: integer *8
!    - nker: integer *8 (must be 1)
!    - wnear: real *8(nker,nquad)
!        precomputed near quadrature for K = alpha*S' + beta*D'
!    - novers: integer *8(npatches)
!    - nptso: integer *8
!    - ixyzso: integer *8(npatches+1)
!    - srcover: real *8(12,nptso)
!    - whtsover: real *8(nptso)
!    - lwork: integer *8 (unused)
!    - work: real *8(lwork) (unused)
!    - idensflag: integer *8 (unused)
!    - ndim_s: integer *8 (must be 1)
!    - sigma: real *8(npts)
!    - ipotflag: integer *8 (unused)
!    - ndim_p: integer *8 (must be 1)
!
!  Output arguments:
!    - pot: real *8(ntarg)
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
      complex *16, intent(in) :: zpars(max(ndz,1))
      integer *8, intent(in) :: ipars(max(ndi,1))

      integer *8, intent(in) :: nnz, nquad
      integer *8, intent(in) :: row_ptr(ntarg+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)

      integer *8, intent(in) :: nker
      real *8, intent(in) :: wnear(nker,nquad)

      integer *8, intent(in) :: nptso
      integer *8, intent(in) :: ixyzso(npatches+1), novers(npatches)
      real *8, intent(in) :: srcover(12,nptso), whtsover(nptso)

      integer *8, intent(in) :: lwork
      real *8, intent(in) :: work(max(lwork,1))

      integer *8, intent(in) :: ndim_s, ndim_p
      integer *8, intent(in) :: idensflag, ipotflag

      real *8, intent(in) :: sigma(npts)
      real *8, intent(out) :: pot(ntarg)

      integer *8 ns
      real *8 alpha, beta
      integer *8 ifcharge, ifdipole

      real *8, allocatable :: sources(:,:), srctmp(:,:), srctmp2(:,:)
      real *8, allocatable :: charges(:), dipvec(:,:), sigmaover(:)
      real *8, allocatable :: ctmp2(:), dtmp2(:,:)
      real *8, allocatable :: pot_aux(:), grad_aux(:,:)

      integer *8 i, j, jpatch, jquadstart, jstart
      integer *8 nss, l, npols, nmax
      real *8 pottmp, gradtmp(3), thresh
      integer *8 ier

      integer *8 nd, ntarg0
      integer *8 int8_3, int8_12
      parameter (nd=1, ntarg0=1)
      int8_3  = 3
      int8_12 = 12

      ns = nptso
      alpha = dpars(1)
      beta  = dpars(2)

      ifcharge = 1
      ifdipole = 1
      if (abs(alpha) .lt. 1.0d-16) ifcharge = 0
      if (abs(beta)  .lt. 1.0d-16) ifdipole = 0

      allocate(sources(3,ns), srctmp(3,ntarg))
      allocate(charges(ns), dipvec(3,ns))
      allocate(sigmaover(ns))
      allocate(pot_aux(ntarg), grad_aux(3,ntarg))

!     oversample density
      call oversample_fun_surf(nd, npatches, norders, ixyzs, iptype, &
        npts, sigma, novers, ixyzso, ns, sigmaover)

!     set FMM source/target arrays
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i = 1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)
        charges(i)   = alpha*sigmaover(i)*whtsover(i)
        dipvec(1,i)  = beta*sigmaover(i)*whtsover(i)*srcover(10,i)
        dipvec(2,i)  = beta*sigmaover(i)*whtsover(i)*srcover(11,i)
        dipvec(3,i)  = beta*sigmaover(i)*whtsover(i)*srcover(12,i)
      enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i = 1,ntarg
        srctmp(1,i)   = targs(1,i)
        srctmp(2,i)   = targs(2,i)
        srctmp(3,i)   = targs(3,i)
        pot(i)        = 0
        pot_aux(i)    = 0
        grad_aux(1,i) = 0
        grad_aux(2,i) = 0
        grad_aux(3,i) = 0
      enddo
!$OMP END PARALLEL DO

!     call FMM: compute gradient of (alpha*S + beta*D)[sigma] at targets
      ier = 0
      if (ifcharge .eq. 1 .and. ifdipole .eq. 1) then
        call lfmm3d_t_cd_g(eps, ns, sources, charges, dipvec, &
          ntarg, srctmp, pot_aux, grad_aux, ier)
      elseif (ifcharge .eq. 1 .and. ifdipole .eq. 0) then
        call lfmm3d_t_c_g(eps, ns, sources, charges, &
          ntarg, srctmp, pot_aux, grad_aux, ier)
      elseif (ifcharge .eq. 0 .and. ifdipole .eq. 1) then
        call lfmm3d_t_d_g(eps, ns, sources, dipvec, &
          ntarg, srctmp, pot_aux, grad_aux, ier)
      endif

!     dot gradient with target normal to get S'/D' contribution
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i = 1,ntarg
        pot(i) = grad_aux(1,i)*targs(10,i) + &
                 grad_aux(2,i)*targs(11,i) + &
                 grad_aux(3,i)*targs(12,i)
      enddo
!$OMP END PARALLEL DO

!     threshold for ignoring near direct interactions
      call get_fmm_thresh(int8_12, ns, srcover, int8_12, ntarg, &
        srctmp, thresh)

!     add precomputed near quadrature
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart) &
!$OMP PRIVATE(jstart,npols,l)
      do i = 1,ntarg
        do j = row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch)
          do l = 1,npols
            pot(i) = pot(i) + wnear(1,jquadstart+l-1)*sigma(jstart+l-1)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

!     subtract near FMM contribution (add-subtract correction)
      call get_near_corr_max(ntarg, row_ptr, nnz, col_ind, npatches, &
        ixyzso, nmax)
      allocate(srctmp2(3,nmax), ctmp2(nmax), dtmp2(3,nmax))

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2) &
!$OMP PRIVATE(ctmp2,dtmp2,nss,l,jstart,pottmp,gradtmp)
      do i = 1,ntarg
        nss = 0
        do j = row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l = ixyzso(jpatch),ixyzso(jpatch+1)-1
            nss = nss+1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)
            if (ifcharge .eq. 1) ctmp2(nss)   = charges(l)
            if (ifdipole .eq. 1) then
              dtmp2(1,nss) = dipvec(1,l)
              dtmp2(2,nss) = dipvec(2,l)
              dtmp2(3,nss) = dipvec(3,l)
            endif
          enddo
        enddo

        pottmp     = 0
        gradtmp(1) = 0
        gradtmp(2) = 0
        gradtmp(3) = 0

        if (ifcharge .eq. 1 .and. ifdipole .eq. 1) then
          call l3ddirectcdg(nd, srctmp2, ctmp2, dtmp2, &
            nss, srctmp(1,i), ntarg0, pottmp, gradtmp, thresh)
        elseif (ifcharge .eq. 1 .and. ifdipole .eq. 0) then
          call l3ddirectcg(nd, srctmp2, ctmp2, &
            nss, srctmp(1,i), ntarg0, pottmp, gradtmp, thresh)
        elseif (ifcharge .eq. 0 .and. ifdipole .eq. 1) then
          call l3ddirectdg(nd, srctmp2, dtmp2, &
            nss, srctmp(1,i), ntarg0, pottmp, gradtmp, thresh)
        endif

        pot(i) = pot(i) - (gradtmp(1)*targs(10,i) + &
                            gradtmp(2)*targs(11,i) + &
                            gradtmp(3)*targs(12,i))
      enddo
!$OMP END PARALLEL DO

      return
      end subroutine lap_comb_cprime_eval_addsub
