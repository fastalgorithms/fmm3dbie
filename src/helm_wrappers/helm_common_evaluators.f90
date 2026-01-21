      subroutine getnearquad_helm_comb_wdd_eval(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)
!
!  This subroutine generates the near field quadrature
!  for the representation:
!
!  u = alpha*S_{k}[\sigma_{1}] + beta*D_{k}[\sigma_{2}]    -   (1)
!
!  This routine is just a wrapper to the combined field
!  layer potential evaluator with two different densities
!  in helm_common_evaluators.f90                                                
!  
!  For targets on surface, this routine returns the principal
!  value part of D. The user is responsible for adding the
!  contribution of the identity term.            
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
!  NOTES:
!    - wnear must be of size (2,nquad) as 2 different layer
!      potentials are returned
!      * the first kernel is S_{k}
!      * the second kernel is D_{k}
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
!    - zpars: complex *16 (1)
!        kernel parameters 
!        zpars(1) = k, Helmholtz wavenumber
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
!        radius parameter for switching to predetermined quadarature
!        rule        
!    - nquad: integer *8
!        number of near field entries corresponding to each source target
!        pair
!
!  Output arguments
!    - wnear: complex *16(2,nquad)
!        The desired near field quadrature
!        wnear(1,:) - stores the quadrature corrections for S_{k}
!        wnear(2,:) - stores the quadrature correction for D_{k}        
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
      complex *16, intent(in) :: zpars(1)
      integer *8, intent(in) :: iquadtype, nnz
      integer *8, intent(in) :: row_ptr(ntarg+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: nquad
      complex *16, intent(out) :: wnear(2,nquad)

!  Temporary variables      
      complex *16 zk, ima, zpars_tmp(3)
      integer *8 ipv, ndd, ndi, ndz, i
      real *8 dpars
      integer *8 ipars
      complex *16, allocatable :: wneartmp(:)
  
      procedure (), pointer :: fker
      external h3d_slp, h3d_dlp

      data ima/(0.0d0,1.0d0)/
    
      ndd = 0
      ndz = 1
      ndi = 0


      allocate(wneartmp(nquad))
    
    
      if (iquadtype.eq.1) then
        ipv = 0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wneartmp(i) = 0 
        enddo
!$OMP END PARALLEL DO
        fker => h3d_slp
        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, &
          zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, &
          rfac0, nquad, wneartmp)

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

        ipv = 1
        fker => h3d_dlp
        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, &
          zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, &
          rfac0, nquad, wneartmp)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wnear(2,i) = wneartmp(i)
        enddo
!$OMP END PARALLEL DO
      endif


      return
      end
!




      subroutine helm_comb_wdd_eval_addsub(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        eps, ndd, dpars, ndz, zpars, ndi, ipars, &
        nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
        nptso, ixyzso, srcover, whtsover, lwork, work, idensflag, &
        ndim_s, sigma, ipotflag, ndim_p, pot)

!
!
!  This subroutine computes the potential u
!  for the representation:
!
!  u = alpha*S_{k}[\sigma_{1}] + beta*D_{k}[\sigma_{2}] 
!                         
! 
!  where the near field is precomputed and stored
!  in the row sparse compressed format.
!        
!  For targets on surface, this routine returns the principal
!  value part of D. The user is responsible for adding the
!  contribution of the identity term.            
!
!  The fmm is used to accelerate the far-field and 
!  near-field interactions are handled via precomputed quadrature
!
!  Using add and subtract - no need to call tree and set fmm parameters
!  can directly call existing fmm library
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
!        leading dimension of target information array
!    - ntarg: integer *8
!        number of targets
!    - targs: real *8(ndtarg,ntarg)
!        target information 
!    - eps: real *8
!        precision requested
!    - ndd: integer *8
!        number of real parameters defining the kernel/
!        integral representation (unused in this routine)
!    - dpars: real *8(ndd)
!        real parameters defining the kernel/
!        integral representation (unused in this routine)
!    - ndz: integer *8
!        number of complex parameters defining the kernel/
!        integral representation (must be three)
!    - zpars: complex *16(ndz)
!        complex parameters defining the kernel/
!        integral represnetation.
!        * zpars(1) = k 
!        * zpars(2) = alpha
!        * zpars(3) = beta    
!    - ndi: integer *8
!        number of integer parameters defining the kernel/
!        integral representation (unused in this routine)
!    - ipars: integer *8(ndi)
!        integer parameters defining the kernel/
!        integral represnetation (unused in this routine)
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
!    - nquad: integer *8
!        number of near field entries corresponding to each source target
!        pair
!    - nker: integer *8
!        number of kernels in quadrature correction, must be 2
!    - wnear: complex *16(nker, nquad)
!        precomputed quadrature corrections            
!    - novers: integer *8(npatches)
!        order of discretization for oversampled sources and
!        density
!    - ixyzso: integer *8(npatches+1)
!        ixyzso(i) denotes the starting location in srcover,
!        corresponding to patch i
!    - nptso: integer *8
!        total number of oversampled points
!    - srcover: real *8 (12,nptso)
!        oversampled set of source information
!    - whtsover: real *8 (nptso)
!        smooth quadrature weights at oversampled nodes
!    - lwork: integer *8
!        size of work array (unused in this routine)
!    - work: real *8(lwork)
!        work array (unused in this routine)
!    - idensflag: integer *8
!        Flag for types of denisties (unused in this case)
!    - ndim_s: integer *8
!        number of densities per point on the surface,
!        must be 2 for this routine
!    - sigma: complex *16(ndim_s, npts)
!        The densities \sigma_{1}, and \sigma_{2}
!    - ipotflag: integer *8
!        Flag for determining output type 
!        ipotflag = 1, only potential is returned
!        ipotflag = 2, potential + gradients are
!        returned (currently unsupported)       
!    - ndim_p: integer *8
!        number of potentials per point
!        ndim_p = 1, if ipotflag = 1
!        ndim_p = 4, if ipotflag = 2                 
!    
!  Output arguments:
!    - pot: complex *16 (ndim_p, ntarg)
!        u above
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
      complex *16, intent(in) :: wnear(nker,nquad)

      integer *8, intent(in) :: nptso
      integer *8, intent(in) :: ixyzso(npatches+1), novers(npatches)
      real *8, intent(in) :: srcover(12,nptso), whtsover(nptso)
    
      integer *8, intent(in) :: lwork
      real *8, intent(in) :: work(lwork)

      integer *8, intent(in) :: ndim_s, ndim_p
      integer *8, intent(in) :: idensflag, ipotflag
  
      complex *16, intent(in) :: sigma(ndim_s, npts)
  
      complex *16, intent(out) :: pot(ndim_p, ntarg)

      real *8, allocatable :: sources(:,:), targtmp(:,:)
      complex *16, allocatable :: charges(:), dipvec(:,:), sigmaover(:,:)
      complex *16, allocatable :: potloc(:)

      integer *8 npols
      integer *8 ns
      integer *8 ifpgh, ifpghtarg

      integer *8 i, j, jpatch, jquadstart, jstart
      complex *16 pottmp

!$    real *8 omp_get_wtime      

      real *8, allocatable :: srctmp2(:,:)
      complex *16, allocatable :: ctmp2(:), dtmp2(:,:)
      real *8 thresh

      integer *8 nss, l, ier
      complex *16 ima, ztmp, ztmp2

      integer *8 nd, ntarg0, nmax

      real *8 done, rr, rtmp
      integer *8 nddtmp, nditmp, ndztmp
      integer *8 int8_12
      complex *16 zpars2(2)
      data ima/(0.0d0,1.0d0)/

      parameter (nd=1,ntarg0=1)

      ns = nptso
      int8_12 = 12

      ifpgh = 0
      ifpghtarg = 1
  
      allocate(sources(3,ns), targtmp(3,ntarg))
      allocate(charges(ns), dipvec(3,ns))
      allocate(sigmaover(2,ns))
!
!  estimate max number of sources in the near field of any target
!
!    
      call get_near_corr_max(ntarg, row_ptr, nnz, col_ind, npatches, &
        ixyzso, nmax)
      allocate(srctmp2(3,nmax), ctmp2(nmax), dtmp2(3,nmax))
! 
!       oversample density
!
      call oversample_fun_surf(2*ndim_s, npatches, norders, ixyzs, &
        iptype, npts, sigma, novers, ixyzso, ns, sigmaover)

!
!  extract source and target info
!
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(rtmp, ztmp, ztmp2)
      do i=1,ns
        rtmp = whtsover(i)
        ztmp = rtmp*zpars(2)
        ztmp2 = rtmp*zpars(3)
    

        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)
        charges(i) = sigmaover(1,i)*ztmp
        dipvec(1,i) = sigmaover(2,i)*srcover(10,i)*ztmp2
        dipvec(2,i) = sigmaover(2,i)*srcover(11,i)*ztmp2
        dipvec(3,i) = sigmaover(2,i)*srcover(12,i)*ztmp2
      enddo
!$OMP END PARALLEL DO      

      allocate(potloc(ntarg))

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        targtmp(1,i) = targs(1,i)
        targtmp(2,i) = targs(2,i)
        targtmp(3,i) = targs(3,i)
        pot(1,i) = 0
        potloc(i) = 0
      enddo

!$OMP END PARALLEL DO

! 
!  Compute S_{k}'
      ier = 0
      call hfmm3d_t_cd_p(eps, zpars(1), ns, sources, charges, dipvec, &
        ntarg, targtmp, potloc, ier)    
  
!  Compute threshold for ignoring local computation
      call get_fmm_thresh(int8_12, ns, srcover, ndtarg, ntarg, targs, thresh)
!
!  Add near field precomputed contribution
! 
  

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, jpatch, jquadstart) &
!$OMP PRIVATE(jstart, npols, l)
      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1) - ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            potloc(i) = potloc(i) + &
              wnear(1,jquadstart+l-1)*sigma(1,jstart+l-1)*zpars(2)
            potloc(i) = potloc(i) + &
              wnear(2,jquadstart+l-1)*sigma(2,jstart+l-1)*zpars(3)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

!
!     Remove near contribution of the FMM
! 

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, jpatch, srctmp2) &
!$OMP PRIVATE(ctmp2, dtmp2, l, jstart, nss, pottmp)
      do i=1,ntarg
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
            nss = nss+1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)
            ctmp2(nss) = charges(l)
            dtmp2(1,nss) = dipvec(1,l)
            dtmp2(2,nss) = dipvec(2,l)
            dtmp2(3,nss) = dipvec(3,l)
          enddo
        enddo

        pottmp = 0

        call h3ddirectcdp(nd, zpars(1), srctmp2, ctmp2, dtmp2, nss, &
          targtmp(1,i), ntarg0, pottmp, thresh)
        potloc(i) = potloc(i) - pottmp
        pot(1,i) = potloc(i)
      enddo
!$OMP END PARALLEL DO      

      return
      end
!
!
!
!
