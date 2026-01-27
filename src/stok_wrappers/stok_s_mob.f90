!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  Single layer representation for the Stokes mobility problem 
!
!  PDE:
!    \Delta u       = \nabla p
!    \Delta \cdot u = 0
!
!  Boundary conditions:
!    u = u_{i} + \omega_{i} \times (x - x_{c,i}) on component \Gamma_{i}
!      where u_{i}, \omega_{i} are unknown
!    \int_{\Gamma_{i}} f dS = F_{i}, f is the surface traction
!    \int_{\Gamma_{i}} (x - x_{c,i}) \times f ds = T_{i}
!      where F_{i}, T_{i} are the forces and torques on body i and
!      are given
!
!
!  Representation:
!    u = S_{stok}[\sigma] + S_{stok}[\sigma_{0}] 
!  
!  where \sigma_{0} = F_{i}/|\Gamma_{i}| + 
!     \tau_{i}^{-1} T_{i} \times (x - x_{c,i})
!
!  with \tau_{i} is the moment of inertia tensor. In addition, 
!
!  \sigma satisfies
!    \int_{\Gamma_{i}} \sigma = 0, 
!    \int_{\Gamma_{i}} (x-x_{c,i}) \times \sigma = 0
!
!  Integral equations obtained by imposing:
!    f^{-} = 0 \,,
!
!  where f^{-} is the interior surface traction
!
!  and is given by:
!    1/2 \sigma + S_{stok}'[\sigma] = -S_{stok}'[\sigma_{0}]
!
!  along with the constraints above. 
!
!  Instead, we use the generalized 1's matrix trick to solve for
!  sigma as
!    1/2 \sigma + S_{stok}'[\sigma] + L[\sigma] = -S_{stok}'[\sigma_{0}]
!
!  where L[\sigma] is given by
!     = 1/|\Gamma_{i}|\int_{\Gamma_{i}} \sigma dS + 
!       \tau_{i}^{-1}(\int_{{\Gamma_{i}}} (y-x_{c,i}) \times \sigma dS) \times (x-x_{c,i})
!  on \Gamma_{i}
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  User callable routines:
!    - stok_s_mob_solver: Given data F_{i}, T_{i}, 
!       this routine returns the solution \sigma, and u_{i}, \omega_{i}
!
!    - stok_s_vel_eval: Given \sigma this routine
!        evaluates the solution at a collection of targets
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  Advanced interfaces:
!    - getnearquad_stok_s_mob: compute the near quadrature correction
!        for constructing the on-surface integral equation for the
!        traction data corresponding to the single layer
!        representation with user-provided near-information prescribed
!        in row-sparse compressed format
!
!    - lpcomp_stok_s_mob_addsub: apply the principal value part
!        of the integeral equation on surface. On input, user provides
!        precomputed near quadrature in row-sparse compressed format,
!        and oversampling information for handling the far part of the
!        computation
!
!    - stok_s_mob_solver_guru: Guru solver routine, where user is
!        responsible for providing precomputed near quadrature
!        information in row-sparse compressed format and oversampling
!        surface information for the far-part
!
!

      subroutine getnearquad_stok_s_mob(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, eps, &
        iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear)
!
!  This subroutine generates the near field quadrature
!  for the representation:
!
!  u = S_{stok}[\sigma] - (1) 
!
!  and returns quantities related to evaluating the surface 
!  traction of u on surface at the surface discretization nodes.
!
!  On imposing the boundary condition, we get the following operator
!
!  u = \sigma/2 + S_{stok}'[\sigma] 
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
!        radius parameter for switching to predetermined quadarature
!        rule        
!    - nquad: integer *8
!        number of near field entries corresponding to each source 
!        target pair
!
!  Output arguments
!    - wnear: real *8(6,nquad)
!        The desired near field quadrature for
!        \alpha S_{stok} + \beta D_{stok}. Since the Stokes 
!        tensor is symmetric, quadratures are only
!        generated for the lower half of the tensor
!        * wnear(1,:) - stores the quadratures for (1,1) entry of the 
!                       tensor 
!        * wnear(2,:) - stores the quadratures for (1,2) entry of the 
!                       tensor 
!        * wnear(3,:) - stores the quadratures for (1,3) entry of the 
!                       tensor 
!        * wnear(4,:) - stores the quadratures for (2,2) entry of the 
!                       tensor 
!        * wnear(5,:) - stores the quadratures for (2,3) entry of the 
!                       tensor 
!        * wnear(6,:) - stores the quadratures for (3,3) entry of the 
!                       tensor 
!
      implicit none 
      integer *8, intent(in) :: npatches, norders(npatches), npts, nquad
      integer *8, intent(in) :: ixyzs(npatches+1), iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts), eps
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: iquadtype
      integer *8, intent(in) :: nnz
      integer *8, intent(in) :: row_ptr(npts+1), col_ind(nnz), iquad(nnz+1)
      real *8, intent(out) :: wnear(6,nquad)

      integer *8 :: ndtarg, ntarg
      integer *8, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)

      real *8 alpha, beta
      real *8, allocatable :: wneartmp(:)
      integer *8 ndd, ndi, ndz
      complex *16 zpars(1)
      integer *8 ipars(2)

      integer *8 i, ipv, j, ii, ijloc(2,6)

      procedure (), pointer :: fker
      external st3d_strac 
      

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

      ndd = 2
      ndi = 2
      ndz = 0
      fker => st3d_strac
      ipv = 1

      ijloc(1,1) = 1
      ijloc(2,1) = 1
      ijloc(1,2) = 1
      ijloc(2,2) = 2
      ijloc(1,3) = 1
      ijloc(2,3) = 3
      ijloc(1,4) = 2
      ijloc(2,4) = 2
      ijloc(1,5) = 2
      ijloc(2,5) = 3
      ijloc(1,6) = 3
      ijloc(2,6) = 3

      allocate(wneartmp(nquad))


      if(iquadtype.eq.1) then

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
      subroutine stok_s_vel_eval(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, sigma, pot)
!
!f2py intent(in) npatches,norders,ixyzs,iptype,npts,srccoefs,srcvals
!f2py intent(in) ndtarg,ntarg,targs,ipatch_id,uvs_targ,eps
!f2py intent(in) sigma
!f2py intent(out) pot
!
!
!------------------------------
!  This subroutine evaluates the layer potential for the 
!  Stokes single layer representation 
!
!  Representation:
!    u = \mathcal{S}_{stok}[\sigma]
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
!    - ipatch_id: integer *8(ntarg)
!        id of patch of target i, id = -1, if target is off-surface
!    - uvs_targ: double precision (2,ntarg)
!        local uv coordinates on patch if on surface, otherwise
!        set to 0 by default
!    - eps: double precision
!        precision requested
!     - sigma: double precision(3,npts)
!         density for layer potential
!
!  Output arguments
!    - pot: double precision(3,ntarg)
!        layer potential evaluated at the target points
!
!-----------------------------------
!
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

      real *8 dpars(2)

      dpars(1) = 1.0d0
      dpars(2) = 0.0d0

      call stok_comb_vel_eval(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, dpars, sigma, pot)

      return
      end
!
!
!  
!
!
      subroutine lpcomp_stok_s_mob_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, &
        novers, nptso, ixyzso, srcover, whtsover, lwork, work, ndim, &
        sigma, pot)
!
!  This subroutine evaluates the surface traction +
!  the generalized ones matrix corresponding to
!  the following integral representation:
!  
!  u = S_{stok}[\sigma]
!
!  pot = S_{stok}'[\sigma] + L[\sigma]
!
!  where L[\sigma] is given by
!     = 1/|\Gamma_{i}|\int_{\Gamma_{i}} \sigma dS + 
!       tau_{i}^{-1}(\int_{{\Gamma_{i}}} (y-x_{c,i}) \times \sigma dS) \times (x-x_{c,i})
!  on \Gamma_{i}
!
!  The near field is precomputed and stored
!  in the row sparse compressed format.
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
!    - eps: real *8
!        precision requested
!    - ndd: integer *8
!        number of real parameters defining the kernel/
!        integral representation (must be 3)
!    - dpars: real *8(ndd)
!        real parameters defining the kernel/
!        integral representation 
!        * dpars(1) = alpha (single layer strength) 
!        * dpars(2) = beta (double layer strength)
!        * dpars(3) = strength of 1s matrix of 
!            n \int \sigma \cdot n
!    - ndz: integer *8
!        number of complex parameters defining the kernel/
!        integral representation (unused in this routine)
!    - zpars: complex *16(ndz)
!        complex parameters defining the kernel/
!        integral represnetation (unused in this routine)
!    - ndi: integer *8
!        must be ncomp + 1, where ncomp is the number of
!        components in the geometry
!    - ipars: integer *8(ndi)
!        ipars(1) = ncomp
!        ipars(i+1) = starting patch index for patches on component
!          i (note that the patches must be ordered by component)
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
!    - nquad: integer *8
!        number of near field entries corresponding to each source target
!        pair
!    - nker: integer *8
!        number of kernels in quadrature correction, must be 6
!    - wnear: real *8(nker, nquad)
!        Precomputed near field quadrature for
!        S_{stok}'. Since the Stokes 
!        tensor is symmetric, quadratures are only
!        generated for the lower half of the tensor
!        * wnear(1,:) - stores the quadratures for (1,1) entry of the 
!                       tensor 
!        * wnear(2,:) - stores the quadratures for (1,2) entry of the 
!                       tensor 
!        * wnear(3,:) - stores the quadratures for (1,3) entry of the 
!                       tensor 
!        * wnear(4,:) - stores the quadratures for (2,2) entry of the 
!                       tensor 
!        * wnear(5,:) - stores the quadratures for (2,3) entry of the 
!                       tensor 
!        * wnear(6,:) - stores the quadratures for (3,3) entry of the 
!                       tensor 
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
!        size of work array (must be npts)
!    - work: real *8(lwork)
!        work array
!        * work(1:npts) stores the coordinates the quadrature
!          weights for integrating smooth functions
!    - ndim: integer *8
!        number of densities per point on the surface,
!        must be 3 for this routine
!    - sigma: real *8(3,npts)
!        The density sigma above                                        
!
!  Output arguments:
!    - pot: real *8(3,npts)
!        u corresponding to representation
!

  
      implicit none
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches),ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)

      real *8, intent(in) :: eps
      
      integer *8, intent(in) :: ndd, ndz, ndi
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      integer *8, intent(in) :: ipars(ndi)

      integer *8, intent(in) :: nnz, nquad
      integer *8, intent(in) :: row_ptr(npts+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      
      integer *8, intent(in) :: nker
      real *8, intent(in) :: wnear(nker,nquad)

      integer *8, intent(in) :: nptso
      integer *8, intent(in) :: ixyzso(npatches+1), novers(npatches)
      real *8, intent(in) :: srcover(12,nptso), whtsover(nptso)
      
      integer *8, intent(in) :: lwork
      real *8, intent(in) :: work(lwork)

      integer *8, intent(in) :: ndim
    
      real *8, intent(in) :: sigma(ndim,npts)
    
      real *8, intent(out) :: pot(ndim,npts)

      integer *8 ndtarg, idensflag, ipotflag, ndim_p, i
      real *8 rint

      integer *8 norder,npols,nover,npolso
      real *8, allocatable :: potsort(:)

      real *8, allocatable :: sources(:,:), targvals(:,:)
      real *8, allocatable :: sigmaover(:,:)
      real *8, allocatable :: stoklet(:,:)
      real *8, allocatable :: pottarg(:,:), gradtarg(:,:,:)
      real *8, allocatable :: pretarg(:)
      integer *8 ns, nt
      real *8 alpha, beta
      integer *8 ifstoklet, ifstrslet
      integer *8 ifppreg, ifppregtarg
      real *8 tmp(10), val
      real *8 over4pi

      real *8 w11, w12, w13, w21, w22, w23, w31, w32, w33
      real *8 sig1, sig2, sig3

      integer *8 istress

      integer *8 i, j, jpatch, jquadstart, jstart

      real *8 pottmp
      real *8, allocatable :: sttmp2(:,:)
      real *8 strstmp2(3), strsvec2(3)

      integer *8 nmax
      real *8 timeinfo(10), t1, t2, omp_get_wtime

      real *8, allocatable :: srctmp2(:,:)
      real *8 thresh, ra
      real *8 rr, rmin
      integer *8 nss, ii, l, npover

      integer *8 nd, ntarg0, ntarg
      integer *8 ier, iper

      integer *8 ndtmp, ndsigma

      real *8 ttot, done, pi

      real *8 potv(3), pv, gradv(3,3)
      real *8 strslet(3), strvec(3)

      real *8 rint(3), rint_torque(3), xyzc(3), xdiff(3), rtmp(3)
      integer *8 ncomp, icomp

      parameter (nd=1,ntarg0=1)
      data over4pi/0.07957747154594767d0/

      ns = nptso
      done = 1
      pi = atan(done)*4

!
!    estimate max number of sources in neear field of 
!    any target
!
      
      
      nmax = 0
      ntarg = npts
      call get_near_corr_max(ntarg, row_ptr, nnz, col_ind, npatches, &
        ixyzso, nmax)
      allocate(srctmp2(3,nmax), sttmp2(3,nmax))
           
           
      ifppreg = 0
      ifppregtarg = 3
      allocate(sources(3,ns), targvals(3,ntarg))
      allocate(stoklet(3,ns), strslet(3,ns), strsvec(3,ns))
      allocate(sigmaover(3,ns))
      allocate(pottarg(3,ntarg), pretarg(ntarg))
      allocate(gradtarg(3,3,ntarg))

! 
!       oversample density
!

      ndsigma = 3
      call oversample_fun_surf(ndsigma, npatches, norders, ixyzs, &
        iptype, npts, sigma, novers, ixyzso, ns, sigmaover)



!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)      
      do i=1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)

        stoklet(1,i) = sigmaover(1,i)*whtsover(i)*over4pi
        stoklet(2,i) = sigmaover(2,i)*whtsover(i)*over4pi
        stoklet(3,i) = sigmaover(3,i)*whtsover(i)*over4pi
      enddo
!$OMP END PARALLEL DO      
      
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i = 1,ntarg
        targvals(1,i) = srcvals(1,i)
        targvals(2,i) = srcvals(2,i)
        targvals(3,i) = srcvals(3,i)
      enddo
!$OMP END PARALLEL DO      
      

      ifstoklet = 1
      ifstrslet = 0

      iper = 0
      ier = 0
!
!       call the fmm
!
      call cpu_time(t1)
!$      t1 = omp_get_wtime()      
      call stfmm3d(nd, eps, ns, sources, ifstoklet, stoklet, &
        ifstrslet, strslet, strsvec, ifppreg, tmp, tmp, tmp, ntarg, &
        targvals, ifppregtarg, pottarg, pretarg, gradtarg, ier)
      call cpu_time(t2)
!$      t2 = omp_get_wtime()

      timeinfo(1) = t2-t1
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i = 1,ntarg
        pot(1,i) = srcvals(10,i)*(-pretarg(i) + 2*gradtarg(1,1,i)) + &
             srcvals(11,i)*(gradtarg(1,2,i) + gradtarg(2,1,i)) + &
             srcvals(12,i)*(gradtarg(1,3,i) + gradtarg(3,1,i))
        pot(2,i) = srcvals(10,i)*(gradtarg(1,2,i) + gradtarg(2,1,i)) + &
             srcvals(11,i)*(-pretarg(i) + 2*gradtarg(2,2,i)) + &
             srcvals(12,i)*(gradtarg(2,3,i) + gradtarg(3,2,i))
        pot(3,i) = srcvals(10,i)*(gradtarg(1,3,i) + gradtarg(3,1,i)) + &
             srcvals(11,i)*(gradtarg(2,3,i) + gradtarg(3,2,i)) + &
             srcvals(12,i)*(-pretarg(i) + 2*gradtarg(3,3,i))
      enddo
!$OMP END PARALLEL DO 

!
!        compute threshold for ignoring local computation
!
      ndtmp = 3 
      call get_fmm_thresh(ndtmp, ns, sources, ndtmp, ntarg, targvals, &
        thresh)
      
!
!
!       add in precomputed quadrature
!

      call cpu_time(t1)
!$      t1 = omp_get_wtime()

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j, jpatch, jquadstart) &
!$OMP PRIVATE(jstart, npols, l, sig1, sig2, sig3, w11, w12, w13, w21) &
!$OMP PRIVATE(w22, w23, w31, w32, w33)
      do i = 1,ntarg
        do j = row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1) - ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            sig1 = sigma(1,jstart+l-1)
            sig2 = sigma(2,jstart+l-1)
            sig3 = sigma(3,jstart+l-1)

            w11 = wnear(1,jquadstart+l-1)
            w12 = wnear(2,jquadstart+l-1)
            w13 = wnear(3,jquadstart+l-1)

            w21 = w12 
            w22 = wnear(4,jquadstart+l-1)
            w23 = wnear(5,jquadstart+l-1)

            w31 = w13 
            w32 = w23 
            w33 = wnear(6,jquadstart+l-1)


            pot(1,i) = pot(1,i) + w11*sig1 + w12*sig2 + w13*sig3 
            pot(2,i) = pot(2,i) + w21*sig1 + w22*sig2 + w23*sig3 
            pot(3,i) = pot(3,i) + w31*sig1 + w32*sig2 + w33*sig3 
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO


!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, jpatch, srctmp2) &
!$OMP PRIVATE(sttmp2, strstmp2, strsvec2, nss, l, jstart, ii, npover) &
!$OMP PRIVATE(potv, pv, gradv, istress)
      do i = 1,ntarg
        nss = 0
        do j = row_ptr(i), row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l = ixyzso(jpatch),ixyzso(jpatch+1)-1
            nss = nss + 1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)

            sttmp2(1,nss) = stoklet(1,l)
            sttmp2(2,nss) = stoklet(2,l)
            sttmp2(3,nss) = stoklet(3,l)
          enddo
        enddo

        potv(1) = 0
        potv(2) = 0
        potv(3) = 0
        pv = 0
        gradv(1:3,1:3) = 0
        
        istress = 0
        call st3ddirectstokstrsg(nd, srctmp2, sttmp2, istress, &
          strstmp2, strsvec2, nss, targvals(1,i), ntarg0, potv, pv, &
          gradv, thresh)
        potv(1) = srcvals(10,i)*(-pv + 2*gradv(1,1)) + &
             srcvals(11,i)*(gradv(1,2) + gradv(2,1)) + &
             srcvals(12,i)*(gradv(1,3) + gradv(3,1))
        potv(2) = srcvals(10,i)*(gradv(1,2) + gradv(2,1)) + &
             srcvals(11,i)*(-pv + 2*gradv(2,2)) + &
             srcvals(12,i)*(gradv(2,3) + gradv(3,2))
        potv(3) = srcvals(10,i)*(gradv(1,3) + gradv(3,1)) + &
             srcvals(11,i)*(gradv(2,3) + gradv(3,2)) + &
             srcvals(12,i)*(-pv + 2*gradv(3,3))

        pot(1,i) = pot(1,i) - potv(1)
        pot(2,i) = pot(2,i) - potv(2)
        pot(3,i) = pot(3,i) - potv(3)

      enddo
      
      call cpu_time(t2)
!$      t2 = omp_get_wtime()     


!
!  Now add in contribution of L[\sigma]
!
      ncomp = ipars(1)
      do icomp = 1,ncomp
        rint(1:3) = 0
        rint_torque(1:3) = 0
        ipstart = ipars(icomp+1)
        ipend = npatches
        if (icomp.ne.ncomp) ipend = ipars(icomp+2)-1

        istart = ixyzs(ipstart)
        iend = ixyzs(ipend)
        nploc = ipend - ipstart + 1
        nptsloc = iend - istart + 1

        do i=1,nploc+1
          ixyzsloc(i) = ixyzs(i+istart-1) - ixyzs(istart)+1
        enddo
        rmoi(1:3,1:3) = 0
        rmoi_inv(1:3,1:3) = 0
        area = 0
        centroid(1:3) = 0
        call get_surf_moments(nploc, norders(ipstart), ixyzsloc, &
          iptype(ipstart), nptsloc, srcvals(1,istart), work(istart), &
          area, centroid, rmoi)
        info = 0
        call dinverse(3, rmoi, info, rmoi_inv) 
        
        do i = istart,iend
          rint(1:3) = rint(1:3) + sigma(1:3,i)*work(i)
          rtmp(1:3) = 0
          xdiff(1:3) = srcvals(1:3,i) - centroid(1:3)
          call cross_prod3d(xdiff, sigma(1:3,i), rtmp)
          rint_torque(1:3) = rint_torque(1:3) + rtmp(1:3)*work(i)
        enddo

        do i = istart,iend
          pot(1:3,i) = pot(1:3,i) + rint(1:3)/area
          xdiff(1:3) = srcvals(1:3,i) - centroid(1:3)
          rtmp(1:3) = 0
          call cross_prod3d(rint_torque, xdiff, rtmp)
          pot(1,i) = pot(1,i) + rmoi_inv(1,1)*rtmp(1) + &
            rmoi_inv(1,2)*rtmp(2) + rmoi_inv(1,3)*rtmp(3)
          pot(2,i) = pot(2,i) + rmoi_inv(2,1)*rtmp(1) + &
            rmoi_inv(2,2)*rtmp(2) + rmoi_inv(2,3)*rtmp(3)
          pot(3,i) = pot(3,i) + rmoi_inv(3,1)*rtmp(1) + &
            rmoi_inv(3,2)*rtmp(2) + rmoi_inv(3,3)*rtmp(3)
        enddo
      enddo

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
      subroutine stok_s_mob_solver(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ncomp, icomps, eps, numit, &
        forces, torques, eps_gmres, niter, errs, rres, soln, &
        trans_vels, rot_vels)

!
!
!  This subroutine solves the Stokes mobility problem
!  where the potential is represented using a 
!  single layer representation.
!
!  This subroutine is the simple interface as opposed to the
!  _solver_guru routine which is called after initialization
!  in this routine.
!
!
!  Representation:
!    u = S_{stok}[\sigma] + S_{stok}[\sigma_{0}] 
!  
!  where \sigma_{0} = F_{i}/|\Gamma_{i}| + 
!     \tau_{i}^{-1} T_{i} \times (x - x_{c,i})
!
!  with \tau_{i} is the moment of inertia tensor. In addition, 
!
!  \sigma satisfies
!    \int_{\Gamma_{i}} \sigma = 0, 
!    \int_{\Gamma_{i}} (x-x_{c,i}) \times \sigma = 0
!
!  Integral equations obtained by imposing:
!    f^{-} = 0 \,,
!
!  where f^{-} is the interior surface traction
!
!  and is given by:
!    1/2 \sigma + S_{stok}'[\sigma] = -S_{stok}'[\sigma_{0}]
!
!  along with the constraints above. 
!
!  Instead, we use the generalized 1's matrix trick to solve for
!  sigma as
!    1/2 \sigma + S_{stok}'[\sigma] + L[\sigma] = -S_{stok}'[\sigma_{0}]
!
!  where L[\sigma] is given by
!     = 1/|\Gamma_{i}|\int_{\Gamma_{i}} \sigma dS + 
!
!  The linear system is solved iteratively using GMRES
!  until a relative residual of eps_gmres is reached
!
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
!    - ncomp: integer *8
!        number of components
!    - icomps: integer *8(ncomp)
!        icomps(i) = starting patch index for patches on component
!          i (note that the patches must be ordered by component)
!    - eps: real *8
!        precision requested for computing quadrature and fmm
!        tolerance
!    - numit: integer *8
!        max number of gmres iterations
!    - ifinout: integer *8
!        ifinout = 0, interior problem
!        ifinout = 1, exterior problem
!    - forces: real *8(3,ncomp)
!        prescribed force on the components
!    - torques: real *8(3,ncomp)
!        prescribed torque on the components
!    - eps_gmres: real *8
!        gmres tolerance requested
!      
!
!  Output arguments:
!    - niter: integer *8
!        number of gmres iterations required for relative residual
!        to converge to eps_gmres
!    - errs:  real *8 (numit+1)
!        relative residual as a function of iteration
!        number (only errs(1:niter) is populated))
!    - rres: real *8
!        relative residual for computed solution
!    - soln: real *8(3,npts)
!        density which solves the mobility problem \sigma
!    - trans_vels: real *8(3,ncomp)
!        translational velocity of the components
!    - rot_vels: real *8(3,ncomp)
!        rotational velocity of the components
!
      implicit none
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: ifinout
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      integer *8, intent(in) :: ncomp, icomps(ncomp)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps, eps_gmres
      real *8, intent(in) :: forces(3,ncomp), torques(3,ncomp) 
      integer *8, intent(in) :: numit
      real *8, intent(out) :: soln(3,npts)
      real *8, intent(out) :: trans_vels(3,ncomp), rot_vels(3,ncmp)
      real *8, intent(out) :: errs(numit+1)
      real *8, intent(out) :: rres
      integer *8, intent(out) :: niter

      integer *8 norder, npols
      real *8, allocatable :: targs(:,:)
      integer *8, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      integer *8 ndtarg, ntarg



      integer *8 nover, npolso, nptso
      integer *8 nnz,nquad
      integer *8, allocatable :: row_ptr(:), col_ind(:), iquad(:)
      real *8, allocatable :: wnear(:,:), wnear_s(:,:)

      real *8, allocatable :: srcover(:,:), wover(:)
      integer *8, allocatable :: ixyzso(:), novers(:)

      real *8, allocatable :: cms(:,:), rads(:), rad_near(:) 
      real *8, allocatable :: rhs(:,:), sigma0(:,:), uvel(:,:)

      integer *8 i, j, jpatch, jquadstart, jstart

      integer *8, allocatable :: ipars(:)
      complex *16 zpars
      real *8 dpars(2)
      real *8 timeinfo(10), t1, t2, omp_get_wtime

      real *8 ttot, done, pi
      real *8 rfac, rfac0
      integer *8 iptype_avg, norder_avg
      integer *8 ikerorder, iquadtype, npts_over

      real *8 xdiff(3), rtmp(3), rtuse(3), ftuse(3), rtmp2(3)
      real *8 rmoi(3,3), rmoi_inv(3,3), centroid(3)
      integer *8, allocatable :: ixyzsloc(:)
      real *8, allocatable :: wts(:)
 
!
!
!       gmres variables
!
      real *8 did, dtmp
      complex *16 ztmp
      integer *8 nker


      done = 1
      pi = atan(done)*4
!
!
!        setup targets as on surface discretization points
! 
      ndtarg = 3
      ntarg = npts
      allocate(targs(ndtarg,npts), uvs_targ(2,ntarg), ipatch_id(ntarg))

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        targs(1,i) = srcvals(1,i)
        targs(2,i) = srcvals(2,i)
        targs(3,i) = srcvals(3,i)
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo
!$OMP END PARALLEL DO   


!
!    initialize patch_id and uv_targ for on surface targets
!
      call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts, &
        ipatch_id, uvs_targ)

!
!
!        this might need fixing
!
      iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
      norder_avg = floor(sum(norders)/(npatches+0.0d0))

      call get_rfacs(norder_avg, iptype_avg, rfac, rfac0)


      allocate(cms(3,npatches), rads(npatches), rad_near(npatches))

      call get_centroid_rads(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, cms, rads)

!$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,npatches
        rad_near(i) = rads(i)*rfac
      enddo
!$OMP END PARALLEL DO      

!
!    find near quadrature correction interactions
!
      print *, "entering find near mem"
      call findnearmem(cms, npatches, rad_near, ndtarg, targs, npts, nnz)

      allocate(row_ptr(npts+1), col_ind(nnz))
      
      call findnear(cms, npatches, rad_near, ndtarg, targs, npts, &
        row_ptr, col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches, ixyzs, npts, nnz, row_ptr, col_ind, &
        iquad)

      ikerorder = -1
      if(abs(dpars(2)).gt.1.0d-16) ikerorder = 0

!
!    estimate oversampling for far-field, and oversample geometry
!

      allocate(novers(npatches), ixyzso(npatches+1))

      print *, "beginning far order estimation"

      ztmp = 0

      call get_far_order(eps, npatches, norders, ixyzs, iptype, cms, &
        rads, npts, srccoefs, ndtarg, npts, targs, ikerorder, ztmp, &
        nnz, row_ptr, col_ind, rfac, novers, ixyzso)

      npts_over = ixyzso(npatches+1) - 1

      allocate(srcover(12,npts_over), wover(npts_over))

      call oversample_geom(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, srcvals, novers, ixyzso, npts_over, srcover)

      call get_qwts(npatches, novers, ixyzso, iptype, npts_over, &
        srcover, wover)

!
!   compute near quadrature correction
!
      nquad = iquad(nnz+1) - 1
      nker = 6
      allocate(wnear(nker,nquad), wnear_s(nker,nquad))
      
!$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,nquad
        wnear(1,i) = 0
        wnear(2,i) = 0
        wnear(3,i) = 0
        wnear(4,i) = 0
        wnear(5,i) = 0
        wnear(6,i) = 0

        wnear_s(1,i) = 0
        wnear_s(2,i) = 0
        wnear_s(3,i) = 0
        wnear_s(4,i) = 0
        wnear_s(5,i) = 0
        wnear_s(6,i) = 0
      enddo
!$OMP END PARALLEL DO    


      iquadtype = 1

      print *, "starting to generate near quadrature"
      call cpu_time(t1)
!$      t1 = omp_get_wtime()     

!
!  quadrature for S'
!
      call getnearquad_stok_s_mob(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, eps, &
        iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear)
!
!  quadrature for S
!
      dpars(1) = 1.0d0
      dpars(2) = 0
      call getnearquad_stok_comb_vel(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, eps, dpars, iquadtype, &
        nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear_s)
      call cpu_time(t2)
!$      t2 = omp_get_wtime()     


!
!  Compute S[\sigma_{0}]
!
      allocate(wts(npts))
      call get_qwts(npatches, novers, ixyzs, iptype, npts, &
        srcvals, wts)
      allocate(sigma0(3,npts), rhs(3,npts), uvel(3,npts))
      allocate(ixyzsloc(npatches+1))
      ndd = 2
      ndz = 0
      ndi = 0
      nker = 6
      lwork = 0
      ndim_s = 3
      idensflag = 0
      ipotflag = 0
      ndim_p = 0

      allocate(ipars(ncomp+1))
      do icomp = 1,ncomp
        ipstart = ipars(icomp+1)
        ipend = npatches
        if (icomp.ne.ncomp) ipend = ipars(icomp+2)-1

        istart = ixyzs(ipstart)
        iend = ixyzs(ipend)
        nploc = ipend - ipstart + 1
        nptsloc = iend - istart + 1

        do i=1,nploc+1
          ixyzsloc(i) = ixyzs(i+istart-1) - ixyzs(istart)+1
        enddo
        rmoi(1:3,1:3) = 0
        rmoi_inv(1:3,1:3) = 0
        area = 0
        centroid(1:3) = 0
        call get_surf_moments(nploc, norders(ipstart), ixyzsloc, &
          iptype(ipstart), nptsloc, srcvals(1,istart), wts(istart), &
          area, centroid, rmoi)
        info = 0
        call dinverse(3, rmoi, info, rmoi_inv) 
        rtuse(1:3) = rtuse(1:3)
        rtuse(1) = rtuse(1) = rmoi_inv(1,1)*torques(1,icomp) + &
           rmoi_inv(1,2)*torques(2,icomp) + &
           rmoi_inv(1,3)*torques(3,icomp)
        rtuse(2) = rtuse(2) = rmoi_inv(2,1)*torques(1,icomp) + &
           rmoi_inv(2,2)*torques(2,icomp) + &
           rmoi_inv(2,3)*torques(3,icomp)
        rtuse(3) = rtuse(3) = rmoi_inv(3,1)*torques(1,icomp) + &
           rmoi_inv(3,2)*torques(2,icomp) + &
           rmoi_inv(3,3)*torques(3,icomp)
        ftuse(1:3) = forces(1:3,icomp)/area
        do i = istart,iend
          xdiff(1:3) = srcvals(1:3,i) - centroid(1:3)
          call cross_prod3d(rtuse, xdiff, rtmp)
          sigma0(1:3,i) = ftuse(1:3) + rtmp(1:3)
        enddo
      enddo

      call stok_comb_vel_eval_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        eps, ndd, dpars, ndz, zpars, ndi, ipars, nnz, row_ptr, &
        col_ind, iquad, nquad, nker, wnear_s, novers, npts_over, &
        ixyzso, srcover, wover, lwork, work, idensflag, ndim_s, &
        sigma0, ipotflag, ndim_p, rhs)
!$OMP PARALLEL DO PRIVATE(i)
      do i=1,npts
        rhs(1:3,i) = -rhs(1:3,i)
      enddo
!$OMP END PARALLEL DO 

      
      print *, "done generating near quadrature, now starting gmres"
      ndi = ncomp+1
      call stok_s_mob_solver_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndi, ipars, numit, &
        rhs, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
        npts_over, ixyzso, srcover, wover, eps_gmres, niter, &
        errs, rres, soln)

      call stok_comb_vel_eval_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        eps, ndd, dpars, ndz, zpars, ndi, ipars, nnz, row_ptr, &
        col_ind, iquad, nquad, nker, wnear_s, novers, npts_over, &
        ixyzso, srcover, wover, lwork, work, idensflag, ndim_s, &
        soln, ipotflag, ndim_p, uvel)
!
!  Compute translational and rotational velocities on each component
!
!  Note the sign flip in rhs as it was -S[\sigma_{0}]
      do i=1,npts
        uvel(1:3,i) = uvel(1:3,i) - rhs(1:3,i)
      enddo

      do icomp = 1,ncomp
        trans_vels(1:3,icomp) = 0
        rot_vels(1:3,icomp) = 0

        ipstart = ipars(icomp+1)
        ipend = npatches
        if (icomp.ne.ncomp) ipend = ipars(icomp+2)-1

        istart = ixyzs(ipstart)
        iend = ixyzs(ipend)
        nploc = ipend - ipstart + 1
        nptsloc = iend - istart + 1

        do i=1,nploc+1
          ixyzsloc(i) = ixyzs(i+istart-1) - ixyzs(istart)+1
        enddo
        rmoi(1:3,1:3) = 0
        rmoi_inv(1:3,1:3) = 0
        area = 0
        centroid(1:3) = 0
        call get_surf_moments(nploc, norders(ipstart), ixyzsloc, &
          iptype(ipstart), nptsloc, srcvals(1,istart), wts(istart), &
          area, centroid, rmoi)
        info = 0
        call dinverse(3, rmoi, info, rmoi_inv) 
        do i = istart,iend
          trans_vels(1:3,icomp) = trans_vels(1:3,icomp) + uvel(1:3,i)*wts(i)
        enddo
        trans_vels(1:3,icomp) = trans_vels(1:3,icomp)/area
        
        rtmp(1:3) = 0
        do i=1,istart,iend
          xdiff(1:3) = srcvals(1:3,i) - centroid(1:3)
          rtuse(1:3) = uvel(1:3,i) - trans_vels(1:3,icomp)
          rtmp2(1:3) = 0
          call cross_prod3d(xdiff, rtuse, rtmp2)
          rtmp(1:3) = rtmp(1:3) + rtmp2(1:3)*wts(i)
        enddo
        rot_vels(1:3,icomp) = rmoi_inv(1:3,1)*rtmp(1) + 
          rmoi_inv(1:3,2)*rtmp(2) + rmoi_inv(1:3,3)*rtmp(3)
      enddo

      
      
      return
      end
!
!
!
!
!

      subroutine stok_s_mob_solver_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndi, ipars, numit, &
        rhs, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
        nptso, ixyzso, srcover, whtsover, eps_gmres, niter, &
        errs, rres, soln)
!
!
!  This subroutine solves the Stokes mobility problem
!  where the potential is represented using a 
!  single layer representation.
!
!  This subroutine is the simple interface as opposed to the
!  _solver_guru routine which is called after initialization
!  in this routine.
!
!
!  Representation:
!    u = S_{stok}[\sigma] + S_{stok}[\sigma_{0}] 
!  
!  where \sigma_{0} = F_{i}/|\Gamma_{i}| + 
!     \tau_{i}^{-1} T_{i} \times (x - x_{c,i})
!
!  with \tau_{i} is the moment of inertia tensor. In addition, 
!
!  \sigma satisfies
!    \int_{\Gamma_{i}} \sigma = 0, 
!    \int_{\Gamma_{i}} (x-x_{c,i}) \times \sigma = 0
!
!  Integral equations obtained by imposing:
!    f^{-} = 0 \,,
!
!  where f^{-} is the interior surface traction
!
!  and is given by:
!    1/2 \sigma + S_{stok}'[\sigma] = -S_{stok}'[\sigma_{0}]
!
!  along with the constraints above. 
!
!  Instead, we use the generalized 1's matrix trick to solve for
!  sigma as
!    1/2 \sigma + S_{stok}'[\sigma] + L[\sigma] = -S_{stok}'[\sigma_{0}]
!
!  where L[\sigma] is given by
!     = 1/|\Gamma_{i}|\int_{\Gamma_{i}} \sigma dS + 
!
!  The linear system is solved iteratively using GMRES
!  until a relative residual of eps_gmres is reached
!
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
!    - eps: real *8
!        precision requested for computing quadrature and fmm
!        tolerance
!    - ndi: integer *8
!        number of ipars, must be atleast ncomp+1, where ncomp
!        is the number of components
!    - ipars: integer *8(ndi)
!        ipars(1) = ncomp
!        ipars(i+1) = starting patch index for patches on component
!          i (note that the patches must be ordered by component)
!    - numit: integer *8
!        max number of gmres iterations
!    - ifinout: integer *8
!        ifinout = 0, interior problem
!        ifinout = 1, exterior problem
!    - rhs: real *8(3,npts)
!        velocity data
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
!    - nquad: integer *8
!        number of near field entries corresponding to each source target
!        pair
!    - nker: integer *8
!        number of kernels in quadrature correction, must be 6
!    - wnear: real *8(nker, nquad)
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
!    - eps_gmres: real *8
!        gmres tolerance requested
!
!  output
!    - niter: integer *8
!        number of gmres iterations required for relative residual
!        to converge to eps_gmres
!    - errs:  real *8 (numit+1)
!        relative residual as a function of iteration
!        number (only errs(1:niter) is populated))
!    - rres: real *8
!        relative residual for computed solution
!    - soln: real *8(3,npts)
!        density which solves the neumann problem \sigma
!				 

      implicit none
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps, eps_gmres
      integer *8, intent(in) :: ndi
      integer *8, intent(in) :: ipars(ndi)
      real *8, intent(in) :: rhs(3,npts)
      integer *8, intent(in) :: numit

      real *8, intent(out) :: errs(numit+1)
      real *8, intent(out) :: rres
      integer *8, intent(out) :: niter
      
      integer *8, intent(in) :: nnz, nquad
      integer *8, intent(in) :: row_ptr(npts+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)

      integer *8, intent(in) :: nker
      real *8, intent(in) :: wnear(nker,nquad)

      integer *8, intent(in) :: nptso
      integer *8, intent(in) :: novers(npatches), ixyzso(npatches+1)
      real *8, intent(in) :: srcover(12,nptso), whtsover(nptso)
      
  
      real *8, intent(out) :: soln(3,npts)

      real *8 did
      real *8 dpars
      integer *8 ndd_use

      procedure (), pointer :: fker
      external lpcomp_stok_comb_vel_addsub

      integer *8 ndd, ndi, ndz, lwork, ndim
      integer *8 nkertmp
      complex *16 zpars

      integer *8 ndtarg
      complex *16 ima
      data ima/(0.0d0,1.0d0)/
      real *8, allocatable :: wts(:)

      integer *8 i

      did = 0.5d0 
      fker => lpcomp_stok_comb_vel_addsub

      ndd = 0
      ndz = 0

      lwork = npts
      ndim = 3
      allocate(wts(npts))

      call get_qwts(npatches, norders, ixyzs, iptype, npts, &
        srcvals, wts)

      call dgmres_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, wts, &
        eps, ndd, dpars, ndz, zpars, ndi, ipars, &
        nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
        nptso, ixyzso, srcover, whtsover, npts, wts, &
        ndim, fker, did, rhs, numit, eps_gmres, niter, errs, &
        rres, soln)

      return
      end

!
