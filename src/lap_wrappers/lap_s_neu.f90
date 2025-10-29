!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!   Single layer representation for Laplace Neumann problem
!
!  PDE:
!    \Delta u = 0
!
!  Boundary conditions:
!    du/dn = f
!
!  Representation:
!    u = S_{0}[\sigma]
!
!  Integral equations obtained by imposing:
!    du/dn = f
!
!  and is given by:
!    d \sigma + S_{0}[\sigma] = f
!
!  where:
!    d =  1/2 for interior problem
!      = -1/2 for exterior problem
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  User callable routines:
!    - lap_s_neu_solver: Given data f this routine returns the 
!        solution \sigma
!
!    - lap_s_neu_eval: Given \sigma this routine evaluates the solution 
!        at a collection of targets
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  Advanced interfaces:
!    - getnearquad_lap_s_neu: compute the near quadrature correction
!        for constructing the on-surface integral equation for the
!        Neumann data corresponding to the single layer 
!        representation with user-provided near-information prescribed
!        in row-sparse compressed format
!
!    - getnearquad_lap_s_neu_eval: compute the near quadrature
!        correction for target points which can be either on-surface or
!        off-surface with user-provided near-information prescribed in
!        row-sparse compressed format
!
!    - lpcomp_lap_s_neu_addsub: apply the principal value part
!        of the integeral equation on surface. On input, user provides
!        precomputed near quadrature in row-sparse compressed format,
!        and oversampling information for handling the far part of the
!        computation
!
!    - lap_s_neu_eval_addsub: compute the solution u at a
!        collection of targets (on-surface or off-surface), given \sigma.
!        On input, user provides precomputed near quadrature in
!        row-sparse compressed format and oversampling surface
!        information for the far-part
!
!    - lap_s_neu_solver_guru: Guru solver routine, where user is
!        responsible for providing precomputed near quadrature
!        information in row-sparse compressed format and oversampling
!        surface information for the far-part
!
!
!

      subroutine getnearquad_lap_s_neu(npatches, norders, &
       ixyzs, iptype, npts, srccoefs, srcvals, eps, iquadtype, nnz, & 
       row_ptr, col_ind, iquad, rfac0, nquad, wnear)
!
!  This subroutine generates the near field quadrature
!  for the representation:
!
!  u = S_{0}[\rho]
!
!  and returns quantities related to evaluating du/dn on surface
!  at the surface discretization nodes
!
!  If values at other nodes is desired then the solution
!  should be reinterpolated to those nodes
!
!
!  On imposing the boundary condition, we get the following operator
!
!  du/dn = \pm I/2 + S_{0}' 

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
!    - rfac0: integer
!        radius parameter for near field
!    - nquad: integer
!        number of near field entries corresponding to each source target
!        pair. 
!
!  Output arguments
!    - wnear: real *8(nquad)
!        The desired near field quadrature
!               
!

      implicit none 
      integer, intent(in) :: npatches, norders(npatches)
      integer, intent(in) :: ixyzs(npatches+1), iptype(npatches)
      integer, intent(in) :: npts
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps
      integer, intent(in) :: iquadtype, nnz
      integer, intent(in) :: row_ptr(npts+1), col_ind(nnz), iquad(nnz+1)
      real *8, intent(in) :: rfac0
      integer, intent(in) :: nquad
      
      real *8, intent(out) :: wnear(nquad)

      integer ndtarg, ntarg, ndd, ndz, ndi
      real *8 dpars(1)
      complex *16 zpars(1)
      integer ipars(1)

      real *8, allocatable :: uvs_targ(:,:)
      integer, allocatable :: ipatch_id(:)

      integer i
      integer ipv

      procedure (), pointer :: fker
      external l3d_sprime


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

      ipv = 1
      fker => l3d_sprime 
      call dgetnearquad_ggq_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, srcvals, &
        ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear)

      return
      end subroutine getnearquad_lap_s_neu
!
!
!
!
!
      subroutine lpcomp_lap_s_neu_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, &
        novers, nptso, ixyzso, srcover, whtsover, lwork, work, ndim, &
        sigma, pot)
!
!  This subroutine evaluates the Neumann data corresponding to
!  the following integral representation:
!  
!  u = S_{0}[\sigma] 
!
!  dudn = S_{0}'[\sigma] + c \int_{\Gamma} \sigma
!
!  where the one's matrix is added for the interior Neumann problem
!
!  where the near field is precomputed and stored
!  in the row sparse compressed format.
!        
!  The fmm is used to accelerate the far-field and 
!  near-field interactions are handled via precomputed quadrature
!
!  Using add and subtract - no need to call tree and set fmm parameters
!  can directly call existing fmm library
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
!        iptype = 12, quadrangular patch discretized with Chebyshev nodes
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
!    - ndd: integer
!        number of real parameters defining the kernel/
!        integral representation (must be 1)
!    - dpars: real *8(ndd)
!        real parameters defining the kernel/
!        integral representation 
!        * dpars(1) = c, the multiple of integral of \sigma
!    - ndz: integer
!        number of complex parameters defining the kernel/
!        integral representation (unused in this routine)
!    - zpars: complex *16(ndz)
!        complex parameters defining the kernel/
!        integral represnetation (unused in this routine)
!    - ndi: integer
!        number of integer parameters defining the kernel/
!        integral representation (unused in this routine)
!    - ipars: integer(ndi)
!        integer parameters defining the kernel/
!        integral represnetation (unused in this routine)
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
!    - nquad: integer
!        number of near field entries corresponding to each source target
!        pair
!    - nker: integer
!        number of kernels in quadrature correction, must be 1
!    - wnear: real *8(nker, nquad)
!        precomputed quadrature corrections            
!    - novers: integer(npatches)
!        order of discretization for oversampled sources and
!        density
!    - ixyzso: integer(npatches+1)
!        ixyzso(i) denotes the starting location in srcover,
!        corresponding to patch i
!    - nptso: integer
!        total number of oversampled points
!    - srcover: real *8 (12,nptso)
!        oversampled set of source information
!    - whtsover: real *8 (nptso)
!        smooth quadrature weights at oversampled nodes
!    - lwork: integer
!        size of work array (unused in this routine)
!    - work: real *8(lwork)
!        work array (unused in this routine)
!    - ndim: integer
!        number of densities per point on the surface,
!        must be 1 for this routine
!    - sigma: real *8(npts)
!        The density sigma above                                        
!
!  Output arguments:
!    - pot: real *8(npts)
!        dudn corresponding to representation
!
      implicit none
      integer, intent(in) :: npatches, npts
      integer, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)

      real *8, intent(in) :: eps
      
      integer, intent(in) :: ndd, ndz, ndi
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      integer, intent(in) :: ipars(ndi)

      integer, intent(in) :: nnz, nquad
      integer, intent(in) :: row_ptr(npts+1), col_ind(nnz)
      integer, intent(in) :: iquad(nnz+1)
      
      integer, intent(in) :: nker
      real *8, intent(in) :: wnear(nker,nquad)

      integer, intent(in) :: nptso
      integer, intent(in) :: ixyzso(npatches+1), novers(npatches)
      real *8, intent(in) :: srcover(12,nptso), whtsover(nptso)
      
      integer, intent(in) :: lwork
      real *8, intent(in) :: work(lwork)

      integer, intent(in) :: ndim
    
      real *8, intent(in) :: sigma(npts)
    
      real *8, intent(out) :: pot(npts)

      real *8, allocatable :: sources(:,:), srctmp(:,:)
      real *8, allocatable :: charges(:), sigmaover(:)
      integer ns, nt, ntarg
      integer ifpgh, ifpghtarg

      integer i, j, jpatch, jquadstart, jstart, npols
      real *8 pottmp, gradtmp(3)
      real *8, allocatable :: pot_aux(:), grad_aux(:,:)

      real *8 timeinfo(10), t1, t2, omp_get_wtime

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      real *8, allocatable :: ctmp2(:)
      real *8 thresh, ra
      real *8 rr, rmin
      real *8 over4pi, rint
      integer nss, ii, l, npover, ier

      integer nd, ntarg0, nmax

      real *8 ttot, done, pi
      data over4pi/0.07957747154594767d0/

      parameter (nd=1, ntarg0=1)

      ns = nptso
      ntarg = npts

      done = 1
      pi = atan(done)*4

      ifpgh = 0
      ifpghtarg = 2
      allocate(sources(3,ns), srctmp(3,npts))
      allocate(charges(ns))
      allocate(sigmaover(ns))
      allocate(pot_aux(npts), grad_aux(3,npts))
!
!  estimate max number of sources in the near field of any target
!
!    
      call get_near_corr_max(npts, row_ptr, nnz, col_ind, npatches, & 
        ixyzso, nmax)
      allocate(srctmp2(3,nmax), ctmp2(nmax))
! 
!       oversample density
!
      call oversample_fun_surf(nd, npatches, norders, ixyzs, iptype, & 
        npts, sigma, novers, ixyzso, ns, sigmaover)

!
!  extract source and target info
!
!$OMP PARALLEL DO DEFAULT(SHARED)
      do i = 1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)
        charges(i) = sigmaover(i)*whtsover(i)*over4pi
      enddo
!$OMP END PARALLEL DO      


!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        srctmp(1,i) = srcvals(1,i)
        srctmp(2,i) = srcvals(2,i)
        srctmp(3,i) = srcvals(3,i)
        pot_aux(i) = 0
        grad_aux(1,i) = 0
        grad_aux(2,i) = 0
        grad_aux(3,i) = 0
      enddo
!$OMP END PARALLEL DO

! 
!  Compute S_{0}'
!
      call lfmm3d_t_c_g(eps, ns, sources, charges, npts, &
        srctmp, pot_aux, grad_aux, ier)

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        pot(i) = grad_aux(1,i)*srcvals(10,i) + &
                 grad_aux(2,i)*srcvals(11,i) + &
                 grad_aux(3,i)*srcvals(12,i)
      enddo
!$OMP END PARALLEL DO      
    
        
!        compute threshold for ignoring local computation
      call get_fmm_thresh(12, ns, srcover, 12, npts, srcvals, thresh)
!
!       Add near field precomputed contribution
!
      call cpu_time(t1)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, jpatch, jquadstart) &
!$OMP PRIVATE(jstart, npols, l)
      do i = 1,npts
        do j = row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1) - ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l = 1,npols
            pot(i) = pot(i) + wnear(1,jquadstart+l-1)*sigma(jstart+l-1)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

!
!     Remove near contribution of the FMM
!

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, jpatch, srctmp2) &
!$OMP PRIVATE(ctmp2, l, jstart, nss, pottmp, gradtmp)
      do i = 1,npts
        nss = 0
        do j = row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l = ixyzso(jpatch),ixyzso(jpatch+1)-1
            nss = nss + 1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)
            ctmp2(nss) = charges(l)
          enddo
        enddo

        pottmp = 0
        gradtmp(1) = 0
        gradtmp(2) = 0
        gradtmp(3) = 0

        call l3ddirectcg(nd, srctmp2, ctmp2, nss, srctmp(1,i), &
          ntarg0, pottmp, gradtmp, thresh)
        pot(i) = pot(i) - gradtmp(1)*srcvals(10,i) - &
          gradtmp(2)*srcvals(11,i) - gradtmp(3)*srcvals(12,i)
      enddo
!$OMP END PARALLEL DO      

      rint = 0
!$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rint)
      do i = 1,nptso
        rint = rint + sigmaover(i)*whtsover(i)
      enddo
!$OMP END PARALLEL DO      
      rint = rint*dpars(1)
!$OMP PARALLEL DO DEFAULT(SHARED)      
      do i = 1,npts
        pot(i) = pot(i) + rint
      enddo
!$OMP END PARALLEL DO

      return
      end subroutine lpcomp_lap_s_neu_addsub
!
!
!
!
!        
      subroutine lap_s_neu_solver(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, numit, ifinout, &
        rhs, eps_gmres, niter, errs, rres, soln)

!
!
!  This subroutine solves the Laplace Neumann problem
!  on the exterior/interior of an object where the potential
!  is represented as a single layer potential 
!
!  This subroutine is the simple interface as opposed to the
!  _solver_guru routine which is called after initialization
!  in this routine.
!
!
!  Representation:
!    u = S_{0}[\sigma]
!
!  Boundary condition:
!    dudn = f
!
!  The linear system is solved iteratively using GMRES
!  until a relative residual of eps_gmres is reached
!
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
!        iptype = 12, quadrangular patch discretized with Chebyshev nodes
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
!        precision requested for computing quadrature and fmm
!        tolerance
!    - numit: integer
!        max number of gmres iterations
!    - ifinout: integer
!        ifinout = 0, interior problem
!        ifinout = 1, exterior problem
!    - rhs: real *8(npts)
!        Neumann data
!    - eps_gmres: real *8
!        gmres tolerance requested
!
!  Output arguments:
!    - niter: integer
!        number of gmres iterations required for relative residual
!        to converge to eps_gmres
!    - errs:  real *8 (numit+1)
!        relative residual as a function of iteration
!        number (only errs(1:niter) is populated))
!    - rres: real *8
!        relative residual for computed solution
!    - soln: real *8(npts)
!        density which solves the Neumann problem \sigma
!				 
      implicit none
      integer, intent(in) :: npatches, npts
      integer, intent(in) :: ifinout
      integer, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps, eps_gmres
      real *8, intent(in) :: rhs(npts)
      integer, intent(in) :: numit
      real *8, intent(out) :: soln(npts)
      real *8, intent(out) :: errs(numit+1)
      real *8, intent(out) :: rres
      integer, intent(out) :: niter

      integer norder, npols
      real *8, allocatable :: targs(:,:)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      integer ndtarg, ntarg



      integer nover, npolso, nptso
      integer nnz,nquad
      integer, allocatable :: row_ptr(:), col_ind(:), iquad(:)
      real *8, allocatable :: wnear(:)

      real *8, allocatable :: srcover(:,:), wover(:)
      integer, allocatable :: ixyzso(:), novers(:)

      real *8, allocatable :: cms(:,:), rads(:), rad_near(:) 

      integer i, j, jpatch, jquadstart, jstart

      integer ipars
      complex *16 zpars
      real *8 timeinfo(10), t1, t2, omp_get_wtime


      real *8 ttot, done, pi
      real *8 rfac, rfac0
      integer iptype_avg, norder_avg
      integer ikerorder, iquadtype, npts_over

!
!
!       gmres variables
!
      real *8 did, dtmp
      complex *16 ztmp
      integer nker


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

      ikerorder = 0

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
      allocate(wnear(nquad))
      
!$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,nquad
        wnear(i) = 0
      enddo
!$OMP END PARALLEL DO    


      iquadtype = 1

      print *, "starting to generate near quadrature"
      call cpu_time(t1)
!$      t1 = omp_get_wtime()      

      call getnearquad_lap_s_neu(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, eps, iquadtype, &
        nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear)
      call cpu_time(t2)
!$      t2 = omp_get_wtime()     

      print *, "done generating near quadrature, now starting gmres"
      
      nker = 1
      call lap_s_neu_solver_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, numit, ifinout, &
        rhs, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
        npts_over, ixyzso, srcover, wover, eps_gmres, niter, &
        errs, rres, soln)
      
      return
      end
!
!
!
!
!

      subroutine lap_s_neu_solver_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, numit, ifinout, &
        rhs, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
        nptso, ixyzso, srcover, whtsover, eps_gmres, niter, &
        errs, rres, soln)
!
!
!  This subroutine solves the Laplace Neumann problem
!  on the interiro/exterior of an object where the potential
!  is represented using a single layer potential 
!
!
!  Representation:
!    u = S_{0}[\sigma] 
!
!  Boundary condition:
!    u = f
!
!  The linear system is solved iteratively using GMRES
!  until a relative residual of eps_gmres is reached
!
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
!        iptype = 12, quadrangular patch discretized with Chebyshev nodes
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
!        precision requested for computing quadrature and fmm
!        tolerance
!    - numit: integer
!        max number of gmres iterations
!    - ifinout: integer
!        ifinout = 0, interior problem
!        ifinout = 1, exterior problem
!    - rhs: real *8(npts)
!        Dirichlet data
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
!    - nquad: integer
!        number of near field entries corresponding to each source target
!        pair
!    - nker: integer
!        number of kernels in quadrature correction
!    - wnear: real *8(nker, nquad)
!        precomputed quadrature corrections 
!    - novers: integer(npatches)
!        order of discretization for oversampled sources and
!        density
!    - ixyzso: integer(npatches+1)
!        ixyzso(i) denotes the starting location in srcover,
!        corresponding to patch i
!    - nptso: integer
!        total number of oversampled points
!    - srcover: real *8 (12,nptso)
!        oversampled set of source information
!    - whtsover: real *8 (nptso)
!        smooth quadrature weights at oversampled nodes
!    - eps_gmres: real *8
!        gmres tolerance requested
!
!  output
!    - niter: integer
!        number of gmres iterations required for relative residual
!        to converge to eps_gmres
!    - errs:  real *8 (numit+1)
!        relative residual as a function of iteration
!        number (only errs(1:niter) is populated))
!    - rres: real *8
!        relative residual for computed solution
!    - soln: real *8(npts)
!        density which solves the neumann problem \sigma
!				 

      implicit none
      integer, intent(in) :: npatches, npts
      integer, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps, eps_gmres
      integer, intent(in) :: ifinout
      real *8, intent(in) :: rhs(npts)
      integer, intent(in) :: numit

      real *8, intent(out) :: errs(numit+1)
      real *8, intent(out) :: rres
      integer, intent(out) :: niter
      
      integer, intent(in) :: nnz, nquad
      integer, intent(in) :: row_ptr(npts+1), col_ind(nnz)
      integer, intent(in) :: iquad(nnz+1)

      integer, intent(in) :: nker
      real *8, intent(in) :: wnear(nker,nquad)

      integer, intent(in) :: nptso
      integer, intent(in) :: novers(npatches), ixyzso(npatches+1)
      real *8, intent(in) :: srcover(12,nptso), whtsover(nptso)
      
      real *8, intent(out) :: soln(npts)

      real *8 did, rint

      procedure (), pointer :: fker
      external lpcomp_lap_s_neu_addsub

      integer ndd, ndi, ndz, lwork, ndim, i
      real *8 work
      integer ipars, nkertmp
      complex *16 zpars
      real *8 dpars

      integer ndtarg
      real *8, allocatable :: wts(:)


      did = 0.5d0*(-1)**(ifinout)
      fker => lpcomp_lap_s_neu_addsub

      ndd = 1
      ndi = 0
      ndz = 0

      lwork = 0
      ndim = 1
      allocate(wts(npts))

      call get_qwts(npatches, norders, ixyzs, iptype, npts, &
        srcvals, wts)
      
      dpars = 0
      if (ifinout.eq.0) then
        rint = 0.0d0
!$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rint)      
        do i = 1,nptso
          rint = rint + whtsover(i) 
        enddo
!$OMP END PARALLEL DO
        dpars = 1.0d0/rint
      endif  

      call dgmres_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, wts, &
        eps, ndd, dpars, ndz, zpars, ndi, ipars, &
        nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
        nptso, ixyzso, srcover, whtsover, lwork, work, &
        ndim, fker, did, rhs, numit, eps_gmres, niter, errs, &
        rres, soln)

      return
      end

!
!
!
!
      subroutine getnearquad_lap_s_neu_eval(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, &
        targs, ipatch_id, uvs_targ, eps, iquadtype, nnz, &
        row_ptr, col_ind, iquad, rfac0, nquad, wnear)
!
!  This subroutine generates the near field quadrature
!  for the representation:
!
!  u = S_{0}[\sigma]
!
!  and returns quantities related to evaluating u both on and off
!  surface.
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
!    - ndtarg: integer
!        leading dimension of target information array
!    - ntarg: integer
!        number of targets
!    - targs: real *8(ndtarg,ntarg)
!        target information 
!    - ipatch_id: integer(ntarg)
!        ipatch_id(i) indicates the patch on which target i
!        is, if it is on surface. ipatch_id(i) should be 0 
!        otherwise
!    - uvs_targ: real *8(2,ntarg)
!        if ipatch_id(i) > 0, then uvs_targ(1:2,i) are the
!        local uv coordinates of the target on the patch,
!    - eps: real *8
!        precision requested
!    - iquadtype: integer
!        quadrature type
!          * iquadtype = 1, use ggq for self + adaptive integration
!            for rest
!    - nnz: integer
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer(ntarg+1)
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
!
!  Output arguments
!    - wnear: real *8(nquad)
!        The desired near field quadrature for
!        S_{0}
!               
!

      implicit none 
      integer, intent(in) :: npatches, norders(npatches), npts, nquad
      integer, intent(in) :: ixyzs(npatches+1), iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts), eps
      real *8, intent(in) :: rfac0
      integer, intent(in) :: iquadtype
      integer, intent(in) :: nnz
      integer, intent(in) :: ndtarg, ntarg
      real *8, intent(in) :: targs(ndtarg, ntarg)
      integer, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      integer, intent(in) :: row_ptr(ntarg+1), col_ind(nnz), iquad(nnz+1)
      real *8, intent(out) :: wnear(nquad)

      real *8 dpars(2)

      dpars(1) = 1.0d0
      dpars(2) = 0.0d0

      call getnearquad_lap_comb_dir_eval(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, &
        targs, ipatch_id, uvs_targ, eps, dpars, iquadtype, nnz, &
        row_ptr, col_ind, iquad, rfac0, nquad, wnear)

      return
      end
!      
!
!
!
!
      subroutine lap_s_neu_eval(npatches, norders, ixyzs, &
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
!  Laplace single layer potential 
!
!  Representation:
!    u = \mathcal{S}_{0}[\sigma] 
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
!        iptype = 12, quadrangular patch discretized with Chebyshev nodes
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
!    - ndtarg: integer
!        leading dimension of target information array
!    - ntarg: integer
!        number of targets
!    - targs: real *8(ndtarg,ntarg)
!        target information 
!    - ipatch_id: integer(ntarg)
!        id of patch of target i, id = -1, if target is off-surface
!    - uvs_targ: double precision (2,ntarg)
!        local uv coordinates on patch if on surface, otherwise
!        set to 0 by default
!    - eps: double precision
!        precision requested
!     - sigma: double precision(npts)
!         density for layer potential
!
!  Output arguments
!    - pot: double precision(ntarg)
!        layer potential evaluated at the target points
!
!-----------------------------------
!
      implicit none
      integer, intent(in) :: npatches, npts
      integer, intent(in) :: ndtarg, ntarg
      integer, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts), eps
      real *8, intent(in) :: targs(ndtarg,ntarg)
      real *8, intent(in) :: sigma(npts)
      integer, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)

      real *8, intent(out) :: pot(ntarg)

      real *8 dpars(2)


      dpars(1) = 1.0d0
      dpars(2) = 0.0d0

      call lap_comb_dir_eval(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, dpars, sigma, pot)

      return
      end
!
!
!
!
!
      subroutine lap_s_neu_eval_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        eps, ndd, dpars, ndz, zpars, ndi, ipars, &
        nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
        nptso, ixyzso, srcover, whtsover, lwork, work, idensflag, &
        ndim_s, sigma, ipotflag, ndim_p, pot)
!
!
!  This subroutine computes the potential u
!  for the representation:
!
!  u = S_{0}[\sigma] 
!
!  where the near field is precomputed and stored
!  in the row sparse compressed format.
!        
!  The fmm is used to accelerate the far-field and 
!  near-field interactions are handled via precomputed quadrature
!
!  Using add and subtract - no need to call tree and set fmm parameters
!  can directly call existing fmm library
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
!        iptype = 12, quadrangular patch discretized with Chebyshev nodes
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
!    - ndtarg: integer
!        leading dimension of target information array
!    - ntarg: integer
!        number of targets
!    - targs: real *8(ndtarg,ntarg)
!        target information 
!    - eps: real *8
!        precision requested
!    - ndd: integer
!        number of real parameters defining the kernel/
!        integral representation (unused in this routine)
!    - dpars: real *8(ndd)
!        real parameters defining the kernel/
!        integral representation (unused in this routine) 
!    - ndz: integer
!        number of complex parameters defining the kernel/
!        integral representation (unused in this routine)
!    - zpars: complex *16(ndz)
!        complex parameters defining the kernel/
!        integral represnetation (unused in this routine)
!    - ndi: integer
!        number of integer parameters defining the kernel/
!        integral representation (unused in this routine)
!    - ipars: integer(ndi)
!        integer parameters defining the kernel/
!        integral represnetation (unused in this routine)
!    - nnz: integer
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer(ntarg+1)
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
!    - nquad: integer
!        number of near field entries corresponding to each source target
!        pair
!    - nker: integer
!        number of kernels in quadrature correction, must be 1
!    - wnear: real *8(nker, nquad)
!        precomputed quadrature corrections            
!    - novers: integer(npatches)
!        order of discretization for oversampled sources and
!        density
!    - ixyzso: integer(npatches+1)
!        ixyzso(i) denotes the starting location in srcover,
!        corresponding to patch i
!    - nptso: integer
!        total number of oversampled points
!    - srcover: real *8 (12,nptso)
!        oversampled set of source information
!    - whtsover: real *8 (nptso)
!        smooth quadrature weights at oversampled nodes
!    - lwork: integer
!        size of work array (unused in this routine)
!    - work: real *8(lwork)
!        work array (unused in this routine)
!    - idensflag: integer
!        Flag for types of denisties (unused in this case)
!    - ndim_s: integer
!        number of densities per point on the surface,
!        must be 1 for this routine
!    - sigma: real *8(ndim_s, npts)
!        The density sigma
!    - ipotflag: integer
!        Flag for determining output type 
!        ipotflag = 1, only potential is returned
!        ipotflag = 2, potential + gradients are
!        returned (currently unsupported)       
!    - ndim_p: integer
!        number of potentials per point
!        ndim_p = 1, if ipotflag = 1
!        ndim_p = 4, if ipotflag = 2        
!  Output arguments:
!    - pot: real *8 (ndim_p, ntarg)
!        u above
!

      implicit none
      integer, intent(in) :: npatches, npts
      integer, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
        
      integer, intent(in) :: ndtarg, ntarg
      real *8, intent(in) :: targs(ndtarg,ntarg)

      real *8, intent(in) :: eps
        
      integer, intent(in) :: ndd, ndz, ndi
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      integer, intent(in) :: ipars(ndi)
  
      integer, intent(in) :: nnz, nquad
      integer, intent(in) :: row_ptr(ntarg+1), col_ind(nnz)
      integer, intent(in) :: iquad(nnz+1)
        
      integer, intent(in) :: nker
      real *8, intent(in) :: wnear(nker,nquad)
  
      integer, intent(in) :: nptso
      integer, intent(in) :: ixyzso(npatches+1), novers(npatches)
      real *8, intent(in) :: srcover(12,nptso), whtsover(nptso)
        
      integer, intent(in) :: lwork
      real *8, intent(in) :: work(lwork)
  
      integer, intent(in) :: ndim_s, ndim_p
      integer, intent(in) :: idensflag, ipotflag
      
      real *8, intent(in) :: sigma(npts)
      
      real *8, intent(out) :: pot(ntarg)

      real *8 dpars_use(2)

      dpars_use(1) = 1.0d0
      dpars_use(2) = 0.0d0

      call lap_comb_dir_eval_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        eps, ndd, dpars_use, ndz, zpars, ndi, ipars, &
        nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
        nptso, ixyzso, srcover, whtsover, lwork, work, idensflag, &
        ndim_s, sigma, ipotflag, ndim_p, pot)

      return
      end

