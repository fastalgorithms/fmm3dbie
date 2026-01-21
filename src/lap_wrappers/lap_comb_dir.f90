!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!   Combined field representation for Laplace dirichlet problem
!
!  PDE:
!    \Delta u = 0
!
!  Boundary conditions:
!    u = f
!
!  Representation:
!    u = (\alpha S_{0}  + \beta D_{0})[\sigma]
!
!  Integral equations obtained by imposing:
!    u = f
!
!  and is given by:
!    d \sigma + (\alpha S_{0}  + \beta D_{0})[\sigma] = f
!
!  where:
!    d = -1/2 beta for interior problem
!      =  1/2 beta for exterior problem
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  User callable routines:
!    - lap_comb_dir_solver: Given data f,and parameters \alpha, \beta, 
!        this routine returns the solution \sigma
!
!    - lap_comb_dir_eval: Given \sigma, and parameters \alpha, \beta, 
!        evaluates the solution at a collection of targets
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  Advanced interfaces:
!    - getnearquad_lap_comb_dir: compute the near quadrature correction
!        for constructing the on-surface integral equation for the
!        Dirichlet data corresponding to the combined field
!        representation with user-provided near-information prescribed
!        in row-sparse compressed format
!
!    - getnearquad_lap_comb_dir_eval: compute the near quadrature
!        correction for target points which can be either on-surface or
!        off-surface with user-provided near-information prescribed in
!        row-sparse compressed format
!
!    - lpcomp_lap_comb_dir_addsub: apply the principal value part
!        of the integeral equation on surface. On input, user provides
!        precomputed near quadrature in row-sparse compressed format,
!        and oversampling information for handling the far part of the
!        computation
!
!    - lap_comb_dir_eval_addsub: compute the solution u at a
!        collection of targets (on-surface or off-surface), given \sigma.
!        On input, user provides precomputed near quadrature in
!        row-sparse compressed format and oversampling surface
!        information for the far-part
!
!    - lap_comb_dir_solver_guru: Guru solver routine, where user is
!        responsible for providing precomputed near quadrature
!        information in row-sparse compressed format and oversampling
!        surface information for the far-part
!
!  Other routines in the file, these are currently in beta mode
!  and should be used at your own risk:
!
!  There are two sets of fast direct solver routines
!  and neither of them currently used oversampling.
!  The fast direct solver routines are currently in beta
!  mode
!
!   - lap_comb_dir_fds_csc_mem: get memory requirements for initialization
!       routine for subsequent calls to fast direct solver
!
!   - lap_comb_dir_fds_csc_init: initialize various arrays to be later
!       used for fast direct solver
!
!   - lap_comb_dir_fds_csc_matgen: query entries of the combined field
!       representation matrix (input indices must be in column
!       sparse compressed format, and must be preceeded by a call
!       to lap_comb_dir_fds_init)
!
!
!

      subroutine getnearquad_lap_comb_dir(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, eps, dpars, &
        iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear)
!
!  This subroutine generates the near field quadrature
!  for the representation:
!
!  u = \alpha S_{0}[\sigma] + \beta D_{0}[\sigma]]    -   (1)
!
!  and returns quantities related to evaluating u on surface
!  at the surface discretization nodes.
!
!  If values at other points on the surface is desired then the 
!  data should be reinterpolated to those points from the
!  discretization nodes
!
!
!  On imposing the boundary condition, we get the following operator
!
!  u = \pm \beta \sigma/2 + \alpha S_{0}[\sigma] + \beta D_{0}[\sigma] 
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
!    - dpars: real *8(2)
!        kernel parameters (Referring to formula (1))
!        dpars(1) = alpha (single layer strength) 
!        dpars(2) = beta (double layer strength)
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
!        number of near field entries corresponding to each source target
!        pair
!
!  Output arguments
!    - wnear: real *8(nquad)
!        The desired near field quadrature for
!        \alpha S_{0} + \beta D_{0}
!               
!
      implicit none 
      integer *8, intent(in) :: npatches, norders(npatches), npts, nquad
      integer *8, intent(in) :: ixyzs(npatches+1), iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts), eps
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: iquadtype
      real *8, intent(in) :: dpars(2)
      integer *8, intent(in) :: nnz
      integer *8, intent(in) :: row_ptr(npts+1), col_ind(nnz), iquad(nnz+1)
      real *8, intent(out) :: wnear(nquad)

      integer *8 :: ndtarg, ntarg
      integer *8, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      integer *8 i

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

      call getnearquad_lap_comb_dir_eval(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, &
        srcvals, ipatch_id, uvs_targ, eps, dpars, iquadtype, nnz, &
        row_ptr, col_ind, iquad, rfac0, nquad, wnear)
      
      return
      end
!
!
!
!
!
      subroutine getnearquad_lap_comb_dir_eval(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, &
        targs, ipatch_id, uvs_targ, eps, dpars, iquadtype, nnz, &
        row_ptr, col_ind, iquad, rfac0, nquad, wnear)
!
!  This subroutine generates the near field quadrature
!  for the representation:
!
!  u = \alpha S_{0}[\sigma] + \beta D_{0}[\sigma]]    -   (1)
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
!    - eps: real *8
!        precision requested
!    - dpars: real *8(2)
!        kernel parameters (Referring to formula (1))
!        dpars(1) = alpha (single layer strength) 
!        dpars(2) = beta (double layer strength)
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
!    - wnear: real *8(nquad)
!        The desired near field quadrature for
!        \alpha S_{0} + \beta D_{0}
!               
!

      implicit none 
      integer *8, intent(in) :: npatches, norders(npatches), npts, nquad
      integer *8, intent(in) :: ixyzs(npatches+1), iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts), eps
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: iquadtype
      real *8, intent(in) :: dpars(2)
      integer *8, intent(in) :: nnz
      integer *8, intent(in) :: ndtarg, ntarg
      real *8, intent(in) :: targs(ndtarg, ntarg)
      integer *8, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      integer *8, intent(in) :: row_ptr(ntarg+1), col_ind(nnz), iquad(nnz+1)
      real *8, intent(out) :: wnear(nquad)

      real *8 alpha, beta
      integer *8 ndd, ndi, ndz
      complex *16 zpars(1)
      integer *8 ipars(1)

      integer *8 i, ipv

      procedure (), pointer :: fker
      external l3d_comb, l3d_slp, l3d_dlp
      
!
!        initialize the appropriate kernel function
!

      alpha = dpars(1)
      beta = dpars(2)

      ndd = 2
      ndi = 0
      ndz = 0
      fker => l3d_comb
      ipv = 1
      if(abs(alpha).ge.1.0d-16.and.abs(beta).lt.1.0d-16) then
        fker=>l3d_slp
        ipv = 0 
      else if(abs(alpha).lt.1.0d-16.and.abs(beta).ge.1.0d-16) then
        fker=>l3d_dlp
      endif

      if(iquadtype.eq.1) then
        call dgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, &
          zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, &
          nquad, wnear)
      endif



      if(abs(alpha).ge.1.0d-16.and.abs(beta).lt.1.0d-16) then
!$OMP PARALLEL DO DEFAULT(SHARED)        
        do i=1,nquad
          wnear(i) = wnear(i)*alpha
        enddo
!$OMP END PARALLEL DO        
      else if(abs(alpha).lt.1.0d-16.and.abs(beta).ge.1.0d-16) then
!$OMP PARALLEL DO DEFAULT(SHARED)        
        do i=1,nquad
          wnear(i) = wnear(i)*beta
        enddo
!$OMP END PARALLEL DO        
      endif

      return
      end
!
!
!
!
!
      subroutine lap_comb_dir_eval(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, dpars, sigma, pot)
!
!f2py intent(in) npatches,norders,ixyzs,iptype,npts,srccoefs,srcvals
!f2py intent(in) ndtarg,ntarg,targs,ipatch_id,uvs_targ,eps,dpars
!f2py intent(in) sigma
!f2py intent(out) pot
!
!
!------------------------------
!  This subroutine evaluates the layer potential for the 
!  Laplace combined field representation 
!
!  Representation:
!    u = (\alpha \mathcal{S}_{0} + \beta \mathcal{D}_{0})[\sigma]
!
!  Note: For targets on the boundary, this routine only computes
!  the principal value part, the identity term corresponding to the jump
!  in the layer potential is not included in the layer potential.
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
!    - dpars: real *8(2)
!        kernel parameters (Referring to formula above)
!        dpars(1) = $\alpha$
!        dpars(2) = $\beta$
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
      integer *8 norder ,npols
      integer *8, allocatable :: row_ptr(:), col_ind(:), iquad(:)
      real *8, allocatable :: wnear(:)

      real *8, allocatable :: srcover(:,:), wover(:)
      integer *8, allocatable :: ixyzso(:), novers(:)

      real *8, allocatable :: cms(:,:), rads(:), rad_near(:) 

      integer *8 i, j, jpatch, jquadstart, jstart

      complex *16 zpars
      integer *8 ipars

      real *8 ttot, done
      real *8 rfac, rfac0
      integer *8 iptype_avg, norder_avg
      integer *8 ikerorder, iquadtype, npts_over

      integer *8 ndd, ndz, ndi, nker, lwork, ndim_s, ndim_p
      integer *8 idensflag, ipotflag
      real *8 work(1)

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

      call findnearmem(cms, npatches, rad_near, ndtarg, targs, ntarg, &
        nnz)

      allocate(row_ptr(ntarg+1), col_ind(nnz))
      
      call findnear(cms, npatches, rad_near, ndtarg, targs, ntarg, &
        row_ptr, col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches, ixyzs, ntarg, nnz, row_ptr, &
        col_ind, iquad)

      ikerorder = -1
      if (abs(dpars(2)).gt.1.0d-16) ikerorder = 0

!
!    estimate oversampling for far-field, and oversample geometry
!

      allocate(novers(npatches), ixyzso(npatches+1))

      zpars = 0
      call get_far_order(eps, npatches, norders, ixyzs, iptype, cms, &
        rads, npts, srccoefs, ndtarg, ntarg, targs, ikerorder, zpars, &
        nnz, row_ptr, col_ind, rfac, novers, ixyzso)

      npts_over = ixyzso(npatches+1)-1

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

      call getnearquad_lap_comb_dir_eval(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, dpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)
!
!
!   compute layer potential
!
      ndd = 2
      ndz = 0
      ndi = 0
      nker = 1
      lwork = 0
      ndim_s = 1
      idensflag = 0
      ipotflag = 0
      ndim_p = 1
!
      call lap_comb_dir_eval_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        eps, ndd, dpars, ndz, zpars, ndi, ipars, nnz, row_ptr, &
        col_ind, iquad, nquad, nker, wnear, novers, npts_over, &
        ixyzso, srcover, wover, lwork, work, idensflag, ndim_s, &
        sigma, ipotflag, ndim_p, pot)

      return
      end
!
!
!  
!
!
      subroutine lpcomp_lap_comb_dir_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, &
        novers, nptso, ixyzso, srcover, whtsover, lwork, work, ndim, &
        sigma, pot)
!
!  This subroutine evaluates the Dirichlet data corresponding to
!  the following integral representation:
!  
!  u = \alpha S_{0}[\sigma] + \beta*D_{0}[\sigma]
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
!        integral representation (must be 2)
!    - dpars: real *8(ndd)
!        real parameters defining the kernel/
!        integral representation 
!        * dpars(1) = alpha (single layer strength) 
!        * dpars(2) = beta (double layer strength)
!    - ndz: integer *8
!        number of complex parameters defining the kernel/
!        integral representation (unused in this routine)
!    - zpars: complex *16(ndz)
!        complex parameters defining the kernel/
!        integral represnetation (unused in this routine)
!    - ndi: integer *8
!        number of integer parameters defining the kernel/
!        integral representation (unused in this routine)
!    - ipars: integer *8(ndi)
!        integer parameters defining the kernel/
!        integral represnetation (unused in this routine)
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
!        number of kernels in quadrature correction, must be 1
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
!    - lwork: integer *8
!        size of work array (unused in this routine)
!    - work: real *8(lwork)
!        work array (unused in this routine)
!    - ndim: integer *8
!        number of densities per point on the surface,
!        must be 1 for this routine
!    - sigma: real *8(npts)
!        The density sigma above                                        
!
!  Output arguments:
!    - pot: real *8(npts)
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
    
      real *8, intent(in) :: sigma(npts)
    
      real *8, intent(out) :: pot(npts)

      integer *8 ndtarg, idensflag, ipotflag, ndim_p

      ndtarg = 12
      idensflag = 0
      ipotflag = 0
      ndim_p = 1

      call lap_comb_dir_eval_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, npts, srcvals, &
        eps, ndd, dpars, ndz, zpars, ndi, ipars, &
        nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
        nptso, ixyzso, srcover, whtsover, lwork, work, idensflag, &
        ndim, sigma, ipotflag, ndim_p, pot)

      return
      end
!
!
!
!
      subroutine lap_comb_dir_eval_addsub(npatches, norders, ixyzs, &
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
!  u = \alpha S_{0}[\sigma] + \beta*D_{0}[\sigma] 
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
!        integral representation (must be 2)
!    - dpars: real *8(ndd)
!        real parameters defining the kernel/
!        integral representation 
!        * dpars(1) = alpha (single layer strength) 
!        * dpars(2) = beta (double layer strength)
!    - ndz: integer *8
!        number of complex parameters defining the kernel/
!        integral representation (unused in this routine)
!    - zpars: complex *16(ndz)
!        complex parameters defining the kernel/
!        integral represnetation (unused in this routine)
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
!        number of kernels in quadrature correction, must be 1
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
!    - lwork: integer *8
!        size of work array (unused in this routine)
!    - work: real *8(lwork)
!        work array (unused in this routine)
!    - idensflag: integer *8
!        Flag for types of denisties (unused in this case)
!    - ndim_s: integer *8
!        number of densities per point on the surface,
!        must be 1 for this routine
!    - sigma: real *8(ndim_s, npts)
!        The density sigma
!    - ipotflag: integer *8
!        Flag for determining output type 
!        ipotflag = 1, only potential is returned
!        ipotflag = 2, potential + gradients are
!        returned (currently unsupported)       
!    - ndim_p: integer *8
!        number of potentials per point
!        ndim_p = 1, if ipotflag = 1
!        ndim_p = 4, if ipotflag = 2        
!  Output arguments:
!    - pot: real *8 (ndim_p, ntarg)
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
      real *8, intent(in) :: wnear(nker,nquad)
  
      integer *8, intent(in) :: nptso
      integer *8, intent(in) :: ixyzso(npatches+1), novers(npatches)
      real *8, intent(in) :: srcover(12,nptso), whtsover(nptso)
        
      integer *8, intent(in) :: lwork
      real *8, intent(in) :: work(lwork)
  
      integer *8, intent(in) :: ndim_s, ndim_p
      integer *8, intent(in) :: idensflag, ipotflag
      
      real *8, intent(in) :: sigma(npts)
      
      real *8, intent(out) :: pot(ntarg)

      integer *8 norder,npols,nover,npolso
      real *8, allocatable :: potsort(:)

      real *8, allocatable :: sources(:,:), targvals(:,:)
      real *8, allocatable :: charges(:), dipvec(:,:), sigmaover(:)
      integer *8 ns, nt
      real *8 alpha, beta
      integer *8 ifcharge, ifdipole
      integer *8 ifpgh, ifpghtarg
      real *8 tmp(10), val

      integer *8 i, j, jpatch, jquadstart, jstart

      real *8 pottmp
      real *8, allocatable :: ctmp2(:), dtmp2(:,:)

      integer *8 nmax
      real *8 timeinfo(10), t1, t2, omp_get_wtime

      real *8, allocatable :: srctmp2(:,:)
      real *8 thresh, ra
      real *8 rr, rmin
      integer *8 nss, ii, l, npover

      integer *8 nd, ntarg0
      integer *8 ier, iper

      integer *8 ndtmp

      real *8 ttot, done

      parameter (nd=1,ntarg0=1)

      ns = nptso
      done = 1

!
!    estimate max number of sources in neear field of 
!    any target
!
      
      
      nmax = 0
      call get_near_corr_max(ntarg, row_ptr, nnz, col_ind, npatches, &
        ixyzso, nmax)
      allocate(srctmp2(3,nmax), ctmp2(nmax), dtmp2(3,nmax))
           
           
      ifpgh = 0
      ifpghtarg = 1
      allocate(sources(3,ns), targvals(3,ntarg))
      allocate(charges(ns), dipvec(3,ns))
      allocate(sigmaover(ns))

! 
!       oversample density
!

      call oversample_fun_surf(nd, npatches, norders, ixyzs, iptype, &
        npts, sigma, novers, ixyzso, ns, sigmaover)


!
!      set relevant parameters for the fmm
!
 
      alpha = dpars(1)
      beta = dpars(2)


!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)      
      do i=1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)

        charges(i) = sigmaover(i)*whtsover(i)*alpha
        dipvec(1,i) = sigmaover(i)*whtsover(i)*srcover(10,i)*beta
        dipvec(2,i) = sigmaover(i)*whtsover(i)*srcover(11,i)*beta
        dipvec(3,i) = sigmaover(i)*whtsover(i)*srcover(12,i)*beta
      enddo
!$OMP END PARALLEL DO      
      
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i = 1,ntarg
        targvals(1,i) = targs(1,i)
        targvals(2,i) = targs(2,i)
        targvals(3,i) = targs(3,i)

        pot(i) = 0
      enddo
!$OMP END PARALLEL DO      
      

      ifcharge = 1
      ifdipole = 1

      if(alpha.eq.0) ifcharge = 0
      if(beta.eq.0) ifdipole = 0

      iper = 0
      ier = 0

!
!
!       call the fmm
!

      call cpu_time(t1)
!$      t1 = omp_get_wtime()      
      call lfmm3d(nd, eps, ns, sources, ifcharge, charges, &
        ifdipole, dipvec, iper, ifpgh, tmp, tmp, tmp, ntarg, targvals, &
        ifpghtarg, pot, tmp, tmp, ier)
      call cpu_time(t2)
!$      t2 = omp_get_wtime()

      timeinfo(1) = t2-t1

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
!$OMP PRIVATE(jstart,pottmp,npols,l)
      do i = 1,ntarg
        do j = row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1) - ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            pot(i) = pot(i) + wnear(1,jquadstart+l-1)*sigma(jstart+l-1)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO



!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, jpatch, srctmp2) &
!$OMP PRIVATE(ctmp2, dtmp2, nss, l, jstart, ii, val, npover)
      do i = 1,ntarg
        nss = 0
        do j = row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l = ixyzso(jpatch),ixyzso(jpatch+1)-1
            nss = nss + 1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)

            if(ifcharge.eq.1) ctmp2(nss) = charges(l)
            if(ifdipole.eq.1) then
              dtmp2(1,nss) = dipvec(1,l)
              dtmp2(2,nss) = dipvec(2,l)
              dtmp2(3,nss) = dipvec(3,l)
            endif
          enddo
        enddo

        val = 0
        if(ifcharge.eq.1.and.ifdipole.eq.0) then
          call l3ddirectcp(nd, srctmp2, ctmp2, &
            nss, targvals(1,i), ntarg0, val, thresh)
        endif

        if(ifcharge.eq.0.and.ifdipole.eq.1) then
          call l3ddirectdp(nd, srctmp2, dtmp2, &
            nss, targvals(1,i), ntarg0, val, thresh)
        endif

        if(ifcharge.eq.1.and.ifdipole.eq.1) then
          call l3ddirectcdp(nd, srctmp2, ctmp2, dtmp2, &
            nss, targvals(1,i), ntarg0, val, thresh)
        endif
        pot(i) = pot(i) - val
      enddo
      
      call cpu_time(t2)
!$      t2 = omp_get_wtime()     

      
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
      subroutine lap_comb_dir_solver(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, dpars, numit, ifinout, &
        rhs, eps_gmres, niter, errs, rres, soln)

!
!
!  This subroutine solves the Laplace Dirichlet problem
!  on the exterior/interior of an object where the potential
!  is represented using a combined field integral representation.
!
!  This subroutine is the simple interface as opposed to the
!  _solver_guru routine which is called after initialization
!  in this routine.
!
!
!  Representation:
!    u = \alpha S_{0}[\sigma]+\beta*D_{0}[\sigma]
!
!  Boundary condition:
!    u = f
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
!    - dpars: real *8(2)
!        kernel parameters (Referring to formula (1))
!          * dpars(1) = \alpha 
!          * dpars(2) = \beta
!    - numit: integer *8
!        max number of gmres iterations
!    - ifinout: integer *8
!        ifinout = 0, interior problem
!        ifinout = 1, exterior problem
!    - rhs: real *8(npts)
!        Dirichlet data
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
!    - soln: real *8(npts)
!        density which solves the Dirichlet problem \sigma
!				 
      implicit none
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: ifinout
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps, eps_gmres
      real *8, intent(in) :: dpars(2)
      real *8, intent(in) :: rhs(npts)
      integer *8, intent(in) :: numit
      real *8, intent(out) :: soln(npts)
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
      real *8, allocatable :: wnear(:)

      real *8, allocatable :: srcover(:,:), wover(:)
      integer *8, allocatable :: ixyzso(:), novers(:)

      real *8, allocatable :: cms(:,:), rads(:), rad_near(:) 

      integer *8 i, j, jpatch, jquadstart, jstart

      integer *8 ipars
      complex *16 zpars
      real *8 timeinfo(10), t1, t2, omp_get_wtime


      real *8 ttot, done
      real *8 rfac, rfac0
      integer *8 iptype_avg, norder_avg
      integer *8 ikerorder, iquadtype, npts_over

!
!
!       gmres variables
!
      real *8 did, dtmp
      complex *16 ztmp
      integer *8 nker


      done = 1
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

      call getnearquad_lap_comb_dir(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, eps, dpars, iquadtype, &
        nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear)
      call cpu_time(t2)
!$      t2 = omp_get_wtime()     

      print *, "done generating near quadrature, now starting gmres"
      
      nker = 1
      call lap_comb_dir_solver_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, dpars, numit, ifinout, &
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

      subroutine lap_comb_dir_solver_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, dpars, numit, ifinout, &
        rhs, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
        nptso, ixyzso, srcover, whtsover, eps_gmres, niter, &
        errs, rres, soln)
!
!
!  This subroutine solves the Laplace Dirichlet problem
!  on the exterior of an object where the potential
!  is represented using the combined field integral representation.
!
!
!  Representation:
!    u = \alpha S_{0}[\sigma] + \beta*D_{0}[\sigma]
!
!  Boundary condition:
!    u = f
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
!    - dpars: real *8(2)
!        kernel parameters (Referring to formula (1))
!          * dpars(1) = \alpha 
!          * dpars(2) = \beta 
!    - numit: integer *8
!        max number of gmres iterations
!    - ifinout: integer *8
!        ifinout = 0, interior problem
!        ifinout = 1, exterior problem
!    - rhs: real *8(npts)
!        Dirichlet data
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
!        number of kernels in quadrature correction
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
!    - soln: real *8(npts)
!        density which solves the neumann problem \sigma
!				 

      implicit none
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps, eps_gmres
      integer *8, intent(in) :: ifinout
      real *8, intent(in) :: dpars(2)
      real *8, intent(in) :: rhs(npts)
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
      
  
      real *8, intent(out) :: soln(npts)

      real *8 did

      procedure (), pointer :: fker
      external lpcomp_lap_comb_dir_addsub

      integer *8 ndd, ndi, ndz, lwork, ndim
      real *8 work
      integer *8 ipars, nkertmp
      complex *16 zpars

      integer *8 ndtarg
      complex *16 ima
      data ima/(0.0d0,1.0d0)/
      real *8, allocatable :: wts(:)

      complex *16 zpars_use(3)

      did = -(-1)**(ifinout)*dpars(2)/2
      fker => lpcomp_lap_comb_dir_addsub

      ndd = 2
      ndi = 0
      ndz = 0

      lwork = 0
      ndim = 1
      allocate(wts(npts))

      call get_qwts(npatches, norders, ixyzs, iptype, npts, &
        srcvals, wts)
!

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
!
!
      subroutine lpcomp_lap_comb_dir_setsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        eps, dpars, nnz, row_ptr, col_ind, iquad, nquad, wnear, &
        sigma, novers, nptso, ixyzso, srcover, whtsover, pot)
!
!
!      this subroutine evaluates the layer potential for
!      the representation u = (\alpha S_{0} + \beta D_{0}) 
!      where the near field is precomputed and stored
!      in the row sparse compressed format.
!
!
!     The fmm is used to accelerate the far-field and 
!     near-field interactions are handled via precomputed quadrature
!
!
!       input:
!         npatches - integer *8
!            number of patches
!
!         norders- integer *8(npatches)
!            order of discretization on each patch 
!
!         ixyzs - integer *8(npatches+1)
!            ixyzs(i) denotes the starting location in srccoefs,
!               and srcvals array corresponding to patch i
!   
!         iptype - integer *8(npatches)
!            type of patch
!             iptype = 1, triangular patch discretized using RV nodes
!
!         npts - integer *8
!            total number of discretization points on the boundary
! 
!         srccoefs - real *8 (9,npts)
!            koornwinder expansion coefficients of xyz, dxyz/du,
!            and dxyz/dv on each patch. 
!            For each point srccoefs(1:3,i) is xyz info
!                           srccoefs(4:6,i) is dxyz/du info
!                           srccoefs(7:9,i) is dxyz/dv info
!
!         srcvals - real *8 (12,npts)
!             xyz(u,v) and derivative info sampled at the 
!             discretization nodes on the surface
!             srcvals(1:3,i) - xyz info
!             srcvals(4:6,i) - dxyz/du info
!             srcvals(7:9,i) - dxyz/dv info
! 
!         ndtarg - integer *8
!            leading dimension of target array
!        
!         ntarg - integer *8
!            number of targets
!
!         targs - real *8 (ndtarg,ntarg)
!            target information
!
!          eps - real *8
!             precision requested
!
!          dpars - real *8 (2)
!              kernel parameters (Referring to formula (1))
!              dpars(2) = alpha
!              dpars(3) = beta
!
!           nnz - integer *8
!             number of source patch-> target interactions in the near field
! 
!           row_ptr - integer *8(ntarg+1)
!              row_ptr(i) is the pointer
!              to col_ind array where list of relevant source patches
!              for target i start
!
!           col_ind - integer *8 (nnz)
!               list of source patches relevant for all targets, sorted
!               by the target number
!
!           iquad - integer *8(nnz+1)
!               location in wnear array where quadrature for col_ind(i)
!               starts
!
!           nquad - integer *8
!               number of entries in wnear
!
!           wnear - real *8(nquad)
!               the near field quadrature correction
!
!           sigma - real *8(npts)
!               density for layer potential
!
!           novers - integer *8(npatches)
!              order of discretization for oversampled sources and
!               density
!
!         ixyzso - integer *8(npatches+1)
!            ixyzso(i) denotes the starting location in srcover,
!               corresponding to patch i
!   
!           nptso - integer *8
!              total number of oversampled points
!
!           srcover - real *8 (12,nptso)
!              oversampled set of source information
!
!           whtsover - real *8 (nptso)
!             smooth quadrature weights at oversampled nodes
!
!           
!               
!
      implicit none
      integer *8 npatches,norder,npols,npts
      integer *8 ndtarg,ntarg
      integer *8 norders(npatches),ixyzs(npatches+1)
      integer *8 ixyzso(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps
      real *8 targs(ndtarg,ntarg)
      real *8 dpars(2)
      integer *8 nnz,row_ptr(ntarg+1),col_ind(nnz),nquad
      integer *8 iquad(nnz+1)
      real *8 wnear(nquad),sigma(npts)
      integer *8 novers(npatches+1)
      integer *8 nover,npolso,nptso
      real *8 srcover(12,nptso),whtsover(nptso)
      real *8 pot(ntarg)
      real *8, allocatable :: potsort(:)

      real *8, allocatable :: sources(:,:),targvals(:,:)
      real *8, allocatable :: charges(:),dipvec(:,:),sigmaover(:)
      integer *8, allocatable :: iboxtarg(:),iboxsrc(:)
      integer *8 ns,nt
      real *8 alpha,beta
      integer *8 ifcharge,ifdipole
      integer *8 ifpgh,ifpghtarg
      real *8 tmp(10),val


      integer *8 i,j,jpatch,jquadstart,jstart


      integer *8 ifaddsub

      integer *8 ltree,ipointer(8)
      integer *8, allocatable :: itree(:)
      integer *8, allocatable :: il1(:),il2(:),ilint(:),il1m2(:),il2m1(:)
      real *8, allocatable :: boxsize(:),centers(:,:)
      integer *8, allocatable :: isrcse(:,:),isrcper(:)
      integer *8, allocatable :: itargse(:,:),itargper(:)

      integer *8, allocatable :: nlist1(:),list1(:,:)
      integer *8, allocatable :: nlist2(:),list2(:,:)
      integer *8, allocatable :: nlist3(:),list3(:,:)
      integer *8, allocatable :: nlist4(:),list4(:,:)

      real *8 expc(3)
      integer *8 ibox,nexpc,idivflag,iert,ifnear,ii,isend,isep,isource
      integer *8 isstart,itarg,itend,itstart,itt,jbox,jpt,mhung,mnbors
      integer *8 iss,l,lpt,mnlist1,mnlist2,mnlist3,mnlist4
      integer *8 n1m2,n2m1,nadd,nbmax,nboxes,nchild,ndiv,nl2,nlevels
      integer *8 nlmax,npover,nl1,ntj
      integer *8 nlmin,iper,ier,ifunif
      
      real *8 ddot,pottmp
      real *8, allocatable :: ctmp1(:),ctmp2(:),dtmp1(:,:),dtmp2(:,:)
      real *8 radexp,epsfmm

      integer *8 ipars
      complex *16 zpars
      real *8 timeinfo(10),t1,t2,omp_get_wtime,timeinfo_fmm(6)

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp1(:,:),srctmp2(:,:)
      real *8 thresh,ra
      integer *8 nss

      integer *8 nd,ntarg0
      integer *8 int8_1

      real *8 ttot,done

      parameter (nd=1,ntarg0=1)

      int8_1 = 1
      ns = nptso
      done = 1

           
      ifpgh = 0
      ifpghtarg = 1
      allocate(sources(3,ns),targvals(3,ntarg))
      allocate(charges(ns),dipvec(3,ns))
      allocate(sigmaover(ns))

! 
!       oversample density
!

      call oversample_fun_surf(int8_1, npatches, norders, ixyzs, iptype, &
         npts,sigma,novers,ixyzso,ns,sigmaover)
!
!      set relevatn parameters for the fmm
!
      alpha = dpars(1)
      beta = dpars(2)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)      
      do i=1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)

        charges(i) = sigmaover(i)*whtsover(i)*alpha
        dipvec(1,i) = sigmaover(i)*whtsover(i)*srcover(10,i)*beta
        dipvec(2,i) = sigmaover(i)*whtsover(i)*srcover(11,i)*beta
        dipvec(3,i) = sigmaover(i)*whtsover(i)*srcover(12,i)*beta
      enddo
!$OMP END PARALLEL DO      

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,ntarg
        targvals(1,i) = targs(1,i)
        targvals(2,i) = targs(2,i)
        targvals(3,i) = targs(3,i)
      enddo
!$OMP END PARALLEL DO      

      ifcharge = 1
      ifdipole = 1

      if(alpha.eq.0) ifcharge = 0
      if(beta.eq.0) ifdipole = 0

!
!       setup tree
!

      isep = 1
      nlmax = 51
      nbmax = 0
      nlevels = 0
      nboxes = 0
      ltree = 0

      idivflag = 0
      mnlist1 = 0
      mnlist2 = 0
      mnlist3 = 0
      mnlist4 = 0


      call lndiv(eps,ns,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg, &
        ndiv,idivflag) 
!
!c      set tree flags
! 
       nlmax = 51
       nlevels = 0
       nboxes = 0
       ltree = 0
       nlmin = 0
       ifunif = 0
       iper = 0


!
!c     memory management code for contructing level restricted tree
      call pts_tree_mem(sources,ns,targvals,ntarg,idivflag,ndiv, &
        nlmin,nlmax,iper,ifunif,nlevels,nboxes,ltree)
      
       allocate(itree(ltree))
       allocate(boxsize(0:nlevels))
       allocate(centers(3,nboxes))

!       Call tree code
      call pts_tree_build(sources,ns,targvals,ntarg,idivflag,ndiv, &
        nlmin,nlmax,iper,ifunif,nlevels,nboxes,ltree,itree,ipointer, &
        centers,boxsize)
      
      

      allocate(isrcse(2,nboxes),itargse(2,nboxes))
      allocate(isrcper(ns),itargper(ntarg))

      call pts_tree_sort(ns,sources,itree,ltree,nboxes,nlevels, &
        ipointer,centers,isrcper,isrcse)

      call pts_tree_sort(ntarg,targvals,itree,ltree,nboxes,nlevels, &
        ipointer,centers,itargper,itargse)

      ifnear = 0

      mnbors = 27
      isep = 1

      call computemnlists(nlevels,nboxes,itree(ipointer(1)),boxsize, &
        centers,itree(ipointer(3)),itree(ipointer(4)), &
        itree(ipointer(5)),isep,itree(ipointer(6)),mnbors, &
        itree(ipointer(7)),iper,mnlist1,mnlist2,mnlist3,mnlist4)
      
      allocate(list1(mnlist1,nboxes),nlist1(nboxes))
      allocate(list2(mnlist2,nboxes),nlist2(nboxes))
      allocate(list3(mnlist3,nboxes),nlist3(nboxes))
      allocate(list4(mnlist4,nboxes),nlist4(nboxes))

      call computelists(nlevels,nboxes,itree(ipointer(1)),boxsize, &
        centers,itree(ipointer(3)),itree(ipointer(4)), &
        itree(ipointer(5)),isep,itree(ipointer(6)),mnbors, &
        itree(ipointer(7)),iper,nlist1,mnlist1,list1,nlist2, &
        mnlist2,list2,nlist3,mnlist3,list3, &
        nlist4,mnlist4,list4)



!
!
!       call the fmm
!

      ier = 0
      iper = 0
      call cpu_time(t1)
!$      t1 = omp_get_wtime()      
      call lfmm3d_ndiv(nd,eps,ns,sources,ifcharge,charges, &
        ifdipole,dipvec,iper,ifpgh,tmp,tmp,tmp,ntarg,targvals,ifpghtarg, &
        pot,tmp,tmp,ndiv,idivflag,ifnear,timeinfo_fmm,ier)
      call cpu_time(t2)
!$      t2 = omp_get_wtime()

      timeinfo(1) = t2-t1

      
!
!
!       add in precomputed quadrature
!

      thresh = 2.0d0**(-51)*boxsize(0)
      call cpu_time(t1)
!$      t1 = omp_get_wtime()

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart) &
!$OMP PRIVATE(jstart,pottmp,npols,l)
      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            pot(i) = pot(i) + wnear(jquadstart+l-1)*sigma(jstart+l-1)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO



!
!
!    work with sorted potentials and unsort them again later
!
      allocate(potsort(ntarg))
      call dreorderf(int8_1,ntarg,pot,potsort,itargper)



!
!    subtract  precomputed near quadrature /setminus list1 
!       also needs to go from pts (targs) -> pts (sources)
! 
!
!    il1 - list of sources in the near field of a target (A)
!    il2 - list of sources in the list1 of the target from fmm
!        perspective (B)
!    il1m2 = A \cap (A \cap B)^{c}
!    il2m1 = B \cap (A \cap B)^{c}
!

     
      allocate(il2(ndiv*mnlist1),il2m1(ndiv*mnlist1))
      allocate(ctmp2(ndiv*mnlist1),dtmp2(3,ndiv*mnlist1))
      allocate(srctmp2(3,ndiv*mnlist1))
      allocate(srctmp1(3,ndiv*mnlist1))
      allocate(ctmp1(ndiv*mnlist1),dtmp1(3,ndiv*mnlist1))
      allocate(il1(ndiv*mnlist1),il1m2(ndiv*mnlist1))

  

      call cpu_time(t1)
!$      t1 = omp_get_wtime()     

!$OMP PARALLEL DO DEFAULT(SHARED) & 
!$OMP PRIVATE(ibox,nchild,nl2) &
!$OMP PRIVATE(i,jbox,isstart,isend,j,isource,il2) &
!$OMP PRIVATE(itstart,itend,itt,itarg,nl1,il1,il1m2,il2m1) &
!$OMP PRIVATE(jpatch,l,jpt,lpt,n1m2,n2m1,ii,val,npover) &
!$OMP PRIVATE(ctmp1,ctmp2,dtmp1,dtmp2,srctmp1,srctmp2) &
!$OMP SCHEDULE(DYNAMIC)
      do ibox = 1,nboxes
        nchild = itree(ipointer(4)+ibox-1)
        if(nchild.eq.0) then

!
!     populate il2
!
          nl2 = 0
          do i=1,nlist1(ibox)
            jbox = list1(i,ibox) 
            isstart = isrcse(1,jbox) 
            isend = isrcse(2,jbox)
            do j=isstart,isend
              isource = isrcper(j) 
              nl2 = nl2 + 1
              il2(nl2) = isource
            enddo
          enddo


!
!    end of populating il2.
!    
!    now loop over targets in this box
!
          itstart = itargse(1,ibox) 
          itend = itargse(2,ibox)
          do itt = itstart,itend
            itarg = itargper(itt) 
            
            nl1 = 0
            do j=row_ptr(itarg),row_ptr(itarg+1)-1
              jpatch = col_ind(j)
              nl1 = nl1 + ixyzso(jpatch+1)-ixyzso(jpatch)
            enddo

!
!    populate il1 
!

            lpt = 0
            do j = row_ptr(itarg),row_ptr(itarg+1)-1
              jpatch = col_ind(j)
              npover = ixyzso(jpatch+1)-ixyzso(jpatch)
              do l=1,npover
                jpt = ixyzso(jpatch)+l-1
                lpt = lpt + 1
                il1(lpt) = jpt
              enddo
            enddo
!
!   end of populating il1. now perform various set subtractions
!
            n1m2 = 0
            n2m1 = 0
            call setsub(il1,nl1,il2,nl2,il1m2,n1m2,il2m1,n2m1)


!
!   subtract off il1m2
!
!   gather step
!
            do i=1,n1m2
              ii = il1m2(i)
              srctmp1(1,i) = srcover(1,ii)
              srctmp1(2,i) = srcover(2,ii)
              srctmp1(3,i) = srcover(3,ii)
            enddo
            if(ifcharge.eq.1) then
              do i=1,n1m2
                ii = il1m2(i)
                ctmp1(i) = charges(ii)
              enddo
            endif

            if(ifdipole.eq.1) then
              do i=1,n1m2
                ii = il1m2(i)
                dtmp1(1,i) = dipvec(1,ii)
                dtmp1(2,i) = dipvec(2,ii)
                dtmp1(3,i) = dipvec(3,ii)
              enddo
            endif

            val = 0
            if(ifcharge.eq.1.and.ifdipole.eq.0) then
              call l3ddirectcp(nd,srctmp1,ctmp1, &
                n1m2,targvals(1,itarg),ntarg0,val,thresh)
            endif

            if(ifcharge.eq.0.and.ifdipole.eq.1) then
              call l3ddirectdp(nd,srctmp1,dtmp1, &
                n1m2,targvals(1,itarg),ntarg0,val,thresh)
            endif

            if(ifcharge.eq.1.and.ifdipole.eq.1) then
              call l3ddirectcdp(nd,srctmp1,ctmp1,dtmp1, &
                n1m2,targvals(1,itarg),ntarg0,val,thresh)
            endif
!
!  scatter step
!
            potsort(itt) = potsort(itt) - val



!
!   add il2m1
!
!
!   gather step
!
            do i=1,n2m1
              ii = il2m1(i)
              srctmp2(1,i) = srcover(1,ii)
              srctmp2(2,i) = srcover(2,ii)
              srctmp2(3,i) = srcover(3,ii)
            enddo
            if(ifcharge.eq.1) then
              do i=1,n2m1
                ii = il2m1(i)
                ctmp2(i) = charges(ii)
              enddo
            endif

            if(ifdipole.eq.1) then
              do i=1,n2m1
                ii = il2m1(i)
                dtmp2(1,i) = dipvec(1,ii)
                dtmp2(2,i) = dipvec(2,ii)
                dtmp2(3,i) = dipvec(3,ii)
              enddo
            endif

            val = 0
            if(ifcharge.eq.1.and.ifdipole.eq.0) then
              call l3ddirectcp(nd,srctmp2,ctmp2, &
                n2m1,targvals(1,itarg),ntarg0,val,thresh)
            endif

            if(ifcharge.eq.0.and.ifdipole.eq.1) then
              call l3ddirectdp(nd,srctmp2,dtmp2, &
                n2m1,targvals(1,itarg),ntarg0,val,thresh)
            endif

            if(ifcharge.eq.1.and.ifdipole.eq.1) then
              call l3ddirectcdp(nd,srctmp2,ctmp2,dtmp2, &
                n2m1,targvals(1,itarg),ntarg0,val,thresh)
            endif
!
! scatter step
!
            potsort(itt) = potsort(itt) + val

          enddo
        endif
      enddo
!$OMP END PARALLEL DO      
      call cpu_time(t2)
!$      t2 = omp_get_wtime()      
      timeinfo(2) = t2-t1

      call dreorderi(int8_1,ntarg,potsort,pot,itargper)

      
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
      subroutine lap_comb_dir_fds_csc_mem(npatches,norders,ixyzs, &
        iptype,npts,srccoefs,srcvals,eps,dpars,nifds,nrfds,nzfds)
!
!f2py   intent(in) npatches,norders,ixyzs,iptype,npts
!f2py   intent(in) srccoefs,srcvals,eps,dpars
!f2py   intent(out) nifds,nrfds,nzfds

!
!----------------------------------
!  This subroutine estimates the memory requirements
!  for the precomputation routine of the fast direct solver
!  for solving Dirichlet problems for Laplace's equation
!  using the combined field integral equation:
!
!  .. math::
!
!     u = \alpha S_{0}[\sigma] + \beta D_{0} [\sigma]
! 
!  The precomputation routine computes an integer array (ifds)
!  a real array (rfds), and complex array (zfds)
! 
!  The following quantities will be computed during the 
!  precomputation phase of the fast direct solver
!
!  - ifds(1): npts_over
!  - ifds(2): nnz
!  - ifds(3): nquad
!  - ifds(4): nximat
!  - ifds(5:6+npts-1): row_ptr (for near quadrature info)
!  - ifds(6+npts:6+npts+nnz-1): col_ind
!  - ifds(6+npts+nnz:6+npts+2*nnz): iquad
!  - ifds(7+npts+2*nnz:7+npts+2*nnz+npatches-1): novers
!  - ifds(7+npts+2*nnz+npatches:7+npts+2*nnz+2*npatches): ixyzso
!  - ifds(8+npts+2*nnz+2*npatches:8+npts+2*nnz+3*npatches-1): iximat
!
!  - rfds(1:2): dpars
!  - rfds(3:12*npts_over+2): srcover
!  - rfds(12*npts_over+3:13*npts_over+2): wover
!  - rfds(13*npts_over+3:13*npts_over+nximat+2): ximats
!  - rfds(13*npts_over+nximat+3:13*npts_over+nximat+nquad+2): wnear
!
!  - zfds: unused 
!
!  Input arguments:
!  
!    - npatches: integer *8
!        number of patches
!    - norders: integer *8(npatches)
!        order of discretization of each patch
!    - ixyzs: integer *8(npatches+1)
!        starting location of points on patch i
!    - iptype: integer *8(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer *8
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!    - eps: double precision
!        precision requested
!    - dpars: double precision(2)
!        weights for single and double layer potential in representation
!        dpars(1) = \alpha
!        dpars(2) = \beta
!
!  Output arguments:
!  
!    - nifds: integer *8
!        size of integer array
!    - nrfds: integer *8
!        size of double precision array
!    - nzfds: integer *8
!        size of double complex array
!     
      implicit none
      integer *8, intent(in) :: npatches,norders(npatches)
      integer *8, intent(in) :: ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches),npts
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
      real *8, intent(in) :: eps
      real *8, intent(in) :: dpars(2)
      integer *8, intent(out) :: nifds,nrfds,nzfds

      integer *8 nnz
      real *8, allocatable :: targs(:,:)
      integer *8 iptype_avg,norder_avg
      integer *8 ntarg,ndtarg,ikerorder
      real *8 rfac,rfac0
      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)
      integer *8, allocatable :: iquad(:),row_ptr(:),col_ind(:)
      integer *8, allocatable :: novers(:),ixyzso(:)
      integer *8 npts_over,nquad,nximat
      

      complex *16 ztmp
      integer *8 i


!
!
!        setup targets as on surface discretization points
! 
      ndtarg = 3
      ntarg = npts
      allocate(targs(ndtarg,npts))

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        targs(1,i) = srcvals(1,i)
        targs(2,i) = srcvals(2,i)
        targs(3,i) = srcvals(3,i)
      enddo
!$OMP END PARALLEL DO   


!
!
!        this might need fixing
!
      iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
      norder_avg = floor(sum(norders)/(npatches+0.0d0))

      call get_rfacs(norder_avg,iptype_avg,rfac,rfac0)


      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))

      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts, & 
        srccoefs,cms,rads)

!$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,npatches
        rad_near(i) = rads(i)*rfac
      enddo
!$OMP END PARALLEL DO      

!
!    find near quadrature correction interactions
!
      call findnearmem(cms,npatches,rad_near,ndtarg,targs,npts,nnz)

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,ndtarg,targs,npts,row_ptr, &
        col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind, &
        iquad)

      ikerorder = -1
      if(abs(dpars(2)).gt.1.0d-16) ikerorder = 0


!
!    estimate oversampling for far-field, and oversample geometry
!

      allocate(novers(npatches),ixyzso(npatches+1))

      ztmp = 0
      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms, &
        rads,npts,srccoefs,ndtarg,npts,targs,ikerorder,ztmp, &
        nnz,row_ptr,col_ind,rfac,novers,ixyzso)


      

      npts_over = ixyzso(npatches+1)-1

      nquad = iquad(nnz+1)-1

      call get_nximat(npatches,ixyzs,ixyzso,nximat)

      nifds = 7+npts+2*nnz+3*npatches
      nrfds = 2+13*npts_over+nximat+nquad
      nzfds = 0
      

      return
      end
!
!
!
!
!

      subroutine lap_comb_dir_fds_csc_init(npatches,norders,ixyzs, &
        iptype,npts,srccoefs,srcvals,eps,dpars,nifds,ifds,nrfds, &
        rfds,nzfds,zfds)
!f2py intent(in) npatches,norders,ixyzs,iptype,npts
!f2py intent(in) srccoefs,srcvals,eps,dpars
!f2py intent(in) nifds,nrfds,nzfds
!f2py intent(out) ifds,rfds,zfds

!
!----------------------------------
!  This subroutine precomputes a few arrays
!  for subsequent calls to the fast direct solver
!  for solving Dirichlet problems for Laplace's equation
!  using the combined field integral equation:
!
!  .. math::
!
!     u = \alpha S_{0}[\sigma] + \beta D_{0} [\sigma]
! 
!  The precomputation routine computes an integer array (ifds)
!  a real array (rfds), and complex array (zfds)
! 
!  The following quantities will be computed during the 
!  precomputation phase of the fast direct solver
!
!  - ifds(1): npts_over
!  - ifds(2): nnz
!  - ifds(3): nquad
!  - ifds(4): nximat
!  - ifds(5:6+npts-1): row_ptr (for near quadrature info)
!  - ifds(6+npts:6+npts+nnz-1): col_ind
!  - ifds(6+npts+nnz:6+npts+2*nnz): iquad
!  - ifds(7+npts+2*nnz:7+npts+2*nnz+npatches-1): novers
!  - ifds(7+npts+2*nnz+npatches:7+npts+2*nnz+2*npatches): ixyzso
!  - ifds(8+npts+2*nnz+2*npatches:8+npts+2*nnz+3*npatches-1): iximat
!
!  - rfds(1:2): dpars
!  - rfds(3:12*npts_over+2): srcover
!  - rfds(12*npts_over+3:13*npts_over+2): wover
!  - rfds(13*npts_over+3:13*npts_over+nximat+2): ximats
!  - rfds(13*npts_over+nximat+3:13*npts_over+nximat+nquad+2): wnear
!
!  - zfds: unused 
!
!  Input arguments:
!  
!    - npatches: integer *8
!        number of patches
!    - norders: integer *8(npatches)
!        order of discretization of each patch
!    - ixyzs: integer *8(npatches+1)
!        starting location of points on patch i
!    - iptype: integer *8(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer *8
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!    - eps: double precision
!        precision requested
!    - dpars: double precision(2)
!        weights for single and double layer potential in representation
!        dpars(1) = \alpha
!        dpars(2) = \beta
!    - nifds: integer *8
!        size of integer array
!    - nrfds: integer *8
!        size of double precision array
!    - nzfds: integer *8
!        size of double complex array
!
!  Output arguments:
!
!    - ifds: integer *8(nifds)
!        precomputed integer array
!    - rfds: double precision(nifds)
!        precomputed double precision array
!    - zfds: double complex(nifds)
!        precomputed double complex array
!  
!     
      implicit none
      integer *8, intent(in) :: npatches,norders(npatches)
      integer *8, intent(in) :: ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches),npts
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
      real *8, intent(in) :: eps
      real *8, intent(in) :: dpars(2)
      integer *8, intent(in) :: nifds,nrfds,nzfds
      integer *8, intent(out) :: ifds(nifds)
      real *8, intent(out) :: rfds(nrfds)
      complex *16, intent(out) :: zfds(nzfds)

      integer *8 nnz
      complex *16 ztmp
      real *8, allocatable :: targs(:,:)
      real *8, allocatable :: uvs_targ(:,:)
      integer *8, allocatable :: ipatch_id(:)
      integer *8 iptype_avg,norder_avg
      integer *8 ntarg,ndtarg,ikerorder
      real *8 rfac,rfac0
      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)
      integer *8 npts_over,nquad

      integer *8 iquadtype,istart,iend,i

      integer *8 irow_ptr,icol_ind,iiquad,inovers,iixyzso,iximat
      integer *8 nximat

!
!
!        setup targets as on surface discretization points
! 
      ndtarg = 3
      ntarg = npts
      allocate(targs(ndtarg,npts))

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        targs(1,i) = srcvals(1,i)
        targs(2,i) = srcvals(2,i)
        targs(3,i) = srcvals(3,i)
      enddo
!$OMP END PARALLEL DO   


!
!
!        this might need fixing
!
      iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
      norder_avg = floor(sum(norders)/(npatches+0.0d0))

      call get_rfacs(norder_avg,iptype_avg,rfac,rfac0)


      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))

      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts, &
        srccoefs,cms,rads)

!$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,npatches
        rad_near(i) = rads(i)*rfac
      enddo
!$OMP END PARALLEL DO      

!
!    find near quadrature correction interactions
!
      call findnearmem(cms,npatches,rad_near,ndtarg,targs,npts,nnz)

      irow_ptr = 5
      icol_ind = irow_ptr+npts+1
      iiquad = icol_ind+nnz
      inovers = iiquad + nnz+1
      iixyzso = inovers+npatches
      iximat = iixyzso+npatches+1

      
      call findnear(cms,npatches,rad_near,ndtarg,targs,npts, &
        ifds(irow_ptr),ifds(icol_ind))

      call get_iquad_rsc(npatches,ixyzs,npts,nnz,ifds(irow_ptr), &
        ifds(icol_ind),ifds(iiquad))


      ikerorder = -1
      if(abs(dpars(2)).gt.1.0d-16) ikerorder = 0


!
!    estimate oversampling for far-field, and oversample geometry
!

      ztmp = 0.0d0
      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms, &
        rads,npts,srccoefs,ndtarg,npts,targs,ikerorder,ztmp, &
        nnz,ifds(irow_ptr),ifds(icol_ind),rfac,ifds(inovers), &
        ifds(iixyzso))

      npts_over = ifds(iixyzso+npatches)-1
      nquad = ifds(iiquad+nnz)-1
      
      rfds(1) = dpars(1)
      rfds(2) = dpars(2)

      call oversample_geom(npatches,norders,ixyzs,iptype,npts, & 
        srccoefs,srcvals,ifds(inovers),ifds(iixyzso),npts_over, &
        rfds(3))

      call get_qwts(npatches,ifds(inovers),ifds(iixyzso),iptype, &
        npts_over,rfds(3),rfds(12*npts_over+3))

      

      ifds(1) = npts_over
      ifds(2) = nnz
      ifds(3) = nquad

      call get_nximat(npatches,ixyzs,ifds(iixyzso),nximat)

      ifds(4) = nximat

      istart = 13*npts_over+3

      call get_ximats(npatches,iptype,norders,ixyzs,ifds(inovers), &
        ifds(iixyzso),nximat,rfds(istart),ifds(iximat))


      iquadtype = 1
      istart = 13*npts_over + nximat + 3

      call getnearquad_lap_comb_dir(npatches,norders, &
        ixyzs,iptype,npts,srccoefs,srcvals,eps,dpars,iquadtype, &
        nnz,ifds(irow_ptr),ifds(icol_ind),ifds(iiquad),rfac0,nquad, &
        rfds(istart))
      
      return
      end
!
!
!
!
!
!

      subroutine lap_comb_dir_fds_csc_matgen(npatches,norders,ixyzs, &
        iptype,npts,srccoefs,srcvals,eps,dpars,nifds,ifds,nrfds, &
        rfds,nzfds,zfds,nent_csc,col_ptr,row_ind,dmatent)
!f2py intent(in) npatches,norders,ixyzs,iptype,npts
!f2py intent(in) srccoefs,srcvals,eps,dpars
!f2py intent(in) nifds,nrfds,nzfds
!f2py intent(in) ifds,rfds,zfds
!f2py intent(in) nent_csc,col_ptr,row_ind
!f2py intent(out) dmatent
   
!----------------------------------
!  This subroutine generates matrix entries queried 
!  for solving Dirichlet problems for Laplace's equation
!  using the combined field integral equation:
!
!  .. math::
!
!     u = \alpha S_{0}[\sigma] + \beta D_{0} [\sigma]
! 
!  The precomputation routine computes an integer array (ifds)
!  a real array (rfds), and complex array (zfds)
! 
!  The following quantities will be computed during the 
!  precomputation phase of the fast direct solver
!
!  - ifds(1): npts_over
!  - ifds(2): nnz
!  - ifds(3): nquad
!  - ifds(4): nximat
!  - ifds(5:6+npts-1): row_ptr (for near quadrature info)
!  - ifds(6+npts:6+npts+nnz-1): col_ind
!  - ifds(6+npts+nnz:6+npts+2*nnz): iquad
!  - ifds(7+npts+2*nnz:7+npts+2*nnz+npatches-1): novers
!  - ifds(7+npts+2*nnz+npatches:7+npts+2*nnz+2*npatches): ixyzso
!  - ifds(8+npts+2*nnz+2*npatches:8+npts+2*nnz+3*npatches-1): iximat
!
!  - rfds(1:2): dpars
!  - rfds(3:12*npts_over+2): srcover
!  - rfds(12*npts_over+3:13*npts_over+2): wover
!  - rfds(13*npts_over+3:13*npts_over+nximat+2): ximats
!  - rfds(13*npts_over+nximat+3:13*npts_over+nximat+nquad+2): wnear
!
!  - zfds: unused
!
!  The list of input entries requested must be provided in 
!  column sparse compressed format. Use the conv_to_csc
!  utility to convert a list of (i,j) entries to column sparse
!  compressed representation.
!
!  Note: if double layer potential is used in representation,
!  then only principal value part of the layer potential
!  is returned. The term corresponding to the jump in layer potential
!  will have to be added manually.
!
!
!  Input arguments:
!  
!    - npatches: integer *8
!        number of patches
!    - norders: integer *8(npatches)
!        order of discretization of each patch
!    - ixyzs: integer *8(npatches+1)
!        starting location of points on patch i
!    - iptype: integer *8(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer *8
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!    - eps: double precision
!        precision requested
!    - dpars: double precision(2)
!        weights for single and double layer potential in representation
!        dpars(1) = \alpha
!        dpars(2) = \beta
!    - nifds: integer *8
!        size of integer array
!    - ifds: integer *8(nifds)
!        precomputed integer array
!    - nrfds: integer *8
!        size of double precision array
!    - rfds: double precision(nifds)
!        precomputed double precision array
!    - nzfds: integer *8
!        size of double complex array
!    - zfds: double complex(nifds)
!        precomputed double complex array
!    - nent_csc: integer *8
!        number of matrix entries requested
!    - col_ptr: integer *8(npts+1)
!        indicates location in row_ind where relevant entries for
!        column i begin
!    - row_ind: integer *8(nent_csc)
!        list of row indices. (row_ind(col_ptr(i):col_ptr(i+1)-1),i)
!        are all the matrix entries requested
!
!  Output arguments:
!
!    - dmatent: double precision(npts)
!        matrix entries requested 
!
!  
!     

      implicit none
      integer *8, intent(in) :: npatches,norders(npatches)
      integer *8, intent(in) :: ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches),npts
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
      real *8, intent(in) :: eps,dpars(2)
      integer *8, intent(in) :: nifds,nrfds,nzfds
      integer *8, intent(in) :: ifds(nifds)
      real *8, intent(in) :: rfds(nrfds)
      complex *16, intent(in) :: zfds(nzfds)
      integer *8, intent(in) :: nent_csc,col_ptr(npts+1),row_ind(nent_csc)
      real *8, intent(out) :: dmatent(nent_csc)

!
!        temporary variables
!
      integer *8 i,j,k,l,ipatch,npols
      integer *8 ndd,ndz,ndi
      integer *8 ilstart,ilend,istart,iend
      integer *8 irow_ptr,icol_ind,iiquad,inovers,iixyzso
      integer *8 npts_over,nnz,nquad

      integer *8, allocatable :: iuni(:),iuniind(:)
      integer *8, allocatable :: col_ptr_src(:),row_ind_src(:),iper(:)
      integer *8, allocatable :: aintb(:),iaintba(:),aintbc(:),iaintbc(:)
      real *8, allocatable :: wquadn(:,:),wquad(:,:)
      real *8, allocatable :: wquadf(:,:),wquadf2(:,:)


      real *8 alpha, beta
      complex *16 zk,zpars
      real *8, allocatable :: srcover(:,:),wtsover(:)
      integer *8 ipars,ipt,iquad,itind,iximat,ixist,j2,jind
      integer *8 jind0,juniind,n2,naintb,naintbc,nmax,nn,npolso
      integer *8 nuni,nximat,iwnear
      integer *8, allocatable :: iaintbb(:)
      integer *8 ifcharge,ifdipole

      real *8 dzero,done

      procedure (), pointer :: fker
      external l3d_slp, l3d_dlp, l3d_comb

      done = 1
      dzero = 0

!
!
!        initialize the appropriate kernel function
!
 
      
      alpha = rfds(1)
      beta = rfds(2)

      fker => l3d_comb
      ndd = 2
      ifcharge = 1
      ifdipole = 1
      if(abs(alpha).ge.1.0d-16.and.abs(beta).lt.1.0d-16) then
        fker=>l3d_slp
        ndd = 0
        ifdipole = 0
      else if(abs(alpha).lt.1.0d-16.and.abs(beta).ge.1.0d-16) then
        fker=>l3d_dlp
        ndd = 0
        ifcharge = 0
      endif


      ndi = 0
      ndz = 0

      npts_over = ifds(1)
      nnz = ifds(2)
      nquad = ifds(3)
      nximat = ifds(4)

      irow_ptr = 5
      icol_ind = irow_ptr+npts+1
      iiquad = icol_ind+nnz
      inovers = iiquad + nnz+1
      iixyzso = inovers+npatches
      iximat = iixyzso+npatches+1

      allocate(col_ptr_src(npatches+1),row_ind_src(nnz),iper(nnz))

      call rsc_to_csc(npatches,npts,nnz,ifds(irow_ptr),ifds(icol_ind), &
        col_ptr_src,row_ind_src,iper)

!
!    estimate max oversampling
!
      nmax = 0
      do ipatch=1,npatches
        npolso = ifds(iixyzso+ipatch)-ifds(iixyzso+ipatch-1)
        if(npolso.gt.nmax) nmax = npolso
      enddo

      allocate(srcover(12,nmax),wtsover(nmax))

      iwnear = 13*npts_over + nximat+2 
     

      do ipatch=1,npatches

!
!
!   combine all list of targets requested for current
!   patch and find unique list of targets
!
      
        istart = ixyzs(ipatch)
        iend = ixyzs(ipatch+1)-1
        npols= iend-istart+1
         
        ilstart = col_ptr(istart)
        ilend = col_ptr(iend+1)-1

        nn = ilend-ilstart+1

        allocate(iuni(nn),iuniind(nn))

        nuni = 0
        call get_iuni1(nn,row_ind(ilstart),nuni,iuni,iuniind)
        
        allocate(aintb(nuni),iaintba(nuni),aintbc(nuni),iaintbc(nuni))
        allocate(iaintbb(nuni))

        n2 = col_ptr_src(ipatch+1)-col_ptr_src(ipatch)

!
!
!    separate list of targets into near field targets for which 
!    quadrature is already computed and far-field targets for which
!    quadrature is required to be computed
!
        naintb = 0
        naintbc = 0

        call setdecomp(nuni,iuni,n2, &
          row_ind_src(col_ptr_src(ipatch)),naintb,aintb,iaintba, &
          iaintbb,naintbc,aintbc,iaintbc)

       
!
!    for the entries in aintb, the quadrature has already been computed
!    so that needs to be extracted and sent to appropriate entries
!  
        allocate(wquad(nuni,npols))
        do i=1,naintb
           jind0 = iaintbb(i)+col_ptr_src(ipatch)-1
           jind = iper(jind0)
           iquad = ifds(iiquad + jind-1)
           wquad(iaintba(i),:) = rfds(iwnear+iquad: &
             iwnear+iquad+npols-1)
        enddo

!
!        compute the oversampled quadrature
!
        istart = ifds(iixyzso+ipatch-1)
        iend = ifds(iixyzso+ipatch)-1
        npolso = iend-istart+1

        allocate(wquadf(npolso,naintbc))
        allocate(wquadf2(npols,naintbc))

!
!      extract srcover, wtsover, and targvals 
!

        do i=1,npolso
          do j=1,12
            srcover(j,i) = rfds(2+12*(istart+i-2)+j)
          enddo
          wtsover(i) = rfds(2+12*npts_over+istart+i-1)
        enddo


        do i=1,naintbc
          itind = aintbc(i)
          do j=1,npolso
            call fker(srcover(1,j),3,srcvals(1,itind),ndd,dpars,ndz, &
              zfds,ndi,ipars,wquadf(j,i))
            wquadf(j,i) = wquadf(j,i)*wtsover(j)
          enddo
        enddo

        if(ifdipole.eq.0) then
          do i=1,naintbc
            do j=1,npolso
              wquadf(j,i) = wquadf(j,i)*alpha
            enddo
          enddo
        endif

        if(ifcharge.eq.0) then
          do i=1,naintbc
            do j=1,npolso
              wquadf(j,i) = wquadf(j,i)*beta
            enddo
          enddo
        endif

!
!      now multiply wquadf by ximat
!
        ixist = ifds(iximat+ipatch-1) + 13*npts_over+2
        call dgemm_guru('t','n',npols,naintbc,npolso,done,rfds(ixist), &
          npolso,wquadf,npolso,dzero,wquadf2,npols)
        
        do i=1,naintbc
          wquad(iaintbc(i),:) = wquadf2(:,i)
        enddo
       
        do i = 1,npols
          ipt = ixyzs(ipatch) + i-1
          do j=col_ptr(ipt),col_ptr(ipt+1)-1
             jind = row_ind(j)

             j2 = j-col_ptr(ixyzs(ipatch))+1
             juniind = iuniind(j2)
             dmatent(j) = wquad(juniind,i)
          enddo
        enddo

        deallocate(iuni,iuniind,aintb,iaintba,iaintbb,aintbc,iaintbc)
        deallocate(wquad,wquadf,wquadf2)
      enddo
      
      return
      end
