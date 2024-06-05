!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!   Combined field representation for Stokes velocity problem
!
!  PDE:
!    \Delta u       = \nabla p
!    \Delta \cdot u = 0
!
!  Boundary conditions:
!    u = f
!
!  Representation:
!    u = (\alpha S_{stok}  + \beta D_{stok})[\sigma]
!
!  Integral equations obtained by imposing:
!    u = f
!
!  and is given by:
!    d \sigma + (\alpha S_{stok}  + \beta D_{stok})[\sigma] = f
!
!  where:
!    d = -1/2 beta for interior problem
!      =  1/2 beta for exterior problem
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  User callable routines:
!    - stok_comb_vel_solver: Given data f,and parameters \alpha, \beta, 
!        this routine returns the solution \sigma
!
!    - stok_comb_vel_eval: Given \sigma, and parameters \alpha, \beta, 
!        evaluates the solution at a collection of targets
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  Advanced interfaces:
!    - getnearquad_stok_comb_vel: compute the near quadrature correction
!        for constructing the on-surface integral equation for the
!        velocity data corresponding to the combined field
!        representation with user-provided near-information prescribed
!        in row-sparse compressed format
!
!    - getnearquad_stok_comb_vel_eval: compute the near quadrature
!        correction for target points which can be either on-surface or
!        off-surface with user-provided near-information prescribed in
!        row-sparse compressed format
!
!    - lpcomp_stok_comb_vel_addsub: apply the principal value part
!        of the integeral equation on surface. On input, user provides
!        precomputed near quadrature in row-sparse compressed format,
!        and oversampling information for handling the far part of the
!        computation
!
!    - stok_comb_vel_eval_addsub: compute the solution u at a
!        collection of targets (on-surface or off-surface), given \sigma.
!        On input, user provides precomputed near quadrature in
!        row-sparse compressed format and oversampling surface
!        information for the far-part
!
!    - stok_comb_vel_solver_guru: Guru solver routine, where user is
!        responsible for providing precomputed near quadrature
!        information in row-sparse compressed format and oversampling
!        surface information for the far-part
!
!  Other routines in the file, these are currently in beta mode
!  and should be used at your own risk:
!
!    - stok_comb_vel_matgen: matrix entry generator for Stokes
!      layer potential operators
!

      subroutine getnearquad_stok_comb_vel(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, eps, dpars, &
        iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear)
!
!  This subroutine generates the near field quadrature
!  for the representation:
!
!  u = \alpha S_{stok}[\sigma] + \beta D_{stok}[\sigma]]    -   (1)
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
!  u = \pm \beta \sigma/2 + \alpha S_{stok}[\sigma] + \beta D_{stok}[\sigma] 
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
!    - eps: real *8
!        precision requested
!    - dpars: real *8(2)
!        kernel parameters (Referring to formula (1))
!        dpars(1) = alpha (single layer strength) 
!        dpars(2) = beta (double layer strength)
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
!    - rfac0: real *8
!        radius parameter for switching to predetermined quadarature
!        rule        
!    - nquad: integer
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
      integer, intent(in) :: npatches, norders(npatches), npts, nquad
      integer, intent(in) :: ixyzs(npatches+1), iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts), eps
      real *8, intent(in) :: rfac0
      integer, intent(in) :: iquadtype
      real *8, intent(in) :: dpars(2)
      integer, intent(in) :: nnz
      integer, intent(in) :: row_ptr(npts+1), col_ind(nnz), iquad(nnz+1)
      real *8, intent(out) :: wnear(6,nquad)

      integer :: ndtarg, ntarg
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      integer i

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

      call getnearquad_stok_comb_vel_eval(npatches, norders, &
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
      subroutine getnearquad_stok_comb_vel_eval(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, &
        targs, ipatch_id, uvs_targ, eps, dpars, iquadtype, nnz, &
        row_ptr, col_ind, iquad, rfac0, nquad, wnear)
!
!  This subroutine generates the near field quadrature
!  for the representation:
!
!  u = \alpha S_{stok}[\sigma] + \beta D_{stok}[\sigma]]    -   (1)
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
!    - dpars: real *8(2)
!        kernel parameters (Referring to formula (1))
!        dpars(1) = alpha (single layer strength) 
!        dpars(2) = beta (double layer strength)
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
      integer, intent(in) :: npatches, norders(npatches), npts, nquad
      integer, intent(in) :: ixyzs(npatches+1), iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts), eps
      real *8, intent(in) :: rfac0
      integer, intent(in) :: iquadtype
      real *8, intent(in) :: dpars(2)
      integer, intent(in) :: nnz
      integer, intent(in) :: ndtarg, ntarg
      real *8, intent(in) :: targs(ndtarg, ntarg)
      integer, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      integer, intent(in) :: row_ptr(ntarg+1), col_ind(nnz), iquad(nnz+1)
      real *8, intent(out) :: wnear(6,nquad)

      real *8 alpha, beta
      real *8, allocatable :: wneartmp(:)
      integer ndd, ndi, ndz
      complex *16 zpars(1)
      integer ipars(2)

      integer i, ipv, j, ii, ijloc(2,6)

      procedure (), pointer :: fker
      external st3d_comb 
      
!
!        initialize the appropriate kernel function
!

      alpha = dpars(1)
      beta = dpars(2)


      ndd = 2
      ndi = 2
      ndz = 0
      fker => st3d_comb
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
      subroutine stok_comb_vel_eval(npatches, norders, ixyzs, &
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
!  Stokes combined field representation 
!
!  Representation:
!    u = (\alpha \mathcal{S}_{stok} + \beta \mathcal{D}_{stok})[\sigma]
!
!  Note: For targets on the boundary, this routine only computes
!  the principal value part, the identity term corresponding to the jump
!  in the layer potential is not included in the layer potential.
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
!    - dpars: real *8(2)
!        kernel parameters (Referring to formula above)
!        dpars(1) = $\alpha$
!        dpars(2) = $\beta$
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
      integer, intent(in) :: npatches, npts
      integer, intent(in) :: ndtarg, ntarg
      integer, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts), eps
      real *8, intent(in) :: targs(ndtarg,ntarg)
      real *8, intent(in) :: dpars(2)
      real *8, intent(in) :: sigma(3,npts)
      integer, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)

      real *8, intent(out) :: pot(3,ntarg)


      integer nptso, nnz, nquad

      integer nover, npolso
      integer norder ,npols
      integer, allocatable :: row_ptr(:), col_ind(:), iquad(:)
      real *8, allocatable :: wnear(:,:)

      real *8, allocatable :: srcover(:,:), wover(:)
      integer, allocatable :: ixyzso(:), novers(:)

      real *8, allocatable :: cms(:,:), rads(:), rad_near(:) 

      integer i, j, jpatch, jquadstart, jstart

      complex *16 zpars
      integer ipars

      real *8 ttot, done, pi
      real *8 rfac, rfac0
      integer iptype_avg, norder_avg
      integer ikerorder, iquadtype, npts_over

      integer ndd, ndz, ndi, nker, lwork, ndim_s, ndim_p
      integer idensflag, ipotflag
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
      allocate(wnear(6,nquad))
      
!$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,nquad
        wnear(1,i) = 0
        wnear(2,i) = 0
        wnear(3,i) = 0
        wnear(4,i) = 0
        wnear(5,i) = 0
        wnear(6,i) = 0
      enddo
!$OMP END PARALLEL DO    


      iquadtype = 1

      call getnearquad_stok_comb_vel_eval(npatches, norders, &
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
      nker = 6
      lwork = 0
      ndim_s = 3
      idensflag = 0
      ipotflag = 0
      ndim_p = 3


!
      call stok_comb_vel_eval_addsub(npatches, norders, ixyzs, &
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
      subroutine lpcomp_stok_comb_vel_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, &
        novers, nptso, ixyzso, srcover, whtsover, lwork, work, ndim, &
        sigma, pot)
!
!  This subroutine evaluates the velocity data corresponding to
!  the following integral representation:
!  
!  u = \alpha S_{stok}[\sigma] + \beta*D_{stok}[\sigma]
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
!        integral representation (must be 3)
!    - dpars: real *8(ndd)
!        real parameters defining the kernel/
!        integral representation 
!        * dpars(1) = alpha (single layer strength) 
!        * dpars(2) = beta (double layer strength)
!        * dpars(3) = strength of 1s matrix of 
!            n \int \sigma \cdot n
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
!        number of kernels in quadrature correction, must be 6
!    - wnear: real *8(nker, nquad)
!        Precomputed near field quadrature for
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
!        size of work array (must be npts)
!    - work: real *8(lwork)
!        work array
!        * work(1:npts) stores the smooth quadrature weights for
!          integration on the discretized surface
!    - ndim: integer
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
      integer, intent(in) :: npatches, npts
      integer, intent(in) :: norders(npatches),ixyzs(npatches+1)
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
    
      real *8, intent(in) :: sigma(ndim,npts)
    
      real *8, intent(out) :: pot(ndim,npts)

      integer ndtarg, idensflag, ipotflag, ndim_p, i
      real *8 rint

      ndtarg = 12
      idensflag = 0
      ipotflag = 0
      ndim_p = 3

      call stok_comb_vel_eval_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, npts, srcvals, &
        eps, ndd, dpars, ndz, zpars, ndi, ipars, &
        nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
        nptso, ixyzso, srcover, whtsover, lwork, work, idensflag, &
        ndim, sigma, ipotflag, ndim_p, pot)
!
!  Add in correction due to 1's matrix
!
      rint = 0
!$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rint)
      do i = 1,npts
        rint = rint + (sigma(1,i)*srcvals(10,i) + &
                       sigma(2,i)*srcvals(11,i) + &
                       sigma(3,i)*srcvals(12,i))*work(i)
      enddo
!$OMP END PARALLEL DO      
      rint = rint*dpars(3)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)      
      do i = 1,npts
        pot(1,i) = pot(1,i) + rint*srcvals(10,i)
        pot(2,i) = pot(2,i) + rint*srcvals(11,i)
        pot(3,i) = pot(3,i) + rint*srcvals(12,i)
      enddo
!$OMP END PARALLEL DO

      return
      end
!
!
!
!
      subroutine stok_comb_vel_eval_addsub(npatches, norders, ixyzs, &
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
!  u = \alpha S_{stok}[\sigma] + \beta*D_{stok}[\sigma] 
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
!        integral representation (must be 3)
!    - dpars: real *8(ndd)
!        real parameters defining the kernel/
!        integral representation 
!        * dpars(1) = alpha (single layer strength) 
!        * dpars(2) = beta (double layer strength)
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
!        number of kernels in quadrature correction, must be 6
!    - wnear: real *8(nker, nquad)
!        Precomputed near field quadrature for
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
!        must be 3 for this routine
!    - sigma: real *8(ndim_s, npts)
!        The density sigma
!    - ipotflag: integer
!        Flag for determining output type 
!        ipotflag = 1, only velocity is returned
!    - ndim_p: integer
!        number of potentials per point
!        ndim_p = 3, if ipotflag = 1
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
      
      real *8, intent(in) :: sigma(ndim_s,npts)
      
      real *8, intent(out) :: pot(ndim_p,ntarg)

      integer norder,npols,nover,npolso
      real *8, allocatable :: potsort(:)

      real *8, allocatable :: sources(:,:), targvals(:,:)
      real *8, allocatable :: sigmaover(:,:)
      real *8, allocatable :: stoklet(:,:), strslet(:,:)
      real *8, allocatable :: pottarg(:,:), strsvec(:,:)
      integer ns, nt
      real *8 alpha, beta
      integer ifstoklet, ifstrslet
      integer ifppreg, ifppregtarg
      real *8 tmp(10), val
      real *8 over4pi

      real *8 w11, w12, w13, w21, w22, w23, w31, w32, w33
      real *8 sig1, sig2, sig3

      integer istress

      integer i, j, jpatch, jquadstart, jstart

      real *8 pottmp
      real *8, allocatable :: sttmp2(:,:), strstmp2(:,:), strsvec2(:,:)

      integer nmax
      real *8 timeinfo(10), t1, t2, omp_get_wtime

      real *8, allocatable :: srctmp2(:,:)
      real *8 thresh, ra
      real *8 rr, rmin
      integer nss, ii, l, npover

      integer nd, ntarg0
      integer ier, iper

      integer ndtmp, ndsigma

      real *8 ttot, done, pi

      real *8 potv(3), pv, gradv(3,3)

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
      call get_near_corr_max(ntarg, row_ptr, nnz, col_ind, npatches, &
        ixyzso, nmax)
      allocate(srctmp2(3,nmax), sttmp2(3,nmax), strstmp2(3,nmax))
      allocate(strsvec2(3,nmax))
           
           
      ifppreg = 0
      ifppregtarg = 1
      allocate(sources(3,ns), targvals(3,ntarg))
      allocate(stoklet(3,ns), strslet(3,ns), strsvec(3,ns))
      allocate(sigmaover(3,ns))
      allocate(pottarg(3,ntarg))

! 
!       oversample density
!

      ndsigma = 3
      call oversample_fun_surf(ndsigma, npatches, norders, ixyzs, &
        iptype, npts, sigma, novers, ixyzso, ns, sigmaover)


!
!      set relevant parameters for the fmm
!
 
      alpha = dpars(1)*over4pi
      beta = dpars(2)*over4pi


!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)      
      do i=1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)

        stoklet(1,i) = sigmaover(1,i)*whtsover(i)*alpha
        stoklet(2,i) = sigmaover(2,i)*whtsover(i)*alpha
        stoklet(3,i) = sigmaover(3,i)*whtsover(i)*alpha
!
!  Stokes double layer is negative of a stresslet with orienatation
!  given by surface normal
!
        strslet(1,i) = -sigmaover(1,i)*whtsover(i)*beta 
        strslet(2,i) = -sigmaover(2,i)*whtsover(i)*beta 
        strslet(3,i) = -sigmaover(3,i)*whtsover(i)*beta

        strsvec(1,i) = srcover(10,i)
        strsvec(2,i) = srcover(11,i)
        strsvec(3,i) = srcover(12,i)

      enddo
!$OMP END PARALLEL DO      
      
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i = 1,ntarg
        targvals(1,i) = targs(1,i)
        targvals(2,i) = targs(2,i)
        targvals(3,i) = targs(3,i)
      enddo
!$OMP END PARALLEL DO      
      

      ifstoklet = 1
      ifstrslet = 1

      if(alpha.eq.0) ifstoklet = 0
      if(beta.eq.0) ifstrslet = 0

      iper = 0
      ier = 0

!
!
!       call the fmm
!

      call cpu_time(t1)
!$      t1 = omp_get_wtime()      
      call stfmm3d(nd, eps, ns, sources, ifstoklet, stoklet, &
        ifstrslet, strslet, strsvec, ifppreg, tmp, tmp, tmp, ntarg, &
        targvals, ifppregtarg, pottarg, tmp, tmp, ier)
      call cpu_time(t2)
!$      t2 = omp_get_wtime()

      timeinfo(1) = t2-t1
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i = 1,ntarg
        pot(1,i) = pottarg(1,i)
        pot(2,i) = pottarg(2,i)
        pot(3,i) = pottarg(3,i)
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
            
            strstmp2(1,nss) = strslet(1,l)
            strstmp2(2,nss) = strslet(2,l)
            strstmp2(3,nss) = strslet(3,l)
            
            strsvec2(1,nss) = strsvec(1,l)
            strsvec2(2,nss) = strsvec(2,l)
            strsvec2(3,nss) = strsvec(3,l)
          enddo
        enddo

        potv(1) = 0
        potv(2) = 0
        potv(3) = 0
        pv = 0
        gradv(1:3,1:3) = 0
        
        istress = 1
        call st3ddirectstokstrsg(nd, srctmp2, sttmp2, istress, &
          strstmp2, strsvec2, nss, targvals(1,i), ntarg0, potv, pv, &
          gradv, thresh)
        
        pot(1,i) = pot(1,i) - potv(1)
        pot(2,i) = pot(2,i) - potv(2)
        pot(3,i) = pot(3,i) - potv(3)

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
      subroutine stok_comb_vel_solver(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, dpars, numit, ifinout, &
        rhs, eps_gmres, niter, errs, rres, soln)

!
!
!  This subroutine solves the Stokes velocity problem
!  on the exterior/interior of an object where the potential
!  is represented using a combined field integral representation.
!
!  This subroutine is the simple interface as opposed to the
!  _solver_guru routine which is called after initialization
!  in this routine.
!
!
!  Representation:
!    u = \alpha S_{stok}[\sigma]+\beta*D_{stok}[\sigma]
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
!    - dpars: real *8(2)
!        kernel parameters (Referring to formula (1))
!          * dpars(1) = \alpha 
!          * dpars(2) = \beta
!    - numit: integer
!        max number of gmres iterations
!    - ifinout: integer
!        ifinout = 0, interior problem
!        ifinout = 1, exterior problem
!    - rhs: real *8(3,npts)
!        velocity data
!    - eps_gmres: real *8
!        gmres tolerance requested
!      
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
!    - soln: real *8(3,npts)
!        density which solves the velocity problem \sigma
!				 
      implicit none
      integer, intent(in) :: npatches, npts
      integer, intent(in) :: ifinout
      integer, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps, eps_gmres
      real *8, intent(in) :: dpars(2)
      real *8, intent(in) :: rhs(3,npts)
      integer, intent(in) :: numit
      real *8, intent(out) :: soln(3,npts)
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
      real *8, allocatable :: wnear(:,:)

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
      allocate(wnear(nker,nquad))
      
!$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,nquad
        wnear(1,i) = 0
        wnear(2,i) = 0
        wnear(3,i) = 0
        wnear(4,i) = 0
        wnear(5,i) = 0
        wnear(6,i) = 0
      enddo
!$OMP END PARALLEL DO    


      iquadtype = 1

      print *, "starting to generate near quadrature"
      call cpu_time(t1)
!$      t1 = omp_get_wtime()      

      call getnearquad_stok_comb_vel(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, eps, dpars, iquadtype, &
        nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear)
      call cpu_time(t2)
!$      t2 = omp_get_wtime()     

      print *, "done generating near quadrature, now starting gmres"
      
      call stok_comb_vel_solver_guru(npatches, norders, ixyzs, &
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

      subroutine stok_comb_vel_solver_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, dpars, numit, ifinout, &
        rhs, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
        nptso, ixyzso, srcover, whtsover, eps_gmres, niter, &
        errs, rres, soln)
!
!
!  This subroutine solves the Stokes velocity problem
!  on the exterior of an object where the potential
!  is represented using the combined field integral representation.
!
!
!  Representation:
!    u = \alpha S_{stok}[\sigma] + \beta*D_{stok}[\sigma]
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
!    - dpars: real *8(2)
!        kernel parameters (Referring to formula (1))
!          * dpars(1) = \alpha 
!          * dpars(2) = \beta 
!    - numit: integer
!        max number of gmres iterations
!    - ifinout: integer
!        ifinout = 0, interior problem
!        ifinout = 1, exterior problem
!    - rhs: real *8(3,npts)
!        velocity data
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
!        number of kernels in quadrature correction, must be 6
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
!    - soln: real *8(3,npts)
!        density which solves the neumann problem \sigma
!				 

      implicit none
      integer, intent(in) :: npatches, npts
      integer, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps, eps_gmres
      integer, intent(in) :: ifinout
      real *8, intent(in) :: dpars(2)
      real *8, intent(in) :: rhs(3,npts)
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
      
  
      real *8, intent(out) :: soln(3,npts)

      real *8 did
      real *8 dpars_use(3), pi, done, rsurf
      integer ndd_use

      procedure (), pointer :: fker
      external lpcomp_stok_comb_vel_addsub

      integer ndd, ndi, ndz, lwork, ndim
      real *8 work
      integer ipars, nkertmp
      complex *16 zpars

      integer ndtarg
      complex *16 ima
      data ima/(0.0d0,1.0d0)/
      real *8, allocatable :: wts(:)

      complex *16 zpars_use(3)
      integer i

      did = -(-1)**(ifinout)*dpars(2)/2
      fker => lpcomp_stok_comb_vel_addsub

      ndd_use = 3
      dpars_use(1) = dpars(1)
      dpars_use(2) = dpars(2)
      dpars_use(3) = 0
      
      ndi = 0
      ndz = 0

      lwork = 0
      ndim = 3
      allocate(wts(npts))

      call get_qwts(npatches, norders, ixyzs, iptype, npts, &
        srcvals, wts)
      rsurf = 0
!$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rsurf)      
      do i=1,npts
        rsurf = rsurf + wts(i)
      enddo
!$OMP END PARALLEL DO      
      done = 1.0d0
      pi = atan(done)*4.0d0
      if(ifinout.eq.0) dpars_use(3) = -2*pi/rsurf
! 

      call dgmres_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, wts, &
        eps, ndd_use, dpars_use, ndz, zpars, ndi, ipars, &
        nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
        nptso, ixyzso, srcover, whtsover, npts, wts, &
        ndim, fker, did, rhs, numit, eps_gmres, niter, errs, &
        rres, soln)

      return
      end

!
!
!
!        
      subroutine stok_comb_vel_matgen(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, dpars, ifinout, &
        xmat)
!
!f2py  intent(in) npatches,norders,ixyzs,iptype
!f2py  intent(in) npts,srccoefs,srcvals,eps,dpars
!f2py  intent(in) ifinout
!f2py  intent(out) xmat
!
!     this subroutine returns the discretization matrix
!     for the on-surface discretization of the stokes combined
!     field layer potential. The unknowns are ordered as:
!
!     \sigma_{1}(x_{1}), \sigma_{2}(x_{1}), \sigma_{3}(x_{1})...
!      \sigma_{1}(x_{2}), \sigma_{2}(x_{2})l \sigma_{3}(x_{2})..
!                           .
!                           .
!                           .
!                           .
!      \sigma_{1}(x_{n}), \sigma_{2}(x_{n}), \sigma_{3}(x_{n})
!
!     And the same ordering for the velocity components as well
!
!
!     Representation:
!        u = alpha S \sigma + beta D \sigma
!     
!
!       input:
!         npatches - integer
!            number of patches
!
!         norders- integer(npatches)
!            order of discretization on each patch 
!
!         ixyzs - integer(npatches+1)
!            ixyzs(i) denotes the starting location in srccoefs,
!               and srcvals array corresponding to patch i
!   
!         iptype - integer(npatches)
!            type of patch
!             iptype = 1, triangular patch discretized using RV nodes
!
!         npts - integer
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
!          eps - real *8
!             precision requested for computing quadrature and fmm
!             tolerance
!
!          dpars - real *8 (2)
!             alpha = dpars(1), beta = dpars(2) in the layer potential
!             representation
!      
!          ifinout - integer
!              flag for interior or exterior problems (normals assumed to 
!                be pointing in exterior of region)
!              ifinout = 0, interior problem
!              ifinout = 1, exterior problem
!
!         output
!           xmat - real *8(3*npts,3*npts)
!              discretization matrix for the stokes boundary value
!              problem
!
!
      implicit none
      integer npatches,norder,npols,npts
      integer ifinout
      integer norders(npatches),ixyzs(npatches+1)
      integer iptype(npatches)
      real *8 srccoefs(9,npts), srcvals(12,npts), eps, eps_gmres
      real *8 dpars(2)
      real *8 xmat(3*npts,3*npts)

      real *8 uint

      real *8, allocatable :: targs(:,:)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      integer ndtarg,ntarg

      real *8 rres,eps2
      integer niter


      integer nover,npolso,nptso
      integer nnz,nquad
      integer, allocatable :: row_ptr(:), col_ind(:), iquad(:)
      real *8, allocatable :: wnear(:,:), wts(:)

      real *8, allocatable :: srcover(:,:), wover(:)
      integer, allocatable :: ixyzso(:), novers(:)

      real *8, allocatable :: cms(:,:), rads(:), rad_near(:) 

      integer i, j, jpatch, jquadstart, jstart

      integer ipars
      complex *16 zpars
      real *8 timeinfo(10), t1, t2, omp_get_wtime


      real *8 ttot, done, pi, rsurf
      real *8 rfac, rfac0, alpha, beta
      integer iptype_avg, norder_avg
      integer ikerorder, iquadtype, npts_over

      real *8 did, ra
      integer jj, l, nmat
      real *8 w11, w12, w13, w21, w22, w23, w31, w32, w33
      

      alpha = dpars(1)
      beta = dpars(2)
      
      nmat = 3*npts
      done = 1
      pi = atan(done)*4


!
!
!        this might need fixing
!
      iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
      norder_avg = floor(sum(norders)/(npatches+0.0d0))

      call get_rfacs(norder_avg, iptype_avg, rfac, rfac0)

      nnz = npts*npatches
      allocate(row_ptr(npts+1), col_ind(nnz))

      do i=1,npts + 1
        row_ptr(i) = (i-1)*npatches + 1
      enddo


      do i=1,npts
        do j=1,npatches
          col_ind((i-1)*npatches + j) = j
        enddo
      enddo

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches, ixyzs, npts, nnz, row_ptr, col_ind, &
        iquad)

      ikerorder = -1
      if(abs(dpars(2)).gt.1.0d-16) ikerorder = 0
      
!
!   compute near quadrature correction
!
      nquad = iquad(nnz+1)-1
      allocate(wnear(6,nquad))
      
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j)
      do i=1,nquad
        do j = 1,6
          wnear(j,i) = 0
        enddo
      enddo
!$OMP END PARALLEL DO    


      iquadtype = 1

      print *, "starting to generate near quadrature"
      call cpu_time(t1)
!$      t1 = omp_get_wtime()      

      call getnearquad_stok_comb_vel(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, &
        eps, dpars, iquadtype, nnz, row_ptr, col_ind, &
        iquad, rfac0, nquad, wnear)
      call cpu_time(t2)
!$      t2 = omp_get_wtime()     

      call prin2('quadrature generation time=*',t2-t1,1)

      allocate(wts(npts))
      call get_qwts(npatches, norders, ixyzs, iptype, npts, &
        srcvals, wts)


!
!  compute scaling for one's matrix correction for interior problem
!

      rsurf = 0
      do i=1,npts
        rsurf = rsurf + wts(i)
      enddo
      
      ra = 0
      if(ifinout.eq.0) ra = -2*pi/rsurf
      

      call cpu_time(t1)
!$      t1 = omp_get_wtime()

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, jpatch, jquadstart) &
!$OMP PRIVATE(jstart, npols, l, jj, w11, w12, w13) &
!$OMP PRIVATE(w21, w22, w23, w31, w32, w33)
      do i = 1,npts
        do j = row_ptr(i), row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
             jj = jstart + l-1
             w11 = wnear(1,jquadstart+l-1)
             w12 = wnear(2,jquadstart+l-1)
             w13 = wnear(3,jquadstart+l-1)
             w21 = w12
             w22 = wnear(4,jquadstart+l-1)
             w23 = wnear(5,jquadstart+l-1)
             w31 = w13
             w32 = w23
             w33 = wnear(6,jquadstart+l-1)
             xmat(3*(i-1)+1,3*(jj-1)+1) = w11 + & 
               srcvals(10,i)*srcvals(10,jj)*ra*wts(jj)
             xmat(3*(i-1)+1,3*(jj-1)+2) = w12 + &
                srcvals(10,i)*srcvals(11,jj)*ra*wts(jj)
             xmat(3*(i-1)+1,3*(jj-1)+3) = w13 + &
                srcvals(10,i)*srcvals(12,jj)*ra*wts(jj)

             xmat(3*(i-1)+2,3*(jj-1)+1) = w21 + &
               srcvals(11,i)*srcvals(10,jj)*ra*wts(jj)
             xmat(3*(i-1)+2,3*(jj-1)+2) = w22 + &
               srcvals(11,i)*srcvals(11,jj)*ra*wts(jj)
             xmat(3*(i-1)+2,3*(jj-1)+3) = w23 + &
               srcvals(11,i)*srcvals(12,jj)*ra*wts(jj)

             xmat(3*(i-1)+3,3*(jj-1)+1) = w31 + &
               srcvals(12,i)*srcvals(10,jj)*ra*wts(jj)
             xmat(3*(i-1)+3,3*(jj-1)+2) = w32 + & 
               srcvals(12,i)*srcvals(11,jj)*ra*wts(jj)
             xmat(3*(i-1)+3,3*(jj-1)+3) = w33 + &
               srcvals(12,i)*srcvals(12,jj)*ra*wts(jj)

          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
     

      did = -(-1)**(ifinout)*beta/2
      do i=1,3*npts
         xmat(i,i) = xmat(i,i) + did
      enddo

     
!     
      return
      end
!     
!
!
