!
!  Non resonant charge current integral equation (NRCCIE)
!  for scattering from perfect electric conductors
!
!  PDE:
!    \nabla \times E =  ik H
!    \nabla \cdot  E =     0
!    \nabla \times H = -ik E
!    \nabla \cdot  H =     0
!
!  Boundary conditions
!    n \times (E + E_in) = 0
!    n \cdot  (E + E_in) = \rho
!
!    n \times (H + H_in) = J
!    n \cdot  (H + H_in) = 0
!  
!
!  Representation:
!    H = \nabla \times S_{k} [J]
!    E = ik S_{k} [J] - \nabla S_{k} [\rho]
!
!  Integral equations solve for [J,\rho] and obtained by imposing
!
!  n \times (H + H_in) - \alpha (n \times n \times (E + E_in)) = J
!  n \cdot  (E + E_in) + \alpha/ik \nabla \cdot E              = \rho
!
!  and are given by
!  J/2 - n \times \nabla S_{k} [J] +
!     \alpha (n \times n \times (ik S_{k} [J] - \nabla S_{k} [\rho])) = 
!     n \times H_in - \alpha n \times n \times E_in      ----  (1)
!
!  \rho/2 + S_{k}'[\rho] - ik n \cdot S_{k} [J] + 
!     \alpha (\nabla \cdot S_{k}[J] -ik S_{k}[\rho]) = n \cdot E_in 
!
!
!  Like all other maxwell formulations, the vector densties and equation
!  are expanded in a local orthonormal basis given by
!  \partial_{u} xyz, and \bn \times \partial_{u} xyz appropriately 
!  normalized. Let $\ru$ and $\rv$ denote those unit tangent vectors, 
!  then the currents are given by
!
!  so J = j_{u} \ru + j_{v} \rv  
!
!  and equation (1) is imposed by dotting it with \ru and \rv
!
!  User callable routines:
!    - em_nrccie_pec_solver: Given incident electric and magnetic 
!        fields E_in, H_in, helmholtz wave number k,
!        and parameter \alpha, this routine returns the current J, and
!        charge \rho
!
!    - em_nrccie_pec_eval: Given J, and \rho,
!        evaluate the electric and magnetic fields at a collection of targets
!
!  Advanced interfaces:
!    - getnearquad_em_nrccie_pec: Computes the quadrature
!        correction for constructing the on-surface integral equation
!        with user-provided near-information prescribed in 
!        row-sparse compressed format
!
!    - lpcomp_em_nrccie_pec_addsub: Apply the principal value part
!        of the integral equation on surface. On input, user
!        provides precomputed near quadrature in row-sparse
!        compressed format, and oversampling information for 
!        handling the far part of the computation
!
!    - em_nrccie_pec_solver_guru: Guru solver routine, 
!        where user is responsible for providing precomputed
!        near quadrature information in row-sparse compressed format
!        and oversampling surface information for the far-part
!
!    - getnearquad_em_nrccie_pec_eval: Generate quadrature
!        for the post-processor, where targets can be in the volume
!        or on surface
!
!
!    - em_nrccie_pec_eval_addsub: Compute the solution
!        E, H at a collection of targets, given J and
!        \rho. On input, user provides precomputed
!        near quadrature in row-sparse compressed format, 
!        and oversampling surface information for handling the
!        far part
!   



      subroutine getnearquad_em_nrccie_pec(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, eps, zpars, &
        iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
        wnear)
!
!  This subroutine generates the near field quadrature
!  for the non-resonant charge current integral equation.        
!
!  The quadrature is computed by the following strategy
!  targets within a sphere of radius rfac0*rs
!  of a chunk centroid is handled using adaptive integration
!  where rs is the radius of the bounding sphere
!  for the patch
!  
!  All other targets in the near field are handled via
!  oversampled quadrature
!
!  The recommended parameter for rfac0 is 1.25d0
!
!  NOTES:
!  This subroutine returns 9 kernels
!  1) \ru \cdot (-n \times \nabla \times S_{k}  + \alpha n \times n \times S_{k})[j_{u} \ru]
!  2) \ru \cdot (-n \times \nabla \times S_{k}  + \alpha n \times n \times S_{k})[j_{v} \rv]
!  3) \ru \cdot (-\alpha n \times n \times \nabla S_{k}[\rho]) 
!  4) \rv \cdot (-n \times \nabla \times S_{k}  + \alpha n \times n \times S_{k})[j_{u} \ru]
!  5) \rv \cdot (-n \times \nabla \times S_{k}  + \alpha n \times n \times S_{k})[j_{v} \rv]
!  6) \rv \cdot (-\alpha n \times n \times \nabla S_{k}[\rho])
!  7) (-ik n \cdot S_{k} + \alpha \nabla \cdot S_{k})[j_{u} \ru]
!  8) (-ik n \cdot S_{k} + \alpha \nabla \cdot S_{k})[j_{v} \rv]
!  9) S_{k}'[\rho]
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
!    - zpars: complex *16 (2)
!        kernel parameters
!          * zpars(1) = k 
!          * zpars(2) = alpha
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
!        pair. The size of wnear is 9*nquad since there are 9 kernels
!        per source target pair
!
!    output
!      wnear - complex *16 (9, nquad) 
!        Near field quadrature corrections for nrccie
!        * wnear(1,:) stores corrections for
!          \ru \cdot (-n \times \nabla \times S_{k}  + 
!          \alpha n \times n \times S_{k})[j_{u} \ru]
!        * wnear(2,:) stores corrections for
!          \ru \cdot (-n \times \nabla \times S_{k}  + 
!          \alpha n \times n \times S_{k})[j_{v} \rv]
!        * wnear(3,:) stores corrections for
!            \ru \cdot (-\alpha n \times n \times \nabla S_{k}[\rho])        
!        * wnear(4,:) stores corrections for
!          \rv \cdot (-n \times \nabla \times S_{k}  + 
!          \alpha n \times n \times S_{k})[j_{u} \ru]
!        * wnear(5,:) stores corrections for
!          \rv \cdot (-n \times \nabla \times S_{k}  + 
!          \alpha n \times n \times S_{k})[j_{v} \rv]
!        * wnear(6,:) stores corrections for
!          \rv \cdot (-\alpha n \times n \times \nabla S_{k}[\rho])
!        * wnear(7,:) stores corrections for
!         (-ik n \cdot S_{k} + \alpha \nabla \cdot S_{k})[j_{u} \ru]
!        * wnear(8,:) stores corrections for
!         (-ik n \cdot S_{k} + \alpha \nabla \cdot S_{k})[j_{v} \rv]
!        * wnear(9,:) stores corrections for S_{k}'[\rho]        
!

      implicit none 
      integer, intent(in) :: npatches, norders(npatches), npts, nquad
      integer, intent(in) :: ixyzs(npatches+1), iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts), eps
      real *8, intent(in) :: rfac0
      
      integer, intent(in) :: iquadtype, nnz
      integer, intent(in) :: row_ptr(npts+1), col_ind(nnz)
      integer, intent(in) :: iquad(nnz+1)
      
      complex *16, intent(in) :: zpars(2)
      complex *16, intent(out) :: wnear(9, nquad)
      
      integer ndtarg
      
      integer ipars(2)
      real *8 dpars(1)
      
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_src(:,:)
      complex *16, allocatable :: wneartmp(:)


      complex *16 alpha, beta
      integer i, j, ndi, ndd, ndz, iuse, idim, l

      integer ipv

      procedure (), pointer :: fker
      external  fker_em_nrccie_pec_s, h3d_sprime

      ndz=2
      ndd=1
      ndi=2
      ipv=1
      ndtarg = 12

      fker =>  fker_em_nrccie_pec_s

      allocate(ipatch_id(npts), uvs_src(2,npts))
      !$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        ipatch_id(i) = -1
        uvs_src(1,i) = 0
        uvs_src(2,i) = 0
      enddo
!$OMP END PARALLEL DO      

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, idim)
      do i=1,nquad
        do idim = 1,9
          wnear(idim,i) = 0
        enddo
      enddo
!$OMP END PARALLEL DO

      allocate(wneartmp(nquad))
      call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts, &
        ipatch_id, uvs_src)

      ipv = 1
      do j=1,3
        do i=1,3
          ipars(1) = j
          ipars(2) = i
          iuse = (j-1)*3 + i
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(l)
          do l=1,nquad
            wneartmp(l) = 0
          enddo
!$OMP END PARALLEL DO                   
          call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
            iptype, npts, srccoefs, srcvals, ndtarg, npts, srcvals, &
            ipatch_id, uvs_src, eps, ipv, fker, ndd, dpars, ndz, &
            zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, &
            nquad, wneartmp)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(l)
          do l=1,nquad
            wnear(iuse,l) = wneartmp(l)
          enddo
!$OMP END PARALLEL DO  
        enddo
      enddo
      
      return
      end subroutine getnearquad_em_nrccie_pec
!
!
!
!
!

      subroutine lpcomp_em_nrccie_pec_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, &
        novers, nptso, ixyzso, srcover, whtsover, lwork, work, &
        ndim, sigma, pot)
!
!
!  This subroutine evaluates the layer potential for
!  the boundary integral equation:
!
!    J/2 - M_{k}[J] + alpha·nxnx(ikS_{k}[J] - grad S_{k}[rho]) = 
!      = nxH_inc + alpha nxnxE_inc
!    rho/2 + S'_{k}[rho] - ikS_{k}[J] + 
!          alpha(divS_{k}[J] - ikS_{k}[rho]) = n·E_inc
!
!  where the near field is precomputed and stored
!  in the row sparse compressed format.
!
!  Note the 4\pi scaling is NOT included as the FMM output was scaled
!  appropriately
!
!  Note: the identity J/2 is not included as the gmres takes care of
!  that
!
!  The fmm is used to accelerate the far-field and 
!  near-field interactions are handled via precomputed quadrature
!
!  Using add and subtract - no need to call tree and set fmm 
!  parameters can directly call existing fmm library
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
!        integral representation (unused in this routine)
!    - dpars: real *8(ndd)
!        real parameters defining the kernel/
!        integral represnetation (unused in this routine)
!    - ndz: integer
!        number of complex parameters defining the kernel/
!        integral representation, must be >=2 
!    - zpars: complex *16(ndz)
!        complex parameters defining the kernel/
!        integral represnetation.
!        * zpars(1) = k 
!        * zpars(2) = alpha
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
!        number of kernels in quadrature correction, must be 9
!    - wnear: complex *16(nker, nquad)
!        Near field quadrature corrections for nrccie
!        * wnear(1,:) stores corrections for
!          \ru \cdot (-n \times \nabla \times S_{k}  + 
!          \alpha n \times n \times S_{k})[j_{u} \ru]
!        * wnear(2,:) stores corrections for
!          \ru \cdot (-n \times \nabla \times S_{k}  + 
!          \alpha n \times n \times S_{k})[j_{v} \rv]
!        * wnear(3,:) stores corrections for
!            \ru \cdot (-\alpha n \times n \times \nabla S_{k}[\rho])        
!        * wnear(4,:) stores corrections for
!          \rv \cdot (-n \times \nabla \times S_{k}  + 
!          \alpha n \times n \times S_{k})[j_{u} \ru]
!        * wnear(5,:) stores corrections for
!          \rv \cdot (-n \times \nabla \times S_{k}  + 
!          \alpha n \times n \times S_{k})[j_{v} \rv]
!        * wnear(6,:) stores corrections for
!          \rv \cdot (-\alpha n \times n \times \nabla S_{k}[\rho])
!        * wnear(7,:) stores corrections for
!         (-ik n \cdot S_{k} + \alpha \nabla \cdot S_{k})[j_{u} \ru]
!        * wnear(8,:) stores corrections for
!         (-ik n \cdot S_{k} + \alpha \nabla \cdot S_{k})[j_{v} \rv]
!        * wnear(9,:) stores corrections for S_{k}'[\rho]        
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
!        number of densities (must be 3)
!    - sigma: complex *16(ndim, npts)
!        * sigma(1,:) is the ru component of J 
!        * sigma(2,:) is the rv component of J
!        * sigma(3,:) is rho 

!
!  Output arguments
!    - pot: complex *16 (ndim, npts)
!        The application of the integral representation 
!        (excluding the identity term)
!        * pot(1,:) is the ru component of 
!           n x H + \alpha n x n x E  equation
!        * pot(2,:) is the rv component of 
!           n x H + \alpha n x n x E  equation
!        * pot(3,:) is the n . E + alpha S_{k}(\div J - ik rho)
!

      implicit none
      integer, intent(in) :: npatches, npts
      integer, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps
      integer, intent(in) :: ndd, ndi, ndz
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      integer, intent(in) :: ipars(ndi)
      integer, intent(in) :: nnz, nquad
      integer, intent(in) :: row_ptr(npts+1), col_ind(nnz)
      integer, intent(in) :: iquad(nnz+1)
      
      integer, intent(in) :: ndim
      complex *16, intent(in) :: sigma(ndim,npts)
      
      integer, intent(in) :: nptso
      integer, intent(in) :: novers(npatches+1)
      integer, intent(in) :: ixyzso(npatches+1)
      real *8, intent(in) :: srcover(12,nptso), whtsover(nptso)
      
      integer, intent(in) :: nker
      complex *16 wnear(nker, nquad)

      integer, intent(in) :: lwork
      real *8, intent(in) :: work(lwork)

      complex *16, intent(out) :: pot(ndim,npts)

!
!  temporary variables
!
      integer nduse 
      real *8, allocatable :: ru(:,:), rv(:,:)
      real *8, allocatable :: ruover(:,:), rvover(:,:)

      complex *16, allocatable :: zjvec(:,:), charges(:,:)
      complex *16, allocatable :: zpot(:,:), zgrad(:,:,:)
      complex *16, allocatable :: sigmaover(:,:)
      integer nmax
      real *8 thresh,ra
      real *8, allocatable :: sources(:,:)

      complex *16 alpha, zk 

      real *8, allocatable :: srctmp(:,:), srctmp2(:,:)
      complex *16, allocatable :: ctmp(:,:)
      complex *16 zpottmp(4), zgradtmp(4,3), pottmp(3)

      integer ns, nt, npols
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      integer i,j,jpatch,jquadstart,jstart

      real *8 timeinfo(10),t1,t2
!$    real *8 omp_get_wtime      

      real *8 over4pi
      integer nss,ii,l,npover, m, n, iind, ier

      integer nd,ntarg0

      real *8 ttot,done,pi

      parameter (nd=1,ntarg0=1)
      data over4pi/0.07957747154594767d0/


      ns = nptso
      done = 1
      pi = atan(done)*4
           
      ifpgh = 0
      ifpghtarg = 1
      nduse = 4
      allocate(sources(3,ns), srctmp(3,npts))
      allocate(charges(nduse,ns), zjvec(3,ns))
      allocate(zpot(nduse,npts), zgrad(nduse,3,npts))
      allocate(sigmaover(3,ns))
      allocate(ru(3,npts), rv(3,npts))
      allocate(ruover(3,nptso), rvover(3,nptso))

      call orthonormalize_all(srcvals(4:6,:), srcvals(10:12,:), ru, &
         rv, npts)
      call orthonormalize_all(srcover(4:6,:), srcover(10:12,:), &
         ruover, rvover, ns)

      call cpu_time(t1)
!$      t1=omp_get_wtime()      
      call get_near_corr_max(npts, row_ptr, nnz, col_ind, npatches, &
        ixyzso, nmax)
      call cpu_time(t2)
!$      t2=omp_get_wtime()     

      call cpu_time(t1)
!$      t1=omp_get_wtime()      

! 
!       oversample density
!

      call oversample_fun_surf(6, npatches, norders, ixyzs, iptype, & 
         npts, sigma, novers, ixyzso, ns, sigmaover)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,ns
         sources(1:3,i) = srcover(1:3,i)
         zjvec(1:3,i) = sigmaover(1,i)*ruover(1:3,i) + &
              sigmaover(2,i)*rvover(1:3,i)
         charges(1:3,i) = zjvec(1:3,i)*whtsover(i)*over4pi
         charges(4,i) = sigmaover(3,i)*whtsover(i)*over4pi
      enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,npts
         srctmp(1:3,i) = srcvals(1:3,i)
      enddo
!$OMP END PARALLEL DO


!
!        compute threshold for ignoring local computation
!
      call get_fmm_thresh(12, ns, srcover, 12, npts, srcvals, thresh)
      call cpu_time(t2)
!$      t2=omp_get_wtime()      
      
      ier = 0
      zk = zpars(1)
      alpha = zpars(2)
      call hfmm3d_t_c_g_vec(nduse, eps, zk, ns, sources, &
        charges, npts, srctmp, zpot, zgrad, ier)
!
! Convert cartesian fmm info into components of integral equation
!

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) 
      do i=1,npts
         call get_nrccie_inteq_comps_from_potgrad(zk, alpha, &
            srcvals(1,i), ru(1,i), rv(1,i), zpot(1,i), zgrad(1,1,i), &
            pot(1,i))
      enddo
!$OMP END PARALLEL DO

      
!
!       Add near field precomputed contribution
!
      call cpu_time(t1)
!$      t1 = omp_get_wtime()

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, jpatch, jquadstart) &
!$OMP PRIVATE(jstart, npols, l, m, n, iind)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            do m=1,3
              do n=1,3
                iind = (m-1)*3 + n
                pot(m,i) = pot(m,i) + wnear(iind,jquadstart+l-1)*  &
                             sigma(n,jstart+l-1)
              enddo
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

      call cpu_time(t2)
!$      t2 = omp_get_wtime()
      timeinfo(2) = t2-t1

!
!     Remove near contribution of the FMM
!
      allocate(srctmp2(3,nmax), ctmp(nduse,nmax))
      call cpu_time(t1)
!$      t1 = omp_get_wtime()
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, nss, j, jpatch, l)  &
!$OMP PRIVATE(srctmp2, ctmp, zpottmp, zgradtmp, pottmp)
      do i=1,npts
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l=ixyzso(jpatch), ixyzso(jpatch+1)-1
            nss = nss + 1
            srctmp2(1:3,nss) = sources(1:3,l)            
            ctmp(1:nduse,nss) = charges(1:nduse,l)
          enddo
        enddo
        zpottmp(1:4) = 0
        zgradtmp(1:4,1:3) = 0
        pottmp(1:3) = 0
        call h3ddirectcg(nduse, zk, srctmp2, ctmp, nss, &
          srctmp(1,i), ntarg0, zpottmp, zgradtmp, thresh)

        call get_nrccie_inteq_comps_from_potgrad(zk, alpha, &
            srcvals(1,i), ru(1,i), rv(1,i), zpottmp, zgradtmp, &
            pottmp)

        pot(1:3,i) = pot(1:3,i) - pottmp(1:3)
      enddo
!$OMP END PARALLEL DO      
      call cpu_time(t2)
!$      t2 = omp_get_wtime()     

      timeinfo(3) = t2-t1
      
      
      ttot = timeinfo(1) + timeinfo(2) +timeinfo(3)


      return
      end subroutine lpcomp_em_nrccie_pec_addsub
!
!
!
!
!

      subroutine get_nrccie_inteq_comps_from_potgrad(zk, alpha, srcvals, &
            ru, rv, zpot, zgrad, pot)
!
!  Given the potential and gradient of the vector helmholtz FMM which computes
!  the cartesian coordinates of S_{k}[J] and S_{k}[\rho], this subroutine
!  returns the components of the nrccie integral equation given by
!
!    -M_{k}[J]+alpha·nxnx(ikS_{k}[J]-gradS_{k}[rho])  -  ru component 
!    -M_{k}[J]+alpha·nxnx(ikS_{k}[J]-gradS_{k}[rho])  -  rv component 
!     S'_{k}[rho]-ikn \cdot S_{k}[J]+alpha(divS_{k}[J]-ikS_{k}[rho]) 
!  
!  where M_{k} = n \times \nabla \times S_{k} is the magnetic field 
!  integral equation operator
!  
!  Input Arguments:
!    - zk: complex *16
!         Wavenumber (k) above
!    - alpha: complex *16
!         the parameter alpha above determining the linear combination
!         of the tangential component of the electric field,
!         and continuity condition to be included in the 
!         integral equation
!    - srcvals: real *8(12)
!         target info, xyz coordinates, d/du (xyz), d/dv(xyz), and 
!         normals
!    - ru: real *8(3)
!         orthonormalized tangent vectors, d/du(xyz)/|d/du(xyz)|
!    - rv: real *8(3)
!         second orthonormalized tangent vector
!         ru \times normal
!    - zpot: complex *16(4)
!         the 4 components of the potential, 
!         S_{k}[J_{1}], S_{k}[J_{2}], S_{k}[J_{3}], S_{k}[\rho]
!    - zgrad: complex *16(4,3)
!         gradients of the potential above
!
!
!  Output arguments:
!    - pot: complex *16(3)
!        the components of the integral equation
!
!
!
      implicit none
      complex *16, intent(in) :: zk, alpha
      real *8, intent(in) :: srcvals(12), ru(3), rv(3)
      complex *16, intent(in) :: zpot(4), zgrad(4,3)
      complex *16, intent(out) :: pot(3)

      complex *16 zcurl(3), zvec(3), zvec2(3), zvec3(3)
      complex *16 ztmp, ztmp2, ima
      data ima/(0.0d0,1.0d0)/


      zcurl(1) = zgrad(3,2) - zgrad(2,3)
      zcurl(2) = zgrad(1,3) - zgrad(3,1)
      zcurl(3) = zgrad(2,1) - zgrad(1,2)

      call dzcross_prod3d(srcvals(10), zcurl, zvec)

      zvec2(1:3) = ima*zk*zpot(1:3) - zgrad(4,1:3)
      ztmp = zvec2(1)*srcvals(10) + zvec2(2)*srcvals(11) + &
               zvec2(3)*srcvals(12)
      zvec3(1:3) = ztmp*srcvals(10:12) - zvec2(1:3) 
      zvec3(1:3) = alpha*zvec3(1:3) - zvec(1:3)

      pot(1) = zvec3(1)*ru(1) + zvec3(2)*ru(2) + zvec3(3)*ru(3)
      pot(2) = zvec3(1)*rv(1) + zvec3(2)*rv(2) + zvec3(3)*rv(3)
         
      ztmp = (zgrad(4,1) - ima*zk*zpot(1))*srcvals(10) + &
             (zgrad(4,2) - ima*zk*zpot(2))*srcvals(11) + &
             (zgrad(4,3) - ima*zk*zpot(3))*srcvals(12)
      ztmp2 = zgrad(1,1) + zgrad(2,2) + zgrad(3,3) - ima*zk*zpot(4)

      pot(3) = ztmp + alpha*ztmp2

      return
      end
!
!
!
!
!
      subroutine em_nrccie_pec_solver(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, numit, einc, &
        hinc, eps_gmres, niter, errs, rres, zjvec, rho)
!
!
!  This subroutine solves the Scattering Maxwell p.e.c. problem.
!  The the equations are:
!
!    curlE=ikH; curlH =-ikE
!
!  Representation:
!
!    H=curlS_{k}[J]
!
!    E=ikS_{k}[J]-gradS_{k}[rho]
!
!  Boundary conditions:
!
!    (1) + alpha nx(2)         b.c.1
!
!    (3) + alpha (4)           b.c.2
!
!    where:
!
!    n x H + n x H_inc = J     (1)
!
!    n x E + n x E_inc = 0     (2)
!
!    n · E + n · E_inc = rho   (3)
!
!    (div E)/ik = 0                 (4)
!
!  The incoming fields must be 'compatible' (solutions to Maxwell's 
!  equations in free space)
!
!  Boundary integral equation:
!
!    J/2 - M_{k}[J] + alpha n x n x (ikS_{k}[J]-gradS_{k}[rho]) = 
!      = nxH_inc - alpha nxnxE_inc
!    rho/2 + S_{k}'[rho] - ik n . S_{k}[J] + alpha (div S_{k}[J]-ikS_{k}[rho]) = 
!      = n·E_inc
!
!  The linear system is solved iteratively using GMRES
!  until a relative residual of eps_gmres is reached
!
!
!  input:
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
!    eps - real *8
!      precision requested for computing quadrature and fmm
!      tolerance
!    zpars - complex *16 (3)
!      kernel parameters (Referring to formula (1))
!      zpars(1) = k 
!      zpars(2) = alpha
!    numit - integer
!      max number of gmres iterations        
!    einc - complex *16(3, npts)
!      incident electric field
!    hinc - complex *16(3, npts)
!      incident magnetic field
!    eps_gmres - real *8
!      gmres tolerance requested
!
!    output
!    - niter: integer
!        number of gmres iterations required for relative residual
!        to converge to eps_gmres
!    - errs:  real *8 (numit+1)
!        relative residual as a function of iteration
!        number (only errs(1:niter) is populated))
!    - rres: real *8
!        relative residual for computed solution
!    - zjvec: complex *16(3,npts)
!        the induced surface current
!    - rho: complex *16(npts)
!        the surface charge density          
!

      implicit none
      integer, intent(in) :: npatches, npts
      integer, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps, eps_gmres
      complex *16, intent(in) :: zpars(2)
      complex *16, intent(in) :: hinc(3,npts), einc(3,npts)
      integer, intent(in) :: numit

      real *8, intent(out) :: errs(numit+1), rres
      integer, intent(out) :: niter
      complex *16, intent(out) :: zjvec(3,npts), rho(npts)

      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_src(:,:)
    

      integer nover, npolso, nptso
      integer nnz, nquad
      integer, allocatable :: row_ptr(:), col_ind(:), iquad(:)

      complex *16, allocatable :: wnear(:,:)

      real *8, allocatable :: srcover(:,:), wover(:)
      integer, allocatable :: ixyzso(:), novers(:)

      real *8, allocatable :: cms(:,:), rads(:), rad_near(:) 

      integer i, j, jpatch, jquadstart, jstart

      integer ndtarg, nker
      real *8 ttot
      real *8 rfac, rfac0
      integer iptype_avg, norder_avg
      integer ikerorder, iquadtype, npts_over

!
!
!       gmres variables
!
      complex *16 zid
      integer iind, k, l
!
!
! 
      allocate(uvs_src(2, npts), ipatch_id(npts))
!C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,npts
        ipatch_id(i) = -1
        uvs_src(1,i) = 0
        uvs_src(2,i) = 0
      enddo
!C$OMP END PARALLEL DO   

!
!    initialize patch_id and uv_targ for on surface targets
!
      call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts, &
        ipatch_id, uvs_src)

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
      ndtarg = 12
      call findnearmem(cms, npatches, rad_near, ndtarg, srcvals, npts, &
        nnz)

      allocate(row_ptr(npts+1), col_ind(nnz))
      
      call findnear(cms, npatches, rad_near, ndtarg, srcvals, npts, &
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

      call get_far_order(eps, npatches, norders, ixyzs, iptype, cms, &
        rads, npts, srccoefs, ndtarg, npts, srcvals, ikerorder, &
        zpars(1), nnz, row_ptr, col_ind, rfac, novers, ixyzso)

      npts_over = ixyzso(npatches+1)-1
      

      allocate(srcover(12,npts_over), wover(npts_over))

      call oversample_geom(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, srcvals, novers, ixyzso, npts_over, srcover)

      call get_qwts(npatches, novers, ixyzso, iptype, npts_over, &
        srcover, wover)


!
!   compute near quadrature correction
!
      nquad = iquad(nnz+1)-1
      print *, "nquad=",nquad
      allocate(wnear(9,nquad))
      
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j)      
      do i=1,nquad
        do j=1,9
          wnear(j,i)=0
        enddo
      enddo
!$OMP END PARALLEL DO          


      iquadtype = 1

      call getnearquad_em_nrccie_pec(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, eps, zpars, &
        iquadtype, nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
        wnear)
 
      
      print *, "done generating near quadrature, now starting gmres"
      nker = 9

      call em_nrccie_pec_solver_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, numit, &
        einc, hinc, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, &
        novers, npts_over, ixyzso, srcover, wover, eps_gmres, niter, &
        errs, rres, zjvec, rho)

!
      return
      end subroutine em_nrccie_pec_solver
!
!
!
!
!                        
      subroutine em_nrccie_pec_solver_guru(npatches, norders, ixyzs, &
      iptype, npts, srccoefs, srcvals, eps, zpars, numit, &
      einc, hinc, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, &
      novers, nptso, ixyzso, srcover, whtsover, eps_gmres, niter, &
      errs, rres, zjvec, rho)
!
!
!  This subroutine solves the Scattering Maxwell p.e.c. problem.
!  The the equations are:
!
!    curlE=ikH; curlH =-ikE
!
!  Representation:
!
!    H=curlS_{k}[J]
!
!    E=ikS_{k}[J]-gradS_{k}[rho]
!
!  Boundary conditions:
!
!    (1) - alpha nx(2)         b.c.1
!
!    (3) + alpha (4)           b.c.2
!
!    where:
!
!    n x H + n x H_inc = J     (1)
!
!    n x E + n x E_inc = 0     (2)
!
!    n · E + n · E_inc = rho   (3)
!
!    (div E)/ik = 0                 (4)
!
!  The incoming fields must be 'compatible' (solutions to Maxwell's 
!  equations in free space)
!
!  Boundary integral equation:
!
!    J/2 - M_{k}[J] + alpha n x n x (ikS_{k}[J]-gradS_{k}[rho]) = 
!      = nxH_inc - alpha nxnxE_inc
!    rho/2 + S_{k}'[rho] - ik n . S_{k}[J] + alpha (div S_{k}[J]-ikS_{k}[rho]) = 
!      = n·E_inc
!
!  The linear system is solved iteratively using GMRES
!  until a relative residual of eps_gmres is reached.
!
!  This is the guru solver routine which on input
!  takes in precomputed near-quadrature corrections, 
!  and oversampled surface information for the far part                            
!  
!
!  input:
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
!    eps - real *8
!      precision requested for computing quadrature and fmm
!      tolerance
!    zpars - complex *16 (3)
!      kernel parameters (Referring to formula (1))
!      zpars(1) = k 
!      zpars(2) = alpha
!    numit - integer
!      max number of gmres iterations
!    einc - complex *16(3, npts)
!      incident electric field              
!    hinc - complex *16(3, npts)
!      incident magnetic field
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
!        Near field quadrature corrections for nrccie
!        * wnear(1,:) stores corrections for
!          \ru \cdot (-n \times \nabla \times S_{k}  + 
!          \alpha n \times n \times S_{k})[j_{u} \ru]
!        * wnear(2,:) stores corrections for
!          \ru \cdot (-n \times \nabla \times S_{k}  + 
!          \alpha n \times n \times S_{k})[j_{v} \rv]
!        * wnear(3,:) stores corrections for
!            \ru \cdot (-\alpha n \times n \times \nabla S_{k}[\rho])        
!        * wnear(4,:) stores corrections for
!          \rv \cdot (-n \times \nabla \times S_{k}  + 
!          \alpha n \times n \times S_{k})[j_{u} \ru]
!        * wnear(5,:) stores corrections for
!          \rv \cdot (-n \times \nabla \times S_{k}  + 
!          \alpha n \times n \times S_{k})[j_{v} \rv]
!        * wnear(6,:) stores corrections for
!          \rv \cdot (-\alpha n \times n \times \nabla S_{k}[\rho])
!        * wnear(7,:) stores corrections for
!         (-ik n \cdot S_{k} + \alpha \nabla \cdot S_{k})[j_{u} \ru]
!        * wnear(8,:) stores corrections for
!         (-ik n \cdot S_{k} + \alpha \nabla \cdot S_{k})[j_{v} \rv]
!        * wnear(9,:) stores corrections for S_{k}'[\rho]        
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
!    eps_gmres - real *8
!      gmres tolerance requested
!
!    output
!    - niter: integer
!        number of gmres iterations required for relative residual
!        to converge to eps_gmres
!    - errs:  real *8 (numit+1)
!        relative residual as a function of iteration
!        number (only errs(1:niter) is populated))
!    - rres: real *8
!        relative residual for computed solution
!    - zjvec: complex *16(3,npts)
!        the induced surface current
!    - rho: complex *16(npts)
!        the surface charge density          
!
      implicit none
      integer, intent(in) :: npatches, npts
      integer, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps
      complex *16, intent(in) :: zpars(2)
      integer, intent(in) :: numit
      complex *16, intent(in) :: hinc(3,npts), einc(3,npts)
      integer, intent(in) :: nnz, nquad, nker
      integer, intent(in) :: row_ptr(npts+1), col_ind(nnz), iquad(nnz+1)
      complex *16, intent(in) :: wnear(nker,nquad)
      integer, intent(in) :: nptso
      integer, intent(in) :: novers(npatches), ixyzso(npatches+1)
      real *8, intent(in) :: srcover(12,nptso), whtsover(nptso)
      real *8, intent(in) :: eps_gmres
      integer, intent(out) :: niter
      real *8, intent(out) :: errs(numit+1), rres
      complex *16, intent(out) :: zjvec(3,npts), rho(npts)
      
      integer ndd, ndi, ndz
      integer ipars
      real *8 dpars, work
      complex *16, allocatable :: rhs(:,:), soln(:,:)
      real *8, allocatable :: ru(:,:), rv(:,:), wts(:)
      complex *16 zvec(3), zvec2(3), zvec3(3), ztmp, zid

      integer i, j, k, l, ndim, idim, lwork

      procedure (), pointer :: fker
      external lpcomp_em_nrccie_pec_addsub

      allocate(rhs(3,npts), soln(3,npts))
      allocate(ru(3,npts), rv(3,npts), wts(npts))

      call orthonormalize_all(srcvals(4:6,:), srcvals(10:12,:), ru, &
         rv, npts)

!
!  Initialize rhs from e and h fields
!                 
   
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(zvec, zvec2, ztmp)
      do i=1,npts

!  Compute n \times hinc        
        call dzcross_prod3d(srcvals(10,i), hinc(1,i), zvec)
       
!  Compute n \times n \times einc        
        ztmp = einc(1,i)*srcvals(10,i) + einc(2,i)*srcvals(11,i) + &
               einc(3,i)*srcvals(12,i)
        zvec2(1:3) = ztmp*srcvals(10:12,i) - einc(1:3,i)
        
        zvec(1:3) = zvec(1:3) - zpars(2)*zvec2(1:3)
        rhs(1,i) = ru(1,i)*zvec(1) + ru(2,i)*zvec(2) + ru(3,i)*zvec(3)
        rhs(2,i) = rv(1,i)*zvec(1) + rv(2,i)*zvec(2) + rv(3,i)*zvec(3)
        rhs(3,i) = einc(1,i)*srcvals(10,i) + einc(2,i)*srcvals(11,i) + &
          einc(3,i)*srcvals(12,i)
      
        soln(1:3,i) = 0
      enddo
!$OMP END PARALLEL DO    
      
      ndd = 0
      ndz = 2
      ndi = 0
      lwork = 0

      fker => lpcomp_em_nrccie_pec_addsub

      ndim = 3

      call get_qwts(npatches, norders, ixyzs, iptype, npts, &
      srcvals, wts)
!
      zid = 0.5d0
      call zgmres_guru(npatches, norders, ixyzs, &
            iptype, npts, srccoefs, srcvals, wts, &
            eps, ndd, dpars, ndz, zpars, ndi, ipars, &
            nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
            nptso, ixyzso, srcover, whtsover, lwork, work, &
            ndim, fker, zid, rhs, numit, eps_gmres, niter, errs, &
            rres, soln)
!
!  Now construct current and charge densities from solution
!                         

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,npts
        zjvec(1:3,i) = soln(1,i)*ru(1:3,i) + soln(2,i)*rv(1:3,i)
        rho(i) = soln(3,i)
      enddo
!$OMP END PARALLEL DO            


      return
      end

!
!
!  Post processing routines      
!
!      

      subroutine getnearquad_em_nrccie_pec_eval(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)
!
!  This subroutine generates the near field quadrature
!  for the non-resonant charge current integral representation:
!
!    H = \nabla \times S_{k} [J]
!    E = ik S_{k} [J] - \nabla S_{k} [\rho]                       
!
!  The quadrature is computed by the following strategy
!  targets within a sphere of radius rfac0*rs
!  of a chunk centroid is handled using adaptive integration
!  where rs is the radius of the bounding sphere
!  for the patch
!  
!  All other targets in the near field are handled via
!  oversampled quadrature
!
!  The recommended parameter for rfac0 is 1.25d0
!
!  NOTES:
!  This subroutine returns 4 kernels
!  1) S_{k}
!  2) \partial_{x} S_{k}
!  3) \partial_{y} S_{k}
!  4) \partial_{z} S_{k}
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
!        uvs_targ(1:2,i) is unused otherwise
!    - eps: real *8
!        precision requested
!    - zpars: complex *16 (1)
!        kernel parameters
!          * zpars(1) = k 
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
!    - rfac0: integer
!        radius parameter for near field
!    - nquad: integer
!        number of near field entries corresponding to each source target
!        pair. The size of wnear is 9*nquad since there are 9 kernels
!        per source target pair
!
!    output
!      wnear - complex *16 (4, nquad) 
!        Near field quadrature corrections for nrccie
!        * wnear(1,:) stores corrections for S_{k}
!        * wnear(2,:) stores corrections for \partial_{x} S_{k}
!        * wnear(3,:) stores corrections for \partial_{y} S_{k}        
!        * wnear(4,:) stores corrections for \partial_{z} S_{k}        
!

      implicit none 
      integer, intent(in) :: npatches, norders(npatches), npts, nquad
      integer, intent(in) :: ixyzs(npatches+1), iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts), eps
      integer, intent(in) :: ndtarg, ntarg
      real *8, intent(in) :: targs(ndtarg, ntarg)
      integer, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      real *8, intent(in) :: rfac0
        
      integer, intent(in) :: iquadtype, nnz
      integer, intent(in) :: row_ptr(ntarg+1), col_ind(nnz)
      integer, intent(in) :: iquad(nnz+1)
        
      complex *16, intent(in) :: zpars(1)
      complex *16, intent(out) :: wnear(4, nquad)
        
        
      integer ipars(2)
      real *8 dpars(1)
        
      complex *16, allocatable :: wneartmp(:)
  
      complex *16 alpha, beta
      integer i, j, ndi, ndd, ndz, iuse, idim
  
      integer ipv
  
      procedure (), pointer :: fker
      external  h3d_slp, h3d_sgradx, h3d_sgrady, h3d_sgradz


!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, idim)
      do i=1,nquad
        do idim = 1,4
          wnear(idim,i) = 0
        enddo
      enddo
!$OMP END PARALLEL DO

      allocate(wneartmp(nquad))

      ndd = 0
      ndz = 1
      ndi = 0
      ipv = 0
      do j=1,4
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wneartmp(i) = 0
        enddo
!$OMP END PARALLEL DO             
        if (j.eq.1) fker => h3d_slp
        if (j.eq.2) fker => h3d_sgradx
        if (j.eq.3) fker => h3d_sgrady
        if (j.eq.4) fker => h3d_sgradz
        if (j.ge.2) ipv = 1
               
        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, &
          zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, &
          nquad, wneartmp)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wnear(j,i) = wneartmp(i)
        enddo
!$OMP END PARALLEL DO                
      enddo


      return
      end
!
!

!
!                  
      subroutine em_nrccie_eval_addsub(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        eps, ndd, dpars, ndz, zpars, ndi, ipars, &
        nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
        nptso, ixyzso, srcover, whtsover, lwork, work, idensflag, &
        ndim_s, sigma, ipotflag, ndim_p, pot)

!
!
!  This subroutine computes the Electric (E) and Magnetic (H) fields
!  for the representation:
!
!    H = \nabla \times S_{k} [J]
!    E = ik S_{k} [J] - \nabla S_{k} [\rho]                       
! 
!  where the near field is precomputed and stored
!  in the row sparse compressed format.
!        
!  For targets on surface, this routine returns the principal
!  value part of the various layer potentials. 
!  The user is responsible for adding the
!  contribution of the identity terms.            
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
!        integral representation 
!    - zpars: complex *16(ndz)
!        complex parameters defining the kernel/
!        integral represnetation.
!        * zpars(1) = k 
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
!        number of kernels in quadrature correction, must be 4
!    - wnear: complex *16(nker, nquad)
!        precomputed quadrature corrections
!        * wnear(1,:) stores corrections for S_{k}
!        * wnear(2,:) stores corrections for \nabla_{x_{1}} S_{k} 
!        * wnear(3,:) stores corrections for \nabla_{x_{2}} S_{k}
!        * wnear(3,:) stores corrections for \nabla_{x_{3}} S_{k}
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
!         unused for this routine
!    - ndim_s: integer
!        number of densities per point on the surface,
!        must be 4 for this routine
!    - sigma: complex *16(ndim_s, npts)
!        The densities
!        * sigma(1,:) stores the x_{1} component of J
!        * sigma(2,:) stores the x_{2} component of J
!        * sigma(3,:) stores the x_{3} component of J
!        * sigma(4,:) stores \rho
!    - ipotflag: integer
!        Flag for determining output type 
!        ipotflag = 1, only E is returned
!        ipotflag = 2, only H is returned
!        ipotflag = 3, both E, and H are returned   
!    - ndim_p: integer
!        number of potentials per point, must be 6
                
!    
!  Output arguments:
!    - pot: complex *16 (ndim_p, ntarg)
!        pot(1:3, :) is E if ipotflag = 1 or 3
!        pot(4:6, :) is H if ipotflag = 2 or 3  
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
      complex *16, intent(in) :: wnear(nker,nquad)

      integer, intent(in) :: nptso
      integer, intent(in) :: ixyzso(npatches+1), novers(npatches)
      real *8, intent(in) :: srcover(12,nptso), whtsover(nptso)
        
      integer, intent(in) :: lwork
      real *8, intent(in) :: work(lwork)

      integer, intent(in) :: ndim_s, ndim_p
      integer, intent(in) :: idensflag, ipotflag

      complex *16, intent(in) :: sigma(ndim_s,npts)

      complex *16, intent(out) :: pot(ndim_p,ntarg)

      complex *16, allocatable :: sigmause(:,:)
      integer idensflaguse

      complex *16 zparsuse(7), ima
      integer ifj, ifm, ifrho, iftau, ndim_s_use, i
      data ima/(0.0d0,1.0d0)/
      
      ifj = 1
      ifrho = 1
      ifm = 0
      iftau = 0
      
      idensflaguse = ifj + 2**(ifm)*ifm + 4**(ifrho)*ifrho + 8**(iftau)*iftau
      zparsuse(1) = zpars(1)
      zparsuse(2) = ima*zpars(1)
      zparsuse(3) = -1.0d0
      zparsuse(4) = 0
      zparsuse(5) = 0
      zparsuse(6) = 0
      zparsuse(7) = 1.0d0

      allocate(sigmause(8,npts))
      ndim_s_use = 8
!$OMP PARALLEL DO DEFAULT(SHARED)
      do i = 1,npts
        sigmause(1,i) = sigma(1,i)
        sigmause(2,i) = sigma(2,i)
        sigmause(3,i) = sigma(3,i)
        sigmause(4,i) = 0
        sigmause(5,i) = 0
        sigmause(6,i) = 0
        sigmause(7,i) = sigma(4,i)
        sigmause(8,i) = 0
      enddo
!$OMP END PARALLEL DO            
      
      call em_debye_eval_addsub(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        eps, ndd, dpars, ndz, zparsuse, ndi, ipars, &
        nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
        nptso, ixyzso, srcover, whtsover, lwork, work, idensflaguse, &
        ndim_s_use, sigmause, ipotflag, ndim_p, pot)

      return
      end
!
!      
!
!
!            
      subroutine em_nrccie_pec_eval(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, zpars, zjvec, rho, ife, e, ifh, h)
!
!
!
!  This subroutine computes the Electric (E) and Magnetic (H) fields
!  for the representation:
!
!    H = \nabla \times S_{k} [J]
!    E = ik S_{k} [J] - \nabla S_{k} [\rho]                       
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
!        ipatch_id(i) indicates the patch on which target i
!        is, if it is on surface. ipatch_id(i) should be 0 
!        otherwise
!    - uvs_targ: real *8(2,ntarg)
!        if ipatch_id(i) > 0, then uvs_targ(1:2,i) are the
!        local uv coordinates of the target on the patch,
!        uvs_targ(1:2,i) is unused otherwise
!    - eps: real *8
!        precision requested
!    - zpars: complex *16(1)
!        complex parameters defining the kernel/
!        integral represnetation.
!        * zpars(1) = k 
!    - zjvec: complex *16(3,npts)
!        current J above
!    - rho: complex *16 (npts)
!        charge rho above
!    - ife: integer
!        flag for computing the e field, electric field will be returned
!        if ife = 1                        
!    - ifh: integer
!        flag for computing the h field, magnetic field will be returned
!        if ifh = 1
!                                        
!  Output arguments:
!    - e: complex *16(3,ntarg)
!        electric field at the target locations, if requested
!    - h: complex *16(3,ntarg)
!        magnetic field at the target locations, if requested				 
!
      implicit none
      integer, intent(in) :: npatches
      integer, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches), npts
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      integer, intent(in) :: ndtarg, ntarg
      real *8, intent(in) :: targs(ndtarg,ntarg)
      integer, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg), eps
      complex *16, intent(in) :: zpars(1), zjvec(3,npts), rho(npts)
      integer, intent(in) :: ife, ifh
      complex *16, intent(out) :: e(3,ntarg), h(3,ntarg)
      complex *16, allocatable :: sigmause(:,:), potuse(:,:)
      integer nnz, nquad
      integer, allocatable :: row_ptr(:), col_ind(:), iquad(:)
  
      complex *16, allocatable :: wnear(:,:)
  
      real *8, allocatable :: srcover(:,:), wover(:)
      integer, allocatable :: ixyzso(:), novers(:)
  
      real *8, allocatable :: cms(:,:), rads(:), rad_near(:) 

      
  
      integer i, j, jstart
  
      integer ipars(1)
      real *8 dpars(1), work(1)
      integer ndim_s, ndim_p
      real *8 t1, t2
  !$   real *8 omp_get_wtime
  
      integer ndd, ndi, ndz, lwork
      integer nker
      real *8 rfac, rfac0
      integer iptype_avg, norder_avg
      integer ikerorder, iquadtype, npts_over
      complex *16 ima
      integer ndtarg0
      integer idensflag, ipotflag

      ima=(0.0d0,1.0d0)
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

      allocate(row_ptr(ntarg+1),col_ind(nnz))

      
      call findnear(cms, npatches, rad_near, ndtarg, targs, ntarg, &
        row_ptr, col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches, ixyzs, ntarg, nnz, row_ptr, col_ind,&
         iquad)
      nquad = iquad(nnz+1) - 1
      nker = 4
      allocate(wnear(nker,nquad))

      iquadtype = 1
    
      call getnearquad_em_nrccie_pec_eval(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)

      ikerorder = 0
!
!    estimate oversampling for far-field, and oversample geometry
!

      allocate(novers(npatches),ixyzso(npatches+1))

      print *, "beginning far order estimation"

      call get_far_order(eps, npatches, norders, ixyzs, iptype, cms, &
       rads, npts, srccoefs, ndtarg, ntarg, targs, ikerorder, zpars(1), &
       nnz, row_ptr, col_ind, rfac, novers, ixyzso)
      
      npts_over = ixyzso(npatches+1)-1

      allocate(srcover(12,npts_over), wover(npts_over))

      call oversample_geom(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, srcvals, novers, ixyzso, npts_over, srcover)

      call get_qwts(npatches, novers, ixyzso, iptype, npts_over, &
        srcover, wover)

      ndz = 1
      ndd = 0 
      ndi = 0  

      if(ife.eq.1.and.ifh.eq.0) ipotflag = 1
      if(ife.eq.0.and.ifh.eq.1) ipotflag = 2
      if(ife.eq.1.and.ifh.eq.1) ipotflag = 3
      ndim_s = 4
      ndim_p = 6
      allocate(sigmause(ndim_s,npts), potuse(ndim_p,ntarg))

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i = 1,npts
        sigmause(1:3,i) = zjvec(1:3,i)
        sigmause(4,i) = rho(i)
      enddo
!$OMP END PARALLEL DO         
      call em_nrccie_eval_addsub(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        eps, ndd, dpars, ndz, zpars, ndi, ipars, &
        nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
        npts_over, ixyzso, srcover, wover, lwork, work, idensflag, &
        ndim_s, sigmause, ipotflag, ndim_p, potuse)

      if(ife.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i = 1,ntarg
          e(1:3,i) = potuse(1:3,i)
        enddo
!$OMP END PARALLEL DO                
      endif

      if(ifh.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i = 1,ntarg
          h(1:3,i) = potuse(4:6,i)
        enddo
!$OMP END PARALLEL DO                
      endif  

      return
      end

!
!
!
!  CONTINUE FROM HERE
!
subroutine fker_em_nrccie_pec_s(srcinfo, ndt,targinfo,ndd, dpars,ndz,zpars,&
 &ndi,ipars,E_val)
implicit none

!  This function provides the near field kernel that will use
!  zgetnearquad_ggq_guru 
!  through getnearquad_em_nrccie_pec

    !List of calling arguments
    integer, intent(in) :: ndt,ndd,ndz,ndi
    real ( kind = 8 ), intent(in) :: srcinfo(12)
    real ( kind = 8 ), intent(in) :: targinfo(ndt)
    integer, intent(in) :: ipars(ndi)
    real ( kind = 8 ), intent(in) :: dpars(ndd)
    complex ( kind = 8 ), intent(in) :: zpars(ndz)
    complex ( kind = 8 ), intent(out) :: E_val

    
    !List of local variables
    complex ( kind = 8 ) E_mat(3,3)
    real ( kind = 8 ) ru_s(3),rv_s(3),n_s(3)
    real ( kind = 8 ) ru_t(3),rv_t(3),n_t(3)
    real ( kind = 8 ) du(3), dv(3),sour(3),targ(3)
    real ( kind = 8 ) r, dr(3),aux_real
    complex ( kind = 8 ) nxcurlSka(2,2),Dkrho,nxnxSkb(2,2)
    complex ( kind = 8 ) nxgradSklambda(2,1),nxnxgradSklambda(2,1),nSkb(1,2)
    complex ( kind = 8 ) Sklambda,divSkb(1,2),ngradSklambda
    real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3),xprod_aux3(3),xprod_aux4(3)
    complex ( kind = 8 ) R1, ima,my_exp, zk,alpha
    real ( kind = 8 ) pi
    
    pi=3.1415926535897932384626433832795028841971d0
    ima=(0.0d0,1.0d0)
    zk=zpars(1)
    alpha=zpars(2)
        
    sour(1)=srcinfo(1)
    sour(2)=srcinfo(2)
    sour(3)=srcinfo(3)
    
    n_s(1)=srcinfo(10)
    n_s(2)=srcinfo(11)
    n_s(3)=srcinfo(12)

    targ(1)=targinfo(1)
    targ(2)=targinfo(2)
    targ(3)=targinfo(3)

    n_t(1)=targinfo(10)
    n_t(2)=targinfo(11)
    n_t(3)=targinfo(12)

    dr(1)=targ(1)-sour(1)
    dr(2)=targ(2)-sour(2)
    dr(3)=targ(3)-sour(3)
    r=sqrt((dr(1))**2+(dr(2))**2+(dr(3))**2)
    
    R1=(ima*zk*r-1.0d0)/r**3*exp(ima*zk*r)/(4.0d0*pi)
    my_exp=exp(ima*zk*r)/(4.0d0*pi)

    call orthonormalize(srcinfo(4:6),n_s,ru_s,rv_s)
    call orthonormalize(targinfo(4:6),n_t,ru_t,rv_t)

    call get_nxcurlSka(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r,nxcurlSka)
    call get_nxnxSkb(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,my_exp,r,nxnxSkb)
    E_mat(1,1)=-nxcurlSka(1,1)+ima*zk*alpha*nxnxSkb(1,1)
    E_mat(1,2)=-nxcurlSka(1,2)+ima*zk*alpha*nxnxSkb(1,2)
    E_mat(2,1)=-nxcurlSka(2,1)+ima*zk*alpha*nxnxSkb(2,1)
    E_mat(2,2)=-nxcurlSka(2,2)+ima*zk*alpha*nxnxSkb(2,2)
    call get_nxnxgradSklambda(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r, &
      nxnxgradSklambda)
    E_mat(1,3)=-alpha*nxnxgradSklambda(1,1)
    E_mat(2,3)=-alpha*nxnxgradSklambda(2,1)
    call get_ngradSklambda(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r, &
      ngradSklambda)
    call get_Sklambda(my_exp,r,Sklambda)
    E_mat(3,3)=ngradSklambda-ima*zk*alpha*Sklambda
    call get_divSkb(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r,divSkb)
    call get_nSkb(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,my_exp,r,nSkb)
    E_mat(3,1)=-ima*zk*nSkb(1,1)+alpha*divSkb(1,1)
    E_mat(3,2)=-ima*zk*nSkb(1,2)+alpha*divSkb(1,2)

    E_val=E_mat(ipars(1),ipars(2))

return
end subroutine fker_em_nrccie_pec_s


subroutine em_nrccie_pec_FMM(eps,zpars,ns,nt,srcvals,targvals,wts,a_u,a_v,&
 &rho_in,AA_u,AA_v,PHI,thresh,ifdir)
implicit none

!
!  This subroutine computes the far field contribution if the NRCCIE 
!  operator via FMM
!
!  Boundary integral equation:
!
!    J/2-M_{k}[J]+alpha·nxnx(ikS_{k}[J]-gradS_{k}[rho]) =
!         = nxH_inc + alpha nxnxE_inc  ->  (AA_u,AA_v)
!    rho/2+S'_{k}[rho]-ikS_{k}[J]+alpha(divS_{k}[J]-ikS_{k}[rho]) =
!         = n·E_inc  ->  PHI
! 
!
!  input:
!    eps - real * 8
!      epsilon for the fmm call
!
!    zpars - complex *16(3)
!      kernel parameters
!      zpars(1) = k 
!      zpars(2) = alpha
!      zpars(3) = - (not used)
!
!    ns - integer
!      number of sources (this is the oversampled set of sources)
!   
!    nt - integer 
!      number of targets (this is the not oversampled set of targets)
!
!    srcvals - real *8(12,ns) (oversampled surface)
!      xyz(u,v) and derivative info sampled at the 
!      discretization nodes on the surface
!      srcvals(1:3,i) - xyz info
!      srcvals(4:6,i) - dxyz/du info
!      srcvals(7:9,i) - dxyz/dv info
!
!    targvals - real *8(3,nt) (not oversampled surface)
!      xyz(u,v) and derivative info sampled at the 
!      discretization nodes on the surface
!      srcvals(1:3,i) - xyz info
!      srcvals(4:6,i) - dxyz/du info
!      srcvals(7:9,i) - dxyz/dv info
!
!    wts - real *8(ns)
!      smooth quadrature weights at oversampled nodes
!
!    a_u,a_v - complex *16(ns)
!      two components of the tangent induced current J on the surface
!      along srcvals(4:6,i) and (srcvals(10:12,i) x srcvals(4:6,i))
!      directions
! 
!    rho_in - complex *16(ns)
!      induced charge on the surface
!
!    thresh - real *8
!      threshold to remove the selfo interaction term
!
!    ifdir - integer
!      flag, ifdir=1 direct calculation N^2 (used to remove teh near terms)
!            ifdir=0 FMM activated
!
!  output:
!    AA_u - complex *16(nt)
!      first component of -M_{k}[J]+alpha·nxnx(ikS_{k}[J]-gradS_{k}[rho])
!      along the srcvals(4:6,i) direction
!    AA_v - complex *16(nt)
!      second component of -M_{k}[J]+alpha·nxnx(ikS_{k}[J]-gradS_{k}[rho])
!      along the (srcvals(10:12,i) x srcvals(4:6,i)) direction
!    PHI - complex *16(nt)
!      scalar component of the operator:
!      S'_{k}[rho]-ikS_{k}[J]+alpha(divS_{k}[J]-ikS_{k}[rho])
!            

    !List of calling arguments
	real ( kind = 8 ), intent(in) :: eps
	complex ( kind = 8 ), intent(in) :: zpars(3)
    integer, intent(in) :: ns,nt
	real ( kind = 8 ), intent(in) :: srcvals(12,ns),wts(ns),targvals(12,nt)
    complex ( kind = 8 ), intent(in) :: a_u(ns),a_v(ns),rho_in(ns)
    complex ( kind = 8 ), intent(out) :: AA_u(nt),AA_v(nt),PHI(nt)
	real ( kind = 8 ), intent(in) :: thresh
	integer, intent(in) :: ifdir 


    !List of local variables
	real ( kind = 8 ), allocatable :: n_vect_s(:,:),u_vect_s(:,:)
	real ( kind = 8 ), allocatable :: v_vect_s(:,:),source(:,:)

	real ( kind = 8 ), allocatable :: n_vect_t(:,:),u_vect_t(:,:)
	real ( kind = 8 ), allocatable :: v_vect_t(:,:),targets(:,:)

    complex ( kind = 8 ), allocatable :: a_vect(:,:),b_vect(:,:),lambda(:)
    complex ( kind = 8 ), allocatable :: E(:,:),curlE(:,:),divE(:),rho(:)
    complex ( kind = 8 ), allocatable :: b_vect_t(:,:)
	complex ( kind = 8 ) ima,zk,alpha

    integer count1,count2,ier
    integer ifa_vect,ifb_vect,iflambda,ifrho,ifE,ifcurlE,ifdivE
	real ( kind = 8 ) pi
	pi=3.1415926535897932384626433832795028841971d0

	ima=(0.0d0,1.0d0)
	zk=zpars(1)
	alpha=zpars(2)

    allocate(a_vect(3,ns))
    allocate(b_vect(3,ns))
    allocate(lambda(ns))
	allocate(rho(ns))
    allocate(E(3,nt))
    allocate(curlE(3,nt))
    allocate(divE(nt))
	allocate(n_vect_s(3,ns))
	allocate(n_vect_t(3,nt))
	allocate(u_vect_s(3,ns))
	allocate(v_vect_s(3,ns))
	allocate(u_vect_t(3,nt))
	allocate(v_vect_t(3,nt))
	allocate(source(3,ns))
	allocate(targets(3,nt))
    allocate(b_vect_t(3,nt))

	do count1=1,ns
     n_vect_s(:,count1)=srcvals(10:12,count1)
     source(:,count1)=srcvals(1:3,count1)
	enddo
	call orthonormalize_all(srcvals(4:6,:),srcvals(10:12,:),u_vect_s,&
	 &v_vect_s,ns)
	
	do count1=1,nt
     n_vect_t(:,count1)=targvals(10:12,count1)
     targets(:,count1)=targvals(1:3,count1)
	enddo
	call orthonormalize_all(targvals(4:6,:),targvals(10:12,:),u_vect_t,&
	 &v_vect_t,nt)
	
	
	do count1=1,ns
     b_vect(1,count1)=(a_u(count1)*u_vect_s(1,count1)+a_v(count1)*&
      &v_vect_s(1,count1))
     b_vect(2,count1)=(a_u(count1)*u_vect_s(2,count1)+a_v(count1)*&
      &v_vect_s(2,count1))
     b_vect(3,count1)=(a_u(count1)*u_vect_s(3,count1)+a_v(count1)*&
      &v_vect_s(3,count1))
    enddo

    !Computing the full operator
    ifa_vect=0
    ifb_vect=1
    iflambda=0
    ifrho=0
    ifE=1
    ifcurlE=1
    ifdivE=1

	call Vector_Helmholtz_targ2(eps,zk,ns,source,wts,ifa_vect,a_vect,&
	 &ifb_vect,b_vect,iflambda,lambda,ifrho,rho,n_vect_s,ifE,E,ifcurlE,&
	 &curlE,ifdivE,divE,nt,targets,thresh,ifdir)

    do count1=1,nt
        b_vect_t(1,count1)=n_vect_t(2,count1)*curlE(3,count1)-&
		 &n_vect_t(3,count1)*curlE(2,count1)
        b_vect_t(2,count1)=n_vect_t(3,count1)*curlE(1,count1)-&
		 &n_vect_t(1,count1)*curlE(3,count1)
        b_vect_t(3,count1)=n_vect_t(1,count1)*curlE(2,count1)-&
		 &n_vect_t(2,count1)*curlE(1,count1)
    enddo

    do count1=1,nt
     AA_u(count1)=-alpha*ima*zk*(E(1,count1)*u_vect_t(1,count1)+E(2,count1)*&
      &u_vect_t(2,count1)+E(3,count1)*u_vect_t(3,count1))
     AA_v(count1)=-alpha*ima*zk*(E(1,count1)*v_vect_t(1,count1)+E(2,count1)*&
      &v_vect_t(2,count1)+E(3,count1)*v_vect_t(3,count1))
    enddo

    do count1=1,nt
     AA_u(count1)=AA_u(count1)-(b_vect_t(1,count1)*u_vect_t(1,count1)+&
      &b_vect_t(2,count1)*u_vect_t(2,count1)&
      &+b_vect_t(3,count1)*u_vect_t(3,count1))
     AA_v(count1)=AA_v(count1)-(b_vect_t(1,count1)*v_vect_t(1,count1)+&
      &b_vect_t(2,count1)*v_vect_t(2,count1)&
      &+b_vect_t(3,count1)*v_vect_t(3,count1))
    enddo

    do count1=1,nt
     PHI(count1)=-ima*zk*(E(1,count1)*n_vect_t(1,count1)+E(2,count1)*&
      &n_vect_t(2,count1)+E(3,count1)*n_vect_t(3,count1))
     PHI(count1)=PHI(count1)+alpha*divE(count1)
	enddo

	do count1=1,ns
     rho(count1)=rho_in(count1)*wts(count1)
	enddo

	do count1=1,nt
     divE(count1)=0.0d0
     E(1,count1)=0.0d0
     E(2,count1)=0.0d0
     E(3,count1)=0.0d0
	enddo

    !Computing the full operator
	if (ifdir.eq.1) then
     call h3ddirectcg(1,zk,source,rho,ns,targets,nt,divE,E&
     &,thresh)
    else
     call hfmm3d_t_c_g_vec(1,eps,zk,ns,source,rho,nt,targets,divE,&
     &E,ier)
	endif		
	
	do count1=1,nt
     E(:,count1)=E(:,count1)/(4.0d0*pi)
     divE(count1)=divE(count1)/(4.0d0*pi)
	enddo

    do count1=1,nt
     AA_u(count1)=AA_u(count1)+alpha*(E(1,count1)*u_vect_t(1,count1)+&
      &E(2,count1)*u_vect_t(2,count1)+E(3,count1)*u_vect_t(3,count1))
     AA_v(count1)=AA_v(count1)+alpha*(E(1,count1)*v_vect_t(1,count1)+&
      &E(2,count1)*v_vect_t(2,count1)+E(3,count1)*v_vect_t(3,count1))
    enddo

    do count1=1,nt
     PHI(count1)=PHI(count1)-alpha*ima*zk*divE(count1)+(E(1,count1)*&
      &n_vect_t(1,count1)+E(2,count1)*n_vect_t(2,count1)+E(3,count1)*&
      &n_vect_t(3,count1))
	enddo

	deallocate(a_vect)
	deallocate(b_vect)
	deallocate(lambda)
	deallocate(rho)
	deallocate(E)
	deallocate(curlE)
	deallocate(divE)
	
	deallocate(u_vect_s)
	deallocate(v_vect_s)
	deallocate(n_vect_s)
	deallocate(source)

	deallocate(u_vect_t)
	deallocate(v_vect_t)
	deallocate(n_vect_t)
	deallocate(targets)
    deallocate(b_vect_t)

return
end subroutine em_nrccie_pec_FMM


subroutine test_accuracy_em_nrccie_pec(eps_FMM,sol,zpars,ns,wts,srcvals,P0,&
 &vf,Pt)
implicit none

!
!  This function test the accuracy of the solution computed in the exterior
!  region by testing the extintion theorem in the interior region.
!
!  input:
!    eps_FMM - real *8
!      epsilon for the fmm call
!
!    sol - complex *16(3*ns)
!      induced charge and current on the surface
!      sol(1:ns) - first component of  J along
!      the srcvals(4:6,i) direction
!      sol(ns+1:2*ns) - second component of J along
!      the (srcvals(10:12,i) x srcvals(4:6,i)) direction
!      sol(2*ns+1:3*ns) - electric charge density rho
!
!    zpars - complex *16(3)
!      kernel parameters
!      zpars(1) = k 
!      zpars(2) = alpha
!      zpars(3) = - (not used)
!
!    ns - integer
!      number of sources
!
!    wts - real *8(ns)
!      smooth quadrature weights at original nodes
!
!    srcvals - real *8(12,ns)
!      xyz(u,v) and derivative info sampled at the 
!      discretization nodes on the surface
!      srcvals(1:3,i) - xyz info
!      srcvals(4:6,i) - dxyz/du info
!      srcvals(7:9,i) - dxyz/dv info
!
!    P0 - real * 8(3)
!      location of the source point at the INTERIOR region
!            
!    vf - complex *16(3)
!      Orientation of the magnetic and electric dipoles located at Pt 
!
!    Pt - real * 8(3)
!      location of the source point at the EXTERIOR region
!

    !List of calling arguments
    complex ( kind = 8 ), intent(in) :: zpars(3)
    integer, intent(in) :: ns
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),eps_FMM
    real ( kind = 8 ), intent(in) :: wts(ns),P0(3),Pt(3)
    complex ( kind = 8 ), intent(in) :: sol(3*ns),vf(3)

!   List of local variables
    complex ( kind = 8 ) a_u00,a_v00,zk,alpha
    complex ( kind = 8 ) ima, R1, R2,Et1(3),Et2(3),aux_cmp,Ht2(3),Ht1(3)
    real ( kind = 8 ) ru_s(3),rv_s(3),n_s(3),sour(3),r,dr(3)
    real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3),error_E,error_H
    real ( kind = 8 ) pi

    integer count1

    ima=(0.0d0,1.0d0)
    pi=3.1415926535897932384626433832795028841971d0
    zk=zpars(1)
    alpha=zpars(2)

    write (*,*) 'P0',P0
    call em_nrccie_pec_FMM_targ(eps_FMM,zk,ns,srcvals,1,P0,wts,sol(1:ns),&
      sol(ns+1:2*ns),sol(2*ns+1:3*ns),Et1,Ht1)
      
    call fieldsED(zk,Pt,P0,1,Et2,Ht2,vf,0)
    call fieldsMD(zk,Pt,P0,1,Et2,Ht2,vf,1)

    call prin2('Et2=*',Et2,6)
    call prin2('Ht2=*',Ht2,6)

!	
!   Here we are testing the extintion theorem,
!   that's why we ADD incoming and scattered fields.
!

    error_E=sqrt(abs(Et1(1)+Et2(1))**2+abs(Et1(2)+Et2(2))**2+&
        abs(Et1(3)+Et2(3))**2)
!	write (*,*) 'Error E: ', error_E
    write (*,*) 'Relative Error E: ', error_E/&
     sqrt(abs(Et2(1))**2+abs(Et2(2))**2+abs(Et2(3))**2)

    error_H=sqrt(abs(Ht1(1)+Ht2(1))**2+abs(Ht1(2)+Ht2(2))**2+&
      abs(Ht1(3)+Ht2(3))**2)
    write (*,*) 'Relative Error H: ', error_H/ &
        sqrt(abs(Ht2(1))**2+abs(Ht2(2))**2+abs(Ht2(3))**2)

  return
end subroutine test_accuracy_em_nrccie_pec
!
!
!
!
!
subroutine em_nrccie_pec_FMM_targ(eps,zk,ns,srcvals,nt,targ,wts,a_u,a_v,&
    rho_in,E,H)
implicit none
!
!  This funciton computes the fields E,H at a given point using the FMM
!  It doesn't contain near field corrections (it's for debugging purposes)   
!  The representation for the potentials is:
!  Representation:
!
!    H=curlS_{k}[J]
!
!    E=ikS_{k}[J]-gradS_{k}[rho]
!
!  it requires as input \rho and \rho_2:=S_{ik}[\rho]
!
!
!  input:
!    eps - real *8
!      epsilon for the fmm call
!
!    zk - complex *16
!      Helmholtz parameter 
!
!    ns - integer
!      number of sources
!   
!    nt - integer
!      number of targets
!
!    srcvals - real *8(12,ns)
!      xyz(u,v) and derivative info sampled at the 
!      discretization nodes on the surface
!      srcvals(1:3,i) - xyz info
!      srcvals(4:6,i) - dxyz/du info
!      srcvals(7:9,i) - dxyz/dv info
!
!    targ - real *8(3,nt)
!      location x,y,z of the target points
!
!    wts - real *8(ns)
!      smooth quadrature weights at original nodes
!
!    a_u,a_v - complex *16(ns)
!      two components of the tangent induced current J on the surface
!      along srcvals(4:6,i) and (srcvals(10:12,i) x srcvals(4:6,i)) directions
! 
!    rho_in - complex *16(ns)
!      induced charge on the surface
!
!  output:
!    E - complex  *16(3,nt)
!      value of the Electric field at the target points
!
!    H - complex  *16(3,nt)
!      value of the Magnetic field at the target points
!

    !List of calling arguments
    real ( kind = 8 ), intent(in) :: eps
    complex ( kind = 8 ), intent(in) :: zk
    integer, intent(in) :: ns,nt
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),targ(3,nt)
    real ( kind = 8 ), intent(in) :: wts(ns)
    complex ( kind = 8 ), intent(in) :: a_u(ns),a_v(ns),rho_in(ns)
    complex ( kind = 8 ), intent(out) :: E(3,nt),H(3,nt)

    !List of local variables
    real ( kind = 8 ), allocatable :: n_vect(:,:),u_vect(:,:),v_vect(:,:)
    real ( kind = 8 ), allocatable :: source(:,:)

    complex ( kind = 8 ), allocatable :: a_vect(:,:),b_vect(:,:)
    complex ( kind = 8 ), allocatable :: lambda(:),rho(:)
    complex ( kind = 8 ), allocatable :: curlE(:,:),gradpot(:,:),divE(:)
    complex ( kind = 8 ) ima

    integer count1,count2,ier
    integer ifa_vect,ifb_vect,iflambda,ifrho,ifE,ifcurlE,ifdivE
    real ( kind = 8 ) pi
    pi=3.1415926535897932384626433832795028841971d0

    ima=(0.0d0,1.0d0)

    allocate(a_vect(3,ns))
    allocate(b_vect(3,ns))
    allocate(lambda(ns))
    allocate(rho(ns))
    allocate(curlE(3,nt))
    allocate(gradpot(3,nt))
    allocate(divE(nt))

    allocate(n_vect(3,ns))
    allocate(u_vect(3,ns))
    allocate(v_vect(3,ns))
    allocate(source(3,ns))

    do count1=1,ns
      n_vect(:,count1)=srcvals(10:12,count1)
      source(:,count1)=srcvals(1:3,count1)
    enddo
    call orthonormalize_all(srcvals(4:6,:),srcvals(10:12,:),u_vect,v_vect,ns)

    do count1=1,ns
     b_vect(1,count1)=a_u(count1)*u_vect(1,count1)+&
      &a_v(count1)*v_vect(1,count1)
     b_vect(2,count1)=a_u(count1)*u_vect(2,count1)+&
      &a_v(count1)*v_vect(2,count1)
     b_vect(3,count1)=a_u(count1)*u_vect(3,count1)+&
      &a_v(count1)*v_vect(3,count1)
    enddo

!  Computing the full operator
    ifa_vect=0
    ifb_vect=1
    iflambda=0
    ifrho=0
    ifE=1
    ifcurlE=1
    ifdivE=0

    call Vector_Helmholtz_targ(eps,zk,ns,source,wts,ifa_vect,a_vect,&
     ifb_vect,b_vect,iflambda,lambda,ifrho,rho,n_vect,ifE,E,ifcurlE,H,&
     ifdivE,divE,nt,targ)

    do count1=1,ns
     rho(count1)=rho_in(count1)*wts(count1)
    enddo


    !Computing the full operator
    call hfmm3d_t_c_g(eps,zk,ns,source,rho,nt,targ,divE,gradpot,ier)


    do count1=1,nt
     gradpot(:,count1)=gradpot(:,count1)/(4.0d0*pi)
     divE(count1)=divE(count1)/(4.0d0*pi)
    enddo

    do count1=1,nt
      E(1,count1)=ima*zk*E(1,count1)-gradpot(1,count1)
      E(2,count1)=ima*zk*E(2,count1)-gradpot(2,count1)
      E(3,count1)=ima*zk*E(3,count1)-gradpot(3,count1)
    enddo

    deallocate(a_vect)
    deallocate(b_vect)
    deallocate(lambda)
    deallocate(curlE)
    deallocate(rho)
    deallocate(gradpot)
    deallocate(divE)

    deallocate(u_vect)
    deallocate(v_vect)
    deallocate(n_vect)
    deallocate(source)

  return
end subroutine em_nrccie_pec_FMM_targ
!
!
!
!
!
!
!
subroutine em_nrccie_pec_FMM_targ_oversamp(eps,zk,npatches,norders, &
    ixyzs,iptype,npts,srcvals,nt,targ,wts,novers,nptso,ixyzso,srcover, &
    whtsover,a_u,a_v,rho_in,E,H)
implicit none
!
!  This funciton computes the fields E,H at a given point using the FMM
!
!  It doesn't contain near field corrections (it's for debugging purposes)   
!
!  It uses an oversampled quadrature rule just to avoid some near evaluation
!  routines although it won't be perfect for near points even after the correction
!
!  The representation for the potentials is:
!  Representation:
!
!    H=curlS_{k}[J]
!
!    E=ikS_{k}[J]-gradS_{k}[rho]
!
!  it requires as input \rho and \rho_2:=S_{ik}[\rho]
!
!
!  input:
!    - eps: real *8
!        precision requested
!    - zk: complex *16 
!        wavenumber
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
!    - npts: integer
!        total number of discretization points on the boundary
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the 
!        discretization nodes on the surface
!          * srcvals(1:3,i) - xyz info
!          * srcvals(4:6,i) - dxyz/du info
!          * srcvals(7:9,i) - dxyz/dv info
!          * srcvals(10:12,i) - normals info
!    - wts: real *8 (npts)
!        quadrature weights for integrating smooth functions
!    - nt: integer
!        number of targets
!    - targs: real *8 (3,nt)
!        target locations
!    - novers: integer(npatches)
!        order of discretization for oversampled sources and
!        density
!    - nptso: integer
!        total number of oversampled points
!    - ixyzso: integer(npatches+1)
!        ixyzso(i) denotes the starting location in srcover,
!        corresponding to patch i
!    - srcover: real *8 (12,nptso)
!        oversampled set of source information
!    - whtsover: real *8 (nptso)
!        smooth quadrature weights at oversampled nodes
!    - a_u,a_v: complex *16 (npts)
!        two components of the tanget induced current J on the surface
!        along the srcvals(4:6,i) and srcvals(10:12,i) \times srcvals(4:6,i)
!        directions
!    - rho_in: complex *16 (npts)
!        induced charge on the surface
!
!  Output arguments:
!    - E: complex *16(3,nt)
!        value of the electric field at the target points
!    - H: complex *16(3,nt)
!        value of the magnetic field at the target points
!

! List of calling arguments
    
    real *8 eps
    complex *16 zk
    integer npatches,norders(npatches),ixyzs(npatches+1),iptype(npatches)
    integer npts
    real *8 srcvals(12,npts),wts(npts)
    integer nt
    real *8 targ(3,nt)
    integer novers(npatches),nptso,ixyzso(npatches+1)
    real *8 srcover(12,nptso),whtsover(nptso)
    complex *16 a_u(npts),a_v(npts),rho_in(npts)
    
    complex *16 E(3,nt),H(3,nt)

!    List of local variables

    real *8, allocatable :: source(:,:)
    complex *16, allocatable :: b_vect(:,:)
    complex *16, allocatable :: a_vect_over(:,:),b_vect_over(:,:),rho_over(:)
    complex *16, allocatable :: lambda_over(:)

    complex *16, allocatable :: curlE(:,:),gradpot(:,:),divE(:)

    real *8, allocatable :: u(:,:),v(:,:),rn(:,:)
    real *8, allocatable :: rnover(:,:)


    complex ( kind = 8 ) ima

    integer count1,count2,ier,i
    integer ifa_vect,ifb_vect,iflambda,ifrho,ifE,ifcurlE,ifdivE
    real *8 pi
    pi=3.1415926535897932384626433832795028841971d0

    ima=(0.0d0,1.0d0)

    allocate(u(3,npts),v(3,npts),rn(3,npts),rnover(3,nptso))
    call orthonormalize_all(srcvals(4:6,:),srcvals(10:12,:),u,v,npts)

    allocate(b_vect(3,npts))

    do i=1,npts
     b_vect(1:3,i)=a_u(i)*u(1:3,i) + a_v(i)*v(1:3,i)
    enddo

    allocate(a_vect_over(3,nptso),b_vect_over(3,nptso))
    allocate(lambda_over(nptso),rho_over(nptso))

    a_vect_over = 0
    lambda_over = 0
    rho_over = 0

    call oversample_fun_surf(6,npatches,norders,ixyzs,iptype,npts, &
      b_vect,novers,ixyzso,nptso,b_vect_over)

    allocate(source(3,nptso))

    do i=1,nptso
      rnover(1:3,i) = srcover(10:12,i)
      source(1:3,i) = srcover(1:3,i)
    enddo


    allocate(curlE(3,nt))
    allocate(gradpot(3,nt))
    allocate(divE(nt))


!  Computing the full operator
    ifa_vect=0
    ifb_vect=1
    iflambda=0
    ifrho=0
    ifE=1
    ifcurlE=1
    ifdivE=0

    call Vector_Helmholtz_targ(eps,zk,nptso,source,whtsover, &
      ifa_vect,a_vect_over,ifb_vect,b_vect_over,iflambda, &
      lambda_over,ifrho,rho_over,rnover,ifE,E,ifcurlE,H,ifdivE, &
      divE,nt,targ)
    
    call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,npts, &
      rho_in,novers,ixyzso,nptso,rho_over)

    do i=1,nptso
     rho_over(i)=rho_over(i)*whtsover(i)
    enddo


! Computing the full operator
    ier = 0 
    call hfmm3d_t_c_g(eps,zk,nptso,source,rho_over,nt,targ,divE,gradpot,ier)


    do count1=1,nt
       gradpot(:,count1)=gradpot(:,count1)/(4.0d0*pi)
       divE(count1)=divE(count1)/(4.0d0*pi)
    enddo

    do count1=1,nt
      E(1,count1)=ima*zk*E(1,count1)-gradpot(1,count1)
      E(2,count1)=ima*zk*E(2,count1)-gradpot(2,count1)
      E(3,count1)=ima*zk*E(3,count1)-gradpot(3,count1)
    enddo



  return
end subroutine em_nrccie_pec_FMM_targ_oversamp
!
!
!
!
!
subroutine 	get_rhs_em_nrccie_pec(P0,vf,alpha,ns,srcvals,zk,RHS)
implicit none
!
!  This function obtains the right hand side for the NRCCIE formulation
!  for the integral boundary equation:
!
!    J/2-M_{k}[J]+alpha·nxnx(ikS_{k}[J]-gradS_{k}[rho]) =
!      = nxH_inc + alpha nxnxE_inc
!    rho/2+S'_{k}[rho]-ikS_{k}[J]+alpha(divS_{k}[J]-ikS_{k}[rho]) = 
!      = n·E_inc
!
!  input:
!    P0 - real * 8 (3)
!      location of the source point at the exterior region
!      WARNING! notice that this formulation uses a representation theorem
!      for the incoming field in the interior region (MFIE) therefore
!      therefore it only works for incoming fields generated by sources in
!      the exterior region (or at infinity like plane waves)
!
!    vf - complex *16(3)
!      Orientation of the magnetic and electric dipoles located at P0 
!
!    alpha - complex *16
!      parameter in the combined formulation
!   
!    ns - integer
!      total number of points on the surface
!
!    srcvals - real *8(12,ns)
!      xyz(u,v) and derivative info sampled at the 
!      discretization nodes on the surface
!      srcvals(1:3,i) - xyz info
!      srcvals(4:6,i) - dxyz/du info
!      srcvals(7:9,i) - dxyz/dv info
!
!    zk - complex *16
!      Helmholtz parameter 
!
!  output:
!    RHS - complex  *16(3*ns)
!      right hand side
!      RHS(1:ns) - first component of  nxH_inc + alpha nxnxE_inc along
!       the srcvals(4:6,i) direction
!      RHS(ns+1:2*ns) - second component of  nxH_inc + alpha nxnxE_inc
!       along the (srcvals(10:12,i) x srcvals(4:6,i)) direction
!      RHS(2*ns+1:3*ns) - normal component of the electric field n·E_inc
!

	!List of calling arguments
	integer ( kind = 4 ), intent(in) :: ns
	real ( kind = 8 ), intent(in) :: P0(3)
	complex ( kind = 8 ), intent(in) :: vf(3)
	real ( kind = 8 ), intent(in) :: srcvals(12,ns)
	complex ( kind = 8 ), intent(in) :: zk,alpha
	complex ( kind = 8 ), intent(out) :: RHS(3*ns)
	
	!List of local variables
	complex ( kind = 8 ), allocatable :: E(:,:), H(:,:)
	integer count1
	real ( kind = 8 ) ru(3),rv(3),cross_aux(3)
		
	allocate(E(3,ns), H(3,ns))

	call fieldsED(zk,P0,srcvals,ns,E,H,vf,0)
	call fieldsMD(zk,P0,srcvals,ns,E,H,vf,1)
	do count1=1,ns	
      call orthonormalize(srcvals(4:6,count1),srcvals(10:12,count1),ru,rv)
      RHS(count1)=-DOT_PRODUCT(rv,H(:,count1))+&
	   &alpha*DOT_PRODUCT(ru,E(:,count1))
      RHS(ns+count1)=DOT_PRODUCT(ru,H(:,count1))+&
       &alpha*DOT_PRODUCT(rv,E(:,count1))
      RHS(2*ns+count1)=DOT_PRODUCT(srcvals(10:12,count1),E(:,count1))
	enddo

return
end subroutine get_rhs_em_nrccie_pec
