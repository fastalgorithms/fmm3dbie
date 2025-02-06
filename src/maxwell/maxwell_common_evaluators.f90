      subroutine em_debye_eval_addsub(npatches, norders, &
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
!  E = \alpha_{E} S_{k} [J] + \beta_{E} \nabla S_{k}[\rho] + 
!      \gamma_{E} \nabla \times S_{k} [M]
!  
!  H = \alpha_{H} S_{k} [M] + \beta_{H} \nabla S_{k}[\tau] + 
!      \gamma_{H} \nabla \times S_{k} [J] 
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
!        integral representation (must be seven)
!    - zpars: complex *16(ndz)
!        complex parameters defining the kernel/
!        integral represnetation.
!        * zpars(1) = k 
!        * zpars(2) = \alpha_{E}
!        * zpars(3) = \beta_{E}
!        * zpars(4) = \gamma_{E}
!        * zpars(5) = \alpha_{H}
!        * zpars(6) = \beta_{H}
!        * zpars(7) = \gamma_{H}
!        Parts of this representation maybe unused 
!        depending on the ipotflag, see below        
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
!        number of kernels in quadrature correction, must be 1 or 4
!        nker is 1 is only allowed if E = S_{k}[J] or H = S_{k}[M]
!        is requested
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
!        Flag for types of denisties
!        idensflag = ifj + 2^(ifm)*ifm + 4^(ifrho)*ifrho +
!                      8^(iftau)*iftau
!        where ifj is 1 if J current is included and 0 otherwise
!              ifm is 1 if M current is included and 0 otherwise
!              ifrho is 1 if \rho charge is included and 0 otherwise
!              iftau is 1 if \tau charge is included and 0 otherwise
!        i.e. the digits of the binary expansion of idensflag
!             are encoded by (iftau, ifrho, ifm, ifj)                
!    - ndim_s: integer
!        number of densities per point on the surface,
!        must be 8 for this routine
!    - sigma: complex *16(ndim_s, npts)
!        The densities
!        * sigma(1,:) stores the x_{1} component of J
!        * sigma(2,:) stores the x_{2} component of J
!        * sigma(3,:) stores the x_{3} component of J
!        * sigma(4,:) stores the x_{1} component of M
!        * sigma(5,:) stores the x_{2} component of M
!        * sigma(6,:) stores the x_{3} component of M
!        * sigma(7,:) stores \rho
!        * sigma(8,:) stores \tau            
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

      complex *16, intent(in) :: sigma(ndim_s, npts)

      complex *16, intent(out) :: pot(ndim_p, ntarg)


      real *8, allocatable :: sources(:,:), targtmp(:,:)
      complex *16, allocatable :: h_current(:,:,:), e_current(:,:,:)
      complex *16, allocatable :: e_charge(:,:)
      complex *16, allocatable :: sigmaover(:,:)
      complex *16, allocatable :: potloc(:,:,:), divpotloc(:,:)
      complex *16, allocatable :: curlpotloc(:,:,:)

      integer ifdens(4), itmp

      integer npols
      integer ns
      integer ifpgh, ifpghtarg


      integer i, j, jpatch, jquadstart, jstart
      complex *16, allocatable :: pottmp(:,:), curlpottmp(:,:)
      complex *16, allocatable :: divpottmp(:)

!$    real *8 omp_get_wtime      

      real *8, allocatable :: srctmp2(:,:)
      complex *16, allocatable :: h_current2(:,:,:), e_current2(:,:,:)
      complex *16, allocatable :: e_charge2(:,:)
      real *8 thresh

      real *8 over4pi
      integer nss, l, ier
      complex *16 ima, ztmp, ztmp2, ztmp3

      integer nd, ntarg0, nmax, idim

      integer i1, i2
      integer ifj, ifm, ifrho, iftau
      integer ife_current, ife_charge, ifh_current
      integer ife, ifdive, ifcurle

      real *8 done, pi, rr, rtmp
      integer nddtmp, nditmp, ndztmp
      complex *16 zpars2(2)
      data ima/(0.0d0,1.0d0)/
      data over4pi/0.07957747154594767d0/

      parameter (ntarg0=1)

      ns = nptso

      ifpgh = 0
      ifpghtarg = 1

      i1 = 0
      i2 = 0
      if(ipotflag.eq.1.or.ipotflag.eq.3) i1 = 1
      if(ipotflag.eq.2.or.ipotflag.eq.3) i2 = 1

      nd = i1 + i2
      


      allocate(sources(3,ns), targtmp(3,ntarg))
      allocate(h_current(nd,3,ns), e_current(nd,3,ns))
      allocate(e_charge(nd,ns))
      allocate(sigmaover(ndim_s,ns))
      allocate(pottmp(nd,3), curlpottmp(nd,3), divpottmp(nd))
!
!  estimate max number of sources in the near field of any target
!
!    
      call get_near_corr_max(ntarg, row_ptr, nnz, col_ind, npatches, &
        ixyzso, nmax)
      allocate(srctmp2(3,nmax), h_current2(nd,3,nmax))
      allocate(e_current2(nd,3,nmax), e_charge2(nd,nmax))
! 
!       oversample density
!
      call oversample_fun_surf(2*ndim_s, npatches, norders, ixyzs, &
        iptype, npts, sigma, novers, ixyzso, ns, sigmaover)

      itmp = idensflag
      do i=1,4
        ifdens(i) = mod(itmp,2)
        itmp = itmp/2
      enddo
      ifj = ifdens(1)
      ifm = ifdens(2)
      ifrho = ifdens(3)
      iftau = ifdens(4)

      ife_current = 0
      ifh_current = 0
      ife_charge = 0
!                  
!

!
!  extract source and target info
!
      
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(rtmp, ztmp, ztmp2, ztmp3, idim)
      do i=1,ns
        idim = 1
        rtmp = whtsover(i)*over4pi
        ztmp = rtmp*zpars(2)
        ztmp2 = rtmp*zpars(3)
        ztmp3 = rtmp*zpars(4)
        

        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)

        e_current(1:nd,1:3,i) = 0
        h_current(1:nd,1:3,i) = 0
        e_charge(1:nd,i) = 0

        if(i1.eq.1) then
          if(ifj.eq.1) then
            e_current(1,1:3,i) = sigmaover(1:3,i)*ztmp
          endif
          if(ifm.eq.1) then
            h_current(1,1:3,i) = sigmaover(4:6,i)*ztmp3
          endif
          if(ifrho.eq.1) then
            e_charge(1,i) = sigmaover(7,i)*ztmp2
          endif
          idim = idim + 1
        endif
        
        ztmp = rtmp*zpars(5)
        ztmp2 = rtmp*zpars(6)
        ztmp3 = rtmp*zpars(7)
        
        if(i2.eq.1) then
          if(ifj.eq.1) then
            h_current(idim,1:3,i) = sigmaover(1:3,i)*ztmp3
          endif
          if(ifm.eq.1) then
            e_current(idim,1:3,i) = sigmaover(4:6,i)*ztmp
          endif
          if(iftau.eq.1) then
            e_charge(idim,i) = sigmaover(8,i)*ztmp2
          endif
        endif
      enddo
!$OMP END PARALLEL DO      

      if(ifj.eq.1.and.i1.eq.1) ife_current = 1
      if(ifm.eq.1.and.i1.eq.1) ifh_current = 1
      if(ifrho.eq.1.and.i1.eq.1) ife_charge = 1

      if(ifj.eq.1.and.i2.eq.1) ifh_current = 1
      if(ifm.eq.1.and.i2.eq.1) ife_current = 1
      if(iftau.eq.1.and.i2.eq.1) ife_charge = 1      
      allocate(potloc(nd,3,ntarg), divpotloc(nd,ntarg))
      allocate(curlpotloc(nd,3,ntarg))

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        targtmp(1,i) = targs(1,i)
        targtmp(2,i) = targs(2,i)
        targtmp(3,i) = targs(3,i)
        pot(1:6,i) = 0
        potloc(1:nd,1:3,i) = 0
        divpotloc(1:nd,i) = 0
        curlpotloc(1:nd,1:3,i) = 0
      enddo

!$OMP END PARALLEL DO

! 
!  Compute the maxwell fmms
      ier = 0
      ife = 1
      ifdive = 0
      ifcurle = 0

      call emfmm3d(nd, eps, zpars(1), ns, sources, ifh_current, &
        h_current, ife_current, e_current, ife_charge, e_charge, &
        ntarg, targtmp, ife, potloc, ifcurle, curlpotloc, ifdive, &
        divpotloc, ier)    
      
      idim = 1
      if(i1.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED)        
        do i=1,ntarg
          pot(1:3,i) = potloc(1,1:3,i)
        enddo
!$OMP END PARALLEL DO        
        idim = idim + 1
      endif

      if(i2.eq.1) then
!$OMP PARALLEL DO DEFAULT(SHARED)        
        do i=1,ntarg
          pot(4:6,i) = potloc(idim,1:3,i)
        enddo
!$OMP END PARALLEL DO
      endif
!        compute threshold for ignoring local computation
      call get_fmm_thresh(12, ns, srcover, ndtarg, ntarg, targs, thresh)

!
!       Add near field precomputed contribution
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
            if(i1.eq.1.and.ifj.eq.1) then
              pot(1:3,i) = pot(1:3,i) + &
                wnear(1,jquadstart+l-1)*sigma(1:3,jstart+l-1)*zpars(2)
            endif
            if(i1.eq.1.and.ifm.eq.1) then
              pot(1,i) = pot(1,i) + &
                (wnear(3,jquadstart+l-1)*sigma(6,jstart+l-1) - &
                  wnear(4,jquadstart+l-1)*sigma(5,jstart+l-1))*zpars(4)
              pot(2,i) = pot(2,i) + &
                (wnear(4,jquadstart+l-1)*sigma(4,jstart+l-1) - &
                  wnear(2,jquadstart+l-1)*sigma(6,jstart+l-1))*zpars(4)
              pot(3,i) = pot(3,i) + &
                (wnear(2,jquadstart+l-1)*sigma(5,jstart+l-1) - &
                  wnear(3,jquadstart+l-1)*sigma(4,jstart+l-1))*zpars(4)
            endif
            if(i1.eq.1.and.ifrho.eq.1) then
              pot(1:3,i) = pot(1:3,i) + &
                wnear(2:4,jquadstart+l-1)*sigma(7,jstart+l-1)*zpars(3)
            endif

            if(i2.eq.1.and.ifm.eq.1) then
              pot(4:6,i) = pot(4:6,i) + &
                wnear(1,jquadstart+l-1)*sigma(4:6,jstart+l-1)*zpars(5)
            endif
            if(i2.eq.1.and.ifj.eq.1) then
              pot(4,i) = pot(4,i) + &
                (wnear(3,jquadstart+l-1)*sigma(3,jstart+l-1) - &
                  wnear(4,jquadstart+l-1)*sigma(2,jstart+l-1))*zpars(7)
              pot(5,i) = pot(5,i) + &
                (wnear(4,jquadstart+l-1)*sigma(1,jstart+l-1) - &
                  wnear(2,jquadstart+l-1)*sigma(3,jstart+l-1))*zpars(7)
              pot(6,i) = pot(6,i) + &
                (wnear(2,jquadstart+l-1)*sigma(2,jstart+l-1) - &
                  wnear(3,jquadstart+l-1)*sigma(3,jstart+l-1))*zpars(7)
            endif
            if(i2.eq.1.and.iftau.eq.1) then
              pot(4:6,i) = pot(4:6,i) + &
                wnear(2:4,jquadstart+l-1)*sigma(8,jstart+l-1)*zpars(6)
            endif

          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

    
!
!     Remove near contribution of the FMM
! 

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, jpatch, srctmp2) &
!$OMP PRIVATE(e_current2, h_current2, e_charge2, l, jstart, nss) &
!$OMP PRIVATE(pottmp, curlpottmp, divpottmp, idim)
      do i=1,ntarg
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
            nss = nss+1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)
            e_current2(1:nd,1:3,nss) = e_current(1:nd,1:3,l)
            h_current2(1:nd,1:3,nss) = h_current(1:nd,1:3,l)
            e_charge2(1:nd,nss) = e_charge(1:nd,l)
          enddo
        enddo

        pottmp = 0

        call em3ddirect(nd, zpars(1), nss, srctmp2, ifh_current, &
          h_current2, ife_current, e_current2, ife_charge, e_charge2, &
          ntarg0, targtmp(1,i), ife, pottmp, ifcurle, curlpottmp, &
          ifdive, divpottmp, thresh)
        idim = 1
        if(i1.eq.1) then
          pot(1:3,i) = pot(1:3,i) - pottmp(1,1:3)
          idim = idim + 1
        endif
        if(i2.eq.1) then
          pot(4:6,i) = pot(4:6,i) - pottmp(idim,1:3)
        endif
      enddo
!$OMP END PARALLEL DO      



return
end
!
!
!
!
