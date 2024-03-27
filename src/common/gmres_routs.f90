!  
!  This file has two user callable routines:
!    * zgmres_guru: low-threshold stagnation free gmres for 
!      complex matrices
!
!    * dgmres_guru: low-threshold stagnation free gmres for 
!      real matrices
!


      

        subroutine zgmres_guru(npatches, norders, ixyzs, &
            iptype, npts, srccoefs, srcvals, wts, &
            eps, ndd, dpars, ndz, zpars, ndi, ipars, &
            nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
            nptso, ixyzso, srcover, whtsover, lwork, work, &
            ndim, fker, zid, rhs, numit, eps_gmres, niter, errs, &
            rres, soln)
!  
!  Low-threshold stagnation free gmres for complex matrices of the form
!    (zI + K)x = y, 
!  where K is available as a matrix vector product through the 
!  subroutine fker with calling sequence
!            
!   subroutine fker(npatches, norders, ixyzs, &
!       iptype, npts, srccoefs, srcvals, &
!       eps, ndd, dpars, ndz, zpars, ndi, ipars, &
!       nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
!       nptso, ixyzso, srcover, whtsover, lwork, work, &
!       ndim, x, y)
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
!    - wts: real *8 (npts)
!        quadrature weights for integrating smooth functions     
!    - eps: real *8
!        precision requested
!    - ndd: integer
!        number of real parameters defining the kernel/
!        integral representation
!    - dpars: real *8(ndd)
!        real parameters defining the kernel/
!        integral representation
!    - ndz: integer
!        number of complex parameters defining the kernel/
!        integral representation
!    - zpars: real *8(ndz)
!        complex parameters defining the kernel/
!        integral representation.
!    - ndi: integer
!        number of integer parameters defining the kernel/
!        integral representation
!    - ipars: real *8(ndi)
!        integer parameters defining the kernel/
!        integral representation
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
!    - wnear: complex *16(nker, nquad)
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
!        number of densities
!    - fker: procedure pointer
!        function handle for evaluating the kernel K
!    - zid: complex *16
!        multiple of identity 'z' to be used
!    - rhs: complex *16(ndim, npts) or complex *16(npts)
!        data y, in case of vector densities at a point,
!        ndim > 1, then different densities at the same
!        point are continuous in memory
!    - numit: integer
!        max number of gmres iterations
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
!    - soln: complex *16(ndim, npts)
!        solution of linear system

      implicit none
      integer, intent(in) :: npatches, npts
      integer, intent(in) :: norders(npatches),ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
      real *8, intent(in) :: wts(npts)
      real *8, intent(in) :: eps
      integer, intent(in) :: ndd, ndi, ndz
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      integer, intent(in) :: ipars(ndi)
      integer, intent(in) :: nnz,nquad
      integer, intent(in) :: row_ptr(npts+1),col_ind(nnz)
      integer, intent(in) :: iquad(nnz+1)
      
      integer, intent(in) :: ndim
      
      
      integer, intent(in) :: nptso
      integer, intent(in) :: novers(npatches+1)
      integer, intent(in) :: ixyzso(npatches+1)
      real *8, intent(in) :: srcover(12,nptso), whtsover(nptso)
      
      integer, intent(in) :: nker
      complex *16 wnear(nker, nquad)

      integer, intent(in) :: lwork
      real *8, intent(in) :: work(lwork)

      
      procedure (), pointer :: fker

      complex *16, intent(in) :: zid
      complex *16, intent(in) :: rhs(ndim,npts)
      integer, intent(in) :: numit
      real *8, intent(in) :: eps_gmres

      integer, intent(out) :: niter
      real *8, intent(out) :: errs(numit+1)
      real *8, intent(out) :: rres

      complex *16, intent(out) :: soln(ndim,npts)

!
!  Temporary variables
!      
      real *8 rb, wnrm2, rmyerr
      complex *16 ztmp, ztmp2
      complex *16, allocatable :: vmat(:,:,:), hmat(:,:)
      complex *16, allocatable :: cs(:), sn(:)
      complex *16, allocatable :: svec(:), yvec(:), wtmp(:,:)
      
      integer i, j, k, l, it, it1, idim, iind 

      allocate(vmat(ndim,npts,numit+1), hmat(numit,numit))
      allocate(cs(numit), sn(numit))
      allocate(wtmp(ndim,npts), svec(numit+1), yvec(numit+1))
      
      
      niter=0
!
!      compute norm of right hand side and initialize v
! 
      rb = 0

      do i=1,numit
        cs(i) = 0
        sn(i) = 0
      enddo

!
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(idim) REDUCTION(+:rb)
      do i=1,npts
        do idim=1,ndim
           rb = rb + abs(rhs(idim,i))**2*wts(i)
        enddo
      enddo
!$OMP END PARALLEL DO      
      rb = sqrt(rb)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(idim)
      do i=1,npts
        do idim=1,ndim
          vmat(idim,i,1) = rhs(idim,i)/rb
        enddo
      enddo
!$OMP END PARALLEL DO      

      svec(1) = rb

      do it=1,numit
        it1 = it + 1

        call fker(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, &
          eps, ndd, dpars, ndz, zpars, ndi, ipars, &
          nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
          nptso, ixyzso, srcover, whtsover, lwork, work, &
          ndim, vmat(1,1,it), wtmp)

        do k=1,it
          ztmp = 0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(idim) REDUCTION(+:ztmp)          
          do j=1,npts
            do idim=1,ndim    
              ztmp = ztmp + wtmp(idim,j)*conjg(vmat(idim,j,k))*wts(j)
            enddo
          enddo
!$OMP END PARALLEL DO  
          hmat(k,it) = ztmp

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(idim)
          do j=1,npts
            do idim=1,ndim
              wtmp(idim,j) = wtmp(idim,j) - hmat(k,it)*vmat(idim,j,k)
            enddo
          enddo
!$OMP END PARALLEL DO          
        enddo
          
        hmat(it,it) = hmat(it,it)+zid
        wnrm2 = 0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(idim) REDUCTION(+:wnrm2)        
        do j=1,npts
          do idim=1,ndim
            wnrm2 = wnrm2 + abs(wtmp(idim,j))**2*wts(j)
          enddo
        enddo
!$OMP END PARALLEL DO        
        wnrm2 = sqrt(wnrm2)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(idim)
        do j=1,npts
          do idim=1,ndim
            vmat(idim,j,it1) = wtmp(idim,j)/wnrm2
          enddo
        enddo
!$OMP END PARALLEL DO        

        do k=1,it-1
          ztmp2 = cs(k)*hmat(k,it) + conjg(sn(k))*hmat(k+1,it)
          hmat(k+1,it) = -sn(k)*hmat(k,it) + cs(k)*hmat(k+1,it)
          hmat(k,it) = ztmp2
        enddo

        ztmp = wnrm2

        call zrotmat_gmres(hmat(it,it), ztmp, cs(it), sn(it))
          
        hmat(it,it) = cs(it)*hmat(it,it) + conjg(sn(it))*wnrm2
        svec(it1) = -sn(it)*svec(it)
        svec(it) = cs(it)*svec(it)
        rmyerr = abs(svec(it1))/rb
        errs(it) = rmyerr
        print *, "iter=",it,errs(it)

        if(rmyerr.le.eps_gmres.or.it.eq.numit) then

!
!            solve the linear system corresponding to
!            upper triangular part of hmat to obtain yvec
!
!            y = triu(H(1:it,1:it))\s(1:it);
!
          do j=1,it
            iind = it-j+1
            yvec(iind) = svec(iind)
            do l=iind+1,it
              yvec(iind) = yvec(iind) - hmat(iind,l)*yvec(l)
            enddo
            yvec(iind) = yvec(iind)/hmat(iind,iind)
          enddo



!
!          estimate x
!
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim) 
          do j=1,npts
            do idim=1,ndim
              soln(idim,j) = 0
              do i=1,it
                soln(idim,j) = soln(idim,j) + yvec(i)*vmat(idim,j,i)
              enddo
            enddo
          enddo
!$OMP END PARALLEL DO          


          rres = 0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(idim)        
          do i=1,npts
            do idim=1,ndim
              wtmp(idim,i) = 0
            enddo
          enddo
!$OMP END PARALLEL DO          

          call fker(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, &
          eps, ndd, dpars, ndz, zpars, ndi, ipars, &
          nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
          nptso, ixyzso, srcover, whtsover, lwork, work, &
          ndim, soln, wtmp)

!$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rres) PRIVATE(idim)
          do i=1,npts
            do idim=1,ndim
              rres = rres + abs(zid*soln(idim,i) + wtmp(idim,i) - & 
                   rhs(idim,i))**2*wts(i)
            enddo
          enddo
!$OMP END PARALLEL DO          
          rres = sqrt(rres)/rb
          niter = it
          return

        endif
      enddo

      return
      end



      subroutine dgmres_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, wts, &
        eps, ndd, dpars, ndz, zpars, ndi, ipars, &
        nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
        nptso, ixyzso, srcover, whtsover, lwork, work, &
        ndim, fker, did, rhs, numit, eps_gmres, niter, errs, &
        rres, soln)
!  
!  Low-threshold stagnation free gmres for real matrices of the form
!    (zI + K)x = y, 
!  where K is available as a matrix vector product through the 
!  subroutine fker with calling sequence
!            
!   subroutine fker(npatches, norders, ixyzs, &
!       iptype, npts, srccoefs, srcvals, &
!       eps, ndd, dpars, ndz, zpars, ndi, ipars, &
!       nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
!       nptso, ixyzso, srcover, whtsover, lwork, work, &
!       ndim, x, y)
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
!    - wts: real *8 (npts)
!        quadrature weights for integrating smooth functions                
!    - eps: real *8
!        precision requested
!    - ndd: integer
!        number of real parameters defining the kernel/
!        integral representation
!    - dpars: real *8(ndd)
!        real parameters defining the kernel/
!        integral representation
!    - ndz: integer
!        number of complex parameters defining the kernel/
!        integral representation
!    - zpars: real *8(ndz)
!        complex parameters defining the kernel/
!        integral representation.
!    - ndi: integer
!        number of integer parameters defining the kernel/
!        integral representation
!    - ipars: real *8(ndi)
!        integer parameters defining the kernel/
!        integral representation
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
!    - lwork: integer
!        size of work array (unused in this routine)
!    - work: real *8(lwork)
!        work array (unused in this routine)
!    - ndim: integer
!        number of densities
!    - fker: procedure pointer
!        function handle for evaluating the kernel K
!    - did: real *8
!        multiple of identity 'z' to be used
!    - rhs: real *8(ndim, npts) or real *8(npts)
!        data y, in case of vector densities at a point,
!        ndim > 1, then different densities at the same
!        point are continuous in memory
!    - numit: integer
!        max number of gmres iterations
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
!    - soln: real *8(ndim, npts)
!        solution of linear system

  implicit none
  integer, intent(in) :: npatches, npts
  integer, intent(in) :: norders(npatches),ixyzs(npatches+1)
  integer, intent(in) :: iptype(npatches)
  real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
  real *8, intent(in) :: wts(npts)
  real *8, intent(in) :: eps
  integer, intent(in) :: ndd, ndi, ndz
  real *8, intent(in) :: dpars(ndd)
  complex *16, intent(in) :: zpars(ndz)
  integer, intent(in) :: ipars(ndi)
  integer, intent(in) :: nnz,nquad
  integer, intent(in) :: row_ptr(npts+1),col_ind(nnz)
  integer, intent(in) :: iquad(nnz+1)
  
  integer, intent(in) :: ndim
  
  
  integer, intent(in) :: nptso
  integer, intent(in) :: novers(npatches+1)
  integer, intent(in) :: ixyzso(npatches+1)
  real *8, intent(in) :: srcover(12,nptso), whtsover(nptso)
  
  integer, intent(in) :: nker
  real *8 wnear(nker, nquad)

  integer, intent(in) :: lwork
  real *8, intent(in) :: work(lwork)

  
  procedure (), pointer :: fker

  real *8, intent(in) :: did
  real *8, intent(in) :: rhs(ndim,npts)
  integer, intent(in) :: numit
  real *8, intent(in) :: eps_gmres

  integer, intent(out) :: niter
  real *8, intent(out) :: errs(numit+1)
  real *8, intent(out) :: rres

  real *8, intent(out) :: soln(ndim,npts)

!
!  Temporary variables
!      
  real *8 rb, wnrm2, rmyerr
  real *8 ztmp, ztmp2
  real *8, allocatable :: vmat(:,:,:), hmat(:,:)
  real *8, allocatable :: cs(:), sn(:)
  real *8, allocatable :: svec(:), yvec(:), wtmp(:,:)
  
  integer i, j, k, l, it, it1, idim, iind 

  allocate(vmat(ndim,npts,numit+1), hmat(numit,numit))
  allocate(cs(numit), sn(numit))
  allocate(wtmp(ndim,npts), svec(numit+1), yvec(numit+1))
  
  
  niter=0
!
!      compute norm of right hand side and initialize v
! 
  rb = 0

  do i=1,numit
    cs(i) = 0
    sn(i) = 0
  enddo

!
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(idim) REDUCTION(+:rb)
  do i=1,npts
    do idim=1,ndim
       rb = rb + abs(rhs(idim,i))**2*wts(i)
    enddo
  enddo
!$OMP END PARALLEL DO      
  rb = sqrt(rb)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(idim)
  do i=1,npts
    do idim=1,ndim
      vmat(idim,i,1) = rhs(idim,i)/rb
    enddo
  enddo
!$OMP END PARALLEL DO      

  svec(1) = rb

  do it=1,numit
    it1 = it + 1

    call fker(npatches, norders, ixyzs, &
      iptype, npts, srccoefs, srcvals, &
      eps, ndd, dpars, ndz, zpars, ndi, ipars, &
      nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
      nptso, ixyzso, srcover, whtsover, lwork, work, &
      ndim, vmat(1,1,it), wtmp)

    do k=1,it
      ztmp = 0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(idim) REDUCTION(+:ztmp)          
      do j=1,npts
        do idim=1,ndim    
          ztmp = ztmp + wtmp(idim,j)*vmat(idim,j,k)*wts(j)
        enddo
      enddo
!$OMP END PARALLEL DO  
      hmat(k,it) = ztmp

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(idim)
      do j=1,npts
        do idim=1,ndim
          wtmp(idim,j) = wtmp(idim,j) - hmat(k,it)*vmat(idim,j,k)
        enddo
      enddo
!$OMP END PARALLEL DO          
    enddo
      
    hmat(it,it) = hmat(it,it)+did
    wnrm2 = 0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(idim) REDUCTION(+:wnrm2)        
    do j=1,npts
      do idim=1,ndim
        wnrm2 = wnrm2 + abs(wtmp(idim,j))**2*wts(j)
      enddo
    enddo
!$OMP END PARALLEL DO        
    wnrm2 = sqrt(wnrm2)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(idim)
    do j=1,npts
      do idim=1,ndim
        vmat(idim,j,it1) = wtmp(idim,j)/wnrm2
      enddo
    enddo
!$OMP END PARALLEL DO        

    do k=1,it-1
      ztmp2 = cs(k)*hmat(k,it) + sn(k)*hmat(k+1,it)
      hmat(k+1,it) = -sn(k)*hmat(k,it) + cs(k)*hmat(k+1,it)
      hmat(k,it) = ztmp2
    enddo

    ztmp = wnrm2

    call rotmat_gmres(hmat(it,it), ztmp, cs(it), sn(it))
      
    hmat(it,it) = cs(it)*hmat(it,it) + sn(it)*wnrm2
    svec(it1) = -sn(it)*svec(it)
    svec(it) = cs(it)*svec(it)
    rmyerr = abs(svec(it1))/rb
    errs(it) = rmyerr
    print *, "iter=",it,errs(it)

    if(rmyerr.le.eps_gmres.or.it.eq.numit) then

!
!            solve the linear system corresponding to
!            upper triangular part of hmat to obtain yvec
!
!            y = triu(H(1:it,1:it))\s(1:it);
!
      do j=1,it
        iind = it-j+1
        yvec(iind) = svec(iind)
        do l=iind+1,it
          yvec(iind) = yvec(iind) - hmat(iind,l)*yvec(l)
        enddo
        yvec(iind) = yvec(iind)/hmat(iind,iind)
      enddo



!
!          estimate x
!
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim) 
      do j=1,npts
        do idim=1,ndim
          soln(idim,j) = 0
          do i=1,it
            soln(idim,j) = soln(idim,j) + yvec(i)*vmat(idim,j,i)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO          


      rres = 0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(idim)        
      do i=1,npts
        do idim=1,ndim
          wtmp(idim,i) = 0
        enddo
      enddo
!$OMP END PARALLEL DO          

      call fker(npatches, norders, ixyzs, &
      iptype, npts, srccoefs, srcvals, &
      eps, ndd, dpars, ndz, zpars, ndi, ipars, &
      nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, &
      nptso, ixyzso, srcover, whtsover, lwork, work, &
      ndim, soln, wtmp)

!$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rres) PRIVATE(idim)
      do i=1,npts
        do idim=1,ndim
          rres = rres + abs(did*soln(idim,i) + wtmp(idim,i) - & 
               rhs(idim,i))**2*wts(i)
        enddo
      enddo
!$OMP END PARALLEL DO          
      rres = sqrt(rres)/rb
      niter = it
      return

    endif
  enddo

  return
  end