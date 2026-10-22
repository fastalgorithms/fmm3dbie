c
c
c   This file contains subroutines for computing the subtraction
c   step corresponding to the computed quadrature correction
c
c   The subtraction is done in place on the near field quadrature,
c   so it improves the CPU time performance of the solver without
c   increasing its memory requirements
c
c
c
c
c     This file has the following user callable routines
c
c        ?getnearquadsub - subtract, in place, the oversampled
c         smooth rule from a near field quadrature, interpolated back
c         to the discretization nodes
c
c        ?getnearquadsub_guru - guru interface for ?getnearquadsub,
c         with user provided oversampled geometry and fmm threshold
c        
c
c         z - complex
c         d - double precision
c
c
c
c
c
c

      subroutine zgetnearquadsub_guru(npatches,norders,ixyzs,
     1   iptype,novers,ixyzso,npts_over,srcover,wover,ndtarg,ntarg,
     2   targvals,thresh,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,
     3   row_ptr,col_ind,iquad,nquad,nker,iker,wnear)
c
c------------------------
c  This subroutine subtracts, in place, the oversampled smooth
c  quadrature rule from a near field quadrature wnear (as computed
c  by zgetnearquad_ggq_guru), for a given kernel, where the near
c  field is specified by the user in row sparse compressed format.
c
c  The smooth rule is mapped back to the discretization nodes on
c  each patch, so that wnear and the result share the same layout
c  (iquad). That is, for target i and source patch j in its near
c  field, the matrix entries
c
c    wnear(iker,iquad+m-1) <- wnear(iker,iquad+m-1) -
c       sum_{l} K(x_i,y_l) wover(l) P_{lm}
c
c  where y_l are the oversampled nodes on patch j, and P is the
c  interpolation from the discretization nodes on patch j to the
c  oversampled nodes. Oversampled nodes within thresh of the target
c  are skipped.
c
c  The interpolation matrices are computed once per unique
c  (norder,nover,iptype).
c
c  Input arguments:
c    - npatches: integer *8
c        number of patches
c    - norders: integer *8(npatches)
c        order of discretization on each patch
c    - ixyzs: integer *8(npatches+1)
c        ixyzs(i) denotes the starting location of patch i in the
c        discretization nodes
c    - iptype: integer *8(npatches)
c        type of patch
c        *  iptype = 1, triangular patch discretized using RV nodes
c        *  iptype = 11, quad patch with GL nodes
c        *  iptype = 12, quad patch with Chebyshev nodes
c    - novers: integer *8(npatches)
c        oversampled order on each patch
c    - ixyzso: integer *8(npatches+1)
c        ixyzso(i) denotes the starting location in srcover,
c        corresponding to patch i
c    - npts_over: integer *8
c        total number of oversampled points on the boundary
c    - srcover: double precision (12,npts_over)
c        xyz(u,v) and derivative info sampled at the
c        oversampled nodes on the surface (see oversample_geom)
c    - wover: double precision (npts_over)
c        quadrature weights at the oversampled nodes (see get_qwts)
c    - ndtarg: integer *8
c        leading dimension of target array
c    - ntarg: integer *8
c        number of targets
c    - targvals: double precision (ndtarg,ntarg)
c        target info. First three components must be x,y,z
c        coordinates
c    - thresh: double precision
c        threshold for ignoring close interactions
c        (see get_fmm_thresh)
c    - fker: procedure pointer
c        function handle for the kernel. Calling sequence
c        * fker(srcinfo,ndtarg,targinfo,ndd,dpars,ndz,zpars,ndi,ipars,val)
c    - ndd: integer *8
c        number of double precision parameters
c    - dpars: double precision(ndd)
c        double precision parameters
c    - ndz: integer *8
c        number of double complex parameters
c    - zpars: double complex(ndz)
c        double complex parameters
c    - ndi: integer *8
c        number of integer *8 parameters
c    - ipars: integer *8(ndi)
c        integer *8 parameters
c    - nnz: integer *8
c        number of source patch-> target interactions in the near field
c    - row_ptr: integer *8(ntarg+1)
c        row_ptr(i) is the pointer
c        to col_ind array where list of relevant source patches
c        for target i start
c    - col_ind: integer *8 (nnz)
c        list of source patches relevant for all targets, sorted
c        by the target number
c    - iquad: integer *8(nnz+1)
c        location in wnear array where quadrature for col_ind(i)
c        starts (on the discretization nodes, not the oversampled
c        ones)
c    - nquad: integer *8
c        number of near field entries corresponding to each source
c        target pair
c    - nker: integer *8
c        number of kernels stored in wnear
c    - iker: integer *8
c        kernel to update, 1 <= iker <= nker
c
c  Input/output parameters:
c    - wnear: double complex(nker,nquad)
c        on input, wnear(iker,:) is the near field quadrature; on
c        output, it is the near field quadrature minus the
c        oversampled smooth rule. Other kernels are not accessed
c----------------------------------------------------
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: npatches,norders(npatches)
      integer *8, intent(in) :: ixyzs(npatches+1),iptype(npatches)
      integer *8, intent(in) :: novers(npatches)
      integer *8, intent(in) :: ixyzso(npatches+1),npts_over
      real *8, intent(in) :: srcover(12,npts_over),wover(npts_over)
      integer *8, intent(in) :: ndtarg,ntarg
      real *8, intent(in) :: targvals(ndtarg,ntarg),thresh
      integer *8, intent(in) :: ndd,ndz,ndi
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      integer *8, intent(in) :: ipars(ndi)
      integer *8, intent(in) :: nnz
      integer *8, intent(in) :: row_ptr(ntarg+1),col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1),nquad
      integer *8, intent(in) :: nker,iker
      complex *16, intent(inout) :: wnear(nker,nquad)

      integer *8, allocatable :: iuni(:,:),iuse(:),ixmat(:)
      real *8, allocatable :: xmats(:),uvs(:,:),uvso(:,:)
      real *8, allocatable :: umat(:,:),vmat(:,:),wts(:),pols(:)
      complex *16, allocatable :: vals(:)
      integer *8 nuni,k,norder,nover,ipt,npols,npolso,npmaxo,lxmat
      integer *8 i,j,l,m,ipatch,istarto,iqstart,ix
      real *8 rr,threshsq
      complex *16 val

      external fker

      threshsq = thresh**2

c
c  unique (norder,nover,iptype) and their interpolation matrices
c  xmat(m,l) = sum_k umat(k,m) p_k(uvso(l)), stored as (npols,npolso)
c
      allocate(iuni(3,npatches),iuse(npatches))
      call get_iuni3(npatches,norders,novers,iptype,nuni,iuni,iuse)

      allocate(ixmat(nuni+1))
      ixmat(1) = 1
      npmaxo = 0
      do k=1,nuni
        call get_npols(iuni(3,k),iuni(1,k),npols)
        call get_npols(iuni(3,k),iuni(2,k),npolso)
        ixmat(k+1) = ixmat(k) + npols*npolso
        if(npolso.gt.npmaxo) npmaxo = npolso
      enddo
      lxmat = ixmat(nuni+1)-1
      allocate(xmats(lxmat))

      do k=1,nuni
        norder = iuni(1,k)
        nover = iuni(2,k)
        ipt = iuni(3,k)
        call get_npols(ipt,norder,npols)
        call get_npols(ipt,nover,npolso)
        allocate(uvs(2,npols),umat(npols,npols),vmat(npols,npols))
        allocate(wts(npols),uvso(2,npolso),pols(npols))
        call get_disc_exps(norder,npols,ipt,uvs,umat,vmat,wts)
        call get_disc_nodes(nover,npolso,ipt,uvso)
        do l=1,npolso
          call get_basis_pols(uvso(1,l),norder,npols,ipt,pols)
          ix = ixmat(k) + (l-1)*npols - 1
          do m=1,npols
            xmats(ix+m) = 0
            do i=1,npols
              xmats(ix+m) = xmats(ix+m) + umat(i,m)*pols(i)
            enddo
          enddo
        enddo
        deallocate(uvs,umat,vmat,wts,uvso,pols)
      enddo

C$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,ipatch,istarto,npolso)
C$OMP$PRIVATE(npols,iqstart,ix,l,m,rr,val,vals)
      allocate(vals(npmaxo))
C$OMP DO SCHEDULE(DYNAMIC)
      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          ipatch = col_ind(j)
          istarto = ixyzso(ipatch)
          npolso = ixyzso(ipatch+1)-ixyzso(ipatch)
          npols = ixyzs(ipatch+1)-ixyzs(ipatch)
          iqstart = iquad(j)-1
          ix = ixmat(iuse(ipatch))-1
          do l=1,npolso
            vals(l) = 0
            rr = (srcover(1,l+istarto-1) - targvals(1,i))**2 +
     1          (srcover(2,l+istarto-1) - targvals(2,i))**2 +
     2          (srcover(3,l+istarto-1) - targvals(3,i))**2
            if(rr.gt.threshsq) then
              call fker(srcover(1,l+istarto-1),ndtarg,targvals(1,i),
     1          ndd,dpars,ndz,zpars,ndi,ipars,val)
              vals(l) = val*wover(l+istarto-1)
            endif
          enddo
          do l=1,npolso
            do m=1,npols
              wnear(iker,iqstart+m) = wnear(iker,iqstart+m) -
     1           xmats(ix+(l-1)*npols+m)*vals(l)
            enddo
          enddo
        enddo
      enddo
C$OMP END DO
      deallocate(vals)
C$OMP END PARALLEL

      return
      end

c
c
c
c
      subroutine dgetnearquadsub_guru(npatches,norders,ixyzs,
     1   iptype,novers,ixyzso,npts_over,srcover,wover,ndtarg,ntarg,
     2   targvals,thresh,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,
     3   row_ptr,col_ind,iquad,nquad,nker,iker,wnear)
c
c------------------------
c  This subroutine subtracts, in place, the oversampled smooth
c  quadrature rule from a near field quadrature wnear (as computed
c  by dgetnearquad_ggq_guru), for a given kernel, where the near
c  field is specified by the user in row sparse compressed format.
c
c  The smooth rule is mapped back to the discretization nodes on
c  each patch, so that wnear and the result share the same layout
c  (iquad). That is, for target i and source patch j in its near
c  field, the matrix entries
c
c    wnear(iker,iquad+m-1) <- wnear(iker,iquad+m-1) -
c       sum_{l} K(x_i,y_l) wover(l) P_{lm}
c
c  where y_l are the oversampled nodes on patch j, and P is the
c  interpolation from the discretization nodes on patch j to the
c  oversampled nodes. Oversampled nodes within thresh of the target
c  are skipped.
c
c  The interpolation matrices are computed once per unique
c  (norder,nover,iptype).
c
c  Input arguments:
c    - npatches: integer *8
c        number of patches
c    - norders: integer *8(npatches)
c        order of discretization on each patch
c    - ixyzs: integer *8(npatches+1)
c        ixyzs(i) denotes the starting location of patch i in the
c        discretization nodes
c    - iptype: integer *8(npatches)
c        type of patch
c        *  iptype = 1, triangular patch discretized using RV nodes
c        *  iptype = 11, quad patch with GL nodes
c        *  iptype = 12, quad patch with Chebyshev nodes
c    - novers: integer *8(npatches)
c        oversampled order on each patch
c    - ixyzso: integer *8(npatches+1)
c        ixyzso(i) denotes the starting location in srcover,
c        corresponding to patch i
c    - npts_over: integer *8
c        total number of oversampled points on the boundary
c    - srcover: double precision (12,npts_over)
c        xyz(u,v) and derivative info sampled at the
c        oversampled nodes on the surface (see oversample_geom)
c    - wover: double precision (npts_over)
c        quadrature weights at the oversampled nodes (see get_qwts)
c    - ndtarg: integer *8
c        leading dimension of target array
c    - ntarg: integer *8
c        number of targets
c    - targvals: double precision (ndtarg,ntarg)
c        target info. First three components must be x,y,z
c        coordinates
c    - thresh: double precision
c        threshold for ignoring close interactions
c        (see get_fmm_thresh)
c    - fker: procedure pointer
c        function handle for the kernel. Calling sequence
c        * fker(srcinfo,ndtarg,targinfo,ndd,dpars,ndz,zpars,ndi,ipars,val)
c    - ndd: integer *8
c        number of double precision parameters
c    - dpars: double precision(ndd)
c        double precision parameters
c    - ndz: integer *8
c        number of double complex parameters
c    - zpars: double complex(ndz)
c        double complex parameters
c    - ndi: integer *8
c        number of integer *8 parameters
c    - ipars: integer *8(ndi)
c        integer *8 parameters
c    - nnz: integer *8
c        number of source patch-> target interactions in the near field
c    - row_ptr: integer *8(ntarg+1)
c        row_ptr(i) is the pointer
c        to col_ind array where list of relevant source patches
c        for target i start
c    - col_ind: integer *8 (nnz)
c        list of source patches relevant for all targets, sorted
c        by the target number
c    - iquad: integer *8(nnz+1)
c        location in wnear array where quadrature for col_ind(i)
c        starts (on the discretization nodes, not the oversampled
c        ones)
c    - nquad: integer *8
c        number of near field entries corresponding to each source
c        target pair
c    - nker: integer *8
c        number of kernels stored in wnear
c    - iker: integer *8
c        kernel to update, 1 <= iker <= nker
c
c  Input/output parameters:
c    - wnear: double precision(nker,nquad)
c        on input, wnear(iker,:) is the near field quadrature; on
c        output, it is the near field quadrature minus the
c        oversampled smooth rule. Other kernels are not accessed
c----------------------------------------------------
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: npatches,norders(npatches)
      integer *8, intent(in) :: ixyzs(npatches+1),iptype(npatches)
      integer *8, intent(in) :: novers(npatches)
      integer *8, intent(in) :: ixyzso(npatches+1),npts_over
      real *8, intent(in) :: srcover(12,npts_over),wover(npts_over)
      integer *8, intent(in) :: ndtarg,ntarg
      real *8, intent(in) :: targvals(ndtarg,ntarg),thresh
      integer *8, intent(in) :: ndd,ndz,ndi
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      integer *8, intent(in) :: ipars(ndi)
      integer *8, intent(in) :: nnz
      integer *8, intent(in) :: row_ptr(ntarg+1),col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1),nquad
      integer *8, intent(in) :: nker,iker
      real *8, intent(inout) :: wnear(nker,nquad)

      integer *8, allocatable :: iuni(:,:),iuse(:),ixmat(:)
      real *8, allocatable :: xmats(:),uvs(:,:),uvso(:,:)
      real *8, allocatable :: umat(:,:),vmat(:,:),wts(:),pols(:)
      real *8, allocatable :: vals(:)
      integer *8 nuni,k,norder,nover,ipt,npols,npolso,npmaxo,lxmat
      integer *8 i,j,l,m,ipatch,istarto,iqstart,ix
      real *8 rr,threshsq
      real *8 val

      external fker

      threshsq = thresh**2

c
c  unique (norder,nover,iptype) and their interpolation matrices
c  xmat(m,l) = sum_k umat(k,m) p_k(uvso(l)), stored as (npols,npolso)
c
      allocate(iuni(3,npatches),iuse(npatches))
      call get_iuni3(npatches,norders,novers,iptype,nuni,iuni,iuse)

      allocate(ixmat(nuni+1))
      ixmat(1) = 1
      npmaxo = 0
      do k=1,nuni
        call get_npols(iuni(3,k),iuni(1,k),npols)
        call get_npols(iuni(3,k),iuni(2,k),npolso)
        ixmat(k+1) = ixmat(k) + npols*npolso
        if(npolso.gt.npmaxo) npmaxo = npolso
      enddo
      lxmat = ixmat(nuni+1)-1
      allocate(xmats(lxmat))

      do k=1,nuni
        norder = iuni(1,k)
        nover = iuni(2,k)
        ipt = iuni(3,k)
        call get_npols(ipt,norder,npols)
        call get_npols(ipt,nover,npolso)
        allocate(uvs(2,npols),umat(npols,npols),vmat(npols,npols))
        allocate(wts(npols),uvso(2,npolso),pols(npols))
        call get_disc_exps(norder,npols,ipt,uvs,umat,vmat,wts)
        call get_disc_nodes(nover,npolso,ipt,uvso)
        do l=1,npolso
          call get_basis_pols(uvso(1,l),norder,npols,ipt,pols)
          ix = ixmat(k) + (l-1)*npols - 1
          do m=1,npols
            xmats(ix+m) = 0
            do i=1,npols
              xmats(ix+m) = xmats(ix+m) + umat(i,m)*pols(i)
            enddo
          enddo
        enddo
        deallocate(uvs,umat,vmat,wts,uvso,pols)
      enddo

C$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,ipatch,istarto,npolso)
C$OMP$PRIVATE(npols,iqstart,ix,l,m,rr,val,vals)
      allocate(vals(npmaxo))
C$OMP DO SCHEDULE(DYNAMIC)
      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          ipatch = col_ind(j)
          istarto = ixyzso(ipatch)
          npolso = ixyzso(ipatch+1)-ixyzso(ipatch)
          npols = ixyzs(ipatch+1)-ixyzs(ipatch)
          iqstart = iquad(j)-1
          ix = ixmat(iuse(ipatch))-1
          do l=1,npolso
            vals(l) = 0
            rr = (srcover(1,l+istarto-1) - targvals(1,i))**2 +
     1          (srcover(2,l+istarto-1) - targvals(2,i))**2 +
     2          (srcover(3,l+istarto-1) - targvals(3,i))**2
            if(rr.gt.threshsq) then
              call fker(srcover(1,l+istarto-1),ndtarg,targvals(1,i),
     1          ndd,dpars,ndz,zpars,ndi,ipars,val)
              vals(l) = val*wover(l+istarto-1)
            endif
          enddo
          do l=1,npolso
            do m=1,npols
              wnear(iker,iqstart+m) = wnear(iker,iqstart+m) -
     1           xmats(ix+(l-1)*npols+m)*vals(l)
            enddo
          enddo
        enddo
      enddo
C$OMP END DO
      deallocate(vals)
C$OMP END PARALLEL

      return
      end
c
c
c
c
      subroutine zgetnearquadsub(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,ndtarg,ntarg,targvals,novers,fker,ndd,dpars,
     2   ndz,zpars,ndi,ipars,nnz,row_ptr,col_ind,iquad,nquad,nker,iker,
     3   wnear)
c
c------------------------
c  This subroutine subtracts, in place, the smooth quadrature rule
c  oversampled to novers from a near field quadrature wnear (as
c  computed by zgetnearquad_ggq_guru), for a given kernel, where
c  the near field is specified by the user in row sparse compressed
c  format.
c
c  The oversampled geometry, smooth quadrature weights, and fmm
c  threshold are computed here; see zgetnearquadsub_guru for
c  details.
c
c  Input arguments:
c    - npatches: integer *8
c        number of patches
c    - norders: integer *8(npatches)
c        order of discretization on each patch
c    - ixyzs: integer *8(npatches+1)
c        ixyzs(i) denotes the starting location in srccoefs,
c        and srcvals array corresponding to patch i
c    - iptype: integer *8(npatches)
c        type of patch
c        *  iptype = 1, triangular patch discretized using RV nodes
c        *  iptype = 11, quad patch with GL nodes
c        *  iptype = 12, quad patch with Chebyshev nodes
c    - npts: integer *8
c        total number of discretization points on the boundary
c    - srccoefs: double precision (9,npts)
c        basis expansion coefficients of xyz, dxyz/du, and dxyz/dv
c        on each patch
c    - srcvals: double precision (12,npts)
c        xyz(u,v) and derivative info sampled at the discretization
c        nodes on the surface
c    - ndtarg: integer *8
c        leading dimension of target array
c    - ntarg: integer *8
c        number of targets
c    - targvals: double precision (ndtarg,ntarg)
c        target info. First three components must be x,y,z
c        coordinates
c    - novers: integer *8(npatches)
c        oversampled order on each patch, must match the one used
c        by the fmm in the evaluator
c    - fker, ndd, dpars, ndz, zpars, ndi, ipars, nnz, row_ptr,
c      col_ind, iquad, nquad, nker, iker:
c        as in zgetnearquadsub_guru
c
c  Input/output parameters:
c    - wnear: double complex(nker,nquad)
c        on input, wnear(iker,:) is the near field quadrature; on
c        output, it is the near field quadrature minus the
c        oversampled smooth rule. Other kernels are not accessed
c----------------------------------------------------
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: npatches,norders(npatches)
      integer *8, intent(in) :: ixyzs(npatches+1),iptype(npatches)
      integer *8, intent(in) :: npts
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
      integer *8, intent(in) :: ndtarg,ntarg
      real *8, intent(in) :: targvals(ndtarg,ntarg)
      integer *8, intent(in) :: novers(npatches)
      integer *8, intent(in) :: ndd,ndz,ndi
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      integer *8, intent(in) :: ipars(ndi)
      integer *8, intent(in) :: nnz
      integer *8, intent(in) :: row_ptr(ntarg+1),col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1),nquad
      integer *8, intent(in) :: nker,iker
      complex *16, intent(inout) :: wnear(nker,nquad)

      integer *8, allocatable :: ixyzso(:)
      real *8, allocatable :: srcover(:,:),wover(:)
      integer *8 i,npolso,nptso,int8_12
      real *8 thresh

      external fker

      int8_12 = 12

      allocate(ixyzso(npatches+1))
      ixyzso(1) = 1
      do i=1,npatches
        call get_npols(iptype(i),novers(i),npolso)
        ixyzso(i+1) = ixyzso(i) + npolso
      enddo
      nptso = ixyzso(npatches+1)-1

      allocate(srcover(12,nptso),wover(nptso))
      call oversample_geom(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,novers,ixyzso,nptso,srcover)
      call get_qwts(npatches,novers,ixyzso,iptype,nptso,srcover,wover)
      call get_fmm_thresh(int8_12,nptso,srcover,ndtarg,ntarg,targvals,
     1   thresh)

      call zgetnearquadsub_guru(npatches,norders,ixyzs,iptype,novers,
     1   ixyzso,nptso,srcover,wover,ndtarg,ntarg,targvals,thresh,fker,
     2   ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,col_ind,iquad,nquad,
     3   nker,iker,wnear)

      return
      end
c
c
c
c
      subroutine dgetnearquadsub(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,ndtarg,ntarg,targvals,novers,fker,ndd,dpars,
     2   ndz,zpars,ndi,ipars,nnz,row_ptr,col_ind,iquad,nquad,nker,iker,
     3   wnear)
c
c------------------------
c  This subroutine subtracts, in place, the smooth quadrature rule
c  oversampled to novers from a near field quadrature wnear (as
c  computed by dgetnearquad_ggq_guru), for a given kernel, where
c  the near field is specified by the user in row sparse compressed
c  format.
c
c  The oversampled geometry, smooth quadrature weights, and fmm
c  threshold are computed here; see dgetnearquadsub_guru for
c  details.
c
c  Input arguments:
c    - npatches: integer *8
c        number of patches
c    - norders: integer *8(npatches)
c        order of discretization on each patch
c    - ixyzs: integer *8(npatches+1)
c        ixyzs(i) denotes the starting location in srccoefs,
c        and srcvals array corresponding to patch i
c    - iptype: integer *8(npatches)
c        type of patch
c        *  iptype = 1, triangular patch discretized using RV nodes
c        *  iptype = 11, quad patch with GL nodes
c        *  iptype = 12, quad patch with Chebyshev nodes
c    - npts: integer *8
c        total number of discretization points on the boundary
c    - srccoefs: double precision (9,npts)
c        basis expansion coefficients of xyz, dxyz/du, and dxyz/dv
c        on each patch
c    - srcvals: double precision (12,npts)
c        xyz(u,v) and derivative info sampled at the discretization
c        nodes on the surface
c    - ndtarg: integer *8
c        leading dimension of target array
c    - ntarg: integer *8
c        number of targets
c    - targvals: double precision (ndtarg,ntarg)
c        target info. First three components must be x,y,z
c        coordinates
c    - novers: integer *8(npatches)
c        oversampled order on each patch, must match the one used
c        by the fmm in the evaluator
c    - fker, ndd, dpars, ndz, zpars, ndi, ipars, nnz, row_ptr,
c      col_ind, iquad, nquad, nker, iker:
c        as in dgetnearquadsub_guru
c
c  Input/output parameters:
c    - wnear: double precision(nker,nquad)
c        on input, wnear(iker,:) is the near field quadrature; on
c        output, it is the near field quadrature minus the
c        oversampled smooth rule. Other kernels are not accessed
c----------------------------------------------------
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(in) :: npatches,norders(npatches)
      integer *8, intent(in) :: ixyzs(npatches+1),iptype(npatches)
      integer *8, intent(in) :: npts
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
      integer *8, intent(in) :: ndtarg,ntarg
      real *8, intent(in) :: targvals(ndtarg,ntarg)
      integer *8, intent(in) :: novers(npatches)
      integer *8, intent(in) :: ndd,ndz,ndi
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      integer *8, intent(in) :: ipars(ndi)
      integer *8, intent(in) :: nnz
      integer *8, intent(in) :: row_ptr(ntarg+1),col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1),nquad
      integer *8, intent(in) :: nker,iker
      real *8, intent(inout) :: wnear(nker,nquad)

      integer *8, allocatable :: ixyzso(:)
      real *8, allocatable :: srcover(:,:),wover(:)
      integer *8 i,npolso,nptso,int8_12
      real *8 thresh

      external fker

      int8_12 = 12

      allocate(ixyzso(npatches+1))
      ixyzso(1) = 1
      do i=1,npatches
        call get_npols(iptype(i),novers(i),npolso)
        ixyzso(i+1) = ixyzso(i) + npolso
      enddo
      nptso = ixyzso(npatches+1)-1

      allocate(srcover(12,nptso),wover(nptso))
      call oversample_geom(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,novers,ixyzso,nptso,srcover)
      call get_qwts(npatches,novers,ixyzso,iptype,nptso,srcover,wover)
      call get_fmm_thresh(int8_12,nptso,srcover,ndtarg,ntarg,targvals,
     1   thresh)

      call dgetnearquadsub_guru(npatches,norders,ixyzs,iptype,novers,
     1   ixyzso,nptso,srcover,wover,ndtarg,ntarg,targvals,thresh,fker,
     2   ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,col_ind,iquad,nquad,
     3   nker,iker,wnear)

      return
      end
