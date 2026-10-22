      subroutine test_near_corr_bel(nsuccess)
c
c------------------------------------------------------------------
c  This test verifies getnearcorr_lap_bel_log and
c  getnearcorr_helm_bel_hank, for the remainder and the parametrix
c  kernels, with on surface targets carrying the mean curvature in
c  component 13.
c
c  For a random density sigma, the near quadrature minus the
c  corrections, applied to sigma, must equal the oversampled smooth
c  rule applied to sigma, i.e. the sum over the oversampled nodes
c  y_l of the near patches of K(x,y_l) w_l sigma_over(l).
c------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)

      integer *8 int8_1, int8_2, int8_12
      real *8 c0(3)

      integer *8, allocatable :: norders(:), iptype(:), ixyzs(:)
      real *8, allocatable :: srcvals(:,:), srccoefs(:,:)
      real *8, allocatable :: cms(:,:), rads(:), rad_near(:)
      real *8, allocatable :: targs(:,:), uvs_targ(:,:)
      integer *8, allocatable :: ipatch_id(:)

      integer *8, allocatable :: row_ptr(:), col_ind(:), iquad(:)
      integer *8, allocatable :: novers(:), ixyzso(:)
      real *8, allocatable :: srcover(:,:), wover(:)

      real *8, allocatable :: dwq(:), dwc(:), dsig(:), dsigover(:)
      complex *16, allocatable :: zwq(:), zwc(:)
      complex *16, allocatable :: zsig(:), zsigover(:)

      real *8 dpars(1)
      complex *16 zk, zk0, zval, zlhs, zrhs
      integer *8 ipars(1)

      call prini(6,13)

      int8_1 = 1
      int8_2 = 2
      int8_12 = 12

      nsuccess = 0

c
c   unit sphere of triangles, on surface targets
c
      a = 1.0d0
      na = 2
      norder = 4
      c0(1) = 0
      c0(2) = 0
      c0(3) = 0
      call get_sphere_npat_mem(a, na, c0, norder, int8_1, npatches,
     1   npts)
      allocate(srcvals(12,npts), srccoefs(9,npts))
      allocate(norders(npatches), iptype(npatches), ixyzs(npatches+1))
      call get_sphere_npat(a, na, c0, norder, int8_1, npatches, npts,
     1   norders, ixyzs, iptype, srccoefs, srcvals)

      ndtarg = 13
      ntarg = npts
      allocate(targs(ndtarg,ntarg), ipatch_id(ntarg))
      allocate(uvs_targ(2,ntarg))
      do i=1,ntarg
        do k=1,12
          targs(k,i) = srcvals(k,i)
        enddo
        targs(13,i) = 1.0d0
      enddo
      call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts,
     1   ipatch_id, uvs_targ)

c
c   near field, patch dependent oversampling, fmm threshold
c
      call get_rfacs(norder, int8_1, rfac, rfac0)
      allocate(cms(3,npatches), rads(npatches), rad_near(npatches))
      call get_centroid_rads(npatches, norders, ixyzs, iptype, npts,
     1   srccoefs, cms, rads)
      do i=1,npatches
        rad_near(i) = rads(i)*rfac
      enddo
      call findnearmem(cms, npatches, rad_near, ndtarg, targs, ntarg,
     1   nnz)
      allocate(row_ptr(ntarg+1), col_ind(nnz), iquad(nnz+1))
      call findnear(cms, npatches, rad_near, ndtarg, targs, ntarg,
     1   row_ptr, col_ind)
      call get_iquad_rsc(npatches, ixyzs, ntarg, nnz, row_ptr, col_ind,
     1   iquad)
      nquad = iquad(nnz+1) - 1

      allocate(novers(npatches), ixyzso(npatches+1))
      ixyzso(1) = 1
      do i=1,npatches
        novers(i) = norders(i) + 2 + mod(i,2)
        call get_npols(iptype(i), novers(i), npolso)
        ixyzso(i+1) = ixyzso(i) + npolso
      enddo
      nptso = ixyzso(npatches+1) - 1
      allocate(srcover(12,nptso), wover(nptso))
      call oversample_geom(npatches, norders, ixyzs, iptype, npts,
     1   srccoefs, srcvals, novers, ixyzso, nptso, srcover)
      call get_qwts(npatches, novers, ixyzso, iptype, nptso, srcover,
     1   wover)
      call get_fmm_thresh(int8_12, nptso, srcover, ndtarg, ntarg, targs,
     1   thresh)

      allocate(dsig(npts), dsigover(nptso))
      allocate(zsig(npts), zsigover(nptso))
      do i=1,npts
        dsig(i) = hkrand(0)
        zsig(i) = hkrand(0) + hkrand(0)*(0.0d0,1.0d0)
      enddo
      call oversample_fun_surf(int8_1, npatches, norders, ixyzs,
     1   iptype, npts, dsig, novers, ixyzso, nptso, dsigover)
      call oversample_fun_surf(int8_2, npatches, norders, ixyzs,
     1   iptype, npts, zsig, novers, ixyzso, nptso, zsigover)

      eps = 1.0d-6
      iquadtype = 1
      ndd = 0
      ndz = 1
      ndi = 0
      zk0 = 0
      zk = 1.3d0
      allocate(dwq(nquad), dwc(nquad), zwq(nquad), zwc(nquad))

c
c   Laplace-Beltrami, remainder (iktype = 1) and parametrix (2)
c
      errmax = 0
      do iktype=1,2
        call getnearquad_lap_bel_log(npatches, norders, ixyzs, iptype,
     1     npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id,
     2     uvs_targ, eps, iktype, iquadtype, nnz, row_ptr, col_ind,
     3     iquad, rfac0, nquad, dwq)
        call getnearcorr_lap_bel_log(npatches, norders, ixyzs, iptype,
     1     npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id,
     2     uvs_targ, eps, iktype, iquadtype, nnz, row_ptr, col_ind,
     3     iquad, rfac0, nquad, novers, dwc)

        err = 0
        rmax = 0
        do i=1,ntarg
          dlhs = 0
          drhs = 0
          do j=row_ptr(i),row_ptr(i+1)-1
            jpatch = col_ind(j)
            npols = ixyzs(jpatch+1) - ixyzs(jpatch)
            do m=1,npols
              dlhs = dlhs + (dwq(iquad(j)+m-1) - dwc(iquad(j)+m-1))*
     1           dsig(ixyzs(jpatch)+m-1)
            enddo
            do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
              rr = sqrt((srcover(1,l)-targs(1,i))**2 +
     1           (srcover(2,l)-targs(2,i))**2 +
     2           (srcover(3,l)-targs(3,i))**2)
              if(rr.le.thresh) cycle
              if(iktype.eq.1) call lap_bel_res(srcover(1,l), ndtarg,
     1           targs(1,i), ndd, dpars, ndz, zk0, ndi, ipars, dval)
              if(iktype.eq.2) call lap_bel_log(srcover(1,l), ndtarg,
     1           targs(1,i), ndd, dpars, ndz, zk0, ndi, ipars, dval)
              drhs = drhs + dval*wover(l)*dsigover(l)
            enddo
          enddo
          err = max(err, abs(dlhs-drhs))
          rmax = max(rmax, abs(drhs))
        enddo
        errmax = max(errmax, err/rmax)
      enddo

      call prin2('laplace-beltrami: relative error in subtraction=*',
     1   errmax, 1)
      if(errmax.lt.1.0d-12) nsuccess = nsuccess + 1

c
c   Helmholtz-Beltrami, remainder (iktype = 1) and parametrix (2)
c
      errmax = 0
      do iktype=1,2
        call getnearquad_helm_bel_hank(npatches, norders, ixyzs, iptype,
     1     npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id,
     2     uvs_targ, eps, zk, iktype, iquadtype, nnz, row_ptr, col_ind,
     3     iquad, rfac0, nquad, zwq)
        call getnearcorr_helm_bel_hank(npatches, norders, ixyzs, iptype,
     1     npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id,
     2     uvs_targ, eps, zk, iktype, iquadtype, nnz, row_ptr, col_ind,
     3     iquad, rfac0, nquad, novers, zwc)

        err = 0
        rmax = 0
        do i=1,ntarg
          zlhs = 0
          zrhs = 0
          do j=row_ptr(i),row_ptr(i+1)-1
            jpatch = col_ind(j)
            npols = ixyzs(jpatch+1) - ixyzs(jpatch)
            do m=1,npols
              zlhs = zlhs + (zwq(iquad(j)+m-1) - zwc(iquad(j)+m-1))*
     1           zsig(ixyzs(jpatch)+m-1)
            enddo
            do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
              rr = sqrt((srcover(1,l)-targs(1,i))**2 +
     1           (srcover(2,l)-targs(2,i))**2 +
     2           (srcover(3,l)-targs(3,i))**2)
              if(rr.le.thresh) cycle
              if(iktype.eq.1) call helm_bel_res(srcover(1,l), ndtarg,
     1           targs(1,i), ndd, dpars, ndz, zk, ndi, ipars, zval)
              if(iktype.eq.2) call helm_bel_hank(srcover(1,l), ndtarg,
     1           targs(1,i), ndd, dpars, ndz, zk, ndi, ipars, zval)
              zrhs = zrhs + zval*wover(l)*zsigover(l)
            enddo
          enddo
          err = max(err, abs(zlhs-zrhs))
          rmax = max(rmax, abs(zrhs))
        enddo
        errmax = max(errmax, err/rmax)
      enddo

      call prin2('helmholtz-beltrami: relative error in subtraction=*',
     1   errmax, 1)
      if(errmax.lt.1.0d-12) nsuccess = nsuccess + 1

      return
      end
