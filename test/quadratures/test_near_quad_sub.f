      subroutine test_near_quad_sub(nsuccess)
c
c------------------------------------------------------------------
c  This test verifies z/dgetnearquadsub_guru on a mesh with
c  triangles and quads, and patch dependent oversampling orders.
c
c  For a random density sigma, the amount subtracted from the
c  near quadrature, applied to sigma, must equal the oversampled
c  smooth rule applied to sigma, i.e. the sum over the oversampled
c  nodes y_l of the near patches of K(x,y_l) w_l sigma_over(l),
c  where sigma_over is sigma interpolated with oversample_fun_surf
c  (as done by the add-and-subtract evaluators). Only row iker of
c  wnear(nker,nquad) may change.
c------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)

      integer *8 int8_1, int8_2, int8_3, int8_11, int8_12
      real *8 c01(3), c02(3)

      integer *8, allocatable :: norders1(:), norders2(:)
      integer *8, allocatable :: iptype1(:), iptype2(:)
      integer *8, allocatable :: ixyzs1(:), ixyzs2(:)
      real *8, allocatable :: srcvals1(:,:), srcvals2(:,:)
      real *8, allocatable :: srccoefs1(:,:), srccoefs2(:,:)

      integer *8, allocatable :: norders(:), iptype(:), ixyzs(:)
      real *8, allocatable :: srcvals(:,:), srccoefs(:,:)
      real *8, allocatable :: cms(:,:), rads(:), rad_near(:)
      real *8, allocatable :: targs(:,:)

      integer *8, allocatable :: row_ptr(:), col_ind(:), iquad(:)
      integer *8, allocatable :: novers(:), ixyzso(:)
      real *8, allocatable :: srcover(:,:), wover(:)

      complex *16, allocatable :: zwnear(:,:), zwnear0(:,:)
      complex *16, allocatable :: zsig(:), zsigover(:)
      real *8, allocatable :: dwnear(:,:), dwnear0(:,:)
      real *8, allocatable :: dsig(:), dsigover(:)

      real *8 dpars(2)
      complex *16 zpars(1), zval, zlhs, zrhs
      integer *8 ipars(1)
      procedure (), pointer :: fker
      external h3d_slp, l3d_comb

      call prini(6,13)

      int8_1 = 1
      int8_2 = 2
      int8_3 = 3
      int8_11 = 11
      int8_12 = 12

      done = 1
      nsuccess = 0

c
c   two nearby spheres: triangles of order 4, and GL quads of
c   order 5, so that near fields contain both patch types
c
      a = 1.0d0
      na = 1
      norder1 = 4
      norder2 = 5

      c01(1) = 0
      c01(2) = 0
      c01(3) = 0

      c02(1) = 2.3d0
      c02(2) = 0
      c02(3) = 0

      call get_sphere_npat_mem(a, na, c01, norder1, int8_1, npatches1,
     1   npts1)
      allocate(norders1(npatches1), iptype1(npatches1))
      allocate(ixyzs1(npatches1+1))
      allocate(srcvals1(12,npts1), srccoefs1(9,npts1))
      call get_sphere_npat(a, na, c01, norder1, int8_1, npatches1,
     1   npts1, norders1, ixyzs1, iptype1, srccoefs1, srcvals1)

      call get_sphere_npat_mem(a, na, c02, norder2, int8_11, npatches2,
     1   npts2)
      allocate(norders2(npatches2), iptype2(npatches2))
      allocate(ixyzs2(npatches2+1))
      allocate(srcvals2(12,npts2), srccoefs2(9,npts2))
      call get_sphere_npat(a, na, c02, norder2, int8_11, npatches2,
     1   npts2, norders2, ixyzs2, iptype2, srccoefs2, srcvals2)

      npatches = npatches1 + npatches2
      npts = npts1 + npts2

      allocate(norders(npatches), iptype(npatches), ixyzs(npatches+1))
      allocate(srcvals(12,npts), srccoefs(9,npts))

      do i=1,npatches1
        norders(i) = norders1(i)
        iptype(i) = iptype1(i)
        ixyzs(i) = ixyzs1(i)
      enddo
      do i=1,npatches2+1
        if(i.le.npatches2) then
          norders(npatches1+i) = norders2(i)
          iptype(npatches1+i) = iptype2(i)
        endif
        ixyzs(npatches1+i) = npts1 + ixyzs2(i)
      enddo

      do i=1,npts1
        srcvals(1:12,i) = srcvals1(1:12,i)
        srccoefs(1:9,i) = srccoefs1(1:9,i)
      enddo
      do i=1,npts2
        srcvals(1:12,npts1+i) = srcvals2(1:12,i)
        srccoefs(1:9,npts1+i) = srccoefs2(1:9,i)
      enddo

c
c   targets: all discretization nodes, and points slightly off
c   every third node
c
      noff = (npts+2)/3
      ntarg = npts + noff
      ndtarg = 3
      allocate(targs(ndtarg,ntarg))
      do i=1,npts
        targs(1:3,i) = srcvals(1:3,i)
      enddo
      do i=1,noff
        j = 3*(i-1) + 1
        targs(1:3,npts+i) = srcvals(1:3,j) + 0.05d0*srcvals(10:12,j)
      enddo

c
c   near field
c
      allocate(cms(3,npatches), rads(npatches), rad_near(npatches))
      call get_centroid_rads(npatches, norders, ixyzs, iptype, npts,
     1   srccoefs, cms, rads)
      do i=1,npatches
        rad_near(i) = 1.5d0*rads(i)
      enddo

      call findnearmem(cms, npatches, rad_near, ndtarg, targs, ntarg,
     1   nnz)
      allocate(row_ptr(ntarg+1), col_ind(nnz), iquad(nnz+1))
      call findnear(cms, npatches, rad_near, ndtarg, targs, ntarg,
     1   row_ptr, col_ind)
      call get_iquad_rsc(npatches, ixyzs, ntarg, nnz, row_ptr, col_ind,
     1   iquad)
      nquad = iquad(nnz+1) - 1

c
c   patch dependent oversampling, and oversampled geometry
c
      allocate(novers(npatches), ixyzso(npatches+1))
      ixyzso(1) = 1
      do i=1,npatches
        novers(i) = norders(i) + 1 + mod(i,int8_3)
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

      nker = 3

c
c   complex: helmholtz single layer, update row 2
c
      allocate(zwnear(nker,nquad), zwnear0(nker,nquad))
      do i=1,nquad
        do k=1,nker
          zwnear(k,i) = hkrand(0) + hkrand(0)*(0.0d0,1.0d0)
          zwnear0(k,i) = zwnear(k,i)
        enddo
      enddo

      iker = 2
      zpars(1) = (1.1d0,0.1d0)
      ndd = 0
      ndz = 1
      ndi = 0
      fker => h3d_slp
      call zgetnearquadsub_guru(npatches, norders, ixyzs, iptype,
     1   novers, ixyzso, nptso, srcover, wover, ndtarg, ntarg, targs,
     2   thresh, fker, ndd, dpars, ndz, zpars, ndi, ipars, nnz,
     3   row_ptr, col_ind, iquad, nquad, nker, iker, zwnear)

      allocate(zsig(npts), zsigover(nptso))
      do i=1,npts
        zsig(i) = hkrand(0) + hkrand(0)*(0.0d0,1.0d0)
      enddo
      call oversample_fun_surf(int8_2, npatches, norders, ixyzs,
     1   iptype, npts, zsig, novers, ixyzso, nptso, zsigover)

      nuntouched = 0
      do i=1,nquad
        if(zwnear(1,i).ne.zwnear0(1,i)) nuntouched = nuntouched + 1
        if(zwnear(3,i).ne.zwnear0(3,i)) nuntouched = nuntouched + 1
      enddo

      err = 0
      rmax = 0
      do i=1,ntarg
        zlhs = 0
        zrhs = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1) - ixyzs(jpatch)
          do m=1,npols
            zlhs = zlhs + (zwnear0(iker,iquad(j)+m-1) -
     1         zwnear(iker,iquad(j)+m-1))*zsig(ixyzs(jpatch)+m-1)
          enddo
          do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
            rr = sqrt((srcover(1,l)-targs(1,i))**2 +
     1         (srcover(2,l)-targs(2,i))**2 +
     2         (srcover(3,l)-targs(3,i))**2)
            if(rr.le.thresh) cycle
            call h3d_slp(srcover(1,l), ndtarg, targs(1,i), ndd, dpars,
     1         ndz, zpars, ndi, ipars, zval)
            zrhs = zrhs + zval*wover(l)*zsigover(l)
          enddo
        enddo
        err = max(err, abs(zlhs-zrhs))
        rmax = max(rmax, abs(zrhs))
      enddo
      err = err/rmax

      call prin2('complex: relative error in subtraction=*', err, 1)
      call prinf('complex: other kernels modified=*', nuntouched, 1)
      if(err.lt.1.0d-12.and.nuntouched.eq.0) nsuccess = nsuccess + 1

c
c   real: laplace combined field, update row 3
c
      allocate(dwnear(nker,nquad), dwnear0(nker,nquad))
      do i=1,nquad
        do k=1,nker
          dwnear(k,i) = hkrand(0)
          dwnear0(k,i) = dwnear(k,i)
        enddo
      enddo

      iker = 3
      dpars(1) = 0.7d0
      dpars(2) = -1.3d0
      ndd = 2
      ndz = 0
      ndi = 0
      fker => l3d_comb
      call dgetnearquadsub_guru(npatches, norders, ixyzs, iptype,
     1   novers, ixyzso, nptso, srcover, wover, ndtarg, ntarg, targs,
     2   thresh, fker, ndd, dpars, ndz, zpars, ndi, ipars, nnz,
     3   row_ptr, col_ind, iquad, nquad, nker, iker, dwnear)

      allocate(dsig(npts), dsigover(nptso))
      do i=1,npts
        dsig(i) = hkrand(0)
      enddo
      call oversample_fun_surf(int8_1, npatches, norders, ixyzs,
     1   iptype, npts, dsig, novers, ixyzso, nptso, dsigover)

      nuntouched = 0
      do i=1,nquad
        if(dwnear(1,i).ne.dwnear0(1,i)) nuntouched = nuntouched + 1
        if(dwnear(2,i).ne.dwnear0(2,i)) nuntouched = nuntouched + 1
      enddo

      err = 0
      rmax = 0
      do i=1,ntarg
        dlhs = 0
        drhs = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1) - ixyzs(jpatch)
          do m=1,npols
            dlhs = dlhs + (dwnear0(iker,iquad(j)+m-1) -
     1         dwnear(iker,iquad(j)+m-1))*dsig(ixyzs(jpatch)+m-1)
          enddo
          do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
            rr = sqrt((srcover(1,l)-targs(1,i))**2 +
     1         (srcover(2,l)-targs(2,i))**2 +
     2         (srcover(3,l)-targs(3,i))**2)
            if(rr.le.thresh) cycle
            call l3d_comb(srcover(1,l), ndtarg, targs(1,i), ndd, dpars,
     1         ndz, zpars, ndi, ipars, dval)
            drhs = drhs + dval*wover(l)*dsigover(l)
          enddo
        enddo
        err = max(err, abs(dlhs-drhs))
        rmax = max(rmax, abs(drhs))
      enddo
      err = err/rmax

      call prin2('real: relative error in subtraction=*', err, 1)
      call prinf('real: other kernels modified=*', nuntouched, 1)
      if(err.lt.1.0d-12.and.nuntouched.eq.0) nsuccess = nsuccess + 1

      return
      end
