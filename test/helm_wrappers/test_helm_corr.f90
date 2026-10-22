      subroutine test_helm_corr(nsuccess)
!
!  Tests of the near quadrature corrections (getnearcorr_*, ifcorr = 1)
!  for the Helmholtz wrappers on a stellarator.
!
!  Applying the operators and solving with corrections must reproduce
!  the add and subtract evaluators to rounding (tests 1-14), and the
!  simple solvers, which use the corrections, must solve exterior
!  problems with a known solution (tests 15-16).
!
!  Output:
!    - nsuccess: number of successful tests out of 16
!
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(out) :: nsuccess

      real *8, allocatable :: srcvals(:,:), srccoefs(:,:), targs(:,:)
      real *8, allocatable :: cms(:,:), rads(:), rad_near(:)
      integer *8, allocatable :: norders(:), ixyzs(:), iptype(:)
      integer *8, allocatable :: row_ptr(:), col_ind(:), iquad(:)
      integer *8, allocatable :: row_ptr_t(:), col_ind_t(:), iquad_t(:)
      integer *8, allocatable :: novers(:), ixyzso(:)
      integer *8, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      real *8, allocatable :: srcover(:,:), wover(:)
      complex *16, allocatable :: wnear(:,:), wcorr(:,:)
      complex *16, allocatable :: wnear_t(:,:), wcorr_t(:,:)
      complex *16, allocatable :: sig(:), p0(:), p1(:)
      complex *16, allocatable :: sig2(:,:), p20(:,:), p21(:,:)
      complex *16, allocatable :: pt0(:), pt1(:)
      complex *16, allocatable :: rhs(:), soln0(:), soln1(:)
      complex *16, allocatable :: sik0(:), sik1(:), zpimp(:), zlams(:)
      real *8, allocatable :: errs(:)
      real *8 dpars(1), work(1), xyz_in(3), xyz_out(3)
      complex *16 zpars(6), zk, ima, pex, pnum(1)
      integer *8 ipars(1), nuv(2), iptype0
      integer *8 int8_0, int8_1, int8_3, int8_12
      integer *8 ipatch_out(1)
      real *8 uvs_out(2,1)
      real *8 zreldiff
      data ima/(0.0d0,1.0d0)/

      int8_0 = 0
      int8_1 = 1
      int8_3 = 3
      int8_12 = 12
      nsuccess = 0

!
!  geometry
!
      norder = 4
      nuv(1) = 4
      nuv(2) = 12
      iptype0 = 1
      call get_stellarator_npat_mem(nuv, norder, iptype0, npatches, &
        npts)
      allocate(srcvals(12,npts), srccoefs(9,npts))
      allocate(norders(npatches), ixyzs(npatches+1), iptype(npatches))
      call get_stellarator_npat(nuv, norder, iptype0, npatches, npts, &
        norders, ixyzs, iptype, srccoefs, srcvals)

      xyz_in(1) = -4.501d0
      xyz_in(2) = 1.7d-3
      xyz_in(3) = 0.00001d0

      xyz_out(1) = -3.5d0
      xyz_out(2) = 3.1d0
      xyz_out(3) = 20.1d0

      eps = 1.0d-6
      zk = 1.1d0

!
!  near field of the on surface targets, oversampling orders
!
      allocate(cms(3,npatches), rads(npatches), rad_near(npatches))
      call get_centroid_rads(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, cms, rads)
      call get_rfacs(norder, iptype0, rfac, rfac0)
      do i = 1,npatches
        rad_near(i) = rads(i)*rfac
      enddo

      call findnearmem(cms, npatches, rad_near, int8_12, srcvals, npts, &
        nnz)
      allocate(row_ptr(npts+1), col_ind(nnz), iquad(nnz+1))
      call findnear(cms, npatches, rad_near, int8_12, srcvals, npts, &
        row_ptr, col_ind)
      call get_iquad_rsc(npatches, ixyzs, npts, nnz, row_ptr, col_ind, &
        iquad)
      nquad = iquad(nnz+1) - 1

      allocate(novers(npatches), ixyzso(npatches+1))
      ixyzso(1) = 1
      do i = 1,npatches
        novers(i) = norders(i) + 3
        call get_npols(iptype(i), novers(i), npolso)
        ixyzso(i+1) = ixyzso(i) + npolso
      enddo
      nptso = ixyzso(npatches+1) - 1
      allocate(srcover(12,nptso), wover(nptso))
      call oversample_geom(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, srcvals, novers, ixyzso, nptso, srcover)
      call get_qwts(npatches, novers, ixyzso, iptype, nptso, srcover, &
        wover)

!
!  off surface targets, at a tenth of the local patch size from
!  every fifth node, on either side of the surface
!
      ntarg = (npts+4)/5
      allocate(targs(3,ntarg), ipatch_id(ntarg), uvs_targ(2,ntarg))
      do i = 1,ntarg
        j = 5*(i-1) + 1
        ipatch = 1
        do k = 1,npatches
          if (j.ge.ixyzs(k).and.j.lt.ixyzs(k+1)) ipatch = k
        enddo
        h = 0.1d0*rads(ipatch)*(-1)**i
        targs(1:3,i) = srcvals(1:3,j) + h*srcvals(10:12,j)
        ipatch_id(i) = -1
        uvs_targ(1:2,i) = 0
      enddo
      call findnearmem(cms, npatches, rad_near, int8_3, targs, ntarg, &
        nnz_t)
      allocate(row_ptr_t(ntarg+1), col_ind_t(nnz_t), iquad_t(nnz_t+1))
      call findnear(cms, npatches, rad_near, int8_3, targs, ntarg, &
        row_ptr_t, col_ind_t)
      call get_iquad_rsc(npatches, ixyzs, ntarg, nnz_t, row_ptr_t, &
        col_ind_t, iquad_t)
      nquad_t = iquad_t(nnz_t+1) - 1

      allocate(wnear(6,nquad), wcorr(6,nquad))
      allocate(wnear_t(2,nquad_t), wcorr_t(2,nquad_t))
      allocate(sig(npts), p0(npts), p1(npts))
      allocate(sig2(2,npts), p20(2,npts), p21(2,npts))
      allocate(pt0(ntarg), pt1(ntarg))
      allocate(rhs(npts), soln0(npts), soln1(npts), sik0(npts))
      allocate(sik1(npts))
      do i = 1,npts
        sig(i) = hkrand(0) - 0.5d0 + ima*(hkrand(0) - 0.5d0)
        sig2(1,i) = hkrand(0) - 0.5d0 + ima*(hkrand(0) - 0.5d0)
        sig2(2,i) = hkrand(0) - 0.5d0 + ima*(hkrand(0) - 0.5d0)
        rhs(i) = hkrand(0) - 0.5d0 + ima*(hkrand(0) - 0.5d0)
      enddo

      iquadtype = 1
      ndd = 0
      ndi = 0
      lwork = 0
      idensflag = 0
      ipotflag = 1
      ndim = 1
      ndim_p = 1
      numit = 200
      eps_gmres = 1.0d-12
      ifinout = 1
      allocate(errs(numit+1))

!
!  1-2. helm_comb_dir: matvec and solver
!
      zpars(1) = zk
      zpars(2) = -ima*zk
      zpars(3) = 0.7d0
      ndz = 3
      nker = 1
      call getnearquad_helm_comb_dir(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)
      call getnearcorr_helm_comb_dir(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, novers, wcorr)
      call lpcomp_helm_comb_dir_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, &
        novers, nptso, ixyzso, srcover, wover, lwork, work, ndim, sig, p0)
      call lpcomp_helm_comb_dir_corr(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wcorr, &
        novers, nptso, ixyzso, srcover, wover, lwork, work, ndim, sig, p1)
      call report_zdiff('helm_comb_dir matvec*', npts, p0, p1, 1.0d-12, &
        nsuccess)
      call helm_comb_dir_solver_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, numit, ifinout, &
        rhs, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, int8_0, &
        novers, nptso, ixyzso, srcover, wover, eps_gmres, niter, errs, &
        rres, soln0)
      call helm_comb_dir_solver_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, numit, ifinout, &
        rhs, nnz, row_ptr, col_ind, iquad, nquad, nker, wcorr, int8_1, &
        novers, nptso, ixyzso, srcover, wover, eps_gmres, niter, errs, &
        rres, soln1)
      call report_zdiff('helm_comb_dir_solver_guru*', npts, soln0, &
        soln1, 1.0d-9, nsuccess)

!
!  3. helm_comb_dir: off surface evaluator
!
      call getnearquad_helm_comb_dir_eval(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, int8_3, ntarg, targs, &
        ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz_t, row_ptr_t, &
        col_ind_t, iquad_t, rfac0, nquad_t, wnear_t)
      call getnearcorr_helm_comb_dir_eval(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, int8_3, ntarg, targs, &
        ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz_t, row_ptr_t, &
        col_ind_t, iquad_t, rfac0, nquad_t, novers, wcorr_t)
      call helm_comb_dir_eval_addsub(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, int8_3, ntarg, targs, eps, ndd, dpars, &
        ndz, zpars, ndi, ipars, nnz_t, row_ptr_t, col_ind_t, iquad_t, &
        nquad_t, nker, wnear_t, int8_0, novers, nptso, ixyzso, srcover, &
        wover, lwork, work, idensflag, ndim, sig, ipotflag, ndim_p, pt0)
      call helm_comb_dir_eval_addsub(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, int8_3, ntarg, targs, eps, ndd, dpars, &
        ndz, zpars, ndi, ipars, nnz_t, row_ptr_t, col_ind_t, iquad_t, &
        nquad_t, nker, wcorr_t, int8_1, novers, nptso, ixyzso, srcover, &
        wover, lwork, work, idensflag, ndim, sig, ipotflag, ndim_p, pt1)
      call report_zdiff('helm_comb_dir off surface eval*', ntarg, pt0, &
        pt1, 1.0d-12, nsuccess)

!
!  4. helm_comb_cprime: on surface evaluator
!
      zpars(2) = 1.0d0
      zpars(3) = 0.5d0
      call getnearquad_helm_comb_cprime(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, &
        row_ptr, col_ind, iquad, rfac0, nquad, wnear)
      call getnearcorr_helm_comb_cprime(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, &
        row_ptr, col_ind, iquad, rfac0, nquad, novers, wcorr)
      call helm_comb_cprime_eval_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, int8_12, npts, srcvals, eps, &
        ndd, dpars, ndz, zpars, ndi, ipars, nnz, row_ptr, col_ind, &
        iquad, nquad, nker, wnear, int8_0, novers, nptso, ixyzso, &
        srcover, wover, lwork, work, idensflag, ndim, sig, ipotflag, &
        ndim_p, p0)
      call helm_comb_cprime_eval_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, int8_12, npts, srcvals, eps, &
        ndd, dpars, ndz, zpars, ndi, ipars, nnz, row_ptr, col_ind, &
        iquad, nquad, nker, wcorr, int8_1, novers, nptso, ixyzso, &
        srcover, wover, lwork, work, idensflag, ndim, sig, ipotflag, &
        ndim_p, p1)
      call report_zdiff('helm_comb_cprime on surface eval*', npts, p0, &
        p1, 1.0d-12, nsuccess)

!
!  5-6. helm_s_neu: matvec and solver
!
      ndz = 1
      call getnearquad_helm_s_neu(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)
      call getnearcorr_helm_s_neu(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, novers, wcorr)
      call lpcomp_helm_s_neu_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, &
        novers, nptso, ixyzso, srcover, wover, lwork, work, ndim, sig, p0)
      call lpcomp_helm_s_neu_corr(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wcorr, &
        novers, nptso, ixyzso, srcover, wover, lwork, work, ndim, sig, p1)
      call report_zdiff('helm_s_neu matvec*', npts, p0, p1, 1.0d-12, &
        nsuccess)
      call helm_s_neu_solver_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, numit, ifinout, &
        rhs, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, int8_0, &
        novers, nptso, ixyzso, srcover, wover, eps_gmres, niter, errs, &
        rres, soln0)
      call helm_s_neu_solver_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, numit, ifinout, &
        rhs, nnz, row_ptr, col_ind, iquad, nquad, nker, wcorr, int8_1, &
        novers, nptso, ixyzso, srcover, wover, eps_gmres, niter, errs, &
        rres, soln1)
      call report_zdiff('helm_s_neu_solver_guru*', npts, soln0, soln1, &
        1.0d-9, nsuccess)

!
!  7. helm_s: off surface evaluator
!
      call getnearquad_helm_s_eval(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, int8_3, ntarg, targs, &
        ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz_t, row_ptr_t, &
        col_ind_t, iquad_t, rfac0, nquad_t, wnear_t)
      call getnearcorr_helm_s_eval(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, int8_3, ntarg, targs, &
        ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz_t, row_ptr_t, &
        col_ind_t, iquad_t, rfac0, nquad_t, novers, wcorr_t)
      call helm_s_eval_addsub(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, int8_3, ntarg, targs, eps, ndd, dpars, &
        ndz, zpars, ndi, ipars, nnz_t, row_ptr_t, col_ind_t, iquad_t, &
        nquad_t, nker, wnear_t, int8_0, novers, nptso, ixyzso, srcover, &
        wover, lwork, work, idensflag, ndim, sig, ipotflag, ndim_p, pt0)
      call helm_s_eval_addsub(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, int8_3, ntarg, targs, eps, ndd, dpars, &
        ndz, zpars, ndi, ipars, nnz_t, row_ptr_t, col_ind_t, iquad_t, &
        nquad_t, nker, wcorr_t, int8_1, novers, nptso, ixyzso, srcover, &
        wover, lwork, work, idensflag, ndim, sig, ipotflag, ndim_p, pt1)
      call report_zdiff('helm_s off surface eval*', ntarg, pt0, pt1, &
        1.0d-12, nsuccess)

!
!  8-9. helm_rpcomb_neu: matvec and solver (solution and S_{i|k|} of it)
!
      zpars(2) = 1.3d0
      ndz = 2
      nker = 4
      call getnearquad_helm_rpcomb_neu(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, &
        row_ptr, col_ind, iquad, rfac0, nquad, wnear)
      call getnearcorr_helm_rpcomb_neu(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, &
        row_ptr, col_ind, iquad, rfac0, nquad, novers, wcorr)
      call lpcomp_helm_rpcomb_neu_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, &
        novers, nptso, ixyzso, srcover, wover, lwork, work, ndim, sig, p0)
      call lpcomp_helm_rpcomb_neu_corr(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wcorr, &
        novers, nptso, ixyzso, srcover, wover, lwork, work, ndim, sig, p1)
      call report_zdiff('helm_rpcomb_neu matvec*', npts, p0, p1, &
        1.0d-12, nsuccess)
      call helm_rpcomb_neu_solver_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, numit, ifinout, &
        rhs, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, int8_0, &
        novers, nptso, ixyzso, srcover, wover, eps_gmres, niter, errs, &
        rres, soln0, sik0)
      call helm_rpcomb_neu_solver_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, numit, ifinout, &
        rhs, nnz, row_ptr, col_ind, iquad, nquad, nker, wcorr, int8_1, &
        novers, nptso, ixyzso, srcover, wover, eps_gmres, niter, errs, &
        rres, soln1, sik1)
      err = max(zreldiff(npts, soln0, soln1), zreldiff(npts, sik0, sik1))
      call prin2('helm_rpcomb_neu_solver_guru*', err, 1)
      if (err.lt.1.0d-9) nsuccess = nsuccess + 1

!
!  10. helm_rpcomb: off surface evaluator, two densities
!
      nker = 2
      ndim_s = 2
      call getnearquad_helm_rpcomb_eval(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, int8_3, ntarg, targs, &
        ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz_t, row_ptr_t, &
        col_ind_t, iquad_t, rfac0, nquad_t, wnear_t)
      call getnearcorr_helm_rpcomb_eval(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, int8_3, ntarg, targs, &
        ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz_t, row_ptr_t, &
        col_ind_t, iquad_t, rfac0, nquad_t, novers, wcorr_t)
      call helm_rpcomb_eval_addsub(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, int8_3, ntarg, targs, eps, ndd, dpars, &
        ndz, zpars, ndi, ipars, nnz_t, row_ptr_t, col_ind_t, iquad_t, &
        nquad_t, nker, wnear_t, int8_0, novers, nptso, ixyzso, srcover, &
        wover, lwork, work, idensflag, ndim_s, sig2, ipotflag, ndim_p, pt0)
      call helm_rpcomb_eval_addsub(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, int8_3, ntarg, targs, eps, ndd, dpars, &
        ndz, zpars, ndi, ipars, nnz_t, row_ptr_t, col_ind_t, iquad_t, &
        nquad_t, nker, wcorr_t, int8_1, novers, nptso, ixyzso, srcover, &
        wover, lwork, work, idensflag, ndim_s, sig2, ipotflag, ndim_p, pt1)
      call report_zdiff('helm_rpcomb off surface eval*', ntarg, pt0, &
        pt1, 1.0d-12, nsuccess)

!
!  11-12. helm_rpcomb_imp: matvec (impedance in zpars(3:)) and solver
!
      nker = 6
      allocate(zpimp(2+npts), zlams(npts))
      zpimp(1) = zk
      zpimp(2) = zpars(2)
      do i = 1,npts
        zlams(i) = 1.0d0 + 0.3d0*hkrand(0) + ima*0.2d0
        zpimp(2+i) = zlams(i)
      enddo
      ndz = 2 + npts
      call getnearquad_helm_rpcomb_imp(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, &
        row_ptr, col_ind, iquad, rfac0, nquad, wnear)
      call getnearcorr_helm_rpcomb_imp(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, &
        row_ptr, col_ind, iquad, rfac0, nquad, novers, wcorr)
      call lpcomp_helm_rpcomb_imp_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpimp, &
        int8_1, ifinout, nnz, row_ptr, col_ind, iquad, nquad, nker, &
        wnear, novers, nptso, ixyzso, srcover, wover, lwork, work, &
        ndim, sig, p0)
      call lpcomp_helm_rpcomb_imp_corr(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpimp, &
        int8_1, ifinout, nnz, row_ptr, col_ind, iquad, nquad, nker, &
        wcorr, novers, nptso, ixyzso, srcover, wover, lwork, work, &
        ndim, sig, p1)
      call report_zdiff('helm_rpcomb_imp matvec*', npts, p0, p1, &
        1.0d-12, nsuccess)
      call helm_rpcomb_imp_solver_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, zlams, numit, &
        ifinout, rhs, nnz, row_ptr, col_ind, iquad, nquad, nker, &
        wnear, int8_0, novers, nptso, ixyzso, srcover, wover, &
        eps_gmres, niter, errs, rres, soln0, sik0)
      call helm_rpcomb_imp_solver_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, zlams, numit, &
        ifinout, rhs, nnz, row_ptr, col_ind, iquad, nquad, nker, &
        wcorr, int8_1, novers, nptso, ixyzso, srcover, wover, &
        eps_gmres, niter, errs, rres, soln1, sik1)
      err = max(zreldiff(npts, soln0, soln1), zreldiff(npts, sik0, sik1))
      call prin2('helm_rpcomb_imp_solver_guru*', err, 1)
      if (err.lt.1.0d-9) nsuccess = nsuccess + 1

!
!  13-14. helm_comb_trans: matvec and solver, two densities. The
!  operator is a difference over the two media, whose cancellation
!  leaves the rounding errors of each medium in the result
!
      zpars(1) = zk
      zpars(2) = 1.0d0
      zpars(3) = 1.0d0
      zpars(4) = 2.1d0
      zpars(5) = 1.0d0
      zpars(6) = 1.3d0
      ndz = 6
      nker = 4
      call getnearquad_helm_comb_trans(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, &
        row_ptr, col_ind, iquad, rfac0, nquad, wnear)
      call getnearcorr_helm_comb_trans(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, &
        row_ptr, col_ind, iquad, rfac0, nquad, novers, wcorr)
      call lpcomp_helm_comb_trans_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, &
        novers, nptso, ixyzso, srcover, wover, lwork, work, 2_8, sig2, &
        p20)
      call lpcomp_helm_comb_trans_corr(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wcorr, &
        novers, nptso, ixyzso, srcover, wover, lwork, work, 2_8, sig2, &
        p21)
      call report_zdiff('helm_comb_trans matvec*', 2*npts, p20, p21, &
        1.0d-10, nsuccess)
      call helm_comb_trans_solver_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, numit, sig2, nnz, &
        row_ptr, col_ind, iquad, nquad, nker, wnear, int8_0, novers, &
        nptso, ixyzso, srcover, wover, eps_gmres, niter, errs, rres, &
        p20)
      call helm_comb_trans_solver_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, numit, sig2, nnz, &
        row_ptr, col_ind, iquad, nquad, nker, wcorr, int8_1, novers, &
        nptso, ixyzso, srcover, wover, eps_gmres, niter, errs, rres, &
        p21)
      call report_zdiff('helm_comb_trans_solver_guru*', 2*npts, p20, &
        p21, 1.0d-9, nsuccess)

!
!  15-16. simple solvers: exterior problems with the field of a point
!  source inside, compared with the exact solution outside. The
!  discretization errors are below 1e-3, the same as with the add and
!  subtract evaluators, an inconsistent correction gives O(1) errors
!
      ipatch_out(1) = -1
      uvs_out(1:2,1) = 0
      zpars(1) = zk
      call h3d_slp(xyz_in, int8_3, xyz_out, int8_0, dpars, int8_1, &
        zpars, int8_0, ipars, pex)

      do i = 1,npts
        call h3d_slp(xyz_in, int8_3, srcvals(1,i), int8_0, dpars, &
          int8_1, zpars, int8_0, ipars, rhs(i))
      enddo
      zpars(2) = -ima*zk
      zpars(3) = 1.0d0
      call helm_comb_dir_solver(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, srcvals, eps, zpars, numit, ifinout, rhs, eps_gmres, &
        niter, errs, rres, soln1)
      call helm_comb_dir_eval(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, srcvals, int8_3, int8_1, xyz_out, ipatch_out, &
        uvs_out, eps, zpars, soln1, pnum)
      err = abs(pnum(1)-pex)/abs(pex)
      call prin2('helm_comb_dir_solver, error at exterior target=*', &
        err, 1)
      if (err.lt.1.0d-2) nsuccess = nsuccess + 1

      do i = 1,npts
        call h3d_sprime(xyz_in, int8_12, srcvals(1,i), int8_0, dpars, &
          int8_1, zpars, int8_0, ipars, rhs(i))
      enddo
      zpars(2) = 1.3d0
      call helm_rpcomb_neu_solver(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, zpars, numit, ifinout, rhs, &
        eps_gmres, niter, errs, rres, soln1, sik1)
      call helm_rpcomb_eval(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, srcvals, int8_3, int8_1, xyz_out, ipatch_out, &
        uvs_out, eps, zpars, soln1, sik1, pnum)
      err = abs(pnum(1)-pex)/abs(pex)
      call prin2('helm_rpcomb_neu_solver, error at exterior target=*', &
        err, 1)
      if (err.lt.1.0d-2) nsuccess = nsuccess + 1

      return
      end
!
!
!
!
!
      real *8 function zreldiff(n, a, b)
!
!  Relative max norm difference of the complex arrays a and b.
!
      implicit none
      integer *8, intent(in) :: n
      complex *16, intent(in) :: a(n), b(n)

      zreldiff = maxval(abs(a-b))/maxval(abs(a))

      return
      end
!
!
!
!
!
      subroutine report_zdiff(str, n, a, b, tol, nsuccess)
!
!  Prints the relative max norm difference of a and b, and
!  increments nsuccess if it is below tol.
!
      implicit none
      character (len=*), intent(in) :: str
      integer *8, intent(in) :: n
      complex *16, intent(in) :: a(n), b(n)
      real *8, intent(in) :: tol
      integer *8, intent(inout) :: nsuccess
      real *8 err, zreldiff

      err = zreldiff(n, a, b)
      call prin2(str, err, 1)
      if (err.lt.tol) nsuccess = nsuccess + 1

      return
      end
