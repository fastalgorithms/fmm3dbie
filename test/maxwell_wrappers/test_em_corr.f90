      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
!
!  Tests of the near quadrature corrections (getnearcorr_*, ifcorr = 1)
!  for the Maxwell wrappers on a sphere.
!
!  Applying the operators and solving with corrections must reproduce
!  the add and subtract evaluators to rounding (tests 1-6), and the
!  simple solvers, which use the corrections, must reproduce the guru
!  solvers with the add and subtract evaluators and the same
!  discretization parameters (tests 7-8).
!
!  The Muller kernels are differences over the two media of
!  hypersingular fields, whose cancellation at close source points
!  amplifies the rounding errors of each medium, so the Muller matvec
!  is compared with a larger tolerance.
!
      real *8, allocatable :: srcvals(:,:), srccoefs(:,:), targs(:,:)
      real *8, allocatable :: cms(:,:), rads(:), rad_near(:)
      integer *8, allocatable :: norders(:), ixyzs(:), iptype(:)
      integer *8, allocatable :: row_ptr(:), col_ind(:), iquad(:)
      integer *8, allocatable :: row_ptr_t(:), col_ind_t(:), iquad_t(:)
      integer *8, allocatable :: novers(:), ixyzso(:)
      integer *8, allocatable :: ipatch_id(:), ipatch_src(:)
      real *8, allocatable :: uvs_targ(:,:), uvs_src(:,:)
      real *8, allocatable :: srcover(:,:), wover(:)
      complex *16, allocatable :: wnear(:,:), wcorr(:,:)
      complex *16, allocatable :: wnear_t(:,:), wcorr_t(:,:)
      complex *16, allocatable :: wnear_e(:,:), wcorr_e(:,:)
      complex *16, allocatable :: wm(:), wmc(:)
      complex *16, allocatable :: sig(:,:), p0(:,:), p1(:,:)
      complex *16, allocatable :: sig4(:,:), pt0(:,:), pt1(:,:)
      complex *16, allocatable :: einc(:,:), hinc(:,:)
      complex *16, allocatable :: zj0(:,:), zj1(:,:), rho0(:), rho1(:)
      complex *16, allocatable :: sm(:), pm0(:), pm1(:), sm0(:), sm1(:)
      real *8, allocatable :: errs(:)
      real *8 dpars(1), work(1), c0(3)
      complex *16 zpars(2), zpm(5), ima, zk, zph
      integer *8 ipars(1)
      real *8 zreldiff
      data ima/(0.0d0,1.0d0)/

      call prini(6,13)
      nsuccess = 0
      ntests = 8

!
!  geometry: sphere of radius 1 with triangles
!
      norder = 4
      iptype0 = 1
      a = 1
      na = 3
      c0(1:3) = 0
      call get_sphere_npat_mem(a, na, c0, norder, iptype0, npatches, &
        npts)
      allocate(srcvals(12,npts), srccoefs(9,npts))
      allocate(norders(npatches), ixyzs(npatches+1), iptype(npatches))
      call get_sphere_npat(a, na, c0, norder, iptype0, npatches, npts, &
        norders, ixyzs, iptype, srccoefs, srcvals)

      eps = 1.0d-6
      zk = 1.1d0
      zpars(1) = zk
      zpars(2) = 1.0d0
!
!  omega, ep0, mu0, ep1, mu1 of the Muller equation, with
!  omega = zk so that both use the same discretization
!
      zpm(1) = zk
      zpm(2) = 1.0d0
      zpm(3) = 1.0d0
      zpm(4) = 1.44d0
      zpm(5) = 1.0d0

!
!  near field and oversampling as in the simple solvers
!
      call get_rfacs(norder, iptype0, rfac, rfac0)
      allocate(cms(3,npatches), rads(npatches), rad_near(npatches))
      call get_centroid_rads(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, cms, rads)
      do i = 1,npatches
        rad_near(i) = rads(i)*rfac
      enddo
      ndtarg = 12
      call findnearmem(cms, npatches, rad_near, ndtarg, srcvals, npts, &
        nnz)
      allocate(row_ptr(npts+1), col_ind(nnz), iquad(nnz+1))
      call findnear(cms, npatches, rad_near, ndtarg, srcvals, npts, &
        row_ptr, col_ind)
      call get_iquad_rsc(npatches, ixyzs, npts, nnz, row_ptr, col_ind, &
        iquad)
      nquad = iquad(nnz+1) - 1

      allocate(novers(npatches), ixyzso(npatches+1))
      ikerorder = 0
      call get_far_order(eps, npatches, norders, ixyzs, iptype, cms, &
        rads, npts, srccoefs, ndtarg, npts, srcvals, ikerorder, &
        zpars(1), nnz, row_ptr, col_ind, rfac, novers, ixyzso)
      nptso = ixyzso(npatches+1) - 1
      allocate(srcover(12,nptso), wover(nptso))
      call oversample_geom(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, srcvals, novers, ixyzso, nptso, srcover)
      call get_qwts(npatches, novers, ixyzso, iptype, nptso, srcover, &
        wover)

      allocate(ipatch_src(npts), uvs_src(2,npts))
      call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts, &
        ipatch_src, uvs_src)

!
!  off surface targets, at a tenth of the local patch size from
!  every fifth node, on either side of the surface, carrying the
!  tangent and normal of that node
!
      ntarg = (npts+4)/5
      allocate(targs(12,ntarg), ipatch_id(ntarg), uvs_targ(2,ntarg))
      do i = 1,ntarg
        j = 5*(i-1) + 1
        ipatch = 1
        do k = 1,npatches
          if (j.ge.ixyzs(k).and.j.lt.ixyzs(k+1)) ipatch = k
        enddo
        h = 0.1d0*rads(ipatch)*(-1)**i
        targs(1:12,i) = srcvals(1:12,j)
        targs(1:3,i) = srcvals(1:3,j) + h*srcvals(10:12,j)
        ipatch_id(i) = -1
        uvs_targ(1:2,i) = 0
      enddo
      call findnearmem(cms, npatches, rad_near, ndtarg, targs, ntarg, &
        nnz_t)
      allocate(row_ptr_t(ntarg+1), col_ind_t(nnz_t), iquad_t(nnz_t+1))
      call findnear(cms, npatches, rad_near, ndtarg, targs, ntarg, &
        row_ptr_t, col_ind_t)
      call get_iquad_rsc(npatches, ixyzs, ntarg, nnz_t, row_ptr_t, &
        col_ind_t, iquad_t)
      nquad_t = iquad_t(nnz_t+1) - 1

      allocate(sig(3,npts), p0(3,npts), p1(3,npts), sig4(4,npts))
      allocate(pt0(6,ntarg), pt1(6,ntarg))
      allocate(einc(3,npts), hinc(3,npts))
      allocate(zj0(3,npts), zj1(3,npts), rho0(npts), rho1(npts))
      allocate(sm(4*npts), pm0(4*npts), pm1(4*npts))
      allocate(sm0(4*npts), sm1(4*npts))
      do i = 1,npts
        do k = 1,3
          sig(k,i) = hkrand(0) - 0.5d0 + ima*(hkrand(0) - 0.5d0)
          sig4(k,i) = sig(k,i)
        enddo
        sig4(4,i) = hkrand(0) - 0.5d0 + ima*(hkrand(0) - 0.5d0)
        do k = 0,3
          sm(i+k*npts) = hkrand(0) - 0.5d0 + ima*(hkrand(0) - 0.5d0)
        enddo
!
!  plane wave incident field
!
        zph = exp(ima*zk*srcvals(3,i))
        einc(1,i) = zph
        einc(2,i) = 0
        einc(3,i) = 0
        hinc(1,i) = 0
        hinc(2,i) = zph
        hinc(3,i) = 0
      enddo

      iquadtype = 1
      ndd = 0
      ndi = 0
      ndz = 2
      lwork = 0
      numit = 200
      eps_gmres = 1.0d-10
      allocate(errs(numit+1))

!
!  1-2. nrccie: matvec and solver
!
      nker = 9
      allocate(wnear(9,nquad), wcorr(9,nquad))
      call getnearquad_em_nrccie_pec(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)
      call getnearcorr_em_nrccie_pec(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, zpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, novers, wcorr)
      ndim = 3
      call lpcomp_em_nrccie_pec_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, &
        novers, nptso, ixyzso, srcover, wover, lwork, work, ndim, sig, &
        p0)
      call lpcomp_em_nrccie_pec_corr(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wcorr, &
        novers, nptso, ixyzso, srcover, wover, lwork, work, ndim, sig, &
        p1)
      call report_zdiff('em_nrccie_pec matvec*', 3*npts, p0, p1, &
        1.0d-12, nsuccess)

      call em_nrccie_pec_solver_guru(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, zpars, numit, einc, hinc, nnz, &
        row_ptr, col_ind, iquad, nquad, nker, wnear, 0_8, novers, &
        nptso, ixyzso, srcover, wover, eps_gmres, niter, errs, rres, &
        zj0, rho0)
      call em_nrccie_pec_solver_guru(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, zpars, numit, einc, hinc, nnz, &
        row_ptr, col_ind, iquad, nquad, nker, wcorr, 1_8, novers, &
        nptso, ixyzso, srcover, wover, eps_gmres, niter, errs, rres, &
        zj1, rho1)
      err = max(zreldiff(3*npts, zj0, zj1), zreldiff(npts, rho0, rho1))
      call prin2('em_nrccie_pec_solver_guru*', err, 1)
      if (err.lt.1.0d-8) nsuccess = nsuccess + 1

!
!  3. nrccie operator at off surface targets
!
      allocate(wnear_t(9,nquad_t), wcorr_t(9,nquad_t))
      call getnearquad_em_nrccie_pec_targ(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz_t, row_ptr_t, &
        col_ind_t, iquad_t, rfac0, nquad_t, wnear_t)
      call getnearcorr_em_nrccie_pec_targ(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz_t, row_ptr_t, &
        col_ind_t, iquad_t, rfac0, nquad_t, novers, wcorr_t)
      call lpcomp_em_nrccie_pec_addsub_targ(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, eps, &
        ndd, dpars, ndz, zpars, ndi, ipars, nnz_t, row_ptr_t, &
        col_ind_t, iquad_t, nquad_t, nker, wnear_t, 0_8, novers, &
        nptso, ixyzso, srcover, wover, lwork, work, ndim, sig, pt0)
      call lpcomp_em_nrccie_pec_addsub_targ(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, eps, &
        ndd, dpars, ndz, zpars, ndi, ipars, nnz_t, row_ptr_t, &
        col_ind_t, iquad_t, nquad_t, nker, wcorr_t, 1_8, novers, &
        nptso, ixyzso, srcover, wover, lwork, work, ndim, sig, pt1)
      call report_zdiff('em_nrccie_pec off surface operator*', &
        3*ntarg, pt0, pt1, 1.0d-12, nsuccess)

!
!  4. nrccie post-processor: E and H at off surface targets
!
      nker = 4
      allocate(wnear_e(4,nquad_t), wcorr_e(4,nquad_t))
      call getnearquad_em_nrccie_pec_eval(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz_t, row_ptr_t, &
        col_ind_t, iquad_t, rfac0, nquad_t, wnear_e)
      call getnearcorr_em_nrccie_pec_eval(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz_t, row_ptr_t, &
        col_ind_t, iquad_t, rfac0, nquad_t, novers, wcorr_e)
      idensflag = 0
      ipotflag = 3
      ndim_s = 4
      ndim_p = 6
      call em_nrccie_eval_addsub(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, ndtarg, ntarg, targs, eps, ndd, &
        dpars, 1_8, zpars, ndi, ipars, nnz_t, row_ptr_t, col_ind_t, &
        iquad_t, nquad_t, nker, wnear_e, 0_8, novers, nptso, ixyzso, &
        srcover, wover, lwork, work, idensflag, ndim_s, sig4, ipotflag, &
        ndim_p, pt0)
      call em_nrccie_eval_addsub(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, ndtarg, ntarg, targs, eps, ndd, &
        dpars, 1_8, zpars, ndi, ipars, nnz_t, row_ptr_t, col_ind_t, &
        iquad_t, nquad_t, nker, wcorr_e, 1_8, novers, nptso, ixyzso, &
        srcover, wover, lwork, work, idensflag, ndim_s, sig4, ipotflag, &
        ndim_p, pt1)
      call report_zdiff('em_nrccie off surface E and H*', 6*ntarg, pt0, &
        pt1, 1.0d-12, nsuccess)

!
!  5-6. Muller: matvec and solver
!
      allocate(wm(16*nquad), wmc(16*nquad))
      call getnearquad_em_muller_trans(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, npts, srcvals, &
        ipatch_src, uvs_src, eps, zpm, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wm)
      call getnearcorr_em_muller_trans(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, npts, srcvals, &
        ipatch_src, uvs_src, eps, zpm, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, novers, wmc)
      call lpcomp_em_muller_trans_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, npts, srcvals, eps, &
        zpm, nnz, row_ptr, col_ind, iquad, nquad, sm, novers, nptso, &
        ixyzso, srcover, wover, pm0, wm)
      call lpcomp_em_muller_trans_corr(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, npts, srcvals, eps, &
        zpm, nnz, row_ptr, col_ind, iquad, nquad, sm, novers, nptso, &
        ixyzso, srcover, wover, pm1, wmc)
      call report_zdiff('em_muller_trans matvec*', 4*npts, pm0, pm1, &
        1.0d-9, nsuccess)

      call em_muller_trans_solver_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpm, numit, sm, nnz, &
        row_ptr, col_ind, iquad, nquad, wm, 0_8, novers, nptso, ixyzso, &
        srcover, wover, eps_gmres, niter, errs, rres, sm0)
      call em_muller_trans_solver_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpm, numit, sm, nnz, &
        row_ptr, col_ind, iquad, nquad, wmc, 1_8, novers, nptso, &
        ixyzso, srcover, wover, eps_gmres, niter, errs, rres, sm1)
      call report_zdiff('em_muller_trans_solver_guru*', 4*npts, sm0, &
        sm1, 1.0d-8, nsuccess)

!
!  7-8. simple solvers, with the same near field and oversampling as
!  above, against the guru solvers with ifcorr = 0
!
      call em_nrccie_pec_solver(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, srcvals, eps, zpars, numit, einc, hinc, eps_gmres, &
        niter, errs, rres, zj1, rho1)
      err = max(zreldiff(3*npts, zj0, zj1), zreldiff(npts, rho0, rho1))
      call prin2('em_nrccie_pec_solver*', err, 1)
      if (err.lt.1.0d-8) nsuccess = nsuccess + 1

      call em_muller_trans_solver(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, zpm, numit, sm, eps_gmres, niter, &
        errs, rres, sm1)
      call report_zdiff('em_muller_trans_solver*', 4*npts, sm0, sm1, &
        1.0d-8, nsuccess)

      call prinf('number of tests passed=*', nsuccess, 1)
      call prinf('number of tests=*', ntests, 1)

      stop
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
