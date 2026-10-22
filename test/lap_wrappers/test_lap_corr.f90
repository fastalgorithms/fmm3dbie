      subroutine test_lap_corr(nsuccess)
!
!  Tests of the near quadrature corrections (getnearcorr_*, ifcorr = 1)
!  for the Laplace wrappers on a stellarator.
!
!  Applying the operators with corrections must reproduce the add and
!  subtract evaluators to rounding (tests 1-7), and the simple solvers,
!  which use the corrections, must solve exterior problems with a
!  known solution (tests 8-9).
!
!  Output:
!    - nsuccess: number of successful tests out of 9
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
      real *8, allocatable :: wnear(:), wcorr(:)
      real *8, allocatable :: wnear_t(:), wcorr_t(:)
      real *8, allocatable :: sigma(:), pot0(:), pot1(:)
      real *8, allocatable :: pott0(:), pott1(:), rhs(:)
      real *8, allocatable :: soln0(:), soln1(:), errs(:)
      real *8 dpars(2), work(1), xyz_in(3), xyz_out(3), pex, pnum(1)
      complex *16 zpars(1)
      integer *8 ipars(1), nuv(2), iptype0
      integer *8 int8_0, int8_3, int8_12
      integer *8 ipatch_out(1)
      real *8 uvs_out(2,1)

      int8_0 = 0
      int8_3 = 3
      int8_12 = 12
      nsuccess = 0
      done = 1
      pi = atan(done)*4

!
!  geometry
!
      norder = 5
      nuv(1) = 6
      nuv(2) = 18
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

      eps = 0.51d-8

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

      allocate(wnear(nquad), wcorr(nquad))
      allocate(wnear_t(nquad_t), wcorr_t(nquad_t))
      allocate(sigma(npts), pot0(npts), pot1(npts))
      allocate(pott0(ntarg), pott1(ntarg))
      do i = 1,npts
        sigma(i) = hkrand(0) - 0.5d0
      enddo

      iquadtype = 1
      ndz = 0
      ndi = 0
      nker = 1
      lwork = 0
      ndim = 1
      idensflag = 0
      ipotflag = 1
      ndim_p = 1

!
!  1. lap_comb_dir: matvec
!
      ndd = 2
      dpars(1) = 1.1d0
      dpars(2) = -0.7d0
      call getnearquad_lap_comb_dir(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, dpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)
      call getnearcorr_lap_comb_dir(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, dpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, novers, wcorr)
      call lpcomp_lap_comb_dir_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, &
        novers, nptso, ixyzso, srcover, wover, lwork, work, ndim, &
        sigma, pot0)
      call lpcomp_lap_comb_dir_corr(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wcorr, &
        novers, nptso, ixyzso, srcover, wover, lwork, work, ndim, &
        sigma, pot1)
      call report_diff('lap_comb_dir matvec*', npts, pot0, pot1, &
        1.0d-12, nsuccess)

!
!  2. lap_comb_dir: off surface evaluator
!
      call getnearquad_lap_comb_dir_eval(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, int8_3, ntarg, targs, &
        ipatch_id, uvs_targ, eps, dpars, iquadtype, nnz_t, row_ptr_t, &
        col_ind_t, iquad_t, rfac0, nquad_t, wnear_t)
      call getnearcorr_lap_comb_dir_eval(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, int8_3, ntarg, targs, &
        ipatch_id, uvs_targ, eps, dpars, iquadtype, nnz_t, row_ptr_t, &
        col_ind_t, iquad_t, rfac0, nquad_t, novers, wcorr_t)
      call lap_comb_dir_eval_addsub(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, int8_3, ntarg, targs, eps, ndd, dpars, &
        ndz, zpars, ndi, ipars, nnz_t, row_ptr_t, col_ind_t, iquad_t, &
        nquad_t, nker, wnear_t, int8_0, novers, nptso, ixyzso, srcover, &
        wover, lwork, work, idensflag, ndim, sigma, ipotflag, ndim_p, &
        pott0)
      call lap_comb_dir_eval_addsub(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, int8_3, ntarg, targs, eps, ndd, dpars, &
        ndz, zpars, ndi, ipars, nnz_t, row_ptr_t, col_ind_t, iquad_t, &
        nquad_t, nker, wcorr_t, 1_8, novers, nptso, ixyzso, srcover, &
        wover, lwork, work, idensflag, ndim, sigma, ipotflag, ndim_p, &
        pott1)
      call report_diff('lap_comb_dir off surface eval*', ntarg, pott0, &
        pott1, 1.0d-12, nsuccess)

!
!  3. lap_comb_cprime: on surface evaluator
!
      dpars(1) = 1.0d0
      dpars(2) = 0.5d0
      call getnearquad_lap_comb_cprime(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, dpars, iquadtype, nnz, &
        row_ptr, col_ind, iquad, rfac0, nquad, wnear)
      call getnearcorr_lap_comb_cprime(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, dpars, iquadtype, nnz, &
        row_ptr, col_ind, iquad, rfac0, nquad, novers, wcorr)
      call lap_comb_cprime_eval_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, int8_12, npts, srcvals, eps, &
        ndd, dpars, ndz, zpars, ndi, ipars, nnz, row_ptr, col_ind, &
        iquad, nquad, nker, wnear, int8_0, novers, nptso, ixyzso, &
        srcover, wover, lwork, work, idensflag, ndim, sigma, ipotflag, &
        ndim_p, pot0)
      call lap_comb_cprime_eval_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, int8_12, npts, srcvals, eps, &
        ndd, dpars, ndz, zpars, ndi, ipars, nnz, row_ptr, col_ind, &
        iquad, nquad, nker, wcorr, 1_8, novers, nptso, ixyzso, &
        srcover, wover, lwork, work, idensflag, ndim, sigma, ipotflag, &
        ndim_p, pot1)
      call report_diff('lap_comb_cprime on surface eval*', npts, pot0, &
        pot1, 1.0d-12, nsuccess)

!
!  4. lap_s_neu: matvec, with the rank one term
!
      ndd = 1
      dpars(1) = 0.3d0
      call getnearquad_lap_s_neu(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)
      call getnearcorr_lap_s_neu(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, novers, wcorr)
      call lpcomp_lap_s_neu_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, &
        novers, nptso, ixyzso, srcover, wover, lwork, work, ndim, &
        sigma, pot0)
      call lpcomp_lap_s_neu_corr(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wcorr, &
        novers, nptso, ixyzso, srcover, wover, lwork, work, ndim, &
        sigma, pot1)
      call report_diff('lap_s_neu matvec*', npts, pot0, pot1, 1.0d-12, &
        nsuccess)

!
!  5. lap_s_neu: off surface evaluator
!
      call getnearquad_lap_s_neu_eval(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, int8_3, ntarg, targs, &
        ipatch_id, uvs_targ, eps, iquadtype, nnz_t, row_ptr_t, &
        col_ind_t, iquad_t, rfac0, nquad_t, wnear_t)
      call getnearcorr_lap_s_neu_eval(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, int8_3, ntarg, targs, &
        ipatch_id, uvs_targ, eps, iquadtype, nnz_t, row_ptr_t, &
        col_ind_t, iquad_t, rfac0, nquad_t, novers, wcorr_t)
      call lap_s_neu_eval_addsub(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, int8_3, ntarg, targs, eps, ndd, dpars, &
        ndz, zpars, ndi, ipars, nnz_t, row_ptr_t, col_ind_t, iquad_t, &
        nquad_t, nker, wnear_t, int8_0, novers, nptso, ixyzso, srcover, &
        wover, lwork, work, idensflag, ndim, sigma, ipotflag, ndim_p, &
        pott0)
      call lap_s_neu_eval_addsub(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, int8_3, ntarg, targs, eps, ndd, dpars, &
        ndz, zpars, ndi, ipars, nnz_t, row_ptr_t, col_ind_t, iquad_t, &
        nquad_t, nker, wcorr_t, 1_8, novers, nptso, ixyzso, srcover, &
        wover, lwork, work, idensflag, ndim, sigma, ipotflag, ndim_p, &
        pott1)
      call report_diff('lap_s_neu off surface eval*', ntarg, pott0, &
        pott1, 1.0d-12, nsuccess)

!
!  6-7. solver gurus: same solution with and without corrections
!
      numit = 200
      eps_gmres = 1.0d-12
      ifinout = 1
      allocate(rhs(npts), soln0(npts), soln1(npts), errs(numit+1))
      do i = 1,npts
        call l3d_slp(xyz_in, int8_3, srcvals(1,i), int8_0, dpars, &
          int8_0, zpars, int8_0, ipars, rhs(i))
      enddo

      dpars(1) = 1.0d0
      dpars(2) = 1.0d0
      call getnearquad_lap_comb_dir(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, dpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)
      call getnearcorr_lap_comb_dir(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, dpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, novers, wcorr)
      call lap_comb_dir_solver_guru(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, dpars, numit, ifinout, rhs, nnz, &
        row_ptr, col_ind, iquad, nquad, nker, wnear, int8_0, novers, &
        nptso, ixyzso, srcover, wover, eps_gmres, niter0, errs, rres, &
        soln0)
      call lap_comb_dir_solver_guru(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, dpars, numit, ifinout, rhs, nnz, &
        row_ptr, col_ind, iquad, nquad, nker, wcorr, 1_8, novers, &
        nptso, ixyzso, srcover, wover, eps_gmres, niter1, errs, rres, &
        soln1)
      call report_diff('lap_comb_dir_solver_guru*', npts, soln0, soln1, &
        1.0d-9, nsuccess)

      do i = 1,npts
        call l3d_sprime(xyz_in, int8_12, srcvals(1,i), int8_0, dpars, &
          int8_0, zpars, int8_0, ipars, rhs(i))
      enddo
      call getnearquad_lap_s_neu(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)
      call getnearcorr_lap_s_neu(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, novers, wcorr)
      call lap_s_neu_solver_guru(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, numit, ifinout, rhs, nnz, &
        row_ptr, col_ind, iquad, nquad, nker, wnear, int8_0, novers, &
        nptso, ixyzso, srcover, wover, eps_gmres, niter0, errs, rres, &
        soln0)
      call lap_s_neu_solver_guru(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, numit, ifinout, rhs, nnz, &
        row_ptr, col_ind, iquad, nquad, nker, wcorr, 1_8, novers, &
        nptso, ixyzso, srcover, wover, eps_gmres, niter1, errs, rres, &
        soln1)
      call report_diff('lap_s_neu_solver_guru*', npts, soln0, soln1, &
        1.0d-9, nsuccess)

!
!  8-9. simple solvers: exterior problems with the field of a point
!  source inside, compared with the exact solution outside. The
!  discretization error is about 1e-6 (Dirichlet) and 2e-5 (Neumann),
!  same as with the add and subtract evaluators, an inconsistent
!  correction gives O(1) errors
!
      ipatch_out(1) = -1
      uvs_out(1:2,1) = 0
      call l3d_slp(xyz_in, int8_3, xyz_out, int8_0, dpars, int8_0, &
        zpars, int8_0, ipars, pex)

      do i = 1,npts
        call l3d_slp(xyz_in, int8_3, srcvals(1,i), int8_0, dpars, &
          int8_0, zpars, int8_0, ipars, rhs(i))
      enddo
      dpars(1) = 1.0d0
      dpars(2) = 1.0d0
      call lap_comb_dir_solver(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, srcvals, eps, dpars, numit, ifinout, rhs, eps_gmres, &
        niter1, errs, rres, soln1)
      call lap_comb_dir_eval(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, srcvals, int8_3, 1_8, xyz_out, ipatch_out, uvs_out, &
        eps, dpars, soln1, pnum)
      err = abs(pnum(1)-pex)/abs(pex)
      call prin2('lap_comb_dir_solver, error at exterior target=*', &
        err, 1)
      if (err.lt.1.0d-4) nsuccess = nsuccess + 1

      do i = 1,npts
        call l3d_sprime(xyz_in, int8_12, srcvals(1,i), int8_0, dpars, &
          int8_0, zpars, int8_0, ipars, rhs(i))
      enddo
      call lap_s_neu_solver(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, srcvals, eps, numit, ifinout, rhs, eps_gmres, niter1, &
        errs, rres, soln1)
      call lap_s_neu_eval(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, srcvals, int8_3, 1_8, xyz_out, ipatch_out, uvs_out, &
        eps, soln1, pnum)
      err = abs(pnum(1)-pex)/abs(pex)
      call prin2('lap_s_neu_solver, error at exterior target=*', err, 1)
      if (err.lt.1.0d-4) nsuccess = nsuccess + 1

      return
      end
!
!
!
!
!
      subroutine report_diff(str, n, a, b, tol, nsuccess)
!
!  Prints the relative max norm difference of a and b, and
!  increments nsuccess if it is below tol.
!
      implicit none
      character (len=*), intent(in) :: str
      integer *8, intent(in) :: n
      real *8, intent(in) :: a(n), b(n), tol
      integer *8, intent(inout) :: nsuccess
      real *8 err

      err = maxval(abs(a-b))/maxval(abs(a))
      call prin2(str, err, 1)
      if (err.lt.tol) nsuccess = nsuccess + 1

      return
      end
