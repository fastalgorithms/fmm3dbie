      subroutine test_stok_corr(nsuccess)
!
!  Tests of the near quadrature corrections (getnearcorr_*, ifcorr = 1)
!  for the Stokes wrappers on a stellarator.
!
!  Applying the operators and solving with corrections must reproduce
!  the add and subtract evaluators to rounding (tests 1-6), and the
!  simple velocity solver, which uses the corrections, must solve an
!  exterior problem with a known solution (test 7).
!
!  Output:
!    - nsuccess: number of successful tests out of 7
!
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, intent(out) :: nsuccess

      real *8, allocatable :: srcvals(:,:), srccoefs(:,:), targs(:,:)
      real *8, allocatable :: cms(:,:), rads(:), rad_near(:), wts(:)
      integer *8, allocatable :: norders(:), ixyzs(:), iptype(:)
      integer *8, allocatable :: row_ptr(:), col_ind(:), iquad(:)
      integer *8, allocatable :: row_ptr_t(:), col_ind_t(:), iquad_t(:)
      integer *8, allocatable :: novers(:), ixyzso(:)
      integer *8, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      real *8, allocatable :: srcover(:,:), wover(:)
      real *8, allocatable :: wnear(:,:), wcorr(:,:)
      real *8, allocatable :: wnear_t(:,:), wcorr_t(:,:)
      real *8, allocatable :: sig(:,:), p0(:,:), p1(:,:)
      real *8, allocatable :: pt0(:,:), pt1(:,:), rhs(:,:)
      real *8, allocatable :: soln0(:,:), soln1(:,:), errs(:)
      real *8 dpars(3), xyz_in(3), xyz_out(3), force(3)
      real *8 uex(3), unum(3,1), gij, uvs_out(2,1)
      complex *16 zpars(1)
      integer *8 ipars(3), nuv(2), iptype0, ij(2)
      integer *8 int8_0, int8_1, int8_2, int8_3, int8_12
      integer *8 ipatch_out(1)

      int8_0 = 0
      int8_1 = 1
      int8_2 = 2
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
      allocate(srcvals(12,npts), srccoefs(9,npts), wts(npts))
      allocate(norders(npatches), ixyzs(npatches+1), iptype(npatches))
      call get_stellarator_npat(nuv, norder, iptype0, npatches, npts, &
        norders, ixyzs, iptype, srccoefs, srcvals)
      call get_qwts(npatches, norders, ixyzs, iptype, npts, srcvals, &
        wts)

      xyz_in(1) = -4.501d0
      xyz_in(2) = 1.7d-3
      xyz_in(3) = 0.00001d0

      xyz_out(1) = -3.5d0
      xyz_out(2) = 3.1d0
      xyz_out(3) = 20.1d0

      eps = 1.0d-6

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
!  every fifth node, on either side of the surface, carrying the
!  normal of that node (for the traction)
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
      call findnearmem(cms, npatches, rad_near, int8_12, targs, ntarg, &
        nnz_t)
      allocate(row_ptr_t(ntarg+1), col_ind_t(nnz_t), iquad_t(nnz_t+1))
      call findnear(cms, npatches, rad_near, int8_12, targs, ntarg, &
        row_ptr_t, col_ind_t)
      call get_iquad_rsc(npatches, ixyzs, ntarg, nnz_t, row_ptr_t, &
        col_ind_t, iquad_t)
      nquad_t = iquad_t(nnz_t+1) - 1

      allocate(wnear(6,nquad), wcorr(6,nquad))
      allocate(wnear_t(6,nquad_t), wcorr_t(6,nquad_t))
      allocate(sig(3,npts), p0(3,npts), p1(3,npts))
      allocate(pt0(3,ntarg), pt1(3,ntarg))
      allocate(rhs(3,npts), soln0(3,npts), soln1(3,npts))
      do i = 1,npts
        do k = 1,3
          sig(k,i) = hkrand(0) - 0.5d0
          rhs(k,i) = hkrand(0) - 0.5d0
        enddo
      enddo

      iquadtype = 1
      ndz = 0
      nker = 6
      ndim = 3
      idensflag = 0
      ipotflag = 1
      numit = 200
      eps_gmres = 1.0d-12
      ifinout = 1
      allocate(errs(numit+1))

!
!  1-3. stok_comb_vel: matvec (with the 1's matrix term), solver, and
!  off surface evaluator
!
      dpars(1) = 1.1d0
      dpars(2) = -0.7d0
      dpars(3) = 0.3d0
      ndd = 3
      ndi = 0
      call getnearquad_stok_comb_vel(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, dpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)
      call getnearcorr_stok_comb_vel(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, dpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, novers, wcorr)
      call lpcomp_stok_comb_vel_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, &
        novers, nptso, ixyzso, srcover, wover, npts, wts, ndim, sig, p0)
      call lpcomp_stok_comb_vel_corr(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wcorr, &
        novers, nptso, ixyzso, srcover, wover, npts, wts, ndim, sig, p1)
      call report_ddiff('stok_comb_vel matvec*', 3*npts, p0, p1, &
        1.0d-12, nsuccess)

      call stok_comb_vel_solver_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, dpars, numit, ifinout, &
        rhs, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, int8_0, &
        novers, nptso, ixyzso, srcover, wover, eps_gmres, niter, errs, &
        rres, soln0)
      call stok_comb_vel_solver_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, dpars, numit, ifinout, &
        rhs, nnz, row_ptr, col_ind, iquad, nquad, nker, wcorr, int8_1, &
        novers, nptso, ixyzso, srcover, wover, eps_gmres, niter, errs, &
        rres, soln1)
      call report_ddiff('stok_comb_vel_solver_guru*', 3*npts, soln0, &
        soln1, 1.0d-9, nsuccess)

      ndd = 2
      call getnearquad_stok_comb_vel_eval(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, int8_12, ntarg, targs, &
        ipatch_id, uvs_targ, eps, dpars, iquadtype, nnz_t, row_ptr_t, &
        col_ind_t, iquad_t, rfac0, nquad_t, wnear_t)
      call getnearcorr_stok_comb_vel_eval(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, int8_12, ntarg, targs, &
        ipatch_id, uvs_targ, eps, dpars, iquadtype, nnz_t, row_ptr_t, &
        col_ind_t, iquad_t, rfac0, nquad_t, novers, wcorr_t)
      call stok_comb_vel_eval_addsub(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, int8_12, ntarg, targs, eps, ndd, &
        dpars, ndz, zpars, ndi, ipars, nnz_t, row_ptr_t, col_ind_t, &
        iquad_t, nquad_t, nker, wnear_t, int8_0, novers, nptso, ixyzso, &
        srcover, wover, int8_0, wts, idensflag, ndim, sig, ipotflag, &
        ndim, pt0)
      call stok_comb_vel_eval_addsub(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, int8_12, ntarg, targs, eps, ndd, &
        dpars, ndz, zpars, ndi, ipars, nnz_t, row_ptr_t, col_ind_t, &
        iquad_t, nquad_t, nker, wcorr_t, int8_1, novers, nptso, ixyzso, &
        srcover, wover, int8_0, wts, idensflag, ndim, sig, ipotflag, &
        ndim, pt1)
      call report_ddiff('stok_comb_vel off surface eval*', 3*ntarg, &
        pt0, pt1, 1.0d-12, nsuccess)

!
!  4-5. stok_s_mob: matvec (with the generalized 1's matrix term for a
!  single body) and solver
!
      ndd = 1
      dpars(1) = 0.5d0
      ndi = 3
      ipars(1) = 1
      ipars(2) = 1
      ipars(3) = npatches + 1
      call getnearquad_stok_s_mob(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)
      call getnearcorr_stok_s_mob(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, novers, wcorr)
      call stok_s_mob_addsub(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, &
        novers, nptso, ixyzso, srcover, wover, npts, wts, ndim, sig, p0)
      call stok_s_mob_corr(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wcorr, &
        novers, nptso, ixyzso, srcover, wover, npts, wts, ndim, sig, p1)
      call report_ddiff('stok_s_mob matvec*', 3*npts, p0, p1, 1.0d-12, &
        nsuccess)

      call stok_s_mob_solver_guru(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, ndi, ipars, numit, rhs, nnz, &
        row_ptr, col_ind, iquad, nquad, nker, wnear, int8_0, novers, &
        nptso, ixyzso, srcover, wover, eps_gmres, niter, errs, rres, &
        soln0)
      call stok_s_mob_solver_guru(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, ndi, ipars, numit, rhs, nnz, &
        row_ptr, col_ind, iquad, nquad, nker, wcorr, int8_1, novers, &
        nptso, ixyzso, srcover, wover, eps_gmres, niter, errs, rres, &
        soln1)
      call report_ddiff('stok_s_mob_solver_guru*', 3*npts, soln0, &
        soln1, 1.0d-9, nsuccess)

!
!  6. stok_sprime: off surface evaluator of the traction
!
      ndd = 0
      ndi = 0
      call getnearquad_stok_sprime(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, int8_12, ntarg, targs, ipatch_id, &
        uvs_targ, eps, iquadtype, nnz_t, row_ptr_t, col_ind_t, iquad_t, &
        rfac0, nquad_t, wnear_t)
      call getnearcorr_stok_sprime(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, int8_12, ntarg, targs, ipatch_id, &
        uvs_targ, eps, iquadtype, nnz_t, row_ptr_t, col_ind_t, iquad_t, &
        rfac0, nquad_t, novers, wcorr_t)
      call stok_sprime_eval_addsub(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, int8_12, ntarg, targs, eps, ndd, &
        dpars, ndz, zpars, ndi, ipars, nnz_t, row_ptr_t, col_ind_t, &
        iquad_t, nquad_t, nker, wnear_t, int8_0, novers, nptso, ixyzso, &
        srcover, wover, int8_0, wts, ndim, sig, pt0)
      call stok_sprime_eval_addsub(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, int8_12, ntarg, targs, eps, ndd, &
        dpars, ndz, zpars, ndi, ipars, nnz_t, row_ptr_t, col_ind_t, &
        iquad_t, nquad_t, nker, wcorr_t, int8_1, novers, nptso, ixyzso, &
        srcover, wover, int8_0, wts, ndim, sig, pt1)
      call report_ddiff('stok_sprime off surface eval*', 3*ntarg, pt0, &
        pt1, 1.0d-12, nsuccess)

!
!  7. simple velocity solver: exterior problem with the flow of a
!  stokeslet inside, compared with the exact flow outside. The
!  discretization error is below 1e-3, the same as with the add and
!  subtract evaluators, an inconsistent correction gives O(1) errors
!
      force(1) = 1.0d0
      force(2) = -0.5d0
      force(3) = 0.3d0
      do i = 1,npts
        rhs(1:3,i) = 0
        do k = 1,3
          do l = 1,3
            ij(1) = k
            ij(2) = l
            call st3d_slp(xyz_in, int8_3, srcvals(1,i), int8_0, dpars, &
              int8_0, zpars, int8_2, ij, gij)
            rhs(k,i) = rhs(k,i) + gij*force(l)
          enddo
        enddo
      enddo
      uex(1:3) = 0
      do k = 1,3
        do l = 1,3
          ij(1) = k
          ij(2) = l
          call st3d_slp(xyz_in, int8_3, xyz_out, int8_0, dpars, int8_0, &
            zpars, int8_2, ij, gij)
          uex(k) = uex(k) + gij*force(l)
        enddo
      enddo

      dpars(1) = 1.0d0
      dpars(2) = 1.0d0
      call stok_comb_vel_solver(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, srcvals, eps, dpars, numit, ifinout, rhs, eps_gmres, &
        niter, errs, rres, soln1)
      ipatch_out(1) = -1
      uvs_out(1:2,1) = 0
      call stok_comb_vel_eval(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, srcvals, int8_3, int8_1, xyz_out, ipatch_out, &
        uvs_out, eps, dpars, soln1, unum)
      err = maxval(abs(unum(1:3,1)-uex))/maxval(abs(uex))
      call prin2('stok_comb_vel_solver, error at exterior target=*', &
        err, 1)
      if (err.lt.1.0d-2) nsuccess = nsuccess + 1

      return
      end
!
!
!
!
!
      subroutine report_ddiff(str, n, a, b, tol, nsuccess)
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
