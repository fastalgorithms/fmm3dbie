!
!  Tests for the name-dispatched near quadrature (src/surf_generic).
!
!  1. d/zgetnearquad_surf_generic reproduces the existing Laplace and
!     Helmholtz wrappers entry for entry.
!  2. The smooth-rule subtraction (ifsub = 1) agrees with applying the
!     oversampled smooth rule directly to an oversampled density.
!  3. Unknown names and real/complex mismatches are reported via ier.
!
      program test_surf_generic
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      real *8, allocatable :: srcvals(:,:), srccoefs(:,:), targs(:,:)
      real *8, allocatable :: cms(:,:), rads(:), rad_near(:)
      integer *8, allocatable :: norders(:), ixyzs(:), iptype(:)
      integer *8, allocatable :: ipatch_id(:), row_ptr(:), col_ind(:)
      integer *8, allocatable :: iquad(:), novers(:), ixyzso(:)
      real *8, allocatable :: uvs_targ(:,:)
      real *8, allocatable :: w1(:), w2(:), wsub(:)
      complex *16, allocatable :: z1(:), z2(:), zsub(:)
      real *8, allocatable :: sigma(:), sigover(:), srcover(:,:)
      real *8, allocatable :: wover(:)
      real *8 dpars(2), c0(3)
      complex *16 zpars(3), zk
      integer *8 ipars(1), nuv(2)
      integer *8 kname(64)
      integer *8 int8_0, int8_1
      character(len=64) nm
      external l3d_comb

      call prini(6,13)
      int8_0 = 0
      int8_1 = 1
      ntests = 0
      nsuccess = 0

      norder = 4
      nuv(1) = 4
      nuv(2) = 8
      iptype0 = 1
      call get_stellarator_npat_mem(nuv, norder, iptype0, npatches, &
        npts)
      allocate(srcvals(12,npts), srccoefs(9,npts))
      allocate(norders(npatches), ixyzs(npatches+1), iptype(npatches))
      call get_stellarator_npat(nuv, norder, iptype0, npatches, npts, &
        norders, ixyzs, iptype, srccoefs, srcvals)

      ndtarg = 12
      allocate(targs(12,npts))
      do i = 1, npts
        targs(1:12,i) = srcvals(1:12,i)
      enddo

      allocate(ipatch_id(npts), uvs_targ(2,npts))
      call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts, &
        ipatch_id, uvs_targ)

      allocate(cms(3,npatches), rads(npatches), rad_near(npatches))
      call get_centroid_rads(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, cms, rads)
      call get_rfacs(norder, iptype0, rfac, rfac0)
      do i = 1, npatches
        rad_near(i) = rads(i)*rfac
      enddo
      call findnearmem(cms, npatches, rad_near, ndtarg, targs, npts, &
        nnz)
      allocate(row_ptr(npts+1), col_ind(nnz))
      call findnear(cms, npatches, rad_near, ndtarg, targs, npts, &
        row_ptr, col_ind)
      allocate(iquad(nnz+1))
      call get_iquad_rsc(npatches, ixyzs, npts, nnz, row_ptr, col_ind, &
        iquad)
      nquad = iquad(nnz+1) - 1
      call prinf('npts=*', npts, 1)
      call prinf('nquad=*', nquad, 1)

      allocate(w1(nquad), w2(nquad), wsub(nquad))
      allocate(z1(nquad), z2(nquad), zsub(nquad))
      allocate(novers(npatches))
      do i = 1, npatches
        novers(i) = norder + 4
      enddo

      eps = 1.0d-7
      iquadtype = 1
      thresh = 1.0d-12

!
!  Laplace
!
      dpars(1) = 1.1d0
      dpars(2) = 0.7d0
      call getnearquad_lap_comb_dir_eval(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, npts, targs, &
        ipatch_id, uvs_targ, eps, dpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, w1)
      call run_d('lap3d_c', 2_8, dpars, w2)
      call check('lap3d_c', maxval(abs(w1-w2))/maxval(abs(w1)))

      dpars(1) = 1.0d0
      dpars(2) = 0.0d0
      call getnearquad_lap_comb_dir_eval(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, npts, targs, &
        ipatch_id, uvs_targ, eps, dpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, w1)
      call run_d('lap3d_s', 2_8, dpars, w2)
      call check('lap3d_s', maxval(abs(w1-w2))/maxval(abs(w1)))

      dpars(1) = 0.6d0
      dpars(2) = 1.3d0
      call getnearquad_lap_comb_cprime_eval(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, npts, targs, &
        ipatch_id, uvs_targ, eps, dpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, w1)
      call run_d('lap3d_cp', 2_8, dpars, w2)
      call check('lap3d_cp', maxval(abs(w1-w2))/maxval(abs(w1)))

      dpars(1) = 1.0d0
      dpars(2) = 0.0d0
      call getnearquad_lap_comb_cprime_eval(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, npts, targs, &
        ipatch_id, uvs_targ, eps, dpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, w1)
      call run_d('lap3d_cp_pv', 2_8, dpars, w2)
      call check('lap3d_cp_pv', maxval(abs(w1-w2))/maxval(abs(w1)))

!
!  Helmholtz
!
      zk = dcmplx(1.3d0, 0.1d0)
      zpars(1) = zk
      zpars(2) = dcmplx(1.0d0, 0.2d0)
      zpars(3) = dcmplx(-0.4d0, 1.0d0)
      call getnearquad_helm_comb_dir_eval(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, npts, targs, &
        ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, z1)
      call run_z('helm3d_c', 3_8, zpars, z2)
      call check('helm3d_c', maxval(abs(z1-z2))/maxval(abs(z1)))

      call getnearquad_helm_comb_cprime_eval(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, npts, targs, &
        ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, z1)
      call run_z('helm3d_cp', 3_8, zpars, z2)
      call check('helm3d_cp', maxval(abs(z1-z2))/maxval(abs(z1)))

!
!  Subtraction: for a smooth density sigma,
!    sum_k wsub(i,k) sigma_k  =  - sum_l K(t_i,s_l) wover_l sigover_l
!  over the near patches of target i, where sigover is sigma
!  oversampled with oversample_fun_surf.
!
      allocate(ixyzso(npatches+1))
      ixyzso(1) = 1
      do i = 1, npatches
        call get_npols(iptype(i), novers(i), npo)
        ixyzso(i+1) = ixyzso(i) + npo
      enddo
      nptso = ixyzso(npatches+1) - 1
      allocate(srcover(12,nptso), wover(nptso), sigma(npts))
      allocate(sigover(nptso))
      call oversample_geom(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, srcvals, novers, ixyzso, nptso, srcover)
      call get_qwts(npatches, novers, ixyzso, iptype, nptso, srcover, &
        wover)
      do i = 1, npts
        sigma(i) = cos(srcvals(1,i)) + srcvals(2,i)*srcvals(3,i)
      enddo
      call oversample_fun_surf(int8_1, npatches, norders, ixyzs, &
        iptype, npts, sigma, novers, ixyzso, nptso, sigover)

      dpars(1) = 1.1d0
      dpars(2) = 0.7d0
      call run_d('lap3d_c', 2_8, dpars, w2)
      call str2codes('lap3d_c', kname, nc)
      call dgetnearquad_surf_generic(nc, kname, npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, npts, targs, &
        ipatch_id, uvs_targ, eps, 2_8, dpars, int8_0, zpars, &
        int8_0, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
        int8_1, novers, wsub, ier)
      errmax = 0
      rmax = 0
      do i = 1, npts
        a1 = 0
        a2 = 0
        do j = row_ptr(i), row_ptr(i+1)-1
          ip = col_ind(j)
          do k = 1, ixyzs(ip+1)-ixyzs(ip)
            a1 = a1 + (wsub(iquad(j)+k-1) - w2(iquad(j)+k-1))* &
              sigma(ixyzs(ip)+k-1)
          enddo
          do l = ixyzso(ip), ixyzso(ip+1)-1
            rr = sqrt(sum((srcover(1:3,l)-targs(1:3,i))**2))
            if (rr .lt. 1.0d-14*max(sqrt(sum(srcover(1:3,l)**2)), &
              1.0d0)) cycle
            call l3d_comb(srcover(1,l), ndtarg, targs(1,i), 2_8, &
              dpars, int8_0, zpars, int8_0, ipars, val)
            a2 = a2 - val*wover(l)*sigover(l)
          enddo
        enddo
        errmax = max(errmax, abs(a1-a2))
        rmax = max(rmax, abs(a2))
      enddo
      call check('lap3d_c subtraction', errmax/rmax)

!
!  Error reporting
!
      call str2codes('no_such_kernel', kname, nc)
      call dgetnearquad_surf_generic(nc, kname, npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, npts, targs, &
        ipatch_id, uvs_targ, eps, 2_8, dpars, int8_0, zpars, &
        int8_0, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
        int8_0, novers, w2, ier)
      ntests = ntests + 1
      if (ier .eq. 1) nsuccess = nsuccess + 1
      call prinf('unknown name ier (expect 1)=*', ier, 1)

      call str2codes('helm3d_s', kname, nc)
      call dgetnearquad_surf_generic(nc, kname, npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, npts, targs, &
        ipatch_id, uvs_targ, eps, 2_8, dpars, 3_8, zpars, &
        int8_0, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
        int8_0, novers, w2, ier)
      ntests = ntests + 1
      if (ier .eq. 2) nsuccess = nsuccess + 1
      call prinf('real call on complex kernel ier (expect 2)=*', ier, 1)

      write(*,'(a,i3,a,i3,a)') 'surf_generic: ', nsuccess, ' of ', &
        ntests, ' tests passed'
      open(unit=33, file='../../print_testres.txt', access='append')
      write(33,'(a,i3,a,i3,a)') 'Successfully completed ', nsuccess, &
        ' out of ', ntests, ' in surf_generic testing suite'
      close(33)

      contains

        subroutine run_d(name, nd, dp, w)
          character(len=*), intent(in) :: name
          integer *8, intent(in) :: nd
          real *8, intent(in) :: dp(nd)
          real *8, intent(out) :: w(nquad)
          integer *8 :: kn(64), ncc, ierr
          call str2codes(name, kn, ncc)
          call dgetnearquad_surf_generic(ncc, kn, npatches, norders, &
            ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, npts, &
            targs, ipatch_id, uvs_targ, eps, nd, dp, int8_0, zpars, &
            int8_0, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
            int8_0, novers, w, ierr)
          if (ierr .ne. 0) call prinf('run_d ier=*', ierr, 1)
        end subroutine run_d

        subroutine run_z(name, nz, zp, w)
          character(len=*), intent(in) :: name
          integer *8, intent(in) :: nz
          complex *16, intent(in) :: zp(nz)
          complex *16, intent(out) :: w(nquad)
          integer *8 :: kn(64), ncc, ierr
          call str2codes(name, kn, ncc)
          call zgetnearquad_surf_generic(ncc, kn, npatches, norders, &
            ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, npts, &
            targs, ipatch_id, uvs_targ, eps, int8_0, dpars, nz, zp, &
            int8_0, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
            int8_0, novers, w, ierr)
          if (ierr .ne. 0) call prinf('run_z ier=*', ierr, 1)
        end subroutine run_z

        subroutine check(name, err)
          character(len=*), intent(in) :: name
          real *8, intent(in) :: err
          ntests = ntests + 1
          if (err .lt. 1.0d-12) nsuccess = nsuccess + 1
          write(*,'(a,a,es10.2)') name, ': rel err = ', err
        end subroutine check

      end program test_surf_generic

      subroutine str2codes(s, codes, n)
      implicit none
      character(len=*), intent(in) :: s
      integer *8, intent(out) :: codes(64), n
      integer *8 :: i
      n = len_trim(s)
      do i = 1, n
        codes(i) = iachar(s(i:i))
      enddo
      return
      end subroutine str2codes
