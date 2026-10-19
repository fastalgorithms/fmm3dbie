!
!  getnearquad_surf_generic.f90
!
!  Near quadrature for scalar kernels selected by registry name (see
!  kernel_registry.f90). The name is passed as ASCII codes.
!
!  If ifsub = 1, the smooth rule oversampled to novers is subtracted,
!  giving the correction to an oversampled smooth quadrature.
!
!  User callable routines:
!    - dgetnearquad_surf_generic, zgetnearquad_surf_generic
!    - dgetnearquad_surf_generic_sub, zgetnearquad_surf_generic_sub
!
!  ier: 0 success, 1 unknown kernel, 2 real/complex mismatch
!
!-----------------------------------------------------------------------
      module surf_generic_over
!
!  Setup shared by the smooth-rule subtraction routines below.
!
      implicit none
      contains

      subroutine surf_generic_over_setup(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, novers, ixyzso, nptso, &
        srcover, wover, ikey, ipmat, pmats)
!
!  Oversampled geometry and weights, and one interpolation matrix per
!  distinct (norder, nover, iptype). Matrix k is stored at pmats(ipmat(k))
!  as pt(npols,npolso).
!
      implicit none
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches), novers(npatches)
      real *8,    intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      integer *8, intent(out):: ixyzso(npatches+1), nptso
      real *8, allocatable, intent(out) :: srcover(:,:), wover(:)
      integer *8, intent(out):: ikey(npatches)
      integer *8, allocatable, intent(out) :: ipmat(:)
      real *8, allocatable, intent(out) :: pmats(:)

      integer *8, allocatable :: keys(:,:)
      real *8, allocatable :: umat(:,:)
      integer *8 :: i, k, nkey, npols, npolso, ntot

      ixyzso(1) = 1
      do i = 1, npatches
        call get_npols(iptype(i), novers(i), npolso)
        ixyzso(i+1) = ixyzso(i) + npolso
      enddo
      nptso = ixyzso(npatches+1) - 1

      allocate(srcover(12,nptso), wover(nptso))
      call oversample_geom(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, srcvals, novers, ixyzso, nptso, srcover)
      call get_qwts(npatches, novers, ixyzso, iptype, nptso, srcover, &
        wover)

      allocate(keys(3,npatches))
      nkey = 0
      do i = 1, npatches
        ikey(i) = 0
        do k = 1, nkey
          if (keys(1,k) .eq. norders(i) .and. keys(2,k) .eq. novers(i) &
            .and. keys(3,k) .eq. iptype(i)) then
            ikey(i) = k
            exit
          endif
        enddo
        if (ikey(i) .eq. 0) then
          nkey = nkey + 1
          keys(1,nkey) = norders(i)
          keys(2,nkey) = novers(i)
          keys(3,nkey) = iptype(i)
          ikey(i) = nkey
        endif
      enddo

      allocate(ipmat(nkey+1))
      ipmat(1) = 1
      do k = 1, nkey
        call get_npols(keys(3,k), keys(1,k), npols)
        call get_npols(keys(3,k), keys(2,k), npolso)
        ipmat(k+1) = ipmat(k) + npols*npolso
      enddo
      ntot = ipmat(nkey+1) - 1
      allocate(pmats(ntot))

      do k = 1, nkey
        call get_npols(keys(3,k), keys(1,k), npols)
        call get_npols(keys(3,k), keys(2,k), npolso)
        allocate(umat(npols,npols))
        umat = 0
        do i = 1, npols
          umat(i,i) = 1
        enddo
        call oversample_fun_guru(npols, keys(1,k), npols, keys(3,k), &
          umat, keys(2,k), npolso, pmats(ipmat(k)))
        deallocate(umat)
      enddo

      return
      end subroutine surf_generic_over_setup

      end module surf_generic_over

!-----------------------------------------------------------------------
      subroutine dgetnearquad_surf_generic(nchar, kername, &
        npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, &
        ndtarg, ntarg, targvals, ipatch_id, uvs_targ, eps, &
        ndd, dpars, ndz, zpars, ndi, ipars, &
        nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
        ifsub, novers, wnear, ier)
!
!  Input arguments:
!    - nchar: integer *8
!        length of kername
!    - kername: integer *8(nchar)
!        ASCII codes of the registry name of the kernel
!    - npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals:
!        source surface, standard fmm3dbie format
!        (srccoefs(9,npts), srcvals(12,npts))
!    - ndtarg, ntarg, targvals(ndtarg,ntarg): targets
!    - ipatch_id(ntarg): patch id of on-surface targets, -1 otherwise
!    - uvs_targ(2,ntarg): local coordinates of on-surface targets
!    - eps: real *8
!        requested precision
!    - ndd, dpars, ndz, zpars, ndi, ipars:
!        kernel parameters, passed unchanged to the evaluator
!    - nnz, row_ptr, col_ind, iquad: near field in row sparse
!        compressed format (see getnear)
!    - rfac0: real *8
!        radius parameter for the near field
!    - nquad: integer *8
!        length of wnear, = iquad(nnz+1)-1
!    - ifsub: integer *8
!        1 to subtract the oversampled smooth rule, 0 otherwise
!    - novers: integer *8(npatches)
!        oversampling order on each patch; only used when ifsub = 1
!
!  Output arguments:
!    - wnear: real *8(nquad)
!        near field quadrature (correction if ifsub = 1)
!    - ier: integer *8
!        error code, see top of file
!
      use kernel_registry, only: lookup_kernel
      implicit none

      integer *8, intent(in) :: nchar
      integer *8, intent(in) :: kername(max(nchar,1_8))
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8,    intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      integer *8, intent(in) :: ndtarg, ntarg
      real *8,    intent(in) :: targvals(ndtarg,ntarg)
      integer *8, intent(in) :: ipatch_id(ntarg)
      real *8,    intent(in) :: uvs_targ(2,ntarg)
      real *8,    intent(in) :: eps
      integer *8, intent(in) :: ndd, ndz, ndi
      real *8,    intent(in) :: dpars(max(ndd,1_8))
      complex *16,intent(in) :: zpars(max(ndz,1_8))
      integer *8, intent(in) :: ipars(max(ndi,1_8))
      integer *8, intent(in) :: nnz
      integer *8, intent(in) :: row_ptr(ntarg+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8,    intent(in) :: rfac0
      integer *8, intent(in) :: nquad
      integer *8, intent(in) :: ifsub
      integer *8, intent(in) :: novers(npatches)
      real *8,    intent(out):: wnear(nquad)
      integer *8, intent(out):: ier

      character(len=64) :: name
      integer *8 :: ipv, ifcomplex, ln, i
      procedure(), pointer :: fker => null()

      do i = 1, nquad
        wnear(i) = 0
      enddo

      call kername_to_fstr(kername, nchar, name, ln)
      call lookup_kernel(trim(name), fker, ipv, ifcomplex, ier)
      if (ier .ne. 0) return
      if (ifcomplex .ne. 0) then
        ier = 2
        return
      endif

      call dgetnearquad_ggq_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targvals, &
        ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear)

      if (ifsub .eq. 1) then
        call dgetnearquad_surf_generic_sub(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targvals, &
          fker, ndd, dpars, ndz, zpars, ndi, ipars, nnz, row_ptr, &
          col_ind, iquad, nquad, novers, wnear)
      endif

      return
      end subroutine dgetnearquad_surf_generic

!-----------------------------------------------------------------------
      subroutine zgetnearquad_surf_generic(nchar, kername, &
        npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, &
        ndtarg, ntarg, targvals, ipatch_id, uvs_targ, eps, &
        ndd, dpars, ndz, zpars, ndi, ipars, &
        nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
        ifsub, novers, wnear, ier)
!
!  Complex valued version of dgetnearquad_surf_generic; the arguments
!  are the same except that wnear is complex *16(nquad).
!
      use kernel_registry, only: lookup_kernel
      implicit none

      integer *8, intent(in) :: nchar
      integer *8, intent(in) :: kername(max(nchar,1_8))
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8,    intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      integer *8, intent(in) :: ndtarg, ntarg
      real *8,    intent(in) :: targvals(ndtarg,ntarg)
      integer *8, intent(in) :: ipatch_id(ntarg)
      real *8,    intent(in) :: uvs_targ(2,ntarg)
      real *8,    intent(in) :: eps
      integer *8, intent(in) :: ndd, ndz, ndi
      real *8,    intent(in) :: dpars(max(ndd,1_8))
      complex *16,intent(in) :: zpars(max(ndz,1_8))
      integer *8, intent(in) :: ipars(max(ndi,1_8))
      integer *8, intent(in) :: nnz
      integer *8, intent(in) :: row_ptr(ntarg+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8,    intent(in) :: rfac0
      integer *8, intent(in) :: nquad
      integer *8, intent(in) :: ifsub
      integer *8, intent(in) :: novers(npatches)
      complex *16,intent(out):: wnear(nquad)
      integer *8, intent(out):: ier

      character(len=64) :: name
      integer *8 :: ipv, ifcomplex, ln, i
      procedure(), pointer :: fker => null()

      do i = 1, nquad
        wnear(i) = 0
      enddo

      call kername_to_fstr(kername, nchar, name, ln)
      call lookup_kernel(trim(name), fker, ipv, ifcomplex, ier)
      if (ier .ne. 0) return
      if (ifcomplex .ne. 1) then
        ier = 2
        return
      endif

      call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targvals, &
        ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, zpars, &
        ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, wnear)

      if (ifsub .eq. 1) then
        call zgetnearquad_surf_generic_sub(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targvals, &
          fker, ndd, dpars, ndz, zpars, ndi, ipars, nnz, row_ptr, &
          col_ind, iquad, nquad, novers, wnear)
      endif

      return
      end subroutine zgetnearquad_surf_generic


!-----------------------------------------------------------------------
      subroutine dgetnearquad_surf_generic_sub(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, &
        targvals, fker, ndd, dpars, ndz, zpars, ndi, ipars, nnz, &
        row_ptr, col_ind, iquad, nquad, novers, wnear)
!
!  Subtract, in place, the oversampled smooth rule for the scalar real
!  kernel fker from the near quadrature wnear. Arguments are as in
!  dgetnearquad_surf_generic; wnear is input/output.
!
      use surf_generic_over, only: surf_generic_over_setup
      implicit none
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8,    intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      integer *8, intent(in) :: ndtarg, ntarg
      real *8,    intent(in) :: targvals(ndtarg,ntarg)
      external fker
      integer *8, intent(in) :: ndd, ndz, ndi
      real *8,    intent(in) :: dpars(*)
      complex *16,intent(in) :: zpars(*)
      integer *8, intent(in) :: ipars(*)
      integer *8, intent(in) :: nnz
      integer *8, intent(in) :: row_ptr(ntarg+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1), nquad
      integer *8, intent(in) :: novers(npatches)
      real *8, intent(inout) :: wnear(nquad)

      integer *8, allocatable :: ixyzso(:), ikey(:), ipmat(:)
      real *8, allocatable :: srcover(:,:), wover(:), pmats(:)
      integer *8 :: nptso, i, j, k, l, ipatch, npols, npolso
      integer *8 :: istarto, ip, iq, ipl
      real *8 :: rr, rs, val

      allocate(ixyzso(npatches+1), ikey(npatches))
      call surf_generic_over_setup(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, novers, ixyzso, nptso, srcover, &
        wover, ikey, ipmat, pmats)

!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(DYNAMIC) &
!$OMP PRIVATE(i,j,k,l,ipatch,npols,npolso,istarto,ip,iq,ipl,rr,rs,val)
      do i = 1, ntarg
        do j = row_ptr(i), row_ptr(i+1)-1
          ipatch  = col_ind(j)
          npols   = ixyzs(ipatch+1) - ixyzs(ipatch)
          istarto = ixyzso(ipatch)
          npolso  = ixyzso(ipatch+1) - istarto
          ip      = ipmat(ikey(ipatch))
          iq      = iquad(j)
          do l = 1, npolso
            rr = (srcover(1,istarto+l-1) - targvals(1,i))**2 + &
                 (srcover(2,istarto+l-1) - targvals(2,i))**2 + &
                 (srcover(3,istarto+l-1) - targvals(3,i))**2
            rs = max(srcover(1,istarto+l-1)**2 + &
                 srcover(2,istarto+l-1)**2 + &
                 srcover(3,istarto+l-1)**2, 1.0d0)
            if (rr .lt. 1.0d-28*rs) cycle
            call fker(srcover(1,istarto+l-1), ndtarg, targvals(1,i), &
              ndd, dpars, ndz, zpars, ndi, ipars, val)
            val = val*wover(istarto+l-1)
            ipl = ip + (l-1)*npols - 1
            do k = 1, npols
              wnear(iq+k-1) = wnear(iq+k-1) - val*pmats(ipl+k)
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

      return
      end subroutine dgetnearquad_surf_generic_sub

!-----------------------------------------------------------------------
      subroutine zgetnearquad_surf_generic_sub(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, &
        targvals, fker, ndd, dpars, ndz, zpars, ndi, ipars, nnz, &
        row_ptr, col_ind, iquad, nquad, novers, wnear)
!
!  Complex valued version of dgetnearquad_surf_generic_sub.
!
      use surf_generic_over, only: surf_generic_over_setup
      implicit none
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8,    intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      integer *8, intent(in) :: ndtarg, ntarg
      real *8,    intent(in) :: targvals(ndtarg,ntarg)
      external fker
      integer *8, intent(in) :: ndd, ndz, ndi
      real *8,    intent(in) :: dpars(*)
      complex *16,intent(in) :: zpars(*)
      integer *8, intent(in) :: ipars(*)
      integer *8, intent(in) :: nnz
      integer *8, intent(in) :: row_ptr(ntarg+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1), nquad
      integer *8, intent(in) :: novers(npatches)
      complex *16, intent(inout) :: wnear(nquad)

      integer *8, allocatable :: ixyzso(:), ikey(:), ipmat(:)
      real *8, allocatable :: srcover(:,:), wover(:), pmats(:)
      integer *8 :: nptso, i, j, k, l, ipatch, npols, npolso
      integer *8 :: istarto, ip, iq, ipl
      real *8 :: rr, rs
      complex *16 :: val

      allocate(ixyzso(npatches+1), ikey(npatches))
      call surf_generic_over_setup(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, novers, ixyzso, nptso, srcover, &
        wover, ikey, ipmat, pmats)

!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(DYNAMIC) &
!$OMP PRIVATE(i,j,k,l,ipatch,npols,npolso,istarto,ip,iq,ipl,rr,rs,val)
      do i = 1, ntarg
        do j = row_ptr(i), row_ptr(i+1)-1
          ipatch  = col_ind(j)
          npols   = ixyzs(ipatch+1) - ixyzs(ipatch)
          istarto = ixyzso(ipatch)
          npolso  = ixyzso(ipatch+1) - istarto
          ip      = ipmat(ikey(ipatch))
          iq      = iquad(j)
          do l = 1, npolso
            rr = (srcover(1,istarto+l-1) - targvals(1,i))**2 + &
                 (srcover(2,istarto+l-1) - targvals(2,i))**2 + &
                 (srcover(3,istarto+l-1) - targvals(3,i))**2
            rs = max(srcover(1,istarto+l-1)**2 + &
                 srcover(2,istarto+l-1)**2 + &
                 srcover(3,istarto+l-1)**2, 1.0d0)
            if (rr .lt. 1.0d-28*rs) cycle
            call fker(srcover(1,istarto+l-1), ndtarg, targvals(1,i), &
              ndd, dpars, ndz, zpars, ndi, ipars, val)
            val = val*wover(istarto+l-1)
            ipl = ip + (l-1)*npols - 1
            do k = 1, npols
              wnear(iq+k-1) = wnear(iq+k-1) - val*pmats(ipl+k)
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

      return
      end subroutine zgetnearquad_surf_generic_sub

!-----------------------------------------------------------------------
      subroutine rsc_to_ij(npatches, ixyzs, ntarg, nnz, row_ptr, &
        col_ind, iquad, nquad, irow, jcol)
!
!  Target and source point index of every entry of a scalar wnear
!  array, so that sparse(irow, jcol, wnear) is the near matrix.
!
      implicit none
      integer *8, intent(in) :: npatches, ixyzs(npatches+1)
      integer *8, intent(in) :: ntarg, nnz
      integer *8, intent(in) :: row_ptr(ntarg+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1), nquad
      integer *8, intent(out):: irow(nquad), jcol(nquad)

      integer *8 :: i, j, k, ipatch

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,ipatch)
      do i = 1, ntarg
        do j = row_ptr(i), row_ptr(i+1)-1
          ipatch = col_ind(j)
          do k = 0, ixyzs(ipatch+1)-ixyzs(ipatch)-1
            irow(iquad(j)+k) = i
            jcol(iquad(j)+k) = ixyzs(ipatch)+k
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

      return
      end subroutine rsc_to_ij
