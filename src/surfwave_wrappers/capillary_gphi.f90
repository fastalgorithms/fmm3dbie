!
!  Near-field quadrature for the capillary wave free-surface Green's
!  function G_phi and its Laplacian Delta G_phi.
!
!  zpars(1:3) = cubic dispersion roots
!  zpars(4:6) = corresponding residues
!

      subroutine getnearquad_capillary_gphi(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, &
        eps, zpars, iquadtype, nnz, row_ptr, col_ind, &
        iquad, rfac0, nquad, wnear)
!
!  Near-field quadrature for the capillary wave free-surface Green's
!  function G_phi (gphihelmkern), evaluated on-surface (targets = sources).
!
!  Input:  zpars(6): cubic roots zpars(1:3) and residues zpars(4:6)
!  Output: wnear(nquad): near-field quadrature corrections for G_phi

      implicit none
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8, intent(in) ::  srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps
      complex *16, intent(in) :: zpars(6)
      integer *8, intent(in) :: iquadtype, nnz
      integer *8, intent(in) :: row_ptr(npts+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: nquad
      complex *16, intent(out) :: wnear(nquad)

      complex *16 zpars_tmp(3)
      integer *8 ipars(2)
      real *8 dpars(1)


      real *8, allocatable :: uvs_targ(:,:)
      integer *8, allocatable :: ipatch_id(:)


      integer *8 ipv, i, ndi, ndd, ndz

      integer *8 ndtarg, ntarg

      procedure (), pointer :: fker
      external gphihelmkern, p3d_log


      ndz=6
      ndd=0
      ndi=0
      ndtarg = 12
      ntarg = npts

      allocate(ipatch_id(npts),uvs_targ(2,npts))
!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo
!$OMP END PARALLEL DO      

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,nquad
        wnear(i) = 0
      enddo
!$OMP END PARALLEL DO

      call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts, &
        ipatch_id, uvs_targ)
      if (iquadtype.eq.1) then
        ipv=0

        fker => gphihelmkern
        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, srcvals, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, zpars, &
          ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
          wnear)

      endif

      return
      end subroutine getnearquad_capillary_gphi
!
!
!
      subroutine getnearquad_capillary_lapgphi(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, &
        eps, zpars, iquadtype, nnz, row_ptr, col_ind, &
        iquad, rfac0, nquad, wnear)
!
!  Near-field quadrature for the Laplacian of the capillary wave
!  free-surface Green's function, Delta G_phi (lapgphihelmkern),
!  evaluated on-surface (targets = sources).
!
!  Input:  zpars(6): cubic roots zpars(1:3) and residues zpars(4:6)
!  Output: wnear(nquad): near-field quadrature corrections for Delta G_phi

      implicit none
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8, intent(in) ::  srccoefs(9,npts), srcvals(12,npts)
      real *8, intent(in) :: eps
      complex *16, intent(in) :: zpars(6)
      integer *8, intent(in) :: iquadtype, nnz
      integer *8, intent(in) :: row_ptr(npts+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: nquad
      complex *16, intent(out) :: wnear(nquad)

      complex *16 zpars_tmp(3)
      integer *8 ipars(2)
      real *8 dpars(1)


      real *8, allocatable :: uvs_targ(:,:)
      integer *8, allocatable :: ipatch_id(:)


      integer *8 ipv, i, ndi, ndd, ndz

      integer *8 ndtarg, ntarg

      procedure (), pointer :: fker
      external lapgphihelmkern, p3d_log


      ndz=6
      ndd=0
      ndi=0
      ndtarg = 12
      ntarg = npts

      allocate(ipatch_id(npts),uvs_targ(2,npts))
!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo
!$OMP END PARALLEL DO      

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,nquad
        wnear(i) = 0
      enddo
!$OMP END PARALLEL DO

      call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts, &
        ipatch_id, uvs_targ)
      if (iquadtype.eq.1) then
        ipv=0

        fker => lapgphihelmkern
        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, srcvals, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, zpars, &
          ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad, &
          wnear)

      endif

      return
      end subroutine getnearquad_capillary_lapgphi
!
!
!
