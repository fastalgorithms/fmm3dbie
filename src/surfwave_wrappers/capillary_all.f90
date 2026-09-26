!
!  Near field quadrature for the capillary surface wave kernels
!

      subroutine getnearquad_capillary_all(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, &
        ndtarg, ntarg, targvals, ipatch_id_targ, uvs_targ_in, &
        eps, zpars, iquadtype, nnz, row_ptr, col_ind, &
        iquad, rfac0, nquad, iker, wnear)
!
!  Near field quadrature for the capillary kernels at arbitrary targets.
!  Targets on the source surface pass their patch id in ipatch_id_targ
!  and local coordinates in uvs_targ_in; other targets pass
!  ipatch_id_targ = -1.
!
!  iker selects the kernel:
!    0 = G_S, 1 = G_phi, 3 = lap G_phi, 5 = S3d G_phi,
!    8 = S' G_S, 9 = S' G_phi
!
      implicit none
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8, intent(in) ::  srccoefs(9,npts), srcvals(12,npts)
      integer *8, intent(in) :: ndtarg, ntarg
      real *8, intent(in) :: targvals(ndtarg,ntarg)
      integer *8, intent(in) :: ipatch_id_targ(ntarg)
      real *8, intent(in) :: uvs_targ_in(2,ntarg)
      real *8, intent(in) :: eps
      complex *16, intent(in) :: zpars(6)
      integer *8, intent(in) :: iquadtype, nnz
      integer *8, intent(in) :: row_ptr(ntarg+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: nquad
      integer *8, intent(in) :: iker
      complex *16, intent(out) :: wnear(nquad)

      integer *8 ipars(2)
      real *8 dpars(1)
      integer *8 ipv, i, ndi, ndd, ndz

      procedure (), pointer :: fker
      external gphihelmkern, p3d_log, gshelmkern, lapgphihelmkern
      external s3dgphihelmkern, gshelm_sp_kern, gphihelm_sp_kern

      ndz = 6
      ndd = 0
      ndi = 0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,nquad
        wnear(i) = 0
      enddo
!$OMP END PARALLEL DO

      if (iquadtype.eq.1) then
        ipv = 0

        if (iker.eq.0) then
          fker => gshelmkern
        elseif (iker.eq.1) then
          fker => gphihelmkern
        elseif (iker.eq.3) then
          fker => lapgphihelmkern
        elseif (iker.eq.5) then
          fker => s3dgphihelmkern
        elseif (iker.eq.8) then
          fker => gshelm_sp_kern
          ipv = 1
        elseif (iker.eq.9) then
          fker => gphihelm_sp_kern
          ipv = 1
        endif

        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targvals, &
          ipatch_id_targ, uvs_targ_in, eps, ipv, fker, ndd, dpars, &
          ndz, zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, &
          rfac0, nquad, wnear)

      endif

      return
      end subroutine getnearquad_capillary_all
