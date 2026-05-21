!
!  Near-field quadrature for the two flexural wave plate boundary conditions.
!
!  Assembles both the moment (BC1) and shear force (BC2) kernels
!  simultaneously, returning wnear(1,:) for flexrepbc1 and
!  wnear(2,:) for flexrepbc2.
!
!  zpars(1:13) as in flexural_all.f90 (roots, residues, alpha, gamma, nu)
!

      subroutine getnearquad_flex_bcs(npatches, norders, &
        ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
        ipatch_id, uvs_targ, eps, zpars, iquadtype, nnz, row_ptr, &
        col_ind, iquad, rfac0, nquad, wnear)
!
!
!  This subroutine generates the near field quadrature
!  for the flexural wave boundary conditions:
!
!  BC1 (moment):    (-gamma/2 * M[G_s] + M[G_phi]) / 2
!  BC2 (shear):     (-gamma/2 * V[G_s] + V[G_phi]) / 2
!
!  where M and V are the second- and third-order boundary operators
!  from flexrepbc1 and flexrepbc2 respectively.
!
!  The quadrature is computed by the following strategy:
!  targets within a sphere of radius rfac0*rs of a patch centroid
!  are handled using adaptive integration (ggq_guru);
!  all other near-field targets use oversampled quadrature.
      implicit none
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8, intent(in) ::  srccoefs(9,npts), srcvals(12,npts)
      integer *8, intent(in) :: ndtarg, ntarg
      real *8, intent(in) :: targs(ndtarg, ntarg)
      integer *8, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      real *8, intent(in) :: eps
      complex *16, intent(in) :: zpars(13)
      integer *8, intent(in) :: iquadtype, nnz
      integer *8, intent(in) :: row_ptr(ntarg+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: nquad
      complex *16, intent(out) :: wnear(2,nquad)

!  Temporary variables      
      complex *16 zk, ima, zpars_tmp(3)
      integer *8 ipv, ndd, ndi, ndz, i
      real *8 dpars
      integer *8 ipars
      complex *16, allocatable :: wneartmp(:)

      procedure (), pointer :: fker
! replace kernels here      
      external flexrepbc1,flexrepbc2 
      
      data ima/(0.0d0,1.0d0)/
    
      ndd = 0
      ndz = 13
      ndi = 0


      allocate(wneartmp(nquad))
    
    
      if (iquadtype.eq.1) then
        ipv = 0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wneartmp(i) = 0 
        enddo
!$OMP END PARALLEL DO
        fker => flexrepbc1
        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, &
          zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, &
          rfac0, nquad, wneartmp)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wnear(1,i) = wneartmp(i)
        enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wneartmp(i) = 0
        enddo
!$OMP END PARALLEL DO
        fker => flexrepbc2
        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
          ipatch_id, uvs_targ, eps, ipv, fker, ndd, dpars, ndz, &
          zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, &
          rfac0, nquad, wneartmp)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nquad
          wnear(2,i) = wneartmp(i)
        enddo
!$OMP END PARALLEL DO
      endif


      return
      end
!


