!
! radpp: kernels built from a radial Green's function G(r) and its
! in-plane derivatives, with G stored as a piecewise polynomial (vpp)
! expansion of its on-axis profiles. The expansion is built in matlab
! by radpp_fit.m (see @kernel3d/radpp.m).
!
! ipars layout
!   ipars(1:5)  vpp pointers [nbin, n, ibins, icfs, nv], see vpp.f90
!   ipars(6)    nderiv, highest derivative order stored (0..3)
!   ipars(7)    ncomp, number of tensor components (1, 3, 6, 10)
!   ipars(8)    iw, start of this kernel's weights in targinfo
!
! The profiles are complex, stored as interlaced real/imaginary parts,
! sampled on the positive x axis:
!   nderiv = 0:  G
!   nderiv = 1:  + G_x
!   nderiv = 2:  + G_xx, G_yy
!   nderiv = 3:  + G_xxx, G_xyy
!
! For a radial G, the derivative tensors in direction (c,s) = x/|x|
! follow by rotation. The tensor components are ordered as in chunkie
!   [G | G_x G_y | H_xx H_xy H_yy | T_xxx T_xxy T_xyy T_yyy]
! with derivatives taken with respect to the target.
!
! The kernel is
!   K(x,y) = sum_c w_c(x) D_c G(x-y)
! with complex weights w_c stored at the target as
!   targinfo(iw+2c-2) + i*targinfo(iw+2c-1),  c = 1..ncomp
!
! radpp_tensors     - tensor components at a displacement (dx,dy,dz)
! radpp_kern        - fmm3dbie kernel, the contraction above
! getnearquad_radpp - near quadrature for nrow such kernels, with
!                     their weights stored one after the other in targs
!


subroutine radpp_tensors(dx,dy,dz,ndd,dpars,ndi,ipars,comps)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: dpars(ndd)
  integer *8 :: ipars(ndi)
  complex *16 :: comps(*)
  complex *16 :: g, gr, p, q, a, b
  real *8 :: w(12)

  nderiv = ipars(6)

  r = sqrt(dx*dx + dy*dy + dz*dz)
  c = 1.0d0
  s = 0.0d0
  if (r .gt. 0) then
     c = dx/r
     s = dy/r
  endif

  call vpp_eval(r,ndd,dpars,ndi,ipars,w)

  g = cmplx(w(1),w(2),kind=8)
  comps(1) = g

  if (nderiv .ge. 1) then
     gr = cmplx(w(3),w(4),kind=8)
     comps(2) = gr*c
     comps(3) = gr*s
  endif

  if (nderiv .ge. 2) then
     p = cmplx(w(5),w(6),kind=8)
     q = cmplx(w(7),w(8),kind=8)
     comps(4) = p*c*c + q*s*s
     comps(5) = (p-q)*c*s
     comps(6) = p*s*s + q*c*c
  endif

  if (nderiv .ge. 3) then
     a = cmplx(w(9),w(10),kind=8)
     b = cmplx(w(11),w(12),kind=8)
     comps(7)  = a*c**3 + 3*b*c*s*s
     comps(8)  = a*c*c*s + b*(s**3 - 2*c*c*s)
     comps(9)  = a*c*s*s + b*(c**3 - 2*c*s*s)
     comps(10) = a*s**3 + 3*b*s*c*c
  endif

  return
end subroutine radpp_tensors


subroutine radpp_kern(src,ndt,targ,ndd,dpars,ndz,zpars,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: src(*), targ(ndt), dpars(ndd)
  integer *8 :: ipars(ndi)
  complex *16 :: zpars(*), val
  complex *16 :: comps(10)

  ncomp = ipars(7)
  iw = ipars(8)

  call radpp_tensors(targ(1)-src(1),targ(2)-src(2),targ(3)-src(3), &
       ndd,dpars,ndi,ipars,comps)

  val = 0
  do i = 1,ncomp
     val = val + cmplx(targ(iw+2*i-2),targ(iw+2*i-1),kind=8)*comps(i)
  enddo

  return
end subroutine radpp_kern


      subroutine getnearquad_radpp(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, &
        uvs_targ, eps, nrow, ipvs, nnz, row_ptr, col_ind, iquad, &
        rfac0, nquad, ndd, dpars, ndi, ipars, wnear)
!
!  Near field quadrature for nrow radpp kernels
!
!    K_j(x,y) = sum_c w_jc(x) D_c G(x-y),  j = 1..nrow
!
!  computed with zgetnearquad_ggq_guru, one call per kernel.
!
!  Input arguments (see zgetnearquad_ggq_guru for the rest):
!    - ndtarg: integer *8
!        leading dimension of targs, 3 + 2*ncomp*nrow
!    - targs: real *8 (ndtarg,ntarg)
!        targs(1:3,i) is the location of target i, followed by the
!        complex weights of kernel 1, then kernel 2, ..., each as
!        2*ncomp interlaced real/imaginary parts
!    - nrow: integer *8
!        number of kernels
!    - ipvs: integer *8 (nrow)
!        ipv for each kernel, passed to the guru
!          * ipv = 0, weakly singular
!          * ipv = 1, principal value / hypersingular
!    - ndd, dpars: vpp expansion of the profiles
!    - ndi, ipars: ipars(1:7) as in the header of this file, ndi >= 7
!
!  Output arguments
!    - wnear: complex *16 (nrow,nquad)
!        near field quadrature, row j for kernel j
!
      implicit none
      integer *8, intent(in) :: npatches, npts
      integer *8, intent(in) :: norders(npatches), ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts), srcvals(12,npts)
      integer *8, intent(in) :: ndtarg, ntarg
      real *8, intent(in) :: targs(ndtarg,ntarg)
      integer *8, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg), eps
      integer *8, intent(in) :: nrow, ipvs(nrow)
      integer *8, intent(in) :: nnz, row_ptr(ntarg+1), col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: nquad, ndd, ndi
      real *8, intent(in) :: dpars(ndd)
      integer *8, intent(in) :: ipars(ndi)
      complex *16, intent(out) :: wnear(nrow,nquad)

      integer *8 ipars_row(8), ndi_row, ndz, ncomp, irow, i
      complex *16 zpars(1)
      complex *16, allocatable :: wtmp(:)
      procedure (), pointer :: fker
      external radpp_kern

      fker => radpp_kern
      ndz = 0
      zpars(1) = 0
      ndi_row = 8
      ipars_row(1:7) = ipars(1:7)
      ncomp = ipars(7)

      allocate(wtmp(nquad))

      do irow = 1,nrow
        ipars_row(8) = 4 + 2*ncomp*(irow-1)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i = 1,nquad
          wtmp(i) = 0
        enddo
!$OMP END PARALLEL DO

        call zgetnearquad_ggq_guru(npatches, norders, ixyzs, &
          iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
          ipatch_id, uvs_targ, eps, ipvs(irow), fker, ndd, dpars, &
          ndz, zpars, ndi_row, ipars_row, nnz, row_ptr, col_ind, &
          iquad, rfac0, nquad, wtmp)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i = 1,nquad
          wnear(irow,i) = wtmp(i)
        enddo
!$OMP END PARALLEL DO
      enddo

      deallocate(wtmp)

      return
      end subroutine getnearquad_radpp
