!
! radpp: kernels built from a radial Green's function G(r) and its
! derivatives (full 3D, or in the (x,y) plane), with G stored as a
! piecewise polynomial (vpp) expansion of its on-axis profiles. The
! expansion is built in matlab by radpp_fit.m (see @kernel3d/radpp.m).
!
! ipars layout
!   ipars(1:5)  vpp pointers [nbin, n, ibins, icfs, nv], see vpp.f90
!   ipars(6)    nderiv, highest derivative order stored (0..3)
!   ipars(7)    ncomp, number of tensor components
!   ipars(8)    dim, 3 (full 3D tensors) or 2 (in-plane tensors only)
!   ipars(9)    iw, start of this kernel's weights in targinfo
!
! The profiles are complex, stored as interlaced real/imaginary parts,
! sampled on the positive x axis:
!   nderiv = 0:  G
!   nderiv = 1:  + G_x
!   nderiv = 2:  + G_xx, G_yy
!   nderiv = 3:  + G_xxx, G_xyy
!
! For a radial G, the derivative tensors follow from the profiles by
! rotation. With derivatives taken with respect to the target, they are
! ordered (ncomp = 1, 4, 10, 20)
!   dim = 3: [G | G_x G_y G_z | H_xx H_xy H_xz H_yy H_yz H_zz |
!             T_xxx T_xxy T_xxz T_xyy T_xyz T_xzz T_yyy T_yyz T_yzz T_zzz]
! or as in chunkie (ncomp = 1, 3, 6, 10)
!   dim = 2: [G | G_x G_y | H_xx H_xy H_yy | T_xxx T_xxy T_xyy T_yyy]
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
  idim = ipars(8)

  if (idim .eq. 3) then
     call radpp_tensors3(dx,dy,dz,ndd,dpars,ndi,ipars,comps)
     return
  endif

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


subroutine radpp_tensors3(dx,dy,dz,ndd,dpars,ndi,ipars,comps)
!
! full 3D tensors, see the header
!   H_ij  = (p-q) u_i u_j + q d_ij
!   T_ijk = (a-3b) u_i u_j u_k + b (d_ij u_k + d_ik u_j + d_jk u_i)
!
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: dpars(ndd)
  integer *8 :: ipars(ndi)
  complex *16 :: comps(*)
  complex *16 :: p, q, a, b
  real *8 :: w(12), u(3), d(3,3)

  nderiv = ipars(6)

  r = sqrt(dx*dx + dy*dy + dz*dz)
  u(1) = 1.0d0
  u(2) = 0.0d0
  u(3) = 0.0d0
  if (r .gt. 0) then
     u(1) = dx/r
     u(2) = dy/r
     u(3) = dz/r
  endif

  call vpp_eval(r,ndd,dpars,ndi,ipars,w)

  comps(1) = cmplx(w(1),w(2),kind=8)
  if (nderiv .lt. 1) return

  p = cmplx(w(3),w(4),kind=8)
  do i = 1,3
     comps(1+i) = p*u(i)
  enddo
  if (nderiv .lt. 2) return

  d = 0
  do i = 1,3
     d(i,i) = 1
  enddo

  p = cmplx(w(5),w(6),kind=8)
  q = cmplx(w(7),w(8),kind=8)
  m = 4
  do i = 1,3
     do j = i,3
        m = m+1
        comps(m) = (p-q)*u(i)*u(j) + q*d(i,j)
     enddo
  enddo
  if (nderiv .lt. 3) return

  a = cmplx(w(9),w(10),kind=8)
  b = cmplx(w(11),w(12),kind=8)
  do i = 1,3
     do j = i,3
        do k = j,3
           m = m+1
           comps(m) = (a-3*b)*u(i)*u(j)*u(k) + &
                b*(d(i,j)*u(k) + d(i,k)*u(j) + d(j,k)*u(i))
        enddo
     enddo
  enddo

  return
end subroutine radpp_tensors3


subroutine radpp_kern(src,ndt,targ,ndd,dpars,ndz,zpars,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: src(*), targ(ndt), dpars(ndd)
  integer *8 :: ipars(ndi)
  complex *16 :: zpars(*), val
  complex *16 :: comps(20)

  ncomp = ipars(7)
  iw = ipars(9)

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
!  Near quadrature for the nrow rows of a radpp kernel, one call to
!  zgetnearquad_ggq_guru per row.
!
!    - targs: real *8 (ndtarg,ntarg), locations then the complex weights
!        of each row, ndtarg = 3 + 2*ncomp*nrow
!    - ipvs: integer *8 (nrow), ipv of each row (0 weakly singular,
!        1 principal value or hypersingular)
!    - wnear: complex *16 (nrow,nquad)
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

      integer *8 ipars_row(9), ndi_row, ndz, ncomp, irow, i
      complex *16 zpars(1)
      complex *16, allocatable :: wtmp(:)
      procedure (), pointer :: fker
      external radpp_kern

      fker => radpp_kern
      ndz = 0
      zpars(1) = 0
      ndi_row = 9
      ipars_row(1:8) = ipars(1:8)
      ncomp = ipars(7)

      allocate(wtmp(nquad))

      do irow = 1,nrow
        ipars_row(9) = 4 + 2*ncomp*(irow-1)

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
