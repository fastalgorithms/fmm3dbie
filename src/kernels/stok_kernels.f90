


!
! the following routines rely on the srcinfo and targinfo arrays
! containing the following information, standardized in the following
! way:
!
!     *info(1:3) = xyz
!     *info(4:6) = tanget vector 1
!     *info(7:9) = tangent vector 2
!     *info(10:12) = normal vector
!
!

!      
!     We take the following conventions for the Stokes kernels
!
!     For a source y and target x, let r_i = x_i-y_i
!     and let r = sqrt(r_1^2 + r_2^2 + r_3^2)
!
!     The Stokeslet, G_{ij}, and its associated pressure tensor, P_j,
!     (without the 1/4pi scaling) are
!
!     G_{ij}(x,y) = (r_i r_j)/(2r^3) + delta_{ij}/(2r)
!     P_j(x,y) = -r_j/r^3
!
!     The (Type I) stresslet, T_{ijk}, and its associated pressure
!     tensor, PI_{jk}, (without the 1/4pi scaling) are
!     
!     T_{ijk}(x,y) = -3 r_i r_j r_k/ r^5
!     PI_{jk} = -2 delta_{jk} + 6 r_j r_k/r^5      
!
!



subroutine st3d_slp_vec(nd,src,ndt,targ,ndd,dpars,ndz,zk,ndi, &
     ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: src(*), targ(ndt),dpars(ndd)
  integer ipars(ndi)
  
  complex *16 :: zk
  real *8 :: val(nd)

!f2py intent(in) nd,src,ndt,targ,ndd,dpars,ndz,zk,ndi,ipars
!f2py intent(out) val

  !
  ! returns the Stokes single layer potential kernel
  !
  ! G_ij = 0.5*(targ_i-src_i)(targ_j-src_j)/|src-targ|^3 +
  !               0.5*delta_ij/|src-targ|
  !
  ! does *not* include the 1/4pi scaling to be consistent
  ! with the FMM. The output is given ordered by standard
  ! linear indexing, ij -> i+(j-1)*3
  
  dx=targ(1)-src(1)
  dy=targ(2)-src(2)
  dz=targ(3)-src(3)

  r=sqrt(dx**2+dy**2+dz**2)

  rinv = 1.0d0/r
  rinv3 = 0.5d0*rinv**3
  rinv = 0.5d0*rinv
  dxdy = dx*dy*rinv3
  dxdz = dx*dz*rinv3
  dydz = dy*dz*rinv3
  
  val(1) = rinv + dx*dx*rinv3
  val(2) = dxdy
  val(3) = dxdz
  val(4) = dxdy
  val(5) = rinv + dy*dy*rinv3
  val(6) = dydz
  val(7) = dxdz
  val(8) = dydz
  val(9) = rinv + dz*dz*rinv3


  return
end subroutine st3d_slp_vec


subroutine st3d_slp(src,ndt,targ,ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: src(*), targ(ndt),dpars(ndd)
  integer ipars(ndi)

  complex *16 :: zk
  real *8 :: val, dr(3)
!f2py intent(in) src,ndt,targ,ndd,dpars,ndz,zk,ndi,ipars
!f2py intent(out) val

  !
  ! returns one entry of the Stokes single layer potential
  ! kernel
  !
  !
  ! G_ij = 0.5*(targ_i-src_i)(targ_j-src_j)/|src-targ|^3 +
  !               0.5*delta_ij/|src-targ|
  !
  ! does *not* include the 1/4pi scaling to be consistent with
  ! the FMM. Returns val=G_ij with i = ipars(1), j = ipars(2)
  

  i = ipars(1)
  j = ipars(2)
  
  dr(1)=targ(1)-src(1)
  dr(2)=targ(2)-src(2)
  dr(3)=targ(3)-src(3)

  dxi = dr(i)
  dxj = dr(j)
  
  r=sqrt(dx**2+dy**2+dz**2)
  rinv = 1.0d0/r
  rinv3 = 0.5d0*rinv**3
  rinv = 0.5d0*rinv
  
  diag = 0.0d0
  if (i .eq. j) diag = rinv

  val = diag + dxi*dxj*rinv3

  return
end subroutine st3d_slp





subroutine st3d_dlp_vec(nd,srcinfo,ndt,targ,ndd,dpars,ndz,zk,ndi, &
     ipars,val)
!f2py intent(in) nd,src,ndt,targ,ndd,dpars,ndz,zk,ndi,ipars
!f2py intent(out) val
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(12), targ(ndt),dpars(ndd)
  integer ipars(ndi)
  complex *16 :: zk
  real *8 :: val(nd)

  real *8 :: src(3), srcnorm(3)
  !
  ! returns the Stokes double layer kernel
  !
  ! D_ij = 3*(sum_k (targ_k-src_k) n_k)
  !               *(targ_j-src_j)(targ_i-src_i)/r^5
  !
  ! where r = |src-targ| and n is the normal vector at 
  ! the source. The output is given ordered by standard
  ! linear indexing, ij -> i+(j-1)*3

  src(1)=srcinfo(1)
  src(2)=srcinfo(2)
  src(3)=srcinfo(3)
  srcnorm(1)=srcinfo(10)
  srcnorm(2)=srcinfo(11)
  srcnorm(3)=srcinfo(12)

!  call prin2('srcinfo=*',srcinfo,12)

  dx=targ(1)-src(1)
  dy=targ(2)-src(2)
  dz=targ(3)-src(3)

  dprod = dx*srcnorm(1) + dy*srcnorm(2) + dz*srcnorm(3)

  r=sqrt(dx**2+dy**2+dz**2)
  rinv5 = 3.0d0*dprod*(1.0d0/r)**5

  dxdy = dx*dy*rinv5
  dxdz = dx*dz*rinv5
  dydz = dy*dz*rinv5

  val(1) = dx*dx*rinv5
  val(2) = dxdy
  val(3) = dxdz
  val(4) = dxdy
  val(5) = dy*dy*rinv5
  val(6) = dydz
  val(7) = dxdz
  val(8) = dydz
  val(9) = dz*dz*rinv5

  return
end subroutine st3d_dlp_vec


subroutine st3d_dlp(srcinfo,ndt,targ,ndd,dpars,ndz, &
     zk,ndi,ipars,val)
!f2py intent(in) src,ndt,targ,ndd,dpars,ndz,zk,ndi,ipars
!f2py intent(out) val
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(12), targ(ndt),dpars(ndd)
  integer ipars(ndi)
  complex *16 :: zk
  real *8 :: val

  real *8 :: src(3), srcnorm(3), dr(3)
  !
  ! returns one entry of the Stokes double layer kernel
  !
  ! D_ij = -T_jik n_k = 3*(sum_k (targ_k-src_k) n_k)
  !               *(targ_j-src_j)(targ_i-src_i)/r^5
  !
  ! where r = |src-targ| and n is the normal vector at
  ! the source. Returns val=D_ij with i = ipars(1), j = ipars(2)

  i = ipars(1)
  j = ipars(2)
  
  src(1)=srcinfo(1)
  src(2)=srcinfo(2)
  src(3)=srcinfo(3)
  srcnorm(1)=srcinfo(10)
  srcnorm(2)=srcinfo(11)
  srcnorm(3)=srcinfo(12)

!  call prin2('srcinfo=*',srcinfo,12)

  dx=targ(1)-src(1)
  dy=targ(2)-src(2)
  dz=targ(3)-src(3)

  dprod = dx*srcnorm(1) + dy*srcnorm(2) + dz*srcnorm(3)

  r=sqrt(dx**2+dy**2+dz**2)
  rinv5 = 3.0d0*dprod*(1.0d0/r)**5

  dr(1) = dx
  dr(2) = dy
  dr(3) = dz

  dxi = dr(i)
  dxj = dr(j)
  val = dxi*dxj*rinv5

  return
end subroutine st3d_dlp


subroutine st3d_comb_vec(nd,srcinfo,ndt,targ,ndd,dpars,ndz,zk,ndi, &
     ipars,val)
!f2py intent(in) nd,src,ndt,targ,ndd,dpars,ndz,zk,ndi,ipars
!f2py intent(out) val
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(12), targ(ndt),dpars(ndd)
  integer ipars(ndi)
  complex *16 :: zk
  real *8 :: val(nd)

  real *8 :: src(3), srcnorm(3), alpha, beta
  !
  ! returns the Stokes combined layer kernel
  !
  ! K_ij = alpha G_ij + beta D_ij
  !
  ! Where alpha = dpars(1), beta = dpars(2), G_ij is the
  ! Stokeslet and :
  !
  ! D_ij = 3*(sum_k (targ_k-src_k) n_k)
  !               *(targ_j-src_j)(targ_i-src_i)/r^5
  !
  ! where r = |src-targ| and n is the normal vector at 
  ! the source. The output is given ordered by standard
  ! linear indexing, ij -> i+(j-1)*3

  alpha = dpars(1)
  beta = dpars(2)
  
  src(1)=srcinfo(1)
  src(2)=srcinfo(2)
  src(3)=srcinfo(3)
  srcnorm(1)=srcinfo(10)
  srcnorm(2)=srcinfo(11)
  srcnorm(3)=srcinfo(12)

!  call prin2('srcinfo=*',srcinfo,12)

  dx=targ(1)-src(1)
  dy=targ(2)-src(2)
  dz=targ(3)-src(3)

  dprod = dx*srcnorm(1) + dy*srcnorm(2) + dz*srcnorm(3)

  r=sqrt(dx**2+dy**2+dz**2)

  rinv = 1.0d0/r
  rinv3 = 0.5d0*rinv**3*alpha
  rinv5 = 3.0d0*dprod*rinv**5*beta
  rinv = 0.5d0*rinv*alpha

  dcomb = rinv3+rinv5
  
  dxdy = dx*dy*dcomb
  dxdz = dx*dz*dcomb
  dydz = dy*dz*dcomb

  val(1) = rinv + dx*dx*dcomb
  val(2) = dxdy
  val(3) = dxdz
  val(4) = dxdy
  val(5) = rinv + dy*dy*dcomb
  val(6) = dydz
  val(7) = dxdz
  val(8) = dydz
  val(9) = rinv + dz*dz*dcomb

  return
end subroutine st3d_comb_vec


subroutine st3d_comb(srcinfo,ndt,targ,ndd,dpars,ndz, &
     zk,ndi,ipars,val)
!f2py intent(in) src,ndt,targ,ndd,dpars,ndz,zk,ndi,ipars
!f2py intent(out) val
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(12), targ(ndt),dpars(ndd)
  integer ipars(ndi)
  complex *16 :: zk
  real *8 :: val

  real *8 :: src(3), srcnorm(3), dr(3)
  !
  ! returns one entry of the Stokes double layer kernel
  !
  ! D_ij = -T_jik n_k = 3*(sum_k (targ_k-src_k) n_k)
  !               *(targ_j-src_j)(targ_i-src_i)/r^5
  !
  ! where r = |src-targ| and n is the normal vector at
  ! the source. Returns val=D_ij with i = ipars(1), j = ipars(2)

  i = ipars(1)
  j = ipars(2)

  alpha = dpars(1)
  beta = dpars(2)

  src(1)=srcinfo(1)
  src(2)=srcinfo(2)
  src(3)=srcinfo(3)
  srcnorm(1)=srcinfo(10)
  srcnorm(2)=srcinfo(11)
  srcnorm(3)=srcinfo(12)

!  call prin2('srcinfo=*',srcinfo,12)

  dx=targ(1)-src(1)
  dy=targ(2)-src(2)
  dz=targ(3)-src(3)

  dprod = dx*srcnorm(1) + dy*srcnorm(2) + dz*srcnorm(3)

  r=sqrt(dx**2+dy**2+dz**2)
  rinv = 1.0d0/r
  rinv3 = 0.5d0*rinv**3*alpha
  rinv5 = 3.0d0*dprod*rinv**5*beta
  rinv = 0.5d0*rinv*alpha

  dcomb = rinv3 + rinv5

  dr(1) = dx
  dr(2) = dy
  dr(3) = dz

  dxi = dr(i)
  dxj = dr(j)
  val = dxi*dxj*dcomb

  if (i .eq. j) val = val + rinv

  return
end subroutine st3d_comb



subroutine st3d_strac_vec(nd,srcinfo,ndt,targinfo, &
     ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(*), targinfo(12),dpars(ndd)
  integer ipars(ndi)
  real *8 :: val(9), targ(3), src(3), targnorm(3)
  complex *16 :: zk
!f2py intent(in) nd,src,ndt,targ,ndd,dpars,ndz,zk,ndi,ipars
!f2py intent(out) val

  !
  ! returns the traction of the Stokes single layer kernel
  !
  ! t(S)_ij = T_ijk n_k = -3*(sum_k (targ_k-src_k) n_k)
  !               *(targ_j-src_j)(targ_i-src_i)/r^5
  !
  ! where r = |src-targ| and n is the normal vector at the
  ! target. The output is given ordered by standard
  ! linear indexing, ij -> i+(j-1)*3

  src(1)=srcinfo(1)
  src(2)=srcinfo(2)
  src(3)=srcinfo(3)
  targ(1)=targinfo(1)
  targ(2)=targinfo(2)
  targ(3)=targinfo(3)
  targnorm(1)=targinfo(10)
  targnorm(2)=targinfo(11)
  targnorm(3)=targinfo(12)
  
!  call prin2('srcinfo=*',srcinfo,12)

  dx=targ(1)-src(1)
  dy=targ(2)-src(2)
  dz=targ(3)-src(3)

  dprod = dx*targnorm(1) + dy*targnorm(2) + dz*targnorm(3)

  r=sqrt(dx**2+dy**2+dz**2)
  rinv5 = -3.0d0*dprod*(1.0d0/r)**5

  dxdy = dx*dy*rinv5
  dxdz = dx*dz*rinv5
  dydz = dy*dz*rinv5

  val(1) = dx*dx*rinv5
  val(2) = dxdy
  val(3) = dxdz
  val(4) = dxdy
  val(5) = dy*dy*rinv5
  val(6) = dydz
  val(7) = dxdz
  val(8) = dydz
  val(9) = dz*dz*rinv5

  return
end subroutine st3d_strac_vec


subroutine st3d_strac(srcinfo,ndt,targinfo,ndd,dpars, &
     ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(*), targinfo(12),dpars(ndd)
  integer ipars(ndi)
  real *8 :: val, targ(3), src(3), targnorm(3), dr(3)
!f2py intent(in) src,ndt,targ,ndd,dpars,ndz,zk,ndi,ipars
!f2py intent(out) val

  !
  ! returns one entry of the traction of the Stokes single
  ! layer kernel
  !
  ! t(S)_ij = T_ijk n_k = -3*(sum_k (targ_k-src_k) n_k)
  !               *(targ_j-src_j)(targ_i-src_i)/r^5
  !
  ! where r = |src-targ| and n is the normal vector at the
  ! target. Returns val=t(S)_ij, with i = ipars(1), j = ipars(2)

  src(1)=srcinfo(1)
  src(2)=srcinfo(2)
  src(3)=srcinfo(3)
  targ(1)=targinfo(1)
  targ(2)=targinfo(2)
  targ(3)=targinfo(3)
  targnorm(1)=targinfo(10)
  targnorm(2)=targinfo(11)
  targnorm(3)=targinfo(12)
  
!  call prin2('srcinfo=*',srcinfo,12)

  dx=targ(1)-src(1)
  dy=targ(2)-src(2)
  dz=targ(3)-src(3)
  
  dprod = dx*targnorm(1) + dy*targnorm(2) + dz*targnorm(3)

  r=sqrt(dx**2+dy**2+dz**2)
  rinv5 = -3.0d0*dprod*(1.0d0/r)**5


  dr(1) = dx
  dr(2) = dy
  dr(3) = dz

  i = ipars(1)
  j = ipars(2)

  dxi = dr(i)
  dxj = dr(j)
  
  val = dxi*dxj*rinv5

  return
end subroutine st3d_strac


