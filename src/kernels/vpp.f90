!
! vpp is a lightweight piecewise polynomial evaluator for
! vector-valued functions in one dimension. to do complex
! valued, simply double the vector length
!
! The expansions are built in matlab by radcheb_fit.m (see
! @kernel3d/radcheb.m); only the evaluators live here.
!
! vpp_eval - evaluates the expansion at a given point r
!
! vpp_kern - evaluates the expansion at r where r is the
!        distance between given source and target points. it has the
!        correct interface for an fmm3dbie kernel
!
! these do no bounds-checking etc.


!
! developer notes:
!
! - [a,b] is subdivided into sub-intervals.
! - function represented by shifted monomial series on each
!  subinterval, with coefficients ordered highest to lowest degree.
! - on [c,d] with coefficients cf, the polynomial representation is
!     cf(1)*(x-c)^(n-1) + cf(2)*(x-c)^(n-2) + ... cf(n)
! - vector-valued function coefs are interlaced, i.e. shape (nv,n)
!  where nv is dimension of vector and n is number of terms in expansion
! - ipars array stores pointers and basic info
! --   ipars(1) = number of subintervals
! --   ipars(2) = number of terms in poly expansions (degree + 1)
! --   ipars(3) = start of subinterval endpoints in dpars 
! --   ipars(4) = start of coefficient storage in dpars
! --   ipars(5) = dimension of vector valued data (1 for scalar)
! - dpars stores subinterval endpoints and coefficients 
 



subroutine vpp_kern(src,ndt,targ,ndd,dpars,ndz,zpars,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: src(*), targ(ndt), dpars(ndd)
  integer *8 ipars(ndi)
  real *8 :: val(*)
  complex *16 :: zpars(ndz)

  dx=targ(1)-src(1)
  dy=targ(2)-src(2)
  dz=targ(3)-src(3)

  r=sqrt(dx**2+dy**2+dz**2)

  call vpp_eval(r,ndd,dpars,ndi,ipars,val)
  
  return
end subroutine vpp_kern

subroutine vpp_eval(r,ndd,dpars,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: dpars(ndd)
  integer *8 ipars(ndi)
  real *8 :: val(*)

  nbin = ipars(1)
  n = ipars(2)
  ibins = ipars(3)
  icfs = ipars(4)
  nv = ipars(5)
  
  call ipp_getbin(r,dpars(ibins),nbin,ibin,a)
  
  call vpp_eval0(r,a,dpars(icfs+n*nv*(ibin-1)),n,nv,val)
  
  return
end subroutine vpp_eval

subroutine vpp_eval0(r,a,dcfs,n,nv,val)
  ! coefficients are ordered highest degree to lowest.
  ! scaled so that polynomial is
  !     cf(1)*(x-a)^(n-1) + cf(2)*(x-a)^(n-2) + ... cf(n)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: dcfs(nv,n), val(*)

  x = r-a
  do j = 1,nv
     val(j) = dcfs(j,1)
  enddo
  do i = 2,n
     do j = 1,nv
        val(j) = dcfs(j,i) + x*val(j)
     enddo
  enddo
  
  return
end subroutine vpp_eval0

subroutine ipp_getbin(r,as,nbin,ibin,a)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: r, as(nbin+1), a
  integer *8 :: nbin, ibin

  ibin = 1
  jbin = nbin

  a = as(ibin)
  a2 = as(jbin)
  if (r .lt. a) then
     ibin = 1
     return
  endif
  if (r .ge. a2) then
     ibin = nbin
     a = a2
     return
  endif

  do i = 1,nbin
     if (jbin .le. ibin + 1) exit
     midbin = (ibin + jbin)/2
     am = as(midbin)
     
     if (r .ge. am) then
        a = am
        ibin = midbin
     else
        a2 = am
        jbin = midbin
     end if
  end do
  
  return
end subroutine ipp_getbin
