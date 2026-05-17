!     end of debugging code
!
! vpp is a lightweight piecewise polynomial evaluator for
! vector-valued functions in one dimension. to do complex
! valued, simply double the vector length
!
! there are 3 user-callable routines
!
! vpp_build - given a function handle for a function of 1 variable
!        (possibly vector-valued) and a range [a,b], this builds
!        a set of subintervals and monomial expansions which
!        approximate the function to the given tolerance
!
! vpp_buildkern - given a subroutine handle for a function of 1 variable
!        (possibly vector-valued) with the usual fmm3dbie kernel calling
!        sequence, and a range [a,b], this builds
!        a set of subintervals and monomial expansions which
!        approximate the function to the given tolerance. The one variable
!        is passed to the subroutine as the x coordinate (first entry)
!        of srcvals
!
! vpp_eval - after the appropriate arrays are created by vpp_build
!        this routine can evaluate the expansion at any given point
!
! vpp_kern - after the appropriate arrays are created by vpp_build
!        this routine will evaluate the expansion at r where r is the
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
 

subroutine vpp_build(fun,nv,a,b,n,tol,maxsub,maxdepth,ietype, &
     ipars,ndd,ldpars,dpars,ier)

  ! ier = 1, maxsub exceeded
  ! ier = 2, maxdepth exceeded
  ! ier = 4, ldpars insufficient   

  implicit real *8 (a-h,o-z)
  real *8 :: dpars(*)
  integer :: ipars(*)

  real *8, allocatable :: subs(:,:), fs(:,:,:), fovers(:,:), &
       amat(:,:), finterp(:,:),x(:),x2(:),val(:),ak(:),as(:), &
       amono(:,:), fk(:,:,:), cfs(:,:,:),work(:)
  integer, allocatable :: ikeep(:), is(:), ik(:), ipiv(:)

  logical :: allres

  external fun

  ier = 0
  
  allocate(subs(2,2*maxsub+100),ikeep(2*maxsub+100), &
       fs(n,nv,2*maxsub+100), &
       fovers(2*n,nv),finterp(2*n,nv),x(n),x2(2*n),amat(2*n,n), &
       val(nv))

  nsubp = 0
  nsub = 1

  subs(1,1) = a
  subs(2,1) = b

  n2 = 2*n
  
  call ppchebpts(x,n)

  do j = 1,n
     x2(j) = -1 + (x(j)+1)/2d0
     x2(j+n) = (x(j)+1)/2d0
  enddo

  call ppchebinterp(n,x2,n2,amat)

  aa = subs(1,1)
  bb = subs(2,1)
  do j = 1,n
     t = aa + (x(j)+1)*(bb-aa)/2d0
     call fun(t,val)
     do k = 1,nv
        fs(j,k,1) = val(k)
     enddo
  enddo
  
  do ii = 1,maxdepth

     istart = nsubp+1
     iend = nsub

     nsubp = nsub

     allres = .true.
     
     do jj = istart,iend
        ! check resolution, divide if necessary
        aa = subs(1,jj)
        bb = subs(2,jj)
        do j = 1,n2
           t = aa + (x2(j)+1)*(bb-aa)/2d0
           call fun(t,val)
           do k = 1,nv
              fovers(j,k) = val(k)
           enddo
        enddo

        call dmatmat(n2,n,amat,nv,fs(1,1,jj),finterp)
        call ppisres(ietype,fovers,finterp,n2*nv,tol,ikeep(jj))

        if (ikeep(jj) .eq. 0) then
           if (nsub + 2 .gt. maxsub) then
              ier = 1
              return
           endif
           amid = (aa+bb)/2
           subs(1,nsub+1) = aa
           subs(2,nsub+1) = amid           
           subs(1,nsub+2) = amid
           subs(2,nsub+2) = bb

           ! save fun evals
           do j = 1,n
              do k = 1,nv
                 fs(j,k,nsub+1) = fovers(j,k)
                 fs(j,k,nsub+2) = fovers(j+n,k)
              enddo
           enddo

           nsub = nsub+2
           allres = .false.

        endif
        
     enddo

     if (allres) exit

  enddo
  
  if (.not. allres) then
     ier = 2
     return
  end if

  nkeep = 0
  do j = 1,nsub
     nkeep = nkeep + ikeep(j)
  enddo

  ndd = nkeep*nv*n + nkeep + 1
  if (ldpars .lt. 0) then
     ! query
     return
  endif

  if (ldpars .lt. ndd) then
     ier = 4
     return
  endif

  allocate(ak(nkeep),ik(nkeep),is(nkeep),as(nkeep),fk(n,nv,nkeep), &
       cfs(n,nv,nkeep),amono(n,n),ipiv(n),work(n))

  ii = 0
  do j = 1,nsub
     if (ikeep(j) .eq. 1) then
        ii = ii + 1
        ak(ii) = subs(1,j)
        ik(ii) = j
     endif
  enddo

  call ppsort(ak,nkeep,as,is)

  ibins = 1
  do j = 1,nkeep
     do k = 1,nv
        do i = 1,n
           fk(i,k,j) = fs(i,k,ik(is(j)))
        enddo
     enddo

     dpars(ibins+j-1) = as(j)
  enddo

  dpars(ibins+nkeep) = subs(2,ik(is(nkeep)))
  
  icfs = nkeep+1+ibins
  ipars(1) = nkeep
  ipars(2) = n
  ipars(3) = ibins
  ipars(4) = icfs
  ipars(5) = nv

  do k = 1,nkeep
     do j = 1,nv
        do i = 1,n
           cfs(i,j,k) = fk(i,j,k)
        enddo
     enddo
  enddo

  call ppmonomat(x,n,amono)
  k = nv*nkeep
  call dgeco_sa(amono,n,n,ipiv,dcond,work)
  job = 0
  do k = 1,nkeep
     do j = 1,nv
        call dgesl_sa(amono,n,n,ipiv,cfs(1,j,k),job)
     enddo
  enddo

  call pp_dcopy_scal(cfs,dpars(icfs),n,nv,nkeep,dpars(ibins))

  return
end subroutine vpp_build

subroutine vpp_buildkern(kern,ndd0,dpars0,ndz0,zpars0,ndi0,ipars0, &
     nv,a,b,n,tol,maxsub,maxdepth,ietype, &
     ipars,ndd,ldpars,dpars,ier)

  ! ier = 1, maxsub exceeded
  ! ier = 2, maxdepth exceeded
  ! ier = 4, ldpars insufficient   

  implicit real *8 (a-h,o-z)
  real *8 :: dpars(*), dpars0(ndd0)
  integer :: ipars(*), ipars0(ndi0)
  complex *16 :: zpars(ndz0)

  real *8 :: src(3), targ(3)
  real *8, allocatable :: subs(:,:), fs(:,:,:), fovers(:,:), &
       amat(:,:), finterp(:,:),x(:),x2(:),val(:),ak(:),as(:), &
       amono(:,:), fk(:,:,:), cfs(:,:,:), work(:)
  integer, allocatable :: ikeep(:), is(:), ik(:), ipiv(:)

  logical :: allres

  external kern

  src(:) = 0
  targ(:) = 0
  
  ndt = 3
  
  ier = 0
  
  allocate(subs(2,2*maxsub+100),ikeep(2*maxsub+100), &
       fs(n,nv,2*maxsub+100), &
       fovers(2*n,nv),finterp(2*n,nv),x(n),x2(2*n),amat(2*n,n), &
       val(nv))

  nsubp = 0
  nsub = 1

  subs(1,1) = a
  subs(2,1) = b

  n2 = 2*n
  
  call ppchebpts(x,n)

  do j = 1,n
     x2(j) = -1 + (x(j)+1)/2d0
     x2(j+n) = (x(j)+1)/2d0
  enddo

  call ppchebinterp(n,x2,n2,amat)

  aa = subs(1,1)
  bb = subs(2,1)
  do j = 1,n
     t = aa + (x(j)+1)*(bb-aa)/2d0
     targ(1) = t
     call kern(src,ndt,targ,ndd0,dpars0,ndz0,zpars0,ndi0,ipars0,val)
     do k = 1,nv
        fs(j,k,1) = val(k)
     enddo
  enddo
  
  do ii = 1,maxdepth

     istart = nsubp+1
     iend = nsub

     nsubp = nsub

     allres = .true.
     
     do jj = istart,iend
        ! check resolution, divide if necessary
        aa = subs(1,jj)
        bb = subs(2,jj)
        do j = 1,n2
           t = aa + (x2(j)+1)*(bb-aa)/2d0
           targ(1) = t
           call kern(src,ndt,targ,ndd0,dpars0,ndz0,zpars0, &
                ndi0,ipars0,val)           
           do k = 1,nv
              fovers(j,k) = val(k)
           enddo
        enddo

        call dmatmat(n2,n,amat,nv,fs(1,1,jj),finterp)
        call ppisres(ietype,fovers,finterp,n2*nv,tol,ikeep(jj))

        if (ikeep(jj) .eq. 0) then
           if (nsub + 2 .gt. maxsub) then
              ier = 1
              return
           endif
           amid = (aa+bb)/2
           subs(1,nsub+1) = aa
           subs(2,nsub+1) = amid           
           subs(1,nsub+2) = amid
           subs(2,nsub+2) = bb

           ! save fun evals
           do j = 1,n
              do k = 1,nv
                 fs(j,k,nsub+1) = fovers(j,k)
                 fs(j,k,nsub+2) = fovers(j+n,k)
              enddo
           enddo

           nsub = nsub+2
           allres = .false.

        endif
        
     enddo

     if (allres) exit

  enddo
  
  if (.not. allres) then
     ier = 2
     return
  end if

  nkeep = 0
  do j = 1,nsub
     nkeep = nkeep + ikeep(j)
  enddo

  ndd = nkeep*nv*n + nkeep + 1
  if (ldpars .lt. 0) then
     ! query
     return
  endif

  if (ldpars .lt. ndd) then
     ier = 4
     return
  endif

  allocate(ak(nkeep),ik(nkeep),is(nkeep),as(nkeep),fk(n,nv,nkeep), &
       cfs(n,nv,nkeep),amono(n,n),ipiv(n),work(n))

  ii = 0
  do j = 1,nsub
     if (ikeep(j) .eq. 1) then
        ii = ii + 1
        ak(ii) = subs(1,j)
        ik(ii) = j
     endif
  enddo

  call ppsort(ak,nkeep,as,is)

  ibins = 1
  do j = 1,nkeep
     do k = 1,nv
        do i = 1,n
           fk(i,k,j) = fs(i,k,ik(is(j)))
        enddo
     enddo

     dpars(ibins+j-1) = as(j)
  enddo

  dpars(ibins+nkeep) = subs(2,ik(is(nkeep)))
  
  icfs = nkeep+1+ibins
  ipars(1) = nkeep
  ipars(2) = n
  ipars(3) = ibins
  ipars(4) = icfs
  ipars(5) = nv

  do k = 1,nkeep
     do j = 1,nv
        do i = 1,n
           cfs(i,j,k) = fk(i,j,k)
        enddo
     enddo
  enddo

  call ppmonomat(x,n,amono)
  call dgeco_sa(amono,n,n,ipiv,dcond,work)
  job = 0
  do k = 1,nkeep
     do j = 1,nv
        call dgesl_sa(amono,n,n,ipiv,cfs(1,j,k),job)
     enddo
  enddo
  
  call pp_dcopy_scal(cfs,dpars(icfs),n,nv,nkeep,dpars(ibins))

  return
end subroutine vpp_buildkern

subroutine pp_dcopy_scal(x,y,n,nv,nkeep,bins)
  implicit real *8 (a-h,o-z)
  real *8 :: x(n,nv,nkeep), y(nv,n,nkeep), bins(nkeep+1)

  do k = 1,nkeep
     h = 1d0/(bins(k+1)-bins(k))
     do j = 1,nv
        hp = 1d0
        do i = 1,n
           y(j,n-i+1,k) = x(n-i+1,j,k)*hp
           hp = hp*h
        enddo
     enddo
  enddo

  return
end subroutine pp_dcopy_scal

subroutine ppisres(ietype,f1,f2,n,tol,ikeep)
  ! in some cases, important that f1 considered "true" value
  implicit real *8 (a-h,o-z)
  real *8 :: f1(*), f2(*)

  ikeep = 1
  
  derr = 0
  drel = 0
  do j = 1,n
     derr = max(derr,abs(f1(j)-f2(j)))
     drel = max(drel,abs(f1(j)))
  enddo

  derr = derr/(max(drel,1d0))
  if (derr .gt. tol) ikeep = 0
  
  return
end subroutine ppisres

subroutine vpp_kern(src,ndt,targ,ndd,dpars,ndz,zpars,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: src(*), targ(ndt), dpars(ndd)
  integer ipars(ndi)
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
  real *8 :: dpars(ndd)
  integer ipars(ndi)
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
  real *8 :: r, as(nbin+1), a
  integer :: nbin, ibin

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


subroutine ppchebpts(x,n)
  implicit real *8 (a-h,o-z)
  real *8 :: x(n)

  pi = 4*atan(1d0)

  do j = 1,n
     x(j) = cos( (2*(j-1) + 1)*pi/(2d0*n))
  enddo

  return
end subroutine ppchebpts

subroutine ppchebbary(x,n)
  implicit real *8 (a-h,o-z)
  real *8 :: x(n)

  pi = 4*atan(1d0)

  isgn = 1
  do j = 1,n
     x(j) = isgn*sin( (2*(j-1) + 1)*pi/(2d0*n))
     isgn = -isgn
  enddo

  return
end subroutine ppchebbary

subroutine ppmonomat(x,n,amat)
  implicit real *8 (a-h,o-z)
  real *8 :: amat(n,n), x(n)

  do j = 1,n
     t = (x(j)+1)/2d0
     tp = 1d0
     do i = 1,n
        amat(j,n-i+1) = tp
        tp = tp*t
     enddo
  enddo
  
  return
end subroutine ppmonomat

subroutine ppsort(a,n,as,is)
  implicit real *8 (a-h,o-z)
  real *8 :: a(n), as(n)
  integer :: is(n)

  real *8, allocatable :: at(:)
  integer, allocatable :: it(:)

  allocate(at(n),it(n))
  
  nb = 1

  do i = 1,n
     as(i) = a(i)
     at(i) = a(i)
     is(i) = i
     it(i) = i
  enddo

  ! lazy merge sort 
  do jjj = 1,1000
     nb2 = 2*nb
     nb2s = (n-1)/nb2 + 1

     do j = 1,nb2s
        i1 = 1
        i2 = 1

        do k = 1,nb2
           j1 = i1 + (j-1)*nb2
           j2 = i2 + (j-1)*nb2 + nb
           l = (j-1)*nb2 + k
           if (l .gt. n) exit
           if (i1 .gt. nb .and. i2 .gt. nb) exit
           if (i1 .gt. nb .or. j1 .gt. n) then
              as(l) = at(j2)
              is(l) = it(j2)
              i2 = i2 + 1
              cycle
           endif
           if (i2 .gt. nb .or. j2 .gt. n) then
              as(l) = at(j1)
              is(l) = it(j1)
              i1 = i1 + 1
              cycle
           endif
           if (at(j2) .lt. at(j1)) then
              as(l) = at(j2)
              is(l) = it(j2)
              i2 = i2 + 1
           else
              as(l) = at(j1)
              is(l) = it(j1)
              i1 = i1 + 1
           endif
        enddo
     enddo

     do j = 1,n
        at(j) = as(j)
        it(j) = is(j)
     enddo
     
     nb = nb2
     if (nb .gt. 2*n) exit

  enddo
  return
end subroutine ppsort

subroutine ppchebinterp(n,x2,n2,amat)
  implicit real *8 (a-h,o-z)
  real *8 :: x2(n2), amat(n2,n)

  real *8, allocatable :: x(:), bw(:)

  allocate(x(n),bw(n))

  call ppchebpts(x,n)
  call ppchebbary(bw,n)

  do i = 1,n2
     xx = x2(i)
     dtot = 0

     ifix = 0
     jfix = 0
     
     do j = 1,n
        amat(i,j) = bw(j)/(xx-x(j))
        dtot = dtot + amat(i,j)
        if (xx .eq. x(j)) then
           ifix = 1
           jfix = j
        endif
     enddo

     do j = 1,n
        amat(i,j) = amat(i,j)/dtot
     enddo

     if (ifix .ne. 0) then
        do j = 1,n
           amat(i,j) = 0
           if (j .eq. jfix) amat(i,j) = 1
        enddo
     endif
  enddo

  return
end subroutine ppchebinterp
      
