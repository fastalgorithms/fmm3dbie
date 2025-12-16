! (c) Mike O'Neil
!
! this file contains basic common-sense wrappers for lapack routines,
! that way you never need to call the complicated routine with
! several parameters...
!
! gaussian elimination
!
!     dgausselim - real matrix, real rhs, real solution
!     dgausselim_zrhs - real matrix, complex rhs, compex solution
!
subroutine dcopy_guru(n,a,incx,b,incy)
  implicit real *8 (a-h,o-z)

  integer *8 n,incx,incy
  integer n1,incx1,incy1

  real *8 a(*),b(*)

  n1 = n
  incx1 = incx
  incy1 = incy
  
  call dcopy(n1, a, incx1, b, incy1)

  return
end subroutine dcopy_guru
!
!
!
!

subroutine zvecvec(n, x, y, cd)
!
!-----------------
!  This subroutine computes the dot product of two complex vectors
!
!  Input arguments:
!  
!    - n: integer
!        length of the complex vectors
!    - x: double complex(n)
!        input vector 1
!    - y: double complex(n)
!        input vector 2
!
!  Output arguments:
!
!    - cd: double complex
!        inner product \sum_{i} x(i) y(i)
!---------------------------        
!
!

  implicit double precision (a-h,o-z)
  integer *8 n
  complex *16, intent(in) :: x(n), y(n)
  complex *16, intent(out) :: cd
  complex *16 zdotu
  integer n1

  n1 = n
  cd = 0
  cd = zdotu(n1,x,1,y,1)
  return
end subroutine zvecvec





subroutine zzero(n, z)
!
!--------------------
!  This subroutine zeros out a complex vector
!
!  Input arguments:
!
!    - n: integer
!        length of the vector
!
!  Output arguments:
!
!    - z: double complex(n)
!        zeroed out vector
!---------------------
!
  implicit double precision (a-h,o-z)
  integer *8 :: n
  complex *16 :: z(n)


  do i = 1,n
    z(i) = 0
  end do
  return
end subroutine zzero





subroutine ztranspose(m, n, a)
!
!---------------------
!  This subroutine returns an inplace transpose
!  of a complex matrix
!
!  Input arguments:
!    - m: integer
!        number of rows
!    - n: integer
!        number of columns
!
!  Inout arguments:
!    - a: double complex(m,n)
!        on input, the matrix a, on output
!        the transpose of matrix a
!-------------------------
!
  implicit double precision (a-h,o-z)
  integer *8, intent(in) :: m,n
  double complex, intent(inout) :: a(*)
  double complex, allocatable :: at(:)
  integer *8 i,j

  integer mn, incx, incy

  mn = m*n

  allocate(at(mn))

  do j = 1,n
    do i = 1,m
      at(j+n*(i-1)) = a(i+m*(j-1))
    end do
  end do

  incx = 1
  incy = 1
  call zcopy(mn, at, incx, a, incy)

  deallocate(at)
  
  return
end subroutine ztranspose



subroutine zdiff(n, x, y, z)
  implicit double precision (a-h,o-z)
  integer *8 n,i
  complex *16 :: x(n), y(n), z(n)

  do i = 1,n
    z(i) = x(i)-y(i)
  end do
  return
end subroutine zdiff





subroutine zratios(n, x, y, z)
  implicit double precision (a-h,o-z)
  integer *8 n,i
  complex *16 :: x(n), y(n), z(n)

  do i = 1,n
    z(i) = x(i)/y(i)
  end do
  return
end subroutine zratios






subroutine zeigs(n, a, info, vleft, zlams, vright)
  implicit double precision (a-h,o-z)
  integer *8 n,i,j
  complex *16 :: a(n,n), vleft(n,n), zlams(n), vright(n,n)

  character :: jobvl, jobvr
  double precision, allocatable :: dwork(:)
  complex *16, allocatable :: work(:), atemp(:,:)

  integer nuse

  jobvl = 'V'
  jobvr = 'V'
  ldvl = n
  ldvr = n
  lwork = 1000000

  nuse = n

  lwork = 10*n
  allocate(work(lwork))
  allocate(dwork(lwork))
  allocate(atemp(n,n))

  do i = 1,n
    do j = 1,n
      atemp(i,j) = a(i,j)
    end do
  end do

  call zgeev(jobvl, jobvr, nuse, atemp, nuse, zlams, vleft, &
      ldvl, vright, ldvr, work, lwork, dwork, info)  

  if (info .ne. 0) then
    !call prinf('in zeigs, info = *', info, 1)
    print *, "in zeigs, info = ", info
    stop
  endif
  
  return
end subroutine zeigs





subroutine zsvd(m, n, a, u, s, vt)
  implicit double precision (a-h,o-z)
  integer *8 m,n,i,j
  double precision :: s(*)
  complex *16 :: a(m,n), u(*), vt(*)

  character :: jobu, jobvt
  double precision, allocatable :: dwork(:)
  complex *16, allocatable :: work(:), atemp(:,:)

  integer muse,nuse
  
  !
  ! note: V^* is returned, not V
  !
  jobu = 'S'
  jobvt = 'S'
  lda = m
  ldu = m
  ldvt = min(m,n)

  lwork = 2*min(m,n) + max(m,n) + 1000
  allocate(work(lwork))
  allocate(dwork(5*min(m,n)+1000))
  allocate(atemp(m,n))

  do i = 1,m
    do j = 1,n
      atemp(i,j) = a(i,j)
    end do
  end do

  muse = m
  nuse = n
  
  call zgesvd(jobu, jobvt, muse, nuse, atemp, lda, s, u, ldu, &
      vt, ldvt, work, lwork, dwork, info)

  if (info .ne. 0) then
    print *, "in zsvd, info = ", info
    !!call prinf('in zsvd, info = *', info, 1)
    stop
  endif
  
  return
end subroutine zsvd





subroutine dsvd(m, n, a, u, s, vt)
  implicit double precision (a-h,o-z)
  integer *8 m,n,i,j
  double precision :: s(*)
  double precision :: a(m,n), u(*), vt(*)

  character :: jobu, jobvt
  double precision, allocatable :: work(:), atemp(:,:)

  integer muse,nuse
  
  !
  ! note: V^* is returned, not V
  !
  jobu = 'S'
  jobvt = 'S'
  lda = m
  ldu = m
  ldvt = min(m,n)

  lwork = 5*max(m,n) + 1000
  allocate(work(lwork))
  allocate(atemp(m,n))

  do i = 1,m
    do j = 1,n
      atemp(i,j) = a(i,j)
    end do
  end do

  muse = m
  nuse = n
  
  call dgesvd(jobu, jobvt, muse, nuse, atemp, lda, s, u, ldu, &
      vt, ldvt, work, lwork, info)

  if (info .ne. 0) then
    print *, "in dsvd, info = ", info
    !!!!call prinf('in dsvd, info = *', info, 1)
    stop
  endif
  
  return
end subroutine dsvd





subroutine dinverse(n, a, info, ainv)
  implicit double precision (a-h,o-z)
  integer *8 n, info
  double precision :: a(n,n), ainv(n,n)

  integer nuse, incx, incy, infouse  
  double precision, allocatable :: work(:)
  integer, allocatable :: ipiv(:)

  !
  ! call the LAPACK routine to compute the inverse of the matrix a
  !
  ! Input:
  !   n - dimension of a
  !   a - matrix to invert, it is not harmed
  !
  ! Output:
  !   info - lapack info
  !   ainv - the inverse of a, computed using Gaussian elimination
  !
  !

  nuse = n
  incx = 1
  incy = 1
  call dcopy(nuse*nuse, a, incx, ainv, incy)
  allocate(ipiv(nuse))
  call dgetrf(nuse, nuse, ainv, nuse, ipiv, infouse)

  lwork = 10*nuse
  allocate(work(lwork))
  call dgetri(nuse, ainv, nuse, ipiv, work, lwork, infouse)

  info = infouse

  return
end subroutine dinverse





subroutine zinverse(n, a, info, ainv)
  implicit double precision (a-h,o-z)
  integer *8 n,info
  complex *16 :: a(n,n), ainv(n,n)

  integer nuse,infouse
  integer, allocatable :: ipiv(:)
  complex *16, allocatable :: work(:)

  !
  ! call the LAPACK routine to compute the inverse of the matrix a
  !
  ! Input:
  !   n - dimension of a
  !   a - matrix to invert, it is not harmed
  !
  ! Output:
  !   info - lapack info
  !   ainv - the inverse of a, computed using Gaussian elimination
  !
  !
  nuse = n
  incx = 1
  incy = 1
  call zcopy(nuse*nuse, a, incx, ainv, incy)

  allocate(ipiv(nuse))
  call zgetrf(nuse, nuse, ainv, nuse, ipiv, infouse)

  lwork = 20*nuse
  allocate(work(lwork))
  call zgetri(nuse, ainv, nuse, ipiv, work, lwork, infouse)

  info = infouse

  return
end subroutine zinverse





subroutine dgausselim(n, a, rhs, info, sol, dcond)
  implicit double precision (a-h,o-z)
  integer *8 n, info
  double precision :: a(n,n), rhs(n), sol(n)
  
  integer nuse, infouse
  integer, allocatable :: ipiv(:), iwork(:)
  double precision, allocatable :: af(:,:), rscales(:), cscales(:), work(:)
  character *1 :: fact, trans, equed
  !
  ! just a wrapper for the lapack routine...
  ! a is untouched
  !
  nuse = n
  allocate(ipiv(nuse))
  allocate(af(nuse,nuse))
  allocate(rscales(nuse))
  allocate(cscales(nuse))
  allocate(work(5*nuse))
  allocate(iwork(nuse))
  
  fact = 'N'
  trans = 'N'
  nrhs = 1
  lda = n
  ldaf = n
  equed = 'N'
  ldb = n
  ldx = n
  
  call dgesvx( fact, trans, nuse, nrhs, a, lda, af, ldaf, ipiv, &
      equed, rscales, cscales, rhs, ldb, sol, ldx, dcond, ferr, berr, &
      work, iwork, infouse )

  dcond = 1/dcond
  info = infouse
  
  return
end subroutine dgausselim


subroutine dgausselim_vec(n, a, k, rhs, info, sol, dcond)
  implicit double precision (a-h,o-z)
  integer *8 n,k, info
  double precision :: a(n,n), rhs(n,k), sol(n,k)
 
  integer nuse,kuse,infouse
  integer, allocatable :: ipiv(:), iwork(:)
  double precision, allocatable :: af(:,:), rscales(:), cscales(:), work(:)
  character *1 :: fact, trans, equed
  !
  ! This a wrapper for double precision Gaussian elimination in
  ! LAPACK for multipole right hand sides
  ! 
  !
  
  nuse = n
  kuse = k
  allocate(ipiv(n))
  allocate(af(n,n))
  allocate(rscales(n))
  allocate(cscales(n))
  allocate(work(5*n))
  allocate(iwork(n))
  
  fact = 'N'
  trans = 'N'
  nrhs = k
  lda = n
  ldaf = n
  equed = 'N'
  ldb = n
  ldx = n
  
  call dgesvx( fact, trans, nuse, nrhs, a, lda, af, ldaf, ipiv, &
      equed, rscales, cscales, rhs, ldb, sol, ldx, dcond, ferr, berr, &
      work, iwork, infouse )

  dcond = 1/dcond
  info = infouse
  
  return
end subroutine dgausselim_vec


subroutine dleastsq(m,n,a,nrhs,rhs,eps,info,sol,irank)
  !
  ! a wrapper for dgelsy. a and rhs are not destroyed
  !
  ! eps is the rank cut-off for solving the system
  ! in general, eps can be set near machine precision
  ! recommended: 1d-15 .lt. eps .lt. 1d-12
  !
  implicit none
  integer *8 m,n,info,nrhs,irank
  real *8 a(m,n), rhs(m,nrhs), sol(n,nrhs), eps
  ! local  
  real *8, allocatable :: work(:), atemp(:,:), rhstemp(:,:)
  integer, allocatable :: ipiv(:)
  real *8 dcond

  integer muse,nuse,infouse,nrhsuse,irankuse
  integer lwork, lda, ldb, i, j, mn, mn2, nb, nb1, nb2, nb3, nb4
  integer lwkmin, lwkopt

  integer ilaenv


  muse = m
  nuse = n
  nrhsuse = nrhs

  dcond = 1d0/eps
  dcond = eps

  dcond = 1d-8
  
  mn = min(m,n)
  mn2 = max(m,n)
  
  if( mn.eq.0 .or. nrhs.eq.0 ) then
     lwkmin = 1
     lwkopt = 1
  else
     nb1 = ilaenv( 1, 'DGEQRF', ' ', m, n, -1, -1 )
     nb2 = ilaenv( 1, 'DGERQF', ' ', m, n, -1, -1 )
     nb3 = ilaenv( 1, 'DORMQR', ' ', m, n, nrhs, -1 )
     nb4 = ilaenv( 1, 'DORMRQ', ' ', m, n, nrhs, -1 )
     nb = max( nb1, nb2, nb3, nb4 )
     lwkmin = mn + max( 2*mn, n + 1, mn + nrhs )
     lwkopt = max( lwkmin, mn + 2*n + nb*( n + 1 ), 2*mn + nb*nrhs )
  endif
  
  lwork = lwkopt
  
  allocate(work(lwork),ipiv(n),atemp(m,n),rhstemp(mn2,nrhs))

  do i = 1,n
     ipiv(i) = 0
     do j = 1,m
        atemp(j,i) = a(j,i)
     enddo
  enddo

  do i = 1,nrhs
     do j = 1,m
        rhstemp(j,i) = rhs(j,i)
     enddo
  enddo

  lda = m
  ldb = mn2

  call dgelsy(muse,nuse,nrhsuse,atemp,lda,rhstemp,ldb,ipiv,dcond, &
       irankuse,work,lwork,infouse)

  do i = 1,nrhs
     do j = 1,n
        sol(j,i) = rhstemp(j,i)
     enddo
  enddo

  irank = irankuse
  info = infouse

  return
end subroutine dleastsq
  


subroutine zgausselim(n, a, rhs, info, sol, dcond)
  implicit double precision (a-h,o-z)
  integer *8 n,info
  complex *16 :: a(n,n), rhs(n), sol(n)
  
  integer nuse,infouse
  integer, allocatable :: ipiv(:)
  double precision, allocatable :: rscales(:), cscales(:), rwork(:)
  complex *16, allocatable :: af(:,:), work(:)
  character *1 :: fact, trans, equed
  !
  ! just a wrapper for the lapack routine...
  ! a is untouched
  !
  !call prinf('n = *', n, 1)
  !call prin2('amat = *', a, n*n*2)
  !call prin2('rhs = *', rhs, 2*n)
  !stop
  
  nuse = n
  allocate(ipiv(n))
  allocate(af(n,n))
  allocate(rscales(n))
  allocate(cscales(n))
  allocate(work(2*n))
  allocate(rwork(2*n))

  
  fact = 'N'
  trans = 'N'
  nrhs = 1
  lda = n
  ldaf = n
  equed = 'N'
  ldb = n
  ldx = n
  
  call zgesvx(fact, trans, nuse, nrhs, a, lda, af, ldaf, ipiv, &
      equed, rscales, cscales, rhs, ldb, sol, ldx, dcond, ferr, berr, &
      work, rwork, infouse)
  dcond = 1/dcond

  info = infouse

  deallocate(ipiv, af, rscales, cscales, work, rwork)
  
  return
end subroutine zgausselim






subroutine dgausselim_zrhs(n, a, rhs, info, sol, dcond)
  implicit double precision (a-h,o-z)
  integer *8 n, info
  double precision :: a(n,n), rhs(2,n), sol(2,n)
  
  integer nuse, infouse
  integer, allocatable :: ipiv(:), iwork(:)
  double precision, allocatable :: af(:,:), rscales(:), cscales(:), work(:)
  double precision, allocatable :: x(:,:), b(:,:)
  character *1 :: fact, trans, equed
  !
  ! just a wrapper for the lapack routine...
  ! a is untouched
  !
  
  nuse = n
  allocate(ipiv(n))
  allocate(af(n,n))
  allocate(rscales(n))
  allocate(cscales(n))
  allocate(b(n,2), x(n,2))
  allocate(work(5*n))
  allocate(iwork(n))
  
  fact = 'N'
  trans = 'N'
  nrhs = 2
  lda = n
  ldaf = n
  equed = 'N'
  ldb = n
  ldx = n

  do i = 1,n
    b(i,1) = rhs(1,i)
    b(i,2) = rhs(2,i)
  end do
  
  call dgesvx( fact, trans, nuse, nrhs, a, lda, af, ldaf, ipiv, &
      equed, rscales, cscales, b, ldb, x, ldx, dcond, ferr, berr, &
      work, iwork, infouse )

  do i = 1,n
    sol(1,i) = x(i,1)
    sol(2,i) = x(i,2)
  end do

  info = infouse
  
  dcond = 1/dcond

  !call prin2('x = *', x, 2*n)
  !stop
  
  return
end subroutine dgausselim_zrhs





subroutine zmatvec(m, n, a, x, y)
  implicit double precision (a-h,o-z)
  integer *8 :: m,n
  complex *16 :: a(m,n), x(n), y(m)
  character *1 :: trans
  complex *16 :: alpha, beta
  integer :: muse,nuse

  muse = m
  nuse = n

  trans = 'N'
  alpha = 1
  incx = 1
  beta = 0
  incy = 1

  if (m.eq.0.or.n.eq.0) return
  call zgemv(trans, muse, nuse, alpha, a, muse, x, incx, beta, y, incy)

  return
end subroutine zmatvec





subroutine dmatvec(m, n, a, x, y)
  implicit double precision (a-h,o-z)
  integer *8 :: m, n
  double precision :: a(m,n), x(n), y(m)
  character *1 :: trans
  double precision :: alpha, beta
  integer :: muse, nuse

  muse = m
  nuse = n
  trans = 'N'
  alpha = 1
  incx = 1
  beta = 0
  incy = 1

  if (m.eq.0.or.n.eq.0) return
  call dgemv(trans, muse, nuse, alpha, a, muse, x, incx, beta, y, incy)

  return
end subroutine dmatvec




subroutine dmatzvec(m, n, a, x, y)
  implicit double precision (a-h,o-z)
  integer *8 :: m, n
  double precision :: a(m,n)
  complex *16 :: x(n), y(m)

  integer :: muse, nuse
  complex *16 :: cd

  muse = m
  nuse = n

  do i = 1,muse
    cd = 0
    do j = 1,nuse
      cd = cd + a(i,j)*x(j)
    end do
    y(i) = cd
  end do

  return
end subroutine dmatzvec





subroutine dmatzvec_saved(m, n, a, x, y)
  implicit double precision (a-h,o-z)
  integer *8 m, n
  double precision :: a(m,n)
  complex *16 :: x(n), y(m)
  character *1 :: trans
  double precision :: alpha, beta

  integer muse,nuse

  double precision, allocatable :: xr(:), xi(:), yr(:), yi(:)
  complex *16 :: ima

  muse = m
  nuse = n

  allocate(xr(n), xi(n))
  allocate(yr(m), yi(m))
  ima = (0,1)
  
  trans = 'N'
  alpha = 1
  incx = 1
  beta = 0
  incy = 1

  do i = 1,n
    xr(i) = dble(x(i))
    xi(i) = imag(x(i))
  end do

  if(m.eq.0.or.n.eq.0) return
  call dgemv(trans, muse, nuse, alpha, a, muse, xr, incx, beta, yr, &
      incy)
  call dgemv(trans, muse, nuse, alpha, a, muse, xi, incx, beta, yi, &
      incy)

  do i = 1,m
    y(i) = yr(i) + ima*yi(i)
  end do

  deallocate(xr, xi, yr, yi)
  
  return
end subroutine dmatzvec_saved





subroutine zmatmat(m, n, a, k, b, c)
  implicit double precision (a-h,o-z)
  integer *8 :: m, n, k
  complex *16 :: a(m,n), b(n,k), c(m,k)
  character *1 :: transa, transb
  complex *16 :: alpha, beta
  integer :: muse, nuse, kuse

  !
  ! note different dimensions than usual...
  !

  muse = m
  nuse = n
  kuse = k
  transa = 'N'
  transb = 'N'
  alpha = 1
  beta = 0
  lda = m
  ldb = n
  ldc = m

  if(m.eq.0.or.n.eq.0.or.k.eq.0) return
  call zgemm(transa, transb, muse, kuse, nuse, alpha, a, lda, b, ldb, &
      beta, c, ldc)

  return
end subroutine zmatmat





subroutine dmatmat(m, n, a, k, b, c)
  implicit double precision (a-h,o-z)
  integer *8 :: m, n, k
  double precision :: a(m,n), b(n,k), c(m,k)
  character *1 :: transa, transb
  double precision :: alpha, beta
  integer :: muse, nuse, kuse

  !
  ! note different dimensions than usual...
  !
  muse = m
  nuse = n
  kuse = k

  transa = 'N'
  transb = 'N'
  alpha = 1
  beta = 0
  lda = m
  ldb = n
  ldc = m

  if(m.eq.0.or.n.eq.0.or.k.eq.0) return
  call dgemm(transa, transb, muse, kuse, nuse, alpha, a, lda, b, ldb, &
      beta, c, ldc)

  return
end subroutine dmatmat










subroutine zrmatmatt(m, n, a, k, b, c)
  implicit double precision (a-h,o-z)
  integer *8 :: m,n,k
  complex *16 :: a(n,m),c(k,m)
  real *8 :: b(n,k)
  complex *16, allocatable :: bz(:,:)
  character *1 :: transa, transb
  complex *16 :: alpha, beta
  integer muse,nuse,kuse
  integer *8 nk8
  integer nk

  nk8 = n*k
  nk = n*k
  
  muse = m
  nuse = n
  kuse = k
  allocate(bz(n,k))
  call zzero(nk8,bz)
  call dcopy(nk,b,1,bz,2)

  !
  ! note different dimensions than usual...
  !
  transa = 'T'
  transb = 'N'
  alpha = 1
  beta = 0

  if(m.eq.0.or.n.eq.0.or.k.eq.0) return
  call zgemm(transa, transb, kuse, muse, nuse, alpha, bz, nuse, a, nuse, &
      beta, c, kuse)

  return
end subroutine zrmatmatt
!
!
!
!
!
!
subroutine dgemm_guru(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
  double precision alpha,beta
  integer *8 k,lda,ldb,ldc,m,n
  integer kuse, ldause, ldbuse, ldcuse, muse, nuse
  character transa,transb

  double precision a(lda,*),b(ldb,*),c(ldc,*)

  if(m.eq.0.or.n.eq.0.or.k.eq.0.or.lda.eq.0.or.ldb.eq.0.or.ldc.eq.0) &
    return

  kuse = k
  ldause = lda
  ldbuse = ldb
  ldcuse = ldc
  muse = m
  nuse = n

  call dgemm(transa,transb,muse,nuse,kuse,alpha,a,ldause,b, &
    ldbuse,beta,c,ldcuse)

  return
end subroutine dgemm_guru
!
!
!
!
!
subroutine zgemm_guru(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
  double complex alpha,beta
  integer *8 k,lda,ldb,ldc,m,n
  character transa,transb
  integer kuse,ldause,ldbuse,ldcuse,muse,nuse

  double complex a(lda,*),b(ldb,*),c(ldc,*)

  if(m.eq.0.or.n.eq.0.or.k.eq.0.or.lda.eq.0.or.ldb.eq.0.or.ldc.eq.0) &
    return

  kuse = k
  ldause = lda
  ldbuse = ldb
  ldcuse = ldc
  muse = m
  nuse = n

  call zgemm(transa,transb,muse,nuse,kuse,alpha,a,ldause,b, &
    ldbuse,beta,c,ldcuse)


  return
end subroutine zgemm_guru
!
!
!
!
!
subroutine dgemv_guru(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
  double precision alpha,beta
  integer *8 incx,incy,lda,m,n
  character trans
  integer incxuse,incyuse,ldause,muse,nuse

  double precision a(lda,*),x(*),y(*)

  if(m.eq.0.or.n.eq.0.or.lda.eq.0) return

  incxuse = incx
  incyuse = incy
  ldause = lda
  muse = m
  nuse = n

  call dgemv(trans,muse,nuse,alpha,a,ldause,x,incxuse, &
    beta,y,incyuse)

  return
end subroutine dgemv_guru
!
!
!
!
!
subroutine zgemv_guru(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
  double complex alpha,beta
  integer *8 incx,incy,lda,m,n
  character trans
  integer incxuse,incyuse,ldause,muse,nuse

  double complex a(lda,*),x(*),y(*)

  if(m.eq.0.or.n.eq.0.or.lda.eq.0) return
  incxuse = incx
  incyuse = incy
  ldause = lda
  muse = m
  nuse = n
  call zgemv(trans,muse,nuse,alpha,a,ldause,x,incxuse, &
    beta,y,incyuse)

  return
end subroutine zgemv_guru

