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
!
subroutine dcopy_guru(n,a,incx,b,incy)
  implicit real *8 (a-h,o-z)

  integer n,incx,incy
  integer *8 n1,incx1,incy1

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
  complex *16, intent(in) :: x(n), y(n)
  complex *16, intent(out) :: cd
  complex *16 zdotu
  integer *8 n1,one

  n1 = n
  one = 1

  cd = 0
  cd = zdotu(n1,x,one,y,one)
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
  integer, intent(in) :: m,n
  double complex, intent(inout) :: a(*)
  double complex, allocatable :: at(:)
  integer *8 mn,incx,incy

  allocate(at(m*n))

  do j = 1,n
    do i = 1,m
      at(j+n*(i-1)) = a(i+m*(j-1))
    end do
  end do

  incx = 1
  incy = 1

  mn = m*n

  call zcopy(mn, at, incx, a, incy)

  deallocate(at)
  
  return
end subroutine ztranspose



subroutine zdiff(n, x, y, z)
  implicit double precision (a-h,o-z)
  complex *16 :: x(n), y(n), z(n)

  do i = 1,n
    z(i) = x(i)-y(i)
  end do
  return
end subroutine zdiff





subroutine zratios(n, x, y, z)
  implicit double precision (a-h,o-z)
  complex *16 :: x(n), y(n), z(n)

  do i = 1,n
    z(i) = x(i)/y(i)
  end do
  return
end subroutine zratios






subroutine zeigs(n, a, info, vleft, zlams, vright)
  implicit double precision (a-h,o-z)
  complex *16 :: a(n,n), vleft(n,n), zlams(n), vright(n,n)

  character :: jobvl, jobvr
  double precision, allocatable :: dwork(:)
  complex *16, allocatable :: work(:), atemp(:,:)
  integer *8 lwork,info1,ldvr,ldvl,n1

  n1 = n

  jobvl = 'V'
  jobvr = 'V'
  ldvl = n
  ldvr = n
  lwork = 1000000

  lwork = 10*n
  allocate(work(lwork))
  allocate(dwork(lwork))
  allocate(atemp(n,n))

  do i = 1,n
    do j = 1,n
      atemp(i,j) = a(i,j)
    end do
  end do

  call zgeev(jobvl, jobvr, n1, atemp, n1, zlams, vleft, &
      ldvl, vright, ldvr, work, lwork, dwork, info1) 
  info = info1

  if (info .ne. 0) then
    !call prinf('in zeigs, info = *', info, 1)
    print *, "in zeigs, info = ", info
    return
  endif
  
  return
end subroutine zeigs





subroutine zsvd(m, n, a, u, s, vt)
  implicit double precision (a-h,o-z)
  double precision :: s(*)
  complex *16 :: a(m,n), u(*), vt(*)

  character :: jobu, jobvt
  double precision, allocatable :: dwork(:)
  complex *16, allocatable :: work(:), atemp(:,:)

  integer *8 info,lwork,ldvt,ldu,lda,m1,n1

  m1 = m
  n1 = n
  
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
  
  call zgesvd(jobu, jobvt, m1, n1, atemp, lda, s, u, ldu, &
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
  double precision :: s(*)
  double precision :: a(m,n), u(*), vt(*)

  character :: jobu, jobvt
  double precision, allocatable :: work(:), atemp(:,:)
  integer *8 info,lwork,ldvt,ldu,lda,m1,n1
  
  !
  ! note: V^* is returned, not V
  !
  m1 = m
  n1 = n
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
  
  call dgesvd(jobu, jobvt, m1, n1, atemp, lda, s, u, ldu, &
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
  double precision :: a(n,n), ainv(n,n)
  
  double precision, allocatable :: work(:)
  integer *8, allocatable :: ipiv(:)
  integer *8 incx,incy,n2,n1,info1,lwork

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

  n1 = n
  n2 = n*n
  incx = 1
  incy = 1
  call dcopy(n2, a, incx, ainv, incy)
  allocate(ipiv(n))
  call dgetrf(n1, n1, ainv, n1, ipiv, info1)

  lwork = 10*n
  allocate(work(lwork))
  call dgetri(n1, ainv, n1, ipiv, work, lwork, info1)
  info  = info1

  return
end subroutine dinverse





subroutine zinverse(n, a, info, ainv)
  implicit double precision (a-h,o-z)
  complex *16 :: a(n,n), ainv(n,n)

  integer, allocatable :: ipiv(:)
  complex *16, allocatable :: work(:)
  integer *8 n1,n2,incx,incy,lwork,info1

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
  n1 = n
  n2 = n*n
  incx = 1
  incy = 1
  call zcopy(n2, a, incx, ainv, incy)

  allocate(ipiv(n))
  call zgetrf(n1, n1, ainv, n1, ipiv, info1)

  lwork = 20*n
  allocate(work(lwork))
  call zgetri(n1, ainv, n1, ipiv, work, lwork, info1)
  info = info1

  return
end subroutine zinverse






subroutine dgausselim(n, a, rhs, info, sol, dcond)
  implicit double precision (a-h,o-z)
  double precision :: a(n,n), rhs(n), sol(n)
  
  integer *8, allocatable :: ipiv(:), iwork(:)
  integer *8 info1,n1,lda,ldafldb,ldx,nrhs
  double precision, allocatable :: af(:,:), rscales(:), cscales(:), work(:)
  character *1 :: fact, trans, equed
  !
  ! just a wrapper for the lapack routine...
  ! a is untouched
  !
  
  allocate(ipiv(n))
  allocate(af(n,n))
  allocate(rscales(n))
  allocate(cscales(n))
  allocate(work(5*n))
  allocate(iwork(n))

  n1 = n
  info1 = 0
  
  fact = 'N'
  trans = 'N'
  nrhs = 1
  lda = n
  ldaf = n
  equed = 'N'
  ldb = n
  ldx = n
  
  call dgesvx( fact, trans, n1, nrhs, a, lda, af, ldaf, ipiv, &
      equed, rscales, cscales, rhs, ldb, sol, ldx, dcond, ferr, berr, &
      work, iwork, info1 )

  dcond = 1.0d0/dcond
  info = info1
  
  return
end subroutine dgausselim


subroutine dgausselim_vec(n, a, k, rhs, info, sol, dcond)
  implicit double precision (a-h,o-z)
  double precision :: a(n,n), rhs(n,k), sol(n,k)
  
  integer *8, allocatable :: ipiv(:), iwork(:)
  double precision, allocatable :: af(:,:), rscales(:), cscales(:), work(:)
  integer *8 n1,lda,ldaf,nrhs,ldb,ldx,info1
  character *1 :: fact, trans, equed
  !
  ! This a wrapper for double precision Gaussian elimination in
  ! LAPACK for multipole right hand sides
  ! 
  !
  
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
  n1 = n
  
  
  call dgesvx( fact, trans, n1, nrhs, a, lda, af, ldaf, ipiv, &
      equed, rscales, cscales, rhs, ldb, sol, ldx, dcond, ferr, berr, &
      work, iwork, info1 )

  dcond = 1.0d0/dcond
  info = info1
  
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
  integer m,n,info,nrhs,irank
  real *8 a(m,n), rhs(m,nrhs), sol(n,nrhs), eps
  ! local  
  real *8, allocatable :: work(:), atemp(:,:), rhstemp(:,:)
  integer *8, allocatable :: ipiv(:)
  real *8 dcond

  integer *8 lwork, lda, ldb, i, j, mn, mn2, nb, nb1, nb2, nb3, nb4
  integer *8 lwkmin, lwkopt

  integer *8 ilaenv,info1,m1,n1,irank1

  dcond = 1d0/eps
  dcond = eps


  m1 = m
  n1 = n

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

  call dgelsy(m1,n1,nrhs,atemp,lda,rhstemp,ldb,ipiv,dcond, &
       irank1,work,lwork,info1)

  do i = 1,nrhs
     do j = 1,n
        sol(j,i) = rhstemp(j,i)
     enddo
  enddo

  info = info1
  irank = irank1

  return
end subroutine dleastsq
  


subroutine zgausselim(n, a, rhs, info, sol, dcond)
  implicit double precision (a-h,o-z)
  complex *16 :: a(n,n), rhs(n), sol(n)
  
  integer *8, allocatable :: ipiv(:)
  double precision, allocatable :: rscales(:), cscales(:), rwork(:)
  complex *16, allocatable :: af(:,:), work(:)
  character *1 :: fact, trans, equed
  integer *8 info1,n1,lda,ldf,nrhs,ldb,ldx
  !
  ! just a wrapper for the lapack routine...
  ! a is untouched
  !
  !call prinf('n = *', n, 1)
  !call prin2('amat = *', a, n*n*2)
  !call prin2('rhs = *', rhs, 2*n)
  !stop
  
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

  n1 = n
  
  call zgesvx(fact, trans, n1, nrhs, a, lda, af, ldaf, ipiv, &
      equed, rscales, cscales, rhs, ldb, sol, ldx, dcond, ferr, berr, &
      work, rwork, info1)
  dcond = 1.0d0/dcond
  info = info1

  deallocate(ipiv, af, rscales, cscales, work, rwork)
  
  return
end subroutine zgausselim






subroutine dgausselim_zrhs(n, a, rhs, info, sol, dcond)
  implicit double precision (a-h,o-z)
  double precision :: a(n,n), rhs(2,n), sol(2,n)
  
  integer *8, allocatable :: ipiv(:), iwork(:)
  double precision, allocatable :: af(:,:), rscales(:), cscales(:), work(:)
  double precision, allocatable :: x(:,:), b(:,:)
  character *1 :: fact, trans, equed
  integer *8 n1,info1,nrhs,lda,ldaf,ldb,ldx
  !
  ! just a wrapper for the lapack routine...
  ! a is untouched
  !
  
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
  n1 = n

  do i = 1,n
    b(i,1) = rhs(1,i)
    b(i,2) = rhs(2,i)
  end do
  
  call dgesvx( fact, trans, n1, nrhs, a, lda, af, ldaf, ipiv, &
      equed, rscales, cscales, b, ldb, x, ldx, dcond, ferr, berr, &
      work, iwork, info1 )

  do i = 1,n
    sol(1,i) = x(i,1)
    sol(2,i) = x(i,2)
  end do
  
  dcond = 1.0d0/dcond
  info = info1

  !call prin2('x = *', x, 2*n)
  !stop
  
  return
end subroutine dgausselim_zrhs





subroutine zmatvec(m, n, a, x, y)
  implicit double precision (a-h,o-z)
  complex *16 :: a(m,n), x(n), y(m)
  character *1 :: trans
  complex *16 :: alpha, beta
  integer *8 incx,incy,m1,n1

  trans = 'N'
  alpha = 1
  incx = 1
  beta = 0
  incy = 1
  m1 = m
  n1 = n
  if (m.eq.0.or.n.eq.0) return
  call zgemv(trans, m1, n1, alpha, a, m1, x, incx, beta, y, incy)

  return
end subroutine zmatvec





subroutine dmatvec(m, n, a, x, y)
  implicit double precision (a-h,o-z)
  double precision :: a(m,n), x(n), y(m)
  character *1 :: trans
  double precision :: alpha, beta
  integer *8 m1,n1,incx,incy

  trans = 'N'
  alpha = 1
  incx = 1
  beta = 0
  incy = 1
  m1 = m
  n1 = n
  if (m.eq.0.or.n.eq.0) return
  call dgemv(trans, m1, n1, alpha, a, m1, x, incx, beta, y, incy)

  return
end subroutine dmatvec




subroutine dmatzvec(m, n, a, x, y)
  implicit double precision (a-h,o-z)
  double precision :: a(m,n)
  complex *16 :: x(n), y(m)

  complex *16 :: cd

  do i = 1,m
    cd = 0
    do j = 1,n
      cd = cd + a(i,j)*x(j)
    end do
    y(i) = cd
  end do

  return
end subroutine dmatzvec





subroutine dmatzvec_saved(m, n, a, x, y)
  implicit double precision (a-h,o-z)
  double precision :: a(m,n)
  complex *16 :: x(n), y(m)
  character *1 :: trans
  double precision :: alpha, beta

  double precision, allocatable :: xr(:), xi(:), yr(:), yi(:)
  complex *16 :: ima
  integer *8 incx,incy,m1,n1

  allocate(xr(n), xi(n))
  allocate(yr(m), yi(m))
  ima = (0,1)
  
  trans = 'N'
  alpha = 1
  incx = 1
  beta = 0
  incy = 1
  m1 = m
  n1 = n

  do i = 1,n
    xr(i) = dble(x(i))
    xi(i) = imag(x(i))
  end do

  if(m.eq.0.or.n.eq.0) return
  call dgemv(trans, m1, n1, alpha, a, m1, xr, incx, beta, yr, incy)
  call dgemv(trans, m1, n1, alpha, a, m1, xi, incx, beta, yi, incy)

  do i = 1,m
    y(i) = yr(i) + ima*yi(i)
  end do

  deallocate(xr, xi, yr, yi)
  
  return
end subroutine dmatzvec_saved





subroutine zmatmat(m, n, a, k, b, c)
  implicit double precision (a-h,o-z)
  complex *16 :: a(m,n), b(n,k), c(m,k)
  character *1 :: transa, transb
  complex *16 :: alpha, beta
  integer *8 lda,ldb,ldc,m1,n1,k1

  !
  ! note different dimensions than usual...
  !
  transa = 'N'
  transb = 'N'
  alpha = 1
  beta = 0
  lda = m
  ldb = n
  ldc = m

  m1 = m
  n1 = n
  k1 = k 

  if(m.eq.0.or.n.eq.0.or.k.eq.0) return
  call zgemm(transa, transb, m1, k1, n1, alpha, a, lda, b, ldb, &
      beta, c, ldc)

  return
end subroutine zmatmat





subroutine dmatmat(m, n, a, k, b, c)
  implicit double precision (a-h,o-z)
  double precision :: a(m,n), b(n,k), c(m,k)
  character *1 :: transa, transb
  double precision :: alpha, beta
  integer *8 m1,n1,k1,lda,ldb,ldc

  !
  ! note different dimensions than usual...
  !
  transa = 'N'
  transb = 'N'
  alpha = 1
  beta = 0
  lda = m
  ldb = n
  ldc = m
  m1 = m
  n1 = n
  k1 = k

  if(m.eq.0.or.n.eq.0.or.k.eq.0) return
  call dgemm(transa, transb, m1, k1, n1, alpha, a, lda, b, ldb, &
      beta, c, ldc)

  return
end subroutine dmatmat










subroutine zrmatmatt(m, n, a, k, b, c)
  implicit double precision (a-h,o-z)
  complex *16 :: a(n,m),c(k,m)
  real *8 :: b(n,k)
  complex *16, allocatable :: bz(:,:)
  character *1 :: transa, transb
  complex *16 :: alpha, beta

  integer *8 nk,incx,incy,m1,n1,k1


  allocate(bz(n,k))
  call zzero(n*k,bz)

  nk = n*k
  incx = 1
  incy = 2
  call dcopy(nk,b,incx,bz,incy)

  !
  ! note different dimensions than usual...
  !
  transa = 'T'
  transb = 'N'
  alpha = 1
  beta = 0
  k1 = k
  n1 = n
  m1 = m

  if(m.eq.0.or.n.eq.0.or.k.eq.0) return
  call zgemm(transa, transb, k1, m1, n1, alpha, bz, n1, a, n1, &
      beta, c, k1)

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
  integer k,lda,ldb,ldc,m,n
  integer *8 k1,lda1,ldb1,ldc1,m1,n1
  character transa,transb

  double precision a(lda,*),b(ldb,*),c(ldc,*)

  k1 = k
  lda1 = lda
  ldb1 = ldb
  ldc1 = ldc
  m1 = m
  n1 = n

  if(m.eq.0.or.n.eq.0.or.k.eq.0.or.lda.eq.0.or.ldb.eq.0.or.ldc.eq.0) &
    return
  call dgemm(transa,transb,m1,n1,k1,alpha,a,lda1,b,ldb1,beta,c,ldc1)

  return
end subroutine dgemm_guru
!
!
!
!
!
subroutine zgemm_guru(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
  double complex alpha,beta
  integer k,lda,ldb,ldc,m,n
  integer *8 k1,lda1,ldb1,ldc1,m1,n1
  character transa,transb

  double complex a(lda,*),b(ldb,*),c(ldc,*)

  k1 = k
  lda1 = lda
  ldb1 = ldb
  ldc1 = ldc
  m1 = m
  n1 = n

  if(m.eq.0.or.n.eq.0.or.k.eq.0.or.lda.eq.0.or.ldb.eq.0.or.ldc.eq.0) &
    return
  call zgemm(transa,transb,m1,n1,k1,alpha,a,lda1,b,ldb1,beta,c,ldc1)

  return
end subroutine zgemm_guru
!
!
!
!
!
subroutine dgemv_guru(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
  double precision alpha,beta
  integer incx,incy,lda,m,n
  integer *8 incx1,incy1,lda1,m1,n1
  character trans
  
  double precision a(lda,*),x(*),y(*)

  lda1 = lda
  incx1 = incx
  incy1 = incy
  m1 = m
  n1 = n

  if(m.eq.0.or.n.eq.0.or.lda.eq.0) return
  call dgemv(trans,m1,n1,alpha,a,lda1,x,incx1,beta,y,incy1)

  return
end subroutine dgemv_guru
!
!
!
!
!
subroutine zgemv_guru(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
  double complex alpha,beta
  integer incx,incy,lda,m,n
  integer *8 incx1,incy1,lda1,m1,n1
  character trans
  
  double complex a(lda,*),x(*),y(*)

  lda1 = lda
  incx1 = incx
  incy1 = incy
  m1 = m
  n1 = n

  if(m.eq.0.or.n.eq.0.or.lda.eq.0) return
  call zgemv(trans,m1,n1,alpha,a,lda1,x,incx1,beta,y,incy1)

  return
end subroutine zgemv_guru


