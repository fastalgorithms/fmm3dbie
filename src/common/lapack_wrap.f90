! (c) Mike O'Neil
!
! A few interfaces by Manas Rachh
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

  cd = 0
  cd = zdotu(n,x,1,y,1)
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

  allocate(at(m*n))

  do j = 1,n
    do i = 1,m
      at(j+n*(i-1)) = a(i+m*(j-1))
    end do
  end do

  incx = 1
  incy = 1
  call zcopy(m*n, at, incx, a, incy)

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

  call zgeev(jobvl, jobvr, n, atemp, n, zlams, vleft, &
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
  double precision :: s(*)
  complex *16 :: a(m,n), u(*), vt(*)

  character :: jobu, jobvt
  double precision, allocatable :: dwork(:)
  complex *16, allocatable :: work(:), atemp(:,:)
  
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
  
  call zgesvd(jobu, jobvt, m, n, atemp, lda, s, u, ldu, &
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
  
  call dgesvd(jobu, jobvt, m, n, atemp, lda, s, u, ldu, &
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
  
  double precision, allocatable :: ipiv(:), work(:)

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
  incx = 1
  incy = 1
  call dcopy(n*n, a, incx, ainv, incy)
  allocate(ipiv(n))
  call dgetrf(n, n, ainv, n, ipiv, info)

  lwork = 10*n
  allocate(work(lwork))
  call dgetri(n, ainv, n, ipiv, work, lwork, info)

  return
end subroutine dinverse





subroutine zinverse(n, a, info, ainv)
  implicit double precision (a-h,o-z)
  complex *16 :: a(n,n), ainv(n,n)

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
  incx = 1
  incy = 1
  call zcopy(n*n, a, incx, ainv, incy)

  allocate(ipiv(n))
  call zgetrf(n, n, ainv, n, ipiv, info)

  lwork = 20*n
  allocate(work(lwork))
  call zgetri(n, ainv, n, ipiv, work, lwork, info)

  return
end subroutine zinverse





subroutine dgausselim(n, a, rhs, info, sol, dcond)
  implicit double precision (a-h,o-z)
  double precision :: a(n,n), rhs(n), sol(n)
  
  integer, allocatable :: ipiv(:), iwork(:)
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
  
  fact = 'N'
  trans = 'N'
  nrhs = 1
  lda = n
  ldaf = n
  equed = 'N'
  ldb = n
  ldx = n
  
  call dgesvx( fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, &
      equed, rscales, cscales, rhs, ldb, sol, ldx, dcond, ferr, berr, &
      work, iwork, info )

  dcond = 1/dcond
  
  return
end subroutine dgausselim


subroutine dgausselim_vec(n, a, k, rhs, info, sol, dcond)
  implicit double precision (a-h,o-z)
  double precision :: a(n,n), rhs(n,k), sol(n,k)
  
  integer, allocatable :: ipiv(:), iwork(:)
  double precision, allocatable :: af(:,:), rscales(:), cscales(:), work(:)
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
  
  call dgesvx( fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, &
      equed, rscales, cscales, rhs, ldb, sol, ldx, dcond, ferr, berr, &
      work, iwork, info )

  dcond = 1/dcond
  
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
  integer, allocatable :: ipiv(:)
  real *8 dcond

  integer lwork, lda, ldb, i, j, mn, mn2, nb, nb1, nb2, nb3, nb4
  integer lwkmin, lwkopt

  integer ilaenv

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

  call dgelsy(m,n,nrhs,atemp,lda,rhstemp,ldb,ipiv,dcond, &
       irank,work,lwork,info)

  do i = 1,nrhs
     do j = 1,n
        sol(j,i) = rhstemp(j,i)
     enddo
  enddo

  return
end subroutine dleastsq
  


subroutine zgausselim(n, a, rhs, info, sol, dcond)
  implicit double precision (a-h,o-z)
  complex *16 :: a(n,n), rhs(n), sol(n)
  
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
  
  call zgesvx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, &
      equed, rscales, cscales, rhs, ldb, sol, ldx, dcond, ferr, berr, &
      work, rwork, info)
  dcond = 1/dcond

  deallocate(ipiv, af, rscales, cscales, work, rwork)
  
  return
end subroutine zgausselim






subroutine dgausselim_zrhs(n, a, rhs, info, sol, dcond)
  implicit double precision (a-h,o-z)
  double precision :: a(n,n), rhs(2,n), sol(2,n)
  
  integer, allocatable :: ipiv(:), iwork(:)
  double precision, allocatable :: af(:,:), rscales(:), cscales(:), work(:)
  double precision, allocatable :: x(:,:), b(:,:)
  character *1 :: fact, trans, equed
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

  do i = 1,n
    b(i,1) = rhs(1,i)
    b(i,2) = rhs(2,i)
  end do
  
  call dgesvx( fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, &
      equed, rscales, cscales, b, ldb, x, ldx, dcond, ferr, berr, &
      work, iwork, info )

  do i = 1,n
    sol(1,i) = x(i,1)
    sol(2,i) = x(i,2)
  end do
  
  dcond = 1/dcond

  !call prin2('x = *', x, 2*n)
  !stop
  
  return
end subroutine dgausselim_zrhs





subroutine zmatvec(m, n, a, x, y)
  implicit double precision (a-h,o-z)
  complex *16 :: a(m,n), x(n), y(m)
  character *1 :: trans
  complex *16 :: alpha, beta

  trans = 'N'
  alpha = 1
  incx = 1
  beta = 0
  incy = 1
  call zgemv(trans, m, n, alpha, a, m, x, incx, beta, y, incy)

  return
end subroutine zmatvec





subroutine dmatvec(m, n, a, x, y)
  implicit double precision (a-h,o-z)
  double precision :: a(m,n), x(n), y(m)
  character *1 :: trans
  double precision :: alpha, beta

  trans = 'N'
  alpha = 1
  incx = 1
  beta = 0
  incy = 1
  call dgemv(trans, m, n, alpha, a, m, x, incx, beta, y, incy)

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

  call dgemv(trans, m, n, alpha, a, m, xr, incx, beta, yr, incy)
  call dgemv(trans, m, n, alpha, a, m, xi, incx, beta, yi, incy)

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

  call zgemm(transa, transb, m, k, n, alpha, a, lda, b, ldb, &
      beta, c, ldc)

  return
end subroutine zmatmat





subroutine dmatmat(m, n, a, k, b, c)
  implicit double precision (a-h,o-z)
  double precision :: a(m,n), b(n,k), c(m,k)
  character *1 :: transa, transb
  double precision :: alpha, beta

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

  call dgemm(transa, transb, m, k, n, alpha, a, lda, b, ldb, &
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


  allocate(bz(n,k))
  call zzero(n*k,bz)
  call dcopy(n*k,b,1,bz,2)

  !
  ! note different dimensions than usual...
  !
  transa = 'T'
  transb = 'N'
  alpha = 1
  beta = 0

  call zgemm(transa, transb, k, m, n, alpha, bz, n, a, n, &
      beta, c, k)

  return
end subroutine zrmatmatt


