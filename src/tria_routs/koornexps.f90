


!
! Orthgonal polynomials on the simplex (Koornwinder pols)
! (c) 2017 Mike O'Neil
! Contact: oneil@cims.nyu.edu
!
! These codes are similar in nature to those in ortho2eva.f and
! ortho2exps.f, however they *only* compute polynomials on the
! simplex, and not the "standard" equilateral triangle or
! Koornwinders original triangle, 0<u<v<1.
!
! Furthermore, consolidated version of Vioreanu-Rokhlin quadrature
! nodes and weight calculations are contained.
!
!
! Relevant routines for the average user:
!
!   koorn_pols - return a bunch of values of orthonormal pols on the
!       simplex
!
!   koorn_ders - return a bunch of values of orthonormal pols, along
!       with their first and second derivatives on the simplex
!
!   koorn_vals2coefs - construct the matrix mapping values to
!       coefficients in an orthogonal pol expansion (works for
!       arbitrary point distributions)
!
!   koorn_coefs2vals - construct the matrix mapping coefficients to
!       values in an orthogonal pol expansion on the simplex
!
!   koorn_evalexp - evaluates an orthgonal polynomial expansion
!
!   koorn_evalexp2 - evaluates an orthogonal polynomial expansion,
!       along with first and second partial derivatives
!
!   vioreanu_quad -
!
!   koornf_init - initialization routine for faster evaluation of
!     polynomials and ders
!
!   koornf_pols - faster koornwinder polynomial evaluation routine
!      after precomputed coeffcients
!
!   koornf_ders - faster koornwinder polynomial +derivative evaluation 
!      routine after precomputed coeffcients
!
! Library dependencies:
! Library dependencies:
!   prini.f
!   LAPACK
!   lapack_wrap.f90
!
!

subroutine get_vioreanu_wts(norder,npols,wts)
!
!f2py intent(in) norder,npols
!f2py intent(out) wts
!
!  This subroutine extracts the precomputed vioreanu quadrature weights
!  
!  Input arguments
!    
!    - norder: integer *8
!        input order. Must be between 0 and 20
!    - npols: integer *8
!        npols = (norder+1)*(norder+2)/2
!
!  Output arguments
!
!    - wts: double precision(npols)
!        vioreanu quadrature weights of order=norder+1
!
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  integer *8, intent(in) :: norder,npols
  real *8, intent(out) :: wts(npols)

  INCLUDE 'koorn-wts-dat.txt'
  
end subroutine get_vioreanu_wts






subroutine get_vioreanu_nodes(norder,npols,uvs)
!
!f2py intent(in) norder,npols
!f2py intent(out) uvs
!
!  This subroutine extracts the precomputed vioreanu quadrature
!  or discretization nodes
!  
!  Input arguments
!    
!    - norder: integer *8
!        input order. Must be between 0 and 20
!    - npols: integer *8
!        npols = (norder+1)*(norder+2)/2
!
!  Output arguments
!
!    - uvs: double precision(2,npols)
!        vioreanu quadrature nodes of order=norder+1
!
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  integer *8, intent(in) :: norder,npols
  real *8, intent(out) :: uvs(2,npols)

  INCLUDE 'koorn-uvs-dat.txt'
  
end subroutine get_vioreanu_nodes






subroutine get_vioreanu_nodes_wts(norder,npols,uvs,wts)
!
!f2py intent(in) norder,npols
!f2py intent(out) uvs,wts
!
!  This subroutine extracts the precomputed vioreanu quadrature
!  or discretization nodes and the corresponding weights
!  
!  Input arguments
!    
!    - norder: integer *8
!        input order. Must be between 0 and 20
!    - npols: integer *8
!        npols = (norder+1)*(norder+2)/2
!
!  Output arguments
!
!    - uvs: double precision(2,npols)
!        vioreanu quadrature nodes of order=norder+1
!    - wts: double precision(npols)
!        vioreanu quadrature weights of order=norder+1
!
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  integer *8, intent(in) :: norder,npols
  real *8, intent(out) :: uvs(2,npols),wts(npols)

  call get_vioreanu_nodes(norder,npols,uvs)
  call get_vioreanu_wts(norder,npols,wts)
  
end subroutine get_vioreanu_nodes_wts





subroutine vioreanu_simplex_quad(norder, npols, uvs, &
    umatr, vmatr, wts)
!
!f2py intent(in) norder,npols
!f2py intent(out) uvs,wts,umatr,vmatr
!
!  This subroutine extracts the precomputed vioreanu quadrature
!  or discretization nodes, the corresponding weights, the coeffs
!  to values matrix, and the values to coeffs matrix
!  
!  Input arguments
!    
!    - norder: integer *8
!        input order. Must be between 0 and 20
!    - npols: integer *8
!        npols = (norder+1)*(norder+2)/2
!
!  Output arguments
!
!    - uvs: double precision(2,npols)
!        vioreanu quadrature nodes of order=norder+1
!    - umatr: double precision (npols,npols)
!        vals to coeffs matrix
!    - vmatr: double precision (npols,npols)
!        coefs to vals matrix
!    - wts: double precision(npols)
!        vioreanu quadrature weights of order=norder+1
!
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  integer *8, intent(in) :: norder,npols
  real *8, intent(out) :: uvs(2,npols),wts(npols)
  real *8, intent(out) :: umatr(npols,npols),vmatr(npols,npols)


  call get_vioreanu_nodes(norder,npols,uvs)
  call get_vioreanu_wts(norder,npols,wts)
  call koorn_vals2coefs_coefs2vals(norder,npols,umatr,vmatr)


  return
end subroutine vioreanu_simplex_quad




subroutine koorn_pols(uv, nmax, npols, pols)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: uv(2), pols(*)

  real *8 :: legpols(0:100), jacpols(0:100,0:100)

  !
  ! This subroutine evalutes a bunch of orthogonal polynomials on the
  ! simplex with vertices (0,0), (1,0), (0,1). The polynomials
  ! computed by this routine are normalized and rotated Koornwinder
  ! polynomials, which are classically defined on the triangle with
  ! vertices (0,0), (1,0), (1,1), and given analytically as:
  !
  !   K_{n,k}(x,y) = P_{n-k}^(0,2k+1) (2x-1)  x^k  P_k(2y/x-1)
  !
  ! After mapping to the uv-simplex via the change of variables
  !
  !      u = y,  v = 1-x
  !
  ! these polynomials are given as
  !
  !   K_{n,k}(u,v) = P_{n-k}^(0,2k+1) (1-2v)  (1-v)^k  \cdot
  !        P_k((2u+v-1)/(1-v))
  !
  ! See Koornwinder 1975 for details, or the NIST handbook.
  !
  ! Input:
  !   uv - a point in the simplex (0,0), (1,0), (0,1)
  !   nmax - maximum degree polynomials to compute, a total of
  !       (nmax+1)(nmax+2)/2 values are returned
  !
  ! Output:
  !   npols - number of pols = (nmax+1)*(nmax+2)/2
  !   pols - values of all the polynomials, ordered in the following
  !       (n,k) manner:
  !
  !            (0,0), (1,0), (1,1), (2,0), (2,1), (2,2), etc.
  !
  !
  done = 1
  if (nmax .ge. 100) then
    !!!!call prinf('nmax too large, nmax = *', nmax, 1)
    print *, "nmax too large, nmae = ", nmax
    stop
  end if
  

  !
  ! first run the recurrence for P_k(z/y)*y^k
  !
  u = uv(1)
  v = uv(2)
  z = 2*u+v-1
  y = 1-v

  legpols(0) = 1
  legpols(1) = z

  do k=1,nmax
    legpols(k+1) = ((2*k+1)*z*legpols(k) - k*legpols(k-1)*y*y)/(k+1)
  end do

  !!call prin2('legpols = *', legpols, nmax+1)
  !!stop
  

  !
  ! now loop over degrees, in reverse order,
  ! and for each run the jacobi recurrence
  !
  x = 1-2*v
  
  do k = 0,nmax
    beta = 2*k+1    
    jacpols(0,k) = 1
    jacpols(1,k) = (-beta + (2+beta)*x)/2
    
    do n = 1,nmax-k-1
      an = (2*n+beta+1)*(2*n+beta+2)/2/(n+1)/(n+beta+1)
      bn = (-beta**2)*(2*n+beta+1)/2/(n+1)/(n+beta+1)/(2*n+beta)
      cn = n*(n+beta)*(2*n+beta+2)/(n+1)/(n+beta+1)/(2*n+beta)
      jacpols(n+1,k) = (an*x + bn)*jacpols(n,k) - cn*jacpols(n-1,k)
    end do

  end do


  !do k = 0,nmax
  !  call prinf('k = *', k, 1)
  !  call prin2('jacpols = *', jacpols(0,k), nmax-k+1)
  !end do

  !
  ! now assemble the ORTHONORMAL koornwinder polynomials
  !
  iii = 0
  do n = 0,nmax
    do k = 0,n
      sc = sqrt(done/(2*k+1)/(2*n+2))
      iii = iii + 1
      pols(iii) = legpols(k)*jacpols(n-k,k)/sc
    end do
  end do

  npols = iii
  
  return
end subroutine koorn_pols





subroutine koorn_ders(uv, nmax, npols, pols, ders)
!-----------------------
!  This subroutine evaluates the koornwinder polynomials and its
!  derivatives upto a total order nmax
!
!  Input arguments:
!
!    - uv: double precision (2)
!        uv coordinates on $T_{0}$ where polynomials and derivatives
!        are to be evaluated
!
!    - nmax: integer *8
!        total order upto which polynomials and derivatives are to
!        be evaluated
!
!    - npols: integer *8
!        total number of polynomials = (norder+1)*(norder+2)/2
!
!  Output parameters:
!   
!    - pols: double precision(npols)
!        Koornwinder polynomials at uv
!
!    - ders: double precision(2,npols)
!        Derivative of koornwinder polynomails at uv
!---------------------
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  integer *8, intent(in) :: nmax,npols
  real *8, intent(in) :: uv(2)
  real *8, intent(out) :: pols(npols), ders(2,npols)

  real *8 :: legpols(0:100), jacpols(0:100,0:100)
  real *8 :: legu(0:100), legv(0:100)
  real *8 :: jacv(0:100,0:100)
  
  
  !
  ! This subroutine evalutes a bunch of orthogonal polynomials (and
  ! their first partial derivatives) on the simplex with
  ! vertices (0,0), (1,0), (0,1). The polynomials computed by this
  ! routine are normalized and rotated Koornwinder polynomials, which
  ! are classically defined on the triangle with vertices (0,0),
  ! (1,0), (1,1), and given analytically as:
  !
  !   K_{n,k}(x,y) = P_{n-k}^(0,2k+1) (2x-1)  x^k  P_k(2y/x-1)
  !
  ! After mapping to the uv-simplex via the change of variables
  !
  !      u = y,  v = 1-x
  !
  ! these polynomials are given as
  !
  !   K_{n,k}(u,v) = P_{n-k}^(0,2k+1) (1-2v)  (1-v)^k  \cdot
  !        P_k((2u+v-1)/(1-v))
  !
  ! See Koornwinder 1975 for details, or the NIST handbook.
  !
  ! Input:
  !   uv - a point in the simplex (0,0), (1,0), (0,1)
  !   nmax - maximum degree polynomials to compute, a total of
  !       (nmax+1)(nmax+2)/2 values are returned
  !
  ! Output:
  !   pols - values of all the polynomials, ordered in the following
  !       (n,k) manner:
  !
  !            (0,0), (1,0), (0,1), (2,0), (1,1), (0,2), etc.
  !
  !   ders - partial derivatives with respect to u and v
  !
  !


  done = 1
  if (nmax .ge. 100) then
    call prinf('nmax too large, nmax = *', nmax, 1)
    stop
  end if
  

  !
  ! first run the recurrence for P_k(z/y)*y^k
  !
  u = uv(1)
  v = uv(2)
  z = 2*u+v-1
  y = 1-v

  legpols(0) = 1
  legpols(1) = z

  legu(0) = 0
  legu(1) = 2

  legv(0) = 0
  legv(1) = 1

  do k=1,nmax
    legpols(k+1) = ((2*k+1)*z*legpols(k) - k*legpols(k-1)*y*y)/(k+1)
    legu(k+1) = ((2*k+1)*(2*legpols(k) + z*legu(k)) - &
        k*legu(k-1)*y*y)/(k+1)
    legv(k+1) = ((2*k+1)*(legpols(k) + z*legv(k)) - &
        k*(legv(k-1)*y*y - 2*legpols(k-1)*y ))/(k+1)
  end do


  !
  ! now loop over degrees, in reverse order,
  ! and for each run the jacobi recurrence
  !
  x = 1-2*v
  
  do k = 0,nmax
    beta = 2*k+1    
    jacpols(0,k) = 1
    jacpols(1,k) = (-beta + (2+beta)*x)/2

    jacv(0,k) = 0
    jacv(1,k) = -(2+beta)

    do n = 1,nmax-k-1
      an = (2*n+beta+1)*(2*n+beta+2)/2/(n+1)/(n+beta+1)
      bn = (-beta**2)*(2*n+beta+1)/2/(n+1)/(n+beta+1)/(2*n+beta)
      cn = n*(n+beta)*(2*n+beta+2)/(n+1)/(n+beta+1)/(2*n+beta)
      jacpols(n+1,k) = (an*x + bn)*jacpols(n,k) - cn*jacpols(n-1,k)
      jacv(n+1,k) = -2*an*jacpols(n,k) + (an*x + bn)*jacv(n,k) &
          - cn*jacv(n-1,k)
    end do

  end do


  !
  ! now assemble the ORTHONORMAL koornwinder polynomials
  !
  iii = 0
  do n = 0,nmax
    do k = 0,n
      sc = sqrt(done/(2*k+1)/(2*n+2)) 
      iii = iii + 1
      pols(iii) = legpols(k)*jacpols(n-k,k)/sc 
      ders(1,iii) = legu(k)*jacpols(n-k,k)/sc 
      ders(2,iii) = (legv(k)*jacpols(n-k,k) + &
           legpols(k)*jacv(n-k,k))/sc 
    end do
  end do


  return
end subroutine koorn_ders





subroutine koorn_vals2coefs(nmax, npols, uvs, amat)
!
!f2py intent(in) nmax,npols,uvs
!f2py intent(out) amat
!
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  integer *8, intent(in) :: nmax,npols
  real *8, intent(in) :: uvs(2,npols)
  real *8, intent(out) :: amat(npols,npols)

  real *8, allocatable :: bmat(:,:)
  
  !
  ! This routine returns a square matrix that maps point values of a
  ! function to coefficients in an orthonormal polynomial expansion on
  ! the simplex (0,0), (1,0), (0,1). The point locations can be
  ! arbitrary, but note that they WILL affect the conditioning of this
  ! matrix.
  !
  ! Input:
  !   nmax, npols - order of expansion, it should be the case that
  !       npols = (nmax+1)(nmax+2)/2
  !   uvs - point locations, npols of them
  !
  ! Output:
  !   amat - matrix such that coefs = amat*vals
  !
  !

  done = 1

  allocate(bmat(npols,npols))
  call koorn_coefs2vals(nmax, npols, uvs, bmat)
  
  !
  ! now construct its inverse
  !
  call dinverse(npols, bmat, info, amat)
  
  return
end subroutine koorn_vals2coefs





subroutine koorn_coefs2vals(nmax, npols, uvs, amat)
!
!f2py intent(in) nmax,npols,uvs
!f2py intent(out) amat
!
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  integer *8, intent(in) :: nmax,npols
  real *8, intent(in) :: uvs(2,npols)
  real *8, intent(out) :: amat(npols,npols)

  real *8 :: pols(2000)
  
  !
  ! This routine returns a square matrix that maps coefficients in an
  ! orthonormal polynomial expansion on the simplex (0,0), (1,0),
  ! (0,1) to function values at the points uvs. The point locations
  ! can be arbitrary, but note that they WILL affect the conditioning
  ! of this matrix.
  !
  ! Input:
  !   nmax, npols - order of expansion, it should be the case that
  !       npols = (nmax+1)(nmax+2)/2
  !   uvs - point locations, npols of them
  !
  ! Output:
  !   amat - matrix such that vals = amat*coefs
  !
  !

  done = 1

  do i = 1,npols
    call koorn_pols(uvs(1,i), nmax, npols2, pols)
    do j = 1,npols
      amat(i,j) = pols(j)
    end do
  end do
 
  return
end subroutine koorn_coefs2vals




subroutine koorn_vals2coefs_coefs2vals(korder,kpols,umatr,vmatr)
!
!f2py intent(in) korder,kpols
!f2py intent(out) umatr,vmatr
  !
  !!   compute vals2coefs and coefs2vals matrix 
  !
  !    input
  !    korder     in: integer *8
  !                 order of rokhlin vioreanu (rv) nodes at which function is
  !                 sampled
  !    kpols      in: integer *8
  !                 number of nodes corresponding to korder rv nodes
  !
  !    output
  !    umatr      out: real *8 (kpols,kpols)
  !               vals2coefs matrix
  ! 
  !    vmatr      out: real *8 (kpols,kpols)
  !               coefs2vals matrix
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  integer *8, intent(in) :: korder,kpols
  real *8, intent(out) :: umatr(kpols,kpols)
  real *8, intent(out) :: vmatr(kpols,kpols)
  real *8 xys(2,kpols)

  call get_vioreanu_nodes(korder,kpols,xys)
  call koorn_coefs2vals(korder,kpols,xys,vmatr)
  call dinverse(kpols, vmatr, info, umatr)


  return
end subroutine koorn_vals2coefs_coefs2vals
!
!
!
!
!
subroutine koorn_coefs2vals_vioreanu(norder, npols, amat)
!
!f2py intent(in) npols,norder
!f2py intent(out) amat(npols,npols)
!
!
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  integer *8, intent(in) :: norder, npols
  real *8, intent(out) :: amat(npols,npols)
  real *8 :: uvs(2,npols)

  INCLUDE 'koorn-uvs-dat.txt'
  call koorn_coefs2vals(norder, npols, uvs, amat)

  return
end subroutine koorn_coefs2vals_vioreanu
!
!
!
!
!


subroutine koorn_oversamp_mat(korder,kpols,norder,npols,interpmat)
  !
  !!   compute ovesampling matrix
  !
  !    input
  !    korder     in: integer *8
  !                 order of rokhlin vioreanu (rv) nodes at which function is
  !                 sampled
  !    kpols      in: integer *8
  !                 number of nodes corresponding to korder rv nodes
  !    norder     in: integer *8
  !                 order of oversampled rv nodes at which function is
  !                 to be oversampled
  !    npols      in: integer *8
  !                 number of nodes corresponding to norder rv nodes
  !
  !    output
  !    interpmat   out: real *8(npols,kpols)
  !                  upsampling matrix from order korder rv nodes
  !                  to order norder rv nodes
  !    

  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  integer *8 korder,kpols,norder,npols
  real *8 interpmat(npols,kpols)
  real *8 umat(kpols,kpols), vmat(kpols,kpols)
  real *8 xys(2,npols),wts(npols),umato(npols,npols),vmato(npols,npols)
  real *8 pmat(npols,kpols),pols(kpols)

  call koorn_vals2coefs_coefs2vals(korder,kpols,umat,vmat)
  call vioreanu_simplex_quad(norder,npols,xys,umato,vmato,wts)
  
  do i=1,npols
     call koorn_pols(xys(1,i),korder,kpols,pols)
     do j=1,kpols
        pmat(i,j) = pols(j)
     enddo
  enddo

  call dmatmat(npols,kpols,pmat,kpols,umat,interpmat)
  
  return
end subroutine koorn_oversamp_mat




subroutine koorn_evalexp(nmax, npols, uv, coefs, val)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: uv(2), coefs(npols)

  real *8 :: pols(2000)

  !
  ! Evaluate the orthgonal polynomial expansion with given
  ! coefficients at the points uv
  !
  ! Input:
  !   nmax, npols - number of terms in expansion,
  !       npols = (nmax+1)*(nmax+2)/2
  !   uv - point in the simplex at which to evaluate
  !   coefs - coefs of expansion, ordered the same as koorn_pols
  !
  ! Output:
  !   val - the value of the expansion at uv
  !
  !

  call koorn_pols(uv, nmax, npols2, pols)

  val = 0
  do i = 1,npols
    val = val + coefs(i)*pols(i)
  end do

  return
end subroutine koorn_evalexp





subroutine koorn_evalexp2(nmax, npols, uv, coefs, val, der)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: uv(2), coefs(npols), der(2)

  real *8 :: pols(2000), ders(2,2000)

  !
  ! Evaluate the orthgonal polynomial expansion with given
  ! coefficients at the points uv
  !
  ! Input:
  !   nmax, npols - number of terms in expansion,
  !       npols = (nmax+1)*(nmax+2)/2
  !   uv - point in the simplex at which to evaluate
  !   coefs - coefs of expansion, ordered the same as koorn_pols
  !
  ! Output:
  !   val - the value of the expansion at uv
  !   der - the partial derivative at uv
  !

  call koorn_ders(uv, nmax, npols2, pols, ders)

  val = 0
  du = 0
  dv = 0
  duu = 0
  duv = 0
  dvv = 0
  do i = 1,npols
    val = val + coefs(i)*pols(i)
    du = du + coefs(i)*ders(1,i)
    dv = dv + coefs(i)*ders(2,i)
  end do

  der(1) = du
  der(2) = dv
  
  return
end subroutine koorn_evalexp2





subroutine koorn_zevalexp2(nmax, npols, uv, coefs, val, der)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: uv(2)
  complex *16 :: coefs(npols), val, der(2)

  real *8 :: pols(2000), ders(2,2000)
  complex *16 :: du, dv, duu, duv, dvv

  !
  ! Evaluate the orthgonal polynomial expansion with given
  ! coefficients at the points uv
  !
  ! Input:
  !   nmax, npols - number of terms in expansion,
  !       npols = (nmax+1)*(nmax+2)/2
  !   uv - point in the simplex at which to evaluate
  !   coefs - coefs of expansion, ordered the same as koorn_pols
  !
  ! Output:
  !   val - the value of the expansion at uv
  !   der - the partial derivative at uv
  !

  call koorn_ders(uv, nmax, npols2, pols, ders)

  val = 0
  du = 0
  dv = 0
  duu = 0
  duv = 0
  dvv = 0
  do i = 1,npols
    val = val + coefs(i)*pols(i)
    du = du + coefs(i)*ders(1,i)
    dv = dv + coefs(i)*ders(2,i)
  end do

  der(1) = du
  der(2) = dv
  
  return
end subroutine koorn_zevalexp2




























  


!
! faster version for real argument 
!


subroutine koornf_init(nmax, rat1, rat2, rsc1)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
!
! Precompute the recurrence coefficients for the fast
! evaluation of koornwinder functions on triangles and their derivatives
!    
! Parameters:
!   nmax                      must be non-negative
!   rat1(1:2,0:nmax)          recurrence coefficients
!   rat2(1:3,0:nmax,0:nmax)   recurrence coefficients
!   rsc1(0:nmax,0:nmax)       scaling factors
!
  real *8 :: rat1(2,0:nmax)
  real *8 :: rat2(3,0:nmax,0:nmax)
  real *8 :: rsc1(0:nmax,0:nmax)

  do k=0,nmax
    do l=0,nmax
      rsc1(l,k) = 0
      rat2(1,l,k) = 0
      rat2(2,l,k) = 0
      rat2(3,l,k) = 0
    enddo
    rat1(1,k) = 0
    rat1(2,k) = 0
  enddo

  do k=1,nmax
     rat1(1,k) = (2*k+1)/(k+1.0d0)
     rat1(2,k) = -k/(k+1.0d0)
  enddo
  
  do k = 0,nmax
     beta = 2*k+1

     if(nmax.gt.0) then
       rat2(1,1,k) = (2+beta)/2
       rat2(2,1,k) = -beta/2
     endif

     do n = 1,nmax-k-1
        an = (2*n+beta+1)*(2*n+beta+2)/2/(n+1.0d0)/(n+beta+1)
        bn = (-beta**2)*(2*n+beta+1)/2/(n+1.0d0)/(n+beta+1)/(2*n+beta)
        cn = n*(n+beta)*(2*n+beta+2)/(n+1.0d0)/(n+beta+1)/(2*n+beta)
        rat2(1,n+1,k) = an
        rat2(2,n+1,k) = bn
        rat2(3,n+1,k) = -cn
     enddo
  enddo

  done=1
  !
  ! now assemble the ORTHONORMAL koornwinder polynomials
  !
  iii = 0
  do n = 0,nmax
    do k = 0,n
      sc = sqrt(done/(2*k+1)/(2*n+2))
      iii = iii + 1
      rsc1(n-k,k) = 1/sc
    end do
  end do
  
  return
end subroutine koornf_init



subroutine koornf_pols(uv, nmax, npols, pols, rat1, rat2, rsc1)
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8 :: uv(2), pols(*)

  real *8 :: rat1(2,0:nmax)
  real *8 :: rat2(3,0:nmax,0:nmax)
  real *8 :: rsc1(0:nmax,0:nmax)
  real *8 :: legpols(0:100), jacpols(0:100,0:100)

  !
  ! This subroutine evalutes a bunch of orthogonal polynomials on the
  ! simplex with vertices (0,0), (1,0), (0,1). The polynomials
  ! computed by this routine are normalized and rotated Koornwinder
  ! polynomials, which are classically defined on the triangle with
  ! vertices (0,0), (1,0), (1,1), and given analytically as:
  !
  !   K_{n,k}(x,y) = P_{n-k}^(0,2k+1) (2x-1)  x^k  P_k(2y/x-1)
  !
  ! After mapping to the uv-simplex via the change of variables
  !
  !      u = y,  v = 1-x
  !
  ! these polynomials are given as
  !
  !   K_{n,k}(u,v) = P_{n-k}^(0,2k+1) (1-2v)  (1-v)^k  \cdot
  !        P_k((2u+v-1)/(1-v))
  !
  ! See Koornwinder 1975 for details, or the NIST handbook.
  !
  ! Input:
  !   uv - a point in the simplex (0,0), (1,0), (0,1)
  !   nmax - maximum degree polynomials to compute, a total of
  !       (nmax+1)(nmax+2)/2 values are returned
  !
  ! Output:
  !   npols - number of pols = (nmax+1)*(nmax+2)/2
  !   pols - values of all the polynomials, ordered in the following
  !       (n,k) manner:
  !
  !            (0,0), (1,0), (0,1), (2,0), (1,1), (0,2), etc.
  !
  !
  done = 1
  if (nmax .ge. 100) then
     call prinf('nmax too large, nmax = *', nmax, 1)
     stop
  end if

  u = uv(1)
  v = uv(2)

  !
  ! first run the recurrence for P_k(z/y)*y^k
  !
  z = 2*u+v-1
  y = 1-v

  legpols(0) = 1
  legpols(1) = z

  do k=1,nmax-1
     legpols(k+1) = rat1(1,k)*z*legpols(k) + rat1(2,k)*legpols(k-1)*y*y
  end do

  !!call prin2('legpols = *', legpols, nmax+1)
  !!stop
  

  !
  ! now loop over degrees, in reverse order,
  ! and for each run the jacobi recurrence
  !
  x = 1-2*v
  
  do k = 0,nmax
     jacpols(0,k) = 1
     if(nmax.gt.0) jacpols(1,k) = rat2(1,1,k)*x + rat2(2,1,k)

     do n = 1,nmax-k-1
        jacpols(n+1,k) = (rat2(1,n+1,k)*x + rat2(2,n+1,k))*jacpols(n,k) + &
             rat2(3,n+1,k)*jacpols(n-1,k)
     end do

  end do


  !do k = 0,nmax
  !  call prinf('k = *', k, 1)
  !  call prin2('jacpols = *', jacpols(0,k), nmax-k+1)
  !end do

  !
  ! now assemble the ORTHONORMAL koornwinder polynomials
  ! compatible with koorn_pols
  !
  iii = 0
  do n = 0,nmax
    do k = 0,n
      iii = iii + 1
      pols(iii) = legpols(k)*jacpols(n-k,k) * rsc1(n-k,k)
    end do
  end do

  npols = iii

  return
end subroutine koornf_pols







subroutine koornf_ders(uv, nmax, npols, pols, ders, rat1, rat2, rsc1)
!-----------------------
!  This subroutine evaluates the koornwinder polynomials and its
!  derivatives upto a total order nmax, where precomputed
!  coeffs for recurrence relations have been precomputed.
!
!  Note that the factors rat1,rat2, and rsc1 must be precomputed
!  with a call to koornf_init
!
!  Input arguments:
!
!    - uv: double precision (2)
!        uv coordinates on $T_{0}$ where polynomials and derivatives
!        are to be evaluated
!
!    - nmax: integer *8
!        total order upto which polynomials and derivatives are to
!        be evaluated
!
!    - npols: integer *8
!        total number of polynomials = (norder+1)*(norder+2)/2
!
!    - rat1: double precision (2,0:nmax)
!        recurrence coefficients 1 from koornf_init
!        
!    - rat2: double precision (3,0:nmax,0:nmax)
!        recurrence coefficients 2 from koornf_init
!        
!    - rsc1: double precision (3,0:nmax,0:nmax)
!        scaling factors from koornf_init
!
!  Output parameters:
!   
!    - pols: double precision(npols)
!        Koornwinder polynomials at uv
!
!    - ders: double precision(2,npols)
!        Derivative of koornwinder polynomails at uv
!---------------------
!
  implicit real *8 (a-h,o-z)
  implicit integer *8 (i-n)
  real *8, intent(in) :: uv(2)
  integer *8, intent(in) :: nmax
  real *8, intent(in) :: rat1(2,0:nmax)
  real *8, intent(in) :: rat2(3,0:nmax,0:nmax)
  real *8, intent(in) :: rsc1(0:nmax,0:nmax)
  real *8, intent(out) :: pols(npols), ders(2,npols)

  real *8, allocatable :: jacpols(:,:),jacv(:,:)
  real *8 :: legpols(0:100) 
  real *8 :: legu(0:100), legv(0:100)

  allocate(jacpols(0:100,0:100),jacv(0:100,0:100))

  
  
  !
  ! This subroutine evalutes a bunch of orthogonal polynomials (and
  ! their first partial derivatives) on the simplex with
  ! vertices (0,0), (1,0), (0,1). The polynomials computed by this
  ! routine are normalized and rotated Koornwinder polynomials, which
  ! are classically defined on the triangle with vertices (0,0),
  ! (1,0), (1,1), and given analytically as:
  !
  !   K_{n,k}(x,y) = P_{n-k}^(0,2k+1) (2x-1)  x^k  P_k(2y/x-1)
  !
  ! After mapping to the uv-simplex via the change of variables
  !
  !      u = y,  v = 1-x
  !
  ! these polynomials are given as
  !
  !   K_{n,k}(u,v) = P_{n-k}^(0,2k+1) (1-2v)  (1-v)^k  \cdot
  !        P_k((2u+v-1)/(1-v))
  !
  ! See Koornwinder 1975 for details, or the NIST handbook.
  !
  !


  done = 1
  if (nmax .ge. 100) then
    call prinf('nmax too large, nmax = *', nmax, 1)
    stop
  end if
  

  !
  ! first run the recurrence for P_k(z/y)*y^k
  !
  u = uv(1)
  v = uv(2)
  z = 2*u+v-1
  y = 1-v

  legpols(0) = 1
  legpols(1) = z

  legu(0) = 0
  legu(1) = 2

  legv(0) = 0
  legv(1) = 1

  do k=1,nmax
     legpols(k+1) = rat1(1,k)*z*legpols(k) + rat1(2,k)*legpols(k-1)*y*y
     legu(k+1) = (rat1(1,k)*(2*legpols(k) + z*legu(k)) + &
          rat1(2,k)*legu(k-1)*y*y)
     legv(k+1) = (rat1(1,k)*(legpols(k) + z*legv(k)) + &
          rat1(2,k)*(legv(k-1)*y*y - 2*legpols(k-1)*y ))
  end do

  !
  ! now loop over degrees, in reverse order,
  ! and for each run the jacobi recurrence
  !
  x = 1-2*v
  
  do k = 0,nmax
     jacpols(0,k) = 1
     jacpols(1,k) = rat2(1,1,k)*x + rat2(2,1,k)

     jacv(0,k) = 0
     jacv(1,k) = -2*rat2(1,1,k)

     do n = 1,nmax-k-1
        jacpols(n+1,k) = (rat2(1,n+1,k)*x + rat2(2,n+1,k))*jacpols(n,k) + &
             rat2(3,n+1,k)*jacpols(n-1,k)
        jacv(n+1,k) = -2*rat2(1,n+1,k)*jacpols(n,k) + &
             (rat2(1,n+1,k)*x + rat2(2,n+1,k))*jacv(n,k) + &
             rat2(3,n+1,k)*jacv(n-1,k)
     end do
  end do

  !
  ! now assemble the ORTHOGONAL koornwinder polynomials
  ! compatible with ortho2eva3
  !
  iii = 0
  do n = 0,nmax
     do k = 0,n
        iii = iii + 1
        pols(iii) = legpols(k)*jacpols(n-k,k) * rsc1(n-k,k)
        ders(1,iii) = legu(k)*jacpols(n-k,k) * rsc1(n-k,k)
        ders(2,iii) = (legv(k)*jacpols(n-k,k) + &
             legpols(k)*jacv(n-k,k))* rsc1(n-k,k)
     end do
  end do

  npols = iii

  return
end subroutine koornf_ders

