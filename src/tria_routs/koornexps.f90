


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
! Library dependencies:
!   prini.f
!   LAPACK
!   lapack_wrap.f90
!
!

subroutine get_vioreanu_wts(norder,npols,wts)
  implicit real *8 (a-h,o-z)
  real *8 wts(npols),xsout(npols),ysout(npols)
  
  call ortho2smexps(itype,norder,npols,xsout,ysout,wts)
  do i=1,npols
    wts(i) = wts(i)/sqrt(3.0d0)/2
  enddo
end subroutine get_vioreanu_wts






subroutine get_vioreanu_nodes(norder,npols,xys)
  implicit real *8 (a-h,o-z)
  real *8 wts(npols),xsout(npols),ysout(npols)
  real *8 xys(2,npols)
  
  call ortho2smexps(itype,norder,npols,xsout,ysout,wts)
  do i=1,npols
    call ortho2_stdtosimplex(xsout(i),ysout(i), xys(1,i), xys(2,i))
  enddo
end subroutine get_vioreanu_nodes




subroutine get_vioreanu_nodes_wts(norder,npols,xys,wts)
  implicit real *8 (a-h,o-z)
  real *8 wts(npols),xsout(npols),ysout(npols)
  real *8 xys(2,npols)
  
  call ortho2smexps(itype,norder,npols,xsout,ysout,wts)
  do i=1,npols
    call ortho2_stdtosimplex(xsout(i),ysout(i), xys(1,i), xys(2,i))
    wts(i) = wts(i)/sqrt(3.0d0)/2
  enddo
end subroutine get_vioreanu_nodes_wts







subroutine vioreanu_simplex_quad(norder, npols, xys, &
    umatr, vmatr, whts)
  implicit real *8 (a-h,o-z)
  real *8 :: xys(2,*), whts(*)
  dimension xsout(1000), ysout(1000), umatr(1),vmatr(1)
  !
  ! This subroutine constructs the Vioreanu-Rokhlin quadrature nodes
  ! and weights for polynomials on the simplex triangle with the
  ! vertices
  !
  !       (0,0), (1,0), (0,1)                                   (1)
  ! 
  ! It also constructs the matrix vmatr converting the coefficients of
  ! the polynomial expansion to its values at the interpolation nodes,
  ! and its inverse umatr, converting the values of a function at the
  ! interpolation nodes into the coefficients of the polynomial
  ! expansion.
  ! 
  !          Input parameters:
  ! 
  !     norder - the order of the quadrature
  ! 
  ! 
  !          Output parameters:
  ! 
  !     npols - the number of the polynomials, i.e. the terms in the
  !           expansion
  !     xsout - the x-coordinates of the interpolation nodes (npols of them)
  !     ysout - the y-coordinates of the interpolation nodes (npols of them)
  !     umatr - the npols*npols matrix, converting the values of the polynomial
  !           at the interpolation nodes into the coefficients of its expansion
  !     vmatr - the npols*npols matrix, converting the coefficients of 
  !           the expansion into its values at the interpolation nodes.
  !     wsout - the corresponding quadrature weigths (npols of them)
  ! 
  ! c
  ! c
  ! c       ... first, get the interpolation nodes for the standard triangle
  ! c
  ! c
  call ortho2smexps(itype,norder,npols,xsout,ysout,whts)

!  print *, "constructed vioreanu quad, umatr and vmatr are wrong!!"
!  print *, "please update using koorn_pols"
  
  !
  ! ... and convert them into simplex coordinates
  !
  do i=1,npols
    call ortho2_stdtosimplex(xsout(i),ysout(i), xys(1,i), xys(2,i))
    whts(i) = whts(i)/sqrt(3.0d0)
    whts(i) = whts(i)*0.5d0
  end do


  !
  !!  temp fix to umatr and matr via calling vals2coefs and coefs2vals
  !

  call koorn_vals2coefs(norder, npols, xys, umatr)
  call koorn_coefs2vals(norder, npols, xys, vmatr)


  return
end subroutine vioreanu_simplex_quad




subroutine koorn_pols(uv, nmax, npols, pols)
  implicit real *8 (a-h,o-z)
  real *8 :: uv(2), pols(*)

  real *8 :: legpols(0:100), jacpols(0:100,0:100)

  !
  ! This subroutine evalutes a bunch of orthogonal polynomials on the
  ! simplex with vertices (0,0), (1,0), (0,1). The polynomials
  ! computed by this routine are normalized and rotated Koornwinder
  ! polynomials, which are classically defined on the triangle with
  ! vertices (0,0), (1,0), (1,1), and given analytically as:
  !
  !   K_{n,k}(x,y) = P_{n-k}^(0,2k+1) (2x+1)  x^k  P_k(2y/x-1)
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





subroutine koorn_ders(uv, nmax, npols, pols, ders, der2s)
  implicit real *8 (a-h,o-z)
  real *8 :: uv(2), pols(*), ders(2,*), der2s(3,*)

  real *8 :: legpols(0:100), jacpols(0:100,0:100)
  real *8 :: legu(0:100), legv(0:100), leguu(0:100)
  real *8 :: leguv(0:100), legvv(0:100)
  real *8 :: jacv(0:100,0:100), jacvv(0:100,0:100)
  
  
  !
  ! This subroutine evalutes a bunch of orthogonal polynomials (and
  ! their first and second partial derivatives) on the simplex with
  ! vertices (0,0), (1,0), (0,1). The polynomials computed by this
  ! routine are normalized and rotated Koornwinder polynomials, which
  ! are classically defined on the triangle with vertices (0,0),
  ! (1,0), (1,1), and given analytically as:
  !
  !   K_{n,k}(x,y) = P_{n-k}^(0,2k+1) (2x+1)  x^k  P_k(2y/x-1)
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
  !            (0,0), (1,0), (1,1), (2,0), (2,1), (2,2), etc.
  !
  !   ders - partial derivatives with respect to u and v
  !   der2s - second partial derivatives, uu, uv, vv
  !
  !


  done = 1
  if (nmax .ge. 100) then
    print *, "nmax too large, nmae = ", nmax
    !!!!call prinf('nmax too large, nmax = *', nmax, 1)
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

  leguu(0) = 0
  leguu(1) = 0
  
  leguv(0) = 0
  leguv(1) = 0
  
  legvv(0) = 0
  legvv(1) = 0

  do k=1,nmax
    legpols(k+1) = ((2*k+1)*z*legpols(k) - k*legpols(k-1)*y*y)/(k+1)
    legu(k+1) = ((2*k+1)*(2*legpols(k) + z*legu(k)) - &
        k*legu(k-1)*y*y)/(k+1)
    legv(k+1) = ((2*k+1)*(legpols(k) + z*legv(k)) - &
        k*(legv(k-1)*y*y - 2*legpols(k-1)*y ))/(k+1)
    leguu(k+1) = ((2*k+1)*(4*legu(k) + z*leguu(k)) - &
        k*leguu(k-1)*y*y)/(k+1)
    leguv(k+1) = ((2*k+1)*(2*legv(k) + legu(k) + z*leguv(k)) - &
        k*(leguv(k-1)*y*y - 2*y*legu(k-1) ))/(k+1)
    legvv(k+1) = ((2*k+1)*(2*legv(k) + z*legvv(k)) - &
        k*(legvv(k-1)*y*y - 4*y*legv(k-1) + 2*legpols(k-1) ))/(k+1)
  end do


  ! !
  ! ! temporarily return these functions to test the partial deriviates
  ! ! via finite differences
  ! !
  ! npols = nmax+1
  ! do i = 1,nmax+1
  !   pols(i) = legpols(i-1)
  !   ders(1,i) = legu(i-1)
  !   ders(2,i) = legv(i-1)
  !   der2s(1,i) = leguu(i-1)
  !   der2s(2,i) = leguv(i-1)
  !   der2s(3,i) = legvv(i-1)
  ! end do

  ! return
  
  
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

    jacv(0,k) = 0
    jacv(1,k) = -(2+beta)

    jacvv(0,k) = 0
    jacvv(1,k) = 0

    do n = 1,nmax-k-1
      an = (2*n+beta+1)*(2*n+beta+2)/2/(n+1)/(n+beta+1)
      bn = (-beta**2)*(2*n+beta+1)/2/(n+1)/(n+beta+1)/(2*n+beta)
      cn = n*(n+beta)*(2*n+beta+2)/(n+1)/(n+beta+1)/(2*n+beta)
      jacpols(n+1,k) = (an*x + bn)*jacpols(n,k) - cn*jacpols(n-1,k)
      jacv(n+1,k) = -2*an*jacpols(n,k) + (an*x + bn)*jacv(n,k) &
          - cn*jacv(n-1,k)
      jacvv(n+1,k) = -4*an*jacv(n,k) + (an*x + bn)*jacvv(n,k) &
          - cn*jacvv(n-1,k)
    end do

  end do



  ! !
  ! ! temporarily return these functions to test the partial deriviates
  ! ! via finite differences
  ! !
  ! i = 0
  ! do k = 0,nmax
  !   do n = 0,nmax-k
  !     i = i + 1
  !     pols(i) = jacpols(n,k)
  !     ders(1,i) = 0
  !     ders(2,i) = jacv(n,k)
  !     der2s(1,i) = 0
  !     der2s(2,i) = 0
  !     der2s(3,i) = jacvv(n,k)
  !   end do
  ! end do

  ! npols = i
  ! return


  

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
      ders(1,iii) = legu(k)*jacpols(n-k,k)/sc
      ders(2,iii) = (legv(k)*jacpols(n-k,k) &
          + legpols(k)*jacv(n-k,k))/sc
      der2s(1,iii) = leguu(k)*jacpols(n-k,k)/sc
      der2s(2,iii) = (leguv(k)*jacpols(n-k,k) &
          + legu(k)*jacv(n-k,k))/sc
      der2s(3,iii) = (legvv(k)*jacpols(n-k,k) &
          + 2*legv(k)*jacv(n-k,k) + legpols(k)*jacvv(n-k,k))/sc
    end do
  end do

  npols = iii

  return
end subroutine koorn_ders





subroutine koorn_vals2coefs(nmax, npols, uvs, amat)
  implicit real *8 (a-h,o-z)
  real *8 :: uvs(2,npols), amat(npols,npols)

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
  implicit real *8 (a-h,o-z)
  real *8 :: uvs(2,npols), amat(npols,npols)

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
  !!   compute vals2coefs and coefs2vals matrix 
  !
  !    input
  !    korder     in: integer
  !                 order of rokhlin vioreanu (rv) nodes at which function is
  !                 sampled
  !    kpols      in: integer
  !                 number of nodes corresponding to korder rv nodes
  !
  !    output
  !    umatr      out: real *8 (kpols,kpols)
  !               vals2coefs matrix
  ! 
  !    vmatr      out: real *8 (kpols,kpols)
  !               coefs2vals matrix
  implicit real *8 (a-h,o-z)
  integer korder
  real *8 umatr(kpols,kpols)
  real *8 vmatr(kpols,kpols)
  real *8 xys(2,kpols)
  real *8 wts(kpols)

  call vioreanu_simplex_quad(korder, kpols, xys, umatr,vmatr,wts)
  call koorn_coefs2vals(korder,kpols,xys,vmatr)
  call dinverse(kpols, vmatr, info, umatr)


  return
end subroutine koorn_vals2coefs_coefs2vals




subroutine koorn_oversamp_mat(korder,kpols,norder,npols,interpmat)
  !
  !!   compute ovesampling matrix
  !
  !    input
  !    korder     in: integer
  !                 order of rokhlin vioreanu (rv) nodes at which function is
  !                 sampled
  !    kpols      in: integer
  !                 number of nodes corresponding to korder rv nodes
  !    norder     in: integer
  !                 order of oversampled rv nodes at which function is
  !                 to be oversampled
  !    npols      in: integer
  !                 number of nodes corresponding to norder rv nodes
  !
  !    output
  !    interpmat   out: real *8(npols,kpols)
  !                  upsampling matrix from order korder rv nodes
  !                  to order norder rv nodes
  !    

  implicit real *8 (a-h,o-z)
  integer korder,kpols,norder,npols
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





subroutine koorn_evalexp2(nmax, npols, uv, coefs, val, der, der2)
  implicit real *8 (a-h,o-z)
  real *8 :: uv(2), coefs(npols), der(2), der2(3)

  real *8 :: pols(2000), ders(2,2000), der2s(3,2000)

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
  !   der2 - the second partial derivatives at uv
  !

  call koorn_ders(uv, nmax, npols2, pols, ders, der2s)

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
    duu = duu + coefs(i)*der2s(1,i)
    duv = duv + coefs(i)*der2s(2,i)
    dvv = dvv + coefs(i)*der2s(3,i)
  end do

  der(1) = du
  der(2) = dv
  der2(1) = duu
  der2(2) = duv
  der2(3) = dvv
  
  return
end subroutine koorn_evalexp2





subroutine koorn_zevalexp2(nmax, npols, uv, coefs, val, der, der2)
  implicit real *8 (a-h,o-z)
  real *8 :: uv(2)
  complex *16 :: coefs(npols), val, der(2), der2(3)

  real *8 :: pols(2000), ders(2,2000), der2s(3,2000)
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
  !   der2 - the second partial derivatives at uv
  !

  call koorn_ders(uv, nmax, npols2, pols, ders, der2s)

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
    duu = duu + coefs(i)*der2s(1,i)
    duv = duv + coefs(i)*der2s(2,i)
    dvv = dvv + coefs(i)*der2s(3,i)
  end do

  der(1) = du
  der(2) = dv
  der2(1) = duu
  der2(2) = duv
  der2(3) = dvv
  
  return
end subroutine koorn_zevalexp2




























  


!
! faster version for real argument 
!


subroutine koornf_init(nmax, rat1, rat2, rsc1)
  implicit real *8 (a-h,o-z)
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

  do k=1,nmax
     rat1(1,k) = (2*k+1)/(k+1.0d0)
     rat1(2,k) = -k/(k+1.0d0)
  enddo
  
  do k = 0,nmax
     beta = 2*k+1
     rat2(1,1,k) = (2+beta)/2
     rat2(2,1,k) = -beta/2

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
  !   K_{n,k}(x,y) = P_{n-k}^(0,2k+1) (2x+1)  x^k  P_k(2y/x-1)
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
     jacpols(1,k) = rat2(1,1,k)*x + rat2(2,1,k)

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

