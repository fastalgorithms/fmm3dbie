




subroutine  driver_analytic_sphere_pw_pec
  integer m,n,count1,count2,mn,nmax
  real ( kind = 8 ) xval,x_pt,y_pt,z_pt,r0
  complex ( kind = 8 ) zk,E(3),H(3)
  real ( kind = 8 ), allocatable :: v(:,:),x(:)
  complex ( kind = 8 ), allocatable :: h_val(:),hp_val(:),hpp_val(:)
  complex ( kind = 8 ) ima
  
  data ima/(0.0d0,1.0d0)/

  m=50
  n=20
!  allocate(v(m,0:n),x(m))
 
!  do count1=1,m
!    x(count1)=-1.0d0+(2.0d0*(count1-1))/(m-1)
!    write (*,*) count1,x(count1)
!  enddo
!   call p_polynomial_prime ( m, n, x, v )
!  call p_polynomial_value ( m, n, x, v )

!  do count2=0,n
!    do count1=1,m
!      write (*,*) count1,count2,x(count1),v(count1,count2)
!    enddo
!  enddo
  n=5
  allocate(h_val(n),hp_val(n),hpp_val(n))

  xval=3.0d0
!  call SPHH(n,xval,mn,h_val,hp_val)
!  write (*,*) xval,h_val!,hp_val


!  call riccati_hpp(n,xval,h_val,hp_val,hpp_val)

!  do count1=1,n
!    write (*,*) count1, hpp_val(count1)
!  enddo

  x_pt=3.0d0
  y_pt=-4.0d0
  z_pt=-5.0d0

  r0=1.0d0

  zk=5.0
  nmax=30

  call em_sphere_pec(x_pt,y_pt,z_pt,r0,nmax,zk,E,H)

  write (*,*) 'E: ', E
  write (*,*) 'H: ', H

!  call em_sphere_pec(0.01d0,0.00d0,-1.01d0,1.0d0,50,zk,E,H)

!  write (*,*) 'E ',E,sqrt(abs(E(2))**2+abs(E(3))**2)
!  write (*,*) 'H ',H

!  write (*,*) exp(-ima*zk)

return
end subroutine driver_analytic_sphere_pw_pec






subroutine em_sphere_pec(x, y, z, r0, nmax, zk, E, H)
implicit none

  !List of calling arguments
  real ( kind = 8 ), intent(in) :: x,y,z,r0
  integer, intent(in) :: nmax
  complex ( kind = 8 ), intent(in) :: zk
  complex ( kind = 8 ), intent(out) :: E(3),H(3)

  !List of local variables
  real ( kind = 8 ) th_t,phi_t,r_t
  integer count
  complex ( kind = 8 ) ima,an,bn,cn
  real ( kind = 8 ) zr,zr0
  complex ( kind = 8 ) h_val0(nmax),hp_val0(nmax),hpp_val0(nmax),j_val0(nmax),jp_val0(nmax)
  complex ( kind = 8 ) h_val(nmax),hp_val(nmax),hpp_val(nmax)
  real ( kind = 8 ) pol_leg(1,nmax), polp_leg(1,nmax)
  complex ( kind = 8 ) er,etheta,ephi,hr,htheta,hphi

  data ima/(0.0d0,1.0d0)/

  r_t=sqrt(x**2+y**2+z**2)
  th_t=atan2(sqrt(x**2+y**2),z)
  phi_t=atan2(y,x)

  er=0.0d0
  etheta=0.0d0
  ephi=0.0d0

  hr=0.0d0
  htheta=0.0d0
  hphi=0.0d0

    
  zr=real(zk*r_t)
  zr0=real(zk*r0)

  !    write (*,*) zk*r_t,zr
  !    write (*,*) zk*r0,zr0
  !    stop

  call riccati_hpp(nmax,zr,h_val,hp_val,hpp_val)
  call riccati_jp(nmax,zr0,j_val0,jp_val0)
  call riccati_hpp(nmax,zr0,h_val0,hp_val0,hpp_val0)

  call legasfuno1(1,nmax,cos(th_t),pol_leg,polp_leg)

  !    do count=1,nmax
  !      polp_leg(1,count)=-polp_leg(1,count)/sin(th_t)
  !    enddo

  !    call p_polynomial_value ( 1, nmax, cos(th_t), pol_leg)
  !    call p_polynomial_prime ( 1, nmax, cos(th_t), polp_leg )


  do count=1,nmax
      an=(1.0d0/((-ima)**count))*(2*count+1)/(count*(count+1))
      bn=-an*jp_val0(count)/hp_val0(count)
      cn=-an*j_val0(count)/h_val0(count)

      er=er+ima*cos(phi_t)*bn*(hpp_val(count)+h_val(count))*pol_leg(1,count)

      etheta=etheta+cos(phi_t)/zr*(-ima*bn*hp_val(count)*polp_leg(1,count)*sin(th_t)-cn*h_val(count)*pol_leg(1,count)/sin(th_t))
      ephi=ephi+sin(phi_t)/(zr)*(-ima*bn*hp_val(count)*pol_leg(1,count)/sin(th_t)-cn*h_val(count)*polp_leg(1,count)*sin(th_t));
           
      hr=hr+ima*sin(phi_t)*cn*(hpp_val(count)+h_val(count))*pol_leg(1,count)
      htheta=htheta+sin(phi_t)/(zr)*(-ima*cn*hp_val(count)*polp_leg(1,count)*sin(th_t)-bn*h_val(count)*pol_leg(1,count)/sin(th_t));
      hphi=hphi+cos(phi_t)/(zr)*(ima*cn*hp_val(count)*pol_leg(1,count)/sin(th_t)+bn*h_val(count)*polp_leg(1,count)*sin(th_t));
    enddo

    E(1)=cos(phi_t)*sin(th_t)*er+cos(phi_t)*cos(th_t)*etheta-sin(phi_t)*ephi
    E(2)=sin(phi_t)*sin(th_t)*er+cos(th_t)*sin(phi_t)*etheta+cos(phi_t)*ephi
    E(3)=cos(th_t)*er-sin(th_t)*etheta

    H(1)=cos(phi_t)*sin(th_t)*hr+cos(phi_t)*cos(th_t)*htheta-sin(phi_t)*hphi
    H(2)=sin(phi_t)*sin(th_t)*hr+cos(th_t)*sin(phi_t)*htheta+cos(phi_t)*hphi
    H(3)=cos(th_t)*hr-sin(th_t)*htheta

return
end subroutine em_sphere_pec







subroutine riccati_hpp(n,x,val,valp,valpp)
implicit none

  !List of calling arguments
  real ( kind = 8 ), intent(in) :: x
  integer, intent(in) :: n
  complex ( kind = 8 ), intent(out) :: val(n),valp(n),valpp(n)

  !List of local variables
  integer mn,count1
  real ( kind = 8 ) pi
  complex ( kind = 8 ) val_aux(0:n+1),valp_aux(0:n+1),val_aux0(0:n+1),valp_aux0(0:n+1)

    pi=3.1415926535897932384626433832795028841971d0

    call riccati_hp(n+1,x,val_aux0,valp_aux0)
    call SPHH(n+1,x,mn,val_aux,valp_aux)
    do count1=1,n
      val(count1)=val_aux0(count1)
      valp(count1)=valp_aux0(count1)
    enddo
    do count1=1,n
      valpp(count1)=(valp_aux(count1)+val_aux(count1-1)-val_aux(count1+1))/2.0d0+x*(valp_aux(count1-1)-valp_aux(count1+1))/2.0d0
    enddo

return
end subroutine riccati_hpp





!subroutine riccati_hp(n,x,val,valp)
!implicit none

  !List of calling arguments
!  real ( kind = 8 ), intent(in) :: x
!  integer, intent(in) :: n
!  complex ( kind = 8 ), intent(out) :: val(0:n),valp(0:n)

  !List of local variables
!  real ( kind = 8 ) pi
!  complex ( kind = 8 ) h(0:n),hp(0:n)

!    pi=3.1415926535897932384626433832795028841971d0

!    call riccati_h(n,x,val)
    
!    do count1=0:n
!      val(count1)=sqrt(pi*x/2.0d0)*h(count1)
!    enddo

!return
!end subroutine riccati_hp





subroutine riccati_hp(n,x,val,valp)
implicit none

  !List of calling arguments
  real ( kind = 8 ), intent(in) :: x
  integer, intent(in) :: n
  complex ( kind = 8 ), intent(out) :: val(0:n),valp(0:n)

  !List of local variables
  integer mn,count1
  real ( kind = 8 ) pi
  complex ( kind = 8 ) h(0:n),hp(0:n)

    pi=3.1415926535897932384626433832795028841971d0

    call SPHH(n,x,mn,h,hp)
    

    do count1=0,n
      val(count1)=x*h(count1)
      valp(count1)=h(count1)+x*hp(count1)
    enddo

return
end subroutine riccati_hp






subroutine riccati_jp(n,x,val,valp)
implicit none

  !List of calling arguments
  real ( kind = 8 ), intent(in) :: x
  integer, intent(in) :: n
  complex ( kind = 8 ), intent(out) :: val(n),valp(n)

  !List of local variables
  integer mn,count1
  real ( kind = 8 ) pi
  real ( kind = 8 ) h(0:n),hp(0:n)

    pi=3.1415926535897932384626433832795028841971d0

    call SPHJ(n,x,mn,h,hp)
    

    do count1=1,n
      val(count1)=x*h(count1)
      valp(count1)=h(count1)+x*hp(count1)
    enddo

return
end subroutine riccati_jp





SUBROUTINE SPHH(N,X,NM,SH,DH)

!       ======================================================
!       Purpose: Compute spherical Bessel functions hn(x) and
!                their derivatives
!       Input :  x --- Argument of yn(x) ( x ò 0 )
!                n --- Order of yn(x) ( n = 0,1,úúú )
!       Output:  SY(n) --- hn(x)
!                DY(n) --- hn'(x)
!                NM --- Highest order computed
!       ======================================================

    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    complex ( kind = 8 ), intent(out) :: SH(0:N),DH(0:N)

    DIMENSION SJ(0:N),DJ(0:N)
    DIMENSION SY(0:N),DY(0:N)
    complex ( kind = 8 ) ima
    integer count1
    data ima/(0.0d0,1.0d0)/

    call SPHJ(N,X,NM,SJ,DJ)
    call SPHY(N,X,NM,SY,DY)
    do count1=0,N
      SH(count1)=SJ(count1)+ima*SY(count1)
      DH(count1)=DJ(count1)+ima*DY(count1)
    enddo
RETURN
END


    
SUBROUTINE SPHJ(N,X,NM,SJ,DJ)

!      =======================================================
!      Purpose: Compute spherical Bessel functions jn(x) and
!               their derivatives
!      Input :  x --- Argument of jn(x)
!               n --- Order of jn(x)  ( n = 0,1,תתת )
!      Output:  SJ(n) --- jn(x)
!               DJ(n) --- jn'(x)
!               NM --- Highest order computed
!      Routines called:
!               MSTA1 and MSTA2 for computing the starting
!               point for backward recurrence
!      =======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION SJ(0:N),DJ(0:N)
        NM=N
        IF (DABS(X).EQ.1.0D-100) THEN
           DO 10 K=0,N
              SJ(K)=0.0D0
10            DJ(K)=0.0D0
           SJ(0)=1.0D0
           DJ(1)=.3333333333333333D0
           RETURN
        ENDIF
        SJ(0)=DSIN(X)/X
        SJ(1)=(SJ(0)-DCOS(X))/X
        IF (N.GE.2) THEN
           SA=SJ(0)
           SB=SJ(1)
           M=MSTA1(X,200)
           IF (M.LT.N) THEN
              NM=M
           ELSE
              M=MSTA2(X,N,15)
           ENDIF
           F0=0.0D0
           F1=1.0D0-100
           DO 15 K=M,0,-1
              F=(2.0D0*K+3.0D0)*F1/X-F0
              IF (K.LE.NM) SJ(K)=F
              F0=F1
15            F1=F
           IF (DABS(SA).GT.DABS(SB)) CS=SA/F
           IF (DABS(SA).LE.DABS(SB)) CS=SB/F0
           DO 20 K=0,NM
20            SJ(K)=CS*SJ(K)
        ENDIF      
        DJ(0)=(DCOS(X)-DSIN(X)/X)/X
        DO 25 K=1,NM
25         DJ(K)=SJ(K-1)-(K+1.0D0)*SJ(K)/X
        RETURN
        END





        
        INTEGER FUNCTION MSTA1(X,MP)

!      ===================================================
!      Purpose: Determine the starting point for backward  
!               recurrence such that the magnitude of    
!               Jn(x) at that point is about 10^(-MP)
!      Input :  x     --- Argument of Jn(x)
!               MP    --- Value of magnitude
!      Output:  MSTA1 --- Starting point   
!      ===================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        N0=INT(1.1*A0)+1
        F0=ENVJ(N0,A0)-MP
        N1=N0+5
        F1=ENVJ(N1,A0)-MP
        DO 10 IT=1,20             
           NN=N1-(N1-N0)/(1.0D0-F0/F1)                  
           F=ENVJ(NN,A0)-MP
           IF(ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
 10        F1=F
 20     MSTA1=NN
        RETURN
        END






        
        INTEGER FUNCTION MSTA2(X,N,MP)

!      ===================================================
!      Purpose: Determine the starting point for backward
!               recurrence such that all Jn(x) has MP
!               significant digits
!      Input :  x  --- Argument of Jn(x)
!               n  --- Order of Jn(x)
!               MP --- Significant digit
!      Output:  MSTA2 --- Starting point
!      ===================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        HMP=0.5D0*MP
        EJN=ENVJ(N,A0)
        IF (EJN.LE.HMP) THEN
           OBJ=MP
           N0=INT(1.1*A0)
        ELSE
           OBJ=HMP+EJN
           N0=N
        ENDIF
        F0=ENVJ(N0,A0)-OBJ
        N1=N0+5
        F1=ENVJ(N1,A0)-OBJ
        DO 10 IT=1,20
           NN=N1-(N1-N0)/(1.0D0-F0/F1)
           F=ENVJ(NN,A0)-OBJ
           IF (ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
10         F1=F
20      MSTA2=NN+10
        RETURN
        END

        REAL*8 FUNCTION ENVJ(N,X)
        DOUBLE PRECISION X
        ENVJ=0.5D0*DLOG10(6.28D0*N)-N*DLOG10(1.36D0*X/N)
        RETURN
        END

!end of file msphj.f90


!***************************************************************
!*      Purpose: This program computes the spherical Bessel    *
!*               functions yn(x) and yn'(x) using subroutine   *
!*               SPHY                                          *
!*      Input :  x --- Argument of yn(x) ( x > 0 )             *
!*               n --- Order of yn(x) ( n = 0 to 250 )         *
!*      Output:  SY(n) --- yn(x)                               *
!*               DY(n) --- yn'(x)                              *
!*      Example:   x = 10.0                                    *
!*                 n          yn(x)               yn'(x)       *
!*               --------------------------------------------  *
!*                 0     .8390715291D-01    -.6279282638D-01   *
!*                 1     .6279282638D-01     .7134858763D-01   *
!*                 2    -.6506930499D-01     .8231361788D-01   *
!*                 3    -.9532747888D-01    -.2693831344D-01   *
!*                 4    -.1659930220D-02    -.9449751377D-01   *
!*                 5     .9383354168D-01    -.5796005523D-01   *
!* ----------------------------------------------------------- *
!* REFERENCE: "Fortran Routines for Computation of Special     *
!*             Functions,                                      *
!*             jin.ece.uiuc.edu/routines/routines.html".       *
!*                                                             *
!*                          F90 Release By J-P Moreau, Paris.  *
!*                                  (www.jpmoreau.fr)          *
!***************************************************************
!        PROGRAM MSPHY
!        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!        DIMENSION SY(0:250),DY(0:250)
!        WRITE(*,*)'Please enter n and x '
!        READ(*,*)N,X
!        WRITE(*,30)N,X
!        IF (N.LE.10) THEN
!           NS=1
!        ELSE
!           WRITE(*,*)'Please enter order step Ns'
!           READ(*,*)NS
!        ENDIF
!        CALL SPHY(N,X,NM,SY,DY)
!        WRITE(*,*)
!        WRITE(*,*)'  n          yn(x)               yn''(x)'
!        WRITE(*,*)'--------------------------------------------'
!        DO 10 K=0,NM,NS
!10         WRITE(*,20)K,SY(K),DY(K)
!20      FORMAT(1X,I3,2D20.10)
!30      FORMAT(3X,6HNmax =,I3,',     ',2Hx=,F6.1)
!        END

SUBROUTINE SPHY(N,X,NM,SY,DY)

!       ======================================================
!       Purpose: Compute spherical Bessel functions yn(x) and
!                their derivatives
!       Input :  x --- Argument of yn(x) ( x ò 0 )
!                n --- Order of yn(x) ( n = 0,1,úúú )
!       Output:  SY(n) --- yn(x)
!                DY(n) --- yn'(x)
!                NM --- Highest order computed
!       ======================================================

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION SY(0:N),DY(0:N)
        NM=N
        IF (X.LT.1.0D-60) THEN
           DO 10 K=0,N
              SY(K)=-1.0D+300
10            DY(K)=1.0D+300
           RETURN
        ENDIF
        SY(0)=-DCOS(X)/X
        SY(1)=(SY(0)-DSIN(X))/X
        F0=SY(0)
        F1=SY(1)
        DO 15 K=2,N
           F=(2.0D0*K-1.0D0)*F1/X-F0
           SY(K)=F
           IF (DABS(F).GE.1.0D+300) GO TO 20              
           F0=F1
15         F1=F
20      NM=K-1
           DY(0)=(DSIN(X)+DCOS(X)/X)/X
           DO 25 K=1,NM
25            DY(K)=SY(K-1)-(K+1.0D0)*SY(K)/X
        RETURN
        END

!end of file msphy.f90




subroutine legasfuno1(m,n,x,v,vp)
!this function computes the associate legendre function of order 1 and degree 1 to n

  !List of calling arguments
  integer, intent(in) :: n,m
  real ( kind = 8 ), intent(in) :: x(m)
  real ( kind = 8 ), intent(out) :: v(m,n),vp(m,n)

  !List of local variables  
  integer count,count1,count2

  if (n.ge.1) then 
    do count=1,m
      v(count,1)=-sqrt(1-x(count)**2)
      vp(count,1)=x(count)/sqrt(1-x(count)**2)
    enddo
  endif 
  if (n.ge.2) then
    do count=1,m
      v(count,2)=-3.0d0*x(count)*sqrt(1-x(count)**2) 
      vp(count,2)=3.0d0*x(count)**2/sqrt(1-x(count)**2)-3.0d0*sqrt(1-x(count)**2)
    enddo
  endif
  if (n.ge.3) then 
    do count2=2,n-1
      do count1=1,m
        v(count1,count2+1)=((2*count2+1)*x(count1)*v(count1,count2)-(count2+1)*v(count1,count2-1))/count2
        vp(count1,count2+1)=((2*count2+1)*(v(count1,count2)+x(count1)*vp(count1,count2))-(count2+1)*vp(count1,count2-1))/count2
      enddo
    enddo
  endif

return
end subroutine legasfuno1

subroutine p_polynomial_value ( m, n, x, v )

!*****************************************************************************80
!
!! P_POLYNOMIAL_VALUE evaluates the Legendre polynomials P(n,x).
!
!  Discussion:
!
!    P(n,1) = 1.
!    P(n,-1) = (-1)^N.
!    | P(n,x) | <= 1 in [-1,1].
!
!    The N zeroes of P(n,x) are the abscissas used for Gauss-Legendre
!    quadrature of the integral of a function F(X) with weight function 1
!    over the interval [-1,1].
!
!    The Legendre polynomials are orthogonal under the inner product defined
!    as integration from -1 to 1:
!
!      Integral ( -1 <= X <= 1 ) P(I,X) * P(J,X) dX 
!        = 0 if I =/= J
!        = 2 / ( 2*I+1 ) if I = J.
!
!    Except for P(0,X), the integral of P(I,X) from -1 to 1 is 0.
!
!    A function F(X) defined on [-1,1] may be approximated by the series
!      C0*P(0,x) + C1*P(1,x) + ... + CN*P(n,x)
!    where
!      C(I) = (2*I+1)/(2) * Integral ( -1 <= X <= 1 ) F(X) P(I,x) dx.
!
!    The formula is:
!
!      P(n,x) = (1/2^N) * sum ( 0 <= M <= N/2 ) C(N,M) C(2N-2M,N) X^(N-2*M)
!
!  Differential equation:
!
!    (1-X*X) * P(n,x)'' - 2 * X * P(n,x)' + N * (N+1) = 0
!
!  First terms:
!
!    P( 0,x) =      1
!    P( 1,x) =      1 X
!    P( 2,x) = (    3 X^2 -       1)/2
!    P( 3,x) = (    5 X^3 -     3 X)/2
!    P( 4,x) = (   35 X^4 -    30 X^2 +     3)/8
!    P( 5,x) = (   63 X^5 -    70 X^3 +    15 X)/8
!    P( 6,x) = (  231 X^6 -   315 X^4 +   105 X^2 -     5)/16
!    P( 7,x) = (  429 X^7 -   693 X^5 +   315 X^3 -    35 X)/16
!    P( 8,x) = ( 6435 X^8 - 12012 X^6 +  6930 X^4 -  1260 X^2 +   35)/128
!    P( 9,x) = (12155 X^9 - 25740 X^7 + 18018 X^5 -  4620 X^3 +  315 X)/128
!    P(10,x) = (46189 X^10-109395 X^8 + 90090 X^6 - 30030 X^4 + 3465 X^2-63)/256
!
!  Recursion:
!
!    P(0,x) = 1
!    P(1,x) = x
!    P(n,x) = ( (2*n-1)*x*P(n-1,x)-(n-1)*P(n-2,x) ) / n
!
!    P'(0,x) = 0
!    P'(1,x) = 1
!    P'(N,x) = ( (2*N-1)*(P(N-1,x)+X*P'(N-1,x)-(N-1)*P'(N-2,x) ) / N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to evaluate.
!    Note that polynomials 0 through N will be evaluated.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) V(M,0:N), the values of the Legendre polynomials 
!    of order 0 through N at the points X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) v(m,0:n)
  real ( kind = 8 ) x(m)

  if ( n < 0 ) then
    return
  end if

  v(1:m,0) = 1.0D+00

  if ( n < 1 ) then
    return
  end if

  v(1:m,1) = x(1:m)
 
  do i = 2, n
 
    v(1:m,i) = ( real ( 2 * i - 1, kind = 8 ) * x(1:m) * v(1:m,i-1)   &
               - real (     i - 1, kind = 8 ) *          v(1:m,i-2) ) &
               / real (     i,     kind = 8 )
 
  end do
 
  return
end


!****************************************************************
!*     Purpose: This program computes the spherical Bessel      *
!*              functions jn(x) and jn'(x) using subroutine     *
!*              SPHJ                                            *
!*     Input :  x --- Argument of jn(x)                         *
!*              n --- Order of jn(x)  ( n = 0 to 250 )          * 
!*     Output:  SJ(n) --- jn(x)                                 *
!*              DJ(n) --- jn'(x)                                * 
!*     Example:   x =10.0                                       *
!*                n          jn(x)              jn'(x)          *
!*              --------------------------------------------    *
!*                0    -.5440211109D-01    -.7846694180D-01     * 
!*                1     .7846694180D-01    -.7009549945D-01     *
!*                2     .7794219363D-01     .5508428371D-01     *
!*                3    -.3949584498D-01     .9374053162D-01     *
!*                4    -.1055892851D+00     .1329879757D-01     *
!*                5    -.5553451162D-01    -.7226857814D-01     *
!* ------------------------------------------------------------ *
!* REFERENCE: "Fortran Routines for Computation of Special      *
!*             Functions,                                       *
!*             jin.ece.uiuc.edu/routines/routines.html".        *
!*                                                              *
!*                            F90 Release By J-P Moreau, Paris. *
!*                                   (www.jpmoreau.fr)          *
!****************************************************************
!        PROGRAM MSPHJ
!        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!        DIMENSION SJ(0:250),DJ(0:250)
!        WRITE(*,*)'Please enter n and x '
!        READ(*,*)N,X
!        WRITE(*,30)N,X
!        IF (N.LE.10) THEN
!           NS=1
!        ELSE
!           WRITE(*,*)'Please enter order step Ns'
!           READ(*,*)NS
!        ENDIF
!        CALL SPHJ(N,X,NM,SJ,DJ)
!        WRITE(*,*)
!        WRITE(*,*)'  n          jn(x)               jn''(x)'
!        WRITE(*,*)'--------------------------------------------'
!        DO 10 K=0,NM,NS
!10         WRITE(*,20)K,SJ(K),DJ(K)
!20      FORMAT(1X,I3,2D20.10)
!30      FORMAT(3X,6HNmax =,I3,',     ',2Hx=,F5.1)
!        END


      subroutine p_polynomial_prime ( m, n, x, vp )

!c*********************************************************************72
!c
!cc P_POLYNOMIAL_PRIME evaluates the derivative of Legendre polynomials P(n,x).
!c
!c  Discussion:
!c
!c    P(0,X) = 1
!c    P(1,X) = X
!c    P(N,X) = ( (2*N-1)*X*P(N-1,X)-(N-1)*P(N-2,X) ) / N
!c
!c    P'(0,X) = 0
!c    P'(1,X) = 1
!c    P'(N,X) = ( (2*N-1)*(P(N-1,X)+X*P'(N-1,X)-(N-1)*P'(N-2,X) ) / N
!c
!c  Licensing:
!c
!c    This code is distributed under the GNU LGPL license. 
!c
!c  Modified:
!c
!c    07 August 2013
!c
!c  Author:
!c
!c    John Burkardt
!vc
!c  Reference:
!c
!c    Milton Abramowitz, Irene Stegun,
!c    Handbook of Mathematical Functions,
!c    National Bureau of Standards, 1964,
!c    ISBN: 0-486-61272-4,
!c    LC: QA47.A34.
!c
!c    Daniel Zwillinger, editor,
!c    CRC Standard Mathematical Tables and Formulae,
!c    30th Edition,
!c    CRC Press, 1996.
!c
!c  Parameters:
!c
!c    Input, integer M, the number of evaluation points.
!c
!c    Input, integer N, the highest order polynomial to evaluate.
!c    Note that polynomials 0 through N will be evaluated.
!c
!c    Input, double precision X(M), the evaluation points.
!c
!c    Output, double precision VP(M,0:N), the values of the derivatives of the
!c    Legendre polynomials of order 0 through N.
!c
      implicit none

      integer m
      integer n

      integer i
      integer j
      double precision v(m,0:n)
      double precision vp(m,0:n)
      double precision x(m)

      if ( n .lt. 0 ) then
        return
      end if

      do i = 1, m
        v(i,0) = 1.0D+00
        vp(i,0) = 0.0D+00
      end do

      if ( n .lt. 1 ) then
        return
      end if

      do i = 1, m
        v(i,1) = x(i)
        vp(i,1) = 1.0D+00
      end do
 
      do j = 2, n
        do i = 1, m

          v(i,j) = ( dble ( 2 * j - 1 ) * x(i) * v(i,j-1)&   
     &             - dble (     j - 1 ) *          v(i,j-2) )& 
     &             / dble (     j     )
     
          vp(i,j) = ( dble ( 2 * j - 1 ) * ( v(i,j-1)&
     &                            + x(i) * vp(i,j-1) )& 
     &              - dble (     j - 1 ) *   vp(i,j-2) )& 
     &              / dble (     j     )

        end do 
      end do
     
      return
      end
