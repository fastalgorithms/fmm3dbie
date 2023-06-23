c
c       This file contains code for generating quadratures for the
c       evaluation of integrals of the form
c
c           \int    K(x,y) \sigma(y) dS(y)                                (1)
c               S
c       where:
c
c           S is a surface element and dS is the corresponding
c           surface measure,
c
c           \sigma is a smooth function.
c           
c           K(x,y) is the double or single layer potential on S, and
c 
c           x is a specified point in S.
c
c       It is assumed that the surface element S is the image under
c       a parameterization
c
c          p: T \to \mathbb{R}^3
c
c       given over a triangle T.  
c
c       The behavior of the Jacobian dp of the parameterization
c       p at the point x has a large influence on the form of the 
c       integrand of (1).  This code proceeds by first composing the 
c       given parameterization p with an appropriate linear mapping 
c
c                 A: \mathbb{R}^2 \to \mathbb{R}^2
c
c       in order to  form a new parameterization p' = p \ocirc A such 
c       that the Jacobian of p' at the point x is conformal.  Then
c       the quadrature rules from radial.f are used to evaluate
c       the resulting integral.
c
c       The quadratures returned by raddiag are formed using precomputed
c       quadrature tables stored on the disk.  These precomputed
c       tables determine the possible orders for the polynomials p
c       and q in (2).  See radial_init for a list of possible orders.
c
c       The following subroutines are user-callable:
c
c
c   self_quadrature_new - 
c       return a quadrature for evaluating an integral of
c       the form (1) for compact kernels
c
c   pv_quadrature_new
c       return a quadrature for evaluating an integral of
c       the form (1) for pv kernels
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine self_quadrature_new(irad, verts0, nv,x0, y0, dr,
     1    nquad, xs, ys, whts)
        implicit double precision (a-h,o-z)
        dimension verts(2,nv),dr(3,2),a(2,2),ainv(2,2),verts0(2,nv)
        dimension xs(*),ys(*),whts(*),b(3,2)
c
c       Return a quadrature formula for evaluating integrals of the
c       form (1).
c
c                            Input Parameters:
c
c    verts - a (2,3) matrix which specifys the vertices of the triangle
c       over which the surface element is parameterized
c    (x0,y0) - the coordinates (in the parameterization variables)
c       of the target point x
c    dr - the (3,2) Jacobian of the parameterization at (x0,y0)
c
c                           Output Parameters:
c
c    ier - an error return code;
c       ier = 0     indicates successful execution
c       ier = 128   means that the quadrature formula could not be
c                   constructed; this usually means 
c
c    nquad - the number of nodes in the resulting quadrature formula
c    xs - this user-supplied array will contain the x coordinates of 
c       the quadrature nodes upon return
c    ys - this user-supplied array will contain the y coordinates of 
c       the quadrature nodes upone return
c    whts - this user-supplied array will contain the quadrature weights
c       upon return
c
c
        ier = 0
c
c
c       Find a linear map A which makes the Jacobian at the target node
c       conformal.  Also, compute the inverse of A.
c
        call find_conf_map(dr,a,ainv)
c
c
c       Apply the inverse of A to the vertices of the triangle T in
c       order to find the vertices of the triangle T_0.
c
        do i=1,nv
          x = verts0(1,i)-x0
          y = verts0(2,i)-y0
c
          xx = ainv(1,1)*x+ainv(1,2)*y
          yy = ainv(2,1)*x+ainv(2,2)*y
c
          verts(1,i) = xx
          verts(2,i) = yy
        enddo
c
c       Fetch a quadrature on T_0 for radially singular functions.
c
        call raddiag(irad,verts,nv,nquad,xs,ys,whts)
c
c       Apply the mapping A to the quadrature formula.
c
        det = abs(a(1,1)*a(2,2)-a(1,2)*a(2,1))
c
        do i=1,nquad
          s   = xs(i)
          t   = ys(i)
          wht = whts(i)
c       
          u = a(1,1)*s + a(1,2)*t
          v = a(2,1)*s + a(2,2)*t
          wht = wht*det
c
c
          xs(i)   = u+x0
          ys(i)   = v+y0
          whts(i) = wht
        enddo
c
        return
        end
c
c
c
c
c
        subroutine pv_self_quadrature_new(irad,verts0,nv,x0,y0,dr,
     1       nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        integer irad
        dimension verts(2,nv),dr(3,2),a(2,2),ainv(2,2),verts0(2,nv)
        dimension xs(*),ys(*),whts(1),b(3,2)
c
c       Return a quadrature formula for evaluating integrals of the
c       form (1).
c
c                            Input Parameters:
c
c    verts - a (2,3) matrix which specifys the vertices of the triangle
c       over which the surface element is parameterized
c    (x0,y0) - the coordinates (in the parameterization variables)
c       of the target point x
c    dr - the (3,2) Jacobian of the parameterization at (x0,y0)
c
c                           Output Parameters:
c
c    ier - an error return code;
c       ier = 0     indicates successful execution
c       ier = 128   means that the quadrature formula could not be
c                   constructed
c
c    nquad - the number of nodes in the resulting quadrature formula
c    xs - this user-supplied array will contain the x coordinates of
c       the quadrature nodes upon return
c    ys - this user-supplied array will contain the y coordinates of
c       the quadrature nodes upone return
c    whts - this user-supplied array will contain the quadrature weights
c       upon return
c
c
        ier = 0
c
c
c       Find a linear map A which makes the Jacobian at the target node
c       conformal.  Also, compute the inverse of A.
c
        ainv=0
        call find_conf_map(dr,a,ainv)

c
c       Apply the inverse of A to the vertices of the triangle T in
c       order to find the vertices of the triangle T_0.
c
        do i=1,nv
          x = verts0(1,i)-x0
          y = verts0(2,i)-y0
c
          xx = ainv(1,1)*x+ainv(1,2)*y
          yy = ainv(2,1)*x+ainv(2,2)*y
c
          verts(1,i) = xx
          verts(2,i) = yy
        enddo
c
c       Fetch a quadrature on T_0 for radially singular functions.
c
        call pv_raddiag(irad,verts,nv,nquad,xs,ys,whts)
c
c       Apply the mapping A to the quadrature formula.
c
        det = abs(a(1,1)*a(2,2)-a(1,2)*a(2,1))
c
c
        do i=1,nquad
          s   = xs(i)
          t   = ys(i)
          wht = whts(i)
c
          u = a(1,1)*s + a(1,2)*t
          v = a(2,1)*s + a(2,2)*t
          wht = wht*det
c
c
          xs(i)   = u+x0
          ys(i)   = v+y0
          whts(i) = wht
        enddo
c
        return
        end
c
c
c
c
c
        subroutine find_conf_map(a,b,binv)
        implicit double precision (a-h,o-z)
        dimension a(3,2),b(2,2),binv(2,2),q(3,2)
c
c       Given a (3,2) matrix A, return a (2,2) matrix B such that
c       A*B has orthonormal columns.  Also, return the inverse
c       BINV of the mapping B.
c
c
c                             --- WARNING ---
c       THIS ROUTINE IS A BIT ROUGH --- IT IS EXPECTED TO FAIL
c       IN THE EVENT THAT THE MATRIX A IS NEARLY SINGULAR.  BUT
c       WHAT ARE YOU DOING USING PARAMETERIZATION WITH JACOBIANS
c       WHICH ARE NEARLY SINGULAR ANYWAY?
c                              --------------
c
c
c                          Input Parameters:
c
c   a - the (3,2) input matrix 
c
c                         Output Parameters:
c
c   b - a (2,2) matrix such that a*b has orthonormal columns
c   binv - the (2,2) inverse of b
c
c
c       Orthonormalize the columns of the matrix a.
c
        sum1 = sqrt(a(1,1)**2 + a(2,1)** 2 + a(3,1)**2)
        q(1,1) = a(1,1) / sum1
        q(2,1) = a(2,1) / sum1
        q(3,1) = a(3,1) / sum1
c
        dip = a(1,2)*q(1,1)+a(2,2)*q(2,1)+a(3,2)*q(3,1)
c
        q(1,2) = a(1,2) - dip * q(1,1)
        q(2,2) = a(2,2) - dip * q(2,1)
        q(3,2) = a(3,2) - dip * q(3,1)
c
        sum2 = sqrt(q(1,2)**2 + q(2,2)** 2 + q(3,2)**2)
        q(1,2) = q(1,2) / sum2
        q(2,2) = q(2,2) / sum2
        q(3,2) = q(3,2) / sum2
c
c       Compute BINV = Q'*A.
c     
        binv(1,1) = q(1,1)*a(1,1) + q(2,1)*a(2,1) + q(3,1)*a(3,1)
        binv(1,2) = q(1,1)*a(1,2) + q(2,1)*a(2,2) + q(3,1)*a(3,2)
        binv(2,1) = q(1,2)*a(1,1) + q(2,2)*a(2,1) + q(3,2)*a(3,1)
        binv(2,2) = q(1,2)*a(1,2) + q(2,2)*a(2,2) + q(3,2)*a(3,2)
c
c       Compute B = (BINV)^(-1).
c
        det = binv(1,1)*binv(2,2) - binv(1,2)*binv(2,1)
        b(1,1) = binv(2,2)/det
        b(2,2) = binv(1,1)/det
        b(1,2) = -binv(1,2)/det
        b(2,1) = -binv(2,1)/det
c
        end

