
c
c
c
c
      subroutine cquadints(eps,istrat,intype,npatches,norder,
     1     ipoly,ttype,npols,
     2     srccoefs,ndtarg,ntarg,xyztarg,ifp,
     3     xyzproxy,itargptr,ntargptr,nporder,nppols,
     3     fker,ndd,dpars,ndz,zpars,ndi,ipars,nqorder,nquadmax,
     4     rfac,cintvals,ifmetric,rn1,n2)
c
c  Task:
c    This subroutine computes the integrals
c
c      \int_{[-1,1]^2} K(x_{i},\rho_{m}(y)) B_{n}(y) J_{m}(y) dy \, , 
c
c    B_{n}(y) are either tensor product legendre polynomials or
c    Chebyshev polynomials on   [-1,1]^2 of either total degree 
c    or full degree 
c
c       J_{m}(y) = det(D\rho_{m}^{T} \cdot D\rho_{m})
c
c    The maps \rho_{m} are stored in the same basis as the 
c    basis functions used in the integration  along with
c    expansions of first derivative information of xyz with 
c    respect to u,v
c
c  Method:
c    This is the guru subroutine. It uses one of three methods.
c    * strategy 1 (istrat=1): use a distance based criterion
c      to create a fixed hierarchy of quad patches for computing 
c      the integrals for each target
c    * startegy 2 (istrat=2): use an adaptive integration
c      method for computing the integrals
c    * strategy 3 (istrat=3): first create a fixed hierarchy
c      of quads based on a distance criterion and then 
c      further refine the evaluation of the integrals
c      using an adaptive method
c
c  Input arguments:
c    - eps: real *8
c        precision/tolerance requested
c    - istrat: integer
c        strategy to be used for computing the integrals
c    - npatches: integer
c        number of patches
c    - norder: integer
c        discretization order of patches
c    - ipoly: integer
c        Type of polynomial to be used
c        * ipoly = 0, Legendre polynomials
c        * ipoly = 1, Chebyshev polynomials
c    - ttype: character
c        type of number of polynomials used
c        * ttype = 'F' or 'f', full degree polynomials
c        * ttype = 'T' or 't', total degree polynomials
c    - npols: integer
c        Number of polynomials 
c        npols = norder*norder if ttype = 'F'
c        npols = (norder+1)*(norder+2)/2 if ttype = 'T'
c    - srccoefs: real *8 (9,npols,npatches)
c        coefficients of basis expansion of xyz coordinates,
c        dxyz/du, and dxyz/dv
c    - ndtarg: integer
c        leading dimension of target info array. Must be at least
c        3, and those components must correspond to the xyz
c        coordinates
c    - ntarg: integer
c        number of targets
c    - xyztarg: real *8 (ndtarg,ntarg)
c        target information
c    - ifp: integer
c        flag for whether a set of proxy targets are used for
c        determining the distance criterion
c    - xyzproxy: real *8 (3,ntarg)
c        location of proxy targets for each target, unused if
c        ifp.ne.1
c    - itargptr: integer(npatches)
c        xyztargs(:,itargptr(i):itargptr(i)+ntargptr(i)-1)
c        are the relevant set of targets for patch i
c    - ntargptr: integer(npatches)
c        ntargptr(i) is the number of relevant targets for patch i
c    - nporder: integer
c        order of basis functions to be integrated
c    - nppols: integer
c        number of basis functions to be integrated
c        * nppols = nporder*nporder, if ttype = 'F'
c        * nppols = (nporder+1)*(nporder+2)/2, if ttype = 'T'
c    - fker: function handle
c        function handle for evaluating the kernel K
c        * expected calling sequence,
c            fker(y,ndtarg,x,ndd,dpars,ndz,zpars,ndi,ipars,f)
c        The output 'f' is complex for this subroutine
c    - ndd: integer
c        number of real parameters
c    - dpars: real *8 (ndd)
c        real parameters for the fker routine
c    - ndz: integer
c        number of complex parameters
c    - zpars: complex *16 (ndz)
c        complex parameters for the fker routine
c    - ndi: integer
c        number of integer parameters
c    - ipars: integer(ndi)
c        integer parameters for the fker routine
c    - rfac: real *8
c        parameter for defining refinement criterion, quadrangles
c        refined if not separated from target by
c        rfac*r_{t} where r_{t} is radius of enclosing sphere
c        in quadrangle
c    - ifmetric: integer
c        flag for storing and priting the metrics, currently unused
c 
c  Output arguments:
c    - cintvals: complex *16 (nppols,ntarg)
c        Integral against all basis functions from patches
c        to targets
c    - rn1,rn2: real *8
c        metrics of performance of the method, currently unused
c
      implicit none

c
cc     calling sequence variables
c
      real *8 eps
      integer istrat
      integer intype
      integer npatches,norder,npols
      integer nporder,nppols
      real *8 srccoefs(9,npols,npatches)
      
      integer ntarg,ndtarg
      real *8 xyztarg(ndtarg,ntarg)
      integer ifp
      real *8 xyzproxy(3,*)
      integer itargptr(npatches)
      integer ntargptr(npatches)
      
      external fker
      integer ndd,ndi,ndz
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer ipars(ndi)
      integer nquadmax

      integer nqorder
      real *8 rfac

      integer ipoly,ier
      character *1 ttype

      integer ifmetric
      real *8 rn1

      integer n2

      complex *16 cintvals(nppols,ntarg)
      


      if(istrat.eq.1) then
          rn1 = 0
          n2 = 0
          call cquadints_dist(eps,intype,npatches,norder,ipoly,
     1        ttype,npols,srccoefs,ndtarg,ntarg,xyztarg,ifp,xyzproxy,
     2        itargptr,ntargptr,nporder,nppols,nquadmax,fker,ndd,
     3        dpars,ndz,zpars,ndi,ipars,nqorder,rfac,cintvals,ier)
      endif

      if(istrat.eq.2) then
          rn1 = 0
          n2 = 0
          call cquadints_adap(eps,intype,npatches,norder,ipoly,ttype,
     1        npols,srccoefs,ndtarg,ntarg,xyztarg,itargptr,ntargptr,
     2        nporder,nppols,nquadmax,fker,ndd,dpars,ndz,zpars,ndi,
     3        ipars,nqorder,cintvals)
      endif

      if(istrat.eq.3) then
          rn1 = 0
          n2 = 0
          call cquadints_comb(eps,intype,npatches,norder,ipoly,ttype,
     1        npols,srccoefs,ndtarg,ntarg,xyztarg,ifp,xyzproxy,itargptr,
     2        ntargptr,nporder,nppols,nquadmax,fker,ndd,dpars,ndz,zpars,
     3        ndi,ipars,nqorder,rfac,cintvals)
      endif
      return
      end




c
c
c
c



