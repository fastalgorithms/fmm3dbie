
c  This file has subroutines for computing the integrals
c
c       \int_{T} K(x_{i},\rho_{m}(y)) K_{n}(y) J_{m}(y) dy \, , (1)
c
c   where $T$ is the standard simplex,
c   K_{n}(y) are the koornwinder polynomials, and 
c   J_{m}(y) = det(D\rho_{m}^{T} \cdot D\rho_{m})
c
c   The maps \rho_{m} are stored as coefficients 
c   of the koornwinder expansion of xyz(u,v) 
c   along with expansions of first derivative 
c   information of xyz with respect to u,v.
c
c   The kernel K may be scalar or vector valued.
c
c
c  (1) This code is meant to be called with a batch of patches
c  over which computation of K_{n}(y) can be amortized.
c  This code is also going to be hard to parallelize
c  and is not going to be designed to be so.
c  Hence an outer wrapper which breaks the problem
c  down into smaller pieces and reassembles the 
c  quadrature matrix is needed. 
c
c  There are three strategies used for computing the integrals 
c
c  Strategy 1:
c    The triangle is subdivided into a collection of smaller
c    triangles, until for each combination of target and triangle
c    the target/or proxy points associated with the target are 
c    separated from the centroid of the triangle
c    by at least rfac*r_{t} where r_{t} is the radius of the 
c    smallest sphere enclosing the triangle centered at the 
c    triangle centroid
c
c    We first loop over all combinations of patches, and relevant
c    targets to find the maximal grid of triangles on 
c    which the densities is to be computed
c
c    For each source triangle, we compute all the geometry
c    information on this grid
c
c    And finally, we identify the subset relevant for a particular
c    target, source triangle combination and compute the 
c    integral
c
c  Strategy 2: 
c    Purely adaptive integration on triangles
c
c
c  Strategy 3:
c    Combination of strategy 1 and strategy 2. First compute
c    the integrals on some coarse hierarchy decided
c    based on the requested accuracy and from that
c    point on, flag all triangles to be in
c    the stack for refinement and use adaptive 
c    integration from that point on
c 
c
c  User callable routines:
c    ctriaints: routine for scalar valued complex kernels
c    ctriaints_vec: routine for vector valued complex kernels
c
      subroutine ctriaints(eps, istrat, intype, npatches, norder, npols,
     1  isd, ndsc, srccoefs, ndtarg, ntarg, xyztarg, ifp, xyzproxy, 
     2  itargptr, ntargptr, nporder, nppols, fker, ndd, dpars, ndz,
     3  zpars, ndi, ipars, nqorder, ntrimax, rfac, cintvals, ifmetric,
     4  rn1, n2)
c     
c  Routine for computing scalar valued complex integrals of the
c  form (1) 
c  
c        
c  Input arguments:
c    - eps: real *8
c        requested tolerance. Currently unsued, nqorder and rfac,
c        should be set based on eps.
c    - istart: integer
c        strategy to be used for computing integrals
c        * istrat = 1, distance based criterion
c        * istrat = 2, adaptive integration
c        * istrat = 3, combination of istrat = 1, and istrat = 2
c    - intype: integer
c        node type
c        * intype = 1, rokhlin vioreanu nodes
c        * intype = 2, xiao gimbutas nodes
c    - npatches: integer
c        number of patches
c    - norder: integer
c        order of discretization on patches
c    - npols: integer
c        number of points/basis functions on the 
c        patch = (norder+1)*(norder+2)/2
c    - isd: integer
c        flag for computing and passing second derivative
c        information in srcvals.
c        if isd = 0, leading dim of srcvals = 12, with 
c        xyz, dxyz/du, dxyz/dv, and normals passed
c        else, leading dim of srcvals = 24, with
c        xyz, dxyz/du, dxyz/dv, normals, d2xyz/du2, d2xyz/duv,
c        d2xyz/dv2, det(g), kap1, kap2 where kap1, kap2
c        are principal curvatures, det(g) is the determinant
c        of the metric tensor
c    - ndsc: integer
c        leading dimension of srccoefs array, should be 9
c        if isd = 0, and 18 otherwise
c    - srccoefs: real *8 (ndsc, npols, npatches)
c        Koornwinder expansion coefficients of xyz, d/du (xyz), 
c        d/dv (xyz), (d2(xyz)/du2, d2(xyz)/duv, d2(xyz)/dv2)
c    - ndtarg: integer
c        leading dimension of target information
c    - ntarg: integer
c        number of targets
c    - xyztarg: real *8 (ndtarg, ntarg)
c        target information
c    - ifp: integer
c        flag for using proxy points for refining based 
c        on distance criterion
c    - xyzproxy: real *8 (3,*)
c        Proxy target points for refining triangles based
c        distance criterion, distance to proxy
c        point will be used instead of distance to target
c        should be of size (3,ntarg), if ifp=1
c    - itargptr: integer(npatches)
c        Pointer array for determining which targets are relevant
c        for patch i. 
c        xyztargs(:,itargptr(i):itargptr(i)+ntargptr(i)-1)
c        are the relevant set of targets for patch i
c    - ntargptr: integer(npatches)
c        number of targets relevant for patch i
c    - nporder: integer
c        order of Koornwinder polynomials to be integrated
c    - nppols: integer
c        number of polynomials corresponding to 
c        nporder = (nporder+1)*(nporder+2)/2
c    - fker: function handle
c        function handle for evaluating the kernel k
c        * expected calling sequence
c            fker(x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
c               
c        In this routine the output is expected to be a complex
c        scalar
c    - ndd: integer
c        number of real *8/double precision parameters
c    - dpars: real *8 (ndd)
c         real *8/ double precision paramters
c    - ndz: integer
c        number of complex *16 parameters
c    - zpars: complex *16(ndz)
c        complex *16 parameters
c    - ndi: integer
c        number of integer parameters
c    - ipars: integer(ndi)
c        integer parameters
c    - nqorder: integer
c        order of quadrature nodes to be used on each triangle
c    - ntrimax: integer
c        maximum number of triangles allowed in heirarchy
c        of triangles. Routine will return without
c        computing anything and an error code, if ntrimax
c        is too small. Recommended value 3000.
c    - rfac: real *8
c        scaling factor for determining when a triangle 
c        needs refinement
c    - ifmetric: integer
c        flag for computing metrics for adaptive integration
c        (Currently not in use)
c 
c  Output arguments:
c    - cintvals: complex *8 (nppols, ntarg)
c        the computed integrals
c    - rn1: real *8
c        metric for number of function evaluations (unused)
c    - n2: integer
c        metric for number of function evaluations (unused)
c
      implicit none

c
cc     calling sequence variables
c
      real *8 eps
      integer *8 istrat
      integer *8 intype
      integer *8 npatches,norder,npols
      integer *8 nporder,nppols
      real *8 srccoefs(9,npols,npatches)
      
      integer *8 ntarg,ndtarg
      real *8 xyztarg(ndtarg,ntarg)
      integer *8 ifp
      real *8 xyzproxy(3,*)
      integer *8 itargptr(npatches)
      integer *8 ntargptr(npatches)
      
      external fker
      integer *8 ndd,ndi,ndz
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer *8 ipars(ndi)
      integer *8 ntrimax

      integer *8 nqorder
      integer *8 isd, ndsc
      real *8 rfac

      integer *8 ifmetric
      real *8 rn1

      real *8, allocatable :: rat1(:,:), rat2(:,:,:), rsc1(:,:)
      real *8, allocatable :: rat1p(:,:), rat2p(:,:,:), rsc1p(:,:)

      integer *8 n2
      integer *8 ier

      complex *16 cintvals(nppols,ntarg)
      
      allocate(rat1(2,0:norder),rat2(3,0:norder,0:norder))
      allocate(rsc1(0:norder,0:norder))
       
      call koornf_init(norder,rat1,rat2,rsc1)


      allocate(rat1p(2,0:nporder),rat2p(3,0:nporder,0:nporder))
      allocate(rsc1p(0:nporder,0:nporder))
       
      call koornf_init(nporder,rat1p,rat2p,rsc1p)



      if(istrat.eq.1) then
          rn1 = 0
          n2 = 0
          ier = 0
          
          call ctriaints_dist(eps, intype, npatches, norder, npols,
     1      isd, ndsc, srccoefs, ndtarg, ntarg, xyztarg, ifp, 
     2      xyzproxy, itargptr, ntargptr, nporder, nppols, ntrimax,
     3      rat1, rat2, rsc1, rat1p, rat2p, rsc1p, fker, ndd, dpars,
     4      ndz, zpars, ndi, ipars, nqorder, rfac, cintvals, ier)
          if(ier.ne.0) then
              print *, "Too few tirangles in ctriaints_dist"
              print *, "Try running with larger value of ntrimax"
              return
          endif
      endif

      if(istrat.eq.2) then
          rn1 = 0
          n2 = 0
          call ctriaints_adap(eps, intype, npatches, norder, npols,
     1      isd, ndsc, srccoefs, ndtarg, ntarg, xyztarg, itargptr,
     2      ntargptr, nporder, nppols, ntrimax, rat1, rat2, rsc1,
     3      rat1p, rat2p, rsc1p, fker, ndd, dpars, ndz, zpars, ndi,
     4      ipars, nqorder, cintvals)
      endif

      if(istrat.eq.3) then
          rn1 = 0
          n2 = 0
          call ctriaints_comb(eps, intype, npatches, norder, npols,
     1     isd, ndsc, srccoefs, ndtarg, ntarg, xyztarg, ifp, xyzproxy,
     2     itargptr, ntargptr, nporder, nppols, ntrimax, rat1, rat2,
     3     rsc1, rat1p, rat2p, rsc1p, fker, ndd, dpars, ndz, zpars, 
     4     ndi, ipars, nqorder, rfac, cintvals)
      endif
      return
      end
c
c
c
c
      subroutine ctriaints_vec(eps, istrat, intype, npatches, norder, 
     1  npols, isd, ndsc, srccoefs, ndtarg, ntarg, xyztarg, ifp, 
     2  xyzproxy, itargptr, ntargptr, nporder, nppols, fker, nd, ndd, 
     3  dpars, ndz, zpars, ndi, ipars, nqorder, ntrimax, rfac, 
     4  cintvals, ifmetric, rn1, n2)
c
c     
c  Routine for computing vector valued complex integrals of the
c  form (1) 
c  
c  Input arguments:
c    - eps: real *8
c        requested tolerance. Currently unsued, nqorder and rfac,
c        should be set based on eps.
c    - istart: integer
c        strategy to be used for computing integrals
c        * istrat = 1, distance based criterion
c        * istrat = 2, adaptive integration
c        * istrat = 3, combination of istrat = 1, and istrat = 2
c    - intype: integer
c        node type
c        * intype = 1, rokhlin vioreanu nodes
c        * intype = 2, xiao gimbutas nodes
c    - npatches: integer
c        number of patches
c    - norder: integer
c        order of discretization on patches
c    - npols: integer
c        number of points/basis functions on the 
c        patch = (norder+1)*(norder+2)/2
c    - isd: integer
c        flag for computing and passing second derivative
c        information in srcvals.
c        if isd = 0, leading dim of srcvals = 12, with 
c        xyz, dxyz/du, dxyz/dv, and normals passed
c        else, leading dim of srcvals = 24, with
c        xyz, dxyz/du, dxyz/dv, normals, d2xyz/du2, d2xyz/duv,
c        d2xyz/dv2, det(g), kap1, kap2 where kap1, kap2
c        are principal curvatures, det(g) is the determinant
c        of the metric tensor
c    - ndsc: integer
c        leading dimension of srccoefs array, should be 9
c        if isd = 0, and 18 otherwise
c    - srccoefs: real *8 (ndsc, npols, npatches)
c        Koornwinder expansion coefficients of xyz, d/du (xyz), 
c        d/dv (xyz), (d2(xyz)/du2, d2(xyz)/duv, d2(xyz)/dv2)
c    - ndtarg: integer
c        leading dimension of target information
c    - ntarg: integer
c        number of targets
c    - xyztarg: real *8 (ndtarg, ntarg)
c        target information
c    - ifp: integer
c        flag for using proxy points for refining based 
c        on distance criterion
c    - xyzproxy: real *8 (3,*)
c        Proxy target points for refining triangles based
c        distance criterion, distance to proxy
c        point will be used instead of distance to target
c        should be of size (3,ntarg), if ifp=1
c    - itargptr: integer(npatches)
c        Pointer array for determining which targets are relevant
c        for patch i. 
c        xyztargs(:,itargptr(i):itargptr(i)+ntargptr(i)-1)
c        are the relevant set of targets for patch i
c    - ntargptr: integer(npatches)
c        number of targets relevant for patch i
c    - nporder: integer
c        order of Koornwinder polynomials to be integrated
c    - nppols: integer
c        number of polynomials corresponding to 
c        nporder = (nporder+1)*(nporder+2)/2
c    - fker: function handle
c        function handle for evaluating the kernel k
c        * expected calling sequence
c            fker(nd,x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
c               
c        In this routine the output is expected to be a complex
c        vector
c    - nd: integer
c        number of components in the kernel
c    - ndd: integer
c        number of real *8/double precision parameters
c    - dpars: real *8 (ndd)
c         real *8/ double precision paramters
c    - ndz: integer
c        number of complex *16 parameters
c    - zpars: complex *16(ndz)
c        complex *16 parameters
c    - ndi: integer
c        number of integer parameters
c    - ipars: integer(ndi)
c        integer parameters
c    - nqorder: integer
c        order of quadrature nodes to be used on each triangle
c    - ntrimax: integer
c        maximum number of triangles allowed in heirarchy
c        of triangles. Routine will return without
c        computing anything and an error code, if ntrimax
c        is too small. Recommended value 3000.
c    - rfac: real *8
c        scaling factor for determining when a triangle 
c        needs refinement
c    - ifmetric: integer
c        flag for computing metrics for adaptive integration
c        (Currently not in use)
c 
c  Output arguments:
c    - cintvals: complex *16 (nd, nppols, ntarg)
c        the computed integrals
c    - rn1: real *8
c        metric for number of function evaluations (unused)
c    - n2: integer
c        metric for number of function evaluations (unused)
c
c
      implicit none

c
cc     calling sequence variables
c
      real *8 eps
      integer *8 istrat
      integer *8 intype
      integer *8 npatches,norder,npols
      integer *8 nporder,nppols
      integer *8 nd
      real *8 srccoefs(9,npols,npatches)
      
      integer *8 ntarg,ndtarg
      real *8 xyztarg(ndtarg,ntarg)
      integer *8 ifp
      real *8 xyzproxy(3,*)
      integer *8 itargptr(npatches)
      integer *8 ntargptr(npatches)
      
      external fker
      integer *8 ndd,ndi,ndz
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer *8 ipars(ndi)
      integer *8 ntrimax

      integer *8 nqorder
      integer *8 isd, ndsc
      real *8 rfac

      integer *8 ifmetric
      real *8 rn1

      real *8, allocatable :: rat1(:,:), rat2(:,:,:), rsc1(:,:)
      real *8, allocatable :: rat1p(:,:), rat2p(:,:,:), rsc1p(:,:)



      integer *8 n2

      complex *16 cintvals(nd,nppols,ntarg)

      allocate(rat1(2,0:norder),rat2(3,0:norder,0:norder))
      allocate(rsc1(0:norder,0:norder))
       
      call koornf_init(norder,rat1,rat2,rsc1) 

      allocate(rat1p(2,0:nporder),rat2p(3,0:nporder,0:nporder))
      allocate(rsc1p(0:nporder,0:nporder))
       
      call koornf_init(nporder,rat1p,rat2p,rsc1p) 


      if(istrat.eq.1) then
          rn1 = 0
          n2 = 0
          call ctriaints_dist_vec(eps, intype, npatches, norder, npols,
     1      isd, ndsc, srccoefs, ndtarg, ntarg, xyztarg, ifp, 
     2      xyzproxy, itargptr, ntargptr, nporder, nppols, ntrimax,
     3      rat1, rat2, rsc1, rat1p, rat2p, rsc1p, fker, nd, ndd, 
     4      dpars, ndz, zpars, ndi, ipars, nqorder, rfac, cintvals)
      endif

      if(istrat.eq.2) then
          rn1 = 0
          n2 = 0
          call ctriaints_adap_vec(eps, intype, npatches, norder, npols,
     1      isd, ndsc, srccoefs, ndtarg, ntarg, xyztarg, itargptr,
     2      ntargptr, nporder, nppols, ntrimax, rat1, rat2, rsc1,
     3      rat1p, rat2p, rsc1p, fker, nd, ndd, dpars, ndz, zpars, ndi,
     4      ipars, nqorder, cintvals)
      endif

      if(istrat.eq.3) then
          rn1 = 0
          n2 = 0
          call ctriaints_comb_vec(eps, intype, npatches, norder, npols,
     1     isd, ndsc, srccoefs, ndtarg, ntarg, xyztarg, ifp, xyzproxy,
     2     itargptr, ntargptr, nporder, nppols, ntrimax, rat1, rat2,
     3     rsc1, rat1p, rat2p, rsc1p, fker, nd, ndd, dpars, ndz, zpars, 
     4     ndi, ipars, nqorder, rfac, cintvals)
      endif


      return
      end



