
c
c
c
c
      subroutine ctriaints(eps,istrat,intype,npatches,norder,npols,
     2     srccoefs,ndtarg,ntarg,xyztarg,ifp,
     3     xyzproxy,itargptr,ntargptr,nporder,nppols,
     3     fker,ndd,dpars,ndz,zpars,ndi,ipars,nqorder,ntrimax,
     4     rfac,cintvals,ifmetric,rn1,n2)
c
c
c       this subroutine computes the integrals
c
c       \int_{T} K(x_{i},\rho_{m}(y)) K_{n}(y) J_{m}(y) dy \, ,
c
c        where $T$ is the standard simplex,
c
c        K_{n}(y) are the koornwinder polynomials
c
c        J_{m}(y) = det(D\rho_{m}^{T} \cdot D\rho_{m})
c
c        The maps \rho_{m} are stored as coefficients 
c        of the koornwinder expansion of xyz(u,v) 
c        along with expansions of first and
c        second derivative information of xyz with respect
c        to u,v. 
c
c
c        (1) This code is meant to be called with a batch of patches
c        over which computation of K_{n}(y) can be amortized.
c        This code is also going to be hard to parallelize
c        and is not going to be designed to be so.
c        Hence an outer wrapper which breaks the problem
c        down into smaller pieces and reassembles the 
c        quadrature matrix is needed. 
c
c        There are two strategies used for adaptive integration
c
c         strategy 1:
c           The triangle is subdivided into a collection of smaller
c        triangles, until for each combination of target and triangle
c        the target/or proxy points associated with the target are 
c        separated from the centroid of the triangle
c        by at least rfac*r_{t} where r_{t} is the radius of the 
c        smallest sphere enclosing the triangle centered at the 
c        triangle centroid
c
c           we first loop over all combinations of patches, and relevant
c        targets to find the maximal grid of triangles on 
c        which the densities is to be computed
c
c          For each source triangle, we compute all the geometry
c        information on this grid
c
c          And finally, we identify the subset relevant for a particular
c        target, source triangle combination and compute the 
c        integral
c
c          strategy 2: 
c             Purely adaptive integration on triangles.
c             if number of triangles exceeds ntrimax, then
c             adaptive integration is restarted after allocating
c             more triangles
c           
c
c
c          strategy 3:
c             combination of strategy 1 and strategy 2. First compute
c             the integrals on some coarse hierarchy decided
c             based on the requested accuracy and from that
c             point on, flag all triangles to be in
c             the stack for refinement and use adaptive 
c             integration from that point on.
c 
c             if number of triangles exceeds ntrimax, then
c             computation for current target is currently
c             halted there and the subroutine returns
c             outside the code
c 
c        
c        input arguments:
c        eps:     requested precision (Currently only
c                  used for adaptive strategy)
c
c        istart:  Adaptive integration strategy
c                 istrat = 1, refinement based on separation
c                  criterion
c
c                 istrat = 2, adpative integration on triangles
c
c        intype:   quadrature node type
c                   intype = 1, rokhlin vioreanu nodes
c                   intype = 2, xiao gimbutas nodes
c   
c
c        npatches: number of patches
c        norder: order of discretization nodes on the patches
c        npols = (norder+1)*(norder+2)/2: number of discretization 
c                   nodes on each patch
c        srccoefs(9,npols,npatches): coefficients
c                 of koornwinder expansion of xyz coordinates,
c                 and that of dxyz/du, dxyz/dv
c
c 
c        ndtarg - dimension of each target point (if istart =1/3,
c          ndtarg must at least be 3, and correspond to location
c          of target points)
c        ntarg - total number of target points
c        xyztarg(ndtarg,ntarg) - 
c                       target vectors  
c                       (currently may have repeats, 
c                        being slightly wasteful in
c                        memory, but see comment (1))
c
c        ifp - flag for proxy target locations for startegy 1
c
c        xyzproxy(3,ntarg) - proxy locations for measuring distances
c
c        itargptr(npatches) - xyztargs(:,itargptr(i):itargptr(i)+ntargptr(i)-1)
c                       are the relevant set of targets for patch i
c        ntargptr(npatches) - number of relevant targets for patch i
c
c        nporder - order of koornwinder polynomials to be integrated
c        nppols - number of koornwinder polynomials to be integrated 
c                  (nppols = (nporder+1)*(nporder+2)/2)
c
c
c        fker - function handle for evaluating the kernel k
c 
c               expected calling sequence
c               fker(x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
c               
c               the output is assumed to be complex for the time
c               being
c
c         ndd - number of real parameters
c         dpars(ndd) - real parameters for the fker routine
c         ndz - number of complex parameters
c         zpars(ndz) - complex parameters for the fker routine
c         ndi - number of integer(8) parameters
c         ipars(ndi) - integer(8) parameters for the fker routine
c         nqorder - order of quadrature nodes on each subtriangle
c                   to be used
c         ntrimax - max number of triangles to be used for 
c                    any of the integrals 
c         rfac - distance criterion for deciding whether a triangle
c                is in the far-field or not (See comment (2))
c
c         note: the ifmetric feature is currently not in use
c
c         output:
c
c         cintvals(nppols,ntarg) - integrals at all targets
c                                  for all koornwinder
c                                  polynomials
c
      implicit none

c
cc     calling sequence variables
c
      real *8 eps
      integer(8) istrat
      integer(8) intype
      integer(8) npatches,norder,npols
      integer(8) nporder,nppols
      real *8 srccoefs(9,npols,npatches)
      
      integer(8) ntarg,ndtarg
      real *8 xyztarg(ndtarg,ntarg)
      integer(8) ifp
      real *8 xyzproxy(3,*)
      integer(8) itargptr(npatches)
      integer(8) ntargptr(npatches)
      
      external fker
      integer(8) ndd,ndi,ndz
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer(8) ipars(ndi)
      integer(8) ntrimax

      integer(8) nqorder
      real *8 rfac

      integer(8) ifmetric
      real *8 rn1

      real *8, allocatable :: rat1(:,:), rat2(:,:,:), rsc1(:,:)
      real *8, allocatable :: rat1p(:,:), rat2p(:,:,:), rsc1p(:,:)



      integer(8) n2

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
          
          call ctriaints_dist(eps,intype,npatches,norder,npols,
     1        srccoefs,ndtarg,ntarg,xyztarg,ifp,xyzproxy,itargptr,
     2        ntargptr,nporder,nppols,ntrimax,rat1,rat2,rsc1,
     3        rat1p,rat2p,rsc1p,fker,ndd,
     3        dpars,ndz,zpars,ndi,ipars,nqorder,rfac,cintvals)
      endif

      if(istrat.eq.2) then
          rn1 = 0
          n2 = 0
          call ctriaints_adap(eps,intype,npatches,norder,npols,
     1        srccoefs,ndtarg,ntarg,xyztarg,itargptr,
     2        ntargptr,nporder,nppols,ntrimax,rat1,rat2,rsc1,
     3        rat1p,rat2p,rsc1p,fker,ndd,dpars,ndz,zpars,ndi,ipars,
     4        nqorder,cintvals)
      endif

      if(istrat.eq.3) then
          rn1 = 0
          n2 = 0
          call ctriaints_comb(eps,intype,npatches,norder,npols,
     1        srccoefs,ndtarg,ntarg,xyztarg,ifp,xyzproxy,itargptr,
     2        ntargptr,nporder,nppols,ntrimax,rat1,rat2,rsc1,rat1p,
     3        rat2p,rsc1p,fker,ndd,dpars,ndz,zpars,ndi,ipars,nqorder,
     4        rfac,cintvals)
      endif
      return
      end




c
c
c
c
      subroutine ctriaints_vec(eps,istrat,intype,npatches,norder,npols,
     2     srccoefs,ndtarg,ntarg,xyztarg,ifp,
     3     xyzproxy,itargptr,ntargptr,nporder,nppols,fker,nd,ndd,dpars,
     3     ndz,zpars,ndi,ipars,nqorder,ntrimax,rfac,cintvals,ifmetric,
     4     rn1,n2)
c
c
c       this subroutine computes the integrals
c
c       \int_{T} K_{\ell}(x_{i},\rho_{m}(y)) P_{n}(y) J_{m}(y) dy \, ,
c
c        where $T$ is the standard simplex,
c
c        P_{n}(y) are the koornwinder polynomials
c
c        J_{m}(y) = det(D\rho_{m}^{T} \cdot D\rho_{m})
c
c        The maps \rho_{m} are stored as coefficients 
c        of the koornwinder expansion of xyz(u,v) 
c        along with expansions of first and
c        second derivative information of xyz with respect
c        to u,v. 
c
c
c        (1) This code is meant to be called with a batch of patches
c        over which computation of P_{n}(y) can be amortized.
c        This code is also going to be hard to parallelize
c        and is not going to be designed to be so.
c        Hence an outer wrapper which breaks the problem
c        down into smaller pieces and reassembles the 
c        quadrature matrix is needed. 
c
c        There are two strategies used for adaptive integration
c
c          Strategy 1:
c
c          The triangle is subdivided into a collection of smaller
c        triangles, until for each combination of target and triangle
c        the target/or proxy points associated with the target are 
c        separated from the centroid of the triangle
c        by at least rfac*r_{t} where r_{t} is the radius of the 
c        smallest sphere enclosing the triangle centered at the 
c        triangle centroid
c
c           we first loop over all combinations of patches, and relevant
c        targets to find the maximal grid of triangles on 
c        which the densities is to be computed
c
c          For each source triangle, we compute all the geometry
c        information on this grid
c
c          And finally, we identify the subset relevant for a particular
c        target, source triangle combination and compute the 
c        integral
c
c          strategy 2: 
c             Purely adaptive integration on triangles
c
c
c          strategy 3:
c             combination of strategy 1 and strategy 2. First compute
c             the integrals on some coarse hierarchy decided
c             based on the requested accuracy and from that
c             point on, flag all triangles to be in
c             the stack for refinement and use adaptive 
c             integration from that point on
c 
c        
c        input arguments:
c        eps:     requested precision (Currently only
c                  used for adaptive strategy)
c
c        istart:  Adaptive integration strategy
c                 istrat = 1, refinement based on separation
c                  criterion
c
c                 istrat = 2, adpative integration on triangles
c
c        intype:   quadrature node type
c                   intype = 1, rokhlin vioreanu nodes
c                   intype = 2, xiao gimbutas nodes
c   
c
c        npatches: number of patches
c        norder: order of discretization nodes on the patches
c        npols = (norder+1)*(norder+2)/2: number of discretization 
c                   nodes on each patch
c        srccoefs(9,npols,npatches): coefficients
c                 of koornwinder expansion of xyz coordinates,
c                 and that of dxyz/du, dxyz/dv
c
c        ndtarg - dimension of each target point (if istart =1/3,
c          ndtarg must at least be 3, and correspond to location
c          of target points)
c        ntarg - total number of target points
c        xyztarg(ndtarg,ntarg) - target vectors 
c                       (currently may have repeats, 
c                        being slightly wasteful in
c                        memory, but see comment (1))
c
c        ifp - flag for proxy target locations for startegy 1
c
c        xyzproxy(3,ntarg) - proxy locations for measuring distances
c
c        itargptr(npatches) - xyztargs(:,itargptr(i):itargptr(i)+ntargptr(i)-1)
c                       are the relevant set of targets for patch i
c        ntargptr(npatches) - number of relevant targets for patch i
c        nporder - order of koornwinder polynomials to be integrated
c        nppols - number of koornwinder polynomials to be integrated 
c                  (nppols = (nporder+1)*(nporder+2)/2)
c
c        fker - function handle for evaluating the kernel k
c 
c               expected calling sequence
c               fker(nd,x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
c               
c               the output is assumed to be nd complex numbers 
c
c         nd - number of outputdimensions
c         ndd - number of real parameters
c         dpars(ndd) - real parameters for the fker routine
c         ndz - number of complex parameters
c         zpars(ndz) - complex parameters for the fker routine
c         ndi - number of integer(8) parameters
c         ipars(ndi) - integer(8) parameters for the fker routine
c         nqorder - order of quadrature nodes on each subtriangle
c                   to be used
c         ntrimax - max number of triangles to be used for 
c                    any of the integrals 
c         rfac - distance criterion for deciding whether a triangle
c                is in the far-field or not (See comment (2))
c
c         output:
c
c         cintvals(nd,nppols,ntarg) - integrals at all targets
c                                  for all koornwinder
c                                  polynomials
c
      implicit none

c
cc     calling sequence variables
c
      real *8 eps
      integer(8) istrat
      integer(8) intype
      integer(8) npatches,norder,npols
      integer(8) nporder,nppols
      integer(8) nd
      real *8 srccoefs(9,npols,npatches)
      
      integer(8) ntarg,ndtarg
      real *8 xyztarg(ndtarg,ntarg)
      integer(8) ifp
      real *8 xyzproxy(3,*)
      integer(8) itargptr(npatches)
      integer(8) ntargptr(npatches)
      
      external fker
      integer(8) ndd,ndi,ndz
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer(8) ipars(ndi)
      integer(8) ntrimax

      integer(8) nqorder
      real *8 rfac

      integer(8) ifmetric
      real *8 rn1

      real *8, allocatable :: rat1(:,:), rat2(:,:,:), rsc1(:,:)
      real *8, allocatable :: rat1p(:,:), rat2p(:,:,:), rsc1p(:,:)



      integer(8) n2

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
          call ctriaints_dist_vec(eps,intype,npatches,norder,npols,
     1        srccoefs,ndtarg,ntarg,xyztarg,ifp,xyzproxy,itargptr,
     2        ntargptr,nporder,nppols,ntrimax,rat1,rat2,rsc1,rat1p,
     3        rat2p,rsc1p,fker,nd,ndd,dpars,ndz,zpars,ndi,ipars,
     3        nqorder,rfac,cintvals)
      endif

      if(istrat.eq.2) then
          rn1 = 0
          n2 = 0
          call ctriaints_adap_vec(eps,intype,npatches,norder,npols,
     1        srccoefs,ndtarg,ntarg,xyztarg,itargptr,
     2        ntargptr,nporder,nppols,ntrimax,rat1,rat2,rsc1,
     3        rat1p,rat2p,rsc1p,fker,nd,ndd,dpars,ndz,zpars,ndi,ipars,
     4        nqorder,cintvals)
      endif

      if(istrat.eq.3) then
          rn1 = 0
          n2 = 0
          call ctriaints_comb_vec(eps,intype,npatches,norder,npols,
     1        srccoefs,ndtarg,ntarg,xyztarg,ifp,xyzproxy,itargptr,
     2        ntargptr,nporder,nppols,ntrimax,rat1,rat2,rsc1,rat1p,
     3        rat2p,rsc1p,fker,nd,ndd,dpars,ndz,zpars,ndi,ipars,nqorder,
     4        rfac,cintvals)
      endif


      return
      end



