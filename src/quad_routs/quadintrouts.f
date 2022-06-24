c
c  Notes: The subroutines in this file assume that all patches
c  are discretized with the same order and same type of
c  polynomials
c
c  Typical usage of this subroutine is with one patch only
c
c  Routines in this file
c   - cquadints_adap: 
c       completely adaptive integration for the quadrangle.
c   - cquadints_wnodes: 
c       compute integrals with prescribed nodes and weights.
c   - cquadints_dist:
c       refine until targets/proxy separated
c       by quadrangle by rfac*r_{t} where r_{t} is the radius of 
c       enclosing sphere for quadrangle.
c   - cquadints_comb: 
c      combination of adaptive integration and distance based refinement.
c      Use distance based refinement to compute the integrals on some 
c      coarse hierarchy of meshes and then use adaptive integration 
c      from that point on.
c           
c  Vectorized routines
c ----------------------------
c   By vectorized version of the routines, we mean here
c   that the kernels to be integrated are vector valued
c   as opposed to scalar valued
c
c         
c  - cquadints_adap_vec: 
c       completely adaptive integration for the quadrangle.
c       vectorized verion of cquadints_adap
c   - cquadints_wnodes_vec: 
c       compute integrals with prescribed nodes and weights.
c       vectorized version of cquadints_wnodes
c  - cquadints_dist_vec: 
c       refine until targets/proxy separated
c       by quadrangle by rfac*r_{t} where r_{t} is the radius of 
c       enclosing sphere for quadrangle.
c       vectorized verison of cquadints_dist
c  - cquadints_comb_vec:
c      combination of adaptive integration and distance based refinement.
c      Use distance based refinement to compute the integrals on some 
c      coarse hierarchy of meshes and then use adaptive integration 
c      from that point on.
c      vectorized version of cquadints_comb
c               
c  We integrate against chebyshev/legendre of full order or 
c  total order polynomials on the standard patch [-1,1]^2
c                       
c
c

      subroutine cquadints_wnodes(npatches,norder,ipoly,ttype,npols,
     1    srccoefs,ndtarg,ntarg,xyztarg,itargptr,ntargptr,
     2    nporder,nppols,fker,ndd,dpars,ndz,zpars,ndi,ipars,
     3    nqpts,qnodes,wts,cintvals)
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
c    The subroutine uses prescribed nodes and weights to compute 
c    the integrals.
c
c
c
c  Input arguments:
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
c    - nqpts: integer
c        number of quadrature points used in the
c        integration routine
c    - qnodes: real *8 (2,nqpts)
c        location of quadrature nodes on (-1,1)^2
c    - wts: real *8 (nqpts)
c        the corresponding quadrature weights
c 
c  Output arguments:
c    - cintvals: complex *16 (nppols,ntarg)
c        Integral against all basis functions from patches
c        to targets
c
c------------------------------------- 
c
      implicit none
      integer, intent(in) :: npatches,norder,ipoly,npols,ndtarg,ntarg
      character *1, intent(in) :: ttype
      real *8, intent(in) :: srcoefs(9,npols,nptaches)
      real *8, intent(in) :: xyztarg(ndtarg,ntarg)
      integer, intent(in) :: itargptr(npatches),ntargptr(npatches)
      integer, intent(in) :: nporder,nppols,ndd,ndz,ndi,ipars(ndi),nqpts
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      real *8, intent(in) :: qnodes(2,nqpts),wts(nqpts)
      complex *16, intent(out) :: cintvals(nppols,ntarg)
      external fker

!
!  temporary variables
!
      complex *16, allocatable :: xkernvals(:,:)
      complex *16, allocatable :: sigvals(:,:)
      real *8, allocatable : srcvals(:,:),qwts(:),rsigvals(:,:)
      real *8, allocatable :: pols_tmp(:)

      integer i,j
      character *1 transa,transb

      allocate(pols_tmp(npols))
      allocate(sigvals(nppols,nqpts),xkernvals(nqpts,ntarg))
      allocate(rsigvals(npols,nqpts))

c
c
c       generate the density values
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,pols_tmp)     
      do i=1,nqpts
        call polytens_pols_2d(ipoly,qnodes(1,i),norder-1,
     1     ttype,rsigvals(1,i))
        call polytens_pols_2d(ipoly,qnodes(1,i),nporder-1,
     1     ttype,pols_tmp) 
        do j=1,nppols
          sigvals(j,i) = pols_tmp(j)
        enddo
      enddo
C$OMP END PARALLEL DO      


      ntargmax = maxval(ntargptr(:))
      allocate(xkernvals(nqpts,ntargmax))

      transa = 'N'
      transb = 'N'
      alpha = 1.0d0
      beta = 0.0d0
      lda = 9
      ldb = npols
      ldc = 12

      da = 1



      do ipatch=1,npatches
        call dgemm_guru(transa,transb,9,nqpts,npols,alpha,
     1    srccoefs(1,1,ipatch),lda,rsigvals,ldb,beta,srcvals,ldc)
        
        call get_norms_qwts_quad(nqpts,wts,srcvals,da,qwts)
        ntarg0 = ntargptr(ipatch)
        do itarg = itargptr(ipatch),itargptr(ipatch)+ntarg0-1
          ii = itarg - itargptr(ipatch)+1
          do j=1,nqpts
            call fker(srcvals(1,j),ndtarg,xyztarg(1,itarg),ndd,dpars,
     1        ndz,zpars,ndi,ipars,fval)
            xkernvals(j,ii) = fval*qwts(j)
          enddo
        enddo
        alpha_c = 1.0d0
        beta_c = 0.0d0
        call zgemm_guru(transa,transb,nppols,ntarg0,nqpts,alpha_c,
     1    sigvals,nppols,xkernvals,nqpts,beta_c,
     2    cintvals(1,itargptr(ipatch)),nppols)
      enddo
        

      deallocate(xkernvals)
      return
      end
c
c
c
c
c
c
      subroutine cquadints_dist(eps,intype,
     1     npatches,norder,ipoly,ttype,npols,srccoefs,ndtarg,
     2     ntarg,xyztarg,ifp,xyzproxy,
     3     itargptr,ntargptr,nporder,nppols,nquadmax,
     4     fker,ndd,dpars,ndz,zpars,ndi,ipars,
     5     nqorder,rfac,cintvals,ier)
c-------------
c  Task:
c    This subroutine computes the integrals
c
c      \int_{[-1,1]^2} K(x_{i},\rho_{m}(y)) B_{n}(y) J_{m}(y) dy \, , 
c
c    B_{n}(y) are either tensor product legendre polynomials or
c    Chebyshev polynomials on [-1,1]^2 of either total degree 
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
c    Refine until targets/proxy separated by quadrangle by 
c    rfac*r_{t} where r_{t} is the radius of enclosing sphere 
c    for quadrangle.
c
c  Notes:
c    There are two parameters which determine the quadrature
c    nodes used on each quadrangle in the heirarchy,
c    intype, and nqorder
c    intype switches between chebyshev, legendre, and xiao gimbutas
c    nodes of full or total degree, and nqorder determines
c    the order of nodes used
c
c  Input arguments:
c    - eps: real *8
c        precision requested (unused in this routine)
c    - intype: integer
c        type of nodes to be used as quadrature nodes for each
c        quadrangle in the heirarchy:
c          intype = 1, use either tensor product chebyshev/legendre
c            nodes of total/full degree based on the input parameters
c            ipoly, and ttype
c          intype = 2, use xiao gimbutas nodes on the square
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
c    - nquadmax: integer
c        max number of quadrangles supported in the heirarchy
c        if total number of quadrangles required exceeds
c        nquadmax, subroutine exits with error code
c        ier = 4
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
c    - nqorder: integer
c        order of quadrature nodes to be used
c    - rfac: real *8
c        parameter for defining refinement criterion, quadrangles
c        refined if not separated from target by
c        rfac*r_{t} where r_{t} is radius of enclosing sphere
c        in quadrangle
c  Output arguments:
c    - cintvals: complex *16 (nppols,ntarg)
c        Integral against all basis functions from patches
c        to targets
c    - ier: integer
c        error code, ier = 0, implies successful execution
c        ier = 4, not enough quadrangles, try rerunning the
c          routine with more quads
c
c-------------------
      implicit none

c
cc     calling sequence variables
c
      real *8, intent(in) :: eps
      integer, intent(in) :: intype,ifp
      integer, intent(in) :: npatches,norder,npols,ipoly
      character *1, intent(in) :: ttype
      integer, intent(in) :: nporder,nppols
      real *8, intent(in) :: srccoefs(9,npols,npatches)
      
      integer, intent(in) :: ntarg,ndtarg
      real *8, intent(in) :: xyztarg(ndtarg,ntarg)
      real *8, intent(in) :: xyzproxy(3,*)
      integer, intent(in) :: itargptr(npatches)
      integer, intent(in) :: ntargptr(npatches)
      
      external fker
      integer, intent(in) :: ndd,ndz,ndi
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      integer, intent(in) :: ipars(ndi)

      integer, intent(in) :: nqorder
      real *8, intent(in) :: rfac

      complex *16, intent(out) :: cintvals(nppols,ntarg)
c
cc      temporary variables
c
      real *8, allocatable :: quadcm(:,:,:)
      real *8, allocatable :: quadrad(:,:)
      integer, allocatable :: iquadreltmp(:,:)
      integer, allocatable :: iquadrel(:,:)
      integer, allocatable :: iquadrelall(:)
      real *8, allocatable :: qverts(:,:,:)
      integer, allocatable :: ichild_start(:)

      integer nquad,nquadmax,nlev,iquad,istart,i,j
      integer ier,itarg,jj,jstart,nlmax,npts
      integer iqquad,ii

      integer npmax

      real *8, allocatable :: uvsq(:,:),wts(:),uvtmp(:,:)
      real *8, allocatable :: umattmp(:,:),vmattmp(:,:)
      real *8, allocatable :: da(:)
      integer nqpols
      real *8, allocatable :: sigvals(:,:)
      real *8, allocatable :: srcvals(:,:),qwts(:)
      real *8, allocatable :: xyztargtmp(:,:)
      complex *16, allocatable :: fkervals(:),sigmatmp(:,:)
      real *8, allocatable :: rsigtmp(:,:)
      real *8 xyztmp(3)
      integer itmp
      complex *16 ima
      complex *16 alpha_c, beta_c

      character *1 transa,transb
      double precision :: alpha,beta
      integer lda,ldb,ldc
      

      data ima/(0.0d0,1.0d0)/





cc      max number of levels
c
      nlmax = 20
      allocate(quadcm(3,npatches,nquadmax),quadrad(npatches,nquadmax))
      allocate(qverts(2,3,nquadmax))
      allocate(iquadreltmp(ntarg,nquadmax))
      allocate(da(nquadmax),ichild_start(nquadmax))

      do i=1,nquadmax
        do j=1,ntarg
          iquadreltmp(j,i) = 0
        enddo
      enddo

      
      nquad = 0
      nlev = 0
      ier = 0

      do i=1,nquadmax
        ichild_start(i) = -1
      enddo

      

      if(ifp.eq.1) then
        call getquadtree(npatches,norder,npols,srccoefs,ntarg,
     1     xyzproxy,itargptr,ntargptr,nquadmax,nlmax,rfac,
     2     nquad,nlev,ichild_start,da,quadcm,quadrad,qverts,
     3     iquadreltmp,ier)
      else
        allocate(xyztargtmp(3,ntarg))

        do i=1,ntarg
          do j=1,3
            xyztargtmp(j,i) = xyztarg(j,i)
          enddo
        enddo
        call getquadtree(npatches,norder,npols,srccoefs,ntarg,
     1     xyztargtmp,itargptr,ntargptr,nquadmax,nlmax,rfac,nquad,
     2     nlev,ichild_start,da,quadcm,quadrad,qverts,iquadreltmp,ier)
        deallocate(xyztargtmp)
      endif


      

      if(ier.ne.0) then
        call prinf('Could not allocate quadrangle tree*',i,0)
        call prinf('Exiting without computing anything*',i,0)
        call prinf('Press 0 to continue*',i,0)
        call prinf('Press 1 to exit routine and continue*',i,0)
        call prinf('Press 2 to stop*',i,0)
        read *, ier0
        if(ier0.eq.2) stop
        if(ier0.eq.1) return
      endif

      allocate(iquadrel(nquad,ntarg),iquadrelall(nquad))

c
c       transpose the iquadrel array
c  

      do iquad=1,nquad
        do j=1,ntarg
          iquadrel(iquad,j) = iquadreltmp(j,iquad)
        enddo
      enddo


c
c       compute all the geometry info and koornwinder polynomials
c       at relevant quadrangle given by iquadrelall
c  

c
c       get quadrature nodes and weights on the base quadrangle
c       based on quadrature type
c

      if(intype.eq.1) then
         
         if(ttype.eq.'t'.or.ttype.eq.'T') then
           nqpols = (nqorder+1)*(nqorder+2)/2
         else if(ttype.eq.'f'.or.ttype.eq.'F') then
           nqpols = (nqorder+1)*(nqorder+1)
         else
           ier = 2
           cintvals = 0
           call prinf('invalid polynomial type*',i,0)
           call prinf('leaving without computing anything, ier=*',ier,1)
           return
         endif
         allocate(uvsq(2,nqpols),wts(nqpols))
         allocate(umattmp(nqpols,nqpols),vmattmp(nqpols,nqpols))
         itype = 1
         
         call polytens_exps_nd(2,ipoly,itype,nqorder+1,ttype,uvsq,
     1     umattmp,1,vmattmp,1,wts)
      
         deallocate(umattmp,vmattmp)
      endif


      if(intype.eq.2) then
        call squarearbq_pts(nqporder,nqpols)
        allocate(uvsq(2,nqpols),wts(nqpols))
        call squarearbq(nqorder,uvsq,wts,nqpols)
      endif

cc      call prinf('nqpols=*',nqpols,1)

      npmax = nquad*nqpols
      allocate(sigvals(npols,npmax))
      allocate(srcvals(12,npmax),qwts(npmax))

      allocate(sigmatmp(nppols,npmax),fkervals(npmax))
      allocate(rsigtmp(nppols,npmax))

      allocate(uvtmp(2,nqpols))
      alpha_c = 1
      beta_c = 0


      do iquad=1,nquad
        call mapuv_quad(qverts(1,1,iquad),nqpols,uvsq,uvtmp)
        istart = (iquad-1)*nqpols
        do i=1,nqpols
          ii = istart+i
          call polytens_pols_2d(ipoly,uvtmp(1,i),norder,ttype,
     1       sigvals(1,ii))
          call polytens_pols_2d(ipoly,uvtmp(1,i),nporder,ttype,
     1       rsigtmp(1,ii))
        enddo
      enddo


      do iquad=1,npatches
c
cc       for the current patch compute all the geometry info
c
        transa = 'N'
        transb = 'N'
        alpha = 1
        beta = 0
        lda = 9
        ldb = npols
        ldc = 12

        call dgemm_guru(transa,transb,9,npmax,npols,alpha,
     1     srccoefs(1,1,iquad),lda,sigvals,ldb,beta,srcvals,ldc)


c
c        compute all the quadrature weights
c

        do i=1,nquad
          istart = (i-1)*nqpols+1
          call get_norms_qwts_quad(nqpols,wts,srcvals(1,istart),
     1        da(i),qwts(istart))
        enddo

        do itarg=itargptr(iquad),itargptr(iquad)+ntargptr(iquad)-1
c
c           extract info of geometry on subquadrangles relevant
c           for this computation
c
          npts = 0
          do iqquad=1,nquad
            jstart = (iqquad-1)*nqpols
            if(iquadrel(iqquad,itarg).eq.1) then
              do i=1,nqpols
                jj = jstart+i
                ii = npts+i
                call fker(srcvals(1,jj),ndtarg,xyztarg(1,itarg),
     1             ndd,dpars,ndz,zpars,ndi,ipars,fkervals(ii))
                fkervals(ii) = fkervals(ii)*qwts(jj)
                do j=1,nppols
                  sigmatmp(j,ii) = rsigtmp(j,jj)
                enddo
              enddo
              npts = npts + nqpols
            endif
          enddo

c
c          TODO:  fix this to call mkl blas with single thread
c
          call zgemv_guru('n',nppols,npts,alpha_c,sigmatmp,nppols,
     1       fkervals,1,beta_c,cintvals(1,itarg),1)
        enddo
      enddo

      return
      end
c
c
c
c
      subroutine cquadints_adap(eps,intype,
     1     npatches,norder,npols,srccoefs,ndtarg,
     2     ntarg,xyztarg,itargptr,ntargptr,nporder,nppols,
     3     nquadmax,rat1,rat2,rsc1,rat1p,rat2p,rsc1p,
     4     fker,ndd,dpars,ndz,zpars,ndi,ipars,nqorder,cintvals)
c-------------
c  Task:
c    This subroutine computes the integrals
c
c      \int_{[-1,1]^2} K(x_{i},\rho_{m}(y)) B_{n}(y) J_{m}(y) dy \, , 
c
c    B_{n}(y) are either tensor product legendre polynomials or
c    Chebyshev polynomials on [-1,1]^2 of either total degree 
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
c    Use adaptive integration, 
c
c   CONTINUE FRO HERE
c    
c
c  Notes:
c    There are two parameters which determine the quadrature
c    nodes used on each quadrangle in the heirarchy,
c    intype, and nqorder
c    intype switches between chebyshev, legendre, and xiao gimbutas
c    nodes of full or total degree, and nqorder determines
c    the order of nodes used
c
c  Input arguments:
c    - eps: real *8
c        precision requested (unused in this routine)
c    - intype: integer
c        type of nodes to be used as quadrature nodes for each
c        quadrangle in the heirarchy:
c          intype = 1, use either tensor product chebyshev/legendre
c            nodes of total/full degree based on the input parameters
c            ipoly, and ttype
c          intype = 2, use xiao gimbutas nodes on the square
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
c    - nquadmax: integer
c        max number of quadrangles supported in the heirarchy
c        if total number of quadrangles required exceeds
c        nquadmax, subroutine exits with error code
c        ier = 4
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
c    - nqorder: integer
c        order of quadrature nodes to be used
c    - rfac: real *8
c        parameter for defining refinement criterion, quadrangles
c        refined if not separated from target by
c        rfac*r_{t} where r_{t} is radius of enclosing sphere
c        in quadrangle
c  Output arguments:
c    - cintvals: complex *16 (nppols,ntarg)
c        Integral against all basis functions from patches
c        to targets
c    - ier: integer
c        error code, ier = 0, implies successful execution
c        ier = 4, not enough quadrangles, try rerunning the
c          routine with more quads
c
c-------------------

c
      implicit none

c
cc     calling sequence variables
c
      real *8 eps
      integer intype
      integer npatches,norder,npols
      integer nporder,nppols
      real *8 srccoefs(9,npols,npatches)
      
      integer ntarg,ndtarg
      real *8 xyztarg(ndtarg,ntarg)
      integer itargptr(npatches)
      integer ntargptr(npatches)
      
      external fker
      integer ndd,ndz,ndi
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer ipars(ndi)

      real *8 rat1(2,0:norder)
      real *8 rat2(3,0:norder,0:norder)
      real *8 rsc1(0:norder,0:norder)

      real *8 rat1p(2,0:nporder)
      real *8 rat2p(3,0:nporder,0:nporder)
      real *8 rsc1p(0:nporder,0:nporder)

      integer nqorder

      complex *16 cintvals(nppols,ntarg)

c
c       tree variables
c
      integer nlmax,ltree
      real *8, allocatable :: tvs(:,:,:),da(:)
      integer, allocatable :: ichild_start(:)
      real *8, allocatable :: tvs2(:,:,:),da2(:)
      integer, allocatable :: ichild_start2(:)

      integer nquad,nquadmax,nlev,iquad,istart,i,j,k
      integer ier,itarg,jj,jstart,npts
      integer iqquad,ii


      integer npmax

      real *8, allocatable :: uvsq(:,:),wts(:),uvtmp(:,:)
      real *8, allocatable :: umattmp(:,:),vmattmp(:,:)
      integer nqpols
      real *8, allocatable :: sigvals(:,:),sigvalsdens(:,:)
      real *8, allocatable :: srcvals(:,:),qwts(:)
      real *8, allocatable :: sigvals2(:,:),sigvalsdens2(:,:)
      real *8, allocatable :: srcvals2(:,:),qwts2(:)
      integer itmp

      character *1 transa,transb
      real *8 alpha,beta
      integer lda,ldb,ldc
      integer nn1,nn2,nn3,nn4,npmax0,ntmaxuse,ntmaxuse0
      
      
c
c      get the tree
c
      nlmax = 20


      ntmaxuse = nquadmax
c
c       for each quadrangle, we just store three pieces
c       of info
c         quadrangle vertices
c         area of quadrangle
c         ichild_start = index for first child, 
c         
c
      allocate(ichild_start(nquadmax),tvs(2,3,nquadmax))
      allocate(da(nquadmax))

      do i=1,nquadmax
        ichild_start(i) = -1
        da(i) = 0
        do j=1,3
          do k=1,2
            tvs(k,j,i) = 0
          enddo
        enddo
      enddo

      da(1) = 1.0d0
      
      tvs(1,1,1) = 0
      tvs(2,1,1) = 0

      tvs(1,2,1) = 1
      tvs(2,2,1) = 0

      tvs(1,3,1) = 0
      tvs(2,3,1) = 1

c
c       get quadrature nodes and weights on the base quadrangle
c       based on quadrature type
c

      if(intype.eq.1) then
         nqpols = (nqorder+1)*(nqorder+2)/2
         allocate(uvsq(2,nqpols),wts(nqpols))
         allocate(umattmp(nqpols,nqpols),vmattmp(nqpols,nqpols))
      
         call vioreanu_simplex_quad(nqorder,nqpols,uvsq,umattmp,
     1      vmattmp,wts)
         
         deallocate(umattmp,vmattmp)
      endif

      if(intype.eq.2) then
        call quadasymq_pts(nqorder,nqpols)
        allocate(uvsq(2,nqpols),wts(nqpols))
        call quadasymq(nqorder,tvs(1,1,1),tvs(1,2,1),tvs(1,3,1),
     1         uvsq,wts,nqpols)    
      endif

      allocate(uvtmp(2,nqpols))
     
      npmax = nquadmax*nqpols
      allocate(sigvals(npols,npmax))
      allocate(sigvalsdens(nppols,npmax))
      allocate(srcvals(12,npmax),qwts(npmax))

c
c      current number of quadrangles in the adaptive structure
c
      nquad = 1
c
c        intialize sigvals for root quadrangle
c


      call mapuv_quad(tvs(1,1,1),nqpols,uvsq,uvtmp)
      do i=1,nqpols
        call koornf_pols(uvtmp(1,i),norder,npols,
     1      sigvals(1,i),rat1,rat2,rsc1)
        call koornf_pols(uvtmp(1,i),nporder,nppols,
     1      sigvalsdens(1,i),rat1p,rat2p,rsc1p)
      enddo


      do iquad=1,npatches
c
cc       for the current patch compute geometry info for base quadrangle 
c
        transa = 'N'
        transb = 'N'
        alpha = 1
        beta = 0
        lda = 9
        ldb = npols
        ldc = 12

        call dgemm_guru(transa,transb,9,nqpols,npols,alpha,
     1     srccoefs(1,1,iquad),lda,sigvals,ldb,beta,srcvals,ldc)
        call get_norms_qwts_quad(nqpols,wts,srcvals,da,qwts)

        do itarg=itargptr(iquad),itargptr(iquad)+ntargptr(iquad)-1


 1111     continue
          ier = 0
          call quadadap(eps,nqorder,nqpols,nlmax,ntmaxuse,nquad,
     1          ichild_start,tvs,da,uvsq,wts, 
     1          norder,npols,srccoefs(1,1,iquad),
     2          npmax,srcvals,
     2          qwts,sigvals,nporder,nppols,sigvalsdens,ndtarg,
     3          xyztarg(1,itarg),
     3          rat1,rat2,rsc1,
     3          rat1p,rat2p,rsc1p,
     3          fker,ndd,dpars,ndz,zpars,ndi,ipars,
     3          cintvals(1,itarg),ier)
           if(ier.eq.4) then
             ntmaxuse0 = ntmaxuse*4
             npmax0 = ntmaxuse0*nqpols
             allocate(sigvals2(npols,npmax))
             allocate(sigvalsdens2(nppols,npmax))
             allocate(srcvals2(12,npmax),qwts2(npmax))
             allocate(ichild_start2(ntmaxuse),tvs2(2,3,ntmaxuse))
             allocate(da2(ntmaxuse))
             nn1 = npols*npmax
             nn2 = nppols*npmax
             nn3 = 12*npmax
             nn4 = nquad*6
             call dcopy(nn1,sigvals,1,sigvals2,1)
             call dcopy(nn2,sigvalsdens,1,sigvalsdens2,1)
             call dcopy(nn3,srcvals,1,srcvals2,1)
             call dcopy(npmax,qwts,1,qwts2,1)
             do ii=1,nquad
               ichild_start2(ii) = ichild_start(ii)
             enddo
             call dcopy(nn4,tvs,1,tvs2,1)
             call dcopy(nquad,da,1,da2,1)


             deallocate(sigvals,sigvalsdens,srcvals,qwts,ichild_start)
             deallocate(tvs,da)


             allocate(sigvals(npols,npmax0))
             allocate(sigvalsdens(nppols,npmax0))
             allocate(srcvals(12,npmax0),qwts(npmax0))
             allocate(ichild_start(ntmaxuse0),tvs(2,3,ntmaxuse0))
             allocate(da(ntmaxuse0))

             do ii=1,ntmaxuse0
               ichild_start(ii) = -1
             enddo

             call dcopy(nn1,sigvals2,1,sigvals,1)
             call dcopy(nn2,sigvalsdens2,1,sigvalsdens,1)
             call dcopy(nn3,srcvals2,1,srcvals,1)
             call dcopy(npmax,qwts2,1,qwts,1)
             do ii=1,nquad
               ichild_start(ii) = ichild_start2(ii)
             enddo
             call dcopy(nn4,tvs2,1,tvs,1)
             call dcopy(nquad,da2,1,da,1)

             npmax = npmax0
             ntmaxuse = ntmaxuse0
             call prinf('restarting adaptive inetgration with nquad=*',
     1             ntmaxuse,1)


             deallocate(sigvals2,sigvalsdens2,srcvals2,ichild_start2)
             deallocate(tvs2,da2,qwts2)
             goto 1111
           endif
        enddo

        do i=1,nquad
          ichild_start(i) = -1
        enddo
      enddo

      return
      end

c
c
c
c
c
      subroutine quadadap(eps,m,kpols,nlmax,ntmax,nquad,
     1             ichild_start,tvs,da,uvsq,wts,
     1             norder,npols,srccoefs,npmax,srcvals,
     2             qwts,sigvals,nporder,nppols,sigvalsdens,
     3             ndtarg,xt,
     3             rat1,rat2,rsc1,
     3             rat1p,rat2p,rsc1p,
     3             fker,ndd,dpars,ndz,zpars,ndi,
     3             ipars,cintall,ier)

c
c       this subroutine adaptively computes the integral
c        of the functions
c   
c        \int_{T} 1/|xt- y(u,v)| K_{n,m}(u,v) |J(u,v)| du dv
c        n,m = 1,2\ldots npols
c
c        K_{n,m}(u,v) are the koornwinder polynomials on the 
c        standard simplex (0,0)-(1,0)-(0,1)
c
c         |J(u,v)| = |dy/du \times dy/dv|
c        
c
c        relevant parameters on a quad refined
c        grid along the way
c
c        DOCUMENTATION NEEDS UPDATING:
c
c
c        IN:
c        eps - precision requested
c        m - quadrature order
c        kpols - number of quarature nodes (look up quadasymq 
c                 for getting kpols(m))
c        nlmax - max level of refinement for geometry
c        ntmax - max number of quadrangles
c        nquad - current number of quadrangles in adaptive structure
c        ichild_start(i) - first child of quadrangle i
c        tvs(2,3,ntmax) - vertices of hierarchy of quadrangles
c        da(ntmax) - area of quadrangles
c        uvsq(kpols) - integration nodes on standard quadrangle
c        wts(kpols) - integration weights on standard quadrangle
c        npols - total number of koornwinder polynomials to be integrated
c        norder - order of discretization of the surface
c        npols - (norder+1)*(norder+2)/2: total number of 
c                koornwinder polynomials to be integrated
c        srccoefs(9,npols) - xyz coefficients of koornwinder expansion 
c                             of current surface + derivative info
c 
c        npmax - max number of points = ntmax*kpols
c        srczvals(12,npmax) - geometry info on heirarchy of meshes
c        qwts(npmax) - quadrature weights 
c        sigvals(npols,npmax) - 
c                   koornwinder polynomials computed along the adaptive grid
c        
c        OUT:
c        cintall(npols) - computed integral
c        ier - error code
c           ier = 0, successful execution
c           ier = 4, too few quadrangles, try with more quadrangles
c
c         

      implicit real *8 (a-h,o-z)
      integer, allocatable :: istack(:)
      integer ichild_start(ntmax)
      real *8 da(ntmax)
      real *8 tvs(2,3,ntmax), uvsq(2,kpols),wts(kpols)
      integer nproclist0, nproclist
      integer idone
      real *8 srccoefs(9,npols)
      real *8 sigvals(npols,npmax)
      real *8 sigvalsdens(nppols,npmax)
      real *8 srcvals(12,*),qwts(npmax)
      complex *16, allocatable :: xkernvals(:)
      real *8 xt(ndtarg)
      complex *16 cintall(nppols),fval
      complex *16, allocatable :: ctmp(:)
      complex *16, allocatable :: cvals(:,:)

      integer ndd,ndz,ndi
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer ipars(ndi)
      integer ier

      real *8 rat1(2,0:norder),rat2(3,0:norder,0:norder)
      real *8 rsc1(0:norder,0:norder)

      external fker

c
c         for historic reasons
c
      ksigpols = nppols
      allocate(ctmp(nppols))
      allocate(istack(2*ntmax))
      allocate(cvals(ksigpols,ntmax))

      nfunev = 0

      do i=1,ksigpols
        cvals(i,1) = 0
      enddo

      allocate(xkernvals(npmax))

c
cc      compute integral at level 0
c
      do i=1,kpols
         call fker(srcvals(1,i),ndtarg,xt,ndd,dpars,ndz,zpars,ndi,
     1      ipars,fval)
         xkernvals(i) = fval*qwts(i)
         do j=1,ksigpols
            cvals(j,1) = cvals(j,1)+xkernvals(i)*sigvalsdens(j,i)
         enddo
      enddo

      nfunev = nfunev + kpols

      
      do i=1,ksigpols
         cintall(i) = cvals(i,1)
      enddo

      nproclist0 = 1
      istack(1) = 1


      call quadadap_main(eps,kpols,nlmax,ntmax,nquad,ichild_start,
     1      tvs,da,uvsq,wts,norder,npols,srccoefs,
     2      npmax,srcvals,qwts,sigvals,nporder,nppols,
     3      sigvalsdens,ndtarg,xt,rat1,rat2,rsc1,
     3      rat1p,rat2p,rsc1p,
     3      fker,ndd,dpars,
     3      ndz,zpars,ndi,ipars,cvals,istack,nproclist0,
     4      xkernvals,cintall,ier)
      
      return
      end
c
c
c
c
c
       
      subroutine quadadap_main(eps,kpols,nlmax,ntmax,nquad,ichild_start,
     1      tvs,da,uvsq,wts,norder,npols,srccoefs,
     2      npmax,srcvals,qwts,sigvals,nporder,nppols,sigvalsdens,
     3      ndtarg,xt,rat1,rat2,rsc1,
     3      rat1p,rat2p,rsc1p,
     3      fker,ndd,dpars,ndz,
     3      zpars,ndi,ipars,cvals,istack,nproclist0,xkernvals,
     4      cintall,ier)
      

      implicit real *8 (a-h,o-z)
      integer istack(*),nproclist0
      integer ichild_start(ntmax)
      integer nporder,nppols
      real *8 da(ntmax)
      real *8 tvs(2,3,ntmax), uvsq(2,kpols),wts(kpols)
      integer  nproclist
      integer idone
      real *8 srccoefs(9,npols)
      real *8 sigvals(npols,npmax)
      real *8 sigvalsdens(nppols,npmax)
      real *8 qwts(npmax)
      complex *16 xkernvals(npmax)
      real *8 xt(ndtarg)
      real *8 srcvals(12,*)
      complex *16 cintall(nppols),fval,ctmp(nppols)
      complex *16 cvals(nppols,ntmax)

      integer ndd,ndz,ndi
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer ipars(ndi)

      real *8 rat1(2,0:norder),rat2(3,0:norder,0:norder)
      real *8 rsc1(0:norder,0:norder)

      real *8 rat1p(2,0:nporder),rat2p(3,0:nporder,0:nporder)
      real *8 rsc1p(0:nporder,0:nporder)


      real *8, allocatable :: uvtmp(:,:)
      character *1 transa,transb
      integer lda,ldb,ldc
      external fker
      
      allocate(uvtmp(2,kpols))

c
c         for historic reasons
c
      ksigpols = nppols
      kfine = 4*kpols

      ier = 0


      do ilev=0,nlmax
        idone = 1
        nproclist = 0
        
        do iproc = 1,nproclist0
          iquad = istack(iproc)

c
c           check to see if quadrangle already has 
c           children, if not, set children
c           and compute necessary info
c
          if(ichild_start(iquad).eq.-1) then

c
c            current quadrangle doesn't have children,
c            compute necessary info
c

            if(nquad+4.gt.ntmax) then
c               print *, "Too many quadrangles in cquadadap"
c               print *, "Exiting without computing anything"
               ier = 4
               return
            endif
            
            ichild_start(iquad) = nquad+1
            call getquadchildren(tvs(1,1,iquad),tvs(1,1,nquad+1),
     1             tvs(1,1,nquad+2),tvs(1,1,nquad+3),tvs(1,1,nquad+4))

            rr = 0.25d0*da(iquad)
            do j=nquad+1,nquad+4
              da(j) = rr
              call mapuv_quad(tvs(1,1,j),kpols,uvsq,uvtmp)
              istart = (j-1)*kpols+1
              do i=1,kpols
                 ii = istart+i-1
                call koornf_pols(uvtmp(1,i),norder,npols,sigvals(1,ii),
     1             rat1,rat2,rsc1)
                call koornf_pols(uvtmp(1,i),nporder,nppols,
     1             sigvalsdens(1,ii),rat1p,rat2p,rsc1p)
              enddo
              
              transa = 'N'
              transb = 'N'
              alpha = 1
              beta = 0
              lda = 9
              ldb = npols
              ldc = 12

              call dgemm_guru(transa,transb,9,kpols,npols,alpha,
     1           srccoefs,lda,sigvals(1,istart),ldb,beta,
     2           srcvals(1,istart),ldc)
              call get_norms_qwts_quad(kpols,wts,srcvals(1,istart),
     1           rr,qwts(istart))

            enddo
            nquad = nquad+4
          endif
        
c
cc           compute xkernvals
c
          iquadc1 = ichild_start(iquad)
          istart = (iquadc1-1)*kpols
          do j=1,kfine
            jj=j+istart
            call fker(srcvals(1,jj),ndtarg,xt,ndd,dpars,
     1         ndz,zpars,ndi,ipars,fval)
            xkernvals(jj) = fval*qwts(jj)
          enddo
cc          call prin2('xkernvals=*',xkernvals(istart+1),kfine)

          nfunev = nfunev + kfine

c
cc         subtract conquadbution of current 
c          quadrangle
          do isig=1,ksigpols
            cintall(isig) = cintall(isig)-cvals(isig,iquad)
            ctmp(isig) = 0
          enddo

c
cc        add in conquadbutions of children quadrangles
c
          do iquadc=iquadc1,iquadc1+3
            do isig=1,ksigpols
              cvals(isig,iquadc) = 0
            enddo

            istart = (iquadc-1)*kpols
            do k=1,kpols
              ii = istart+k
              do isig=1,ksigpols
                cvals(isig,iquadc) = cvals(isig,iquadc)+xkernvals(ii)*
     1                                  sigvalsdens(isig,ii)
              enddo
            enddo
              
            do isig=1,ksigpols
              cintall(isig) = cintall(isig) + cvals(isig,iquadc)
              ctmp(isig) = ctmp(isig)+cvals(isig,iquadc)
            enddo
          enddo

c
cc        compare integral of children to integral of parent
c         to determine if the children need to be refined further
c
          errmax = 0
          do isig=1,ksigpols
            if(abs(ctmp(isig)-cvals(isig,iquad)).gt.errmax) 
     1          errmax = abs(ctmp(isig)-cvals(isig,iquad))
          enddo

          if(errmax.gt.eps) then
            idone = 0
            
            do j=1,4
              istack(nproclist0+nproclist+j) = iquadc1+j-1
            enddo
            nproclist = nproclist+4
          endif
c
cc        end of looping over all quadrangles at current stage
        enddo
cc         if idone is still 1, that means that no more refinement
c          is needed
         if(idone.eq.1) goto 1111
        do i=1,nproclist
          istack(i) = istack(nproclist0+i)
        enddo
        nproclist0 = nproclist
      enddo
 1111 continue


      return
      end
c
c
c
c
c
c
c
      subroutine cquadints_comb(eps,intype,
     1     npatches,norder,npols,srccoefs,ndtarg,
     2     ntarg,xyztarg,ifp,xyzproxy,
     3     itargptr,ntargptr,nporder,nppols,nquadmax,rat1,rat2,rsc1,
     3     rat1p,rat2p,rsc1p,
     3     fker,ndd,dpars,ndz,zpars,ndi,ipars,
     4     nqorder,rfac,cintvals)
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
c        this code is the same as cquadints on all counts except
c        for choosing when to subdivide a quadrangle
c        
c        here a quadrangle is subdivided when xyzproxy associated
c        with the quadrangle is separated from the centroid
c        of the quadrangle by at least rfac*r_{t} where r_{t}
c        is the radius of the smallest sphere  
c
c        input arguments:
c        eps:     requested precision (Currently not being used
c                 in any way, later, nqorder, and rfac will
c                 be set internally based on eps)
c
c        intype:   quadrature node type
c                   intype = 1, rokhlin vioreanu nodes
c                   intype = 2, xiao gimbutas nodes
c   
c        npatches: number of patches
c        norder: order of discretization nodes on the patches
c        npols = (norder+1)*(norder+2)/2: number of discretization 
c                   nodes on each patch
c        srccoefs(3,npols,npatches): coefficients
c                 of koornwinder expansion of xyz coordinates
c                 and  of dxyz/du and dxyz/dv
c
c        ndtarg - dimension of each target point (
c          ndtarg must at least be 3, and correspond to location
c          of target points)
c        ntarg - total number of target points
c        xyztarg(ndtarg,ntarg) - location of target points 
c                       (currently may have repeats, 
c                        being slightly wasteful in
c                        memory, but see comment (1))
c
c        ifp - flag for using proxy target locations
c               proxy target locations will be used if ifp = 1
c               else original target locations will be used
c        xyzproxy(3,*) - proxy locations for measuring distances
c
c        itargptr(npatches) - xyztargs(:,itargptr(i):itargptr(i)+ntargptr(i)-1)
c                       are the relevant set of targets for patch i
c        ntargptr(npatches) - number of relevant targets for patch i
c        nquadmax - max number of quadrangles to be used on base quadrangle
c        fker - function handle for evaluating the kernel k
c 
c               
c               expected calling sequence
c               fker(x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
c               the output is assumed to be complex for the time
c               being
c
c         dpars(*) - real parameters for the fker routine
c         zpars(*) - complex parameters for the fker routine
c         ipars(*) - integer parameters for the fker routine
c         nqorder - order of quadrature nodes on each subquadrangle
c                   to be used
c         rfac - distance criterion for deciding whether a quadrangle
c                is in the far-field or not
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
      integer intype,ifp
      integer npatches,norder,npols
      integer nporder,nppols
      real *8 srccoefs(9,npols,npatches)
      
      integer ntarg,ndtarg
      real *8 xyztarg(ndtarg,ntarg)
      real *8 xyzproxy(3,*)
      integer itargptr(npatches)
      integer ntargptr(npatches)
      
      external fker
      integer ndd,ndz,ndi
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer ipars(ndi)

      integer nqorder
      real *8 rfac

      real *8 rat1(2,0:norder),rat2(3,0:norder,0:norder)
      real *8 rsc1(0:norder,0:norder)

      real *8 rat1p(2,0:nporder),rat2p(3,0:nporder,0:nporder)
      real *8 rsc1p(0:nporder,0:nporder)

      complex *16 cintvals(nppols,ntarg)
      

c
cc      temporary variables
c
      real *8, allocatable :: quadcm(:,:,:)
      real *8, allocatable :: quadrad(:,:)
      integer, allocatable :: iquadreltmp(:,:)
      integer, allocatable :: iquadrel(:,:)
      integer, allocatable :: iquadrelall(:)
      real *8, allocatable :: tverts(:,:,:)
      integer, allocatable :: ichild_start(:)
      integer, allocatable :: ichild_start0(:)

      integer nquad,nquadmax,nlev,iquad,istart,i,j
      integer ier,itarg,jj,jstart,nlmax,npts
      integer iqquad,ii

      integer npmax

      real *8, allocatable :: uvsq(:,:),wts(:),uvtmp(:,:)
      real *8, allocatable :: umattmp(:,:),vmattmp(:,:)
      real *8, allocatable :: da(:)
      integer nqpols
      real *8, allocatable :: sigvals(:,:)
      real *8, allocatable :: sigvalsdens(:,:)
      real *8, allocatable :: srcvals(:,:)
      real *8, allocatable :: qwts(:)
      complex *16, allocatable :: sigmatmp(:,:)
      real *8, allocatable :: xyztargtmp(:,:)
      real *8 xyztmp(3)
      integer itmp

      complex *16, allocatable :: xkernvals(:)
      integer, allocatable :: istack(:)

      complex *16, allocatable :: cvals(:,:)
      
      integer nproclist0,nproclist,nquad0,npts0
      complex *16 zz

      
      complex *16 ima

      character *1 transa,transb
      
      real *8 alpha,beta
      integer lda,ldb,ldc

      data ima/(0.0d0,1.0d0)/

      allocate(cvals(nppols,nquadmax))
      allocate(istack(2*nquadmax))





cc      max number of levels
c
      nlmax = 20
      allocate(quadcm(3,npatches,nquadmax),
     1   quadrad(npatches,nquadmax),tverts(2,3,nquadmax))
      allocate(iquadreltmp(ntarg,nquadmax))
      allocate(da(nquadmax),ichild_start(nquadmax))
      allocate(ichild_start0(nquadmax))

      do i=1,nquadmax
        do j=1,ntarg
          iquadreltmp(j,i) = 0
        enddo
        
      enddo

      
      nquad = 0
      nlev = 0
      ier = 0

      do i=1,nquadmax
        ichild_start(i) = -1
      enddo

      

      if(ifp.eq.1) then
        call getquadtree(npatches,norder,npols,srccoefs,ntarg,xyzproxy,
     1     itargptr,ntargptr,nquadmax,nlmax,rfac,nquad,nlev,ichild_start,
     2     da,quadcm,quadrad,tverts,iquadreltmp,ier)
      else
        allocate(xyztargtmp(3,ntarg))
        do i=1,ntarg
          do j=1,3
            xyztargtmp(j,i) = xyztarg(j,i)
          enddo
        enddo
        call getquadtree(npatches,norder,npols,srccoefs,ntarg,
     1     xyztargtmp,
     1     itargptr,ntargptr,nquadmax,nlmax,rfac,nquad,nlev,ichild_start,
     2     da,quadcm,quadrad,tverts,iquadreltmp,ier)
        deallocate(xyztargtmp)
      endif

      nquad0 = nquad

      do i=1,nquad0
        ichild_start0(i) = ichild_start(i)
      enddo


      

      if(ier.ne.0) then
        call prinf('Could not allocate quadrangle tree*',i,0)
        call prinf('Exiting without computing anything*',i,0)
        call prinf('Press 0 to continue*',i,0)
        call prinf('Press 1 to exit routine and continue*',i,0)
        call prinf('Press 2 to stop*',i,0)
        read *, ier
        if(ier.eq.2) stop
        if(ier.eq.1) return
      endif

      allocate(iquadrel(nquad0,ntarg),iquadrelall(nquad0))

c
c       transpose the iquadrel array
c  

      do iquad=1,nquad0
        do j=1,ntarg
          iquadrel(iquad,j) = iquadreltmp(j,iquad)
        enddo
      enddo


c
c       compute all the geometry info and koornwinder polynomials
c       at relevant quadrangle given by iquadrelall
c  

c
c       get quadrature nodes and weights on the base quadrangle
c       based on quadrature type
c

      if(intype.eq.1) then
         nqpols = (nqorder+1)*(nqorder+2)/2
         allocate(uvsq(2,nqpols),wts(nqpols))
         allocate(umattmp(nqpols,nqpols),vmattmp(nqpols,nqpols))
      
         call vioreanu_simplex_quad(nqorder,nqpols,uvsq,umattmp,
     1      vmattmp,wts)
         
         deallocate(umattmp,vmattmp)
      endif


      if(intype.eq.2) then
        call quadasymq_pts(nqorder,nqpols)
        allocate(uvsq(2,nqpols),wts(nqpols))
        call quadasymq(nqorder,tverts(1,1,1),tverts(1,2,1),tverts(1,3,1),
     1         uvsq,wts,nqpols)    
      endif

cc      call prinf('nqpols=*',nqpols,1)

      npts0 = nquad0*nqpols
      npmax = nquadmax*nqpols
      allocate(sigvals(npols,npmax))
      allocate(sigvalsdens(nppols,npmax))
      allocate(srcvals(12,npmax),qwts(npmax))

      allocate(xkernvals(npmax))

      allocate(uvtmp(2,nqpols))


      do iquad=1,nquad0
        call mapuv_quad(tverts(1,1,iquad),nqpols,uvsq,uvtmp)
        istart = (iquad-1)*nqpols
        do i=1,nqpols
          ii = istart+i
          call koornf_pols(uvtmp(1,i),norder,npols,sigvals(1,ii),
     1           rat1,rat2,rsc1)
          call koornf_pols(uvtmp(1,i),nporder,nppols,sigvalsdens(1,ii),
     1           rat1p,rat2p,rsc1p)
        enddo
      enddo


      do iquad=1,npatches
        nquad = nquad0
c
cc       for the current patch compute all the geometry info
c

        transa = 'N'
        transb = 'N'
        alpha = 1.0d0
        beta = 0.0d0
        lda = 9
        ldb = npols
        ldc = 12

        call dgemm_guru(transa,transb,9,npts0,npols,alpha,
     1     srccoefs(1,1,iquad),lda,sigvals,ldb,beta,srcvals,ldc)



c
c        compute all the quadrature weights
c

        do i=1,nquad
          istart = (i-1)*nqpols+1
          call get_norms_qwts_quad(nqpols,wts,srcvals(1,istart),
     1        da(i),qwts(istart))
        enddo

        do itarg=itargptr(iquad),itargptr(iquad)+ntargptr(iquad)-1
          do i=1,nppols
            cintvals(i,itarg) = 0
          enddo
c
c           extract info of geometry on subquadrangles relevant
c           for this computation
c
          nproclist0 = 0
          do iqquad=1,nquad0
            jstart = (iqquad-1)*nqpols
            if(iquadrel(iqquad,itarg).eq.1) then
              nproclist0 = nproclist0 + 1
              istack(nproclist0) = iqquad
c
c
c              compute the integral due to iquad
c

              do i=1,nppols
                cvals(i,iqquad) = 0
              enddo
              do i=1,nqpols
                jj = jstart+i
                call fker(srcvals(1,jj),ndtarg,xyztarg(1,itarg),
     1              ndd,dpars,ndz,zpars,ndi,ipars,xkernvals(jj))
                xkernvals(jj) = xkernvals(jj)*qwts(jj)
                
                do j=1,nppols
                  zz = xkernvals(jj)*sigvalsdens(j,jj)
                  cvals(j,iqquad) = cvals(j,iqquad)+zz
                  cintvals(j,itarg)=cintvals(j,itarg)+zz
                enddo
              enddo
            endif
          enddo
c
c           done with initial computation of 
c

          ier = 0
          call quadadap_main(eps,nqpols,nlmax,nquadmax,nquad,ichild_start,
     1      tverts,da,uvsq,wts,norder,npols,srccoefs(1,1,iquad),
     2      npmax,srcvals,qwts,sigvals,nporder,nppols,sigvalsdens,
     3      ndtarg,xyztarg(1,itarg),
     3      rat1,rat2,rsc1,rat1p,rat2p,rsc1p,
     3      fker,ndd,dpars,ndz,zpars,ndi,ipars,cvals,istack,nproclist0,
     4      xkernvals,cintvals(1,itarg),ier)

      
        enddo

        do i=1,nquad0
          ichild_start(i) = ichild_start0(i)
        enddo
        do i=nquad0+1,nquad
          ichild_start(i) = -1
        enddo
      enddo

      return
      end
c
c
c
c
c
c
c
c
      subroutine cquadints_dist_vec(eps,intype,
     1     npatches,norder,npols,srccoefs,ndtarg,
     2     ntarg,xyztarg,ifp,xyzproxy,
     3     itargptr,ntargptr,nporder,nppols,nquadmax,rat1,rat2,rsc1,
     3     rat1p,rat2p,rsc1p,
     3     fker,nd,ndd,dpars,ndz,zpars,ndi,ipars,
     4     nqorder,rfac,cintvals)
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
c        this code is the same as cquadints on all counts except
c        for choosing when to subdivide a quadrangle
c        
c        here a quadrangle is subdivided when xyzproxy associated
c        with the quadrangle is separated from the centroid
c        of the quadrangle by at least rfac*r_{t} where r_{t}
c        is the radius of the smallest sphere  
c
c        input arguments:
c        eps:     requested precision (Currently not being used
c                 in any way, later, nqorder, and rfac will
c                 be set internally based on eps)
c
c        intype:   quadrature node type
c                   intype = 1, rokhlin vioreanu nodes
c                   intype = 2, xiao gimbutas nodes
c   
c        npatches: number of patches
c        norder: order of discretization nodes on the patches
c        npols = (norder+1)*(norder+2)/2: number of discretization 
c                   nodes on each patch
c        srccoefs(9,npols,npatches): coefficients
c                 of koornwinder expansion of xyz coordinates
c                 and dxyz/du, dxyz/dv
c
c        ndtarg - dimension of each target point (
c          ndtarg must at least be 3, and correspond to location
c          of target points)
c        ntarg - total number of target points
c        xyztarg(ndtarg,ntarg) - location of target points 
c                       (currently may have repeats, 
c                        being slightly wasteful in
c                        memory, but see comment (1))
c
c        ifp - flag for using proxy target locations
c               proxy target locations will be used if ifp = 1
c               else original target locations will be used
c        xyzproxy(3,*) - proxy locations for measuring distances
c
c        itargptr(npatches) - xyztargs(:,itargptr(i):itargptr(i)+ntargptr(i)-1)
c                       are the relevant set of targets for patch i
c        ntargptr(npatches) - number of relevant targets for patch i
c        nquadmax - max number of quadrangles to be used on base quadrangle
c        fker - function handle for evaluating the kernel k
c 
c               expected calling sequence
c               fker(x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
c               
c               the output is assumed to be nd complex numbers
c
c         nd - number of kernels to be integrated
c         dpars(*) - real parameters for the fker routine
c         zpars(*) - complex parameters for the fker routine
c         ipars(*) - integer parameters for the fker routine
c         nqorder - order of quadrature nodes on each subquadrangle
c                   to be used
c         rfac - distance criterion for deciding whether a quadrangle
c                is in the far-field or not
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
      integer intype,ifp,nd
      integer npatches,norder,npols
      integer nporder,nppols
      real *8 srccoefs(9,npols,npatches)
      
      integer ntarg,ndtarg
      real *8 xyztarg(ndtarg,ntarg)
      real *8 xyzproxy(3,*)
      integer itargptr(npatches)
      integer ntargptr(npatches)
      
      external fker
      integer ndd,ndi,ndz
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer ipars(ndi)

      integer nqorder
      real *8 rfac
      real *8 rat1(2,0:norder),rat2(3,0:norder,0:norder)
      real *8 rsc1(0:norder,0:norder)

      real *8 rat1p(2,0:nporder),rat2p(3,0:nporder,0:nporder)
      real *8 rsc1p(0:nporder,0:nporder)

      complex *16 cintvals(nd,nppols,ntarg)
      

c
cc      temporary variables
c
      real *8, allocatable :: quadcm(:,:,:)
      real *8, allocatable :: quadrad(:,:)
      integer, allocatable :: iquadreltmp(:,:)
      integer, allocatable :: iquadrel(:,:)
      integer, allocatable :: iquadrelall(:)
      real *8, allocatable :: tverts(:,:,:)
      integer, allocatable :: ichild_start(:)

      integer nquad,nquadmax,nlev,iquad,istart,i,j
      integer ier,itarg,jj,jstart,nlmax,npts
      integer iqquad,ii

      integer npmax

      real *8, allocatable :: uvsq(:,:),wts(:),uvtmp(:,:)
      real *8, allocatable :: umattmp(:,:),vmattmp(:,:)
      real *8, allocatable :: da(:)
      integer nqpols
      real *8, allocatable :: sigvals(:,:)
      real *8, allocatable :: srcvals(:,:),qwts(:)
      real *8, allocatable :: xyztargtmp(:,:)
      complex *16, allocatable :: fkervals(:,:),sigmatmp(:,:)
      real *8, allocatable :: rsigtmp(:,:)
      real *8 xyztmp(3)
      integer itmp,idim
      complex *16 ima

      character *1 transa,transb
      double precision :: alpha,beta
      complex *16 alpha_c, beta_c
      integer lda,ldb,ldc
      

      data ima/(0.0d0,1.0d0)/





cc      max number of levels
c
      nlmax = 20
      allocate(quadcm(3,npatches,nquadmax),
     1   quadrad(npatches,nquadmax),tverts(2,3,nquadmax))
      allocate(iquadreltmp(ntarg,nquadmax))
      allocate(da(nquadmax),ichild_start(nquadmax))

      do i=1,nquadmax
        do j=1,ntarg
          iquadreltmp(j,i) = 0
        enddo
      enddo

      
      nquad = 0
      nlev = 0
      ier = 0

      do i=1,nquadmax
        ichild_start(i) = -1
      enddo

      

      if(ifp.eq.1) then
        call getquadtree(npatches,norder,npols,srccoefs,ntarg,xyzproxy,
     1     itargptr,ntargptr,nquadmax,nlmax,rfac,nquad,nlev,ichild_start,
     2     da,quadcm,quadrad,tverts,iquadreltmp,ier)
      else
        allocate(xyztargtmp(3,ntarg))
        do i=1,ntarg
          do j=1,3
            xyztargtmp(j,i) = xyztarg(j,i)
          enddo
        enddo

        call getquadtree(npatches,norder,npols,srccoefs,ntarg,xyztargtmp,
     1     itargptr,ntargptr,nquadmax,nlmax,rfac,nquad,nlev,ichild_start,
     2     da,quadcm,quadrad,tverts,iquadreltmp,ier)
        deallocate(xyztargtmp)
      endif


      

      if(ier.ne.0) then
        call prinf('Could not allocate quadrangle tree*',i,0)
        call prinf('Exiting without computing anything*',i,0)
        call prinf('Press 0 to continue*',i,0)
        call prinf('Press 1 to exit routine and continue*',i,0)
        call prinf('Press 2 to stop*',i,0)
        read *, ier
        if(ier.eq.2) stop
        if(ier.eq.1) return
      endif

      allocate(iquadrel(nquad,ntarg),iquadrelall(nquad))

c
c       transpose the iquadrel array
c  

      do iquad=1,nquad
        do j=1,ntarg
          iquadrel(iquad,j) = iquadreltmp(j,iquad)
        enddo
      enddo


c
c       compute all the geometry info and koornwinder polynomials
c       at relevant quadrangle given by iquadrelall
c  

c
c       get quadrature nodes and weights on the base quadrangle
c       based on quadrature type
c

      if(intype.eq.1) then
         nqpols = (nqorder+1)*(nqorder+2)/2
         allocate(uvsq(2,nqpols),wts(nqpols))
         allocate(umattmp(nqpols,nqpols),vmattmp(nqpols,nqpols))
      
         call vioreanu_simplex_quad(nqorder,nqpols,uvsq,umattmp,
     1      vmattmp,wts)
         
         deallocate(umattmp,vmattmp)
      endif


      if(intype.eq.2) then
        call quadasymq_pts(nqorder,nqpols)
        allocate(uvsq(2,nqpols),wts(nqpols))
        call quadasymq(nqorder,tverts(1,1,1),tverts(1,2,1),tverts(1,3,1),
     1         uvsq,wts,nqpols)    
      endif

cc      call prinf('nqpols=*',nqpols,1)

      npmax = nquad*nqpols
      allocate(sigvals(npols,npmax))
      allocate(srcvals(12,npmax),qwts(npmax))

      allocate(sigmatmp(nppols,npmax),fkervals(nd,npmax))
      allocate(rsigtmp(nppols,npmax))

      allocate(uvtmp(2,nqpols))


      do iquad=1,nquad
        call mapuv_quad(tverts(1,1,iquad),nqpols,uvsq,uvtmp)
        istart = (iquad-1)*nqpols
        do i=1,nqpols
          ii = istart+i
          call koornf_pols(uvtmp(1,i),norder,npols,sigvals(1,ii),
     1        rat1,rat2,rsc1)
          call koornf_pols(uvtmp(1,i),nporder,nppols,rsigtmp(1,ii),
     1        rat1p,rat2p,rsc1p)
        enddo
      enddo



      do iquad=1,npatches
c
cc       for the current patch compute all the geometry info
c
        transa = 'N'
        transb = 'N'
        alpha = 1
        beta = 0
        lda = 9
        ldb = npols
        ldc = 12

        call dgemm_guru(transa,transb,9,npmax,npols,alpha,
     1     srccoefs(1,1,iquad),lda,sigvals,ldb,beta,srcvals,ldc)


c
c        compute all the quadrature weights
c

        do i=1,nquad
          istart = (i-1)*nqpols+1
          call get_norms_qwts_quad(nqpols,wts,srcvals(1,istart),
     1        da(i),qwts(istart))
        enddo

        do itarg=itargptr(iquad),itargptr(iquad)+ntargptr(iquad)-1
c
c           extract info of geometry on subquadrangles relevant
c           for this computation
c
          npts = 0
          do iqquad=1,nquad
            jstart = (iqquad-1)*nqpols
            if(iquadrel(iqquad,itarg).eq.1) then
              do i=1,nqpols
                jj = jstart+i
                ii = npts+i
                call fker(nd,srcvals(1,jj),ndtarg,xyztarg(1,itarg),
     1             ndd,dpars,ndz,zpars,ndi,ipars,fkervals(1,ii))
                do idim=1,nd
                  fkervals(idim,ii) = fkervals(idim,ii)*qwts(jj)
                enddo
                do j=1,nppols
                  sigmatmp(j,ii) = rsigtmp(j,jj)
                enddo
              enddo
              npts = npts + nqpols
            endif
          enddo

          transa = 'n'
          transb = 't'
          alpha_c = 1
          beta_c = 0
          call zgemm_guru(transa,transb,nd,nppols,npts,alpha_c,
     1      fkervals,nd,sigmatmp,nppols,beta_c,cintvals(1,1,itarg),
     2      nd)

        enddo
      enddo

      return
      end

c
c
c
c
c
c
c
      subroutine cquadints_adap_vec(eps,intype,
     1     npatches,norder,npols,srccoefs,ndtarg,
     2     ntarg,xyztarg,itargptr,ntargptr,nporder,nppols,nquadmax,
     3     rat1,rat2,rsc1,rat1p,rat2p,rsc1p,
     3     fker,nd,ndd,dpars,ndz,zpars,ndi,ipars,nqorder,cintvals)

c
      implicit none

c
cc     calling sequence variables
c
      real *8 eps
      integer intype,nd
      integer npatches,norder,npols
      integer nporder,nppols
      real *8 srccoefs(9,npols,npatches)
      
      integer ntarg,ndtarg
      real *8 xyztarg(ndtarg,ntarg)
      integer itargptr(npatches)
      integer ntargptr(npatches)
      
      external fker
      integer ndd,ndi,ndz
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer ipars(ndi)

      real *8 rat1(2,0:norder)
      real *8 rat2(3,0:norder,0:norder)
      real *8 rsc1(0:norder,0:norder)

      real *8 rat1p(2,0:norder)
      real *8 rat2p(3,0:norder,0:norder)
      real *8 rsc1p(0:norder,0:norder)

      integer nqorder

      complex *16 cintvals(nd,nppols,ntarg)

c
c       tree variables
c
      integer nlmax,ltree
      real *8, allocatable :: tvs(:,:,:),da(:)
      integer, allocatable :: ichild_start(:)
      real *8, allocatable :: tvs2(:,:,:),da2(:)
      integer, allocatable :: ichild_start2(:)

      integer nquad,nquadmax,nlev,iquad,istart,i,j,k
      integer ier,itarg,jj,jstart,npts
      integer iqquad,ii


      integer npmax

      real *8, allocatable :: uvsq(:,:),wts(:),uvtmp(:,:)
      real *8, allocatable :: umattmp(:,:),vmattmp(:,:)
      integer nqpols
      real *8, allocatable :: sigvals(:,:),sigvalsdens(:,:)
      real *8, allocatable :: srcvals(:,:),qwts(:)
      real *8, allocatable :: sigvals2(:,:),sigvalsdens2(:,:)
      real *8, allocatable :: srcvals2(:,:),qwts2(:)
      integer itmp

      character *1 transa,transb
      real *8 alpha,beta
      integer lda,ldb,ldc
      integer nn1,nn2,nn3,nn4,npmax0,ntmaxuse,ntmaxuse0
      
      
c
c      get the tree
c
      nlmax = 20
      ntmaxuse = nquadmax

      ltree = 2*(nlmax+1) + 6*nquad

c
c       for each quadrangle, we just store three pieces
c       of info
c         quadrangle vertices
c         area of quadrangle
c         ichild_start = index for first child, 
c         
c
      allocate(ichild_start(nquadmax),tvs(2,3,nquadmax))
      allocate(da(nquadmax))

      do i=1,nquadmax
        ichild_start(i) = -1
        da(i) = 0
        do j=1,3
          do k=1,2
            tvs(k,j,i) = 0
          enddo
        enddo
      enddo

      da(1) = 1.0d0
      
      tvs(1,1,1) = 0
      tvs(2,1,1) = 0

      tvs(1,2,1) = 1
      tvs(2,2,1) = 0

      tvs(1,3,1) = 0
      tvs(2,3,1) = 1


c
c       get quadrature nodes and weights on the base quadrangle
c       based on quadrature type
c

      if(intype.eq.1) then
         nqpols = (nqorder+1)*(nqorder+2)/2
         allocate(uvsq(2,nqpols),wts(nqpols))
         allocate(umattmp(nqpols,nqpols),vmattmp(nqpols,nqpols))
      
         call vioreanu_simplex_quad(nqorder,nqpols,uvsq,umattmp,
     1      vmattmp,wts)
         
         deallocate(umattmp,vmattmp)
      endif

      if(intype.eq.2) then
        call quadasymq_pts(nqorder,nqpols)
        allocate(uvsq(2,nqpols),wts(nqpols))
        call quadasymq(nqorder,tvs(1,1,1),tvs(1,2,1),tvs(1,3,1),
     1         uvsq,wts,nqpols)    
      endif

      allocate(uvtmp(2,nqpols))
     
      npmax = nquadmax*nqpols
      allocate(sigvals(npols,npmax),sigvalsdens(nppols,npmax))
      allocate(srcvals(12,npmax),qwts(npmax))


c
c      current number of quadrangles in the adaptive structure
c
      nquad = 1
c
c        intialize sigvals for root quadrangle
c


      call mapuv_quad(tvs(1,1,1),nqpols,uvsq,uvtmp)
      do i=1,nqpols
        call koornf_pols(uvtmp(1,i),norder,npols,
     1      sigvals(1,i),rat1,rat2,rsc1)
        call koornf_pols(uvtmp(1,i),nporder,nppols,
     1      sigvalsdens(1,i),rat1p,rat2p,rsc1p)
      enddo


      do iquad=1,npatches
c
cc       for the current patch compute geometry info for base quadrangle 
c
        transa = 'N'
        transb = 'N'
        alpha = 1
        beta = 0
        lda = 9
        ldb = npols
        ldc = 12

        call dgemm_guru(transa,transb,9,nqpols,npols,alpha,
     1     srccoefs(1,1,iquad),lda,sigvals,ldb,beta,srcvals,ldc)
        call get_norms_qwts_quad(nqpols,wts,srcvals,da,qwts)

        do itarg=itargptr(iquad),itargptr(iquad)+ntargptr(iquad)-1

 1111     continue

          call quadadap_vec(eps,nqorder,nqpols,nlmax,ntmaxuse,nquad,
     1          ichild_start,tvs,da,uvsq,wts, 
     1          norder,npols,srccoefs(1,1,iquad),
     2          npmax,srcvals,
     2          qwts,sigvals,nporder,nppols,sigvalsdens,ndtarg,
     3          xyztarg(1,itarg),
     3          rat1,rat2,rsc1,rat1p,rat2p,rsc1p,
     3          fker,nd,ndd,dpars,ndz,zpars,ndi,ipars,
     3          cintvals(1,1,itarg),ier)
           if(ier.eq.4) then
             ntmaxuse0 = ntmaxuse*4
             npmax0 = ntmaxuse0*nqpols
             allocate(sigvals2(npols,npmax))
             allocate(sigvalsdens2(nppols,npmax))
             allocate(srcvals2(12,npmax),qwts2(npmax))
             allocate(ichild_start2(ntmaxuse),tvs2(2,3,ntmaxuse))
             allocate(da2(ntmaxuse))
             nn1 = npols*npmax
             nn2 = nppols*npmax
             nn3 = 12*npmax
             nn4 = nquad*6
             call dcopy(nn1,sigvals,1,sigvals2,1)
             call dcopy(nn2,sigvalsdens,1,sigvalsdens2,1)
             call dcopy(nn3,srcvals,1,srcvals2,1)
             call dcopy(npmax,qwts,1,qwts2,1)
             do ii=1,nquad
               ichild_start2(ii) = ichild_start(ii)
             enddo
             call dcopy(nn4,tvs,1,tvs2,1)
             call dcopy(nquad,da,1,da2,1)


             deallocate(sigvals,sigvalsdens,srcvals,qwts,ichild_start)
             deallocate(tvs,da)


             allocate(sigvals(npols,npmax0))
             allocate(sigvalsdens(nppols,npmax0))
             allocate(srcvals(12,npmax0),qwts(npmax0))
             allocate(ichild_start(ntmaxuse0),tvs(2,3,ntmaxuse0))
             allocate(da(ntmaxuse0))

             do ii=1,ntmaxuse0
               ichild_start(ii) = -1
             enddo

             call dcopy(nn1,sigvals2,1,sigvals,1)
             call dcopy(nn2,sigvalsdens2,1,sigvalsdens,1)
             call dcopy(nn3,srcvals2,1,srcvals,1)
             call dcopy(npmax,qwts2,1,qwts,1)
             do ii=1,nquad
               ichild_start(ii) = ichild_start2(ii)
             enddo
             call dcopy(nn4,tvs2,1,tvs,1)
             call dcopy(nquad,da2,1,da,1)

             npmax = npmax0
             ntmaxuse = ntmaxuse0
c             call prinf('restarting adaptive inetgration with nquad=*',
c     1             ntmaxuse,1)


             deallocate(sigvals2,sigvalsdens2,srcvals2,ichild_start2)
             deallocate(tvs2,da2,qwts2)
             goto 1111
           endif
        enddo

        do i=1,nquad
          ichild_start(i) = -1
        enddo
      enddo

      return
      end

c
c
c
c
c
      subroutine quadadap_vec(eps,m,kpols,nlmax,ntmax,nquad,
     1             ichild_start,tvs,da,uvsq,wts,
     1             norder,npols,srccoefs,npmax,srcvals,
     2             qwts,sigvals,nporder,nppols,sigvalsdens,ndtarg,xt,
     3             rat1,rat2,rsc1,rat1p,rat2p,rsc1p,
     3             fker,nd,ndd,dpars,ndz,zpars,ndi,
     3             ipars,cintall,ier)

c
c       this subroutine adaptively computes the integral
c        of the functions
c   
c        \int_{T} K_{\ell}(x,y(u,v) P_{n,m}(u,v) |J(u,v)| du dv
c        n,m = 1,2\ldots npols 
c
c        \ell = 1,2,\ldots nd
c
c        P_{n,m}(u,v) are the koornwinder polynomials on the 
c        standard simplex (0,0)-(1,0)-(0,1)
c
c         |J(u,v)| = |dy/du \times dy/dv|
c        
c
c        relevant parameters on a quad refined
c        grid along the way
c
c        IN:
c        eps - precision requested
c        m - quadrature order
c        kpols - number of quarature nodes (look up quadasymq 
c                 for getting kpols(m))
c        nlmax - max level of refinement for geometry
c        ntmax - max number of quadrangles
c        nquad - current number of quadrangles in adaptive structure
c        ichild_start(i) - first child of quadrangle i
c        tvs(2,3,ntmax) - vertices of hierarchy of quadrangles
c        da(ntmax) - area of quadrangles
c        uvsq(kpols) - integration nodes on standard quadrangle
c        wts(kpols) - integration weights on standard quadrangle
c        npols - total number of koornwinder polynomials to be integrated
c        norder - order of discretization of the surface
c        npols - (norder+1)*(norder+2)/2: total number of 
c                koornwinder polynomials to be integrated
c        srccoefs(9,npols) - xyz coefficients of koornwinder expansion 
c                             of current surface + derivative info
c 
c        npmax - max number of points = ntmax*kpols
c        srczvals(12,npmax) - geometry info on heirarchy of meshes
c        qwts(npmax) - quadrature weights 
c        sigvals(npols,npmax) - 
c                   koornwinder polynomials computed along the adaptive grid
c        
c        OUT:
c        cintall(nd,npols) - computed integral 
c
c         

      implicit real *8 (a-h,o-z)
      integer, allocatable :: istack(:)
      integer ichild_start(ntmax),nd
      integer nporder,nppols
      real *8 da(ntmax)
      real *8 tvs(2,3,ntmax), uvsq(2,kpols),wts(kpols)
      integer nproclist0, nproclist
      integer idone
      real *8 srccoefs(9,npols)
      real *8 sigvals(npols,npmax)
      real *8 sigvalsdens(nppols,npmax)
      real *8 srcvals(12,*),qwts(npmax)
      complex *16, allocatable :: xkernvals(:,:)
      real *8 xt(ndtarg)
      complex *16 cintall(nd,npols),fval(nd)
      complex *16, allocatable :: cvals(:,:,:)

      integer ndd,ndz,ndi
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer ipars(ndi),idim

      real *8 rat1(2,0:norder),rat2(3,0:norder,0:norder)
      real *8 rsc1(0:norder,0:norder)

      real *8 rat1p(2,0:nporder),rat2p(3,0:nporder,0:nporder)
      real *8 rsc1p(0:nporder,0:nporder)

      external fker

c
c         for historic reasons
c
      ksigpols = nppols
      allocate(istack(2*ntmax))
      allocate(cvals(nd,ksigpols,ntmax))

      nfunev = 0

      do i=1,ksigpols
        do idim=1,nd
          cvals(idim,i,1) = 0
        enddo
      enddo

      allocate(xkernvals(nd,npmax))

c
cc      compute integral at level 0
c
      do i=1,kpols
         call fker(nd,srcvals(1,i),ndtarg,xt,ndd,dpars,ndz,zpars,
     1        ndi,ipars,fval)
         do idim=1,nd
           xkernvals(idim,i) = fval(idim)*qwts(i)
         enddo
         do j=1,ksigpols
           do idim=1,nd
             cvals(idim,j,1) = cvals(idim,j,1)+
     1           xkernvals(idim,i)*sigvalsdens(j,i)
           enddo
         enddo
      enddo

      nfunev = nfunev + kpols

      
      do i=1,ksigpols
        do idim=1,nd
          cintall(idim,i) = cvals(idim,i,1)
        enddo
      enddo

      nproclist0 = 1
      istack(1) = 1


      call quadadap_main_vec(eps,kpols,nlmax,ntmax,nquad,ichild_start,
     1      tvs,da,uvsq,wts,norder,npols,srccoefs,
     2      npmax,srcvals,qwts,sigvals,nporder,nppols,sigvalsdens,
     3      ndtarg,xt,rat1,rat2,rsc1,rat1p,rat2p,rsc1p,
     3      fker,nd,ndd,dpars,ndz,
     3      zpars,ndi,ipars,cvals,istack,nproclist0,
     4      xkernvals,cintall,ier)
      
      return
      end
c
c
c
c
c
       
      subroutine quadadap_main_vec(eps,kpols,nlmax,ntmax,nquad,
     1      ichild_start,tvs,da,uvsq,wts,norder,npols,srccoefs,
     2      npmax,srcvals,qwts,sigvals,nporder,nppols,sigvalsdens,
     3      ndtarg,xt,rat1,rat2,rsc1,rat1p,rat2p,rsc1p,
     3      fker,nd,ndd,dpars,ndz,
     3      zpars,ndi,ipars,cvals,istack,nproclist0,xkernvals,
     4      cintall,ier)
      

      implicit real *8 (a-h,o-z)
      integer istack(*),nproclist0,nd
      integer ichild_start(ntmax)
      integer nporder,nppols
      real *8 da(ntmax)
      real *8 tvs(2,3,ntmax), uvsq(2,kpols),wts(kpols)
      integer  nproclist
      integer idone
      real *8 srccoefs(9,npols)
      real *8 sigvals(npols,npmax)
      real *8 sigvalsdens(nppols,npmax)
      real *8 qwts(npmax)
      complex *16 xkernvals(nd,npmax)
      real *8 xt(ndtarg)
      real *8 srcvals(12,*)
      complex *16 cintall(nd,nppols),fval(nd),ctmp(nd,nppols)
      complex *16 cvals(nd,nppols,ntmax)

      integer ndd,ndz,ndi
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer ipars(ndi)

      real *8 rat1(2,0:norder),rat2(3,0:norder,0:norder)
      real *8 rsc1(0:norder,0:norder)

      real *8 rat1p(2,0:nporder),rat2p(3,0:nporder,0:nporder)
      real *8 rsc1p(0:nporder,0:nporder)

      real *8, allocatable :: uvtmp(:,:)
      character *1 transa,transb
      integer lda,ldb,ldc
      external fker
      
      allocate(uvtmp(2,kpols))


c
c         for historic reasons
c
      ksigpols = nppols
      kfine = 4*kpols
      ier = 0


      do ilev=0,nlmax
        idone = 1
        nproclist = 0
        
        do iproc = 1,nproclist0
          iquad = istack(iproc)

c
c           check to see if quadrangle already has 
c           children, if not, set children
c           and compute necessary info
c
          if(ichild_start(iquad).eq.-1) then

c
c            current quadrangle doesn't have children,
c            compute necessary info
c

            if(nquad+4.gt.ntmax) then
c               print *, "Too many quadrangles in cquadadap"
c               print *, "Exiting without computing anything"
               ier = 4
               return
            endif
            
            ichild_start(iquad) = nquad+1
            call getquadchildren(tvs(1,1,iquad),tvs(1,1,nquad+1),
     1             tvs(1,1,nquad+2),tvs(1,1,nquad+3),tvs(1,1,nquad+4))

            rr = 0.25d0*da(iquad)
            do j=nquad+1,nquad+4
              da(j) = rr
              call mapuv_quad(tvs(1,1,j),kpols,uvsq,uvtmp)
              istart = (j-1)*kpols+1
              do i=1,kpols
                 ii = istart+i-1
                call koornf_pols(uvtmp(1,i),norder,npols,sigvals(1,ii),
     1             rat1,rat2,rsc1)
                call koornf_pols(uvtmp(1,i),nporder,nppols,
     1             sigvalsdens(1,ii),rat1p,rat2p,rsc1p)
              enddo
              
              transa = 'N'
              transb = 'N'
              alpha = 1
              beta = 0
              lda = 9
              ldb = npols
              ldc = 12

              call dgemm_guru(transa,transb,9,kpols,npols,alpha,
     1           srccoefs,lda,sigvals(1,istart),ldb,beta,
     2           srcvals(1,istart),ldc)
              call get_norms_qwts_quad(kpols,wts,srcvals(1,istart),
     1           rr,qwts(istart))

            enddo
            nquad = nquad+4
          endif
        
c
cc           compute xkernvals
c
          iquadc1 = ichild_start(iquad)
          istart = (iquadc1-1)*kpols
          do j=1,kfine
            jj=j+istart
            call fker(nd,srcvals(1,jj),ndtarg,xt,ndd,dpars,
     1         ndz,zpars,ndi,ipars,fval)
            do idim=1,nd
              xkernvals(idim,jj) = fval(idim)*qwts(jj)
            enddo
          enddo
cc          call prin2('xkernvals=*',xkernvals(istart+1),kfine)

          nfunev = nfunev + kfine

c
cc         subtract conquadbution of current 
c          quadrangle
          do isig=1,ksigpols
            do idim=1,nd
              cintall(idim,isig) = cintall(idim,isig)-
     1            cvals(idim,isig,iquad)
              ctmp(idim,isig) = 0
            enddo
          enddo

c
cc        add in conquadbutions of children quadrangles
c
          do iquadc=iquadc1,iquadc1+3
            do isig=1,ksigpols
              do idim=1,nd
                cvals(idim,isig,iquadc) = 0
              enddo
            enddo

            istart = (iquadc-1)*kpols
            do k=1,kpols
              ii = istart+k
              do isig=1,ksigpols
                do idim=1,nd
                  cvals(idim,isig,iquadc) = cvals(idim,isig,iquadc)+
     1               xkernvals(idim,ii)*sigvalsdens(isig,ii)
                enddo
              enddo
            enddo
              
            do isig=1,ksigpols
              do idim=1,nd
                cintall(idim,isig) = cintall(idim,isig) + 
     1              cvals(idim,isig,iquadc)
                ctmp(idim,isig) = ctmp(idim,isig)+
     1              cvals(idim,isig,iquadc)
              enddo
            enddo
          enddo

c
cc        compare integral of children to integral of parent
c         to determine if the children need to be refined further
c
          errmax = 0
          do isig=1,ksigpols
            do idim=1,nd
              if(abs(ctmp(idim,isig)-cvals(idim,isig,iquad)).gt.errmax) 
     1          errmax = abs(ctmp(idim,isig)-cvals(idim,isig,iquad))
            enddo
          enddo

          if(errmax.gt.eps) then
            idone = 0
            
            do j=1,4
              istack(nproclist0+nproclist+j) = iquadc1+j-1
            enddo
            nproclist = nproclist+4
          endif
c
cc        end of looping over all quadrangles at current stage
        enddo
cc         if idone is still 1, that means that no more refinement
c          is needed
         if(idone.eq.1) goto 1111
        do i=1,nproclist
          istack(i) = istack(nproclist0+i)
        enddo
        nproclist0 = nproclist
      enddo
 1111 continue


      return
      end
c
c
c
c
c
c
c
      subroutine cquadints_comb_vec(eps,intype,
     1     npatches,norder,npols,srccoefs,ndtarg,
     2     ntarg,xyztarg,ifp,xyzproxy,
     3     itargptr,ntargptr,nporder,nppols,nquadmax,rat1,rat2,rsc1,
     3     rat1p,rat2p,rsc1p,
     3     fker,nd,ndd,dpars,ndz,zpars,ndi,ipars,
     4     nqorder,rfac,cintvals)
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
c        this code is the same as cquadints on all counts except
c        for choosing when to subdivide a quadrangle
c        
c        here a quadrangle is subdivided when xyzproxy associated
c        with the quadrangle is separated from the centroid
c        of the quadrangle by at least rfac*r_{t} where r_{t}
c        is the radius of the smallest sphere  
c
c        input arguments:
c        eps:     requested precision (Currently not being used
c                 in any way, later, nqorder, and rfac will
c                 be set internally based on eps)
c
c        intype:   quadrature node type
c                   intype = 1, rokhlin vioreanu nodes
c                   intype = 2, xiao gimbutas nodes
c   
c        npatches: number of patches
c        norder: order of discretization nodes on the patches
c        npols = (norder+1)*(norder+2)/2: number of discretization 
c                   nodes on each patch
c        srccoefs(3,npols,npatches): coefficients
c                 of koornwinder expansion of xyz coordinates
c                 and  of dxyz/du and dxyz/dv
c
c        ndtarg - dimension of each target point (
c          ndtarg must at least be 3, and correspond to location
c          of target points)
c        ntarg - total number of target points
c        xyztarg(ndtarg,ntarg) - location of target points 
c                       (currently may have repeats, 
c                        being slightly wasteful in
c                        memory, but see comment (1))
c
c        ifp - flag for using proxy target locations
c               proxy target locations will be used if ifp = 1
c               else original target locations will be used
c        xyzproxy(3,*) - proxy locations for measuring distances
c
c        itargptr(npatches) - xyztargs(:,itargptr(i):itargptr(i)+ntargptr(i)-1)
c                       are the relevant set of targets for patch i
c        ntargptr(npatches) - number of relevant targets for patch i
c        nquadmax - max number of quadrangles to be used on base quadrangle
c        fker - function handle for evaluating the kernel k
c 
c               expected calling sequence
c               fker(x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
c               
c               the output is assumed to be complex for the time
c               being
c
c         nd - number of kernels to be integrated
c         dpars(*) - real parameters for the fker routine
c         zpars(*) - complex parameters for the fker routine
c         ipars(*) - integer parameters for the fker routine
c         nqorder - order of quadrature nodes on each subquadrangle
c                   to be used
c         rfac - distance criterion for deciding whether a quadrangle
c                is in the far-field or not
c
c         output:
c
c         cintvals(nd,npols,ntarg) - integrals at all targets
c                                  for all koornwinder
c                                  polynomials
c
      implicit none

c
cc     calling sequence variables
c
      real *8 eps
      integer intype,ifp,nd
      integer npatches,norder,npols
      integer nporder,nppols
      real *8 srccoefs(9,npols,npatches)
      
      integer ntarg,ndtarg
      real *8 xyztarg(ndtarg,ntarg)
      real *8 xyzproxy(3,*)
      integer itargptr(npatches)
      integer ntargptr(npatches)
      
      external fker
      integer ndd,ndz,ndi
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer ipars(ndi)

      integer nqorder
      real *8 rfac

      real *8 rat1(2,0:norder),rat2(3,0:norder,0:norder)
      real *8 rsc1(0:norder,0:norder)

      real *8 rat1p(2,0:nporder),rat2p(3,0:nporder,0:nporder)
      real *8 rsc1p(0:nporder,0:nporder)

      complex *16 cintvals(nd,nppols,ntarg)
      

c
cc      temporary variables
c
      real *8, allocatable :: quadcm(:,:,:)
      real *8, allocatable :: quadrad(:,:)
      integer, allocatable :: iquadreltmp(:,:)
      integer, allocatable :: iquadrel(:,:)
      integer, allocatable :: iquadrelall(:)
      real *8, allocatable :: tverts(:,:,:)
      integer, allocatable :: ichild_start(:)
      integer, allocatable :: ichild_start0(:)

      real *8, allocatable :: xyztargtmp(:,:)

      integer nquad,nquadmax,nlev,iquad,istart,i,j
      integer ier,itarg,jj,jstart,nlmax,npts
      integer iqquad,ii

      integer npmax
      integer idim

      real *8, allocatable :: uvsq(:,:),wts(:),uvtmp(:,:)
      real *8, allocatable :: umattmp(:,:),vmattmp(:,:)
      real *8, allocatable :: da(:)
      integer nqpols
      real *8, allocatable :: sigvals(:,:)
      real *8, allocatable :: sigvalsdens(:,:)
      real *8, allocatable :: srcvals(:,:)
      real *8, allocatable :: qwts(:)
      complex *16, allocatable :: sigmatmp(:,:)
      real *8 xyztmp(3)
      integer itmp

      complex *16, allocatable :: xkernvals(:,:)
      integer, allocatable :: istack(:)

      complex *16, allocatable :: cvals(:,:,:)
      
      integer nproclist0,nproclist,nquad0,npts0
      complex *16 zz

      
      complex *16 ima

      character *1 transa,transb
      
      real *8 alpha,beta
      integer lda,ldb,ldc

      data ima/(0.0d0,1.0d0)/

      allocate(cvals(nd,nppols,nquadmax))
      allocate(istack(2*nquadmax))





cc      max number of levels
c
      nlmax = 20
      allocate(quadcm(3,npatches,nquadmax),
     1   quadrad(npatches,nquadmax),tverts(2,3,nquadmax))
      allocate(iquadreltmp(ntarg,nquadmax))
      allocate(da(nquadmax),ichild_start(nquadmax))
      allocate(ichild_start0(nquadmax))

      do i=1,nquadmax
        do j=1,ntarg
          iquadreltmp(j,i) = 0
        enddo
        
      enddo

      
      nquad = 0
      nlev = 0
      ier = 0

      do i=1,nquadmax
        ichild_start(i) = -1
      enddo

      

      if(ifp.eq.1) then
        call getquadtree(npatches,norder,npols,srccoefs,ntarg,xyzproxy,
     1     itargptr,ntargptr,nquadmax,nlmax,rfac,nquad,nlev,ichild_start,
     2     da,quadcm,quadrad,tverts,iquadreltmp,ier)
      else
        allocate(xyztargtmp(3,ntarg))
        do i=1,ntarg
          do j=1,3
            xyztargtmp(j,i) = xyztarg(j,i)
          enddo
        enddo
        call getquadtree(npatches,norder,npols,srccoefs,ntarg,
     1     xyztargtmp,
     1     itargptr,ntargptr,nquadmax,nlmax,rfac,nquad,nlev,ichild_start,
     2     da,quadcm,quadrad,tverts,iquadreltmp,ier)
        deallocate(xyztargtmp)
      endif

      nquad0 = nquad

      do i=1,nquad0
        ichild_start0(i) = ichild_start(i)
      enddo


      

      if(ier.ne.0) then
        call prinf('Could not allocate quadrangle tree*',i,0)
        call prinf('Exiting without computing anything*',i,0)
        call prinf('Press 0 to continue*',i,0)
        call prinf('Press 1 to exit routine and continue*',i,0)
        call prinf('Press 2 to stop*',i,0)
        read *, ier
        if(ier.eq.2) stop
        if(ier.eq.1) return
      endif

      allocate(iquadrel(nquad0,ntarg),iquadrelall(nquad0))

c
c       transpose the iquadrel array
c  

      do iquad=1,nquad0
        do j=1,ntarg
          iquadrel(iquad,j) = iquadreltmp(j,iquad)
        enddo
      enddo


c
c       compute all the geometry info and koornwinder polynomials
c       at relevant quadrangle given by iquadrelall
c  

c
c       get quadrature nodes and weights on the base quadrangle
c       based on quadrature type
c

      if(intype.eq.1) then
         nqpols = (nqorder+1)*(nqorder+2)/2
         allocate(uvsq(2,nqpols),wts(nqpols))
         allocate(umattmp(nqpols,nqpols),vmattmp(nqpols,nqpols))
      
         call vioreanu_simplex_quad(nqorder,nqpols,uvsq,umattmp,
     1      vmattmp,wts)
         
         deallocate(umattmp,vmattmp)
      endif


      if(intype.eq.2) then
        call quadasymq_pts(nqorder,nqpols)
        allocate(uvsq(2,nqpols),wts(nqpols))
        call quadasymq(nqorder,tverts(1,1,1),tverts(1,2,1),tverts(1,3,1),
     1         uvsq,wts,nqpols)    
      endif

cc      call prinf('nqpols=*',nqpols,1)

      npts0 = nquad0*nqpols
      npmax = nquadmax*nqpols
      allocate(sigvals(npols,npmax))
      allocate(sigvalsdens(nppols,npmax))
      allocate(srcvals(12,npmax),qwts(npmax))

      allocate(xkernvals(nd,npmax))

      allocate(uvtmp(2,nqpols))


      do iquad=1,nquad0
        call mapuv_quad(tverts(1,1,iquad),nqpols,uvsq,uvtmp)
        istart = (iquad-1)*nqpols
        do i=1,nqpols
          ii = istart+i
          call koornf_pols(uvtmp(1,i),norder,npols,sigvals(1,ii),
     1           rat1,rat2,rsc1)
          call koornf_pols(uvtmp(1,i),nporder,nppols,sigvalsdens(1,ii),
     1           rat1p,rat2p,rsc1p)
        enddo
      enddo


      do iquad=1,npatches
        nquad = nquad0
c
cc       for the current patch compute all the geometry info
c

        transa = 'N'
        transb = 'N'
        alpha = 1.0d0
        beta = 0.0d0
        lda = 9
        ldb = npols
        ldc = 12

        call dgemm_guru(transa,transb,9,npts0,npols,alpha,
     1     srccoefs(1,1,iquad),lda,sigvals,ldb,beta,srcvals,ldc)



c
c        compute all the quadrature weights
c

        do i=1,nquad
          istart = (i-1)*nqpols+1
          call get_norms_qwts_quad(nqpols,wts,srcvals(1,istart),
     1        da(i),qwts(istart))
        enddo

        do itarg=itargptr(iquad),itargptr(iquad)+ntargptr(iquad)-1
          do i=1,nppols
            do idim=1,nd
              cintvals(idim,i,itarg) = 0
            enddo
          enddo
c
c           extract info of geometry on subquadrangles relevant
c           for this computation
c
          nproclist0 = 0
          do iqquad=1,nquad0
            jstart = (iqquad-1)*nqpols
            if(iquadrel(iqquad,itarg).eq.1) then
              nproclist0 = nproclist0 + 1
              istack(nproclist0) = iqquad
c
c
c              compute the integral due to iquad
c

              do i=1,nppols
                do idim=1,nd
                  cvals(idim,i,iqquad) = 0
                enddo
              enddo
              do i=1,nqpols
                jj = jstart+i
                call fker(nd,srcvals(1,jj),ndtarg,xyztarg(1,itarg),
     1              ndd,dpars,ndz,zpars,ndi,ipars,xkernvals(1,jj))
                do idim=1,nd
                  xkernvals(idim,jj) = xkernvals(idim,jj)*qwts(jj)
                enddo
                
                do j=1,nppols
                  do idim=1,nd
                    zz = xkernvals(idim,jj)*sigvalsdens(j,jj)
                    cvals(idim,j,iqquad) = cvals(idim,j,iqquad)+zz
                    cintvals(idim,j,itarg)=cintvals(idim,j,itarg)+zz
                  enddo
                enddo
              enddo
            endif
          enddo
c
c           done with initial computation of 
c


          call quadadap_main_vec(eps,nqpols,nlmax,nquadmax,nquad,
     1      ichild_start,
     1      tverts,da,uvsq,wts,norder,npols,srccoefs(1,1,iquad),
     2      npmax,srcvals,qwts,sigvals,nporder,nppols,sigvalsdens,
     3      ndtarg,xyztarg(1,itarg),
     3      rat1,rat2,rsc1,rat1p,rat2p,rsc1p,
     3      fker,nd,ndd,dpars,ndz,zpars,ndi,ipars,cvals,istack,
     3      nproclist0,
     4      xkernvals,cintvals(1,1,itarg),ier)

      
        enddo

        do i=1,nquad0
          ichild_start(i) = ichild_start0(i)
        enddo
        do i=nquad0+1,nquad
          ichild_start(i) = -1
        enddo
      enddo

      return
      end
c
c
c
c
c
c
