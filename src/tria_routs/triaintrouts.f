c  This file contains routines for computing integrals of the form
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
c
c  Routines in this file
c  ----------------------       
c  * ctriaints_dist: refine until targets/proxy separated
c                    by triangle by rfac*r_{t}
c                    where r_{t} is the radius of enclosing
c                    sphere for triangle
c 
c  * ctriaints_adap: completely adaptive integration
c                    for the triangle
c           
c  * ctriaints_comb: combination of adaptive integration
c                    and distance based refinement.
c                    Use distance based refinement to compute
c                    the integrals on some coarse hierarchy
c                    of meshes and then use 
c                    adaptive integration from that
c                    point on
c
c  * ctriaints_wnodes: compute integral using prescribed nodes and 
c                      weights
c                       
c
c  Vectorized routines
c  --------------------
c         
c * ctriaints_dist_vec: refine until targets/proxy separated
c                       by triangle by rfac*r_{t}
c                       where r_{t} is the radius of enclosing
c                       sphere for triangle
c                       vectorized verison of ctriaints_dist
c
c 
c * ctriaints_adap_vec: completely adaptive integration
c                       for the triangle
c                       vectorized verion of ctriaints_adap
c
c           
c * ctriaints_comb_vec: combination of adaptive integration
c                       and distance based refinement.
c                       Use distance based refinement to compute
c                       the integrals on some coarse hierarchy
c                       of meshes and then use 
c                       adaptive integration from that
c                       point on
c                       vectorized version of ctriaints_comb
c
* ctriaints_wnodes_vec: compute integral using prescribed nodes and 
c                       weights. 
c                       Vectorized version of ctriaints_wnodes

c
c
      subroutine ctriaints_wnodes(npatches, norder, npols,
     1   isd, ndsc, srccoefs, ndtarg, ntarg, xyztarg, itargptr, 
     2   ntargptr, nporder, nppols, fker, ndd, dpars, ndz, zpars, 
     3   ndi, ipars, nqpts, qnodes, wts, cintvals)
c
c  Compute the integrals in (1) defined at the top of the file
c  using a prescribed set of nodes and weights.
c
c  Input arguments:
c    - npatches: integer *8
c        number of patches
c    - norder: integer *8
c        order of discretization on patches
c    - npols: integer *8
c        number of points/basis functions on the 
c        patch = (norder+1)*(norder+2)/2
c    - isd: integer *8
c        flag for computing and passing second derivative
c        information in srcvals.
c        if isd = 0, leading dim of srcvals = 12, with 
c        xyz, dxyz/du, dxyz/dv, and normals passed
c        else, leading dim of srcvals = 30, with
c        xyz, dxyz/du, dxyz/dv, normals, d2xyz/du2, d2xyz/duv,
c        d2xyz/dv2, det(g), (11,12,22) entries of I, (11,12,22)
c        entries of II, kap1, kap2 where kap1, kap2
c        are principal curvatures, det(g) is the determinant
c        of the metric tensor
c    - ndsc: integer *8
c        leading dimension of srccoefs array, should be 9
c        if isd = 0, and 18 otherwise
c    - srccoefs: real *8 (ndsc, npols, npatches)
c        Koornwinder expansion coefficients of xyz, d/du (xyz), 
c        d/dv (xyz), (d2(xyz)/du2, d2(xyz)/duv, d2(xyz)/dv2)
c    - ndtarg: integer *8
c        leading dimension of target information
c    - ntarg: integer *8
c        number of targets
c    - xyztarg: real *8 (ndtarg, ntarg)
c        target information
c    - itargptr: integer *8(npatches)
c        Pointer array for determining which targets are relevant
c        for patch i. 
c        xyztargs(:,itargptr(i):itargptr(i)+ntargptr(i)-1)
c        are the relevant set of targets for patch i
c    - ntargptr: integer *8(npatches)
c        number of targets relevant for patch i
c    - nporder: integer *8
c        order of Koornwinder polynomials to be integrated
c    - nppols: integer *8
c        number of polynomials corresponding to 
c        nporder = (nporder+1)*(nporder+2)/2
c    - fker: function handle
c        function handle for evaluating the kernel k
c        * expected calling sequence
c            fker(x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
c               
c        In this routine the output is expected to be a complex
c        scalar
c    - ndd: integer *8
c        number of real *8/double precision parameters
c    - dpars: real *8 (ndd)
c         real *8/ double precision paramters
c    - ndz: integer *8
c        number of complex *16 parameters
c    - zpars: complex *16(ndz)
c        complex *16 parameters
c    - ndi: integer *8
c        number of integer parameters
c    - ipars: integer *8(ndi)
c        integer parameters
c    - nqpts: integer *8
c        number of quadrature points
c    - qnodes: real *8 (2,nqpts)
c        u,v location of the quadrature nodes
c    - wts: real *8 nqpts
c        quadrature weights
c 
c  Output arguments:
c    - cintvals: complex *16 (nppols, ntarg)
c        the computed integrals
c 

      implicit none
      integer *8, intent(in) :: npatches, norder, npols, ndtarg
      integer *8, intent(in) :: nporder, nppols
      integer *8, intent(in) :: ndsc
      real *8, intent(in) :: srccoefs(ndsc,npols,npatches)
      real *8, intent(in) :: xyztarg(ndtarg,ntarg)
      integer *8, intent(in) :: ntarg
      integer *8, intent(in) :: itargptr(npatches), ntargptr(npatches)
      integer *8, intent(in) :: isd
      integer *8, intent(in) :: ndd, ndz, ndi
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      integer *8, intent(in) :: ipars(ndi),nqpts
      real *8, intent(in) :: qnodes(2,nqpts),wts(nqpts)

      complex *16, intent(out) :: cintvals(nppols,ntarg)

      real *8, allocatable :: rat1(:,:), rat2(:,:,:), rsc1(:,:)
      real *8, allocatable :: rat1p(:,:), rat2p(:,:,:), rsc1p(:,:)
      real *8, allocatable :: srcvals(:,:),qwts(:)
      real *8, allocatable :: rsigvals(:,:),rsigtmp(:)
      complex *16, allocatable :: sigvals(:,:)
      complex *16, allocatable :: xkernvals(:,:)

      integer *8 i,ipatch,j,lda,ldb,itarg,ldc,ntarg0,ii
      integer *8 nds

      complex *16 fval
      real *8 da,ra

      character *1 transa,transb
      real *8 alpha,beta
      complex *16 alpha_c,beta_c

      external fker

c
c       initialize koornwinder polynomials
c

      
      allocate(rat1(2,0:norder),rat2(3,0:norder,0:norder))
      allocate(rsc1(0:norder,0:norder))
      call koornf_init(norder,rat1,rat2,rsc1) 

      allocate(rat1p(2,0:nporder),rat2p(3,0:nporder,0:nporder))
      allocate(rsc1p(0:nporder,0:nporder))
      call koornf_init(nporder,rat1p,rat2p,rsc1p) 

      allocate(rsigtmp(nppols))
      allocate(sigvals(nppols,nqpts),rsigvals(npols,nqpts))
      
      nds = 12
      if(isd.ne.0) nds = 30
      allocate(srcvals(nds,nqpts),qwts(nqpts))

      do i=1,nqpts
        call koornf_pols(qnodes(1,i),norder,npols,rsigvals(1,i),
     1        rat1,rat2,rsc1)
        call koornf_pols(qnodes(1,i),nporder,nppols,rsigtmp,rat1p,
     2        rat2p,rsc1p)
        do j=1,nppols
          sigvals(j,i) = rsigtmp(j)
        enddo
      enddo




      do ipatch=1,npatches
c
c        compute srcvals
c
        transa = 'N'
        transb = 'N'
        alpha = 1
        beta = 0
        lda = ndsc
        ldb = npols
        ldc = nds

        da = 1

        call dgemm_guru(transa,transb,ndsc,nqpts,npols,alpha,
     1     srccoefs(1,1,ipatch),lda,rsigvals,ldb,beta,srcvals,ldc)

        call get_srcvals_auxinfo_tri(nqpts, wts, isd, nds, srcvals, 
     1    da, qwts)
c
c
c          compute the kernel values for all the targets
c
        
        ntarg0 = ntargptr(ipatch)
        allocate(xkernvals(ntarg0, nqpts))
        do j=1,nqpts
          do itarg=itargptr(ipatch),itargptr(ipatch)+ntarg0-1
            ii = itarg - itargptr(ipatch)+1
            call fker(srcvals(1,j),ndtarg,xyztarg(1,itarg),ndd,dpars,
     1         ndz,zpars,ndi,ipars,fval)
            xkernvals(ii,j) = fval*qwts(j)
          enddo
        enddo


        transa = 'n'
        transb = 't'
        alpha_c = 1
        beta_c = 0

      call zgemm_guru(transa, transb, nppols, ntarg0, nqpts,
     1   alpha_c, sigvals, nppols, xkernvals, ntarg0, beta_c, 
     2   cintvals(1,itargptr(ipatch)), nppols)
      
        deallocate(xkernvals)
c

      enddo

      

      return
      end


      subroutine ctriaints_dist(eps, intype, npatches, norder, npols,
     1     isd, ndsc, srccoefs, ndtarg, ntarg, xyztarg, ifp, xyzproxy,
     2     itargptr, ntargptr, nporder, nppols, ntrimax, rat1, rat2, 
     3     rsc1, rat1p, rat2p, rsc1p, fker, ndd, dpars, ndz, zpars,
     4     ndi, ipars, nqorder, rfac, cintvals, ier)
c
c  Compute the integrals in (1) defined at the top of the file
c  using a distance based criterion to subdivide a triangle.
c  A triangle is subdivided when xyzproxy associated
c  with the triangle is separated from the centroid
c  of the triangle by at least rfac*r_{t} where r_{t}
c  is the radius of the smallest sphere  
c
c
c  Input arguments:
c    - eps: real *8
c        requested tolerance. Currently unsued, nqorder and rfac,
c        should be set based on eps.
c    - intype: integer *8
c        node type
c        * intype = 1, rokhlin vioreanu nodes
c        * intype = 2, xiao gimbutas nodes
c    - npatches: integer *8
c        number of patches
c    - norder: integer *8
c        order of discretization on patches
c    - npols: integer *8
c        number of points/basis functions on the 
c        patch = (norder+1)*(norder+2)/2
c    - isd: integer *8
c        flag for computing and passing second derivative
c        information in srcvals.
c        if isd = 0, leading dim of srcvals = 12, with 
c        xyz, dxyz/du, dxyz/dv, and normals passed
c        else, leading dim of srcvals = 30, with
c        xyz, dxyz/du, dxyz/dv, normals, d2xyz/du2, d2xyz/duv,
c        d2xyz/dv2, det(g), (11,12,22) entries of I, (11,12,22)
c        entries of II, kap1, kap2 where kap1, kap2
c        are principal curvatures, det(g) is the determinant
c        of the metric tensor
c    - ndsc: integer *8
c        leading dimension of srccoefs array, should be 9
c        if isd = 0, and 18 otherwise
c    - srccoefs: real *8 (ndsc, npols, npatches)
c        Koornwinder expansion coefficients of xyz, d/du (xyz), 
c        d/dv (xyz), (d2(xyz)/du2, d2(xyz)/duv, d2(xyz)/dv2)
c    - ndtarg: integer *8
c        leading dimension of target information
c    - ntarg: integer *8
c        number of targets
c    - xyztarg: real *8 (ndtarg, ntarg)
c        target information
c    - ifp: integer *8
c        flag for using proxy points for refining based 
c        on distance criterion
c    - xyzproxy: real *8 (3,*)
c        Proxy target points for refining triangles based
c        distance criterion, distance to proxy
c        point will be used instead of distance to target
c        should be of size (3,ntarg), if ifp=1
c    - itargptr: integer *8(npatches)
c        Pointer array for determining which targets are relevant
c        for patch i. 
c        xyztargs(:,itargptr(i):itargptr(i)+ntargptr(i)-1)
c        are the relevant set of targets for patch i
c    - ntargptr: integer *8(npatches)
c        number of targets relevant for patch i
c    - nporder: integer *8
c        order of Koornwinder polynomials to be integrated
c    - nppols: integer *8
c        number of polynomials corresponding to 
c        nporder = (nporder+1)*(nporder+2)/2
c    - ntrimax: integer *8
c        maximum number of triangles allowed in heirarchy
c        of triangles. Routine will return without
c        computing anything and an error code, if ntrimax
c        is too small. Recommended value 3000.
c    - rat1: real *8 (2, 0:norder)
c        recurrence coefficients for Legendre polynomials
c        with discretization order = norder
c    - rat2: real *8 (3, 0:norder, 0:norder)
c        recurrence coefficients for Jacobi polynomials
c        with discretization order = norder
c    - rsc1: real *8 (0:norder, 0:norder)
c        scaling parameters for the Koornwinder polynomials
c        with discretization order = norder
c    - rat1p: real *8 (2, 0:nporder)
c        recurrence coefficients for Legendre polynomials
c        with discretization order = nporder
c    - rat2p: real *8 (3, 0:nporder, 0:nporder)
c        recurrence coefficients for Jacobi polynomials
c        with discretization order = nporder
c    - rsc1p: real *8 (0:nporder, 0:nporder)
c        scaling parameters for the Koornwinder polynomials
c        with discretization order = nporder
c    - fker: function handle
c        function handle for evaluating the kernel k
c        * expected calling sequence
c            fker(x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
c               
c        In this routine the output is expected to be a complex
c        scalar
c    - ndd: integer *8
c        number of real *8/double precision parameters
c    - dpars: real *8 (ndd)
c         real *8/ double precision paramters
c    - ndz: integer *8
c        number of complex *16 parameters
c    - zpars: complex *16(ndz)
c        complex *16 parameters
c    - ndi: integer *8
c        number of integer parameters
c    - ipars: integer *8(ndi)
c        integer parameters
c    - nqorder: integer *8
c        order of quadrature nodes to be used on each triangle
c    - rfac: real *8
c        scaling factor for determining when a triangle 
c        needs refinement
c 
c  Output arguments:
c    - cintvals: complex *16 (nppols, ntarg)
c        the computed integrals
c    - ier: integer *8
c        error code, ier = 0 means successful execution
c        ier = 2, means that too few triangles, code should be run
c        with a larger value of ntrimax
c      
c
      implicit none

c
cc     calling sequence variables
c
      real *8 eps
      integer *8 intype,ifp
      integer *8 npatches,norder,npols
      integer *8 nporder,nppols
      integer *8 isd, ndsc
      real *8 srccoefs(ndsc,npols,npatches)
      
      integer *8 ntarg,ndtarg
      real *8 xyztarg(ndtarg,ntarg)
      real *8 xyzproxy(3,*)
      integer *8 itargptr(npatches)
      integer *8 ntargptr(npatches)
      
      external fker
      integer *8 ndd,ndz,ndi
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer *8 ipars(ndi)

      integer *8 nqorder
      real *8 rfac
      real *8 rat1(2,0:norder),rat2(3,0:norder,0:norder)
      real *8 rsc1(0:norder,0:norder)

      real *8 rat1p(2,0:nporder),rat2p(3,0:nporder,0:nporder)
      real *8 rsc1p(0:nporder,0:nporder)

      complex *16 cintvals(nppols,ntarg)
      

c
cc      temporary variables
c
      real *8, allocatable :: tricm(:,:,:)
      real *8, allocatable :: trirad(:,:)
      integer *8, allocatable :: itrireltmp(:,:)
      integer *8, allocatable :: itrirel(:,:)
      integer *8, allocatable :: itrirelall(:)
      real *8, allocatable :: tverts(:,:,:)
      integer *8, allocatable :: ichild_start(:)

      integer *8 ntri,ntrimax,nlev,itri,istart,i,j
      integer *8 ier,itarg,jj,jstart,nlmax,npts
      integer *8 iqtri,ii

      integer *8 npmax

      real *8, allocatable :: uvsq(:,:),wts(:),uvtmp(:,:)
      real *8, allocatable :: umattmp(:,:),vmattmp(:,:)
      real *8, allocatable :: da(:)
      integer *8 nqpols
      real *8, allocatable :: sigvals(:,:)
      real *8, allocatable :: srcvals(:,:),qwts(:)
      real *8, allocatable :: xyztargtmp(:,:)
      complex *16, allocatable :: fkervals(:),sigmatmp(:,:)
      real *8, allocatable :: rsigtmp(:,:)
      real *8 xyztmp(3)
      integer *8 itmp,nds
      complex *16 ima
      complex *16 alpha_c, beta_c

      character *1 transa,transb
      double precision :: alpha,beta
      integer *8 lda,ldb,ldc
      integer *8 int8_1
      

      data ima/(0.0d0,1.0d0)/





      int8_1 = 1
cc      max number of levels
c
      nlmax = 20
      allocate(tricm(3,npatches,ntrimax),
     1   trirad(npatches,ntrimax),tverts(2,3,ntrimax))
      allocate(itrireltmp(ntarg,ntrimax))
      allocate(da(ntrimax),ichild_start(ntrimax))

      do i=1,ntrimax
        do j=1,ntarg
          itrireltmp(j,i) = 0
        enddo
      enddo

      
      ntri = 0
      nlev = 0
      ier = 0

      do i=1,ntrimax
        ichild_start(i) = -1
      enddo

      

      if(ifp.eq.1) then
        call gettritree(npatches, norder, npols, ndsc, srccoefs, ntarg,
     1     xyzproxy, itargptr, ntargptr, ntrimax, nlmax, rfac, ntri, 
     2     nlev, ichild_start, da, tricm, trirad, tverts, itrireltmp,
     3     ier)
      else
        allocate(xyztargtmp(3,ntarg))

        do i=1,ntarg
          do j=1,3
            xyztargtmp(j,i) = xyztarg(j,i)
          enddo
        enddo
        call gettritree(npatches, norder, npols, ndsc, srccoefs, ntarg,
     1      xyztargtmp, itargptr, ntargptr, ntrimax, nlmax, rfac, ntri, 
     2      nlev, ichild_start, da, tricm, trirad, tverts, itrireltmp,
     3      ier)        
         deallocate(xyztargtmp)
      endif


      

      if(ier.ne.0) then
        ier = 2
        return
      endif

      allocate(itrirel(ntri,ntarg),itrirelall(ntri))

c
c       transpose the itrirel array
c  

      do itri=1,ntri
        do j=1,ntarg
          itrirel(itri,j) = itrireltmp(j,itri)
        enddo
      enddo


c
c       compute all the geometry info and koornwinder polynomials
c       at relevant triangle given by itrirelall
c  

c
c       get quadrature nodes and weights on the base triangle
c       based on quadrature type
c

      if(intype.eq.1) then
         nqpols = (nqorder+1)*(nqorder+2)/2
         allocate(uvsq(2,nqpols),wts(nqpols))
         allocate(umattmp(nqpols,nqpols),vmattmp(nqpols,nqpols))
      
         call vioreanu_simplex_quad(nqorder, nqpols, uvsq, umattmp,
     1      vmattmp, wts)
         
         deallocate(umattmp,vmattmp)
      endif


      if(intype.eq.2) then
        call triasymq_pts(nqorder,nqpols)
        allocate(uvsq(2,nqpols),wts(nqpols))
        call triasymq(nqorder, tverts(1,1,1), tverts(1,2,1), 
     1    tverts(1,3,1), uvsq, wts, nqpols)    
      endif

cc      call prinf('nqpols=*',nqpols,1)

      npmax = ntri*nqpols
      allocate(sigvals(npols,npmax))
      nds = 12
      if(isd.ne.0) nds = 30
      allocate(srcvals(nds,npmax),qwts(npmax))
      
      allocate(sigmatmp(nppols,npmax),fkervals(npmax))
      allocate(rsigtmp(nppols,npmax))

      allocate(uvtmp(2,nqpols))
      alpha_c = 1
      beta_c = 0


      do itri=1,ntri
        call mapuv_tri(tverts(1,1,itri), nqpols, uvsq, uvtmp)
        istart = (itri-1)*nqpols
        do i=1,nqpols
          ii = istart+i
          call koornf_pols(uvtmp(1,i), norder, npols, sigvals(1,ii),
     1       rat1, rat2, rsc1)
          call koornf_pols(uvtmp(1,i), nporder, nppols, rsigtmp(1,ii),
     1        rat1p, rat2p, rsc1p)
        enddo
      enddo


      do itri=1,npatches
c
cc       for the current patch compute all the geometry info
c
        transa = 'N'
        transb = 'N'
        alpha = 1
        beta = 0
        lda = ndsc
        ldb = npols
        ldc = nds

        call dgemm_guru(transa, transb, ndsc, npmax, npols, alpha,
     1     srccoefs(1,1,itri), lda, sigvals, ldb, beta, srcvals,
     2     ldc)

c
c        compute all the quadrature weights
c

        do i=1,ntri
          istart = (i-1)*nqpols+1
          call get_srcvals_auxinfo_tri(nqpols, wts, isd, 
     1       nds, srcvals(1,istart), da(i), qwts(istart))
        enddo

        do itarg=itargptr(itri),itargptr(itri)+ntargptr(itri)-1
c
c           extract info of geometry on subtriangles relevant
c           for this computation
c
          npts = 0
          do iqtri=1,ntri
            jstart = (iqtri-1)*nqpols
            if(itrirel(iqtri,itarg).eq.1) then
              do i=1,nqpols
                jj = jstart+i
                ii = npts+i
                call fker(srcvals(1,jj), ndtarg, xyztarg(1,itarg),
     1             ndd, dpars, ndz, zpars, ndi, ipars, fkervals(ii))
                fkervals(ii) = fkervals(ii)*qwts(jj)
                do j=1,nppols
                  sigmatmp(j,ii) = rsigtmp(j,jj)
                enddo
              enddo
              npts = npts + nqpols
            endif
          enddo

          call zgemv_guru('n', nppols, npts, alpha_c, sigmatmp,
     1      nppols, fkervals, int8_1, beta_c, cintvals(1,itarg),
     2      int8_1)
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
      subroutine ctriaints_adap(eps, intype,
     1     npatches, norder, npols, isd, ndsc, srccoefs, ndtarg,
     2     ntarg, xyztarg, itargptr, ntargptr, nporder, nppols,
     3     ntrimax, rat1, rat2, rsc1, rat1p, rat2p, rsc1p,
     4     fker, ndd, dpars, ndz, zpars, ndi, ipars, nqorder,
     5     cintvals)
c
c  Compute the integrals in (1) defined at the top of the file
c  using adaptive integration.
c  The refinement stops when the integral computed using 
c  the children of the triangle agree with the integral
c  on the triangle to a tolerance \eps
c
c
c  Input arguments:
c    - eps: real *8
c        requested tolerance. Currently unsued, nqorder and rfac,
c        should be set based on eps.
c    - intype: integer *8
c        node type
c        * intype = 1, rokhlin vioreanu nodes
c        * intype = 2, xiao gimbutas nodes
c    - npatches: integer *8
c        number of patches
c    - norder: integer *8
c        order of discretization on patches
c    - npols: integer *8
c        number of points/basis functions on the 
c        patch = (norder+1)*(norder+2)/2
c    - isd: integer
c        flag for computing and passing second derivative
c        information in srcvals.
c        if isd = 0, leading dim of srcvals = 12, with 
c        xyz, dxyz/du, dxyz/dv, and normals passed
c        else, leading dim of srcvals = 30, with
c        xyz, dxyz/du, dxyz/dv, normals, d2xyz/du2, d2xyz/duv,
c        d2xyz/dv2, det(g), (11,12,22) entries of I, (11,12,22)
c        entries of II, kap1, kap2 where kap1, kap2
c        are principal curvatures, det(g) is the determinant
c        of the metric tensor
c    - ndsc: integer
c        leading dimension of srccoefs array, should be 9
c        if isd = 0, and 18 otherwise
c    - srccoefs: real *8 (ndsc, npols, npatches)
c        Koornwinder expansion coefficients of xyz, d/du (xyz), 
c        d/dv (xyz), (d2(xyz)/du2, d2(xyz)/duv, d2(xyz)/dv2)
c    - ndtarg: integer *8
c        leading dimension of target information
c    - ntarg: integer *8
c        number of targets
c    - xyztarg: real *8 (ndtarg, ntarg)
c        target information
c    - itargptr: integer *8(npatches)
c        Pointer array for determining which targets are relevant
c        for patch i. 
c        xyztargs(:,itargptr(i):itargptr(i)+ntargptr(i)-1)
c        are the relevant set of targets for patch i
c    - ntargptr: integer *8(npatches)
c        number of targets relevant for patch i
c    - nporder: integer *8
c        order of Koornwinder polynomials to be integrated
c    - nppols: integer *8
c        number of polynomials corresponding to 
c        nporder = (nporder+1)*(nporder+2)/2
c    - ntrimax: integer *8
c        Initial guess for maximum number of triangles 
c        allowed in heirarchy of triangles. 
c        Recommended value 3000.
c    - rat1: real *8 (2, 0:norder)
c        recurrence coefficients for Legendre polynomials
c        with discretization order = norder
c    - rat2: real *8 (3, 0:norder, 0:norder)
c        recurrence coefficients for Jacobi polynomials
c        with discretization order = norder
c    - rsc1: real *8 (0:norder, 0:norder)
c        scaling parameters for the Koornwinder polynomials
c        with discretization order = norder
c    - rat1p: real *8 (2, 0:nporder)
c        recurrence coefficients for Legendre polynomials
c        with discretization order = nporder
c    - rat2p: real *8 (3, 0:nporder, 0:nporder)
c        recurrence coefficients for Jacobi polynomials
c        with discretization order = nporder
c    - rsc1p: real *8 (0:nporder, 0:nporder)
c        scaling parameters for the Koornwinder polynomials
c        with discretization order = nporder
c    - fker: function handle
c        function handle for evaluating the kernel k
c        * expected calling sequence
c            fker(x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
c               
c        In this routine the output is expected to be a complex
c        scalar
c    - ndd: integer *8
c        number of real *8/double precision parameters
c    - dpars: real *8 (ndd)
c         real *8/ double precision paramters
c    - ndz: integer *8
c        number of complex *16 parameters
c    - zpars: complex *16(ndz)
c        complex *16 parameters
c    - ndi: integer *8
c        number of integer parameters
c    - ipars: integer *8(ndi)
c        integer parameters
c    - nqorder: integer *8
c        order of quadrature nodes to be used on each triangle
c 
c  Output arguments:
c    - cintvals: complex *16 (nppols,ntarg)
c        the computed integrals

c
      implicit none

c
cc     calling sequence variables
c
      real *8 eps
      integer *8 intype
      integer *8 npatches,norder,npols
      integer *8 nporder,nppols
      integer *8 isd,ndsc
      real *8 srccoefs(ndsc,npols,npatches)
      
      integer *8 ntarg,ndtarg
      real *8 xyztarg(ndtarg,ntarg)
      integer *8 itargptr(npatches)
      integer *8 ntargptr(npatches)
      
      external fker
      integer *8 ndd,ndz,ndi
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer *8 ipars(ndi)

      real *8 rat1(2,0:norder)
      real *8 rat2(3,0:norder,0:norder)
      real *8 rsc1(0:norder,0:norder)

      real *8 rat1p(2,0:nporder)
      real *8 rat2p(3,0:nporder,0:nporder)
      real *8 rsc1p(0:nporder,0:nporder)

      integer *8 nqorder

      complex *16 cintvals(nppols,ntarg)

c
c       tree variables
c
      integer *8 nlmax,ltree
      real *8, allocatable :: tvs(:,:,:),da(:)
      integer *8, allocatable :: ichild_start(:)
      real *8, allocatable :: tvs2(:,:,:),da2(:)
      integer *8, allocatable :: ichild_start2(:)

      integer *8 ntri,ntrimax,nlev,itri,istart,i,j,k
      integer *8 ier,itarg,jj,jstart,npts
      integer *8 iqtri,ii


      integer *8 npmax

      real *8, allocatable :: uvsq(:,:),wts(:),uvtmp(:,:)
      real *8, allocatable :: umattmp(:,:),vmattmp(:,:)
      integer *8 nqpols
      real *8, allocatable :: sigvals(:,:),sigvalsdens(:,:)
      real *8, allocatable :: srcvals(:,:),qwts(:)
      real *8, allocatable :: sigvals2(:,:),sigvalsdens2(:,:)
      real *8, allocatable :: srcvals2(:,:),qwts2(:)
      integer *8 itmp

      character *1 transa,transb
      real *8 alpha,beta
      integer *8 lda,ldb,ldc,nds
      integer *8 nn1,nn2,nn3,nn4,npmax0,ntmaxuse,ntmaxuse0
      integer *8 int8_1
      
      
      int8_1 = 1
c
c      get the tree
c
      nlmax = 20


      ntmaxuse = ntrimax
c
c       for each triangle, we just store three pieces
c       of info
c         triangle vertices
c         area of triangle
c         ichild_start = index for first child, 
c         
c
      allocate(ichild_start(ntrimax),tvs(2,3,ntrimax))
      allocate(da(ntrimax))

      do i=1,ntrimax
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
c       get quadrature nodes and weights on the base triangle
c       based on quadrature type
c

      if(intype.eq.1) then
         nqpols = (nqorder+1)*(nqorder+2)/2
         allocate(uvsq(2,nqpols),wts(nqpols))
         allocate(umattmp(nqpols,nqpols),vmattmp(nqpols,nqpols))
      
         call vioreanu_simplex_quad(nqorder, nqpols, uvsq, umattmp,
     1      vmattmp, wts)
         
         deallocate(umattmp,vmattmp)
      endif

      if(intype.eq.2) then
        call triasymq_pts(nqorder,nqpols)
        allocate(uvsq(2,nqpols),wts(nqpols))
        call triasymq(nqorder, tvs(1,1,1), tvs(1,2,1), tvs(1,3,1),
     1         uvsq, wts, nqpols)    
      endif

      allocate(uvtmp(2,nqpols))
     
      npmax = ntrimax*nqpols
      allocate(sigvals(npols,npmax))
      allocate(sigvalsdens(nppols,npmax))
      nds = 12
      if(isd.ne.0) nds = 30
      allocate(srcvals(nds,npmax),qwts(npmax))
       
c
c      current number of triangles in the adaptive structure
c
      ntri = 1
c
c        intialize sigvals for root triangle
c


      call mapuv_tri(tvs(1,1,1), nqpols, uvsq, uvtmp)
      do i=1,nqpols
        call koornf_pols(uvtmp(1,i), norder, npols,
     1      sigvals(1,i), rat1, rat2, rsc1)
        call koornf_pols(uvtmp(1,i), nporder, nppols,
     1      sigvalsdens(1,i), rat1p, rat2p, rsc1p)
      enddo


      do itri=1,npatches
c
cc       for the current patch compute geometry info for base triangle 
c
        transa = 'N'
        transb = 'N'
        alpha = 1
        beta = 0
        lda = ndsc
        ldb = npols
        ldc = nds 

        call dgemm_guru(transa, transb, ndsc, nqpols, npols, alpha,
     1     srccoefs(1,1,itri), lda, sigvals, ldb, beta, srcvals,
     2     ldc)
        call get_srcvals_auxinfo_tri(nqpols, wts, isd, nds, srcvals, 
     1    da, qwts)

        do itarg=itargptr(itri),itargptr(itri)+ntargptr(itri)-1


 1111     continue
          ier = 0
          call triaadap(eps, nqorder, nqpols, nlmax, ntmaxuse,
     1       ntri, ichild_start, tvs, da, uvsq,wts, 
     2       norder, npols, isd, ndsc, srccoefs(1,1,itri), npmax, 
     3       nds, srcvals, qwts, sigvals, nporder, nppols, sigvalsdens, 
     4       ndtarg, xyztarg(1,itarg), rat1, rat2, rsc1, rat1p, rat2p, 
     5       rsc1p, fker, ndd, dpars, ndz, zpars, ndi, ipars,
     6       cintvals(1,itarg), ier)
           if(ier.eq.4) then
             ntmaxuse0 = ntmaxuse*4
             npmax0 = ntmaxuse0*nqpols
             allocate(sigvals2(npols,npmax))
             allocate(sigvalsdens2(nppols,npmax))
             allocate(srcvals2(nds,npmax),qwts2(npmax))
             allocate(ichild_start2(ntmaxuse),tvs2(2,3,ntmaxuse))
             allocate(da2(ntmaxuse))
             nn1 = npols*npmax
             nn2 = nppols*npmax
             nn3 = nds*npmax
             nn4 = ntri*6
             call dcopy_guru(nn1, sigvals, int8_1, sigvals2, int8_1)
             call dcopy_guru(nn2, sigvalsdens, int8_1, sigvalsdens2,
     1          int8_1)
             call dcopy_guru(nn3, srcvals, int8_1, srcvals2, int8_1)
             call dcopy_guru(npmax, qwts, int8_1, qwts2, int8_1)
             do ii=1,ntri
               ichild_start2(ii) = ichild_start(ii)
             enddo
             call dcopy_guru(nn4,tvs,int8_1,tvs2,int8_1)
             call dcopy_guru(ntri,da,int8_1,da2,int8_1)


             deallocate(sigvals,sigvalsdens,srcvals,qwts,ichild_start)
             deallocate(tvs,da)


             allocate(sigvals(npols,npmax0))
             allocate(sigvalsdens(nppols,npmax0))
             allocate(srcvals(nds,npmax0),qwts(npmax0))
             allocate(ichild_start(ntmaxuse0),tvs(2,3,ntmaxuse0))
             allocate(da(ntmaxuse0))

             do ii=1,ntmaxuse0
               ichild_start(ii) = -1
             enddo

             call dcopy_guru(nn1, sigvals2, int8_1, sigvals, int8_1)
             call dcopy_guru(nn2, sigvalsdens2, int8_1, sigvalsdens,
     1          int8_1)
             call dcopy_guru(nn3, srcvals2, int8_1, srcvals, int8_1)
             call dcopy_guru(npmax, qwts2, int8_1, qwts, int8_1)
             do ii=1,ntri
               ichild_start(ii) = ichild_start2(ii)
             enddo
             call dcopy_guru(nn4, tvs2, int8_1, tvs, int8_1)
             call dcopy_guru(ntri, da2, int8_1, da, int8_1)

             npmax = npmax0
             ntmaxuse = ntmaxuse0
             call prinf('restarting adaptive inetgration with ntri=*',
     1             ntmaxuse,1)


             deallocate(sigvals2,sigvalsdens2,srcvals2,ichild_start2)
             deallocate(tvs2,da2,qwts2)
             goto 1111
           endif
        enddo

        do i=1,ntri
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
      subroutine triaadap(eps, m, kpols, nlmax, ntmax, ntri,
     1             ichild_start, tvs, da, uvsq, wts,
     2             norder, npols, isd, ndsc, srccoefs, npmax, nds, 
     3             srcvals, qwts, sigvals, nporder, nppols, 
     4             sigvalsdens, ndtarg, xt, rat1, rat2, rsc1,
     5             rat1p, rat2p, rsc1p,
     6             fker, ndd, dpars, ndz, zpars, ndi,
     7             ipars, cintall, ier)

c
c  Compute the integrals in (1) defined at the top of the file
c  using adaptive integration.
c  This is an intermediate routine which initializes the guru
c  adaptive integration routine
c
c
c  Input arguments:
c    - eps: real *8
c        requested tolerance. Currently unsued, nqorder and rfac,
c        should be set based on eps.
c    - m: integer *8
c        order for quadrature nodes 
c    - kpols: integer *8
c        number of quadrature nodes 
c    - nlmax: integer *8
c        max number of levels
c    - ntmax: integer *8
c        max number of triangles
c    - ntri: integer *8
c        current number of triangles in adaptive integration
c    - ichild_start: integer *8(ntmax)
c        ichild_start(i) is the first child of traingle i
c    - tvs: real *8(2,3,ntmax)
c        vertices of hierarchy of triangles
c    - da: real *8(ntmax)
c        area of triangles
c    - uvsq: real *8(kpols)
c        integration nodes on standard triangle
c    - wts: real *8(kpols)
c        integration weights on standard triangle
c    - norder: integer *8
c        order of discretization on patches
c    - npols: integer *8
c        number of points/basis functions on the 
c        patch = (norder+1)*(norder+2)/2
c    - isd: integer *8
c        flag for computing and passing second derivative
c        information in srcvals.
c        if isd = 0, leading dim of srcvals = 12, with 
c        xyz, dxyz/du, dxyz/dv, and normals passed
c        else, leading dim of srcvals = 30, with
c        xyz, dxyz/du, dxyz/dv, normals, d2xyz/du2, d2xyz/duv,
c        d2xyz/dv2, det(g), (11,12,22) entries of I, (11,12,22)
c        entries of II, kap1, kap2 where kap1, kap2
c        are principal curvatures, det(g) is the determinant
c        of the metric tensor
c    - ndsc: integer *8
c        leading dimension of srccoefs array, should be 9
c        if isd = 0, and 18 otherwise
c    - srccoefs: real *8 (ndsc, npols)
c        Koornwinder expansion coefficients of xyz, d/du (xyz), 
c        d/dv (xyz), (d2(xyz)/du2, d2(xyz)/duv, d2(xyz)/dv2)
c    - npmax: integer *8
c        max number of points = ntmax*kpols
c    - nds: integer *8
c        leading dimension of srcvals array, if isd = 0,
c        then nds must be 12, and if isd = 1, then
c        nds must be 30
c    - srcvals: real *8(nds,npmax)
c        geometry info on heirarchy of meshes
c    - qwts: real *8(npmax)
c        quadrature weights 
c    - sigvals: real *8(npols,npmax) - 
c        koornwinder polynomials computed along the adaptive grid
c        of order = norder
c    - nporder: integer *8
c        order of Koornwinder polynomials to be integrated
c    - nppols: integer *8
c        number of polynomials corresponding to 
c        nporder = (nporder+1)*(nporder+2)/2
c    - sigvalsdens: real *8(nppols,npmax) - 
c        koornwinder polynomials computed along the adaptive grid
c        of order = nporder
c    - ndtarg: integer *8
c        leading dimension of target information
c    - xt: real *8 (ndtarg)
c        target information
c    - rat1: real *8 (2, 0:norder)
c        recurrence coefficients for Legendre polynomials
c        with discretization order = norder
c    - rat2: real *8 (3, 0:norder, 0:norder)
c        recurrence coefficients for Jacobi polynomials
c        with discretization order = norder
c    - rsc1: real *8 (0:norder, 0:norder)
c        scaling parameters for the Koornwinder polynomials
c        with discretization order = norder
c    - rat1p: real *8 (2, 0:nporder)
c        recurrence coefficients for Legendre polynomials
c        with discretization order = nporder
c    - rat2p: real *8 (3, 0:nporder, 0:nporder)
c        recurrence coefficients for Jacobi polynomials
c        with discretization order = nporder
c    - rsc1p: real *8 (0:nporder, 0:nporder)
c        scaling parameters for the Koornwinder polynomials
c        with discretization order = nporder
c    - fker: function handle
c        function handle for evaluating the kernel k
c        * expected calling sequence
c            fker(x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
c               
c        In this routine the output is expected to be a complex
c        scalar
c    - ndd: integer *8
c        number of real *8/double precision parameters
c    - dpars: real *8 (ndd)
c         real *8/ double precision paramters
c    - ndz: integer *8
c        number of complex *16 parameters
c    - zpars: complex *16(ndz)
c        complex *16 parameters
c    - ndi: integer *8
c        number of integer parameters
c    - ipars: integer *8(ndi)
c        integer parameters
c    - nqorder: integer *8
c        order of quadrature nodes to be used on each triangle
c 
c  Output arguments:
c    - cintall: complex *16 (nppols)
c        the computed integrals
c    - ier: integer *8
c        error code
c        * ier = 0, successful execution
c        * ier = 4, too few triangles, try with more triangles
c         

      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, allocatable :: istack(:)
      integer *8 ichild_start(ntmax)
      real *8 da(ntmax)
      real *8 tvs(2,3,ntmax), uvsq(2,kpols),wts(kpols)
      integer *8 nproclist0, nproclist
      integer *8 idone
      integer *8 isd, ndsc, nds
      real *8 srccoefs(ndsc,npols)
      real *8 sigvals(npols,npmax)
      real *8 sigvalsdens(nppols,npmax)
      real *8 srcvals(nds,*),qwts(npmax)
      complex *16, allocatable :: xkernvals(:)
      real *8 xt(ndtarg)
      complex *16 cintall(nppols),fval
      complex *16, allocatable :: ctmp(:)
      complex *16, allocatable :: cvals(:,:)

      integer *8 ndd,ndz,ndi
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer *8 ipars(ndi)
      integer *8 ier

      real *8 rat1(2,0:norder),rat2(3,0:norder,0:norder)
      real *8 rsc1(0:norder,0:norder)

      real *8 rat1p(2,0:nporder)
      real *8 rat2p(3,0:nporder,0:nporder)
      real *8 rsc1p(0:nporder,0:nporder)


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
         call fker(srcvals(1,i), ndtarg, xt, ndd, dpars, ndz, zpars, 
     1     ndi, ipars, fval)
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


      call triaadap_main(eps, kpols, nlmax, ntmax, ntri, ichild_start,
     1      tvs, da, uvsq, wts, norder, npols, isd, ndsc, srccoefs,
     2      npmax, nds, srcvals, qwts, sigvals, nporder, nppols,
     3      sigvalsdens, ndtarg, xt, rat1, rat2, rsc1,
     3      rat1p, rat2p, rsc1p, fker, ndd, dpars,
     3      ndz, zpars, ndi, ipars, cvals, istack, nproclist0,
     4      xkernvals, cintall, ier)
      
      return
      end
c
c
c
c
c
       
      subroutine triaadap_main(eps, kpols, nlmax, ntmax, 
     1    ntri, ichild_start, tvs, da, uvsq, wts, norder, npols,
     2    isd, ndsc, srccoefs, npmax, nds, srcvals, qwts, sigvals, 
     3    nporder, nppols, sigvalsdens, ndtarg, xt, rat1, rat2, rsc1,
     4    rat1p, rat2p, rsc1p, fker, ndd, dpars, ndz,
     3    zpars, ndi, ipars, cvals, istack, nproclist0, xkernvals,
     4    cintall, ier)
      

      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8 istack(*),nproclist0
      integer *8 ichild_start(ntmax)
      integer *8 nporder,nppols
      real *8 da(ntmax)
      real *8 tvs(2,3,ntmax), uvsq(2,kpols),wts(kpols)
      integer *8 nproclist
      integer *8 idone
      integer *8 isd, ndsc, nds
      real *8 srccoefs(ndsc,npols)
      real *8 sigvals(npols,npmax)
      real *8 sigvalsdens(nppols,npmax)
      real *8 qwts(npmax)
      complex *16 xkernvals(npmax)
      real *8 xt(ndtarg)
      real *8 srcvals(nds,*)
      complex *16 cintall(nppols),fval,ctmp(nppols)
      complex *16 cvals(nppols,ntmax)

      integer *8 ndd,ndz,ndi
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer *8 ipars(ndi)

      real *8 rat1(2,0:norder),rat2(3,0:norder,0:norder)
      real *8 rsc1(0:norder,0:norder)

      real *8 rat1p(2,0:nporder),rat2p(3,0:nporder,0:nporder)
      real *8 rsc1p(0:nporder,0:nporder)


      real *8, allocatable :: uvtmp(:,:)
      character *1 transa,transb
      integer *8 lda,ldb,ldc
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
          itri = istack(iproc)

c
c           check to see if triangle already has 
c           children, if not, set children
c           and compute necessary info
c
          if(ichild_start(itri).eq.-1) then

c
c            current triangle doesn't have children,
c            compute necessary info
c

            if(ntri+4.gt.ntmax) then
c               print *, "Too many triangles in ctriaadap"
c               print *, "Exiting without computing anything"
               ier = 4
               return
            endif
            
            ichild_start(itri) = ntri+1
            call gettrichildren(tvs(1,1,itri), tvs(1,1,ntri+1),
     1             tvs(1,1,ntri+2), tvs(1,1,ntri+3), tvs(1,1,ntri+4))

            rr = 0.25d0*da(itri)
            do j=ntri+1,ntri+4
              da(j) = rr
              call mapuv_tri(tvs(1,1,j), kpols, uvsq, uvtmp)
              istart = (j-1)*kpols+1
              do i=1,kpols
                 ii = istart+i-1
                call koornf_pols(uvtmp(1,i), norder, npols, 
     1            sigvals(1,ii), rat1, rat2, rsc1)
                call koornf_pols(uvtmp(1,i), nporder, nppols,
     1             sigvalsdens(1,ii), rat1p, rat2p, rsc1p)
              enddo
              
              transa = 'N'
              transb = 'N'
              alpha = 1
              beta = 0
              lda = ndsc
              ldb = npols
              ldc = nds

              call dgemm_guru(transa, transb, ndsc, kpols, npols, alpha,
     1           srccoefs, lda, sigvals(1,istart), ldb, beta,
     2           srcvals(1,istart), ldc)
              call get_srcvals_auxinfo_tri(kpols, wts, isd, 
     1          nds, srcvals(1,istart), rr, qwts(istart)) 
            enddo
            ntri = ntri+4
          endif
        
c
cc           compute xkernvals
c
          itric1 = ichild_start(itri)
          istart = (itric1-1)*kpols
          do j=1,kfine
            jj=j+istart
            call fker(srcvals(1,jj), ndtarg, xt, ndd, dpars,
     1         ndz, zpars, ndi, ipars, fval)
            xkernvals(jj) = fval*qwts(jj)
          enddo

          nfunev = nfunev + kfine

c
cc         subtract contribution of current 
c          triangle
          do isig=1,ksigpols
            cintall(isig) = cintall(isig)-cvals(isig,itri)
            ctmp(isig) = 0
          enddo

c
cc        add in contributions of children triangles
c
          do itric=itric1,itric1+3
            do isig=1,ksigpols
              cvals(isig,itric) = 0
            enddo

            istart = (itric-1)*kpols
            do k=1,kpols
              ii = istart+k
              do isig=1,ksigpols
                cvals(isig,itric) = cvals(isig,itric)+xkernvals(ii)*
     1                                  sigvalsdens(isig,ii)
              enddo
            enddo
              
            do isig=1,ksigpols
              cintall(isig) = cintall(isig) + cvals(isig,itric)
              ctmp(isig) = ctmp(isig)+cvals(isig,itric)
            enddo
          enddo

c
cc        compare integral of children to integral of parent
c         to determine if the children need to be refined further
c
          errmax = 0
          do isig=1,ksigpols
            if(abs(ctmp(isig)-cvals(isig,itri)).gt.errmax) 
     1          errmax = abs(ctmp(isig)-cvals(isig,itri))
          enddo

          if(errmax.gt.eps) then
            idone = 0
            
            do j=1,4
              istack(nproclist0+nproclist+j) = itric1+j-1
            enddo
            nproclist = nproclist+4
          endif
c
cc        end of looping over all triangles at current stage
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
      subroutine ctriaints_comb(eps, intype,
     1     npatches, norder, npols, isd, ndsc, srccoefs, ndtarg,
     2     ntarg, xyztarg, ifp, xyzproxy,
     3     itargptr, ntargptr, nporder, nppols, ntrimax,
     4     rat1, rat2, rsc1, rat1p, rat2p, rsc1p,
     3     fker, ndd, dpars, ndz, zpars, ndi, ipars,
     4     nqorder, rfac, cintvals)

c
c  Compute the integrals in (1) defined at the top of the file
c  using a combination of distance based criterion
c  followed by adaptive integration.
c
c  A triangle is subdivided when xyzproxy associated
c  with the triangle is separated from the centroid
c  of the triangle by at least rfac*r_{t} where r_{t}
c  is the radius of the smallest sphere  
c  The refinement stops when the integral computed using 
c  the children of the triangle agree with the integral
c  on the triangle to a tolerance \eps
c
c
c  Input arguments:
c    - eps: real *8
c        requested tolerance. Currently unsued, nqorder and rfac,
c        should be set based on eps.
c    - intype: integer *8
c        node type
c        * intype = 1, rokhlin vioreanu nodes
c        * intype = 2, xiao gimbutas nodes
c    - npatches: integer *8
c        number of patches
c    - norder: integer *8
c        order of discretization on patches
c    - npols: integer *8
c        number of points/basis functions on the 
c        patch = (norder+1)*(norder+2)/2
c    - isd: integer *8
c        flag for computing and passing second derivative
c        information in srcvals.
c        if isd = 0, leading dim of srcvals = 12, with 
c        xyz, dxyz/du, dxyz/dv, and normals passed
c        else, leading dim of srcvals = 30, with
c        xyz, dxyz/du, dxyz/dv, normals, d2xyz/du2, d2xyz/duv,
c        d2xyz/dv2, det(g), (11,12,22) entries of I, (11,12,22)
c        entries of II, kap1, kap2 where kap1, kap2
c        are principal curvatures, det(g) is the determinant
c        of the metric tensor
c    - ndsc: integer *8
c        leading dimension of srccoefs array, should be 9
c        if isd = 0, and 18 otherwise
c    - srccoefs: real *8 (ndsc, npols, npatches)
c        Koornwinder expansion coefficients of xyz, d/du (xyz), 
c        d/dv (xyz), (d2(xyz)/du2, d2(xyz)/duv, d2(xyz)/dv2)
c    - ndtarg: integer *8
c        leading dimension of target information
c    - ntarg: integer *8
c        number of targets
c    - xyztarg: real *8 (ndtarg, ntarg)
c        target information
c    - ifp: integer *8
c        flag for using proxy points for refining based 
c        on distance criterion
c    - xyzproxy: real *8 (3,*)
c        Proxy target points for refining triangles based
c        distance criterion, distance to proxy
c        point will be used instead of distance to target
c        should be of size (3,ntarg), if ifp=1
c    - itargptr: integer *8(npatches)
c        Pointer array for determining which targets are relevant
c        for patch i. 
c        xyztargs(:,itargptr(i):itargptr(i)+ntargptr(i)-1)
c        are the relevant set of targets for patch i
c    - ntargptr: integer *8(npatches)
c        number of targets relevant for patch i
c    - nporder: integer *8
c        order of Koornwinder polynomials to be integrated
c    - nppols: integer *8
c        number of polynomials corresponding to 
c        nporder = (nporder+1)*(nporder+2)/2
c    - ntrimax: integer *8
c        Initial guess for maximum number of triangles 
c        allowed in heirarchy of triangles. 
c        Recommended value 3000.
c    - rat1: real *8 (2, 0:norder)
c        recurrence coefficients for Legendre polynomials
c        with discretization order = norder
c    - rat2: real *8 (3, 0:norder, 0:norder)
c        recurrence coefficients for Jacobi polynomials
c        with discretization order = norder
c    - rsc1: real *8 (0:norder, 0:norder)
c        scaling parameters for the Koornwinder polynomials
c        with discretization order = norder
c    - rat1p: real *8 (2, 0:nporder)
c        recurrence coefficients for Legendre polynomials
c        with discretization order = nporder
c    - rat2p: real *8 (3, 0:nporder, 0:nporder)
c        recurrence coefficients for Jacobi polynomials
c        with discretization order = nporder
c    - rsc1p: real *8 (0:nporder, 0:nporder)
c        scaling parameters for the Koornwinder polynomials
c        with discretization order = nporder
c    - fker: function handle
c        function handle for evaluating the kernel k
c        * expected calling sequence
c            fker(x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
c               
c        In this routine the output is expected to be a complex
c        scalar
c    - ndd: integer *8
c        number of real *8/double precision parameters
c    - dpars: real *8 (ndd)
c         real *8/ double precision paramters
c    - ndz: integer *8
c        number of complex *16 parameters
c    - zpars: complex *16(ndz)
c        complex *16 parameters
c    - ndi: integer *8
c        number of integer parameters
c    - ipars: integer *8(ndi)
c        integer parameters
c    - nqorder: integer *8
c        order of quadrature nodes to be used on each triangle
c    - rfac: real *8
c        scaling factor for determining when a triangle 
c        needs refinement
c 
c  Output arguments:
c    - cintvals: complex *16 (nppols, ntarg)
c        the computed integrals

      implicit none

c
cc     calling sequence variables
c
      real *8 eps
      integer *8 intype,ifp
      integer *8 npatches,norder,npols
      integer *8 nporder,nppols
      integer *8 isd, ndsc
      real *8 srccoefs(ndsc,npols,npatches)
      
      integer *8 ntarg,ndtarg
      real *8 xyztarg(ndtarg,ntarg)
      real *8 xyzproxy(3,*)
      integer *8 itargptr(npatches)
      integer *8 ntargptr(npatches)
      
      external fker
      integer *8 ndd,ndz,ndi
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer *8 ipars(ndi)

      integer *8 nqorder
      real *8 rfac

      real *8 rat1(2,0:norder),rat2(3,0:norder,0:norder)
      real *8 rsc1(0:norder,0:norder)

      real *8 rat1p(2,0:nporder),rat2p(3,0:nporder,0:nporder)
      real *8 rsc1p(0:nporder,0:nporder)

      complex *16 cintvals(nppols,ntarg)
      

c
cc      temporary variables
c
      real *8, allocatable :: tricm(:,:,:)
      real *8, allocatable :: trirad(:,:)
      integer *8, allocatable :: itrireltmp(:,:)
      integer *8, allocatable :: itrirel(:,:)
      integer *8, allocatable :: itrirelall(:)
      real *8, allocatable :: tverts(:,:,:)
      integer *8, allocatable :: ichild_start(:)
      integer *8, allocatable :: ichild_start0(:)

      integer *8 ntri,ntrimax,nlev,itri,istart,i,j
      integer *8 ier,itarg,jj,jstart,nlmax,npts
      integer *8 iqtri,ii

      integer *8 npmax

      real *8, allocatable :: uvsq(:,:),wts(:),uvtmp(:,:)
      real *8, allocatable :: umattmp(:,:),vmattmp(:,:)
      real *8, allocatable :: da(:)
      integer *8 nqpols
      real *8, allocatable :: sigvals(:,:)
      real *8, allocatable :: sigvalsdens(:,:)
      real *8, allocatable :: srcvals(:,:)
      real *8, allocatable :: qwts(:)
      complex *16, allocatable :: sigmatmp(:,:)
      real *8, allocatable :: xyztargtmp(:,:)
      real *8 xyztmp(3)
      integer *8 itmp, nds

      complex *16, allocatable :: xkernvals(:)
      integer *8, allocatable :: istack(:)

      complex *16, allocatable :: cvals(:,:)
      
      integer *8 nproclist0,nproclist,ntri0,npts0
      complex *16 zz

      
      complex *16 ima

      character *1 transa,transb
      
      real *8 alpha,beta
      integer *8 lda,ldb,ldc

      data ima/(0.0d0,1.0d0)/

      allocate(cvals(nppols,ntrimax))
      allocate(istack(2*ntrimax))





cc      max number of levels
c
      nlmax = 20
      allocate(tricm(3,npatches,ntrimax),
     1   trirad(npatches,ntrimax),tverts(2,3,ntrimax))
      allocate(itrireltmp(ntarg,ntrimax))
      allocate(da(ntrimax),ichild_start(ntrimax))
      allocate(ichild_start0(ntrimax))

      do i=1,ntrimax
        do j=1,ntarg
          itrireltmp(j,i) = 0
        enddo
        
      enddo

      
      ntri = 0
      nlev = 0
      ier = 0

      do i=1,ntrimax
        ichild_start(i) = -1
      enddo

      

      if(ifp.eq.1) then
        call gettritree(npatches, norder, npols, ndsc, srccoefs, ntarg,
     1     xyzproxy, itargptr, ntargptr, ntrimax, nlmax, rfac, 
     2     ntri, nlev, ichild_start, da, tricm, trirad, tverts,
     3     itrireltmp, ier)
      else
        allocate(xyztargtmp(3,ntarg))
        do i=1,ntarg
          do j=1,3
            xyztargtmp(j,i) = xyztarg(j,i)
          enddo
        enddo
        call gettritree(npatches, norder, npols, ndsc, srccoefs, ntarg,
     1     xyztargtmp, itargptr, ntargptr, ntrimax, nlmax, rfac,
     2     ntri, nlev, ichild_start, da, tricm, trirad, tverts,
     3     itrireltmp, ier)
        deallocate(xyztargtmp)
      endif

      ntri0 = ntri

      do i=1,ntri0
        ichild_start0(i) = ichild_start(i)
      enddo


      

      if(ier.ne.0) then
        call prinf('Could not allocate triangle tree*',i,0)
        call prinf('Exiting without computing anything*',i,0)
        call prinf('Press 0 to continue*',i,0)
        call prinf('Press 1 to exit routine and continue*',i,0)
        call prinf('Press 2 to stop*',i,0)
        read *, ier
        if(ier.eq.2) stop
        if(ier.eq.1) return
      endif

      allocate(itrirel(ntri0,ntarg),itrirelall(ntri0))

c
c       transpose the itrirel array
c  

      do itri=1,ntri0
        do j=1,ntarg
          itrirel(itri,j) = itrireltmp(j,itri)
        enddo
      enddo


c
c       compute all the geometry info and koornwinder polynomials
c       at relevant triangle given by itrirelall
c  

c
c       get quadrature nodes and weights on the base triangle
c       based on quadrature type
c

      if(intype.eq.1) then
         nqpols = (nqorder+1)*(nqorder+2)/2
         allocate(uvsq(2,nqpols),wts(nqpols))
         allocate(umattmp(nqpols,nqpols),vmattmp(nqpols,nqpols))
      
         call vioreanu_simplex_quad(nqorder, nqpols, uvsq, umattmp,
     1      vmattmp, wts)
         
         deallocate(umattmp,vmattmp)
      endif


      if(intype.eq.2) then
        call triasymq_pts(nqorder,nqpols)
        allocate(uvsq(2,nqpols),wts(nqpols))
        call triasymq(nqorder, tverts(1,1,1), tverts(1,2,1),
     1     tverts(1,3,1), uvsq, wts, nqpols)    
      endif

cc      call prinf('nqpols=*',nqpols,1)

      npts0 = ntri0*nqpols
      npmax = ntrimax*nqpols
      allocate(sigvals(npols,npmax))
      allocate(sigvalsdens(nppols,npmax))
      allocate(qwts(npmax))

      nds = 12
      if(isd.gt.0) nds = 30
      allocate(srcvals(nds,npmax))

      

      allocate(xkernvals(npmax))

      allocate(uvtmp(2,nqpols))


      do itri=1,ntri0
        call mapuv_tri(tverts(1,1,itri),nqpols,uvsq,uvtmp)
        istart = (itri-1)*nqpols
        do i=1,nqpols
          ii = istart+i
          call koornf_pols(uvtmp(1,i), norder, npols, sigvals(1,ii),
     1           rat1, rat2, rsc1)
          call koornf_pols(uvtmp(1,i), nporder, nppols, 
     1       sigvalsdens(1,ii), rat1p, rat2p, rsc1p)
        enddo
      enddo


      do itri=1,npatches
        ntri = ntri0
c
cc       for the current patch compute all the geometry info
c

        transa = 'N'
        transb = 'N'
        alpha = 1.0d0
        beta = 0.0d0
        lda = ndsc
        ldb = npols
        ldc = nds

        call dgemm_guru(transa, transb, ndsc, npts0, npols, alpha,
     1     srccoefs(1,1,itri), lda, sigvals, ldb, beta, srcvals, ldc)



c
c        compute all the quadrature weights
c

        do i=1,ntri
          istart = (i-1)*nqpols+1
          call get_srcvals_auxinfo_tri(nqpols, wts, isd, 
     1      nds, srcvals(1,istart), da(i), qwts(istart)) 
        enddo

        do itarg=itargptr(itri),itargptr(itri)+ntargptr(itri)-1
          do i=1,nppols
            cintvals(i,itarg) = 0
          enddo
c
c           extract info of geometry on subtriangles relevant
c           for this computation
c
          nproclist0 = 0
          do iqtri=1,ntri0
            jstart = (iqtri-1)*nqpols
            if(itrirel(iqtri,itarg).eq.1) then
              nproclist0 = nproclist0 + 1
              istack(nproclist0) = iqtri
c
c
c              compute the integral due to itri
c

              do i=1,nppols
                cvals(i,iqtri) = 0
              enddo
              do i=1,nqpols
                jj = jstart+i
                call fker(srcvals(1,jj), ndtarg, xyztarg(1,itarg),
     1              ndd, dpars, ndz, zpars, ndi, ipars, xkernvals(jj))
                xkernvals(jj) = xkernvals(jj)*qwts(jj)
                
                do j=1,nppols
                  zz = xkernvals(jj)*sigvalsdens(j,jj)
                  cvals(j,iqtri) = cvals(j,iqtri)+zz
                  cintvals(j,itarg)=cintvals(j,itarg)+zz
                enddo
              enddo
            endif
          enddo
c
c           done with initial computation of 
c

          ier = 0
          call triaadap_main(eps, nqpols, nlmax, ntrimax, ntri, 
     1     ichild_start, tverts, da, uvsq, wts, norder, npols,
     2     isd, ndsc, srccoefs(1,1,itri), npmax, nds, srcvals, qwts, 
     3     sigvals, nporder, nppols, sigvalsdens, ndtarg, 
     4     xyztarg(1,itarg), rat1, rat2, rsc1, rat1p, rat2p, rsc1p,
     5     fker, ndd, dpars, ndz, zpars, ndi, ipars, cvals,
     6     istack, nproclist0, xkernvals, cintvals(1,itarg), ier)

      
        enddo

        do i=1,ntri0
          ichild_start(i) = ichild_start0(i)
        enddo
        do i=ntri0+1,ntri
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
      subroutine ctriaints_wnodes_vec(npatches, norder, npols, isd, 
     1   ndsc, srccoefs, ndtarg, ntarg, xyztarg, itargptr, ntargptr,
     2   nporder, nppols, fker, nd, ndd, dpars, ndz, zpars, ndi, ipars,
     3   nqpts, qnodes, wts, cintvals)
c
c  Compute the integrals in (1) defined at the top of the file
c  using a prescribed set of nodes and weights.
c
c  Input arguments:
c    - npatches: integer *8
c        number of patches
c    - norder: integer *8
c        order of discretization on patches
c    - npols: integer *8
c        number of points/basis functions on the 
c        patch = (norder+1)*(norder+2)/2
c    - isd: integer *8
c        flag for computing and passing second derivative
c        information in srcvals.
c        if isd = 0, leading dim of srcvals = 12, with 
c        xyz, dxyz/du, dxyz/dv, and normals passed
c        else, leading dim of srcvals = 30, with
c        xyz, dxyz/du, dxyz/dv, normals, d2xyz/du2, d2xyz/duv,
c        d2xyz/dv2, det(g), (11,12,22) entries of I, (11,12,22)
c        entries of II, kap1, kap2 where kap1, kap2
c        are principal curvatures, det(g) is the determinant
c        of the metric tensor
c    - ndsc: integer *8
c        leading dimension of srccoefs array, should be 9
c        if isd = 0, and 18 otherwise
c    - srccoefs: real *8 (ndsc, npols, npatches)
c        Koornwinder expansion coefficients of xyz, d/du (xyz), 
c        d/dv (xyz), (d2(xyz)/du2, d2(xyz)/duv, d2(xyz)/dv2)
c    - ndtarg: integer *8
c        leading dimension of target information
c    - ntarg: integer *8
c        number of targets
c    - xyztarg: real *8 (ndtarg, ntarg)
c        target information
c    - itargptr: integer *8(npatches)
c        Pointer array for determining which targets are relevant
c        for patch i. 
c        xyztargs(:,itargptr(i):itargptr(i)+ntargptr(i)-1)
c        are the relevant set of targets for patch i
c    - ntargptr: integer *8(npatches)
c        number of targets relevant for patch i
c    - nporder: integer *8
c        order of Koornwinder polynomials to be integrated
c    - nppols: integer *8
c        number of polynomials corresponding to 
c        nporder = (nporder+1)*(nporder+2)/2
c    - fker: function handle
c        function handle for evaluating the kernel k
c        * expected calling sequence
c            fker(nd,x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
c               
c        In this routine the output is expected to be a complex
c        vector
c    - nd: integer *8
c        number of kernels
c    - ndd: integer *8
c        number of real *8/double precision parameters
c    - dpars: real *8 (ndd)
c         real *8/ double precision paramters
c    - ndz: integer *8
c        number of complex *16 parameters
c    - zpars: complex *16(ndz)
c        complex *16 parameters
c    - ndi: integer *8
c        number of integer parameters
c    - ipars: integer *8(ndi)
c        integer parameters
c    - nqpts: integer *8
c        number of quadrature points
c    - qnodes: real *8 (2,nqpts)
c        u,v location of the quadrature nodes
c    - wts: real *8 nqpts
c        quadrature weights
c 
c  Output arguments:
c    - cintvals: complex *16 (nppols, ntarg)
c        the computed integrals
c 

      implicit none
      integer *8, intent(in) :: npatches, norder, npols, ndtarg
      integer *8, intent(in) :: nporder, nppols
      integer *8, intent(in) :: isd, ndsc
      real *8, intent(in) :: srccoefs(ndsc,npols,npatches)
      real *8, intent(in) :: xyztarg(ndtarg,ntarg)
      integer *8, intent(in) :: ntarg
      integer *8, intent(in) :: itargptr(npatches), ntargptr(npatches)
      integer *8, intent(in) :: nd, ndd, ndz, ndi
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      integer *8, intent(in) :: ipars(ndi),nqpts
      real *8, intent(in) :: qnodes(2,nqpts),wts(nqpts)

      complex *16, intent(out) :: cintvals(nd,nppols,ntarg)

      real *8, allocatable :: rat1(:,:), rat2(:,:,:), rsc1(:,:)
      real *8, allocatable :: rat1p(:,:), rat2p(:,:,:), rsc1p(:,:)
      real *8, allocatable :: srcvals(:,:),qwts(:)
      real *8, allocatable :: rsigvals(:,:),rsigtmp(:)
      complex *16, allocatable :: sigvals(:,:)
      complex *16, allocatable :: xkernvals(:,:,:)
      complex *16, allocatable :: cinttmp(:,:,:)

      integer *8 i,ipatch,j,lda,ldb,itarg,ldc,ntarg0,ii, idim
      integer *8 nds

      complex *16 fval(nd)
      real *8 da,ra

      character *1 transa,transb
      real *8 alpha,beta
      complex *16 alpha_c,beta_c

      external fker

c
c       initialize koornwinder polynomials
c

      
      allocate(rat1(2,0:norder),rat2(3,0:norder,0:norder))
      allocate(rsc1(0:norder,0:norder))
      call koornf_init(norder,rat1,rat2,rsc1) 

      allocate(rat1p(2,0:nporder),rat2p(3,0:nporder,0:nporder))
      allocate(rsc1p(0:nporder,0:nporder))
      call koornf_init(nporder,rat1p,rat2p,rsc1p) 

      allocate(rsigtmp(nppols))

      allocate(sigvals(nppols,nqpts),rsigvals(npols,nqpts))

      nds = 12
      if(isd.ne.0) nds = 30
      allocate(srcvals(nds,nqpts),qwts(nqpts))

      do i=1,nqpts
        call koornf_pols(qnodes(1,i),norder,npols,rsigvals(1,i),
     1        rat1,rat2,rsc1)
        call koornf_pols(qnodes(1,i),nporder,nppols,rsigtmp,rat1p,
     2        rat2p,rsc1p)
        do j=1,nppols
          sigvals(j,i) = rsigtmp(j)
        enddo
      enddo




      do ipatch=1,npatches
c
c        compute srcvals
c
        transa = 'N'
        transb = 'N'
        alpha = 1
        beta = 0
        lda = ndsc
        ldb = npols
        ldc = nds

        da = 1

        call dgemm_guru(transa, transb, ndsc, nqpts, npols, alpha,
     1     srccoefs(1,1,ipatch), lda, rsigvals, ldb, beta,
     2     srcvals, ldc)

   
        call get_srcvals_auxinfo_tri(nqpts, wts, isd, nds, srcvals, 
     1    da, qwts)
c
c
c          compute the kernel values for all the targets
c
        
        ntarg0 = ntargptr(ipatch)
        allocate(xkernvals(nd, ntarg0, nqpts))
        allocate(cinttmp(nppols, nd, ntarg0))
        do j=1,nqpts
          do itarg=itargptr(ipatch),itargptr(ipatch)+ntarg0-1
            ii = itarg - itargptr(ipatch)+1
            call fker(nd, srcvals(1,j), ndtarg, xyztarg(1,itarg),
     1         ndd, dpars, ndz, zpars, ndi, ipars, fval)
            xkernvals(:,ii,j) = fval(:)*qwts(j)
          enddo
        enddo


        transa = 'n'
        transb = 't'
        alpha_c = 1
        beta_c = 0

        call zgemm_guru(transa, transb, nppols, nd*ntarg0, nqpts,
     1     alpha_c, sigvals, nppols, xkernvals, nd*ntarg0, beta_c, 
     2     cinttmp, nppols)

        do itarg=itargptr(ipatch), itargptr(ipatch)+ntarg0-1
           do idim=1,nd
             do j=1,nppols
                ii = itarg -itargptr(ipatch)+1
                cintvals(idim,j,itarg) = cinttmp(j,idim,ii)
             enddo
          enddo 
        enddo
      
        deallocate(xkernvals)
c

      enddo

      

      return
      end


      subroutine ctriaints_dist_vec(eps, intype,
     1     npatches, norder, npols, isd, ndsc, srccoefs, ndtarg,
     2     ntarg, xyztarg, ifp, xyzproxy,
     3     itargptr, ntargptr, nporder, nppols, ntrimax, 
     4     rat1, rat2, rsc1, rat1p, rat2p, rsc1p,
     3     fker, nd, ndd, dpars, ndz, zpars, ndi, ipars,
     4     nqorder, rfac, cintvals)
c
c  Compute the integrals in (1) defined at the top of the file
c  using a distance based criterion to subdivide a triangle.
c  A triangle is subdivided when xyzproxy associated
c  with the triangle is separated from the centroid
c  of the triangle by at least rfac*r_{t} where r_{t}
c  is the radius of the smallest sphere.
c
c  Vectorized version of ctriaints_dist 
c
c  Input arguments:
c    - eps: real *8
c        requested tolerance. Currently unsued, nqorder and rfac,
c        should be set based on eps.
c    - intype: integer *8
c        node type
c        * intype = 1, rokhlin vioreanu nodes
c        * intype = 2, xiao gimbutas nodes
c    - npatches: integer *8
c        number of patches
c    - norder: integer *8
c        order of discretization on patches
c    - npols: integer *8
c        number of points/basis functions on the 
c        patch = (norder+1)*(norder+2)/2
c    - isd: integer *8
c        flag for computing and passing second derivative
c        information in srcvals.
c        if isd = 0, leading dim of srcvals = 12, with 
c        xyz, dxyz/du, dxyz/dv, and normals passed
c        else, leading dim of srcvals = 30, with
c        xyz, dxyz/du, dxyz/dv, normals, d2xyz/du2, d2xyz/duv,
c        d2xyz/dv2, det(g), (11,12,22) entries of I, (11,12,22)
c        entries of II, kap1, kap2 where kap1, kap2
c        are principal curvatures, det(g) is the determinant
c        of the metric tensor
c    - ndsc: integer *8
c        leading dimension of srccoefs array, should be 9
c        if isd = 0, and 18 otherwise
c    - srccoefs: real *8 (ndsc, npols, npatches)
c        Koornwinder expansion coefficients of xyz, d/du (xyz), 
c        d/dv (xyz), (d2(xyz)/du2, d2(xyz)/duv, d2(xyz)/dv2)
c    - ndtarg: integer *8
c        leading dimension of target information
c    - ntarg: integer *8
c        number of targets
c    - xyztarg: real *8 (ndtarg, ntarg)
c        target information
c    - ifp: integer *8
c        flag for using proxy points for refining based 
c        on distance criterion
c    - xyzproxy: real *8 (3,*)
c        Proxy target points for refining triangles based
c        distance criterion, distance to proxy
c        point will be used instead of distance to target
c        should be of size (3,ntarg), if ifp=1
c    - itargptr: integer *8(npatches)
c        Pointer array for determining which targets are relevant
c        for patch i. 
c        xyztargs(:,itargptr(i):itargptr(i)+ntargptr(i)-1)
c        are the relevant set of targets for patch i
c    - ntargptr: integer *8(npatches)
c        number of targets relevant for patch i
c    - nporder: integer *8
c        order of Koornwinder polynomials to be integrated
c    - nppols: integer *8
c        number of polynomials corresponding to 
c        nporder = (nporder+1)*(nporder+2)/2
c    - ntrimax: integer *8
c        maximum number of triangles allowed in heirarchy
c        of triangles. Routine will return without
c        computing anything and an error code, if ntrimax
c        is too small. Recommended value 3000.
c    - rat1: real *8 (2, 0:norder)
c        recurrence coefficients for Legendre polynomials
c        with discretization order = norder
c    - rat2: real *8 (3, 0:norder, 0:norder)
c        recurrence coefficients for Jacobi polynomials
c        with discretization order = norder
c    - rsc1: real *8 (0:norder, 0:norder)
c        scaling parameters for the Koornwinder polynomials
c        with discretization order = norder
c    - rat1p: real *8 (2, 0:nporder)
c        recurrence coefficients for Legendre polynomials
c        with discretization order = nporder
c    - rat2p: real *8 (3, 0:nporder, 0:nporder)
c        recurrence coefficients for Jacobi polynomials
c        with discretization order = nporder
c    - rsc1p: real *8 (0:nporder, 0:nporder)
c        scaling parameters for the Koornwinder polynomials
c        with discretization order = nporder
c    - fker: function handle
c        function handle for evaluating the kernel k
c        * expected calling sequence
c            fker(nd,x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
c               
c        In this routine the output is expected to be a complex
c        vector
c    - nd: integer *8
c        number of kernels
c    - ndd: integer *8
c        number of real *8/double precision parameters
c    - dpars: real *8 (ndd)
c         real *8/ double precision paramters
c    - ndz: integer *8
c        number of complex *16 parameters
c    - zpars: complex *16(ndz)
c        complex *16 parameters
c    - ndi: integer *8
c        number of integer parameters
c    - ipars: integer *8(ndi)
c        integer parameters
c    - nqorder: integer *8
c        order of quadrature nodes to be used on each triangle
c    - rfac: real *8
c        scaling factor for determining when a triangle 
c        needs refinement
c    - cintvals: complex *16 (nd, nppols, ntarg)
c        the computed integrals


      implicit none

c
cc     calling sequence variables
c
      real *8 eps
      integer *8 intype,ifp,nd
      integer *8 npatches,norder,npols
      integer *8 nporder,nppols
      integer *8 isd, ndsc
      real *8 srccoefs(ndsc,npols,npatches)
      
      integer *8 ntarg,ndtarg
      real *8 xyztarg(ndtarg,ntarg)
      real *8 xyzproxy(3,*)
      integer *8 itargptr(npatches)
      integer *8 ntargptr(npatches)
      
      external fker
      integer *8 ndd,ndi,ndz
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer *8 ipars(ndi)

      integer *8 nqorder
      real *8 rfac
      real *8 rat1(2,0:norder),rat2(3,0:norder,0:norder)
      real *8 rsc1(0:norder,0:norder)

      real *8 rat1p(2,0:nporder),rat2p(3,0:nporder,0:nporder)
      real *8 rsc1p(0:nporder,0:nporder)

      complex *16 cintvals(nd,nppols,ntarg)
      

c
cc      temporary variables
c
      real *8, allocatable :: tricm(:,:,:)
      real *8, allocatable :: trirad(:,:)
      integer *8, allocatable :: itrireltmp(:,:)
      integer *8, allocatable :: itrirel(:,:)
      integer *8, allocatable :: itrirelall(:)
      real *8, allocatable :: tverts(:,:,:)
      integer *8, allocatable :: ichild_start(:)

      integer *8 ntri,ntrimax,nlev,itri,istart,i,j
      integer *8 ier,itarg,jj,jstart,nlmax,npts
      integer *8 iqtri,ii

      integer *8 npmax

      real *8, allocatable :: uvsq(:,:),wts(:),uvtmp(:,:)
      real *8, allocatable :: umattmp(:,:),vmattmp(:,:)
      real *8, allocatable :: da(:)
      integer *8 nqpols
      real *8, allocatable :: sigvals(:,:)
      real *8, allocatable :: srcvals(:,:),qwts(:)
      real *8, allocatable :: xyztargtmp(:,:)
      complex *16, allocatable :: fkervals(:,:),sigmatmp(:,:)
      real *8, allocatable :: rsigtmp(:,:)
      real *8 xyztmp(3)
      integer *8 itmp,idim,nds
      complex *16 ima

      character *1 transa,transb
      double precision :: alpha,beta
      complex *16 alpha_c, beta_c
      integer *8 lda,ldb,ldc
      

      data ima/(0.0d0,1.0d0)/





cc      max number of levels
c
      nlmax = 20
      allocate(tricm(3,npatches,ntrimax),
     1   trirad(npatches,ntrimax),tverts(2,3,ntrimax))
      allocate(itrireltmp(ntarg,ntrimax))
      allocate(da(ntrimax),ichild_start(ntrimax))

      do i=1,ntrimax
        do j=1,ntarg
          itrireltmp(j,i) = 0
        enddo
      enddo

      
      ntri = 0
      nlev = 0
      ier = 0

      do i=1,ntrimax
        ichild_start(i) = -1
      enddo

      

      if(ifp.eq.1) then
        call gettritree(npatches, norder, npols, ndsc, srccoefs, ntarg,
     1    xyzproxy, itargptr, ntargptr, ntrimax, nlmax, rfac, ntri,
     2    nlev, ichild_start, da, tricm, trirad, tverts, itrireltmp,
     3    ier)
      else
        allocate(xyztargtmp(3,ntarg))
        do i=1,ntarg
          do j=1,3
            xyztargtmp(j,i) = xyztarg(j,i)
          enddo
        enddo

        call gettritree(npatches, norder, npols, ndsc, srccoefs, ntarg,
     1     xyztargtmp, itargptr, ntargptr, ntrimax, nlmax, rfac, 
     2     ntri, nlev, ichild_start, da, tricm, trirad, tverts,
     3     itrireltmp, ier)
        deallocate(xyztargtmp)
      endif


      

      if(ier.ne.0) then
        call prinf('Could not allocate triangle tree*',i,0)
        call prinf('Exiting without computing anything*',i,0)
        call prinf('Press 0 to continue*',i,0)
        call prinf('Press 1 to exit routine and continue*',i,0)
        call prinf('Press 2 to stop*',i,0)
        read *, ier
        if(ier.eq.2) stop
        if(ier.eq.1) return
      endif

      allocate(itrirel(ntri,ntarg),itrirelall(ntri))

c
c       transpose the itrirel array
c  

      do itri=1,ntri
        do j=1,ntarg
          itrirel(itri,j) = itrireltmp(j,itri)
        enddo
      enddo


c
c       compute all the geometry info and koornwinder polynomials
c       at relevant triangle given by itrirelall
c  

c
c       get quadrature nodes and weights on the base triangle
c       based on quadrature type
c

      if(intype.eq.1) then
         nqpols = (nqorder+1)*(nqorder+2)/2
         allocate(uvsq(2,nqpols),wts(nqpols))
         allocate(umattmp(nqpols,nqpols),vmattmp(nqpols,nqpols))
      
         call vioreanu_simplex_quad(nqorder, nqpols, uvsq, umattmp,
     1      vmattmp, wts)
         
         deallocate(umattmp,vmattmp)
      endif


      if(intype.eq.2) then
        call triasymq_pts(nqorder,nqpols)
        allocate(uvsq(2,nqpols),wts(nqpols))
        call triasymq(nqorder, tverts(1,1,1), tverts(1,2,1),
     1    tverts(1,3,1), uvsq, wts, nqpols)    
      endif

cc      call prinf('nqpols=*',nqpols,1)

      npmax = ntri*nqpols
      allocate(sigvals(npols,npmax))
      nds = 12
      if(isd.ne.0) nds = 30
      allocate(srcvals(nds,npmax),qwts(npmax))

      allocate(sigmatmp(nppols,npmax),fkervals(nd,npmax))
      allocate(rsigtmp(nppols,npmax))

      allocate(uvtmp(2,nqpols))


      do itri=1,ntri
        call mapuv_tri(tverts(1,1,itri),nqpols,uvsq,uvtmp)
        istart = (itri-1)*nqpols
        do i=1,nqpols
          ii = istart+i
          call koornf_pols(uvtmp(1,i), norder, npols, sigvals(1,ii),
     1        rat1, rat2, rsc1)
          call koornf_pols(uvtmp(1,i), nporder, nppols, rsigtmp(1,ii),
     1        rat1p, rat2p, rsc1p)
        enddo
      enddo



      do itri=1,npatches
c
cc       for the current patch compute all the geometry info
c
        transa = 'N'
        transb = 'N'
        alpha = 1
        beta = 0
        lda = ndsc
        ldb = npols
        ldc = nds

        call dgemm_guru(transa, transb, ndsc, npmax, npols, alpha,
     1     srccoefs(1,1,itri), lda, sigvals, ldb, beta, srcvals, ldc)


c
c        compute all the quadrature weights
c

        do i=1,ntri
          istart = (i-1)*nqpols+1
          call get_srcvals_auxinfo_tri(nqpols, wts, isd, 
     1       nds, srcvals(1,istart), da(i), qwts(istart))
        enddo

        do itarg=itargptr(itri),itargptr(itri)+ntargptr(itri)-1
c
c           extract info of geometry on subtriangles relevant
c           for this computation
c
          npts = 0
          do iqtri=1,ntri
            jstart = (iqtri-1)*nqpols
            if(itrirel(iqtri,itarg).eq.1) then
              do i=1,nqpols
                jj = jstart+i
                ii = npts+i
                call fker(nd, srcvals(1,jj), ndtarg, xyztarg(1,itarg),
     1             ndd, dpars, ndz, zpars, ndi, ipars, fkervals(1,ii))
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
          call zgemm_guru(transa, transb, nd, nppols, npts, alpha_c,
     1      fkervals, nd, sigmatmp, nppols, beta_c, cintvals(1,1,itarg),
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
      subroutine ctriaints_adap_vec(eps, intype,
     1     npatches, norder, npols, isd, ndsc, srccoefs, ndtarg,
     2     ntarg, xyztarg, itargptr, ntargptr, nporder, nppols, 
     3     ntrimax, rat1, rat2, rsc1, rat1p, rat2p, rsc1p,
     4     fker, nd, ndd, dpars, ndz, zpars, ndi, ipars, nqorder,
     5     cintvals)
c
c  Compute the integrals in (1) defined at the top of the file
c  using adaptive integration.
c  The refinement stops when the integral computed using 
c  the children of the triangle agree with the integral
c  on the triangle to a tolerance \eps.
c
c  Vectorized version of ctriaints_adap
c
c
c  Input arguments:
c    - eps: real *8
c        requested tolerance. Currently unsued, nqorder and rfac,
c        should be set based on eps.
c    - intype: integer *8
c        node type
c        * intype = 1, rokhlin vioreanu nodes
c        * intype = 2, xiao gimbutas nodes
c    - npatches: integer *8
c        number of patches
c    - norder: integer *8
c        order of discretization on patches
c    - npols: integer *8
c        number of points/basis functions on the 
c        patch = (norder+1)*(norder+2)/2
c    - isd: integer *8
c        flag for computing and passing second derivative
c        information in srcvals.
c        if isd = 0, leading dim of srcvals = 12, with 
c        xyz, dxyz/du, dxyz/dv, and normals passed
c        else, leading dim of srcvals = 30, with
c        xyz, dxyz/du, dxyz/dv, normals, d2xyz/du2, d2xyz/duv,
c        d2xyz/dv2, det(g), (11,12,22) entries of I, (11,12,22)
c        entries of II, kap1, kap2 where kap1, kap2
c        are principal curvatures, det(g) is the determinant
c        of the metric tensor
c    - ndsc: integer *8
c        leading dimension of srccoefs array, should be 9
c        if isd = 0, and 18 otherwise
c    - srccoefs: real *8 (ndsc, npols, npatches)
c        Koornwinder expansion coefficients of xyz, d/du (xyz), 
c        d/dv (xyz), (d2(xyz)/du2, d2(xyz)/duv, d2(xyz)/dv2)
c    - ndtarg: integer *8
c        leading dimension of target information
c    - ntarg: integer *8
c        number of targets
c    - xyztarg: real *8 (ndtarg, ntarg)
c        target information
c    - itargptr: integer *8(npatches)
c        Pointer array for determining which targets are relevant
c        for patch i. 
c        xyztargs(:,itargptr(i):itargptr(i)+ntargptr(i)-1)
c        are the relevant set of targets for patch i
c    - ntargptr: integer *8(npatches)
c        number of targets relevant for patch i
c    - nporder: integer *8
c        order of Koornwinder polynomials to be integrated
c    - nppols: integer *8
c        number of polynomials corresponding to 
c        nporder = (nporder+1)*(nporder+2)/2
c    - ntrimax: integer *8
c        Initial guess for maximum number of triangles 
c        allowed in heirarchy of triangles. 
c        Recommended value 3000.
c    - rat1: real *8 (2, 0:norder)
c        recurrence coefficients for Legendre polynomials
c        with discretization order = norder
c    - rat2: real *8 (3, 0:norder, 0:norder)
c        recurrence coefficients for Jacobi polynomials
c        with discretization order = norder
c    - rsc1: real *8 (0:norder, 0:norder)
c        scaling parameters for the Koornwinder polynomials
c        with discretization order = norder
c    - rat1p: real *8 (2, 0:nporder)
c        recurrence coefficients for Legendre polynomials
c        with discretization order = nporder
c    - rat2p: real *8 (3, 0:nporder, 0:nporder)
c        recurrence coefficients for Jacobi polynomials
c        with discretization order = nporder
c    - rsc1p: real *8 (0:nporder, 0:nporder)
c        scaling parameters for the Koornwinder polynomials
c        with discretization order = nporder
c    - fker: function handle
c        function handle for evaluating the kernel k
c        * expected calling sequence
c            fker(nd, x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
c               
c        In this routine the output is expected to be a complex
c        vector
c    - nd: integer *8
c        number of kernels
c    - ndd: integer *8
c        number of real *8/double precision parameters
c    - dpars: real *8 (ndd)
c         real *8/ double precision paramters
c    - ndz: integer *8
c        number of complex *16 parameters
c    - zpars: complex *16(ndz)
c        complex *16 parameters
c    - ndi: integer *8
c        number of integer parameters
c    - ipars: integer *8(ndi)
c        integer parameters
c    - nqorder: integer *8
c        order of quadrature nodes to be used on each triangle
c 
c  Output arguments:
c    - cintvals: complex *16 (nd,nppols,ntarg)
c        the computed integrals

c
      implicit none

c
cc     calling sequence variables
c
      real *8 eps
      integer *8 intype,nd
      integer *8 npatches,norder,npols
      integer *8 nporder,nppols
      integer *8 isd, ndsc
      real *8 srccoefs(9,npols,npatches)
      
      integer *8 ntarg,ndtarg
      real *8 xyztarg(ndtarg,ntarg)
      integer *8 itargptr(npatches)
      integer *8 ntargptr(npatches)
      
      external fker
      integer *8 ndd,ndi,ndz
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer *8 ipars(ndi)

      real *8 rat1(2,0:norder)
      real *8 rat2(3,0:norder,0:norder)
      real *8 rsc1(0:norder,0:norder)

      real *8 rat1p(2,0:norder)
      real *8 rat2p(3,0:norder,0:norder)
      real *8 rsc1p(0:norder,0:norder)

      integer *8 nqorder

      complex *16 cintvals(nd,nppols,ntarg)

c
c       tree variables
c
      integer *8 nlmax,ltree
      real *8, allocatable :: tvs(:,:,:),da(:)
      integer *8, allocatable :: ichild_start(:)
      real *8, allocatable :: tvs2(:,:,:),da2(:)
      integer *8, allocatable :: ichild_start2(:)

      integer *8 ntri,ntrimax,nlev,itri,istart,i,j,k
      integer *8 ier,itarg,jj,jstart,npts
      integer *8 iqtri,ii


      integer *8 npmax

      real *8, allocatable :: uvsq(:,:),wts(:),uvtmp(:,:)
      real *8, allocatable :: umattmp(:,:),vmattmp(:,:)
      integer *8 nqpols
      real *8, allocatable :: sigvals(:,:),sigvalsdens(:,:)
      real *8, allocatable :: srcvals(:,:),qwts(:)
      real *8, allocatable :: sigvals2(:,:),sigvalsdens2(:,:)
      real *8, allocatable :: srcvals2(:,:),qwts2(:)
      integer *8 itmp

      character *1 transa,transb
      real *8 alpha,beta
      integer *8 lda,ldb,ldc,nds
      integer *8 nn1,nn2,nn3,nn4,npmax0,ntmaxuse,ntmaxuse0
      integer *8 int8_1
      
      
      int8_1 = 1
c
c      get the tree
c
      nlmax = 20
      ntmaxuse = ntrimax

      ltree = 2*(nlmax+1) + 6*ntri

c
c       for each triangle, we just store three pieces
c       of info
c         triangle vertices
c         area of triangle
c         ichild_start = index for first child, 
c         
c
      allocate(ichild_start(ntrimax),tvs(2,3,ntrimax))
      allocate(da(ntrimax))

      do i=1,ntrimax
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
c       get quadrature nodes and weights on the base triangle
c       based on quadrature type
c

      if(intype.eq.1) then
         nqpols = (nqorder+1)*(nqorder+2)/2
         allocate(uvsq(2,nqpols),wts(nqpols))
         allocate(umattmp(nqpols,nqpols),vmattmp(nqpols,nqpols))
      
         call vioreanu_simplex_quad(nqorder, nqpols, uvsq, umattmp,
     1      vmattmp, wts)
         
         deallocate(umattmp,vmattmp)
      endif

      if(intype.eq.2) then
        call triasymq_pts(nqorder,nqpols)
        allocate(uvsq(2,nqpols),wts(nqpols))
        call triasymq(nqorder, tvs(1,1,1), tvs(1,2,1), tvs(1,3,1),
     1         uvsq, wts, nqpols)    
      endif

      allocate(uvtmp(2,nqpols))
     
      npmax = ntrimax*nqpols
      allocate(sigvals(npols,npmax),sigvalsdens(nppols,npmax))
      nds = 12
      if(isd.ne.0) nds = 30
      allocate(srcvals(nds,npmax),qwts(npmax))


c
c      current number of triangles in the adaptive structure
c
      ntri = 1
c
c        intialize sigvals for root triangle
c


      call mapuv_tri(tvs(1,1,1),nqpols,uvsq,uvtmp)
      do i=1,nqpols
        call koornf_pols(uvtmp(1,i), norder, npols,
     1      sigvals(1,i), rat1, rat2, rsc1)
        call koornf_pols(uvtmp(1,i), nporder, nppols,
     1      sigvalsdens(1,i), rat1p, rat2p, rsc1p)
      enddo


      do itri=1,npatches
c
cc       for the current patch compute geometry info for base triangle 
c
        transa = 'N'
        transb = 'N'
        alpha = 1
        beta = 0
        lda = ndsc
        ldb = npols
        ldc = nds

        call dgemm_guru(transa, transb, ndsc, nqpols, npols, alpha,
     1     srccoefs(1,1,itri), lda, sigvals, ldb, beta, srcvals,
     2     ldc)
        call get_srcvals_auxinfo_tri(nqpols, wts, isd, nds, srcvals, da, 
     1     qwts)

        do itarg=itargptr(itri),itargptr(itri)+ntargptr(itri)-1

 1111     continue

          call triaadap_vec(eps, nqorder, nqpols, nlmax, ntmaxuse,
     1       ntri, ichild_start, tvs, da, uvsq, wts, 
     2       norder, npols, isd, ndsc, srccoefs(1,1,itri), npmax, nds, 
     2       srcvals, qwts, sigvals, nporder, nppols, sigvalsdens, 
     3       ndtarg, xyztarg(1,itarg), rat1, rat2, rsc1, rat1p, rat2p, 
     3       rsc1p, fker, nd, ndd, dpars, ndz, zpars, ndi, ipars,
     3       cintvals(1,1,itarg), ier)
           if(ier.eq.4) then
             ntmaxuse0 = ntmaxuse*4
             npmax0 = ntmaxuse0*nqpols
             allocate(sigvals2(npols,npmax))
             allocate(sigvalsdens2(nppols,npmax))
             allocate(srcvals2(nds,npmax),qwts2(npmax))
             allocate(ichild_start2(ntmaxuse),tvs2(2,3,ntmaxuse))
             allocate(da2(ntmaxuse))
             nn1 = npols*npmax
             nn2 = nppols*npmax
             nn3 = nds*npmax
             nn4 = ntri*6
             call dcopy_guru(nn1, sigvals, int8_1, sigvals2, int8_1)
             call dcopy_guru(nn2, sigvalsdens, int8_1, sigvalsdens2,
     1          int8_1)
             call dcopy_guru(nn3, srcvals, int8_1, srcvals2, int8_1)
             call dcopy_guru(npmax, qwts, int8_1, qwts2, int8_1)
             do ii=1,ntri
               ichild_start2(ii) = ichild_start(ii)
             enddo
             call dcopy_guru(nn4, tvs, int8_1, tvs2, int8_1)
             call dcopy_guru(ntri, da, int8_1, da2, int8_1)


             deallocate(sigvals,sigvalsdens,srcvals,qwts,ichild_start)
             deallocate(tvs,da)


             allocate(sigvals(npols,npmax0))
             allocate(sigvalsdens(nppols,npmax0))
             allocate(srcvals(nds,npmax0),qwts(npmax0))
             allocate(ichild_start(ntmaxuse0),tvs(2,3,ntmaxuse0))
             allocate(da(ntmaxuse0))

             do ii=1,ntmaxuse0
               ichild_start(ii) = -1
             enddo

             call dcopy_guru(nn1, sigvals2, int8_1, sigvals, int8_1)
             call dcopy_guru(nn2, sigvalsdens2, int8_1, sigvalsdens,
     1          int8_1)
             call dcopy_guru(nn3, srcvals2, int8_1, srcvals, int8_1)
             call dcopy_guru(npmax, qwts2, int8_1, qwts, int8_1)
             do ii=1,ntri
               ichild_start(ii) = ichild_start2(ii)
             enddo
             call dcopy_guru(nn4, tvs2, int8_1, tvs, int8_1)
             call dcopy_guru(ntri, da2, int8_1, da, int8_1)

             npmax = npmax0
             ntmaxuse = ntmaxuse0
c             call prinf('restarting adaptive inetgration with ntri=*',
c     1             ntmaxuse,1)


             deallocate(sigvals2,sigvalsdens2,srcvals2,ichild_start2)
             deallocate(tvs2,da2,qwts2)
             goto 1111
           endif
        enddo

        do i=1,ntri
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
      subroutine triaadap_vec(eps, m, kpols, nlmax, ntmax, ntri,
     1             ichild_start, tvs, da, uvsq, wts,
     1             norder, npols, isd, ndsc, srccoefs, npmax, nds, 
     2             srcvals, qwts, sigvals, nporder, nppols, sigvalsdens, 
     3             ndtarg, xt, rat1, rat2, rsc1, rat1p, rat2p, rsc1p,
     3             fker, nd, ndd, dpars, ndz, zpars, ndi,
     3             ipars, cintall, ier)
c
c  Compute the integrals in (1) defined at the top of the file
c  using adaptive integration.
c  This is an intermediate routine which initializes the guru
c  adaptive integration routine.
c
c  Vectorized version of triaadap.
c
c
c  Input arguments:
c    - eps: real *8
c        requested tolerance. Currently unsued, nqorder and rfac,
c        should be set based on eps.
c    - m: integer *8
c        order for quadrature nodes 
c    - kpols: integer *8
c        number of quadrature nodes 
c    - nlmax: integer *8
c        max number of levels
c    - ntmax: integer *8
c        max number of triangles
c    - ntri: integer *8
c        current number of triangles in adaptive integration
c    - ichild_start: integer *8(ntmax)
c        ichild_start(i) is the first child of traingle i
c    - tvs: real *8(2,3,ntmax)
c        vertices of hierarchy of triangles
c    - da: real *8(ntmax)
c        area of triangles
c    - uvsq: real *8(kpols)
c        integration nodes on standard triangle
c    - wts: real *8(kpols)
c        integration weights on standard triangle
c    - norder: integer *8
c        order of discretization on patches
c    - npols: integer *8
c        number of points/basis functions on the 
c        patch = (norder+1)*(norder+2)/2
c    - isd: integer *8
c        flag for computing and passing second derivative
c        information in srcvals.
c        if isd = 0, leading dim of srcvals = 12, with 
c        xyz, dxyz/du, dxyz/dv, and normals passed
c        else, leading dim of srcvals = 30, with
c        xyz, dxyz/du, dxyz/dv, normals, d2xyz/du2, d2xyz/duv,
c        d2xyz/dv2, det(g), (11,12,22) entries of I, (11,12,22)
c        entries of II, kap1, kap2 where kap1, kap2
c        are principal curvatures, det(g) is the determinant
c        of the metric tensor
c    - ndsc: integer *8
c        leading dimension of srccoefs array, should be 9
c        if isd = 0, and 18 otherwise
c    - srccoefs: real *8 (ndsc, npols)
c        Koornwinder expansion coefficients of xyz, d/du (xyz), 
c        d/dv (xyz), (d2(xyz)/du2, d2(xyz)/duv, d2(xyz)/dv2)
c    - npmax: integer *8
c        max number of points = ntmax*kpols
c    - nds: integer *8
c        leading dimension of srcvals array
c    - srcvals: real *8(nds,npmax)
c        geometry info on heirarchy of meshes
c    - qwts: real *8(npmax)
c        quadrature weights 
c    - sigvals: real *8(npols,npmax) - 
c        koornwinder polynomials computed along the adaptive grid
c        of order = norder
c    - nporder: integer *8
c        order of Koornwinder polynomials to be integrated
c    - nppols: integer *8
c        number of polynomials corresponding to 
c        nporder = (nporder+1)*(nporder+2)/2
c    - sigvalsdens: real *8(nppols,npmax) - 
c        koornwinder polynomials computed along the adaptive grid
c        of order = nporder
c    - ndtarg: integer *8
c        leading dimension of target information
c    - xt: real *8 (ndtarg)
c        target information
c    - rat1: real *8 (2, 0:norder)
c        recurrence coefficients for Legendre polynomials
c        with discretization order = norder
c    - rat2: real *8 (3, 0:norder, 0:norder)
c        recurrence coefficients for Jacobi polynomials
c        with discretization order = norder
c    - rsc1: real *8 (0:norder, 0:norder)
c        scaling parameters for the Koornwinder polynomials
c        with discretization order = norder
c    - rat1p: real *8 (2, 0:nporder)
c        recurrence coefficients for Legendre polynomials
c        with discretization order = nporder
c    - rat2p: real *8 (3, 0:nporder, 0:nporder)
c        recurrence coefficients for Jacobi polynomials
c        with discretization order = nporder
c    - rsc1p: real *8 (0:nporder, 0:nporder)
c        scaling parameters for the Koornwinder polynomials
c        with discretization order = nporder
c    - fker: function handle
c        function handle for evaluating the kernel k
c        * expected calling sequence
c            fker(nd, x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
c               
c        In this routine the output is expected to be a complex
c        vector
c    - nd: integer *8
c        number of kernels
c    - ndd: integer *8
c        number of real *8/double precision parameters
c    - dpars: real *8 (ndd)
c         real *8/ double precision paramters
c    - ndz: integer *8
c        number of complex *16 parameters
c    - zpars: complex *16(ndz)
c        complex *16 parameters
c    - ndi: integer *8
c        number of integer parameters
c    - ipars: integer *8(ndi)
c        integer parameters
c    - nqorder: integer *8
c        order of quadrature nodes to be used on each triangle
c 
c  Output arguments:
c    - cintall: complex *16 (nd,nppols)
c        the computed integrals
c    - ier: integer *8
c        error code
c        * ier = 0, successful execution
c        * ier = 4, too few triangles, try with more triangles


      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, allocatable :: istack(:)
      integer *8 ichild_start(ntmax),nd
      integer *8 nporder,nppols
      real *8 da(ntmax)
      real *8 tvs(2,3,ntmax), uvsq(2,kpols),wts(kpols)
      integer *8 nproclist0, nproclist
      integer *8 idone
      integer *8 isd, ndsc, nds
      real *8 srccoefs(ndsc,npols)
      real *8 sigvals(npols,npmax)
      real *8 sigvalsdens(nppols,npmax)
      real *8 srcvals(nds,*),qwts(npmax)
      complex *16, allocatable :: xkernvals(:,:)
      real *8 xt(ndtarg)
      complex *16 cintall(nd,npols),fval(nd)
      complex *16, allocatable :: cvals(:,:,:)

      integer *8 ndd,ndz,ndi
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer *8 ipars(ndi),idim

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
         call fker(nd,srcvals(1,i), ndtarg, xt, ndd, dpars, ndz, zpars,
     1        ndi, ipars, fval)
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


      call triaadap_main_vec(eps, kpols, nlmax, ntmax, ntri, 
     1       ichild_start, tvs, da, uvsq, wts, norder, npols,
     2       isd, ndsc, srccoefs, npmax, nds, srcvals, qwts, sigvals, 
     3       nporder, nppols, sigvalsdens, ndtarg, xt, rat1, rat2, rsc1,
     4       rat1p, rat2p, rsc1p, fker, nd, ndd, dpars, ndz,
     5       zpars, ndi, ipars, cvals, istack, nproclist0,
     6       xkernvals, cintall, ier)
      
      return
      end
c
c
c
c
c
       
      subroutine triaadap_main_vec(eps, kpols, nlmax, ntmax, ntri,
     1      ichild_start, tvs, da, uvsq, wts, norder, npols, 
     2      isd, ndsc, srccoefs, npmax, nds, srcvals, qwts, sigvals, 
     3      nporder, nppols, sigvalsdens, ndtarg, xt, rat1, rat2, rsc1, 
     4      rat1p, rat2p, rsc1p, fker, nd, ndd, dpars, ndz, zpars, ndi, 
     5      ipars, cvals, istack, nproclist0, xkernvals, cintall, ier)
      

      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8 istack(*),nproclist0,nd
      integer *8 ichild_start(ntmax)
      integer *8 nporder,nppols
      real *8 da(ntmax)
      real *8 tvs(2,3,ntmax), uvsq(2,kpols),wts(kpols)
      integer *8 nproclist
      integer *8 idone
      integer *8 isd, ndsc, nds
      real *8 srccoefs(ndsc,npols)
      real *8 sigvals(npols,npmax)
      real *8 sigvalsdens(nppols,npmax)
      real *8 qwts(npmax)
      complex *16 xkernvals(nd,npmax)
      real *8 xt(ndtarg)
      real *8 srcvals(nds,*)
      complex *16 cintall(nd,nppols),fval(nd),ctmp(nd,nppols)
      complex *16 cvals(nd,nppols,ntmax)

      integer *8 ndd,ndz,ndi
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer *8 ipars(ndi)

      real *8 rat1(2,0:norder),rat2(3,0:norder,0:norder)
      real *8 rsc1(0:norder,0:norder)

      real *8 rat1p(2,0:nporder),rat2p(3,0:nporder,0:nporder)
      real *8 rsc1p(0:nporder,0:nporder)

      real *8, allocatable :: uvtmp(:,:)
      character *1 transa,transb
      integer *8 lda,ldb,ldc
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
          itri = istack(iproc)

c
c           check to see if triangle already has 
c           children, if not, set children
c           and compute necessary info
c
          if(ichild_start(itri).eq.-1) then

c
c            current triangle doesn't have children,
c            compute necessary info
c

            if(ntri+4.gt.ntmax) then
c               print *, "Too many triangles in ctriaadap"
c               print *, "Exiting without computing anything"
               ier = 4
               return
            endif
            
            ichild_start(itri) = ntri+1
            call gettrichildren(tvs(1,1,itri), tvs(1,1,ntri+1),
     1             tvs(1,1,ntri+2), tvs(1,1,ntri+3), tvs(1,1,ntri+4))

            rr = 0.25d0*da(itri)
            do j=ntri+1,ntri+4
              da(j) = rr
              call mapuv_tri(tvs(1,1,j), kpols, uvsq, uvtmp)
              istart = (j-1)*kpols+1
              do i=1,kpols
                 ii = istart+i-1
                call koornf_pols(uvtmp(1,i), norder, npols, 
     1           sigvals(1,ii), rat1, rat2, rsc1)
                call koornf_pols(uvtmp(1,i), nporder, nppols,
     1             sigvalsdens(1,ii), rat1p, rat2p, rsc1p)
              enddo
              
              transa = 'N'
              transb = 'N'
              alpha = 1
              beta = 0
              lda = ndsc
              ldb = npols
              ldc = nds

              call dgemm_guru(transa, transb, ndsc, kpols, npols, alpha,
     1           srccoefs, lda, sigvals(1,istart), ldb, beta,
     2           srcvals(1,istart), ldc)
              call get_srcvals_auxinfo_tri(kpols, wts, isd, 
     1          nds, srcvals(1,istart), rr, qwts(istart))

            enddo
            ntri = ntri+4
          endif
        
c
cc           compute xkernvals
c
          itric1 = ichild_start(itri)
          istart = (itric1-1)*kpols
          do j=1,kfine
            jj=j+istart
            call fker(nd,srcvals(1,jj), ndtarg, xt, ndd, dpars,
     1         ndz, zpars, ndi, ipars, fval)
            do idim=1,nd
              xkernvals(idim,jj) = fval(idim)*qwts(jj)
            enddo
          enddo
cc          call prin2('xkernvals=*',xkernvals(istart+1),kfine)

          nfunev = nfunev + kfine

c
cc         subtract contribution of current 
c          triangle
          do isig=1,ksigpols
            do idim=1,nd
              cintall(idim,isig) = cintall(idim,isig)-
     1            cvals(idim,isig,itri)
              ctmp(idim,isig) = 0
            enddo
          enddo

c
cc        add in contributions of children triangles
c
          do itric=itric1,itric1+3
            do isig=1,ksigpols
              do idim=1,nd
                cvals(idim,isig,itric) = 0
              enddo
            enddo

            istart = (itric-1)*kpols
            do k=1,kpols
              ii = istart+k
              do isig=1,ksigpols
                do idim=1,nd
                  cvals(idim,isig,itric) = cvals(idim,isig,itric)+
     1               xkernvals(idim,ii)*sigvalsdens(isig,ii)
                enddo
              enddo
            enddo
              
            do isig=1,ksigpols
              do idim=1,nd
                cintall(idim,isig) = cintall(idim,isig) + 
     1              cvals(idim,isig,itric)
                ctmp(idim,isig) = ctmp(idim,isig)+
     1              cvals(idim,isig,itric)
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
              if(abs(ctmp(idim,isig)-cvals(idim,isig,itri)).gt.errmax) 
     1          errmax = abs(ctmp(idim,isig)-cvals(idim,isig,itri))
            enddo
          enddo

          if(errmax.gt.eps) then
            idone = 0
            
            do j=1,4
              istack(nproclist0+nproclist+j) = itric1+j-1
            enddo
            nproclist = nproclist+4
          endif
c
cc        end of looping over all triangles at current stage
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
      subroutine ctriaints_comb_vec(eps, intype,
     1     npatches, norder, npols, isd, ndsc, srccoefs, ndtarg,
     2     ntarg, xyztarg, ifp, xyzproxy,
     3     itargptr, ntargptr, nporder, nppols, ntrimax,
     4     rat1, rat2, rsc1, rat1p, rat2p, rsc1p,
     5     fker, nd, ndd, dpars, ndz, zpars, ndi, ipars,
     6     nqorder, rfac, cintvals)
c
c
c  Compute the integrals in (1) defined at the top of the file
c  using a combination of distance based criterion
c  followed by adaptive integration.
c
c  A triangle is subdivided when xyzproxy associated
c  with the triangle is separated from the centroid
c  of the triangle by at least rfac*r_{t} where r_{t}
c  is the radius of the smallest sphere  
c  The refinement stops when the integral computed using 
c  the children of the triangle agree with the integral
c  on the triangle to a tolerance \eps.
c
c  Vectorized version of ctriaints_comb
c
c
c  Input arguments:
c    - eps: real *8
c        requested tolerance. Currently unsued, nqorder and rfac,
c        should be set based on eps.
c    - intype: integer *8
c        node type
c        * intype = 1, rokhlin vioreanu nodes
c        * intype = 2, xiao gimbutas nodes
c    - npatches: integer *8
c        number of patches
c    - norder: integer *8
c        order of discretization on patches
c    - npols: integer *8
c        number of points/basis functions on the 
c        patch = (norder+1)*(norder+2)/2
c    - isd: integer *8
c        flag for computing and passing second derivative
c        information in srcvals.
c        if isd = 0, leading dim of srcvals = 12, with 
c        xyz, dxyz/du, dxyz/dv, and normals passed
c        else, leading dim of srcvals = 30, with
c        xyz, dxyz/du, dxyz/dv, normals, d2xyz/du2, d2xyz/duv,
c        d2xyz/dv2, det(g), (11,12,22) entries of I, (11,12,22)
c        entries of II, kap1, kap2 where kap1, kap2
c        are principal curvatures, det(g) is the determinant
c        of the metric tensor
c    - ndsc: integer *8
c        leading dimension of srccoefs array, should be 9
c        if isd = 0, and 18 otherwise
c    - srccoefs: real *8 (ndsc, npols, npatches)
c        Koornwinder expansion coefficients of xyz, d/du (xyz), 
c        d/dv (xyz), (d2(xyz)/du2, d2(xyz)/duv, d2(xyz)/dv2)
c    - ndtarg: integer *8
c        leading dimension of target information
c    - ntarg: integer *8
c        number of targets
c    - xyztarg: real *8 (ndtarg, ntarg)
c        target information
c    - ifp: integer *8
c        flag for using proxy points for refining based 
c        on distance criterion
c    - xyzproxy: real *8 (3,*)
c        Proxy target points for refining triangles based
c        distance criterion, distance to proxy
c        point will be used instead of distance to target
c        should be of size (3,ntarg), if ifp=1
c    - itargptr: integer *8(npatches)
c        Pointer array for determining which targets are relevant
c        for patch i. 
c        xyztargs(:,itargptr(i):itargptr(i)+ntargptr(i)-1)
c        are the relevant set of targets for patch i
c    - ntargptr: integer *8(npatches)
c        number of targets relevant for patch i
c    - nporder: integer *8
c        order of Koornwinder polynomials to be integrated
c    - nppols: integer *8
c        number of polynomials corresponding to 
c        nporder = (nporder+1)*(nporder+2)/2
c    - ntrimax: integer *8
c        Initial guess for maximum number of triangles 
c        allowed in heirarchy of triangles. 
c        Recommended value 3000.
c    - rat1: real *8 (2, 0:norder)
c        recurrence coefficients for Legendre polynomials
c        with discretization order = norder
c    - rat2: real *8 (3, 0:norder, 0:norder)
c        recurrence coefficients for Jacobi polynomials
c        with discretization order = norder
c    - rsc1: real *8 (0:norder, 0:norder)
c        scaling parameters for the Koornwinder polynomials
c        with discretization order = norder
c    - rat1p: real *8 (2, 0:nporder)
c        recurrence coefficients for Legendre polynomials
c        with discretization order = nporder
c    - rat2p: real *8 (3, 0:nporder, 0:nporder)
c        recurrence coefficients for Jacobi polynomials
c        with discretization order = nporder
c    - rsc1p: real *8 (0:nporder, 0:nporder)
c        scaling parameters for the Koornwinder polynomials
c        with discretization order = nporder
c    - fker: function handle
c        function handle for evaluating the kernel k
c        * expected calling sequence
c            fker(nd,x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
c               
c        In this routine the output is expected to be a complex
c        vector
c    - nd: integer *8
c        number of kernels
c    - ndd: integer *8
c        number of real *8/double precision parameters
c    - dpars: real *8 (ndd)
c         real *8/ double precision paramters
c    - ndz: integer *8
c        number of complex *16 parameters
c    - zpars: complex *16(ndz)
c        complex *16 parameters
c    - ndi: integer *8
c        number of integer parameters
c    - ipars: integer *8(ndi)
c        integer parameters
c    - nqorder: integer *8
c        order of quadrature nodes to be used on each triangle
c    - rfac: real *8
c        scaling factor for determining when a triangle 
c        needs refinement
c 
c  Output arguments:
c    - cintvals: complex *16 (nppols, ntarg)
c        the computed integrals

      implicit none

c
cc     calling sequence variables
c
      real *8 eps
      integer *8 intype,ifp,nd
      integer *8 npatches,norder,npols
      integer *8 nporder,nppols
      integer *8 ndsc, isd
      real *8 srccoefs(ndsc,npols,npatches)
      
      integer *8 ntarg,ndtarg
      real *8 xyztarg(ndtarg,ntarg)
      real *8 xyzproxy(3,*)
      integer *8 itargptr(npatches)
      integer *8 ntargptr(npatches)
      
      external fker
      integer *8 ndd,ndz,ndi
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer *8 ipars(ndi)

      integer *8 nqorder
      real *8 rfac

      real *8 rat1(2,0:norder),rat2(3,0:norder,0:norder)
      real *8 rsc1(0:norder,0:norder)

      real *8 rat1p(2,0:nporder),rat2p(3,0:nporder,0:nporder)
      real *8 rsc1p(0:nporder,0:nporder)

      complex *16 cintvals(nd,nppols,ntarg)
      

c
cc      temporary variables
c
      real *8, allocatable :: tricm(:,:,:)
      real *8, allocatable :: trirad(:,:)
      integer *8, allocatable :: itrireltmp(:,:)
      integer *8, allocatable :: itrirel(:,:)
      integer *8, allocatable :: itrirelall(:)
      real *8, allocatable :: tverts(:,:,:)
      integer *8, allocatable :: ichild_start(:)
      integer *8, allocatable :: ichild_start0(:)

      real *8, allocatable :: xyztargtmp(:,:)

      integer *8 ntri,ntrimax,nlev,itri,istart,i,j
      integer *8 ier,itarg,jj,jstart,nlmax,npts
      integer *8 iqtri,ii

      integer *8 npmax
      integer *8 idim

      real *8, allocatable :: uvsq(:,:),wts(:),uvtmp(:,:)
      real *8, allocatable :: umattmp(:,:),vmattmp(:,:)
      real *8, allocatable :: da(:)
      integer *8 nqpols
      real *8, allocatable :: sigvals(:,:)
      real *8, allocatable :: sigvalsdens(:,:)
      real *8, allocatable :: srcvals(:,:)
      real *8, allocatable :: qwts(:)
      complex *16, allocatable :: sigmatmp(:,:)
      real *8 xyztmp(3)
      integer *8 itmp

      complex *16, allocatable :: xkernvals(:,:)
      integer *8, allocatable :: istack(:)

      complex *16, allocatable :: cvals(:,:,:)
      
      integer *8 nproclist0,nproclist,ntri0,npts0
      complex *16 zz

      
      complex *16 ima

      character *1 transa,transb
      
      real *8 alpha,beta
      integer *8 lda,ldb,ldc,nds

      data ima/(0.0d0,1.0d0)/

      allocate(cvals(nd,nppols,ntrimax))
      allocate(istack(2*ntrimax))





cc      max number of levels
c
      nlmax = 20
      allocate(tricm(3,npatches,ntrimax),
     1   trirad(npatches,ntrimax),tverts(2,3,ntrimax))
      allocate(itrireltmp(ntarg,ntrimax))
      allocate(da(ntrimax),ichild_start(ntrimax))
      allocate(ichild_start0(ntrimax))

      do i=1,ntrimax
        do j=1,ntarg
          itrireltmp(j,i) = 0
        enddo
        
      enddo

      
      ntri = 0
      nlev = 0
      ier = 0

      do i=1,ntrimax
        ichild_start(i) = -1
      enddo

      

      if(ifp.eq.1) then
        call gettritree(npatches, norder, npols, ndsc, srccoefs, ntarg,
     1    xyzproxy, itargptr, ntargptr, ntrimax, nlmax, rfac, ntri,
     2    nlev, ichild_start, da, tricm, trirad, tverts, itrireltmp,
     3    ier)
      else
        allocate(xyztargtmp(3,ntarg))
        do i=1,ntarg
          do j=1,3
            xyztargtmp(j,i) = xyztarg(j,i)
          enddo
        enddo
        call gettritree(npatches, norder, npols, ndsc, srccoefs, ntarg,
     1     xyztargtmp, itargptr, ntargptr, ntrimax, nlmax, rfac, ntri,
     2     nlev, ichild_start, da, tricm, trirad, tverts, itrireltmp,
     3     ier)
        deallocate(xyztargtmp)
      endif

      ntri0 = ntri

      do i=1,ntri0
        ichild_start0(i) = ichild_start(i)
      enddo


      

      if(ier.ne.0) then
        call prinf('Could not allocate triangle tree*',i,0)
        call prinf('Exiting without computing anything*',i,0)
        call prinf('Press 0 to continue*',i,0)
        call prinf('Press 1 to exit routine and continue*',i,0)
        call prinf('Press 2 to stop*',i,0)
        read *, ier
        if(ier.eq.2) stop
        if(ier.eq.1) return
      endif

      allocate(itrirel(ntri0,ntarg),itrirelall(ntri0))

c
c       transpose the itrirel array
c  

      do itri=1,ntri0
        do j=1,ntarg
          itrirel(itri,j) = itrireltmp(j,itri)
        enddo
      enddo


c
c       compute all the geometry info and koornwinder polynomials
c       at relevant triangle given by itrirelall
c  

c
c       get quadrature nodes and weights on the base triangle
c       based on quadrature type
c

      if(intype.eq.1) then
         nqpols = (nqorder+1)*(nqorder+2)/2
         allocate(uvsq(2,nqpols),wts(nqpols))
         allocate(umattmp(nqpols,nqpols),vmattmp(nqpols,nqpols))
      
         call vioreanu_simplex_quad(nqorder, nqpols, uvsq, umattmp,
     1      vmattmp, wts)
         
         deallocate(umattmp,vmattmp)
      endif


      if(intype.eq.2) then
        call triasymq_pts(nqorder,nqpols)
        allocate(uvsq(2,nqpols),wts(nqpols))
        call triasymq(nqorder, tverts(1,1,1), tverts(1,2,1), 
     1   tverts(1,3,1), uvsq, wts, nqpols)    
      endif

cc      call prinf('nqpols=*',nqpols,1)

      npts0 = ntri0*nqpols
      npmax = ntrimax*nqpols
      allocate(sigvals(npols,npmax))
      allocate(sigvalsdens(nppols,npmax))
      nds = 12
      if(isd.eq.0) nds = 30
      allocate(srcvals(nds,npmax),qwts(npmax))

      allocate(xkernvals(nd,npmax))

      allocate(uvtmp(2,nqpols))


      do itri=1,ntri0
        call mapuv_tri(tverts(1,1,itri), nqpols, uvsq, uvtmp)
        istart = (itri-1)*nqpols
        do i=1,nqpols
          ii = istart+i
          call koornf_pols(uvtmp(1,i), norder, npols, sigvals(1,ii),
     1           rat1, rat2, rsc1)
          call koornf_pols(uvtmp(1,i), nporder, nppols,
     1    sigvalsdens(1,ii), rat1p, rat2p, rsc1p)
        enddo
      enddo


      do itri=1,npatches
        ntri = ntri0
c
cc       for the current patch compute all the geometry info
c

        transa = 'N'
        transb = 'N'
        alpha = 1.0d0
        beta = 0.0d0
        lda = ndsc
        ldb = npols
        ldc = nds

        call dgemm_guru(transa, transb, ndsc, npts0, npols, alpha,
     1     srccoefs(1,1,itri), lda, sigvals, ldb, beta, srcvals, ldc)



c
c        compute all the quadrature weights
c

        do i=1,ntri
          istart = (i-1)*nqpols+1
          call get_srcvals_auxinfo_tri(nqpols, wts, isd, 
     1      nds, srcvals(1,istart), da(i), qwts(istart))
        enddo

        do itarg=itargptr(itri),itargptr(itri)+ntargptr(itri)-1
          do i=1,nppols
            do idim=1,nd
              cintvals(idim,i,itarg) = 0
            enddo
          enddo
c
c           extract info of geometry on subtriangles relevant
c           for this computation
c
          nproclist0 = 0
          do iqtri=1,ntri0
            jstart = (iqtri-1)*nqpols
            if(itrirel(iqtri,itarg).eq.1) then
              nproclist0 = nproclist0 + 1
              istack(nproclist0) = iqtri
c
c
c              compute the integral due to itri
c

              do i=1,nppols
                do idim=1,nd
                  cvals(idim,i,iqtri) = 0
                enddo
              enddo
              do i=1,nqpols
                jj = jstart+i
                call fker(nd,srcvals(1,jj), ndtarg, xyztarg(1,itarg),
     1              ndd, dpars, ndz, zpars, ndi, ipars, 
     2              xkernvals(1,jj))
                do idim=1,nd
                  xkernvals(idim,jj) = xkernvals(idim,jj)*qwts(jj)
                enddo
                
                do j=1,nppols
                  do idim=1,nd
                    zz = xkernvals(idim,jj)*sigvalsdens(j,jj)
                    cvals(idim,j,iqtri) = cvals(idim,j,iqtri)+zz
                    cintvals(idim,j,itarg)=cintvals(idim,j,itarg)+zz
                  enddo
                enddo
              enddo
            endif
          enddo
c
c           done with initial computation of 
c


          call triaadap_main_vec(eps, nqpols, nlmax, ntrimax, ntri,
     1      ichild_start, tverts, da, uvsq, wts, norder, npols, 
     2      isd, ndsc, srccoefs(1,1,itri), npmax, nds, srcvals, qwts, 
     2      sigvals,
     3      nporder, nppols, sigvalsdens, ndtarg, xyztarg(1,itarg),
     3      rat1, rat2, rsc1, rat1p, rat2p, rsc1p,
     3      fker, nd, ndd, dpars, ndz, zpars, ndi, ipars, cvals, 
     3      istack, nproclist0, xkernvals, cintvals(1,1,itarg), ier)

      
        enddo

        do i=1,ntri0
          ichild_start(i) = ichild_start0(i)
        enddo
        do i=ntri0+1,ntri
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
