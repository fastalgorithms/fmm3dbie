c
c
c        routines in this file
c         
c           dtriaints_dist - refine until targets/proxy separated
c                      by triangle by rfac*r_{t}
c                      where r_{t} is the radius of enclosing
c                      sphere for triangle
c 
c           dtriaints_adap - completely adaptive integration
c                       for the triangle
c           
c           dtriaints_comb - combination of adaptive integration
c                       and distance based refinement.
c                       Use distance based refinement to compute
c                       the integrals on some coarse hierarchy
c                       of meshes and then use 
c                       adaptive integration from that
c                       point on
c
c            dtriaints_wnodes - compute integral using 
c              prescribed nodes and weights
c                       
c
c
c
c
c
      subroutine dtriaints_wnodes(npatches,norder,npols,
     1   srccoefs,ndtarg,ntarg,xyztarg,itargptr,ntargptr,
     2   nporder,nppols,fker,ndd,dpars,ndz,zpars,ndi,ipars,nqpts,
     3   qnodes,wts,cintvals)
c
c
c       this subroutine computes the integral of
c       a kernel against koornwinder polynomials
c       given a prescribed set of nodes and weights
c

      implicit none
      integer *8 npatches,norder,npols,ndtarg
      integer *8 nporder,nppols
      real *8 srccoefs(9,npols,npatches),xyztarg(ndtarg,ntarg)
      integer *8 ntarg
      integer *8 itargptr(npatches),ntargptr(npatches)
      integer *8 ndd,ndz,ndi
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer *8 ipars(ndi),nqpts
      real *8 qnodes(2,nqpts),wts(nqpts)

      real *8 cintvals(nppols,ntarg)
      real *8, allocatable :: rat1(:,:), rat2(:,:,:), rsc1(:,:)
      real *8, allocatable :: rat1p(:,:), rat2p(:,:,:), rsc1p(:,:)
      real *8, allocatable :: srcvals(:,:),qwts(:)
      real *8, allocatable :: rsigvals(:,:),rsigtmp(:)
      real *8, allocatable :: sigvals(:,:)
      real *8, allocatable :: xkernvals(:,:)

      integer *8 i,ipatch,j,lda,ldb,itarg,ldc,ntarg0,ii

      real *8 fval
      real *8 da,ra

      character *1 transa,transb
      real *8 alpha,beta
      real *8 alpha_c,beta_c

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


c
c    
c
c
      allocate(sigvals(nppols,nqpts),rsigvals(npols,nqpts))

      allocate(srcvals(12,nqpts),qwts(nqpts))

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
        lda = 9
        ldb = npols
        ldc = 12

        da = 1

        call dgemm_guru(transa,transb,int(9,8),nqpts,npols,alpha,
     1     srccoefs(1,1,ipatch),lda,rsigvals,ldb,beta,srcvals,ldc)

   
        call get_norms_qwts_tri(nqpts,wts,srcvals,
     1        da,qwts)
c
c
c          compute the kernel values for all the targets
c
        
        ntarg0 = ntargptr(ipatch)
        allocate(xkernvals(nqpts,ntarg0))
        do itarg=itargptr(ipatch),itargptr(ipatch)+ntarg0-1
          ii = itarg - itargptr(ipatch)+1
          do j=1,nqpts
            call fker(srcvals(1,j),ndtarg,xyztarg(1,itarg),ndd,dpars,
     1         ndz,zpars,ndi,ipars,fval)
            xkernvals(j,ii) = fval*qwts(j)
          enddo
        enddo


        transa = 'n'
        transb = 'n'
        alpha_c = 1
        beta_c = 0

      call dgemm_guru(transa,transb,nppols,ntarg0,nqpts,alpha_c,sigvals,
     1   nppols,xkernvals,nqpts,beta_c,cintvals(1,itargptr(ipatch)),
     2   nppols)
      
        deallocate(xkernvals)
c

      enddo

      

      return
      end


      subroutine dtriaints_dist(eps,intype,
     1     npatches,norder,npols,srccoefs,ndtarg,
     2     ntarg,xyztarg,ifp,xyzproxy,
     3     itargptr,ntargptr,nporder,nppols,ntrimax,rat1,rat2,rsc1,
     3     rat1p,rat2p,rsc1p,fker,ndd,dpars,ndz,zpars,ndi,ipars,
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
c        this code is the same as dtriaints on all counts except
c        for choosing when to subdivide a triangle
c        
c        here a triangle is subdivided when xyzproxy associated
c        with the triangle is separated from the centroid
c        of the triangle by at least rfac*r_{t} where r_{t}
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
c        xyztarg(ndtarg,ntarg) - 
c                       target vectors
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
c        nporder - order of koornwinder polynomials to be integrated
c        nppols - number of koornwinder polynomials to be integrated 
c                  (nppols = (nporder+1)*(nporder+2)/2)
c
c        ntrimax - max number of triangles to be used on base triangle
c        fker - function handle for evaluating the kernel k
c 
c               expected calling sequence
c               fker(x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
c               
c
c         ndd - number of real parameters
c         dpars(ndd) - real parameters for the fker routine
c         ndz - number of complex parameters
c         zpars(ndz) - complex parameters for the fker routine
c         ndi - number of integer *8 parameters
c         ipars(ndi) - integer *8 parameters for the fker routine
c         nqorder - order of quadrature nodes on each subtriangle
c                   to be used
c         rfac - distance criterion for deciding whether a triangle
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
      integer *8 intype,ifp
      integer *8 npatches,norder,npols
      integer *8 nporder,nppols
      real *8 srccoefs(9,npols,npatches)
      
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

      real *8 cintvals(nppols,ntarg)
      

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
      real *8, allocatable :: fkervals(:),sigmatmp(:,:)
      real *8, allocatable :: rsigtmp(:,:)
      real *8 xyztmp(3)
      integer *8 itmp
      real *8 alpha_c, beta_c

      character *1 transa,transb
      double precision :: alpha,beta
      integer *8 lda,ldb,ldc
      





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
        call gettritree(npatches,norder,npols,srccoefs,ntarg,
     1     xyzproxy,
     1     itargptr,ntargptr,ntrimax,nlmax,rfac,ntri,nlev,ichild_start,
     2     da,tricm,trirad,tverts,itrireltmp,ier)
      else
        allocate(xyztargtmp(3,ntarg))

        do i=1,ntarg
          do j=1,3
            xyztargtmp(j,i) = xyztarg(j,i)
          enddo
        enddo
        call gettritree(npatches,norder,npols,srccoefs,ntarg,
     1     xyztargtmp,
     1     itargptr,ntargptr,ntrimax,nlmax,rfac,ntri,nlev,ichild_start,
     2     da,tricm,trirad,tverts,itrireltmp,ier)
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
      
         call vioreanu_simplex_quad(nqorder,nqpols,uvsq,umattmp,
     1      vmattmp,wts)
         
         deallocate(umattmp,vmattmp)
      endif


      if(intype.eq.2) then
        call triasymq_pts(nqorder,nqpols)
        allocate(uvsq(2,nqpols),wts(nqpols))
        call triasymq(nqorder,tverts(1,1,1),tverts(1,2,1),tverts(1,3,1),
     1         uvsq,wts,nqpols)    
      endif

cc      call prinf('nqpols=*',nqpols,1)

      npmax = ntri*nqpols
      allocate(sigvals(npols,npmax))
      allocate(srcvals(12,npmax),qwts(npmax))

      allocate(sigmatmp(nppols,npmax),fkervals(npmax))
      allocate(rsigtmp(nppols,npmax))

      allocate(uvtmp(2,nqpols))
      alpha_c = 1
      beta_c = 0


      do itri=1,ntri
        call mapuv_tri(tverts(1,1,itri),nqpols,uvsq,uvtmp)
        istart = (itri-1)*nqpols
        do i=1,nqpols
          ii = istart+i
          call koornf_pols(uvtmp(1,i),norder,npols,sigvals(1,ii),
     1        rat1,rat2,rsc1)
          call koornf_pols(uvtmp(1,i),nporder,nppols,rsigtmp(1,ii),
     1        rat1p,rat2p,rsc1p)
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
        lda = 9
        ldb = npols
        ldc = 12

        call dgemm_guru(transa,transb,int(9,8),npmax,npols,alpha,
     1     srccoefs(1,1,itri),lda,sigvals,ldb,beta,srcvals,ldc)


c
c        compute all the quadrature weights
c

        do i=1,ntri
          istart = (i-1)*nqpols+1
          call get_norms_qwts_tri(nqpols,wts,srcvals(1,istart),
     1        da(i),qwts(istart))
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

          call dgemv_guru('n',nppols,npts,alpha_c,sigmatmp,nppols,
     1      fkervals,int(1,8),beta_c,cintvals(1,itarg),int(1,8))
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
      subroutine dtriaints_adap(eps,intype,
     1     npatches,norder,npols,srccoefs,ndtarg,
     2     ntarg,xyztarg,itargptr,ntargptr,nporder,nppols,
     3     ntrimax,rat1,rat2,rsc1,rat1p,rat2p,rsc1p,
     4     fker,ndd,dpars,ndz,zpars,ndi,ipars,nqorder,cintvals)

c
      implicit none

c
cc     calling sequence variables
c
      real *8 eps
      integer *8 intype
      integer *8 npatches,norder,npols
      integer *8 nporder,nppols
      real *8 srccoefs(9,npols,npatches)
      
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

      real *8 cintvals(nppols,ntarg)

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
      integer *8 lda,ldb,ldc
      integer *8 nn1,nn2,nn3,nn4,npmax0,ntmaxuse,ntmaxuse0
      
      
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
      
         call vioreanu_simplex_quad(nqorder,nqpols,uvsq,umattmp,
     1      vmattmp,wts)
         
         deallocate(umattmp,vmattmp)
      endif

      if(intype.eq.2) then
        call triasymq_pts(nqorder,nqpols)
        allocate(uvsq(2,nqpols),wts(nqpols))
        call triasymq(nqorder,tvs(1,1,1),tvs(1,2,1),tvs(1,3,1),
     1         uvsq,wts,nqpols)    
      endif

      allocate(uvtmp(2,nqpols))
     
      npmax = ntrimax*nqpols
      allocate(sigvals(npols,npmax))
      allocate(sigvalsdens(nppols,npmax))
      allocate(srcvals(12,npmax),qwts(npmax))

c
c      current number of triangles in the adaptive structure
c
      ntri = 1
c
c        intialize sigvals for root triangle
c


      call mapuv_tri(tvs(1,1,1),nqpols,uvsq,uvtmp)
      do i=1,nqpols
        call koornf_pols(uvtmp(1,i),norder,npols,
     1      sigvals(1,i),rat1,rat2,rsc1)
        call koornf_pols(uvtmp(1,i),nporder,nppols,
     1      sigvalsdens(1,i),rat1p,rat2p,rsc1p)
      enddo


      do itri=1,npatches
c
cc       for the current patch compute geometry info for base triangle 
c
        transa = 'N'
        transb = 'N'
        alpha = 1
        beta = 0
        lda = 9
        ldb = npols
        ldc = 12

        call dgemm_guru(transa,transb,int(9,8),nqpols,npols,alpha,
     1     srccoefs(1,1,itri),lda,sigvals,ldb,beta,srcvals,ldc)
        call get_norms_qwts_tri(nqpols,wts,srcvals,da,qwts)

        do itarg=itargptr(itri),itargptr(itri)+ntargptr(itri)-1


 1111     continue
          call dtriaadap(eps,nqorder,nqpols,nlmax,ntmaxuse,ntri,
     1          ichild_start,tvs,da,uvsq,wts, 
     1          norder,npols,srccoefs(1,1,itri),
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
             nn4 = ntri*6
             call dcopy_guru(nn1,sigvals,1,sigvals2,1)
             call dcopy_guru(nn2,sigvalsdens,1,sigvalsdens2,1)
             call dcopy_guru(nn3,srcvals,1,srcvals2,1)
             call dcopy_guru(npmax,qwts,1,qwts2,1)
             do ii=1,ntri
               ichild_start2(ii) = ichild_start(ii)
             enddo
             call dcopy_guru(nn4,tvs,1,tvs2,1)
             call dcopy_guru(ntri,da,1,da2,1)


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

             call dcopy_guru(nn1,sigvals2,1,sigvals,1)
             call dcopy_guru(nn2,sigvalsdens2,1,sigvalsdens,1)
             call dcopy_guru(nn3,srcvals2,1,srcvals,1)
             call dcopy_guru(npmax,qwts2,1,qwts,1)
             do ii=1,ntri
               ichild_start(ii) = ichild_start2(ii)
             enddo
             call dcopy_guru(nn4,tvs2,1,tvs,1)
             call dcopy_guru(ntri,da2,1,da,1)

             npmax = npmax0
             ntmaxuse = ntmaxuse0
c             print *, "restrating adaptive integration with ntri=",
c     1          ntmaxuse


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
      subroutine dtriaadap(eps,m,kpols,nlmax,ntmax,ntri,
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
c        kpols - number of quarature nodes (look up triasymq 
c                 for getting kpols(m))
c        nlmax - max level of refinement for geometry
c        ntmax - max number of triangles
c        ntri - current number of triangles in adaptive structure
c        ichild_start(i) - first child of triangle i
c        tvs(2,3,ntmax) - vertices of hierarchy of triangles
c        da(ntmax) - area of triangles
c        uvsq(kpols) - integration nodes on standard triangle
c        wts(kpols) - integration weights on standard triangle
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
c
c         

      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, allocatable :: istack(:)
      integer *8 ichild_start(ntmax)
      real *8 da(ntmax)
      real *8 tvs(2,3,ntmax), uvsq(2,kpols),wts(kpols)
      integer *8 nproclist0, nproclist
      integer *8 idone
      real *8 srccoefs(9,npols)
      real *8 sigvals(npols,npmax)
      real *8 sigvalsdens(nppols,npmax)
      real *8 srcvals(12,*),qwts(npmax)
      real *8, allocatable :: xkernvals(:)
      real *8 xt(ndtarg)
      real *8 cintall(nppols),fval
      real *8, allocatable :: ctmp(:)
      real *8, allocatable :: cvals(:,:)

      integer *8 ndd,ndz,ndi
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer *8 ipars(ndi)

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


      call dtriaadap_main(eps,kpols,nlmax,ntmax,ntri,ichild_start,
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
       
      subroutine dtriaadap_main(eps,kpols,nlmax,ntmax,ntri,ichild_start,
     1      tvs,da,uvsq,wts,norder,npols,srccoefs,
     2      npmax,srcvals,qwts,sigvals,nporder,nppols,sigvalsdens,
     3      ndtarg,xt,rat1,rat2,rsc1,
     3      rat1p,rat2p,rsc1p,
     3      fker,ndd,dpars,ndz,
     3      zpars,ndi,ipars,cvals,istack,nproclist0,xkernvals,
     4      cintall,ier)
      

      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8 istack(*),nproclist0
      integer *8 ichild_start(ntmax)
      integer *8 nporder,nppols
      real *8 da(ntmax)
      real *8 tvs(2,3,ntmax), uvsq(2,kpols),wts(kpols)
      integer *8  nproclist
      integer *8 idone
      real *8 srccoefs(9,npols)
      real *8 sigvals(npols,npmax)
      real *8 sigvalsdens(nppols,npmax)
      real *8 qwts(npmax)
      real *8 xkernvals(npmax)
      real *8 xt(ndtarg)
      real *8 srcvals(12,*)
      real *8 cintall(nppols),fval,ctmp(nppols)
      real *8 cvals(nppols,ntmax)

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
      
      ier = 0
      allocate(uvtmp(2,kpols))

c
c         for historic reasons
c
      ksigpols = nppols
      kfine = 4*kpols


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
c               print *, "Too many triangles in dtriaadap"
c               print *, "Exiting without computing anything"
               ier = 4

               return
            endif
            
            ichild_start(itri) = ntri+1
            call gettrichildren(tvs(1,1,itri),tvs(1,1,ntri+1),
     1             tvs(1,1,ntri+2),tvs(1,1,ntri+3),tvs(1,1,ntri+4))

            rr = 0.25d0*da(itri)
            do j=ntri+1,ntri+4
              da(j) = rr
              call mapuv_tri(tvs(1,1,j),kpols,uvsq,uvtmp)
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

              call dgemm_guru(transa,transb,int(9,8),kpols,npols,alpha,
     1           srccoefs,lda,sigvals(1,istart),ldb,beta,
     2           srcvals(1,istart),ldc)
              call get_norms_qwts_tri(kpols,wts,srcvals(1,istart),
     1           rr,qwts(istart))

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
            call fker(srcvals(1,jj),ndtarg,xt,ndd,dpars,
     1         ndz,zpars,ndi,ipars,fval)
            xkernvals(jj) = fval*qwts(jj)
          enddo
cc          call prin2('xkernvals=*',xkernvals(istart+1),kfine)

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
c
c
      subroutine dtriaints_comb(eps,intype,
     1     npatches,norder,npols,srccoefs,ndtarg,
     2     ntarg,xyztarg,ifp,xyzproxy,
     3     itargptr,ntargptr,nporder,nppols,ntrimax,rat1,rat2,rsc1,
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
c        this code is the same as dtriaints on all counts except
c        for choosing when to subdivide a triangle
c        
c        here a triangle is subdivided when xyzproxy associated
c        with the triangle is separated from the centroid
c        of the triangle by at least rfac*r_{t} where r_{t}
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
c        ntrimax - max number of triangles to be used on base triangle
c        fker - function handle for evaluating the kernel k
c 
c               
c               expected calling sequence
c               fker(x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
c
c         dpars(*) - real parameters for the fker routine
c         zpars(*) - complex parameters for the fker routine
c         ipars(*) - integer *8 parameters for the fker routine
c         nqorder - order of quadrature nodes on each subtriangle
c                   to be used
c         rfac - distance criterion for deciding whether a triangle
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
      integer *8 intype,ifp
      integer *8 npatches,norder,npols
      integer *8 nporder,nppols
      real *8 srccoefs(9,npols,npatches)
      
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

      real *8 cintvals(nppols,ntarg)
      

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
      real *8, allocatable :: sigmatmp(:,:)
      real *8, allocatable :: xyztargtmp(:,:)
      real *8 xyztmp(3)
      integer *8 itmp

      real *8, allocatable :: xkernvals(:)
      integer *8, allocatable :: istack(:)

      real *8, allocatable :: cvals(:,:)
      
      integer *8 nproclist0,nproclist,ntri0,npts0
      real *8 zz

      

      character *1 transa,transb
      
      real *8 alpha,beta
      integer *8 lda,ldb,ldc

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
        call gettritree(npatches,norder,npols,srccoefs,ntarg,xyzproxy,
     1     itargptr,ntargptr,ntrimax,nlmax,rfac,ntri,nlev,ichild_start,
     2     da,tricm,trirad,tverts,itrireltmp,ier)
      else
        allocate(xyztargtmp(3,ntarg))
        do i=1,ntarg
          do j=1,3
            xyztargtmp(j,i) = xyztarg(j,i)
          enddo
        enddo
        call gettritree(npatches,norder,npols,srccoefs,ntarg,
     1     xyztargtmp,
     1     itargptr,ntargptr,ntrimax,nlmax,rfac,ntri,nlev,ichild_start,
     2     da,tricm,trirad,tverts,itrireltmp,ier)
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
      
         call vioreanu_simplex_quad(nqorder,nqpols,uvsq,umattmp,
     1      vmattmp,wts)
         
         deallocate(umattmp,vmattmp)
      endif


      if(intype.eq.2) then
        call triasymq_pts(nqorder,nqpols)
        allocate(uvsq(2,nqpols),wts(nqpols))
        call triasymq(nqorder,tverts(1,1,1),tverts(1,2,1),tverts(1,3,1),
     1         uvsq,wts,nqpols)    
      endif

cc      call prinf('nqpols=*',nqpols,1)

      npts0 = ntri0*nqpols
      npmax = ntrimax*nqpols
      allocate(sigvals(npols,npmax))
      allocate(sigvalsdens(nppols,npmax))
      allocate(srcvals(12,npmax),qwts(npmax))

      allocate(xkernvals(npmax))

      allocate(uvtmp(2,nqpols))


      do itri=1,ntri0
        call mapuv_tri(tverts(1,1,itri),nqpols,uvsq,uvtmp)
        istart = (itri-1)*nqpols
        do i=1,nqpols
          ii = istart+i
          call koornf_pols(uvtmp(1,i),norder,npols,sigvals(1,ii),
     1           rat1,rat2,rsc1)
          call koornf_pols(uvtmp(1,i),nporder,nppols,sigvalsdens(1,ii),
     1           rat1p,rat2p,rsc1p)
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
        lda = 9
        ldb = npols
        ldc = 12

        call dgemm_guru(transa,transb,int(9,8),npts0,npols,alpha,
     1     srccoefs(1,1,itri),lda,sigvals,ldb,beta,srcvals,ldc)



c
c        compute all the quadrature weights
c

        do i=1,ntri
          istart = (i-1)*nqpols+1
          call get_norms_qwts_tri(nqpols,wts,srcvals(1,istart),
     1        da(i),qwts(istart))
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
                call fker(srcvals(1,jj),ndtarg,xyztarg(1,itarg),
     1              ndd,dpars,ndz,zpars,ndi,ipars,xkernvals(jj))
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

          call dtriaadap_main(eps,nqpols,nlmax,ntrimax,ntri,
     1      ichild_start,
     1      tverts,da,uvsq,wts,norder,npols,srccoefs(1,1,itri),
     2      npmax,srcvals,qwts,sigvals,nporder,nppols,sigvalsdens,
     3      ndtarg,xyztarg(1,itarg),
     3      rat1,rat2,rsc1,rat1p,rat2p,rsc1p,
     3      fker,ndd,dpars,ndz,zpars,ndi,ipars,cvals,istack,nproclist0,
     4      xkernvals,cintvals(1,itarg),ier)

      
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
c        vectorized routines
c ----------------------------
c         
c           dtriaints_dist_vec - refine until targets/proxy separated
c                      by triangle by rfac*r_{t}
c                      where r_{t} is the radius of enclosing
c                      sphere for triangle
c            vectorized verison of dtriaints_dist
c
c 
c           dtriaints_adap_vec - completely adaptive integration
c                       for the triangle
c             vectorized verion of dtriaints_adap
c
c           
c           dtriaints_comb_vec - combination of adaptive integration
c                       and distance based refinement.
c                       Use distance based refinement to compute
c                       the integrals on some coarse hierarchy
c                       of meshes and then use 
c                       adaptive integration from that
c                       point on
c               vectorized version of dtriaints_comb
c               
c                       
c
c
c
c
c
      subroutine dtriaints_dist_vec(eps,intype,
     1     npatches,norder,npols,srccoefs,ndtarg,
     2     ntarg,xyztarg,ifp,xyzproxy,
     3     itargptr,ntargptr,nporder,nppols,ntrimax,rat1,rat2,rsc1,
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
c        this code is the same as dtriaints on all counts except
c        for choosing when to subdivide a triangle
c        
c        here a triangle is subdivided when xyzproxy associated
c        with the triangle is separated from the centroid
c        of the triangle by at least rfac*r_{t} where r_{t}
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
c        ntrimax - max number of triangles to be used on base triangle
c        fker - function handle for evaluating the kernel k
c 
c               expected calling sequence
c               fker(x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
c               
c               the output is assumed to be nd real numbers
c
c         nd - number of kernels to be integrated
c         dpars(*) - real parameters for the fker routine
c         zpars(*) - complex parameters for the fker routine
c         ipars(*) - integer *8 parameters for the fker routine
c         nqorder - order of quadrature nodes on each subtriangle
c                   to be used
c         rfac - distance criterion for deciding whether a triangle
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
      integer *8 intype,ifp,nd
      integer *8 npatches,norder,npols
      integer *8 nporder,nppols
      real *8 srccoefs(9,npols,npatches)
      
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

      real *8 cintvals(nd,nppols,ntarg)
      

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
      real *8, allocatable :: fkervals(:,:),sigmatmp(:,:)
      real *8, allocatable :: rsigtmp(:,:)
      real *8 xyztmp(3)
      integer *8 itmp,idim

      character *1 transa,transb
      double precision :: alpha,beta
      real *8 alpha_c, beta_c
      integer *8 lda,ldb,ldc
      





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
        call gettritree(npatches,norder,npols,srccoefs,ntarg,xyzproxy,
     1     itargptr,ntargptr,ntrimax,nlmax,rfac,ntri,nlev,ichild_start,
     2     da,tricm,trirad,tverts,itrireltmp,ier)
      else
        allocate(xyztargtmp(3,ntarg))
        do i=1,ntarg
          do j=1,3
            xyztargtmp(j,i) = xyztarg(j,i)
          enddo
        enddo

        call gettritree(npatches,norder,npols,srccoefs,ntarg,xyztargtmp,
     1     itargptr,ntargptr,ntrimax,nlmax,rfac,ntri,nlev,ichild_start,
     2     da,tricm,trirad,tverts,itrireltmp,ier)
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
      
         call vioreanu_simplex_quad(nqorder,nqpols,uvsq,umattmp,
     1      vmattmp,wts)
         
         deallocate(umattmp,vmattmp)
      endif


      if(intype.eq.2) then
        call triasymq_pts(nqorder,nqpols)
        allocate(uvsq(2,nqpols),wts(nqpols))
        call triasymq(nqorder,tverts(1,1,1),tverts(1,2,1),tverts(1,3,1),
     1         uvsq,wts,nqpols)    
      endif

cc      call prinf('nqpols=*',nqpols,1)

      npmax = ntri*nqpols
      allocate(sigvals(npols,npmax))
      allocate(srcvals(12,npmax),qwts(npmax))

      allocate(sigmatmp(nppols,npmax),fkervals(nd,npmax))
      allocate(rsigtmp(nppols,npmax))

      allocate(uvtmp(2,nqpols))


      do itri=1,ntri
        call mapuv_tri(tverts(1,1,itri),nqpols,uvsq,uvtmp)
        istart = (itri-1)*nqpols
        do i=1,nqpols
          ii = istart+i
          call koornf_pols(uvtmp(1,i),norder,npols,sigvals(1,ii),
     1        rat1,rat2,rsc1)
          call koornf_pols(uvtmp(1,i),nporder,nppols,rsigtmp(1,ii),
     1        rat1p,rat2p,rsc1p)
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
        lda = 9
        ldb = npols
        ldc = 12

        call dgemm_guru(transa,transb,int(9,8),npmax,npols,alpha,
     1     srccoefs(1,1,itri),lda,sigvals,ldb,beta,srcvals,ldc)


c
c        compute all the quadrature weights
c

        do i=1,ntri
          istart = (i-1)*nqpols+1
          call get_norms_qwts_tri(nqpols,wts,srcvals(1,istart),
     1        da(i),qwts(istart))
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

c
c          TODO:  fix this to call mkl blas with single thread
c
          transa = 'n'
          transb = 't'
          alpha_c = 1
          beta_c = 0
          call dgemm_guru(transa,transb,nd,nppols,npts,alpha_c,
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
      subroutine dtriaints_adap_vec(eps,intype,
     1     npatches,norder,npols,srccoefs,ndtarg,
     2     ntarg,xyztarg,itargptr,ntargptr,nporder,nppols,ntrimax,
     3     rat1,rat2,rsc1,rat1p,rat2p,rsc1p,
     3     fker,nd,ndd,dpars,ndz,zpars,ndi,ipars,nqorder,cintvals)

c
      implicit none

c
cc     calling sequence variables
c
      real *8 eps
      integer *8 intype,nd
      integer *8 npatches,norder,npols
      integer *8 nporder,nppols
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

      real *8 cintvals(nd,nppols,ntarg)

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
      integer *8 lda,ldb,ldc
      integer *8 nn1,nn2,nn3,nn4,npmax0,ntmaxuse,ntmaxuse0
      
      
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
      
         call vioreanu_simplex_quad(nqorder,nqpols,uvsq,umattmp,
     1      vmattmp,wts)
         
         deallocate(umattmp,vmattmp)
      endif

      if(intype.eq.2) then
        call triasymq_pts(nqorder,nqpols)
        allocate(uvsq(2,nqpols),wts(nqpols))
        call triasymq(nqorder,tvs(1,1,1),tvs(1,2,1),tvs(1,3,1),
     1         uvsq,wts,nqpols)    
      endif

      allocate(uvtmp(2,nqpols))
     
      npmax = ntrimax*nqpols
      allocate(sigvals(npols,npmax),sigvalsdens(nppols,npmax))
      allocate(srcvals(12,npmax),qwts(npmax))


c
c      current number of triangles in the adaptive structure
c
      ntri = 1
c
c        intialize sigvals for root triangle
c


      call mapuv_tri(tvs(1,1,1),nqpols,uvsq,uvtmp)
      do i=1,nqpols
        call koornf_pols(uvtmp(1,i),norder,npols,
     1      sigvals(1,i),rat1,rat2,rsc1)
        call koornf_pols(uvtmp(1,i),nporder,nppols,
     1      sigvalsdens(1,i),rat1p,rat2p,rsc1p)
      enddo


      do itri=1,npatches
c
cc       for the current patch compute geometry info for base triangle 
c
        transa = 'N'
        transb = 'N'
        alpha = 1
        beta = 0
        lda = 9
        ldb = npols
        ldc = 12

        call dgemm_guru(transa,transb,int(9,8),nqpols,npols,alpha,
     1     srccoefs(1,1,itri),lda,sigvals,ldb,beta,srcvals,ldc)
        call get_norms_qwts_tri(nqpols,wts,srcvals,da,qwts)

        do itarg=itargptr(itri),itargptr(itri)+ntargptr(itri)-1

 1111     continue
          call dtriaadap_vec(eps,nqorder,nqpols,nlmax,ntmaxuse,ntri,
     1          ichild_start,tvs,da,uvsq,wts, 
     1          norder,npols,srccoefs(1,1,itri),
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
             nn4 = ntri*6
             call dcopy_guru(nn1,sigvals,1,sigvals2,1)
             call dcopy_guru(nn2,sigvalsdens,1,sigvalsdens2,1)
             call dcopy_guru(nn3,srcvals,1,srcvals2,1)
             call dcopy_guru(npmax,qwts,1,qwts2,1)
             do ii=1,ntri
               ichild_start2(ii) = ichild_start(ii)
             enddo
             call dcopy_guru(nn4,tvs,1,tvs2,1)
             call dcopy_guru(ntri,da,1,da2,1)


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

             call dcopy_guru(nn1,sigvals2,1,sigvals,1)
             call dcopy_guru(nn2,sigvalsdens2,1,sigvalsdens,1)
             call dcopy_guru(nn3,srcvals2,1,srcvals,1)
             call dcopy_guru(npmax,qwts2,1,qwts,1)
             do ii=1,ntri
               ichild_start(ii) = ichild_start2(ii)
             enddo
             call dcopy_guru(nn4,tvs2,1,tvs,1)
             call dcopy_guru(ntri,da2,1,da,1)

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
      subroutine dtriaadap_vec(eps,m,kpols,nlmax,ntmax,ntri,
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
c        kpols - number of quarature nodes (look up triasymq 
c                 for getting kpols(m))
c        nlmax - max level of refinement for geometry
c        ntmax - max number of triangles
c        ntri - current number of triangles in adaptive structure
c        ichild_start(i) - first child of triangle i
c        tvs(2,3,ntmax) - vertices of hierarchy of triangles
c        da(ntmax) - area of triangles
c        uvsq(kpols) - integration nodes on standard triangle
c        wts(kpols) - integration weights on standard triangle
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
      implicit integer *8 (i-n)
      integer *8, allocatable :: istack(:)
      integer *8 ichild_start(ntmax),nd
      integer *8 nporder,nppols
      real *8 da(ntmax)
      real *8 tvs(2,3,ntmax), uvsq(2,kpols),wts(kpols)
      integer *8 nproclist0, nproclist
      integer *8 idone
      real *8 srccoefs(9,npols)
      real *8 sigvals(npols,npmax)
      real *8 sigvalsdens(nppols,npmax)
      real *8 srcvals(12,*),qwts(npmax)
      real *8, allocatable :: xkernvals(:,:)
      real *8 xt(ndtarg)
      real *8 cintall(nd,npols),fval(nd)
      real *8, allocatable :: cvals(:,:,:)

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


      call dtriaadap_main_vec(eps,kpols,nlmax,ntmax,ntri,ichild_start,
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
       
      subroutine dtriaadap_main_vec(eps,kpols,nlmax,ntmax,ntri,
     1      ichild_start,tvs,da,uvsq,wts,norder,npols,srccoefs,
     2      npmax,srcvals,qwts,sigvals,nporder,nppols,sigvalsdens,
     3      ndtarg,xt,rat1,rat2,rsc1,rat1p,rat2p,rsc1p,
     3      fker,nd,ndd,dpars,ndz,
     3      zpars,ndi,ipars,cvals,istack,nproclist0,xkernvals,
     4      cintall,ier)
      

      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8 istack(*),nproclist0,nd
      integer *8 ichild_start(ntmax)
      integer *8 nporder,nppols
      real *8 da(ntmax)
      real *8 tvs(2,3,ntmax), uvsq(2,kpols),wts(kpols)
      integer *8  nproclist
      integer *8 idone
      real *8 srccoefs(9,npols)
      real *8 sigvals(npols,npmax)
      real *8 sigvalsdens(nppols,npmax)
      real *8 qwts(npmax)
      real *8 xkernvals(nd,npmax)
      real *8 xt(ndtarg)
      real *8 srcvals(12,*)
      real *8 cintall(nd,nppols),fval(nd),ctmp(nd,nppols)
      real *8 cvals(nd,nppols,ntmax)

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
c               print *, "Too many triangles in dtriaadap"
c               print *, "Exiting without computing anything"
               ier = 4
               return
            endif
            
            ichild_start(itri) = ntri+1
            call gettrichildren(tvs(1,1,itri),tvs(1,1,ntri+1),
     1             tvs(1,1,ntri+2),tvs(1,1,ntri+3),tvs(1,1,ntri+4))

            rr = 0.25d0*da(itri)
            do j=ntri+1,ntri+4
              da(j) = rr
              call mapuv_tri(tvs(1,1,j),kpols,uvsq,uvtmp)
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

              call dgemm_guru(transa,transb,int(9,8),kpols,npols,alpha,
     1           srccoefs,lda,sigvals(1,istart),ldb,beta,
     2           srcvals(1,istart),ldc)
              call get_norms_qwts_tri(kpols,wts,srcvals(1,istart),
     1           rr,qwts(istart))

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
            call fker(nd,srcvals(1,jj),ndtarg,xt,ndd,dpars,
     1         ndz,zpars,ndi,ipars,fval)
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
      subroutine dtriaints_comb_vec(eps,intype,
     1     npatches,norder,npols,srccoefs,ndtarg,
     2     ntarg,xyztarg,ifp,xyzproxy,
     3     itargptr,ntargptr,nporder,nppols,ntrimax,rat1,rat2,rsc1,
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
c        this code is the same as dtriaints on all counts except
c        for choosing when to subdivide a triangle
c        
c        here a triangle is subdivided when xyzproxy associated
c        with the triangle is separated from the centroid
c        of the triangle by at least rfac*r_{t} where r_{t}
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
c        ntrimax - max number of triangles to be used on base triangle
c        fker - function handle for evaluating the kernel k
c 
c               expected calling sequence
c               fker(x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
c               
c
c         nd - number of kernels to be integrated
c         dpars(*) - real parameters for the fker routine
c         zpars(*) - complex parameters for the fker routine
c         ipars(*) - integer *8 parameters for the fker routine
c         nqorder - order of quadrature nodes on each subtriangle
c                   to be used
c         rfac - distance criterion for deciding whether a triangle
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
      integer *8 intype,ifp,nd
      integer *8 npatches,norder,npols
      integer *8 nporder,nppols
      real *8 srccoefs(9,npols,npatches)
      
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

      real *8 cintvals(nd,nppols,ntarg)
      

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
      real *8, allocatable :: sigmatmp(:,:)
      real *8 xyztmp(3)
      integer *8 itmp

      real *8, allocatable :: xkernvals(:,:)
      integer *8, allocatable :: istack(:)

      real *8, allocatable :: cvals(:,:,:)
      
      integer *8 nproclist0,nproclist,ntri0,npts0
      real *8 zz

      

      character *1 transa,transb
      
      real *8 alpha,beta
      integer *8 lda,ldb,ldc

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
        call gettritree(npatches,norder,npols,srccoefs,ntarg,xyzproxy,
     1     itargptr,ntargptr,ntrimax,nlmax,rfac,ntri,nlev,ichild_start,
     2     da,tricm,trirad,tverts,itrireltmp,ier)
      else
        allocate(xyztargtmp(3,ntarg))
        do i=1,ntarg
          do j=1,3
            xyztargtmp(j,i) = xyztarg(j,i)
          enddo
        enddo
        call gettritree(npatches,norder,npols,srccoefs,ntarg,
     1     xyztargtmp,
     1     itargptr,ntargptr,ntrimax,nlmax,rfac,ntri,nlev,ichild_start,
     2     da,tricm,trirad,tverts,itrireltmp,ier)
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
      
         call vioreanu_simplex_quad(nqorder,nqpols,uvsq,umattmp,
     1      vmattmp,wts)
         
         deallocate(umattmp,vmattmp)
      endif


      if(intype.eq.2) then
        call triasymq_pts(nqorder,nqpols)
        allocate(uvsq(2,nqpols),wts(nqpols))
        call triasymq(nqorder,tverts(1,1,1),tverts(1,2,1),tverts(1,3,1),
     1         uvsq,wts,nqpols)    
      endif

cc      call prinf('nqpols=*',nqpols,1)

      npts0 = ntri0*nqpols
      npmax = ntrimax*nqpols
      allocate(sigvals(npols,npmax))
      allocate(sigvalsdens(nppols,npmax))
      allocate(srcvals(12,npmax),qwts(npmax))

      allocate(xkernvals(nd,npmax))

      allocate(uvtmp(2,nqpols))


      do itri=1,ntri0
        call mapuv_tri(tverts(1,1,itri),nqpols,uvsq,uvtmp)
        istart = (itri-1)*nqpols
        do i=1,nqpols
          ii = istart+i
          call koornf_pols(uvtmp(1,i),norder,npols,sigvals(1,ii),
     1           rat1,rat2,rsc1)
          call koornf_pols(uvtmp(1,i),nporder,nppols,sigvalsdens(1,ii),
     1           rat1p,rat2p,rsc1p)
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
        lda = 9
        ldb = npols
        ldc = 12

        call dgemm_guru(transa,transb,int(9,8),npts0,npols,alpha,
     1     srccoefs(1,1,itri),lda,sigvals,ldb,beta,srcvals,ldc)



c
c        compute all the quadrature weights
c

        do i=1,ntri
          istart = (i-1)*nqpols+1
          call get_norms_qwts_tri(nqpols,wts,srcvals(1,istart),
     1        da(i),qwts(istart))
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
                call fker(nd,srcvals(1,jj),ndtarg,xyztarg(1,itarg),
     1              ndd,dpars,ndz,zpars,ndi,ipars,xkernvals(1,jj))
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


          call dtriaadap_main_vec(eps,nqpols,nlmax,ntrimax,ntri,
     1      ichild_start,
     1      tverts,da,uvsq,wts,norder,npols,srccoefs(1,1,itri),
     2      npmax,srcvals,qwts,sigvals,nporder,nppols,sigvalsdens,
     3      ndtarg,xyztarg(1,itarg),
     3      rat1,rat2,rsc1,rat1p,rat2p,rsc1p,
     3      fker,nd,ndd,dpars,ndz,zpars,ndi,ipars,cvals,istack,
     3      nproclist0,
     4      xkernvals,cintvals(1,1,itarg),ier)

      
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
