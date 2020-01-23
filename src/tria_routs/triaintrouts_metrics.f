
c
c
c
c
      subroutine ctriaints2_withmetrics(eps,
     1     npatches,norder,npols,xyzcoefs,
     2     dxyzcoefs,hxyzcoefs,ntarg,xyztarg,xyzproxy,
     3     itargptr,ntargptr,fker,dpars,zpars,ipars,nqorder,rfac,
     4     cintvals,rn1,n2)
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
c        this code is the same as ctriaints on all counts except
c        for choosing when to subdivide a triangle
c        
c        here a triangle is subdivided when xyzproxy associated
c        with the triangle is separated from the centroid
c        of the triangle by at least rfac*r_{t} where r_{t}
c        is the radius of the smallest sphere  (see point 2 of
c        ctriaints)
c
c        input arguments:
c        eps:     requested precision (Currently not being used
c                 in any way, later, nqorder, and rfac will
c                 be set internally based on eps)
c
c        npatches: number of patches
c        norder: order of discretization nodes on the patches
c        npols = (norder+1)*(norder+2)/2: number of discretization 
c                   nodes on each patch
c        xyzcoefs(3,npols,npatches): coefficients
c                 of koornwinder expansion of xyz coordinates
c        dxyzcoefs(6,npols,npatches):
c                 coefficients of koornwinder expansion 
c                 of dxyz/duv
c        hxyzcoefs(9,npols,npatches):
c                 coefficients of koornwinder expansion
c                 of d2xyz/d2uv, 
c                 1 = d2x/duu
c                 2 = d2x/duv
c                 3 = d2x/dvv
c                 4 = d2y/duu
c                      .
c                      .
c                 9 = d2z/dvv
c
c        ntarg - total number of target points
c        xyztarg(3,ntarg) - location of target points 
c                       (currently may have repeats, 
c                        being slightly wasteful in
c                        memory, but see comment (1))
c
c        xyzproxy(3,ntarg) - proxy locations for measuring distances
c
c        itargptr(npatches) - xyztargs(:,itargptr(i):itargptr(i)+ntargptr(i)-1)
c                       are the relevant set of targets for patch i
c        ntargptr(npatches) - number of relevant targets for patch i
c        fker - function handle for evaluating the kernel k
c 
c               expected calling sequence
c               fker(x,y,dpars,zpars,ipars,f)
c               
c               the output is assumed to be complex for the time
c               being
c
c         dpars(*) - real parameters for the fker routine
c         zpars(*) - complex parameters for the fker routine
c         ipars(*) - integer parameters for the fker routine
c         nqorder - order of quadrature nodes on each subtriangle
c                   to be used
c         rfac - distance criterion for deciding whether a triangle
c                is in the far-field or not (See comment (2))
c
c         output:
c
c         cintvals(npols,ntarg) - integrals at all targets
c                                  for all koornwinder
c                                  polynomials
c
      implicit none

c
cc     calling sequence variables
c
      real *8 eps
      integer npatches,norder,npols
      real *8 xyzcoefs(3,npols,npatches)
      real *8 dxyzcoefs(6,npols,npatches)
      real *8 hxyzcoefs(9,npols,npatches)
      
      integer ntarg
      real *8 xyztarg(3,ntarg)
      real *8 xyzproxy(3,ntarg)
      integer itargptr(npatches)
      integer ntargptr(npatches)
      
      external fker
      real *8 dpars(*)
      complex *16 zpars(*)
      integer ipars(*)

      integer nqorder
      real *8 rfac
      real *8 rn1
      integer n2

      complex *16 cintvals(npols,ntarg)

c
cc      temporary variables
c
      real *8, allocatable :: tricm(:,:,:)
      real *8, allocatable :: trirad(:,:)
      integer, allocatable :: itrireltmp(:,:)
      integer, allocatable :: itrirel(:,:)
      integer, allocatable :: itrirelall(:)
      real *8, allocatable :: tverts(:,:,:)

      integer laddr(2,0:200)
      integer ntri,ntrimax,nlev,itri,istart,i,j
      integer ier,itarg,jj,jstart,nlmax,npts
      integer iqtri,ii

      integer npmax

      real *8 uvsq(2,10000),wts(10000)
      real *8 uvtmp(2,1000)
      integer nqpols
      real *8, allocatable :: sigvals(:,:)
      real *8, allocatable :: pcoefs(:,:),pvals(:,:),qwts(:)
      complex *16, allocatable :: fkervals(:),sigmatmp(:,:)
      real *8 rtmp(npols,2)
      real *8 xyztmp(3)
      complex *16 ima
      integer n1,nnt

      data ima/(0.0d0,1.0d0)/

      ntrimax = 20000
c
cc      max number of levels
c
      nlmax = 200
      allocate(tricm(3,npatches,ntrimax),
     1   trirad(npatches,ntrimax),tverts(2,3,ntrimax))
      allocate(itrireltmp(ntarg,ntrimax))

      do i=1,ntrimax
        do j=1,ntarg
          itrireltmp(j,i) = 0
        enddo
      enddo
      
      ntri = 0
      nlev = 0
      ier = 0

      call gettritree(npatches,norder,npols,xyzcoefs,ntarg,xyzproxy,
     1       itargptr,ntargptr,ntrimax,nlmax,rfac,ntri,nlev,laddr,
     2       tricm,trirad,tverts,itrireltmp,ier)

cc      call prinf('ntri=*',ntri,1)
cc      call prinf('nlev=*',nlev,1)



      allocate(itrirel(ntri,ntarg),itrirelall(ntri))

c
c       transpose the itrirel array
c  

c
cc       later have the option to generate statistics 
c        for the difference between ntri 
c        and actual  number of triangles needed 
c        for the computation.
c        Currently not using that optimization 
c        in order to avoid if statements

      do itri=1,ntri
cc        itrirelall(itri) = 0
        do j=1,ntarg
cc          if(itrireltmp(j,itri).eq.1) itrirelall(itri) = 1
          itrirel(itri,j) = itrireltmp(j,itri)
        enddo
      enddo
c
c       compute all the geometry info and koornwinder polynomials
c       at relevant triangle given by itrirelall
c  

      call triasymq(nqorder,tverts(1,1,1),tverts(1,2,1),tverts(1,3,1),
     1       uvsq,wts,nqpols)    


      npmax = ntri*nqpols
      allocate(sigvals(npols,npmax))
      allocate(pcoefs(9,npols))
      allocate(pvals(9,npmax),qwts(npmax))

      allocate(sigmatmp(npols,npmax),fkervals(npmax))


      do itri=1,ntri
        call mapuv(tverts(1,1,itri),nqpols,uvsq,uvtmp)
        istart = (itri-1)*nqpols
        do i=1,nqpols
          ii = istart+i
          call koorn_pols(uvtmp(1,i),norder,npols,sigvals(1,ii))
        enddo
      enddo

      n1 = 0
      n2 = 0
      do itri=1,npatches
c
cc       for the current patch compute all the geometry info
c
        do i=1,npols
          do j=1,3
            pcoefs(j,i) = xyzcoefs(j,i,itri)
          enddo

          do j=1,6
            pcoefs(j+3,i) = dxyzcoefs(j,i,itri)
          enddo
        enddo

        call dmatmat(9,npols,pcoefs,npmax,sigvals,pvals)
        call getqwts(nqpols,nlev,npmax,laddr,wts,pvals,qwts)
        do itarg=itargptr(itri),itargptr(itri)+ntargptr(itri)-1
          dpars(1) = xyzproxy(1,itarg)
          dpars(2) = xyzproxy(2,itarg)
          dpars(3) = xyzproxy(3,itarg)

          ipars(1) = itarg

c
c           extract info of geometry on subtriangles relevant
c           for this computation
c
          npts = 0
          nnt = 0
          do iqtri=1,ntri
            jstart = (iqtri-1)*nqpols
            if(itrirel(iqtri,itarg).eq.1) then
              nnt = nnt + 1
              do i=1,nqpols
                jj = jstart+i
                ii = npts+i
                xyztmp(1) = pvals(1,jj)
                xyztmp(2) = pvals(2,jj)
                xyztmp(3) = pvals(3,jj)
cc                call prin2('xyztmp=*',xyztmp,3)
cc                call prinf('ipars=*',ipars,4)
cc                call prin2('zpars=*',zpars,12)
                call fker(xyztarg(1,itarg),xyztmp,dpars,zpars,ipars,
     1                  fkervals(ii))
                fkervals(ii) = fkervals(ii)*qwts(jj)
                do j=1,npols
                  sigmatmp(j,ii) = sigvals(j,jj)
                enddo
              enddo
              npts = npts + nqpols
            endif
          enddo
cc          call prinf('npts=*',npts,1)
cc          call prin2('fkervals=*',fkervals,2*npts)
          call zmatvec2(npols,npts,sigmatmp,fkervals,cintvals(1,itarg))
          n1 = n1 + nnt
          if(nnt.gt.n2) n2 = nnt
        enddo
      enddo

      rn1 = (n1+0.0d0)/(ntarg+0.0d0)

      return
      end
c
c
c
c
c
c
      subroutine ctriaints3_withmetrics(eps,
     1     npatches,norder,npols,xyzcoefs,
     2     dxyzcoefs,hxyzcoefs,ntarg,xyztarg,xyzproxy,
     3     itargptr,ntargptr,fker,dpars,zpars,ipars,nqorder,rfac,
     4     cintvals,rn1,n2)

c
      implicit none

c
cc     calling sequence variables
c
      real *8 eps
      integer npatches,norder,npols
      real *8 xyzcoefs(3,npols,npatches)
      real *8 dxyzcoefs(6,npols,npatches)
      real *8 hxyzcoefs(9,npols,npatches)
      
      integer ntarg
      real *8 xyztarg(3,ntarg)
      real *8 xyzproxy(3,ntarg)
      integer itargptr(npatches)
      integer ntargptr(npatches)
      
      external fker
      real *8 dpars(*)
      complex *16 zpars(*)
      integer ipars(*)

      integer nqorder
      real *8 rfac,rn1
      integer n2

      complex *16 cintvals(npols,ntarg)

c
c       tree variables
c
      integer nlmax,ltree
      real *8, allocatable :: tverts(:,:,:),da(:)
      integer, allocatable :: itree(:)
      integer iptr(6)

      integer ntri,ntrimax,nlev,itri,istart,i,j
      integer ier,itarg,jj,jstart,npts
      integer iqtri,ii

      integer npmax

      real *8 uvsq(2,10000),wts(10000)
      real *8 uvtmp(2,1000)
      integer nqpols
      real *8, allocatable :: sigvals(:,:)
      real *8, allocatable :: pcoefs(:,:),pvals(:,:),qwts(:)

      integer n1,nnt
      
c
c      get the tree
c 
      nlmax = 6
      ntri = (4**(nlmax+1)-1)/3.0d0

      ltree = 2*(nlmax+1) + 6*ntri

      allocate(tverts(2,3,ntri),da(0:nlmax))
      allocate(itree(ltree))

      call triatreeunifini(nlmax,ntri,ltree,tverts,da,itree,iptr)

c
c       compute all the geometry info and koornwinder polynomials
c  

      call triasymq(nqorder,tverts(1,1,1),tverts(1,2,1),tverts(1,3,1),
     1       uvsq,wts,nqpols)    

cc      call prinf('nqpols=*',nqpols,1)

      npmax = ntri*nqpols
      allocate(sigvals(npols,npmax))
      allocate(pcoefs(9,npols))
      allocate(pvals(9,npmax),qwts(npmax))

      do itri=1,ntri
        call mapuv(tverts(1,1,itri),nqpols,uvsq,uvtmp)
        istart = (itri-1)*nqpols
        do i=1,nqpols
          ii = istart+i
          call koorn_pols(uvtmp(1,i),norder,npols,sigvals(1,ii))
        enddo
      enddo

      n1 = 0
      n2 = 0
      do itri=1,npatches
c
cc       for the current patch compute all the geometry info
c
        do i=1,npols
          do j=1,3
            pcoefs(j,i) = xyzcoefs(j,i,itri)
          enddo

          do j=1,6
            pcoefs(j+3,i) = dxyzcoefs(j,i,itri)
          enddo
        enddo

        call dmatmat(9,npols,pcoefs,npmax,sigvals,pvals)
        call getqwts(nqpols,nlmax,npmax,itree(1),wts,pvals,qwts)
        do itarg=itargptr(itri),itargptr(itri)+ntargptr(itri)-1

          dpars(1) = xyzproxy(1,itarg)
          dpars(2) = xyzproxy(2,itarg)
          dpars(3) = xyzproxy(3,itarg)

          ipars(1) = itarg

          call triaadap_withmetrics(eps,nqorder,nqpols,nlmax,
     1           ntri,npmax,
     1          itree,ltree,iptr,pvals,qwts,npols,sigvals,
     2          xyztarg(1,itarg),fker,dpars,zpars,ipars,
     3          cintvals(1,itarg),nnt)
          n1 = n1 + nnt
          if(n2.lt.nnt) n2 = nnt
          
cc          call prin2('cintvals=*',cintvals,3)
        enddo
      enddo
      rn1 = (n1+0.0d0)/(ntarg+0.0d0)

      return
      end

c
c
c
c
c
      subroutine triaadap_withmetrics(eps,m,kpols,nlmax,ntmax,
     1             npmax,itree,
     1             ltree,iptr,pvals,qwts,ksigpols,
     2             sigvals,xt,fker,dpars,zpars,ipars,cintall,nnt)

c
c       this subroutine adaptively computes the integral
c        of the functions
c   
c        \int_{T} 1/|xt- y(u,v)| K_{n,m}(u,v) |J(u,v)| du dv
c        n,m = 1,2\ldots ksigpols
c
c        K_{n,m}(u,v) are the koornwinder polynomials on the 
c        standard simplex (0,0)-(1,0)-(0,1)
c
c         |J(u,v)| = |dy/du \times dy/dv|
c        
c
c        relevant parameters on a dyadically refined
c        grid along the way
c
c        IN:
c        eps - precision requested
c        m - quadrature order
c        kpols - number of quarature nodes (look up triasymq 
c                 for getting kpols(m))
c        nlmax - max level of refinement for geometry
c        ntmax - max number of triangles = (4**(nlmax+1)-1)/3
c        npmax - max number of points = (4**(nlmax+1)-1)/3*kpols
c        itree - quad tree for storing the what geometry information
c                is already computed
c        ltree - length of itree
c        iptr - pointer for itree array
c               iptr(1) = laddr
c               iptr(2) = ilevel
c               iptr(3) = iparent
c               iptr(4) = ichild
c        pvals(9,npmax) - parameters computed along the adaptive grid
c        qwts(npmax) - quadrature weights 
c        ksigpols - total number of koornwinder polynomials to be
c                   integrated = (ksigorder+1)(ksigorder+2)/2
c        sigvals(ksigpols,npmax) - 
c                   koornwinder polynomials computed along the adaptive grid
c        
c        OUT:
c        cintall(ksigpols) - computed integral 
c
c         
      implicit real *8 (a-h,o-z)
      integer itree(ltree),iptr(5)
      integer, allocatable :: irefineflag(:)
      integer idone
      real *8 sigvals(ksigpols,npmax)
      real *8 pvals(9,npmax),qwts(npmax)
      complex *16, allocatable :: xkernvals(:)
      real *8 xt(3),xs(3)
      complex *16 cintall(ksigpols),fval,ctmp(ksigpols)
      complex *16, allocatable :: cvals(:,:)

      real *8 dpars(*)
      complex *16 zpars(*)
      integer ipars(*)

      external fker
 
      allocate(irefineflag(ntmax),cvals(ksigpols,ntmax))
      kfine = 4*kpols

      do i=1,ntmax
        irefineflag(i) = 0
      enddo

      nnt = 0

      irefineflag(1) = 1
      do i=1,ksigpols
         cvals(i,1) = 0
      enddo
      allocate(xkernvals(npmax))

c
cc      compute integral at level 0
c
      do i=1,kpols
         call fker(xt,pvals(1,i),dpars,zpars,ipars,fval)
         xkernvals(i) = fval*qwts(i)
         do j=1,ksigpols
            cvals(j,1) = cvals(j,1)+xkernvals(i)*sigvals(j,i)
         enddo
      enddo

      nnt = nnt + 1

      
      do i=1,ksigpols
         cintall(i) = cvals(i,1)
      enddo


      do ilev=0,nlmax-1
        idone = 1
        do itri = itree(2*ilev+1),itree(2*ilev+2)

c
cc       check if triangle needs to be refined
c
          if(irefineflag(itri).eq.1) then

            irefineflag(itri) = 0
            itric1 = itree(iptr(4)+4*(itri-1))
            do i=1,4
              irefineflag(itric1+i-1) = 0
            enddo
c
cc           compute xkernvals
c
            istart = (itric1-1)*kpols
            do j=1,kfine
              jj=j+istart
              call fker(xt,pvals(1,jj),dpars,zpars,ipars,fval)
              xkernvals(jj) = fval*qwts(jj)
            enddo

            nnt = nnt + 4

c
cc           subtract contribution of current 
c            triangle
            do isig=1,ksigpols
              cintall(isig) = cintall(isig)-cvals(isig,itri)
              ctmp(isig) = 0
            enddo

c
cc           add in contributions of children triangles
c
            do i=1,4
              itric = itree(iptr(4)+4*(itri-1)+i-1)
              do isig=1,ksigpols
                cvals(isig,itric) = 0
              enddo

              istart = (itric-1)*kpols
cc              call prinf('istart=*',istart,1)
              do k=1,kpols
                ii = istart+k
                do isig=1,ksigpols
                  cvals(isig,itric) = cvals(isig,itric)+xkernvals(ii)*
     1                                  sigvals(isig,ii)
                enddo
              enddo
              
              do isig=1,ksigpols
                cintall(isig) = cintall(isig) + cvals(isig,itric)
                ctmp(isig) = ctmp(isig)+cvals(isig,itric)
              enddo
            enddo

c
cc           compare integral of children to integral of parent
c            to determine if the children need to be refined further
c
            errmax = 0
            do isig=1,ksigpols
              if(abs(ctmp(isig)-cvals(isig,itri)).gt.errmax) 
     1          errmax = abs(ctmp(isig)-cvals(isig,itri))
            enddo

            if(errmax.gt.eps) then
              idone = 0
              do j=1,4
                itric = itree(iptr(4)+4*(itri-1)+j-1)
                irefineflag(itric) = 1
              enddo
            endif
          endif
c
cc        end of looping over all triangles
        enddo
cc         if idone is still 1, that means that no more refinement
c          is needed
         if(idone.eq.1) goto 1111
      enddo
 1111 continue


      return
      end
      
c
c
c
cc
c
c
      subroutine ctriaints4_withmetrics(eps,
     1     npatches,norder,npols,xyzcoefs,
     2     dxyzcoefs,hxyzcoefs,ntarg,xyztarg,xyzproxy,
     3     itargptr,ntargptr,fker,dpars,zpars,ipars,nqorder,rfac,
     4     cintvals,rn1,n2)
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
c        this code is the same as ctriaints on all counts except
c        for choosing when to subdivide a triangle
c        
c        here a triangle is subdivided when xyzproxy associated
c        with the triangle is separated from the centroid
c        of the triangle by at least rfac*r_{t} where r_{t}
c        is the radius of the smallest sphere  (see point 2 of
c        ctriaints)
c
c        input arguments:
c        eps:     requested precision (Currently not being used
c                 in any way, later, nqorder, and rfac will
c                 be set internally based on eps)
c
c        npatches: number of patches
c        norder: order of discretization nodes on the patches
c        npols = (norder+1)*(norder+2)/2: number of discretization 
c                   nodes on each patch
c        xyzcoefs(3,npols,npatches): coefficients
c                 of koornwinder expansion of xyz coordinates
c        dxyzcoefs(6,npols,npatches):
c                 coefficients of koornwinder expansion 
c                 of dxyz/duv
c        hxyzcoefs(9,npols,npatches):
c                 coefficients of koornwinder expansion
c                 of d2xyz/d2uv, 
c                 1 = d2x/duu
c                 2 = d2x/duv
c                 3 = d2x/dvv
c                 4 = d2y/duu
c                      .
c                      .
c                 9 = d2z/dvv
c
c        ntarg - total number of target points
c        xyztarg(3,ntarg) - location of target points 
c                       (currently may have repeats, 
c                        being slightly wasteful in
c                        memory, but see comment (1))
c
c        xyzproxy(3,ntarg) - proxy locations for measuring distances
c
c        itargptr(npatches) - xyztargs(:,itargptr(i):itargptr(i)+ntargptr(i)-1)
c                       are the relevant set of targets for patch i
c        ntargptr(npatches) - number of relevant targets for patch i
c        fker - function handle for evaluating the kernel k
c 
c               expected calling sequence
c               fker(x,y,dpars,zpars,ipars,f)
c               
c               the output is assumed to be complex for the time
c               being
c
c         dpars(*) - real parameters for the fker routine
c         zpars(*) - complex parameters for the fker routine
c         ipars(*) - integer parameters for the fker routine
c         nqorder - order of quadrature nodes on each subtriangle
c                   to be used
c         rfac - distance criterion for deciding whether a triangle
c                is in the far-field or not (See comment (2))
c
c         output:
c
c         cintvals(npols,ntarg) - integrals at all targets
c                                  for all koornwinder
c                                  polynomials
c
      implicit none

c
cc     calling sequence variables
c
      real *8 eps
      integer npatches,norder,npols
      real *8 xyzcoefs(3,npols,npatches)
      real *8 dxyzcoefs(6,npols,npatches)
      real *8 hxyzcoefs(9,npols,npatches)
      
      integer ntarg
      real *8 xyztarg(3,ntarg)
      real *8 xyzproxy(3,ntarg)
      integer itargptr(npatches)
      integer ntargptr(npatches)
      
      external fker
      real *8 dpars(*)
      complex *16 zpars(*)
      integer ipars(*)

      integer nqorder
      real *8 rfac,rn1

      complex *16 cintvals(npols,ntarg)
      integer n2

c
cc      temporary variables
c
      real *8, allocatable :: tricm(:,:,:)
      real *8, allocatable :: trirad(:,:)
      integer, allocatable :: itrireltmp(:,:)
      integer, allocatable :: itrirel(:,:)
      integer, allocatable :: itrirelall(:)
      real *8, allocatable :: tverts(:,:,:)

      integer laddr(2,0:200)
      integer ntri,ntrimax,nlev,itri,istart,i,j
      integer ier,itarg,jj,jstart,nlmax,npts
      integer iqtri,ii

      integer npmax

      real *8 uvsq(2,10000),wts(10000)
      real *8, allocatable :: umatrq(:,:),vmatrq(:,:)
      real *8 uvtmp(2,1000)
      integer nqpols
      real *8, allocatable :: sigvals(:,:)
      real *8, allocatable :: pcoefs(:,:),pvals(:,:),qwts(:)
      complex *16, allocatable :: fkervals(:),sigmatmp(:,:)
      real *8 rtmp(npols,2)
      real *8 xyztmp(3)
      complex *16 ima
      integer n1,nnt

      data ima/(0.0d0,1.0d0)/


      ntrimax = 20000
c
cc      max number of levels
c
      nlmax = 200
      allocate(tricm(3,npatches,ntrimax),
     1   trirad(npatches,ntrimax),tverts(2,3,ntrimax))
      allocate(itrireltmp(ntarg,ntrimax))

      do i=1,ntrimax
        do j=1,ntarg
          itrireltmp(j,i) = 0
        enddo
      enddo
      
      ntri = 0
      nlev = 0
      ier = 0

cc      call prin2('rfac=*',rfac,1)

      call gettritree(npatches,norder,npols,xyzcoefs,ntarg,xyzproxy,
     1       itargptr,ntargptr,ntrimax,nlmax,rfac,ntri,nlev,laddr,
     2       tricm,trirad,tverts,itrireltmp,ier)



      allocate(itrirel(ntri,ntarg),itrirelall(ntri))

c
c       transpose the itrirel array
c  

c
cc       later have the option to generate statistics 
c        for the difference between ntri 
c        and actual  number of triangles needed 
c        for the computation.
c        Currently not using that optimization 
c        in order to avoid if statements


      do itri=1,ntri
        do j=1,ntarg
          itrirel(itri,j) = itrireltmp(j,itri)
        enddo
      enddo
c
c       compute all the geometry info and koornwinder polynomials
c       at relevant triangle given by itrirelall
c  
      nqpols = (nqorder+1)*(nqorder+2)/2
      allocate(umatrq(nqpols,nqpols),vmatrq(nqpols,nqpols))
      
      call vioreanu_simplex_quad(nqorder,nqpols,uvsq,umatrq,vmatrq,wts)

      npmax = ntri*nqpols
      allocate(sigvals(npols,npmax))
      allocate(pcoefs(9,npols))
      allocate(pvals(9,npmax),qwts(npmax))

      allocate(sigmatmp(npols,npmax),fkervals(npmax))

      n1 = 0
      n2 = 0
      do itri=1,ntri
        call mapuv(tverts(1,1,itri),nqpols,uvsq,uvtmp)
        istart = (itri-1)*nqpols
        do i=1,nqpols
          ii = istart+i
          call koorn_pols(uvtmp(1,i),norder,npols,sigvals(1,ii))
        enddo
      enddo


      do itri=1,npatches
c
cc       for the current patch compute all the geometry info
c
        do i=1,npols
          do j=1,3
            pcoefs(j,i) = xyzcoefs(j,i,itri)
          enddo

          do j=1,6
            pcoefs(j+3,i) = dxyzcoefs(j,i,itri)
          enddo
        enddo

        call dmatmat(9,npols,pcoefs,npmax,sigvals,pvals)
        call getqwts(nqpols,nlev,npmax,laddr,wts,pvals,qwts)
        do itarg=itargptr(itri),itargptr(itri)+ntargptr(itri)-1
          dpars(1) = xyzproxy(1,itarg)
          dpars(2) = xyzproxy(2,itarg)
          dpars(3) = xyzproxy(3,itarg)

          ipars(1) = itarg
c
c           extract info of geometry on subtriangles relevant
c           for this computation
c
          npts = 0
          nnt = 0
          do iqtri=1,ntri
            jstart = (iqtri-1)*nqpols
            if(itrirel(iqtri,itarg).eq.1) then
              nnt = nnt + 1
              do i=1,nqpols
                jj = jstart+i
                ii = npts+i
                xyztmp(1) = pvals(1,jj)
                xyztmp(2) = pvals(2,jj)
                xyztmp(3) = pvals(3,jj)
cc                call prin2('xyztmp=*',xyztmp,3)
cc                call prinf('ipars=*',ipars,4)
cc                call prin2('zpars=*',zpars,12)
                call fker(xyztarg(1,itarg),xyztmp,dpars,zpars,ipars,
     1                  fkervals(ii))
                fkervals(ii) = fkervals(ii)*qwts(jj)
                do j=1,npols
                  sigmatmp(j,ii) = sigvals(j,jj)
                enddo
              enddo
              npts = npts + nqpols
            endif
          enddo
cc          call prinf('npts=*',npts,1)
cc          call prin2('fkervals=*',fkervals,2*npts)
          call zmatvec2(npols,npts,sigmatmp,fkervals,cintvals(1,itarg))
          n1 = n1 + nnt
          if(nnt.gt.n2) n2 = nnt
        enddo
      enddo

      rn1 = (n1+0.0d0)/(ntarg+0.0d0)

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
c
c
      subroutine ctriaints3_withtrip(eps,
     1     npatches,norder,npols,xyzcoefs,
     2     dxyzcoefs,hxyzcoefs,ntarg,xyztarg,xyzproxy,
     3     itargptr,ntargptr,fker,dpars,zpars,ipars,nqorder,rfac,
     4     cintvals)

c
      implicit none

c
cc     calling sequence variables
c
      real *8 eps
      integer npatches,norder,npols
      real *8 xyzcoefs(3,npols,npatches)
      real *8 dxyzcoefs(6,npols,npatches)
      real *8 hxyzcoefs(9,npols,npatches)
      
      integer ntarg
      real *8 xyztarg(3,ntarg)
      real *8 xyzproxy(3,ntarg)
      integer itargptr(npatches)
      integer ntargptr(npatches)

      external fker
      real *8 dpars(*)
      complex *16 zpars(*)
      integer ipars(*)

      integer nqorder
      real *8 rfac

      complex *16 cintvals(npols,ntarg)

c
c       tree variables
c
      integer nlmax,ltree
      real *8, allocatable :: tverts(:,:,:),da(:)
      integer, allocatable :: itree(:)
      integer iptr(6)

      integer ntri,ntrimax,nlev,itri,istart,i,j
      integer ier,itarg,jj,jstart,npts
      integer iqtri,ii

      integer npmax

      real *8 uvsq(2,10000),wts(10000)
      real *8 uvtmp(2,1000)
      integer nqpols
      real *8, allocatable :: sigvals(:,:)
      real *8, allocatable :: pcoefs(:,:),pvals(:,:),qwts(:)

      integer, allocatable :: iscur(:)
      integer iw
      
c
c      get the tree
c 
      nlmax = 6
      ntri = (4**(nlmax+1)-1)/3.0d0
      ltree = 2*(nlmax+1) + 6*ntri

      allocate(tverts(2,3,ntri),da(0:nlmax))
      allocate(itree(ltree))

      call triatreeunifini(nlmax,ntri,ltree,tverts,da,itree,iptr)

c
c       compute all the geometry info and koornwinder polynomials
c  

      call triasymq(nqorder,tverts(1,1,1),tverts(1,2,1),tverts(1,3,1),
     1       uvsq,wts,nqpols)    

cc      call prinf('nqpols=*',nqpols,1)

      npmax = ntri*nqpols
      allocate(sigvals(npols,npmax))
      allocate(pcoefs(9,npols))
      allocate(pvals(9,npmax),qwts(npmax))


      do itri=1,ntri
        call mapuv(tverts(1,1,itri),nqpols,uvsq,uvtmp)
        istart = (itri-1)*nqpols
        do i=1,nqpols
          ii = istart+i
          call koorn_pols(uvtmp(1,i),norder,npols,sigvals(1,ii))
        enddo
      enddo

      allocate(iscur(ntri))



      do itri=1,npatches
c
cc       for the current patch compute all the geometry info
c
        do i=1,npols
          do j=1,3
            pcoefs(j,i) = xyzcoefs(j,i,itri)
          enddo

          do j=1,6
            pcoefs(j+3,i) = dxyzcoefs(j,i,itri)
          enddo
        enddo

        call dmatmat(9,npols,pcoefs,npmax,sigvals,pvals)
        call getqwts(nqpols,nlmax,npmax,itree(1),wts,pvals,qwts)
        do itarg=itargptr(itri),itargptr(itri)+ntargptr(itri)-1

          do ii=1,ntri
            iscur(ii) = 0
          enddo

          dpars(1) = xyzproxy(1,itarg)
          dpars(2) = xyzproxy(2,itarg)
          dpars(3) = xyzproxy(3,itarg)

          ipars(1) = itarg

          call triaadap_withtrip(eps,nqorder,nqpols,nlmax,ntri,npmax,
     1          itree,ltree,iptr,pvals,qwts,npols,sigvals,
     2          xyztarg(1,itarg),fker,dpars,zpars,ipars,iscur,
     3          cintvals(1,itarg))

          iw = 33 + itarg
          call printtri(iw,ntri,tverts,iscur)

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
      subroutine triaadap_withtrip(eps,m,kpols,nlmax,ntmax,npmax,itree,
     1             ltree,iptr,pvals,qwts,ksigpols,
     2             sigvals,xt,fker,dpars,zpars,ipars,iscur,cintall)

c
c       this subroutine adaptively computes the integral
c        of the functions
c   
c        \int_{T} 1/|xt- y(u,v)| K_{n,m}(u,v) |J(u,v)| du dv
c        n,m = 1,2\ldots ksigpols
c
c        K_{n,m}(u,v) are the koornwinder polynomials on the 
c        standard simplex (0,0)-(1,0)-(0,1)
c
c         |J(u,v)| = |dy/du \times dy/dv|
c        
c
c        relevant parameters on a dyadically refined
c        grid along the way
c
c        IN:
c        eps - precision requested
c        m - quadrature order
c        kpols - number of quarature nodes (look up triasymq 
c                 for getting kpols(m))
c        nlmax - max level of refinement for geometry
c        ntmax - max number of triangles = (4**(nlmax+1)-1)/3
c        npmax - max number of points = (4**(nlmax+1)-1)/3*kpols
c        itree - quad tree for storing the what geometry information
c                is already computed
c        ltree - length of itree
c        iptr - pointer for itree array
c               iptr(1) = laddr
c               iptr(2) = ilevel
c               iptr(3) = iparent
c               iptr(4) = ichild
c        pvals(9,npmax) - parameters computed along the adaptive grid
c        qwts(npmax) - quadrature weights 
c        ksigpols - total number of koornwinder polynomials to be
c                   integrated = (ksigorder+1)(ksigorder+2)/2
c        sigvals(ksigpols,npmax) - 
c                   koornwinder polynomials computed along the adaptive grid
c        
c        OUT:
c        cintall(ksigpols) - computed integral 
c
c         
      implicit real *8 (a-h,o-z)
      integer itree(ltree),iptr(5)
      integer, allocatable :: irefineflag(:)
      integer idone
      real *8 sigvals(ksigpols,npmax)
      real *8 pvals(9,npmax),qwts(npmax)
      complex *16, allocatable :: xkernvals(:)
      real *8 xt(3),xs(3)
      complex *16 cintall(ksigpols),fval,ctmp(ksigpols)
      complex *16, allocatable :: cvals(:,:)

      real *8 dpars(*)
      complex *16 zpars(*)
      integer ipars(*)
      integer iscur(*)

      external fker
 
      allocate(irefineflag(ntmax),cvals(ksigpols,ntmax))
      kfine = 4*kpols

      nfunev = 0

      do i=1,ntmax
        irefineflag(i) = 0
      enddo

     

      irefineflag(1) = 1
      do i=1,ksigpols
         cvals(i,1) = 0
      enddo
      allocate(xkernvals(npmax))


c
cc      compute integral at level 0
c

      iscur(1) = 1
      do i=1,kpols
         call fker(xt,pvals(1,i),dpars,zpars,ipars,fval)
         xkernvals(i) = fval*qwts(i)
         do j=1,ksigpols
            cvals(j,1) = cvals(j,1)+xkernvals(i)*sigvals(j,i)
         enddo
      enddo

      nfunev = nfunev + kpols

      
      do i=1,ksigpols
         cintall(i) = cvals(i,1)
      enddo


      do ilev=0,nlmax-1
        idone = 1
        do itri = itree(2*ilev+1),itree(2*ilev+2)

c
cc       check if triangle needs to be refined
c
          if(irefineflag(itri).eq.1) then

            irefineflag(itri) = 0
            itric1 = itree(iptr(4)+4*(itri-1))
            do i=1,4
              irefineflag(itric1+i-1) = 0
            enddo
c
cc           compute xkernvals
c
            istart = (itric1-1)*kpols
            do j=1,kfine
              jj=j+istart
              call fker(xt,pvals(1,jj),dpars,zpars,ipars,fval)
              xkernvals(jj) = fval*qwts(jj)
            enddo

            nfunev = nfunev + kfine

c
cc           subtract contribution of current 
c            triangle
            do isig=1,ksigpols
              cintall(isig) = cintall(isig)-cvals(isig,itri)
              ctmp(isig) = 0
            enddo

c
cc           add in contributions of children triangles
c
            do i=1,4
              itric = itree(iptr(4)+4*(itri-1)+i-1)
              iscur(itric) = 1

              do isig=1,ksigpols
                cvals(isig,itric) = 0
              enddo

              istart = (itric-1)*kpols
              do k=1,kpols
                ii = istart+k
                do isig=1,ksigpols
                  cvals(isig,itric) = cvals(isig,itric)+xkernvals(ii)*
     1                                  sigvals(isig,ii)
                enddo
              enddo
              
              do isig=1,ksigpols
                cintall(isig) = cintall(isig) + cvals(isig,itric)
                ctmp(isig) = ctmp(isig)+cvals(isig,itric)
              enddo
            enddo

c
cc           compare integral of children to integral of parent
c            to determine if the children need to be refined further
c
            errmax = 0
            do isig=1,ksigpols
              if(abs(ctmp(isig)-cvals(isig,itri)).gt.errmax) 
     1          errmax = abs(ctmp(isig)-cvals(isig,itri))
            enddo

            if(errmax.gt.eps) then
              idone = 0
              do j=1,4
                itric = itree(iptr(4)+4*(itri-1)+j-1)
                irefineflag(itric) = 1
              enddo
            endif
          endif
c
cc        end of looping over all triangles
        enddo
cc         if idone is still 1, that means that no more refinement
c          is needed
         if(idone.eq.1) goto 1111
      enddo
 1111 continue


      return
      end
c
c
c---------------------------------
      subroutine printtri(iw,ntri,tvs,iprint)
      implicit real *8 (a-h,o-z)
      real *8 tvs(2,3,ntri)
      integer iprint(ntri)

 2100 format(4(2x,e11.5))
      do i=1,ntri
        if(iprint(i).eq.1) then
           write(iw,2100) tvs(1,1,i),tvs(1,2,i),tvs(2,1,i),tvs(2,2,i)
           write(iw,2100) tvs(1,2,i),tvs(1,3,i),tvs(2,2,i),tvs(2,3,i)
           write(iw,2100) tvs(1,3,i),tvs(1,1,i),tvs(2,3,i),tvs(2,1,i)
        endif
      enddo

      return
      end

c---------------------------------      
