c
c  Notes: The subroutines in this file assume that all patches
c  are discretized with the same order and same type of
c  polynomials
c
c  Typical usage of this subroutine is with one patch only
c
c  Routines in this file
c   - cquadints_adap: completely adaptive integration
c                       for the quad patch
c
c   - cquadints_wnodes: compute integration with prescribed nodes
c                       and weights
c           
c
c  We integrate against chebyshev/legendre of full order or 
c  total order polynomials on the standard patch [-1,1]^2
c                       
c
c

      subroutine cquadints_wnodes(npatches,norder,ipoly,ttype,npols,
     1    srccoefs,ndtarg,ntarg,xyztarg,itargptr,ntargptr,
     2    nporder,nppols,fker,ndd,dpars,ndz,zpars,ndi,ipars,
     3    nqpts,qnodes,qwts,cintvals)
c
c  This subroutine computes the integrals
c
c    \int_{[-1,1]^2} K(x_{i},\rho_{m}(y)) B_{n}(y) J_{m}(y) dy \, , 
c
c  B_{n}(y) are either tensor product legendre polynomials or
c  Chebyshev polynomials on   [-1,1]^2 of either total degree 
c  or full degree 
c
c     J_{m}(y) = det(D\rho_{m}^{T} \cdot D\rho_{m})
c
c  The maps \rho_{m} are stored in the same basis as the 
c  basis functions used in the integration  along with
c  expansions of first derivative information of xyz with 
c  respect to u,v
c       
c  The subroutine uses prescribed nodes and weights to compute 
c  the integrals.
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
c    - qwts: real *8 (nqpts)
c        the corresponding quadrature weights
c 
c  Output arguments:
c    - cintvals: complex *16 (npols,ntarg)
c        Integral against all basis functions from patches
c        to targets
c
c------------------------------------- 
c
      implicit none
      integer norder,npols,ntarg,ipars(*),nqpts
      character *1 ttype
      real *8 xyztarg(3,ntarg),dpars(*),qnodes(2,nqpts),qwts(nqpts)
      complex *16 zpars(*),cintvals(npols,ntarg),alpha,beta,fval

      complex *16, allocatable :: xkernvals(:,:)
      complex *16, allocatable :: sigvals(:,:)
      real *8, allocatable :: pols_tmp(:)

      integer i,j
      character *1 transa,transb
      external fker

      allocate(pols_tmp(npols))
      allocate(sigvals(npols,nqpts),xkernvals(nqpts,ntarg))

c
c
c       generate the density values
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,pols_tmp)     
      do i=1,nqpts
        call legetens_pols_2d(qnodes(1,i),norder-1,type,pols_tmp)
        do j=1,npols
          sigvals(j,i) = pols_tmp(j)
        enddo
      enddo
C$OMP END PARALLEL DO      


C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,fval)
      do i=1,ntarg
        do j=1,nqpts
          call fker(qnodes(1,j),xyztarg(1,i),dpars,zpars,ipars,fval)
          xkernvals(j,i) = fval*qwts(j)
        enddo
      enddo
C$OMP END PARALLEL DO      

      
      transa = 'n'
      transb = 'n'
      alpha = 1
      beta = 0

cc      print *, "Starting blas"
      call zgemm(transa,transb,npols,ntarg,nqpts,alpha,sigvals,npols,
     1   xkernvals,nqpts,beta,cintvals,npols)
      

      return
      end
c
c
c
c
c
c
      subroutine cquadints_adap(eps,intype,norder,ttype,npols,
     1  ntarg,xyztarg,nquadmax,fker,dpars,zpars,ipars,nqorder,
     2  iflg,cintvals)
c
c       this subroutine computes the integrals
c
c       \int_{[-1,1]^2} K(x_{i},y) P_{n}(y_{1}) P_{m}(y_{2}) dy \, ,
c
c        P_{n}(y) are Legendre polynomials on [-1,1]
c
c        using adaptive integration
c
c
c        input arguments:
c        eps:     requested precision 
c
c        intype:   quadrature node type
c                   intype = 1, tensor product gauss-legendre nodes
c                   intype = 2, xiao gimbutas nodes
c   
c        norder: order of polynomials on the patch 
c        npols = norder*norder, number of polynomials to be integrated 
c        ntarg - total number of target points
c        xyztarg(3,ntarg) - location of target points 
c        nqmax - max number of quads to be used on base quad
c        fker - function handle for evaluating the kernel k
c 
c               expected calling sequence
c               fker(x,y,ynorms,dpars,zpars,ipars,f)
c
c               x \in \mathbb{R}^{3}, y \in \mathbb{R}^{2}
c               
c               the output is assumed to be complex for the time
c               being
c
c         dpars(*) - real parameters for the fker routine
c         zpars(*) - complex parameters for the fker routine
c         ipars(*) - integer parameters for the fker routine
c         nqorder - order of quadrature nodes on each subquad
c                   to be used
c         iflg  - flag for determining if precomputed
c            tensor product legendre polynomial values
c            are resued across different targets
c            iflg = 1 => legendre polynomial values reused
c            iflg = 2 => legendre polynomila values not reused
c
c         output:
c
c         cintvals(npols,ntarg) - integrals at all targets
c                                  for all tensor product
c                                  chebyshev polynomials
c
c
c
      implicit none

c
cc     calling sequence variables
c
      real *8 eps
      integer intype
      integer norder,npols
      
      integer ntarg
      real *8 xyztarg(3,ntarg)
      
      external fker
      real *8 dpars(*)
      complex *16 zpars(*)
      integer ipars(*)

      integer nqorder

      integer nquadmax

      complex *16 cintvals(npols,ntarg)

      character ttype
      integer iflg


      if(iflg.eq.1) then
        call cquadints_adap_precomp(eps,intype,norder,ttype,
     1      npols,ntarg,xyztarg,nquadmax,fker,dpars,zpars,
     2      ipars,nqorder,cintvals)
      else
        call cquadints_adap_noprecomp(eps,intype,norder,ttype,
     1      npols,ntarg,xyztarg,nquadmax,fker,dpars,zpars,
     2      ipars,nqorder,cintvals)

      endif


      return
      end




      subroutine cquadints_adap_precomp(eps,intype,
     1     norder,type,npols,ntarg,xyztarg,nquadmax,
     3     fker,dpars,zpars,ipars,nqorder,cintvals)

c
c       this subroutine computes the integrals
c
c       \int_{[-1,1]^2} K(x_{i},y) P_{n}(y_{1}) P_{m}(y_{2}) dy \, ,
c
c        P_{n}(y) are Legendre polynomials on [-1,1]
c
c        using adaptive integration
c
c
c        input arguments:
c        eps:     requested precision 
c
c        intype:   quadrature node type
c                   intype = 1, tensor product gauss-legendre nodes
c                   intype = 2, xiao gimbutas nodes
c   
c        norder: order of polynomials on the patch 
c        npols = norder*norder, number of polynomials to be integrated 
c        ntarg - total number of target points
c        xyztarg(3,ntarg) - location of target points 
c        nqmax - max number of quads to be used on base quad
c        fker - function handle for evaluating the kernel k
c 
c               expected calling sequence
c               fker(x,y,ynorms,dpars,zpars,ipars,f)
c
c               x \in \mathbb{R}^{3}, y \in \mathbb{R}^{2}
c               
c               the output is assumed to be complex for the time
c               being
c
c         dpars(*) - real parameters for the fker routine
c         zpars(*) - complex parameters for the fker routine
c         ipars(*) - integer parameters for the fker routine
c         nqorder - order of quadrature nodes on each subquad
c                   to be used
c
c         output:
c
c         cintvals(npols,ntarg) - integrals at all targets
c                                  for all tensor product
c                                  chebyshev polynomials
c
c
c
      implicit none

c
cc     calling sequence variables
c
      real *8 eps
      integer intype
      integer norder,npols
      
      integer ntarg
      real *8 xyztarg(3,ntarg)
      
      external fker
      real *8 dpars(*)
      complex *16 zpars(*)
      integer ipars(*)

      integer nqorder

      integer nquadmax

      complex *16 cintvals(npols,ntarg)

      character type

c
c       tree variables
c
      integer nlmax,ltree
      real *8, allocatable :: tvs(:,:,:),da(:)
      integer, allocatable :: ichild_start(:)

      integer nquad,nlev,iquad,istart,i,j,k
      integer ier,itarg,jj,jstart,npts
      integer iqquad,ii


      integer npmax

      real *8, allocatable :: uvsq(:,:),wts(:),uvtmp(:,:)
      real *8, allocatable :: umattmp(:,:),vmattmp(:,:)
      integer nqpols
      real *8, allocatable :: sigvals(:,:)
      real *8, allocatable :: uvvals(:,:),qwts(:)
      integer itmp

      character *1 transa,transb
      real *8 alpha,beta,ra
      integer lda,ldb,ldc
      real *8 u, v
      integer ldu, ldv, itype,ndeg
      
c
c       for each quad, we just store three pieces
c       of info
c         quad vertices
c         area of quad
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
c
c
c        quad vertices nomenclature
c
c       v3
c        ________ 
c        |       |
c        |       |
c        |       |
c        ---------
c        v1       v2
c
c
c
      tvs(1,1,1) = -1
      tvs(2,1,1) = -1

      tvs(1,2,1) = 1
      tvs(2,2,1) = -1

      tvs(1,3,1) = -1
      tvs(2,3,1) = 1

c
c       get quadrature nodes and weights on the base quad
c       based on quadrature type
c

      if(intype.eq.1) then
        nqpols = nqorder*nqorder
        allocate(uvsq(2,nqpols),wts(nqpols))

        ldu = 1
        ldv = 1
        itype = 1
        call legetens_exps_2d(itype,nqorder,type,uvsq,u,ldu,v,ldv,wts)        
c        call get_tens_leg_nodes_2d(nqorder,nqpols,uvsq,wts)
      endif

      if(intype.eq.2) then
        call squarearbq_pts(nqorder,nqpols)
        allocate(uvsq(2,nqpols),wts(nqpols))
        call squarearbq(nqorder,uvsq,wts,nqpols)
      endif


      allocate(uvtmp(2,nqpols))
     
      npmax = nquadmax*nqpols
      allocate(sigvals(npols,npmax))
      allocate(uvvals(2,npmax),qwts(npmax))

c
c      current number of quads in the adaptive structure
c
      nquad = 1
c
c        intialize sigvals for root quad
c

      ndeg = norder-1
      call mapuv_quad(tvs(1,1,1),nqpols,uvsq,uvvals)
      do i=1,nqpols
c     call tensleg_pols(uvvals(1,i),norder,npols,sigvals(1,i))
         call legetens_pols_2d(uvvals(1,i),ndeg,type,sigvals(1,i))
        qwts(i) = wts(i)
      enddo



c
cc       for the current patch compute geometry info for base quad 
c

       
      nlmax = 20 
      do itarg=1,ntarg
        
        call quadadap(eps,nqorder,nqpols,nlmax,nquadmax,nquad,
     1    ichild_start,tvs,da,uvsq,wts, 
     1    norder,type,npols,npmax,uvvals,qwts,sigvals,xyztarg(1,itarg),
     3    fker,dpars,zpars,ipars,cintvals(1,itarg))
      enddo


      return
      end

c
c
c
c
c
      subroutine quadadap(eps,m,kpols,nlmax,nqmax,nquad,
     1             ichild_start,tvs,da,uvsq,wts,
     1             norder,type,npols,npmax,uvvals,qwts,
     2             sigvals,xt,fker,dpars,zpars,
     3             ipars,cintall)

c
c       this subroutine adaptively computes the integral
c        of the functions
c   
c        \int_{(-1,1)^2} 1/|xt- y(u,v)| P_{n}(u)P_{m}(v) du dv
c        n,m = 1,2\ldots npols
c
c        P_{n}(u) are the Legendre polynomials on (-1,1) 
c
c
c        relevant parameters on a quad refined
c        grid along the way
c
c        IN:
c        eps - precision requested
c        m - quadrature order
c        kpols - number of quarature nodes (look up squarearbq 
c                 for getting kpols(m))
c        nlmax - max level of refinement for geometry
c        nqmax - max number of quads
c        nquad - current number of quads in adaptive structure
c        ichild_start(i) - first child of quad i
c        tvs(2,3,nqmax) - vertices of hierarchy of quads
c        da(nqmax) - area of quads
c        uvsq(kpols) - integration nodes on standard quad
c        wts(kpols) - integration weights on standard quad
c        npols - total number of koornwinder polynomials to be integrated
c        norder - order of discretization of the surface
c        npols - (norder+1)*(norder+2)/2: total number of 
c                koornwinder polynomials to be integrated
c        npmax - max number of points = nqmax*kpols
c        uvvals(2,npmax) - geometry info on heirarchy of meshes
c        qwts(npmax) - quadrature weights 
c        sigvals(npols,npmax) - 
c                   tensor product GL polynomials computed along the adaptive grid
c        
c        OUT:
c        cintall(npols) - computed integral 
c
c         

      implicit real *8 (a-h,o-z)
      integer, allocatable :: istack(:)
      integer ichild_start(nqmax)
      real *8 da(nqmax)
      real *8 tvs(2,3,nqmax), uvsq(2,kpols),wts(kpols)
      integer nproclist0, nproclist
      integer idone
      real *8 sigvals(npols,npmax)
      real *8 uvvals(2,*),qwts(*)
      complex *16, allocatable :: xkernvals(:)
      real *8 xt(3),xs(2)
      complex *16 cintall(npols),fval,ctmp(npols)
      complex *16, allocatable :: cvals(:,:)

      real *8 dpars(*)
      complex *16 zpars(*)
      integer ipars(*)

      character type
      
      external fker

c
c         for historic reasons
c
      ksigpols = npols
      allocate(istack(2*nqmax))
      allocate(cvals(ksigpols,nqmax))

      do i=1,ksigpols
         cvals(i,1) = 0
      enddo

      allocate(xkernvals(npmax))

c
cc      compute integral at level 0
c
      do i=1,kpols
         call fker(uvvals(1,i),xt,dpars,zpars,ipars,fval)
         xkernvals(i) = fval*qwts(i)
         do j=1,ksigpols
            cvals(j,1) = cvals(j,1)+xkernvals(i)*sigvals(j,i)
         enddo
      enddo

      
      do i=1,ksigpols
         cintall(i) = cvals(i,1)
      enddo

      nproclist0 = 1
      istack(1) = 1


      call quadadap_main(eps,kpols,nlmax,nqmax,nquad,ichild_start,
     1      tvs,da,uvsq,wts,norder,type,npols,
     2      npmax,uvvals,qwts,sigvals,xt,fker,dpars,
     3      zpars,ipars,cvals,istack,nproclist0,
     4      xkernvals,cintall)

      
      return
      end
c
c
c
c
c
       
      subroutine quadadap_main(eps,kpols,nlmax,nqmax,nquad,
     1    ichild_start,tvs,da,uvsq,wts,norder,type,npols,
     2    npmax,uvvals,qwts,sigvals,xt,fker,dpars,
     3    zpars,ipars,cvals,istack,nproclist0,xkernvals,
     4    cintall)
      

      implicit real *8 (a-h,o-z)
      integer istack(*),nproclist0
      integer ichild_start(nqmax)
      real *8 da(nqmax)
      real *8 tvs(2,3,nqmax), uvsq(2,kpols),wts(kpols)
      integer  nproclist
      integer idone
      real *8 sigvals(npols,npmax)
      complex *16 xkernvals(npmax)
      real *8 xt(3)
      real *8 uvvals(2,*),qwts(*)
      complex *16 cintall(npols),fval,ctmp(npols)
      complex *16 cvals(npols,nqmax)

      character type

      real *8 dpars(*)
      complex *16 zpars(*)
      integer ipars(*)

      real *8, allocatable :: uvtmp(:,:)
      character *1 transa,transb
      integer lda,ldb,ldc
      external fker
      
      allocate(uvtmp(2,kpols))


c
c         for historic reasons
c
      ksigpols = npols
      kfine = 4*kpols


      do ilev=0,nlmax
        idone = 1
        nproclist = 0
        
        do iproc = 1,nproclist0
          iquad = istack(iproc)


c
c           check to see if quad already has 
c           children, if not, set children
c           and compute necessary info
c
          if(ichild_start(iquad).eq.-1) then

c
c            current quad doesn't have children,
c            compute necessary info
c

            if(nquad+4.gt.nqmax) then
               print *, "Too many quads in cquadadap"
               print *, "Exiting before convergence"

               return
            endif
            
            ichild_start(iquad) = nquad+1
            call getquadchildren(tvs(1,1,iquad),tvs(1,1,nquad+1),
     1             tvs(1,1,nquad+2),tvs(1,1,nquad+3),tvs(1,1,nquad+4))
            

            ndeg = norder-1

            rr = 0.25d0*da(iquad)
            do j=nquad+1,nquad+4
              da(j) = rr
              istart = (j-1)*kpols+1
              call mapuv_quad(tvs(1,1,j),kpols,uvsq,uvvals(1,istart))
              do i=1,kpols
                ii = istart+i-1
c                call tensleg_pols(uvvals(1,ii),norder,npols,
c     1              sigvals(1,ii))
                call legetens_pols_2d(uvvals(1,ii),ndeg,type,
     1               sigvals(1,ii))
                qwts(ii) = rr*wts(i)
              enddo
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
            call fker(uvvals(1,jj),xt,dpars,
     1         zpars,ipars,fval)
            xkernvals(jj) = fval*qwts(jj)
          enddo


c
cc         subtract contribution of current quad
c          
          do isig=1,ksigpols
            cintall(isig) = cintall(isig)-cvals(isig,iquad)
            ctmp(isig) = 0
          enddo

c
cc        add in contributions of quads 
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
     1                                  sigvals(isig,ii)
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
cc        end of looping over all quads at current stage
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

      subroutine cquadints_adap_noprecomp(eps,intype,
     1     norder,ttype,npols,ntarg,xyztarg,nquadmax,
     3     fker,dpars,zpars,ipars,nqorder,cintvals)

c
c       this subroutine computes the integrals
c
c       \int_{[-1,1]^2} K(x_{i},y) P_{n}(y_{1}) P_{m}(y_{2}) dy \, ,
c
c        P_{n}(y) are Legendre polynomials on [-1,1]
c
c        using adaptive integration without any precomputation
c
c
c        input arguments:
c        eps:     requested precision 
c
c        intype:   quadrature node type
c                   intype = 1, tensor product gauss-legendre nodes
c                   intype = 2, xiao gimbutas nodes
c   
c        norder: order of polynomials on the patch 
c        npols = norder*norder, number of polynomials to be integrated 
c        ntarg - total number of target points
c        xyztarg(3,ntarg) - location of target points 
c        nqmax - max number of quads to be used on base quad
c        fker - function handle for evaluating the kernel k
c 
c               expected calling sequence
c               fker(x,y,ynorms,dpars,zpars,ipars,f)
c
c               x \in \mathbb{R}^{3}, y \in \mathbb{R}^{2}
c               
c               the output is assumed to be complex for the time
c               being
c
c         dpars(*) - real parameters for the fker routine
c         zpars(*) - complex parameters for the fker routine
c         ipars(*) - integer parameters for the fker routine
c         nqorder - order of quadrature nodes on each subquad
c                   to be used
c
c         output:
c
c         cintvals(npols,ntarg) - integrals at all targets
c                                  for all tensor product
c                                  chebyshev polynomials
c
c
c
      implicit none

c
cc     calling sequence variables
c
      real *8 eps
      integer intype
      integer norder,npols
      
      integer ntarg
      real *8 xyztarg(3,ntarg)
      
      external fker
      real *8 dpars(*)
      complex *16 zpars(*)
      integer ipars(*)

      integer nqorder

      integer nquadmax

      complex *16 cintvals(npols,ntarg)

      character ttype

c
c       tree variables
c
      integer nlmax,ltree
      real *8, allocatable :: tvs(:,:,:),da(:)
      integer, allocatable :: ichild_start(:)

      integer nquad,nlev,iquad,istart,i,j,k
      integer ier,itarg,jj,jstart,npts
      integer iqquad,ii


      integer npmax

      real *8, allocatable :: uvsq(:,:),wts(:),uvtmp(:,:)
      real *8, allocatable :: umattmp(:,:),vmattmp(:,:)
      integer nqpols
      real *8, allocatable :: sigvals(:,:)
      real *8, allocatable :: uvvals(:,:),qwts(:)
      integer itmp

      character *1 transa,transb
      real *8 alpha,beta,ra
      integer lda,ldb,ldc
      real *8 u, v
      integer ldu, ldv, itype,ndeg
      
c
c       for each quad, we just store three pieces
c       of info
c         quad vertices
c         area of quad
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
c
c
c        quad vertices nomenclature
c
c       v3
c        ________ 
c        |       |
c        |       |
c        |       |
c        ---------
c        v1       v2
c
c
c
      tvs(1,1,1) = -1
      tvs(2,1,1) = -1

      tvs(1,2,1) = 1
      tvs(2,2,1) = -1

      tvs(1,3,1) = -1
      tvs(2,3,1) = 1

c
c       get quadrature nodes and weights on the base quad
c       based on quadrature type
c

      if(intype.eq.1) then
        nqpols = nqorder*nqorder
        allocate(uvsq(2,nqpols),wts(nqpols))

        ldu = 1
        ldv = 1
        itype = 1
        call legetens_exps_2d(itype,nqorder,ttype,uvsq,u,ldu,v,ldv,wts)        
      endif

      if(intype.eq.2) then
        call squarearbq_pts(nqorder,nqpols)
        allocate(uvsq(2,nqpols),wts(nqpols))
        call squarearbq(nqorder,uvsq,wts,nqpols)
      endif


      allocate(uvtmp(2,nqpols))
     
      npmax = nquadmax*nqpols
      allocate(sigvals(npols,nqpols))
      allocate(uvvals(2,nqpols),qwts(nqpols))

c
c      current number of quads in the adaptive structure
c
      nquad = 1
c
c        intialize sigvals for root quad
c


c
cc       for the current patch compute geometry info for base quad 
c

       
      nlmax = 20 
      do itarg=1,ntarg
        nquad = 1
        ndeg = norder-1
        call mapuv_quad(tvs(1,1,1),nqpols,uvsq,uvvals)
        do i=1,nqpols
          call legetens_pols_2d(uvvals(1,i),ndeg,ttype,sigvals(1,i))
          qwts(i) = wts(i)
        enddo
        
        call quadadap_noprecomp(eps,nqorder,nqpols,nlmax,nquadmax,nquad,
     1    ichild_start,tvs,da,uvsq,wts, 
     1    norder,ttype,npols,npmax,uvvals,qwts,sigvals,xyztarg(1,itarg),
     3    fker,dpars,zpars,ipars,cintvals(1,itarg))
      enddo


      return
      end

c
c
c
c
c
      subroutine quadadap_noprecomp(eps,m,kpols,nlmax,nqmax,nquad,
     1             ichild_start,tvs,da,uvsq,wts,
     1             norder,ttype,npols,npmax,uvvals,qwts,
     2             sigvals,xt,fker,dpars,zpars,
     3             ipars,cintall)

c
c       this subroutine adaptively computes the integral
c        of the functions
c   
c        \int_{(-1,1)^2} 1/|xt- y(u,v)| P_{n}(u)P_{m}(v) du dv
c        n,m = 1,2\ldots npols
c
c        P_{n}(u) are the Legendre polynomials on (-1,1) 
c
c
c        relevant parameters on a quad refined
c        grid along the way
c
c        IN:
c        eps - precision requested
c        m - quadrature order
c        kpols - number of quarature nodes (look up squarearbq 
c                 for getting kpols(m))
c        nlmax - max level of refinement for geometry
c        nqmax - max number of quads
c        nquad - current number of quads in adaptive structure
c        ichild_start(i) - first child of quad i
c        tvs(2,3,nqmax) - vertices of hierarchy of quads
c        da(nqmax) - area of quads
c        uvsq(kpols) - integration nodes on standard quad
c        wts(kpols) - integration weights on standard quad
c        npols - total number of koornwinder polynomials to be integrated
c        norder - order of discretization of the surface
c        npols - (norder+1)*(norder+2)/2: total number of 
c                koornwinder polynomials to be integrated
c        npmax - max number of points = nqmax*kpols
c        uvvals(2,npmax) - geometry info on heirarchy of meshes
c        qwts(npmax) - quadrature weights 
c        sigvals(npols,npmax) - 
c                   tensor product GL polynomials computed along the adaptive grid
c        
c        OUT:
c        cintall(npols) - computed integral 
c
c         

      implicit real *8 (a-h,o-z)
      integer, allocatable :: istack(:)
      integer ichild_start(nqmax)
      real *8 da(nqmax)
      real *8 tvs(2,3,nqmax), uvsq(2,kpols),wts(kpols)
      integer nproclist0, nproclist
      integer idone
      real *8 sigvals(npols,kpols)
      real *8 uvvals(2,kpols),qwts(kpols)
      complex *16, allocatable :: xkernvals(:)
      real *8 xt(3),xs(2)
      complex *16 cintall(npols),fval,ctmp(npols)
      complex *16, allocatable :: cvals(:,:)

      real *8 dpars(*)
      complex *16 zpars(*)
      integer ipars(*)

      character ttype
      
      external fker

c
c         for historic reasons
c
      ksigpols = npols
      allocate(istack(2*nqmax))
      allocate(cvals(ksigpols,nqmax))
      nproclist0 = 1
      istack(1) = 1


      do i=1,ksigpols
         cvals(i,1) = 0
      enddo

      allocate(xkernvals(kpols))

c
cc      compute integral at level 0
c
      do i=1,kpols
         call fker(uvvals(1,i),xt,dpars,zpars,ipars,fval)
         xkernvals(i) = fval*qwts(i)
         do j=1,ksigpols
            cvals(j,1) = cvals(j,1)+xkernvals(i)*sigvals(j,i)
         enddo
      enddo

      
      do i=1,ksigpols
         cintall(i) = cvals(i,1)
      enddo


      call quadadap_main_noprecomp(eps,kpols,nlmax,nqmax,nquad,
     1      ichild_start,
     1      tvs,da,uvsq,wts,norder,ttype,npols,
     2      npmax,uvvals,qwts,sigvals,xt,fker,dpars,
     3      zpars,ipars,cvals,istack,nproclist0,
     4      xkernvals,cintall)

      
      return
      end
c
c
c
c
c
       
      subroutine quadadap_main_noprecomp(eps,kpols,nlmax,nqmax,nquad,
     1    ichild_start,tvs,da,uvsq,wts,norder,ttype,npols,
     2    npmax,uvvals,qwts,sigvals,xt,fker,dpars,
     3    zpars,ipars,cvals,istack,nproclist0,xkernvals,
     4    cintall)
      

      implicit real *8 (a-h,o-z)
      integer istack(*),nproclist0
      integer ichild_start(nqmax)
      real *8 da(nqmax)
      real *8 tvs(2,3,nqmax), uvsq(2,kpols),wts(kpols)
      integer  nproclist
      integer idone
      real *8 sigvals(npols,kpols)
      complex *16 xkernvals(kpols)
      real *8 xt(3)
      real *8 uvvals(2,kpols),qwts(kpols)
      complex *16 cintall(npols),fval,ctmp(npols)
      complex *16 cvals(npols,nqmax)

      character ttype

      real *8 dpars(*)
      complex *16 zpars(*)
      integer ipars(*)

      character *1 transa,transb
      integer lda,ldb,ldc
      external fker
      

c
c         for historic reasons
c
      ksigpols = npols
      kfine = 4*kpols


      do ilev=0,nlmax
        idone = 1
        nproclist = 0
        
        do iproc = 1,nproclist0
          iquad = istack(iproc)
          if(nquad+4.gt.nqmax) then
            print *, "Too many quads in cquadadap"
            print *, "Exiting before convergence"

            return
          endif
            
          ichild_start(iquad) = nquad+1
          call getquadchildren(tvs(1,1,iquad),tvs(1,1,nquad+1),
     1         tvs(1,1,nquad+2),tvs(1,1,nquad+3),tvs(1,1,nquad+4))
            

          ndeg = norder-1

c
cc         subtract contribution of current quad
c          
          do isig=1,ksigpols
            cintall(isig) = cintall(isig)-cvals(isig,iquad)
            ctmp(isig) = 0
          enddo

          rr = 0.25d0*da(iquad)
          do j=nquad+1,nquad+4
            da(j) = rr
            call mapuv_quad(tvs(1,1,j),kpols,uvsq,uvvals)
            do i=1,kpols
              call legetens_pols_2d(uvvals(1,i),ndeg,ttype,
     1               sigvals(1,i))
              call fker(uvvals(1,i),xt,dpars,zpars,ipars,fval)
              xkernvals(i) = fval*rr*wts(i)
            enddo
            do isig=1,ksigpols
              cvals(isig,j) = 0
            enddo

            do k=1,kpols
              do isig=1,ksigpols
                cvals(isig,j) = cvals(isig,j) + xkernvals(k)*
     1                      sigvals(isig,k)             
              enddo
            enddo
            do isig=1,ksigpols
              cintall(isig) = cintall(isig) + cvals(isig,j)
              ctmp(isig) = ctmp(isig) + cvals(isig,j)
            enddo
          enddo

          errmax = 0
          do isig=1,ksigpols
            if(abs(ctmp(isig)-cvals(isig,iquad)).gt.errmax) 
     1          errmax = abs(ctmp(isig)-cvals(isig,iquad))
          enddo

          if(errmax.gt.eps) then
            idone = 0
            do j=1,4
              istack(nproclist0+nproclist+j) = nquad+j
            enddo
            nproclist = nproclist+4
          endif
          nquad = nquad+4
cc       end of looping ovr all quads at current level          
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
c--------------------------------------------------------------------------------
        
      subroutine mapuv_quad(verts,kpols,uvs,uvout)
      implicit real *8 (a-h,o-z)
      integer kpols
      real *8 verts(2,3),uvs(2,kpols),uvout(2,kpols)

      dx = verts(1,2)-verts(1,1)
      dy = verts(2,3)-verts(2,1) 

      do i=1,kpols
        uvout(1,i) = verts(1,1) + dx*(uvs(1,i)+1)/2
        uvout(2,i) = verts(2,1) + dy*(uvs(2,i)+1)/2
      enddo

      return
      end
c-----------------------------------------      
