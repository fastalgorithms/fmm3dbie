c
c
c    generate level restricted quad tree based on resolving a 
c    function to desired precision
c
c
c    The function handle is of the form
c    call fun(nd,xy,dpars,zpars,ipars,f)
c
c      where xy is the location in (-L/2,L/2)^2
c
c    dpars are a collection of real parameters, zpars are complex
c    parameters and ipars are integer parameters, and f is a 
c    real array of size nd
c     
c
c    For the Helmholtz/Maxwell tree, the boxes are refined until 
c     Re(zk)*boxsize<5, beyond which the function resolution criterion
c    kicks in. 
c
c    A function is said to be resolved if its interpolant at the 8
c    child nodes, agrees with the function values at those nodes
c    up to the user specified tolerance. 
c    The error is scaled by h**(eta)
c    where eta is user specified and h is the boxsize. If there is 
c    any confusion, the user should seet \eta to 0
c    Let \tilde{f} denote the interpolant of f, then
c    the refinement criterion is 
c      \int_{B_{j} |\tilde{f}-f|^{p} *h^{\eta} < 
c        \varepsilon V_{j}^{1/p}/V_{0}^{1/p}/(\int_{B_{0}}|f|^{p})^{1/p}
c    This implies that
c       \int_{B_{0}} |\tilde{f}-f|^{p} =
c          \sum_{j} \int_{B_{j}} |\tilde{f}-f|^{p}
c          \leq \sum_{j} \eps^{p}*h^{\eta p} 
c                  V_{j}/V_{0}/(\int_{B_{0}} |f|^p)
c      If \eta = 0,
c          \leq \eps^{p}/(\int_{B_{0}} |f|^{p})
c
c    i.e., this strategy guarantees that the interpolated function
c      approximates the function with relative lp accuracy of \eps
c      
c    This code has 2 main user callable routines
c      make_vol_tree_mem -> Returns the memory requirements, 
c          tree length, number of boxes, number of levels
c      make_vol_tree -> Makes the actual tree, returns centers of boxes,
c          colleague info, function values on leaf boxes
c       
c          
c          iptr(1) - laddr
c          iptr(2) - ilevel
c          iptr(3) - iparent
c          iptr(4) - nchild
c          iptr(5) - ichild
c          iptr(6) - ncoll
c          iptr(7) - coll
c          iptr(8) - ltree
c




      subroutine surf_tree_mem(nsplit,eps,zk,boxlen,norder,iptype,
     1     eta,fun,nd,dpars,zpars,ipars,nlevels,nboxes,ltree,rintl,
     2     nboxesleaf)
c
c      get memory requirements for the tree
c
c     input parameters:
c        nsplit - integer(2)
c           nsplit(1): split [-boxlen/2,boxlen/2]^2 into nsplit(1) subdivisions at root level in the x-direction,
c           nsplit(1) = 1 means no subdivisions in the x-direction
c           nsplit(2): split [-boxlen/2,boxlen/2]^2 into nsplit(2) subdivisions at root level in the y-direction
c           nsplit(2) = 1 means no subdivisions in the y-direction
c        eps - double precision
c           precision requested
c        zk - double complex
c           Helmholtz parameter
c        boxlen - double precision
c           length of box in which volume is contained, 
c           if boxlen = L, then volume is [-L/2,L/2]^2
c        norder - integer
c           order of discretization
c        eta - double precision
c           scaling parameter for error
c        fun - function handle
c           function to evalute it everywhere in the volume
c        nd - integer
c           number of real functions returned by fun
c        dpars - double precision
c           real parameters for function evaluation
c        zpars - double complex
c           complex parameters for function evaluation
c        ipars - integer
c           integer parameters for function evaluation
c 
c        output:
c           nlevels - integer
c             number of levels
c           nboxes - integer
c             number of boxes
c           ltree - integer
c             length of tree
c           rintl(0:nlevels) - real *8
c             lp norm to scale the functions by
c             (on input rintl should be of size(0:200)
c           nboxesleaf - integer
c             number of leaf boxes      
c
c      
c

      implicit none
      real *8 eps,boxlen,eta,dpars(*),delta
      real *8, allocatable :: fvals(:,:,:)
      complex *16 zpars(*),zk
      integer nd,ipars(*),iptype,nsplit(2),nroot,ipoly
      integer nlevels,nboxes,ltree,norder,nboxesleaf

      external fun

      integer, allocatable :: laddr(:,:),ilevel(:),iparent(:),nchild(:)
      integer, allocatable :: ichild(:,:),ncoll(:),icoll(:,:)
      real *8, allocatable :: centers(:,:)
      integer, allocatable :: nbors(:,:),nnbors(:)

      integer, allocatable :: ilevel2(:),iparent2(:),nchild2(:),
     1    ichild2(:,:)
      real *8, allocatable :: centers2(:,:),fvals2(:,:,:)

      integer nbmax,nlmax,npbox,npc
      real *8, allocatable :: grid(:,:),qwts(:)
      real *8, allocatable:: xq(:),wts(:),umat(:,:),vmat(:,:)
      real *8 xy(2)
      real *8, allocatable :: wts2(:),xref2(:,:),umat2(:,:)
      real *8 rintl(0:200)
      real *8 rint
      real *8, allocatable :: rintbs(:),rintbs2(:)
      integer i,itype,j,npols,k

      real *8, allocatable :: fval1(:,:,:,:),centerstmp(:,:,:)
      real *8, allocatable :: boxsize(:,:)
      integer, allocatable :: irefinebox(:)

      real *8 rsc,ra
      integer nbloc,nbctr,nbadd,irefine,ilev,ifirstbox,ilastbox
      integer nbtot,iii,idim,iper
      

      nbmax = 100000
      nlmax = 200

      allocate(boxsize(2,0:nlmax))

      
      allocate(laddr(2,0:nlmax),ilevel(nbmax),iparent(nbmax))
      allocate(nchild(nbmax),ichild(4,nbmax))

      allocate(fvals(nd,norder**2,nbmax),centers(2,nbmax))

      allocate(rintbs(nbmax))

c
c      set tree info for level 0
c
      nroot = nsplit(1)*nsplit(2)
      laddr(1,0) = 1
      laddr(2,0) = nroot
      
      do j=1,nroot
         ilevel(j) = 0
         iparent(j) = -1
         nchild(j) = 0
         rintbs(j) = 0
         do i=1,4
            ichild(i,j) = -1
         enddo
      enddo
c
c
c      

      if(nroot.eq.1) then
         boxsize(1,0) = boxlen
         boxsize(2,0) = boxlen
         centers(1,1) = 0
         centers(2,1) = 0
      else if(nroot.gt.1) then
         boxsize(1,0) = boxlen/nsplit(1)
         boxsize(2,0) = boxlen/nsplit(2)
      endif

      k = 1
      do i=1,nsplit(1)    
         do j=1,nsplit(2)
            centers(1,k)=-boxlen/2.0d0
     1        + boxlen/(nsplit(1))/2.0d0
     2        + boxlen/(nsplit(1))*(i-1)     
            centers(2,k)=-boxlen/2.0d0
     1           + boxlen/(nsplit(2))/2.0d0
     2           + boxlen/(nsplit(2))*(j-1)           
            k = k + 1           
         enddo

      enddo
c
c
      npbox = norder**2
      allocate(grid(2,npbox))
c
c     Generate a grid on the box [-1/2,1/2]^2
c
      allocate(xq(norder),wts(norder))
      allocate(umat(norder,norder),vmat(norder,norder))
      
c
      itype = 2
      call chebexps(itype,norder,xq,umat,vmat,wts)
      do i=1,norder
        xq(i) = xq(i)/2
      enddo

      call surf_mesh2d(xq,norder,xq,norder,grid)


      npols = norder**2
      allocate(wts2(npbox),xref2(2,npbox))
      ipoly = 1
      itype = 1
      call polytens_exps_2d(ipoly,itype,norder,'f',xref2,umat2,
     1     npols,vmat,1,wts2)
c      call chebtens_exps_2d(itype,norder,'f',xref2,umat2,npols,
c     1   vmat,1,wts2)
      
c
c       compute fvals at the grid
c

      rint = 0


c
c   note extra factor of 4 since wts2 are on [-1,1]^2 
c   as opposed to [-1/2,1/2]^2
c
      rsc = boxlen**2/4
      rsc = boxsize(1,0)*boxsize(2,0)/(4.0)
      do k=1,nroot

         do i=1,npbox            
            xy(1) = centers(1,k) + grid(1,i)*boxsize(1,0)
            xy(2) = centers(2,k) + grid(2,i)*boxsize(2,0)
            call fun(nd,xy,dpars,zpars,ipars,fvals(1,i,k))
            if(iptype.eq.0) then
               do idim=1,nd
                  if(abs(fvals(idim,i,1)).gt.rintbs(1)) rintbs(1) = 
     1                 abs(fvals(idim,i,1))
               enddo
            endif

            if(iptype.eq.1) then
               do idim=1,nd
                  rintbs(1) = rintbs(1)+abs(fvals(idim,i,1))*wts2(i)*rsc
               enddo
            endif

            if(iptype.eq.2) then
               do idim=1,nd
                  rintbs(1) = rintbs(1) + fvals(idim,i,1)**2*wts2(i)*rsc
               enddo
            endif
         enddo
      enddo

      if(iptype.eq.0.or.iptype.eq.1) rint = rintbs(1)
      if(iptype.eq.2) rint = sqrt(rintbs(1))

      rintl(0) = rint

      nbctr = nroot
     
      
      do ilev=0,nlmax-1
        irefine = 0

        ifirstbox = laddr(1,ilev) 
        ilastbox = laddr(2,ilev)

        nbloc = ilastbox-ifirstbox+1

        allocate(fval1(nd,norder**2,4,nbloc))
        allocate(centerstmp(2,4,nbloc))
        allocate(irefinebox(nbloc))

        
        if(iptype.eq.2) rsc = sqrt(1.0d0/boxlen**2)
        if(iptype.eq.1) rsc = (1.0d0/boxlen**2)
        if(iptype.eq.0) rsc = 1.0d0

        rsc = rsc*rint

        print *, ilev, rint,rsc

        call surf_tree_find_box_refine(nd,iptype,eta,eps,zk,norder,
     1       npbox,fvals,npols,umat,boxsize(1,ilev),nbmax,ifirstbox,
     2       nbloc,rsc,irefinebox,irefine)
     
c
c
c          figure out if current set of boxes is sufficient
c

        nbadd = 0 
        do i=1,nbloc
          if(irefinebox(i).eq.1) nbadd = nbadd+4
        enddo

        nbtot = nbctr+nbadd

c
c         if current memory is not sufficient reallocate
c
        if(nbtot.gt.nbmax) then
          print *, "Reallocating"
          allocate(centers2(2,nbmax),ilevel2(nbmax),iparent2(nbmax))
          allocate(nchild2(nbmax),ichild2(4,nbmax))
          allocate(fvals2(nd,npbox,nbmax),rintbs2(nbmax))

          call surf_tree_copy(nd,nbctr,npbox,centers,ilevel,iparent,
     1            nchild,ichild,fvals,centers2,ilevel2,iparent2,nchild2,
     2            ichild2,fvals2)
          call dcopy(nbctr,rintbs,1,rintbs2,1)

          deallocate(centers,ilevel,iparent,nchild,ichild,fvals,rintbs)

          nbmax = nbtot
          allocate(centers(2,nbmax),ilevel(nbmax),iparent(nbmax))
          allocate(nchild(nbmax),ichild(4,nbmax),fvals(nd,npbox,nbmax))
          allocate(rintbs(nbmax))

          call surf_tree_copy(nd,nbctr,npbox,centers2,ilevel2,iparent2,
     1            nchild2,ichild2,fvals2,centers,ilevel,iparent,nchild,
     2            ichild,fvals)
          call dcopy(nbctr,rintbs2,1,rintbs,1)

          deallocate(centers2,ilevel2,iparent2,nchild2,ichild2,fvals2)
          deallocate(rintbs2)
        endif


        if(irefine.eq.1) then
           boxsize(1,ilev+1) = boxsize(1,ilev)/2
           boxsize(2,ilev+1) = boxsize(2,ilev)/2
           laddr(1,ilev+1) = nbctr+1

          call surf_tree_refine_boxes(irefinebox,nd,npbox,fvals,
     1      fun,dpars,zpars,ipars,grid,nbmax,ifirstbox,nbloc,centers,
     2      boxsize(1,ilev+1),nbctr,ilev+1,ilevel,iparent,nchild,ichild)


          rsc = boxsize(1,ilev+1)*boxsize(2,ilev+1)/4
          call surf_update_rints(nd,npbox,nbmax,fvals,ifirstbox,nbloc,
     1       iptype,nchild,ichild,wts2,rsc,rintbs,rint)
          
          rintl(ilev+1) = rint

          
          laddr(2,ilev+1) = nbctr
        else
          exit
        endif

        deallocate(fval1,irefinebox,centerstmp)
      enddo

      nboxes = nbctr
      nlevels = ilev

      if(nlevels.ge.2) then

        nbtot = 8*nboxes
        if(nbtot.gt.nbmax) then
          print *, "Reallocating"
          allocate(centers2(2,nbmax),ilevel2(nbmax),iparent2(nbmax))
          allocate(nchild2(nbmax),ichild2(4,nbmax))
          allocate(fvals2(nd,npbox,nbmax),rintbs2(nbmax))

          call surf_tree_copy(nd,nboxes,npbox,centers,ilevel,iparent,
     1       nchild,ichild,fvals,centers2,ilevel2,iparent2,nchild2,
     2       ichild2,fvals2)
          call dcopy(nboxes,rintbs,1,rintbs2,1)

          deallocate(centers,ilevel,iparent,nchild,ichild,fvals,rintbs)

          nbmax = nbtot
          allocate(centers(2,nbmax),ilevel(nbmax),iparent(nbmax))
          allocate(nchild(nbmax),ichild(4,nbmax),fvals(nd,npbox,nbmax))
          allocate(rintbs(nbmax))

          call surf_tree_copy(nd,nboxes,npbox,centers2,ilevel2,iparent2,
     1          nchild2,ichild2,fvals2,centers,ilevel,iparent,nchild,
     2          ichild,fvals)
          call dcopy(nboxes,rintbs2,1,rintbs,1)

          deallocate(centers2,ilevel2,iparent2,nchild2,ichild2,fvals2)
          deallocate(rintbs2)
        endif

        allocate(nnbors(nbmax))
        allocate(nbors(9,nbmax))

        do i=1,nboxes
          nnbors(i) = 0
          do j=1,9
            nbors(j,i) = -1
          enddo
        enddo

        iper = 0

        call surf_computecoll2d(nlevels,nboxes,laddr,boxsize,centers,
     1        iparent,nchild,ichild,iper,nnbors,nbors)

        if(nlevels.ge.2) then
c         call vol_tree_fix_lr(fun,nd,dpars,zpars,ipars,norder,npbox,
c     1         fvals,grid,centers,nlevels,nboxes,boxsize,nbmax,nlmax,
c     2         laddr,ilevel,iparent,nchild,ichild,nnbors,nbors)
        endif
      endif

      ltree = 17*nboxes + 2*(nlevels+1)

      nboxesleaf = 0
      do i =1,nboxes
         if(nchild(i).eq.0) then
c     is leaf
            nboxesleaf = nboxesleaf + 1         
                  
         endif
      end do

      return
      end
c
c
c
c
c

      subroutine surf_tree_build(nsplit,eps,zk,boxlen,norder,iptype,eta,
     1  fun,nd,dpars,zpars,ipars,nlevels,nboxes,ltree,rintl,itree,iptr,
     2  fvals,centers,boxsize)
c
c      construct the tree
c
c     input parameters:
c        nsplit - integer(2)
c           nsplit(1): split [-boxlen/2,boxlen/2]^2 into nsplit(1) subdivisions at root level in the x-direction,
c           nsplit(1) = 1 means no subdivisions in the x-direction
c           nsplit(2): split [-boxlen/2,boxlen/2]^2 into nsplit(2) subdivisions at root level in the y-direction
c           nsplit(2) = 1 means no subdivisions in the y-direction
c        eps - double precision
c           precision requested
c        zk - double complex
c           Helmholtz parameter
c        boxlen - double precision
c           length of box in which volume is contained, 
c           if boxlen = L, then volume is [-L/2,L/2]^2
c        norder - integer
c           order of discretization
c        iptype - integer
c           error norm
c           iptype = 0 - linf
c           iptype = 1 - l1
c           iptype = 2 - l2
c        eta - double precision
c           scaling parameter for error
c        fun - function handle
c           function to evalute it everywhere in the volume
c        nd - integer
c           number of real functions returned by fun
c        dpars - double precision
c           real parameters for function evaluation
c        zpars - double complex
c           complex parameters for function evaluation
c        ipars - integer
c           integer parameters for function evaluation
c        nlevels - integer
c          number of levels
c        nboxes - integer
c          number of boxes
c        ltree - integer
c          length of tree = 2*(nlevels+1)+17*nboxes
c        rintl - real *8 (0:nlevels)
c          estimate of lp norm for scaling the errors
c          at various levels. 
c          We require the estimate at each level to make sure
c          that the memory estimate code is consitent
c          with the build code else there could be potential
c          memory issues 
c         
c
c      output:
c        itree - integer(ltree)
c          tree info
c        iptr - integer(8)
c          iptr(1) - laddr
c          iptr(2) - ilevel
c          iptr(3) - iparent
c          iptr(4) - nchild
c          iptr(5) - ichild
c          iptr(6) - ncoll
c          iptr(7) - coll
c          iptr(8) - ltree
c        fvals - double precision (nd,norder**2,nboxes)
c          function values at discretization nodes
c        centers - double precision (2,nboxes)
c          xyz coordinates of box centers in the quad tree
c        boxsize - double precision (0:nlevels)
c          size of box at each of the levels
c

      implicit none
      real *8 eps,boxlen,eta,dpars(*)
      complex *16 zk,zpars(*)
      integer nd,ipars(*),iptype
      integer nlevels,nboxes,ltree,norder,nroot
      integer iptr(8),nsplit(2)
      integer itree(ltree),ier
      real *8 fvals(nd,norder**2,nboxes),centers(2,nboxes)
      real *8, allocatable :: fval1(:,:,:,:),centerstmp(:,:,:)
      integer, allocatable :: irefinebox(:)
      real *8 boxsize(2,0:nlevels)
c      real *8 xq(norder),wts(norder),umat(norder,norder)
c      real *8 vmat(norder,norder)
      real *8, allocatable :: xq(:),wts(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: grid(:,:),ximat(:,:),qwts(:)
      real *8 rintl(0:nlevels)
      real *8 xy(2)

      integer i,ilev,irefine,itype,nbmax,nlmax,npbox,npc,ii,k
      integer ifirstbox,ilastbox,nbctr,nbloc
      real *8 rsc

      real *8 ra
      integer j,nboxes0,npols,iper

      external fun

      allocate(xq(norder),wts(norder))
      allocate(umat(norder,norder),vmat(norder,norder))
c
      iptr(1) = 1
      iptr(2) = 2*(nlevels+1)+1
      iptr(3) = iptr(2) + nboxes
      iptr(4) = iptr(3) + nboxes
      iptr(5) = iptr(4) + nboxes
      iptr(6) = iptr(5) + 4*nboxes
      iptr(7) = iptr(6) + nboxes
      iptr(8) = iptr(7) + 9*nboxes



      nroot = nsplit(1)*nsplit(2)
      itree(1) = 1
      itree(2) = nroot

c
c      set tree info for level 0
c

      
      do j=1,nroot
         itree(iptr(2)+j-1) = 0
         itree(iptr(3)+j-1) = -1
         itree(iptr(4)+j-1) = 0
         do i=1,4
            itree(iptr(5)+(j-1)*4 +i-1) = -1
         enddo
      enddo

      if(nroot.eq.1) then
         boxsize(1,0) = boxlen
         boxsize(2,0) = boxlen
         centers(1,1) = 0
         centers(2,1) = 0
      else if(nroot.gt.1) then
         boxsize(1,0) = boxlen/nsplit(1)
         boxsize(2,0) = boxlen/nsplit(2)
      endif


      k = 1
      do i=1,nsplit(1)    
         do j=1,nsplit(2)
            centers(1,k)=-boxlen/2.0d0
     1           + boxlen/(nsplit(1))/2.0d0
     2           + boxlen/(nsplit(1))*(i-1)     
            centers(2,k)=-boxlen/2.0d0
     1           + boxlen/(nsplit(2))/2.0d0
     2           + boxlen/(nsplit(2))*(j-1)
            k = k + 1
         enddo
      enddo
      
c
c
      npbox = norder**2
      allocate(grid(2,npbox))
c
c     Generate a grid on the box [-1/2,1/2]^2
c
c
      itype = 2
      call chebexps(itype,norder,xq,umat,vmat,wts)
      do i=1,norder
        xq(i) = xq(i)/2
      enddo

      call surf_mesh2d(xq,norder,xq,norder,grid)
      npols = norder**2
      
      do k=1,nroot
         do i=1,npbox           
            xy(1) = centers(1,k) + grid(1,i)*boxsize(1,0)
            xy(2) = centers(2,k) + grid(2,i)*boxsize(2,0)
            call fun(nd,xy,dpars,zpars,ipars,fvals(1,i,k))
         enddo
      enddo
      

c
c       Reset nlevels, nboxes
c     
      nbctr = nroot
      
      do ilev=0,nlevels-1
         irefine = 0
         
         ifirstbox = itree(2*ilev+1) 
         ilastbox = itree(2*ilev+2)
         
         nbloc = ilastbox-ifirstbox+1

         allocate(fval1(nd,norder**2,4,nbloc))
         allocate(centerstmp(2,4,nbloc))
         allocate(irefinebox(nbloc))

         
         if(iptype.eq.2) rsc = sqrt(1.0d0/boxlen**2)
         if(iptype.eq.1) rsc = (1.0d0/boxlen**2)
         if(iptype.eq.0) rsc = 1.0d0
         rsc = rsc*rintl(ilev)
         call surf_tree_find_box_refine(nd,iptype,eta,eps,zk,norder,
     1        npbox,fvals,npols,umat,boxsize(1,ilev),nboxes,ifirstbox,
     2        nbloc,rsc,irefinebox,irefine)
         

         if(irefine.eq.1) then
            boxsize(1,ilev+1) = boxsize(1,ilev)/2
            boxsize(2,ilev+1) = boxsize(2,ilev)/2
            itree(2*ilev+3) = nbctr+1
            call surf_tree_refine_boxes(irefinebox,nd,npbox,fvals,
     1           fun,dpars,zpars,ipars,grid,nboxes,ifirstbox,nbloc,
     2           centers,boxsize(1,ilev+1),nbctr,ilev+1,itree(iptr(2)),
     3           itree(iptr(3)),itree(iptr(4)),itree(iptr(5)))
            itree(2*ilev+4) = nbctr
         else
            exit
         endif

         deallocate(fval1,irefinebox,centerstmp)
      enddo
      nboxes0 = nbctr
      nlevels = ilev

      do i=1,nboxes0
         itree(iptr(6)+i-1) = 0
         do j=1,9
            itree(iptr(7)+9*(i-1)+j-1) = -1
         enddo
      enddo

      iper = 0

      call surf_computecoll2d(nlevels,nboxes0,itree(iptr(1)),boxsize,
     1        centers,itree(iptr(3)),itree(iptr(4)),itree(iptr(5)),iper,
     2        itree(iptr(6)),itree(iptr(7)))

      if(nlevels.ge.2) then
c         call vol_tree_fix_lr(fun,nd,dpars,zpars,ipars,norder,npbox,
c     1       fvals,grid,centers,nlevels,nboxes0,boxsize,nboxes,nlevels,
c     2       itree(iptr(1)),itree(iptr(2)),itree(iptr(3)),
c     3       itree(iptr(4)),itree(iptr(5)),itree(iptr(6)),
c     4       itree(iptr(7)))
      endif

cccc      call prinf('nboxes0=*',nboxes0,1)
cccc      call prinf('nlevels=*',nlevels,1)
      return
      end
c
c
c      
c
c
      subroutine surf_tree_find_box_refine(nd,iptype,eta,eps,zk,norder,
     1  npbox,fvals,npols,umat,
     2  boxsize,nboxes,ifirstbox,nbloc,rsc,irefinebox,irefine)
      implicit none
      integer nd,npc,npbox,norder,iptype,npols
      integer nboxes,nbloc
      real *8 eta,eps,fvals(nd,npbox,nboxes)
      real *8 umat(norder,norder)
      real *8 rsc
      real *8, allocatable :: fcoefs(:,:,:),rmask(:)
      integer, allocatable :: iind2p(:,:)
      real *8 alpha,beta,boxsize(2),rsum
      integer irefinebox(nbloc),xind(4),yind(4)
      complex *16 zk
      integer ifirstbox

      integer irefine

      integer i,j,k,l,ibox,ifunif,i1
      real *8 rscale2,err,bs,bs2
      data xind/-1,1,-1,1/
      data yind/-1,-1,1,1/

      character *1 transa,transb

      external fun

      ifunif = 0

      transa = 'n'
      transb = 't'

      
      allocate(fcoefs(nd,npols,nbloc))
      allocate(iind2p(2,npols))
      call polytens_ind2pow_2d(norder-1,'f',iind2p)

      
      allocate(rmask(npols))

      rsum = 0
      do i=1,npols
        rmask(i) = 0.0d0
        i1=iind2p(1,i)+iind2p(2,i)
        if(i1.eq.norder-1) then
          rmask(i) = 1.0d0
          rsum = rsum + 1
        endif
      enddo

      if(iptype.eq.2) rsum = sqrt(rsum)
      if(iptype.eq.0) rsum = 1
      



      alpha = 1
      beta = 0

c      bs = boxsize/4.0d0     
      bs = boxsize(1)*boxsize(2)/4
      bs2 = 2*bs
      rscale2 = bs2**eta

      if(real(zk)*boxsize(1).gt.5) then
         do i=1,nbloc
            irefinebox(i) = 1
         enddo
         goto 1000
      endif

      if(real(zk)*boxsize(2).gt.5) then
         do i=1,nbloc
            irefinebox(i) = 1
         enddo
         goto 1000
      endif


C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ibox,err)
      do i=1,nbloc
        irefinebox(i) = 0

        ibox = ifirstbox + i-1

        call legval2coefs_2d(nd,norder,fvals(1,1,ibox),
     1      fcoefs(1,1,i),umat)

        call surf_fun_err(nd,npols,fcoefs(1,1,i),rmask,
     1     iptype,rscale2,err)

        err = err/rsum
		  
        if(err.gt.eps*rsc) then
          irefinebox(i) = 1
        endif
      enddo
C$OMP END PARALLEL DO     

 1000 continue
      irefine = maxval(irefinebox(1:nbloc))


c
c       make tree uniform
c

      if(ifunif.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,nbloc
          irefinebox(i) = irefine
        enddo
C$OMP END PARALLEL DO      
      endif

      return
      end
c
c
c
c
c


      subroutine surf_tree_refine_boxes(irefinebox,nd,npbox,fvals,
     1  fun,dpars,zpars,ipars,grid,nboxes,ifirstbox,nbloc,centers,
     2  bs,nbctr,nlctr,ilevel,iparent,nchild,ichild)
      implicit none
      integer nd,npbox
      integer ipars(*)
      real *8 dpars(*)
      complex *16 zpars(*)
      real *8 fvals(nd,npbox,nboxes)
      integer nboxes,nbloc,nbctr,nlctr
      real *8 grid(2,npbox),bs(2),xy(2)
      real *8 centers(2,nboxes)
      integer ilevel(nboxes),iparent(nboxes)
      integer ichild(4,nboxes),nchild(nboxes)
      integer irefinebox(nbloc)
      integer ifirstbox
      integer, allocatable :: isum(:)

      integer i,ibox,nel0,j,l,jbox,nel1,nbl
      integer xind(4),yind(4)

      real *8 bshx,bshy
      data xind/-1,1,-1,1/
      data yind/-1,-1,1,1/

      external fun

  
      allocate(isum(nbloc))
      call cumsum(nbloc,irefinebox,isum)

c      bsh = bs/2
      bshx = bs(1)/2
      bshy = bs(2)/2

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ibox,nbl,j,jbox,l,xy)
      do i = 1,nbloc
        ibox = ifirstbox + i-1
        if(irefinebox(i).eq.1) then
          nbl = nbctr + (isum(i)-1)*4
          nchild(ibox) = 4
          do j=1,4
             jbox = nbl+j             
             centers(1,jbox) = centers(1,ibox)+xind(j)*bshx
             centers(2,jbox) = centers(2,ibox)+yind(j)*bshy
            do l=1,npbox
              xy(1) = centers(1,jbox) + grid(1,l)*bs(1)
              xy(2) = centers(2,jbox) + grid(2,l)*bs(2)
              call fun(nd,xy,dpars,zpars,ipars,fvals(1,l,jbox))
            enddo
            iparent(jbox) = ibox
            nchild(jbox) = 0
            do l=1,4
              ichild(l,jbox) = -1
            enddo
            ichild(j,ibox) = jbox
            ilevel(jbox) = nlctr 
          enddo
        endif
      enddo
C$OMP END PARALLEL DO      

      nbctr = nbctr + isum(nbloc)*4


      return
      end
c
c
c
c
c
      subroutine surf_fun_err(nd,n,fcoefs,rmask,iptype,rscale,err)
c       this subroutine estimates the error based on the expansion
c       coefficients in a given basis
c       
c       input
c        nd - integer
c          number of functions
c        n -  integer
c           number of points at which function is tabulated
c        fcoefs: double precision(nd,n) 
c          tensor product legendre coeffs of func
c        rmask: double precision(n)
c           coefficients which will be accounted for in the error
c        iptype: integer
c           type of error to be computed
c           iptype = 0 - linf error
c           iptype = 1 - l1 error
c           iptype = 2 - l2 error
c        rscale: double precision
c          scaling factor
c
c       output
c         err: double precision
c           max scaled error in the functions
c
        implicit none
        integer n,i,iptype,idim,nd
        real *8 rscale,err
        real *8 fcoefs(nd,n),rmask(n)
        real *8, allocatable :: errtmp(:),ftmp(:,:)
        real *8 alpha, beta

        allocate(errtmp(nd),ftmp(nd,n))

        alpha = 1.0d0
        beta = 0.0d0
   
        err = 0
        do idim=1,nd
          errtmp(idim) = 0
        enddo
        if(iptype.eq.0) then
          do i=1,n
            if(rmask(i).gt.0.5d0) then
              do idim = 1,nd
                if(errtmp(idim).lt.abs(fcoefs(idim,i))) 
     1             errtmp(idim)=abs(fcoefs(idim,i))
              enddo
            endif
          enddo
        endif

        if(iptype.eq.1) then
          do i=1,n 
            do idim=1,nd
              ftmp(idim,i) = abs(fcoefs(idim,i))
            enddo
          enddo
          call dgemv('n',nd,n,alpha,ftmp,nd,rmask,1,beta,errtmp,1)
        endif


        if(iptype.eq.2) then
          do i=1,n 
            do idim=1,nd
              ftmp(idim,i) = fcoefs(idim,i)**2
            enddo
          enddo
          call dgemv('n',nd,n,alpha,ftmp,nd,rmask,1,beta,errtmp,1)

          do idim=1,nd
            errtmp(idim) = sqrt(errtmp(idim))
          enddo
        endif

        err = 0
        do idim=1,nd
          if(errtmp(idim).gt.err) err = errtmp(idim)
        enddo

        err = err*rscale

        return
        end

c
c
c
c       
c
      subroutine surf_update_rints(nd,npbox,nbmax,fvals,ifirstbox,nbloc,
     1       iptype,nchild,ichild,wts,rsc,rintbs,rint)
c
c------------------------
c  This subroutine updates the integrals of the function to
c  be resolved on the computational domain. It subtracts the
c  integral of the boxes which have been refined and adds
c  in the integrals corresponding to the function values 
c  tabulated at the children
c
c  Input arguments:
c  
c    - nd: integer
c        number of functions
c    - npbox: integer
c        number of points per box where the function is tabulated
c    - nbmax: integer
c        max number of boxes
c    - fvals: real *8 (nd,npbox,nbmax)
c        tabulated function values
c    - ifirstbox: integer
c        first box in the list of boxes to be processed
c    - nbloc: integer
c        number of boxes to be processed
c    - iptype: integer
c        Lp version of the scheme
c        * iptype = 0, linf
c        * iptype = 1, l1
c        * iptype = 2, l2
c    - nchild: integer(nbmax)
c        number of children 
c    - ichild: integer(8,nbmax)
c        list of children
c    - wts: real *8 (npbox)
c        quadrature weights for intgegrating functions 
c    - rsc: real *8
c        scaling parameter for computing integrals
c  
c  Inout arguemnts:
c
c     - rintbs: real *8(nbmax)
c         the integral for the new boxes cretated will be updated
c     - rint: real *8
c         the total integral will be updated
c    
c  
c      
      implicit real *8 (a-h,o-z)
      integer, intent(in) :: nd,npbox,nbmax
      real *8, intent(in) :: fvals(nd,npbox,nbmax)
      integer, intent(in) :: ifirstbox,nbloc,iptype
      integer, intent(in) :: nchild(nbmax),ichild(4,nbmax)
      real *8, intent(in) :: wts(npbox),rsc
      real *8, intent(inout) :: rintbs(nbmax),rint

c
c
c      compute the integrals for the newly formed boxes
c   and update the overall integral
c
      if(iptype.eq.0) then
         do i=1,nbloc  
            ibox = ifirstbox+i-1
            if(nchild(ibox).gt.0) then
               do j=1,4
                  jbox = ichild(j,ibox)
                  rintbs(jbox) = maxval(fvals(1:nd,1:npbox,jbox))
                  if(rintbs(jbox).gt.rint) rint = rintbs(jbox)
               enddo
            endif
         enddo
      endif

      if(iptype.eq.1) then
        do i=1,nbloc
          ibox=ifirstbox+i-1
          if(nchild(ibox).gt.0) then
c     subtract contribution of ibox from rint
            rint = rint - rintbs(ibox) 
          endif
        enddo

c
c     add back contribution of children
c
        do i=1,nbloc
          ibox = ifirstbox+i-1
          if(nchild(ibox).gt.0) then
            do j=1,4
              jbox = ichild(j,ibox)
              rintbs(jbox) = 0
              do l=1,npbox
                do idim=1,nd
                  rintbs(jbox) = rintbs(jbox) + 
     1               abs(fvals(idim,l,jbox))*wts(l)*rsc
                enddo
              enddo
              rint = rint + rintbs(jbox)
            enddo
          endif
        enddo
      endif

      if(iptype.eq.2) then
        rintsq = rint**2
        do i=1,nbloc
          ibox=ifirstbox+i-1
          if(nchild(ibox).gt.0) then
c
c    note that if iptype = 2, then rintbs stores squares
c    of the integral on the box
c
             rintsq = rintsq - rintbs(ibox)
          endif
        enddo

        do i=1,nbloc
          ibox = ifirstbox+i-1
          if(nchild(ibox).gt.0) then
            do j=1,4
              jbox = ichild(j,ibox)
              rintbs(jbox) = 0
              do l=1,npbox
                do idim=1,nd
                  rintbs(jbox) = rintbs(jbox) + 
     1               fvals(idim,l,jbox)**2*wts(l)*rsc
                enddo
              enddo
              rintsq = rintsq + rintbs(jbox)
            enddo
          endif
        enddo
        rint = sqrt(rintsq)
      endif
          

      return
      end
c
c
c
c
c
       subroutine surf_tree_copy(nd,nb,npb,centers,ilevel,iparent,
     1            nchild,ichild,fvals,centers2,ilevel2,iparent2,nchild2,
     2            ichild2,fvals2)

       implicit none
       integer nd,nb,npb
       real *8 centers(2,nb),centers2(2,nb)
       integer ilevel(nb),ilevel2(nb)
       integer iparent(nb),iparent2(nb)
       integer nchild(nb),nchild2(nb)
       integer ichild(4,nb),ichild2(4,nb)
       real *8 fvals(nd,npb,nb),fvals2(nd,npb,nb)

       integer i,j,nel

       nel = nd*npb*nb
       call dcopy(nel,fvals,1,fvals2,1)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
       do i=1,nb
         centers2(1,i) = centers(1,i)
         centers2(2,i) = centers(2,i)
         ilevel2(i) = ilevel(i)
         iparent2(i) = iparent(i)
         nchild2(i) = nchild(i)
         do j=1,4
           ichild2(j,i) = ichild(j,i)
         enddo
       enddo
C$OMP END PARALLEL DO       
       

       return
       end
c
c
c
c
c
c-------------------------------------------------------------      
      subroutine surf_tree_fix_lr(fun,nd,dpars,zpars,ipars,norder,
     1       npbox,fvals,grid,centers,nlevels,nboxes,boxsize,
     2       nbmax,nlmax,laddr,ilevel,iparent,nchild,ichild,nnbors,
     3       nbors)
c
c
c       convert an adaptive tree into a level restricted tree
c
      implicit none
      integer nd,ipars(*),norder,npbox,nlevels,nboxes,nlmax
      integer nbmax
      real *8 dpars(*),fvals(nd,npbox,nbmax),grid(2,npbox)
      real *8 centers(2,nbmax),boxsize(0:nlmax)
      complex *16 zpars(*)
      integer laddr(2,0:nlmax),ilevel(nbmax),iparent(nbmax)
      integer nchild(nbmax),ichild(4,nbmax),nnbors(nbmax)
      integer nbors(9,nbmax)
      integer laddrtail(2,0:nlmax),isum
      integer, allocatable :: iflag(:)

      integer i,j,k,l,ibox,jbox,kbox,ilev,idad,igranddad
      integer nbloc,ict,iper
      real *8 xdis,ydis,distest

      external fun

      allocate(iflag(nbmax))

c     Initialize flag array
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,nboxes
         iflag(i) = 0
      enddo
C$OMP END PARALLEL DO     



c     Flag boxes that violate level restriction by "1"
c     Violatioin refers to any box that is directly touching
c     a box that is more than one level finer
c
c     Method:
c     1) Carry out upward pass. For each box B, look at
c     the colleagues of B's grandparent
c     2) See if any of those colleagues are childless and in
c     contact with B.
c
c     Note that we only need to get up to level two, as
c     we will not find a violation at level 0 and level 1
c
c     For such boxes, we set iflag(i) = 1
c
      do ilev=nlevels,2,-1
c        This is the distance to test if two boxes separated
c        by two levels are touching
         distest = 1.05d0*(boxsize(ilev-1) + boxsize(ilev-2))/2.0d0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,idad,igranddad,i,jbox)         
C$OMP$ PRIVATE(ict,xdis,ydis)
         do ibox = laddr(1,ilev),laddr(2,ilev) 
            idad = iparent(ibox)
            igranddad = iparent(idad)
            
c           Loop over colleagues of granddad            
            do i=1,nnbors(igranddad)
               jbox = nbors(i,igranddad)
c              Check if the colleague of grandad
c              is a leaf node. This automatically
c              eliminates the granddad
               if(nchild(jbox).eq.0.and.iflag(jbox).eq.0) then
                   xdis = centers(1,jbox) - centers(1,idad)
                   ydis = centers(2,jbox) - centers(2,idad)
                   ict = 0
                   if(abs(xdis).le.distest) ict = ict + 1
                   if(abs(ydis).le.distest) ict = ict + 1
                   if(ict.eq.2) then
                      iflag(jbox) = 1
                   endif
               endif
c              End of checking criteria for the colleague of
c              granddad
            enddo
c           End of looping over colleagues of
c           granddad
         enddo
c        End of looping over boxes at ilev         
C$OMP END PARALLEL DO
      enddo
c     End of looping over levels and flagging boxes


c     Find all boxes that need to be given a flag+
c     A flag+ box will be denoted by setting iflag(box) = 2
c     This refers to any box that is not already flagged and
c     is bigger than and is contacting a flagged box
c     or another box that has already been given a flag +.
c     It is found by performing an upward pass and looking
c     at the flagged box's parents colleagues and a flag+
c     box's parents colleagues and seeing if they are
c     childless and present the case where a bigger box 
c     is contacting a flagged or flag+ box.

      do ilev = nlevels,1,-1
c        This is the distance to test if two boxes separated
c        by one level are touching
         distest = 1.05d0*(boxsize(ilev) + boxsize(ilev-1))/2.0d0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,idad,i,jbox,xdis,ydis)
C$OMP$PRIVATE(ict)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            if(iflag(ibox).eq.1.or.iflag(ibox).eq.2) then
               idad = iparent(ibox)
c              Loop over dad's colleagues               
               do i=1,nnbors(idad)
                  jbox = nbors(i,idad)
c                 Check if the colleague of dad
c                 is a leaf node. This automatically
c                 eliminates the dad
                  if(nchild(jbox).eq.0.and.iflag(jbox).eq.0) then
                     xdis = centers(1,jbox) - centers(1,ibox)
                     ydis = centers(2,jbox) - centers(2,ibox)
                     ict = 0
                     if(abs(xdis).le.distest) ict = ict + 1
                     if(abs(ydis).le.distest) ict = ict + 1
                     if(ict.eq.2) then
                        iflag(jbox) = 2
                     endif
                  endif
c                 End of checking criteria for the colleague of
c                dad
               enddo
c              End of looping over dad's colleagues               
            endif
c           End of checking if current box is relevant for
c           flagging flag+ boxes
         enddo
c        End of looping over boxes at ilev        
C$OMP END PARALLEL DO 
      enddo
c     End of looping over levels

c     Subdivide all flag and flag+ boxes. Flag all the children
c     of flagged boxes as flag++. Flag++ boxes are denoted
c     by setting iflag(box) = 3. The flag++ boxes need 
c     to be checked later to see which of them need further
c     refinement. While creating new boxes, we will
c     need to update all the tree structures as well.
c     Note that all the flagged boxes live between
c     levels 1 and nlevels - 2. We process the boxes via a
c     downward pass. We first determine the number of boxes
c     that are going to be subdivided at each level and 
c     everything else accordingly
      do ilev = 0,nlevels
         laddrtail(1,ilev) = 0
         laddrtail(2,ilev) = -1
      enddo
 
      do ilev = 1,nlevels-2
c        First subdivide all the flag and flag+
c        boxes with boxno nboxes+1, nboxes+ 2
c        and so on. In the second step, we reorganize
c        all the structures again to bring it back
c        in the standard format

         laddrtail(1,ilev+1) = nboxes+1

         nbloc = laddr(2,ilev)-laddr(1,ilev)+1

         call surf_tree_refine_boxes_flag(iflag,nd,npbox,fvals,
     1    fun,dpars,zpars,ipars,grid,nbmax,laddr(1,ilev),nbloc,
     2    centers,boxsize(ilev+1),nboxes,ilev,ilevel,iparent,nchild,
     3    ichild)

         laddrtail(2,ilev+1) = nboxes
      enddo
c     Reorganize the tree to get it back in the standard format

      call surf_tree_reorg(nboxes,nd,npbox,centers,nlevels,laddr,
     1          laddrtail,ilevel,iparent,nchild,ichild,fvals,iflag)

c       Compute colleague information again      

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
      do i=1,nboxes
         nnbors(i) = 0
         do j=1,9
            nbors(j,i) = -1
         enddo
      enddo
C$OMP END PARALLEL DO     
      iper = 0
      call surf_computecoll2d(nlevels,nboxes,laddr, boxsize,
     1                   centers,iparent,nchild,
     2                   ichild,iper,nnbors,nbors)

c     Processing of flag and flag+ boxes is done
c     Start processing flag++ boxes. We will use a similar
c     strategy as before. We keep checking the flag++
c     boxes that require subdivision if they still
c     violate the level restriction criterion, create
c     the new boxes, append them to the end of the list to begin
c     with and in the end reorganize the tree structure.
c     We shall accomplish this via a downward pass
c     as new boxes that get added in the downward pass
c     will also be processed simultaneously.
c     We shall additionally also need to keep on updating
c     the colleague information as we proceed in the 
c     downward pass

c     Reset the flags array to remove all the flag and flag+
c     cases. This is to ensure reusability of the subdivide
c     _flag routine to handle the flag++ case

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox)
      do ibox=1,nboxes
         if(iflag(ibox).ne.3) iflag(ibox) = 0
      enddo
C$OMP END PARALLEL DO      
 
      do ilev = 0,nlevels
         laddrtail(1,ilev) = 0
         laddrtail(2,ilev) = -1
      enddo

      do ilev = 2,nlevels-2
c     Step 1: Determine which of the flag++ boxes need
c     further division. In the even a flag++ box needs
c     further subdivision then flag the box with iflag(box) = 1
c     This will again ensure that the subdivide_flag routine
c     will take care of handling the flag++ case
         call surf_updateflags(ilev,nboxes,nlevels,laddr,nchild,ichild,
     1                    nnbors,nbors,centers,boxsize,iflag)

         call surf_updateflags(ilev,nboxes,nlevels,laddrtail,nchild,
     1          ichild,nnbors,nbors,centers,boxsize,iflag)
         
c      Step 2: Subdivide all the boxes that need subdivision
c      in the laddr set and the laddrtail set as well
         laddrtail(1,ilev+1) = nboxes + 1

         nbloc = laddr(2,ilev)-laddr(1,ilev)+1
         call surf_tree_refine_boxes_flag(iflag,nd,npbox,fvals,
     1    fun,dpars,zpars,ipars,grid,nbmax,laddr(1,ilev),nbloc,
     2    centers,boxsize(ilev+1),nboxes,ilev,ilevel,iparent,nchild,
     3    ichild)

         nbloc = laddrtail(2,ilev)-laddrtail(1,ilev)+1

         call surf_tree_refine_boxes_flag(iflag,nd,npbox,fvals,
     1    fun,dpars,zpars,ipars,grid,nbmax,laddrtail(1,ilev),nbloc,
     2    centers,boxsize(ilev+1),nboxes,ilev,ilevel,iparent,nchild,
     3    ichild)

          laddrtail(2,ilev+1) = nboxes         
c      Step 3: Update the colleague information for the newly
c      created boxes

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,i,idad,jbox,j,kbox)
          do ibox = laddrtail(1,ilev+1),laddrtail(2,ilev+1)
            nnbors(ibox) = 0
c           Find the parent of the current box         
            idad = iparent(ibox)
c           Loop over the neighbors of the parent box
c           to find out colleagues
            do i=1,nnbors(idad)
                jbox = nbors(i,idad)
                do j=1,4
c               ichild(j,jbox) is one of the children of the
c               neighbors of the parent of the current
c               box
                   kbox = ichild(j,jbox)
                   if(kbox.gt.0) then
c               Check if kbox is a nearest neighbor or in list 2
                      if((abs(centers(1,kbox)-centers(1,ibox)).le.
     1                   1.05*boxsize(ilev+1)).and.
     2                   (abs(centers(2,kbox)-centers(2,ibox)).le.
     3                   1.05*boxsize(ilev+1))) then
                     
                         nnbors(ibox) = nnbors(ibox)+1
                         nbors(nnbors(ibox),ibox) = kbox
                      endif
                   endif
                enddo
            enddo
c           End of computing colleagues of box i
         enddo
C$OMP END PARALLEL DO         
      enddo

c     Reorganize tree once again and we are all done      
      call surf_tree_reorg(nboxes,nd,npbox,centers,nlevels,laddr,
     1          laddrtail,ilevel,iparent,nchild,ichild,fvals,iflag)

c     Compute colleague information again      

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
      do i=1,nboxes
         nnbors(i) = 0
         do j=1,9
            nbors(j,i) = -1
         enddo
      enddo
C$OMP END PARALLEL DO    

      call surf_computecoll2d(nlevels,nboxes,laddr, boxsize,
     1                   centers,iparent,nchild,
     2                   ichild,iper,nnbors,nbors)
      

      return
      end
      

c-------------------------------------------------------------      
      subroutine surf_tree_reorg(nboxes,nd,npbox,centers,nlevels,laddr,
     1     laddrtail,ilevel,iparent,nchild,ichild,fvals,iflag)

c    This subroutine reorganizes the current data in all the tree
c    arrays to rearrange them in the standard format.
c    The boxes on input are assumed to be arranged in the following
c    format
c    boxes on level i are the boxes from laddr(1,i) to 
c    laddr(2,i) and also from laddrtail(1,i) to laddrtail(2,i)
c
c    At the end of the sorting, the boxes on level i
c    are arranged from laddr(1,i) to laddr(2,i)  
c
c    INPUT/OUTPUT arguments
c    nboxes         in: integer
c                   number of boxes
c
c    nd             in: integer
c                   number of real value functions
c
c    npbox          in: integer
c                   number of grid points per function
c
c    centers        in/out: double precision(3,nboxes)
c                   x and y coordinates of the center of boxes
c
c    nlevels        in: integer
c                   Number of levels in the tree
c
c    laddr          in/out: integer(2,0:nlevels)
c                   boxes at level i are numbered between
c                   laddr(1,i) to laddr(2,i)
c
c    laddrtail      in: integer(2,0:nlevels)
c                   new boxes to be added to the tree
c                   structure are numbered from
c                   laddrtail(1,i) to laddrtail(2,i)
c
c    ilevel      in/out: integer(nboxes)
c                ilevel(i) is the level of box i
c
c    iparent     in/out: integer(nboxes)
c                 iparent(i) is the parent of box i
c
c    nchild      in/out: integer(nboxes)
c                nchild(i) is the number of children 
c                of box i
c
c    ichild       in/out: integer(8,nboxes)
c                 ichild(j,i) is the jth child of box i
c
c    iflag        in/out: integer(nboxes)
c                 iflag(i) is a flag for box i required to generate
c                 level restricted tree from adaptive tree

      implicit none
c     Calling sequence variables and temporary variables
      integer nboxes,nlevels,npbox,nd
      double precision centers(2,nboxes)
      integer laddr(2,0:nlevels), tladdr(2,0:nlevels)
      integer laddrtail(2,0:nlevels)
      integer ilevel(nboxes)
      integer iparent(nboxes)
      integer nchild(nboxes)
      integer ichild(4,nboxes)
      integer iflag(nboxes)
      double precision fvals(nd,npbox,nboxes)
      
      integer, allocatable :: tilevel(:),tiparent(:),tnchild(:)
      integer, allocatable :: tichild(:,:),tiflag(:)
      integer, allocatable :: iboxtocurbox(:),ilevptr(:),ilevptr2(:)

      double precision, allocatable :: tfvals(:,:,:),tcenters(:,:)



c     Temporary variables
      integer i,j,k,l
      integer ibox,ilev, curbox,idim,nblev

      allocate(tilevel(nboxes),tiparent(nboxes),tnchild(nboxes))
      allocate(tichild(4,nboxes),tiflag(nboxes),iboxtocurbox(nboxes))
      allocate(tfvals(nd,npbox,nboxes),tcenters(2,nboxes))

      do ilev = 0,nlevels
         tladdr(1,ilev) = laddr(1,ilev)
         tladdr(2,ilev) = laddr(2,ilev)
      enddo
      call surf_tree_copy(nd,nboxes,npbox,centers,ilevel,iparent,nchild,
     1            ichild,fvals,tcenters,tilevel,tiparent,tnchild,
     2            tichild,tfvals)

      do ibox=1,nboxes
         tiflag(ibox) = iflag(ibox)
      enddo
     
c     Rearrange old arrays now

      do ilev = 0,1
         do ibox = laddr(1,ilev),laddr(2,ilev)
           iboxtocurbox(ibox) = ibox
         enddo
      enddo

      allocate(ilevptr(nlevels+1),ilevptr2(nlevels))

      ilevptr(2) = laddr(1,2)


      do ilev=2,nlevels
        nblev = laddr(2,ilev)-laddr(1,ilev)+1
        ilevptr2(ilev) = ilevptr(ilev) + nblev
        nblev = laddrtail(2,ilev)-laddrtail(1,ilev)+1
        ilevptr(ilev+1) = ilevptr2(ilev) + nblev
      enddo

      curbox = laddr(1,2)
      do ilev=2,nlevels
         laddr(1,ilev) = curbox
         do ibox = tladdr(1,ilev),tladdr(2,ilev)
            ilevel(curbox) = tilevel(ibox)
            nchild(curbox) = tnchild(ibox)
            centers(1,curbox) = tcenters(1,ibox)
            centers(2,curbox) = tcenters(2,ibox)
            do i=1,npbox
              do idim=1,nd
                fvals(idim,i,curbox) = tfvals(idim,i,ibox)
              enddo
            enddo
            iflag(curbox) = tiflag(ibox)
            iboxtocurbox(ibox) = curbox

            curbox = curbox + 1
         enddo
         do ibox = laddrtail(1,ilev),laddrtail(2,ilev)
            ilevel(curbox) = tilevel(ibox)
            centers(1,curbox) = tcenters(1,ibox)
            centers(2,curbox) = tcenters(2,ibox)
            nchild(curbox) = tnchild(ibox)
            do i=1,npbox
              do idim=1,nd
                fvals(idim,i,curbox) = tfvals(idim,i,ibox)
              enddo
            enddo
            iflag(curbox) = tiflag(ibox)
            iboxtocurbox(ibox) = curbox

            curbox = curbox + 1
         enddo
         laddr(2,ilev) = curbox-1
      enddo

c     Handle the parent children part of the tree 
c     using the mapping iboxtocurbox

      do ibox=1,nboxes
         if(tiparent(ibox).eq.-1) iparent(iboxtocurbox(ibox)) = -1
         if(tiparent(ibox).gt.0) 
     1    iparent(iboxtocurbox(ibox)) = iboxtocurbox(tiparent(ibox))
         do i=1,4
            if(tichild(i,ibox).eq.-1) ichild(i,iboxtocurbox(ibox)) = -1
            if(tichild(i,ibox).gt.0) 
     1      ichild(i,iboxtocurbox(ibox)) = iboxtocurbox(tichild(i,ibox))
         enddo
      enddo

      return
      end
c--------------------------------------------------------------------      
      subroutine surf_updateflags(curlev,nboxes,nlevels,laddr,nchild,
     1      ichild,nnbors,nbors,centers,boxsize,iflag)

c      This subroutine is to check the boxes flagged as flag++
c      and determine which of the boxes need refinement. The flag
c      of the box which need refinement is updated to iflag(box)=1
c      and that of the boxes which do not need refinement is
c      updated to iflag(box) = 0
c
c      INPUT arguments
c      curlev         in: integer
c                     the level for which boxes need to be processed
c
c      nboxes         in: integer
c                     total number of boxes
c
c      nlevels        in: integer
c                     total number of levels
c
c      laddr          in: integer(2,0:nlevels)
c                     boxes from laddr(1,ilev) to laddr(2,ilev)
c                     are at level ilev
c
c      nchild         in: integer(nboxes)
c                     nchild(ibox) is the number of children
c                     of box ibox
c
c      ichild         in: integer(4,nboxes)
c                     ichild(j,ibox) is the box id of the jth
c                     child of box ibox
c
c      nnbors         in: integer(nboxes)
c                     nnbors(ibox) is the number of colleagues
c                     of box ibox
c
c      nbors          in: integer(9,nboxes)
c                     nbors(j,ibox) is the jth colleague of box
c                     ibox
c
c      centers        in: double precision(3,nboxes)
c                     x and y coordinates of the box centers
c
c      boxsize        in: double precision(0:nlevels)
c                     boxsize(i) is the size of the box at level i
c
c      iflag          in/out: integer(nboxes)
c                     iflag(ibox)=3 if it is flag++. iflag(ibox) =1
c                     or 0 at the end of routine depending on
c                     whether box needs to be subdivided or not
c
      implicit none
c     Calling sequence variables
      integer curlev, nboxes, nlevels
      integer laddr(2,0:nlevels),nchild(nboxes),ichild(4,nboxes)
      integer nnbors(nboxes), nbors(9,nboxes)
      integer iflag(nboxes)
      double precision centers(2,nboxes),boxsize(0:nlevels)

c     Temporary variables
      integer i,j,k,l,ibox,jbox,kbox,lbox, ict
      double precision distest,xdis,ydis

      distest = 1.05d0*(boxsize(curlev) + boxsize(curlev+1))/2.0d0
c     Loop over all boxes at the current level     

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,i,jbox,j,kbox,xdis,ydis)
C$OMP$PRIVATE(ict)
      do ibox = laddr(1,curlev),laddr(2,curlev)
         if(iflag(ibox).eq.3) then
            iflag(ibox) = 0
c           Loop over colleagues of the current box      
            do i=1,nnbors(ibox)
c              Loop over colleagues of flag++ box        
               jbox = nbors(i,ibox)
              
c              Loop over the children of the colleague box
c              Note we do not need to exclude self from
c              the list of colleagues as a self box which
c              is flag++ does not have any children 
c              and will not enter the next loop
               do j=1,4
                  kbox = ichild(j,jbox)
                  if(kbox.gt.0) then
                     if(nchild(kbox).gt.0) then
                        xdis = centers(1,kbox) - centers(1,ibox)
                        ydis = centers(2,kbox) - centers(2,ibox)
                        ict = 0
                        if(abs(xdis).le.distest) ict = ict + 1
                        if(abs(ydis).le.distest) ict = ict + 1
                        if(ict.eq.2) then
                           iflag(ibox) = 1
                           goto 1111
                        endif
                     endif
                  endif
c                 End of looping over the children of the child
c                 of the colleague box
               enddo
c              End of looping over the children of the colleague box       
            enddo
c           End of looping over colleagues            
 1111       continue        
         endif
c        End of testing if the current box needs to checked for         
      enddo
c     End of looping over boxes at the current level      
C$OMP END PARALLEL DO      

      return
      end
c
c
c
c
c


      subroutine surf_tree_refine_boxes_flag(iflag,nd,npbox,fvals,
     1  fun,dpars,zpars,ipars,grid,nboxes,ifirstbox,nbloc,centers,
     2  bs,nbctr,nlctr,ilevel,iparent,nchild,ichild)
      implicit none
      integer nd,npbox,nboxes
      real *8 fvals(nd,npbox,nboxes)
      integer nbloc,nbctr,nlctr
      real *8 centers(2,nboxes),bs,grid(2,npbox),xy(2)
      integer ilevel(nboxes),iparent(nboxes)
      integer ichild(4,nboxes),nchild(nboxes)
      integer iflag(nboxes)
      integer ifirstbox,ilastbox
      integer, allocatable :: isum(:)
      real *8 dpars(*)
      integer ipars(*)
      complex *16 zpars(*)

      integer i,ibox,nel,j,l,jbox,nbl,ii
      integer xind(4),yind(4)

      real *8 bsh
      data xind/-1,1,-1,1/
      data yind/-1,-1,1,1/

      external fun


      ilastbox = ifirstbox+nbloc-1
      bsh = bs/2

      allocate(isum(nbloc))

      call surf_cumsum_nz(nbloc,iflag(ifirstbox),isum)
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,j,jbox,nbl,l)
C$OMP$PRIVATE(xy)
      do ibox = ifirstbox,ilastbox
        if(iflag(ibox).gt.0) then
          nchild(ibox) = 4
          nbl = nbctr + (isum(ibox-ifirstbox+1)-1)*4
          do j=1,4
            jbox = nbl+j
            centers(1,jbox) = centers(1,ibox)+xind(j)*bsh
            centers(2,jbox) = centers(2,ibox)+yind(j)*bsh
            do l=1,npbox
              xy(1) = centers(1,jbox) + grid(1,l)*bs
              xy(2) = centers(2,jbox) + grid(2,l)*bs
              call fun(nd,xy,dpars,zpars,ipars,fvals(1,l,jbox))
            enddo
            iparent(jbox) = ibox
            nchild(jbox) = 0
            do l=1,4
              ichild(l,jbox) = -1
            enddo
            ichild(j,ibox) = jbox
            ilevel(jbox) = nlctr+1 
            if(iflag(ibox).eq.1) iflag(jbox) = 3
            if(iflag(ibox).eq.2) iflag(jbox) = 0
          enddo
        endif
      enddo
C$OMP END PARALLEL DO
      
      if(nbloc.gt.0) nbctr = nbctr + isum(nbloc)*4


      return
      end
c
c
c     
c
      subroutine surf_cumsum_nz(n,a,b)
c
c        this subroutine computes the cumulative sum of postive
c        entries of an array
c
c
c       TODO: need to make this routine openmp
c
c       input:
c         n - number of elements
c         a - integer(n)
c              input array
c
c       output:
c         b - integer(n)
c            b(i) = sum_{j=1}^{i} I_{a(j)>0}
c
      implicit none
      integer a(n),b(n),n,i,isum

      isum = 0
      do i=1,n
        if(a(i).gt.0) isum = isum+1
        b(i) = isum
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
c
c
c

      subroutine surf_computecoll2d(nlevels,nboxes,laddr,boxsize,
     1                       centers,iparent,nchild,ichild,iper,
     2                       nnbors,nbors)

c     This subroutine computes the colleagues for an adaptive
c     pruned tree. box j is a colleague of box i, if they share a
c     vertex or an edge and the two boxes are at the same
c     level in the tree
c
c     INPUT arguments
c     nlevels     in: integer
c                 Number of levels
c
c     nboxes      in: integer
c                 Total number of boxes
c
c     laddr       in: integer(2,0:nlevels)
c                 indexing array providing access to boxes at
c                 each level. 
c                 the first box on level i is laddr(1,i)
c                 the last box on level i is laddr(2,i)
c
c     boxsize     in: double precision(0:nlevels)
c                 Array of boxsizes
c 
c     centers     in: double precision(2,nboxes)
c                 array of centers of boxes
c   
c     iparent     in: integer(nboxes)
c                 iparent(i) is the box number of the parent of
c                 box i
c
c     nchild      in: integer(nboxes)
c                 nchild(i) is the number of children of box i
c
c     ichild      in: integer(4,nboxes)
c                 ichild(j,i) is the box id of the jth child of
c                 box i
c
c     iper        in: integer
c                 flag for periodic implementations. 
c                 Currently not used. Feature under construction.
c
c----------------------------------------------------------------
c     OUTPUT
c     nnbors      out: integer(nboxes)
c                 nnbors(i) is the number of colleague boxes of
c                 box i
c
c     nbors       out: integer(9,nboxes)
c                 nbors(j,i) is the box id of the jth colleague
c                 box of box i
c---------------------------------------------------------------
      implicit none
      integer nlevels,nboxes
      integer iper
      integer laddr(2,0:nlevels)
      double precision boxsize(2,0:nlevels)
      double precision centers(2,nboxes)
      integer iparent(nboxes), nchild(nboxes), ichild(4,nboxes)
      integer nnbors(nboxes)
      integer nbors(9,nboxes)

c     Temp variables
      integer ilev,ibox,jbox,kbox,dad
      integer i,j,ifirstbox,ilastbox


c     Setting parameters for level = 0
      nnbors(1) = 1
      nbors(1,1) = 1
      do ilev = 1,nlevels
c        Find the first and the last box at level ilev      
         ifirstbox = laddr(1,ilev)
         ilastbox = laddr(2,ilev)
c        Loop over all boxes to evaluate neighbors, list1 and updating
c        hunglists of targets

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,dad,i,jbox,j,kbox)
         do ibox = ifirstbox,ilastbox
c           Find the parent of the current box         
            dad = iparent(ibox)
c           Loop over the neighbors of the parent box
c           to find out list 1 and list 2
            do i=1,nnbors(dad)
                jbox = nbors(i,dad)
                do j=1,4
c               ichild(j,jbox) is one of the children of the
c               neighbors of the parent of the current
c               box
                   kbox = ichild(j,jbox)
                   if(kbox.gt.0) then
c               Check if kbox is a nearest neighbor or in list 2
                      if((abs(centers(1,kbox)-centers(1,ibox)).le.
     1                   1.05*boxsize(1,ilev)).and.
     2                   (abs(centers(2,kbox)-centers(2,ibox)).le.
     3                   1.05*boxsize(2,ilev))) then
                     
                         nnbors(ibox) = nnbors(ibox)+1
                         nbors(nnbors(ibox),ibox) = kbox
                      endif
                   endif
                enddo
            enddo
c           End of computing colleagues of box i
         enddo
C$OMP END PARALLEL DO         
      enddo

      return
      end
c
c
c
c
      subroutine surf_mesh2d(x,nx,y,ny,xy)
      implicit real *8 (a-h,o-z)
      dimension x(*), y(*), xy(2,*)

      ind = 0
      do iy = 1,ny
         do ix = 1,nx
            ind = ind+1
            xy(1,ind) = x(ix)
            xy(2,ind) = y(iy)
         enddo
      enddo

      return
      end
c
c
c
c

c      
c     end of file
c
c
c
c
      
      
