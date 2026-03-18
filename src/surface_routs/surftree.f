c     
c     
c     generate oct tree based on resolving a 
c     function to desired precision. Observe that the tree is not
c     level restricted
c     
c     
c     The function handle is of the form
c     call fun(nd,xy,dpars,zpars,ipars,f)
c     
c     where xy is the location in (-L/2,L/2)^2
c     
c     dpars are a collection of real parameters, zpars are complex
c     parameters and ipars are integer *8 parameters, and f is a 
c     real array of size nd
c     
c     
c     For the Helmholtz/Maxwell tree, the boxes are refined until 
c     Re(zk)*boxsize<5, beyond which the function resolution criterion
c     kicks in. 
c     
c     A function is said to be resolved if its interpolant at the 8
c     child nodes, agrees with the function values at those nodes
c     up to the user specified tolerance. 
c     The error is scaled by h**(eta)
c     where eta is user specified and h is the boxsize. If there is 
c     any confusion, the user should seet \eta to 0
c     Let \tilde{f} denote the interpolant of f, then
c     the refinement criterion is 
c     \int_{B_{j} |\tilde{f}-f|^{p} *h^{\eta} < 
c     \varepsilon V_{j}^{1/p}/V_{0}^{1/p}/(\int_{B_{0}}|f|^{p})^{1/p}
c     This implies that
c     \int_{B_{0}} |\tilde{f}-f|^{p} =
c     \sum_{j} \int_{B_{j}} |\tilde{f}-f|^{p}
c     \leq \sum_{j} \eps^{p}*h^{\eta p} 
c     V_{j}/V_{0}/(\int_{B_{0}} |f|^p)
c     If \eta = 0,
c     \leq \eps^{p}/(\int_{B_{0}} |f|^{p})
c     
c     i.e., this strategy guarantees that the interpolated function
c     approximates the function with relative lp accuracy of \eps
c     
c     This code has 2 main user callable routines
c     make_vol_tree_mem -> Returns the memory requirements, 
c     tree length, number of boxes, number of levels
c     make_vol_tree -> Makes the actual tree, returns centers of boxes,
c     colleague info, function values on leaf boxes
c     
c     
c     iptr(1) - laddr
c     iptr(2) - ilevel
c     iptr(3) - iparent
c     iptr(4) - nchild
c     iptr(5) - ichild
c     iptr(6) - ncoll
c     iptr(7) - coll
c     iptr(8) - ltree
c     




      subroutine surf_tree_mem(nuv,eps,zk,boxlen,norder,iptype,
     1     eta,fun_dummyint,nd,dpars,zpars,ipars,nlevels,nboxes,ltree,
     2     rintlin,nboxesleaf)
c
c      get memory requirements for the tree
c
c     input parameters:
c        nuv - integer *8(2)
c           nuv(1): split [-boxlen/2,boxlen/2]^2 into nuv(1) subdivisions at root level in the x-direction,
c           nuv(1) = 1 means no subdivisions in the x-direction
c           nuv(2): split [-boxlen/2,boxlen/2]^2 into nuv(2) subdivisions at root level in the y-direction
c           nuv(2) = 1 means no subdivisions in the y-direction
c        eps - double precision
c           precision requested
c        zk - double complex
c           Helmholtz parameter
c        boxlen - double precision
c           length of box in which volume is contained, 
c           if boxlen = L, then volume is [-L/2,L/2]^2
c        norder - integer *8
c           order of discretization
c        eta - double precision
c           scaling parameter for error
c        fun_dummyint - DUMMY INT .function handle
c           function to evalute it everywhere in the volume
c        nd - integer *8
c           number of real functions returned by fun
c        dpars - double precision
c           real parameters for function evaluation
c        zpars - double complex
c           complex parameters for function evaluation
c        ipars - integer *8
c           integer *8 parameters for function evaluation. ipars(1) is
c           reserved for function selection in fun
c 
c        output:
c           nlevels - integer *8
c             number of levels
c           nboxes - integer *8
c             number of boxes
c           ltree - integer *8
c             length of tree
c           rintlin(1:nlevels+1) - real *8
c             lp norm to scale the functions by
c             (on input rintl should be of size(1:201)
c           nboxesleaf - integer *8
c             number of leaf boxes      
c
c      
c

      implicit none

      integer *8 i,itype,j,npols,k,nbloc,nbctr,nbadd,irefine,idim,iper
      integer *8 nbmax,nlmax,npbox,npc,ilev,ifirstbox,ilastbox,iii
      integer *8 nd,ipars(*),iptype,nuv(2),nroot,fun_dummyint,nbtot
      integer *8 nlevels,nboxes,ltree,norder,nboxesleaf

      real *8 eps,boxlen,eta,dpars(*),delta,xy(2),rint
      real *8 rintl(0:200),rintlin(201),rsc,ra

      complex *16 zpars(*),zk
      
      external fun

      integer *8, allocatable :: laddr(:,:),ilevel(:),iparent(:)
      integer *8, allocatable :: ichild(:,:),ncoll(:),icoll(:,:)
      integer *8, allocatable :: nbors(:,:),nnbors(:),nchild(:)
      integer *8, allocatable :: ilevel2(:),iparent2(:),nchild2(:)
      integer *8, allocatable :: irefinebox(:),ichild2(:,:)
      
      real *8, allocatable :: centers(:,:)
      real *8, allocatable :: centers2(:,:),fvals2(:,:,:)
      real *8, allocatable :: grid(:,:),qwts(:)
      real *8, allocatable :: xq(:),wts(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: rintbs(:),rintbs2(:)
      real *8, allocatable :: fvals(:,:,:)
      real *8, allocatable :: wts2(:),xref2(:,:),umat2(:,:),vmat2(:,:)
      real *8, allocatable :: boxsize(:,:)
      
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
      nroot = nuv(1)*nuv(2)
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
         boxsize(1,0) = boxlen/nuv(1)
         boxsize(2,0) = boxlen/nuv(2)
      endif

      k = 1
      do i=1,nuv(1)    
         do j=1,nuv(2)
            centers(1,k)=-boxlen/2.0d0
     1        + boxlen/(nuv(1))/2.0d0
     2        + boxlen/(nuv(1))*(i-1)     
            centers(2,k)=-boxlen/2.0d0
     1           + boxlen/(nuv(2))/2.0d0
     2           + boxlen/(nuv(2))*(j-1)           
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
c     call legeexps(itype,norder,xq,umat,vmat,wts)
      call chebexps(itype,norder,xq,umat,vmat,wts)
      do i=1,norder
        xq(i) = xq(i)/2
      enddo

      call mesh2d(xq,norder,xq,norder,grid)

      k = 0
      do j = 1,norder
         do i = 1,norder
            k = k+1
            grid(1,k) = xq(i)
            grid(2,k) = xq(j)
         enddo
      enddo

      npols = norder**2
      allocate(wts2(npbox),xref2(2,npbox))
      allocate(umat2(npols,npols),vmat2(npols,npols))
      itype = 1
c      call legetens_exps_2d(itype,norder,'f',xref2,umat2,npols,
c     1     vmat,1,wts2)
      call chebtens_exps_2d(itype,norder,'f',xref2,umat2,npols,
     1     vmat2,1,wts2)
      
c
c       compute fvals at the grid
c

      rint = 0

c
c   note extra factor of 4 since wts2 are on [-1,1]^2 
c   as opposed to [-1/2,1/2]^2
c
c      rsc = boxlen**2/4
      rsc = boxsize(1,0)*boxsize(2,0)/(4.0)
      do k=1,nroot

         do i=1,npbox            
            xy(1) = centers(1,k) + grid(1,i)*boxsize(1,0)
            xy(2) = centers(2,k) + grid(2,i)*boxsize(2,0)
            call fun(nd,xy,dpars,zpars,ipars,fvals(1,i,k))
            if(iptype.eq.0) then
               do idim=1,nd
                  if(abs(fvals(idim,i,k)).gt.rintbs(1)) rintbs(1) = 
     1                 abs(fvals(idim,i,k))
               enddo
            endif

            if(iptype.eq.1) then
               do idim=1,nd
                  rintbs(1) = rintbs(1)+abs(fvals(idim,i,k))*wts2(i)*rsc
               enddo
            endif

            if(iptype.eq.2) then
               do idim=1,nd
                  rintbs(1) = rintbs(1) + fvals(idim,i,k)**2*wts2(i)*rsc
               enddo
            endif
         enddo
      enddo

      if(iptype.eq.0.or.iptype.eq.1) rint = rintbs(1)
      if(iptype.eq.2) rint = sqrt(rintbs(1))

      rintl(0) = rint
      rintlin(1) = rint

      nbctr = nroot
     
      
      do ilev=0,nlmax-1
        irefine = 0

        ifirstbox = laddr(1,ilev) 
        ilastbox = laddr(2,ilev)

        nbloc = ilastbox-ifirstbox+1

        allocate(irefinebox(nbloc))

        
        if(iptype.eq.2) rsc = sqrt(1.0d0/boxlen**2)
        if(iptype.eq.1) rsc = (1.0d0/boxlen**2)
        if(iptype.eq.0) rsc = 1.0d0

        rsc = rsc*rint

        print *, ilev, rint, rsc
        
        call surf_tree_find_box_refine(nd,iptype,eta,eps,zk,norder,
     1       npbox,fvals,npols,umat,boxsize(1,ilev),nbmax,ifirstbox,
     2       nbloc,rsc,irefinebox,irefine)

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

          call update_rints(nd,npbox,nbmax,fvals,ifirstbox,nbloc,
     1       iptype,nchild,ichild,wts2,rsc,rintbs,rint)

          rintl(ilev+1) = rint
          
          rintlin(ilev+2) = rint
          
          laddr(2,ilev+1) = nbctr
        else
          exit
        endif
        deallocate(irefinebox)
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

        call computecoll2d(nlevels,nboxes,laddr,boxsize,centers,
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

      subroutine surf_tree_build(nuv,eps,zk,boxlen,norder,iptype,eta,
     1  fun_dummyint,nd,dpars,zpars,ipars,nlevels,nboxes,ltree,rintlin,
     2  itree,iptr,fvals,centers,boxsizein)
c
c      construct the tree
c
c     input parameters:
c        nuv - integer *8(2)
c           nuv(1): split [-boxlen/2,boxlen/2]^2 into nuv(1) subdivisions at root level in the x-direction,
c           nuv(1) = 1 means no subdivisions in the x-direction
c           nuv(2): split [-boxlen/2,boxlen/2]^2 into nuv(2) subdivisions at root level in the y-direction
c           nuv(2) = 1 means no subdivisions in the y-direction
c        eps - double precision
c           precision requested
c        zk - double complex
c           Helmholtz parameter
c        boxlen - double precision
c           length of box in which volume is contained, 
c           if boxlen = L, then volume is [-L/2,L/2]^2
c        norder - integer *8
c           order of discretization
c        iptype - integer *8
c           error norm
c           iptype = 0 - linf
c           iptype = 1 - l1
c           iptype = 2 - l2
c        eta - double precision
c           scaling parameter for error
c        fun - function handle
c           function to evalute it everywhere in the volume
c        nd - integer *8
c           number of real functions returned by fun
c        dpars - double precision
c           real parameters for function evaluation
c        zpars - double complex
c           complex parameters for function evaluation
c        ipars - integer *8
c           integer *8 parameters for function evaluation
c        nlevels - integer *8
c          number of levels
c        nboxes - integer *8
c          number of boxes
c        ltree - integer *8
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
c        itree - integer *8 (ltree)
c          tree info
c        iptr - integer *8(8)
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
c          xyz coordinates of box centers in the oct tree
c        boxsize - double precision (0:nlevels)
c          size of box at each of the levels
c

      implicit none
      real *8 eps,boxlen,eta,dpars(*)
      complex *16 zk,zpars(*)
      integer *8 nd,ipars(*),iptype,fun_dummyint
      integer *8 nlevels,nboxes,ltree,norder,nroot
      integer *8 iptr(8),nuv(2)
      integer *8 itree(ltree),ier
      real *8 fvals(nd,norder**2,nboxes),centers(2,nboxes)
      integer *8, allocatable :: irefinebox(:)
      real *8 boxsize(2,0:nlevels)
      real *8 boxsizein(2,nlevels+1)
      real *8, allocatable :: xq(:),wts(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: grid(:,:),ximat(:,:),qwts(:)
      real *8 rintl(0:nlevels), rintlin(nlevels+1)
      real *8 xy(2)

      integer *8 i,ilev,irefine,itype,nbmax,nlmax,npbox,npc,ii,k
      integer *8 ifirstbox,ilastbox,nbctr,nbloc
      real *8 rsc

      real *8 ra
      integer *8 j,nboxes0,nlevels0,npols,iper

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

      nroot = nuv(1)*nuv(2)
      itree(1) = 1
      itree(2) = nroot

      nbmax = 100000
      
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
         boxsize(1,0) = boxlen/nuv(1)
         boxsize(2,0) = boxlen/nuv(2)
      endif


      k = 1
      do i=1,nuv(1)    
         do j=1,nuv(2)
            centers(1,k)=-boxlen/2.0d0
     1           + boxlen/(nuv(1))/2.0d0
     2           + boxlen/(nuv(1))*(i-1)     
            centers(2,k)=-boxlen/2.0d0
     1           + boxlen/(nuv(2))/2.0d0
     2           + boxlen/(nuv(2))*(j-1)
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
c     call legeexps(itype,norder,xq,umat,vmat,wts)
      call chebexps(itype,norder,xq,umat,vmat,wts)
      do i=1,norder
        xq(i) = xq(i)/2
      enddo

c      call mesh2d(xq,norder,xq,norder,grid)
      k = 0
      do j = 1,norder
         do i = 1,norder
            k = k+1
            grid(1,k) = xq(i)
            grid(2,k) = xq(j)
         enddo
      enddo

      npols = norder**2
      
      do k=1,nroot
         do i=1,npbox           
            xy(1) = centers(1,k) + grid(1,i)*boxsize(1,0)
            xy(2) = centers(2,k) + grid(2,i)*boxsize(2,0)
            call fun(nd,xy,dpars,zpars,ipars,fvals(1,i,k))
         enddo
      enddo

      do i=0,nlevels
         rintl(i) = rintlin(i+1)
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
         
         allocate(irefinebox(nbloc))

         if(iptype.eq.2) rsc = sqrt(1.0d0/boxlen**2)
         if(iptype.eq.1) rsc = (1.0d0/boxlen**2)
         if(iptype.eq.0) rsc = 1.0d0
         rsc = rsc*rintl(ilev)
ccc custom         

         nbmax = nboxes
         
         call surf_tree_find_box_refine(nd,iptype,eta,eps,zk,norder,
     1        npbox,fvals,npols,umat,boxsize(1,ilev),nbmax,ifirstbox,
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
        deallocate(irefinebox)
      enddo

      

      do i=1,nboxes
         itree(iptr(6)+i-1) = 0
         do j=1,9
            itree(iptr(7)+9*(i-1)+j-1) = -1
         enddo
      enddo

      iper = 0

      call computecoll2d(nlevels,nboxes,itree(iptr(1)),boxsize,
     1        centers,itree(iptr(3)),itree(iptr(4)),itree(iptr(5)),iper,
     2     itree(iptr(6)),itree(iptr(7)))

c     populate output array to matlab
      do i=0,nlevels
         boxsizein(1,i+1) = boxsize(1,i)
         boxsizein(2,i+1) = boxsize(2,i)
      enddo
      
      if(nlevels.ge.2) then
c         call vol_tree_fix_lr(fun,nd,dpars,zpars,ipars,norder,npbox,
c     1       fvals,grid,centers,nlevels,nboxes0,boxsize,nboxes,nlevels,
c     2       itree(iptr(1)),itree(iptr(2)),itree(iptr(3)),
c     3       itree(iptr(4)),itree(iptr(5)),itree(iptr(6)),
c     4       itree(iptr(7)))
      endif
      
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
      integer *8 nd,npc,npbox,norder,iptype,npols
      integer *8 nboxes,nbloc
      real *8 eta,eps,fvals(nd,npbox,nboxes)
      real *8 umat(norder,norder)
      real *8 rsc
      real *8, allocatable :: fcoefs(:,:,:),rmask(:)
      integer *8, allocatable :: iind2p(:,:)
      real *8 alpha,beta,boxsize(2),rsum
      integer *8 irefinebox(nbloc),xind(4),yind(4)
      complex *16 zk
      integer *8 ifirstbox

      integer *8 irefine

      integer *8 i,j,k,l,ibox,ifunif,i1
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
      
      call legetens_ind2pow_2d(norder-1,'f',iind2p)
      
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

      bs = boxsize(1)*boxsize(2)/4.0d0
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
     1       fcoefs(1,1,i),umat)

        call fun_err(nd,npols,fcoefs(1,1,i),rmask,
     1       iptype,rscale2,err)
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
      integer *8 nd,npbox
      integer *8 ipars(*)
      real *8 dpars(*)
      complex *16 zpars(*)
      real *8 fvals(nd,npbox,nboxes)
      integer *8 nboxes,nbloc,nbctr,nlctr
      real *8 grid(2,npbox),bs(2),xy(2)
      real *8 centers(2,nboxes)
      integer *8 ilevel(nboxes),iparent(nboxes)
      integer *8 ichild(4,nboxes),nchild(nboxes)
      integer *8 irefinebox(nbloc)
      integer *8 ifirstbox
      integer *8, allocatable :: isum(:)

      integer *8 i,ibox,nel0,j,l,jbox,nel1,nbl
      integer *8 xind(4),yind(4)

      real *8 bsh(2)
      data xind/-1,1,-1,1/
      data yind/-1,-1,1,1/

      external fun

  
      allocate(isum(nbloc))
      call cumsum_surf(nbloc,irefinebox,isum)

c      bsh = bs/2
      bsh(1) = bs(1)/2.0d0
      bsh(2) = bs(2)/2.0d0

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ibox,nbl,j,jbox,l,xy)
      do i = 1,nbloc
        ibox = ifirstbox + i-1
        if(irefinebox(i).eq.1) then
          nbl = nbctr + (isum(i)-1)*4
          nchild(ibox) = 4
          do j=1,4
             jbox = nbl+j             
             centers(1,jbox) = centers(1,ibox)+xind(j)*bsh(1)
             centers(2,jbox) = centers(2,ibox)+yind(j)*bsh(2)
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
      subroutine fun_err(nd,n,fcoefs,rmask,iptype,rscale,err)
c       this subroutine estimates the error based on the expansion
c       coefficients in a given basis
c       
c       input
c        nd - integer *8
c          number of functions
c        n -  integer *8
c           number of points at which function is tabulated
c        fcoefs: double precision(nd,n) 
c          tensor product legendre coeffs of func
c        rmask: double precision(n)
c           coefficients which will be accounted for in the error
c        iptype: integer *8
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
        integer *8 n,i,iptype,idim,nd,int8_1
        real *8 rscale,err
        real *8 fcoefs(nd,n),rmask(n)
        real *8, allocatable :: errtmp(:),ftmp(:,:)
        real *8 alpha, beta

        allocate(errtmp(nd),ftmp(nd,n))

        alpha = 1.0d0
        beta = 0.0d0

        int8_1 = 1
   
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
          call dgemv_guru('n',nd,n,alpha,ftmp,nd,rmask,int8_1,beta,
     1         errtmp,int8_1)
        endif
        if(iptype.eq.2) then
          do i=1,n 
            do idim=1,nd
              ftmp(idim,i) = fcoefs(idim,i)**2
            enddo
         enddo
         call dgemv_guru('n',nd,n,alpha,ftmp,nd,rmask,int8_1,beta,
     1        errtmp,int8_1)
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
      subroutine update_rints(nd,npbox,nbmax,fvals,ifirstbox,nbloc,
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
c    - nd: integer *8
c        number of functions
c    - npbox: integer *8
c        number of points per box where the function is tabulated
c    - nbmax: integer *8
c        max number of boxes
c    - fvals: real *8 (nd,npbox,nbmax)
c        tabulated function values
c    - ifirstbox: integer *8
c        first box in the list of boxes to be processed
c    - nbloc: integer *8
c        number of boxes to be processed
c    - iptype: integer *8
c        Lp version of the scheme
c        * iptype = 0, linf
c        * iptype = 1, l1
c        * iptype = 2, l2
c    - nchild: integer *8(nbmax)
c        number of children 
c    - ichild: integer *8(8,nbmax)
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
      implicit integer *8 (i-n)
      integer *8, intent(in) :: nd,npbox,nbmax
      real *8, intent(in) :: fvals(nd,npbox,nbmax)
      integer *8, intent(in) :: ifirstbox,nbloc,iptype
      integer *8, intent(in) :: nchild(nbmax),ichild(4,nbmax)
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
      subroutine get_children_qwts(norder,npc,wts,qwts)
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      real *8 wts(norder),qwts(npc)
      
      ipt = 1
      do i=1,norder
        do j=1,norder
           qwts(ipt) = wts(i)*wts(j)/4
           ipt = ipt+1
        enddo
      enddo

      npbox = norder**2
      do i=1,3
        do j=1,npbox
          qwts(i*npbox+j) = qwts(j)
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
       subroutine surf_tree_copy(nd,nb,npb,centers,ilevel,iparent,
     1            nchild,ichild,fvals,centers2,ilevel2,iparent2,nchild2,
     2            ichild2,fvals2)

       implicit none
       integer *8 nd,nb,npb
       real *8 centers(2,nb),centers2(2,nb)
       integer *8 ilevel(nb),ilevel2(nb)
       integer *8 iparent(nb),iparent2(nb)
       integer *8 nchild(nb),nchild2(nb)
       integer *8 ichild(4,nb),ichild2(4,nb)
       real *8 fvals(nd,npb,nb),fvals2(nd,npb,nb)

       integer *8 i,j,nel

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
      integer *8 nd,ipars(*),norder,npbox,nlevels,nboxes,nlmax
      integer *8 nbmax
      real *8 dpars(*),fvals(nd,npbox,nbmax),grid(2,npbox)
      real *8 centers(2,nbmax),boxsize(2,0:nlmax)
      complex *16 zpars
      integer *8 laddr(2,0:nlmax),ilevel(nbmax),iparent(nbmax)
      integer *8 nchild(nbmax),ichild(4,nbmax),nnbors(nbmax)
      integer *8 nbors(9,nbmax)
      integer *8 laddrtail(2,0:nlmax),isum
      integer *8, allocatable :: iflag(:)

      integer *8 i,j,k,l,ibox,jbox,kbox,ilev,idad,igranddad
      integer *8 nbloc,ict,iper
      real *8 xdis,ydis,xdistest,ydistest

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
         xdistest = 1.05d0*(boxsize(1,ilev-1) + boxsize(1,ilev-2))/2.0d0
         ydistest = 1.05d0*(boxsize(2,ilev-1) + boxsize(2,ilev-2))/2.0d0
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
                   if(abs(xdis).le.xdistest) ict = ict + 1
                   if(abs(ydis).le.ydistest) ict = ict + 1
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
         xdistest = 1.05d0*(boxsize(1,ilev) + boxsize(1,ilev-1))/2.0d0
         ydistest = 1.05d0*(boxsize(2,ilev) + boxsize(2,ilev-1))/2.0d0
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
                     if(abs(xdis).le.xdistest) ict = ict + 1
                     if(abs(ydis).le.ydistest) ict = ict + 1
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
     2    centers,boxsize(1,ilev+1),nboxes,ilev,ilevel,iparent,nchild,
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
      call computecoll2d(nlevels,nboxes,laddr, boxsize,
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
     2    centers,boxsize(1,ilev+1),nboxes,ilev,ilevel,iparent,nchild,
     3    ichild)

         nbloc = laddrtail(2,ilev)-laddrtail(1,ilev)+1

         call surf_tree_refine_boxes_flag(iflag,nd,npbox,fvals,
     1    fun,dpars,zpars,ipars,grid,nbmax,laddrtail(1,ilev),nbloc,
     2    centers,boxsize(1,ilev+1),nboxes,ilev,ilevel,iparent,nchild,
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
     1                   1.05*boxsize(1,ilev+1)).and.
     2                   (abs(centers(2,kbox)-centers(2,ibox)).le.
     3                   1.05*boxsize(2,ilev+1))) then
                     
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

      call computecoll2d(nlevels,nboxes,laddr, boxsize,
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
c    nboxes         in: integer *8
c                   number of boxes
c
c    nd             in: integer *8
c                   number of real value functions
c
c    npbox          in: integer *8
c                   number of grid points per function
c
c    centers        in/out: double precision(2,nboxes)
c                   x and y coordinates of the center of boxes
c
c    nlevels        in: integer *8
c                   Number of levels in the tree
c
c    laddr          in/out: integer *8(2,0:nlevels)
c                   boxes at level i are numbered between
c                   laddr(1,i) to laddr(2,i)
c
c    laddrtail      in: integer *8(2,0:nlevels)
c                   new boxes to be added to the tree
c                   structure are numbered from
c                   laddrtail(1,i) to laddrtail(2,i)
c
c    ilevel      in/out: integer *8(nboxes)
c                ilevel(i) is the level of box i
c
c    iparent     in/out: integer *8(nboxes)
c                 iparent(i) is the parent of box i
c
c    nchild      in/out: integer *8(nboxes)
c                nchild(i) is the number of children 
c                of box i
c
c    ichild       in/out: integer *8(8,nboxes)
c                 ichild(j,i) is the jth child of box i
c
c    iflag        in/out: integer *8(nboxes)
c                 iflag(i) is a flag for box i required to generate
c                 level restricted tree from adaptive tree

      implicit none
c     Calling sequence variables and temporary variables
      integer *8 nboxes,nlevels,npbox,nd
      double precision centers(2,nboxes)
      integer *8 laddr(2,0:nlevels), tladdr(2,0:nlevels)
      integer *8 laddrtail(2,0:nlevels)
      integer *8 ilevel(nboxes)
      integer *8 iparent(nboxes)
      integer *8 nchild(nboxes)
      integer *8 ichild(4,nboxes)
      integer *8 iflag(nboxes)
      double precision fvals(nd,npbox,nboxes)
      
      integer *8, allocatable :: tilevel(:),tiparent(:),tnchild(:)
      integer *8, allocatable :: tichild(:,:),tiflag(:)
      integer *8, allocatable :: iboxtocurbox(:),ilevptr(:),ilevptr2(:)

      double precision, allocatable :: tfvals(:,:,:),tcenters(:,:)



c     Temporary variables
      integer *8 i,j,k,l
      integer *8 ibox,ilev, curbox,idim,nblev

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
c      curlev         in: integer *8
c                     the level for which boxes need to be processed
c
c      nboxes         in: integer *8
c                     total number of boxes
c
c      nlevels        in: integer *8
c                     total number of levels
c
c      laddr          in: integer *8(2,0:nlevels)
c                     boxes from laddr(1,ilev) to laddr(2,ilev)
c                     are at level ilev
c
c      nchild         in: integer *8(nboxes)
c                     nchild(ibox) is the number of children
c                     of box ibox
c
c      ichild         in: integer *8(4,nboxes)
c                     ichild(j,ibox) is the box id of the jth
c                     child of box ibox
c
c      nnbors         in: integer *8(nboxes)
c                     nnbors(ibox) is the number of colleagues
c                     of box ibox
c
c      nbors          in: integer *8(9,nboxes)
c                     nbors(j,ibox) is the jth colleague of box
c                     ibox
c
c      centers        in: double precision(2,nboxes)
c                     x and y coordinates of the box centers
c
c      boxsize        in: double precision(0:nlevels)
c                     boxsize(i) is the size of the box at level i
c
c      iflag          in/out: integer *8(nboxes)
c                     iflag(ibox)=3 if it is flag++. iflag(ibox) =1
c                     or 0 at the end of routine depending on
c                     whether box needs to be subdivided or not
c
      implicit none
c     Calling sequence variables
      integer *8 curlev, nboxes, nlevels
      integer *8 laddr(2,0:nlevels),nchild(nboxes),ichild(4,nboxes)
      integer *8 nnbors(nboxes), nbors(9,nboxes)
      integer *8 iflag(nboxes)
      double precision centers(2,nboxes),boxsize(2,0:nlevels)

c     Temporary variables
      integer *8 i,j,k,l,ibox,jbox,kbox,lbox, ict
      double precision xdistest,ydistest,xdis,ydis

      xdistest = 1.05d0*(boxsize(1,curlev) + boxsize(1,curlev+1))/2.0d0
      ydistest = 1.05d0*(boxsize(2,curlev) + boxsize(2,curlev+1))/2.0d0
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
                        if(abs(xdis).le.xdistest) ict = ict + 1
                        if(abs(ydis).le.ydistest) ict = ict + 1
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
      integer *8 nd,npbox,nboxes
      real *8 fvals(nd,npbox,nboxes)
      integer *8 nbloc,nbctr,nlctr
      real *8 centers(2,nboxes),bs(2),grid(2,npbox),xy(2)
      integer *8 ilevel(nboxes),iparent(nboxes)
      integer *8 ichild(4,nboxes),nchild(nboxes)
      integer *8 iflag(nboxes)
      integer *8 ifirstbox,ilastbox
      integer *8, allocatable :: isum(:)
      real *8 dpars(*)
      integer *8 ipars(*)
      complex *16 zpars

      integer *8 i,ibox,nel,j,l,jbox,nbl,ii
      integer *8 xind(4),yind(4)

      real *8 bsh(2)
      data xind/-1,1,-1,1/
      data yind/-1,-1,1,1/

      external fun


      ilastbox = ifirstbox+nbloc-1
      bsh(1) = bs(1)/2
      bsh(2) = bs(2)/2

      allocate(isum(nbloc))

      call cumsum_nz_surf(nbloc,iflag(ifirstbox),isum)
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,j,jbox,nbl,l)
C$OMP$PRIVATE(xy)
      do ibox = ifirstbox,ilastbox
        if(iflag(ibox).gt.0) then
          nchild(ibox) = 4
          nbl = nbctr + (isum(ibox-ifirstbox+1)-1)*4
          do j=1,4
            jbox = nbl+j
            centers(1,jbox) = centers(1,ibox)+xind(j)*bsh(1)
            centers(2,jbox) = centers(2,ibox)+yind(j)*bsh(2)
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
c
      subroutine cumsum_surf(n,a,b)
c
c        this subroutine computes the cumulative sum of an array
c
c
c       TODO: need to make this routine openmp
c
c       input:
c         n - number of elements
c         a - integer *8(n)
c              input array
c
c       output:
c         b - integer *8(n)
c            b(i) = sum_{j=1}^{i} a(j)
c
      implicit none
      integer *8 a(n),b(n),n,i,isum
      isum = 0


      do i=1,n
        isum = isum + a(i)
        b(i) = isum
      enddo
      
      return
      end
c
c
c
c
c
      subroutine cumsum_nz_surf(n,a,b)
c
c        this subroutine computes the cumulative sum of postive
c        entries of an array
c
c
c       TODO: need to make this routine openmp
c
c       input:
c         n - number of elements
c         a - integer *8(n)
c              input array
c
c       output:
c         b - integer *8(n)
c            b(i) = sum_{j=1}^{i} I_{a(j)>0}
c
      implicit none
      integer *8 a(n),b(n),n,i,isum

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

      subroutine computecoll2d(nlevels,nboxes,laddr,boxsize,
     1                       centers,iparent,nchild,ichild,iper,
     2                       nnbors,nbors)

c     This subroutine computes the colleagues for an adaptive
c     pruned tree. box j is a colleague of box i, if they share a
c     vertex or an edge and the two boxes are at the same
c     level in the tree
c
c     INPUT arguments
c     nlevels     in: integer *8
c                 Number of levels
c
c     nboxes      in: integer *8
c                 Total number of boxes
c
c     laddr       in: integer *8(2,0:nlevels)
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
c     iparent     in: integer *8(nboxes)
c                 iparent(i) is the box number of the parent of
c                 box i
c
c     nchild      in: integer *8(nboxes)
c                 nchild(i) is the number of children of box i
c
c     ichild      in: integer *8(4,nboxes)
c                 ichild(j,i) is the box id of the jth child of
c                 box i
c
c     iper        in: integer *8
c                 flag for periodic implementations. 
c                 Currently not used. Feature under construction.
c
c----------------------------------------------------------------
c     OUTPUT
c     nnbors      out: integer *8(nboxes)
c                 nnbors(i) is the number of colleague boxes of
c                 box i
c
c     nbors       out: integer *8(9,nboxes)
c                 nbors(j,i) is the box id of the jth colleague
c                 box of box i
c---------------------------------------------------------------
      implicit none
      integer *8 nlevels,nboxes
      integer *8 iper
      integer *8 laddr(2,0:nlevels)
      double precision boxsize(2,0:nlevels)
      double precision centers(2,nboxes)
      integer *8 iparent(nboxes), nchild(nboxes), ichild(4,nboxes)
      integer *8 nnbors(nboxes)
      integer *8 nbors(9,nboxes)

c     Temp variables
      integer *8 ilev,ibox,jbox,kbox,dad
      integer *8 i,j,ifirstbox,ilastbox


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
c     End of computing colleagues of box i
          enddo
C$OMP END PARALLEL DO         
       enddo

       return
       end
c     
c     
c     
c     
c     


      subroutine fun(nd,xy,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8 nd,ipars(100)
      real *8 dpars(1000),f(nd),xy(2)
      complex *16 zpars
      
      call fungeo(nd,xy,dpars,zpars,ipars,f)

      end

      subroutine fungeo(nd,xy,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8 nd,ipars(100),igeom
      real *8 dpars(1000),f(nd),xy(2)
      complex *16 zpars

      igeom = ipars(1)
      
      if(igeom.eq.0) then
         call funellipsoid(nd,xy,dpars,zpars,ipars,f)
      elseif(igeom.eq.1) then
         call funwtorus(nd,xy,dpars,zpars,ipars,f)
      elseif(igeom.eq.2) then
c     call funstell(nd,xy,dpars,zpars,ipars,f)
         call funtensorfourier_mean(nd,xy,dpars,zpars,ipars,f,ipars(4))
      elseif(igeom.eq.3) then
         call fun_torus_curve_tube(nd,xy,dpars,zpars,ipars,f)
      elseif(igeom.eq.4) then
         call fun_fourier_curve_tube(nd,xy,dpars,zpars,ipars,f)
      elseif(igeom.eq.5) then
c     call fun_cruller(nd,xy,dpars,zpars,ipars,f)
         call fun_cruller_mean(nd,xy,dpars,zpars,ipars,f)         
      elseif(igeom.eq.6) then
         call funaxisymstar(nd,xy,dpars,zpars,ipars,f)
      elseif(igeom.eq.7) then
         call funaxisymsaw(nd,xy,dpars,zpars,ipars,f)
      else
c         print *, "No geometry selected"
      endif
      end


      subroutine funwtorus(nd,st,dpars,zpars,ipars,f)

      implicit none
      
      integer *8 nd,ipars(*)
      real *8 dpars(*),f(nd),st(2),radii(3),s,t,pi
      real *8 xyz(3),nosc,dxyzdst(3,2),scales(3)
      complex *16 zpars

      pi = 4*datan(1.0d0)
            
c     assume that xy is in [-1/2,1/2]^2
c     scale to s,t in [0,2\pi]\times[0,2\pi]
      
      s = (st(1) + 0.5d0)*pi*2.0d0
      t = (st(2) + 0.5d0)*pi*2.0d0
      
c     rminor = dpars(1)

c     rmajor = dpars(2)

      radii(1) = dpars(1)
      radii(2) = dpars(2)
      radii(3) = dpars(3)
      scales(1) = dpars(4)
      scales(2) = dpars(5)
      scales(3) = dpars(6)
      nosc = dpars(7)
     
      call wtorus_eval(s, t, radii, scales, nosc, xyz, dxyzdst)
      
      f(1) = xyz(1)
      f(2) = xyz(2)
      f(3) = xyz(3)

      if (nd.gt.3) then
         f(4) = dxyzdst(1,1)
         f(5) = dxyzdst(2,1)
         f(6) = dxyzdst(3,1)
         
         f(7) = dxyzdst(1,2)
         f(8) = dxyzdst(2,2)
         f(9) = dxyzdst(3,2)
         
      endif      
      
      return
      
      end

      
      subroutine funtensorfourier(nd,st,dpars,zpars,ipars,f,m)

      implicit none
      
      integer *8  nd,ipars(*),m,nfp,nij,k,i,j,idx
      real *8 st(2),s,t,dpars(*),f(nd),xyz(3),dxyzdst(3,2)
      real *8 coefs(2*m+1,2*m+1,3),scales(3),pi
      complex *16 zpars

      pi = 4*datan(1.0d0)

!     scale from [-1/2,1/2]^2 to [0,2*pi]^2. Will be removed.
      s = (st(1) + 0.5d0)*pi*2.0d0
      t = (st(2) + 0.5d0)*pi*2.0d0
      
      nfp = ipars(3)
      nij = 2*m+1

      scales(1) = dpars(nij*nij*3+1)
      scales(2) = dpars(nij*nij*3+2)
      scales(3) = dpars(nij*nij*3+3)
     
      call  xyz_tensor_fourier_eval(s, t, dpars, m, nfp, scales, 
     1     xyz, dxyzdst)

      f(1) = xyz(1)
      f(2) = xyz(2)
      f(3) = xyz(3)

      if (nd.gt.3) then
         
         f(4) = dxyzdst(1,1)
         f(5) = dxyzdst(2,1)
         f(6) = dxyzdst(3,1)
         
         f(7) = dxyzdst(1,2)
         f(8) = dxyzdst(2,2)
         f(9) = dxyzdst(3,2)
         
      endif

      return
      
      end

      subroutine fun_torus_curve_tube(nd,st,dpars,zpars,ipars,f)

      implicit none
      
      integer *8  nd,ipars(*)
      real *8 st(2),s,t,dpars(*),f(nd),xyz(3),dxyzdst(3,2)
      real *8 pi
      complex *16 zpars

      pi = 4*datan(1.0d0)

!     scale from [-1/2,1/2]^2 to [0,2*pi]^2. Will be removed.
      s = (st(1) + 0.5d0)*pi*2.0d0
      t = (st(2) + 0.5d0)*pi*2.0d0
      
      call  torus_curve_tube_eval(s, t, dpars, xyz, dxyzdst)

      f(1) = xyz(1)
      f(2) = xyz(2)
      f(3) = xyz(3)

      if (nd.gt.3) then
         
         f(4) = dxyzdst(1,1)
         f(5) = dxyzdst(2,1)
         f(6) = dxyzdst(3,1)
         
         f(7) = dxyzdst(1,2)
         f(8) = dxyzdst(2,2)
         f(9) = dxyzdst(3,2)
         
      endif

      return
      
      end

      subroutine fun_fourier_curve_tube(nd,st,dpars,zpars,ipars,f)

      implicit none
      
      integer *8  nd,ipars(*)
      real *8 st(2),s,t,dpars(*),f(nd),xyz(3),dxyzdst(3,2)
      real *8 pi
      complex *16 zpars

      pi = 4*datan(1.0d0)

!     scale from [-1/2,1/2]^2 to [0,2*pi]^2. Will be removed.
      s = (st(1) + 0.5d0)*pi*2.0d0
      t = (st(2) + 0.5d0)*pi*2.0d0
      
      call  fourier_curve_tube_eval(s, t, dpars, xyz, dxyzdst)

      f(1) = xyz(1)
      f(2) = xyz(2)
      f(3) = xyz(3)

      if (nd.gt.3) then
         
         f(4) = dxyzdst(1,1)
         f(5) = dxyzdst(2,1)
         f(6) = dxyzdst(3,1)
         
         f(7) = dxyzdst(1,2)
         f(8) = dxyzdst(2,2)
         f(9) = dxyzdst(3,2)
         
      endif

      return
      
      end

      subroutine fun_cruller_mean(nd,st,dpars,zpars,ipars,f)

      implicit none
      
      integer *8 nd,ipars(*)
      real *8 dpars(*),f(nd),st(2),radii(5),s,t,pi
      real *8 xyz(3),dxyzdst(3,2),scales(3)
      real *8 rmaj,rmean,ramp,freqs,freqt
      real *8 cs,ss,ct,stt,phi,h,rr
      real *8 dhds,dhdt,drrds,drrdt
      real *8 d2hdss,d2hdst,d2hdtt
      real *8 d2rrdss,d2rrdst,d2rrdtt
      real *8 xs(3),xt(3),xss(3),xst(3),xtt(3)
      real *8 nvec(3),nrm
      real *8 e1,f1,g1,e2,f2,g2,hmean,denom
      complex *16 zpars(*)

      pi = 4*datan(1.0d0)
            
c     assume that st is in [-1/2,1/2]^2
c     scale to s,t in [0,2*pi] x [0,2*pi]
      
      s = (st(1) + 0.5d0)*pi*2.0d0
      t = (st(2) + 0.5d0)*pi*2.0d0
      
      radii(1) = dpars(1)
      radii(2) = dpars(2)
      radii(3) = dpars(3)
      radii(4) = dpars(4)
      radii(5) = dpars(5)
      
      scales(1) = dpars(6)
      scales(2) = dpars(7)
      scales(3) = dpars(8)
           
      call cruller_eval(s,t,radii,scales,xyz,dxyzdst)
      
      f(1) = xyz(1)
      f(2) = xyz(2)
      f(3) = xyz(3)

      if (nd.ge.9) then
         f(4) = dxyzdst(1,1)
         f(5) = dxyzdst(2,1)
         f(6) = dxyzdst(3,1)
         
         f(7) = dxyzdst(1,2)
         f(8) = dxyzdst(2,2)
         f(9) = dxyzdst(3,2)
      endif

      if (nd.ge.10) then
         rmaj  = radii(1)
         rmean = radii(2)
         ramp  = radii(3)
         freqs = radii(4)
         freqt = radii(5)

         cs  = cos(s)
         ss  = sin(s)
         ct  = cos(t)
         stt = sin(t)

         phi = freqs*s + freqt*t

         h  = rmean + ramp*cos(phi)
         rr = rmaj + h*ct

         dhds = -ramp*freqs*sin(phi)
         dhdt = -ramp*freqt*sin(phi)

         drrds = dhds*ct
         drrdt = dhdt*ct - h*stt

         xs(1) = dxyzdst(1,1)
         xs(2) = dxyzdst(2,1)
         xs(3) = dxyzdst(3,1)

         xt(1) = dxyzdst(1,2)
         xt(2) = dxyzdst(2,2)
         xt(3) = dxyzdst(3,2)

         d2hdss = -ramp*freqs*freqs*cos(phi)
         d2hdst = -ramp*freqs*freqt*cos(phi)
         d2hdtt = -ramp*freqt*freqt*cos(phi)

         d2rrdss = d2hdss*ct
         d2rrdst = d2hdst*ct - dhds*stt
         d2rrdtt = d2hdtt*ct - 2.0d0*dhdt*stt - h*ct

         xss(1) = scales(1)*(d2rrdss*cs - 2.0d0*drrds*ss - rr*cs)
         xss(2) = scales(2)*(d2rrdss*ss + 2.0d0*drrds*cs - rr*ss)
         xss(3) = scales(3)*(d2hdss*stt)

         xst(1) = scales(1)*(d2rrdst*cs - drrdt*ss)
         xst(2) = scales(2)*(d2rrdst*ss + drrdt*cs)
         xst(3) = scales(3)*(d2hdst*stt + dhds*ct)

         xtt(1) = scales(1)*(d2rrdtt*cs)
         xtt(2) = scales(2)*(d2rrdtt*ss)
         xtt(3) = scales(3)*(d2hdtt*stt + 2.0d0*dhdt*ct - h*stt)

         e1 = xs(1)*xs(1) + xs(2)*xs(2) + xs(3)*xs(3)
         f1 = xs(1)*xt(1) + xs(2)*xt(2) + xs(3)*xt(3)
         g1 = xt(1)*xt(1) + xt(2)*xt(2) + xt(3)*xt(3)

         nvec(1) = xs(2)*xt(3) - xs(3)*xt(2)
         nvec(2) = xs(3)*xt(1) - xs(1)*xt(3)
         nvec(3) = xs(1)*xt(2) - xs(2)*xt(1)

         nrm = sqrt(nvec(1)*nvec(1) + nvec(2)*nvec(2) +
     1              nvec(3)*nvec(3))

         nvec(1) = nvec(1)/nrm
         nvec(2) = nvec(2)/nrm
         nvec(3) = nvec(3)/nrm

         e2 = xss(1)*nvec(1) + xss(2)*nvec(2) + xss(3)*nvec(3)
         f2 = xst(1)*nvec(1) + xst(2)*nvec(2) + xst(3)*nvec(3)
         g2 = xtt(1)*nvec(1) + xtt(2)*nvec(2) + xtt(3)*nvec(3)

         denom = 2.0d0*(e1*g1 - f1*f1)
         hmean = (e2*g1 - 2.0d0*f2*f1 + g2*e1)/denom

         f(10) = hmean
      endif
      
      return
      
      end

      
      subroutine fun_cruller(nd,st,dpars,zpars,ipars,f)

      implicit none
      
      integer *8 nd,ipars(*)
      real *8 dpars(*),f(nd),st(2),radii(5),s,t,pi
      real *8 xyz(3),nosc,dxyzdst(3,2),scales(3)
      complex *16 zpars

      pi = 4*datan(1.0d0)
            
c     assume that xy is in [-1/2,1/2]^2
c     scale to s,t in [0,2\pi]\times[0,2\pi]
      
      s = (st(1) + 0.5d0)*pi*2.0d0
      t = (st(2) + 0.5d0)*pi*2.0d0
      
c     rminor = dpars(1)

c     rmajor = dpars(2)

      radii(1) = dpars(1)
      radii(2) = dpars(2)
      radii(3) = dpars(3)
      radii(4) = dpars(4)
      radii(5) = dpars(5)
      
      scales(1) = dpars(6)
      scales(2) = dpars(7)
      scales(3) = dpars(8)
           
      call cruller_eval(s, t, radii, scales, xyz, dxyzdst)
      
      f(1) = xyz(1)
      f(2) = xyz(2)
      f(3) = xyz(3)

      if (nd.eq.9) then
         f(4) = dxyzdst(1,1)
         f(5) = dxyzdst(2,1)
         f(6) = dxyzdst(3,1)
         
         f(7) = dxyzdst(1,2)
         f(8) = dxyzdst(2,2)
         f(9) = dxyzdst(3,2)         
      endif      
      
      return
      
      end

     
      subroutine fungeo2(nd,xy,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8 nd,ipars(100),igeom
      real *8 dpars(1000),f(nd),xy(2)
      complex *16 zpars

      igeom = ipars(1)
      
      if(igeom.eq.0) then
         call funellipsoid(nd,xy,dpars,zpars,ipars,f)
      elseif(igeom.eq.1) then
         call funwtorus2(nd,xy,dpars,zpars,ipars,f)
      elseif(igeom.eq.2) then
c     call funstell(nd,xy,dpars,zpars,ipars,f)
         call funtensorfourier2(nd,xy,dpars,zpars,ipars,f)
      elseif(igeom.eq.3) then
         call funtrefoil(nd,xy,dpars,zpars,ipars,f)
      elseif(igeom.eq.4) then
         call funcruller2(nd,xy,dpars,zpars,ipars,f)
      elseif(igeom.eq.5) then
         call funaxisymstar(nd,xy,dpars,zpars,ipars,f)
      elseif(igeom.eq.6) then
         call funaxisymsaw(nd,xy,dpars,zpars,ipars,f)
      else
c         print *, "No geometry selected"
      endif
      end

      subroutine funellipsoid(nd,xy,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8 nd,ipars(100)
      real *8 dpars(1000),f(nd),xy(2),s,t,pi
      complex *16 zpars

      pi = 4*datan(1.0d0)
      
c     assume that xy is in [-1/2,1/2]^2
c     scale to s,t in [0,2\pi]\times[-pi/2,\pi/2]
      
      s = (xy(1) + 0.5d0)*pi*2.0d0
      t = (xy(2) + 0.5d0)*pi-pi/2

      f(1) =  dpars(1)*dsin(t)*dcos(s)+dpars(4)
      f(2) =  dpars(1)*dsin(t)*dsin(s)+dpars(5)
      f(3) =  dpars(1)*dcos(t)+dpars(6)
      
      end      
      
      
      subroutine funwtorus2(nd,xy,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8 nd,ipars(100),nosc
      real *8 dpars(1000),f(nd),xy(2),s,t,pi
      real *8 rx,ry,rz,cx,cy,cz,x,y,z
      complex *16 zpars

      pi = 4*datan(1.0d0)

      print *,"IN HERE W TORUS"
      
c     assume that xy is in [-1/2,1/2]^2
c     scale to s,t in [0,2\pi]\times[0,2\pi]
      
      s = (xy(1) + 0.5d0)*pi*2.0d0
      t = (xy(2) + 0.5d0)*pi*2.0d0

c     rminor = dpars(1)

c     rmajor = dpars(2)
      rminor = dpars(1)
      rmajor = dpars(2)
      rwave = dpars(3)
      
      a = dpars(4)
      b = dpars(5)
      c = dpars(6)

      cx = dpars(7)
      cy = dpars(8)
      cz = dpars(9)

      rx = dpars(10)
      ry = dpars(11)
      rz = dpars(12)

c     nosc = dpars(13)
      nosc = ipars(3)
      
      rr = rmajor+rminor*cos(t)+rwave*cos(nosc*s)

      x = a*rr*cos(s)
      y = b*rr*sin(s)
      z = c*rminor*sin(t)

      f(1) = x 
      f(2) = y
      f(3) = z

      print *,"ND = ", nd
      print *,"ZPARS = ", zpars
      if (nd.gt.3) then
         
         drrds = -nosc*rwave*sin(nosc*s)
         drrdt = -rminor*sin(t)

         dxds = a*(drrds*cos(s) - rr*sin(s))
         dyds = b*(drrds*sin(s) + rr*cos(s))
         dzds = 0.0d0      
         
         dxdt = a*drrdt*cos(s)
         dydt = b*drrdt*sin(s)
         dzdt = c*rminor*cos(t)

         f(4) = dxds
         f(5) = dyds
         f(6) = dzds
         f(7) = dxdt
         f(8) = dydt
         f(9) = dzdt
         
      endif
      
      return
      
      end


      subroutine funstell(nd,xy,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8 nd,ipars(100),m,n
      real *8 dpars(1000),f(nd),xy(2),s,t,pi
c     real *8, allocatable :: deltas(:,:)
      complex *16 zpars

      pi = 4*datan(1.0d0)

      m = ipars(2)
      n = ipars(3)

c     assume that xy is in [-1/2,1/2]^2
c     scale to s,t in [0,2\pi]\times[0,2\pi]
      
      s = (xy(1) + 0.5d0)*pi*2.0d0
      t = (xy(2) + 0.5d0)*pi*2.0d0


      f(1) = 0
      f(2) = 0
      f(3) = 0

      
      ct = cos(t)
      st = sin(t) 
      do i = -1,m
         do j = -1,n
            cst = cos((1.0d0-i)*s + j*t)
            sst = sin((1.0d0-i)*s + j*t)
            f(1) = f(1) + ct*dpars((m+2)*(j+1)+i+2)*cst
            f(2) = f(2) + st*dpars((m+2)*(j+1)+i+2)*cst
            f(3) = f(3) + dpars((m+2)*(j+1)+i+2)*sst

         end do
      end do



      end

      subroutine funtensorfourier_mean(nd,st,dpars,zpars,ipars,f,m)

      implicit none
      
      integer *8 nd,ipars(*),m,nfp,nij,k,i,j,idx
      real *8 st(2),s,t,dpars(*),f(nd),xyz(3),dxyzdst(3,2)
      real *8 scales(3),pi,hmean
      real *8 hatxyz(3),dhatxyzds(3),dhatxyzdt(3)
      real *8 d2hatxyzdss(3),d2hatxyzdst(3),d2hatxyzdtt(3)
      real *8 bis(2*m+1),bjs(2*m+1),bdis(2*m+1),bdjs(2*m+1)
      real *8 bddis(2*m+1),bddjs(2*m+1)
      real *8 xs(3),xt(3),xss(3),xst(3),xtt(3)
      real *8 nvec(3),nrm,e1,f1,g1,e2,f2,g2,denom
      real *8 ct,stt,c4t,s4t
      complex *16 zpars(*),zmuls,zmult,ima,zfacs,zfact

      data ima/(0.0d0,1.0d0)/

      pi = 4*datan(1.0d0)

c     scale from [-1/2,1/2]^2 to [0,2*pi]^2. Will be removed.
      s = (st(1) + 0.5d0)*pi*2.0d0
      t = (st(2) + 0.5d0)*pi*2.0d0
      
      nfp = ipars(3)
      nij = 2*m+1

      scales(1) = dpars(nij*nij*3+1)
      scales(2) = dpars(nij*nij*3+2)
      scales(3) = dpars(nij*nij*3+3)
     
      call xyz_tensor_fourier_eval(s,t,dpars,m,nfp,scales,
     1     xyz,dxyzdst)

      f(1) = xyz(1)
      f(2) = xyz(2)
      f(3) = xyz(3)

      if (nd.gt.3) then
         f(4) = dxyzdst(1,1)
         f(5) = dxyzdst(2,1)
         f(6) = dxyzdst(3,1)
         
         f(7) = dxyzdst(1,2)
         f(8) = dxyzdst(2,2)
         f(9) = dxyzdst(3,2)
      endif

      if (nd.gt.9) then

         hatxyz(1) = 0.0d0
         hatxyz(2) = 0.0d0
         hatxyz(3) = 0.0d0

         dhatxyzds(1) = 0.0d0
         dhatxyzds(2) = 0.0d0
         dhatxyzds(3) = 0.0d0

         dhatxyzdt(1) = 0.0d0
         dhatxyzdt(2) = 0.0d0
         dhatxyzdt(3) = 0.0d0

         d2hatxyzdss(1) = 0.0d0
         d2hatxyzdss(2) = 0.0d0
         d2hatxyzdss(3) = 0.0d0

         d2hatxyzdst(1) = 0.0d0
         d2hatxyzdst(2) = 0.0d0
         d2hatxyzdst(3) = 0.0d0

         d2hatxyzdtt(1) = 0.0d0
         d2hatxyzdtt(2) = 0.0d0
         d2hatxyzdtt(3) = 0.0d0

         ct = cos(t)
         stt = sin(t)

         c4t = cos(nfp*t)
         s4t = sin(nfp*t)

         zmuls = cos(s) + ima*sin(s)
         zmult = c4t + ima*s4t

         zfacs = zmuls
         zfact = zmult

         bis(1) = 1.0d0
         bdis(1) = 0.0d0
         bddis(1) = 0.0d0

         bjs(1) = 1.0d0
         bdjs(1) = 0.0d0
         bddjs(1) = 0.0d0

         do i=1,m
            bis(i+1) = real(zmuls)
            bis(i+m+1) = imag(zmuls)

            bdis(i+1) = -i*imag(zmuls)
            bdis(i+m+1) = i*real(zmuls)

            bddis(i+1) = -i*i*real(zmuls)
            bddis(i+m+1) = -i*i*imag(zmuls)

            bjs(i+1) = real(zmult)
            bjs(i+m+1) = imag(zmult)

            bdjs(i+1) = -nfp*i*imag(zmult)
            bdjs(i+m+1) = nfp*i*real(zmult)

            bddjs(i+1) = -(nfp*i)*(nfp*i)*real(zmult)
            bddjs(i+m+1) = -(nfp*i)*(nfp*i)*imag(zmult)

            zmuls = zmuls*zfacs
            zmult = zmult*zfact
         enddo

         do j=1,2*m+1
            do i=1,2*m+1

               idx = i + (j-1)*(2*m+1)

               hatxyz(1) = hatxyz(1) +
     1            dpars(idx)*bis(i)*bjs(j)
               hatxyz(2) = hatxyz(2) +
     1            dpars(idx+nij*nij)*bis(i)*bjs(j)
               hatxyz(3) = hatxyz(3) +
     1            dpars(idx+2*nij*nij)*bis(i)*bjs(j)

               dhatxyzds(1) = dhatxyzds(1) +
     1            dpars(idx)*bdis(i)*bjs(j)
               dhatxyzds(2) = dhatxyzds(2) +
     1            dpars(idx+nij*nij)*bdis(i)*bjs(j)
               dhatxyzds(3) = dhatxyzds(3) +
     1            dpars(idx+2*nij*nij)*bdis(i)*bjs(j)

               dhatxyzdt(1) = dhatxyzdt(1) +
     1            dpars(idx)*bis(i)*bdjs(j)
               dhatxyzdt(2) = dhatxyzdt(2) +
     1            dpars(idx+nij*nij)*bis(i)*bdjs(j)
               dhatxyzdt(3) = dhatxyzdt(3) +
     1            dpars(idx+2*nij*nij)*bis(i)*bdjs(j)

               d2hatxyzdss(1) = d2hatxyzdss(1) +
     1            dpars(idx)*bddis(i)*bjs(j)
               d2hatxyzdss(2) = d2hatxyzdss(2) +
     1            dpars(idx+nij*nij)*bddis(i)*bjs(j)
               d2hatxyzdss(3) = d2hatxyzdss(3) +
     1            dpars(idx+2*nij*nij)*bddis(i)*bjs(j)

               d2hatxyzdst(1) = d2hatxyzdst(1) +
     1            dpars(idx)*bdis(i)*bdjs(j)
               d2hatxyzdst(2) = d2hatxyzdst(2) +
     1            dpars(idx+nij*nij)*bdis(i)*bdjs(j)
               d2hatxyzdst(3) = d2hatxyzdst(3) +
     1            dpars(idx+2*nij*nij)*bdis(i)*bdjs(j)

               d2hatxyzdtt(1) = d2hatxyzdtt(1) +
     1            dpars(idx)*bis(i)*bddjs(j)
               d2hatxyzdtt(2) = d2hatxyzdtt(2) +
     1            dpars(idx+nij*nij)*bis(i)*bddjs(j)
               d2hatxyzdtt(3) = d2hatxyzdtt(3) +
     1            dpars(idx+2*nij*nij)*bis(i)*bddjs(j)

            enddo
         enddo

         xs(1) = dxyzdst(1,1)
         xs(2) = dxyzdst(2,1)
         xs(3) = dxyzdst(3,1)

         xt(1) = dxyzdst(1,2)
         xt(2) = dxyzdst(2,2)
         xt(3) = dxyzdst(3,2)

         xss(1) = (d2hatxyzdss(1)*ct - d2hatxyzdss(2)*stt)
     1      *scales(1)
         xss(2) = (d2hatxyzdss(1)*stt + d2hatxyzdss(2)*ct)
     1      *scales(2)
         xss(3) = d2hatxyzdss(3)*scales(3)

         xst(1) = -(dhatxyzds(1)*stt + dhatxyzds(2)*ct)
         xst(1) = xst(1) + (d2hatxyzdst(1)*ct - d2hatxyzdst(2)*stt)
         xst(1) = xst(1)*scales(1)

         xst(2) = dhatxyzds(1)*ct - dhatxyzds(2)*stt
         xst(2) = xst(2) + (d2hatxyzdst(1)*stt + d2hatxyzdst(2)*ct)
         xst(2) = xst(2)*scales(2)

         xst(3) = d2hatxyzdst(3)*scales(3)

         xtt(1) = d2hatxyzdtt(1)*ct - d2hatxyzdtt(2)*stt
         xtt(1) = xtt(1) - 2.0d0*(dhatxyzdt(1)*stt + dhatxyzdt(2)*ct)
         xtt(1) = xtt(1) - (hatxyz(1)*ct - hatxyz(2)*stt)
         xtt(1) = xtt(1)*scales(1)

         xtt(2) = d2hatxyzdtt(1)*stt + d2hatxyzdtt(2)*ct
         xtt(2) = xtt(2) + 2.0d0*(dhatxyzdt(1)*ct - dhatxyzdt(2)*stt)
         xtt(2) = xtt(2) - (hatxyz(1)*stt + hatxyz(2)*ct)
         xtt(2) = xtt(2)*scales(2)

         xtt(3) = d2hatxyzdtt(3)*scales(3)

         e1 = xs(1)*xs(1) + xs(2)*xs(2) + xs(3)*xs(3)
         f1 = xs(1)*xt(1) + xs(2)*xt(2) + xs(3)*xt(3)
         g1 = xt(1)*xt(1) + xt(2)*xt(2) + xt(3)*xt(3)

         nvec(1) = xs(2)*xt(3) - xs(3)*xt(2)
         nvec(2) = xs(3)*xt(1) - xs(1)*xt(3)
         nvec(3) = xs(1)*xt(2) - xs(2)*xt(1)

         nrm = sqrt(nvec(1)*nvec(1) + nvec(2)*nvec(2) +
     1      nvec(3)*nvec(3))

         nvec(1) = nvec(1)/nrm
         nvec(2) = nvec(2)/nrm
         nvec(3) = nvec(3)/nrm

         e2 = xss(1)*nvec(1) + xss(2)*nvec(2) + xss(3)*nvec(3)
         f2 = xst(1)*nvec(1) + xst(2)*nvec(2) + xst(3)*nvec(3)
         g2 = xtt(1)*nvec(1) + xtt(2)*nvec(2) + xtt(3)*nvec(3)

         denom = 2.0d0*(e1*g1 - f1*f1)
         hmean = (e2*g1 - 2.0d0*f2*f1 + g2*e1)/denom

         f(10) = hmean

      endif

      return
      
      end
      subroutine funtensorfourier2(nd,xy,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer, parameter :: m = 2
      integer *8 nd,ipars(100)
      real *8 dpars(1000),f(nd),xy(2),s,t,pi
c      real *8, allocatable:: bis(:), bjs(:),coefs(:,:,:)
c     real *8, allocatable:: bdis(:), bdjs(:)
      real*8 bis(2*m+1), bjs(2*m+1)
      real*8 bdis(2*m+1), bdjs(2*m+1)
      real*8 coefs(2*m+1,2*m+1,3)
      complex *16 zpars
      complex *16 zmuls, zmult, ima, zfacs, zfact
      integer *8 nfp,ni,nj,idx
      data ima/(0.0d0,1.0d0)/

      real *8 :: dxyzdst(3,2)
      real *8 :: hatxyz(3), dhatxyzds(3), dhatxyzdt(3)
      real *8 :: scales(3)

!
!  This subroutine is the chart for a double fourier toroidal
!  discretization
!
!  \hat(x) = \sum_{i=1}^{2m+1} \sum_{j=1} x_{ij} b_{i}(s) b_{j}(nfp*t)
!  \hat(y) = \sum_{i=1}^{2m+1} \sum_{j=1} y_{ij} b_{i}(s) b_{j}(nfp*t)
!  \hat(z) = \sum_{i=1}^{2m+1} \sum_{j=1} z_{ij} b_{i}(s) b_{j}(nfp*t)
!
!  x(s,t) = (\hat(x) \cos(t) - \hat(y) \sin(t))*scales(1)
!  y(s,t) = (\hat(x) \sin(t) + \hat(y) \cos(t))*scales(2)
!  z(s,t) = \hat(z)*scales(3)
!
!  it returns r = (x,y,z), dr/ds and dr/dt
      
      pi = 4*datan(1.0d0)
      nfp = ipars(3) 
c      m = ipars(4)

      print *,"m = ", ipars(4)
      
c      allocate(bis(2*m+1), bjs(2*m+1), bdis(2*m+1), bdjs(2*m+1))
c      allocate(coefs(2*m+1,2*m+1,3))
c     assume that xy is in [-1/2,1/2]^2
c     scale to s,t in [0,2\pi]\times[0,2\pi]
      
      s = (xy(1) + 0.5d0)*pi*2.0d0
      t = (xy(2) + 0.5d0)*pi*2.0d0

      hatxyz(1:3) = 0

      dxyzdst(1,1) = 0
      dxyzdst(2,1) = 0
      dxyzdst(3,1) = 0

      dhatxyzds(1:3) = 0

      dxyzdst(1,2) = 0
      dxyzdst(2,2) = 0
      dxyzdst(3,2) = 0

      dhatxyzdt(1:3) = 0
  
      ct = cos(t)
      st = sin(t)
      
      c4t = cos(nfp*t)
      s4t = sin(nfp*t)
      zmuls = cos(s) + ima*sin(s)
      zmult = c4t + ima*s4t


      zfacs = zmuls
      zfact = zmult
  
      bis(1) = 1
      bdis(1) = 0
  
      bjs(1) = 1
      bdjs(1) = 0

      
      do i=1,m
        bis(i+1) = real(zmuls)
        bis(i+1+m) = imag(zmuls)

        bdis(i+1) = -i*imag(zmuls)
        bdis(i+m+1) = i*real(zmuls)
     
        bjs(i+1) = real(zmult)
        bjs(i+1+m) = imag(zmult)

        bdjs(i+1) = -nfp*i*imag(zmult)
        bdjs(i+1+m) = nfp*i*real(zmult)


        zmuls = zmuls*zfacs
        zmult = zmult*zfact
      enddo
!
!
      print *,"firstloop done"
c      print *,"ZPARs = ", zpars
      print *,"nd = ", nd
      ni = 2*m+1
      nj = 2*m+1

      idx = 0
      do k = 1,3
         do j = 1,nj
            do i = 1,ni
!     idx = i + (j-1)*ni + (k-1)*ni*nj
               idx = idx + 1
               coefs(i,j,k) = dpars(idx)
c               print *,"coefs(",i,",",j,"k",k,") = ", coefs(i,j,k)
            enddo
         enddo
      enddo
      print *,"IDX = ", idx
      scales(1) = dpars(ni*nj*3+1)
      scales(2) = dpars(ni*nj*3+2)
      scales(3) = dpars(ni*nj*3+3)      

      print *,"SCALES = ", scales
      do j = 1,nj
         do i = 1,ni
            idx = i + (j-1)*ni

            hatxyz(1) = hatxyz(1) + dpars(idx)*bis(i)*bjs(j)
            hatxyz(2) = hatxyz(2) + dpars(idx + ni*nj)*bis(i)*bjs(j)
            hatxyz(3) = hatxyz(3) + dpars(idx + 2*ni*nj)*bis(i)*bjs(j)

c          hatxyz(1) = hatxyz(1) + coefs(i,j,1)*bis(i)*bjs(j)
c          hatxyz(2) = hatxyz(2) + coefs(i,j,2)*bis(i)*bjs(j)
c     hatxyz(3) = hatxyz(3) + coefs(i,j,3)*bis(i)*bjs(j)

          

          dhatxyzds(1) = dhatxyzds(1) + coefs(i,j,1)*bdis(i)*bjs(j)
          dhatxyzds(2) = dhatxyzds(2) + coefs(i,j,2)*bdis(i)*bjs(j)
          dhatxyzds(3) = dhatxyzds(3) + coefs(i,j,3)*bdis(i)*bjs(j)

          dhatxyzdt(1) = dhatxyzdt(1) + coefs(i,j,1)*bis(i)*bdjs(j)
          dhatxyzdt(2) = dhatxyzdt(2) + coefs(i,j,2)*bis(i)*bdjs(j)
          dhatxyzdt(3) = dhatxyzdt(3) + coefs(i,j,3)*bis(i)*bdjs(j)

       enddo
      enddo


       f(1) = (hatxyz(1)*ct - hatxyz(2)*st)*scales(1)
       f(2) = (hatxyz(1)*st + hatxyz(2)*ct)*scales(2)
       f(3) = hatxyz(3)*scales(3)

       print *,"f(1) f(2) f(3)", f(1),f(2),f(3)
       
       if (nd.gt.3) then
          
          f(4) = (dhatxyzds(1)*ct - dhatxyzds(2)*st)*scales(1) 
          f(5) = (dhatxyzds(1)*st + dhatxyzds(2)*ct)*scales(2) 
          f(6) = dhatxyzds(3)*scales(3)

          f(7) = -(hatxyz(1)*st + hatxyz(2)*ct)
          f(7) = f(7) + (dhatxyzdt(1)*ct - dhatxyzdt(2)*st)
          f(7) = f(7)*scales(1)

          f(8) = (hatxyz(1)*ct - hatxyz(2)*st)
          f(8) = f(8) + (dhatxyzdt(1)*st + dhatxyzdt(2)*ct)
          f(8) = f(8)*scales(2)

          f(9) = dhatxyzdt(3)*scales(3)
       endif
       
c     xyz(1) = (hatxyz(1)*ct - hatxyz(2)*st)*scales(1)
c      xyz(2) = (hatxyz(1)*st + hatxyz(2)*ct)*scales(2)
c      xyz(3) = hatxyz(3)*scales(3)
  
c      dxyzdst(1,1) = (dhatxyzds(1)*ct - dhatxyzds(2)*st)*scales(1) 
c      dxyzdst(2,1) = (dhatxyzds(1)*st + dhatxyzds(2)*ct)*scales(2) 
c      dxyzdst(3,1) = dhatxyzds(3)*scales(3)

c      dxyzdst(1,2) = -(hatxyz(1)*st + hatxyz(2)*ct)
c      dxyzdst(1,2) = dxyzdst(1,2) + (dhatxyzdt(1)*ct - dhatxyzdt(2)*st)
c      dxyzdst(1,2) = dxyzdst(1,2)*scales(1)

c      dxyzdst(2,2) = (hatxyz(1)*ct - hatxyz(2)*st)
c      dxyzdst(2,2) = dxyzdst(2,2) + (dhatxyzdt(1)*st + dhatxyzdt(2)*ct)
c      dxyzdst(2,2) = dxyzdst(2,2)*scales(2)

c      dxyzdst(3,2) = dhatxyzdt(3)*scales(3)

c       deallocate(bis,bjs,bdis,bdjs,coefs)
       
      end

      
      subroutine funtrefoil(nd,xy,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8 nd,ipars(100)
      real *8 dpars(1000),f(nd),xy(2),s,t,pi,r
c     real *8, allocatable :: deltas(:,:)
      complex *16 zpars

      pi = 4*datan(1.0d0)

      r = dpars(1)
      a = dpars(2)
      b = dpars(3)
      c = dpars(4)
      d = dpars(5)
      
c     assume that xy is in [-1/2,1/2]^2
c     scale to s,t in [0,2\pi]\times[0,2\pi]
      
      s = (xy(1) + 0.5d0)*pi*2.0d0
      t = (xy(2) + 0.5d0)*pi*2.0d0

      f(1) = (a + b*dcos(c*s))*dcos(d*s) + r*dcos(c*s)*dcos(d*s)*dcos(t)
     1     + (dsqrt(2.0d0)*r*(d*(a + b*dcos(c*s))*dcos(d*s)*dsin(c*s)
     2     - b*c*dsin(d*s))*dsin(t))/dsqrt(2.0d0*(a**2.0d0)*(d**2.0d0)
     3     + (b**2.0d0)*(2.0d0*(c**2.0d0) + d**2.0d0)
     4     + b*(d**2.0d0)*(4*a*dcos(c*s) + b*dcos(2.0d0*c*s)))
      
      f(2) = (a + b*dcos(c*s))*dsin(d*s) + r*dcos(c*s)*dcos(t)*dsin(d*s)
     1     + (dsqrt(2.0d0)*r*(b*c*dcos(d*s)
     2     + d*(a + b*dcos(c*s))*dsin(c*s)*dsin(d*s))*dsin(t))
     3     /dsqrt(2.0d0*a**2.0d0*d**2.0d0 + b**2.0d0*(2.0d0*c**2.0d0
     4     + d**2.0d0) + b*d**2.0d0*(4*a*dcos(c*s) + b*dcos(2.0d0*c*s)))
      
      f(3) = b*dsin(c*s) + r*dcos(t)*dsin(c*s)
     1     - (dsqrt(2.0d0)*d*r*dcos(c*s)*(a + b*dcos(c*s))*dsin(t))
     2     /dsqrt(2.0d0*a**2.0d0*d**2.0d0 + b**2.0d0*(2.0d0*c**2.0d0
     3     + d**2.0d0) + b*d**2.0d0*(4*a*dcos(c*s) + b*dcos(2.0d0*c*s)))



      end

      subroutine funcruller2(nd,xy,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8 nd,ipars(100)
      real *8 dpars(1000),f(nd),xy(2),s,t,pi
c     real *8, allocatable :: deltas(:,:)
      complex *16 zpars

      pi = 4*datan(1.0d0)
      
c     assume that xy is in [-1/2,1/2]^2
c     scale to s,t in [0,2\pi]\times[0,2\pi]
      
      s = (xy(1) + 0.5d0)*pi*2.0d0
      t = (xy(2) + 0.5d0)*pi*2.0d0

      h = dpars(2) + dpars(3)*dcos(dpars(4)*s+dpars(5)*t)
      f(1) =  (dpars(1)+h*dcos(t))*dcos(s)
      f(2) =  (dpars(1)+h*dcos(t))*dsin(s)
      f(3) =  h*dsin(t)

      end


      subroutine funaxisymstar(nd,xy,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8 nd,ipars(100)
      real *8 dpars(1000),f(nd),xy(2),s,t,pi
c     real *8, allocatable :: deltas(:,:)
      complex *16 zpars

      pi = 4*datan(1.0d0)
      
c     assume that xy is in [-1/2,1/2]^2
c     scale to s,t in [0,2\pi]\times[0,2\pi]
      
      s = (xy(1) + 0.5d0)*pi*2.0d0
      t = (xy(2) + 0.5d0)*pi*2.0d0

      r = dpars(1)
      a = dpars(2)
      b = dpars(3)
      c = dpars(4)

      h = r*(1+a*dcos(-b*s))
      f(1) = (dcos(s)*h + c)*dcos(-t)
      f(2) = (dcos(s)*h + c)*dsin(-t)
      f(3) = dsin(s)*h

      end

      subroutine funaxisymsaw(nd,xy,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8 nd,ipars(100)
      real *8 dpars(1000),f(nd),xy(2),s,t,pi
c     real *8, allocatable :: deltas(:,:)
      complex *16 zpars

      pi = 4*datan(1.0d0)
      
c     assume that xy is in [-1/2,1/2]^2
c     scale to s,t in [0,2\pi]\times[0,2\pi]
      
      s = (xy(1) + 0.5d0)*pi*2.0d0
      t = (xy(2) + 0.5d0)*pi*2.0d0

      r = dpars(1)
      a = dpars(2)
      b = dpars(3)
      c = dpars(4)
      N = dpars(5)
      or = dpars(6)

      f(1) = dcos(-t)*(c+
     1     r*((b+a*dsin(or*N*s))*dcos(or*s+a*dsin(or*N*s))))
      f(2) = dsin(-t)*(c+
     1     r*((b+a*dsin(or*N*s))*dcos(or*s+a*dsin(or*N*s))))
      f(3) = r*((b+a*dsin(or*N*s))*dsin(or*s+a*dsin(or*N*s)))

      end      


      subroutine legval2coefs_2d(nd,norder,fvals,fcoefs,umat)
c     
c     converts function values at 2d Legendre tensor product grid
c     to Legendre expansion coefficients
c     
c     
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      real *8 fvals(nd,norder,norder)
      real *8 fcoefs(nd,norder,norder),umat(norder,norder)
      real *8, allocatable:: fcv(:,:,:)

      allocate(fcv(nd,norder,norder))
      
      do j=1,norder
         do k=1,norder
            do ind=1,nd
               dd=0
               do k1=1,norder
                  dd=dd+umat(k,k1)*fvals(ind,k1,j)
               enddo
               fcv(ind,k,j)=dd
            enddo
         enddo
      enddo

      do j=1,norder
         do k=1,norder
            do ind=1,nd
               dd=0
               do j1=1,norder
                  dd=dd+umat(j,j1)*fcv(ind,k,j1)
               enddo
               fcoefs(ind,k,j)=dd
            enddo
         enddo
      enddo

      return
      end      


      subroutine legetens_ind2pow_2d(ndeg,type,iind2p)
      integer *8 ndeg, iind2p(2,*)
      character type

      integer *8 i, j, ipol


      if (type .eq. 'f' .or. type .eq. 'F') then
         ipol = 0
         do i = 1,ndeg+1
            do j = 1,ndeg+1
               ipol = ipol+1
               iind2p(1,ipol) = j-1
               iind2p(2,ipol) = i-1                
            enddo
         enddo
      else if (type .eq. 't' .or. type .eq. 'T') then
         ipol = 0
         do i = 1,ndeg+1
            do j = 1,ndeg+1+1-i
               ipol = ipol+1
               iind2p(1,ipol) = j-1
               iind2p(2,ipol) = i-1
            enddo
         enddo
      endif

      return
      end
c     


      subroutine legetens_exps_2d(itype,n,type,x,u,ldu,v,ldv,w)
c     input parameters:
c     
c     itype - the type of the calculation to be performed
c     itype=0 means that only the gaussian nodes are 
c     to be constructed. 
c     itype=1 means that only the nodes and the weights 
c     are to be constructed
c     itype=2 means that the nodes, the weights, and
c     the matrices u, v are to be constructed
c     itype=3 only construct u
c     itype=4 only construct v
      
      implicit none
      integer *8 itype, n, ldu, ldv
      character type
      real *8 x(2,*),w(*)
      real *8 u(ldu,*), v(ldv,*)
      real *8 x1d(n), w1d(n), u1d(n,n), v1d(n,n)
      integer *8 i,j,ipt,itype1d,io, jo,ipol
      
      itype1d = 0
      if (itype .ge. 1) then
         itype1d = 1
      endif
      if (itype .ge. 2) then
         itype1d = 2
      endif
      
      call legeexps(itype1d,n,x1d,u1d,v1d,w1d)

      ipt = 0
      do i=1,n
         do j=1,n
            ipt = ipt + 1
            x(1,ipt) = x1d(j)
            x(2,ipt) = x1d(i)               
         enddo
      enddo

      if (itype .ge. 1) then
         ipt = 0
         do i=1,n
            do j=1,n
               ipt = ipt + 1
               w(ipt) = w1d(i)*w1d(j)
            enddo
         enddo
      endif


      if (itype .eq. 2 .or. itype .eq. 3) then
c     construct u from 1d u
         if (type .eq. 'f' .or. type .eq. 'F') then         
            ipt = 0
            do io = 1,n
               do jo = 1,n
                  ipt = ipt + 1
                  ipol = 0
                  do i=1,n
                     do j=1,n
                        ipol = ipol + 1
                        u(ipol,ipt) =  u1d(i,io)*u1d(j,jo)
                     enddo
                  enddo
               enddo
            enddo

         else if (type .eq. 't' .or. type .eq. 'T') then

            ipt = 0
            do io = 1,n
               do jo = 1,n
                  ipt = ipt + 1
                  ipol = 0
                  do i=1,n
                     do j=1,n+1-i
                        ipol = ipol + 1
                        u(ipol,ipt) = u1d(i,io)*u1d(j,jo)
                     enddo
                  enddo
               enddo
            enddo
         endif
      endif
      
      if (itype .eq. 2 .or. itype .eq. 4) then
c     construct v from 1d v
         if (type .eq. 'f' .or. type .eq. 'F') then         
            ipol = 0
            do io = 1,n
               do jo = 1,n
                  ipol = ipol + 1
                  ipt = 0
                  do i=1,n
                     do j=1,n
                        ipt = ipt + 1
                        v(ipt,ipol) =  v1d(i,io)*v1d(j,jo)
                     enddo
                  enddo
               enddo
            enddo

         else if (type .eq. 't' .or. type .eq. 'T') then

            ipol = 0
            do io = 1,n
               do jo = 1,n+1-io
                  ipol = ipol + 1
                  ipt = 0
                  do i=1,n
                     do j=1,n
                        ipt = ipt + 1
                        v(ipt,ipol) = v1d(i,io)*v1d(j,jo)
                     enddo
                  enddo
               enddo
            enddo
         endif
      endif         

      return
      end


c     
c     
c     
c     
c     
      subroutine chebtens_exps_2d(itype,n,type,x,u,ldu,v,ldv,w)
c     input parameters:
c     
c     itype - the type of the calculation to be performed
c     itype=0 means that only the gaussian nodes are 
c     to be constructed. 
c     itype=1 means that only the nodes and the weights 
c     are to be constructed
c     itype=2 means that the nodes, the weights, and
c     the matrices u, v are to be constructed
c     itype=3 only construct u
c     itype=4 only construct v
      
      implicit none
      integer *8 itype, n, ldu, ldv
      character type
      real *8 x(2,*),w(*)
      real *8 u(ldu,*), v(ldv,*)
      real *8 x1d(n), w1d(n), u1d(n,n), v1d(n,n)
      integer *8 i,j,ipt,itype1d,io, jo,ipol
      
      itype1d = 0
      if (itype .ge. 1) then
         itype1d = 1
      endif
      if (itype .ge. 2) then
         itype1d = 2
      endif
      
      call chebexps(itype1d,n,x1d,u1d,v1d,w1d)

      ipt = 0
      do i=1,n
         do j=1,n
            ipt = ipt + 1
            x(1,ipt) = x1d(j)
            x(2,ipt) = x1d(i)               
         enddo
      enddo

      if (itype .ge. 1) then
         ipt = 0
         do i=1,n
            do j=1,n
               ipt = ipt + 1
               w(ipt) = w1d(i)*w1d(j)
            enddo
         enddo
      endif


      if (itype .eq. 2 .or. itype .eq. 3) then
c     construct u from 1d u
         if (type .eq. 'f' .or. type .eq. 'F') then         
            ipt = 0
            do io = 1,n
               do jo = 1,n
                  ipt = ipt + 1
                  ipol = 0
                  do i=1,n
                     do j=1,n
                        ipol = ipol + 1
                        u(ipol,ipt) =  u1d(i,io)*u1d(j,jo)
                     enddo
                  enddo
               enddo
            enddo

         else if (type .eq. 't' .or. type .eq. 'T') then

            ipt = 0
            do io = 1,n
               do jo = 1,n
                  ipt = ipt + 1
                  ipol = 0
                  do i=1,n
                     do j=1,n+1-i
                        ipol = ipol + 1
                        u(ipol,ipt) = u1d(i,io)*u1d(j,jo)
                     enddo
                  enddo
               enddo
            enddo
         endif
      endif
      
      if (itype .eq. 2 .or. itype .eq. 4) then
c     construct v from 1d v
         if (type .eq. 'f' .or. type .eq. 'F') then         
            ipol = 0
            do io = 1,n
               do jo = 1,n
                  ipol = ipol + 1
                  ipt = 0
                  do i=1,n
                     do j=1,n
                        ipt = ipt + 1
                        v(ipt,ipol) =  v1d(i,io)*v1d(j,jo)
                     enddo
                  enddo
               enddo
            enddo

         else if (type .eq. 't' .or. type .eq. 'T') then

            ipol = 0
            do io = 1,n
               do jo = 1,n+1-io
                  ipol = ipol + 1
                  ipt = 0
                  do i=1,n
                     do j=1,n
                        ipt = ipt + 1
                        v(ipt,ipol) = v1d(i,io)*v1d(j,jo)
                     enddo
                  enddo
               enddo
            enddo
         endif
      endif         

      return
      end
c     
c     
      subroutine mesh2d(x,nx,y,ny,xy)
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
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

cccccccccccccccccccccccccccccccccccccc
c     Routines for generating surface
cccccccccccccccccccccccccccccccccccccc

      subroutine surf_mesh_build(norder,iptype,nd,dpars,zpars,ipars,
     1     nlevels,nboxes,ltree,itree,iptr,centers,boxsizein,npatches,
     2     npts,srcvals,srccoefs,norders,ixyzs,iptypes)


      implicit none
      
      integer  ifplot

      integer *8 ibox,ilev,npatches,npmax,norder,nd
      integer *8 nlevels,nboxes,ltree,npts,iptype
      integer *8 itree(ltree),iptr(8)
      integer *8 norders(npatches),ixyzs(npatches+1),iptypes(npatches)
      integer *8 i,npols,ifun
      integer *8, target :: mmax,nmax,ipars(100)
      integer *8, pointer :: iptr1,iptr2,iptr3,iptr4
      
      
      real *8 done,pi,umin,umax,vmin,vmax,centers(2,nboxes)
      real *8 srcvals(12,npts),srccoefs(9,npts)
      real *8 boxsize(2,0:nlevels),boxsizein(2,nlevels+1)

      real *8, target :: dpars(1000)
      real *8, target :: p1(10),p2(10),p3(10),p4(10),nosc
      real *8, pointer :: dptr1,dptr2,dptr3,dptr4

      real *8, target :: dummy_real

      real *8, allocatable, target :: triaskel(:,:,:)
      real *8, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts(:)

      
      complex *16, target :: zpars
      complex *16, pointer :: zptr1
      
      character *100 fname
      character*20 tmp
      
      procedure (), pointer :: xtri_geometry

      external xtri_geo_fun

      ifplot = ipars(5)
      done = 1
      pi = atan(done)*4
      do i=0,nlevels
         boxsize(1,i) = boxsizein(1,i+1)
         boxsize(2,i) = boxsizein(2,i+1)
      enddo


         
      if (iptype.ne.1) then
         print *, "Only iptype=1 is implemented. Exiting"
         return
      endif

        
        allocate(triaskel(3,3,npatches))
        npmax = npatches
        
        done = 1
        pi = atan(done)*4

        umin = 0.0d0
        umax = 2.0d0*pi
        vmin = 0.0d0
        vmax = 2.0d0*pi

        call xtri_surftree_rectmesh(umin,umax,vmin,vmax,
     1       nlevels,nboxes,itree(iptr(2)),itree(iptr(4)),
     2       centers,boxsize,triaskel,npmax)


        dummy_real = 0.0d0
        dptr3 => dummy_real

        xtri_geometry => xtri_geo_fun
        dptr1 => triaskel(1,1,1)
        dptr2 => dpars(1)        
        iptr4 => ipars(1)

c     change to zpars
c        zptr1 => zpars
        
        if(ifplot.eq.1) then
           fname = 'tri_surf_geometry_'
           write(tmp,'(I0)') ipars(1)
           fname = trim(fname)//trim(tmp)// '.vtk'
           call xtri_vtk_surf(fname,npatches,xtri_geometry, dptr1,
     1          dptr2,dptr3,iptr4,norder,
     2          'Triangulated surface of input geometry')
        endif


        npols = (norder + 1) * (norder + 2)/2

        allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols))
        allocate(wts(npols))

        call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)

        call getgeominfo(npatches,xtri_geometry,dptr1,dptr2,dptr3,iptr4,
     1     npols,uvs,umatr,srcvals,srccoefs)

         
         do i=1,npatches
            norders(i) = norder
            ixyzs(i) = 1 +(i-1)*npols
            iptypes(i) = iptype
         enddo

         ixyzs(npatches+1) = 1 + npts
         
      return
      
      end





      
      subroutine xtri_surftree_rectmesh(umin,umax,vmin,vmax,
     1     nlevels,nboxes,ilevel,nchild,
     2     centers,boxsize,triaskel,maxtri)
      
      implicit none

      integer *8 ntri,maxtri,ibox,ilev,nlevels,nboxes
      integer *8 ilevel(nboxes), nchild(nboxes)
      
      real *8 centers(2,nboxes),boxsize(2,0:nlevels)
      real *8 umin,umax,vmin,vmax,width,height
      real *8 triaskel(3,3,maxtri)
      real *8 u0,u1,v0,v1
!     
!     PLEASE NOTE: This code assumes that points in the quad-tree
!     are mapped from [-1/2,1/2]x[-1/2,1/2].
!     
!     this routine returns a mesh consisting of triangles, created by
!     splitting the leaf level rectangles [umin,umax] \times [vmin,vmax]
!     in a quad-tree. Each leaf square is split into two triangels.
!     --- keep in mind the
!     triangles are stored with u,v,w components, with w=0 (in order to
!     be consistent with other skeleton meshes)
!     
!     Input:
!     umin,umax,vmin,vmax -  dimensions of the rectangle
!     [umin,umax]x[vmin,vmax]
!     maxtri - maximum number of allowable triangles
!     
!     Output:
!     ntri - number of triangles created
!     triaskel - skeleton mesh info, basically just vertices of each
!     triangle
!

      width = umax-umin
      height = vmax-vmin
      
      ntri = 0
!     loop over all boxes, if the box is a leaf, split it
      do ibox = 1,nboxes         
         if(nchild(ibox).eq.0) then
!     is leaf, split
            ilev = ilevel(ibox)
            u0 = (centers(1,ibox) - boxsize(1,ilev)/2 + 0.5d0)
     1           *width+umin
            u1 = (centers(1,ibox) + boxsize(1,ilev)/2 + 0.5d0)
     1           *width + umin
            v0 = (centers(2,ibox) - boxsize(2,ilev)/2 + 0.5d0)
     1           *height + vmin
            v1 = (centers(2,ibox) + boxsize(2,ilev)/2 + 0.5d0)
     1           *height + vmin
            call xtri_rectmesh0(u0, u1, v0, v1, triaskel(1,1,ntri+1))
            ntri= ntri + 2
         endif
      enddo
      
      return
      end



      
