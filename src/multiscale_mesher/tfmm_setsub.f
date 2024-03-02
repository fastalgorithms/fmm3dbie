c
c
c
c
      subroutine tfmm3d_setsub(eps, ns, src, srcnorm, wts, nt, targs, 
     1  sigma, sigma_grad, pottarg, gradtarg) 
c----------------------------------------------------
c  This subroutine evaluates the integral 
c    \int_{S} \psi_{\sigma}(x,x') dx 
c
c  where \psi_{\sigma}(x,x') is given by
c
c    -\frac{n(x') \cdot (x-x')}{4\pi} \left(erf(r)/|x-x'|^3 - 
c       \sqrt{2/pi} \exp(-r^2)/(\sigma(x) |x-x'|^2 \right) 
c 
c  with r = |x-x'|/(\sqrt(2) \sigma).

c  When |x-x'|>= 8 \sigma (x) then  the integrand is just the 
c  Laplace double layer potential.
c
c  The integrand is discretized using a smooth quadrature rule,
c  and the far-field is accelerated using an fmm
c
c  Input arguments:
c    - eps: real *8
c        tolerance requested
c    - ns: integer
c        number of sources
c    - src: real *8 (3,ns)
c        discretization points on the surface
c    - srcnorm: real *8 (3,ns)
c        normal vectors at the discretization points
c    - wts: real *8 (ns)
c        quadrature rule for integrating smooth functions
c    - nt: integer
c        number of targets
c    - targs: real *8 (3,nt)
c        target locations
c    - sigma: real *8 (nt)
c        sigma function at the target locations
c    - sigma_grad: real *8 (3,nt)
c        derivative of the sigma function at the target locations
c 
c  Output arguments:
c    - pottarg: real *8 (nt)
c        potential at the target locations
c    - gradtarg: real *8 (3,nt)
c        gradient at the target locations
c
c           
c               
c
      implicit none
      integer ns, npols, nt, ndtarg
      real *8 src(3,ns), srcnorm(3,ns), wts(ns)
      real *8 targs(3,nt), sigma(nt), sigma_grad(3,nt)
      real *8 pottarg(nt), gradtarg(3,nt)
      real *8, allocatable :: potsort(:), gradsort(:,:)
      real *8 eps

      real *8, allocatable :: charges(:), dipvec(:,:)
      integer, allocatable :: iboxtarg(:),iboxsrc(:)
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      real *8 tmp(10),val,vgrad(3)

      integer i,j,jpatch,jquadstart,jstart
      integer ifaddsub

      integer *8 ltree,ipointer(8)
      integer, allocatable :: itree(:)
      integer, allocatable :: il1(:),il2(:),ilint(:),il1m2(:),il2m1(:)
      real *8, allocatable :: boxsize(:),centers(:,:)
      integer, allocatable :: isrcse(:,:),isrcper(:)
      integer, allocatable :: itargse(:,:),itargper(:)

      real *8 expc(3)
      integer ibox,nexpc,idivflag,iert,ifnear,ii,isend,isep,isource
      integer isstart,itarg,itend,itstart,itt,jbox,jpt,mhung,mnbors
      integer iss,l,lpt,mnlist1,mnlist2,mnlist3,mnlist4
      integer n1m2,n2m1,nadd,nbmax,nboxes,nchild,ndiv,nl2,nlevels
      integer nlmax,npover,nl1,ntj
      integer ier,nlmin,iper,ifunif

      integer, allocatable :: nlist1(:),list1(:,:)
      integer, allocatable :: nlist2(:),list2(:,:)
      integer, allocatable :: nlist3(:),list3(:,:)
      integer, allocatable :: nlist4(:),list4(:,:)
      
      complex *16 zdotu
      real *8 pottmp
      real *8, allocatable :: ctmp1(:),ctmp2(:),dtmp1(:,:),
     1   dtmp2(:,:)
      real *8 radexp,epsfmm

      integer ipars
      real *8 dpars,timeinfo(10),t1,t2,omp_get_wtime
      real *8 timeinfo_fmm(6)

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp1(:,:),srctmp2(:,:)
      real *8 thresh,ra
      real *8 over4pi
      integer nss

      integer nd,ntarg0

      real *8 ttot,done,pi
      real *8 distest, rr, xdis, ydis, zdis, rtest
      integer ilev, ilevup, jlev, isrc, lbox, llev, n1, ncoll
      data over4pi/0.07957747154594767d0/

      parameter (nd=1,ntarg0=1)

      done = 1
      pi = atan(done)*4

           
      ifpgh = 0
      ifpghtarg = 2

      allocate(charges(ns), dipvec(3,ns))

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)      
      do i=1,ns
        charges(i) = 0 
        dipvec(1,i) = -srcnorm(1,i)*wts(i)*over4pi
        dipvec(2,i) = -srcnorm(2,i)*wts(i)*over4pi
        dipvec(3,i) = -srcnorm(3,i)*wts(i)*over4pi
      enddo
C$OMP END PARALLEL DO      

      ifcharge = 0
      ifdipole = 1

c
c       setup tree
c
c

      isep = 1
      nlmax = 51
      nbmax = 0
      nlevels = 0
      nboxes = 0
      ltree = 0

      idivflag = 0
      mnlist1 = 0
      mnlist2 = 0
      mnlist3 = 0
      mnlist4 = 0

      call lndiv(eps,ns,nt,ifcharge,ifdipole,ifpgh,ifpghtarg,
     1   ndiv,idivflag) 
c
cc      set tree flags
c 
      nlmax = 51
      nlevels = 0
      nboxes = 0
      ltree = 0
      nlmin = 0 
      iper = 0
      ifunif = 0
      ifnear = 0


c 
cc     memory management code for contructing level restricted tree
      call pts_tree_mem(src, ns, targs, nt, idivflag, ndiv,
     1  nlmin, nlmax, iper, ifunif, nlevels, nboxes, ltree)
      
       allocate(itree(ltree))
       allocate(boxsize(0:nlevels))
       allocate(centers(3,nboxes))

c       Call tree code
      call pts_tree_build(src, ns, targs, nt, idivflag, ndiv,
     1  nlmin, nlmax, iper, ifunif, nlevels, nboxes, ltree, itree,
     2  ipointer, centers, boxsize)
      
      

      allocate(isrcse(2,nboxes),itargse(2,nboxes))
      allocate(isrcper(ns),itargper(nt))

      call pts_tree_sort(ns,src,itree,ltree,nboxes,nlevels,
     1   ipointer,centers,isrcper,isrcse)

      call pts_tree_sort(nt,targs,itree,ltree,nboxes,nlevels,
     1   ipointer,centers,itargper,itargse)


      mnbors = 27
      isep = 1

      call computemnlists(nlevels,nboxes,itree(ipointer(1)),boxsize,
     1  centers,itree(ipointer(3)),itree(ipointer(4)),
     2  itree(ipointer(5)),isep,itree(ipointer(6)),mnbors,
     2  itree(ipointer(7)),iper,mnlist1,mnlist2,mnlist3,mnlist4)
      
      allocate(list1(mnlist1,nboxes),nlist1(nboxes))
      allocate(list2(mnlist2,nboxes),nlist2(nboxes))
      allocate(list3(mnlist3,nboxes),nlist3(nboxes))
      allocate(list4(mnlist4,nboxes),nlist4(nboxes))

      call computelists(nlevels,nboxes,itree(ipointer(1)),boxsize,
     1  centers,itree(ipointer(3)),itree(ipointer(4)),
     2  itree(ipointer(5)),isep,itree(ipointer(6)),mnbors,
     3  itree(ipointer(7)),iper,nlist1,mnlist1,list1,nlist2,
     4  mnlist2,list2,nlist3,mnlist3,list3,
     4  nlist4,mnlist4,list4)


c
c
c       call the fmm
c

      call cpu_time(t1)
C$      t1 = omp_get_wtime()      
      call lfmm3d_ndiv(nd, eps, ns, src, ifcharge, charges,
     1  ifdipole, dipvec, iper, ifpgh, tmp, tmp, tmp, nt, targs,
     1  ifpghtarg, pottarg, gradtarg, tmp, ndiv, idivflag, 
     1  ifnear, timeinfo_fmm, ier)
      call cpu_time(t2)
C$      t2 = omp_get_wtime()

      timeinfo(1) = t2-t1
      
c
c
c    work with sorted potentials and unsort them again later
c
      allocate(potsort(nt))
      allocate(gradsort(3,nt))
      call dreorderf(1,nt,pottarg,potsort,itargper)
      call dreorderf(3,nt,gradtarg,gradsort,itargper)


      thresh = 1.0d-14



c
c    subtract  precomputed near quadrature /setminus list1 
c       also needs to go from pts (targs) -> pts (sources)
c 
c
c    il1 - list of sources in the near field of a target (A)
c    il2 - list of sources in the list1 of the target from fmm
c        perspective (B)
c    il1m2 = A \cap (A \cap B)^{c}
c    il2m1 = B \cap (A \cap B)^{c}
c

     
      allocate(il2(ndiv*mnlist1),il2m1(ndiv*mnlist1))
      allocate(dtmp2(3,ndiv*mnlist1))
      allocate(srctmp2(3,ndiv*mnlist1))
      allocate(srctmp1(3,ns))
      allocate(dtmp1(3,ns))
      allocate(il1(ns),il1m2(ns))


  

      call cpu_time(t1)
C$      t1 = omp_get_wtime()     

C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP$PRIVATE(ibox,nchild,nl2)
C$OMP$PRIVATE(i,jbox,isstart,isend,j,isource,il2)
C$OMP$PRIVATE(itstart,itend,itt,itarg,nl1,il1,il1m2,il2m1)
C$OMP$PRIVATE(jpatch,l,jpt,lpt,n1m2,n2m1,ii,val,vgrad,npover)
C$OMP$PRIVATE(dtmp1,dtmp2,srctmp1,srctmp2)
C$OMP$PRIVATE(isrc,iss,ilev,rtest,ilevup,ncoll,lbox,xdis,ydis)
C$OMP$PRIVATE(zdis,llev,distest,rr)
C$OMP$SCHEDULE(DYNAMIC)
      do ibox = 1,nboxes
        nchild = itree(ipointer(4)+ibox-1)
        ilev = itree(ipointer(2)+ibox-1)
        if(nchild.eq.0) then

c
c     populate il2
c
          nl2 = 0
          do i=1,nlist1(ibox)
            jbox = list1(i,ibox) 
            isstart = isrcse(1,jbox) 
            isend = isrcse(2,jbox)
            do j=isstart,isend
              isource = isrcper(j) 
              nl2 = nl2 + 1
              il2(nl2) = isource
            enddo
          enddo


c
c    end of populating il2.
c    
c    now loop over targets in this box
c
          itstart = itargse(1,ibox)
          itend = itargse(2,ibox) 
          do itt = itstart,itend
            itarg = itargper(itt) 
            
c
c    populate il1 
c
            nl1 = 0
            rtest = 64*sigma(itarg)**2
            ilevup = ceiling(log(8*sigma(itarg)*2/boxsize(ilev))
     1          /log(2.0d0))
            jbox = ibox
            ilevup = min(ilevup,ilev)

            do i=1,ilevup
              jbox = itree(ipointer(3)+jbox-1)
            enddo
c
c   loop over colleagues
c
            ncoll = itree(ipointer(6)+jbox-1)
            do l=1,ncoll
              lbox = itree(ipointer(7)+mnbors*(jbox-1)+l-1)
c
c  check if lbox does not intersect with region of interest
c 
              xdis = abs(centers(1,lbox)-targs(1,itarg))
              ydis = abs(centers(2,lbox)-targs(2,itarg))
              zdis = abs(centers(3,lbox)-targs(3,itarg))

              llev = itree(ipointer(2)+lbox-1)
              distest = 8*sigma(itarg)+boxsize(llev)/2
              if((xdis.le.distest).and.(ydis.le.distest).and.
     1            (zdis.le.distest)) then
                 
                do iss=isrcse(1,lbox),isrcse(2,lbox)
                   isrc = isrcper(iss)
                   rr = (src(1,isrc)-targs(1,itarg))**2 + 
     1                (src(2,isrc)-targs(2,itarg))**2 + 
     2                (src(3,isrc)-targs(3,itarg))**2
                   if(rr.le.rtest) then
                     nl1 = nl1 + 1
                     il1(nl1) = isrc
                   endif
                 enddo
               endif
             enddo

c
c   loop over list1 boxes which are larger
c
            jlev = itree(ipointer(2)+jbox-1)
            do l=1,nlist1(jbox)
              lbox = list1(l,jbox)
              llev = itree(ipointer(2)+lbox-1)
              if(llev.lt.jlev) then
                distest = 8*sigma(itarg)+boxsize(llev)/2
                if((xdis.le.distest).and.(ydis.le.distest).and.
     1              (zdis.le.distest)) then
                 
                  do iss=isrcse(1,lbox),isrcse(2,lbox)
                    isrc = isrcper(iss)
                    rr = (src(1,isrc)-targs(1,itarg))**2 + 
     1                 (src(2,isrc)-targs(2,itarg))**2 + 
     2                 (src(3,isrc)-targs(3,itarg))**2
                    if(rr.le.rtest) then
                      nl1 = nl1 + 1
                      il1(nl1) = isrc
                     endif
                  enddo
                endif
              endif
            enddo

c
c   end of populating il1. now perform various set subtractions
c
            n1m2 = 0
            n2m1 = 0
            call setsub(il1,nl1,il2,nl2,il1m2,n1m2,il2m1,n2m1)

            do i=1,nl1
              isrc = il1(i)
              val = 0
              vgrad(1:3) = 0
              call nearkernel(src(1,isrc), dipvec(1,isrc), 
     1           targs(1,itarg), sigma(itarg), sigma_grad(1,itarg), 
     2           val, vgrad)
              potsort(itt) = potsort(itt) + val
              gradsort(1:3,itt) = gradsort(1:3,itt) + vgrad(1:3)

            enddo




c
c   subtract off il1m2
c
c   gather step
c
            do i=1,n1m2
              ii = il1m2(i)
              srctmp1(1,i) = src(1,ii)
              srctmp1(2,i) = src(2,ii)
              srctmp1(3,i) = src(3,ii)
            enddo

            do i=1,n1m2
              ii = il1m2(i)
              dtmp1(1,i) = dipvec(1,ii)
              dtmp1(2,i) = dipvec(2,ii)
              dtmp1(3,i) = dipvec(3,ii)
            enddo

            val = 0
            vgrad(1:3) = 0

            call l3ddirectdg(nd,srctmp1,dtmp1,
     1        n1m2,targs(1,itarg),ntarg0,val,vgrad,thresh)

c
c  scatter step
c
            potsort(itt) = potsort(itt) - val
            gradsort(1:3,itt) = gradsort(1:3,itt) - vgrad(1:3)

c
c   add il2m1
c
c
c   gather step
c
            do i=1,n2m1
              ii = il2m1(i)
              srctmp2(1,i) = src(1,ii)
              srctmp2(2,i) = src(2,ii)
              srctmp2(3,i) = src(3,ii)
            enddo

            do i=1,n2m1
              ii = il2m1(i)
              dtmp2(1,i) = dipvec(1,ii)
              dtmp2(2,i) = dipvec(2,ii)
              dtmp2(3,i) = dipvec(3,ii)
            enddo

            val = 0
            vgrad(1:3) = 0
            call l3ddirectdg(nd,srctmp2,dtmp2,
     1          n2m1,targs(1,itarg),ntarg0,val,vgrad,thresh)
c
c  scatter step
c
            potsort(itt) = potsort(itt) + val
            gradsort(1:3,itt) = gradsort(1:3,itt) + vgrad(1:3)

          enddo
        endif
      enddo
C$OMP END PARALLEL DO      
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      
      timeinfo(2) = t2-t1

      call dreorderi(1,nt,potsort,pottarg,itargper)
      call dreorderi(3,nt,gradsort,gradtarg,itargper)

      
      ttot = timeinfo(1) + timeinfo(2)
cc      call prin2('timeinfo in setsub=*',timeinfo,2)

        
      
      return
      end

c
c
c
c
      subroutine nearkernel(src, dipvec, targ, sigma, grad_sigma, pot,
     1    grad)
      implicit real *8 (a-h,o-z)
      real *8 src(3), targ(3), dipvec(3), sigma, grad_sigma(3), pot
      real *8 grad(3)

      real *8 r, h, hh, hhh, prod
      real *8 dxyz(3)

      dxyz(1) = src(1) - targ(1)
      dxyz(2) = src(2) - targ(2)
      dxyz(3) = src(3) - targ(3)
      r = sqrt(dxyz(1)*dxyz(1) + dxyz(2)*dxyz(2) + dxyz(3)*dxyz(3))
      call surf_smooth_ker(r, sigma, h, hh, hhh)
      prod = dxyz(1)*dipvec(1) + dxyz(2)*dipvec(2) + dxyz(3)*dipvec(3)
      pot = -h*prod

      
      grad(1:3) = -prod*(-hh*dxyz(1:3) + hhh*grad_sigma(1:3)) + 
     1    h*dipvec(1:3)

      return
      end



      subroutine surf_smooth_ker(r, sgma, h, hh, hhh)
      implicit real *8 (a-h,o-z)
      real *8 dmexp
      real *8, parameter :: pi=3.141592653589793238462643383d0
      real *8, parameter :: fourpi = 12.5663706143592d0
      real *8, parameter :: foursqrtpi3 = 22.2733119873268d0
      real *8, parameter :: threesqrtpi = 5.31736155271655d0
      real *8, parameter :: sqpi = 1.7724538509055159d0
      real *8, parameter :: sq2 = 1.4142135623730951d0

    
      dmexp=exp(-r**2/(2*sgma**2))
      hhh=-sqrt(2/pi)*(1/sgma**4)*dmexp

      uuu = r/(sq2*sgma)
      if (r==0.0d0) then
        h = sq2/(3.0d0*sgma**3*sqpi)
        hh = -sq2/(5.0d0*sgma**5.0d0*sqpi)
      else if (uuu<0.01d0) then
        uuu2 = uuu*uuu
        uuu4 = uuu2*uuu2
        uuu6 = uuu2*uuu4
        uuu8 = uuu2*uuu6
        h = 0.0d0
        hh = 0.0d0
        h = uuu8/132.0d0 - uuu6/27.0d0 + uuu4/7.0d0 - 
     1      2*uuu2/5.0d0 + 2.0d0/3.0d0
        h = h/(sq2*sqpi*(sgma**3))
        hh = -uuu8/78.0d0 + 2*uuu6/33.0d0 - 2*uuu4/9 + 4*uuu2/7 -
     1     4.0d0/5.0d0
        hh = hh/((sq2*sgma)**5)
      else
        dmerf = erf(uuu) 
        r2 = r*r
        r3 = r2*r
        r5 = r2*r3
      
        h = dmerf/r3 - sq2*dmexp/(sqpi*sgma*r2)
        denom = sqpi*r5*sgma**3
        hh = (sq2*r3*dmexp - threesqrtpi*sgma**3*dmerf 
     1          + 3*sq2*r*sgma**2*dmexp)/denom
      endif


      return
      end
