      implicit real *8 (a-h,o-z)

      call test_get_uni(isuccess0,isuccess1)
      call test_setdecomp(isuccess)

      stop
      end




      subroutine test_get_uni(isuccess0,isuccess1)

      implicit real *8 (a-h,o-z)
      integer, allocatable :: a(:),iuni(:),iuniind(:)
      integer, allocatable :: b(:),c(:)
      integer, allocatable :: iuni3(:,:)


      call prini(6,13)

      n = 37
      allocate(a(n),iuni(n),iuniind(n),iuni3(3,n),b(n),c(n))

      nuni = 0

      nmax = n/4

      ifprint = 0

      do i=1,n
        r = hkrand(0)*nmax
        a(i) = (r+0.0d0)
        r = hkrand(0) 
        b(i) = 0
        if(r.gt.0.8d0) b(i) = 1
        r = hkrand(0)
        c(i) = 0
        if(r.gt.0.6d0) c(i) = 1
      enddo

      if(ifprint.ge.1) call prinf('a=*',a,n)
      call get_iuni1_omp(n,a,nuni,iuni,iuniind)

      if(ifprint.ge.1) then
        call prinf('nuni=*',nuni,1)
        call prinf('iuni=*',iuni,nuni)
        call prinf('iuniind=*',iuniind,n)
      endif


      isuccess0=1
      do i=1,n
        if(abs(a(i)-iuni(iuniind(i))).ge.1.0d-16) then 
          print *, "problem in iuni1 with index i=",i
          isuccess0 = 0
        endif
      enddo

      ifprint = 0
      if(ifprint.ge.1) then
        call prinf('a=*',a,n)
        call prinf('b=*',b,n)
        call prinf('c=*',c,n)
      endif

      
      call get_iuni3(n,b,c,a,nuni,iuni3,iuniind)

      if(ifprint.ge.1) then
        call prinf('nuni=*',nuni,1)
        call prinf('iuni3=*',iuni3,3*nuni)
        call prinf('iuniind=*',iuniind,n)
      endif

      isuccess1 = 1
      do i=1,n
        if(abs(b(i)-iuni3(1,iuniind(i))).ge.1.0d-16.or.
     1     abs(c(i)-iuni3(2,iuniind(i))).ge.1.0d-16.or.
     2     abs(a(i)-iuni3(3,iuniind(i))).ge.1.0d-16) then 
          print *, "problem in iuni3 with index i=",i
          isuccess1 = 0
        endif
      enddo


      if(isuccess0.eq.1) call prinf('test iuni1 passed*',i,0)
      if(isuccess1.eq.1) call prinf('test iuni3 passed*',i,0)
      return
      end


      
      subroutine test_setdecomp(isuccess)

      implicit real *8 (a-h,o-z)
      integer, allocatable :: a(:),iauni(:),auni(:)
      integer, allocatable :: b(:),ibuni(:),buni(:)
      integer, allocatable :: aintb(:),iaintb(:),aintbc(:)
      integer, allocatable :: iaintbc(:),binta(:),ibinta(:)
      integer, allocatable :: bintac(:),ibaintac(:)
      integer, allocatable :: itmp1(:),itmp2(:),aunisort(:),bunisort(:)


      call prini(6,13)

      n = 31
      m = 37
      allocate(a(n),auni(n),iauni(n),b(m),buni(m),ibuni(m))
      allocate(aunisort(n),bunisort(m))

      nuni = 0

      nmax = n*1.2
      mmax = m*1.2

      ifprint = 1

      do i=1,n
        a(i) = hkrand(0)*nmax
      enddo

      do i=1,m
        b(i) = hkrand(0)*mmax
      enddo

      call get_iuni1_omp(n,a,nauni,aunisort,iauni)
      call get_iuni1_omp(m,b,nbuni,bunisort,ibuni)


c
c
c    randomly scramble aunisort,bunisort arrays
c

      allocate(itmp1(nauni),itmp2(nbuni))
      call sorti(nauni,a,itmp1)
      call sorti(nbuni,b,itmp2)

      do i=1,nauni
        auni(i) = aunisort(itmp1(i))
      enddo

      do i=1,nbuni
        buni(i) = bunisort(itmp2(i))
      enddo

      if(ifprint.ge.1) then
        call prinf('nauni=*',nauni,1)
        call prinf('auni=*',auni,nauni)
        call prinf('nbuni=*',nbuni,1)
        call prinf('buni=*',buni,nbuni)
      endif


      allocate(aintb(nauni),iaintb(nauni),aintbc(nauni),iaintbc(nauni))
      call setdecomp(nauni,auni,nbuni,buni,naintb,aintb,iaintb,naintbc,
     1      aintbc,iaintbc)
      
      if(ifprint.ge.1) then
        call prinf('naintb=*',naintb,1)
        call prinf('aintb=*',aintb,naintb)
        call prinf('iaintb=*',iaintb,naintb)

        call prinf('naintbc=*',naintbc,1)
        call prinf('aintbc=*',aintbc,naintbc)
        call prinf('iaintbc=*',iaintbc,naintbc)
      endif

      isuccess = 1
      do i=1,naintb
         
c
c        first check location in a array
c
        if(auni(iaintb(i)).ne.aintb(i)) then
         isuccess = -1
         print *, 'location in a array of aintb incorrect for i=',i
        endif

        imatch = 0
        do j=1,nbuni
          if(buni(j).eq.aintb(i)) imatch = 1
        enddo

        if(imatch.eq.0) then
          print *, 'aintb not found in b for i=',i 
          isuccess = 0
        endif
      enddo

      do i=1,naintbc
         
c
c        first check location in a array
c
        if(auni(iaintbc(i)).ne.aintbc(i)) then
         isuccess = 1
         print *, 'location in a array of aintbc incorrect for i=',i
        endif

        imatch = 0
        do j=1,nbuni
          if(buni(j).eq.aintbc(i)) imatch = 1
        enddo

        if(imatch.eq.1) then
          print *, 'aintbc found in b for i=',i 
          isuccess = 2
        endif
      enddo

      if(isuccess.eq.1) call prinf('test setdecomp passed*',i,0)

      stop
      end
