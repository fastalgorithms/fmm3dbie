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
      stop
      end
