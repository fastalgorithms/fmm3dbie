      implicit real *8 (a-h,o-z)

      ntests = 6
      call test_get_uni(i1,i2)
      call test_setdecomp(i3)
      call test_cumsum(i4)
      call test_sort(i5)
      call test_conv_csc(i6)

      nsuccess = i1+i2+i3+i4+i5+i6

      open(unit=33,file='../../print_testres.txt')
      write(33,'(a,i1,a,i1,a)') 'Successfully completed ',nsuccess,
     1  ' out of ',ntests,' in common testing suite'
      close(33)
      

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

      ifprint = 1

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
      integer, allocatable :: aintb(:),iaintba(:),iaintbb(:),aintbc(:)
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


      allocate(aintb(nauni),iaintba(nauni),aintbc(nauni),iaintbc(nauni))
      allocate(iaintbb(nauni))
      call setdecomp(nauni,auni,nbuni,buni,naintb,aintb,iaintba,
     1     iaintbb,naintbc,aintbc,iaintbc)
      
      if(ifprint.ge.1) then
        call prinf('naintb=*',naintb,1)
        call prinf('aintb=*',aintb,naintb)
        call prinf('iaintba=*',iaintba,naintb)
        call prinf('iaintbb=*',iaintbb,naintb)

        call prinf('naintbc=*',naintbc,1)
        call prinf('aintbc=*',aintbc,naintbc)
        call prinf('iaintbc=*',iaintbc,naintbc)
      endif

      isuccess = 1
      do i=1,naintb
         
c
c        first check location in a array
c
        if(auni(iaintba(i)).ne.aintb(i)) then
         isuccess = -1
         print *, 'location in a array of aintb incorrect for i=',i
        endif

        if(buni(iaintbb(i)).ne.aintb(i)) then
         isuccess = 0
         print *, 'location in b array of aintb incorrect for i=',i
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

      return
      end


            
      subroutine test_cumsum(isuccess)
c$    use omp_lib
      implicit real *8 (a-h,o-z)
      integer, allocatable :: a(:), b(:), b2(:)
      integer :: ns(10), idiff(10), nn, nmax, i, nrange, ntimes

      call prini(6,13)

      ns(1) = 240000
      ns(2) = 270000
      ns(3) = 1000000
      ns(4) = 1500000
      ns(5) = 10000000
      nn = 5
      
      nmax = 0
      do i = 1,nn
         nmax = max(nmax,ns(i))
      enddo


      allocate(a(nmax),b(nmax),b2(nmax))

      nrange = nmax*1.2

      ifprint = 1

      isuccess = 1
      ntimes = 100
      do i = 1,nn
         do j=1,ns(i)
            a(j) = hkrand(0)*nrange
         enddo
         call cpu_time(t1)
C$       t1 = omp_get_wtime()
         do j = 1,ntimes
            call cumsum1(ns(i),a,b)
         enddo
         call cpu_time(t2)
C$       t2 = omp_get_wtime()

         call cumsum(ns(i),a,b2)
         call cpu_time(t3)
C$       t3 = omp_get_wtime()
         do j = 1,ntimes
            call cumsum(ns(i),a,b2)
         enddo
         call cpu_time(t4)
C$       t4 = omp_get_wtime()


         print *, 'speed-up ',(t2-t1)/(t4-t3), ' n = ', ns(i)

         idiff(i) = 0
         do j = 1,ns(i)
            idiff(i) = idiff(i) + abs(b(j)-b2(j))
         enddo
         if (idiff(i) .gt. 0) then
            isuccess = 0
            print *, 'para cumsum not equal to direct, ns=',ns(i)
         endif
      enddo


      if(isuccess.eq.1) call prinf('test cumsum passed*',i,0)
      
      return
      end


      
      subroutine test_sort(isuccess)
c$    use omp_lib
      implicit none
      integer isuccess
      integer, allocatable :: data(:),ind1(:),ind2(:)
      integer i,j,n,nn,ns(5)
      double precision t1,t2,t3,t4
      integer,parameter :: seed = 86456

      isuccess = 1
      ns(1) = 100001
      ns(2) = 1000003
      ns(3) = 10000005
      ns(4) = 100000006
      ns(5) = 1000000007
      nn=3
      call srand(seed)
      do i=1,nn
        n = ns(i)
        allocate(data(n),ind1(n),ind2(n))
        do j=1,n
          data(n-j+1) = irand()
        enddo
        call cpu_time(t1)
C$      t1 = omp_get_wtime()
        call sorti(n,data,ind1)
        call cpu_time(t2)
C$      t2 = omp_get_wtime()

        call cpu_time(t3)
C$      t3 = omp_get_wtime()
        call sorti_para(n,data,ind2)
        call cpu_time(t4)
C$      t4 = omp_get_wtime()
        write(*,*) "serial time: ",t2-t1,"para time: ",t4-t3
        write(*,*) "speed up: ",(t2-t1)/(t4-t3),"n=",n
        do j=2,n
          if(data(ind2(j)) < data(ind2(j-1))) then
            isuccess=0
            write(*,*) "sort error",j,data(ind2(j-1)),data(ind2(j))
          endif
        enddo
        deallocate(data,ind1,ind2)
      enddo

      if(isuccess.eq.1) call prinf('test sort passed*',i,0)

      return
      end



      subroutine test_conv_csc(isuccess)
c$    use omp_lib
      implicit none
      integer isuccess,nent,m,i,ns(5),ms(5),j,nn
      integer, allocatable :: iind(:),jind(:)
      integer, allocatable :: col_ptr1(:),row_ind1(:)
      integer, allocatable :: col_ptr2(:),row_ind2(:)
      double precision t1,t2,t3,t4
      integer,parameter :: seed = 86456

      isuccess = 1
      ns(1) = 100001
      ns(2) = 1000003
      ns(3) = 10000005
      ns(4) = 100000006
      ns(5) = 1000000007

      ms(1) = 1000
      ms(2) = 10000
      ms(3) = 100000
      ms(4) = 1000000
      ms(5) = 10000000

      nn=3
      do j=1,nn
        m = ms(j)
        nent = ns(j)
        allocate(iind(nent),jind(nent),col_ptr1(m+1),row_ind1(nent))
        allocate(col_ptr2(m+1),row_ind2(nent))

        call srand(seed)
        do i=1,nent
          iind(i) = mod(irand(),m)
          jind(i) = mod(irand(),m)
        enddo

        call cpu_time(t1)
C$      t1 = omp_get_wtime()
        call conv_to_csc(nent,m,iind,jind,col_ptr1,row_ind1)
        call cpu_time(t2)
C$      t2 = omp_get_wtime()

        call cpu_time(t3)
C$      t3 = omp_get_wtime()
        call conv_to_csc_para(nent,m,iind,jind,col_ptr2,row_ind2)
        call cpu_time(t4)
C$      t4 = omp_get_wtime()
         
        write(*,*) "serial time: ",t2-t1,"para time: ",t4-t3
        write(*,*) "speed up: ",(t2-t1)/(t4-t3),"nent=",nent,"m=",m

        do i=1,m
          if(col_ptr1(i).ne.col_ptr2(i)) isuccess = 0
        enddo

        deallocate(iind,jind,col_ptr1,col_ptr2,row_ind1,row_ind2)
      enddo

      if(isuccess.eq.1) call prinf('test convert to csc passed*',i,0)

      return
      end
