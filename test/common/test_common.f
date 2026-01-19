      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)

      ntests = 8
      call test_get_uni(i1,i2)
      call test_setdecomp(i3)
      call test_cumsum(i4)
      call test_sort(i5)
      call test_conv_csc(i6)
      call test_findnear(i7)
      call test_rsc_to_csc(i8)

      nsuccess = i1+i2+i3+i4+i5+i6+i7+i8

      open(unit=33,file='../../print_testres.txt')
      write(33,'(a,i1,a,i1,a)') 'Successfully completed ',nsuccess,
     1  ' out of ',ntests,' in common testing suite'
      close(33)
      

      stop
      end



      subroutine test_get_uni(isuccess0,isuccess1)

      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, allocatable :: a(:),iuni(:),iuniind(:)
      integer *8, allocatable :: b(:),c(:)
      integer *8, allocatable :: iuni3(:,:)
      integer *8 n, nuni


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
      implicit integer *8 (i-n)
      integer *8, allocatable :: a(:),iauni(:),auni(:)
      integer *8, allocatable :: b(:),ibuni(:),buni(:)
      integer *8, allocatable :: aintb(:),iaintba(:)
      integer *8, allocatable :: iaintbb(:),aintbc(:)
      integer *8, allocatable :: iaintbc(:),binta(:),ibinta(:)
      integer *8, allocatable :: bintac(:),ibaintac(:)
      integer *8, allocatable :: itmp1(:),itmp2(:)
      integer *8, allocatable :: aunisort(:),bunisort(:)
      integer *8 nauni, nbuni, m, n


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
      implicit integer *8 (i-n)
      integer *8, allocatable :: a(:), b(:), b2(:)
      integer *8 :: ns(10), idiff(10), nn, nmax, i, nrange, ntimes

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
      integer *8 isuccess
      integer *8, allocatable :: data(:),ind1(:),ind2(:)
      integer *8 i,j,n,nn,ns(5)
      double precision t1,t2,t3,t4,hkrand
      integer *8,parameter :: seed = 86456
      integer *8 :: irand,nmax

      isuccess = 1
      ns(1) = 100001
      ns(2) = 1000003
      ns(3) = 10000005
      ns(4) = 100000006
      ns(5) = 1000000007
      nn=2
      nmax = 100000000
      do i=1,nn
        n = ns(i)
        allocate(data(n),ind1(n),ind2(n))
        do j=1,n
          data(n-j+1) = hkrand(0)*nmax 
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
      integer *8 isuccess,nent,m,i,ns(5),ms(5),j,nn
      integer *8, allocatable :: iind(:),jind(:)
      integer *8, allocatable :: col_ptr1(:),row_ind1(:)
      integer *8, allocatable :: col_ptr2(:),row_ind2(:)
      double precision t1,t2,t3,t4,hkrand
      integer *8,parameter :: seed = 86456
      integer *8 :: ii,jj

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

      nn=2
      do j=1,nn
        m = ms(j)
        nent = ns(j)
        allocate(iind(nent),jind(nent),col_ptr1(m+1),row_ind1(nent))
        allocate(col_ptr2(m+1),row_ind2(nent))

c        call srand(seed)
        do i=1,nent
          ii = hkrand(0)*m*2
          jj = hkrand(0)*m*2
          iind(i) = mod(ii,m)
          jind(i) = mod(jj,m)
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
c
c
c
c
c
      subroutine test_findnear(isuccess)
      implicit none
      real *8, allocatable :: src(:,:),targ(:,:),rads(:)
      real *8 hkrand
      integer *8, allocatable :: row_ptr(:),row_ptr2(:)
      integer *8, allocatable :: col_ind(:),col_ind2(:)
      integer *8, allocatable :: isort(:),wsort(:),isort2(:),wsort2(:)
      integer *8 n1, n2, isuccess, ns, nt, erra
      integer *8 i, j, nnz, nnz2, ndt

      call prini(6,13)


      ns = 1001
      nt = 10000
      isuccess = 1
      
      allocate(src(3,ns),targ(3,nt),rads(ns))
      allocate(row_ptr(nt+1),row_ptr2(nt+1))
      allocate(isort(nt),wsort(nt),isort2(nt),wsort2(nt))

      do i=1,ns
        src(1,i) = hkrand(0)
        src(2,i) = hkrand(0)
        src(3,i) = hkrand(0)
        rads(i) = 2.0d0**(-10*hkrand(0))/3.0d0
      enddo

      do i=1,nt
        targ(1,i) = hkrand(0)
        targ(2,i) = hkrand(0)
        targ(3,i) = hkrand(0)
      enddo

      nnz = 0
      nnz2 = 0
      ndt = 3
      call findnearslowmem(src,ns,rads,ndt,targ,nt,nnz2)
      call findnearmem(src,ns,rads,ndt,targ,nt,nnz)


      if(nnz.ne.nnz2) then
        call prinf('number of non zero elements dont match*',i,0)
        isuccess = 0
        return
      endif

      allocate(col_ind(nnz),col_ind2(nnz))
      call findnear(src,ns,rads,ndt,targ,nt,row_ptr,col_ind)
      call findnearslow(src,ns,rads,ndt,targ,nt,row_ptr2,col_ind2)

      do i=1,nt
        n1 = row_ptr(i+1)-row_ptr(i)
        n2 = row_ptr2(i+1)-row_ptr2(i)
        if(n1.ne.n2) then
          call prinf('number of sources doesnt match for target i=*',
     1       i,1)
          isuccess = 0
          return
        endif

        call sorti(n1,col_ind(row_ptr(i)),wsort)
        call sorti(n2,col_ind2(row_ptr2(i)),wsort2)

        erra = 0
        do j=1,n1
          isort(j) = col_ind(row_ptr(i)+wsort(j)-1)
          isort2(j) = col_ind2(row_ptr2(i)+wsort2(j)-1)
          erra = erra + abs(isort(j)-isort2(j))
        enddo

        if(erra.ne.0) then
          call prinf('list of sources dont match for target i=*',i,1)
          call prinf('correct sources=*',isort2,n2)
          call prinf('computed sources=*',isort,n1)
          isuccess = 0
          return
        endif
      enddo
      
      print *, "find near test successfully completed"



      return
      end



      subroutine test_rsc_to_csc(isuccess)
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8, allocatable :: row_ptr(:),col_ind(:),col_ptr(:),
     1   row_ind(:),iper(:)
      integer *8 ns, nt, nnz

      
      ns = 301
      nt = 101

      isuccess = 1

      allocate(row_ptr(nt+1),col_ptr(ns+1))
      

      nnz_max = ns*nt
      allocate(col_ind(nnz_max),iper(nnz_max),row_ind(nnz_max))

      
      row_ptr(1) = 1
      do i=1,nt
        nrel = hkrand(0)*0.5*ns
        do j=1,nrel
          col_ind(row_ptr(i)+j-1) = ceiling(hkrand(0)*ns)
        enddo
        row_ptr(i+1) = row_ptr(i) + nrel
      enddo
      nnz = row_ptr(nt+1)-1

      call rsc_to_csc(ns,nt,nnz,row_ptr,col_ind,col_ptr,row_ind,iper)

      do i=1,ns
        do j=col_ptr(i),col_ptr(i+1)-1
           itarg = row_ind(j)
           ii1 = iper(j)
           isrc = col_ind(ii1)
           if(ii1.lt.row_ptr(itarg).or.ii1.ge.row_ptr(itarg+1)) then
              call prinf('error in permutation array*',i,0)
              print *, ii1,row_ptr(itarg),row_ptr(itarg+1)
              isuccess = 0
              return
           endif

           if(isrc.ne.i) then
             call prinf('error in converting rsc to csc*',i,0)
             print *, "No match found for interaction between ",isrc,
     1          itarg
             isuccess = 0
             return

           endif
        enddo
      enddo

      print *, "rsc_to_csc completed successfully"

      return
      end

