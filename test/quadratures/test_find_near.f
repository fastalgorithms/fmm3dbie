      subroutine test_find_near(nsuccess)

      implicit real *8 (a-h,o-z)
      real *8, allocatable :: src(:,:),targ(:,:),rads(:)
      integer, allocatable :: row_ptr(:),row_ptr2(:)
      integer, allocatable :: col_ind(:),col_ind2(:)
      integer, allocatable :: isort(:),wsort(:),isort2(:),wsort2(:)

      call prini(6,13)


      ns = 101
      nt = 1000
      i0 = 1
      i1 = 1
      i2 = 1
      
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
      call findnearslowmem(src,ns,rads,3,targ,nt,nnz2)
      call findnearmem(src,ns,rads,3,targ,nt,nnz)

      if(nnz.ne.nnz2) then
        call prinf('number of non zero elements dont match*',i,0)
        i0 = 0
      endif

      allocate(col_ind(nnz),col_ind2(nnz))
      call findnear(src,ns,rads,3,targ,nt,row_ptr,col_ind)
      call findnearslow(src,ns,rads,3,targ,nt,row_ptr2,col_ind2)

      do i=1,nt
        n1 = row_ptr(i+1)-row_ptr(i)
        n2 = row_ptr2(i+1)-row_ptr2(i)
        if(n1.ne.n2) then
          call prinf('number of sources doesnt match for target i=*',
     1       i,1)
          i1 = 0
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
          i2 = 0
        endif
      enddo

      isuc = 0
      if(i0+i1+i2.eq.3) isuc = 1
      nsuccess = isuc


      return
      end
