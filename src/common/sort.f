c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sortr(n,data,index)
C===================================================================
C
C This routine performs an sort of the first n elements
C of array data, returning into array index the indices of elements of
C data arranged in ascending order.  Thus,
C
C    data(index(1)) will be the smallest number in array data;
C    data(index(n)) will be the largest number in data.
C
C The original data is not physically rearranged.  The original order
C of equal input values is not necessarily preserved.
C
C===================================================================
      implicit none
      integer i,n,index(n)
      real *8 data(n)

      do i=1,n
        index(i) = i
      enddo

      call sortr_help(n,n,data,index)
      end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sortr_help(m,n,data,index)
C===================================================================
C
C
C===================================================================
      implicit none
      integer n,m,n_thresh,p,i,j,pivot_index,index(n)
      real *8 pivot,data(m)
      parameter (n_thresh = 16)
      if (n<=n_thresh) then
          call insertion_sortr(m,n,data,index)
      else
c         choose the pivot
          p = (1+n)/2
          pivot_index = index(p)
          pivot = data(pivot_index)
          if (data(index(1)) > pivot) then
              index(p) = index(1)
              index(1) = pivot_index
              pivot_index = index(p)
              pivot = data(pivot_index)
          endif
          if (pivot > data(index(n))) then
              if (data(index(1)) > data(index(n))) then
                  index(p) = index(1)
                  index(1) = index(n)
              else
                  index(p) = index(n)
              endif
              index(n) = pivot_index
              pivot_index = index(p)
              pivot = data(pivot_index)
          endif
c         end of choose the pivot
c         start to do swaps and partition the array using pivot
          i = 1
          j = n
          do while (i<j)
            i = i+1
            do while (data(index(i)) < pivot)
              i = i+1
            enddo
            j = j-1
            do while (data(index(j)) > pivot)
              j = j-1
            enddo
            if (i<j) then
                p = index(i)
                index(i) = index(j)
                index(j) = p
            endif
          enddo
c         done partition
c         sort two subarrays
          if (i-1 > 1) then
              call sortr_help(m, i-1, data, index(1))
          endif
          if (n-j > 1) then
              call sortr_help(m, n-j, data, index(j+1))
          endif
      endif
      end
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine insertion_sortr(m,n,data,index)
C===================================================================
C
C This routine performs an sort of the first N elements
C of array DATA, returning into array INDEX the indices of elements of
C DATA arranged in ascending order.  Thus,
C
C    DATA(INDEX(1)) will be the smallest number in array DATA;
C    DATA(INDEX(N)) will be the largest number in DATA.
C
C The original data is not physically rearranged.  The original order
C of equal input values is not necessarily preserved.
C This subroutine uses insertion sort.
C
C===================================================================
      implicit none
      integer n,m,i,j,key_index,index(n)
      real *8 key,data(m)
      do i = 2,n
        if(data(index(i-1)) > data(index(i))) then
          key_index = index(i)
          key = data(key_index)
          j = i - 1
          do while (data(index(j))>key .and. j>0)
            index(j+1) = index(j)
            j = j - 1
          enddo
          index(j+1) = key_index
        endif
      enddo
      end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sorti(n,data,index)
C===================================================================
C
C This routine performs an sort of the first n elements
C of array data, returning into array index the indices of elements of
C data arranged in ascending order.  Thus,
C
C    data(index(1)) will be the smallest number in array data;
C    data(index(n)) will be the largest number in data.
C
C The original data is not physically rearranged.  The original order
C of equal input values is not necessarily preserved.
C
C===================================================================
      implicit none
      integer i,n,index(n)
      integer data(n)

      do i=1,n
        index(i) = i
      enddo

      call sorti_help(n,n,data,index)
      end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sorti_help(m,n,data,index)
C===================================================================
C
C
C===================================================================
      implicit none
      integer n,m,n_thresh,p,i,j,pivot_index,index(n)
      integer pivot,data(m)
      parameter (n_thresh = 16)
      if (n<=n_thresh) then
          call insertion_sorti(m,n,data,index)
      else
c         choose the pivot
          p = (1+n)/2
          pivot_index = index(p)
          pivot = data(pivot_index)
          if (data(index(1)) > pivot) then
              index(p) = index(1)
              index(1) = pivot_index
              pivot_index = index(p)
              pivot = data(pivot_index)
          endif
          if (pivot > data(index(n))) then
              if (data(index(1)) > data(index(n))) then
                  index(p) = index(1)
                  index(1) = index(n)
              else
                  index(p) = index(n)
              endif
              index(n) = pivot_index
              pivot_index = index(p)
              pivot = data(pivot_index)
          endif
c         end of choose the pivot
c         start to do swaps and partition the array using pivot
          i = 1
          j = n
          do while (i<j)
            i = i+1
            do while (data(index(i)) < pivot)
              i = i+1
            enddo
            j = j-1
            do while (data(index(j)) > pivot)
              j = j-1
            enddo
            if (i<j) then
                p = index(i)
                index(i) = index(j)
                index(j) = p
            endif
          enddo
c         done partition
c         sort two subarrays
          if (i-1 > 1) then
              call sorti_help(m, i-1, data, index(1))
          endif
          if (n-j > 1) then
              call sorti_help(m, n-j, data, index(j+1))
          endif
      endif
      end
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine insertion_sorti(m,n,data,index)
C===================================================================
C
C This routine performs an sort of the first N elements
C of array DATA, returning into array INDEX the indices of elements of
C DATA arranged in ascending order.  Thus,
C
C    DATA(INDEX(1)) will be the smallest number in array DATA;
C    DATA(INDEX(N)) will be the largest number in DATA.
C
C The original data is not physically rearranged.  The original order
C of equal input values is not necessarily preserved.
C This subroutine uses insertion sort.
C
C===================================================================
      implicit none
      integer n,m,i,j,key_index,index(n)
      integer key,data(m)
      do i = 2,n
        if(data(index(i-1)) > data(index(i))) then
          key_index = index(i)
          key = data(key_index)
          j = i - 1
          do while (data(index(j))>key .and. j>0)
            index(j+1) = index(j)
            j = j - 1
          enddo
          index(j+1) = key_index
        endif
      enddo
      end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sorti_para(n,data,index)
C===================================================================
C
C     sorti_para -- sort, integer input, index output
C
C
C     Input:  n     integer
C             data  integer
C
C     Output: index integer (dimension n)
C
C This openmp threaded routine performs an sort of the first n elements
C of array data, returning into array index the indices of elements of
C data arranged in ascending order.  Thus,
C
C    data(index(1)) will be the smallest number in array data;
C    data(index(n)) will be the largest number in data.
C
C The original data is not physically rearranged.  The original order
C of equal input values is not necessarily preserved.
C===================================================================

      implicit none
      integer n,data(n),index(n)
      integer nthreads,nelems,ielems,istart,iend,i,j,swap
      integer(8) ithread
      integer istart2,iend2,nelems2
      integer, allocatable :: split(:)
      integer, allocatable :: index2(:)
      integer, external :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
      integer, external :: OMP_GET_MAX_THREADS

      nthreads=1
C$    nthreads=OMP_GET_MAX_THREADS()
c      write(*,*) "max num threads: ", nthreads

      if(N<2*nthreads .or. nthreads==1) then
c        write(*,*) "call sorti"
        call sorti(n,data,index)
        return
      endif

c      write(*,*) "call omp sort"
ccc   split array evenly
      allocate(split(nthreads+1))
      allocate(index2(n))
      split(nthreads+1)=n+1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ithread)
      do ithread = 1,nthreads
        split(ithread) = ((ithread-1)*n)/nthreads+1
      enddo
C$OMP END PARALLEL DO
c      write(*,*) nthreads+1,split(1),split(2),split(3)

ccc   sort each part
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ithread,nelems,istart,iend)
      do ithread = 1,nthreads
        istart = split(ithread)
        iend = split(ithread+1)
        nelems = iend-istart
c        write(*,*) "nelems: ", nelems
        call sorti(nelems,data(istart),index(istart))
        do ielems = istart,iend-1
          index(ielems) = index(ielems) + istart - 1
        enddo
      enddo
C$OMP END PARALLEL DO

ccc   merge two parts at a time
      swap=1
      j=1
      do while (j<nthreads)
        i=1
        do while (i<nthreads+1)
          if(i+j<nthreads+1) then
            istart=split(i)
            iend=split(i+j)
            istart2=split(i+j)
            if(i+2*j>nthreads+1) then
              iend2=split(nthreads+1)
            else
              iend2=split(i+2*j)
            endif
            nelems=iend-istart
            nelems2=iend2-istart2
            if(swap==1) then
              call merge_arr_para(data,index(istart),index(istart2),
     1             nelems,nelems2,index2(istart),nthreads)
            else
              call merge_arr_para(data,index2(istart),index2(istart2),
     1             nelems,nelems2,index(istart),nthreads)
            endif
          else
            if(swap==0) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(istart)
              do istart=split(i),split(nthreads+1)-1
                index(istart) = index2(istart)
              enddo
C$OMP END PARALLEL DO
            else
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(istart)
              do istart=split(i),split(nthreads+1)-1
                index2(istart) = index(istart)
              enddo
C$OMP END PARALLEL DO
            endif
          endif
          i=i+2*j
        end do
        j=j*2
        if(swap==1) then
          swap=0
        else
          swap=1
        endif
      end do

      if(swap==0) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(istart)
        do istart=1,n
          index(istart) = index2(istart)
        enddo
C$OMP END PARALLEL DO
      endif

      deallocate(split)
      deallocate(index2)
      end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine merge_arr_para(data,ind1,ind2,n1,n2,indout,nthreads)
C===================================================================
C
C This openmp threaded routine performs merge of two arrays.
C ind1 contains the index of array1 of data, ind2 contains the index
C of array2 of data.
C The merged index set is in indout of increaing order, it contains
C the index respect to data, the original data is not physically
C rearranged.
C
C===================================================================
      implicit none
      integer n1,n2,nthreads,i,j,k,ns,indx,first,val,split1,split2
      integer(8) tmpk
      integer split_size1,split_size2,n1tmp,n2tmp
      integer istart1,istart2,istart3
      integer data(*)
      integer ind1(n1),ind2(n2),indout(n1+n2)
      integer, allocatable :: split(:), split_size(:)
      integer, allocatable :: split_ind1(:), split_ind2(:),indtmp(:)
      if(n1==0 .and. n2==0) return
      if(n1==0) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,k,tmpk)
        do k=0,nthreads-1
          tmpk = k
          i = (tmpk*n2)/nthreads+1
          j = ((tmpk+1)*n2)/nthreads
          indout(i:j) = ind2(i:j)
        enddo
C$OMP END PARALLEL DO
        return
      endif

      if(n2==0) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,k,tmpk)
        do k=0,nthreads-1
          tmpk = k
          i = (tmpk*n1)/nthreads+1
          j = ((tmpk+1)*n1)/nthreads
          indout(i:j) = ind1(i:j)
        enddo
C$OMP END PARALLEL DO
        return
      endif

      ns=10
      allocate(split(nthreads*ns*2))
      allocate(split_size(nthreads*ns*2))
      allocate(indtmp(nthreads*ns*2))
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,k,indx,first,val,tmpk)
      do i=0,nthreads-1
        do j=0,ns-1
          k = i*ns+j
          tmpk = k
          indtmp(k+1)=k+1
          indx = (tmpk*n1)/(nthreads*ns)
          split(k+1) = indx+1
          val = data(ind1(split(k+1)))
          call lower_bound(data,ind2,n2,val,first)
          split_size(k+1) = indx+first-1

          indx = (tmpk*n2)/(nthreads*ns)
          k = k+nthreads*ns
          split(k+1) = indx+1
          val = data(ind2(split(k+1)))
          call lower_bound(data,ind1,n1,val,first)
          split_size(k+1) = indx+first-1
        enddo
      enddo
C$OMP END PARALLEL DO

      allocate(split_ind1(nthreads+1),split_ind2(nthreads+1))
      split_ind1(1)=1
      split_ind2(1)=1
      split_ind1(nthreads+1)=n1+1
      split_ind2(nthreads+1)=n2+1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,k,split1,split2,split_size1,split_size2,val,first)
C$OMP$PRIVATE(tmpk)
      do i=1,nthreads-1
        tmpk = i
        k = (tmpk*(n1+n2))/nthreads

        call lower_bound(split_size(1),indtmp,nthreads*ns,
     1       k,first)
        j=first
        if(j>nthreads*ns) j=nthreads*ns
        split1=split(j)
        split_size1=split_size(j)

        call lower_bound(split_size(1+nthreads*ns),indtmp,nthreads*ns,
     1       k,first)
        j=first+nthreads*ns
        if(j>2*nthreads*ns) j=2*nthreads*ns
        if(abs(k-split_size(j)) < abs(k-split_size1)) then
          val = data(ind2(split(j)))
        else
          val = data(ind1(split1))
        endif
        call lower_bound(data,ind1,n1,val,first)
        split_ind1(i+1) = first
        call lower_bound(data,ind2,n2,val,first)
        split_ind2(i+1) = first
      enddo
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,n1tmp,n2tmp,istart1,istart2,istart3)
      do i=1,nthreads
        n1tmp = split_ind1(i+1) - split_ind1(i)
        n2tmp = split_ind2(i+1) - split_ind2(i)
        istart1 = split_ind1(i)
        istart2 = split_ind2(i)
        istart3 = split_ind1(i)+split_ind2(i)-1
c        write(*,*) "i:",i,"n1:",n1tmp,"n2:",n2tmp
        call merge_arr(data,ind1(istart1),ind2(istart2),n1tmp,n2tmp,
     1       indout(istart3),nthreads)
      enddo
C$OMP END PARALLEL DO

      return
      end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine merge_arr(data,ind1,ind2,n1,n2,indout,nthreads)
      implicit none
      integer n1,n2,nthreads,i,j,k
      integer data(*)
      integer ind1(n1),ind2(n2),indout(n1+n2)

      i=1
      j=1
      k=1

      do while(i<=n1.and.j<=n2)
        if(data(ind1(i))<data(ind2(j))) then
          indout(k)=ind1(i)
          i=i+1
        else
          indout(k)=ind2(j)
          j=j+1
        endif
        k=k+1
      enddo

      do while(i<=n1)
        indout(k) = ind1(i)
        i=i+1
        k=k+1
      enddo

      do while(j<=n2)
        indout(k) = ind2(j)
        j=j+1
        k=k+1
      enddo

      return
      end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine lower_bound(data,ind,n,val,first)
        integer data(*),ind(*),n
        integer val
        integer it,first
        integer cnt,step

        first=1
        cnt=n
        do while(cnt>0)
          it = first
          step = cnt/2
          it = it + step
          if(data(ind(it))<val) then
            first = it+1
            cnt = cnt-step-1
          else
            cnt = step
          endif
        enddo

      return
      end
