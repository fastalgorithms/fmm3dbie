c-------------------------------------------
c> @brief
c> Compute the cumulative sum of a vector:
c>
c> \f{equation*}{b_{i} = \sum_{j=1}^{i} a_{j} \f}
c>
c> @param[in] n: length of array
c> @param[in] a: input array
c> @param[out] b: cumulative sum of input array
c------------------------------------------------
      subroutine cumsum(n,a,b)
      implicit none
      integer, intent(in) :: n,a(n)
      integer, intent(out) :: b(n)
      
      integer lend, nsmall, maxth, offset
      parameter (lend = 200)
      integer d(lend)
      integer, allocatable :: d2(:)
c$    integer omp_get_max_threads
      data nsmall / 5000 /

c     if small problem, no parallel
      if (n .lt. nsmall) goto 300

c     if not a ton of processors, use d on stack
      maxth = 1
c$    maxth = omp_get_max_threads()
      if (maxth .le. lend) goto 200

c     if tons of processors, allocate d2
      allocate(d2(maxth))
      call cumsum_para2(n,a,b,maxth,d2)
      return

 200  continue

      call cumsum_para2(n,a,b,lend,d)
      return
      
 300  continue
            
      call cumsum1(n,a,b)
      return
      
      end
c
c
      subroutine cumsum1(n,a,b)
      implicit none
      integer, intent(in) :: n,a(n)
      integer, intent(out) :: b(n)

      integer i,isum

      isum = 0

      do i=1,n
        isum = isum + a(i)
        b(i) = isum
      enddo
      
      return
      end
c

c      subroutine cumsum_para(n,a,b,nd,d)
cc--------------------------------------------
cc      
c      implicit none
c      integer, intent(in) :: n,a(n),nd
c      integer, intent(out) :: b(n),d(nd)
c
c      integer i,isum,offset,istart,iend
c      integer nth, id, nbatch, inext
cc$    integer omp_get_thread_num, omp_get_num_threads      
c
c
cc$OMP PARALLEL DEFAULT(none) SHARED(a,b,n,nth,d,nbatch)
cc$OMP$ PRIVATE(i,isum,offset,istart,iend,id,inext)
c
c
cc     get some info on parallel set up (serial)
c      
cc$OMP SINGLE      
c      nth = 1
cc$    nth = omp_get_num_threads()
c      nbatch = (n-1)/nth+1
cc$OMP END SINGLE      
c
c      
cc     partial sums (parallel)
c
c      id = 1
cc$    id = omp_get_thread_num()+1
c
c      istart = 1+(id-1)*nbatch
c      iend = min(istart + nbatch-1,n)
c
c      isum = 0
c      do i = istart,iend
c         isum = isum + a(i)
c         b(i) = isum
c      enddo
c      d(id) = isum
c
cc$OMP BARRIER
c      
cc     accumulate over ends of partial sums (serial)
c
cc$OMP SINGLE      
c      isum = 0
c      do i = 1,nth
c         inext = d(i)
c         d(i) = isum
c         isum = isum+inext
c      enddo
cc$OMP END SINGLE
c
cc     broadcast these accumulations (parallel)
c
c      offset = d(id)
c      do i = istart,iend
c         b(i) = b(i) + offset
c      enddo
c      
cc$OMP END PARALLEL
c      
c      return
c      end
cc


      
      subroutine cumsum_para2(n,a,b,nd,d)
c----------------------------------------------------------------------
c
c
c      
c     WARNING: this subroutine depends on the properties
c     of the "static" schedule in the OPENMP
c     specification. Specifically, the static schedule
c     with no chunk_size specified assigns at most one
c     block of iterates to each thread in thread order.
c     The assignment of indices is the same for any do
c     loops of the same cardinality within a parallel
c     region. Changing the schedule will make the code non
c     conforming.
c
c----------------------------------------------------------------------      
      
      implicit none
      integer, intent(in) :: n,a(n),nd
      integer, intent(out) :: b(n),d(nd)

      integer i,isum,offset
      integer nth, id, inext
c$    integer omp_get_thread_num, omp_get_num_threads      


c$OMP PARALLEL DEFAULT(none) SHARED(a,b,n,d)
c$OMP$ PRIVATE(i,isum,offset,id,inext,nth)

      id = 1
c$    id = omp_get_thread_num()+1

c     compute cumulative sums of portions (parallel)
      
      
      isum = 0
c$OMP DO SCHEDULE(static)       
      do i = 1,n
         isum = isum+a(i)
         b(i) = isum
      enddo
c$OMP END DO no wait      
      d(id) = isum

c$OMP BARRIER

c     accumulate the ends of the partial sums (each thread)
      
      offset = 0
      do i = 1,id-1
         offset = offset + d(i)
      enddo

c$OMP DO SCHEDULE(static)       
      do i = 1,n
         b(i) = b(i) + offset
      enddo
c$OMP END DO no wait      
      
c$OMP END PARALLEL
      
      return
      end
c
      
