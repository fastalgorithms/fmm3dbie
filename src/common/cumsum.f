c
c
c        compute the cumulative sum of an integer array
c
      subroutine cumsum(n,a,b)
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
c            b(i) = sum_{j=1}^{i} a(j)
c            
c
      implicit none
      integer a(n),b(n),n,i,isum

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
