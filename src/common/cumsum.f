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

      integer i,isum

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
