      subroutine dmatmatt(m,n,a,k,b,c)
      implicit real *8 (a-h,o-z)
      real *8 a(m,n), b(k,n), c(k,m)

c
c        this routine computes 
c         c = (ab^{t})^{t}
c

      do i=1,m
        do j=1,k
          c(j,i) = 0
        enddo
      enddo

      do j=1,n
        do i=1,m
          do l=1,k
            c(l,i) = c(l,i) + a(i,j)*b(l,j)
          enddo
        enddo
      enddo
      

      return
      end

      

      subroutine zrmatmatt_slow(m,n,k,a,b,c)
      implicit real *8 (a-h,o-z)
      complex *16 a(n,m),c(k,m)
      real *8 b(n,k)

      do i=1,m
        do j=1,k
          c(j,i) = 0
          do l=1,n
            c(j,i) = c(j,i) + a(l,i)*b(l,j)
          enddo
        enddo
      enddo

      return
      end

