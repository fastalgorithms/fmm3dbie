
      subroutine dapplyrigidbody3d(r0,v,omega,nr,r,u)
c
c     set
c
c     u(1:3,i) = v + omega \times (r(1:3,i)-r0)
c
c     for each i (in the least squares sense)
c
c     input:
c
c     v - real *8 (3), the translation part of the
c                      formula above
c     omega - real *8 (3), the vector determining the
c                       infinitesimal rotation in the
c                       formula above
c     r0 - real *8 (3), reference point in formula above     
c     nr - integer *8, number of points
c     r - real *8 (3,nr), array of point locations
c
c     output:
c      
c     u - real *8 (3,nr), array of velocities
c
c      
      implicit none
      integer *8 nr
      real *8 u(3,nr),r(3,nr),r0(3),v(3),omega(3)
c     local
      real *8 rtemp(3)
      integer *8 i
      
      do i = 1,nr
         rtemp(1) = r(1,i)-r0(1)
         rtemp(2) = r(2,i)-r0(2) 
         rtemp(3) = r(3,i)-r0(3)

         call dcross(omega,rtemp,u(1,i))

         u(1,i) = u(1,i) + v(1)
         u(2,i) = u(2,i) + v(2)
         u(3,i) = u(3,i) + v(3)
      enddo

      return
      end

            subroutine dfindrigidbody3d(nu,u,r,r0,v,omega)
c
c     solve for v and omega so that
c
c     u(1:3,i) = v + omega \times (r(1:3,i)-r0)
c
c     for each i (in the least squares sense)
c
c     input:
c
c     nu - integer *8, number of velocities and points
c     u - real *8 (3,nu), array of velocities
c     r - real *8 (3,nu), array of point locations
c     r0 - real *8 (3), reference point in formula above
c
c     output:
c
c     v - real *8 (3), the translation part of the
c                      formula above
c     omega - real *8 (3), the vector determining the
c                       infinitesimal rotation in the
c                       formula above
c      
      implicit none
      integer *8 nu
      real *8 u(3,nu),r(3,nu),r0(3),v(3),omega(3)
c     local

      real *8, allocatable :: amat(:,:), rhs(:), sol(:)
      real *8 bmat(3,3), epsleast, rtemp(3)
      integer *8 m, n, i, j, nrhs, info, irank, istart
      
      m = nu*4
      n = 6
      nrhs = 1
      

      allocate(amat(m,n),rhs(m),sol(n))

      do j = 1,n
         do i = 1,m
            amat(i,j) = 0
         enddo
      enddo

      do i = 1,nu
         rtemp(1) = r(1,i)-r0(1)
         rtemp(2) = r(2,i)-r0(2)
         rtemp(3) = r(3,i)-r0(3)         
         rhs(i) = u(1,i)*rtemp(1) + u(2,i)*rtemp(2)
     1        + u(3,i)*rtemp(3)
      enddo

      do i = 1,3
         do j = 1,nu
            amat(j,i) = r(i,j)
         enddo
      enddo

      do i = 1,nu
         istart = 3*(i-1)+nu
         rhs(istart+1) = -u(1,i)
         rhs(istart+2) = -u(2,i)
         rhs(istart+3) = -u(3,i)

         amat(istart+1,1) = -1
         amat(istart+2,2) = -1
         amat(istart+3,3) = -1
         
         rtemp(1) = r(1,i)-r0(1)
         rtemp(2) = r(2,i)-r0(2)
         rtemp(3) = r(3,i)-r0(3)         
         
         call dcrossvectomat(rtemp,bmat)
         amat(istart+1,4) = bmat(1,1)
         amat(istart+2,4) = bmat(2,1)
         amat(istart+3,4) = bmat(3,1)         
         amat(istart+1,5) = bmat(1,2)
         amat(istart+2,5) = bmat(2,2)
         amat(istart+3,5) = bmat(3,2)         
         amat(istart+1,6) = bmat(1,3)
         amat(istart+2,6) = bmat(2,3)
         amat(istart+3,6) = bmat(3,3)         
      enddo
         
      epsleast = 1d-14
      call dleastsq(m,n,amat,nrhs,rhs,epsleast,info,sol,irank)

c      call prinf('in dfindrigidbody3d, info *',info,1)
c      call prinf('in dfindrigidbody3d, irank *',irank,1)

      do i = 1,3
         v(i) = sol(i)
         omega(i) = sol(i+3)
      enddo


      return
      end

      
      subroutine dcross(a,b,c)
c
c     c = a \times b
c
      implicit none
      real *8 a(3),b(3),c(3),a1,a2,a3,b1,b2,b3

      a1 = a(1)
      a2 = a(2)
      a3 = a(3)
      b1 = b(1)
      b2 = b(2)
      b3 = b(3)
      
      c(1) = a2*b3-a3*b2
      c(2) = a3*b1-a1*b3
      c(3) = a1*b2-a2*b1
      
      return
      end
      
      subroutine dcrossvectomat(r,a)
c
c     return the matrix corresponding to
c     the cross product for the given vector,
c     i.e. return a such that
c
c     a x = r \times x
c     
      implicit none
      real *8 r(3), a(3,3), r1, r2, r3

      r1 = r(1)
      r2 = r(2)
      r3 = r(3)

      a(1,1) = 0
      a(2,1) = r3
      a(3,1) = -r2
      a(1,2) = -r3
      a(2,2) = 0
      a(3,2) = r1
      a(1,3) = r2
      a(2,3) = -r1
      a(3,3) = 0

      return
      end

     
      
