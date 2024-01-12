c
c
c
c     this file contains routines for working with tensor
c     product gauss-legendre nodes and polynomials in 2 and
c     3 dimensions; also provides some related spectral
c     differentiation operators
c
c--------------------------------------------------------------
c
c     CONVENTIONS
c
c     ndeg refers to polynomial degree
c      
c     n refers to order of approximation and number of nodes
c     used in each direction (n = ndeg+1)
c
c     type character specifies the set of tensor
c     polynomials to use. current options
c
c     type = 'F', full degree polynomials (T_i(x)T_j(y)T_k(z) where
c                 each of i, j, and k goes from 0 to ndeg)
c     type = 'T', total degree polynomials (T_i(x)T_j(y)T_k(z) where
c                 sum of i,j,k>=0 is less than or equal to ndeg)
c
c     TODO: implement "Euclidean degree"
c
c     a tensor grid of points is traversed with x on the inner
c        loop, y on the loop above that, and z above that.
c        e.g. for a 2x2x2 grid, we have the order:
c            (x1,y1,z1), (x2,y1,z1), (x1,y2,z1), (x2,y2,z1)
c     (x1,y1,z2), (x2,y1,z2), (x1,y2,z2), (x2,y2,z2)
c      
c     tensor product polynomials are numbered analagously
c        e.g. for full degree at ndeg = 3, we would have the order
c        1, T_1(x), T_2(x), T_1(y), T_1(x)T_1(y), T_2(x)T_1(y),
c        T_2(y), T_1(x)T_2(y), T_2(x)T_2(y), T_1(z), T_1(x)T_1(z),
c        T_2(x) T_1(z), T_1(y)T_1(z), T_1(x)T_1(y)T_1(z), T_2(x)
c        T_1(y)T_1(z), T_2(y) T_1(z), T_1(x)T_2(y) T_1(z),
c        T_2(x) T_2(y) T_1(z), T_2(z), T_1(x)T_2(z), T_2(x) T_2(z),
c        T_1(y)T_2(z), T_1(x)T_1(y)T_2(z), T_2(x) T_1(y) T_2(z),
c        T_2(y) T_2(z), T_1(x)T_2(y) T_2(z), T_2(x) T_2(y) T_2(z)
c 
c        e.g. for total degree, we would have the order
c        1, T_1(x), T_2(x), T_1(y), T_1(x)T_1(y), T_2(y),
c        T_1(z), T_1(x)T_1(z), T_1(y)T_1(z), T_2(z)
c
c-------------------------------------------------------------     
c
c     ROUTINES
c      
c     legetens_npol_*d - total number of polynomials up
c                        to a given degree of specified
c                        type in 2 or 3 dimensions
c     legetens_exps_*d - get nodes, weights, and, optionally,
c                        projection and evaluation matrices      
c                        in 2 or 3 dimensions
c     legecoeff_dmat - rectangular 1d spectral diff mat
c     legecoeff_d2mat - rectangular 1d spectral 2nd deriv mat
c     legetens_pols_*d - evaluate the tensor legendre polys at
c                        a given point      
c     legetens_lapmat_*d - rectangular spectral laplacian in
c                          2 and 3 dims
c     legetens_eyemat_*d - rectangular "identity" matrix between
c                          different polynomial orders in
c                          2 and 3 dims 

      subroutine polytens_npol_2d(ndeg,type,npol)
c
c     return the number of polynomials of a given type
c     up to degree ndeg
c      
      integer *8 n, npol, ndeg
      character type

      n = ndeg + 1
      
      if (type .eq. 'f' .or. type .eq. 'F') then
         npol = n**2
      else if (type .eq. 't' .or. type .eq. 'T') then
         npol = n*(n+1)/2
      endif

      return
      end

c
c
c
c
      subroutine polytens_npol_3d(ndeg,type,npol)
c
c     return the number of polynomials of a given type
c     up to degree ndeg
c      
      integer *8 n, npol, ndeg
      character type

      n = ndeg + 1
      
      if (type .eq. 'f' .or. type .eq. 'F') then
         npol = n**3
      else if (type .eq. 't' .or. type .eq. 'T') then
         npol = n*(n+1)*(n+2)/6
      endif

      return
      end
c
c
c
c


      subroutine polytens_pow2ind_3d(ndeg,type,ip2ind)
      integer *8 ndeg, ip2ind(ndeg+1,ndeg+1,ndeg+1)
      character type

      integer *8 i, j, k, ipol

      do i = 1,ndeg+1
         do j = 1,ndeg+1
            do k = 1,ndeg+1
               ip2ind(k,j,i) = -1
            enddo
         enddo
      enddo

      if (type .eq. 'f' .or. type .eq. 'F') then
         ipol = 0
         do i = 1,ndeg+1
            do j = 1,ndeg+1
               do k = 1,ndeg+1
                  ipol = ipol+1
                  ip2ind(k,j,i) = ipol
               enddo
            enddo
         enddo
      else if (type .eq. 't' .or. type .eq. 'T') then
         ipol = 0
         do i = 1,ndeg+1
            do j = 1,ndeg+1+1-i
               do k = 1,ndeg+1+2-i-j
                  ipol = ipol+1
                  ip2ind(k,j,i) = ipol
               enddo
            enddo
         enddo
      endif

      return
      end
c
c
c
c
      subroutine polytens_ind2pow(ndim,ndeg,type,iind2p)
      integer *8 ndeg, iind2p(ndim,*)
      character type

      integer *8 i

      if (ndim.eq.1) then
         do i=1,ndeg+1
            iind2p(1,i)=i-1
         enddo
      elseif (ndim.eq.2) then
         call polytens_ind2pow_2d(ndeg,type,iind2p)
      elseif (ndim.eq.3) then
         call polytens_ind2pow_3d(ndeg,type,iind2p)
      endif
      
      return
      end
c
c
c
c
      subroutine polytens_ind2pow_3d(ndeg,type,iind2p)
      integer *8 ndeg, iind2p(3,*)
      character type

      integer *8 i, j, k, ipol


      if (type .eq. 'f' .or. type .eq. 'F') then
         ipol = 0
         do i = 1,ndeg+1
            do j = 1,ndeg+1
               do k = 1,ndeg+1
                  ipol = ipol+1
                  iind2p(1,ipol) = k-1
                  iind2p(2,ipol) = j-1
                  iind2p(3,ipol) = i-1                
               enddo
            enddo
         enddo
      else if (type .eq. 't' .or. type .eq. 'T') then
         ipol = 0
         do i = 1,ndeg+1
            do j = 1,ndeg+1+1-i
               do k = 1,ndeg+1+2-i-j
                  ipol = ipol+1
                  iind2p(1,ipol) = k-1
                  iind2p(2,ipol) = j-1
                  iind2p(3,ipol) = i-1
               enddo
            enddo
         enddo
      endif

      return
      end
c
c
c
c
      subroutine polytens_pow2ind_2d(ndeg,type,ip2ind)
      integer *8 ndeg, ip2ind(ndeg+1,ndeg+1)
      character type

      integer *8 i, j, ipol

      do i = 1,ndeg+1
         do j = 1,ndeg+1
            ip2ind(j,i) = -1
         enddo
      enddo
      
      if (type .eq. 'f' .or. type .eq. 'F') then
         ipol = 0
         do i = 1,ndeg+1
            do j = 1,ndeg+1
               ipol = ipol+1
               ip2ind(j,i) = ipol
            enddo
         enddo
      else if (type .eq. 't' .or. type .eq. 'T') then
         ipol = 0
         do i = 1,ndeg+1
            do j = 1,ndeg+1+1-i
               ipol = ipol+1
               ip2ind(j,i) = ipol
            enddo
         enddo
      endif

      return
      end
c
c
c
c
c
      subroutine polytens_ind2pow_2d(ndeg,type,iind2p)
      integer *8 ndeg, iind2p(2,*)
      character type

      integer *8 i, j, ipol


      if (type .eq. 'f' .or. type .eq. 'F') then
         ipol = 0
         do i = 1,ndeg+1
            do j = 1,ndeg+1
               ipol = ipol+1
               iind2p(1,ipol) = j-1
               iind2p(2,ipol) = i-1                
            enddo
         enddo
      else if (type .eq. 't' .or. type .eq. 'T') then
         ipol = 0
         do i = 1,ndeg+1
            do j = 1,ndeg+1+1-i
               ipol = ipol+1
               iind2p(1,ipol) = j-1
               iind2p(2,ipol) = i-1
            enddo
         enddo
      endif

      return
      end
c
c
      subroutine polytens_exps_nd(ndim,ipoly,itype,n,type,x,
     1    u,ldu,v,ldv,w)
c                 input parameters:
c
c  itype - the type of the calculation to be performed
c          itype=0 means that only the gaussian nodes are 
c                  to be constructed. 
c          itype=1 means that only the nodes and the weights 
c                  are to be constructed
c          itype=2 means that the nodes, the weights, and
c                  the matrices u, v are to be constructed
c          itype=3 only construct u
c          itype=4 only construct v
      
      implicit none
      integer *8 ndim,ipoly, itype, n, ldu, ldv
      character type
      real *8 x(ndim,*),w(*)
      real *8 u(ldu,*), v(ldv,*)

      if (ndim.eq.1) then
         if (ipoly.eq.0) then
            call legeexps(itype,n,x,u,v,w)
         elseif (ipoly.eq.1) then
            call chebexps(itype,n,x,u,v,w)
         endif
      elseif (ndim.eq.2) then
         call polytens_exps_2d(ipoly,itype,n,type,x,u,ldu,v,ldv,w)
      elseif (ndim.eq.3) then
         call polytens_exps_3d(ipoly,itype,n,type,x,u,ldu,v,ldv,w)
      endif
      
      return
      end
c
c
      subroutine polytens_exps_2d(ipoly,itype,n,type,x,u,ldu,v,ldv,w)
c                 input parameters:
c
c  itype - the type of the calculation to be performed
c          itype=0 means that only the gaussian nodes are 
c                  to be constructed. 
c          itype=1 means that only the nodes and the weights 
c                  are to be constructed
c          itype=2 means that the nodes, the weights, and
c                  the matrices u, v are to be constructed
c          itype=3 only construct u
c          itype=4 only construct v
      
      implicit none
      integer *8 ipoly, itype, n, ldu, ldv
      character type
      real *8 x(2,*),w(*)
      real *8 u(ldu,*), v(ldv,*)
      real *8 x1d(n), w1d(n), u1d(n,n), v1d(n,n)
      integer *8 i,j,ipt,itype1d,ii, jo,ipol
      
      itype1d = 0
      if (itype .ge. 1) then
         itype1d = 1
      endif
      if (itype .ge. 2) then
         itype1d = 2
      endif
      
      if (ipoly.eq.0) call legeexps(itype1d,n,x1d,u1d,v1d,w1d)
      if (ipoly.eq.1) call chebexps(itype1d,n,x1d,u1d,v1d,w1d)

      ipt = 0
      do i=1,n
         do j=1,n
            ipt = ipt + 1
            x(1,ipt) = x1d(j)
            x(2,ipt) = x1d(i)               
        enddo
      enddo

      if (itype .ge. 1) then
         ipt = 0
         do i=1,n
            do j=1,n
               ipt = ipt + 1
               w(ipt) = w1d(i)*w1d(j)
            enddo
         enddo
      endif



      if (itype .eq. 2 .or. itype .eq. 3) then
c     construct u from 1d u
         if (type .eq. 'f' .or. type .eq. 'F') then        
            ipt = 0
            do ii = 1,n
            do jo = 1,n
               ipt = ipt + 1
               ipol = 0
               do i=1,n
               do j=1,n
                  ipol = ipol + 1
                  u(ipol,ipt) =  u1d(i,ii)*u1d(j,jo)
               enddo
               enddo
            enddo
            enddo


         else if (type .eq. 't' .or. type .eq. 'T') then

            ipt = 0
            do ii = 1,n
            do jo = 1,n
               ipt = ipt + 1
               ipol = 0
               do i=1,n
               do j=1,n+1-i
                  ipol = ipol + 1
                  u(ipol,ipt) = u1d(i,ii)*u1d(j,jo)
               enddo
               enddo
            enddo
            enddo
         endif
      endif
      
      if (itype .eq. 2 .or. itype .eq. 4) then
c     construct v from 1d v
         if (type .eq. 'f' .or. type .eq. 'F') then         
            ipol = 0
            do ii = 1,n
            do jo = 1,n
               ipol = ipol + 1
               ipt = 0
               do i=1,n
               do j=1,n
                  ipt = ipt + 1
                  v(ipt,ipol) =  v1d(i,ii)*v1d(j,jo)
               enddo
               enddo
            enddo
            enddo

         else if (type .eq. 't' .or. type .eq. 'T') then

            ipol = 0
            do ii = 1,n
            do jo = 1,n+1-ii
               ipol = ipol + 1
               ipt = 0
               do i=1,n
               do j=1,n
                  ipt = ipt + 1
                  v(ipt,ipol) = v1d(i,ii)*v1d(j,jo)
               enddo
               enddo
            enddo
            enddo
         endif
      endif         

      return
      end
c
c
      subroutine polytens_exps_3d(ipoly,itype,n,type,x,u,ldu,v,ldv,w)
c                 input parameters:
c
c  itype - the type of the calculation to be performed
c          itype=0 means that only the gaussian nodes are 
c                  to be constructed. 
c          itype=1 means that only the nodes and the weights 
c                  are to be constructed
c          itype=2 means that the nodes, the weights, and
c                  the matrices u, v are to be constructed
c          itype=3 only construct x,w, and u
c          itype=4 only construct x,w, and v
      
      implicit none
      integer *8 ipoly, itype, n, ldu, ldv
      character type
      real *8 x(3,*),w(*)
      real *8 u(ldu,*), v(ldv,*)
      real *8 x1d(n), w1d(n), u1d(n,n), v1d(n,n)
      integer *8 i,j,ipt,itype1d, k, io, jo, ko, ipol
      
      itype1d = 0
      if (itype .ge. 1) then
         itype1d = 1
      endif
      if (itype .ge. 2) then
         itype1d = 2
      endif
      
      if (ipoly.eq.0) call legeexps(itype1d,n,x1d,u1d,v1d,w1d)
      if (ipoly.eq.1) call chebexps(itype1d,n,x1d,u1d,v1d,w1d)

      ipt = 0
      do i=1,n
         do j=1,n
            do k = 1,n
               ipt = ipt + 1
               x(1,ipt) = x1d(k)
               x(2,ipt) = x1d(j)               
               x(3,ipt) = x1d(i)
            enddo
        enddo
      enddo

      if (itype .ge. 1) then
         ipt = 0
         do i=1,n
            do j=1,n
               do k = 1,n
                  ipt = ipt + 1
                  w(ipt) = w1d(i)*w1d(j)*w1d(k) 
               enddo
            enddo
         enddo
      endif


      if (itype .eq. 2 .or. itype .eq. 3) then
c     construct u from 1d u
         if (type .eq. 'f' .or. type .eq. 'F') then         
            ipt = 0
            do io = 1,n
            do jo = 1,n
            do ko = 1,n
               ipt = ipt + 1
               ipol = 0
               do i=1,n
               do j=1,n
               do k = 1,n
                  ipol = ipol + 1
                  u(ipol,ipt) =  u1d(i,io)*u1d(j,jo)*u1d(k,ko)
               enddo
               enddo
               enddo
            enddo
            enddo
            enddo

         else if (type .eq. 't' .or. type .eq. 'T') then

            ipt = 0
            do io = 1,n
            do jo = 1,n
            do ko = 1,n
               ipt = ipt + 1
               ipol = 0
               do i=1,n
               do j=1,n+1-i
               do k = 1,n+2-i-j
                  ipol = ipol + 1
                  u(ipol,ipt) = u1d(i,io)*u1d(j,jo)*u1d(k,ko)
               enddo
               enddo
               enddo
            enddo
            enddo
            enddo
         endif
      endif
      
      if (itype .eq. 2 .or. itype .eq. 4) then
c     construct v from 1d v
         if (type .eq. 'f' .or. type .eq. 'F') then         
            ipol = 0
            do io = 1,n
            do jo = 1,n
            do ko = 1,n
               ipol = ipol + 1
               ipt = 0
               do i=1,n
               do j=1,n
               do k=1,n
                  ipt = ipt + 1
                  v(ipt,ipol) =  v1d(i,io)*v1d(j,jo)*v1d(k,ko)
               enddo
               enddo
               enddo
            enddo
            enddo
            enddo

         else if (type .eq. 't' .or. type .eq. 'T') then

            ipol = 0
            do io = 1,n
            do jo = 1,n+1-io
            do ko = 1,n+2-io-jo
               ipol = ipol + 1
               ipt = 0
               do i=1,n
               do j=1,n
               do k=1,n
                  ipt = ipt + 1
                  v(ipt,ipol) = v1d(i,io)*v1d(j,jo)*v1d(k,ko)
               enddo
               enddo
               enddo
            enddo
            enddo
            enddo
         endif
      endif         

      return
      end
c
c

      subroutine polytens_pols_2d(ipoly,x,ndeg,type,pols)
c
c     type = 'F' full degree polynomials ((ndeg+1)**3)
c     type = 'T' total degree polynomials ((ndeg+1)*(ndeg+2)*(ndeg+3)/6)
c      
      implicit none
      integer *8 ipoly,ndeg,npols
      real *8 x(2),pols(*),px(ndeg+1),py(ndeg+1)
      integer *8 i,j,ipol,n
      character type

      n = ndeg + 1

      if (ipoly.eq.0) then
         call legepols(x(1),ndeg,px)
         call legepols(x(2),ndeg,py)
      elseif (ipoly.eq.1) then
         call chebpols(x(1),ndeg,px)
         call chebpols(x(2),ndeg,py)
      endif
      
      if (type .eq. 'f' .or. type .eq. 'F') then
         ipol = 0
         do i=1,n
            do j =1,n
               ipol = ipol + 1
               pols(ipol) = px(j)*py(i)
            enddo
         enddo
      else if (type .eq. 't' .or. type .eq. 'T') then
         ipol = 0
         do i=1,n
            do j=1,n+1-i
               ipol = ipol + 1
               pols(ipol) = px(j)*py(i)
            enddo
         enddo
      endif
      
      return
      end

      subroutine polytens_pols_3d(ipoly,x,ndeg,type,pols)
c
c     type = 'F' full degree polynomials ((ndeg+1)**3)
c     type = 'T' total degree polynomials ((ndeg+1)*(ndeg+2)*(ndeg+3)/6)
c      
      implicit none
      integer *8 ipoly,ndeg,npols
      real *8 x(3),pols(*),px(ndeg+1),py(ndeg+1),pz(ndeg+1)
      integer *8 i,j,ipol, k, n
      character type

      n = ndeg + 1

      if (ipoly.eq.0) then
         call legepols(x(1),ndeg,px)
         call legepols(x(2),ndeg,py)
         call legepols(x(3),ndeg,pz)      
      elseif (ipoly.eq.1) then
         call chebpols(x(1),ndeg,px)
         call chebpols(x(2),ndeg,py)
         call chebpols(x(3),ndeg,pz)      
      endif
      
      if (type .eq. 'f' .or. type .eq. 'F') then
         ipol = 0
         do i=1,n
            do j=1,n
               do k = 1,n
                  ipol = ipol + 1
                  pols(ipol) = px(k)*py(j)*pz(i)
               enddo
            enddo
         enddo
      else if (type .eq. 't' .or. type .eq. 'T') then
         ipol = 0
         do i=1,n
            do j=1,n+1-i
               do k = 1,n+2-i-j
                  ipol = ipol + 1
                  pols(ipol) = px(k)*py(j)*pz(i)
               enddo
            enddo
         enddo
      endif
      
      return
      end

      
      subroutine polytens_ders_2d(ipoly,x,ndeg,type,pols,ders)
c
c     type = 'F' full degree polynomials ((ndeg+1)**3)
c     type = 'T' total degree polynomials ((ndeg+1)*(ndeg+2)*(ndeg+3)/6)
c      
      implicit none
      integer *8 ipoly,ndeg,npols,idim
      real *8 x(2),ders(2,*),px(ndeg+1),py(ndeg+1),tmp(ndeg+1)
      real *8 pxder(ndeg+1),pyder(ndeg+1),pols(*)
      integer *8 i,j,ipol, k, n
      character type

      n = ndeg + 1

      if (ipoly.eq.0) then
         call legepolders(x(1),px,pxder,ndeg)
         call legepolders(x(2),py,pyder,ndeg)
      elseif (ipoly.eq.1) then
         call chebpolders(x(1),px,pxder,ndeg)
         call chebpolders(x(2),py,pyder,ndeg)
      endif

      if (type .eq. 'f' .or. type .eq. 'F') then
         ipol = 0
         do i=1,n
            do j=1,n
               ipol = ipol + 1
               pols(ipol) = px(j)*py(i)
               ders(1,ipol) = pxder(j)*py(i)
               ders(2,ipol) = px(j)*pyder(i)
            enddo
         enddo
      else if (type .eq. 't' .or. type .eq. 'T') then
         ipol = 0
         do i=1,n
            do j=1,n+1-i
               ipol = ipol + 1
               pols(ipol) = px(j)*py(i)
               ders(1,ipol) = pxder(j)*py(i)
               ders(2,ipol) = px(j)*pyder(i)
            enddo
         enddo
      endif
      
      return
      end

      
      subroutine polycoeff_dmat(ipoly,ndeg,dmat,lddmat)
c
c     
c      

      implicit none
      integer *8 ipoly, ndeg, lddmat
      real *8 dmat(lddmat,ndeg+1)
c     local
      real *8 polin(ndeg+1), polout(ndeg+1)
      integer *8 i, j

      do i = 1,ndeg+1
         do j = 1,ndeg+1
            polin(j) = 0.0d0
            if (j .eq. i) polin(j) = 1.0d0
         enddo
         if (ipoly.eq.0) call legediff(polin,ndeg,polout)
         if (ipoly.eq.1) call chebdiff(polin,ndeg,polout)
         do j = 1,ndeg
            dmat(j,i) = polout(j)
         enddo
         do j = ndeg+1,lddmat
            dmat(j,i) = 0
         enddo
      enddo

      return
      end
      

      subroutine polycoeff_d2mat(ipoly,ndeg,dmat,lddmat)
      implicit none
      integer *8 ipoly, ndeg, lddmat
      real *8 dmat(lddmat,ndeg+1)
c     local
      real *8 polin(ndeg+1), polout(ndeg+1)
      integer *8 i, j

      do i = 1,ndeg+1
         do j = 1,ndeg+1
            polin(j) = 0.0d0
            if (j .eq. i) polin(j) = 1.0d0
         enddo
         if (ipoly.eq.0) then
            call legediff(polin,ndeg,polout)
            call legediff(polout,ndeg-1,polin)         
         elseif (ipoly.eq.1) then
            call chebdiff(polin,ndeg,polout)
            call chebdiff(polout,ndeg-1,polin)
         endif
         do j = 1,ndeg-1
            dmat(j,i) = polin(j)
         enddo
         do j = ndeg,lddmat
            dmat(j,i) = 0
         enddo
      enddo

      return
      end

      subroutine polytens_lapmat_2d(ipoly,ndegin,ndegout,type,
     1     lapmat,ldlapmat)
      implicit none
      integer *8 ipoly, ndegin, ndegout, ldlapmat
      character type
      real *8 lapmat(ldlapmat,*)
c     local
      real *8 lap1d(ndegout+1,ndegin+1)
      integer *8 ldlap1d, io, ii, jo, ji, iopol, iipol

      do ii = 1,ndegin+1
         do io = 1,ndegout+1
            lap1d(io,ii) = 0.0d0
         enddo
      enddo
      
      
      ldlap1d = ndegout+1

      call polycoeff_d2mat(ipoly,ndegin,lap1d,ldlap1d)

      if (type .eq. 'f' .or. type .eq. 'F') then
         iipol = 0
         do ii = 1,ndegin+1
            do ji = 1,ndegin+1
               iipol = iipol + 1
               iopol = 0
               do io = 1,ndegout+1
                  do jo = 1,ndegout+1
                     iopol = iopol + 1
                     lapmat(iopol,iipol) = 0.0d0
                     if (io .eq. ii) then
                        lapmat(iopol,iipol) =
     1                       lapmat(iopol,iipol) + lap1d(jo,ji)
                     endif
                     if (jo .eq. ji) then
                        lapmat(iopol,iipol) =
     1                       lapmat(iopol,iipol) + lap1d(io,ii)
                     endif
                  enddo
               enddo
            enddo
         enddo
      else if (type .eq. 't' .or. type .eq. 'T') then
         iipol = 0
         do ii = 1,ndegin+1
            do ji = 1,ndegin+1+1-ii
               iipol = iipol + 1
               iopol = 0
               do io = 1,ndegout+1
                  do jo = 1,ndegout+1+1-io
                     iopol = iopol + 1
                     lapmat(iopol,iipol) = 0.0d0
                     if (io .eq. ii) then
                        lapmat(iopol,iipol) =
     1                       lapmat(iopol,iipol) + lap1d(jo,ji)
                     endif
                     if (jo .eq. ji) then
                        lapmat(iopol,iipol) =
     1                       lapmat(iopol,iipol) + lap1d(io,ii)
                     endif
                     
                  enddo
               enddo
            enddo
         enddo
      endif
      
      return
      end

      subroutine polytens_lapmat_3d(ipoly,ndegin,ndegout,type,
     1     lapmat,ldlapmat)
      implicit none
      integer *8 ipoly, ndegin, ndegout, ldlapmat
      character type
      real *8 lapmat(ldlapmat,*)
c     local
      real *8 lap1d(ndegout+1,ndegin+1)
      integer *8 ldlap1d, io, ii, jo, ji, ko, ki, iopol, iipol

      ldlap1d = ndegout+1

      do ii = 1,ndegin+1
         do io = 1,ndegout+1
            lap1d(io,ii) = 0.0d0
         enddo
      enddo
      
      call polycoeff_d2mat(ipoly,ndegin,lap1d,ldlap1d)

      if (type .eq. 'f' .or. type .eq. 'F') then
         iipol = 0
         do ii = 1,ndegin+1
            do ji = 1,ndegin+1
               do ki = 1,ndegin+1
                  iipol = iipol + 1
                  iopol = 0
                  do io = 1,ndegout+1
                     do jo = 1,ndegout+1
                        do ko = 1,ndegout+1
                           iopol = iopol + 1
                           lapmat(iopol,iipol) = 0.0d0
                           if (jo .eq. ji .and. io .eq. ii) then
                              lapmat(iopol,iipol) =
     1                             lapmat(iopol,iipol) + lap1d(ko,ki)
                           endif
                           if (ko .eq. ki .and. io .eq. ii) then
                              lapmat(iopol,iipol) =
     1                             lapmat(iopol,iipol) + lap1d(jo,ji)
                           endif
                           if (jo .eq. ji .and. ko .eq. ki) then
                              lapmat(iopol,iipol) =
     1                             lapmat(iopol,iipol) + lap1d(io,ii)
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      else if (type .eq. 't' .or. type .eq. 'T') then
         iipol = 0
         do ii = 1,ndegin+1
            do ji = 1,ndegin+1+1-ii
               do ki = 1,ndegin+1+2-ii-ji
                  iipol = iipol + 1
                  iopol = 0
                  do io = 1,ndegout+1
                     do jo = 1,ndegout+1+1-io
                        do ko = 1,ndegout+1+2-io-jo
                           iopol = iopol + 1
                           lapmat(iopol,iipol) = 0.0d0
                           if (jo .eq. ji .and. io .eq. ii) then
                              lapmat(iopol,iipol) =
     1                             lapmat(iopol,iipol) + lap1d(ko,ki)
                           endif
                           if (ko .eq. ki .and. io .eq. ii) then
                              lapmat(iopol,iipol) =
     1                             lapmat(iopol,iipol) + lap1d(jo,ji)
                           endif
                           if (jo .eq. ji .and. ko .eq. ki) then
                              lapmat(iopol,iipol) =
     1                             lapmat(iopol,iipol) + lap1d(io,ii)
                           endif

                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
      
      return
      end


      subroutine polytens_eyemat_2d(ndegin,ndegout,type,
     1     eyemat,ldeyemat)
      implicit none
      integer *8 ndegin, ndegout, ldeyemat
      character type
      real *8 eyemat(ldeyemat,*)
c     local
      integer *8 io, ii, jo, ji, iopol, iipol

      if (type .eq. 'f' .or. type .eq. 'F') then
         iipol = 0
         do ii = 1,ndegin+1
            do ji = 1,ndegin+1
               iipol = iipol + 1
               iopol = 0
               do io = 1,ndegout+1
                  do jo = 1,ndegout+1
                     iopol = iopol + 1
                     eyemat(iopol,iipol) = 0.0d0
                     if (jo .eq. ji .and. io .eq. ii) then
                        eyemat(iopol,iipol) = 1.0d0
                     endif
                  enddo
               enddo
            enddo
         enddo
      else if (type .eq. 't' .or. type .eq. 'T') then
         iipol = 0
         do ii = 1,ndegin+1
            do ji = 1,ndegin+1+1-ii
               iipol = iipol + 1
               iopol = 0
               do io = 1,ndegout+1
                  do jo = 1,ndegout+1+1-io
                     iopol = iopol + 1
                     eyemat(iopol,iipol) = 0.0d0
                     if (jo .eq. ji .and. io .eq. ii) then
                        eyemat(iopol,iipol) = 1.0d0
                     endif
                  enddo
               enddo
            enddo
         enddo
      endif

      return
      end
      
      subroutine polytens_eyemat_3d(ndegin,ndegout,type,
     1     eyemat,ldeyemat)
      implicit none
      integer *8 ndegin, ndegout, ldeyemat
      character type
      real *8 eyemat(ldeyemat,*)
c     local
      integer *8 io, ii, jo, ji, ko, ki, iopol, iipol

      if (type .eq. 'f' .or. type .eq. 'F') then
         iipol = 0
         do ii = 1,ndegin+1
            do ji = 1,ndegin+1
               do ki = 1,ndegin+1
                  iipol = iipol + 1
                  iopol = 0
                  do io = 1,ndegout+1
                     do jo = 1,ndegout+1
                        do ko = 1,ndegout+1
                           iopol = iopol + 1
                           eyemat(iopol,iipol) = 0.0d0
                           if (jo .eq. ji .and. io .eq. ii .and.
     1                          ko .eq. ki) then
                              eyemat(iopol,iipol) = 1.0d0
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      else if (type .eq. 't' .or. type .eq. 'T') then
         iipol = 0
         do ii = 1,ndegin+1
            do ji = 1,ndegin+1+1-ii
               do ki = 1,ndegin+1+2-ii-ji
                  iipol = iipol + 1
                  iopol = 0
                  do io = 1,ndegout+1
                     do jo = 1,ndegout+1+1-io
                        do ko = 1,ndegout+1+2-io-jo
                           iopol = iopol + 1
                           eyemat(iopol,iipol) = 0.0d0
                           if (jo .eq. ji .and. io .eq. ii .and.
     1                          ko .eq. ki) then
                              eyemat(iopol,iipol) = 1.0d0
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
      
      return
      end
      

      subroutine polytens_slicezcoeffs_3d(ipoly,ndeg,type,idim,val,cin,
     1     ldcin,nfun,cout,ldcout,ifdiff,dcout,lddcout)
c
c     complex valued coeffs, for real see legetens_slicedcoeffs...
c      
c     given a set of polynomials defined by the coefficients in
c     cin, get the 2d coefficients of the slice of that polynomial
c     where x_idim = val
c     
      implicit none
      integer *8 ipoly,ndeg, idim, ldcin, nfun, ldcout, ifdiff, lddcout
      character type
      real *8 val
      complex *16 cin(ldcin,*), cout(ldcout,*), dcout(lddcout,*)
c     local
      integer *8 iipol, iopol, ii, ji, ki, io, jo, ko, n, npol2, jpol
      integer *8 ilist(ndeg+1,ndeg+1), ifun
      real *8 pols(ndeg+1), ders(ndeg+1)
      complex *16 cv, zero, dv
      data zero / (0.0d0,0.0d0) /

      if (ipoly.eq.0) call legepolders(val,pols,ders,ndeg)
      if (ipoly.eq.1) call chebpolders(val,pols,ders,ndeg)

      call polytens_npol_2d(ipoly,ndeg,type,npol2)

      n = ndeg+1

      if (type .eq. 'f' .or. type .eq. 'F') then
         
         iopol = 0
         do io = 1,n
            do jo = 1,n
               iopol = iopol+1
               ilist(jo,io) = iopol
            enddo
         enddo

         do ifun = 1,nfun

            do jpol = 1,npol2
               cout(jpol,ifun) = zero
            enddo

            if (ifdiff .eq. 1) then
               do jpol = 1,npol2
                  dcout(jpol,ifun) = zero
               enddo
            endif
            
            iipol = 0
            do ii = 1,n
               do ji = 1,n
                  do ki = 1,n
                     iipol = iipol+1
                     if (idim .eq. 1) then
                        io = ii
                        jo = ji
                        cv = pols(ki)*cin(iipol,ifun)
                        dv = ders(ki)*cin(iipol,ifun)
                     else if (idim .eq. 2) then
                        io = ii
                        jo = ki
                        cv = pols(ji)*cin(iipol,ifun)
                        dv = ders(ji)*cin(iipol,ifun)                        
                     else if (idim .eq. 3) then
                        io = ji
                        jo = ki
                        cv = pols(ii)*cin(iipol,ifun)
                        dv = ders(ii)*cin(iipol,ifun)
                     endif
                     iopol = ilist(jo,io)
                     cout(iopol,ifun) = cout(iopol,ifun)+cv
                     if (ifdiff .eq. 1) then
                        dcout(iopol,ifun) = dcout(iopol,ifun)+dv
                     endif
                  enddo
               enddo
            enddo
         enddo
         
      endif
      
      if (type .eq. 't' .or. type .eq. 'T') then
         
         iopol = 0
         do io = 1,n
            do jo = 1,n+1-io
               iopol = iopol+1
               ilist(jo,io) = iopol
            enddo
         enddo

         do ifun = 1,nfun

            do jpol = 1,npol2
               cout(jpol,ifun) = zero
            enddo

            if (ifdiff .eq. 1) then
               do jpol = 1,npol2
                  dcout(jpol,ifun) = zero
               enddo
            endif
            
            iipol = 0
            do ii = 1,n
               do ji = 1,n+1-ii
                  do ki = 1,n+2-ii-ji
                     iipol = iipol+1
                     if (idim .eq. 1) then
                        io = ii
                        jo = ji
                        cv = pols(ki)*cin(iipol,ifun)
                        dv = ders(ki)*cin(iipol,ifun)                        
                     else if (idim .eq. 2) then
                        io = ii
                        jo = ki
                        cv = pols(ji)*cin(iipol,ifun)
                        dv = ders(ji)*cin(iipol,ifun)                        
                     else if (idim .eq. 3) then
                        io = ji
                        jo = ki
                        cv = pols(ii)*cin(iipol,ifun)
                        dv = ders(ii)*cin(iipol,ifun)                        
                     endif
                     iopol = ilist(jo,io)
                     cout(iopol,ifun) = cout(iopol,ifun)+cv
                     if (ifdiff .eq. 1) then
                        dcout(iopol,ifun) = dcout(iopol,ifun)+dv
                     endif
                  enddo
               enddo
            enddo
         enddo
         
      endif
      
      return
      end



      subroutine polytens_lape_3d(ipoly,nd,ndeg,type,polin,polout)
c
c     compute the Laplacian of the input coefficients, returning
c     in the output coefficients
c
c     input:
c     
c     nd - dimensionality of the coefficients
c     polin(nd,*) - real *8 array. nd sets of tensor Legendre
c       coefficients in the legetens ordering
c     ndeg - integer *8, degree of polynomial basis
c     type - character, 'T' gives total degree polynomials
c       'F' gives full degree
c
c     output:
c
c     polout(nd,*) - real *8 array. nd sets of tensor Legendre
c       coefficients in the legetens ordering giving the
c       Laplacian of the input coefficients (of same degree)
c      

      implicit none
      integer *8 :: ipoly, nd, ndeg
      character :: type
      real *8 :: polin(nd,*), polout(nd,*)
c     local
      real *8 :: dp, d2i, d2j, d2k
      real *8, allocatable :: dxxmat(:,:)
      integer *8 :: n, i, j, npol, iii, ip, jp, kp, id, jd, kd, ind
      integer *8 :: jjj
      integer *8, allocatable :: ip2ind(:,:,:), ind2p(:,:)

      n = ndeg+1
      allocate(dxxmat(0:ndeg,0:ndeg))
      call polycoeff_d2mat(ipoly,ndeg,dxxmat,n)

      call polytens_npol_3d(ipoly,ndeg,type,npol)

      n = ndeg + 1
      
      allocate(ip2ind(0:ndeg,0:ndeg,0:ndeg),ind2p(3,npol))
      
      call polytens_pow2ind_3d(ndeg,type,ip2ind)
      call polytens_ind2pow_3d(ndeg,type,ind2p)

      do i = 1,npol
         do j = 1,nd
            polout(j,i) = 0.0d0
         enddo
      enddo

      if (ndeg .lt. 2) return

c     grab appropriate entries
      
      do iii = 1,npol
         ip = ind2p(1,iii)
         jp = ind2p(2,iii)
         kp = ind2p(3,iii)

         do kd = 0,kp-2
            d2k = dxxmat(kd,kp)
            ind = ip2ind(ip,jp,kd)
            do jjj = 1,nd
               dp = polin(jjj,iii)
               polout(jjj,ind) = polout(jjj,ind) + dp*d2k
            enddo
         enddo
         
         do jd = 0,jp-2
            d2j = dxxmat(jd,jp)
            ind = ip2ind(ip,jd,kp)
            do jjj = 1,nd
               dp = polin(jjj,iii)
               polout(jjj,ind) = polout(jjj,ind) + dp*d2j
            enddo
         enddo
         
         do id = 0,ip-2
            d2i = dxxmat(id,ip)
            ind = ip2ind(id,jp,kp)
            do jjj = 1,nd
               dp = polin(jjj,iii)
               polout(jjj,ind) = polout(jjj,ind) + dp*d2i
            enddo
         enddo
      enddo

      return
      end
      
    

      subroutine polydiff_3d(ipoly,ndeg,ttype,polin,idir,dmat,polout)
c
c       This subroutine differentiates a tensor
c       product legendre expansion expressed in 
c       its coefficient basis and returns
c       the tensor product legendre expansion
c       coeffs of the derivative
c 
c       input arguments:
c
c     ndeg refers to polynomial degree
c      
c     n refers to order of approximation and number of nodes
c     used in each direction (n = ndeg+1)
c
c     type character specifies the set of tensor
c     polynomials to use. current options
c
c     ttype = 'F', full degree polynomials (T_i(x)T_j(y)T_k(z) where
c                 each of i, j, and k goes from 0 to ndeg)
c     ttype = 'T', total degree polynomials (T_i(x)T_j(y)T_k(z) where
c                 sum of i,j,k>=0 is less than or equal to ndeg)
c     
c     polin - real *8 (*)
c       tensor product legendre expansion coeffs of 
c       input polynomial
c
c     dmat - real *8 (ndeg+1,ndeg+1)
c       differentiation martix in 1d
c
c     idir - integer *8
c       whether to compute x,y,or z derivative of expansion
c       idir = 1, compute x derivative
c       idir = 2, compute y derivative
c       idir = 3, compute z derivative
c
c
c
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      real *8 dmat(ndeg+1,ndeg+1),polin(*),polout(*)
      real *8 coef1(ndeg+1,ndeg+1,ndeg+1),coef2(ndeg+1,ndeg+1,ndeg+1)
      integer *8 iind2pow(3,(ndeg+1)**3)
      character ttype


      call polytens_ind2pow_3d(ndeg,ttype,iind2pow)
      call polytens_npol_3d(ndeg,ttype,npol)
c
c   extract coef mat
c
c

      coef1 = 0
      coef2 = 0
      do i=1,npol
        i1 = iind2pow(1,i)+1
        i2 = iind2pow(2,i)+1
        i3 = iind2pow(3,i)+1
        coef1(i1,i2,i3) = polin(i)
      enddo

      if(idir.eq.1) then
        do j=1,ndeg+1
          do k=1,ndeg+1
            do l=1,ndeg+1
              do i=1,ndeg+1
                coef2(l,k,j) = coef2(l,k,j) + dmat(l,i)*coef1(i,k,j)
              enddo
            enddo
          enddo
        enddo
      endif

      if(idir.eq.2) then
        do j=1,ndeg+1
          do k=1,ndeg+1
            do l=1,ndeg+1
              do i=1,ndeg+1
                coef2(k,l,j) = coef2(k,l,j) + dmat(l,i)*coef1(k,i,j)
              enddo
            enddo
          enddo
        enddo
      endif

      if(idir.eq.3) then
        do j=1,ndeg+1
          do k=1,ndeg+1
            do l=1,ndeg+1
              do i=1,ndeg+1
                coef2(k,j,l) = coef2(k,j,l) + dmat(l,i)*coef1(k,j,i)
              enddo
            enddo
          enddo
        enddo
      endif

      do i=1,npol
        i1 = iind2pow(1,i)+1
        i2 = iind2pow(2,i)+1
        i3 = iind2pow(3,i)+1
        polout(i) = coef2(i1,i2,i3) 
      enddo


      return
      end
c
c
c
c
      subroutine tens_prod_trans_nd(ndim,nd,norder,fvals,fcoefs,umat)
c
c     user interface for tensor product transformation in n dimensions
c
c
c
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      real *8 fvals(nd,norder**ndim)
      real *8 fcoefs(nd,norder**ndim),umat(norder,norder)

      if (ndim.eq.1) then
         call tens_prod_trans_1d(nd,norder,fvals,fcoefs,umat)
      elseif (ndim.eq.2) then
         call tens_prod_trans_2d(nd,norder,fvals,fcoefs,umat)
      elseif (ndim.eq.3) then
         call tens_prod_trans_3d(nd,norder,fvals,fcoefs,umat)
      endif
      
      return
      end
c
c
c
c
      subroutine tens_prod_trans_1d(nd,norder,fvals,fcoefs,umat)
c
c     1D tensor product transformation
c
c
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      real *8 fvals(nd,norder)
      real *8 fcoefs(nd,norder),umat(norder,norder)
      
      do j=1,norder
         do ind=1,nd
            dd=0
            do j1=1,norder
               dd=dd+umat(j,j1)*fvals(ind,j1)
            enddo
            fcoefs(ind,j)=dd
         enddo
      enddo

      return
      end
c
c
c
c
      subroutine tens_prod_trans_2d(nd,norder,fvals,fcoefs,umat)
c
c     2D tensor product transformation
c
c
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      real *8 fvals(nd,norder,norder)
      real *8 fcoefs(nd,norder,norder),umat(norder,norder)
      real *8, allocatable:: fcv(:,:,:)

      allocate(fcv(nd,norder,norder))
      
      do j=1,norder
         do k=1,norder
            do ind=1,nd
               dd=0
               do k1=1,norder
                  dd=dd+umat(k,k1)*fvals(ind,k1,j)
               enddo
               fcv(ind,k,j)=dd
            enddo
         enddo
      enddo

      do j=1,norder
         do k=1,norder
            do ind=1,nd
               dd=0
               do j1=1,norder
                  dd=dd+umat(j,j1)*fcv(ind,k,j1)
               enddo
               fcoefs(ind,k,j)=dd
            enddo
         enddo
      enddo

      return
      end
c
c
c
c
      subroutine tens_prod_trans_3d(nd,norder,fvals,fcoefs,umat)
c
c     3D tensor product transformation
c
c
c
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      real *8 fvals(nd,norder,norder,norder)
      real *8 fcoefs(nd,norder,norder,norder),umat(norder,norder)
      real *8, allocatable:: fcvv(:,:,:,:),fccv(:,:,:,:)

      allocate(fcvv(nd,norder,norder,norder))
      allocate(fccv(nd,norder,norder,norder))
      
      do i=1,norder
         do j=1,norder
            do k=1,norder
               do ind=1,nd
                  dd=0
                  do k1=1,norder
                     dd=dd+umat(k,k1)*fvals(ind,k1,j,i)
                  enddo
                  fcvv(ind,k,j,i)=dd
               enddo
            enddo
         enddo
      enddo

      do i=1,norder
         do j=1,norder
            do k=1,norder
               do ind=1,nd
                  dd=0
                  do j1=1,norder
                     dd=dd+umat(j,j1)*fcvv(ind,k,j1,i)
                  enddo
                  fccv(ind,k,j,i)=dd
               enddo
            enddo
         enddo
      enddo

      do i=1,norder
         do j=1,norder
            do k=1,norder
               do ind=1,nd
                  dd=0
                  do i1=1,norder
                     dd=dd+umat(i,i1)*fccv(ind,k,j,i1)
                  enddo
                  fcoefs(ind,k,j,i)=dd
               enddo
            enddo
         enddo
      enddo

      return
      end
c
c
c
c
      subroutine tens_prod_get_rmask(ndim,iptype,norder,npols,
     1    rmask,rsum)
c
c     get the tail indices for estimating the function error
c
c
c
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8 iind2p(ndim,npols)
      real *8 rmask(npols)

      call polytens_ind2pow(ndim,norder-1,'f',iind2p)
      
      rsum = 0
      do i=1,npols
        rmask(i) = 0.0d0
        i1=0
        do k=1,ndim
           i1=i1+iind2p(k,i)
        enddo
        if(i1.eq.norder-1) then
          rmask(i) = 1.0d0
          rsum = rsum + 1
        endif
      enddo

      if(iptype.eq.2) rsum = sqrt(rsum)
      if(iptype.eq.0) rsum = 1

      return
      end
c
c
c
c
      
