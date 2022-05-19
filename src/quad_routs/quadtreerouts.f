
      
c
c
c

      subroutine getquadchildren(v0,v1,v2,v3,v4)
c  
cc       given the three vertices of a quad v0,
c        this subroutine returns the vertices of 4 
c        smaller quads constructed using the
c        midpoints of the quads
c 
c        input:
c        v0 - real *8 (2,3)  
c              vertices of parent quad
c
c        output:
c        v1,v2,v3,v4 - real *8 (2,3)
c                 vertices of children quads
      
      implicit real *8 (a-h,o-z)
      real *8 v0(2,3),v1(2,3),v2(2,3),v3(2,3),v4(2,3),vm(2,5)


      vm(1,1) = (v0(1,1)+v0(1,2))/2
      vm(2,1) = v0(2,1)

      vm(1,2) = (v0(1,3)+v0(1,2))/2
      vm(2,2) = (v0(2,3)+v0(2,2))/2

      vm(1,3) = v0(1,1)
      vm(2,3) = (v0(2,1)+v0(2,3))/2

      vm(1,4) = v0(1,2)
      vm(2,4) = (v0(2,1) + v0(2,3))/2

      vm(1,5) = (v0(1,1)+v0(1,2))/2
      vm(2,5) = v0(2,3)

c
cc     first quad
c
      v1(1,1) = v0(1,1)
      v1(2,1) = v0(2,1)

      v1(1,2) = vm(1,1)
      v1(2,2) = vm(2,1)

      v1(1,3) = vm(1,3)
      v1(2,3) = vm(2,3)
c
cc      second quad
c
      v2(1,1) = vm(1,1)
      v2(2,1) = vm(2,1)

      v2(1,2) = v0(1,2)
      v2(2,2) = v0(2,2)

      v2(1,3) = vm(1,2)
      v2(2,3) = vm(2,2)

c
cc      third quad
c
      v3(1,1) = vm(1,3)
      v3(2,1) = vm(2,3)

      v3(1,2) = vm(1,2)
      v3(2,2) = vm(2,2)

      v3(1,3) = v0(1,3)
      v3(2,3) = v0(2,3)

c
cc      fourth quad
c
      v4(1,1) = vm(1,2)
      v4(2,1) = vm(2,2)

      v4(1,2) = vm(1,4)
      v4(2,2) = vm(2,4)

      v4(1,3) = vm(1,5)
      v4(2,3) = vm(2,5)

      return
      end

c----------------------------------

      subroutine gen_xg_uniftree_nodes(nqorder,nnodes,nu,npts,qnodes,
     1   qwts)
c
c
c        this subroutine generates quadrature nodes
c        and weights on a nu \times nu uniform quads on [-1,1]^2 
c
c      input
c        nqorder - order of xiao-gimbutas nodes to be used
c        nnodes - number of xiao-gimbutas nodes of order nqorder
c        nu - number of quads in each direction
c        npts - total number of points (nu*nu*nnodes)
c      
c      output
c        qnodes(2,npts) - quadrature nodes
c        qwts(npts) - quadrature weights
c

      implicit real *8 (a-h,o-z)
      real *8 qnodes(2,npts),qwts(npts)
      real *8 qnodes0(2,nnodes),qwts0(nnodes)

      call squarearbq(nqorder,qnodes0,qwts0,nnodes)

      ra = 0
      do i=1,nnodes
        ra = ra + qwts0(i)
      enddo

      
      bs = 2.0d0/nu
      do iquad = 1,nu
        do jquad = 1,nu
          xc = -1 + (jquad-1)*bs + bs/2
          yc = -1 + (iquad-1)*bs + bs/2

          do i=1,nnodes
            ipt = ((iquad-1)*nu + jquad-1)*nnodes + i
            qnodes(1,ipt) = xc + qnodes0(1,i)*bs/2
            qnodes(2,ipt) = yc + qnodes0(2,i)*bs/2

            qwts(ipt) = qwts0(i)/4*bs*bs
          enddo
        enddo
      enddo


      return
      end

