
c
cc
cc
c
      subroutine getquadtree(npatches,norder,ipoly,ttype,npols,srccoefs,
     1  ntarg,xyztarg,itargptr,ntargptr,nquadmax,nlmax,rfac,nquad,nlev,
     2  ichild_start,da,quadcm,quadrad,qverts,iquadrel,ier)
c-------------
c  Task:
c    Construct a tree of quads on a high-order quadrangle patch,
c    such that each, for every target, there is decomposition
c    of the quadrangle patch, such that each subpatch
c    is well-separated from the target by a
c    user-prescribed factor
c
c  Algorithm:
c    densely checking whether each quad needs further subdivision.
c    This algorithm is supposed to be used when the number of targets
c    per patch is O(1).
c
c  Input arguments:
c    - npatches: integer *8
c        number of patches
c    - norder: integer *8
c        discretization order of patches
c    - ipoly: integer *8
c        Type of polynomial to be used
c        * ipoly = 0, Legendre polynomials
c        * ipoly = 1, Chebyshev polynomials
c    - ttype: character
c        type of number of polynomials used
c        * ttype = 'F' or 'f', full degree polynomials
c        * ttype = 'T' or 't', total degree polynomials
c    - npols: integer *8
c        Number of polynomials 
c    - srccoefs: real *8 (9,npols,npatches)
c        coefficients of basis expansion of xyz coordinates,
c        dxyz/du, and dxyz/dv
c    - ntarg: integer *8
c        number of targets
c    - xyztarg: real *8 (ndtarg,ntarg)
c        target information
c    - itargptr: integer *8(npatches)
c        xyztargs(:,itargptr(i):itargptr(i)+ntargptr(i)-1)
c        are the relevant set of targets for patch i
c    - ntargptr: integer *8(npatches)
c        ntargptr(i) is the number of relevant targets for patch i
c    - nquadmax: integer *8
c        max number of quadrangles supported in the heirarchy
c        if total number of quadrangles required exceeds
c        nquadmax, subroutine exits with error code
c        ier = 4
c    - nlmax: integer *8
c        max depth of tree
c    - rfac: real *8
c        parameter for defining refinement criterion, quadrangles
c        refined if not separated from target by
c        rfac*r_{t} where r_{t} is radius of enclosing sphere
c        in quadrangle
c  Output arguments:
c    - nquad: integer *8
c        number of quadrangles used
c    - nlev: integer *8
c        depth of tree
c    - ichild_start: integer *8(nquadmax)
c        location in set of quads where children for quad i start
c    - da: real *8(nquadmax)
c        da(i) = denotes 2*area of quadrangle i
c    - quadcm: real *8(3,npatches,nquadmax)
c        quadcm(:,i,j) are the xyz coordinates of centroid
c             of subquadrangle j for patch i
c    - quadrad: real *8(npatches,nquadmax)
c        is the square of the radius of smallest enclosing sphere 
c        for subquadrangle j for patch i
c    - qverts: real *8(2,3,nquadmax)
c        uv coordinates of the three vertices of each subquadrangle
c    - iquadrel: integer *8(ntarg,nquadmax)
c        iquadrel(j,i) = 1 if subquadrangle i will be used in 
c          the computation at target j
c    - ier: error code
c        * ier = 0, normal execution
c        * ier = 4, nquad>nquadmax, reallocate all arrays of 
c            appropriate sizes
c---------
c        
      implicit none
      integer *8, intent(in) :: npatches,npols,norder
      real *8, intent(in) :: srccoefs(9,npols,npatches)
      integer *8, intent(in) :: ntarg
      integer *8, intent(in) :: ipoly
      character *1, intent(in) :: ttype
      
      real *8, intent(in) :: xyztarg(3,ntarg)
      integer *8, intent(in) :: itargptr(npatches),ntargptr(npatches)
      
      integer *8, intent(in) :: nquadmax,nlmax
      real *8, intent(in) :: rfac
      real *8, intent(out) :: quadcm(3,npatches,nquadmax)
      real *8, intent(out) :: quadrad(npatches,nquadmax)
      real *8, intent(out) :: qverts(2,3,nquadmax)
      integer *8, intent(out) :: iquadrel(ntarg,nquadmax)
      integer *8, intent(out) :: ichild_start(nquadmax)
      real *8, intent(out) :: da(nquadmax)
      integer *8, intent(out) :: ier

      real *8 rfac2
      integer *8 laddr(2,0:nlmax)
      
      integer *8 nquad,nlev
      integer *8 irefineall
      integer *8, allocatable :: irefinequad(:)


c
cc       temporary variables
c
      integer *8 i,itarg,j,ii
      integer *8 iquad,ilev
      real *8 rr

      rfac2 = rfac**2

      ier = 0

      allocate(irefinequad(nquadmax))

      laddr(1,0) = 1
      laddr(2,0) = 1
      nlev = 0

      nquad = 1
c
c  check whether da needs to be rescaled
c
      da(1) = 1.0d0


      qverts(1,1,1) = -1
      qverts(2,1,1) = -1

      qverts(1,2,1) = 1
      qverts(2,2,1) = -1
      
      qverts(1,3,1) = -1
      qverts(2,3,1) = 1


      call getquadcms(npatches,norder,ipoly,ttype,npols,srccoefs,
     1  qverts(1,1,1),quadcm(1,1,1),quadrad(1,1))



      irefineall = 0
      irefinequad(1) = 0


      do i=1,npatches
         do itarg = itargptr(i),itargptr(i)+ntargptr(i)-1
           iquadrel(itarg,1) = 1
c
cc          compute distance between targ, and relevant quadcm
c
           rr = (xyztarg(1,itarg)-quadcm(1,i,1))**2 +
     1          (xyztarg(2,itarg)-quadcm(2,i,1))**2 + 
     2          (xyztarg(3,itarg)-quadcm(3,i,1))**2

           if(rr.lt.rfac2*quadrad(i,1)) then
             iquadrel(itarg,1) = -1 
             irefinequad(1) = 1
             irefineall = 1
           endif
         enddo
      enddo


      do ilev=1,nlmax
        if(irefineall.eq.1) then
c
c          loop over quadrangles at the previous level
c

          irefineall = 0
          laddr(1,ilev) = nquad+1
          do iquad=laddr(1,ilev-1),laddr(2,ilev-1)
c
c            determine if current quadrangle needs to be refined
c
             if(irefinequad(iquad).eq.1) then

c
cc                check if total number of quadrangles
c                 exceeds max quadrangles
                if(nquad+4.gt.nquadmax) then
                  ier = 4
                  return
                endif
                ichild_start(iquad) = nquad+1

c
c               get vertices of the children quadrangle
c

                call getquadchildren(qverts(1,1,iquad),
     1            qverts(1,1,nquad+1),qverts(1,1,nquad+2),
     2            qverts(1,1,nquad+3),qverts(1,1,nquad+4))
                

c
c               get cms of children quadrangle
c
                rr = da(iquad)*0.25d0
                do i=1,4
                  call getquadcms(npatches,norder,ipoly,ttype,npols,
     1              srccoefs,qverts(1,1,nquad+i),quadcm(1,1,nquad+i),
     2              quadrad(1,nquad+i))
                  irefinequad(nquad+i) = 0
                  da(nquad+i) = rr
                enddo

c
c               decide if the children quadrangle are relevant for 
c               a target or not, and if children quadrangles need
c               to be further refined
c
                do i=1,npatches
                  do itarg = itargptr(i),itargptr(i)+ntargptr(i)-1
                    if(iquadrel(itarg,iquad).eq.-1) then
                      do j=1,4
                        iquadrel(itarg,nquad+j) = 1
c
cc          compute distance between targ, and relevant quadcm
c
                        rr = (xyztarg(1,itarg)-quadcm(1,i,nquad+j))**2 +
     1                       (xyztarg(2,itarg)-quadcm(2,i,nquad+j))**2 +
     2                       (xyztarg(3,itarg)-quadcm(3,i,nquad+j))**2
                        if(rr.lt.rfac2*quadrad(i,nquad+j)) then
                          iquadrel(itarg,nquad+j) = -1 
                          irefinequad(nquad+j) = 1
                          irefineall = 1
                        endif
                      enddo
                    endif
                  enddo
                enddo

c
c               update number of quadrangles
               nquad = nquad + 4
             endif
          enddo
          laddr(2,ilev) = nquad
          nlev = nlev+1

        else
          exit
        endif
      enddo

      return
      end
c
c      
c
c
c
      subroutine getquadcms(npatches,norder,ipoly,ttype,npols,srccoefs,
     1   qverts,quadcm,quadrad)
c
c-------------
c  Task:
c    Compute the centroid and radius of enclosing sphere of mapped
c    quadrangle with vertices qverts through their basis function
c    expansions
c
c  Input arguments:
c    - npatches: integer *8
c        number of patches
c    - norder: integer *8
c        discretization order of patches
c    - ipoly: integer *8
c        Type of polynomial to be used
c        * ipoly = 0, Legendre polynomials
c        * ipoly = 1, Chebyshev polynomials
c    - ttype: character
c        type of number of polynomials used
c        * ttype = 'F' or 'f', full degree polynomials
c        * ttype = 'T' or 't', total degree polynomials
c    - npols: integer *8
c        Number of polynomials 
c    - srccoefs: real *8 (9,npols,npatches)
c        coefficients of basis expansion of xyz coordinates,
c        dxyz/du, and dxyz/dv
c    - qverts: real *8(2,3,nquadmax)
c        uv coordinates of the three vertices of each subquadrangle
c  Output arguments:
c    - quadcm: real *8(3,npatches)
c        quadcm(:,i) are the xyz coordinates of centroid
c             of patch i
c    - quadrad: real *8(npatches)
c        is the square of the radius of smallest enclosing sphere 
c        for patch i
c-----------------
c
      implicit none
      integer *8 npatches,npols,norder,ipoly
      character *1 ttype
      real *8 srccoefs(9,npols,npatches)
      real *8 qverts(2,3)
      real *8 quadcm(3,npatches)
      real *8 quadrad(npatches)
      real *8 vs(3,4)
      real *8, allocatable :: pols(:,:)
      real *8 qvuse(2,4)
      real *8 rr

      integer *8 i,j,k

      allocate(pols(npols,4))
      qvuse(1:2,1:3) = qverts
      qvuse(1,4) = qverts(1,2)
      qvuse(2,4) = qverts(2,3)

      call polytens_pols_2d(ipoly,qvuse(1,1),norder,ttype,pols(1,1))
      call polytens_pols_2d(ipoly,qvuse(1,2),norder,ttype,pols(1,2))
      call polytens_pols_2d(ipoly,qvuse(1,3),norder,ttype,pols(1,3))
      call polytens_pols_2d(ipoly,qvuse(1,4),norder,ttype,pols(1,4))


      do i=1,npatches
        do j=1,3
          do k=1,4
            vs(j,k) = 0
          enddo
        enddo
        do j=1,4
          do k=1,npols
            vs(1,j) = vs(1,j)+srccoefs(1,k,i)*pols(k,j)
            vs(2,j) = vs(2,j)+srccoefs(2,k,i)*pols(k,j)
            vs(3,j) = vs(3,j)+srccoefs(3,k,i)*pols(k,j)
          enddo
        enddo
c
cc        compute the centroid
c
        do j=1,3
          quadcm(j,i) = 0
          do k=1,4
            quadcm(j,i) = quadcm(j,i) + vs(j,k)
          enddo
          quadcm(j,i) = quadcm(j,i)/4
        enddo

c
cc         compute the radius
c
        quadrad(i) = 0
        do j=1,4
          rr = (quadcm(1,i) - vs(1,j))**2 + 
     1         (quadcm(2,i) - vs(2,j))**2 +
     2          (quadcm(3,i) - vs(3,j))**2
          if(rr.gt.quadrad(i)) quadrad(i) = rr
        enddo
      enddo

      return
      end
c
c
c
c
c
c

      subroutine getquadchildren(v0,v1,v2,v3,v4)
c-----------
c  Task:
c    Given the three vertices of a quad v0, compute
c    the vertices of the 4 smaller quads constructed
c    using the midpoints of the quads
c  
c  Input:
c    v0: real *8 (2,3)
c      vertices of the quad to be subdivided
c  Output:
c    v1,v2,v3,v4: real *8 (2,3)
c      vertices of children quads
c----------------      
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
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
c
c
c
c
      subroutine gen_xg_unif_nodes_quad(nlev,nqorder,nnodes,npts,qnodes,
     1   qwts)
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      real *8, allocatable :: tvs(:,:,:)
      real *8, allocatable :: qnodes0(:,:),qwts0(:)
      real *8 qnodes(2,npts),qwts(npts)

      

      nquad = (4**(nlev+1)-1)/3

      allocate(tvs(2,3,nquad))

      tvs(1,1,1) = -1
      tvs(2,1,1) = -1

      tvs(1,2,1) = 1
      tvs(2,2,1) = -1

      tvs(1,3,1) = -1
      tvs(2,3,1) = 1

      allocate(qnodes0(2,nnodes),qwts0(nnodes))
      nnodes0 = nnodes
      call squarearbq(nqorder,qnodes0,qwts0,nnodes0)

      do i=0,nlev-1
        istart = (4**(i)-1)/3 + 1
        nb = 4**i
        iend = istart + nb-1

        do iquadp = istart,iend
          iquadc1 = (iquadp-istart)*4+iend
c
cc   compute the area element and the location of vertices
c    of the children
c
          call getquadchildren(tvs(1,1,iquadp),tvs(1,1,iquadc1+1), 
     1       tvs(1,1,iquadc1+2),tvs(1,1,iquadc1+3), 
     2       tvs(1,1,iquadc1+4))   
        enddo
      enddo


      istart = nquad-4**nlev+1
      iend = nquad


      do iquad=istart,iend
        instart = (iquad-istart)*nnodes + 1
        call mapuv_quad(tvs(1,1,iquad),nnodes,qnodes0,qnodes(1,instart))
        do i=1,nnodes
          qwts(instart+i-1) = qwts0(i)/4**nlev
        enddo
      enddo

      return
      end
c
c
c
c--------------------------------------------------------------------------------
        
      subroutine mapuv_quad(verts,kpols,uvs,uvout)
c--------
c  Task:
c    Map discretization nodes on [-1,1]^2 to an arbitrary
c    quadrangle aligned along the coordinate axes
c  
c  Input arguments:
c    - verts: real *8 (2,3)
c        Vertices defining the sub-quadrangle
c    - kpols: integer *8
c        number of input discretization points
c    - uvs: real *8 (2,kpols)
c        discretization points on [-1,1]^2
c  Output argument:
c    - uvout: real *8 (2,kpols)
c        discretization points on the mapped quadrangle
c------------

      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8 kpols
      real *8 verts(2,3),uvs(2,kpols),uvout(2,kpols)

      dx = verts(1,2)-verts(1,1)
      dy = verts(2,3)-verts(2,1) 

      do i=1,kpols
        uvout(1,i) = verts(1,1) + dx*(uvs(1,i)+1)/2
        uvout(2,i) = verts(2,1) + dy*(uvs(2,i)+1)/2
      enddo

      return
      end
c-----------------------------------------      
        subroutine get_norms_qwts_quad(kpols,whts,srcvals,da,
     1       qwts)
c----------------
c  Task: 
c    Compute normals and quadrature weights on a quadrangle
c  
c  Input arguments:
c    - kpols: integer *8
c        number of points on the quadrangle
c    - whts: real *8 (kpols)
c        quadrature weights for integrating smooth functions 
c        on [-1,1]^2
c    - srcvals: real *8 (12,kpols)
c        xyz coorindates, dxyz/du, dxyz/dv and normals
c        at the discretization nodes on the patch
c    - da: real *8
c        surface area of the subquad
c   
c  Output arguments:
c    - qwts: real *8 (kpols)
c        quadrature weights for integrating smooth functions
c        on the patch
c-------------
c
        implicit real *8 (a-h,o-z)
        implicit integer *8 (i-n)
        real *8 srcvals(12,kpols),qwts(kpols),whts(kpols)
        real *8 tmp(3),da


        do i=1,kpols
          call cross_prod3d(srcvals(4,i),srcvals(7,i),srcvals(10,i))
          rr = sqrt(srcvals(10,i)**2 + srcvals(11,i)**2 + 
     1        srcvals(12,i)**2)
          qwts(i) = rr*da*whts(i)
          srcvals(10,i) = srcvals(10,i)/rr
          srcvals(11,i) = srcvals(11,i)/rr
          srcvals(12,i) = srcvals(12,i)/rr
        enddo

        return
        end
c
c
c
c
c-------------------------------------------------

