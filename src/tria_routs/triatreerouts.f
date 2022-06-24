c
cc
cc
c
      subroutine gettritree(npatches,norder,npols,srccoefs,ntarg,
     1     xyztarg,itargptr,ntargptr,ntrimax,nlmax,rfac,ntri,nlev,
     2     ichild_start,da,tricm,trirad,tverts,itrirel,ier)
c
c
cc      this subroutine identifies the relevant quad tree on
c       the standard simplex relevant for all targets
c
c       and also identifies a bunch of other information along
c       the way
c
c       input arguments
c        npatches: number of patches
c        norder  : order of discretization
c        npols   : number of discretization 
c                   nodes on each patch
c        srccoefs(9,npols,npatches): coefficients
c                 of koornwinder expansion of xyz coordinates
c                 and derivative info
c        ntarg - total number of target points
c        xyztarg(3,ntargs) - location of target points 
c        itargptr(npatches) - xyztargs(:,itargptr(i):itargptr(i)+ntargptr(i))
c                       are the relevant set of targets for patch i
c        ntargptr(npatches) - number of relevant targets for patch i
c        ntrimax - max number of triangles expected
c        nlmax - max number of levels
c        rfac - distance criterion for deciding whether a triangle is in
c               in the far-field of a target or not
c               
c               target t, is considered in the far field of triangle T
c               if |t-c| > rfac*r
c               where c is the centroid of the triangle T, and r
c               is the radius of the enclosing sphere
c             
c        
c        output arguments
c         ntri - total number of triangles used.
c         nlev - total number of levels used
c         ichild_start - ichild_start(i) denotes the
c              starting index for the children of 
c              triangle i. If the triangle
c              is not subdivided then it stays at -1
c         da - denotes the 2*area of triangle i
c              
c         tricm(3,npatches,ntrimax) - tricm(:,i,j) 
c                          are the xyz coordinates of centroid
c                          of subtriangle j for patch i
c         trirad(npatches,ntrimax) - is the square of the
c                                  radius of smallest
c                                  enclosing sphere for subtriangle j 
c                                  for patch i
c         tverts(2,3,ntrimax) - uv coordinates of the three vertices
c                                of each subtriangle
c         itrirel(ntarg,ntrimax) - itrirel(j,i) = 1 if subtriangle
c                                    i will be used in the computation
c                                    at target j
c          ier  - error code
c                 ier  = 0, normal execution
c                 ier = 1, ntri>ntrimax, reallocate all arrays
c                          of appropriate sizes
c                                       
c
c        
      implicit none
      integer npatches,npols,norder
      real *8 srccoefs(9,npols,npatches)
      integer ntarg
      
      real *8 xyztarg(3,ntarg)
      integer itargptr(npatches),ntargptr(npatches)
      
      integer ntrimax,nlmax
      real *8 rfac,rfac2
      real *8 tricm(3,npatches,ntrimax)
      real *8 trirad(npatches,ntrimax)
      real *8 tverts(2,3,ntrimax)
      integer itrirel(ntarg,ntrimax)
      integer ichild_start(ntrimax)
      real *8 da(ntrimax)

      integer laddr(2,0:nlmax)

      integer ier
      
      integer ntri,nlev
      integer irefineall
      integer, allocatable :: irefinetri(:)


c
cc       temporary variables
c
      integer i,itarg,j,ii
      integer itri,ilev
      real *8 rr

      rfac2 = rfac**2

      ier = 0

      allocate(irefinetri(ntrimax))

      laddr(1,0) = 1
      laddr(2,0) = 1
      nlev = 0

      ntri = 1
      da(1) = 1.0d0
      tverts(1,1,1) = 0
      tverts(2,1,1) = 0

      tverts(1,2,1) = 1
      tverts(2,2,1) = 0
      
      tverts(1,3,1) = 0
      tverts(2,3,1) = 1


      call gettricms(npatches,norder,npols,srccoefs,tverts(1,1,1),
     1      tricm(1,1,1),trirad(1,1))



      irefineall = 0
      irefinetri(1) = 0


      do i=1,npatches
         do itarg = itargptr(i),itargptr(i)+ntargptr(i)-1
           itrirel(itarg,1) = 1
c
cc          compute distance between targ, and relevant tricm
c
           rr = (xyztarg(1,itarg)-tricm(1,i,1))**2 +
     1          (xyztarg(2,itarg)-tricm(2,i,1))**2 + 
     2          (xyztarg(3,itarg)-tricm(3,i,1))**2

           if(rr.lt.rfac2*trirad(i,1)) then
             itrirel(itarg,1) = -1 
             irefinetri(1) = 1
             irefineall = 1
           endif
         enddo
      enddo


      do ilev=1,nlmax
        if(irefineall.eq.1) then
c
c          loop over triangles at the previous level
c

          irefineall = 0
          laddr(1,ilev) = ntri+1
          do itri=laddr(1,ilev-1),laddr(2,ilev-1)
c
c            determine if current triangle needs to be refined
c
             if(irefinetri(itri).eq.1) then

c
cc                check if total number of triangles
c                 exceeds max triangles
                if(ntri+4.gt.ntrimax) then
                  ier = 1
                  return
                endif
                ichild_start(itri) = ntri+1

c
c               get vertices of the children triangle
c

                call gettrichildren(tverts(1,1,itri),
     1            tverts(1,1,ntri+1),tverts(1,1,ntri+2),
     2            tverts(1,1,ntri+3),tverts(1,1,ntri+4))
                

c
c               get cms of children triangle
c
                rr = da(itri)*0.25d0
                do i=1,4
                  call gettricms(npatches,norder,npols,srccoefs,
     1             tverts(1,1,ntri+i),tricm(1,1,ntri+i),
     2             trirad(1,ntri+i))
                  irefinetri(ntri+i) = 0
                  da(ntri+i) = rr
                enddo

c
c               decide if the children triangle are relevant for 
c               a target or not, and if children triangles need
c               to be further refined
c
                do i=1,npatches
                  do itarg = itargptr(i),itargptr(i)+ntargptr(i)-1
                    if(itrirel(itarg,itri).eq.-1) then
                      do j=1,4
                        itrirel(itarg,ntri+j) = 1
c
cc          compute distance between targ, and relevant tricm
c
                        rr = (xyztarg(1,itarg)-tricm(1,i,ntri+j))**2 +
     1                       (xyztarg(2,itarg)-tricm(2,i,ntri+j))**2 + 
     2                       (xyztarg(3,itarg)-tricm(3,i,ntri+j))**2
                        if(rr.lt.rfac2*trirad(i,ntri+j)) then
                          itrirel(itarg,ntri+j) = -1 
                          irefinetri(ntri+j) = 1
                          irefineall = 1
                        endif
                      enddo
                    endif
                  enddo
                enddo

c
c               update number of triangles
               ntri = ntri + 4
             endif
          enddo
          laddr(2,ilev) = ntri
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
      subroutine gettricms(npatches,norder,npols,srccoefs,tverts,tricm,
     1   trirad)
c
cc     this subroutine computes the centroid and radius of enclosing
c      sphere of the mapped triangle with vertices tverts
c      through the koornwinder expansion srccoefs
c
      implicit none
      integer npatches,npols,norder
      real *8 srccoefs(9,npols,npatches)
      real *8 tverts(2,3)
      real *8 tricm(3,npatches)
      real *8 trirad(npatches)
      real *8 vs(3,3)
      real *8, allocatable :: pols(:,:)
      real *8 rr

      integer i,j,k

      allocate(pols(npols,3))

      call koorn_pols(tverts(1,1),norder,npols,pols(1,1)) 
      call koorn_pols(tverts(1,2),norder,npols,pols(1,2)) 
      call koorn_pols(tverts(1,3),norder,npols,pols(1,3))


      do i=1,npatches
        do j=1,3
          do k=1,3
            vs(j,k) = 0
          enddo
        enddo
        do j=1,3
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
          tricm(j,i) = 0
          do k=1,3
            tricm(j,i) = tricm(j,i) + vs(j,k)
          enddo
          tricm(j,i) = tricm(j,i)/3
        enddo

c
cc         compute the radius
c
        trirad(i) = 0
        do j=1,3
          rr = (tricm(1,i) - vs(1,j))**2 + 
     1         (tricm(2,i) - vs(2,j))**2 +
     2          (tricm(3,i) - vs(3,j))**2
          if(rr.gt.trirad(i)) trirad(i) = rr
        enddo
      enddo

      return
      end
c
c
c
c
c

      subroutine gettrichildren(v0,v1,v2,v3,v4)
c  
cc       given the three vertices of a triangle v0,
c        this subroutine returns the vertices of 4 
c        smaller triangles constructed using the
c        midpoints of the triangles
c 
c        input:
c        v0 - real *8 (2,3)  
c              vertices of parent triangle
c
c        output:
c        v1,v2,v3,v4 - real *8 (2,3)
c                 vertices of children triangle
      
      implicit real *8 (a-h,o-z)
      real *8 v0(2,3),v1(2,3),v2(2,3),v3(2,3),v4(2,3),vm(2,3)

      vm(1,1) = (v0(1,1)+v0(1,2))/2
      vm(2,1) = (v0(2,1)+v0(2,2))/2

      vm(1,2) = (v0(1,3)+v0(1,2))/2
      vm(2,2) = (v0(2,3)+v0(2,2))/2

      vm(1,3) = (v0(1,1)+v0(1,3))/2
      vm(2,3) = (v0(2,1)+v0(2,3))/2

c
cc     first triangle
c
      v1(1,1) = v0(1,1)
      v1(2,1) = v0(2,1)

      v1(1,2) = vm(1,1)
      v1(2,2) = vm(2,1)

      v1(1,3) = vm(1,3)
      v1(2,3) = vm(2,3)
c
cc      second triangle
c
      v2(1,1) = vm(1,1)
      v2(2,1) = vm(2,1)

      v2(1,2) = v0(1,2)
      v2(2,2) = v0(2,2)

      v2(1,3) = vm(1,2)
      v2(2,3) = vm(2,2)

c
cc      third triangle
c
      v3(1,1) = vm(1,3)
      v3(2,1) = vm(2,3)

      v3(1,2) = vm(1,2)
      v3(2,2) = vm(2,2)

      v3(1,3) = v0(1,3)
      v3(2,3) = v0(2,3)

c
cc      fourth triangle
c
      v4(1,1) = vm(1,2)
      v4(2,1) = vm(2,2)

      v4(1,2) = vm(1,3)
      v4(2,2) = vm(2,3)

      v4(1,3) = vm(1,1)
      v4(2,3) = vm(2,1)

      return
      end

c----------------------------------

c
c
c
      subroutine gen_xg_unif_nodes_tri(nlev,nqorder,nnodes,npts,qnodes,
     1   qwts)
      implicit real *8 (a-h,o-z)
      real *8, allocatable :: tvs(:,:,:)
      real *8, allocatable :: qnodes0(:,:),qwts0(:)
      real *8 qnodes(2,npts),qwts(npts)

      

      ntri = (4**(nlev+1)-1)/3

      allocate(tvs(2,3,ntri))

      tvs(1,1,1) = 0
      tvs(2,1,1) = 0

      tvs(1,2,1) = 1
      tvs(2,2,1) = 0

      tvs(1,3,1) = 0
      tvs(2,3,1) = 1

      allocate(qnodes0(2,nnodes),qwts0(nnodes))
      call triasymq(nqorder,tvs(1,1,1),tvs(1,2,1),tvs(1,3,1),
     1       qnodes0,qwts0,nnodes0)

      do i=0,nlev-1
        istart = (4**(i)-1)/3 + 1
        nb = 4**i
        iend = istart + nb-1

        do itrip = istart,iend
          itric1 = (itrip-istart)*4+iend
c
cc   compute the area element and the location of vertices
c    of the children
c
          call gettrichildren(tvs(1,1,itrip),tvs(1,1,itric1+1), 
     1       tvs(1,1,itric1+2),tvs(1,1,itric1+3), 
     2       tvs(1,1,itric1+4))   
        enddo
      enddo


      istart = ntri-4**nlev+1
      iend = ntri


      do itri=istart,iend
        instart = (itri-istart)*nnodes + 1
        call mapuv_tri(tvs(1,1,itri),nnodes,qnodes0,qnodes(1,instart))
        do i=1,nnodes
          qwts(instart+i-1) = qwts0(i)/4**nlev
        enddo
      enddo

      return
      end
c
c
c
c
      subroutine writetritree(iw,ntri,tverts,iprint)
      implicit none
      integer ntri,iw,iprint(ntri)
      real *8 tverts(2,3,ntri)

      integer i,j,k,itri

 1100 format(6(2x,e11.5))
      do itri=1,ntri
        if(iprint(itri).ne.0) then
           write(iw,1100) tverts(1,1,itri),tverts(1,2,itri),
     1              tverts(1,3,itri),tverts(2,1,itri),
     2              tverts(2,2,itri),tverts(2,3,itri)
        endif
      enddo

      return
      end

c
c
c
c
c
c
c--------------------------------------------------------------------------------
        
      subroutine mapuv_tri(verts,kpols,uvs,uvout)
      implicit real *8 (a-h,o-z)
      integer kpols
      real *8 verts(2,3),uvs(2,kpols),uvout(2,kpols)

      dx = verts(1,2)-verts(1,1)
      dy = verts(2,3)-verts(2,1) 

      do i=1,kpols
        uvout(1,i) = verts(1,1) + dx*uvs(1,i)
        uvout(2,i) = verts(2,1) + dy*uvs(2,i)
      enddo

      return
      end
c-----------------------------------------      
        subroutine get_norms_qwts_tri(kpols,whts,srcvals,da,
     1       qwts)
        implicit real *8 (a-h,o-z)
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

