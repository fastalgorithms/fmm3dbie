      subroutine get_fmm_thresh(nds,ns,src,ndt,nt,trg,thresh)
c---------------------------------------
c  This subroutine returns the threshold used by an fmm
c  for ignoring self interactions for a given collection of sources
c  and targets.
c
c  The rotuine identifies the size of the smallest bounding cube
c  and sets the threshold to be 2**(-51)*(length of bounding cube edge)
c
c
c  Input arguments:
c  
c    - nds: integer(8)
c        leading dimension for sources array (must be at least 3)
c    - ns: integer(8)
c        number of sources
c    - src: real *8 (nds,ns)
c        source info array (src(1:3,i) should contain the x,y,z components)
c    - ndt: integer(8)
c        leading dimension for targets array (must be at least 3)
c    - nt: integer(8)
c        number of targets
c    - trg: real *8 (ndt,nt)c
c        targ info array (targ(1:3,i) should contain the x,y,z components)
c
c  Output arguments:
c  
c    - thresh: real *8
c        fmm threshold
c
      implicit real *8 (a-h,o-z)
      integer(8), intent(in) :: nds,ns,ndt,nt
      real *8, intent(in) :: src(nds,ns),trg(ndt,nt)
      real *8, intent(out) :: thresh

      integer(8) i
      real *8 xmin,xmax,ymin,ymax,zmin,zmax,bsizex,bsizey,bsizez,bsize

      xmin = src(1,1)
      xmax = src(1,1)
      ymin = src(2,1)
      ymax = src(2,1)
      zmin = src(3,1)
      zmax = src(3,1)
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$REDUCTION(min:xmin,ymin,zmin)
C$OMP$REDUCTION(max:xmax,ymax,zmax)
C$OMP$PRIVATE(i)
      do i=1,ns
        if(src(1,i) .lt. xmin) xmin=src(1,i)
        if(src(1,i) .gt. xmax) xmax=src(1,i)
        if(src(2,i) .lt. ymin) ymin=src(2,i)
        if(src(2,i) .gt. ymax) ymax=src(2,i)
        if(src(3,i) .lt. zmin) zmin=src(3,i)
        if(src(3,i) .gt. zmax) zmax=src(3,i)
      enddo
C$OMP END PARALLEL DO    

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$REDUCTION(min:xmin,ymin,zmin)
C$OMP$REDUCTION(max:xmax,ymax,zmax)
C$OMP$PRIVATE(i)
      do i=1,nt
        if(trg(1,i) .lt. xmin) xmin=trg(1,i)
        if(trg(1,i) .gt. xmax) xmax=trg(1,i)
        if(trg(2,i) .lt. ymin) ymin=trg(2,i)
        if(trg(2,i) .gt. ymax) ymax=trg(2,i)
        if(trg(3,i) .lt. zmin) zmin=trg(3,i)
        if(trg(3,i) .gt. zmax) zmax=trg(3,i)
      enddo
C$OMP END PARALLEL DO    

      bsize = xmax-xmin
      bsizey = ymax-ymin
      bsizez = zmax-zmin
      if(bsizey.gt.bsize) bsize = bsizey
      if(bsizez.gt.bsize) bsize = bsizez

      thresh = 2.0d0**(-51)*bsize

      return
      end
c
c
c
c
c
c
c
      subroutine get_bsize(nds,ns,src,ndt,nt,trg,bsize)
c---------------------------------------
c  This subroutine returns the size of the smallest bounding cube 
c  for a given collection of sources and targets.
c
c
c  Input arguments:
c  
c    - nds: integer(8)
c        leading dimension for sources array (must be at least 3)
c    - ns: integer(8)
c        number of sources
c    - src: real *8 (nds,ns)
c        source info array (src(1:3,i) should contain the x,y,z components)
c    - ndt: integer(8)
c        leading dimension for targets array (must be at least 3)
c    - nt: integer(8)
c        number of targets
c    - trg: real *8 (ndt,nt)
c        targ info array (targ(1:3,i) should contain the x,y,z components)
c
c  Output arguments:
c  
c    - bsize: real *8
c        edge length of smallest bounding cube
c
      implicit real *8 (a-h,o-z)
      integer(8), intent(in) :: nds,ns,ndt,nt
      real *8, intent(in) :: src(nds,ns),trg(ndt,nt)
      real *8, intent(out) :: bsize

      integer(8) i
      real *8 xmin,xmax,ymin,ymax,zmin,zmax,bsizex,bsizey,bsizez

      xmin = src(1,1)
      xmax = src(1,1)
      ymin = src(2,1)
      ymax = src(2,1)
      zmin = src(3,1)
      zmax = src(3,1)
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$REDUCTION(min:xmin,ymin,zmin)
C$OMP$REDUCTION(max:xmax,ymax,zmax)
C$OMP$PRIVATE(i)
      do i=1,ns
        if(src(1,i) .lt. xmin) xmin=src(1,i)
        if(src(1,i) .gt. xmax) xmax=src(1,i)
        if(src(2,i) .lt. ymin) ymin=src(2,i)
        if(src(2,i) .gt. ymax) ymax=src(2,i)
        if(src(3,i) .lt. zmin) zmin=src(3,i)
        if(src(3,i) .gt. zmax) zmax=src(3,i)
      enddo
C$OMP END PARALLEL DO    

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$REDUCTION(min:xmin,ymin,zmin)
C$OMP$REDUCTION(max:xmax,ymax,zmax)
C$OMP$PRIVATE(i)
      do i=1,nt
        if(trg(1,i) .lt. xmin) xmin=trg(1,i)
        if(trg(1,i) .gt. xmax) xmax=trg(1,i)
        if(trg(2,i) .lt. ymin) ymin=trg(2,i)
        if(trg(2,i) .gt. ymax) ymax=trg(2,i)
        if(trg(3,i) .lt. zmin) zmin=trg(3,i)
        if(trg(3,i) .gt. zmax) zmax=trg(3,i)
      enddo
C$OMP END PARALLEL DO    

      bsize = xmax-xmin
      bsizey = ymax-ymin
      bsizez = zmax-zmin
      if(bsizey.gt.bsize) bsize = bsizey
      if(bsizez.gt.bsize) bsize = bsizez


      return
      end
c
c
c
c
c
c
c
c
      subroutine get_bsize_rect(nds,ns,src,ndt,nt,trg,bsizex,bsizey,
     1   bsizez)
c---------------------------------------
c  This subroutine returns the size of the smallest bounding rect
c  parallelopiped
c  for a given collection of sources and targets.
c
c
c  Input arguments:
c  
c    - nds: integer(8)
c        leading dimension for sources array (must be at least 3)
c    - ns: integer(8)
c        number of sources
c    - src: real *8 (nds,ns)
c        source info array (src(1:3,i) should contain the x,y,z components)
c    - ndt: integer(8)
c        leading dimension for targets array (must be at least 3)
c    - nt: integer(8)
c        number of targets
c    - trg: real *8 (ndt,nt)
c        targ info array (targ(1:3,i) should contain the x,y,z components)
c
c  Output arguments:
c  
c    - bsize: real *8
c        edge length of smallest bounding cube
c
      implicit real *8 (a-h,o-z)
      integer(8), intent(in) :: nds,ns,ndt,nt
      real *8, intent(in) :: src(nds,ns),trg(ndt,nt)
      real *8, intent(out) :: bsizex,bsizey,bsizez

      integer(8) i
      real *8 xmin,xmax,ymin,ymax,zmin,zmax

      xmin = src(1,1)
      xmax = src(1,1)
      ymin = src(2,1)
      ymax = src(2,1)
      zmin = src(3,1)
      zmax = src(3,1)
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$REDUCTION(min:xmin,ymin,zmin)
C$OMP$REDUCTION(max:xmax,ymax,zmax)
C$OMP$PRIVATE(i)
      do i=1,ns
        if(src(1,i) .lt. xmin) xmin=src(1,i)
        if(src(1,i) .gt. xmax) xmax=src(1,i)
        if(src(2,i) .lt. ymin) ymin=src(2,i)
        if(src(2,i) .gt. ymax) ymax=src(2,i)
        if(src(3,i) .lt. zmin) zmin=src(3,i)
        if(src(3,i) .gt. zmax) zmax=src(3,i)
      enddo
C$OMP END PARALLEL DO    

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$REDUCTION(min:xmin,ymin,zmin)
C$OMP$REDUCTION(max:xmax,ymax,zmax)
C$OMP$PRIVATE(i)
      do i=1,nt
        if(trg(1,i) .lt. xmin) xmin=trg(1,i)
        if(trg(1,i) .gt. xmax) xmax=trg(1,i)
        if(trg(2,i) .lt. ymin) ymin=trg(2,i)
        if(trg(2,i) .gt. ymax) ymax=trg(2,i)
        if(trg(3,i) .lt. zmin) zmin=trg(3,i)
        if(trg(3,i) .gt. zmax) zmax=trg(3,i)
      enddo
C$OMP END PARALLEL DO    

      bsizex = xmax-xmin
      bsizey = ymax-ymin
      bsizez = zmax-zmin


      return
      end
c
c
c
c
c
c
c
