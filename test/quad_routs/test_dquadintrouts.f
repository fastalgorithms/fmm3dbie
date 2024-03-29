
      subroutine dquadintrouts_testing(istrat,ifp,intype,ipoly,ttype,
     1   isuccess)
      implicit real *8 (a-h,o-z)
      real *8 xyztarg(3,3),xyztarg0(3,3)
      real *8, allocatable :: srcvals(:,:,:)
      real *8, allocatable :: srccoefs(:,:,:)
      real *8, allocatable :: umatr(:,:),vmatr(:,:),uvs(:,:),wts(:)
      real *8, allocatable :: slp(:,:)
      real *8, allocatable :: slp_ex(:,:)
      integer ipars(10)
      integer, allocatable :: iprint(:,:)
      integer itargptr(2),ntargptr(2)
      character *1 ttype


      real *8 dpars(5)
      complex *16 zpars(3),ima

      external lslp

      character type

      ima = dcmplx(0.0d0,1.0d0)
      
      done = 1
      pi = atan(done)*4


      call prini(6,13)


      norder = 8
      npatches = 2


      npts = (norder+1)*(norder+1)
      if(ttype.eq.'t'.or.ttype.eq.'T') then
        npols = (norder+1)*(norder+2)/2
      else if(ttype.eq.'f'.or.ttype.eq.'F') then
        npols = npts
      else
        ier = 2
        call prinf('invalid polynomial type*',i,0)
        isuccess = 0
        return
      endif
      
      ntarg = 3
      allocate(umatr(npols,npts),vmatr(npts,npols))
      allocate(uvs(2,npts),wts(npts))
      print *, "norder=",norder
      print *, "npts=",npts
      print *, "npols=",npols
      itype = 2
      call polytens_exps_nd(2,ipoly,itype,norder+1,ttype,uvs,
     1     umatr,npols,vmatr,npts,wts)
      thet = pi/4
      rr = 1.0d-1

      allocate(srccoefs(9,npols,npatches),srcvals(9,npols,npatches))
      do i=1,npols
        srcvals(1,i,1) = uvs(1,i)*rr
        srcvals(2,i,1) = uvs(2,i)*rr
        srcvals(3,i,1) = 0
        
        srcvals(4,i,1) = rr
        srcvals(5,i,1) = 0
        srcvals(6,i,1) = 0

        srcvals(7,i,1) = 0
        srcvals(8,i,1) = rr
        srcvals(9,i,1) = 0

        srcvals(1,i,2) = cos(thet)*uvs(1,i) + sin(thet)*uvs(2,i)+2.1d0
        srcvals(2,i,2) = -sin(thet)*uvs(1,i) + cos(thet)*uvs(2,i)+1.1d0
        srcvals(3,i,2) = 0
        
        srcvals(4,i,2) = cos(thet)
        srcvals(5,i,2) = -sin(thet)
        srcvals(6,i,2) = 0

        srcvals(7,i,2) = sin(thet)
        srcvals(8,i,2) = cos(thet)
        srcvals(9,i,2) = 0
      enddo

      do ipatch=1,npatches
        do i=1,npols
          do j=1,9
            srccoefs(j,i,ipatch) = 0
            do l=1,npts
               srccoefs(j,i,ipatch) = srccoefs(j,i,ipatch) + umatr(i,l)*
     1            srcvals(j,l,ipatch)
            enddo
          enddo
        enddo
      enddo



c
cc       note the exact integrals are hardcoded for the
c        the following geometry
c
c        update documentation needed
c
c        see mathematica notebook quadintegrals.nb
c        
      

      zk = 1.1d0

      zpars(1) = zk
      xyztarg(1,1) = 0
      xyztarg(2,1) = 0.1d0*rr
      xyztarg(3,1) = 0.019d0*rr

      xyztarg(1,2) = 1.2d0*rr
      xyztarg(2,2) = 0.1d0*rr
      xyztarg(3,2) = 0.07d0*rr

      xyztarg(1,3) = cos(thet)*0.3d0 + sin(thet)*1.4d0+2.1d0
      xyztarg(2,3) = -sin(thet)*0.3d0 + cos(thet)*1.4d0+1.1d0
      xyztarg(3,3) = -3.0d0/7.0d0


      allocate(slp(npols,ntarg),slp_ex(npols,ntarg))

      do i=1,ntarg
        do j=1,npols
          slp(j,i) = 0
          slp_ex(j,i) = 0
        enddo
      enddo

      eps = 1.0d-7
      nqorder = 20
      nquadmax = 5000
      nn = norder + 2
c
c  fix exact integrals
c
      if(ipoly.eq.0) then

      slp_ex(1,1) = 6.918455037474758d0*rr 
      slp_ex(2,1) = 0
      slp_ex(3,1) = -1.169804038261973d0*rr
      slp_ex(nn,1) = 0.3395839474103037d0*rr 

      slp_ex(1,2) = 3.647562052924359d0*rr 
      slp_ex(2,2) = 0.9569212230520274d0*rr
      slp_ex(3,2) = 0.2829375669797821d0*rr
      slp_ex(nn,2)= 0.09916135325699458d0*rr

      slp_ex(1,3) = 2.787747194509981d0 
      slp_ex(2,3) = 0.1395798678743148d0 
      slp_ex(3,3) = -0.09564700500083024d0
      slp_ex(nn,3) = 0.5439248539409640d0

      elseif(ipoly.eq.1) then

      slp_ex(1,1) = 6.918455037474758d0*rr 
      slp_ex(2,1) = 0
      slp_ex(3,1) = -3.865890396840883d0*rr
      slp_ex(nn,1) = 0.3395839474103037d0*rr 

      slp_ex(1,2) = 3.647562052924359d0*rr 
      slp_ex(2,2) = 0.9569212230520274d0*rr
      slp_ex(3,2) = -0.8386039283350767d0*rr
      slp_ex(nn,2)= 0.09916135325699458d0*rr

      slp_ex(1,3) = 2.787747194509981d0 
      slp_ex(2,3) = 0.1395798678743148d0 
      slp_ex(3,3) = -1.056778404837767d0
      slp_ex(nn,3) = 0.5439248539409640d0


      endif

      rfac = 3.0d0
      ifmetric = 0
      rn1 = 0
      n2 = 0
      itargptr(1) = 1
      itargptr(2) = 3

      ntargptr(1) = 2
      ntargptr(2) = 1
      call dquadints(eps,istrat,intype,npatches,norder,ipoly,ttype,
     1  npols,srccoefs,3,ntarg,xyztarg,ifp,xyztarg,itargptr,ntargptr,
     2  norder,npols,lslp,0,dpars,1,zpars,0,ipars,nqorder,nquadmax,
     3  rfac,slp,ifmetric,rn1,n2)

      
      errmax = 0.0d0
      do i=1,ntarg
        print *, "itarg =", i
        print *, ""
        print *, ""

        err1s = abs(slp(1,i)-slp_ex(1,i))
        err2s = abs(slp(2,i)-slp_ex(2,i))
        err3s = abs(slp(nn,i)-slp_ex(nn,i))
        err4s = abs(slp(3,i)-slp_ex(3,i))


        if(err1s.gt.errmax) errmax = err1s
        if(err2s.gt.errmax) errmax = err2s
        if(err3s.gt.errmax) errmax = err3s
        if(err4s.gt.errmax) errmax = err4s

        write(*,'(a,2(2x,e11.5))'), "0,0 int=", slp(1,i),slp_ex(1,i)
        write(*,'(a,2(2x,e11.5))'), "1,0 int=", slp(2,i),slp_ex(2,i)
        write(*,'(a,2(2x,e11.5))'), "2,0 int=", slp(3,i),slp_ex(3,i)
        write(*,'(a,2(2x,e11.5))'), "0,1 int=", slp(nn,i),slp_ex(nn,i)

        print *, "0,0 int=",err1s
        print *, "1,0 int=",err2s
        print *, "0,1 int=",err3s
        print *, "2,0 int=",err4s
        print *, ""
        print *, ""
      enddo

      
      i1 = 0
      if(errmax.lt.eps) i1 = 1

      isuccess = i1

      return
      end
c
c
c
c
c
c

      subroutine lslp(x,ndt,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(3),y(ndt),dpars(ndd)
      complex *16 zpars(ndz),ima
      data ima/(0.0d0,1.0d0)/
      integer ipars(ndi)
      real *8 f

      rr = sqrt((x(1)-y(1))**2 + (x(2)-y(2))**2 + (x(3)-y(3))**2)

      f = 1.0d0/rr

      return
      end

c      
c      
c
c

