      implicit real *8 (a-h,o-z)
      real *8 uvs(2,10000),wts(10000)
      real *8, allocatable :: umatr(:,:),vmatr(:,:)
      real *8, allocatable :: xyzvals(:,:,:)
      real *8, allocatable :: srccoefs(:,:,:)
      real *8, allocatable :: dxyzvals(:,:,:)
      real *8 xyztarg(3,1000)
      complex *16 ima
      complex *16, allocatable :: cintvals(:,:,:)
      complex *16, allocatable :: cintex(:,:,:)
      integer itargptr(1000)
      integer ntargptr(1000)
      integer ipars(10)
      integer, allocatable :: iprint(:,:)

      data ima/(0.0d0,1.0d0)/

      real *8 dpars(5)

      complex * 16 zpars(3)

      external vslp

      done = 1
      pi = atan(done)*4


      call prini(6,13)
      !call prinf('enter n*',n,0)
      !read *, n


      npatches = 2
      nd = 2

      norder = 8
      npols = (norder+1)*(norder+2)/2


      nporder = 10
      nppols = (nporder+1)*(nporder+2)/2
      allocate(srccoefs(9,npols,npatches),xyzvals(3,npols,npatches))

      allocate(dxyzvals(6,npols,npatches))

      allocate(umatr(npols,npols),vmatr(npols,npols))
      
c
cc      generate xyzcoefs for the two patches
c
      call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)
      call koorn_vals2coefs(norder,npols,uvs,umatr)
      call koorn_coefs2vals(norder,npols,uvs,vmatr)

      thet = pi/4

      do i=1,npols
        xyzvals(1,i,1) = uvs(1,i) 
        xyzvals(2,i,1) = uvs(2,i)
        xyzvals(3,i,1) = 0

        dxyzvals(1,i,1) = 1
        dxyzvals(2,i,1) = 0
        dxyzvals(3,i,1) = 0

        dxyzvals(4,i,1) = 0
        dxyzvals(5,i,1) = 1
        dxyzvals(6,i,1) = 0

        xyzvals(1,i,2) = cos(thet)*uvs(1,i) + sin(thet)*uvs(2,i)+2.1d0
        xyzvals(2,i,2) = -sin(thet)*uvs(1,i) + cos(thet)*uvs(2,i)+1.1d0
        xyzvals(3,i,2) = 0

        dxyzvals(1,i,2) = cos(thet)
        dxyzvals(2,i,2) = -sin(thet)
        dxyzvals(3,i,2) = 0

        dxyzvals(4,i,2) = sin(thet)
        dxyzvals(5,i,2) = cos(thet)
        dxyzvals(6,i,2) = 0
      enddo

      npts = npatches*npols


      do ipatch =1,npatches
        do i=1,npols
          do j=1,3
            srccoefs(j,i,ipatch) = 0
            do l=1,npols
              srccoefs(j,i,ipatch) = srccoefs(j,i,ipatch)+umatr(i,l)*
     1                                 xyzvals(j,l,ipatch)
            enddo
          enddo

          do j=1,6
            srccoefs(j+3,i,ipatch) = 0
            do l=1,npols
              srccoefs(j+3,i,ipatch)=srccoefs(j+3,i,ipatch)+umatr(i,l)*
     1                                   dxyzvals(j,l,ipatch)
            enddo
          enddo
        enddo
      enddo

cc      call prin2('xyzcoefs=*',xyzcoefs,3*npatches*npols)
cc      call prin2('dxyzcoefs=*',dxyzcoefs,6*npatches*npols)

      ntarg = 3

c
cc       note the exact integrals are hardcoded for the
c        the following geometry
c
c        patch 1 is the standard simplex, the first
c        two targets are associated with patch 1 and are
c        located at (0,0.1,0.4), (0.019,-0.019,0)
c
c        patch 2 is the standard simplex rotated
c        by pi/4 in the anti-clockwise direction and shifted
c        by (2.1,1.1,0), and the target associated with this
c        patch is located at (2.817,1,0)
c
c        see mathematica notebook triangleintegrals.nb
c        
      
      
      xyztarg(1,1) = 0
      xyztarg(2,1) = 0.1d0
      xyztarg(3,1) = 0.4d0

      xyztarg(1,2) = 0.019d0
      xyztarg(2,2) = -0.019d0
      xyztarg(3,2) = 0

      
      xyztarg(1,3) = 2.817d0
      xyztarg(2,3) = 1.0d0
      xyztarg(3,3) = 0

      itargptr(1) = 1
      itargptr(2) = 3
      
      ntargptr(1) = 2
      ntargptr(2) = 1

      allocate(cintvals(2,nppols,ntarg),cintex(2,nppols,ntarg))

      ntrimax = 3000

      do i=1,ntarg
        do j=1,nppols
          do idim=1,nd
            cintvals(idim,j,i) = 0
            cintex(idim,j,i) = 0
          enddo
        enddo
      enddo

cc      call prin2('xyztarg=*',xyztarg,9)

      rfac= 3.0d0
      cintex(1,1,1) = 1.156374590854795d0 
      cintex(1,2,1) = 0.0942702604583580d0 
      cintex(1,3,1) = -0.239019699324178d0

      cintex(1,1,2) = 1.7158925715165438d0 
      cintex(1,2,2) = 0.588248265036294074d0 
      cintex(1,3,2) = -0.957021021711130743d0 

      cintex(1,1,3) = 2.341262814070374d0 
      cintex(1,2,3) = -0.05030309850820814d0 
      cintex(1,3,3) = 0.8574902331922906d0

      cintex(2,1,1) = 0.8936785872359208d0 + 0.7114527818057991d0*ima 
      cintex(2,2,1) = 0.11034267767644536d0 + 0.00711422324933908d0*ima
      cintex(2,3,1) = -0.2860965856806985d0 - 0.0241919666010641d0*ima

      cintex(2,1,2) = 1.493516353671125d0 + 0.726877282085253d0*ima
      cintex(2,2,2) = 0.6294936116618000d0 + 0.0161762785701602d0*ima
      cintex(2,3,2) = -1.0198449639415986d0 - 0.0234575057934004d0*ima

      cintex(2,1,3) = 2.175031826546217d0 + 0.749807030379924d0*ima
      cintex(2,2,3) = -0.04227033492325212d0 + 0.00373575741028573d0*ima
      cintex(2,3,3) = 0.8955868465640032d0 + 0.0111093163537682d0*ima

      eps = 1.0d-6



      nqorder = 16
c
c
c    
c
      write(*,*) ' '
      write(*,*) ' '
      write(*,*) '======================='
      write(*,*) 'Testing strategy 1, no proxy points, RV nodes'
      write(*,*) ' '

      istrat = 1
      intype = 1
      ifp = 0
      ifmetric = 0

      

      call ctriaints_vec(eps,istrat,intype,npatches,norder,npols,
     1      srccoefs,3,ntarg,xyztarg,ifp,tmp,itargptr,
     2      ntargptr,nporder,nppols,vslp,nd,dpars,zpars,ipars,nqorder,
     3      ntrimax,rfac,cintvals,ifmetric,rn1,n2)

      erra = 0
      ra = 0
      do i=1,ntarg
        do j=1,3
          do idim=1,nd
            ra = ra + abs(cintex(idim,j,i))**2
            erra = erra + abs(cintex(idim,j,i)-cintvals(idim,j,i))**2
          enddo
        enddo

      enddo

      erra = sqrt(erra/ra)
      call prin2('error in integrals = *',erra,1)

c
c    
c
      write(*,*) ' '
      write(*,*) ' '
      write(*,*) '======================='
      write(*,*) 'Testing strategy 1, no proxy points, XG nodes'
      write(*,*) ' '

      istrat = 1
      intype = 2
      ifp = 0
      ifmetric = 0

      

      call ctriaints_vec(eps,istrat,intype,npatches,norder,npols,
     1      srccoefs,3,ntarg,xyztarg,ifp,tmp,itargptr,
     2      ntargptr,nporder,nppols,vslp,nd,dpars,zpars,ipars,nqorder,
     3      ntrimax,rfac,cintvals,
     3      ifmetric,rn1,n2)

      erra = 0
      ra = 0
      do i=1,ntarg
        do j=1,3
          do idim=1,nd
            ra = ra + abs(cintex(idim,j,i))**2
            erra = erra + abs(cintex(idim,j,i)-cintvals(idim,j,i))**2
          enddo
        enddo

      enddo

      erra = sqrt(erra/ra)
      call prin2('error in integrals = *',erra,1)




      

 2000 continue

c
c    
c
      write(*,*) ' '
      write(*,*) ' '
      write(*,*) '======================='
      write(*,*) 'Testing strategy 2, no proxy points, RV nodes'
      write(*,*) ' '

      eps = 1.0d-6
      istrat = 2
      intype = 1
      ifp = 0
      ifmetric = 0


      call ctriaints_vec(eps,istrat,intype,npatches,norder,npols,
     1      srccoefs,3,ntarg,xyztarg,ifp,tmp,itargptr,
     2      ntargptr,nporder,nppols,vslp,nd,dpars,zpars,ipars,nqorder,
     3      ntrimax,rfac,cintvals,
     3      ifmetric,rn1,n2)

      erra = 0
      ra = 0
      do i=1,ntarg
        do j=1,3
          do idim=1,nd
            ra = ra + abs(cintex(idim,j,i))**2
            erra = erra + abs(cintex(idim,j,i)-cintvals(idim,j,i))**2
          enddo
        enddo
      enddo

      erra = sqrt(erra/ra)
      call prin2('error in integrals = *',erra,1)


c
c    
c
      write(*,*) ' '
      write(*,*) ' '
      write(*,*) '======================='
      write(*,*) 'Testing strategy 2, no proxy points, XG nodes'
      write(*,*) ' '

      istrat = 2
      intype = 2
      ifp = 0
      ifmetric = 0



      call ctriaints_vec(eps,istrat,intype,npatches,norder,npols,
     1      srccoefs,3,ntarg,xyztarg,ifp,tmp,itargptr,
     2      ntargptr,nporder,nppols,vslp,nd,dpars,zpars,ipars,nqorder,
     3      ntrimax,rfac,cintvals,ifmetric,rn1,n2)

      erra = 0
      ra = 0
      do i=1,ntarg
        do j=1,3
          do idim=1,nd
            ra = ra + abs(cintex(idim,j,i))**2
            erra = erra + abs(cintex(idim,j,i)-cintvals(idim,j,i))**2
          enddo
        enddo
      enddo

      erra = sqrt(erra/ra)
      call prin2('error in integrals = *',erra,1)


 3000 continue


      rfac = 0.8d0
      nqorder = 6
      eps = 1.0d-6

c
c    
c
      write(*,*) ' '
      write(*,*) ' '
      write(*,*) '======================='
      write(*,*) 'Testing strategy 3, no proxy points, RV nodes'
      write(*,*) ' '

      istrat = 3
      intype = 1
      ifp = 0
      ifmetric = 0

      call ctriaints_vec(eps,istrat,intype,npatches,norder,npols,
     1      srccoefs,3,ntarg,xyztarg,ifp,tmp,itargptr,
     2      ntargptr,nporder,nppols,vslp,nd,dpars,zpars,ipars,nqorder,
     3      ntrimax,rfac,cintvals,ifmetric,rn1,n2)

      erra = 0
      ra = 0
      do i=1,ntarg
        do j=1,3
          do idim=1,nd
            ra = ra + abs(cintex(idim,j,i))**2
            erra = erra + abs(cintex(idim,j,i)-cintvals(idim,j,i))**2
          enddo
        enddo
      enddo

      erra = sqrt(erra/ra)
      call prin2('error in integrals = *',erra,1)




c
c    
c
      write(*,*) ' '
      write(*,*) ' '
      write(*,*) '======================='
      write(*,*) 'Testing strategy 3, no proxy points, XG nodes'
      write(*,*) ' '

      istrat = 3
      intype = 2
      ifp = 0
      ifmetric = 0

      nqorder = 10


      call ctriaints_vec(eps,istrat,intype,npatches,norder,npols,
     1      srccoefs,3,ntarg,xyztarg,ifp,tmp,itargptr,
     2      ntargptr,nporder,nppols,vslp,nd,dpars,zpars,ipars,nqorder,
     3      ntrimax,rfac,cintvals,ifmetric,rn1,n2)

      erra = 0
      ra = 0
      do i=1,ntarg
        do j=1,3
          do idim=1,nd
            ra = ra + abs(cintex(idim,j,i))**2
            erra = erra + abs(cintex(idim,j,i)-cintvals(idim,j,i))**2
          enddo
        enddo
      enddo

      erra = sqrt(erra/ra)
      call prin2('error in integrals = *',erra,1)




      stop
      end
c
c
c
c
c

      subroutine vslp(nd,x,y,dpars,zpars,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(3),y(3),dpars(*)
      complex *16 zpars(*)
      integer ipars(*)
      complex *16 f(2),ima
      data ima/(0.0d0,1.0d0)/
      
      rr = (x(1)-y(1))**2 + (x(2)-y(2))**2 + (x(3)-y(3))**2
      rr = sqrt(rr)
      f(1) = 1/rr
      f(2) = exp(ima*1.1d0*rr)/rr

      return
      end


      
      


