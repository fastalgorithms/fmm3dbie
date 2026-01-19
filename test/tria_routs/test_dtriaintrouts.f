      subroutine test_dtriarouts(nsuccess)
      
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      real *8 uvs(2,10000),wts(10000)
      real *8, allocatable :: umatr(:,:),vmatr(:,:)
      real *8, allocatable :: xyzvals(:,:,:)
      real *8, allocatable :: srccoefs(:,:,:)
      real *8, allocatable :: dxyzvals(:,:,:)
      real *8, allocatable :: qnodes(:,:),qwts(:)
      real *8 xyztarg(3,1000)
      real *8, allocatable :: cintvals(:,:)
      real *8, allocatable :: cintex(:,:)
      integer *8 itargptr(1000)
      integer *8 ntargptr(1000)
      integer *8 ipars(10)
      integer *8, allocatable :: iprint(:,:)

      real *8 dpars(5)

      complex * 16 zpars(3)

      external dlslp

      done = 1
      pi = atan(done)*4
      ndim = 3


      call prini(6,13)

      nsuccess = 0


      npatches = 2

      norder = 8
      npols = (norder+1)*(norder+2)/2
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

      rr = 1.0d-2

      do i=1,npols
        xyzvals(1,i,1) = uvs(1,i)*rr 
        xyzvals(2,i,1) = uvs(2,i)*rr
        xyzvals(3,i,1) = 0

        dxyzvals(1,i,1) = rr
        dxyzvals(2,i,1) = 0
        dxyzvals(3,i,1) = 0

        dxyzvals(4,i,1) = 0
        dxyzvals(5,i,1) = rr
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
      xyztarg(2,1) = 0.1d0*rr
      xyztarg(3,1) = 0.4d0*rr

      xyztarg(1,2) = 0.019d0*rr
      xyztarg(2,2) = -0.019d0*rr
      xyztarg(3,2) = 0

      
      xyztarg(1,3) = 2.817d0
      xyztarg(2,3) = 1.0d0
      xyztarg(3,3) = 0

      itargptr(1) = 1
      itargptr(2) = 3
      
      ntargptr(1) = 2
      ntargptr(2) = 1

      nporder = 10
      nppols = (nporder+1)*(nporder+2)/2

      allocate(cintvals(nppols,ntarg),cintex(nppols,ntarg))

      do i=1,ntarg
        do j=1,nppols
          cintvals(j,i) = 0
          cintex(j,i) = 0
        enddo
      enddo

      rfac= 3.0d0
      cintex(1,1) = 1.156374590854795d0*rr 
      cintex(2,1) = 0.0942702604583580d0*rr 
      cintex(3,1) = -0.239019699324178d0*rr

      cintex(1,2) = 1.7158925715165438d0*rr 
      cintex(2,2) = 0.588248265036294074d0*rr 
      cintex(3,2) = -0.957021021711130743d0*rr

      cintex(1,3) = 2.341262814070374d0 
      cintex(2,3) = -0.05030309850820814d0 
      cintex(3,3) = 0.8574902331922906d0

      eps = 1.0d-6

      ndd = 1
      ndi = 1
      ndz = 1


      nqorder = 16
c
c
c    
c
      write(*,*) '======================='
      write(*,*) 'Testing strategy 1, no proxy points, RV nodes'
      write(*,*) ' '
      write(*,*) ' '

      istrat = 1
      intype = 1
      ifp = 0
      ifmetric = 0
      ntrimax = 3000

      

      call dtriaints(eps,istrat,intype,npatches,norder,npols,srccoefs,
     1      ndim,ntarg,xyztarg,ifp,tmp,itargptr,
     2      ntargptr,nporder,nppols,dlslp,ndd,dpars,ndz,zpars,ndi,ipars,
     3      nqorder,
     3      ntrimax,rfac,cintvals,ifmetric,rn1,n2)

      erra = 0
      ra = 0
      do i=1,ntarg
        do j=1,3
          ra = ra + abs(cintex(j,i))**2
          erra = erra + abs(cintex(j,i)-cintvals(j,i))**2
        enddo

      enddo

      erra = sqrt(erra/ra)
      if(erra.lt.eps) nsuccess = nsuccess+1
      call prin2('error in integrals = *',erra,1)

c
c    
c
      write(*,*) '======================='
      write(*,*) 'Testing strategy 1, no proxy points, XG nodes'
      write(*,*) ' '
      write(*,*) ' '

      istrat = 1
      intype = 2
      ifp = 0
      ifmetric = 0

      

      call dtriaints(eps,istrat,intype,npatches,norder,npols,srccoefs,
     1      ndim,ntarg,xyztarg,ifp,tmp,itargptr,
     2      ntargptr,nporder,nppols,dlslp,ndd,dpars,ndz,zpars,ndi,ipars,
     3      nqorder,ntrimax,rfac,cintvals,ifmetric,rn1,n2)

      erra = 0
      ra = 0
      do i=1,ntarg
        do j=1,3
          ra = ra + abs(cintex(j,i))**2
          erra = erra + abs(cintex(j,i)-cintvals(j,i))**2
        enddo

      enddo

      erra = sqrt(erra/ra)
      if(erra.lt.eps) nsuccess = nsuccess + 1
      call prin2('error in integrals = *',erra,1)


 2000 continue

c
c    
c
      write(*,*) '======================='
      write(*,*) 'Testing strategy 2, no proxy points, RV nodes'
      write(*,*) ' '
      write(*,*) ' '

      istrat = 2
      intype = 1
      ifp = 0
      ifmetric = 0


      call dtriaints(eps,istrat,intype,npatches,norder,npols,srccoefs,
     1      ndim,ntarg,xyztarg,ifp,tmp,itargptr,
     2      ntargptr,nporder,nppols,dlslp,ndd,dpars,ndz,zpars,ndi,ipars,
     3      nqorder,ntrimax,rfac,cintvals,ifmetric,rn1,n2)

      erra = 0
      ra = 0
      ifprint = 0
      do i=1,ntarg
        if(ifprint.ge.1) call prinf('itarg=*',i,1)
        if(ifprint.ge.1) call prin2('cintex=*',cintex(1,i),6)
        if(ifprint.ge.1) call prin2('cintvals=*',cintvals(1,i),6)
        do j=1,3
          ra = ra + abs(cintex(j,i))**2
          erra = erra + abs(cintex(j,i)-cintvals(j,i))**2
        enddo
      enddo


      erra = sqrt(erra/ra)
      if(erra.lt.eps) nsuccess = nsuccess + 1
      call prin2('error in integrals = *',erra,1)

c
c    
c
      write(*,*) '======================='
      write(*,*) 'Testing strategy 2, no proxy points, XG nodes'
      write(*,*) ' '
      write(*,*) ' '

      istrat = 2
      intype = 2
      ifp = 0
      ifmetric = 0


      call dtriaints(eps,istrat,intype,npatches,norder,npols,srccoefs,
     1      ndim,ntarg,xyztarg,ifp,tmp,itargptr,
     2      ntargptr,nporder,nppols,dlslp,ndd,dpars,ndz,zpars,ndi,ipars,
     3      nqorder,ntrimax,rfac,cintvals,ifmetric,rn1,n2)

      erra = 0
      ra = 0
      do i=1,ntarg
        if(ifprint.ge.1) call prinf('itarg=*',i,1)
        if(ifprint.ge.1) call prin2('cintex=*',cintex(1,i),6)
        if(ifprint.ge.1) call prin2('cintvals=*',cintvals(1,i),6)
        do j=1,3
          ra = ra + abs(cintex(j,i))**2
          erra = erra + abs(cintex(j,i)-cintvals(j,i))**2
        enddo
      enddo


      erra = sqrt(erra/ra)
      if(erra.lt.eps) nsuccess = nsuccess + 1
      call prin2('error in integrals = *',erra,1)


 3000 continue


      rfac = 0.8d0
      nqorder = 6
      eps = 1.0d-6

c
c    
c
      write(*,*) '======================='
      write(*,*) 'Testing strategy 3, no proxy points, RV nodes'
      write(*,*) ' '
      write(*,*) ' '

      istrat = 3
      intype = 1
      ifp = 0
      ifmetric = 0


      call dtriaints(eps,istrat,intype,npatches,norder,npols,srccoefs,
     1      ndim,ntarg,xyztarg,ifp,tmp,itargptr,
     2      ntargptr,nporder,nppols,
     3      dlslp,ndd,dpars,ndz,zpars,ndi,ipars,nqorder,ntrimax,
     3      rfac,cintvals,ifmetric,rn1,n2)

      erra = 0
      ra = 0
      ifprint = 0
      do i=1,ntarg
        if(ifprint.ge.1) call prinf('itarg=*',i,1)
        if(ifprint.ge.1) call prin2('cintex=*',cintex(1,i),6)
        if(ifprint.ge.1) call prin2('cintvals=*',cintvals(1,i),6)
        do j=1,3
          ra = ra + abs(cintex(j,i))**2
          erra = erra + abs(cintex(j,i)-cintvals(j,i))**2
        enddo
      enddo


      erra = sqrt(erra/ra)
      if(erra.lt.eps) nsuccess = nsuccess + 1
      call prin2('error in integrals = *',erra,1)

c
c    
c
      write(*,*) '======================='
      write(*,*) 'Testing strategy 3, no proxy points, XG nodes'
      write(*,*) ' '
      write(*,*) ' '

      istrat = 3
      intype = 2
      ifp = 0
      ifmetric = 0

      nqorder = 10


      call dtriaints(eps,istrat,intype,npatches,norder,npols,srccoefs,
     1      ndim,ntarg,xyztarg,ifp,tmp,itargptr,
     2      ntargptr,nporder,nppols,
     3      dlslp,ndd,dpars,ndz,zpars,ndi,ipars,nqorder,ntrimax,rfac,
     3      cintvals,ifmetric,rn1,n2)

      erra = 0
      ra = 0

      ifprint = 0
      do i=1,ntarg
        if(ifprint.ge.1) call prinf('itarg=*',i,1)
        if(ifprint.ge.1) call prin2('cintex=*',cintex(1,i),6)
        if(ifprint.ge.1) call prin2('cintvals=*',cintvals(1,i),6)
        do j=1,3
          ra = ra + abs(cintex(j,i))**2
          erra = erra + abs(cintex(j,i)-cintvals(j,i))**2
        enddo
      enddo



      erra = sqrt(erra/ra)
      if(erra.lt.eps) nsuccess = nsuccess + 1
      call prin2('error in integrals = *',erra,1)
c
c
c
c
 4000 continue
      write(*,*) '======================='
      write(*,*) 'Testing wnodes XG nodes'
      write(*,*) ' '
      write(*,*) ' '

      nqorder = 20
      nlev = 3
      
      call triasymq_pts(nqorder,nnodes)
      
      ntri = 4**nlev
      npts = ntri*nnodes
      allocate(qnodes(2,npts),qwts(npts))

      call gen_xg_unif_nodes_tri(nlev,nqorder,nnodes,npts,qnodes,qwts)



      call dtriaints_wnodes(npatches,norder,npols,srccoefs,
     1      ndim,ntarg,xyztarg,itargptr,
     2      ntargptr,nporder,nppols,
     3      dlslp,ndd,dpars,ndz,zpars,ndi,ipars,npts,qnodes,
     3      qwts,cintvals)

      erra = 0
      ra = 0

      do i=1,ntarg
        if(ifprint.ge.1) call prinf('itarg=*',i,1)
        if(ifprint.ge.1) call prin2('cintex=*',cintex(1,i),6)
        if(ifprint.ge.1) call prin2('cintvals=*',cintvals(1,i),6)
        do j=1,3
          ra = ra + abs(cintex(j,i))**2
          erra = erra + abs(cintex(j,i)-cintvals(j,i))**2
        enddo
      enddo


      erra = sqrt(erra/ra)
      if(erra.lt.1.0d-3) nsuccess = nsuccess + 1
      call prin2('error in integrals = *',erra,1)

      return
      end
c
c
c
c
      subroutine test_dtriarouts_vec(nsuccess)
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      real *8 uvs(2,10000),wts(10000)
      real *8, allocatable :: umatr(:,:),vmatr(:,:)
      real *8, allocatable :: xyzvals(:,:,:)
      real *8, allocatable :: srccoefs(:,:,:)
      real *8, allocatable :: dxyzvals(:,:,:)
      real *8 xyztarg(3,1000)
      real *8, allocatable :: cintvals(:,:,:)
      real *8, allocatable :: cintex(:,:,:)
      integer *8 itargptr(1000)
      integer *8 ntargptr(1000)
      integer *8 ipars(10)
      integer *8, allocatable :: iprint(:,:)

      real *8 dpars(5)

      complex * 16 zpars(3)

      external dvslp

      done = 1
      pi = atan(done)*4
      ndim = 3


      call prini(6,13)
      nsuccess = 0


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

      cintex(2,1,1) = 0.8936785872359208d0 
      cintex(2,2,1) = 0.11034267767644536d0 
      cintex(2,3,1) = -0.2860965856806985d0 

      cintex(2,1,2) = 1.493516353671125d0 
      cintex(2,2,2) = 0.6294936116618000d0 
      cintex(2,3,2) = -1.0198449639415986d0 

      cintex(2,1,3) = 2.175031826546217d0 
      cintex(2,2,3) = -0.04227033492325212d0 
      cintex(2,3,3) = 0.8955868465640032d0 

      eps = 1.0d-6

      ndd = 1
      ndi = 1
      ndz = 1


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

      

      call dtriaints_vec(eps,istrat,intype,npatches,norder,npols,
     1      srccoefs,ndim,ntarg,xyztarg,ifp,tmp,itargptr,
     2      ntargptr,nporder,nppols,dvslp,nd,ndd,dpars,ndz,zpars,ndi,
     3      ipars,nqorder,ntrimax,rfac,cintvals,ifmetric,rn1,n2)

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
      if(erra.lt.eps) nsuccess = nsuccess+1
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

      

      call dtriaints_vec(eps,istrat,intype,npatches,norder,npols,
     1      srccoefs,ndim,ntarg,xyztarg,ifp,tmp,itargptr,
     2      ntargptr,nporder,nppols,dvslp,nd,ndd,dpars,ndz,zpars,ndi,
     3      ipars,nqorder,ntrimax,rfac,cintvals,
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
      if(erra.lt.eps) nsuccess = nsuccess+1
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


      call dtriaints_vec(eps,istrat,intype,npatches,norder,npols,
     1      srccoefs,ndim,ntarg,xyztarg,ifp,tmp,itargptr,
     2      ntargptr,nporder,nppols,dvslp,nd,ndd,dpars,ndz,zpars,ndi,
     3      ipars,nqorder,
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
      if(erra.lt.eps) nsuccess = nsuccess+1
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



      call dtriaints_vec(eps,istrat,intype,npatches,norder,npols,
     1      srccoefs,ndim,ntarg,xyztarg,ifp,tmp,itargptr,
     2      ntargptr,nporder,nppols,dvslp,nd,ndd,dpars,ndz,zpars,ndi,
     3      ipars,nqorder,ntrimax,rfac,cintvals,ifmetric,rn1,n2)

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
      if(erra.lt.eps) nsuccess = nsuccess+1
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

      call dtriaints_vec(eps,istrat,intype,npatches,norder,npols,
     1      srccoefs,ndim,ntarg,xyztarg,ifp,tmp,itargptr,
     2      ntargptr,nporder,nppols,dvslp,nd,ndd,dpars,ndz,zpars,ndi,
     3      ipars,nqorder,ntrimax,rfac,cintvals,ifmetric,rn1,n2)

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
      if(erra.lt.eps) nsuccess = nsuccess+1
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


      call dtriaints_vec(eps,istrat,intype,npatches,norder,npols,
     1      srccoefs,ndim,ntarg,xyztarg,ifp,tmp,itargptr,
     2      ntargptr,nporder,nppols,dvslp,nd,ndd,dpars,ndz,zpars,ndi,
     3      ipars,nqorder,ntrimax,rfac,cintvals,ifmetric,rn1,n2)

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
      if(erra.lt.eps) nsuccess = nsuccess+1
      call prin2('error in integrals = *',erra,1)


      return
      end
c
c
c
c
c

      subroutine dvslp(nd,x,ndt,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      real *8 x(*),y(ndt),dpars(ndd)
      complex *16 zpars(ndz)
      integer *8 ipars(ndi)
      integer *8 ndd,ndi,ndz
      real *8 f(2)
      
      rr = (x(1)-y(1))**2 + (x(2)-y(2))**2 + (x(3)-y(3))**2
      rr = sqrt(rr)
      f(1) = 1/rr
      f(2) = cos(1.1d0*rr)/rr

      return
      end


      
      


c

      subroutine dlslp(x,ndt,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      real *8 x(*),y(ndt),dpars(ndd)
      complex *16 zpars(ndz)
      integer *8 ipars(ndi)
      real *8 f
      
      rr = (x(1)-y(1))**2 + (x(2)-y(2))**2 + (x(3)-y(3))**2
      f = 1/sqrt(rr)

      return
      end


      
      


