      implicit real *8 (a-h,o-z)
      implicit integer(8) (i-n)

      ntests = 2
      nsuccess = 0
      call test_surf_lap(nsuccess)


      open(unit=33,file='../../print_testres.txt',access='append')
      write(33,'(a,i1,a,i1,a)') 'Successfully completed ',nsuccess,
     1  ' out of ',ntests,' in surface routs testing suite'
      close(33)

      stop
      end
c
c
c
c
c
      subroutine test_surf_lap(nsuccess)
      implicit real *8 (a-h,o-z) 
      implicit integer(8) (i-n)
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:)
      integer(8) ipars(2)

      integer(8), allocatable :: norders(:),ixyzs(:),iptype(:)

      real *8 xyz_out(3),xyz_in(3)
      real *8, allocatable :: ffform(:,:,:),ffformex(:,:,:)
      real *8, allocatable :: ffforminv(:,:,:),ffformexinv(:,:,:)

      complex *16, allocatable :: sigma(:),rhs(:),drhs(:,:)
      complex *16, allocatable :: drhs_cart_ex(:,:)
      complex *16, allocatable :: drhs_cart(:,:)
      complex *16, allocatable :: sigma2(:),rhs2(:)
      real *8, allocatable :: errs(:)
      real *8 thet,phi,eps_gmres
      complex * 16 zpars(3)
      integer(8) numit,niter
      character *100 title,dirname
      character *300 fname

      integer(8) ipatch_id
      real *8 uvs_targ(2)
      real *8, allocatable :: w(:,:)

      logical isout0,isout1

      complex *16 pot,potex,ztmp,ima

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4


      zk = 4.4d0+ima*0.0d0
      zpars(1) = zk 
      zpars(2) = -ima*zk
      zpars(3) = 2.0d0

      
      xyz_in(1) = 0.11d0
      xyz_in(2) = 0.0d-5
      xyz_in(3) = 0.37d0

      xyz_out(1) = -3.5d0
      xyz_out(2) = 3.1d0
      xyz_out(3) = 20.1d0

      igeomtype = 1
      ipars(1) = 4
      npatches=12*(4**ipars(1))

      norder = 8 
      npols = (norder+1)*(norder+2)/2

      npts = npatches*npols
      allocate(srcvals(12,npts),srccoefs(9,npts))
      ifplot = 0

      call setup_geom(igeomtype,norder,npatches,ipars, 
     1       srcvals,srccoefs,ifplot,fname)

      allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))

      do i=1,npatches
        norders(i) = norder
        ixyzs(i) = 1 +(i-1)*npols
        iptype(i) = 1
      enddo

      ixyzs(npatches+1) = 1+npols*npatches

c      call surf_quadratic_msh_vtk_plot(npatches,norders,ixyzs,iptype,
c     1  npts,srccoefs,srcvals,'sph-msh.vtk','msh')
      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)


      allocate(sigma(npts),rhs(npts),rhs2(npts))
      allocate(ffform(2,2,npts))

c
c       define rhs to be one of the ynm's
c
      nn = 2
      mm = 1
      nmax = nn
      allocate(w(0:nmax,0:nmax))
      call l3getsph(nmax,mm,nn,int(12,8),srcvals,rhs,npts,w)

      allocate(drhs(2,npts),drhs_cart(3,npts),drhs_cart_ex(3,npts))

      call get_surf_grad(int(2,8),npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,rhs,drhs)
       
c
c      estimate correct scaling
c
      erra = 0
      do i=1,npts
        call cart2polar(srcvals(1,i),r,theta,phi)
        ztmp = -sqrt(15.0d0/2.0d0)*exp(ima*phi)
        ct = cos(theta)
        st = sin(theta)
        cp = cos(phi)
        sp = sin(phi)
        c2t = cos(2*theta)
        drhs_cart_ex(1,i) = ztmp*(c2t*ct*cp - ima*ct*sp) 
        drhs_cart_ex(2,i) = ztmp*(c2t*ct*sp + ima*ct*cp)
        drhs_cart_ex(3,i) = -ztmp*st*c2t


        drhs_cart(1,i)=drhs(1,i)*srcvals(4,i) + drhs(2,i)*srcvals(7,i)
        drhs_cart(2,i)=drhs(1,i)*srcvals(5,i) + drhs(2,i)*srcvals(8,i)
        drhs_cart(3,i)=drhs(1,i)*srcvals(6,i) + drhs(2,i)*srcvals(9,i)

        err1 = abs(drhs_cart(1,i)-drhs_cart_ex(1,i))**2
        err2 = abs(drhs_cart(2,i)-drhs_cart_ex(2,i))**2
        err3 = abs(drhs_cart(3,i)-drhs_cart_ex(3,i))**2
        
        erra = erra + err1 + err2 + err3
      enddo
      erra = sqrt(erra/npts)

      print *, "error in surface gradient=",erra
      if(erra.lt.1.0d-8) nsuccess=nsuccess+1

      call get_surf_div(int(2,8),npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,drhs,rhs2)

      erra = 0
      do i=1,npts
        erra = erra + abs(rhs2(i)+rhs(i)*nn*(nn+1))**2
      enddo

      erra=  sqrt(erra/npts)
      print *, "error in surface divergence=",erra
      if(erra.lt.1.0d-6) nsuccess=nsuccess+1

      return
      end




      subroutine setup_geom(igeomtype,norder,npatches,ipars, 
     1    srcvals,srccoefs,ifplot,fname)
      implicit real *8 (a-h,o-z)
      implicit integer(8) (i-n)
      integer(8) igeomtype,norder,npatches,ipars(*),ifplot
      character (len=*) fname
      real *8 srcvals(12,*), srccoefs(9,*)
      real *8, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts(:)

      real *8, pointer :: ptr1,ptr2,ptr3,ptr4
      integer(8), pointer :: iptr1,iptr2,iptr3,iptr4
      real *8, target :: p1(10),p2(10),p3(10),p4(10)
      real *8, allocatable, target :: triaskel(:,:,:)
      real *8, allocatable, target :: deltas(:,:)
      integer(8), allocatable :: isides(:)
      integer(8), target :: nmax,mmax

      procedure (), pointer :: xtri_geometry


      external xtri_stell_eval,xtri_sphere_eval
      
      npols = (norder+1)*(norder+2)/2
      allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols))
      allocate(wts(npols))

      call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)

      if(igeomtype.eq.1) then
        itype = 2
        allocate(triaskel(3,3,npatches))
        allocate(isides(npatches))
        npmax = npatches
        ntri = 0
        call xtri_platonic(itype, ipars(1), npmax, ntri, 
     1      triaskel, isides)

        xtri_geometry => xtri_sphere_eval
        ptr1 => triaskel(1,1,1)
        ptr2 => p2(1)
        ptr3 => p3(1)
        ptr4 => p4(1)


        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, 
     1         ptr3,ptr4, norder,'Triangulated surface of the sphere')
        endif


        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4,
     1     npols,uvs,umatr,srcvals,srccoefs)
      endif

      if(igeomtype.eq.2) then
        done = 1
        pi = atan(done)*4
        umin = 0
        umax = 2*pi
        vmin = 2*pi
        vmax = 0
        allocate(triaskel(3,3,npatches))
        nover = 0
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2),
     1     nover,npatches,npatches,triaskel)

        mmax = 2
        nmax = 1
        xtri_geometry => xtri_stell_eval

        allocate(deltas(-1:mmax,-1:nmax))
        deltas(-1,-1) = 0.17d0
        deltas(0,-1) = 0
        deltas(1,-1) = 0
        deltas(2,-1) = 0

        deltas(-1,0) = 0.11d0
        deltas(0,0) = 1
        deltas(1,0) = 4.5d0
        deltas(2,0) = -0.25d0

        deltas(-1,1) = 0
        deltas(0,1) = 0.07d0
        deltas(1,1) = 0
        deltas(2,1) = -0.45d0

        ptr1 => triaskel(1,1,1)
        ptr2 => deltas(-1,-1)
        iptr3 => mmax
        iptr4 => nmax

        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, 
     1         iptr3,iptr4, norder,
     2         'Triangulated surface of the stellarator')
        endif

        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,iptr3,iptr4,
     1     npols,uvs,umatr,srcvals,srccoefs)
      endif
      
      return  
      end


      subroutine test_exterior_pt(npatches,norder,npts,srcvals,
     1   srccoefs,wts,xyzout,isout)
c
c
c  this subroutine tests whether the pt xyzin, is
c  in the exterior of a surface, and also estimates the error
c  in representing e^{ir/2}/r and \grad e^{ir/2}/r \cdot n
c  centered at the interior point. Whether a point 
c  is in the interior or not is tested using Gauss' 
c  identity for the flux due to a point charge
c
c
c  input:
c    npatches - integer(8)
c       number of patches
c    norder - integer(8)
c       order of discretization
c    npts - integer(8)
c       total number of discretization points on the surface
c    srccoefs - real *8 (9,npts)
c       koornwinder expansion coefficients of geometry info
c    xyzout -  real *8 (3)
c       point to be tested
c
c  output: 
c    isout - boolean
c      whether the target is in the interior or not
c

      implicit none
      integer(8) npatches,norder,npts,npols
      real *8 srccoefs(9,npts),srcvals(12,npts),xyzout(3),wts(npts)
      real *8 tmp(3)
      real *8 dpars,done,pi
      real *8, allocatable :: rsurf(:),err_p(:,:) 
      integer(8) ipars,norderhead,nd
      complex *16, allocatable :: sigma_coefs(:,:), sigma_vals(:,:)
      complex *16 zk,val

      integer(8) ipatch,j,i
      real *8 ra,ds
      logical isout

      done = 1
      pi = atan(done)*4

      npols = (norder+1)*(norder+2)/2


      zk = 0

      ra = 0



      do ipatch=1,npatches
        do j=1,npols
          i = (ipatch-1)*npols + j
          call h3d_sprime(xyzout,srcvals(1,i),dpars,zk,ipars,val)
          call cross_prod3d(srcvals(4,i),srcvals(7,i),tmp)
          ds = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)
          ra = ra + real(val)*wts(i)
        enddo
      enddo

      if(abs(ra+4*pi).le.1.0d-3) isout = .false.
      if(abs(ra).le.1.0d-3) isout = .true.

      return
      end

   



      subroutine l3getsph(nmax,mm,nn,ndx,xyzs,ynms,npts,ynm)
      implicit real *8 (a-h,o-z)
      implicit integer(8) (i-n)
      real *8 :: xyzs(ndx,npts)
      complex *16 ynms(npts),ima
      real *8 rat1(10000),rat2(10000)
      real *8 ynm(0:nmax,0:nmax)
      data ima/(0.0d0,1.0d0)/
  
      call ylgndrini(nmax, rat1, rat2)
  
      do i=1,npts
        x=xyzs(1,i)
        y=xyzs(2,i)
        z=xyzs(3,i)
        r=sqrt(x**2+y**2+z**2)
        call cart2polar(xyzs(1,i),r,theta,phi)
        ctheta = cos(theta)
        call ylgndrf(nmax, ctheta, ynm, rat1, rat2)
        ynms(i) = ynm(nn,abs(mm))*exp(ima*mm*phi)        
      enddo
       
      return
      end






