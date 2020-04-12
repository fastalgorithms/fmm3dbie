      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:)
      integer ipars(2)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)

      real *8 xyz_out(3),xyz_in(3)
      complex *16, allocatable :: sigma(:),rhs(:)
      complex *16, allocatable :: sigma2(:),rhs2(:)
      real *8, allocatable :: errs(:)
      real *8 thet,phi,eps_gmres
      complex * 16 zpars(3)
      integer numit,niter
      character *100 title,fname

      integer ipatch_id
      real *8 uvs_targ(2)

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

      
      xyz_in(1) = 2.01d0
      xyz_in(2) = 0.0d-5
      xyz_in(3) = 0.5d0

      xyz_out(1) = -3.5d0
      xyz_out(2) = 3.1d0
      xyz_out(3) = 20.1d0

      fname = '../../geometries/A380_Final_o03_r00.go3'
      
      call open_gov3_geometry_mem(fname,npatches,npts)

      call prinf('npatches=*',npatches,1)
      call prinf('npts=*',npts,1)

      allocate(srcvals(12,npts),srccoefs(9,npts))
      allocate(ixyzs(npatches+1),iptype(npatches),norders(npatches))
      allocate(wts(npts))

      call open_gov3_geometry(fname,npatches,norders,ixyzs,
     1   iptype,npts,srcvals,srccoefs,wts)
      
      fname = 'a380.vtk'
      title = 'a380'
      call surf_vtk_plot(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,fname,title)

      norder = norders(1)
      call test_exterior_pt(npatches,norder,npts,srcvals,srccoefs,
     1  wts,xyz_in,isout0)
      
      call test_exterior_pt(npatches,norder,npts,srcvals,srccoefs,
     1   wts,xyz_out,isout1)

       print *, isout0,isout1


      allocate(sigma(npts),rhs(npts))
      allocate(sigma2(npts),rhs2(npts))
      ifinout = 1

      thet = hkrand(0)*pi
      phi = hkrand(0)*2*pi
      do i=1,npts
        if(ifinout.eq.0) 
     1     call h3d_slp(xyz_out,srcvals(1,i),dpars,zpars,ipars,rhs(i))
        if(ifinout.eq.1) 
     1     call h3d_slp(xyz_in,srcvals(1,i),dpars,zpars,ipars,rhs(i))
        rhs(i) = rhs(i)
        x = srcvals(1,i)*cos(thet)
        y = srcvals(2,i)*sin(thet)*cos(phi)
        z = srcvals(3,i)*sin(thet)*sin(phi)
        rhs2(i) = exp(ima*zk*(x+y+z))
        sigma(i) = 0
        sigma2(i) = 0
      enddo


      numit = 200
      niter = 0
      allocate(errs(numit+1))

      eps = 0.51d-6
      eps_gmres = eps

      call helm_comb_dir_solver(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,eps,zpars,numit,ifinout,rhs,eps_gmres,
     2  niter,errs,rres,sigma)

      call prinf('niter=*',niter,1)
      call prin2('rres=*',rres,1)
      call prin2('errs=*',errs,niter)


c
c       test solution at interior point
c
      call h3d_slp(xyz_out,xyz_in,dpars,zpars,ipars,potex)
      pot = 0
      do i=1,npts
        if(ifinout.eq.0) 
     1     call h3d_comb(srcvals(1,i),xyz_in,dpars,zpars,ipars,ztmp)
        if(ifinout.eq.1) 
     1     call h3d_comb(srcvals(1,i),xyz_out,dpars,zpars,ipars,ztmp)
        pot = pot + sigma(i)*wts(i)*ztmp
      enddo

      call prin2('potex=*',potex,2)
      call prin2('pot=*',pot,2)
      erra = abs(pot-potex)/abs(potex)
      call prin2('relative error=*',erra,1)

      allocate(rsigma(npts))
      do i=1,npts
        rsigma(i) = real(sigma(i))*100
      enddo

      fname = 'a380_rsigma_ext.vtk'
      title = 'a380 - real part'
      call surf_vtk_plot_scalar(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,rsigma,fname,title)

      call helm_comb_dir_solver(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,eps,zpars,numit,ifinout,rhs2,niter,errs,
     2  rres,sigma2)

      do i=1,npts
        rsigma(i) = real(sigma2(i))*100
      enddo

      fname = 'a380_rsigma2_ext.vtk'
      title = 'a380 - real part - pw data'
      call surf_vtk_plot_scalar(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,rsigma,fname,title)

      stop
      end




      subroutine setup_geom(igeomtype,norder,npatches,ipars, 
     1    srcvals,srccoefs,ifplot,fname)
      implicit real *8 (a-h,o-z)
      integer igeomtype,norder,npatches,ipars(*),ifplot
      character (len=*) fname
      real *8 srcvals(12,*), srccoefs(9,*)
      real *8, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts(:)

      real *8, pointer :: ptr1,ptr2,ptr3,ptr4
      integer, pointer :: iptr1,iptr2,iptr3,iptr4
      real *8, target :: p1(10),p2(10),p3(10),p4(10)
      real *8, allocatable, target :: triaskel(:,:,:)
      real *8, allocatable, target :: deltas(:,:)
      integer, allocatable :: isides(:)
      integer, target :: nmax,mmax

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
        vmin = 0
        vmax = 2*pi
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
c    npatches - integer
c       number of patches
c    norder - integer
c       order of discretization
c    npts - integer
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
      integer npatches,norder,npts,npols
      real *8 srccoefs(9,npts),srcvals(12,npts),xyzout(3),wts(npts)
      real *8 tmp(3)
      real *8 dpars,done,pi
      real *8, allocatable :: rsurf(:),err_p(:,:) 
      integer ipars,norderhead,nd
      complex *16, allocatable :: sigma_coefs(:,:), sigma_vals(:,:)
      complex *16 zk,val

      integer ipatch,j,i
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

   




