      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:)
      character *100 fname
      integer ipars(2)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)

      real *8 xyz_out(3),xyz_in(3),xyz_src(3),xyz_targ(3),xyz_all(3,2)
      complex *16, allocatable :: sigma(:),rhs(:),sigma1(:)
      real *8, allocatable :: errs(:)
      real *8 eps_gmres
      complex * 16 zpars(6),zpars2(3)
      complex *16 omega,ep0,ep1,mu0,mu1,zk0,zk1,ztmp,ztmp2
      complex *16 u0,dudn0,u1,dudn1
      complex *16 alpha0,alpha1,beta0,beta1
      integer numit,niter

      integer ipatch_id(2)
      real *8 uvs_targ(2,2)
      complex *16 potall(2)

      logical isout0,isout1

      complex *16 pot,potex,ima

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4


c
c       select geometry type
c       igeomtype = 1 => sphere
c 
      igeomtype = 1
      ipars(1) = 1
      npatches = 12*(4**ipars(1))

      omega = 0.33d0
      ep0 = 1.0d0
      mu0 = 1.0d0
      ep1 = 1.1d0
      mu1 = 1.2d0
      zk0 = omega*sqrt(ep0*mu0)
      zk1 = omega*sqrt(ep1*mu1)

      alpha0 = 1.0d0
      alpha1 = 1.0d0
      beta0 = 1.0d0
      beta1 = 1.0d0

      zpars(1) = zk0
      zpars(2) = alpha0
      zpars(3) = beta0
      zpars(4) = zk1
      zpars(5) = alpha1
      zpars(6) = beta1

      xyz_out(1) = 3.17d0
      xyz_out(2) = -0.03d0
      xyz_out(3) = 3.15d0

      xyz_in(1) = 0.17d0
      xyz_in(2) = 0.23d0
      xyz_in(3) = -0.11d0

      norder = 4
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
      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)

      allocate(sigma(2*npts),rhs(2*npts))

c
c  get boundary data
c
      do i=1,npts
        call h3d_slp(xyz_in,12,srcvals(1,i),0,dpars,1,zk0,0,
     1     ipars,u0)
        call h3d_sprime(xyz_in,12,srcvals(1,i),0,dpars,1,zk0,0,
     1     ipars,dudn0)
        call h3d_slp(xyz_out,12,srcvals(1,i),0,dpars,1,zk1,0,
     1     ipars,u1)
        call h3d_sprime(xyz_out,12,srcvals(1,i),0,dpars,1,zk1,0,
     1     ipars,dudn1)
        rhs(i) = alpha0*u0 - alpha1*u1
        rhs(npts+i) = beta0*dudn0 - beta1*dudn1
      enddo


      numit = 200
      niter = 0
      allocate(errs(numit+1))

      eps = 0.51d-4

      eps_gmres = 0.5d-6

      call helm_comb_trans_solver(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,eps,zpars,numit,rhs,eps_gmres,
     2  niter,errs,rres,sigma)

      call prinf('niter=*',niter,1)
      call prin2('rres=*',rres,1)
      call prin2('errs=*',errs,niter)



      ndtarg = 3
      ntarg = 2
      ipatch_id(1) = -1
      ipatch_id(2) = -1
      uvs_targ(1,1) = 0
      uvs_targ(2,1) = 0
      uvs_targ(1,2) = 0
      uvs_targ(2,2) = 0
      xyz_all(:,1) = xyz_in
      xyz_all(:,2) = xyz_out
      call helm_comb_trans_eval(npatches,norders,ixyzs,iptype,
     1  npts,srccoefs,srcvals,ndtarg,ntarg,xyz_all,ipatch_id,
     2  uvs_targ,eps,zpars,sigma,potall)

c
c       test solution at interior point
c
      pot = potall(1)
      call h3d_slp(xyz_out,3,xyz_in,0,dpars,1,zk1,0,ipars,potex)
      call prin2('potex=*',potex,2)
      call prin2('pot=*',pot,2)
      erra = abs(pot-potex)/abs(potex)
      call prin2('relative error at interior target=*',erra,1)
c
c       test solution at exterior point
c
      pot = potall(2)
      call h3d_slp(xyz_in,3,xyz_out,0,dpars,1,zk0,0,ipars,potex)
      call prin2('potex=*',potex,2)
      call prin2('pot=*',pot,2)
      erra = abs(pot-potex)/abs(potex)
      call prin2('relative error at exterior target=*',erra,1)

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
          call h3d_sprime(xyzout,12,srcvals(1,i),0,dpars,1,zk,0,ipars,
     1       val)

          call cross_prod3d(srcvals(4,i),srcvals(7,i),tmp)
          ds = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)
          ra = ra + real(val)*wts(i)
        enddo
      enddo

      if(abs(ra+4*pi).le.1.0d-1) isout = .false.
      if(abs(ra).le.1.0d-1) isout = .true.

      return
      end

   




