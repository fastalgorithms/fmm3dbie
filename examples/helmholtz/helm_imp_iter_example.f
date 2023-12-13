      implicit real *8 (a-h,o-z) 
      implicit integer(8) (i-n)
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:)
      character *100 fname
      integer(8) ipars(2)

      integer(8), allocatable :: norders(:),ixyzs(:),iptype(:)

      real *8 xyz_out(3),xyz_in(3),xyz_src(3),xyz_targ(3)
      complex *16, allocatable :: sigma(:),rhs(:),sigma1(:)
      complex *16, allocatable :: zlams(:)
      real *8, allocatable :: errs(:)
      real *8 eps_gmres
      complex * 16 zpars(3)
      integer(8) numit,niter

      integer(8) ipatch_id
      real *8 uvs_targ(2)

      logical isout0,isout1

      complex *16 pot,potex,ztmp,ima

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4


c
c       select geometry type
c       igeomtype = 1 => sphere
c 
      igeomtype = 1
c
c    
c     iimp_fun = 1, then impedance = constant
c     iimp_fun = 2, then impednace = (2+z)/5
c     iimp_fun = 3, then impedance = (1+y)
c
      iimp_fun = 3

      ipars(1) = 1
      npatches = 12*(4**ipars(1))

      zk = 1.11d0+ima*0.1d0
      zpars(1) = zk 
      zpars(2) = 1.0d0

      xyz_out(1) = 3.17d0
      xyz_out(2) = -0.03d0
      xyz_out(3) = 3.15d0

      xyz_in(1) = 0.17d0
      xyz_in(2) = 0.23d0
      xyz_in(3) = -0.11d0

      norder = 4 
      npols = (norder+1)*(norder+2)/2

      ifinout = 0
      if(ifinout.eq.0) then
        xyz_src(1) = xyz_out(1)
        xyz_src(2) = xyz_out(2)
        xyz_src(3) = xyz_out(3)

        xyz_targ(1) = xyz_in(1)
        xyz_targ(2) = xyz_in(2)
        xyz_targ(3) = xyz_in(3)
      endif

      if(ifinout.eq.1) then
        xyz_src(1) = xyz_in(1)
        xyz_src(2) = xyz_in(2)
        xyz_src(3) = xyz_in(3)

        xyz_targ(1) = xyz_out(1)
        xyz_targ(2) = xyz_out(2)
        xyz_targ(3) = xyz_out(3)
      endif
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

      allocate(sigma(npts),rhs(npts),sigma1(npts))
      allocate(zlams(npts))

      if(iimp_fun.eq.1) then
        do i=1,npts
          zlams(i) = 1.0d0
        enddo
      endif

      if(iimp_fun.eq.2) then
        do i=1,npts
          zlams(i) = (2+srcvals(3,i))/5.0d0
        enddo
      endif

      if(iimp_fun.eq.3) then
        do i=1,npts
          zlams(i) = 1+ srcvals(2,i)
        enddo

      endif
     

      do i=1,npts
        call h3d_sprime(xyz_src,int(12,8),srcvals(1,i),int(0,8),dpars,
     1     int(1,8),zpars,int(0,8),ipars,rhs(i))
        
        call h3d_slp(xyz_src,int(12,8),srcvals(1,i),int(0,8),dpars,
     1     int(1,8),zpars,int(0,8),ipars,ztmp)
        rhs(i) = rhs(i) + ima*zpars(1)*zlams(i)*ztmp
        sigma(i) = 0
        sigma1(i) = 0
      enddo


      numit = 200
      niter = 0
      allocate(errs(numit+1))

      eps = 0.51d-4

      eps_gmres = 0.5d-6

      call helm_rpcomb_imp_solver_memest(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,eps,zpars,numit,rmem)
      
      call prin2('estimated memory=*',rmem,1)


      call helm_rpcomb_imp_solver(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,eps,zpars,zlams,numit,ifinout,rhs,eps_gmres,
     2  niter,errs,rres,sigma,sigma1)

      call prinf('niter=*',niter,1)
      call prin2('rres=*',rres,1)
      call prin2('errs=*',errs,niter)


c
c       test solution at interior point
c
      call h3d_slp(xyz_targ,int(3,8),xyz_src,int(0,8),dpars,int(1,8),
     1             zpars,int(0,8),ipars,potex)

      ndtarg = 3
      ntarg = 1
      ipatch_id = -1
      uvs_targ(1) = 0
      uvs_targ(2) = 0
      call lpcomp_helm_rpcomb_dir(npatches,norders,ixyzs,iptype,
     1  npts,srccoefs,srcvals,ndtarg,ntarg,xyz_targ,ipatch_id,
     2  uvs_targ,eps,zpars,sigma,sigma1,pot)

      call prin2('potex=*',potex,2)
      call prin2('pot=*',pot,2)
      erra = abs(pot-potex)/abs(potex)
      call prin2('relative error=*',erra,1)



      stop
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
          call h3d_sprime(xyzout,int(12,8),srcvals(1,i),int(0,8),dpars,
     1       int(1,8),zk,int(0,8),ipars,val)

          call cross_prod3d(srcvals(4,i),srcvals(7,i),tmp)
          ds = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)
          ra = ra + real(val)*wts(i)
        enddo
      enddo

      if(abs(ra+4*pi).le.1.0d-1) isout = .false.
      if(abs(ra).le.1.0d-1) isout = .true.

      return
      end

   




