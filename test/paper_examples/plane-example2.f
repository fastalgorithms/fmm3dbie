      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:)
      integer ipars(2),ndims(3)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)
      integer, allocatable :: isout(:)

      real *8, allocatable :: xyz_in(:,:)

      complex *16, allocatable :: sigma(:),rhs(:),charges(:)
      complex *16, allocatable :: sigmatmp(:)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: targs(:,:),uvs_targ(:,:)
      complex *16, allocatable :: pvals(:),pvalsex(:)
      real *8, allocatable :: evals(:),evals2(:)
      real *8, allocatable :: errs(:)
      real *8 thet,phi
      real *8 xyz_start(3),dxyz(3)
      complex * 16 zpars(3),zpars_tmp(3)
      integer numit,niter
      character *100 title,fname



      complex *16 pot,potex,ztmp,ima

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4

c
c   simulation for plane 50 wavelengths in size
c
      zk = 28.56d0/100.0d0
      zpars(1) = zk 
      zpars(2) = -ima*zk
      zpars(3) = 2.0d0

      
      fname = '../../geometries/A380_Final_o03_r00.go3'
      
      call open_gov3_geometry_mem(fname,npatches,npts)

      call prinf('npatches=*',npatches,1)
      call prinf('npts=*',npts,1)

      allocate(srcvals(12,npts),srccoefs(9,npts))
      allocate(ixyzs(npatches+1),iptype(npatches),norders(npatches))
      allocate(wts(npts))

      call open_gov3_geometry(fname,npatches,norders,ixyzs,
     1   iptype,npts,srcvals,srccoefs,wts)

      open(unit=39,file='tail-sources.txt')
      open(unit=40,file='fuselage-sources.txt')

      read(39,*) ntail
      read(40,*) nfuse


      ntot = nfuse+ntail
      allocate(xyz_in(3,ntot),charges(ntot))
      do i=1,ntail
        read(39,*) xyz_in(1,i),xyz_in(2,i),xyz_in(3,i)
      enddo

      do i=1,nfuse
        ii = i+ntail
        read(40,*) xyz_in(1,ii),xyz_in(2,ii),xyz_in(3,ii)
      enddo

      close(39)
      close(40)

      open(unit=39,file='charges.txt')
      do i=1,ntot
        read(39,*) rtmp
        charges(i) = rtmp
      enddo

      

      allocate(sigma(npts),rhs(npts))
      ifinout = 1

      thet = hkrand(0)*pi
      phi = hkrand(0)*2*pi
      do i=1,npts
        rhs(i) = 0
        do j=1,ntot
          call h3d_slp(xyz_in(1,j),srcvals(1,i),dpars,zpars,ipars,pot)
          rhs(i) = rhs(i) + charges(j)*pot
        enddo
        sigma(i) = 0
      enddo


      numit = 200
      niter = 0
      allocate(errs(numit+1))

      eps = 0.51d-6


      call cpu_time(t1)
C$      t1 = omp_get_wtime()      
      call helm_comb_dir_solver(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,eps,zpars,numit,ifinout,rhs,niter,errs,
     2  rres,sigma)

      call prinf('niter=*',niter,1)
      call prin2('rres=*',rres,1)
      call prin2('errs=*',errs,niter)

      call cpu_time(t2)
C$       t2 = omp_get_wtime()
      call prin2('analytic solve time=*',t2-t1,1)

      open(unit=33,file='plane-res/sigma-analyticmulti-50l.dat')
      do i=1,npts
        write(33,*) real(sigma(i)),imag(sigma(i))
      enddo
      close(33)

      open(unit=33,file='plane-res/rhs-analyticmulti-50l.dat')
      do i=1,npts
        write(33,*) real(rhs(i)),imag(rhs(i))
      enddo
      close(33)

      nlat = 301
      ntarg = nlat*nlat
      allocate(targs(3,ntarg),ipatch_id(ntarg),uvs_targ(2,ntarg))
      allocate(isout(ntarg))
      xyz_start(1) = -0.5d0
      xyz_start(2) = -6.0d0
      xyz_start(3) = -3.0d0

      dxyz(1) = 0.04d0
      dxyz(2) = 0.04d0
      dxyz(3) = 1.0d0

      do i=1,nlat
        do j=1,nlat
          ipt = (i-1)*nlat+j
          targs(1,ipt) = xyz_start(1) + dxyz(1)*(j-1)
          targs(2,ipt) = xyz_start(2) + dxyz(2)*(i-1) 
          targs(3,ipt) = 0
          ipatch_id(ipt) = -1
          uvs_targ(1,ipt) = 0.0d0
          uvs_targ(2,ipt) = 0.0d0
          isout(ipt) = 0
        enddo
      enddo

      allocate(sigmatmp(npts))
      zpars_tmp(1) = 1.0d-6
      zpars_tmp(2) = 0.0d0
      zpars_tmp(3) = 1.0d0

      allocate(pvals(ntarg),pvalsex(ntarg))
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,npts
        sigmatmp(i) = 1
      enddo
C$OMP END PARALLEL DO      


      call lpcomp_helm_comb_dir(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,3,ntarg,targs,ipatch_id,uvs_targ,eps,zpars_tmp,
     2  sigmatmp,pvalsex)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,pot)
      do i=1,ntarg
        if(abs(pvalsex(i)).le.1.0d-3*4*pi) isout(i) = 1
        pvalsex(i) = 0
        do j=1,ntot
          call h3d_slp(xyz_in(1,j),targs(1,i),dpars,zpars,ipars,pot)
          pvalsex(i) = pvalsex(i) + pot*charges(j)
        enddo
      enddo
C$OMP END PARALLEL DO 

      call lpcomp_helm_comb_dir(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,3,ntarg,targs,ipatch_id,uvs_targ,eps,zpars,
     2  sigma,pvals)
      
      allocate(evals(ntarg),evals2(ntarg))

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(rtmp)
      do i=1,ntarg
        if(isout(i).eq.0) evals(i) = -6
        if(isout(i).eq.0) evals2(i) = -6
        if(isout(i).eq.1) then
          rtmp = abs(pvals(i)-pvalsex(i))
          evals(i) = log(rtmp)/log(10.0d0)
          rtmp = rtmp/abs(pvalsex(i))
          evals2(i) = log(rtmp)/log(10.0d0) 
        endif
      enddo
C$OMP END PARALLEL DO

      allocate(rsigma(npts))
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,npts
        rsigma(i) = real(sigma(i))
      enddo
C$OMP END PARALLEL DO      

      fname = 'plane-res/a380_rsigma_analyticmulti_ext.vtk'
      title = 'a380 - real part'
      call surf_vtk_plot_scalar(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,rsigma,fname,title)

      ndims(1) = nlat
      ndims(2) = nlat
      ndims(3) = 1
      fname = 'plane-res/a380_err_analyticmulti_ext.vtk'
      title = 'abc'
      call vtk_write_plane(ndims,ntarg,xyz_start,dxyz,evals,title,
     1   fname)

      fname = 'plane-res/a380_errrel_analyticmulti_ext.vtk'
      title = 'abc'
      xyz_start(3) = -6
      call vtk_write_plane(ndims,ntarg,xyz_start,dxyz,evals2,title,
     1   fname)


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

   




