      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:)
      integer ipars(2),ndims(3)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)
      integer, allocatable :: isout(:)

      real *8, allocatable :: xyz_in(:,:),cms(:,:),rads(:)
      real *8, allocatable :: rad_near(:)

      complex *16, allocatable :: sigma(:),rhs(:),charges(:)
      complex *16, allocatable :: rhs_coefs(:)
      complex *16, allocatable :: sigmatmp(:)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: targs(:,:),uvs_targ(:,:),pdis(:)

      complex *16, allocatable :: pvals(:),pvalsex(:)
      real *8, allocatable :: evals(:),evals2(:)
      real *8, allocatable :: errs(:)
      real *8, allocatable :: errp_dens(:,:),rsurf(:),errp_surf(:,:)
      real *8 errm_surf(100),errm_dens(100),eps_gmres
      real *8 thet,phi
      real *8 xyz_start(3),dxyz(3)
      complex * 16 zpars(3),zpars_tmp(3)
      integer numit,niter
      logical isin(1000)
      character *100 title,fname



      complex *16 pot,potex,ztmp,ima

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4

c
c   simulation for plane 50 wavelengths in size
c
      zk = 28.56d0
      zpars(1) = zk 
      zpars(2) = -ima*zk
      zpars(3) = 1.0d0

      
      fname = '../../geometries/A380_Final_o03_r02.go3'
      
      call open_gov3_geometry_mem(fname,npatches,npts)

      call prinf('npatches=*',npatches,1)
      call prinf('npts=*',npts,1)

      allocate(srcvals(12,npts),srccoefs(9,npts))
      allocate(ixyzs(npatches+1),iptype(npatches),norders(npatches))
      allocate(wts(npts))

      call open_gov3_geometry(fname,npatches,norders,ixyzs,
     1   iptype,npts,srcvals,srccoefs,wts)

      nmax = 2000
      allocate(xyz_in(3,nmax),charges(nmax))
      open(unit=33,file='fuselage-sources-8.txt')
      read(33,*) nfuse

      ntail = 0
      open(unit=34,file='tail-sources-3.5.txt')
      read(34,*) ntail

      ntot = nfuse + ntail
      do i=1,nfuse
        read(33,*) xyz_in(1,i),xyz_in(2,i),xyz_in(3,i)
      enddo

      do i=1,ntail
        ii = nfuse+i
        read(34,*) xyz_in(1,ii),xyz_in(2,ii),xyz_in(3,ii)
      enddo
      close(33)
      close(34)

      do i=1,ntot
        charges(i) = (hkrand(0)-0.5d0) + ima*(hkrand(0)-0.5d0)
      enddo


      call prin2('xyz_in=*',xyz_in,3*ntot)
      call prin2('charges=*',charges,2*ntot)
      open(unit=33,file='plane-res/chargerhs-analyticmulti-50l-6.dat')
      write(33,*) ntot
      do i=1,ntot
        write(33,*) xyz_in(1,i),xyz_in(2,i),xyz_in(3,i),
     1     real(charges(i)),imag(charges(i))
      enddo
      close(33)

      do i=1,ntot
        call test_exterior_pt(npatches,norders,npts,srcvals,srccoefs,
     1  wts,xyz_in(1,i),isin(i))
      enddo

      allocate(sigma(npts),rhs(npts))
      ifinout = 1

      thet = hkrand(0)*pi
      phi = hkrand(0)*2*pi
      rhs_m = 0
      do i=1,npts
        rhs(i) = 0
        do j=1,ntot
          call h3d_slp(xyz_in(1,j),srcvals(1,i),dpars,zpars,ipars,pot)
          rhs(i) = rhs(i) + charges(j)*pot
        enddo
        if(abs(rhs(i)).gt.rhs_m) rhs_m = abs(rhs(i))
        sigma(i) = 0
      enddo

      allocate(errp_dens(2,npatches),errp_surf(9,npatches))

      allocate(rsurf(npatches))
      rmax = 0
      rmin = 100.0d0
      imax = 1
      imin = 1
      do i=1,npatches
        rsurf(i) = 0
        istart = ixyzs(i)
        npols = ixyzs(i+1)-ixyzs(i)
        do j=1,npols
          jpt = istart+j-1
          rsurf(i) = rsurf(i) + wts(jpt)
        enddo
        if(rsurf(i).gt.rmax) then
          rmax = rsurf(i)
          imax = i
        endif
        if(rsurf(i).lt.rmin) then
          rmin = rsurf(i)
          imin = i
        endif
      enddo

      print *, "rmax,rmin=",rmax,rmin
      print *, "imax,imin=",imax,imin

      call prin2('rsurf=*',rsurf,24)
      allocate(rhs_coefs(npts),pdis(npatches))
      call surf_vals_to_coefs(2,npatches,norders,ixyzs,iptype,npts,
     1  rhs,rhs_coefs)
      
      call surf_fun_error(2,npatches,norders,ixyzs,iptype,npts,rsurf,
     1  rhs_coefs,errp_dens,errm_dens)
      
      call surf_fun_error(9,npatches,norders,ixyzs,iptype,npts,rsurf,
     1  srccoefs,errp_surf,errm_surf)

      call get_patch_distortion(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,wts,pdis)

      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))
      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts, 
     1     srccoefs,cms,rads)
      do i=1,npatches
        rad_near(i) = 3.5d0*rads(i)
      enddo
      print *, ntot
      call findnearmem(cms,npatches,rad_near,xyz_in,ntot,nnz)
      print *, nnz

      do i=1,npatches
        errp0 = maxval(errp_dens(:,i))
        errs0 = maxval(errp_surf(:,i))
        write(77,*) i,rsurf(i),rads(i),pdis(i),errs0,errp0
      enddo
      call prin2('errm_dens=*',errm_dens,2)
      print *, "rhs_m=",rhs_m
      


      numit = 150
      niter = 0
      allocate(errs(numit+1))

      eps = 0.51d-6
      eps_gmres = eps


      call cpu_time(t1)
C$      t1 = omp_get_wtime()      
      call helm_comb_dir_solver(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,eps,zpars,numit,ifinout,rhs,eps_gmres,
     2  niter,errs,rres,sigma)

      call prinf('niter=*',niter,1)
      call prin2('rres=*',rres,1)
      call prin2('errs=*',errs,niter)

      call cpu_time(t2)
C$       t2 = omp_get_wtime()
      call prin2('analytic solve time=*',t2-t1,1)

      open(unit=33,file='plane-res/sigma-analyticmulti-50l-6.dat')
      do i=1,npts
        write(33,*) real(sigma(i)),imag(sigma(i))
      enddo
      close(33)

      open(unit=33,file='plane-res/rhs-analyticmulti-50l-6.dat')
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
ccC$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,npts
        sigmatmp(i) = 1
      enddo
ccC$OMP END PARALLEL DO      


      call lpcomp_helm_comb_dir(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,3,ntarg,targs,ipatch_id,uvs_targ,eps,zpars_tmp,
     2  sigmatmp,pvalsex)

      do i=1,ntarg
        if(abs(pvalsex(i)).le.4*pi*1.0d-3) isout(i) = 1
      enddo

      print *, "ntot=",ntot
      call prin2('xyz_in=*',xyz_in,3*ntot)
      call prin2('charges=*',charges,2*ntot)

      do i=1,ntarg
        pvalsex(i) = 0
        do j=1,ntot
          pot = 0
          call h3d_slp(xyz_in(1,j),targs(1,i),dpars,zpars,ipars,pot)
          pvalsex(i) = pvalsex(i) + pot*charges(j)
        enddo
      enddo

      call lpcomp_helm_comb_dir(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,3,ntarg,targs,ipatch_id,uvs_targ,eps,zpars,
     2  sigma,pvals)
      
      allocate(evals(ntarg),evals2(ntarg))
      emax = -100
      emax2 = -100

      rmax = 0
      do i=1,ntarg
        if(isout(i).eq.1) then
          if(abs(pvalsex(i)).gt.rmax) rmax = abs(pvalsex(i))
        endif
      enddo

      open(unit=33,file='plane-res/pottarg-analyticmulti-50l-6.dat')
ccC$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(rtmp)
      do i=1,ntarg
        if(isout(i).eq.0) evals(i) = -6
        if(isout(i).eq.0) evals2(i) = -6
        if(isout(i).eq.1) then
          rtmp = abs(pvals(i)-pvalsex(i))
          evals(i) = log(rtmp)/log(10.0d0)
          rtmp = rtmp/rmax
          evals2(i) = log(rtmp)/log(10.0d0) 
          if(evals(i).gt.emax) emax = evals(i)
          if(evals2(i).gt.emax2) emax2 = evals2(i)
        endif
        write(33,*) real(pvals(i)),imag(pvals(i)),real(pvalsex(i)),
     1    imag(pvalsex(i)),evals(i),evals2(i)
      enddo
ccC$OMP END PARALLEL DO
      close(33)

      print *, "max errors=",emax,emax2
      stop

      allocate(rsigma(npts))
ccC$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,npts
        rsigma(i) = real(sigma(i))
      enddo
ccC$OMP END PARALLEL DO      

      fname = 'plane-res/a380_rsigma_analytic_halfl_ext.vtk'
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

      fname = 'plane-res/a380_errrel_analytic1_halfl_ext.vtk'
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

   




