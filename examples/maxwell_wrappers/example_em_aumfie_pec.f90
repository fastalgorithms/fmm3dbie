      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:)
      integer ipars(2)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)

      real *8 xyz_out(3),xyz_in(3)
      complex *16, allocatable :: sigma(:),rhs(:)
      complex *16, allocatable :: sigma2(:),rhs2(:)

      complex ( kind = 8 ), allocatable :: a_vect(:),RHS_vect(:),rhs_nE(:),rhs_nE_aux(:),a_s(:)
      complex *16 vf(3)


      real *8, allocatable :: errs(:)
      real *8 thet,phi
      complex * 16 zpars(3)
      integer numit,niter
      character *200 title,fname,fname1,fname2

      integer ipatch_id
      real *8 uvs_targ(2)

      logical isout0,isout1

      complex *16 pot,potex,ztmp,ima,zk
      complex *16 alpha_rhs

      integer count1

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4

!
!   simulation for plane 50 wavelengths in size
!
      zk = 2.0d-1
	  zk=8.d-1
      zpars(1) = zk 
      zpars(2) = 1.0d0
      zpars(3) = 0.0d0

      
      xyz_in(1) = 0.00d0
      xyz_in(2) = 0.0d0
      xyz_in(3) = 1.5d0

      xyz_out(1) = 50.0d0
      xyz_out(2) = 60.0d0
      xyz_out(3) = 70.0d0
	  
!!      fname = '../../geometries/A380_Final_o03_r12.go3'
      fname = '../../../../Geometries_go3/' // &
     & 'simplest_cube_quadratic_v4_o04_r01.go3'
            
      call open_gov3_geometry_mem(fname,npatches,npts)

      call prinf('npatches=*',npatches,1)
      call prinf('npts=*',npts,1)

      allocate(srcvals(12,npts),srccoefs(9,npts))
      allocate(ixyzs(npatches+1),iptype(npatches),norders(npatches))
      allocate(wts(npts))

      call open_gov3_geometry(fname,npatches,norders,ixyzs,&
     &iptype,npts,srcvals,srccoefs,wts)


!      fname = 'plane-res/a380.vtk'
!      title = 'a380'
!      call surf_vtk_plot(npatches,norders,ixyzs,iptype,npts,
!     1   srccoefs,srcvals,fname,title)

      norder = norders(1)
!      call test_exterior_pt(npatches,norder,npts,srcvals,srccoefs,&
!     &wts,xyz_in,isout0)
      
!      call test_exterior_pt(npatches,norder,npts,srcvals,srccoefs,&
!     &wts,xyz_out,isout1)

       print *, isout0,isout1 


      allocate(sigma(npts),rhs(npts))
      allocate(sigma2(npts),rhs2(npts))
	  
	  allocate(a_vect(2*npts),RHS_vect(2*npts),rhs_nE(npts),rhs_nE_aux(npts),a_s(npts))
      ifinout = 1
	  
      thet = hkrand(0)*pi
      phi = hkrand(0)*2*pi
	  
	  vf(1)=1.0d0
	  vf(2)=2.0d0
	  vf(3)=3.0d0
!	  
!	  call get_RHS_EMFIE(xyz_in,vf,npts,srcvals,zk,RHS_vect)
	  call get_rhs_em_mfie_pec(xyz_out,vf,zpars(2),npts,srcvals,zpars(1),RHS_vect)

!	  do count1=1,10
!	      write (*,*) count1,RHS_vect(count1),RHS_vect(npts+count1)
!	  enddo
	  
!      do i=1,npts
!        if(ifinout.eq.0)& 
!     &call h3d_slp(xyz_out,srcvals(1,i),dpars,zpars,ipars,rhs(i))
!        if(ifinout.eq.1)& 
!     &call h3d_slp(xyz_in,srcvals(1,i),dpars,zpars,ipars,rhs(i))
!        rhs(i) = rhs(i)
!        x = srcvals(1,i)*cos(thet)
!        y = srcvals(2,i)*sin(thet)*cos(phi)
!        z = srcvals(3,i)*sin(thet)*sin(phi)
!        rhs2(i) = exp(ima*zk*(x+y+z))
!        sigma(i) = 0
!        sigma2(i) = 0
!      enddo


      numit = 400
      niter = 0
      allocate(errs(numit+1))

      eps = 1d-6
	  eps_gmres=1d-10

      call cpu_time(t1)
!C$      t1 = omp_get_wtime()      
      call em_mfie_solver(npatches,norders,ixyzs,iptype,npts,&
     &srccoefs,srcvals,eps,zpars,numit,ifinout,RHS_vect,eps_gmres,niter,errs,&
     &rres,a_vect,rhs_nE)
	 
	call test_accuracy_em_mfie_pec(eps,a_vect,zpars,npts,wts,srcvals,xyz_in,vf,xyz_out)

      call get_rhs_em_aumfie_pec(xyz_out,vf,zpars(2),npts,srcvals,zpars(1),rhs_nE_aux)
	  
	  do count1=1,npts
        rhs_nE(count1)=rhs_nE(count1)+rhs_nE_aux(count1)
      enddo
      call prinf('niter=*',niter,1)
      call prin2('rres=*',rres,1)
      call prin2('errs=*',errs,niter)

      call cpu_time(t2)
!C$       t2 = omp_get_wtime()
      call prin2('analytic solve time=*',t2-t1,1)

!      open(unit=33,file='plane-res/sigma-analytic-50l.dat')
!      do i=1,npts
!        write(33,*) real(sigma(i)),imag(sigma(i))
!      enddo
!      close(33)


      call em_aumfie_solver(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,eps,zpars,numit,ifinout,&
     &rhs_nE,eps_gmres,niter,errs,rres,a_s)

call test_accuracy_em_aumfie_pec(eps,a_vect,a_s,zpars,npts,wts,srcvals,xyz_in,vf,xyz_out)

!
!       test solution at interior point
!
      
	  

      stop
      end




      subroutine setup_geom(igeomtype,norder,npatches,ipars,& 
     &srcvals,srccoefs,ifplot,fname)
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
        call xtri_platonic(itype, ipars(1), npmax, ntri,& 
     &triaskel, isides)

        xtri_geometry => xtri_sphere_eval
        ptr1 => triaskel(1,1,1)
        ptr2 => p2(1)
        ptr3 => p3(1)
        ptr4 => p4(1)


        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2,& 
     &ptr3,ptr4, norder,'Triangulated surface of the sphere')
        endif


        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4,&
     &npols,uvs,umatr,srcvals,srccoefs)
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
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2),&
     &nover,npatches,npatches,triaskel)

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
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2,& 
     &iptr3,iptr4, norder,'Triangulated surface of the stellarator')
        endif

        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,iptr3,iptr4,&
     &npols,uvs,umatr,srcvals,srccoefs)
      endif
      
      return  
      end


      subroutine test_exterior_pt(npatches,norder,npts,srcvals,&
     &srccoefs,wts,xyzout,isout)
!
!
!  this subroutine tests whether the pt xyzin, is
!  in the exterior of a surface, and also estimates the error
!  in representing e^{ir/2}/r and \grad e^{ir/2}/r \cdot n
!  centered at the interior point. Whether a point 
!  is in the interior or not is tested using Gauss' 
!  identity for the flux due to a point charge
!
!
!  input:
!    npatches - integer
!       number of patches
!    norder - integer
!       order of discretization
!    npts - integer
!       total number of discretization points on the surface
!    srccoefs - real *8 (9,npts)
!       koornwinder expansion coefficients of geometry info
!    xyzout -  real *8 (3)
!       point to be tested
!
!  output: 
!    isout - boolean
!      whether the target is in the interior or not
!

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

   




