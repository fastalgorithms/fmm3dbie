      implicit real *8 (a-h,o-z) 
      implicit integer *8 (i-n)
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:)
      integer *8 ipars(2)

      integer *8, allocatable :: norders(:),ixyzs(:),iptype(:)

      real *8 xyz_out(3),xyz_in(3)
      complex *16, allocatable :: sigma(:),rhs(:)
      complex *16, allocatable :: sigma2(:),rhs2(:)

      real *8, allocatable :: errs(:)
      real *8 thet,phi
	  complex * 16  zpars(5), omega, ep0,mu0,ep1,mu1

      integer *8 numit,niter
      character *200 title,fname,fname1,fname2

      integer *8 ipatch_id
      real *8 uvs_targ(2)

      logical isout0,isout1

      complex *16 pot,potex,ztmp,ima,zk
      complex *16 alpha_rhs

      integer *8 count1

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4

!
! Here we define the helmholtz parameters:    
!
	  
	omega=0.330d0
	ep0=1.00000d0
	mu0=1.0d0
	ep1=(1.1d0,0.0d0)
	mu1=1.200d0

	zpars(1) = omega
	zpars(2) = ep0
	zpars(3) = mu0
	zpars(4) = ep1
	zpars(5) = mu1
	  
!	write (*,*) 'Electrical length for the cube (simplest_cube_quadratic): ', 3.0d0*omega*sqrt(ep1*mu1)/(2*3.1415927d0)


      xyz_out(1) = 50.0d0
      xyz_out(2) = 60.0d0
      xyz_out(3) = 70.0d0

      fname = '../../../../Geometries_go3/' // &
     & 'simplest_cube_quadratic_v4_o04_r01.go3'
      xyz_in(1) = 0.00d0
      xyz_in(2) = 0.0d0
      xyz_in(3) = 1.5d0


!       fname =  '../../../../Geometries_go3/' // &
!      & 'Round_2_o04_r01.go3'
!      xyz_in(1) = 0.00d0
!      xyz_in(2) = 0.0d0
!      xyz_in(3) = 0.0d0

      write (*,*) 'Nombre fichero:', fname
      call open_gov3_geometry_mem(fname,npatches,npts)

      call prinf('npatches=*',npatches,1)
      call prinf('npts=*',npts,1)

      allocate(srcvals(12,npts),srccoefs(9,npts))
      allocate(ixyzs(npatches+1),iptype(npatches),norders(npatches))
      allocate(wts(npts))

      call open_gov3_geometry(fname,npatches,norders,ixyzs,&
     &iptype,npts,srcvals,srccoefs,wts )


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


      allocate(sigma(2*npts),rhs(2*npts))
!      allocate(sigma2(npts),rhs2(npts))
	  
!	  allocate(a_vect(npts),RHS_vect(npts))
      ifinout = 1
	  
      thet = hkrand(0)*pi
      phi = hkrand(0)*2*pi
	  
!     Get right hand side (f and g)

      call get_RHS_transmission(xyz_in,xyz_out,npts,srcvals,zpars,rhs)

!     Start solver (GMRES)

      numit = 400
      niter = 0
      allocate(errs(numit+1))

      eps = 1d-9
	  eps_gmres=1d-14

      call cpu_time(t1)
!C$      t1 = omp_get_wtime()      
      call h_transmission_solver(npatches,norders,ixyzs,iptype,npts,&
     &srccoefs,srcvals,eps,zpars,numit,ifinout,rhs,eps_gmres,niter,errs,&
     &rres,sigma)

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


!
!       test solution at a point in the interior and exterior region
!

		    
		call test_accuracy_transmission(eps,sigma,zpars,npts,wts,srcvals,xyz_in,xyz_out)
	  
      stop
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
!    npatches - integer *8
!       number of patches
!    norder - integer *8
!       order of discretization
!    npts - integer *8
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
      integer *8 npatches,norder,npts,npols
      real *8 srccoefs(9,npts),srcvals(12,npts),xyzout(3),wts(npts)
      real *8 tmp(3)
      real *8 dpars,done,pi
      real *8, allocatable :: rsurf(:),err_p(:,:) 
      integer *8 ipars,norderhead,nd
      complex *16, allocatable :: sigma_coefs(:,:), sigma_vals(:,:)
      complex *16 zk,val

      integer *8 ipatch,j,i
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

   




