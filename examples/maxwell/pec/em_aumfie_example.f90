program em_aumfie_example

  implicit double precision (a-h,o-z) 
  double precision, allocatable :: srcvals(:,:),srccoefs(:,:)
  double precision, allocatable :: wts(:),rsigma(:)
  integer :: ipars(2)

  integer, allocatable :: norders(:),ixyzs(:),iptype(:)

  double precision xyz_out(3),xyz_in(3)
  complex *16, allocatable :: sigma(:),rhs(:)
  complex *16, allocatable :: sigma2(:),rhs2(:)
  
  complex *16, allocatable :: a_vect(:),RHS_vect(:),rhs_nE(:),rhs_nE_aux(:),a_s(:)
  complex *16 :: vf(3), e0(3), kvec(3)
  double complex, allocatable :: ein(:,:), hin(:,:), jvec(:,:)
  double complex :: uvec(3), vvec(3)


  double precision, allocatable :: errs(:)
  complex * 16 zpars(3)
  integer numit,niter
  character *200 title,fname,fname1,fname2

  integer ipatch_id
  double precision uvs_targ(2)

  logical isout0,isout1

  complex *16 pot,potex,ztmp,ima,zk
  complex *16 alpha_rhs

  integer count1

  data ima/(0.0d0,1.0d0)/


  call prini(6,13)

  done = 1
  pi = atan(done)*4

  ! set some parameters, wavelength, wavenumber, etc.
  dlam = 1.0d0
  zk = 2*pi/dlam

  zpars(1) = zk 
  zpars(2) = 1.0d0
  zpars(3) = 0.0d0

  xyz_in(1) = 0.11d0
  xyz_in(2) = 0.0d-5
  xyz_in(3) = 0.37d0

  xyz_out(1) = -3.5d0
  xyz_out(2) = 3.1d0
  xyz_out(3) = 20.1d0
      
  fname = '../../geometries/sphere_192_o03.go3'
  print *
  print *, 'loading file ', trim(fname)
  call open_gov3_geometry_mem(fname, npatches, npts)

  call prinf('npatches =*',npatches,1)
  call prinf('npts =*',npts,1)

  allocate(srcvals(12,npts),srccoefs(9,npts))
  allocate(ixyzs(npatches+1),iptype(npatches),norders(npatches))
  allocate(wts(npts))

  ! open the geometry and plot it
  call open_gov3_geometry(fname,npatches,norders,ixyzs,&
       iptype,npts,srcvals,srccoefs,wts)

  call surf_vtk_plot(npatches, norders, ixyzs, iptype, &
       npts, srccoefs, srcvals, 'aumfie_surface.vtk', 'the scatterer')

  norder = norders(1)
  ! allocate some room for the solution
  allocate(sigma(npts),rhs(npts))
  allocate(sigma2(npts),rhs2(npts))
  allocate(a_vect(2*npts),rhs_vect(2*npts),rhs_nE(npts), &
       rhs_nE_aux(npts),a_s(npts))
  

  print *
  print *, '. . . generating the right hand side for MFIE'

  vf(1)=1.0d0
  vf(2)=2.0d0
  vf(3)=3.0d0

  e0(1) = 1
  e0(2) = 0
  e0(3) = 0
  kvec(1) = 0
  kvec(2) = 1
  kvec(3) = 0
  allocate( ein(3,npts), hin(3,npts) )
  call em_planewave(e0, zk, kvec, npts, srcvals, ein, hin)







  call get_rhs_em_mfie_pec(xyz_out, vf, zpars(2), npts, srcvals, &
       zpars(1), rhs_vect)

    
  ! set some parameters and call the solver
  numit = 400
  niter = 0
  allocate(errs(numit+100))

  eps = 1d-3
  eps_gmres=1d-10
  ifinout = 1
  
  call cpu_time(t1)
  !$  t1 = omp_get_wtime()      
  call em_mfie_solver(npatches, norders, ixyzs, iptype, npts, &
       srccoefs, srcvals, eps, zpars, numit, ifinout, rhs_vect, &
       eps_gmres, &
       niter, errs, rres, a_vect, rhs_nE)

  print *
  call prinf('from mfie solve, gmres niter=*',niter,1)
  call prin2('from mfie solve, rres=*',rres,1)
  call prin2('from mfie solve, errs=*',errs,niter)

  ! convert the solution a_vect into regular cartesian 3-vectors and plot it
  allocate( jvec(3,npts) )
  do i = 1,npts
     uvec(1)
     jvec(1,i) = a_vect(i)*srcvals(4,i)
  end do


  
  ! estimate the accuracy in the mfie solve
  call test_accuracy_em_mfie_pec(eps, a_vect, zpars, npts, wts, &
       srcvals, xyz_in, vf, xyz_out)

  
  call get_rhs_em_aumfie_pec(xyz_out,vf,zpars(2),npts,srcvals,&
       zpars(1), rhs_nE_aux)
  
  do count1=1,npts
     rhs_nE(count1)=rhs_nE(count1)+rhs_nE_aux(count1)
  enddo

  
  call em_aumfie_solver(npatches,norders,ixyzs,&
       &iptype,npts,srccoefs,srcvals,eps,zpars,numit,ifinout,&
       &rhs_nE,eps_gmres,niter,errs,rres,a_s)

  call cpu_time(t2)
  !$       t2 = omp_get_wtime()
  print *
  call prin2('analytic solve time=*',t2-t1,1)


  !
  ! Test accuracy of solution
  !
  call test_accuracy_em_aumfie_pec(eps,a_vect,a_s,zpars,npts,wts,&
       srcvals,xyz_in,vf,xyz_out)

      stop
end program em_aumfie_example





subroutine em_planewave(e0, zk, kvec, npts, srcvals, efield, hfield)
  implicit double precision (a-h,o-z)
  !
  ! Author: Mike O'Neil
  ! Email:  moneil@flatironinstitute.org
  ! Date:   May 19, 2023
  !
  ! This subroutine evaluates the electromagnetic linearly polarized planewave
  ! described by the electric polarization vector e0 and wavevector kvec. In
  ! order to satisfy Maxwell's equations:
  !
  !   curl(E) = ikH, curl(H) = -ikE
  !    div(E) = 0,    div(H) = 0
  !
  ! it must be the case that magnetic polarization vector h0 = kvec X e0
  !
  ! Input:
  !   e0 - the electric field linear polarization direction
  !   zk - the complex valued wavenumber
  !   kvec - the wave propagation direction, up in the exponent
  !   npts - number of points at which to evaluate the planewave
  !   srcvals(12,npts) - the points at which to evaluate the planewave,
  !       it's assumed that the first 3 components contain the x,y,z coordinates
  !
  ! Output:
  !   efield - the electric field at all the srcvals
  !   hfield - the magnetic field at all the srcvals
  !

  ! calling sequence variables
  complex *16 :: e0(3), zk, kvec(3), efield(3,npts), hfield(3,npts)
  double precision :: srcvals(12,npts)

  ! local variables
  complex *16 :: h0(3), ima, zp, cd

  ! check that e0 is orthgonal to kvec first, and that kvec has length 1,
  ! otherwise something will go wrong
  dd = 0
  dd = e0(1)*dconjg(kvec(1)) + e0(2)*dconjg(kvec(2)) + e0(3)*dconjg(kvec(3))
  dnorm = abs(e0(1))**2 + abs(e0(2))**2 + abs(e0(3))**2
  dnorm = sqrt(dnorm)
  if (abs(dd)/dnorm .gt. 1.0d-12) then
     call prin2('inner product in em_planewave= *', dd, 1)
     stop
  end if

  ! compute the magnetic polarization vector, h0 = kvec X e0
  h0(1) = kvec(2)*e0(3) - kvec(3)*e0(2)
  h0(2) = kvec(3)*e0(1) - kvec(1)*e0(3)
  h0(3) = kvec(1)*e0(2) - kvec(2)*e0(1)

  ima = (0,1)
  do i = 1,npts
     zp = kvec(1)*srcvals(1,i) + kvec(2)*srcvals(2,i) + kvec(3)*srcvals(3,i)
     cd = exp(ima*zk*zp)
     efield(1,i) = e0(1)*cd
     efield(2,i) = e0(2)*cd
     efield(3,i) = e0(3)*cd
     hfield(1,i) = h0(1)*cd
     hfield(2,i) = h0(2)*cd
     hfield(3,i) = h0(3)*cd
  end do

  return
end subroutine em_planewave






      subroutine setup_geom(igeomtype,norder,npatches,ipars,& 
     &srcvals,srccoefs,ifplot,fname)
      implicit double precision (a-h,o-z)
      integer igeomtype,norder,npatches,ipars(*),ifplot
      character (len=*) fname
      double precision srcvals(12,*), srccoefs(9,*)
      double precision, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts(:)

      double precision, pointer :: ptr1,ptr2,ptr3,ptr4
      integer, pointer :: iptr1,iptr2,iptr3,iptr4
      double precision, target :: p1(10),p2(10),p3(10),p4(10)
      double precision, allocatable, target :: triaskel(:,:,:)
      double precision, allocatable, target :: deltas(:,:)
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
!    srccoefs - double precision (9,npts)
!       koornwinder expansion coefficients of geometry info
!    xyzout -  double precision (3)
!       point to be tested
!
!  output: 
!    isout - boolean
!      whether the target is in the interior or not
!

      implicit none
      integer npatches,norder,npts,npols
      double precision srccoefs(9,npts),srcvals(12,npts),xyzout(3),wts(npts)
      double precision tmp(3)
      double precision dpars,done,pi
      double precision, allocatable :: rsurf(:),err_p(:,:) 
      integer ipars,norderhead,nd
      complex *16, allocatable :: sigma_coefs(:,:), sigma_vals(:,:)
      complex *16 zk,val

      integer ipatch,j,i
      double precision ra,ds
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

   




