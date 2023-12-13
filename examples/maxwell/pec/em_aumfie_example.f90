program em_aumfie_example

  implicit double precision (a-h,o-z) 
  implicit integer(8) (i-n)
  double precision, allocatable :: srcvals(:,:),srccoefs(:,:)
  double precision, allocatable :: wts(:),rsigma(:)
  integer(8) :: ipars(2)

  integer(8), allocatable :: norders(:),ixyzs(:),iptype(:)

  double precision xyz_out(3),xyz_in(3)
  complex *16, allocatable :: sigma(:),rhs(:)

  complex *16, allocatable :: a_vect(:),RHS_vect(:),rhs_nE(:),rhs_nE_aux(:),a_s(:)
  complex *16 :: vf(3), e0(3), kvec(3)
  double complex, allocatable :: ein(:,:), hin(:,:), jvec(:,:)
  double complex :: uvec(3), vvec(3)


  double precision, allocatable :: errs(:)
  complex * 16 :: zpars(3)
  integer(8) :: numit,niter
  character *200 :: title,fname,fname1,fname2

  integer(8) :: ipatch_id
  double precision :: uvs_targ(2), ru(3), rv(3)

  logical isout0,isout1

  complex *16 pot,potex,ztmp,ima,zk
  complex *16 alpha_rhs

  integer(8) count1

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
  allocate( rhs(2*npts) )
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


  print *, 'generating incoming field . . .'
  allocate( ein(3,npts), hin(3,npts) )

  call em_field_electric_dipole(zk, xyz_out, vf, npts, srcvals(1:3,:), &
       ein, hin, int(0,8))
  call em_field_magnetic_dipole(zk, xyz_out, vf, npts, srcvals, ein, hin, int(1,8))

  ! now compute - n \times H
  !do j = 1,npts
  !   call orthonormalize(srcvals(4:6,j), srcvals(10:12,j), ru, rv)
  !   call dzdot_prod3d(rv, hin(:,j), rhs_vect(j))
  !   rhs_vect(j) = -rhs_vect(j)
  !   call dzdot_prod3d(ru, hin(:,j), rhs_vect(npts+j))
  !enddo

  !call prin2('and again, rhs_vect = *', rhs_vect, 30)

    
  ! set some parameters and call the solver
  numit = 400
  niter = 0
  allocate(errs(numit+100))

  eps = 1d-3
  eps_gmres=1d-8
  ifinout = 1
  
  call cpu_time(t1)
  !$  t1 = omp_get_wtime()      
  call em_mfie_solver(npatches, norders, ixyzs, iptype, npts, &
       srccoefs, srcvals, eps, zpars, numit, ifinout, &
       ein, hin, rhs_vect, &
       eps_gmres, &
       niter, errs, rres, a_vect, rhs_nE)

  print *
  call prinf('from mfie solve, gmres niter=*',niter,1)
  call prin2('from mfie solve, rres=*',rres,1)
  call prin2('from mfie solve, errs=*',errs,niter)

  ! convert the solution a_vect into regular cartesian 3-vectors and plot it
  !allocate( jvec(3,npts) )
  !do i = 1,npts
  !   !call orthonormalize(srcvals(4,i),srcvals(), uvec, vvec)
  !   jvec(1,i) = a_vect(i)*uvec(1) + a_vect(npts+i)*vvec(1)
  !   jvec(2,i) = a_vect(i)*uvec(2) + a_vect(npts+i)*vvec(2)
  !   jvec(3,i) = a_vect(i)*uvec(3) + a_vect(npts+i)*vvec(3)
  !end do


  
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





subroutine get_rhs_em_mfie_pec777(p0, vf, alpha, ns, srcvals, zk,rhs)
  implicit none
  !
  !  This function obtains the right hand side for the MFIE
  !  formulation for the integral boundary equation:
  !
  !         J/2 - M_{k}[J] = nxH_inc
  !
  !  input:
  !    P0 - real * 8 (3)
  !      location of the source point at the exterior region
  !      WARNING! notice that this formulation uses a representation theorem
  !      for the incoming field in the interior region (MFIE) therefore
  !      therefore it only works for incoming fields generated by sources in
  !      the exterior region (or at infinity like plane waves)
  !
  !    vf - complex *16(3)
  !      Orientation of the magnetic and electric dipoles located at P0 
  !
  !    alpha - complex *16
  !      parameter in the combined formulation
  !   
  !    ns - integer
  !      total number of points on the surface
  !
  !    srcvals - real *8(12,ns)
  !      xyz(u,v) and derivative info sampled at the 
  !      discretization nodes on the surface
  !      srcvals(1:3,i) - xyz info
  !      srcvals(4:6,i) - dxyz/du info
  !      srcvals(7:9,i) - dxyz/dv info
  !
  !    zk - complex *16
  !      Helmholtz parameter 
  !
  !    output:
  !      RHS - complex  *16(2*ns)
  !        right hand side
  !          RHS(1:ns) - first component of  nxH_inc along
  !          the srcvals(4:6,i) direction
  !          RHS(ns+1:2*ns) - second component of nxH_inc along
  !          the (srcvals(10:12,i) x srcvals(4:6,i)) direction
  !
  
  !List of calling arguments
  integer ( kind = 8 ), intent(in) :: ns
  real ( kind = 8 ), intent(in) :: P0(3)
  complex ( kind = 8 ), intent(in) :: vf(3)
  real ( kind = 8 ), intent(in) :: srcvals(12,ns)
  complex ( kind = 8 ), intent(in) :: zk,alpha
  complex ( kind = 8 ), intent(out) :: RHS(2*ns)

  !List of local variables
  complex ( kind = 8 ), allocatable :: efield(:,:), hfield(:,:)
  integer(8) :: count1, lda
  real ( kind = 8 ) :: ru(3),rv(3),cross_aux(3), dnorm
  
  allocate(efield(3,ns), hfield(3,ns))
  
  call em_field_electric_dipole(zk, P0, vf, ns, srcvals(1:3,:), &
       efield, hfield, int(0,8))
  call em_field_magnetic_dipole(zk, P0, vf, ns, srcvals, efield, hfield, int(1,8))

  ! now compute - n \times H
  do count1=1,ns

     call orthonormalize(srcvals(4:6,count1),srcvals(10:12,count1),ru,rv)

     call dzdot_prod3d(rv, hfield(:,count1), rhs(count1))
     rhs(count1) = -rhs(count1)

     call dzdot_prod3d(ru, hfield(:,count1), rhs(ns+count1))

  enddo




  !call prin2('rhs = *', rhs, 20)
  !stop

  return
end subroutine get_rhs_em_mfie_pec777






subroutine em_field_electric_dipole(zk, p0, dipvec, n, points, &
     efield, hfield, init)
  implicit none
  !
  ! Generate the Maxwell field due to an electric dipole, i.e. the fields are
  ! generated as
  !
  !     H = curl(A)
  !     E = 1/(-ima*zk) * curl(curl(A))
  !
  ! where  A = dipvec*exp(ima*zk*r)/r
  ! and r = (points - P0)
  !

  ! Input:
  !   zk - the wavenumber
  !   p0 - location of the dipole
  !   dipvec - dipole direction, no assumption of magnitude, complex-valued
  !   n - number of targets at which to evaluate field
  !   points - evaluation locations, assumed to be dimension(3,n)
  !   init - initialization switch,
  !       0: set E,H to 0
  !       1: increment current E,H values
  !
  ! Output:
  !   E, H - the electric and magnetic ields
  !

  ! List of calling arguments
  integer(8), intent(in) :: n, init
  double precision, intent(in) :: P0(3), points(3,n)
  complex *16, intent(in) :: dipvec(3), zk
  complex *16, intent(out) :: efield(3,n), hfield(3,n)

  ! List of local variables
  double precision :: dx,dy,dz,r
  complex *16 :: R1,R2,au1,au2,ima, cd
  integer(8) :: i

  ima = (0,1)

  if (init .eq. 0) then
     do i = 1,n
        efield(1,i)=0
        efield(2,i)=0
        efield(3,i)=0
        hfield(1,i)=0
        hfield(2,i)=0
        hfield(3,i)=0
     end do
  end if
  
  do i=1,n

     dx=points(1,i)-P0(1)
     dy=points(2,i)-P0(2)
     dz=points(3,i)-P0(3)
     r=sqrt(dx**2+dy**2+dz**2)

     R1=(ima*zk/r**2-1/r**3)
     R2=((ima*zk)**2/r**3-3*ima*zk/r**4+3/r**5)

     cd = exp(ima*zk*r)
     au1=(zk**2/r+R1)*cd/(-ima*zk)
     au2=dx*dipvec(1)+dy*dipvec(2)
     au2=au2+dz*dipvec(3)
     au2=au2*R2*cd/(-ima*zk)

     efield(1,i) = efield(1,i)+(dipvec(1)*au1+dx*au2)
     efield(2,i) = efield(2,i)+(dipvec(2)*au1+dy*au2)
     efield(3,i) = efield(3,i)+(dipvec(3)*au1+dz*au2)

     hfield(1,i) = hfield(1,i)+(dy*dipvec(3)-dz*dipvec(2))*R1*cd
     hfield(2,i) = hfield(2,i)+(dz*dipvec(1)-dx*dipvec(3))*R1*cd
     hfield(3,i) = hfield(3,i)+(dx*dipvec(2)-dy*dipvec(1))*R1*cd
  enddo

  return
end subroutine em_field_electric_dipole






subroutine em_field_magnetic_dipole(zk,P0, vf, n, points, E, H, initial)
implicit none
! a small magnetic dipole; that is: F=vf*exp(ima*zk*r)/r
! E=curlF
! H=1/(ima*zk)*curlcurlF

	!List of calling arguments
	integer(8), intent(in) :: n,initial
	real ( kind = 8 ), intent(in) :: P0(3),points(12,n)
	complex ( kind = 8 ), intent(in) :: vf(3),zk
	complex ( kind = 8 ), intent(out) :: E(3,n),H(3,n)

	!List of local variables
	real ( kind = 8 ) dx,dy,dz,r
	complex ( kind = 8 ) R1,R2,au1,au2,ima
	integer(8) i

	ima = (0.0d0,1.0d0)

	do i=1,n
		if (initial .eq. 0) then
			E(1,i)=0
			E(2,i)=0
			E(3,i)=0
			H(1,i)=0
			H(2,i)=0
			H(3,i)=0
       endif
		dx=points(1,i)-P0(1)
		dy=points(2,i)-P0(2)
		dz=points(3,i)-P0(3)
		r=sqrt(dx**2+dy**2+dz**2)
		R1=(ima*zk/r**2-1/r**3)
		R2=((ima*zk)**2/r**3-3*ima*zk/r**4+3/r**5)
		au1=(zk**2/r+R1)*exp(ima*zk*r)/(ima*zk)
		au2=dx*vf(1)+dy*vf(2)
		au2=au2+dz*vf(3)
		au2=au2*R2*exp(ima*zk*r)/(ima*zk)
		H(1,i)=H(1,i)+(vf(1)*au1+dx*au2)
		H(2,i)=H(2,i)+(vf(2)*au1+dy*au2)
		H(3,i)=H(3,i)+(vf(3)*au1+dz*au2)
		E(1,i)=E(1,i)+(dy*vf(3)-dz*vf(2))*R1*exp(ima*zk*r)
		E(2,i)=E(2,i)+(dz*vf(1)-dx*vf(3))*R1*exp(ima*zk*r)
		E(3,i)=E(3,i)+(dx*vf(2)-dy*vf(1))*R1*exp(ima*zk*r)
	enddo

return
end subroutine em_field_magnetic_dipole










subroutine em_planewave(e0, zk, kvec, npts, srcvals, efield, hfield)
  implicit double precision (a-h,o-z)
  implicit integer(8) (i-n)
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








!       subroutine setup_geom(igeomtype,norder,npatches,ipars,& 
!      &srcvals,srccoefs,ifplot,fname)
!       implicit double precision (a-h,o-z)
!       integer igeomtype,norder,npatches,ipars(*),ifplot
!       character (len=*) fname
!       double precision srcvals(12,*), srccoefs(9,*)
!       double precision, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts(:)

!       double precision, pointer :: ptr1,ptr2,ptr3,ptr4
!       integer, pointer :: iptr1,iptr2,iptr3,iptr4
!       double precision, target :: p1(10),p2(10),p3(10),p4(10)
!       double precision, allocatable, target :: triaskel(:,:,:)
!       double precision, allocatable, target :: deltas(:,:)
!       integer, allocatable :: isides(:)
!       integer, target :: nmax,mmax

!       procedure (), pointer :: xtri_geometry


!       external xtri_stell_eval,xtri_sphere_eval
      
!       npols = (norder+1)*(norder+2)/2
!       allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols))
!       allocate(wts(npols))

!       call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)

!       if(igeomtype.eq.1) then
!         itype = 2
!         allocate(triaskel(3,3,npatches))
!         allocate(isides(npatches))
!         npmax = npatches
!         ntri = 0
!         call xtri_platonic(itype, ipars(1), npmax, ntri,& 
!      &triaskel, isides)

!         xtri_geometry => xtri_sphere_eval
!         ptr1 => triaskel(1,1,1)
!         ptr2 => p2(1)
!         ptr3 => p3(1)
!         ptr4 => p4(1)


!         if(ifplot.eq.1) then
!            call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2,& 
!      &ptr3,ptr4, norder,'Triangulated surface of the sphere')
!         endif


!         call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4,&
!      &npols,uvs,umatr,srcvals,srccoefs)
!       endif

!       if(igeomtype.eq.2) then
!         done = 1
!         pi = atan(done)*4
!         umin = 0
!         umax = 2*pi
!         vmin = 0
!         vmax = 2*pi
!         allocate(triaskel(3,3,npatches))
!         nover = 0
!         call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2),&
!      &nover,npatches,npatches,triaskel)

!         mmax = 2
!         nmax = 1
!         xtri_geometry => xtri_stell_eval

!         allocate(deltas(-1:mmax,-1:nmax))
!         deltas(-1,-1) = 0.17d0
!         deltas(0,-1) = 0
!         deltas(1,-1) = 0
!         deltas(2,-1) = 0

!         deltas(-1,0) = 0.11d0
!         deltas(0,0) = 1
!         deltas(1,0) = 4.5d0
!         deltas(2,0) = -0.25d0

!         deltas(-1,1) = 0
!         deltas(0,1) = 0.07d0
!         deltas(1,1) = 0
!         deltas(2,1) = -0.45d0

!         ptr1 => triaskel(1,1,1)
!         ptr2 => deltas(-1,-1)
!         iptr3 => mmax
!         iptr4 => nmax

!         if(ifplot.eq.1) then
!            call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2,& 
!      &iptr3,iptr4, norder,'Triangulated surface of the stellarator')
!         endif

!         call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,iptr3,iptr4,&
!      &npols,uvs,umatr,srcvals,srccoefs)
!       endif
      
!       return  
!       end




      

!       subroutine test_exterior_pt(npatches,norder,npts,srcvals,&
!      &srccoefs,wts,xyzout,isout)
! !
! !
! !  this subroutine tests whether the pt xyzin, is
! !  in the exterior of a surface, and also estimates the error
! !  in representing e^{ir/2}/r and \grad e^{ir/2}/r \cdot n
! !  centered at the interior point. Whether a point 
! !  is in the interior or not is tested using Gauss' 
! !  identity for the flux due to a point charge
! !
! !
! !  input:
! !    npatches - integer
! !       number of patches
! !    norder - integer
! !       order of discretization
! !    npts - integer
! !       total number of discretization points on the surface
! !    srccoefs - double precision (9,npts)
! !       koornwinder expansion coefficients of geometry info
! !    xyzout -  double precision (3)
! !       point to be tested
! !
! !  output: 
! !    isout - boolean
! !      whether the target is in the interior or not
! !

!       implicit none
!       integer npatches,norder,npts,npols
!       double precision srccoefs(9,npts),srcvals(12,npts),xyzout(3),wts(npts)
!       double precision tmp(3)
!       double precision dpars,done,pi
!       double precision, allocatable :: rsurf(:),err_p(:,:) 
!       integer ipars,norderhead,nd
!       complex *16, allocatable :: sigma_coefs(:,:), sigma_vals(:,:)
!       complex *16 zk,val

!       integer ipatch,j,i
!       double precision ra,ds
!       logical isout

!       done = 1
!       pi = atan(done)*4

!       npols = (norder+1)*(norder+2)/2


!       zk = 0

!       ra = 0



!       do ipatch=1,npatches
!         do j=1,npols
!           i = (ipatch-1)*npols + j
!           call h3d_sprime(xyzout,srcvals(1,i),dpars,zk,ipars,val)
!           call cross_prod3d(srcvals(4,i),srcvals(7,i),tmp)
!           ds = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)
!           ra = ra + real(val)*wts(i)
!         enddo
!       enddo

!       if(abs(ra+4*pi).le.1.0d-3) isout = .false.
!       if(abs(ra).le.1.0d-3) isout = .true.

!       return
!       end

   




