program efie_rwg_open

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
  integer :: ipatch_id
  integer :: numit,niter

  character *200 title,fname,fname1,fname2,fname3

  real *8 uvs_targ(2)

  logical isout0,isout1

  double complex :: pot,potex,ztmp,ima,zk
  complex *16 alpha_rhs

  integer nnn,mmm
  real ( kind = 8 ) xyz_0(3),r0

  integer count1

  integer norder, n, npatches_aux,npts_aux,npoints_est
  integer npatches,npts,nrwg,n_vert_flat,n_tri_flat
  real ( kind = 8 ), allocatable :: vert_flat(:,:)
  integer, allocatable :: tri_flat(:,:),rwg(:,:)

  data ima/(0.0d0,1.0d0)/


  !
  ! the driver code begins here... -------------------------------------------
  !


  call prini(6,13)

  done = 1
  pi = atan(done)*4

  zk = 2.0d-1
  zk = 1.0d0
  !      zk=4.4934d0 this frequency is interesting to test the resonance
  ! of a mfie on a sphere of radius 1
  !https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=4334247

  zpars(1) = zk
  zpars(2) = 0.0d0
  zpars(3) = -1.0d0

  !zpars(2) and zpars(3) must be real and with different sign to avoid resonances

  xyz_in(1) = 0.11d0
  xyz_in(2) = 0.0d-5
  xyz_in(3) = 0.37d0

  xyz_out(1) = -3.5d0
  xyz_out(2) = 3.1d0
  xyz_out(3) = 20.1d0
      
  fname = '../../geometries/sphere_768_o03.go3'
      
  eps = 1d-4
  eps_gmres=1d-5
  n = 2

  write (*,*) 'fname: ',fname


  call open_gov3_geometry_mem(fname, npatches_aux, npts_aux)
  call prinf('npatches=*',npatches_aux,1)
  call prinf('npts=*',npts_aux,1)

  norder=3
  npatches=(n**2)*npatches_aux
  npoints_est=npatches*3
  allocate(rwg(4,3*npatches))

  npts=npatches*(norder+1)*(norder+2)/2
  allocate(srcvals(12,npts),srccoefs(9,npts))
  allocate(ixyzs(npatches+1),iptype(npatches),norders(npatches))
  allocate(wts(npts))
     
  allocate(vert_flat(3,npoints_est))
  allocate(tri_flat(3,npatches))

  call open_go3_4_rwg(fname, npts_aux, npatches_aux, norder, npatches, &
       norders, ixyzs, &
       iptype, npts, srcvals, srccoefs, wts, n, npoints_est, nrwg, rwg, &
       vert_flat,&
       tri_flat, n_vert_flat, n_tri_flat)

  norder = norders(1)
  allocate(sigma(npts),rhs(npts))
  allocate(sigma2(npts),rhs2(npts))
  allocate(a_vect(nrwg),RHS_vect(nrwg),rhs_nE(npts), &
       rhs_nE_aux(npts),a_s(npts))

  ifinout = 1
  
  thet = hkrand(0)*pi
  phi = hkrand(0)*2*pi
  
  vf(1)=1.0d0
  vf(2)=2.0d0
  vf(3)=3.0d0
  write (*,*) 'about to get rhs'


  call get_rhs_em_cfie_rwg_pec_sphere(xyz_out,vf,npts,srcvals,zpars,RHS_vect,&
       rwg,nrwg,vert_flat,tri_flat,n_vert_flat,n_tri_flat,norder,wts)



  numit = 400
  niter = 0
  allocate(errs(numit+1))

  call cpu_time(t1)
  !$ t1 = omp_get_wtime()

  call em_cfie_rwg_solver(npatches,norders,ixyzs,iptype,npts,&
       srccoefs,srcvals,eps,zpars,numit,ifinout,RHS_vect,eps_gmres,&
       niter,errs,rres,a_vect,rhs_nE,&
       rwg,nrwg,vert_flat,tri_flat,n_vert_flat,n_tri_flat,norder,wts)
       
  call cpu_time(t2)
  !$ t2 = omp_get_wtime()
  call prin2('time taken in efie solver=*',t2-t1,1)
 
  call test_accuracy_em_cfie_rwg_pec_sphere(eps,a_vect,zpars,npts,wts, &
       srcvals,xyz_in,vf,xyz_out,&
       rwg,nrwg,vert_flat,tri_flat,n_vert_flat,n_tri_flat,norder)

  stop
end program efie_rwg_open











subroutine setup_geom(igeomtype,norder,npatches,ipars,&
     srcvals,srccoefs,ifplot,fname)
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

  if(igeomtype .eq. 1) then
     itype = 2
     allocate(triaskel(3,3,npatches))
     allocate(isides(npatches))
     npmax = npatches
     ntri = 0
     call xtri_platonic(itype, ipars(1), npmax, ntri,&
          triaskel, isides)

     xtri_geometry => xtri_sphere_eval
     ptr1 => triaskel(1,1,1)
     ptr2 => p2(1)
     ptr3 => p3(1)
     ptr4 => p4(1)

     if(ifplot.eq.1) then
        call xtri_vtk_surf(fname, npatches,xtri_geometry, ptr1,ptr2, &
             ptr3,ptr4, norder,'Triangulated surface of the sphere')
     endif


     call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4,&
          npols,uvs,umatr,srcvals,srccoefs)
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
end subroutine setup_geom






subroutine test_exterior_pt(npatches,norder,npts,srcvals,&
     srccoefs,wts,xyzout,isout)
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

   




