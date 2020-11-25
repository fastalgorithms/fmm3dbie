      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:),targs(:,:)
      real *8, allocatable :: wts(:)
      
      real *8 ts(2), rres
      real *8, allocatable :: rfacs(:,:), errs(:)
      character *100 fname
      integer ipars(2)
      integer, allocatable :: row_ptr(:),col_ind(:)
      integer, allocatable :: iquad(:)
      real *8, allocatable :: srcover(:,:),wover(:)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)
      integer, allocatable :: ixyzso(:),nfars(:)

      integer, allocatable :: ipatch_id(:),inode_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      real *8 xyz_out(3),xyz_in(3),stracmat(3,3),smat(3,3), dmat(3,3)
      real *8 xyz_src(3),xyz_targ(3)
      real *8 velgrad(3,3), vel(3), pre, tractemp(3)
      real *8 sigout(3), uin(3), uintest(3), dpars(2), st1(3), du1(3)
      real *8 udir(3), uneu(3,10), uavecomp(3), uavetest(3)
      real *8 st2(3), du2(3), uconst(3)
      real *8 v(3), omega(3), r0(3), udiff(3,10), udiff2(3,10)      
      real *8, allocatable :: uval(:,:), tracval(:,:), soln(:,:)
      complex * 16 zpars

      call prini(6,13)

      done = 1
      pi = atan(done)*4


c
c       select geometry type
c       igeomtype = 1 => sphere
c       igeomtype = 2 => stellarator
c 
      igeomtype = 2
      if(igeomtype.eq.1) ipars(1) = 1
      if(igeomtype.eq.2) ipars(1) = 5*2

      if(igeomtype.eq.1) then
        npatches = 12*(4**ipars(1))
      endif
      if(igeomtype.eq.2) then
        ipars(2) = ipars(1)*3
        npatches = 2*ipars(1)*ipars(2)
      endif


      if(igeomtype.eq.1) then
        xyz_out(1) = 3.17d0
        xyz_out(2) = -0.03d0
        xyz_out(3) = 3.15d0

        xyz_in(1) = 0.17d0
        xyz_in(2) = 0.23d0
        xyz_in(3) = -0.11d0
      endif

      if(igeomtype.eq.2) then
        xyz_in(1) = -4.501d0
        xyz_in(2) = 1.7d-3
        xyz_in(3) = 0.00001d0

        xyz_out(1) = -3.5d0
        xyz_out(2) = 3.1d0
        xyz_out(3) = 20.1d0
      endif
      ifinout = 1

      if(ifinout.eq.0) then
        xyz_src(1) = xyz_out(1) 
        xyz_src(2) = xyz_out(2) 
        xyz_src(3) = xyz_out(3)

        xyz_targ(1) = xyz_in(1)
        xyz_targ(2) = xyz_in(2)
        xyz_targ(3) = xyz_in(3)
      else
        xyz_src(1) = xyz_in(1) 
        xyz_src(2) = xyz_in(2) 
        xyz_src(3) = xyz_in(3)

        xyz_targ(1) = xyz_out(1)
        xyz_targ(2) = xyz_out(2)
        xyz_targ(3) = xyz_out(3)
      endif

      norder = 4
      npols = (norder+1)*(norder+2)/2

      npts = npatches*npols
      allocate(srcvals(12,npts),srccoefs(9,npts))
      allocate(targs(3,npts))
      ifplot = 0



      call setup_geom(igeomtype,norder,npatches,ipars, 
     1       srcvals,srccoefs,ifplot,fname)

      allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))
      allocate(ixyzso(npatches+1),nfars(npatches))

      do i=1,npatches
        norders(i) = norder
        ixyzs(i) = 1 +(i-1)*npols
        iptype(i) = 1
      enddo

      ixyzs(npatches+1) = 1+npols*npatches
      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)

      allocate(uval(3,npts),tracval(3,npts))

      sigout(1) = 1.1d0
      sigout(2) = -0.27d0
      sigout(3) = .31d0
      
      do i=1,npts
         call st3d_slp_vec(9,xyz_src,3,srcvals(1,i),0,dpars,0,zpars,0,
     1        ipars,smat)
         uval(1,i) = smat(1,1)*sigout(1) + smat(1,2)*sigout(2)
     1        + smat(1,3)*sigout(3)
         uval(2,i) = smat(2,1)*sigout(1) + smat(2,2)*sigout(2)
     1        + smat(2,3)*sigout(3)
         uval(3,i) = smat(3,1)*sigout(1) + smat(3,2)*sigout(2)
     1        + smat(3,3)*sigout(3)
      enddo



      call st3d_slp_vec(9,xyz_src,3,xyz_targ,0,dpars,0,zpars,0,
     1     ipars,smat)
      
       uintest(1) = smat(1,1)*sigout(1) + smat(1,2)*sigout(2)
     1        + smat(1,3)*sigout(3)
       uintest(2) = smat(2,1)*sigout(1) + smat(2,2)*sigout(2)
     1        + smat(2,3)*sigout(3)
       uintest(3) = smat(3,1)*sigout(1) + smat(3,2)*sigout(2)
     1        + smat(3,3)*sigout(3)


      npt1 = 1
      allocate(ipatch_id(npt1),uvs_targ(2,npt1))
      do i=1,npt1
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo
      
c
c     solve dirichlet (velocity) problem
c

      alpha = 1
      beta = 1
      
      dpars(1) = alpha
      dpars(2) = beta
      numit = 200

      allocate(soln(3,npts),errs(numit+1))

      eps_gmres = 1d-7
      eps = 1d-7
      
      call stok_comb_vel_solver(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,eps,dpars,numit,
     2     ifinout,uval,eps_gmres,niter,errs,rres,soln)


      call prin2('gmres errs *',errs,niter)
      call prin2('gmres rres *',rres,1)

      ndt_in = 3
      nt_in = 1
      call lpcomp_stok_comb_vel(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndt_in,nt_in,xyz_targ,
     2     ipatch_id,uvs_targ,eps,dpars,soln,udir)

      udir(1) = udir(1) 
      udir(2) = udir(2) 
      udir(3) = udir(3) 

      sum = 0
      sumrel = 0

      do i = 1,3
         sum = sum + (udir(i)-uintest(i))**2
         sumrel = sumrel + (uintest(i))**2
      enddo
      
      call prin2('rel err in velocity *',sqrt(sum/sumrel),1)

c
c  solve scattering problem
c
      do i=1,npts
        uval(1,i) = -1
        uval(2,i) = 0
        uval(3,i) = 0
        soln(1,i) = 0
        soln(2,i) = 0
        soln(3,i) = 0
      enddo
      
      call stok_comb_vel_solver(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,eps,dpars,numit,
     2     ifinout,uval,eps_gmres,niter,errs,rres,soln)


      call prin2('gmres errs *',errs,niter)
      call prin2('gmres rres *',rres,1)

      call surf_vtk_plot_vec(npatches,norders,ixyzs,iptype,
     1   npts,srccoefs,srcvals,soln,'stell-stok-scat-soln-ref1.vtk','a')

      
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

