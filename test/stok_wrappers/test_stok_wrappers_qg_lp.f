      implicit real *8 (a-h,o-z) 
      implicit integer(8) (i-n)
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:),targs(:,:)
      real *8, allocatable :: wts(:)
      
      real *8 ts(2), rres
      real *8, allocatable :: rfacs(:,:), errs(:)
      character *100 fname
      integer(8) ipars(2)
      integer(8), allocatable :: row_ptr(:),col_ind(:)
      integer(8), allocatable :: iquad(:)
      real *8, allocatable :: srcover(:,:),wover(:)

      integer(8), allocatable :: norders(:),ixyzs(:),iptype(:)
      integer(8), allocatable :: ixyzso(:),nfars(:)

      integer(8), allocatable :: ipatch_id(:),inode_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      real *8 xyz_out(3),xyz_in(3,10),stracmat(3,3),smat(3,3), dmat(3,3)
      real *8 velgrad(3,3), vel(3), pre, tractemp(3)
      real *8 sigout(3), uin(3), uintest(3,10), dpars(2), st1(3), du1(3)
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
      igeomtype = 1
      if(igeomtype.eq.1) ipars(1) = 1
      if(igeomtype.eq.2) ipars(1) = 20

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

        xyz_in(1,1) = 0.17d0
        xyz_in(2,1) = 0.23d0
        xyz_in(3,1) = -0.11d0

        do i = 2,10
           xyz_in(1,i) = 0.2d0*cos(100.0d0*i+5.0d0)
           xyz_in(2,i) = 0.2d0*cos(211.0d0*i+5.0d0)
           xyz_in(3,i) = 0.2d0*cos(357.0d0*i+5.0d0)
        enddo

      endif

      if(igeomtype.eq.2) then
        xyz_in(1,1) = -4.501d0
        xyz_in(2,1) = 1.7d-3
        xyz_in(3,1) = 0.00001d0

        xyz_out(1) = -3.5d0
        xyz_out(2) = 3.1d0
        xyz_out(3) = 20.1d0
      endif

      norder = 3
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
         call st3d_slp_vec(9,xyz_out,3,srcvals(1,i),0,dpars,0,zpars,0,
     1        ipars,smat)
         call st3d_strac_vec(9,xyz_out,12,srcvals(1,i),0,dpars,0,zpars,
     1        0,ipars,stracmat)
         uval(1,i) = smat(1,1)*sigout(1) + smat(1,2)*sigout(2)
     1        + smat(1,3)*sigout(3)
         uval(2,i) = smat(2,1)*sigout(1) + smat(2,2)*sigout(2)
     1        + smat(2,3)*sigout(3)
         uval(3,i) = smat(3,1)*sigout(1) + smat(3,2)*sigout(2)
     1        + smat(3,3)*sigout(3)
         tracval(1,i) = stracmat(1,1)*sigout(1)+ stracmat(1,2)*sigout(2)
     1        + stracmat(1,3)*sigout(3)
         tracval(2,i) = stracmat(2,1)*sigout(1)+ stracmat(2,2)*sigout(2)
     1        + stracmat(2,3)*sigout(3)
         tracval(3,i) = stracmat(3,1)*sigout(1)+ stracmat(3,2)*sigout(2)
     1        + stracmat(3,3)*sigout(3)
         
      enddo



      nin = 9
      
      do i = 1,nin
         call st3d_slp_vec(9,xyz_out,3,xyz_in(1,i),0,dpars,0,zpars,0,
     1        ipars,smat)
      
         uintest(1,i) = smat(1,1)*sigout(1) + smat(1,2)*sigout(2)
     1        + smat(1,3)*sigout(3)
         uintest(2,i) = smat(2,1)*sigout(1) + smat(2,2)*sigout(2)
     1        + smat(2,3)*sigout(3)
         uintest(3,i) = smat(3,1)*sigout(1) + smat(3,2)*sigout(2)
     1        + smat(3,3)*sigout(3)

      enddo


c
c     basic kernel check: green's i.d. with smooth quad at 1 target
c     

      uin(1) = 0
      uin(2) = 0
      uin(3) = 0

      st2(1) = 0
      st2(2) = 0
      st2(3) = 0

      du2(1) = 0
      du2(2) = 0
      du2(3) = 0

      do i = 1,npts
         call st3d_slp_vec(9,srcvals(1,i),3,xyz_in,0,dpars,0,zpars,0,
     1        ipars,smat)
         call st3d_dlp_vec(9,srcvals(1,i),3,xyz_in,0,dpars,0,zpars,
     1        0,ipars,dmat)

         u1 = uval(1,i)
         u2 = uval(2,i)
         u3 = uval(3,i)
         t1 = tracval(1,i)
         t2 = tracval(2,i)
         t3 = tracval(3,i)         

         uin(1) = uin(1) + wts(i)*(
     1        smat(1,1)*t1+smat(1,2)*t2+smat(1,3)*t3
     2        -dmat(1,1)*u1-dmat(1,2)*u2-dmat(1,3)*u3)
         uin(2) = uin(2) + wts(i)*(
     1        smat(2,1)*t1+smat(2,2)*t2+smat(2,3)*t3
     2        -dmat(2,1)*u1-dmat(2,2)*u2-dmat(2,3)*u3)
         uin(3) = uin(3) + wts(i)*(
     1        smat(3,1)*t1+smat(3,2)*t2+smat(3,3)*t3
     2        -dmat(3,1)*u1-dmat(3,2)*u2-dmat(3,3)*u3)

         st2(1) = st2(1) + wts(i)*(
     1        smat(1,1)*t1+smat(1,2)*t2+smat(1,3)*t3)
         st2(2) = st2(2) + wts(i)*(
     1        smat(2,1)*t1+smat(2,2)*t2+smat(2,3)*t3)
         st2(3) = st2(3) + wts(i)*(
     1        smat(3,1)*t1+smat(3,2)*t2+smat(3,3)*t3)

         du2(1) = du2(1) + wts(i)*(
     2        dmat(1,1)*u1+dmat(1,2)*u2+dmat(1,3)*u3)
         du2(2) = du2(2) + wts(i)*(
     2        dmat(2,1)*u1+dmat(2,2)*u2+dmat(2,3)*u3)
         du2(3) = du2(3) + wts(i)*(
     2        dmat(3,1)*u1+dmat(3,2)*u2+dmat(3,3)*u3)


      enddo

      uin(1) = uin(1)/(4*pi)
      uin(2) = uin(2)/(4*pi)
      uin(3) = uin(3)/(4*pi)      

      sum = 0
      sumrel = 0

      do i = 1,3
         sum = sum + (uin(i)-uintest(i,1))**2
         sumrel = sumrel + uintest(i,1)**2
      enddo

      i1 = 0
      errl2 = sqrt(sum/sumrel)
      if (errl2 .lt. 1d-3) i1 = 1
      
      call prin2('rel err in Gauss ID --- direct, smooth quad *',
     1     errl2,1)

c
c     same test but calling fmm interface
c      


      npt1 = nin
      allocate(ipatch_id(npt1),uvs_targ(2,npt1))
      do i=1,npt1
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo


      ndt_in = 3
      nt_in = 1
      eps = 1d-12
      dpars(1) = 1
      dpars(2) = 0
      st1(1) = 0
      st1(2) = 0
      st1(3) = 0      
      call lpcomp_stok_comb_vel(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndt_in,nt_in,xyz_in,
     2     ipatch_id,uvs_targ,eps,dpars,tracval,st1)

      dpars(1) = 0
      dpars(2) = 1
      du1(1) = 0
      du1(2) = 0
      du1(3) = 0      
      call lpcomp_stok_comb_vel(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndt_in,nt_in,xyz_in,
     2     ipatch_id,uvs_targ,eps,dpars,uval,du1)



      uin(1) = (st1(1)-du1(1))/(4*pi)
      uin(2) = (st1(2)-du1(2))/(4*pi)
      uin(3) = (st1(3)-du1(3))/(4*pi)

      sum = 0
      sumrel = 0
      do i = 1,3
         sum = sum + (uin(i)-uintest(i,1))**2
         sumrel = sumrel + uintest(i,1)**2
      enddo

      i2= 0
      errl2 = sqrt(sum/sumrel)
      if (errl2 .lt. 1d-3) i2 = 2 

      call prin2('rel err in Gauss ID --- FMM, adaptive quadrature *',
     1     errl2,1)

      ntests=2
      nsuccess = i1+i2

      open(unit=33,file='../../print_testres.txt',access='append')
      write(33,'(a,i1,a,i1,a)') 'Successfully completed ',nsuccess,
     1     ' out of ',ntests,' in stok_wrappers testing suite'
      close(33)
      
      
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

