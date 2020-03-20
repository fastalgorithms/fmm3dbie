      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:),targs(:,:)
      real *8, allocatable :: wts(:)
      real *8 errs(6),ts(2)
      character *100 fname
      integer ipars(2)
      complex *16, allocatable :: uval(:),dudnval(:)
      complex *16, allocatable :: pot(:),potslp(:),potdlp(:)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)

      integer, allocatable :: ipatch_id(:),inode_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      real *8 xyz_out(3)
      complex * 16 zpars(3)
      integer norder_list(5)


      external h3d_ggq_comb,h3d_ggq_slp

      call prini(6,13)

      done = 1
      pi = atan(done)*4


      norder_list(1) = 2
      norder_list(2) = 3
      norder_list(3) = 4
      norder_list(4) = 6
      norder_list(5) = 8

c
c       select geometry type
c       igeomtype = 1 => sphere
c       igeomtype = 2 => stellarator
c 
      igeomtype = 2

      zk = 1.0d0
      eps = 0.50001d-6

      if(igeomtype.eq.1) then
        xyz_out(1) = 3.17d0
        xyz_out(2) = -0.03d0
        xyz_out(3) = 3.15d0
      endif

      if(igeomtype.eq.2) then
        xyz_out(1) = -3.5d0
        xyz_out(2) = 3.1d0
        xyz_out(3) = 20.1d0
      endif

      do iref = 1,3

        if(igeomtype.eq.1) ipars(1) = 2+iref
        if(igeomtype.eq.2) ipars(1) = 10*2**(iref-1)

        if(igeomtype.eq.1) then
          npatches = 12*(4**ipars(1))
        endif
        if(igeomtype.eq.2) then
          ipars(2) = ipars(1)*3
          npatches = 2*ipars(1)*ipars(2)
        endif


        zpars(1) = zk
        zpars(2) = 1.0d0
        zpars(3) = 0.0d0

        allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))

        do iorder_list = 1,1

          norder = norder_list(iorder_list)-1 
          print *, norder,npatches
          npols = (norder+1)*(norder+2)/2

          npts = npatches*npols
          allocate(srcvals(12,npts),srccoefs(9,npts))
          allocate(targs(3,npts))
          ifplot = 0



          call setup_geom(igeomtype,norder,npatches,ipars, 
     1       srcvals,srccoefs,ifplot,fname)


          do i=1,npatches
            norders(i) = norder
            ixyzs(i) = 1 +(i-1)*npols
            iptype(i) = 1
          enddo

          ixyzs(npatches+1) = 1+npols*npatches
          allocate(wts(npts))
          call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)


          allocate(pot(npts),potslp(npts),potdlp(npts))
          allocate(uval(npts),dudnval(npts))

          do i=1,npts
            call h3d_slp(xyz_out,srcvals(1,i),dpars,zpars,ipars,
     1         uval(i))
            call h3d_sprime(xyz_out,srcvals(1,i),dpars,zpars,ipars,
     1         dudnval(i))
          enddo

          ndtarg = 3
     
          do i=1,npts
            targs(1,i) = srcvals(1,i)
            targs(2,i) = srcvals(2,i)
            targs(3,i) = srcvals(3,i)
          enddo

          allocate(ipatch_id(npts),uvs_targ(2,npts))
          do i=1,npts
            ipatch_id(i) = -1
            uvs_targ(1,i) = 0
            uvs_targ(2,i) = 0
          enddo

          call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts, 
     1         ipatch_id,uvs_targ)


          call cpu_time(t1)
          zpars(1) = zk
          zpars(2) = 1.0d0
          zpars(3) = 0.0d0
          call lpcomp_helm_comb_dir(npatches,norders,ixyzs,
     1      iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,ipatch_id,
     2      uvs_targ,eps,zpars,dudnval,potslp)
      

          zpars(2) = 0.0d0
          zpars(3) = 1.0d0
          call lpcomp_helm_comb_dir(npatches,norders,ixyzs,
     1      iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,ipatch_id,
     2      uvs_targ,eps,zpars,uval,potdlp)
      
c
c
c      compute error
c
          errl2 = 0
          rl2 = 0
          do i=1,npts
            pot(i) = (potslp(i) - potdlp(i))/2/pi
            errl2 = errl2 + abs(uval(i)-pot(i))**2*wts(i)
            rl2 = rl2 + abs(uval(i))**2*wts(i)
          enddo
          err = sqrt(errl2/rl2)

          call prin2('error in greens identity=*',err,1)
          open(unit=33,file='stell_ooc_green_mac.txt',
     1       access='append')
          write(33,*) npatches,norder,err
          deallocate(srcvals,srccoefs,targs,wts,uval,dudnval)
          deallocate(pot,potslp,potdlp,ipatch_id,uvs_targ)
        enddo
        deallocate(norders,ixyzs,iptype)
      enddo


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

