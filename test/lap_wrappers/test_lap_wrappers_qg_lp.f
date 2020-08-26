      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:),targs(:,:)
      real *8, allocatable :: wts(:)
      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)
      real *8 errs(6),ts(2)
      real *8, allocatable :: rfacs(:,:)
      complex *16 zk
      character *100 fname
      integer ipars(2)
      integer, allocatable :: row_ptr(:),col_ind(:)
      integer, allocatable :: iquad(:)
      real *8, allocatable :: srcover(:,:),wover(:)
      real *8, allocatable :: uval(:),dudnval(:)
      real *8, allocatable :: sigmaover(:),slp_near(:),dlp_near(:)
      real *8, allocatable :: pot(:),potslp(:),potdlp(:)
      real *8, allocatable :: potslp2(:)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)
      integer, allocatable :: ixyzso(:),nfars(:)

      integer, allocatable :: ipatch_id(:),inode_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      real *8 xyz_out(3),xyz_in(3)
      real *8, allocatable :: sigma(:)
      real *8 dpars(2)


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


      dpars(1) = 1.0d0
      dpars(2) = 0.0d0

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



      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))
      allocate(pot(npts),potslp(npts),potdlp(npts))

      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts, 
     1     srccoefs,cms,rads)

      allocate(sigma(npts),uval(npts),dudnval(npts))

      do i=1,npts
        call l3d_slp(xyz_out,3,srcvals(1,i),0,dpars,0,zpars,0,
     1     ipars,uval(i))
        call l3d_sprime(xyz_out,12,srcvals(1,i),0,dpars,0,zpars,0,
     1     ipars,dudnval(i))
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

 
c
c    find near field
c
      iptype = 1
      call get_rfacs(norder,iptype,rfac,rfac0)
      do i=1,npatches 
        rad_near(i) = rads(i)*rfac
      enddo
      

      call findnearmem(cms,npatches,rad_near,ndtarg,targs,npts,nnz)

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,ndtarg,targs,npts,row_ptr, 
     1        col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind,
     1         iquad)

      nquad = iquad(nnz+1)-1
      allocate(slp_near(nquad),dlp_near(nquad))


      ndtarg = 3

      eps = 0.50001d-6

      ikerorder = -1
      zk = 0


      call cpu_time(t1)
      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,
     1    rads,npts,srccoefs,ndtarg,npts,targs,ikerorder,zk,
     2    nnz,row_ptr,col_ind,rfac,nfars,ixyzso)
      call cpu_time(t2)
      tfar = t2-t1


      npts_over = ixyzso(npatches+1)-1

      print *, "npts_over=",npts_over


      allocate(srcover(12,npts_over),sigmaover(npts_over),
     1         wover(npts_over))

          
      call oversample_geom(npatches,norders,ixyzs,iptype,npts, 
     1   srccoefs,srcvals,nfars,ixyzso,npts_over,srcover)

      call get_qwts(npatches,nfars,ixyzso,iptype,npts_over,
     1        srcover,wover)


      do i=1,nquad
        slp_near(i) = 0
        dlp_near(i) = 0
      enddo



      call cpu_time(t1)

      dpars(1) = 1.0d0
      dpars(2) = 0.0d0

      iquadtype = 1

cc      goto 1111

      call getnearquad_lap_comb_dir(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     1      ipatch_id,uvs_targ,eps,dpars,iquadtype,nnz,row_ptr,col_ind,
     1      iquad,rfac0,nquad,slp_near)

      
      dpars(1) = 0.0d0
      dpars(2) = 1.0d0
      call getnearquad_lap_comb_dir(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     1      ipatch_id,uvs_targ,eps,dpars,iquadtype,
     1      nnz,row_ptr,col_ind,iquad,
     1      rfac0,nquad,dlp_near)
      
      call cpu_time(t2)
      tquadgen = t2-t1



      ifinout = 1     

      dpars(1) = 1.0d0
      dpars(2) = 0.0d0


      call cpu_time(t1)

      call lpcomp_lap_comb_dir_addsub(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     2  eps,dpars,nnz,row_ptr,col_ind,iquad,nquad,slp_near,
     3  dudnval,nfars,npts_over,ixyzso,srcover,wover,potslp)


      dpars(1) = 0.0d0
      dpars(2) = 1.0d0


      call lpcomp_lap_comb_dir_setsub(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     2  eps,dpars,nnz,row_ptr,col_ind,iquad,nquad,dlp_near,
     3  uval,nfars,npts_over,ixyzso,srcover,wover,potdlp)


      call cpu_time(t2)
      tlpcomp = t2-t1


c
c
c      compute error
c
      errl2 = 0
      rl2 = 0
      do i=1,npts
        pot(i) = (potslp(i) - potdlp(i))*2
        errl2 = errl2 + abs(uval(i)-pot(i))**2*wts(i)
        rl2 = rl2 + abs(uval(i))**2*wts(i)
      enddo


      err = sqrt(errl2/rl2)

      call prin2('error in greens identity=*',err,1)

      i1 = 0
      if(err.lt.1.0d-2) i1 = 1

      allocate(potslp2(npts))

      dpars(1) = 1.0d0
      dpars(2) = 0.0d0
      
      call lpcomp_lap_comb_dir(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,ipatch_id,
     2  uvs_targ,eps,dpars,dudnval,potslp2)


      errl2 = 0
      rl2 = 0
      do i=1,npts

        errl2 = errl2 + abs(potslp(i)-potslp2(i))**2*wts(i)
        rl2 = rl2 + abs(potslp(i))**2*wts(i) 
      enddo
      errl2 = sqrt(errl2/rl2)

      call prin2('error in simpler calling interface for lp eval=*',
     1   errl2,1)

      i2 = 0
      if(errl2.lt.1.0d-12) i2 = 1
      ntests = 2

      nsuccess = i1+i2

      open(unit=33,file='../../print_testres.txt',access='append')
      write(33,'(a,i1,a,i1,a)') 'Successfully completed ',nsuccess,
     1  ' out of ',ntests,' in lap_wrappers testing suite'
      close(33)
      
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

