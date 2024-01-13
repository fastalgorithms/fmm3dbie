      implicit real *8 (a-h,o-z) 
      implicit integer *8 (i-n)
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:),targs(:,:)
      real *8, allocatable :: wts(:)
      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)
      real *8 errs(6),ts(2)
      real *8, allocatable :: rfacs(:,:)
      character *100 fname
      integer *8 ipars(2)
      integer *8, allocatable :: row_ptr(:),col_ind(:)
      integer *8, allocatable :: iquad(:)
      real *8, allocatable :: srcover(:,:),wover(:)
      complex *16, allocatable :: uval(:),dudnval(:)
      complex *16, allocatable :: sigmaover(:),slp_near(:),dlp_near(:)
      complex *16, allocatable :: pot(:),potslp(:),potdlp(:)
      complex *16, allocatable :: potslp2(:)

      complex *16 zk

      integer *8, allocatable :: norders(:),ixyzs(:),iptype(:)
      integer *8, allocatable :: ixyzso(:),nfars(:)

      integer *8, allocatable :: ipatch_id(:),inode_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      real *8 xyz_out(3),xyz_in(3)
      complex *16, allocatable :: sigma(:)
      complex * 16 zpars(3)
      integer *8 int8_0,int8_1,int8_3,int8_12


      external h3d_ggq_comb,h3d_ggq_slp

      int8_0 = 0
      int8_1 = 1
      int8_3 = 3
      int8_12 = 12

      call prini(6,13)

      done = 1
      pi = atan(done)*4


      igeomtype = 2
      iasp = 4
      iref = 0
      iprec = 0


      if(iasp.eq.1) then
        ipars(1) = 5
        ipars(2) = 15
      endif

      if(iasp.eq.2) then
        ipars(1) = 6
        ipars(2) = 12
      endif

      if(iasp.eq.3) then
        ipars(1) = 9
        ipars(2) = 9
      endif

      if(iasp.eq.4) then
        ipars(1) = 12
        ipars(2) = 6
      endif

      if(iasp.eq.5) then
        ipars(1) = 15
        ipars(2) = 5
      endif

      ipars(1) = ipars(1)*2**(iref)
      ipars(2) = ipars(2)*2**(iref)

      npatches = 2*ipars(1)*ipars(2)

      norder = 4 
      npols = (norder+1)*(norder+2)/2
      write(fname,'(a,i1,a,i1,a,i1,a)') 'data/stell_aquad_iasp_',iasp,
     1  '_iref_',iref,'_norder_',norder,'.dat'


      zk = 1.0d0
      zpars(1) = zk
      zpars(2) = 1.0d0
      zpars(3) = 0.0d0

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
        call h3d_slp(xyz_out,int8_3,srcvals(1,i),int8_0,dpars,int8_1,
     1     zpars,int8_0,ipars,uval(i))
        call h3d_sprime(xyz_out,int8_12,srcvals(1,i),int8_0,dpars,
     1     int8_1,zpars,int8_0,ipars,dudnval(i))
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
      rfac = 2.75d0
      rfac0 = 1.25d0
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

      eps = 0.50001d-12

      ikerorder = -1


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

      zpars(1) = zk
      zpars(2) = 1.0d0
      zpars(3) = 0.0d0

      iquadtype = 1

cc      goto 1111

      call getnearquad_helm_comb_dir(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     1      ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,
     1      iquad,
     1      rfac0,nquad,slp_near)

      
      zpars(2) = 0.0d0
      zpars(3) = 1.0d0
      call getnearquad_helm_comb_dir(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     1      ipatch_id,uvs_targ,eps,zpars,iquadtype,
     1      nnz,row_ptr,col_ind,iquad,
     1      rfac0,nquad,dlp_near)
      
      call cpu_time(t2)
      tquadgen = t2-t1



      ifinout = 1     

      zpars(2) = 1.0d0
      zpars(3) = 0.0d0


      call cpu_time(t1)

      call lpcomp_helm_comb_dir_addsub(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     2  eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,slp_near,
     3  dudnval,nfars,npts_over,ixyzso,srcover,wover,potslp)


      zpars(2) = 0.0d0
      zpars(3) = 1.0d0


      call lpcomp_helm_comb_dir_setsub(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     2  eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,dlp_near,
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


      open(unit=23,file=fname,action='readwrite',form='unformatted',
     1  access='stream')
      write(unit=23) slp_near
      write(unit=23) dlp_near
      
      close(23)

      stop
      end




      subroutine setup_geom(igeomtype,norder,npatches,ipars, 
     1    srcvals,srccoefs,ifplot,fname)
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8 igeomtype,norder,npatches,ipars(*),ifplot
      character (len=*) fname
      real *8 srcvals(12,*), srccoefs(9,*)
      real *8, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts(:)

      real *8, pointer :: ptr1,ptr2,ptr3,ptr4
      integer *8, pointer :: iptr1,iptr2,iptr3,iptr4
      real *8, target :: p1(10),p2(10),p3(10),p4(10)
      real *8, allocatable, target :: triaskel(:,:,:)
      real *8, allocatable, target :: deltas(:,:)
      integer *8, allocatable :: isides(:)
      integer *8, target :: nmax,mmax

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

