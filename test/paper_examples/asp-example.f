      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:),targs(:,:)
      real *8, allocatable :: wts(:)
      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)
      real *8 errs(6),ts(2)
      real *8, allocatable :: rfacs(:,:)
      character *100 fname
      integer ipars(2)
      integer, allocatable :: row_ptr(:),col_ind(:)
      integer, allocatable :: iquad(:)
      real *8, allocatable :: srcover(:,:),wover(:)
      complex *16, allocatable :: sigmaover(:),slp_near(:)
      complex *16, allocatable :: pot(:)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)
      integer, allocatable :: ixyzso(:),nfars(:)

      real *8, allocatable :: pdis(:)

      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      real *8 xyz_out(3)
      complex *16, allocatable :: sigma(:)
      complex * 16 zpars(3)

      real *8 epss(4)
      integer norder_list(5),iasps(2,5)

      call prini(6,13)

      allocate(rfacs(2,6))

      norder = 7
      eps = 0.50001d-6

 1110 format(2x,i2,2x,e11.5,3(2x,i6),2x,e11.5,2x,i9,7(2x,e11.5))
 1120 format(2x,i2,2(2x,e11.5))

      igeomtype = 2
      iasps(1,1) = 20
      iasps(2,1) = 60

      iasps(1,2) = 24
      iasps(2,2) = 48
      
      iasps(1,3) = 35
      iasps(2,3) = 35

      iasps(1,4) = 48
      iasps(2,4) = 24

      iasps(1,5) = 60
      iasps(2,5) = 20
      

      ifaddsub = 1
      
      igeomtype = 2
      iquadtype = 1

      zk = 1.0d0
      zpars(1) = zk
      zpars(2) = 1.0d0
      zpars(3) = 0.0d0

      xyz_out(1) = -3.5d0
      xyz_out(2) = 3.1d0
      xyz_out(3) = 20.1d0

      do iasp = 5,5
        ipars(1) = iasps(1,iasp)
        ipars(2) = iasps(2,iasp)
        npatches = 2*ipars(1)*ipars(2)

        npols = (norder+1)*(norder+2)/2

        npts = npatches*npols
        allocate(srcvals(12,npts),srccoefs(9,npts))
        allocate(targs(3,npts))
        ifplot = 0


        call setup_geom(igeomtype,norder,npatches,ipars, 
     1     srcvals,srccoefs,ifplot,fname)


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
        allocate(pot(npts))

        call get_centroid_rads(npatches,norders,ixyzs,iptype,npts, 
     1     srccoefs,cms,rads)


        allocate(sigma(npts))

        iquadtype = 1

        do i=1,npts
          call h3d_slp(xyz_out,srcvals(1,i),dpars,zpars,ipars,sigma(i))
        enddo


     
        do i=1,npts
          targs(1,i) = srcvals(1,i)
          targs(2,i) = srcvals(2,i)
          targs(3,i) = srcvals(3,i)
        enddo

        allocate(ipatch_id(npts),uvs_targ(2,npts),pdis(npatches))
        do i=1,npts
          ipatch_id(i) = -1
          uvs_targ(1,i) = 0
          uvs_targ(2,i) = 0
        enddo

        call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts, 
     1           ipatch_id,uvs_targ)
        
        call get_patch_distortion(npatches,norders,ixyzs,iptype,npts,
     1    srccoefs,srcvals,wts,pdis)
        
        
        pavg = 0
        pmax = 0
        do i=1,npatches
          pavg = pavg + pdis(i)
          if(pdis(i).gt.pmax) pmax = pdis(i)
        enddo
        pavg = pavg/(npatches+0.0d0)

        open(unit=34,file='stell20_asp_order8.txt',
     1       access='append')
        write(34,1120) norder,pavg,pmax
        close(34)

        call get_rfacs(norder,iptype,rfac,rfac0)

        do i=1,npatches 
          rad_near(i) = rads(i)*rfac
        enddo
      

        call findnearmem(cms,npatches,rad_near,targs,npts,nnz)

        allocate(row_ptr(npts+1),col_ind(nnz))
      
        call findnear(cms,npatches,rad_near,targs,npts,row_ptr, 
     1        col_ind)

        allocate(iquad(nnz+1)) 
        call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind,
     1         iquad)

        nquad = iquad(nnz+1)-1
        allocate(slp_near(nquad))


        ndtarg = 3

        write(*,*) " "
        write(*,*) " "
        write(*,*) " "
        write(*,*) " "
        write(*,*) " "
        write(*,*) " "
        write(*,*) " "
        write(*,*) norder,rfac,rfac0,eps
        write(*,*) " "
        write(*,*) " "
        write(*,*) " "
        write(*,*) " "

        ikerorder = -1

        call cpu_time(t1)
        call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,
     1    rads,npts,srccoefs,ndtarg,npts,targs,ikerorder,zk,
     2    nnz,row_ptr,col_ind,rfac,nfars,ixyzso)
        call cpu_time(t2)
        tfar = t2-t1

        npts_over = ixyzso(npatches+1)-1

        print *, npts_over

        allocate(srcover(12,npts_over),sigmaover(npts_over),
     1     wover(npts_over))

          
        call oversample_geom(npatches,norders,ixyzs,iptype,npts, 
     1    srccoefs,nfars,ixyzso,npts_over,srcover)
        call get_qwts(npatches,nfars,ixyzso,iptype,npts_over,
     1    srcover,wover)


        nd = 2
        call oversample_fun_surf(nd,npatches,norders,ixyzs,iptype, 
     1     npts,sigma,nfars,ixyzso,npts_over,sigmaover)


        do i=1,nquad
          slp_near(i) = 0
        enddo


        call cpu_time(t1)
C$        t1 = omp_get_wtime()          

        call getnearquad_helm_comb_dir(npatches,norders,
     1    ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     1    ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,
     1    iquad,rfac0,nquad,slp_near)

        call cpu_time(t2)
C$        t2 = omp_get_wtime()          
        tquadgen = t2-t1

        call cpu_time(t1)
C$        t1 = omp_get_wtime()          
        call lpcomp_helm_comb_dir_addsub(npatches,norders,ixyzs,
     1   iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     2   eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,slp_near,
     3   sigma,nfars,npts_over,ixyzso,srcover,wover,pot)

        call cpu_time(t2)
C$        t2 = omp_get_wtime()          
        tlpcomp = t2-t1

        call cpu_time(t1)
C$        t1 = omp_get_wtime()          
        call hfmm3d_s_c_p(eps,zk,npts,targs,sigma,pot)

        call cpu_time(t2)
C$        t2 = omp_get_wtime()          
        tlpcompref = t2-t1

        deallocate(srcover,sigmaover,wover)

        rmem = (nquad+0.0d0)/(npts+0.0d0)
        rquadgen = (npts+0.0d0)/tquadgen

        rlpcomp = (npts+0.0d0)/tlpcomp
        rlpcompref = (npts+0.0d0)/tlpcompref

        rfar = (npts+0.0d0)/tfar

        rover = (npts_over+0.0d0)/(npts+0.0d0)



        open(unit=33,file='stell20_asp_mac_order8_iprec2.txt',
     1     access='append')
        write(33,1110) norder,pavg,npatches,npts,
     1       npts_over,rover,nquad,rmem,tquadgen,rquadgen,
     2       tlpcomp,rlpcomp,tlpcompref,rlpcompref
        close(33)
        deallocate(row_ptr,col_ind,iquad,slp_near)
        deallocate(srcvals,srccoefs,targs,wts,pdis)
        deallocate(norders,ixyzs,iptype,uvs_targ,ipatch_id)
        deallocate(ixyzso,nfars,cms,rads,rad_near,sigma,pot)
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

