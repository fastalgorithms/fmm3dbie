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
      complex *16, allocatable :: uval(:),dudnval(:)
      complex *16, allocatable :: sigmaover(:),slp_near(:)
      complex *16, allocatable :: pot(:)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)
      integer, allocatable :: ixyzso(:),nfars(:)

      integer, allocatable :: ipatch_id(:),inode_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      real *8 xyz_out(3)
      complex *16, allocatable :: sigma(:)
      complex * 16 zpars(3)

      real *8 epss(4)
      integer norder_list(5)

      external h3d_ggq_comb,h3d_ggq_slp

      call prini(6,13)


c
c       select geometry type
c       igeomtype = 1 => sphere
c       igeomtype = 2 => stellarator
c 
      igeomtype = 2
      ipars(1) = 20

      if(igeomtype.eq.1) then
        npatches = 12*(4**ipars(1))
      endif
      if(igeomtype.eq.2) then
        ipars(2) = ipars(1)*3
        npatches = 2*ipars(1)*ipars(2)
      endif


      zk = 1.0d0
      zpars(1) = zk
      zpars(2) = 1.0d0
      zpars(3) = 0.0d0

      call get_exterior_pt(igeomtype,xyz_out)

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


      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))
      allocate(pot(npts))

      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts, 
     1     srccoefs,cms,rads)


      allocate(sigma(npts))

      do i=1,npts
        call h3d_slp(xyz_out,srcvals(1,i),dpars,zpars,ipars,sigma(i))
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
        inode_id(i) = i
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

      eps = 0.50001d-6

      ikerorder = -1

      call cpu_time(t1)
      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,
     1    rads,npts,srccoefs,ndtarg,npts,targs,ikerorder,zk,
     2    nnz,row_ptr,col_ind,rfac,nfars,ixyzso)
      call cpu_time(t2)
      tfar = t2-t1

cc          call prinf('ixyzso=*',ixyzso,npatches+1)


      npts_over = ixyzso(npatches+1)-1


      allocate(srcover(12,npts_over),sigmaover(npts_over),
     1         wover(npts_over))

          
      call oversample_geom(npatches,norders,ixyzs,iptype,npts, 
     1   srccoefs,nfars,ixyzso,npts_over,srcover)

      call get_qwts(npatches,nfars,ixyzso,iptype,npts_over,
     1        srcover,wover)


      nd = 2
      call oversample_fun_surf(nd,npatches,norders,ixyzs,iptype, 
     1      npts,sigma,nfars,ixyzso,npts_over,sigmaover)


      do i=1,nquad
        slp_near(i) = 0
      enddo


      call cpu_time(t1)

      zpars(1) = zk
      zpars(2) = 1.0d0
      zpars(3) = 0.0d0

      call getnearquad_helm_comb_dir_bdry(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     1      ipatch_id,uvs_targ,eps,zpars,nnz,row_ptr,col_ind,iquad,
     1      rfac0,nquad,slp_near)

            
      call cpu_time(t2)
      tquadgen = t2-t1

      ifinout = 1      


      call cpu_time(t1)
      call lpcomp_helm_comb_dir_bdry(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,inode_id,
     2  eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,slp_near,
     2  sigma,nfars,npts_over,ixyzso,srcover,wover,ifinout,pot)
      call cpu_time(t2)
      tlpcomp = t2-t1


      stop
      end

