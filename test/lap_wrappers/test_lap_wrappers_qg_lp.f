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
      integer nuv(2), iptype0
      integer ndd, ndz, ndi, nker, lwork, ndim
      real *8 work(1)


      call prini(6,13)

      done = 1
      pi = atan(done)*4


      norder = 3
      nuv(1) = 4
      nuv(2) = 12

      iptype0 = 1
      call get_stellarator_npat_mem(nuv, norder, iptype0, npatches,
     1  npts)

      allocate(srcvals(12,npts),srccoefs(9,npts))
      allocate(targs(3,npts))
      allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))
      allocate(ixyzso(npatches+1),nfars(npatches))

      call get_stellarator_npat(nuv, norder, iptype0, npatches, npts,
     1  norders, ixyzs, iptype, srccoefs, srcvals)

      xyz_in(1) = -4.501d0
      xyz_in(2) = 1.7d-3
      xyz_in(3) = 0.00001d0

      xyz_out(1) = -3.5d0
      xyz_out(2) = 3.1d0
      xyz_out(3) = 20.1d0

      allocate(wts(npts))
      call get_qwts(npatches, norders, ixyzs, iptype, npts, srcvals,
     1  wts)

      allocate(cms(3,npatches), rads(npatches), rad_near(npatches))
      allocate(pot(npts), potslp(npts), potdlp(npts))

      call get_centroid_rads(npatches, norders, ixyzs, iptype, npts, 
     1     srccoefs, cms, rads)

      allocate(sigma(npts), uval(npts), dudnval(npts))

      do i=1,npts
        call l3d_slp(xyz_out, 3, srcvals(1,i), 0, dpars, 0, zpars, 0,
     1    ipars, uval(i))
        call l3d_sprime(xyz_out, 12, srcvals(1,i), 0, dpars, 0, zpars,
     1    0, ipars, dudnval(i))
      enddo

     
      ndtarg = 3
      do i=1,npts
        targs(1,i) = srcvals(1,i)
        targs(2,i) = srcvals(2,i)
        targs(3,i) = srcvals(3,i)
      enddo

      allocate(ipatch_id(npts), uvs_targ(2,npts))
      do i=1,npts
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo

      call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts, 
     1         ipatch_id, uvs_targ)
 
c
c    find near field
c
      iptype = 1
      call get_rfacs(norder, iptype, rfac, rfac0)
      do i=1,npatches 
        rad_near(i) = rads(i)*rfac
      enddo
      

      call findnearmem(cms, npatches, rad_near, ndtarg, targs, npts,  
     1  nnz)

      allocate(row_ptr(npts+1), col_ind(nnz))
      
      call findnear(cms, npatches, rad_near, ndtarg, targs, npts,
     1  row_ptr, col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches, ixyzs, npts, nnz, row_ptr, col_ind,
     1  iquad)

      nquad = iquad(nnz+1)-1
      allocate(slp_near(nquad), dlp_near(nquad))



      eps = 0.50001d-6

      ikerorder = -1
      zk = 0

      call cpu_time(t1)
      call get_far_order(eps, npatches, norders, ixyzs, iptype, cms,
     1  rads, npts, srccoefs, ndtarg, npts, targs, ikerorder, zk,
     2  nnz, row_ptr, col_ind, rfac, nfars, ixyzso)
      call cpu_time(t2)
      tfar = t2-t1


      npts_over = ixyzso(npatches+1)-1

      print *, "npts_over=",npts_over


      allocate(srcover(12,npts_over), sigmaover(npts_over),
     1         wover(npts_over))

          
      call oversample_geom(npatches, norders, ixyzs, iptype, npts, 
     1  srccoefs, srcvals, nfars, ixyzso, npts_over, srcover)

      call get_qwts(npatches, nfars, ixyzso, iptype, npts_over, 
     1  srcover, wover)


      do i = 1,nquad
        slp_near(i) = 0
        dlp_near(i) = 0
      enddo

!      goto 1111

      call cpu_time(t1)

      dpars(1) = 1.0d0
      dpars(2) = 0.0d0

      iquadtype = 1

      call getnearquad_lap_comb_dir(npatches, norders,
     1  ixyzs, iptype, npts, srccoefs, srcvals, eps, dpars, iquadtype, 
     1  nnz, row_ptr, col_ind, iquad, rfac0, nquad, slp_near)

      
      dpars(1) = 0.0d0
      dpars(2) = 1.0d0

      call getnearquad_lap_comb_dir(npatches, norders,
     1  ixyzs, iptype, npts, srccoefs, srcvals, eps, dpars, iquadtype, 
     1  nnz, row_ptr, col_ind, iquad, rfac0, nquad, dlp_near)

      call cpu_time(t2)
      tquadgen = t2-t1

      ifinout = 1     

      dpars(1) = 1.0d0
      dpars(2) = 0.0d0

      ndd = 2
      ndz = 0
      ndi = 0
      nker = 1
      ndim = 1
      lwork = 0

      call cpu_time(t1)

      call lpcomp_lap_comb_dir_addsub(npatches, norders, ixyzs,
     1  iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars,
     2  ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, slp_near,
     3  nfars, npts_over, ixyzso, srcover, wover, lwork, work, ndim, 
     4  dudnval, potslp)


      dpars(1) = 0.0d0
      dpars(2) = 1.0d0

      call lpcomp_lap_comb_dir_addsub(npatches, norders, ixyzs,
     1  iptype, npts, srccoefs, srcvals, eps, ndd, dpars, ndz, zpars,
     2  ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, dlp_near,
     3  nfars, npts_over, ixyzso, srcover, wover, lwork, work, ndim, 
     4  uval, potdlp)

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
        errl2 = errl2 + abs(uval(i) - pot(i))**2*wts(i)
        rl2 = rl2 + abs(uval(i))**2*wts(i)
      enddo

      err = sqrt(errl2/rl2)

      call prin2('error in greens identity=*',err,1)

      i1 = 0
      if(err.lt.1.0d-2) i1 = 1

 1111 continue
      allocate(potslp2(npts))

      dpars(1) = 1.0d0
      dpars(2) = 0.0d0

      
      call lap_comb_dir_eval(npatches, norders, ixyzs,
     1  iptype, npts, srccoefs, srcvals, ndtarg, npts, targs, ipatch_id,
     2  uvs_targ, eps, dpars, dudnval, potslp2)


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



