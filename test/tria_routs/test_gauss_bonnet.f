      implicit real *8 (a-h,o-z) 
      implicit integer *8 (i-n)
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:),targs(:,:)
      real *8, allocatable :: wts(:)
      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)
      real *8 errs(6),ts(2)
      real *8, allocatable :: rfacs(:,:)
      complex *16 zk
      character *100 fname
      integer *8 ipars(2)

      integer *8, allocatable :: norders(:),ixyzs(:),iptype(:)

      integer *8, allocatable :: ipatch_id(:),inode_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      real *8 xyz_out(3),xyz_in(3)
      real *8 dpars(2)
      integer *8 nuv(2), iptype0
      integer *8 ndd, ndz, ndi, nker, lwork, ndim
      real *8 work(1)

      real *8, allocatable :: srcvals2(:,:), srccoefs2(:,:)
      real *8, allocatable :: sfform(:,:,:)


      call prini(6,13)

      done = 1
      pi = atan(done)*4


      norder = 8
      nuv(1) = 8
      nuv(2) = 24

      iptype0 = 1
      call get_stellarator_npat_mem(nuv, norder, iptype0, npatches,
     1  npts)

      allocate(srcvals(12,npts),srccoefs(9,npts),wts(npts))
      allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))

      call get_stellarator_npat(nuv, norder, iptype0, npatches, npts,
     1  norders, ixyzs, iptype, srccoefs, srcvals)
      call get_qwts(npatches, norders, ixyzs, iptype, npts, srcvals, 
     1  wts)

      ntarg = 1
      xyz_in(1) = -4.501d0
      xyz_in(2) = 1.7d-3
      xyz_in(3) = 0.00001d0

c
c  test 
c
      allocate(srccoefs2(18,npts), srcvals2(30,npts))
      allocate(sfform(2,2,npts))
      
      call get_second_derivative_surfinfo(npatches, norders, 
     1  ixyzs, iptype, npts, srccoefs, srcvals, srccoefs2, srcvals2)

      call get_second_fundamental_form(npatches, norders, ixyzs,
     1  iptype, npts, srccoefs, srcvals, sfform)
      rint = 0
      rint2 = 0
      do i=1,npts
        rint = rint + srcvals2(29,i)*srcvals2(30,i)*wts(i)
      enddo

      print *, "rint = ",rint

      rint = 0
      eps = 1.0d-7
      istrat = 1
      intype = 1
      np0 = 1


      do ipat = 1,npatches
        npols = ixyzs(ipat+1) - ixyzs(ipat)
        ndsc = 18
c        call dtriaints(eps, istart, intype, np0, norders(ipat), 
c     1    npols, ndsc, )

      enddo
      
      stop
      end



