      implicit real *8 (a-h,o-z) 
      implicit integer *8 (i-n)
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:),targs(:,:)
      real *8, allocatable :: wts(:)
      
      real *8 ts(2), rres
      real *8, allocatable :: rfacs(:,:), errs(:)
      character *100 fname
      integer *8 ipars(2)
      integer *8, allocatable :: row_ptr(:),col_ind(:)
      integer *8, allocatable :: iquad(:)
      real *8, allocatable :: srcover(:,:),wover(:)

      integer *8, allocatable :: norders(:),ixyzs(:),iptype(:)

      integer *8, allocatable :: ipatch_id(:),inode_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      real *8 c0(3)
      real *8 area, centroid(3), rmoi(3,3)
      real *8 forces1(3), torques1(3)
      real *8 trans_vels1(3), rot_vels1(3)
      real *8, allocatable :: uval(:,:), tracval(:,:), soln(:,:)
      complex * 16 zpars
      integer *8 int8_0,int8_3,int8_9
      integer *8 ncomp
      integer *8, allocatable :: icomps(:)

      call prini(6,13)

      int8_0 = 0
      int8_3 = 3
      int8_9 = 9
      done = 1
      pi = atan(done)*4

      a = 1.0d0
      na = 4
      c0(1) = 1.1d0
      c0(2) = 0.3d0
      c0(3) = 0.2d0

      iptype0 = 1

      norder = 4
      npols = (norder+1)*(norder+2)/2

      call get_sphere_npat_mem(a, na, c0, norder, iptype0, npatches, &
        npts)
      call prinf('npatches=*',npatches,1)
      call prinf('npts=*',npts,1)
      
      allocate(srcvals(12,npts), srccoefs(9,npts), ixyzs(npatches+1), &
        iptype(npatches), norders(npatches))
      allocate(wts(npts))
      call get_sphere_npat(a, na, c0, norder, iptype0, npatches, &
        npts, norders, ixyzs, iptype, srccoefs, srcvals)

      call get_qwts(npatches, norders, ixyzs, iptype, npts, srcvals, &
        wts)

      call get_surf_moments(npatches, norders, ixyzs, iptype, npts, &
        srcvals, wts, area, centroid, rmoi)
      
      call prin2('area=*',area,1)
      call prin2('centroid=*',centroid,3)
      call prin2('rmoi=*',rmoi,9)

      err_moi = rmoi(1,1) - 4.0d0/3.0d0*2*pi*(a**4)
      call prin2('err_moi=*', err_moi,1)

      forces1(1) = 1.0d0
      forces1(2) = 0.0d0
      forces1(3) = 0.0d0
      
      torques1(1) = 0.0d0
      torques1(2) = 0.0d0
      torques1(3) = 0.0d0

      numit = 100
      allocate(errs(numit+1))
      allocate(soln(3,npts))
      ncomp = 1
      allocate(icomps(ncomp+1))
      icomps(1) = 1
      icomps(2) = npatches+1

      eps = 1.0d-8


      call stok_s_mob_solver(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, srcvals, ncomp, icomps, eps, numit, forces1, &
        torques1, eps, niter, errs, rres, soln, trans_vels1, rot_vels1)
      
      call prin2('trans_vels1=*',trans_vels1,3)
      call prin2('rot_vels1=*',rot_vels1,3)

      erra = abs(trans_vels1(1) - 1.0d0/6/pi)
      call prin2('error in mobility matrix=*',erra,1)
      
      
      


      
      stop
      end


