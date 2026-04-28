      implicit real *8 (a-h,o-z) 
      implicit integer *8 (i-n)
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:)
      integer *8, allocatable :: norders(:),ixyzs(:),iptype(:)
      
      real *8, allocatable :: srcvals1(:,:),srccoefs1(:,:)
      real *8, allocatable :: srccoefs2(:,:)
      real *8, allocatable :: wts1(:)
      integer *8, allocatable :: norders1(:),ixyzs1(:),iptype1(:)

      real *8 ts(2), rres
      real *8, allocatable :: errs(:)
      character *100 fname

      real *8 c0(3), abc(3)
      real *8 nabc(3)
      real *8, allocatable :: forces(:,:), torques(:,:)
      real *8, allocatable :: trans_vels_ex(:,:), rot_vels_ex(:,:)
      real *8, allocatable :: trans_vels(:,:), rot_vels(:,:)
      real *8, allocatable :: centroids(:,:), shifts(:,:), rmois(:,:,:)
      real *8, allocatable :: soln_mob(:,:), soln_res(:,:), rhs_res(:,:)
      complex * 16 zpars
      integer *8 int8_0,int8_3,int8_9
      integer *8 ncomp
      integer *8, allocatable :: icomps(:)
      real *8 xdiff(3), rtmp(3), dpars(2) 

      call prini(6,13)

      int8_0 = 0
      int8_3 = 3
      int8_9 = 9
      done = 1
      pi = atan(done)*4

      nlatx = 2
      nlaty = 1
      nlatz = 1

      dx = 4.0d0
      dy = 4.0d0
      dz = 4.0d0

      dh = 0.01d0*0

      ncomp = nlatx*nlaty*nlatz

      allocate(shifts(3,ncomp))
      allocate(trans_vels_ex(3,ncomp), rot_vels_ex(3,ncomp))
      do ilatz = 1,nlatz
        do ilaty = 1,nlaty
          do ilatx = 1,nlatx
            icomp = (ilatz-1)*nlatx*nlaty + (ilaty-1)*nlatx + ilatx
            shifts(1,icomp) = (ilatx-1)*dx + 2*dh*(hkrand(0)-0.5d0)
            shifts(2,icomp) = (ilaty-1)*dy + 2*dh*(hkrand(0)-0.5d0)
            shifts(3,icomp) = (ilatz-1)*dz + 2*dh*(hkrand(0)-0.5d0)

            trans_vels_ex(1,icomp) = hkrand(0)-0.5d0
            trans_vels_ex(2,icomp) = (hkrand(0)-0.5d0)
            trans_vels_ex(3,icomp) = (hkrand(0)-0.5d0)

            rot_vels_ex(1,icomp) = (hkrand(0)-0.5d0)
            rot_vels_ex(2,icomp) = (hkrand(0)-0.5d0)
            rot_vels_ex(3,icomp) = (hkrand(0)-0.5d0)
          enddo
        enddo
      enddo

      call prin2('shifts=*',shifts,3*ncomp)

!
!  extract shape of geometry to be extracted
!
      a = 1.0d0
      na = 2
      c0(1) = 0
      c0(2) = 0
      c0(3) = 0

      iptype0 = 1

      norder = 7
      npols = (norder+1)*(norder+2)/2
      
      npatches1 = 0
      npts1 = 0
      call get_sphere_npat_mem(a, na, c0, norder, iptype0, npatches1, &
        npts1)
      call prinf('npatches=*',npatches1,1)
      call prinf('npts=*',npts1,1)
      
      allocate(srcvals1(12,npts1), srccoefs1(9,npts1), &
        ixyzs1(npatches1+1), iptype1(npatches1), norders1(npatches1))
      call get_sphere_npat(a, na, c0, norder, iptype0, npatches1, &
        npts1, norders1, ixyzs1, iptype1, srccoefs1, srcvals1)

      allocate(icomps(ncomp+1))
      icomps(1) = 1
      icomps(2) = npatches+1
      npatches = ncomp*npatches1
      npts = ncomp*npts1
!
!
!
      do icomp=1,ncomp+1
        icomps(icomp) = npatches1*(icomp-1)+1
      enddo
!
!  Define ixyzs
!
      allocate(ixyzs(npatches+1), iptype(npatches), norders(npatches))
      allocate(srcvals(12,npts),srccoefs(9,npts))
      do icomp = 1,ncomp
        ipstart = icomps(icomp)
        do j=1,npatches1
          ixyzs(ipstart+j-1) = ixyzs1(j) + (icomp-1)*npts1
          norders(ipstart+j-1) = norders1(j)
          iptype(ipstart+j-1) = iptype1(j)
        enddo
        istart = (icomp-1)*npts1 + 1
        do j=1,npts1
           srcvals(1:12,istart+j-1) = srcvals1(1:12,j) 
           srcvals(1:3,istart+j-1) = srcvals(1:3,istart+j-1) + &
                     shifts(1:3,icomp)
        enddo
      enddo


      ixyzs(npatches+1) = npts+1
      allocate(srccoefs2(9,npts))
      int8_9 = 9
      call surf_vals_to_coefs(int8_9, npatches, norders, ixyzs, iptype, &
       npts, srcvals(1:9,1:npts), srccoefs)
      allocate(wts(npts))
      call get_qwts(npatches, norders, ixyzs, iptype, npts, srcvals, &
        wts)


      allocate(rhs_res(3,npts),soln_res(3,npts))
      allocate(centroids(3,ncomp),rmois(3,3,ncomp))
      do icomp=1,ncomp
        ipstart = icomps(icomp)
        ipend = icomps(icomp+1)-1
        istart = ixyzs(ipstart)
        iend = ixyzs(ipend+1)-1
        nptsloc = iend - istart + 1
        call get_surf_moments(npatches1, norders(ipstart), ixyzs1, &
          iptype(ipstart), npts1, srcvals(1,istart), wts(istart), &
          area, centroids(1,icomp), rmois(1,1,icomp))
        do i=istart,iend
          xdiff(1:3) = srcvals(1:3,i) - centroids(1:3,icomp)
          call cross_prod3d(xdiff, rot_vels_ex(1:3,icomp), rtmp)
          rhs_res(1:3,i) = trans_vels_ex(1:3,icomp) + rtmp(1:3)
        enddo
      enddo
      call prin2('centroids=*',centroids,3*ncomp)
      call prin2('shifts=*',shifts,3*ncomp)
      call prin2('rmois=*',rmois,9*ncomp)

      numit = 100
      allocate(errs(numit+1))
      eps = 1.0d-8
      eps_gmres = 1.0d-8
      ifinout = 1
      dpars(1) = 1.0d0
      dpars(2) = 2.0d0

      call stok_comb_vel_solver(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, srcvals, eps, dpars, numit, ifinout, rhs_res, &
        eps_gmres, niter, errs, rres, soln_res)
      call prin2('soln_res=*',soln_res,24)
      allocate(forces(3,ncomp), torques(3,ncomp))
      allocate(trans_vels(3,ncomp), rot_vels(3,ncomp))
      
      do icomp=1,ncomp
        ipstart = icomps(icomp)
        ipend = icomps(icomp+1)-1
        istart = ixyzs(ipstart)
        iend = ixyzs(ipend+1)-1
        nptsloc = iend - istart + 1
        forces(1:3,icomp) = 0
        torques(1:3,icomp) = 0
        trans_vels(1:3,icomp) = 0
        rot_vels(1:3,icomp) = 0
        do i=istart,iend
          forces(1:3,icomp) = forces(1:3,icomp) + &
            soln_res(1:3,i)*wts(i)
          xdiff(1:3) = srcvals(1:3,i) - centroids(1:3,icomp)
          call cross_prod3d(xdiff, soln_res(1,i), rtmp)
          torques(1:3,icomp) = torques(1:3,icomp) + rtmp(1:3)*wts(i)
        enddo
      enddo
      call prin2('forces=*',forces,3*ncomp)
      call prin2('torques=*',torques,3*ncomp)



      allocate(soln_mob(3,npts))


      call stok_s_mob_solver(npatches, norders, ixyzs, iptype, npts, &
        srccoefs, srcvals, ncomp, icomps, eps, numit, forces, &
        torques, eps, niter, errs, rres, soln_mob, trans_vels, rot_vels)
      
      print *, "============================="
      call prin2('trans_vels1=*',trans_vels,3*ncomp)
      call prin2('trans_vels_ex=*',trans_vels_ex,3*ncomp)
      print *, " "
      print *, "============================="
      print *, " "
      call prin2('rot_vels1=*',rot_vels,3*ncomp)
      call prin2('rot_vels_ex=*',rot_vels_ex,3*ncomp)
      print *, " "
      print *, "============================="


      
      stop
      end


