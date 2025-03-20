      implicit real *8 (a-h,o-z)
      real *8 rmid
      integer npars(3)
      real *8, allocatable :: ptinfo(:,:), srcvals(:,:), srccoefs(:,:)
      real *8, allocatable :: ptcoefs(:,:)
      integer, allocatable :: ixyzs(:), iptype(:), norders(:)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_src(:,:)
      complex *16, allocatable :: fints(:)
      complex *16, allocatable :: amat(:,:)
      complex *16 zk
      integer, allocatable :: row_ptr(:), col_ind(:), iquad(:)
      external h3d_slp_disk
      
      call prini(6,13)
      rmid = 0.3d0
      npars(1) = 3
      npars(2) = 5
      npars(3) = 5
      
      iort = 1
      

      norder = 8
      iptype0 = 11
      call mesh_circle_pts_mem(rmid, npars, iort, norder, iptype0, &
        npatches, npts)
      call prinf('npatches=*',npatches,1)
      call prinf('npts=*',npts,1)
      
      allocate(ptinfo(6,npts), srcvals(12,npts), srccoefs(9,npts))
      allocate(ixyzs(npatches+1), iptype(npatches), norders(npatches))
      call mesh_circle_pts(rmid, npars, iort, norder, iptype0, &
        npatches, npts, ixyzs, ptinfo)
      do i=1,npatches
        iptype(i) = iptype0
        norders(i) = norder
      enddo
      allocate(ptcoefs(6,npts))
      call surf_vals_to_coefs(6, npatches, norders, ixyzs, iptype, npts, &
        ptinfo, ptcoefs)

      do i=1,npts
        srcvals(1,i) = ptinfo(1,i)
        srcvals(2,i) = ptinfo(2,i)
        srcvals(3,i) = 0

        srcvals(4,i) = ptinfo(3,i)
        srcvals(5,i) = ptinfo(4,i)
        srcvals(6,i) = 0

        srcvals(7,i) = ptinfo(5,i)
        srcvals(8,i) = ptinfo(6,i)
        srcvals(9,i) = 0

        srcvals(10,i) = 0
        srcvals(11,i) = 0
        srcvals(12,i) = 1
        
        srccoefs(1,i) = ptcoefs(1,i)
        srccoefs(2,i) = ptcoefs(2,i)
        srccoefs(3,i) = 0

        srccoefs(4,i) = ptcoefs(3,i)
        srccoefs(5,i) = ptcoefs(4,i)
        srccoefs(6,i) = 0

        srccoefs(7,i) = ptcoefs(5,i)
        srccoefs(8,i) = ptcoefs(6,i)
        srccoefs(9,i) = 0
      enddo

      call prin2('srcvals=*',srcvals,24)
      call prin2('srccoefs=*',srccoefs,24)


!
!
!
      nptuse = 24
      allocate(ipatch_id(npts), uvs_src(2,npts))
      allocate(fints(nptuse))
      call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts, &
        ipatch_id, uvs_src)
      nnz = npatches*nptuse
      allocate(row_ptr(nptuse+1), col_ind(nnz), iquad(nnz+1))
      do i=1,nptuse+1
        row_ptr(i) = (i-1)*npatches+1
      enddo

      do i=1,nptuse
        do j=1,npatches
           ii = (i-1)*npatches+j
           col_ind(ii) = j
        enddo
      enddo

      npols = ixyzs(2) - ixyzs(1)
      do i=1,nnz+1
        iquad(i) = (i-1)*npols+1
      enddo
      nquad = iquad(nnz+1)-1
      call prinf('nptuse=*',nptuse,1)



      ndtarg = 12
      eps = 1.0d-8
      allocate(amat(npts,nptuse))
      zk = 2.1d0
      ndz = 1
      ndi = 0 
      ndd = 0
      rfac0 = 1.25d0
      ipv = 0
      call zgetnearquad_ggq_guru(npatches, norders, ixyzs, iptype, &
        npts, srcoefs, srcvals, ndtarg, nptuse, srcvals, ipatch_id, &
        uvs_src, eps, ipv, h3d_slp_disk, ndd, dpars, ndz, zk, ndi, ipars, &
        nnz, row_ptr, col_ind, iquad, rfac0, nquad, amat)
      
      do i=1,nptuse
        fints(i) = 0
        do j=1,npts
          fints(i) = fints(i) + amat(j,i)
        enddo
        write(33,*) real(fints(i)), imag(fints(i))
      enddo






      stop
      end
