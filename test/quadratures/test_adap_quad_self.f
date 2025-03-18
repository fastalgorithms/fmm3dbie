      subroutine test_adap_quad_self(nsuccess)
      implicit real *8 (a-h,o-z)
      real *8 c0(3)
      integer, allocatable :: norders(:), iptype(:), ixyzs(:)
      real *8, allocatable :: srcvals(:,:), srccoefs(:,:)
      complex *16 zk
      complex *16, allocatable :: znear_ggq(:), znear_adap(:)
      real *8, allocatable :: dnear_ggq(:), dnear_adap(:)

      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_src(:,:)

      integer row_ptr(2), col_ind(1), iquad(1)

      procedure (), pointer :: fker

      external h3d_slp, l3d_slp

      call prini(6,13)
      
      a = 1.0d0
      na = 1
      c0(1) = 0.0d0
      c0(2) = 0.0d0
      c0(3) = 0.0d0

      norder = 6
      iptype0 = 11

      npatches = 0
      npts = 0
      call get_sphere_npat_mem(a, na, c0, norder, iptype0, npatches, 
     1   npts)
      allocate(norders(npatches), iptype(npatches), ixyzs(npatches+1))
      allocate(srcvals(12,npts), srccoefs(9,npts))

      call get_sphere_npat(a, na, c0, norder, iptype0, npatches, npts,
     1  norders, ixyzs, iptype, srccoefs, srcvals)
      
      allocate(ipatch_id(npts), uvs_src(2,npts))

      call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts, 
     1  ipatch_id, uvs_src)

      npols = ixyzs(2) - ixyzs(1)
      allocate(znear_ggq(npols), znear_adap(npols))
      allocate(dnear_ggq(npols), dnear_adap(npols))

      npatches_use = 1
      ndd = 0
      ndz = 1
      ndi = 0

      ndtarg = 12
      ntarg = 1
      eps = 1.0d-9
      ipv = 0
      fker => h3d_slp
      zk = 1.1d0
      nnz = 1
      row_ptr(1) = 1
      row_ptr(2) = 2
      col_ind(1) = 1
      iquad(1) = 1
      rfac0 = 1.25d0
      call prinf('ipatch_id=*',ipatch_id,1)
      call prinf('iptype=*', iptype,1)
      call prin2('uvs_src=*',uvs_src,2)
      call zgetnearquad_ggq_guru(npatches_use, norders, ixyzs, 
     1  iptype, npts, srccoefs, srcvals, ndtarg, ntarg, srcvals,
     2  ipatch_id, uvs_src, eps, ipv, fker, ndd, dpars, ndz, zk, 
     3  ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, npols,
     4  znear_ggq)
      
      call zgetnearquad_adap_guru(npatches_use, norders, ixyzs, 
     1  iptype, npts, srccoefs, srcvals, ndtarg, ntarg, srcvals,
     2  ipatch_id, uvs_src, eps, fker, ndd, dpars, ndz, zk, 
     3  ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, npols, 
     4  znear_adap)
      

      call prin2('znear_ggq=*',znear_ggq,2*npols)
      call prin2('znear_adap=*',znear_adap,2*npols)
      
      erra = 0
      ra = 0
      do i=1,npols
        ra = ra + abs(znear_ggq(i))**2
        erra = erra + abs(znear_ggq(i) - znear_adap(i))**2
      enddo

      erra_abs = sqrt(erra)
      erra = sqrt(erra/ra)
      call prin2('relative error in adaptive integration =*',erra,1)
      call prin2('absolute error in adaptive integration =*',erra_abs,1)

      return
      end
