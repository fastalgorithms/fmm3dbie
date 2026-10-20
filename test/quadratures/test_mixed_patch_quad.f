      subroutine test_mixed_patch_quad(nsuccess)
c
c------------------------------------------------------------------
c  This test verifies the "unique (norder,iptype) caching"
c  optimization in ?getnearquad_ggq_guru (see
c  src/quadratures/ggq-quads.f, get_ptype_uni/fill_disc_exps_uni)
c  by exercising it on a genuinely mixed-type, mixed-order mesh.
c
c
c  Strategy:
c    * Build three small sphere meshes of different (order,type):
c        block 1: iptype = 1  (triangle, RV nodes),   order n1
c        block 2: iptype = 11 (quad, GL nodes),        order n2
c        block 3: iptype = 12 (quad, cheb nodes),       order n1
c    * Concatenate the three blocks into a single srccoefs/srcvals/
c      norders/iptype/ixyzs mesh with npatches = npatches1+2+3.
c    * Build a near field with exactly one target per patch
c    * Call z/dgetnearquad_ggq_guru once on the combined mesh 
c      and compare to the result of separately calling the same 
c      routine on each patch in isolation 
c------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)

      integer *8 int8_1, int8_11, int8_12
      real *8 c01(3), c02(3), c03(3)

      integer *8, allocatable :: norders1(:), norders2(:), norders3(:)
      integer *8, allocatable :: iptype1(:), iptype2(:), iptype3(:)
      integer *8, allocatable :: ixyzs1(:), ixyzs2(:), ixyzs3(:)
      real *8, allocatable :: srcvals1(:,:), srcvals2(:,:)
      real *8, allocatable :: srcvals3(:,:)
      real *8, allocatable :: srccoefs1(:,:), srccoefs2(:,:)
      real *8, allocatable :: srccoefs3(:,:)

      integer *8, allocatable :: norders(:), iptype(:), ixyzs(:)
      real *8, allocatable :: srcvals(:,:), srccoefs(:,:)

      integer *8, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_src(:,:)

      integer *8, allocatable :: row_ptr(:), col_ind(:), iquad(:)
      integer *8 row_ptr1(2), col_ind1(1), iquad1(1)
      integer *8 ipatch_id1(1), ixyzs_loc(2), norders_loc(1)
      integer *8 iptype_loc(1)
      real *8, allocatable :: srccoefs_loc(:,:), srcvals_loc(:,:)
      real *8 targvals_loc(12), uvs_targ_loc(2)

      real *8, allocatable :: targvals(:,:)
      integer *8, allocatable :: ipatch_id_t(:)
      real *8, allocatable :: uvs_targ_t(:,:)

      complex *16, allocatable :: znear_all(:), znear_solo(:)
      real *8, allocatable :: dnear_all(:), dnear_solo(:)

      complex *16 zk
      procedure (), pointer :: fker
      external h3d_slp, l3d_slp

      call prini(6,13)

      int8_1 = 1
      int8_11 = 11
      int8_12 = 12

c     radii/centers chosen simply so the three blocks are disjoint
c     pieces of geometry -- the near field is specified explicitly
c     below (self targets only) and does not depend on the blocks
c     being adjacent in space

      a1 = 1.0d0
      a2 = 1.3d0
      a3 = 0.7d0

      n1 = 4
      n2 = 6

      na = 1

      c01(1) = 0.0d0
      c01(2) = 0.0d0
      c01(3) = 0.0d0

      c02(1) = 5.0d0
      c02(2) = 0.0d0
      c02(3) = 0.0d0

      c03(1) = -5.0d0
      c03(2) = 0.0d0
      c03(3) = 0.0d0

c
c   block 1: triangles, order n1
c
      call get_sphere_npat_mem(a1, na, c01, n1, int8_1, npatches1,
     1   npts1)
      allocate(norders1(npatches1), iptype1(npatches1))
      allocate(ixyzs1(npatches1+1))
      allocate(srcvals1(12,npts1), srccoefs1(9,npts1))
      call get_sphere_npat(a1, na, c01, n1, int8_1, npatches1, npts1,
     1   norders1, ixyzs1, iptype1, srccoefs1, srcvals1)

c
c   block 2: GL quads, order n2
c
      call get_sphere_npat_mem(a2, na, c02, n2, int8_11, npatches2,
     1   npts2)
      allocate(norders2(npatches2), iptype2(npatches2))
      allocate(ixyzs2(npatches2+1))
      allocate(srcvals2(12,npts2), srccoefs2(9,npts2))
      call get_sphere_npat(a2, na, c02, n2, int8_11, npatches2, npts2,
     1   norders2, ixyzs2, iptype2, srccoefs2, srcvals2)

c
c   block 3: cheb quads, order n1 (same order as block 1, different
c   type -- exercises a third, previously-absent cache key)
c
      call get_sphere_npat_mem(a3, na, c03, n1, int8_12, npatches3,
     1   npts3)
      allocate(norders3(npatches3), iptype3(npatches3))
      allocate(ixyzs3(npatches3+1))
      allocate(srcvals3(12,npts3), srccoefs3(9,npts3))
      call get_sphere_npat(a3, na, c03, n1, int8_12, npatches3, npts3,
     1   norders3, ixyzs3, iptype3, srccoefs3, srcvals3)

      print *, "npatches1 (tri, order n1)      =", npatches1
      print *, "npatches2 (GL quad, order n2)  =", npatches2
      print *, "npatches3 (cheb quad, order n1)=", npatches3

      npatches = npatches1 + npatches2 + npatches3
      npts = npts1 + npts2 + npts3

      allocate(norders(npatches), iptype(npatches))
      allocate(ixyzs(npatches+1))
      allocate(srcvals(12,npts), srccoefs(9,npts))

c
c    concatenate the three blocks into one mesh
c
      do i=1,npatches1
        norders(i) = norders1(i)
        iptype(i) = iptype1(i)
      enddo
      do i=1,npatches2
        norders(npatches1+i) = norders2(i)
        iptype(npatches1+i) = iptype2(i)
      enddo
      do i=1,npatches3
        norders(npatches1+npatches2+i) = norders3(i)
        iptype(npatches1+npatches2+i) = iptype3(i)
      enddo

      do i=1,npatches1+1
        ixyzs(i) = ixyzs1(i)
      enddo
      ioff = ixyzs(npatches1+1) - 1
      do i=1,npatches2+1
        ixyzs(npatches1+i) = ioff + ixyzs2(i)
      enddo
      ioff = ixyzs(npatches1+npatches2+1) - 1
      do i=1,npatches3+1
        ixyzs(npatches1+npatches2+i) = ioff + ixyzs3(i)
      enddo

      do i=1,npts1
        do j=1,12
          srcvals(j,i) = srcvals1(j,i)
        enddo
        do j=1,9
          srccoefs(j,i) = srccoefs1(j,i)
        enddo
      enddo

      do i=1,npts2
        do j=1,12
          srcvals(j,npts1+i) = srcvals2(j,i)
        enddo
        do j=1,9
          srccoefs(j,npts1+i) = srccoefs2(j,i)
        enddo
      enddo

      do i=1,npts3
        do j=1,12
          srcvals(j,npts1+npts2+i) = srcvals3(j,i)
        enddo
        do j=1,9
          srccoefs(j,npts1+npts2+i) = srccoefs3(j,i)
        enddo
      enddo

      allocate(ipatch_id(npts), uvs_src(2,npts))
      call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts,
     1  ipatch_id, uvs_src)

c
c    build a near field with one target per patch: the target is
c    the patch's own first discretization node (a self-interaction,
c    forcing use of the self-quadrature branch), so ntarg = npatches
c    and nnz = npatches (row i only sees column i).
c
      ntarg = npatches
      nnz = npatches

      allocate(row_ptr(ntarg+1), col_ind(nnz), iquad(nnz+1))
      do i=1,ntarg
        row_ptr(i) = i
        col_ind(i) = i
      enddo
      row_ptr(ntarg+1) = ntarg+1

      iquad(1) = 1
      do i=1,nnz
        ip = col_ind(i)
        npols = ixyzs(ip+1)-ixyzs(ip)
        iquad(i+1) = iquad(i) + npols
      enddo
      nquad = iquad(nnz+1)-1

      allocate(znear_all(nquad), dnear_all(nquad))
      allocate(znear_solo(nquad), dnear_solo(nquad))

      ndd = 0
      ndz = 1
      ndi = 0
      ipv = 0
      eps = 1.0d-9
      rfac0 = 1.25d0
      zk = 1.1d0

      ndtarg = 12

c
c    targets: self node of each patch, i.e. srcvals(:,ixyzs(ip)),
c    which is exactly what get_patch_id_uvs/ixyzs already index, so
c    we can just reuse srcvals(1,ixyzs(ip)) as the target list --
c    build a compact targvals array with one column per patch.
c
      allocate(targvals(ndtarg,ntarg))
      allocate(ipatch_id_t(ntarg), uvs_targ_t(2,ntarg))
      do ip=1,npatches
        istart = ixyzs(ip)
        do j=1,ndtarg
          targvals(j,ip) = srcvals(j,istart)
        enddo
        ipatch_id_t(ip) = ipatch_id(istart)
        uvs_targ_t(1,ip) = uvs_src(1,istart)
        uvs_targ_t(2,ip) = uvs_src(2,istart)
      enddo

      fker => h3d_slp
      call zgetnearquad_ggq_guru(npatches, norders, ixyzs, iptype,
     1  npts, srccoefs, srcvals, ndtarg, ntarg, targvals,
     2  ipatch_id_t, uvs_targ_t, eps, ipv, fker, ndd, dpars, ndz, zk,
     3  ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad,
     4  znear_all)

      fker => l3d_slp
      call dgetnearquad_ggq_guru(npatches, norders, ixyzs, iptype,
     1  npts, srccoefs, srcvals, ndtarg, ntarg, targvals,
     2  ipatch_id_t, uvs_targ_t, eps, ipv, fker, ndd, dpars, ndz, zk,
     3  ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, nquad,
     4  dnear_all)

c
c    now recompute the same self-quadrature entries one patch at a
c    time (npatches_use = 1 slices, exactly as test_adap_quad_self.f
c    does).  Each such call only ever builds nuni = 1 cache entry,
c    so it cannot suffer from any multi-entry indexing bug -- this
c    is the reference to compare against.
c
      row_ptr1(1) = 1
      row_ptr1(2) = 2
      col_ind1(1) = 1
      iquad1(1) = 1
      ipatch_id1(1) = 1
      ixyzs_loc(1) = 1

c     each "solo" call below is built on freshly-sized local buffers
c     holding only patch ip's own points, renumbered to start at
c     local index 1 -- this avoids any pointer-arithmetic aliasing
c     between the global mesh's point offsets and the small
c     (npatches=1) array shapes the callee expects, so the only
c     thing that can differ between this call and the combined call
c     above is whether the vals-to-coefs cache had 1 or nuni=3
c     entries when norder/iptype for this exact patch was looked up.

      do ip=1,npatches
        npols = ixyzs(ip+1)-ixyzs(ip)
        istart = ixyzs(ip)

        norders_loc(1) = norders(ip)
        iptype_loc(1) = iptype(ip)
        ixyzs_loc(2) = npols+1

        allocate(srccoefs_loc(9,npols), srcvals_loc(12,npols))
        do i=1,npols
          do j=1,9
            srccoefs_loc(j,i) = srccoefs(j,istart+i-1)
          enddo
          do j=1,12
            srcvals_loc(j,i) = srcvals(j,istart+i-1)
          enddo
        enddo

        do j=1,ndtarg
          targvals_loc(j) = targvals(j,ip)
        enddo
        uvs_targ_loc(1) = uvs_targ_t(1,ip)
        uvs_targ_loc(2) = uvs_targ_t(2,ip)

        fker => h3d_slp
        call zgetnearquad_ggq_guru(int8_1, norders_loc, ixyzs_loc,
     1    iptype_loc, npols, srccoefs_loc, srcvals_loc,
     2    ndtarg, int8_1, targvals_loc, ipatch_id1,
     3    uvs_targ_loc, eps, ipv, fker, ndd, dpars, ndz, zk,
     4    ndi, ipars, int8_1, row_ptr1, col_ind1, iquad1,
     5    rfac0, npols, znear_solo(iquad(ip)))

        fker => l3d_slp
        call dgetnearquad_ggq_guru(int8_1, norders_loc, ixyzs_loc,
     1    iptype_loc, npols, srccoefs_loc, srcvals_loc,
     2    ndtarg, int8_1, targvals_loc, ipatch_id1,
     3    uvs_targ_loc, eps, ipv, fker, ndd, dpars, ndz, zk,
     4    ndi, ipars, int8_1, row_ptr1, col_ind1, iquad1,
     5    rfac0, npols, dnear_solo(iquad(ip)))

        deallocate(srccoefs_loc, srcvals_loc)
      enddo

      erra_z = 0
      ra_z = 0
      erra_d = 0
      ra_d = 0
      do i=1,nquad
        erra_z = erra_z + abs(znear_all(i)-znear_solo(i))**2
        ra_z = ra_z + abs(znear_solo(i))**2
        erra_d = erra_d + (dnear_all(i)-dnear_solo(i))**2
        ra_d = ra_d + dnear_solo(i)**2
      enddo
      erra_z = sqrt(erra_z/ra_z)
      erra_d = sqrt(erra_d/ra_d)

      call prin2('relative error, complex mixed-type near quad =*',
     1   erra_z, 1)
      call prin2('relative error, real mixed-type near quad =*',
     1   erra_d, 1)

      nsuccess = 0
      if(erra_z.lt.1.0d-12) nsuccess = nsuccess + 1
      if(erra_d.lt.1.0d-12) nsuccess = nsuccess + 1

      return
      end
