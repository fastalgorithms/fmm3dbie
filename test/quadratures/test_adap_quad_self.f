      subroutine test_adap_quad_self(nsuccess)
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8 int8_1, int8_11, int8_12

      int8_1 = 1
      int8_11 = 11
      int8_12 = 12

      na = 1
      norder = 6

      print *, "================"
      print *, "starting triangle test"
      call zcompute_errors(int8_1, na, norder, erra_abs_tri1)
      print *, ""
      print *, ""
      print *, "================"
      print *, "starting quadrangle test with Legendre points"
      call zcompute_errors(int8_11, na, norder, erra_abs_lege_poly1)
      print *, ""
      print *, ""
      print *, "================"
      print *, "starting quadrangle test with Chebyshev points"
      call zcompute_errors(int8_12, na, norder, erra_abs_cheb_poly1)
      print *, ""
      print *, ""


      na = 2
      print *, "================"
      print *, "starting triangle test"
      call zcompute_errors(int8_1, na, norder, erra_abs_tri2)
      print *, ""
      print *, ""
      print *, "================"
      print *, "starting quadrangle test with Legendre points"
      call zcompute_errors(int8_11, na, norder, erra_abs_lege_poly2)
      print *, ""
      print *, ""
      print *, "================"
      print *, "starting quadrangle test with Chebyshev points"
      call zcompute_errors(int8_12, na, norder, erra_abs_cheb_poly2)
      print *, ""
      print *, ""
      print *, "================"

      rat_tri = erra_abs_tri1/erra_abs_tri2
      rat_lege = erra_abs_lege_poly1/erra_abs_lege_poly2
      rat_cheb = erra_abs_cheb_poly1/erra_abs_cheb_poly2

      na = 1
      print *, "================"
      print *, "starting dtriangle test"
      call dcompute_errors(int8_1, na, norder, erra_abs_tri1)

      print *, ""
      print *, ""
      print *, "================"
      print *, "starting dquadrangle test with Legendre points"
      call dcompute_errors(int8_11, na, norder, erra_abs_lege_poly1)
      print *, ""
      print *, ""
      print *, "================"
      print *, "starting dquadrangle test with Chebyshev points"
      call dcompute_errors(int8_12, na, norder, erra_abs_cheb_poly1)
      print *, ""
      print *, ""


      na = 2
      print *, "================"
      print *, "starting dtriangle test"
      call dcompute_errors(int8_1, na, norder, erra_abs_tri2)
      print *, ""
      print *, ""
      print *, "================"
      print *, "starting dquadrangle test with Legendre points"
      call dcompute_errors(int8_11, na, norder, erra_abs_lege_poly2)
      print *, ""
      print *, ""
      print *, "================"
      print *, "starting dquadrangle test with Chebyshev points"
      call dcompute_errors(int8_12, na, norder, erra_abs_cheb_poly2)
      print *, ""
      print *, ""
      print *, "================"

      drat_tri = erra_abs_tri1/erra_abs_tri2
      drat_lege = erra_abs_lege_poly1/erra_abs_lege_poly2
      drat_cheb = erra_abs_cheb_poly1/erra_abs_cheb_poly2

      print *, "ratio tri =", rat_tri
      print *, "ratio lege =", rat_lege
      print *, "ratio cheb =", rat_cheb

      print *, "d ratio tri =", drat_tri
      print *, "d ratio lege =", drat_lege
      print *, "d ratio cheb =", drat_cheb

      nsuccess = 0
      if(rat_tri.ge.2**(norder-1)) nsuccess = nsuccess + 1
      if(rat_lege.ge.2**(norder-1)) nsuccess = nsuccess + 1
      if(rat_cheb.ge.2**(norder-1)) nsuccess = nsuccess + 1
      if(drat_tri.ge.2**(norder-1)) nsuccess = nsuccess + 1
      if(drat_lege.ge.2**(norder-1)) nsuccess = nsuccess + 1
      if(drat_cheb.ge.2**(norder-1)) nsuccess = nsuccess + 1
      print *, "nsuccess=", nsuccess



      return
      end



      subroutine zcompute_errors(iptype0, na, norder, erra_abs) 
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      real *8 c0(3)
      integer *8, allocatable :: norders(:), iptype(:), ixyzs(:)
      real *8, allocatable :: srcvals(:,:), srccoefs(:,:)
      complex *16 zk
      complex *16, allocatable :: znear_ggq(:), znear_adap(:)

      integer *8, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_src(:,:)
      real *8, allocatable :: pols(:,:)
      complex *16, allocatable :: fints(:), fints_ex(:)

      integer *8 row_ptr(2), col_ind(1), iquad(1)
      integer *8, allocatable :: iind2p(:,:)

      procedure (), pointer :: fker

      external h3d_slp

      call prini(6,13)
      
      a = 1.0d0
      c0(1) = 0.0d0
      c0(2) = 0.0d0
      c0(3) = 0.0d0

      norder = 6
      npowmax = (norder+1)*(norder+1)
      allocate(iind2p(2,npowmax))

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

      ipt = 5


      call zgetnearquad_ggq_guru(npatches_use, norders, ixyzs, 
     1  iptype, npts, srccoefs, srcvals, ndtarg, ntarg, srcvals(1,ipt),
     2  ipatch_id(ipt), uvs_src(1,ipt), eps, ipv, fker, ndd, dpars, 
     3  ndz, zk, ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, npols,
     4  znear_ggq)
      
      call zgetnearquad_adap_guru(npatches_use, norders, ixyzs, 
     1  iptype, npts, srccoefs, srcvals, ndtarg, ntarg, srcvals(1,ipt),
     2  ipatch_id(ipt), uvs_src(1,ipt), eps, fker, ndd, dpars, ndz, zk, 
     3  ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, npols, 
     4  znear_adap)
      
      
      if(iptype0.eq.1) then
      
        erra = 0
        ra = 0
        do i=1,npols
          erra = erra + abs(znear_ggq(i) - znear_adap(i))**2
        enddo

        erra_abs = sqrt(erra)
      endif

      if(iptype0.eq.11.or.iptype0.eq.12) then
        call polytens_ind2pow_2d(norder, "F", iind2p)
        allocate(pols(npols, npols))
        do i=1,npols
          call get_basis_pols(uvs_src(1,i), norder, npols, iptype0,
     1       pols(1,i))
        enddo

        allocate(fints(npols), fints_ex(npols))

        do ipol = 1,npols
          fints(ipol) = 0
          fints_ex(ipol) = 0
          do j=1,npols
            fints(ipol) = fints(ipol) + pols(ipol,j)*znear_adap(j)
            fints_ex(ipol) = fints_ex(ipol) + pols(ipol,j)*znear_ggq(j)
          enddo
        enddo

        erra_abs = 0.0d0
        do i=1,npols
          ideg = iind2p(1,i) + iind2p(2,i)
          if(ideg.le.norder) then 
             erra = abs(fints(i) - fints_ex(i)) 
             if(erra.gt.erra_abs) erra_abs = erra
          endif
        enddo
      endif
      call prin2('absolute error in adaptive integration =*',
     1     erra_abs,1)


      return
      end




      subroutine dcompute_errors(iptype0, na, norder, erra_abs) 
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      real *8 c0(3)
      integer *8, allocatable :: norders(:), iptype(:), ixyzs(:)
      real *8, allocatable :: srcvals(:,:), srccoefs(:,:)
      complex *16 zk
      real *8, allocatable :: dnear_ggq(:), dnear_adap(:)

      integer *8, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_src(:,:)
      real *8, allocatable :: pols(:,:)
      real *8, allocatable :: fints(:), fints_ex(:)

      integer *8 row_ptr(2), col_ind(1), iquad(1)
      integer *8, allocatable :: iind2p(:,:)

      procedure (), pointer :: fker
      character *100 fname

      external l3d_slp

      call prini(6,13)
      
      a = 1.0d0
      c0(1) = 0.0d0
      c0(2) = 0.0d0
      c0(3) = 0.0d0

      norder = 6
      npowmax = (norder+1)*(norder+1)
      allocate(iind2p(2,npowmax))

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
      allocate(dnear_ggq(npols), dnear_adap(npols))

      npatches_use = 1
      ndd = 0
      ndz = 0
      ndi = 0

      ndtarg = 12
      ntarg = 1
      eps = 1.0d-9
      ipv = 0
      fker => l3d_slp
      zk = 1.1d0
      nnz = 1
      row_ptr(1) = 1
      row_ptr(2) = 2
      col_ind(1) = 1
      iquad(1) = 1
      rfac0 = 1.25d0

      ipt = 5


      call dgetnearquad_ggq_guru(npatches_use, norders, ixyzs, 
     1  iptype, npts, srccoefs, srcvals, ndtarg, ntarg, srcvals(1,ipt),
     2  ipatch_id(ipt), uvs_src(1,ipt), eps, ipv, fker, ndd, dpars, 
     3  ndz, zk, ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, npols,
     4  dnear_ggq)
      
      call dgetnearquad_adap_guru(npatches_use, norders, ixyzs, 
     1  iptype, npts, srccoefs, srcvals, ndtarg, ntarg, srcvals(1,ipt),
     2  ipatch_id(ipt), uvs_src(1,ipt), eps, fker, ndd, dpars, ndz, zk, 
     3  ndi, ipars, nnz, row_ptr, col_ind, iquad, rfac0, npols, 
     4  dnear_adap)
      
      
      if(iptype0.eq.1) then
      
        erra = 0
        ra = 0
        do i=1,npols
          erra = erra + abs(dnear_ggq(i) - dnear_adap(i))**2
        enddo

        erra_abs = sqrt(erra)
      endif

      if(iptype0.eq.11.or.iptype0.eq.12) then
        call polytens_ind2pow_2d(norder, "F", iind2p)
        allocate(pols(npols, npols))
        do i=1,npols
          call get_basis_pols(uvs_src(1,i), norder, npols, iptype0,
     1       pols(1,i))
        enddo

        allocate(fints(npols), fints_ex(npols))

        do ipol = 1,npols
          fints(ipol) = 0
          fints_ex(ipol) = 0
          do j=1,npols
            fints(ipol) = fints(ipol) + pols(ipol,j)*dnear_adap(j)
            fints_ex(ipol) = fints_ex(ipol) + pols(ipol,j)*dnear_ggq(j)
          enddo
        enddo

        erra_abs = 0.0d0
        do i=1,npols
          ideg = iind2p(1,i) + iind2p(2,i)
          if(ideg.le.norder) then 
             erra = abs(fints(i) - fints_ex(i)) 
             if(erra.gt.erra_abs) erra_abs = erra
          endif
        enddo
      endif
      call prin2('absolute error in adaptive integration =*',
     1     erra_abs,1)


      return
      end
