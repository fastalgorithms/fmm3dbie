
      implicit real *8 (a-h,o-z) 

      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      character *100 fname
      integer ipars(2)

      integer, allocatable :: norders(:), ixyzs(:), iptype(:)

      real *8 abc(3), c0(3)
      integer nabc(3)

      integer ipatch_id
      real *8 uvs_targ(2)

      call prini(6,13)


      done = 1
      pi = atan(done)*4

      iptype0 = 1
      dlam = 0.5d0

      abc(1) = 2.1d0
      abc(2) = 1.0d0
      abc(3) = 4.0d0

      nabc(1) = 5
      nabc(2) = 3
      nabc(3) = 10

      c0(1) = 0
      c0(2) = 5
      c0(3) = -3

      npatches = 0
      npts = 0

      norder = 4

      call get_ellipsoid_npat_mem(abc, nabc, c0, norder, iptype0, 
     1  npatches, npts)
      print *, "npatches tri=", npatches

      allocate(srcvals(12,npts), srccoefs(9,npts))
      allocate(norders(npatches), ixyzs(npatches+1), iptype(npatches))

      call get_ellipsoid_npat(abc, nabc, c0, norder, iptype0, 
     1  npatches, npts, norders, ixyzs, iptype, srccoefs, srcvals)

      call plot_surface_info_all(dlam, npatches, norders, ixyzs, iptype,
     1  npts, srccoefs, srcvals, 'ellip_tri.vtk','a')

      call surf_quadratic_msh_vtk_plot(npatches, norders, ixyzs, iptype,
     1   npts, srccoefs, srcvals, 'ellip_tri_msh.vtk','a')


      deallocate(srcvals, srccoefs, norders, ixyzs, iptype)

      npatches = 0
      npts = 0
      iptype0 = 11

      call get_ellipsoid_npat_mem(abc, nabc, c0, norder, iptype0, 
     1  npatches, npts) 
      print *, "npatches quad=", npatches

      allocate(srcvals(12,npts), srccoefs(9,npts))
      allocate(norders(npatches), ixyzs(npatches+1), iptype(npatches))

      call get_ellipsoid_npat(abc, nabc, c0, norder, iptype0, 
     1  npatches, npts, norders, ixyzs, iptype, srccoefs, srcvals)

      call plot_surface_info_all(dlam, npatches, norders, ixyzs, iptype,
     1  npts, srccoefs, srcvals, 'ellip_quad.vtk','a')

      stop
      end


   




