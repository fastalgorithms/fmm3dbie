      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:), srccoefs(:,:)
      real *8, allocatable :: wts(:)
      character *100 fname
      integer ipars(2)

      integer, allocatable :: norders(:), ixyzs(:), iptype(:)

      real *8 xyz_out(3), xyz_in(3)
      real *8, allocatable :: sigma(:), rhs(:)
      real *8, allocatable :: errs(:)
      real *8 eps_gmres
      real *8 dpars(2)
      real *8 c0(3)
      integer nuv(2)
      complex * 16 zpars(3)
      integer numit,niter

      integer ipatch_id
      real *8 uvs_targ(2), dtmp
      integer iptype0

      real *8 pot,potex
      character *100, igeom

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4

      norder = 4 

      dpars(1) = -1.0d0
      dpars(2) = 1.0d0 

c
c  patch type, iptype0 = 1, for triangles with RV nodes
c                      = 11, for quadrangles with GL nodes
c                      = 12, for quadrangles with Cheb nodes      
c      
      iptype0 = 1

      igeom = 'stellarator'
      if (trim(igeom).eq.'sphere') then
        a = 1
        na = 4
        c0(1) = 0
        c0(2) = 0
        c0(3) = 0
        call get_sphere_npat_mem(a, na, c0, norder, iptype0, 
     1    npatches, npts)

        allocate(srcvals(12,npts), srccoefs(9,npts))
        allocate(ixyzs(npatches+1), iptype(npatches))
        allocate(norders(npatches))

        call get_sphere_npat(a, na, c0, norder, iptype0,
     1    npatches, npts, norders, ixyzs, iptype, srccoefs, srcvals)

        xyz_out(1) = 3.17d0
        xyz_out(2) = -0.03d0
        xyz_out(3) = 3.15d0

        xyz_in(1) = 0.17d0
        xyz_in(2) = 0.23d0
        xyz_in(3) = -0.11d0
      elseif (trim(igeom).eq.'stellarator') then

        nuv(1) = 10
        nuv(2) = 30
        
        call get_stellarator_npat_mem(nuv, norder, iptype0, 
     1    npatches, npts)

        allocate(srcvals(12,npts), srccoefs(9,npts))
        allocate(ixyzs(npatches+1), iptype(npatches))
        allocate(norders(npatches))

        call get_stellarator_npat(nuv, norder, iptype0,
     1    npatches, npts, norders, ixyzs, iptype, srccoefs, srcvals)


        xyz_in(1) = -4.501d0
        xyz_in(2) = 1.7d-3
        xyz_in(3) = 0.00001d0

        xyz_out(1) = -3.5d0
        xyz_out(2) = 3.1d0
        xyz_out(3) = 20.1d0
      endif



      allocate(wts(npts))
      call get_qwts(npatches, norders, ixyzs, iptype, npts, srcvals, 
     1   wts)


      allocate(sigma(npts),rhs(npts))

      do i=1,npts
        call l3d_slp(xyz_out, 3, srcvals(1,i), 0, dpars, 1, zpars, 0,
     1     ipars, rhs(i))
        sigma(i) = 0 
      enddo

      numit = 200
      ifinout = 0
      niter = 0
      allocate(errs(numit+1))

      eps = 0.51d-4

      eps_gmres = 0.5d-6

      call lap_comb_dir_solver(npatches, norders, ixyzs, iptype, npts,
     1  srccoefs, srcvals, eps, dpars, numit, ifinout, rhs, eps_gmres,
     2  niter, errs, rres, sigma)

      call prinf('niter=*',niter,1)
      call prin2('rres=*',rres,1)
      call prin2('errs=*',errs,niter)

c
c       test solution at interior point
c
      call l3d_slp(xyz_out, 3, xyz_in, 0, dpars, 1, zpars, 0,
     1   ipars,potex)
      pot = 0
      do i=1,npts
        call l3d_comb(srcvals(1,i), 3, xyz_in, 0, dpars, 3, zpars,
     1    0, ipars, dtmp)
        pot = pot + sigma(i)*wts(i)*dtmp
      enddo

      call prin2('potex=*',potex,1)
      call prin2('pot=*',pot,1)
      erra = abs(pot-potex)/abs(potex)
      call prin2('relative error after direct computation=*',erra,1)

      ndtarg = 3
      ntarg = 1
      ipatch_id = -1
      uvs_targ(1) = 0
      uvs_targ(2) = 0
      call lap_comb_dir_eval(npatches, norders, ixyzs, iptype,
     1  npts, srccoefs, srcvals, ndtarg, ntarg, xyz_in, ipatch_id,
     2  uvs_targ, eps, dpars, sigma, pot)

      call prin2('potex=*',potex,1)
      call prin2('pot=*',pot,1)
      erra = abs(pot-potex)/abs(potex)
      call prin2('relative error after fmm call=*',erra,1)



      stop
      end


