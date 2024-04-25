      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:), srccoefs(:,:)
      integer, allocatable :: norders(:), ixyzs(:), iptype(:)
      real *8, allocatable :: wts(:), rsigma(:)
      integer ipars(2)


      real *8 xyz_out(3), xyz_in(3)

      complex *16, allocatable :: einc(:,:), hinc(:,:)
      complex *16, allocatable :: zjvec(:,:), rho(:)
      complex *16 zedips(3), zhdips(3) 
      complex *16 e_ex(3), h_ex(3)
      complex *16 e_comp(3), h_comp(3)


      real *8, allocatable :: errs(:)
      real *8 thet,phi
      complex * 16 zpars(3)
      integer numit,niter
      character *200 title, fname, fname1, fname2

      integer ipatch_id
      real *8 uvs_targ(2)

      logical isout0, isout1

      complex *16 pot, potex, ztmp, ima, zk
      complex *16 alpha_rhs
      real *8 ptinfo_out(12)
      real *8 c0(3)

      integer count1

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4


      norder = 4 
!
!  patch type, iptype0 = 1, for triangles with RV nodes
!                      = 11, for quadrangles with GL nodes
!                      = 12, for quadrangles with Cheb nodes      
!      
      iptype0 = 1

      a = 1
      na = 4
      c0(1) = 0
      c0(2) = 0
      c0(3) = 0

      call get_sphere_npat_mem(a, na, c0, norder, iptype0, &
       npatches, npts)

      allocate(srcvals(12,npts), srccoefs(9,npts))
      allocate(ixyzs(npatches+1), iptype(npatches))
      allocate(norders(npatches))

      call get_sphere_npat(a, na, c0, norder, iptype0, &
        npatches, npts, norders, ixyzs, iptype, srccoefs, srcvals)

      xyz_out(1) = 3.17d0
      xyz_out(2) = -0.03d0
      xyz_out(3) = 3.15d0

      xyz_in(1) = 0.17d0
      xyz_in(2) = 0.23d0
      xyz_in(3) = -0.11d0


      zk = 1.1d0 
      zpars(1) = zk 
      zpars(2) = 1.0d0
      zpars(3) = 0.0d0

      allocate(wts(npts))


      allocate(einc(3,npts), hinc(3,npts))
      allocate(zjvec(3,npts), rho(npts))

      zedips(1) = 3.1d0
      zedips(2) = 1.2d0
      zedips(3) = -0.1d0 + ima*3.3d0

      zhdips(1) = 1
      zhdips(2) = ima
      zhdips(3) = 1.1d0
      
  
!
!   Get fields due to magnetic and electric dipoles
!     
      call em_eval_dipoles(zpars(1), 1, xyz_out, zedips, zhdips, &
        12, npts, srcvals, einc, hinc)
      call prin2('einc=*', einc, 24)
      call prin2('hinc=*', hinc, 24)
            

      numit = 400
      niter = 0
      allocate(errs(numit+1))

      eps = 1d-6
      eps_gmres = 1d-10

      call cpu_time(t1)
!$      t1 = omp_get_wtime()          

      call em_nrccie_pec_solver(npatches, norders, ixyzs, iptype, &
        npts, srccoefs, srcvals, eps, zpars, numit, einc, hinc, &
        eps_gmres, niter, errs, rres, zjvec, rho)

      call prinf('niter=*',niter,1)
      call prin2('rres=*',rres,1)
      call prin2('errs=*',errs,niter)

      call cpu_time(t2)
!$       t2 = omp_get_wtime()
      call prin2('analytic solve time=*',t2-t1,1)

      
!
!       test solution at exterior point
!

      ntarg = 1
      ndtarg = 3

      call em_eval_dipoles(zpars(1), 1, xyz_out, zedips, zhdips, &
        ndtarg, ntarg, xyz_in, e_ex, h_ex)

      ife = 1
      ifh = 1
      ipatch_id = -1
      uvs_targ(1) = 0
      uvs_targ(2) = 0
      e_comp(1:3) = 0
      h_comp(1:3) = 0

      call em_nrccie_pec_eval(npatches, norders, ixyzs, &
            iptype, npts, srccoefs, srcvals, ndtarg, ntarg, xyz_in, &
            ipatch_id, uvs_targ, eps, zpars, zjvec, rho, ife, e_comp, &
            ifh, h_comp)

      erra = 0
      ra = 0
      erra = erra + abs(e_comp(1) + e_ex(1))**2
      erra = erra + abs(e_comp(2) + e_ex(2))**2
      erra = erra + abs(e_comp(3) + e_ex(3))**2

      erra = erra + abs(h_comp(1) + h_ex(1))**2
      erra = erra + abs(h_comp(2) + h_ex(2))**2
      erra = erra + abs(h_comp(3) + h_ex(3))**2

      ra = ra + abs(e_ex(1))*2 + abs(e_ex(2))**2 + abs(e_ex(3))**2
      ra = ra + abs(h_ex(1))*2 + abs(h_ex(2))**2 + abs(h_ex(3))**2

      erra = sqrt(erra/ra)
      call prin2('e_comp=*',e_comp,6)
      call prin2('e_ex=*',e_ex,6)

      call prin2('h_comp=*',h_comp,6)
      call prin2('h_ex=*',h_ex,6)
      print *, ""
      print *, ""
      print *, "error in computed fields=",erra

      stop
      end





