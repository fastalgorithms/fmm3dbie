      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:), srccoefs(:,:)
      real *8, allocatable :: wts(:)
      character *100 fname
      integer ipars(2)

      integer, allocatable :: norders(:), ixyzs(:), iptype(:)

      real *8 xyz_out(3), xyz_in(3), xyz_src(3), xyz_targ(3)
      complex *16, allocatable :: sigma(:), rhs(:), sigma1(:)
      real *8, allocatable :: errs(:)
      real *8 eps_gmres
      complex * 16 zpars(3)
      integer numit, niter

      integer ipatch_id
      real *8 uvs_targ(2)

      logical isout0, isout1

      complex *16 pot, potex, ztmp, ima
      real *8 c0(3)
      character *100 rep

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4

      zk = 1.1d0 + ima*0.01d0
      zpars(1) = zk
      zpars(2) = 1.0d0

      norder = 4

c
c  patch type, iptype0 = 1, for triangles with RV nodes
c                      = 11, for quadrangles with GL nodes
c                      = 12, for quadrangles with Cheb nodes      
c      
      iptype0 = 1

c
c  representation, rep = 'rpcomb', right preconditioned combined
c                                  field representation             
c                      = 's', single layer representation
      rep = 's'      

      a = 1
      na = 2
      c0(1) = 0
      c0(2) = 0
      c0(3) = 0
      call get_sphere_npat_mem(a, na, c0, norder, iptype0, 
     1  npatches, npts)

      allocate(srcvals(12,npts), srccoefs(9,npts))
      allocate(ixyzs(npatches+1), iptype(npatches))
      allocate(norders(npatches))

      call get_sphere_npat(a, na, c0, norder, iptype0,
     1  npatches, npts, norders, ixyzs, iptype, srccoefs, srcvals)

      ra = 1.05d0
      thet = 0.8d0*pi
      phi = 1.13d0*2*pi

      xyz_out(1) = ra*sin(thet)*cos(phi) 
      xyz_out(2) = ra*sin(thet)*sin(phi) 
      xyz_out(3) = ra*cos(thet) 


      xyz_in(1) = 0.17d0
      xyz_in(2) = 0.23d0
      xyz_in(3) = -0.11d0

      ifinout = 1
      if(ifinout.eq.0) then
        xyz_src(1) = xyz_out(1)
        xyz_src(2) = xyz_out(2)
        xyz_src(3) = xyz_out(3)

        xyz_targ(1) = xyz_in(1)
        xyz_targ(2) = xyz_in(2)
        xyz_targ(3) = xyz_in(3)
      endif

      if(ifinout.eq.1) then
        xyz_src(1) = xyz_in(1)
        xyz_src(2) = xyz_in(2)
        xyz_src(3) = xyz_in(3)

        xyz_targ(1) = xyz_out(1)
        xyz_targ(2) = xyz_out(2)
        xyz_targ(3) = xyz_out(3)
      endif

      allocate(wts(npts))
      call get_qwts(npatches, norders, ixyzs, iptype, npts, srcvals, 
     1   wts)

      allocate(sigma(npts),rhs(npts),sigma1(npts))

      do i=1,npts
        call h3d_sprime(xyz_src, 12, srcvals(1,i), 0, dpars, 1, 
     1    zpars, 0, ipars, rhs(i))
        sigma(i) = 0
        sigma1(i) = 0
      enddo


      numit = 200
      niter = 0
      allocate(errs(numit+1))

      eps = 0.51d-4

      eps_gmres = 0.5d-6

      if (trim(rep).eq.'rpcomb') then

        call helm_rpcomb_neu_solver_memest(npatches, norders, ixyzs,
     1    iptype, npts, srccoefs, srcvals, eps, zpars, numit, rmem)
      
        call helm_rpcomb_neu_solver(npatches, norders, ixyzs, iptype, 
     1    npts, srccoefs, srcvals, eps, zpars, numit, ifinout, rhs, 
     2    eps_gmres, niter, errs, rres, sigma, sigma1)

      elseif (trim(rep).eq.'s') then

        call helm_s_neu_solver_memest(npatches, norders, ixyzs,
     1    iptype, npts, srccoefs, srcvals, eps, zpars, numit, rmem)
      
        call helm_s_neu_solver(npatches, norders, ixyzs, iptype, 
     1    npts, srccoefs, srcvals, eps, zpars, numit, ifinout, rhs, 
     2    eps_gmres, niter, errs, rres, sigma)

      endif

      call prin2('estimated memory=*',rmem,1)
      call prinf('niter=*',niter,1)
      call prin2('rres=*',rres,1)
      call prin2('errs=*',errs,niter)


c
c       test solution at interior point
c
      call h3d_slp(xyz_targ, 3, xyz_src, 0, dpars, 1, zpars, 0, 
     1  ipars, potex)

      ndtarg = 3
      ntarg = 1
      ipatch_id = -1
      uvs_targ(1) = 0
      uvs_targ(2) = 0
      
      if (trim(rep).eq.'rpcomb') then

        call helm_rpcomb_eval(npatches, norders, ixyzs, iptype,
     1    npts, srccoefs, srcvals, ndtarg, ntarg, xyz_targ, ipatch_id,
     2    uvs_targ, eps, zpars, sigma, sigma1, pot)

      elseif (trim(rep).eq.'s') then

        call helm_s_eval(npatches, norders, ixyzs, iptype,
     1    npts, srccoefs, srcvals, ndtarg, ntarg, xyz_targ, ipatch_id,
     2    uvs_targ, eps, zpars, sigma, pot)

      endif

      call prin2('potex=*',potex,2)
      call prin2('pot=*',pot,2)
      erra = abs(pot-potex)/abs(potex)
      call prin2('relative error=*',erra,1)



      stop
      end




      subroutine test_exterior_pt(npatches,norder,npts,srcvals,
     1   srccoefs,wts,xyzout,isout)
c
c
c  this subroutine tests whether the pt xyzin, is
c  in the exterior of a surface, and also estimates the error
c  in representing e^{ir/2}/r and \grad e^{ir/2}/r \cdot n
c  centered at the interior point. Whether a point 
c  is in the interior or not is tested using Gauss' 
c  identity for the flux due to a point charge
c
c
c  input:
c    npatches - integer
c       number of patches
c    norder - integer
c       order of discretization
c    npts - integer
c       total number of discretization points on the surface
c    srccoefs - real *8 (9,npts)
c       koornwinder expansion coefficients of geometry info
c    xyzout -  real *8 (3)
c       point to be tested
c
c  output: 
c    isout - boolean
c      whether the target is in the interior or not
c

      implicit none
      integer npatches,norder,npts,npols
      real *8 srccoefs(9,npts),srcvals(12,npts),xyzout(3),wts(npts)
      real *8 tmp(3)
      real *8 dpars,done,pi
      real *8, allocatable :: rsurf(:),err_p(:,:) 
      integer ipars,norderhead,nd
      complex *16, allocatable :: sigma_coefs(:,:), sigma_vals(:,:)
      complex *16 zk,val

      integer ipatch,j,i
      real *8 ra,ds
      logical isout

      done = 1
      pi = atan(done)*4

      npols = (norder+1)*(norder+2)/2


      zk = 0

      ra = 0



      do ipatch=1,npatches
        do j=1,npols
          i = (ipatch-1)*npols + j
          call h3d_sprime(xyzout,12,srcvals(1,i),0,dpars,1,zk,0,ipars,
     1       val)

          call cross_prod3d(srcvals(4,i),srcvals(7,i),tmp)
          ds = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)
          ra = ra + real(val)*wts(i)
        enddo
      enddo

      if(abs(ra+4*pi).le.1.0d-1) isout = .false.
      if(abs(ra).le.1.0d-1) isout = .true.

      return
      end

   




