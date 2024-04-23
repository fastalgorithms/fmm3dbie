      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:)
      character *100 fname
      integer ipars(2)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)

      real *8 xyz_out(3),xyz_in(3)
      complex *16, allocatable :: sigma(:),rhs(:)
      complex *16 zk
      real *8, allocatable :: errs(:)
      real *8 eps_gmres
      complex * 16 zpars(3)
      integer numit,niter

      integer ipatch_id
      real *8 uvs_targ(2)

      logical isout0,isout1

      complex *16 pot,potex,ztmp,ima
      real *8 c0(3)
      integer nuv(2)

      character *100, igeom

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4

      norder = 4 
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



      zk = 1.11d0+ima*0.0d0
      zpars(1) = zk 
      zpars(2) = ima*zk
      zpars(3) = 1.0d0 

      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)

      isout0 = .false.
      isout1 = .false.
      call test_exterior_pt(npatches,norder,npts,srcvals,srccoefs,
     1  wts,xyz_in,isout0)
      
      call test_exterior_pt(npatches,norder,npts,srcvals,srccoefs,
     1   wts,xyz_out,isout1)

       print *, isout0,isout1


      allocate(sigma(npts),rhs(npts))

      do i=1,npts
        call h3d_slp(xyz_out,3,srcvals(1,i),0,dpars,1,zpars,0,
     1     ipars,rhs(i))
        sigma(i) = 0 
      enddo



      numit = 200
      ifinout = 0
      niter = 0
      allocate(errs(numit+1))

      eps = 0.51d-4

      eps_gmres = 0.5d-6

      call helm_comb_dir_solver(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,eps,zpars,numit,ifinout,rhs,eps_gmres,
     2  niter,errs,rres,sigma)

      call prinf('niter=*',niter,1)
      call prin2('rres=*',rres,1)
      call prin2('errs=*',errs,niter)

 2000 continue

c
c       test solution at interior point
c
      call h3d_slp(xyz_out,3,xyz_in,0,dpars,1,zpars,0,ipars,potex)
      pot = 0
      do i=1,npts
        call h3d_comb(srcvals(1,i),3,xyz_in,0,dpars,3,zpars,0,ipars,
     1     ztmp)
        pot = pot + sigma(i)*wts(i)*ztmp
      enddo

      call prin2('potex=*',potex,2)
      call prin2('pot=*',pot,2)
      erra = abs(pot-potex)/abs(potex)
      call prin2('relative error=*',erra,1)

      ndtarg = 3
      ntarg = 1
      ipatch_id = -1
      uvs_targ(1) = 0
      uvs_targ(2) = 0
      call helm_comb_dir_eval(npatches,norders,ixyzs,iptype,
     1  npts,srccoefs,srcvals,ndtarg,ntarg,xyz_in,ipatch_id,
     2  uvs_targ,eps,zpars,sigma,pot)

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

   




