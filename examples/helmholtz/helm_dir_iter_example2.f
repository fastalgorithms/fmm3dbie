      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:)
      integer ipars(2)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)

      real *8 xyz_out(3),xyz_in(3)
      complex *16, allocatable :: sigma(:),rhs(:)
      complex *16, allocatable :: sigma2(:),rhs2(:)
      real *8, allocatable :: errs(:)
      real *8 thet,phi,eps_gmres
      complex * 16 zpars(3)
      integer numit,niter
      character *100 title,fname

      integer ipatch_id
      real *8 uvs_targ(2)

      logical isout0,isout1

      complex *16 pot,potex,ztmp,ima

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4


      zk = 4.4d0+ima*0.0d0
      zpars(1) = zk 
      zpars(2) = -ima*zk
      zpars(3) = 2.0d0

      
      xyz_in(1) = 0.11d0
      xyz_in(2) = 0.0d-5
      xyz_in(3) = 0.37d0

      xyz_out(1) = -3.5d0
      xyz_out(2) = 3.1d0
      xyz_out(3) = 20.1d0

      fname = '../../geometries/sphere_192_o03.go3'
      
      call open_gov3_geometry_mem(fname,npatches,npts)

      call prinf('npatches=*',npatches,1)
      call prinf('npts=*',npts,1)

      allocate(srcvals(12,npts),srccoefs(9,npts))
      allocate(ixyzs(npatches+1),iptype(npatches),norders(npatches))
      allocate(wts(npts))

      call open_gov3_geometry(fname,npatches,norders,ixyzs,
     1   iptype,npts,srcvals,srccoefs,wts)
      

      allocate(sigma(npts),rhs(npts))
      ifinout = 1

      do i=1,npts
        if(ifinout.eq.0) 
     1     call h3d_slp(xyz_out,3,srcvals(1,i),0,dpars,1,zpars,0,
     2       ipars,rhs(i))
        if(ifinout.eq.1) 
     1     call h3d_slp(xyz_in,3,srcvals(1,i),0,dpars,1,zpars,0,
     1     ipars,rhs(i))
        rhs(i) = rhs(i)
        sigma(i) = 0
      enddo


      numit = 200
      niter = 0
      allocate(errs(numit+1))

      eps = 0.51d-6
      eps_gmres = eps
c
c  estimate memory required
c

      call cpu_time(t1)
C$      t1 = omp_get_wtime()      
      call helm_comb_dir_solver_memest(npatches,norders,ixyzs,iptype,
     1  npts,srccoefs,srcvals,eps,zpars,numit,rmem)
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      
      call prin2('memory required in GB=*',rmem,1)
      call prin2('time taken in memory estimation code=*',t2-t1,1)


      call helm_comb_dir_solver(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,eps,zpars,numit,ifinout,rhs,eps_gmres,
     2  niter,errs,rres,sigma)

      call prinf('niter=*',niter,1)
      call prin2('rres=*',rres,1)
      call prin2('errs=*',errs,niter)


c
c       test solution at interior point
c
      call h3d_slp(xyz_out,3,xyz_in,0,dpars,1,zpars,0,ipars,potex)
      pot = 0
      do i=1,npts
        if(ifinout.eq.0) 
     1     call h3d_comb(srcvals(1,i),3,xyz_in,0,dpars,3,zpars,1,ipars,
     2      ztmp)
        if(ifinout.eq.1) 
     1     call h3d_comb(srcvals(1,i),3,xyz_out,0,dpars,3,zpars,1,ipars,
     1     ztmp)
        pot = pot + sigma(i)*wts(i)*ztmp
      enddo

      call prin2('potex=*',potex,2)
      call prin2('pot=*',pot,2)
      erra = abs(pot-potex)/abs(potex)
      call prin2('relative error=*',erra,1)

      stop
      end





