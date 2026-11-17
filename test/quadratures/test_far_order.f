      subroutine test_far_order(nsuccess)
c
c  Test get_far_order_guru on patches up to and above the largest
c  far-field rule (order 16)
c
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8 iptypes(3),norders(5),nfars(10)
      real *8 cm(3),targs(3,30)
      real *8, allocatable :: uvs(:,:),umat(:,:),vmat(:,:),wts(:)
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      complex *16 zk

      data iptypes /1,11,12/
      data norders /16,17,20,30,40/
      data nfars /1,2,3,4,6,8,10,12,14,16/

      eps = 1.0d-8
      zk = 1.0d0
      ndtarg = 3
      ntarg = 30
      rad = 3.0d0

      nsuccess = 1
      do it=1,3
        iptype = iptypes(it)
        do io=1,5
          norder = norders(io)
          if(iptype.eq.1.and.norder.gt.20) cycle

          call get_npols(iptype,norder,npols)
          allocate(uvs(2,npols),umat(npols,npols),vmat(npols,npols))
          allocate(wts(npols),srcvals(9,npols),srccoefs(9,npols))
          call get_disc_exps(norder,npols,iptype,uvs,umat,vmat,wts)
c
c  curved patch z = 0.1 u^2, and its expansion coefficients
c
          do i=1,npols
            u = uvs(1,i)
            v = uvs(2,i)
            srcvals(1,i) = u
            srcvals(2,i) = v
            srcvals(3,i) = 0.1d0*u**2
            srcvals(4,i) = 1
            srcvals(5,i) = 0
            srcvals(6,i) = 0.2d0*u
            srcvals(7,i) = 0
            srcvals(8,i) = 1
            srcvals(9,i) = 0
          enddo
          srccoefs = matmul(srcvals,transpose(umat))

          cm(1) = sum(uvs(1,:)*wts)/sum(wts)
          cm(2) = sum(uvs(2,:)*wts)/sum(wts)
          cm(3) = 0
          do j=1,ntarg
            r = 0.5d0 + 0.1d0*j
            targs(1,j) = cm(1) + r*cos(1.0d0*j)
            targs(2,j) = cm(2) + r*sin(1.0d0*j)
            targs(3,j) = cm(3) + 0.3d0*r
          enddo

          do ikerorder=-1,1
            call get_far_order_guru(eps,norder,npols,iptype,cm,rad,
     1        srccoefs,ikerorder,zk,ndtarg,ntarg,targs,nfar,npolsf)
            call get_npols(iptype,nfar,npolsf0)
            if(.not.any(nfars.eq.nfar).or.npolsf.ne.npolsf0) then
              call prinf('bad far order, iptype=*',iptype,1)
              call prinf('norder=*',norder,1)
              call prinf('nfar=*',nfar,1)
              nsuccess = 0
            endif
          enddo

          deallocate(uvs,umat,vmat,wts,srcvals,srccoefs)
        enddo
      enddo

c
c  repeat on axisymmetric spheres, whose polar caps go through
c  the x{tri,quad}_axissym_*circ_eval routines
c
      call test_far_order_axissym(nsuccess)

      return
      end
c
c
c
c
      subroutine test_far_order_axissym(nsuccess)
c
c  Build the unit sphere as an axisymmetric surface, with the
c  generating curve given both as a function handle
c  (get_axissym_fcurve_npat) and as legendre chunks
c  (get_axissym_npat), at orders above the old fixed work
c  array sizes. On each surface check the area, surf_fun_error,
c  fun_surf_interp_coefs, and get_far_order_guru.
c
c  nsuccess is set to 0 on failure, and left unchanged otherwise
c
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8 iptypes(3),nords(3),nrts(2,3),ixys2d(4),iptype2d(3)
      integer *8 nfars(10)
      real *8 tchse(4),pars(2),targs(3,30),errm
      real *8 xs(20),umat2d(20,20),vmat2d(20,20),ws2d(20)
      real *8 srcvals2d(8,60),srccoefs2d(6,60)
      integer *8, allocatable :: norders(:),ixyzs(:),iptype(:)
      integer *8, allocatable :: ipatches(:)
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:),wts(:)
      real *8, allocatable :: cms(:,:),rads(:),errp(:)
      real *8, allocatable :: uvs(:,:),uvs_targ(:,:),uover(:,:)
      complex *16 zk
      external fcurve_sphere_tfo

      data iptypes /1,11,12/
      data nords /20,25,25/
      data nfars /1,2,3,4,6,8,10,12,14,16/

      done = 1
      pi = atan(done)*4
      eps = 1.0d-8
      zk = 1.0d0
      ndtarg = 3
      ntarg = 30

      nch2d = 3
      do i=1,nch2d+1
        tchse(i) = (i-1.0d0)/nch2d
      enddo
      do i=1,nch2d
        nrts(1,i) = 2
        nrts(2,i) = 4
      enddo
      rmid = 0.5d0
      nmid = 2
      iort = 1
      np = 2
      pars(1) = 0
      pars(2) = 0
c
c  legendre chunks for r = sin(pi t), z = cos(pi t), t in [0,1],
c  derivatives taken with respect to the local chunk variable
c
      k = 20
      itype = 2
      call legeexps(itype,k,xs,umat2d,vmat2d,ws2d)
      npts2d = nch2d*k
      do ich=1,nch2d
        ixys2d(ich) = (ich-1)*k + 1
        iptype2d(ich) = 1
        h = tchse(ich+1) - tchse(ich)
        do j=1,k
          ipt = ixys2d(ich) + j - 1
          t = tchse(ich) + (xs(j)+1)/2*h
          call fcurve_sphere_tfo(t,np,pars,r,z,drdt,dzdt,
     1       d2rdt2,d2zdt2)
          srcvals2d(1,ipt) = r
          srcvals2d(2,ipt) = z
          srcvals2d(3,ipt) = drdt*h/2
          srcvals2d(4,ipt) = dzdt*h/2
          srcvals2d(5,ipt) = d2rdt2*h*h/4
          srcvals2d(6,ipt) = d2zdt2*h*h/4
          ds = sqrt(srcvals2d(3,ipt)**2 + srcvals2d(4,ipt)**2)
          srcvals2d(7,ipt) = srcvals2d(4,ipt)/ds
          srcvals2d(8,ipt) = -srcvals2d(3,ipt)/ds
        enddo
        do j=1,k
          ipt = ixys2d(ich) + j - 1
          do l=1,6
            srccoefs2d(l,ipt) = 0
          enddo
          do jj=1,k
            jpt = ixys2d(ich) + jj - 1
            do l=1,6
              srccoefs2d(l,ipt) = srccoefs2d(l,ipt) +
     1           umat2d(j,jj)*srcvals2d(l,jpt)
            enddo
          enddo
        enddo
      enddo
      ixys2d(nch2d+1) = npts2d + 1

      do it=1,3
        iptype0 = iptypes(it)
        norder = nords(it)
        do imode=1,2
          npatches = 0
          npts = 0
          if(imode.eq.1) then
            call get_axissym_fcurve_npat_mem(nch2d,tchse,
     1         fcurve_sphere_tfo,np,pars,nrts,rmid,nmid,iort,
     2         norder,iptype0,npatches,npts)
          else
            call get_axissym_npat_mem(nch2d,npts2d,ixys2d,iptype2d,
     1         srcvals2d,srccoefs2d,nrts,rmid,nmid,iort,
     2         norder,iptype0,npatches,npts)
          endif

          allocate(srcvals(12,npts),srccoefs(9,npts),wts(npts))
          allocate(norders(npatches),ixyzs(npatches+1))
          allocate(iptype(npatches))

          if(imode.eq.1) then
            call get_axissym_fcurve_npat(nch2d,tchse,
     1         fcurve_sphere_tfo,np,pars,nrts,rmid,nmid,iort,
     2         norder,iptype0,npatches,npts,norders,ixyzs,iptype,
     3         srccoefs,srcvals)
          else
            call get_axissym_npat(nch2d,npts2d,ixys2d,iptype2d,
     1         srcvals2d,srccoefs2d,nrts,rmid,nmid,iort,
     2         norder,iptype0,npatches,npts,norders,ixyzs,iptype,
     3         srccoefs,srcvals)
          endif
c
c  area of the unit sphere
c
          call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,
     1       wts)
          ra = sum(wts)
          erra = abs(ra - 4*pi)/(4*pi)
c
c  tail error of xyz, should be resolved
c
          allocate(errp(npatches))
          nd = 3
          call surf_fun_error(nd,npatches,norders,ixyzs,iptype,npts,
     1       srcvals(1:3,1:npts),wts,errp,errm)
c
c  interpolate xyz back to the discretization nodes
c
          call get_npols(iptype0,norder,npols)
          allocate(uvs(2,npols),ipatches(npts),uvs_targ(2,npts))
          allocate(uover(9,npts))
          call get_disc_nodes(norder,npols,iptype0,uvs)
          do ip=1,npatches
            do j=1,npols
              ipt = ixyzs(ip) + j - 1
              ipatches(ipt) = ip
              uvs_targ(1,ipt) = uvs(1,j)
              uvs_targ(2,ipt) = uvs(2,j)
            enddo
          enddo
          nd = 9
          call fun_surf_interp_coefs(nd,npatches,norders,ixyzs,
     1       iptype,npts,srccoefs,npts,ipatches,uvs_targ,uover)
          erri = 0
          do i=1,npts
            do l=1,3
              erri = max(erri,abs(uover(l,i)-srcvals(l,i)))
            enddo
          enddo

          if(erra.gt.1.0d-9.or.errm.gt.1.0d-6.or.
     1       erri.gt.1.0d-10) then
            call prinf('axissym failed, iptype=*',iptype0,1)
            call prinf('norder=*',norder,1)
            call prinf('imode=*',imode,1)
            call prin2('area err=*',erra,1)
            call prin2('surf_fun_error=*',errm,1)
            call prin2('interp err=*',erri,1)
            nsuccess = 0
          endif
c
c  far order on a polar cap patch and on a patch from the middle
c
          allocate(cms(3,npatches),rads(npatches))
          call get_centroid_rads(npatches,norders,ixyzs,iptype,npts,
     1       srccoefs,cms,rads)
          do iip=1,2
            ip = 1
            if(iip.eq.2) ip = npatches/2
            rad = 3*rads(ip)
            do j=1,ntarg
              r = (0.5d0 + 0.1d0*j)*rads(ip)
              targs(1,j) = cms(1,ip) + r*cos(1.0d0*j)
              targs(2,j) = cms(2,ip) + r*sin(1.0d0*j)
              targs(3,j) = cms(3,ip) + 0.3d0*r
            enddo
            istart = ixyzs(ip)
            ikerorder = 0
            call get_far_order_guru(eps,norder,npols,iptype0,
     1         cms(1,ip),rad,srccoefs(1,istart),ikerorder,zk,ndtarg,
     2         ntarg,targs,nfar,npolsf)
            call get_npols(iptype0,nfar,npolsf0)
            if(.not.any(nfars.eq.nfar).or.npolsf.ne.npolsf0) then
              call prinf('axissym bad far order, iptype=*',iptype0,1)
              call prinf('norder=*',norder,1)
              call prinf('nfar=*',nfar,1)
              nsuccess = 0
            endif
          enddo

          deallocate(srcvals,srccoefs,wts,norders,ixyzs,iptype)
          deallocate(errp,uvs,ipatches,uvs_targ,uover,cms,rads)
        enddo
      enddo

      return
      end
c
c
c
c
      subroutine fcurve_sphere_tfo(t,np,pars,r,z,drdt,dzdt,
     1   d2rdt2,d2zdt2)
c
c  generating curve of the unit sphere, t in [0,1]
c
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      real *8 pars(np)

      done = 1
      pi = atan(done)*4
      r = sin(pi*t)
      z = cos(pi*t)
      drdt = pi*cos(pi*t)
      dzdt = -pi*sin(pi*t)
      d2rdt2 = -pi*pi*sin(pi*t)
      d2zdt2 = -pi*pi*cos(pi*t)

      return
      end
