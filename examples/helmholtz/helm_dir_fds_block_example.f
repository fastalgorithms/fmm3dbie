      implicit real *8 (a-h,o-z) 
      implicit integer *8 (i-n)
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:)
      character *100 fname
      integer *8 ipars(2)

      real *8, allocatable :: targs(:,:),uvs_targ(:,:)
      integer *8, allocatable :: ipatch_id(:)

      integer *8, allocatable :: norders(:),ixyzs(:),iptype(:)

      integer *8, allocatable :: ifds(:)
      real *8, allocatable :: rfds(:)
      complex *16, allocatable :: zfds(:)

      real *8 xyz_out(3),xyz_in(3)
      complex *16, allocatable :: sigma(:),rhs(:),sigma2(:)
      complex *16 zid
      real *8, allocatable :: errs(:)
      complex * 16 zpars(3)
      integer *8, allocatable :: irand(:),isort(:),isum(:)
      integer *8, allocatable :: jrand(:),jsort(:)
      real *8, allocatable :: rrand(:)
      complex *16, allocatable :: xmat(:,:),xmattmp(:,:)
      integer *8, allocatable :: itarg(:),jsrc(:)
      complex *16, allocatable :: zent(:)
      real *8 c0(3)
      integer *8 numit,niter
      integer *8 ifwrite
      integer *8, allocatable :: row_ind(:),col_ind(:)
      integer *8 int8_0,int8_1,int8_3

      complex *16 pot,potex,ztmp,ima
      
      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      int8_0 = 0
      int8_1 = 1
      int8_3 = 3
      done = 1
      pi = atan(done)*4

      norder = 3 
c
c  patch type, iptype0 = 1, for triangles with RV nodes
c                      = 11, for quadrangles with GL nodes
c                      = 12, for quadrangles with Cheb nodes      
c      
      iptype0 = 1

      a = 1
      na = 4
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

      xyz_out(1) = 3.17d0
      xyz_out(2) = -0.03d0
      xyz_out(3) = 3.15d0

      xyz_in(1) = 0.17d0
      xyz_in(2) = 0.23d0
      xyz_in(3) = -0.11d0


      zk = 1.11d0+ima*0.0d0
      zpars(1) = zk 
      zpars(2) = -3.0d0
      zpars(3) = 1.0d0

      allocate(wts(npts))
      allocate(targs(3,npts),ipatch_id(npts),uvs_targ(2,npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)


      ndtarg = 3
     
      do i=1,npts
        targs(1,i) = srcvals(1,i)
        targs(2,i) = srcvals(2,i)
        targs(3,i) = srcvals(3,i)
      enddo

      do i=1,npts
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo

      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts, 
     1         ipatch_id,uvs_targ)


      allocate(sigma(npts),rhs(npts),sigma2(npts))

      do i=1,npts
        rhs(i) = 0
      enddo


      do i=1,npts
        call h3d_slp(xyz_out,int8_3,srcvals(1,i),int8_0,dpars,
     1     int8_1,zpars,int8_0,ipars,rhs(i))
      enddo


      eps = 0.51d-3

      call helm_comb_dir_fds_block_mem(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,eps,zpars,ifwrite,nifds,nrfds,
     1  nzfds)

      call prinf('nifds=*',nifds,1)
      call prinf('nrfds=*',nrfds,1)
      call prinf('nzfds=*',nzfds,1)

      allocate(ifds(nifds),rfds(nrfds),zfds(nzfds))

      
      call helm_comb_dir_fds_block_init(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,eps,zpars,nifds,ifds,nzfds,zfds)
      
      allocate(irand(npts),isort(npts))
      allocate(jrand(npts),jsort(npts))

      do i=1,npts
        irand(i) = hkrand(0)*npts
        jrand(i) = hkrand(0)*npts
      enddo


      call sorti(npts,irand,isort)
      call sorti(npts,jrand,jsort)

      allocate(xmat(npts,npts))
      allocate(row_ind(npts),col_ind(npts))

      do i=1,npts
        do j=1,npts
          xmat(j,i) = 0
        enddo
      enddo


      nbat = 30
      npbat = (npts+0.0d0)/(nbat+0.0d0)
      do ibat=1,nbat
        istart = (ibat-1)*npbat+1
        iend = min(npts,ibat*npbat)
        nrow = iend-istart+1
        do i=1,nrow
          row_ind(i) = isort(istart+i-1)
        enddo

        do jbat=1,nbat
          jstart = (jbat-1)*npbat+1
          jend = min(npts,jbat*npbat)
          ncol = jend-jstart+1
          do j=1,ncol
            col_ind(j) = jsort(jstart+j-1)
          enddo

          allocate(xmattmp(nrow,ncol))

          call helm_comb_dir_fds_block_matgen(npatches,norders,
     1     ixyzs,iptype,npts,srccoefs,srcvals,wts,eps,zpars,nifds,
     2     ifds,nzfds,zfds,nrow,row_ind,ncol,col_ind,ifwrite,xmattmp)
          

          
          do j=1,ncol
            do i=1,nrow
              xmat(row_ind(i),col_ind(j))=xmattmp(i,j)
            enddo
          enddo
          deallocate(xmattmp)
        enddo
      enddo

      do i=1,npts
        sigma(i) = 0
        sigma2(i) = 0
        do j=1,npts
          sigma(i) = sigma(i) + xmat(i,j)*rhs(j)
        enddo
      enddo



      call helm_comb_dir_eval(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,ipatch_id,
     2  uvs_targ,eps,zpars,rhs,sigma2)
      
      err = 0
      ra = 0
      do i=1,npts
        err = err + abs(sigma(i)-sigma2(i))**2
        ra = ra + abs(sigma2(i))**2
      enddo
      err = sqrt(err/ra)
      call prin2('error in matvec=*',err,1)


      ifinout = 0 
      zid = -(-1)**(ifinout)*zpars(3)/2.0d0
      do i=1,npts
        xmat(i,i) = xmat(i,i) + zid 
      enddo


      call zgausselim(npts,xmat,rhs,info,sigma,dcond)


c
c       test solution at interior point
c
      call h3d_slp(xyz_out,int8_3,xyz_in,int8_0,dpars,int8_1,
     1             zpars,int8_0,ipars,potex)
      pot = 0
      do i=1,npts
        call h3d_comb(srcvals(1,i),int8_3,xyz_in,int8_0,dpars,
     1     int8_3,zpars,int8_0,ipars,ztmp)
        pot = pot + sigma(i)*wts(i)*ztmp
      enddo

      erra = abs(pot-potex)/abs(potex)
      call prin2('relative error in solve=*',erra,1)

      stop
      end


