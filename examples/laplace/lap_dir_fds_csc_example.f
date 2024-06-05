      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:)
      character *100 fname
      integer ipars(2)

      real *8, allocatable :: targs(:,:),uvs_targ(:,:)
      integer, allocatable :: ipatch_id(:)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)

      integer, allocatable :: ifds(:)
      real *8, allocatable :: rfds(:)
      complex *16, allocatable :: zfds(:)

      real *8 xyz_out(3),xyz_in(3)
      real *8, allocatable :: sigma(:),rhs(:),sigma2(:)
      real *8 did
      real *8, allocatable :: errs(:)
      real *8 dpars(2)
      integer, allocatable :: irand(:),isort(:),isum(:)
      real *8, allocatable :: rrand(:)
      real *8, allocatable :: xmat(:,:)
      integer, allocatable :: itarg(:),jsrc(:)
      integer, allocatable :: col_ptr(:),row_ind(:)
      real *8, allocatable :: zent(:)
      integer numit,niter
      real *8 c0(3)

      real *8 pot,potex
      complex *16 ztmp,ima

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

      xyz_out(1) = 3.17d0
      xyz_out(2) = -0.03d0
      xyz_out(3) = 3.15d0

      xyz_in(1) = 0.17d0
      xyz_in(2) = 0.23d0
      xyz_in(3) = -0.11d0


      dpars(1) = -1.1d0
      dpars(2) = 1.2d0


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
        call l3d_slp(xyz_out,3,srcvals(1,i),0,dpars,0,zpars,0,ipars,
     1     rhs(i))
      enddo

      

      eps = 0.51d-6

      print *, dpars(1),dpars(2)

      call lap_comb_dir_fds_csc_mem(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,eps,dpars,nifds,nrfds,nzfds)

      call prinf('nifds=*',nifds,1)
      call prinf('nrfds=*',nrfds,1)
      call prinf('nzfds=*',nzfds,1)

      allocate(ifds(nifds),rfds(nrfds),zfds(nzfds))

      
      call lap_comb_dir_fds_csc_init(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,eps,dpars,nifds,ifds,nrfds,rfds,nzfds,zfds)




c
c       create a permutation of indices from 1 to npts**2
c       and divide them into nbat random sized partitions
c
      nn = npts**2
      allocate(irand(nn),isort(nn))
      do i=1,nn
        irand(i) = hkrand(0)*nn
      enddo

      call sorti(nn,irand,isort)

      nbat = 50
      allocate(rrand(nbat))
      ra = 0
      do i=1,nbat
        rrand(i) = hkrand(0)
        ra = ra + rrand(i)
      enddo

   
      do i=1,nbat
         rrand(i) = rrand(i)/ra
         irand(i) = rrand(i)*nn
      enddo

      allocate(isum(nbat+1))

      isum(1) = 1
      call cumsum(nbat,irand,isum(2))
      

      if(isum(nbat+1).lt.nn+1) isum(nbat+1) = nn+1


      allocate(col_ptr(npts+1))
      allocate(xmat(npts,npts))

      do i=1,npts
        do j=1,npts
          xmat(i,j) = 0
        enddo
      enddo

      do ibat=1,nbat
        nent = isum(ibat+1)-isum(ibat)

        allocate(zent(nent),row_ind(nent),itarg(nent),jsrc(nent))
        do i=1,nent
          iind = isort(isum(ibat)+i-1)

          idiv = iind/npts
          idiv = idiv+1
          irem = iind - (idiv-1)*npts
          if(irem.eq.0) then
            idiv = idiv - 1
            irem = npts
          endif

          jsrc(i) = irem
          itarg(i) = idiv

        enddo


        call conv_to_csc(nent,npts,itarg,jsrc,col_ptr,row_ind)


        

        
        
        
        call lap_comb_dir_fds_csc_matgen(npatches,norders,ixyzs,iptype,
     1     npts,srccoefs,srcvals,eps,dpars,nifds,ifds,nrfds,rfds,nzfds,
     2     zfds,nent,col_ptr,row_ind,zent)

        
       
cc        call prinf('col_ptr=*',col_ptr,npts+1)

        do i=1,npts
          do j=col_ptr(i),col_ptr(i+1)-1
            xmat(row_ind(j),i) = zent(j)
          enddo
        enddo

        deallocate(zent,row_ind,itarg,jsrc)
      enddo

      do i=1,npts
        sigma(i) = 0
        sigma2(i) = 0
        do j=1,npts
          sigma(i) = sigma(i) + xmat(i,j)*rhs(j)
        enddo
      enddo


      call lap_comb_dir_eval(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,ipatch_id,
     2  uvs_targ,eps,dpars,rhs,sigma2)

      err = 0
      ra = 0
      do i=1,npts
        err = err + abs(sigma(i)-sigma2(i))**2
        ra = ra + abs(sigma2(i))**2
      enddo
      err = sqrt(err/ra)
      call prin2('error in matvec=*',err,1)

      ifinout = 0
      did = -(-1)**(ifinout)*dpars(2)/2
      do i=1,npts
        xmat(i,i) = xmat(i,i) + did 
      enddo


      call dgausselim(npts,xmat,rhs,info,sigma,dcond)


c
c       test solution at interior point
c
      call l3d_slp(xyz_out,3,xyz_in,0,dpars,0,zpars,0,ipars,potex)
      pot = 0
      do i=1,npts
        call l3d_comb(srcvals(1,i),3,xyz_in,2,dpars,0,zpars,0,ipars,
     1     ztmp)
        pot = pot + sigma(i)*wts(i)*ztmp
      enddo

      erra = abs(pot-potex)/abs(potex)
      call prin2('relative error in solve=*',erra,1)

      stop
      end


