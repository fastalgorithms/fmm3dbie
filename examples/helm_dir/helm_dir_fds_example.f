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
      complex *16, allocatable :: sigma(:),rhs(:),sigma2(:)
      complex *16 zid
      real *8, allocatable :: errs(:)
      complex * 16 zpars(3)
      integer, allocatable :: irand(:),isort(:),isum(:)
      real *8, allocatable :: rrand(:)
      complex *16, allocatable :: xmat(:,:)
      integer, allocatable :: itarg(:),jsrc(:)
      integer, allocatable :: col_ptr(:),row_ind(:)
      complex *16, allocatable :: zent(:)
      integer numit,niter

      complex *16 pot,potex,ztmp,ima

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4


c
c       select geometry type
c       igeomtype = 1 => sphere
c       igeomtype = 2 => stellarator
c 
      igeomtype = 2
      if(igeomtype.eq.1) ipars(1) = 1
      if(igeomtype.eq.2) ipars(1) = 10

      if(igeomtype.eq.1) then
        npatches = 12*(4**ipars(1))
      endif
      if(igeomtype.eq.2) then
        ipars(2) = ipars(1)*3
        npatches = 2*ipars(1)*ipars(2)
      endif


      zk = 1.11d0+ima*0.0d0
      zpars(1) = zk 
      zpars(2) = -3.0d0
      zpars(3) = 0.0d0

      if(igeomtype.eq.1) then
        xyz_out(1) = 3.17d0
        xyz_out(2) = -0.03d0
        xyz_out(3) = 3.15d0

        xyz_in(1) = 0.17d0
        xyz_in(2) = 0.23d0
        xyz_in(3) = -0.11d0
      endif

      if(igeomtype.eq.2) then
        xyz_out(1) = -3.5d0
        xyz_out(2) = 3.1d0
        xyz_out(3) = 20.1d0
      endif

      norder = 4 
      npols = (norder+1)*(norder+2)/2

      npts = npatches*npols
      allocate(srcvals(12,npts),srccoefs(9,npts))
      ifplot = 0


      call setup_geom(igeomtype,norder,npatches,ipars, 
     1       srcvals,srccoefs,ifplot,fname)

      allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))

      do i=1,npatches
        norders(i) = norder
        ixyzs(i) = 1 +(i-1)*npols
        iptype(i) = 1
      enddo

      ixyzs(npatches+1) = 1+npols*npatches
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
        call h3d_slp(xyz_out,3,srcvals(1,i),0,dpars,1,zpars,0,ipars,
     1     rhs(i))
      enddo

      

      eps = 0.51d-6

      call helm_comb_dir_fds_mem(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,eps,zpars,nifds,nrfds,nzfds)

      call prinf('nifds=*',nifds,1)
      call prinf('nrfds=*',nrfds,1)
      call prinf('nzfds=*',nzfds,1)

      allocate(ifds(nifds),rfds(nrfds),zfds(nzfds))

      
      call helm_comb_dir_fds_init(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,eps,zpars,nifds,ifds,nrfds,rfds,nzfds,zfds)




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


        

        
        
        
        call helm_comb_dir_fds_matgen(npatches,norders,ixyzs,iptype,
     1     npts,srccoefs,srcvals,eps,zpars,nifds,ifds,nrfds,rfds,nzfds,
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


      call lpcomp_helm_comb_dir(npatches,norders,ixyzs,
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



      zid = -(-1)**(ifinout)*2*pi*zpars(3)
      do i=1,npts
        xmat(i,i) = xmat(i,i) + zid 
      enddo


      call zgausselim(npts,xmat,rhs,info,sigma,dcond)


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

      erra = abs(pot-potex)/abs(potex)
      call prin2('relative error in solve=*',erra,1)

      stop
      end




      subroutine setup_geom(igeomtype,norder,npatches,ipars, 
     1    srcvals,srccoefs,ifplot,fname)
      implicit real *8 (a-h,o-z)
      integer igeomtype,norder,npatches,ipars(*),ifplot
      character (len=*) fname
      real *8 srcvals(12,*), srccoefs(9,*)
      real *8, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts(:)

      real *8, pointer :: ptr1,ptr2,ptr3,ptr4
      integer, pointer :: iptr1,iptr2,iptr3,iptr4
      real *8, target :: p1(10),p2(10),p3(10),p4(10)
      real *8, allocatable, target :: triaskel(:,:,:)
      real *8, allocatable, target :: deltas(:,:)
      integer, allocatable :: isides(:)
      integer, target :: nmax,mmax

      procedure (), pointer :: xtri_geometry


      external xtri_stell_eval,xtri_sphere_eval
      
      npols = (norder+1)*(norder+2)/2
      allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols))
      allocate(wts(npols))

      call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)

      if(igeomtype.eq.1) then
        itype = 2
        allocate(triaskel(3,3,npatches))
        allocate(isides(npatches))
        npmax = npatches
        ntri = 0
        call xtri_platonic(itype, ipars(1), npmax, ntri, 
     1      triaskel, isides)

        xtri_geometry => xtri_sphere_eval
        ptr1 => triaskel(1,1,1)
        ptr2 => p2(1)
        ptr3 => p3(1)
        ptr4 => p4(1)


        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, 
     1         ptr3,ptr4, norder,'Triangulated surface of the sphere')
        endif


        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4,
     1     npols,uvs,umatr,srcvals,srccoefs)
      endif

      if(igeomtype.eq.2) then
        done = 1
        pi = atan(done)*4
        umin = 0
        umax = 2*pi
        vmin = 2*pi
        vmax = 0
        allocate(triaskel(3,3,npatches))
        nover = 0
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2),
     1     nover,npatches,npatches,triaskel)

        mmax = 2
        nmax = 1
        xtri_geometry => xtri_stell_eval

        allocate(deltas(-1:mmax,-1:nmax))
        deltas(-1,-1) = 0.17d0
        deltas(0,-1) = 0
        deltas(1,-1) = 0
        deltas(2,-1) = 0

        deltas(-1,0) = 0.11d0
        deltas(0,0) = 1
        deltas(1,0) = 4.5d0
        deltas(2,0) = -0.25d0

        deltas(-1,1) = 0
        deltas(0,1) = 0.07d0
        deltas(1,1) = 0
        deltas(2,1) = -0.45d0

        ptr1 => triaskel(1,1,1)
        ptr2 => deltas(-1,-1)
        iptr3 => mmax
        iptr4 => nmax

        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, 
     1         iptr3,iptr4, norder,
     2         'Triangulated surface of the stellarator')
        endif

        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,iptr3,iptr4,
     1     npols,uvs,umatr,srcvals,srccoefs)
      endif
      
      return  
      end

