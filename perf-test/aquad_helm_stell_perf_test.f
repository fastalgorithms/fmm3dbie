      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:),targs(:,:)
      real *8, allocatable :: wts(:)
      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)
      real *8 errs(6),ts(2)
      real *8, allocatable :: rfacs(:,:)
      character *100 fname
      integer ipars(2)
      integer, allocatable :: row_ptr(:),col_ind(:)
      integer, allocatable :: iquad(:)
      integer, allocatable :: col_ptr(:),row_ind(:)
      integer, allocatable :: iper(:)
      real *8, allocatable :: srcover(:,:),wover(:)
      complex *16, allocatable :: uval(:),dudnval(:)
      complex *16, allocatable :: sigmaover(:),slp_near(:),dlp_near(:)
      complex *16, allocatable :: slp_near_ex(:),dlp_near_ex(:)
      complex *16, allocatable :: pot(:),potslp(:),potdlp(:)
      complex *16, allocatable :: potslp2(:)

      complex *16 zk

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)
      integer, allocatable :: ixyzso(:),nfars(:)

      integer, allocatable :: ipatch_id(:),inode_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      real *8 xyz_out(3),xyz_in(3)
      complex *16, allocatable :: sigma(:)
      complex * 16 zpars(3)


      external h3d_ggq_comb,h3d_ggq_slp

      call prini(6,13)

      done = 1
      pi = atan(done)*4


      igeomtype = 2
      iasp = 3
      iref = 0
      iprec = 4


      if(iasp.eq.1) then
        ipars(1) = 5
        ipars(2) = 15
      endif

      if(iasp.eq.2) then
        ipars(1) = 6
        ipars(2) = 12
      endif

      if(iasp.eq.3) then
        ipars(1) = 9
        ipars(2) = 9
      endif

      if(iasp.eq.4) then
        ipars(1) = 12
        ipars(2) = 6
      endif

      if(iasp.eq.5) then
        ipars(1) = 15
        ipars(2) = 5
      endif

      if(iprec.eq.0) eps = 0.5001d-2
      if(iprec.eq.1) eps = 0.5001d-3
      if(iprec.eq.2) eps = 0.5001d-6
      if(iprec.eq.3) eps = 0.5001d-9
      if(iprec.eq.4) eps = 0.5001d-12

      ipars(1) = ipars(1)*2**(iref)
      ipars(2) = ipars(2)*2**(iref)

      npatches = 2*ipars(1)*ipars(2)

      norder = 8 
      npols = (norder+1)*(norder+2)/2
      write(fname,'(a,i1,a,i1,a,i1,a)') 'data/stell_aquad_iasp_',iasp,
     1  '_iref_',iref,'_norder_',norder,'.dat'


      zk = 1.0d0
      zpars(1) = zk
      zpars(2) = 1.0d0
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
        xyz_in(1) = -4.501d0
        xyz_in(2) = 1.7d-3
        xyz_in(3) = 0.00001d0

        xyz_out(1) = -3.5d0
        xyz_out(2) = 3.1d0
        xyz_out(3) = 20.1d0
      endif


      npts = npatches*npols
      allocate(srcvals(12,npts),srccoefs(9,npts))
      allocate(targs(3,npts))
      ifplot = 0



      call setup_geom(igeomtype,norder,npatches,ipars, 
     1       srcvals,srccoefs,ifplot,fname)

      allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))
      allocate(ixyzso(npatches+1),nfars(npatches))

      do i=1,npatches
        norders(i) = norder
        ixyzs(i) = 1 +(i-1)*npols
        iptype(i) = 1
      enddo

      ixyzs(npatches+1) = 1+npols*npatches
      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)



      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))
      allocate(pot(npts),potslp(npts),potdlp(npts))

      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts, 
     1     srccoefs,cms,rads)

      allocate(sigma(npts),uval(npts),dudnval(npts))

      do i=1,npts
        call h3d_slp(xyz_out,3,srcvals(1,i),0,dpars,1,zpars,0,
     1     ipars,uval(i))
        call h3d_sprime(xyz_out,12,srcvals(1,i),0,dpars,1,zpars,0,
     1     ipars,dudnval(i))
      enddo

      ndtarg = 3
     
      do i=1,npts
        targs(1,i) = srcvals(1,i)
        targs(2,i) = srcvals(2,i)
        targs(3,i) = srcvals(3,i)
      enddo

      allocate(ipatch_id(npts),uvs_targ(2,npts))
      do i=1,npts
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo

      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts, 
     1         ipatch_id,uvs_targ)

 
c
c    find near field
c
      iptype = 1
      call get_rfacs(norder,iptype,rfac,rfac0)
      rfac = 2.75d0
      rfac0 = 1.25d0
      do i=1,npatches 
        rad_near(i) = rads(i)*rfac
      enddo
      

      call findnearmem(cms,npatches,rad_near,ndtarg,targs,npts,nnz)

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,ndtarg,targs,npts,row_ptr, 
     1        col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind,
     1         iquad)

      nquad = iquad(nnz+1)-1
      allocate(slp_near(nquad),dlp_near(nquad))
      allocate(slp_near_ex(nquad),dlp_near_ex(nquad))


      
      ndtarg = 3
      ikerorder = -1



      do i=1,nquad
        slp_near(i) = 0
        dlp_near(i) = 0
      enddo

      print *, "done reading1"



      call cpu_time(t1)
C$       t1 = omp_get_wtime()      

      zpars(1) = zk
      zpars(2) = 1.0d0
      zpars(3) = 0.0d0

      iquadtype = 1
      call prin2('eps=*',eps,1)
      call prinf('npatches=*',npatches,1)

cc      goto 1111
      call getnearquad_helm_comb_dir(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     1      ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,
     1      iquad,
     1      rfac0,nquad,slp_near)
      call prin2('done with slp*',i,0)

      
      zpars(2) = 0.0d0
      zpars(3) = 1.0d0
      call getnearquad_helm_comb_dir(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     1      ipatch_id,uvs_targ,eps,zpars,iquadtype,
     1      nnz,row_ptr,col_ind,iquad,
     1      rfac0,nquad,dlp_near)
      
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      
      tquadgen = t2-t1
      print *, "done quadgen"
      call prin2('tquadgen=*',tquadgen,1)


      allocate(col_ptr(npatches+1),row_ind(nnz),iper(nnz))
      call rsc_to_csc(npatches,npts,nnz,row_ptr,col_ind,col_ptr,
     1   row_ind,iper)
      call prinf('col_ptr=*',col_ptr,20)
      call prinf('row_ind=*',row_ind,20)

      open(unit=23,file=fname,action='readwrite',form='unformatted',
     1  access='stream')
      read(unit=23) slp_near_ex
      read(unit=23) dlp_near_ex
      call prin2('slp_near=*',slp_near,24)
      call prin2('slp_near_ex=*',slp_near_ex,24)

c
c  now compute accuracy
c
      rrat = 0.0d0
      do ipatch=1,npatches
        
        rsc = wts(ixyzs(ipatch)) 
        do j=ixyzs(ipatch),ixyzs(ipatch+1)-1
          if(wts(j).lt.rsc) rsc = wts(j)
        enddo
        rsc = sqrt(rsc)
        npols = ixyzs(ipatch+1)-ixyzs(ipatch)
        eps_test = eps*rsc
        eps_test = max(eps_test,1.0d-14)
        do j=row_ind(ipatch),row_ind(ipatch+1)-1
          jquadstart = iquad(iper(j))
          do l=1,npols
            erra = abs(slp_near_ex(jquadstart+l-1) - 
     1          slp_near(jquadstart+l-1))
            errb = abs(dlp_near_ex(jquadstart+l-1) - 
     1          dlp_near(jquadstart+l-1))
            if(errb.gt.erra) erra = errb
            rrat0 = erra/eps_test
            if(rrat0.gt.rrat) rrat = rrat0 
          enddo
        enddo
      enddo

      call get_quadparams_adap(eps,nqorder,eps_adap,nlev,nqorder_f)

      open(unit=33,file='data/res-sum.dat',access='append')
 1231 format(2x,i1,2x,e11.5,2x,i2,2x,e11.5,2(2x,i2),2(2x,e11.5))      
      write(33,1231) iasp,eps,nqorder,eps_adap,nlev,nqorder_f,tquadgen,
     1     rrat 
      write(*,1231) iasp,eps,nqorder,eps_adap,nlev,nqorder_f,tquadgen,
     1     rrat 






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

