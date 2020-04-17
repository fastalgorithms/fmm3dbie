
c
c
c     This file has the following user callable routines
c
c        getnearquad_ggq_compact_guru - guru interface for 
c         computing near field quadrature for a compact kernel
c         (for e.g. single layer potentials,
c            double layer potentials for smooth domains)
c        
c        getnearquad_ggq_pv_guru - guru interface for computing near
c          field quadrature for principal value kernels
c          (for e.g. tangential derivatives of single layer potentials)
c
c
c
c
c

      subroutine getnearquad_ggq_compact_guru(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targvals,
     2   ipatch_id,uvs_targ,eps,fker,ndd,dpars,ndz,zpars,ndi,ipars,
     3   nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear)
c
c       this subroutine generates the near field quadrature
c       for the given kernel which is assumed to be
c       a compact integral operator
c       where the near field is specified by the user 
c       in row sparse compressed format.
c
c
c       The quadrature is computed by the following strategy
c        targets within a sphere of radius rfac0*rs
c        of a chunk centroid is handled using adaptive integration
c        where rs is the radius of the bounding sphere
c        for the patch
c  
c       All other targets in the near field are handled via
c        oversampled quadrature
c
c        
c
c       input:
c         npatches - integer
c            number of patches
c
c         norders - integer(npatches)
c            order of discretization on each patch
c
c         ixyzs - integer(npatches+1)
c            ixyzs(i) denotes the starting location in srccoefs,
c               and srcvals array corresponding to patch i
c   
c         iptype - integer(npatches)
c            type of patch
c             iptype = 1, triangular patch discretized using RV nodes
c
c         npts - integer
c            total number of discretization points on the boundary
c 
c         srccoefs - real *8 (9,npts)
c            koornwinder expansion coefficients of xyz, dxyz/du,
c            and dxyz/dv on each patch. 
c            For each point srccoefs(1:3,i) is xyz info
c                           srccoefs(4:6,i) is dxyz/du info
c                           srccoefs(7:9,i) is dxyz/dv info
c
c          srcvals - real *8 (12,npts)
c             xyz(u,v) and derivative info sampled at the 
c             discretization nodes on the surface
c             srcvals(1:3,i) - xyz info
c             srcvals(4:6,i) - dxyz/du info
c             srcvals(7:9,i) - dxyz/dv info
c
c
c          ndtarg - integer
c             leading dimension of target array
c
c          ntarg - integer
c             number of targets
c
c          targvals - real *8 (ndtarg,ntarg)
c             boundary target info
c
c          ipatch_id - integer(ntarg)
c             patch id of target, patch_id = -1, if target off-surface
c         
c          uvs_targ - real *8 (2,ntarg)
c            local uv coordinates on patch if target on surface
c
c          eps - real *8
c             precision requested
c
c          zpars - complex *16 (*)
c
c           nnz - integer *8
c             number of source patch-> target interactions in the near field
c 
c           row_ptr - integer(npts+1)
c              row_ptr(i) is the pointer
c              to col_ind array where list of relevant source patches
c              for target i start
c
c           col_ind - integer (nnz)
c               list of source patches relevant for all targets, sorted
c               by the target number
c
c           iquad - integer(nnz+1)
c               location in wnear array where quadrature for col_ind(i)
c               starts
c
c           rfac0 - integer
c               radius parameter for near field
c
c           nquad - integer
c               number of entries in wnear
c
c        output
c            wnear - complex *16(nquad)
c               the desired near field quadrature
c               
c

      implicit real *8 (a-h,o-z)
      integer ipars(*)
      integer ndtarg
      integer ntrimax
      integer npatches,norders(npatches),npts
      integer ixyzs(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps
      integer ntarg,ipatch_id(ntarg)
      real *8 uvs_targ(2,ntarg)
      real *8 targvals(ndtarg,ntarg)
      real *8, allocatable :: cms(:,:),rads(:)
      real *8, allocatable :: targ_near(:,:),targ_far(:,:)
      integer, allocatable :: iind_near(:),iind_far(:)
      complex *16 zpars(3)
      integer nnz
      integer row_ptr(ntarg+1),col_ind(nnz),iquad(nnz+1)
      complex *16 wnear(nquad)

      real *8, allocatable :: umatr(:,:),vmatr(:,:),uvs(:,:),wts(:)

c
c        temporary variables
c
      integer, allocatable :: col_ptr(:),row_ind(:),iper(:)
      complex *16, allocatable :: sints_n(:,:),svtmp_n(:,:)
      complex *16, allocatable :: sints_f(:,:),svtmp_f(:,:)
      complex *16, allocatable :: svtmp2(:,:)

      real *8, allocatable :: rad(:),xyztarg2(:,:)
      integer ilrad(3),nlrad(3),irad

      real *8, allocatable :: qnodes(:,:),qwts(:)
      real *8 ra

      complex *16 ff1,ff2,cra1,cra2
      integer nlev, nqorder_f
      real *8 rfac0

      character *200 dirname
      character *300 fname


      external fker
c      
c       determine number of discretization per node
c       and get the relevant interpolation matrices
c



c
c
c       transpose the row sparse compressed format
c       to the column sparse compressed format needed for the
c       the adaptive integration routines
c
      allocate(row_ind(nnz),col_ptr(npatches+1),iper(nnz))

      call rsc_to_csc(npatches,ntarg,nnz,row_ptr,col_ind,col_ptr,
     1   row_ind,iper)


      allocate(cms(3,npatches),rads(npatches))


c
c        find chunk radisu and centroid again
c      note: this cost can be saved if needed
c      by passing it in as arguments
c

      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts,
     1     srccoefs,cms,rads)

      
c
c
c       get adaptive integration + bremer self quadrature parameters
c
      norder0 = norder
      nqorder = 8
      eps_adap = eps


c      suppress computation of 
c      some relevant metrics for analyzing adaptive
c      quadrature performance
c       
      ifmetric = 0
      rn1 = 0
      n2 = 0

c       use xiao gimbutas nodes      
      intype = 2
c       use adaptive integration (see documentation of ctriaints)     
      istrat = 2
c       relevant parameter for istrat =1/3      
      ifp = 0
c      number of source patches to be processed per batch      
      ntest0 = 1
      itargptr = 1

      ntrimax = 3000

      rfac = 1.0d0
c
c
c    Near quadrature params
c      nqorder is the order of quadrature nodes
c      used on each triangle for handling the adaptive
c      integration. Depends on order of discretization
c      and accuracy requested
c
c      norder0 is the smallest integer multiple of 4
c      greater than norder, needed for bremer 
c      self-quadrature scheme
c
c      eps_adap is the precision set to adaptive 
c      integration routines, it is slightly greater
c      than eps requested, and that depends nqorder
c
c
c   Far quadrature params
c      nlev is the number of refinements of the standard simplex
c
c      nqorder_f is the order of XG nodes on each of the smaller
c      triangles
c       
c      


      call get_quadparams_adap(eps,nqorder,eps_adap,nlev,nqorder_f)

      call triasymq_pts(nqorder_f,nnodes)
      
      ntri_f = 4**nlev
      npts_f = ntri_f*nnodes
      allocate(qnodes(2,npts_f),qwts(npts_f))

      call gen_xg_unif_nodes(nlev,nqorder_f,nnodes,npts_f,qnodes,qwts)


c
c       setup self quadrature for bremer
c
      
      
      call getenv("S3D_INSTALL_DIR",dirname)
      fname = trim(dirname)//
     1   '/src/quadratures/ggq-self-quads/radial4-8-12.bin'
      open(unit=171,file=fname,action='read',
     1     form='unformatted',access='stream',status='old')
      read(unit=171) ilrad
      read(unit=171) nlrad
      read(unit=171) lrad
      allocate(rad(lrad))
      read(unit=171) rad
      close(171)

      t1 = second()
C$        t1 = omp_get_wtime()

C$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(DYNAMIC)
C$OMP$PRIVATE(ipatch,ntarg2,xyztarg2,sints_n,sints_f,svtmp_n,svtmp_f)
C$OMP$PRIVATE(ii,ii2,i,jpatch,svtmp2,iiif,l,ntarg2m)
C$OMP$PRIVATE(j,iii,istart,itarg2,iqstart,jpt,jtarg)
C$OMP$PRIVATE(targ_near,targ_far,iind_near,iind_far,rr)
C$OMP$PRIVATE(ntarg_f,ntarg_n,npols,norder,irad)
C$OMP$PRIVATE(uvs,umatr,vmatr,wts)

      do ipatch=1,npatches
        
        npols = ixyzs(ipatch+1)-ixyzs(ipatch)
        norder = norders(ipatch)

c
c         
c
        if(norder.le.4) then
          irad = 1
        else if(norder.le.8) then
          irad = 2
        else
          irad = 3
        endif

        allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols))
        allocate(wts(npols))
        if(iptype(ipatch).eq.1) 
     1     call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)

        ntarg2m = col_ptr(ipatch+1)-col_ptr(ipatch) 
        allocate(xyztarg2(ndtarg,ntarg2m),svtmp2(npols,ntarg2m))
        allocate(targ_near(ndtarg,ntarg2m),iind_near(ntarg2m))
        allocate(targ_far(ndtarg,ntarg2m),iind_far(ntarg2m))

        ii = 1
        ii2 = 1
        do i=col_ptr(ipatch),col_ptr(ipatch+1)-1
          jtarg = row_ind(i)
          jpatch = ipatch_id(jtarg)

          if(ipatch.ne.jpatch) then
            do iiif=1,ndtarg
               xyztarg2(iiif,ii2) = targvals(iiif,jtarg)
            enddo
            ii2 = ii2 + 1
          endif
        enddo
        ntarg2 = ii2-1

c
c          split near far targs 
c
       rr = rfac0*rads(ipatch)
       call get_near_far_split_pt(ndtarg,ntarg2,xyztarg2,rr,
     1        cms(1,ipatch),ntarg_f,ntarg_n,targ_far,targ_near,
     2        iind_far,iind_near)

        allocate(sints_n(npols,ntarg_n))
        allocate(svtmp_n(npols,ntarg_n))

        allocate(sints_f(npols,ntarg_f))
        allocate(svtmp_f(npols,ntarg_f))

        istart = ixyzs(ipatch)

c
c       fill out near part of single layer
c

        if(iptype(ipatch).eq.1) 
     1      call ctriaints(eps_adap,istrat,intype,ntest0,norder,npols,
     1      srccoefs(1,istart),ndtarg,ntarg_n,targ_near,ifp,xyztarg2,
     2      itargptr,ntarg_n,norder,npols,fker,ndd,dpars,ndz,zpars,ndi,
     3      ipars,
     3      nqorder,ntrimax,rfac,sints_n,ifmetric,
     3      rn1,n2)
        
        call zrmatmatt_slow(ntarg_n,npols,npols,sints_n,umatr,svtmp_n)
c
c
c       fill out far part of single layer
c

        if(iptype(ipatch).eq.1) 
     1     call ctriaints_wnodes(ntest0,norder,npols,
     1     srccoefs(1,istart),ndtarg,ntarg_f,targ_far,
     2      itargptr,ntarg_f,norder,npols,fker,ndd,dpars,ndz,zpars,
     3      ndi,ipars,
     3      npts_f,qnodes,qwts,sints_f)
        
        call zrmatmatt_slow(ntarg_f,npols,npols,sints_f,umatr,svtmp_f)
c
c       combine svtmp_f, svtmp_n to fill out svtmp2
c
c

        do i=1,ntarg_f
          ii = iind_far(i)
          do j=1,npols
            svtmp2(j,ii) = svtmp_f(j,i)
          enddo
        enddo

        do i=1,ntarg_n
          ii = iind_near(i)
          do j=1,npols
            svtmp2(j,ii) = svtmp_n(j,i)
          enddo
        enddo

c
c
c       fill out relevant sections of wnear and
c
        itarg2 = 1
        do i=col_ptr(ipatch),col_ptr(ipatch+1)-1
          jtarg = row_ind(i)
          jpatch = ipatch_id(jtarg)
          iqstart = iquad(iper(i))-1

          if(ipatch.eq.jpatch) then
            call get_ggq_self_quad_pt(norder,npols,uvs_targ(1,jtarg),
     1      umatr,srccoefs(1,istart),ndtarg,targvals(1,jtarg),
     2      rad(ilrad(irad)),nlrad(irad),fker,ndd,dpars,ndz,zpars,
     2      ndi,pars,
     2      wnear(iqstart+1))
          endif

          if(ipatch.ne.jpatch) then
            do l=1,npols
              wnear(iqstart+l) = svtmp2(l,itarg2)
            enddo
            itarg2 = itarg2 + 1
          endif
        enddo

        deallocate(xyztarg2,svtmp2,svtmp_f,svtmp_n,sints_f,sints_n)
        deallocate(targ_near,targ_far,iind_near,iind_far)
        deallocate(uvs,umatr,vmatr,wts)
      enddo
C$OMP END PARALLEL DO      

      t2 = second()
C$      t2 = omp_get_wtime()     


      return
      end
c
c
c
c
c
c
c
c
c
c

      subroutine getnearquad_ggq_pv_guru(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targvals,
     2   ipatch_id,uvs_targ,eps,fker,ndd,dpars,ndz,zpars,ndi,ipars,
     3   nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear)
c
c       this subroutine generates the near field quadrature
c       for the given kernel which is assumed to be
c       a compact integral operator
c       where the near field is specified by the user 
c       in row sparse compressed format.
c
c
c       The quadrature is computed by the following strategy
c        targets within a sphere of radius rfac0*rs
c        of a chunk centroid is handled using adaptive integration
c        where rs is the radius of the bounding sphere
c        for the patch
c  
c       All other targets in the near field are handled via
c        oversampled quadrature
c
c        
c
c       input:
c         npatches - integer
c            number of patches
c
c         norders - integer(npatches)
c            order of discretization on each patch
c
c         ixyzs - integer(npatches+1)
c            ixyzs(i) denotes the starting location in srccoefs,
c               and srcvals array corresponding to patch i
c   
c         iptype - integer(npatches)
c            type of patch
c             iptype = 1, triangular patch discretized using RV nodes
c
c         npts - integer
c            total number of discretization points on the boundary
c 
c         srccoefs - real *8 (9,npts)
c            koornwinder expansion coefficients of xyz, dxyz/du,
c            and dxyz/dv on each patch. 
c            For each point srccoefs(1:3,i) is xyz info
c                           srccoefs(4:6,i) is dxyz/du info
c                           srccoefs(7:9,i) is dxyz/dv info
c
c          srcvals - real *8 (12,npts)
c             xyz(u,v) and derivative info sampled at the 
c             discretization nodes on the surface
c             srcvals(1:3,i) - xyz info
c             srcvals(4:6,i) - dxyz/du info
c             srcvals(7:9,i) - dxyz/dv info
c
c
c          ndtarg - integer
c             leading dimension of target array
c
c          ntarg - integer
c             number of targets
c
c          targvals - real *8 (ndtarg,ntarg)
c             boundary target info
c
c          ipatch_id - integer(ntarg)
c             patch id of target, patch_id = -1, if target off-surface
c         
c          uvs_targ - real *8 (2,ntarg)
c            local uv coordinates on patch if target on surface
c
c          eps - real *8
c             precision requested
c
c          zpars - complex *16 (*)
c
c           nnz - integer *8
c             number of source patch-> target interactions in the near field
c 
c           row_ptr - integer(npts+1)
c              row_ptr(i) is the pointer
c              to col_ind array where list of relevant source patches
c              for target i start
c
c           col_ind - integer (nnz)
c               list of source patches relevant for all targets, sorted
c               by the target number
c
c           iquad - integer(nnz+1)
c               location in wnear array where quadrature for col_ind(i)
c               starts
c
c           rfac0 - integer
c               radius parameter for near field
c
c           nquad - integer
c               number of entries in wnear
c
c        output
c            wnear - complex *16(nquad)
c               the desired near field quadrature
c               
c

      implicit real *8 (a-h,o-z)
      integer ipars(*)
      integer ndtarg
      integer ntrimax
      integer npatches,norders(npatches),npts
      integer ixyzs(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps
      integer ntarg,ipatch_id(ntarg)
      real *8 uvs_targ(2,ntarg)
      real *8 targvals(ndtarg,ntarg)
      real *8, allocatable :: cms(:,:),rads(:)
      real *8, allocatable :: targ_near(:,:),targ_far(:,:)
      integer, allocatable :: iind_near(:),iind_far(:)
      complex *16 zpars(3)
      integer nnz
      integer row_ptr(ntarg+1),col_ind(nnz),iquad(nnz+1)
      complex *16 wnear(nquad)

      real *8, allocatable :: umatr(:,:),vmatr(:,:),uvs(:,:),wts(:)
      character *200 dirname
      character *100 fname

c
c        temporary variables
c
      integer, allocatable :: col_ptr(:),row_ind(:),iper(:)
      complex *16, allocatable :: sints_n(:,:),svtmp_n(:,:)
      complex *16, allocatable :: sints_f(:,:),svtmp_f(:,:)
      complex *16, allocatable :: svtmp2(:,:)

      real *8, allocatable :: rad(:),xyztarg2(:,:)
      integer ilrad(3),nlrad(3),irad

      real *8, allocatable :: qnodes(:,:),qwts(:)
      real *8 ra

      complex *16 ff1,ff2,cra1,cra2
      integer nlev, nqorder_f
      real *8 rfac0


      external fker
c      
c       determine number of discretization per node
c       and get the relevant interpolation matrices
c



c
c
c       transpose the row sparse compressed format
c       to the column sparse compressed format needed for the
c       the adaptive integration routines
c
      allocate(row_ind(nnz),col_ptr(npatches+1),iper(nnz))

      call rsc_to_csc(npatches,ntarg,nnz,row_ptr,col_ind,col_ptr,
     1   row_ind,iper)


      allocate(cms(3,npatches),rads(npatches))


c
c        find chunk radisu and centroid again
c      note: this cost can be saved if needed
c      by passing it in as arguments
c

      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts,
     1     srccoefs,cms,rads)

      
c
c
c       get adaptive integration + bremer self quadrature parameters
c
      norder0 = norder
      nqorder = 8
      eps_adap = eps


c      suppress computation of 
c      some relevant metrics for analyzing adaptive
c      quadrature performance
c       
      ifmetric = 0
      rn1 = 0
      n2 = 0

c       use xiao gimbutas nodes      
      intype = 2
c       use adaptive integration (see documentation of ctriaints)     
      istrat = 2
c       relevant parameter for istrat =1/3      
      ifp = 0
c      number of source patches to be processed per batch      
      ntest0 = 1
      itargptr = 1

      ntrimax = 3000

      rfac = 1.0d0
c
c
c    Near quadrature params
c      nqorder is the order of quadrature nodes
c      used on each triangle for handling the adaptive
c      integration. Depends on order of discretization
c      and accuracy requested
c
c      norder0 is the smallest integer multiple of 4
c      greater than norder, needed for bremer 
c      self-quadrature scheme
c
c      eps_adap is the precision set to adaptive 
c      integration routines, it is slightly greater
c      than eps requested, and that depends nqorder
c
c
c   Far quadrature params
c      nlev is the number of refinements of the standard simplex
c
c      nqorder_f is the order of XG nodes on each of the smaller
c      triangles
c       
c      


      call get_quadparams_adap(eps,nqorder,eps_adap,nlev,nqorder_f)

      call triasymq_pts(nqorder_f,nnodes)
      
      ntri_f = 4**nlev
      npts_f = ntri_f*nnodes
      allocate(qnodes(2,npts_f),qwts(npts_f))

      call gen_xg_unif_nodes(nlev,nqorder_f,nnodes,npts_f,qnodes,qwts)



      
c
c       setup self quadrature for bremer
c
      
      call getenv("S3D_INSTALL_DIR",dirname)
      fname = trim(dirname)//
     1   '/src/quadratures/ggq-self-quads/pvradial4-8-12.bin'
      open(unit=171,file=fname,action='read',
     1     form='unformatted',access='stream',status='OLD')
      read(unit=171) ilrad
      read(unit=171) nlrad
      read(unit=171) lrad
      allocate(rad(lrad))
      read(unit=171) rad
      close(171)



      t1 = second()
C$        t1 = omp_get_wtime()

C$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(DYNAMIC)
C$OMP$PRIVATE(ipatch,ntarg2,xyztarg2,sints_n,sints_f,svtmp_n,svtmp_f)
C$OMP$PRIVATE(ii,ii2,i,jpatch,svtmp2,iiif,l,ntarg2m)
C$OMP$PRIVATE(j,iii,istart,itarg2,iqstart,jpt,jtarg)
C$OMP$PRIVATE(targ_near,targ_far,iind_near,iind_far,rr)
C$OMP$PRIVATE(ntarg_f,ntarg_n,npols,norder,irad)
C$OMP$PRIVATE(uvs,umatr,vmatr,wts)

      do ipatch=1,npatches
        
        npols = ixyzs(ipatch+1)-ixyzs(ipatch)
        norder = norders(ipatch)

c
c         
c
        if(norder.le.4) then
          irad = 1
        else if(norder.le.8) then
          irad = 2
        else
          irad = 3
        endif

        allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols))
        allocate(wts(npols))
        if(iptype(ipatch).eq.1) 
     1     call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)

        ntarg2m = col_ptr(ipatch+1)-col_ptr(ipatch) 
        allocate(xyztarg2(ndtarg,ntarg2m),svtmp2(npols,ntarg2m))
        allocate(targ_near(ndtarg,ntarg2m),iind_near(ntarg2m))
        allocate(targ_far(ndtarg,ntarg2m),iind_far(ntarg2m))

        ii = 1
        ii2 = 1
        do i=col_ptr(ipatch),col_ptr(ipatch+1)-1
          jtarg = row_ind(i)
          jpatch = ipatch_id(jtarg)

          if(ipatch.ne.jpatch) then
            do iiif=1,ndtarg
               xyztarg2(iiif,ii2) = targvals(iiif,jtarg)
            enddo
            ii2 = ii2 + 1
          endif
        enddo
        ntarg2 = ii2-1

c
c          split near far targs 
c
       rr = rfac0*rads(ipatch)
       call get_near_far_split_pt(ndtarg,ntarg2,xyztarg2,rr,
     1        cms(1,ipatch),ntarg_f,ntarg_n,targ_far,targ_near,
     2        iind_far,iind_near)

        allocate(sints_n(npols,ntarg_n))
        allocate(svtmp_n(npols,ntarg_n))

        allocate(sints_f(npols,ntarg_f))
        allocate(svtmp_f(npols,ntarg_f))

        istart = ixyzs(ipatch)

c
c       fill out near part of single layer
c

        if(iptype(ipatch).eq.1) 
     1      call ctriaints(eps_adap,istrat,intype,ntest0,norder,npols,
     1      srccoefs(1,istart),ndtarg,ntarg_n,targ_near,ifp,xyztarg2,
     2      itargptr,ntarg_n,norder,npols,fker,ndd,dpars,ndz,zpars,
     3      ndi,ipars,
     3      nqorder,ntrimax,rfac,sints_n,ifmetric,
     3      rn1,n2)
        
        call zrmatmatt_slow(ntarg_n,npols,npols,sints_n,umatr,svtmp_n)
c
c
c       fill out far part of single layer
c

        if(iptype(ipatch).eq.1) 
     1     call ctriaints_wnodes(ntest0,norder,npols,
     1     srccoefs(1,istart),ndtarg,ntarg_f,targ_far,
     2      itargptr,ntarg_f,norder,npols,fker,ndd,dpars,ndz,zpars,
     3      ndi,ipars,
     3      npts_f,qnodes,qwts,sints_f)
        
        call zrmatmatt_slow(ntarg_f,npols,npols,sints_f,umatr,svtmp_f)
c
c       combine svtmp_f, svtmp_n to fill out svtmp2
c
c

        do i=1,ntarg_f
          ii = iind_far(i)
          do j=1,npols
            svtmp2(j,ii) = svtmp_f(j,i)
          enddo
        enddo

        do i=1,ntarg_n
          ii = iind_near(i)
          do j=1,npols
            svtmp2(j,ii) = svtmp_n(j,i)
          enddo
        enddo

c
c
c       fill out relevant sections of wnear and
c
        itarg2 = 1
        do i=col_ptr(ipatch),col_ptr(ipatch+1)-1
          jtarg = row_ind(i)
          jpatch = ipatch_id(jtarg)
          iqstart = iquad(iper(i))-1

          if(ipatch.eq.jpatch) then
            call get_ggq_self_quad_pv_pt(norder,npols,
     1      uvs_targ(1,jtarg),
     1      umatr,srccoefs(1,istart),ndtarg,targvals(1,jtarg),
     2      rad(ilrad(irad)),nlrad(irad),fker,ndd,dpars,ndz,zpars,ndi,
     2      ipars,wnear(iqstart+1))
          endif

          if(ipatch.ne.jpatch) then
            do l=1,npols
              wnear(iqstart+l) = svtmp2(l,itarg2)
            enddo
            itarg2 = itarg2 + 1
          endif
        enddo

        deallocate(xyztarg2,svtmp2,svtmp_f,svtmp_n,sints_f,sints_n)
        deallocate(targ_near,targ_far,iind_near,iind_far)
        deallocate(uvs,umatr,vmatr,wts)
      enddo
C$OMP END PARALLEL DO      

      t2 = second()
C$      t2 = omp_get_wtime()     


      return
      end
c
c
c
c
c
c
c

      subroutine get_ggq_self_quad(norder,npols,uvs,umat,vmat,
     1   srccoefs,rad,lrad,
     1   fker,ndd,dpars,ndz,zpars,ndi,ipars,zquad)
c
c
c       this is the replacement for xtri_zself 
c       for implementing bremer quadrature in ytri. 
c       
c       the main difference is the calling sequence
c       of the subroutine for getting the kernel
c       this is made consistent with adaptive integration
c       on triangles.
c
c
c        The calling sequence for fker is of the form
c          fker(srcinfo,targinfo,
c
c


      implicit real *8(a-h,o-z)
      integer norder,npols,lrad
      real *8 rad(*)
      real *8 srccoefs(9,npols)
      real *8 uvs(2,npols),umat(npols,npols)
      real *8 vmat(npols,npols)
      complex *16 zquad(*)
      integer ndd,ndz,ndi
      integer ipars(ndi)
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      real *8, allocatable :: srcvals(:,:)
      real *8, allocatable :: srctmp(:,:)
      real *8, allocatable :: qwts(:),sigvals(:,:)
      real *8, allocatable :: xs(:),ys(:),ws(:)
      real *8 uv(2),verts(2,3)
      real *8 alpha,beta
      complex *16 fval
      complex *16, allocatable :: fint(:,:),fcoefs(:),fvals(:)
      character *1 transa,transb

      external fker

      nmax = 3000
      allocate(ws(nmax),xs(nmax),ys(nmax))
      allocate(srcvals(9,npols),fcoefs(npols),fvals(npols))
      verts(1,1) = 0
      verts(2,1) = 0
      verts(1,2) = 1
      verts(2,2) = 0
      verts(1,3) = 0
      verts(2,3) = 1
      allocate(fint(npols,npols))

cc      call prin2('zpars=*',zpars,2)
cc
cc      call prin2('srccoefs=*',srccoefs,9*npols)
c
c
c        check for optimizations in this part of the code
c
      do i=1,npols
        do k=1,9
          srcvals(k,i) = 0
        enddo

        do j=1,npols
          do k=1,9
            srcvals(k,i) = srcvals(k,i) + vmat(i,j)*srccoefs(k,j)
          enddo
        enddo
      enddo



      do inode=1,npols
        do j=1,npols
          fint(j,inode) = 0
          fcoefs(j) = 0
        enddo
        u = uvs(1,inode)
        v = uvs(2,inode)
        ns = 0
        call self_quadrature_new(ier,rad,verts,uvs(1,inode),
     1     uvs(2,inode),srcvals(4,inode),ns,xs,ys,ws)
     
        allocate(srctmp(12,ns),qwts(ns))
        allocate(sigvals(npols,ns))

        do i=1,ns
          uv(1) = xs(i)
          uv(2) = ys(i)

          call koorn_pols(uv,norder,npols,sigvals(1,i))
        enddo

        transa = 'N'
        transb = 'N'
        alpha = 1
        beta = 0
        lda = 9
        ldb = npols
        ldc = 12

        call dgemm(transa,transb,9,ns,npols,alpha,
     1        srccoefs,lda,sigvals,ldb,beta,srctmp,ldc)

cc        call dmatmat(3,npols,xyzcoefs,ns,sigvals,xyztmp)
cc        call dmatmat(6,npols,xyzcoefs,ns,sigvals,dtmp)


c
        do i=1,ns
          call cross_prod3d(srctmp(4,i),srctmp(7,i),srctmp(10,i))
          rr = sqrt(srctmp(10,i)**2+srctmp(11,i)**2+srctmp(12,i)**2)
          qwts(i) = rr*ws(i)
          srctmp(10,i) = srctmp(10,i)/rr
          srctmp(11,i) = srctmp(11,i)/rr
          srctmp(12,i) = srctmp(12,i)/rr

          call fker(srctmp(1,i),ndtarg,srcvals(1,inode),ndd,dpars,
     2       ndz,zpars,ndi,ipars,fval)
           
          do j=1,npols
            fint(j,inode) = fint(j,inode) + fval*sigvals(j,i)*qwts(i)
          enddo
        enddo

        deallocate(srctmp,qwts,sigvals)
      enddo

      call zrmatmatt_slow(npols,npols,npols,fint,umat,zquad)


      return
      end
c
c
c
c
c
c

      subroutine get_ggq_self_quad_pt(norder,npols,uvs,umat,
     1   srccoefs,ndtarg,targ,rad,lrad,
     1   fker,ndd,dpars,ndz,zpars,ndi,ipars,zquad)
c
c
c       this is the replacement for xtri_zself 
c       for implementing bremer quadrature in ytri. 
c       
c       the main difference is the calling sequence
c       of the subroutine for getting the kernel
c       this is made consistent with adaptive integration
c       on triangles.
c
c
c        The calling sequence for fker is of the form
c          fker(srcinfo,targinfo,
c
c


      implicit real *8(a-h,o-z)
      integer norder,npols,lrad,ndtarg
      real *8 rad(*)
      real *8 srccoefs(9,npols),srcvals(9)
      real *8, allocatable :: sigvalstmp(:)
      real *8 targ(ndtarg)
      real *8 uvs(2),umat(npols,npols)
      complex *16 zquad(*)
      integer ndi,ndd,ndz
      integer ipars(ndi)
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      real *8, allocatable :: srctmp(:,:)
      real *8, allocatable :: qwts(:),sigvals(:,:)
      real *8, allocatable :: xs(:),ys(:),ws(:)
      real *8 uv(2),verts(2,3)
      real *8 alpha,beta
      complex *16 fval
      complex *16, allocatable :: fint(:),fvals(:)
      character *1 transa,transb

      external fker

      nmax = 10000
      allocate(ws(nmax),xs(nmax),ys(nmax))
      allocate(fvals(npols))
      verts(1,1) = 0
      verts(2,1) = 0
      verts(1,2) = 1
      verts(2,2) = 0
      verts(1,3) = 0
      verts(2,3) = 1
      allocate(fint(npols),sigvalstmp(npols))

c
c       compute all source info at target point on patch
c
      call koorn_pols(uvs,norder,npols,sigvalstmp)

      alpha = 1.0d0
      beta = 0.0d0
      call dgemv('n',9,npols,alpha,srccoefs,9,sigvalstmp,1,beta,
     1      srcvals,1)
      

      do j=1,npols
        fint(j) = 0
      enddo
      ns = 0
      call self_quadrature_new(ier,rad,verts,uvs(1),
     1  uvs(2),srcvals(4),ns,xs,ys,ws)
      
     
      allocate(srctmp(12,ns),qwts(ns))
      allocate(sigvals(npols,ns))

      do i=1,ns
        uv(1) = xs(i)
        uv(2) = ys(i)

        call koorn_pols(uv,norder,npols,sigvals(1,i))
      enddo

      transa = 'N'
      transb = 'N'
      alpha = 1
      beta = 0
      lda = 9
      ldb = npols
      ldc = 12

      call dgemm(transa,transb,9,ns,npols,alpha,
     1      srccoefs,lda,sigvals,ldb,beta,srctmp,ldc)


c
      do i=1,ns
        call cross_prod3d(srctmp(4,i),srctmp(7,i),srctmp(10,i))
        rr = sqrt(srctmp(10,i)**2+srctmp(11,i)**2+srctmp(12,i)**2)
        qwts(i) = rr*ws(i)
        srctmp(10,i) = srctmp(10,i)/rr
        srctmp(11,i) = srctmp(11,i)/rr
        srctmp(12,i) = srctmp(12,i)/rr

        call fker(srctmp(1,i),ndtarg,targ,ndd,dpars,
     2         ndz,zpars,ndi,ipars,fval)
           
        do j=1,npols
          fint(j) = fint(j) + fval*sigvals(j,i)*qwts(i)
        enddo
      enddo

      deallocate(srctmp,qwts,sigvals)

      call zrmatmatt_slow(1,npols,npols,fint,umat,zquad)


      return
      end
c
c
c
c
c
c
c
c

      subroutine get_ggq_self_quad_pt_vec(norder,npols,uvs,umat,
     1   srccoefs,ndtarg,targ,rad,lrad,
     1   fker,nd,ndd,dpars,ndz,zpars,ndi,ipars,zquad)
c
c
c       Vectorized version of get_ggq_self_quad_pt_vec
c
c
c        The calling sequence for fker is of the form
c          fker(nd,srcinfo,targinfo,
c
c


      implicit real *8(a-h,o-z)
      integer norder,npols,lrad,ndtarg,nd
      real *8 rad(*)
      real *8 srccoefs(9,npols),srcvals(9)
      real *8, allocatable :: sigvalstmp(:)
      real *8 targ(ndtarg)
      real *8 uvs(2),umat(npols,npols)
      complex *16 zquad(nd,*)
      integer ipars(ndi)
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      real *8, allocatable :: srctmp(:,:)
      real *8, allocatable :: qwts(:),sigvals(:,:)
      real *8, allocatable :: xs(:),ys(:),ws(:)
      real *8 uv(2),verts(2,3)
      real *8 alpha,beta
      complex *16, allocatable :: fint(:,:),fvals(:,:),fval(:)
      integer idim
      character *1 transa,transb

      external fker

      nmax = 3000
      allocate(ws(nmax),xs(nmax),ys(nmax))
      allocate(fvals(nd,npols),fval(nd))
      verts(1,1) = 0
      verts(2,1) = 0
      verts(1,2) = 1
      verts(2,2) = 0
      verts(1,3) = 0
      verts(2,3) = 1
      allocate(fint(nd,npols),sigvalstmp(npols))

c
c       compute all source info at target point on patch
c
      call koorn_pols(uvs,norder,npols,sigvalstmp)

      alpha = 1.0d0
      beta = 0.0d0
      call dgemv('n',9,npols,alpha,srccoefs,9,sigvalstmp,1,beta,
     1      srcvals,1)
      

      do j=1,npols
        do idim=1,nd
          fint(idim,j) = 0
        enddo
      enddo
      ns = 0
      call self_quadrature_new(ier,rad,verts,uvs(1),
     1  uvs(2),srcvals(4),ns,xs,ys,ws)
     
      allocate(srctmp(12,ns),qwts(ns))
      allocate(sigvals(npols,ns))

      do i=1,ns
        uv(1) = xs(i)
        uv(2) = ys(i)

        call koorn_pols(uv,norder,npols,sigvals(1,i))
      enddo

      transa = 'N'
      transb = 'N'
      alpha = 1
      beta = 0
      lda = 9
      ldb = npols
      ldc = 12

      call dgemm(transa,transb,9,ns,npols,alpha,
     1      srccoefs,lda,sigvals,ldb,beta,srctmp,ldc)


c
      do i=1,ns
        call cross_prod3d(srctmp(4,i),srctmp(7,i),srctmp(10,i))
        rr = sqrt(srctmp(10,i)**2+srctmp(11,i)**2+srctmp(12,i)**2)
        qwts(i) = rr*ws(i)
        srctmp(10,i) = srctmp(10,i)/rr
        srctmp(11,i) = srctmp(11,i)/rr
        srctmp(12,i) = srctmp(12,i)/rr

        call fker(nd,srctmp(1,i),ndtarg,targ,ndd,dpars,
     2         ndz,zpars,ndi,ipars,fval)
           
        do j=1,npols
          do idim=1,nd
            fint(idim,j) = fint(idim,j) + fval(idim)*
     1         sigvals(j,i)*qwts(i)
          enddo
        enddo
      enddo

      deallocate(srctmp,qwts,sigvals)

      call zrmatmatt_slow(nd,npols,npols,fint,umat,zquad)


      return
      end
c
c
c
c
c
c
c

      subroutine get_ggq_self_quad_pv_pt(norder,npols,uvs,umat,
     1   srccoefs,ndtarg,targ,rad,lrad,
     1   fker,ndd,dpars,ndz,zpars,ndi,ipars,zquad)
c
c
c       this is the replacement for xtri_zself 
c       for implementing bremer quadrature in ytri. 
c       
c       the main difference is the calling sequence
c       of the subroutine for getting the kernel
c       this is made consistent with adaptive integration
c       on triangles.
c
c
c        The calling sequence for fker is of the form
c          fker(srcinfo,targinfo,
c
c


      implicit real *8(a-h,o-z)
      integer norder,npols,lrad,ndtarg
      real *8 rad(*)
      real *8 srccoefs(9,npols),srcvals(9)
      real *8, allocatable :: sigvalstmp(:)
      real *8 targ(ndtarg)
      real *8 uvs(2),umat(npols,npols)
      complex *16 zquad(*)
      integer ipars(ndi)
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      real *8, allocatable :: srctmp(:,:)
      real *8, allocatable :: qwts(:),sigvals(:,:)
      real *8, allocatable :: xs(:),ys(:),ws(:)
      real *8 uv(2),verts(2,3)
      real *8 alpha,beta
      complex *16 fval
      complex *16, allocatable :: fint(:),fvals(:)
      character *1 transa,transb

      external fker

      nmax = 10000
      allocate(ws(nmax),xs(nmax),ys(nmax))
      allocate(fvals(npols))
      verts(1,1) = 0
      verts(2,1) = 0
      verts(1,2) = 1
      verts(2,2) = 0
      verts(1,3) = 0
      verts(2,3) = 1
      allocate(fint(npols),sigvalstmp(npols))

c
c       compute all source info at target point on patch
c
      call koorn_pols(uvs,norder,npols,sigvalstmp)

      alpha = 1.0d0
      beta = 0.0d0
      call dgemv('n',9,npols,alpha,srccoefs,9,sigvalstmp,1,beta,
     1      srcvals,1)
      

      do j=1,npols
        fint(j) = 0
      enddo
      ns = 0
      call pv_self_quadrature_new(ier,rad,verts,uvs(1),
     1  uvs(2),srcvals(4),ns,xs,ys,ws)
      
     
      allocate(srctmp(12,ns),qwts(ns))
      allocate(sigvals(npols,ns))

      do i=1,ns
        uv(1) = xs(i)
        uv(2) = ys(i)

        call koorn_pols(uv,norder,npols,sigvals(1,i))
      enddo

      transa = 'N'
      transb = 'N'
      alpha = 1
      beta = 0
      lda = 9
      ldb = npols
      ldc = 12

      call dgemm(transa,transb,9,ns,npols,alpha,
     1      srccoefs,lda,sigvals,ldb,beta,srctmp,ldc)


c
      do i=1,ns
        call cross_prod3d(srctmp(4,i),srctmp(7,i),srctmp(10,i))
        rr = sqrt(srctmp(10,i)**2+srctmp(11,i)**2+srctmp(12,i)**2)
        qwts(i) = rr*ws(i)
        srctmp(10,i) = srctmp(10,i)/rr
        srctmp(11,i) = srctmp(11,i)/rr
        srctmp(12,i) = srctmp(12,i)/rr

        call fker(srctmp(1,i),ndtarg,targ,ndd,dpars,
     2         ndz,zpars,ndi,ipars,fval)
           
        do j=1,npols
          fint(j) = fint(j) + fval*sigvals(j,i)*qwts(i)
        enddo
      enddo

      deallocate(srctmp,qwts,sigvals)

      call zrmatmatt_slow(1,npols,npols,fint,umat,zquad)


      return
      end
c
c
c
c
c
c
c
c

      subroutine get_ggq_self_quad_pv_pt_vec(norder,npols,uvs,umat,
     1   srccoefs,ndtarg,targ,rad,lrad,
     1   fker,nd,ndd,dpars,ndz,zpars,ndi,ipars,zquad)
c
c
c       Vectorized version of get_ggq_self_quad_pv_pt
c
c
c        The calling sequence for fker is of the form
c          fker(nd,srcinfo,targinfo,
c
c


      implicit real *8(a-h,o-z)
      integer norder,npols,lrad,ndtarg,nd
      real *8 rad(*)
      real *8 srccoefs(9,npols),srcvals(9)
      real *8, allocatable :: sigvalstmp(:)
      real *8 targ(ndtarg)
      real *8 uvs(2),umat(npols,npols)
      complex *16 zquad(nd,*)
      integer ipars(ndi)
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      real *8, allocatable :: srctmp(:,:)
      real *8, allocatable :: qwts(:),sigvals(:,:)
      real *8, allocatable :: xs(:),ys(:),ws(:)
      real *8 uv(2),verts(2,3)
      real *8 alpha,beta
      complex *16, allocatable :: fint(:,:),fvals(:,:),fval(:)
      integer idim
      character *1 transa,transb

      external fker

      nmax = 3000
      allocate(ws(nmax),xs(nmax),ys(nmax))
      allocate(fvals(nd,npols),fval(nd))
      verts(1,1) = 0
      verts(2,1) = 0
      verts(1,2) = 1
      verts(2,2) = 0
      verts(1,3) = 0
      verts(2,3) = 1
      allocate(fint(nd,npols),sigvalstmp(npols))

c
c       compute all source info at target point on patch
c
      call koorn_pols(uvs,norder,npols,sigvalstmp)

      alpha = 1.0d0
      beta = 0.0d0
      call dgemv('n',9,npols,alpha,srccoefs,9,sigvalstmp,1,beta,
     1      srcvals,1)
      

      do j=1,npols
        do idim=1,nd
          fint(idim,j) = 0
        enddo
      enddo
      ns = 0
      call pv_self_quadrature_new(ier,rad,verts,uvs(1),
     1  uvs(2),srcvals(4),ns,xs,ys,ws)
     
      allocate(srctmp(12,ns),qwts(ns))
      allocate(sigvals(npols,ns))

      do i=1,ns
        uv(1) = xs(i)
        uv(2) = ys(i)

        call koorn_pols(uv,norder,npols,sigvals(1,i))
      enddo

      transa = 'N'
      transb = 'N'
      alpha = 1
      beta = 0
      lda = 9
      ldb = npols
      ldc = 12

      call dgemm(transa,transb,9,ns,npols,alpha,
     1      srccoefs,lda,sigvals,ldb,beta,srctmp,ldc)


c
      do i=1,ns
        call cross_prod3d(srctmp(4,i),srctmp(7,i),srctmp(10,i))
        rr = sqrt(srctmp(10,i)**2+srctmp(11,i)**2+srctmp(12,i)**2)
        qwts(i) = rr*ws(i)
        srctmp(10,i) = srctmp(10,i)/rr
        srctmp(11,i) = srctmp(11,i)/rr
        srctmp(12,i) = srctmp(12,i)/rr

        call fker(nd,srctmp(1,i),ndtarg,targ,ndd,dpars,
     2         ndz,zpars,ndi,ipars,fval)
           
        do j=1,npols
          do idim=1,nd
            fint(idim,j) = fint(idim,j) + fval(idim)*
     1         sigvals(j,i)*qwts(i)
          enddo
        enddo
      enddo

      deallocate(srctmp,qwts,sigvals)

      call zrmatmatt_slow(nd,npols,npols,fint,umat,zquad)


      return
      end
c
