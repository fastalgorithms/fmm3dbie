c
c
c
c
c     This file has the following user callable routines
c
c        ?getnearquad_ggq_guru - guru interface for 
c         computing near field quadrature for a compact/pv kernel
c        
c
c         z - complex
c         d - double precision
c
c  TODO:
c    * Fix behavior of getting quadrature nodes for different
c      patch types for "far-near" part, currently
c      nodes for both triangles are quads are generated. This needs
c      to be changed to figuring out unique sets of patch types
c      and storing the appropirate nodes and weights and then
c      accessing them as needed in loop
c
c

      subroutine zgetnearquad_ggq_guru(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targvals,
     2   ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,
     3   ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear)
c
c------------------------
c  This subroutine generates the near field quadrature
c  for a given kernel which is assumed to be
c  a compact/principal value integral operator
c  where the near field is specified by the user 
c  in row sparse compressed format.
c
c
c  The quadrature is computed by the following strategy
c  targets within a sphere of radius rfac0*rs
c  of a chunk centroid is handled using adaptive integration
c  where rs is the radius of the bounding sphere
c  for the patch
c  
c  All other targets in the near field are handled via
c  oversampled quadrature
c
c
c  Input arguments:
c
c    - npatches: integer
c        number of patches
c    - norders: integer(npatches)
c        order of discretization on each patch
c    - ixyzs: integer(npatches+1)
c        ixyzs(i) denotes the starting location in srccoefs,
c        and srcvals array corresponding to patch i
c    - iptype: integer(npatches)
c        type of patch
c        *  iptype = 1, triangular patch discretized using RV nodes
c        *  iptype = 11, quadrangular patch discretized using GL nodes,
c                        full degree polynomials
c        *  iptype = 12, quadrangular patch discretized using cheb nodes,
c                        full degree polynomials
c    - npts: integer
c        total number of discretization points on the boundary
c    - srccoefs: double precision (9,npts)
c        basis expansion coefficients of xyz, dxyz/du,
c        and dxyz/dv on each patch. For each patch
c        * srccoefs(1:3,i) is xyz info
c        * srccoefs(4:6,i) is dxyz/du info
c        * srccoefs(7:9,i) is dxyz/dv info
c    - srcvals: double precision (12,npts)
c        xyz(u,v) and derivative info sampled at the 
c        discretization nodes on the surface
c        * srcvals(1:3,i) - xyz info
c        * srcvals(4:6,i) - dxyz/du info
c        * srcvals(7:9,i) - dxyz/dv info
c        * srcvals(10:12,i) - normals info
c    - ndtarg: integer
c        leading dimension of target array
c    - ntarg: integer
c        number of targets
c    - targvals: double precision (ndtarg,ntarg)
c        target info. First three components must be x,y,z
c        coordinates
c    - ipatch_id: integer(ntarg)
c        patch id of target, patch_id = -1, if target off-surface
c    - uvs_targ: double precision (2,ntarg)
c        local uv coordinates on patch if target on surface
c    - eps: double precision
c        precision requested
c    - ipv: integer
c        Flag for choosing type of self-quadrature
c        * ipv = 0, for compact/weakly singular operators
c        * ipv = 1, for singular operators
c    - fker: procedure pointer
c        function handle for the kernel. Calling sequence 
c        * fker(srcinfo,ndtarg,targinfo,ndd,dpars,ndz,zpars,ndi,ipars,val)
c    - ndd: integer
c        number of double precision parameters
c    - dpars: double precision(ndd)
c        double precision parameters
c    - ndz: integer
c        number of double complex parameters
c    - zpars: double complex(ndz)
c        double complex parameters
c    - ndi: integer
c        number of integer parameters
c    - ipars: integer(ndi)
c        integer parameters
c    - nnz: integer
c        number of source patch-> target interactions in the near field
c    - row_ptr: integer(ntarg+1)
c        row_ptr(i) is the pointer
c        to col_ind array where list of relevant source patches
c        for target i start
c    - col_ind: integer (nnz)
c        list of source patches relevant for all targets, sorted
c        by the target number
c    - iquad: integer(nnz+1)
c        location in wnear array where quadrature for col_ind(i)
c        starts
c    - rfac0: integer
c        radius parameter for near field
c    - nquad: integer
c        number of entries in wnear
c
c  Output parameters:
c
c    - wnear: double complex(nquad)
c        near field quadrature corrections
c----------------------------------------------------               
c
      implicit real *8 (a-h,o-z)
      integer, intent(in) :: ndi,ndd,ndz,ipv
      integer, intent(in) :: ipars(ndi)
      integer, intent(in) :: ndtarg
      integer, intent(in) :: npatches,norders(npatches),npts
      integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      integer, intent(in) :: ntarg,ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      real *8, intent(in) :: targvals(ndtarg,ntarg)
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      integer, intent(in) :: nnz
      integer, intent(in) :: row_ptr(ntarg+1),col_ind(nnz),iquad(nnz+1)
      complex *16, intent(out) :: wnear(nquad)

      integer ntrimax
      real *8, allocatable :: cms(:,:),rads(:)
      real *8, allocatable :: targ_near(:,:),targ_far(:,:)
      integer, allocatable :: iind_near(:),iind_far(:)
      real *8, allocatable :: umatr(:,:),vmatr(:,:),uvs(:,:),wts(:)

c
c        temporary variables
c
      integer, allocatable :: col_ptr(:),row_ind(:),iper(:)
      complex *16, allocatable :: sints_n(:,:),svtmp_n(:,:)
      complex *16, allocatable :: sints_f(:,:),svtmp_f(:,:)
      complex *16, allocatable :: svtmp2(:,:)

      real *8, allocatable :: xyztarg2(:,:)
      integer irad

      real *8, allocatable :: qnodes_tri(:,:),qwts_tri(:)
      real *8, allocatable :: qnodes_quad(:,:),qwts_quad(:)
      real *8 ra

      complex *16 ff1,ff2,cra1,cra2
      integer nlev, nqorder_f,norder_avg
      real *8 rfac0,rsc,rr,tmp(3),epsp

      integer ipoly
      character *1 ttype

      external fker
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

      npmax = 3000
      norder_avg = floor(sum(norders)/(npatches+0.0d0))
      if(norder_avg.ge.8) npmax = 6000

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
       

c
c    Get triangle parameters
c
      call get_quadparams_adap(eps,1,nqorder,eps_adap,nlev,nqorder_f)

      call triasymq_pts(nqorder_f,nnodes)
      
      ntri_f = 4**nlev
      npts_f_tri = ntri_f*nnodes
      allocate(qnodes_tri(2,npts_f_tri),qwts_tri(npts_f_tri))

      call gen_xg_unif_nodes_tri(nlev,nqorder_f,nnodes,npts_f_tri,
     1   qnodes_tri,qwts_tri)

c
c   Get quad parameters
c
      
      call get_quadparams_adap(eps,11,nqorder,eps_adap,nlev,nqorder_f)

      call squarearbq_pts(nqorder_f,nnodes)
      
      nquad_f = 4**nlev
      npts_f_quad = nquad_f*nnodes
      allocate(qnodes_quad(2,npts_f_quad),qwts_quad(npts_f_quad))

      call gen_xg_unif_nodes_quad(nlev,nqorder_f,nnodes,npts_f_quad,
     1   qnodes_quad,qwts_quad)
      ra = sum(qwts_quad)



      t1 = second()
C$        t1 = omp_get_wtime()

C$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(DYNAMIC)
C$OMP$PRIVATE(ipatch,ntarg2,xyztarg2,sints_n,sints_f,svtmp_n,svtmp_f)
C$OMP$PRIVATE(ii,ii2,i,jpatch,svtmp2,iiif,l,ntarg2m)
C$OMP$PRIVATE(j,iii,istart,itarg2,iqstart,jpt,jtarg)
C$OMP$PRIVATE(targ_near,targ_far,iind_near,iind_far,rr)
C$OMP$PRIVATE(ntarg_f,ntarg_n,npols,norder,irad)
C$OMP$PRIVATE(uvs,umatr,vmatr,wts)
C$OMP$PRIVATE(epsp,rsc,tmp,ipoly,ttype)

      do ipatch=1,npatches
        
        npols = ixyzs(ipatch+1)-ixyzs(ipatch)
        norder = norders(ipatch)
        allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols))
        allocate(wts(npols))
        if(iptype(ipatch).eq.11) ipoly = 0
        if(iptype(ipatch).eq.12) ipoly = 1
        ttype = "F"

c        print *, "ipatch=",ipatch
c        print *, "npols=",npols
c        print *, "iptype(ipatch)=",iptype(ipatch)
        

        call get_disc_exps(norder,npols,iptype(ipatch),uvs,
     1     umatr,vmatr,wts)

        if(norder.le.4) then
          irad = 1
        else if(norder.le.8) then
          irad = 2
        else
          irad = 3
        endif

c
c  estimate rescaling of epsp needed to account for the scale of the
c    patch 
c

        ii = ixyzs(ipatch)
        call cross_prod3d(srcvals(4,ii),srcvals(7,ii),tmp)
        rsc = (tmp(1)**2 + tmp(2)**2 + tmp(3)**2)*wts(1)**2
        do i=2,npols
          call cross_prod3d(srcvals(4,ii+i-1),srcvals(7,ii+i-1),tmp)
          rr = (tmp(1)**2 + tmp(2)**2 + tmp(3)**2)*wts(i)**2
          if(rr.lt.rsc) rsc = rr
        enddo
        epsp = eps_adap*rsc**0.25d0

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
c       fill out near part of layer potential
c

        sints_n = 0
        sints_f = 0
        if(iptype(ipatch).eq.1) 
     1      call ctriaints(epsp,istrat,intype,ntest0,norder,npols,
     1      srccoefs(1,istart),ndtarg,ntarg_n,targ_near,ifp,xyztarg2,
     2      itargptr,ntarg_n,norder,npols,fker,ndd,dpars,ndz,zpars,ndi,
     3      ipars,nqorder,npmax,rfac,sints_n,ifmetric,rn1,n2)
        
        if(iptype(ipatch).eq.11.or.iptype(ipatch).eq.12) 
     1      call cquadints(epsp,istrat,intype,ntest0,norder,ipoly,
     1      ttype,npols,srccoefs(1,istart),ndtarg,ntarg_n,targ_near,
     1      ifp,xyztarg2,itargptr,ntarg_n,norder,npols,fker,ndd,dpars,
     1      ndz,zpars,ndi,ipars,nqorder,npmax,rfac,sints_n,ifmetric,
     1      rn1,n2)
        

        call zrmatmatt(ntarg_n,npols,sints_n,npols,umatr,svtmp_n)
cc
cc       fill out far part of layer potential
c
        if(iptype(ipatch).eq.1) 
     1     call ctriaints_wnodes(ntest0,norder,npols,
     1      srccoefs(1,istart),ndtarg,ntarg_f,targ_far,
     2      itargptr,ntarg_f,norder,npols,fker,ndd,dpars,ndz,zpars,
     3      ndi,ipars,npts_f_tri,qnodes_tri,qwts_tri,sints_f)
        if(iptype(ipatch).eq.11.or.iptype(ipatch).eq.12) 
     1     call cquadints_wnodes(ntest0,norder,ipoly,ttype,npols,
     1      srccoefs(1,istart),ndtarg,ntarg_f,targ_far,
     2      itargptr,ntarg_f,norder,npols,fker,ndd,dpars,ndz,zpars,
     3      ndi,ipars,npts_f_quad,qnodes_quad,qwts_quad,sints_f)
        
        
        call zrmatmatt(ntarg_f,npols,sints_f,npols,umatr,svtmp_f)
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
            call zget_ggq_self_quad_pt(ipv,norder,npols,
     1      uvs_targ(1,jtarg),iptype(ipatch),
     2      umatr,srccoefs(1,istart),ndtarg,
     2      targvals(1,jtarg),irad,fker,ndd,dpars,ndz,zpars,ndi,
     3      ipars,wnear(iqstart+1))
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

      subroutine dgetnearquad_ggq_guru(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targvals,
     2   ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,
     3   ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear)
c
c------------------------
c  This subroutine generates the near field quadrature
c  for a given kernel which is assumed to be
c  a compact/principal value integral operator
c  where the near field is specified by the user 
c  in row sparse compressed format.
c
c
c  The quadrature is computed by the following strategy
c  targets within a sphere of radius rfac0*rs
c  of a chunk centroid is handled using adaptive integration
c  where rs is the radius of the bounding sphere
c  for the patch
c  
c  All other targets in the near field are handled via
c  oversampled quadrature
c
c
c  Input arguments:
c
c    - npatches: integer
c        number of patches
c    - norders: integer(npatches)
c        order of discretization on each patch
c    - ixyzs: integer(npatches+1)
c        ixyzs(i) denotes the starting location in srccoefs,
c        and srcvals array corresponding to patch i
c    - iptype: integer(npatches)
c        type of patch
c        *  iptype = 1, triangular patch discretized using RV nodes
c        *  iptype = 11, quadrangular patch discretized using GL nodes,
c                        full degree polynomials
c        *  iptype = 12, quadrangular patch discretized using cheb nodes,
c                        full degree polynomials
c    - npts: integer
c        total number of discretization points on the boundary
c    - srccoefs: double precision (9,npts)
c        basis expansion coefficients of xyz, dxyz/du,
c        and dxyz/dv on each patch. For each patch
c        * srccoefs(1:3,i) is xyz info
c        * srccoefs(4:6,i) is dxyz/du info
c        * srccoefs(7:9,i) is dxyz/dv info
c    - srcvals: double precision (12,npts)
c        xyz(u,v) and derivative info sampled at the 
c        discretization nodes on the surface
c        * srcvals(1:3,i) - xyz info
c        * srcvals(4:6,i) - dxyz/du info
c        * srcvals(7:9,i) - dxyz/dv info
c        * srcvals(10:12,i) - normals info
c    - ndtarg: integer
c        leading dimension of target array
c    - ntarg: integer
c        number of targets
c    - targvals: double precision (ndtarg,ntarg)
c        target info. First three components must be x,y,z
c        coordinates
c    - ipatch_id: integer(ntarg)
c        patch id of target, patch_id = -1, if target off-surface
c    - uvs_targ: double precision (2,ntarg)
c        local uv coordinates on patch if target on surface
c    - eps: double precision
c        precision requested
c    - ipv: integer
c        Flag for choosing type of self-quadrature
c        * ipv = 0, for compact/weakly singular operators
c        * ipv = 1, for singular operators
c    - fker: procedure pointer
c        function handle for the kernel. Calling sequence 
c        * fker(srcinfo,ndtarg,targinfo,ndd,dpars,ndz,zpars,ndi,ipars,val)
c    - ndd: integer
c        number of double precision parameters
c    - dpars: double precision(ndd)
c        double precision parameters
c    - ndz: integer
c        number of double complex parameters
c    - zpars: double complex(ndz)
c        double complex parameters
c    - ndi: integer
c        number of integer parameters
c    - ipars: integer(ndi)
c        integer parameters
c    - nnz: integer
c        number of source patch-> target interactions in the near field
c    - row_ptr: integer(ntarg+1)
c        row_ptr(i) is the pointer
c        to col_ind array where list of relevant source patches
c        for target i start
c    - col_ind: integer (nnz)
c        list of source patches relevant for all targets, sorted
c        by the target number
c    - iquad: integer(nnz+1)
c        location in wnear array where quadrature for col_ind(i)
c        starts
c    - rfac0: integer
c        radius parameter for near field
c    - nquad: integer
c        number of entries in wnear
c
c  Output parameters:
c
c    - wnear: double precision(nquad)
c        near field quadrature corrections
c----------------------------------------------------               
c
      implicit real *8 (a-h,o-z)
      integer, intent(in) :: ndi,ndd,ndz,ipv
      integer, intent(in) :: ipars(ndi)
      integer, intent(in) :: ndtarg
      integer, intent(in) :: npatches,norders(npatches),npts
      integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      integer, intent(in) :: ntarg,ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      real *8, intent(in) :: targvals(ndtarg,ntarg)
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      integer, intent(in) :: nnz
      integer, intent(in) :: row_ptr(ntarg+1),col_ind(nnz),iquad(nnz+1)
      real *8, intent(out) :: wnear(nquad)

      integer ntrimax
      real *8, allocatable :: cms(:,:),rads(:)
      real *8, allocatable :: targ_near(:,:),targ_far(:,:)
      integer, allocatable :: iind_near(:),iind_far(:)
      real *8, allocatable :: umatr(:,:),vmatr(:,:),uvs(:,:),wts(:)

c
c        temporary variables
c
      integer, allocatable :: col_ptr(:),row_ind(:),iper(:)
      real *8, allocatable :: sints_n(:,:),svtmp_n(:,:)
      real *8, allocatable :: sints_f(:,:),svtmp_f(:,:)
      real *8, allocatable :: svtmp2(:,:)

      real *8, allocatable :: xyztarg2(:,:)
      integer irad

      real *8, allocatable :: qnodes_tri(:,:),qwts_tri(:)
      real *8, allocatable :: qnodes_quad(:,:),qwts_quad(:)
      real *8 ra

      real *8 ff1,ff2,cra1,cra2
      integer nlev, nqorder_f
      real *8 rfac0,rsc,rr,tmp(3),epsp
      real *8 done,dzero
      integer norder_avg
      integer ipoly
      character *1 ttype


      external fker

      done = 1
      dzero = 0
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

      npmax = 3000
      norder_avg = floor(sum(norders)/(npatches+0.0d0))
      if(norder_avg.ge.8) npmax = 6000

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
c
c    Get triangle parameters
c
      call get_quadparams_adap(eps,1,nqorder,eps_adap,nlev,nqorder_f)

      call triasymq_pts(nqorder_f,nnodes)
      
      ntri_f = 4**nlev
      npts_f_tri = ntri_f*nnodes
      allocate(qnodes_tri(2,npts_f_tri),qwts_tri(npts_f_tri))

      call gen_xg_unif_nodes_tri(nlev,nqorder_f,nnodes,npts_f_tri,
     1   qnodes_tri,qwts_tri)

c
c   Get quad parameters
c
      
      call get_quadparams_adap(eps,11,nqorder,eps_adap,nlev,nqorder_f)

      call squarearbq_pts(nqorder_f,nnodes)
      
      nquad_f = 4**nlev
      npts_f_quad = nquad_f*nnodes
      allocate(qnodes_quad(2,npts_f_quad),qwts_quad(npts_f_quad))

      call gen_xg_unif_nodes_quad(nlev,nqorder_f,nnodes,npts_f_quad,
     1   qnodes_quad,qwts_quad)


      t1 = second()
C$        t1 = omp_get_wtime()

C$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(DYNAMIC)
C$OMP$PRIVATE(ipatch,ntarg2,xyztarg2,sints_n,sints_f,svtmp_n,svtmp_f)
C$OMP$PRIVATE(ii,ii2,i,jpatch,svtmp2,iiif,l,ntarg2m)
C$OMP$PRIVATE(j,iii,istart,itarg2,iqstart,jpt,jtarg)
C$OMP$PRIVATE(targ_near,targ_far,iind_near,iind_far,rr)
C$OMP$PRIVATE(ntarg_f,ntarg_n,npols,norder,irad)
C$OMP$PRIVATE(uvs,umatr,vmatr,wts)
C$OMP$PRIVATE(rsc,tmp,epsp,ipoly,ttype)

      do ipatch=1,npatches
        
        npols = ixyzs(ipatch+1)-ixyzs(ipatch)
        norder = norders(ipatch)

        if(norder.le.4) then
          irad = 1
        else if(norder.le.8) then
          irad = 2
        else
          irad = 3
        endif

        allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols))
        allocate(wts(npols))

        call get_disc_exps(norder,npols,iptype(ipatch),uvs,umatr,
     1     vmatr,wts)

c
c  estimate rescaling of epsp needed to account for the scale of the
c    patch 
c

        ii = ixyzs(ipatch)
        call cross_prod3d(srcvals(4,ii),srcvals(7,ii),tmp)
        rsc = (tmp(1)**2 + tmp(2)**2 + tmp(3)**2)*wts(1)**2
        do i=2,npols
          call cross_prod3d(srcvals(4,ii+i-1),srcvals(7,ii+i-1),tmp)
          rr = (tmp(1)**2 + tmp(2)**2 + tmp(3)**2)*wts(i)**2
          if(rr.lt.rsc) rsc = rr
        enddo
        epsp = eps_adap*rsc**0.25d0

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
c       fill out near part of layer potential
c

        if(iptype(ipatch).eq.1) 
     1      call dtriaints(epsp,istrat,intype,ntest0,norder,npols,
     1      srccoefs(1,istart),ndtarg,ntarg_n,targ_near,ifp,xyztarg2,
     2      itargptr,ntarg_n,norder,npols,fker,ndd,dpars,ndz,zpars,ndi,
     3      ipars,nqorder,npmax,rfac,sints_n,ifmetric,rn1,n2)
        
        if(iptype(ipatch).eq.11.or.iptype(ipatch).eq.12) 
     1      call dquadints(epsp,istrat,intype,ntest0,norder,npols,
     1      ipoly,ttype,srccoefs(1,istart),ndtarg,ntarg_n,targ_near,
     1      ifp,xyztarg2,itargptr,ntarg_n,norder,npols,fker,ndd,dpars,
     1      ndz,zpars,ndi,ipars,nqorder,npmax,rfac,sints_n,ifmetric,
     1      rn1,n2)
        

        call dgemm_guru('t','n',npols,ntarg_n,npols,done,umatr,npols,
     1     sints_n,npols,dzero,svtmp_n,npols)
c
c       fill out far part of layer potential
c
        if(iptype(ipatch).eq.1) 
     1     call dtriaints_wnodes(ntest0,norder,npols,
     1      srccoefs(1,istart),ndtarg,ntarg_f,targ_far,
     2      itargptr,ntarg_f,norder,npols,fker,ndd,dpars,ndz,zpars,
     3      ndi,ipars,npts_f_tri,qnodes_tri,qwts_tri,sints_f)
        if(iptype(ipatch).eq.11.or.iptype(ipatch).eq.12) 
     1     call dquadints_wnodes(ntest0,norder,ipoly,ttype,npols,
     1      srccoefs(1,istart),ndtarg,ntarg_f,targ_far,
     2      itargptr,ntarg_f,norder,npols,fker,ndd,dpars,ndz,zpars,
     3      ndi,ipars,npts_f_quad,qnodes_quad,qwts_quad,sints_f)
        
        call dgemm_guru('t','n',npols,ntarg_f,npols,done,umatr,npols,
     1     sints_f,npols,dzero,svtmp_f,npols)
c
c       combine svtmp_f, svtmp_n to fill out svtmp2
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
            call dget_ggq_self_quad_pt(ipv,norder,npols,
     1      uvs_targ(1,jtarg),iptype(ipatch),umatr,srccoefs(1,istart),
     2      ndtarg,targvals(1,jtarg),irad,fker,ndd,dpars,ndz,zpars,ndi,
     3      ipars,wnear(iqstart+1))
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

      subroutine zget_ggq_self_quad_pt(ipv,norder,npols,uvs,iptype,umat,
     1   srccoefs,ndtarg,targ,irad,fker,ndd,dpars,ndz,zpars,ndi,
     2   ipars,zquad)
c
c
c------------------
c  This subroutine evaluates the integral of an integral
c  operator at a target on the interior of a triangluar patch 
c  using the generalized
c  Gaussian quadratures developed by Bremer and Gimbutas.
c
c  The quadrature currently cannot handle targets on the boundary
c  of the triangle
c
c  Input arguments:
c    - ipv: integer
c        Flag for choosing type of self-quadrature
c        * ipv = 0, for compact/weakly singular operators
c        * ipv = 1, for singular operators
c    - norder: integer
c        order of patch discretization
c    - npols: integer
c        number of discretization nodes on patch
c    - uvs: double precision(2)
c        local u,v coordinates of target
c    - iptype: integer(npatches)
c        type of patch
c        *  iptype = 1, triangular patch discretized using RV nodes
c        *  iptype = 11, quadrangular patch discretized using GL nodes,
c                        full degree polynomials
c        *  iptype = 12, quadrangular patch discretized using cheb nodes,
c                        full degree polynomials
c    - umat: double precision(npols,npols)
c        values to coefficient matrix
c    - srccoefs: double precision (9,npts)
c        basis expansion coefficients of xyz, dxyz/du,
c        and dxyz/dv on each patch. For each patch
c        * srccoefs(1:3,i) is xyz info
c        * srccoefs(4:6,i) is dxyz/du info
c        * srccoefs(7:9,i) is dxyz/dv info
c    - ndtarg: integer
c        leading dimension of target array
c    - targ: double precision (ndtarg)
c        target info. First three components must be x,y,z
c        coordinates
c    - irad: integer
c        which order of self quadrature to use
c        * for norder<=4, irad=1
c        * for 4<norder<=8, irad=2
c        * for 8<norder<=12, irad=3
c        * irad = 3, for higher orders, but there might be some loss
c          in precision
c    - fker: procedure pointer
c        function handle for the kernel. Calling sequence 
c        * fker(srcinfo,ndtarg,targinfo,ndd,dpars,ndz,zpars,ndi,ipars,val)
c    - ndd: integer
c        number of double precision parameters
c    - dpars: double precision(ndd)
c        double precision parameters
c    - ndz: integer
c        number of double complex parameters
c    - zpars: double complex(ndz)
c        double complex parameters
c    - ndi: integer
c        number of integer parameters
c    - ipars: integer(ndi)
c        integer parameters
c
c  Output parameters:
c    
c    - zquad: double complex(npols)
c        near quadrature correction
c    
c
c


      implicit real *8(a-h,o-z)
      integer, intent(in) :: norder,npols,ndtarg
      real *8, intent(in) :: srccoefs(9,npols)
      real *8, intent(in) :: targ(ndtarg)
      real *8, intent(in) :: uvs(2),umat(npols,npols)
      integer, intent(in) :: ndi,ndd,ndz
      integer, intent(in) :: ipars(ndi)
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      complex *16, intent(out) :: zquad(npols)

      real *8 srcvals(9)
      real *8, allocatable :: sigvalstmp(:)
      real *8, allocatable :: srctmp(:,:)
      real *8, allocatable :: qwts(:),sigvals(:,:)
      real *8, allocatable :: xs(:),ys(:),ws(:)
      real *8 uv(2),verts(2,4)
      real *8 alpha,beta
      complex *16 fval
      complex *16, allocatable :: fint(:),fvals(:)
      character *1 transa,transb

      external fker

      nmax = 20000
      allocate(ws(nmax),xs(nmax),ys(nmax))
      allocate(fvals(npols))
      if(iptype.eq.1) then
        nv = 3
        verts(1,1) = 0
        verts(2,1) = 0
        verts(1,2) = 1
        verts(2,2) = 0
        verts(1,3) = 0
        verts(2,3) = 1
      endif
      if(iptype.eq.11.or.iptype.eq.12) then
         nv = 4
         verts(1,1) = -1
         verts(2,1) = -1

         verts(1,2) = 1
         verts(2,2) = -1

         verts(1,3) = 1
         verts(2,3) = 1
         
         verts(1,4) = -1
         verts(2,4) = 1
      endif
      allocate(fint(npols),sigvalstmp(npols))

c
c       compute all source info at target point on patch
c
      call get_basis_pols(uvs,norder,npols,iptype,sigvalstmp)

      alpha = 1.0d0
      beta = 0.0d0
      call dgemv_guru('n',9,npols,alpha,srccoefs,9,sigvalstmp,1,beta,
     1      srcvals,1)
      

      do j=1,npols
        fint(j) = 0
      enddo
      ns = 0
      if(ipv.eq.0) then
        call self_quadrature_new(irad,verts,nv,uvs(1),
     1    uvs(2),srcvals(4),ns,xs,ys,ws)
      endif

      if(ipv.eq.1) then
        call pv_self_quadrature_new(irad,verts,nv,uvs(1),
     1    uvs(2),srcvals(4),ns,xs,ys,ws)
      endif
     
      allocate(srctmp(12,ns),qwts(ns))
      allocate(sigvals(npols,ns))

      do i=1,ns
        uv(1) = xs(i)
        uv(2) = ys(i)
        call get_basis_pols(uv,norder,npols,iptype,sigvals(1,i))
      enddo

      transa = 'N'
      transb = 'N'
      alpha = 1
      beta = 0
      lda = 9
      ldb = npols
      ldc = 12

      call dgemm_guru(transa,transb,9,ns,npols,alpha,
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

      call zrmatmatt(1,npols,fint,npols,umat,zquad)

      return
      end
c
c
c
c
c
c

      subroutine dget_ggq_self_quad_pt(ipv,norder,npols,uvs,iptype,
     1   umat,srccoefs,ndtarg,targ,irad,fker,ndd,dpars,ndz,zpars,ndi,
     2   ipars,dquad)
c
c
c------------------
c  This subroutine evaluates the integral of an integral
c  operator at a target on the interior of a triangluar patch 
c  using the generalized
c  Gaussian quadratures developed by Bremer and Gimbutas.
c
c  The quadrature currently cannot handle targets on the boundary
c  of the triangle
c
c  Input arguments:
c    - ipv: integer
c        Flag for choosing type of self-quadrature
c        * ipv = 0, for compact/weakly singular operators
c        * ipv = 1, for singular operators
c    - norder: integer
c        order of patch discretization
c    - npols: integer
c        number of discretization nodes on patch
c    - uvs: double precision(2)
c        local u,v coordinates of target
c    - iptype: integer(npatches)
c        type of patch
c        *  iptype = 1, triangular patch discretized using RV nodes
c        *  iptype = 11, quadrangular patch discretized using GL nodes,
c                        full degree polynomials
c        *  iptype = 12, quadrangular patch discretized using cheb nodes,
c                        full degree polynomials
c    - umat: double precision(npols,npols)
c        values to coefficient matrix
c    - srccoefs: double precision (9,npts)
c        basis expansion coefficients of xyz, dxyz/du,
c        and dxyz/dv on each patch. For each patch
c        * srccoefs(1:3,i) is xyz info
c        * srccoefs(4:6,i) is dxyz/du info
c        * srccoefs(7:9,i) is dxyz/dv info
c    - ndtarg: integer
c        leading dimension of target array
c    - targ: double precision (ndtarg)
c        target info. First three components must be x,y,z
c        coordinates
c    - irad: integer
c        which order of self quadrature to use
c        * for norder<=4, irad=1
c        * for 4<norder<=8, irad=2
c        * for 8<norder<=12, irad=3
c        * irad = 3, for higher orders, but there might be some loss
c          in precision
c    - fker: procedure pointer
c        function handle for the kernel. Calling sequence 
c        * fker(srcinfo,ndtarg,targinfo,ndd,dpars,ndz,zpars,ndi,ipars,val)
c    - ndd: integer
c        number of double precision parameters
c    - dpars: double precision(ndd)
c        double precision parameters
c    - ndz: integer
c        number of double complex parameters
c    - zpars: double complex(ndz)
c        double complex parameters
c    - ndi: integer
c        number of integer parameters
c    - ipars: integer(ndi)
c        integer parameters
c
c  Output parameters:
c    
c    - dquad: double precision(npols)
c        near quadrature correction
c    
c
c


      implicit real *8(a-h,o-z)
      integer, intent(in) :: norder,npols,ndtarg
      real *8, intent(in) :: srccoefs(9,npols)
      real *8, intent(in) :: targ(ndtarg)
      real *8, intent(in) :: uvs(2),umat(npols,npols)
      integer, intent(in) :: ndi,ndd,ndz
      integer, intent(in) :: ipars(ndi)
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      real *8, intent(out) :: dquad(npols)

      real *8 srcvals(9)
      real *8, allocatable :: sigvalstmp(:)
      real *8, allocatable :: srctmp(:,:)
      real *8, allocatable :: qwts(:),sigvals(:,:)
      real *8, allocatable :: xs(:),ys(:),ws(:)
      real *8 uv(2),verts(2,4)
      real *8 alpha,beta
      real *8 fval
      real *8, allocatable :: fint(:),fvals(:)
      character *1 transa,transb
      real *8 done,dzero 

      external fker

      done = 1
      dzero = 0

      nmax = 10000
      allocate(ws(nmax),xs(nmax),ys(nmax))
      allocate(fvals(npols))
      if(iptype.eq.1) then
        nv = 3
        verts(1,1) = 0
        verts(2,1) = 0
        verts(1,2) = 1
        verts(2,2) = 0
        verts(1,3) = 0
        verts(2,3) = 1
      endif
      if(iptype.eq.11.or.iptype.eq.12) then
         nv = 4
         verts(1,1) = -1
         verts(2,1) = -1

         verts(1,2) = 1
         verts(2,2) = -1

         verts(1,3) = 1
         verts(2,3) = 1
         
         verts(1,4) = -1
         verts(2,4) = 1
      endif
      allocate(fint(npols),sigvalstmp(npols))

c
c       compute all source info at target point on patch
c
      call get_basis_pols(uvs,norder,npols,iptype,sigvalstmp)

      alpha = 1.0d0
      beta = 0.0d0
      call dgemv_guru('n',9,npols,alpha,srccoefs,9,sigvalstmp,1,beta,
     1      srcvals,1)
      

      do j=1,npols
        fint(j) = 0
      enddo
      ns = 0
      if(ipv.eq.0) then
        call self_quadrature_new(irad,verts,nv,uvs(1),
     1    uvs(2),srcvals(4),ns,xs,ys,ws)
      endif

      if(ipv.eq.1) then
        call pv_self_quadrature_new(irad,verts,nv,uvs(1),
     1    uvs(2),srcvals(4),ns,xs,ys,ws)
      endif
     
      allocate(srctmp(12,ns),qwts(ns))
      allocate(sigvals(npols,ns))

      do i=1,ns
        uv(1) = xs(i)
        uv(2) = ys(i)
        call get_basis_pols(uv,norder,npols,iptype,sigvals(1,i))
      enddo

      transa = 'N'
      transb = 'N'
      alpha = 1
      beta = 0
      lda = 9
      ldb = npols
      ldc = 12

      call dgemm_guru(transa,transb,9,ns,npols,alpha,
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

      call dgemv_guru('t',npols,npols,done,umat,npols,fint,1,dzero,
     1   dquad,1)

      return
      end
c
c
c
c
c
c
