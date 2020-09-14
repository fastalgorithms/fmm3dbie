c
c
c     This file contains the following user callable
c     routines: 
c 
c       getnearquad_stok_s_trac - generates the near
c        field quadrature for the traction of
c        the single layer
c        representation (targ arrays must have normal
c        vector info)
c
c       lpcomp_stok_s_trac 
c          simpler version of Stokes single layer traction
c          potential evaluator
c          only geometry, targets,
c          and density sampled at discretization required on input,
c          output is the layer potential evaluated at the target points
c          Because this routine evaluates the traction, the targets
c          array must have the source array format, in particular
c          normal information should be in targs(10:12,i)
c
c       stok_s_trac_solver - solves the interior/exterior traction
c         problem for the Stokes equation using the single layer
c         representation
c
c    Advanced user interfaces: 
c*****************************
c     Note for developers: One of the two advanced user interfaces
c     is necessary for the easy user interface and
c     efficient iterative solvers. It seems that the add and subtract
c     version tends to have better CPU-time performance, but we expect
c     the setsub version to be more numerically stable
c**************************************
c       lpcomp_stok_s_trac_addsub 
c         compute traction of the layer potential
c         corresponding to the combined field representation
c         using add and subtract
c 
c
c
c



      subroutine getnearquad_stok_s_trac(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2   ipatch_id,uvs_targ,eps,iquadtype,nnz,row_ptr,col_ind,
     3   iquad,rfac0,nquad,wnear)
c
c       this subroutine generates the near field quadrature
c       for the *traction* of the representation u = 4\pi S
c       where S is the Stokes single layer potential.
c       the near field is specified by the user 
c       in row sparse compressed format.
c
c
c        Note: the 4 \pi scaling is included to be consistent with the FMM
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
c       The recommended parameter for rfac0 is 1.25d0
c
c       Note that all the parameters are real and the quadrature
c       correction returned is real. Currently implemented in
c       a hacky manner by calling the complex routine and then
c       converting to real routines
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
c            starting location of data on patch i
c  
c         iptype - integer(npatches)
c           type of patch
c           iptype = 1 -> triangular patch discretized with RV nodes
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
c             srcvals(10:12,i) - normals info
c 
c         ndtarg - integer
c            leading dimension of target array
c        
c         ntarg - integer
c            number of targets
c
c         targs - real *8 (ndtarg,ntarg)
c            target information
c
c         ipatch_id - integer(ntarg)
c            id of patch of target i, id = -1, if target is off-surface
c
c         uvs_targ - real *8 (2,ntarg)
c            local uv coordinates on patch if on surface, otherwise
c            set to 0 by default
c            
c          eps - real *8
c             precision requested
c
c           iquadtype - integer
c              quadrature type
c              iquadtype = 1, use ggq for self + adaptive integration
c                 for rest
c 
c
c           nnz - integer
c             number of source patch-> target interactions in the near field
c 
c           row_ptr - integer(ntarg+1)
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
c               number of near field entries corresponding to
c               each source-target pair 
c
c        output
c            wnear - real *8(6,nquad)
c               the desired near field quadrature
c               
c

      implicit none 
      integer, intent(in) :: npatches,norders(npatches),npts,nquad
      integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      real *8, intent(in) :: rfac0
      integer, intent(in) :: ndtarg,ntarg
      integer, intent(in) :: iquadtype
      real *8, intent(in) :: targs(ndtarg,ntarg)
      integer, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      integer, intent(in) :: nnz
      integer, intent(in) :: row_ptr(ntarg+1),col_ind(nnz),iquad(nnz+1)
      real *8, intent(out) :: wnear(6,nquad)

      real *8, allocatable :: wnear1(:)

      integer ipars(2), ijloc(2,6)
      integer ndd,ndz,ndi
      complex *16 zpars

      real *8 dpars
      integer i,j,ii,l
      integer ipv

      procedure (), pointer :: fker
      external st3d_strac

c
c
c        initialize the appropriate kernel function
c

      allocate(wnear1(nquad))

      ndd = 0
      ndi = 2
      ndz = 0

      ijloc(1,1) = 1
      ijloc(2,1) = 1
      ijloc(1,2) = 1
      ijloc(2,2) = 2
      ijloc(1,3) = 1
      ijloc(2,3) = 3
      ijloc(1,4) = 2
      ijloc(2,4) = 2
      ijloc(1,5) = 2
      ijloc(2,5) = 3
      ijloc(1,6) = 3
      ijloc(2,6) = 3
      
      do l = 1,6
         i = ijloc(1,l)
         j = ijloc(2,l)
         ipars(1) = i
         ipars(2) = j
         fker => st3d_strac
         ipv = 1
         
         call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     1        ipatch_id,uvs_targ,
     1        eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,
     1        row_ptr,col_ind,iquad,rfac0,nquad,wnear1)
         
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ii)       
         do ii=1,nquad
            wnear(l,ii) = wnear1(ii)
         enddo
C$OMP END PARALLEL DO        
      enddo

      return
      end
c
c
c
c
c
      subroutine lpcomp_stok_s_trac(npatches,norders,ixyzs,
     1   iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2   ipatch_id,uvs_targ,eps,sigma,pot)
c
cf2py intent(in) npatches,norders,ixyzs,iptype,npts,srccoefs,srcvals
cf2py intent(in) ndtarg,ntarg,targs,ipatch_id,uvs_targ,eps
cf2py intent(in) sigma,ndpot
cf2py intent(out) pot
c
c
c------------------------------
c     This subroutine evaluates the traction of the layer potential for
c     the representation 
c
c
c  .. math ::
c  
c      u = 4 \pi \mathcal{S} 
c
c  Note: the $4\pi$ scaling is included to be consistent with the FMM3D
c     libraries.
c
c     NOTE: TARGETS MUST BE ON BOUDNARY, i.e. targs array must be
c     length 12 with surface normal information in targs(10:12,i)
c     This routine only computes the principal value part, the identity
c     term corresponding to the jump in the layer potential is not included
c     in the layer potential.
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
c        and srcvals array where information for patch i begins
c    - iptype: integer(npatches)
c        type of patch
c    - npts: integer
c        total number of discretization points on the boundary
c    - srccoefs: double precision (9,npts)
c        koornwinder expansion coefficients of x, $\partial_{u} x$,
c        and $\partial_{v} x$. 
c    - srcvals: double precision (12,npts)
c        x, $\partial_{u} x$, $\partial_{v} x$, and $n$ sampled at
c        discretization nodes
c    - ndtarg: integer
c        leading dimension of target array
c    - ntarg: integer
c        number of targets
c    - targs: double precision (ndtarg,ntarg)
c        target information
c    - ipatch_id: integer(ntarg)
c        id of patch of target i, id = -1, if target is off-surface
c    - uvs_targ: double precision (2,ntarg)
c        local uv coordinates on patch if on surface, otherwise
c        set to 0 by default
c    - eps: double precision
c        precision requested
c     - sigma: double precision(ndsigma,npts)
c        density for layer potential
c
c  Output arguments
c    - pot: double precision(ndpot,ntarg)
c        traction of layer potential evaluated at the target points
c
c-----------------------------------
c
      implicit none
      integer, intent(in) :: npatches,npts
      integer, intent(in) :: ndtarg,ntarg
      integer, intent(in) :: norders(npatches),ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      real *8, intent(in) :: targs(ndtarg,ntarg)
      real *8, intent(in) :: sigma(3,npts)
      integer, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)

      real *8, intent(out) :: pot(3,ntarg)


      integer nptso,nnz,nquad


      integer nover,npolso
      integer norder,npols
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      real *8, allocatable :: wnear(:,:)

      real *8, allocatable :: srcover(:,:),wover(:)
      integer, allocatable :: ixyzso(:),novers(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)

      integer i,j,jpatch,jquadstart,jstart

      complex *16 zpars
      integer ipars
      real *8 timeinfo(10),t1,t2,omp_get_wtime


      real *8 ttot,done,pi
      real *8 rfac,rfac0
      integer iptype_avg,norder_avg
      integer ikerorder, iquadtype,npts_over

c
c
c        this might need fixing
c
      iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
      norder_avg = floor(sum(norders)/(npatches+0.0d0))

      call get_rfacs(norder_avg,iptype_avg,rfac,rfac0)


      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))
      
      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts, 
     1     srccoefs,cms,rads)

C$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,npatches
        rad_near(i) = rads(i)*rfac
      enddo
C$OMP END PARALLEL DO


c     
c    find near quadrature correction interactions
c
      call findnearmem(cms,npatches,rad_near,ndtarg,targs,ntarg,nnz)

      allocate(row_ptr(ntarg+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,ndtarg,targs,ntarg,row_ptr, 
     1        col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,ntarg,nnz,row_ptr,col_ind,
     1         iquad)

      
      ikerorder = 0

c
c    estimate oversampling for far-field, and oversample geometry
c

      allocate(novers(npatches),ixyzso(npatches+1))

      zpars = 0
      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,
     1    rads,npts,srccoefs,ndtarg,ntarg,targs,ikerorder,zpars,
     2    nnz,row_ptr,col_ind,rfac,novers,ixyzso)

      npts_over = ixyzso(npatches+1)-1

      allocate(srcover(12,npts_over),wover(npts_over))

      call oversample_geom(npatches,norders,ixyzs,iptype,npts, 
     1   srccoefs,srcvals,novers,ixyzso,npts_over,srcover)

      call get_qwts(npatches,novers,ixyzso,iptype,npts_over,
     1        srcover,wover)


c
c   compute near quadrature correction
c

c     dimension of kernel (3x3)
      
      nquad = iquad(nnz+1)-1
      allocate(wnear(6,nquad))
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j)   
      do i=1,nquad
         do j = 1,6
            wnear(j,i) = 0
         enddo
      enddo
C$OMP END PARALLEL DO    


      iquadtype = 1

      call getnearquad_stok_s_trac(npatches,norders,
     1     ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     1     ipatch_id,uvs_targ,eps,iquadtype,nnz,row_ptr,col_ind,
     1     iquad,rfac0,nquad,wnear)
      

c
c     
c     compute layer potential
c     
      call lpcomp_stok_s_trac_addsub(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2     eps,nnz,row_ptr,col_ind,iquad,nquad,wnear,
     3     sigma,novers,npts_over,ixyzso,srcover,wover,pot)
      


      return
      end
c
c
c
c
c
      subroutine lpcomp_stok_s_trac_addsub(npatches,norders,ixyzs,
     1   iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2     eps,nnz,row_ptr,col_ind,iquad,nquad,wnear,sigma,
     3     novers,nptso,ixyzso,srcover,whtsover,pot)
c
c
c      this subroutine evaluates the layer potential for
c      the traction of the representation u = 4\pi S 
c      where S is the single layer potential for the Stokes operator.
c      the near field is precomputed and stored
c      in the row sparse compressed format.
c
c       Note the 4\pi scaling is included to be consistent with the FMM
c
c
c     The fmm is used to accelerate the far-field and 
c     near-field interactions are handled via precomputed quadrature
c
c
c     Using add and subtract - no need to call tree and set fmm parameters
c      can directly call existing fmm library
c
c
c       input:
c         npatches - integer
c            number of patches
c
c         norders- integer(npatches)
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
c         srcvals - real *8 (12,npts)
c             xyz(u,v) and derivative info sampled at the 
c             discretization nodes on the surface
c             srcvals(1:3,i) - xyz info
c             srcvals(4:6,i) - dxyz/du info
c             srcvals(7:9,i) - dxyz/dv info
c             srcvals(10:12,i) - normals info
c 
c         ndtarg - integer
c            leading dimension of target array
c        
c         ntarg - integer
c            number of targets
c
c         targs - real *8 (ndtarg,ntarg)
c            target information
c
c          eps - real *8
c             precision requested
c
c           nnz - integer *8
c             number of source patch-> target interactions in the near field
c 
c           row_ptr - integer(ntarg+1)
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
c           nquad - integer
c               number of entries in wnear
c
c           wnear - real *8(6,nquad)
c               the near field quadrature correction
c
c           sigma - real *8(3,npts)
c               density for layer potential
c
c           novers - integer(npatches)
c              order of discretization for oversampled sources and
c               density
c
c         ixyzso - integer(npatches+1)
c            ixyzso(i) denotes the starting location in srcover,
c               corresponding to patch i
c   
c           nptso - integer
c              total number of oversampled points
c
c           srcover - real *8 (12,nptso)
c              oversampled set of source information
c
c           whtsover - real *8 (nptso)
c             smooth quadrature weights at oversampled nodes
c
c
c         output
c           pot - real *8(3,npts)
c              layer potential evaluated at the target points
c
c           
c               
c
      implicit none
      integer, intent(in) :: npatches,npts
      integer, intent(in) :: ndtarg,ntarg
      integer, intent(in) :: norders(npatches),ixyzs(npatches+1)
      integer, intent(in) :: ixyzso(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      real *8, intent(in) :: targs(ndtarg,ntarg)
      integer, intent(in) :: nnz,row_ptr(ntarg+1),col_ind(nnz),nquad
      integer, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: wnear(6,nquad),sigma(3,npts)
      integer, intent(in) :: novers(npatches+1)
      integer, intent(in) :: nptso
      real *8, intent(in) :: srcover(12,nptso),whtsover(nptso)
      real *8, intent(out) :: pot(3,ntarg)

      integer norder,npols,nover,npolso,ndsigma
      real *8, allocatable :: potsort(:)

      real *8, allocatable :: sources(:,:),targvals(:,:)
      real *8, allocatable :: sigmaover(:,:), stoklet(:,:)
      real *8, allocatable :: pottarg(:,:), pretarg(:)
      real *8, allocatable :: gradtarg(:,:,:)
      integer ns,nt

      integer ifstoklet,ifstrslet
      integer ifppreg,ifppregtarg
      real *8 tmp(10),val

      real *8 xmin,xmax,ymin,ymax,zmin,zmax,sizey,sizez,boxsize


      integer i,j,jpatch,jquadstart,jstart


      integer ifaddsub

      integer ntj
      
      real *8 ddot,pottmp
      real *8, allocatable :: sttmp2(:,:)
      real *8 radexp,epsfmm

      integer ipars
      complex *16 zpars
      real *8 timeinfo(10),t1,t2,omp_get_wtime

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      real *8 thresh
      real *8 rr,rmin
      integer nss,ii,l,npover

      integer ntarg0, nd

      real *8 ttot,done,pi

      real *8 e11,e12,e13,e21,e22,e23,e31,e32,e33,p,dn1,dn2,dn3
      real *8 w11,w12,w13,w21,w22,w23,w31,w32,w33,sig1,sig2,sig3
      real *8 tmp3(3), gradv(3,3)

      parameter (ntarg0=1)

      ns = nptso
      done = 1
      pi = atan(done)*4

           
      allocate(sources(3,ns),targvals(3,ntarg))
      allocate(pottarg(3,ntarg),pretarg(ntarg),gradtarg(3,3,ntarg))

      allocate(sigmaover(3,ns),stoklet(3,ns))

c 
c       oversample density
c

      ndsigma = 3
      
      call oversample_fun_surf(ndsigma,npatches,norders,ixyzs,iptype, 
     1    npts,sigma,novers,ixyzso,ns,sigmaover)




c
c       set relevatn parameters for the fmm
c

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)      
      do i=1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)

        stoklet(1,i) = sigmaover(1,i)*whtsover(i)
        stoklet(2,i) = sigmaover(2,i)*whtsover(i)
        stoklet(3,i) = sigmaover(3,i)*whtsover(i)
      enddo
C$OMP END PARALLEL DO      

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,ntarg
        targvals(1,i) = targs(1,i)
        targvals(2,i) = targs(2,i)
        targvals(3,i) = targs(3,i)
      enddo
C$OMP END PARALLEL DO      

      ifstoklet = 1
      ifstrslet = 0

c
c       call the fmm
c

      nd = 1
      ifppreg = 0
      ifppregtarg = 3
      call cpu_time(t1)
C$      t1 = omp_get_wtime()      
      call stfmm3d(nd,eps,ns,sources,ifstoklet,stoklet,
     1     ifstrslet,tmp,tmp,ifppreg,tmp,tmp,tmp,
     2     ntarg,targvals,ifppregtarg,pottarg,pretarg,gradtarg)
      call cpu_time(t2)
C$      t2 = omp_get_wtime()

      timeinfo(1) = t2-t1

c$OMP PARALLEL DO DEFAULT(SHARED)
c$OMP& PRIVATE(i,e11,e12,e13,e21,e22,e23,e31,e32,e33,dn1,dn2,dn3,p)
      do i = 1,ntarg
         e11 = gradtarg(1,1,i)*2
         e12 = gradtarg(1,2,i) + gradtarg(2,1,i)
         e13 = gradtarg(1,3,i) + gradtarg(3,1,i)
         e21 = e12
         e22 = gradtarg(2,2,i)*2
         e23 = gradtarg(2,3,i) + gradtarg(3,2,i)
         e31 = e13
         e32 = e23
         e33 = gradtarg(3,3,i)*2
         dn1 = targs(10,i)
         dn2 = targs(11,i)
         dn3 = targs(12,i)
         p = pretarg(i)
         
         pot(1,i) = e11*dn1 + e12*dn2 + e13*dn3 - p*dn1
         pot(2,i) = e21*dn1 + e22*dn2 + e23*dn3 - p*dn2
         pot(3,i) = e31*dn1 + e32*dn2 + e33*dn3 - p*dn3
      enddo
c$OMP END PARALLEL DO      

c
c        compute threshold for ignoring local computation
c
      
      xmin = sources(1,1)
      xmax = sources(1,1)
      ymin = sources(2,1)
      ymax = sources(2,1)
      zmin = sources(3,1)
      zmax = sources(3,1)

      do i=1,ns
        if(sources(1,i).lt.xmin) xmin = sources(1,i)
        if(sources(1,i).gt.xmax) xmax = sources(1,i)
        if(sources(2,i).lt.ymin) ymin = sources(2,i)
        if(sources(2,i).gt.ymax) ymax = sources(2,i)
        if(sources(3,i).lt.zmin) zmin = sources(3,i)
        if(sources(3,i).gt.zmax) zmax = sources(3,i)
      enddo

      do i=1,ntarg
        if(targvals(1,i).lt.xmin) xmin = targvals(1,i)
        if(targvals(1,i).gt.xmax) xmax = targvals(1,i)
        if(targvals(2,i).lt.ymin) ymin = targvals(2,i)
        if(targvals(2,i).gt.ymax) ymax = targvals(2,i)
        if(targvals(3,i).lt.zmin) zmin = targvals(3,i)
        if(targvals(3,i).gt.zmax) zmax = targvals(3,i)
      enddo
      
      boxsize = xmax-xmin
      sizey = ymax - ymin
      sizez = zmax - zmin

      if(sizey.gt.boxsize) boxsize = sizey
      if(sizez.gt.boxsize) boxsize = sizez

      thresh = 2.0d0**(-51)*boxsize

      
c
c
c       add in precomputed quadrature
c



      call cpu_time(t1)
C$      t1 = omp_get_wtime()

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart)
C$OMP$PRIVATE(jstart,pottmp,npols,w11,w12,w13,w21,w22,w23)
C$OMP$PRIVATE(w31,w32,w33,sig1,sig2,sig3)
      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
             sig1 = sigma(1,jstart+l-1)
             sig2 = sigma(2,jstart+l-1)
             sig3 = sigma(3,jstart+l-1)
             w11 = wnear(1,jquadstart+l-1)
             w12 = wnear(2,jquadstart+l-1)
             w13 = wnear(3,jquadstart+l-1)
             w21 = w12
             w22 = wnear(4,jquadstart+l-1)
             w23 = wnear(5,jquadstart+l-1)
             w31 = w13
             w32 = w23
             w33 = wnear(6,jquadstart+l-1)
             pot(1,i) = pot(1,i) + w11*sig1+w12*sig2+w13*sig3
             pot(2,i) = pot(2,i) + w21*sig1+w22*sig2+w23*sig3
             pot(3,i) = pot(3,i) + w31*sig1+w32*sig2+w33*sig3             
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO


c

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
C$OMP$PRIVATE(sttmp2,nss,l,jstart,ii,npover,gradv,tmp3,
C$OMP$   e11,e12,e13,e21,e22,e23,e31,e32,e33,dn1,dn2,dn3,p)
      do i=1,ntarg
         nss = 0
         do j=row_ptr(i),row_ptr(i+1)-1
            jpatch = col_ind(j)
            nss = nss + ixyzso(jpatch+1)-ixyzso(jpatch)
         enddo
         allocate(srctmp2(3,nss),sttmp2(3,nss))

         rmin = 1.0d6
         ii = 0
         do j=row_ptr(i),row_ptr(i+1)-1
            jpatch = col_ind(j)
            jstart = ixyzso(jpatch)-1
            npover = ixyzso(jpatch+1)-ixyzso(jpatch)
            do l=1,npover
               ii = ii+1
               srctmp2(1,ii) = srcover(1,jstart+l)
               srctmp2(2,ii) = srcover(2,jstart+l)
               srctmp2(3,ii) = srcover(3,jstart+l)

               sttmp2(1,ii) = stoklet(1,jstart+l)
               sttmp2(2,ii) = stoklet(2,jstart+l)
               sttmp2(3,ii) = stoklet(3,jstart+l)
            enddo
         enddo

         p=0
         gradv(1,1) = 0
         gradv(2,1) = 0
         gradv(3,1) = 0
         gradv(1,2) = 0
         gradv(2,2) = 0
         gradv(3,2) = 0
         gradv(1,3) = 0
         gradv(2,3) = 0
         gradv(3,3) = 0
         
         call st3ddirectstokg(nd,srctmp2,sttmp2,nss,
     1        targvals(1,i),ntarg0,tmp3,p,gradv,thresh)

         e11 = gradv(1,1)*2
         e12 = gradv(1,2) + gradv(2,1)
         e13 = gradv(1,3) + gradv(3,1)
         e21 = e12
         e22 = gradv(2,2)*2
         e23 = gradv(2,3) + gradv(3,2)
         e31 = e13
         e32 = e23
         e33 = gradv(3,3)*2
         dn1 = targs(10,i)
         dn2 = targs(11,i)
         dn3 = targs(12,i)

         pot(1,i) = pot(1,i) - (e11*dn1 + e12*dn2 + e13*dn3 - p*dn1)
         pot(2,i) = pot(2,i) - (e21*dn1 + e22*dn2 + e23*dn3 - p*dn2)
         pot(3,i) = pot(3,i) - (e31*dn1 + e32*dn2 + e33*dn3 - p*dn3)

         deallocate(srctmp2,sttmp2)
      enddo
      
      call cpu_time(t2)
C$      t2 = omp_get_wtime()     

      timeinfo(2) = t2-t1


cc      call prin2('quadrature time=*',timeinfo,2)
      
      ttot = timeinfo(1) + timeinfo(2)
cc      call prin2('time in lpcomp=*',ttot,1)

      
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
c
c
c
c        
      subroutine stok_s_trac_solver(npatches,norders,ixyzs,
     1    iptype,npts,srccoefs,srcvals,eps,numit,ifinout,
     2     rhs,eps_gmres,niter,errs,rres,soln)
c
c
c        this subroutine solves the Stokes traction problem
c     on the interior or exterior of an object where the potential
c     is represented as a single layer potential.
c
c
c     Representation:
c        u = S \sigma
c      
c     Equation:
c
c     1/2 \sigma(x) + t(S \sigma)[x]
c         + n(x) \int n(y) \cdot \sigma \, dS
c
c     the integral
c
c     \int n(y) \cdot \sigma \, dS
c
c     should be zero, so it does not affect the representation.
c     In practice, this removes a nullspace
c      
c     The linear system is solved iteratively using GMRES
c
c
c       input:
c         npatches - integer
c            number of patches
c
c         norders- integer(npatches)
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
c         srcvals - real *8 (12,npts)
c             xyz(u,v) and derivative info sampled at the 
c             discretization nodes on the surface
c             srcvals(1:3,i) - xyz info
c             srcvals(4:6,i) - dxyz/du info
c             srcvals(7:9,i) - dxyz/dv info
c 
c          eps - real *8
c             precision requested for computing quadrature and fmm
c             tolerance
c
c          ifinout - integer
c              flag for interior or exterior problems (normals assumed to 
c                be pointing in exterior of region)
c              ifinout = 0, interior problem
c              ifinout = 1, exterior problem
c
c           rhs - real *8(3,npts)
c              right hand side
c
c           eps_gmres - real *8
c                gmres tolerance requested
c
c           numit - integer
c              max number of gmres iterations
c
c         output
c           niter - integer
c              number of gmres iterations required for relative residual
c          
c           errs(1:iter) - relative residual as a function of iteration
c              number
c 
c           rres - real *8
c              relative residual for computed solution
c              
c           soln - real *8(3,npts)
c              density which solves the traction problem
c
c
      implicit none
      integer npatches,norder,npols,npts
      integer ifinout
      integer norders(npatches),ixyzs(npatches+1)
      integer iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,eps_gmres
      real *8 dpars(2)
      real *8 rhs(3*npts)
      real *8 soln(3*npts)

      real *8, allocatable :: targs(:,:)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:), wts(:)
      integer ndtarg,ntarg

      real *8 errs(numit+1)
      real *8 rres,eps2
      integer niter


      integer nover,npolso,nptso
      integer nnz,nquad
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      real *8, allocatable :: wnear(:,:)

      real *8, allocatable :: srcover(:,:),wover(:)
      integer, allocatable :: ixyzso(:),novers(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:) 

      integer i,j,jpatch,jquadstart,jstart

      integer ipars
      complex *16 zpars
      real *8 timeinfo(10),t1,t2,omp_get_wtime


      real *8 ttot,done,pi
      real *8 rfac,rfac0,sdotn
      integer iptype_avg,norder_avg
      integer ikerorder, iquadtype,npts_over

c
c
c       gmres variables
c
      integer nmat
      real *8 did,dtmp
      real *8 rb,wnrm2
      integer numit,it,iind,it1,k,l
      real *8 rmyerr
      real *8 temp
      real *8, allocatable :: vmat(:,:),hmat(:,:)
      real *8, allocatable :: cs(:),sn(:)
      real *8, allocatable :: svec(:),yvec(:),wtmp(:)

      complex *16 ztmp


      nmat = 3*npts
      allocate(vmat(nmat,numit+1),hmat(numit,numit))
      allocate(cs(numit),sn(numit))
      allocate(wtmp(nmat),svec(numit+1),yvec(numit+1))


      done = 1
      pi = atan(done)*4


c
c
c        setup targets as on surface discretization points
c 
      ndtarg = 12
      ntarg = npts
      allocate(targs(ndtarg,npts),uvs_targ(2,ntarg),ipatch_id(ntarg))

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        targs(1,i) = srcvals(1,i)
        targs(2,i) = srcvals(2,i)
        targs(3,i) = srcvals(3,i)
        targs(4,i) = srcvals(4,i)
        targs(5,i) = srcvals(5,i)
        targs(6,i) = srcvals(6,i)
        targs(7,i) = srcvals(7,i)
        targs(8,i) = srcvals(8,i)
        targs(9,i) = srcvals(9,i)
        targs(10,i) = srcvals(10,i)
        targs(11,i) = srcvals(11,i)
        targs(12,i) = srcvals(12,i)
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo
C$OMP END PARALLEL DO   


c
c    initialize patch_id and uv_targ for on surface targets
c
      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,
     1  ipatch_id,uvs_targ)

c
c
c        this might need fixing
c
      iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
      norder_avg = floor(sum(norders)/(npatches+0.0d0))

      call get_rfacs(norder_avg,iptype_avg,rfac,rfac0)


      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))

      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts, 
     1     srccoefs,cms,rads)

C$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,npatches
        rad_near(i) = rads(i)*rfac
      enddo
C$OMP END PARALLEL DO      

c
c    find near quadrature correction interactions
c
      print *, "entering find near mem"
      call findnearmem(cms,npatches,rad_near,ndtarg,targs,npts,nnz)
      print *, "nnz=",nnz

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,ndtarg,targs,npts,row_ptr, 
     1        col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind,
     1         iquad)

      ikerorder = -1
      if(abs(dpars(2)).gt.1.0d-16) ikerorder = 0


c
c    estimate oversampling for far-field, and oversample geometry
c

      allocate(novers(npatches),ixyzso(npatches+1))

      print *, "beginning far order estimation"

      ztmp = 0

      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,
     1    rads,npts,srccoefs,ndtarg,npts,targs,ikerorder,ztmp,
     2    nnz,row_ptr,col_ind,rfac,novers,ixyzso)

      npts_over = ixyzso(npatches+1)-1
      print *, "npts_over=",npts_over

      allocate(srcover(12,npts_over),wover(npts_over),wts(npts))

      call oversample_geom(npatches,norders,ixyzs,iptype,npts, 
     1   srccoefs,srcvals,novers,ixyzso,npts_over,srcover)

      call get_qwts(npatches,norders,ixyzs,iptype,npts,
     1        srcvals,wts)

      call get_qwts(npatches,novers,ixyzso,iptype,npts_over,
     1        srcover,wover)


c
c   compute near quadrature correction
c
      nquad = iquad(nnz+1)-1
      print *, "nquad=",nquad
      allocate(wnear(6,nquad))
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j)
      do i=1,nquad
         do j = 1,6
            wnear(j,i) = 0
         enddo
      enddo
C$OMP END PARALLEL DO    


      iquadtype = 1

      print *, "starting to generate near quadrature"
      call cpu_time(t1)
C$      t1 = omp_get_wtime()      

      call getnearquad_stok_s_trac(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     1      ipatch_id,uvs_targ,eps,iquadtype,nnz,row_ptr,col_ind,
     1      iquad,rfac0,nquad,wnear)
      call cpu_time(t2)
C$      t2 = omp_get_wtime()     

      call prin2('quadrature generation time=*',t2-t1,1)
      
      print *, "done generating near quadrature, now starting gmres"


c
c
c     start gmres code here
c
c     NOTE: matrix equation should be of the form (z*I + K)x = y
c       the identity scaling (z) is defined via did below,
c       and K represents the action of the principal value 
c       part of the matvec
c
      did = (-1)**(ifinout)*2*pi


      niter=0

c
c      compute norm of right hand side and initialize v
c 
      rb = 0

      do i=1,numit
        cs(i) = 0
        sn(i) = 0
      enddo


c
      do i=1,nmat
         rb = rb + abs(rhs(i))**2
      enddo
      rb = sqrt(rb)

      do i=1,nmat
        vmat(i,1) = rhs(i)/rb
      enddo

      svec(1) = rb

      do it=1,numit
         it1 = it + 1
         
c     
c     NOTE:
c     replace this routine by appropriate layer potential
c     evaluation routine  
c     


         call lpcomp_stok_s_trac_addsub(npatches,norders,ixyzs,
     1        iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     2        eps,nnz,row_ptr,col_ind,iquad,nquad,wnear,
     3        vmat(1,it),novers,npts_over,ixyzso,srcover,wover,
     4        wtmp)


         sdotn = 0
         do j = 1,npts
            sdotn = sdotn + vmat((j-1)*3+1,it)*srcvals(10,j)*wts(j)
            sdotn = sdotn + vmat((j-1)*3+2,it)*srcvals(11,j)*wts(j)
            sdotn = sdotn + vmat((j-1)*3+3,it)*srcvals(12,j)*wts(j)
         enddo

         do j = 1,npts
            wtmp((j-1)*3+1) = wtmp((j-1)*3+1) + srcvals(10,j)*sdotn
            wtmp((j-1)*3+2) = wtmp((j-1)*3+2) + srcvals(11,j)*sdotn
            wtmp((j-1)*3+3) = wtmp((j-1)*3+3) + srcvals(12,j)*sdotn
         enddo

         do k=1,it
            hmat(k,it) = 0
            do j=1,nmat
               hmat(k,it) = hmat(k,it) + wtmp(j)*vmat(j,k)
            enddo

            do j=1,nmat
               wtmp(j) = wtmp(j)-hmat(k,it)*vmat(j,k)
            enddo
         enddo
         
        hmat(it,it) = hmat(it,it)+did
        wnrm2 = 0
        do j=1,nmat
           wnrm2 = wnrm2 + abs(wtmp(j))**2
        enddo
        wnrm2 = sqrt(wnrm2)

        do j=1,nmat
           vmat(j,it1) = wtmp(j)/wnrm2
        enddo

        do k=1,it-1
           temp = cs(k)*hmat(k,it)+sn(k)*hmat(k+1,it)
           hmat(k+1,it) = -sn(k)*hmat(k,it)+cs(k)*hmat(k+1,it)
           hmat(k,it) = temp
        enddo

        dtmp = wnrm2

        call rotmat_gmres(hmat(it,it),dtmp,cs(it),sn(it))
          
        hmat(it,it) = cs(it)*hmat(it,it)+sn(it)*wnrm2
        svec(it1) = -sn(it)*svec(it)
        svec(it) = cs(it)*svec(it)
        rmyerr = abs(svec(it1))/rb
        errs(it) = rmyerr
        print *, "iter=",it,errs(it)

        if(rmyerr.le.eps_gmres.or.it.eq.numit) then

c     
c     solve the linear system corresponding to
c     upper triangular part of hmat to obtain yvec
c     
c     y = triu(H(1:it,1:it))\s(1:it);
c     
           do j=1,it
              iind = it-j+1
              yvec(iind) = svec(iind)
              do l=iind+1,it
                 yvec(iind) = yvec(iind) - hmat(iind,l)*yvec(l)
              enddo
              yvec(iind) = yvec(iind)/hmat(iind,iind)
           enddo

           
           
c     
c     estimate x
c
           do j=1,nmat
              soln(j) = 0
              do i=1,it
                 soln(j) = soln(j) + yvec(i)*vmat(j,i)
              enddo
           enddo
           

           rres = 0
           do i=1,nmat
              wtmp(i) = 0
           enddo
c     
c     NOTE:
c     replace this routine by appropriate layer potential
c     evaluation routine  
c     

           
           call lpcomp_stok_s_trac_addsub(npatches,norders,ixyzs,
     1          iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     2          eps,nnz,row_ptr,col_ind,iquad,nquad,wnear,
     3          soln,novers,npts_over,ixyzso,srcover,wover,
     4          wtmp)


           sdotn = 0
           do j = 1,npts
              sdotn = sdotn + soln((j-1)*3+1)*srcvals(10,j)*wts(j)
              sdotn = sdotn + soln((j-1)*3+2)*srcvals(11,j)*wts(j)
              sdotn = sdotn + soln((j-1)*3+3)*srcvals(12,j)*wts(j)
           enddo

           do j = 1,npts
              wtmp((j-1)*3+1) = wtmp((j-1)*3+1) + srcvals(10,j)*sdotn
              wtmp((j-1)*3+2) = wtmp((j-1)*3+2) + srcvals(11,j)*sdotn
              wtmp((j-1)*3+3) = wtmp((j-1)*3+3) + srcvals(12,j)*sdotn
           enddo
           
           do i=1,nmat
              rres = rres + abs(did*soln(i) + wtmp(i)-rhs(i))**2
           enddo
           rres = sqrt(rres)/rb
           niter = it
           return
        endif
      enddo
c     
      return
      end
c     
c
c
