c
c
c     This file contains the following user callable
c     routines: 
c 
c       getnearquad_helm_comb_dir - generates the near
c        field quadrature for the Dirichlet data
c        corresponding to the combined field
c        representation 
c
c       lpcomp_helm_comb_dir 
c          simpler version of helmholtz layer potential evaluator
c          only geometry, targets, representation parameters (alpha,beta,k)
c          and density sampled at discretization required on input,
c          output is the layer potential evaluated at the target points
c          (note that the identity term is not included for targets on
c           surface)
c
c       helm_comb_dir_solver - solves the interior/exterior Dirichlet
c         problem for Helmholtz equation using the combined field
c         representation
c
c       helm_comb_dir_fds_mem - get memory requirements for initialization
c          routine for subsequent calls to fast direct solver
c
c       helm_comb_dir_fds_init - initialize various arrays to be later 
c          used for fast direct solver
c
c       helm_comb_dir_fds_matgen - query entries of the combined field
c          representation matrix (input indices must be in column 
c          sparse compressed format, and must be preceeded by a call
c          to helm_comb_dir_fds_init)
c        
c
c    Advanced user interfaces: 
c*****************************
c       Note for developers: One of the two advanced user interfaces
c         is necessary for the easy user interface and
c         efficient iterative solvers. It seems that the add and subtract
c         version tends to have better CPU-time performance, but we expect
c         the setsub version to be more numerically stable
c**************************************
c       lpcomp_helm_comb_dir_addsub 
c         compute layer potential for the Dirichlet
c         data corresponding to the combined field representation
c         using add and subtract
c 
c
c       lpcomp_helm_comb_dir_setsub 
c          compute layer potential for Dirichlet data corresponding
c          to the combined field representation using 
c          set subtraction and turning off list 1
c
c
c
c



      subroutine getnearquad_helm_comb_dir(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2   ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,
     3   iquad,rfac0,nquad,wnear)
c
c       this subroutine generates the near field quadrature
c       for the representation u = 4\pi (\alpha S_{k}  + \beta D_{k}) ---(1)
c       where the near field is specified by the user 
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
c          zpars - complex *16 (3)
c              kernel parameters (Referring to formula (1))
c              zpars(1) = k 
c              zpars(2) = alpha
c              zpars(3) = beta
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
c               number of entries in wnear
c
c        output
c            wnear - complex *16(nquad)
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
      complex *16, intent(in) :: zpars(3)
      integer, intent(in) :: nnz
      integer, intent(in) :: row_ptr(ntarg+1),col_ind(nnz),iquad(nnz+1)
      complex *16, intent(out) :: wnear(nquad)


      integer ipars
      integer ndd,ndz,ndi
      real *8 dpars

      complex *16 alpha,beta
      integer i,j
      integer ipv

      procedure (), pointer :: fker
      external h3d_slp, h3d_dlp, h3d_comb

c
c
c        initialize the appropriate kernel function
c

      alpha = zpars(2)
      beta = zpars(3)

      ndz = 3
      ndi = 0
      ndd = 0
      if(iquadtype.eq.1) then
        fker => h3d_comb
        ipv = 1
        if(abs(alpha).ge.1.0d-16.and.abs(beta).lt.1.0d-16) then
          fker=>h3d_slp
          ipv = 0 
        else if(abs(alpha).lt.1.0d-16.and.abs(beta).ge.1.0d-16) then
          fker=>h3d_dlp
        endif


        call zgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     1     ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,
     1     ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear)
      endif


      if(abs(alpha).ge.1.0d-16.and.abs(beta).lt.1.0d-16) then
C$OMP PARALLEL DO DEFAULT(SHARED)        
        do i=1,nquad
          wnear(i) = wnear(i)*alpha
        enddo
C$OMP END PARALLEL DO        
      else if(abs(alpha).lt.1.0d-16.and.abs(beta).ge.1.0d-16) then
C$OMP PARALLEL DO DEFAULT(SHARED)        
        do i=1,nquad
          wnear(i) = wnear(i)*beta
        enddo
C$OMP END PARALLEL DO        
      endif

      return
      end
c
c
c
c
c
      subroutine lpcomp_helm_comb_dir(npatches,norders,ixyzs,
     1   iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2   ipatch_id,uvs_targ,eps,zpars,sigma,pot)
c
cf2py intent(in) npatches,norders,ixyzs,iptype,npts,srccoefs,srcvals
cf2py intent(in) ndtarg,ntarg,targs,ipatch_id,uvs_targ,eps,zpars
cf2py intent(in) sigma
cf2py intent(out) pot
c
c
c------------------------------
c  This subroutine evaluates the layer potential for the representation 
c
c
c  .. math ::
c  
c      u = 4 \pi(\alpha \mathcal{S}_{k} + \beta \mathcal{D}_{k})
c
c  Note: the $4\pi$ scaling is included to be consistent with the FMM3D
c  libraries. For targets on the boundary, this routine only computes
c  the principal value part, the identity term corresponding to the jump
c  in the layer potential is not included in the layer potential.
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
c    - zpars: double complex (3)
c        kernel parameters (Referring to formula above)
c        zpars(1) = k 
c        zpars(2) = $\alpha$
c        zpars(3) = $\beta$
c     - sigma: double complex(npts)
c         density for layer potential
c
c  Output arguments
c    - pot: double complex(ntarg)
c        layer potential evaluated at the target points
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
      complex *16, intent(in) :: zpars(3)
      complex *16, intent(in) :: sigma(npts)
      integer, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)

      complex *16, intent(out) :: pot(ntarg)


      integer nptso,nnz,nquad


      integer nover,npolso
      integer norder,npols
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      complex *16, allocatable :: wnear(:)

      real *8, allocatable :: srcover(:,:),wover(:)
      integer, allocatable :: ixyzso(:),novers(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:) 

      integer i,j,jpatch,jquadstart,jstart

      integer ipars
      real *8 dpars,timeinfo(10),t1,t2,omp_get_wtime


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
      call findnearmem(cms,npatches,rad_near,targs,ntarg,nnz)

      allocate(row_ptr(ntarg+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,targs,ntarg,row_ptr, 
     1        col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,ntarg,nnz,row_ptr,col_ind,
     1         iquad)

      ikerorder = -1
      if(abs(zpars(3)).gt.1.0d-16) ikerorder = 0

c
c    estimate oversampling for far-field, and oversample geometry
c

      allocate(novers(npatches),ixyzso(npatches+1))

      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,
     1    rads,npts,srccoefs,ndtarg,ntarg,targs,ikerorder,zpars(1),
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
      nquad = iquad(nnz+1)-1
      allocate(wnear(nquad))
      
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,nquad
        wnear(i) = 0
      enddo
C$OMP END PARALLEL DO    


      iquadtype = 1

      call getnearquad_helm_comb_dir(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     1      ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,
     1      iquad,rfac0,nquad,wnear)


c
c
c   compute layer potential
c
      call lpcomp_helm_comb_dir_addsub(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2  eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,
     3  sigma,novers,npts_over,ixyzso,srcover,wover,pot)



      return
      end
c
c
c
c
c
      subroutine lpcomp_helm_comb_dir_addsub(npatches,norders,ixyzs,
     1   iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2   eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,sigma,novers,
     3   nptso,ixyzso,srcover,whtsover,pot)
c
c
c      this subroutine evaluates the layer potential for
c      the representation u = 4\pi (\alpha S_{k} + \beta D_{k}) 
c      where the near field is precomputed and stored
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
c          zpars - complex *16 (3)
c              kernel parameters (Referring to formula (1))
c              zpars(1) = k 
c              zpars(2) = alpha
c              zpars(3) = beta
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
c           wnear - complex *16(nquad)
c               the near field quadrature correction
c
c           sigma - complex *16(npts)
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
c           pot - complex *16(npts)
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
      complex *16, intent(in) :: zpars(3)
      integer, intent(in) :: nnz,row_ptr(ntarg+1),col_ind(nnz),nquad
      integer, intent(in) :: iquad(nnz+1)
      complex *16, intent(in) :: wnear(nquad),sigma(npts)
      integer, intent(in) :: novers(npatches+1)
      integer, intent(in) :: nptso
      real *8, intent(in) :: srcover(12,nptso),whtsover(nptso)
      complex *16, intent(out) :: pot(ntarg)

      integer norder,npols,nover,npolso
      complex *16, allocatable :: potsort(:)

      real *8, allocatable :: sources(:,:),targvals(:,:)
      complex *16, allocatable :: charges(:),dipvec(:,:),sigmaover(:)
      integer ns,nt
      complex *16 alpha,beta
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      complex *16 tmp(10),val

      real *8 xmin,xmax,ymin,ymax,zmin,zmax,sizey,sizez,boxsize


      integer i,j,jpatch,jquadstart,jstart


      integer ifaddsub

      integer ntj
      
      complex *16 zdotu,pottmp
      complex *16, allocatable :: ctmp2(:),dtmp2(:,:)
      real *8 radexp,epsfmm

      integer ipars
      real *8 dpars,timeinfo(10),t1,t2,omp_get_wtime

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      real *8 thresh,ra
      real *8 rr,rmin
      integer nss,ii,l,npover

      integer nd,ntarg0

      real *8 ttot,done,pi

      parameter (nd=1,ntarg0=1)

      ns = nptso
      done = 1
      pi = atan(done)*4

           
      ifpgh = 0
      ifpghtarg = 1
      allocate(sources(3,ns),targvals(3,ntarg))
      allocate(charges(ns),dipvec(3,ns))
      allocate(sigmaover(ns))

c 
c       oversample density
c

      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype, 
     1    npts,sigma,novers,ixyzso,ns,sigmaover)


      ra = 0


c
c       set relevatn parameters for the fmm
c
      alpha = zpars(2)
      beta = zpars(3)
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)      
      do i=1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)

        charges(i) = sigmaover(i)*whtsover(i)*alpha
        dipvec(1,i) = sigmaover(i)*whtsover(i)*srcover(10,i)*beta
        dipvec(2,i) = sigmaover(i)*whtsover(i)*srcover(11,i)*beta
        dipvec(3,i) = sigmaover(i)*whtsover(i)*srcover(12,i)*beta
      enddo
C$OMP END PARALLEL DO      

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,ntarg
        targvals(1,i) = targs(1,i)
        targvals(2,i) = targs(2,i)
        targvals(3,i) = targs(3,i)
      enddo
C$OMP END PARALLEL DO      

      ifcharge = 1
      ifdipole = 1

      if(alpha.eq.0) ifcharge = 0
      if(beta.eq.0) ifdipole = 0

c
c
c       call the fmm
c

      call cpu_time(t1)
C$      t1 = omp_get_wtime()      
      call hfmm3d(nd,eps,zpars(1),ns,sources,ifcharge,charges,
     1  ifdipole,dipvec,ifpgh,tmp,tmp,tmp,ntarg,targvals,ifpghtarg,
     1  pot,tmp,tmp)
      call cpu_time(t2)
C$      t2 = omp_get_wtime()

      timeinfo(1) = t2-t1


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
C$OMP$PRIVATE(jstart,pottmp,npols)
      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            pot(i) = pot(i) + wnear(jquadstart+l-1)*sigma(jstart+l-1)
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO


c

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
C$OMP$PRIVATE(ctmp2,dtmp2,nss,l,jstart,ii,val,npover)
      do i=1,ntarg
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          nss = nss + ixyzso(jpatch+1)-ixyzso(jpatch)
        enddo
        allocate(srctmp2(3,nss),ctmp2(nss),dtmp2(3,nss))

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

            if(ifcharge.eq.1) ctmp2(ii) = charges(jstart+l)
            if(ifdipole.eq.1) then
              dtmp2(1,ii) = dipvec(1,jstart+l)
              dtmp2(2,ii) = dipvec(2,jstart+l)
              dtmp2(3,ii) = dipvec(3,jstart+l)
            endif
          enddo
        enddo

        val = 0
        if(ifcharge.eq.1.and.ifdipole.eq.0) then
          call h3ddirectcp(nd,zpars(1),srctmp2,ctmp2,
     1        nss,targvals(1,i),ntarg0,val,thresh)
        endif

        if(ifcharge.eq.0.and.ifdipole.eq.1) then
          call h3ddirectdp(nd,zpars(1),srctmp2,dtmp2,
     1          nss,targvals(1,i),ntarg0,val,thresh)
        endif

        if(ifcharge.eq.1.and.ifdipole.eq.1) then
          call h3ddirectcdp(nd,zpars(1),srctmp2,ctmp2,dtmp2,
     1          nss,targvals(1,i),ntarg0,val,thresh)
        endif
        pot(i) = pot(i) - val

        deallocate(srctmp2,ctmp2,dtmp2)
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
      subroutine lpcomp_helm_comb_dir_setsub(npatches,norders,ixyzs,
     1   iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2   eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,sigma,novers,
     3   nptso,ixyzso,srcover,whtsover,pot)
c
c
c      this subroutine evaluates the layer potential for
c      the representation u = (\alpha S_{k} + \beta D_{k}) 4 \pi 
c      where the near field is precomputed and stored
c      in the row sparse compressed format.
c
c
c          Note the 4\pi scaling is included to be consistent with fmm
c
c
c     The fmm is used to accelerate the far-field and 
c     near-field interactions are handled via precomputed quadrature
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
c          zpars - complex *16 (3)
c              kernel parameters (Referring to formula (1))
c              zpars(1) = k 
c              zpars(2) = alpha
c              zpars(3) = beta
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
c           wnear - complex *16(nquad)
c               the near field quadrature correction
c
c           sigma - complex *16(npts)
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
c               
c
      implicit none
      integer npatches,norder,npols,npts
      integer ndtarg,ntarg
      integer norders(npatches),ixyzs(npatches+1)
      integer ixyzso(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps
      real *8 targs(ndtarg,ntarg)
      complex *16 zpars(3)
      integer nnz,row_ptr(ntarg+1),col_ind(nnz),nquad
      integer iquad(nnz+1)
      complex *16 wnear(nquad),sigma(npts)
      integer novers(npatches+1)
      integer nover,npolso,nptso
      real *8 srcover(12,nptso),whtsover(nptso)
      complex *16 pot(ntarg)
      complex *16, allocatable :: potsort(:)

      real *8, allocatable :: sources(:,:),targvals(:,:)
      complex *16, allocatable :: charges(:),dipvec(:,:),sigmaover(:)
      integer, allocatable :: iboxtarg(:),iboxsrc(:)
      integer ns,nt
      complex *16 alpha,beta
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      complex *16 tmp(10),val


      integer i,j,jpatch,jquadstart,jstart


      integer ifaddsub

      integer *8 ltree,ipointer(32)
      integer, allocatable :: itree(:)
      integer, allocatable :: il1(:),il2(:),ilint(:),il1m2(:),il2m1(:)
      real *8, allocatable :: boxsize(:),centers(:,:)

      real *8 expc(3)
      integer ibox,nexpc,idivflag,iert,ifnear,ii,isend,isep,isource
      integer isstart,itarg,itend,itstart,itt,jbox,jpt,mhung,mnbors
      integer iss,l,lpt,mnlist1,mnlist2,mnlist3,mnlist4
      integer n1m2,n2m1,nadd,nbmax,nboxes,nchild,ndiv,nl2,nlevels
      integer nlist1,nlmax,npover,nl1,ntj
      
      complex *16 zdotu,pottmp
      complex *16, allocatable :: ctmp1(:),ctmp2(:),dtmp1(:,:),
     1   dtmp2(:,:)
      real *8 radexp,epsfmm

      integer ipars
      real *8 dpars,timeinfo(10),t1,t2,omp_get_wtime

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp1(:,:),srctmp2(:,:)
      real *8 thresh,ra
      integer nss

      integer nd,ntarg0

      real *8 ttot,done,pi

      parameter (nd=1,ntarg0=1)

      ns = nptso
      done = 1
      pi = atan(done)*4

           
      ifpgh = 0
      ifpghtarg = 1
      allocate(sources(3,ns),targvals(3,ntarg))
      allocate(charges(ns),dipvec(3,ns))
      allocate(sigmaover(ns))

c 
c       oversample density
c

      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype, 
     1    npts,sigma,novers,ixyzso,ns,sigmaover)
      call prinf('inside lpcomp, done oversampling density*',i,0)


      ra = 0


c
c       set relevatn parameters for the fmm
c
      alpha = zpars(2)
      beta = zpars(3)
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)      
      do i=1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)

        charges(i) = sigmaover(i)*whtsover(i)*alpha
        dipvec(1,i) = sigmaover(i)*whtsover(i)*srcover(10,i)*beta
        dipvec(2,i) = sigmaover(i)*whtsover(i)*srcover(11,i)*beta
        dipvec(3,i) = sigmaover(i)*whtsover(i)*srcover(12,i)*beta
      enddo
C$OMP END PARALLEL DO      

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,ntarg
        targvals(1,i) = targs(1,i)
        targvals(2,i) = targs(2,i)
        targvals(3,i) = targs(3,i)
      enddo
C$OMP END PARALLEL DO      

      ifcharge = 1
      ifdipole = 1

      if(alpha.eq.0) ifcharge = 0
      if(beta.eq.0) ifdipole = 0


c
c       setup tree
c
c

      isep = 1
      nlmax = 200
      nbmax = 0
      nlevels = 0
      nboxes = 0
      mhung = 0
      ltree = 0

      nexpc = 0
      nadd = 0
      ntj = 0
      
      idivflag = 0


      mnlist1 = 0
      mnlist2 = 0
      mnlist3 = 0
      mnlist4 = 0

      allocate(radsrc(ns))
      
      do i=1,ns
        radsrc(i) = 0
      enddo

      radexp = 0
      iert = 0

cc      ndiv = ns + ntarg

       if(eps.ge.0.5d-0) then
         ndiv = 300
       else if(eps.ge.0.5d-1) then
         ndiv = 300
       else if(eps.ge.0.5d-2) then
         ndiv = 300
       else if(eps.ge.0.5d-3) then
         ndiv = 300
       else if(eps.ge.0.5d-6) then
         ndiv = 1000
       else if(eps.ge.0.5d-9) then
         ndiv = 1000
       else if(eps.ge.0.5d-12) then
         ndiv = 1000
       else if(eps.ge.0.5d-15) then
         ndiv = 1000
       else
         ndiv = ns+ntarg
       endif

       ndiv = ndiv/4



      call mklraptreemem(iert,sources,ns,radsrc,targvals,ntarg,
     1   expc,nexpc,radexp,idivflag,ndiv,isep,nlmax,nbmax,
     2   nlevels,nboxes,mnbors,mnlist1,mnlist2,mnlist3,mnlist4,
     3   mhung,ltree)

      allocate(itree(ltree),boxsize(0:nlevels),centers(3,nboxes))

      call mklraptree(sources,ns,radsrc,targvals,ntarg,expc,nexpc,
     1  radexp,idivflag,ndiv,isep,mhung,mnbors,mnlist1,mnlist2,mnlist3,
     2  mnlist4,nlevels,nboxes,centers,boxsize,itree,ltree,ipointer)


      ifnear = 0


c
c
c       call the fmm
c

      call cpu_time(t1)
C$      t1 = omp_get_wtime()      
      call hfmm3d_ndiv(nd,eps,zpars(1),ns,sources,ifcharge,charges,
     1  ifdipole,dipvec,ifpgh,tmp,tmp,tmp,ntarg,targvals,ifpghtarg,
     1  pot,tmp,tmp,ndiv,idivflag,ifnear)
      call cpu_time(t2)
C$      t2 = omp_get_wtime()

      timeinfo(1) = t2-t1

      
c
c
c       add in precomputed quadrature
c

      thresh = 2.0d0**(-51)*boxsize(0)
      call cpu_time(t1)
C$      t1 = omp_get_wtime()

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart)
C$OMP$PRIVATE(jstart,pottmp,npols)
      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            pot(i) = pot(i) + wnear(jquadstart+l-1)*sigma(jstart+l-1)
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO



c
c
c    work with sorted potentials and unsort them again later
c
      allocate(potsort(ntarg))
      call dreorderf(2,ntarg,pot,potsort,itree(ipointer(6)))



c
c    subtract  precomputed near quadrature /setminus list1 
c       also needs to go from pts (targs) -> pts (sources)
c 
c
c    il1 - list of sources in the near field of a target (A)
c    il2 - list of sources in the list1 of the target from fmm
c        perspective (B)
c    il1m2 = A \cap (A \cap B)^{c}
c    il2m1 = B \cap (A \cap B)^{c}
c

     
      allocate(il2(ndiv*mnlist1),il2m1(ndiv*mnlist1))
      allocate(ctmp2(ndiv*mnlist1),dtmp2(3,ndiv*mnlist1))
      allocate(srctmp2(3,ndiv*mnlist1))

  

      call cpu_time(t1)
C$      t1 = omp_get_wtime()     

C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP$PRIVATE(ibox,nchild,nl2)
C$OMP$PRIVATE(nlist1,i,jbox,isstart,isend,j,isource,il2)
C$OMP$PRIVATE(itstart,itend,itt,itarg,nl1,il1,il1m2,il2m1)
C$OMP$PRIVATE(jpatch,l,jpt,lpt,n1m2,n2m1,ii,val,npover)
C$OMP$PRIVATE(ctmp1,ctmp2,dtmp1,dtmp2,srctmp1,srctmp2)
C$OMP$SCHEDULE(DYNAMIC)
      do ibox = 1,nboxes
        nchild = itree(ipointer(3)+ibox-1)
        if(nchild.eq.0) then

c
c     populate il2
c
          nl2 = 0
          nlist1 = itree(ipointer(20)+ibox-1)
          do i=1,nlist1
            jbox = itree(ipointer(21) + mnlist1*(ibox-1)+i-1)
            isstart = itree(ipointer(10)+jbox-1)
            isend = itree(ipointer(11)+jbox-1)
            do j=isstart,isend
              isource = itree(ipointer(5)+j-1)
              nl2 = nl2 + 1
              il2(nl2) = isource
            enddo
          enddo


c
c    end of populating il2.
c    
c    now loop over targets in this box
c
          itstart = itree(ipointer(12)+ibox-1)
          itend = itree(ipointer(13)+ibox-1)
          do itt = itstart,itend
            itarg = itree(ipointer(6)+itt-1)
            
            nl1 = 0
            do j=row_ptr(itarg),row_ptr(itarg+1)-1
              jpatch = col_ind(j)
              nl1 = nl1 + ixyzso(jpatch+1)-ixyzso(jpatch)
            enddo

            allocate(il1(nl1),il1m2(nl1),ctmp1(nl1),dtmp1(3,nl1))
            allocate(srctmp1(3,nl1))
c
c    populate il1 
c

            lpt = 0
            do j = row_ptr(itarg),row_ptr(itarg+1)-1
              jpatch = col_ind(j)
              npover = ixyzso(jpatch+1)-ixyzso(jpatch)
              do l=1,npover
                jpt = ixyzso(jpatch)+l-1
                lpt = lpt + 1
                il1(lpt) = jpt
              enddo
            enddo

cc            call prinf('il1=*',il1,nl1)
c
c   end of populating il1. now perform various set subtractions
c
            n1m2 = 0
            n2m1 = 0
            call setsub(il1,nl1,il2,nl2,il1m2,n1m2,il2m1,n2m1)


c
c   subtract off il1m2
c
c   gather step
c
            do i=1,n1m2
              ii = il1m2(i)
              srctmp1(1,i) = srcover(1,ii)
              srctmp1(2,i) = srcover(2,ii)
              srctmp1(3,i) = srcover(3,ii)
            enddo
            if(ifcharge.eq.1) then
              do i=1,n1m2
                ii = il1m2(i)
                ctmp1(i) = charges(ii)
              enddo
            endif

            if(ifdipole.eq.1) then
              do i=1,n1m2
                ii = il1m2(i)
                dtmp1(1,i) = dipvec(1,ii)
                dtmp1(2,i) = dipvec(2,ii)
                dtmp1(3,i) = dipvec(3,ii)
              enddo
            endif

            val = 0
            if(ifcharge.eq.1.and.ifdipole.eq.0) then
              call h3ddirectcp(nd,zpars(1),srctmp1,ctmp1,
     1          n1m2,targvals(1,itarg),ntarg0,val,thresh)
            endif

            if(ifcharge.eq.0.and.ifdipole.eq.1) then
              call h3ddirectdp(nd,zpars(1),srctmp1,dtmp1,
     1          n1m2,targvals(1,itarg),ntarg0,val,thresh)
            endif

            if(ifcharge.eq.1.and.ifdipole.eq.1) then
              call h3ddirectcdp(nd,zpars(1),srctmp1,ctmp1,dtmp1,
     1          n1m2,targvals(1,itarg),ntarg0,val,thresh)
            endif
c
c  scatter step
c
            potsort(itt) = potsort(itt) - val



c
c   add il2m1
c
c
c   gather step
c
            do i=1,n2m1
              ii = il2m1(i)
              srctmp2(1,i) = srcover(1,ii)
              srctmp2(2,i) = srcover(2,ii)
              srctmp2(3,i) = srcover(3,ii)
            enddo
            if(ifcharge.eq.1) then
              do i=1,n2m1
                ii = il2m1(i)
                ctmp2(i) = charges(ii)
              enddo
            endif

            if(ifdipole.eq.1) then
              do i=1,n2m1
                ii = il2m1(i)
                dtmp2(1,i) = dipvec(1,ii)
                dtmp2(2,i) = dipvec(2,ii)
                dtmp2(3,i) = dipvec(3,ii)
              enddo
            endif

            val = 0
            if(ifcharge.eq.1.and.ifdipole.eq.0) then
              call h3ddirectcp(nd,zpars(1),srctmp2,ctmp2,
     1          n2m1,targvals(1,itarg),ntarg0,val,thresh)
            endif

            if(ifcharge.eq.0.and.ifdipole.eq.1) then
              call h3ddirectdp(nd,zpars(1),srctmp2,dtmp2,
     1          n2m1,targvals(1,itarg),ntarg0,val,thresh)
            endif

            if(ifcharge.eq.1.and.ifdipole.eq.1) then
              call h3ddirectcdp(nd,zpars(1),srctmp2,ctmp2,dtmp2,
     1          n2m1,targvals(1,itarg),ntarg0,val,thresh)
            endif
c
c  scatter step
c
            potsort(itt) = potsort(itt) + val

            deallocate(il1,il1m2,ctmp1,dtmp1,srctmp1)
          enddo
        endif
      enddo
C$OMP END PARALLEL DO      
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      
      timeinfo(2) = t2-t1

      call dreorderi(2,ntarg,potsort,pot,itree(ipointer(6)))

cc      call prin2('quadrature time=*',timeinfo,2)
      
      ttot = timeinfo(1) + timeinfo(2)
cc      call prin2('time in lpcomp=*',ttot,1)

cc      call prin2('at end of lpcomp*',i,0)
cc      call prin2('pot=*',pot,24)
        
      
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
      subroutine helm_comb_dir_solver(npatches,norders,ixyzs,
     1    iptype,npts,srccoefs,srcvals,eps,zpars,numit,ifinout,
     2    rhs,eps_gmres,niter,errs,rres,soln)
c
c
c        this subroutine solves the helmholtz dirichlet problem
c     on the interior or exterior of an object where the potential
c     is represented as a combined field integral equation.
c
c
c     Representation:
c        u = \alpha S_{k} + \beta D_{k}
c     
c     The linear system is solved iteratively using GMRES
c     until a relative residual of 1e-15 is reached
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
c          zpars - complex *16 (3)
c              kernel parameters (Referring to formula (1))
c              zpars(1) = k 
c              zpars(2) = alpha
c              zpars(3) = beta
c
c          ifinout - integer
c              flag for interior or exterior problems (normals assumed to 
c                be pointing in exterior of region)
c              ifinout = 0, interior problem
c              ifinout = 1, exterior problem
c
c           rhs - complex *16(npts)
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
c           soln - complex *16(npts)
c              density which solves the dirichlet problem
c
c
      implicit none
      integer npatches,norder,npols,npts
      integer ifinout
      integer norders(npatches),ixyzs(npatches+1)
      integer iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,eps_gmres
      complex *16 zpars(3)
      complex *16 rhs(npts)
      complex *16 soln(npts)

      real *8, allocatable :: targs(:,:)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      integer ndtarg,ntarg

      real *8 errs(numit+1)
      real *8 rres,eps2
      integer niter


      integer nover,npolso,nptso
      integer nnz,nquad
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      complex *16, allocatable :: wnear(:)

      real *8, allocatable :: srcover(:,:),wover(:)
      integer, allocatable :: ixyzso(:),novers(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:) 

      integer i,j,jpatch,jquadstart,jstart

      integer ipars
      real *8 dpars,timeinfo(10),t1,t2,omp_get_wtime


      real *8 ttot,done,pi
      real *8 rfac,rfac0
      integer iptype_avg,norder_avg
      integer ikerorder, iquadtype,npts_over

c
c
c       gmres variables
c
      complex *16 zid,ztmp
      real *8 rb,wnrm2
      integer numit,it,iind,it1,k,l
      real *8 rmyerr
      complex *16 temp
      complex *16, allocatable :: vmat(:,:),hmat(:,:)
      complex *16, allocatable :: cs(:),sn(:)
      complex *16, allocatable :: svec(:),yvec(:),wtmp(:)


      allocate(vmat(npts,numit+1),hmat(numit,numit))
      allocate(cs(numit),sn(numit))
      allocate(wtmp(npts),svec(numit+1),yvec(numit+1))


      done = 1
      pi = atan(done)*4


c
c
c        setup targets as on surface discretization points
c 
      ndtarg = 3
      ntarg = npts
      allocate(targs(ndtarg,npts),uvs_targ(2,ntarg),ipatch_id(ntarg))

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        targs(1,i) = srcvals(1,i)
        targs(2,i) = srcvals(2,i)
        targs(3,i) = srcvals(3,i)
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
      call findnearmem(cms,npatches,rad_near,targs,npts,nnz)
      print *, "nnz=",nnz

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,targs,npts,row_ptr, 
     1        col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind,
     1         iquad)

      ikerorder = -1
      if(abs(zpars(3)).gt.1.0d-16) ikerorder = 0


c
c    estimate oversampling for far-field, and oversample geometry
c

      allocate(novers(npatches),ixyzso(npatches+1))

      print *, "beginning far order estimation"

      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,
     1    rads,npts,srccoefs,ndtarg,npts,targs,ikerorder,zpars(1),
     2    nnz,row_ptr,col_ind,rfac,novers,ixyzso)

      npts_over = ixyzso(npatches+1)-1
      print *, "npts_over=",npts_over
      call prinf('novers=*',novers,100)

      allocate(srcover(12,npts_over),wover(npts_over))

      call oversample_geom(npatches,norders,ixyzs,iptype,npts, 
     1   srccoefs,srcvals,novers,ixyzso,npts_over,srcover)

      call get_qwts(npatches,novers,ixyzso,iptype,npts_over,
     1        srcover,wover)


c
c   compute near quadrature correction
c
      nquad = iquad(nnz+1)-1
      print *, "nquad=",nquad
      allocate(wnear(nquad))
      
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,nquad
        wnear(i) = 0
      enddo
C$OMP END PARALLEL DO    


      iquadtype = 1

cc      eps2 = 1.0d-8

      print *, "starting to generate near quadrature"
      call cpu_time(t1)
C$      t1 = omp_get_wtime()      

      call getnearquad_helm_comb_dir(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     1      ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,
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
c       the identity scaling (z) is defined via zid below,
c       and K represents the action of the principal value 
c       part of the matvec
c
      zid = -(-1)**(ifinout)*2*pi*zpars(3)


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
      do i=1,npts
        rb = rb + abs(rhs(i))**2
      enddo
      rb = sqrt(rb)

      do i=1,npts
        vmat(i,1) = rhs(i)/rb
      enddo

      svec(1) = rb

      do it=1,numit
        it1 = it + 1

c
c        NOTE:
c        replace this routine by appropriate layer potential
c        evaluation routine  
c


        call lpcomp_helm_comb_dir_addsub(npatches,norders,ixyzs,
     1    iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     2    eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,
     3    vmat(1,it),novers,npts_over,ixyzso,srcover,wover,wtmp)

        do k=1,it
          hmat(k,it) = 0
          do j=1,npts
            hmat(k,it) = hmat(k,it) + wtmp(j)*conjg(vmat(j,k))
          enddo

          do j=1,npts
            wtmp(j) = wtmp(j)-hmat(k,it)*vmat(j,k)
          enddo
        enddo
          
        hmat(it,it) = hmat(it,it)+zid
        wnrm2 = 0
        do j=1,npts
          wnrm2 = wnrm2 + abs(wtmp(j))**2
        enddo
        wnrm2 = sqrt(wnrm2)

        do j=1,npts
          vmat(j,it1) = wtmp(j)/wnrm2
        enddo

        do k=1,it-1
          temp = cs(k)*hmat(k,it)+sn(k)*hmat(k+1,it)
          hmat(k+1,it) = -sn(k)*hmat(k,it)+cs(k)*hmat(k+1,it)
          hmat(k,it) = temp
        enddo

        ztmp = wnrm2

        call zrotmat_gmres(hmat(it,it),ztmp,cs(it),sn(it))
          
        hmat(it,it) = cs(it)*hmat(it,it)+sn(it)*wnrm2
        svec(it1) = -sn(it)*svec(it)
        svec(it) = cs(it)*svec(it)
        rmyerr = abs(svec(it1))/rb
        errs(it) = rmyerr
        print *, "iter=",it,errs(it)

        if(rmyerr.le.eps_gmres.or.it.eq.numit) then

c
c            solve the linear system corresponding to
c            upper triangular part of hmat to obtain yvec
c
c            y = triu(H(1:it,1:it))\s(1:it);
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
c          estimate x
c
          do j=1,npts
            soln(j) = 0
            do i=1,it
              soln(j) = soln(j) + yvec(i)*vmat(j,i)
            enddo
          enddo


          rres = 0
          do i=1,npts
            wtmp(i) = 0
          enddo
c
c        NOTE:
c        replace this routine by appropriate layer potential
c        evaluation routine  
c


          call lpcomp_helm_comb_dir_addsub(npatches,norders,ixyzs,
     1      iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     2      eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,
     3      soln,novers,npts_over,ixyzso,srcover,wover,wtmp)

            
          do i=1,npts
            rres = rres + abs(zid*soln(i) + wtmp(i)-rhs(i))**2
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
c
      subroutine helm_comb_dir_fds_mem(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,eps,zpars,nifds,nrfds,nzfds)
c
cf2py   intent(in) npatches,norders,ixyzs,iptype,npts
cf2py   intent(in) srccoefs,srcvals,eps,zpars
cf2py   intent(out) nifds,nrfds,nzfds

c
c       This subroutine estimates the memory requirements
c       for the precomputation routine of the fast direct solver
c 
c       The precomputation routine computes an integer array (ifds)
c       a real array (rfds), and complex array (zfds)
c 
c       The following quantities will be computed during the 
c       precomputation phase of the fast direct solver
c
c          ifds(1) - npts_over
c          ifds(2) - nnz
c          ifds(3) - nquad
c          ifds(4) - nximat
c          ifds(5:6+npts-1) - row_ptr (for near quadrature info)
c          ifds(6+npts:6+npts+nnz-1) - col_ind
c          ifds(6+npts+nnz:6+npts+2*nnz) - iquad
c          ifds(7+npts+2*nnz:7+npts+2*nnz+npatches-1) - novers
c          ifds(7+npts+2*nnz+npatches:7+npts+2*nnz+2*npatches) - ixyzso
c          ifds(8+npts+2*nnz+2*npatches:8+npts+2*nnz+3*npatches-1) -
c             iximat
c
c          rfds(1:12*npts_over) - srcover
c          rfds(12*npts_over+1:13*npts_over) - wover
c          rfds(13*npts_over+1:13*npts_over+nximat) - ximats)
c
c          zfds(1:3) - zpars(1:3)
c          zfds(4:4+nquad-1) - wnear
c
c        Thus this subroutine on output returns 
c          nifds = 6+npts+2*nnz+2*npatches
c          nrfds = 13*npts_over
c          nzfds = 3+nquad
c     
      implicit none
      integer npatches,norders(npatches),ixyzs(npatches+1)
      integer iptype(npatches),npts
      real *8 srccoefs(9,npts),srcvals(12,npts)
      real *8 eps
      complex *16 zpars(3)
      integer nifds,nrfds,nzfds
      integer nnz

      real *8, allocatable :: targs(:,:)
      integer iptype_avg,norder_avg
      integer ntarg,ndtarg,ikerorder
      real *8 rfac,rfac0
      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)
      integer, allocatable :: iquad(:),row_ptr(:),col_ind(:)
      integer, allocatable :: novers(:),ixyzso(:)
      integer npts_over,nquad,nximat

      integer i


c
c
c        setup targets as on surface discretization points
c 
      ndtarg = 3
      ntarg = npts
      allocate(targs(ndtarg,npts))

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        targs(1,i) = srcvals(1,i)
        targs(2,i) = srcvals(2,i)
        targs(3,i) = srcvals(3,i)
      enddo
C$OMP END PARALLEL DO   


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
      call findnearmem(cms,npatches,rad_near,targs,npts,nnz)

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,targs,npts,row_ptr, 
     1        col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind,
     1         iquad)

      ikerorder = -1
      if(abs(zpars(3)).gt.1.0d-16) ikerorder = 0


c
c    estimate oversampling for far-field, and oversample geometry
c

      allocate(novers(npatches),ixyzso(npatches+1))

      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,
     1    rads,npts,srccoefs,ndtarg,npts,targs,ikerorder,zpars(1),
     2    nnz,row_ptr,col_ind,rfac,novers,ixyzso)


      

      npts_over = ixyzso(npatches+1)-1

      nquad = iquad(nnz+1)-1

      call get_nximat(npatches,ixyzs,ixyzso,nximat)

      nifds = 7+npts+2*nnz+3*npatches
      nrfds = 13*npts_over+nximat
      nzfds = 3 + nquad
      

      return
      end
c
c
c
c
c

      subroutine helm_comb_dir_fds_init(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,eps,zpars,nifds,ifds,nrfds,
     2     rfds,nzfds,zfds)
cf2py   intent(in) npatches,norders,ixyzs,iptype,npts
cf2py   intent(in) srccoefs,srcvals,eps,zpars
cf2py   intent(in) nifds,nrfds,nzfds
cf2py   intent(out) ifds,rfds,zfds

c
c       This subroutine is the precomputation routine of the fast direct solver
c 
c       The following quantities will be computed during the 
c       precomputation phase of the fast direct solver
c
c
c          ifds(1) - npts_over
c          ifds(2) - nnz
c          ifds(3) - nquad
c          ifds(4) - nximat
c          ifds(5:6+npts-1) - row_ptr (for near quadrature info)
c          ifds(6+npts:6+npts+nnz-1) - col_ind
c          ifds(6+npts+nnz:6+npts+2*nnz) - iquad
c          ifds(7+npts+2*nnz:7+npts+2*nnz+npatches-1) - novers
c          ifds(7+npts+2*nnz+npatches:7+npts+2*nnz+2*npatches) - ixyzso
c          ifds(8+npts+2*nnz+2*npatches:8+npts+2*nnz+3*npatches-1) -
c             iximat
c
c          rfds(1:12*npts_over) - srcover
c          rfds(12*npts_over+1:13*npts_over) - wover
c          rfds(13*npts_over+1:13*npts_over+nximat) - ximats
c
c          zfds(1:3) - zpars(1:3)
c          zfds(4:4+nquad-1) - wnear
c
c        Thus this subroutine on output returns 
c          nifds = 6+npts+2*nnz+2*npatches
c          nrfds = 13*npts_over
c          nzfds = 3+nquad
c     
      implicit none
      integer npatches,norders(npatches),ixyzs(npatches+1)
      integer iptype(npatches),npts
      real *8 srccoefs(9,npts),srcvals(12,npts)
      real *8 eps
      complex *16 zpars(3)
      integer nnz
      integer nifds,nrfds,nzfds
      integer ifds(nifds)
      real *8 rfds(nrfds)
      complex *16 zfds(nzfds)

      real *8, allocatable :: targs(:,:)
      real *8, allocatable :: uvs_targ(:,:)
      integer, allocatable :: ipatch_id(:)
      integer iptype_avg,norder_avg
      integer ntarg,ndtarg,ikerorder
      real *8 rfac,rfac0
      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)
      integer npts_over,nquad

      integer iquadtype,istart,iend,i

      integer irow_ptr,icol_ind,iiquad,inovers,iixyzso,iximat
      integer nximat




c
c
c        setup targets as on surface discretization points
c 
      ndtarg = 3
      ntarg = npts
      allocate(targs(ndtarg,npts))

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        targs(1,i) = srcvals(1,i)
        targs(2,i) = srcvals(2,i)
        targs(3,i) = srcvals(3,i)
      enddo
C$OMP END PARALLEL DO   


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
      call findnearmem(cms,npatches,rad_near,targs,npts,nnz)

      irow_ptr = 5
      icol_ind = irow_ptr+npts+1
      iiquad = icol_ind+nnz
      inovers = iiquad + nnz+1
      iixyzso = inovers+npatches
      iximat = iixyzso+npatches+1

      
      call findnear(cms,npatches,rad_near,targs,npts,ifds(irow_ptr), 
     1        ifds(icol_ind))

      call get_iquad_rsc(npatches,ixyzs,npts,nnz,ifds(irow_ptr),
     1   ifds(icol_ind),ifds(iiquad))


      ikerorder = -1
      if(abs(zpars(3)).gt.1.0d-16) ikerorder = 0


c
c    estimate oversampling for far-field, and oversample geometry
c


      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,
     1    rads,npts,srccoefs,ndtarg,npts,targs,ikerorder,zpars(1),
     2    nnz,ifds(irow_ptr),ifds(icol_ind),rfac,ifds(inovers),
     3    ifds(iixyzso))

      npts_over = ifds(iixyzso+npatches)-1
      nquad = ifds(iiquad+nnz)-1

      call oversample_geom(npatches,norders,ixyzs,iptype,npts, 
     1   srccoefs,srcvals,ifds(inovers),ifds(iixyzso),npts_over,rfds)

      call get_qwts(npatches,ifds(inovers),ifds(iixyzso),iptype,
     1      npts_over,rfds,rfds(12*npts_over+1))

      

      ifds(1) = npts_over
      ifds(2) = nnz
      ifds(3) = nquad

      zfds(1:3) = zpars

c
c    initialize patch_id and uv_targ for on surface targets
c
      allocate(ipatch_id(ntarg),uvs_targ(2,ntarg))
      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,
     1  ipatch_id,uvs_targ)

      iquadtype = 1

      call getnearquad_helm_comb_dir(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     1      ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,ifds(irow_ptr),
     1      ifds(icol_ind),ifds(iiquad),rfac0,nquad,zfds(4))
      

      call get_nximat(npatches,ixyzs,ifds(iixyzso),nximat)

      ifds(4) = nximat


      call get_ximats(npatches,iptype,norders,ixyzs,ifds(inovers),
     1  ifds(iixyzso),nximat,rfds(13*npts_over+1),ifds(iximat))

      return
      end
c
c
c
c
c
c

      subroutine helm_comb_dir_fds_matgen(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,eps,zpars,nifds,ifds,nrfds,
     2     rfds,nzfds,zfds,nent_csc,col_ptr,row_ind,zmatent)
cf2py   intent(in) npatches,norders,ixyzs,iptype,npts
cf2py   intent(in) srccoefs,srcvals,eps,zpars
cf2py   intent(in) nifds,nrfds,nzfds
cf2py   intent(in) ifds,rfds,zfds
cf2py   intent(in) nent_csc,col_ptr,row_ind
cf2py   intent(out) zmatent
      
      implicit none
      integer npatches,norders(npatches),ixyzs(npatches+1)
      integer iptype(npatches),npts
      real *8 srccoefs(9,npts),srcvals(12,npts)
      real *8 eps
      complex *16 zpars(3)
      integer nifds,nrfds,nzfds
      integer ifds(nifds)
      real *8 rfds(nrfds)
      complex *16 zfds(nzfds)
      integer nent_csc,col_ptr(npts+1),row_ind(nent_csc)
      complex *16 zmatent(nent_csc)

c
c        temporary variables
c
      complex *16 zid
      integer i,j,k,l,ipatch,npols
      integer ndd,ndz,ndi
      integer ilstart,ilend,istart,iend
      integer irow_ptr,icol_ind,iiquad,inovers,iixyzso
      integer npts_over,nnz,nquad

      integer, allocatable :: iuni(:),iuniind(:)
      integer, allocatable :: col_ptr_src(:),row_ind_src(:),iper(:)
      integer, allocatable :: aintb(:),iaintba(:),aintbc(:),iaintbc(:)
      complex *16, allocatable :: wquadn(:,:),wquad(:,:)
      complex *16, allocatable :: wquadf(:,:),wquadf2(:,:)


      complex *16 alpha, beta,zk
      real *8, allocatable :: srcover(:,:),wtsover(:)
      real *8 dpars
      integer ipars,ipt,iquad,itind,iximat,ixist,j2,jind
      integer jind0,juniind,n2,naintb,naintbc,nmax,nn,npolso
      integer nuni,nximat
      integer, allocatable :: iaintbb(:)
      integer ifcharge,ifdipole

      procedure (), pointer :: fker
      external h3d_slp, h3d_dlp, h3d_comb

c
c
c        initialize the appropriate kernel function
c
 
      
      alpha = zfds(2)
      beta = zfds(3)
      fker => h3d_comb
      ndz = 3
      ifcharge = 1
      ifdipole = 1
      if(abs(alpha).ge.1.0d-16.and.abs(beta).lt.1.0d-16) then
        fker=>h3d_slp
        ndz = 1
        ifdipole = 0
      else if(abs(alpha).lt.1.0d-16.and.abs(beta).ge.1.0d-16) then
        fker=>h3d_dlp
        ndz = 1
        ifcharge = 0
      endif

      ndd = 0
      ndi = 0

      npts_over = ifds(1)
      nnz = ifds(2)
      nquad = ifds(3)
      nximat = ifds(4)


      irow_ptr = 5
      icol_ind = irow_ptr+npts+1
      iiquad = icol_ind+nnz
      inovers = iiquad + nnz+1
      iixyzso = inovers+npatches
      iximat = iixyzso+npatches+1

      allocate(col_ptr_src(npatches+1),row_ind_src(nnz),iper(nnz))

      call rsc_to_csc(npatches,npts,nnz,ifds(irow_ptr),ifds(icol_ind),
     1  col_ptr_src,row_ind_src,iper)

c
c    estimate max oversampling
c
      nmax = 0
      do ipatch=1,npatches
        npolso = ifds(iixyzso+ipatch)-ifds(iixyzso+ipatch-1)
        if(npolso.gt.nmax) nmax = npolso
      enddo

      allocate(srcover(12,nmax),wtsover(nmax))
     

      do ipatch=1,npatches

c
c
c   combine all list of targets requested for current
c   patch and find unique list of targets
c
      
        istart = ixyzs(ipatch)
        iend = ixyzs(ipatch+1)-1
        npols= iend-istart+1
         
        ilstart = col_ptr(istart)
        ilend = col_ptr(iend+1)-1

        nn = ilend-ilstart+1

        allocate(iuni(nn),iuniind(nn))

        call get_iuni1(nn,row_ind(ilstart),nuni,iuni,iuniind)
        
        allocate(aintb(nuni),iaintba(nuni),aintbc(nuni),iaintbc(nuni))
        allocate(iaintbb(nuni))

        n2 = col_ptr_src(ipatch+1)-col_ptr_src(ipatch)

c
c
c    separate list of targets into near field targets for which 
c    quadrature is already computed and far-field targets for which
c    quadrature is required to be computed
c
        naintb = 0
        naintbc = 0

        call setdecomp(nuni,iuni,n2,
     1     row_ind_src(col_ptr_src(ipatch)),naintb,aintb,iaintba,
     2     iaintbb,naintbc,aintbc,iaintbc)

       
c
c    for the entries in aintb, the quadrature has already been computed
c    so that needs to be extracted and sent to appropriate entries
c  
        allocate(wquad(nuni,npols))
        do i=1,naintb
           jind0 = iaintbb(i)+col_ptr_src(ipatch)-1
           jind = iper(jind0)
           iquad = ifds(iiquad + jind-1)
           wquad(iaintba(i),:) = zfds(3+iquad:3+iquad+npols-1)
        enddo

c
c        compute the oversampled quadrature
c
        istart = ifds(iixyzso+ipatch-1)
        iend = ifds(iixyzso+ipatch)-1
        npolso = iend-istart+1

        allocate(wquadf(npolso,naintbc))
        allocate(wquadf2(npols,naintbc))

c
c      extract srcover, wtsover, and targvals 
c

        do i=1,npolso
          do j=1,12
            srcover(j,i) = rfds(12*(istart+i-2)+j)
          enddo
          wtsover(i) = rfds(12*npts_over+istart+i-1)
        enddo

cc        call prin2('srcover=*',srcover,12*npolso)
cc        call prin2('wtsover=*',wtsover,npolso)


cc        call prin2('zfds=*',zfds,6)

        do i=1,naintbc
          itind = aintbc(i)
          do j=1,npolso
            call fker(srcover(1,j),3,srcvals(1,itind),ndd,dpars,ndz,
     1        zfds,ndi,ipars,wquadf(j,i))
            wquadf(j,i) = wquadf(j,i)*wtsover(j)
          enddo
        enddo

        if(ifdipole.eq.0) then
          do i=1,naintbc
            do j=1,npolso
              wquadf(j,i) = wquadf(j,i)*alpha
            enddo
          enddo
        endif

        if(ifcharge.eq.0) then
          do i=1,naintbc
            do j=1,npolso
              wquadf(j,i) = wquadf(j,i)*beta
            enddo
          enddo
        endif

c
c      now multiply wquadf by ximat
c
        ixist = ifds(iximat+ipatch-1) + 13*npts_over
        call zrmatmatt_slow(naintbc,npolso,npols,wquadf,rfds(ixist),
     1        wquadf2)
        
        do i=1,naintbc
          wquad(iaintbc(i),:) = wquadf2(:,i)
        enddo
       
        do i = 1,npols
          ipt = ixyzs(ipatch) + i-1
          do j=col_ptr(ipt),col_ptr(ipt+1)-1
             jind = row_ind(j)

             j2 = j-col_ptr(ixyzs(ipatch))+1
             juniind = iuniind(j2)
             zmatent(j) = wquad(juniind,i)
          enddo
        enddo

        deallocate(iuni,iuniind,aintb,iaintba,iaintbb,aintbc,iaintbc)
        deallocate(wquad,wquadf,wquadf2)
      enddo
      
      
      


      return
      end
