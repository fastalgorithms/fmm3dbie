c
c
c     This file contains the following user callable
c     routines: 
c 
c       getnearquad_lap_comb_dir - generates the near
c        field quadrature for the Dirichlet data
c        corresponding to the combined field
c        representation 
c
c       lpcomp_lap_comb_dir 
c          simpler version of Laplace layer potential evaluator
c          only geometry, targets, representation parameters (alpha,beta)
c          and density sampled at discretization required on input,
c          output is the layer potential evaluated at the target points
c          (note that the identity term is not included for targets on
c           surface)
c
c       lap_comb_dir_solver - solves the interior/exterior Dirichlet
c         problem for Laplace's equation using the combined field
c         representation
c
c       lap_comb_cap_solver - solver for the capacitance problem
c          for Laplace equation using modified Dirichlet approach 
c
c    Advanced user interfaces: 
c*****************************
c       Note for developers: One of the two advanced user interfaces
c         is necessary for the easy user interface and
c         efficient iterative solvers. It seems that the add and subtract
c         version tends to have better CPU-time performance, but we expect
c         the setsub version to be more numerically stable
c**************************************
c       lpcomp_lap_comb_dir_addsub 
c         compute layer potential for the Dirichlet
c         data corresponding to the combined field representation
c         using add and subtract
c 
c
c
c



      subroutine getnearquad_lap_comb_dir(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2   ipatch_id,uvs_targ,eps,dpars,iquadtype,nnz,row_ptr,col_ind,
     3   iquad,rfac0,nquad,wnear)
c
c       this subroutine generates the near field quadrature
c       for the representation u = (\alpha S_{0}  + \beta D_{0}) ---(1)
c       where the near field is specified by the user 
c       in row sparse compressed format.
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
c       input:
c         npatches - integer *8
c            number of patches
c
c         norders - integer *8(npatches)
c            order of discretization on each patch 
c
c         ixyzs - integer *8(npatches+1)
c            starting location of data on patch i
c  
c         iptype - integer *8(npatches)
c           type of patch
c           iptype = 1 -> triangular patch discretized with RV nodes
c
c         npts - integer *8
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
c         ndtarg - integer *8
c            leading dimension of target array
c        
c         ntarg - integer *8
c            number of targets
c
c         targs - real *8 (ndtarg,ntarg)
c            target information
c
c         ipatch_id - integer *8(ntarg)
c            id of patch of target i, id = -1, if target is off-surface
c
c         uvs_targ - real *8 (2,ntarg)
c            local uv coordinates on patch if on surface, otherwise
c            set to 0 by default
c            
c          eps - real *8
c             precision requested
c
c          dpars - double precision (2)
c              kernel parameters (Referring to formula (1))
c              dpars(1) = alpha
c              dpars(2) = beta
c
c           iquadtype - integer *8
c              quadrature type
c              iquadtype = 1, use ggq for self + adaptive integration
c                 for rest
c 
c
c           nnz - integer *8
c             number of source patch-> target interactions in the near field
c 
c           row_ptr - integer *8(ntarg+1)
c              row_ptr(i) is the pointer
c              to col_ind array where list of relevant source patches
c              for target i start
c
c           col_ind - integer *8 (nnz)
c               list of source patches relevant for all targets, sorted
c               by the target number
c
c           iquad - integer *8(nnz+1)
c               location in wnear array where quadrature for col_ind(i)
c               starts
c
c           rfac0 - integer *8
c               radius parameter for near field
c
c           nquad - integer *8
c               number of entries in wnear
c
c        output
c            wnear - real *16(nquad)
c               the desired near field quadrature
c               
c

      implicit none 
      integer *8, intent(in) :: npatches,norders(npatches),npts,nquad
      integer *8, intent(in) :: ixyzs(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: ndtarg,ntarg
      integer *8, intent(in) :: iquadtype
      real *8, intent(in) :: targs(ndtarg,ntarg)
      integer *8, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      real *8, intent(in) :: dpars(2)
      integer *8, intent(in) :: nnz
      integer *8, intent(in) :: row_ptr(ntarg+1),col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(out) :: wnear(nquad)


      integer *8 ipars
      integer *8 ndd,ndz,ndi
      complex *16 zpars

      real *8 alpha,beta
      integer *8 i,j
      integer *8 ipv

      procedure (), pointer :: fker
      external l3d_slp, l3d_dlp, l3d_comb

c
c
c        initialize the appropriate kernel function
c

      alpha = dpars(1)
      beta = dpars(2)

      ndd = 2
      ndi = 0
      ndz = 0
      if(iquadtype.eq.1) then
        fker => l3d_comb
        ipv = 1
        if(abs(alpha).ge.1.0d-16.and.abs(beta).lt.1.0d-16) then
          fker=>l3d_slp
          ipv = 0 
        else if(abs(alpha).lt.1.0d-16.and.abs(beta).ge.1.0d-16) then
          fker=>l3d_dlp
        endif

        call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     1     ipatch_id,uvs_targ,
     1     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     1     col_ind,iquad,
     1     rfac0,nquad,wnear)
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
      subroutine lpcomp_lap_comb_dir(npatches,norders,ixyzs,
     1   iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2   ipatch_id,uvs_targ,eps,dpars,sigma,pot)
c
cf2py intent(in) npatches,norders,ixyzs,iptype,npts,srccoefs,srcvals
cf2py intent(in) ndtarg,ntarg,targs,ipatch_id,uvs_targ,eps,dpars
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
c      u = (\alpha \mathcal{S}_{k} + \beta \mathcal{D}_{k})
c
c  Note: For targets on the boundary, this routine only computes
c  the principal value part, the identity term corresponding to the jump
c  in the layer potential is not included in the layer potential.
c
c
c  Input arguments:
c
c    - npatches: integer *8
c        number of patches
c    - norders: integer *8(npatches)
c        order of discretization on each patch 
c    - ixyzs: integer *8(npatches+1)
c        ixyzs(i) denotes the starting location in srccoefs,
c        and srcvals array where information for patch i begins
c    - iptype: integer *8(npatches)
c        type of patch
c    - npts: integer *8
c        total number of discretization points on the boundary
c    - srccoefs: double precision (9,npts)
c        koornwinder expansion coefficients of x, $\partial_{u} x$,
c        and $\partial_{v} x$. 
c    - srcvals: double precision (12,npts)
c        x, $\partial_{u} x$, $\partial_{v} x$, and $n$ sampled at
c        discretization nodes
c    - ndtarg: integer *8
c        leading dimension of target array
c    - ntarg: integer *8
c        number of targets
c    - targs: double precision (ndtarg,ntarg)
c        target information
c    - ipatch_id: integer *8(ntarg)
c        id of patch of target i, id = -1, if target is off-surface
c    - uvs_targ: double precision (2,ntarg)
c        local uv coordinates on patch if on surface, otherwise
c        set to 0 by default
c    - eps: double precision
c        precision requested
c    - dpars: double complex (2)
c        kernel parameters (Referring to formula above)
c        dpars(1) = $\alpha$
c        dpars(2) = $\beta$
c     - sigma: double precision(npts)
c         density for layer potential
c
c  Output arguments
c    - pot: double precision(ntarg)
c        layer potential evaluated at the target points
c
c-----------------------------------
c
      implicit none
      integer *8, intent(in) :: npatches,npts
      integer *8, intent(in) :: ndtarg,ntarg
      integer *8, intent(in) :: norders(npatches),ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      real *8, intent(in) :: targs(ndtarg,ntarg)
      real *8, intent(in) :: dpars(2)
      real *8, intent(in) :: sigma(npts)
      integer *8, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)

      real *8, intent(out) :: pot(ntarg)


      integer *8 nptso,nnz,nquad


      integer *8 nover,npolso
      integer *8 norder,npols
      integer *8, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      real *8, allocatable :: wnear(:)

      real *8, allocatable :: srcover(:,:),wover(:)
      integer *8, allocatable :: ixyzso(:),novers(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:) 

      integer *8 i,j,jpatch,jquadstart,jstart

      complex *16 zpars
      integer *8 ipars
      real *8 timeinfo(10),t1,t2,omp_get_wtime


      real *8 ttot,done,pi
      real *8 rfac,rfac0
      integer *8 iptype_avg,norder_avg
      integer *8 ikerorder, iquadtype,npts_over

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
      call findnearmem(cms,npatches,rad_near,ndtarg,targs,
     1   ntarg,nnz)

      allocate(row_ptr(ntarg+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,ndtarg,targs,ntarg,
     1   row_ptr,col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,ntarg,nnz,row_ptr,col_ind,
     1         iquad)

      ikerorder = -1
      if(abs(dpars(2)).gt.1.0d-16) ikerorder = 0


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
      nquad = iquad(nnz+1)-1
      allocate(wnear(nquad))
      
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,nquad
        wnear(i) = 0
      enddo
C$OMP END PARALLEL DO    


      iquadtype = 1

      call getnearquad_lap_comb_dir(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     1      ipatch_id,uvs_targ,eps,dpars,iquadtype,nnz,row_ptr,col_ind,
     1      iquad,rfac0,nquad,wnear)


c
c
c   compute layer potential
c
      call lpcomp_lap_comb_dir_addsub(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2  eps,dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,
     3  sigma,novers,npts_over,ixyzso,srcover,wover,pot)



      return
      end
c
c
c
c
c
      subroutine lpcomp_lap_comb_dir_addsub(npatches,norders,ixyzs,
     1   iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2   eps,dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,sigma,novers,
     3   nptso,ixyzso,srcover,whtsover,pot)
c
c
c      this subroutine evaluates the layer potential for
c      the representation u = (\alpha S_{0} + \beta D_{0}) 
c      where the near field is precomputed and stored
c      in the row sparse compressed format.
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
c         npatches - integer *8
c            number of patches
c
c         norders- integer *8(npatches)
c            order of discretization on each patch 
c
c         ixyzs - integer *8(npatches+1)
c            ixyzs(i) denotes the starting location in srccoefs,
c               and srcvals array corresponding to patch i
c   
c         iptype - integer *8(npatches)
c            type of patch
c             iptype = 1, triangular patch discretized using RV nodes
c
c         npts - integer *8
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
c         ndtarg - integer *8
c            leading dimension of target array
c        
c         ntarg - integer *8
c            number of targets
c
c         targs - real *8 (ndtarg,ntarg)
c            target information
c
c          eps - real *8
c             precision requested
c
c          dpars - real *8 (2)
c              kernel parameters (Referring to formula (1))
c              dpars(1) = alpha
c              dpars(2) = beta
c
c           nnz - integer *8 8
c             number of source patch-> target interactions in the near field
c 
c           row_ptr - integer *8(ntarg+1)
c              row_ptr(i) is the pointer
c              to col_ind array where list of relevant source patches
c              for target i start
c
c           col_ind - integer *8 (nnz)
c               list of source patches relevant for all targets, sorted
c               by the target number
c
c           iquad - integer *8(nnz+1)
c               location in wnear array where quadrature for col_ind(i)
c               starts
c
c           nquad - integer *8
c               number of entries in wnear
c
c           wnear - real *8(nquad)
c               the near field quadrature correction
c
c           sigma - real *8(npts)
c               density for layer potential
c
c           novers - integer *8(npatches)
c              order of discretization for oversampled sources and
c               density
c
c         ixyzso - integer *8(npatches+1)
c            ixyzso(i) denotes the starting location in srcover,
c               corresponding to patch i
c   
c           nptso - integer *8
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
c           pot - real *8(npts)
c              layer potential evaluated at the target points
c
c           
c               
c
      implicit none
      integer *8, intent(in) :: npatches,npts
      integer *8, intent(in) :: ndtarg,ntarg
      integer *8, intent(in) :: norders(npatches),ixyzs(npatches+1)
      integer *8, intent(in) :: ixyzso(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      real *8, intent(in) :: targs(ndtarg,ntarg)
      real *8, intent(in) :: dpars(2)
      integer *8, intent(in) :: nnz,row_ptr(ntarg+1),col_ind(nnz),nquad
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: wnear(nquad),sigma(npts)
      integer *8, intent(in) :: novers(npatches+1)
      integer *8, intent(in) :: nptso
      real *8, intent(in) :: srcover(12,nptso),whtsover(nptso)
      real *8, intent(out) :: pot(ntarg)

      integer *8 norder,npols,nover,npolso
      real *8, allocatable :: potsort(:)

      real *8, allocatable :: sources(:,:),targvals(:,:)
      real *8, allocatable :: charges(:),dipvec(:,:),sigmaover(:)
      integer *8 ns,nt
      real *8 alpha,beta
      integer *8 ifcharge,ifdipole
      integer *8 ifpgh,ifpghtarg
      real *8 tmp(10),val
      real *8 over4pi

      real *8 xmin,xmax,ymin,ymax,zmin,zmax,sizey,sizez,boxsize


      integer *8 i,j,jpatch,jquadstart,jstart


      integer *8 ifaddsub

      integer *8 ntj
      
      real *8 ddot,pottmp
      real *8, allocatable :: ctmp2(:),dtmp2(:,:)
      real *8 radexp,epsfmm

      integer *8 ipars,nmax
      complex *16 zpars
      real *8 timeinfo(10),t1,t2,omp_get_wtime

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      real *8 thresh,ra
      real *8 rr,rmin
      integer *8 nss,ii,l,npover

      integer *8 nd,ntarg0
      integer *8 ier,iper

      real *8 ttot,done,pi

      parameter (nd=1,ntarg0=1)
      data over4pi/0.07957747154594767d0/

      ns = nptso
      done = 1
      pi = atan(done)*4

c
c    estimate max number of sources in neear field of 
c    any target
c
      nmax = 0
      call get_near_corr_max(ntarg,row_ptr,nnz,col_ind,npatches,
     1  ixyzso,nmax)
      allocate(srctmp2(3,nmax),ctmp2(nmax),dtmp2(3,nmax))
           
           
      ifpgh = 0
      ifpghtarg = 1
      allocate(sources(3,ns),targvals(3,ntarg))
      allocate(charges(ns),dipvec(3,ns))
      allocate(sigmaover(ns))

c 
c       oversample density
c

      call oversample_fun_surf(int(1,8),npatches,norders,ixyzs,iptype, 
     1    npts,sigma,novers,ixyzso,ns,sigmaover)


c
c       set relevatn parameters for the fmm
c
      alpha = dpars(1)*over4pi
      beta = dpars(2)*over4pi
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

      iper = 0
      ier = 0

c
c
c       call the fmm
c

      call cpu_time(t1)
C$      t1 = omp_get_wtime()      
      call lfmm3d(nd,eps,ns,sources,ifcharge,charges,
     1  ifdipole,dipvec,iper,ifpgh,tmp,tmp,tmp,ntarg,targvals,ifpghtarg,
     1  pot,tmp,tmp,ier)
      call cpu_time(t2)
C$      t2 = omp_get_wtime()

      timeinfo(1) = t2-t1


c
c        compute threshold for ignoring local computation
c
      
      call get_fmm_thresh(int(3,8),ns,sources,int(3,8),ntarg,targvals,
     1     thresh)
      
c
c
c       add in precomputed quadrature
c

      call cpu_time(t1)
C$      t1 = omp_get_wtime()

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart)
C$OMP$PRIVATE(jstart,pottmp,npols,l)
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
          do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
            nss = nss + 1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)

            if(ifcharge.eq.1) ctmp2(nss) = charges(l)
            if(ifdipole.eq.1) then
              dtmp2(1,nss) = dipvec(1,l)
              dtmp2(2,nss) = dipvec(2,l)
              dtmp2(3,nss) = dipvec(3,l)
            endif
          enddo
        enddo

        val = 0
        if(ifcharge.eq.1.and.ifdipole.eq.0) then
          call l3ddirectcp(nd,srctmp2,ctmp2,
     1        nss,targvals(1,i),ntarg0,val,thresh)
        endif

        if(ifcharge.eq.0.and.ifdipole.eq.1) then
          call l3ddirectdp(nd,srctmp2,dtmp2,
     1          nss,targvals(1,i),ntarg0,val,thresh)
        endif

        if(ifcharge.eq.1.and.ifdipole.eq.1) then
          call l3ddirectcdp(nd,srctmp2,ctmp2,dtmp2,
     1          nss,targvals(1,i),ntarg0,val,thresh)
        endif
        pot(i) = pot(i) - val
      enddo
      
      call cpu_time(t2)
C$      t2 = omp_get_wtime()     

      timeinfo(2) = t2-t1


cc      call prin2('quadrature time=*',timeinfo,2)
      
      ttot = timeinfo(1) + timeinfo(2)

      
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
      subroutine lpcomp_lap_comb_dir_setsub(npatches,norders,ixyzs,
     1   iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2   eps,dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,sigma,novers,
     3   nptso,ixyzso,srcover,whtsover,pot)
c
c
c      this subroutine evaluates the layer potential for
c      the representation u = (\alpha S_{0} + \beta D_{0}) 
c      where the near field is precomputed and stored
c      in the row sparse compressed format.
c
c
c     The fmm is used to accelerate the far-field and 
c     near-field interactions are handled via precomputed quadrature
c
c
c       input:
c         npatches - integer *8
c            number of patches
c
c         norders- integer *8(npatches)
c            order of discretization on each patch 
c
c         ixyzs - integer *8(npatches+1)
c            ixyzs(i) denotes the starting location in srccoefs,
c               and srcvals array corresponding to patch i
c   
c         iptype - integer *8(npatches)
c            type of patch
c             iptype = 1, triangular patch discretized using RV nodes
c
c         npts - integer *8
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
c         ndtarg - integer *8
c            leading dimension of target array
c        
c         ntarg - integer *8
c            number of targets
c
c         targs - real *8 (ndtarg,ntarg)
c            target information
c
c          eps - real *8
c             precision requested
c
c          dpars - real *8 (2)
c              kernel parameters (Referring to formula (1))
c              dpars(2) = alpha
c              dpars(3) = beta
c
c           nnz - integer *8 8
c             number of source patch-> target interactions in the near field
c 
c           row_ptr - integer *8(ntarg+1)
c              row_ptr(i) is the pointer
c              to col_ind array where list of relevant source patches
c              for target i start
c
c           col_ind - integer *8 (nnz)
c               list of source patches relevant for all targets, sorted
c               by the target number
c
c           iquad - integer *8(nnz+1)
c               location in wnear array where quadrature for col_ind(i)
c               starts
c
c           nquad - integer *8
c               number of entries in wnear
c
c           wnear - real *8(nquad)
c               the near field quadrature correction
c
c           sigma - real *8(npts)
c               density for layer potential
c
c           novers - integer *8(npatches)
c              order of discretization for oversampled sources and
c               density
c
c         ixyzso - integer *8(npatches+1)
c            ixyzso(i) denotes the starting location in srcover,
c               corresponding to patch i
c   
c           nptso - integer *8
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
      integer *8 npatches,norder,npols,npts
      integer *8 ndtarg,ntarg
      integer *8 norders(npatches),ixyzs(npatches+1)
      integer *8 ixyzso(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps
      real *8 targs(ndtarg,ntarg)
      real *8 dpars(2)
      integer *8 nnz,row_ptr(ntarg+1),col_ind(nnz),nquad
      integer *8 iquad(nnz+1)
      real *8 wnear(nquad),sigma(npts)
      integer *8 novers(npatches+1)
      integer *8 nover,npolso,nptso
      real *8 srcover(12,nptso),whtsover(nptso)
      real *8 pot(ntarg)
      real *8, allocatable :: potsort(:)

      real *8, allocatable :: sources(:,:),targvals(:,:)
      real *8, allocatable :: charges(:),dipvec(:,:),sigmaover(:)
      integer *8, allocatable :: iboxtarg(:),iboxsrc(:)
      integer *8 ns,nt
      real *8 alpha,beta
      integer *8 ifcharge,ifdipole
      integer *8 ifpgh,ifpghtarg
      real *8 tmp(10),val


      integer *8 i,j,jpatch,jquadstart,jstart


      integer *8 ifaddsub

      integer *8 ltree,ipointer(8)
      integer *8, allocatable :: itree(:)
      integer *8, allocatable :: il1(:),il2(:),ilint(:),il1m2(:)
      integer *8, allocatable :: il2m1(:)
      real *8, allocatable :: boxsize(:),centers(:,:)
      integer *8, allocatable :: isrcse(:,:),isrcper(:)
      integer *8, allocatable :: itargse(:,:),itargper(:)

      integer *8, allocatable :: nlist1(:),list1(:,:)
      integer *8, allocatable :: nlist2(:),list2(:,:)
      integer *8, allocatable :: nlist3(:),list3(:,:)
      integer *8, allocatable :: nlist4(:),list4(:,:)

      real *8 expc(3)
      integer *8 ibox,nexpc,idivflag,iert,ifnear,ii,isend,isep,isource
      integer *8 isstart,itarg,itend,itstart,itt,jbox,jpt,mhung,mnbors
      integer *8 iss,l,lpt,mnlist1,mnlist2,mnlist3,mnlist4
      integer *8 n1m2,n2m1,nadd,nbmax,nboxes,nchild,ndiv,nl2,nlevels
      integer *8 nlmax,npover,nl1,ntj
      integer *8 nlmin,iper,ier,ifunif
      
      real *8 ddot,pottmp
      real *8, allocatable :: ctmp1(:),ctmp2(:),dtmp1(:,:),
     1   dtmp2(:,:)
      real *8 radexp,epsfmm

      integer *8 ipars
      complex *16 zpars
      real *8 timeinfo(10),t1,t2,omp_get_wtime,timeinfo_fmm(6)

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp1(:,:),srctmp2(:,:)
      real *8 thresh,ra
      real *8 over4pi
      integer *8 nss

      integer *8 nd,ntarg0

      real *8 ttot,done,pi

      parameter (nd=1,ntarg0=1)
      data over4pi/0.07957747154594767d0/

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

      call oversample_fun_surf(int(1,8),npatches,norders,ixyzs,iptype, 
     1    npts,sigma,novers,ixyzso,ns,sigmaover)
c
c       set relevatn parameters for the fmm
c
      alpha = dpars(1)*over4pi
      beta = dpars(2)*over4pi
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

      isep = 1
      nlmax = 51
      nbmax = 0
      nlevels = 0
      nboxes = 0
      ltree = 0

      idivflag = 0
      mnlist1 = 0
      mnlist2 = 0
      mnlist3 = 0
      mnlist4 = 0


      call lndiv(eps,ns,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg,
     1   ndiv,idivflag) 
c
cc      set tree flags
c 
       nlmax = 51
       nlevels = 0
       nboxes = 0
       ltree = 0
       nlmin = 0
       ifunif = 0
       iper = 0


c
cc     memory management code for contructing level restricted tree
      call pts_tree_mem(sources,ns,targvals,ntarg,idivflag,ndiv,
     1  nlmin,nlmax,iper,ifunif,nlevels,nboxes,ltree)
      
       allocate(itree(ltree))
       allocate(boxsize(0:nlevels))
       allocate(centers(3,nboxes))

c       Call tree code
      call pts_tree_build(sources,ns,targvals,ntarg,idivflag,ndiv,
     1  nlmin,nlmax,iper,ifunif,nlevels,nboxes,ltree,itree,ipointer,
     2  centers,boxsize)
      
      

      allocate(isrcse(2,nboxes),itargse(2,nboxes))
      allocate(isrcper(ns),itargper(ntarg))

      call pts_tree_sort(ns,sources,itree,ltree,nboxes,nlevels,
     1   ipointer,centers,isrcper,isrcse)

      call pts_tree_sort(ntarg,targvals,itree,ltree,nboxes,nlevels,
     1   ipointer,centers,itargper,itargse)

      ifnear = 0

      mnbors = 27
      isep = 1

      call computemnlists(nlevels,nboxes,itree(ipointer(1)),boxsize,
     1  centers,itree(ipointer(3)),itree(ipointer(4)),
     2  itree(ipointer(5)),isep,itree(ipointer(6)),mnbors,
     2  itree(ipointer(7)),iper,mnlist1,mnlist2,mnlist3,mnlist4)
      
      allocate(list1(mnlist1,nboxes),nlist1(nboxes))
      allocate(list2(mnlist2,nboxes),nlist2(nboxes))
      allocate(list3(mnlist3,nboxes),nlist3(nboxes))
      allocate(list4(mnlist4,nboxes),nlist4(nboxes))

      call computelists(nlevels,nboxes,itree(ipointer(1)),boxsize,
     1  centers,itree(ipointer(3)),itree(ipointer(4)),
     2  itree(ipointer(5)),isep,itree(ipointer(6)),mnbors,
     3  itree(ipointer(7)),iper,nlist1,mnlist1,list1,nlist2,
     4  mnlist2,list2,nlist3,mnlist3,list3,
     4  nlist4,mnlist4,list4)



c
c
c       call the fmm
c

      ier = 0
      iper = 0
      call cpu_time(t1)
C$      t1 = omp_get_wtime()      
      call lfmm3d_ndiv(nd,eps,ns,sources,ifcharge,charges,
     1  ifdipole,dipvec,iper,ifpgh,tmp,tmp,tmp,ntarg,targvals,ifpghtarg,
     1  pot,tmp,tmp,ndiv,idivflag,ifnear,timeinfo_fmm,ier)
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
C$OMP$PRIVATE(jstart,pottmp,npols,l)
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
      call dreorderf(int(1,8),ntarg,pot,potsort,itargper)



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
      allocate(srctmp1(3,ndiv*mnlist1))
      allocate(ctmp1(ndiv*mnlist1),dtmp1(3,ndiv*mnlist1))
      allocate(il1(ndiv*mnlist1),il1m2(ndiv*mnlist1))

  

      call cpu_time(t1)
C$      t1 = omp_get_wtime()     

C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP$PRIVATE(ibox,nchild,nl2)
C$OMP$PRIVATE(i,jbox,isstart,isend,j,isource,il2)
C$OMP$PRIVATE(itstart,itend,itt,itarg,nl1,il1,il1m2,il2m1)
C$OMP$PRIVATE(jpatch,l,jpt,lpt,n1m2,n2m1,ii,val,npover)
C$OMP$PRIVATE(ctmp1,ctmp2,dtmp1,dtmp2,srctmp1,srctmp2)
C$OMP$SCHEDULE(DYNAMIC)
      do ibox = 1,nboxes
        nchild = itree(ipointer(4)+ibox-1)
        if(nchild.eq.0) then

c
c     populate il2
c
          nl2 = 0
          do i=1,nlist1(ibox)
            jbox = list1(i,ibox) 
            isstart = isrcse(1,jbox) 
            isend = isrcse(2,jbox)
            do j=isstart,isend
              isource = isrcper(j) 
              nl2 = nl2 + 1
              il2(nl2) = isource
            enddo
          enddo


c
c    end of populating il2.
c    
c    now loop over targets in this box
c
          itstart = itargse(1,ibox) 
          itend = itargse(2,ibox)
          do itt = itstart,itend
            itarg = itargper(itt) 
            
            nl1 = 0
            do j=row_ptr(itarg),row_ptr(itarg+1)-1
              jpatch = col_ind(j)
              nl1 = nl1 + ixyzso(jpatch+1)-ixyzso(jpatch)
            enddo

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
              call l3ddirectcp(nd,srctmp1,ctmp1,
     1          n1m2,targvals(1,itarg),ntarg0,val,thresh)
            endif

            if(ifcharge.eq.0.and.ifdipole.eq.1) then
              call l3ddirectdp(nd,srctmp1,dtmp1,
     1          n1m2,targvals(1,itarg),ntarg0,val,thresh)
            endif

            if(ifcharge.eq.1.and.ifdipole.eq.1) then
              call l3ddirectcdp(nd,srctmp1,ctmp1,dtmp1,
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
              call l3ddirectcp(nd,srctmp2,ctmp2,
     1          n2m1,targvals(1,itarg),ntarg0,val,thresh)
            endif

            if(ifcharge.eq.0.and.ifdipole.eq.1) then
              call l3ddirectdp(nd,srctmp2,dtmp2,
     1          n2m1,targvals(1,itarg),ntarg0,val,thresh)
            endif

            if(ifcharge.eq.1.and.ifdipole.eq.1) then
              call l3ddirectcdp(nd,srctmp2,ctmp2,dtmp2,
     1          n2m1,targvals(1,itarg),ntarg0,val,thresh)
            endif
c
c  scatter step
c
            potsort(itt) = potsort(itt) + val

          enddo
        endif
      enddo
C$OMP END PARALLEL DO      
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      
      timeinfo(2) = t2-t1

      call dreorderi(int(1,8),ntarg,potsort,pot,itargper)

cc      call prin2('quadrature time=*',timeinfo,2)
      
      ttot = timeinfo(1) + timeinfo(2)
      
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
      subroutine lap_comb_dir_solver(npatches,norders,ixyzs,
     1    iptype,npts,srccoefs,srcvals,eps,dpars,numit,ifinout,
     2    rhs,eps_gmres,niter,errs,rres,soln)
c
c
c        this subroutine solves the Laplace dirichlet problem
c     on the interior or exterior of an object where the potential
c     is represented as a combined field integral equation.
c
c
c     Representation:
c        u = \alpha S_{0} + \beta D_{0}
c     
c     The linear system is solved iteratively using GMRES
c
c
c       input:
c         npatches - integer *8
c            number of patches
c
c         norders- integer *8(npatches)
c            order of discretization on each patch 
c
c         ixyzs - integer *8(npatches+1)
c            ixyzs(i) denotes the starting location in srccoefs,
c               and srcvals array corresponding to patch i
c   
c         iptype - integer *8(npatches)
c            type of patch
c             iptype = 1, triangular patch discretized using RV nodes
c
c         npts - integer *8
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
c          dpars - real *8 (2)
c              kernel parameters (Referring to formula (1))
c              dpars(1) = alpha
c              dpars(2) = beta
c
c          ifinout - integer *8
c              flag for interior or exterior problems (normals assumed to 
c                be pointing in exterior of region)
c              ifinout = 0, interior problem
c              ifinout = 1, exterior problem
c
c           rhs - real *8(npts)
c              right hand side
c
c           eps_gmres - real *8
c                gmres tolerance requested
c
c           numit - integer *8
c              max number of gmres iterations
c
c         output
c           niter - integer *8
c              number of gmres iterations required for relative residual
c          
c           errs(1:iter) - relative residual as a function of iteration
c              number
c 
c           rres - real *8
c              relative residual for computed solution
c              
c           soln - real *8(npts)
c              density which solves the dirichlet problem
c
c
      implicit none
      integer *8 npatches,norder,npols,npts
      integer *8 ifinout
      integer *8 norders(npatches),ixyzs(npatches+1)
      integer *8 iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,eps_gmres
      real *8 dpars(2)
      real *8 rhs(npts)
      real *8 soln(npts)

      real *8, allocatable :: targs(:,:)
      integer *8, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      integer *8 ndtarg,ntarg

      real *8 errs(numit+1)
      real *8 rres,eps2
      integer *8 niter


      integer *8 nover,npolso,nptso
      integer *8 nnz,nquad
      integer *8, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      real *8, allocatable :: wnear(:)

      real *8, allocatable :: srcover(:,:),wover(:)
      integer *8, allocatable :: ixyzso(:),novers(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:) 

      integer *8 i,j,jpatch,jquadstart,jstart

      integer *8 ipars
      complex *16 zpars
      real *8 timeinfo(10),t1,t2,omp_get_wtime


      real *8 ttot,done,pi
      real *8 rfac,rfac0
      integer *8 iptype_avg,norder_avg
      integer *8 ikerorder, iquadtype,npts_over

c
c
c       gmres variables
c
      real *8 did,dtmp
      real *8 rb,wnrm2
      integer *8 numit,it,iind,it1,k,l
      real *8 rmyerr
      real *8 temp
      real *8, allocatable :: vmat(:,:),hmat(:,:)
      real *8, allocatable :: cs(:),sn(:)
      real *8, allocatable :: svec(:),yvec(:),wtmp(:)

      complex *16 ztmp


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

      print *, "starting to generate near quadrature"
      call cpu_time(t1)
C$      t1 = omp_get_wtime()      

      call getnearquad_lap_comb_dir(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     1      ipatch_id,uvs_targ,eps,dpars,iquadtype,nnz,row_ptr,col_ind,
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
      did = -(-1)**(ifinout)*dpars(2)/2


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
C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rb)
      do i=1,npts
        rb = rb + abs(rhs(i))**2
      enddo
C$OMP END PARALLEL DO      
      rb = sqrt(rb)

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        vmat(i,1) = rhs(i)/rb
      enddo
C$OMP END PARALLEL DO      

      svec(1) = rb

      do it=1,numit
        it1 = it + 1

c
c        NOTE:
c        replace this routine by appropriate layer potential
c        evaluation routine  
c


        call lpcomp_lap_comb_dir_addsub(npatches,norders,ixyzs,
     1    iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     2    eps,dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,
     3    vmat(1,it),novers,npts_over,ixyzso,srcover,wover,wtmp)

        do k=1,it
          dtmp = 0
C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:dtmp)          
          do j=1,npts
            dtmp = dtmp + wtmp(j)*vmat(j,k)
          enddo
C$OMP END PARALLEL DO          
          hmat(k,it) = dtmp

C$OMP PARALLEL DO DEFAULT(SHARED) 
          do j=1,npts
            wtmp(j) = wtmp(j)-hmat(k,it)*vmat(j,k)
          enddo
C$OMP END PARALLEL DO          
        enddo
          
        hmat(it,it) = hmat(it,it)+did
        wnrm2 = 0
C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:wnrm2)        
        do j=1,npts
          wnrm2 = wnrm2 + abs(wtmp(j))**2
        enddo
C$OMP END PARALLEL DO        
        wnrm2 = sqrt(wnrm2)

C$OMP PARALLEL DO DEFAULT(SHARED) 
        do j=1,npts
          vmat(j,it1) = wtmp(j)/wnrm2
        enddo
C$OMP END PARALLEL DO        

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
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
          do j=1,npts
            soln(j) = 0
            do i=1,it
              soln(j) = soln(j) + yvec(i)*vmat(j,i)
            enddo
          enddo
C$OMP END PARALLEL DO          


          rres = 0
C$OMP PARALLEL DO DEFAULT(SHARED)          
          do i=1,npts
            wtmp(i) = 0
          enddo
C$OMP END PARALLEL DO          
c
c        NOTE:
c        replace this routine by appropriate layer potential
c        evaluation routine  
c


          call lpcomp_lap_comb_dir_addsub(npatches,norders,ixyzs,
     1      iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     2      eps,dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,
     3      soln,novers,npts_over,ixyzso,srcover,wover,wtmp)

            
C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rres)            
          do i=1,npts
            rres = rres + abs(did*soln(i) + wtmp(i)-rhs(i))**2
          enddo
C$OMP END PARALLEL DO          
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
c
      subroutine lap_comb_dir_fds_mem(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,eps,dpars,nifds,nrfds,nzfds)
c
cf2py   intent(in) npatches,norders,ixyzs,iptype,npts
cf2py   intent(in) srccoefs,srcvals,eps,dpars
cf2py   intent(out) nifds,nrfds,nzfds

c
c----------------------------------
c  This subroutine estimates the memory requirements
c  for the precomputation routine of the fast direct solver
c  for solving Dirichlet problems for Laplace's equation
c  using the combined field integral equation:
c
c  .. math::
c
c     u = \alpha S_{0}[\sigma] + \beta D_{0} [\sigma]
c 
c  The precomputation routine computes an integer *8 array (ifds)
c  a real array (rfds), and complex array (zfds)
c 
c  The following quantities will be computed during the 
c  precomputation phase of the fast direct solver
c
c  - ifds(1): npts_over
c  - ifds(2): nnz
c  - ifds(3): nquad
c  - ifds(4): nximat
c  - ifds(5:6+npts-1): row_ptr (for near quadrature info)
c  - ifds(6+npts:6+npts+nnz-1): col_ind
c  - ifds(6+npts+nnz:6+npts+2*nnz): iquad
c  - ifds(7+npts+2*nnz:7+npts+2*nnz+npatches-1): novers
c  - ifds(7+npts+2*nnz+npatches:7+npts+2*nnz+2*npatches): ixyzso
c  - ifds(8+npts+2*nnz+2*npatches:8+npts+2*nnz+3*npatches-1): iximat
c
c  - rfds(1:2): dpars
c  - rfds(3:12*npts_over+2): srcover
c  - rfds(12*npts_over+3:13*npts_over+2): wover
c  - rfds(13*npts_over+3:13*npts_over+nximat+2): ximats
c  - rfds(13*npts_over+nximat+3:13*npts_over+nximat+nquad+2): wnear
c
c  - zfds: unused 
c
c  Input arguments:
c  
c    - npatches: integer *8
c        number of patches
c    - norders: integer *8(npatches)
c        order of discretization of each patch
c    - ixyzs: integer *8(npatches+1)
c        starting location of points on patch i
c    - iptype: integer *8(npatches)
c        type of patch
c        iptype = 1, triangle discretized using RV nodes
c    - npts: integer *8
c        total number of points on the surface
c    - srccoefs: double precision (9,npts)
c        koornwinder expansion coefs of geometry info
c    - srcvals: double precision (12,npts)
c        xyz, dxyz/du,dxyz/dv, normals at all nodes
c    - eps: double precision
c        precision requested
c    - dpars: double precision(2)
c        weights for single and double layer potential in representation
c        dpars(1) = \alpha
c        dpars(2) = \beta
c
c  Output arguments:
c  
c    - nifds: integer *8
c        size of integer *8 array
c    - nrfds: integer *8
c        size of double precision array
c    - nzfds: integer *8
c        size of double complex array
c     
      implicit none
      integer *8, intent(in) :: npatches,norders(npatches)
      integer *8, intent(in) :: ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches),npts
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
      real *8, intent(in) :: eps
      real *8, intent(in) :: dpars(2)
      integer *8, intent(out) :: nifds,nrfds,nzfds

      integer *8 nnz
      real *8, allocatable :: targs(:,:)
      integer *8 iptype_avg,norder_avg
      integer *8 ntarg,ndtarg,ikerorder
      real *8 rfac,rfac0
      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)
      integer *8, allocatable :: iquad(:),row_ptr(:),col_ind(:)
      integer *8, allocatable :: novers(:),ixyzso(:)
      integer *8 npts_over,nquad,nximat
      

      complex *16 ztmp
      integer *8 i


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
      call findnearmem(cms,npatches,rad_near,ndtarg,targs,npts,nnz)

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

      ztmp = 0
      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,
     1    rads,npts,srccoefs,ndtarg,npts,targs,ikerorder,ztmp,
     2    nnz,row_ptr,col_ind,rfac,novers,ixyzso)


      

      npts_over = ixyzso(npatches+1)-1

      nquad = iquad(nnz+1)-1

      call get_nximat(npatches,ixyzs,ixyzso,nximat)

      nifds = 7+npts+2*nnz+3*npatches
      nrfds = 2+13*npts_over+nximat+nquad
      nzfds = 0
      

      return
      end
c
c
c
c
c

      subroutine lap_comb_dir_fds_init(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,eps,dpars,nifds,ifds,nrfds,
     2     rfds,nzfds,zfds)
cf2py   intent(in) npatches,norders,ixyzs,iptype,npts
cf2py   intent(in) srccoefs,srcvals,eps,dpars
cf2py   intent(in) nifds,nrfds,nzfds
cf2py   intent(out) ifds,rfds,zfds

c
c----------------------------------
c  This subroutine precomputes a few arrays
c  for subsequent calls to the fast direct solver
c  for solving Dirichlet problems for Laplace's equation
c  using the combined field integral equation:
c
c  .. math::
c
c     u = \alpha S_{0}[\sigma] + \beta D_{0} [\sigma]
c 
c  The precomputation routine computes an integer *8 array (ifds)
c  a real array (rfds), and complex array (zfds)
c 
c  The following quantities will be computed during the 
c  precomputation phase of the fast direct solver
c
c  - ifds(1): npts_over
c  - ifds(2): nnz
c  - ifds(3): nquad
c  - ifds(4): nximat
c  - ifds(5:6+npts-1): row_ptr (for near quadrature info)
c  - ifds(6+npts:6+npts+nnz-1): col_ind
c  - ifds(6+npts+nnz:6+npts+2*nnz): iquad
c  - ifds(7+npts+2*nnz:7+npts+2*nnz+npatches-1): novers
c  - ifds(7+npts+2*nnz+npatches:7+npts+2*nnz+2*npatches): ixyzso
c  - ifds(8+npts+2*nnz+2*npatches:8+npts+2*nnz+3*npatches-1): iximat
c
c  - rfds(1:2): dpars
c  - rfds(3:12*npts_over+2): srcover
c  - rfds(12*npts_over+3:13*npts_over+2): wover
c  - rfds(13*npts_over+3:13*npts_over+nximat+2): ximats
c  - rfds(13*npts_over+nximat+3:13*npts_over+nximat+nquad+2): wnear
c
c  - zfds: unused 
c
c  Input arguments:
c  
c    - npatches: integer *8
c        number of patches
c    - norders: integer *8(npatches)
c        order of discretization of each patch
c    - ixyzs: integer *8(npatches+1)
c        starting location of points on patch i
c    - iptype: integer *8(npatches)
c        type of patch
c        iptype = 1, triangle discretized using RV nodes
c    - npts: integer *8
c        total number of points on the surface
c    - srccoefs: double precision (9,npts)
c        koornwinder expansion coefs of geometry info
c    - srcvals: double precision (12,npts)
c        xyz, dxyz/du,dxyz/dv, normals at all nodes
c    - eps: double precision
c        precision requested
c    - dpars: double precision(2)
c        weights for single and double layer potential in representation
c        dpars(1) = \alpha
c        dpars(2) = \beta
c    - nifds: integer *8
c        size of integer *8 array
c    - nrfds: integer *8
c        size of double precision array
c    - nzfds: integer *8
c        size of double complex array
c
c  Output arguments:
c
c    - ifds: integer *8(nifds)
c        precomputed integer *8 array
c    - rfds: double precision(nifds)
c        precomputed double precision array
c    - zfds: double complex(nifds)
c        precomputed double complex array
c  
c     
      implicit none
      integer *8, intent(in) :: npatches,norders(npatches)
      integer *8, intent(in) :: ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches),npts
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
      real *8, intent(in) :: eps
      real *8, intent(in) :: dpars(2)
      integer *8, intent(in) :: nifds,nrfds,nzfds
      integer *8, intent(out) :: ifds(nifds)
      real *8, intent(out) :: rfds(nrfds)
      complex *16, intent(out) :: zfds(nzfds)

      integer *8 nnz
      complex *16 ztmp
      real *8, allocatable :: targs(:,:)
      real *8, allocatable :: uvs_targ(:,:)
      integer *8, allocatable :: ipatch_id(:)
      integer *8 iptype_avg,norder_avg
      integer *8 ntarg,ndtarg,ikerorder
      real *8 rfac,rfac0
      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)
      integer *8 npts_over,nquad

      integer *8 iquadtype,istart,iend,i

      integer *8 irow_ptr,icol_ind,iiquad,inovers,iixyzso,iximat
      integer *8 nximat

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
      call findnearmem(cms,npatches,rad_near,ndtarg,targs,npts,nnz)

      irow_ptr = 5
      icol_ind = irow_ptr+npts+1
      iiquad = icol_ind+nnz
      inovers = iiquad + nnz+1
      iixyzso = inovers+npatches
      iximat = iixyzso+npatches+1

      
      call findnear(cms,npatches,rad_near,ndtarg,targs,npts,
     1   ifds(irow_ptr),ifds(icol_ind))

      call get_iquad_rsc(npatches,ixyzs,npts,nnz,ifds(irow_ptr),
     1   ifds(icol_ind),ifds(iiquad))


      ikerorder = -1
      if(abs(dpars(2)).gt.1.0d-16) ikerorder = 0


c
c    estimate oversampling for far-field, and oversample geometry
c

      ztmp = 0.0d0
      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,
     1    rads,npts,srccoefs,ndtarg,npts,targs,ikerorder,ztmp,
     2    nnz,ifds(irow_ptr),ifds(icol_ind),rfac,ifds(inovers),
     3    ifds(iixyzso))

      npts_over = ifds(iixyzso+npatches)-1
      nquad = ifds(iiquad+nnz)-1
      
      rfds(1) = dpars(1)
      rfds(2) = dpars(2)

      call oversample_geom(npatches,norders,ixyzs,iptype,npts, 
     1   srccoefs,srcvals,ifds(inovers),ifds(iixyzso),npts_over,
     2   rfds(3))

      call get_qwts(npatches,ifds(inovers),ifds(iixyzso),iptype,
     1      npts_over,rfds(3),rfds(12*npts_over+3))

      

      ifds(1) = npts_over
      ifds(2) = nnz
      ifds(3) = nquad

      call get_nximat(npatches,ixyzs,ifds(iixyzso),nximat)

      ifds(4) = nximat

      istart = 13*npts_over+3

      call get_ximats(npatches,iptype,norders,ixyzs,ifds(inovers),
     1  ifds(iixyzso),nximat,rfds(istart),ifds(iximat))


c
c    initialize patch_id and uv_targ for on surface targets
c
      allocate(ipatch_id(ntarg),uvs_targ(2,ntarg))
      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,
     1  ipatch_id,uvs_targ)

      iquadtype = 1
      istart = 13*npts_over + nximat + 3

      call getnearquad_lap_comb_dir(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     1      ipatch_id,uvs_targ,eps,dpars,iquadtype,nnz,ifds(irow_ptr),
     1      ifds(icol_ind),ifds(iiquad),rfac0,nquad,rfds(istart))
      
      return
      end
c
c
c
c
c
c

      subroutine lap_comb_dir_fds_matgen(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,eps,dpars,nifds,ifds,nrfds,
     2     rfds,nzfds,zfds,nent_csc,col_ptr,row_ind,dmatent)
cf2py   intent(in) npatches,norders,ixyzs,iptype,npts
cf2py   intent(in) srccoefs,srcvals,eps,dpars
cf2py   intent(in) nifds,nrfds,nzfds
cf2py   intent(in) ifds,rfds,zfds
cf2py   intent(in) nent_csc,col_ptr,row_ind
cf2py   intent(out) dmatent
   
c----------------------------------
c  This subroutine generates matrix entries queried 
c  for solving Dirichlet problems for Laplace's equation
c  using the combined field integral equation:
c
c  .. math::
c
c     u = \alpha S_{0}[\sigma] + \beta D_{0} [\sigma]
c 
c  The precomputation routine computes an integer *8 array (ifds)
c  a real array (rfds), and complex array (zfds)
c 
c  The following quantities will be computed during the 
c  precomputation phase of the fast direct solver
c
c  - ifds(1): npts_over
c  - ifds(2): nnz
c  - ifds(3): nquad
c  - ifds(4): nximat
c  - ifds(5:6+npts-1): row_ptr (for near quadrature info)
c  - ifds(6+npts:6+npts+nnz-1): col_ind
c  - ifds(6+npts+nnz:6+npts+2*nnz): iquad
c  - ifds(7+npts+2*nnz:7+npts+2*nnz+npatches-1): novers
c  - ifds(7+npts+2*nnz+npatches:7+npts+2*nnz+2*npatches): ixyzso
c  - ifds(8+npts+2*nnz+2*npatches:8+npts+2*nnz+3*npatches-1): iximat
c
c  - rfds(1:2): dpars
c  - rfds(3:12*npts_over+2): srcover
c  - rfds(12*npts_over+3:13*npts_over+2): wover
c  - rfds(13*npts_over+3:13*npts_over+nximat+2): ximats
c  - rfds(13*npts_over+nximat+3:13*npts_over+nximat+nquad+2): wnear
c
c  - zfds: unused
c
c  The list of input entries requested must be provided in 
c  column sparse compressed format. Use the conv_to_csc
c  utility to convert a list of (i,j) entries to column sparse
c  compressed representation.
c
c  Note: if double layer potential is used in representation,
c  then only principal value part of the layer potential
c  is returned. The term corresponding to the jump in layer potential
c  will have to be added manually.
c
c
c  Input arguments:
c  
c    - npatches: integer *8
c        number of patches
c    - norders: integer *8(npatches)
c        order of discretization of each patch
c    - ixyzs: integer *8(npatches+1)
c        starting location of points on patch i
c    - iptype: integer *8(npatches)
c        type of patch
c        iptype = 1, triangle discretized using RV nodes
c    - npts: integer *8
c        total number of points on the surface
c    - srccoefs: double precision (9,npts)
c        koornwinder expansion coefs of geometry info
c    - srcvals: double precision (12,npts)
c        xyz, dxyz/du,dxyz/dv, normals at all nodes
c    - eps: double precision
c        precision requested
c    - dpars: double precision(2)
c        weights for single and double layer potential in representation
c        dpars(1) = \alpha
c        dpars(2) = \beta
c    - nifds: integer *8
c        size of integer *8 array
c    - ifds: integer *8(nifds)
c        precomputed integer *8 array
c    - nrfds: integer *8
c        size of double precision array
c    - rfds: double precision(nifds)
c        precomputed double precision array
c    - nzfds: integer *8
c        size of double complex array
c    - zfds: double complex(nifds)
c        precomputed double complex array
c    - nent_csc: integer *8
c        number of matrix entries requested
c    - col_ptr: integer *8(npts+1)
c        indicates location in row_ind where relevant entries for
c        column i begin
c    - row_ind: integer *8(nent_csc)
c        list of row indices. (row_ind(col_ptr(i):col_ptr(i+1)-1),i)
c        are all the matrix entries requested
c
c  Output arguments:
c
c    - dmatent: double precision(npts)
c        matrix entries requested 
c
c  
c     

      implicit none
      integer *8, intent(in) :: npatches,norders(npatches)
      integer *8, intent(in) :: ixyzs(npatches+1)
      integer *8, intent(in) :: iptype(npatches),npts
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
      real *8, intent(in) :: eps,dpars(2)
      integer *8, intent(in) :: nifds,nrfds,nzfds
      integer *8, intent(in) :: ifds(nifds)
      real *8, intent(in) :: rfds(nrfds)
      complex *16, intent(in) :: zfds(nzfds)
      integer *8, intent(in) :: nent_csc,col_ptr(npts+1)
      integer *8, intent(in) :: row_ind(nent_csc)
      real *8, intent(out) :: dmatent(nent_csc)

c
c        temporary variables
c
      integer *8 i,j,k,l,ipatch,npols
      integer *8 ndd,ndz,ndi
      integer *8 ilstart,ilend,istart,iend
      integer *8 irow_ptr,icol_ind,iiquad,inovers,iixyzso
      integer *8 npts_over,nnz,nquad

      integer *8, allocatable :: iuni(:),iuniind(:)
      integer *8, allocatable :: col_ptr_src(:),row_ind_src(:),iper(:)
      integer *8, allocatable :: aintb(:),iaintba(:),aintbc(:)
      integer *8, allocatable :: iaintbc(:)
      real *8, allocatable :: wquadn(:,:),wquad(:,:)
      real *8, allocatable :: wquadf(:,:),wquadf2(:,:)


      real *8 alpha, beta
      complex *16 zk,zpars
      real *8, allocatable :: srcover(:,:),wtsover(:)
      integer *8 ipars,ipt,iquad,itind,iximat,ixist,j2,jind
      integer *8 jind0,juniind,n2,naintb,naintbc,nmax,nn,npolso
      integer *8 nuni,nximat,iwnear
      integer *8, allocatable :: iaintbb(:)
      integer *8 ifcharge,ifdipole

      real *8 dzero,done

      procedure (), pointer :: fker
      external l3d_slp, l3d_dlp, l3d_comb

      done = 1
      dzero = 0

c
c
c        initialize the appropriate kernel function
c
 
      
      alpha = rfds(1)
      beta = rfds(2)

      fker => l3d_comb
      ndd = 2
      ifcharge = 1
      ifdipole = 1
      if(abs(alpha).ge.1.0d-16.and.abs(beta).lt.1.0d-16) then
        fker=>l3d_slp
        ndd = 0
        ifdipole = 0
      else if(abs(alpha).lt.1.0d-16.and.abs(beta).ge.1.0d-16) then
        fker=>l3d_dlp
        ndd = 0
        ifcharge = 0
      endif


      ndi = 0
      ndz = 0

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

      iwnear = 13*npts_over + nximat+2 
     

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

        nuni = 0
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
           wquad(iaintba(i),:) = rfds(iwnear+iquad:
     1        iwnear+iquad+npols-1)
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
            srcover(j,i) = rfds(2+12*(istart+i-2)+j)
          enddo
          wtsover(i) = rfds(2+12*npts_over+istart+i-1)
        enddo

cc        call prin2('srcover=*',srcover,12*npolso)
cc        call prin2('wtsover=*',wtsover,npolso)


cc        call prin2('zfds=*',zfds,6)

        do i=1,naintbc
          itind = aintbc(i)
          do j=1,npolso
            call fker(srcover(1,j),int(3,8),srcvals(1,itind),ndd,dpars,
     1      ndz,zfds,ndi,ipars,wquadf(j,i))
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
        ixist = ifds(iximat+ipatch-1) + 13*npts_over+2
        call dgemm_guru('t','n',npols,naintbc,npolso,done,rfds(ixist),
     1    npolso,wquadf,npolso,dzero,wquadf2,npols)
        
        do i=1,naintbc
          wquad(iaintbc(i),:) = wquadf2(:,i)
        enddo
       
        do i = 1,npols
          ipt = ixyzs(ipatch) + i-1
          do j=col_ptr(ipt),col_ptr(ipt+1)-1
             jind = row_ind(j)

             j2 = j-col_ptr(ixyzs(ipatch))+1
             juniind = iuniind(j2)
             dmatent(j) = wquad(juniind,i)
          enddo
        enddo

        deallocate(iuni,iuniind,aintb,iaintba,iaintbb,aintbc,iaintbc)
        deallocate(wquad,wquadf,wquadf2)
      enddo
      
      
      


      return
      end
