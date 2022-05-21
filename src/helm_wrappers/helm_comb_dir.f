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
c       getnearquadsub_helm_comb_dir - generates the near
c        field fmm substraction for the Dirichlet data
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
c       helm_comb_dir_solver_memest - memory estimation code
c         for determining how much memory will be used by
c         the solver (including the memory used by the fmm
c         code per iterate)
c
c
c       There are two sets of fast direct solver routines
c       and neither of them currently used oversampling.
c       The fast direct solver routines are currently in beta
c       mode
c
c       helm_comb_dir_fds_csc_mem - get memory requirements for initialization
c          routine for subsequent calls to fast direct solver
c
c       helm_comb_dir_fds_csc_init - initialize various arrays to be later 
c          used for fast direct solver
c
c       helm_comb_dir_fds_csc_matgen - query entries of the combined field
c          representation matrix (input indices must be in column 
c          sparse compressed format, and must be preceeded by a call
c          to helm_comb_dir_fds_init)
c
c       helm_comb_dir_fds_block_mem,
c       helm_comb_dir_fds_block_init,
c       helm_comb_dir_fds_block_matgen
c           (routines analogous to the three routines above but which
c             do not include any oversampling, also on input
c             are a collection of row indices and column indices and on
c             output the routine returns a complex matrix 
c             of nrows \times ncols)
c
c       REMARK:
c         The fast direct solver routines need to be optimized
c         for performance. Improvements coming shortly
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
c       for the representation u = (\alpha S_{k}  + \beta D_{k}) ---(1)
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
      subroutine getnearquadsub_helm_comb_dir(npatches,norders,
     1   ixyzso,iptype,npts_over,srcover,wover,ndtarg,ntarg,targs,
     2   ipatch_id,uvs_targ,thresh,zpars,iquadtype,nnz,row_ptr,col_ind,
     3   iquad,nquad,wnear)
c
c       this subroutine generates the near field fmm substraction
c       for the representation u = (\alpha S_{k}  + \beta D_{k}) ---(1)
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
      integer, intent(in) :: npatches,norders(npatches),npts_over,nquad
      integer, intent(in) :: ixyzso(npatches+1),iptype(npatches)
      real *8, intent(in) :: srcover(12,npts_over),thresh
      real *8, intent(in) :: wover(npts_over)
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


        call zgetnearquadsub_guru(npatches,norders,ixyzso,
     1     iptype,npts_over,srcover,wover,ndtarg,ntarg,targs,
     1     ipatch_id,uvs_targ,thresh,ipv,fker,ndd,dpars,ndz,zpars,
     1     ndi,ipars,nnz,row_ptr,col_ind,iquad,nquad,wnear)
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
      subroutine getnearquadsub_helm_comb_dir_low_mem(npatches,norders,
     1   ixyzso,iptype,nptso,srcover,wover,ndtarg,ntarg,targs,
     2   ipatch_id,uvs_targ,thresh,zpars,iquadtype,nnz,row_ptr,col_ind,
     3   iquad,nquad,wnear,ixyzs,novers,npts)
c
c       this subroutine generates the near field fmm substraction
c       for the representation u = (\alpha S_{k}  + \beta D_{k}) ---(1)
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
      integer, intent(in) :: npatches,norders(npatches),nptso,nquad
      integer, intent(in) :: novers(npatches),npts
      integer, intent(in) :: ixyzso(npatches+1),iptype(npatches)
      integer, intent(in) :: ixyzs(npatches+1)
      real *8, intent(in) :: srcover(12,nptso),thresh
      real *8, intent(in) :: wover(nptso)
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
      real *8 dpars,rr,threshsq

      complex *16 alpha,beta,coef,val
      integer i,j,ipatch,iquadstart,istart,l
      integer ipv,npols,npolso
      integer nthd,ithd,nmax,nmaxo
      integer omp_get_max_threads,omp_get_thread_num
      complex *16, allocatable :: wnear_sub(:,:)
      complex *16, allocatable :: wnear_omp(:,:)

      procedure (), pointer :: fker
      external h3d_slp, h3d_dlp, h3d_comb

c
c
c        initialize the appropriate kernel function
c
      nthd = 1
C$    nthd=omp_get_max_threads()
      alpha = zpars(2)
      beta = zpars(3)
      threshsq = thresh**2
      nmax = 0
      nmaxo = 0
      do i = 1,npatches
        if(norders(i).gt.nmax) nmax = norders(i)
        if(novers(i).gt.nmaxo) nmaxo = novers(i)
      enddo
      nmax = (nmax+1)*(nmax+2)/2
      nmaxo = (nmaxo+1)*(nmaxo+2)/2
      allocate(wnear_sub(nmaxo,nthd),wnear_omp(nmax,nthd))

      if(abs(alpha).ge.1.0d-16.and.abs(beta).lt.1.0d-16) then
        coef = alpha
      else if(abs(alpha).lt.1.0d-16.and.abs(beta).ge.1.0d-16) then
        coef = beta
      endif

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

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,ipatch)
C$OMP$PRIVATE(istart,npols,npolso,iquadstart,l,val,rr,ithd)
        do i=1,ntarg
          ithd = 0
C$        ithd = omp_get_thread_num()
          ithd = ithd + 1
          do j=row_ptr(i),row_ptr(i+1)-1
            ipatch = col_ind(j)
            istart = ixyzso(ipatch)
            npolso = ixyzso(ipatch+1)-ixyzso(ipatch)
            npols = ixyzs(ipatch+1)-ixyzs(ipatch)
            iquadstart = iquad(j)
            do l=1,npolso
              wnear_sub(l, ithd) = 0
              rr = (srcover(1,l+istart-1) - targs(1,i))**2 +
     1             (srcover(2,l+istart-1) - targs(2,i))**2 +
     2             (srcover(3,l+istart-1) - targs(3,i))**2
              if(rr.gt.threshsq) then
                call fker(srcover(1,l+istart-1),ndtarg,targs(1,i),
     1                    ndd,dpars,ndz,zpars,ndi,ipars,val)
                wnear_sub(l, ithd) = coef * val * wover(l+istart-1)
              endif
            enddo
            call oversample_fun_tri_transpose(2,
     1           norders(ipatch),npols,
     2           wnear_sub(1,ithd),
     3           novers(ipatch),npolso,wnear_omp(1,ithd))
            do l=1,npols
              wnear(iquadstart+l-1) = wnear(iquadstart+l-1) -
     1                                wnear_omp(l,ithd)
            enddo
          enddo
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
c      u = (\alpha \mathcal{S}_{k} + \beta \mathcal{D}_{k})
c
c  Note: For targets on the boundary, this routine only computes
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


      integer nptso,nnz,nquad,nquadsub


      integer nover,npolso
      integer norder,npols
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      integer, allocatable :: iquadsub(:)
      complex *16, allocatable :: wnear(:),wnearsub(:)

      real *8, allocatable :: srcover(:,:),wover(:)
      integer, allocatable :: ixyzso(:),novers(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:) 

      integer i,j,jpatch,jquadstart,jstart

      integer ipars
      real *8 dpars,timeinfo(10),t1,t2,omp_get_wtime


      real *8 ttot,done,pi
      real *8 rfac,rfac0
      real *8 over4pi
      integer iptype_avg,norder_avg
      integer ikerorder, iquadtype,npts_over
      data over4pi/0.07957747154594767d0/


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

      ikerorder = -1
      if(abs(zpars(3)).gt.1.0d-16) ikerorder = 0

c
c    estimate oversampling for far-field, and oversample geometry
c

      allocate(novers(npatches),ixyzso(npatches+1))

      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,
     1     rads,npts,srccoefs,ndtarg,ntarg,targs,ikerorder,zpars(1),
     2     nnz,row_ptr,col_ind,rfac,novers,ixyzso)

     
      allocate(iquadsub(nnz+1))
      call get_iquad_rsc(npatches,ixyzso,ntarg,nnz,row_ptr,col_ind,
     1     iquadsub)


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

      nquadsub = iquadsub(nnz+1)-1
      allocate(wnearsub(nquadsub))
      
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,nquadsub
        wnearsub(i) = 0
      enddo
C$OMP END PARALLEL DO    

      iquadtype = 1

      call getnearquad_helm_comb_dir(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     1      ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,
     1      iquad,rfac0,nquad,wnear)

c      call getnearquadsub_helm_comb_dir(npatches,norders,
c     1      ixyzso,iptype,npts_over,srcover,wover,ndtarg,ntarg,targs,
c     1      ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,
c     1      iquadsub,nquadsub,wnearsub)

      call getnearquadsub_helm_comb_dir_low_mem(npatches,norders,
     1   ixyzso,iptype,npts_over,srcover,wover,ndtarg,ntarg,targs,
     2   ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,
     3   iquad,nquad,wnear,ixyzs,novers,npts)
c
c
c   compute layer potential
c
      call lpcomp_helm_comb_dir_addsub2(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2  eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,iquadsub,nquadsub,
     3  wnear,wnearsub,sigma,novers,npts_over,ixyzso,srcover,
     4  wover,pot)

c      call lpcomp_helm_comb_dir_addsub(npatches,norders,ixyzs,
c     1  iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
c     2  eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,
c     3  sigma,novers,npts_over,ixyzso,srcover,wover,pot)



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
c      the representation u = (\alpha S_{k} + \beta D_{k}) 
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
      real *8 radexp,epsfmm

      integer ipars
      real *8 dpars,timeinfo(10),t1,t2,omp_get_wtime

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      complex *16, allocatable :: ctmp2(:),dtmp2(:,:)
      real *8 thresh,ra
      real *8 rr,rmin
      real *8 over4pi
      integer nss,ii,l,npover
      integer nmax,ier,iper

      integer nd,ntarg0

      real *8 ttot,done,pi
      data over4pi/0.07957747154594767d0/

      parameter (nd=1,ntarg0=1)

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

      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype, 
     1    npts,sigma,novers,ixyzso,ns,sigmaover)


      ra = 0


c
c       set relevatn parameters for the fmm
c
      alpha = zpars(2)*over4pi
      beta = zpars(3)*over4pi
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
     1  ifdipole,dipvec,iper,ifpgh,tmp,tmp,tmp,ntarg,targvals,ifpghtarg,
     1  pot,tmp,tmp,ier)
      call cpu_time(t2)
C$      t2 = omp_get_wtime()

      timeinfo(1) = t2-t1


c
c        compute threshold for ignoring local computation
c
      call get_fmm_thresh(3,ns,sources,3,ntarg,targvals,thresh)

c
c       add in precomputed quadrature

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
            nss = nss+1
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
      subroutine lpcomp_helm_comb_dir_addsub2(npatches,norders,ixyzs,
     1   iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2   eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,iquadsub,nquadsub,
     3   wnear,wnearsub,sigma,novers,nptso,ixyzso,srcover,whtsover,pot)
c
c
c      this subroutine evaluates the layer potential for
c      the representation u = (\alpha S_{k} + \beta D_{k}) 
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
      integer, intent(in) :: nquadsub
      integer, intent(in) :: iquad(nnz+1),iquadsub(nnz+1)
      complex *16, intent(in) :: wnear(nquad),sigma(npts)
      complex *16, intent(in) :: wnearsub(nquadsub)
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
      real *8 radexp,epsfmm

      integer ipars
      real *8 dpars,timeinfo(10),t1,t2,omp_get_wtime

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      complex *16, allocatable :: ctmp2(:),dtmp2(:,:)
      real *8 thresh,ra
      real *8 rr,rmin
      real *8 over4pi
      integer nss,ii,l,npover
      integer nmax,ier,iper

      integer nd,ntarg0

      real *8 ttot,done,pi
      data over4pi/0.07957747154594767d0/

      parameter (nd=1,ntarg0=1)

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

      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype, 
     1    npts,sigma,novers,ixyzso,ns,sigmaover)


      ra = 0


c
c       set relevatn parameters for the fmm
c
      alpha = zpars(2)*over4pi
      beta = zpars(3)*over4pi
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
     1  ifdipole,dipvec,iper,ifpgh,tmp,tmp,tmp,ntarg,targvals,ifpghtarg,
     1  pot,tmp,tmp,ier)
      call cpu_time(t2)
C$      t2 = omp_get_wtime()

      timeinfo(1) = t2-t1


c
c        compute threshold for ignoring local computation
c
      call get_fmm_thresh(3,ns,sources,3,ntarg,targvals,thresh)

c
c       add in precomputed quadrature

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
          npols = ixyzso(jpatch+1)-ixyzso(jpatch)
          jquadstart = iquadsub(j)
          jstart = ixyzso(jpatch)
          do l=1,npols
            pot(i) = pot(i) -
     1               0*wnearsub(jquadstart+l-1)*sigmaover(jstart+l-1)
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO


c

cC$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
cC$OMP$PRIVATE(ctmp2,dtmp2,nss,l,jstart,ii,val,npover)
c      do i=1,ntarg
c        nss = 0
c        do j=row_ptr(i),row_ptr(i+1)-1
c          jpatch = col_ind(j)
c          do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
c            nss = nss+1
c            srctmp2(1,nss) = srcover(1,l)
c            srctmp2(2,nss) = srcover(2,l)
c            srctmp2(3,nss) = srcover(3,l)
c
c            if(ifcharge.eq.1) ctmp2(nss) = charges(l)
c            if(ifdipole.eq.1) then
c              dtmp2(1,nss) = dipvec(1,l)
c              dtmp2(2,nss) = dipvec(2,l)
c              dtmp2(3,nss) = dipvec(3,l)
c            endif
c          enddo
c        enddo
c
c        val = 0
c        if(ifcharge.eq.1.and.ifdipole.eq.0) then
c          call h3ddirectcp(nd,zpars(1),srctmp2,ctmp2,
c     1        nss,targvals(1,i),ntarg0,val,thresh)
c        endif
c
c        if(ifcharge.eq.0.and.ifdipole.eq.1) then
c          call h3ddirectdp(nd,zpars(1),srctmp2,dtmp2,
c     1          nss,targvals(1,i),ntarg0,val,thresh)
c        endif
c
c        if(ifcharge.eq.1.and.ifdipole.eq.1) then
c          call h3ddirectcdp(nd,zpars(1),srctmp2,ctmp2,dtmp2,
c     1          nss,targvals(1,i),ntarg0,val,thresh)
c        endif
c        pot(i) = pot(i) - val
c      enddo
      
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
c      the representation u = (\alpha S_{k} + \beta D_{k})  
c      where the near field is precomputed and stored
c      in the row sparse compressed format.
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

      integer *8 ltree,ipointer(8)
      integer, allocatable :: itree(:)
      integer, allocatable :: il1(:),il2(:),ilint(:),il1m2(:),il2m1(:)
      real *8, allocatable :: boxsize(:),centers(:,:)
      integer, allocatable :: isrcse(:,:),isrcper(:)
      integer, allocatable :: itargse(:,:),itargper(:)

      real *8 expc(3)
      integer ibox,nexpc,idivflag,iert,ifnear,ii,isend,isep,isource
      integer isstart,itarg,itend,itstart,itt,jbox,jpt,mhung,mnbors
      integer iss,l,lpt,mnlist1,mnlist2,mnlist3,mnlist4
      integer n1m2,n2m1,nadd,nbmax,nboxes,nchild,ndiv,nl2,nlevels
      integer nlmax,npover,nl1,ntj
      integer ier,nlmin,iper,ifunif

      integer, allocatable :: nlist1(:),list1(:,:)
      integer, allocatable :: nlist2(:),list2(:,:)
      integer, allocatable :: nlist3(:),list3(:,:)
      integer, allocatable :: nlist4(:),list4(:,:)
      
      complex *16 zdotu,pottmp
      complex *16, allocatable :: ctmp1(:),ctmp2(:),dtmp1(:,:),
     1   dtmp2(:,:)
      real *8 radexp,epsfmm

      integer ipars
      real *8 dpars,timeinfo(10),t1,t2,omp_get_wtime
      real *8 timeinfo_fmm(6)

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp1(:,:),srctmp2(:,:)
      real *8 thresh,ra
      real *8 over4pi
      integer nss

      integer nd,ntarg0

      real *8 ttot,done,pi
      data over4pi/0.07957747154594767d0/

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


c
c       set relevatn parameters for the fmm
c
      alpha = zpars(2)*over4pi
      beta = zpars(3)*over4pi
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

      call hndiv(eps,ns,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg,
     1   ndiv,idivflag) 
c
cc      set tree flags
c 
       nlmax = 51
       nlevels = 0
       nboxes = 0
       ltree = 0
       nlmin = 0 
       iper = 0
       ifunif = 0

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

      call cpu_time(t1)
C$      t1 = omp_get_wtime()      
      call hfmm3d_ndiv(nd,eps,zpars(1),ns,sources,ifcharge,charges,
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
      call dreorderf(2,ntarg,pot,potsort,itargper)



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

          enddo
        endif
      enddo
C$OMP END PARALLEL DO      
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      
      timeinfo(2) = t2-t1

      call dreorderi(2,ntarg,potsort,pot,itargper)

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
c     until a relative residual of eps_gmres is reached
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
      integer nnz,nquad,nquadsub
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      integer, allocatable :: iquadsub(:)
      complex *16, allocatable :: wnear(:),wnearsub(:)

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
      call findnearmem(cms,npatches,rad_near,ndtarg,targs,npts,nnz)
      print *, "nnz=",nnz

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,ndtarg,targs,npts,row_ptr, 
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

      allocate(iquadsub(nnz+1))
      call get_iquad_rsc(npatches,ixyzso,ntarg,nnz,row_ptr,col_ind,
     1     iquadsub)

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

      nquadsub = iquadsub(nnz+1)-1
      allocate(wnearsub(nquadsub))
      
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,nquadsub
        wnearsub(i) = 0
      enddo
C$OMP END PARALLEL DO    

      iquadtype = 1

      print *, "starting to generate near quadrature"
      call cpu_time(t1)
C$      t1 = omp_get_wtime()      

      call getnearquad_helm_comb_dir(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     1      ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,
     1      iquad,rfac0,nquad,wnear)

c      call getnearquadsub_helm_comb_dir(npatches,norders,
c     1      ixyzso,iptype,npts_over,srcover,wover,ndtarg,npts,targs,
c     1      ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,
c     1      iquadsub,nquadsub,wnearsub)
      call getnearquadsub_helm_comb_dir_low_mem(npatches,norders,
     1   ixyzso,iptype,npts_over,srcover,wover,ndtarg,ntarg,targs,
     2   ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,
     3   iquad,nquad,wnear,ixyzs,novers,npts)

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
      zid = -(-1)**(ifinout)*zpars(3)/2


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

        call lpcomp_helm_comb_dir_addsub2(npatches,norders,ixyzs,
     1    iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     2    eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,iquadsub,nquadsub,
     3    wnear,wnearsub,vmat(1,it),novers,npts_over,ixyzso,srcover,
     4    wover,wtmp)

c        call lpcomp_helm_comb_dir_addsub(npatches,norders,ixyzs,
c     1    iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
c     2    eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,
c     3    vmat(1,it),novers,npts_over,ixyzso,srcover,wover,wtmp)

        do k=1,it
          ztmp = 0
C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:ztmp)          
          do j=1,npts
            ztmp = ztmp + wtmp(j)*conjg(vmat(j,k))
          enddo
C$OMP END PARALLEL DO          
          hmat(k,it) = ztmp

C$OMP PARALLEL DO DEFAULT(SHARED) 
          do j=1,npts
            wtmp(j) = wtmp(j)-hmat(k,it)*vmat(j,k)
          enddo
C$OMP END PARALLEL DO          
        enddo
          
        hmat(it,it) = hmat(it,it)+zid
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

          call lpcomp_helm_comb_dir_addsub2(npatches,norders,ixyzs,
     1      iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     2      eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,iquadsub,nquadsub,
     3      wnear,wnearsub,soln,novers,npts_over,ixyzso,srcover,
     4      wover,wtmp)

c          call lpcomp_helm_comb_dir_addsub(npatches,norders,ixyzs,
c     1      iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
c     2      eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,
c     3      soln,novers,npts_over,ixyzso,srcover,wover,wtmp)

C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rres)            
          do i=1,npts
            rres = rres + abs(zid*soln(i) + wtmp(i)-rhs(i))**2
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
c        
      subroutine helm_comb_dir_solver_memest(npatches,norders,ixyzs,
     1    iptype,npts,srccoefs,srcvals,eps,zpars,numit,
     2    rmem)
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
c     until a relative residual of eps_gmres is reached
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
c           numit - integer
c              max number of gmres iterations
c
c         output
c           rmem - real *8
c              estimated memory required by code in GB. Note that
c              this is meant to serve as an estimate only. 
c              The exact memory usage might be between (0.75,1.25)*rmem
c              The memory estimate may not be reliable for a
c              very small number of points
c 
c
c
      implicit none
      integer npatches,norder,npols,npts
      integer ifinout
      integer norders(npatches),ixyzs(npatches+1)
      integer iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps
      complex *16 zpars(3)
      real *8 rmem

      integer *8 lmem8,bigint
      real *8 rmemfmm
      real *8, allocatable :: targs(:,:)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      integer ndtarg,ntarg

      real *8 rres,eps2
      integer niter


      integer nover,npolso,nptso
      integer nnz,nquad
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      complex *16, allocatable :: wnear(:)

      real *8, allocatable :: srcover(:,:),wover(:),sources(:,:)
      integer, allocatable :: ixyzso(:),novers(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:) 

      integer i,j,jpatch,jquadstart,jstart

      integer ipars,ifcharge,ifdipole,ifpgh,ifpghtarg
      real *8 dpars,timeinfo(10),t1,t2,omp_get_wtime


      real *8 ttot,done,pi
      real *8 rfac,rfac0
      integer iptype_avg,norder_avg
      integer ikerorder, iquadtype,npts_over,iper

c
c
c       gmres variables
c
      integer numit,k,l
      complex *16 temp
      
      bigint = numit+1
      bigint = bigint*npts*2
      lmem8 = lmem8 + bigint



      lmem8 = lmem8 + numit*(numit+5)*2 + npts*2


      done = 1
      pi = atan(done)*4


c
c
c        setup targets as on surface discretization points
c 
      ndtarg = 3
      ntarg = npts
      allocate(targs(ndtarg,npts),uvs_targ(2,ntarg),ipatch_id(ntarg))
      lmem8 = lmem8 + ndtarg*npts + 3*ntarg 

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
      lmem8 = lmem8 + 5*npatches

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
      lmem8 = lmem8 + npts+1+ 2*nnz
      call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind,
     1         iquad)

      ikerorder = -1
      if(abs(zpars(3)).gt.1.0d-16) ikerorder = 0


c
c    estimate oversampling for far-field, and oversample geometry
c

      allocate(novers(npatches),ixyzso(npatches+1))
      lmem8 = lmem8 + 2*npatches + 1

      print *, "beginning far order estimation"

      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,
     1    rads,npts,srccoefs,ndtarg,npts,targs,ikerorder,zpars(1),
     2    nnz,row_ptr,col_ind,rfac,novers,ixyzso)

      npts_over = ixyzso(npatches+1)-1

      allocate(srcover(12,npts_over),wover(npts_over))
      lmem8 = lmem8 + 15*npts_over

      call oversample_geom(npatches,norders,ixyzs,iptype,npts, 
     1   srccoefs,srcvals,novers,ixyzso,npts_over,srcover)

      call get_qwts(npatches,novers,ixyzso,iptype,npts_over,
     1        srcover,wover)


c
c   compute near quadrature correction
c
      nquad = iquad(nnz+1)-1
      allocate(wnear(nquad))
      lmem8 = lmem8 + nquad*2
      rmem = lmem8*8/1024/1024/1024

      ifcharge = 0
      ifdipole = 0
      if(abs(zpars(2)).gt.1.0d-16) ifcharge = 1
      if(abs(zpars(3)).gt.1.0d-16) ifdipole = 1

      ifpgh = 0
      ifpghtarg = 1
      allocate(sources(3,npts_over))
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,npts_over
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)
      enddo
C$OMP END PARALLEL DO      

      iper = 0
      rmemfmm = 0
      call hfmm3d_memest(1,eps,zpars(1),npts_over,sources,ifcharge,
     1   ifdipole,iper,ifpgh,npts,targs,ifpghtarg,rmemfmm)
      rmem = rmem + rmemfmm
      
c
      return
      end
c
c
c
c
c
c
      subroutine helm_comb_dir_fds_csc_mem(npatches,norders,ixyzs,
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
      call findnearmem(cms,npatches,rad_near,ndtarg,targs,npts,nnz)

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,ndtarg,targs,npts,row_ptr, 
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
c
c  note oversampling turned off temporarily
c

      do i=1,npatches
        novers(i) = norders(i)
        ixyzso(i) = ixyzs(i)
      enddo
      ixyzso(npatches+1) = ixyzs(npatches+1)

cc      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,
cc     1    rads,npts,srccoefs,ndtarg,npts,targs,ikerorder,zpars(1),
cc     2    nnz,row_ptr,col_ind,rfac,novers,ixyzso)
      

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

      subroutine helm_comb_dir_fds_csc_init(npatches,norders,ixyzs,
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
      if(abs(zpars(3)).gt.1.0d-16) ikerorder = 0


c
c    estimate oversampling for far-field, and oversample geometry
c

c
c  note oversampling turned off currently
c   
       do i=1,npatches
         ifds(inovers+i-1) = norders(i)
         ifds(iixyzso+i-1) = ixyzs(i)
       enddo
       ifds(iixyzso+npatches) = ixyzs(npatches+1)

c      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,
c     1    rads,npts,srccoefs,ndtarg,npts,targs,ikerorder,zpars(1),
c     2    nnz,ifds(irow_ptr),ifds(icol_ind),rfac,ifds(inovers),
c     3    ifds(iixyzso))

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

      subroutine helm_comb_dir_fds_csc_matgen(npatches,norders,ixyzs,
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
      real *8 over4pi

      procedure (), pointer :: fker
      external h3d_slp, h3d_dlp, h3d_comb
      data over4pi/0.07957747154594767d0/

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

c
c    if no relevant targets for this source do nothing
c
        if(nn.le.0) goto 1111

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
        call zrmatmatt(naintbc,npolso,wquadf,npols,rfds(ixist),
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
 1111   continue        
      enddo
      
      
      


      return
      end
c
c
c
c
c
      subroutine helm_comb_dir_fds_block_mem(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,eps,zpars,ifwrite,
     2     nifds,nrfds,nzfds)
c
cf2py   intent(in) npatches,norders,ixyzs,iptype,npts
cf2py   intent(in) srccoefs,srcvals,eps,zpars,ifwrite
cf2py   intent(out) nifds,nrfds,nzfds

c
c       This subroutine estimates the memory requirements
c       for the precomputation routine of the fast direct solver
c       without oversampling
c 
c       The precomputation routine computes an integer array (ifds)
c       a real array (rfds), and complex array (zfds)
c 
c       The following quantities will be computed during the 
c       precomputation phase of the fast direct solver
c
c          ifds(1) - nnz
c          ifds(2) - nquad
c          ifds(3:4+npatches-1) - col_ptr (for near quadrature info)
c          ifds(4+npatches:4+npatches+nnz-1) - row_ind
c          ifds(4+npatches+nnz:4+npatches+2*nnz-1) - iper 
c          ifds(4+npatches+2*nnz:4+npatches+3*nnz) - iquad
c          zfds(1:3) - zpars(1:3)
c          zfds(4:4+nquad-1) - wnear
c
c        Thus this subroutine on output returns 
c          nifds = 4+npatches+3*nnz
c          nrfds = 0
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
      integer ifwrite

      real *8, allocatable :: targs(:,:)
      integer iptype_avg,norder_avg
      integer ntarg,ndtarg
      real *8 rfac,rfac0
      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)
      integer, allocatable :: iquad(:),row_ptr(:),col_ind(:)
      integer npts_over,nquad

      character *100 fname

      integer i


      if(ifwrite.eq.1) then
        write(fname,'(a,i5.5,a)') 'nrowcol-wtorus-',npatches,'.dat'
        open(unit=33,file=trim(fname),status='replace')
        write(33,*) "beginning nrow"
        close(33)
      
        write(fname,'(a,i5.5,a)') 'irowcol-wtorus-',npatches,'.dat'
        open(unit=33,file=trim(fname),status='replace')
        write(33,*) "Beginning irow"
        close(33)
      endif


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

      nquad = iquad(nnz+1)-1

      nifds = 4+npatches+3*nnz
      nrfds = 0 
      nzfds = 3 + nquad
      

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

      subroutine helm_comb_dir_fds_block_init(npatches,norders,
     1     ixyzs,iptype,npts,srccoefs,srcvals,eps,zpars,nifds,ifds,
     2     nzfds,zfds)
cf2py   intent(in) npatches,norders,ixyzs,iptype,npts
cf2py   intent(in) srccoefs,srcvals,eps,zpars
cf2py   intent(in) nifds,nzfds
cf2py   intent(out) ifds,zfds

c
c       This subroutine is the precomputation routine of the fast direct solver
c 
c       The following quantities will be computed during the 
c       precomputation phase of the fast direct solver
c
c          ifds(1) - nnz
c          ifds(2) - nquad
c          ifds(3:4+npatches-1) - col_ptr (for near quadrature info)
c          ifds(4+npatches:4+npatches+nnz-1) - row_ind
c          ifds(4+npatches+nnz:4+npatches+2*nnz-1) - iper 
c          ifds(4+npatches+2*nnz:4+npatches+3*nnz) - iquad
c
c          zfds(1:3) - zpars(1:3)
c          zfds(4:4+nquad-1) - wnear
c
c     
      implicit none
      integer npatches,norders(npatches),ixyzs(npatches+1)
      integer iptype(npatches),npts
      real *8 srccoefs(9,npts),srcvals(12,npts)
      real *8 eps
      complex *16 zpars(3)
      integer nnz
      integer nifds,nzfds
      integer ifds(nifds)
      complex *16 zfds(nzfds)

      real *8, allocatable :: targs(:,:)
      real *8, allocatable :: uvs_targ(:,:)
      integer, allocatable :: ipatch_id(:)
      integer, allocatable :: row_ptr(:),col_ind(:)
      integer iptype_avg,norder_avg
      integer ntarg,ndtarg
      real *8 rfac,rfac0
      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)
      integer npts_over,nquad

      integer iquadtype,istart,iend,i

      integer icol_ptr,irow_ind,iiquad,iiper




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

      icol_ptr = 3
      irow_ind = icol_ptr+npatches+1
      iiper = irow_ind+nnz
      iiquad = iiper+nnz

      allocate(row_ptr(npts+1),col_ind(nnz))

      
      call findnear(cms,npatches,rad_near,ndtarg,targs,npts,row_ptr,
     1   col_ind) 

      call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind,
     1   ifds(iiquad))

      call rsc_to_csc(npatches,npts,nnz,row_ptr,col_ind,
     1  ifds(icol_ptr),ifds(irow_ind),ifds(iiper))

      nquad = ifds(iiquad+nnz)-1

      ifds(1) = nnz
      ifds(2) = nquad

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
     1      ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,
     1      col_ind,ifds(iiquad),rfac0,nquad,zfds(4))
      
      return
      end
c
c
c
c
c
c
c

      subroutine helm_comb_dir_fds_block_matgen(npatches,norders,
     1     ixyzs,iptype,npts,srccoefs,srcvals,wts,eps,zpars,nifds,ifds,
     2     nzfds,zfds,nrow,row_ind,ncol,col_ind,ifwrite,zmat)
cf2py   intent(in) npatches,norders,ixyzs,iptype,npts
cf2py   intent(in) srccoefs,srcvals,eps,zpars,wts
cf2py   intent(in) nifds,nzfds
cf2py   intent(in) ifds,zfds,ifwrite
cf2py   intent(in) nrow,row_ind,ncol,col_ind 
cf2py   intent(out) zmat
      
      implicit real *8 (a-h,o-z)
      integer npatches,norders(npatches),ixyzs(npatches+1)
      integer iptype(npatches),npts
      real *8 srccoefs(9,npts),srcvals(12,npts)
      real *8 eps
      real *8 wts(npts)
      complex *16 zpars(3)
      integer nifds,nzfds
      integer ifds(nifds)
      complex *16 zfds(nzfds)
      integer nrow,ncol,row_ind(nrow),col_ind(ncol)
      integer ifwrite
      complex *16 zmat(nrow,ncol)

c
c        temporary variables
c
      integer ifcharge,ifdipole
      real *8, allocatable :: src(:,:),targ(:,:)
      complex *16, allocatable :: zcharge(:),zdipvec(:,:)
      integer, allocatable :: patch_ind(:),node_ind(:)
      integer, allocatable :: iuni_patch(:),iuni_patch_ind(:)
      integer, allocatable :: aintb(:,:),iaintba(:,:),iaintbb(:,:)
      integer, allocatable :: aintbc(:,:),iaintbc(:,:)
      integer, allocatable :: naintb(:),naintbc(:)

      integer icol_ptr,iiper,irow_ind
      real *8 over4pi

      procedure (), pointer :: fker
      complex *16 alpha,beta
      external h3d_slp, h3d_dlp, h3d_comb
      character *100 fname
      data over4pi/0.07957747154594767d0/


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
c
c       extract relevant source and target locations
c
c
      allocate(src(3,ncol),targ(3,nrow))
      allocate(zcharge(ncol),zdipvec(3,ncol))

      if(ifwrite.eq.1) then
        write(fname,'(a,i5.5,a)') 'nrowcol-wtorus-',npatches,'.dat'
        open(unit=33,file=trim(fname),access='append')
        write(33,*) nrow,ncol
        close(33)
      
        write(fname,'(a,i5.5,a)') 'irowcol-wtorus-',npatches,'.dat'
        open(unit=33,file=trim(fname),access='append')
        write(33,*) row_ind(1:nrow)
        write(33,*) col_ind(1:ncol)
        close(33)
      endif

C$OMP PARALLEL DO
      do i=1,ncol
        src(1,i) = srcvals(1,col_ind(i))
        src(2,i) = srcvals(2,col_ind(i))
        src(3,i) = srcvals(3,col_ind(i))
      enddo
C$OMP END PARALLEL DO      


      if(ifcharge.eq.1) then
C$OMP PARALLEL DO      
        do i=1,ncol
          zcharge(i) = wts(col_ind(i))*alpha*over4pi
        enddo
C$OMP END PARALLEL DO        
      endif

      if(ifdipole.eq.1) then
C$OMP PARALLEL DO      
        do i=1,ncol
          zdipvec(1,i) = srcvals(10,col_ind(i))*wts(col_ind(i))*beta*
     1       over4pi
          zdipvec(2,i) = srcvals(11,col_ind(i))*wts(col_ind(i))*beta*
     1       over4pi
          zdipvec(3,i) = srcvals(12,col_ind(i))*wts(col_ind(i))*beta*
     1       over4pi
        enddo
C$OMP END PARALLEL DO        
      endif

C$OMP PARALLEL DO
      do i=1,nrow
        targ(1,i) = srcvals(1,row_ind(i))
        targ(2,i) = srcvals(2,row_ind(i))
        targ(3,i) = srcvals(3,row_ind(i))
      enddo
C$OMP END PARALLEL DO      


c
c    note that thresh does not matter since near entries will
c    be corrected anyway
c
      thresh = 0
      nd = 1
c
c
c       fill out uncorrected matrix
c

C$OMP PARALLEL DO PRIVATE(i,j)
      do i=1,ncol
        do j=1,nrow
          zmat(j,i) = 0
        enddo

        if(ifcharge.eq.1.and.ifdipole.eq.0) then
          call h3ddirectcp(nd,zpars(1),src(1,i),zcharge(i),1,targ,
     1       nrow,zmat(1,i),thresh)
        endif
        if(ifcharge.eq.0.and.ifdipole.eq.1) then
          call h3ddirectdp(nd,zpars(1),src(1,i),zdipvec(1,i),1,targ,
     1       nrow,zmat(1,i),thresh)
        endif
        if(ifcharge.eq.1.and.ifdipole.eq.1) then
          call h3ddirectcdp(nd,zpars(1),src(1,i),zcharge(i),
     1      zdipvec(1,i),1,targ,nrow,zmat(1,i),thresh)
        endif
      enddo
C$OMP END PARALLEL DO      

      nnz = ifds(1)
      nquad = ifds(2)

      allocate(patch_ind(ncol),node_ind(ncol))

      call col_ind_to_patch_node_ind(npatches,ixyzs,ncol,col_ind,
     1       patch_ind,node_ind)

      icol_ptr = 3
      irow_ind = icol_ptr+npatches+1
      iiper = irow_ind+nnz
      iiquad = iiper+nnz

      allocate(iuni_patch(ncol),iuni_patch_ind(ncol))

      nuni = 0
      call get_iuni1(ncol,patch_ind,nuni,iuni_patch,iuni_patch_ind)


      allocate(aintb(nrow,nuni),iaintba(nrow,nuni),iaintbb(nrow,nuni))
      allocate(aintbc(nrow,nuni),iaintbc(nrow,nuni),naintb(nuni))
      allocate(naintbc(nuni))


C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ipatch,m,istart)
      do i=1,nuni
        ipatch = iuni_patch(i)
c
c     collect set b
c
        m = ifds(icol_ptr+ipatch)-ifds(icol_ptr+ipatch-1)
        istart = ifds(icol_ptr+ipatch-1)
        naintb(i) = 0
        naintbc(i) = 0
        call setdecomp(nrow,row_ind,m,ifds(irow_ind+istart-1),
     1    naintb(i),aintb(1,i),iaintba(1,i),iaintbb(1,i),
     2    naintbc(i),aintbc(1,i),iaintbc(1,i))
      enddo
C$OMP END PARALLEL DO      

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(icol,iuni,inode,ipatch)
C$OMP$PRIVATE(j,jind,jrow,jquad,jquadloc)
      do icol=1,ncol
        iuni = iuni_patch_ind(icol)
        inode = node_ind(icol)
        ipatch = patch_ind(icol)
        do j=1,naintb(iuni)
          jind = iaintba(j,iuni)
          jrow = iaintbb(j,iuni)+ ifds(icol_ptr+ipatch-1)-1
          jquadloc = ifds(iiper+jrow-1)
          jquad = ifds(iiquad+jquadloc-1)
          zmat(jind,icol) = zfds(3+jquad+inode-1)
        enddo
      enddo
C$OMP END PARALLEL DO      




      return
      end
