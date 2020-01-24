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
c       lpcomp_helm_comb_dir_addsub (to be tested) -
c         compute layer potential for the Dirichlet
c         data corresponding to the combined field representation
c         using add and subtract
c 
c
c       lpcomp_helm_comb_dir_setsub (to be tested) (OPTIONAL)
c          compute layer potential for Dirichlet data corresponding
c          to the combined field representation using 
c          set subtraction and turning off list 1
c
c     TO WRITE:
c       setparams0_fmm_helm_comb_dir
c
c       setparams_fmm_helm_comb_dir
c
c       helm_comb_dir_mv
c
c       helm_comb_dir_solver
c
c       setparams0_fds_helm_comb_dir
c 
c       setparams_fds_helm_comb_dir
c
c       helm_comb_dir_matgen
c
c


      subroutine getnearquad_helm_comb_dir(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2   ipatch_id,uvs_targ,eps,zpars,nnz,row_ptr,col_ind,
     3   iquad,rfac0,nquad,wnear)
c
c       this subroutine generates the near field quadrature
c       for the representation u = \alpha S_{k}  + \beta D_{k} ---(1)
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
c          (maybe better to find closest uv on patch using
c            newton)
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
c           nnz - integer
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

      implicit none 
      integer npatches,norders(npatches),npts,nquad
      integer ixyzs(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,rfac0
      integer ndtarg,ntarg
      real *8 targs(ndtarg,ntarg)
      complex *16 zpars(3)
      integer nnz,ndtarg,ipars
      real *8 dpars
      integer row_ptr(npts+1),col_ind(nnz),iquad(nnz+1)
      complex *16 wnear(nquad)

      integer ipatch_id(ntarg)
      real *8 uvs_targ(2,ntarg)

      complex *16 alpha,beta
      integer i,j

      procedure (), pointer :: fker
      external h3d_ggq_slp, h3d_ggq_dlp, h3d_ggq_comb

c
c
c        initialize the appropriate kernel function
c
      alpha = zpars(2)
      beta = zpars(3)
      fker => h3d_ggq_comb
      if(abs(alpha).ge.1.0d-16.and.abs(beta).lt.1.0d-16) then
        fker=>h3d_ggq_slp
      else if(abs(alpha).lt.1.0d-16.and.abs(beta).ge.1.0d-16) then
        fker=>h3d_ggq_dlp
      endif

      call getnearquad_compact_guru(npatches,norders,ixyzs,
     1   iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     1   ipatch_id,uvs_targ,
     1   eps,fker,dpars,zpars,ipars,nnz,row_ptr,col_ind,iquad,
     1   rfac0,nquad,wnear)

      return
      end
c
c
c
c
c
c
c
      subroutine lpcomp_helm_comb_dir_addsub(npatches,norders,ixyzs,
     1   iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,inode_id,
     2   eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,sigma,novers,
     3   nptso,ixyzso,srcover,whtsover,ifinout,pot)
c
c
c      this subroutine evaluates the layer potential for
c      the representation u = \alpha S_{k} + \beta D_{k} 
c      where the near field is precomputed and stored
c      in the row sparse compressed format.
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
c         ipt_id - integer(ntarg)
c            id of node of target i, id = -1, if target is off-surface
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
c           ifinout - integer
c              flag for computing interior/exterior limit of layer
c              potential at boundary target points
c
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
      integer inode_id(ntarg)
      complex *16 zpars(3)
      integer nnz,row_ptr(npts+1),col_ind(nnz),nquad
      integer iquad(nnz+1)
      complex *16 wnear(nquad),sigma(npts)
      integer novers(npatches+1)
      integer nover,npolso,nptso
      real *8 srcover(12,nptso),whtsover(nptso)
      complex *16 pot(npts)
      complex *16, allocatable :: potsort(:)

      real *8, allocatable :: sources(:,:),targvals(:,:)
      complex *16, allocatable :: charges(:),dipvec(:,:),sigmaover(:)
      integer ns,nt
      complex *16 alpha,beta
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      complex *16 tmp(10),val


      integer i,j,jpatch,jquadstart,jstart


      integer ifaddsub

      integer ntarg,ntj
      
      complex *16 zdotu,pottmp
      complex *16, allocatable :: ctmp2(:),dtmp2(:,:)
      real *8 radexp,epsfmm

      integer ipars
      real *8 dpars,timeinfo(10),t1,t2,omp_get_wtime

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
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
c
c       add in precomputed quadrature
c

      thresh = 1.0d-16
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
c
c
c         check sign of identity term
c
        if(inode_id(i).gt.0) pot(i) = pot(i) - 
     1      (-1)**(ifinout)*sigma(inode_id(i))/2/pi*zpars(3)
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


      call prin2('quadrature time=*',timeinfo,2)
      
      ttot = timeinfo(1) + timeinfo(2)
      call prin2('time in lpcomp=*',ttot,1)

      
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
     1   iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,inode_id,
     2   eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,sigma,novers,
     3   nptso,ixyzso,srcover,whtsover,ifinout,pot)
c
c
c      this subroutine evaluates the layer potential for
c      the representation u = \alpha S_{k} + \beta D_{k} 
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
c         ipt_id - integer(ntarg)
c            id of node of target i, id = -1, if target is off-surface
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
c           ifinout - integer
c              flag for computing interior/exterior limit of layer
c              potential at boundary target points
c
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
      integer inode_id(ntarg)
      complex *16 zpars(3)
      integer nnz,row_ptr(npts+1),col_ind(nnz),nquad
      integer iquad(nnz+1)
      complex *16 wnear(nquad),sigma(npts)
      integer novers(npatches+1)
      integer nover,npolso,nptso
      real *8 srcover(12,nptso),whtsover(nptso)
      complex *16 pot(npts)
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
      integer nlist1,nlmax,npover,nl1,ntarg,ntj
      
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

      thresh = 1.0d-16
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
c
c
c         check sign of identity term
c
        if(inode_id(i).gt.0) pot(i) = pot(i) - 
     1      (-1)**(ifinout)*sigma(inode_id(i))/2/pi*zpars(3)
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

      call prin2('quadrature time=*',timeinfo,2)
      
      ttot = timeinfo(1) + timeinfo(2)
      call prin2('time in lpcomp=*',ttot,1)

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
