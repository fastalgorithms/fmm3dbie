c
c     Helmholtz Beltrami solvers, based on the representations
c     This file contains the following user callable
c     routines: 
c 
c       getnearquad_helm_bel_<rep> - 
c        generates the near  field quadrature required for
c        representation
c
c
c
c     For the quadrature generation routines, rep could also be
c     log: implements u = \int_{\Gamma} H_0(zk|x-y|) \sigma(y) dy
c
c     The integral reprensentation for this case then becomes
c        \int_{\Gamma} -zk^2 H_2(_) (n(x). (x-y)/|x-y|)^2 
c             + zk 2H(x) H_1(_) n(x).(x-y)/|x-y|
c
c     where H(x) is the mean curvature
c
c  
c     Note that the addsub routines might have different input and
c      output structures, since one of them requires the additional
c      input of S[1]
c
c     In the solvers, we implement a low stagnation threshold stable
c     version of GMRES which benefits the two integral equations
c     where the identitiy term is explicitly pulled out, and
c     we also use a weighted Krylov space to account for the 
c     the piecewise spectral representation of the density
c     on the surface
c

c
c
      subroutine getnearquad_helm_bel_hank(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,eps,
     2   iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,
     1   meancrvs,nquad,zpars,iktype,wnear)
c
c
c  This subroutine generates the near field quadrature
c  for the integral equation:
c
c        \int_{\Gamma} (n(x). (x-y)/|x-y|)^2 - 2H(x) n(x).(x-y)/|x-y|^2
c
c  where H(x) is the mean curvature
c  
c  The quadrature is computed by the following strategy
c  targets within a sphere of radius rfac0*rs
c  of a patch centroid is handled using adaptive integration
c  where rs is the radius of the bounding sphere
c  for the patch
c  
c  All other targets in the near field are handled via
c  oversampled quadrature
c
c  The recommended parameter for rfac0 is 1.25d0
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
c        iptype = 1, triangular patch discretized using RV nodes
c    - npts: integer
c        total number of discretization points on the boundary
c    - srccoefs: real *8 (9,npts)
c        koornwinder expansion coefficients of xyz, dxyz/du,
c        and dxyz/dv on each patch. 
c        For each point 
c          * srccoefs(1:3,i) is xyz info
c          * srccoefs(4:6,i) is dxyz/du info
c          * srccoefs(7:9,i) is dxyz/dv info
c    - srcvals: real *8 (12,npts)
c        xyz(u,v) and derivative info sampled at the 
c        discretization nodes on the surface
c          * srcvals(1:3,i) - xyz info
c          * srcvals(4:6,i) - dxyz/du info
c          * srcvals(7:9,i) - dxyz/dv info
c          * srcvals(10:12,i) - normals info
c    - eps: real *8
c        precision requested
c    - iquadtype: integer
c        quadrature type
c          * iquadtype = 1, use ggq for self + adaptive integration
c            for rest
c    - nnz: integer
c        number of source patch-> target interactions in the near field
c    - row_ptr: integer(npts+1)
c        row_ptr(i) is the pointer
c        to col_ind array where list of relevant source patches
c        for target i start
c    - col_ind: integer (nnz)
c        list of source patches relevant for all targets, sorted
c        by the target number
c    - iquad: integer(nnz+1)
c        location in wnear_ij array where quadrature for col_ind(i)
c        starts for a single kernel. In this case the different kernels
c        are matrix entries are located at (m-1)*nquad+iquad(i), where
c        m is the kernel number
c    - rfac0: integer
c        radius parameter for near field
c    - nquad: integer
c        number of near field entries corresponding to each source target
c        pair. The size of wnear is (nquad,3) since there are 4 kernels
c        per source target pair
c
c  Output arguments
c    - wnear: real *8(nquad)
c        The desired near field quadrature
c               
c

      implicit none 
      integer *8, intent(in) :: npatches,norders(npatches),npts,nquad
      integer *8, intent(in) :: ixyzs(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: iquadtype
      integer *8, intent(in) :: nnz
      integer *8, intent(in) :: row_ptr(npts+1),col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      complex *16, intent(out) :: wnear(nquad)
      real *8, intent(in) :: meancrvs(npts)


      integer *8, intent(in) :: iktype
      complex *16, intent(in) :: zpars

      integer *8, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_src(:,:), targvals(:,:)
      integer *8 ipars(1)
      integer *8 ndd,ndz,ndi
      real *8 dpars

      integer *8 ndtarg

      real *8 done,pi
      integer *8 i,j
      integer *8 ipv

      procedure (), pointer :: fker
      external helm_bel_hank,helm_bel_res

      done = 1
      pi = atan(done)*4

      allocate(ipatch_id(npts),uvs_src(2,npts))
      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,
     1   ipatch_id,uvs_src)

c
c  get mean curvature 
c
ccc      allocate(meancrvs(npts))
ccc      call get_mean_curvature(npatches, norders, ixyzs, iptype, 
ccc     1  npts, srccoefs, srcvals, meancrvs)
      
      allocate(targvals(13,npts))
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
      do i=1,npts
        do j=1,12
          targvals(j,i) = srcvals(j,i)
        enddo
        targvals(13,i) = meancrvs(i)
      enddo
C$OMP END PARALLEL DO      

c
c
c        initialize the appropriate kernel function
c

      ndd = 0
      ndi = 0
      ndz = 1
      ndtarg = 13
      if(iquadtype.eq.1) then
        if (iktype .eq. 1) then
            ipv = 0
            fker=>helm_bel_res
        elseif (iktype .eq. 2) then
            ipv = 0
            fker=>helm_bel_hank
        endif
        call zgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,targvals,
     1     ipatch_id,uvs_src,
     1     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     1     col_ind,iquad,
     1     rfac0,nquad,wnear)
      endif


      return
      end
c
c
c
c
c

      subroutine getnearquad_lap_bel_log(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,eps,
     2   iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,
     1   meancrvs,nquad,zpars,iktype,wnear)
c
c
c  This subroutine generates the near field quadrature
c  for the integral equation:
c
c        \int_{\Gamma} (n(x). (x-y)/|x-y|)^2 - 2H(x) n(x).(x-y)/|x-y|^2
c
c  where H(x) is the mean curvature
c  
c  The quadrature is computed by the following strategy
c  targets within a sphere of radius rfac0*rs
c  of a patch centroid is handled using adaptive integration
c  where rs is the radius of the bounding sphere
c  for the patch
c  
c  All other targets in the near field are handled via
c  oversampled quadrature
c
c  The recommended parameter for rfac0 is 1.25d0
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
c        iptype = 1, triangular patch discretized using RV nodes
c    - npts: integer
c        total number of discretization points on the boundary
c    - srccoefs: real *8 (9,npts)
c        koornwinder expansion coefficients of xyz, dxyz/du,
c        and dxyz/dv on each patch. 
c        For each point 
c          * srccoefs(1:3,i) is xyz info
c          * srccoefs(4:6,i) is dxyz/du info
c          * srccoefs(7:9,i) is dxyz/dv info
c    - srcvals: real *8 (12,npts)
c        xyz(u,v) and derivative info sampled at the 
c        discretization nodes on the surface
c          * srcvals(1:3,i) - xyz info
c          * srcvals(4:6,i) - dxyz/du info
c          * srcvals(7:9,i) - dxyz/dv info
c          * srcvals(10:12,i) - normals info
c    - eps: real *8
c        precision requested
c    - iquadtype: integer
c        quadrature type
c          * iquadtype = 1, use ggq for self + adaptive integration
c            for rest
c    - nnz: integer
c        number of source patch-> target interactions in the near field
c    - row_ptr: integer(npts+1)
c        row_ptr(i) is the pointer
c        to col_ind array where list of relevant source patches
c        for target i start
c    - col_ind: integer (nnz)
c        list of source patches relevant for all targets, sorted
c        by the target number
c    - iquad: integer(nnz+1)
c        location in wnear_ij array where quadrature for col_ind(i)
c        starts for a single kernel. In this case the different kernels
c        are matrix entries are located at (m-1)*nquad+iquad(i), where
c        m is the kernel number
c    - rfac0: integer
c        radius parameter for near field
c    - nquad: integer
c        number of near field entries corresponding to each source target
c        pair. The size of wnear is (nquad,3) since there are 4 kernels
c        per source target pair
c
c  Output arguments
c    - wnear: real *8(nquad)
c        The desired near field quadrature
c               
c

      implicit none 
      integer *8, intent(in) :: npatches,norders(npatches),npts,nquad
      integer *8, intent(in) :: ixyzs(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      real *8, intent(in) :: rfac0
      integer *8, intent(in) :: iquadtype
      integer *8, intent(in) :: nnz
      integer *8, intent(in) :: row_ptr(npts+1),col_ind(nnz)
      integer *8, intent(in) :: iquad(nnz+1)
      real *8, intent(out) :: wnear(nquad)
      real *8, intent(in) :: meancrvs(npts)


      integer *8, intent(in) :: iktype
      complex *16, intent(in) :: zpars

      integer *8, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_src(:,:), targvals(:,:)
      integer *8 ipars(1)
      integer *8 ndd,ndz,ndi
      real *8 dpars

      integer *8 ndtarg

      real *8 done,pi
      integer *8 i,j
      integer *8 ipv

      procedure (), pointer :: fker
      external lap_bel_log,lap_bel_res

      done = 1
      pi = atan(done)*4

      allocate(ipatch_id(npts),uvs_src(2,npts))
      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,
     1   ipatch_id,uvs_src)

c
c  get mean curvature 
c
ccc      allocate(meancrvs(npts))
ccc      call get_mean_curvature(npatches, norders, ixyzs, iptype, 
ccc     1  npts, srccoefs, srcvals, meancrvs)
      
      allocate(targvals(13,npts))
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
      do i=1,npts
        do j=1,12
          targvals(j,i) = srcvals(j,i)
        enddo
        targvals(13,i) = meancrvs(i)
      enddo
C$OMP END PARALLEL DO      

c
c
c        initialize the appropriate kernel function
c

      ndd = 0
      ndi = 0
      ndz = 1
      ndtarg = 13
      if(iquadtype.eq.1) then
        if (iktype .eq. 1) then
            ipv = 0
            fker=>lap_bel_res
        elseif (iktype .eq. 2) then
            ipv = 0
            fker=>lap_bel_log
        endif
        call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,targvals,
     1     ipatch_id,uvs_src,
     1     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     1     col_ind,iquad,
     1     rfac0,nquad,wnear)
      endif


      return
      end
c
c
c
c
c

