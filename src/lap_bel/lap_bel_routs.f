c
c     Laplace Beltrami solvers, based on the representations
c     S (\Delta_{\Gamma}+ \int_{\Gamma} )S \sigma = f
c
c     with either using the form
c     (\nabla \cdot S [\nabla_{\Gamma} S]  + S W S) \sigma = f
c
c     or using Calderon identities
c     (-1/4 + D^2 - S(S'' + D' + 2HS')) [\sigma] = S[f]
c
c
c     and also using the representation
c     (\Delta_{\Gamma} S^2 + WS^2) \sigma = f 
c     
c     using Calderon identities
c     -1/4 + S'^2 - (S'' + D' + 2HS') + WS^2 = f
c     
c
c     This file contains the following user callable
c     routines: 
c 
c       getnearquad_lap_bel_<rep> - 
c        generates the near  field quadrature required for
c        representation
c
c       lpcomp_lap_bel_<rep>_addsub 
c          Apply integral representation using add and subtract 
c          
c
c       lap_bel_solver - Solves the laplace beltrami integral
c         equation, using one the three representations
c         above
c
c     For all of these routines, rep could be
c       s_lap_s: implements S(\Delta_{\Gamma} + W)S using Calderon
c              identities
c       s_lap_s_noc: implements S(\Delta_{\Gamma} + W) S without
c                    Calderon identities
c       lap_s2: implements (\Delta_{\Gamma} + S^2) \sigma
c                using Calderon identities
c
c
c     For the quadrature generation routines, rep could also be
c     log: implements u = \int_{\Gamma} \log|x-y| \sigma(y) dy
c
c     The integral reprensentation for this case then becomes
c        \int_{\Gamma} (n(x). (x-y)/|x-y|)^2 - 2H(x) n(x).(x-y)/|x-y|^2
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

      subroutine getnearquad_lap_bel_s_lap_s(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,eps,
     2   iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear)
c
c
c  This subroutine generates the near field quadrature
c  for the integral equation:
c
c			-sigma/4+D^2 - S(S''+ D' +2HS')=S[f]
c
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
c  NOTES:
c    - wnear must be of size (nquad,4)
c      * The first kernel is S
c      * The second kernel is D
c      * The third kernel is S'
c      * The fourth kernel is S'' + D'
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
c        pair. The size of wnear is (nquad,4) since there are 4 kernels
c        per source target pair
c
c  Output arguments
c    - wnear: real *8(nquad,4)
c        The desired near field quadrature
c        * The first kernel is S
c        * The second kernel is D
c        * The third kernel is S'
c        * The fourth kernel is S'' + D'
c               
c

      implicit none 
      integer, intent(in) :: npatches,norders(npatches),npts,nquad
      integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      real *8, intent(in) :: rfac0
      integer, intent(in) :: iquadtype
      integer, intent(in) :: nnz
      integer, intent(in) :: row_ptr(npts+1),col_ind(nnz),iquad(nnz+1)
      real *8, intent(out) :: wnear(nquad,4)


      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_src(:,:)
      integer ipars
      integer ndd,ndz,ndi
      complex *16 zpars
      real *8 dpars
      real *8 eps_spp_dp

      integer ndtarg

      real *8 alpha,beta,done,pi
      integer i,j
      integer ipv

      procedure (), pointer :: fker
      external l3d_slp, l3d_dlp, l3d_spp_sum_dp, l3d_sprime

      done = 1
      pi = atan(done)*4

      allocate(ipatch_id(npts),uvs_src(2,npts))
      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,
     1   ipatch_id,uvs_src)

c
c
c        initialize the appropriate kernel function
c


      ndd = 0
      ndi = 0
      ndz = 0
      ndtarg = 12
      if(iquadtype.eq.1) then
        ipv = 1
        fker=>l3d_slp
        call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,srcvals,
     1     ipatch_id,uvs_src,
     1     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     1     col_ind,iquad,
     1     rfac0,nquad,wnear(1,1))
        print *, "done with kernel 1"

        fker=>l3d_dlp
        call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,srcvals,
     1     ipatch_id,uvs_src,
     1     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     1     col_ind,iquad,
     1     rfac0,nquad,wnear(1,2))
        print *, "done with kernel 2"

        fker=>l3d_sprime
        call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,srcvals,
     1     ipatch_id,uvs_src,
     1     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     1     col_ind,iquad,
     1     rfac0,nquad,wnear(1,3))
        print *, "done with kernel 3"

        eps_spp_dp = 0.5d-8
        fker=>l3d_spp_sum_dp
        call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,srcvals,
     1     ipatch_id,uvs_src,
     1     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,
     1     row_ptr,col_ind,iquad,rfac0,nquad,wnear(1,4))
        print *, "done with kernel 4"

      endif


      return
      end
c
c
c
c
c
c
c
      subroutine getnearquad_lap_bel_s_lap_s_noc(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,eps,
     2   iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear)
c
c
c  This subroutine generates the near field quadrature
c  for the integral equation:
c
c  (\nabla \cdot S([\nabla_{\Gamma} S]) + SWS) sigma = S[f]
c
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
c  NOTES:
c    - wnear must be of size (nquad,4)
c      * The first kernel is S
c      * The second kernel is \nabla_{x} S
c      * The third kernel is \nabla_{y} S
c      * The fourth kernel is \nabla_{z} S
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
c        pair. The size of wnear is (nquad,4) since there are 4 kernels
c        per source target pair
c
c  Output arguments
c    - wnear: real *8(nquad,4)
c        The desired near field quadrature
c        * The first kernel is S
c        * The second kernel is \nabla_{x} S
c        * The third kernel is \nabla_{y} S
c        * The fourth kernel is \nabla_{z} S
c               
c

      implicit none 
      integer, intent(in) :: npatches,norders(npatches),npts,nquad
      integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      real *8, intent(in) :: rfac0
      integer, intent(in) :: iquadtype
      integer, intent(in) :: nnz
      integer, intent(in) :: row_ptr(npts+1),col_ind(nnz),iquad(nnz+1)
      real *8, intent(out) :: wnear(nquad,4)


      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_src(:,:)
      integer ipars
      integer ndd,ndz,ndi
      complex *16 zpars
      real *8 dpars

      integer ndtarg

      real *8 alpha,beta,done,pi
      integer i,j
      integer ipv

      procedure (), pointer :: fker
      external l3d_slp, l3d_sgradx, l3d_sgrady, l3d_sgradz

      done = 1
      pi = atan(done)*4

      allocate(ipatch_id(npts),uvs_src(2,npts))
      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,
     1   ipatch_id,uvs_src)

c
c
c        initialize the appropriate kernel function
c


      ndd = 0
      ndi = 0
      ndz = 0
      ndtarg = 12
      if(iquadtype.eq.1) then
        ipv = 1
        fker=>l3d_slp
        call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,srcvals,
     1     ipatch_id,uvs_src,
     1     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     1     col_ind,iquad,
     1     rfac0,nquad,wnear(1,1))
        print *, "done with kernel 1"

        fker=>l3d_sgradx
        call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,srcvals,
     1     ipatch_id,uvs_src,
     1     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     1     col_ind,iquad,
     1     rfac0,nquad,wnear(1,2))
        print *, "done with kernel 2"

        fker=>l3d_sgrady
        call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,srcvals,
     1     ipatch_id,uvs_src,
     1     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     1     col_ind,iquad,
     1     rfac0,nquad,wnear(1,3))
        print *, "done with kernel 3"

        fker=>l3d_sgradz
        call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,srcvals,
     1     ipatch_id,uvs_src,
     1     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     1     col_ind,iquad,
     1     rfac0,nquad,wnear(1,4))
        print *, "done with kernel 4"

      endif


      return
      end
c
c
c
c
c
c
c
      subroutine getnearquad_lap_bel_lap_s2(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,eps,
     2   iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear)
c
c
c  This subroutine generates the near field quadrature
c  for the integral equation:
c
c			-sigma/4+S'^2 - (S''+ D' +2HS')S = f
c
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
c  NOTES:
c    - wnear must be of size (nquad,3)
c      * The first kernel is S
c      * The second kernel is S'
c      * The third kernel is S'' + D' 
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
c    - wnear: real *8(nquad,3)
c        The desired near field quadrature
c        * The first kernel is S
c        * The second kernel is S'
c        * The third kernel is S'' + D'
c               
c

      implicit none 
      integer, intent(in) :: npatches,norders(npatches),npts,nquad
      integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      real *8, intent(in) :: rfac0
      integer, intent(in) :: iquadtype
      integer, intent(in) :: nnz
      integer, intent(in) :: row_ptr(npts+1),col_ind(nnz),iquad(nnz+1)
      real *8, intent(out) :: wnear(nquad,3)


      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_src(:,:)
      integer ipars
      integer ndd,ndz,ndi
      complex *16 zpars
      real *8 dpars

      integer ndtarg

      real *8 alpha,beta,done,pi
      integer i,j
      integer ipv

      procedure (), pointer :: fker
      external l3d_slp, l3d_sprime, l3d_spp_sum_dp

      done = 1
      pi = atan(done)*4

      allocate(ipatch_id(npts),uvs_src(2,npts))
      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,
     1   ipatch_id,uvs_src)

c
c
c        initialize the appropriate kernel function
c


      ndd = 0
      ndi = 0
      ndz = 0
      ndtarg = 12
      if(iquadtype.eq.1) then
        ipv = 1
        fker=>l3d_slp
        call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,srcvals,
     1     ipatch_id,uvs_src,
     1     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     1     col_ind,iquad,
     1     rfac0,nquad,wnear(1,1))
        print *, "done with kernel 1"

        fker=>l3d_sprime
        call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,srcvals,
     1     ipatch_id,uvs_src,
     1     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     1     col_ind,iquad,
     1     rfac0,nquad,wnear(1,2))
        print *, "done with kernel 2"

        fker=>l3d_spp_sum_dp
        call dgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,srcvals,
     1     ipatch_id,uvs_src,
     1     eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     1     col_ind,iquad,
     1     rfac0,nquad,wnear(1,3))
        print *, "done with kernel 3"

      endif


      return
      end
c
c
c
c
c
c
      subroutine getnearquad_lap_bel_log(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,eps,
     2   iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear)
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
      integer, intent(in) :: npatches,norders(npatches),npts,nquad
      integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      real *8, intent(in) :: rfac0
      integer, intent(in) :: iquadtype
      integer, intent(in) :: nnz
      integer, intent(in) :: row_ptr(npts+1),col_ind(nnz),iquad(nnz+1)
      real *8, intent(out) :: wnear(nquad)


      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_src(:,:), meancrvs(:), targvals(:,:)
      integer ipars
      integer ndd,ndz,ndi
      complex *16 zpars
      real *8 dpars

      integer ndtarg

      real *8 alpha,beta,done,pi
      integer i,j
      integer ipv

      procedure (), pointer :: fker
      external lap_bel_log

      done = 1
      pi = atan(done)*4

      allocate(ipatch_id(npts),uvs_src(2,npts))
      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,
     1   ipatch_id,uvs_src)

c
c  get mean curvature 
c
      allocate(meancrvs(npts))
      call get_mean_curvature(npatches, norders, ixyzs, iptype, 
     1  npts, srccoefs, srcvals, meancrvs)
      
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
      ndz = 0
      ndtarg = 13
      if(iquadtype.eq.1) then
        ipv = 1
        fker=>lap_bel_log
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
c
      subroutine lpcomp_lap_bel_s_lap_s_addsub(npatches,norders,ixyzs,
     1   iptype,npts,srccoefs,srcvals,eps,nnz,row_ptr,col_ind,iquad,
     2   nquad,wnear,sigma,novers,nptso,ixyzso,srcover,whtsover,pot,
     3   s_pot)
c
c
c
c  This subroutine evaluates the layer potential for the integral
c  representation:
c
c			-sigma/4+D^2 - S(S''+ D' +2HS')\sigma + SWS \sigma
c
c  where the near field is precomputed and stored in the row
c  sparse compressed format.
c
c  The subroutine also returns S[\sigma] which is the representation
c  for the solution
c
c  The fmm is used to accelerate the far-field with two calls to the
c  vector fmm, and the near-field interactions are 
c  handled via precomputed quadrature.
c
c  Using add and subtract: no need to call tree and set fmm parameters,
c  can call existing fmm library directly
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
c    - nquad: integer
c        number of near field entries corresponding to each source target
c        pair. The size of wnear is (nquad,4) since there are 4 kernels
c        per source target pair
c    - wnear: real *8(nquad,4)
c        The desired near field quadrature
c        * The first kernel is S
c        * The second kernel is D
c        * The third kernel is S'
c        * The fourth kernel is S'' + D'
c    - sigma: real *8(npts)
c        density for the layer potential
c    - novers: integer(npatches)
c        order of discretization for oversampled sources and
c        density
c    - ixyzso: integer(npatches+1)
c        ixyzso(i) denotes the starting location in srcover,
c        corresponding to patch i
c    - nptso: integer
c        total number of oversampled points
c    - srcover: real *8 (12,nptso)
c        oversampled set of source information
c    - whtsover: real *8 (nptso)
c        smooth quadrature weights at oversampled nodes
c
c  Output arguments
c    - pot: real *8 (npts)
c        (S \Delta_{\Gamma} S + SWS) \sigma applied using Calderon
c        identities
c------------------------

      implicit none
      integer, intent(in) :: npatches,npts
      integer, intent(in) :: norders(npatches),ixyzs(npatches+1)
      integer, intent(in) :: ixyzso(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      integer, intent(in) :: nnz,row_ptr(npts+1),col_ind(nnz),nquad
      integer, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: wnear(nquad,4),sigma(npts)
      integer, intent(in) :: novers(npatches+1)
      integer, intent(in) :: nptso
      real *8, intent(in) :: srcover(12,nptso),whtsover(nptso)
      real *8, intent(out) :: pot(npts),s_pot(npts)

      integer norder,npols,nover,npolso

      integer ndtarg,ntarg

      real *8, allocatable :: sources(:,:),targvals(:,:)
      real *8, allocatable :: charges(:),dipvec(:,:),sigmaover(:)

      real *8, allocatable :: pot2(:)
      real *8, allocatable :: pot3(:)
      real *8, allocatable :: curv(:)
      real *8, allocatable :: abc(:,:),pot_tmp(:,:),grad_tmp(:,:,:)
      real *8, allocatable :: hess_tmp(:,:,:),charges_tmp(:,:)
      real *8, allocatable :: dipvec_tmp(:,:,:),abcover(:,:)
      integer ns,nt
      real *8 alpha,beta
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      real *8 tmp(10),val(2),vgrad(2,3),nvgrad,nvhess
      real *8 vhess(2,6)
      real *8 xmin,xmax,ymin,ymax,zmin,zmax,sizey,sizez,boxsize
      integer i,j,jpatch,jquadstart,jstart

      integer ifaddsub

      integer ntj
      
      real *8 ddot,pottmp
      real *8, allocatable :: ctmp2(:,:),dtmp2(:,:,:)
      real *8, allocatable :: wts(:)
      real *8 radexp,epsfmm

      integer ipars
      complex *16 zpars
      real *8 timeinfo(10),t1,t2,omp_get_wtime
      real *8 dpars(2)

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      real *8 thresh,ra
      real *8 rr,rmin,rint,rsurf
      integer nss,ii,l,npover

      integer nd,ntarg0,nmax
      integer ier,iper

      real *8 ttot,done,pi,over4pi

      parameter (ntarg0=1)



      ns = nptso
      done = 1.0d0
      pi = atan(done)*4.0d0
      over4pi = 1.0d0/4.0d0/pi

      nmax = 0
      call get_near_corr_max(npts,row_ptr,nnz,col_ind,npatches,
     1   ixyzso,nmax)
      allocate(srctmp2(3,nmax),ctmp2(2,nmax),dtmp2(2,3,nmax))

      nd = 2
      allocate(abc(2,npts),charges_tmp(2,ns),dipvec_tmp(2,3,ns))
      allocate(pot_tmp(2,npts),grad_tmp(2,3,npts),hess_tmp(2,6,npts))
      allocate(abcover(2,ns),sigmaover(ns),curv(npts))

c
c   Set source info and targinfo for the fmm
c
c
      allocate(sources(3,ns),targvals(3,npts))

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)
      enddo
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,npts
        targvals(1,i) = srcvals(1,i)
        targvals(2,i) = srcvals(2,i)
        targvals(3,i) = srcvals(3,i)
      enddo
C$OMP END PARALLEL DO

c
c  get mean curvature
c
c

      call get_mean_curvature(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,curv)

c
c
c  oversample the density sigma
c

      call oversample_fun_surf(1,npatches,norders,ixyzs,iptype, 
     1    npts,sigma,novers,ixyzso,ns,sigmaover)
c
c  Set the charge and dipole densities for evaluating
c   
c
c     abc(1,:) = D[\sigma]
c     abc(2,:) =  ((S''+D')-2HS'-WS)[\sigma]
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,ns
        charges_tmp(1,i) = sigmaover(i)*whtsover(i)*over4pi
        charges_tmp(2,i) = 0

        dipvec_tmp(1,1,i) = 0
        dipvec_tmp(1,2,i) = 0
        dipvec_tmp(1,3,i) = 0

        dipvec_tmp(2,1,i) = 
     1     sigmaover(i)*whtsover(i)*srcover(10,i)*over4pi 
        dipvec_tmp(2,2,i) = 
     1     sigmaover(i)*whtsover(i)*srcover(11,i)*over4pi 
        dipvec_tmp(2,3,i) = 
     1     sigmaover(i)*whtsover(i)*srcover(12,i)*over4pi 
      enddo
C$OMP END PARALLEL DO


c
c   call the FMM 
c
c

      iper = 0
      ier = 0
      ifpgh = 0
      ifpghtarg = 3
      ifdipole = 1
      ifcharge = 1


      call lfmm3d(nd,eps,ns,sources,ifcharge,charges_tmp,
     1  ifdipole,dipvec_tmp,iper,ifpgh,tmp,tmp,tmp,npts,targvals,
     1  ifpghtarg,pot_tmp,grad_tmp,hess_tmp,ier)

c
c  Now start assembling components
c  s_pot will hold S[\sigma]
c  pot2 will hold 2*H*S'[\sigma]
c  pot3 will hold (S'' + D')[\sigma]
c  abc(1,:) will hold D[\sigma]
c 
      
      allocate(pot2(npts),pot3(npts))
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)       
      do i=1,npts
         s_pot(i) = pot_tmp(1,i)
         abc(1,i) = pot_tmp(2,i)
         pot2(i) = grad_tmp(1,1,i)*srcvals(10,i) + 
     1             grad_tmp(1,2,i)*srcvals(11,i) + 
     2             grad_tmp(1,3,i)*srcvals(12,i)


         pot3(i) = grad_tmp(2,1,i)*srcvals(10,i) + 
     1             grad_tmp(2,2,i)*srcvals(11,i) + 
     2             grad_tmp(2,3,i)*srcvals(12,i)

         pot3(i) = pot3(i) + 
     1        hess_tmp(1,1,i)*srcvals(10,i)*srcvals(10,i) + 
     1        hess_tmp(1,4,i)*srcvals(11,i)*srcvals(10,i) + 
     2        hess_tmp(1,5,i)*srcvals(12,i)*srcvals(10,i) +
     1        hess_tmp(1,4,i)*srcvals(10,i)*srcvals(11,i) +
     1        hess_tmp(1,2,i)*srcvals(11,i)*srcvals(11,i) +
     1        hess_tmp(1,6,i)*srcvals(12,i)*srcvals(11,i) +
     1        hess_tmp(1,5,i)*srcvals(10,i)*srcvals(12,i) +
     1        hess_tmp(1,6,i)*srcvals(11,i)*srcvals(12,i) +
     1        hess_tmp(1,3,i)*srcvals(12,i)*srcvals(12,i) 
      enddo
C$OMP END PARALLEL DO     
c
c  Fix quadrature corrections
c

      call cpu_time(t1)
C$      t1 = omp_get_wtime()

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart)
C$OMP$PRIVATE(jstart,npols,l)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch)
          do l=1,npols
            s_pot(i) = s_pot(i) + 
     1          wnear(jquadstart+l-1,1)*sigma(jstart+l-1)
            abc(1,i) = abc(1,i) + 
     1          wnear(jquadstart+l-1,2)*sigma(jstart+l-1)
            pot2(i) = pot2(i) + 
     1          wnear(jquadstart+l-1,3)*sigma(jstart+l-1)
            pot3(i) = pot3(i) + 
     1          wnear(jquadstart+l-1,4)*sigma(jstart+l-1)
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      

c      print *, "Near quad addition done"

      thresh = 0.0d0
      call get_fmm_thresh(3,ns,sources,3,npts,targvals,thresh)

c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
C$OMP$PRIVATE(ctmp2,dtmp2,nss,l,jstart,ii,val,vgrad,vhess,npover)
      do i=1,npts

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
            

            ctmp2(1,ii) = charges_tmp(1,jstart+l)
            ctmp2(2,ii) = charges_tmp(2,jstart+l)

            dtmp2(1,1,ii) = dipvec_tmp(1,1,jstart+l)
            dtmp2(1,2,ii) = dipvec_tmp(1,2,jstart+l)
            dtmp2(1,3,ii) = dipvec_tmp(1,3,jstart+l)

            dtmp2(2,1,ii) = dipvec_tmp(2,1,jstart+l)
            dtmp2(2,2,ii) = dipvec_tmp(2,2,jstart+l)
            dtmp2(2,3,ii) = dipvec_tmp(2,3,jstart+l)

          enddo
        enddo
        nss = ii

        val(1) = 0
        val(2) = 0

        vgrad(1:2,1:3) = 0
        vhess(1:2,1:6) = 0


        call l3ddirectcdh(nd,srctmp2,ctmp2,dtmp2,
     1        nss,targvals(1,i),ntarg0,val,vgrad,vhess,thresh)
       
        s_pot(i) = s_pot(i) - val(1)
        abc(1,i) = abc(1,i) - val(2)
        pot2(i) = pot2(i) - vgrad(1,1)*srcvals(10,i) -
     1      vgrad(1,2)*srcvals(11,i) - vgrad(1,3)*srcvals(12,i)
        pot3(i) = pot3(i) - vgrad(2,1)*srcvals(10,i) -
     1      vgrad(2,2)*srcvals(11,i) - vgrad(2,3)*srcvals(12,i)
        pot3(i) = pot3(i) - 
     1        (vhess(1,1)*srcvals(10,i)*srcvals(10,i) + 
     1        vhess(1,4)*srcvals(11,i)*srcvals(10,i) + 
     2        vhess(1,5)*srcvals(12,i)*srcvals(10,i) +
     1        vhess(1,4)*srcvals(10,i)*srcvals(11,i) +
     1        vhess(1,2)*srcvals(11,i)*srcvals(11,i) +
     1        vhess(1,6)*srcvals(12,i)*srcvals(11,i) +
     1        vhess(1,5)*srcvals(10,i)*srcvals(12,i) +
     1        vhess(1,6)*srcvals(11,i)*srcvals(12,i) +
     1        vhess(1,3)*srcvals(12,i)*srcvals(12,i)) 
      enddo
C$OMP END PARALLEL DO     

c      print *, "Subtraction done"
c
c    s_pot = S[\sigma]
c    abc(1,:) = D[\sigma]
c    pot2 = S'[\sigma]
c    pot3 = (S'' + D')[\sigma]
c
      ra = 0
      rsurf = 0
      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) REDUCTION(+:ra,rsurf)
      do i=1,npts
        ra = ra + s_pot(i)*wts(i)
        rsurf = rsurf + wts(i)
      enddo
C$OMP END PARALLEL DO
c
c  Now assemble
c     abc(2,:) =  -((S''+D')+2HS'-WS)[\sigma]
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,npts
        abc(2,i) = -(pot3(i) + 2*curv(i)*pot2(i) - ra/rsurf)
      enddo
C$OMP END PARALLEL DO


      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype, 
     1    npts,abc,novers,ixyzso,ns,abcover)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,ns
        charges_tmp(1,i) = abcover(2,i)*whtsover(i)*over4pi
        charges_tmp(2,i) = 0

        dipvec_tmp(1,1,i) = 0
        dipvec_tmp(1,2,i) = 0
        dipvec_tmp(1,3,i) = 0

        dipvec_tmp(2,1,i) = 
     1     abcover(1,i)*whtsover(i)*srcover(10,i)*over4pi 
        dipvec_tmp(2,2,i) = 
     1     abcover(1,i)*whtsover(i)*srcover(11,i)*over4pi 
        dipvec_tmp(2,3,i) = 
     1     abcover(1,i)*whtsover(i)*srcover(12,i)*over4pi 
      enddo
C$OMP END PARALLEL DO

c
c   call the FMM 
c
c

      iper = 0
      ier = 0
      ifpgh = 0
      ifpghtarg = 1
      ifdipole = 1
      ifcharge = 1
      call lfmm3d(nd,eps,ns,sources,ifcharge,charges_tmp,
     1  ifdipole,dipvec_tmp,iper,ifpgh,tmp,tmp,tmp,npts,targvals,
     1  ifpghtarg,pot_tmp,grad_tmp,hess_tmp,ier)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)       
      do i=1,npts
         pot(i) = pot_tmp(1,i) + pot_tmp(2,i)
      enddo
C$OMP END PARALLEL DO      

c
c  Fix quadrature corrections
c

      call cpu_time(t1)
C$      t1 = omp_get_wtime()

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart)
C$OMP$PRIVATE(jstart,l,npols)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            pot(i) = pot(i) + 
     1          wnear(jquadstart+l-1,1)*abc(2,jstart+l-1)
            pot(i) = pot(i) + 
     1          wnear(jquadstart+l-1,2)*abc(1,jstart+l-1)
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO
      

      call cpu_time(t2)
C$      t2 = omp_get_wtime()      

c      print *, "Near quad addition done"



c

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
C$OMP$PRIVATE(ctmp2,dtmp2,nss,l,jstart,ii,val,vgrad,vhess,npover)
      do i=1,npts

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
            
            ctmp2(1,ii) = charges_tmp(1,jstart+l)
            ctmp2(2,ii) = charges_tmp(2,jstart+l)

            dtmp2(1,1,ii) = dipvec_tmp(1,1,jstart+l)
            dtmp2(1,2,ii) = dipvec_tmp(1,2,jstart+l)
            dtmp2(1,3,ii) = dipvec_tmp(1,3,jstart+l)

            dtmp2(2,1,ii) = dipvec_tmp(2,1,jstart+l)
            dtmp2(2,2,ii) = dipvec_tmp(2,2,jstart+l)
            dtmp2(2,3,ii) = dipvec_tmp(2,3,jstart+l)

          enddo
        enddo
        nss = ii

        val(1) = 0
        val(2) = 0

        call l3ddirectcdp(nd,srctmp2,ctmp2,dtmp2,
     1        nss,targvals(1,i),ntarg0,val,thresh)
        pot(i) = pot(i)-val(1)-val(2)
      enddo
C$OMP END PARALLEL DO


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
      subroutine lpcomp_lap_bel_s_lap_s_noc_addsub(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,eps,nnz,row_ptr,col_ind,
     2   iquad,nquad,wnear,sigma,novers,nptso,ixyzso,srcover,whtsover,
     3   s_one,pot,s_pot)
c
c
c
c  This subroutine evaluates the layer potential for the integral
c  representation:
c
c			\nabla \cdot S \nabla_{\Gamma} S [\sigma] + SWS \sigma
c
c  where the near field is precomputed and stored in the row
c  sparse compressed format.
c
c  The subroutine also returns S[\sigma] which is the representation
c  for the solution
c
c  The fmm is used to accelerate the far-field with two calls to the
c  vector fmm, and the near-field interactions are 
c  handled via precomputed quadrature.
c
c  Using add and subtract: no need to call tree and set fmm parameters,
c  can call existing fmm library directly
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
c    - nquad: integer
c        number of near field entries corresponding to each source target
c        pair. The size of wnear is (nquad,4) since there are 4 kernels
c        per source target pair
c    - wnear: real *8(nquad,4)
c        The desired near field quadrature
c        * The first kernel is S
c        * The second kernel is \nabla_{1} S
c        * The third kernel is \nabla_{2} S
c        * The fourth kernel is \nabla_{3} S 
c    - sigma: real *8(npts)
c        density for the layer potential
c    - novers: integer(npatches)
c        order of discretization for oversampled sources and
c        density
c    - ixyzso: integer(npatches+1)
c        ixyzso(i) denotes the starting location in srcover,
c        corresponding to patch i
c    - nptso: integer
c        total number of oversampled points
c    - srcover: real *8 (12,nptso)
c        oversampled set of source information
c    - whtsover: real *8 (nptso)
c        smooth quadrature weights at oversampled nodes
c    - s_one: real *8 (npts)
c        S[1] on surface
c
c  Output arguments
c    - pot: real *8 (npts)
c        (S \Delta_{\Gamma} S + SWS) \sigma applied using Calderon
c        identities
c------------------------

      implicit none
      integer, intent(in) :: npatches,npts
      integer, intent(in) :: norders(npatches),ixyzs(npatches+1)
      integer, intent(in) :: ixyzso(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      integer, intent(in) :: nnz,row_ptr(npts+1),col_ind(nnz),nquad
      integer, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: wnear(nquad,4),sigma(npts)
      integer, intent(in) :: novers(npatches+1)
      integer, intent(in) :: nptso
      real *8, intent(in) :: srcover(12,nptso),whtsover(nptso)
      real *8, intent(in) :: s_one(npts)
      real *8, intent(out) :: pot(npts),s_pot(npts)

      integer norder,npols,nover,npolso

      integer ndtarg,ntarg

      real *8, allocatable :: sources(:,:),targvals(:,:)
      real *8, allocatable :: charges(:),dipvec(:,:),sigmaover(:)

      real *8, allocatable :: curv(:)
      real *8, allocatable :: abc(:,:),pot_tmp(:,:),grad_tmp(:,:,:)
      real *8, allocatable :: hess_tmp(:,:,:),charges_tmp(:,:)
      real *8, allocatable :: dipvec_tmp(:,:,:),abcover(:,:)
      real *8, allocatable :: ffforminv(:,:,:)
      real *8, allocatable :: wtmp1(:,:),wtmp2(:,:)
      integer ns,nt
      real *8 alpha,beta
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      real *8 tmp(10),val(3),grad(3),nvgrad,nvhess
      real *8 hess(6),vgrad(3,3),vhess(3,6)
      real *8 xmin,xmax,ymin,ymax,zmin,zmax,sizey,sizez,boxsize
      integer i,j,jpatch,jquadstart,jstart

      integer ifaddsub

      integer ntj
      
      real *8 ddot,pottmp
      real *8, allocatable :: ctmp2(:,:),dtmp2(:,:,:)
      real *8, allocatable :: wts(:)
      real *8 radexp,epsfmm

      integer ipars
      complex *16 zpars
      real *8 timeinfo(10),t1,t2,omp_get_wtime
      real *8 dpars(2)
      real *8 w1,w2,w3,u1,u2,uf,vf

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      real *8 thresh,ra
      real *8 rr,rmin,rint,rsurf
      integer nss,ii,l,npover

      integer nd,ntarg0,nmax
      integer ier,iper

      real *8 ttot,done,pi,over4pi

      parameter (ntarg0=1)



      ns = nptso
      done = 1.0d0
      pi = atan(done)*4.0d0
      over4pi = 1.0d0/4.0d0/pi

      nmax = 0
      call get_near_corr_max(npts,row_ptr,nnz,col_ind,npatches,
     1   ixyzso,nmax)
      allocate(srctmp2(3,nmax))

      nd = 1
      allocate(abc(3,npts),charges_tmp(nd,ns),dipvec_tmp(nd,3,ns))
      allocate(abcover(3,ns),sigmaover(ns),curv(npts))
      allocate(wtmp1(3,npts),wtmp2(3,npts),ffforminv(2,2,npts))

      allocate(ctmp2(nd,nmax),dtmp2(nd,3,nmax))
      allocate(pot_tmp(nd,npts),grad_tmp(nd,3,npts),hess_tmp(nd,6,npts))

c
c   Set source info and targinfo for the fmm
c
c
      allocate(sources(3,ns),targvals(3,npts))

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)
      enddo
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,npts
        targvals(1,i) = srcvals(1,i)
        targvals(2,i) = srcvals(2,i)
        targvals(3,i) = srcvals(3,i)
      enddo
C$OMP END PARALLEL DO

c
c  get mean curvature
c
c

      call get_mean_curvature(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,curv)
      call get_inv_first_fundamental_form(npatches,norders,ixyzs, 
     1   iptype,npts,srccoefs,srcvals,ffforminv)
      

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
         wtmp1(1:3,i) = ffforminv(1,1,i)*srcvals(4:6,i) + 
     1      ffforminv(1,2,i)*srcvals(7:9,i)
         wtmp2(1:3,i) = ffforminv(2,1,i)*srcvals(4:6,i) + 
     1      ffforminv(2,2,i)*srcvals(7:9,i)
      enddo
C$OMP END PARALLEL DO     
      

c
c
c  oversample the density sigma
c

      call oversample_fun_surf(1,npatches,norders,ixyzs,iptype, 
     1    npts,sigma,novers,ixyzso,ns,sigmaover)
c
c  Set the charge and dipole densities for evaluating
c   
c
c     abc(1:3,:) = \nabla_{\Gamma} S[\sigma]
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,ns
        charges_tmp(1,i) = sigmaover(i)*whtsover(i)*over4pi

        dipvec_tmp(1,1,i) = 0
        dipvec_tmp(1,2,i) = 0
        dipvec_tmp(1,3,i) = 0

      enddo
C$OMP END PARALLEL DO


c
c   call the FMM 
c
c

      iper = 0
      ier = 0
      ifpgh = 0
      ifpghtarg = 2
      ifcharge = 1
      ifdipole = 0


      call lfmm3d(nd,eps,ns,sources,ifcharge,charges_tmp,
     1  ifdipole,dipvec_tmp,iper,ifpgh,tmp,tmp,tmp,npts,targvals,
     1  ifpghtarg,pot_tmp,grad_tmp,hess_tmp,ier)

c
c  Now start assembling components
c  s_pot will hold S[\sigma]
c  abc(1:3,:) will hold \nabla_{\Gamma} S[\sigma]
c 
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,u1,u2)       
      do i=1,npts
         s_pot(i) = pot_tmp(1,i)
         call dot_prod3d(grad_tmp(1,1:3,i),srcvals(4,i),u1)
         call dot_prod3d(grad_tmp(1,1:3,i),srcvals(7,i),u2)
         abc(1:3,i) = u1*wtmp1(1:3,i) + u2*wtmp2(1:3,i)
      enddo
C$OMP END PARALLEL DO     
c
c  Fix quadrature corrections
c

      call cpu_time(t1)
C$      t1 = omp_get_wtime()

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart)
C$OMP$PRIVATE(jstart,npols,l,w1,w2,w3,uf,vf)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch)
          do l=1,npols
            s_pot(i) = s_pot(i) + 
     1          wnear(jquadstart+l-1,1)*sigma(jstart+l-1)
            w1 = wnear(jquadstart+l-1,2)
            w2 = wnear(jquadstart+l-1,3)
            w3 = wnear(jquadstart+l-1,4)
            uf = w1*srcvals(4,i) + w2*srcvals(5,i) + w3*srcvals(6,i)
            vf = w1*srcvals(7,i) + w2*srcvals(8,i) + w3*srcvals(9,i)
            abc(1:3,i) = abc(1:3,i) + 
     1         (uf*wtmp1(1:3,i)+vf*wtmp2(1:3,i))*sigma(jstart+l-1)
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      

c      print *, "Near quad addition done"

      thresh = 0.0d0
      call get_fmm_thresh(3,ns,sources,3,npts,targvals,thresh)

c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
C$OMP$PRIVATE(ctmp2,dtmp2,nss,l,jstart,ii,val,grad,u1,u2,npover)
      do i=1,npts

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
            

            ctmp2(1,ii) = charges_tmp(1,jstart+l)

          enddo
        enddo
        nss = ii

        val(1) = 0
        grad(1:3) = 0


        call l3ddirectcg(nd,srctmp2,ctmp2,
     1        nss,targvals(1,i),ntarg0,val,grad,thresh)
        call dot_prod3d(grad,srcvals(4,i),u1)
        call dot_prod3d(grad,srcvals(7,i),u2)

       
        s_pot(i) = s_pot(i) - val(1)
        abc(1:3,i) = abc(1:3,i) - (u1*wtmp1(1:3,i) + u2*wtmp2(1:3,i)) 
      enddo
C$OMP END PARALLEL DO     

c      print *, "Subtraction done"
c
c    s_pot = S[\sigma]
c    abc(1:3,:) = \nabla_{\Gamma} S[\sigma]
c
      ra = 0
      rsurf = 0
      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) REDUCTION(+:ra,rsurf)
      do i=1,npts
        ra = ra + s_pot(i)*wts(i)
        rsurf = rsurf + wts(i)
      enddo
C$OMP END PARALLEL DO

      nd = 3
      deallocate(ctmp2,dtmp2)
      deallocate(pot_tmp,grad_tmp,hess_tmp)
      deallocate(charges_tmp,dipvec_tmp)

      allocate(ctmp2(nd,nmax),dtmp2(nd,3,nmax))
      allocate(charges_tmp(nd,ns),dipvec_tmp(nd,3,ns))
      allocate(pot_tmp(nd,npts),grad_tmp(nd,3,npts),hess_tmp(nd,6,npts))


      call oversample_fun_surf(nd,npatches,norders,ixyzs,iptype, 
     1    npts,abc,novers,ixyzso,ns,abcover)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,ns
        charges_tmp(1:3,i) = abcover(1:3,i)*whtsover(i)*over4pi
        dipvec_tmp(1:3,1:3,i) = 0
      enddo
C$OMP END PARALLEL DO

c
c   call the FMM 
c
c

      iper = 0
      ier = 0
      ifpgh = 0
      ifpghtarg = 2
      ifcharge = 1
      ifdipole = 0
      call lfmm3d(nd,eps,ns,sources,ifcharge,charges_tmp,
     1  ifdipole,dipvec_tmp,iper,ifpgh,tmp,tmp,tmp,npts,targvals,
     1  ifpghtarg,pot_tmp,grad_tmp,hess_tmp,ier)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)       
      do i=1,npts
         pot(i) = grad_tmp(1,1,i) + grad_tmp(2,2,i) + grad_tmp(3,3,i)
      enddo
C$OMP END PARALLEL DO      

c
c  Fix quadrature corrections
c

      call cpu_time(t1)
C$      t1 = omp_get_wtime()

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart)
C$OMP$PRIVATE(jstart,l,npols,w1,w2,w3)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            w1 = wnear(jquadstart+l-1,2)
            w2 = wnear(jquadstart+l-1,3)
            w3 = wnear(jquadstart+l-1,4)
            pot(i) = pot(i) + w1*abc(1,jstart+l-1) +
     1         w2*abc(2,jstart+l-1) + w3*abc(3,jstart+l-1) 
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO
      

      call cpu_time(t2)
C$      t2 = omp_get_wtime()      

c      print *, "Near quad addition done"



c

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
C$OMP$PRIVATE(ctmp2,dtmp2,nss,l,jstart,ii,val,vgrad,npover)
      do i=1,npts

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
            
            ctmp2(1,ii) = charges_tmp(1,jstart+l)
            ctmp2(2,ii) = charges_tmp(2,jstart+l)
            ctmp2(3,ii) = charges_tmp(3,jstart+l)

          enddo
        enddo
        nss = ii

        val(1) = 0
        val(2) = 0
        val(3) = 0
        vgrad(1:3,1:3) = 0

        call l3ddirectcg(nd,srctmp2,ctmp2,
     1        nss,targvals(1,i),ntarg0,val,vgrad,thresh)
        pot(i) = pot(i)-vgrad(1,1)-vgrad(2,2)-vgrad(3,3)
c
c  add SWS now
c
        pot(i) = pot(i) + s_one(i)*ra/rsurf
      enddo
C$OMP END PARALLEL DO

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
      subroutine lpcomp_lap_bel_lap_s2_addsub(npatches,norders,ixyzs,
     1   iptype,npts,srccoefs,srcvals,eps,nnz,row_ptr,col_ind,iquad,
     2   nquad,wnear,sigma,novers,nptso,ixyzso,srcover,whtsover,pot,
     3   s2_pot)

c
c
c
c  This subroutine evaluates the layer potential for the integral
c  representation:
c
c			-sigma/4+(S')^2 - (S''+ D' +2HS')S\sigma + WS^2 \sigma
c
c  where the near field is precomputed and stored in the row
c  sparse compressed format.
c
c  The subroutine also returns S[\sigma] which is the representation
c  for the solution
c
c  The fmm is used to accelerate the far-field with two calls to the
c  vector fmm, and the near-field interactions are 
c  handled via precomputed quadrature.
c
c  Using add and subtract: no need to call tree and set fmm parameters,
c  can call existing fmm library directly
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
c    - nquad: integer
c        number of near field entries corresponding to each source target
c        pair. The size of wnear is (nquad,4) since there are 4 kernels
c        per source target pair
c    - wnear: real *8(nquad,4)
c        The desired near field quadrature
c        * The first kernel is S
c        * The second kernel is S'
c        * The third kernel is S'' + D'
c    - sigma: real *8(npts)
c        density for the layer potential
c    - novers: integer(npatches)
c        order of discretization for oversampled sources and
c        density
c    - ixyzso: integer(npatches+1)
c        ixyzso(i) denotes the starting location in srcover,
c        corresponding to patch i
c    - nptso: integer
c        total number of oversampled points
c    - srcover: real *8 (12,nptso)
c        oversampled set of source information
c    - whtsover: real *8 (nptso)
c        smooth quadrature weights at oversampled nodes
c
c  Output arguments
c    - pot: real *8 (npts)
c        (S \Delta_{\Gamma} S + SWS) \sigma applied using Calderon
c        identities
c    - s2_pot: real *8 (npts)
c        The solution u:=S^2 \sigma 
c------------------------

      implicit none
      integer, intent(in) :: npatches,npts
      integer, intent(in) :: norders(npatches),ixyzs(npatches+1)
      integer, intent(in) :: ixyzso(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      integer, intent(in) :: nnz,row_ptr(npts+1),col_ind(nnz),nquad
      integer, intent(in) :: iquad(nnz+1)
      real *8, intent(in) :: wnear(nquad,3),sigma(npts)
      integer, intent(in) :: novers(npatches+1)
      integer, intent(in) :: nptso
      real *8, intent(in) :: srcover(12,nptso),whtsover(nptso)
      real *8, intent(out) :: pot(npts),s2_pot(npts)

      integer norder,npols,nover,npolso

      integer ndtarg,ntarg

      real *8, allocatable :: sources(:,:),targvals(:,:)
      real *8, allocatable :: charges(:),dipvec(:,:),sigmaover(:)

      real *8, allocatable :: pot2(:)
      real *8, allocatable :: pot3(:)
      real *8, allocatable :: curv(:)
      real *8, allocatable :: abc(:,:),pot_tmp(:,:),grad_tmp(:,:,:)
      real *8, allocatable :: hess_tmp(:,:,:),charges_tmp(:,:)
      real *8, allocatable :: dipvec_tmp(:,:,:),abcover(:,:)
      integer ns,nt
      real *8 alpha,beta
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      real *8 tmp(10),val(3),grad(3),nvgrad,nvhess
      real *8 hess(6),vgrad(3,3),vhess(3,6)
      real *8 xmin,xmax,ymin,ymax,zmin,zmax,sizey,sizez,boxsize
      integer i,j,jpatch,jquadstart,jstart

      integer ifaddsub

      integer ntj
      
      real *8 ddot,pottmp
      real *8, allocatable :: ctmp2(:,:),dtmp2(:,:,:)
      real *8, allocatable :: wts(:)
      real *8 radexp,epsfmm

      integer ipars
      complex *16 zpars
      real *8 timeinfo(10),t1,t2,omp_get_wtime
      real *8 dpars(2)

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      real *8 thresh,ra
      real *8 rr,rmin,rint,rsurf
      integer nss,ii,l,npover

      integer nd,ntarg0,nmax
      integer ier,iper

      real *8 ttot,done,pi,over4pi

      parameter (ntarg0=1)



      ns = nptso
      done = 1.0d0
      pi = atan(done)*4.0d0
      over4pi = 1.0d0/4.0d0/pi

      nmax = 0
      call get_near_corr_max(npts,row_ptr,nnz,col_ind,npatches,
     1   ixyzso,nmax)
      allocate(srctmp2(3,nmax))

      nd = 1
      allocate(abc(2,npts),abcover(2,ns))
      allocate(sigmaover(ns),curv(npts))

      allocate(charges_tmp(nd,ns),dipvec_tmp(nd,3,ns))
      allocate(pot_tmp(nd,npts),grad_tmp(nd,3,npts),hess_tmp(nd,6,npts))

c
c   Set source info and targinfo for the fmm
c
c
      allocate(sources(3,ns),targvals(3,npts))

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)
      enddo
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,npts
        targvals(1,i) = srcvals(1,i)
        targvals(2,i) = srcvals(2,i)
        targvals(3,i) = srcvals(3,i)
      enddo
C$OMP END PARALLEL DO

c
c  get mean curvature
c
c

      call get_mean_curvature(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,curv)

c
c
c  oversample the density sigma
c

      call oversample_fun_surf(1,npatches,norders,ixyzs,iptype, 
     1    npts,sigma,novers,ixyzso,ns,sigmaover)
c
c  Set the charge and dipole densities for evaluating
c   
c
c     abc(1,:) = S[\tau] 
c     abc(2,:) = S'[\tau] 
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,ns
        charges_tmp(1,i) = sigmaover(i)*whtsover(i)*over4pi

        dipvec_tmp(1,1,i) = 0
        dipvec_tmp(1,2,i) = 0
        dipvec_tmp(1,3,i) = 0
      enddo
C$OMP END PARALLEL DO


c
c   call the FMM 
c
c

      iper = 0
      ier = 0
      ifpgh = 0
      ifpghtarg = 2
      ifdipole = 0
      ifcharge = 1


      call lfmm3d(nd,eps,ns,sources,ifcharge,charges_tmp,
     1  ifdipole,dipvec_tmp,iper,ifpgh,tmp,tmp,tmp,npts,targvals,
     1  ifpghtarg,pot_tmp,grad_tmp,hess_tmp,ier)

c
c  Now start assembling components
c  abc(1,:) will hold S[\sigma]
c  abc(2,:) will hold S'[\sigma]
c 
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)       
      do i=1,npts
         abc(1,i) = pot_tmp(1,i)
         abc(2,i) = grad_tmp(1,1,i)*srcvals(10,i) + 
     1             grad_tmp(1,2,i)*srcvals(11,i) + 
     2             grad_tmp(1,3,i)*srcvals(12,i)
      enddo
C$OMP END PARALLEL DO     
c
c  Fix quadrature corrections
c

      call cpu_time(t1)
C$      t1 = omp_get_wtime()

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart)
C$OMP$PRIVATE(jstart,npols,l)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch)
          do l=1,npols
            abc(1,i) = abc(1,i) + 
     1          wnear(jquadstart+l-1,1)*sigma(jstart+l-1)
            abc(2,i) = abc(2,i) + 
     1          wnear(jquadstart+l-1,2)*sigma(jstart+l-1)
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      

c      print *, "Near quad addition done"

      thresh = 0.0d0
      call get_fmm_thresh(3,ns,sources,3,npts,targvals,thresh)
      allocate(ctmp2(nd,nmax),dtmp2(nd,3,nmax))

c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
C$OMP$PRIVATE(ctmp2,dtmp2,nss,l,jstart,ii,val,grad,npover)
      do i=1,npts

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
            

            ctmp2(1,ii) = charges_tmp(1,jstart+l)

          enddo
        enddo
        nss = ii

        val(1) = 0
        val(2) = 0
        val(3) = 0

        grad(1:3) = 0


        call l3ddirectcg(nd,srctmp2,ctmp2,
     1        nss,targvals(1,i),ntarg0,val,grad,thresh)
       
        abc(1,i) = abc(1,i) - val(1)
        abc(2,i) = abc(2,i) - (grad(1)*srcvals(10,i) + 
     1      grad(2)*srcvals(11,i)+grad(3)*srcvals(12,i))
      enddo
C$OMP END PARALLEL DO     

c      print *, "Subtraction done"
c
c    abc(1,:) = S[\sigma]
c    abc(2,:) = S'[\sigma]
c

      deallocate(charges_tmp,dipvec_tmp,pot_tmp,grad_tmp,hess_tmp)
      deallocate(ctmp2,dtmp2)

      nd = 3
      allocate(charges_tmp(nd,ns),dipvec_tmp(nd,3,ns))
      allocate(pot_tmp(nd,npts),grad_tmp(nd,3,npts),hess_tmp(nd,6,npts))
      allocate(ctmp2(nd,nmax),dtmp2(nd,3,nmax))


      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype, 
     1    npts,abc,novers,ixyzso,ns,abcover)
      

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,ns
        charges_tmp(1,i) = abcover(2,i)*whtsover(i)*over4pi
        charges_tmp(2,i) = abcover(1,i)*whtsover(i)*over4pi
        charges_tmp(3,i) = 0

        dipvec_tmp(1,1,i) = 0
        dipvec_tmp(1,2,i) = 0
        dipvec_tmp(1,3,i) = 0

        dipvec_tmp(2,1,i) = 0
        dipvec_tmp(2,2,i) = 0
        dipvec_tmp(2,3,i) = 0

        dipvec_tmp(3,1,i) = 
     1     abcover(1,i)*whtsover(i)*srcover(10,i)*over4pi 
        dipvec_tmp(3,2,i) = 
     1     abcover(1,i)*whtsover(i)*srcover(11,i)*over4pi 
        dipvec_tmp(3,3,i) = 
     1     abcover(1,i)*whtsover(i)*srcover(12,i)*over4pi 
      enddo
C$OMP END PARALLEL DO

c
c   call the FMM 
c
c

      iper = 0
      ier = 0
      ifpgh = 0
      ifpghtarg = 3
      ifdipole = 1
      ifcharge = 1
      call lfmm3d(nd,eps,ns,sources,ifcharge,charges_tmp,
     1  ifdipole,dipvec_tmp,iper,ifpgh,tmp,tmp,tmp,npts,targvals,
     1  ifpghtarg,pot_tmp,grad_tmp,hess_tmp,ier)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)       
      do i=1,npts

c
c   add contribution of S'2[\sigma]
c
         pot(i) = grad_tmp(1,1,i)*srcvals(10,i) + 
     1       grad_tmp(1,2,i)*srcvals(11,i) +
     2       grad_tmp(1,3,i)*srcvals(12,i)
c
c   add contribution of 2HS'S\sigma
c
         pot(i) = pot(i) - 2*curv(i)*(grad_tmp(2,1,i)*srcvals(10,i) + 
     1        grad_tmp(2,2,i)*srcvals(11,i) + 
     2        grad_tmp(2,3,i)*srcvals(12,i))
c
c  add contribution of -S''S\sigma - D' S\sigma
c
c

         pot(i) = pot(i) - (grad_tmp(3,1,i)*srcvals(10,i) + 
     1             grad_tmp(3,2,i)*srcvals(11,i) + 
     2             grad_tmp(3,3,i)*srcvals(12,i))

         pot(i) = pot(i) - 
     1        (hess_tmp(2,1,i)*srcvals(10,i)*srcvals(10,i) + 
     1        hess_tmp(2,4,i)*srcvals(11,i)*srcvals(10,i) + 
     2        hess_tmp(2,5,i)*srcvals(12,i)*srcvals(10,i) +
     1        hess_tmp(2,4,i)*srcvals(10,i)*srcvals(11,i) +
     1        hess_tmp(2,2,i)*srcvals(11,i)*srcvals(11,i) +
     1        hess_tmp(2,6,i)*srcvals(12,i)*srcvals(11,i) +
     1        hess_tmp(2,5,i)*srcvals(10,i)*srcvals(12,i) +
     1        hess_tmp(2,6,i)*srcvals(11,i)*srcvals(12,i) +
     1        hess_tmp(2,3,i)*srcvals(12,i)*srcvals(12,i))


         s2_pot(i) = pot_tmp(2,i)

      enddo
C$OMP END PARALLEL DO      

c
c  Fix quadrature corrections
c

      call cpu_time(t1)
C$      t1 = omp_get_wtime()

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart)
C$OMP$PRIVATE(jstart,l,npols)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            pot(i) = pot(i) + 
     1          wnear(jquadstart+l-1,2)*abc(2,jstart+l-1)
            pot(i) = pot(i) - 
     1          2*curv(i)*wnear(jquadstart+l-1,2)*abc(1,jstart+l-1)
            pot(i) = pot(i) - wnear(jquadstart+l-1,3)*abc(1,jstart+l-1)
            s2_pot(i) = s2_pot(i) +
     1         wnear(jquadstart+l-1,1)*abc(1,jstart+l-1)
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO
      

      call cpu_time(t2)
C$      t2 = omp_get_wtime()      

c      print *, "Near quad addition done"



c

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
C$OMP$PRIVATE(ctmp2,dtmp2,nss,l,jstart,ii,val,vgrad,vhess,npover)
      do i=1,npts

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
            
            ctmp2(1,ii) = charges_tmp(1,jstart+l)
            ctmp2(2,ii) = charges_tmp(2,jstart+l)
            ctmp2(3,ii) = charges_tmp(3,jstart+l)

            dtmp2(1,1,ii) = dipvec_tmp(1,1,jstart+l)
            dtmp2(1,2,ii) = dipvec_tmp(1,2,jstart+l)
            dtmp2(1,3,ii) = dipvec_tmp(1,3,jstart+l)

            dtmp2(2,1,ii) = dipvec_tmp(2,1,jstart+l)
            dtmp2(2,2,ii) = dipvec_tmp(2,2,jstart+l)
            dtmp2(2,3,ii) = dipvec_tmp(2,3,jstart+l)

            dtmp2(3,1,ii) = dipvec_tmp(3,1,jstart+l)
            dtmp2(3,2,ii) = dipvec_tmp(3,2,jstart+l)
            dtmp2(3,3,ii) = dipvec_tmp(3,3,jstart+l)

          enddo
        enddo
        nss = ii

        val(1) = 0
        val(2) = 0
        val(3) = 0
        vgrad(1:3,1:3) = 0
        vhess(1:3,1:6) = 0

        call l3ddirectcdh(nd,srctmp2,ctmp2,dtmp2,
     1        nss,targvals(1,i),ntarg0,val,vgrad,vhess,thresh)
        

        pot(i) = pot(i) - vgrad(1,1)*srcvals(10,i) -
     1      vgrad(1,2)*srcvals(11,i) - vgrad(1,3)*srcvals(12,i)
        
        pot(i) = pot(i) + 2*curv(i)*(vgrad(2,1)*srcvals(10,i) +
     1      vgrad(2,2)*srcvals(11,i) + vgrad(2,3)*srcvals(12,i))
        
        pot(i) = pot(i) + vgrad(3,1)*srcvals(10,i) +
     1      vgrad(3,2)*srcvals(11,i) + vgrad(3,3)*srcvals(12,i)

        pot(i) = pot(i) + 
     1        (vhess(2,1)*srcvals(10,i)*srcvals(10,i) + 
     1        vhess(2,4)*srcvals(11,i)*srcvals(10,i) + 
     2        vhess(2,5)*srcvals(12,i)*srcvals(10,i) +
     1        vhess(2,4)*srcvals(10,i)*srcvals(11,i) +
     1        vhess(2,2)*srcvals(11,i)*srcvals(11,i) +
     1        vhess(2,6)*srcvals(12,i)*srcvals(11,i) +
     1        vhess(2,5)*srcvals(10,i)*srcvals(12,i) +
     1        vhess(2,6)*srcvals(11,i)*srcvals(12,i) +
     1        vhess(2,3)*srcvals(12,i)*srcvals(12,i))
        s2_pot(i) = s2_pot(i) - val(2)
      enddo
C$OMP END PARALLEL DO

c
c  now subtract WS^2 from pot
c
      ra = 0
      rsurf = 0
      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) REDUCTION(+:ra,rsurf)
      do i=1,npts
        ra = ra + s2_pot(i)*wts(i)
        rsurf = rsurf + wts(i)
      enddo
C$OMP END PARALLEL DO
  
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,npts
        pot(i) = pot(i) + ra/rsurf
      enddo
C$OMP END PARALLEL DO

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





      subroutine lpcomp_lap_bel(npatches,norders,ixyzs,
     1    iptype,npts,srccoefs,srcvals,eps,sigma,irep,
     2    pot,u,t_mv_fmm)
c
c  This subroutine evaluates the Laplace Beltrami problem
c  using one of the following three integral representations:
c
c  irep=1, s_lap_s: pot = S (\Delta_{\Gamma} + W)S[\sigma], where
c                      S \Delta_{\Gamma}S [\sigma] is expanded out using 
c                      Calderon identities, and u=S[\sigma]
c  irep=2, s_lap_s_noc: Same as above but no Calderon identities
c            used
c  irep=3, lap_s2: pot  = (\Delta_{\Gamma} + W) S^2 [\sigma], where
c                    \Delta_{\Gamma}S^2 is expanded out
c                    using Calderon identities, u = S^2[\sigma]
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
c    - sigma: real *8 (npts)
c        the density sigma
c    - irep: integer
c        representation for solving the Laplace Beltrami problem
c        * irep = 1, s_lap_s
c        * irep = 2, s_lap_s_noc
c        * irep = 3, lap_s2
c
c  Output arguments:
c    - pot: real *8 (npts)
c        the potential pot above
c    - u: real *8 (npts)
c        the potential u above
c    - t_mv_fmm: real *8
c        time taken just for the fmm part
c        
c
      implicit none

      integer, intent(in) :: npatches,npts
      integer, intent(in) :: norders(npatches),ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      integer, intent(in) :: irep
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      real *8, intent(in) :: sigma(npts)
      real *8, intent(out) :: pot(npts),u(npts),t_mv_fmm
      

      real *8, allocatable :: targs(:,:)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      integer ndtarg,ntarg

      integer norder,npols
      integer nover,npolso,nptso
      integer nnz,nquad
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      real *8, allocatable :: wnear(:,:)
      real *8, allocatable :: s_one(:),dens_one(:)

      real *8, allocatable :: srcover(:,:),wover(:)
      integer, allocatable :: ixyzso(:),novers(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:) 

      integer i,j,jpatch,jquadstart,jstart

      integer ipars
      complex *16 zpars
      real *8 timeinfo(10),t1,t2,omp_get_wtime


      real *8 ttot,done,pi
      real *8 rfac,rfac0
      real *8 dpars(2),erra,ra
      integer iptype_avg,norder_avg
      integer ikerorder, iquadtype,npts_over

c
c
c       gmres variables
c
      real *8 did,dtmp
      integer it,iind,it1,k,l,ndquad
      real *8 rmyerr
      real *8 temp
      real *8, allocatable :: wts(:)

      complex *16 ztmp


      done = 1
      pi = atan(done)*4

      if(irep.lt.1.or.irep.gt.3) then
        print *, "invalid argument for irep, returning"
        return
      endif


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
      call findnearmem(cms,npatches,rad_near,3,targs,npts,nnz)
      print *, "nnz=",nnz

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,3,targs,npts,row_ptr, 
     1        col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind,
     1         iquad)


c
c    estimate oversampling for far-field, and oversample geometry
c

      allocate(novers(npatches),ixyzso(npatches+1))

      print *, "beginning far order estimation"

      ztmp = 0
      ikerorder = 0
      print *, eps

      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,
     1    rads,npts,srccoefs,ndtarg,npts,targs,ikerorder,ztmp,
     2    nnz,row_ptr,col_ind,rfac,novers,ixyzso)

      npts_over = ixyzso(npatches+1)-1
      print *, "npts_over=",npts_over

      allocate(srcover(12,npts_over),wover(npts_over))

      call oversample_geom(npatches,norders,ixyzs,iptype,npts, 
     1   srccoefs,srcvals,novers,ixyzso,npts_over,srcover)

      allocate(wts(npts))

      call get_qwts(npatches,norders,ixyzs,iptype,npts,
     1        srcvals,wts)

      call get_qwts(npatches,novers,ixyzso,iptype,npts_over,
     1        srcover,wover)

c
c   compute near quadrature correction
c
      nquad = iquad(nnz+1)-1
      print *, "nquad=",nquad
      if(irep.eq.1.or.irep.eq.2) then
        allocate(wnear(nquad,4))
        ndquad = 4
      else
        allocate(wnear(nquad,3))
        ndquad = 3
      endif
      
      do j=1,ndquad
C$OMP PARALLEL DO DEFAULT(SHARED)      
        do i=1,nquad
          wnear(i,j) = 0
        enddo
C$OMP END PARALLEL DO    
      enddo


      iquadtype = 1

      print *, "starting to generate near quadrature"
      call cpu_time(t1)
C$      t1 = omp_get_wtime()      
      if(irep.eq.1) then
        call getnearquad_lap_bel_s_lap_s(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,eps,
     2   iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear)
      else if(irep.eq.2) then
        call getnearquad_lap_bel_s_lap_s_noc(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,eps,
     2   iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear)
      else
        call getnearquad_lap_bel_lap_s2(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,eps,
     2   iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear)
      endif
      call cpu_time(t2)
C$      t2 = omp_get_wtime()     

      call prin2('quadrature generation time=*',t2-t1,1)
      print *, "done generating near quadrature"

      print *, "Starting to generate right hand side for linear system"


      if(irep.eq.2) then
        allocate(s_one(npts),dens_one(npts))

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,npts
          dens_one(i) = 1.0d0
        enddo
C$OMP END PARALLEL DO

        dpars(1) = 1.0d0
        dpars(2) = 0
        call lpcomp_lap_comb_dir_addsub(npatches,norders,ixyzs,
     1    iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,eps,
     2    dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear(1,1),
     3    dens_one,novers,npts_over,ixyzso,srcover,wover,s_one)
        
      endif

      print *, "Evaluate layer potential now"

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,npts
        u(i) = 0
        pot(i) = 0
      enddo
C$OMP END PARALLEL DO


      if(irep.eq.1.or.irep.eq.3) then
        did = -0.25d0
      else
        did = 0.0d0
      endif

      call cpu_time(t1)
C$      t1 = omp_get_wtime()      

      if(irep.eq.1) then
        call lpcomp_lap_bel_s_lap_s_addsub(npatches,norders,ixyzs,
     1   iptype,npts,srccoefs,srcvals,eps,nnz,row_ptr,col_ind,iquad,
     2   nquad,wnear,sigma,novers,npts_over,ixyzso,
     3   srcover,wover,pot,u)
      endif
      if(irep.eq.2) then
        call lpcomp_lap_bel_s_lap_s_noc_addsub(npatches,norders,ixyzs,
     1    iptype,npts,srccoefs,srcvals,eps,nnz,row_ptr,col_ind,iquad,
     2    nquad,wnear,sigma,novers,npts_over,ixyzso,
     3    srcover,wover,s_one,pot,u)
      endif
      if(irep.eq.3) then
        call lpcomp_lap_bel_lap_s2_addsub(npatches,norders,ixyzs,
     1    iptype,npts,srccoefs,srcvals,eps,nnz,row_ptr,col_ind,iquad,
     2    nquad,wnear,sigma,novers,npts_over,ixyzso,
     3    srcover,wover,pot,u)
      endif
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,npts
         pot(i) = pot(i) + did*sigma(i)
      enddo
C$OMP END PARALLEL DO

      call cpu_time(t2)
C$      t2 = omp_get_wtime()      
      t_mv_fmm = t2-t1

c
      return
      end
c
c
c
c
c
c

      subroutine lap_bel_solver(npatches,norders,ixyzs,
     1    iptype,npts,srccoefs,srcvals,eps,numit,rhs,irep,
     2    eps_gmres,niter,errs,rres,soln,u)
c
c  This subroutine solves the Laplace Beltrami problem
c  using one of the following three integral representations:
c
c  irep=1, s_lap_s: S (\Delta_{\Gamma} + W)S[\sigma] = S[f], where
c                      S \Delta_{\Gamma}S [\sigma] is expanded out using 
c                      Calderon identities
c  irep=2, s_lap_s_noc: Same as above but no Calderon identities
c            used
c  irep=3, lap_s2: \Delta_{\Gamma} + W S^2 [\sigma], where
c                    \Delta_{\Gamma}S^2 is expanded out
c                    using Calderon identities
c
c  The linear system is solved iteratively using GMRES.
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
c    - numit: integer
c        maximum number of GMRES iterations
c    - rhs: real *8 (npts)
c        the surface data f
c    - irep: integer
c        representation for solving the Laplace Beltrami problem
c        * irep = 1, s_lap_s
c        * irep = 2, s_lap_s_noc
c        * irep = 3, lap_s2
c    - eps_gmres: real *8
c        relative residual tolerance for GMRES
c
c  Output arguments:
c    - niter: integer
c        number of GMRES iterations used
c    - errs: real *8 (niter+1)
c        On input must be of size (numit+1), gmres residual
c        as a function of iteration number
c    - rres: real *8
c        relative residual to which the linear system is solved
c    - soln: real *8 (npts)
c        the solution sigma
c    - u: real *8 (npts)
c        The solution to the Laplace Beltrami problem
c
      implicit none

      integer, intent(in) :: npatches,npts
      integer, intent(in) :: norders(npatches),ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      integer, intent(in) :: numit,irep
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      real *8, intent(in) :: eps_gmres
      real *8, intent(in) :: rhs(npts)
      real *8, intent(out) :: soln(npts),u(npts)
      real *8, intent(out) :: errs(numit+1)
      real *8, intent(out) :: rres
      integer, intent(out) ::  niter
      

      real *8, allocatable :: rhs_use(:)
      real *8, allocatable :: targs(:,:)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      integer ndtarg,ntarg



      integer norder,npols
      integer nover,npolso,nptso
      integer nnz,nquad
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      real *8, allocatable :: wnear(:,:)
      real *8, allocatable :: s_one(:),dens_one(:)

      real *8, allocatable :: srcover(:,:),wover(:)
      integer, allocatable :: ixyzso(:),novers(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:) 

      integer i,j,jpatch,jquadstart,jstart

      integer ipars
      complex *16 zpars
      real *8 timeinfo(10),t1,t2,omp_get_wtime


      real *8 ttot,done,pi
      real *8 rfac,rfac0
      real *8 dpars(2),erra,ra
      integer iptype_avg,norder_avg
      integer ikerorder, iquadtype,npts_over

c
c
c       gmres variables
c
      real *8 did,dtmp
      real *8 rb,wnrm2
      integer it,iind,it1,k,l,ndquad
      real *8 rmyerr
      real *8 temp
      real *8, allocatable :: vmat(:,:),hmat(:,:)
      real *8, allocatable :: cs(:),sn(:)
      real *8, allocatable :: svec(:),yvec(:),wtmp(:)
      real *8, allocatable :: wts(:)

      complex *16 ztmp


      allocate(vmat(npts,numit+1),hmat(numit,numit))
      allocate(cs(numit),sn(numit))
      allocate(wtmp(npts),svec(numit+1),yvec(numit+1))


      done = 1
      pi = atan(done)*4

      if(irep.lt.1.or.irep.gt.3) then
        print *, "invalid argument for irep, returning"
        return
      endif


c
c
c        setup targets as on surface discretization points
c 
      ndtarg = 3
      ntarg = npts
      allocate(targs(ndtarg,npts),uvs_targ(2,ntarg),ipatch_id(ntarg))

      allocate(rhs_use(npts))

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
      call findnearmem(cms,npatches,rad_near,3,targs,npts,nnz)
      print *, "nnz=",nnz

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,3,targs,npts,row_ptr, 
     1        col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind,
     1         iquad)


c
c    estimate oversampling for far-field, and oversample geometry
c

      allocate(novers(npatches),ixyzso(npatches+1))

      print *, "beginning far order estimation"

      ztmp = 0
      ikerorder = 0
      print *, eps

      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,
     1    rads,npts,srccoefs,ndtarg,npts,targs,ikerorder,ztmp,
     2    nnz,row_ptr,col_ind,rfac,novers,ixyzso)

      npts_over = ixyzso(npatches+1)-1
      print *, "npts_over=",npts_over

      allocate(srcover(12,npts_over),wover(npts_over))

      call oversample_geom(npatches,norders,ixyzs,iptype,npts, 
     1   srccoefs,srcvals,novers,ixyzso,npts_over,srcover)

      allocate(wts(npts))

      call get_qwts(npatches,norders,ixyzs,iptype,npts,
     1        srcvals,wts)

      call get_qwts(npatches,novers,ixyzso,iptype,npts_over,
     1        srcover,wover)

c
c   compute near quadrature correction
c
      nquad = iquad(nnz+1)-1
      print *, "nquad=",nquad
      if(irep.eq.1.or.irep.eq.2) then
        allocate(wnear(nquad,4))
        ndquad = 4
      else
        allocate(wnear(nquad,3))
        ndquad = 3
      endif
      
      do j=1,ndquad
C$OMP PARALLEL DO DEFAULT(SHARED)      
        do i=1,nquad
          wnear(i,j) = 0
        enddo
C$OMP END PARALLEL DO    
      enddo


      iquadtype = 1

      print *, "starting to generate near quadrature"
      call cpu_time(t1)
C$      t1 = omp_get_wtime()      
      if(irep.eq.1) then
        call getnearquad_lap_bel_s_lap_s(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,eps,
     2   iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear)
      else if(irep.eq.2) then
        call getnearquad_lap_bel_s_lap_s_noc(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,eps,
     2   iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear)
      else
        call getnearquad_lap_bel_lap_s2(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,eps,
     2   iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear)
      endif
      call cpu_time(t2)
C$      t2 = omp_get_wtime()     

      call prin2('quadrature generation time=*',t2-t1,1)
      print *, "done generating near quadrature"

      print *, "Starting to generate right hand side for linear system"

c
c
c     compute the right hand side S[f], if irep=1, or irep=2
c     note that
c
      if(irep.eq.1.or.irep.eq.2) then
        dpars(1) = 1.0d0
        dpars(2) = 0
        call lpcomp_lap_comb_dir_addsub(npatches,norders,ixyzs,
     1    iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,eps,
     2    dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear(1,1),
     3    rhs,novers,npts_over,ixyzso,srcover,wover,rhs_use)
      else
        print  *, "here1"
C$OMP PARALLEL DO DEFAULT(SHARED)
        do i=1,npts
          rhs_use(i) = rhs(i)
        enddo
C$OMP END PARALLEL DO
      endif
      print *, "Done generating right hand side for linear system"

      if(irep.eq.2) then
        allocate(s_one(npts),dens_one(npts))

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,npts
          dens_one(i) = 1.0d0
        enddo
C$OMP END PARALLEL DO

        dpars(1) = 1.0d0
        dpars(2) = 0
        call lpcomp_lap_comb_dir_addsub(npatches,norders,ixyzs,
     1    iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,eps,
     2    dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear(1,1),
     3    dens_one,novers,npts_over,ixyzso,srcover,wover,s_one)
        

      endif

      print *, "Start gmres now"
c
c
c     start gmres code here
c
c     NOTE: matrix equation should be of the form (z*I + K)x = y
c       the identity scaling (z) is defined via did below,
c       and K represents the action of the principal value 
c       part of the matvec
c

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,npts
        u(i) = 0
      enddo
C$OMP END PARALLEL DO


      if(irep.eq.1.or.irep.eq.3) then
        did = -0.25d0
      else
        did = 0.0d0
      endif


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
        rb = rb + abs(rhs_use(i))**2*wts(i)
      enddo
      rb = sqrt(rb)

      do i=1,npts
        vmat(i,1) = rhs_use(i)/rb
      enddo

      svec(1) = rb

      do it=1,numit
        it1 = it + 1

c
c        NOTE:
c        replace this routine by appropriate layer potential
c        evaluation routine  
c
        if(irep.eq.1) then
          call lpcomp_lap_bel_s_lap_s_addsub(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,eps,nnz,row_ptr,col_ind,iquad,
     2     nquad,wnear,vmat(1,it),novers,npts_over,ixyzso,
     3     srcover,wover,wtmp,u)
        endif
        if(irep.eq.2) then
          call lpcomp_lap_bel_s_lap_s_noc_addsub(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,eps,nnz,row_ptr,col_ind,iquad,
     2     nquad,wnear,vmat(1,it),novers,npts_over,ixyzso,
     3     srcover,wover,s_one,wtmp,u)
        endif
        if(irep.eq.3) then
          call lpcomp_lap_bel_lap_s2_addsub(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,eps,nnz,row_ptr,col_ind,iquad,
     2     nquad,wnear,vmat(1,it),novers,npts_over,ixyzso,
     3     srcover,wover,wtmp,u)
        endif


        do k=1,it
          hmat(k,it) = 0
          do j=1,npts
            hmat(k,it) = hmat(k,it) + wtmp(j)*vmat(j,k)*wts(j)
          enddo

          do j=1,npts
            wtmp(j) = wtmp(j)-hmat(k,it)*vmat(j,k)
          enddo
        enddo
          
        hmat(it,it) = hmat(it,it)+did
        wnrm2 = 0
        do j=1,npts
          wnrm2 = wnrm2 + abs(wtmp(j))**2*wts(j)
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
          if(irep.eq.1) then
            call lpcomp_lap_bel_s_lap_s_addsub(npatches,norders,ixyzs,
     1       iptype,npts,srccoefs,srcvals,eps,nnz,row_ptr,col_ind,iquad,
     2       nquad,wnear,soln,novers,npts_over,ixyzso,
     3       srcover,wover,wtmp,u)
          endif
          if(irep.eq.2) then
            call lpcomp_lap_bel_s_lap_s_noc_addsub(npatches,norders,
     1       ixyzs,iptype,npts,srccoefs,srcvals,eps,nnz,row_ptr,col_ind,
     2       iquad,nquad,wnear,soln,novers,npts_over,ixyzso,
     3       srcover,wover,s_one,wtmp,u)
          endif

          if(irep.eq.3) then
            call lpcomp_lap_bel_lap_s2_addsub(npatches,norders,ixyzs,
     1       iptype,npts,srccoefs,srcvals,eps,nnz,row_ptr,col_ind,iquad,
     2       nquad,wnear,soln,novers,npts_over,ixyzso,
     3       srcover,wover,wtmp,u)
          endif

          do i=1,npts
            rres = rres + abs(did*soln(i) + wtmp(i)-rhs(i))**2*wts(i)
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
c

      subroutine lap_bel_solver_wquad(npatches,norders,ixyzs,
     1    iptype,npts,srccoefs,srcvals,eps,numit,rhs,irep,
     2    nnz,row_ptr,col_ind,iquad,ndquad,nquad,wnear,
     3    eps_gmres,niter,errs,rres,soln,u)

c
c  This subroutine solves the Laplace Beltrami problem
c  using one of the following three integral representations:
c
c  irep=1, s_lap_s: S (\Delta_{\Gamma} + W)S[\sigma] = S[f], where
c                      S \Delta_{\Gamma}S [\sigma] is expanded out using 
c                      Calderon identities
c  irep=2, s_lap_s_noc: Same as above but no Calderon identities
c            used
c  irep=3, lap_s2: \Delta_{\Gamma} + W S^2 [\sigma], where
c                    \Delta_{\Gamma}S^2 is expanded out
c                    using Calderon identities
c
c  The linear system is solved iteratively using GMRES.
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
c    - numit: integer
c        maximum number of GMRES iterations
c    - rhs: real *8 (npts)
c        the surface data f
c    - irep: integer
c        representation for solving the Laplace Beltrami problem
c        * irep = 1, s_lap_s
c        * irep = 2, s_lap_s_noc
c        * irep = 3, lap_s2
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
c    - ndquad: integer
c        trailing dimension of quadrature correction array
c    - nquad: integer
c        number of near field entries corresponding to each source target
c        pair. The size of wnear is (nquad,4) since there are 4 kernels
c        per source target pair
c    - wnear: real *8(nquad,ndquad)
c        The desired near field quadrature
c    - eps_gmres: real *8
c        relative residual tolerance for GMRES
c
c  Output arguments:
c    - niter: integer
c        number of GMRES iterations used
c    - errs: real *8 (niter+1)
c        On input must be of size (numit+1), gmres residual
c        as a function of iteration number
c    - rres: real *8
c        relative residual to which the linear system is solved
c    - soln: real *8 (npts)
c        the solution sigma
c    - u: real *8 (npts)
c        The solution to the Laplace Beltrami problem
c
      implicit none

      integer, intent(in) :: npatches,npts
      integer, intent(in) :: norders(npatches),ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      integer, intent(in) :: numit,irep
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      real *8, intent(in) :: eps_gmres
      real *8, intent(in) :: rhs(npts)

      integer, intent(in) :: nnz,ndquad,nquad
      integer, intent(in) :: row_ptr(npts+1),col_ind(nnz),iquad(nnz+1)
      real *8, intent(in) :: wnear(nquad,ndquad)

      real *8, intent(out) :: soln(npts),u(npts)
      real *8, intent(out) :: errs(numit+1)
      real *8, intent(out) :: rres
      integer, intent(out) ::  niter
      

      real *8, allocatable :: rhs_use(:)
      real *8, allocatable :: targs(:,:)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      integer ndtarg,ntarg



      integer norder,npols
      integer nover,npolso,nptso
      real *8, allocatable :: s_one(:),dens_one(:)

      real *8, allocatable :: srcover(:,:),wover(:)
      integer, allocatable :: ixyzso(:),novers(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:) 

      integer i,j,jpatch,jquadstart,jstart

      integer ipars
      complex *16 zpars
      real *8 timeinfo(10),t1,t2,omp_get_wtime


      real *8 ttot,done,pi
      real *8 rfac,rfac0
      real *8 dpars(2),erra,ra
      integer iptype_avg,norder_avg
      integer ikerorder, iquadtype,npts_over

c
c
c       gmres variables
c
      real *8 did,dtmp
      real *8 rb,wnrm2
      integer it,iind,it1,k,l
      real *8 rmyerr
      real *8 temp
      real *8, allocatable :: vmat(:,:),hmat(:,:)
      real *8, allocatable :: cs(:),sn(:)
      real *8, allocatable :: svec(:),yvec(:),wtmp(:)
      real *8, allocatable :: wts(:)

      complex *16 ztmp


      allocate(vmat(npts,numit+1),hmat(numit,numit))
      allocate(cs(numit),sn(numit))
      allocate(wtmp(npts),svec(numit+1),yvec(numit+1))


      done = 1
      pi = atan(done)*4

      if(irep.lt.1.or.irep.gt.3) then
        print *, "invalid argument for irep, returning"
        return
      endif


c
c
c        setup targets as on surface discretization points
c 
      ndtarg = 3
      ntarg = npts
      allocate(targs(ndtarg,npts),uvs_targ(2,ntarg),ipatch_id(ntarg))

      allocate(rhs_use(npts))

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
c    estimate oversampling for far-field, and oversample geometry
c

      allocate(novers(npatches),ixyzso(npatches+1))

      print *, "beginning far order estimation"

      ztmp = 0
      ikerorder = 0
      print *, eps

      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,
     1    rads,npts,srccoefs,ndtarg,npts,targs,ikerorder,ztmp,
     2    nnz,row_ptr,col_ind,rfac,novers,ixyzso)

      npts_over = ixyzso(npatches+1)-1
      print *, "npts_over=",npts_over

      allocate(srcover(12,npts_over),wover(npts_over))

      call oversample_geom(npatches,norders,ixyzs,iptype,npts, 
     1   srccoefs,srcvals,novers,ixyzso,npts_over,srcover)

      allocate(wts(npts))

      call get_qwts(npatches,norders,ixyzs,iptype,npts,
     1        srcvals,wts)

      call get_qwts(npatches,novers,ixyzso,iptype,npts_over,
     1        srcover,wover)


      print *, "Starting to generate right hand side for linear system"

c
c
c     compute the right hand side S[f], if irep=1, or irep=2
c     note that
c
      if(irep.eq.1.or.irep.eq.2) then
        dpars(1) = 1.0d0
        dpars(2) = 0
        call lpcomp_lap_comb_dir_addsub(npatches,norders,ixyzs,
     1    iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,eps,
     2    dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear(1,1),
     3    rhs,novers,npts_over,ixyzso,srcover,wover,rhs_use)
      else
        print  *, "here1"
C$OMP PARALLEL DO DEFAULT(SHARED)
        do i=1,npts
          rhs_use(i) = rhs(i)
        enddo
C$OMP END PARALLEL DO
      endif
      print *, "Done generating right hand side for linear system"

      if(irep.eq.2) then
        allocate(s_one(npts),dens_one(npts))

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
        do i=1,npts
          dens_one(i) = 1.0d0
        enddo
C$OMP END PARALLEL DO

        dpars(1) = 1.0d0
        dpars(2) = 0
        call lpcomp_lap_comb_dir_addsub(npatches,norders,ixyzs,
     1    iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,eps,
     2    dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear(1,1),
     3    dens_one,novers,npts_over,ixyzso,srcover,wover,s_one)
        

      endif

      print *, "Start gmres now"
c
c
c     start gmres code here
c
c     NOTE: matrix equation should be of the form (z*I + K)x = y
c       the identity scaling (z) is defined via did below,
c       and K represents the action of the principal value 
c       part of the matvec
c

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,npts
        u(i) = 0
      enddo
C$OMP END PARALLEL DO


      if(irep.eq.1.or.irep.eq.3) then
        did = -0.25d0
      else
        did = 0
      endif


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
        rb = rb + abs(rhs_use(i))**2*wts(i)
      enddo
      rb = sqrt(rb)

      do i=1,npts
        vmat(i,1) = rhs_use(i)/rb
      enddo

      svec(1) = rb

      do it=1,numit
        it1 = it + 1

c
c        NOTE:
c        replace this routine by appropriate layer potential
c        evaluation routine  
c
        if(irep.eq.1) then
          call lpcomp_lap_bel_s_lap_s_addsub(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,eps,nnz,row_ptr,col_ind,iquad,
     2     nquad,wnear,vmat(1,it),novers,npts_over,ixyzso,
     3     srcover,wover,wtmp,u)
        endif
        if(irep.eq.2) then
          call lpcomp_lap_bel_s_lap_s_noc_addsub(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,eps,nnz,row_ptr,col_ind,iquad,
     2     nquad,wnear,vmat(1,it),novers,npts_over,ixyzso,
     3     srcover,wover,s_one,wtmp,u)
        endif
        if(irep.eq.3) then
          call lpcomp_lap_bel_lap_s2_addsub(npatches,norders,ixyzs,
     1     iptype,npts,srccoefs,srcvals,eps,nnz,row_ptr,col_ind,iquad,
     2     nquad,wnear,vmat(1,it),novers,npts_over,ixyzso,
     3     srcover,wover,wtmp,u)
        endif


        do k=1,it
          hmat(k,it) = 0
          do j=1,npts
            hmat(k,it) = hmat(k,it) + wtmp(j)*vmat(j,k)*wts(j)
          enddo

          do j=1,npts
            wtmp(j) = wtmp(j)-hmat(k,it)*vmat(j,k)
          enddo
        enddo
          
        hmat(it,it) = hmat(it,it)+did
        wnrm2 = 0
        do j=1,npts
          wnrm2 = wnrm2 + abs(wtmp(j))**2*wts(j)
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
          if(irep.eq.1) then
            call lpcomp_lap_bel_s_lap_s_addsub(npatches,norders,ixyzs,
     1       iptype,npts,srccoefs,srcvals,eps,nnz,row_ptr,col_ind,iquad,
     2       nquad,wnear,soln,novers,npts_over,ixyzso,
     3       srcover,wover,wtmp,u)
          endif
          if(irep.eq.2) then
            call lpcomp_lap_bel_s_lap_s_noc_addsub(npatches,norders,
     1       ixyzs,iptype,npts,srccoefs,srcvals,eps,nnz,row_ptr,col_ind,
     2       iquad,nquad,wnear,soln,novers,npts_over,ixyzso,
     3       srcover,wover,s_one,wtmp,u)
          endif

          if(irep.eq.3) then
            call lpcomp_lap_bel_lap_s2_addsub(npatches,norders,ixyzs,
     1       iptype,npts,srccoefs,srcvals,eps,nnz,row_ptr,col_ind,iquad,
     2       nquad,wnear,soln,novers,npts_over,ixyzso,
     3       srcover,wover,wtmp,u)
          endif

          do i=1,npts
            rres = rres + abs(did*soln(i) + wtmp(i)-rhs(i))**2*wts(i)
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
c
c

      subroutine get_hodge_decomposition(npatches,norders,ixyzs,
     1    iptype,npts,srccoefs,srcvals,eps,numit,rhs,irep,
     2    eps_gmres,niter,errs,rres,cfree,dfree,vharm)
c
c  This subroutine computes the hodge decomposition of a given
c  vector field v = rhs.
c 
c  First, we compute its tangential projection given by
c      v -> v - (v.n)n 
c  
c  Then the hodge decomposition of the tangential vector field
c  is given by
c      v = \nabla_{\Gamma} \alpha + n\times \nabla_{\Gamma} \beta +
c             h
c  where \nabla_{\Gamma} \alpha is the curl free component (cfree)
c  n \times \nabla_{\Gamma} \beta is the divergence free component (dfree)
c  and h is the harmonic component (vharm)
c
c  The components are computed by solving two Laplace beltrami problems
c  \Delta_{\Gamma} \alpha = \nabla_{\Gamma} \cdot v, and 
c  \Delta_{\Gamma} \beta = -\nabla_{\Gamma} \cdot n \times v.
c
c  Finally h is given by
c   h = v - \nabla_{\Gamma} \alpha - n \times \nabla_{\Gamma} \beta
c  
c  The Laplace Beltrami problems are solved
c  using one of the following three integral representations:
c
c  irep=1, s_lap_s: S (\Delta_{\Gamma} + W)S[\sigma] = S[f], where
c                      S \Delta_{\Gamma}S [\sigma] is expanded out using 
c                      Calderon identities
c  irep=2, s_lap_s_noc: Same as above but no Calderon identities
c            used
c  irep=3, lap_s2: \Delta_{\Gamma} + W S^2 [\sigma], where
c                    \Delta_{\Gamma}S^2 is expanded out
c                    using Calderon identities
c
c  The linear system is solved iteratively using GMRES.
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
c    - numit: integer
c        maximum number of GMRES iterations
c    - rhs: real *8 (3,npts)
c        the input vector field
c    - irep: integer
c        representation for solving the Laplace Beltrami problem
c        * irep = 1, s_lap_s
c        * irep = 2, s_lap_s_noc
c        * irep = 3, lap_s2
c    - eps_gmres: real *8
c        relative residual tolerance for GMRES
c
c  Output arguments:
c    - niter: integer(2)
c        number of GMRES iterations used for the two solves
c    - errs: real *8 (niterj+1,2)
c        On input must be of size (numit+1,2), gmres residual
c        as a function of iteration number
c    - rres: real *8 (2)
c        relative residual to which the linear system is solved
c    - cfree: real *8 (3,npts)
c        the curl free component (\nabla_{\Gamma} \alpha)
c    - dfree: real *8 (3,npts)
c        the divergence free component (n \times \nabla_{\Gamma} \beta)
c    - vharm: real *8 (3,npts)
c        the harmonic component
c
      implicit none

      integer, intent(in) :: npatches,npts
      integer, intent(in) :: norders(npatches),ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches)
      integer, intent(in) :: numit,irep
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      real *8, intent(in) :: eps_gmres
      real *8, intent(in) :: rhs(3,npts)
      real *8, intent(out) :: cfree(3,npts),dfree(3,npts),vharm(3,npts)
      real *8, intent(out) :: errs(numit+1,2)
      real *8, intent(out) :: rres(2)
      integer, intent(out) ::  niter(2)
      

      real *8, allocatable :: rhs_use(:),rhs_tan(:,:),rhs_tan2(:,:)
      real *8, allocatable :: targs(:,:)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      integer ndtarg,ntarg



      integer norder,npols
      integer nover,npolso,nptso
      integer nnz,nquad
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      real *8, allocatable :: wnear(:,:)
      real *8, allocatable :: s_one(:),dens_one(:)

      real *8, allocatable :: srcover(:,:),wover(:)
      integer, allocatable :: ixyzso(:),novers(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)

      real *8, allocatable :: alpha(:),beta(:),sgbeta(:,:)
      real *8, allocatable :: alpha_sig(:),beta_sig(:)

      integer i,j,jpatch,jquadstart,jstart

      integer ipars
      complex *16 zpars
      real *8 timeinfo(10),t1,t2,omp_get_wtime


      real *8 ttot,done,pi
      real *8 rfac,rfac0,rn
      real *8 dpars(2),erra,ra
      integer iptype_avg,norder_avg
      integer ikerorder, iquadtype,npts_over

c
c
c       gmres variables
c
      real *8 did,dtmp
      real *8 rb,wnrm2
      integer it,iind,it1,k,l,ndquad
      real *8 rmyerr
      real *8 temp
      real *8, allocatable :: wts(:)

      complex *16 ztmp


      done = 1
      pi = atan(done)*4

      if(irep.lt.1.or.irep.gt.3) then
        print *, "invalid argument for irep, returning"
        return
      endif


c
c
c        setup targets as on surface discretization points
c 
      ndtarg = 3
      ntarg = npts
      allocate(targs(ndtarg,npts),uvs_targ(2,ntarg),ipatch_id(ntarg))

      allocate(rhs_tan(3,npts),rhs_use(npts),rhs_tan2(3,npts))

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(rn)
      do i=1,ntarg
        targs(1,i) = srcvals(1,i)
        targs(2,i) = srcvals(2,i)
        targs(3,i) = srcvals(3,i)
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
        rn = rhs(1,i)*srcvals(10,i) + rhs(2,i)*srcvals(11,i) + 
     1     rhs(3,i)*srcvals(12,i)
        rhs_tan(1,i) = rhs(1,i) - rn*srcvals(10,i)
        rhs_tan(2,i) = rhs(2,i) - rn*srcvals(11,i)
        rhs_tan(3,i) = rhs(3,i) - rn*srcvals(12,i)
        call cross_prod3d(srcvals(10,i),rhs_tan(1,i),rhs_tan2(1,i))
      enddo
C$OMP END PARALLEL DO   


c
c    initialize patch_id and uv_targ for on surface targets
c
      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,
     1  ipatch_id,uvs_targ)

c
c
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
      call findnearmem(cms,npatches,rad_near,3,targs,npts,nnz)
      print *, "nnz=",nnz

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,3,targs,npts,row_ptr, 
     1        col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind,
     1         iquad)


      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,
     1        srcvals,wts)

c
c   compute near quadrature correction
c
      nquad = iquad(nnz+1)-1
      print *, "nquad=",nquad
      if(irep.eq.1.or.irep.eq.2) then
        allocate(wnear(nquad,4))
        ndquad = 4
      else
        allocate(wnear(nquad,3))
        ndquad = 3
      endif
      
      do j=1,ndquad
C$OMP PARALLEL DO DEFAULT(SHARED)      
        do i=1,nquad
          wnear(i,j) = 0
        enddo
C$OMP END PARALLEL DO    
      enddo


      iquadtype = 1

      print *, "starting to generate near quadrature"
      call cpu_time(t1)
C$      t1 = omp_get_wtime()      
      if(irep.eq.1) then
        call getnearquad_lap_bel_s_lap_s(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,eps,
     2   iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear)
      else if(irep.eq.2) then
        call getnearquad_lap_bel_s_lap_s_noc(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,eps,
     2   iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear)
      else
        call getnearquad_lap_bel_lap_s2(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,eps,
     2   iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear)
      endif
      call cpu_time(t2)
C$      t2 = omp_get_wtime()     

      call prin2('quadrature generation time=*',t2-t1,1)
      print *, "done generating near quadrature"
      

      allocate(alpha(npts),alpha_sig(npts))
      allocate(beta(npts),beta_sig(npts))

      call get_surf_div_cartesian(1,npatches,norders,ixyzs,iptype,
     1  npts,srccoefs,srcvals,rhs_tan,rhs_use)
      

      call lap_bel_solver_wquad(npatches,norders,ixyzs,
     1    iptype,npts,srccoefs,srcvals,eps,numit,rhs_use,irep,
     2    nnz,row_ptr,col_ind,iquad,ndquad,nquad,wnear,
     3    eps_gmres,niter(1),errs(1,1),rres(1),alpha_sig,alpha)
      
      call get_surf_grad_cartesian(1,npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,alpha,cfree)


      call get_surf_div_cartesian(1,npatches,norders,ixyzs,iptype,
     1  npts,srccoefs,srcvals,rhs_tan2,rhs_use)
      
c
c
c
C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
         rhs_use(i) = -rhs_use(i)
      enddo
C$OMP END PARALLEL DO      
      

      call lap_bel_solver_wquad(npatches,norders,ixyzs,
     1    iptype,npts,srccoefs,srcvals,eps,numit,rhs_use,irep,
     2    nnz,row_ptr,col_ind,iquad,ndquad,nquad,wnear,
     3    eps_gmres,niter(2),errs(1,2),rres(2),beta_sig,beta)
      allocate(sgbeta(3,npts))
      
      call get_surf_grad_cartesian(1,npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,beta,sgbeta)
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)      
      do i=1,npts
        call cross_prod3d(srcvals(10,i),sgbeta(1,i),dfree(1,i))
        vharm(1,i) = rhs_tan(1,i) - cfree(1,i) - dfree(1,i)
        vharm(2,i) = rhs_tan(2,i) - cfree(2,i) - dfree(2,i)
        vharm(3,i) = rhs_tan(3,i) - cfree(3,i) - dfree(3,i)
      enddo
C$OMP END PARALLEL DO      
c


c
c
c
c
      return
      end
c
c
c
c
c
