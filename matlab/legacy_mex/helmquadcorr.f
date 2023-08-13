

      subroutine getnearquad_helm_comb_dir_spmat(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2   ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,
     3   iquad,rfac0,nquad,wnear,irowind,icolind)
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
c            irowind - integer(nquad)
c              row indices corresponding to near quadrature
c            icolind - integer(nquad)
c              column indices corresponding to near quadrature
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
      integer, intent(out) :: irowind(nquad), icolind(nquad)
      complex *16, intent(out) :: wnear(nquad)


      integer ipars
      integer ndd,ndz,ndi
      real *8 dpars

      complex *16 alpha,beta
      integer i,j,jpatch,jstart,jquadstart,l,npols
      integer ipv
      real *8 t1,t2
      
      call cpu_time(t1)
      call getnearquad_helm_comb_dir(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2   ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,
     3   iquad,rfac0,nquad,wnear)
      call cpu_time(t2)

ccC$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,npols,jquadstart)
ccC$OMP$PRIVATE(jstart,l)
      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch)
          do l=1,npols
            irowind(jquadstart+l-1) = i
            icolind(jquadstart+l-1) = jstart+l-1
          enddo
        enddo
      enddo
ccC$OMP END PARALLEL DO      

      return
      end



      subroutine getnearquad_helm_rpcomb_neu_spmat(npatches,norders, 
     1  ixyzs,iptype,npts,srccoefs,srcvals, 
     1  eps,zpars,iquadtype,nnz,row_ptr,col_ind,
     1  iquad,rfac0,nquad,wnear,irowind,icolind)
c
c  This subroutine generates the near field quadrature
c  for the representation:
c
c  u = S_{k}[\rho]+i*alpha*D_{k}[S_{i|k|}[\rho]]
c
c  and returns quantities related to evaluating du/dn on surface
c  at the surface discretization nodes
c
c  If values at other nodes is desired then the solution
c  should be reinterpolated to those nodes
c
c
c  On imposing the boundary condition, we get the following operator
c
c  du/dn = -I/2 + S_{k}' + i \alpha (D_{k}'-D_{i|k|}') S_{i|k|}
c    - i \alpha I/4 + i \alpha (S_{i|k|}')^2
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
c    - wnear must be of size 4*nquad as 4 different layer
c      potentials are returned
c      * the first kernel is S_{k}'
c      * the second kernel is S_{i|k|}
c      * the third kernel is S_{i|k|}'
c      * the fourth kernel is D_{k}'-D_{i|k|}'
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
c    - zpars: complex *16 (2)
c        kernel parameters (Referring to formula (1))
c        zpars(1) = k 
c        zpars(2) = alpha
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
c        pair. The size of wnear is 4*nquad since there are 4 kernels
c        per source target pair
c
c  Output arguments
c    - wnear: complex *16(4*nquad)
c        The desired near field quadrature
c    - irowind: integer(nquad)
c        row indices corresponding to near quadrature
c    - icolind: integer(nquad)
c        column indices corresponding to near quadrature
c               
c

      implicit none 
      integer npatches,norders(npatches),npts,nquad
      integer ixyzs(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,rfac0
      integer ndtarg,ntarg
      integer iquadtype
      complex *16 zpars(2)
      complex *16 zpars_tmp(3)
      integer nnz,ipars(2)
      real *8 dpars(1)
      integer row_ptr(npts+1),col_ind(nnz),iquad(nnz+1)
      complex *16 wnear(nquad,4)
      integer irowind(nquad),icolind(nquad)

      real *8, allocatable :: uvs_targ(:,:)
      integer, allocatable :: ipatch_id(:)


      complex *16 alpha,beta,ima,zk
      integer i,j,ndi,ndd,ndz

      integer ipv
      integer jpatch,jquadstart,jstart,l,npols

      procedure (), pointer :: fker
      external h3d_sprime,h3d_slp,h3d_dprime_diff

      data ima/(0.0d0,1.0d0)/

      ndz=1
      ndd=1
      ndi=2
      ndtarg = 12
      ntarg = npts
      zk = zpars(1)

      allocate(ipatch_id(npts),uvs_targ(2,npts))
cc$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo
cc$OMP END PARALLEL DO      

      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,
     1   ipatch_id,uvs_targ)
      ipv=1
      fker => h3d_sprime 
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1   iptype,npts,srccoefs,srcvals,ndtarg,ntarg,srcvals,
     1   ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,
     1   ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(1,1))

      zpars_tmp(1) = ima*abs(zk)
      fker => h3d_slp
      ipv = 1
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs, 
     1   iptype,npts,srccoefs,srcvals,ndtarg,ntarg,srcvals, 
     1   ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars_tmp,
     1   ndi,ipars,nnz,row_ptr,col_ind,iquad, 
     1   rfac0,nquad,wnear(1,2))
      
      fker => h3d_sprime
      ipv = 1
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,ndtarg,ntarg,srcvals,
     1  ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars_tmp,ndi,
     1  ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(1,3))

      zpars_tmp(1) = zk
      zpars_tmp(2) = ima*abs(zk)
      ndz = 2
      fker => h3d_dprime_diff
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,ndtarg,ntarg,srcvals,
     1  ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars_tmp,ndi,
     1  ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(1,4))

ccC$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,npols,jquadstart)
ccC$OMP$PRIVATE(jstart,l)
      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch)
          do l=1,npols
            irowind(jquadstart+l-1) = i
            icolind(jquadstart+l-1) = jstart+l-1
          enddo
        enddo
      enddo
ccC$OMP END PARALLEL DO      


      return
      end subroutine getnearquad_helm_rpcomb_neu_spmat




      subroutine getnearquad_helm_comb_trans_spmat(npatches,norders,
     1  ixyzs,iptype,npts,srccoefs,srcvals, 
     1  eps,zpars,iquadtype,nnz,row_ptr,col_ind,
     1  iquad,rfac0,nquad,wnear,irowind,icolind)
c
c  This subroutine generates the near field quadrature
c  for the representations:
c
c  u1 = ep1^2 S_{k1}[\lambda]+ep1*D_{k1}[\rho] (interior representation)
c  u0 = ep0^2 S_{k0}[\lambda]+ep0*D_{k0}[\rho] (exterior representation)
c
c  and returns quantities related to u0-u1 , and 1/ep0 u0' - 1/ep1 u1'
c  on surface
c
c  On imposing the boundary condition, we get the following 
c  sets of operators
c
c  u0-u1 = (ep0+ep1)/2 \rho + (ep0 D_{k0} - \ep1 D_{k1})[\rho] +
c    (ep0^2 S_{k1} - ep1^2 S_{k0})[\lambda]
c  
c  1/ep0 u0' - 1/ep1 u1' = -(ep0 + ep1)/2 \lambda + 
c    (D_{k0}'-D_{k1}')[\rho] + (ep0 S_{k0}' - ep1 S_{k1}')[\lambda]
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
c    - wnear must be of size 4*nquad as 4 different layer
c      potentials are returned
c      * the first kernel is ep0^2 S_{k0} - ep1^2 S_{k1}
c      * the second kernel is ep0 D_{k0} - ep1 D_{k1}
c      * the third kernel is ep0 S_{k0}' - ep1 S_{k1}'
c      * the fourth kernel is D_{k0}'-D_{k1}'
c    - The parameters that must be provided zpars(5), and how they
c      are related to the problem parameters
c        zpars(1) = \omega
c        zpars(2) = ep0
c        zpars(3) = mu0
c        zpars(4) = ep1
c        zpars(5) = mu1
c        
c        k0  = \omega*sqrt(ep0*mu0)
c        k1 = \omega*sqrt(ep1*mu1)
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
c    - zpars: complex *16 (5)
c        kernel parameters (See notes above)
c        zpars(1) = \omega 
c        zpars(2) = ep0
c        zpars(3) = mu0
c        zpars(4) = ep1
c        zpars(5) = mu1
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
c        pair. The size of wnear is 4*nquad since there are 4 kernels
c        per source target pair
c
c  Output arguments
c    - wnear: complex *16(nquad,4)
c        The desired near field quadrature
c        * First kernel is difference of S
c        * Second kernel is difference of D
c        * third kernel is difference of S'
c        * fourth kernel is difference of D'
c                
c    - irowind: integer(nquad)
c        row indices corresponding to near quadrature
c    - icolind: integer(nquad)
c        column indices corresponding to near quadrature
c               
c

      implicit none 
      integer npatches,norders(npatches),npts,nquad
      integer ixyzs(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,rfac0
      integer, intent(out) :: irowind(nquad), icolind(nquad)
      integer ndtarg,ntarg
      integer iquadtype
      complex *16 zpars(5)
      complex *16 zpars_tmp(5)
      complex *16 zk0,zk1
      integer nnz,ipars(2)
      real *8 dpars(1)
      integer row_ptr(npts+1),col_ind(nnz),iquad(nnz+1)
      complex *16 wnear(nquad,4)

      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      integer jpatch,jquadstart,l,npols,jstart


      complex *16 alpha,beta,ima,zk
      integer i,j,ndi,ndd,ndz

      integer ipv

      procedure (), pointer :: fker
      external h3d_sprime_diff,h3d_slp_diff,h3d_dlp_diff
      external h3d_dprime_diff

      data ima/(0.0d0,1.0d0)/

      ndz=4
      ndd=0
      ndi=0
      ndtarg = 12
      ntarg = npts
      zk0 = zpars(1)*sqrt(zpars(2)*zpars(3))
      zk1 = zpars(1)*sqrt(zpars(4)*zpars(5))
      zpars_tmp(1) = zk0
      zpars_tmp(2) = zk1

      allocate(ipatch_id(npts),uvs_targ(2,npts))
c$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo
c$OMP END PARALLEL DO      

      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,
     1   ipatch_id,uvs_targ)
      ipv=0
      fker => h3d_slp_diff
      zpars_tmp(3) = zpars(2)**2
      zpars_tmp(4) = zpars(4)**2
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1   iptype,npts,srccoefs,srcvals,ndtarg,ntarg,srcvals,
     2   ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars_tmp,
     2   ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(1,1))

      fker => h3d_dlp_diff
      zpars_tmp(3) = zpars(2)
      zpars_tmp(4) = zpars(4)
      ipv = 1
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs, 
     1   iptype,npts,srccoefs,srcvals,ndtarg,ntarg,srcvals, 
     2   ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars_tmp, 
     2   ndi,ipars,nnz,row_ptr,col_ind,iquad, 
     2   rfac0,nquad,wnear(1,2))
      
      fker => h3d_sprime_diff
      ipv = 1
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,ndtarg,ntarg,srcvals,
     2  ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars_tmp,ndi,
     2  ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(1,3))

      ndz = 2
      fker => h3d_dprime_diff
      ipv = 0
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,ndtarg,ntarg,srcvals,
     2  ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars_tmp,ndi,
     2  ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(1,4))


ccC$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,npols,jquadstart)
ccC$OMP$PRIVATE(jstart,l)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch)
          do l=1,npols
            irowind(jquadstart+l-1) = i
            icolind(jquadstart+l-1) = jstart+l-1
          enddo
        enddo
      enddo
ccC$OMP END PARALLEL DO      


      return
      end subroutine getnearquad_helm_comb_trans_spmat
c
c
c
c
