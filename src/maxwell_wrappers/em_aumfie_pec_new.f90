!
!
! This file contains the following suer callable routines
! for solving Maxwell's equations with perfect
! electric conductor boundary conditions
! using the magnetic field integral equation (mfie) 
!
! PDE: 
!   curl E = ik H;   curl H = - ik E
!
! Representation:
!   H = curl S_{k}[J] 
!   E =  ik S_{k} [J] - \nabla S_{k} [\rho]
!
! Boundary condtions imposed:
!   n \times H = J
!   n \dot E = \rho
!
! Boundary integral equation
!   J/2 - n \times curl S_{k} [J] = n \times H_{inc}
!   \rho/2 + S'_{k}[\rho] = ik S_{k}[J] + n \cdot E_{inc}
! 
!
! em_aumfie_solver: an FMM accelerated GMRES-solver 
! for finding the current density of the magnetic field integral
! equation
!
! lpcomp_em_aumfie_EH: a layer potential evaluator for
! computing the electric and magnetic fields corresponding to the
! augmented magnetic field integral equation, using both charge and current
! density
!
!
! The file also contains the following advanced routines
!
! getnearquad_em_mfie_pec: generate near quadrature correction for 
!   mfie
!
      subroutine getnearquad_em_aumfie_pec(npatches,norders,ixyzs, &
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,ipatch_id, &
        uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,iquad,rfac0, &
        nquad,wnear)
!
!  This subroutine generates the near field quadrature
!  for the representation:
!
!  H = curl S_{k}[J]
!  E = ik S_{k} [J] - \nabla S_{k} [\rho]
!
!
!  The quadrature is computed by the following strategy
!  targets within a sphere of radius rfac0*rs
!  of a patch centroid is handled using adaptive integration
!  where rs is the radius of the bounding sphere
!  for the patch
!  
!  All other targets in the near field are handled via
!  oversampled quadrature
!
!  The recommended parameter for rfac0 is 1.25d0
!  
!  NOTES:
!    - wnear must be of size nquad,4 as 4 different layer
!      potentials are returned
!      * the first kernel is S_{k}
!      * the second kernel is \partial_{x} S_{k}
!      * the third kernel is \partial_{y} S_{k}
!      * the fourth kernel is \partial_{z} S_{k}
! 
!  Input arguments:
! 
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization on each patch 
!    - ixyzs: integer(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!    - npts: integer
!        total number of discretization points on the boundary
!    - srccoefs: real *8 (9,npts)
!        koornwinder expansion coefficients of xyz, dxyz/du,
!        and dxyz/dv on each patch. 
!        For each point 
!          * srccoefs(1:3,i) is xyz info
!          * srccoefs(4:6,i) is dxyz/du info
!          * srccoefs(7:9,i) is dxyz/dv info
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the 
!        discretization nodes on the surface
!          * srcvals(1:3,i) - xyz info
!          * srcvals(4:6,i) - dxyz/du info
!          * srcvals(7:9,i) - dxyz/dv info
!          * srcvals(10:12,i) - normals info
!    - ndtarg: integer
!        leading dimension of target array, 
!        must be at least 3
!    - ntarg: integer
!        number of targets
!    - targs: real *8 (ndtarg,ntarg)
!         target info, the first three components targs(1:3,i) must
!         be the xyz location of the target
!    - ipatch_id: integer(ntarg)
!        id of patch if target i is on surface, id=-1 if target if off-surface
!    - uvs_targ: real *8 (2,ntarg)
!        local uv coordinates of patch of target, if it is on surface, 
!        otherwise unused
!    - eps: real *8
!        precision requested
!    - zpars: complex *16 (1)
!        kernel parameters 
!        zpars(1) = k 
!    - iquadtype: integer
!        quadrature type
!          * iquadtype = 1, use ggq for self + adaptive integration
!            for rest
!    - nnz: integer
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer(ntarg+1)
!        row_ptr(i) is the pointer
!        to col_ind array where list of relevant source patches
!        for target i start
!    - col_ind: integer (nnz)
!        list of source patches relevant for all targets, sorted
!        by the target number
!    - iquad: integer(nnz+1)
!        location in wnear_ij array where quadrature for col_ind(i)
!        starts for a single kernel. In this case the different kernels
!        are matrix entries are located at (m-1)*nquad+iquad(i), where
!        m is the kernel number
!    - rfac0: integer
!        radius parameter for near field
!    - nquad: integer
!        number of near field entries corresponding to each source target
!        pair. The size of wnear is 4*nquad since there are 4 kernels
!        per source target pair
!
!  Output arguments
!    - wnear: complex *16(nquad,4)
!        The desired near field quadrature
!               
!
      implicit none
      integer, intent(in) :: npatches,norders(npatches)
      integer, intent(in) :: ixyzs(npatches+1),iptype(npatches),npts
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
      integer, intent(in) :: ndtarg,ntarg,ipatch_id(ntarg)
      real *8, intent(in) :: targs(ndtarg,ntarg),uvs_targ(2,ntarg),eps
      complex *16, intent(in) :: zpars(1)
      integer, intent(in) :: iquadtype,nnz,row_ptr(ntarg+1),col_ind(nnz)
      integer, intent(in) :: iquad(nnz+1),nquad
      real *8, intent(in) :: rfac0
      complex *16, intent(out) :: wnear(nquad,4)

      integer ndz,ndd,ndi
      integer ipv
      real *8 dpars
      integer ipars
      procedure (), pointer :: fker
      external h3d_slp,h3d_sgradx,h3d_sgrady,h3d_sgradz

      ndd = 0
      ndz = 1
      ndi = 0
      dpars = 0
      ipars = 0

      ipv = 0
      fker => h3d_slp 
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,iptype,&
       npts,srccoefs,srcvals,ndtarg,ntarg,targs,ipatch_id, &
       uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz, &
       row_ptr,col_ind,iquad,rfac0,wnear(1,1))

      ipv = 1
      fker => h3d_sgradx
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,iptype,&
       npts,srccoefs,srcvals,ndtarg,ntarg,targs,ipatch_id, &
       uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz, &
       row_ptr,col_ind,iquad,rfac0,wnear(1,2))

      fker => h3d_sgrady
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,iptype,&
       npts,srccoefs,srcvals,ndtarg,ntarg,targs,ipatch_id, &
       uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz, &
       row_ptr,col_ind,iquad,rfac0,wnear(1,3))

      fker => h3d_sgradz
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,iptype,&
       npts,srccoefs,srcvals,ndtarg,ntarg,targs,ipatch_id, &
       uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz, &
       row_ptr,col_ind,iquad,rfac0,wnear(1,4))
      
      return
      end subroutine getnearquad_em_aumfie_pec
!
!
!
!
!
!
      subroutine lpcomp_em_aumfie_EH_addsub(npatches,norders,ixyzs, &
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,eps,zpars,nnz, &
        row_ptr,col_ind,iquad,nquad,wnear,zjvec,zrhovec,novers,nptso, &
        ixyzso,srcover,whtsover,E,H)
!
!
!
!  This subroutine evaluates the magnetic field and the electric 
!  field corresponding to the following integral representation
!
!  H = curl S_{k}[J]
!  E = ik S_{k} [J] - \nabla S_{k} [\rho]
!
!  where the near field is precomputed and stored in the
!  row sparse compressed format.
!
!  The identity term for on surface targets
!  is not included in the representation
!  
!  The fmm is used to accelerate the far-field and 
!  the near-field interactions are corrected using the precomputed
!  quadrature
!
!  Using add and subtract - no need to call tree and set fmm
!  parameters. Can directly call existing fmm library.
!  
!  NOTES:
!    - the routine assumes that the current is tangential
!    - the routine also on input assumes that the current 
!    is provided in the orthonormal frame with unit vectors
!    in the direction given by 
!    srcvals(4:6,i), srcvals(4:6,i) \times srcvals(10:12,i)   
!
! 
!  Input arguments:
! 
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization on each patch 
!    - ixyzs: integer(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!    - npts: integer
!        total number of discretization points on the boundary
!    - srccoefs: real *8 (9,npts)
!        koornwinder expansion coefficients of xyz, dxyz/du,
!        and dxyz/dv on each patch. 
!        For each point 
!          * srccoefs(1:3,i) is xyz info
!          * srccoefs(4:6,i) is dxyz/du info
!          * srccoefs(7:9,i) is dxyz/dv info
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the 
!        discretization nodes on the surface
!          * srcvals(1:3,i) - xyz info
!          * srcvals(4:6,i) - dxyz/du info
!          * srcvals(7:9,i) - dxyz/dv info
!          * srcvals(10:12,i) - normals info
!    - ndtarg: integer
!        leading dimension of target array, 
!        must be at least 3
!    - ntarg: integer
!        number of targets
!    - targs: real *8 (ndtarg,ntarg)
!         target info, the first three components targs(1:3,i) must
!         be the xyz location of the target
!    - eps: real *8
!        precision requested
!    - zpars: complex *16 (1)
!        kernel parameters 
!        zpars(1) = k 
!    - iquadtype: integer
!        quadrature type
!          * iquadtype = 1, use ggq for self + adaptive integration
!            for rest
!    - nnz: integer
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer(ntarg+1)
!        row_ptr(i) is the pointer
!        to col_ind array where list of relevant source patches
!        for target i start
!    - col_ind: integer (nnz)
!        list of source patches relevant for all targets, sorted
!        by the target number
!    - iquad: integer(nnz+1)
!        location in wnear_ij array where quadrature for col_ind(i)
!        starts for a single kernel. In this case the different kernels
!        are matrix entries are located at (m-1)*nquad+iquad(i), where
!        m is the kernel number
!    - nquad: integer
!        number of near field entries corresponding to each source target
!        pair. The size of wnear is 4*nquad since there are 4 kernels
!        per source target pair
!    - wnear: complex *16(nquad,4)
!        Near field quadrature corrections
!        * the first kernel is S_{k}
!        * the second kernel is \partial_{x} S_{k}
!        * the third kernel is \partial_{y} S_{k}
!        * the fourth kernel is \partial_{z} S_{k}
!    - zjvec: complex *16 (2,npts)
!        components of current density J in the orthonormal
!        frame defined in the notes above
!    - zrhovec: complex *16 (npts)
!        the charge density $\rho$
!    - novers: integer(npatches)
!        order of discretization for oversampled sources and
!        density
!    - ixyzso: integer(npatches+1)
!        ixyzso(i) denotes the starting location in srcover,
!        corresponding to patch i
!    - nptso: integer
!        total number of oversampled points
!    - srcover: real *8 (12,nptso)
!        oversampled set of source information
!    - whtsover: real *8 (nptso)
!        smooth quadrature weights at oversampled nodes
!
!  Output arguments
!    - E: complex *16 (3,ntarg)
!         Electric field at the target locations
!    - H: complex *16 (3,ntarg)
!         Magnetic field at the target locations
!               
!
!
!
!
      implicit real *8 (a-h,o-z)
!  calling sequence variables      
      integer, intent(in) :: npatches,norders(npatches)
      integer, intent(in) :: ixyzs(npatches+1),iptype(npatches),npts
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      integer, intent(in) :: ndtarg,ntarg
      real *8, intent(in) :: targs(ndtarg,ntarg)
      complex *16, intent(in) :: zpars
      integer, intent(in) :: nnz,row_ptr(ntarg+1),col_ind(nnz)
      integer, intent(in) :: iquad(nnz+1),nquad
      complex *16, intent(in) :: wnear(nquad,4),zjvec(2,npts),rho(npts)
      integer, intent(in) :: novers(npatches),nptso,ixyzso(npatches+1)
      real *8, intent(in) :: srcover(12,nptso),whtsover(nptso)
      complex *16, intent(out) :: E(3,npts),H(3,npts)


!  Temporary variables
      real *8, allocatable :: uvec(:,:),vvec(:,:)
      complex *16, allocatable :: zjuse(:,:),zcurlj(:,:),abc0(:,:)
      complex *16, allocatable :: sigmaover(:,:),charges(:,:)
      complex *16, allocatable :: pot_aux(:,:),grad_aux(:,:,:)
      
      real *8, allocatable :: sources(:,:),targtmp(:,:)
      complex *16, allocatable :: pottmp(:),gradtmp(:,:),ctmp0(:,:)
      real *8, allocatable :: srctmp2(:,:)

      real *8 done,pi,over4pi
      complex *16 ima
      complex *16 w0,w1,w2,w3,w4,ztmp(3)
      


      data over4pi/0.07957747154594767d0/
      data ima/(0.0d0,1.0d0)/

      ns = nptso

      done = 1.0d0
      pi = atan(done)*4.0d0

      allocate(uvec(3,npts),vvec(3,npts))
      call orthonormalize_all(srcvals(4:6,1:npts),srcvals(10:12,1:npts),&
          uvec,vvec)

      ifpgh = 0
      ifpghtarg = 2
      allocate(sources(3,ns),targtmp(3,ntarg))
      nmax = 0
!
!  estimate max number of sources in the near field of any target
!
!    
      call get_near_corr_max(npts,row_ptr,nnz,col_ind,npatches,ixyzso,&
        nmax)

!
!  Allocate various densities
!

      allocate(sigmaover(4,ns),zjuse(3,npts),abc0(4,npts))

!
!  extract source and target info

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)
      enddo
!$OMP END PARALLEL DO      


!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        zjuse(1,i) = uvec(1,i)*zjvec(1,i) + vvec(1,i)*zjvec(2,i)
        zjuse(2,i) = uvec(2,i)*zjvec(1,i) + vvec(2,i)*zjvec(2,i)
        zjuse(3,i) = uvec(2,i)*zjvec(1,i) + vvec(3,i)*zjvec(2,i)
        abc0(1:3,i) = zjuse(1:3,i)
        abc0(4,i) = zrhovec(i) 
      enddo
!$OMP END PARALLEL DO


      nd = 4
      allocate(charges0(nd,ns))

      nd2 = nd*2
      call oversample_fun_surf(nd2,npatches,norders,ixyzs,iptype,& 
          npts,zjuse,novers,ixyzso,ns,sigmaover)

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        targtmp(1:3,i) = targs(1:3,i)
      enddo
!$OMP END PARALLEL DO

        
!
!
!$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,ns
        charges0(1:4,i) = sigmaover(1:4,i)*whtsover(i)*over4pi
      enddo
!$OMP END PARALLEL DO      

      allocate(pot_aux(nd,ntarg),grad_aux(nd,3,ntarg))

!      print *, "before fmm"

      ier = 0
      call hfmm3d_t_c_g_vec(nd,zpars,eps,ns,sources,charges0,npts, &
        targtmp,pot_aux,grad_aux,ier)

!$OMP PARALLEL DO DEFAULT(SHARED)         
      do i=1,ntarg
        H(1,i) = grad_aux(3,2,i) - grad_aux(2,3,i)
        H(2,i) = grad_aux(1,3,i) - grad_aux(3,1,i) 
        H(3,i) = grad_aux(2,1,i) - grad_aux(1,2,i)
        E(1:3,i) = pot_aux(1:3,i)*ima*zpars - grad_aux(4,1:3,i)
      enddo
!$OMP END PARALLEL DO      

      
!  E = ik S_{k} [J] - \nabla S_{k} [\rho]
      
!
!  Add near quadrature correction
!

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart) &
!$OMP PRIVATE(jstart,npols,l,w1,w2,w3,w4)
      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            w1 = wnear(jquadstart+l-1,1)
            w2 = wnear(jquadstart+l-1,2)
            w3 = wnear(jquadstart+l-1,3)
            w4 = wnear(jquadstart+l-1,4)
            H(1,i) = H(1,i) + w3*abc0(3,jstart+l-1) - &
              w4*abc0(2,jstart+l-1)
            H(2,i) = H(2,i) + w4*abc0(1,jstart+l-1) - &
              w2*abc0(3,jstart+l-1)
            H(3,i) = H(3,i) + w2*abc0(2,jstart+l-1) - &
              w3*abc0(1,jstart+l-1)
            E(1:3,i) = E(1:3,i) +ima*zpars*w1*abc0(1:3,jstart+l-1)
            E(1,i) = E(1,i) - w2*abc0(4,jstart+l-1)
            E(2,i) = E(2,i) - w3*abc0(4,jstart+l-1)
            E(3,i) = E(3,i) - w4*abc0(4,jstart+l-1)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO     
      

      call get_fmm_thresh(12,ns,srcover,12,npts,srcvals,thresh)

!
! Subtract near contributions computed via fmm
!
      allocate(pottmp(nd),gradtmp(nd,3))
      allocate(ctmp0(nd,nmax),srctmp2(3,nmax))
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2) &
!$OMP PRIVATE(ctmp0,l,jstart,nss,pottmp,gradtmp)
      do i=1,ntarg
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
            nss = nss+1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)

            ctmp0(1:3,nss)=charges0(1:3,l)
          enddo
        enddo

        pottmp = 0
        gradtmp = 0

        call h3ddirectcg(nd,zpars(1),srctmp2,ctmp0,nss,targtmp(1,i), &
          ntarg0,pottmp,gradtmp,thresh)
        H(1,i) = H(1,i) - (gradtmp(3,2)-gradtmp(2,3))
        H(2,i) = H(2,i) - (gradtmp(1,3)-gradtmp(3,1))
        H(3,i) = H(3,i) - (gradtmp(2,1)-gradtmp(1,2))

        E(1:3,i) = E(1:3,i) - ima*zpars*pottmp(1:3)
        E(1:3,i) = E(1:3,i) + gradtmp(4,1:3)
        
      enddo
!$OMP END PARALLEL DO      

      

      return
      end subroutine lpcomp_em_aumfie_EH_addsub
!
!
!
