
      subroutine getnearquad_helm_comb_trans(npatches,norders,&
       ixyzs,iptype,npts,srccoefs,srcvals, &
       eps,zpars,iquadtype,nnz,row_ptr,col_ind,&
       iquad,rfac0,nquad,wnear)
!
!  This subroutine generates the near field quadrature
!  for the representations:
!
!  u1 = ep1^2 S_{k1}[\lambda]+ep1*D_{k1}[\rho] (interior representation)
!  u0 = ep0^2 S_{k0}[\lambda]+ep0*D_{k0}[\rho] (exterior representation)
!
!  and returns quantities related to u0-u1 , and 1/ep0 u0' - 1/ep1 u1'
!  on surface
!
!  On imposing the boundary condition, we get the following 
!  sets of operators
!
!  u0-u1 = (ep0+ep1)/2 \rho + (ep0 D_{k0} - \ep1 D_{k1})[\rho] +
!    (ep0^2 S_{k1} - ep1^2 S_{k0})[\lambda]
!  
!  1/ep0 u0' - 1/ep1 u1' = -(ep0 + ep1)/2 \lambda + 
!    (D_{k0}'-D_{k1}')[\rho] + (ep0 S_{k0}' - ep1 S_{k1}')[\lambda]
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
!    - wnear must be of size 4*nquad as 4 different layer
!      potentials are returned
!      * the first kernel is ep0^2 S_{k0} - ep1^2 S_{k1}
!      * the second kernel is ep0 D_{k0} - ep1 D_{k1}
!      * the third kernel is ep0 S_{k0}' - ep1 S_{k1}'
!      * the fourth kernel is D_{k0}'-D_{k1}'
!    - The parameters that must be provided zpars(5), and how they
!      are related to the problem parameters
!        zpars(1) = \omega
!        zpars(2) = ep0
!        zpars(3) = mu0
!        zpars(4) = ep1
!        zpars(5) = mu1
!        
!        k0  = \omega*sqrt(ep0*mu0)
!        k1 = \omega*sqrt(ep1*mu1)
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
!    - eps: real *8
!        precision requested
!    - zpars: complex *16 (5)
!        kernel parameters (See notes above)
!        zpars(1) = \omega 
!        zpars(2) = ep0
!        zpars(3) = mu0
!        zpars(4) = ep1
!        zpars(5) = mu1
!    - iquadtype: integer
!        quadrature type
!          * iquadtype = 1, use ggq for self + adaptive integration
!            for rest
!    - nnz: integer
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer(npts+1)
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
!    - wnear: complex *16(4*nquad)
!        The desired near field quadrature
!               
!

      implicit none 
      integer npatches,norders(npatches),npts,nquad
      integer ixyzs(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,rfac0
      integer ndtarg,ntarg
      integer iquadtype
      complex *16 zpars(5)
      complex *16 zpars_tmp(5)
      complex *16 zk0,zk1
      integer nnz,ipars(2)
      real *8 dpars(1)
      integer row_ptr(npts+1),col_ind(nnz),iquad(nnz+1)
      complex *16 wnear(4*nquad)

      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)


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
!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo
!$OMP END PARALLEL DO      

      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,&
        ipatch_id,uvs_targ)
      ipv=0
      fker => h3d_slp_diff
      zpars_tmp(3) = zpars(2)**2
      zpars_tmp(4) = zpars(4)**2
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,srcvals,&
        ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars_tmp,&
        ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear)

      fker => h3d_dlp_diff
      zpars_tmp(3) = zpars(2)
      zpars_tmp(4) = zpars(4)
      ipv = 1
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs, &
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,srcvals, &
        ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars_tmp, &
        ndi,ipars,nnz,row_ptr,col_ind,iquad, &
        rfac0,nquad,wnear(nquad+1))
      
      fker => h3d_sprime_diff
      ipv = 1
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
       iptype,npts,srccoefs,srcvals,ndtarg,ntarg,srcvals,&
       ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars_tmp,ndi,&
       ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(2*nquad+1))

      ndz = 2
      fker => h3d_dprime_diff
      ipv = 0
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
       iptype,npts,srccoefs,srcvals,ndtarg,ntarg,srcvals,&
       ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars_tmp,ndi,&
       ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(3*nquad+1))

      return
      end subroutine getnearquad_helm_comb_trans
!
!
!
!
!
!

      subroutine lpcomp_helm_comb_trans_addsub(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,eps,zpars,nnz,row_ptr,col_ind, &
        iquad,nquad,wnear,sigma,novers,nptso,ixyzso,srcover,whtsover, &
        pot)

!
!  This subroutine evaluates the transmission data corresponding to
!  the following integral representation:
!
!  u1 = ep1^2 S_{k1}[\lambda]+ep1*D_{k1}[\rho] (interior representation)
!  u0 = ep0^2 S_{k0}[\lambda]+ep0*D_{k0}[\rho] (exterior representation)
!
!  and returns quantities related to u0-u1 , and 1/ep0 u0' - 1/ep1 u1'
!  on surface
!
!  On imposing the boundary condition, we get the following 
!  sets of operators
!
!  u0-u1 = (ep0+ep1)/2 \rho + (ep0 D_{k0} - \ep1 D_{k1})[\rho] +
!    (ep0^2 S_{k1} - ep1^2 S_{k0})[\lambda]
!  
!  1/ep0 u0' - 1/ep1 u1' = -(ep0 + ep1)/2 \lambda + 
!    (D_{k0}'-D_{k1}')[\rho] + (ep0 S_{k0}' - ep1 S_{k1}')[\lambda]
!
!  where the near field is precomputed and stored
!  in the row sparse compressed format.
!
!  NOTES: 
!    - on output, the identity terms are not included
!    - We return 1/ep1 u1' - 1/ep0 u0' to align the identity
!      terms in the two sets of equations.
!    - The parameters that must be provided zpars(5), and how they
!      are related to the problem parameters
!        zpars(1) = \omega
!        zpars(2) = ep0
!        zpars(3) = mu0
!        zpars(4) = ep1
!        zpars(5) = mu1
!        
!        k0  = \omega*sqrt(ep0*mu0)
!        k1 = \omega*sqrt(ep1*mu1)
! 
!  Input arguments:
!
!  The fmm is used to accelerate the far-field and 
!  near-field interactions are handled via precomputed quadrature
!
!  Using add and subtract - no need to call tree and set fmm parameters
!  can directly call existing fmm library
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
!    - eps: real *8
!        precision requested
!    - zpars: complex *16 (5)
!        kernel parameters (See notes above)
!        zpars(1) = \omega 
!        zpars(2) = ep0
!        zpars(3) = mu0
!        zpars(4) = ep1
!        zpars(5) = mu1
!    - nnz: integer *8
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer(npts+1)
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
!    - wnear: complex *16(4*nquad)
!        Precomputed near field quadrature
!          * the first kernel is ep0^2 S_{k0} - ep1^2 S_{k1}
!          * the second kernel is ep0 D_{k0} - ep1 D_{k1}
!          * the third kernel is ep0 S_{k0}' - ep1 S_{k1}'
!          * the fourth kernel is D_{k0}'-D_{k1}'
!    - sigma: complex *16(2*npts)
!        * sigma(1:npts) is the density \lambda
!        * sigma(npts+1:2*npts) is the density \rho
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
!  Output arguments:
!    - pot: complex *16 (2*npts)
!        * pot(1:npts) is the data u0-u1
!        * pot(npts+1:2*npts) is the pv part of the 
!          data 1/ep1 u1' - 1/ep0u0'
!

      implicit none
      integer npatches,norder,npols,npts
      integer ndtarg,ntarg
      integer norders(npatches),ixyzs(npatches+1)
      integer ixyzso(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps
      complex *16 zpars(5),zpars_tmp(4)
      complex *16 zk0,zk1
      integer nnz,row_ptr(npts+1),col_ind(nnz),nquad
      integer iquad(nnz+1)
      complex *16 sigma(2*npts)
      complex *16 wnear(4*nquad)

      integer novers(npatches)
      integer nover,npolso,nptso
      real *8 srcover(12,nptso),whtsover(nptso)
      complex *16 pot(2*npts)

      real *8, allocatable :: sources(:,:),srctmp(:,:)
      complex *16, allocatable :: charges0(:),dipvec0(:,:),sigmaover(:)
      complex *16, allocatable :: charges1(:),dipvec1(:,:)
      integer ns,nt
      complex *16 alpha,beta
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      complex *16 tmp(10),val,E(4)

      integer i,j,jpatch,jquadstart,jstart,count1,count2
      complex *16 zdotu,pottmp,gradtmp(3)
      complex *16, allocatable :: pot_aux(:),grad_aux(:,:)
      complex *16 ep0,ep1,ep0sq,ep1sq,ep0inv,ep1inv

      real *8 radexp,epsfmm

      integer ipars(2)
      real *8 dpars(1),timeinfo(10),t1,t2,omp_get_wtime

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      complex *16, allocatable :: ctmp0(:),dtmp0(:,:)
      complex *16, allocatable :: ctmp1(:),dtmp1(:,:)
      real *8 thresh,ra
      real *8 rr,rmin
      real *8 over4pi
      integer nss,ii,l,npover,ier
      complex *16 ima,ztmp

      integer nd,ntarg0,nmax
      integer ndd,ndz,ndi

      real *8 ttot,done,pi
      data ima/(0.0d0,1.0d0)/
      data over4pi/0.07957747154594767d0/

      parameter (nd=1,ntarg0=1)

      ns = nptso
      ntarg = npts
      done = 1
      pi = atan(done)*4

!
!  compute interior and exteiror wave numbers
!
      zk0 = zpars(1)*sqrt(zpars(2)*zpars(3))
      zk1 = zpars(1)*sqrt(zpars(4)*zpars(5))
      ep0 = zpars(2)
      ep1 = zpars(4)
      ep0sq = ep0*ep0
      ep1sq = ep1*ep1
      ep0inv = 1.0d0/ep0
      ep1inv = 1.0d0/ep1
           
      ifpgh = 0
      ifpghtarg = 2
      allocate(sources(3,ns),srctmp(3,npts))
      allocate(charges0(ns),dipvec0(3,ns))
      allocate(charges1(ns),dipvec1(3,ns))
      allocate(sigmaover(2*ns))
      allocate(pot_aux(npts),grad_aux(3,npts))
!
!  estimate max number of sources in the near field of any target
!
!    
      call get_near_corr_max(npts,row_ptr,nnz,col_ind,npatches,ixyzso,&
        nmax)
      allocate(srctmp2(3,nmax),ctmp0(nmax),dtmp0(3,nmax))
      allocate(ctmp1(nmax),dtmp1(3,nmax))
! 
!       oversample densities
!
      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
        npts,sigma,novers,ixyzso,ns,sigmaover)
      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
        npts,sigma(npts+1),novers,ixyzso,ns,sigmaover(ns+1))

!
!  extract source and target info

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)
        charges0(i) = sigmaover(ns+i)*whtsover(i)*over4pi*ep0sq
        charges1(i) = sigmaover(ns+i)*whtsover(i)*over4pi*ep1sq

        dipvec0(1,i) = sigmaover(i)*whtsover(i)*over4pi* &
           ep0*srcover(10,i)
        dipvec0(2,i) = sigmaover(i)*whtsover(i)*over4pi* &
           ep0*srcover(11,i)
        dipvec0(3,i) = sigmaover(i)*whtsover(i)*over4pi* &
           ep0*srcover(12,i)

        dipvec1(1,i) = sigmaover(i)*whtsover(i)*over4pi* &
           ep1*srcover(10,i)
        dipvec1(2,i) = sigmaover(i)*whtsover(i)*over4pi* &
           ep1*srcover(11,i)
        dipvec1(3,i) = sigmaover(i)*whtsover(i)*over4pi* &
           ep1*srcover(12,i)
      enddo
!$OMP END PARALLEL DO      


!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        srctmp(1,i) = srcvals(1,i)
        srctmp(2,i) = srcvals(2,i)
        srctmp(3,i) = srcvals(3,i)
        pot_aux(i) = 0
        grad_aux(1,i) = 0
        grad_aux(2,i) = 0
        grad_aux(3,i) = 0
        pot(i) = 0
        pot(npts+i) = 0
      enddo

!$OMP END PARALLEL DO

! 
!  Compute ep0^2*S_{k0}[\lambda] + ep0*D_{k0}[\rho], 
!  ep0^2*S_{k0}'[\lambda] + ep0*D_{k0}'[\rho]

      call hfmm3d_t_cd_g(eps,zk0,ns,sources,charges0,dipvec0,npts, &
        srctmp,pot_aux,grad_aux,ier)

!
!   Add ep0^2*S_{k0}[\lambda]+ep0 D_{k0}[\rho] to pot and 
!     1/ep0*(ep0^2*S_{k0}' + ep0 D_{k0}[\rho] to pot(npts+:)
!
!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        pot(i) = pot_aux(i) 
        pot(npts+i) = ep0inv*(grad_aux(1,i)*srcvals(10,i)+ &
          grad_aux(2,i)*srcvals(11,i)+ &
          grad_aux(3,i)*srcvals(12,i))
      enddo
!$OMP END PARALLEL DO      
    
! 
!  Compute ep1^2*S_{k1}[\lambda] + ep1*D_{k1}[\rho], 
!  ep1^2*S_{k1}'[\lambda] + ep1*D_{k1}'[\rho]
!
      call hfmm3d_t_cd_g(eps,zk1,ns,sources,charges1,dipvec1,npts, &
        srctmp,pot_aux,grad_aux,ier)

!
!   subtract ep1^2*S_{k1}[\lambda]+ep1 D_{k1}[\rho] from pot and 
!     1/ep1*(ep1^2*S_{k1}' + ep0 D_{k1}[\rho] from pot(npts+:)
!
!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        pot(i) = pot(i)- pot_aux(i) 
        pot(npts+i) = pot(npts+i)-ep1inv*(grad_aux(1,i)*srcvals(10,i)+&
          grad_aux(2,i)*srcvals(11,i)+ &
          grad_aux(3,i)*srcvals(12,i))
      enddo
!$OMP END PARALLEL DO      

!        compute threshold for ignoring local computation
      call get_fmm_thresh(12,ns,srcover,12,npts,srcvals,thresh)
!
!       Add near field precomputed contribution
!
      call cpu_time(t1)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart) &
!$OMP PRIVATE(jstart,npols,l)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            pot(i) = pot(i) + wnear(jquadstart+l-1)*sigma(jstart+l-1+npts)
            pot(i) = pot(i) + wnear(nquad+jquadstart+l-1)* &
              sigma(jstart+l-1)
            pot(i+npts) = pot(i+npts) + wnear(2*nquad+jquadstart+l-1)* &
              sigma(jstart+l-1+npts)
            pot(i+npts) = pot(i+npts) + wnear(3*nquad+jquadstart+l-1)* &
              sigma(jstart+l-1)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

!     Remove near contribution of the FMM
!

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2) &
!$OMP PRIVATE(ctmp0,dtmp0,ctmp1,dtmp1,l,jstart,nss,pottmp,gradtmp)
      do i=1,npts
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
            nss = nss+1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)

            ctmp0(nss)=charges0(l)
            ctmp1(nss)=charges1(l)

            dtmp0(1,nss)=dipvec0(1,l)
            dtmp0(2,nss)=dipvec0(2,l)
            dtmp0(3,nss)=dipvec0(3,l)

            dtmp1(1,nss)=dipvec1(1,l)
            dtmp1(2,nss)=dipvec1(2,l)
            dtmp1(3,nss)=dipvec1(3,l)
          enddo
        enddo

        pottmp = 0
        gradtmp(1) = 0
        gradtmp(2) = 0
        gradtmp(3) = 0

        call h3ddirectcdg(nd,zk0,srctmp2,ctmp0,dtmp0,nss,srctmp(1,i), &
          ntarg0,pottmp,gradtmp,thresh)
        pot(i) = pot(i)  - pottmp   
        pot(npts+i) = pot(npts+i) - (gradtmp(1)*srcvals(10,i)+  &
          gradtmp(2)*srcvals(11,i)+gradtmp(3)*srcvals(12,i))*ep0inv

        pottmp = 0
        gradtmp(1) = 0
        gradtmp(2) = 0
        gradtmp(3) = 0

        call h3ddirectcdg(nd,zk1,srctmp2,ctmp1,dtmp1,nss,srctmp(1,i), &
          ntarg0,pottmp,gradtmp,thresh)
        pot(i) = pot(i)  + pottmp   
        pot(npts+i) = pot(npts+i) + (gradtmp(1)*srcvals(10,i)+  &
          gradtmp(2)*srcvals(11,i)+gradtmp(3)*srcvals(12,i))*ep1inv
!
!  flip sign of pot(npts+i)
! 
        pot(npts+i) = -pot(npts+i)
      enddo
!$OMP END PARALLEL DO      


      return
      end subroutine lpcomp_helm_comb_trans_addsub
!
!
!
!
!
!

      subroutine helm_comb_trans_solver(npatches,norders,ixyzs, &
        iptype,npts,srccoefs,srcvals,eps,zpars,numit, &
        rhs,eps_gmres,niter,errs,rres,soln)
!
!
!  This subroutine solves the Helmholtz tranmission problem
!  where the potential is represented as an appropriate
!  combined field integral representation in each region
!
!  PDE:
!    \Delta u1 + k1^2 u1 = 0 (interior region)
!    \Delta u0 + k0^2 u0 = 0 (exterior region)
!
!  Representations:
!    u1 = ep1^2 S_{k1}[\lambda]+ep1*D_{k1}[\rho] 
!    u0 = ep0^2 S_{k0}[\lambda]+ep0*D_{k0}[\rho] 
!
!  Boundary condition:
!    u0-u1=f
!    (1/ep0)u0' - (1/ep1)u1' = g
!
!  The linear system is solved iteratively using GMRES
!  until a relative residual of eps_gmres is reached
!
!  NOTES:
!    - The parameters that must be provided zpars(5), and how they
!      are related to the problem parameters
!        zpars(1) = \omega
!        zpars(2) = ep0
!        zpars(3) = mu0
!        zpars(4) = ep1
!        zpars(5) = mu1
!        
!        k0  = \omega*sqrt(ep0*mu0)
!        k1 = \omega*sqrt(ep1*mu1)
!
!  Input:
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
!    - eps: real *8
!        precision requested for computing quadrature and fmm
!        tolerance
!    - zpars: complex *16 (5)
!        kernel parameters (See notes above)
!        zpars(1) = \omega 
!        zpars(2) = ep0
!        zpars(3) = mu0
!        zpars(4) = ep1
!        zpars(5) = mu1
!    - numit: integer
!        max number of gmres iterations
!    - rhs: complex *16(2*npts)
!        boundary data
!          * rhs(1:npts) = f
!          * rhs(npts+1:2*npts) = g
!    - eps_gmres: real *8
!        gmres tolerance requested
!      
!
!  output
!    - niter: integer
!        number of gmres iterations required for relative residual
!        to converge to eps_gmres
!    - errs:  real *8 (numit)
!        relative residual as a function of iteration
!        number (only errs(1:niter) is populated))
!    - rres: real *8
!        relative residual for computed solution
!    - soln: complex *16(2*npts)
!        densities which solve the transmission problem
!        soln(1:npts) = \lambda
!        soln(npts+1:2*npts) = \rho
!
      implicit none
      integer npatches,norder,npols,npts
      integer norders(npatches),ixyzs(npatches+1)
      integer iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,eps_gmres
      complex *16 zpars(5)
      complex *16 rhs(2*npts)
      complex *16 soln(2*npts)

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

      integer ipars(2)
      real *8 dpars,timeinfo(10),t1,t2,omp_get_wtime
      complex *16 zk0,zk1,zkuse


      real *8 ttot,done,pi
      real *8 rfac,rfac0
      integer iptype_avg,norder_avg
      integer ikerorder, iquadtype,npts_over
      integer n_var

!
!
!       gmres variables
!
      complex *16 zid,ztmp
      real *8 rb,wnrm2
      integer numit,it,iind,it1,k,l,count1
      real *8 rmyerr
      complex *16 temp
      complex *16, allocatable :: vmat(:,:),hmat(:,:)
      complex *16, allocatable :: cs(:),sn(:)
      complex *16, allocatable :: svec(:),yvec(:),wtmp(:)
      complex *16 ima

      ima=(0.0d0,1.0d0)

!
!   n_var is the number of unknowns in the linear system.
!   in this case n_var=2*npts
!

      n_var=npts*2


      allocate(vmat(n_var,numit+1),hmat(numit,numit))
      allocate(cs(numit),sn(numit))
      allocate(wtmp(n_var),svec(numit+1),yvec(numit+1))


      done = 1
      pi = atan(done)*4


!
!
!        setup targets as on surface discretization points
! 
      ndtarg = 12  
      ntarg = npts
      allocate(targs(ndtarg,npts),uvs_targ(2,ntarg),ipatch_id(ntarg))
!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        targs(:,i)=srcvals(:,i)
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo
!$OMP END PARALLEL DO   


!
!    initialize patch_id and uv_targ for on surface targets
!
      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,&
        ipatch_id,uvs_targ)

!
!
!        this might need fixing
!
      iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
      norder_avg = floor(sum(norders)/(npatches+0.0d0))

      call get_rfacs(norder_avg,iptype_avg,rfac,rfac0)

      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))

      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts,& 
        srccoefs,cms,rads)

!$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,npatches
        rad_near(i) = rads(i)*rfac
      enddo
!$OMP END PARALLEL DO      

!
!    find near quadrature correction interactions
!
      print *, "entering find near mem"
      call findnearmem(cms,npatches,rad_near,ndtarg,targs,npts,nnz)
      print *, "nnz=",nnz

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,ndtarg,targs,npts,row_ptr,&
         col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind,&
         iquad)

      ikerorder = 0
!
!    estimate oversampling for far-field, and oversample geometry
!

      allocate(novers(npatches),ixyzso(npatches+1))

      print *, "beginning far order estimation"

      zk0 = zpars(1)*sqrt(zpars(2)*zpars(3))
      zk1 = zpars(1)*sqrt(zpars(4)*zpars(5))
      
      zkuse = real(zk0)
      if(real(zk1).gt.real(zk0)) zkuse = real(zk1)

      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,&
       rads,npts,srccoefs,ndtarg,npts,targs,ikerorder,zkuse,&
       nnz,row_ptr,col_ind,rfac,novers,ixyzso)

      npts_over = ixyzso(npatches+1)-1
      print *, "npts_over=",npts_over

      allocate(srcover(12,npts_over),wover(npts_over))

      call oversample_geom(npatches,norders,ixyzs,iptype,npts,&
        srccoefs,srcvals,novers,ixyzso,npts_over,srcover)

      call get_qwts(npatches,novers,ixyzso,iptype,npts_over,&
        srcover,wover)
!
!   compute near quadrature correction
!
      nquad = iquad(nnz+1)-1
      print *, "nquad=",nquad
      allocate(wnear(4*nquad))
      
!$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,4*nquad
        wnear(i)=0
      enddo
!$OMP END PARALLEL DO    

      iquadtype = 1

      print *, "starting to generate near quadrature"
      call cpu_time(t1)
!$      t1 = omp_get_wtime()      

      call getnearquad_helm_comb_trans(npatches,norders,&
       ixyzs,iptype,npts,srccoefs,srcvals, &
       eps,zpars,iquadtype,nnz,row_ptr,col_ind,&
       iquad,rfac0,nquad,wnear)


      call cpu_time(t2)
!$      t2 = omp_get_wtime()     

      call prin2('quadrature generation time=*',t2-t1,1)
      
      print *, "done generating near quadrature, now starting gmres"

!
!
!     start gmres code here
!
!     NOTE: matrix equation should be of the form (z*I + K)x = y
!       the identity scaling (z) is defined via zid below,
!       and K represents the action of the principal value 
!       part of the matvec
!
  
      zid= (zpars(2)+zpars(4))/2.0d0


      niter=0
!
!      compute norm of right hand side and initialize v
! 
      rb = 0

      do i=1,numit
        cs(i) = 0
        sn(i) = 0
      enddo

!
!$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rb)
      do i=1,n_var
        rb = rb + abs(rhs(i))**2
      enddo
!$OMP END PARALLEL DO      
      rb = sqrt(rb)

!
!  Note: sign in g needs to be flipped since we are instead imposing
!  1/ep1 u1' - 1/ep0 u0' = -g for the identity terms to line up
!

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        vmat(i,1) = rhs(i)/rb
        vmat(i+npts,1)=-rhs(i+npts)/rb
      enddo
!$OMP END PARALLEL DO      

      svec(1) = rb

      do it=1,numit
        it1 = it + 1

!
!        NOTE:
!        replace this routine by appropriate layer potential
!        evaluation routine  
!

      call lpcomp_helm_comb_trans_addsub(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,eps,zpars,nnz,row_ptr,col_ind, &
        iquad,nquad,wnear,vmat(1,it),novers,npts_over,ixyzso,srcover,&
        wover,wtmp)

        do k=1,it
          ztmp = 0
!$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:ztmp)          
          do j=1,n_var      
            ztmp = ztmp + wtmp(j)*conjg(vmat(j,k))
          enddo
!$OMP END PARALLEL DO  
          hmat(k,it) = ztmp

!$OMP PARALLEL DO DEFAULT(SHARED)
          do j=1,n_var
            wtmp(j) = wtmp(j)-hmat(k,it)*vmat(j,k)
          enddo
!$OMP END PARALLEL DO          
        enddo
          
        hmat(it,it) = hmat(it,it)+zid
        wnrm2 = 0
!$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:wnrm2)        
        do j=1,n_var
          wnrm2 = wnrm2 + abs(wtmp(j))**2
        enddo
!$OMP END PARALLEL DO        
        wnrm2 = sqrt(wnrm2)

!$OMP PARALLEL DO DEFAULT(SHARED)
        do j=1,n_var
          vmat(j,it1) = wtmp(j)/wnrm2
        enddo
!$OMP END PARALLEL DO        

        do k=1,it-1
          temp = cs(k)*hmat(k,it)+conjg(sn(k))*hmat(k+1,it)
          hmat(k+1,it) = -sn(k)*hmat(k,it)+cs(k)*hmat(k+1,it)
          hmat(k,it) = temp
        enddo

        ztmp = wnrm2

        call zrotmat_gmres(hmat(it,it),ztmp,cs(it),sn(it))
          
        hmat(it,it) = cs(it)*hmat(it,it)+conjg(sn(it))*wnrm2
        svec(it1) = -sn(it)*svec(it)
        svec(it) = cs(it)*svec(it)
        rmyerr = abs(svec(it1))/rb
        errs(it) = rmyerr
        print *, "iter=",it,errs(it)

        if(rmyerr.le.eps_gmres.or.it.eq.numit) then

!
!            solve the linear system corresponding to
!            upper triangular part of hmat to obtain yvec
!
!            y = triu(H(1:it,1:it))\s(1:it);
!
          do j=1,it
            iind = it-j+1
            yvec(iind) = svec(iind)
            do l=iind+1,it
              yvec(iind) = yvec(iind) - hmat(iind,l)*yvec(l)
            enddo
            yvec(iind) = yvec(iind)/hmat(iind,iind)
          enddo



!
!          estimate x
!
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
          do j=1,n_var
            soln(j) = 0
            do i=1,it
              soln(j) = soln(j) + yvec(i)*vmat(j,i)
            enddo
          enddo
!$OMP END PARALLEL DO          


          rres = 0
!$OMP PARALLEL DO DEFAULT(SHARED)          
          do i=1,n_var
            wtmp(i) = 0
          enddo
!$OMP END PARALLEL DO          
!
!        NOTE:
!        replace this routine by appropriate layer potential
!        evaluation routine  
!

      call lpcomp_helm_comb_trans_addsub(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,eps,zpars,nnz,row_ptr,col_ind, &
        iquad,nquad,wnear,soln,novers,npts_over,ixyzso,srcover,&
        wover,wtmp)

!$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rres)
          do i=1,npts
            rres = rres + abs(zid*soln(i) + wtmp(i)-rhs(i))**2
            rres = rres + abs(zid*soln(npts+i) + &
              wtmp(npts+i)+rhs(npts+i))**2
          enddo
!$OMP END PARALLEL DO          
          rres = sqrt(rres)/rb
          niter = it
          return

        endif
      enddo
!
      return
      end subroutine helm_comb_trans_solver
!
!
!
      subroutine lpcomp_helm_comb_split_dir(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,ipatch_id, &
        uvs_targ,eps,zpars,sigma,pot)
!
!  This subroutine evaluates the dirichlet data for
!  the combined field representation where the density of the single and
!  double layer potential are allowed to be different
!
!  Representation:
!    u = \alpha S_{k}[\lambda]+beta*D_{k}[\rho]
!
!  v1: we implement using two calls to the combined field dirichlet 
!  representation. 
!
!  Note: while the number of fmm calls can be reduced from 2->1,
!  the bulk of the computation cost is the separate near quadarture
!  generation of S_{k} and D_{k}.
!
!  A faster version of this routine can be implemented by separating
!  the near quadrature generation, and then just using 1 fmm
!
!  Input:
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
!        leading dimension of target array
!    - ntarg: integer
!        number of targets
!    - targs: real *8 (ndtarg,ntarg)
!        target info, the first three coordinates must be
!        the xyz components of the target
!    - ipatch_id: integer(ntarg)
!        patch on which target is on if on-surface, =-1/0 if 
!        target is off-surface
!    - uvs_targ: real *8 (2,ntarg)
!        local uv coordinates on patch, if target is on surface
!    - eps: real *8
!        precision requested for computing quadrature and fmm
!        tolerance
!    - zpars: complex *16 (3)
!        kernel parameters (Referring to formula (1))
!          * zpars(1) = k 
!          * zpars(2) = alpha
!          * zpars(3) = beta
!    - sigma: complex *16(2*npts)
!        density sigma above
!
!  Output arguments:
!    - pot: complex *16(ntarg)
!        potential at the target locations
!				 
!
      implicit none
      integer, intent(in) :: npatches
      integer, intent(in) :: norders(npatches),ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches),npts
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
      integer, intent(in) :: ndtarg,ntarg
      real *8, intent(in) :: targs(ndtarg,ntarg)
      integer, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg),eps
      complex *16, intent(in) :: zpars(3),sigma(2*npts)
      complex *16, intent(out) :: pot(ntarg)

      complex *16, allocatable :: pottmp2(:)
      complex *16 zpars_tmp(3),ima
      integer, allocatable :: ipatch_id_src(:)
      real *8, allocatable :: uvs_src(:,:)
      integer i,ndtarg0

      data ima/(0.0d0,1.0d0)/


      allocate(pottmp2(ntarg))

!!!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        pottmp2(i) = 0
        pot(i) = 0
      enddo
!!!$OMP END PARALLEL DO
      
!
!  compute \alpha S_{k} [\lambda]
!
      zpars_tmp(1) = zpars(1)
      zpars_tmp(2) = zpars(2) 
      zpars_tmp(3) = 0

      
      call lpcomp_helm_comb_dir(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,ipatch_id, &
        uvs_targ,eps,zpars_tmp,sigma(npts+1),pot)
      
      
!
!  compute \beta D_{k} [\rho]
!
      zpars_tmp(1) = zpars(1)
      zpars_tmp(2) = 0
      zpars_tmp(3) = zpars(3)

      
      call lpcomp_helm_comb_dir(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,ipatch_id, &
        uvs_targ,eps,zpars_tmp,sigma,pottmp2)


!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        pot(i) = pot(i) + pottmp2(i)
      enddo
!$OMP END PARALLEL DO

     
      return
      end subroutine lpcomp_helm_comb_split_dir

