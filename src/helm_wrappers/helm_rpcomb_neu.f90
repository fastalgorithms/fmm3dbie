

      subroutine getnearquad_helm_rpcomb_neu(npatches,norders,&
       ixyzs,iptype,npts,srccoefs,srcvals, &
       eps,zpars,iquadtype,nnz,row_ptr,col_ind,&
       iquad,rfac0,nquad,wnear)
!
!  This subroutine generates the near field quadrature
!  for the representation:
!
!  u = S_{k}[\rho]+i*alpha*D_{k}[S_{i|k|}[\rho]]
!
!  and returns quantities related to evaluating du/dn on surface
!  at the surface discretization nodes
!
!  If values at other nodes is desired then the solution
!  should be reinterpolated to those nodes
!
!
!  On imposing the boundary condition, we get the following operator
!
!  du/dn = -I/2 + S_{k}' + i \alpha (D_{k}'-D_{i|k|}') S_{i|k|}
!    - i \alpha I/4 + i \alpha (S_{i|k|}')^2
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
!    - wnear must be of size 4*nquad as 4 different layer
!      potentials are returned
!      * the first kernel is S_{k}'
!      * the second kernel is S_{i|k|}
!      * the third kernel is S_{i|k|}'
!      * the fourth kernel is D_{k}'-D_{i|k|}'
! 
!  Input arguments:
! 
!    - npatches: integer(8)
!        number of patches
!    - norders: integer(8)(npatches)
!        order of discretization on each patch 
!    - ixyzs: integer(8)(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(8)(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!    - npts: integer(8)
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
!    - zpars: complex *16 (2)
!        kernel parameters (Referring to formula (1))
!        zpars(1) = k 
!        zpars(2) = alpha
!    - iquadtype: integer(8)
!        quadrature type
!          * iquadtype = 1, use ggq for self + adaptive integration
!            for rest
!    - nnz: integer(8)
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer(8)(npts+1)
!        row_ptr(i) is the pointer
!        to col_ind array where list of relevant source patches
!        for target i start
!    - col_ind: integer(8) (nnz)
!        list of source patches relevant for all targets, sorted
!        by the target number
!    - iquad: integer(8)(nnz+1)
!        location in wnear_ij array where quadrature for col_ind(i)
!        starts for a single kernel. In this case the different kernels
!        are matrix entries are located at (m-1)*nquad+iquad(i), where
!        m is the kernel number
!    - rfac0: integer(8)
!        radius parameter for near field
!    - nquad: integer(8)
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
      integer(8) npatches,norders(npatches),npts,nquad
      integer(8) ixyzs(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,rfac0
      integer(8) ndtarg,ntarg
      integer(8) iquadtype
      complex *16 zpars(2)
      complex *16 zpars_tmp(3)
      integer(8) nnz,ipars(2)
      real *8 dpars(1)
      integer(8) row_ptr(npts+1),col_ind(nnz),iquad(nnz+1)
      complex *16 wnear(4*nquad)

      real *8, allocatable :: uvs_targ(:,:)
      integer(8), allocatable :: ipatch_id(:)


      complex *16 alpha,beta,ima,zk
      integer(8) i,j,ndi,ndd,ndz

      integer(8) ipv

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
!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo
!$OMP END PARALLEL DO      

      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,&
        ipatch_id,uvs_targ)
      ipv=1
      fker => h3d_sprime 
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,srcvals,&
        ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
        ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear)

      zpars_tmp(1) = ima*abs(zk)
      fker => h3d_slp
      ipv = 1
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs, &
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,srcvals, &
        ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars_tmp, &
        ndi,ipars,nnz,row_ptr,col_ind,iquad, &
        rfac0,nquad,wnear(nquad+1))
      
      fker => h3d_sprime
      ipv = 1
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
       iptype,npts,srccoefs,srcvals,ndtarg,ntarg,srcvals,&
       ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars_tmp,ndi,&
       ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(2*nquad+1))

      zpars_tmp(1) = zk
      zpars_tmp(2) = ima*abs(zk)
      ndz = 2
      fker => h3d_dprime_diff
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
       iptype,npts,srccoefs,srcvals,ndtarg,ntarg,srcvals,&
       ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars_tmp,ndi,&
       ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(3*nquad+1))

      return
      end subroutine getnearquad_helm_rpcomb_neu
!
!
!
!
!
!

      subroutine lpcomp_helm_rpcomb_neu_addsub(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,eps,zpars,nnz,row_ptr,col_ind, &
        iquad,nquad,wnear,sigma,novers,nptso,ixyzso,srcover,whtsover, &
        pot,phi1)

!
!  This subroutine evaluates the neumann data corresponding to
!  the following integral representation:
!  
!  u = S_{k}[\rho]+i*alpha*D_{k}[S_{i|k|}[\rho]]
!
!  du/dn = S_{k}'[\rho] + i\alpha S_{i|k|}'^2 [\rho] + 
!    i \alpha (D_{k}' - D_{i|k|}') S_{i|k|}
!
!
!  where the near field is precomputed and stored
!  in the row sparse compressed format.
!
!
!  The fmm is used to accelerate the far-field and 
!  near-field interactions are handled via precomputed quadrature
!
!
!  Using add and subtract - no need to call tree and set fmm parameters
!  can directly call existing fmm library
!
!  Input arguments:
! 
!    - npatches: integer(8)
!        number of patches
!    - norders: integer(8)(npatches)
!        order of discretization on each patch 
!    - ixyzs: integer(8)(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(8)(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!    - npts: integer(8)
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
!    - zpars: complex *16 (2)
!        kernel parameters (Referring to formula (1))
!        zpars(1) = k 
!        zpars(2) = alpha
!    - nnz: integer(8)
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer(8)(npts+1)
!        row_ptr(i) is the pointer
!        to col_ind array where list of relevant source patches
!        for target i start
!    - col_ind: integer(8) (nnz)
!        list of source patches relevant for all targets, sorted
!        by the target number
!    - iquad: integer(8)(nnz+1)
!        location in wnear_ij array where quadrature for col_ind(i)
!        starts for a single kernel. In this case the different kernels
!        are matrix entries are located at (m-1)*nquad+iquad(i), where
!        m is the kernel number
!    - nquad: integer(8)
!        number of near field entries corresponding to each source target
!        pair. The size of wnear is 4*nquad since there are 4 kernels
!        per source target pair
!    - wnear: complex *16(4*nquad)
!        Precomputed near field quadrature
!          * the first kernel is S_{k}'
!          * the second kernel is S_{i|k|}
!          * the third kernel is S_{i|k|}'
!          * the fourth kernel is D_{k}'-D_{i|k|}'
!    - sigma: complex *16(npts)
!        density for layer potential
!    - novers: integer(8)(npatches)
!        order of discretization for oversampled sources and
!        density
!    - ixyzso: integer(8)(npatches+1)
!        ixyzso(i) denotes the starting location in srcover,
!        corresponding to patch i
!    - nptso: integer(8)
!        total number of oversampled points
!    - srcover: real *8 (12,nptso)
!        oversampled set of source information
!    - whtsover: real *8 (nptso)
!        smooth quadrature weights at oversampled nodes
!
!  Output arguments:
!    - pot: complex *16 (npts)
!        du/dn corresponding to representation
!    - phi1: complex *16 (npts)
!        sik(pot) to be used for computing solution later
!

      implicit none
      integer(8) npatches,norder,npols,npts
      integer(8) ndtarg,ntarg
      integer(8) norders(npatches),ixyzs(npatches+1)
      integer(8) ixyzso(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps
      complex *16 zpars(2),zpars2(3)
      integer(8) nnz,row_ptr(npts+1),col_ind(nnz),nquad
      integer(8) iquad(nnz+1)
      complex *16 sigma(npts)
      complex *16 wnear(4*nquad)

      integer(8) novers(npatches)
      integer(8) nover,npolso,nptso
      real *8 srcover(12,nptso),whtsover(nptso)
      complex *16 pot(npts),phi1(npts)

      real *8, allocatable :: sources(:,:),srctmp(:,:)
      complex *16, allocatable :: charges(:),dipvec(:,:),sigmaover(:)
      complex *16, allocatable :: sigmaover_aux(:),sigma_aux(:)
      integer(8) ns,nt
      complex *16 alpha,beta
      integer(8) ifcharge,ifdipole
      integer(8) ifpgh,ifpghtarg
      complex *16 tmp(10),val,E(1)


      integer(8) i,j,jpatch,jquadstart,jstart,count1,count2
      complex *16 zdotu,pottmp,gradtmp(3)
      complex *16, allocatable :: pot_aux(:),grad_aux(:,:)
      complex *16, allocatable :: phi2(:)

      real *8 radexp,epsfmm

      integer(8) ipars(2)
      real *8 dpars(1),timeinfo(10),t1,t2,omp_get_wtime

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      complex *16, allocatable :: ctmp2(:),dtmp2(:,:)
      real *8 thresh,ra
      real *8 rr,rmin
      real *8 over4pi
      integer(8) nss,ii,l,npover,ier
      complex *16 ima,ztmp

      integer(8) nd,ntarg0,nmax
      integer(8) ndd,ndz,ndi

      real *8 ttot,done,pi
      data ima/(0.0d0,1.0d0)/
      data over4pi/0.07957747154594767d0/

      parameter (nd=1,ntarg0=1)

      ns = nptso
      ntarg = npts
      done = 1
      pi = atan(done)*4


           
      ifpgh = 0
      ifpghtarg = 2
      allocate(sources(3,ns),srctmp(3,npts))
      allocate(charges(ns),dipvec(3,ns))
      allocate(sigmaover(ns))
      allocate(pot_aux(npts),grad_aux(3,npts))
      allocate(phi2(npts))
      allocate(sigma_aux(npts))
!
!  estimate max number of sources in the near field of any target
!
!    
      call get_near_corr_max(npts,row_ptr,nnz,col_ind,npatches,ixyzso,&
        nmax)
      allocate(srctmp2(3,nmax),ctmp2(nmax),dtmp2(3,nmax))
! 
!       oversample density
!
      call oversample_fun_surf(int(2,8),npatches,norders,ixyzs,iptype,& 
        npts,sigma,novers,ixyzso,ns,sigmaover)



!
!  extract source and target info
!
!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)
        charges(i) = sigmaover(i)*whtsover(i)*over4pi
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
      enddo

!$OMP END PARALLEL DO

! 
!  Compute S_{k}'
!
      call hfmm3d_t_c_g(eps,zpars(1),ns,sources,charges,npts, &
        srctmp,pot_aux,grad_aux,ier)

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        pot(i) = grad_aux(1,i)*srcvals(10,i)+ &
          grad_aux(2,i)*srcvals(11,i)+ &
          grad_aux(3,i)*srcvals(12,i)
      enddo
!$OMP END PARALLEL DO      
    
        
!        compute threshold for ignoring local computation
      call get_fmm_thresh(int(12,8),ns,srcover,int(12,8),npts,srcvals,thresh)
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
            pot(i) = pot(i) + wnear(jquadstart+l-1)*sigma(jstart+l-1)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

!
!     Remove near contribution of the FMM
!
      ifcharge = 1
      ifdipole = 0 

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2) &
!$OMP PRIVATE(ctmp2,l,jstart,nss,pottmp,gradtmp)
      do i=1,npts
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
            nss = nss+1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)
            ctmp2(nss)=charges(l)
          enddo
        enddo

        pottmp = 0
        gradtmp(1) = 0
        gradtmp(2) = 0
        gradtmp(3) = 0

        call h3ddirectcg(nd,zpars(1),srctmp2,ctmp2,nss,srctmp(1,i), &
          ntarg0,pottmp,gradtmp,thresh)
        pot(i) = pot(i) - gradtmp(1)*srcvals(10,i) - &
          gradtmp(2)*srcvals(11,i)-gradtmp(3)*srcvals(12,i)
      enddo
!$OMP END PARALLEL DO      

!
!
!    Now handle the computation of S_{i|k|}[\rho] and S_{i|k|}'[\rho] 
!    and hold them in separate arrays phi1 and phi2
!
!
      ztmp = ima*abs(zpars(1))
      call hfmm3d_t_c_g(eps,ztmp,ns,sources,charges,npts, &
        srctmp,pot_aux,grad_aux,ier)

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        phi1(i) = pot_aux(i)
        phi2(i) = grad_aux(1,i)*srcvals(10,i)+ &
          grad_aux(2,i)*srcvals(11,i)+ &
          grad_aux(3,i)*srcvals(12,i)
      enddo
!$OMP END PARALLEL DO      
    

!
!       Add near field precomputed contribution
!

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart) &
!$OMP PRIVATE(jstart,pottmp,npols,l)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            phi1(i) = phi1(i) + &
              wnear(nquad+jquadstart+l-1)*sigma(jstart+l-1)
            phi2(i) = phi2(i) + &
              wnear(2*nquad+jquadstart+l-1)*sigma(jstart+l-1)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

!
!     Remove near contribution of the FMM
!

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2) &
!$OMP PRIVATE(ctmp2,nss,l,pottmp,gradtmp)
      do i=1,npts
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
            nss = nss+1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)
            ctmp2(nss)= charges(l) 
          enddo
        enddo
        pottmp = 0
        gradtmp(1) = 0
        gradtmp(2) = 0
        gradtmp(3) = 0
        call h3ddirectcg(nd,ztmp,srctmp2,ctmp2,nss,srctmp(1,i), &
          ntarg0,pottmp,gradtmp,thresh)
        phi1(i) = phi1(i) - pottmp  
        phi2(i) = phi2(i) - gradtmp(1)*srcvals(10,i) - &
          gradtmp(2)*srcvals(11,i)-gradtmp(3)*srcvals(12,i)
      enddo
!$OMP END PARALLEL DO      
!
!  End of computing phi1 = S_{i|k|}[\sigma]
!  and phi2 = S_{i|k|}'[\sigma]
!

!
!  Now compoute S_{i|k|}'[\phi2] and add i\alpha S_{i|k|}'[phi2]
!  to pot
!
      call oversample_fun_surf(int(2,8),npatches,norders,ixyzs,iptype,& 
        npts,phi2,novers,ixyzso,ns,sigmaover)

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ns
        charges(i) = sigmaover(i)*whtsover(i)*over4pi 
      enddo
!$OMP END PARALLEL DO

      call hfmm3d_t_c_g(eps,ztmp,ns,sources,charges,npts, &
        srctmp,pot_aux,grad_aux,ier)

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        pot_aux(i) = grad_aux(1,i)*srcvals(10,i)+ &
          grad_aux(2,i)*srcvals(11,i)+&
          grad_aux(3,i)*srcvals(12,i)
      enddo
!$OMP END PARALLEL DO      
    

!
!       Add near field precomputed contribution
!

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart) &
!$OMP PRIVATE(jstart,npols,l)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            pot_aux(i) = pot_aux(i) + &
              wnear(2*nquad+jquadstart+l-1)*phi2(jstart+l-1)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

!
!     Remove near contribution of the FMM
!

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2) &
!$OMP PRIVATE(ctmp2,nss,l,jstart,pottmp,gradtmp)
      do i=1,npts
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
            nss = nss+1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)
            ctmp2(nss)= charges(l) 
          enddo
        enddo
        pottmp = 0
        gradtmp(1) = 0
        gradtmp(2) = 0
        gradtmp(3) = 0

        call h3ddirectcg(nd,ztmp,srctmp2,ctmp2,nss,srctmp(1,i), &
          ntarg0,pottmp,gradtmp,thresh)
        pot_aux(i) = pot_aux(i) - gradtmp(1)*srcvals(10,i) - &
          gradtmp(2)*srcvals(11,i)-gradtmp(3)*srcvals(12,i)
        pot(i) = pot(i)+ ima*zpars(2)*pot_aux(i)
      enddo
!$OMP END PARALLEL DO      

!
! End of adding i \alpha S_{ik}'^2[\sigma] to pot
!
!

!
!  Begin computation of D_{k}'[\phi1] 
!   we will not handle to subtraction of the near correction
!   until compouting D_{i|k|}'[\phi1] and take care of
!   the total subtraction together
!
!  the array phi2 is no longer neeeded, so we will reuse it
!  for temporary storage
!

      call oversample_fun_surf(int(2,8),npatches,norders,ixyzs,iptype,& 
        npts,phi1,novers,ixyzso,ns,sigmaover)

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ns
        dipvec(1,i) = sigmaover(i)*whtsover(i)*srcover(10,i)*over4pi 
        dipvec(2,i) = sigmaover(i)*whtsover(i)*srcover(11,i)*over4pi 
        dipvec(3,i) = sigmaover(i)*whtsover(i)*srcover(12,i)*over4pi 
      enddo
!$OMP END PARALLEL DO

      call hfmm3d_t_d_g(eps,zpars(1),ns,sources,dipvec,npts, &
        srctmp,pot_aux,grad_aux,ier)

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        phi2(i) = grad_aux(1,i)*srcvals(10,i) + &
          grad_aux(2,i)*srcvals(11,i) + grad_aux(3,i)*srcvals(12,i)
      enddo
!$OMP END PARALLEL DO
      

      call hfmm3d_t_d_g(eps,ztmp,ns,sources,dipvec,npts, &
        srctmp,pot_aux,grad_aux,ier)

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        pot_aux(i) = phi2(i) - grad_aux(1,i)*srcvals(10,i) - &
          grad_aux(2,i)*srcvals(11,i) - grad_aux(3,i)*srcvals(12,i)
      enddo
!$OMP END PARALLEL DO
      

!
!       Add near field precomputed contribution
!

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart) &
!$OMP PRIVATE(jstart,pottmp,npols,l)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            pot_aux(i) = pot_aux(i) + &
              wnear(3*nquad+jquadstart+l-1)*phi1(jstart+l-1)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

!
!     Remove near contribution of the FMM
!

      zpars2(1) = zpars(1)
      zpars2(2) = ztmp
      ndd = 0
      ndi = 0
      ndz = 2
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,l,jpatch,pottmp,rr)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
            pottmp = 0
            rr = sqrt((srcover(1,l)-srcvals(1,i))**2 + &
              (srcover(2,l)-srcvals(2,i))**2 + &
              (srcover(3,l)-srcvals(3,i))**2)
            if(rr.lt.thresh) goto 1311
            call h3d_dprime_diff(srcover(1,l),int(12,8),srcvals(1,i),ndd, &
              dpars,ndz,zpars2,ndi,ipars,pottmp)
            pot_aux(i) = pot_aux(i) - pottmp*sigmaover(l)*whtsover(l)
 1311       continue            
          enddo
        enddo
        pot(i) = pot(i)+ ima*zpars(2)*pot_aux(i)
      enddo
!$OMP END PARALLEL DO      

!
! End of adding i \alpha (D_{k}' - D_{ik}')S_{ik}'[\sigma] to pot



      return
      end subroutine lpcomp_helm_rpcomb_neu_addsub
!
!
!
!
!
!

      subroutine helm_rpcomb_neu_solver(npatches,norders,ixyzs, &
        iptype,npts,srccoefs,srcvals,eps,zpars,numit,ifinout, &
        rhs,eps_gmres,niter,errs,rres,soln,siksoln)
!
!
!  This subroutine solves the Helmholtz Neumann problem
!  on the exterior of an object where the potential
!  is represented as a right preconditioned 
!  combined field integral representation.
!
!
!  Representation:
!    u = S_{k}[\rho]+i*alpha*D_{k}[S_{i|k|}[\rho]]
!
!  Boundary condition:
!    u'=f
!
!  The linear system is solved iteratively using GMRES
!  until a relative residual of eps_gmres is reached
!
!
!  Input:
!
!    - npatches: integer(8)
!        number of patches
!    - norders: integer(8)(npatches)
!        order of discretization on each patch 
!    - ixyzs: integer(8)(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(8)(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!    - npts: integer(8)
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
!    - zpars: complex *16 (2)
!        kernel parameters (Referring to formula (1))
!          * zpars(1) = k 
!          * zpars(2) = alpha
!    - numit: integer(8)
!        max number of gmres iterations
!    - ifinout: integer(8)
!        ifinout = 0, interior problem
!        ifinout = 1, exterior problem
!      
!    - rhs: complex *16(npts)
!        right hand side
!    - eps_gmres: real *8
!        gmres tolerance requested
!      
!
!  output
!    - niter: integer(8)
!        number of gmres iterations required for relative residual
!        to converge to eps_gmres
!    - errs:  real *8 (numit)
!        relative residual as a function of iteration
!        number (only errs(1:niter) is populated))
!    - rres: real *8
!        relative residual for computed solution
!    - soln: complex *16(npts)
!        density which solves the neumann problem \rho
!    - siksoln: complx *16(npts)
!        sik[\rho] which can be used for far field
!        computations later
!				 
!
!
      implicit none
      integer(8) npatches,norder,npols,npts
      integer(8) norders(npatches),ixyzs(npatches+1)
      integer(8) iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,eps_gmres
      integer(8) ifinout
      complex *16 zpars(2)
      complex *16 rhs(npts)
      complex *16 soln(npts),siksoln(npts)

      real *8, allocatable :: targs(:,:)
      integer(8), allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      integer(8) ndtarg,ntarg

      real *8 errs(numit+1)
      real *8 rres,eps2
      integer(8) niter


      integer(8) nover,npolso,nptso
      integer(8) nnz,nquad
      integer(8), allocatable :: row_ptr(:),col_ind(:),iquad(:)

      complex *16, allocatable :: wnear(:)

      real *8, allocatable :: srcover(:,:),wover(:)
      integer(8), allocatable :: ixyzso(:),novers(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:) 

      integer(8) i,j,jpatch,jquadstart,jstart

      integer(8) ipars(2)
      real *8 dpars,timeinfo(10),t1,t2,omp_get_wtime


      real *8 ttot,done,pi
      real *8 rfac,rfac0
      integer(8) iptype_avg,norder_avg
      integer(8) ikerorder, iquadtype,npts_over
      integer(8) n_var

!
!
!       gmres variables
!
      complex *16 zid,ztmp
      real *8 rb,wnrm2
      integer(8) numit,it,iind,it1,k,l,count1
      real *8 rmyerr
      complex *16 temp
      complex *16, allocatable :: vmat(:,:),hmat(:,:)
      complex *16, allocatable :: cs(:),sn(:)
      complex *16, allocatable :: svec(:),yvec(:),wtmp(:)
      complex *16 ima

      ima=(0.0d0,1.0d0)

!
!   n_var is the number of unknowns in the linear system.
!   in this case n_var=npts
!

      n_var=npts


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

      ikerorder = -1
!
!    estimate oversampling for far-field, and oversample geometry
!

      allocate(novers(npatches),ixyzso(npatches+1))

      print *, "beginning far order estimation"

      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,&
       rads,npts,srccoefs,ndtarg,npts,targs,ikerorder,zpars(1),&
       nnz,row_ptr,col_ind,rfac,novers,ixyzso)

      npts_over = ixyzso(npatches+1)-1

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

      call getnearquad_helm_rpcomb_neu(npatches,norders,&
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
!      In this case it is a bit messy because of the identity that
!      comes from the calderon identity
  
      zid=((-1)**(ifinout)/2.0d0-ima*zpars(2)/4)


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

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,n_var
        vmat(i,1) = rhs(i)/rb
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

      call lpcomp_helm_rpcomb_neu_addsub(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,eps,zpars,nnz,row_ptr,col_ind, &
        iquad,nquad,wnear,vmat(1,it),novers,npts_over,ixyzso,srcover,&
        wover,wtmp,siksoln)

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

      call lpcomp_helm_rpcomb_neu_addsub(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,eps,zpars,nnz,row_ptr,col_ind, &
        iquad,nquad,wnear,soln,novers,npts_over,ixyzso,srcover,&
        wover,wtmp,siksoln)

!$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rres)
          do i=1,npts
            rres = rres + abs(zid*soln(i) + wtmp(i)-rhs(i))**2
          enddo
!$OMP END PARALLEL DO          
          rres = sqrt(rres)/rb
          niter = it
          return

        endif
      enddo
!
      return
      end subroutine helm_rpcomb_neu_solver
!
!
!
!        
      subroutine helm_rpcomb_neu_solver_memest(npatches,norders,ixyzs, &
        iptype,npts,srccoefs,srcvals,eps,zpars,numit,rmem)
!
!
!  This subroutine is the memory estimation routine for
!  the Helmholtz Neumann problem
!  solver where the potential
!  is represented as a right preconditioned 
!  combined field integral representation.
!
!
!  Representation:
!    u = S_{k}[\rho]+i*alpha*D_{k}[S_{i|k|}[\rho]]
!
!  Boundary condition:
!    u'=f
!
!  The linear system is solved iteratively using GMRES
!  until a relative residual of eps_gmres is reached
!
!
!  Input:
!
!    - npatches: integer(8)
!        number of patches
!    - norders: integer(8)(npatches)
!        order of discretization on each patch 
!    - ixyzs: integer(8)(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(8)(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!    - npts: integer(8)
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
!    - zpars: complex *16 (2)
!        kernel parameters (Referring to formula (1))
!          * zpars(1) = k 
!          * zpars(2) = alpha
!    - numit: integer(8)
!        max number of gmres iterations
!      
!
!  output
!    - rmem: double precision
!        estimated memory required by code in GB. Note that
!        this is meant to serve as an estimate only. 
!        The exact memory usage might be between (0.75,1.25)*rmem
!        The memory estimate may not be reliable for a
!        very small number of points
! 
!
!
      implicit none
      integer(8) npatches,norder,npols,npts
      integer(8) ifinout
      integer(8) norders(npatches),ixyzs(npatches+1)
      integer(8) iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps
      complex *16 zpars(3)
      real *8 rmem

      integer(8) lmem8,bigint
      real *8 rmemfmm
      real *8, allocatable :: targs(:,:)
      integer(8), allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      integer(8) ndtarg,ntarg

      real *8 rres,eps2
      integer(8) niter


      integer(8) nover,npolso,nptso
      integer(8) nnz,nquad
      integer(8), allocatable :: row_ptr(:),col_ind(:),iquad(:)
      complex *16, allocatable :: wnear(:)

      real *8, allocatable :: srcover(:,:),wover(:),sources(:,:)
      integer(8), allocatable :: ixyzso(:),novers(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:) 

      integer(8) i,j,jpatch,jquadstart,jstart

      integer(8) ipars,ifcharge,ifdipole,ifpgh,ifpghtarg
      real *8 dpars,timeinfo(10),t1,t2,omp_get_wtime


      real *8 ttot,done,pi
      real *8 rfac,rfac0
      integer(8) iptype_avg,norder_avg
      integer(8) ikerorder, iquadtype,npts_over,iper

!
!
!       gmres variables
!
      integer(8) numit,k,l
      complex *16 temp
      
      bigint = numit+1
      bigint = bigint*npts*2
      lmem8 = lmem8 + bigint



      lmem8 = lmem8 + numit*(numit+5)*2 + npts*2


      done = 1
      pi = atan(done)*4


!
!
!        setup targets as on surface discretization points
! 
      ndtarg = 3
      ntarg = npts
      allocate(targs(ndtarg,npts),uvs_targ(2,ntarg),ipatch_id(ntarg))
      lmem8 = lmem8 + ndtarg*npts + 3*ntarg 

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        targs(1,i) = srcvals(1,i)
        targs(2,i) = srcvals(2,i)
        targs(3,i) = srcvals(3,i)
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo
!$OMP END PARALLEL DO   


!
!    initialize patch_id and uv_targ for on surface targets
!
      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts, &
       ipatch_id,uvs_targ)

!
!
!        this might need fixing
!
      iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
      norder_avg = floor(sum(norders)/(npatches+0.0d0))

      call get_rfacs(norder_avg,iptype_avg,rfac,rfac0)


      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))
      lmem8 = lmem8 + 5*npatches

      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts, &
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
      
      call findnear(cms,npatches,rad_near,ndtarg,targs,npts,row_ptr, &
             col_ind)

      allocate(iquad(nnz+1)) 
      lmem8 = lmem8 + npts+1+ 2*nnz
      call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind, &
              iquad)

      ikerorder = 0


!
!    estimate oversampling for far-field, and oversample geometry
!

      allocate(novers(npatches),ixyzso(npatches+1))
      lmem8 = lmem8 + 2*npatches + 1

      print *, "beginning far order estimation"

      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms, &
         rads,npts,srccoefs,ndtarg,npts,targs,ikerorder,zpars(1), &
         nnz,row_ptr,col_ind,rfac,novers,ixyzso)

      npts_over = ixyzso(npatches+1)-1

      allocate(srcover(12,npts_over),wover(npts_over))
      lmem8 = lmem8 + 15*npts_over

      call oversample_geom(npatches,norders,ixyzs,iptype,npts, &
        srccoefs,srcvals,novers,ixyzso,npts_over,srcover)

      call get_qwts(npatches,novers,ixyzso,iptype,npts_over, &
             srcover,wover)


!
!   compute near quadrature correction
!
      nquad = iquad(nnz+1)-1
      allocate(wnear(nquad))
      lmem8 = lmem8 + nquad*2*4
      rmem = lmem8*8/1024/1024/1024

      ifcharge = 1
      ifdipole = 0
      if(abs(zpars(2)).gt.1.0d-16) ifdipole = 1

      ifpgh = 0
      ifpghtarg = 2
      allocate(sources(3,npts_over))
!$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,npts_over
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)
      enddo
!$OMP END PARALLEL DO      

      iper = 0
      rmemfmm = 0
      call hfmm3d_memest(int(1,8),eps,zpars(1),npts_over,sources,ifcharge, &
        ifdipole,iper,ifpgh,npts,targs,ifpghtarg,rmemfmm)
      rmem = rmem + rmemfmm
      
!
      return
      end
!
!
!
!
!
!
!
      subroutine lpcomp_helm_rpcomb_dir(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,ipatch_id, &
        uvs_targ,eps,zpars,sigma,sigma1,pot)
!
!
!  This subroutine evaluates the dirichlet data for
!  the right preconditioned
!  combined field representation.
!
!  Representation:
!    u = S_{k}[\sigma]+i*alpha*D_{k}[S_{i|k|}[\sigma]]
!
!
!  Input:
!
!    - npatches: integer(8)
!        number of patches
!    - norders: integer(8)(npatches)
!        order of discretization on each patch 
!    - ixyzs: integer(8)(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(8)(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!    - npts: integer(8)
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
!    - ndtarg: integer(8)
!        leading dimension of target array
!    - ntarg: integer(8)
!        number of targets
!    - targs: real *8 (ndtarg,ntarg)
!        target info, the first three coordinates must be
!        the xyz components of the target
!    - ipatch_id: integer(8)(ntarg)
!        patch on which target is on if on-surface, =-1/0 if 
!        target is off-surface
!    - uvs_targ: real *8 (2,ntarg)
!        local uv coordinates on patch, if target is on surface
!    - eps: real *8
!        precision requested for computing quadrature and fmm
!        tolerance
!    - zpars: complex *16 (2)
!        kernel parameters (Referring to formula (1))
!          * zpars(1) = k 
!          * zpars(2) = alpha
!    - sigma: complex *16(npts)
!        density sigma above
!    - sigma1: complex *16 (npts)
!        sik(sigma)
!
!  Output arguments:
!    - pot: complex *16(ntarg)
!        potential at the target locations
!				 
!
      implicit none
      integer(8), intent(in) :: npatches
      integer(8), intent(in) :: norders(npatches),ixyzs(npatches+1)
      integer(8), intent(in) :: iptype(npatches),npts
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
      integer(8), intent(in) :: ndtarg,ntarg
      real *8, intent(in) :: targs(ndtarg,ntarg)
      integer(8), intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg),eps
      complex *16, intent(in) :: zpars(2),sigma(npts),sigma1(npts)
      complex *16, intent(out) :: pot(ntarg)

      complex *16, allocatable :: pottmp2(:)
      complex *16 zpars_tmp(3),ima
      integer(8), allocatable :: ipatch_id_src(:)
      real *8, allocatable :: uvs_src(:,:)
      integer(8) i,ndtarg0

      data ima/(0.0d0,1.0d0)/


      allocate(pottmp2(ntarg))

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        pottmp2(i) = 0
        pot(i) = 0
      enddo
!$OMP END PARALLEL DO
      
      
!
!  compute S_{k} [sigma]
!
      zpars_tmp(1) = zpars(1)
      zpars_tmp(2) = 1.0d0 
      zpars_tmp(3) = 0

      
      call lpcomp_helm_comb_dir(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,ipatch_id, &
        uvs_targ,eps,zpars_tmp,sigma,pot)
      
      
!
!  compute D_{k} [pottmp]
!
      zpars_tmp(1) = zpars(1)
      zpars_tmp(2) = 0
      zpars_tmp(3) = 1.0d0

      
      call lpcomp_helm_comb_dir(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,ipatch_id, &
        uvs_targ,eps,zpars_tmp,sigma1,pottmp2)


!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        pot(i) = pot(i) + pottmp2(i)*zpars(2)*ima
      enddo
!$OMP END PARALLEL DO

     
      return
      end subroutine lpcomp_helm_rpcomb_dir

