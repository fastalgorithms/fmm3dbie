!
!  This file contains the following user callable routines
!  for solving Maxwell equations using the mfie integral
!  represenation:
!    getnearquad_em_mfie_pec: routine for generating the near
!     field quadrature correction  
!
!
!





subroutine getnearquad_em_mfie_pec(npatches,norders,&
  ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
  ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,&
  iquad,rfac0,nquad,wnear)
!
!  This subroutine generates the near field quadrature
!  for the MFIE integral equation:
!
!  J/2-M_{k}[J]= n \times H_inc
!
!  and the additional operator ikn·S_{k}[J] to get the charge therm
!  from n·E_inc in a separate solver if E is needed.
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
!    - wnear must be of size 6*nquad as 6 different layer
!      potentials are returned
!      * the uu component of the matrix kernel -M_{k}[J]
!      * the uv component of the matrix kernel -M_{k}[J]
!      * the vu component of the matrix kernel -M_{k}[J]
!      * the vv component of the matrix kernel -M_{k}[J]
!      * the u component of the matrix kernel ikn·S_{k}[J]
!      * the v component of the matrix kernel ikn·S_{k}[J]
!
!
!  Input arguments:
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
!    - zpars: complex *16 (2)
!        kernel parameters (Referring to formula (1))
!        zpars(1) = k 
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
!    - wnear: complex *16(6*nquad)
!        The desired near field quadrature
!              

      implicit none 
      integer, intent(in) :: npatches,norders(npatches),npts,nquad
      integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),eps
      real *8, intent(in) :: rfac0
      integer, intent(in) :: ndtarg,ntarg
      integer, intent(in) :: iquadtype
      real *8, intent(in) :: targs(ndtarg,ntarg)
      complex *16, intent(in) :: zpars
      integer, intent(in) :: nnz
      integer, intent(in) :: row_ptr(ntarg+1),col_ind(nnz),iquad(nnz+1)
      integer, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      complex *16, intent(out) :: wnear(6*nquad)
      integer ipars(2)
      real *8 dpars(1)

      complex *16 alpha,beta
      integer i,j,ndi,ndd,ndz

      integer ipv

      procedure (), pointer :: fker
      external  fker_em_mfie_pec

     ndz=3
     ndd=1
     ndi=2
     ipv=0

     fker =>  fker_em_mfie_pec
     ipars(1)=1
     ipars(2)=1
     call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
     &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(1:nquad))

      ipars(1)=1
      ipars(2)=2
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,&
     &nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(nquad+1:2*nquad))
 
      ipars(1)=2
      ipars(2)=1
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,&
     &nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(2*nquad+1:3*nquad))

      ipars(1)=2
      ipars(2)=2
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,&
     &nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(3*nquad+1:4*nquad))

      ipars(1)=3
      ipars(2)=1
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,&
     &nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(4*nquad+1:5*nquad))

      ipars(1)=3
      ipars(2)=2
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,&
     &nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(5*nquad+1:6*nquad))

      return
      end subroutine getnearquad_em_mfie_pec


      subroutine lpcomp_em_mfie_pec_addsub(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,sigma,novers,&
     &nptso,ixyzso,srcover,whtsover,pot,wnear)
!
!
!  this subroutine evaluates the layer potential for -M{k}[J] the MFIE
!  operator
!  the boundary integral equation is (MFIE):
!
!			J/2-M_{k}[J] = nxH_inc 
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
!    - zpars: complex *16 (2)
!        kernel parameters (Referring to formula (1))
!        zpars(1) = k 
!        zpars(2) = alpha
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
!    - wnear: complex *16(6*nquad)
!        Precomputed near field quadrature
!      * the uu component of the matrix kernel -M_{k}[J]
!      * the uv component of the matrix kernel -M_{k}[J]
!      * the vu component of the matrix kernel -M_{k}[J]
!      * the vv component of the matrix kernel -M_{k}[J]
!      * the u component of the matrix kernel ikn·S_{k}[J]
!      * the v component of the matrix kernel ikn·S_{k}[J]
!    - sigma - complex *16(2*ns)
!        induced physical current on the surface
!        sigma(1:ns) - first component of  J along
!          the srcvals(4:6,i) direction
!        sigma(ns+1:2*ns) - second component of J along
!          the (srcvals(10:12,i) x srcvals(4:6,i)) direction
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
!    - pot(1:ntarg) - complex *16
!        first component of -M_{k}[J] along
!        the srcvals(4:6,i) direction
!    - pot(ntarg+1,2*ntarg) - complex *16
!        second component of -M_{k}[J] along
!        the (srcvals(10:12,i) x srcvals(4:6,i)) direction

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
      complex *16 sigma(2*npts)
      complex *16 pot(2*ntarg)
	  
	  complex *16 wnear(6*nquad)

      integer novers(npatches+1)
      integer nover,npolso,nptso
      real *8 srcover(12,nptso),whtsover(nptso)
      complex *16, allocatable :: potsort(:)

      real *8, allocatable :: sources(:,:),targvals(:,:),wtmp2(:)
      complex *16, allocatable :: charges(:),dipvec(:,:),sigmaover(:)
      integer ns,nt
      complex *16 alpha,beta  
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      complex *16 tmp(10),val,E(2)

      real *8 xmin,xmax,ymin,ymax,zmin,zmax,sizey,sizez,boxsize


      integer i,j,jpatch,jquadstart,jstart


      integer ifaddsub,ifdir

      integer ntj
      
      complex *16 zdotu,pottmp
      complex *16, allocatable :: ctmp2_u(:),ctmp2_v(:),dtmp2(:,:)
      real *8 radexp,epsfmm

      integer ipars(2)
      real *8 dpars(1),timeinfo(10),t1,t2,omp_get_wtime

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
      allocate(sigmaover(3*ns))
! 
!       oversample density
!

      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
     &npts,sigma(1:npts),novers,ixyzso,ns,sigmaover(1:ns))
      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
     &npts,sigma(npts+1:2*npts),novers,ixyzso,ns,sigmaover(ns+1:2*ns))


      ra = 0

!
!        compute threshold for ignoring local computation
!
      call get_fmm_thresh(12,ns,srcover,12,npts,srcvals,thresh)

!
!       fmm call
!
      ifdir=0
	  call em_mfie_pec_FMM(eps,zpars,ns,npts,srcover,targs,whtsover,&
	  &sigmaover(1:ns),sigmaover(ns+1:2*ns),pot(1:npts)&
	  &,pot(npts+1:2*npts),thresh,ifdir)

      
!
!       Add near field precomputed contribution
!

      call cpu_time(t1)
!C$      t1 = omp_get_wtime()

!C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart)
!C$OMP$PRIVATE(jstart,pottmp,npols)
      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            pot(i) = pot(i) + wnear(jquadstart+l-1)*sigma(jstart+l-1)
            pot(i) = pot(i) + wnear(nquad+jquadstart+l-1)*&
			&sigma(jstart+l-1+npts)

            pot(i+npts) = pot(i+npts) + wnear(2*nquad+jquadstart+l-1)&
			&*sigma(jstart+l-1)
            pot(i+npts) = pot(i+npts) + wnear(3*nquad+jquadstart+l-1)&
			&*sigma(jstart+l-1+npts)
          enddo
        enddo
      enddo
!C$OMP END PARALLEL DO

!
!     Remove near contribution of the FMM
!

!C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
!C$OMP$PRIVATE(ctmp2_u,ctmp2_v,wtmp2,nss,l,jstart,ii,E,npover)
      ifdir=1
      do i=1,ntarg
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          nss = nss + ixyzso(jpatch+1)-ixyzso(jpatch)
        enddo
        allocate(srctmp2(12,nss),ctmp2_u(nss),ctmp2_v(nss),wtmp2(nss))

        rmin = 1.0d6
        ii = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          jstart = ixyzso(jpatch)-1
          npover = ixyzso(jpatch+1)-ixyzso(jpatch)
          do l=1,npover
			ii = ii+1
			srctmp2(:,ii) = srcover(:,jstart+l)            
			ctmp2_u(ii)=sigmaover(jstart+l)
			ctmp2_v(ii)=sigmaover(jstart+l+ns)
			wtmp2(ii)=whtsover(jstart+l)
          enddo
        enddo
	  call em_mfie_pec_FMM(eps,zpars,nss,ntarg0,srctmp2,targs(:,i),wtmp2,&
	  &ctmp2_u,ctmp2_v,E(1),E(2),thresh,ifdir)
	    pot(i) = pot(i) - E(1)
	    pot(i+ntarg) = pot(i+ntarg) - E(2)

        deallocate(srctmp2,ctmp2_u,ctmp2_v,wtmp2)
      enddo
      
      call cpu_time(t2)
!C$      t2 = omp_get_wtime()     

      timeinfo(2) = t2-t1


!!      call prin2('quadrature time=*',timeinfo,2)
      
      ttot = timeinfo(1) + timeinfo(2)
!!      call prin2('time in lpcomp=*',ttot,1)

      return
      end subroutine lpcomp_em_mfie_pec_addsub


      subroutine lpcomp_em_mfie_pec_addsub_last(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,sigma,novers,&
     &nptso,ixyzso,srcover,whtsover,pot,wnear,rhs_nE)
!
!  Same as lpcomp_em_mfie_pec_addsub but also provides 
!  rhs_nE:=ikn·S_{k}[J]
!  This is only used for the last iteration of GMRES.
!  The value of rhs_nE is used later outside this subroutine
!  and n·E^{inc} is added to get the RHS of the auxiliar 
!  scalar equation corresponding to the charge term
!
!
!  this subroutine evaluates the layer potential for -M_{k}[J] the MFIE
!  operator
!  the boundary integral equation is (MFIE):
!
!			J/2-M_{k}[J] = nxH_inc 
!
!  and also computes:
!
!          ikn·S_{k}[J]
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
!    - zpars: complex *16 (2)
!        kernel parameters (Referring to formula (1))
!        zpars(1) = k 
!        zpars(2) = alpha
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
!    - wnear: complex *16(6*nquad)
!        Precomputed near field quadrature
!      * the uu component of the matrix kernel -M_{k}[J]
!      * the uv component of the matrix kernel -M_{k}[J]
!      * the vu component of the matrix kernel -M_{k}[J]
!      * the vv component of the matrix kernel -M_{k}[J]
!      * the u component of the matrix kernel ikn·S_{k}[J]
!      * the v component of the matrix kernel ikn·S_{k}[J]
!    - sigma - complex *16(2*ns)
!        induced physical current on the surface
!        sigma(1:ns) - first component of  J along
!          the srcvals(4:6,i) direction
!        sigma(ns+1:2*ns) - second component of J along
!          the (srcvals(10:12,i) x srcvals(4:6,i)) direction
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
!    - pot(1:ntarg) - complex *16
!        first component of -M_{k}[J] along
!        the srcvals(4:6,i) direction
!    - pot(ntarg+1,2*ntarg) - complex *16
!        second component of -M_{k}[J] along
!        the (srcvals(10:12,i) x srcvals(4:6,i)) direction
!    - rhs_nE(1:ntarg) - complex *16
!        value of ikn·S_{k}[J]
!

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
      complex *16 sigma(2*npts)
      complex *16 pot(2*ntarg)
	  complex *16 rhs_nE(ntarg)
	  
	  complex *16 wnear(6*nquad)

      integer novers(npatches+1)
      integer nover,npolso,nptso
      real *8 srcover(12,nptso),whtsover(nptso)
      complex *16, allocatable :: potsort(:)

      real *8, allocatable :: sources(:,:),targvals(:,:),wtmp2(:)
      complex *16, allocatable :: charges(:),dipvec(:,:),sigmaover(:)
      integer ns,nt
      complex *16 alpha,beta
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      complex *16 tmp(10),val,E(3)

      real *8 xmin,xmax,ymin,ymax,zmin,zmax,sizey,sizez,boxsize


      integer i,j,jpatch,jquadstart,jstart,ifdir


      integer ifaddsub

      integer ntj
      
      complex *16 zdotu,pottmp
      complex *16, allocatable :: ctmp2_u(:),ctmp2_v(:),dtmp2(:,:)
      real *8 radexp,epsfmm

      integer ipars(2)
      real *8 dpars(1),timeinfo(10),t1,t2,omp_get_wtime

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
      allocate(sigmaover(3*ns))

! 
!       oversample density
!

      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
     &npts,sigma(1:npts),novers,ixyzso,ns,sigmaover(1:ns))
      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
     &npts,sigma(npts+1:2*npts),novers,ixyzso,ns,sigmaover(ns+1:2*ns))


      ra = 0

!
!        compute threshold for ignoring local computation
!

      call get_fmm_thresh(12,ns,srcover,12,npts,srcvals,thresh)

!
!       fmm call
!
	 ifdir=0
	  call em_mfie_pec_FMM_last(eps,zpars,ns,npts,srcover,targs,&
	 &whtsover,sigmaover(1:ns),sigmaover(ns+1:2*ns),pot(1:npts),&
	 &pot(npts+1:2*npts),rhs_nE,thresh,ifdir)

      
!
!       Add near field precomputed contribution
!



      call cpu_time(t1)
!C$      t1 = omp_get_wtime()

!C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart)
!C$OMP$PRIVATE(jstart,pottmp,npols)
      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            pot(i) = pot(i) + wnear(jquadstart+l-1)*sigma(jstart+l-1)
            pot(i) = pot(i) + wnear(nquad+jquadstart+l-1)*&
			&sigma(jstart+l-1+npts)

            pot(i+npts) = pot(i+npts) + wnear(2*nquad+jquadstart+l-1)&
			&*sigma(jstart+l-1)
            pot(i+npts) = pot(i+npts) + wnear(3*nquad+jquadstart+l-1)&
			&*sigma(jstart+l-1+npts)
			
            rhs_nE(i) = rhs_nE(i) + wnear(4*nquad+jquadstart+l-1)&
			&*sigma(jstart+l-1)
            rhs_nE(i) = rhs_nE(i) + wnear(5*nquad+jquadstart+l-1)&
			
			&*sigma(jstart+l-1+npts)			
			
          enddo
        enddo
      enddo
!C$OMP END PARALLEL DO

!
!     Remove near contribution of the FMM
!

!C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
!C$OMP$PRIVATE(ctmp2_u,ctmp2_v,wtmp2,nss,l,jstart,ii,E,npover)
	ifdir=1
      do i=1,ntarg
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          nss = nss + ixyzso(jpatch+1)-ixyzso(jpatch)
        enddo
        allocate(srctmp2(12,nss),ctmp2_u(nss),ctmp2_v(nss),wtmp2(nss))

        rmin = 1.0d6
        ii = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          jstart = ixyzso(jpatch)-1
          npover = ixyzso(jpatch+1)-ixyzso(jpatch)
          do l=1,npover
			ii = ii+1
			srctmp2(:,ii) = srcover(:,jstart+l)            
			ctmp2_u(ii)=sigmaover(jstart+l)
			ctmp2_v(ii)=sigmaover(jstart+l+ns)
			wtmp2(ii)=whtsover(jstart+l)
          enddo
        enddo
		call em_mfie_pec_FMM_last(eps,zpars,nss,ntarg0,srctmp2,targs(:,i),&
	 &wtmp2,ctmp2_u,ctmp2_v,E(1),&
	 &E(2),E(3),thresh,ifdir)

!		call fker_em_mfie_pec_last(nss,ntarg0,srctmp2, targs(:,i),&
!	   & dpars,zpars,ipars,ctmp2_u,ctmp2_v,wtmp2,E,thresh)

		!!sigmaover(i)*whtsover(i)
	    pot(i) = pot(i) - E(1)
	    pot(i+ntarg) = pot(i+ntarg) - E(2)
	    rhs_nE(i) = rhs_nE(i) - E(3)

        deallocate(srctmp2,ctmp2_u,ctmp2_v,wtmp2)
      enddo
      
      call cpu_time(t2)
!C$      t2 = omp_get_wtime()     

      timeinfo(2) = t2-t1


!!      call prin2('quadrature time=*',timeinfo,2)
      
      ttot = timeinfo(1) + timeinfo(2)
!!      call prin2('time in lpcomp=*',ttot,1)

      return
      end subroutine lpcomp_em_mfie_pec_addsub_last


      subroutine em_mfie_solver(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,eps,zpars,numit,ifinout,&
     &rhs,eps_gmres,niter,errs,rres,soln,rhs_nE)
!
!
!		This subroutine solves the Scattering Maxwell p.e.c. problem.
!		The the equations are:
!
!			curlE=ikH; curlH =-ikE
!
!		Representation:
!
!			H=curlS_{k}[J]
!
!		Boundary conditions:
!
!			nxH = J
!
!		Boundary integral equation:
!
!			J/2-M_{k}[J] = nxH_inc
!
!     This subroutine also outputs the value rhs_nE=ikn·S_{k}[J]
!
!     The linear system is solved iteratively using GMRES
!     until a relative residual of 1e-15 is reached
!
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
!    - zpars: complex *16 (2)
!        kernel parameters (Referring to formula (1))
!          * zpars(1) = k 
!          * zpars(2) = alpha
!    - numit: integer
!        max number of gmres iterations
!    - ifinout: integer
!        ifinout = 0, interior problem
!        ifinout = 1, exterior problem
!      
!    - rhs: complex *16(2*npts)
!        right hand side
!         rhs(1:npts)
!        rhs(1:ns) - first component of nxH_inc along
!          the srcvals(4:6,i) direction (u direction)
!        rhs(ns+1:2*ns) - second component of nxH_inc along
!          the (srcvals(10:12,i) x srcvals(4:6,i)) direction (v direct.)
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
!        soln(1:ns) - first component of J along
!          the srcvals(4:6,i) direction (u direction)
!        soln(ns+1:2*ns) - second component of J along
!          the (srcvals(10:12,i) x srcvals(4:6,i)) direction (v direct.)
!    - rhs_nE: complx *16(npts)
!        ikn·S_{k}[J] which can be used for computing the charge term
!        from n·E_inc using the em_aumfie_pec solver
!				 

      implicit none
      integer npatches,norder,npols,npts
      integer ifinout
      integer norders(npatches),ixyzs(npatches+1)
      integer iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,eps_gmres
      complex *16 zpars(3)
      complex *16 rhs(2*npts)
      complex *16 soln(2*npts)
      complex *16 rhs_nE(npts)

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


!
!   n_var is the number of unknowns in the linear system.
!   as we have one vector unknown: J , 
!   we need n_var=2*npts
!

      n_var=2*npts

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
!C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
		targs(:,i)=srcvals(:,i)
		ipatch_id(i) = -1
		uvs_targ(1,i) = 0
		uvs_targ(2,i) = 0
      enddo
!C$OMP END PARALLEL DO   


!
!    initialize patch_id and uv_targ for on surface targets
!
      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,&
     &ipatch_id,uvs_targ)

!
!
!        this might need fixing
!
      iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
      norder_avg = floor(sum(norders)/(npatches+0.0d0))

      call get_rfacs(norder_avg,iptype_avg,rfac,rfac0)


      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))

      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts,& 
     &srccoefs,cms,rads)

!C$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,npatches
        rad_near(i) = rads(i)*rfac
      enddo
!C$OMP END PARALLEL DO      

!
!    find near quadrature correction interactions
!
      print *, "entering find near mem"
      call findnearmem(cms,npatches,rad_near,ndtarg,targs,npts,nnz)
      print *, "nnz=",nnz

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,ndtarg,targs,npts,row_ptr,&
     &col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind,&
     &iquad)

      ikerorder = -1
      if(abs(zpars(3)).gt.1.0d-16) ikerorder = 0


!
!    estimate oversampling for far-field, and oversample geometry
!

      allocate(novers(npatches),ixyzso(npatches+1))

      print *, "beginning far order estimation"

      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,&
     &rads,npts,srccoefs,ndtarg,npts,targs,ikerorder,zpars(1),&
     &nnz,row_ptr,col_ind,rfac,novers,ixyzso)

      npts_over = ixyzso(npatches+1)-1
      print *, "npts_over=",npts_over

      allocate(srcover(12,npts_over),wover(npts_over))

      call oversample_geom(npatches,norders,ixyzs,iptype,npts,&
     &srccoefs,srcvals,novers,ixyzso,npts_over,srcover)

      call get_qwts(npatches,novers,ixyzso,iptype,npts_over,&
     &srcover,wover)


!
!   compute near quadrature correction
!
      nquad = iquad(nnz+1)-1
      print *, "nquad=",nquad
      allocate(wnear(6*nquad))
      
!C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,6*nquad
		wnear(i)=0.0d0
	  enddo
!C$OMP END PARALLEL DO    


      iquadtype = 1

!!      eps2 = 1.0d-8

      print *, "starting to generate near quadrature"
      call cpu_time(t1)
!C$      t1 = omp_get_wtime()      

      call getnearquad_em_mfie_pec(npatches,norders,&
     &ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,&
     &ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,&
     &iquad,rfac0,nquad,wnear)

      call cpu_time(t2)
!C$      t2 = omp_get_wtime()     

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
!      zid = 0.5
      zid=0.5d0


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
      do i=1,n_var
        rb = rb + abs(rhs(i))**2
      enddo
      rb = sqrt(rb)

      do i=1,n_var
        vmat(i,1) = rhs(i)/rb
      enddo

      svec(1) = rb

      do it=1,numit
        it1 = it + 1

!
!        NOTE:
!        replace this routine by appropriate layer potential
!        evaluation routine  
!

        call lpcomp_em_mfie_pec_addsub(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,&
     &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,&
     &vmat(1,it),novers,npts_over,ixyzso,srcover,wover,wtmp,&
	 &wnear)

        do k=1,it
          hmat(k,it) = 0
          do j=1,n_var      
            hmat(k,it) = hmat(k,it) + wtmp(j)*conjg(vmat(j,k))
          enddo

          do j=1,n_var
            wtmp(j) = wtmp(j)-hmat(k,it)*vmat(j,k)
          enddo
        enddo
          
        hmat(it,it) = hmat(it,it)+zid
        wnrm2 = 0
        do j=1,n_var
          wnrm2 = wnrm2 + abs(wtmp(j))**2
        enddo
        wnrm2 = sqrt(wnrm2)

        do j=1,n_var
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
          do j=1,n_var
            soln(j) = 0
            do i=1,it
              soln(j) = soln(j) + yvec(i)*vmat(j,i)
            enddo
          enddo


          rres = 0
          do i=1,n_var
            wtmp(i) = 0
          enddo
!
!        NOTE:
!        replace this routine by appropriate layer potential
!        evaluation routine  
!



          call lpcomp_em_mfie_pec_addsub_last(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,&
     &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,&
     &soln,novers,npts_over,ixyzso,srcover,wover,wtmp,&
	 &wnear,rhs_nE)

            
          do i=1,npts
            rres = rres + abs(zid*soln(i) + wtmp(i)-rhs(i))**2
          enddo
          rres = sqrt(rres)/rb
          niter = it
          return
        endif
      enddo
!
      return
      end subroutine em_mfie_solver










subroutine fker_em_mfie_pec(srcinfo,ndt,targinfo,ndd,dpars,ndz,zpars,&
 &ndi,ipars,E_val)
implicit none

! this function provides the near field kernel that will use 
! zgetnearquad_ggq_guru through getnearquad_em_mfie_pec

    !List of calling arguments
	integer, intent(in) :: ndt,ndd,ndz,ndi
	real ( kind = 8 ), intent(in) :: srcinfo(12)
	real ( kind = 8 ), intent(in) :: targinfo(ndt)
	integer, intent(in) :: ipars(ndi)
	real ( kind = 8 ), intent(in) :: dpars(ndd)
	complex ( kind = 8 ), intent(in) :: zpars(ndz)
	complex ( kind = 8 ), intent(out) :: E_val

	
	!List of local variables
	complex ( kind = 8 ) E_mat(3,3)
	real ( kind = 8 ) ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ) ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ) du(3), dv(3),sour(3),targ(3)
	real ( kind = 8 ) r, dr(3),aux_real
	
	complex ( kind = 8 ) nxcurlSka(2,2),nxnxSkb(2,2)
	complex ( kind = 8 ) ngradSklambda,nSkb(1,2)

	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3),xprod_aux3(3),xprod_aux4(3)
	complex ( kind = 8 ) R1, ima,my_exp, zk,alpha
	real ( kind = 8 ) pi
	
	pi=3.1415926535897932384626433832795028841971d0
	ima=(0.0d0,1.0d0)
	zk=zpars(1)
	alpha=zpars(2)
	    
	sour(1)=srcinfo(1)
	sour(2)=srcinfo(2)
	sour(3)=srcinfo(3)
	
	n_s(1)=srcinfo(10)
	n_s(2)=srcinfo(11)
	n_s(3)=srcinfo(12)	

	targ(1)=targinfo(1)
	targ(2)=targinfo(2)
	targ(3)=targinfo(3)

	n_t(1)=targinfo(10)
	n_t(2)=targinfo(11)
	n_t(3)=targinfo(12)

	dr(1)=targ(1)-sour(1)
	dr(2)=targ(2)-sour(2)
	dr(3)=targ(3)-sour(3)
	r=sqrt((dr(1))**2+(dr(2))**2+(dr(3))**2)
	
	  R1=(ima*zk*r-1.0d0)/r**3*exp(ima*zk*r)/(4.0d0*pi)
	  my_exp=exp(ima*zk*r)/(4.0d0*pi)
		
	  call orthonormalize(srcinfo(4:6),n_s,ru_s,rv_s)
	  call orthonormalize(targinfo(4:6),n_t,ru_t,rv_t)

	  call get_nxcurlSka(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r,&
	   &nxcurlSka)
	  E_mat(1,1)=-nxcurlSka(1,1)
	  E_mat(1,2)=-nxcurlSka(1,2)
	  E_mat(2,1)=-nxcurlSka(2,1)
	  E_mat(2,2)=-nxcurlSka(2,2)
				
	  call get_nSkb(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,my_exp,r,nSkb)
	  E_mat(3,1)=+ima*zk*nSkb(1,1)
	  E_mat(3,2)=+ima*zk*nSkb(1,2)	
	
	E_val=E_mat(ipars(1),ipars(2))

return
end subroutine fker_em_mfie_pec


subroutine em_mfie_pec_FMM(eps,zpars,ns,nt,srcvals,targvals,wts,a_u,a_v,&
 &AA_u,AA_v,thresh,ifdir)
implicit none

!
!  This subroutine computes the far field contribution if the MFIE operator
!  via FMM
!
!  Boundary integral equation:
!
!	 J/2-M_{k}[J] = nxH_inc  ->  (AA_u,AA_v)
!
!  input:
!    eps - real * 8
!      epsilon for the fmm call
!
!    zpars - complex *16(3)
!      kernel parameters
!      zpars(1) = k 
!      zpars(2) = - (not used)	
!      zpars(3) = - (not used)	
!
!    ns - integer
!      number of sources (this is the oversampled set of sources)
!   
!    nt - integer 
!      number of targets (this is the not oversampled set of targets)
!
!    srcvals - real *8(12,ns) (oversampled surface)
!      xyz(u,v) and derivative info sampled at the 
!      discretization nodes on the surface
!      srcvals(1:3,i) - xyz info
!      srcvals(4:6,i) - dxyz/du info
!      srcvals(7:9,i) - dxyz/dv info
!
!    targvals - real *8(3,nt) (not oversampled surface)
!      xyz(u,v) and derivative info sampled at the 
!      discretization nodes on the surface
!      srcvals(1:3,i) - xyz info
!      srcvals(4:6,i) - dxyz/du info
!      srcvals(7:9,i) - dxyz/dv info
!
!    wts - real *8(ns)
!      smooth quadrature weights at oversampled nodes
!
!    a_u,a_v - complex *16(ns)
!      two components of the tangent induced current J on the surface
!      along srcvals(4:6,i) and (srcvals(10:12,i) x srcvals(4:6,i)) 
!      directions
!
!    thresh - real *8
!      threshold to remove the selfo interaction term
!
!    ifdir - integer
!      flag, ifdir=1 direct calculation N^2 (used to remove teh near terms)
!            ifdir=0 FMM activated
!
!  output:
!    AA_u - complex *16(nt)
!      first component of -M_{k}[J] along the srcvals(4:6,i) direction
!        
!    AA_v - complex *16(nt)
!      second component of -M_{k}[J] along
!      the (srcvals(10:12,i) x srcvals(4:6,i)) direction
!            


    !List of calling arguments
	real ( kind = 8 ), intent(in) :: eps
	complex ( kind = 8 ), intent(in) :: zpars(3)
    integer, intent(in) :: ns,nt
	real ( kind = 8 ), intent(in) :: srcvals(12,ns),wts(ns),targvals(12,nt)
    complex ( kind = 8 ), intent(in) :: a_u(ns),a_v(ns)
    complex ( kind = 8 ), intent(out) :: AA_u(nt),AA_v(nt)
	real ( kind = 8 ), intent(in) :: thresh
	integer, intent(in) :: ifdir 


    !List of local variables
	real ( kind = 8 ), allocatable :: n_vect_s(:,:),u_vect_s(:,:)
	real ( kind = 8 ), allocatable :: v_vect_s(:,:),source(:,:)

	real ( kind = 8 ), allocatable :: n_vect_t(:,:),u_vect_t(:,:)
	real ( kind = 8 ), allocatable :: v_vect_t(:,:),targets(:,:)

    complex ( kind = 8 ), allocatable :: a_vect(:,:),b_vect(:,:),lambda(:)
    complex ( kind = 8 ), allocatable :: E(:,:),curlE(:,:),divE(:),rho(:)
	complex ( kind = 8 ) ima,zk,alpha

    integer count1,count2
    integer ifa_vect,ifb_vect,iflambda,ifrho,ifE,ifcurlE,ifdivE
	real ( kind = 8 ) pi
	pi=3.1415926535897932384626433832795028841971d0

	ima=(0.0d0,1.0d0)
	zk=zpars(1)
!	alpha=zpars(2)

    allocate(a_vect(3,ns))
    allocate(b_vect(3,ns))
    allocate(lambda(ns))
	allocate(rho(ns))
    allocate(E(3,nt))
    allocate(curlE(3,nt))
    allocate(divE(nt))
	allocate(n_vect_s(3,ns))
	allocate(n_vect_t(3,nt))
	allocate(u_vect_s(3,ns))
	allocate(v_vect_s(3,ns))
	allocate(u_vect_t(3,nt))
	allocate(v_vect_t(3,nt))
	allocate(source(3,ns))
	allocate(targets(3,nt))

	do count1=1,ns
	  n_vect_s(:,count1)=srcvals(10:12,count1)
	  source(:,count1)=srcvals(1:3,count1)
	enddo
	call orthonormalize_all(srcvals(4:6,:),srcvals(10:12,:),&
	 &u_vect_s,v_vect_s,ns)
	
	do count1=1,nt
	  n_vect_t(:,count1)=targvals(10:12,count1)
	  targets(:,count1)=targvals(1:3,count1)
	enddo
	call orthonormalize_all(targvals(4:6,:),targvals(10:12,:),&
     &u_vect_t,v_vect_t,nt)
		
	do count1=1,ns
	  b_vect(1,count1)=(a_u(count1)*u_vect_s(1,count1)+&
	   &a_v(count1)*v_vect_s(1,count1))
	  b_vect(2,count1)=(a_u(count1)*u_vect_s(2,count1)+&
       &a_v(count1)*v_vect_s(2,count1))
	  b_vect(3,count1)=(a_u(count1)*u_vect_s(3,count1)+&
	   &a_v(count1)*v_vect_s(3,count1))
    enddo

    !Computing the full operator
    ifa_vect=0
    ifb_vect=1
    iflambda=0
    ifrho=0
    ifE=0
    ifcurlE=1
    ifdivE=0

	call Vector_Helmholtz_targ2(eps,zk,ns,source,wts,ifa_vect,a_vect,&
	 &ifb_vect,b_vect,iflambda,lambda,ifrho,rho,n_vect_s,ifE,E,ifcurlE&
	 &,curlE,ifdivE,divE,nt,targets,thresh,ifdir)
    do count1=1,nt
	  b_vect(1,count1)=n_vect_t(2,count1)*curlE(3,count1)-&
	   &n_vect_t(3,count1)*curlE(2,count1)
	  b_vect(2,count1)=n_vect_t(3,count1)*curlE(1,count1)-&
	   &n_vect_t(1,count1)*curlE(3,count1)
	  b_vect(3,count1)=n_vect_t(1,count1)*curlE(2,count1)-&
	   &n_vect_t(2,count1)*curlE(1,count1)
    enddo

    do count1=1,nt
	  AA_u(count1)=-(b_vect(1,count1)*u_vect_t(1,count1)+b_vect(2,count1)&
	   &*u_vect_t(2,count1)+b_vect(3,count1)*u_vect_t(3,count1))
	  AA_v(count1)=-(b_vect(1,count1)*v_vect_t(1,count1)+b_vect(2,count1)&
	   &*v_vect_t(2,count1)+b_vect(3,count1)*v_vect_t(3,count1))
    enddo

	deallocate(a_vect)
	deallocate(b_vect)
	deallocate(lambda)
	deallocate(rho)
	deallocate(E)
	deallocate(curlE)
	deallocate(divE)
	
	deallocate(u_vect_s)
	deallocate(v_vect_s)
	deallocate(n_vect_s)
	deallocate(source)

	deallocate(u_vect_t)
	deallocate(v_vect_t)
	deallocate(n_vect_t)
	deallocate(targets)

return
end subroutine em_mfie_pec_FMM


subroutine em_mfie_pec_FMM_last(eps,zpars,ns,nt,srcvals,targvals,wts,&
 &a_u,a_v,AA_u,AA_v,PHI,thresh,ifdir)
implicit none

!
!  This subroutine computes the far field contribution if the MFIE operator
!  via FMM
!
!  Boundary integral equation:
!
!	 J/2-M_{k}[J] = nxH_inc  ->  (AA_u,AA_v)
!
!  It also provides the value of PHI := ikn·S_{k}[J] which will be used
!  later in a subsequent solver to find tha charge density
!
!  input:
!    eps - real * 8
!      epsilon for the fmm call
!
!    zpars - complex *16(3)
!      kernel parameters
!      zpars(1) = k 
!      zpars(2) = - (not used)	
!      zpars(3) = - (not used)	
!
!    ns - integer
!      number of sources (this is the oversampled set of sources)
!   
!    nt - integer 
!      number of targets (this is the not oversampled set of targets)
!
!    srcvals - real *8(12,ns) (oversampled surface)
!      xyz(u,v) and derivative info sampled at the 
!      discretization nodes on the surface
!      srcvals(1:3,i) - xyz info
!      srcvals(4:6,i) - dxyz/du info
!      srcvals(7:9,i) - dxyz/dv info
!
!    targvals - real *8(3,nt) (not oversampled surface)
!      xyz(u,v) and derivative info sampled at the 
!      discretization nodes on the surface
!      srcvals(1:3,i) - xyz info
!      srcvals(4:6,i) - dxyz/du info
!      srcvals(7:9,i) - dxyz/dv info
!
!    wts - real *8(ns)
!      smooth quadrature weights at oversampled nodes
!
!    a_u,a_v - complex *16(ns)
!      two components of the tangent induced current J on the surface
!      along srcvals(4:6,i) and (srcvals(10:12,i) x srcvals(4:6,i)) 
!      directions
!
!    thresh - real *8
!      threshold to remove the selfo interaction term
!
!    ifdir - integer
!      flag, ifdir=1 direct calculation N^2 (used to remove teh near terms)
!            ifdir=0 FMM activated
!
!  output:
!    AA_u - complex *16(nt)
!      first component of -M_{k}[J] along the srcvals(4:6,i) direction
!        
!    AA_v - complex *16(nt)
!      second component of -M_{k}[J] along
!      the (srcvals(10:12,i) x srcvals(4:6,i)) direction
!            
!    PHI - complex *16(nt) ikn·S_{k}[J]
!

    !List of calling arguments
	real ( kind = 8 ), intent(in) :: eps
	complex ( kind = 8 ), intent(in) :: zpars(3)
    integer, intent(in) :: ns,nt
	real ( kind = 8 ), intent(in) :: srcvals(12,ns),wts(ns),targvals(12,nt)
    complex ( kind = 8 ), intent(in) :: a_u(ns),a_v(ns)
    complex ( kind = 8 ), intent(out) :: AA_u(nt),AA_v(nt),PHI(nt)
	real ( kind = 8 ), intent(in) :: thresh
	integer, intent(in) :: ifdir 


    !List of local variables
	real ( kind = 8 ), allocatable :: n_vect_s(:,:),u_vect_s(:,:)
    real ( kind = 8 ), allocatable :: v_vect_s(:,:),source(:,:)

	real ( kind = 8 ), allocatable :: n_vect_t(:,:),u_vect_t(:,:)
	real ( kind = 8 ), allocatable :: v_vect_t(:,:),targets(:,:)

    complex ( kind = 8 ), allocatable :: a_vect(:,:),b_vect(:,:),lambda(:)
    complex ( kind = 8 ), allocatable :: E(:,:),curlE(:,:),divE(:),rho(:)
	complex ( kind = 8 ) ima,zk,alpha

    integer count1,count2
    integer ifa_vect,ifb_vect,iflambda,ifrho,ifE,ifcurlE,ifdivE
	real ( kind = 8 ) pi
	pi=3.1415926535897932384626433832795028841971d0

	ima=(0.0d0,1.0d0)
	zk=zpars(1)
	alpha=zpars(2)

    allocate(a_vect(3,ns))
    allocate(b_vect(3,ns))
    allocate(lambda(ns))
    allocate(rho(ns))
    allocate(E(3,nt))
    allocate(curlE(3,nt))
    allocate(divE(nt))
    allocate(n_vect_s(3,ns))
    allocate(n_vect_t(3,nt))
    allocate(u_vect_s(3,ns))
    allocate(v_vect_s(3,ns))
    allocate(u_vect_t(3,nt))
    allocate(v_vect_t(3,nt))
    allocate(source(3,ns))
    allocate(targets(3,nt))

	do count1=1,ns
	  n_vect_s(:,count1)=srcvals(10:12,count1)
	  source(:,count1)=srcvals(1:3,count1)
	enddo
	call orthonormalize_all(srcvals(4:6,:),srcvals(10:12,:),&
	 &u_vect_s,v_vect_s,ns)
	
	do count1=1,nt
	  n_vect_t(:,count1)=targvals(10:12,count1)
	  targets(:,count1)=targvals(1:3,count1)
	enddo
	call orthonormalize_all(targvals(4:6,:),targvals(10:12,:),&
	 &u_vect_t,v_vect_t,nt)
	
	
	do count1=1,ns
	  b_vect(1,count1)=(a_u(count1)*u_vect_s(1,count1)+&
	   &a_v(count1)*v_vect_s(1,count1))
	  b_vect(2,count1)=(a_u(count1)*u_vect_s(2,count1)+&
	   &a_v(count1)*v_vect_s(2,count1))
	  b_vect(3,count1)=(a_u(count1)*u_vect_s(3,count1)+&
	   &a_v(count1)*v_vect_s(3,count1))
    enddo

    !Computing the full operator
    ifa_vect=0
    ifb_vect=1
    iflambda=0
    ifrho=0
    ifE=1
    ifcurlE=1
    ifdivE=0


	call Vector_Helmholtz_targ2(eps,zk,ns,source,wts,ifa_vect,a_vect,&
	 &ifb_vect,b_vect,iflambda,lambda,ifrho,rho,n_vect_s,ifE,E,ifcurlE&
	 &,curlE,ifdivE,divE,nt,targets,thresh,ifdir)

    do count1=1,nt
	  b_vect(1,count1)=n_vect_t(2,count1)*curlE(3,count1)-&
	   &n_vect_t(3,count1)*curlE(2,count1)
	  b_vect(2,count1)=n_vect_t(3,count1)*curlE(1,count1)-&
	   &n_vect_t(1,count1)*curlE(3,count1)
	  b_vect(3,count1)=n_vect_t(1,count1)*curlE(2,count1)-&
	   &n_vect_t(2,count1)*curlE(1,count1)
    enddo

    do count1=1,nt
	  AA_u(count1)=-(b_vect(1,count1)*u_vect_t(1,count1)+b_vect(2,count1)&
	   &*u_vect_t(2,count1)+b_vect(3,count1)*u_vect_t(3,count1))
	  AA_v(count1)=-(b_vect(1,count1)*v_vect_t(1,count1)+b_vect(2,count1)&
	   &*v_vect_t(2,count1)+b_vect(3,count1)*v_vect_t(3,count1))
    enddo

    do count1=1,nt
	  PHI(count1)=+ima*zk*(E(1,count1)*n_vect_t(1,count1)+E(2,count1)&
	   &*n_vect_t(2,count1)+E(3,count1)*n_vect_t(3,count1))
	enddo

	deallocate(a_vect)
	deallocate(b_vect)
	deallocate(lambda)
	deallocate(rho)
	deallocate(E)
	deallocate(curlE)
	deallocate(divE)
	
	deallocate(u_vect_s)
	deallocate(v_vect_s)
	deallocate(n_vect_s)
	deallocate(source)

	deallocate(u_vect_t)
	deallocate(v_vect_t)
	deallocate(n_vect_t)
	deallocate(targets)

return
end subroutine em_mfie_pec_FMM_last



subroutine test_accuracy_em_mfie_pec(eps_FMM,sol,zpars,ns,wts,srcvals,&
 &P0,vf,Pt)
implicit none

!
!  This function test the accuracy of the solution computed in the exterior
!  region by testing the extintion theorem in the interior region.
!
!  input:
!    eps_FMM - real *8
!      epsilon for the fmm call
!
!    sol - complex *16(2*ns)
!      induced charge and current on the surface
!      sol(1:ns) - first component of  J along
!      the srcvals(4:6,i) direction
!      sol(ns+1:2*ns) - second component of J along
!      the (srcvals(10:12,i) x srcvals(4:6,i)) direction
!
!    zpars - complex *16(3)
!      kernel parameters
!      zpars(1) = k 
!      zpars(2) = - (not used)
!      zpars(3) = - (not used)
!
!    ns - integer
!      number of sources
!
!    wts - real *8(ns)
!      smooth quadrature weights at original nodes	
!
!    srcvals - real *8(12,ns)
!      xyz(u,v) and derivative info sampled at the 
!      discretization nodes on the surface
!      srcvals(1:3,i) - xyz info
!      srcvals(4:6,i) - dxyz/du info
!      srcvals(7:9,i) - dxyz/dv info
!
!    P0 - real * 8(3)
!      location of the source point at the INTERIOR region
!            
!    vf - complex *16(3)
!      Orientation of the magnetic and electric dipoles located at Pt 
!
!    Pt - real * 8(3)
!      location of the source point at the EXTERIOR region
!
	

    !List of calling arguments
	complex ( kind = 8 ), intent(in) :: zpars(3)
    integer, intent(in) :: ns
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),eps_FMM
    real ( kind = 8 ), intent(in) :: wts(ns),P0(3),Pt(3)
	complex ( kind = 8 ), intent(in) :: sol(2*ns),vf(3)
	
    !List of local variables
	complex ( kind = 8 ) a_u00,a_v00,zk
	complex ( kind = 8 ) ima, R1, R2,Et1(3),Et2(3),aux_cmp,Ht2(3),Ht1(3)
	real ( kind = 8 ) ru_s(3),rv_s(3),n_s(3),sour(3),r,dr(3)
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3),error_E,error_H
	real ( kind = 8 ) pi

	integer count1
	
	ima=(0.0d0,1.0d0)
	pi=3.1415926535897932384626433832795028841971d0
	zk=zpars(1)
	
	write (*,*) 'P0',P0
	call em_mfie_pec_FMM_targ(eps_FMM,zk,ns,srcvals,1,P0,wts,sol(1:ns),&
	 &sol(ns+1:2*ns),Ht1)
		
	call fieldsED(zk,Pt,P0,1,Et2,Ht2,vf,0)
	call fieldsMD(zk,Pt,P0,1,Et2,Ht2,vf,1)

!	
!   Here we are testing the extintion theorem, 
!   that's why we ADD incoming and scattered fields.
!
		
	error_H=sqrt(abs(Ht1(1)+Ht2(1))**2+abs(Ht1(2)+Ht2(2))**2+&
	 &abs(Ht1(3)+Ht2(3))**2)
!	write (*,*) 'Error H: ', error_H
	write (*,*) 'Relative Error H: ', error_H/&
	 &sqrt(abs(Ht2(1))**2+abs(Ht2(2))**2+abs(Ht2(3))**2)
		
return
end subroutine test_accuracy_em_mfie_pec



subroutine em_mfie_pec_FMM_targ(eps,zk,ns,srcvals,nt,targ,wts,a_u,a_v,H)
implicit none
!
!  This funciton computes the fields E,H at a given point using the FMM
!  It doesn't contain near field corrections (it's for debugging purposes)   
!  The representation for the potentials is:
!    Representation:
!
!	    H=curlS_{k}[J]
!
!  input:
!    eps - real *8
!      epsilon for the fmm call
!
!    zk - complex *16
!      Helmholtz parameter 
!
!    ns - integer
!      number of sources
!   
!    nt - integer
!      number of targets
!
!    srcvals - real *8(12,ns)
!      xyz(u,v) and derivative info sampled at the 
!      discretization nodes on the surface
!      srcvals(1:3,i) - xyz info
!      srcvals(4:6,i) - dxyz/du info
!      srcvals(7:9,i) - dxyz/dv info
!
!    targ - real *8(3,nt)
!      location x,y,z of the target points
!
!    wts - real *8(ns)
!      smooth quadrature weights at original nodes
!
!    a_u,a_v - complex *16(ns)
!      two components of the tangent induced current J on the surface
!      along srcvals(4:6,i) and (srcvals(10:12,i) x srcvals(4:6,i)) 
!      directions
! 
!  output:
!    H - complex  *16(3,nt)
!      value of the Magnetic field at the target points


    !List of calling arguments
	real ( kind = 8 ), intent(in) :: eps
	complex ( kind = 8 ), intent(in) :: zk
    integer, intent(in) :: ns,nt
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),targ(3,nt)
    real ( kind = 8 ), intent(in) :: wts(ns)
    complex ( kind = 8 ), intent(in) :: a_u(ns),a_v(ns)
    complex ( kind = 8 ), intent(out) :: H(3,nt)

    !List of local variables
    real ( kind = 8 ), allocatable :: n_vect(:,:),u_vect(:,:),v_vect(:,:)
    real ( kind = 8 ), allocatable :: source(:,:)
    complex ( kind = 8 ), allocatable :: a_vect(:,:),b_vect(:,:)
    complex ( kind = 8 ), allocatable :: lambda(:),rho(:)
    complex ( kind = 8 ), allocatable :: curlE(:,:),gradpot(:,:),divE(:)
	complex ( kind = 8 ) ima
	complex ( kind = 8 ), allocatable :: E(:,:)


    integer count1,count2
    integer ifa_vect,ifb_vect,iflambda,ifrho,ifE,ifcurlE,ifdivE
	real ( kind = 8 ) pi
	pi=3.1415926535897932384626433832795028841971d0

	ima=(0.0d0,1.0d0)

    allocate(a_vect(3,ns))
    allocate(b_vect(3,ns))
    allocate(lambda(ns))
    allocate(rho(ns))
    allocate(curlE(3,nt))
    allocate(gradpot(3,nt))
	allocate(divE(nt))
	allocate(E(3,nt))
	allocate(n_vect(3,ns))
	allocate(u_vect(3,ns))
	allocate(v_vect(3,ns))
	allocate(source(3,ns))

	do count1=1,ns
      n_vect(:,count1)=srcvals(10:12,count1)
	  source(:,count1)=srcvals(1:3,count1)
	enddo
	call orthonormalize_all(srcvals(4:6,:),srcvals(10:12,:),u_vect,&
	 &v_vect,ns)

    do count1=1,ns
	  b_vect(1,count1)=a_u(count1)*u_vect(1,count1)+a_v(count1)*&
	   &v_vect(1,count1)
	  b_vect(2,count1)=a_u(count1)*u_vect(2,count1)+a_v(count1)*&
	   &v_vect(2,count1)
	  b_vect(3,count1)=a_u(count1)*u_vect(3,count1)+a_v(count1)*&
	   &v_vect(3,count1)
    enddo
	
    !Computing the full operator
    ifa_vect=0
    ifb_vect=1
    iflambda=0
    ifrho=0
    ifE=0
    ifcurlE=1
    ifdivE=0

	call Vector_Helmholtz_targ(eps,zk,ns,source,wts,ifa_vect,a_vect,&
	 &ifb_vect,b_vect,iflambda,lambda,ifrho,rho,n_vect,ifE,E,ifcurlE,&
	 &H,ifdivE,divE,nt,targ)

	deallocate(a_vect)
	deallocate(b_vect)
	deallocate(lambda)
	deallocate(curlE)
	deallocate(rho)
	deallocate(gradpot)
	deallocate(divE)
	deallocate(u_vect)
	deallocate(v_vect)
	deallocate(n_vect)
	deallocate(source)

return
end subroutine em_mfie_pec_FMM_targ



subroutine 	get_rhs_em_mfie_pec(P0,vf,alpha,ns,srcvals,zk,RHS)
implicit none

!
!  This function obtains the right hand side for the MFIE formulation
!  for the integral boundary equation:
!
!			J/2 - M_{k}[J] = nxH_inc 
!
!  input:
!    P0 - real * 8 (3)
!      location of the source point at the exterior region
!      WARNING! notice that this formulation uses a representation theorem
!      for the incoming field in the interior region (MFIE) therefore
!      therefore it only works for incoming fields generated by sources in
!      the exterior region (or at infinity like plane waves)
!
!    vf - complex *16(3)
!      Orientation of the magnetic and electric dipoles located at P0 
!
!    alpha - complex *16
!      parameter in the combined formulation
!   
!    ns - integer
!      total number of points on the surface
!
!    srcvals - real *8(12,ns)
!      xyz(u,v) and derivative info sampled at the 
!      discretization nodes on the surface
!      srcvals(1:3,i) - xyz info
!      srcvals(4:6,i) - dxyz/du info
!      srcvals(7:9,i) - dxyz/dv info
!
!    zk - complex *16
!      Helmholtz parameter 
!
!    output:
!      RHS - complex  *16(2*ns)
!        right hand side
!          RHS(1:ns) - first component of  nxH_inc along
!          the srcvals(4:6,i) direction
!          RHS(ns+1:2*ns) - second component of nxH_inc along
!          the (srcvals(10:12,i) x srcvals(4:6,i)) direction
!

	!List of calling arguments
	integer ( kind = 4 ), intent(in) :: ns
	real ( kind = 8 ), intent(in) :: P0(3)
	complex ( kind = 8 ), intent(in) :: vf(3)
	real ( kind = 8 ), intent(in) :: srcvals(12,ns)
	complex ( kind = 8 ), intent(in) :: zk,alpha
	complex ( kind = 8 ), intent(out) :: RHS(2*ns)
	
	!List of local variables
	complex ( kind = 8 ), allocatable :: E(:,:), H(:,:)
	integer count1
	real ( kind = 8 ) ru(3),rv(3),cross_aux(3)
		
	allocate(E(3,ns), H(3,ns))

	call fieldsED(zk,P0,srcvals,ns,E,H,vf,0)
	call fieldsMD(zk,P0,srcvals,ns,E,H,vf,1)
	do count1=1,ns	
	  call orthonormalize(srcvals(4:6,count1),srcvals(10:12,count1),ru,rv)		
	  RHS(count1)=-DOT_PRODUCT(rv,H(:,count1))
	  RHS(ns+count1)=DOT_PRODUCT(ru,H(:,count1))
	enddo

return
end subroutine get_rhs_em_mfie_pec




