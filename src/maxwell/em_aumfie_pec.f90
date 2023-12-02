      subroutine getnearquad_em_aumfie_pec(npatches,norders,&
     &ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,&
     &iquad,rfac0,nquad,wnear)
 
!
!  This subroutine generates the near field quadrature
!  for the auMFIE integral equation:
!
!			rho/2+S'_{k}[rho]=+ikS_{k}[J] + n·E_inc
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
!    - wnear must be of size nquad, corresponding to the kernel S'_{k}[rho]
!
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
!        zpars(2) = not in use
!        zpars(3) = not in use
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
!    - wnear: complex *16(nquad)
!        The desired near field quadrature


      implicit none 
      integer(8) npatches,norders(npatches),npts,nquad
      integer(8) ixyzs(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,rfac0
      integer(8) ndtarg,ntarg
      integer(8) iquadtype
      real *8 targs(ndtarg,ntarg)
      complex *16 zpars(3)
      integer(8) nnz,ipars(2)
      real *8 dpars(1)
      integer(8) row_ptr(ntarg+1),col_ind(nnz),iquad(nnz+1)
      complex *16 wnear(nquad)

      integer(8) ipatch_id(ntarg)
      real *8 uvs_targ(2,ntarg)

      complex *16 alpha,beta
      integer(8) i,j,ndi,ndd,ndz

      integer(8) ipv

!      procedure (), pointer :: fker
!      external h3d_slp, h3d_dlp, h3d_comb

      procedure (), pointer :: fker
	  external  fker_em_aumfie_pec

     ndz=3
	 ndd=1
	 ndi=2
	 ipv=0

      fker =>  fker_em_aumfie_pec
	  ipars(1)=1
	  ipars(2)=1
	  call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
     &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear)

      return
      end subroutine getnearquad_em_aumfie_pec
!
!
!
!

      subroutine lpcomp_em_aumfie_pec_addsub(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,sigma,novers,&
     &nptso,ixyzso,srcover,whtsover,pot,wnear)

!
!  this subroutine evaluates  the S'_{k}[rho] operator
!  the boundary integral equation is (auMFIE):
!
!		rho/2+S'_{k}[rho] = RHS:= +ikS_{k}[J] + n·E_inc
!
!  where the near field is precomputed and stored
!  in the row sparse compressed format.
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
!    - nnz: integer(8) *8
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
!    - wnear: complex *16(nquad)
!        Precomputed near field quadrature S'_{k}[sigma]
!    - sigma - complex *16(ns)
!        induced physical charge density on the surface
!        sigma(1:ns) 
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
!    - pot(1:ntarg) - complex *16 S'_{k}[sigma]


      implicit none
      integer(8) npatches,norder,npols,npts
      integer(8) ndtarg,ntarg
      integer(8) norders(npatches),ixyzs(npatches+1)
      integer(8) ixyzso(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps
      real *8 targs(ndtarg,ntarg)
      complex *16 zpars(3)
      integer(8) nnz,row_ptr(ntarg+1),col_ind(nnz),nquad
      integer(8) iquad(nnz+1)
      complex *16 sigma(npts)
      complex *16 pot(ntarg)
	  
	  complex *16 wnear(nquad)

      integer(8) novers(npatches+1)
      integer(8) nover,npolso,nptso
      real *8 srcover(12,nptso),whtsover(nptso)
      complex *16, allocatable :: potsort(:)

      real *8, allocatable :: sources(:,:),targvals(:,:),wtmp2(:)
      complex *16, allocatable :: charges(:),dipvec(:,:),sigmaover(:)
      integer(8) ns,nt
      complex *16 alpha,beta
      integer(8) ifcharge,ifdipole
      integer(8) ifpgh,ifpghtarg
      complex *16 tmp(10),val,E(1)

      real *8 xmin,xmax,ymin,ymax,zmin,zmax,sizey,sizez,boxsize


      integer(8) i,j,jpatch,jquadstart,jstart


      integer(8) ifaddsub,ifdir

      integer(8) ntj
      
      complex *16 zdotu,pottmp
      complex *16, allocatable :: ctmp2_s(:),dtmp2(:,:)
      real *8 radexp,epsfmm

      complex *16 total_charge

      integer(8) ipars(2)
      real *8 dpars(1),timeinfo(10),t1,t2,omp_get_wtime

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      real *8 thresh,ra
      real *8 rr,rmin
      integer(8) nss,ii,l,npover

      integer(8) nd,ntarg0

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

! 
!       oversample density
!

      call oversample_fun_surf(int(2,8),npatches,norders,ixyzs,iptype,& 
     &npts,sigma(1:npts),novers,ixyzso,ns,sigmaover(1:ns))


      ra = 0

!
!        compute threshold for ignoring local computation
!
      call get_fmm_thresh(int(12,8),ns,srcover,int(12,8),npts,srcvals,thresh)

!
!       fmm call
!
	ifdir=0
      call em_aumfie_pec_FMM(eps,zpars,ns,npts,srcover,targs,whtsover,&
	 &sigmaover(1:ns),pot(1:npts),thresh,ifdir)

      
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
          enddo
        enddo
      enddo
!C$OMP END PARALLEL DO

!
!     Remove near contribution of the FMM
!

      ifdir=1
!C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
!C$OMP$PRIVATE(ctmp2_u,ctmp2_v,wtmp2,nss,l,jstart,ii,E,npover)
      do i=1,ntarg
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          nss = nss + ixyzso(jpatch+1)-ixyzso(jpatch)
        enddo
        allocate(srctmp2(12,nss),ctmp2_s(nss),wtmp2(nss))

        rmin = 1.0d6
        ii = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          jstart = ixyzso(jpatch)-1
          npover = ixyzso(jpatch+1)-ixyzso(jpatch)
          do l=1,npover
			ii = ii+1
			srctmp2(:,ii) = srcover(:,jstart+l)            
			ctmp2_s(ii)=sigmaover(jstart+l)
			wtmp2(ii)=whtsover(jstart+l)
          enddo
        enddo
      call em_aumfie_pec_FMM(eps,zpars,nss,ntarg0,srctmp2,targs(:,i),wtmp2,&
	 &ctmp2_s,E(1),thresh,ifdir)

!		call fker_em_aumfie_pec(nss,ntarg0,srctmp2, targs(:,i), dpars,&
!	   &zpars,ipars,ctmp2_s,wtmp2,E,thresh)
		pot(i) = pot(i) - E(1)

        deallocate(srctmp2,ctmp2_s,wtmp2)
      enddo
      
      call cpu_time(t2)
!C$      t2 = omp_get_wtime()     

      timeinfo(2) = t2-t1


!!      call prin2('quadrature time=*',timeinfo,2)
      
      ttot = timeinfo(1) + timeinfo(2)
!!      call prin2('time in lpcomp=*',ttot,1)

!!!!! This must be used in low frequency regime (to kill the nullspace)
      total_charge=0.0d0
      do i=1,ns
	       total_charge=total_charge+sigmaover(i)*whtsover(i)
	  enddo
	  do i=1,npts
		 pot(i)=pot(i)+total_charge
	  enddo
!!!!! This must be used in low frequency regime (to kill the nullspace)

      return
      end subroutine lpcomp_em_aumfie_pec_addsub



      


subroutine em_aumfie_solver(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,eps,zpars,numit,ifinout,&
     &rhs,eps_gmres,niter,errs,rres,soln)
	 
!
!		This subroutine solves the Scattering Maxwell p.e.c. problem.
!		The the equations are:
!
!			curlE=ikH; curlH =-ikE
!
!		Representation:
!
!			E=ikS_{k}[J]-gradS_{k}[rho] (J is already known)
!
!		Boundary conditions:
!
!			n·E+n·E_inc = rho
!
!		Boundary integral equation:
!
!			rho/2 + S'_{k}[rho] = +ikS_{k}[J] + n·E_inc
!
!     The linear system is solved iteratively using GMRES
!     until a relative residual of 1e-15 is reached
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
!          * zpars(2,3) = not used
!    - numit: integer(8)
!        max number of gmres iterations
!    - ifinout: integer(8)
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
!    - niter: integer(8)
!        number of gmres iterations required for relative residual
!        to converge to eps_gmres
!    - errs:  real *8 (numit)
!        relative residual as a function of iteration
!        number (only errs(1:niter) is populated))
!    - rres: real *8
!        relative residual for computed solution
!    - soln - complex *16(3*npts)
!        soln(1:npts) value of rho
!				 

	 
!
!
!		This subroutine solves the Scattering Maxwell p.e.c. problem.
!		The the equations are:
!
!			curlE=ikH; curlH =-ikE
!
!		Representation:
!
!			E=ikS_{k}[J]-gradS_{k}[rho] (J is already known)
!
!		Boundary conditions:
!
!			n·E+n·E_inc = rho
!
!		The incoming fields must be 'compatible' 
! (solutions to Maxwell's equations in free space)
!
!		Boundary integral equation:
!
!			rho/2 + S'_{k}[rho] = +ikS_{k}[J] + n·E_inc
!
!     The linear system is solved iteratively using GMRES
!     until a relative residual of 1e-15 is reached
!
!
!       input:
!         npatches - integer(8)
!            number of patches
!
!         norders- integer(8)(npatches)
!            order of discretization on each patch 
!
!         ixyzs - integer(8)(npatches+1)
!            ixyzs(i) denotes the starting location in srccoefs,
!               and srcvals array corresponding to patch i
!   
!         iptype - integer(8)(npatches)
!            type of patch
!             iptype = 1, triangular patch discretized using RV nodes
!
!         npts - integer(8)
!            total number of discretization points on the boundary
! 
!         srccoefs - real *8 (9,npts)
!            koornwinder expansion coefficients of xyz, dxyz/du,
!            and dxyz/dv on each patch. 
!            For each point srccoefs(1:3,i) is xyz info
!                           srccoefs(4:6,i) is dxyz/du info
!                           srccoefs(7:9,i) is dxyz/dv info
!
!         srcvals - real *8 (12,npts)
!             xyz(u,v) and derivative info sampled at the 
!             discretization nodes on the surface
!             srcvals(1:3,i) - xyz info
!             srcvals(4:6,i) - dxyz/du info
!             srcvals(7:9,i) - dxyz/dv info
! 
!          eps - real *8
!             precision requested for computing quadrature and fmm
!             tolerance
!
!          zpars - complex *16 (3)
!              kernel parameters (Referring to formula (1))
!              zpars(1) = k 
!              zpars(2) = alpha
!              zpars(3) = - (not used)
!
!          ifinout - integer(8)
!              flag for interior or exterior problems (normals assumed to 
!                be pointing in exterior of region)
!              ifinout = 0, interior problem
!              ifinout = 1, exterior problem
!
!           rhs - complex *16(npts)
!              right hand side
!
!           eps_gmres - real *8
!                gmres tolerance requested
!
!           numit - integer(8)
!              max number of gmres iterations
!
!         output
!           niter - integer(8)
!              number of gmres iterations required for relative residual
!               to converge to 1e-15
!          
!           errs(1:iter) - relative residual as a function of iteration
!              number
! 
!           rres - real *8
!              relative residual for computed solution
!              
!           soln - complex *16(3*npts)
!              soln(1:npts) value of rho
!

      implicit none
      integer(8) npatches,norder,npols,npts
      integer(8) ifinout
      integer(8) norders(npatches),ixyzs(npatches+1)
      integer(8) iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,eps_gmres
      complex *16 zpars(3)
      complex *16 rhs(npts)
      complex *16 soln(npts)

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


!
!   n_var is the number of unknowns in the linear system.
!   as there is only one scalar unknown
!

      n_var=1*npts

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
      allocate(wnear(nquad))
      
!C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,nquad
		wnear(i)=0
      enddo
!C$OMP END PARALLEL DO    


      iquadtype = 1

!!      eps2 = 1.0d-8

      print *, "starting to generate near quadrature"
      call cpu_time(t1)
!C$      t1 = omp_get_wtime()      

      call getnearquad_em_aumfie_pec(npatches,norders,&
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
!

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

        call lpcomp_em_aumfie_pec_addsub(npatches,norders,ixyzs,&
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
!        replace this routine bvimy appropriate layer potential
!        evaluation routine  
!



          call lpcomp_em_aumfie_pec_addsub(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,&
     &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,&
     &soln,novers,npts_over,ixyzso,srcover,wover,wtmp,&
	 &wnear)

            
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
      end subroutine em_aumfie_solver










subroutine fker_em_aumfie_pec(srcinfo, ndt,targinfo,ndd, &
 &dpars,ndz,zpars,ndi,ipars,E_val)
implicit none

!  this function provides the near field kernel that will use
!  zgetnearquad_ggq_guru through getnearquad_em_aumfie_pec

    !List of calling arguments
	integer(8), intent(in) :: ndt,ndd,ndz,ndi
	real ( kind = 8 ), intent(in) :: srcinfo(12)
	real ( kind = 8 ), intent(in) :: targinfo(ndt)
	integer(8), intent(in) :: ipars(ndi)
	real ( kind = 8 ), intent(in) :: dpars(ndd)
	complex ( kind = 8 ), intent(in) :: zpars(ndz)
	complex ( kind = 8 ), intent(out) :: E_val

	
	!List of local variables
!	complex ( kind = 8 ) E_mat(1,1)
	real ( kind = 8 ) ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ) ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ) du(3), dv(3),sour(3),targ(3)
	real ( kind = 8 ) r, dr(3),aux_real
!	complex ( kind = 8 ) nxcurlSka(2,2)
	
	complex ( kind = 8 ) ngradSklambda

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
				
		call get_ngradSklambda(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,&
		&R1,my_exp,r,ngradSklambda)
		E_val=ngradSklambda

return
end subroutine fker_em_aumfie_pec



subroutine em_aumfie_pec_FMM(eps,zpars,ns,nt,srcvals,targvals,wts,rho_in,PHI,&
 &thresh,ifdir)
implicit none

!
!  This subroutine computes the far field contribution if the MFIE operator
!  via FMM
!
!  Boundary integral equation:
!
!			rho/2+S'_{k}[rho]=+ikS_{k}[J] + n·E_inc
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
!    ns - integer(8)
!      number of sources (this is the oversampled set of sources)
!   
!    nt - integer(8) 
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
!    rho_in - complex *16(ns)
!      charge density
!
!    thresh - real *8
!      threshold to remove the selfo interaction term
!
!    ifdir - integer(8)
!      flag, ifdir=1 direct calculation N^2 (used to remove teh near terms)
!            ifdir=0 FMM activated
!
!
!  output:
!    PHI - complex *16(nt) S'_{k}[rho_in]
!

    !List of calling arguments
	real ( kind = 8 ), intent(in) :: eps
	complex ( kind = 8 ), intent(in) :: zpars(3)
    integer(8), intent(in) :: ns,nt
	real ( kind = 8 ), intent(in) :: srcvals(12,ns),wts(ns),targvals(12,nt)
    complex ( kind = 8 ), intent(in) :: rho_in(ns)
    complex ( kind = 8 ), intent(out) :: PHI(nt)
	real ( kind = 8 ), intent(in) :: thresh
	integer(8), intent(in) :: ifdir 

    !List of local variables
	real ( kind = 8 ), allocatable :: n_vect_s(:,:),u_vect_s(:,:),v_vect_s(:,:),source(:,:)
	real ( kind = 8 ), allocatable :: n_vect_t(:,:),u_vect_t(:,:),v_vect_t(:,:),targets(:,:)

    complex ( kind = 8 ), allocatable :: E(:,:),divE(:),rho(:)
	complex ( kind = 8 ) ima,zk,alpha

    integer(8) count1,count2,ier
    integer(8) ifa_vect,ifb_vect,iflambda,ifrho,ifE,ifcurlE,ifdivE
	real ( kind = 8 ) pi
	pi=3.1415926535897932384626433832795028841971d0

	ima=(0.0d0,1.0d0)
	zk=zpars(1)
	alpha=zpars(2)

    allocate(rho(ns))
    allocate(E(3,nt))
    allocate(divE(nt))
    allocate(n_vect_s(3,ns))
    allocate(u_vect_s(3,ns))
    allocate(v_vect_s(3,ns))
    allocate(n_vect_t(3,nt))
    allocate(u_vect_t(3,nt))
    allocate(v_vect_t(3,nt))
    allocate(source(3,ns))
    allocate(targets(3,nt))

    do count1=1,nt
        divE(count1)=0.0d0
        E(1,count1)=0.0d0
        E(2,count1)=0.0d0
        E(3,count1)=0.0d0
    enddo

    do count1=1,ns
      n_vect_s(:,count1)=srcvals(10:12,count1)
      source(:,count1)=srcvals(1:3,count1)
    enddo
    call orthonormalize_all(srcvals(4:6,:),srcvals(10:12,:),u_vect_s,v_vect_s,ns)

    do count1=1,nt
      n_vect_t(:,count1)=targvals(10:12,count1)
      targets(:,count1)=targvals(1:3,count1)
    enddo
    call orthonormalize_all(targvals(4:6,:),targvals(10:12,:),u_vect_t,v_vect_t,nt)

    do count1=1,ns
       rho(count1)=rho_in(count1)*wts(count1)
    enddo

    !Computing the full operator
    if (ifdir.eq.1) then
       call h3ddirectcg(int(1,8),zk,source,rho,ns,targets,nt,divE,E,thresh)
    else
       call hfmm3d_t_c_g(eps,zk,ns,source,rho,nt,targets,divE,E,ier)
    endif
    do count1=1,nt
       E(:,count1)=E(:,count1)/(4.0d0*pi)
    enddo

    do count1=1,nt
      PHI(count1)=&
       +(E(1,count1)*n_vect_t(1,count1)+E(2,count1)*n_vect_t(2,count1)+E(3,count1)*n_vect_t(3,count1))
    enddo

	deallocate(rho)
	deallocate(E)
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
end subroutine em_aumfie_pec_FMM


subroutine test_accuracy_em_aumfie_pec(eps_FMM,sol_v,sol_s,zpars,ns,wts,srcvals,P0,vf,Pt)
implicit none


!
!  This function test the accuracy of the solution computed in the exterior
!  region by testing the extintion theorem in the interior region.
!  Tested E and H fields
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
!      sol(2*ns+1:3*ns) - electric charge density rho
!
!    zpars - complex *16(3)
!      kernel parameters
!      zpars(1) = k 
!      zpars(2) = - (not used)
!      zpars(3) = - (not used)
!
!    ns - integer(8)
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
    integer(8), intent(in) :: ns
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),eps_FMM
    real ( kind = 8 ), intent(in) :: wts(ns),P0(3),Pt(3)
	complex ( kind = 8 ), intent(in) :: sol_v(2*ns),sol_s(ns),vf(3)
	
    !List of local variables
	complex ( kind = 8 ) a_u00,a_v00,zk,alpha
	complex ( kind = 8 ) ima, R1, R2,Et1(3),Et2(3),aux_cmp,Ht2(3),Ht1(3)
	real ( kind = 8 ) ru_s(3),rv_s(3),n_s(3),sour(3),r,dr(3)
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3),error_E,error_H
	real ( kind = 8 ) pi

	integer(8) count1
	
	ima=(0.0d0,1.0d0)
	pi=3.1415926535897932384626433832795028841971d0
	zk=zpars(1)
	alpha=zpars(2)
	
	write (*,*) 'P0',P0
	call em_aumfie_pec_FMM_targ(eps_FMM,zk,ns,srcvals,1,P0,wts,sol_v(1:ns),sol_v(ns+1:2*ns),sol_s,Et1,Ht1)
		
	call fieldsED(zk,Pt,P0,int(1,8),Et2,Ht2,vf,int(0,8))
	call fieldsMD(zk,Pt,P0,int(1,8),Et2,Ht2,vf,int(1,8))

!	
!   Here we are testing the extintion theorem, that's why we ADD incoming and scattered fields.
!

	error_E=sqrt(abs(Et1(1)+Et2(1))**2+abs(Et1(2)+Et2(2))**2+&
	 &abs(Et1(3)+Et2(3))**2)
	write (*,*) 'Relative Error E: ', error_E/&
	 &sqrt(abs(Et2(1))**2+abs(Et2(2))**2+abs(Et2(3))**2)
		
	error_H=sqrt(abs(Ht1(1)+Ht2(1))**2+abs(Ht1(2)+Ht2(2))**2+&
	 &abs(Ht1(3)+Ht2(3))**2)
!	write (*,*) 'Error H: ', error_H
	write (*,*) 'Relative Error H: ', error_H/&
	 &sqrt(abs(Ht2(1))**2+abs(Ht2(2))**2+abs(Ht2(3))**2)
		
return
end subroutine test_accuracy_em_aumfie_pec



subroutine em_aumfie_pec_FMM_targ(eps,zk,ns,srcvals,nt,targ,wts,a_u,a_v,&
 &rho_in,E,H)
implicit none

!
!  This funciton computes the fields E,H at a given point using the FMM
!  It doesn't contain near field corrections (it's for debugging purposes)   
!  The representation for the potentials is:
!    Representation:
!
!	    H=curlS_{k}[J]
!		E=ikS_{k}[J]-gradS_{k}[rho]
!
!  input:
!    eps - real *8
!      epsilon for the fmm call
!
!    zk - complex *16
!      Helmholtz parameter 
!
!    ns - integer(8)
!      number of sources
!   
!    nt - integer(8)
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
!    rho_in - complex *16(ns)
!      induced charge on the surface
! 
!  output:
!    E - complex  *16(3,nt)
!      value of the Electric field at the target points
!    H - complex  *16(3,nt)
!      value of the Magnetic field at the target points

    !List of calling arguments
	real ( kind = 8 ), intent(in) :: eps
	complex ( kind = 8 ), intent(in) :: zk
    integer(8), intent(in) :: ns,nt
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),targ(3,nt)
    real ( kind = 8 ), intent(in) :: wts(ns)
    complex ( kind = 8 ), intent(in) :: a_u(ns),a_v(ns),rho_in(ns)
    complex ( kind = 8 ), intent(out) :: E(3,nt),H(3,nt)

    !List of local variables
	real ( kind = 8 ), allocatable :: n_vect(:,:),u_vect(:,:),v_vect(:,:),source(:,:)
    complex ( kind = 8 ), allocatable :: a_vect(:,:),b_vect(:,:),lambda(:),rho(:)
    complex ( kind = 8 ), allocatable :: curlE(:,:),gradpot(:,:),divE(:)
	complex ( kind = 8 ) ima

    integer(8) count1,count2,ier
    integer(8) ifa_vect,ifb_vect,iflambda,ifrho,ifE,ifcurlE,ifdivE
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
	
	allocate(n_vect(3,ns))
	allocate(u_vect(3,ns))
	allocate(v_vect(3,ns))
	allocate(source(3,ns))

	do count1=1,ns
		n_vect(:,count1)=srcvals(10:12,count1)
		source(:,count1)=srcvals(1:3,count1)
	enddo
	call orthonormalize_all(srcvals(4:6,:),srcvals(10:12,:),u_vect,v_vect,ns)

    do count1=1,ns
        b_vect(1,count1)=a_u(count1)*u_vect(1,count1)+a_v(count1)*v_vect(1,count1)
        b_vect(2,count1)=a_u(count1)*u_vect(2,count1)+a_v(count1)*v_vect(2,count1)
        b_vect(3,count1)=a_u(count1)*u_vect(3,count1)+a_v(count1)*v_vect(3,count1)
    enddo
	
    !Computing the full operator
    ifa_vect=0
    ifb_vect=1
    iflambda=0
    ifrho=0
    ifE=1
    ifcurlE=1
    ifdivE=0

	call Vector_Helmholtz_targ(eps,zk,ns,source,wts,ifa_vect,a_vect,ifb_vect,&
   &b_vect,iflambda,lambda,ifrho,rho,n_vect,ifE,E,ifcurlE,H,ifdivE,divE,nt,targ)

	do count1=1,ns
		rho(count1)=rho_in(count1)*wts(count1)
	enddo


    !Computing the full operator
    call hfmm3d_t_c_g(eps,zk,ns,source,rho,nt,targ,divE,gradpot,ier)

	do count1=1,nt
		gradpot(:,count1)=gradpot(:,count1)/(4.0d0*pi)
		divE(count1)=divE(count1)/(4.0d0*pi)
	enddo
	
	do count1=1,nt
		E(1,count1)=ima*zk*E(1,count1)-gradpot(1,count1)
		E(2,count1)=ima*zk*E(2,count1)-gradpot(2,count1)
		E(3,count1)=ima*zk*E(3,count1)-gradpot(3,count1)
	enddo

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
end subroutine em_aumfie_pec_FMM_targ



subroutine 	get_rhs_em_aumfie_pec(P0,vf,alpha,ns,srcvals,zk,RHS)
implicit none

!
!  This function obtains the right hand side for the auMFIE formulation
!  for the integral boundary equation:
!
!			rho/2 + S'_{k}[rho] = RHS := +ikS_{k}[J] + n·E_inc
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
!    ns - integer(8)
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
!      RHS - complex  *16(ns)
!        right hand side
!          RHS(1:ns) - first component of  nxH_inc along
!          the srcvals(4:6,i) direction
!

	!List of calling arguments
	integer(8), intent(in) :: ns
	real ( kind = 8 ), intent(in) :: P0(3)
	complex ( kind = 8 ), intent(in) :: vf(3)
	real ( kind = 8 ), intent(in) :: srcvals(12,ns)
	complex ( kind = 8 ), intent(in) :: zk,alpha
	complex ( kind = 8 ), intent(out) :: RHS(ns)
	
	!List of local variables
	complex ( kind = 8 ), allocatable :: E(:,:), H(:,:)
	integer(8) count1
	real ( kind = 8 ) ru(3),rv(3),cross_aux(3)
		
	allocate(E(3,ns), H(3,ns))

	call fieldsED(zk,P0,srcvals,ns,E,H,vf,int(0,8))
	call fieldsMD(zk,P0,srcvals,ns,E,H,vf,int(1,8))
	do count1=1,ns	
	  call orthonormalize(srcvals(4:6,count1),srcvals(10:12,count1),ru,rv)		
	  RHS(count1)=DOT_PRODUCT(srcvals(10:12,count1),E(:,count1))
	enddo

return
end subroutine get_rhs_em_aumfie_pec
