      subroutine getnearquad_em_dfie_trans(npatches,norders,&
     &ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,&
     &iquad,rfac0,nquad,wnear)
!
!
!  This subroutine generates the near field quadrature
!  for the integral DFIE:
!
!  Note: the 4 \pi scaling is NOT!! included as the output of the FMM
!  has been rescaled.
!
!
!  The quadrature is computed by the following strategy
!  targets within a sphere of radius rfac0*rs
!  of a chunk centroid is handled using adaptive integration
!  where rs is the radius of the bounding sphere
!  for the patch
!  
!  All other targets in the near field are handled via
!  oversampled quadrature
!
!  The recommended parameter for rfac0 is 1.25d0
!
!        
!
!  input:
!    npatches - integer
!      number of patches
!
!    norders - integer(npatches)
!      order of discretization on each patch 
!
!    ixyzs - integer(npatches+1)
!      starting location of data on patch i
!  
!    iptype - integer(npatches)
!      type of patch
!      iptype = 1 -> triangular patch discretized with RV nodes
!
!    npts - integer
!      total number of discretization points on the boundary
!
!    srccoefs - real *8 (9,npts)
!      koornwinder expansion coefficients of xyz, dxyz/du,
!      and dxyz/dv on each patch. 
!      For each point srccoefs(1:3,i) is xyz info
!        srccoefs(4:6,i) is dxyz/du info
!        srccoefs(7:9,i) is dxyz/dv info
!
!    srcvals - real *8 (12,npts)
!      xyz(u,v) and derivative info sampled at the 
!      discretization nodes on the surface
!      srcvals(1:3,i) - xyz info
!      srcvals(4:6,i) - dxyz/du info
!      srcvals(7:9,i) - dxyz/dv info
! 
!    ndtarg - integer
!      leading dimension of target array
!        
!    ntarg - integer
!      number of targets
!
!    targs - real *8 (ndtarg,ntarg)
!      target information
!
!    ipatch_id - integer(ntarg)
!      id of patch of target i, id = -1, if target is off-surface
!
!         uvs_targ - real *8 (2,ntarg)
!            local uv coordinates on patch if on surface, otherwise
!            set to 0 by default
!          (maybe better to find closest uv on patch using
!            newton)
!            
!          eps - real *8
!             precision requested
!
!          zpars - complex *16(3)
!              kernel parameters
!              zpars(1) = omega 
!              zpars(2) = ep0
!              zpars(3) = mu0
!              zpars(4) = ep1
!              zpars(5) = mu1
!
!           iquadtype - integer
!              quadrature type
!              iquadtype = 1, use ggq for self + adaptive integration
!                 for rest
! 
!
!           nnz - integer
!             number of source patch-> target interactions in the near
!             field
! 
!           row_ptr - integer(ntarg+1)
!              row_ptr(i) is the pointer
!              to col_ind array where list of relevant source patches
!              for target i start
!
!           col_ind - integer (nnz)
!               list of source patches relevant for all targets, sorted
!               by the target number
!
!           iquad - integer(nnz+1)
!               location in wnear array where quadrature for col_ind(i)
!               starts
!
!           rfac0 - integer
!               radius parameter for near field
!
!           nquad - integer
!               number of entries in wnear
!
!        output
!           wnear - complex *16(36*nquad)
!               the desired near field quadrature
!

      implicit none 
      integer npatches,norders(npatches),npts,nquad
      integer ixyzs(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,rfac0
      integer ndtarg,ntarg
      integer iquadtype
      real *8 targs(ndtarg,ntarg)
      complex *16 zpars(5)
      integer nnz,ipars(2)
      real *8 dpars(1)
      integer row_ptr(ntarg+1),col_ind(nnz),iquad(nnz+1)
      complex *16 wnear(36*nquad)
	  
      integer ipatch_id(ntarg)
      real *8 uvs_targ(2,ntarg)

      complex *16 alpha,beta
      integer i,j,ndi,ndd,ndz,count1,count2,icount

      integer ipv
      integer :: t1, t2,clock_rate, clock_max

      procedure (), pointer :: fker
	  external  fker_em_dfie_trans
	  external  em_dfie_trans

     ndz=5
	 ndd=1
	 ndi=2
	 ipv=1

!!      fker => fker_em_dfie_trans
	  fker => em_dfie_trans
!    call system_clock ( t1, clock_rate, clock_max )


	  icount=0
	  do count1=1,6
	    do count2=1,6
		  ipars(1)=count1
		  ipars(2)=count2
		  call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
		   &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
		   &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
		   &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
		   &wnear(icount*nquad+1:(icount+1)*nquad))
		  icount=icount+1
        enddo
	  enddo

!call system_clock ( t2, clock_rate, clock_max )
!write (*,*) 'time quadratures: ',real(t2-t1)/real(clock_rate)
!stop
	  return
      end subroutine getnearquad_em_dfie_trans


      subroutine lpcomp_em_dfie_trans_addsub(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,sigma,novers,&
     &nptso,ixyzso,srcover,whtsover,pot,wnear)

!
!  This subroutine evaluates the layer potential for
!  the DFIE boundary integral equation:
!
!  where the near field is precomputed and stored
!  in the row sparse compressed format.
!
!  Note the 4\pi scaling is NOT included as the FMM output was scaled
!  appropriately
!
!  Note: the identities are not included as the gmres takes care of that
!
!
!  The fmm is used to accelerate the far-field and 
!  near-field interactions are handled via precomputed quadrature
!
!
!  Using add and subtract - no need to call tree and set fmm parameters
!  can directly call existing fmm library
!
!
!  input:
!    npatches - integer
!      number of patches
!
!    norders- integer(npatches)
!      order of discretization on each patch 
!
!    ixyzs - integer(npatches+1)
!      ixyzs(i) denotes the starting location in srccoefs,
!      and srcvals array corresponding to patch i
!   
!    iptype - integer(npatches)
!      type of patch
!      iptype = 1, triangular patch discretized using RV nodes
!
!    npts - integer
!      total number of discretization points on the boundary
! 
!    srccoefs - real *8 (9,npts)
!      koornwinder expansion coefficients of xyz, dxyz/du,
!      and dxyz/dv on each patch. 
!      For each point srccoefs(1:3,i) is xyz info
!        srccoefs(4:6,i) is dxyz/du info
!        srccoefs(7:9,i) is dxyz/dv info
!
!    srcvals - real *8 (12,npts)
!      xyz(u,v) and derivative info sampled at the 
!      discretization nodes on the surface
!      srcvals(1:3,i) - xyz info
!      srcvals(4:6,i) - dxyz/du info
!      srcvals(7:9,i) - dxyz/dv info
! 
!    ndtarg - integer
!      leading dimension of target array
!        
!    ntarg - integer
!      number of targets
!
!    targs - real *8 (ndtarg,ntarg)
!      target information
!
!    eps - real *8
!      precision requested
!
!    zpars - complex *16(3)
!      kernel parameters
!      zpars(1) = omega 
!      zpars(2) = ep0
!      zpars(3) = mu0
!      zpars(4) = ep1
!      zpars(5) = mu1
!
!    nnz - integer *8
!      number of source patch-> target interactions in the near field
! 
!    row_ptr - integer(ntarg+1)
!      row_ptr(i) is the pointer
!      to col_ind array where list of relevant source patches
!      for target i start
!
!    col_ind - integer (nnz)
!      list of source patches relevant for all targets, sorted
!      by the target number
!
!    iquad - integer(nnz+1)
!      location in wnear array where quadrature for col_ind(i)
!      starts
!
!    nquad - integer
!      number of entries in wnear
!
!    wnear  - complex *16(36*nquad)
!      near field precomputed quadrature
!
!    sigma - complex *16(2*ns)
!      induced charge and current on the surface
!      sigma(1:ns) - first component of 'a' along
!        the srcvals(4:6,i) direction
!      sigma(ns+1:2*ns) - second component of 'a' along
!        the (srcvals(10:12,i) x srcvals(4:6,i)) direction
!      sigma(2*ns+1:3*ns) - scalar sigma on the surface
!      sigma(3*ns+1:4*ns) - first component of 'b' along
!        the srcvals(4:6,i) direction
!      sigma(4*ns+1:5*ns) - second component of 'b' along
!        the (srcvals(10:12,i) x srcvals(4:6,i)) direction
!      sigma(5*ns+1:6*ns) - scalar rho on the surface
!
!    novers - integer(npatches)
!      order of discretization for oversampled sources and density
!
!    ixyzso - integer(npatches+1)
!      ixyzso(i) denotes the starting location in srcover,
!      corresponding to patch i
!   
!    nptso - integer
!      total number of oversampled points
!
!    srcover - real *8 (12,nptso)
!      oversampled set of source information
!
!    whtsover - real *8 (nptso)
!      smooth quadrature weights at oversampled nodes
!

      implicit none
      integer npatches,norder,npols,npts
      integer ndtarg,ntarg
      integer norders(npatches),ixyzs(npatches+1)
      integer ixyzso(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps
      real *8 targs(ndtarg,ntarg)
      complex *16 zpars(5)
      integer nnz,row_ptr(ntarg+1),col_ind(nnz),nquad
      integer iquad(nnz+1)
      complex *16 sigma(6*npts),sigma2(npts)
	  
      complex *16 wnear(36*nquad)
	  
      integer novers(npatches+1)
      integer nover,npolso,nptso
      real *8 srcover(12,nptso),whtsover(nptso)
      complex *16 pot(6*ntarg)
      complex *16, allocatable :: potsort(:)

      real *8, allocatable :: sources(:,:),targvals(:,:),wtmp2(:)
      complex *16, allocatable :: charges(:),dipvec(:,:),sigmaover(:)
      integer ns,nt
      complex *16 alpha,beta
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      complex *16 tmp(10),val,E(6)

      real *8 xmin,xmax,ymin,ymax,zmin,zmax,sizey,sizez,boxsize


      integer i,j,jpatch,jquadstart,jstart,count1,count2


      integer ifaddsub,ifdir

      integer ntj
      
      complex *16 zdotu,pottmp
      complex *16, allocatable :: ctmp2_s(:),dtmp2(:,:)
      complex *16, allocatable :: ctmp2_u(:),ctmp2_v(:)
	  complex *16, allocatable :: ctmp2_a_u(:),ctmp2_a_v(:)
	  complex *16, allocatable :: ctmp2_b_u(:),ctmp2_b_v(:)
	  complex *16, allocatable :: ctmp2_rho(:),ctmp2_lambda(:)

	  complex *16, allocatable :: pot_aux(:)

      real *8 radexp,epsfmm

      integer ipars(2)
      real *8 dpars(1),timeinfo(10),t1,t2,omp_get_wtime

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      real *8 thresh,ra
      real *8 rr,rmin
      integer nss,ii,l,npover
      complex *16 ima
	  complex *16 omega,ep0,mu0,ep1,mu1

      integer nd,ntarg0

      real *8 ttot,done,pi

      parameter (nd=1,ntarg0=1)

      ns = nptso
      done = 1
      pi = atan(done)*4
	  ima=(0.0d0,1.0d0)

           
      ifpgh = 0
      ifpghtarg = 1
      allocate(sources(3,ns),targvals(3,ntarg))
      allocate(charges(ns),dipvec(3,ns))
      allocate(sigmaover(6*ns))
	  allocate(pot_aux(6*ntarg))

! 
!       oversample density
    
    call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,&
	&npts,sigma(1:npts),novers,ixyzso,ns,sigmaover(1:ns))
	       
    call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,&
	&npts,sigma(npts+1:2*npts),novers,ixyzso,ns,sigmaover(ns+1:2*ns))

	call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
	&npts,sigma(2*npts+1:3*npts),novers,ixyzso,ns,sigmaover(2*ns+1:3*ns))

	call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
	&npts,sigma(3*npts+1:4*npts),novers,ixyzso,ns,sigmaover(3*ns+1:4*ns))

    call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
    &npts,sigma(4*npts+1:5*npts),novers,ixyzso,ns,sigmaover(4*ns+1:5*ns))

    call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
    &npts,sigma(5*npts+1:6*npts),novers,ixyzso,ns,sigmaover(5*ns+1:6*ns))

      ra = 0


!
!       fmm
!
                                                                        !OUT






		call get_thresh(srcover,ns,targs,ntarg,thresh)
       ifdir=0
	  
		!Calculate the far_field with FMM		
    call em_dfie_trans_FMM(eps,zpars,ns,npts,srcover,targs,whtsover,&
    &sigmaover(1:ns),sigmaover(ns+1:2*ns),sigmaover(2*ns+1:3*ns),&
    &sigmaover(3*ns+1:4*ns),sigmaover(4*ns+1:5*ns),&
    &sigmaover(5*ns+1:6*ns),pot_aux(1:ntarg),pot_aux(ntarg+1:2*ntarg),&
    &pot_aux(2*ntarg+1:3*ntarg),pot_aux(3*ntarg+1:4*ntarg),&
    &pot_aux(4*ntarg+1:5*ntarg),pot_aux(5*ntarg+1:6*ntarg),thresh,ifdir)
                                                                             		

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
		    do count1=0,5
			  do count2=0,5
			    pot_aux(i+count1*npts) = pot_aux(i+count1*npts) + &
				 &wnear((count1*6+count2)*nquad+jquadstart+l-1)*&
				 &sigma(jstart+l-1+npts*count2)
			  enddo
			enddo
          enddo
        enddo
      enddo

!C$OMP END PARALLEL DO
!C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
!C$OMP$PRIVATE(ctmp2_u,ctmp2_v,wtmp2,nss,l,jstart,ii,E,npover)
!		ipars(1)=1
!		ipars(2)=1
	ifdir=1
      do i=1,ntarg
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          nss = nss + ixyzso(jpatch+1)-ixyzso(jpatch)
        enddo
        allocate(srctmp2(12,nss),wtmp2(nss))
        allocate(ctmp2_a_u(nss),ctmp2_a_v(nss))
	    allocate(ctmp2_b_u(nss),ctmp2_b_v(nss))
	    allocate(ctmp2_rho(nss),ctmp2_lambda(nss))

        rmin = 1.0d6
        ii = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          jstart = ixyzso(jpatch)-1
          npover = ixyzso(jpatch+1)-ixyzso(jpatch)
          do l=1,npover
			ii = ii+1
			srctmp2(:,ii) = srcover(:,jstart+l)
			ctmp2_a_u(ii)=sigmaover(jstart+l)
			ctmp2_a_v(ii)=sigmaover(jstart+l+ns)			
			ctmp2_rho(ii)=sigmaover(jstart+l+2*ns)
			ctmp2_b_u(ii)=sigmaover(jstart+l+3*ns)
			ctmp2_b_v(ii)=sigmaover(jstart+l+4*ns)
			ctmp2_lambda(ii)=sigmaover(jstart+l+5*ns)
			
			wtmp2(ii)=whtsover(jstart+l)
          enddo
        enddo
        call em_dfie_trans_FMM(eps,zpars,nss,ntarg0,srctmp2,targs(:,i),&
        &wtmp2,ctmp2_a_u,ctmp2_a_v,ctmp2_rho,ctmp2_b_u,ctmp2_b_v,&
        &ctmp2_lambda,E(1),E(2),E(3),E(4),E(5),E(6),thresh,ifdir)
	 
        do j=0,5
          pot_aux(i+j*ntarg) = pot_aux(i+j*ntarg) - E(j+1)
        enddo
!pot_aux(i) = pot_aux(i) - E(1)
!pot_aux(i+ntarg) = pot_aux(i+ntarg) - E(2)		
!pot_aux(i+2*ntarg) = pot_aux(i+2*ntarg) - E(3)
!pot_aux(i+3*ntarg) = pot_aux(i+3*ntarg) - E(4)
!pot_aux(i+4*ntarg) = pot_aux(i+4*ntarg) - E(5)
!pot_aux(i+5*ntarg) = pot_aux(i+5*ntarg) - E(6)
		
        deallocate(srctmp2,wtmp2)
        deallocate(ctmp2_a_u,ctmp2_a_v,ctmp2_b_u,ctmp2_b_v)
        deallocate(ctmp2_rho,ctmp2_lambda)

      enddo

      omega = zpars(1)
	  ep0 = zpars(2)
      mu0 = zpars(3)
	  ep1 = zpars(4)
	  mu1 = zpars(5)

	  do i=1,ntarg
		pot_aux(i)=pot_aux(i)/(mu0+mu1)
		pot_aux(i+ntarg)=pot_aux(i+ntarg)/(mu0+mu1)
		pot_aux(i+2*ntarg)=pot_aux(i+2*ntarg)/(ep0+ep1)
		pot_aux(i+3*ntarg)=pot_aux(i+3*ntarg)/(ep0+ep1)
		pot_aux(i+4*ntarg)=pot_aux(i+4*ntarg)/(mu0+mu1)
		pot_aux(i+5*ntarg)=-pot_aux(i+5*ntarg)/(ep0+ep1)
	  enddo
      
	  
	  do i=1,ntarg
		pot(i)=pot_aux(i)
		pot(i+ntarg)=pot_aux(i+ntarg)
		pot(i+2*ntarg)=pot_aux(i+4*ntarg)
		pot(i+3*ntarg)=pot_aux(i+2*ntarg)
		pot(i+4*ntarg)=pot_aux(i+3*ntarg)
		pot(i+5*ntarg)=pot_aux(i+5*ntarg)
	  enddo
	  
	  return
      end subroutine lpcomp_em_dfie_trans_addsub



subroutine em_dfie_trans_solver(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,eps,zpars,numit,ifinout,&
     &rhs,eps_gmres,niter,errs,rres,soln)
!
!  This subroutine solves the Scattering Maxwell homogeneous dielectric
!  problem.
!  The the equations are:
!
!    curlE0=ik0H0; curlH0 =-ik0E0   (exterior region)
!    curlE1=ik1H1; curlH1 =-ik1E1   (interior region)
!
!  parameters: k0=omega*sqrt(ep0*mu0); k1=omega*sqrt(ep1*mu1)
!
!  Representation: (1)
!
!    E0 = mu0curlS_{k0}[a] - mu0S_{k0}[n·sigma] + mu0ep0S_{k0}[b] + 
!       + gradS_{k0}[rho]
!    E1 = mu1curlS_{k1}[a] - mu1S_{k1}[n·sigma] + mu1ep1S_{k1}[b] + 
!       + gradS_{k1}[rho]
!
!  Boundary conditions: (2)
!
!                    nxE0-nxE1 = -nxE_inc
!                  divE0-divE1 = 0
!    nxcurlE0/mu0-nxcurlE1/mu1 = -nxcurlE_inc/mu0
!              n·(ep0E0-ep1E1) = -n·ep0·E_inc
!
!
!  The incoming fields must be 'compatible' 
!  (solutions to Maxwell's equations in free space)
!
!  Boundary integral equation:
!
!    equation (86)of paper https://doi.org/10.1080/03605302.2018.1446447
!    is the result of hitting the representation (1) with the boundary 
!    conditions (2) and rescale the the equations with 1/(mu0+mu1) 
!    and 1/(ep0+ep1) to get 1/2 in the diagonal
!
!  The linear system is solved iteratively using GMRES
!  until a relative residual of 1e-15 is reached
!
!  input:
!    npatches - integer
!      number of patches
!
!    norders- integer(npatches)
!      order of discretization on each patch 
!
!    ixyzs - integer(npatches+1)
!      ixyzs(i) denotes the starting location in srccoefs,
!      and srcvals array corresponding to patch i
!   
!    iptype - integer(npatches)
!      type of patch
!      iptype = 1, triangular patch discretized using RV nodes
!
!    npts - integer
!      total number of discretization points on the boundary
! 
!    srccoefs - real *8 (9,npts)
!      koornwinder expansion coefficients of xyz, dxyz/du,
!      and dxyz/dv on each patch. 
!      For each point srccoefs(1:3,i) is xyz info
!        srccoefs(4:6,i) is dxyz/du info
!        srccoefs(7:9,i) is dxyz/dv info
!
!    srcvals - real *8 (12,npts)
!      xyz(u,v) and derivative info sampled at the 
!      discretization nodes on the surface
!      srcvals(1:3,i) - xyz info
!      srcvals(4:6,i) - dxyz/du info
!      srcvals(7:9,i) - dxyz/dv info
! 
!    eps - real *8
!      precision requested for computing quadrature and fmm
!      tolerance
!
!    zpars - complex *16 (5)
!      kernel parameters 
!      zpars(1) = omega 
!      zpars(2) = ep0
!      zpars(3) = mu0
!      zpars(4) = ep1
!      zpars(5) = mu1
!
!    ifinout - integer
!      flag for interior or exterior problems (normals assumed to 
!        be pointing in exterior of region)
!      ifinout = 0, interior problem
!      ifinout = 1, exterior problem
!
!    rhs - complex *16(npts)
!      right hand side
!
!    eps_gmres - real *8
!      gmres tolerance requested
!
!    numit - integer
!      max number of gmres iterations
!
!    output
!      niter - integer
!      number of gmres iterations required for relative residual
!      to converge to 1e-15
!          
!    errs(1:iter) - relative residual as a function of iteration
!      number
! 
!    rres - real *8
!      relative residual for computed solution
!              
!    soln - complex *16(3*npts)
!      soln(1:npts) component of the tangent induced current a on the 
!        surface along srcvals(4:6,i) direction
!      soln(npts+1:2*npts) component of the tangent induced current a 
!        on the surface along (srcvals(10:12,i) x srcvals(4:6,i)) 
!        direction
!      soln(2*npts+1:3*npts) rho on the surface 
!      soln(3*npts+1:4*npts) component of the tangent induced current
!        a on the surface along srcvals(4:6,i) direction
!      soln(4*npts+1:5*npts) component of the tangent induced current
!        a on the surface along (srcvals(10:12,i) x srcvals(4:6,i))
!        direction
!      soln(5*npts+1:6*npts) rho on the surface 
!

      implicit none
      integer npatches,norder,npols,npts
      integer ifinout
      integer norders(npatches),ixyzs(npatches+1)
      integer iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,eps_gmres
      complex *16 zpars(5)
      complex *16 rhs(6*npts)
      complex *16 soln(6*npts)

      real *8, allocatable :: targs(:,:)
      real *8, allocatable :: targs_aux(:,:)
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
      complex *16 ima

      ima=(0.0d0,1.0d0)
	  
!
!   n_var is the number of unknowns in the linear system.
!   as we have tow vector unknown a,b and two scalars rho,sigma
!   we need n_var=6*npts
!

      n_var=6*npts

      allocate(vmat(n_var,numit+1),hmat(numit,numit))
      allocate(cs(numit),sn(numit))
      allocate(wtmp(n_var),svec(numit+1),yvec(numit+1))


      done = 1
      pi = atan(done)*4


!
!        setup targets as on surface discretization points
! 
      ndtarg = 12
      ntarg = npts
      allocate(targs(ndtarg,npts),uvs_targ(2,ntarg),ipatch_id(ntarg))
      allocate(targs_aux(3,npts))
!C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
	    targs_aux(1,i) = srcvals(1,i)
		targs_aux(2,i) = srcvals(2,i)
		targs_aux(3,i) = srcvals(3,i)
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
      call findnearmem(cms,npatches,rad_near,targs_aux,npts,nnz)
      print *, "nnz=",nnz

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,targs_aux,npts,row_ptr,&
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
      call prinf('novers=*',novers,100)
	write (*,*) 'inside no collapse'

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
      allocate(wnear(36*nquad))
	  
!C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,36*nquad
	    wnear(i)=0
      enddo
!C$OMP END PARALLEL DO    

      iquadtype = 1

!!      eps2 = 1.0d-8

      print *, "starting to generate near quadrature"
      call cpu_time(t1)
!C$      t1 = omp_get_wtime()      
	 
	 	 write (*,*) 'generation doing'

      call getnearquad_em_dfie_trans(npatches,norders,&
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


        call lpcomp_em_dfie_trans_addsub(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,&
     &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,&
     &vmat(1,it),novers,npts_over,ixyzso,srcover,wover,wtmp,wnear)

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


          call lpcomp_em_dfie_trans_addsub(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,&
     &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,&
     &soln,novers,npts_over,ixyzso,srcover,wover,wtmp,wnear)

            
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
      end subroutine em_dfie_trans_solver




subroutine fker_em_dfie_trans(srcinfo, ndt,targinfo,ndd, dpars,ndz,&
 &zpars,ndi,ipars,E_val)
implicit none

!  THIS FUNCTION IS OBSOLETE, NOT USED, subroutine em_dfie_trans is ued
!  instead (much faster), but very hard to read..
!  this function provides the near field kernel that will use 
!  zgetnearquad_ggq_guru 
!  through getnearquad_DFIE

    !List of calling arguments
	integer, intent(in) :: ndt,ndd,ndz,ndi
	real ( kind = 8 ), intent(in) :: srcinfo(12)
	real ( kind = 8 ), intent(in) :: targinfo(ndt)
	integer, intent(in) :: ipars(ndi)
	real ( kind = 8 ), intent(in) :: dpars(ndd)
	complex ( kind = 8 ), intent(in) :: zpars(ndz)
	complex ( kind = 8 ), intent(out) :: E_val
	
	!List of local variables
	real ( kind = 8 ) ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ) ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ) du(3), dv(3),sour(3),targ(3)
	real ( kind = 8 ) r, dr(3),aux_real
	complex ( kind = 8 ) nxcurlSk0a(2,2),nxcurlSk1a(2,2),nxcurlcurlSk0a(2,2)
	complex ( kind = 8 ) nxcurlcurlSk1a(2,2),ncurlSk1a(1,2),ncurlSk0a(1,2)
	complex ( kind = 8 ) nxSk0b(2,2),nxSk1b(2,2),divSk0b(1,2),divSk1b(1,2)
	complex ( kind = 8 ) nxSk0nrho(2,1),nxSk1nrho(2,1),nxcurlSk0nrho(2,1)
	complex ( kind = 8 ) nxcurlSk1nrho(2,1),Dk0rho(1,1),Dk1rho(1,1)
	complex ( kind = 8 ) nxgradSk0lambda(2,1),nxgradSk1lambda(2,1),Sk0lambda(1,1)
	complex ( kind = 8 ) Sk1lambda(1,1),ngradSk0lambda(1,1),nSk1nrho(1,1)
	complex ( kind = 8 ) ngradSk1lambda(1,1),nSk0b(1,2),nSk1b(1,2),nSk0nrho(1,1)
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3),xprod_aux3(3),xprod_aux4(3)
	complex ( kind = 8 ) R1_0,R1_1,R2_0,R2_1,ima,my_exp_0,my_exp_1,omega,ep0	
	complex ( kind = 8 ) mu0,ep1,mu1,zk0,zk1
	real ( kind = 8 ) pi
	integer count1,count2
	complex ( kind = 8 )  E_mat(6,6)

	
	pi=3.1415926535897932384626433832795028841971d0
	ima=(0.0d0,1.0d0)
	omega=zpars(1)
	ep0=zpars(2)
	mu0=zpars(3)
	ep1=zpars(4)
	mu1=zpars(5)
	
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
	zk0=omega*sqrt(ep0*mu0)
	zk1=omega*sqrt(ep1*mu1)

	R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
	R2_0=((ima*zk0)**2/r**3-3.0d0*ima*zk0/r**4+3.0d0/r**5)*exp(ima*zk0*r)/&
	 &(4.0d0*pi)
	my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
	R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
	R2_1=((ima*zk1)**2/r**3-3.0d0*ima*zk1/r**4+3.0d0/r**5)*exp(ima*zk1*r)/&
	 &(4.0d0*pi)
	my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
	
	call orthonormalize(srcinfo(4:6),n_s,ru_s,rv_s)
	call orthonormalize(targinfo(4:6),n_t,ru_t,rv_t)

	call get_nxcurlSka(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_0,my_exp_0,r,nxcurlSk0a)		
	call get_nxcurlcurlSka(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_0,R2_0,zk0,&
	 &my_exp_0,r,nxcurlcurlSk0a)
		
	call get_nxcurlSka(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_1,my_exp_1,r,nxcurlSk1a)		
	call get_nxcurlcurlSka(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_1,R2_1,zk1,&
	 &my_exp_1,r,nxcurlcurlSk1a)

	call get_nxSkb(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,my_exp_0,r,nxSk0b)
	call get_nxSkb(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,my_exp_1,r,nxSk1b)
		
	call get_divSkb(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_0,my_exp_0,r,divSk0b)
	call get_divSkb(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_1,my_exp_1,r,divSk1b)
		
	call get_nSkb(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,my_exp_0,r,nSk0b)
	call get_nSkb(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,my_exp_1,r,nSk1b)
		
	call get_ncurlSka(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_0,my_exp_0,r,ncurlSk0a)
	call get_ncurlSka(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_1,my_exp_1,r,ncurlSk1a)

	call get_nxSknrho(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,my_exp_0,r,nxSk0nrho)
	call get_nxSknrho(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,my_exp_1,r,nxSk1nrho)
		
	call get_nxcurlSknrho(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_0,my_exp_0,r,&
	 &nxcurlSk0nrho)
	call get_nxcurlSknrho(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_1,my_exp_1,r,&
	 &nxcurlSk1nrho)


	call get_Dkrho(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_0,my_exp_0,r,Dk0rho)
	call get_Dkrho(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_1,my_exp_1,r,Dk1rho)
		
	call get_nSknrho(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,my_exp_0,r,nSk0nrho)
	call get_nSknrho(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,my_exp_1,r,nSk1nrho)


	call get_nxgradSklambda(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_0,my_exp_0,r,&
	 &nxgradSk0lambda)	
	call get_nxgradSklambda(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_1,my_exp_1,r,&
	 &nxgradSk1lambda)	
		
	call get_Sklambda(my_exp_0,r,Sk0lambda)
	call get_Sklambda(my_exp_1,r,Sk1lambda)

	call get_ngradSklambda(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_0,my_exp_0,r,&
	 &ngradSk0lambda)
	call get_ngradSklambda(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_1,my_exp_1,r,&
	 &ngradSk1lambda)

	E_mat(1,1)=mu0*nxcurlSk0a(1,1)-mu1*nxcurlSk1a(1,1)
	E_mat(1,2)=mu0*nxcurlSk0a(1,2)-mu1*nxcurlSk1a(1,2)
	E_mat(2,1)=mu0*nxcurlSk0a(2,1)-mu1*nxcurlSk1a(2,1)
	E_mat(2,2)=mu0*nxcurlSk0a(2,2)-mu1*nxcurlSk1a(2,2)

	E_mat(3,4)=ep0*nxcurlSk0a(1,1)-ep1*nxcurlSk1a(1,1)
	E_mat(3,5)=ep0*nxcurlSk0a(1,2)-ep1*nxcurlSk1a(1,2)
	E_mat(4,4)=ep0*nxcurlSk0a(2,1)-ep1*nxcurlSk1a(2,1)
	E_mat(4,5)=ep0*nxcurlSk0a(2,2)-ep1*nxcurlSk1a(2,2)

	E_mat(3,1)=(nxcurlcurlSk0a(1,1)-nxcurlcurlSk1a(1,1))
	E_mat(3,2)=(nxcurlcurlSk0a(1,2)-nxcurlcurlSk1a(1,2))
	E_mat(4,1)=(nxcurlcurlSk0a(2,1)-nxcurlcurlSk1a(2,1))
	E_mat(4,2)=(nxcurlcurlSk0a(2,2)-nxcurlcurlSk1a(2,2))

	E_mat(1,4)=ep0*mu0*nxSk0b(1,1)-ep1*mu1*nxSk1b(1,1)
	E_mat(1,5)=ep0*mu0*nxSk0b(1,2)-ep1*mu1*nxSk1b(1,2)
	E_mat(2,4)=ep0*mu0*nxSk0b(2,1)-ep1*mu1*nxSk1b(2,1)
	E_mat(2,5)=ep0*mu0*nxSk0b(2,2)-ep1*mu1*nxSk1b(2,2)

	E_mat(5,1)=0.0d0
	E_mat(5,2)=0.0d0

	E_mat(3,6)=0.0d0
	E_mat(4,6)=0.0d0
		
	E_mat(5,4)=ep0*mu0*divSk0b(1,1)-ep1*mu1*divSk1b(1,1)
	E_mat(5,5)=ep0*mu0*divSk0b(1,2)-ep1*mu1*divSk1b(1,2)
		
	E_mat(6,4)=mu0*ep0**2*nSk0b(1,1)-mu1*ep1**2*nSk1b(1,1)
	E_mat(6,5)=mu0*ep0**2*nSk0b(1,2)-mu1*ep1**2*nSk1b(1,2)
		
	E_mat(6,1)=ep0*mu0*ncurlSk0a(1,1)-ep1*mu1*ncurlSk1a(1,1)
	E_mat(6,2)=ep0*mu0*ncurlSk0a(1,2)-ep1*mu1*ncurlSk1a(1,2)
		
	E_mat(1,3)=-mu0*nxSk0nrho(1,1)+mu1*nxSk1nrho(1,1)
	E_mat(2,3)=-mu0*nxSk0nrho(2,1)+mu1*nxSk1nrho(2,1)
		
	E_mat(3,3)=-nxcurlSk0nrho(1,1)+nxcurlSk1nrho(1,1)
	E_mat(4,3)=-nxcurlSk0nrho(2,1)+nxcurlSk1nrho(2,1)
		
	E_mat(5,3)=mu0*Dk0rho(1,1)-mu1*Dk1rho(1,1)
	E_mat(6,3)=-ep0*mu0*nSk0nrho(1,1)+ep1*mu1*nSk1nrho(1,1)

	E_mat(1,6)=nxgradSk0lambda(1,1)-nxgradSk1lambda(1,1)
	E_mat(2,6)=nxgradSk0lambda(2,1)-nxgradSk1lambda(2,1)
		
	E_mat(5,6)=-zk0**2*Sk0lambda(1,1)+zk1**2*Sk1lambda(1,1)
	E_mat(6,6)=ep0*ngradSk0lambda(1,1)-ep1*ngradSk1lambda(1,1)
		
	E_val=E_mat(ipars(1),ipars(2))


return
end subroutine fker_em_dfie_trans



subroutine em_dfie_trans_FMM(eps,zpars,ns,nt,srcvals,targvals,wts,&
 &a_u,a_v,rho_in,b_u,b_v,lambda_in,AA_u,AA_v,BB_u,BB_v,PHI,PSI,&
 &thresh,ifdir)
implicit none

    !List of calling arguments
    real ( kind = 8 ), intent(in) :: eps
    complex ( kind = 8 ), intent(in) :: zpars(5)
    integer, intent(in) :: ns,nt
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),wts(ns)
    real ( kind = 8 ), intent(in) :: targvals(12,nt)
    complex ( kind = 8 ), intent(in) :: a_u(ns),a_v(ns),rho_in(ns)
    complex ( kind = 8 ), intent(in) :: b_u(ns),b_v(ns),lambda_in(ns)
    complex ( kind = 8 ), intent(out) :: AA_u(nt),AA_v(nt),PHI(nt)
    complex ( kind = 8 ), intent(out) :: BB_u(nt),BB_v(nt),PSI(nt)
    real ( kind = 8 ), intent(in) :: thresh
    integer, intent(in) :: ifdir 


    !List of local variables
	real ( kind = 8 ), allocatable :: u_vect_s(:,:),v_vect_s(:,:)
	real ( kind = 8 ), allocatable :: source(:,:),n_vect_s(:,:)

	real ( kind = 8 ), allocatable :: n_vect_t(:,:),u_vect_t(:,:)
	real ( kind = 8 ), allocatable :: targets(:,:),v_vect_t(:,:)

    complex ( kind = 8 ), allocatable :: a_vect(:,:),b_vect(:,:)
    complex ( kind = 8 ), allocatable :: lambda(:),rho(:)
    complex ( kind = 8 ), allocatable :: E(:,:),curlE(:,:),divE(:)
    complex ( kind = 8 ) ima,zk0,zk1

    complex ( kind = 8 ) omega,ep0,mu0,ep1,mu1

    integer count1,count2
    integer ifa_vect,ifb_vect,iflambda,ifrho,ifE,ifcurlE,ifdivE

	ima=(0.0d0,1.0d0)

	omega=zpars(1)
	ep0=zpars(2)
	mu0=zpars(3)
	ep1=zpars(4)
	mu1=zpars(5)
	
	zk0=omega*sqrt(ep0*mu0)
	zk1=omega*sqrt(ep1*mu1)


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
	call orthonormalize_all(srcvals(4:6,:),srcvals(10:12,:),u_vect_s,&
	 &v_vect_s,ns)
	
	do count1=1,nt
      n_vect_t(:,count1)=targvals(10:12,count1)
      targets(:,count1)=targvals(1:3,count1)
	enddo
	call orthonormalize_all(targvals(4:6,:),targvals(10:12,:),u_vect_t,&
	 &v_vect_t,nt)
	
	do count1=1,ns
      a_vect(1,count1)=mu0*(a_u(count1)*u_vect_s(1,count1)+a_v(count1)*&
       &v_vect_s(1,count1))
      a_vect(2,count1)=mu0*(a_u(count1)*u_vect_s(2,count1)+a_v(count1)*&
       &v_vect_s(2,count1))
      a_vect(3,count1)=mu0*(a_u(count1)*u_vect_s(3,count1)+a_v(count1)*&
       &v_vect_s(3,count1))
				
      b_vect(1,count1)=mu0*ep0*(b_u(count1)*u_vect_s(1,count1)+&
	   &b_v(count1)*v_vect_s(1,count1))
      b_vect(2,count1)=mu0*ep0*(b_u(count1)*u_vect_s(2,count1)+&
	   &b_v(count1)*v_vect_s(2,count1))
      b_vect(3,count1)=mu0*ep0*(b_u(count1)*u_vect_s(3,count1)+&
	   &b_v(count1)*v_vect_s(3,count1))
		
      rho(count1)=-mu0*rho_in(count1)
      lambda(count1)=lambda_in(count1)
	enddo


    !Computing the full operator
    ifa_vect=1
    ifb_vect=1
    iflambda=1
    ifrho=1
    ifE=1
    ifcurlE=1
    ifdivE=1

	call Vector_Helmholtz_targ2(eps,zk0,ns,source,wts,ifa_vect,a_vect,&
	 &ifb_vect,b_vect,iflambda,lambda,ifrho,rho,n_vect_s,ifE,E,ifcurlE,&
	 &curlE,ifdivE,divE,nt,targets,thresh,ifdir)

!	call Vector_Helmholtz_targ(eps,izk,ns,source,wts,ifa_vect,a_vect,ifb_vect,&
!	 &b_vect,iflambda,lambda,ifrho,rho,n_vect_s,ifE,E,ifcurlE,curlE,ifdivE,divE,nt,targets)

    do count1=1,nt
      b_vect(1,count1)=n_vect_t(2,count1)*E(3,count1)-&
	   &n_vect_t(3,count1)*E(2,count1)
      b_vect(2,count1)=n_vect_t(3,count1)*E(1,count1)-&
	   &n_vect_t(1,count1)*E(3,count1)
      b_vect(3,count1)=n_vect_t(1,count1)*E(2,count1)-&
	   &n_vect_t(2,count1)*E(1,count1)
    enddo

    do count1=1,nt
      AA_u(count1)=b_vect(1,count1)*u_vect_t(1,count1)+&
	   &b_vect(2,count1)*u_vect_t(2,count1)+&
	   &b_vect(3,count1)*u_vect_t(3,count1)
      AA_v(count1)=b_vect(1,count1)*v_vect_t(1,count1)+&
	   &b_vect(2,count1)*v_vect_t(2,count1)+&
	   &b_vect(3,count1)*v_vect_t(3,count1)
    enddo



    do count1=1,nt
      b_vect(1,count1)=n_vect_t(2,count1)*curlE(3,count1)-&
	   &n_vect_t(3,count1)*curlE(2,count1)
      b_vect(2,count1)=n_vect_t(3,count1)*curlE(1,count1)-&
	   &n_vect_t(1,count1)*curlE(3,count1)
      b_vect(3,count1)=n_vect_t(1,count1)*curlE(2,count1)-&
	   &n_vect_t(2,count1)*curlE(1,count1)
    enddo

    do count1=1,nt
      BB_u(count1)=(b_vect(1,count1)*u_vect_t(1,count1)+&
	   &b_vect(2,count1)*u_vect_t(2,count1)+&
	   &b_vect(3,count1)*u_vect_t(3,count1))/mu0
      BB_v(count1)=(b_vect(1,count1)*v_vect_t(1,count1)+&
	   &b_vect(2,count1)*v_vect_t(2,count1)+&
	   &b_vect(3,count1)*v_vect_t(3,count1))/mu0
    enddo
		
	do count1=1,nt
      PHI(count1)=divE(count1)
      PSI(count1)=ep0*(n_vect_t(1,count1)*E(1,count1)+&
	   &n_vect_t(2,count1)*E(2,count1)+n_vect_t(3,count1)*E(3,count1))
	enddo
!!!!!


	do count1=1,ns
      a_vect(1,count1)=mu1*(a_u(count1)*u_vect_s(1,count1)+a_v(count1)*&
	   &v_vect_s(1,count1))
      a_vect(2,count1)=mu1*(a_u(count1)*u_vect_s(2,count1)+a_v(count1)*&
	   &v_vect_s(2,count1))
      a_vect(3,count1)=mu1*(a_u(count1)*u_vect_s(3,count1)+a_v(count1)*&
	   &v_vect_s(3,count1))
				
      b_vect(1,count1)=mu1*ep1*(b_u(count1)*u_vect_s(1,count1)+&
	   &b_v(count1)*v_vect_s(1,count1))
      b_vect(2,count1)=mu1*ep1*(b_u(count1)*u_vect_s(2,count1)+&
	   &b_v(count1)*v_vect_s(2,count1))
      b_vect(3,count1)=mu1*ep1*(b_u(count1)*u_vect_s(3,count1)+&
	   &b_v(count1)*v_vect_s(3,count1))
		
      rho(count1)=-mu1*rho_in(count1)
      lambda(count1)=lambda_in(count1)
	enddo
	
	
	  !Computing the full operator
    ifa_vect=1
    ifb_vect=1
    iflambda=1
    ifrho=1
    ifE=1
    ifcurlE=1
    ifdivE=1

	call Vector_Helmholtz_targ2(eps,zk1,ns,source,wts,ifa_vect,a_vect,&
	 &ifb_vect,b_vect,iflambda,lambda,ifrho,rho,n_vect_s,ifE,E,ifcurlE,&
	 &curlE,ifdivE,divE,nt,targets,thresh,ifdir)

	
	
	do count1=1,nt
      b_vect(1,count1)=n_vect_t(2,count1)*E(3,count1)-&
	   &n_vect_t(3,count1)*E(2,count1)
      b_vect(2,count1)=n_vect_t(3,count1)*E(1,count1)-&
	   &n_vect_t(1,count1)*E(3,count1)
      b_vect(3,count1)=n_vect_t(1,count1)*E(2,count1)-&
	   &n_vect_t(2,count1)*E(1,count1)
    enddo

    do count1=1,nt
      AA_u(count1)=AA_u(count1)-(b_vect(1,count1)*u_vect_t(1,count1)+&
	   &b_vect(2,count1)*u_vect_t(2,count1)+&
	   &b_vect(3,count1)*u_vect_t(3,count1))
      AA_v(count1)=AA_v(count1)-(b_vect(1,count1)*v_vect_t(1,count1)+&
	   &b_vect(2,count1)*v_vect_t(2,count1)+&
	   &b_vect(3,count1)*v_vect_t(3,count1))
    enddo

    do count1=1,nt
      b_vect(1,count1)=n_vect_t(2,count1)*curlE(3,count1)-&
	   &n_vect_t(3,count1)*curlE(2,count1)
      b_vect(2,count1)=n_vect_t(3,count1)*curlE(1,count1)-&
	   &n_vect_t(1,count1)*curlE(3,count1)
      b_vect(3,count1)=n_vect_t(1,count1)*curlE(2,count1)-&
	   &n_vect_t(2,count1)*curlE(1,count1)
    enddo

    do count1=1,nt
      BB_u(count1)=BB_u(count1)-(b_vect(1,count1)*u_vect_t(1,count1)+&
	   &b_vect(2,count1)*u_vect_t(2,count1)+b_vect(3,count1)*&
	   &u_vect_t(3,count1))/mu1
      BB_v(count1)=BB_v(count1)-(b_vect(1,count1)*v_vect_t(1,count1)+&
	   &b_vect(2,count1)*v_vect_t(2,count1)+b_vect(3,count1)*&
	   &v_vect_t(3,count1))/mu1
    enddo
		
	do count1=1,nt
      PHI(count1)=PHI(count1)-divE(count1)
      PSI(count1)=PSI(count1)-ep1*(n_vect_t(1,count1)*E(1,count1)&
       &+n_vect_t(2,count1)*E(2,count1)+n_vect_t(3,count1)*E(3,count1))
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
end subroutine em_dfie_trans_FMM



subroutine 	get_rhs_em_dfie_trans(P0,Pt,vf,ns,srcvals,zpars,RHS_out)
implicit none
!
!  This function obtains the right hand side for the DFIE formulation
!
!  input:
!    P0 - real * 8 (3)
!      location of the source point at the exterior region
!      WARNING! notice that this formulation uses a representation
!      theorem for the incoming field in the interior region (MFIE)
!      therefore therefore it only works for incoming fields generated
!      by sources in the exterior region (or at infinity like plane
!      waves)
!
!    vf - complex *16(3)
!      Orientation of the magnetic and electric dipoles located at P0 
!      and Pt
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
!  output:
!    RHS - complex  *16(6*ns)
!      right hand side
!      RHS(1:ns) - first component of -nxE_inc along the srcvals(4:6,i)
!       direction scaled with 1/(mu0+mu1) to get 1/2 as identity
!      RHS(ns+1:2*ns) - second component of -nxE_inc along the
!       (srcvals(10:12,i) x srcvals(4:6,i)) direction scaled with 
!       1/(mu0+mu1) to get 1/2 as identity
!      RHS(2*ns+1:3*ns) - divergence of E_inc (should be zero) scaled 
!       with 1/(mu0+mu1) to get 1/2 as identity
!      RHS(3*ns+1:4*ns) - first component of -nxE_inc along the 
!       srcvals(4:6,i) direction scaled with 1/(ep0+ep1) to get 1/2 as 
!       identity
!      RHS(4*ns+1:5*ns) - second component of -nxE_inc along the
!       (srcvals(10:12,i) x srcvals(4:6,i)) direction  scaled 
!       with 1/(ep0+ep1) to get 1/2 as identity
!      RHS(5*ns+1:6*ns) - normal component of the incoming 
!       field -n·ep0·E_inc scaled with 1/(ep0+ep1) to get 1/2 as 
!       identity
!


	!List of calling arguments
	integer ns
	real ( kind = 8 ), intent(in) :: P0(3),Pt(3)
	complex ( kind = 8 ), intent(in) :: vf(3)
	real ( kind = 8 ), intent(in) :: srcvals(12,ns)
	complex ( kind = 8 ), intent(in) :: zpars(5)
	complex ( kind = 8 ), intent(out) :: RHS_out(6*ns)
	
	!List of local variables
	complex ( kind = 8 ) omega, ep0,mu0,ep1,mu1
	complex ( kind = 8 ), allocatable :: E(:,:), H(:,:)
	complex ( kind = 8 ), allocatable :: RHS(:)

	integer count1,i
	complex ( kind = 8 ) ima,vf_m(3),vf2(3),vf2_m(3)
	real ( kind = 8 ) ru(3),rv(3),cross_aux(3)
	allocate(E(3,ns), H(3,ns), RHS(6*ns))

	ima=(0.0d0,1.0d0)
	omega = zpars(1)
	ep0 = zpars(2)
	mu0 = zpars(3)
	ep1 = zpars(4)
	mu1 = zpars(5)
	
	vf_m(1)=-1.0d0*vf(1)
	vf_m(2)=-1.0d0*vf(2)
	vf_m(3)=-1.0d0*vf(3)
	
	vf2(1)=ep0*vf(1)
	vf2(2)=ep0*vf(2)
	vf2(3)=ep0*vf(3)
	
	vf2_m(1)=-ep1*vf(1)
	vf2_m(2)=-ep1*vf(2)
	vf2_m(3)=-ep1*vf(3)
	
	call fieldsEDomega(omega,ep0,mu0,P0,srcvals,ns,E,H,vf,0)
	call fieldsMDomega(omega,ep0,mu0,P0,srcvals,ns,E,H,vf,1)

	call fieldsEDomega(omega,ep1,mu1,Pt,srcvals,ns,E,H,vf_m,1)
	call fieldsMDomega(omega,ep1,mu1,Pt,srcvals,ns,E,H,vf_m,1)
!	read (*,*)
	
	do count1=1,ns	
      call orthonormalize(srcvals(4:6,count1),srcvals(10:12,count1),&
	   &ru,rv)		
      RHS(count1)=-DOT_PRODUCT(rv,E(:,count1))
      RHS(ns+count1)=DOT_PRODUCT(ru,E(:,count1))
      RHS(2*ns+count1)=-DOT_PRODUCT(rv,H(:,count1))*ima*omega
	  RHS(3*ns+count1)=DOT_PRODUCT(ru,H(:,count1))*ima*omega
      RHS(4*ns+count1)=0.0d0
	enddo
	call fieldsEDomega(omega,ep0,mu0,P0,srcvals,ns,E,H,vf2,0)
	call fieldsMDomega(omega,ep0,mu0,P0,srcvals,ns,E,H,vf2,1)

	call fieldsEDomega(omega,ep1,mu1,Pt,srcvals,ns,E,H,vf2_m,1)
	call fieldsMDomega(omega,ep1,mu1,Pt,srcvals,ns,E,H,vf2_m,1)
	
	do count1=1,ns
      RHS(5*ns+count1)=DOT_PRODUCT(srcvals(10:12,count1),E(:,count1))
	enddo
	
	do count1=1,ns
      RHS(count1)=RHS(count1)/(mu0+mu1)
      RHS(count1+ns)=RHS(count1+ns)/(mu0+mu1)
      RHS(count1+2*ns)=RHS(count1+2*ns)/(ep0+ep1)
      RHS(count1+3*ns)=RHS(count1+3*ns)/(ep0+ep1)
      RHS(count1+4*ns)=RHS(count1+4*ns)/(mu0+mu1)
      RHS(count1+5*ns)=-RHS(count1+5*ns)/(ep0+ep1)
	enddo


    do i=1,ns
      RHS_out(i)=RHS(i)
      RHS_out(i+ns)=RHS(i+ns)
      RHS_out(i+2*ns)=RHS(i+4*ns)
      RHS_out(i+3*ns)=RHS(i+2*ns)
      RHS_out(i+4*ns)=RHS(i+3*ns)
      RHS_out(i+5*ns)=RHS(i+5*ns)
    enddo

return
end subroutine get_rhs_em_dfie_trans


subroutine test_accuracy_em_dfie_trans(eps_FMM,sol,zpars,ns,wts,&
 &srcvals,P0,vf,Pt)
implicit none
!
!  This function test the accuracy of the solution computed in the
!  exterior region
!
!  input:
!    eps_FMM - real *8
!      epsilon for the fmm call
!
!    sol - complex *16(6*ns)
!      solution of the DFIE equation
!	
!    zpars - complex *16(3)
!      kernel parameters
!      zpars(1) = omega 
!      zpars(2) = ep0
!      zpars(3) = mu0
!      zpars(4) = ep1
!      zpars(5) = mu1
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
	complex ( kind = 8 ), intent(in) :: zpars(5)
    integer, intent(in) :: ns
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),eps_FMM
    real ( kind = 8 ), intent(in) :: wts(ns),P0(3),Pt(3)
	complex ( kind = 8 ), intent(in) :: sol(6*ns),vf(3)
	
    !List of local variables
	complex ( kind = 8 ) ima, R1, R2,Et1(3),Et2(3),Ht1(3),Ht2(3)
	complex ( kind = 8 ) E01(3),E02(3),H01(3),H02(3)
	real ( kind = 8 ) ru_s(3),rv_s(3),n_s(3),sour(3),r,dr(3)
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3)
	real ( kind = 8 ) error_E,rel_err_E,error_H,rel_err_H
	real ( kind = 8 ) pi
	complex ( kind = 8 ) omega,ep0,mu0,ep1,mu1

	integer count1

	ima=(0.0d0,1.0d0)
	pi=3.1415926535897932384626433832795028841971d0
	omega = zpars(1)
	ep0 = zpars(2)
	mu0 = zpars(3)
	ep1 = zpars(4)
	mu1 = zpars(5)

	write (*,*) 'solution:' ,sol(1:10)
	
	call em_dfie_trans_FMM_targ(eps_FMM,omega,ep0,mu0,ns,srcvals,1,Pt,&
	 &wts,sol(1:ns),sol(ns+1:2*ns),sol(2*ns+1:3*ns),sol(3*ns+1:4*ns),&
	 &sol(4*ns+1:5*ns),sol(5*ns+1:6*ns),Et1)
	call em_dfie_trans_FMM_targ(eps_FMM,omega,ep1,mu1,ns,srcvals,1,P0,&
	 &wts,sol(1:ns),sol(ns+1:2*ns),sol(2*ns+1:3*ns),sol(3*ns+1:4*ns),&
	 &sol(4*ns+1:5*ns),sol(5*ns+1:6*ns),E01)
		
	call fieldsEDomega(omega,ep0,mu0,P0,Pt,1,Et2,Ht2,vf,0)
	call fieldsMDomega(omega,ep0,mu0,P0,Pt,1,Et2,Ht2,vf,1)
		
	call fieldsEDomega(omega,ep1,mu1,Pt,P0,1,E02,H02,vf,0)
	call fieldsMDomega(omega,ep1,mu1,Pt,P0,1,E02,H02,vf,1)
	
	write (*,*) 'Errors at the EXTERIOR region:'
	error_E=sqrt(abs(Et1(1)-Et2(1))**2+abs(Et1(2)-Et2(2))**2+&
	 &abs(Et1(3)-Et2(3))**2)
!	write (*,*) 'Error E: ', error_E
	write (*,*) 'Relative Error in E: ', error_E/sqrt(abs(Et2(1))**2+&
	 &abs(Et2(2))**2+abs(Et2(3))**2)
		
	write (*,*) 'Errors at the INTERIOR region:'
	error_E=sqrt(abs(E01(1)-E02(1))**2+abs(E01(2)-E02(2))**2+&
	 &abs(E01(3)-E02(3))**2)
!	write (*,*) 'Error E: ', error_E
	write (*,*) 'Relative Error in E: ', error_E/sqrt(abs(E02(1))**2+&
	 &abs(E02(2))**2+abs(E02(3))**2)

return
end subroutine test_accuracy_em_dfie_trans


subroutine em_dfie_trans_FMM_targ(eps,omega,ep,mu,ns,srcvals,nt,targ,&
 &wts,a_u,a_v,rho_in,b_u,b_v,lambda_in,E)
implicit none
!
!  This funciton computes the fields E0,E1 at a given point using the
!  FMM
!  It doesn't contain near field corrections (it's for debugging
!  purposes)
!
!  Representation:
!
!    E=mucurlS_{k}[a]-muS_{k}[n·sigma]+muepS_{k}[b]+gradS_{k}[rho]
!
!      where k=omega*sqrt(ep*mu)
!
!  input:
!    eps - real *8
!      epsilon for the fmm call
!
!    omega - complex *16
!      frequency (pulsation)
!
!    ep - complex *16
!      electric constant 
!
!    mu - complex *16
!      magnetic constant
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
!    rho_in - complex *16(ns)
!      iscalar source on the surface
!
!    b_u,b_v - complex *16(ns)
!      two components of the tangent induced current J on the surface
!      along srcvals(4:6,i) and (srcvals(10:12,i) x srcvals(4:6,i))
!      directions
! 
!    lambda_in - complex *16(ns)
!      scalar source on the surface
!
!
!  output:
!    E0 - complex  *16(3,nt)
!      value of the Electric field at the target points in the exterior
!      region
!
!    E1 - complex  *16(3,nt)
!      value of the Electric field at the target points in the interior
!      region
!

    !List of calling arguments
	real ( kind = 8 ), intent(in) :: eps
	complex ( kind = 8 ), intent(in) :: omega,ep,mu
    integer, intent(in) :: ns,nt
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),targ(3,nt)
    real ( kind = 8 ), intent(in) :: wts(ns)
    complex ( kind = 8 ), intent(in) :: a_u(ns),a_v(ns),rho_in(ns)
    complex ( kind = 8 ), intent(in) :: b_u(ns),b_v(ns),lambda_in(ns)
    complex ( kind = 8 ), intent(out) :: E(3,nt)

    !List of local variables
    real ( kind = 8 ), allocatable :: n_vect(:,:),u_vect(:,:)	
    real ( kind = 8 ), allocatable :: v_vect(:,:),source(:,:)
    complex ( kind = 8 ), allocatable :: a_vect(:,:),b_vect(:,:)
    complex ( kind = 8 ), allocatable :: lambda(:),rho(:)
    complex ( kind = 8 ), allocatable :: curlE(:,:),divE(:)
    complex ( kind = 8 ) ima,zk

    integer count1,count2
    integer ifa_vect,ifb_vect,iflambda,ifrho,ifE,ifcurlE,ifdivE

	ima=(0.0d0,1.0d0)
	zk=omega*sqrt(ep*mu)

    allocate(a_vect(3,ns))
    allocate(b_vect(3,ns))
    allocate(lambda(ns))
	allocate(rho(ns))
    allocate(curlE(3,nt))
	allocate(divE(nt))
	
	allocate(n_vect(3,ns))
	allocate(u_vect(3,ns))
	allocate(v_vect(3,ns))
	allocate(source(3,ns))

	do count1=1,ns
      n_vect(:,count1)=srcvals(10:12,count1)
      source(:,count1)=srcvals(1:3,count1)
	enddo
	call orthonormalize_all(srcvals(4:6,:),srcvals(10:12,:),&
	 &u_vect,v_vect,ns)


    do count1=1,ns
      a_vect(1,count1)=mu*(a_u(count1)*u_vect(1,count1)+a_v(count1)*&
	   &v_vect(1,count1))
      a_vect(2,count1)=mu*(a_u(count1)*u_vect(2,count1)+a_v(count1)*&
	   &v_vect(2,count1))
      a_vect(3,count1)=mu*(a_u(count1)*u_vect(3,count1)+a_v(count1)*&
	   &v_vect(3,count1))
				
      b_vect(1,count1)=mu*ep*(b_u(count1)*u_vect(1,count1)+b_v(count1)*&
	   &v_vect(1,count1))
      b_vect(2,count1)=mu*ep*(b_u(count1)*u_vect(2,count1)+b_v(count1)*&
	   &v_vect(2,count1))
      b_vect(3,count1)=mu*ep*(b_u(count1)*u_vect(3,count1)+b_v(count1)*&
	   &v_vect(3,count1))
		
      rho(count1)=-mu*rho_in(count1)
      lambda(count1)=lambda_in(count1)
	enddo

    !Computing the full operator
    ifa_vect=1
    ifb_vect=1
    iflambda=1
    ifrho=1
    ifE=1
    ifcurlE=1
    ifdivE=1

    call Vector_Helmholtz_targ(eps,zk,ns,source,wts,ifa_vect,a_vect,&
	 &ifb_vect,b_vect,iflambda,lambda,ifrho,rho,n_vect,ifE,E,ifcurlE,&
	 &curlE,ifdivE,divE,nt,targ)

	deallocate(a_vect)
	deallocate(b_vect)
	deallocate(lambda)
	deallocate(curlE)
	deallocate(rho)
	
	deallocate(u_vect)
	deallocate(v_vect)
	deallocate(n_vect)
	deallocate(source)
	deallocate(divE)

return
end subroutine em_dfie_trans_FMM_targ





subroutine em_dfie_trans(srcinfo, ndt,targinfo,ndd, dpars,ndz,zpars,&
 &ndi,ipars,E_val)
implicit none

! this function provides the near field kernel that will use
! zgetnearquad_ggq_guru 
! through getnearquad_DFIE

    !List of calling arguments
	integer, intent(in) :: ndt,ndd,ndz,ndi
	real ( kind = 8 ), intent(in) :: srcinfo(12)
	real ( kind = 8 ), intent(in) :: targinfo(ndt)
	integer, intent(in) :: ipars(ndi)
	real ( kind = 8 ), intent(in) :: dpars(ndd)
	complex ( kind = 8 ), intent(in) :: zpars(ndz)
	complex ( kind = 8 ), intent(out) :: E_val
	
	!List of local variables
	real ( kind = 8 ) ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ) ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ) du(3), dv(3),sour(3),targ(3)
	real ( kind = 8 ) r, dr(3),aux_real
    real ( kind = 8 ) c_aux,c_aux_A,c_aux_B,c_aux_C,c_aux_D
    real ( kind = 8 ) xprod_aux3(3),xprod_aux4(3)	
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3)
	complex ( kind = 8 ) R1_0,R1_1,R2_0,R2_1,ima	
	complex ( kind = 8 ) my_exp_0,my_exp_1,omega,ep0	
	complex ( kind = 8 ) mu0,ep1,mu1,zk0,zk1,z_aux_0,z_aux_1
	real ( kind = 8 ) pi
	integer count1,count2
	
	pi=3.1415926535897932384626433832795028841971d0
	ima=(0.0d0,1.0d0)
	omega=zpars(1)
	ep0=zpars(2)
	mu0=zpars(3)
	ep1=zpars(4)
	mu1=zpars(5)
	
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
	
	zk0=omega*sqrt(ep0*mu0)
	zk1=omega*sqrt(ep1*mu1)

!R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
!R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
	
	
!R2_0=((ima*zk0)**2/r**3-3.0d0*ima*zk0/r**4+3.0d0/r**5)*exp(ima*zk0*r)/&
!&(4.0d0*pi)
!R2_1=((ima*zk1)**2/r**3-3.0d0*ima*zk1/r**4+3.0d0/r**5)*exp(ima*zk1*r)/&
!&(4.0d0*pi)
	 
!my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
!my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
	
	call orthonormalize(srcinfo(4:6),n_s,ru_s,rv_s)
	call orthonormalize(targinfo(4:6),n_t,ru_t,rv_t)

    if (ipars(1).eq.1) then
      if (ipars(2).eq.1) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
	  	call my_cross_v2(dr,ru_s,xprod_aux1)		
		c_aux=-DOT_PRODUCT(xprod_aux1,rv_t)
		E_val=c_aux*(mu0*R1_0-mu1*R1_1)
!		write (*,*) 'dentro', E_val
!E_mat(1,1)=mu0*nxcurlSk0a(1,1)-mu1*nxcurlSk1a(1,1)
	  elseif (ipars(2).eq.2) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
    	call my_cross_v2(dr,rv_s,xprod_aux2)
	  	c_aux=-DOT_PRODUCT(xprod_aux2,rv_t)
		E_val=c_aux*(mu0*R1_0-mu1*R1_1)
!E_mat(1,2)=mu0*nxcurlSk0a(1,2)-mu1*nxcurlSk1a(1,2)
	  elseif (ipars(2).eq.3) then
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
	    call my_cross_v2(n_t,n_s,xprod_aux1)
	    c_aux=DOT_PRODUCT(ru_t,xprod_aux1)
	    E_val=c_aux/r*(-mu0*my_exp_0+mu1*my_exp_1)	  
!E_mat(1,3)=-mu0*nxSk0nrho(1,1)+mu1*nxSk1nrho(1,1)
	  elseif (ipars(2).eq.4) then
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        c_aux=DOT_PRODUCT(rv_t,ru_s)
        E_val=-c_aux/r*(ep0*mu0*my_exp_0-ep1*mu1*my_exp_1)
!E_mat(1,4)=ep0*mu0*nxSk0b(1,1)-ep1*mu1*nxSk1b(1,1)
	  elseif (ipars(2).eq.5) then
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        c_aux=DOT_PRODUCT(rv_t,rv_s)
        e_val=-c_aux/r*(ep0*mu0*my_exp_0-ep1*mu1*my_exp_1)
!E_mat(1,5)=ep0*mu0*nxSk0b(1,2)-ep1*mu1*nxSk1b(1,2)
	  else
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
        c_aux=DOT_PRODUCT(rv_t,dr)
        E_val=-c_aux*(R1_0-R1_1)
!E_mat(1,6)=nxgradSk0lambda(1,1)-nxgradSk1lambda(1,1)
	  endif
	elseif (ipars(1).eq.2) then
      if (ipars(2).eq.1) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
		call my_cross_v2(dr,ru_s,xprod_aux1)
	  	c_aux=DOT_PRODUCT(xprod_aux1,ru_t)
		E_val=c_aux*(mu0*R1_0-mu1*R1_1)
!E_mat(2,1)=mu0*nxcurlSk0a(2,1)-mu1*nxcurlSk1a(2,1)
	  elseif (ipars(2).eq.2) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
		call my_cross_v2(dr,rv_s,xprod_aux2)
	  	c_aux=DOT_PRODUCT(xprod_aux2,ru_t)
		E_val=c_aux*(mu0*R1_0-mu1*R1_1)	  
!E_mat(2,2)=mu0*nxcurlSk0a(2,2)-mu1*nxcurlSk1a(2,2)
	  elseif (ipars(2).eq.3) then
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
	  	call my_cross_v2(n_t,n_s,xprod_aux1)
	    c_aux=DOT_PRODUCT(rv_t,xprod_aux1)
	    E_val=c_aux/r*(-mu0*my_exp_0+mu1*my_exp_1)
!E_mat(2,3)=-mu0*nxSk0nrho(2,1)+mu1*nxSk1nrho(2,1)
	  elseif (ipars(2).eq.4) then
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        c_aux=DOT_PRODUCT(ru_t,ru_s)
	    E_val=c_aux/r*(ep0*mu0*my_exp_0-ep1*mu1*my_exp_1)
!E_mat(2,4)=ep0*mu0*nxSk0b(2,1)-ep1*mu1*nxSk1b(2,1)
	  elseif (ipars(2).eq.5) then
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
	  	c_aux=DOT_PRODUCT(ru_t,rv_s)
		E_val=c_aux/r*(ep0*mu0*my_exp_0-ep1*mu1*my_exp_1)
!E_mat(2,5)=ep0*mu0*nxSk0b(2,2)-ep1*mu1*nxSk1b(2,2)
	  else
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
        c_aux=DOT_PRODUCT(ru_t,dr)
        E_val=c_aux*(R1_0-R1_1)
!E_mat(2,6)=nxgradSk0lambda(2,1)-nxgradSk1lambda(2,1)
	  endif
	elseif (ipars(1).eq.3) then
      if (ipars(2).eq.1) then
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        R1_0=(ima*zk0*r-1.0d0)/r**3*my_exp_0
        R1_1=(ima*zk1*r-1.0d0)/r**3*my_exp_1
        R2_0=((ima*zk0)**2/r**3-3.0d0*ima*zk0/r**4+3.0d0/r**5)*&
		 &exp(ima*zk0*r)/(4.0d0*pi)
        R2_1=((ima*zk1)**2/r**3-3.0d0*ima*zk1/r**4+3.0d0/r**5)*&
		 &exp(ima*zk1*r)/(4.0d0*pi)
	    c_aux_A=DOT_PRODUCT(rv_t,ru_s)
	    c_aux_B=DOT_PRODUCT(rv_t,ru_s)
	    c_aux_C=DOT_PRODUCT(rv_t,dr)
	    c_aux_D=DOT_PRODUCT(ru_s,-dr)
	    z_aux_0=zk0**2*(-c_aux_A*my_exp_0/r)+(-c_aux_B*R1_0)-&
		 &(-c_aux_C*c_aux_D*R2_0)
	    z_aux_1=zk1**2*(-c_aux_A*my_exp_1/r)+(-c_aux_B*R1_1)-&
		 &(-c_aux_C*c_aux_D*R2_1)
	    E_val=z_aux_0-z_aux_1
!E_mat(3,1)=(nxcurlcurlSk0a(1,1)-nxcurlcurlSk1a(1,1))
	  elseif (ipars(2).eq.2) then
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        R1_0=(ima*zk0*r-1.0d0)/r**3*my_exp_0
        R1_1=(ima*zk1*r-1.0d0)/r**3*my_exp_1
        R2_0=((ima*zk0)**2/r**3-3.0d0*ima*zk0/r**4+3.0d0/r**5)*&
		 &exp(ima*zk0*r)/(4.0d0*pi)
        R2_1=((ima*zk1)**2/r**3-3.0d0*ima*zk1/r**4+3.0d0/r**5)*&
		 &exp(ima*zk1*r)/(4.0d0*pi)
	  	c_aux_A=DOT_PRODUCT(rv_t,rv_s)
	    c_aux_B=DOT_PRODUCT(rv_t,rv_s)
	    c_aux_C=DOT_PRODUCT(rv_t,dr)
	    c_aux_D=DOT_PRODUCT(rv_s,-dr)
	    z_aux_0=zk0**2*(-c_aux_A*my_exp_0/r)+(-c_aux_B*R1_0)-&
		 &(-c_aux_C*c_aux_D*R2_0)
        z_aux_1=zk1**2*(-c_aux_A*my_exp_1/r)+(-c_aux_B*R1_1)-&
		 &(-c_aux_C*c_aux_D*R2_1)
        E_val=z_aux_0-z_aux_1
!E_mat(3,2)=(nxcurlcurlSk0a(1,2)-nxcurlcurlSk1a(1,2))
	  elseif (ipars(2).eq.3) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
        call my_cross_v2(dr, n_s, xprod_aux1)
        c_aux=DOT_PRODUCT(rv_t,xprod_aux1)
        E_val=-c_aux*(-R1_0+R1_1)	
!E_mat(3,3)=-nxcurlSk0nrho(1,1)+nxcurlSk1nrho(1,1)
	  elseif (ipars(2).eq.4) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
		call my_cross_v2(dr,ru_s,xprod_aux1)
		c_aux=-DOT_PRODUCT(xprod_aux1,rv_t)
		E_val=c_aux*(ep0*R1_0-ep1*R1_1)
!E_mat(3,4)=ep0*nxcurlSk0a(1,1)-ep1*nxcurlSk1a(1,1)
	  elseif (ipars(2).eq.5) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
		call my_cross_v2(dr,rv_s,xprod_aux2)
	  	c_aux=-DOT_PRODUCT(xprod_aux2,rv_t)
		E_val=c_aux*(ep0*R1_0-ep1*R1_1)
!E_mat(3,5)=ep0*nxcurlSk0a(1,2)-ep1*nxcurlSk1a(1,2)
	  else
	    E_val=0.0d0
	  endif
	elseif (ipars(1).eq.4) then
      if (ipars(2).eq.1) then
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        R1_0=(ima*zk0*r-1.0d0)/r**3*my_exp_0
        R1_1=(ima*zk1*r-1.0d0)/r**3*my_exp_1
        R2_0=((ima*zk0)**2/r**3-3.0d0*ima*zk0/r**4+3.0d0/r**5)*&
		 &exp(ima*zk0*r)/(4.0d0*pi)
        R2_1=((ima*zk1)**2/r**3-3.0d0*ima*zk1/r**4+3.0d0/r**5)*&
		 &exp(ima*zk1*r)/(4.0d0*pi)
        c_aux_A=DOT_PRODUCT(ru_t,ru_s)
        c_aux_B=DOT_PRODUCT(ru_t,ru_s)
        c_aux_C=DOT_PRODUCT(ru_t,dr)
        c_aux_D=DOT_PRODUCT(ru_s,-dr)
        z_aux_0=zk0**2*(c_aux_A*my_exp_0/r)+(+c_aux_B*R1_0)-&
		 &(+c_aux_C*c_aux_D*R2_0)
        z_aux_1=zk1**2*(c_aux_A*my_exp_1/r)+(+c_aux_B*R1_1)-&
		 &(+c_aux_C*c_aux_D*R2_1)
        E_val=z_aux_0-z_aux_1
!E_mat(4,1)=(nxcurlcurlSk0a(2,1)-nxcurlcurlSk1a(2,1))
	  elseif (ipars(2).eq.2) then	
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        R1_0=(ima*zk0*r-1.0d0)/r**3*my_exp_0
        R1_1=(ima*zk1*r-1.0d0)/r**3*my_exp_1
        R2_0=((ima*zk0)**2/r**3-3.0d0*ima*zk0/r**4+3.0d0/r**5)*&
		 &exp(ima*zk0*r)/(4.0d0*pi)
        R2_1=((ima*zk1)**2/r**3-3.0d0*ima*zk1/r**4+3.0d0/r**5)*&
		 &exp(ima*zk1*r)/(4.0d0*pi)
	  	c_aux_A=DOT_PRODUCT(ru_t,rv_s)
	    c_aux_B=DOT_PRODUCT(ru_t,rv_s)
	    c_aux_C=DOT_PRODUCT(ru_t,dr)
	    c_aux_D=DOT_PRODUCT(rv_s,-dr)
	    z_aux_0=zk0**2*(c_aux_A*my_exp_0/r)+(+c_aux_B*R1_0)-&
		 &(+c_aux_C*c_aux_D*R2_0)
	    z_aux_1=zk1**2*(c_aux_A*my_exp_1/r)+(+c_aux_B*R1_1)-&
		 &(+c_aux_C*c_aux_D*R2_1)
        E_val=z_aux_0-z_aux_1
!E_mat(4,2)=(nxcurlcurlSk0a(2,2)-nxcurlcurlSk1a(2,2))
	  elseif (ipars(2).eq.3) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
        call my_cross_v2(dr, n_s, xprod_aux1)
        c_aux=DOT_PRODUCT(ru_t,xprod_aux1)
        E_val=c_aux*(-R1_0+R1_1)
!E_mat(4,3)=-nxcurlSk0nrho(2,1)+nxcurlSk1nrho(2,1)
	  elseif (ipars(2).eq.4) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
		call my_cross_v2(dr,ru_s,xprod_aux1)
	  	c_aux=DOT_PRODUCT(xprod_aux1,ru_t)
		E_val=c_aux*(ep0*R1_0-ep1*R1_1)
!E_mat(4,4)=ep0*nxcurlSk0a(2,1)-ep1*nxcurlSk1a(2,1)
	  elseif (ipars(2).eq.5) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
		call my_cross_v2(dr,rv_s,xprod_aux2)
	  	c_aux=DOT_PRODUCT(xprod_aux2,ru_t)
		E_val=c_aux*(ep0*R1_0-ep1*R1_1)	  
!E_mat(4,5)=ep0*nxcurlSk0a(2,2)-ep1*nxcurlSk1a(2,2)
	  else
	    E_val=0.0d0
	  endif
	elseif (ipars(1).eq.5) then
      if (ipars(2).eq.1) then
	    E_val=0.0d0
	  elseif (ipars(2).eq.2) then
        E_val=0.0d0
	  elseif (ipars(2).eq.3) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
        c_aux=DOT_PRODUCT(n_s,dr)
        E_val=-c_aux*(mu0*R1_0-mu1*R1_1)
!E_mat(5,3)=mu0*Dk0rho(1,1)-mu1*Dk1rho(1,1)
      elseif (ipars(2).eq.4) then	  
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
        c_aux=DOT_PRODUCT(dr,ru_s)
        E_val=c_aux*(ep0*mu0*R1_0-ep1*mu1*R1_1)
!E_mat(5,4)=ep0*mu0*divSk0b(1,1)-ep1*mu1*divSk1b(1,1)
	  elseif (ipars(2).eq.5) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
        c_aux=DOT_PRODUCT(dr,rv_s)
        E_val=c_aux*(ep0*mu0*R1_0-ep1*mu1*R1_1)
!E_mat(5,5)=ep0*mu0*divSk0b(1,2)-ep1*mu1*divSk1b(1,2)
	  else
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
	    E_val=(-zk0**2*my_exp_0+zk1**2*my_exp_1)/r
!E_mat(5,6)=-zk0**2*Sk0lambda(1,1)+zk1**2*Sk1lambda(1,1)
	  endif
	else
      if (ipars(2).eq.1) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
        call my_cross_v2(dr, ru_s, xprod_aux1)
        c_aux=DOT_PRODUCT(n_t,xprod_aux1)
        E_val=c_aux*(ep0*mu0*R1_0-ep1*mu1*R1_1)	  
!E_mat(6,1)=ep0*mu0*ncurlSk0a(1,1)-ep1*mu1*ncurlSk1a(1,1)
	  elseif (ipars(2).eq.2) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
        call my_cross_v2(dr, rv_s, xprod_aux2)
        c_aux=DOT_PRODUCT(n_t,xprod_aux2)
        E_val=c_aux*(ep0*mu0*R1_0-ep1*mu1*R1_1)
!E_mat(6,2)=ep0*mu0*ncurlSk0a(1,2)-ep1*mu1*ncurlSk1a(1,2)
	  elseif (ipars(2).eq.3) then
        c_aux=DOT_PRODUCT(n_t,n_s)
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        E_val=c_aux/r*(-ep0*mu0*my_exp_0+ep1*mu1*my_exp_1)
!E_mat(6,3)=-ep0*mu0*nSk0nrho(1,1)+ep1*mu1*nSk1nrho(1,1)
	  elseif (ipars(2).eq.4) then
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        c_aux=DOT_PRODUCT(n_t,ru_s)
        E_val=c_aux/r*(mu0*ep0**2*my_exp_0-mu1*ep1**2*my_exp_1)
!E_mat(6,4)=mu0*ep0**2*nSk0b(1,1)-mu1*ep1**2*nSk1b(1,1)
	  elseif (ipars(2).eq.5) then
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        c_aux=DOT_PRODUCT(n_t,rv_s)
        E_val=c_aux/r*(mu0*ep0**2*my_exp_0-mu1*ep1**2*my_exp_1)
!E_mat(6,5)=mu0*ep0**2*nSk0b(1,2)-mu1*ep1**2*nSk1b(1,2)
	  else
        c_aux=DOT_PRODUCT(n_t,dr)
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
        E_val=c_aux*(ep0*R1_0-ep1*R1_1)
!E_mat(6,6)=ep0*ngradSk0lambda(1,1)-ep1*ngradSk1lambda(1,1)
	  endif
	endif

return
end subroutine em_dfie_trans
