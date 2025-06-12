      subroutine getnearquad_em_muller_trans(npatches,norders,&
     &ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,&
     &iquad,rfac0,nquad,wnear)
!
!
!  This subroutine generates the near field quadrature
!  for the Muller integral equation:
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
!           wnear - complex *16(16*nquad)
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
      complex *16 wnear(16*nquad)
	  
      integer ipatch_id(ntarg)
      real *8 uvs_targ(2,ntarg)

      complex *16 alpha,beta
      integer i,j,ndi,ndd,ndz,count1,count2,icount

      integer ipv
      integer :: t1, t2,clock_rate, clock_max

      procedure (), pointer :: fker
      external  fker_em_muller_trans
      external  em_muller_trans

     ndz=5
	   ndd=1
	   ndi=2
	   ipv=1

!!  fker => fker_em_muller_trans
	  fker => em_muller_trans
!    call system_clock ( t1, clock_rate, clock_max )


!	  icount=0
!	  do count1=1,4
!	    do count2=1,4
!		  ipars(1)=count1
!      ipars(2)=count2
!      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
!		   &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
!		   &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
!		   &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
!		   &wnear(icount*nquad+1:(icount+1)*nquad))
!		   icount=icount+1
!        enddo
!	  enddo



		  ipars(1)=1
      ipars(2)=1
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
		   &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
		   &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
		   &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
		   &wnear(1:nquad))

		  ipars(1)=1
      ipars(2)=2
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
		   &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
		   &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
		   &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
		   &wnear(nquad+1:2*nquad))

		  ipars(1)=1
      ipars(2)=3
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
		   &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
		   &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
		   &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
		   &wnear(2*nquad+1:3*nquad))

		  ipars(1)=1
      ipars(2)=4
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
		   &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
		   &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
		   &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
		   &wnear(3*nquad+1:4*nquad))


		  ipars(1)=2
      ipars(2)=1
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
		   &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
		   &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
		   &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
		   &wnear(4*nquad+1:5*nquad))

		  ipars(1)=2
      ipars(2)=2
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
		   &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
		   &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
		   &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
		   &wnear(5*nquad+1:6*nquad))

		  ipars(1)=2
      ipars(2)=3
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
		   &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
		   &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
		   &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
		   &wnear(6*nquad+1:7*nquad))

		  ipars(1)=2
      ipars(2)=4
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
		   &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
		   &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
		   &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
		   &wnear(7*nquad+1:8*nquad))

		  ipars(1)=3
      ipars(2)=1
       wnear(8*nquad+1:9*nquad)=-wnear(2*nquad+1:3*nquad)

		  ipars(1)=3
      ipars(2)=2
       wnear(9*nquad+1:10*nquad)=-wnear(3*nquad+1:4*nquad)

		  ipars(1)=3
      ipars(2)=3
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
		   &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
		   &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
		   &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
		   &wnear(10*nquad+1:11*nquad))

		  ipars(1)=3
      ipars(2)=4
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
		   &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
		   &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
		   &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
		   &wnear(11*nquad+1:12*nquad))

		  ipars(1)=4
      ipars(2)=1
       wnear(12*nquad+1:13*nquad)=-wnear(6*nquad+1:7*nquad)

		  ipars(1)=4
      ipars(2)=2
       wnear(13*nquad+1:14*nquad)=-wnear(7*nquad+1:8*nquad)

		  ipars(1)=4
      ipars(2)=3
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
		   &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
		   &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
		   &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
		   &wnear(14*nquad+1:15*nquad))

		  ipars(1)=4
      ipars(2)=4
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
		   &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
		   &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
		   &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
		   &wnear(15*nquad+1:16*nquad))


!call system_clock ( t2, clock_rate, clock_max )
!write (*,*) 'time quadratures: ',real(t2-t1)/real(clock_rate)
!stop
    return
    end subroutine getnearquad_em_muller_trans


      subroutine lpcomp_em_muller_trans_addsub(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,sigma,novers,&
     &nptso,ixyzso,srcover,whtsover,pot,wnear)

!
!  This subroutine evaluates the layer potential for
!  the Muller boundary integral equation:
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
!    wnear  - complex *16(16*nquad)
!      near field precomputed quadrature
!
!    sigma - complex *16(2*ns)
!      induced charge and current on the surface
!      sigma(1:ns) - first component of 'a' along
!        the srcvals(4:6,i) direction
!      sigma(ns+1:2*ns) - second component of 'a' along
!        the (srcvals(10:12,i) x srcvals(4:6,i)) direction
!      sigma(2*ns+1:3*ns) - first component of 'b' along
!        the srcvals(4:6,i) direction
!      sigma(3*ns+1:4*ns) - second component of 'b' along
!        the (srcvals(10:12,i) x srcvals(4:6,i)) direction
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
      complex *16 sigma(4*npts),sigma2(npts)
	  
      complex *16 wnear(16*nquad)
	  
      integer novers(npatches+1)
      integer nover,npolso,nptso
      real *8 srcover(12,nptso),whtsover(nptso)
      complex *16 pot(4*ntarg)
      complex *16, allocatable :: potsort(:)

      real *8, allocatable :: sources(:,:),targvals(:,:),wtmp2(:)
      complex *16, allocatable :: charges(:),dipvec(:,:),sigmaover(:)
      integer ns,nt
      complex *16 alpha,beta
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      complex *16 tmp(10),val,E(4)

      real *8 xmin,xmax,ymin,ymax,zmin,zmax,sizey,sizez,boxsize


      integer i,j,jpatch,jquadstart,jstart,count1,count2


      integer ifaddsub,ifdir

      integer ntj
      
      complex *16 zdotu,pottmp
      complex *16, allocatable :: ctmp2_s(:),dtmp2(:,:)
      complex *16, allocatable :: ctmp2_u(:),ctmp2_v(:)
      complex *16, allocatable :: ctmp2_a_u(:),ctmp2_a_v(:)
      complex *16, allocatable :: ctmp2_b_u(:),ctmp2_b_v(:)

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
      allocate(sigmaover(4*ns))
	    allocate(pot_aux(4*ntarg))

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

      ra = 0


!
!       fmm
!

		call get_fmm_thresh(12,ns,srcover,ndtarg,ntarg,targs,thresh)
       ifdir=0
	  
		!Calculate the far_field with FMM		
    call em_muller_trans_FMM2(eps,zpars,ns,npts,srcover,targs,whtsover,&
    &sigmaover(1:ns),sigmaover(ns+1:2*ns),sigmaover(2*ns+1:3*ns),&
    &sigmaover(3*ns+1:4*ns),pot_aux(1:ntarg),pot_aux(ntarg+1:2*ntarg),&
    &pot_aux(2*ntarg+1:3*ntarg),pot_aux(3*ntarg+1:4*ntarg),&
    &thresh,ifdir)

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
		        do count1=0,3
			       do count2=0,3
			         pot_aux(i+count1*npts) = pot_aux(i+count1*npts) + &
				      &wnear((count1*4+count2)*nquad+jquadstart+l-1)*&
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
			      ctmp2_b_u(ii)=sigmaover(jstart+l+2*ns)
			      ctmp2_b_v(ii)=sigmaover(jstart+l+3*ns)
			      wtmp2(ii)=whtsover(jstart+l)
          enddo
        enddo
        call em_muller_trans_FMM2(eps,zpars,nss,ntarg0,srctmp2,targs(:,i),&
        &wtmp2,ctmp2_a_u,ctmp2_a_v,ctmp2_b_u,ctmp2_b_v,&
        &E(1),E(2),E(3),E(4),thresh,ifdir)
	 
        do j=0,3
          pot_aux(i+j*ntarg) = pot_aux(i+j*ntarg) - E(j+1)
        enddo
		
        deallocate(srctmp2,wtmp2)
        deallocate(ctmp2_a_u,ctmp2_a_v,ctmp2_b_u,ctmp2_b_v)

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
	  enddo
      
	  
	  do i=1,ntarg
		  pot(i)=pot_aux(i)
		  pot(i+ntarg)=pot_aux(i+ntarg)
		  pot(i+2*ntarg)=pot_aux(i+2*ntarg)
		  pot(i+3*ntarg)=pot_aux(i+3*ntarg)
	  enddo
	  
	  return
    end subroutine lpcomp_em_muller_trans_addsub



subroutine em_muller_trans_solver(npatches, norders, ixyzs, &
      iptype, npts, srccoefs, srcvals, eps, zpars, numit, &
      rhs, eps_gmres, niter, errs, rres, soln)
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
!  Representation: (1), Standard Muller
!
!    E0 = mu0curlS_{k0}[a] + curlcurlS_{k0}[b]/(-i\omega)
!    E1 = mu1curlS_{k1}[a] + curlcurlS_{k1}[b]/(-i\omega)
!
!    H0 = ep0curlS_{k0}[b] + curlcurlS_{k0}[a]/(-i\omega)
!    H1 = ep1curlS_{k1}[b] + curlcurlS_{k1}[a]/(-i\omega)

!
!  Boundary conditions: (2), Standard EM transmission
!
!                    nxE0-nxE1 = -nxE_inc
!                    nxH0-nxH1 = -nxH_inc
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
!    rhs - complex *16(4*npts)
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
!      soln(2*npts+1:3*npts) component of the tangent induced current
!        a on the surface along srcvals(4:6,i) direction
!      soln(3*npts+1:4*npts) component of the tangent induced current
!        a on the surface along (srcvals(10:12,i) x srcvals(4:6,i))
!        direction
!

      implicit none
      integer npatches,norder,npols,npts
      integer norders(npatches),ixyzs(npatches+1)
      integer iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,eps_gmres
      complex *16 zpars(5)
      complex *16 rhs(4*npts)
      complex *16 soln(4*npts)

      real *8, allocatable :: targs(:,:)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      integer ndtarg,ntarg

      real *8 errs(numit+1)
      real *8 rres,eps2
      integer niter, numit


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
      complex *16 ima

      ima=(0.0d0,1.0d0)
  
!
!   n_var is the number of unknowns in the linear system.
!   as we have two vector unknown a and b 
!   we need n_var=4*npts
!

      n_var=4*npts


      done = 1
      pi = atan(done)*4


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

      ikerorder = 0

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
      allocate(wnear(16*nquad))
	  
!C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,16*nquad
	    wnear(i)=0
      enddo
!C$OMP END PARALLEL DO    

      iquadtype = 1

!      goto 1111


      print *, "starting to generate near quadrature"
      call cpu_time(t1)
!C$      t1 = omp_get_wtime()      

       call getnearquad_em_muller_trans(npatches,norders,&
     &ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,&
     &ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,&
     &iquad,rfac0,nquad,wnear)
      call cpu_time(t2)
!C$      t2 = omp_get_wtime()     

      call prin2('quadrature generation time=*',t2-t1,1)
     
 1111 continue     
      print *, "done generating near quadrature, now starting gmres"

      call em_muller_trans_solver_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, numit, rhs, &
        nnz, row_ptr, col_ind, iquad, nquad, wnear, novers, npts_over, &
        ixyzso, srcover, wover, eps_gmres, niter, errs, rres, soln)

      return
      end
!
!
!
!
!

      subroutine em_muller_trans_solver_guru(npatches, norders, ixyzs, &
        iptype, npts, srccoefs, srcvals, eps, zpars, numit, rhs, &
        nnz, row_ptr, col_ind, iquad, nquad, wnear, novers, npts_over, &
        ixyzso, srcover, wover, eps_gmres, niter, errs, rres, soln)
!
!  Guru subroutine for em_muller_trans_solver.
!
!  Documentation coming soon
!
!
      implicit none

      integer npatches,norder,npols,npts
      integer norders(npatches),ixyzs(npatches+1)
      integer iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,eps_gmres
      complex *16 zpars(5)
      complex *16 rhs(4*npts)
      complex *16 soln(4*npts)

      real *8, allocatable :: targs(:,:)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      integer ndtarg,ntarg

      real *8 errs(numit+1)
      real *8 rres,eps2
      integer niter


      integer nover,npolso,nptso
      integer nnz,nquad
      integer row_ptr(npts+1), col_ind(nnz), iquad(nnz+1)

      complex *16 wnear(16*nquad)

      real *8 srcover(12,npts_over),wover(npts_over)
      integer ixyzso(npatches+1),novers(npatches)


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

!
!   n_var is the number of unknowns in the linear system.
!   as we have two vector unknown a and b 
!   we need n_var=4*npts
!

      n_var=4*npts

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
      call get_patch_id_uvs(npatches, norders, ixyzs, iptype, npts, &
        ipatch_id, uvs_targ)




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


        call lpcomp_em_muller_trans_addsub(npatches,norders,ixyzs,&
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
!        replace this routine by appropriate layer potential
!        evaluation routine  
!


          call lpcomp_em_muller_trans_addsub(npatches,norders,ixyzs,&
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
      end subroutine em_muller_trans_solver_guru




subroutine fker_em_muller_trans(srcinfo, ndt,targinfo,ndd, dpars,ndz,&
 &zpars,ndi,ipars,E_val)
implicit none

!  THIS FUNCTION IS OBSOLETE, NOT USED, subroutine em_dfie_trans is ued
!  instead (much faster), but very hard to read..
!  this function provides the near field kernel that will use 
!  zgetnearquad_ggq_guru 
!  through getnearquad_Muller

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


		E_mat(1,1)=mu0*nxcurlSk0a(1,1)-mu1*nxcurlSk1a(1,1)
		E_mat(1,2)=mu0*nxcurlSk0a(1,2)-mu1*nxcurlSk1a(1,2)
		E_mat(2,1)=mu0*nxcurlSk0a(2,1)-mu1*nxcurlSk1a(2,1)
		E_mat(2,2)=mu0*nxcurlSk0a(2,2)-mu1*nxcurlSk1a(2,2)

		E_mat(3,3)=ep0*nxcurlSk0a(1,1)-ep1*nxcurlSk1a(1,1)
		E_mat(3,4)=ep0*nxcurlSk0a(1,2)-ep1*nxcurlSk1a(1,2)
		E_mat(4,3)=ep0*nxcurlSk0a(2,1)-ep1*nxcurlSk1a(2,1)
		E_mat(4,4)=ep0*nxcurlSk0a(2,2)-ep1*nxcurlSk1a(2,2)

		E_mat(1,3)=(nxcurlcurlSk0a(1,1)-nxcurlcurlSk1a(1,1))/(-ima*omega)
		E_mat(1,4)=(nxcurlcurlSk0a(1,2)-nxcurlcurlSk1a(1,2))/(-ima*omega)
		E_mat(2,3)=(nxcurlcurlSk0a(2,1)-nxcurlcurlSk1a(2,1))/(-ima*omega)
		E_mat(2,4)=(nxcurlcurlSk0a(2,2)-nxcurlcurlSk1a(2,2))/(-ima*omega)
		
		E_mat(3,1)=-E_mat(1,3)
		E_mat(3,2)=-E_mat(1,4)
		E_mat(4,1)=-E_mat(2,3)
		E_mat(4,2)=-E_mat(2,4)
		
	E_val=E_mat(ipars(1),ipars(2))


return
end subroutine fker_em_muller_trans

subroutine em_muller_trans(srcinfo,ndt,targinfo,ndd,dpars,ndz,zpars,&
 &ndi,ipars,E_val)
implicit none

! this function provides the near field kernel that will use
! zgetnearquad_ggq_guru 
! through getnearquad_Muller

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
	
!  ep0=targinfo(13)+ima*targinfo(14)
!  mu0=targinfo(15)+ima*targinfo(16)
!	ep1=targinfo(17)+ima*targinfo(18)
!  mu1=targinfo(19)+ima*targinfo(20)


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
!E_mat(1,1)=mu0*nxcurlSk0a(1,1)-mu1*nxcurlSk1a(1,1)      OK
	  elseif (ipars(2).eq.2) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
    	call my_cross_v2(dr,rv_s,xprod_aux2)
	  	c_aux=-DOT_PRODUCT(xprod_aux2,rv_t)
		E_val=c_aux*(mu0*R1_0-mu1*R1_1)
!E_mat(1,2)=mu0*nxcurlSk0a(1,2)-mu1*nxcurlSk1a(1,2)       OK
	  elseif (ipars(2).eq.3) then
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
	    E_val=(z_aux_0-z_aux_1)/(-ima*omega)
!E_mat(1,3)=(nxcurlcurlSk0a(1,1)-nxcurlcurlSk1a(1,1))/(-ima*omega)     OK
	  else
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
        E_val=(z_aux_0-z_aux_1)/(-ima*omega)
!E_mat(1,4)=(nxcurlcurlSk0a(1,2)-nxcurlcurlSk1a(1,2))/(-ima*omega)   OK
	  endif
	elseif (ipars(1).eq.2) then
      if (ipars(2).eq.1) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
		call my_cross_v2(dr,ru_s,xprod_aux1)
	  	c_aux=DOT_PRODUCT(xprod_aux1,ru_t)
		E_val=c_aux*(mu0*R1_0-mu1*R1_1)
!E_mat(2,1)=mu0*nxcurlSk0a(2,1)-mu1*nxcurlSk1a(2,1)     OK
	  elseif (ipars(2).eq.2) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
		call my_cross_v2(dr,rv_s,xprod_aux2)
	  	c_aux=DOT_PRODUCT(xprod_aux2,ru_t)
		E_val=c_aux*(mu0*R1_0-mu1*R1_1)	  
!E_mat(2,2)=mu0*nxcurlSk0a(2,2)-mu1*nxcurlSk1a(2,2)     OK
	  elseif (ipars(2).eq.3) then
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
        E_val=(z_aux_0-z_aux_1)/(-ima*omega)
!E_mat(2,3)=(nxcurlcurlSk0a(2,1)-nxcurlcurlSk1a(2,1))/(-ima*omega)   OK
	  else
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
        E_val=(z_aux_0-z_aux_1)/(-ima*omega)
!E_mat(2,4)=(nxcurlcurlSk0a(2,2)-nxcurlcurlSk1a(2,2))/(-ima*omega)
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
	    E_val=-(z_aux_0-z_aux_1)/(-ima*omega)
!E_mat(3,1)=-(nxcurlcurlSk0a(1,1)-nxcurlcurlSk1a(1,1))/(-ima*omega)     OK
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
        E_val=-(z_aux_0-z_aux_1)/(-ima*omega)
!E_mat(3,2)=-(nxcurlcurlSk0a(1,2)-nxcurlcurlSk1a(1,2))/(-ima*omega)   OK
	  elseif (ipars(2).eq.3) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
		    call my_cross_v2(dr,ru_s,xprod_aux1)
		    c_aux=-DOT_PRODUCT(xprod_aux1,rv_t)
		    E_val=c_aux*(ep0*R1_0-ep1*R1_1)
!E_mat(3,3)=ep0*nxcurlSk0a(1,1)-ep1*nxcurlSk1a(1,1)
	  else
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
		    call my_cross_v2(dr,rv_s,xprod_aux2)
	      c_aux=-DOT_PRODUCT(xprod_aux2,rv_t)
		    E_val=c_aux*(ep0*R1_0-ep1*R1_1)
!E_mat(3,4)=ep0*nxcurlSk0a(1,2)-ep1*nxcurlSk1a(1,2)
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
        E_val=-(z_aux_0-z_aux_1)/(-ima*omega)
!E_mat(4,1)=-(nxcurlcurlSk0a(2,1)-nxcurlcurlSk1a(2,1))/(-ima*omega)   OK
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
        E_val=-(z_aux_0-z_aux_1)/(-ima*omega)
!E_mat(4,2)=-(nxcurlcurlSk0a(2,2)-nxcurlcurlSk1a(2,2))/(-ima*omega)
	  elseif (ipars(2).eq.3) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
		call my_cross_v2(dr,ru_s,xprod_aux1)
	  	c_aux=DOT_PRODUCT(xprod_aux1,ru_t)
		E_val=c_aux*(ep0*R1_0-ep1*R1_1)
!E_mat(4,3)=ep0*nxcurlSk0a(2,1)-ep1*nxcurlSk1a(2,1)
	  else
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
		call my_cross_v2(dr,rv_s,xprod_aux2)
	  c_aux=DOT_PRODUCT(xprod_aux2,ru_t)
		E_val=c_aux*(ep0*R1_0-ep1*R1_1)	  
!E_mat(4,4)=ep0*nxcurlSk0a(2,2)-ep1*nxcurlSk1a(2,2)
	  endif
	endif
return
end subroutine em_muller_trans


subroutine em_muller_trans_FMM(eps,zpars,ns,nt,srcvals,targvals,wts,&
 &a_u,a_v,b_u,b_v,AA_u,AA_v,BB_u,BB_v,thresh,ifdir)
implicit none

!Not under use now, we use em_muller_trans_FMM2 because it uses one single
! helmholtz fmm calls per wave number with 6 scalar sources

    !List of calling arguments
    real ( kind = 8 ), intent(in) :: eps
    complex ( kind = 8 ), intent(in) :: zpars(5)
    integer, intent(in) :: ns,nt
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),wts(ns)
    real ( kind = 8 ), intent(in) :: targvals(12,nt)
    complex ( kind = 8 ), intent(in) :: a_u(ns),a_v(ns)
    complex ( kind = 8 ), intent(in) :: b_u(ns),b_v(ns)
    complex ( kind = 8 ), intent(out) :: AA_u(nt),AA_v(nt)
    complex ( kind = 8 ), intent(out) :: BB_u(nt),BB_v(nt)
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
		  a_vect(1,count1)=(b_u(count1)*u_vect_s(1,count1)+b_v(count1)*v_vect_s(1,count1))/(-ima*omega)
		  a_vect(2,count1)=(b_u(count1)*u_vect_s(2,count1)+b_v(count1)*v_vect_s(2,count1))/(-ima*omega)
		  a_vect(3,count1)=(b_u(count1)*u_vect_s(3,count1)+b_v(count1)*v_vect_s(3,count1))/(-ima*omega)
				
      b_vect(1,count1)=(a_u(count1)*u_vect_s(1,count1)+a_v(count1)*v_vect_s(1,count1))*mu0
      b_vect(2,count1)=(a_u(count1)*u_vect_s(2,count1)+a_v(count1)*v_vect_s(2,count1))*mu0
      b_vect(3,count1)=(a_u(count1)*u_vect_s(3,count1)+a_v(count1)*v_vect_s(3,count1))*mu0
	  enddo

    !Computing the full operator
    ifa_vect=1
    ifb_vect=1
    iflambda=0
    ifrho=0
    ifE=0
    ifcurlE=1
    ifdivE=0

	call Vector_Helmholtz_targ2(eps,zk0,ns,source,wts,ifa_vect,a_vect,&
	 &ifb_vect,b_vect,iflambda,lambda,ifrho,rho,n_vect_s,ifE,E,ifcurlE,&
	 &curlE,ifdivE,divE,nt,targets,thresh,ifdir)

!	call Vector_Helmholtz_targ(eps,izk,ns,source,wts,ifa_vect,a_vect,ifb_vect,&
!	 &b_vect,iflambda,lambda,ifrho,rho,n_vect_s,ifE,E,ifcurlE,curlE,ifdivE,divE,nt,targets)

    do count1=1,nt
        b_vect(1,count1)=n_vect_t(2,count1)*curlE(3,count1)-n_vect_t(3,count1)*curlE(2,count1)
        b_vect(2,count1)=n_vect_t(3,count1)*curlE(1,count1)-n_vect_t(1,count1)*curlE(3,count1)
        b_vect(3,count1)=n_vect_t(1,count1)*curlE(2,count1)-n_vect_t(2,count1)*curlE(1,count1)
    enddo

    do count1=1,nt
        AA_u(count1)=b_vect(1,count1)*u_vect_t(1,count1)+b_vect(2,count1)*u_vect_t(2,count1)&
		   &+b_vect(3,count1)*u_vect_t(3,count1)
        AA_v(count1)=b_vect(1,count1)*v_vect_t(1,count1)+b_vect(2,count1)*v_vect_t(2,count1)&
		   &+b_vect(3,count1)*v_vect_t(3,count1)
    enddo

	  do count1=1,ns
		    b_vect(1,count1)=(b_u(count1)*u_vect_s(1,count1)+b_v(count1)*v_vect_s(1,count1))*ep0
		    b_vect(2,count1)=(b_u(count1)*u_vect_s(2,count1)+b_v(count1)*v_vect_s(2,count1))*ep0
		    b_vect(3,count1)=(b_u(count1)*u_vect_s(3,count1)+b_v(count1)*v_vect_s(3,count1))*ep0
				
        a_vect(1,count1)=(a_u(count1)*u_vect_s(1,count1)+a_v(count1)*v_vect_s(1,count1))/(ima*omega)
        a_vect(2,count1)=(a_u(count1)*u_vect_s(2,count1)+a_v(count1)*v_vect_s(2,count1))/(ima*omega)
        a_vect(3,count1)=(a_u(count1)*u_vect_s(3,count1)+a_v(count1)*v_vect_s(3,count1))/(ima*omega)
	  enddo

    !Computing the full operator
    ifa_vect=1
    ifb_vect=1
    iflambda=0
    ifrho=0
    ifE=0
    ifcurlE=1
    ifdivE=0

	  call Vector_Helmholtz_targ2(eps,zk0,ns,source,wts,ifa_vect,a_vect,&
	   &ifb_vect,b_vect,iflambda,lambda,ifrho,rho,n_vect_s,ifE,E,ifcurlE,&
	   &curlE,ifdivE,divE,nt,targets,thresh,ifdir)

    do count1=1,nt
        b_vect(1,count1)=n_vect_t(2,count1)*curlE(3,count1)-n_vect_t(3,count1)*curlE(2,count1)
        b_vect(2,count1)=n_vect_t(3,count1)*curlE(1,count1)-n_vect_t(1,count1)*curlE(3,count1)
        b_vect(3,count1)=n_vect_t(1,count1)*curlE(2,count1)-n_vect_t(2,count1)*curlE(1,count1)
    enddo

    do count1=1,nt
        BB_u(count1)=b_vect(1,count1)*u_vect_t(1,count1)+b_vect(2,count1)*u_vect_t(2,count1)&
         &+b_vect(3,count1)*u_vect_t(3,count1)
        BB_v(count1)=b_vect(1,count1)*v_vect_t(1,count1)+b_vect(2,count1)*v_vect_t(2,count1)&
         &+b_vect(3,count1)*v_vect_t(3,count1)
    enddo


!!! now zk1


    do count1=1,ns
		    a_vect(1,count1)=(b_u(count1)*u_vect_s(1,count1)+b_v(count1)*v_vect_s(1,count1))/(-ima*omega)
		    a_vect(2,count1)=(b_u(count1)*u_vect_s(2,count1)+b_v(count1)*v_vect_s(2,count1))/(-ima*omega)
		    a_vect(3,count1)=(b_u(count1)*u_vect_s(3,count1)+b_v(count1)*v_vect_s(3,count1))/(-ima*omega)
				
        b_vect(1,count1)=(a_u(count1)*u_vect_s(1,count1)+a_v(count1)*v_vect_s(1,count1))*mu1
        b_vect(2,count1)=(a_u(count1)*u_vect_s(2,count1)+a_v(count1)*v_vect_s(2,count1))*mu1
        b_vect(3,count1)=(a_u(count1)*u_vect_s(3,count1)+a_v(count1)*v_vect_s(3,count1))*mu1
    enddo

    ifa_vect=1
    ifb_vect=1
    iflambda=0
    ifrho=0
    ifE=0
    ifcurlE=1
    ifdivE=0

	  call Vector_Helmholtz_targ2(eps,zk1,ns,source,wts,ifa_vect,a_vect,&
	   &ifb_vect,b_vect,iflambda,lambda,ifrho,rho,n_vect_s,ifE,E,ifcurlE,&
	   &curlE,ifdivE,divE,nt,targets,thresh,ifdir)


    do count1=1,nt
        b_vect(1,count1)=n_vect_t(2,count1)*curlE(3,count1)-n_vect_t(3,count1)*curlE(2,count1)
        b_vect(2,count1)=n_vect_t(3,count1)*curlE(1,count1)-n_vect_t(1,count1)*curlE(3,count1)
        b_vect(3,count1)=n_vect_t(1,count1)*curlE(2,count1)-n_vect_t(2,count1)*curlE(1,count1)
    enddo

    do count1=1,nt
        AA_u(count1)=AA_u(count1)-(b_vect(1,count1)*u_vect_t(1,count1)+b_vect(2,count1)*u_vect_t(2,count1)&
		     &+b_vect(3,count1)*u_vect_t(3,count1))
        AA_v(count1)=AA_v(count1)-(b_vect(1,count1)*v_vect_t(1,count1)+b_vect(2,count1)*v_vect_t(2,count1)&
		     &+b_vect(3,count1)*v_vect_t(3,count1))
    enddo


	  do count1=1,ns
		    b_vect(1,count1)=(b_u(count1)*u_vect_s(1,count1)+b_v(count1)*v_vect_s(1,count1))*ep1
		    b_vect(2,count1)=(b_u(count1)*u_vect_s(2,count1)+b_v(count1)*v_vect_s(2,count1))*ep1
		    b_vect(3,count1)=(b_u(count1)*u_vect_s(3,count1)+b_v(count1)*v_vect_s(3,count1))*ep1
				
        a_vect(1,count1)=(a_u(count1)*u_vect_s(1,count1)+a_v(count1)*v_vect_s(1,count1))/(ima*omega)
        a_vect(2,count1)=(a_u(count1)*u_vect_s(2,count1)+a_v(count1)*v_vect_s(2,count1))/(ima*omega)
        a_vect(3,count1)=(a_u(count1)*u_vect_s(3,count1)+a_v(count1)*v_vect_s(3,count1))/(ima*omega)
	  enddo


    ifa_vect=1
    ifb_vect=1
    iflambda=0
    ifrho=0
    ifE=0
    ifcurlE=1
    ifdivE=0

	  call Vector_Helmholtz_targ2(eps,zk1,ns,source,wts,ifa_vect,a_vect,&
	   &ifb_vect,b_vect,iflambda,lambda,ifrho,rho,n_vect_s,ifE,E,ifcurlE,&
	   &curlE,ifdivE,divE,nt,targets,thresh,ifdir)

    do count1=1,nt
        b_vect(1,count1)=n_vect_t(2,count1)*curlE(3,count1)-n_vect_t(3,count1)*curlE(2,count1)
        b_vect(2,count1)=n_vect_t(3,count1)*curlE(1,count1)-n_vect_t(1,count1)*curlE(3,count1)
        b_vect(3,count1)=n_vect_t(1,count1)*curlE(2,count1)-n_vect_t(2,count1)*curlE(1,count1)
    enddo

    do count1=1,nt
        BB_u(count1)=BB_u(count1)-(b_vect(1,count1)*u_vect_t(1,count1)+b_vect(2,count1)*u_vect_t(2,count1)&
		     &+b_vect(3,count1)*u_vect_t(3,count1))
        BB_v(count1)=BB_v(count1)-(b_vect(1,count1)*v_vect_t(1,count1)+b_vect(2,count1)*v_vect_t(2,count1)&
		     &+b_vect(3,count1)*v_vect_t(3,count1))
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
end subroutine em_muller_trans_FMM


subroutine 	get_rhs_em_muller_trans(P0,Pt,vf,ns,srcvals,zpars,RHS)
implicit none

	!List of calling arguments
	integer ( kind = 4 ), intent(in) :: ns
	real ( kind = 8 ), intent(in) :: P0(3),Pt(3)
	complex ( kind = 8 ), intent(in) :: vf(3)
	real ( kind = 8 ), intent(in) :: srcvals(12,ns)
  complex ( kind = 8 ), intent(in) :: zpars(5)
	complex ( kind = 8 ), intent(out) :: RHS(4*ns)
	
	!List of local variables
  complex ( kind = 8 ) omega, ep0,mu0,ep1,mu1
	complex ( kind = 8 ), allocatable :: E(:,:), H(:,:)
	integer count1
	complex ( kind = 8 ) vf_minus(3)
	
	real ( kind = 8 ) ru(3),rv(3),cross_aux(3)
	allocate(E(3,ns), H(3,ns))
		
	omega = zpars(1)
	ep0 = zpars(2)
	mu0 = zpars(3)
	ep1 = zpars(4)
	mu1 = zpars(5)

	vf_minus(1)=-1.0d0*vf(1)
	vf_minus(2)=-1.0d0*vf(2)
	vf_minus(3)=-1.0d0*vf(3)
	write (*,*) 'P0, Pt',P0,Pt
	

	call fieldsEDomega(omega,ep0,mu0,P0,srcvals,ns,E,H,vf,0)
	call fieldsMDomega(omega,ep0,mu0,P0,srcvals,ns,E,H,vf,1)

	call fieldsEDomega(omega,ep1,mu1,Pt,srcvals,ns,E,H,vf_minus,1)
	call fieldsMDomega(omega,ep1,mu1,Pt,srcvals,ns,E,H,vf_minus,1)
!	read (*,*)
	
	do count1=1,ns	
		call orthonormalize(srcvals(4:6,count1),srcvals(10:12,count1),ru,rv)		
		RHS(count1)=-DOT_PRODUCT(rv,E(:,count1))
		RHS(ns+count1)=DOT_PRODUCT(ru,E(:,count1))	
		RHS(2*ns+count1)=-DOT_PRODUCT(rv,H(:,count1))
		RHS(3*ns+count1)=DOT_PRODUCT(ru,H(:,count1))
	enddo

	do count1=1,ns
      RHS(count1)=RHS(count1)/(mu0+mu1)
      RHS(count1+ns)=RHS(count1+ns)/(mu0+mu1)
      RHS(count1+2*ns)=RHS(count1+2*ns)/(ep0+ep1)
      RHS(count1+3*ns)=RHS(count1+3*ns)/(ep0+ep1)
	enddo


return
end subroutine get_rhs_em_muller_trans


      subroutine em_muller_trans_eval_oneside(eps_FMM, sol, zpars, ns, & 
        wts, srcvals, ntarg, targs, iside, E, H)
!
!  Postprocess muller solution at a given collection of targets
!
!  documentation tbd
!
!
      implicit real *8 (a-h,o-z)
      real *8 eps_FMM
      complex *16 sol(4*npts), zpars(5)
      integer ns, ntarg
      real *8 wts(ns), srcvals(12,npts), targs(3,ntarg)
      integer iside
      complex *16 E(3,ntarg), H(3,ntarg)
      complex *16 omega, ep, mu

      omega = zpars(1)

      if(iside.eq.0) then
        ep = zpars(2)
        mu = zpars(3)
      else
        ep = zpars(4)
        mu = zpars(5)
      endif

      call em_muller_trans_FMM_targ(eps_FMM, omega, ep, mu, ns, &
      srcvals, ntarg, targs, wts, sol(1), sol(ns+1), sol(2*ns+1), &
      sol(3*ns+1:4*ns), E, H)

      return
      end



subroutine test_accuracy_em_muller_trans(eps_FMM,sol,zpars,ns,wts,srcvals,P0,vf,Pt)
implicit none

    !List of calling arguments
    complex ( kind = 8 ), intent(in) :: zpars(5)
    integer, intent(in) :: ns
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),eps_FMM
    real ( kind = 8 ), intent(in) :: wts(ns),P0(3),Pt(3)
	  complex ( kind = 8 ), intent(in) :: sol(4*ns),vf(3)
	
    !List of local variables
    complex ( kind = 8 ) omega,ep0,mu0,ep1,mu1
	  complex ( kind = 8 ) ima, R1, R2,Et1(3),Et2(3),Ht1(3),Ht2(3)
	  complex ( kind = 8 ) E01(3),E02(3),H01(3),H02(3)
	  real ( kind = 8 ) ru_s(3),rv_s(3),n_s(3),sour(3),r,dr(3)
	  real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3),error_E,rel_err_E,error_H,rel_err_H
	  real ( kind = 8 ) pi
	  integer count1, iside


		ima=(0.0d0,1.0d0)
		pi=3.1415926535897932384626433832795028841971d0

    omega = zpars(1)
	  ep0 = zpars(2)
 	  mu0 = zpars(3)
	  ep1 = zpars(4)
	  mu1 = zpars(5)

      iside  = 0
      call em_muller_trans_eval_oneside(eps_FMM, sol, zpars, ns, & 
        wts, srcvals, 1, Pt, iside, Et1, Ht1)
      iside = 1
      call em_muller_trans_eval_oneside(eps_FMM, sol, zpars, ns, & 
        wts, srcvals, 1, P0, iside, E01, H01)

!      call em_muller_trans_FMM_targ(eps_FMM,omega,ep0,mu0,ns,srcvals,1,Pt,wts,&
!     &sol(1:ns),sol(ns+1:2*ns),sol(2*ns+1:3*ns),sol(3*ns+1:4*ns),Et1,Ht1)
!		call em_muller_trans_FMM_targ(eps_FMM,omega,ep1,mu1,ns,srcvals,1,P0,wts,&
!     &sol(1:ns),sol(ns+1:2*ns),sol(2*ns+1:3*ns),sol(3*ns+1:4*ns),E01,H01)

        call fieldsEDomega(omega,ep0,mu0,P0,Pt,1,Et2,Ht2,vf,0)
        call fieldsMDomega(omega,ep0,mu0,P0,Pt,1,Et2,Ht2,vf,1)

        call fieldsEDomega(omega,ep1,mu1,Pt,P0,1,E02,H02,vf,0)
        call fieldsMDomega(omega,ep1,mu1,Pt,P0,1,E02,H02,vf,1)

        write (*,*) 'Errors at the EXTERIOR region:'
        error_E=sqrt(abs(Et1(1)-Et2(1))**2+abs(Et1(2)-Et2(2))**2+abs(Et1(3)-Et2(3))**2)
        call prin2('Et1=*',Et1,6)
        call prin2('Et2=*',Et2,6)
        write (*,*) 'Absolute Error in E: ', error_E
        write (*,*) 'Relative Error in E: ', error_E/sqrt(abs(Et2(1))**2+abs(Et2(2))**2+abs(Et2(3))**2)

        error_H=sqrt(abs(Ht1(1)-Ht2(1))**2+abs(Ht1(2)-Ht2(2))**2+abs(Ht1(3)-Ht2(3))**2)
        write (*,*) 'Absolute Error in H: ', error_H
        write (*,*) 'Relative Error in H: ', error_H/sqrt(abs(Ht2(1))**2+abs(Ht2(2))**2+abs(Ht2(3))**2)

        write (*,*) 'Errors at the INTERIOR region:'
        call prin2('Et1=*',E01,6)
        call prin2('Et2=*',E02,6)
        error_E=sqrt(abs(E01(1)-E02(1))**2+abs(E01(2)-E02(2))**2+abs(E01(3)-E02(3))**2)
        write (*,*) 'Absolute Error in E: ', error_E
        write (*,*) 'Relative Error in E: ', error_E/sqrt(abs(E02(1))**2+abs(E02(2))**2+abs(E02(3))**2)

        error_H=sqrt(abs(H01(1)-H02(1))**2+abs(H01(2)-H02(2))**2+abs(H01(3)-H02(3))**2)
        write (*,*) 'Absolute Error in H: ', error_H
        write (*,*) 'Relative Error in H: ', error_H/sqrt(abs(H02(1))**2+abs(H02(2))**2+abs(H02(3))**2)

return
end subroutine test_accuracy_em_muller_trans


subroutine em_muller_trans_FMM_targ(eps,omega,ep,mu,ns,srcvals,nt,targ,wts,a_u,a_v,b_u,b_v,E,H)
implicit none

    !List of calling arguments
	real ( kind = 8 ), intent(in) :: eps
	complex ( kind = 8 ), intent(in) :: omega,ep,mu
    integer, intent(in) :: ns,nt
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),targ(3,nt)
    real ( kind = 8 ), intent(in) :: wts(ns)
    complex ( kind = 8 ), intent(in) :: a_u(ns),a_v(ns),b_u(ns),b_v(ns)
    complex ( kind = 8 ), intent(out) :: E(3,nt),H(3,nt)

    !List of local variables
	real ( kind = 8 ), allocatable :: n_vect(:,:),u_vect(:,:),v_vect(:,:),source(:,:)
    complex ( kind = 8 ), allocatable :: a_vect(:,:),b_vect(:,:),lambda(:),rho(:)
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
	call orthonormalize_all(srcvals(4:6,:),srcvals(10:12,:),u_vect,v_vect,ns)

    do count1=1,ns
		a_vect(1,count1)=(b_u(count1)*u_vect(1,count1)+b_v(count1)*v_vect(1,count1))/(-ima*omega)
		a_vect(2,count1)=(b_u(count1)*u_vect(2,count1)+b_v(count1)*v_vect(2,count1))/(-ima*omega)
		a_vect(3,count1)=(b_u(count1)*u_vect(3,count1)+b_v(count1)*v_vect(3,count1))/(-ima*omega)
				
        b_vect(1,count1)=(a_u(count1)*u_vect(1,count1)+a_v(count1)*v_vect(1,count1))*mu
        b_vect(2,count1)=(a_u(count1)*u_vect(2,count1)+a_v(count1)*v_vect(2,count1))*mu
        b_vect(3,count1)=(a_u(count1)*u_vect(3,count1)+a_v(count1)*v_vect(3,count1))*mu
	enddo

    !Computing the full operator
    ifa_vect=1
    ifb_vect=1
    iflambda=0
    ifrho=0
    ifE=0
    ifcurlE=1
    ifdivE=0

call Vector_Helmholtz_targ(eps,zk,ns,source,wts,ifa_vect,a_vect,ifb_vect,&
 &b_vect,iflambda,lambda,ifrho,rho,n_vect,ifE,curlE,ifcurlE,E,ifdivE,divE,nt,targ)	!watch out!! the E and curlE are fliped

    do count1=1,ns
		a_vect(1,count1)=(a_u(count1)*u_vect(1,count1)+a_v(count1)*v_vect(1,count1))/(ima*omega)
		a_vect(2,count1)=(a_u(count1)*u_vect(2,count1)+a_v(count1)*v_vect(2,count1))/(ima*omega)
		a_vect(3,count1)=(a_u(count1)*u_vect(3,count1)+a_v(count1)*v_vect(3,count1))/(ima*omega)
				
        b_vect(1,count1)=(b_u(count1)*u_vect(1,count1)+b_v(count1)*v_vect(1,count1))*ep
        b_vect(2,count1)=(b_u(count1)*u_vect(2,count1)+b_v(count1)*v_vect(2,count1))*ep
        b_vect(3,count1)=(b_u(count1)*u_vect(3,count1)+b_v(count1)*v_vect(3,count1))*ep
	enddo

    !Computing the full operator
    ifa_vect=1
    ifb_vect=1
    iflambda=0
    ifrho=0
    ifE=0
    ifcurlE=1
    ifdivE=0

call Vector_Helmholtz_targ(eps,zk,ns,source,wts,ifa_vect,a_vect,ifb_vect,&
 &b_vect,iflambda,lambda,ifrho,rho,n_vect,ifE,curlE,ifcurlE,H,ifdivE,divE,nt,targ)	!watch out!! the E and curlE are fliped


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
end subroutine em_muller_trans_FMM_targ





subroutine em_muller_trans_FMM2(eps,zpars,ns,nt,srcvals,targvals,wts,&
 &a_u,a_v,b_u,b_v,AA_u,AA_v,BB_u,BB_v,thresh,ifdir)
implicit none

    !List of calling arguments
    real ( kind = 8 ), intent(in) :: eps
    complex ( kind = 8 ), intent(in) :: zpars(5)
    integer, intent(in) :: ns,nt
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),wts(ns)
    real ( kind = 8 ), intent(in) :: targvals(12,nt)
    complex ( kind = 8 ), intent(in) :: a_u(ns),a_v(ns)
    complex ( kind = 8 ), intent(in) :: b_u(ns),b_v(ns)
    complex ( kind = 8 ), intent(out) :: AA_u(nt),AA_v(nt)
    complex ( kind = 8 ), intent(out) :: BB_u(nt),BB_v(nt)
    real ( kind = 8 ), intent(in) :: thresh
    integer, intent(in) :: ifdir 


    !List of local variables
	  real ( kind = 8 ), allocatable :: u_vect_s(:,:),v_vect_s(:,:)
	  real ( kind = 8 ), allocatable :: source(:,:),n_vect_s(:,:)

	  real ( kind = 8 ), allocatable :: n_vect_t(:,:),u_vect_t(:,:)
	  real ( kind = 8 ), allocatable :: targets(:,:),v_vect_t(:,:)

    complex ( kind = 8 ), allocatable :: a_vect(:,:,:),b_vect(:,:,:)
    complex ( kind = 8 ), allocatable :: b_vect_t(:,:,:)
    complex ( kind = 8 ), allocatable :: lambda(:,:),rho(:,:)
    complex ( kind = 8 ), allocatable :: E(:,:,:),curlE(:,:,:),divE(:,:)
    complex ( kind = 8 ) ima,zk0,zk1

    complex ( kind = 8 ) omega,ep0,mu0,ep1,mu1

    integer count1,count2,nd
    integer ifa_vect,ifb_vect,iflambda,ifrho,ifE,ifcurlE,ifdivE

    ima=(0.0d0,1.0d0)

    omega=zpars(1)
    ep0=zpars(2)
    mu0=zpars(3)
    ep1=zpars(4)
    mu1=zpars(5)

    zk0=omega*sqrt(ep0*mu0)
    zk1=omega*sqrt(ep1*mu1)

    nd=2
    allocate(a_vect(nd,3,ns))
    allocate(b_vect(nd,3,ns))
    allocate(b_vect_t(nd,3,nt))
    allocate(lambda(nd,ns))
    allocate(rho(nd,ns))
    allocate(E(nd,3,nt))
    allocate(curlE(nd,3,nt))
    allocate(divE(nd,nt))
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
      a_vect(1,1,count1)=(b_u(count1)*u_vect_s(1,count1) + &
         b_v(count1)*v_vect_s(1,count1))/(-ima*omega)
      a_vect(1,2,count1)=(b_u(count1)*u_vect_s(2,count1) + &
         b_v(count1)*v_vect_s(2,count1))/(-ima*omega)
      a_vect(1,3,count1)=(b_u(count1)*u_vect_s(3,count1) + &
         b_v(count1)*v_vect_s(3,count1))/(-ima*omega)
      b_vect(1,1,count1)=(a_u(count1)*u_vect_s(1,count1) + &
        a_v(count1)*v_vect_s(1,count1))*mu0
      b_vect(1,2,count1)=(a_u(count1)*u_vect_s(2,count1) + &
        a_v(count1)*v_vect_s(2,count1))*mu0
      b_vect(1,3,count1)=(a_u(count1)*u_vect_s(3,count1) + &
        a_v(count1)*v_vect_s(3,count1))*mu0


      b_vect(2,1,count1)=(b_u(count1)*u_vect_s(1,count1) + &
        b_v(count1)*v_vect_s(1,count1))*ep0
      b_vect(2,2,count1)=(b_u(count1)*u_vect_s(2,count1) + &
        b_v(count1)*v_vect_s(2,count1))*ep0
      b_vect(2,3,count1)=(b_u(count1)*u_vect_s(3,count1) + &
        b_v(count1)*v_vect_s(3,count1))*ep0

      a_vect(2,1,count1)=(a_u(count1)*u_vect_s(1,count1) + &
        a_v(count1)*v_vect_s(1,count1))/(ima*omega)
      a_vect(2,2,count1)=(a_u(count1)*u_vect_s(2,count1) + &
        a_v(count1)*v_vect_s(2,count1))/(ima*omega)
      a_vect(2,3,count1)=(a_u(count1)*u_vect_s(3,count1) + &
        a_v(count1)*v_vect_s(3,count1))/(ima*omega)

     enddo

    !Computing the full operator
    ifa_vect=1
    ifb_vect=1
    iflambda=0
    ifrho=0
    ifE=0
    ifcurlE=1
    ifdivE=0

    call Vector_Helmholtz_targ2_vect(nd,eps,zk0,ns,source,wts, &
      ifa_vect,a_vect,&
      ifb_vect,b_vect,iflambda,lambda,ifrho,rho,n_vect_s,ifE,E,ifcurlE,&
      curlE,ifdivE,divE,nt,targets,thresh,ifdir)

    do count1=1,nt
        b_vect_t(1,1,count1)=n_vect_t(2,count1)*curlE(1,3,count1)-&
          n_vect_t(3,count1)*curlE(1,2,count1)
        b_vect_t(1,2,count1)=n_vect_t(3,count1)*curlE(1,1,count1)-&
          n_vect_t(1,count1)*curlE(1,3,count1)
        b_vect_t(1,3,count1)=n_vect_t(1,count1)*curlE(1,2,count1)-&
          n_vect_t(2,count1)*curlE(1,1,count1)

        b_vect_t(2,1,count1)=n_vect_t(2,count1)*curlE(2,3,count1)-&
          n_vect_t(3,count1)*curlE(2,2,count1)
        b_vect_t(2,2,count1)=n_vect_t(3,count1)*curlE(2,1,count1)-&
          n_vect_t(1,count1)*curlE(2,3,count1)
        b_vect_t(2,3,count1)=n_vect_t(1,count1)*curlE(2,2,count1)-&
          n_vect_t(2,count1)*curlE(2,1,count1)
    enddo

    do count1=1,nt
      AA_u(count1)=b_vect_t(1,1,count1)*u_vect_t(1,count1) + &
          b_vect_t(1,2,count1)*u_vect_t(2,count1) + &
          b_vect_t(1,3,count1)*u_vect_t(3,count1)
      AA_v(count1)=b_vect_t(1,1,count1)*v_vect_t(1,count1) + &
          b_vect_t(1,2,count1)*v_vect_t(2,count1) +  &
          b_vect_t(1,3,count1)*v_vect_t(3,count1)

      BB_u(count1)=b_vect_t(2,1,count1)*u_vect_t(1,count1) + &
          b_vect_t(2,2,count1)*u_vect_t(2,count1) + &
          b_vect_t(2,3,count1)*u_vect_t(3,count1)
      BB_v(count1)=b_vect_t(2,1,count1)*v_vect_t(1,count1) + &
          b_vect_t(2,2,count1)*v_vect_t(2,count1) + &
          b_vect_t(2,3,count1)*v_vect_t(3,count1)
    enddo


!!! now zk1


    do count1=1,ns
      a_vect(1,1,count1)=(b_u(count1)*u_vect_s(1,count1) + &
         b_v(count1)*v_vect_s(1,count1))/(-ima*omega)
      a_vect(1,2,count1)=(b_u(count1)*u_vect_s(2,count1) + &
         b_v(count1)*v_vect_s(2,count1))/(-ima*omega)
      a_vect(1,3,count1)=(b_u(count1)*u_vect_s(3,count1) + &
            b_v(count1)*v_vect_s(3,count1))/(-ima*omega)

      b_vect(1,1,count1)=(a_u(count1)*u_vect_s(1,count1) + &
        a_v(count1)*v_vect_s(1,count1))*mu1
      b_vect(1,2,count1)=(a_u(count1)*u_vect_s(2,count1) + &
        a_v(count1)*v_vect_s(2,count1))*mu1
      b_vect(1,3,count1)=(a_u(count1)*u_vect_s(3,count1) + &
        a_v(count1)*v_vect_s(3,count1))*mu1

      b_vect(2,1,count1)=(b_u(count1)*u_vect_s(1,count1) + &
        b_v(count1)*v_vect_s(1,count1))*ep1
      b_vect(2,2,count1)=(b_u(count1)*u_vect_s(2,count1) + &
        b_v(count1)*v_vect_s(2,count1))*ep1
      b_vect(2,3,count1)=(b_u(count1)*u_vect_s(3,count1) + &
        b_v(count1)*v_vect_s(3,count1))*ep1

      a_vect(2,1,count1)=(a_u(count1)*u_vect_s(1,count1) + &
        a_v(count1)*v_vect_s(1,count1))/(ima*omega)
      a_vect(2,2,count1)=(a_u(count1)*u_vect_s(2,count1) + &
        a_v(count1)*v_vect_s(2,count1))/(ima*omega)
      a_vect(2,3,count1)=(a_u(count1)*u_vect_s(3,count1) + &
        a_v(count1)*v_vect_s(3,count1))/(ima*omega)
    enddo

    ifa_vect=1
    ifb_vect=1
    iflambda=0
    ifrho=0
    ifE=0
    ifcurlE=1
    ifdivE=0

    call Vector_Helmholtz_targ2_vect(nd,eps,zk1,ns,source,wts, &
      ifa_vect,a_vect,ifb_vect,b_vect,iflambda,lambda,ifrho,rho, &
      n_vect_s,ifE,E,ifcurlE,curlE,ifdivE,divE,nt,targets,thresh, &
      ifdir)


    do count1=1,nt
      b_vect_t(1,1,count1)=n_vect_t(2,count1)*curlE(1,3,count1)- &
        n_vect_t(3,count1)*curlE(1,2,count1)
      b_vect_t(1,2,count1)=n_vect_t(3,count1)*curlE(1,1,count1)-&
        n_vect_t(1,count1)*curlE(1,3,count1)
      b_vect_t(1,3,count1)=n_vect_t(1,count1)*curlE(1,2,count1)-&
        n_vect_t(2,count1)*curlE(1,1,count1)

      b_vect_t(2,1,count1)=n_vect_t(2,count1)*curlE(2,3,count1)-&
        n_vect_t(3,count1)*curlE(2,2,count1)
      b_vect_t(2,2,count1)=n_vect_t(3,count1)*curlE(2,1,count1)-&
        n_vect_t(1,count1)*curlE(2,3,count1)
      b_vect_t(2,3,count1)=n_vect_t(1,count1)*curlE(2,2,count1)-&
        n_vect_t(2,count1)*curlE(2,1,count1)
    enddo

    do count1=1,nt
        AA_u(count1)=AA_u(count1) - &
          (b_vect_t(1,1,count1)*u_vect_t(1,count1) + &
           b_vect_t(1,2,count1)*u_vect_t(2,count1) + &
           b_vect_t(1,3,count1)*u_vect_t(3,count1))
        AA_v(count1)=AA_v(count1) - &
          (b_vect_t(1,1,count1)*v_vect_t(1,count1) + &
           b_vect_t(1,2,count1)*v_vect_t(2,count1) + &
           b_vect_t(1,3,count1)*v_vect_t(3,count1))

        BB_u(count1)=BB_u(count1) - &
          (b_vect_t(2,1,count1)*u_vect_t(1,count1) + &
           b_vect_t(2,2,count1)*u_vect_t(2,count1) + &
           b_vect_t(2,3,count1)*u_vect_t(3,count1))

        BB_v(count1)=BB_v(count1) - &
         (b_vect_t(2,1,count1)*v_vect_t(1,count1) + &
          b_vect_t(2,2,count1)*v_vect_t(2,count1) + &
          b_vect_t(2,3,count1)*v_vect_t(3,count1))
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
end subroutine em_muller_trans_FMM2
