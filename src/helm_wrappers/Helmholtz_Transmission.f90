      subroutine h_transmission_solver(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,eps,zpars,numit,ifinout,&
     &rhs,eps_gmres,niter,errs,rres,soln)
!
!
!      This subroutine solves the helmholtz transmission problem:
!
!       \Delta u1 + zk1^2 u1 = 0  (interior region)
!       \Delta u0 + zk0^2 u0 = 0  (exterior region)
!
!      with junp conditions:
!	     u0-u1 = f
!		 (1/ep0)u0'-(1/ep1)u1' = g
!
!      Representation:
!        u1 = ep1 D_{k}\rho + ep1^2 S_{k}\lambda
!        u0 = ep0 D_{k}\rho + ep0^2 S_{k}\lambda
!
!	   The parameters verify:
! 		 zk0=omega*sqrt(ep0*mu0)
! 		 zk1=omega*sqrt(ep1*mu1)
!
!     The linear system is solved iteratively using GMRES
!     until a relative residual of 1e-15 is reached
!
!
!       input:
!         npatches - integer
!            number of patches
!
!         norders- integer(npatches)
!            order of discretization on each patch 
!
!         ixyzs - integer(npatches+1)
!            ixyzs(i) denotes the starting location in srccoefs,
!               and srcvals array corresponding to patch i
!   
!         iptype - integer(npatches)
!            type of patch
!             iptype = 1, triangular patch discretized using RV nodes
!
!         npts - integer
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
!          zpars - complex *16 (5)
!              kernel parameters
!              zpars(1) = omega 
!              zpars(2) = ep0
!              zpars(3) = mu0
!              zpars(4) = ep
!              zpars(5) = mu
!
!          ifinout - integer
!              flag for interior or exterior problems (normals assumed to 
!                be pointing in exterior of region)
!              ifinout = 0, interior problem
!              ifinout = 1, exterior problem
!
!           rhs - complex *16(2*npts)
!              right hand side
!				rhs(1:npts)=f
!				rhs(npts+1:2*npts)=g
!
!           eps_gmres - real *8
!                gmres tolerance requested
!
!           numit - integer
!              max number of gmres iterations
!
!         output:
!           niter - integer
!              number of gmres iterations required for relative residual
!               to converge to 1e-15
!          
!           errs(1:iter) - relative residual as a function of iteration
!              number
! 
!           rres - real *8
!              relative residual for computed solution
!              
!           soln - complex *16(2*npts)
!              densities \rho,\lambda which solve the dirichlet problem:
!			       soln(1:npts)=\rho
!			       soln(npts+1:npts)=\lambda
!
!
      implicit none
      integer npatches,norder,npols,npts
      integer ifinout
      integer norders(npatches),ixyzs(npatches+1)
      integer iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,eps_gmres
      complex *16 zpars(5)
      complex *16 rhs(2*npts)
      complex *16 soln(2*npts)

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

      complex *16, allocatable :: wnear11(:),wnear12(:)
      complex *16, allocatable :: wnear21(:),wnear22(:)

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
!   as we have two scalar unknowns: \rho,\lambda, we need n_var=2*npts
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
!    this might need fixing
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

!move this inside get_far_irder.. or better, ask Manas..

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
      allocate(wnear11(nquad),wnear12(nquad))
      allocate(wnear21(nquad),wnear22(nquad))
      
!C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,nquad
	    wnear11(i)=0
	    wnear12(i)=0
	    wnear21(i)=0
	    wnear22(i)=0
      enddo
!C$OMP END PARALLEL DO    


      iquadtype = 1

!!      eps2 = 1.0d-8

      print *, "starting to generate near quadrature"
      call cpu_time(t1)
!C$      t1 = omp_get_wtime()      
	 
      call getnearquad_h_transmission(npatches,norders,&
     &ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,&
     &ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,&
     &iquad,rfac0,nquad,wnear11,wnear12,wnear21,wnear22)
	 
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
!      zid = (ep0+ep1)/2.0d0
	  
	  zid=(zpars(2)+zpars(4))/2.0d0

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
!        evaluation routine  (in that case helmholtz transmission)
!


        call lpcomp_h_transmission_addsub(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,&
     &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,&
     &vmat(1,it),novers,npts_over,ixyzso,srcover,wover,wtmp,&
	 &wnear11,wnear12,wnear21,wnear22)

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


          call lpcomp_h_transmission_addsub(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,&
     &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,&
     &soln,novers,npts_over,ixyzso,srcover,wover,wtmp,&
	 &wnear11,wnear12,wnear21,wnear22)

            
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
      end subroutine h_transmission_solver




      subroutine getnearquad_h_transmission(npatches,norders,&
     &ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,&
     &iquad,rfac0,nquad,wnear11,wnear12,wnear21,wnear22)
!
!       this subroutine generates the near field quadrature
!       for the representation:
!        u1 = ep1 D_{k}\rho + ep1^2 S_{k}\lambda
!        u0 = ep0 D_{k}\rho + ep0^2 S_{k}\lambda
!
!
!        Note: the 4 \pi scaling is NOT!! included as the output of the FMM
!        has been rescaled.
!
!
!       The quadrature is computed by the following strategy
!        targets within a sphere of radius rfac0*rs
!        of a chunk centroid is handled using adaptive integration
!        where rs is the radius of the bounding sphere
!        for the patch
!  
!       All other targets in the near field are handled via
!        oversampled quadrature
!
!       The recommended parameter for rfac0 is 1.25d0
!
!        
!
!       input:
!         npatches - integer
!            number of patches
!
!         norders - integer(npatches)
!            order of discretization on each patch 
!
!         ixyzs - integer(npatches+1)
!            starting location of data on patch i
!  
!         iptype - integer(npatches)
!           type of patch
!           iptype = 1 -> triangular patch discretized with RV nodes
!
!         npts - integer
!            total number of discretization points on the boundary
!
!         srccoefs - real *8 (9,npts)
!            koornwinder expansion coefficients of xyz, dxyz/du,
!            and dxyz/dv on each patch. 
!            For each point srccoefs(1:3,i) is xyz info
!                           srccoefs(4:6,i) is dxyz/du info
!                           srccoefs(7:9,i) is dxyz/dv info
!
!          srcvals - real *8 (12,npts)
!             xyz(u,v) and derivative info sampled at the 
!             discretization nodes on the surface
!             srcvals(1:3,i) - xyz info
!             srcvals(4:6,i) - dxyz/du info
!             srcvals(7:9,i) - dxyz/dv info
! 
!         ndtarg - integer
!            leading dimension of target array
!        
!         ntarg - integer
!            number of targets
!
!         targs - real *8 (ndtarg,ntarg)
!            target information
!
!         ipatch_id - integer(ntarg)
!            id of patch of target i, id = -1, if target is off-surface
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
!          zpars - complex *16 (5)
!              kernel parameters (Referring to formula (1))
!              zpars(1) = omega 
!              zpars(2) = ep0
!              zpars(3) = mu0
!              zpars(4) = ep
!              zpars(5) = mu
!
!           iquadtype - integer
!              quadrature type
!              iquadtype = 1, use ggq for self + adaptive integration
!                 for rest
! 
!
!           nnz - integer
!             number of source patch-> target interactions in the near field
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
!               location in wnear_ij array where quadrature for col_ind(i)
!               starts
!
!           rfac0 - integer
!               radius parameter for near field
!
!           nquad - integer
!               number of entries in wnear_ij
!
!        output
!            wnear11,wnear12,wnear21,wnear22 - complex *16(nquad)
!               the desired near field quadrature
!               
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
      complex *16 wnear11(nquad),wnear12(nquad)
	  complex *16 wnear21(nquad),wnear22(nquad)

      integer ipatch_id(ntarg)
      real *8 uvs_targ(2,ntarg)

      complex *16 alpha,beta
      integer i,j,ndi,ndd,ndz

      integer ipv

      procedure (), pointer :: fker
	  external  fker_h_transmissions

     ndz=5
	 ndd=1
	 ndi=2
	 ipv=1

      fker =>  fker_h_transmissions
	  ipars(1)=1
	  ipars(2)=1
	  call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
     &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear11)
      ipars(1)=1
	  ipars(2)=2
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,&
     &eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,col_ind,iquad,&
     &rfac0,nquad,wnear12)
      ipars(1)=2
	  ipars(2)=1
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,&
     &eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,col_ind,iquad,&
     &rfac0,nquad,wnear21)
      ipars(1)=2
	  ipars(2)=2
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,&
     &eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,col_ind,iquad,&
     &rfac0,nquad,wnear22)

      return
      end subroutine getnearquad_h_transmission


subroutine fker_h_transmissions(srcinfo, ndt,targinfo,ndd, dpars,ndz,zpars,ndi,ipars,E_val)
implicit none

! this function provides the near field kernel that will use zgetnearquad_ggq_guru inside
! getnearquad_h_neumann

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
	complex ( kind = 8 ) Dk1rho(1,1),Dk0rho(1,1),Sk1lambda(1,1),Sk0lambda(1,1)
	complex ( kind = 8 ) ngradSk1lambda(1,1),ngradSk0lambda(1,1)
	complex ( kind = 8 ) ngradDk1rho(1,1),ngradDk0rho(1,1),E_mat(2,2)
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3),xprod_aux3(3),xprod_aux4(3)
	complex ( kind = 8 ) R1_0,R1_1,R2_0,R2_1,ima,my_exp_0,my_exp_1, omega,ep0,mu0,ep1,mu1,zk0,zk1
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

	if (ipars(1).eq.1) then
	  if (ipars(2).eq.1) then
	  	R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
    	R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
	    my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
	    my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
	  	call get_Dkrho(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_0,my_exp_0,r,Dk0rho)
	    call get_Dkrho(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_1,my_exp_1,r,Dk1rho)
		E_val=ep0*Dk0rho(1,1)-ep1*Dk1rho(1,1)
      else	
	    my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
	    my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
	    call get_Sklambda(my_exp_0,r,Sk0lambda)
	    call get_Sklambda(my_exp_1,r,Sk1lambda)
		E_val=ep0**2*Sk0lambda(1,1)-ep1**2*Sk1lambda(1,1)
	  endif
	else	
	  if (ipars(2).eq.1) then
		R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
    	R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
		R2_0=((ima*zk0)**2/r**3-3.0d0*ima*zk0/r**4+3.0d0/r**5)*exp(ima*zk0*r)/(4.0d0*pi)
		R2_1=((ima*zk1)**2/r**3-3.0d0*ima*zk1/r**4+3.0d0/r**5)*exp(ima*zk1*r)/(4.0d0*pi)
	    my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
	    my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
	    call get_ngradDkrho(n_s,n_t,dr,R1_0,R2_0,zk0,my_exp_0,r,ngradDk0rho)
	    call get_ngradDkrho(n_s,n_t,dr,R1_1,R2_1,zk1,my_exp_1,r,ngradDk1rho)
		E_val=ngradDk0rho(1,1)-ngradDk1rho(1,1)
      else		  
	  	R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
    	R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
	    my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
	    my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
	  	call get_ngradSklambda(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_0,my_exp_0,r,ngradSk0lambda)
    	call get_ngradSklambda(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_1,my_exp_1,r,ngradSk1lambda)
		E_val=ep0*ngradSk0lambda(1,1)-ep1*ngradSk1lambda(1,1)
	  endif	
	endif

return
end subroutine fker_h_transmissions

subroutine fker_h_transmissions_not_used(srcinfo, ndt,targinfo,ndd, dpars,ndz,zpars,ndi,ipars,E_val)
implicit none

! this function provides the near field kernel that will use zgetnearquad_ggq_guru inside
! getnearquad_h_transmission

! This subroutine is not used. It does the same as fker_h_neumanns but is prepared for adaptive integration
! of different kernels in parallel (in case it end up being the most effective way)


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
	complex ( kind = 8 ) Dk1rho(1,1),Dk0rho(1,1),Sk1lambda(1,1),Sk0lambda(1,1)
	complex ( kind = 8 ) ngradSk1lambda(1,1),ngradSk0lambda(1,1)
	complex ( kind = 8 ) ngradDk1rho(1,1),ngradDk0rho(1,1),E_mat(2,2)
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3),xprod_aux3(3),xprod_aux4(3)
	complex ( kind = 8 ) R1_0,R1_1,R2_0,R2_1,ima,my_exp_0,my_exp_1, omega,ep0,mu0,ep1,mu1,zk0,zk1
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

		R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
		R2_0=((ima*zk0)**2/r**3-3.0d0*ima*zk0/r**4+3.0d0/r**5)*exp(ima*zk0*r)/(4.0d0*pi)
		my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
		R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
		R2_1=((ima*zk1)**2/r**3-3.0d0*ima*zk1/r**4+3.0d0/r**5)*exp(ima*zk1*r)/(4.0d0*pi)
		my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
				
		call orthonormalize(srcinfo(4:6),n_s,ru_s,rv_s)
		call orthonormalize(targinfo(4:6),n_t,ru_t,rv_t)

		call get_Sklambda(my_exp_0,r,Sk0lambda)
		call get_Sklambda(my_exp_1,r,Sk1lambda)
		call get_Dkrho(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_0,my_exp_0,r,Dk0rho)
		call get_Dkrho(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_1,my_exp_1,r,Dk1rho)
		call get_ngradSklambda(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_0,my_exp_0,r,ngradSk0lambda)
		call get_ngradSklambda(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_1,my_exp_1,r,ngradSk1lambda)

		call get_ngradDkrho(n_s,n_t,dr,R1_0,R2_0,zk0,my_exp_0,r,ngradDk0rho)
		call get_ngradDkrho(n_s,n_t,dr,R1_1,R2_1,zk1,my_exp_1,r,ngradDk1rho)

		E_mat(1,1)=ep0*Dk0rho(1,1)-ep1*Dk1rho(1,1)
		E_mat(1,2)=ep0**2*Sk0lambda(1,1)-ep1**2*Sk1lambda(1,1)
		E_mat(2,1)=ngradDk0rho(1,1)-ngradDk1rho(1,1)
		E_mat(2,2)=ep0*ngradSk0lambda(1,1)-ep1*ngradSk1lambda(1,1)
		E_val=E_mat(ipars(1),ipars(2))


return
end subroutine fker_h_transmissions_not_used


subroutine fker_h_transmission(ns,nt,srcinfo, targinfo, dpars,zpars,ipars,sigma,rho,wts,E,thresh)
implicit none

! this function does something very similar to fker_h_transmissions but oriented to 
! remove the near interaction of the FMM 


    !List of calling arguments
	integer, intent(in) :: ns, nt
	real ( kind = 8 ), intent(in) :: srcinfo(12,ns)
	real ( kind = 8 ), intent(in) :: targinfo(12,nt)
	integer, intent(in) :: ipars(2)
	real ( kind = 8 ), intent(in) :: dpars(1)
	complex ( kind = 8 ), intent(in) :: zpars(5),sigma(ns),rho(ns)
	complex ( kind = 8 ), intent(out) :: E(2)
    real ( kind = 8 ), intent(in) :: thresh,wts(ns)
	
	!List of local variables
	real ( kind = 8 ) ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ) ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ) du(3), dv(3),sour(3),targ(3)
	real ( kind = 8 ) r, dr(3),aux_real
	complex ( kind = 8 ) Dk1rho(1,1),Dk0rho(1,1),Sk1lambda(1,1),Sk0lambda(1,1)
	complex ( kind = 8 ) ngradSk1lambda(1,1),ngradSk0lambda(1,1)
	complex ( kind = 8 ) ngradDk1rho(1,1),ngradDk0rho(1,1),E_mat(2,2)
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3),xprod_aux3(3),xprod_aux4(3)
	complex ( kind = 8 ) R1_0,R1_1,R2_0,R2_1,ima,my_exp_0,my_exp_1, omega,ep0,mu0,ep1,mu1,zk0,zk1
	real ( kind = 8 ) pi
	integer count1,count2
	
	pi=3.1415926535897932384626433832795028841971d0
	ima=(0.0d0,1.0d0)
	omega=zpars(1)
	ep0=zpars(2)
	mu0=zpars(3)
	ep1=zpars(4)
	mu1=zpars(5)
	do count1=1,nt
	  E(1)=0.0d0
	  E(2)=0.0d0
	  
	  do count2=1,ns

		sour(1)=srcinfo(1,count2)
		sour(2)=srcinfo(2,count2)
		sour(3)=srcinfo(3,count2)
		
		n_s(1)=srcinfo(10,count2)
		n_s(2)=srcinfo(11,count2)
		n_s(3)=srcinfo(12,count2)	

		targ(1)=targinfo(1,count1)
		targ(2)=targinfo(2,count1)
		targ(3)=targinfo(3,count1)

		n_t(1)=targinfo(10,count1)
		n_t(2)=targinfo(11,count1)
		n_t(3)=targinfo(12,count1)

		dr(1)=targ(1)-sour(1)
		dr(2)=targ(2)-sour(2)
		dr(3)=targ(3)-sour(3)
		
		r=sqrt((dr(1))**2+(dr(2))**2+(dr(3))**2)
	    if (r.ge.thresh) then
			zk0=omega*sqrt(ep0*mu0)
			zk1=omega*sqrt(ep1*mu1)

			R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
			R2_0=((ima*zk0)**2/r**3-3.0d0*ima*zk0/r**4+3.0d0/r**5)*exp(ima*zk0*r)/(4.0d0*pi)
			my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
			R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
			R2_1=((ima*zk1)**2/r**3-3.0d0*ima*zk1/r**4+3.0d0/r**5)*exp(ima*zk1*r)/(4.0d0*pi)
			my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
					
			call orthonormalize(srcinfo(4:6,count2),n_s,ru_s,rv_s)
			call orthonormalize(targinfo(4:6,count1),n_t,ru_t,rv_t)

			call get_Sklambda(my_exp_0,r,Sk0lambda)
			call get_Sklambda(my_exp_1,r,Sk1lambda)
			call get_Dkrho(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_0,my_exp_0,r,Dk0rho)
			call get_Dkrho(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_1,my_exp_1,r,Dk1rho)
			call get_ngradSklambda(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_0,my_exp_0,r,ngradSk0lambda)
			call get_ngradSklambda(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_1,my_exp_1,r,ngradSk1lambda)

			call get_ngradDkrho(n_s,n_t,dr,R1_0,R2_0,zk0,my_exp_0,r,ngradDk0rho)
			call get_ngradDkrho(n_s,n_t,dr,R1_1,R2_1,zk1,my_exp_1,r,ngradDk1rho)

			E_mat(1,1)=ep0*Dk0rho(1,1)-ep1*Dk1rho(1,1)
			E_mat(1,2)=ep0**2*Sk0lambda(1,1)-ep1**2*Sk1lambda(1,1)
			E_mat(2,1)=ngradDk0rho(1,1)-ngradDk1rho(1,1)
			E_mat(2,2)=ep0*ngradSk0lambda(1,1)-ep1*ngradSk1lambda(1,1)

		else

			E_mat(1,1)=0.0d0
			E_mat(1,2)=0.0d0
			E_mat(2,1)=0.0d0
			E_mat(2,2)=0.0d0

		endif
		
		E(1)=E(1)+E_mat(1,1)*sigma(count2)*wts(count2)+E_mat(1,2)*rho(count2)*wts(count2)
		E(2)=E(2)+E_mat(2,1)*sigma(count2)*wts(count2)+E_mat(2,2)*rho(count2)*wts(count2)		

	  enddo
	  enddo
return
end subroutine fker_h_transmission


      subroutine lpcomp_h_transmission_addsub(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,sigma,novers,&
     &nptso,ixyzso,srcover,whtsover,pot,wnear11,wnear12,&
	 &wnear21,wnear22)

!
!      this subroutine evaluates the layer potential (jump in u and u') for
!      the representation:
!        u1 = ep1 D_{k}\rho + ep1^2 S_{k}\lambda
!        u0 = ep0 D_{k}\rho + ep0^2 S_{k}\lambda
!      where the near field is precomputed and stored
!      in the row sparse compressed format.
!
!       Note the 4\pi scaling is NOT included as the FMM output was scaled
!       appropriately
!
!
!     The fmm is used to accelerate the far-field and 
!     near-field interactions are handled via precomputed quadrature
!
!
!     Using add and subtract - no need to call tree and set fmm parameters
!      can directly call existing fmm library
!
!
!       input:
!         npatches - integer
!            number of patches
!
!         norders- integer(npatches)
!            order of discretization on each patch 
!
!         ixyzs - integer(npatches+1)
!            ixyzs(i) denotes the starting location in srccoefs,
!               and srcvals array corresponding to patch i
!   
!         iptype - integer(npatches)
!            type of patch
!             iptype = 1, triangular patch discretized using RV nodes
!
!         npts - integer
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
!         ndtarg - integer
!            leading dimension of target array
!        
!         ntarg - integer
!            number of targets
!
!         targs - real *8 (ndtarg,ntarg)
!            target information
!
!          eps - real *8
!             precision requested
!
!          zpars - complex *16 (5)
!              kernel parameters (Referring to formula (1))
!              zpars(1) = omega 
!              zpars(2) = ep0
!              zpars(3) = mu0
!              zpars(4) = ep
!              zpars(5) = mu
!
!           nnz - integer *8
!             number of source patch-> target interactions in the near field
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
!           nquad - integer
!               number of entries in wnear
!
!           wnear - complex *16(nquad)
!               the near field quadrature correction
!
!           sigma - complex *16(2*npts)
!               density for layer potentials (\rho and \lambda)
!
!           novers - integer(npatches)
!              order of discretization for oversampled sources and
!               density
!
!         ixyzso - integer(npatches+1)
!            ixyzso(i) denotes the starting location in srcover,
!               corresponding to patch i
!   
!           nptso - integer
!              total number of oversampled points
!
!           srcover - real *8 (12,nptso)
!              oversampled set of source information
!
!           whtsover - real *8 (nptso)
!             smooth quadrature weights at oversampled nodes
!
!
!       output:
!         pot- complex *16(2*ntarg)
!            jump in u and (1/ep)u' 
!           
!               
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
      complex *16 wnear(nquad),sigma(2*npts),sigma2(npts)
	  
	  complex *16 wnear11(nquad),wnear12(nquad)
	  complex *16 wnear21(nquad),wnear22(nquad)
	  
      integer novers(npatches+1)
      integer nover,npolso,nptso
      real *8 srcover(12,nptso),whtsover(nptso)
      complex *16 pot(2*ntarg)
      complex *16, allocatable :: potsort(:)

      real *8, allocatable :: sources(:,:),targvals(:,:),wtmp2(:)
      complex *16, allocatable :: charges(:),dipvec(:,:),sigmaover(:)
      integer ns,nt
      complex *16 alpha,beta
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      complex *16 tmp(10),val,E(2)

      real *8 xmin,xmax,ymin,ymax,zmin,zmax,sizey,sizez,boxsize


      integer i,j,jpatch,jquadstart,jstart,count1,count2


      integer ifaddsub

      integer ntj
      
      complex *16 zdotu,pottmp
      complex *16, allocatable :: ctmp2_u(:),ctmp2_v(:),ctmp2_s(:),dtmp2(:,:)
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
      allocate(sigmaover(2*ns))

! 
!       oversample density \rho and \lambda
!

      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
     &npts,sigma(1:npts),novers,ixyzso,ns,sigmaover(1:ns))
	       
	call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
     &npts,sigma(npts+1:2*npts),novers,ixyzso,ns,sigmaover(ns+1:2*ns))



      ra = 0


!
!       fmm
!
	  
		!Calculate the far_field with FMM
		call transmission_FMM(eps,zpars,npts,targs,srcover,whtsover,&
	   &sigmaover(1:ns),sigmaover(ns+1:2*ns),pot(1:ntarg),pot(ntarg+1:2*ntarg),ns)
!        compute threshold for ignoring local computation
		call get_thresh(srcover,ns,targs,ntarg,thresh)

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
		  	pot(i) = pot(i) + wnear11(jquadstart+l-1)*sigma(jstart+l-1)
			pot(i) = pot(i) + wnear12(jquadstart+l-1)*sigma(jstart+l-1+npts)
			pot(i+npts) = pot(i+npts) + wnear21(jquadstart+l-1)*sigma(jstart+l-1)
			pot(i+npts) = pot(i+npts) + wnear22(jquadstart+l-1)*sigma(jstart+l-1+npts)
          enddo
        enddo
      enddo
	  
!C$OMP END PARALLEL DO

!
!     Remove near contribution of the FMM
!

!C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
!C$OMP$PRIVATE(ctmp2_u,ctmp2_v,wtmp2,nss,l,jstart,ii,E,npover)
!		ipars(1)=1
!		ipars(2)=1
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

     
		call fker_h_transmission(nss,ntarg0,srctmp2, targs(:,i), dpars,zpars&
	   &,ipars,ctmp2_u,ctmp2_v,wtmp2,E,thresh)
	    pot(i) = pot(i) - E(1)
	    pot(i+ntarg) = pot(i+ntarg) - E(2)
		deallocate(srctmp2,ctmp2_u,ctmp2_v,wtmp2)
      enddo
		do j=1,ntarg
				pot(j+ntarg)=-pot(j+ntarg)
		enddo

      return
      end subroutine lpcomp_h_transmission_addsub





subroutine transmission_FMM(eps,zpars,ns,srcvals,srcover,wts,rho_in,lambda_in,PHI,PSI,ns_over)
implicit none

    !List of calling arguments
	real ( kind = 8 ), intent(in) :: eps
	complex ( kind = 8 ), intent(in) :: zpars(5)
    integer, intent(in) :: ns,ns_over
    real ( kind = 8 ), intent(in) :: srcvals(12,ns)
    real ( kind = 8 ), intent(in) :: wts(ns_over),srcover(12,ns_over)
    complex ( kind = 8 ), intent(in) :: rho_in(ns_over),lambda_in(ns_over)
    complex ( kind = 8 ), intent(out) :: PHI(ns),PSI(ns)

    !List of local variables
    complex ( kind = 8 ), allocatable :: charge(:),dipvec(:,:)
    complex ( kind = 8 ), allocatable :: pot(:),gradpot(:,:)
	complex ( kind = 8 ) ima,zk0,zk1,ep0,mu0,ep1,mu1,omega

    integer count1,count2
	real ( kind = 8 ) pi
	pi=3.1415926535897932384626433832795028841971d0

	
	ima=(0.0d0,1.0d0)
    omega=zpars(1)
	ep0=zpars(2)
	mu0=zpars(3)
	ep1=zpars(4)
	mu1=zpars(5)

	zk0=omega*sqrt(ep0*mu0)
	zk1=omega*sqrt(ep1*mu1)
	
    allocate(charge(ns_over))
	allocate(dipvec(3,ns_over))
    allocate(gradpot(3,ns))
    allocate(pot(ns))

    do count1=1,ns_over
		charge(count1)=lambda_in(count1)*ep0**2*wts(count1)
		dipvec(1,count1)=srcover(10,count1)*ep0*rho_in(count1)*wts(count1)
		dipvec(2,count1)=srcover(11,count1)*ep0*rho_in(count1)*wts(count1)
		dipvec(3,count1)=srcover(12,count1)*ep0*rho_in(count1)*wts(count1)
	enddo

    !Computing the full operator
	call hfmm3d_t_cd_g(eps,zk0,ns_over,srcover(1:3,:),charge,dipvec,ns,srcvals(1:3,:),pot,gradpot)

	do count1=1,ns
		PHI(count1)=pot(count1)
		PSI(count1)=(srcvals(10,count1)*gradpot(1,count1)+srcvals(11,count1)&
		 &*gradpot(2,count1)+srcvals(12,count1)*gradpot(3,count1))/ep0
	enddo
	
    do count1=1,ns_over
		charge(count1)=lambda_in(count1)*ep1**2*wts(count1)
		dipvec(1,count1)=srcover(10,count1)*ep1*rho_in(count1)*wts(count1)
		dipvec(2,count1)=srcover(11,count1)*ep1*rho_in(count1)*wts(count1)
		dipvec(3,count1)=srcover(12,count1)*ep1*rho_in(count1)*wts(count1)
	enddo

    !Computing the full operator
	call hfmm3d_t_cd_g(eps,zk1,ns_over,srcover(1:3,:),charge,dipvec,ns,srcvals(1:3,:),pot,gradpot)

	do count1=1,ns
		PHI(count1)=PHI(count1)-pot(count1)
		PSI(count1)=PSI(count1)-(srcvals(10,count1)*gradpot(1,count1)+&
		 &srcvals(11,count1)*gradpot(2,count1)+srcvals(12,count1)*gradpot(3,count1))/ep1
	enddo


	do count1=1,ns
		PHI(count1)=PHI(count1)/(4.0d0*pi)
		PSI(count1)=PSI(count1)/(4.0d0*pi)	
	enddo


	deallocate(charge)
	deallocate(dipvec)
	deallocate(pot)
	deallocate(gradpot)

return
end subroutine transmission_FMM







subroutine test_accuracy_transmission(eps_FMM,sol,zpars,ns,wts,srcvals,P0,Pt)
implicit none

!
!	 This function test the accuracy in the interior and exterior region.
!

    !List of calling arguments
	complex ( kind = 8 ), intent(in) :: zpars(5)
    integer, intent(in) :: ns
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),eps_FMM
    real ( kind = 8 ), intent(in) :: wts(ns),P0(3),Pt(3)
	complex ( kind = 8 ), intent(in) :: sol(2*ns)
	
    !List of local variables
	complex ( kind = 8 ) omega,ep0,mu0,ep1,mu1
	complex ( kind = 8 ) ima, pot_t1(1),pot_t2(1),pot_01(1), pot_02(1),gradpot_aux(3),zk0,zk1,qv
	real ( kind = 8 ) pi
	integer count1


		ima=(0.0d0,1.0d0)
		pi=3.1415926535897932384626433832795028841971d0
		omega = zpars(1)
	    ep0 = zpars(2)
	    mu0 = zpars(3)
	    ep1 = zpars(4)
	    mu1 = zpars(5)
		zk0=omega*sqrt(ep0*mu0)
		zk1=omega*sqrt(ep1*mu1)
		call transmission_FMM_targ(eps_FMM,omega,ep0,mu0,ns,srcvals,1,Pt,wts,sol(1:ns),sol(ns+1:2*ns),pot_t1)
		call transmission_FMM_targ(eps_FMM,omega,ep1,mu1,ns,srcvals,1,P0,wts,sol(1:ns),sol(ns+1:2*ns),pot_01)
		
			 qv=1.0d0
		call point_source_scalar_helmholtz2(P0,qv,1,Pt,zk0,pot_t2,gradpot_aux,0)
		call point_source_scalar_helmholtz2(Pt,qv,1,P0,zk1,pot_02,gradpot_aux,0)
		

		write (*,*) 'Errors at the EXTERIOR region:'
		write (*,*) 'Relative Error in E: ', abs(pot_t1(1)-pot_t2(1))/abs(pot_t2(1))
		
		write (*,*) 'Errors at the Interior region:'
		write (*,*) 'Relative Error in E: ', abs(pot_01(1)-pot_02(1))/abs(pot_02(1))

return
end subroutine test_accuracy_transmission






subroutine transmission_FMM_targ(eps,omega,ep,mu,ns,srcvals,nt,targ,wts,rho_in,lambda_in,pot)
implicit none

	!
	!	This funciton computes the potential u (u0 or u1) at a given point using the FMM
	!   It doesn't contain near field corrections (it's for debugging purposes)   
	!   The representation for the potentials is:
	!      Representation:
	!        u1 = ep1 D_{k}\rho + ep1^2 S_{k}\lambda
	!        u0 = ep0 D_{k}\rho + ep0^2 S_{k}\lambda
	!   It's up to the user to evaluate each representation in the appropriate region
	!

    !List of calling arguments
	real ( kind = 8 ), intent(in) :: eps
	complex ( kind = 8 ), intent(in) :: omega,ep,mu
    integer, intent(in) :: ns,nt
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),targ(3,nt)
    real ( kind = 8 ), intent(in) :: wts(ns)
    complex ( kind = 8 ), intent(in) :: rho_in(ns),lambda_in(ns)
    complex ( kind = 8 ), intent(out) :: pot(nt)

    !List of local variables
	real ( kind = 8 ), allocatable :: n_vect(:,:),source(:,:)
    complex ( kind = 8 ), allocatable :: charge(:),dipvec(:,:)
	complex ( kind = 8 ) ima,zk

    integer count1,count2
    integer ifa_vect,ifb_vect,iflambda,ifrho,ifE,ifcurlE,ifdivE
	real ( kind = 8 ) pi
	pi=3.1415926535897932384626433832795028841971d0


	ima=(0.0d0,1.0d0)
	zk=omega*sqrt(ep*mu)

	allocate(charge(ns))
	allocate(dipvec(3,ns))
	
	allocate(n_vect(3,ns))
	allocate(source(3,ns))


	do count1=1,ns
		n_vect(:,count1)=srcvals(10:12,count1)
		source(:,count1)=srcvals(1:3,count1)
	enddo
	
	do count1=1,ns
		charge(count1)=lambda_in(count1)*ep**2*wts(count1)
		dipvec(1,count1)=n_vect(1,count1)*ep*rho_in(count1)*wts(count1)
		dipvec(2,count1)=n_vect(2,count1)*ep*rho_in(count1)*wts(count1)
		dipvec(3,count1)=n_vect(3,count1)*ep*rho_in(count1)*wts(count1)
	enddo

    !Computing the full operator
    call hfmm3d_t_cd_p(eps,zk,ns,source,charge,dipvec,nt,targ,pot)
	do count1=1,nt
		pot(count1)=pot(count1)/(4*pi)
	enddo


	deallocate(charge)
	deallocate(dipvec)
	deallocate(n_vect)
	deallocate(source)

return
end subroutine transmission_FMM_targ




subroutine 	get_RHS_transmission(P0,Pt,ns,srcvals,zpars,RHS)
implicit none

	!List of calling arguments
	integer ( kind = 4 ), intent(in) :: ns
	real ( kind = 8 ), intent(in) :: P0(3),Pt(3)
	real ( kind = 8 ), intent(in) :: srcvals(12,ns)
	complex (kind = 8 ), intent(in) :: zpars(5)
	complex ( kind = 8 ), intent(out) :: RHS(2*ns)
	
	!List of local variables
	complex ( kind = 8 ) omega, ep0,mu0,ep1,mu1
	complex ( kind = 8 ), allocatable :: gradpot(:,:),pot(:)
	integer count1
	complex ( kind = 8 ) ima,zk0,zk1,qv
	real ( kind = 8 ) ru(3),rv(3),cross_aux(3)
		allocate(gradpot(3,ns),pot(ns))

	ima=(0.0d0,1.0d0)
	
	omega = zpars(1)
	ep0 = zpars(2)
	mu0 = zpars(3)
	ep1 = zpars(4)
	mu1 = zpars(5)

	zk0=omega*sqrt(ep0*mu0)
	zk1=omega*sqrt(ep1*mu1)
		
	qv=1.0d0
	call point_source_scalar_helmholtz2(P0,qv,ns,srcvals(1:3,:),zk0,pot,gradpot,0)	
	
	do count1=1,ns
		RHS(count1)=pot(count1)
		RHS(ns+count1)=(1.0d0/ep0)*DOT_PRODUCT(srcvals(10:12,count1),gradpot(:,count1))
	enddo
	
	call point_source_scalar_helmholtz2(Pt,qv,ns,srcvals(1:3,:),zk1,pot,gradpot,0)
	
	do count1=1,ns
		RHS(count1)=RHS(count1)-pot(count1)	
	RHS(ns+count1)=RHS(ns+count1)-(1.0d0/ep1)*DOT_PRODUCT(srcvals(10:12,count1),gradpot(:,count1))
		RHS(ns+count1)=-RHS(ns+count1)
	enddo
	
	
!	do count1=1,ns
!		RHS(ns+count1)=RHS(ns+count1)/omega
!	enddo
	
	
	deallocate(pot)
	deallocate(gradpot)

return
end subroutine get_RHS_transmission



