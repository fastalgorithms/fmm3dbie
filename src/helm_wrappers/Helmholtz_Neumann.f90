      subroutine h_neumann_solver(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,eps,zpars,numit,ifinout,&
     &rhs,eps_gmres,niter,errs,rres,soln,sigma2)
!
!
!        this subroutine solves the helmholtz Neumann problem
!     on the exterior of an object where the potential
!     is represented as a combined field integral equation.
!
!
!     Representation:
!        u = S_{k}[\rho]+i*alpha*D_{k}[S_{ik}[\rho]]
!
!     Boundary condition:
!        u'=f
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
!          zpars - complex *16 (3)
!              kernel parameters (Referring to formula (1))
!              zpars(1) = k 
!              zpars(2) = alpha
!              zpars(3) =  - (not used)
!
!          ifinout - integer
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
!           numit - integer
!              max number of gmres iterations
!
!         output
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
!           soln - complex *16(npts)
!              density which solves the dirichlet problem \rho
!
!			sigma2 - complex *16(npts)
!			   sigma2=S_{ik}[\rho] this is also output to make live 
!			     easyer when evalutating the far field from the sources
!				 
!
!
      implicit none
      integer npatches,norder,npols,npts
      integer ifinout
      integer norders(npatches),ixyzs(npatches+1)
      integer iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,eps_gmres
      complex *16 zpars(3)
      complex *16 rhs(npts)
      complex *16 soln(npts),sigma2(npts)

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

      call getnearquad_h_neumann(npatches,norders,&
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
!      In this case it is a bit messy because of the identity that
!      comes from the calderon identity
	  
	  zid=(-1.0d0/2.0d0-ima*zpars(2)/4.0d0)
	  
!     write (*,*) 'zid', zid,ifinout,zpars(3)


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


        call lpcomp_h_neumann_addsub(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,&
     &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,&
     &vmat(1,it),novers,npts_over,ixyzso,srcover,wover,wtmp,sigma2,&
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


          call lpcomp_h_neumann_addsub(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,&
     &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,&
     &soln,novers,npts_over,ixyzso,srcover,wover,wtmp,sigma2,&
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
      end subroutine h_neumann_solver








      subroutine getnearquad_h_neumann(npatches,norders,&
     &ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,&
     &iquad,rfac0,nquad,wnear11,wnear12,wnear21,wnear22)
!
!       this subroutine generates the near field quadrature
!       for the representation:
!        u = S_{k}[\rho]+i*alpha*D_{k}[S_{ik}[\rho]]
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
!          zpars - complex *16 (3)
!              kernel parameters (Referring to formula (1))
!              zpars(1) = k 
!              zpars(2) = alpha
!              zpars(3) =  - (not used)
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
      complex *16 zpars(3)
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
	  external  fker_h_neumanns

     ndz=3
	 ndd=1
	 ndi=2
	 ipv=1

      fker =>  fker_h_neumanns
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
      end subroutine getnearquad_h_neumann


subroutine fker_h_neumanns(srcinfo, ndt,targinfo,ndd, dpars,ndz,zpars,ndi,ipars,E_val)
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
	complex ( kind = 8 ) Dk1rho(1,1),Dk0rho(1,1),Sk1lambda(1,1),Sk0lambda(1,1),ngradSk1lambda(1,1),ngradSk0lambda(1,1)
	complex ( kind = 8 ) ngradDk1rho(1,1),ngradDk0rho(1,1),E_mat(2,2)
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3),xprod_aux3(3),xprod_aux4(3)
	complex ( kind = 8 ) R1_0,R1_1,R2_0,R2_1,ima,my_exp_0,my_exp_1, zk,alpha,zk0,zk1
	complex ( kind = 8 ) R1_0_stab,R1_1_stab,R2_0_stab,R2_1_stab
	real ( kind = 8 ) pi
	integer count1,count2
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
	
	call orthonormalize(srcinfo(4:6),n_s,ru_s,rv_s)
	call orthonormalize(targinfo(4:6),n_t,ru_t,rv_t)
	
	if (ipars(1).eq.1) then
	  if (ipars(2).eq.1) then
		zk0=zk
	  	R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
		my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
	    call get_ngradSklambda(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_0,my_exp_0,r,ngradSk0lambda)
	    E_val=ngradSk0lambda(1,1)
	  else
		zk1=zk*ima
	    my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
		call get_Sklambda(my_exp_1,r,Sk1lambda)
	    E_val=Sk1lambda(1,1)
	  endif
	else
	  if (ipars(2).eq.1) then
		zk1=zk*ima
	  	R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
		my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
		call get_ngradSklambda(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_1,my_exp_1,r,ngradSk1lambda)
		E_val=ngradSk1lambda(1,1)
	  else
        zk0=zk
        zk1=zk*ima
        my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
        my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
        R1_0_stab=(1.0d0/(r**3*4.0d0*pi))*((ima*zk0*r-1.0d0)*(exp(ima*zk0*r/2.0d0)*2*ima*sin(zk0*r/2.0d0))+ima*zk0*r)
        R2_0_stab=(1.0d0/(4.0d0*pi*r**5))*(((ima*zk0*r)**2-3.0d0*ima*zk0*r+3.0d0)*(exp(ima*zk0*r/2.0d0)*&
        &2*ima*sin(zk0*r/2.0d0))+(ima*zk0*r)**2-3.0d0*ima*zk0*r)
        R1_1_stab=(1.0d0/(r**3*4.0d0*pi))*((ima*zk1*r-1.0d0)*(exp(ima*zk1*r/2.0d0)*2*ima*sin(zk1*r/2.0d0))+ima*zk1*r)
        R2_1_stab=(1.0d0/(4.0d0*pi*r**5))*(((ima*zk1*r)**2-3.0d0*ima*zk1*r+3.0d0)*(exp(ima*zk1*r/2.0d0)*&
        &2*ima*sin(zk1*r/2.0d0))+(ima*zk1*r)**2-3.0d0*ima*zk1*r)
        call get_ngradDkrho(n_s,n_t,dr,R1_0_stab,R2_0_stab,zk0,my_exp_0,r,ngradDk0rho)
        call get_ngradDkrho(n_s,n_t,dr,R1_1_stab,R2_1_stab,zk1,my_exp_1,r,ngradDk1rho)
        E_val=ngradDk0rho(1,1)-ngradDk1rho(1,1)
	  endif	
	endif

return
end subroutine fker_h_neumanns


subroutine fker_h_neumanns_not_used(srcinfo, ndt,targinfo,ndd, dpars,ndz,zpars,ndi,ipars,E_val)
implicit none

! this function provides the near field kernel that will use zgetnearquad_ggq_guru inside
! getnearquad_h_neumann

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
	complex ( kind = 8 ) Dk1rho(1,1),Dk0rho(1,1),Sk1lambda(1,1),Sk0lambda(1,1),ngradSk1lambda(1,1),ngradSk0lambda(1,1)
	complex ( kind = 8 ) ngradDk1rho(1,1),ngradDk0rho(1,1),E_mat(2,2)
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3),xprod_aux3(3),xprod_aux4(3)
	complex ( kind = 8 ) R1_0,R1_1,R2_0,R2_1,ima,my_exp_0,my_exp_1, zk,alpha,zk0,zk1
	complex ( kind = 8 ) R1_0_stab,R1_1_stab,R2_0_stab,R2_1_stab
	real ( kind = 8 ) pi
	integer count1,count2
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
	
	zk0=zk
	zk1=zk*ima
	
	R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
	R2_0=((ima*zk0)**2/r**3-3.0d0*ima*zk0/r**4+3.0d0/r**5)*exp(ima*zk0*r)/(4.0d0*pi)
	my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
	R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
	R2_1=((ima*zk1)**2/r**3-3.0d0*ima*zk1/r**4+3.0d0/r**5)*exp(ima*zk1*r)/(4.0d0*pi)
	my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
			
	call orthonormalize(srcinfo(4:6),n_s,ru_s,rv_s)
	call orthonormalize(targinfo(4:6),n_t,ru_t,rv_t)

	call get_ngradSklambda(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_0,my_exp_0,r,ngradSk0lambda)
	E_mat(1,1)=ngradSk0lambda(1,1)
	call get_ngradSklambda(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_1,my_exp_1,r,ngradSk1lambda)
	E_mat(2,1)=ngradSk1lambda(1,1)
	call get_Sklambda(my_exp_1,r,Sk1lambda)
	E_mat(1,2)=Sk1lambda(1,1)
			
			
	R1_0_stab=(1.0d0/(r**3*4.0d0*pi))*((ima*zk0*r-1.0d0)*(exp(ima*zk0*r/2.0d0)*2*ima*sin(zk0*r/2.0d0))+ima*zk0*r)
	R2_0_stab=(1.0d0/(4.0d0*pi*r**5))*(((ima*zk0*r)**2-3.0d0*ima*zk0*r+3.0d0)*(exp(ima*zk0*r/2.0d0)*&
	&2*ima*sin(zk0*r/2.0d0))+(ima*zk0*r)**2-3.0d0*ima*zk0*r)
						
	R1_1_stab=(1.0d0/(r**3*4.0d0*pi))*((ima*zk1*r-1.0d0)*(exp(ima*zk1*r/2.0d0)*2*ima*sin(zk1*r/2.0d0))+ima*zk1*r)
	R2_1_stab=(1.0d0/(4.0d0*pi*r**5))*(((ima*zk1*r)**2-3.0d0*ima*zk1*r+3.0d0)*(exp(ima*zk1*r/2.0d0)*&
	&2*ima*sin(zk1*r/2.0d0))+(ima*zk1*r)**2-3.0d0*ima*zk1*r)

	call get_ngradDkrho(n_s,n_t,dr,R1_0_stab,R2_0_stab,zk0,my_exp_0,r,ngradDk0rho)
	call get_ngradDkrho(n_s,n_t,dr,R1_1_stab,R2_1_stab,zk1,my_exp_1,r,ngradDk1rho)
			

	E_mat(2,2)=ngradDk0rho(1,1)-ngradDk1rho(1,1)
	E_val=E_mat(ipars(1),ipars(2))
	

return
end subroutine fker_h_neumanns_not_used


subroutine fker_h_neumann(ns,nt,srcinfo, targinfo, dpars,zpars,ipars,sigma,wts,E_val,thresh)
implicit none

! this function does something very similar to fker_h_neumanns but oriented to 
! remove the near interaction of the FMM 

    !List of calling arguments
	integer, intent(in) :: ns, nt
	real ( kind = 8 ), intent(in) :: srcinfo(12,ns)
	real ( kind = 8 ), intent(in) :: targinfo(12,nt)
	integer, intent(in) :: ipars(2)
	real ( kind = 8 ), intent(in) :: dpars(1)
	complex ( kind = 8 ), intent(in) :: zpars(3),sigma(ns)
	complex ( kind = 8 ), intent(out) :: E_val(nt)
    real ( kind = 8 ), intent(in) :: thresh,wts(ns)
	
	!List of local variables
	real ( kind = 8 ) ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ) ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ) du(3), dv(3),sour(3),targ(3)
	real ( kind = 8 ) r, dr(3),aux_real
	complex ( kind = 8 ) E_mat(2,2)
	complex ( kind = 8 ) Dk1rho(1,1),Dk0rho(1,1),Sk1lambda(1,1),Sk0lambda(1,1),ngradSk1lambda(1,1),ngradSk0lambda(1,1)
	complex ( kind = 8 ) ngradDk1rho(1,1),ngradDk0rho(1,1)
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3),xprod_aux3(3),xprod_aux4(3)
	complex ( kind = 8 ) R1_0,R1_1,R2_0,R2_1,ima,my_exp_0,my_exp_1, zk,alpha,zk0,zk1
	complex ( kind = 8 ) R1_0_stab,R1_1_stab,R2_0_stab,R2_1_stab
	real ( kind = 8 ) pi
	integer count1,count2
	
	pi=3.1415926535897932384626433832795028841971d0
	ima=(0.0d0,1.0d0)
	zk=zpars(1)
	alpha=zpars(2)
	do count1=1,nt
	  E_val(count1)=0.0d0
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
		  zk0=zk
		  zk1=zk*ima

		  R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
		  R2_0=((ima*zk0)**2/r**3-3.0d0*ima*zk0/r**4+3.0d0/r**5)*exp(ima*zk0*r)/(4.0d0*pi)
		  my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
		  R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
		  R2_1=((ima*zk1)**2/r**3-3.0d0*ima*zk1/r**4+3.0d0/r**5)*exp(ima*zk1*r)/(4.0d0*pi)
		  my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
			
		  call orthonormalize(srcinfo(4:6,count2),n_s,ru_s,rv_s)
		  call orthonormalize(targinfo(4:6,count1),n_t,ru_t,rv_t)

		  call get_ngradSklambda(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_0,my_exp_0,r,ngradSk0lambda)
		  E_mat(1,1)=ngradSk0lambda(1,1)
		  call get_ngradSklambda(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_1,my_exp_1,r,ngradSk1lambda)
		  E_mat(2,1)=ngradSk1lambda(1,1)
		  call get_Sklambda(my_exp_1,r,Sk1lambda)
		  E_mat(1,2)=Sk1lambda(1,1)
			
!	R1_0_stab=(1.0d0/(r**3*4.0d0*pi))*((ima*zk0*r-1.0d0)*(exp(ima*zk0*r/2.0d0)*2*ima*sin(zk0*r/2.0d0))+ima*zk0*r)
!			R2_0_stab=(1.0d0/(4.0d0*pi*r**5))*(((ima*zk0*r)**2-3.0d0*ima*zk0*r+3.0d0)*(exp(ima*zk0*r/2.0d0)*&
!			&2*ima*sin(zk0*r/2.0d0))+(ima*zk0*r)**2-3.0d0*ima*zk0*r)
						
!	R1_1_stab=(1.0d0/(r**3*4.0d0*pi))*((ima*zk1*r-1.0d0)*(exp(ima*zk1*r/2.0d0)*2*ima*sin(zk1*r/2.0d0))+ima*zk1*r)
!			R2_1_stab=(1.0d0/(4.0d0*pi*r**5))*(((ima*zk1*r)**2-3.0d0*ima*zk1*r+3.0d0)*(exp(ima*zk1*r/2.0d0)*&
!			&2*ima*sin(zk1*r/2.0d0))+(ima*zk1*r)**2-3.0d0*ima*zk1*r)

!		call get_ngradDkrho(n_s,n_t,dr,R1_0_stab,R2_0_stab,zk0,my_exp_0,r,ngradDk0rho)
!		call get_ngradDkrho(n_s,n_t,dr,R1_1_stab,R2_1_stab,zk1,my_exp_1,r,ngradDk1rho)
			
!!!			Here is where I compute D'k alone so that I can remove the effect of the FMM near contribution
			call get_ngradDkrho(n_s,n_t,dr,R1_0,R2_0,zk0,my_exp_0,r,ngradDk0rho)
!			call get_ngradDkrho(n_s,n_t,dr,R1_1,R2_1,zk1,my_exp_1,r,ngradDk1rho)
			E_mat(2,2)=ngradDk0rho(1,1)!-ngradDk1rho(1,1)

		else

			E_mat(1,1)=0.0d0
			E_mat(1,2)=0.0d0
			E_mat(2,1)=0.0d0
			E_mat(2,2)=0.0d0

		endif
	    E_val(count1)=E_val(count1)+E_mat(ipars(1),ipars(2))*sigma(count2)*wts(count2)
	  enddo
	enddo
return
end subroutine fker_h_neumann



      subroutine lpcomp_h_neumann_addsub(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,sigma,novers,&
     &nptso,ixyzso,srcover,whtsover,pot,sigma2,wnear11,wnear12,&
	 &wnear21,wnear22)

!
!      this subroutine evaluates the layer potential for
!      the representation:
!        u = S_{k}[\rho]+i*alpha*D_{k}[S_{ik}[\rho]]
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
!          zpars - complex *16 (3)
!              kernel parameters (Referring to formula (1))
!              zpars(1) = k 
!              zpars(2) = alpha
!              zpars(3) =  - (not used)
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
!               location in wnear_ij array where quadrature for col_ind(i)
!               starts
!
!           nquad - integer
!               number of entries in wnear_ij
!
!           wnear_ij - complex *16(nquad)
!               the near field quadrature correction for the different
!               operators that apear in the formulation
!
!           sigma - complex *16(npts)
!               density for layer potential
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
      complex *16 zpars(3),zpars2(3)
      integer nnz,row_ptr(ntarg+1),col_ind(nnz),nquad
      integer iquad(nnz+1)
      complex *16 sigma(npts),sigma2(npts)
	  
	  complex *16 wnear11(nquad),wnear12(nquad)
	  complex *16 wnear21(nquad),wnear22(nquad)
	  
      integer novers(npatches+1)
      integer nover,npolso,nptso
      real *8 srcover(12,nptso),whtsover(nptso)
      complex *16 pot(ntarg)
      complex *16, allocatable :: potsort(:)

      real *8, allocatable :: sources(:,:),targvals(:,:),wtmp2(:)
      complex *16, allocatable :: charges(:),dipvec(:,:),sigmaover(:)
	  complex *16, allocatable :: sigmaover_aux(:),sigma_aux(:)
      integer ns,nt
      complex *16 alpha,beta
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      complex *16 tmp(10),val,E(1)

      real *8 xmin,xmax,ymin,ymax,zmin,zmax,sizey,sizez,boxsize


      integer i,j,jpatch,jquadstart,jstart,count1,count2


      integer ifaddsub

      integer ntj
      
      complex *16 zdotu,pottmp
      complex *16, allocatable :: ctmp2_u(:),ctmp2_v(:),ctmp2_s(:),dtmp2(:,:)
	  complex *16, allocatable :: pot_aux(:),pot_aux2(:)

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
      allocate(sigmaover(ns),sigmaover_aux(ns))
	  allocate(pot_aux(ntarg),pot_aux2(ntarg),sigma_aux(ntarg))



! 
!       oversample density
!
      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
     &npts,sigma(1:npts),novers,ixyzso,ns,sigmaover(1:ns))


      ra = 0


!
!       fmm
!
	  
		!NOW WE COMPUTE Sk'(rho) and add to PHI
		!Calculate the far_field with FMM
		call h_neumann_FMM_11(eps,zpars(1),npts,targs,srcover,whtsover,sigmaover,pot,ns)
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
          enddo
        enddo
      enddo
!C$OMP END PARALLEL DO

!
!     Remove near contribution of the FMM
!

!C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
!C$OMP$PRIVATE(ctmp2_u,ctmp2_v,wtmp2,nss,l,jstart,ii,E,npover)
		ipars(1)=1
		ipars(2)=1
      do i=1,ntarg
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          nss = nss + ixyzso(jpatch+1)-ixyzso(jpatch)
        enddo
        allocate(srctmp2(12,nss),ctmp2_u(nss),wtmp2(nss))

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
			wtmp2(ii)=whtsover(jstart+l)
          enddo
        enddo

call fker_h_neumann(nss,ntarg0,srctmp2, targs(:,i), dpars,zpars,ipars,ctmp2_u,wtmp2,E,thresh)

	    pot(i) = pot(i) - E(1)
		deallocate(srctmp2,ctmp2_u,wtmp2)
      enddo

      call cpu_time(t2 )
!C$      t2 = omp_get_wtime()     

      timeinfo(2) = t2-t1
	  

!!      call prin2('quadrature time=*',timeinfo,2)
      
      ttot = timeinfo(1) + timeinfo(2)
!!      call prin2('time in lpcomp=*',ttot,1)


		!NOW WE COMPUTE Sik'(rho) and add to PHI
		!Calculate the far_field with FMM
call h_neumann_FMM_11(eps,ima*zpars(1),npts,targs,srcover,whtsover,sigmaover,pot_aux,ns)

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
		    pot_aux(i) = pot_aux(i) + wnear21(jquadstart+l-1)*sigma(jstart+l-1)
          enddo
        enddo
      enddo
!C$OMP END PARALLEL DO

!
!     Remove near contribution of the FMM
!

!C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
!C$OMP$PRIVATE(ctmp2_u,ctmp2_v,wtmp2,nss,l,jstart,ii,E,npover)

		ipars(1)=2
		ipars(2)=1
      do i=1,ntarg
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          nss = nss + ixyzso(jpatch+1)-ixyzso(jpatch)
        enddo
        allocate(srctmp2(12,nss),ctmp2_u(nss),wtmp2(nss))

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
			wtmp2(ii)=whtsover(jstart+l)
          enddo
        enddo

call fker_h_neumann(nss,ntarg0,srctmp2, targs(:,i), dpars,zpars,ipars,ctmp2_u,wtmp2,E,thresh)

	    pot_aux(i) = pot_aux(i) - E(1)
		deallocate(srctmp2,ctmp2_u,wtmp2)
      enddo
	  
      
      call cpu_time(t2)
!C$      t2 = omp_get_wtime()     

      timeinfo(2) = t2-t1


!!      call prin2('quadrature time=*',timeinfo,2)
      
      ttot = timeinfo(1) + timeinfo(2)
!!      call prin2('time in lpcomp=*',ttot,1)


      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
     &npts,pot_aux(1:npts),novers,ixyzso,ns,sigmaover_aux(1:ns))

	 do i=1,ntarg
		sigma_aux(i)=pot_aux(i)
	 enddo


		!NOW WE COMPUTE Sik'(rho) and add to PHI
		!Calculate the far_field with FMM
call h_neumann_FMM_11(eps,ima*zpars(1),npts,targs,srcover,whtsover,sigmaover_aux,pot_aux,ns)
!        compute threshold for ignoring local computation
!      call get_thresh(srcover,ns,targs,ntarg,thresh)
!       add in precomputed quadrature
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
		    pot_aux(i) = pot_aux(i) + wnear21(jquadstart+l-1)*sigma_aux(jstart+l-1)
          enddo
        enddo
      enddo
!C$OMP END PARALLEL DO
!C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
!C$OMP$PRIVATE(ctmp2_u,ctmp2_v,wtmp2,nss,l,jstart,ii,E,npover)
		ipars(1)=2
		ipars(2)=1

      do i=1,ntarg
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          nss = nss + ixyzso(jpatch+1)-ixyzso(jpatch)
        enddo
        allocate(srctmp2(12,nss),ctmp2_u(nss),wtmp2(nss))

        rmin = 1.0d6
        ii = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          jstart = ixyzso(jpatch)-1
          npover = ixyzso(jpatch+1)-ixyzso(jpatch)
          do l=1,npover
			ii = ii+1
			srctmp2(:,ii) = srcover(:,jstart+l)
			ctmp2_u(ii)=sigmaover_aux(jstart+l)
			wtmp2(ii)=whtsover(jstart+l)
          enddo
        enddo

call fker_h_neumann(nss,ntarg0,srctmp2, targs(:,i), dpars,zpars,ipars,ctmp2_u,wtmp2,E,thresh)

	    pot_aux(i) = pot_aux(i) - E(1)
		deallocate(srctmp2,ctmp2_u,wtmp2)
      enddo
      
!      call cpu_time(t2)
!C$      t2 = omp_get_wtime()     

!      timeinfo(2) = t2-t1



!!      call prin2('quadrature time=*',timeinfo,2)
      
!      ttot = timeinfo(1) + timeinfo(2)
!!      call prin2('time in lpcomp=*',ttot,1)


		do count1=1,ntarg
			pot(count1)=pot(count1)+ima*zpars(2)*pot_aux(count1)
			pot_aux(count1)=0.0d0
		enddo


		!NOW WE COMPUTE Sik(rho) and put into pot_aux->rho_aux
		!Calculate the far_field with FMM
call h_neumann_FMM_12(eps,zpars(1),npts,targs,srcover,whtsover,sigmaover,pot_aux,ns)
		
!       add in precomputed quadrature
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
		    pot_aux(i) = pot_aux(i) + wnear12(jquadstart+l-1)*sigma(jstart+l-1)
          enddo
        enddo
      enddo
!C$OMP END PARALLEL DO
!C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
!C$OMP$PRIVATE(ctmp2_u,ctmp2_v,wtmp2,nss,l,jstart,ii,E,npover)
		ipars(1)=1
		ipars(2)=2
      do i=1,ntarg
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          nss = nss + ixyzso(jpatch+1)-ixyzso(jpatch)
        enddo
        allocate(srctmp2(12,nss),ctmp2_u(nss),wtmp2(nss))

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
			wtmp2(ii)=whtsover(jstart+l)
          enddo
        enddo

call fker_h_neumann(nss,ntarg0,srctmp2, targs(:,i), dpars,zpars,ipars,ctmp2_u,wtmp2,E,thresh)

	    pot_aux(i) = pot_aux(i) - E(1)
		deallocate(srctmp2,ctmp2_u,wtmp2)
      enddo

      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
     &npts,pot_aux(1:npts),novers,ixyzso,ns,sigmaover_aux(1:ns))

	 do i=1,ntarg
		sigma_aux(i)=pot_aux(i)
		sigma2(i)=pot_aux(i)
	 enddo

		!NOW WE COMPUTE Dk'(rho) and add to PHI
		!Calculate the far_field with FMM
		call h_neumann_FMM_22(eps,zpars(1),npts,targs,srcover,whtsover,sigmaover_aux,pot_aux,ns)

!C$OMP END PARALLEL DO
!C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
!C$OMP$PRIVATE(ctmp2_u,ctmp2_v,wtmp2,nss,l,jstart,ii,E,npover)
		ipars(1)=2
		ipars(2)=2
      do i=1,ntarg
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          nss = nss + ixyzso(jpatch+1)-ixyzso(jpatch)
        enddo
        allocate(srctmp2(12,nss),ctmp2_u(nss),wtmp2(nss))

        rmin = 1.0d6
        ii = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          jstart = ixyzso(jpatch)-1
          npover = ixyzso(jpatch+1)-ixyzso(jpatch)
          do l=1,npover
			ii = ii+1
			srctmp2(:,ii) = srcover(:,jstart+l)
			ctmp2_u(ii)=sigmaover_aux(jstart+l)
			wtmp2(ii)=whtsover(jstart+l)
          enddo
        enddo

call fker_h_neumann(nss,ntarg0,srctmp2, targs(:,i), dpars,zpars,ipars,ctmp2_u,wtmp2,E,thresh)

	    pot_aux(i) = pot_aux(i) - E(1)
		deallocate(srctmp2,ctmp2_u,wtmp2)
      enddo

		call h_neumann_FMM_22(eps,ima*zpars(1),npts,targs,srcover,whtsover,sigmaover_aux,pot_aux2,ns)
		
		!C$OMP END PARALLEL DO
!C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
!C$OMP$PRIVATE(ctmp2_u,ctmp2_v,wtmp2,nss,l,jstart,ii,E,npover)
		ipars(1)=2
		ipars(2)=2
      do i=1,ntarg
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          nss = nss + ixyzso(jpatch+1)-ixyzso(jpatch)
        enddo
        allocate(srctmp2(12,nss),ctmp2_u(nss),wtmp2(nss))

        rmin = 1.0d6
        ii = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          jstart = ixyzso(jpatch)-1
          npover = ixyzso(jpatch+1)-ixyzso(jpatch)
          do l=1,npover
			ii = ii+1
			srctmp2(:,ii) = srcover(:,jstart+l)
			ctmp2_u(ii)=sigmaover_aux(jstart+l)
			wtmp2(ii)=whtsover(jstart+l)
          enddo
        enddo

		zpars2(1)=zpars(1)*ima
		zpars2(2)=zpars(2)
		zpars2(3)=zpars(3)

		 call fker_h_neumann(nss,ntarg0,srctmp2, targs(:,i), dpars,zpars2,ipars,ctmp2_u,wtmp2,E,thresh)

	    pot_aux2(i) = pot_aux2(i) - E(1)
		deallocate(srctmp2,ctmp2_u,wtmp2)
      enddo


do i=1,ntarg
	pot_aux(i)=pot_aux(i)-pot_aux2(i)
enddo

!       add in precomputed quadrature

!      call cpu_time(t1)
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
		    pot_aux(i) = pot_aux(i) + wnear22(jquadstart+l-1)*sigma_aux(jstart+l-1)
          enddo
        enddo
      enddo
      
!      call cpu_time(t2)
!C$      t2 = omp_get_wtime()     

      timeinfo(2) = t2-t1

!!      call prin2('quadrature time=*',timeinfo,2)
      
      ttot = timeinfo(1) + timeinfo(2)
!!      call prin2('time in lpcomp=*',ttot,1)


		do count1=1,ntarg
			pot(count1)=pot(count1)+ima*zpars(2)*pot_aux(count1)
		enddo

      return
      end subroutine lpcomp_h_neumann_addsub



subroutine h_neumann_FMM_11(eps,zk,ns,srcvals,srcover,wts,rho_in,PHI,ns_over)
implicit none

!
!	One of the FMM operators needed in the formulation
!

    !List of calling arguments
	real ( kind = 8 ), intent(in) :: eps
	complex ( kind = 8 ), intent(in) :: zk
    integer, intent(in) :: ns, ns_over
    real ( kind = 8 ), intent(in) :: srcvals(12,ns), srcover(12,ns_over)
    real ( kind = 8 ), intent(in) :: wts(ns_over)
    complex ( kind = 8 ), intent(in) :: rho_in(ns_over)
    complex ( kind = 8 ), intent(out) :: PHI(ns)

    !List of local variables
    complex ( kind = 8 ), allocatable :: charge(:),dipvec(:,:)
    complex ( kind = 8 ), allocatable :: pot(:),gradpot(:,:)
	complex ( kind = 8 ) ima,zk0

    integer count1,count2
	real ( kind = 8 ) pi
	pi=3.1415926535897932384626433832795028841971d0

	ima=(0.0d0,1.0d0)
	zk0=zk
    allocate(charge(ns_over))
    allocate(gradpot(3,ns))
    allocate(pot(ns))


    do count1=1,ns_over
		charge(count1)=rho_in(count1)*wts(count1)
	enddo

    !Computing the full operator
	call hfmm3d_t_c_g(eps,zk,ns_over,srcover(1:3,:),charge,ns,srcvals(1:3,:),pot,gradpot)
	
	do count1=1,ns
		PHI(count1)=(srcvals(10,count1)*gradpot(1,count1)+srcvals(11,count1)&
		 &*gradpot(2,count1)+srcvals(12,count1)*gradpot(3,count1))/(4.0d0*pi)
	enddo

	deallocate(charge)
	deallocate(pot)
	deallocate(gradpot)

return
end subroutine h_neumann_FMM_11




subroutine h_neumann_FMM_12(eps,zk,ns,srcvals,srcover,wts,rho_in,PHI,ns_over)
implicit none

!
!	One of the FMM operators needed in the formulation
!

    !List of calling arguments
	real ( kind = 8 ), intent(in) :: eps
	complex ( kind = 8 ), intent(in) :: zk
    integer, intent(in) :: ns,ns_over
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),srcover(12,ns_over)
    real ( kind = 8 ), intent(in) :: wts(ns_over)
    complex ( kind = 8 ), intent(in) :: rho_in(ns_over)
    complex ( kind = 8 ), intent(out) :: PHI(ns)

    !List of local variables
!	real ( kind = 8 ), allocatable :: n_vect(:,:),source(:,:)
    complex ( kind = 8 ), allocatable :: charge(:),dipvec(:,:)
    complex ( kind = 8 ), allocatable :: pot(:),gradpot(:,:)
	complex ( kind = 8 ) ima,zk0,zk1

    integer count1,count2
	real ( kind = 8 ) pi
	pi=3.1415926535897932384626433832795028841971d0


	ima=(0.0d0,1.0d0)
	zk0=zk
	zk1=ima*zk

    allocate(charge(ns_over))
    allocate(pot(ns))

    do count1=1,ns_over
		charge(count1)=rho_in(count1)*wts(count1)
	enddo

    !Computing the full operator
	call hfmm3d_t_c_p(eps,zk1 ,ns_over,srcover(1:3,:),charge,ns,srcvals(1:3,:),pot)
	
	do count1=1,ns
		PHI(count1)=pot(count1)/(4.0d0*pi)
	enddo

	deallocate(charge)
	deallocate(pot)

return
end subroutine h_neumann_FMM_12



subroutine h_neumann_FMM_22(eps,zk,ns,srcvals,srcover,wts,rho_in,PHI,ns_over)
implicit none

!
!	One of the FMM operators needed in the formulation
!

    !List of calling arguments
	real ( kind = 8 ), intent(in) :: eps
	complex ( kind = 8 ), intent(in) :: zk
    integer, intent(in) :: ns,ns_over
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),srcover(12,ns_over)
    real ( kind = 8 ), intent(in) :: wts(ns_over)
    complex ( kind = 8 ), intent(in) :: rho_in(ns_over)
    complex ( kind = 8 ), intent(out) :: PHI(ns)

    !List of local variables
    complex ( kind = 8 ), allocatable :: dipvec(:,:)
    complex ( kind = 8 ), allocatable :: pot(:),gradpot(:,:)
	complex ( kind = 8 ) ima,zk0,zk1

    integer count1,count2
	real ( kind = 8 ) pi
	pi=3.1415926535897932384626433832795028841971d0

	ima=(0.0d0,1.0d0)
	zk0=zk

	allocate(dipvec(3,ns_over))
	allocate(gradpot(3,ns))
	allocate(pot(ns))

    do count1=1,ns_over
		dipvec(1,count1)=srcover(10,count1)*rho_in(count1)*wts(count1)
		dipvec(2,count1)=srcover(11,count1)*rho_in(count1)*wts(count1)
		dipvec(3,count1)=srcover(12,count1)*rho_in(count1)*wts(count1)
	enddo

    !Computing the full operator
	call hfmm3d_t_d_g(eps,zk0,ns_over,srcover(1:3,:),dipvec,ns,srcvals(1:3,:),pot,gradpot)
	
	do count1=1,ns
		PHI(count1)=(srcvals(10,count1)*gradpot(1,count1)+srcvals(11,count1)&
		 &*gradpot(2,count1)+srcvals(12,count1)*gradpot(3,count1))/(4.0d0*pi)
	enddo

	deallocate(dipvec)
	deallocate(pot)
	deallocate(gradpot)

return
end subroutine h_neumann_FMM_22

subroutine test_accuracy_h_neumann(eps_FMM,sol,sol2,zpars,ns,wts,srcvals,P0,Pt)
implicit none

!
!	 This function test the accuracy in the exterior region.
!

    !List of calling arguments
	complex ( kind = 8 ), intent(in) :: zpars(2)
    integer, intent(in) :: ns
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),eps_FMM
    real ( kind = 8 ), intent(in) :: wts(ns),P0(3),Pt(3)
	complex ( kind = 8 ), intent(in) :: sol(ns),sol2(ns)
	
    !List of local variables
	complex ( kind = 8 ) ima, pot_t1(1),pot_t2(1),pot_01(1), pot_02(1),gradpot_aux(3)
	complex (kind = 8 ) zk0,zk1,qv,zk,alpha
	real ( kind = 8 ) pi
	integer count1
		zk=zpars(1)
		alpha=zpars(2)
		ima=(0.0d0,1.0d0)
		pi=3.1415926535897932384626433832795028841971d0
		
		call h_neumann_FMM_targ(eps_FMM,zk,alpha,ns,srcvals,1,Pt,wts,sol(1:ns)&
		&,sol2(1:ns),pot_t1)
		
		qv=1.0d0
		call point_source_scalar_helmholtz2(P0,qv,1,Pt,zk,pot_t2,gradpot_aux,0)
		!write (*,*) pot_t1,pot_t2
		write (*,*) 'Relative Error in POT: ', abs(pot_t1(1)-pot_t2(1))/abs(pot_t2(1))

return
end subroutine test_accuracy_h_neumann






subroutine h_neumann_FMM_targ(eps,zk,alpha,ns,srcvals,nt,targ,wts,rho_in,rho_in2,pot)
implicit none
	!
	!	This funciton computes the potential u at a given point using the FMM
	!   It doesn't contain near field corrections (it's for debugging purposes)   
	!   The representation for the potentials is:
	!      Representation:
	!        u = S_{k}[\rho]+i*alpha*D_{k}[S_{ik}[\rho]]
	!
	! it requires as input \rho and \rho_2:=S_{ik}[\rho]

    !List of calling arguments
	real ( kind = 8 ), intent(in) :: eps
	complex ( kind = 8 ), intent(in) :: zk, alpha
    integer, intent(in) :: ns,nt
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),targ(3,nt)
    real ( kind = 8 ), intent(in) :: wts(ns)
    complex ( kind = 8 ), intent(in) :: rho_in(ns),rho_in2(ns)
    complex ( kind = 8 ), intent(out) :: pot(nt)


    !List of local variables
    complex ( kind = 8 ), allocatable :: charge(:),dipvec(:,:),pot_aux(:),rho_over(:)
	complex ( kind = 8 ) ima

    integer count1,count2,npols,ntri,npols_over
    integer ifa_vect,ifb_vect,iflambda,ifrho,ifE,ifcurlE,ifdivE
	real ( kind = 8 ) pi
	pi=3.1415926535897932384626433832795028841971d0


	ima=(0.0d0,1.0d0)

	allocate(charge(ns))
	allocate(pot_aux(ns))
	allocate(dipvec(3,ns))
	
	do count1=1,ns
		charge(count1)=rho_in(count1)*wts(count1)
		dipvec(1,count1)=srcvals(10,count1)*rho_in2(count1)*wts(count1)*ima*alpha
		dipvec(2,count1)=srcvals(11,count1)*rho_in2(count1)*wts(count1)*ima*alpha
		dipvec(3,count1)=srcvals(12,count1)*rho_in2(count1)*wts(count1)*ima*alpha
	enddo

    !Computing the full operator
    call hfmm3d_t_cd_p(eps,zk,ns,srcvals(1:3,:),charge,dipvec,nt,targ,pot)
	
	do count1=1,nt
		pot(count1)=pot(count1)/(4.0d0*pi)
	enddo

	deallocate(charge)
	deallocate(dipvec)

return
end subroutine h_neumann_FMM_targ


subroutine 	get_RHS_h_neumann(P0,ns,srcvals,zk,RHS)
implicit none

!
!   It obtains the right hand side of the problem, that is u'_inc
!

	!List of calling arguments
	integer, intent(in) :: ns
	real ( kind = 8 ), intent(in) :: P0(3)
	real ( kind = 8 ), intent(in) :: srcvals(12,ns)
	complex ( kind = 8 ), intent(in) :: zk
	complex ( kind = 8 ), intent(out) :: RHS(ns)
	
	!List of local variables
	complex ( kind = 8 ), allocatable :: gradpot(:,:),pot(:)
	integer count1
	complex ( kind = 8 ) ima,zk0,zk1,qv
	real ( kind = 8 ) ru(3),rv(3),cross_aux(3)
	
	allocate(gradpot(3,ns),pot(ns))

	ima=(0.0d0,1.0d0)
	
	qv=1.0d0
	call point_source_scalar_helmholtz2(P0,qv,ns,srcvals(1:3,:),zk,pot,gradpot,0)	
	
	do count1=1,ns
		RHS(count1)=DOT_PRODUCT(srcvals(10:12,count1),gradpot(:,count1))
	enddo
	
	deallocate(pot)
	deallocate(gradpot)

return
end subroutine get_RHS_h_neumann

