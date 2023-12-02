!                                                                       1
      subroutine getnearquad_em_sdpie_pec(npatches,norders,&
     &ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,&
     &iquad,rfac0,nquad,wnear)
	
!
!  This subroutine generates the near field quadrature
!  for the integral equation:
!
!    J/2+nxcurlS_{k}[J]+i·alpha·nxcurlcurlS_{k}[nxS_{ik}[J]]= -nxE_inc
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
!  input:
!    npatches - integer(8)
!      number of patches
!
!    norders - integer(8)(npatches)
!      order of discretization on each patch 
!
!    ixyzs - integer(8)(npatches+1)
!      starting location of data on patch i
!  
!    iptype - integer(8)(npatches)
!      type of patch
!      iptype = 1 -> triangular patch discretized with RV nodes
!
!    npts - integer(8)
!      total number of discretization points on the boundary
!
!      srccoefs - real *8 (9,npts)
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
!    ndtarg - integer(8)
!      leading dimension of target array
!        
!    ntarg - integer(8)
!      number of targets
!
!    targs - real *8 (ndtarg,ntarg)
!      target information
!
!    ipatch_id - integer(8)(ntarg)
!      id of patch of target i, id = -1, if target is off-surface
!
!    uvs_targ - real *8 (2,ntarg)
!      local uv coordinates on patch if on surface, otherwise
!      set to 0 by default
!      (maybe better to find closest uv on patch using
!      newton)
!            
!    eps - real *8
!      precision requested
!
!    zpars - complex *16 (3)
!      kernel parameters
!      zpars(1) = k 
!      zpars(2) = alpha
!      zpars(3) = - (not used)
!
!    iquadtype - integer(8)
!      quadrature type
!      iquadtype = 1, use ggq for self + adaptive integration
!      for rest
!
!    nnz - integer(8)
!      number of source patch-> target interactions in the near field
! 
!    row_ptr - integer(8)(ntarg+1)
!      row_ptr(i) is the pointer
!      to col_ind array where list of relevant source patches
!      for target i start
!
!    col_ind - integer(8) (nnz)
!      list of source patches relevant for all targets, sorted
!      by the target number
!
!    iquad - integer(8)(nnz+1)
!      location in wnear array where quadrature for col_ind(i)
!      starts
!
!    rfac0 - integer(8)
!      radius parameter for near field
!
!    nquad - integer(8)
!      number of entries in wnear
!
!  output
!    wnear  - complex *16(16*nquad)
!    the desired near field quadrature
!         
!

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
      complex *16 wnear(2*nquad)
      integer(8) ipatch_id(ntarg),count1,count2,icount
      real *8 uvs_targ(2,ntarg)

      complex *16 alpha,beta
      integer(8) i,j,ndi,ndd,ndz

      integer(8) ipv

!      procedure (), pointer :: fker
!      external h3d_slp, h3d_dlp, h3d_comb

      procedure (), pointer :: fker
      external  fker_em_sdpie_pec
      ndz=3
      ndd=1
      ndi=2
      ipv=1

!     write (*,*) 'Datos Entrada'
      fker =>  fker_em_sdpie_pec
!	 write(*,*) ndtarg,ntarg
!	 write (*,*) 'Fin Entrada'

!		 read (*,*)

	 write (*,*) 'generating near'
      icount=1
      ipars(1)=1
      ipars(2)=1
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
     &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
     &wnear((icount-1)*nquad+1:(icount)*nquad))
      icount=2
      ipars(1)=2
      ipars(2)=1
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
     &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
     &wnear((icount-1)*nquad+1:(icount)*nquad))
      return
      end subroutine getnearquad_em_sdpie_pec


subroutine em_sdpie_pec_solver(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,eps,zpars,numit,ifinout,&
     &rhs,eps_gmres,niter,errs,rres,soln,&
	 &ncomponents,idcomponent)

!
!
!  This subroutine solves the Scattering Maxwell p.e.c. problem.
!  The the equations are:
!
!    curlE=ikH; curlH =-ikE
!
!  This subroutine only solves E, there is an auxiliar scalar equation
!  to obtain H
!
!  Representation:
!
!    E=curlS_{k}[J]+i·alpha·curlcurlS_{k}[nxS_{ik}[J]]
!
!  Note: in order to make live simpler to evaluate E, the output will
!  be
!    J and both nxS_{ik}[J] and n·curlS_{ik}[J];
!    J -> soln and (nxS_{ik}[J],n·curlS_{ik}[J]) -> sigma2
!
!  Note: the representation for H is:
!    H=ik·S_{k}[J]+grad(PSI)+i·alpha·k·curlS_{k}[nxS_{ik}[J]]
!    therefore, an extra auxiliary scalar Neumann problem has to be
!    solved to obtain the magentic field H
!
!  Boundary conditions:
!
!    nxE = -nxE_inc
!
!  The incoming fields can be any vector function tangent to the surface 
!
!  Boundary integral equation:
!
!    J/2+nxcurlS_{k}[J]+i·alpha·nxcurlcurlS_{k}[nxS_{ik}[J]]= -nxE_inc
!
!  The linear system is solved iteratively using GMRES
!  until a relative residual of 1e-15 is reached
!
!
!  input:
!    npatches - integer(8)
!      number of patches
!
!    norders- integer(8)(npatches)
!      order of discretization on each patch 
!
!    ixyzs - integer(8)(npatches+1)
!      ixyzs(i) denotes the starting location in srccoefs,
!      and srcvals array corresponding to patch i
!
!    iptype - integer(8)(npatches)
!      type of patch
!      iptype = 1, triangular patch discretized using RV nodes
!
!    npts - integer(8)
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
!    zpars - complex *16 (3)
!      kernel parameters (Referring to formula (1))
!      zpars(1) = k 
!      zpars(2) = alpha
!      zpars(3) = - (not used)
!
!    ifinout - integer(8)
!      flag for interior or exterior problems (normals assumed to 
!      be pointing in exterior of region)
!        ifinout = 0, interior problem
!        ifinout = 1, exterior problem
!
!    rhs - complex *16(npts)
!      right hand side
!
!    eps_gmres - real *8
!      gmres tolerance requested
!
!    numit - integer(8)
!      max number of gmres iterations
!
!  output
!    niter - integer(8)
!      number of gmres iterations required for relative residual
!      to converge to 1e-15
!          
!    errs(1:iter) - relative residual as a function of iteration
!      number
! 
!    rres - real *8
!      relative residual for computed solution
!              
!    soln - complex *16(2*npts)
!      soln(1:npts) component of the tangent induced current J on the
!        surface along srcvals(4:6,i) direction
!      soln(npts+1:2*npts) component of the tangent induced current J 
!        on the surface along (srcvals(10:12,i) x srcvals(4:6,i))
!        direction
!
!    sigma2 - complex *16(3*npts)
!      sigma2(1:npts) component of nxS_{ik][J] on the surface
!        along srcvals(4:6,i) direction
!      sigma2(npts+1:2*npts) component of nxS_{ik][J]  on the surface
!        along (srcvals(10:12,i) x srcvals(4:6,i)) direction
!      sigma2(2*npts+1:3*npts) component of n·curlS_{ik][J] on the
!        surface along srcvals(10:12,i) direction
!
!    scalar_RHS - complex *16(npts)     !! (VECT2SCAL)
!      contains n·ik·S_{k}[J]+i·alpha·k·n·curlS_{k}[nxS_{ik}[J]]
!      which is used to solve the auxiliary scalar Neumann problem 
!      to obtain H
!

      implicit none
      integer(8) ipars,ncomponents
      integer(8) npatches,norder,npols,npts
      integer(8) ifinout
      integer(8) norders(npatches),ixyzs(npatches+1),idcomponent(npatches)
      integer(8) iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,eps_gmres
      complex *16 zpars(3)
      complex *16 rhs(npts+ncomponents)
      complex *16 soln(npts+ncomponents)

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
      real *8, allocatable :: srcover(:,:),wover(:),wts(:)
      integer(8), allocatable :: ixyzso(:),novers(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:) 

      integer(8) i,j,jpatch,jquadstart,jstart

      real *8 dpars(1),timeinfo(10),t1,t2,omp_get_wtime


      real *8 ttot,done,pi
      real *8 rfac,rfac0
      integer(8) iptype_avg,norder_avg
      integer(8) ikerorder, iquadtype,npts_over
	  integer(8) n_var

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
!   as we have one vector unknown J we need n_var=2*npts
!

      n_var=npts+ncomponents
	  !CAMBIAR
!      n_var=npts

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
      call prinf('novers=*',novers,100)

      allocate(srcover(12,npts_over),wover(npts_over),wts(npts))

      call oversample_geom(npatches,norders,ixyzs,iptype,npts,&
     &srccoefs,srcvals,novers,ixyzso,npts_over,srcover)

      call get_qwts(npatches,novers,ixyzso,iptype,npts_over,&
     &srcover,wover)

      call get_qwts(npatches,norders,ixyzs,iptype,npts,&
     &srcvals,wts)

!
!   compute near quadrature correction
!
      nquad = iquad(nnz+1)-1
      print *, "nquad=",nquad
      allocate(wnear(2*nquad))
      
!C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,2*nquad
		wnear(i)=0
      enddo
!C$OMP END PARALLEL DO    


      iquadtype = 1

!!      eps2 = 1.0d-8

      print *, "starting to generate near quadrature"
      call cpu_time(t1)
!C$      t1 = omp_get_wtime()      


      call getnearquad_em_sdpie_pec(npatches,norders,&
     &ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,&
     &ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,&
     &iquad,rfac0,nquad,wnear)
	 
!	 do count1=1,100
!		write(*,*) count1,wnear34(count1)
!	 enddo
!stop
!	 read (*,*)
	 
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
!      zid = -(-1)**(ifinout)*2*pi!*zpars(3)
	  zid=0.5d0
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
        call lpcomp_sdpie_addsub(npatches,norders,ixyzs,&
        &iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,&
        &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,&
        &vmat(1,it),novers,npts_over,ixyzso,srcover,wover,wtmp,&
        &wnear,ncomponents,idcomponent,wts)

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
	 
          call lpcomp_sdpie_addsub(npatches,norders,ixyzs,&
          &iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,&
          &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,&
          &soln,novers,npts_over,ixyzso,srcover,wover,wtmp,&
          &wnear,ncomponents,idcomponent,wts)

            
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
      end subroutine em_sdpie_pec_solver


subroutine fker_em_sdpie_pec(srcinfo, ndt,targinfo,ndd, dpars,ndz,&
 &zpars,ndi,ipars,E_val)
implicit none

! this function provides the near field kernel that will use 
! zgetnearquad_ggq_guru through getnearquad_CKi

    !List of calling arguments
    integer(8), intent(in) :: ndt,ndd,ndz,ndi
    real ( kind = 8 ), intent(in) :: srcinfo(12)
    real ( kind = 8 ), intent(in) :: targinfo(ndt)
    integer(8), intent(in) :: ipars(ndi)
    real ( kind = 8 ), intent(in) :: dpars(ndd)
    complex ( kind = 8 ), intent(in) :: zpars(ndz)
    complex ( kind = 8 ), intent(out) :: E_val
	
    !List of local variables
	complex ( kind = 8 ) zpars_aux(2)
	real ( kind = 8 ) ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ) ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ) du(3), dv(3),sour(3),targ(3)
	real ( kind = 8 ) r, dr(3),aux_real
	complex ( kind = 8 ) Dkrho(1,1),Sklambda(1,1),ngradDkrho(1,1),ngradDk0rho(1,1)
	complex ( kind = 8 ) ngradDk_diff(1,1),ngradSklambda(1,1)
	complex ( kind = 8 ) R1,R2, ima,my_exp, zk,alpha,E_mat(2,1),R1_0,R2_0,zk0,my_exp_0,D_diff
	real ( kind = 8 ) pi
	pi=3.1415926535897932384626433832795028841971d0
	ima=(0.0d0,1.0d0)
	zk=zpars(1)
	zk0=0.0d0
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
	
    call get_Dkrho(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r,Dkrho)
    call get_Sklambda(my_exp,r,Sklambda)
	
	call get_ngradSklambda(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r,ngradSklambda)
	
	zpars_aux(1)=zpars(1)
	zpars_aux(2)=0.0d0
	
!	R2=((ima*zk)**2/r**3-3.0d0*ima*zk/r**4+3.0d0/r**5)*exp(ima*zk*r)/&
!	 &(4.0d0*pi)
!	 call get_ngradDkrho(n_s,n_t,dr,R1,R2,zk,my_exp,r,ngradDkrho)
	 
!	R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
!	my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
!	R2_0=((ima*zk0)**2/r**3-3.0d0*ima*zk0/r**4+3.0d0/r**5)*exp(ima*zk0*r)/&
!	 &(4.0d0*pi)
	
!	call get_ngradDkrho(n_s,n_t,dr,R1_0,R2_0,zk0,my_exp_0,r,ngradDk0rho)
	
	call h3d_dprime_diff(srcinfo, ndt,targinfo, ndd,dpars,ndz,zpars_aux,ndi, &
  ipars,D_diff)
!  D_diff=D_diff
!write(*,*) 'diff1',ngradDkrho(1,1)-ngradDk0rho(1,1),D_diff
!read (*,*)
	E_mat(1,1)=Dkrho(1,1)-ima*alpha*Sklambda(1,1)
!	E_mat(2,1)=(ngradDkrho(1,1)-ngradDk0rho(1,1))-ima*alpha*ngradSklambda(1,1)
	E_mat(2,1)=(D_diff)-ima*alpha*ngradSklambda(1,1)
!write (*,*) ngradDkrho(1,1)-ngradDk0rho(1,1),D_diff
!read (*,*)
    E_val=E_mat(ipars(1),ipars(2))

return
end subroutine fker_em_sdpie_pec

      subroutine lpcomp_sdpie_addsub(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,sigma,novers,&
     &nptso,ixyzso,srcover,whtsover,pot,&
     &wnear,ncomponents,idcomponent,wts)
!
!  This subroutine evaluates the layer potential for
!  the boundary integral equation:
!
!    J/2+nxcurlS_{k}[J]+i·alpha·nxcurlcurlS_{k}[nxS_{ik}[J]]= -nxE_inc
!
!  where the near field is precomputed and stored
!  in the row sparse compressed format.
!
!  Note the 4\pi scaling is NOT included as the FMM output was scaled
!  appropriately
!
!  Note: the identity J/2 is not included as the gmres takes care of
!  that
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
!    npatches - integer(8)
!      number of patches
!
!    norders- integer(8)(npatches)
!      order of discretization on each patch 
!
!    ixyzs - integer(8)(npatches+1)
!      ixyzs(i) denotes the starting location in srccoefs,
!      and srcvals array corresponding to patch i
!   
!    iptype - integer(8)(npatches)
!      type of patch
!      iptype = 1, triangular patch discretized using RV nodes
!
!    npts - integer(8)
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
!    ndtarg - integer(8)
!      leading dimension of target array
!        
!    ntarg - integer(8)
!      number of targets
!
!    targs - real *8 (ndtarg,ntarg)
!      target information
!
!    eps - real *8
!      precision requested
!
!    zpars - complex *16 (3)
!      kernel parameters (Referring to formula (1))
!      zpars(1) = k 
!      zpars(2) = alpha
!      zpars(3) = - ( not used )
!
!    nnz - integer(8) *8
!      number of source patch-> target interactions in the near field
! 
!    row_ptr - integer(8)(ntarg+1)
!      row_ptr(i) is the pointer
!      to col_ind array where list of relevant source patches
!      for target i start
!
!    col_ind - integer(8) (nnz)
!      list of source patches relevant for all targets, sorted
!      by the target number
!
!    iquad - integer(8)(nnz+1)
!      location in wnear array where quadrature for col_ind(i)
!      starts
!
!    nquad - integer(8)
!      number of entries in wnear
!
!    wnear  - complex *16(12*nquad)
!      the desired near field quadrature
!
!    sigma - complex *16(2*ns)
!      induced charge and current on the surface
!      sigma(1:ns) - first component of  J along
!      the srcvals(4:6,i) direction
!      sigma(ns+1:2*ns) - second component of J along
!      the (srcvals(10:12,i) x srcvals(4:6,i)) direction
!
!    novers - integer(8)(npatches)
!      order of discretization for oversampled sources and
!      density
!
!    ixyzso - integer(8)(npatches+1)
!      ixyzso(i) denotes the starting location in srcover,
!      corresponding to patch i
!   
!    nptso - integer(8)
!      total number of oversampled points
!
!    srcover - real *8 (12,nptso)
!      oversampled set of source information
!
!    whtsover - real *8 (nptso)
!      smooth quadrature weights at oversampled nodes
!

      implicit none
      integer(8) npatches,norder,npols,npts
      integer(8) ndtarg,ntarg
      integer(8) norders(npatches),ixyzs(npatches+1)
      integer(8) ixyzso(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps
      real *8 targs(ndtarg,ntarg),wts(npts)
      complex *16 zpars(3),alpha,zk
      integer(8) nnz,row_ptr(ntarg+1),col_ind(nnz),nquad
      integer(8) iquad(nnz+1),ncomponents,idcomponent(npatches)
      complex *16 sigma(npts+ncomponents)
!      complex *16 sigma(npts)
	  !CAMBIAR
	  
      complex *16 wnear(2*nquad)

      integer(8) novers(npatches+1)
      integer(8) nover,npolso,nptso
      real *8 srcover(12,nptso),whtsover(nptso)
      complex *16 pot(ntarg+ncomponents)
	  !CAMBIAR
!	  complex *16 pot(ntarg)
	  
      complex *16, allocatable :: potsort(:)

      real *8, allocatable :: sources(:,:),targvals(:,:),wtmp2(:)
      complex *16, allocatable :: charges(:),dipvec(:,:),sigmaover(:)
	  complex *16, allocatable :: sigma_aux(:),sigmaover_aux(:)
	  complex *16, allocatable :: pot_aux(:)

      integer(8) ns,nt
      complex *16 beta
      integer(8) ifcharge,ifdipole
      integer(8) ifpgh,ifpghtarg
      complex *16 tmp(10),val,E(2),zzzt

      real *8 xmin,xmax,ymin,ymax,zmin,zmax,sizey,sizez,boxsize

      integer(8) i,j,jpatch,jquadstart,jstart

      integer(8) ifaddsub,ifdir

      integer(8) ntj
      
      complex *16 zdotu,pottmp
      complex *16, allocatable :: dtmp2(:,:),ctmp2_b_u(:),ctmp2_b_v(:)
	  complex *16, allocatable :: ctmp2_u(:),ctmp2_v(:),ctmp2_s(:),ctmp2_b_s(:)
      real *8 radexp,epsfmm

      integer(8) ipars(2)
      real *8 dpars(1),timeinfo(10),t1,t2,omp_get_wtime

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      real *8 thresh,ra
      real *8 rr,rmin
      integer(8) nss,ii,l,npover

      integer(8) nd,ntarg0,count1,count2,icount,naux

      real *8 ttot,done,pi

      complex *16 sigma_x,sigma_y,sigma_z
      complex *16 pot_x,pot_y,pot_z,pot2_x,pot2_y,pot2_z
      real *8 u_vect(3),v_vect(3),n_vect(3)
	  complex *16 ima

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
	    allocate(sigma_aux(ntarg+ncomponents),pot_aux(2*ntarg))

! 
!       oversample density
!
		  zk=zpars(1)
		  alpha=zpars(2)

      call oversample_fun_surf(int(2,8),npatches,norders,ixyzs,iptype,& 
     &npts,sigma(1:npts),novers,ixyzso,ns,sigmaover(1:ns))

      ra = 0
      
      !
      !     FMM call to compute nxSkia and -n·curlSkia
      !
	   call get_thresh(srcover,ns,targs,ntarg,thresh)
	  
       ifdir=0
	   call em_sdpie_pec_FMM(eps,zpars,ns,npts,srcover,targs,whtsover,sigmaover(1:ns),&
	   &pot_aux(1:npts),pot_aux(npts+1:2*npts),thresh,ifdir)	  
!
!        compute threshold for ignoring local computation
!



!
!       Add near field precomputed contribution
!
      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch)
          do l=1,npols
            pot_aux(i)=pot_aux(i)+wnear(jquadstart+l-1)*sigma(jstart+l-1)
!			write (*,*) pot_aux(i),pot_aux(i+npts),wnear(nquad+jquadstart+l-1),sigma(jstart+l-1)
!			read (*,*)
            pot_aux(i+npts)=pot_aux(i+npts)+wnear(nquad+jquadstart+l-1)*sigma(jstart+l-1)
          enddo
        enddo
      enddo
!
!     Remove near contribution of the FMM
!

      ifdir=1
	  ipars(1)=1
	  ipars(2)=2
	  do i=1,ntarg
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          nss = nss + ixyzso(jpatch+1)-ixyzso(jpatch)
        enddo
        allocate(srctmp2(12,nss),ctmp2_u(nss))
        allocate(wtmp2(nss))

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
		
		call em_sdpie_pec_FMM(eps,zpars,nss,ntarg0,srctmp2,targs(:,i),wtmp2,ctmp2_u,E(1),E(2),thresh,ifdir)	  
!
!        call em_sdpie_pec_FMM(eps,zpars(1),nss,ntarg0,srctmp2,&
!        &targs(:,i),wtmp2,ctmp2_u,E(1),E(2),thresh,&
!        &ifdir)
	    pot_aux(i) = pot_aux(i) - E(1)
	    pot_aux(i+ntarg) = pot_aux(i+ntarg) - E(2)

        deallocate(srctmp2,ctmp2_u,wtmp2)
      enddo

       do count1=1,npts
    	   pot(count1)=pot_aux(count1)
	   enddo
	   !CAMBIAR
       do count1=1,ncomponents
	       pot(npts+count1)=0.0d0
	   enddo
		icount=1
!        zzzt=0.0d0
		do count1=1,npatches
		  naux=(norders(count1)+1)*(norders(count1)+2)/2
		  do count2=1,naux
			pot(ntarg+idcomponent(count1))=pot(ntarg+idcomponent(count1))+&
			&(pot_aux(ntarg+icount)+ima*alpha*.5d0*sigma(icount))*wts(icount)
			
			
!			zzzt=zzzt+&
!			&(pot_aux(ntarg+icount)+ima*alpha*.5d0*sigma(icount))*wts(icount)
!           write (*,*) zzzt,pot_aux(ntarg+icount),sigma(icount),wts(icount),ima,alpha
!		   read (*,*)
			pot(icount)=pot(icount)+sigma(npts+idcomponent(count1))
			icount=icount+1
		  enddo
		enddo
       do count1=1,ncomponents
!	   write (*,*) 'pot_value: ',pot(npts+count1),sigma(npts+count1)
	       pot(npts+count1)=pot(npts+count1)-.5d0*sigma(npts+count1)
	   enddo
!write (*,*) 'valor del flujo: ',zzzt
      return
      end subroutine lpcomp_sdpie_addsub



subroutine em_sdpie_pec_FMM(eps,zpars,ns,nt,srcvals,targvals,wts,rho_in,PHI,nGRAd,thresh,ifdir)
implicit none

!
!	One of the FMM operators needed in the formulation
!

    !List of calling arguments
	  real ( kind = 8 ), intent(in) :: eps
	  complex ( kind = 8 ), intent(in) :: zpars(3)
    integer(8), intent(in) :: ns,nt
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),wts(ns)
    real ( kind = 8 ), intent(in) :: targvals(12,nt)
    complex ( kind = 8 ), intent(in) :: rho_in(ns)
    complex ( kind = 8 ), intent(out) :: PHI(nt),nGRAD(nt)
    real ( kind = 8 ), intent(in) :: thresh
    integer(8), intent(in) :: ifdir 


    !List of local variables
    complex ( kind = 8 ), allocatable :: dipvec(:,:),rho(:)
    complex ( kind = 8 ), allocatable :: pot(:),gradpot(:,:)
    real ( kind = 8 ), allocatable :: source(:,:),targets(:,:)

	  complex ( kind = 8 ) ima,zk

    integer(8) count1,count2,nd,ier
	  real ( kind = 8 ) pi
	  pi=3.1415926535897932384626433832795028841971d0

	  ima=(0.0d0,1.0d0)
	  zk=zpars(1)
    nd=1
	  allocate(rho(ns))
	  allocate(dipvec(3,ns))
	  allocate(gradpot(3,nt))
	  allocate(pot(nt))
    allocate(source(3,ns))
    allocate(targets(3,nt))

    do count1=1,ns
      source(:,count1)=srcvals(1:3,count1)
    enddo
	
    do count1=1,nt
      targets(:,count1)=targvals(1:3,count1)
    enddo

    do count1=1,ns
		  rho(count1)=-ima*zpars(2)*rho_in(count1)*wts(count1)
		  dipvec(1,count1)=srcvals(10,count1)*rho_in(count1)*wts(count1)
		  dipvec(2,count1)=srcvals(11,count1)*rho_in(count1)*wts(count1)
		  dipvec(3,count1)=srcvals(12,count1)*rho_in(count1)*wts(count1)
	  enddo

    !Computing the full operator
		if (ifdir.eq.1) then
      do count1=1,nt
		    PHI(count1)=0.0d0
		    gradpot(1,count1)=0.0d0
		    gradpot(2,count1)=0.0d0
		    gradpot(3,count1)=0.0d0
	    enddo
      call h3ddirectcdg(nd,zk,source,rho,dipvec,ns&
	     &,targets,nt,PHI,gradpot,thresh)
    else
      call hfmm3d_t_cd_g_vec(nd,eps,zk,ns,source,rho,dipvec&
	     &,nt,targets,PHI,gradpot,ier)
    endif

	  do count1=1,nt
      PHI(count1)=PHI(count1)/(4.0d0*pi)
      nGRAD(count1)=(targvals(10,count1)*gradpot(1,count1)+targvals(11,count1)&
      &*gradpot(2,count1)+targvals(12,count1)*gradpot(3,count1))/(4.0d0*pi)
	  enddo
	
	  zk=1.0d-15
    if (ifdir.eq.1) then
      do count1=1,nt
		  pot(count1)=0.0d0
		  gradpot(1,count1)=0.0d0
		  gradpot(2,count1)=0.0d0
		  gradpot(3,count1)=0.0d0
    enddo
	  nd=2
	  call l3ddirectdg(nd,source,dipvec,ns,targets,nt,pot,&
	   &gradpot,thresh)
    else
	  nd=2
      call lfmm3d_t_d_g_vec(nd,eps,ns,source,dipvec,nt,targets&
	   &,pot,gradpot)
    endif

	  do count1=1,nt
		 nGRAD(count1)=nGRAD(count1)-(targvals(10,count1)*gradpot(1,count1)+targvals(11,count1)&
		 &*gradpot(2,count1)+targvals(12,count1)*gradpot(3,count1))/(4.0d0*pi)
	  enddo

	deallocate(rho)
	deallocate(dipvec)
	deallocate(pot)
	deallocate(gradpot)
	deallocate(source)
	deallocate(targets)

return
end subroutine em_sdpie_pec_FMM





subroutine 	get_rhs_sdpie_pec(P0,vf,ns,srcvals,zk,RHS,&
&npatches,norders,ncomponents,idcomponent,ixyzs,iptype)
implicit none

!
!  This function obtains the right hand side for the CKi formulation
!  for the integral boundary equation:
!
!    J/2+nxcurlS_{k}[J]+i·alpha·nxcurlcurlS_{k}[nxS_{ik}[J]]= -nxE_inc
!
!  input:
!    P0 - real * 8 (3)
!      location of the source point at the exterior region
!      WARNING! notice that this formulation uses a representation 
!      theorem for the incoming field in the interior region (MFIE)
!      therefore it only works for incoming fields generated by sources
!      in the exterior region (or at infinity like plane waves)
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
!  output:
!    RHS - complex  *16(3*ns)
!      right hand side
!      RHS(1:ns) - first component of -nxE_inc along
!      the srcvals(4:6,i) direction
!      RHS(ns+1:2*ns) - second component of  nxE_inc along
!      the (srcvals(10:12,i) x srcvals(4:6,i)) direction
!


	!List of calling arguments
	integer(8), intent(in) :: ns
	real ( kind = 8 ), intent(in) :: P0(3)
	complex ( kind = 8 ), intent(in) :: vf(3)
	real ( kind = 8 ), intent(in) :: srcvals(12,ns)
	complex ( kind = 8 ), intent(in) :: zk
	integer(8), intent(in) :: npatches,ncomponents
	integer(8), intent(in) :: idcomponent(npatches),norders(npatches)
	integer(8), intent(in) :: ixyzs(npatches+1),iptype(npatches)
	complex ( kind = 8 ), intent(out) :: RHS(ns+ncomponents)
	
	!List of local variables
	complex ( kind = 8 ), allocatable :: E(:,:), H(:,:),dpot_dnormal(:),pot(:)
	integer(8) count1,count2,icount,naux
	real ( kind = 8 ) ru(3),rv(3),cross_aux(3)
	real ( kind = 8 ), allocatable :: wts(:)

	allocate(E(3,ns), H(3,ns))
	allocate(pot(ns),dpot_dnormal(ns))
	allocate(wts(ns))
!!	call fieldsMD(zk,P0,srcvals(1:3,:),ns,curlE,E,vf,0)
!	write (*,*) 'P0',P0,ns
	call point_source_scalar_helmholtz(P0,ns,srcvals(1:3,:),srcvals(10:12,:),zk,pot,dpot_dnormal)
!	read (*,*)
	do count1=1,ns	
		RHS(count1)=pot(count1)		
	enddo
    call get_qwts(npatches,norders,ixyzs,iptype,ns,srcvals,wts)
    do count1=1,ncomponents
		RHS(ns+count1)=0.0d0
	enddo

	icount=1
    do count1=1,npatches
      naux=(norders(count1)+1)*(norders(count1)+2)/2
	  do count2=1,naux
	    RHS(ns+idcomponent(count1))=RHS(ns+idcomponent(count1))+dpot_dnormal(icount)*wts(icount)
	    icount=icount+1
	  enddo
    enddo
return 
end subroutine get_rhs_sdpie_pec




subroutine test_accuracy_sdpie(eps_FMM,sol,zk,alpha,ns,wts,srcvals,P0,vf,Pt)
implicit none

    !List of calling arguments
	complex ( kind = 8 ), intent(in) :: zk,alpha
    integer(8), intent(in) :: ns
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),eps_FMM
    real ( kind = 8 ), intent(in) :: wts(ns),P0(3),Pt(3)
	complex ( kind = 8 ), intent(in) :: sol(3*ns),vf(3)
	
    !List of local variables
	complex ( kind = 8 ) a_u00,a_v00
	complex ( kind = 8 ) ima, R1, R2,Et1(3),Et2(3),aux_cmp,Ht2(3),phi1(1),E_aux(3),phi2
	real ( kind = 8 ) ru_s(3),rv_s(3),n_s(3),sour(3),r,dr(3)
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3),error,rel_err
	real ( kind = 8 ) pi

	integer(8) count1
		ima=(0.0d0,1.0d0)
		pi=3.1415926535897932384626433832795028841971d0

call sdpie_FMM_targ(eps_FMM,zk,alpha,ns,srcvals,1,Pt,wts,sol(1:ns),phi1)
		
		
!		call fieldsED(zk,P0,Pt,1,Et2,Ht2,vf,0)
!		call fieldsMD(zk,P0,Pt,1,Et2,Ht2,vf,1)
		!call point_source_vector_helmholtz(1,P0,vf,Pt,zk,Et2,E_aux,phi2)
		!!	call point_source_scalar_helmholtz(P0,ns,points,normals,zk,pot,dpot_dnormal)
!		call point_source_scalar_helmholtz(P0,ns,srcvals(1:3,:),srcvals(10:12,:),zk,pot,dpot_dnormal)
        call point_source_scalar_helmholtz(P0,1,Pt,Pt,zk,phi2,Et2(1))

!		write (*,*) Et2

!		error=sqrt(abs(Et1(1)-Et2(1))**2+abs(Et1(2)-Et2(2))**2+abs(Et1(3)-Et2(3))**2)
!		write (*,*) 'Error: ', error
!		write (*,*) 'Relative Error: ', error/sqrt(abs(Et2(1))**2+abs(Et2(2))**2+abs(Et2(3))**2)
		write (*,*) 'relative error in phi: ',abs(phi1-phi2)/abs(phi2)
!		write (*,*) 'phi1: ', phi1, 'phi2: ', phi2

return
end subroutine test_accuracy_sdpie


subroutine sdpie_FMM_targ(eps,zk,alpha,ns,srcvals,nt,targ,wts,rho_in,pot)
implicit none

  !List of calling arguments
	real ( kind = 8 ), intent(in) :: eps
	complex ( kind = 8 ), intent(in) :: zk,alpha
  integer(8), intent(in) :: ns,nt
  real ( kind = 8 ), intent(in) :: srcvals(12,ns),targ(3,nt)
  real ( kind = 8 ), intent(in) :: wts(ns)
  complex ( kind = 8 ), intent(in) :: rho_in(ns)
  complex ( kind = 8 ), intent(out) :: pot(nt)

  !List of local variables
	real ( kind = 8 ), allocatable :: n_vect(:,:),u_vect(:,:),v_vect(:,:),source(:,:)
  complex ( kind = 8 ), allocatable :: dipvec(:,:),rho(:)
	complex ( kind = 8 ) ima

  integer(8) count1,count2,ier
	real ( kind = 8 ) pi
	pi=3.1415926535897932384626433832795028841971d0

	ima=(0.0d0,1.0d0)

  allocate(dipvec(3,ns))
	allocate(rho(ns))
	
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
		rho(count1)=-ima*alpha*rho_in(count1)*wts(count1)
		dipvec(1,count1)=n_vect(1,count1)*rho_in(count1)*wts(count1)
		dipvec(2,count1)=n_vect(2,count1)*rho_in(count1)*wts(count1)
		dipvec(3,count1)=n_vect(3,count1)*rho_in(count1)*wts(count1)
	enddo

  !Computing the full operator
	
  call hfmm3d_t_cd_p_vec(int(1,8),eps,zk,ns,source,rho,dipvec&
	 &,nt,targ,pot,ier)
	do count1=1,nt
	  pot(count1)=pot(count1)/(4.0d0*pi)
  enddo
	deallocate(dipvec)
	deallocate(rho)
	
	deallocate(u_vect)
	deallocate(v_vect)
	deallocate(n_vect)
	deallocate(source)

return
end subroutine sdpie_FMM_targ



subroutine sdpie_grad_FMM_targ(eps,zk,alpha,ns,srcvals,nt,targ,wts,rho_in,pot,gradpot)
implicit none

  !List of calling arguments
	real ( kind = 8 ), intent(in) :: eps
	complex ( kind = 8 ), intent(in) :: zk,alpha
  integer(8), intent(in) :: ns,nt
  real ( kind = 8 ), intent(in) :: srcvals(12,ns),targ(3,nt)
  real ( kind = 8 ), intent(in) :: wts(ns)
  complex ( kind = 8 ), intent(in) :: rho_in(ns)
  complex ( kind = 8 ), intent(out) :: pot(nt),gradpot(3,nt)

  !List of local variables
	real ( kind = 8 ), allocatable :: n_vect(:,:),u_vect(:,:),v_vect(:,:),source(:,:)
  complex ( kind = 8 ), allocatable :: dipvec(:,:),rho(:)
	complex ( kind = 8 ) ima

  integer(8) count1,count2,ier
	real ( kind = 8 ) pi
	pi=3.1415926535897932384626433832795028841971d0

	ima=(0.0d0,1.0d0)

  allocate(dipvec(3,ns))
	allocate(rho(ns))
	
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
		rho(count1)=-ima*alpha*rho_in(count1)*wts(count1)
		dipvec(1,count1)=n_vect(1,count1)*rho_in(count1)*wts(count1)
		dipvec(2,count1)=n_vect(2,count1)*rho_in(count1)*wts(count1)
		dipvec(3,count1)=n_vect(3,count1)*rho_in(count1)*wts(count1)
	enddo

  !Computing the full operator
	
  call hfmm3d_t_cd_g_vec(int(1,8),eps,zk,ns,source,rho,dipvec&
	  &,nt,targ,pot,gradpot,ier)

	do count1=1,nt
	  pot(count1)=pot(count1)/(4.0d0*pi)
    gradpot(1,count1)=gradpot(1,count1)/(4.0d0*pi)
    gradpot(2,count1)=gradpot(2,count1)/(4.0d0*pi)
    gradpot(3,count1)=gradpot(3,count1)/(4.0d0*pi)
  enddo

	deallocate(dipvec)
	deallocate(rho)
	
	deallocate(u_vect)
	deallocate(v_vect)
	deallocate(n_vect)
	deallocate(source)

return
end subroutine sdpie_grad_FMM_targ
