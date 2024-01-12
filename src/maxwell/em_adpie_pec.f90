!                                                                       1
      subroutine getnearquad_em_adpie_pec(npatches,norders,&
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
!    npatches - integer *8
!      number of patches
!
!    norders - integer *8(npatches)
!      order of discretization on each patch 
!
!    ixyzs - integer *8(npatches+1)
!      starting location of data on patch i
!  
!    iptype - integer *8(npatches)
!      type of patch
!      iptype = 1 -> triangular patch discretized with RV nodes
!
!    npts - integer *8
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
!    ndtarg - integer *8
!      leading dimension of target array
!        
!    ntarg - integer *8
!      number of targets
!
!    targs - real *8 (ndtarg,ntarg)
!      target information
!
!    ipatch_id - integer *8(ntarg)
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
!    iquadtype - integer *8
!      quadrature type
!      iquadtype = 1, use ggq for self + adaptive integration
!      for rest
!
!    nnz - integer *8
!      number of source patch-> target interactions in the near field
! 
!    row_ptr - integer *8(ntarg+1)
!      row_ptr(i) is the pointer
!      to col_ind array where list of relevant source patches
!      for target i start
!
!    col_ind - integer *8 (nnz)
!      list of source patches relevant for all targets, sorted
!      by the target number
!
!    iquad - integer *8(nnz+1)
!      location in wnear array where quadrature for col_ind(i)
!      starts
!
!    rfac0 - integer *8
!      radius parameter for near field
!
!    nquad - integer *8
!      number of entries in wnear
!
!  output
!    wnear  - complex *16(16*nquad)
!    the desired near field quadrature
!         
!

      implicit none 
      integer *8 npatches,norders(npatches),npts,nquad
      integer *8 ixyzs(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,rfac0
      integer *8 ndtarg,ntarg
      integer *8 iquadtype
      real *8 targs(ndtarg,ntarg)
      complex *16 zpars(3)
      integer *8 nnz,ipars(2)
      real *8 dpars(1)
      integer *8 row_ptr(ntarg+1),col_ind(nnz),iquad(nnz+1)
      complex *16 wnear(12*nquad)
      integer *8 ipatch_id(ntarg),count1,count2,icount
      real *8 uvs_targ(2,ntarg)

      complex *16 alpha,beta
      integer *8 i,j,ndi,ndd,ndz

      integer *8 ipv

!      procedure (), pointer :: fker
!      external h3d_slp, h3d_dlp, h3d_comb

      procedure (), pointer :: fker
      external  fker_em_adpie_pec
      ndz=3
      ndd=1
      ndi=2
      ipv=1

!     write (*,*) 'Datos Entrada'
      fker =>  fker_em_adpie_pec
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
      ipars(1)=1
      ipars(2)=2
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
     &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
     &wnear((icount-1)*nquad+1:(icount)*nquad))
	 
      icount=3     
      ipars(1)=1
      ipars(2)=3
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
     &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
     &wnear((icount-1)*nquad+1:(icount)*nquad))
	 
      icount=4
      ipars(1)=2
      ipars(2)=1
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
     &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
     &wnear((icount-1)*nquad+1:(icount)*nquad))
	 
	 
      icount=5	 
      ipars(1)=2
      ipars(2)=2
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
     &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
     &wnear((icount-1)*nquad+1:(icount)*nquad))
	 
      icount=6
      ipars(1)=2
      ipars(2)=3
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
     &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
     &wnear((icount-1)*nquad+1:(icount)*nquad))
	 
      icount=7
      ipars(1)=3
      ipars(2)=1
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
     &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
     &wnear((icount-1)*nquad+1:(icount)*nquad))

      icount=8
      ipars(1)=3
      ipars(2)=2
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,& 
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
     &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
     &wnear((icount-1)*nquad+1:(icount)*nquad))

      icount=9
      ipars(1)=3
      ipars(2)=3
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,& 
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
     &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
     &wnear((icount-1)*nquad+1:(icount)*nquad))

      icount=10
      ipars(1)=4
      ipars(2)=1
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,& 
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
     &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
     &wnear((icount-1)*nquad+1:(icount)*nquad))

      icount=11
      ipars(1)=4
      ipars(2)=2
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,& 
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
     &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
     &wnear((icount-1)*nquad+1:(icount)*nquad))

      icount=12
      ipars(1)=4
      ipars(2)=3
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,& 
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
     &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
     &wnear((icount-1)*nquad+1:(icount)*nquad))
      return
      end subroutine getnearquad_em_adpie_pec


subroutine em_adpie_pec_solver(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,eps,zpars,numit,ifinout,&
     &rhs,eps_gmres,niter,errs,rres,soln,sigma2,&
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
!    npatches - integer *8
!      number of patches
!
!    norders- integer *8(npatches)
!      order of discretization on each patch 
!
!    ixyzs - integer *8(npatches+1)
!      ixyzs(i) denotes the starting location in srccoefs,
!      and srcvals array corresponding to patch i
!
!    iptype - integer *8(npatches)
!      type of patch
!      iptype = 1, triangular patch discretized using RV nodes
!
!    npts - integer *8
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
!    ifinout - integer *8
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
!    numit - integer *8
!      max number of gmres iterations
!
!  output
!    niter - integer *8
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
      integer *8 ipars,ncomponents
      integer *8 npatches,norder,npols,npts
      integer *8 ifinout
      integer *8 norders(npatches),ixyzs(npatches+1),idcomponent(npatches)
      integer *8 iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,eps_gmres
      complex *16 zpars(3)
      complex *16 rhs(3*npts+ncomponents)
      complex *16 soln(3*npts+ncomponents),sigma2(3*npts+ncomponents)

      real *8, allocatable :: targs(:,:)
      integer *8, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      integer *8 ndtarg,ntarg

      real *8 errs(numit+1)
      real *8 rres,eps2
      integer *8 niter


      integer *8 nover,npolso,nptso
      integer *8 nnz,nquad
      integer *8, allocatable :: row_ptr(:),col_ind(:),iquad(:)

      complex *16, allocatable :: wnear(:)
      real *8, allocatable :: srcover(:,:),wover(:),wts(:)
      integer *8, allocatable :: ixyzso(:),novers(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:) 

      integer *8 i,j,jpatch,jquadstart,jstart

      real *8 dpars(1),timeinfo(10),t1,t2,omp_get_wtime


      real *8 ttot,done,pi
      real *8 rfac,rfac0
      integer *8 iptype_avg,norder_avg
      integer *8 ikerorder, iquadtype,npts_over
	  integer *8 n_var

!
!       gmres variables
!

      complex *16 zid,ztmp
      real *8 rb,wnrm2
      integer *8 numit,it,iind,it1,k,l,count1
      real *8 rmyerr
      complex *16 temp
      complex *16, allocatable :: vmat(:,:),hmat(:,:)
      complex *16, allocatable :: cs(:),sn(:)
      complex *16, allocatable :: svec(:),yvec(:),wtmp(:)

!
!   n_var is the number of unknowns in the linear system.
!   as we have one vector unknown J we need n_var=2*npts
!

      n_var=3*npts+ncomponents

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
      allocate(wnear(12*nquad))
      
!C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,12*nquad
		wnear(i)=0
      enddo
!C$OMP END PARALLEL DO    


      iquadtype = 1

!!      eps2 = 1.0d-8

      print *, "starting to generate near quadrature"
      call cpu_time(t1)
!C$      t1 = omp_get_wtime()      


      call getnearquad_em_adpie_pec(npatches,norders,&
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
        call lpcomp_adpie_addsub(npatches,norders,ixyzs,&
        &iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,&
        &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,&
        &vmat(1,it),novers,npts_over,ixyzso,srcover,wover,wtmp,sigma2,&
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
	 
          call lpcomp_adpie_addsub(npatches,norders,ixyzs,&
          &iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,&
          &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,&
          &soln,novers,npts_over,ixyzso,srcover,wover,wtmp,sigma2,&
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
      end subroutine em_adpie_pec_solver


subroutine fker_em_adpie_pec(srcinfo, ndt,targinfo,ndd, dpars,ndz,&
 &zpars,ndi,ipars,E_val)
implicit none

! this function provides the near field kernel that will use 
! zgetnearquad_ggq_guru through getnearquad_CKi

    !List of calling arguments
    integer *8, intent(in) :: ndt,ndd,ndz,ndi
    real ( kind = 8 ), intent(in) :: srcinfo(12)
    real ( kind = 8 ), intent(in) :: targinfo(ndt)
    integer *8, intent(in) :: ipars(ndi)
    real ( kind = 8 ), intent(in) :: dpars(ndd)
    complex ( kind = 8 ), intent(in) :: zpars(ndz)
    complex ( kind = 8 ), intent(out) :: E_val
	
    !List of local variables
	real ( kind = 8 ) ru_s(3),rv_s(3),n_s(3)
	real ( kind = 8 ) ru_t(3),rv_t(3),n_t(3)
	real ( kind = 8 ) du(3), dv(3),sour(3),targ(3)
	real ( kind = 8 ) r, dr(3),aux_real
	complex ( kind = 8 ) nxcurlSka(2,2),Dkrho,nxSknrho(2,1),nxSknxb(2,2),nxgradSklambda(2,1)
	complex ( kind = 8 ) Sklambda,divSknxb(1,2),ncurlSka(1,2),nSknxb(1,2),nSknrho(1,1),ngradSklambda(1,1)
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3),xprod_aux3(3),xprod_aux4(3)
	complex ( kind = 8 ) R1, ima,my_exp, zk,alpha,E_mat(4,3)
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

	call get_nxcurlSka(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r,nxcurlSka)
	call get_nxSknxb(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,my_exp,r,nxSknxb)
	E_mat(1,1)=nxcurlSka(1,1)+ima*alpha*nxSknxb(1,1)
	E_mat(1,2)=nxcurlSka(1,2)+ima*alpha*nxSknxb(1,2)
	E_mat(2,1)=nxcurlSka(2,1)+ima*alpha*nxSknxb(2,1)
	E_mat(2,2)=nxcurlSka(2,2)+ima*alpha*nxSknxb(2,2)
	call get_nxSknrho(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,my_exp,r,nxSknrho)
	call get_nxgradSklambda(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r,nxgradSklambda)
	E_mat(1,3)=0.0d0-nxSknrho(1,1)+ima*alpha*nxgradSklambda(1,1)
	E_mat(2,3)=0.0d0-nxSknrho(2,1)+ima*alpha*nxgradSklambda(2,1)
	call get_Dkrho(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r,Dkrho)
	call get_Sklambda(my_exp,r,Sklambda)
	E_mat(3,3)=Dkrho+ima*alpha*Sklambda*(-zk**2)
	call get_divSknxb(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r,divSknxb)
	E_mat(3,1)=0.0d0+ima*alpha*divSknxb(1,1)
	E_mat(3,2)=0.0d0+ima*alpha*divSknxb(1,2)	

	call get_ncurlSka(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r,ncurlSka)	
	call get_nSknxb(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,my_exp,r,nSknxb)
	E_mat(4,1)=ncurlSka(1,1)+ima*alpha*nSknxb(1,1)
	E_mat(4,2)=ncurlSka(1,2)+ima*alpha*nSknxb(1,2)


    call get_nSknrho(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,my_exp,r,nSknrho)
    call get_ngradSklambda(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1,my_exp,r,ngradSklambda)
    E_mat(4,3)=-nSknrho(1,1)+ima*alpha*ngradSklambda(1,1)

    E_val=E_mat(ipars(1),ipars(2))

return
end subroutine fker_em_adpie_pec

      subroutine lpcomp_adpie_addsub(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,sigma,novers,&
     &nptso,ixyzso,srcover,whtsover,pot,sigma2,&
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
!    npatches - integer *8
!      number of patches
!
!    norders- integer *8(npatches)
!      order of discretization on each patch 
!
!    ixyzs - integer *8(npatches+1)
!      ixyzs(i) denotes the starting location in srccoefs,
!      and srcvals array corresponding to patch i
!   
!    iptype - integer *8(npatches)
!      type of patch
!      iptype = 1, triangular patch discretized using RV nodes
!
!    npts - integer *8
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
!    ndtarg - integer *8
!      leading dimension of target array
!        
!    ntarg - integer *8
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
!    nnz - integer *8 *8
!      number of source patch-> target interactions in the near field
! 
!    row_ptr - integer *8(ntarg+1)
!      row_ptr(i) is the pointer
!      to col_ind array where list of relevant source patches
!      for target i start
!
!    col_ind - integer *8 (nnz)
!      list of source patches relevant for all targets, sorted
!      by the target number
!
!    iquad - integer *8(nnz+1)
!      location in wnear array where quadrature for col_ind(i)
!      starts
!
!    nquad - integer *8
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
!    novers - integer *8(npatches)
!      order of discretization for oversampled sources and
!      density
!
!    ixyzso - integer *8(npatches+1)
!      ixyzso(i) denotes the starting location in srcover,
!      corresponding to patch i
!   
!    nptso - integer *8
!      total number of oversampled points
!
!    srcover - real *8 (12,nptso)
!      oversampled set of source information
!
!    whtsover - real *8 (nptso)
!      smooth quadrature weights at oversampled nodes
!

      implicit none
      integer *8 npatches,norder,npols,npts
      integer *8 ndtarg,ntarg
      integer *8 norders(npatches),ixyzs(npatches+1)
      integer *8 ixyzso(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps
      real *8 targs(ndtarg,ntarg),wts(npts)
      complex *16 zpars(3),alpha,zk
      integer *8 nnz,row_ptr(ntarg+1),col_ind(nnz),nquad
      integer *8 iquad(nnz+1),ncomponents,idcomponent(npatches)
      complex *16 sigma(3*npts+ncomponents),sigma2(3*npts)
	  
      complex *16 wnear(12*nquad)

      integer *8 novers(npatches+1)
      integer *8 nover,npolso,nptso
      real *8 srcover(12,nptso),whtsover(nptso)
      complex *16 pot(3*ntarg+ncomponents)
      complex *16, allocatable :: potsort(:)

      real *8, allocatable :: sources(:,:),targvals(:,:),wtmp2(:)
      complex *16, allocatable :: charges(:),dipvec(:,:),sigmaover(:)
	  complex *16, allocatable :: sigma_aux(:),sigmaover_aux(:)
	  complex *16, allocatable :: pot_aux(:)

      integer *8 ns,nt
      complex *16 beta
      integer *8 ifcharge,ifdipole
      integer *8 ifpgh,ifpghtarg
      complex *16 tmp(10),val,E(4)

      real *8 xmin,xmax,ymin,ymax,zmin,zmax,sizey,sizez,boxsize

      integer *8 i,j,jpatch,jquadstart,jstart

      integer *8 ifaddsub,ifdir

      integer *8 ntj
      
      complex *16 zdotu,pottmp
      complex *16, allocatable :: dtmp2(:,:),ctmp2_b_u(:),ctmp2_b_v(:)
	  complex *16, allocatable :: ctmp2_u(:),ctmp2_v(:),ctmp2_s(:),ctmp2_b_s(:)
      real *8 radexp,epsfmm

      integer *8 ipars(2)
      real *8 dpars(1),timeinfo(10),t1,t2,omp_get_wtime

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      real *8 thresh,ra
      real *8 rr,rmin
      integer *8 nss,ii,l,npover

      integer *8 nd,ntarg0,count1,count2,icount,naux

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
      allocate(sigmaover(3*ns),sigmaover_aux(3*ns))
	  allocate(sigma_aux(3*ntarg+ncomponents),pot_aux(4*ntarg))

! 
!       oversample density
!
		zk=zpars(1)
		alpha=zpars(2)

      call oversample_fun_surf(int(2,8),npatches,norders,ixyzs,iptype,& 
     &npts,sigma(1:npts),novers,ixyzso,ns,sigmaover(1:ns))
      call oversample_fun_surf(int(2,8),npatches,norders,ixyzs,iptype,& 
     &npts,sigma(npts+1:2*npts),novers,ixyzso,ns,sigmaover(ns+1:2*ns))
      call oversample_fun_surf(int(2,8),npatches,norders,ixyzs,iptype,& 
     &npts,sigma(2*npts+1:3*npts),novers,ixyzso,ns,sigmaover(2*ns+1:3*ns))

      ra = 0
      
      !
      !     FMM call to compute nxSkia and -n·curlSkia
      !
	  call get_thresh(srcover,ns,targs,ntarg,thresh)
	  
       ifdir=0
	  call em_adpie_pec_FMM(eps,zpars,ns,npts,srcover,targs,&
      &whtsover,sigmaover(1:ns),sigmaover(ns+1:2*ns),sigmaover(2*ns+1:3*ns),&
	  &pot_aux(1:npts),pot_aux(npts+1:2*npts),pot_aux(2*npts+1:3*npts),&
	  &pot_aux(3*npts+1:4*npts),thresh,ifdir)
	  
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
            pot_aux(i)=pot_aux(i)+wnear(nquad+jquadstart+l-1)*sigma(jstart+l-1+npts)
            pot_aux(i)=pot_aux(i)+wnear(2*nquad+jquadstart+l-1)*sigma(jstart+l-1+2*npts)
			

            pot_aux(i+npts)=pot_aux(i+npts)+wnear(3*nquad+jquadstart+l-1)*sigma(jstart+l-1)
            pot_aux(i+npts)=pot_aux(i+npts)+wnear(4*nquad+jquadstart+l-1)*sigma(jstart+l-1+npts)
            pot_aux(i+npts)=pot_aux(i+npts)+wnear(5*nquad+jquadstart+l-1)*sigma(jstart+l-1+2*npts)

            pot_aux(i+2*npts)=pot_aux(i+2*npts)+wnear(6*nquad+jquadstart+l-1)*sigma(jstart+l-1)
            pot_aux(i+2*npts)=pot_aux(i+2*npts)+wnear(7*nquad+jquadstart+l-1)*sigma(jstart+l-1+npts)
            pot_aux(i+2*npts)=pot_aux(i+2*npts)+wnear(8*nquad+jquadstart+l-1)*sigma(jstart+l-1+2*npts)

            pot_aux(i+3*npts)=pot_aux(i+3*npts)+wnear(9*nquad+jquadstart+l-1)*sigma(jstart+l-1)
            pot_aux(i+3*npts)=pot_aux(i+3*npts)+wnear(10*nquad+jquadstart+l-1)*sigma(jstart+l-1+npts)
            pot_aux(i+3*npts)=pot_aux(i+3*npts)+wnear(11*nquad+jquadstart+l-1)*sigma(jstart+l-1+2*npts)
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
        allocate(ctmp2_v(nss),ctmp2_s(nss),wtmp2(nss))

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
			ctmp2_s(ii)=sigmaover(jstart+l+2*ns)
			wtmp2(ii)=whtsover(jstart+l)
          enddo
        enddo
        call em_adpie_pec_FMM(eps,zpars(1),nss,ntarg0,srctmp2,&
        &targs(:,i),wtmp2,ctmp2_u,ctmp2_v,ctmp2_s,E(1),E(2),E(3),E(4),thresh,&
        &ifdir)
	
		
	    pot_aux(i) = pot_aux(i) - E(1)
	    pot_aux(i+ntarg) = pot_aux(i+ntarg) - E(2)
	    pot_aux(i+2*ntarg) = pot_aux(i+2*ntarg) - E(3)
	    pot_aux(i+3*ntarg) = pot_aux(i+3*ntarg) - E(4)

        deallocate(srctmp2,ctmp2_u,ctmp2_v,ctmp2_s,wtmp2)
      enddo

       do count1=1,3*npts
    	   pot(count1)=pot_aux(count1)
	   enddo
       do count1=1,ncomponents
	       pot(3*npts+count1)=0.0d0
	   enddo
		icount=1
		do count1=1,npatches
		  naux=(norders(count1)+1)*(norders(count1)+2)/2
		  do count2=1,naux
			pot(3*ntarg+idcomponent(count1))=pot(3*ntarg+idcomponent(count1))+&
			&(pot_aux(3*ntarg+icount)-ima*alpha*.5d0*sigma(2*npts+icount))*wts(icount)
			pot(2*ntarg+icount)=pot(2*ntarg+icount)+sigma(3*npts+idcomponent(count1))
			icount=icount+1
		  enddo
		enddo
!		write (*,*) 'flujo gmres: ', pot(3*ntarg+1)
!		write (*,*) 'constante gmres: ', sigma(3*npts+1)

       do count1=1,ncomponents
	       pot(3*npts+count1)=pot(3*npts+count1)-.5d0*sigma(3*npts+count1)
	   enddo

      return
      end subroutine lpcomp_adpie_addsub



subroutine em_adpie_pec_FMM(eps,zpars,ns,nt,srcvals,targvals,wts,a_u,a_v,a_s&
 &,AA_u,AA_v,PHI,nE,thresh,ifdir)
implicit none

!
!  This subroutine computes the far field contribution if the Cki 
!  operator via FMM
!
!  operator:
!
!    nxS_{k}[J]
!
!  input:
!    eps - real * 8
!      epsilon for the fmm call
!
!    zpars - complex *16(3)
!      kernel parameters
!      zpars(1) = k 
!      zpars(2) = alpha
!      zpars(3) = - (not used)
!
!    ns - integer *8
!      number of sources (this is the oversampled set of sources)
!   
!    nt - integer *8 
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
!
!  output:
!    AA_u - complex *16(nt)
!      first component of nxS_{k}[J] along
!      the srcvals(4:6,i) direction
!    AA_v - complex *16(nt)
!      second component of nxS_{k}[J] along
!            the (srcvals(10:12,i) x srcvals(4:6,i)) direction
!            


    !List of calling arguments
    real ( kind = 8 ), intent(in) :: eps
    complex ( kind = 8 ), intent(in) :: zpars(3)
    integer *8, intent(in) :: ns,nt
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),wts(ns)
    real ( kind = 8 ), intent(in) :: targvals(12,nt)
    complex ( kind = 8 ), intent(in) :: a_u(ns),a_v(ns),a_s(ns)
    complex ( kind = 8 ), intent(out) :: AA_u(nt),AA_v(nt),PHI(nt),nE(nt)
    real ( kind = 8 ), intent(in) :: thresh
    integer *8, intent(inout) :: ifdir 

    !List of local variables
	  real ( kind = 8 ), allocatable :: n_vect_s(:,:),u_vect_s(:,:)
	  real ( kind = 8 ), allocatable :: v_vect_s(:,:),source(:,:)
	  real ( kind = 8 ), allocatable :: n_vect_t(:,:),u_vect_t(:,:)
	  real ( kind = 8 ), allocatable :: v_vect_t(:,:),targets(:,:)

    complex ( kind = 8 ), allocatable :: a_vect(:,:),b_vect(:,:),b_vect_aux(:,:)
    complex ( kind = 8 ), allocatable :: lambda(:),rho(:),sigma(:),dipvect(:,:)
    complex ( kind = 8 ), allocatable :: E(:,:),curlE(:,:),divE(:)
    complex ( kind = 8 ) ima,zk,alpha

    integer *8 count1,count2,nd
    integer *8 ifa_vect,ifb_vect,iflambda,ifrho,ifE,ifcurlE,ifdivE
    real ( kind = 8 ) pi
	
    pi=3.141592653589793238462643383d0
  
	  ima=(0.0d0,1.0d0)
    zk=zpars(1)
	  alpha=zpars(2)
!     ifdir=1
    allocate(a_vect(3,ns))
    allocate(b_vect(3,ns))
    allocate(b_vect_aux(3,nt))
    allocate(sigma(ns))
    allocate(dipvect(3,ns))
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
        a_vect(1,count1)=a_u(count1)*u_vect_s(1,count1)+a_v(count1)*v_vect_s(1,count1)
        a_vect(2,count1)=a_u(count1)*u_vect_s(2,count1)+a_v(count1)*v_vect_s(2,count1)
        a_vect(3,count1)=a_u(count1)*u_vect_s(3,count1)+a_v(count1)*v_vect_s(3,count1)
		rho(count1)=-a_s(count1)
		lambda(count1)=ima*alpha*a_s(count1)
		
		b_vect(1,count1)=ima*alpha*(a_u(count1)*v_vect_s(1,count1)-a_v(count1)*u_vect_s(1,count1))
		b_vect(2,count1)=ima*alpha*(a_u(count1)*v_vect_s(2,count1)-a_v(count1)*u_vect_s(2,count1))
		b_vect(3,count1)=ima*alpha*(a_u(count1)*v_vect_s(3,count1)-a_v(count1)*u_vect_s(3,count1))
    enddo

    !Computing the full operator
    ifa_vect=1
    ifb_vect=1
    iflambda=1
    ifrho=1
    ifE=1
    ifcurlE=0
    ifdivE=1

	call Vector_Helmholtz_targ2(eps,zk,ns,source,wts,ifa_vect,a_vect,&
	 &ifb_vect,b_vect,iflambda,lambda,ifrho,rho,n_vect_s,ifE,E,ifcurlE,&
	 &curlE,ifdivE,PHI,nt,targets,thresh,ifdir)


!   do count1=1,nt
!    write (*,*) 'first: ',count1,PHI(count1)
!   enddo
!	call Vector_Helmholtz_targ(eps,zk,ns,source,wts,ifa_vect,a_vect,ifb_vect,&
!	 &b_vect,iflambda,rho_in,ifrho,rho,n_vect_s,ifE,E,ifcurlE,curlE,ifdivE,divE,nt,targets)


    do count1=1,nt
      b_vect_aux(1,count1)=n_vect_t(2,count1)*E(3,count1)-&
      &n_vect_t(3,count1)*E(2,count1)
      b_vect_aux(2,count1)=n_vect_t(3,count1)*E(1,count1)-&
      &n_vect_t(1,count1)*E(3,count1)
      b_vect_aux(3,count1)=n_vect_t(1,count1)*E(2,count1)-&
      &n_vect_t(2,count1)*E(1,count1)
    enddo

    do count1=1,nt
      AA_u(count1)=b_vect_aux(1,count1)*u_vect_t(1,count1)+&
      &b_vect_aux(2,count1)*u_vect_t(2,count1)+b_vect_aux(3,count1)*&
      &u_vect_t(3,count1)
      AA_v(count1)=b_vect_aux(1,count1)*v_vect_t(1,count1)+&
      &b_vect_aux(2,count1)*v_vect_t(2,count1)+b_vect_aux(3,count1)*&
      &v_vect_t(3,count1)
	  
	  nE(count1)=n_vect_t(1,count1)*E(1,count1)+&
	  &n_vect_t(2,count1)*E(2,count1)+n_vect_t(3,count1)*E(3,count1)
    enddo


        !Computing the div vector part
    ifa_vect=0
    ifb_vect=1
    iflambda=0
    ifrho=0
    ifE=1
    ifcurlE=0
    ifdivE=1

!!	call Vector_Helmholtz_targ2(eps,zk,ns,source,wts,ifa_vect,a_vect,&
!!	 &ifb_vect,b_vect,iflambda,lambda,ifrho,rho,n_vect_s,ifE,E,ifcurlE,&
!!	 &curlE,ifdivE,divE,nt,targets,thresh,ifdir)


!!   do count1=1,ns
!!    sigma(count1)=-ima*alpha*zk**2*a_s(count1)*wts(count1)
!!    dipvect(1,count1)=n_vect_s(1,count1)*a_s(count1)*wts(count1)
!!    dipvect(2,count1)=n_vect_s(2,count1)*a_s(count1)*wts(count1)
!!    dipvect(3,count1)=n_vect_s(3,count1)*a_s(count1)*wts(count1)
!!  enddo

!!  do count1=1,nt
!!    PHI(count1)=0.0d0
!!  enddo

   nd=1
   if (ifdir.eq.1) then
!!    call h3ddirectcdp(nd,zk,source,sigma,dipvect,ns&
!!   &,targets,nt,PHI,thresh)
  else
!!  call hfmm3d_t_cd_p_vec(nd,eps,zk,ns,source,sigma,dipvect&
!!  &,nt,targets,PHI)
  endif


  do count1=1,nt
!!    PHI(count1)=PHI(count1)/(4.0d0*pi)+divE(count1)
  enddo


!  do count1=1,nt
!    write (*,*) 'second: ',count1,PHI(count1)
!  enddo

!  write (*,*) 'this is nt: ',nt

!    read (*,*)

    deallocate(a_vect)
    deallocate(b_vect)
    deallocate(b_vect_aux)
    deallocate(lambda)
    deallocate(rho)
    deallocate(E)
    deallocate(curlE)
    deallocate(divE)
    deallocate(sigma)
    deallocate(dipvect)

    deallocate(u_vect_s)
    deallocate(v_vect_s)
    deallocate(n_vect_s)
    deallocate(source)

    deallocate(u_vect_t)
    deallocate(v_vect_t)
    deallocate(n_vect_t)
    deallocate(targets)

return
end subroutine em_adpie_pec_FMM

subroutine 	get_rhs_adpie_pec(P0,vf,ns,srcvals,zk,RHS,&
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
!    ns - integer *8
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
	integer *8, intent(in) :: ns
	real ( kind = 8 ), intent(in) :: P0(3)
	complex ( kind = 8 ), intent(in) :: vf(3)
	real ( kind = 8 ), intent(in) :: srcvals(12,ns)
	complex ( kind = 8 ), intent(in) :: zk
	integer *8, intent(in) :: npatches,ncomponents
	integer *8, intent(in) :: idcomponent(npatches),norders(npatches)
	integer *8, intent(in) :: ixyzs(npatches+1),iptype(npatches)
	complex ( kind = 8 ), intent(out) :: RHS(3*ns+ncomponents)
	
	!List of local variables
	complex ( kind = 8 ), allocatable :: E(:,:), H(:,:),divE(:),nE(:)
	integer *8 count1,count2,icount,naux
	real ( kind = 8 ) ru(3),rv(3),cross_aux(3)
	real ( kind = 8 ), allocatable :: wts(:)

	allocate(E(3,ns), H(3,ns))
	allocate(divE(ns),nE(ns))
	allocate(wts(ns))
!	call fieldsMD(zk,P0,srcvals(1:3,:),ns,curlE,E,vf,0)
!	write (*,*) 'P0',P0,ns
!!	call fieldsED(zk,P0,srcvals(1:3,:),ns,E,H,vf,0)
!!	call fieldsMD(zk,P0,srcvals(1:3,:),ns,E,H,vf,1)
	call point_source_vector_helmholtz(ns,P0,vf,srcvals(1:3,:),zk,E,H,divE)	
!	read (*,*)
	do count1=1,ns	
		call orthonormalize(srcvals(4:6,count1),srcvals(10:12,count1),ru,rv)		
		RHS(count1)=-DOT_PRODUCT(rv,E(:,count1))
		RHS(ns+count1)=DOT_PRODUCT(ru,E(:,count1))
		!RHS(2*ns+count1)=0.0d0
		RHS(2*ns+count1)=divE(count1)		
		nE(count1)=DOT_PRODUCT(srcvals(10:12,count1),E(:,count1))
	enddo
    call get_qwts(npatches,norders,ixyzs,iptype,ns,srcvals,wts)
    do count1=1,ncomponents
		RHS(3*ns+count1)=0.0d0
	enddo

	icount=1
    do count1=1,npatches
      naux=(norders(count1)+1)*(norders(count1)+2)/2
	  do count2=1,naux
	    RHS(3*ns+idcomponent(count1))=RHS(3*ns+idcomponent(count1))+nE(icount)*wts(icount)
	    icount=icount+1
	  enddo
    enddo
!	write (*,*) 'valor esperado del flujo: ', RHS(3*ns+1)
return 
end subroutine get_rhs_adpie_pec


subroutine 	get_rhs_edpie_pec(P0,vf,ns,srcvals,zk,RHS,&
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
!    ns - integer *8
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
	integer *8, intent(in) :: ns
	real ( kind = 8 ), intent(in) :: P0(3)
	complex ( kind = 8 ), intent(in) :: vf(3)
	real ( kind = 8 ), intent(in) :: srcvals(12,ns)
	complex ( kind = 8 ), intent(in) :: zk
	integer *8, intent(in) :: npatches,ncomponents
	integer *8, intent(in) :: idcomponent(npatches),norders(npatches)
	integer *8, intent(in) :: ixyzs(npatches+1),iptype(npatches)
	complex ( kind = 8 ), intent(out) :: RHS(3*ns+ncomponents)
	
	!List of local variables
	complex ( kind = 8 ), allocatable :: E(:,:), H(:,:),divE(:),nE(:)
  complex ( kind = 8 ), allocatable :: A(:,:), gradpot(:,:),pot(:)
  
	integer *8 count1,count2,icount,naux
	real ( kind = 8 ) ru(3),rv(3),cross_aux(3)
	real ( kind = 8 ), allocatable :: wts(:)

	allocate(E(3,ns), H(3,ns))
  allocate(A(3,ns),gradpot(3,ns),pot(ns))
	allocate(divE(ns),nE(ns))
	allocate(wts(ns))

  call fieldsEDdpie(zk,P0,srcvals,ns,int(12,8),vf,int(0,8),E,H,A,pot,gradpot)
  call fieldsMDdpie(zk,P0,srcvals,ns,int(12,8),vf,int(1,8),E,H,A,pot,gradpot)
	do count1=1,ns	
		call orthonormalize(srcvals(4:6,count1),srcvals(10:12,count1),ru,rv)		
		RHS(count1)=-DOT_PRODUCT(rv,E(:,count1))
		RHS(ns+count1)=DOT_PRODUCT(ru,E(:,count1))
		RHS(2*ns+count1)=0.0d0
!		RHS(2*ns+count1)=divE(count1)		
		nE(count1)=DOT_PRODUCT(srcvals(10:12,count1),E(:,count1))
	enddo
  call get_qwts(npatches,norders,ixyzs,iptype,ns,srcvals,wts)
  do count1=1,ncomponents
		RHS(3*ns+count1)=0.0d0
	enddo

	icount=1
  do count1=1,npatches
    naux=(norders(count1)+1)*(norders(count1)+2)/2
	  do count2=1,naux
	    RHS(3*ns+idcomponent(count1))=RHS(3*ns+idcomponent(count1))+nE(icount)*wts(icount)
	    icount=icount+1
	  enddo
  enddo
!	write (*,*) 'valor esperado del flujo: ', RHS(3*ns+1)
return 
end subroutine get_rhs_edpie_pec


subroutine 	get_rhs_dpie_pec(P0,vf,ns,srcvals,zk,RHS,rhs_s,&
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
!    ns - integer *8
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
	integer *8, intent(in) :: ns
	real ( kind = 8 ), intent(in) :: P0(3)
	complex ( kind = 8 ), intent(in) :: vf(3)
	real ( kind = 8 ), intent(in) :: srcvals(12,ns)
	complex ( kind = 8 ), intent(in) :: zk
	integer *8, intent(in) :: npatches,ncomponents
	integer *8, intent(in) :: idcomponent(npatches),norders(npatches)
	integer *8, intent(in) :: ixyzs(npatches+1),iptype(npatches)
	complex ( kind = 8 ), intent(out) :: RHS(3*ns+ncomponents),rhs_s(ns+ncomponents)
	
	!List of local variables
	complex ( kind = 8 ), allocatable :: E(:,:), H(:,:),divE(:),nA(:)
  complex ( kind = 8 ), allocatable :: A(:,:), gradpot(:,:),pot(:)
  complex ( kind = 8 ), allocatable :: dpot_normal(:)
  complex ( kind = 8 ) ima  
	integer *8 count1,count2,icount,naux
	real ( kind = 8 ) ru(3),rv(3),cross_aux(3)
	real ( kind = 8 ), allocatable :: wts(:)

	ima=(0.0d0,1.0d0)

	allocate(E(3,ns), H(3,ns))
  allocate(A(3,ns),gradpot(3,ns),pot(ns))
	allocate(divE(ns),nA(ns),dpot_normal(ns))
	allocate(wts(ns))

  call fieldsEDdpie(zk,P0,srcvals,ns,int(12,8),vf,int(0,8),E,H,A,pot,gradpot)
  call fieldsMDdpie(zk,P0,srcvals,ns,int(12,8),vf,int(1,8),E,H,A,pot,gradpot)
	do count1=1,ns	
		call orthonormalize(srcvals(4:6,count1),srcvals(10:12,count1),ru,rv)		
		RHS(count1)=-DOT_PRODUCT(rv,A(:,count1))
		RHS(ns+count1)=DOT_PRODUCT(ru,A(:,count1))
		RHS(2*ns+count1)=-ima*zk*pot(count1)
!		RHS(2*ns+count1)=divE(count1)
		nA(count1)=DOT_PRODUCT(srcvals(10:12,count1),A(:,count1))
    dpot_normal(count1)=DOT_PRODUCT(srcvals(10:12,count1),gradpot(:,count1))
    rhs_s(count1)=pot(count1)
	enddo
  call get_qwts(npatches,norders,ixyzs,iptype,ns,srcvals,wts)
  do count1=1,ncomponents
		RHS(3*ns+count1)=0.0d0
    rhs_s(ns+count1)=0.0d0
	enddo

	icount=1
  do count1=1,npatches
    naux=(norders(count1)+1)*(norders(count1)+2)/2
	  do count2=1,naux
	    RHS(3*ns+idcomponent(count1))=RHS(3*ns+idcomponent(count1))+nA(icount)*wts(icount)
      rhs_s(ns+idcomponent(count1))=rhs_s(ns+idcomponent(count1))+dpot_normal(icount)*wts(icount)
      icount=icount+1
	  enddo
  enddo
!	write (*,*) 'valor esperado del flujo: ', RHS(3*ns+1),rhs_s(ns+1)
return 
end subroutine get_rhs_dpie_pec

subroutine test_accuracy_adpie(eps_FMM,sol,zk,alpha,ns,wts,srcvals,P0,vf,Pt)
implicit none

    !List of calling arguments
	complex ( kind = 8 ), intent(in) :: zk,alpha
    integer *8, intent(in) :: ns
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),eps_FMM
    real ( kind = 8 ), intent(in) :: wts(ns),P0(3),Pt(3)
	complex ( kind = 8 ), intent(in) :: sol(3*ns),vf(3)
	
    !List of local variables
	complex ( kind = 8 ) a_u00,a_v00
	complex ( kind = 8 ) ima, R1, R2,Et1(3),Et2(3),aux_cmp,Ht2(3),phi1(1),E_aux(3),phi2
	real ( kind = 8 ) ru_s(3),rv_s(3),n_s(3),sour(3),r,dr(3)
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3),error,rel_err
	real ( kind = 8 ) pi

	integer *8 count1
		ima=(0.0d0,1.0d0)
		pi=3.1415926535897932384626433832795028841971d0

!		write (*,*) 'P0',P0
call EDPIE1_FMM_targ(eps_FMM,zk,alpha,ns,srcvals,int(1,8),Pt,wts,sol(1:ns),sol(ns+1:2*ns),sol(2*ns+1:3*ns),Et1,phi1)
		
!		call fieldsED(zk,P0,Pt,1,Et2,Ht2,vf,0)
!		call fieldsMD(zk,P0,Pt,1,Et2,Ht2,vf,1)
!	    write (*,*) 'Antes: ',P0,vf,Pt,zk
		call point_source_vector_helmholtz(int(1,8),P0,vf,Pt,zk,Et2,E_aux,phi2)
!		write (*,*) Et2

		error=sqrt(abs(Et1(1)-Et2(1))**2+abs(Et1(2)-Et2(2))**2+abs(Et1(3)-Et2(3))**2)
!		write (*,*) 'Error: ', error
		write (*,*) 'Relative Error: ', error/sqrt(abs(Et2(1))**2+abs(Et2(2))**2+abs(Et2(3))**2)
		write (*,*) 'relative error in phi: ',abs(phi1-phi2)/abs(phi2)
!		write (*,*) 'Et1: ', Et1
!		write (*,*) 'Et2: ', Et2
!		write (*,*) 'phi1: ', phi1, 'phi2: ', phi2

return
end subroutine test_accuracy_adpie




subroutine test_accuracy_edpie(eps_FMM,sol,zk,alpha,ns,wts,srcvals,P0,vf,Pt)
implicit none

  !List of calling arguments
	complex ( kind = 8 ), intent(in) :: zk,alpha
  integer *8, intent(in) :: ns
  real ( kind = 8 ), intent(in) :: srcvals(12,ns),eps_FMM
  real ( kind = 8 ), intent(in) :: wts(ns),P0(3),Pt(3)
	complex ( kind = 8 ), intent(in) :: sol(3*ns),vf(3)
	
  !List of local variables
	complex ( kind = 8 ) a_u00,a_v00
	complex ( kind = 8 ) ima, R1, R2,Et1(3),Et2(3),aux_cmp,Ht2(3),phi1(1),E_aux(3),phi2
  complex ( kind = 8 ) At2(3),pott2,gradpott2(3)
	real ( kind = 8 ) ru_s(3),rv_s(3),n_s(3),sour(3),r,dr(3)
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3),error,rel_err
	real ( kind = 8 ) pi

	integer *8 count1
	ima=(0.0d0,1.0d0)
	pi=3.1415926535897932384626433832795028841971d0

	write (*,*) 'P0',P0
  call EDPIE1_FMM_targ(eps_FMM,zk,alpha,ns,srcvals,int(1,8),Pt,wts,sol(1:ns),sol(ns+1:2*ns),sol(2*ns+1:3*ns),Et1,phi1)
		
!		call fieldsED(zk,P0,Pt,1,Et2,Ht2,vf,0)
!		call fieldsMD(zk,P0,Pt,1,Et2,Ht2,vf,1)
!	write (*,*) 'Antes: ',P0,vf,Pt,zk

  call fieldsEDdpie(zk,P0,Pt,int(1,8),int(3,8),vf,int(0,8),Et2,Ht2,At2,pott2,gradpott2)
  call fieldsMDdpie(zk,P0,Pt,int(1,8),int(3,8),vf,int(1,8),Et2,Ht2,At2,pott2,gradpott2)
!	call point_source_vector_helmholtz(1,P0,vf,Pt,zk,Et2,E_aux,phi2)
!	write (*,*) 'Et2: ', Et2
!	write (*,*) 'Et1: ', Et1

	error=sqrt(abs(Et1(1)-Et2(1))**2+abs(Et1(2)-Et2(2))**2+abs(Et1(3)-Et2(3))**2)
	write (*,*) 'Relative Error: ', error/sqrt(abs(Et2(1))**2+abs(Et2(2))**2+abs(Et2(3))**2)
!	write (*,*) 'relative error in phi: ',abs(phi1-phi2)/abs(phi2)

return
end subroutine test_accuracy_edpie

subroutine test_accuracy_dpie(eps_FMM,sol,sol_s,zk,alpha,ns,wts,srcvals,P0,vf,Pt)
implicit none

  !List of calling arguments
	complex ( kind = 8 ), intent(in) :: zk,alpha
  integer *8, intent(in) :: ns
  real ( kind = 8 ), intent(in) :: srcvals(12,ns),eps_FMM
  real ( kind = 8 ), intent(in) :: wts(ns),P0(3),Pt(3)
	complex ( kind = 8 ), intent(in) :: sol(3*ns),sol_s(ns),vf(3)
	
  !List of local variables
	complex ( kind = 8 ) a_u00,a_v00
	complex ( kind = 8 ) ima, R1, R2

  complex ( kind = 8 ) Et1(3),Et2(3),Ht1(3),Ht2(3),At1(3),At2(3)
  complex ( kind = 8 ) pott1(1),pott2(1),gradpott1(3),gradpott2(3)
	
  real ( kind = 8 ) ru_s(3),rv_s(3),n_s(3),sour(3),r,dr(3)
	real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3),error,rel_err
	real ( kind = 8 ) pi

	integer *8 count1
	ima=(0.0d0,1.0d0)
	pi=3.1415926535897932384626433832795028841971d0

	write (*,*) 'P0',P0
  call dpie_FMM_targ(eps_FMM,zk,alpha,ns,srcvals,int(1,8),Pt,wts,sol(1:ns),sol(ns+1:2*ns),sol(2*ns+1:3*ns),At1,pott1,Ht1)
		
  call sdpie_grad_FMM_targ(eps_FMM,zk,alpha,ns,srcvals,int(1,8),Pt,wts,sol_s,pott1,gradpott1)

!		call fieldsED(zk,P0,Pt,1,Et2,Ht2,vf,0)
!		call fieldsMD(zk,P0,Pt,1,Et2,Ht2,vf,1)
!	write (*,*) 'Antes: ',P0,vf,Pt,zk

  call fieldsEDdpie(zk,P0,Pt,int(1,8),int(3,8),vf,int(0,8),Et2,Ht2,At2,pott2,gradpott2)
  call fieldsMDdpie(zk,P0,Pt,int(1,8),int(3,8),vf,int(1,8),Et2,Ht2,At2,pott2,gradpott2)
!	write (*,*) 'Et2: ', Et2
!	write (*,*) 'Et1: ', Et1
  Et1(1)=ima*zk*At1(1)+gradpott1(1)
  Et1(2)=ima*zk*At1(2)+gradpott1(2)
  Et1(3)=ima*zk*At1(3)+gradpott1(3)


	error=sqrt(abs(Et1(1)-Et2(1))**2+abs(Et1(2)-Et2(2))**2+abs(Et1(3)-Et2(3))**2)
	write (*,*) 'Relative Error in E: ', error/sqrt(abs(Et2(1))**2+abs(Et2(2))**2+abs(Et2(3))**2)

	error=sqrt(abs(Ht1(1)-Ht2(1))**2+abs(Ht1(2)-Ht2(2))**2+abs(Ht1(3)-Ht2(3))**2)
	write (*,*) 'Relative Error in H: ', error/sqrt(abs(Ht2(1))**2+abs(Ht2(2))**2+abs(Ht2(3))**2)

  error=sqrt(abs(At1(1)-At2(1))**2+abs(At1(2)-At2(2))**2+abs(At1(3)-At2(3))**2)
	write (*,*) 'Relative Error in A: ', error/sqrt(abs(At2(1))**2+abs(At2(2))**2+abs(At2(3))**2)

	write (*,*) 'relative error in phi: ',abs(pott1-pott2)/abs(pott2)

  error=sqrt(abs(gradpott1(1)-gradpott2(1))**2+abs(gradpott1(2)-gradpott2(2))&
   &**2+abs(gradpott1(3)-gradpott2(3))**2)
	
  write (*,*) 'Relative Error in gradpot: ', error/sqrt(abs(gradpott2(1))**2+&
   &abs(gradpott2(2))**2+abs(gradpott2(3))**2)

return
end subroutine test_accuracy_dpie


subroutine EDPIE1_FMM_targ(eps,zk,alpha,ns,srcvals,nt,targ,wts,a_u,a_v,rho_in,E,phi)
implicit none

    !List of calling arguments
	real ( kind = 8 ), intent(in) :: eps
	complex ( kind = 8 ), intent(in) :: zk,alpha
    integer *8, intent(in) :: ns,nt
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),targ(3,nt)
    real ( kind = 8 ), intent(in) :: wts(ns)
    complex ( kind = 8 ), intent(in) :: a_u(ns),a_v(ns),rho_in(ns)
    complex ( kind = 8 ), intent(out) :: E(3,nt),phi(nt)

    !List of local variables
	real ( kind = 8 ), allocatable :: n_vect(:,:),u_vect(:,:),v_vect(:,:),source(:,:)
    complex ( kind = 8 ), allocatable :: a_vect(:,:),b_vect(:,:),lambda(:),rho(:)
    complex ( kind = 8 ), allocatable :: curlE(:,:)
	complex ( kind = 8 ) ima

    integer *8 count1,count2
    integer *8 ifa_vect,ifb_vect,iflambda,ifrho,ifE,ifcurlE,ifdivE

	ima=(0.0d0,1.0d0)

    allocate(a_vect(3,ns))
    allocate(b_vect(3,ns))
    allocate(lambda(ns))
	allocate(rho(ns))
    allocate(curlE(3,nt))
	
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
        a_vect(1,count1)=a_u(count1)*u_vect(1,count1)+a_v(count1)*v_vect(1,count1)
        a_vect(2,count1)=a_u(count1)*u_vect(2,count1)+a_v(count1)*v_vect(2,count1)
        a_vect(3,count1)=a_u(count1)*u_vect(3,count1)+a_v(count1)*v_vect(3,count1)
		rho(count1)=-rho_in(count1)
		lambda(count1)=ima*alpha*rho_in(count1)
		
		b_vect(1,count1)=ima*alpha*(a_u(count1)*v_vect(1,count1)-a_v(count1)*u_vect(1,count1))
		b_vect(2,count1)=ima*alpha*(a_u(count1)*v_vect(2,count1)-a_v(count1)*u_vect(2,count1))
		b_vect(3,count1)=ima*alpha*(a_u(count1)*v_vect(3,count1)-a_v(count1)*u_vect(3,count1))
    enddo


    !Computing the full operator
    ifa_vect=1
    ifb_vect=1
    iflambda=1
    ifrho=1
    ifE=1
    ifcurlE=0
    ifdivE=1

call Vector_Helmholtz_targ(eps,zk,ns,source,wts,ifa_vect,a_vect,ifb_vect,&
 &b_vect,iflambda,lambda,ifrho,rho,n_vect,ifE,E,ifcurlE,curlE,ifdivE,phi,nt,targ)

	deallocate(a_vect)
	deallocate(b_vect)
	deallocate(lambda)
	deallocate(curlE)
	deallocate(rho)
	
	deallocate(u_vect)
	deallocate(v_vect)
	deallocate(n_vect)
	deallocate(source)

return
end subroutine EDPIE1_FMM_targ



subroutine dpie_FMM_targ(eps,zk,alpha,ns,srcvals,nt,targ,wts,a_u,a_v,rho_in,A,pot,H)
implicit none

  !List of calling arguments
	real ( kind = 8 ), intent(in) :: eps
	complex ( kind = 8 ), intent(in) :: zk,alpha
  integer *8, intent(in) :: ns,nt
  real ( kind = 8 ), intent(in) :: srcvals(12,ns),targ(3,nt)
  real ( kind = 8 ), intent(in) :: wts(ns)
  complex ( kind = 8 ), intent(in) :: a_u(ns),a_v(ns),rho_in(ns)
  complex ( kind = 8 ), intent(out) :: A(3,nt),pot(nt),H(3,nt)

  !List of local variables
	real ( kind = 8 ), allocatable :: n_vect(:,:),u_vect(:,:),v_vect(:,:),source(:,:)
  complex ( kind = 8 ), allocatable :: a_vect(:,:),b_vect(:,:),lambda(:),rho(:)
  complex ( kind = 8 ), allocatable :: curlE(:,:)
	complex ( kind = 8 ) ima

  integer *8 count1,count2
  integer *8 ifa_vect,ifb_vect,iflambda,ifrho,ifE,ifcurlE,ifdivE

	  ima=(0.0d0,1.0d0)

    allocate(a_vect(3,ns))
    allocate(b_vect(3,ns))
    allocate(lambda(ns))
	  allocate(rho(ns))
    allocate(curlE(3,nt))
	
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
        a_vect(1,count1)=a_u(count1)*u_vect(1,count1)+a_v(count1)*v_vect(1,count1)
        a_vect(2,count1)=a_u(count1)*u_vect(2,count1)+a_v(count1)*v_vect(2,count1)
        a_vect(3,count1)=a_u(count1)*u_vect(3,count1)+a_v(count1)*v_vect(3,count1)
		rho(count1)=-rho_in(count1)
		lambda(count1)=ima*alpha*rho_in(count1)
		
		b_vect(1,count1)=ima*alpha*(a_u(count1)*v_vect(1,count1)-a_v(count1)*u_vect(1,count1))
		b_vect(2,count1)=ima*alpha*(a_u(count1)*v_vect(2,count1)-a_v(count1)*u_vect(2,count1))
		b_vect(3,count1)=ima*alpha*(a_u(count1)*v_vect(3,count1)-a_v(count1)*u_vect(3,count1))
    enddo


    !Computing the full operator
    ifa_vect=1
    ifb_vect=1
    iflambda=1
    ifrho=1
    ifE=1
    ifcurlE=1
    ifdivE=1

    call Vector_Helmholtz_targ(eps,zk,ns,source,wts,ifa_vect,a_vect,ifb_vect,&
     &b_vect,iflambda,lambda,ifrho,rho,n_vect,ifE,A,ifcurlE,H,ifdivE,pot,nt,targ)

    do count1=1,nt
      pot(count1)=pot(count1)/(-ima*zk)
      write (*,*) count1,pot(count1),zk,ima
    enddo

	  deallocate(a_vect)
	  deallocate(b_vect)
	  deallocate(lambda)
	  deallocate(curlE)
	  deallocate(rho)
	
	  deallocate(u_vect)
	  deallocate(v_vect)
	  deallocate(n_vect)
	  deallocate(source)

return
end subroutine dpie_FMM_targ
