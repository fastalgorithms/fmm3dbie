!                                                                       1
      subroutine getnearquad_CKi(npatches,norders,&
       ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
       ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,&
       iquad,rfac0,nquad,wnear)
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
!    iquadtype - integer
!      quadrature type
!      iquadtype = 1, use ggq for self + adaptive integration
!      for rest
!
!    nnz - integer
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
!    rfac0 - integer
!      radius parameter for near field
!
!    nquad - integer
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
      external  fker_em_auCKi_pec
      external  em_auCKi_pec
      ndz=3
      ndd=1
      ndi=2
      ipv=1

!     write (*,*) 'Datos Entrada'
      fker =>  fker_em_auCKi_pec
      fker =>  em_auCKi_pec
!	  	 write(*,*) ndtarg,ntarg
!	 write (*,*) 'Fin Entrada'

!		 read (*,*)


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
      ipars(1)=2
      ipars(2)=1
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
     &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
     &wnear((icount-1)*nquad+1:(icount)*nquad))
	 
      icount=4
      ipars(1)=2
      ipars(2)=2
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
     &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
     &wnear((icount-1)*nquad+1:(icount)*nquad))
	 
	 
      icount=5	 
      ipars(1)=3
      ipars(2)=3
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
     &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
     &wnear((icount-1)*nquad+1:(icount)*nquad))
	 
      icount=6
      ipars(1)=4
      ipars(2)=3
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
     &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
     &wnear((icount-1)*nquad+1:(icount)*nquad))
	 
      icount=7
      ipars(1)=3
      ipars(2)=4
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
     &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
     &wnear((icount-1)*nquad+1:(icount)*nquad))

      icount=8
      ipars(1)=4
      ipars(2)=4
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,& 
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
     &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
     &wnear((icount-1)*nquad+1:(icount)*nquad))

      icount=9
      ipars(1)=1
      ipars(2)=4
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,& 
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
     &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
     &wnear((icount-1)*nquad+1:(icount)*nquad))

      icount=10
      ipars(1)=2
      ipars(2)=4
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,& 
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
     &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
     &wnear((icount-1)*nquad+1:(icount)*nquad))
	 
      icount=11
      ipars(1)=1
      ipars(2)=3
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,& 
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
     &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
     &wnear((icount-1)*nquad+1:(icount)*nquad))
	 
      icount=12
      ipars(1)=3
      ipars(2)=1
      call zgetnearquad_ggq_guru(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,& 
     &ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
     &ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,&
     &wnear((icount-1)*nquad+1:(icount)*nquad))
	
      return
      end subroutine getnearquad_CKi


subroutine em_auCKi_pec_solver(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,eps,zpars,numit,ifinout,&
     &rhs,eps_gmres,niter,errs,rres,soln,sigma2,scalar_RHS)

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
!    zpars - complex *16 (3)
!      kernel parameters (Referring to formula (1))
!      zpars(1) = k 
!      zpars(2) = alpha
!      zpars(3) = - (not used)
!
!    ifinout - integer
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
!    numit - integer
!      max number of gmres iterations
!
!  output
!    niter - integer
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
      integer *8 npatches,norder,npols,npts
      integer *8 ifinout
      integer *8 norders(npatches),ixyzs(npatches+1)
      integer *8 iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,eps_gmres
      complex *16 zpars(3)
      complex *16 rhs(2*npts)
      complex *16 soln(2*npts),sigma2(3*npts),scalar_RHS(npts)

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
      real *8, allocatable :: srcover(:,:),wover(:)
      integer *8, allocatable :: ixyzso(:),novers(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:) 

      integer *8 i,j,jpatch,jquadstart,jstart

      integer *8 ipars
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


      call getnearquad_CKi(npatches,norders,&
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

        call lpcomp_CKi_addsub(npatches,norders,ixyzs,&
        &iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,&
        &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,&
        &vmat(1,it),novers,npts_over,ixyzso,srcover,wover,wtmp,sigma2,&
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
!        replace this routine by appropriate layer potential
!        evaluation routine  
!
	 
          call lpcomp_CKi2RHS_addsub(npatches,norders,ixyzs,&
          &iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,&
          &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,&
          &soln,novers,npts_over,ixyzso,srcover,wover,wtmp,sigma2,&
          &scalar_RHS,wnear)

            
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
      end subroutine em_auCKi_pec_solver


subroutine fker_em_auCKi_pec(srcinfo, ndt,targinfo,ndd, dpars,ndz,&
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
    complex ( kind = 8 ) nxcurlSk0a(2,2),nxcurlSk1a(2,2)
    complex ( kind = 8 ) nxcurlcurlSk0a(2,2)
    complex ( kind = 8 ) nxgradSk0lambda(2,1),ncurlSk1a(1,2)
    complex ( kind = 8 ) ncurlSk0a(1,2)
    complex ( kind = 8 ) Sk1lambda(1,1),Sk0lambda(1,1),E_mat(4,4)
    real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3),xprod_aux3(3)
    real ( kind = 8 ) xprod_aux4(3)
    complex ( kind = 8 ) R1_0,R1_1,R2_0,R2_1,ima,my_exp_0,my_exp_1,zk
    complex ( kind = 8 ) alpha,zk0,zk1,nxSk0b(2,2),nxSk1b(2,2)
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
    zk0=zk
    zk1=ima*zk
    R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
    my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
    R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
    my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)

    call orthonormalize(srcinfo(4:6),n_s,ru_s,rv_s)
    call orthonormalize(targinfo(4:6),n_t,ru_t,rv_t)

    call get_nxcurlSka(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_0,my_exp_0,r,&
    &nxcurlSk0a)		
    E_mat(1,1)=nxcurlSk0a(1,1)
    E_mat(1,2)=nxcurlSk0a(1,2)
    E_mat(2,1)=nxcurlSk0a(2,1)
    E_mat(2,2)=nxcurlSk0a(2,2)

	call get_Sklambda(my_exp_0,r,Sk0lambda)
		
	E_mat(3,1)=Sk0lambda(1,1)		
		
	call get_Sklambda(my_exp_1,r,Sk1lambda)
		
	E_mat(1,3)=Sk1lambda(1,1)

	call get_nxgradSklambda(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_0,my_exp_0,r,nxgradSk0lambda)

	E_mat(3,3)=nxgradSk0lambda(1,1)
	E_mat(4,3)=nxgradSk0lambda(2,1)

	call get_ncurlSka(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_1,my_exp_1,r,ncurlSk1a)
		
	E_mat(3,4)=-ncurlSk1a(1,1)
	E_mat(4,4)=-ncurlSk1a(1,2)

	call get_ncurlSka(ru_s,rv_s,n_s,ru_t,rv_t,n_t,dr,R1_0,my_exp_0,r,ncurlSk0a)
		
	E_mat(1,4)=ncurlSk0a(1,1)
	E_mat(2,4)=ncurlSk0a(1,2)

	E_val=E_mat(ipars(1),ipars(2))

return
end subroutine fker_em_auCKi_pec

      subroutine lpcomp_CKi_addsub(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,sigma,novers,&
     &nptso,ixyzso,srcover,whtsover,pot,sigma2,&
     &wnear)
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
!    zpars - complex *16 (3)
!      kernel parameters (Referring to formula (1))
!      zpars(1) = k 
!      zpars(2) = alpha
!      zpars(3) = - ( not used )
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
!    novers - integer(npatches)
!      order of discretization for oversampled sources and
!      density
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
      integer *8 npatches,norder,npols,npts
      integer *8 ndtarg,ntarg
      integer *8 norders(npatches),ixyzs(npatches+1)
      integer *8 ixyzso(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps
      real *8 targs(ndtarg,ntarg)
      complex *16 zpars(3),alpha,zk
      integer *8 nnz,row_ptr(ntarg+1),col_ind(nnz),nquad
      integer *8 iquad(nnz+1)
      complex *16 sigma(2*npts),sigma2(3*npts)
	  
      complex *16 wnear(12*nquad)

      integer *8 novers(npatches+1)
      integer *8 nover,npolso,nptso
      real *8 srcover(12,nptso),whtsover(nptso)
      complex *16 pot(2*ntarg)
      complex *16, allocatable :: potsort(:)

      real *8, allocatable :: sources(:,:),targvals(:,:),wtmp2(:)
      complex *16, allocatable :: charges(:),dipvec(:,:),sigmaover(:)
	  complex *16, allocatable :: sigma_aux(:),sigmaover_aux(:)
	  complex *16, allocatable :: pot_aux(:)

      integer *8 ns,nt
      complex *16 beta
      integer *8 ifcharge,ifdipole
      integer *8 ifpgh,ifpghtarg
      complex *16 tmp(10),val,E(3)

      real *8 xmin,xmax,ymin,ymax,zmin,zmax,sizey,sizez,boxsize

      integer *8 i,j,jpatch,jquadstart,jstart

      integer *8 ifaddsub,ifdir

      integer *8 ntj
      
      complex *16 zdotu,pottmp
      complex *16, allocatable :: dtmp2(:,:),ctmp2_b_u(:),ctmp2_b_v(:)
	  complex *16, allocatable :: ctmp2_a_u(:),ctmp2_a_v(:),ctmp2_b_s(:)
      real *8 radexp,epsfmm

      integer *8 ipars(2)
      real *8 dpars(1),timeinfo(10),t1,t2,omp_get_wtime

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      real *8 thresh,ra
      real *8 rr,rmin
      integer *8 nss,ii,l,npover

      integer *8 nd,ntarg0

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
	  allocate(sigma_aux(3*ntarg),pot_aux(3*ntarg))

! 
!       oversample density
!
		zk=zpars(1)
		alpha=zpars(2)

      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
     &npts,sigma(1:npts),novers,ixyzso,ns,sigmaover(1:ns))
      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
     &npts,sigma(npts+1:2*npts),novers,ixyzso,ns,sigmaover(ns+1:2*ns))

      ra = 0
      
      !
      !     FMM call to compute nxSkia and -n·curlSkia
      !
	  call get_fmm_thresh(12,ns,srcover,ndtarg,ntarg,targs,thresh)
	  
       ifdir=0
	   
	  call em_auCKi_pec_reg_FMM(eps,zpars(1),ns,npts,srcover,targs,&
      &whtsover,sigmaover(1:ns),sigmaover(ns+1:2*ns),pot_aux(1:npts),&
      &pot_aux(npts+1:2*npts),pot_aux(2*npts+1:3*npts),thresh,ifdir)
	  
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
            pot_aux(i+2*npts)=pot_aux(i+2*npts)+&
            &wnear(6*nquad+jquadstart+l-1)*sigma(jstart+l-1)
            pot_aux(i+2*npts)=pot_aux(i+2*npts)+&
            &wnear(7*nquad+jquadstart+l-1)*sigma(jstart+l-1+npts)		  
		  
            call orthonormalize(srcvals(4:6,jstart+l-1),&
            &srcvals(10:12,jstart+l-1),u_vect,v_vect)
            sigma_x=(sigma(jstart+l-1)*u_vect(1)+sigma(jstart+l-1+npts)*&
            &v_vect(1))
            sigma_y=(sigma(jstart+l-1)*u_vect(2)+sigma(jstart+l-1+npts)*&
            &v_vect(2))
            sigma_z=(sigma(jstart+l-1)*u_vect(3)+sigma(jstart+l-1+npts)*&
            &v_vect(3))

            pot_x=wnear(10*nquad+jquadstart+l-1)*sigma_x
            pot_y=wnear(10*nquad+jquadstart+l-1)*sigma_y
            pot_z=wnear(10*nquad+jquadstart+l-1)*sigma_z
		  
            n_vect=srcvals(10:12,i)
            pot2_x=n_vect(2)*pot_z-n_vect(3)*pot_y
            pot2_y=n_vect(3)*pot_x-n_vect(1)*pot_z
            pot2_z=n_vect(1)*pot_y-n_vect(2)*pot_x
            call orthonormalize(srcvals(4:6,i),srcvals(10:12,i),&
            &u_vect,v_vect)

            pot_aux(i) = pot_aux(i) + pot2_x*u_vect(1)+pot2_y*u_vect(2)+&
            &pot2_z*u_vect(3)
            pot_aux(i+npts) = pot_aux(i+npts) + pot2_x*v_vect(1)+pot2_y*&
            &v_vect(2)+pot2_z*v_vect(3)
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
        allocate(srctmp2(12,nss),ctmp2_a_u(nss))
        allocate(ctmp2_a_v(nss),wtmp2(nss))

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
			wtmp2(ii)=whtsover(jstart+l)
          enddo
        enddo

        call em_auCKi_pec_reg_FMM(eps,zpars(1),nss,ntarg0,srctmp2,&
        &targs(:,i),wtmp2,ctmp2_a_u,ctmp2_a_v,E(1),E(2),E(3),thresh,&
        &ifdir)
		
	    pot_aux(i) = pot_aux(i) - E(1)
	    pot_aux(i+ntarg) = pot_aux(i+ntarg) - E(2)
	    pot_aux(i+2*ntarg) = pot_aux(i+2*ntarg) - E(3)

        deallocate(srctmp2,ctmp2_a_u,ctmp2_a_v,wtmp2)
      enddo
	  	  

      do i=1,ntarg
	    sigma_aux(i)=pot_aux(i)
	    sigma_aux(i+ntarg)=pot_aux(i+ntarg)
	    sigma_aux(i+2*ntarg)=pot_aux(i+2*ntarg)

	    pot_aux(i)=0.0d0
	    pot_aux(i+ntarg)=0.0d0
	    pot_aux(i+2*ntarg)=0.0d0

	    sigma2(i)=sigma_aux(i)
	    sigma2(i+ntarg)=sigma_aux(i+ntarg)
	    sigma2(i+2*ntarg)=sigma_aux(i+2*ntarg)
	  
	    sigma_aux(i)=ima*alpha*zk**2*sigma_aux(i)
	    sigma_aux(i+ntarg)=ima*alpha*zk**2*sigma_aux(i+ntarg)
	    sigma_aux(i+2*ntarg)=ima*alpha*sigma_aux(i+2*ntarg)
      enddo


      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
     &npts,sigma_aux(1:npts),novers,ixyzso,ns,sigmaover_aux(1:ns))
      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
     &npts,sigma_aux(npts+1:2*npts),novers,ixyzso,ns,&
     &sigmaover_aux(ns+1:2*ns))
      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
     &npts,sigma_aux(2*npts+1:3*npts),novers,ixyzso,ns,&
     &sigmaover_aux(2*ns+1:3*ns))

		! FMM call, Here I compute nxSkb
		ifdir=0
	  call CKif_good_FMM(eps,zpars(1),ns,npts,srcover,targs,whtsover,&
     &sigmaover(1:ns),sigmaover(ns+1:2*ns),sigmaover_aux(1:ns),&
     &sigmaover_aux(ns+1:2*ns),sigmaover_aux(2*ns+1:3*ns),&
     &pot_aux(1:npts),pot_aux(npts+1:2*npts),thresh,ifdir)
	  

      do i=1,ntarg

        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch)
          do l=1,npols		  
            call orthonormalize(srcvals(4:6,jstart+l-1),&
            &srcvals(10:12,jstart+l-1),u_vect,v_vect)
            sigma_x=(sigma_aux(jstart+l-1)*u_vect(1)+&
            &sigma_aux(jstart+l-1+npts)*v_vect(1))
            sigma_y=(sigma_aux(jstart+l-1)*u_vect(2)+&
            &sigma_aux(jstart+l-1+npts)*v_vect(2))
            sigma_z=(sigma_aux(jstart+l-1)*u_vect(3)+&
            &sigma_aux(jstart+l-1+npts)*v_vect(3))

            pot_x=wnear(11*nquad+jquadstart+l-1)*sigma_x
            pot_y=wnear(11*nquad+jquadstart+l-1)*sigma_y
            pot_z=wnear(11*nquad+jquadstart+l-1)*sigma_z
		  
            n_vect=srcvals(10:12,i)
            pot2_x=n_vect(2)*pot_z-n_vect(3)*pot_y
            pot2_y=n_vect(3)*pot_x-n_vect(1)*pot_z
            pot2_z=n_vect(1)*pot_y-n_vect(2)*pot_x
            call orthonormalize(srcvals(4:6,i),srcvals(10:12,i),&
            &u_vect,v_vect)

            pot_aux(i) = pot_aux(i) +&
            & pot2_x*u_vect(1)+pot2_y*u_vect(2)+pot2_z*u_vect(3)
            pot_aux(i+npts) = pot_aux(i+npts) +&
            & pot2_x*v_vect(1)+pot2_y*v_vect(2)+pot2_z*v_vect(3)
			
            pot_aux(i) = pot_aux(i) +&
            & wnear(4*nquad+jquadstart+l-1)*sigma_aux(jstart+l-1+2*npts)
            pot_aux(i+npts) = pot_aux(i+npts) +&
            & wnear(5*nquad+jquadstart+l-1)*sigma_aux(jstart+l-1+2*npts)

            pot_aux(i) = pot_aux(i) + &
            &wnear(jquadstart+l-1)*sigma(jstart+l-1)
            pot_aux(i) = pot_aux(i) +&
            & wnear(1*nquad+jquadstart+l-1)*sigma(jstart+l-1+npts)
            pot_aux(i+npts) = pot_aux(i+npts) +&
            & wnear(2*nquad+jquadstart+l-1)*sigma(jstart+l-1)
            pot_aux(i+npts) = pot_aux(i+npts) +&
            & wnear(3*nquad+jquadstart+l-1)*sigma(jstart+l-1+npts)


          enddo
        enddo
      enddo

	ifdir=1

	  do i=1,ntarg
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          nss = nss + ixyzso(jpatch+1)-ixyzso(jpatch)
        enddo
        allocate(srctmp2(12,nss),ctmp2_b_u(nss),ctmp2_b_v(nss))
        allocate(ctmp2_b_s(nss),ctmp2_a_u(nss),ctmp2_a_v(nss))
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
			ctmp2_b_u(ii)=sigmaover_aux(jstart+l)
			ctmp2_b_v(ii)=sigmaover_aux(jstart+l+ns)
			ctmp2_b_s(ii)=sigmaover_aux(jstart+l+2*ns)
			
			
			ctmp2_a_u(ii)=sigmaover(jstart+l)
			ctmp2_a_v(ii)=sigmaover(jstart+l+ns)

			wtmp2(ii)=whtsover(jstart+l)
          enddo
        enddo
		
        call CKif_good_FMM(eps,zpars(1),nss,ntarg0,srctmp2,targs(:,i),&
        &wtmp2,ctmp2_a_u,ctmp2_a_v,ctmp2_b_u,ctmp2_b_v,ctmp2_b_s,&
        &E(1),E(2),thresh,ifdir)
	    pot_aux(i) = pot_aux(i) - E(1)
	    pot_aux(i+ntarg) = pot_aux(i+ntarg) - E(2)

        deallocate(srctmp2,ctmp2_a_u,ctmp2_a_v)
        deallocate(ctmp2_b_s,ctmp2_b_u,ctmp2_b_v,wtmp2)
      enddo
	  
	  do i=1,npts
	    pot(i)=pot_aux(i)
	    pot(i+npts)=pot_aux(i+npts)
	  enddo


      return
      end subroutine lpcomp_CKi_addsub



subroutine CKif_good_FMM(eps,zk,ns,nt,srcvals,targvals,wts,a_u,a_v,&
 &b_u,b_v,rho_in,AA_u,AA_v,thresh,ifdir)
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
    complex ( kind = 8 ), intent(in) :: zk
    integer *8, intent(in) :: ns,nt
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),wts(ns)
    real ( kind = 8 ), intent(in) :: targvals(12,nt)
    complex ( kind = 8 ), intent(in) :: a_u(ns),a_v(ns),b_u(ns),b_v(ns)
    complex ( kind = 8 ), intent(in) :: rho_in(ns)
    complex ( kind = 8 ), intent(out) :: AA_u(nt),AA_v(nt)
    real ( kind = 8 ), intent(in) :: thresh
    integer *8, intent(in) :: ifdir 

    !List of local variables
	real ( kind = 8 ), allocatable :: n_vect_s(:,:),u_vect_s(:,:)
	real ( kind = 8 ), allocatable :: v_vect_s(:,:),source(:,:)
	real ( kind = 8 ), allocatable :: n_vect_t(:,:),u_vect_t(:,:)
	real ( kind = 8 ), allocatable :: v_vect_t(:,:),targets(:,:)

    complex ( kind = 8 ), allocatable :: a_vect(:,:),b_vect(:,:)
    complex ( kind = 8 ), allocatable :: b_vect_t(:,:)
    complex ( kind = 8 ), allocatable :: lambda(:),rho(:)
    complex ( kind = 8 ), allocatable :: E(:,:),curlE(:,:),divE(:)
    complex ( kind = 8 ) ima

    integer *8 count1,count2
    integer *8 ifa_vect,ifb_vect,iflambda,ifrho,ifE,ifcurlE,ifdivE

	ima=(0.0d0,1.0d0)


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
    allocate(b_vect_t(3,nt))

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
      a_vect(1,count1)=(a_u(count1)*u_vect_s(1,count1)+&
      &a_v(count1)*v_vect_s(1,count1))
      a_vect(2,count1)=(a_u(count1)*u_vect_s(2,count1)+&
      &a_v(count1)*v_vect_s(2,count1))
      a_vect(3,count1)=(a_u(count1)*u_vect_s(3,count1)+&
      &a_v(count1)*v_vect_s(3,count1))
    enddo
    do count1=1,ns	
      b_vect(1,count1)=(b_u(count1)*u_vect_s(1,count1)+&
      &b_v(count1)*v_vect_s(1,count1))
      b_vect(2,count1)=(b_u(count1)*u_vect_s(2,count1)+&
      &b_v(count1)*v_vect_s(2,count1))
      b_vect(3,count1)=(b_u(count1)*u_vect_s(3,count1)+&
      &b_v(count1)*v_vect_s(3,count1))
    enddo


    !Computing the full operator
    ifa_vect=1
    ifb_vect=1
    iflambda=1
    ifrho=0
    ifE=1
    ifcurlE=0
    ifdivE=0

	call Vector_Helmholtz_targ2(eps,zk,ns,source,wts,ifa_vect,a_vect,&
	 &ifb_vect,b_vect,iflambda,rho_in,ifrho,rho,n_vect_s,ifE,E,ifcurlE,&
	 &curlE,ifdivE,divE,nt,targets,thresh,ifdir)

!	call Vector_Helmholtz_targ(eps,zk,ns,source,wts,ifa_vect,a_vect,ifb_vect,&
!	 &b_vect,iflambda,rho_in,ifrho,rho,n_vect_s,ifE,E,ifcurlE,curlE,ifdivE,divE,nt,targets)

    do count1=1,nt
      b_vect_t(1,count1)=n_vect_t(2,count1)*E(3,count1)-&
      &n_vect_t(3,count1)*E(2,count1)
      b_vect_t(2,count1)=n_vect_t(3,count1)*E(1,count1)-&
      &n_vect_t(1,count1)*E(3,count1)
      b_vect_t(3,count1)=n_vect_t(1,count1)*E(2,count1)-&
      &n_vect_t(2,count1)*E(1,count1)
    enddo

    do count1=1,nt
      AA_u(count1)=b_vect_t(1,count1)*u_vect_t(1,count1)+&
      &b_vect_t(2,count1)*u_vect_t(2,count1)+b_vect_t(3,count1)*&
      &u_vect_t(3,count1)
      AA_v(count1)=b_vect_t(1,count1)*v_vect_t(1,count1)+&
      &b_vect_t(2,count1)*v_vect_t(2,count1)+b_vect_t(3,count1)*&
      &v_vect_t(3,count1)		
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
    deallocate(b_vect_t)

return
end subroutine CKif_good_FMM


subroutine 	get_RHS_CKif(P0,vf,ns,srcvals,zk,alpha,RHS,normal_H)
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
!    RHS - complex  *16(3*ns)
!      right hand side
!      RHS(1:ns) - first component of -nxE_inc along
!      the srcvals(4:6,i) direction
!      RHS(ns+1:2*ns) - second component of  nxE_inc along
!      the (srcvals(10:12,i) x srcvals(4:6,i)) direction
!

    !List of calling arguments
    integer ( kind = 8 ), intent(in) :: ns
    real ( kind = 8 ), intent(in) :: P0(3)
    complex ( kind = 8 ), intent(in) :: vf(3)
    real ( kind = 8 ), intent(in) :: srcvals(12,ns)
    complex ( kind = 8 ), intent(in) :: zk,alpha
    complex ( kind = 8 ), intent(out) :: RHS(2*ns),normal_H(ns)
	
    !List of local variables
    complex ( kind = 8 ), allocatable :: E(:,:), H(:,:)
    integer *8 count1
    real ( kind = 8 ) ru(3),rv(3),cross_aux(3)

    allocate(E(3,ns), H(3,ns))
    call fieldsED(zk,P0,srcvals,ns,E,H,vf,0)
    call fieldsMD(zk,P0,srcvals,ns,E,H,vf,1)

    do count1=1,ns	
      call orthonormalize(srcvals(4:6,count1),srcvals(10:12,count1),&
      &ru,rv)		
      RHS(count1)=-DOT_PRODUCT(rv,E(:,count1))
      RHS(ns+count1)=DOT_PRODUCT(ru,E(:,count1))	
      normal_H(count1)=DOT_PRODUCT(srcvals(10:12,count1),H(:,count1))
    enddo
return
end subroutine get_RHS_CKif





subroutine em_auCKi_pec_reg_FMM(eps,zk,ns,nt,srcvals,targvals,wts,a_u,&
&a_v,AA_u,AA_v,PHI,thresh,ifdir)
implicit none

!
!  This subroutine computes the far field contribution if 
!  the Cki operator via FMM
!
!  operators:
!
!    nxS_{ik}[J] and -n·curlS_{ik}[J]
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
!
!  output:
!    AA_u - complex *16(nt)
!      first component of nxS_{ik}[J] along
!      the srcvals(4:6,i) direction
!    AA_v - complex *16(nt)
!      second component of nxS_{ik}[J] along
!      the (srcvals(10:12,i) x srcvals(4:6,i)) direction
!    PHI - complex *16(nt)
!      scalar component of the operator:
!      -n·curlS_{ik}[J]
!            

    !List of calling arguments
    real ( kind = 8 ), intent(in) :: eps
    complex ( kind = 8 ), intent(in) :: zk
    integer *8, intent(in) :: ns,nt
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),wts(ns)
    real ( kind = 8 ), intent(in) :: targvals(12,nt)
    complex ( kind = 8 ), intent(in) :: a_u(ns),a_v(ns)
    complex ( kind = 8 ), intent(out) :: AA_u(nt),AA_v(nt),PHI(nt)
    real ( kind = 8 ), intent(in) :: thresh
    integer *8, intent(in) :: ifdir

    !List of local variables
    real ( kind = 8 ), allocatable :: n_vect_s(:,:),u_vect_s(:,:)
    real ( kind = 8 ), allocatable :: v_vect_s(:,:),source(:,:)
    real ( kind = 8 ), allocatable :: n_vect_t(:,:),u_vect_t(:,:)
    real ( kind = 8 ), allocatable :: v_vect_t(:,:),targets(:,:)

    complex ( kind = 8 ), allocatable :: a_vect(:,:),b_vect(:,:)
    complex ( kind = 8 ), allocatable :: b_vect_t(:,:)
    complex ( kind = 8 ), allocatable :: lambda(:),rho(:)
    complex ( kind = 8 ), allocatable :: E(:,:),curlE(:,:),divE(:)
    complex ( kind = 8 ) ima,izk

    integer *8 count1,count2
    integer *8 ifa_vect,ifb_vect,iflambda,ifrho,ifE,ifcurlE,ifdivE

    ima=(0.0d0,1.0d0)
    izk=zk*ima

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
    allocate(b_vect_t(3,nt))

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
      b_vect(1,count1)=(a_u(count1)*u_vect_s(1,count1)+a_v(count1)*&
      &v_vect_s(1,count1))
      b_vect(2,count1)=(a_u(count1)*u_vect_s(2,count1)+a_v(count1)*&
      &v_vect_s(2,count1))
      b_vect(3,count1)=(a_u(count1)*u_vect_s(3,count1)+a_v(count1)*&
      &v_vect_s(3,count1))
    enddo

    !Computing the full operator
    ifa_vect=0
    ifb_vect=1
    iflambda=0
    ifrho=0
    ifE=1
    ifcurlE=1
    ifdivE=0


    call Vector_Helmholtz_targ2(eps,izk,ns,source,wts,ifa_vect,a_vect,&
    &ifb_vect,b_vect,iflambda,lambda,ifrho,rho,n_vect_s,ifE,E,ifcurlE,&
    &curlE,ifdivE,divE,nt,targets,thresh,ifdir)

!	call Vector_Helmholtz_targ(eps,izk,ns,source,wts,ifa_vect,a_vect,ifb_vect,&
!	 &b_vect,iflambda,lambda,ifrho,rho,n_vect_s,ifE,E,ifcurlE,curlE,ifdivE,divE,nt,targets)

    do count1=1,nt
      b_vect_t(1,count1)=n_vect_t(2,count1)*E(3,count1)-&
      &n_vect_t(3,count1)*E(2,count1)
      b_vect_t(2,count1)=n_vect_t(3,count1)*E(1,count1)-&
      &n_vect_t(1,count1)*E(3,count1)
      b_vect_t(3,count1)=n_vect_t(1,count1)*E(2,count1)-&
      &n_vect_t(2,count1)*E(1,count1)
    enddo

    do count1=1,nt
      AA_u(count1)=b_vect_t(1,count1)*u_vect_t(1,count1)+b_vect_t(2,count1)&
      &*u_vect_t(2,count1)+b_vect_t(3,count1)*u_vect_t(3,count1)
      AA_v(count1)=b_vect_t(1,count1)*v_vect_t(1,count1)+b_vect_t(2,count1)&
      &*v_vect_t(2,count1)+b_vect_t(3,count1)*v_vect_t(3,count1)		
      PHI(count1)=-(n_vect_t(1,count1)*curlE(1,count1)+&
      &n_vect_t(2,count1)*curlE(2,count1)&
      &+n_vect_t(3,count1)*curlE(3,count1))
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
    deallocate(b_vect_t)
	
return
end subroutine em_auCKi_pec_reg_FMM




subroutine em_auCKi_pec_vect2rhs_FMM(eps,zpars,ns,nt,srcvals,targvals,&
&wts,a_u,a_v,b_u,b_v,PHI,thresh,ifdir)
implicit none

!
!  This subroutine computes the far field contribution if the Cki 
!  operator via FMM
!
!  operator:
!
!    n·ik·S_{k}[J]+i·alpha·k·n·curlS_{k}[nxS_{ik}[J]]
!
!  input:
!
!    eps - real * 8
!      epsilon for the fmm call
!
!    zpars - complex *16(3)
!      kernel parameters
!      zpars(1) = k 
!      zpars(2) = alpha
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
!    b_u,b_v - complex *16(ns)
!      two components of the tangent induced current nxS_{ik}[J] on 
!      the surface along srcvals(4:6,i) and 
!      (srcvals(10:12,i) x srcvals(4:6,i)) directions
!
!  output:
!    PHI - complex *16(nt)
!      Normal component of the field
!            


    !List of calling arguments
    real ( kind = 8 ), intent(in) :: eps
    complex ( kind = 8 ), intent(in) :: zpars(3)
    integer *8, intent(in) :: ns,nt
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),wts(ns)
    real ( kind = 8 ), intent(in) :: targvals(12,nt)
    complex ( kind = 8 ), intent(in) :: a_u(ns),a_v(ns)
    complex ( kind = 8 ), intent(in) :: b_u(ns),b_v(ns)
    complex ( kind = 8 ), intent(out) :: PHI(nt)
    real ( kind = 8 ), intent(in) :: thresh
    integer *8, intent(in) :: ifdir


    !List of local variables
    real ( kind = 8 ), allocatable :: n_vect_s(:,:),u_vect_s(:,:)
    real ( kind = 8 ), allocatable :: v_vect_s(:,:),source(:,:)
    real ( kind = 8 ), allocatable :: n_vect_t(:,:),u_vect_t(:,:)
    real ( kind = 8 ), allocatable :: v_vect_t(:,:),targets(:,:)

    complex ( kind = 8 ), allocatable :: a_vect(:,:),b_vect(:,:)
    complex ( kind = 8 ), allocatable :: lambda(:)
    complex ( kind = 8 ), allocatable :: E(:,:),curlE(:,:),divE(:),rho(:)
    complex ( kind = 8 ) ima,izk,zk,alpha

    integer *8 count1,count2
    integer *8 ifa_vect,ifb_vect,iflambda,ifrho,ifE,ifcurlE,ifdivE

    ima=(0.0d0,1.0d0)
    zk=zpars(1)
    alpha=zpars(2)
    izk=zk*ima

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
      b_vect(1,count1)=-ima*zk*(a_u(count1)*u_vect_s(1,count1)+&
      &a_v(count1)*v_vect_s(1,count1))
      b_vect(2,count1)=-ima*zk*(a_u(count1)*u_vect_s(2,count1)+&
      &a_v(count1)*v_vect_s(2,count1))
      b_vect(3,count1)=-ima*zk*(a_u(count1)*u_vect_s(3,count1)+&
      &a_v(count1)*v_vect_s(3,count1))
    enddo
    do count1=1,ns
      a_vect(1,count1)=zk*alpha*(b_u(count1)*u_vect_s(1,count1)+&
      &b_v(count1)*v_vect_s(1,count1))
      a_vect(2,count1)=zk*alpha*(b_u(count1)*u_vect_s(2,count1)+&
      &b_v(count1)*v_vect_s(2,count1))
      a_vect(3,count1)=zk*alpha*(b_u(count1)*u_vect_s(3,count1)+&
      &b_v(count1)*v_vect_s(3,count1))
    enddo

!setting flags
    ifa_vect=1
    ifb_vect=1
    iflambda=0
    ifrho=0
    ifE=1
    ifcurlE=0
    ifdivE=0


    call Vector_Helmholtz_targ2(eps,zk,ns,source,wts,ifa_vect,a_vect,&
    &ifb_vect,b_vect,iflambda,lambda,ifrho,rho,n_vect_s,ifE,E,ifcurlE,&
    &curlE,ifdivE,divE,nt,targets,thresh,ifdir)

    do count1=1,nt
      PHI(count1)=(n_vect_t(1,count1)*E(1,count1)+n_vect_t(2,count1)*&
      &E(2,count1)+n_vect_t(3,count1)*E(3,count1))
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
end subroutine em_auCKi_pec_vect2rhs_FMM

    
subroutine test_accuracy_auCKi(eps_FMM,a_vect,a_vect2,sigma,sigma2,zk,&
&alpha,ns,wts,srcvals,P0,vf,Pt)
implicit none

!
!  This function test the accuracy of the solution computed in the 
!  exterior region
!
!  input:
!    eps_FMM - real *8
!      epsilon for the fmm call
!
!    a_vect - complex *16(2*ns)
!      induced current on the surface
!      a_vect(1:ns) - first component of  J along
!      the srcvals(4:6,i) direction
!      a_vect(ns+1:2*ns) - second component of J along
!      the (srcvals(10:12,i) x srcvals(4:6,i)) direction
!	
!    a_vect2 - complex *16(2*ns)
!      induced current on the surface
!      a_vect2(1:ns) - first component of nxS_{ik}[J] along
!      the srcvals(4:6,i) direction
!      a_vect2(ns+1:2*ns) - second component of nxS_{ik}[J] along
!      the (srcvals(10:12,i) x srcvals(4:6,i)) direction
!
!    sigma - complex *16(ns)
!      induced charge on the surface solution of the auxiliary
!      Neumann problem
!	
!    sigma2 - complex *16(ns)
!      value of S_{ik}[sigma] on the surface
!
!    zpars - complex *16(3)
!      kernel parameters
!      zpars(1) = k 
!      zpars(2) = alpha
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
    complex ( kind = 8 ), intent(in) :: zk,alpha
    integer *8, intent(in) :: ns
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),eps_FMM
    real ( kind = 8 ), intent(in) :: wts(ns),P0(3),Pt(3)
    complex ( kind = 8 ), intent(in) :: a_vect(2*ns),vf(3),a_vect2(3*ns)
    complex ( kind = 8 ), intent(in) :: sigma(ns),sigma2(ns)

    !List of local variables
    complex ( kind = 8 ) ima, R1, R2,Et1(3),Et2(3),Ht1(3),Ht2(3)
    complex ( kind = 8 ) E01(3),E02(3),H01(3),H02(3)
    real ( kind = 8 ) ru_s(3),rv_s(3),n_s(3),sour(3),r,dr(3)
    real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3)
    real ( kind = 8 ) error_E,rel_err_E,error_H,rel_err_H
    real ( kind = 8 ) pi
    integer *8 count1

    ima=(0.0d0,1.0d0)
    pi=3.1415926535897932384626433832795028841971d0
!	write (*,*) 'alpha value at accu: ', alpha
    write (*,*) 'AUGMENTED ALGORITHM:'

    call auCKif_FMM_targ(eps_FMM,zk,alpha,ns,srcvals,1,Pt,wts,&
    &a_vect(1:ns),a_vect(ns+1:2*ns),a_vect2(1:ns),a_vect2(ns+1:2*ns),&
    &a_vect2(2*ns+1:3*ns),Et1,Ht1)
!	call CKif_FMM_targ_v2(eps_FMM,zk,alpha,ns,srcvals,1,Pt,wts,sol(1:ns),sol(ns+1:2*ns),Et1,Ht1,nnz,norder,row_ptr,col_ind,CKif_near)
!	write (*,*) Et1
	
    call fieldsED(zk,P0,Pt,1,Et2,Ht2,vf,0)
    call fieldsMD(zk,P0,Pt,1,Et2,Ht2,vf,1)

    write (*,*) 'Errors at the EXTERIOR region:'
    error_E=sqrt(abs(Et1(1)-Et2(1))**2+abs(Et1(2)-Et2(2))**2+&
    &abs(Et1(3)-Et2(3))**2)
!		write (*,*) 'Error E: ', error_E
    write (*,*) 'Relative Error in E: ', error_E/&
    &sqrt(abs(Et2(1))**2+abs(Et2(2))**2+abs(Et2(3))**2)
				
    call h_neumann_Hfield_FMM_targ(eps_FMM,zk,alpha,ns,srcvals,1,Pt,&
    &wts,sigma,sigma2,Ht1)

    error_H=sqrt(abs(Ht1(1)-Ht2(1))**2+abs(Ht1(2)-Ht2(2))**2+&
    &abs(Ht1(3)-Ht2(3))**2)
		
    write (*,*) 'Error H: ', error_H
    write (*,*) 'Relative Error in H: ', error_H/sqrt(abs(Ht2(1))**2+&
    &abs(Ht2(2))**2+abs(Ht2(3))**2)

return
end subroutine test_accuracy_auCKi


subroutine h_neumann_Hfield_FMM_targ(eps,zk,alpha,ns,srcvals,nt,targ,&
&wts,rho_in,sigma2,H)
implicit none

!
!  This funciton computes grad PHI term in the representation of H 
!  using the FMM
!  It doesn't contain near field corrections (it's for debugging 
!  purposes)   
!      
!  Representation:
!
!    H=ik·S_{k}[J]+(...)+i·alpha·k·curlS_{k}[nxS_{ik}[J]]
!    (...)=grad(S_{k}[sigma]+i·alpha·D_{k}[S_{ik}[sigma]]) 
!          <-THIS IS THE TERM COMPUTED HERE
!
!  it requires as input rho_in=sigma and sigma2=S_{ik}[sigma]
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
!    rho_in - complex *16(ns)
!      induced density on the surface
! 
!    sigma2 - complex *16(ns)
!      sigma2=S_{ik}[rho_in]
!
!  output:
!    E - complex  *16(3,nt)
!      value of the Electric field at the target points
!
!    H - complex  *16(3,nt)
!      part of the value of the Magnetic field at the target 
!      points (without the term (...)=grad PHI ) 
!

    !List of calling arguments
    real ( kind = 8 ), intent(in) :: eps
    complex ( kind = 8 ), intent(in) :: zk, alpha
    integer *8, intent(in) :: ns,nt
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),targ(3,nt)
    real ( kind = 8 ), intent(in) :: wts(ns)
    complex ( kind = 8 ), intent(in) :: rho_in(ns),sigma2(ns)
    complex ( kind = 8 ), intent(inout) :: H(3,nt)

    !List of local variables
    real ( kind = 8 ), allocatable :: n_vect(:,:),source(:,:)
    complex ( kind = 8 ), allocatable :: charge(:),dipvec(:,:)
    complex ( kind = 8 ), allocatable :: pot_aux(:),H_aux(:,:)
    complex ( kind = 8 ), allocatable :: pot_aux2(:)
    complex ( kind = 8 ) ima

    integer *8 count1,count2,npols,ntri,ier
    integer *8 ifa_vect,ifb_vect,iflambda,ifrho,ifE,ifcurlE,ifdivE
    real ( kind = 8 ) pi
    pi=3.1415926535897932384626433832795028841971d0

    ima=(0.0d0,1.0d0)
    allocate(charge(ns))
    allocate(pot_aux(ns))
    allocate(dipvec(3,ns))
    allocate(H_aux(3,nt))
    allocate(pot_aux2(nt))
    allocate(n_vect(3,ns))
    allocate(source(3,ns))

    do count1=1,ns
      n_vect(:,count1)=srcvals(10:12,count1)
      source(:,count1)=srcvals(1:3,count1)
    enddo

    do count1=1,ns
      charge(count1)=rho_in(count1)*wts(count1)
      dipvec(1,count1)=n_vect(1,count1)*sigma2(count1)*&
      &wts(count1)*ima*alpha
      dipvec(2,count1)=n_vect(2,count1)*sigma2(count1)*&
      &wts(count1)*ima*alpha
      dipvec(3,count1)=n_vect(3,count1)*sigma2(count1)*&
      &wts(count1)*ima*alpha
    enddo


    !Computing the full operator
    call hfmm3d_t_cd_g(eps,zk,ns,source,charge,dipvec,nt,targ,pot_aux2,&
    &H_aux,ier)
!    call hfmm3d_t_cd_p(eps,zk,ns,source,charge,dipvec,nt,targ,pot)
			

    do count1=1,nt
      H_aux(1,count1)=H_aux(1,count1)/(4.0d0*pi)
      H_aux(2,count1)=H_aux(2,count1)/(4.0d0*pi)
      H_aux(3,count1)=H_aux(3,count1)/(4.0d0*pi)
    enddo
    do count1=1,nt
      H(1,count1)=H(1,count1)+H_aux(1,count1)
      H(2,count1)=H(2,count1)+H_aux(2,count1)
      H(3,count1)=H(3,count1)+H_aux(3,count1)
    enddo

    deallocate(charge)
    deallocate(dipvec)
    deallocate(n_vect)
    deallocate(source)

return
end subroutine h_neumann_Hfield_FMM_targ


subroutine auCKif_FMM_targ(eps,zk,alpha,ns,srcvals,nt,targ,wts,a_u,a_v,&
&b_u,b_v,sigma_s,E,H)
implicit none

!
!  This funciton computes the fields E,H at a given point using the FMM
!  It doesn't contain near field corrections (it's for debugging 
!  purposes)   
!
!  Representation:
!
!    E=curlS_{k}[J]+i·alpha·curlcurlS_{k}[nxS_{ik}[J]]
!    H=ik·S_{k}[J]+(...)+i·alpha·k·curlS_{k}[nxS_{ik}[J]]
!    (...)=grad(S_{k}[sigma]+i·alpha·D_{k}[S_{ik}[sigma]])
!
!  it requires as input J, nxS_{ik}[J] and n·curlS_{ik}[J]
!  the gradPHI term is not computed here, it will be added later
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
!    b_u,b_v - complex *16(ns)
!      two components of nxS_{ik}[J] on the surface
!      along srcvals(4:6,i) and (srcvals(10:12,i) x srcvals(4:6,i))
!      directions
!      rsigma_s - complex *16(ns)
!      the scalar quantity -n·curlS_{ik}[J] on the surface
!
!  output:
!    E - complex  *16(3,nt)
!      value of the Electric field at the target points
!
!    H - complex  *16(3,nt)
!      part of the value of the Magnetic field at the target points
!      (without the term (...)=grad PHI ) 
!

    !List of calling arguments
    real ( kind = 8 ), intent(in) :: eps
    complex ( kind = 8 ), intent(in) :: zk,alpha
    integer *8, intent(in) :: ns,nt
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),targ(3,nt)
    real ( kind = 8 ), intent(in) :: wts(ns)
    complex ( kind = 8 ), intent(in) :: a_u(ns),a_v(ns)
    complex ( kind = 8 ), intent(in) :: b_u(ns),b_v(ns),sigma_s(ns)
    complex ( kind = 8 ), intent(out) :: E(3,nt),H(3,nt)

    !List of local variables
    real ( kind = 8 ), allocatable :: n_vect(:,:),u_vect(:,:)	
    real ( kind = 8 ), allocatable :: v_vect(:,:),source(:,:)
    complex ( kind = 8 ), allocatable :: a_vect(:,:),b_vect(:,:)
    complex ( kind = 8 ), allocatable :: lambda(:),rho(:)
    complex ( kind = 8 ), allocatable :: curlE(:,:),divE(:)
    complex ( kind = 8 ) ima

    integer *8 count1,count2
    integer *8 ifa_vect,ifb_vect,iflambda,ifrho,ifE,ifcurlE,ifdivE

    ima=(0.0d0,1.0d0)

    allocate(a_vect(3,ns))
    allocate(b_vect(3,ns))
    allocate(lambda(ns))
    allocate(rho(ns))
!	allocate(b_u(ns))
!	allocate(b_v(ns))
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
    call orthonormalize_all(srcvals(4:6,:),srcvals(10:12,:),u_vect,&
	&v_vect,ns)

    do count1=1,ns
      a_vect(1,count1)=(a_u(count1)*u_vect(1,count1)+a_v(count1)*&
	  &v_vect(1,count1))
      a_vect(2,count1)=(a_u(count1)*u_vect(2,count1)+a_v(count1)*&
	  &v_vect(2,count1))
      a_vect(3,count1)=(a_u(count1)*u_vect(3,count1)+a_v(count1)*&
	  &v_vect(3,count1))
				
      b_vect(1,count1)=(b_u(count1)*u_vect(1,count1)+b_v(count1)*&
	  &v_vect(1,count1))*zk**2*ima*alpha
      b_vect(2,count1)=(b_u(count1)*u_vect(2,count1)+b_v(count1)*&
	  &v_vect(2,count1))*zk**2*ima*alpha
      b_vect(3,count1)=(b_u(count1)*u_vect(3,count1)+b_v(count1)*&
	  &v_vect(3,count1))*zk**2*ima*alpha
		
      lambda(count1)=sigma_s(count1)*ima*alpha
    enddo

    !Computing the full operator
    ifa_vect=1
    ifb_vect=1
    iflambda=1
    ifrho=0
    ifE=1
    ifcurlE=0   !Cambia esto para el cálculo de H estable!!
    ifdivE=0

    call Vector_Helmholtz_targ(eps,zk,ns,source,wts,ifa_vect,a_vect,&
    &ifb_vect,b_vect,iflambda,lambda,ifrho,rho,n_vect,ifE,E,ifcurlE,&
    &H,ifdivE,divE,nt,targ)

    do count1=1,ns
      a_vect(1,count1)=(b_u(count1)*u_vect(1,count1)+b_v(count1)*&
      &v_vect(1,count1))*alpha*zk
      a_vect(2,count1)=(b_u(count1)*u_vect(2,count1)+b_v(count1)*&
      &v_vect(2,count1))*alpha*zk
      a_vect(3,count1)=(b_u(count1)*u_vect(3,count1)+b_v(count1)*&
      &v_vect(3,count1))*alpha*zk
				
      b_vect(1,count1)=(a_u(count1)*u_vect(1,count1)+a_v(count1)*&
      &v_vect(1,count1))*(-ima*zk)
      b_vect(2,count1)=(a_u(count1)*u_vect(2,count1)+a_v(count1)*&
      &v_vect(2,count1))*(-ima*zk)
      b_vect(3,count1)=(a_u(count1)*u_vect(3,count1)+a_v(count1)*&
      &v_vect(3,count1))*(-ima*zk)
    enddo

!!!This part is for the magnetic field

    !Computing the full operator
    ifa_vect=1
    ifb_vect=1
    iflambda=0
    ifrho=0
    ifE=1
    ifcurlE=0
    ifdivE=0

    call Vector_Helmholtz_targ(eps,zk,ns,source,wts,ifa_vect,a_vect,&
    &ifb_vect,b_vect,iflambda,lambda,ifrho,rho,n_vect,ifE,H,ifcurlE,&
    &curlE,ifdivE,divE,nt,targ) !! WATCH OUT!!, E and H are switched
!	from the usual place!!!!


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
end subroutine auCKif_FMM_targ




subroutine test_accuracy_CKi(eps_FMM,a_vect,a_vect2,zk,alpha,ns,wts,&
&srcvals,P0,vf,Pt)
implicit none

!
!  This function test the accuracy of the solution computed in the 
!  exterior region
!  Same as test_accuracy_auCKi but computing H as curlE/(ik) 
!  (in a unstable way, tol compare with the augmented version here)
!
!  input:
!    eps_FMM - real *8
!      epsilon for the fmm call
!
!    a_vect - complex *16(2*ns)
!      induced current on the surface
!      a_vect(1:ns) - first component of  J along
!      the srcvals(4:6,i) direction
!      a_vect(ns+1:2*ns) - second component of J along
!      the (srcvals(10:12,i) x srcvals(4:6,i)) direction
!	
!    a_vect2 - complex *16(3*ns)
!      induced current on the surface
!      a_vect2(1:ns) - first component of nxS_{ik}[J] along
!      the srcvals(4:6,i) direction
!      a_vect2(ns+1:2*ns) - second component of nxS_{ik}[J] along
!      the (srcvals(10:12,i) x srcvals(4:6,i)) direction
!      a_vect2(2*ns+1:3*ns) - normal component, -ncurlS_{ik}[J] along
!      the direction srcvals(10:12,i)
!
!    zpars - complex *16(3)
!      kernel parameters
!      zpars(1) = k 
!      zpars(2) = alpha
!      zpars(3) = - (not used)
!
!    ns - integer *8
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
    complex ( kind = 8 ), intent(in) :: zk,alpha
    integer *8, intent(in) :: ns
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),eps_FMM
    real ( kind = 8 ), intent(in) :: wts(ns),P0(3),Pt(3)
    complex ( kind = 8 ), intent(in) :: a_vect(2*ns),vf(3),a_vect2(3*ns)


    !List of local variables
    complex ( kind = 8 ) ima, R1, R2,Et1(3),Et2(3),Ht1(3),Ht2(3)
    complex ( kind = 8 ) E01(3),E02(3),H01(3),H02(3)
    real ( kind = 8 ) ru_s(3),rv_s(3),n_s(3),sour(3),r,dr(3)
    real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3)
    real ( kind = 8 ) error_E,rel_err_E,error_H,rel_err_H
    real ( kind = 8 ) pi
    integer *8 count1

    ima=(0.0d0,1.0d0)
    pi=3.1415926535897932384626433832795028841971d0
!		write (*,*) 'alpha value at accu: ', alpha
    call CKif_FMM_targ(eps_FMM,zk,alpha,ns,srcvals,1,Pt,wts,&
    &a_vect(1:ns),a_vect(ns+1:2*ns),a_vect2(1:ns),a_vect2(ns+1:2*ns),&
    &a_vect2(2*ns+1:3*ns),Et1,Ht1)
!	call CKif_FMM_targ_v2(eps_FMM,zk,alpha,ns,srcvals,1,Pt,wts,sol(1:ns),sol(ns+1:2*ns),Et1,Ht1,nnz,norder,row_ptr,col_ind,CKif_near)
!	write (*,*) Et1
	
    call fieldsED(zk,P0,Pt,1,Et2,Ht2,vf,0)
    call fieldsMD(zk,P0,Pt,1,Et2,Ht2,vf,1)

    write (*,*) 'Errors at the EXTERIOR region:'
    error_E=sqrt(abs(Et1(1)-Et2(1))**2+abs(Et1(2)-Et2(2))**2+&
    &abs(Et1(3)-Et2(3))**2)
!		write (*,*) 'Error E: ', error_E
    write (*,*) 'Relative Error in E: ', error_E/&
    &sqrt(abs(Et2(1))**2+abs(Et2(2))**2+abs(Et2(3))**2)
!		write (*,*) Et1
!		write (*,*) Et2
		
!call h_neumann_Hfield_FMM_targ(eps_FMM,zk,alpha,ns,srcvals,1,Pt,wts,sigma,sigma2)

    error_H=sqrt(abs(Ht1(1)-Ht2(1))**2+abs(Ht1(2)-Ht2(2))**2+&
    &abs(Ht1(3)-Ht2(3))**2)
    write (*,*) 'Error H: ', error_H
    write (*,*) 'Relative Error in H: ', error_H/&
    &sqrt(abs(Ht2(1))**2+abs(Ht2(2))**2+abs(Ht2(3))**2)

return
end subroutine test_accuracy_CKi



subroutine CKif_FMM_targ(eps,zk,alpha,ns,srcvals,nt,targ,wts,a_u,a_v,&
&b_u,b_v,sigma_s,E,H)
implicit none

!
!  This funciton computes the fields E,H at a given point using the FMM
!  It doesn't contain near field corrections (it's for debugging 
!  purposes)   
!
!  Representation:
!
!    E=curlS_{k}[J]+i·alpha·curlcurlS_{k}[nxS_{ik}[J]]
!    H=curl E /(ik)   -> this is known to be unstable in low frequency
!    You should use instead auCKif_FMM_targ together 
!    with h_neumann_Hfield_FMM_targ
!
!  it requires as input J, nxS_{ik}[J] and n·curlS_{ik}[J]
!  the gradPHI term is not computed here, it will be added later
!
!  input:
!    eps - real *8
!      epsilon for the fmm call
!
!    zk - complex *16
!      Helmholtz parameter 
!
!    ns - integer *8
!      number of sources
!   
!    nt - integer *8
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
!    b_u,b_v - complex *16(ns)
!      two components of nxS_{ik}[J] on the surface
!      along srcvals(4:6,i) and (srcvals(10:12,i) x srcvals(4:6,i))
!      directions
!
!    sigma_s - complex *16(ns)
!      the scalar quantity -n·curlS_{ik}[J] on the surface
!
!  output:
!    E - complex  *16(3,nt)
!      value of the Electric field at the target points
!
!    H - complex  *16(3,nt)
!      part of the value of the Magnetic field at the target points
!

    !List of calling arguments
    real ( kind = 8 ), intent(in) :: eps
    complex ( kind = 8 ), intent(in) :: zk,alpha
    integer *8, intent(in) :: ns,nt
    real ( kind = 8 ), intent(in) :: srcvals(12,ns),targ(3,nt)
    real ( kind = 8 ), intent(in) :: wts(ns)
    complex ( kind = 8 ), intent(in) :: a_u(ns),a_v(ns),b_u(ns),b_v(ns)
    complex ( kind = 8 ), intent(in) :: sigma_s(ns)
    complex ( kind = 8 ), intent(out) :: E(3,nt),H(3,nt)

    !List of local variables
    real ( kind = 8 ), allocatable :: n_vect(:,:),u_vect(:,:)
    real ( kind = 8 ), allocatable :: v_vect(:,:),source(:,:)
    complex ( kind = 8 ), allocatable :: a_vect(:,:),b_vect(:,:)
    complex ( kind = 8 ), allocatable :: lambda(:),rho(:)
    complex ( kind = 8 ), allocatable :: curlE(:,:),divE(:)
    complex ( kind = 8 ) ima

    integer *8 count1,count2
    integer *8 ifa_vect,ifb_vect,iflambda,ifrho,ifE,ifcurlE,ifdivE

    ima=(0.0d0,1.0d0)

    allocate(a_vect(3,ns))
    allocate(b_vect(3,ns))
    allocate(lambda(ns))
    allocate(rho(ns))
!	allocate(b_u(ns))
!	allocate(b_v(ns))
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
    call orthonormalize_all(srcvals(4:6,:),srcvals(10:12,:),u_vect,&
    &v_vect,ns)

    do count1=1,ns
      a_vect(1,count1)=(a_u(count1)*u_vect(1,count1)+a_v(count1)*&
      &v_vect(1,count1))
      a_vect(2,count1)=(a_u(count1)*u_vect(2,count1)+a_v(count1)*&
      &v_vect(2,count1))
      a_vect(3,count1)=(a_u(count1)*u_vect(3,count1)+a_v(count1)*&
      &v_vect(3,count1))

      b_vect(1,count1)=(b_u(count1)*u_vect(1,count1)+b_v(count1)*&
      &v_vect(1,count1))*zk**2*ima*alpha
      b_vect(2,count1)=(b_u(count1)*u_vect(2,count1)+b_v(count1)*&
      &v_vect(2,count1))*zk**2*ima*alpha
      b_vect(3,count1)=(b_u(count1)*u_vect(3,count1)+b_v(count1)*&
      &v_vect(3,count1))*zk**2*ima*alpha

      lambda(count1)=sigma_s(count1)*ima*alpha
    enddo

    !Computing the full operator
    ifa_vect=1
    ifb_vect=1
    iflambda=1
    ifrho=0
    ifE=1
    ifcurlE=1   !Cambia esto para el cálculo de H estable!!
    ifdivE=0

    call Vector_Helmholtz_targ(eps,zk,ns,source,wts,ifa_vect,a_vect,&
    &ifb_vect,b_vect,iflambda,lambda,ifrho,rho,n_vect,ifE,E,ifcurlE,H,&
    &ifdivE,divE,nt,targ)

    do count1=1,nt
      H(1,count1)=H(1,count1)/(ima*zk)
      H(2,count1)=H(2,count1)/(ima*zk)
      H(3,count1)=H(3,count1)/(ima*zk)
    enddo

    deallocate(a_vect)
    deallocate(b_vect)
    deallocate(lambda)
    deallocate(curlE)
    deallocate(rho)
!	deallocate(b_u)
!	deallocate(b_v)
	
    deallocate(u_vect)
    deallocate(v_vect)
    deallocate(n_vect)
    deallocate(source)
    deallocate(divE)

return
end subroutine CKif_FMM_targ


      subroutine lpcomp_CKi2RHS_addsub(npatches,norders,ixyzs,&
     &iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
     &eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,sigma,novers,&
     &nptso,ixyzso,srcover,whtsover,pot,sigma2,scalar_RHS,&
	 &wnear)
!
!  This subroutine evaluates the layer potential for
!  the boundary integral equation:
!
!    J/2+nxcurlS_{k}[J]+i·alpha·nxcurlcurlS_{k}[nxS_{ik}[J]]= -nxE_inc
!
!  where the near field is precomputed and stored in the row sparse 
!  compressed format.
!
!  it also provides the calculation:
!
!    scalar_RHS=n·ik·S_{k}[J]+i·alpha·k·n·curlS_{k}[nxS_{ik}[J]] 
!
!  Note the 4\pi scaling is NOT included as the FMM output was scaled
!  appropriately
!
!  Note: the identity J/2 is not included as the gmres takes care 
!  of that
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
!        srcvals(1:3,i) - xyz info
!        srcvals(4:6,i) - dxyz/du info
!        srcvals(7:9,i) - dxyz/dv info
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
!      location in wnear array where quadrature for col_ind(i) starts
!
!    nquad - integer *8
!      number of entries in wnear
!
!    wnear  - complex *16(nquad)
!      the desired near field quadrature
!
!    sigma - complex *16(2*ns)
!      induced charge and current on the surface
!      sol(1:ns) - first component of  J along
!        the srcvals(4:6,i) direction
!      sol(ns+1:2*ns) - second component of J along
!        the (srcvals(10:12,i) x srcvals(4:6,i)) direction
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
      real *8 targs(ndtarg,ntarg)
      complex *16 zpars(3),alpha,zk
      integer *8 nnz,row_ptr(ntarg+1),col_ind(nnz),nquad
      integer *8 iquad(nnz+1)
      complex *16 sigma(2*npts),sigma2(3*npts)
	  
      complex *16 wnear(12*nquad)

      integer *8 novers(npatches+1)
      integer *8 nover,npolso,nptso
      real *8 srcover(12,nptso),whtsover(nptso)
      complex *16 pot(2*ntarg),scalar_RHS(npts)
      complex *16, allocatable :: potsort(:)

      real *8, allocatable :: sources(:,:),targvals(:,:),wtmp2(:)
      complex *16, allocatable :: charges(:),dipvec(:,:),sigmaover(:)
	  complex *16, allocatable :: sigma_aux(:),sigmaover_aux(:)
	  complex *16, allocatable :: pot_aux(:)

      integer *8 ns,nt
      complex *16 beta
      integer *8 ifcharge,ifdipole
      integer *8 ifpgh,ifpghtarg
      complex *16 tmp(10),val,E(3)

      real *8 xmin,xmax,ymin,ymax,zmin,zmax,sizey,sizez,boxsize

      integer *8 i,j,jpatch,jquadstart,jstart

      integer *8 ifaddsub,ifdir

      integer *8 ntj
      
      complex *16 zdotu,pottmp
      complex *16, allocatable :: dtmp2(:,:),ctmp2_b_u(:),ctmp2_b_v(:)
	  complex *16, allocatable :: ctmp2_a_u(:),ctmp2_a_v(:),ctmp2_b_s(:)
      complex *16, allocatable :: ctmp2_u(:),ctmp2_v(:),ctmp2_s(:)
      real *8 radexp,epsfmm

      integer *8 ipars(2)
      real *8 dpars(1),timeinfo(10),t1,t2,omp_get_wtime

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      real *8 thresh,ra
      real *8 rr,rmin
      integer *8 nss,ii,l,npover

      integer *8 nd,ntarg0

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
	  allocate(sigma_aux(3*ntarg),pot_aux(3*ntarg))

! 
!       oversample density
!
		zk=zpars(1)
		alpha=zpars(2)

      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
     &npts,sigma(1:npts),novers,ixyzso,ns,sigmaover(1:ns))
      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
     &npts,sigma(npts+1:2*npts),novers,ixyzso,ns,sigmaover(ns+1:2*ns))

      ra = 0
      
      !
      !     FMM call to compute nxSkia and -n·curlSkia
      !
	  call get_fmm_thresh(12,ns,srcover,ndtarg,ntarg,targs,thresh)
	  
       ifdir=0
	   
	  call em_auCKi_pec_reg_FMM(eps,zpars(1),ns,npts,srcover,targs,&
      &whtsover,sigmaover(1:ns),sigmaover(ns+1:2*ns),pot_aux(1:npts),&
      &pot_aux(npts+1:2*npts),pot_aux(2*npts+1:3*npts),thresh,ifdir)
	  
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
            pot_aux(i+2*npts)=pot_aux(i+2*npts)+&
            &wnear(6*nquad+jquadstart+l-1)*sigma(jstart+l-1)
            pot_aux(i+2*npts)=pot_aux(i+2*npts)+&
            &wnear(7*nquad+jquadstart+l-1)*sigma(jstart+l-1+npts)		  
		  
            call orthonormalize(srcvals(4:6,jstart+l-1),&
            &srcvals(10:12,jstart+l-1),u_vect,v_vect)
            sigma_x=(sigma(jstart+l-1)*u_vect(1)+sigma(jstart+l-1+npts)*&
            &v_vect(1))
            sigma_y=(sigma(jstart+l-1)*u_vect(2)+sigma(jstart+l-1+npts)*&
            &v_vect(2))
            sigma_z=(sigma(jstart+l-1)*u_vect(3)+sigma(jstart+l-1+npts)*&
            &v_vect(3))

            pot_x=wnear(10*nquad+jquadstart+l-1)*sigma_x
            pot_y=wnear(10*nquad+jquadstart+l-1)*sigma_y
            pot_z=wnear(10*nquad+jquadstart+l-1)*sigma_z
		  
            n_vect=srcvals(10:12,i)
            pot2_x=n_vect(2)*pot_z-n_vect(3)*pot_y
            pot2_y=n_vect(3)*pot_x-n_vect(1)*pot_z
            pot2_z=n_vect(1)*pot_y-n_vect(2)*pot_x
            call orthonormalize(srcvals(4:6,i),srcvals(10:12,i),&
            &u_vect,v_vect)

            pot_aux(i) = pot_aux(i) + pot2_x*u_vect(1)+pot2_y*u_vect(2)+&
            &pot2_z*u_vect(3)
            pot_aux(i+npts) = pot_aux(i+npts) + pot2_x*v_vect(1)+pot2_y*&
            &v_vect(2)+pot2_z*v_vect(3)
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
        allocate(srctmp2(12,nss),ctmp2_a_u(nss))
        allocate(ctmp2_a_v(nss),wtmp2(nss))

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
			wtmp2(ii)=whtsover(jstart+l)
          enddo
        enddo

        call em_auCKi_pec_reg_FMM(eps,zpars(1),nss,ntarg0,srctmp2,&
        &targs(:,i),wtmp2,ctmp2_a_u,ctmp2_a_v,E(1),E(2),E(3),thresh,&
        &ifdir)
		
	    pot_aux(i) = pot_aux(i) - E(1)
	    pot_aux(i+ntarg) = pot_aux(i+ntarg) - E(2)
	    pot_aux(i+2*ntarg) = pot_aux(i+2*ntarg) - E(3)

        deallocate(srctmp2,ctmp2_a_u,ctmp2_a_v,wtmp2)
      enddo
	  	  

      do i=1,ntarg
	    sigma_aux(i)=pot_aux(i)
	    sigma_aux(i+ntarg)=pot_aux(i+ntarg)
	    sigma_aux(i+2*ntarg)=pot_aux(i+2*ntarg)

	    pot_aux(i)=0.0d0
	    pot_aux(i+ntarg)=0.0d0
	    pot_aux(i+2*ntarg)=0.0d0

	    sigma2(i)=sigma_aux(i)
	    sigma2(i+ntarg)=sigma_aux(i+ntarg)
	    sigma2(i+2*ntarg)=sigma_aux(i+2*ntarg)
	  
	    sigma_aux(i)=ima*alpha*zk**2*sigma_aux(i)
	    sigma_aux(i+ntarg)=ima*alpha*zk**2*sigma_aux(i+ntarg)
	    sigma_aux(i+2*ntarg)=ima*alpha*sigma_aux(i+2*ntarg)
      enddo


      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
     &npts,sigma_aux(1:npts),novers,ixyzso,ns,sigmaover_aux(1:ns))
      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
     &npts,sigma_aux(npts+1:2*npts),novers,ixyzso,ns,&
     &sigmaover_aux(ns+1:2*ns))
      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
     &npts,sigma_aux(2*npts+1:3*npts),novers,ixyzso,ns,&
     &sigmaover_aux(2*ns+1:3*ns))

		! FMM call, Here I compute nxSkb
		ifdir=0
	  call CKif_good_FMM(eps,zpars(1),ns,npts,srcover,targs,whtsover,&
     &sigmaover(1:ns),sigmaover(ns+1:2*ns),sigmaover_aux(1:ns),&
     &sigmaover_aux(ns+1:2*ns),sigmaover_aux(2*ns+1:3*ns),&
     &pot_aux(1:npts),pot_aux(npts+1:2*npts),thresh,ifdir)
	  

      do i=1,ntarg

        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch)
          do l=1,npols		  
            call orthonormalize(srcvals(4:6,jstart+l-1),&
            &srcvals(10:12,jstart+l-1),u_vect,v_vect)
            sigma_x=(sigma_aux(jstart+l-1)*u_vect(1)+&
            &sigma_aux(jstart+l-1+npts)*v_vect(1))
            sigma_y=(sigma_aux(jstart+l-1)*u_vect(2)+&
            &sigma_aux(jstart+l-1+npts)*v_vect(2))
            sigma_z=(sigma_aux(jstart+l-1)*u_vect(3)+&
            &sigma_aux(jstart+l-1+npts)*v_vect(3))

            pot_x=wnear(11*nquad+jquadstart+l-1)*sigma_x
            pot_y=wnear(11*nquad+jquadstart+l-1)*sigma_y
            pot_z=wnear(11*nquad+jquadstart+l-1)*sigma_z
		  
            n_vect=srcvals(10:12,i)
            pot2_x=n_vect(2)*pot_z-n_vect(3)*pot_y
            pot2_y=n_vect(3)*pot_x-n_vect(1)*pot_z
            pot2_z=n_vect(1)*pot_y-n_vect(2)*pot_x
            call orthonormalize(srcvals(4:6,i),srcvals(10:12,i),&
            &u_vect,v_vect)

            pot_aux(i) = pot_aux(i) +&
            & pot2_x*u_vect(1)+pot2_y*u_vect(2)+pot2_z*u_vect(3)
            pot_aux(i+npts) = pot_aux(i+npts) +&
            & pot2_x*v_vect(1)+pot2_y*v_vect(2)+pot2_z*v_vect(3)
			
            pot_aux(i) = pot_aux(i) +&
            & wnear(4*nquad+jquadstart+l-1)*sigma_aux(jstart+l-1+2*npts)
            pot_aux(i+npts) = pot_aux(i+npts) +&
            & wnear(5*nquad+jquadstart+l-1)*sigma_aux(jstart+l-1+2*npts)

            pot_aux(i) = pot_aux(i) + &
            &wnear(jquadstart+l-1)*sigma(jstart+l-1)
            pot_aux(i) = pot_aux(i) +&
            & wnear(1*nquad+jquadstart+l-1)*sigma(jstart+l-1+npts)
            pot_aux(i+npts) = pot_aux(i+npts) +&
            & wnear(2*nquad+jquadstart+l-1)*sigma(jstart+l-1)
            pot_aux(i+npts) = pot_aux(i+npts) +&
            & wnear(3*nquad+jquadstart+l-1)*sigma(jstart+l-1+npts)


          enddo
        enddo
      enddo

	ifdir=1

	  do i=1,ntarg
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          nss = nss + ixyzso(jpatch+1)-ixyzso(jpatch)
        enddo
        allocate(srctmp2(12,nss),ctmp2_b_u(nss),ctmp2_b_v(nss))
        allocate(ctmp2_b_s(nss),ctmp2_a_u(nss),ctmp2_a_v(nss))
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
			ctmp2_b_u(ii)=sigmaover_aux(jstart+l)
			ctmp2_b_v(ii)=sigmaover_aux(jstart+l+ns)
			ctmp2_b_s(ii)=sigmaover_aux(jstart+l+2*ns)
			
			
			ctmp2_a_u(ii)=sigmaover(jstart+l)
			ctmp2_a_v(ii)=sigmaover(jstart+l+ns)

			wtmp2(ii)=whtsover(jstart+l)
          enddo
        enddo
		
        call CKif_good_FMM(eps,zpars(1),nss,ntarg0,srctmp2,targs(:,i),&
        &wtmp2,ctmp2_a_u,ctmp2_a_v,ctmp2_b_u,ctmp2_b_v,ctmp2_b_s,&
        &E(1),E(2),thresh,ifdir)
	    pot_aux(i) = pot_aux(i) - E(1)
	    pot_aux(i+ntarg) = pot_aux(i+ntarg) - E(2)

        deallocate(srctmp2,ctmp2_a_u,ctmp2_a_v)
        deallocate(ctmp2_b_s,ctmp2_b_u,ctmp2_b_v,wtmp2)
      enddo
	  
	  do i=1,npts
	    pot(i)=pot_aux(i)
	    pot(i+npts)=pot_aux(i+npts)
	  enddo


!this is the scalar part scalar_RHS:
	do i=1,3*ntarg
	 sigma_aux(i)=sigma2(i)
	  pot_aux(i)=0.0d0
	enddo

      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
     &npts,sigma_aux(1:npts),novers,ixyzso,ns,sigmaover_aux(1:ns))
      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
     &npts,sigma_aux(npts+1:2*npts),novers,ixyzso,ns,&
     &sigmaover_aux(ns+1:2*ns))
      call oversample_fun_surf(2,npatches,norders,ixyzs,iptype,& 
     &npts,sigma_aux(2*npts+1:3*npts),novers,ixyzso,ns,&
     &sigmaover_aux(2*ns+1:3*ns))




      ifdir=0

	  call em_auCKi_pec_vect2rhs_FMM(eps,zpars,ns,npts,srcover,&
      &targs,whtsover,sigmaover(1:ns),sigmaover(ns+1:2*ns),&
      &sigmaover_aux(1:ns),sigmaover_aux(ns+1:2*ns),pot_aux(1:npts),&
      &thresh,ifdir)
	  
!     write (*,*) sigmaover(1:5)
!	  write (*,*) sigmaover(ns+1:ns+5)

!      write (*,*) sigmaover_aux(1:5)
!	  write (*,*) sigmaover_aux(ns+1:ns+5)

!write (*,*) pot_aux(1:10)
!stop

      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch)
          do l=1,npols
            pot_aux(i)=pot_aux(i)+wnear(8*nquad+jquadstart+l-1)*&
            &sigma_aux(jstart+l-1)*zpars(1)*zpars(2)
            pot_aux(i)=pot_aux(i)+wnear(9*nquad+jquadstart+l-1)*&
            &sigma_aux(jstart+l-1+npts)*zpars(1)*zpars(2)		  
		  
            call orthonormalize(srcvals(4:6,jstart+l-1),srcvals(10:12,&
            &jstart+l-1),u_vect,v_vect)
            sigma_x=(sigma(jstart+l-1)*u_vect(1)+sigma(jstart+l-1+npts)*&
            &v_vect(1))
            sigma_y=(sigma(jstart+l-1)*u_vect(2)+sigma(jstart+l-1+npts)*&
            &v_vect(2))
            sigma_z=(sigma(jstart+l-1)*u_vect(3)+sigma(jstart+l-1+npts)*&
            &v_vect(3))

            pot_x=wnear(11*nquad+jquadstart+l-1)*sigma_x
            pot_y=wnear(11*nquad+jquadstart+l-1)*sigma_y
            pot_z=wnear(11*nquad+jquadstart+l-1)*sigma_z
		  
            n_vect=srcvals(10:12,i)
            pot_aux(i)=pot_aux(i)+(n_vect(1)*pot_x+n_vect(2)*pot_y+&
            &n_vect(3)*pot_z)*(-ima*zpars(1))
          enddo
        enddo
      enddo
ifdir=1
	  ipars(1)=3
	  ipars(2)=1
	  do i=1,ntarg
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          nss = nss + ixyzso(jpatch+1)-ixyzso(jpatch)
        enddo
        allocate(srctmp2(12,nss),ctmp2_a_u(nss),ctmp2_a_v(nss))
        allocate(ctmp2_b_u(nss),ctmp2_b_v(nss),wtmp2(nss))

        rmin = 1.0d6
        ii = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          jstart = ixyzso(jpatch)-1
          npover = ixyzso(jpatch+1)-ixyzso(jpatch)
          do l=1,npover
			ii = ii+1
			srctmp2(:,ii) = srcover(:,jstart+l)
			ctmp2_b_u(ii)=sigmaover_aux(jstart+l)
			ctmp2_b_v(ii)=sigmaover_aux(jstart+l+ns)
			ctmp2_a_u(ii)=sigmaover(jstart+l)
			ctmp2_a_v(ii)=sigmaover(jstart+l+ns)			
			wtmp2(ii)=whtsover(jstart+l)
          enddo
        enddo

	  call em_auCKi_pec_vect2rhs_FMM(eps,zpars,nss,ntarg0,srctmp2,&
      &targs(:,i),wtmp2,ctmp2_a_u,ctmp2_a_v,ctmp2_b_u,ctmp2_b_v,&
	  &E(1),thresh,ifdir)

	    pot_aux(i) = pot_aux(i) - E(1)


        deallocate(srctmp2,ctmp2_a_u,ctmp2_a_v,wtmp2)
		deallocate(ctmp2_b_u,ctmp2_b_v)
      enddo
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
        deallocate(srctmp2,ctmp2_u,ctmp2_v,wtmp2)
      enddo

		do i=1,npts
			scalar_RHS(i)=pot_aux(i)
		enddo

      return
      end subroutine lpcomp_CKi2RHS_addsub





subroutine em_auCKi_pec(srcinfo,ndt,targinfo,ndd,dpars,&
&ndz,zpars,ndi,ipars,E_val)
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
    real ( kind = 8 ) c_aux,c_aux_A,c_aux_B,c_aux_C,c_aux_D
    real ( kind = 8 ) xprod_aux3(3),xprod_aux4(3)	
    real ( kind = 8 ) xprod_aux1(3),xprod_aux2(3)
    complex ( kind = 8 ) R1_0,R1_1,ima,my_exp_0,my_exp_1,zk
    complex ( kind = 8 ) alpha,zk0,zk1,nxSk0b(2,2),nxSk1b(2,2)
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
    zk0=zk
    zk1=ima*zk

    call orthonormalize(srcinfo(4:6),n_s,ru_s,rv_s)
    call orthonormalize(targinfo(4:6),n_t,ru_t,rv_t)


    if (ipars(1).eq.1) then
      if (ipars(2).eq.1) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
	  	call my_cross_v2(dr,ru_s,xprod_aux1)		
		c_aux=-DOT_PRODUCT(xprod_aux1,rv_t)
		E_val=c_aux*(R1_0)
!E_mat(1,1)=mu0*nxcurlSk0a(1,1)-mu1*nxcurlSk1a(1,1)
	  elseif (ipars(2).eq.2) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
    	call my_cross_v2(dr,rv_s,xprod_aux2)
	  	c_aux=-DOT_PRODUCT(xprod_aux2,rv_t)
		E_val=c_aux*(R1_0)
!E_mat(1,2)=mu0*nxcurlSk0a(1,2)-mu1*nxcurlSk1a(1,2)
	  elseif (ipars(2).eq.3) then
		my_exp_1=exp(ima*zk1*r)/(4.0d0*pi)
	    E_val=(my_exp_1)/r
!E_mat(5,6)=-zk0**2*Sk0lambda(1,1)+zk1**2*Sk1lambda(1,1)
	  else
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        call my_cross_v2(dr, ru_s, xprod_aux1)
        c_aux=DOT_PRODUCT(n_t,xprod_aux1)
        E_val=c_aux*(R1_0)	  
!E_mat(6,1)=ep0*mu0*ncurlSk0a(1,1)-ep1*mu1*ncurlSk1a(1,1)
	  endif
	elseif (ipars(1).eq.2) then
      if (ipars(2).eq.1) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
		call my_cross_v2(dr,ru_s,xprod_aux1)
	  	c_aux=DOT_PRODUCT(xprod_aux1,ru_t)
		E_val=c_aux*(R1_0)
!E_mat(2,1)=mu0*nxcurlSk0a(2,1)-mu1*nxcurlSk1a(2,1)
	  elseif (ipars(2).eq.2) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
		call my_cross_v2(dr,rv_s,xprod_aux2)
	  	c_aux=DOT_PRODUCT(xprod_aux2,ru_t)
		E_val=c_aux*(R1_0)	  
!E_mat(2,2)=mu0*nxcurlSk0a(2,2)-mu1*nxcurlSk1a(2,2)
	  elseif (ipars(2).eq.3) then
	  else
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        call my_cross_v2(dr, rv_s, xprod_aux2)
        c_aux=DOT_PRODUCT(n_t,xprod_aux2)
        E_val=c_aux*(R1_0)
!E_mat(6,2)=ep0*mu0*ncurlSk0a(1,2)-ep1*mu1*ncurlSk1a(1,2)
	  endif
	elseif (ipars(1).eq.3) then
      if (ipars(2).eq.1) then	  
	    my_exp_0=exp(ima*zk0*r)/(4.0d0*pi)
	    E_val=(my_exp_0)/r
!E_mat(5,6)=-zk0**2*Sk0lambda(1,1)+zk1**2*Sk1lambda(1,1)
	  elseif (ipars(2).eq.2) then
	  elseif (ipars(2).eq.3) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        c_aux=DOT_PRODUCT(rv_t,dr)
        E_val=-c_aux*(R1_0)
!E_mat(1,6)=nxgradSk0lambda(1,1)-nxgradSk1lambda(1,1)
      else
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
        call my_cross_v2(dr, ru_s, xprod_aux1)
        c_aux=DOT_PRODUCT(n_t,xprod_aux1)
        E_val=c_aux*(-R1_1)	  
!E_mat(6,1)=ep0*mu0*ncurlSk0a(1,1)-ep1*mu1*ncurlSk1a(1,1)
      endif
	else
      if (ipars(2).eq.1) then
	  elseif (ipars(2).eq.2) then
	  elseif (ipars(2).eq.3) then
        R1_0=(ima*zk0*r-1.0d0)/r**3*exp(ima*zk0*r)/(4.0d0*pi)
        c_aux=DOT_PRODUCT(ru_t,dr)
        E_val=c_aux*(R1_0)
!E_mat(2,6)=nxgradSk0lambda(2,1)-nxgradSk1lambda(2,1)
	  else
        R1_1=(ima*zk1*r-1.0d0)/r**3*exp(ima*zk1*r)/(4.0d0*pi)
        call my_cross_v2(dr, rv_s, xprod_aux2)
        c_aux=DOT_PRODUCT(n_t,xprod_aux2)
        E_val=c_aux*(-R1_1)
!E_mat(6,2)=ep0*mu0*ncurlSk0a(1,2)-ep1*mu1*ncurlSk1a(1,2)
	  endif
	endif


return
end subroutine em_auCKi_pec
